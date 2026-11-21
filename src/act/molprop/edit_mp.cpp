/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#include "actpre.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/basics/allmols.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/particletype.h"
#include "act/forces/forcecomputer.h"
#include "act/molprop/experiment.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_json.h"
#include "act/molprop/molprop_sqlite3.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/molprop/molpropobservable.h"
#include "act/molprop/topologyentry.h"
#include "act/statistics/statistics.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_low.h"
#include "act/utility/jsontree.h"
#include "act/utility/output_environment.h"
#include "act/utility/units.h"
#include "edit_mp.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

//! \brief Map to count the occurence of strings, e.g. atom names.
typedef std::map<const std::string, int> stringCount;

//! \brief Dump important parts of a MolProp to the msghandler
static bool dump_molecule(MsgHandler        *msghandler,
                          ForceComputer     *forceComp,
                          stringCount       *atomTypeCount,
                          stringCount       *bccTypeCount,
                          ForceField        &pd,
                          MolProp           *mp)
{
    alexandria::ACTMol actmol;
    actmol.Merge(mp);
    actmol.GenerateTopology(msghandler, &pd, missingParameters::Error);
    if (msghandler->ok())
    {
        std::vector<gmx::RVec> coords = actmol.xOriginal();
        //! \todo check whether this is needed.
        //actmol.symmetrizeCharges(pd, qsymm, nullptr);
        std::map<MolPropObservable, iqmType> iqm = {
            { MolPropObservable::CHARGE, iqmType::QM }
        };
        actmol.getExpProps(msghandler, &pd, iqm, 0.0, 100);
        auto fhandler = actmol.fragmentHandler();
        if (fhandler->topologies().size() == 1)
        {
            const std::vector<double> dummy;
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            actmol.generateCharges(msghandler, &pd, forceComp,
                                   pd.chargeGenerationAlgorithm(), &coords, &forces);
        }
    }
    if (!msghandler->ok())
    {
        return false;
    }
    else
    {
        std::map<MolPropObservable, iqmType> iqm = {
            { MolPropObservable::DIPOLE, iqmType::Both },
            { MolPropObservable::QUADRUPOLE, iqmType::Both },
            { MolPropObservable::POLARIZABILITY, iqmType::Both },
        };
        // Fetch TextWriter object
        auto tw = msghandler->tw();
        tw->writeStringFormatted("Molecule: %s\n", actmol.getMolname().c_str());
        actmol.getExpProps(msghandler, &pd, iqm, -1);
        actmol.Dump(tw);
        // Atoms!
        auto &atoms = actmol.topology()->atoms();
        std::vector<Identifier> atomId;
        for (size_t i = 0; i < atoms.size(); i++)
        {
            const auto &atype = atoms[i].ffType();
            tw->writeStringFormatted("atom: %2lu  %5s  %5s", i+1, 
                    atoms[i].name().c_str(), atype.c_str());
            const Identifier pid(atype);
            atomId.push_back(pid);
            if (pd.hasParticleType(pid))
            {
                auto pIter = pd.findParticleType(pid);
                // Each particle type has its own charge distribution
                if (pIter->hasParameter(zetaName))
                {
                    const auto zid = pIter->id();
                    tw->writeStringFormatted("  %s", zid.id().c_str());
                    const auto atypeMap = atomTypeCount->find(zid.id());
                    if (atypeMap == atomTypeCount->end())
                    {
                        atomTypeCount->insert(std::pair<const std::string, int>(zid.id(), 1));
                    }
                    else
                    {
                        atypeMap->second += 1;
                    }
                }
            }
            tw->writeStringFormatted("\n");
        }
        // Bonds! We cannot use the bonds in the molprop anymore
        // since there may be vsites or shells and the numbering
        // is messed up.
        auto bctype = InteractionType::BONDCORRECTIONS;
        const auto &mybonds = actmol.topology()->entry(InteractionType::BONDS);
        for (size_t i = 0; i < mybonds.size(); i++)
        {
            const auto mybond = static_cast<const Bond &>(*mybonds[i]->self());
            const auto &aa = mybond.atomIndices();
            if (aa.size() != 2)
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Bond with %zu atoms", aa.size()).c_str()));
            }
            int    ai = aa[0];
            int    aj = aa[1];
            double bo = mybond.bondOrders()[0];
            tw->writeStringFormatted("bcc: %3d  %3d  %5g", ai+1, aj+1, bo);
            if (pd.hasParticleType(atomId[ai]) && pd.hasParticleType(atomId[aj]))
            {
                auto pidI = pd.findParticleType(atomId[ai]);
                auto pidJ = pd.findParticleType(atomId[aj]);
                // The bond corrections are identified by the ACM types, see MolGen
                auto acmtype = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
                // A force field without bond charge corrections, such as one with
                // fixed charges, has none to count
                if (pd.interactionPresent(bctype) &&
                    pidI->hasInteractionType(acmtype) && pidJ->hasInteractionType(acmtype))
                {
                    const auto zidI = pidI->interactionTypeToIdentifier(acmtype);
                    const auto zidJ = pidJ->interactionTypeToIdentifier(acmtype);
                    Identifier mybond({ zidI.id(), zidJ.id()}, { bo }, CanSwap::No);
                    auto btypeMap   = bccTypeCount->find(mybond.id());
                    bool bondExists = false;
                    const auto &fs         = pd.findForcesConst(bctype);
                    if (fs.parameterExists(mybond))
                    {
                        tw->writeStringFormatted("  %s", mybond.id().c_str());
                        bondExists = true;
                    }
                    else
                    {
                        const Identifier mybond2({ zidJ.id(), zidI.id()}, { bo }, CanSwap::No);
                        mybond = mybond2;
                        btypeMap   = bccTypeCount->find(mybond.id());
                        const auto &fs = pd.findForcesConst(bctype);
                        if (fs.parameterExists(mybond))
                        {
                            tw->writeStringFormatted("  %s", mybond.id().c_str());
                            bondExists = true;
                        }
                        else
                        {
                            tw->writeStringFormatted("  N/A");
                        }
                    }
                    if (bondExists)
                    {
                        if (btypeMap == bccTypeCount->end())
                        {
                            bccTypeCount->insert(std::pair<const std::string, int>(mybond.id(), 1));
                        }
                        else
                        {
                            btypeMap->second += 1;
                        }
                    }
                }
            }
            tw->writeStringFormatted("\n");
        }
    }
    return true;
}

/*! \brief Check the contents of molprop structures
 * Making the topologies takes the experiments from the compounds, so
 * they should not be used afterwards.
 * \param[in]    msghandler For the results
 * \param[in]    ffname     The force field file
 * \param[inout] mp         The compounds
 * \return the names of the compounds for which no topology could be made
 */
static std::set<std::string> check_mp(MsgHandler           *msghandler,
                                      const char           *ffname,
                                      std::vector<MolProp> *mp)
{
    std::set<std::string> failed;
    alexandria::ForceField pd;
    try
    {
        alexandria::readForceField(ffname, &pd, msghandler);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    stringCount atomTypeCount;
    stringCount bccTypeCount;

    auto forceComp = new ForceComputer();

    if (msghandler)
    {
        auto tw = msghandler->tw();
        tw->writeStringFormatted("Force field file %s\n", ffname);
        int numberOk = 0, numberFailed = 0;
        for (auto m = mp->begin(); m < mp->end(); ++m)
        {
            typedef struct
            {
                std::string name;
                rvec        mu;
            } name_mu;
            std::string basis, method;
            std::vector<name_mu> mus;
            for (auto &ci : m->experimentConst())
            {
                int nH = 0, nC = 0;
                for (auto &cai : ci.calcAtomConst())
                {
                    const std::string &name = cai.getName();
                    if (name.compare("H") == 0)
                    {
                        nH++;
                    }
                    else if (name.compare("C") == 0)
                    {
                        nC++;
                    }
                }
                if (nC > 0 && nH == 0)
                {
                     tw->writeStringFormatted("%s #C %d #H %d\n",
                                              ci.getDatafile().c_str(), 
                                              nC, nH);
                }
                if (ci.NAtom() > 0)
                {
                    method = ci.getMethod();
                    basis  = ci.getBasisset();
                }
                double T = 0;
                // Dimers, for instance, have no dipole or potential of an optimized structure
                if (m->hasQMProperty(MolPropObservable::DIPOLE, T, JobType::OPT))
                {
                    auto &gp = m->qmProperty(MolPropObservable::DIPOLE, T, JobType::OPT);
                    if (gp)
                    {
                        std::vector<double> mu = gp->getVector();
                        const name_mu nmu = { ci.getDatafile(), { static_cast<real>(mu[XX]), static_cast<real>(mu[YY]), static_cast<real>(mu[ZZ]) } };
                        mus.push_back(nmu);
                    }
                }
                if (!m->hasQMProperty(MolPropObservable::POTENTIAL, T, JobType::OPT))
                {
                    continue;
                }
                auto Xcalc = ci.getCoordinates();
                auto &gp = m->qmProperty(MolPropObservable::POTENTIAL, T, JobType::OPT);
                if (gp)
                {
                    const ElectrostaticPotential *ep = static_cast<ElectrostaticPotential *>(gp.get());
                    auto xyz = ep->xyz();
                    if (xyz.size() >= Xcalc.size() && Xcalc.size() > 1)
                    {
                        double msd = 0;
                        const auto xunit = ep->getXYZunit();
                        double fac = convertToGromacs(1.0, xunit);
                        for(size_t i = 0; i < Xcalc.size(); i++)
                        {
                            msd += (gmx::square(Xcalc[i][XX]-fac*xyz[i][XX])+
                                    gmx::square(Xcalc[i][YY]-fac*xyz[i][YY])+
                                    gmx::square(Xcalc[i][ZZ]-fac*xyz[i][ZZ]));
                        }
                        double rmsd = std::sqrt(msd/Xcalc.size());
                        if (rmsd != 0)
                        {
                            tw->writeStringFormatted("%s RMSD coordinates between ESP and QM %g\n",
                                                     m->getMolname().c_str(), rmsd);
                        }
                        if (rmsd > 1e-3)
                        {
                            for(size_t i = 0; i < Xcalc.size(); i++)
                            {
                                tw->writeStringFormatted("%2zu %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                                                         i+1,
                                                         Xcalc[i][XX], Xcalc[i][YY], Xcalc[i][ZZ],
                                                         fac*xyz[i][XX], fac*xyz[i][YY], fac*xyz[i][ZZ]);
                            }
                        }
                    }
                }
            }
            // Check dipoles
            if (msghandler->debug())
            {
                for(const auto &mi : mus)
                {
                    msghandler->msg(ACTStatus::Debug,
                                    gmx::formatString("%s %s %.2f %.2f %.2f\n", m->getMolname().c_str(),
                                                      mi.name.c_str(), mi.mu[XX], mi.mu[YY], mi.mu[ZZ]));
                }
            }
            
            if (dump_molecule(msghandler, forceComp, &atomTypeCount,
                              &bccTypeCount, pd, &(*m)))
            {
                numberOk++;
            }
            else
            {
                numberFailed++;
                failed.insert(m->getMolname());
            }
        }
        tw->writeStringFormatted("Succeed making %d topologies, failed for %d compounds\n",
                                 numberOk, numberFailed);
        tw->writeStringFormatted("Statistics\n");
        for(auto &atc : atomTypeCount)
        {
            tw->writeStringFormatted("atom: %-6s  %5d\n", atc.first.c_str(), atc.second);
        }
        for(auto &bcc : bccTypeCount)
        {
            tw->writeStringFormatted("bcc: %-12s  %5d\n", bcc.first.c_str(), bcc.second);
        }
    }
    return failed;
}

//! \brief Generate energy histogram.
static void gen_ehist(MsgHandler                 *msghandler,
                      const std::vector<MolProp> *mpt,
                      gmx_output_env_t           *oenv,
                      real                        ewarnLow,
                      real                        ewarnHi)
{
    const std::set<MolPropObservable> mpset = { MolPropObservable::DELTAE0, MolPropObservable::INTERACTIONENERGY };
    auto tw = msghandler->tw();
    for(auto mp = mpt->begin(); mp < mpt->end(); ++mp)
    {
        std::map<MolPropObservable, gmx_stats> histo;
        for(auto mps : mpset)
        {
            histo[mps] = gmx_stats();
            for(auto &eee : mp->experimentConst())
            {
                if (eee.hasProperty(mps))
                {
                    for (auto &prop : eee.propertyConst(mps))
                    {
                        histo[mps].add_point(prop->getValue());
                    }
                }
            }
            if (histo[mps].get_npoints() > 0)
            {
                std::vector<double> x, y;
                real binwidth = 1;
                int  nbins    = 0;
                histo[mps].make_histogram(binwidth, &nbins, eHisto::Y, false, &x, &y);
                const auto filename = gmx::formatString("%s-%s.xvg", mp->getMolname().c_str(), mpo_name(mps));
                FILE *fp     = xvgropen(filename.c_str(), "Energy distribution in molprop", "Energy (kJ/mol)", "(a.u.)", oenv);
                bool warnLow = false;
                bool warnHi  = false;
                for(size_t i = 0; i < x.size(); i++)
                {
                    fprintf(fp, "%10g  %10g\n", x[i], y[i]);
                    if (x[i] < ewarnLow)
                    {
                        warnLow = true;
                    }
                    else if (x[i] > ewarnHi)
                    {
                        warnHi = true;
                    }
                }
                xvgrclose(fp);
                if (warnLow)
                {
                    tw->writeStringFormatted("Warning: low energies encountered for %s\n", mp->getMolname().c_str());
                }
                if (warnHi)
                {
                    tw->writeStringFormatted("Warning: high energies encountered for %s\n", mp->getMolname().c_str());
                }
                tw->writeStringFormatted("Range of energies for %s : %g - %g%s\n",
                                         mp->getMolname().c_str(), x[0], x[x.size()-1],
                                         x[0] > 0 ? " ALL-POSITIVE" : "");
            }
        }
    }
}

/*! \brief Update internal molecule information.
 * This will use an external data file, alexandria.json, containing
 * information about > 9000 molecules.
 */
static void updateMolInfo(MsgHandler           *msghandler,
                          std::vector<MolProp> *mpt)
{
    const AlexandriaMols amols;
    auto tw = msghandler->tw();
    for(auto mp = mpt->begin(); mp < mpt->end(); ++mp)
    {
        auto frags = mp->fragments();
        if (frags.size() == 1)
        {
            auto amol = amols.findInChi(frags[0].inchi());
            if (amol)
            {
                if (tw)
                {
                    tw->writeStringFormatted("Updating %s\n", mp->getMolname().c_str());
                }
                mp->SetIupac(amol->iupac);
                mp->SetCas(amol->cas);
                mp->SetCid(std::to_string(amol->csid));
                mp->SetInchi(amol->inchi);
            }
        }
    }
}

/*! \brief Function to modify molprop files in a number of ways.
 * \param[in] argc Number of command line arguments
 * \param[in] argv The command line arguments
 * \return 0 if ok, another value otherwise.
 */
void EditMPOptions::check() const
{
    if (inputFiles.empty())
    {
        GMX_THROW(gmx::InvalidInputError("There is no molprop file to read (-mp)"));
    }
    if (outputFile.empty())
    {
        GMX_THROW(gmx::InvalidInputError("There is no molprop file to write (-o)"));
    }
    if (!std::isfinite(temperature) || temperature <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -temperature should be more than zero",
                                                           temperature)));
    }
    if (!std::isfinite(energyWarnLow))
    {
        GMX_THROW(gmx::InvalidInputError("The value of -ewarnLow is not a number"));
    }
    if (!std::isfinite(energyWarnHigh))
    {
        GMX_THROW(gmx::InvalidInputError("The value of -ewarnHi is not a number"));
    }
}

void runEditMP(MsgHandler *msghandler, const EditMPOptions &options, gmx_output_env_t *oenv)
{
    options.check();
    std::vector<MolProp> mpt;
    const auto          &fns = options.inputFiles;
    merge_xml(msghandler, fns, &mpt);
    int mptsize = mpt.size();
    if (msghandler->ok())
    {
        ReadSqlite3(msghandler, options.databaseFile.empty() ? nullptr : options.databaseFile.c_str(),
                    &mpt, options.temperature);

        printf("Read %d molecules\n", mptsize);
    }
    else
    {
        fprintf(stderr, "Too many warnings, not generating output\n");
        mpt.clear();
    }
    const char *ffname = options.forceFieldFile.empty() ? nullptr : options.forceFieldFile.c_str();

    auto tw = msghandler->tw();
    if (options.molinfo)
    {
        updateMolInfo(msghandler, &mpt);
    }
    const char *molpropout = options.outputFile.c_str();
    if (ffname)
    {
        printf("Since you provided a force field file and a log file name, I will now check the compounds.\n");
        // The check takes the experiments from the compounds, so it uses a
        // copy, the compounds for which it fails are not written.
        std::vector<MolProp> checked;
        merge_xml(msghandler, fns, &checked);
        auto failed = check_mp(msghandler, ffname, &checked);
        mpt.erase(std::remove_if(mpt.begin(), mpt.end(),
                                 [&failed](const MolProp &mp) { return failed.count(mp.getMolname()) > 0; }),
                  mpt.end());
    }
    // The result of the validation for the JSON output
    JsonTree jvalidation("validation");
    jvalidation.addBool("checked", false);
    if (options.validate)
    {
        // The problems are errors of the message handler, but the molecules
        // are written anyway, without the experiments that were removed
        std::vector<ACTRecordedMessage> messages;
        msghandler->record(&messages, ACTStatus::Warning);
        const auto validation = validateMolProps(msghandler, &mpt);
        msghandler->record(nullptr);
        jvalidation   = molPropValidationToJson(validation, messages);
        auto nproblem = validation.nproblem;
        printf("Validation of %zu molecules: %zu problem(s)%s\n", mpt.size(), nproblem,
               nproblem > 0 ? gmx::formatString(", see %s", msghandler->filename().c_str()).c_str() : "");
    }
    if (options.energyHistogram)
    {
        // A program without a command line has no output environment for the histograms
        const OutputEnvironment output(oenv);
        gen_ehist(msghandler, &mpt, output.get(), options.energyWarnLow, options.energyWarnHigh);
    }
    if (tw)
    {
        printf("Please check %s for analyses.\n", msghandler->filename().c_str());
    }
    printf("Will store %zu molprops in output file %s\n", mpt.size(), molpropout);
    MolPropWrite(molpropout, mpt, options.compress);
    if (!options.jsonFile.empty())
    {
        // What was written, for programs that show it
        JsonTree jtree("edit_mp");
        jtree.addObject(molPropSummaryToJson(mpt));
        JsonTree jmols = JsonTree::array("molecules");
        for (const auto &mp : mpt)
        {
            jmols.addObject(molPropToJson(mp));
        }
        jtree.addObject(jmols);
        jtree.addObject(jvalidation);
        jtree.write(options.jsonFile, true);
    }
}

int edit_mp(int argc, char *argv[])
{
    static const char               *desc[] =
    {
        "edit_mp manipulates molprop files. It can read multiple molprop files and merges",
        "the molecule descriptions into a single new file. By specifying the [TT]-db[TT]", 
        "option additional experimental information will be read from a SQLite3 database.[PAR]",
        "The program can run in paralllel to read a file on one processer, then send it over",
        "an MPI connection to one or more other processors to write. In this manner the MPI transfer",
        "software in ACT can be tested.[PAR]",
        "edit_mp can check calculations for missing hydrogens and inconsistent dipoles if a force field file is given.", 
        "It also can try to make a topology and reports errors doing this. Output is to a log file.[PAR]",
        "With [TT]-validate[tt] the molecules are checked for consistency: all experiments",
        "of a molecule need the atoms in the same order, those that differ are removed;",
        "each compound (a fragment, identified by its InChI) needs its atoms in the same",
        "order in all molecules; the atoms of each fragment have to match its formula; the",
        "fragments of a dimer may not be bonded to each other; and a dimer A#B may not occur",
        "as B#A as well, nor have its fragments in the order of B#A. The problems are",
        "printed and written to the log file.[PAR]",
        "With [TT]-json[tt] a description of the molecules that are written goes to",
        "the given file in json format, edit_mp.json by default: for each molecule",
        "its fragments, the number of experiments and the number of values of each",
        "property, the totals, and the result of [TT]-validate[tt] with the names of",
        "the molecules that have problems and the messages."

    };
    std::vector<t_filenm> fnm =
    {
        { efXML, "-mp",  "data",    ffOPTRDMULT },
        { efXML, "-o",   "allmols", ffWRITE     },
        { efXML, "-ff",  "aff",     ffOPTRD     },
        { efDAT, "-db",  "sqlite",  ffOPTRD     },
        { efJSON, "-json", "edit_mp", ffOPTWR    }
    };
    bool     compress    = false;
    real     temperature = 298.15;
    bool     energyHisto = false;
    bool     molinfo     = true;
    real     ewarnLow    = -20;
    real     ewarnHi     = 100;
    bool     validate    = false;
    std::vector<t_pargs> pa =
    {
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML files" },
        { "-temperature", FALSE, etREAL, {&temperature},
          "Temperature for properties to extract from the SQLite database" },
        { "-molinfo", FALSE, etBOOL, {&molinfo},
          "Update iupac, cas etc. using the Alexandria Molecule Database" },
        { "-ehisto", FALSE, etBOOL, {&energyHisto},
          "Make a histogram of the energy distribution per molecule or complex." },
        { "-ewarnLow", FALSE, etREAL, {&ewarnLow},
          "Print a warning if energy is lower than this number (kJ/mol)" },
        { "-ewarnHi", FALSE, etREAL, {&ewarnHi},
          "Print a warning if energy is higher than this number (kJ/mol)" },
        { "-validate", FALSE, etBOOL, {&validate},
          "Check the order of the molecules and atoms in all experiments and compounds, see above." }
    };
    MsgHandler msghandler;
    gmx_output_env_t     *oenv;
    msghandler.addOptions(&pa, &fnm, "edit_mp.log");

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, fnm.size(), fnm.data(),
                           pa.size(), pa.data(), sizeof(desc)/sizeof(desc[0]), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }
    CommunicationRecord cr;
    cr.init(cr.size());
    msghandler.optionsFinished(fnm, &cr);

    // Only the master edits the compounds
    if (!cr.isMaster())
    {
        return 0;
    }
    try
    {
        EditMPOptions options;
        const auto    fns = opt2fns("-mp", fnm.size(), fnm.data());
        options.inputFiles.assign(fns.begin(), fns.end());
        options.outputFile  = opt2fn("-o", fnm.size(), fnm.data());
        const auto optional = [&fnm](const char *option)
        {
            const char *name = opt2fn_null(option, fnm.size(), fnm.data());
            return std::string(name ? name : "");
        };
        options.forceFieldFile  = optional("-ff");
        options.databaseFile    = optional("-db");
        options.jsonFile        = optional("-json");
        options.compress        = compress;
        options.temperature     = temperature;
        options.molinfo         = molinfo;
        options.energyHistogram = energyHisto;
        options.energyWarnLow   = ewarnLow;
        options.energyWarnHigh  = ewarnHi;
        options.validate        = validate;
        runEditMP(&msghandler, options, oenv);
    }
    catch (const gmx::InvalidInputError &ex)
    {
        fprintf(stderr, "%s\n", ex.what());
        return 1;
    }
    return 0;
}

} // namespace alexandria
