/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024,2025
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
    
#include "compound_reader.h"

#include "act/alexandria/babel_io.h"
#include "act/alexandria/fetch_charges.h"
#include "act/basics/msg_handler.h"
#include "act/molprop/molprop_xml.h"

namespace alexandria
{

static std::vector<const char *> crDesc = {
    "It is highly recommended to provide the [TT]-charges[tt] option",
    "with a molprop file that will be used to generate charges for the",
    "file specified with the [TT]-f[tt] option.[PAR]"
    "Alternatively, you can use the [TT]-db 'molecule(s)'[tt]",
    "option to extract one or more compounds from the molprop file.[PAR]",
};

void CompoundReader::addOptions(std::vector<t_pargs>      *pargs,
                                std::vector<t_filenm>     *filenm,
                                std::vector<const char *> *desc)
{
    std::vector<t_filenm>     fnm = {
        { efXML, "-charges", "molprop",    ffOPTRD }
    };
    for (const auto &fn : fnm)
    {
        filenm->push_back(fn);
    }
    
    for (const auto &cr : crDesc)
    {
        desc->push_back(cr);
    }

    std::vector<t_pargs> mypargs = {
        { "-f",      FALSE, etSTR,  {&filename_},
          "Molecular structure file in e.g. pdb format" },
        { "-generateCharges", FALSE, etBOOL, {&genCharges_},
          "Generate charges for your compound(s) using the Alexandria Charge Model." },
        { "-db",     FALSE, etSTR,  {&dbname_},
          "Read one or more molecules from the database rather than from a file. To specify multiple molecules please use quotes, e.g. [TT]-db[tt] 'water methane ammonia'." },
        { "-oneH", FALSE, etBOOL, {&oneH_},
          "Map all different hydrogen atom types back to H, mainly for debugging." },
        { "-qtot",   FALSE, etREAL, {&qtot_},
          "Combined charge of the molecule(s). This will be taken from the input file by default, but that is not always reliable. If the -qcustom flag is supplied, that will be used instead." },
        { "-qqm",    FALSE, etSTR,  {&qqm_},
          "Use a method from quantum mechanics that needs to be present in the input file. Either ESP, Hirshfeld, CM5 or Mulliken may be available." },
        { "-qcustom", FALSE, etSTR, {&qcustom_}, 
          "Here a quoted string of custom charges can be provided such that a third party source can be used. It is then possible to generate multipoles and compare the ESP to a quantum chemistry result. The number of charges provided must match the number of particles (including shells if present in the force field used)." }
    };
    for(const auto &mp : mypargs)
    {
        pargs->push_back(mp);
    }
}

void CompoundReader::optionsOK(MsgHandler                  *msghandler,
                               const std::vector<t_filenm> &filenm)
{
    if (strlen(filename_) == 0 && strlen(dbname_) == 0)
    {
        msghandler->msg(ACTStatus::Error,
                        "Please provide either a filename (-f flag) with a compound/dimer structure to optimize or simulate, or a compound/dimer from the charges database (-db flag).");
        return;
    }
    if (strlen(filename_) != 0 && strlen(dbname_) != 0)
    {
        msghandler->msg(ACTStatus::Error,
                        "Please provide just one of the flags -f and -db.");
        return;
    }
    const char *qfn = opt2fn_null("-charges", filenm.size(), filenm.data());
    if (qfn)
    {
        qmapfn_.assign(qfn);
    }
    if (!qmapfn_.empty() && (strlen(qcustom_) > 0 || strlen(qqm_) > 0 || genCharges_))
    {
        msghandler->msg(ACTStatus::Error,
                        "If you provide a charge map please do not provide a custom charge string or a QM charge selection or the generateCharges flag at the same time.");
        return;
    }
    if (strlen(dbname_) > 0 && qmapfn_.empty())
    {
        msghandler->msg(ACTStatus::Error,
                        "Please provide -charges flag in conjunction with -db.");
        return;
    }
    if (strlen(qcustom_) > 0 && strlen(qqm_) > 0)
    {
        msghandler->msg(ACTStatus::Error,
                        "Please do not provide both the custom charges and the QM charge type to read.");
        return;
    }
    if (genCharges_ && (strlen(qcustom_) > 0 || strlen(qqm_) > 0))
    {
        msghandler->msg(ACTStatus::Error,
                        "Please do not provide both the generateCharges flag and either custom charges or the QM charge type to read.");
        return;
    }
}

void CompoundReader::setCharges(MsgHandler          *msghandler,
                                ForceField          &pd,
                                ACTMol              *mol,
                                const chargeMap     &qmap,
                                const ForceComputer *forceComp,
                                bool                 warnQtot)
{
    if (mol->totalCharge() != qtot_ && warnQtot)
    {
        msghandler->msg(ACTStatus::Warning,
                        gmx::formatString("Detected total charge %d, command line says %g.",
                                          mol->totalCharge(), qtot_));
    }
    
    mol->GenerateTopology(msghandler, &pd, missingParameters::Error);

    if (msghandler->ok())
    {
        std::vector<gmx::RVec> coords = mol->xOriginal();
        auto fragments  = mol->fragmentHandler();
        if (!qmap.empty())
        {
            fragments->setCharges(msghandler, qmap);
            if (msghandler->ok())
            {
                // Copy charges to the high-level topology as well
                fragments->fetchCharges(mol->atoms());
            }
            else
            {
                msghandler->msg(ACTStatus::Error,
                                gmx::formatString("CompoundReader: not all compounds present in the charge map %s", qmapfn_.c_str()));
                return;
            }
        }
        else
        {
            std::vector<gmx::RVec> forces(mol->atomsConst().size());

            std::vector<double> myq;
            auto alg   = pd.chargeGenerationAlgorithm();
            auto qtype = qType::Calc;
            if (strlen(qcustom_) > 0)
            {
                auto mycharges = gmx::splitString(qcustom_);
                for(auto &q : mycharges)
                {
                    myq.push_back(my_atof(q.c_str(), "custom q"));
                }
                alg = ChargeGenerationAlgorithm::Custom;
            }
            else if (strlen(qqm_) > 0)
            {
                alg   = ChargeGenerationAlgorithm::Read;
                qtype = stringToQtype(qqm_);
            }
            if (!genCharges_)
            {
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("Using %s to generate charges. It is recommended to use a charge database instead of this option.\n", chargeGenerationAlgorithmName(alg).c_str()));
            }
            mol->GenerateCharges(msghandler, &pd, forceComp, alg, qtype, myq, &coords, &forces,
                                 msghandler->verbose());
        }
    }
}

void CompoundReader::readFile(MsgHandler *msghandler,
                              ForceField &pd,
                              ACTMol     *mol)
{
    matrix box;
    clear_mat(box);
    if (strlen(filename_) == 0)
    {
        msghandler->msg(ACTStatus::Error, "Empty filename");
        return;
    }
    std::vector<MolProp> mps;
    double               qtot_babel = qtot_;
    std::string          method, basis;
    int                  maxpot = 100;
    int                  nsymm  = 1;
    bool                 addHydrogen = false;
    if (!readBabel(&pd, filename_, &mps, molnm_, molnm_, "", &method,
                   &basis, maxpot, nsymm, "Opt", userQtot(), &qtot_babel,
                   addHydrogen, box, oneH_))
    {
        msghandler->msg(ACTStatus::Error,
                        gmx::formatString("Reading %s failed.\n", filename_));
        return;
    }

    if (mps.size() > 1)
    {
        msghandler->msg(ACTStatus::Warning,
                        gmx::formatString("will only use the first compound (out of %zu) in %s\n",
                                          mps.size(), filename_));
    }
    if (mps.size() == 0)
    {
        msghandler->msg(ACTStatus::Error,
                        gmx::formatString("Failed to import coordinate file %s using OpenBabel", filename_));
        return;
    }
    mol->Merge(&mps[0]);
}

std::vector<ACTMol> CompoundReader::read(MsgHandler          *msghandler,
                                         ForceField          &pd,
                                         const ForceComputer *forceComp)
{
    std::vector<ACTMol>   mols;
    std::set<std::string> lookup;
    // Try reading from a file first
    bool readCoordinates = strlen(filename_) > 0;
    if (readCoordinates)
    {
        ACTMol mol;
        readFile(msghandler, pd, &mol);
        if (msghandler->ok())
        {
            auto fp = mol.fragmentPtr();
            if (fp)
            {
                for(auto ic = fp->begin(); ic < fp->end(); ++ic)
                {
                    lookup.insert(ic->iupac());
                }
            }
            else
            {
                lookup.insert(mol.getMolname());
            }
            mols.push_back(mol);
        }
    }
    else if (strlen(dbname_) > 0)
    {
        // Determine what compounds the usesr selected
        for(const auto &mymol : split(dbname_, ' '))
        {
            lookup.insert(mymol);
        }
    }
    if (lookup.empty())
    {
        msghandler->msg(ACTStatus::Info,
                        gmx::formatString("CompoundReader will include all compounds from %s\n.", qmapfn_.c_str()));
    }
    else
    {
        std::string msg("CompoundReader found the following compounds:");
        for(const auto &lu : lookup)
        {
            msg += " " + lu;
        }
        if (strlen(dbname_) > 0)
        {
            msg += " in the charges molprop.";
        }
        else
        {
            msg += " in ";
            msg += filename_;
        }
        msghandler->msg(ACTStatus::Info, msg);
    }
    chargeMap qmap;
    if (!qmapfn_.empty())
    {
        std::vector<MolProp> mps;
        MolPropRead(qmapfn_.c_str(), &mps);
        qmap = fetchChargeMap(msghandler, &pd, forceComp, mps, lookup);
        msghandler->msg(ACTStatus::Info,
                        gmx::formatString("CompoundReader read %lu out of %lu entries into charge map from %s\n",
                                          qmap.size(), lookup.size(), qmapfn_.c_str()));
        
        // Throw away those compounds that are not in the selection
        if (!lookup.empty() && !readCoordinates)
        {
            for(auto mp = mps.begin(); mp < mps.end(); )
            {
                if (lookup.find(mp->getMolname()) == lookup.end() &&
                    lookup.find(mp->getIupac()) == lookup.end())
                {
                    mp = mps.erase(mp);
                }
                else
                {
                    msghandler->msg(ACTStatus::Info,
                                    gmx::formatString("Keeping %s (%s) from molprop file %s\n",
                                                      mp->getMolname().c_str(), mp->getIupac().c_str(), qmapfn_.c_str()));
                    ++mp;
                }
            }
        }
        if (!readCoordinates)
        {
            for(auto mp : mps)
            {
                ACTMol mm;
                mm.Merge(&mp);
                mols.push_back(mm);
            }
        }
    }
    // Did we find any molecule?
    if (mols.empty())
    {
        msghandler->msg(ACTStatus::Error, "Couldn't find any molecule!");
    }
    else
    {
        // Only warn about the total charge for monomers
        // since we cannot know how charge is divided between molecules.
        bool warnQtot = mols.size() == 1;
        for(auto mol = mols.begin(); mol < mols.end(); )
        {
            setCharges(msghandler, pd, &(*mol), qmap, forceComp, warnQtot);
            // Load all the other properties of the compounds as well
            mol->getExpProps(msghandler, &pd, {});
            if (!msghandler->ok())
            {
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("CompoundReader could not determine charges for '%s' from '%s'\n",
                                                  mol->getMolname().c_str(), filename_));
                // Prevent false positives, delete compound and reset status
                msghandler->resetStatus();
                mol = mols.erase(mol);
            }
            else
            {
                mol++;
            }
        }
    }
    return mols;
}

}
