/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/actmol.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/forces/forcecomputer.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_sqlite3.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

namespace alexandria
{

typedef std::map<const std::string, int> stringCount;

static bool dump_molecule(FILE                                        *fp,
                          ForceComputer                               *forceComp,
                          CommunicationRecord                         *cr,
                          const gmx::MDLogger                         &mdlog,
                          stringCount                                 *atomTypeCount,
                          stringCount                                 *bccTypeCount,
                          const ForceField                               &pd,
                          MolProp                                     *mp,
                          std::map<std::string, std::vector<double> > *qmap,
                          t_inputrec                                  *inputrec)
{
    alexandria::ACTMol actmol;
    actmol.Merge(mp);
    actmol.setInputrec(inputrec);
    auto imm = actmol.GenerateTopology(fp, &pd, missingParameters::Error,
                                      false);
    if (immStatus::OK == imm)
    {
        std::vector<gmx::RVec> coords = actmol.xOriginal();
        //actmol.symmetrizeCharges(pd, qsymm, nullptr);
        actmol.initQgenResp(&pd, coords, 0.0, 100);
        auto fhandler = actmol.fragmentHandler();
        if (fhandler->topologies().size() == 1)
        {
            std::vector<double> dummy;
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            imm = actmol.GenerateCharges(&pd, forceComp, mdlog, cr,
                                        pd.chargeGenerationAlgorithm(),
                                        qType::ACM,
                                        dummy, &coords, &forces);
            if (immStatus::OK == imm)
            {
                // Add ACM charges
                auto myexp = mp->findExperiment(JobType::OPT);
                if (nullptr == myexp)
                {
                    myexp = mp->findExperiment(JobType::TOPOLOGY);
                }
                if (nullptr != myexp)
                {
                    auto ca       = myexp->calcAtom();
                    auto topatoms = actmol.topology()->atoms();
                    size_t index  = 0;
                    std::vector<double> newq;
                    for(auto atom = ca->begin(); atom < ca->end(); atom++)
                    {
                        // TODO this is not general!
                        double qq = topatoms[index++].charge();
                        if (pd.polarizable())
                        {
                            qq += topatoms[index++].charge();
                        }
                        newq.push_back(qq);
                        atom->AddCharge(qType::ACM, qq);
                    }
                    qmap->insert({fhandler->ids()[0], newq});
                }
            }
        }
    }
    if (immStatus::OK != imm)
    {
        fprintf(fp, "Failed to generate topology for %s. Outcome: %s\n",
                actmol.getMolname().c_str(), immsg(imm));
        return false;
    }
    else
    {
        std::map<MolPropObservable, iqmType> iqm = {
            { MolPropObservable::DIPOLE, iqmType::Both },
            { MolPropObservable::QUADRUPOLE, iqmType::Both },
            { MolPropObservable::POLARIZABILITY, iqmType::Both },
        };
       
        fprintf(fp, "Molecule: %s\n", actmol.getMolname().c_str());
        for(const auto &f : actmol.fragments())
        {
            f.dump(fp);
        }
        actmol.getExpProps(iqm, -1);
        actmol.Dump(fp);
        // Atoms!
        auto &atoms = actmol.topology()->atoms();
        std::vector<Identifier> atomId;
        auto ztype = InteractionType::COULOMB;
        for (size_t i = 0; i < atoms.size(); i++)
        {
            const auto &atype = atoms[i].ffType();
            fprintf(fp, "atom: %2lu  %5s  %5s", i+1, 
                    atoms[i].name().c_str(), atype.c_str());
            Identifier pid(atype);
            atomId.push_back(pid);
            if (pd.hasParticleType(pid))
            {
                auto pIter = pd.findParticleType(pid);
                if (pIter->hasInteractionType(ztype))
                {
                    auto zid = pIter->interactionTypeToIdentifier(ztype);
                    fprintf(fp, "  %s", zid.id().c_str());
                    auto atypeMap = atomTypeCount->find(zid.id());
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
            fprintf(fp, "\n");
        }
        // Bonds!
        auto bctype = InteractionType::BONDCORRECTIONS;
        for (const auto &b : mp->bondsConst())
        {
            int ai = b.aI();
            int aj = b.aJ();
            fprintf(fp, "bcc: %3d  %3d  %5g", ai+1, aj+1, b.bondOrder());
            if (pd.hasParticleType(atomId[ai]) && pd.hasParticleType(atomId[aj]))
            {
                auto pidI = pd.findParticleType(atomId[ai]);
                auto pidJ = pd.findParticleType(atomId[aj]);
                if (pidI->hasInteractionType(ztype) && pidJ->hasInteractionType(ztype))
                {
                    auto zidI = pidI->interactionTypeToIdentifier(ztype);
                    auto zidJ = pidJ->interactionTypeToIdentifier(ztype);
                    Identifier mybond({ zidI.id(), zidJ.id()}, { b.bondOrder() }, CanSwap::No);
                    auto btypeMap   = bccTypeCount->find(mybond.id());
                    bool bondExists = false;
                    auto fs         = pd.findForcesConst(bctype);
                    if (fs.parameterExists(mybond))
                    {
                        fprintf(fp, "  %s", mybond.id().c_str());
                        bondExists = true;
                    }
                    else
                    {
                        Identifier mybond2({ zidJ.id(), zidI.id()}, { b.bondOrder() }, CanSwap::No);
                        mybond = mybond2;
                        btypeMap   = bccTypeCount->find(mybond.id());
                        auto fs = pd.findForcesConst(bctype);
                        if (fs.parameterExists(mybond))
                        {
                            fprintf(fp, "  %s", mybond.id().c_str());
                            bondExists = true;
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
            fprintf(fp, "\n");
        }
    }
    return true;
}

static void monomer2cluster(FILE                                              *fp,
                            std::vector<MolProp>                              *mp,
                            const ForceField                                  &pd,
                            const std::map<std::string, std::vector<double> > &qmap,
                            t_inputrec                                        *inputrec)
{
    // Premature optimization and all that
    for(auto mm = mp->begin(); mm < mp->end(); mm++)
    {
         alexandria::ACTMol actmol;
         actmol.Merge(&(*mm));
         actmol.setInputrec(inputrec);
         auto imm = actmol.GenerateTopology(fp, &pd, missingParameters::Error,
                                           false);
         if (immStatus::OK == imm)
         {
             std::vector<gmx::RVec> coords = actmol.xOriginal();
             //actmol.symmetrizeCharges(pd, qsymm, nullptr);
             actmol.initQgenResp(&pd, coords, 0.0, 100);
             auto fhandler = actmol.fragmentHandler();
             if (fhandler->topologies().size() > 1)
             {
                 // Add ACM charges
                 auto allexp = mm->experiment();
                 if (allexp->empty())
                 {
                     if (fp)
                     {
                         fprintf(fp, "No experiment for %s, removing it.\n",
                                 mm->getMolname().c_str());
                     }
                     mp->erase(mm);
                 }
                 for(auto myexp = allexp->begin(); myexp < allexp->end(); ++myexp)
                 {
                     auto ca    = myexp->calcAtom();
                     size_t idx = 0;
                     bool found = true;
                     for(const auto &id : fhandler->ids())
                     {
                         auto qm = qmap.find(id);
                         found   = qmap.end() != qm;
                         if (!found)
                         {
                             break;
                         }
                         for(size_t iq = 0; iq < qm->second.size(); iq++)
                         {
                             if (idx == ca->size())
                             {
                                 GMX_THROW(gmx::InternalError(gmx::formatString("Size of charge array (%lu) does not match number of atoms (%d) for %s", qm->second.size(), myexp->NAtom(), mm->getMolname().c_str()).c_str()));
                             }
                             (*ca)[idx].AddCharge(qType::ACM, qm->second[iq]);
                             idx++;
                         }
                     }
                     if (!found)
                     {
                         // No support for this dimer, remove it.
                         if (fp)
                         {
                             fprintf(fp, "No charge support for all fragments in %s, removing it.\n",
                                     mm->getMolname().c_str());
                         }
                         mp->erase(mm);
                     }
                 }
             }
         }
    }
}

static void check_mp(const char           *ffname,
                     const char           *logname,
                     bool                  genCharges,
                     std::vector<MolProp> *mp)
{
    alexandria::ForceField pd;
    try
    {
        alexandria::readForceField(ffname, &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    auto inputrec  = new t_inputrec();
    fill_inputrec(inputrec);

    stringCount atomTypeCount;
    stringCount bccTypeCount;

    CommunicationRecord commRec;
    auto forceComp = new ForceComputer();
    auto mdlog     = gmx::MDLogger {};
    std::map<std::string, std::vector<double> > qmap;

    FILE *mylog = gmx_ffopen(logname, "w");
    fprintf(mylog, "Force field file %s\n", ffname);
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
                std::string name = cai.getName();
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
                fprintf(mylog, "%s #C %d #H %d\n",
                        ci.getDatafile().c_str(), 
                        nC, nH);
            }
            if (ci.NAtom() > 0)
            {
                method = ci.getMethod();
                basis  = ci.getBasisset();
            }
            double T = 0;
            auto gp = m->qmProperty(MolPropObservable::DIPOLE, T, JobType::OPT);
            if (gp)
            {
                std::vector<double> mu = gp->getVector();
                name_mu nmu = { ci.getDatafile(), { mu[XX], mu[YY], mu[ZZ] } };
                mus.push_back(nmu);
            }
        
            auto Xcalc = ci.getCoordinates();
            auto Esp   = ci.electrostaticPotentialConst();
            if (Esp.size() >= Xcalc.size() && Xcalc.size() > 1)
            {
                double msd = 0;
                auto xunit = Esp[0].getXYZunit();
                double fac = convertToGromacs(1.0, xunit);
                for(size_t i = 0; i < Xcalc.size(); i++)
                {
                    msd += (gmx::square(Xcalc[i][XX]-fac*Esp[i].getX())+
                            gmx::square(Xcalc[i][YY]-fac*Esp[i].getY())+
                            gmx::square(Xcalc[i][ZZ]-fac*Esp[i].getZ()));
                }
                double rmsd = std::sqrt(msd/Xcalc.size());
                if (rmsd != 0)
                {
                    fprintf(mylog, "%s RMSD coordinates between ESP and QM %g\n",
                            m->getMolname().c_str(), rmsd);
                }
                if (rmsd > 1e-3)
                {
                    for(size_t i = 0; i < Xcalc.size(); i++)
                    {
                        fprintf(mylog, "%2d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                                static_cast<int>(i+1),
                                Xcalc[i][XX], Xcalc[i][YY], Xcalc[i][ZZ],
                                fac*Esp[i].getX(), fac*Esp[i].getY(),
                                fac*Esp[i].getZ());
                    }
                }
            }
        }
        // Check dipoles
        if (debug)
        {
            for(const auto &mi : mus)
            {
                fprintf(debug, "%s %s %.2f %.2f %.2f\n", m->getMolname().c_str(),
                        mi.name.c_str(), mi.mu[XX], mi.mu[YY], mi.mu[ZZ]);
            }
        }

        if (dump_molecule(mylog, forceComp, &commRec, mdlog, &atomTypeCount,
                          &bccTypeCount, pd, &(*m), &qmap, inputrec))
        {
            numberOk++;
        }
        else
        {
            numberFailed++;
            mp->erase(m);
        }
    }
    fprintf(mylog, "Succeed making %d topologies, failed for %d compounds\n",
            numberOk, numberFailed);
    fprintf(mylog, "Statistics\n");
    for(auto &atc : atomTypeCount)
    {
        fprintf(mylog, "atom: %-6s  %5d\n", atc.first.c_str(), atc.second);
    }
    for(auto &bcc : bccTypeCount)
    {
        fprintf(mylog, "bcc: %-12s  %5d\n", bcc.first.c_str(), bcc.second);
    }
    if (genCharges)
    {
        // Copy charges from monomers
        monomer2cluster(mylog, mp, pd, qmap, inputrec);
    }
    fclose(mylog);

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
        "It also can try to make a topology and reports errors doing this. Output is to a log file.",
        "If the optional [TT]-charges[tt] flag is given, this program will generate",
        "charges for all the monomeric compounds in the molprop file,",
        "according to the input force field file and store them",
        "in the output molprop file as ACM (Alexandria Charge Model) charges.",
        "If there are clusters in the input molprop file they will get charges",
        "from the monomers.[PAR]"

    };
    std::vector<t_filenm> fnm =
    {
        { efXML, "-mp",  "data",    ffOPTRDMULT },
        { efXML, "-o",   "allmols", ffWRITE     },
        { efXML, "-ff",  "aff",     ffOPTRD     },
        { efLOG, "-g",   "check",   ffOPTWR     },
        { efDAT, "-db",  "sqlite",  ffOPTRD     }
    };
    int      compress    = 1;
    real     temperature = 298.15;
    bool     forceMerge  = false;
    gmx_bool bcast       = false;
    bool     genCharges  = false;
    int      maxwarn     = 0;
    int      writeNode   = 0;
    t_pargs pa[] =
    {
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML files" },
        { "-maxwarn", FALSE, etINT, {&maxwarn},
          "Will only write output if number of warnings is at most this." },
        { "-temperature", FALSE, etREAL, {&temperature},
          "Temperature for properties to extract from the SQLite database" },
        { "-force", FALSE, etBOOL, {&forceMerge},
          "Force merging compounds with the same name even if not the formula matches" },
        { "-bcast", FALSE, etBOOL, {&bcast},
          "Use broadcast instead of send/receive when running in parallel" },
        { "-wn", FALSE, etINT, {&writeNode},
          "Processor ID to write from if in parallel." },
        { "-charges", FALSE, etBOOL, {&genCharges},
          "Compute charges based on monomers and store them in the output. If there are cluster, e.g. dimers they will get the same charges." }
    };
    std::vector<MolProp>  mpt;
    ForceField            pd;

    gmx_output_env_t     *oenv;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, fnm.size(), fnm.data(),
                           asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    auto fns = opt2fns("-mp", fnm.size(), fnm.data());

    CommunicationRecord cr;
    cr.init(cr.size());
    auto comm = MPI_COMM_WORLD;
    int root  = 0;
    if (cr.isMaster())
    {
        int nwarn   = merge_xml(fns, &mpt, nullptr, nullptr, nullptr, forceMerge);
        int mptsize = mpt.size();
        if (nwarn <= maxwarn)
        {
            ReadSqlite3(opt2fn_null("-db", fnm.size(), fnm.data()), &mpt, temperature);

            printf("Read %d molecules\n", mptsize);
        }
        else
        {
            printf("Too many warnings (%d), not generating output\n", nwarn);
            mpt.clear();
        }
            
        if (bcast)
        {
            cr.bcast(&mptsize, comm);
            for(auto &mm : mpt)
            {
                mm.BroadCast(&cr, root, comm);
            }
        }
        else
        {
            for(int dest = 1; dest < cr.size(); dest++)
            {
                cr.send_int(dest, mptsize);
                for(const auto &mm : mpt)
                {
                    mm.Send(&cr, dest);
                }
            }
        }
    }
    else
    {
        if (bcast)
        {
            int nmpt;
            cr.bcast(&nmpt, comm);
            for(int i = 0; i < nmpt; i++)
            {
                MolProp mp;
                mp.BroadCast(&cr, root, comm);
                mpt.push_back(mp);
            }
        }
        else
        {
            int nmpt = cr.recv_int(cr.superior());
            for(int i = 0; i < nmpt; i++)
            {
                MolProp mp;
                mp.Receive(&cr, cr.superior());
                mpt.push_back(mp);
            }
        }
    }
    auto ffname = opt2fn_null("-ff", fnm.size(), fnm.data());
    if (ffname)
    {
        check_mp(ffname, opt2fn("-g", fnm.size(), fnm.data()), genCharges, &mpt);
    }
    if (writeNode == cr.rank())
    {
        MolPropWrite(opt2fn("-o", fnm.size(), fnm.data()), mpt, compress);
    }
    return 0;
}

} // namespace alexandria
