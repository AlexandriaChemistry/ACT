/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024-2026
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

#include "act/alexandria/molselect.h"
#include "act/basics/msg_handler.h"
#include "act/import/fetch_charges.h"
#include "act/import/import.h"
#include "act/molprop/molprop_xml.h"
#include "act/utility/stringutil.h"

namespace alexandria
{

//! \brief Description of the compound reader for the command line help
static std::vector<const char *> crDesc = {
    "It is highly recommended to provide the [TT]-charges[tt] option",
    "with a molprop file that will be used to generate charges for the",
    "file specified with the [TT]-f[tt] option.",
    "This file can also be used to read charges from, if present",
    "with the [TT]-qqm[tt] flag.[PAR]", 
    "Alternatively, you can use the [TT]-db 'molecule(s)'[tt]",
    "option to extract one or more compounds from the molprop file.[PAR]",
    "A selection of molecules into a training set and a test set (or ignore set)",
    "can be made using option [TT]-sel[tt]. The format of this file is:[PAR]",
    "iupac|Train[PAR]",
    "iupac|Test[PAR]",
    "iupac|Ignore[PAR]",
    "and you should ideally have a line for each molecule in the molecule database",
    "([TT]-mp[tt] option). Missing molecules will be ignored. Selection files",
    "can be generated using the [TT]molselect[tt] script."
};

void CompoundReader::addOptions(std::vector<t_pargs>      *pargs,
                                std::vector<t_filenm>     *filenm,
                                std::vector<const char *> *desc)
{
    std::vector<t_filenm>     fnm = {
        { efXML, "-charges", "molprop",    ffOPTRD },
        { efDAT, "-sel",  "molselect",     ffOPTRD }
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
          "Molecular structure file in e.g. pdb format. Is ignored with train_ff." },
        { "-db",     FALSE, etSTR,  {&dbname_},
          "Read one or more molecules from the database rather than from a file. To specify multiple molecules please use quotes, e.g. [TT]-db[tt] 'water methane ammonia'. Is ignored with train_ff." },
        { "-oneH", FALSE, etBOOL, {&oneH_},
          "Map all different hydrogen atom types back to H, mainly for debugging." },
        { "-qtot",   FALSE, etREAL, {&qtot_},
          "Combined charge of the molecule(s). This will be taken from the input file by default, but that is not always reliable. If the -qcustom flag is supplied, that will be used instead." },
        { "-qalg",   FALSE, etSTR,  {&qalgorithm_},
          "Algorithm to generate charges. Can be either of NONE (read from force field), EEM or SQE (depending on force field parameter availability), ESP (fit to electrostatic potential which then must be present in the molprop file, -charges option). Flag can be Custom (see flag -qcustom below) or Read (see flag -qqm below). When fitting parameters to the ESP, do not provide this flag and not the charges file either. Sorry for the confusion." },
        { "-qqm",    FALSE, etSTR,  {&qqm_},
          "Use a method from quantum mechanics that needs to be present in the input file. For instance, qESP, qHirshfeld, qCM5 or qMulliken may be available but check your input." },
        { "-qcustom", FALSE, etSTR, {&qcustom_}, 
          "Here a quoted string of custom charges can be provided such that a third party source can be used. It is then possible to generate multipoles and compare the ESP to a quantum chemistry result. The number of charges provided must match the number of particles (including shells if present in the force field used)." }
    };
    for(const auto &mp : mypargs)
    {
        pargs->push_back(mp);
    }
}

void CompoundReader::optionsFinished(MsgHandler                  *msghandler,
                                     const std::vector<t_filenm> &filenm)
{
    // molselect?
    auto mfn = opt2fn_null("-sel", filenm.size(), filenm.data());
    if (mfn != nullptr)
    {
        if (strlen(filename_) != 0 || strlen(dbname_) != 0)
        {
            msghandler->msg(ACTStatus::Error,
                            "Please provide just one of the flags -f and -db and -sel.");
            return;
        }
        molselect_.read(mfn);
    }
    else if (strlen(filename_) != 0 && strlen(dbname_) != 0)
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
    if (strlen(qalgorithm_) > 0)
    {
        qAlgorithm_ = nameToChargeGenerationAlgorithm(qalgorithm_);
    }
    else
    {
        if (strlen(qcustom_) > 0)
        {
            qAlgorithm_ = ChargeGenerationAlgorithm::Custom;
        }
        else if (strlen(qqm_) > 0)
        {
            qAlgorithm_ = ChargeGenerationAlgorithm::Read;
        }
    }
    if (strlen(dbname_) > 0 && qmapfn_.empty())
    {
        msghandler->msg(ACTStatus::Error,
                        "Please provide -charges flag in conjunction with -db.");
    }
    else if (strlen(qcustom_) > 0 && strlen(qqm_) > 0)
    {
        msghandler->msg(ACTStatus::Error,
                        "Please do not provide both the custom charges and the QM charge type to read.");
    }
    else if (strlen(qalgorithm_) > 0 && 
             ((strlen(qcustom_) > 0 && qAlgorithm_ != ChargeGenerationAlgorithm::Custom) || 
              (strlen(qqm_) > 0 && qAlgorithm_ != ChargeGenerationAlgorithm::Custom)))
    {
        msghandler->msg(ACTStatus::Warning,
                        "Please do not provide both a charge algorithm and either custom charges or the QM charge type to read.");
    }
}

void CompoundReader::setCharges(MsgHandler          *msghandler,
                                ForceField          &pd,
                                ACTMol              *mol,
                                const ForceComputer *forceComp,
                                bool                 warnQtot)
{
    if (mol->totalCharge() != qtot_ && warnQtot)
    {
        msghandler->msg(ACTStatus::Warning,
                        gmx::formatString("Detected total charge %d, command line says %g.",
                                          mol->totalCharge(), qtot_));
    }
    
    if (msghandler->ok())
    {
        std::vector<gmx::RVec> coords = mol->xOriginal();
        std::vector<gmx::RVec> forces(mol->atomsConst().size());
        if (qAlgorithm_ == ChargeGenerationAlgorithm::Custom)
        {
            std::vector<double> qcustom;
            auto mycharges = gmx::splitString(qcustom_);
            for(auto &q : mycharges)
            {
                qcustom.push_back(my_atof(q.c_str(), "custom q"));
            }
            mol->setCharges(qcustom);
        }
        else if (qAlgorithm_ == ChargeGenerationAlgorithm::NONE)
        {
            mol->setCharges(msghandler, &pd);
        }
        else if (strlen(filename_) > 0 && qmap_.empty())
        {
            mol->generateCharges(msghandler, &pd, forceComp, qAlgorithm_, &coords, &forces,
                                 msghandler->verbose());
        }
        else
        {
            mol->setCharges(msghandler, qmap_);
        }
        mol->minimizeShells(msghandler, &pd, forceComp, &coords, &forces);
        if (msghandler->debug())
        {
            auto msg = gmx::formatString("CHARGES for %s", mol->getMolname().c_str());
            auto ac = mol->atomsConst();
            for(size_t i = 0; i < ac.size(); i++)
            {
                msg += gmx::formatString(" %s-%zu %g", ac[i].name().c_str(), i+1, ac[i].charge());
            }
            msghandler->writeDebug(msg);
        }
    }
}

void CompoundReader::readFile(MsgHandler          *msghandler,
                              ForceField          &pd,
                              std::vector<ACTMol> *mols)
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
    importFile(msghandler, &pd, filename_, &mps, "",
               JobType::OPT, userQtot(), &qtot_babel,
               oneH_);
    if (!msghandler->ok())
    {
        msghandler->msg(ACTStatus::Error,
                        gmx::formatString("Reading %s failed.\n", filename_));
        return;
    }
    if (mps.size() == 0)
    {
        msghandler->msg(ACTStatus::Error,
                        gmx::formatString("Failed to import coordinate file %s using OpenBabel", filename_));
        return;
    }
    for(auto &mp : mps)
    {
        ACTMol amol;
        amol.Merge(std::move(&mp));
        mols->push_back(std::move(amol));
    }
}

void CompoundReader::read(MsgHandler          *msghandler,
                          ForceField          &pd,
                          const ForceComputer *forceComp,
                          std::vector<ACTMol> *mols)
{
    bool                  readCoordinates = strlen(filename_) > 0;
    std::set<std::string> lookup;
    std::string           lookupSource;
    if (qmapfn_.empty() && !readCoordinates)
    {
        return;
    }
    if (molselect_.nMol() > 0)
    {
        for(const auto &ims: molselect_.imolSelect())
        {
            const auto cccs = split(ims.iupac(), '#');
            for(const auto &ccc : cccs)
            {
                lookup.insert(ccc);
            }
        }
        // For printing further down
        lookupSource.assign("selection file");
    }
    else
    {
        // Try reading from a file first
        if (readCoordinates)
        {
            readFile(msghandler, pd, mols);
            if (msghandler->ok())
            {
                for(auto mol = mols->begin(); mol < mols->end(); ++mol)
                {
                    auto fp = mol->fragmentPtr();
                    if (fp)
                    {
                        for(auto ic = fp->begin(); ic < fp->end(); ++ic)
                        {
                            if (!ic->iupac().empty())
                            {
                                lookup.insert(ic->iupac());
                            }
                            else
                            {
                                lookup.insert(ic->inchi());
                            }
                        }
                    }
                    else
                    {
                        lookup.insert(mol->getMolname());
                    }
                }
            }
            // For printing further down
            lookupSource.assign(" from file ");
            lookupSource += filename_;
        }
        else if (strlen(dbname_) > 0)
        {
            // Determine what compounds the usesr selected
            for(const auto &mymol : split(dbname_, ' '))
            {
                lookup.insert(mymol);
            }
            // For printing further down
            lookupSource.assign(" from the command line flag -db");
        }
    }
    if (lookup.empty())
    {
        if (!qmapfn_.empty())
        {
            msghandler->msg(ACTStatus::Info,
                            gmx::formatString("CompoundReader will include all compounds from %s\n.", qmapfn_.c_str()));
        }
        else
        {
            // No lookup and no file.
            return;
        }
    }
    else
    {
        std::string msg("CompoundReader found the following compounds:");
        for(const auto &lu : lookup)
        {
            msg += gmx::formatString(" '%s'", lu.c_str());
        }
        msg += lookupSource;
        msghandler->msg(ACTStatus::Info, msg);
    }
    if (!qmapfn_.empty())
    {
        if (readCoordinates)
        {
            // If we read molecules from an input file, we do not want to overwrite them
            std::vector<ACTMol> temp_mols;
            qmap_ = fetchChargeMap(msghandler, &pd, forceComp, qmapfn_.c_str(),
                                   &temp_mols, lookup, qAlgorithm_, qqm_);
        }
        else
        {
            qmap_ = fetchChargeMap(msghandler, &pd, forceComp, qmapfn_.c_str(),
                                   mols, lookup, qAlgorithm_, qqm_);
        }
        msghandler->msg(ACTStatus::Info,
                        gmx::formatString("CompoundReader read %lu out of %lu entries into charge map from %s\n",
                                          qmap_.size(), lookup.size(), qmapfn_.c_str()));
    }
    // Did we find any molecule?
    if (mols->empty())
    {
        msghandler->msg(ACTStatus::Error, "Couldn't find any molecule!");
    }
    else
    {
        // Only warn about the total charge for monomers
        // since we cannot know how charge is divided between molecules.
        bool warnQtot = mols->size() == 1;
        for(auto mol = mols->begin(); mol != mols->end(); )
        {
            mol->GenerateTopology(msghandler, &pd, missingParameters::Error);
            if (msghandler->ok())
            {
                // Load all the other properties of the compounds
                mol->getExpProps(msghandler, &pd, {{ MolPropObservable::POTENTIAL, iqmType::QM }});
            }
            if (msghandler->ok())
            {
                setCharges(msghandler, pd, &(*mol), forceComp, warnQtot);
            }
            if (!msghandler->ok())
            {
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("CompoundReader could not determine charges for '%s' from '%s'\n",
                                                  mol->getMolname().c_str(), filename_));
                // Prevent false positives, delete compound and reset status
                msghandler->resetStatus();
                mol = mols->erase(mol);
            }
            else
            {
                mol++;
            }
        }
    }
}

}
