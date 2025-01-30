/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2025
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
#include "actpre.h" 

#include "msg_handler.h"

#include <cstdio>
#include <cstdlib>

#include "gromacs/commandline/filenm.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

std::map<ACTMessage, const char *> ACTMessages = {
    { ACTMessage::Unknown,                  "Unknown status" },
    { ACTMessage::OK,                       "OK" },
    { ACTMessage::NoAtoms,                  "No Atoms" },
    { ACTMessage::ZeroDip,                  "Zero Dipole" },
    { ACTMessage::NoQuad,                   "No Quadrupole" },
    { ACTMessage::Charged,                  "Charged" },
    { ACTMessage::AtomTypes,                "Atom type problem" },
    { ACTMessage::AtomNumber,               "Atom number problem" },
    { ACTMessage::MolpropConv,              "Converting from molprop" },
    { ACTMessage::Multiplicity,             "Number of electrons does not match the multiplicity. Is the total charge correct?" },
    { ACTMessage::BondOrder,                "Determining bond order" },
    { ACTMessage::RespInit,                 "RESP Initialization" },
    { ACTMessage::ChargeGeneration,         "Charge generation" },
    { ACTMessage::MissingChargeGenerationParameters, "Parameters for charge generation missing" },
    { ACTMessage::ShellMinimization,        "Shell minimization" },
    { ACTMessage::Topology,                 "No input to generate a topology" },
    { ACTMessage::FragmentHandler,          "Fragment Handler could not make topologies" },
    { ACTMessage::QMInconsistency,          "QM Inconsistency (ESP dipole does not match Electronic)" },
    { ACTMessage::Test,                     "Compound not in training set" },
    { ACTMessage::NoData,                   "No data" },
    { ACTMessage::NoMolpropCharges,         "No charges in the molprop file" },
    { ACTMessage::GenShells,                "Generating shells" },
    { ACTMessage::GenBonds,                 "Generating bonds" },
    { ACTMessage::CommProblem,              "Communicating MolProp" },
    { ACTMessage::ZeroZeta,                 "Charge distribution width zeta is zero unexpectedly" },
    { ACTMessage::InsufficientDATA,         "The number of data is lower than mindata" },
    { ACTMessage::NoDipole,                 "No Dipole moment" },
    { ACTMessage::NotSupportedBond,         "NotSupportedBond" },
    { ACTMessage::NotSupportedAngle,        "NotSupportedAngle" },
    { ACTMessage::NotSupportedLinearAngle,  "NotSupportedLinearAngle" },
    { ACTMessage::NotSupportedDihedral,     "NotSupportedDihedral" },
    { ACTMessage::MissingFFParameter,       "Missing parameter in force field" },
    { ACTMessage::MinimizationFailed,       "Minimization failed" },
    { ACTMessage::Info,                     "Information" }
};

std::map<ACTStatus, const char *> statnm = {
    { ACTStatus::Fatal, "Fatal" },
    { ACTStatus::Error, "Error" },
    { ACTStatus::Warning, "Warning" },
    { ACTStatus::Verbose, "Verbose" },
    { ACTStatus::Debug, "Debug" }
};

const char *actMessage(ACTMessage actm)
{
    auto m = ACTMessages.find(actm);
    if (ACTMessages.end() == m)
    {
        return ACTMessages[ACTMessage::Unknown];
    }
    else
    {
        return m->second;
    }
}
    
void MsgHandler::print(ACTStatus   level,
                       ACTMessage  actm,
                       const char *msg) const
{
    auto mymsg = gmx::formatString("%s - %s: %s", 
                                   statnm[level],
                                   ACTMessages[actm], msg);
    if (level <= printLevel_)
    {
        if (fp_)
        {
            fprintf(fp_, "%s\n", mymsg.c_str());
            if (flush_)
            {
                fflush(fp_);
            }
        }
        else
        {
            if (level <= ACTStatus::Error)
            {
                fprintf(stderr, "%s\n", mymsg.c_str());
                if (flush_)
                {
                    fflush(stderr);
                }
            }
            else
            {
                fprintf(stdout, "%s\n", mymsg.c_str());
                if (flush_)
                {
                    fflush(stdout);
                }
            }
        }
    }
}

MsgHandler::MsgHandler()
{
    for(const auto actm: ACTMessages)
    {
        wcount_.insert( { actm.first, 0 });
    }
}

MsgHandler::~MsgHandler()
{
    if (!filename_.empty() && fp_)
    {
        gmx_ffclose(fp_);
    }
}

void MsgHandler::addOptions(std::vector<t_pargs>      *pargs,
                            std::vector<t_filenm>     *filenm,
                            const std::string         &defaultLogName)
{
    pargs->push_back( { "-v", FALSE, etINT, { &ilevel_ },
            "Verbosity level: 0 (Fatal), 1 (Error), 2 (Warning), 3 (Info), 4 (Debug)" } );
    pargs->push_back( { "-flush", FALSE, etBOOL, {&flush_},
            "Flush output immediately rather than letting the OS buffer it. Don't use for production simulations."} );

    filenm->push_back( { efLOG, "-g", defaultLogName.c_str(), ffWRITE }); 
}

void MsgHandler::optionsFinished(const std::vector<t_filenm> &filenm)
{
    if (ilevel_ == 0)
    {
        printLevel_ = ACTStatus::Fatal;
    }
    else if (ilevel_ == 1)
    {
        printLevel_ = ACTStatus::Error;
    }
    else if (ilevel_ == 2)
    {
        printLevel_ = ACTStatus::Warning;
    }
    else if (ilevel_ == 3)
    {
        printLevel_ = ACTStatus::Verbose;
    }
    else
    {
        printLevel_ = ACTStatus::Debug;
    }
    std::string fn(opt2fn("-g", filenm.size(), filenm.data()));
    setFileName(fn);
}

void MsgHandler::setFileName(const std::string &fn)
{
    fp_ = gmx_ffopen(fn.c_str(), "w");
    if (!fp_)
    {
        fatal(ACTMessage::Info, gmx::formatString("Cannot open file '%s' for writing", fn.c_str()).c_str());
    }
}

void MsgHandler::fatal(ACTMessage  actm,
                       const char *msg)
{
    status_ = ACTStatus::Fatal;
    print(status_, actm, msg);
    if (fp_)
    {
        fclose(fp_);
    }
    GMX_THROW(gmx::InvalidInputError("See message above"));
}

void MsgHandler::msg(ACTStatus         level,
                     ACTMessage        actm,
                     const std::string &msg)
{
    if (level == ACTStatus::Fatal)
    {
        fatal(actm, msg.c_str());
    }
    if (level < status_)
    {
        status_ = level;
    }
    last_ = actm;
    print(level, actm, msg.c_str());
    wcount_.find(actm)->second++;
}

void MsgHandler::summary() const
{
    for(const auto &wc : wcount_)
    {
        if (fp_)
        {
            fprintf(fp_, "%s : %d\n", ACTMessages[wc.first], wc.second);
        }
    }
}

}
