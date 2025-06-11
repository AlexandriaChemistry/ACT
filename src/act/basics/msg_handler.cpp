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

#include "act/basics/version.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

std::map<ACTMessage, const char *> ACTMessages = {
    { ACTMessage::Silent,                   nullptr },
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
    { ACTMessage::MinimizationFailed,       "Minimization failed" }
};

std::map<ACTStatus, const char *> statnm = {
    { ACTStatus::Fatal, "Fatal" },
    { ACTStatus::Error, "Error" },
    { ACTStatus::Warning, "Warning" },
    { ACTStatus::Info, "Info" },
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
    std::string mymsg;

    if (actm == ACTMessage::Silent)
    {
        mymsg = gmx::formatString("%s: %s", statnm[level], msg);
    }
    else
    {
        mymsg = gmx::formatString("%s - %s: %s", 
                                  statnm[level],
                                  ACTMessages[actm], msg);
    }
    if (level <= printLevel_)
    {
        if (level == ACTStatus::Debug)
        {
            writeDebug(mymsg);
        }
        else
        {
            write(mymsg);
        }
        if (level <= ACTStatus::Error)
        {
            fprintf(stderr, "%s\n", mymsg.c_str());
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
    if (!filename_.empty())
    {
        printf("\nPlease check output in file %s.\n\n", filename_.c_str());
    }
    if (printLevel_ == ACTStatus::Debug)
    {
        printf("\nPlease check debug statements in the actXXXXX.debug files.\n\n");
    }
    if (tw_)
    {
        tw_->writeLine("Program finished, please check the output in this file.");
        tw_->writeLine();
        write(act_goodbye());
        tw_->close();
        delete tw_;
    }
    if (twdebug_)
    {
        twdebug_->close();
        delete twdebug_;
    }
}

void MsgHandler::flush()
{
    if (tw_)
    {
        tw_->flush();
    }
}

void MsgHandler::addOptions(std::vector<t_pargs>      *pargs,
                            std::vector<t_filenm>     *filenm,
                            const std::string         &defaultLogName)
{
    defaultLogName_ = defaultLogName;
    pargs->push_back( { "-v", FALSE, etINT, { &ilevel_ },
            "Verbosity level: 0 (Fatal), 1 (Error), 2 (Warning), 3 (Info), 4 (Verbose), 5 (Debug)" } );

    filenm->push_back( { efLOG, "-g", defaultLogName_.c_str(), ffWRITE }); 
}

void MsgHandler::optionsFinished(const std::vector<t_filenm> &filenm,
                                 const CommunicationRecord   *cr)
{
    std::map<int, ACTStatus> i2s = {
        { 0, ACTStatus::Fatal },
        { 1, ACTStatus::Error },
        { 2, ACTStatus::Warning },
        { 3, ACTStatus::Info },
        { 4, ACTStatus::Verbose },
        { 5, ACTStatus::Debug }
    };

    auto i2 = i2s.find(ilevel_);
    if (i2s.end() != i2)
    {
        printLevel_ = i2->second;
    }
    else
    {
        printLevel_ = ACTStatus::Info;
    }
    if (cr->isMaster())
    {
        filename_.assign(opt2fn("-g", filenm.size(), filenm.data()));
        tw_ = new gmx::TextWriter(filename_);
        tw_->writeLine(act_welcome());
        tw_->writeLineFormatted("Verbosity level %d (%s or more serious messages are printed).",
                                ilevel_, statnm[printLevel_]);
        time_t my_t;
        time(&my_t);
        tw_->writeStringFormatted("# This file was created %s", ctime(&my_t));
        tw_->writeLine();
    }
    if (printLevel_ == ACTStatus::Debug)
    {
        std::string debugfn = gmx::formatString("act%05d.debug", cr->rank());
        twdebug_ = new gmx::TextWriter(debugfn);
    }
}

void MsgHandler::fatal(ACTMessage  actm,
                       const char *msg)
{
    status_ = ACTStatus::Fatal;
    print(status_, actm, msg);
    if (tw_)
    {
        tw_->close();
        delete tw_;
        tw_ = nullptr;
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

void MsgHandler::msg(ACTStatus         level,
                     const std::string &msg)
{
    auto actm = ACTMessage::Silent;
    if (level == ACTStatus::Fatal)
    {
        fatal(actm, msg.c_str());
    }
    if (level < status_)
    {
        status_ = level;
    }
    print(level, actm, msg.c_str());
    wcount_.find(actm)->second++;
}

void MsgHandler::summary() const
{
    for(const auto &wc : wcount_)
    {
        if (tw_)
        {
            tw_->writeStringFormatted("%s : %d\n", ACTMessages[wc.first], wc.second);
        }
    }
}

}
