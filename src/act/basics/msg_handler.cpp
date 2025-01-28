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

#include <cstdlib>

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
    { ACTMessage::NoData,                   "No experimental data" },
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
    { ACTMessage::NotSupportedDihedral,     "NotSupportedDihedral" }
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
    
    static std::string combine(ACTMessage  actm,
                               const char *msg)
    {
        return gmx::formatString("%s: %s", ACTMessages[actm], msg);
    }

    MsgHandler::MsgHandler()
    {
        for(const auto actm: ACTMessages)
        {
            wcount_.insert( { actm.first, 0 });
        }
    }    

    void MsgHandler::fatal(ACTMessage  actm,
                           const char *msg) const
    {
        if (fp_)
        {
            fprintf(fp_, "%s\n", combine(actm, msg).c_str());
            
            fclose(fp_);
        }
        fprintf(stderr, "%s\n", combine(actm, msg).c_str());
        throw;
    }

    void MsgHandler::warning(ACTMessage  actm,
                             const std::string &msg)
    {
        if (verbose_)
        {
            if (fp_)
            {
                fprintf(fp_, "%s\n", combine(actm, msg.c_str()).c_str());
            }
            else
            {
                fprintf(stdout, "%s\n", combine(actm, msg.c_str()).c_str());
            }
        }
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
