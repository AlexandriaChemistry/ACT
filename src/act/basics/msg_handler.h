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

#ifndef ACT_MSGHANDLER_H
#define ACT_MSGHANDLER_H

#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/textwriter.h"

namespace alexandria
{

class CommunicationRecord;

enum class ACTMessage
    {
        //! Silent, do not print anything
        Silent,
        //! Unknown error, this should not happen
        Unknown,
        //! No error, all OK.
        OK,
        //! No atoms foumd
        NoAtoms,
        //! Zero dipole found in a molecule
        ZeroDip,
        //! No Quadrupole present
        NoQuad,
        //! The compound is charged
        Charged,
        //! Cannot determine the atom types for one or more atoms
        AtomTypes,
        //! Cannot determine the atom number for one or more atoms
        AtomNumber,
        //! Incorrect multiplicity
        Multiplicity,
        //! Problem converting data from the underlying molprop structure
        MolpropConv,
        //! Incorrect or unknown bondorder
        BondOrder,
        //! Initializing the Restrained Electrostatic Potential algorithm
        RespInit,
        //! No charges in the input
        NoMolpropCharges,
        //! Problem generating charges
        ChargeGeneration,
        //! Missing charge generation parameters
        MissingChargeGenerationParameters,
        //! Shell minimization did not converge
        ShellMinimization,
        //! QM Inconsistency (ESP dipole does not match Electronic)
        QMInconsistency,
        //! No input to generate topology
        Topology,
        //! Something in FragmentHandler
        FragmentHandler,
        //! Compound not in training set
        Test,
        //! No experimental data
        NoData,
        //! Problem generating shells
        GenShells,
        //! Problem generating bonds
        GenBonds,
        //! Problem communicating MolProp between processors
        CommProblem,
        //! Charge distribution width zeta is zero unexpectedly
        ZeroZeta,
        //! The number of data is lower than mindata
        InsufficientDATA,
        //! No dipole moment
        NoDipole,
        //! Not a supported bond
        NotSupportedBond,
        //! A not supported angle was found
        NotSupportedAngle,
        //! A not supported LinearAngle was found
        NotSupportedLinearAngle,
        //! A not supported Dihedral was found
        NotSupportedDihedral,
        //! Missing FF parameter
        MissingFFParameter,
        //! Minimization failed
        MinimizationFailed
    };

extern std::map<ACTMessage, const char *> ACTMessages;

/*! \brief Return error message corresponding to code
 * \param[in] actm The code
 * \return The corresponding message
 */
const char *actMessage(ACTMessage actm);

/*! \brief enum to determine how verbose we are going to be
 * Each level includes the previous levels as well
 */
enum class ACTStatus {
    //! Only print fatal errors
    Fatal = 0,
    //! Errors that are not fatal at once
    Error = 1,
    //! Also print warnings
    Warning = 2,
    //! Just right amount of output
    Info = 3,
    //! Give even more output
    Verbose = 4,
    //! Print debugging messages as well
    Debug = 5
};

/*! \brief Simple class to print message to a file or stdout
 * ACTStatus level can be set.
 * Messages are counted and can be summarized
 */
class MsgHandler
    {
    private:
        //! TextWriter for normal output
        gmx::TextWriter *tw_          = nullptr;
        //! TextWriter for debug output
        gmx::TextWriter *twdebug_     = nullptr;
        //! File name of log file if known
        std::string      filename_;
        //! Default log file name
        std::string      defaultLogName_;
        //! ACTStatus level, messages at this or lower level are printed
        ACTStatus        printLevel_  = ACTStatus::Fatal;
        //! Level integer for command line selection of print level
        int              ilevel_      = 2;
        //! Lowest level of error encountered (lower is more severe)
        ACTStatus        status_      = ACTStatus::Debug;
        //! Last message reported
        ACTMessage       last_;
        //! Warning count per type
        std::map<ACTMessage, unsigned int> wcount_;
        /*! \brief Print at different levels of seriousness
         * \param[in] level Seriousness
         * \param[in] actm  The message type
         * \param[in] msg   Additional information to provide
         */
        void print(ACTStatus   level,
                   ACTMessage  actm,
                   const char *msg) const;

    public:
        //! Constructor
        MsgHandler();

        /*! Destructor
         * If we are responsible for opening the file, close it as well.
         */
        ~MsgHandler();

        /*! \brief Add command-line options
         * \param[in] pargs  Command line flags
         * \param[in] filenm Filenames
         * \param[in] defaultLogName Variable name says it all
         */
        void addOptions(std::vector<t_pargs>      *pargs,
                        std::vector<t_filenm>     *filenm,
                        const std::string         &defaultLogName);

        //! \brief Check and evaluate command line options
        void optionsFinished(const std::vector<t_filenm> &filenm,
                             const CommunicationRecord   *cr);

        //! \return the file name (may be empty)
        std::string filename() const { return filename_; }

        /*! Set verbosity level for printing
         * \param[in] verbose Whether or not to print a lot
         */
        void setPrintLevel(ACTStatus level) { printLevel_ = level; }

        //! \return the current verbosity level
        ACTStatus printLevel() const { return printLevel_; }

        //! \brief Reset severity level
        void resetStatus() { status_ = ACTStatus::Debug; }

        /*! \brief Fatal error message, will throw a fatal error
         * \param[in] actm The message type
         * \param[in] msg  Additional information to provide
         */
        void fatal(ACTMessage  actm,
                   const char *msg);

        /*! \brief Fatal error message, will throw a fatal error
         * \param[in] actm The message type
         * \param[in] msg  Additional information to provide
         */
        void fatal(const std::string &msg) { fatal(ACTMessage::Silent, msg.c_str()); }

        /*! \brief Message, will only print if verbosity level
         * is at least what has been configure but type will be logged
         * \param[in] level The verbosity of this message
         * \param[in] actm  The message type
         * \param[in] msg   Additional information to provide
         */
        void msg(ACTStatus          level,
                 ACTMessage         actm,
                 const std::string &msg);

        /*! \brief Simple message, will only print if verbosity level
         * is at least what has been configure but type will be logged
         * \param[in] level The verbosity of this message
         * \param[in] msg   Additional information to provide
         */
        void msg(ACTStatus          level,
                 const std::string &msg);

        //! Return pointer to internal TextWriter object
        gmx::TextWriter *tw() { return tw_; }

        //! Return pointer to internal TextWriter object for debugging
        gmx::TextWriter *twDebug() { return twdebug_; }

        //! \brief Flush the output if possible
        void flush();

        /*! \brief Just write a string
         * \param[in] str String
         */
        void write(const std::string &s) const { if (tw_) tw_->writeLine(s); }

        /*! \brief Just write a string
         * \param[in] str String
         */
        void writeDebug(const std::string &s) const { if (twdebug_) twdebug_->writeLine(s); }

        //! \return Level of severity encountered
        ACTStatus status() const { return status_; }

        //! \return whether the status is above error level 
        bool ok() const { return status_ > ACTStatus::Error; }

        //! \return whether we are in a info mode
        bool info() const { return printLevel_ >= ACTStatus::Info; }

        //! \return whether we are in a verbose mode
        bool verbose() const { return printLevel_ >= ACTStatus::Verbose; }

        //! \return whether we are in debug mode
        bool debug() const { return printLevel_ == ACTStatus::Debug; }

        //! \return ID of last message
        ACTMessage last() const { return last_; }

        /*! \brief Print summary of warnings to file pointer or stdout
         * independent of verbosity.
         */
        void summary() const;

        //! \return warning count
        unsigned int warningCount(ACTMessage actm) { return wcount_[actm]; }
    };

} // namespace

#endif
