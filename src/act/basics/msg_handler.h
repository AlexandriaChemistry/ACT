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

namespace alexandria
{
enum class ACTMessage
    {
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
        NotSupportedDihedral
    };

    extern std::map<ACTMessage, const char *> ACTMessages;

    /*! \brief Return error message corresponding to code
     * \param[in] actm The code
     * \return The corresponding message
     */
    const char *actMessage(ACTMessage actm);

/*! \brief Simple class to print message to a file or stdout
 * Verbosity level can be set.
 * Messages are counted and can be summarized
 */
class MsgHandler
    {
    private:
        //! File pointer to write stuff to
        FILE         *fp_      = nullptr;
        //! Verbosity level
        bool          verbose_ = false;
        //! Warning count per type
        std::map<ACTMessage, unsigned int> wcount_;
    public:
        //! Constructor
        MsgHandler();

        /*! \brief Set the file pointer
         * \param[in] fp the file pointer
         */
        void setFilePointer(FILE *fp)
        {
            fp_ = fp;
        }
        /*! Set verbosity level
         * \param[in] verbose Whether or not to print a lot
         */
        void setVerbosity(bool verbose) { verbose_ = verbose; }
        //! \return the current verbosity level
        bool verbose() const { return verbose_; }
        /*! \brief Fatal error message, will throw a fatal error
         * \param[in] actm The message type
         * \param[in] msg  Additional information to provide
         */
        void fatal(ACTMessage  actm,
                   const char *msg) const;
        /*! \brief Fatal error message, will throw a fatal error
         * \param[in] actm The message type
         * \param[in] msg  Additional information to provide
         */
        void fatal(ACTMessage         actm,
                   const std::string &msg) const { fatal(actm, msg.c_str()); }
        /*! \brief Warning message, will only print if verbose, but type will be logged
         * \param[in] actm The message type
         * \param[in] msg  Additional information to provide
         */
        void warning(ACTMessage         actm,
                     const std::string &msg);
        //! \brief Print summary of warnings to file pointer or stdout, independent of verbosity
        void summary() const;
        //! \return warning count
        unsigned int warningCount(ACTMessage actm) { return wcount_[actm]; }
    };

} // namespace

#endif
