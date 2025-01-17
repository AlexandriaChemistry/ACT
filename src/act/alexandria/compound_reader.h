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
#ifndef COMPOUND_READER_H
#define COMPOUND_READER_H

#include <cstdio>
    
#include <vector>
    
#include "act/alexandria/actmol.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"

namespace alexandria
{
    /*! Class to read compounds from files and turn them into full-blown
     * ACTmol objects by adding charges and topologies.
     */
    class CompoundReader
    {
    private:
        //! Total charge of all compounds
        double  qtot_       = 0;
        //! Charge method to select from QM methods available in the molprop file
        char   *qqm_        = (char *)"";
        //! String of custom charges
        char   *qcustom_    = (char *)"";
        //! Whether to generate charges from the force field
        bool    genCharges_ = false;
        //! File name for reading structure from
        char   *filename_   = (char *)"";
        //! Molecule name to use if there is none in the input file
        char   *molnm_      = (char *)"";
        //! List of compounds/dimers to extract from charges files
        char   *dbname_     = (char *)"";
        //! Map back hydrogen atoms to one type
        bool               oneH_                 = false;
        //! Name of the charge map file
        std::string qmapfn_;
        //! File pointer for debug messages
        FILE   *logFile_    = stderr;
        /*! Read molecule from a single file
         * \param[in] pd   The force field
         * \param[out] mol The molecule
         * \return true if successful
         */
        bool readFile(ForceField &pd,
                      ACTMol     *mol);
        /*! Set charges for a single molecule
         * \param[in]  pd        The force field
         * \param[out] mol       The molecule
         * \param[in]  qmap      A charge map
         * \param[in]  forceComp A force computer
         * \param[in]  warnQtot  Print a warning when qtot does not match the input
         * \return true if successful
         */
        bool setCharges(ForceField          &pd,
                        ACTMol              *mol,
                        const chargeMap     &qmap,
                        const ForceComputer *forceComp,
                        bool                 warnQtot);
    public:
        // Constructor
        CompoundReader() {}

        /*! \brief Add command-line options
         * \param[inout] pargs  List of command line flags
         * \param[inout] filenm List of file options
         * \param[inout] desc   Additional help text
         */
        void addOptions(std::vector<t_pargs>      *pargs,
                        std::vector<t_filenm>     *filenm,
                        std::vector<const char *> *desc);

        /*! Check whether passed options make sense
         * \param[in] filenm The filenames after processing
         * \return true if all is fine
         */
        bool optionsOK(const std::vector<t_filenm> &filenm);

        /*! \brief Set a new log file
         * \param[in] logfile The new file pointer to use
         */
        void setLogfile(FILE *logfile) { logFile_ = logfile; }

        //! \return whether there is one H only
        bool oneH() const { return oneH_; }

        /*! \brief Do the actual reading and processing
         * \param[in] pd        The force field
         * \param[in] forceComp A force computer
         * \return A vector of zero or more compounds
         */
        std::vector<ACTMol> read(ForceField          &pd,
                                 const ForceComputer *forceComp);
        
        //! Return true if user provided charges
        bool userQtot() const { return strlen(qcustom_) > 0; }
    };
    
}
#endif
