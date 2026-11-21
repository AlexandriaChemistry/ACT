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
#ifndef COMPOUND_READER_H
#define COMPOUND_READER_H

    
#include <cstring>
#include <string>
#include <vector>

#include "act/alexandria/molselect.h"
#include "act/basics/chargemodel.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/import/chargemap.h"
#include "act/topology/actmol.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/real.h"

namespace alexandria
{
    /*! Class to read compounds from files and turn them into full-blown
     * ACTmol objects by adding charges and topologies.
     */
    class CompoundReader
    {
    private:
        /*! \brief The text options as the command line parser sets them
         *
         * They point into the command line, and are only read by
         * optionsFinished, which copies them. All other code uses the
         * strings that the reader owns, such that it can be copied.
         */
        struct ParsedText
        {
            //! -qqm
            const char *qqm        = "";
            //! -qcustom
            const char *qcustom    = "";
            //! -qalg
            const char *qalgorithm = "";
            //! -f
            const char *filename   = "";
            //! -smiles
            const char *smiles     = "";
            //! -db
            const char *dbname     = "";
        };
        //! What the parser set, see ParsedText
        ParsedText  parsed_;
        //! Total charge of all compounds
        real        qtot_       = 0;
        //! Charge method to select from QM methods available in the molprop file
        std::string qqm_;
        //! String of custom charges
        std::string qcustom_;
        //! How to generate charges, either from the force field, by reading, or fitting ESP
        ChargeGenerationAlgorithm qAlgorithm_ = ChargeGenerationAlgorithm::SQE;
        //! File name for reading structures from
        std::string filename_;
        //! SMILES string to make a structure from
        std::string smiles_;
        //! List of compounds/dimers to extract from charges files
        std::string dbname_;
        //! Map back hydrogen atoms to one type
        bool        oneH_       = false;
        //! Name of the molprop file that contains charges or reference structures for generating charges
        std::string qmapfn_;
        //! Data structure for charges
        ChargeMap   qmap_;
        //! MolSelect structre if data was provided by the user
        MolSelect   molselect_;

    public:
        //! Constructor
        CompoundReader() = default;

        /*! \brief Add command-line options
         * \param[inout] pargs  List of command line flags
         * \param[inout] filenm List of file options
         * \param[inout] desc   Additional help text
         */
        void addOptions(std::vector<t_pargs>      *pargs,
                        std::vector<t_filenm>     *filenm,
                        std::vector<const char *> *desc);

        /*! Check whether passed options make sense
         * \param[in] msghandler Error and message handler
         * \param[in] filenm The filenames after processing
         */
        void optionsFinished(MsgHandler                  *msghandler,
                             const std::vector<t_filenm> &filenm);

        /*! \brief Read the structure from a file, as the option -f does
         * For a program that uses this class without a command line. The
         * SMILES string, if any, is forgotten.
         * \param[in] fileName The file
         */
        void setFileName(const std::string &fileName);

        /*! \brief Make the structure from a SMILES string, as the option -smiles does
         * For a program that uses this class without a command line. The
         * file name, if any, is forgotten.
         * \param[in] smiles The SMILES string
         */
        void setSmiles(const std::string &smiles);

        /*! \brief Set the total charge of the compound, as the option -qtot does
         * \param[in] qtot The charge
         */
        void setTotalCharge(double qtot) { qtot_ = qtot; }

        //! \return whether there is one H only
        bool oneH() const { return oneH_; }

        //! Set the oneH variable
        void setOneH(bool oneH) { oneH_ = oneH; }

        //! \return the charge to read
        const char *qread() const { return qqm_.c_str(); }

        //! \return the charge algorithm
        ChargeGenerationAlgorithm algorithm() const { return qAlgorithm_; }

        /*! \brief Set the charge algorithm
         * \param[in] qA The new algorithm
         */
        void setAlgorithm(ChargeGenerationAlgorithm qA) { qAlgorithm_ = qA; }

        //! \return charge map
        const ChargeMap &chargeMapConst() const { return qmap_; }

        //! \return the file of the monomers for the charges (-charges), empty if none
        const std::string &chargeFile() const { return qmapfn_; }

        //! \return charge map for editing
        ChargeMap *chargeMap() { return &qmap_; }

        //! \return the SMILES string to make a molecule from (-smiles), empty if none
        const char *smiles() const { return smiles_.c_str(); }

        //! \return the file to read the structure from (-f), empty if none
        const char *fileName() const { return filename_.c_str(); }

        //! \return the compound to take from the charge file (-db), empty if none
        const char *dbName() const { return dbname_.c_str(); }

        //! \return the custom charges (-qcustom), empty if none
        const char *customCharges() const { return qcustom_.c_str(); }

        //! \return the total charge of the compounds (-qtot)
        double totalCharge() const { return qtot_; }

        /*! Read molecule(s) from a single file, or make one from the SMILES string
         * \param[in]  msghandler Error and message handler, check whether ok after returning
         * \param[in]  pd         The force field
         * \param[out] mols       The molecule
         */
        void readFile(MsgHandler          *msghandler,
                      ForceField          &pd,
                      std::vector<ACTMol> *mols);
        /*! Set charges for a single molecule
         * \param[in]  msghandler Error and message handler, check whether ok after returning
         * \param[in]  pd         The force field
         * \param[out] mol        The molecule
         * \param[in]  forceComp  A force computer
         * \param[in]  warnQtot   Print a warning when qtot does not match the input
         */
        void setCharges(MsgHandler          *msghandler,
                        ForceField          &pd,
                        ACTMol              *mol,
                        const ForceComputer *forceComp,
                        bool                 warnQtot);
        /*! \brief Do the actual reading and processing
         * \param[in]  msghandler Error and message handler
         * \param[in]  pd         The force field
         * \param[in]  forceComp  A force computer
         * \param[out] mols       The ACTMol structures
         */
        void read(MsgHandler          *msghandler,
                  ForceField          &pd,
                  const ForceComputer *forceComp,
                  std::vector<ACTMol> *mols);

        //! \return the molselect structure (may be empty!)
        const MolSelect &molselect() const { return molselect_; }

        //! \return true if user provided charges
        bool userQtot() const { return !qcustom_.empty(); }
    };
    
}
#endif
