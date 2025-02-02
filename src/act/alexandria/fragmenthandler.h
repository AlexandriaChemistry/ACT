/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2025
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
#ifndef ACT_ALEXANDRIA_FRAGMENTHANDLER_H
#define ACT_ALEXANDRIA_FRAGMENTHANDLER_H
#include <vector>

#include "act/alexandria/actmol_low.h"
#include "act/alexandria/fetch_charges.h"
#include "act/qgen/qgen_acm.h"
#include "act/molprop/fragment.h"

namespace alexandria
{

    class Topology;
    
    /*! Class to control charge generation for fragments
     */
    class FragmentHandler
    {
    private:
        //! We need one ACM structure for each fragment
        std::vector<QgenAcm>               QgenAcm_;
        //! What algorithm do we use for generating charges
        ChargeGenerationAlgorithm          algorithm_ = ChargeGenerationAlgorithm::EEM;
        //! A complete topology for each fragment is needed to compute energies
        std::vector<Topology *>            topologies_;
        //! And a vector of bonds
        std::vector<std::vector<Bond> >    bonds_;
        //! Fragment InChi identifiers
        std::vector<std::string>           ids_;
        //! Array denoting where the atoms start in the global system
        std::vector<size_t>                atomStart_;
        //! Pointer to copy of the fragments
        std::vector<double>                qtotal_;
        //! Whether charges are fixed or not. They are when set from a charge map.
        bool                               fixedQ_ = false;
        //! Total number of atoms
        size_t                             natoms_ = 0;
    public:
        /*! Constructor
         * \param[in] msghandler    Message handler
         * \param[in] pd            Force field data
         * \param[in] coordinates   The atomic coordinates, may be extended if shells are present.
         * \param[in] atoms         The ActAtoms
         * \param[in] bonds         The bonds
         * \param[in] fragments     The fragmentation information
         * \param[in] missing       How to deal with missing parameters
         */
        FragmentHandler(MsgHandler                   *msghandler,
                        ForceField                   *pd,
                        const std::vector<gmx::RVec> &coordinates,
                        const std::vector<ActAtom>   &atoms,
                        const std::vector<Bond>      &bonds,
                        const std::vector<Fragment>  *fragments,
                        missingParameters             missing);

        /*! \brief Fetch charges for all atoms
         * \param[out] qq Vector that will be reinitialized at correct length
         */
        void fetchCharges(std::vector<double> *qq);
        
        //! Return the fragment ids
        const std::vector<std::string> &ids() const { return ids_; }

        //! \return the atomStart_ vector, containing one extra index more than the number of fragments.
        const std::vector<size_t> atomStart() const { return atomStart_; }

        //! \return the vector of Topology structures
        const std::vector<Topology *> topologies() const { return topologies_; }

        //! \return the vector of Topology structures
        std::vector<Topology *> &topologiesPtr() { return topologies_; }

        /*! \brief Generate charges for all fragments
         * \param[in]  fp      Debug file pointer, may be nullptr
         * \param[in]  molname Molecule name for printing
         * \param[in]  x       Atomic coordinates
         * \param[in]  pd      Force field file
         * \param[out] atoms   Atoms to store the charges in
         * \param[in]  symmetric_charges Information on which charges should be symmetrized according to FF
         * \return Error code or OK if all is fine.
         */
        eQgen generateCharges(FILE                         *fp,
                              const std::string            &molname,
                              const std::vector<gmx::RVec> &x,
                              const ForceField             *pd,
                              std::vector<ActAtom>         *atoms,
                              const std::vector<int>       &symmetric_charges);

        /*! \brief Copy charges from the atoms to the fragments
         * \param[in] atoms The atoms from the complete topology for all fragments
         */
        void setCharges(const std::vector<ActAtom> &atoms);
        /*! \brief Copy charges from fragments back to atoms
         * \param[inout] atoms The atoms to copy to.
         * \return whether atoms->size() equals the sum of the number of atoms in the fragments.
         */
        bool fetchCharges(std::vector<ActAtom> *atoms);
        /*! \brief Copy charges from a charge map
         * \param[in] q The charges
         */
        void setCharges(const std::vector<double> &q);
        /*! \brief Copy charges from a charge map. Check msghandler status on return
         * \param[in] msghandler For status and messaging
         * \param[in] qmap       The charge map
         */
        void setCharges(MsgHandler      *msghandler,
                        const chargeMap &qmap);
        //! \return whether charges are fixed or not
        bool fixedCharges() const { return fixedQ_; }
        /*! \brief Set the charge generation algorithm to use
         * \param[in] alg The algorithm to use. Only Read or EEM/SQE are supported.
         */
        void setChargeGenerationAlgorithm(ChargeGenerationAlgorithm alg) { algorithm_ = alg; }
    };
} // namespace alexandria

#endif
