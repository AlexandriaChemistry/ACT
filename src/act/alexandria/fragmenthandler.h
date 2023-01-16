/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
#include <vector>

#include "act/qgen/qgen_acm.h"
#include "act/molprop/fragment.h"

namespace alexandria
{

    /*! Class to control charge generation for fragments
     */
    class FragmentHandler
    {
    private:
        //! We need one ACM structure for each fragment
        std::vector<QgenAcm>               QgenAcm_;
        //! A complete topology for each fragment is needed to compute energies
        std::vector<Topology>              topologies_;
        //! And a vector of bonds
        std::vector<std::vector<Bond> >    bonds_;
        //! Fragment identifiers
        std::vector<std::string>           ids_;
        //! Array denoting where the atoms start in the global
        std::vector<size_t>                atomStart_;
        //! Total number of atoms
        size_t                             natoms_ = 0;
    public:
        /*! Constructor
         * \param[in] pd            Force field data
         * \param[in] coordinates   The atomic coordinates
         * \param[in] atoms         The atoms
         * \param[in] bonds         The bonds
         * \param[in] fragments     The fragmentation information
         * \param[in] shellRenumber Info on renumbering atoms because of shells
         * \param[in] missing       How to deal with missing parameters
         */
        FragmentHandler(const Poldata                *pd,
                        const std::vector<gmx::RVec> &coordinates,
                        const std::vector<ActAtom>   &atoms,
                        const std::vector<Bond>      &bonds,
                        const std::vector<Fragment>  *fragments,
                        const std::vector<int>       &shellRenumber,
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
        const std::vector<Topology> topologies() const { return topologies_; }

        /*! \brief Generate charges for all fragments
         * \param[in]  fp      Debug file pointer, may be nullptr
         * \param[in]  molname Molecule name for printing
         * \param[in]  x       Atomic coordinates
         * \param[in]  pd      Force field file
         * \param[out] atoms   Atoms to store the charges in
         * \return Error code or OK if all is fine.
         */
        eQgen generateCharges(FILE                         *fp,
                              const std::string            &molname,
                              const std::vector<gmx::RVec> &x,
                              const Poldata                *pd,
                              std::vector<ActAtom>         *atoms);
    };
} // namespace alexandria
