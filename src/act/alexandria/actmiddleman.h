/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */
#ifndef ACT_ACTMIDDLEMAN_H
#define ACT_ACTMIDDLEMAN_H

#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/ga/mutator.h"
#include "gromacs/fileio/oenv.h"

namespace gmx
{
class TextWriter;
}
    
namespace alexandria
{
    // Declare some classes rather than including headers.
    class StaticIndividualInfo;
    class ACMFitnessComputer;
    class ACMIndividual;
    class MolGen;
   
    /*! \brief Class that manages the deviation calculations
     */
    class ACTMiddleMan
    {
    private:
        //! Fitness computer
        ACMFitnessComputer    fitComp_;
        //! Mutator
        ga::Mutator          *mutator_ = nullptr;
        //! Individual
        ACMIndividual        *ind_     = nullptr;
        //! Config handler for GA
        GAConfigHandler      *gach_    = nullptr;
        //! Force computer
        ForceComputer         forceComp_;
        //! My ID
        int                   id_;

        //! \brief Stop my helpers
        void stopHelpers();

    public:
        /*! \brief Constructor
         * \param[in] msghandler Message and status handler
         * \param[in] mg         Molecule info
         * \param[in] sii        The individual info
         * \param[in] gach       GA Config handler
         * \param[in] bch        Bayes Config handler
         * \param[in] oenv       GROMACS output environment
         * \param[in] algorithm  The charge generation algorithm
         * \param[in] openConvFiles Whether or not to create convergence files
         */
        ACTMiddleMan(MsgHandler                *msghandler,
                     MolGen                    *mg,
                     StaticIndividualInfo      *sii,
                     GAConfigHandler           *gach,
                     BayesConfigHandler        *bch,
                     gmx_output_env_t          *oenv,
                     ChargeGenerationAlgorithm  algorithm,
                     bool                       openConvFiles);

        //! \brief Destructor
        ~ACTMiddleMan();

        /*! \brief Run the helper process
         * \param[in] msghandler Message and status handler
         */
        void run(MsgHandler *msghandler);

        //! Print the MCMC statistics if appropriate
        void printStatistics(gmx::TextWriter *tw);

    };
    
} // namespace alexandria

#endif
