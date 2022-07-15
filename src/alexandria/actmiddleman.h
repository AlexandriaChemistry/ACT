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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#ifndef ACT_ACTMIDDLEMAN_H
#define ACT_ACTMIDDLEMAN_H

#include "gromacs/fileio/oenv.h"
#include "gromacs/mdtypes/commrec.h"

#include "act/ga//Mutator.h"

#include "confighandler.h"
#include "staticindividualinfo.h"

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
        ACMFitnessComputer   *fitComp_;
        //! Mutator
        ga::Mutator          *mutator_;
        //! Individual
        ACMIndividual        *ind_;
        //! Config handler for GA
        GAConfigHandler      *gach_;
        //! Force computer
        ForceComputer        *forceComp_;
        //! SharedIndividualInfo
        StaticIndividualInfo *sii_;
        //! My ID
        int                   id_;

        //! \brief Stop my helpers
        void stopHelpers();

    public:
        /*! \brief Constructor
         * \param[in] mg         Molecule info
         * \param[in] sii        The individual info
         * \param[in] gach       GA Config handler
         * \param[in] bch        Bayes Config handler
         * \param[in] flush      Whether or not to to flush file output directly instead of letting the OS handle it.
         * \param[in] oenv       GROMACS output environment
         */
        ACTMiddleMan(MolGen               *mg,
                     StaticIndividualInfo *sii,
                     GAConfigHandler      *gach,
                     BayesConfigHandler   *bch,
                     bool                  flush,
                     gmx_output_env_t     *oenv);
        
        //! \brief Run the middleman process
        void run();

        //! Print the MCMC statistics if appropriate
        void printStatistics(FILE *logFile);

    };
    
} // namespace alexandria

#endif
