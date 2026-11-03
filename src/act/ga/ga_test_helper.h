/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
 * Implements testing of the genetic algorithm
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef _GA_TEST_HELPER
#define _GA_TEST_HELPER

#include <string>
#include <vector>

#include "act/alexandria/molgen.h"
#include "act/import/compound_reader.h"
#include "act/utility/communicationrecord.h"
 
namespace alexandria
{

    class ACMInitializer;
    class ACMFitnessComputer;
    class ForceComputer;
    class StaticIndividualInfo;
    class MolGen;
    class MsgHandler;

/*! \brief Set up the components of a genetic algorithm for the tests
 * Reads the ACS-g force field and the alcohol test molecules and creates
 * the initializer and fitness computer.
 */
class GaTestHelper
{
 public:
    //! The communication record
    CommunicationRecord   cr;
    //! Initializer of individuals
    ACMInitializer       *initializer;
    //! Computes the fitness of individuals
    ACMFitnessComputer   *fitnessComputer;
    //! Computes energies and forces
    ForceComputer        *forceComputer;
    //! Information shared by all individuals, including the force field
    StaticIndividualInfo *sii;
    //! The molecules for training
    MolGen               *molgen;
    //! Reader of the molecule input
    CompoundReader        compR;
    //! Message handler
    MsgHandler           *msghandler;

    /*! \brief Constructor
     * \param[in] nmiddlemen The number of individuals
     * \param[in] fitstrings Names of the parameters to fit
     * \param[in] erms       The fitting targets
     */
    GaTestHelper(int                                  nmiddlemen,
                 const std::vector<std::string>      &fitstrings,
                 const std::vector<alexandria::eRMS> &erms);
};

} // namespace

#endif
