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

class GaTestHelper
{
 public:
    CommunicationRecord   cr;
    ACMInitializer       *initializer;
    ACMFitnessComputer   *fitnessComputer;
    ForceComputer        *forceComputer;
    StaticIndividualInfo *sii;
    MolGen               *molgen;
    CompoundReader        compR;
    // Message Handler
    MsgHandler           *msghandler;

    /*! \brief Constructor
     * \param[in] nmiddlemen The number of individuals
     * \param[in] fitstring  Info for fitting
     * \param[in] erms       Stuff to fit
     */
    GaTestHelper(int                                  nmiddlemen,
                 const std::vector<std::string>      &fitstrings,
                 const std::vector<alexandria::eRMS> &erms);
};

} // namespace

#endif
