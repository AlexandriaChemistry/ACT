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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#ifndef ALEXANDRIA_ACMFITNESSCOMPUTER_H
#define ALEXANDRIA_ACMFITNESSCOMPUTER_H

#include <vector>

#include "act/ga/fitness_computer.h"
#include "act/forces/forcecomputer.h"
#include "acmindividual.h"
#include "bayes.h"
#include "devcomputer.h"
#include "molgen.h"

namespace gmx
{
class TextWriter;
}

namespace alexandria
{


/*!
 * \brief Computes \f$ \chi^2 \f$ of an individual as its fitness
 */
class ACMFitnessComputer : public ga::FitnessComputer
{

private: 
    //! \brief A pointer to the BoundsDevComputer.
    BoundsDevComputer         *bdc_  = nullptr;
    //! \brief A vector of devComputers.
    std::vector<DevComputer*>  devComputers_;
    //! \brief StaticIndividualInfo pointer
    StaticIndividualInfo      *sii_ = nullptr;
    //! The force computer
    ForceComputer             *forceComp_;
    //! \brief MolGen pointer
    MolGen *molgen_;
    //! \brief Whether or not to remove molecules that fail to converge in the shell minimization
    bool removeMol_;
    //! \brief Amount of times calcDeviation() has been called
    int numberCalcDevCalled_ = 0;

    /*! \brief Compute multipole moments (if needed), for a given molecule
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param actmol     the molecule
     */
    void computeMultipoles(std::map<eRMS, FittingTarget> *targets,
                           ACTMol                        *actmol);

    /*!
     * \brief Fill the devComputers vector according to the needs of the user
     * \param[in] msghandler                  Message Handler
     * \param[in] zetadiff                    Allowed difference in zeta between cores and shells (if both have distributed charges)
     * \param[in] haveInductionCorrectionData Whether or not this energy term is present in the input data
     */
    void fillDevComputers(MsgHandler *msghandler,
                          double      zetaDiff,
                          bool        haveInductionCorrectionData);

public:

    /*!
     * Constructor
     * \param[in] msghandler Message Handler
     * \param[in] sii        pointer to StaticIndividualInfo
     * \param[in] mg         pointer to molgen
     * \param[in] removeMol  Whether or not to remove molecules that fail to converge in the shell minimization
     */
    ACMFitnessComputer(MsgHandler            *msghandler,
                       StaticIndividualInfo  *sii,
                       MolGen                *molgen,
                       const bool             removeMol,
                       ForceComputer         *forceComp)
        : sii_(sii), forceComp_(forceComp), molgen_(molgen), removeMol_(removeMol)
    {
        fillDevComputers(msghandler, molgen->zetaDiff(),
                         molgen->hasMolPropObservable(MolPropObservable::INDUCTIONCORRECTION));
    }

    /*! \brief Do the actual computation
     * \param[in] msghandler Message Handler
     * \param[in] genome    The genome
     * \param[in] trgtFit   The selection to compute
     * \param[in] forceComp The force computer
     */
    void compute(MsgHandler *msghandler,
                 ga::Genome *genome,
                 iMolSelect  trgtFit);

    /*! \brief Distributes the parameters from middlemen to helpers
     * \param[in] params   The force field parameters
     * \param[in] changed  Indication of which parameters have changed and for which
     *                     the forcefield should be updated. If empty, all parameters will
     *                     be updated.
     */
    void distributeParameters(const std::vector<double> *params,
                              const std::set<int>       &changed);
                     
    /*! \brief Distribute the work
     * \param[in] task Tell the middlemen what task to distribute
     * \return what a poor helper is to do    
     */
    CalcDev distributeTasks(CalcDev task);
    
    /*! \brief Computes deviation from target
     * \param[in] msghandler Message Handler
     * \param[in] task       The task at hand
     * \param[in] ims        The dataset to do computations on
     * \return the square deviation
     */
    double calcDeviation(MsgHandler *msghandler,
                         CalcDev     task,
                         iMolSelect  ims);

    //! \return the number of devComputers
    size_t numDevComputers() const { return devComputers_.size(); }
};


} // namespace alexandria


#endif // ALEXANDRIA_ACMFITNESSCOMPUTER_H
