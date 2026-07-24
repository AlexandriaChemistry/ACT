/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
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
#ifndef ACT_FORCECOMPUTERSTATISTICS_H
#define ACT_FORCECOMPUTERSTATISTICS_H

#include <string>

namespace alexandria
{
class ForceComputer;
class CommunicationRecord;

/*! \brief Utitility to compile force computer statistics
 * Note that helper nodes are not included in the statistics.
 */
class ForceComputerStatistics
{
private:
    //! Total number of energy/force evaluations performed (calls to computeOnce)
    size_t totalEval_      = 0;
    //! Total number of shell-minimization iterations performed across all compute() calls
    size_t totalShellIter_ = 0;
    //! Total number of failed shell-minimizations across all compute() calls
    size_t totalShellConvergeFailed_ = 0;
     
public:
    //! Constructor
    ForceComputerStatistics() {}

    /*! \brief Collects statistics from other nodes.
     * and returns a string summarizing it.
     * \param[in] cr        Communication Infrastructure
     * \param[in] forceComp The local force computer. Values will be reset.
     * \param[in] nElites   The number of elite middlemen to ignore
     * \param[in] aggregate Sum over all generations / iterations
     * \return Summarizing string on master node, empty string on others
     */
    std::string statistics(const CommunicationRecord *cr,
                           const ForceComputer       *forceComp,
                           int                        nElites,
                           bool                       aggregate);
};

} // namespace

#endif
