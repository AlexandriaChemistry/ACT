/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
#ifndef ACT_POLINFO_H
#define ACT_POLINFO_H

#include <vector>

namespace alexandria
{
class Topology;
class ForceField;

class PolInfo
{
private:
    // Is this particle a shell?
    std::vector<bool>              isShell_;
    // One over force constant for this particle
    std::vector<double>            fcShell1_;
    // Hyperforce constant
    std::vector<double>            fcHyper_;
    // Number of shells
    size_t                         nShell_ = 0;
public:
    //! Constructor
    PolInfo() {}
    /*! \brief Store polarization information
     * \param[in] pd  The force field
     * \param[in] top A topology
     */
    void init(const ForceField *pd,
              const Topology   *top);
    //! \return the isShell data
    const std::vector<bool> isShell() const { return isShell_; }
    //! \return the fcShell1 data
    const std::vector<double> fcShell1() const { return fcShell1_; }
    //! \return hyperforce constant
    const std::vector<double> &fcHyper() const { return fcHyper_; }
    //! \return number of shells
    size_t nShell() const { return nShell_; }
};

#endif
