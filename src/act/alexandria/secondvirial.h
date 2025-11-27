/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024
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

#ifndef ACT_SECONDVIRIAL_H
#define ACT_SECONDVIRIAL_H

#include <random>
#include <vector>

#include "act/alexandria/actmol.h"
#include "act/alexandria/dimergenerator.h"
#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/utility/jsontree.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/statistics/statistics.h"

namespace alexandria
{

/*! \brief Print a force field summary
 * \param[in] jtree JSon tree
 * \param[in] pd    The force field
 */
void forceFieldSummary(JsonTree      *jtree,
                       const ForceField *pd);

//! \brief Components of the second virial
enum class b2Type { Classical, Force, Torque1, Torque2, Total };

//! Map from enum to string for b2Type
extern const std::map<b2Type, std::string> b2Type2str;

/*! \brief Convert b2Type to string
 * \param[in] b2t The b2Type
 * \return The string
 */
const std::string &b2TypeToString(b2Type b2t);

class ReRunner
{
private:
    //! The force computer
    ForceComputer      *forceComp_   = nullptr;
    //! The dimer generator
    DimerGenerator     *gendimers_   = nullptr;
    //! GROMACS output stuff
    gmx_output_env_t   *oenv_        = nullptr;
               
    //! Temperature to start
    double              T1_          = 300;
    //! Final temperature
    double              T2_          = 400;
    //! Temperature step
    double              deltaT_      = 10;
    //! Whether to compute interaction energies and potentially B2
    bool                eInter_      = true;
    //! Whether to compute the second virial
    bool                computeB2_   = false;
    //! Plot only the total B2
    bool                totalOnly_   = false;
    //! The temperature array
    std::vector<double> Temperatures_;
    //! Second virial as a function of T for all components (see function temperatures)
    std::map<b2Type, std::vector<double>> b2t_;
    //! Generate temperature series if needed and return it
    const std::vector<double> &temperatures();
    /*! \brief Generate plot with Second virial as a function of temperature
     * \param[in] b2file    The B2(T) output file name
     */
    void plotB2temp(const char *b2file);


public:
    //! \brief Constructor
    ReRunner(bool computeB2) : computeB2_(computeB2) {}
    
    /*! \brief Make copies of utilities
     * \param[in] forceComp The force computer
     * \param[in] gendimers Dimer generator
     * \param[in] oenv      GROMACS plotting stuff
     */
    void setFunctions(ForceComputer    *forceComp,
                      DimerGenerator   *gendimers,
                      gmx_output_env_t *oenv)
    {
        forceComp_ = forceComp;
        gendimers_ = gendimers;
        oenv_      = oenv;
    }
    
    /*! \brief Add command line options
     * \param[inout] pargs  Regular flags
     * \param[inout] filenm File options
     */
    void addOptions(std::vector<t_pargs>  *pargs,
                    std::vector<t_filenm> *filenm);
    
    //! \return whether or not we will compute interaction energies                
    bool eInteraction() const { return eInter_; }

    //! \brief Manually set the temperatures
    void setTemperatures(double T1, double T2, double deltaT) { T1_ = T1; T2_ = T2; deltaT_ = deltaT; }
    
    //! \brief Manually set a temperature array
    void setTemperatures(const std::vector<double> &T) { Temperatures_ = T; }
    
    //! \brief Reset temperatures
    void resetTemperatures() { Temperatures_.clear(); }
    
    /*! \brief Set the interaction energy flag
     * \param[in] eInter The value
     */
    void setEInteraction(bool eInter) { eInter_ = eInter; }

    /*! \brief Do the rerunning with different options
     * \param[in] msghandler The MsgHandler for printing and status
     * \param[in] pd             Force field structure
     * \param[in] actmol         Structure with molecule info
     * \param[in] verbose        Whether or not to print a lot
     */
    void rerun(MsgHandler       *msghandler,
               const ForceField *pd,
               const ACTMol     *actmol,
               bool              verbose);

    /*! \brief Compute second virial coefficient including QM corrections
     * \param[in] cr         Communication Record for parallel calcs
     * \param[in] msghandler The MsgHandler for printing and status
     * \param[in] pd         Force field structure
     * \param[in] actmol     Structure with molecule info
     * \param[in] maxdimer   Number of dimers to generate (if any)
     * \param[in] fnm        The filenames
     */
    void runB2(CommunicationRecord         *cr,
               MsgHandler                  *msghandler,
               const ForceField            *pd,
               const ACTMol                *actmol,
               int                          maxdimer,
               const std::vector<t_filenm> &fnm);

    //! \return the second virial as a function of T.
    const std::vector<double> &b2Temp(b2Type b2t) const { return b2t_.find(b2t)->second; }

};

} // namespace alexandria

#endif // ACT_SECONDVIRIAL_H
