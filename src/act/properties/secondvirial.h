/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024,2026
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

#include <map>
#include <string>
#include <vector>

#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/import/compound_reader.h"
#include "act/properties/dimergenerator.h"
#include "act/topology/actmol.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/jsontree.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/real.h"

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

/*! \brief Recompute energies and forces for (generated or read) dimer
 * structures, and optionally compute the second virial coefficient
 * as a function of temperature
 */
class ReRunner
{
private:
    //! The force computer
    ForceComputer      *forceComp_   = nullptr;
    //! The dimer generator
    DimerGenerator      gendimers_;
    //! GROMACS output stuff
    gmx_output_env_t   *oenv_        = nullptr;
               
    //! Temperature to start
    real                T1_          = 300;
    //! Final temperature
    real                T2_          = 400;
    //! Temperature step
    real                deltaT_      = 10;
    //! Whether to compute interaction energies and potentially B2
    bool                eInter_      = true;
    //! Whether to compute the second virial
    bool                computeB2_   = false;
    //! Plot only the total B2
    bool                totalOnly_   = false;
    //! Dimers with a lower energy (kJ/mol) are treated as overlapping atoms
    real                eMin_        = -1000;
    //! Largest number of distances of an orientation computed together, 0 means all
    int                 blockSize_   = 264;
    //! The file for the Mayer functions (-eh), none if empty
    std::string         mayerFile_;
    //! The file for the second virial as a function of the temperature (-b2), none if empty
    std::string         b2File_;
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
                      gmx_output_env_t *oenv)
    {
        forceComp_ = forceComp;
        oenv_      = oenv;
    }

    //! \return the dimer generator
    DimerGenerator *gendimers() { return &gendimers_; }
    //! \return the dimer generator
    const DimerGenerator &gendimersConst() const { return gendimers_; }

    /*! \brief Add command line options
     * \param[inout] pargs  Regular flags
     * \param[inout] filenm File options
     * \param[inout] desc   Descriptive text lines
     */
    void addOptions(std::vector<t_pargs>      *pargs,
                    std::vector<t_filenm>     *filenm,
                    std::vector<const char *> *desc);
    
    /*! \brief Take the settings of the command line
     * The output files, then see DimerGenerator::finishOptions and validate.
     * \param[in] filenm The file options after parsing
     * \throws gmx::InvalidInputError for settings that make no sense
     */
    void finishOptions(const std::vector<t_filenm> &filenm);

    /*! \brief Check the settings
     * \throws gmx::InvalidInputError that names the option of a setting that is out of range
     */
    void validate() const;

    /*! \brief Set the temperatures, as -T1, -T2 and -dT
     * \param[in] T1     The lowest temperature in K
     * \param[in] T2     The highest temperature in K
     * \param[in] deltaT The step in K
     */
    void setTemperatures(double T1, double T2, double deltaT)
    {
        T1_     = T1;
        T2_     = T2;
        deltaT_ = deltaT;
        Temperatures_.clear();
    }

    /*! \brief Set the energy below which dimers are taken as overlapping, as -emin
     * \param[in] eMin The energy in kJ/mol
     */
    void setEMin(double eMin) { eMin_ = eMin; }

    /*! \brief Set the largest number of distances computed together, as -block
     * \param[in] blockSize The number, 0 for all
     */
    void setBlockSize(int blockSize) { blockSize_ = blockSize; }

    /*! \brief Set the output files, as -eh and -b2
     * \param[in] mayerFile The Mayer functions, none if empty
     * \param[in] b2File    The second virial as a function of the temperature, none if empty
     */
    void setOutputFiles(const std::string &mayerFile, const std::string &b2File)
    {
        mayerFile_ = mayerFile;
        b2File_    = b2File;
    }

    /*! \brief Set the interaction energy flag
     * \param[in] eInter The value
     */
    void setEInteraction(bool eInter) { eInter_ = eInter; }

    /*! \brief Do the rerunning with different options
     * \param[in] msghandler The MsgHandler for printing and status
     * \param[in] pd             Force field structure
     * \param[inout] actmol      Structure with molecule info, the structures
     *                           are stored in it as its conformations
     * \param[in] verbose        Whether or not to print a lot
     */
    void rerun(MsgHandler       *msghandler,
               const ForceField *pd,
               ACTMol           *actmol,
               bool              verbose);

    /*! \brief Compute second virial coefficient including QM corrections
     * \param[in] cr         Communication Record for parallel calcs
     * \param[in] msghandler The MsgHandler for printing and status
     * \param[in] pd         Force field structure
     * \param[inout] actmol   Structure with molecule info, the dimers are
     *                        stored in it as its conformations
     * \param[in] maxdimer   Number of dimers to generate (if any)
     */
    void runB2(CommunicationRecord *cr,
               MsgHandler          *msghandler,
               const ForceField    *pd,
               ACTMol              *actmol,
               int                  maxdimer);

    /*! \brief Add the second virial coefficients to a JSON tree
     * For each temperature the components and the total in cm^3/mol, with
     * all digits, such that results can be compared exactly, as in the
     * tests. Nothing is added if the second virial was not computed.
     * \param[inout] jtree The tree
     */
    void addB2(JsonTree *jtree);

};

/*! \brief What alexandria b2 does, as its command line gives it
 *
 * The program fills it from the command line, a program without one sets
 * the members. runSecondVirial checks it, see check.
 */
struct B2Options
{
    //! The force field file (-ff)
    std::string    forceFieldFile;
    //! The file for part of the output in JSON, with the second virials (-json), none if empty
    std::string    jsonFile;
    //! The number of orientations of the dimers, unless they come from a trajectory or a grid (-maxdimer)
    int            maxDimers      = 1024;
    //! The tolerance for relaxing the shells (-shelltoler)
    double         shellTolerance = 1e-6;
    //! The largest number of iterations for relaxing the shells (-shellMaxIter)
    int            shellMaxIter   = 25;
    //! The temperatures, the dimers, and the output files of the second virial (-T1, -randalg, -b2 and so on)
    ReRunner       rerun { true };
    //! The dimer and the charges of its compounds (-f, -charges, -qalg and so on)
    CompoundReader compounds;

    /*! \brief Check that the options make sense
     * \throws gmx::InvalidInputError that says which option is wrong
     */
    void check() const;
};

/*! \brief The second virial coefficient of a dimer, as alexandria b2 computes it
 * \param[in] cr         The processes, which share the orientations
 * \param[in] msghandler For output and errors
 * \param[in] options    What to do
 * \param[in] oenv       The output environment for the xvg files; one with
 *                       the defaults is made if it is nullptr
 * \return 0 if it went well, 1 if reading the input failed
 * \throws gmx::InvalidInputError for options that make no sense, see
 * B2Options::check
 */
int runSecondVirial(CommunicationRecord *cr,
                    MsgHandler          *msghandler,
                    const B2Options     &options,
                    gmx_output_env_t    *oenv = nullptr);

} // namespace alexandria

#endif // ACT_SECONDVIRIAL_H
