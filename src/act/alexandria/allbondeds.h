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

/*! \internal \brief
 * Implements part of the alexandria program.
 * In particular, this is the core of the bastat program which
 * does an analysis of bond and angle statistics.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ALLBONDEDS_H
#define ALLBONDEDS_H

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "act/alexandria/molselect.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield.h"
#include "act/molprop/molprop.h"
#include "act/statistics/normality.h"
#include "act/statistics/statistics.h"
#include "act/topology/actmol.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

namespace gmx
{
class TextWriter;
}

namespace alexandria
{

class JsonTree;

    /*! \brief Where a measured bond length, angle or dihedral comes from
     */
    struct GeometrySample
    {
        //! The value, in pm or degrees
        double           value      = 0;
        //! The molecule, its place in the list of AllBondeds
        int              molecule   = 0;
        //! The calculation, its place in the experiments of the molecule
        int              experiment = 0;
        //! The atoms, their places in the topology of the molecule
        std::vector<int> atoms;
    };

    /*! \brief Base class for managing bonds, angles and dihedrals
     */
    class OneBonded
    {
    private:
        //! Atoms and bondorder
        Identifier id_;
        //! LSQ fitting structure
        gmx_stats  lsq_;
        //! The values, in the order they were added
        std::vector<double>         values_;
        //! Where they come from, if that is kept
        std::vector<GeometrySample> samples_;
    public:
        /*! Constructor
         * \param[in] id The identifier
         */
        OneBonded(const Identifier &id) : id_(id) {}
    
        //! \return my id
        const Identifier &id() const { return id_; }
    
        /*! Add data point to histogram
         * \param[in] x The point to add
         */
        void addPoint(double x);
        /*! Add where a data point comes from, after addPoint
         * \param[in] sample The point and where it comes from
         */
        void addSample(GeometrySample sample) { samples_.push_back(std::move(sample)); }
        //! \return the values, in the order they were added
        const std::vector<double> &values() const { return values_; }
        //! \return where they come from, empty if that was not kept
        const std::vector<GeometrySample> &samples() const { return samples_; }
    
        /*! Write a histogram
         * \param[in] fn_prefix
         * \param[in] xaxis
         * \param[in] oenv
         * \param[in] binwidth
         */
        void writeHistogram(const char             *fn_prefix,
                            const char             *xaxis,
                            const gmx_output_env_t *oenv,
                            double                  binwidth);
        /*! \brief Extract statistics from this parameter
         *! \todo Implement support for extracting the median value
         * \param[out] average The average value    
         * \param[out] sigma   The standard deviation
         * \param[out] N       The number of data points
         * \return status from statistics library
         */
        eStats getAverageSigmaN(real   *average,
                                real   *sigma,
                                size_t *N);
    };

    /*! Class to analyze structures and derive geometric statistics
     */
    class AllBondeds
    {
    private:
        //! Map from interaction type to list of bondeds
        std::map<InteractionType, std::vector<OneBonded> > bondeds_;
        //! Default dissociation energy
        real                                               De_        = 400;
        //! Default steepness for Morse potentials
        real                                               beta_      = 20;
        //! Default force constant for harmonic angles
        real                                               kt_        = 300;
        //! Default force constant for harmonic bonds
        real                                               kb_        = 300000;
        //! Default force constant for harmonic linear angles
        real                                               klin_      = 20000;
        //! Default force constant for dihedral angles
        real                                               kp_        = 5;
        //! Default force constant for improper dihedrals
        real                                               kimp_      = 1;
        //! Default force constant for Urey-Bradley
        real                                               kub_       = 30000;
        //! Tolerance for warning about large sigma in bond-lengths (pm)
        real                                               bond_tol_  = 8;
        //! Tolerance for warning about large sigma in angles (degrees)
        real                                               angle_tol_ = 8;
        //! Scaling factor for setting min and max for force parameters
        real                                               factor_    = 0.5;
        //! Scaling factor for setting min and max for bond lengths
        real                                               bfactor_   = 0.95;
        //! Scaling factor for setting min and max for bond angles
        real                                               afactor_   = 0.95;
        //! Spacing for bond histograms in pm
        real                                               bspacing_  = 1;
        //! Spacing for angle histograms in degrees
        real                                               aspacing_  = 0.5;
        //! Relative number for linear angles
        real                                               laspacing_ = 0.000001;
        //! Spacing for dihedrals in degrees
        real                                               dspacing_  = 1;
        //! Whether to test the distributions of bond lengths and angles for normality
        gmx_bool                                           analyse_   = false;
        //! Significance level of the test for normality
        real                                               alpha_     = 0.01;
        //! What is known about the distribution of the values of one identifier
        struct Distribution
        {
            //! The identifier
            std::string     id;
            //! The average and the standard deviation as they are in the log file
            double          average = 0;
            double          sigma   = 0;
            //! Whether the standard deviation is above the tolerance
            bool            wide    = false;
            //! The moments of the values and the test for normality
            NormalityResult normality;
        };
        //! The distributions that were analysed, in the order of the log file
        std::map<InteractionType, std::vector<Distribution> > distributions_;
        //! Whether to keep where each value comes from
        bool                                               keepSamples_ = false;
        //! A molecule from which values come: its name, the data files of
        //! its calculations and the names of its atoms
        struct SampledMolecule
        {
            std::string              name;
            std::vector<std::string> datafiles;
            std::vector<std::string> atomNames;
        };
        //! The molecules from which values come, if that is kept
        std::vector<SampledMolecule>                       molecules_;
        /*! \brief Add bonds etc. for one molecule and list of atoms
         * \param[in] msghandler Message handler
         * \param[in] iType   InteractionType
         * \param[in] mmi     Molecule structure
         * \param[in] bondId  The bond identifier
         * \param[in] atomid  List of atoms involved in the interaction
         * \param[in] x       The coordinates
         * \param[in] experiment The calculation of the coordinates, its place in
         *                       the experiments of the molecule
          */
        void addBonded(MsgHandler                   *msghandler,
                       InteractionType               iType,
                       const ACTMol                 &mmi,
                       const Identifier             &bondId,         
                       const std::vector<int>       &atomid,
                       const std::vector<gmx::RVec> &x,
                       int                           experiment = 0);

        /*! \brief Add bonds etc. for one molecule and list of atoms
         * \param[in] msghandler Message handler
         * \param[in] iType   InteractionType
         * \param[in] mmi     Molecule structure
         * \param[in] bondId  The bond identifier
         * \param[in] atomid  List of atoms involved in the interaction
          */
        void addBondeds(MsgHandler             *msghandler,
                        InteractionType         iType,
                        const ACTMol           &mmi,
                        const Identifier       &bondId,         
                        const std::vector<int> &atomid);

        /*! \brief Analyse the distribution of the values of one identifier, with -analyse
         * \param[in] iType     InteractionType
         * \param[in] ob        The values
         * \param[in] average   The average as it is in the log file
         * \param[in] sigma     The standard deviation as it is in the log file
         * \param[in] tolerance The largest standard deviation without a warning
         * \param[in] test      Whether to test for normality
         * \return what to add to the line of the identifier in the log file:
         *         a warning if the values do not come from a normal
         *         distribution, nothing otherwise or without -analyse
         */
        std::string analyseDistribution(InteractionType  iType,
                                        const OneBonded &ob,
                                        double           average,
                                        double           sigma,
                                        double           tolerance,
                                        bool             test);

    public:
        //! Constructor
        AllBondeds() = default;
        /*! \brief Keep where each value comes from, see histogramsToJson
         * \param[in] keep Whether to keep it
         */
        void setKeepSamples(bool keep) { keepSamples_ = keep; }
        /*! \brief The histograms of the values, in JSON
         * For each interaction type a list with, for each identifier, its
         * unit, the width of a bin, the number of values and the bins that
         * have values: the middle of the bin and the number of values in
         * it. With samples kept, a bin lists up to maxSamples of its values
         * with the molecule, the data file of the calculation and the atoms
         * (numbers from one and names with numbers) they come from, and the
         * number of the others (more).
         * \param[in] maxSamples The most values that a bin lists
         * \return the histograms, with the key histograms
         */
        JsonTree histogramsToJson(size_t maxSamples) const;

        //! \return whether the distributions are analysed (-analyse)
        bool analyse() const { return analyse_; }

        //! \return the significance level of the test for normality (-alpha)
        double alpha() const { return alpha_; }

        /*! \brief How many distributions are not normal, after updateForceField
         * \return the number of identifiers, of all interaction types, of
         *         which the values do not come from a normal distribution
         */
        size_t notNormal() const;

        /*! \brief The analysis of the distributions of bond lengths and angles, in JSON
         * The test that was used, the significance level (alpha) and the
         * fewest values that are tested (minimum), and for each interaction
         * type its unit, the tolerance for the standard deviation, the number
         * of distributions, how many are wide (standard deviation above
         * the tolerance), how many were tested and how many are not normal,
         * and for each identifier (entries) the number of values (count),
         * average, sigma, wide, tested and, if tested, normal, the statistic
         * A2, pvalue, skewness and kurtosis.
         * \return the analysis, with the key analysis
         */
        JsonTree analysisToJson() const;
    
        /*! \brief Add command line options for class variables
         * \param[inout] pargs The command line parameters
         */
        void addOptions(std::vector<t_pargs> *pargs);

        /*! \brief Write histograms for bond length distributions etc.
         * \param[in] oenv GROMACS output environment structure
         */
        void writeHistogram(const gmx_output_env_t *oenv);

        /*! \brief How many different bonds, angles etc. were extracted
         * \return for each interaction type the number of identifiers
         */
        std::map<InteractionType, size_t> counts() const;
            
        /*! \brief Store the bond lengths etc. in the force field
         * With -analyse the distributions of the bond lengths and angles
         * are analysed as well, see analysisToJson.
         * \param[in]  tw  To print information to
         * \param[out] pd  The force field structure to update
         */
        void updateForceField(gmx::TextWriter *tw,
                              ForceField      *pd);
             
        /*! \brief Extract bond lengths, angles etc. from molecules
         * \param[in]  msghandler Message handler
         * \param[in]  mp     MolProp array (will be destroyed)
         * \param[out] actmols ACTMol array will be filled here
         * \param[in]  pd     Force field structure
         * \param[in]  gms    Selection of compounds
         */                          
        void extractGeometries(MsgHandler           *msghandler,
                               std::vector<MolProp> *mp,
                               std::vector<ACTMol>  *actmols,
                               ForceField           *pd,
                               const MolSelect      &gms);

        /*! \brief Write how many bonds etc. were found, and with -analyse
         * how many distributions are not normal
         * \param[in] tw Text writer to write to
         */
        void writeSummary(gmx::TextWriter *tw);

    };
    
}

#endif
