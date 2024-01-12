/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#ifndef ALEXANDRIA_TRAIN_FF_H
#define ALEXANDRIA_TRAIN_FF_H

#include <cstdio>

#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "molgen.h"
#include "bayes.h"
#include "staticindividualinfo.h"
#include "acmfitnesscomputer.h"
#include "acminitializer.h"
#include "act/utility/communicationrecord.h"
#include "act/forces/forcecomputer.h"
#include "act/ga/mutator.h"
#include "act/ga/genetic_algorithm.h"


namespace alexandria
{

/*! \brief Wrapper for closing a file
 * Will print a warning if something is wrong when closing.
 * \param[in] fp The file pointer
 */
void my_fclose(FILE *fp);

/*! \brief FIXME: The class that does all the optimization work.
 * This class inherits the MolGen class that holds molecule data and
 * properties, and the Bayes class, that does the Monte Carlo steps.
 *
 * The class can run as a Master process or as a Helper process.
 */
class OptACM : public ConfigHandler
{

private:
    //! Whether or not to remove molecules that fail to converge in the shell minimization
    bool bRemoveMol_ = true;
    //! Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
    bool flush_ = false;
    //! Print extra stuff during the optimization
    bool verbose_ = false;
    //! Pointer to log file
    gmx::unique_cptr<FILE, my_fclose> fplog_ = nullptr;
    //! ACT Communication data structure
    CommunicationRecord   commRec_;
    //! GROMACS output environment
    gmx_output_env_t     *oenv_ = nullptr;
    //! MolGen instance
    MolGen                mg_;
    //! BayesConfigHandler instance
    BayesConfigHandler    bch_;
    //! GAConfigHandler instance
    GAConfigHandler       gach_;
    //! StaticIndividualInfo instance
    StaticIndividualInfo *sii_;
    //! Base name for FF output
    std::string           baseOutputFileName_;

    // This is for MASTER node
    //! GeneticAlgorithm instance
    ga::GeneticAlgorithm *ga_          = nullptr;

    //! Pointer to ForceComputer
    ForceComputer        *forceComp_   = nullptr;
    //! Pointer to ACMFitnessComputer since it will be initialized later
    ACMFitnessComputer   *fitComp_     = nullptr;
    //! Storing the mutator used
    ga::Mutator          *mutator_     = nullptr;
    
    /*!
     * \brief Print to log file (if it exists), the estimated number of times
     * we will call calcDeviation per dataset
     */
    void printNumCalcDevEstimate();

    /*!
     * \brief Print a table about the parameters in a best genome(s) and a population
     *        to the logfile, if it exists
     * \param[in] genome the best Genome object (for different datasets)
     * \param[in] pop    the population to define summary statistics
     */
    void printGenomeTable(const std::map<iMolSelect, ga::Genome> &genome,
                          const ga::GenePool                     &pop);

public:

    //! Constructor
    OptACM() : mg_(&commRec_)
    {
        sii_ = new StaticIndividualInfo(&commRec_);
    }

    /*! \brief Add my options to the command line
     * \param[inout] pargs   Flags
     * \param[inout] fnms  File names
     */
    virtual void add_options(std::vector<t_pargs>  *pargs,
                             std::vector<t_filenm> *fnms);
    
    /*! \brief Evaluate arguments after parsing.
     */
    void check_pargs();

    /*! \brief Routine to be called after processing options
     * \param[in] filenames The file names
     */
    void optionsFinished(const std::vector<t_filenm> &filenames);

    /*! \brief Routine that opens a log file
     * \param[in] filenms The file names
     */
    void openLogFile(const std::vector<t_filenm> &filenms);

    //! \return a file pointer to the open logfile
    FILE *logFile();

    /*! \brief Initialize charge generation
     * \param[in] ims The data set to do the work for
     */
    void initChargeGeneration(iMolSelect ims);

    /*! \brief Do the actual optimization.
     * \param[in] optimize    If true an optimization will be done
     * \param[in] sensitivity If true, a sensitivity analysis will be done
     * \return true if better parameters were found.
     */
    bool runMaster(bool        optimize,
                   bool        sensitivity);

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Initializing stuff                *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief Initialize the main components of the
     * Genetic Algorithm, just on the master.
     * \param[in] fitnessFile Name of the genetic algorithm files
     */
    void initMaster(const char *fitnessFile);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Initializing stuff                  *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return whether or not we remove problematic compounds
    bool removeMol() const { return bRemoveMol_; }

    //! \return whether or not we are in verbose mode
    bool verbose() const { return verbose_; }

    const CommunicationRecord *commRec() const { return &commRec_; }
    
    /*! \brief Set the output environment pointer \p oenv_
     * \param[in] oenv the reference pointer
     */
    void setOenv(gmx_output_env_t *oenv) { oenv_ = oenv; }

    //! \brief Get the output environment \p oenv_ pointer
    gmx_output_env_t *oenv() { return oenv_; }

    //! \brief Get the MolGen \p mg_ pointer
    MolGen *mg() { return &mg_; }

    //! \brief Get the BayesConfigHandler \p bch_ pointer
    BayesConfigHandler *bch() { return &bch_; }

    //! \brief Get the GAConfigHandler \p gach_ pointer
    GAConfigHandler *gach() { return &gach_; }

    //! \brief Get the StaticIndividualInfo \p sii_ pointer
    StaticIndividualInfo *sii() { return sii_; }

    //! \brief Get the GeneticAlgorithm \p ga_ pointer
    ga::GeneticAlgorithm *ga() { return ga_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};

}


#endif //ALEXANDRIA_TRAIN_EEM_H
