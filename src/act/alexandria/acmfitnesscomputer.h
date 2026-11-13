/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#include <array>
#include <map>
#include <vector>

#include "act/ga/fitness_computer.h"
#include "act/forces/forcecomputer.h"
#include "acmindividual.h"
#include "devcomputer.h"
#include "molgen.h"
#include "sensitivity.h"

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
    //! A vector of devComputers.
    std::vector<DevComputer*>  devComputers_;
    //! StaticIndividualInfo pointer
    StaticIndividualInfo      *sii_ = nullptr;
    //! The force computer
    ForceComputer             *forceComp_;
    //! MolGen pointer
    MolGen *molgen_;
    //! Whether or not to remove molecules that fail to converge in the shell minimization
    bool removeMol_;
    //! Charge generation algorithm
    ChargeGenerationAlgorithm  algorithm_;
    //! Amount of times calcDeviation() has been called
    int numberCalcDevCalled_ = 0;
    /*! For each molecule, the version of the parameters of each interaction
     * type that its topologies were filled with, see
     * StaticIndividualInfo::parameterVersion.
     */
    std::vector<std::map<InteractionType, size_t>> filledVersion_;
    //! For each molecule, the parameter version that its charges were generated with
    std::vector<size_t>                            chargeVersion_;
    //! The parameter version for which only the fixed charges were set, as the GPU generates the others
    std::vector<size_t>                            fixedChargeVersion_;
    //! The timer of each DevComputer, see Timings
    std::vector<size_t>                            devTimers_;
    //! The evaluations compared and those that differed, see ACT_GPU_DEVIATION_CHECK in finishDeviation
    int                                            compactChecks_ = 0, compactDiffer_ = 0;
    //! The same for the deviations computed on the GPU
    int                                            gpuChecks_ = 0, gpuDiffer_ = 0;
    //! The optimized parameters that each molecule on the GPU depends on, see buildDependencies
    std::vector<std::vector<int>>                  molDeps_;
    //! The molecules of molDeps_ and the version of the GPU system, see ForceComputer::gpuVersion
    std::vector<const ACTMol *>                    depsMols_;
    size_t                                         depsVersion_ = 0;
    //! The deviations of a molecule for a set of the parameters it depends on
    struct MolCacheEntry
    {
        bool                  valid = false;
        std::vector<double>   key;
        GpuMoleculeDeviations dev;
        //! The penalty of the charges, see BoundsDevComputer::moleculeBound
        double                bound = 0;
    };
    //! For each replica and molecule, the deviations of the last two sets of its parameters
    std::vector<std::vector<std::array<MolCacheEntry, 2>>> molCache_;
    //! The version of the deviations on the GPU of molCache_, see ForceComputer::gpuDeviationsVersion
    size_t                                         cacheVersion_ = 0;
    //! The evaluations the chain on the GPU was compared with and those that differed, see checkChain
    int                                            chainChecks_ = 0, chainDiffer_ = 0;

    //! The number of molecules computed and skipped, for the log
    size_t                                         nComputedMol_ = 0, nSkippedMol_ = 0;
    //! Whether ForceComputer::gpuResidentWhyNot was written to the log
    bool                                           gpuWhyNotReported_ = false;
    //! For each dataset, whether its molecules have data of the targets of the charge properties, see needChargeProperties
    std::map<iMolSelect, std::map<eRMS, bool>>     qPropsData_;


    /*! \brief Update the topology and the charges of a molecule for the parameters now
     * As in prepareDeviation.
     * \param[in]    msghandler Message Handler
     * \param[inout] actmol     The molecule
     * \param[in]    imol       Its index in the molecules of MolGen
     * \param[in]    itUpdate   The interaction types that are optimized
     * \param[in]    gpuCharges Whether the GPU generates the charges, see prepareDeviation
     * \param[in]    needQProps Whether the charge properties are needed
     * \param[in]    allTypes   Whether to fill the parameters of all interaction
     *                          types in the topology, or only of those that the
     *                          processor needs when the GPU makes the rows of
     *                          its parameter tables, see hostNeedsParameters.
     *                          The others then keep their version, such that
     *                          they are filled the next time all are.
     */
    void updateMolecule(MsgHandler                      *msghandler,
                        ACTMol                          *actmol,
                        size_t                           imol,
                        const std::set<InteractionType> &itUpdate,
                        bool                             gpuCharges,
                        bool                             needQProps,
                        bool                             allTypes = true);

    /*! \brief Find which molecules on the GPU each optimized parameter changes
     * Each parameter is changed in turn, and the molecules of which the
     * parameters on the GPU change depend on it, see
     * ForceComputer::gpuMoleculeSignatures. Only those molecules have to be
     * computed again when it changes. The parameters and molecules are
     * restored.
     * \param[in] msghandler Message Handler
     * \param[in] selected   The molecules on the GPU
     * \param[in] itUpdate   The interaction types that are optimized
     * \param[in] needQProps Whether the charge properties are needed
     */
    void buildDependencies(MsgHandler                      *msghandler,
                           const std::vector<ACTMol *>     &selected,
                           const std::set<InteractionType> &itUpdate,
                           bool                             needQProps);

    //! \return the values of the optimized parameters in the force field now
    std::vector<double> parameterValues();
    //! The loss function
    LossFunction lossFunction_;

    /*! \brief Compute multipole moments (if needed), for a given molecule
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param actmol     the molecule
     */
    void computeMultipoles(std::map<eRMS, FittingTarget> *targets,
                           ACTMol                        *actmol);

    /*!
     * \brief Fill the devComputers vector according to the needs of the user
     * \param[in] msghandler                  Message Handler
     * \param[in] haveInductionCorrectionData Whether or not this energy term is present in the input data
     */
    void fillDevComputers(MsgHandler *msghandler,
                          bool        haveInductionCorrectionData);

public:

    /*! \brief Whether to compare ForceEnergyDevComputer::calcDeviationCompact with calcDeviation
     * For developers, with the environment variable ACT_GPU_DEVIATION_CHECK.
     */
    static bool compactCheck();

    /*! \brief Whether the charge properties of the molecules are needed
     * They are needed for the targets of the ESP and the multipoles that
     * have a weight, if the molecules of the dataset have data for them. A
     * target without data adds nothing to the objective function, but would
     * keep the charges on the processor, see prepareDeviation, which is said
     * once.
     * \param[in] msghandler For the warning
     * \param[in] targets    The fitting targets of the dataset
     * \param[in] ims        The dataset
     * \return whether they are needed
     */
    bool needChargeProperties(MsgHandler                          *msghandler,
                              const std::map<eRMS, FittingTarget> &targets,
                              iMolSelect                           ims);

    //! \brief Constructor
    ACMFitnessComputer() {}

    //! \brief Destructor
    virtual ~ACMFitnessComputer();

    /*! Initalize stuff
     * \param[in] msghandler   Message Handler
     * \param[in] sii          pointer to StaticIndividualInfo
     * \param[in] molgen       pointer to molgen
     * \param[in] removeMol    Whether or not to remove molecules that fail to converge in the shell minimization
     * \param[in] forceComp    Pointer to ForceComputer
     * \param[in] algorithm    Charge generation algorithm
     * \param[in] lossFunction The loss function to use
     */
    void init(MsgHandler                *msghandler,
              StaticIndividualInfo      *sii,
              MolGen                    *molgen,
              const bool                 removeMol,
              ForceComputer             *forceComp,
              ChargeGenerationAlgorithm  algorithm,
              LossFunction               lossFunction)
    {
        sii_          = sii;
        forceComp_    = forceComp;
        if (compactCheck())
        {
            forceComp_->setKeepGpuEnergies(true);
        }
        molgen_       = molgen;
        removeMol_    = removeMol;
        algorithm_    = algorithm;
        lossFunction_ = lossFunction;
        if (molgen)
        {
            fillDevComputers(msghandler,
                             molgen->hasMolPropObservable(MolPropObservable::INDUCTIONCORRECTION));
        }
    }

    /*! \brief Do the actual computation
     * \param[in] msghandler Message Handler
     * \param[in] genome    The genome
     * \param[in] trgtFit   The selection to compute
     */
    void compute(MsgHandler                *msghandler,
                 ga::Genome                *genome,
                 iMolSelect                 trgtFit);

    /*! \brief Distributes the parameters from middlemen to helpers
     * \param[in] msghandler Message Handler
     * \param[in] params   The force field parameters
     * \param[in] changed  Indication of which parameters have changed and for which
     *                     the forcefield should be updated. If empty, all parameters will
     *                     be updated.
     */
    void distributeParameters(MsgHandler                *msghandler,
                              const std::vector<double> *params,
                              const std::set<int>       &changed);
                     
    /*! \brief Distribute the work
     * \param[in] task Tell the middlemen what task to distribute
     * \return what a poor helper is to do    
     */
    CalcDev distributeTasks(CalcDev task);
    
    /*! \brief Whether the steps of the Monte Carlo chain can run on the GPU
     * After an evaluation of the training set on this node alone, whose
     * deviations the GPU computed and whose molecules may be skipped, see
     * runChain. ACT_NO_GPU_CHAIN turns it off, and so does ACT_GPU_CHAIN_CHECK.
     */
    bool chainReady() const { return chainOk_; }

    /*! \brief Run steps of the Monte Carlo chain on the GPU, see MCMCMutator::mutate
     * Each step changes a parameter, evaluates the training set as
     * calcDeviation and applies the Metropolis criterion, all on the GPU. The
     * force field gets the parameters at the end.
     * \param[in]  msghandler For messages
     * \param[in]  params     The parameters, of chi2
     * \param[in]  steps      The steps, with their random numbers
     * \param[in]  chi2       The chi squared of params, from the evaluation before
     * \param[in]  haveBest   Whether there is a lowest chi squared
     * \param[in]  bestChi2   The lowest chi squared
     * \param[out] results    The results of the steps
     */
    void runChain(MsgHandler                         *msghandler,
                  const std::vector<double>          &params,
                  const std::vector<gpu::ChainStep>  &steps,
                  double                              chi2,
                  bool                                haveBest,
                  double                              bestChi2,
                  gpu::ChainResults                  *results);

    /*! \brief Computes deviation from target
     * As prepareDeviation, computeEnergies and finishDeviation.
     * \param[in] msghandler   Message Handler
     * \param[in] task         The task at hand
     * \param[in] ims          The dataset to do computations on
     * \return the deviation
     */
    double calcDeviation(MsgHandler   *msghandler,
                         CalcDev       task,
                         iMolSelect    ims);

    //! The state of a deviation calculation between its steps, see calcDeviation
    struct DeviationTask
    {
        //! The task at hand
        CalcDev                        task    = CalcDev::Compute;
        //! The dataset
        iMolSelect                     ims     = iMolSelect::Train;
        //! The fitting targets of the dataset
        std::map<eRMS, FittingTarget> *targets = nullptr;
        //! The molecules of the dataset that this node computes
        std::vector<ACTMol *>          selected;
        //! The number of molecules in the dataset
        int                            ntrain  = 0;
        //! The number of these computed here
        int                            nlocal  = 0;
        //! Whether the GPU computes the energies and the charges, see ForceComputer::gpuChargesOnly
        bool                           gpuCharges = false;
        //! Whether the GPU computed the deviations of the energies, and not the energies of the molecules
        bool                           gpuDeviations = false;
        //! The replica, for several sets of parameters at once, see calcDeviations
        size_t                         replica = 0;
        //! Whether molecules whose parameters did not change are skipped, see prepareDeviation
        bool                           skipping = false;
        //! For each molecule: whether it is computed, the values of the parameters
        //! it depends on, and the entry of molCache_ with its deviations if it is not computed
        std::vector<char>                active;
        std::vector<std::vector<double>> keys;
        std::vector<int>                 hit;
        //! Whether the keys are known, such that the deviations can be kept
        bool                             haveKeys = false;
        //! Whether the chain on the GPU checks this evaluation, see checkChain
        bool                             chainCheck = false;
    };

    /*! \brief Make the Monte Carlo chain on the GPU for the molecules of the dependencies
     * \param[in] dt The evaluation, after prepareDeviation
     * \return Whether the GPU has it
     */
    bool setupChain(const DeviationTask &dt);

    /*! \brief The cache of a replica as the chain on the GPU has it, see GpuTrainingSystem::chainCache
     * \param[in]  replica The replica
     * \param[out] values  The values of the sets
     * \param[out] valid   Whether each is valid
     */
    void chainCacheOf(size_t replica, std::vector<double> *values, std::vector<char> *valid) const;

    /*! \brief For developers, compare the chain on the GPU with this evaluation
     * With ACT_GPU_CHAIN_CHECK set: after prepareDeviation the molecules
     * skipped, after finishDeviation the chi squared and the sets kept.
     * \param[in] msghandler For messages
     * \param[in] dt         The evaluation
     * \param[in] chi2       The chi squared, after finishDeviation, otherwise negative
     */
    void checkChain(MsgHandler *msghandler, DeviationTask *dt, double chi2);

    /*! \brief The sets kept of a replica from those of the chain on the GPU, see chainCacheOf
     * \param[in] replica The replica
     * \param[in] values  The values of the sets
     * \param[in] valid   Whether each is valid
     */
    void chainCacheTo(size_t replica, const std::vector<double> &values, const std::vector<char> &valid);

    //! Whether the evaluation before allows the steps of the chain on the GPU, see chainReady
    bool                  chainOk_ = false;
    //! The molecules of that evaluation
    std::vector<ACTMol *> chainMols_;

    /*! \brief Whether an evaluation allows the steps of the chain on the GPU
     * \param[in] dt The evaluation, after finishDeviation
     * \return Whether it does
     */
    bool chainAllowed(const DeviationTask &dt);

    /*! \brief Compute the deviations of several sets of parameters
     * As distributeParameters and calcDeviation for each set, at once on
     * the GPU if the force computer has a replica for each, see
     * ForceComputer::setGpuReplicas, and only energies are needed. The
     * force field has the parameters of one of the sets afterwards.
     * \param[in] msghandler Message Handler
     * \param[in] params     The parameters of each set, all of them
     * \param[in] ims        The dataset to do computations on
     * \return the deviation of each set
     */
    std::vector<double> calcDeviations(MsgHandler                               *msghandler,
                                       const std::vector<const std::vector<double> *> &params,
                                       iMolSelect                                ims);

    /*! \brief The first step of calcDeviation, before computing energies
     * Receives or sends the dataset, resets the chi squared, updates the
     * topologies and charges of the molecules with the newest parameters
     * and computes their multipoles.
     * \param[in]  msghandler Message Handler
     * \param[in]  task       The task at hand
     * \param[in]  ims        The dataset to do computations on
     * \param[out] dt         The state for the next steps
     * \return false if there are no targets for the dataset, the deviation is zero then
     */
    bool prepareDeviation(MsgHandler    *msghandler,
                          CalcDev        task,
                          iMolSelect     ims,
                          DeviationTask *dt);

    /*! \brief The energy calculation that finishDeviation needs
     * \param[in] dt The state of prepareDeviation
     * \return what to compute for the molecules of dt, nothing without energy targets
     */
    ForceEnergyDevComputer::EnergyTask energyTask(const DeviationTask &dt) const;

    /*! \brief The second step of calcDeviation, the energies of the molecules
     * Energies of the molecules of several calculations can be computed at
     * once with ForceComputer::computeRelaxedEnergies, as energyTask says,
     * instead of this.
     * \param[in]    msghandler Message Handler
     * \param[inout] dt         The state of prepareDeviation, which tells whether the
     *                          GPU computed only the deviations, see useGpuDeviations
     */
    void computeEnergies(MsgHandler *msghandler, DeviationTask *dt);

    /*! \brief Whether the GPU can compute the deviations of the energies of these molecules
     * Instead of the energies of the molecules, see
     * ForceEnergyDevComputer::setupGpuDeviations, for evaluations in
     * training, not for ACT_GPU_DEVIATION_CHECK, which compares them, nor
     * with ACT_NO_GPU_DEVIATIONS, for developers.
     * \param[in] dt The state of prepareDeviation
     * \return whether the deviations are set up on the GPU
     */
    bool useGpuDeviations(const DeviationTask &dt);

    /*! \brief For developers, let the GPU compute the energy deviations of these molecules
     * With ACT_GPU_DEVIATION_CHECK, see ForceEnergyDevComputer::setupGpuDeviations.
     * \param[in] dt The state of prepareDeviation
     */
    void setupGpuDeviations(const DeviationTask &dt);

    /*! \brief The last step of calcDeviation, from the energies computed
     * \param[in]    msghandler Message Handler
     * \param[inout] dt         The state of prepareDeviation
     * \return the deviation
     */
    double finishDeviation(MsgHandler *msghandler, DeviationTask *dt);
    
    //! \return the number of devComputers
    size_t numDevComputers() const { return devComputers_.size(); }

    //! \return forceComputer
    const ForceComputer *forceComputer() const { return forceComp_; }
};


} // namespace alexandria


#endif // ALEXANDRIA_ACMFITNESSCOMPUTER_H
