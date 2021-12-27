#include "GeneticAlgorithm.h"

#include <cstdio>

#include "Initializer.h"
#include "FitnessComputer.h"
#include "ProbabilityComputer.h"
#include "Crossover.h"
#include "Mutator.h"
#include "Terminator.h"


#include "alexandria/acmindividual.h"
#include "alexandria/mcmcmutator.h"


namespace ga
{


void GeneticAlgorithm::fprintPop() const
{
    fprintf(logfile_, "Population:\n");
    for (Individual *ind : oldPop_) ind->fprintSelf(logfile_);
}

void GeneticAlgorithm::fprintBestInd() const
{
    fprintf(logfile_, "Overall Best Individual:\n");
    bestInd_->fprintSelf(logfile_);
}

void GeneticAlgorithm::fprintBestIndInPop() const
{
    fprintf(logfile_, "Best Individual in current population:\n");
    oldPop_[findBestIndex()]->fprintSelf(logfile_);
}

int GeneticAlgorithm::findBestIndex() const
{
    int index = 0;
    double bestFitness = oldPop_[index]->fitnessTrain();
    for (int i = 1; i < oldPop_.size(); i++)
    {
        if (oldPop_[i]->fitnessTrain() < bestFitness)
        {
            index = i;
            bestFitness = oldPop_[i]->fitnessTrain();
        }
    }
    return index;
}

void GeneticAlgorithm::fprintProbability() const
{
    fprintf(logfile_, "Probability: [ ");
    for (Individual *ind : oldPop_) fprintf(logfile_, "%f ", ind->probability());
    fprintf(logfile_, "]\n");
}

void GeneticAlgorithm::fprintFitness() const
{
    for (int i = 0; i < oldPop_.size() - 1; i++)
    {
        fprintf(fileFitnessTrain_, "%lf ", oldPop_[i]->fitnessTrain());
        fprintf(fileFitnessTest_, "%lf ", oldPop_[i]->fitnessTest());
    }
    fprintf(fileFitnessTrain_, "%lf\n", oldPop_[oldPop_.size()-1]->fitnessTrain());
    fprintf(fileFitnessTest_, "%lf\n", oldPop_[oldPop_.size()-1]->fitnessTest());
}

void GeneticAlgorithm::evolveMCMC()
{

    // Simplify syntax
    using alexandria::ACMIndividual;
    using alexandria::MCMCMutator;

    // Initialize population/s
    for (int i = 0; i < gach_->popSize(); i++)
    {
        initializer_->initialize(&(oldPop_[i]));
    }

    // Cast each individual to ACMIndividual for easier evolution
    std::vector<ACMIndividual*> acmPop;
    for (Individual *ind : oldPop_) acmPop.push_back(static_cast<ACMIndividual*>(ind));

    // Open files of each individual
    for (ACMIndividual *ind : acmPop)
    {
        ind->openParamConvFiles(oenv_);
        ind->openChi2ConvFile(oenv_, bch_->evaluateTestset());
    }

    // Evolve each individual
    MCMCMutator *acmMut = static_cast<MCMCMutator*>(mutator_);
    for (ACMIndividual *ind : acmPop) acmMut->MCMC(ind, bch_->evaluateTestset());

    // Collect results into the best individual
    int bestIndex = 0;
    double bestFitness = oldPop_[0]->fitnessTrain();
    for (int i = 0; i < gach_->popSize(); i++)
    {
        if (acmPop[i]->fitnessTrain() < bestFitness)
        {
            bestFitness = acmPop[i]->fitnessTrain();
            bestIndex = i;
        }
    }
    bestInd_ = new ACMIndividual(acmPop[bestIndex]);

    // Close files of each individual
    for (ACMIndividual *ind : acmPop) ind->closeConvFiles();

}

void GeneticAlgorithm::evolveGA()
{
    
    fprintf(logfile_, "\nStarting GA/HYBRID evolution\n");

    // Open surveillance files for fitness
    fileFitnessTrain_ = fopen(filenameFitnessTrain_, "w");
    fileFitnessTest_  = fopen(filenameFitnessTest_, "w");
    GMX_RELEASE_ASSERT(fileFitnessTrain_ != NULL && fileFitnessTest_ != NULL,
                       "GeneticAlgorithm: error opening the fitness output files.");

    // Random number generation
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Iteration variables
    int i, k;

    // Indices for parents
    int parent1;
    int parent2;

    // Generations
    int generation = 0;
    fprintf(logfile_, "\nGeneration %i\n", generation);

    // Initialize the population and compute fitness
    fprintf(logfile_, "Initializing individuals and computing initial fitness...\n");
    for (i = 0; i < oldPop_.size(); i++)
    {
        initializer_->initialize(&(oldPop_[i]));
        fitComputer_->compute(oldPop_[i], Target::Train);
    }
    fprintFitness();

    // FIXME: THIS IS NOT GENERAL. Open files of each individual
    for (Individual *ind : oldPop)
    {
        ACMIndividual *tmpInd = static_cast<ACMIndividual*>(ind);
        tmpInd->openParamConvFiles(oenv_);
        tmpInd->openChi2ConvFile(oenv_, bch_->evaluateTestset());
    }

    // Copy individuals into newPop_
    for (i = 0; i < oldPop_.size(); i++) newPop_[i] = oldPop_[i]->clone();

    // Initialize best individual
    bestInd_ = oldPop_[findBestIndex()]->clone();

    fprintPop();
    fprintBestIndInPop();
    fprintBestInd();

    // Iterate and create new generation
    do
    {

        // Increase generation counter
        generation++;
        fprintf(logfile_, "\nGeneration %i\n", generation);

        // Sort individuals based on fitness
        fprintf(logfile_, "Sorting... (if needed)\n");
        sorter_->sort(&oldPop_);
        fprintPop();

        // Normalize the fitness into a probability
        fprintf(logfile_, "Computing probabilities...\n");
        probComputer_->compute(&oldPop_);
        fprintProbability();

        // Move the "nElites" best individuals (unchanged) into the new population (assuming population is sorted)
        fprintf(logfile_, "Moving the %i best individual(s) into the new population...\n", gach_->nElites());
        for (i = 0; i < gach_->nElites(); i++) newPop_[i]->copyGenome(oldPop_[i]);

        // Generate new population after the elitism
        fprintf(logfile_, "Generating the rest of the new population...\n");
        for (i = gach_->nElites(); i < oldPop_.size(); i += 2)
        {
            fprintf(logfile_, "i = %i, %i\n", i, i + 1);

            // Select parents
            parent1 = selector_->select(oldPop_);
            parent2 = selector_->select(oldPop_);
            fprintf("parent1: %i; parent2: %i\n", parent1, parent2);

            // Do crossover
            fprintf(logfile_, "Before crossover\n");
            fprintf(logfile_, "Parent 1:\n");
            oldPop_[parent1]->fprintSelf(logfile_);
            fprintf(logfile_, "Parent 2:\n");
            oldPop_[parent2]->fprintSelf(logfile_);
            if (dis(gen) <= gach_->prCross())  // If crossover is to be performed
            {
                fprintf(logfile_, "Doing crossover...\n");
                crossover_->offspring(oldPop_[parent1], oldPop_[parent2], newPop_[i], newPop_[i+1]);
            }
            else
            {
                fprintf(logfile_, "Omitting crossover...\n");
                newPop_[i]->copyGenome(oldPop_[parent1]);
                newPop_[i+1]->copyGenome(oldPop_[parent2]);
            }
            fprintf(logfile_, "Child 1:\n");
            newPop_[i]->fprintSelf(logfile_);
            fprintf(logfile_, "Child 2:\n");
            newPop_[i+1]->fprintSelf(logfile_);

            // Do mutation in each child
            fprintf(logfile_, "Doing mutation...\n");
            for (k = 0; k < 2; k++)
            {
                mutator_->mutate(newPop_[i + k], gach_->prMut());
            }
            fprintf(logfile_, "Child 1:\n");
            newPop_[i]->fprintSelf(logfile_);
            fprintf(logfile_, "Child 2:\n");
            newPop_[i+1]->fprintSelf(logfile_);

        }

        // Swap oldPop and newPop
        fprintf(logfile_, "Swapping oldPop and newPop...\n");
        tmpPop_ = oldPop_;
        oldPop_ = newPop_;
        newPop_ = tmpPop_;

        // Compute fitness
        fprintf(logfile_, "Computing fitness of new generation...\n");
        for (i = 0; i < oldPop_.size(); i++)
        {
            fitComputer_->compute(oldPop_[i], Target::Train);
        }
        fprintFitness();

        fprintPop();
        fprintBestIndInPop();

        // Check if a better individual was found, and update if so
        Individual *tmpBest = oldPop_[findBestIndex()];
        if (tmpBest->fitnessTrain() < bestInd_->fitnessTrain())
        {
            fprintf(logfile_, "A new best individual has been found!\n");
            fprintf(logfile_, "Previous best:\n");
            bestInd_->fprintSelf(logfile_);
            fprintf(logfile_, "New best:\n");
            bestInd_ = tmpBest->clone();
        }

        fprintf(logfile_, "Checking termination conditions...\n");

    } while (!terminator_->terminate(oldPop_, generation));

    // FIXME: THIS IS NOT GENERAL. Close files of each individual
    for (Individual *ind : oldPop_) static_cast<ACMIndividual*>(ind)->closeConvFiles();

    // Close surveillance files for fitness
    fclose(fileFitnessTrain_); fileFitnessTrain_ = NULL;
    fclose(fileFitnessTest_); fileFitnessTest_ = NULL;

    fprintf(logfile_, "\nGA/HYBRID Evolution is done!\n");
    fprintBestInd();

}

void GeneticAlgorithm::evolve()
{

    if (strcmp(gach_->optimizer(), "MCMC") == 0)
    {
        evolveMCMC();
    }
    else
    {
        evolveGA();
    }

}


}  //namespace ga
