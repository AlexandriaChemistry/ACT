/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "GeneticAlgorithm.h"

#include <cstdio>

#include "Initializer.h"
#include "FitnessComputer.h"
#include "ProbabilityComputer.h"
#include "Crossover.h"
#include "Mutator.h"
#include "Terminator.h"

#include "gromacs/fileio/gmxfio.h"

namespace ga
{

GeneticAlgorithm::GeneticAlgorithm(FILE                                *logFile,
                                   struct gmx_output_env_t             *oenv,
                                   Initializer                         *initializer,
                                   FitnessComputer                     *fitnessComputer,
                                   Sorter                              *sorter,
                                   ProbabilityComputer                 *probComputer,
                                   Selector                            *selector,
                                   Crossover                           *crossover,
                                   Mutator                             *mutator,
                                   Terminator                          *terminator,
                                   int                                  popSize,
                                   bool                                 evaluateTestSet)
    : popSize_(popSize), evaluateTestSet_(evaluateTestSet), logfile_(logFile), oenv_(oenv),
      newPop_(popSize),
      initializer_(initializer),
      fitComputer_(fitnessComputer), sorter_(sorter), probComputer_(probComputer),
      selector_(selector), crossover_(crossover), mutator_(mutator), terminator_(terminator)
{
    // Create directories for each individual
    for (int i = 0; i <= popSize; i++)
    {
        const std::string command = "mkdir ind" + std::to_string(i);
        system(command.c_str());
        // std::filesystem::create_directory(dirName.c_str());
    }
}

void GeneticAlgorithm::setBestIndividual(Individual *ind)
{
    bestInd_ = ind;
}

void GeneticAlgorithm::fprintPop() const
{
    if (logfile_)
    {
        fprintf(logfile_, "Population:\n");
        for (Individual *ind : oldPop_)
        {
            ind->fprintSelf(logfile_);
        }
    }
}

void GeneticAlgorithm::fprintBestInd() const
{
    if (logfile_)
    {
        fprintf(logfile_, "Overall Best Individual:\n");
        bestInd_->fprintSelf(logfile_);
    }
}

void GeneticAlgorithm::fprintBestIndInPop() const
{
    if (logfile_)
    {
        fprintf(logfile_, "Best Individual in current population:\n");
        oldPop_[findBestIndex()]->fprintSelf(logfile_);
    }
}

int GeneticAlgorithm::findBestIndex() const
{
    int index = 0;
    double bestFitness = oldPop_[index]->fitnessTrain();
    for (size_t i = 1; i < oldPop_.size(); i++)
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
    if (logfile_)
    {
        fprintf(logfile_, "Probability: [ ");
        for (Individual *ind : oldPop_)
        {
            fprintf(logfile_, "%f ", ind->probability());
        }
        fprintf(logfile_, "]\n");
    }
}

void GeneticAlgorithm::fprintFitness() const
{
    for (size_t i = 0; i < oldPop_.size() - 1; i++)
    {
        fprintf(fileFitnessTrain_, "%lf ", oldPop_[i]->fitnessTrain());
        fprintf(fileFitnessTest_, "%lf ", oldPop_[i]->fitnessTest());
    }
    fprintf(fileFitnessTrain_, "%lf\n", oldPop_[oldPop_.size()-1]->fitnessTrain());
    fprintf(fileFitnessTest_, "%lf\n", oldPop_[oldPop_.size()-1]->fitnessTest());
}

void GeneticAlgorithm::swapOldNewPopulations()
{
    auto tmpPop = oldPop_;
    oldPop_ = newPop_;
    newPop_ = tmpPop;
}

void GeneticAlgorithm::copyOldToNewPopulations()
{
    newPop_.clear();
    for (int i = 0; i < populationSize(); i++)
    {
        newPop_.push_back(oldPop_[i]->clone());
    }
}

void GeneticAlgorithm::openFitnessFiles()
{
    fileFitnessTrain_ = gmx_fio_fopen("ga_fitness_train.txt", "w");
    fileFitnessTest_  = gmx_fio_fopen("ga_fitness_test.txt", "w");
    GMX_RELEASE_ASSERT(fileFitnessTrain_ != NULL && fileFitnessTest_ != NULL, "Could not open files");
}

void GeneticAlgorithm::closeFitnessFiles()
{
    gmx_fio_fclose(fileFitnessTrain_);
    gmx_fio_fclose(fileFitnessTest_);
}

}  //namespace ga
