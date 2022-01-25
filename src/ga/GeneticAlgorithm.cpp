/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "GeneticAlgorithm.h"

#include <cstdio>

#include "Crossover.h"
#include "FitnessComputer.h"
//#include "GenePool.h"
#include "Initializer.h"
#include "Mutator.h"
#include "ProbabilityComputer.h"
#include "Terminator.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/exceptions.h"

namespace ga
{

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

void GeneticAlgorithm::fprintFitness(const GenePool &pool)
{
    for (const auto &genome : pool.genePool())
    {
        fprintf(fileFitnessTrain_, "%f ", genome.fitness(iMolSelect::Train));
        fprintf(fileFitnessTest_, "%f ", genome.fitness(iMolSelect::Test));
    }
    fprintf(fileFitnessTrain_, "\n");
    fprintf(fileFitnessTest_, "\n");
}

}  //namespace ga
