/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "GeneticAlgorithm.h"

#include <cstdio>

#include "Crossover.h"
#include "FitnessComputer.h"
#include "Initializer.h"
#include "Mutator.h"
#include "ProbabilityComputer.h"
#include "Terminator.h"

#include "act/basics/dataset.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace ga
{

void GeneticAlgorithm::openFitnessFiles()
{
    for(const auto &im : iMolSelectNames())
    {
        std::string fn = gmx::formatString("ga_fitness_%s.txt", im.second);
        fileFitness_.insert({im.first, gmx_fio_fopen(fn.c_str(), "w")});
        GMX_RELEASE_ASSERT(fileFitness_[im.first] != NULL, "Could not open file");
    }
}

void GeneticAlgorithm::closeFitnessFiles()
{
    for(const auto &ff : fileFitness_)
    {
        gmx_fio_fclose(ff.second);
    }
}

void GeneticAlgorithm::fprintFitness(const GenePool &pool)
{
    for (const auto &ff : fileFitness_)
    {
        for (const auto &genome : pool.genePool())
        {
            if (genome.hasFitness(ff.first))
            {
                fprintf(ff.second, "%f ", genome.fitness(ff.first));
            }
        }
        fprintf(ff.second, "\n");
    }
}

}  //namespace ga
