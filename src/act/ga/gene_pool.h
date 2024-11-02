﻿/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
#ifndef ACT_GA_GENEPOOL_H
#define ACT_GA_GENEPOOL_H

#include <algorithm>
#include <vector>

#include "gromacs/utility/gmxassert.h"

#include "genome.h"

namespace ga
{

/*! \brief Class that holds all the genomes for the whole population
 */
class GenePool
{
private:
    //! The size of each genome, there cannot be deletions or insertions
    size_t              genomeSize_ = 0;
    //! The gene pool is just a vector of Genomes
    std::vector<Genome> genomes_;
public:
    //! Constructor
    GenePool(size_t genomeSize) : genomeSize_(genomeSize) {}
    
    //! Default constructor
    GenePool() {}

    //! \return the genome size
    size_t genomeSize() const { return genomeSize_; } 

    /*! Add a genome
     * \param[in] genome  The parameter values
     * \param[in] fitness The fitness value (optiona)
     * \throws if genome does not have the expected length
     */
    void addGenome(const std::vector<double> &genome, double fitness = 0.0);
    
    /*! Add a genome
     * \param[in] genome  The complete genome
     * \throws if genome does not have the expected length
     */
    void addGenome(const Genome &genome);
    
    //! Replace a genome
    void replaceGenome(size_t index, const Genome &genome) { genomes_[index] = genome; }
    
    //! Return a genome
    const Genome &genome(size_t index) const
    {
        GMX_RELEASE_ASSERT(index < genomes_.size(), "No such genome");
        return genomes_[index]; 
    }
    
    //! Return a mutable genome
    Genome *genomePtr(size_t index) { return &genomes_[index]; }
    
    //! Return the gene pool
    const std::vector<Genome> &genePool() const { return genomes_; }
    
    //! Return the mutable gene pool
    std::vector<Genome> *genePoolPtr() { return &genomes_; }

    /*! \brief Read a genepool from a space and new-line separated file
     * \param[in] fileName The name of the file to read from
     */
    void read(const std::string &fileName);

    /*! \brief Write a genepool to a space and new-line separated file
     * \param[in] fileName The name of the file to write to
     */
    void write(const std::string &fileName) const;

    //! Return the population size
    size_t popSize() const { return genomes_.size(); }
    
    /*! \brief Sort the GenePool according to a data set
     * \param[in] ims The data set to use for the fitness
     */
    void sort(iMolSelect ims);
    
    /*! \return the index of the genome with the lowest fitness.
     * \param[in] ims the dataset to grab the fitness from
     */
    size_t findBestIndex(const iMolSelect ims) const;

    /*! \return a const referene of the genome with the lowest fitness.
     * \param[in] ims the dataset to grab the fitness from
     */
    const Genome &getBest(const iMolSelect ims) const;

    /*! Print the Gene Pool
     * \param[in] fp The file to print to
     */
    void print(FILE *fp) const;
    
    //! \return the minimum value per parameter in this pool
    std::vector<double> min() const;

    //! \return the maximum value per parameter in this pool
    std::vector<double> max() const;

    //! \return the mean value per parameter in this pool
    std::vector<double> mean() const;

    /*!
     * \param[in] mean the mean per parameter in this pool
     * \return the standard deviation per parameter in this pool
     */ 
    std::vector<double> stdev(const std::vector<double> &mean) const;

    /*!
     * \return the standard deviation per parameter in this pool
     */ 
    std::vector<double> stdev() const { return stdev(mean()); }

    //! \return the median value per parameter in this pool
    std::vector<double> median() const;

};

} // namespace ga

#endif
