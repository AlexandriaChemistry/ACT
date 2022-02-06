/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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
#ifndef ACT_GA_GENOME_H
#define ACT_GA_GENOME_H

#include <map>
#include <vector>

#include "act/utility/communicationrecord.h"

#include "Dataset.h"

typedef std::map<iMolSelect, double> FitnessMap;

namespace ga
{

class Genome
{
private:
    //! The values of each "base"
    std::vector<double> genome_;
    //! And the fitness corresponding to this genome
    FitnessMap          fitness_;
    //! Probability of training selection
    double              probability_ = 0.0;
public:
    //! Constructor
    Genome(const std::vector<double> genome,
           FitnessMap                fitness) : 
    genome_(genome), fitness_(fitness) {}
    
    //! Empty constructor
    Genome() {}
    
    //! \return whether this fitness term is present
    bool hasFitness(iMolSelect ims) const { return fitness_.end() != fitness_.find(ims); }
    
    /*! Return the fitness in the data set
     * \param[in] ims The data set
     * \return the fitness
     * \throws if not found
     */
    double fitness(iMolSelect ims) const;

    /*! Return a pointer to the fitness in the data set
     * \param[in] ims The data set
     * \return the fitness
     * \throws if not found
     */
    double *fitnessPtr(iMolSelect ims) { return &(fitness_.find(ims)->second); }

    /*! \brief Set a fitness
     * \param[in] ims     The data set
     * \param[in] fitness The fitness term
     */
    void setFitness(iMolSelect ims, double fitness);
    
    /*! \brief Unset a fitness, i.e. remove it from the map.
     * \param[in] ims     The data set
     */
    void unsetFitness(iMolSelect ims);
    
    //! \return the probability of the training set
    double probability() const { return probability_; }
    
    /*! \brief Set the probability of the training set
     * \param[in] p The probability
     */
    void setProbability(double p) { probability_ = p; }
    
    /*! \brief Add a base to the genome
     * \param[in] base The value
     */
    void addBase(double base) { genome_.push_back(base); }

    /*! \brief Set one base
     * \param[in] index The position in the genome
     * \param[in] value The new value
     */
    void setBase(size_t index, double value);
    
    /*! \brief Get the value of one base
     * \param[in] index The position in the gene
     * \return the value
     */
    double base(size_t index) const;

    //! \return lenght of the genome
    size_t nBase() const { return genome_.size(); }
    
    /*! \brief Get the values of the whole genome
     * \return the value vector
     */
    const std::vector<double> &bases() const { return genome_; }
    
    /*! \brief Get the mutable values of the whole genome
     * \return the value vector
     */
    std::vector<double> *basesPtr() { return &genome_; }
    
    /*! \brief Print my content
     * \param[in] fp File pointer to print to
     */
    void print(FILE *fp) const;
    
    /*! \brief Print my content
     * \param[in] name String to start the line with
     * \param[in] fp   File pointer to print to. If nullptr nothing is printed.
     */
    void print(const char *name, FILE *fp) const;
    
    /*! \brief Send to another processor
     * \param[in] cr   The communication record
     * \param[in] dest The destination processor
     */
    void Send(const alexandria::CommunicationRecord *cr, int dest) const;
    
    /*! \brief Receive from another processor
     * \param[in] cr  The communication record
     * \param[in] src The source processor
     */
    void Receive(const alexandria::CommunicationRecord *cr, int src);
};

} // namespace ga

#endif
