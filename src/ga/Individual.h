/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_INDIVIDUAL_H
#define GA_INDIVIDUAL_H


#include <cstdio>


namespace ga
{


//! \brief Abstract individual for genetic algorithms
class Individual
{

protected:

    //! Fitness for training set
    double fitnessTrain_ = 0.0;
    //! Fitness for test set
    double fitnessTest_ = 0.0;
    //! Probability of selection
    double probability_ = 0.0;

    //! \brief Default constructor
    Individual() {}

    /*!
     * \brief Property constructor
     * \param[in] fitnessTrain  the fitness for training set
     * \param[in] fitnessTest   the fitness for test set
     * \param[in] probability   the probability of selection
     */
    Individual(const double fitnessTrain,
               const double fitnessTest,
               const double probability)
    : fitnessTrain_(fitnessTrain), fitnessTest_(fitnessTest), probability_(probability) {}

public:

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Printing                          *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Print information about this Individual to file
     * \param[in] fp file pointer
     */
    virtual void fprintSelf(FILE *fp) = 0;

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Printing                            *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Cloning                           *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return a copy of this Individual
    virtual Individual *clone() = 0;

    /*!
     * \brief Copy the genome from another Individual
     * \param[in] other pointer to another Individual
     */
    virtual void copyGenome(Individual *other) = 0;

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Cloning                             *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and Setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return the fitness in training set
    double fitnessTrain() const { return fitnessTrain_; }

    //! \return a pointer to the fitness in training set
    double *fitnessTrainPtr() { return &fitnessTrain_; }

    /*!
     * \brief Set the fitness in training set
     * \param[in] fitnessTrain the fitness
     */
    void setFitnessTrain(const double fitnessTrain) { fitnessTrain_ = fitnessTrain; }

    //! \return the fitness in test set
    double fitnessTest() const { return fitnessTest_; }

    //! \return a pointer to the fitness in test set
    double *fitnessTestPtr() { return &fitnessTest_; }

    /*!
     * \brief Set the fitness in test set
     * \param[in] fitnessTest the fitness
     */
    void setFitnessTest(const double fitnessTest) { fitnessTest_ = fitnessTest; }

    //! \return the selection probability
    double probability() const { return probability_; }

    //! \return a pointer to the selection probability
    double *probabilityPtr() { return &probability_; }

    /*!
     * \brief Set the selection probability
     * \param[in] probability the selection probability
     */
    void setProbability(const double probability) { probability_ = probability; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace ga


#endif //GA_INDIVIDUAL_H