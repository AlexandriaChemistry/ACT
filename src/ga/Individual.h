/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julianramon.marradesfurquet.8049@student.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef GA_INDIVIDUAL_H
#define GA_INDIVIDUAL_H


namespace ga
{


//! Abstract individual for genetic algorithms!
class Individual
{

protected:

    //! Fitness for training set
    double fitnessTrain_ = 0.0;
    //! Fitness for test set
    double fitnessTest_ = 0.0;
    //! Probability of selection
    double probability_ = 0.0;

    //! Default constructor FIXME: maybe we have to make it public
    Individual() {}

    /*!
     * Property constructor
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
    * BEGIN: Getters and Setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * Get the fitness of the individual in training set
     * @returns the fitness
     */
    double fitnessTrain() const { return fitnessTrain_; }

    /*!
     * Get a pointer to \p fitnessTrain_
     * @returns the pointer
     */
    double *fitnessTrainPtr() { return &fitnessTrain_; }

    /*!
     * Set the fitness of the individual in training set
     * @param fitnessTrain the fitness
     */
    void setFitnessTrain(const double fitnessTrain) { fitnessTrain_ = fitnessTrain; }

    /*!
     * Get the fitness of the individual in test set
     * @returns the fitness
     */
    double fitnessTest() const { return fitnessTest_; }

    /*!
     * Get a pointer to \p fitnessTest_
     * @returns the pointer
     */
    double *fitnessTestPtr() { return &fitnessTest_; }

    /*!
     * Set the fitness of the individual in test set
     * @param fitnessTest the fitness
     */
    void setFitnessTest(const double fitnessTest) { fitnessTest_ = fitnessTest; }

    /*!
     * Get the selection probability of the individual
     * @returns the selection probability
     */
    double probability() const { return probability_; }

    /*!
     * Get a pointer to \p probability_
     * @returns the pointer
     */
    double *probabilityPtr() { return &probability_; }

    /*!
     * Set the selection probability
     * @param probability the selection probability
     */
    void setProbability(const double probability) { probability_ = probability; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace ga


#endif //GA_INDIVIDUAL_H