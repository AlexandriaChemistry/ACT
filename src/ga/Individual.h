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

    //! Default constructor
    Individual() {}

    /*!
     * Property constructor
     * \param[in] fitnessTrain  the fitness for training set
     * \param[in] fitnessTest  the fitness for test set
     */
    Individual(const double fitnessTrain,
               const double fitnessTest)
    : fitnessTrain_(fitnessTrain), fitnessTest_(fitnessTest) {}

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

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace ga


#endif //GA_INDIVIDUAL_H