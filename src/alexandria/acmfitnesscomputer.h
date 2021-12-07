#ifndef ALEXANDRIA_ACMFITNESSCOMPUTER_H
#define ALEXANDRIA_ACMFITNESSCOMPUTER_H

#include "devcomputer.h"
#include "aliases.h"


namespace alexandria
{


class ACMFitnessComputer : public ga::FitnessComputer
{

private: 

    //! \brief The filepointer to the log file.
    FILE *logfile_;

    //! \brief A pointer to the BoundsDevComputer.
    BoundsDecComputer* bdc_;

    //! \brief A vector of devComputers.
    vector<DevComputer> devComputers_;

    /*! \brief Computes the dipole and quadropole for a given molecule.
     * @param[in] targets   Computates the dipole and quadropole for a given molecule.
     * @param[in] mymol     A pointer to the molecule.
     */
    void computeDiQuad(std::map<eRMS, FittingTarget> *targets,
                       MyMol                         *mymol);

public:

    //! \brief The flag which determines the amount of explaining text for the output.
    bool verbose_;

    /*! \brief Compute the desired entities.
     * @param[in] individual    The pointer to the individual to compute for
     */
    virtual void compute(Individual *individual);

    /*! \brief Computes deviation from target
     * \param[in] individual    The pointer to the individual to initialize
     * \param[in] calcDev       The type of calculation to do
     * \param[in] ims           The data set to do computations on
     * \return the square deviation
     */
    double calcDeviation(Individual   *individual,
                         CalcDev      calcDev,
                         iMolSelect   ims);

    //! \brief Fill the devComputers vector according to the needs of the user.
    void fillDevComputers();
};


} // namespace alexandria


#endif // ALEXANDRIA_ACMFITNESSCOMPUTER_H