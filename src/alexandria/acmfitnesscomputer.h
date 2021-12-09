#ifndef ALEXANDRIA_ACMFITNESSCOMPUTER_H
#define ALEXANDRIA_ACMFITNESSCOMPUTER_H


#include <vector>

#include "devcomputer.h"
#include "aliases.h"
#include "acmindividual.h"


namespace alexandria
{


class ACMFitnessComputer : public ga::FitnessComputer
{

private: 

    //! \brief The filepointer to the log file.
    FILE *logfile_;

    //! \brief A pointer to the BoundsDevComputer.
    BoundsDevComputer *bdc_;

    //! \brief A vector of devComputers.
    std::vector<DevComputer*> devComputers_;

    //! \brief SharedIndividualInfo pointer
    SharedIndividualInfo *sii_;

    //! \brief MolGen pointer
    MolGen *mg_;

    //! \brief Whether or not to remove molecules that fail to converge in the shell minimization
    bool removeMol_;

    //! \brief Whether we are in verbose mode or not
    bool verbose_;

    //! \brief Amount of times calcDeviation() has been called
    int numberCalcDevCalled_ = 0;

    /*! \brief Compute dipole and quadrupole moments (if needed), for a given molecule
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param mymol     the molecule
     */
    void computeDiQuad(std::map<eRMS, FittingTarget> *targets,
                       MyMol                         *mymol);

public:

    //! \brief The flag which determines the amount of explaining text for the output.
    bool verbose_;

    /*! \brief Compute the desired entities.
     * @param[in] ind    The pointer to the individual to compute for
     */
    virtual void compute(ga::Individual *ind);

    /*! \brief Computes deviation from target
     * \param[in] ind    pointer to individual
     * \param[in] verbose       Whether or not to print a lot
     * \param[in] cr            communication record
     * \param[in] calcDev       The type of calculation to do
     * \param[in] ims           The data set to do computations on
     * \return the square deviation
     */
    double calcDeviation(      ACMIndividual   *ind,
                               t_commrec       *cr,
                         const bool             verbose,
                               CalcDev          calcDev,
                               iMolSelect       ims);

    //! \brief Fill the devComputers vector according to the needs of the user
    void fillDevComputers();  // TODO: maybe call this in the constructor?
};


} // namespace alexandria


#endif // ALEXANDRIA_ACMFITNESSCOMPUTER_H