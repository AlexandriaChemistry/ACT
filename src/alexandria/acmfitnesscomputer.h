/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julianramon.marradesfurquet.8049@student.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef ALEXANDRIA_ACMFITNESSCOMPUTER_H
#define ALEXANDRIA_ACMFITNESSCOMPUTER_H


#include <vector>

#include "devcomputer.h"
#include "acmindividual.h"
#include "ga/FitnessComputer.h"
#include "bayes.h"


namespace alexandria
{


class ACMFitnessComputer : public ga::FitnessComputer
{

private: 

    //! Communications record
    t_commrec *cr_;
    //! \brief The filepointer to the log file.
    FILE *logfile_;
    //! \brief A pointer to the BoundsDevComputer.
    BoundsDevComputer *bdc_ = nullptr;
    //! \brief A vector of devComputers.
    std::vector<DevComputer*> devComputers_;
    //! \brief SharedIndividualInfo pointer
    SharedIndividualInfo *sii_;
    //! \brief MolGen pointer
    MolGen *mg_;
    //! \brief Whether or not to remove molecules that fail to converge in the shell minimization
    bool removeMol_;
    //! \brief Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
    bool verbose_;
    //! Whether we consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization
    bool fullQuadrupole_;
    //! \brief Amount of times calcDeviation() has been called
    int numberCalcDevCalled_ = 0;

    /*! \brief Compute dipole and quadrupole moments (if needed), for a given molecule
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param mymol     the molecule
     */
    void computeDiQuad(std::map<eRMS, FittingTarget> *targets,
                       MyMol                         *mymol);

    //! \brief Fill the devComputers vector according to the needs of the user
    void fillDevComputers();

public:

    /*!
     * Constructor
     * \param[in] cr                communcations record
     * \param[in] logfile           pointer to logfile
     * \param[in] sii               pointer to SharedIndividualInfo
     * \param[in] mg                pointer to molgen
     * \param[in] removeMol         Whether or not to remove molecules that fail to converge in the shell minimization
     * \param[in] verbose           Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
     * \param[in] fullQuadrupole    Whether we consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization
     */
    ACMFitnessComputer(      t_commrec             *cr,
                             FILE                  *logfile,
                             SharedIndividualInfo  *sii,
                             MolGen                *mg,
                       const bool                   removeMol,
                       const bool                   verbose,
                       const bool                   fullQuadrupole)
    : cr_(cr), logfile_(logfile), sii_(sii), mg_(mg), removeMol_(removeMol), verbose_(verbose), fullQuadrupole_(fullQuadrupole)
    {
        fillDevComputers();
    }

    virtual void compute(ga::Individual *ind,
                         ga::Target trgtFit);

    /*! \brief Computes deviation from target
     * \param[in] ind           pointer to individual
     * \param[in] verbose       Whether or not to print a lot (for when this gets called from outside the compute() routine)
     * \param[in] calcDev       The type of calculation to do
     * \param[in] ims           The data set to do computations on
     * \return the square deviation
     */
    double calcDeviation(      ACMIndividual   *ind,
                         const bool             verbose,
                               CalcDev          calcDev,
                               iMolSelect       ims);

};


} // namespace alexandria


#endif // ALEXANDRIA_ACMFITNESSCOMPUTER_H