/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef ALEXANDRIA_ACMFITNESSCOMPUTER_H
#define ALEXANDRIA_ACMFITNESSCOMPUTER_H

#include <vector>

#include "act/ga//FitnessComputer.h"

#include "acmindividual.h"
#include "bayes.h"
#include "devcomputer.h"
#include "molgen.h"

namespace alexandria
{


/*!
 * \brief Computes \f$ \chi^2 \f$ of an individual as its fitness
 */
class ACMFitnessComputer : public ga::FitnessComputer
{

private: 
    //! \brief The filepointer to the log file.
    FILE *logfile_;
    //! \brief A pointer to the BoundsDevComputer.
    BoundsDevComputer         *bdc_  = nullptr;
    //! \brief A vector of devComputers.
    std::vector<DevComputer*>  devComputers_;
    //! \brief StaticIndividualInfo pointer
    StaticIndividualInfo      *sii_ = nullptr;
    //! \brief MolGen pointer
    MolGen *molgen_;
    //! \brief Whether or not to remove molecules that fail to converge in the shell minimization
    bool removeMol_;
    //! \brief Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
    bool verbose_;
    //! \brief Amount of times calcDeviation() has been called
    int numberCalcDevCalled_ = 0;

    /*! \brief Compute multipole moments (if needed), for a given molecule
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param mymol     the molecule
     */
    void computeMultipoles(std::map<eRMS, FittingTarget> *targets,
                           MyMol                         *mymol);

    //! \brief Fill the devComputers vector according to the needs of the user
    void fillDevComputers();

public:

    /*!
     * Constructor
     * \param[in] logfile           pointer to logfile
     * \param[in] sii               pointer to StaticIndividualInfo
     * \param[in] mg                pointer to molgen
     * \param[in] removeMol         Whether or not to remove molecules that fail to converge in the shell minimization
     * \param[in] verbose           Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
     */
    ACMFitnessComputer(      FILE                  *logfile,
                             StaticIndividualInfo  *sii,
                             MolGen                *molgen,
                       const bool                   removeMol,
                       const bool                   verbose)
    : logfile_(logfile), sii_(sii), molgen_(molgen), removeMol_(removeMol), verbose_(verbose)
    {
        fillDevComputers();
    }

    virtual void compute(ga::Genome *genome,
                         iMolSelect  trgtFit
                         bool        verbose) override;

    /*! \brief Computes deviation from target
     * \param[in] params   The force field parameters
     * \param[in] calcDev  The type of calculation to do
     * \param[in] ims      The dataset to do computations on
     * \return the square deviation
     */
    double calcDeviation(std::vector<double> *params,
                         CalcDev              calcDev,
                         iMolSelect           ims);

};


} // namespace alexandria


#endif // ALEXANDRIA_ACMFITNESSCOMPUTER_H
