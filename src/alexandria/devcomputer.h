#ifndef ALEXANDRIA_DEVCOMPUTER_H
#define ALEXANDRIA_DEVCOMPUTER_H

#include <cstdio>
#include <map>
#include <vector>

#include "gromacs/mdtypes/commrec.h"

#include "mymol.h"
#include "molgen.h"


namespace alexandria
{


/*!
 * Abstract class for chi-squared deviation computation
 * Each DevComputer handles a particular component of chi-squared
 */
class DevComputer
{

protected:

    //! Pointer to log file
    FILE *logfile_;
    //! Whether we are in verbose mode
    bool verbose_;

    /*! \brief Create a new DevComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     */
    DevComputer(      FILE *logfile,
                const bool  verbose)
    {
        logfile_ = logfile;
        verbose_ = verbose;
    }

public:

    /*! \brief Computes a component of the chi-squared deviation and place it into the appropriate FittingTarget
     * @param mymol     pointer to the molecule to compute the deviation for
     * @param targets   map between different components of chi-squared and their fitting targets
     * @param poldata   pointer to Poldata structure
     * @param param     the current force field parameter vector
     * @param commrec   pointer to communications record
     */
    virtual void calcDeviation(      MyMol                             *mymol,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     Poldata                           *poldata,
                               const std::vector<double>               &param,
                                     t_commrec                         *commrec) = 0;

};


}

#endif //ALEXANDRIA_DEVCOMPUTER_H