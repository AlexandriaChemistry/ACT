#ifndef ALEXANDRIA_PARAMHANDLER_H
#define ALEXANDRIA_PARAMHANDLER_H

#include <vector>

#include "gromacs/commandline/pargs.h"

namespace alexandria
{

/*!
 * Abstract class to handle parameters for methods.
 * It adds its command-line arguments to a given parameter vector
 * and later checks its validity.
 */
class ParamHandler
{

public:

    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     */
    virtual void add_pargs(std::vector<t_pargs> *pargs) = 0;

    /*!
     * \brief Check the validity of the provided arguments
     * Throw an exception if an invalid combination of arguments is encountered
     */
    virtual void check_pargs() = 0;

};

}

#endif //ALEXANDRIA_PARAMHANDLER_H