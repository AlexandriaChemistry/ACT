#ifndef ACTMOL_UTIL_H
#define ACTMOL_UTIL_H
    
#include <vector>

struct t_inputrec;

namespace alexandria
{
    class  ForceField;
    class  ForceComputer;
    class  ACTMol;
    
    /*! \brief Read a file from the test directories and produce a vector of actmols.
     * \param[in]  molname  The name of the compound
     * \param[in]  pd       The force field
     * \param[in]  fcomp    The force computer
     * \param[in]  inputrec GROMACS input record
     * \param[out] mps      The ACTMol structures
     */
    void initACTMol(const char          *molname, 
                    const ForceField    *pd,
                    ForceComputer       *fcomp,
                    t_inputrec          *inputrec,
                    std::vector<ACTMol> *mps);

}

#endif
