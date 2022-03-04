#include <string>
#include <vector>

#include <cstdio>

#include "act/molprop/molpropobservable.h"

namespace alexandria
{

extern std::array<MolPropObservable, 4> mpoMultiPoles;

/*! \brief Return index in array
 * \param[in] id The string, e.g. "xxx" will be interpreted both for its length
 *               and its content. Length will tell which moment is requested.
 * \param[out] index The index in the values array
 * \return true if found, false if not.
 */
bool multipoleIndex(const std::string &id, int *index);
    
/*! \brief Make a pretty formatted multipole
 * \param[in] mpo    Should be a multipole
 * \param[in] values The numerical components
 * \return a vector of strings
 */
std::vector<std::string> formatMultipole(MolPropObservable          mpo,
                                         const std::vector<double> &values);

/*! \brief Print a pretty formatted multipole
 * \param[in] fp     File pointer to print to
 * \param[in] mpo    Should be a multipole
 * \param[in] values The numerical components
 */
void printMultipole(FILE                      *fp,
                    MolPropObservable          mpo,
                    const std::vector<double> &values);

/*! \brief Return the name of a multipole index.
 * The length of m determines the order of the multipole:
 * 1 for dipole, 2 for quadrupole, etc.
 * \param[in] m The indices
 * \return string such as e.g. "xyz" for an octupole term
 */
const std::string &multipoleName(const std::vector<int> &m);

/*! \brief Return the names of a multipole category
 * \param[in] mpo The type of multipole
 * \return string such as e.g. "xyz" for an octupole term
 */
const std::vector<std::string> &multipoleNames(MolPropObservable mpo);

/*! \brief Return the index of an multipole component
 * \param[in] m See above
 * \return the index
 * \throws if input incorrect
 */
int multipoleIndex(const std::vector<int> &m);

} // namespace alexandria
