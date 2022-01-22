#ifndef GA_DATASET_H
#define GA_DATASET_H

#include <map>
#include <string>

//! Class to distinguish data sets
enum class iMolSelect {
    //! Training data set
    Train,
    //! Testing data set
    Test,
    //! Data set to ignore
    Ignore
};

/*! \brief Return string corresponding to data set
 * \param[in] ims The data set
 * \return a string
 */
const char *iMolSelectName(iMolSelect ims);

//! \return map of all names of data sets
const std::map<iMolSelect, const char *> &iMolSelectNames();

/*! \brief Look up data set name
 * \param[in]  name The data set name
 * \param[out] ims  The type of data set
 * \return true if found, false otherwise
 */
bool name2molselect(const std::string &name, iMolSelect *ims);

#endif
