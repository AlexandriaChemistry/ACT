#ifndef ALEXANDRIA_SHAREDINDIVIDUALINFO_H
#define ALEXANDRIA_SHAREDINDIVIDUALINFO_H


#include <vector>


namespace alexandria
{


class SharedIndividualInfo
{

private:

    //! Optimization index for each parameter
    std::vector<OptimizationIndex> optIndex_;

public:

    //! \brief Get a constant \p optIndex_ reference
    const std::vector<OptimizationIndex> &optIndex() const { return optIndex_; }

};


} //namespace alexandria


#endif //ALEXANDRIA_SHAREDINDIVIDUALINFO_H