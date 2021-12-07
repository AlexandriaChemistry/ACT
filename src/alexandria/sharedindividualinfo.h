#ifndef ALEXANDRIA_SHAREDINDIVIDUALINFO_H
#define ALEXANDRIA_SHAREDINDIVIDUALINFO_H


#include <vector>
#include <string>


namespace alexandria
{


class SharedIndividualInfo
{

private:

    //! Training steps per parameter
    std::vector<int> ntrain_;
    //! Lower bound per parameter
    std::vector<double> lowerBound_;
    //! Upper bound per parameter
    std::vector<double> upperBound_;
    //! Mutability per parameter
    std::vector<Mutability> mutability_; 
    //! Weighted temperature
    std::vector<double> weightedTemperature_;
    //! Optimization index for each parameter
    std::vector<OptimizationIndex> optIndex_;
    //! Class-name for each parameter
    std::vector<std::string> paramNames_;
    //! Classes in the system
    std::vector<std::string> paramClass_;
    //! Class-index for each parameter
    std::vector<std::string> paramClassIndex_;
    //! Base name for parameter convergence files
    std::string xgvconv_;
    //! Base name for Chi2 convergence file
    std::string xvgepot_;
    //! Base name for Force field output file
    std::string outputFile_;

public:

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and setters stuff         *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Get a constant \p optIndex_ reference
    const std::vector<OptimizationIndex> &optIndex() const { return optIndex_; }
    
    //! \brief Get a constant \p paramNames_ reference
    const std::vector<std::string> &paramNames() const { return paramNames_; }

    //! \brief Get a constant \p paramClass_ reference
    const std::vector<std::string> &paramClass() const { return paramClass_; }

    //! \brief Return xvg file for convergence information
    const std::string &xvgConv() const { return xvgconv_; }

    //! \brief Return xvg file for epot information
    const std::string &xvgEpot() const { return xvgepot_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and setters stuff           *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace alexandria


#endif //ALEXANDRIA_SHAREDINDIVIDUALINFO_H