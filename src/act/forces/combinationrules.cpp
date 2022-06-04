#include "combinationrules.h"

#include <cmath>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

void CombineLJ(int     CombinationRule,
               double  sigmaI,
               double  sigmaJ,
               double  epsilonI,
               double  epsilonJ,
               double *c6,
               double *c12)
{
    switch (CombinationRule)
    {
    case eCOMB_GEOMETRIC:
        {
            double sig  = std::sqrt(sigmaI * sigmaJ);
            double eps  = std::sqrt(epsilonI * epsilonJ);
            double sig6 = std::pow(sig, 6.0);
            *c6  = 4*eps*sig6;
            *c12 = *c6 * sig6;
        }
        break;
    case eCOMB_ARITHMETIC:
        {
            double sig  = 0.5 * (sigmaI + sigmaJ);
            double eps  = 0.5 * (epsilonI + epsilonJ);
            double sig6 = std::pow(sig, 6.0);
            *c6  = 4*eps*sig6;
            *c12 = *c6 * sig6;
        }
        break;
    case eCOMB_LORENTZ_BERTHELOT:
        {
            double sig  = 0.5 * (sigmaI + sigmaJ);
            double eps  = std::sqrt(epsilonI * epsilonJ);
            double sig6 = std::pow(sig, 6.0);
            *c6  = 4*eps*sig6;
            *c12 = *c6 * sig6;
        }
        break;
    default:
        gmx_fatal(FARGS, "Unsupported combination rule %d for Lennard Jones", CombinationRule);
    }
}

void CombineBham(int     CombinationRule,
                 double  sigmaI,
                 double  sigmaJ,
                 double  epsilonI,
                 double  epsilonJ,
                 double  gammaI,
                 double  gammaJ,
                 double *sigmaIJ,
                 double *epsilonIJ,
                 double *gammaIJ)
{
    switch (CombinationRule)
    {
        case eCOMB_GEOMETRIC:
            *sigmaIJ = std::sqrt(sigmaI * sigmaJ);
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ);
            *gammaIJ = std::sqrt(gammaI * gammaJ);
            break;
        case eCOMB_ARITHMETIC:
            *sigmaIJ = 0.5 * (sigmaI + sigmaJ);
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ);
            *gammaIJ = 0.5 * (gammaI + gammaJ);
            break;
        case eCOMB_KONG_MASON:
            *sigmaIJ = std::sqrt(sigmaI * sigmaJ);
            *epsilonIJ = 2.0 * (epsilonI * epsilonJ)/(epsilonI + epsilonJ);
            *gammaIJ = *sigmaIJ * (0.5*((gammaI/sigmaI)+(gammaJ/sigmaJ)));
            break;
        case eCOMB_HOGERVORST: // Hogervorst, Physica, Volume: 51, Page: 77, Year: 1971. Combination rules for Buckingham.
            {
                *gammaIJ = 0.5 * (gammaI + gammaJ);  
                *epsilonIJ = (2.0 * epsilonI * epsilonJ)/(epsilonI + epsilonJ);
                double itmp = (epsilonI*gammaI*std::pow(sigmaI,6.0))/(gammaI-6.0);
                double jtmp = (epsilonJ*gammaJ*std::pow(sigmaJ,6.0))/(gammaJ-6.0);
                *sigmaIJ = std::pow(std::sqrt(itmp*jtmp)*(*gammaIJ - 6.0)/(*epsilonIJ * *gammaIJ), (1.0/6.0));
            }
            break;   
        case eCOMB_YANG: // Yang, JPhysChemA, Volume: 122, Page: 1672, Year: 2018. Combination rules for Morse.
            *sigmaIJ = ((sigmaI * sigmaJ) *  (sigmaI + sigmaJ))/ (pow(sigmaI,2) + pow(sigmaJ,2));
            *epsilonIJ = (2.0 * epsilonI * epsilonJ)/(epsilonI + epsilonJ);
            *gammaIJ = ((gammaI * gammaJ) *  (gammaI + gammaJ))/ (pow(gammaI,2) + pow(gammaJ,2));   
            break;
        case eCOMB_QI: // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. Combination rules for Buf-14-7. Cubic-mean for sigma, and Waldman-Hagler for epsilon. 
            *sigmaIJ = (pow(sigmaI,3) + pow(sigmaJ,3))/(pow(sigmaI,2) + pow(sigmaJ,2));
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ) * ((2.0 * pow(epsilonI,3) * pow(epsilonJ,3))/(pow(epsilonI,6) + pow(epsilonJ,6)));
            *gammaIJ = 0.5 * (gammaI + gammaJ);
            break;
        case eCOMB_WALDMAN_HAGLER: // Waldman & Hagler, J. Comp. Chem., Year: 1993. 
            *sigmaIJ = pow(((pow(sigmaI,6.0)+pow(sigmaJ,6.0))/2.0),(1.0/6.0));
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ) * ((2.0 * pow(epsilonI,3.0) * pow(epsilonJ,3.0))/(pow(epsilonI,6.0) + pow(epsilonJ,6.0)));
            *gammaIJ = 0.5 * (gammaI + gammaJ);;
            break;            
        case eCOMB_NONE:
        case eCOMB_NR:
            gmx_fatal(FARGS, "Unsupported combination rule %d for Buckingham", CombinationRule);
    }
}

} // namespace alexandria
