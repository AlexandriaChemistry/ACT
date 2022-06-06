#include "combinationrules.h"

#include <cmath>

#include "act/basics/mutability.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
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

static int getCombinationRule(const ForceFieldParameterList &vdw)
{
    auto combRule = vdw.optionValue("combination_rule");
    int  i;
    for(i = 0; i < eCOMB_NR; i++)
    {
        if (combRule.compare(ecomb_names[i]) == 0)
        {
            break;
        }
    }
    GMX_RELEASE_ASSERT(i < eCOMB_NR, gmx::formatString("Cannot find combination rule %s in GROMACS",
                                                       combRule.c_str()).c_str());
    return i;
}

static void generateVdwParameterPairs(Poldata *pd)
{
    auto forcesVdw = pd->findForces(InteractionType::VDW);
    auto ftypeVdW  = forcesVdw->fType();
    int  comb_rule = getCombinationRule(*forcesVdw);
    
    // We temporarily store the new parameters here
    ForceFieldParameterList newParams;
    
    // Fudge unit
    std::string unit("kJ/mol");
    
    // We use dependent mutability to show these are not independent params
    auto mut = Mutability::Dependent;

    // Now do the double loop
    for (auto &ivdw : *forcesVdw->parameters())
    {
        auto iid    = ivdw.first;
        // Check whether this is a single atom parameter
        if (iid.atoms().size() > 1)
        {
            continue;
        }
        auto iparam = ivdw.second;
        double igamma   = 0;
        double isigma   = ivdw.second["sigma"].internalValue();
        double iepsilon = ivdw.second["epsilon"].internalValue();
        if (ftypeVdW == F_BHAM)
        {
            igamma = ivdw.second["gamma"].internalValue();
        }
        for (auto &jvdw : *forcesVdw->parameters())
        {
            auto jid    = jvdw.first;
            // Check whether this is a single atom parameter and
            // whether this is is larger or equal to iid.
            if (jid.atoms().size() > 1 || jid.id() < iid.id())
            {
                continue;
            }
            auto jparam = jvdw.second;
            double jgamma   = 0;
            double jsigma   = jvdw.second["sigma"].internalValue();
            double jepsilon = jvdw.second["epsilon"].internalValue();
            if (ftypeVdW == F_BHAM)
            {
                jgamma = jvdw.second["gamma"].internalValue();
            }
            Identifier pairID({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes);
            switch (ftypeVdW)
            {
            case F_LJ:
                {
                    double c6 = 0, c12 = 0;
                    CombineLJ(comb_rule, isigma, jsigma,
                              iepsilon, jepsilon, &c6, &c12);
                    ForceFieldParameter c6parm(unit, c6, 0, 1, c6, c6, 
                                               mut, true, true);
                    ForceFieldParameter c12parm(unit, c12, 0, 1, c12, c12, 
                                                mut, true, true);
                    newParams.addParameter(pairID, "c6_ij", c6parm);
                    newParams.addParameter(pairID, "c12_ij", c12parm);
                }
                break;
            case F_BHAM:
                {
                    double sigmaij = 0, epsilonij = 0, gammaij = 0;
                    CombineBham(comb_rule, isigma, jsigma,
                                iepsilon, jepsilon, 
                                igamma, jgamma, &sigmaij,
                                &epsilonij, &gammaij);
                    ForceFieldParameter sigparm(unit, sigmaij, 0, 1,
                                                sigmaij, sigmaij,
                                                mut, true, true);
                    ForceFieldParameter epsparm(unit, epsilonij, 0, 1,
                                                epsilonij, epsilonij,
                                                mut, true, true);
                    ForceFieldParameter gamparm(unit, gammaij, 0, 1,
                                                gammaij, gammaij,
                                                mut, true, true);
                    newParams.addParameter(pairID, "sigma_ij", sigparm);
                    newParams.addParameter(pairID, "epsilon_ij", epsparm);
                    newParams.addParameter(pairID, "gamma_ij", gamparm);
                }
                break;
            default:
                fprintf(stderr, "Invalid van der waals type %s\n",
                        interaction_function[ftypeVdW].longname);
            }
        }
    }
    // Finally add the new parameters to the exisiting list
    auto fold = forcesVdw->parameters();
    for(const auto &np : newParams.parametersConst())
    {
        // Remove old copy if it exists
        auto oldfp = fold->find(np.first);
        if (oldfp != fold->end())
        {
            fold->erase(oldfp);
        }
        // Now add the new one
        fold->insert({ np.first, np.second });
    }
    // Phew, we're done!
}

static void generateShellForceConstants(Poldata *pd)
{
    if (!pd->polarizable())
    {
        return;
    }
    auto itype = InteractionType::POLARIZATION;
    auto ffpl  = pd->findForces(itype)->parameters();
    // Loop over particles
    for(const auto &part : pd->particleTypesConst())
    {
        if (part.hasOption("poltype"))
        {
            auto shellType = part.optionValue("poltype");
            if (ffpl->find(shellType) == ffpl->end())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Missing polarization term for %s", shellType.c_str()).c_str()));
            }
            auto  &parms   = ffpl->find(shellType)->second;
            auto   alpha   = parms.find("alpha")->second.internalValue();
            auto   qshell  = pd->findParticleType(shellType)->charge();
            double kshell  = 0;
            if (alpha > 0 && qshell != 0)
            {
                kshell = gmx::square(qshell)*ONE_4PI_EPS0/alpha;
            }
            std::string fc_name("kshell");
            auto fc_parm = parms.find(fc_name);
            if (parms.end() == fc_parm)
            {
                ForceFieldParameter fc_new("kJ/mol nm2", kshell, 0, 1,
                                           kshell, kshell,
                                           Mutability::Dependent, true, true);
                parms.insert({ fc_name, fc_new });
            }
            else
            {
                fc_parm->second.setMutability(Mutability::Free);
                fc_parm->second.setMinimum(kshell);
                fc_parm->second.setMaximum(kshell);
                fc_parm->second.setValue(kshell);
                fc_parm->second.setMutability(Mutability::Dependent);
            }
        }
    }
}

void generateDependentParameter(Poldata *pd)
{
    generateVdwParameterPairs(pd);
    generateShellForceConstants(pd);
}

} // namespace alexandria
