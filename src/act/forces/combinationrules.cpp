/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */
#include "combinationrules.h"

#include <cmath>

#include "act/basics/mutability.h"
#include "act/forcefield/forcefield_parametername.h"
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
    case eCOMB_NONE:
        break;
    default:
        gmx_fatal(FARGS, "Unsupported combination rule %d for Lennard Jones", CombinationRule);
    }
}

void CombineLJ_86(int     CombinationRule,
               double  sigmaI,
               double  sigmaJ,
               double  epsilonI,
               double  epsilonJ,
               double *c6,
               double *c8)
{
    switch (CombinationRule)
    {
    case eCOMB_GEOMETRIC:
        {
            double sig  = std::sqrt(sigmaI * sigmaJ);
            double eps  = std::sqrt(epsilonI * epsilonJ);
	    double sig6 = std::pow(sig, 6.0);
            double sig2 = std::pow(sig, 2.0);
            *c6  = 4*eps*sig6;
            *c8 = *c6 * sig2;
        }
        break;
    case eCOMB_ARITHMETIC:
        {
            double sig  = 0.5 * (sigmaI + sigmaJ);
            double eps  = 0.5 * (epsilonI + epsilonJ);
            double sig2 = std::pow(sig, 2.0);
	    double sig6 = std::pow(sig, 6.0);
            *c6  = 4*eps*sig6;
            *c8 = *c6 * sig2;
        }
        break;
    case eCOMB_LORENTZ_BERTHELOT:
        {
            double sig  = 0.5 * (sigmaI + sigmaJ);
            double eps  = std::sqrt(epsilonI * epsilonJ);
            double sig2 = std::pow(sig, 2.0);
	    double sig6 = std::pow(sig, 6.0);
            *c6  = 4*eps*sig6;
            *c8 = *c6 * sig2;
        }
        break;
    case eCOMB_NONE:
        break;
    default:
        gmx_fatal(FARGS, "Unsupported combination rule %d for Lennard Jones", CombinationRule);
    }
}

void CombineLJ_147(int     CombinationRule,
                 double  sigmaI,
                 double  sigmaJ,
                 double  epsilonI,
                 double  epsilonJ,
                 double  gammaI,
                 double  gammaJ,
		 double  deltaI,
		 double  deltaJ,
                 double *sigmaIJ,
                 double *epsilonIJ,
                 double *gammaIJ,
                 double *deltaIJ)

{
    switch (CombinationRule)
    {
        case eCOMB_GEOMETRIC:
            *sigmaIJ = std::sqrt(sigmaI * sigmaJ);
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ);
            *gammaIJ = std::sqrt(gammaI * gammaJ);
            *deltaIJ = std::sqrt(deltaI * deltaJ);
            break;
        case eCOMB_ARITHMETIC:
            *sigmaIJ = 0.5 * (sigmaI + sigmaJ);
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ);
            *gammaIJ = 0.5 * (gammaI + gammaJ);
            *deltaIJ = 0.5 * (deltaI + deltaJ);
            break;
        case eCOMB_QI: // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. Combination rules for Buf-14-7. Cubic-mean for sigma, and Waldman-Hagler for epsilon. 
            *sigmaIJ = (pow(sigmaI,3) + pow(sigmaJ,3))/(pow(sigmaI,2) + pow(sigmaJ,2));
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ) * ((2.0 * pow(sigmaI,3) * pow(sigmaJ,3))/(pow(sigmaI,6) + pow(sigmaJ,6)));
            *gammaIJ = 0.5 * (gammaI + gammaJ);
            *deltaIJ = std::sqrt(deltaI * deltaJ);
            break;
        case eCOMB_QI_2: // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. Combination rules for Buf-14-7. Cubic-mean for sigma, and Waldman-Hagler for epsilon. 2023 testing, Kriz. is almost the same asi Qi. but has WH rule for gamma, as is the winner for WBH
            *sigmaIJ = (pow(sigmaI,3) + pow(sigmaJ,3))/(pow(sigmaI,2) + pow(sigmaJ,2));
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ) * ((2.0 * pow(sigmaI,3) * pow(sigmaJ,3))/(pow(sigmaI,6) + pow(sigmaJ,6)));
            *gammaIJ = pow(((pow(gammaI,6.0)+pow(gammaJ,6.0))/2.0),(1.0/6.0));
            *deltaIJ = std::sqrt(deltaI * deltaJ);
            break;
        case eCOMB_QKmQG: // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. The best combination rules for Buf-14-7. Cubic-mean for sigma, and Waldman-Hagler for epsilon. Qi /WH for epsilon, KM for gamma (but with geometric sigmaIJ), qi for sigma and geometric for delta
            *sigmaIJ = (pow(sigmaI,3) + pow(sigmaJ,3))/(pow(sigmaI,2) + pow(sigmaJ,2));
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ) * ((2.0 * pow(sigmaI,3) * pow(sigmaJ,3))/(pow(sigmaI,6) + pow(sigmaJ,6)));
            *gammaIJ = std::sqrt(sigmaI * sigmaJ) * (0.5*((gammaI/sigmaI)+(gammaJ/sigmaJ)));
            *deltaIJ = std::sqrt(deltaI * deltaJ);
            break;	    
    case eCOMB_NONE:
        break;
    default:
        gmx_fatal(FARGS, "Unsupported combination rule %d for 14_7 Lennard Jones. Use Geometric, Arithmetic, Qi or Qi_2 (special for 14-7 with Kong-Mason gamma)", CombinationRule);
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
        case eCOMB_KONG_MASON: // Kong, C. L. Combining Rules for Intermolecular Potential Parameters. II. Rules for the Lennard-Jones (12−6) Potential and the Morse Potential. J. Chem. Phys. 1973, 59.
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
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ) * ((2.0 * pow(sigmaI,3) * pow(sigmaJ,3))/(pow(sigmaI,6) + pow(sigmaJ,6)));
            *gammaIJ = 0.5 * (gammaI + gammaJ);
            break;
        case eCOMB_QI_2: // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. Combination rules for Buf-14-7. Cubic-mean for sigma, and Waldman-Hagler for epsilon. 2023 testing, Kriz. is almost the same asi Qi
            *sigmaIJ = (pow(sigmaI,3) + pow(sigmaJ,3))/(pow(sigmaI,2) + pow(sigmaJ,2));
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ) * ((2.0 * pow(sigmaI,3) * pow(sigmaJ,3))/(pow(sigmaI,6) + pow(sigmaJ,6)));
            *gammaIJ = pow(((pow(gammaI,6.0)+pow(gammaJ,6.0))/2.0),(1.0/6.0));
            break;    
        case eCOMB_WALDMAN_HAGLER: // Waldman & Hagler, J. Comp. Chem., Year: 1993. 
            *sigmaIJ = pow(((pow(sigmaI,6.0)+pow(sigmaJ,6.0))/2.0),(1.0/6.0));
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ) * ((2.0 * pow(sigmaI,3.0) * pow(sigmaJ,3.0))/(pow(sigmaI,6.0) + pow(sigmaJ,6.0)));
            *gammaIJ = 0.5 * (gammaI + gammaJ);;
            break;
	case eCOMB_QYQY: // Waldman & Hagler, J. Comp. Chem., Year: 1993. kriz changing to the best rule for GBHAM Qi, Yang, Qi (with the "yang" for delta)
            *sigmaIJ = (pow(sigmaI,3) + pow(sigmaJ,3))/(pow(sigmaI,2) + pow(sigmaJ,2));
            *epsilonIJ = std::sqrt(epsilonI * epsilonJ) * ((2.0 * pow(sigmaI,3.0) * pow(sigmaJ,3.0))/(pow(sigmaI,6.0) + pow(sigmaJ,6.0)));
            *gammaIJ = ((gammaI * gammaJ) *  (gammaI + gammaJ))/ (pow(gammaI,2) + pow(gammaJ,2));
            break;
        case eCOMB_NONE:
            break;
        case eCOMB_NR:
            gmx_fatal(FARGS, "Unsupported combination rule %d for Buckingham", CombinationRule);
    }
}

void CombineGBham(int     CombinationRule,
                  double  rminI,
                  double  rminJ,
                  double  epsilonI,
                  double  epsilonJ,
                  double  gammaI,
                  double  gammaJ,
                  double  deltaI,
                  double  deltaJ,
                  double *rminIJ,
                  double *epsilonIJ,
                  double *gammaIJ,
                  double *deltaIJ)
{
    // This is just a quick hack!
    CombineBham(CombinationRule, rminI, rminJ, epsilonI, epsilonJ,
                gammaI, gammaJ, rminIJ, epsilonIJ, gammaIJ);
//    *deltaIJ = (deltaI+deltaJ)/2;
      *deltaIJ = ((deltaI * deltaJ) *  (deltaI + deltaJ))/ (pow(deltaI,2) + pow(deltaJ,2));
}

int getCombinationRule(const ForceFieldParameterList &vdw)
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

static void generateVdwParameterPairs(ForceField *pd)
{
    auto forcesVdw = pd->findForces(InteractionType::VDW);
    auto ftypeVdW  = forcesVdw->gromacsType();
    int  comb_rule = getCombinationRule(*forcesVdw);

    // We temporarily store the new parameters here
    ForceFieldParameterListMap *parm = forcesVdw->parameters();;
    
    // Fudge unit
    std::string unit("kJ/mol");
    
    // We use dependent mutability to show these are not independent params
    auto mutd = Mutability::Dependent;

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
            // File the parameters, potential dependent
            ForceFieldParameterMap pmap;
            switch (ftypeVdW)
            {
            case F_LJ:
                {
                    auto csigma     = lj_name[ljSIGMA];
                    auto cepsilon   = lj_name[ljEPSILON];
                    double isigma   = ivdw.second[csigma].internalValue();
                    double iepsilon = ivdw.second[cepsilon].internalValue();
                    double jsigma   = jvdw.second[csigma].internalValue();
                    double jepsilon = jvdw.second[cepsilon].internalValue();
                    double c6 = 0, c12 = 0;
                    CombineLJ(comb_rule, isigma, jsigma,
                              iepsilon, jepsilon, &c6, &c12);
                    pmap.insert({lj_name[ljC6_IJ], ForceFieldParameter(unit, c6, 0, 1, c6, c6, 
                                                                       mutd, true, true)});
                    pmap.insert({lj_name[ljC12_IJ], ForceFieldParameter(unit, c12, 0, 1, c12, c12, 
                                                                        mutd, true, true)});
                    // Add some dummy parameters
                    pmap.insert({csigma, ForceFieldParameter(unit, 0.3, 0, 1, 0.3, 0.3, 
                                                            mutd, true, true)});
                    pmap.insert({cepsilon, ForceFieldParameter(unit, 0, 0, 1, 0, 0, 
                                                              mutd, true, true)});
                }
                break;	
            case F_LJ_86:
                {
                    auto csigma     = lj_name[lj_86SIGMA];
                    auto cepsilon   = lj_name[lj_86EPSILON];
                    double isigma   = ivdw.second[csigma].internalValue();
                    double iepsilon = ivdw.second[cepsilon].internalValue();
                    double jsigma   = jvdw.second[csigma].internalValue();
                    double jepsilon = jvdw.second[cepsilon].internalValue();
                    double c6 = 0, c8 = 0;
                    CombineLJ_86(comb_rule, isigma, jsigma,
                              iepsilon, jepsilon, &c6, &c8);
                    pmap.insert({lj_name[lj_C6_IJ], ForceFieldParameter(unit, c6, 0, 1, c6, c6,
                                                                        mutd, true, true)});
                    pmap.insert({lj_name[lj_C8_IJ], ForceFieldParameter(unit, c8, 0, 1, c8, c8,
                                                                        mutd, true, true)});
                    // Add some dummy parameters
                    pmap.insert({csigma, ForceFieldParameter(unit, 0.3, 0, 1, 0.3, 0.3, 
                                                            mutd, true, true)});
                    pmap.insert({cepsilon, ForceFieldParameter(unit, 0, 0, 1, 0, 0, 
                                                              mutd, true, true)});
                }
                break;	
            case F_BHAM:
                {
                    auto csigma     = wbh_name[wbhSIGMA];
                    auto cepsilon   = wbh_name[wbhEPSILON];
                    auto cgamma     = wbh_name[wbhGAMMA];
                    double isigma   = ivdw.second[csigma].internalValue();
                    double iepsilon = ivdw.second[cepsilon].internalValue();
                    double igamma   = ivdw.second[cgamma].internalValue();
                    double jsigma   = jvdw.second[csigma].internalValue();
                    double jepsilon = jvdw.second[cepsilon].internalValue();
                    double jgamma   = jvdw.second[cgamma].internalValue();
                    double sigmaij = 0, epsilonij = 0, gammaij = 0;
                    CombineBham(comb_rule, isigma, jsigma,
                                iepsilon, jepsilon, 
                                igamma, jgamma, &sigmaij,
                                &epsilonij, &gammaij);
                    pmap.insert({wbh_name[wbhSIGMA_IJ], ForceFieldParameter(unit, sigmaij, 0, 1,
                                                                           sigmaij, sigmaij,
                                                                            mutd, true, true)});
                    pmap.insert({wbh_name[wbhEPSILON_IJ], ForceFieldParameter(unit, epsilonij, 0, 1,
                                                                              epsilonij, epsilonij,
                                                                              mutd, true, true)});
                    pmap.insert({wbh_name[wbhGAMMA_IJ], ForceFieldParameter(unit, gammaij, 0, 1,
                                                                            gammaij, gammaij,
                                                                            mutd, true, true)});
                    // Add some dummy parameters
                    pmap.insert({csigma, ForceFieldParameter(unit, 0.3, 0, 1, 0.3, 0.3, 
                                                             mutd, true, true)});
                    pmap.insert({cepsilon, ForceFieldParameter(unit, 0, 0, 1, 0, 0, 
                                                               mutd, true, true)});
                    pmap.insert({cgamma, ForceFieldParameter(unit, 10, 0, 1, 10, 10, 
                                                             mutd, true, true)});
                }
                break;
            case F_LJ_147:
                {
                    auto csigma     = lj_147_name[lj_147SIGMA];
                    auto cepsilon   = lj_147_name[lj_147EPSILON];
                    auto cgamma     = lj_147_name[lj_147GAMMA];
                    auto cdelta     = lj_147_name[lj_147DELTA];
                    double isigma   = ivdw.second[csigma].internalValue();
                    double iepsilon = ivdw.second[cepsilon].internalValue();
                    double igamma   = ivdw.second[cgamma].internalValue();
                    double idelta   = ivdw.second[cdelta].internalValue();
                    double jsigma   = jvdw.second[csigma].internalValue();
                    double jepsilon = jvdw.second[cepsilon].internalValue();
                    double jgamma   = jvdw.second[cgamma].internalValue();
<<<<<<< HEAD
                    double jdelta   = ivdw.second[cdelta].internalValue();
=======
		    double jdelta   = jvdw.second[cdelta].internalValue();
>>>>>>> 011a0e919 (147  commit openmm)
                    double sigmaij = 0, epsilonij = 0, gammaij = 0, deltaij = 0;
                    CombineLJ_147(comb_rule, isigma, jsigma,
                                iepsilon, jepsilon,
                                igamma, jgamma, idelta, jdelta, &sigmaij,
                                &epsilonij, &gammaij, &deltaij);
                    pmap.insert({lj_147_name[lj_147SIGMA_IJ], ForceFieldParameter(unit, sigmaij, 0, 1,
                                                                                  sigmaij, sigmaij,
                                                                                  mutd, true, true)});
                    pmap.insert({lj_147_name[lj_147EPSILON_IJ], ForceFieldParameter(unit, epsilonij, 0, 1,
                                                                                    epsilonij, epsilonij,
                                                                                    mutd, true, true)});
                    pmap.insert({lj_147_name[lj_147GAMMA_IJ], ForceFieldParameter(unit, gammaij, 0, 1,
                                                                                  gammaij, gammaij,
                                                                                  mutd, true, true)});
                    pmap.insert({lj_147_name[lj_147DELTA_IJ], ForceFieldParameter(unit, deltaij, 0, 1,
                                                                                  deltaij, deltaij,
                                                                                  mutd, true, true)});

                    // Add some dummy parameters
                    pmap.insert({csigma, ForceFieldParameter(unit, 0.3, 0, 1, 0.3, 0.3, 
                                                            mutd, true, true)});
                    pmap.insert({cepsilon, ForceFieldParameter(unit, 0, 0, 1, 0, 0, 
                                                               mutd, true, true)});
                    pmap.insert({cgamma, ForceFieldParameter(unit, 10, 0, 1, 10, 10, 
                                                             mutd, true, true)});
                    pmap.insert({cdelta, ForceFieldParameter(unit, 6, 0, 1, 6, 6, 
                                                             mutd, true, true)});
                }
                break;
            case F_GBHAM:
                {
                    auto crmin      = gbh_name[gbhRMIN];
                    auto cepsilon   = gbh_name[gbhEPSILON];
                    auto cgamma     = gbh_name[gbhGAMMA];
                    auto cdelta     = gbh_name[gbhDELTA];
                    double irmin    = ivdw.second[crmin].internalValue();
                    double iepsilon = ivdw.second[cepsilon].internalValue();
                    double igamma   = ivdw.second[cgamma].internalValue();
                    double idelta   = ivdw.second[cdelta].internalValue();
                    double jrmin    = jvdw.second[crmin].internalValue();
                    double jepsilon = jvdw.second[cepsilon].internalValue();
                    double jgamma   = jvdw.second[cgamma].internalValue();
                    double jdelta   = jvdw.second[cdelta].internalValue();
                    double rminij = 0, epsilonij = 0, gammaij = 0, deltaij = 0;
                    CombineGBham(comb_rule, irmin, jrmin, iepsilon, jepsilon, 
                                 igamma, jgamma, idelta, jdelta, &rminij,
                                 &epsilonij, &gammaij, &deltaij);
                    pmap.insert({gbh_name[gbhRMIN_IJ], ForceFieldParameter(unit, rminij, 0, 1,
                                                                          rminij, rminij,
                                                                           mutd, true, true)});
                    pmap.insert({gbh_name[gbhEPSILON_IJ],ForceFieldParameter(unit, epsilonij, 0, 1,
                                                                             epsilonij, epsilonij,
                                                                             mutd, true, true)});
                    pmap.insert({gbh_name[gbhGAMMA_IJ], ForceFieldParameter(unit, gammaij, 0, 1,
                                                                            gammaij, gammaij,
                                                                            mutd, true, true)});
                    pmap.insert({gbh_name[gbhDELTA_IJ], ForceFieldParameter(unit, deltaij, 0, 1,
                                                                            deltaij, deltaij,
                                                                            mutd, true, true)});
                    // Add some dummy parameters
                    pmap.insert({crmin, ForceFieldParameter(unit, 0.3, 0, 1, 0.3, 0.3, 
                                                            mutd, true, true)});
                    pmap.insert({cepsilon, ForceFieldParameter(unit, 0, 0, 1, 0, 0, 
                                                               mutd, true, true)});
                    pmap.insert({cgamma, ForceFieldParameter(unit, 10, 0, 1, 10, 10, 
                                                             mutd, true, true)});
                    pmap.insert({cdelta, ForceFieldParameter(unit, 6, 0, 1, 6, 6, 
                                                             mutd, true, true)});
                }
                break;
            default:
                fprintf(stderr, "Invalid van der waals type %s\n",
                        interaction_function[ftypeVdW].longname);
            }
            parm->insert_or_assign(Identifier({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes),
                                   std::move(pmap));
            
        }
    }
    // Finally add the new parameters to the existing list
    //auto fold = forcesVdw->parameters();
    //for(const auto &np : newParams.parametersConst())
    // {
    //   // Remove old copy if it exists
    //   auto oldfp = fold->find(np.first);
    //   if (oldfp != fold->end())
    //   {
    //       fold->erase(oldfp);
    //   }
    //   // Now add the new one
    //   fold->insert({ np.first, np.second });
    //}
    // Phew, we're done!
}

static void generateCoulombParameterPairs(ForceField *pd)
{
    auto forcesCoul = pd->findForces(InteractionType::COULOMB);
    
    // We temporarily store the new parameters here
    ForceFieldParameterList newParams;
    
    // Fudge unit
    std::string unit("kJ/mol");
    
    // We use dependent mutability to show these are not independent params
    auto mutd = Mutability::Dependent;

    // Now do the double loop
    for (auto &icoul : *forcesCoul->parameters())
    {
        auto iid    = icoul.first;
        // Check whether this is a single atom parameter
        if (iid.atoms().size() > 1)
        {
            continue;
        }
        auto iparam = icoul.second;
        double izeta = icoul.second["zeta"].internalValue();
        for (auto &jcoul : *forcesCoul->parameters())
        {
            auto jid    = jcoul.first;
            // Check whether this is a single atom parameter and
            // whether this is is larger or equal to iid.
            if (jid.atoms().size() > 1 || jid.id() < iid.id())
            {
                continue;
            }
            auto jparam = jcoul.second;
            double jzeta  = jcoul.second["zeta"].internalValue();
            Identifier pairID({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes);
            ForceFieldParameter pi(unit, izeta, 0, 1, izeta, izeta, mutd, true, true);
            ForceFieldParameter pj(unit, jzeta, 0, 1, jzeta, jzeta, mutd, true, true);
            newParams.addParameter(pairID, coul_name[coulZETAI], pi);
            newParams.addParameter(pairID, coul_name[coulZETAJ], pj);
        }
    }
    // Finally add the new parameters to the exisiting list
    auto fold = forcesCoul->parameters();
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

static void generateShellForceConstants(ForceField *pd)
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
        if (part.second.hasOption("poltype"))
        {
            auto shellType = part.second.optionValue("poltype");
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

void generateDependentParameter(ForceField *pd)
{
    generateVdwParameterPairs(pd);
    generateCoulombParameterPairs(pd);
    generateShellForceConstants(pd);
}

} // namespace alexandria
