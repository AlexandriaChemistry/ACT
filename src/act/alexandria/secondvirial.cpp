/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022,2023
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

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "secondvirial.h"

#include <cctype>
#include <cstdlib>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/atype_mapping.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/fetch_charges.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/actmol.h"
#include "act/alexandria/princ.h"
#include "act/alexandria/train_utility.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/jsontree.h"
#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

const std::map<b2Type, std::string> b2Type2str = {
    { b2Type::Classical, "B2(Classical)" },
    { b2Type::Force, "B2(Force)" },
    { b2Type::Torque, "B2(Torque)" },
    { b2Type::Total, "B2(Total)" }
};                    

const std::string &b2TypeToString(b2Type b2t)
{
    return b2Type2str.find(b2t)->second;
}

void forceFieldSummary(JsonTree      *jtree,
                       const ForceField *pd)
{
    jtree->addObject(JsonTree("Force field file", pd->filename()));
    jtree->addObject(JsonTree("Created", pd->timeStamp()));
    jtree->addObject(JsonTree("Checksum", pd->checkSum()));
    jtree->addObject(JsonTree("Polarizable", yesno_names[pd->polarizable()]));
    jtree->addObject(JsonTree("Charge generation", 
                              chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str()));
    int nexclvdw;
    if (!ffOption(*pd, InteractionType::VDW, 
                  "nexcl", &nexclvdw))
    {
        nexclvdw = 0;
    }
    jtree->addObject(JsonTree("# vanderwaals exclusions", gmx_itoa(nexclvdw)));
    int nexclqq;
    if (!ffOption(*pd, InteractionType::COULOMB, 
                  "nexcl", &nexclqq))
    {
        nexclqq = 0;
    }
    jtree->addObject(JsonTree("# coulomb exclusions", gmx_itoa(nexclqq)));
    double epsilonr;
    if (!ffOption(*pd, InteractionType::COULOMB, 
                  "epsilonr", &epsilonr))
    {
        epsilonr = 1;
    }
    jtree->addObject(JsonTree("Relative dielectric constant epsilon_r",
                              gmx_ftoa(epsilonr)));
    jtree->addObject(JsonTree("# particle types", gmx_itoa(pd->getNatypes())));
    
    JsonTree ftree("InteractionTypes");
    for(const auto &fs : pd->forcesConst())
    {
        auto itype = fs.first;
        auto &ffpl = fs.second;
        if (!ffpl.parametersConst().empty())
        {
            JsonTree fftree(interactionTypeToString(itype));
            if (!ffpl.function().empty())
            {
                fftree.addObject(JsonTree("Force function", ffpl.function()));
            }
            fftree.addObject(JsonTree("# entries",
                                      gmx_itoa(ffpl.parametersConst().size())));
            ftree.addObject(fftree);
        }
    }
    if (!ftree.empty())
    {
        jtree->addObject(ftree);
    }
}

double sphereIntegrator(double r1, double r2, double val1, double val2)
{
    // Approximate trapezium by y = ax + b (a == slope)
    // then integrate that multiplied by x^2 to get
    // a/4 x^4 + b/3 x^3
    // insert old and new point
    double a        = (val2-val1)/(r2-r1);
    double b        = val1 - a*r1;
    double integral = ((a/4)*(r2*r2*r2*r2 - r1*r1*r1*r1) + 
                       (b/3)*(r2*r2*r2 - r1*r1*r1));
    return 4*M_PI*integral;
}

void ReRunner::addOptions(std::vector<t_pargs>  *pargs,
                          std::vector<t_filenm> *filenm)
{
    std::vector<t_pargs> pa = {
        { "-traj",   FALSE, etSTR,  {&trajname_},
          "Trajectory or series of structures of the same compound for which the energies will be computed. If this option is present, no simulation will be performed." },
        { "-T1",     FALSE, etREAL, {&T1_},
          "Starting temperature for second virial calculations." },
        { "-T2",     FALSE, etREAL, {&T2_},
          "Starting temperature for second virial calculations." },
        { "-dT",     FALSE, etREAL, {&deltaT_},
          "Temperature increment for calculation of second virial." },
        { "-nbootstrap", FALSE, etINT, {&nbootStrap_},
          "Number of times to iterate B2 calculation based on the same dimer orientations and energies" },
        { "-optimizedB2", FALSE, etBOOL, {&optimizedB2_},
          "Optimize the bootstrapping by pre-calculating stuff" }
    };
    pargs->push_back(pa[0]);
    if (computeB2_)
    {
        for(size_t i = 1; i < pa.size(); i++)
        {
            pargs->push_back(pa[i]);
        }
    
        std::vector<t_filenm> fnm = {
            { efXVG, "-eh", "mayer",      ffOPTWR },
            { efXVG, "-b2", "B2T",        ffOPTWR }
        };
        for(const auto &fn : fnm)
        {
            filenm->push_back(fn);
        }
    }
}

const std::vector<double> &ReRunner::temperatures()
{
    if (Temperatures_.empty())
    {
        Temperatures_.push_back(T1_);
        if (deltaT_ > 0)
        {
            double T = T1_+deltaT_;
            while (T <= T2_)
            {
                Temperatures_.push_back(T);
                T += deltaT_;
            }
        }
        else
        {
            Temperatures_.push_back(T2_);
        }
    }
    return Temperatures_;
}

void ReRunner::plotMayer(const char                             *ehisto,
                         const std::vector<std::vector<double>> &mayer)
{
    if (!ehisto || strlen(ehisto) == 0)
    {
        return;
    }

    FILE *fp = xvgropen(ehisto, "Mayer function", "r (nm)",
                        "< exp[-U12/kBT]-1 >", oenv_);
    std::vector<std::string> label;
    for(auto T : temperatures())
    {
        if (T != 0)
        {
            label.push_back(gmx::formatString("T = %g K", T));
        }
    }
    xvgrLegend(fp, label, oenv_);
    for(size_t i = 0; i < mayer[0].size(); i++)
    {
        for(size_t j = 0; j < mayer.size(); j++)
        {
            fprintf(fp, "  %10g", mayer[j][i]);
        }
        fprintf(fp, "\n");
    }
    xvgrclose(fp);
}

void ReRunner::plotB2temp(const char *b2file)
{
    if (!b2file || strlen(b2file) == 0)
    {
        return;
    }
    auto T = temperatures();
    if (b2t_.empty() || b2t_[b2Type::Total].size() != T.size())
    {
        fprintf(stderr, "Internal inconsistency. There are %zu temperatures and %zu B2 values.\n",
                T.size(), b2t_.size());
        return;
    }
    FILE *b2p = xvgropen(b2file, "Second virial coefficient",
                         "Temperature (K)", "B2(T) cm^3/mol", oenv_);
    std::vector<std::string> legend = {
        "Total", "Classical", "Force", "Torque"
    };
    xvgrLegend(b2p, legend, oenv_);
    for(size_t ii = 0; ii < T.size(); ii++)
    {
        fprintf(b2p, "%10g  %10g  %10g  %10g  %10g\n", T[ii], 
                b2t_[b2Type::Total][ii], b2t_[b2Type::Classical][ii],
                b2t_[b2Type::Force][ii], b2t_[b2Type::Torque][ii]);
    }
    xvgrclose(b2p);
}

void ReRunner::computeB2(FILE                                      *logFile,
                         gmx_stats                                  edist,
                         int                                        ndist,
                         const std::vector<double>                 &mass,
                         const std::vector<gmx::RVec>              &inertia,
                         const std::vector<std::vector<gmx::RVec>> &forceMol,
                         const std::vector<std::vector<gmx::RVec>> &torqueMol,
                         const std::vector<t_filenm>               &fnm)
{
    auto                      N = edist.get_npoints();
    if (N <= 0)
    {
        return;
    }
    const std::vector<double> x = edist.getX();
    const std::vector<double> y = edist.getY();
    double xmin                 = *std::min_element(x.begin(), x.end());
    double xmax                 = *std::max_element(x.begin(), x.end());

    if (N > 2 && xmax > xmin)
    {
        // Default bin width
        double binWidth = 0.0025; // nm
        if (ndist > 1)
        {
            binWidth = 2*(xmax-xmin)/(ndist-1);
        }
        // Bins start from zero for proper integration
        size_t nbins    = 1+std::round(xmax/binWidth);
        // Temp array to store distance and Mayer functions
        auto Temperature = temperatures();
        std::vector<std::vector<double> > mayer(1+Temperature.size());
        // Array for total second virial as a function of T
        size_t iTemp = 1;
        // Will be used to obtain a seed for the random number engine
        std::random_device                 bsRand;  
        //Standard mersenne_twister_engine seeded with rd()
        std::mt19937                       bsGen(bsRand());
        std::uniform_int_distribution<int> bsDistr(0, x.size()-1);
        for(auto T : Temperature)
        {
            if (T == 0)
            {
                fprintf(stderr, "Please provide a finite temperature to compute second virial.\n");
                continue;
            }
            if (nbootStrap_ < 1)
            {
                nbootStrap_ = 1;
            }
            double                      beta = 1.0/(BOLTZ*T);
            std::map<b2Type, gmx_stats> b2BootStrap;
            for(auto &b2b : b2Type2str)
            {
                gmx_stats gs;
                b2BootStrap.insert({ b2b.first, std::move(gs) });
            }
            // Store data temporarily in a structure for memory performance
            typedef struct
            {
                int       index;
                double    g0;
                double    g0_f2[2];
                gmx::RVec g0_tau[2];
            } b2temp_t;
            std::vector<b2temp_t> b2temp(x.size());
            if (optimizedB2_)
            {
                for(size_t jj = 0; jj < x.size(); jj++)
                {
                    double rindex = x[jj]/binWidth;
                    b2temp[jj].index  = rindex;
                    b2temp[jj].g0     = std::exp(-y[jj]*beta);
                    for(int kk = 0; kk < 2; kk++)
                    {
                        b2temp[jj].g0_f2[kk] = b2temp[jj].g0*iprod(forceMol[kk][jj], forceMol[kk][jj]);
                        for(int m = 0; m < DIM; m++)
                        {
                            b2temp[jj].g0_tau[kk][m] += b2temp[jj].g0*torqueMol[kk][jj][m]*torqueMol[kk][jj][m];
                        }
                    }
                }
            }
            
            for(int nb = 0; nb < nbootStrap_; nb++)
            {
                // Temporary arrays for weighted properties.
                std::vector<double>    exp_U12(nbins, 0.0);
                std::vector<double>    exp_F2[2];
                std::vector<gmx::RVec> exp_tau[2];
                for(int kk = 0; kk < 2; kk++)
                {
                    exp_F2[kk].resize(nbins, 0);
                    exp_tau[kk].resize(nbins, { 0.0, 0.0, 0.0 });
                }
                std::vector<int>       n_U12(nbins, 0);
                for(size_t jj = 0; jj < x.size(); jj++)
                {
                    size_t ii = jj;
                    if (nbootStrap_ > 1)
                    {
                        ii = bsDistr(bsGen);
                    }
                    size_t index;
                    if (optimizedB2_)
                    {
                        index = b2temp[ii].index;
                        double g0_12 = b2temp[ii].g0;
                        exp_U12[index] += g0_12-1;
                        for(int kk = 0; kk < 2; kk++)
                        {
                            exp_F2[kk][index] += b2temp[ii].g0_f2[kk];
                            for(int m = 0; m < DIM; m++)
                            {
                                // Gray and Gubbins Eqn. 3.282
                                exp_tau[kk][index][m] += b2temp[ii].g0_tau[kk][m];
                            }
                        }
                    }
                    else
                    {
                        double rindex = x[ii]/binWidth;
                        index  = rindex;
                        // Gray and Gubbins Eqn. 3.261
                        double g0_12 = std::exp(-y[ii]*beta);
                        // Gray and Gubbins Eqn. 3.272
                        exp_U12[index] += g0_12-1;
                        for(int kk = 0; kk < 2; kk++)
                        {
                            // Gray and Gubbins Eqn. 3.281
                            exp_F2[kk][index] += g0_12*iprod(forceMol[kk][ii], forceMol[kk][ii]);
                            for(int m = 0; m < DIM; m++)
                            {
                                // Gray and Gubbins Eqn. 3.282
                                exp_tau[kk][index][m] += g0_12*torqueMol[kk][ii][m]*torqueMol[kk][ii][m];
                            }
                        }
                    }
                    n_U12[index]   += 1;
                }
                double    Bclass       =  0;
                double    BqmForce     =  0;
                double    BqmTorque[2] =  { 0, 0 };
                // We start in the origin even if there is no data.
                double    r1        =  0;
                // Starting energy, all values until first data entry
                double    Uprev     = -1;
                int jj = 0;
                while(jj*binWidth < xmin)
                {
                    exp_U12[jj] = -1;
                    n_U12[jj]   = 1;
                    jj += 1;
                }
                // Starting force
                double    Fprev[2] =  { 0, 0 };
                // Starting torque
                gmx::RVec Tprev[2] = { { 0, 0, 0 }, { 0, 0, 0 } };
                double    hbarfac  = beta*gmx::square(beta*PLANCK/(2*M_PI))/24;
                if (iTemp == 1)
                {
                    // Store distance first time around only
                    mayer[0].push_back(r1);
                }
                mayer[iTemp].push_back(Uprev);
                for(size_t ii = 1; ii < nbins; ii++)
                {
                    double r2 = ii*binWidth;
                    if (n_U12[ii] > 0)
                    {
                        double Unew = exp_U12[ii]/n_U12[ii];
                        if (iTemp == 1)
                        {
                            // Store distance first time around only
                            mayer[0].push_back(r2);
                        }
                        mayer[iTemp].push_back(Unew);
                        auto dB       = sphereIntegrator(r1, r2, Uprev, Unew);
                        // TODO: There is factor 0.5 here
                        Bclass       -= 0.5*dB;
                        Uprev         = Unew;
                        for(int kk = 0; kk < 2; kk++)
                        {
                            // Weighted square force
                            // We follow Eqn. 9 in Schenter, JCP 117 (2002) 6573
                            double Fnew  = exp_F2[kk][ii]/(mass[kk]*n_U12[ii]);
                            BqmForce    += 0.5*hbarfac*sphereIntegrator(r1, r2, Fprev[kk], Fnew);
                            Fprev[kk]    = Fnew;
                            // Contributions from torque
                            for(int m = 0; m < DIM; m++)
                            {
                                if (inertia[kk][m] > 0)
                                {
                                    double Tnew    = exp_tau[kk][ii][m]/(n_U12[ii]*inertia[kk][m]);
                                    BqmTorque[kk] += hbarfac*sphereIntegrator(r1, r2, Tprev[kk][m], Tnew);
                                    Tprev[kk][m]   = Tnew;
                                }
                            }
                        }
                    }
                    r1 = r2;
                }
                // Conversion to regular units cm^3/mol.
                double fac  = AVOGADRO*1e-21;
                double bqt  = (BqmTorque[0]+BqmTorque[1])*0.5;
                double Btot = (Bclass + BqmForce + bqt)*fac;
                // Add to bootstrapping statistics
                b2BootStrap[b2Type::Classical].add_point(nb, Bclass*fac, 0, 0);
                b2BootStrap[b2Type::Force].add_point(nb, BqmForce*fac, 0, 0);
                b2BootStrap[b2Type::Torque].add_point(nb, bqt*fac, 0, 0);
                b2BootStrap[b2Type::Total].add_point(nb, Btot, 0, 0);
            }
            // Done bootstrapping, store the result.
            for(const auto &b2b : b2Type2str)
            {
                real aver, sigma, error;
                b2BootStrap[b2b.first].get_ase(&aver, &sigma, &error);
                b2t_[b2b.first].push_back(aver);
                b2tError_[b2b.first].push_back(sigma);
            }
            if (logFile)
            {
                fprintf(logFile, "T = %g K. ", T);
                for(const auto &b2b : b2Type2str)
                {
                    fprintf(logFile, " %s %8.1f (%5.1f)", b2b.second.c_str(),
                            b2t_[b2b.first].back(), b2tError_[b2b.first].back());
                }
                fprintf(logFile, "\n");
            }
            iTemp += 1;
        }
        if (!fnm.empty())
        {
            // Store Mayer functions if requested.
            plotMayer(opt2fn_null("-eh", fnm.size(), fnm.data()), mayer);
            // Print B2(T) if requested
            plotB2temp(opt2fn("-b2", fnm.size(), fnm.data()));
        }
    }
}

void ReRunner::rerun(FILE                        *logFile,
                     const ForceField            *pd,
                     const ACTMol                *actmol,
                     bool                         userqtot,
                     double                       qtot,
                     bool                         verbose,
                     const std::vector<t_filenm> &fnm)
{
    std::vector<std::vector<gmx::RVec> > dimers;
    std::string          method, basis;
    int                  maxpot = 100;
    int                  nsymm  = 1;
    int                  ndist  = 0;
    const char          *molnm  = "";
    if (verbose && debug)
    {
        print_memory_usage(debug);
    }
    if (trajname_ && strlen(trajname_) > 0)
    {
        std::vector<MolProp> mps;
        std::string tname(trajname_);
        auto pos = tname.find(".xml");
        if (pos != std::string::npos && tname.size() == pos+4)
        {
            // Assume this is a molprop file
            MolPropRead(trajname_, &mps);
        }
        else
        {
            // Read compounds if we have a trajectory file
            matrix box;
            if (!readBabel(pd, trajname_, &mps, molnm, molnm, "", &method,
                           &basis, maxpot, nsymm, "Opt", userqtot, &qtot, false, box))
            {
                fprintf(stderr, "Could not read compounds from %s\n", trajname_);
                return;
            }
        }
        for(size_t i = 0; i < mps.size(); i++)
        {
            auto exper = mps[i].experimentConst();
            for(const auto &ep : exper)
            {
                std::vector<gmx::RVec> xx;
                for(const auto &epx: ep.getCoordinates())
                {
                    xx.push_back(epx);
                    if (pd->polarizable())
                    {
                        xx.push_back(epx);
                    }
                }
                dimers.push_back(xx);
            }
        }
        if (logFile)
        {
            fprintf(logFile, "Doing energy calculation for %zu structures from %s\n",
                    dimers.size(), trajname_);
            fflush(logFile);
        }
    }
    else
    {
        // Generate compounds
        gendimers_->generate(logFile, actmol, &dimers, opt2fn_null("-ox", fnm.size(), fnm.data()));
        ndist = gendimers_->ndist();
    }
    if (verbose && debug)
    {
        print_memory_usage(debug);
    }
    std::map<InteractionType, double> energies;
    int mp_index = 0;
    gmx_stats edist;
    std::vector<std::vector<gmx::RVec>> forceMol;
    std::vector<std::vector<gmx::RVec>> torqueMol;
    forceMol.resize(2);
    torqueMol.resize(2);
    if (eInter_)
    {
        for(int kk = 0; kk < 2; kk++)
        {
            torqueMol[kk].resize(dimers.size());
            forceMol[kk].resize(dimers.size());
        }
    }
    if (verbose && debug)
    {
        print_memory_usage(debug);
    }
    std::vector<gmx::RVec> inertia = { { 0, 0, 0 }, { 0, 0, 0 } };
    // Loop over molecules
    const auto &atoms = actmol->atomsConst();
    for (size_t idim = 0; idim < dimers.size(); idim++)
    {
        std::vector<gmx::RVec> coords;
        if (dimers[idim].size() == atoms.size())
        {
            // Assume there are shells in the input
            coords = dimers[idim];
        }
        else
        {
            size_t index = 0;
            for(size_t i = 0; i < atoms.size(); i++)
            {
                if (index < dimers[idim].size())
                {
                    coords.push_back(dimers[idim][index]);
                }
                else
                {
                    GMX_THROW(gmx::InvalidInputError(gmx::formatString("Number of (generated) coordinates in trajectory (%zu) does not match molecule file (%zu)", dimers[idim].size(), atoms.size()).c_str()));
                }
                if (atoms[i].pType() == eptAtom)
                {
                    index++;
                }
            }
        }
        std::vector<gmx::RVec> forces(coords.size());
        if (verbose)
        {
            fprintf(logFile, "%5d", mp_index);
        }
        if (eInter_)
        {
            std::map<InteractionType, double> einter;
            actmol->calculateInteractionEnergy(pd, forceComp_, &einter, &forces, &coords);
            auto atomStart  = actmol->fragmentHandler()->atomStart();
            std::vector<gmx::RVec> f    = { { 0, 0, 0 }, { 0, 0, 0 } };
            std::vector<gmx::RVec> com        = { { 0, 0, 0 }, { 0, 0, 0 } };
            std::vector<double>    mtot       = { 0, 0 };
            std::vector<gmx::RVec> torque     = { { 0, 0, 0 }, { 0, 0, 0 } };
            std::vector<gmx::RVec> torqueRot  = { { 0, 0, 0 }, { 0, 0, 0 } };
            auto      tops  = actmol->fragmentHandler()->topologies();
            for(int kk = 0; kk < 2; kk++)
            {
                auto natom = tops[kk]->atoms().size();
                for(size_t i = atomStart[kk]; i < atomStart[kk]+natom; i++)
                {
                    gmx::RVec mr1;
                    auto      mi = atoms[i].mass();
                    svmul(mi, coords[i], mr1);
                    mtot[kk] += mi;
                    // Compute center of mass of compound kk
                    rvec_inc(com[kk], mr1);
                    // Compute total force on compound kk
                    rvec_inc(f[kk], forces[i]);
                }
                GMX_RELEASE_ASSERT(mtot[kk] > 0, "Zero mass");
                for(size_t m = 0; m < DIM; m++)
                {
                    // Normalize
                    com[kk][m] /= mtot[kk];
                }
                copy_rvec(f[kk], forceMol[kk][mp_index]);

                // Compute the coordinates relative to the center of mass
                std::vector<real>      mass;
                std::vector<int>       index;
                std::vector<gmx::RVec> x_com;
                for(size_t i = atomStart[kk]; i < atomStart[kk]+natom; i++)
                {
                    // Store atom index relative to molecule start
                    index.push_back(i-atomStart[kk]);
                    // Store mass
                    mass.push_back(atoms[i].mass());
                    // Subtract COM and store
                    gmx::RVec ri;
                    rvec_sub(coords[i], com[kk], ri);
                    x_com.push_back(ri);
                }
                // Compute moments of inertia and transformation matrix
                gmx::RVec inertia1;
                clear_rvec(inertia1);
                matrix trans;
                principal_comp(index, mass, x_com, &trans, &inertia1);

                // Move to inertial frame (only well-defined for
                // rigid molecules).
                // The trans matrix should convert that coordinate to the inertial frame,
                // but what about the force on the atoms? It likely has to be rotated in the
                // same manner. After that, the torque can be computed.
                // Before doing the rotations, the sum of the torque vectors is zero.
                // Since the rotations are different for both molecules, this does not
                // hold after the rotations are done.
                // TODO: write out the math.
                for(size_t i = atomStart[kk]; i < atomStart[kk]+natom; i++)
                {
                    gmx::RVec ri, fi, ti;
                    cprod(x_com[i-atomStart[kk]], forces[i], ti);
                    rvec_inc(torque[kk], ti);
                    // Rotate coordinates
                    mvmul(trans, x_com[i-atomStart[kk]], ri);
                    // Rotate force vector
                    mvmul(trans, forces[i], fi);
                    // Compute torque on this atom
                    cprod(ri, fi, ti);
                    // Update total torque
                    rvec_inc(torqueRot[kk], ti);
                }
                rvec_inc(inertia[kk], inertia1);
                torqueMol[kk][mp_index] = torqueRot[kk];
            }
            if (verbose && debug)
            {
                print_memory_usage(debug);
            }
            gmx::RVec dcom;
            rvec_sub(com[0], com[1], dcom);
            double rcom = norm(dcom);
            fprintf(logFile, " r %g", rcom);
            for (auto &EE: einter)
            {
                fprintf(logFile, " %s %g", interactionTypeToString(EE.first).c_str(), EE.second);
            }
            if (verbose)
            {
                fprintf(logFile, " Force %g %g %g Torque[0] %g %g %g Torque[1] %g %g %g Rotated Torque[0] %g %g %g Rotated Torque[1] %g %g %g",
                        f[0][XX], f[0][YY], f[0][ZZ],
                        torque[0][XX], torque[0][YY], torque[0][ZZ],
                        torque[1][XX], torque[1][YY], torque[1][ZZ],
                        torqueRot[0][XX], torqueRot[0][YY], torqueRot[0][ZZ],
                        torqueRot[1][XX], torqueRot[1][YY], torqueRot[1][ZZ]);
            }
            else
            {
                fprintf(logFile, "\n");
            }
            edist.add_point(rcom, einter[InteractionType::EPOT], 0, 0);
        }
        else
        {
            forceComp_->compute(pd, actmol->topology(),
                                &coords, &forces, &energies);
            for(const auto &ee : energies)
            {
                fprintf(logFile, "  %s %8g", interactionTypeToString(ee.first).c_str(), ee.second);
            }
            fprintf(logFile, "\n");
        }
        if (verbose)
        {
            fprintf(logFile, "\n");
        }
        mp_index++;
    }
    if (verbose && debug)
    {
        print_memory_usage(debug);
    }
    if (computeB2_)
    {
        auto info = gmx::formatString("Done with energy calculations, now time for second virial.");
        printf("%s\n", info.c_str());
        if (logFile)
        {
            fprintf(logFile, "%s\n", info.c_str());
        }
        for(int kk = 0; kk < 2; kk++)
        {
            for(int m = 0; m < DIM; m++)
            {
                inertia[kk][m] /= edist.get_npoints();
            }
        }
        // Compute the relative mass
        std::vector<double> masses = {
            actmol->fragmentHandler()->topologies()[0]->mass(),
            actmol->fragmentHandler()->topologies()[1]->mass()
        };
        computeB2(logFile, edist, ndist, masses, inertia,
                  forceMol, torqueMol, fnm);
    }
    if (verbose && debug)
    {
        print_memory_usage(debug);
    }
}

int b2(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "alexandria b2 will generate conformations of dimers. The",
        "corresponding charges (molprop) file, used for generating the topology",
        "needs to contain information about",
        "the compounds in the dimer. Based on dimer energies the second",
        "virial coefficient will be estimated. Mayer curves can optionally",
        "be plotted and the second virial can be plotted as a function",
        "of temperature."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff",      "aff",     ffREAD  },
        { efXML, "-charges", "charges", ffOPTRD },
        { efLOG, "-g",       "b2",      ffWRITE } 
    };
    gmx_output_env_t         *oenv;
    static char              *molnm      = (char *)"";
    static char              *qqm        = (char *)"";
    static char              *filename   = (char *)"";
    double                    qtot       = 0;
    double                    shellToler = 1e-6;
    bool                      verbose    = false;
    bool                      json       = false;
    std::vector<t_pargs>      pa = {
        { "-f",      FALSE, etSTR,  {&filename},
           "Input file name" },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule." },
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Combined charge of the molecule(s). This will be taken from the input file by default, but that is not always reliable." },
        { "-qqm",    FALSE, etSTR,  {&qqm},
          "Use a method from quantum mechanics that needs to be present in the input file. Either ESP, Hirshfeld, CM5 or Mulliken may be available." },
        { "-v", FALSE, etBOOL, {&verbose},
          "Print more information to the log file." },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)." },
        { "-json", FALSE, etBOOL, {&json},
          "Print part of the output in json format." }
    };
    DimerGenerator gendimers;
    gendimers.addOptions(&pa, &fnm);
    ReRunner       rerun(true);
    rerun.addOptions(&pa, &fnm);
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        status = 1;
        return status;
    }
    gendimers.finishOptions();
    
    ForceField        pd;
    try
    {
        readForceField(opt2fn("-ff", fnm.size(), fnm.data()), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    (void) pd.verifyCheckSum(stderr);
    const char *logFileName = opt2fn("-g", fnm.size(),fnm.data());
    FILE *logFile   = gmx_ffopen(logFileName, "w");
    auto  forceComp = new ForceComputer(shellToler, 100);
    print_header(logFile, pa, fnm);
    
    JsonTree jtree("SecondVirialCoefficient");
    if (verbose)
    {
        forceFieldSummary(&jtree, &pd);
    }

    ACTMol                actmol;
    std::map<std::string, std::vector<double> > qmap;
    auto qfn       = opt2fn_null("-charges", fnm.size(), fnm.data());
    if (qfn)
    {
        qmap = fetchChargeMap(&pd, forceComp, qfn);
        fprintf(logFile, "\nRead %lu entries into charge map from %s\n", qmap.size(), qfn);
    }
    if (strlen(molnm) == 0)
    {
        molnm = (char *)"MOL";
    }
    {
        std::vector<MolProp>  mps;
        double qtot_babel = qtot;
        int maxpot = 100;
        int nsymm  = 1;
        std::string method, basis;
        const char *conf = "";
        const char *jobtype = (char *)"Opt";
        matrix box;
        bool   userqtot = opt2parg_bSet("-qtot", pa.size(), pa.data());
        if (readBabel(&pd, filename, &mps, molnm, molnm, conf, &method, &basis,
                      maxpot, nsymm, jobtype, userqtot, &qtot_babel, false, box))
        {
            if (mps.size() > 1)
            {
                fprintf(stderr, "Warning: will only use the first dimer in %s\n", filename);
            }
            actmol.Merge(&mps[0]);
        }
        else
        {
            gmx_fatal(FARGS, "No input file has been specified.");
        }
    }
    
    immStatus imm = immStatus::OK;
    if (status == 0)
    {
        imm = actmol.GenerateTopology(logFile, &pd, missingParameters::Error);
    }
    std::vector<gmx::RVec> coords = actmol.xOriginal();
    if (immStatus::OK == imm && status == 0)
    {
        auto fragments  = actmol.fragmentHandler();
        if (fragments->setCharges(qmap))
        {
            // Copy charges to the high-level topology as well
            fragments->fetchCharges(actmol.atoms());
        }
        else
        {
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());

            std::vector<double> myq;
            auto alg   = pd.chargeGenerationAlgorithm();
            auto qtype = qType::Calc;
            if (strlen(qqm) > 0)
            {
                alg   = ChargeGenerationAlgorithm::Read;
                qtype = stringToQtype(qqm);
            }
            fprintf(logFile, "WARNING: No information in charge map. Will generate charges using %s algorithm\n", chargeGenerationAlgorithmName(alg).c_str());
            imm    = actmol.GenerateCharges(&pd, forceComp, alg, qtype, myq, &coords, &forces);
        }
    }
    if (immStatus::OK == imm && status == 0)
    {
        if (debug)
        {
            actmol.topology()->dump(debug);
        }
        rerun.setFunctions(forceComp, &gendimers, oenv);
        bool userqtot = opt2parg_bSet("-qtot", pa.size(), pa.data());
        rerun.rerun(logFile, &pd, &actmol, userqtot, qtot, verbose, fnm);
    }
    if (json)
    {
        jtree.write("simulate.json", json);
    }
    else
    {
        jtree.fwrite(logFile, json);
    }
    gmx_ffclose(logFile);
    return status;
}

} // namespace alexandria
