/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

#include <random>

#include <ctype.h>
#include <stdlib.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/futil.h"

#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/poldata/poldata_xml.h"
#include "act/utility/jsontree.h"
#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "alexandria/alex_modules.h"
#include "alexandria/atype_mapping.h"
#include "alexandria/babel_io.h"
#include "alexandria/confighandler.h"
#include "alexandria/molhandler.h"
#include "alexandria/mymol.h"
#include "alexandria/princ.h"
#include "alexandria/tuning_utility.h"

namespace alexandria
{

void forceFieldSummary(JsonTree      *jtree,
                       const Poldata *pd)
{
    jtree->addObject(JsonTree("Force field file", pd->filename()));
    jtree->addObject(JsonTree("Created", pd->timeStamp()));
    jtree->addObject(JsonTree("Checksum", pd->checkSum()));
    jtree->addObject(JsonTree("Polarizable", yesno_names[pd->polarizable()]));
    jtree->addObject(JsonTree("Charge generation", 
                              chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str()));
    jtree->addObject(JsonTree("# exclusions", gmx_itoa(pd->getNexcl())));
    jtree->addObject(JsonTree("Relative dielectric constant epsilon_r",
                              gmx_ftoa(pd->getEpsilonR())));
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
                          std::vector<t_filenm> *filenm,
                          bool                   b2code)
{
    std::vector<t_pargs> pa = {
        { "-traj",   FALSE, etSTR,  {&trajname_},
          "Trajectory or series of structures of the same compound for which the energies will be computed. If this option is present, no simulation will be performed." },
        { "-T1",     FALSE, etREAL, {&T1_},
          "Starting temperature for second virial calculations." },
        { "-T2",     FALSE, etREAL, {&T2_},
          "Starting temperature for second virial calculations." },
        { "-dT",     FALSE, etREAL, {&deltaT_},
          "Temperature increment for calculation of second virial." }
    };
    pargs->push_back(pa[0]);
    if (b2code)
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
    if (b2t_.empty() || b2t_.size() != T.size())
    {
        fprintf(stderr, "Internal inconsistency. There are %zu temperatures and %zu B2 values.\n",
                T.size(), b2t_.size());
        return;
    }
    FILE *b2p = xvgropen(b2file, "Second virial coefficient",
                         "Temperature (K)", "B2(T) cm^3/mol", oenv_);
    for(size_t ii = 0; ii < T.size(); ii++)
    {
        fprintf(b2p, "%10g  %10g\n", T[ii], b2t_[ii]);
    }
    xvgrclose(b2p);
}

void ReRunner::computeB2(FILE                                      *logFile,
                         gmx_stats                                  edist,
                         double                                     mass,
                         const gmx::RVec                            inertia[2],
                         const std::vector<gmx::RVec>              &force1,
                         const std::vector<std::vector<gmx::RVec>> &torqueMol,
                         const std::vector<t_filenm>               &fnm)
{
    auto                      N = edist.get_npoints();
    const std::vector<double> x = edist.getX();
    const std::vector<double> y = edist.getY();
    double xmin                 = *std::min_element(x.begin(), x.end());
    double xmax                 = *std::max_element(x.begin(), x.end());

    if (N > 2 && xmax > xmin)
    {
        real   binwidth = 0.01; // nm
        // Bins start from zero for proper integration
        size_t nbins    = 1+std::round(xmax/binwidth);
        binwidth        = xmax/(nbins-1);
        // Temp array to store distance and Mayer functions
        auto Temperature = temperatures();
        std::vector<std::vector<double> > mayer(1+Temperature.size());
        // Array for total second virial as a function of T
        size_t iTemp = 1;
        for(auto T : Temperature)
        {
            if (T == 0)
            {
                fprintf(stderr, "Please provide a finite temperature to compute second virial.\n");
                continue;
            }
            // Temporary arrays for weighted properties.
            std::vector<double>    exp_U12(nbins, 0.0);
            std::vector<double>    exp_F2(nbins, 0.0);
            std::vector<gmx::RVec> exp_tau[2];
            for(int kk = 0; kk < 2; kk++)
            {
                exp_tau[kk].resize(nbins, { 0.0, 0.0, 0.0 });
            }
            std::vector<int>       n_U12(nbins, 0);
            double beta = 1.0/(BOLTZ*T);
            for(size_t ii = 0; ii < x.size(); ii++)
            {
                double rindex = x[ii]/binwidth;
                size_t index  = rindex;
                // Gray and Gubbins Eqn. 3.261
                double g0_12 = std::exp(-y[ii]*beta);
                // Gray and Gubbins Eqn. 3.272
                exp_U12[index] += g0_12-1;
                // Gray and Gubbins Eqn. 3.281
                exp_F2[index]  += g0_12*iprod(force1[ii], force1[ii]);
                for(int m = 0; m < DIM; m++)
                {
                    // Gray and Gubbins Eqn. 3.282
                    for(int kk = 0; kk < 2; kk++)
                    {
                        exp_tau[kk][index][m] += g0_12*torqueMol[kk][ii][m]*torqueMol[kk][ii][m];
                    }
                }
                n_U12[index]   += 1;
            }
            double    Bclass    =  0;
            double    BqmForce  =  0;
            double    BqmTorque =  0;
            // We start in the origin even if there is no data.
            double    r1        =  0;
            // Starting energy, all values until first data entry
            double    Uprev     = -1;
            int jj = 0;
            while(jj*binwidth < xmin)
            {
                exp_U12[jj] = -1;
                n_U12[jj]   = 1;
                jj += 1;
            }
            // Starting force
            double    Fprev     =  0;
            // Starting torque
            gmx::RVec Tprev[2]  = { { 0, 0, 0 }, { 0, 0, 0 } };
            double hbarfac      = beta*gmx::square(PLANCK*beta/(2*M_PI))/24;
            if (iTemp == 1)
            {
                // Store distance first time around only
                mayer[0].push_back(r1);
            }
            mayer[iTemp].push_back(Uprev);
            for(size_t ii = 1; ii < nbins; ii++)
            {
                double r2 = ii*binwidth;
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
                    // TODO: There is factor 0.5 here, but results are off by a factor two.
                    Bclass       -= 0.5*dB;
                    Uprev         = Unew;
                    // Weighted square force
                    // We follow Eqn. 9 in Schenter, JCP 117 (2002) 6573
                    double Fnew   = 0.5*exp_F2[ii]/(mass*n_U12[ii]);
                    BqmForce     += hbarfac*sphereIntegrator(r1, r2, Fprev, Fnew);
                    Fprev         = Fnew;
                    // Contributions from torque
                    for(int m = 0; m < DIM; m++)
                    {
                        for(int kk = 0; kk < 2; kk++)
                        {
                            if (inertia[kk][m] > 0)
                            {
                                double Tnew  = exp_tau[kk][ii][m]/(n_U12[ii]*inertia[kk][m]);
                                BqmTorque   += 0.5*hbarfac*sphereIntegrator(r1, r2, Tprev[kk][m], Tnew);
                                Tprev[kk][m] = Tnew;
                            }
                        }
                    }
                    r1 = r2;
                }
            }
            // Conversion to regular units cm^3/mol.
            double fac  = AVOGADRO*1e-21;
            double Btot = (Bclass + BqmForce + BqmTorque)*fac;
            b2t_.push_back(Btot);
            if (logFile)
            {
                fprintf(logFile, "T = %g K. Classical second virial coefficient B2cl %g BqmForce %g BqmTorque %g Total %g cm^3/mol\n", T,
                        Bclass*fac, BqmForce*fac, BqmTorque*fac, Btot);
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

class Rotator
{
private:
    std::random_device                     rd_;
    std::mt19937                           gen_;
    std::uniform_real_distribution<double> dis_;
public:
    Rotator(int seed) : gen_(rd_()), dis_(std::uniform_real_distribution<double>(0.0, 1.0))
    {
        if (seed > 0)
        {
            gen_.seed(seed);
        }
    }

    void random(std::vector<gmx::RVec> *coords)
    {
        // Distribution is 0-M_PI, multiply by two to get to 2*M_PI
        double alpha = dis_(gen_) * 2 * M_PI;
        double beta  = dis_(gen_) * 2 * M_PI; 
        double gamma = std::acos(2*dis_(gen_)-1);
        double cosa  = std::cos(alpha);
        double sina  = std::sin(alpha);
        double cosb  = std::cos(beta);
        double sinb  = std::sin(beta);
        double cosg  = std::cos(gamma);
        double sing  = std::sin(gamma);
        
        matrix A;
        A[0][0] = cosb * cosg;
        A[0][1] =-cosb * sing;
        A[0][2] = sinb;
        
        A[1][0] = sina * sinb * cosg + cosa * sing;
        A[1][1] =-sina * sinb * sing + cosa * cosg;
        A[1][2] =-sina * cosb;
        
        A[2][0] =-cosa * sinb * cosg + sina * sing;
        A[2][1] = cosa * sinb * sing + sina * cosg;
        A[2][2] = cosa * cosb;
        
        auto oldX = *coords;
        for(size_t i = 0; i < oldX.size(); i++)
        {
            gmx::RVec newx;
            mvmul(A, oldX[i], newx);
            copy_rvec(newx, (*coords)[i]);
        }
    }
};

void DimerGenerator::addOptions(std::vector<t_pargs>  *pa,
                                std::vector<t_filenm> *fnm)
{
    std::vector<t_pargs> mypa = {
        { "-maxdimer", FALSE, etINT, {&maxdimers_},
          "Number of dimer orientations to generate if you do not provide a trajectory. For each of these a distance scan will be performed." },
        { "-ndist", FALSE, etINT, {&ndist_},
          "Number of distances to use for computing interaction energies and forces. Total number of dimers is the product of maxdimer and ndist." },
        { "-mindist", FALSE, etREAL, {&mindist_},
          "Minimum com-com distance to generate dimers for." },
        { "-maxdist", FALSE, etREAL, {&maxdist_},
          "Maximum com-com distance to generate dimers for." },
        { "-seed", FALSE, etINT, {&seed_},
          "Random number seed to generate monomer orientations, applied if seed is larger than 0. If not, the built-in default will be used." }
    };
    for(auto &pp : mypa)
    {
        pa->push_back(pp);
    }
    std::vector<t_filenm>  myfnm = {
        { efSTX, "-ox", "dimers",     ffOPTWR  }
    };
    for(auto &ff : myfnm)
    {
        fnm->push_back(ff);
    }
}
    
void DimerGenerator::generate(FILE                                *logFile,
                              const MyMol                         *mymol,
                              std::vector<std::vector<gmx::RVec>> *coords,
                              const char                          *outcoords)
{
    auto fragptr = mymol->fragmentHandler();
    if (fragptr->topologies().size() == 2)
    {
        // Random number generation
        Rotator rot(seed_);
        
        // Copy original coordinates
        auto xorig     = mymol->xOriginal();
        // Split the coordinates into two fragments
        auto atomStart = fragptr->atomStart();
        std::vector<gmx::RVec> xmOrig[2];
        for(int m = 0; m < 2; m++)
        {
            for(size_t j = atomStart[m]; j < atomStart[m+1]; j++)
            {
                xmOrig[m].push_back(xorig[j]);
            } 
        }
        // Topologies
        auto tops = fragptr->topologies();
        // Move molecules to their respective COM
        gmx::RVec com[2];
        for(int m = 0; m < 2; m++)
        {
            // Compute center of mass
            clear_rvec(com[m]);
            auto   atoms   = tops[m].atoms();
            double totmass = 0;
            for(size_t j = 0; j < atoms.size(); j++)
            {
                gmx::RVec mx;
                svmul(atoms[j].mass(), xmOrig[m][j], mx);
                rvec_inc(com[m], mx);
                totmass += atoms[j].mass();
            }
            for(int n = 0; n < DIM; n++)
            {
                com[m][n] /= totmass;
            }
            // Subtract center of mass
            for(size_t j = 0; j < atoms.size(); j++)
            {
                rvec_sub(xmOrig[m][j], com[m], xmOrig[m][j]);
            }
        }
        // Loop over orientations
        size_t nmp = maxdimers_*ndist_;
        // Initiate all the MolProps with a copy of the input molecule
        MolProp tmp = *mymol;
        auto exper  = tmp.experiment();
        // Do some cleaning
        exper->clear();
        tmp.clearCategory();
        // Then initiate the big array
        coords->resize(nmp);
        if (logFile)
        {
            print_memory_usage(logFile);
        }
        // Shortcut to the current molecules
        // MolProp *mp = &((*mps)[mp_index++]);
        // exper = mp->experiment();
        for(int ndim = 0; ndim < maxdimers_; ndim++)
        {
            // Copy the coordinates and rotate them
            std::vector<gmx::RVec> xrand[2];
            for(int m = 0; m < 2; m++)
            {
                xrand[m] = xmOrig[m];
                // Random rotation
                rot.random(&xrand[m]);
            }
            // Loop over distances from mindist to maxdist
            double range = maxdist_ - mindist_;
            for(int idist = 0; idist < ndist_; idist++)
            {
                double    dist  = mindist_ + range*((1.0*idist)/(ndist_-1));
                gmx::RVec trans = { 0, 0, dist };
                auto      atoms = tops[1].atoms();
                for(size_t j = 0; j < atoms.size(); j++)
                {
                    rvec_inc(xrand[1][j], trans);
                }
                size_t idim = ndim*ndist_+idist;
                (*coords)[idim].resize(atomStart[2]-atomStart[0]);
                for(int m = 0; m < 2; m++)
                {
                    for(size_t j = atomStart[m]; j < atomStart[m+1]; j++)
                    {
                        auto jindex = j - atomStart[m];
                        copy_rvec(xrand[m][jindex], (*coords)[idim][j]);
                    }
                }
                // Put the coordinates back!
                for(size_t j = 0; j < atoms.size(); j++)
                {
                    rvec_dec(xrand[1][j], trans);
                }
            }
            if (logFile)
            {
                print_memory_usage(logFile);
            }
        }
    }
}

void ReRunner::rerun(FILE                        *logFile,
                     const Poldata               *pd,
                     const MyMol                 *mymol,
                     double                       qtot,
                     bool                         verbose,
                     const std::vector<t_filenm> &fnm)
{
    std::vector<std::vector<gmx::RVec> > dimers;
    std::string          method, basis;
    int                  maxpot = 100;
    int                  nsymm  = 1;
    const char          *molnm  = "";
    if (verbose)
    {
        print_memory_usage(logFile);
    }
    if (trajname_ && strlen(trajname_) > 0)
    {
        std::vector<MolProp> mps;
        // Read compounds
        if (!readBabel(trajname_, &mps, molnm, molnm, "", &method,
                       &basis, maxpot, nsymm, "Opt", &qtot, false))
        {
            fprintf(stderr, "Could not read compounds from %s\n", trajname_);
            return;
        }
        if (logFile)
        {
            fprintf(logFile, "Doing energy calculation for %zu structures from %s\n",
                    mps.size(), trajname_);
        }
        dimers.resize(mps.size());
        for(size_t i = 0; i < mps.size(); i++)
        {
            auto exper = mps[i].experimentConst();
            dimers[i] = exper[0].getCoordinates();
        }
    }
    else
    {
        // Generate compounds
        gendimers_->generate(verbose ? logFile : nullptr, mymol, &dimers, nullptr);
        if (logFile)
        {
            fprintf(logFile, "Doing energy calculation for %zu randomly oriented structures generated at %d distances\n",
                    dimers.size(), gendimers_->ndist());
        }
    }
    if (verbose)
    {
        print_memory_usage(logFile);
    }
    std::map<InteractionType, double> energies;
    int mp_index = 0;
    gmx_stats edist;
    std::vector<gmx::RVec> force1;
    std::vector<std::vector<gmx::RVec>> torqueMol;
    torqueMol.resize(2);
    if (eInter_)
    {
        for(int kk = 0; kk < 2; kk++)
        {
            torqueMol[kk].resize(dimers.size());
        }
        force1.resize(dimers.size());
    }
    if (verbose)
    {
        print_memory_usage(logFile);
    }
    gmx::RVec inertia[2]  = { { 0, 0, 0 }, { 0, 0, 0 } };
    // Loop over molecules
    const auto &atoms = mymol->atomsConst();
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
                    GMX_THROW(gmx::InvalidInputError("Number of generated coordinates in trajectory does not match molecule file"));
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
            auto EE         = mymol->calculateInteractionEnergy(pd, forceComp_, &forces, &coords);
            auto atomStart  = mymol->fragmentHandler()->atomStart();
            if (atomStart.size() != 3)
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("This is not a dimer, there are %zu fragments instead of 2", atomStart.size()-1).c_str()));
            }
            gmx::RVec f[2]    = { { 0, 0, 0 }, { 0, 0, 0 } };
            gmx::RVec com[2]  = { { 0, 0, 0 }, { 0, 0, 0 } };
            double    mtot[2] = { 0, 0 };
            auto      tops  = mymol->fragmentHandler()->topologies();
            for(int kk = 0; kk < 2; kk++)
            {
                auto atoms = tops[kk].atoms();
                for(size_t i = atomStart[kk]; i < atomStart[kk+1]; i++)
                {
                    gmx::RVec mr1;
                    auto      mi = atoms[i-atomStart[kk]].mass();
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
            }
            copy_rvec(f[0], force1[mp_index]);
            gmx::RVec torque[2]     = { { 0, 0, 0 }, { 0, 0, 0 } };
            gmx::RVec torqueRot[2]  = { { 0, 0, 0 }, { 0, 0, 0 } };
            for(int kk = 0; kk < 2; kk++)
            {
                // Compute the coordinates relative to the center of mass
                std::vector<real>      mass;
                std::vector<int>       index;
                std::vector<gmx::RVec> x_com;
                for(size_t i = atomStart[kk]; i < atomStart[kk+1]; i++)
                {
                    // Store atom index relative to molecule start
                    index.push_back(i-atomStart[kk]);
                    // Store mass
                    mass.push_back(atoms[i].mass());
                    // Subtract COM and store
                    gmx::RVec ri;
                    rvec_sub(coords[i], com[0], ri);
                    x_com.push_back(ri);
                }
                // Compute moments of inertia and transformation matrix
                rvec   inertia1;
                matrix trans;
                principal_comp(index.size(), index.data(), mass.data(), 
                               as_rvec_array(x_com.data()),
                               trans, inertia1);
                // Move to inertial frame (only well-defined for
                // rigid molecules).
                // The trans matrix should convert that coordinate to the inertial frame,
                // but what about the force on the atoms? It likely has to be rotated in the
                // same manner. After that, the torque can be computed.
                // Before doing the rotations, the sum of the torque vectors is zero.
                // Since the rotations are different for both molecules, this does not
                // hold after the rotations are done.
                // TODO: write out the math.
                for(size_t i = atomStart[kk]; i < atomStart[kk+1]; i++)
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
            if (verbose)
            {
                print_memory_usage(logFile);
            }
            gmx::RVec dcom;
            rvec_sub(com[0], com[1], dcom);
            double rcom = norm(dcom);
            if (verbose)
            {
                fprintf(logFile, " r %g Einter %g Force %g %g %g Torque[0] %g %g %g Torque[1] %g %g %g Rotated Torque[0] %g %g %g Rotated Torque[1] %g %g %g",
                        rcom, EE, f[0][XX], f[0][YY], f[0][ZZ],
                        torque[0][XX], torque[0][YY], torque[0][ZZ],
                        torque[1][XX], torque[1][YY], torque[1][ZZ],
                        torqueRot[0][XX], torqueRot[0][YY], torqueRot[0][ZZ],
                        torqueRot[1][XX], torqueRot[1][YY], torqueRot[1][ZZ]);
            }
            edist.add_point(rcom, EE, 0, 0);
        }
        else
        {
            forceComp_->compute(pd, mymol->topology(),
                                &coords, &forces, &energies);
            for(const auto &ee : energies)
            {
                fprintf(logFile, "  %s %8g", interactionTypeToString(ee.first).c_str(), ee.second);
            }
        }
        if (verbose)
        {
            fprintf(logFile, "\n");
        }
        mp_index++;
    }
    if (verbose)
    {
        print_memory_usage(logFile);
    }
    if (eInter_)
    {
        for(int kk = 0; kk < 2; kk++)
        {
            for(int m = 0; m < DIM; m++)
            {
                inertia[kk][m] /= edist.get_npoints();
            }
        }
        // Compute the relative mass
        double m0 = mymol->fragmentHandler()->topologies()[0].mass();
        double m1 = mymol->fragmentHandler()->topologies()[1].mass();
        double mm = m0*m1/(m0+m1);
        computeB2(logFile, edist, mm, inertia, force1, torqueMol, fnm);
    }
    print_memory_usage(stdout);
}

int b2(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "alexandria b2 will either read a trajectory",
        "of dimers as input for energy calculations or generate",
        "conformations of dimers. The",
        "corresponding molecule file, used for generating the topology",
        "needs to be a molprop (xml) file and contain information about",
        "the compounds in the dimer. Based on dimer energies the second",
        "virial coefficient will be estimated. Mayer curves can optionally",
        "be plotted and the second virial can be plotted as a function",
        "of temperature."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff", "gentop",     ffREAD  },
        { efXML, "-mp", "molprop",    ffOPTRD },
        { efLOG, "-g",  "b2",         ffWRITE } 
    };
    gmx_output_env_t         *oenv;
    static char              *molnm      = (char *)"";
    static char              *qqm        = (char *)"";
    double                    qtot       = 0;
    double                    shellToler = 1e-6;
    bool                      verbose    = false;
    bool                      json       = false;
    std::vector<t_pargs>      pa = {
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
    ReRunner       rerun;
    rerun.addOptions(&pa, &fnm, true);
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        status = 1;
        return status;
    }
    Poldata        pd;
    try
    {
        readPoldata(opt2fn("-ff", fnm.size(), fnm.data()), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    (void) pd.verifyCheckSum(stderr);
    const char *logFileName = opt2fn("-g", fnm.size(),fnm.data());
    FILE *logFile   = gmx_ffopen(logFileName, "w");
    auto  forceComp = new ForceComputer(shellToler, 100);
    print_header(logFile, pa);
    
    JsonTree jtree("SecondVirialCoefficient");
    if (verbose)
    {
        forceFieldSummary(&jtree, &pd);
    }

    MyMol                mymol;
    {
        std::vector<MolProp> mps;
        MolPropRead(opt2fn("-mp", fnm.size(), fnm.data()), &mps);
        if (mps.size() > 1)
        {
            fprintf(stderr, "Warning: will only use the first compound (out of %zu) in %s\n", mps.size(), 
                    opt2fn("-mp", fnm.size(), fnm.data()));
        }
        mymol.Merge(&mps[0]);
    }
    immStatus imm = immStatus::OK;
    if (status == 0)
    {
        imm = mymol.GenerateTopology(logFile, &pd, missingParameters::Error, false);
    }
    std::vector<gmx::RVec> coords = mymol.xOriginal();
    if (immStatus::OK == imm && status == 0)
    {
        CommunicationRecord cr;
        gmx::MDLogger  mdlog {};
        std::vector<gmx::RVec> forces(mymol.atomsConst().size());

        std::vector<double> myq;
        auto alg   = pd.chargeGenerationAlgorithm();
        auto qtype = qType::Calc;
        if (strlen(qqm) > 0)
        {
            alg   = ChargeGenerationAlgorithm::Read;
            qtype = stringToQtype(qqm);
        }
        imm    = mymol.GenerateCharges(&pd, forceComp, mdlog, &cr, alg, qtype, myq, &coords, &forces);
    }
    if (immStatus::OK == imm && status == 0)
    {
        if (debug)
        {
            mymol.topology()->dump(debug);
        }
        rerun.setFunctions(forceComp, &gendimers, oenv);
        rerun.rerun(logFile, &pd, &mymol, qtot, verbose, fnm);
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
