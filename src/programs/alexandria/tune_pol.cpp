/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <random>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/strconvert.h"

#include "alex_modules.h"
#include "composition.h"
#include "molprop.h"
#include "molprop_tables.h"
#include "molprop_xml.h"
#include "molselect.h"
#include "poldata.h"
#include "poldata_xml.h"

class pType
{
    private:
    //! Poltype name
        std::string         name_;
        bool                bUse_;
    //! Number of data points corresponding to this poltype
        int                 nCopies_;
     //!Vector of polarizability values resulted from bootstrapping
        std::vector<double> polvalues_;
    //! Statistics structure for all data points
        gmx_stats_t         fullStats_;
    //! Statistics structure for data points between the %95 percentile intervals
        gmx_stats_t         improvedStats_;
    public:
    //! Constructor
        pType(std::string name, const bool bUse, const int nCopies);
    //! Is this one used
        bool bUse() { return bUse_; }
    //! Set the flag for being used
        void setUse(bool bUse) { bUse_ = bUse; }
    //! Check whether we have sufficient data to use this one
        void checkUse(int mindata) { bUse_ = (nCopies_ > mindata); }
    //! Return the name of the poltype
        std::string name() { return name_; }
    //! Return the number of this poltype in all molecules
        int nCopies() { return nCopies_; }
    //! Set the number of copies to 0
        void resetCopies() { nCopies_ = 0; }
    //! Increase the number of copies for this poltype
        void incCopies() { nCopies_++; }
    //! Return the statistics structure (which is a pointer)
        gmx_stats_t fullStats() { return fullStats_; }
        gmx_stats_t improvedStats() { return improvedStats_; }
    //! Add polvalue
        void addPolValue(double v) {polvalues_.push_back(std::move(v));}
    //! Improve stats by Percentile Bootstrap Intervals
        void improveStats();
};

pType::pType(std::string name, const bool bUse, const int nCopies)
{
    name_          = name;
    bUse_          = bUse;
    nCopies_       = nCopies;
    fullStats_     = gmx_stats_init();
    improvedStats_ = gmx_stats_init();
}

void pType::improveStats()
{   
    /*Percentile Bootstrap Intervals
      Bootstrapping Regression Models.
      Chapter 21.*/
      
    int lower = 0.025 * polvalues_.size();
    int upper = 0.975 * polvalues_.size();  
       
    std::sort(polvalues_.begin(), polvalues_.end());          
    for (int i = lower; i <= upper; i++)
    {
        gmx_stats_add_point_ydy(improvedStats_, polvalues_[i], 0);
    }
}

static bool check_matrix(double             **a, 
                         double              *x, 
                         size_t               nrow,
                         std::vector<pType>  &ptypes)
{
    int nrownew = nrow;
    for (size_t i = 0; i < ptypes.size(); i++)
    {
        for (size_t j = i+1; j < ptypes.size(); j++)
        {
            bool bSame = true;
            for (size_t k = 0; bSame && (k < nrow); k++)
            {
                bSame = (a[k][i] == a[k][j]);
            }
            if (bSame)
            {
                return false;
                gmx_fatal(FARGS, "Columns %d (%s) and %d (%s) are linearly dependent",
                          static_cast<int>(i), ptypes[i].name().c_str(), 
                          static_cast<int>(j), ptypes[j].name().c_str());
            }
        }
    }
    //! Check diagonal
    if (nrow == ptypes.size())
    {
        for (size_t i = 0; i < ptypes.size(); i++)
        {
            if (a[i][i] == 0)
            {
                return false;
                gmx_fatal(FARGS, "a[%d][%d] = 0. Pol type = %s",
                          static_cast<int>(i),  static_cast<int>(i), ptypes[i].name().c_str());
            }
        }
    }
    return true;
    for (auto i = 0; i < nrownew; i++)
    {
        for (auto j = i+1; j < nrownew; j++)
        {
            bool bSame = true;
            for (size_t k = 0; bSame && (k < ptypes.size()); k++)
            {
                bSame = (a[i][k] == a[j][k]);
            }
            if (bSame)
            {
                fprintf(stderr, "Rows %d and %d are linearly dependent. Removing %d.\n",
                        i, j, j);
                if (j < nrownew-1)
                {
                    for (size_t k = 0; k < ptypes.size(); k++)
                    {
                        a[j][k] = a[nrownew-1][k];
                    }
                    x[j] = x[nrownew-1];
                }
                nrownew--;
            }
        }
    }
    return nrownew;
}

static bool bZeroPol(const char *ptype, std::vector<std::string> zeropol)
{
    for (auto si = zeropol.begin(); (si < zeropol.end()); si++)
    {
        if (si->compare(ptype) == 0)
        {
            return true;
        }
    }
    return false;
}

static void fill_matrices_and_dump_csv(const alexandria::Poldata        &pd,
                                       std::vector<alexandria::MolProp> &mp,
                                       const alexandria::MolSelect      &gms,
                                       std::vector<pType>               &ptypes,
                                       int                               nusemol,
                                       double                            x[],
                                       double                          **a,
                                       double                          **at)
{
    alexandria::CompositionSpecs cs;
    FILE *csv = gmx_ffopen("tune_pol.csv", "w");

    fprintf(csv, "\"molecule\",\"formula\",");
    for (size_t i = 0; i < ptypes.size(); i++)
    {
        fprintf(csv, "\"%zu %s\",", i, ptypes[i].name().c_str());
    }
    fprintf(csv, "\"polarizability\"\n");    
    int nn = 0;
    int j  = 0;
    for (auto mpi = mp.begin(); mpi < mp.end(); mpi++, j++)
    {
        iMolSelect ims = gms.status(mpi->getIupac());
        if (imsIgnore != ims)
        {
            auto mci = mpi->SearchMolecularComposition(cs.searchCS(alexandria::iCalexandria)->name());            
            fprintf(csv, "\"%d %s\",\"%s\",", nn, mpi->getMolname().c_str(), mpi->formula().c_str());
            
            std::vector<int> count;
            count.resize(ptypes.size());
            for (auto ani = mci->BeginAtomNum(); ani < mci->EndAtomNum(); ani++)
            {
                std::string  ptype;
                if (pd.atypeToPtype(ani->getAtom(), ptype))
                {
                    size_t i;
                    for (i = 0; i < ptypes.size(); i++)
                    {
                        if (strcmp(ptype.c_str(), ptypes[i].name().c_str()) == 0)
                        {
                            break;
                        }
                    }
                    if (i < ptypes.size())
                    {
                        count[i] += ani->getNumber();
                    }
                }
                else if (nullptr != debug)
                {
                    fprintf(debug, "Supported molecule %s has unsupported or zeropol atom %s (ptype %s)",
                            mpi->getMolname().c_str(), ani->getAtom().c_str(), ptype.c_str());
                }
            }
            for (size_t i = 0; i < ptypes.size(); i++)
            {
                a[nn][i] = at[i][nn] = count[i];
                fprintf(csv, "%d,", count[i]);
            }
            fprintf(csv, "%.3f\n", x[nn]);
            nn++;
        }
    }
    fclose(csv);
    if (nusemol != nn)
    {
        gmx_fatal(FARGS, "Consistency error: nusemol = %d, nn = %d", nusemol, nn);
    }
}

static bool ptypes_still_have_support(const alexandria::Poldata        &pd,
                                      std::vector<alexandria::MolProp> &mp,
                                      std::vector<pType>               &ptypes,
                                      const char                       *alexandria,
                                      int                               mindata)
{
    int    found    = 0;
    size_t nsupport = 0;
    for(auto &ptype : ptypes)
    {
        found = 0;
        for (auto mpi = mp.begin(); (mpi < mp.end() && found < mindata); ++mpi)
        {
            auto mci = mpi->SearchMolecularComposition(alexandria);
            for (auto ani = mci->BeginAtomNum(); ani < mci->EndAtomNum(); ++ani)
            {
                std::string p;
                if (pd.atypeToPtype(ani->getAtom(), p))
                {
                    if (p == ptype.name())
                    {
                       found++;
                       break;
                    }
                }
            }
        }
        if (found >= mindata)
        {
            nsupport++;
        }
    }
    if (nsupport == ptypes.size())
    {
        return true;
    }
    else
    {
        return false;
    }
}

static int decompose_frag(FILE                             *fplog,
                          const char                       *hisfn,
                          alexandria::Poldata              &pd,
                          std::vector<alexandria::MolProp> &mp,
                          gmx_bool                          bQM, 
                          const char                       *lot,
                          int                               mindata,
                          const alexandria::MolSelect      &gms,
                          gmx_bool                          bZero,
                          gmx_bool                          bForceFit,
                          int                               nBootStrap, 
                          std::vector<std::string>          zeropol,
                          const gmx_output_env_t           *oenv,
                          const char                       *exp_type)
{

    int                          j         = 0;
    int                          niter     = 0;
    int                          row       = 0;
    int                          nusemol   = 0;
    int                          nAccepted = 0;
    double                       poltot    = 0;
    double                       a0        = 0;
    double                       da0       = 0;
    double                       ax        = 0;
    double                       chi2      = 0;
    double                     **a;
    double                     **at;
    double                     **ata; 
    std::vector<double>          x, atx, fpp;
    std::vector<pType>           ptypes;
    
    alexandria::CompositionSpecs cs;
    const char                  *alexandria = cs.searchCS(alexandria::iCalexandria)->name();

    x.resize(mp.size()+1);
    
    // Copy all the polarizability types into array.
    for (auto ptype = pd.getPtypeBegin(); ptype != pd.getPtypeEnd(); ptype++)
    {       
        if (!bZeroPol(ptype->getType().c_str(), zeropol))
        {
            ptypes.push_back(pType(ptype->getType().c_str(), false, 0));
        }
    }
    size_t ptstart = ptypes.size();
    size_t ptsize;
    int    iter      = 1;
    size_t nmol_orig = mp.size();
    // Remove the molecule that do not support none of the poltypes. 
    do
    {
        ptsize = ptypes.size();
        if (nullptr != fplog)
        {
            fprintf(fplog, "iter %d %zu ptypes left\n", iter++, ptsize);
        }
        nusemol  = 0;
        poltot   = 0;
        for (auto mpi = mp.begin(); mpi < mp.end(); )
        {
            iMolSelect ims  = gms.status(mpi->getIupac());
            double T        = -1;
            double pol      = 0;
            double sig_pol  = 0;
            bool   bPol     = mpi->getProp(MPO_POLARIZABILITY,
                                           bQM ? iqmBoth : iqmExp,
                                           lot, 
                                           "", 
                                           exp_type,
                                           &pol, 
                                           &sig_pol, 
                                           &T);
                                           
            auto mci              = mpi->SearchMolecularComposition(alexandria);
            bool bHaveComposition = (mci != mpi->EndMolecularComposition());
            bool bUseMol          = ((imsIgnore != ims) && bPol && (pol > 0) && bHaveComposition);
            if (nullptr != fplog)
            {
                fprintf(fplog, "%s nExper %d pol %g Use: %s\n", 
                        mpi->getMolname().c_str(), mpi->NExperiment(),
                        pol, gmx::boolToString(bUseMol));
            }

            // Check whether all atoms have supported ptypes.
            size_t npolarizable = 0;
            for (auto ani = mci->BeginAtomNum(); (bUseMol && (ani < mci->EndAtomNum())); ++ani)
            {
                const char *atomname = ani->getAtom().c_str();
                std::string ptype;
                if (!pd.atypeToPtype(atomname, ptype))
                {
                    if (nullptr != fplog)
                    {
                        fprintf(fplog, "No ptype for atom %s in molecule %s\n",
                                atomname, mpi->getMolname().c_str());
                    }
                    bUseMol = false;
                }
                else if (!bZeroPol(ptype.c_str(), zeropol))
                {
                    npolarizable++;
                }
            }
            if (nullptr != fplog)
            {
                fprintf(fplog, "Mol: %s, IUPAC: %s, ims: %s, bPol: %s, pol: %g - %s\n",
                        mpi->getMolname().c_str(),
                        mpi->getIupac().c_str(),
                        iMolSelectName(ims), gmx::boolToString(bPol), pol,
                        bUseMol ? "Used" : "Not used");
            }

            if (bUseMol && (npolarizable > 0))
            {
                x[nusemol] = pol;
                poltot    += pol;
                nusemol++;
                mpi++;
            }
            else
            {
                fprintf(fplog, "Removing %s. bPol = %s pol = %g composition found = %s npolarizable = %zu\n",
                        mpi->getMolname().c_str(), gmx::boolToString(bPol), pol,
                        (mci == mpi->EndMolecularComposition()) ? "true" : "false", npolarizable);
                mpi = mp.erase(mpi);
            }
        }
        // Check whether the ptypes still have support in experimental or QM polarizabilities.
        for (auto pi = ptypes.begin(); pi < ptypes.end();)
        {
            pi->resetCopies();
            double pol, sig_pol;
            if ((pd.getPtypePol(pi->name(), &pol, &sig_pol)) &&
                ((pol == 0) || bForceFit))
            {
                for (auto &mpi : mp)
                {
                    auto mci = mpi.SearchMolecularComposition(alexandria);
                    for (auto ani = mci->BeginAtomNum(); ani < mci->EndAtomNum(); ++ani)
                    {
                        std::string p;
                        if (pd.atypeToPtype(ani->getAtom(), p))
                        {
                            if (p == pi->name())
                            {
                                for (int i = 0; i < ani->getNumber(); i++)
                                {
                                    pi->incCopies();
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                if (nullptr != fplog)
                {
                    fprintf(fplog, "Polarizability for %s is %g. Not optimizing this one.\n",
                            pi->name().c_str(), pol);
                }
            }
            pi->checkUse(mindata);
            if (pi->bUse())
            {
                if (nullptr != fplog)
                {
                    fprintf(fplog, "Will optimize polarizability for ptype %s with %d copies\n",
                            pi->name().c_str(), pi->nCopies());
                }
                pi++;
            }
            else
            {            
                pi = ptypes.erase(pi);
            }
        }
    }
    while (ptypes.size() < ptsize);
    if (mp.size() < nmol_orig)
    {
        printf("Reduced number of molecules from %zu to %zu\n", nmol_orig, mp.size());
    }
    if (ptypes.size() == 0)
    {
        printf("No polarization types to optimize, but before compaction %zu.\n", ptstart);
        return 0;
    }
    
    //Fill the matrix 'a' and the transposed matrix 'at'.
    a      = alloc_matrix(nusemol, ptypes.size());
    at     = alloc_matrix(ptypes.size(), nusemol);
    if (nullptr != fplog)
    {
        fprintf(fplog, "There are %d different atomtypes to optimize the polarizabilities\n", (int)ptypes.size());
        fprintf(fplog, "There are %d (experimental) reference polarizabilities.\n", nusemol);
        fflush(fplog);
    }
    fill_matrices_and_dump_csv(pd, mp, gms, ptypes, nusemol, x.data(), a, at);
    if (!check_matrix(a, x.data(), nusemol, ptypes))
    {
        fprintf(stderr, "Matrix is linearly dependent. Sorry!.\n");
    }

    // Bootstrap
    std::random_device              rd;
    std::mt19937                    gen(rd());
    std::uniform_int_distribution<> dis(0, nusemol-1);   
    for (int kk = 0; kk < nBootStrap; kk++)
    {
        fprintf(stderr, "\rBootStrap %d", 1+kk);
        ata    = alloc_matrix(ptypes.size(), ptypes.size());

        double            **a_copy  = alloc_matrix(nusemol, ptypes.size());
        double            **at_copy = alloc_matrix(ptypes.size(), nusemol);
        
        std::vector<alexandria::MolProp> mp_copy;
        std::vector<double> x_copy;
        x_copy.resize(nusemol);       
        do
        {
            mp_copy.clear();
            for (int ii = 0; ii < nusemol; ii++)
            {
                // Pick random molecule uu out of the stack.
                int uu = dis(gen);
                for (size_t jj = 0; jj < ptypes.size(); jj++)
                {
                    // Make row ii equal to uu in the original matrix
                    a_copy[ii][jj] = at_copy[jj][ii] = a[uu][jj];
                }
                // Make experimenal (QM) value equal to the original
                x_copy[ii] = x[uu]; 
                mp_copy.push_back(mp[uu]);
            }
        }
        while(!ptypes_still_have_support(pd, mp_copy, ptypes, alexandria, mindata));
               
        matrix_multiply(debug, nusemol, ptypes.size(), a_copy, at_copy, ata);        
        if (check_matrix(ata, x_copy.data(), ptypes.size(), ptypes) &&
            ((row = matrix_invert(debug, ptypes.size(), ata)) == 0))
        {
            atx.resize(ptypes.size());
            fpp.resize(ptypes.size());
            a0 = 0;
            nAccepted++;
            do
            {
                for (size_t i = 0; i < ptypes.size(); i++)
                {
                    atx[i] = 0;
                    for (j = 0; j < nusemol; j++)
                    {
                        atx[i] += at_copy[i][j]*(x_copy[j]-a0);
                    }
                }
                for (size_t i = 0; i < ptypes.size(); i++)
                {
                    fpp[i] = 0;
                    for (size_t j = 0; j < ptypes.size(); j++)
                    {
                        fpp[i] += ata[i][j]*atx[j];
                    }
                }
                da0  = 0;
                chi2 = 0;
                if (bZero)
                {
                    for (j = 0; j < nusemol; j++)
                    {
                        ax = a0;
                        for (size_t i = 0; i < ptypes.size(); i++)
                        {
                            ax += fpp[i]*a_copy[j][i];
                        }
                        da0  += (x_copy[j]-ax);
                        chi2 += gmx::square(x_copy[j]-ax);
                    }
                    da0 /= nusemol;
                    a0  += da0;
                    niter++;
                    printf("iter: %d <pol> = %g, a0 = %g, chi2 = %g\n",
                           niter, poltot/nusemol, a0, chi2/nusemol);
                }
            }
            while (bZero && (fabs(da0) > 1e-5) && (niter < 1000));           
            for (size_t i = 0; i < ptypes.size(); i++)
            {
                ptypes[i].addPolValue(fpp[i]);
                gmx_stats_add_point_ydy(ptypes[i].fullStats(), fpp[i], 0);
            }
        }
        free_matrix(a_copy);
        free_matrix(at_copy);
        free_matrix(ata);
    }
    fprintf(stderr, "\n");

    //Make the histogram 
    FILE *xp = xvgropen(hisfn, "Polarizability distribution", "alpha (A\\S3\\N)", "", oenv);          
    std::vector<std::string> leg;
    for (auto p = ptypes.begin(); p < ptypes.end(); ++p)
    {
        size_t      pos   = p->name().find("p_");
        std::string ptype = p->name();
        if (pos != std::string::npos)
        {
            ptype = ptype.substr(pos+2);
        }
        leg.push_back(ptype);
    }
    std::vector<const char*> legend;
    for (auto l = leg.begin(); l < leg.end(); ++l)
    {
        legend.push_back(l->c_str());
    }
    xvgr_legend(xp, legend.size(), legend.data(), oenv);
    for (auto ptype = ptypes.begin(); ptype < ptypes.end(); ++ptype)
    {
        real aver, sigma;
        ptype->improveStats();
        if ((estatsOK == gmx_stats_get_average(ptype->improvedStats(), &aver)) &&
            (estatsOK == gmx_stats_get_sigma(ptype->improvedStats(), &sigma)))
        {
            int   nbins = 1+sqrt(nAccepted);
            real *x, *y;           
            pd.setPtypePolarizability(ptype->name().c_str(), aver, sigma);
            fprintf(fplog, "%-5s  %8.3f +/- %.3f\n", ptype->name().c_str(), aver, sigma);           
            if (estatsOK == gmx_stats_make_histogram(ptype->fullStats(), 0, &nbins, ehistoY, 1, &x, &y))
            {
                fprintf(xp, "@type xy\n");
                for (int ll = 0; ll < nbins; ll++)
                {
                    fprintf(xp, "%.3f  %.3f\n", x[ll], y[ll]);
                }
                fprintf(xp, "&\n");
            }
            else
            {
                fprintf(stderr, "Could not make histogram for %s\n", ptype->name().c_str());
            }
            free(x);
            free(y);
        }
        else
        {
            fprintf(stderr, "Could not determine polarizability for %s\n", ptype->name().c_str());
        }
        gmx_stats_free(ptype->fullStats());
        gmx_stats_free(ptype->improvedStats());
    }
    fclose(xp);
    if (bZero)
    {
        const char *null = (const char *)"0";
        pd.addPtype(null, nullptr, null, a0, 0);
    }
    free_matrix(a);
    free_matrix(at);
    return ptypes.size();
}

int alex_tune_pol(int argc, char *argv[])
{
    static const char                     *desc[] =
    {
        "tune_pol optimizes atomic polarizabilities that together build",
        "an additive model for polarization. The set of atomtypes used",
        "is determined by the input force field file ([TT]-di[tt] option). All",
        "atomtypes for which the polarizability is zero, and for which",
        "there is sufficient data in the input data file ([TT]-f[tt] option)",
        "will be used in the least-squares fit (done by Singular Value Decomposition).[PAR]"
        "Bootstrapping is used to estimate the error in the resulting parameters",
        "and this will be stored in the resulting file such that programs",
        "using the parameters can estimate an error in the polarizability.[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored.[PAR]",
        "tune_pol produces a table and a figure to describe the data."
    };
    t_filenm                               fnm[] =
    {
        { efDAT, "-f",     "data",      ffRDMULT },
        { efDAT, "-o",     "allmols",   ffOPTWR  },
        { efDAT, "-di",    "gentop",    ffOPTRD  },
        { efDAT, "-do",    "tune_pol",  ffWRITE  },
        { efDAT, "-sel",   "molselect", ffREAD   },
        { efLOG, "-g",     "tune_pol",  ffWRITE  },
        { efTEX, "-atype", "atomtypes", ffWRITE  },
        { efXVG, "-his",   "polhisto",  ffWRITE  }
    };
    int                                    NFILE       = asize(fnm);
    static char                           *sort[]      = {nullptr, (char *)"molname", (char *)"formula", (char *)"composition", nullptr};
    static char                           *exp_type    = (char *)"Polarizability";
    static char                           *zeropol     = (char *) nullptr;
    static gmx_bool                        bQM         = false;
    static int                             mindata     = 1;
    static int                             nBootStrap  = 1000;
    static int                             maxwarn     = 0;
    static int                             seed;
    static char                           *lot         = (char *)"B3LYP/aug-cc-pVTZ";
    static real                            sigma       = 0;
    static gmx_bool                        bZero       = false;
    static gmx_bool                        bForceFit   = false;
    static gmx_bool                        bCompress   = true;
    t_pargs                                pa[]        =
    {
        { "-sort",   FALSE, etENUM, {sort},
          "Key to sort the final data file on." },
        { "-maxwarn", FALSE, etINT, {&maxwarn},
          "Will only write output if number of warnings is at most this." },
        { "-qm",     FALSE, etBOOL, {&bQM},
          "Use QM data for optimizing the empirical polarizabilities as well." },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory" },
        { "-mindata", FALSE, etINT, {&mindata},
          "Minimum number of data points to optimize a polarizability value" },
        { "-sigma",  FALSE, etREAL, {&sigma},
          "Assumed standard deviation (relative) in the experimental data points" },
        { "-zero",    FALSE, etBOOL, {&bZero},
          "Use a zero term (like in Bosque and Sales)" },
        { "-zeropol", FALSE, etSTR,  {&zeropol},
          "Polarizability for these poltypes (polarizability types) is fixed to zero. Note that all poltypes in gentop.dat that are non zero are fixed as well." },
        { "-force",   FALSE, etBOOL, {&bForceFit},
          "Reset all polarizablities to zero in order to re-optimize based on a gentop.dat file with previous values" },
        { "-nBootStrap", FALSE, etINT, {&nBootStrap},
          "Number of trials for bootstrapping" },
        { "-seed", FALSE, etINT, {&seed},
          "Seed for random numbers in bootstrapping. If <= 0 a seed will be generated." },
        { "-exp_type", FALSE, etSTR, {&exp_type},
          "[HIDDEN]The experimental property type that all the stuff above has to be compared to." },
        { "-compress", FALSE, etBOOL, {&bCompress},
          "Compress output XML files" }
    };
    
    FILE                                  *fplog;
    FILE                                  *tp;
    int                                    i, nalexandria_atypes;
    int                                    nwarn = 0;
    std::vector<alexandria::MolProp>       mp;
    std::vector<std::string>               zpol;
    MolPropSortAlgorithm                   mpsa;
    alexandria::MolSelect                  gms;
    gmx_atomprop_t                         ap;
    alexandria::Poldata                    pd;
    gmx_output_env_t *                     oenv;
    const char                            *pdout;
    const char                            *atype;
    const char                            *mpfn;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    if (nBootStrap <= 1)
    {
        gmx_fatal(FARGS, "nBootStrap should be >= 1");
    }
    
    ap = gmx_atomprop_init();
    
    try
    {
        alexandria::readPoldata(opt2fn_null("-di", NFILE, fnm), pd, ap);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    auto fns = opt2fns("-f", NFILE, fnm);
    nwarn = merge_xml(fns, &mp, nullptr, nullptr, (char *)"double_dip.dat", ap, pd, true);
    
    if (nwarn > maxwarn)
    {
        gmx_fatal(FARGS, "Too many warnings (%d). Terminating.\n", nwarn);
    }   
    if (nullptr != debug)
    {
        for (auto mpi = mp.begin(); mpi < mp.end(); mpi++)
        {
            fprintf(debug, "%s", mpi->getMolname().c_str());
            for (auto mci = mpi->BeginMolecularComposition(); mci < mpi->EndMolecularComposition(); ++mci)
            {
                fprintf(debug, "  %s", mci->getCompName().c_str());
            }
            fprintf(debug, "\n");
        }
    }
    
    gms.read(opt2fn("-sel", NFILE, fnm));
    
    fplog  = opt2FILE("-g", NFILE, fnm, "w");   
    if (nullptr != zeropol)
    {
        zpol = gmx::splitString(zeropol);
    }

    nalexandria_atypes = decompose_frag(fplog, 
                                        opt2fn("-his", NFILE, fnm),
                                        pd, 
                                        mp, 
                                        bQM, 
                                        lot, 
                                        mindata,
                                        gms, 
                                        bZero, 
                                        bForceFit,
                                        nBootStrap,
                                        zpol, 
                                        oenv, 
                                        exp_type);
                                        
    fprintf(fplog, "There are %d alexandria atom types\n", nalexandria_atypes);

    pdout = opt2fn("-do", NFILE, fnm);
    fprintf(fplog, "Now writing force field file %s\n", pdout);
    alexandria::writePoldata(pdout, &pd,  bCompress);

    atype = opt2fn("-atype", NFILE, fnm);
    fprintf(fplog, "Now writing LaTeX description of force field to %s\n", atype);
    tp = fopen(atype, "w");
    std::vector<alexandria::Poldata> pds = { pd };
    alexandria_molprop_atomtype_table(tp, true, pds, mp, lot, exp_type);
    fclose(tp);
    gmx_ffclose(fplog);

    mpfn = opt2fn_null("-o", NFILE, fnm);
    if (nullptr != mpfn)
    {
        mpsa = MPSA_NR;
        if (opt2parg_bSet("-sort", asize(pa), pa))
        {
            for (i = 0; (i < MPSA_NR); i++)
            {
                if (strcasecmp(sort[0], sort[i+1]) == 0)
                {
                    mpsa = (MolPropSortAlgorithm) i;
                    break;
                }
            }
        }
        if (mpsa != MPSA_NR)
        {
            alexandria::MolSelect gms;
            MolPropSort(&mp, mpsa, ap, gms);
        }

        MolPropWrite(mpfn, &mp, bCompress);
    }
    return 0;
}
