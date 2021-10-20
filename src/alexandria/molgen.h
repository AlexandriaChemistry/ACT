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
 
#ifndef MOLGEN_H
#define MOLGEN_H

#include <map>

#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/real.h"

#include "molselect.h"
#include "mymol.h"

typedef struct {
    int       n;
    int       nopt;
    int       nconst;
    int       nopt_c;
    int      *tot_count;
    int      *count;
    char    **name;
    gmx_bool *bConst;
} t_index_count;

extern char *opt_index_count(t_index_count *ic);

namespace alexandria
{

enum class eRMS { BOUNDS, MU, QUAD, CHARGE, CM5, ESP, EPOT, Force2, Polar, TOT };

//! \brief Return string corresponding to eRMS
const char *rmsName(eRMS e);

class FittingTarget 
{
private:
    eRMS       erms_; 
    iMolSelect ims_;
    real       weight_             = 0;
    int        numberOfDatapoints_ = 0;
    real       chiSquared_         = 0;
public:
    FittingTarget(eRMS e, iMolSelect ims) : erms_(e), ims_(ims) {};
  
    void setWeight(real w) { weight_ = w; }
    
    real weight() const { return weight_; }
    
    real *weightPtr() { return &weight_; }
    
    iMolSelect ims() const { return ims_; }
    
    void setNumberOfDatapoints(int n) { numberOfDatapoints_ = n; }
    
    int numberOfDatapoints() const { return numberOfDatapoints_; }
    
    void setChiSquared(real chi2) { chiSquared_ = chi2; }
    
    real chiSquared() const { return chiSquared_; }
    
    real chiSquaredWeighted() const 
    { 
        if (chiSquared_ > 0 && numberOfDatapoints_ > 0)
        { 
            return chiSquared_*weight_/numberOfDatapoints_;
        }
        else
        {
            return 0;
        }
    } 

    /*! \brief
     * Increase the chi2 by delta. 
     * \param[in] ndata Number of data points
     * \param[in] delta Change in chiSquared
     */
    void increase(int ndata, real delta)
    {
        numberOfDatapoints_ += ndata;
        chiSquared_         += delta;
    }
    
    void reset()
    {
        numberOfDatapoints_ = 0;
        chiSquared_         = 0;
    }
  
    /*! \brief Print if non zero
     * \param[in] fp   File pointer, print only if non null
     */
    void print(FILE *fp) const;
};
  
enum eTune {
    etuneEEM, etuneFC, etuneNone
};

class AtomIndex
{
    private:
        std::string  name_;
        int          count_;
        bool         const_;
    public:
        AtomIndex(const std::string name, bool bConst) :
            name_(name), count_(0), const_(bConst) {}

        const std::string name() const { return name_; }

        int count() const { return count_; }

        bool isConst() const { return const_; }
        
        void increment() { count_++; }

        void decrement()
        {
            if (count_ > 0)
            {
                count_--;
            }
            else
            {
                fprintf(stderr, "Trying to decrease number of atoms %s below zero\n",
                        name().c_str());
            }
        }
};

using AtomIndexIterator      = typename std::vector<AtomIndex>::iterator;
using AtomIndexConstIterator = typename std::vector<AtomIndex>::const_iterator;
using RmsFittingTarget       = typename std::map<eRMS, FittingTarget>;

/*! \brief Convenience storage of parameters to optimize
 */
class OptimizationIndex
{
 private:
    //! The particle type. If empty, an interaction type is used.
    std::string     particleType_;
    //! The interaction type
    InteractionType iType_ = InteractionType::CHARGE;
    //! The name of the parameter matching poldata
    Identifier      parameterId_;
    //! The type of parameter, eg. sigma, epsilon
    std::string     parameterType_;
 public:
    /*! \brief Constructor
     * \param[in] iType         The interaction type
     * \param[in] parameterId   The identifier
     * \param[in] parameterType The type
     */
    OptimizationIndex(InteractionType    iType,
                      Identifier         parameterId,
                      const std::string &parameterType) :
        iType_(iType), parameterId_(parameterId), parameterType_(parameterType) {}
    
    /*! \brief Constructor
     * \param[in] pType         The particle type
     * \param[in] parameterType The type
     */
    OptimizationIndex(const std::string  pType,
                      const std::string &parameterType) :
        particleType_(pType), parameterType_(parameterType) {}
    
    //! Return the interaction type
    InteractionType iType() const { return iType_; }
    
    //! Return the id
    Identifier id() const { return parameterId_; }
    
    //! Return particle type
    const std::string &particleType() const { return particleType_; }
    
    //! Return the type
    const std::string &parameterType() const { return parameterType_; }

    //! Return a compound string representing the index
    std::string name() const
    {
        if (InteractionType::CHARGE == iType_)
        {
            return gmx::formatString("%s-%s",
                                     particleType_.c_str(),
                                     parameterType_.c_str());
        }
        else
        {
            return gmx::formatString("%s-%s",
                                     parameterId_.id().c_str(),
                                     parameterType_.c_str());
        }
    }
    
    //! Return a compound string representing the atomindex
    std::string parameterName() const
    {
        return parameterId_.id().c_str();
    }
};
    
class MolGen
{
    private:
        int                             mindata_;
        int                             maxESP_;
        int                             nexcl_;
        int                             nexcl_orig_;
        real                            watoms_;
        eTune                           etune_;
        real                            qtol_;
        int                             qcycle_;
        std::map<iMolSelect, RmsFittingTarget> targets_;
        std::map<iMolSelect, size_t>    targetSize_;
        gmx_bool                        bQM_;
        gmx_bool                        bDone_;
        gmx_bool                        bGenVsite_;
        gmx_bool                        qsymm_;
        gmx_bool                        constrain_; 
        Poldata                         pd_;
        t_commrec                      *cr_;
        gmx::MDLogger                   mdlog_;
        t_inputrec                     *inputrec_;
        gmx_hw_info_t                  *hwinfo_;
        gmx_atomprop_t                  atomprop_;
        //! String for command line to harvest the options to fit
        char                            *fitString_ = nullptr;
        //! Map to determine whether or not to  fit a parameter type
        std::map<std::string, bool>     fit_;
        gmx::MDModules                  mdModules_;
        std::vector<alexandria::MyMol>  mymol_;
        const char                     *lot_;
        /*! \brief Check that we have enough data 
         * Check that we have enough data for all parameters to optimize
         * in this molecule.
         * \param[in] fp File to print logging information to. May be nullptr.
         */
        void checkDataSufficiency(FILE *fp);
        /*! \brief Generate optIndex
         * \param[in] fp File to print logging information to. May be nullptr.
         */
        void generateOptimizationIndex(FILE *fp);

        //! Compute amount of compounds in each group
        void countTargetSize();
        
        size_t nTargetSize(iMolSelect ims) const
        {
            auto ts = targetSize_.find(ims);
            if (ts != targetSize_.end())
            {
                return ts->second;
            }
            return 0;
        }

        //! \brief Fill the  iOpt_ map
        void fillIopt();
        //! Tell us whether this interaction type needs optimizing
        std::map<InteractionType, bool> iOpt_;

    public:
        std::vector<OptimizationIndex>  optIndex_;
        /*! \brief 
         * Constructor of MolGen class.
         */ 
        MolGen();

        /*! \brief 
         * Deconstructor of MolGen class.
         */
        ~MolGen();

        //! \brief Add options to the command line
        void addOptions(std::vector<t_pargs> *pargs, eTune etune);

        /*! \brief Process options after parsing
         */
        void optionsFinished();

        //! \brief Return the poldata as const variable
        const Poldata *poldata() const { return &pd_; }
        
        //! \brief Return the poldata
        Poldata *poldata() { return &pd_; }

        //! \brief Return the atomprop structure
        gmx_atomprop_t atomprop() const { return atomprop_; }

        //! \brief Return the const vector of molecules
        const std::vector<MyMol> &mymols() const { return mymol_; }

        //! \brief Return the mutable vector of molecules
        std::vector<MyMol> &mymols() { return mymol_; }

     	/*! \brief Return size of data set
         * \param[in] ims The data set 
         * \return The size of it
         */
        size_t iMolSelectSize(iMolSelect ims) const
        {
            auto ts = targetSize_.find(ims);
            if (targetSize_.end() == ts)
            {
                return 0;
            }
            else
            {
                return ts->second; 
            }
        }

        gmx::MDLogger  mdlog()  const {return mdlog_; }

        //! Tell the user whether this parameter needs to be fitted
        bool fit(const std::string &type) const { return fit_.find(type) != fit_.end(); }
        //! Return all the type to  fit
        std::map<std::string, bool> typesToFit() const { return fit_; }
        
        //! Tell  the user whether this interaction type needs optimization
        bool optimize(InteractionType itype) const { return iOpt_.find(itype) !=  iOpt_.end(); }
        
        //! Return the whole optimization map
        const std::map<InteractionType, bool> iopt() const { return iOpt_; }
        
        int qcycle() const { return qcycle_;}
        
        int maxPot() const { return maxESP_;}
        
        double qtol() const { return qtol_;}
        
        gmx_bool bQsym() const { return qsymm_;}
        
        gmx_bool bConstrain() const { return constrain_;}
        
        //! \brief Return ESP weighting factor for atoms
        real watoms() const { return watoms_; }

        //! \brief Return minimum amount of data needed
        int mindata() const { return mindata_; }

        //! \brief Return const communication record
        const t_commrec *commrec() const { return cr_; }

        //! \brief Return non-const communication record
        t_commrec *commrec() { return cr_; }
        
        gmx_hw_info_t *hwinfo() {return hwinfo_;}
        
        //! \brief Are we using QM only?
        bool bQM() const { return bQM_; }

        //! \brief Return level of theory
        const char *lot() const { return lot_; }

        /*! \brief Return the fitting targets for editing
         * \param[in] ims The selection to return
         * \return The map of fittingtargets or nullptr
         */
        std::map<eRMS, FittingTarget> *fittingTargets(iMolSelect ims)
        {
            auto tt = targets_.find(ims);
            if (targets_.end() == tt)
            {
                return nullptr;
            }
            else
            {
                return &tt->second;
            }
        }
        
        const std::map<eRMS, FittingTarget> &fittingTargetsConst(iMolSelect ims) const
        {
            auto tt = targets_.find(ims);
            GMX_RELEASE_ASSERT(targets_.end() != tt, gmx::formatString("Cannot find selection %s", iMolSelectName(ims)).c_str());
            return tt->second;
        }
        
        /*! \brief return appropriate fitting target
         * \param[in] ims The selection
         * \param[in] rms The contributor to the chi squared
         * \return FittingTarget class or nullptr if not found
         */
        FittingTarget *target(iMolSelect ims, eRMS rms);
	
        /*! \brief Set the chiSquared to zero.
         * \param[in] ims The selection to reset
         */
        void resetChiSquared(iMolSelect ims)
        {
            auto fts = fittingTargets(ims);
            if (fts != nullptr)
            {
                for (auto &ft : *fts)
                {
                    ft.second.reset();
                }
            }
        }

        /*! \brief 
         * Sum over the energies of the cores if desired.
         * Also multiplies the terms by the weighting factors.
         * \param[in] parallel Whether or not to sum in parallel
         * \param[in] ims      The selection to sum
         */
        void sumChiSquared(bool parallel, iMolSelect ims);

        /*! \brief Print the chiSquared components.
         * \param[in] fp  File pointer to print to, may be nullptr
         * \param[in] ims The selection to print
         */  
        void printChiSquared(FILE *fp, iMolSelect ims) const;

        /*! \brief Read the molecular property data file to generate molecules.
         * \param[in] fp      File pointer for printing information
         * \param[in] fn      Filename for molecules
         * \param[in] pd_fn   Filename for force field file
         * \param[in] bZero   Use compounds with zero dipole
         * \param[in] gms     The molecule selection
         * \param[in] bZPE    Use Zero point energy
         * \param[in] bDHform Use delta H formation
         * \param[in] tabfn   Table function for gromacs
         * \param[in] verbose Whether or not to print stuff
         * \return number of molecules read and processed correctly
         */
        size_t Read(FILE            *fp,
                  const char      *fn,
                  const char      *pd_fn,
                  gmx_bool         bZero,
                  const MolSelect &gms,
                  bool             bZPE,
                  bool             bDHform,
                  const char      *tabfn,
                  bool             verbose);

};

}

#endif
