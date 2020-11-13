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

#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/real.h"

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

enum eRMS {
    ermsBOUNDS  = 0,
    ermsMU      = 1,
    ermsQUAD    = 2,
    ermsCHARGE  = 3,
    ermsESP     = 4,
    ermsEPOT    = 5,
    ermsForce2  = 6,
    ermsPolar   = 7,
    ermsTOT     = 9,
    ermsNR      = 10
};

//! \brief Return string corresponding to eRMS
const char *rmsName(int e);

namespace alexandria
{

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

/*! \brief Convenience storage of parameters to optimize
 */
class OptimizationIndex
{
 private:
    //! The interaction type
    InteractionType iType_;
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
        iType_(iType),parameterId_(parameterId), parameterType_(parameterType) {}
    
    //! Return the interaction type
    InteractionType iType() const { return iType_; }
    
    //! Return the id
    Identifier id() const { return parameterId_; }
    
    //! Return the type
    const std::string type() const { return parameterType_; }

    //! Return a compound string representing the index
    std::string name() const
    {
        return gmx::formatString("%s-%s", //interactionTypeToString(iType_).c_str(),
                                 parameterId_.id().c_str(), parameterType_.c_str()); 
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
        int                             nmol_support_;
        int                             mindata_;
        int                             maxESP_;
        int                             nexcl_;
        int                             nexcl_orig_;
        real                            watoms_;
        eTune                           etune_;
        real                            qtol_;
        int                             qcycle_;
        real                            ener_[ermsNR] = { 0 };
        real                            fc_[ermsNR] = { 0 };
        gmx_bool                        bQM_;
        gmx_bool                        bDone_;
        gmx_bool                        bFinal_;
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

        //! \brief Process options after parsing
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
        
        //! \brief Is this the last calculation?
        bool final() const { return bFinal_; }

        //! \brief Are we using QM only?
        bool bQM() const { return bQM_; }

        //! \brief Return level of theory
        const char *lot() const { return lot_; }

        void setFinal() { bFinal_ = true; }
        
        //! \brief Return the number of compounds in the data set
        int nMolSupport() const { return nmol_support_; }

        //! \brief Set the energy value for the corresponding rms.
        void setEnergy(int rms, real ener)
        {
            ener_[rms] = ener;
        }
        
        //! \brief Return the energy of the corresponding rms. 
        double energy(int rms) const 
        { 
            return ener_[rms]; 
        }
        
        //! \brief Return the weighting factor of the energy.
        double weight(int rms) const 
        { 
            return fc_[rms]; 
        }
        
        //! \brief Set all the energies to zero.
        void resetEnergies()
        {
            for (int rms = 0; rms < ermsNR; rms++)
            {
                ener_[rms] = 0;
            }
        }

        //! \brief Increase the rms component of the energy vector by delta. 
        void increaseEnergy(int rms, real delta)
        {
            ener_[rms] += delta;
        }

        /*! \brief Sum over the energies of the cores.
         * Also multiplies the terms by the weighting factors.
         */
        void sumEnergies();

        /*! \brief Print the energy components.
         * \param[in] fp File pointer to print to, maybe nullptr
         */  
        void printEnergies(FILE *fp) const;

        /*! \brief Read the molecular property data file to generate molecues.
         * TODO: update comments
         */
        void Read(FILE            *fp,
                  const char      *fn,
                  const char      *pd_fn,
                  gmx_bool         bZero,
                  const MolSelect &gms,
                  bool             bZPE,
                  bool             bDHform,
                  const char      *tabfn,
                  iMolSelect       SelectType);

};

}

#endif
