/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#ifndef MOLGEN_H
#define MOLGEN_H

#include <map>

#include "gromacs/commandline/pargs.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/real.h"

#include "act/utility/communicationrecord.h"
#include "molselect.h"
#include "actmol.h"

namespace alexandria
{

//! Different components of the chi-squared used in fitting
enum class eRMS { 
    //! Parameters going out of bound
    BOUNDS,
    //! Unphysical parameters
    UNPHYSICAL,
    //! Charge
    CHARGE,
    //! Molecular dipole
    MU,
    //! Molecular quadrupole
    QUAD,
    //! Molecular octupole
    OCT,
    //! Molecular hexadecapole
    HEXADEC,
    //! Vibrational frequencies
    FREQUENCY,
    //! Infrared intensities
    INTENSITY,
    //! Deviation from CM5 charge
    CM5,
    //! Deviation from Electrostatic potential
    ESP,
    //! Potential energy deviation
    EPOT,
    //! Interaction energy deviation
    Interaction,
    //! Mean square force on the atoms
    Force2,
    //! Deviation of polarizability components
    Polar,
    //! Total root mean square deviation
    TOT
};

//! \return map from each eRMS to its string name
const std::map<eRMS, const char *> &geteRMSNames();

//! \brief Return string corresponding to eRMS
const char *rmsName(eRMS e);

/*! \brief Class to manage one of the fitting targets
 */
class FittingTarget 
{

private:

    //! The actual term stored here
    eRMS       erms_;
    //! The data set used
    iMolSelect ims_;
    //! The weighting factor in the total chi squared
    real       weight_             = 0;
    //! The number of data points for this property
    int        numberOfDatapoints_ = 0;
    //! The unweighted chi squared for this property
    real       chiSquared_         = 0;

public:

    /*! \brief Constructor
     * \param[in] e   The fitting target
     * \param[in] ims The data set used
     */
    FittingTarget(eRMS e, iMolSelect ims) : erms_(e), ims_(ims) {};

    //! \return the eRMS component covered by this FittingTarget
    eRMS erms() const { return erms_; }

    /*! \brief Set the weight factor
     * \param[in] w The new weighting factor
     */
    void setWeight(real w) { weight_ = w; }
    
    //! \return the weighting factor
    real weight() const { return weight_; }
    
    //! \return a pointer to the address of the weighting factor
    real *weightPtr() { return &weight_; }
    
    //! \return the data set used
    iMolSelect ims() const { return ims_; }
    
    /*! \brief Set number of data points
     * \param[in] n The new number of data points
     */
    void setNumberOfDatapoints(int n) { numberOfDatapoints_ = n; }
    
    //! \return the number of data points
    int numberOfDatapoints() const { return numberOfDatapoints_; }
    
    /*! \brief Set new value for chi squared
     * \param[in] chi2 The new chi squared
     */
    void setChiSquared(real chi2) { chiSquared_ = chi2; }
    
    //! \return the present chi squared
    real chiSquared() const { return chiSquared_; }
    
    //! \return the weighted chi squared
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
    
    //! \brief Rest the number of data points and chi squared to zero.
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

class AtomIndex
{
private:
    //! Atom type name
    std::string  name_;
    //! Atom type count
    int          count_;
    //! Whether or not this atom type can be changed
    bool         const_;
public:
    /*! \brief Constructor
     * \param[in] name   The atom type name
     * \param[in] bConst Whether or not this can be changed
     */
    AtomIndex(const std::string &name, bool bConst) :
        name_(name), count_(0), const_(bConst) {}

    //! \return the atom type name
    const std::string name() const { return name_; }

    //! \return the number of copies
    int count() const { return count_; }

    //! \return whether or not this is constant
    bool isConst() const { return const_; }
        
    //! Increment the count
    void increment() { count_++; }

    //! Reduced count if larger than zero
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

//! Map from RMS type to FittingTarget structure
// using RmsFittingTarget       = typename std::map<eRMS, FittingTarget>;

/* \brief Class to manage the molecules in a force field optimization
 */
class MolGen
{

private:

    //! Communication record
    const CommunicationRecord      *cr_;
    //! Minimum number of data points to consider a parameter
    int                             mindata_    = 1;
    //! Map that holds the number of compounds in each data set
    std::map<iMolSelect, size_t>    targetSize_;
    //! Tell us whether this interaction type needs optimizing
    std::map<InteractionType, bool> iOpt_;
    //! Whether or not to use charge symmetry
    bool                            qsymm_      = false;
    //! GROMACS logger structure
    gmx::MDLogger                   mdlog_;
    //! GROMACS MD parameter structure
    t_inputrec                     *inputrec_   = nullptr;
    //! String for command line to harvest the options to fit
    char                           *fitString_ = nullptr;
    //! Map to determine whether or not to fit a parameter type
    std::map<std::string, bool>     fit_;
    //! GROMACS structure containing optional MD modules, used for electric fields
    gmx::MDModules                  mdModules_;
    //! Which charge type to use
    char                           *chargeMethod_ = nullptr;
    //! The molecules used in the optimization
    std::vector<alexandria::ACTMol>  actmol_;
    //! Whether or not to load balance the inputs
    bool                            loadBalance_ = true;
    
    /*! \brief Check that we have enough data 
     * Check that we have enough data for all parameters to optimize
     * in this molecule.
     * \param[in] fp File to print logging information to. May be nullptr.
     * \param[in] pd Pointer to forcefield object
     */
    void checkDataSufficiency(FILE     *fp,
                              ForceField  *pd);
    
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
    
    
public:
    
    //! Default constructor
    MolGen() {}

    /*! \brief 
     * Constructor of MolGen class.
     * \param[in] cr communication record
     */ 
    MolGen(const CommunicationRecord *cr);
    
    /*! \brief Add options to the command line
     * \param[in] pargs   Vector of command line arguments
     * \param[in] targets pointer to targets map for train dataset. Comes from sharedIndividualInfo
     */
    void addOptions(std::vector<t_pargs> *pargs, std::map<eRMS, FittingTarget> *targets);
    
    /*! \brief Process options after parsing
     */
    void optionsFinished();

    /*! \brief Manually add fitting option
     * \param[in] opt The string, e.g. sigma or alpha
     */
    void addFitOption(const std::string &opt) { fit_.insert({ opt, true }); }

    /*! \brief Fill the  iOpt_ map
     * \param[in] pd Pointer to forcefield
     */
    void fillIopt(ForceField *pd);
    
    //! \brief Return the const vector of molecules
    const std::vector<ACTMol> &actmols() const { return actmol_; }

    //! \brief Return the mutable vector of molecules
    std::vector<ACTMol> &actmols() { return actmol_; }

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
    
    //! Return the mdlogger structure
    const gmx::MDLogger &mdlog()  const {return mdlog_; }

    /*! Tell the user whether this parameter needs to be fitted
     * \param[in] type The parameter type
     * \return whether or not this parameter is part of the fitting target
     */
    bool fit(const std::string &type) const
    { 
        return fit_.find(type) != fit_.end(); 
    }
    
    //! Return all the parameter types to fit
    std::map<std::string, bool> typesToFit() const { return fit_; }
        
    //! Tell the user whether this interaction type needs optimization
    bool optimize(InteractionType itype) const
    {
        return iOpt_.find(itype) !=  iOpt_.end(); 
    }
    
    //! Return which chargemethod was selectod
    const char *chargeMethod() const { return chargeMethod_; }
        
    //! Return the whole optimization map
    const std::map<InteractionType, bool> iopt() const { return iOpt_; }
    
    //! Whether or not to use charge symmetry
    gmx_bool bQsym() const { return qsymm_;}
        
    //! \brief Return minimum amount of data needed
    int mindata() const { return mindata_; }
  
    /*! \brief Read the molecular property data file to generate molecules.
     * \param[in] fp      File pointer for printing information, may be nullptr
     * \param[in] fn      Filename for to read molecules from
     * \param[in] pd      Pointer to ForceField object
     * \param[in] gms     The molecule selection
     * \param[in] verbose Whether or not to print extra information
     * \return number of molecules read and processed correctly
     */
    size_t Read(FILE                                *fp,
                const char                          *fn,
                ForceField                             *pd,
                const MolSelect                     &gms,
                const std::map<eRMS, FittingTarget> &targets,
                bool                                 verbose);

};

} //namespace alexandria

#endif //MOLGEN_H
