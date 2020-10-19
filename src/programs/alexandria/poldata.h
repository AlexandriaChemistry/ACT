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
 
#ifndef POLDATA_H
#define POLDATA_H

#include <algorithm>
#include <vector>

#include "gromacs/utility/stringutil.h"

#include "chargemodel.h"
#include "forcefieldparameter.h"
#include "forcefieldparameterlist.h"
#include "interactiontype.h"
#include "poldata_low.h"
#include "stringutil.h"

/* This source code file is part of the Alexandria project */

struct t_commrec;

namespace alexandria
{

class Poldata
{
    public:

        Poldata() {};

        /*! \brief
         * Set the file name gentop.dat
         *
         */
        void  setFilename(const std::string &fn2);

        /*! \brief
         * Set the force field version
         * \param[in] version Force field version
         */
        void setVersion(const std::string &version)
        {
            alexandriaVersion_ = version;
        }

        /*! \brief
         * Get the force field version
         * \return Force field version
         */
        const std::string &getVersion() const
        {
            return alexandriaVersion_;
        }

        //! Return whether or not the original Rappe & Goddard is used
        bool rappe() const;
        
        //! Return whether Yang & Sharp is used
        bool yang() const;
        
        //! Return whether we use a point charge for cores
        bool corePointCharge() const { return false; }

        /*! \brief
         * Set relative dielectric constant
         *
         */
        void setEpsilonR(double epsilonR) { gtEpsilonR_ = epsilonR; }

        /*! \brief
         * Set the number of exclusion
         */
        void setNexcl(int nexcl) { nexcl_ = nexcl; }

        /*! \brief
         * Add the atom types used in Alexandria FF
         *
         **\param[in] elem          The element of the atom type
         **\param[in] desc          The description of the atom type
         **\param[in] atype         The atom type defined for this element in Alexandria FF
         **\param[in] ptype         The polarizability type defined for this element in Alexandria FF
         **\param[in] btype         The bond type defined for elem in Alexandria FF
         **\param[in] ztype         The zeta type defined for elem in Alexandria FF
         **\param[in] ref_enthalpy  The reference enthalpy for elem
         */
        void  addAtype(const std::string &elem,
                       const std::string &desc,
                       const std::string &atype,
                       const std::string &ptype,
                       const std::string &btype,
                       const std::string &ztype,
                       const std::string &ref_enthalpy);

        /*! \brief Add an atom type
         * \param[in] atp The new atom type
         */
        void  addAtype(Ffatype atp) { alexandria_.push_back(atp); }
        /*! \brief
         *  Add the polarizability types
         *
         **\param[in] ptype           The name specifying the polarizability type in Alexandria FF
         **\param[in] miller          Miller polarizability type
         **\param[in] bosque          Bosque polarizability type
         **\param[in] polarizability  The calulated value of polarizability
         **\param[in] sigPol          The uncertainty of the calculated polarizability
         */
        void  addPtype(const std::string &ptype,
                       const std::string &miller,
                       const std::string &bosque,
                       double             polarizability,
                       double             sigPol);

        void  addVsite(const std::string &atype,
                       const std::string &type,
                       int                number,
                       double             distance,
                       double             angle,
                       int                ncontrolatoms);

        /*! \brief
         * Set the vsite angle unit.
         */
        void setVsite_angle_unit(const std::string &angle_unit)
        {
            vsite_angle_unit_ = angle_unit;
        }

        /*! \brief
         * Set the vsite angle unit.
         */
        void setVsite_length_unit(const std::string &length_unit)
        {
            vsite_length_unit_ = length_unit;
        }

        std::vector<Vsite> &getVsite() {return vsite_; }

        int getNexcl() const { return nexcl_; }

        size_t getNatypes() const { return alexandria_.size(); }

        /*! \brief
         * Return the reference enthalpy for the given atom type
         *
         **\param[in] atype Atom type
         **\param[ou] Href  Reference enthalpy
         */
        bool getAtypeRefEnthalpy(const std::string &atype,
                                 double            *Href) const;

        /*! \brief
         * Return the relative dielectric constant
         */
        double getEpsilonR() const { return gtEpsilonR_; }

        std::string  getGeometry(std::string gtBrule);

        /*! \brief
         * Return the discription corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         */
        const std::string &getDesc(const std::string &atype) const;

        /*! \brief
         * Return the element name corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         */
        const std::string &getElem(const std::string &atype) const;
        
        /*! \brief
         * Return the element name corresponding to the zeta type
         *
         * \param[in] ztype  zeta Type
         */
        const std::string &ztype2elem(const std::string &ztype) const;

        /*! \brief
         * Return the atom type corresponding to the zeta type
         *
         * \param[in] ztype  zeta Type
         */
        //const std::string &ztype2atype(const std::string &ztype) const;

        //std::vector<std::string> ztype_names() const;

        /*! \brief
         * Return the charge corresponding to the atyom type
         * from the gentop.dat file
         *
         * \param[in] atype  Atom type
         */
        std::string  getCharge(  std::string atype);
        
        std::vector<Ffatype> &getAtypes() {return alexandria_; }

        FfatypeIterator getAtypeBegin() { return alexandria_.begin(); }

        FfatypeIterator getAtypeEnd() { return alexandria_.end(); }

        FfatypeConstIterator getAtypeBegin() const { return alexandria_.begin(); }

        FfatypeConstIterator getAtypeEnd() const { return alexandria_.end(); }

        /*! \brief
         * Return the iterator corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         * \throw if atom type not  found.
         */
        FfatypeIterator findAtype(const std::string &atype)
        {
            auto atp = std::find_if(alexandria_.begin(), alexandria_.end(),
                                    [atype](Ffatype const &f)
                                    { return (atype == f.getType()); });
            if (atp == alexandria_.end())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such atom ype %s", atype.c_str()).c_str()));
            }
            return atp;
        }

        /*! \brief
         * Return the const_iterator corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         */
        FfatypeConstIterator findAtype(const std::string &atype) const
        {
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [atype](Ffatype const &f)
                                { return (atype == f.getType()); });
        }

        VsiteIterator getVsiteBegin()  { return vsite_.begin(); }

        VsiteConstIterator getVsiteBegin()  const { return vsite_.begin(); }

        VsiteIterator getVsiteEnd() { return vsite_.end(); }

        VsiteConstIterator getVsiteEnd() const { return vsite_.end(); }

        VsiteIterator findVsite(std::string  atype)
        {

            return std::find_if(vsite_.begin(), vsite_.end(),
                                [atype](const Vsite &vs)
                                {
                                    return (atype == vs.atype());
                                });
        }

        VsiteConstIterator findVsite(std::string atype) const
        {

            return std::find_if(vsite_.begin(), vsite_.end(),
                                [atype](const Vsite &vs)
                                {
                                    return (atype == vs.atype());
                                });
        }

        /*! \brief
         * Return the poltype corresponding to atype and true if successful
         *
         * \param[in]  atype  Atom type
         * \param[out] ptype  Polarizability type.
         */
        bool atypeToPtype(const std::string &atype,
                          std::string       *ptype) const;


        /*! \brief
         * Return the bond type corresponding to atom type and true if successful
         *
         * \param[in]  atype  Atom type
         * \param[out] btype  Polarizability type.
         */
        bool atypeToBtype(const std::string &atype,
                          std::string       *btype) const;

        /*! \brief
         * Return the zeta type corresponding to atom type and true if successful
         *
         * \param[in]  atype  Atom type
         * \param[out] ztype  Zeta type.
         */
        bool atypeToZtype(const std::string &atype,
                          std::string       *ztype) const;
 
        /*! \brief Convenience function to create a zeta Identifier
         *
         * Routine will convert atom types to a zeta identifier
         * \param[in] atoms Vector of atom types
         * \return BCC identifier
         * \throws if something wrong
         */
        const Identifier atomtypesToZetaIdentifier(const std::vector<std::string> atoms) const;

        /*! \brief
         * Add a new force list. The routine checks for duplicate itypes
         * and will not add a second block of a certain type is one is
         * present already.
         * \param[in] interaction The type of interaction being modeled
         * \param[in] forces      The data structure
         */
        void addForces(const std::string             &interaction,
                       const ForceFieldParameterList &forces);

        size_t nforces() const { return forces_.size(); }

        std::map<InteractionType, ForceFieldParameterList> *forces() { return &forces_; }
    
        const std::map<InteractionType, ForceFieldParameterList> &forcesConst() const { return forces_; }

        /*! \brief Test whether an interaction is present
         * \param[in] iType The interaction type
         * \return true if present
         */
        bool interactionPresent(InteractionType iType) const { return forces_.find(iType) != forces_.end(); }
        
        /*! \brief Return forces
         * \param[in] iType The interaction type
         * \return parameter list
         * \throws if the list does not exist
         */
        ForceFieldParameterList *findForces(InteractionType iType)
        {
            auto force = forces_.find(iType);
            if (force == forces_.end())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("No such interaction type %s (there are %d interaction types)",
                                                               interactionTypeToString(iType).c_str(), 
                                                               static_cast<int>(forces_.size())).c_str()));
            }
            return &force->second;
        }

        const ForceFieldParameterList &findForcesConst(InteractionType iType) const
        {
            const auto force = forces_.find(iType);
            if (force == forces_.end())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("No such interaction type %s (there are %d interaction types)",
                                                               interactionTypeToString(iType).c_str(),
                                                               static_cast<int>(forces_.size())).c_str()));
            }
            return force->second;
        }

        const std::string &getVsite_angle_unit() const { return vsite_angle_unit_; }

        const std::string &getVsite_length_unit() const { return vsite_length_unit_; }

        void addSymcharges(const std::string &central,
                           const std::string &attached,
                           int                numattach);

        SymchargesIterator getSymchargesBegin() { return symcharges_.begin(); }

        SymchargesIterator getSymchargesEnd() { return symcharges_.end(); }

        SymchargesConstIterator getSymchargesBegin() const { return symcharges_.begin(); }

        SymchargesConstIterator getSymchargesEnd() const { return symcharges_.end(); }

        //! Return whether there is polarizability support for atype
        int havePolSupport(const std::string &atype) const;

        const char *getOpts(const std::string &name) const;

        //! Return the charge distribution type used
        ChargeType chargeType() const { return ChargeType_; }
        
        //! Set the charge distribution type used
        void setChargeType(ChargeType eqdType) { ChargeType_ = eqdType; }
        
        //! Set the charge distribution model used
        void setChargeType(const std::string &eqdType)
        { 
            ChargeType_ = name2ChargeType(eqdType);
        }
        
        //! Return the charge generation algorithm used
        ChargeGenerationAlgorithm chargeGenerationAlgorithm() const
        { return ChargeGenerationAlgorithm_; }
        
        //! Set the charge generation algorithm used
        void setChargeGenerationAlgorithm(ChargeGenerationAlgorithm eqgAlgorithm)
        { ChargeGenerationAlgorithm_ = eqgAlgorithm; }
        
        //! Set the charge generation algorithm used
        void setChargeGenerationAlgorithm(const std::string &eqgAlgorithm)
        {
            ChargeGenerationAlgorithm_ = name2ChargeGenerationAlgorithm(eqgAlgorithm);
        }

        //! Return whether the model used is polarizable.     
        bool polarizable() const { return polarizable_; }

        //! Turn polarizability on or off.
        void setPolarizable(bool polar) { polarizable_ = polar; }

        //! Turn polarizability on or off automatically
        void checkForPolarizability();

        //! Spread from master to slave nodes
        void  broadcast(const t_commrec *cr);

        //! Spread eemprop from master to slave nodes
        void broadcast_eemprop(const t_commrec *cr);
        
        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

        //! \brief Check internal consistency of data structures
        void checkConsistency(FILE *fplog) const;
    private:
        std::string                           filename_;
        std::vector<Ffatype>                  alexandria_;
        std::vector<Vsite>                    vsite_;
        std::vector<std::string>              btype_;
        std::string                           alexandriaVersion_;
        std::string                           vsite_angle_unit_;
        std::string                           vsite_length_unit_;
        int                                   nexcl_ = 0;
        double                                gtEpsilonR_ = 1.0;
        std::map<InteractionType, ForceFieldParameterList> forces_;
        std::vector<Symcharges>               symcharges_;
        bool                                  polarizable_ = false;
        ChargeType                            ChargeType_  = eqtPoint;
        ChargeGenerationAlgorithm             ChargeGenerationAlgorithm_ = eqgEEM;

        void addBtype(const std::string &btype);

        gmx_bool strcasestrStart(std::string needle, std::string haystack);

        template<class Type>
        int indexOfPointInVector(Type * pointer, std::vector<Type> vector)
        {
            return (pointer - &(vector[0]));
        }
};

}
#endif
