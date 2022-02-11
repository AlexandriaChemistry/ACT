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
#include <map>
#include <vector>

#include "gromacs/utility/stringutil.h"

#include "act/basics/chargemodel.h"
#include "act/basics/interactiontype.h"
#include "act/poldata/forcefieldparameter.h"
#include "act/poldata/forcefieldparameterlist.h"
#include "act/poldata/particletype.h"
#include "act/poldata/poldata_low.h"
#include "act/poldata/vsite.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/stringutil.h"

namespace alexandria
{

class Poldata
{
    public:

        //! Default constructor
        Poldata() {};

        /*! \brief
         * Set the file name gentop.dat
         *
         */
        void  setFilename(const std::string &fn2);

        //! Return the actual filename
        const std::string &filename() const { return filename_; }

        /*! \brief
         * Set the force field checksum
         * \param[in] version Force field checksum
         */
        void setCheckSum(const std::string &checkSum)
        {
            checkSum_ = checkSum;
        }
        
        /*! \brief Verify the checksum and print message if required
         * \param[in] fp       A file to write messages to, may be nullptr
         * \param[in] checkSum The checksum computed for the Poldata
         * \return true if checkSum matches.
         */
        bool verifyCheckSum(FILE              *fp,
                            const std::string &checkSum);
    
        /*! \brief Verify the checksum and print message if required
         * \param[in] fp       A file to write messages to, may be nullptr
         * \return true if checkSum is correct.
         */
        bool verifyCheckSum(FILE *fp);
    
        /*! \brief Update the checkSum
         * Compute a new checkSum and store it in the appropriate data field.
         */
        void updateCheckSum();
        /*! \brief
         * Get the force field checkSum
         * \return Force field checkSum
         */
        const std::string &checkSum() const
        {
            return checkSum_;
        }

        /*! \brief Generate a new time stamp with the current time.
         */
        void updateTimeStamp();

        //! \return the time stamp for this force field file
        const std::string &timeStamp() const { return timeStamp_; }

        /*! \brief 
         * Set the time stamp for this force field file
         * \param[in] timeStamp the new value
         */
        void setTimeStamp(const std::string &timeStamp) { timeStamp_ = timeStamp; }
        //! Return whether or not the original Rappe & Goddard is used
        bool rappe() const;
        
        //! Return whether Yang & Sharp is used
        bool yang() const;
        
        /*! \brief
         * Set relative dielectric constant
         *
         */
        void setEpsilonR(double epsilonR) { gtEpsilonR_ = epsilonR; }

        /*! \brief
         * Set the number of exclusion
         */
        void setNexcl(int nexcl) { nexcl_ = nexcl; }

        /*! \brief Add an atom type
         * \param[in] ptp The new particle type
         * \throw if this particle type exists already
         */
        void  addParticleType(const ParticleType &ptp);
        /*! \brief
         *  Add a virtual site type
         *
         * \param[in] atype           The name specifying the vsite type
         * \param[in] type            Unknown
         * \param[in] number
         * \param[in] distance
         * \param[in] angle
         * \param[in] ncontrolatoms
         */
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

        std::vector<Vsite> &getVsite() { return vsite_; }

        const std::vector<Vsite> &getVsiteConst() const { return vsite_; }

        int getNexcl() const { return nexcl_; }

        size_t getNatypes() const { return alexandria_.size(); }

        /*! \brief
         * Return the relative dielectric constant
         */
        double getEpsilonR() const { return gtEpsilonR_; }

        /*! \brief Compute dependent parameters
         *
         * The force field may contain a number of parameters
         * that depend on othe parameters. For instance, the
         * distance r_{ik} in an angle triplet r_{ijk} depends
         * on the distances r_{ij} and r_{kj} as well as the 
         * angle between the vectors. Those dependent properties
         * are computed here in order to get the highest accuracy
         * in the numbers. This is important when computing energies
         * based on these distances.
         */
        void calcDependent();

        /*! \brief
         * Return the element name corresponding to the zeta type
         *
         * \param[in] ztype  zeta Type
         */
        const std::string ztype2elem(const std::string &ztype) const;

    /*! \brief Check whether particle type exists
     * \param[in] id The identifier of the particle type
     * \return true if ptype exists
     */
    bool hasParticleType(const Identifier &id) const
    {
        auto atp = std::find_if(alexandria_.begin(), alexandria_.end(),
                                [id](ParticleType const &f)
                                { return (id == f.id()); });
        return (atp != alexandria_.end());
    }
    /*! \brief Check whether particle type exists
     * \param[in] ptype The name of the particle type
     * \return true if ptype exists
     */
    bool hasParticleType(const std::string &ptype) const
    {
        Identifier id(ptype);
        return hasParticleType(id);
    }
        
        ParticleTypeIterator findParticleType(const Identifier &id)
        {
            auto atp = std::find_if(alexandria_.begin(), alexandria_.end(),
                                    [id](ParticleType const &f)
                                    { return (id == f.id()); });
            if (atp == alexandria_.end())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such atom type %s", id.id().c_str()).c_str()));
            }
            return atp;
        }
        ParticleTypeConstIterator findParticleType(const Identifier &id) const
        {
            auto atp = std::find_if(alexandria_.begin(), alexandria_.end(),
                                    [id](ParticleType const &f)
                                    { return (id == f.id()); });
            if (atp == alexandria_.end())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such atom ype %s", id.id().c_str()).c_str()));
            }
            return atp;
        }
        /*! \brief
         * Return the iterator corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         * \throw if atom type not  found.
         */
        ParticleTypeIterator findParticleType(const std::string &atype)
        {
            Identifier id(atype);
            return findParticleType(id);
        }

        /*! \brief
         * Return the const_iterator corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         */
        ParticleTypeConstIterator findParticleType(const std::string &atype) const
        {
            Identifier id(atype);
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [id](ParticleType const &f)
                                { return (id == f.id()); });
        }
        
        /*! \brief Return mutable vector \p alexandria_
         */
        std::vector<ParticleType> *particleTypes() { return &alexandria_; }

        //! \return the number of particle types
        int nParticleTypes() const { return alexandria_.size(); }
        
        /*! \brief Return const vector \p alexandria_
         */
        const std::vector<ParticleType> &particleTypesConst() const { return alexandria_; }

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
    //const Identifier atomtypesToZetaIdentifier(const std::vector<std::string> atoms) const;

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
        
        /*! \brief Return forces in a mutable form
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
            const auto &force = forces_.find(iType);
            if (force == forces_.end())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("No such interaction type %s (there are %d interaction types)",
                                                               interactionTypeToString(iType).c_str(),
                                                               static_cast<int>(forces_.size())).c_str()));
            }
            return force->second;
        }
        
        /*! Check whether an interaction type corresponding to parameter 
         *  type exists, and return it.
         * \param[in]  type  The string describing e.g. "zeta"
         * \param[out] itype The itype
         * \return whether or not the type was found.
         */
        bool typeToInteractionType(const std::string &type, InteractionType *itype);

        const std::string &getVsite_angle_unit() const { return vsite_angle_unit_; }

        const std::string &getVsite_length_unit() const { return vsite_length_unit_; }

        void addSymcharges(const std::string &central,
                           const std::string &attached,
                           int                numattach);

        //! \return a constant reference of \p symcharges_
        const std::vector<Symcharges> &getSymcharges() const { return symcharges_; }

        SymchargesIterator getSymchargesBegin() { return symcharges_.begin(); }

        SymchargesIterator getSymchargesEnd() { return symcharges_.end(); }

        SymchargesConstIterator getSymchargesBegin() const { return symcharges_.begin(); }

        SymchargesConstIterator getSymchargesEnd() const { return symcharges_.end(); }

        //! Return the charge generation algorithm used
        ChargeGenerationAlgorithm chargeGenerationAlgorithm() const;
        
        //! Return whether the model used is polarizable.     
        bool polarizable() const { return polarizable_; }

        //! Turn polarizability on or off.
        void setPolarizable(bool polar) { polarizable_ = polar; }

        //! Turn polarizability on or off automatically
        void checkForPolarizability();

        /*! Spread from master to other nodes
         * \param[in] cr Communication record
         */
        void sendToHelpers(const CommunicationRecord *cr);

        /*! Spread eemprop to a helper node
         * \param[in] cr   Communication record
         * \param[in] dest Destination node
         */
        void sendEemprops(const CommunicationRecord *cr, int dest);
    
        /*! Receive eemprop from someone
         * \param[in] cr  Communication record
         * \param[in] src Source node
         */
        void receiveEemprops(const CommunicationRecord *cr, int src);
        
        /*! Spread mutable particle properties to a helper node
         * \param[in] cr   Communication record
         * \param[in] dest Destination node
         */
        void sendParticles(const CommunicationRecord *cr, int dest);
        
        /*! Receive particles from someone
         * \param[in] cr  Communication record
         * \param[in] src Source node
         */
        void receiveParticles(const CommunicationRecord *cr, int src);
        
        CommunicationStatus Send(const CommunicationRecord *cr, int dest);

        CommunicationStatus Receive(const CommunicationRecord *cr, int src);

        //! \brief Check internal consistency of data structures
        void checkConsistency(FILE *fplog) const;

        //! \return a constant \p type2Itype_ reference
        const std::map<std::string, InteractionType> &type2Itype() const { return type2Itype_; }

    private:
        std::string                           checkSum_;
        std::string                           timeStamp_;
        std::map<std::string, InteractionType> type2Itype_;
        std::string                           filename_;
        std::vector<ParticleType>             alexandria_;
        std::vector<Vsite>                    vsite_;
        std::string                           vsite_angle_unit_;
        std::string                           vsite_length_unit_;
        int                                   nexcl_ = 0;
        double                                gtEpsilonR_ = 1.0;
        std::map<InteractionType, ForceFieldParameterList> forces_;
        std::vector<Symcharges>               symcharges_;
        bool                                  polarizable_ = false;
        ChargeGenerationAlgorithm             ChargeGenerationAlgorithm_ = ChargeGenerationAlgorithm::EEM;

        gmx_bool strcasestrStart(std::string needle, std::string haystack);

        template<class Type>
        int indexOfPointInVector(Type * pointer, std::vector<Type> vector)
        {
            return (pointer - &(vector[0]));
        }

};

}
#endif
