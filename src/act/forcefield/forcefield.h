/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
 
#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include "gromacs/utility/stringutil.h"

#include "act/basics/chargemodel.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/forcefield/particletype.h"
#include "act/forcefield/symcharges.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/stringutil.h"

namespace alexandria
{

class MsgHandler;

/*! \brief The complete force field
 *
 * Holds the particle types, the parameter lists for all interaction
 * types, the symmetric charges and general properties such as the
 * charge generation algorithm and polarizability.
 */
class ForceField
{
    public:

        //! Default constructor
        ForceField() {};

        /*! Return information about the force field
         * \return array of strings with info
         */
        std::vector<std::string> info() const;

        /*! \brief
         * Set the file name gentop.dat
         *
         * \param[in] fn2 The new file name, must not be empty
         */
        void  setFilename(const std::string &fn2);

        //! \return the actual filename
        const std::string &filename() const { return filename_; }

        /*! \brief
         * Set the force field checksum
         * \param[in] checkSum New Force field checksum
         */
        void setCheckSum(const std::string &checkSum)
        {
            checkSum_ = checkSum;
        }
        
        /*! \brief Verify the checksum and print message if required
         * \param[in] msgHandler A MsgHandler to write messages to, may be nullptr
         * \param[in] checkSum The checksum computed for the ForceField
         * \return true if checkSum matches.
         */
        bool verifyCheckSum(MsgHandler        *msgHandler,
                            const std::string &checkSum);
    
        /*! \brief Verify the checksum and print message if required
         * \param[in] msgHandler A MsgHandler to write messages to, may be nullptr
         * \return true if checkSum is correct.
         */
        bool verifyCheckSum(MsgHandler *msgHandler);
    
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

        /*! \brief Add an atom type
         * \param[in] ptp The new particle type
         * \throw if this particle type exists already
         */
        void  addParticleType(const ParticleType &ptp);

        //! \return the number of particle types
        size_t getNatypes() const { return alexandria_.size(); }

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
         * \return the element name, fatal error if not found
         */
        const std::string ztype2elem(const std::string &ztype) const;

    /*! \brief Check whether particle type exists
     * \param[in] id The identifier of the particle type
     * \return true if ptype exists
     */
    bool hasParticleType(const Identifier &id) const
    {
        return alexandria_.end() != alexandria_.find(id);
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
        
    /*! \brief Find a particle type for editing
     * \param[in] id The identifier of the particle type
     * \return pointer to the particle type
     * \throw if the particle type is not found
     */
    ParticleType *findParticleType(const Identifier &id)
    {
        auto atp = alexandria_.find(id);
        if (atp == alexandria_.end())
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such atom type %s", id.id().c_str()).c_str()));
        }
        return &atp->second;
    }
    
    /*! \brief Find a particle type
     * \param[in] id The identifier of the particle type
     * \return const pointer to the particle type
     * \throw if the particle type is not found
     */
    const ParticleType *findParticleType(const Identifier &id) const
    {
        auto atp = alexandria_.find(id);
        if (atp == alexandria_.end())
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such atom ype %s", id.id().c_str()).c_str()));
        }
        return &atp->second;
    }
    
    /*! \brief
     * Return the iterator corresponding to the atom type
     *
     * \param[in] atype  Atom Type
     * \return pointer to the particle type
     * \throw if atom type not  found.
     */
    ParticleType *findParticleType(const std::string &atype)
    {
        Identifier id(atype);
        return findParticleType(id);
    }

    /*! \brief
     * Return the const_iterator corresponding to the atom type
     *
     * \param[in] atype  Atom Type
     * \return const pointer to the particle type, the atom type must exist
     */
    const ParticleType *findParticleType(const std::string &atype) const
    {
        Identifier id(atype);
        auto iter = alexandria_.find(id);
        return &iter->second;
    }
        
    /*! \brief Return mutable vector \p alexandria_
     * \return pointer to the map of particle types
     */
    std::map<Identifier, ParticleType> *particleTypes() { return &alexandria_; }

    //! \return the number of particle types
    int nParticleTypes() const { return alexandria_.size(); }
    
    /*! \brief Return const vector \p alexandria_
     * \return the map of particle types
     */
    const std::map<Identifier, ParticleType> &particleTypesConst() const { return alexandria_; }
    
    /*! \brief
     * Return the poltype corresponding to atype and true if successful
     *
     * \param[in]  atype  Atom type
     * \param[out] ptype  Polarizability type.
     * \return true if successful
     */
    bool atypeToPtype(const std::string &atype,
                      std::string       *ptype) const;
    
    
    /*! \brief
     * Return the bond type corresponding to atom type and true if successful
     *
     * \param[in]  atype  Atom type
     * \param[out] btype  Bond type.
     * \return true if successful
     */
    bool atypeToBtype(const std::string &atype,
                      std::string       *btype) const;
    
    /*! \brief
     * Return the zeta type corresponding to atom type and true if successful
     *
     * \param[in]  atype  Atom type
     * \param[out] ztype  Zeta type.
     * \return true if successful
     */
    bool atypeToZtype(const std::string &atype,
                      std::string       *ztype) const;
    
    /*! \brief
     * Add a new force list. The routine checks for duplicate itypes
     * and will not add a second block of a certain type is one is
     * present already.
     * \param[in] iType  The type of interaction being modeled
     * \param[in] forces The data structure
     */
    void addForces(InteractionType                iType,
                   const ForceFieldParameterList &forces);
    
    //! \return the number of interaction types with parameters
    size_t nforces() const { return forces_.size(); }
    
    //! \return pointer to the map of parameter lists per interaction type
    std::map<InteractionType, ForceFieldParameterList> *forces() { return &forces_; }
    
    //! \return the map of parameter lists per interaction type
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
            GMX_THROW(gmx::InternalError(gmx::formatString("No such interaction type %s (there are %zu interaction types)",
                                                           interactionTypeToString(iType).c_str(), 
                                                           forces_.size()).c_str()));
        }
        return &force->second;
    }
    
    /*! \brief Return forces in a const form
     * \param[in] iType The interaction type
     * \return parameter list
     * \throws if the list does not exist
     */
    const ForceFieldParameterList &findForcesConst(InteractionType iType) const
    {
        const auto &force = forces_.find(iType);
        if (force == forces_.end())
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("No such interaction type %s (there are %zu interaction types)",
                                                           interactionTypeToString(iType).c_str(),
                                                           forces_.size()).c_str()));
        }
        return force->second;
    }
    
    /*! Check whether an interaction type corresponding to parameter 
     *  type exists, and return it.
     * \param[in]  type  The string describing e.g. "zeta"
     * \param[out] itype The itype
     * \return whether or not the type was found.
     */
    bool typeToInteractionType(const std::string &type, InteractionType *itype) const;
    
    /*! \brief Add a symmetric charge entry unless an identical one exists
     * \param[in] central   The central atom type
     * \param[in] attached  The attached atom type
     * \param[in] numattach The number of attached atoms
     */
    void addSymcharges(const std::string &central,
                       const std::string &attached,
                       int                numattach);
    
    //! \return a constant reference of \p symcharges_
    const std::vector<Symcharges> &getSymcharges() const { return symcharges_; }
    
    //! \return iterator to the start of \p symcharges_
    SymchargesIterator getSymchargesBegin() { return symcharges_.begin(); }
    
    //! \return iterator to the end of \p symcharges_
    SymchargesIterator getSymchargesEnd() { return symcharges_.end(); }
    
    //! \return const iterator to the start of \p symcharges_
    SymchargesConstIterator getSymchargesBegin() const { return symcharges_.begin(); }
    
    //! \return const iterator to the end of \p symcharges_
    SymchargesConstIterator getSymchargesEnd() const { return symcharges_.end(); }
    
    //! \return the charge generation algorithm used
    ChargeGenerationAlgorithm chargeGenerationAlgorithm() const { return ChargeGenerationAlgorithm_; };

    /*! \brief Set the charge generation algorithm
     * \param cga The charge generation algorithm to use
     */
    void setChargeGenerationAlgorithm(ChargeGenerationAlgorithm cga)
    {
        ChargeGenerationAlgorithm_ = cga;
    }

    //! \brief Guess the algorithm based on force field content
    void guessChargeGenerationAlgorithm();

    //! \return whether the model used is polarizable.     
    bool polarizable() const { return polarizable_; }
    
    /*! \brief Turn polarizability on or off.
     * \param[in] polar Whether the model is polarizable
     */
    void setPolarizable(bool polar) { polarizable_ = polar; }
    
    //! Turn polarizability on or off automatically
    void checkForPolarizability();

    /*! \brief Replace virtual sites on top of atoms by a second charge
     * Force field files may contain virtual sites of type VSITE1 that only
     * carry a charge, on top of an atom. These are replaced by a second
     * charge of the atom type, see ParticleType::chargeVS, which takes the
     * charge, the width of the charge distribution and the row of the
     * virtual site. The virtual site types and their parameters are removed.
     * Virtual sites that have other interactions, for instance Van der Waals
     * parameters, are kept.
     * \param[in] msghandler For messages, may be nullptr
     */
    void convertVsite1(MsgHandler *msghandler);
    
    /*! Spread from master to other nodes
     * \param[in] cr    Communication record
     * \param[in] root  The root of this communication
     * \param[in] bcast Whether or not to use boradcast functionality
     * \param[in] msgHandler MsgHandler for debug output
     */
    void sendToHelpers(const CommunicationRecord *cr, int root, bool bcast=true,
                       MsgHandler *msgHandler=nullptr);
    
    /*! Spread eemprop to a helper node
     * \param[in] cr   Communication record
     * \param[in] dest Destination node
     * \param[in] msgHandler MsgHandler for debug output
     */
    void sendEemprops(const CommunicationRecord *cr, int dest,
                      MsgHandler *msgHandler=nullptr);
    
    /*! Receive eemprop from someone
     * \param[in] cr  Communication record
     * \param[in] src Source node
     * \param[in] msgHandler MsgHandler for debug output
     */
    void receiveEemprops(const CommunicationRecord *cr, int src,
                         MsgHandler *msgHandler=nullptr);
    
    /*! Spread mutable particle properties to a helper node
     * \param[in] cr   Communication record
     * \param[in] dest Destination node
     * \param[in] msgHandler MsgHandler for debug output
     */
    void sendParticles(const CommunicationRecord *cr, int dest,
                       MsgHandler *msgHandler=nullptr);
    
    /*! Receive particles from someone
     * \param[in] cr  Communication record
     * \param[in] src Source node
     * \param[in] msgHandler MsgHandler for debug output
     */
    void receiveParticles(const CommunicationRecord *cr, int src,
                          MsgHandler *msgHandler=nullptr);
    
    /*! \brief BroadCast a whole force field
     * \param[in] cr  Communication data structure
     * \param[in] root The root of this communication
     * \param[in] comm MPI communication structure
     * \param[in] msgHandler MsgHandler for debug output
     * \return The status of the whole thing
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm,
                                  MsgHandler                *msgHandler=nullptr);
        
    /*! \brief Send a whole force field to another processor
     * \param[in] cr   Communication data structure
     * \param[in] dest Processor id to send the data to
     * \return The status of the communication
     */
    CommunicationStatus Send(const CommunicationRecord *cr, int dest);
    
    /*! \brief Receive a whole force field from another processor
     * \param[in] cr  Communication data structure
     * \param[in] src Processor id to receive the data from
     * \param[in] msgHandler MsgHandler for debug output
     * \return The status of the communication
     */
    CommunicationStatus Receive(const CommunicationRecord *cr, int src,
                                MsgHandler *msgHandler=nullptr);
    
    /*! \brief Check internal consistency of data structures
     * \param[in] msgHandler MsgHandler for warnings, may be nullptr
     */
    void checkConsistency(MsgHandler *msgHandler) const;
    
    //! \return a constant \p type2Itype_ reference
    const std::map<InteractionType, std::set<std::string>> &type2Itype() const { return type2Itype_; }
    
private:
    //! The force field checksum
    std::string                           checkSum_;
    //! The time stamp of the force field
    std::string                           timeStamp_;
    //! Parameter names per interaction type, filled on demand by typeToInteractionType
    mutable std::map<InteractionType, std::set<std::string>> type2Itype_;
    //! The force field file name
    std::string                           filename_;
    //! The particle types
    std::map<Identifier, ParticleType>    alexandria_;
    //! The parameter lists per interaction type
    std::map<InteractionType, ForceFieldParameterList> forces_;
    //! The symmetric charges
    std::vector<Symcharges>               symcharges_;
    //! Whether the model is polarizable
    bool                                  polarizable_ = false;
    //! The charge generation algorithm
    ChargeGenerationAlgorithm             ChargeGenerationAlgorithm_ = ChargeGenerationAlgorithm::NONE;
    
    /*! \brief Check whether a string starts with another, ignoring case
     * \param[in] needle   The string to look for
     * \param[in] haystack The string to search in
     * \return true if \p needle is found at the start of \p haystack
     */
    gmx_bool strcasestrStart(const std::string &needle, const std::string &haystack);

    /*! \brief Compute the index of an element in a vector
     * \param[in] pointer Pointer to an element of \p vector
     * \param[in] vector  The vector
     * \return the index of \p pointer in \p vector
     */
    template<class Type>
    int indexOfPointInVector(const Type *pointer, const std::vector<Type> &vector)
    {
        return (pointer - &(vector[0]));
    }
};

template<typename T>
/*! \brief Extract option value from a force field component
 * \param[in]  pd    The force field
 * \param[in]  itype The InteractionType
 * \param[in]  name  The name of the option
 * \param[out] value The value found
 * \return true if found, false otherwise
 */
bool ffOption(const ForceField  &pd,
              InteractionType    itype,
              const std::string &name,
              T                 *value);

}
#endif
