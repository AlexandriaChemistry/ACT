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
#ifndef TUNE_FC_UTILS_H
#define TUNE_FC_UTILS_H

#include <map>
#include <string>
#include <vector>

#include "gromacs/utility/stringutil.h"

#include "communication.h"
#include "identifier.h"
#include "mymol.h"
#include "plistwrapper.h"
#include "poldata.h"

struct t_commrec;

namespace alexandria
{

/*! \brief Helper class storing bond/angle/dihedral names
 *
 * For one bond/angle/dihedral here the name of the bondtypes
 * are stored as in e.g. c c h for an angle, along with the number
 * of occurrences in the force field.
 */
class ParameterNames
{
    public:

    //! Empty constructor
    ParameterNames () {}

    /*! \brief Constructor with initiation
     *
     * \param[in] ncopies  Number of copies of this interaction
     * \param[in] ftype    The GROMACS function type
     * \param[in] params   String with all the parameters
     * \param[in] index    Integer identifier
     */
        ParameterNames(int                        ncopies,
                       int                        ftype,
                       const std::vector<double> &params,
                       int                        index) :
              ncopies_(ncopies),
              ftype_(ftype),
              params_(params),
              poldataIndex_(index)
        {}

        void inc() { ncopies_++; }

        int nCopies() const { return ncopies_; }

        /*! \brief Update the parameters for this interaction
         * \param[in] params The new parameters
         */
        void setParameterValues(const std::vector<double> &params)
        {
            params_ = params;
        }

        int poldataIndex() const { return poldataIndex_; }

        const std::vector<double> &parameterValues() const { return params_; }

        size_t nParams() const { return params_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:
        //! Number of copies in the molecule data set
        int                 ncopies_;
        //! Function type for this particular bond
        int                 ftype_;
        //! Vector containing all the parameters
        std::vector<double> params_;
        //! Index in Poldata structure
        int                 poldataIndex_;
};

using ParameterNamesIterator = typename std::vector<ParameterNames>::iterator;

/*! \brief Class holding for one type of interactions all names
 *
 * Class holding the OptNames for each interaction type.
 */
class ForceConstants
{

    public:

        ForceConstants () {}

    ForceConstants(int             ftype, 
                   InteractionType itype, 
                   bool            bOpt)
            :
              ftype_(ftype),
              itype_(itype),
              bOpt_(bOpt)
        {
        }

        /*! \brief Add a force constant to the map
         *
         * \param[in] identifier A string to identify the parameter
         * \param[in] bn         ParameterNames structure
         */
        void addForceConstant(const Identifier &identifier, ParameterNames bn)
        { 
            bn_.insert({identifier, std::move(bn)}); 
        }

        /*! \brief
         * Extract information from the idef structure about parameters
         * \param[in] mm All the molecule structures
         * \param[in] pd The Poldata structure
         */
        void analyzeIdef(const std::vector<MyMol> &mm,
                         const Poldata            *pd);

        /*! \brief Make reverse index from Poldata to ParameterNames
         *
         * The ParameterNames structure stores the Poldata index for
         * all interactions. This routine makes an index to convert
         * the Poldata index to the index in ParameterNames.
         */
        void makeReverseIndex();

        int reverseIndex(int poldataIndex)
        {
            GMX_RELEASE_ASSERT(poldataIndex >= 0 && poldataIndex < static_cast<int>(reverseIndex_.size()), "Incorrect poldataIndex");
            GMX_RELEASE_ASSERT(reverseIndex_[poldataIndex] != -1, "The reverseIndex is incorrect");

            return reverseIndex_[poldataIndex];
        }

        int ftype() const { return ftype_; }

        InteractionType interactionType() const { return itype_; }

        void dump(FILE *fp) const;

        /*! \brief Find a particular bond name
         * 
         * \param[in] identifier String identifying the bond name
         * \throw if identifier not found
         */
        const ParameterNames &bondNamesConst(const Identifier &identifier) const;

        /*! \brief Test whether bondname exists
         * \param[in] identifier identifying the bond name
         * \return existence of the bond name
         */
        bool bondNameExists(const Identifier &identifier) const
        {
            return bn_.find(identifier) != bn_.end();
        }

        /*! \brief Return the whole ParameterNames map const version
         */
        const std::map<Identifier, ParameterNames> &bondNamesConst() const { return bn_; }
        /*! \brief Return the whole ParameterNames map
         */
        std::map<Identifier, ParameterNames> &bondNames() { return bn_; }

        size_t nbad() const { return bn_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:
    //! GROMASCS function type
        int                                  ftype_;
    //! ACT interaction type
        InteractionType                      itype_;
    //! Whether this parameter is to be optimized (changed)
        bool                                 bOpt_;
    //! Map from 
        std::map<Identifier, ParameterNames> bn_;
        std::vector<int>                     reverseIndex_;
};


class PoldataUpdate
{
public:
    /*! Constructor without parameters 
     */
    PoldataUpdate() {}
    /*! Empty constructor
     * \param[in] iType           Interaction type, e.g. ietBONDS
     * \param[in] identifier      Atom / Bond etc. indentifier
     * \param[in] parameterValues Parameter values
     */
    PoldataUpdate(InteractionType            iType,
                  const Identifier           identifier,
                  const std::vector<double> &parameterValues) :
        iType_(iType), identifier_(identifier), 
        parameterValues_(parameterValues)
    {}
    
    /*! \brief
     * Implement the changes in Poldata.
     *
     * \param[inout] pd The Poldata structure
     */
    void execute(Poldata *pd);
    /*! \brief
     * Dump the contents of the structure to a file.
     * \param[in] fp File pointer to dumpt to if not nullptr
     */
    void dump(FILE *fp) const;
    
    CommunicationStatus Send(t_commrec *cr, int dest);
    
    CommunicationStatus Receive(t_commrec *cr, int src);

private:
    //! Interaction type, e.g. ietBONDS
    InteractionType     iType_;
    //! Atom / Bond etc. indentifier
    Identifier          identifier_;
    //! String containing the parameter values
    std::vector<double> parameterValues_;
};


} // namespace alexandria
#endif
