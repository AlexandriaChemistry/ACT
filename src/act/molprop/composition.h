/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022,2026
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

#ifndef COMPOSITION_H
#define COMPOSITION_H

#include <string.h>

#include <map>
#include <string>
#include <vector>

#include "act/utility/communicationrecord.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

/*! \brief
 * Specifies the name of an atom type and the number in a molecular composition
 *
 * \inpublicapi
 * \ingroup group_molprop
 */
class AtomNum
{
    private:
        /*! \brief
         * Atom name
         */
        std::string catom_;
        /*! \brief
         * Atom number
         */
        int         cnumber_;
    public:
        //! Default constructor
        AtomNum() = default;

        /*! \brief
         * Creates a new AtomNum object.
         *
         * \param[in] catom   Atom name
         * \param[in] cnumber Number of copies of this atom
         */
        AtomNum(const char *catom, int cnumber) { SetAtom(catom); SetNumber(cnumber); }

        /*! \brief
         * Creates a new AtomNum object.
         *
         * \param[in] catom   Atom name
         * \param[in] cnumber Number of copies of this atom
         */
        AtomNum(const std::string &catom, int cnumber)
        {
            SetAtom(catom);
            SetNumber(cnumber);
        }

        /*! \brief
         * Return the name of the atom for this AtomNum
         */
        const std::string &getAtom() const { return catom_; }

        /*! \brief
         * Set the name of the atom for this AtomNum
         */
        void SetAtom(const std::string &catom) { catom_ = catom; }

        /*! \brief
         * Set the name of the atom for this AtomNum
         */
        void SetAtom(const char *catom) { catom_.assign(catom); }

        /*! \brief
         * Return the number of atoms for this AtomNum
         */
        int getNumber() const { return cnumber_; }

        /*! \brief
         * Set the number of atoms for this AtomNum
         */
        void SetNumber(int cnumber) { cnumber_ = cnumber; }

};

//! Iterator over a vector of AtomNum
using  AtomNumIterator      = typename std::vector<AtomNum>::iterator;
//! Const iterator over a vector of AtomNum
using  AtomNumConstIterator = typename std::vector<AtomNum>::const_iterator;

/*! \brief
 * Contains the molecular composition in terms of atoms
 *
 * \inpublicapi
 * \ingroup group_molprop
 */
class MolecularComposition
{
    private:
        /*! \brief
         * Composition name
         */
        std::string          compname_;
        /*! \brief
         * A vector of AtomNum object
         */
        std::vector<AtomNum> atomnum_;
    public:
        //! Defult constructor
        MolecularComposition() = default;

        /*! \brief
         * Creates a new MolecularComposition object.
         *
         * \param[in] compname  Name of the composition type
         */
        MolecularComposition(const char *compname) { compname_.assign(compname); }

        /*! \brief
         * Creates a new MolecularComposition object.
         *
         * \param[in] compname  Name of the composition type
         */
        MolecularComposition(const std::string &compname)
        {
            compname_ = compname;
        }

        /*! \brief
         * Return the composition name
         */
        const std::string &getCompName() const { return compname_; }

        /*! \brief
         * Set the composition name
         */
        void SetCompName(const std::string &compname)
        {
            compname_ = compname;
        }

        /*! \brief
         * Add an AtomNum object to the composition
         *
         * \param[in] an  Atom number
         */
        void AddAtom(const AtomNum &an);

        /*! \brief
         * Remove the atom with name catom from the composition
         *
         * \param[in] catom Atom name
         */
        void DeleteAtom(const std::string &catom);

        /*! \brief
         * Replace the oldatom by newatom
         *
         * \param[in] oldatom   Name of the old atom
         * \param[in] newatom   Name of the new atom
         */
        void ReplaceAtom(const std::string &oldatom,
                         const std::string &newatom);

        //! Return const vector of AtomNum
        const std::vector<AtomNum> atomNumConst() const { return atomnum_; }
    
        /*! \brief
         * Return const iterator pointing to a specific atom or EndAtomNum if not found
         *
         * \param[in] an Atom number
         */
        AtomNumConstIterator searchAtomConst(const std::string &an) const;

        /*! \brief
         * Return iterator pointing to a specific atom or EndAtomNum if not found
         *
         * \param[in] an Atom number
         */
        AtomNumIterator searchAtom(const std::string &an);

        /*! \brief
         * Return the number of atoms of a certain type
         *
         **\param[in] atom Atom name
         */
        int CountAtoms(const std::string &atom) const;

        /*! \brief
         * Return the total number of atoms
         */
        int CountAtoms() const;

};
//! Iterator over MolecularComposition items
using MolecularCompositionIterator      = typename std::vector<MolecularComposition>::iterator;
//! Const iterator over MolecularComposition items
using MolecularCompositionConstIterator = typename std::vector<MolecularComposition>::const_iterator;

/*! \brief
 * Contains data on an atom based on a calculation.
 * Coordinates and Forces are store in the specified units.
 * This class coordinates, name, atom type and an array of
 * Charge values based on different methods.
 *
 * \inpublicapi
 * \ingroup group_molprop
 */
class CalcAtom
{
private:
    //! Atom name
    std::string             name_;
    //! OpenBabel atom type
    std::string             obType_;
    //! Unit of the coordinates
    std::string             coord_unit_;
    //! Unit of the forces
    std::string             force_unit_;
    //! Residue name
    std::string             residueName_;
    //! X coordinate, in coord_unit_
    double                  x_ = 0;
    //! Y coordinate, in coord_unit_
    double                  y_ = 0;
    //! Z coordinate, in coord_unit_
    double                  z_ = 0;
    //! X component of the force, in force_unit_
    double                  fx_ = 0;
    //! Y component of the force, in force_unit_
    double                  fy_ = 0;
    //! Z component of the force, in force_unit_
    double                  fz_ = 0;
    //! Atom id
    int                     atomID_ = 0;
    //! Residue number
    int                     residueNumber_ = 0;
    //! Chain id
    int                     chainId_ = 0;
    //! Chain identifier
    char                    chain_ = ' ';
    //! Atomic charges (e) indexed by charge type
    std::map<std::string, double> q_;

    //! Return chain ID
    int chainId() const { return chainId_; }

public:
    //! Default constructor
    CalcAtom() = default;
    
    //! Constructor initiating the name, type and atomid
    CalcAtom(const std::string &name,
             const std::string &obtype,
             int                atomid)
    {
        name_ = name; obType_ = obtype; atomID_ = atomid;
    };
    
    //! Function returning true if the two atoms are equal
    bool Equal(const CalcAtom &ca);
    
    /*! \brief Add an AtomicCharge element to the atom
     * If a charge of the type is present already, it will be
     * overwritten.
     * \param[in] type Charge type
     * \param[in] q    The charge
     */
    void AddCharge(const std::string &type, double q)
    {
        q_.insert_or_assign(type, q);
    }
    
    /*! \brief Return whether the charge type is present
     * \param[in] type Charge type
     * \return true if found
     */
    bool hasCharge(const std::string &type) const
    {
        return q_.find(type) != q_.end();
    }
    
    /*! \brief Return charge of a certain type
     * \param[in] type Charge type
     * \return charge if found
     * \throws if not found
     */
    double charge(const std::string &type) const
    {
        const auto qq = q_.find(type);
        if (q_.end() == qq)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("No charge of type %s", type.c_str()).c_str()));
        }
        return qq->second;
    }
    
    //! \brief Return the whole charge map
    const std::map<std::string, double> &chargesConst() const { return q_; }
    
    //! Return the atom id of the atom
    int getAtomid() const { return atomID_; }
    
    //! Set the atom id of the atom
    void setAtomid(int atomId) { atomID_ = atomId; }
    
    //! Return the name of the atom
    const std::string &getName() const { return name_; }
    
    //! Set the name of the atom
    void setName(const std::string &name) { name_ = name; }
    
    //! Return the OpenBabel type of the atom
    const std::string &getObtype() const { return obType_; }
    
    //! Set the OpenBabel type of the atom
    void setObtype(const std::string &obtype) { obType_ = obtype; }
    
    //! Return the unit of the coordinates of the atom
    const std::string &coordUnit() const { return coord_unit_; }
    
    //! Return the unit of the coordinates of the atom
    const std::string &forceUnit() const { return force_unit_; }
    
    //! Return the name of residue
    const std::string &residueName() const { return residueName_; }
    
    //! Return the residue number
    int residueNumber() const { return residueNumber_; }

    //! Set the unit of the coordinates of the atom
    void setCoordUnit(const std::string &unit);
    
    //! Set the unit of the forces on the atom
    void setForceUnit(const std::string &unit);
    
    //! Set the residue name for the atom
    void setResidue(const std::string &residueName, int residueNumber)
    {
        residueName_   = residueName;
        residueNumber_ = residueNumber;
    }
    
    //! Return chain
    char chain() const { return chain_; }
    
    //! Set the coordinates of the atom
    void setCoords(double x, double y, double z) { x_ = x; y_ = y; z_ = z; }
    
    //! Return all the coordinates of the atom
    void coords(double *x, double *y, double *z) const
    { *x = x_; *y = y_; *z = z_; }

    //! Set the forces on the atom
    void setForces(double fx, double fy, double fz) { fx_ = fx; fy_ = fy; fz_ = fz; }
    
    //! Return all the coordinates of the atom
    void forces(double *fx, double *fy, double *fz) const
    { *fx = fx_; *fy = fy_; *fz = fz_; }

    //! Return the X coordinate of the atom
    double getX() const { return x_; }
    
    //! Return the Y coordinate of the atom
    double getY() const { return y_; }
    
    //! Return the Z coordinate of the atom
    double getZ() const { return z_; }
    
    /*! \brief
     * Broadcast this object over an MPI connection
     *
     * \param[in] cr   data structure for MPI communication
     * \param[in] root The root of this communication
     * \param[in] comm MPI communicator
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);

};

//! Iterator over CalcAtom items
using  CalcAtomIterator      = typename std::vector<CalcAtom>::iterator;

//! Const iterator over CalcAtom items
using  CalcAtomConstIterator = typename std::vector<CalcAtom>::const_iterator;

} // namespace alexandria

#endif
