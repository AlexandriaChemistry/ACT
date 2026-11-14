/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
 * Compact storage of the non-bonded pairs of a topology, see NonBondedPairs.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_TOPOLOGY_NONBONDEDPAIRS_H
#define ACT_TOPOLOGY_NONBONDEDPAIRS_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <vector>

#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/topology/actmol_low.h"

namespace alexandria
{

class ActAtom;
class ForceField;
class MsgHandler;

/*! \brief A non-bonded pair as the force routines use it
 * It offers the same functions as a TopologyEntry of a pair, with the
 * parameters and the identifier of the pair of particle types.
 */
class PairView
{
public:
    /*! \brief Constructor
     * \param[in] atoms  The two particles
     * \param[in] params The parameters of the type pair
     * \param[in] id     The identifier of the type pair
     * \param[in] type   The index of the type pair, see PairTypes
     */
    PairView(const std::array<int, 2> *atoms, const std::vector<double> *params, const Identifier *id,
             uint32_t type = 0)
        : atoms_(atoms), params_(params), id_(id), type_(type) {}

    //! \return the two particles
    const std::array<int, 2> &atomIndices() const { return *atoms_; }

    //! \return particle k of the pair, 0 or 1
    int atomIndex(size_t k) const { return (*atoms_)[k]; }

    //! \return the parameters, in internal units, empty for Coulomb
    const std::vector<double> &params() const { return *params_; }

    //! \return the identifier of the pair of particle types, empty for Coulomb
    const Identifier &id() const { return *id_; }

    //! \return the index of the type pair, see PairTypes, zero for Coulomb
    uint32_t typeIndex() const { return type_; }

    //! \return this pair, such that a range element can be used like a BaseContainer<TopologyEntry>
    const PairView *operator->() const { return this; }

private:
    //! The particles
    const std::array<int, 2>  *atoms_;
    //! The parameters
    const std::vector<double> *params_;
    //! The identifier
    const Identifier          *id_;
    //! The type pair
    uint32_t                   type_;
};

/*! \brief The pairs of particle types of one interaction type, see NonBondedPairs
 * The identifier depends only on the particle types of the two particles,
 * and so do the parameters, which the force field computes with the
 * combination rules.
 */
struct PairTypes
{
    //! The identifier of each type pair, empty if the particle types have none
    std::vector<Identifier>          ids;
    //! The parameters of each type pair in internal units, empty if the force field has none
    std::vector<std::vector<double>> params;
};

class NonBondedPairs;

/*! \brief A part of the pairs of an interaction type, as a TopologyEntryRange
 * The elements are PairView objects.
 */
class PairRange
{
public:
    //! An empty range
    PairRange() = default;

    /*! \brief Constructor
     * \param[in] pairs All pairs of the topology
     * \param[in] pos   The positions in pairs of the pairs of this interaction type
     * \param[in] type  Their type pairs
     * \param[in] types The type pairs, nullptr for Coulomb
     * \param[in] begin The first pair of this range in pos
     * \param[in] end   One beyond the last
     */
    PairRange(const std::array<int, 2> *pairs, const uint32_t *pos, const uint32_t *type,
              const PairTypes *types, size_t begin, size_t end)
        : pairs_(pairs), pos_(pos), type_(type), types_(types), begin_(begin), end_(end) {}

    //! \return the number of pairs
    size_t size() const { return end_ - begin_; }

    //! \return whether there are no pairs
    bool empty() const { return end_ == begin_; }

    //! \return pair k of this range
    PairView operator[](size_t k) const
    {
        size_t e = begin_ + k;
        if (types_)
        {
            return PairView(pairs_ + pos_[e], &types_->params[type_[e]], &types_->ids[type_[e]], type_[e]);
        }
        return PairView(pairs_ + pos_[e], &noParams_, &noId_);
    }

    /*! \brief The pairs from an offset on, as std::span::subspan
     * \param[in] offset The first pair
     * \return the range
     */
    PairRange subspan(size_t offset) const
    {
        PairRange r(*this);
        r.begin_ = begin_ + offset;
        return r;
    }

    /*! \brief The first pairs, as std::span::first
     * \param[in] n The number of pairs
     * \return the range
     */
    PairRange first(size_t n) const
    {
        PairRange r(*this);
        r.end_ = begin_ + n;
        return r;
    }

    //! An iterator over the pairs
    class Iterator
    {
    public:
        /*! \brief Constructor
         * \param[in] range The range
         * \param[in] k     The position in the range
         */
        Iterator(const PairRange *range, size_t k) : range_(range), k_(k) {}
        //! \return the pair
        PairView operator*() const { return (*range_)[k_]; }
        //! \return the next position
        Iterator &operator++() { k_ += 1; return *this; }
        //! \return whether two iterators differ
        bool operator!=(const Iterator &other) const { return k_ != other.k_; }
    private:
        //! The range
        const PairRange *range_;
        //! The position
        size_t           k_;
    };

    //! \return the first pair
    Iterator begin() const { return Iterator(this, 0); }

    //! \return beyond the last pair
    Iterator end() const { return Iterator(this, size()); }

private:
    //! All pairs of the topology
    const std::array<int, 2> *pairs_ = nullptr;
    //! The positions of the pairs of the interaction type
    const uint32_t           *pos_   = nullptr;
    //! Their type pairs
    const uint32_t           *type_  = nullptr;
    //! The type pairs
    const PairTypes          *types_ = nullptr;
    //! The first pair of the range in pos_
    size_t                    begin_ = 0;
    //! One beyond the last
    size_t                    end_   = 0;
    //! The parameters of Coulomb pairs, which take them from the particles
    static const std::vector<double> noParams_;
    //! The identifier of Coulomb pairs
    static const Identifier          noId_;
};

/*! \brief The non-bonded pairs of a topology
 * All non-bonded interaction types (Coulomb, Van der Waals and their
 * corrections) have the same pairs, those of the particles that are not
 * excluded. They are stored once, as two particle indices each, and each
 * interaction type has the positions of its pairs in that list, which are
 * all of them unless the force field has no parameters for some, with the
 * pair of particle types of each. The parameters are stored per pair of
 * particle types, see PairTypes. The pairs of an interaction type are in
 * the same order as the TopologyEntry objects of the type.
 */
class NonBondedPairs
{
public:
    /*! \brief Make the pairs of the particles that are not excluded
     * \param[in] natom      The number of particles
     * \param[in] exclusions For each particle the ones it has no pair with
     */
    void make(size_t natom, const std::vector<std::set<size_t>> &exclusions);

    /*! \brief Add an interaction type with all pairs
     * The identifiers of the type pairs are found as for TopologyEntry
     * objects, see Topology::setEntryIdentifiers.
     * \param[in] pd    The force field, or nullptr, then the type pairs have no identifiers
     * \param[in] itype The interaction type
     * \param[in] atoms The particles
     */
    void addType(const ForceField *pd, InteractionType itype, const std::vector<ActAtom> &atoms);

    /*! \brief Fill the parameters of the type pairs of an interaction type
     * As Topology::fillParameters does for TopologyEntry objects: pairs of
     * which the force field has no parameters are removed, and missing
     * parameters of the others are zero.
     * \param[in] msghandler For messages about missing parameters, once per type pair
     * \param[in] pd         The force field, which has the interaction type
     * \param[in] itype      The interaction type
     * \param[in] missing    What to do with missing parameters
     * \return whether the parameters or the pairs changed
     */
    bool fillParameters(MsgHandler       *msghandler,
                        const ForceField *pd,
                        InteractionType   itype,
                        missingParameters missing);

    /*! \brief Put the pairs for which keepFirst is true first
     * The order is kept otherwise, as by std::stable_partition, which is
     * done for the TopologyEntry objects as well.
     * \param[in] keepFirst Which pairs to put first
     * \param[in] shells    Whether this separates the pairs with shells
     *                      (Topology::separateShells), which replaces the
     *                      separation of the intermolecular pairs, or the
     *                      intermolecular ones (Topology::separateIntermolecular)
     */
    void separate(const std::function<bool(const std::array<int, 2> &)> &keepFirst, bool shells);

    /*! \brief Set the parameters of all type pairs of an interaction type
     * For instance for tests. Pairs that were removed since they had no
     * parameters stay removed.
     * \param[in] itype  The interaction type, which has type pairs
     * \param[in] params The parameters, in internal units
     */
    void setParameters(InteractionType itype, const std::vector<double> &params);

    //! \return whether the interaction type has pairs
    bool hasType(InteractionType itype) const { return types_.find(itype) != types_.end(); }

    //! Remove the pairs of an interaction type
    void removeType(InteractionType itype) { types_.erase(itype); }

    //! \return the interaction types that have pairs
    std::vector<InteractionType> types() const;

    //! \return all pairs of an interaction type
    PairRange range(InteractionType itype) const;

    /*! \brief The intermolecular pairs of an interaction type, see Topology::intermolecular
     * \param[in] itype The interaction type
     * \return the pairs, all of them if they were not separated
     */
    PairRange intermolecular(InteractionType itype) const;

    /*! \brief The pairs with or without shells, see Topology::shellEntries
     * \param[in] itype  The interaction type
     * \param[in] shells Whether those with shells or those without
     * \return the pairs
     */
    PairRange shellEntries(InteractionType itype, bool shells) const;

    //! \return the number of pairs of the topology
    size_t npair() const { return pairs_.size(); }

    //! \return the type pairs of an interaction type, nullptr if it has none
    const PairTypes *pairTypes(InteractionType itype) const;

    /*! \brief Leave out the pairs of an interaction type with particles that do not take part
     * \param[in] itype     The interaction type
     * \param[in] particles For each particle whether its pairs are left out
     */
    void leaveOut(InteractionType itype, const std::vector<char> &particles);

    //! \return the memory used, in bytes, approximately
    size_t memory() const;

private:
    //! The pairs of one interaction type
    struct TypePairs
    {
        //! Positions in pairs_ of the pairs of this type, increasing
        std::vector<uint32_t> pos;
        //! The type pair of each, an index in types
        std::vector<uint32_t> type;
        //! The type pairs, without for Coulomb
        PairTypes             types;
        //! Whether the pairs have type pairs, false for Coulomb
        bool                  hasTypes = false;
        //! The first intermolecular pair in pos, see separate
        size_t                intermolecularStart = 0;
        //! The first pair with a shell in pos, see separate
        size_t                shellStart = 0;
        /*! \brief The force field entry of each type pair, nullptr if it has none
         * Found once for a structure of the parameter list, see
         * ForceFieldParameterList::structureVersion, which is unique
         * for all lists, such that fillParameters only copies values.
         */
        std::vector<const ForceFieldParameterMap *> entries;
        //! The structure version of the parameter list that entries are of, 0 if none
        size_t                entriesVersion = 0;
        //! Whether pairs without parameters may still have to be left out, see fillParameters
        bool                  mayRemove = true;
        /*! \brief Whether each type pair is left out since it does not interact and cannot
         * Its epsilon is zero, and so it stays since the epsilons of its
         * particle types are fixed, see fillParameters. Its pairs are
         * removed once, after which droppedApplied is true.
         */
        std::vector<char>     dropped;
        bool                  droppedApplied = false;
    };
    //! The particles of each pair
    std::vector<std::array<int, 2>>        pairs_;
    //! The pairs of each interaction type
    std::map<InteractionType, TypePairs>   types_;
    //! Whether the intermolecular pairs were separated
    bool                                   intermolecular_ = false;
    //! Whether the pairs with shells were separated
    bool                                   shells_         = false;
    //! \return a range of pairs of an interaction type
    PairRange makeRange(const TypePairs &tp, size_t begin, size_t end) const;
};

} // namespace alexandria

#endif
