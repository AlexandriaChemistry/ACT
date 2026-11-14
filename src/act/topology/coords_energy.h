/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_COORDS_ENERGY_H
#define ACT_COORDS_ENERGY_H

#include <concepts>
#include <map>
#include <string>
#include <vector>

#include "act/basics/coordinate.h"
#include "act/basics/interactiontype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

/*! \brief Vectors of the particles of one conformation, such as its
 * coordinates or forces, without owning them
 * Component m of particle a is at data[(DIM*a+m)*stride], such that the
 * vectors can be those of one conformation in arrays of many (see
 * ConformationArrays) or contiguous, with a stride of one. Particle a
 * is indexed as v[a][m], like a std::vector<gmx::RVec>.
 * \tparam T The type of the components, const for read-only access
 */
template <typename T>
class ParticleVectors
{
public:
    //! The vector of one particle
    class Vector
    {
    public:
        /*! \brief Constructor
         * \param[in] p      The first component
         * \param[in] stride Distance between the components
         */
        Vector(T *p, size_t stride) : p_(p), stride_(stride) {}
        /*! \brief Return a component
         * \param[in] m The dimension
         * \return the component
         */
        T &operator[](int m) const { return p_[m*stride_]; }
        //! \return the vector, for vector functions
        gmx::RVec rvec() const { return { static_cast<real>(p_[0]), static_cast<real>(p_[stride_]), static_cast<real>(p_[2*stride_]) }; }
    private:
        //! The first component
        T     *p_;
        //! Distance between the components
        size_t stride_;
    };

    //! Constructor for no particles
    ParticleVectors() = default;

    /*! \brief Constructor
     * \param[in] data   The first component of the first particle
     * \param[in] stride Distance between components
     * \param[in] natom  Number of particles
     */
    ParticleVectors(T *data, size_t stride, size_t natom) : data_(data), stride_(stride), natom_(natom) {}

    //! Read-only access to the same vectors
    operator ParticleVectors<const T>() const { return { data_, stride_, natom_ }; }

    //! \return the number of particles
    size_t size() const { return natom_; }

    //! \return whether there are no particles
    bool empty() const { return natom_ == 0; }

    //! \return the distance between components
    size_t stride() const { return stride_; }

    //! \return the first component of the first particle
    T *data() const { return data_; }

    /*! \brief Return the vector of a particle
     * \param[in] a The particle
     * \return the vector
     */
    Vector operator[](size_t a) const { return Vector(data_ + DIM*a*stride_, stride_); }

    /*! \brief Set the vectors of all particles
     * \param[in] v The vectors, as many as there are particles
     */
    void assign(const std::vector<gmx::RVec> &v) const
    {
        GMX_RELEASE_ASSERT(v.size() == natom_, "Wrong number of particles");
        for(size_t a = 0; a < natom_; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                (*this)[a][m] = v[a][m];
            }
        }
    }

    //! \return the vectors of all particles
    std::vector<gmx::RVec> toRVec() const
    {
        std::vector<gmx::RVec> v(natom_);
        for(size_t a = 0; a < natom_; a++)
        {
            v[a] = (*this)[a].rvec();
        }
        return v;
    }

private:
    //! The first component of the first particle
    T     *data_   = nullptr;
    //! Distance between components
    size_t stride_ = 1;
    //! Number of particles
    size_t natom_  = 0;
};

/*! \brief Integer coordinates of the particles of one conformation, without owning them
 * As ParticleVectors, with the components stored as integers in units of
 * coordinateUnit, see coordinate.h. A component reads as a value in nm,
 * a value in nm assigned to it is rounded to the nearest unit, it has to be
 * within the range, see toCoordinate. Assigning
 * one component to another copies the integer, which is exact, and the
 * difference of two vectors is computed from the integers.
 * \tparam C Coordinate, or const Coordinate for read-only access
 */
template <typename C>
class CoordinateVectors
{
public:
    //! One component
    class Component
    {
    public:
        /*! \brief Constructor
         * \param[in] p The component
         */
        explicit Component(C *p) : p_(p) {}
        //! Copying a component refers to the same component
        Component(const Component &) = default;
        //! \return the component in nm
        operator double() const { return fromCoordinate(*p_); }
        //! \return the component in units of coordinateUnit
        Coordinate units() const { return *p_; }
        /*! \brief Set the component
         * \param[in] v The value in nm
         * \return this component
         */
        const Component &operator=(double v) const { *p_ = toCoordinate(v); return *this; }
        /*! \brief Copy the integer of another component
         * \param[in] o The other component
         * \return this component
         */
        const Component &operator=(const Component &o) const { *p_ = *o.p_; return *this; }
        //! \copydoc operator=(const Component &) const
        template <typename O>
        requires requires(const O &o) { { o.units() } -> std::same_as<Coordinate>; }
        const Component &operator=(const O &o) const { *p_ = o.units(); return *this; }
        /*! \brief Move the component
         * A move beyond the range of the coordinates ends at the range, see
         * moveCoordinate, since particles such as shells may run away.
         * \param[in] v The displacement in nm
         * \return this component
         */
        const Component &operator+=(double v) const { *p_ = moveCoordinate(*p_, v); return *this; }
        /*! \brief Move the component the other way, see operator+=
         * \param[in] v The displacement in nm
         * \return this component
         */
        const Component &operator-=(double v) const { *p_ = moveCoordinate(*p_, -v); return *this; }
    private:
        //! The component
        C *p_;
    };

    //! The vector of one particle
    class Vector
    {
    public:
        /*! \brief Constructor
         * \param[in] p      The first component
         * \param[in] stride Distance between the components
         */
        Vector(C *p, size_t stride) : p_(p), stride_(stride) {}
        /*! \brief Return a component
         * \param[in] m The dimension
         * \return the component
         */
        Component operator[](int m) const { return Component(p_ + m*stride_); }
        //! \return the vector in nm, for vector functions
        gmx::RVec rvec() const
        {
            return { fromCoordinateAs<real>(p_[0]), fromCoordinateAs<real>(p_[stride_]),
                     fromCoordinateAs<real>(p_[2*stride_]) };
        }
        /*! \brief The difference with another vector, from the integers
         * \param[in] o The other vector
         * \return this vector minus o in nm
         */
        template <typename O>
        requires requires(const O &o) { { o.units(0) } -> std::same_as<Coordinate>; }
        gmx::RVec operator-(const O &o) const
        {
            gmx::RVec d;
            for(int m = 0; m < DIM; m++)
            {
                d[m] = fromCoordinateAs<real>(coordinateDifference(p_[m*stride_], o.units(m)));
            }
            return d;
        }
        /*! \return a component in units of coordinateUnit
         * \param[in] m The dimension
         */
        Coordinate units(int m) const { return p_[m*stride_]; }
    private:
        //! The first component
        C     *p_;
        //! Distance between the components
        size_t stride_;
    };

    //! Constructor for no particles
    CoordinateVectors() = default;

    /*! \brief Constructor
     * \param[in] data   The first component of the first particle
     * \param[in] stride Distance between components
     * \param[in] natom  Number of particles
     */
    CoordinateVectors(C *data, size_t stride, size_t natom) : data_(data), stride_(stride), natom_(natom) {}

    //! Read-only access to the same vectors
    operator CoordinateVectors<const C>() const { return { data_, stride_, natom_ }; }

    //! \return the number of particles
    size_t size() const { return natom_; }

    //! \return whether there are no particles
    bool empty() const { return natom_ == 0; }

    //! \return the distance between components
    size_t stride() const { return stride_; }

    //! \return the first component of the first particle
    C *data() const { return data_; }

    /*! \brief Return the vector of a particle
     * \param[in] a The particle
     * \return the vector
     */
    Vector operator[](size_t a) const { return Vector(data_ + DIM*a*stride_, stride_); }

    /*! \brief Set the vectors of all particles
     * \param[in] v The vectors in nm, as many as there are particles
     */
    void assign(const std::vector<gmx::RVec> &v) const
    {
        GMX_RELEASE_ASSERT(v.size() == natom_, "Wrong number of particles");
        for(size_t a = 0; a < natom_; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                (*this)[a][m] = v[a][m];
            }
        }
    }

    //! \return the vectors of all particles in nm
    std::vector<gmx::RVec> toRVec() const
    {
        std::vector<gmx::RVec> v(natom_);
        for(size_t a = 0; a < natom_; a++)
        {
            v[a] = (*this)[a].rvec();
        }
        return v;
    }

private:
    //! The first component of the first particle
    C     *data_   = nullptr;
    //! Distance between components
    size_t stride_ = 1;
    //! Number of particles
    size_t natom_  = 0;
};

/*! \brief Integer coordinates of the particles of one structure, owning them
 * The components are stored contiguously, as CoordinateVectors with a
 * stride of one, in units of coordinateUnit, see coordinate.h. Particle a
 * is indexed as x[a][m], like a std::vector<gmx::RVec>, a component reads
 * as a value in nm and the difference of two vectors is exact.
 */
class Coordinates
{
public:
    //! Constructor for no particles
    Coordinates() = default;

    /*! \brief Constructor for particles at the origin
     * \param[in] natom The number of particles
     */
    explicit Coordinates(size_t natom) : x_(DIM*natom, 0) {}

    /*! \brief Constructor from coordinates in nm
     * \param[in] v The coordinates, rounded to the unit
     */
    explicit Coordinates(const std::vector<gmx::RVec> &v) : x_(DIM*v.size(), 0) { view().assign(v); }

    /*! \brief Constructor copying the integers of a view
     * \param[in] v The coordinates
     */
    explicit Coordinates(const CoordinateVectors<const Coordinate> &v) : x_(DIM*v.size(), 0)
    {
        for(size_t a = 0; a < v.size(); a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                x_[DIM*a+m] = v[a].units(m);
            }
        }
    }

    //! \return the number of particles
    size_t size() const { return x_.size()/DIM; }

    //! \return whether there are no particles
    bool empty() const { return x_.empty(); }

    //! \return a view of the coordinates
    CoordinateVectors<Coordinate> view() { return { x_.data(), 1, size() }; }

    //! \return a read-only view of the coordinates
    CoordinateVectors<const Coordinate> view() const { return { x_.data(), 1, size() }; }

    /*! \brief Return the vector of a particle
     * \param[in] a The particle
     * \return the vector
     */
    CoordinateVectors<Coordinate>::Vector operator[](size_t a) { return view()[a]; }

    //! \copydoc operator[](size_t)
    CoordinateVectors<const Coordinate>::Vector operator[](size_t a) const { return view()[a]; }

    /*! \brief Set the coordinates of all particles
     * \param[in] v The coordinates in nm, the number of particles changes to theirs
     */
    void assign(const std::vector<gmx::RVec> &v)
    {
        x_.resize(DIM*v.size());
        view().assign(v);
    }

    //! \return the coordinates of all particles in nm
    std::vector<gmx::RVec> toRVec() const { return view().toRVec(); }

private:
    //! The components, that of particle a in dimension m at DIM*a+m
    std::vector<Coordinate> x_;
};

/*! \brief Class to hold corresponding QM (reference) and ACT energies
 */
class ACTEnergy
{
private:
    //! Id of the experiment
    int    experiment_ = 0;
    //! Copy of the datafile name
    std::string datafile_;
    //! QM energy
    double eqm_        = 0;
    //! ACT energy
    double eact_       = 0;
    //! ACT energy backup, for electrostatics calculations
    double eactBefore_ = 0;
    //! QM present
    bool haveQM_       = false;
    //! ACT present
    bool haveACT_      = false;
public:
    //! Empty constructor
    ACTEnergy() = default;
    /*! \brief Constructor
     * \param[in] experiment The id
     * \param[in] datafile   The data file name
     * \param[in] eqm        The QM energy
     * \param[in] eact       The ACT energy
     */
    ACTEnergy(int experiment, const std::string &datafile, double eqm, double eact) : experiment_(experiment), datafile_(datafile)
    {
        setQM(eqm);
        setACT(eact);
    }
    /*! \brief Constructor
     * \param[in] experiment The id
     * \param[in] datafile   The data file name
     */
    ACTEnergy(int experiment, const std::string &datafile) : experiment_(experiment), datafile_(datafile) {}

    //! \return The experiment ID
    int id() const { return experiment_; }

    //! \return The data file name
    const std::string &datafile() const { return datafile_; }

    //! \return whether there is a QM energy
    bool haveQM() const { return haveQM_; }

    //! \return The QM energy    
    double eqm() const
    {
        if (!haveQM_)
        {
            GMX_THROW(gmx::InvalidInputError("Trying to extract a QM energy when there is none"));
        }
        return eqm_;
    }

    /*! \brief Set the QM energy and mark it as present
     * \param[in] eqm The QM energy
     */
    void setQM(double eqm)
    {
        eqm_ = eqm;
        haveQM_ = true;
    }

    //! \return whether there is an ACT energy
    bool haveACT() const { return haveACT_; }

    //! \return The ACT energy
    double eact() const
    {
        if (!haveACT_)
        {
            GMX_THROW(gmx::InternalError("Trying to extract an ACT energy when there is none"));
        }
        return eact_;
    }

    /*! \brief Set the ACT energy and mark it as present
     * \param[in] eact The ACT energy
     */
    void setACT(double eact)
    {
        eact_ = eact;
        haveACT_ = true;
    }
    //! Backup ACT energy
    void backup()
    {
        eactBefore_ = eact_;
    }
    //! Reset backup ACT energy
    void restoreBackup()
    {
        eact_       = eactBefore_;
        eactBefore_ = 0;
    }
};

//! \brief Map to match interaction type to energies
typedef std::map<InteractionType, ACTEnergy> ACTEnergyMap;
//! \brief Vector of the above
typedef std::vector<ACTEnergyMap>            ACTEnergyMapVector;

//! \brief Coordinates, forces and energies of one structure of a molecule
typedef struct
{
    //! Indicates that the shells have converged in this structure
    bool                   done;
    //! Pointer to datafile copied from Experiment
    std::string            dataFile;
    //! Relaxed coordinates, the coordinates are those of ConformationData::x
    std::vector<gmx::RVec> xrelax;
    //! The shells of the previous relaxation, three components per shell, see ForceComputer::setKeepShells
    std::vector<Coordinate> xkept;
    //! Observed forces
    std::vector<gmx::RVec> fobs;
    //! The QM energies, the ACT energies are in ConformationData::energies
    ACTEnergyMap           eMap;
    //! Lowest energy for use in Boltzmann weighting
    double                 eMin;
} CoordsEnergy;

} // namespace

#endif
