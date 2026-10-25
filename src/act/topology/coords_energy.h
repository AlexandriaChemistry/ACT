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

namespace alexandria
{
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
    // Empty Constructor
    ACTEnergy() {}
    /*! \brief Constructor
     * \param[in] experiment The id
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
     */
    ACTEnergy(int experiment, const std::string &datafile) : experiment_(experiment), datafile_(datafile) {}

    //! \return The experiment ID
    int id() const { return experiment_; }

    //! \return The data file name
    const std::string &datafile() const { return datafile_; }

    //! Is there a QM energy
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

    //! set QM value
    void setQM(double eqm)
    {
        eqm_ = eqm;
        haveQM_ = true;
    }

    //! Is there an ACT energy
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

    //! set ACT value
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
    //! \return backup energy
    double eactBefore() const { return eactBefore_; }
};

//! \brief Map to match interaction type to energies
typedef std::map<InteractionType, ACTEnergy> ACTEnergyMap;
//! \brief Vector of the above
typedef std::vector<ACTEnergyMap>            ACTEnergyMapVector;

typedef struct
{
    //! Indicates that the shells have converged in this structure
    bool                   done;
    //! Pointer to datafile copied from Experiment
    std::string            dataFile;
    //! Coordinates
    std::vector<gmx::RVec> x;
    //! Relaxed coordinates
    std::vector<gmx::RVec> xrelax;
    //! Computed forces
    std::vector<gmx::RVec> fcalc;
    //! Observed forces
    std::vector<gmx::RVec> fobs;
    //! Energy map comparing QM to ACT energies
    ACTEnergyMap           eMap;
    //! Lowest energy for use in Boltzmann weighting
    double                 eMin;
} CoordsEnergy;

} // namespace

#endif
