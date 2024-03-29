/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

#include "act/molprop/molpropobservable.h"

#include <map>
#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "act/molprop/multipole_names.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/units.h"

namespace alexandria
{

std::map<MolPropObservable, const char *> mpo_name_ =
{
    { MolPropObservable::POTENTIAL, "potential" }, 
    { MolPropObservable::FREQUENCY, "frequency" },
    { MolPropObservable::INTENSITY, "intensity" },
    { MolPropObservable::DIPOLE, "dipole" }, 
    { MolPropObservable::QUADRUPOLE, "quadrupole" },
    { MolPropObservable::OCTUPOLE, "octupole" },
    { MolPropObservable::HEXADECAPOLE, "hexadecapole" },
    { MolPropObservable::POLARIZABILITY, "polarizability" }, 
    { MolPropObservable::HF, "HF" }, 
    { MolPropObservable::DELTAE0, "DeltaE0" }, 
    { MolPropObservable::INTERACTIONENERGY, "InteractionEnergy" },
    { MolPropObservable::ELECTROSTATICS, "Electrostatics" },
    { MolPropObservable::INDUCTION, "Induction" },
    { MolPropObservable::EXCHANGE, "Exchange" },
    { MolPropObservable::DISPERSION, "Dispersion" },
    { MolPropObservable::DHFORM, "DeltaHform" }, 
    { MolPropObservable::DGFORM, "DeltaGform" },
    { MolPropObservable::DSFORM, "DeltaSform" },
    { MolPropObservable::ZPE, "ZPE" }, 
    { MolPropObservable::ENTROPY, "S0" },
    { MolPropObservable::STRANS, "Strans" },
    { MolPropObservable::SROT, "Srot" },
    { MolPropObservable::SVIB, "Svib" },
    { MolPropObservable::CP, "cp" },
    { MolPropObservable::CV, "cv" },
    { MolPropObservable::CHARGE, "charge" }
};

std::map<MolPropObservable, const char *> mpo_unit_ =
{
    { MolPropObservable::POTENTIAL,         "kJ/mol e" }, 
    { MolPropObservable::FREQUENCY,         "cm^-1" },
    { MolPropObservable::INTENSITY,         "nm/mol" },
    { MolPropObservable::DIPOLE,            "D" }, 
    { MolPropObservable::QUADRUPOLE,        "B" },
    { MolPropObservable::OCTUPOLE,          "D.Angstrom2" },
    { MolPropObservable::HEXADECAPOLE,      "D.Angstrom3" },
    { MolPropObservable::POLARIZABILITY,    "Angstrom3" }, 
    { MolPropObservable::HF,                "kJ/mol" }, 
    { MolPropObservable::DELTAE0,           "kJ/mol" }, 
    { MolPropObservable::INTERACTIONENERGY, "kJ/mol" },
    { MolPropObservable::ELECTROSTATICS,    "kJ/mol" },
    { MolPropObservable::INDUCTION,         "kJ/mol" },
    { MolPropObservable::EXCHANGE,          "kJ/mol" },
    { MolPropObservable::DISPERSION,        "kJ/mol" },
    { MolPropObservable::DHFORM,            "kJ/mol" }, 
    { MolPropObservable::DGFORM,            "kJ/mol" }, 
    { MolPropObservable::DSFORM,            "J/mol K" }, 
    { MolPropObservable::ZPE,               "kJ/mol" }, 
    { MolPropObservable::ENTROPY,           "J/mol K" },
    { MolPropObservable::STRANS,            "J/mol K" },
    { MolPropObservable::SROT,              "J/mol K" },
    { MolPropObservable::SVIB,              "J/mol K" },
    { MolPropObservable::CP,                "J/mol K" },
    { MolPropObservable::CV,                "J/mol K" },
    { MolPropObservable::CHARGE,            "e" }
};

const char *mpo_name(MolPropObservable MPO)
{
    if (mpo_name_.find(MPO) == mpo_name_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Incorrect MolPropObservable %d", int(MPO)).c_str()));
    }
    
    return mpo_name_[MPO];
}

const char *mpo_unit2(MolPropObservable MPO)
{
    return mpo_unit_[MPO];
}

bool stringToMolPropObservable(const std::string &str, MolPropObservable *mpo)
{
    for (const auto &mn : mpo_name_)
    {
        if (strcasecmp(mn.second, str.c_str()) == 0)
        {
            *mpo = mn.first;
            return true;
        }
    }
    return false;
}

static void print_values(FILE *fp, const std::vector<double> values, double factor)
{
    fprintf(fp, "  Values:");
    for(const auto v : values)
    {
        fprintf(fp, " %g", v*factor);
    }
    fprintf(fp, "\n");
}

void GenericProperty::Dump(FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "  Type: %s Unit: %s T: %g K Phase: %s\n",
            getType(), inputUnit_.c_str(), getTemperature(),
            phase2string(getPhase()).c_str());
}

CommunicationStatus GenericProperty::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;

    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        std::string mpo_type(mpo_name(mpo_));
        cr->send_str(dest, &mpo_type);
        cr->send_str(dest, &type_);
        cr->send_str(dest, &inputUnit_);
        cr->send_str(dest, &unit_);
        cr->send_double(dest, T_);
        cr->send_int(dest, (int) eP_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send GenericProperty, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus GenericProperty::BroadCast(const CommunicationRecord *cr,
                                               gmx_unused int             root,
                                               MPI_Comm                   comm)
{
    CommunicationStatus cs = cr->bcast_data(comm);

    if (CommunicationStatus::OK == cs)
    {
        const char *mn = mpo_name(mpo_);
        std::string mpo_type(mn);
        cr->bcast(&mpo_type, comm);
        if (!stringToMolPropObservable(mpo_type, &mpo_))
        {
            gmx_fatal(FARGS, "Unknown observable %s", mpo_type.c_str());
        }
        cr->bcast(&type_, comm);
        cr->bcast(&inputUnit_, comm);
        cr->bcast(&unit_, comm);
        cr->bcast(&T_, comm);
        int ep = static_cast<int>(eP_);
        cr->bcast(&ep, comm);
        eP_ = static_cast<ePhase>(ep);
        
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive GenericProperty, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus GenericProperty::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;

    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        std::string mpo_type;
        cr->recv_str(src, &mpo_type);
        if (!stringToMolPropObservable(mpo_type, &mpo_))
        {
            gmx_fatal(FARGS, "Unknown observable %s", mpo_type.c_str());
        }
        cr->recv_str(src, &type_);
        cr->recv_str(src, &inputUnit_);
        cr->recv_str(src, &unit_);
        T_   = cr->recv_double(src);
        eP_  = (ePhase) cr->recv_int(src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive GenericProperty, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

MolecularMultipole::MolecularMultipole(const std::string &type,
                                       const std::string &inputUnit,
                                       double             T,
                                       MolPropObservable  mpo) :
    GenericProperty(mpo, type, inputUnit, T, ePhase::GAS)
{
    size_t nvalues = multipoleNames(mpo).size();
    values_.resize(nvalues, 0.0);
    setUnit(gromacsUnit(inputUnit));
}

bool MolecularMultipole::hasId(const std::string &myid)
{
    int index;
    
    return (multipoleIndex(myid, &index) ||
            myid.compare("average") == 0 ||
            myid.compare("error") == 0);
}

void MolecularMultipole::Dump(FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    GenericProperty::Dump(fp);
    fprintf(fp, "  Average %g Error %g\n", average_, error_);
    double factor = convertFromGromacs(1, getInputUnit());
    print_values(fp, values_, factor);
}

void MolecularMultipole::setValue(const std::string &myid, double value)
{
    int index; 
    
    double actValue = convertToGromacs(value, getInputUnit());
    if (multipoleIndex(myid, &index))
    {
        range_check(index, 0, values_.size());
        values_[index] = actValue;
    }
    else if (myid.compare("average") == 0)
    {
        average_ = actValue;
    }
    else if (myid.compare("error") == 0)
    {
        error_ = actValue;
    }
    else if (debug)
    {
        fprintf(debug, "Unknown id %s with value %g\n", myid.c_str(), value);
    }
}

CommunicationStatus MolecularMultipole::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_int(dest, values_.size());
        for(auto &m : values_)
        {
            cr->send_double(dest, m);
        }
        cr->send_double(dest, average_);
        cr->send_double(dest, error_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send MolecularMultipole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularMultipole::BroadCast(const CommunicationRecord *cr,
                                                  int                        root,
                                                  MPI_Comm                   comm)
{
    CommunicationStatus cs = GenericProperty::BroadCast(cr, root, comm);

    if (CommunicationStatus::OK == cs)
    {
        int nvalue = values_.size();
        cr->bcast(&nvalue, comm);
        if (cr->rank() != root)
        {
            values_.resize(nvalue);
        }
        cr->bcast(&values_, comm);
        cr->bcast(&average_, comm);
        cr->bcast(&error_, comm);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive MolecularMultipole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularMultipole::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        values_.clear();
        size_t N = cr->recv_int(src);
        for(size_t m = 0; m < N; m++)
        {
            double x = cr->recv_double(src);
            values_.push_back(x);
        }
        average_ = cr->recv_double(src);
        error_ = cr->recv_double(src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive MolecularMultipole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

Harmonics::Harmonics(const std::string &inputUnit,
                     double             T,
                     MolPropObservable  mpo) :
    GenericProperty(mpo, "electronic", inputUnit, T, ePhase::GAS)
{
    setUnit(gromacsUnit(inputUnit));
    if (MolPropObservable::FREQUENCY != mpo &&
        MolPropObservable::INTENSITY != mpo)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Passed incorrect observable %s to Harmonics constructor", mpo_name(mpo)).c_str()));
    }
}

void Harmonics::addValue(double value)
{
    values_.push_back(convertToGromacs(value, getInputUnit()));
}

void Harmonics::Dump(FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    GenericProperty::Dump(fp);
    double factor = convertFromGromacs(1, getInputUnit());
    print_values(fp, values_, factor);
}

CommunicationStatus Harmonics::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_int(dest, values_.size());
        for(auto &m : values_)
        {
            cr->send_double(dest, m);
        }
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send Harmonics, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Harmonics::BroadCast(const CommunicationRecord *cr,
                                         int                        root,
                                         MPI_Comm                   comm)
{
    CommunicationStatus cs = GenericProperty::BroadCast(cr, root, comm);

    if (CommunicationStatus::OK == cs)
    {
        int nvalue = values_.size();
        cr->bcast(&nvalue, comm);
        if (cr->rank() != root)
        {
            values_.clear();
        }
        cr->bcast(&values_, comm);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive MolecularMultipole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Harmonics::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        values_.clear();
        size_t N = cr->recv_int(src);
        for(size_t m = 0; m < N; m++)
        {
            double x = cr->recv_double(src);
            values_.push_back(x);
        }
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive MolecularMultipole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

MolecularPolarizability::MolecularPolarizability(const std::string &type,
                                                 const std::string &inputUnit,
                                                 double T,
                                                 double xx, double yy, double zz,
                                                 double xy, double xz, double yz,
                                                 double average, double error) :
    GenericProperty(MolPropObservable::POLARIZABILITY,  type, inputUnit, T, ePhase::GAS)
{ 
    Set(xx, yy, zz, xy, xz, yz);
    average_ = convertToGromacs(average, inputUnit);
    error_   = convertToGromacs(error, inputUnit);
    setUnit(gromacsUnit(inputUnit));
}

double MolecularPolarizability::getValue() const
{
    if (average_ != 0)
    {
        return average_;
    }
    else
    {
        auto alpha = getTensor();
        return (alpha[XX][XX] + alpha[YY][YY] + alpha[ZZ][ZZ])/3.0;
    }
}

void MolecularPolarizability::Set(double xx, double yy, double zz, double xy, double xz, double yz)
{
    auto unit = getInputUnit();
    alpha_[XX][XX] = convertToGromacs(xx, unit);
    alpha_[YY][YY] = convertToGromacs(yy, unit);
    alpha_[ZZ][ZZ] = convertToGromacs(zz, unit);
    alpha_[XX][YY] = convertToGromacs(xy, unit);
    alpha_[XX][ZZ] = convertToGromacs(xz, unit); 
    alpha_[YY][ZZ] = convertToGromacs(yz, unit);
};

void MolecularPolarizability::Dump(FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    GenericProperty::Dump(fp);
    double factor = convertFromGromacs(1, getInputUnit());
    fprintf(fp, "  Average %g Error %g\n", average_*factor, error_*factor);
    fprintf(fp, "  Alpha:\n");
    for(int m = 0; m < DIM; m++)
    {
        fprintf(fp, "  ");
        for(int n = 0; n < DIM; n++)
        {
            fprintf(fp, " %10g", alpha_[m][n]*factor);
        }
        fprintf(fp, "\n");
    }
}

CommunicationStatus MolecularPolarizability::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        for(int m = 0; m < DIM; m++)
        {
            for(int n = 0; n < DIM; n++)
            {
                cr->send_double(dest, alpha_[m][n]);
            }
        }
    }
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_double(dest, average_);
        cr->send_double(dest, error_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send Polarizability, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularPolarizability::BroadCast(const CommunicationRecord *cr,
                                                       int                        root,
                                                       MPI_Comm                   comm)
{
    CommunicationStatus cs = GenericProperty::BroadCast(cr, root, comm);

    if (CommunicationStatus::OK == cs)
    {
        for(int m = 0; m < DIM; m++)
        {
            for(int n = 0; n < DIM; n++)
            {
                cr->bcast(&alpha_[m][n], comm);
            }
        }
        cr->bcast(&average_, comm);
        cr->bcast(&error_, comm);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to broadcast Polarizability, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularPolarizability::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        for(int m = 0; m < DIM; m++)
        {
            for(int n = 0; n < DIM; n++)
            {
                alpha_[m][n] = cr->recv_double(src);
            }
        }
    }
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        average_ = cr->recv_double(src);
        error_   = cr->recv_double(src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to broadcast Polarizability, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

MolecularEnergy::MolecularEnergy(MolPropObservable mpo,
                                 const std::string &type,
                                 const std::string &inputUnit,
                                 double T,
                                 ePhase ep,
                                 double average,
                                 double error)
    : GenericProperty(mpo, type, inputUnit, T, ep)
{ 
    Set(convertToGromacs(average, inputUnit), convertToGromacs(error, inputUnit));
    setUnit(gromacsUnit(inputUnit));
}

void MolecularEnergy::Dump(FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    GenericProperty::Dump(fp);
    fprintf(fp, "  Average %g Error %g\n", average_, error_);
}

CommunicationStatus MolecularEnergy::BroadCast(const CommunicationRecord *cr,
                                               int                        root,
                                               MPI_Comm                   comm)
{
    CommunicationStatus cs = GenericProperty::BroadCast(cr, root, comm);

    if (CommunicationStatus::OK == cs)
    {
        cr->bcast(&average_, comm);
        cr->bcast(&error_, comm);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive MolecularEnergy, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularEnergy::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        average_ = cr->recv_double(src);
        error_   = cr->recv_double(src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive MolecularEnergy, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularEnergy::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CommunicationStatus::OK == cs &&
        CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_double(dest, average_);
        cr->send_double(dest, error_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send MolecularEnergy, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

ElectrostaticPotential::ElectrostaticPotential(const std::string &xyzInputUnit,
                                               const std::string &vInputUnit)
{
    xyzInputUnit_ = xyzInputUnit;
    vInputUnit_   = vInputUnit;
    xyzUnit_      = gromacsUnit(xyzInputUnit);
    vUnit_        = gromacsUnit(vInputUnit);
}

void ElectrostaticPotential::addPoint(int    espid, 
                                      double x,
                                      double y,
                                      double z,
                                      double V)
{
    double xfac = convertToGromacs(1.0, xyzInputUnit_);
    double Vfac = convertToGromacs(1.0, vInputUnit_);
    espID_.push_back(espid);
    gmx::RVec xyz = { x*xfac, y*xfac, z*xfac };
    xyz_.push_back(xyz);
    V_.push_back(V*Vfac);
}

void ElectrostaticPotential::Dump(FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    GenericProperty::Dump(fp);
    double factor = convertFromGromacs(1, getInputUnit());
    fprintf(fp, "  %4s  %8s  %8s  %8s  %12s\n", "ID", "x", "y", "z", "V");
    for(size_t i = 0; i < espID_.size(); i++)
    {
        fprintf(fp, "  %4d  %8.3f  %8.3f  %8.3f  %12.3f\n", espID_[i],xyz_[i][XX],
                xyz_[i][YY], xyz_[i][ZZ], factor*V_[i]);
    }
}

CommunicationStatus ElectrostaticPotential::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;

    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv_str(src, &xyzInputUnit_);
        cr->recv_str(src, &vInputUnit_);
        cr->recv_str(src, &xyzUnit_);
        cr->recv_str(src, &vUnit_);
        int nesp = cr->recv_int(src);
        espID_.clear();
        xyz_.clear();
        V_.clear();
        for(int i = 0; i < nesp; i++)
        {
            espID_.push_back(cr->recv_int(src));
            std::vector<double> xyz;
            cr->recv_double_vector(src, &xyz);
            gmx::RVec xxx = { xyz[XX], xyz[YY], xyz[ZZ] };
            xyz_.push_back(xxx);
            V_.push_back(cr->recv_double(src));
        }
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive ElectrostaticPotential, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus ElectrostaticPotential::BroadCast(const CommunicationRecord *cr,
                                                      int                        root,
                                                      MPI_Comm                   comm)
{
    CommunicationStatus cs = GenericProperty::BroadCast(cr, root, comm);

    if (CommunicationStatus::OK == cs)
    {
        cr->bcast(&xyzInputUnit_, comm);
        cr->bcast(&vInputUnit_, comm);
        cr->bcast(&xyzUnit_, comm);
        cr->bcast(&vUnit_, comm);
        int nesp = espID_.size();
        cr->bcast(&nesp, comm);
        if (cr->rank() != root)
        {
            espID_.resize(nesp);
            xyz_.resize(nesp);
            V_.resize(nesp);
        }

        for(size_t i = 0; i < espID_.size(); i++)
        {
            cr->bcast(&espID_[i], comm);
            cr->bcast(&xyz_[i][XX], comm);
            cr->bcast(&xyz_[i][YY], comm);
            cr->bcast(&xyz_[i][ZZ], comm);
            cr->bcast(&V_[i], comm);
        }
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive ElectrostaticPotential, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus ElectrostaticPotential::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;

    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_str(dest, &xyzInputUnit_);
        cr->send_str(dest, &vInputUnit_);
        cr->send_str(dest, &xyzUnit_);
        cr->send_str(dest, &vUnit_);
        cr->send_int(dest, espID_.size());
        for(size_t i = 0; i < espID_.size(); i++)
        {
            cr->send_int(dest, espID_[i]);
            std::vector<double> xyz = { xyz_[i][XX], xyz_[i][YY], xyz_[i][ZZ] };
            cr->send_double_vector(dest, &xyz);
            cr->send_double(dest, V_[i]);
        }
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send ElectrostaticPotential, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

} // namespace alexandria
