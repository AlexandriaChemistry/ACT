/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
    { MolPropObservable::DIPOLE, "dipole" }, 
    { MolPropObservable::QUADRUPOLE, "quadrupole" },
    { MolPropObservable::OCTUPOLE, "octupole" },
    { MolPropObservable::HEXADECAPOLE, "hexadecapole" },
    { MolPropObservable::POLARIZABILITY, "polarizability" }, 
    { MolPropObservable::HF, "HF" }, 
    { MolPropObservable::DELTAE0, "DeltaE0" }, 
    { MolPropObservable::DHFORM, "DeltaHform" }, 
    { MolPropObservable::DGFORM, "DeltaGform" },
    { MolPropObservable::DSFORM, "DeltaSform" },
    { MolPropObservable::ZPE, "ZPE" }, 
    { MolPropObservable::EMOL, "emol" }, 
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
    { MolPropObservable::POTENTIAL,      "e/nm" }, 
    { MolPropObservable::DIPOLE,         "D" }, 
    { MolPropObservable::QUADRUPOLE,     "B" },
    { MolPropObservable::OCTUPOLE,       "D.Angstrom2" },
    { MolPropObservable::HEXADECAPOLE,   "D.Angstrom3" },
    { MolPropObservable::POLARIZABILITY, "Angstrom3" }, 
    { MolPropObservable::HF,             "kJ/mol" }, 
    { MolPropObservable::DELTAE0,        "kJ/mol" }, 
    { MolPropObservable::DHFORM,         "kJ/mol" }, 
    { MolPropObservable::DGFORM,         "kJ/mol" }, 
    { MolPropObservable::DSFORM,         "J/mol K" }, 
    { MolPropObservable::ZPE,            "kJ/mol" }, 
    { MolPropObservable::EMOL,           "kJ/mol" }, 
    { MolPropObservable::ENTROPY,        "J/mol K" },
    { MolPropObservable::STRANS,         "J/mol K" },
    { MolPropObservable::SROT,           "J/mol K" },
    { MolPropObservable::SVIB,           "J/mol K" },
    { MolPropObservable::CP,             "J/mol K" },
    { MolPropObservable::CV,             "J/mol K" },
    { MolPropObservable::CHARGE,         "e" }
};

const char *mpo_name(MolPropObservable MPO)
{
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
        fprintf(debug, "Trying to received Polarizability, status %s\n", cs_name(cs));
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

CommunicationStatus ElectrostaticPotential::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;

    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv_str(src, &xyzInputUnit_);
        cr->recv_str(src, &vInputUnit_);
        cr->recv_str(src, &xyzUnit_);
        cr->recv_str(src, &vUnit_);
        espID_ = cr->recv_int(src);
        x_     = cr->recv_double(src);
        y_     = cr->recv_double(src);
        z_     = cr->recv_double(src);
        V_     = cr->recv_double(src);
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
        cr->send_int(dest, espID_);
        cr->send_double(dest, x_);
        cr->send_double(dest, y_);
        cr->send_double(dest, z_);
        cr->send_double(dest, V_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send ElectrostaticPotential, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

void ElectrostaticPotential::set(const std::string &xyzInputUnit,
                                 const std::string &vInputUnit,
                                 int                espid, 
                                 double             x,
                                 double             y,
                                 double             z,
                                 double             V)
{ 
    xyzInputUnit_ = xyzInputUnit;
    vInputUnit_   = vInputUnit;
    espID_ = espid;
    x_ = convertToGromacs(x, xyzInputUnit);
    y_ = convertToGromacs(y, xyzInputUnit);
    z_ = convertToGromacs(z, xyzInputUnit);
    V_ = convertToGromacs(V, vInputUnit); 
    xyzUnit_.assign("nm");
    vUnit_.assign("kJ/mol e");
}

void ElectrostaticPotential::get(std::string *xyz_unit,
                                 std::string *V_unit,
                                 int         *espid,
                                 double      *x,
                                 double      *y,
                                 double      *z,
                                 double *V) const
{
    xyz_unit->assign(xyzUnit_);
    V_unit->assign(vUnit_);
    *espid = espID_;
    *x = x_;
    *y = y_;
    *z = z_;
    *V = V_;
}

} // namespace alexandria
