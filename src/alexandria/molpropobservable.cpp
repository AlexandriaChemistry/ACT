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

#include "molpropobservable.h"

#include <map>
#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "communicationrecord.h"

namespace alexandria
{

std::map<MolPropObservable, const char *> mpo_name_ =
{
    { MolPropObservable::POTENTIAL, "potential" }, 
    { MolPropObservable::DIPOLE, "dipole" }, 
    { MolPropObservable::QUADRUPOLE, "quadrupole" },
    { MolPropObservable::POLARIZABILITY, "polarizability" }, 
    { MolPropObservable::HF, "HF" }, 
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
    { MolPropObservable::POTENTIAL, "e/nm" }, 
    { MolPropObservable::DIPOLE, "D" }, 
    { MolPropObservable::QUADRUPOLE, "B" },
    { MolPropObservable::POLARIZABILITY, "\\AA$^3$" }, 
    { MolPropObservable::HF, "kJ/mol" }, 
    { MolPropObservable::DHFORM, "kJ/mol" }, 
    { MolPropObservable::DGFORM, "kJ/mol" }, 
    { MolPropObservable::DSFORM, "J/mol K" }, 
    { MolPropObservable::ZPE, "kJ/mol" }, 
    { MolPropObservable::EMOL, "kJ/mol" }, 
    { MolPropObservable::ENTROPY, "J/mol K" },
    { MolPropObservable::STRANS, "J/mol K" },
    { MolPropObservable::SROT, "J/mol K" },
    { MolPropObservable::SVIB, "J/mol K" },
    { MolPropObservable::CP, "J/mol K" },
    { MolPropObservable::CV, "J/mol K" },
    { MolPropObservable::CHARGE, "e" }
};

const char *mpo_name(MolPropObservable MPO)
{
    return mpo_name_[MPO];
}

const char *mpo_unit(MolPropObservable MPO)
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
    CommunicationStatus cs;

    cs = cr->send_data(dest);
    if (CS_OK == cs)
    {
        std::string type(mpo_name(mpo_));
        cr->send_str(dest, &type);
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
    CommunicationStatus cs;

    cs = cr->recv_data(src);
    if (CS_OK == cs)
    {
        std::string type;
        cr->recv_str(src, &type);
        if (!stringToMolPropObservable(type, &mpo_))
        {
            gmx_fatal(FARGS, "Unknown observable %s", type.c_str());
        }
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

CommunicationStatus MolecularDipole::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = cr->send_data(dest);
    }
    if (CS_OK == cs)
    {
        cr->send_int(dest, mu_.size());
        for(auto &m : mu_)
        {
            cr->send_double(dest, m);
        }
        cr->send_double(dest, average_);
        cr->send_double(dest, error_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send MolecularDipole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularDipole::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = cr->recv_data(src);
    }
    if (CS_OK == cs)
    {
        mu_.clear();
        size_t N = cr->recv_int(src);
        for(size_t m = 0; m < N; m++)
        {
            double x = cr->recv_double(src);
            mu_.push_back(x);
        }
        average_  = cr->recv_double(src);
        error_    = cr->recv_double(src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive MolecularDipole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularQuadrupole::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = cr->send_data(dest);
    }
    if (CS_OK == cs)
    {
        for(int m = 0; m < DIM; m++)
        {
            for(int n = 0; n < DIM; n++)
            {
                cr->send_double(dest, quad_[m][n]);
            }
        }
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send MolecularQuadrupole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularQuadrupole::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = cr->recv_data(src);
    }
    if (CS_OK == cs)
    {
        for(int m = 0; m < DIM; m++)
        {
            for(int n = 0; n < DIM; n++)
            {
                quad_[m][n] = cr->recv_double(src);
            }
        }
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to received MolecularQuadrupole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
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

CommunicationStatus MolecularPolarizability::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = cr->send_data(dest);
    }
    if (CS_OK == cs)
    {
        for(int m = 0; m < DIM; m++)
        {
            for(int n = 0; n < DIM; n++)
            {
                cr->send_double(dest, alpha_[m][n]);
            }
        }
    }
    if (CS_OK == cs)
    {
        cs = cr->send_data(dest);
    }
    if (CS_OK == cs)
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
    if (CS_OK == cs)
    {
        cs = cr->recv_data(src);
    }
    if (CS_OK == cs)
    {
        for(int m = 0; m < DIM; m++)
        {
            for(int n = 0; n < DIM; n++)
            {
                alpha_[m][n] = cr->recv_double(src);
            }
        }
    }
    if (CS_OK == cs)
    {
        cs = cr->recv_data(src);
    }
    if (CS_OK == cs)
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

CommunicationStatus MolecularEnergy::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = cr->recv_data(src);
    }
    if (CS_OK == cs)
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
    if (CS_OK == cs)
    {
        cs = cr->send_data(dest);
    }
    if (CS_OK == cs)
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
    CommunicationStatus cs;

    cs = cr->recv_data(src);
    if (CS_OK == cs)
    {
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
    CommunicationStatus cs;

    cs = cr->send_data(dest);
    if (CS_OK == cs)
    {
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

void ElectrostaticPotential::set(const std::string &xyz_unit,
                                 const std::string &V_unit,
                                 int                espid, 
                                 double             x,
                                 double             y,
                                 double             z,
                                 double             V)
{ 
    xyzUnit_ = xyz_unit;
    vUnit_ = V_unit;
    espID_ = espid;
    x_ = x;
    y_ = y;
    z_ = z;
    V_ = V; 
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
