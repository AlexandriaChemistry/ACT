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

#include "molprop.h"

#include <cmath>

#include <string>
#include <vector>

#include "gromacs/math/utilities.h"
#include "gromacs/utility/fatalerror.h"

#include "communication.h"
#include "composition.h"
#include "gmx_simple_comm.h"
#include "units.h"

const char *mpo_name[MPO_NR] =
{
    "potential", "dipole", "quadrupole", "polarizability", "energy", "entropy", "charge"
};

const char *mpo_unit[MPO_NR] =
{
    "e/nm", "D", "B", "\\AA$^3$", "kJ/mol", "J/mol K"
};

namespace alexandria
{

const char *dataSourceName(DataSource ds)
{
    switch (ds)
    {
        case dsExperiment:
            return "Experiment";
        case dsTheory:
            return "Theory";
    }
    return nullptr;
}

DataSource dataSourceFromName(const std::string &name)
{
    if (strcasecmp(dataSourceName(dsExperiment), name.c_str()) == 0)
    {
        return dsExperiment;
    }
    else if (strcasecmp(dataSourceName(dsTheory), name.c_str()) == 0)
    {
        return dsTheory;
    }
    gmx_fatal(FARGS, "No data source corresponding to %s", name.c_str());
}

void GenericProperty::SetType(const std::string &type)
{
    if ((type_.size() == 0) && (type.size() > 0))
    {
        type_ = type;
    }
    else
    {
        if (type_.size() == 0)
        {
            fprintf(stderr, "Replacing GenericProperty type '%s' by '%s'\n", type_.c_str(), type.c_str());
        }
    }
}

void GenericProperty::SetUnit(const std::string &unit)
{
    if ((unit_.size() == 0) && (unit.size() > 0))
    {
        unit_ = unit;
    }
    else
    {
        if (unit_.size() == 0)
        {
            fprintf(stderr, "Replacing GenericProperty unit '%s' by '%s'\n", unit_.c_str(), unit.c_str());
        }
    }
}

CommunicationStatus GenericProperty::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &type_);
        gmx_send_str(cr, dest, &unit_);
        gmx_send_double(cr, dest, T_);
        gmx_send_int(cr, dest, (int) eP_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send GenericProperty, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus GenericProperty::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &type_);
        gmx_recv_str(cr, src, &unit_);
        T_  = gmx_recv_double(cr, src);
        eP_ = (ePhase) gmx_recv_int(cr, src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive GenericProperty, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

void CalcAtom::SetUnit(const std::string &unit)
{
    if ((unit_.size() == 0) && (unit.size() > 0))
    {
        unit_ = unit;
    }
    else
    {
        if (unit_.size() == 0)
        {
            fprintf(stderr, "Replacing CalcAtom unit '%s' by '%s'\n", unit_.c_str(), unit.c_str());
        }
    }
}

void CalcAtom::AddCharge(AtomicCharge q)
{
    AtomicChargeIterator aci;

    auto acc = atomicChargeConst();
    for (aci = acc.begin(); aci < acc.end(); ++aci)
    {
        if ((aci->getType().compare(q.getType()) == 0) &&
            (aci->getUnit().compare(q.getUnit()) == 0) &&
            (aci->getTemperature() == q.getTemperature()) &&
            (aci->getQ() == q.getQ()))
        {
            break;
        }
    }
    if (aci == acc.end())
    {
        q_.push_back(q);
    }
}

bool CalcAtom::Equal(CalcAtom ca)
{
    return !((name_.compare(ca.getName()) != 0) ||
             (obType_.compare(ca.getObtype()) != 0) ||
             (x_ != ca.getX()) ||
             (y_ != ca.getY()) ||
             (z_ != ca.getZ()) ||
             (atomID_ != ca.getAtomid()));
}

CommunicationStatus CalcAtom::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    int                 Ncharge;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &name_);
        gmx_recv_str(cr, src, &obType_);
        gmx_recv_str(cr, src, &residueName_);
        residueNumber_ = gmx_recv_int(cr, src);
        atomID_ = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &unit_);
        x_      = gmx_recv_double(cr, src);
        y_      = gmx_recv_double(cr, src);
        z_      = gmx_recv_double(cr, src);
        Ncharge = gmx_recv_int(cr, src);

        for (int n = 0; (CS_OK == cs) && (n < Ncharge); n++)
        {
            AtomicCharge aq;
            cs = aq.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddCharge(aq);
            }
        }
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Received CalcAtom, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus CalcAtom::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus  cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &name_);
        gmx_send_str(cr, dest, &obType_);
        gmx_send_str(cr, dest, &residueName_);
        gmx_send_int(cr, dest, residueNumber_);
        gmx_send_int(cr, dest, atomID_);
        gmx_send_str(cr, dest, &unit_);
        gmx_send_double(cr, dest, x_);
        gmx_send_double(cr, dest, y_);
        gmx_send_double(cr, dest, z_);
        gmx_send_int(cr, dest, q_.size());

        for (auto qi : atomicChargeConst())
        {
            cs = qi.Send(cr, dest);
            if (CS_OK != cs)
            {
                break;
            }
        }
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Sent CalcAtom, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Bond::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, ai_);
        gmx_send_int(cr, dest, aj_);
        gmx_send_int(cr, dest, bondorder_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send Bond, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Bond::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        ai_        = gmx_recv_int(cr, src);
        aj_        = gmx_recv_int(cr, src);
        bondorder_ = gmx_recv_int(cr, src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive Bond, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus AtomNum::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &catom_);
        gmx_send_int(cr, dest, cnumber_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent AtomNum %s %d, status %s\n",
                    catom_.c_str(), cnumber_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus AtomNum::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &catom_);
        cnumber_ = gmx_recv_int(cr, src);
        if (nullptr != debug)
        {
            fprintf(debug, "Received AtomNum %s %d, status %s\n",
                    catom_.c_str(), cnumber_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

void MolecularComposition::AddAtom(AtomNum an)
{
    AtomNumIterator mci = searchAtom(an.getAtom());
    if (mci == atomnum_.end())
    {
        atomnum_.push_back(an);
    }
    else
    {
        mci->SetNumber(mci->getNumber()+an.getNumber());
    }
}

void MolecularComposition::DeleteAtom(const std::string &catom)
{
    AtomNumIterator ani;

    if ((ani = searchAtom(catom)) != atomnum_.end())
    {
        atomnum_.erase(ani);
    }
}

AtomNumConstIterator MolecularComposition::searchAtomConst(const std::string &an) const
{
    for (auto ani = atomnum_.begin(); ani < atomnum_.end(); ++ani)
    {
        if (an.compare(ani->getAtom()) == 0)
        {
            return ani;
        }
    }
    return atomnum_.end();
}

AtomNumIterator MolecularComposition::searchAtom(const std::string &an)
{
    for (auto ani = atomnum_.begin(); ani < atomnum_.end(); ++ani)
    {
        if (an.compare(ani->getAtom()) == 0)
        {
            return ani;
        }
    }
    return atomnum_.end();
}

void MolecularComposition::ReplaceAtom(const std::string &oldatom,
                                       const std::string &newatom)
{

    for (auto &i : atomnum_)
    {
        if (oldatom.compare(i.getAtom()) == 0)
        {
            i.SetAtom(newatom);
            break;
        }
    }
}

int MolecularComposition::CountAtoms(const std::string &atom) const
{
    for (auto &i : atomnum_)
    {
        if (atom.compare(i.getAtom()) == 0)
        {
            return i.getNumber();
        }
    }
    return 0;
}

int MolecularComposition::CountAtoms() const
{
    int             nat = 0;

    for (auto &i : atomnum_)
    {
        nat += i.getNumber();
    }
    return nat;
}

CommunicationStatus MolecularComposition::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, atomnum_.size());
        gmx_send_str(cr, dest, &compname_);
        for (auto &ani : atomnum_)
        {
            cs = ani.Send(cr, dest);
            if (CS_OK != cs)
            {
                break;
            }
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Sent MolecularComposition %s, status %s\n",
                    compname_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus MolecularComposition::Receive(t_commrec *cr, int src)
{
    int                 Natomnum;
    CommunicationStatus cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        Natomnum = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &compname_);
        CommunicationStatus cs2;
        for (int n = 0; n < Natomnum; n++)
        {
            AtomNum an;
            cs2 = an.Receive(cr, src);
            if (CS_OK == cs2)
            {
                AddAtom(an);
            }
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Received MolecularComposition %s, status %s\n",
                    compname_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus ElectrostaticPotential::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &xyzUnit_);
        gmx_recv_str(cr, src, &vUnit_);
        espID_ = gmx_recv_int(cr, src);
        x_     = gmx_recv_double(cr, src);
        y_     = gmx_recv_double(cr, src);
        z_     = gmx_recv_double(cr, src);
        V_     = gmx_recv_double(cr, src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive ElectrostaticPotential, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus ElectrostaticPotential::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &xyzUnit_);
        gmx_send_str(cr, dest, &vUnit_);
        gmx_send_int(cr, dest, espID_);
        gmx_send_double(cr, dest, x_);
        gmx_send_double(cr, dest, y_);
        gmx_send_double(cr, dest, z_);
        gmx_send_double(cr, dest, V_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send ElectrostaticPotential, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus AtomicCharge::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);

    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        q_ = gmx_recv_double(cr, src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive AtomicCharge, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus AtomicCharge::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);

    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, q_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send AtomicCharge, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularDipole::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, _x);
        gmx_send_double(cr, dest, _y);
        gmx_send_double(cr, dest, _z);
        gmx_send_double(cr, dest, _aver);
        gmx_send_double(cr, dest, error_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send MolecularDipole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularDipole::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        _x     = gmx_recv_double(cr, src);
        _y     = gmx_recv_double(cr, src);
        _z     = gmx_recv_double(cr, src);
        _aver  = gmx_recv_double(cr, src);
        error_ = gmx_recv_double(cr, src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive MolecularDipole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularQuadrupole::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, xx_);
        gmx_send_double(cr, dest, yy_);
        gmx_send_double(cr, dest, zz_);
        gmx_send_double(cr, dest, xy_);
        gmx_send_double(cr, dest, xz_);
        gmx_send_double(cr, dest, yz_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send MolecularQuadrupole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularQuadrupole::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        xx_    = gmx_recv_double(cr, src);
        yy_    = gmx_recv_double(cr, src);
        zz_    = gmx_recv_double(cr, src);
        xy_    = gmx_recv_double(cr, src);
        xz_    = gmx_recv_double(cr, src);
        yz_    = gmx_recv_double(cr, src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to received MolecularQuadrupole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

void MolecularPolarizability::Set(double xx, double yy, double zz,
                                  double xy, double xz, double yz,
                                  double average, double error)
{
    xx_      = xx; yy_ = yy; zz_ = zz;
    xy_      = xy; xz_ = xz; yz_ = yz;
    average_ = average; error_ = error;
    if (average_ == 0)
    {
        // Compute average as the 1/3 the trace of the diagonal
        average_ = (xx_ + yy_ + zz_)/3.0;
    }
    else if ((xx_ == 0) && (yy_ == 0) && (zz_ == 0))
    {
        // Estimate tensor as the 1/3 the trace of the diagonal
        xx_ = yy_ = zz_ = average_;
    }
}

CommunicationStatus MolecularPolarizability::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, xx_);
        gmx_send_double(cr, dest, yy_);
        gmx_send_double(cr, dest, zz_);
        gmx_send_double(cr, dest, xy_);
        gmx_send_double(cr, dest, xz_);
        gmx_send_double(cr, dest, yz_);
        gmx_send_double(cr, dest, average_);
        gmx_send_double(cr, dest, error_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send MolecularQuadrupole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularPolarizability::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        xx_      = gmx_recv_double(cr, src);
        yy_      = gmx_recv_double(cr, src);
        zz_      = gmx_recv_double(cr, src);
        xy_      = gmx_recv_double(cr, src);
        xz_      = gmx_recv_double(cr, src);
        yz_      = gmx_recv_double(cr, src);
        average_ = gmx_recv_double(cr, src);
        error_   = gmx_recv_double(cr, src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to received MolecularQuadrupole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularEnergy::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        _value = gmx_recv_double(cr, src);
        error_ = gmx_recv_double(cr, src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive MolecularEnergy, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularEnergy::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, _value);
        gmx_send_double(cr, dest, error_);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send MolecularEnergy, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

Experiment::Experiment(const std::string &program,
                       const std::string &method,
                       const std::string &basisset,
                       const std::string &reference,
                       const std::string &conformation,
                       const std::string &datafile,
                       jobType            jtype)
    :
      dataSource_(dsTheory),
      reference_(reference),
      conformation_(conformation),
      program_(program),
      method_(method),
      basisset_(basisset),
      datafile_(datafile),
      jobtype_(jtype)

{}

void Experiment::Dump(FILE *fp) const
{
    if (nullptr != fp)
    {
        fprintf(fp, "Experiment %s %d polar %d dipole\n",
                dataSourceName(dataSource()), NPolar(), NDipole());
        if (dsExperiment == dataSource())
        {
            fprintf(fp, "reference    = %s\n", reference_.c_str());
            fprintf(fp, "conformation = %s\n", conformation_.c_str());
        }
        else
        {
            fprintf(fp, "program    = %s\n", program_.c_str());
            fprintf(fp, "method     = %s\n", method_.c_str());
            fprintf(fp, "basisset   = %s\n", basisset_.c_str());
            fprintf(fp, "datafile   = %s\n", datafile_.c_str());
            for (auto &cai : calcAtomConst())
            {
                double   x, y, z;
                cai.getCoords(&x, &y, &z);
                fprintf(fp, "%-3s  %-3s  %3d  %10.3f  %10.3f  %10.3f\n",
                        cai.getName().c_str(), cai.getObtype().c_str(),
                        cai.getAtomid(), x, y, z);
            }
        }
    }
}

bool Experiment::getHF(double *value) const
{
    double      T    = -1;
    double      error;
    rvec        vec;
    bool        done = false;
    tensor      quad;

    if (getVal((char *)"HF", MPO_ENERGY, value, &error, &T, vec, quad))
    {
        done = true;
    }
    return done;
}

int Experiment::Merge(const Experiment *src)
{
    int nwarn = 0;

    for (auto &mei : src->molecularEnergyConst())
    {
        alexandria::MolecularEnergy me(mei.getType(),
                                       mei.getUnit(),
                                       mei.getTemperature(),
                                       mei.getPhase(),
                                       mei.getValue(),
                                       mei.getError());
        AddEnergy(me);
    }

    for (auto &dpi : src->dipoleConst())
    {
        alexandria::MolecularDipole dp(dpi.getType(),
                                       dpi.getUnit(),
                                       dpi.getTemperature(),
                                       dpi.getX(), dpi.getY(), dpi.getZ(),
                                       dpi.getAver(), dpi.getError());
        AddDipole(dp);
    }

    for (auto &mpi : src->polarizabilityConst())
    {
        alexandria::MolecularPolarizability mp(mpi.getType(),
                                               mpi.getUnit(),
                                               mpi.getTemperature(),
                                               mpi.getXX(), mpi.getYY(), mpi.getZZ(),
                                               mpi.getXY(), mpi.getXZ(), mpi.getYZ(),
                                               mpi.getAverage(), mpi.getError());
        AddPolar(mp);
    }

    for (auto &mqi : src->quadrupoleConst())
    {
        alexandria::MolecularQuadrupole mq(mqi.getType(), mqi.getUnit(),
                                           mqi.getTemperature(),
                                           mqi.getXX(), mqi.getYY(), mqi.getZZ(),
                                           mqi.getXY(), mqi.getXZ(), mqi.getYZ());
        AddQuadrupole(mq);
    }

    for (auto &cai : src->calcAtomConst())
    {
        double   x, y, z;
        CalcAtom caa(cai.getName(), cai.getObtype(), cai.getAtomid());

        cai.getCoords(&x, &y, &z);
        caa.SetCoords(x, y, z);
        caa.SetUnit(cai.getUnit());
        caa.SetResidue(cai.ResidueName(), cai.ResidueNumber());
        for (auto &aci : cai.atomicChargeConst())
        {
            AtomicCharge aq(aci.getType(), aci.getUnit(),
                            aci.getTemperature(), aci.getQ());
            caa.AddCharge(aq);
        }
        AddAtom(caa);
    }

    for (auto &mep : src->electrostaticPotentialConst())
    {
        alexandria::ElectrostaticPotential ep(mep.getXYZunit(), mep.getVunit(),
                                              mep.getEspid(),
                                              mep.getX(), mep.getY(),
                                              mep.getZ(), mep.getV());
        AddPotential(ep);
    }

    return nwarn;
}

CalcAtomIterator Experiment::searchAtom(CalcAtom ca)
{
    CalcAtomIterator cai;
    for (auto cai = catom_.begin(); (cai < catom_.end()); ++cai)
    {
        if (cai->Equal(ca))
        {
            return cai;
        }
    }
    return catom_.end();
}

void Experiment::AddAtom(CalcAtom ca)
{
    CalcAtomIterator cai = searchAtom(ca);

    if (cai == catom_.end())
    {
        gmx::RVec x;
        auto      unit = ca.getUnit();
        x[XX]          = convertToGromacs(ca.getX(), unit);
        x[YY]          = convertToGromacs(ca.getY(), unit);
        x[ZZ]          = convertToGromacs(ca.getZ(), unit);
        coordinates_.push_back(x);
        catom_.push_back(ca);
    }
    else
    {
        printf("Trying to add identical atom %s (%s) twice. N = %d\n",
               ca.getName().c_str(), ca.getObtype().c_str(),
               (int)catom_.size());
    }
}

bool Experiment::getVal(const std::string &type,
                        MolPropObservable  mpo,
                        double            *value,
                        double            *error,
                        double            *T,
                        rvec               vec,
                        tensor             quad_polar) const
{
    bool   done = false;
    double x, y, z;
    double Told = *T;

    switch (mpo)
    {
        case MPO_ENERGY:
        case MPO_ENTROPY:
            for (auto &mei : molecularEnergyConst())
            {
                if (((type.size() == 0) || (type.compare(mei.getType()) == 0)) &&
                    bCheckTemperature(Told, mei.getTemperature()))
                {
                    mei.get(value, error);
                    *T   = mei.getTemperature();
                    done = true;
                    break;
                }
            }
            break;
        case MPO_DIPOLE:
            for (auto &mdp : dipoleConst())
            {
                if (((type.size() == 0) || (type.compare(mdp.getType()) == 0))  &&
                    bCheckTemperature(Told, mdp.getTemperature()))
                {
                    mdp.get(&x, &y, &z, value, error);
                    vec[XX] = x;
                    vec[YY] = y;
                    vec[ZZ] = z;
                    *T      = mdp.getTemperature();
                    done    = true;
                    break;
                }
            }
            break;
        case MPO_POLARIZABILITY:
        {
            for (auto &mdp : polarizabilityConst())
            {
                if (((type.size() == 0) || (type.compare(mdp.getType()) == 0)) &&
                    bCheckTemperature(Told, mdp.getTemperature()))
                {
                    double xx, yy, zz, xy, xz, yz;
                    mdp.get(&xx, &yy, &zz, &xy, &xz, &yz, value, error);
                    quad_polar[XX][XX] = xx;
                    quad_polar[XX][YY] = xy;
                    quad_polar[XX][ZZ] = xz;
                    quad_polar[YY][XX] = 0;
                    quad_polar[YY][YY] = yy;
                    quad_polar[YY][ZZ] = yz;
                    quad_polar[ZZ][XX] = 0;
                    quad_polar[ZZ][YY] = 0;
                    quad_polar[ZZ][ZZ] = zz;
                    *T                 = mdp.getTemperature();
                    done               = true;
                    break;
                }
            }
        }
        break;
        case MPO_QUADRUPOLE:
            for (auto &mqi : quadrupoleConst())
            {
                if (((type.size() == 0) || (type.compare(mqi.getType()) == 0)) &&
                    bCheckTemperature(Told, mqi.getTemperature()))
                {
                    double xx, yy, zz, xy, xz, yz;
                    mqi.get(&xx, &yy, &zz, &xy, &xz, &yz);
                    quad_polar[XX][XX] = xx;
                    quad_polar[XX][YY] = xy;
                    quad_polar[XX][ZZ] = xz;
                    quad_polar[YY][XX] = 0;
                    quad_polar[YY][YY] = yy;
                    quad_polar[YY][ZZ] = yz;
                    quad_polar[ZZ][XX] = 0;
                    quad_polar[ZZ][YY] = 0;
                    quad_polar[ZZ][ZZ] = zz;
                    *T                 = mqi.getTemperature();
                    done               = true;
                    break;
                }
            }
            break;
        case MPO_CHARGE:
        {
            int i = 0;
            for (auto &mai : calcAtomConst())
            {
                for (auto &q : mai.atomicChargeConst())
                {
                    if (((type.size() == 0) || (type.compare(q.getType()) == 0)) &&
                        bCheckTemperature(Told, q.getTemperature()))
                    {
                        vec[i] = q.getQ();
                        i++;
                    }
                }
            }
            if (i == NAtom())
            {
                done = true;
            }
        }
        break;
    default:
        break;
    }
    return done;
}

CommunicationStatus Experiment::Receive(t_commrec *cr, int src)
{
    CalcAtomIterator               cai;
    CommunicationStatus            cs;
    ElectrostaticPotentialIterator epi;
    std::string                    jobtype;
    int                            Npolar, Ndipole, Nenergy, Npotential, Natom, Nquadrupole;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        dataSource_ = static_cast<DataSource>(gmx_recv_int(cr, src));
        gmx_recv_str(cr, src, &reference_);
        gmx_recv_str(cr, src, &conformation_);
        gmx_recv_str(cr, src, &program_);
        gmx_recv_str(cr, src, &method_);
        gmx_recv_str(cr, src, &basisset_);
        gmx_recv_str(cr, src, &datafile_);
        gmx_recv_str(cr, src, &jobtype);
        jobtype_    = string2jobType(jobtype);
        Npolar      = gmx_recv_int(cr, src);
        Ndipole     = gmx_recv_int(cr, src);
        Nquadrupole = gmx_recv_int(cr, src);
        Nenergy     = gmx_recv_int(cr, src);
        Npotential  = gmx_recv_int(cr, src);
        Natom       = gmx_recv_int(cr, src);

        //! Receive Polarizabilities
        for (int n = 0; (CS_OK == cs) && (n < Npolar); n++)
        {
            MolecularPolarizability mp;
            cs = mp.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddPolar(mp);
            }
        }

        //! Receive Dipoles
        for (int n = 0; (CS_OK == cs) && (n < Ndipole); n++)
        {
            MolecularDipole md;
            cs = md.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddDipole(md);
            }
        }

        //! Receive Quadrupoles
        for (int n = 0; (CS_OK == cs) && (n < Nquadrupole); n++)
        {
            MolecularQuadrupole mq;
            cs = mq.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddQuadrupole(mq);
            }
        }

        //! Receive Energies
        for (int n = 0; (CS_OK == cs) && (n < Nenergy); n++)
        {
            MolecularEnergy me;
            cs = me.Receive(cr, src);
            if  (CS_OK == cs)
            {
                AddEnergy(me);
            }
        }

        //! Receive Potentials
        for (int n = 0; (CS_OK == cs) && (n < Npotential); n++)
        {
            ElectrostaticPotential ep;
            cs = ep.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddPotential(ep);
            }
        }

        //! Receive Atoms
        for (int n = 0; (CS_OK == cs) && (n < Natom); n++)
        {
            CalcAtom ca;
            cs = ca.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddAtom(ca);
            }
        }
    }

    if ((CS_OK != cs) && (nullptr != debug))
    {
        fprintf(debug, "Trying to receive Experiment, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Experiment::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;
    std::string         jobtype;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, static_cast<int>(dataSource_));
        gmx_send_str(cr, dest, &reference_);
        gmx_send_str(cr, dest, &conformation_);
        gmx_send_str(cr, dest, &program_);
        gmx_send_str(cr, dest, &method_);
        gmx_send_str(cr, dest, &basisset_);
        gmx_send_str(cr, dest, &datafile_);
        jobtype.assign(jobType2string(jobtype_));
        gmx_send_str(cr, dest, &jobtype);
        gmx_send_int(cr, dest, polar_.size());
        gmx_send_int(cr, dest, dipole_.size());
        gmx_send_int(cr, dest, quadrupole_.size());
        gmx_send_int(cr, dest, energy_.size());
        gmx_send_int(cr, dest, potential_.size());
        gmx_send_int(cr, dest, catom_.size());

        //! Send Polarizabilities
        for (auto &dpi : polarizabilityConst())
        {
            cs = dpi.Send(cr, dest);
            if (CS_OK != cs)
            {
                break;
            }
        }

        //! Send Dipoles
        if (CS_OK == cs)
        {
            for (auto &dpi : dipoleConst())
            {
                cs = dpi.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }

        //! Send Quadrupoles
        if (CS_OK == cs)
        {
            for (auto &dqi : quadrupoleConst())
            {
                cs = dqi.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }
        
        //! Send Energies
        if (CS_OK == cs)
        {
            for (auto &mei : molecularEnergyConst())
            {
                cs = mei.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }

        //! Send Potentials
        if (CS_OK == cs)
        {
            for (auto &epi : electrostaticPotentialConst())
            {
                cs = epi.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }

        //! Send Atoms
        if (CS_OK == cs)
        {
            for (auto &cai : calcAtomConst())
            {
                cs = cai.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }
    }

    if ((CS_OK != cs) && (nullptr != debug))
    {
        fprintf(debug, "Trying to send Experiment, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

int MolProp::NAtom()
{
    if (mol_comp_.size() > 0)
    {
        int nat = mol_comp_[0].CountAtoms();
        return nat;
    }
    return 0;
}

void MolProp::AddBond(Bond b)
{
    BondConstIterator bi;
    bool              bFound = false;

    for (bi = bondConst().begin(); bi < bondConst().end(); ++bi)
    {
        bFound = (((bi->getAi() == b.getAi()) && (bi->getAj() == b.getAj())) ||
                  ((bi->getAi() == b.getAj()) && (bi->getAj() == b.getAi())));
        if (bFound)
        {
            break;
        }
    }
    if (!bFound)
    {
        bond_.push_back(b);
    }
    else if ((nullptr != debug) && (bi->getBondOrder() != b.getBondOrder()))
    {
        fprintf(debug, "Different bond orders in molecule %s\n", getMolname().c_str());
        fflush(debug);
    }
}

void MolProp::CheckConsistency()
{
}

bool MolProp::SearchCategory(const std::string &catname) const
{
    for (auto &i : category_)
    {
        if (strcasecmp(i.c_str(), catname.c_str()) == 0)
        {
            return true;
        }
    }
    return false;
}

void MolProp::DeleteComposition(const std::string &compname)
{
    std::remove_if(mol_comp_.begin(), mol_comp_.end(),
                   [compname](MolecularComposition mc)
                   {
                       return (compname.compare(mc.getCompName()) == 0);
                   });
}

void MolProp::AddComposition(MolecularComposition mc)
{
    MolecularCompositionIterator mci = SearchMolecularComposition(mc.getCompName());
    if (mci == mol_comp_.end())
    {
        mol_comp_.push_back(mc);
    }
}

bool MolProp::BondExists(Bond b)
{
    for (auto &bi : bondConst())
    {
        if (((bi.getAi() == b.getAi()) && (bi.getAj() == b.getAj())) ||
            ((bi.getAi() == b.getAj()) && (bi.getAj() == b.getAi())))
        {
            return true;
        }
    }
    return false;
}

int MolProp::Merge(const MolProp *src)
{
    double      q = 0, sq = 0;
    std::string stmp;
    int         nwarn = 0;

    for (auto &si : src->categoryConst())
    {
        AddCategory(si);
    }
    SetFormula(src->formula());
    SetMass(src->getMass());
    SetIndex(src->getIndex());
    if (getMultiplicity() <= 1)
    {
        SetMultiplicity(src->getMultiplicity());
    }
    else
    {
        int smult = src->getMultiplicity();
        if ((nullptr != debug) && (smult != getMultiplicity()))
        {
            fprintf(debug, "Not overriding multiplicity to %d when merging since it is %d (%s)\n",
                    smult, getMultiplicity(), src->getMolname().c_str());
            fflush(debug);
        }
    }
    q = getCharge();
    if (q == 0)
    {
        SetCharge(src->getCharge());
    }
    else
    {
        sq = src->getCharge();
        if ((nullptr != debug) && (sq != q))
        {
            fprintf(debug, "Not overriding charge to %g when merging since it is %g (%s)\n",
                    sq, q, getMolname().c_str());
            fflush(debug);
        }
    }

    stmp = src->getMolname();
    if ((getMolname().size() == 0) && (stmp.size() != 0))
    {
        SetMolname(stmp);
    }
    stmp = src->getIupac();
    if ((getIupac().size() == 0) && (stmp.size() != 0))
    {
        SetIupac(stmp);
    }
    stmp = src->getCas();
    if ((getCas().size() == 0) && (stmp.size() != 0))
    {
        SetCas(stmp);
    }
    stmp = src->getCid();
    if ((getCid().size() == 0) && (stmp.size() != 0))
    {
        SetCid(stmp);
    }
    stmp = src->getInchi();
    if ((getInchi().size() == 0) && (stmp.size() != 0))
    {
        SetInchi(stmp);
    }
    if (NBond() == 0)
    {
        for (auto &bi : src->bondConst())
        {
            alexandria::Bond bb(bi.getAi(), bi.getAj(), bi.getBondOrder());
            AddBond(bb);
        }
    }
    else
    {
        for (auto &bi : src->bondConst())
        {
            alexandria::Bond bb(bi.getAi(), bi.getAj(), bi.getBondOrder());
            if (!BondExists(bb))
            {
                fprintf(stderr, "WARNING bond %d-%d not present in %s\n",
                        bi.getAi(), bi.getAj(), getMolname().c_str());
                nwarn++;
            }
        }
    }

    for (auto &ei : src->experimentConst())
    {
        if (dsExperiment == ei.dataSource())
        {
            Experiment ex(ei.getReference(), ei.getConformation());
            nwarn += ex.Merge(&ei);
            AddExperiment(ex);
        }
        else
        {
            auto jtype = ei.getJobtype();
            if (jtype != JOB_NR)
            {
                Experiment ca(ei.getProgram(), ei.getMethod(),
                              ei.getBasisset(), ei.getReference(),
                              ei.getConformation(), ei.getDatafile(),
                              jtype);
                nwarn += ca.Merge(&ei);
                AddExperiment(ca);
            }
        }
    }
    for (auto &mci : src->molecularCompositionConst())
    {
        alexandria::MolecularComposition mc(mci.getCompName());

        for (auto &ani : mci.atomNumConst())
        {
            AtomNum an(ani.getAtom(), ani.getNumber());
            mc.AddAtom(an);
        }
        AddComposition(mc);
    }
    return nwarn;
}

MolecularCompositionIterator MolProp::SearchMolecularComposition(const std::string &str)
{
    return std::find_if(mol_comp_.begin(), mol_comp_.end(),
                        [str](MolecularComposition const &mc)
                        {
                            return (str.compare(mc.getCompName()) == 0);
                        });
}

MolecularCompositionConstIterator MolProp::SearchMolecularComposition(const std::string &str) const
{
    return std::find_if(mol_comp_.begin(), mol_comp_.end(),
                        [str](MolecularComposition const &mc)
                        {
                            return (str.compare(mc.getCompName()) == 0);
                        });
}

void MolProp::Dump(FILE *fp) const
{
    if (fp)
    {
        fprintf(fp, "formula:      %s\n", formula().c_str());
        fprintf(fp, "molname:      %s\n", getMolname().c_str());
        fprintf(fp, "iupac:        %s\n", getIupac().c_str());
        fprintf(fp, "CAS:          %s\n", getCas().c_str());
        fprintf(fp, "cis:          %s\n", getCid().c_str());
        fprintf(fp, "InChi:        %s\n", getInchi().c_str());
        fprintf(fp, "mass:         %g\n", getMass());
        fprintf(fp, "charge:       %d\n", getCharge());
        fprintf(fp, "multiplicity: %d\n", getMultiplicity());
        fprintf(fp, "category:    ");
        for (auto &si : categoryConst())
        {
            fprintf(fp, " '%s'", si.c_str());
        }
        fprintf(fp, "\n");
        for (auto &ei : experimentConst())
        {
            ei.Dump(fp);
        }
    }
}

bool MolProp::GenerateComposition()
{
    ExperimentIterator   ci;
    CalcAtomIterator     cai;
    CompositionSpecs     cs;
    MolecularComposition mci_alexandria(cs.searchCS(iCalexandria)->name());

    // Why was this again?
    mol_comp_.clear();

    int natoms = 0;
    for (auto &ci : experimentConst())
    {
        /* This assumes we have either all atoms or none.
         * A consistency check could be
         * to compare the number of atoms to the formula */
        int nat = 0;
        for (auto &cai : ci.calcAtomConst())
        {
            nat++;
            AtomNum ans(cai.getObtype(), 1);
            mci_alexandria.AddAtom(ans);
        }
        natoms = std::max(natoms, nat);
        if (mci_alexandria.CountAtoms() > 0)
        {
            break;
        }
    }

    if (natoms == mci_alexandria.CountAtoms())
    {
        AddComposition(mci_alexandria);
        if (nullptr != debug)
        {
            fprintf(debug, "LO_COMP: ");
            for (auto &ani : mci_alexandria.atomNumConst())
            {
                fprintf(debug, " %s:%d",
                        ani.getAtom().c_str(), ani.getNumber());
            }
            fprintf(debug, "\n");
            fflush(debug);
        }
        return true;
    }
    return false;
}

static void add_element_toformula_(const char *elem, int number, char *formula, char *texform)
{
    if (number > 0)
    {
        strcat(formula, elem);
        strcat(texform, elem);
        if (number > 1)
        {
            char cnumber[32];

            sprintf(cnumber, "%d", number);
            strcat(formula, cnumber);
            sprintf(cnumber, "$_{%d}$", number);
            strcat(texform, cnumber);
        }
    }
}

bool MolProp::GenerateFormula(gmx_atomprop_t ap)
{
    char             myform[1280], texform[2560];
    std::vector<int> ncomp;
    alexandria::MolecularCompositionIterator mci;

    ncomp.resize(110, 0);
    myform[0]  = '\0';
    texform[0] = '\0';
    mci        = SearchMolecularComposition("bosque");
    if (mci != mol_comp_.end())
    {
        for (auto &ani : mci->atomNumConst())
        {
            int         cnumber, an;
            real        value;
            std::string catom = ani.getAtom();
            cnumber = ani.getNumber();
            if (gmx_atomprop_query(ap, epropElement, "???", catom.c_str(), &value))
            {
                an = std::lround(value);
                range_check(an, 0, 110);
                if (an > 0)
                {
                    ncomp[an] += cnumber;
                }
            }
        }
    }
    add_element_toformula_("C", ncomp[6], myform, texform);
    add_element_toformula_("H", ncomp[1], myform, texform);
    ncomp[6] = ncomp[1] = 0;

    for (int j = 109; (j >= 1); j--)
    {
        add_element_toformula_(gmx_atomprop_element(ap, j), ncomp[j], myform, texform);
    }
    std::string mform = formula();
    if (strlen(myform) > 0)
    {
        if (debug)
        {
            if ((mform.size() > 0) && (strcasecmp(myform, mform.c_str()) != 0))
            {
                fprintf(debug, "Formula '%s' does match '%s' based on composition for %s.\n",
                        mform.c_str(), myform, getMolname().c_str());
                fflush(debug);
            }
        }
        SetFormula(myform);
        SetTexFormula(texform);
    }
    else if ((mform.size() == 0) && debug)
    {
        fprintf(debug, "Empty composition and formula for %s\n",
                getMolname().c_str());
        fflush(debug);
    }

    return (strlen(myform) > 0);
}

bool MolProp::HasComposition(const std::string &composition) const
{
    return std::find_if(mol_comp_.begin(), mol_comp_.end(),
                        [composition](const MolecularComposition &mi)
                        {
                            return mi.getCompName().compare(composition) == 0;
                        }) != mol_comp_.end();
}


bool bCheckTemperature(double Tref, double T)
{
    return (Tref < 0) || (fabs(T - Tref) < 0.05);
}

static bool stringEqual(const std::string &a, const std::string &b)
{
    unsigned int sz = a.size();
    if (b.size() != sz)
    {
        return false;
    }
    for (unsigned int i = 0; i < sz; ++i)
    {
        if (tolower(a[i]) != tolower(b[i]))
        {
            return false;
        }
    }
    return true;
}

bool MolProp::getPropRef(MolPropObservable mpo, iqmType iQM,
                         const std::string &method,
                         const std::string &basis,
                         const std::string &conf,
                         const std::string &type,
                         double *value, double *error, double *T,
                         std::string *ref, std::string *mylot,
                         rvec vec, tensor quad_polar)
{
    bool   done = false;
    double Told = *T;

    if (iQM == iqmBoth)
    {
        for (auto &ei : experimentConst())
        {
            if ((conf.size() == 0) ||
                stringEqual(ei.getConformation(), conf))
            {
                if (ei.getVal(type, mpo, value, error, T, vec, quad_polar) &&
                    bCheckTemperature(Told, *T))
                {
                    ref->assign(ei.getReference());
                    mylot->assign("Experiment");
                    done = true;
                    break;
                }
            }
        }
    }
    if (iQM == iqmExp)
    {
        for (auto &ei : experimentConst())
        {
            if (dsExperiment != ei.dataSource())
            {
                continue;
            }
            if ((conf.size() == 0) ||
                stringEqual(ei.getConformation(), conf))
            {
                if (ei.getVal(type, mpo, value, error, T, vec, quad_polar) &&
                    bCheckTemperature(Told, *T))
                {
                    ref->assign(ei.getReference());
                    mylot->assign("Experiment");
                    done = true;
                    break;
                }
            }
        }
    }
    else if (iQM == iqmQM)
    {
        for (auto &ci : experimentConst())
        {
            if (dsExperiment == ci.dataSource())
            {
                continue;
            }
            if (((method.size() == 0 || method == ci.getMethod()) &&
                 (basis.size() == 0  || basis == ci.getBasisset())) &&
                (conf.size() == 0    || conf == ci.getConformation()))
            {
                if  (ci.getVal(type.c_str(), mpo, value, error, T, vec, quad_polar) &&
                     bCheckTemperature(Told, *T))
                {
                    ref->assign(ci.getReference());
                    mylot->assign(method + "/" + basis);
                    done = true;
                    break;
                }
            }
        }
    }
    return done;
}

bool MolProp::getOptHF(double *value)
{
    bool done = false;

    for (auto &ei : experimentConst())
    {
        if (ei.getJobtype() == JOB_OPT)
        {
            if (ei.getHF(value))
            {
                done = true;
            }
        }
    }
    return done;
}

int MolProp::NOptSP()
{
    int n = 0;

    for (auto &ei : experimentConst())
    {
        if (ei.getJobtype() == JOB_OPT ||
            ei.getJobtype() == JOB_SP)
        {
            n++;
        }
    }
    return n;
}

bool MolProp::getProp(MolPropObservable mpo, iqmType iQM,
                      const std::string &method,
                      const std::string &basis,
                      const std::string &conf,
                      const std::string &type,
                      double *value, double *error, double *T)
{
    double      myerror;
    rvec        vec;
    tensor      quad;
    bool        bReturn;
    std::string myref, mylot;

    bReturn = getPropRef(mpo, iQM, method, basis,
                         conf, type, value, &myerror, T,
                         &myref, &mylot, vec, quad);
    if (nullptr != error)
    {
        *error = myerror;
    }
    return bReturn;
}


ExperimentIterator MolProp::getCalcPropType(const std::string &method,
                                            const std::string &basis,
                                            std::string       *mylot,
                                            MolPropObservable  mpo,
                                            const char        *type)
{
    ExperimentIterator ci;

    for (ci = exper_.begin(); (ci < exper_.end()); ci++)
    {
        if ((method.size() == 0 || strcasecmp(method.c_str(), ci->getMethod().c_str()) == 0) &&
            (basis.size() == 0  || strcasecmp(basis.c_str(), ci->getBasisset().c_str()) == 0))
        {
            bool done = false;
            switch (mpo)
            {
                case MPO_POTENTIAL:
                    done = ci->NPotential() > 0;
                    break;
                case MPO_DIPOLE:
                    for (auto &mdp : ci->dipoleConst())
                    {
                        done = ((nullptr == type) ||
                                (strcasecmp(type, mdp.getType().c_str()) == 0));
                        if (done)
                        {
                            break;
                        }
                    }
                    break;
                case MPO_QUADRUPOLE:
                    for (auto &mdp : ci->quadrupoleConst())
                    {
                        done = ((nullptr == type) ||
                                (strcasecmp(type, mdp.getType().c_str()) == 0));
                        if (done)
                        {
                            break;
                        }
                    }
                    break;
                case MPO_POLARIZABILITY:
                    for (auto &mdp : ci->polarizabilityConst())
                    {
                        done = ((nullptr == type) ||
                                (strcasecmp(type, mdp.getType().c_str()) == 0));
                        if (done)
                        {
                            break;
                        }
                    }
                    break;
                case MPO_ENERGY:
                case MPO_ENTROPY:
                    for (auto &mdp : ci->molecularEnergyConst())
                    {
                        done = ((nullptr == type) ||
                                (strcasecmp(type, mdp.getType().c_str()) == 0));
                        if (done)
                        {
                            break;
                        }
                    }
                    break;
                default:
                    break;
            }
            if (done)
            {
                if (nullptr != mylot)
                {
                    mylot->assign(ci->getMethod() + "/" + ci->getBasisset());
                }
                break;
            }
        }
    }
    return ci;
}

ExperimentIterator MolProp::getCalc(const std::string &method,
                                    const std::string &basis,
                                    std::string       *mylot)
{
    ExperimentIterator ci;
    for (ci = exper_.begin(); (ci < exper_.end()); ++ci)
    {
        if ((method.size() == 0 || strcasecmp(method.c_str(), ci->getMethod().c_str()) == 0) &&
            (basis.size() == 0  || strcasecmp(basis.c_str(), ci->getBasisset().c_str()) == 0))
        {
            if (nullptr != mylot)
            {
                mylot->assign(ci->getMethod() + "/" + ci->getBasisset());
            }
            break;
        }
    }
    return ci;
}

CommunicationStatus MolProp::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus                cs;
    BondIterator                       bi;
    MolecularCompositionIterator       mci;
    std::vector<std::string>::iterator si;
    ExperimentIterator                 ei;

    /* Generic stuff */
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, mass_);
        gmx_send_int(cr, dest, index_);
        gmx_send_int(cr, dest, charge_);
        gmx_send_int(cr, dest, multiplicity_);
        gmx_send_str(cr, dest, &formula_);
        gmx_send_str(cr, dest, &molname_);
        gmx_send_str(cr, dest, &iupac_);
        gmx_send_str(cr, dest, &cas_);
        gmx_send_str(cr, dest, &cid_);
        gmx_send_str(cr, dest, &inchi_);
        gmx_send_int(cr, dest, bond_.size());
        gmx_send_int(cr, dest, mol_comp_.size());
        gmx_send_int(cr, dest, category_.size());
        gmx_send_int(cr, dest, exper_.size());

        /* Send Bonds */
        for (auto &bi : bondConst())
        {
            cs = bi.Send(cr, dest);
            if (CS_OK != cs)
            {
                break;
            }
        }

        /* Send Composition */
        if (CS_OK == cs)
        {
            for (auto &mci : molecularCompositionConst())
            {
                cs = mci.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }

        /* send Categories */
        if (CS_OK == cs)
        {
            for (auto &si : categoryConst())
            {
                cs = gmx_send_data(cr, dest);
                if (CS_OK == cs)
                {
                    std::string sii = si.c_str();
                    gmx_send_str(cr, dest, &sii);
                    if (nullptr != debug)
                    {
                        fprintf(debug, "Sent category %s\n", si.c_str());
                        fflush(debug);
                    }
                }
            }
        }

        /* Send Experiments */
        if (CS_OK == cs)
        {
            for (auto &ei : experimentConst())
            {
                cs = ei.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Sent MolProp %s, status %s\n",
                    getMolname().c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus MolProp::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    int                 Nbond, Nmol_comp, Ncategory, Nexper;

    /* Generic stuff */
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        //! Receive mass and more
        mass_         = gmx_recv_double(cr, src);
        index_        = gmx_recv_int(cr, src);
        charge_       = gmx_recv_int(cr, src);
        multiplicity_ = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &formula_);
        gmx_recv_str(cr, src, &molname_);
        gmx_recv_str(cr, src, &iupac_);
        gmx_recv_str(cr, src, &cas_);
        gmx_recv_str(cr, src, &cid_);
        gmx_recv_str(cr, src, &inchi_);
        Nbond     = gmx_recv_int(cr, src);
        Nmol_comp = gmx_recv_int(cr, src);
        Ncategory = gmx_recv_int(cr, src);
        Nexper    = gmx_recv_int(cr, src);

        if (nullptr != debug)
        {
            fprintf(debug, "Got molname %s\n", getMolname().c_str());
        }
        //! Receive Bonds
        for (int n = 0; (CS_OK == cs) && (n < Nbond); n++)
        {
            Bond b;
            cs = b.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddBond(b);
            }
        }

        //! Receive Compositions
        for (int n = 0; (CS_OK == cs) && (n < Nmol_comp); n++)
        {
            MolecularComposition mc;
            cs = mc.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddComposition(mc);
            }
        }

        //! Receive Categories
        for (int n = 0; (CS_OK == cs) && (n < Ncategory); n++)
        {
            cs = gmx_recv_data(cr, src);
            if (CS_OK == cs)
            {
                std::string str;
                gmx_recv_str(cr, src, &str);
                if (!str.empty())
                {
                    AddCategory(str);
                    if (nullptr != debug)
                    {
                        fprintf(debug, "Received a category %s\n", str.c_str());
                        fflush(debug);
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "A category was promised but I got a nullptr pointer");
                }
            }
        }

        //! Receive Experiments
        for (int n = 0; (CS_OK == cs) && (n < Nexper); n++)
        {
            Experiment ex;
            cs = ex.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddExperiment(ex);
            }
        }

        if (nullptr != debug)
        {
            fprintf(debug, "Received %d experiments from %d for mol %s\n",
                    NExperiment(), src, getMolname().c_str());
            fprintf(debug, "Received MolProp %s, status %s\n",
                    getMolname().c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

const std::string &MolProp::getTexFormula() const
{
    if (texform_.size() > 0)
    {
        return texform_;
    }
    else
    {
        return formula_;
    }
}

}
