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

#include "experiment.h"

#include "gromacs/utility/fatalerror.h"

#include "units.h"

namespace alexandria
{

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

Experiment::Experiment(const std::string &program,
                       const std::string &method,
                       const std::string &basisset,
                       const std::string &reference,
                       const std::string &conformation,
                       const std::string &datafile,
                       JobType            jtype)
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
        fprintf(fp, "Experiment %s\n", dataSourceName(dataSource()));
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
        for(const auto &p : property_)
        {
            fprintf(fp, "Property %s\n", mpo_name(p.first));
            for (const auto &gp : p.second)
            {
                fprintf(fp, "Type: %s Unit: %s T: %g K Phase: %s\n",
                        gp->getType(), gp->getUnit(),
                        gp->getTemperature(),
                        phase2string(gp->getPhase()).c_str());
            }
        }
    }
}

void Experiment::addProperty(MolPropObservable mpo, GenericProperty *gp)
{
    if (property_.find(mpo) == property_.end())
    {
        std::vector<GenericProperty *> gpnew;
        property_.insert(std::pair<MolPropObservable, std::vector<GenericProperty *>>(mpo, gpnew));
    }
    property_.find(mpo)->second.push_back(std::move(gp));
}
        
int Experiment::Merge(const Experiment *src)
{
    int nwarn = 0;

    for (auto &prop : src->property_)
    {
        auto mpo = prop.first;
        for (auto &gp : prop.second)
        {
            addProperty(mpo, gp);
        }
    }

    for (auto &cai : src->calcAtomConst())
    {
        double   x, y, z;
        CalcAtom caa(cai.getName(), cai.getObtype(), cai.getAtomid());

        cai.getCoords(&x, &y, &z);
        caa.SetCoords(x, y, z);
        caa.SetUnit(cai.getUnit());
        caa.SetResidue(cai.ResidueName(), cai.ResidueNumber());
        for (const auto &aci : cai.chargesConst())
        {
            caa.AddCharge(aci.first, aci.second);
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

void Experiment::getReferenceLot(std::string *reference,
                                 std::string *lot) const
{
    if (nullptr != reference)
    {
        reference->assign(getReference());
    }
    if (nullptr != lot)
    {
        if (DataSource() == dsExperiment)
        {
            lot->assign("Experiment");
        }
        else
        {
            lot->assign(gmx::formatString("%s/%s", method_.c_str(), basisset_.c_str()));
        }
    }
}

bool Experiment::getCharges(std::vector<double> *q,
                            qType                qtype,
                            std::string         *reference,
                            std::string         *lot) const
{
    q->resize(NAtom(), 0.0);
    int i = 0;
    for (auto &mai : calcAtomConst())
    {
        if (mai.hasCharge(qtype))
        {
            (*q)[i] = mai.charge(qtype);
            i++;
        }
    }
    getReferenceLot(reference, lot);
    return i == NAtom();
}

CommunicationStatus Experiment::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs;
    std::string         jobtype;

    cs = cr->recv_data(src);
    if (CS_OK == cs)
    {
        dataSource_ = static_cast<DataSource>(cr->recv_int(src));
        cr->recv_str(src, &reference_);
        cr->recv_str(src, &conformation_);
        cr->recv_str(src, &program_);
        cr->recv_str(src, &method_);
        cr->recv_str(src, &basisset_);
        cr->recv_str(src, &datafile_);
        cr->recv_str(src, &jobtype);
        jobtype_    = string2jobType(jobtype);

        int nmpo = cr->recv_int(src);
        
        //! Receive Properties
        for (int i = 0; i < nmpo; i++)
        {
            MolPropObservable mpo;
            std::string       mpo_str;
            cr->recv_str(src, &mpo_str);
            if (!stringToMolPropObservable(mpo_str, &mpo))
            {
                gmx_fatal(FARGS, "Unknown observable %s", mpo_str.c_str());
            }
            int  ngp = cr->recv_int(src);
            for (int n = 0; n < ngp; n++)
            {
                GenericProperty *gp;
                switch (mpo)
                {
                case MolPropObservable::DIPOLE:
                    {
                        gp = new MolecularDipole;
                        break;
                    }
                case MolPropObservable::QUADRUPOLE:
                    {
                       gp = new MolecularQuadrupole;
                       break;
                    }
                case MolPropObservable::POLARIZABILITY:
                    {
                        gp = new MolecularPolarizability;
                        break;
                    }
                case MolPropObservable::HF:
                case MolPropObservable::DHFORM:
                case MolPropObservable::DGFORM:
                case MolPropObservable::DSFORM:
                case MolPropObservable::ENTROPY:
                case MolPropObservable::STRANS:
                case MolPropObservable::SROT:
                case MolPropObservable::SVIB:
                case MolPropObservable::CP:
                case MolPropObservable::CV:
                case MolPropObservable::ZPE:
                case MolPropObservable::EMOL:
                    {
                        gp = new MolecularEnergy;
                        break;
                    }
                case MolPropObservable::POTENTIAL:
                case MolPropObservable::CHARGE:
                case MolPropObservable::COORDINATES:
                    {
                        gmx_fatal(FARGS, "Don't know what to do...");
                    }
                }
                gp->Receive(cr, src);
                addProperty(mpo, gp);

            }
        } 
        
        //! Receive Potentials
        int Npotential = cr->recv_int(src);
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
        int Natom = cr->recv_int(src);
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

CommunicationStatus Experiment::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs;
    std::string         jobtype;

    cs = cr->send_data(dest);
    if (CS_OK == cs)
    {
        cr->send_int(dest, static_cast<int>(dataSource_));
        cr->send_str(dest, &reference_);
        cr->send_str(dest, &conformation_);
        cr->send_str(dest, &program_);
        cr->send_str(dest, &method_);
        cr->send_str(dest, &basisset_);
        cr->send_str(dest, &datafile_);
        jobtype.assign(jobType2string(jobtype_));
        cr->send_str(dest, &jobtype);
        
        cr->send_int(dest, property_.size());
        for(const auto &prop : property_)
        {
            std::string mpo_str(mpo_name(prop.first));
            cr->send_str(dest, &mpo_str);
            cr->send_int(dest, prop.second.size());
            for (const auto &p : prop.second)
            {
                p->Send(cr, dest);
            }
        }
        
        //! Send Potentials
        if (CS_OK == cs)
        {
            cr->send_int(dest, electrostaticPotentialConst().size());
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
            cr->send_int(dest, calcAtomConst().size());
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

} // namespace alexandria
