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

#include "experiment.h"

#include "gromacs/utility/fatalerror.h"

#include "act/utility/units.h"

namespace alexandria
{

std::map<JobType, const char *> job_name =
{
    { JobType::OPT,      "Opt"      },
    { JobType::TOPOLOGY, "Topology" },
    { JobType::SP,       "SP"       },
    { JobType::UNKNOWN,  "unknown"  }
};

const char *jobType2string(JobType jType)

{
    return job_name[jType];
}

JobType string2jobType(const std::string &str)
{
    if (!str.empty())
    {
        for (const auto &s2j : job_name)
        {
            if (str.compare(s2j.second) == 0)
            {
                return s2j.first;
            }
        }
        auto buf = gmx::formatString("Invalid job type %s", str.c_str());
        GMX_THROW(gmx::InvalidInputError(buf.c_str()));
    }
    return JobType::UNKNOWN;
}

CalcAtomIterator Experiment::searchAtom(CalcAtom ca)
{
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

#ifdef EXPERIMENT_DESTRUCTOR
Experiment::~Experiment()
{
    for(auto &mp : property_)
    {
        for(auto pp = mp.second.begin(); pp < mp.second.end(); ++pp)
        {
            if (*pp)
            {
                delete *pp;
            }
        }
    }
}

Experiment::Experiment(const Experiment& other)// = default;
{
    auto oprops = other.propertiesConst();
    for(auto &op : oprops)
    {
        
        property_.insert({ op.first, {} });
        for(auto &gp : op.second)
        {
            property_[op.first].push_back(gp);
        }
    }
}

Experiment::Experiment(Experiment &&other)
{
    auto oprops = other.properties();
    for(auto &op : *oprops)
    {
        property_.insert({ op.first, {} });
        for(auto &gp : op.second)
        {
            property_[op.first].push_back(gp);
            gp = nullptr;
        }
    }
}

#endif

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
                double   x, y, z, fx, fy, fz;
                cai.coords(&x, &y, &z);
                cai.forces(&fx, &fy, &fz);
                fprintf(fp, "%-3s  %-3s  %3d X: %10.3f  %10.3f  %10.3f F: %10.3f  %10.3f  %10.3f\n",
                        cai.getName().c_str(), cai.getObtype().c_str(),
                        cai.getAtomid(), x, y, z, fx, fy, fz);
            }
        }
        for(const auto &p : property_)
        {
            fprintf(fp, "Property %s.\n", mpo_name(p.first));
            for (const auto &gp : p.second)
            {
                gp->Dump(fp);
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

        cai.coords(&x, &y, &z);
        caa.setCoordUnit(cai.coordUnit());
        caa.setCoords(x, y, z);
        cai.forces(&x, &y, &z);
        caa.setForceUnit(cai.forceUnit());
        caa.setForces(x, y, z);
        caa.setResidue(cai.residueName(), cai.residueNumber());
        for (const auto &aci : cai.chargesConst())
        {
            caa.AddCharge(aci.first, aci.second);
        }
        AddAtom(caa);
    }
    return nwarn;
}

void Experiment::setCoordinates()
{
    if (coordinates_.empty())
    {
        gmx::RVec x = { 0, 0, 0 };
        coordinates_.resize(catom_.size(), x);
        size_t j = 0;
        for(auto ca: catom_)
        {
            auto cunit      = ca.coordUnit();
            x[XX]           = convertToGromacs(ca.getX(), cunit);
            x[YY]           = convertToGromacs(ca.getY(), cunit);
            x[ZZ]           = convertToGromacs(ca.getZ(), cunit);
            coordinates_[j] = x;
            j++;
        }
    }
}

void Experiment::AddAtom(CalcAtom ca)
{
    CalcAtomIterator cai = searchAtom(ca);

    if (cai == catom_.end())
    {
        gmx::RVec x, f;
        auto      cunit = ca.coordUnit();
        x[XX]           = convertToGromacs(ca.getX(), cunit);
        x[YY]           = convertToGromacs(ca.getY(), cunit);
        x[ZZ]           = convertToGromacs(ca.getZ(), cunit);
        coordinates_.push_back(x);
        auto      funit = ca.forceUnit();
        if (!funit.empty())
        {
            ca.forces(&f[XX], &f[YY], &f[ZZ]);
            for (int m = 0; m < DIM; m++)
            {
                f[m] = convertToGromacs(f[m], funit);
            }
            forces_.push_back(f);
        }
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

CommunicationStatus Experiment::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;
    std::string         jobtype;

    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send(dest, static_cast<int>(dataSource_));
        cr->send(dest, reference_);
        cr->send(dest, conformation_);
        cr->send(dest, program_);
        cr->send(dest, method_);
        cr->send(dest, basisset_);
        cr->send(dest, datafile_);
        jobtype.assign(jobType2string(jobtype_));
        cr->send(dest, jobtype);
        cr->send(dest, property_.size());
        for(const auto &prop : property_)
        {
            std::string mpo_str(mpo_name(prop.first));
            cr->send(dest, mpo_str);
            cr->send(dest, prop.second.size());
            for (const auto &p : prop.second)
            {
                p->Send(cr, dest);
            }
        }
        
        //! Send Atoms
        if (CommunicationStatus::OK == cs)
        {
            cr->send(dest, calcAtomConst().size());
            for (auto &cai : calcAtomConst())
            {
                cs = cai.Send(cr, dest);
                if (CommunicationStatus::OK != cs)
                {
                    break;
                }
            }
        }
    }

    if ((CommunicationStatus::OK != cs) && (nullptr != debug))
    {
        fprintf(debug, "Trying to send Experiment, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Experiment::BroadCast(const CommunicationRecord *cr,
                                          int                        root,
                                          MPI_Comm                   comm)
{
    CommunicationStatus cs = cr->bcast_data(comm);
    if (CommunicationStatus::OK == cs)
    {
        int         dt = static_cast<int>(dataSource_);
        cr->bcast(&dt, comm);
        dataSource_ = static_cast<DataSource>(dt);
        cr->bcast(&reference_, comm);
        cr->bcast(&conformation_, comm);
        cr->bcast(&program_, comm);
        cr->bcast(&method_, comm);
        cr->bcast(&basisset_, comm);
        cr->bcast(&datafile_, comm);
        const char *jt = jobType2string(jobtype_);
        std::string jobtype(jt);
        cr->bcast(&jobtype, comm);
        jobtype_  = string2jobType(jobtype);
        // BroadCast number of properties
        size_t nprop = property_.size();
        cr->bcast(&nprop, comm);
        if (debug && cr->rank() == root)
        {
            for (auto &p : property_)
            {
                std::string mpo_str(mpo_name(p.first));
                fprintf(debug, "Will broadcast mpo_name '%s'\n", mpo_str.c_str());
            }
        }
        auto thisProp = property_.begin();
        for(size_t n = 0; n < nprop; n++)
        {
            size_t np = 0;
            if (root == cr->rank())
            {
                np = thisProp->second.size();
            }
            cr->bcast(&np, comm);
            if (0 == np)
            {
                continue;
            }
            std::string mpo_str;
            if (root == cr->rank())
            {
                mpo_str = mpo_name(thisProp->first);
            }
            cr->bcast(&mpo_str, comm);
            if (debug && root != cr->rank())
            {
                fprintf(debug, "Received broadcast string '%s' n = %zu\n", mpo_str.c_str(), n);
            }
            MolPropObservable mpo;
            if (!stringToMolPropObservable(mpo_str, &mpo))
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Received unknown string '%s' for a MolPropObservable", mpo_str.c_str()).c_str()));
            }
            for(size_t propid = 0; propid < np; propid++)
            {
                GenericProperty *gp = nullptr;
                if (root == cr->rank())
                {
                    gp = thisProp->second[propid];
                }
                else
                {
                    switch (mpo)
                    {
                    case MolPropObservable::DIPOLE:
                    case MolPropObservable::QUADRUPOLE:
                    case MolPropObservable::OCTUPOLE:
                    case MolPropObservable::HEXADECAPOLE:
                        {
                            gp = new MolecularMultipole;
                            break;
                        }
                    case MolPropObservable::POLARIZABILITY:
                        {
                            gp = new MolecularPolarizability;
                            break;
                    }
                    case MolPropObservable::FREQUENCY:
                    case MolPropObservable::INTENSITY:
                        {
                            gp = new Harmonics;
                            break;
                        }
                    case MolPropObservable::HF:
                    case MolPropObservable::DELTAE0:
                    case MolPropObservable::INTERACTIONENERGY:
                    case MolPropObservable::ELECTROSTATICS:
                    case MolPropObservable::EXCHANGE:
                    case MolPropObservable::DISPERSION:
                    case MolPropObservable::CHARGETRANSFER:
                    case MolPropObservable::INDUCTION:
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
                        {
                            gp = new MolecularEnergy;
                            break;
                        }
                    case MolPropObservable::POTENTIAL:
                        {
                            gp = new ElectrostaticPotential;
                            break;
                        }
                    case MolPropObservable::CHARGE:
                    case MolPropObservable::COORDINATES:
                        {
                            gmx_fatal(FARGS, "Don't know what to do...");
                        }
                    }
                }
                if (gp)
                {
                    gp->BroadCast(cr, root, comm);
                    if (root != cr->rank())
                    {
                        addProperty(mpo, gp);
                    }
                }
            }
            // Find next property to broadcast from the root.
            if (root == cr->rank() && thisProp != property_.end())
            {
                thisProp++;
            }
        }
        
        //! Broadcast Atoms
        int Natom = catom_.size();
        cr->bcast(&Natom, comm);
        for (int n = 0; (CommunicationStatus::OK == cs) && (n < Natom); n++)
        {
            if (cr->rank() == root)
            {
                catom_[n].BroadCast(cr, root, comm);
            }
            else
            {
                CalcAtom ca;
                cs = ca.BroadCast(cr, root, comm);
                if (CommunicationStatus::OK == cs)
                {
                    AddAtom(ca);
                }
            }
        }
    }

    if ((CommunicationStatus::OK != cs) && (nullptr != debug))
    {
        fprintf(debug, "Trying to receive Experiment, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Experiment::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    std::string         jobtype;

    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        int i;
        cr->recv(src, &i);
        dataSource_ = static_cast<DataSource>(i);
        cr->recv(src, &reference_);
        cr->recv(src, &conformation_);
        cr->recv(src, &program_);
        cr->recv(src, &method_);
        cr->recv(src, &basisset_);
        cr->recv(src, &datafile_);
        cr->recv(src, &jobtype);
        jobtype_    = string2jobType(jobtype);
        size_t nmpo;
        cr->recv(src, &nmpo);
        
        //! Receive Properties
        for (size_t i = 0; i < nmpo; i++)
        {
            MolPropObservable mpo;
            std::string       mpo_str;
            cr->recv(src, &mpo_str);
            if (!stringToMolPropObservable(mpo_str, &mpo))
            {
                gmx_fatal(FARGS, "Unknown observable %s", mpo_str.c_str());
            }
            size_t  ngp;
            cr->recv(src, &ngp);
            for (size_t n = 0; n < ngp; n++)
            {
                GenericProperty *gp = nullptr;
                switch (mpo)
                {
                case MolPropObservable::DIPOLE:
                case MolPropObservable::QUADRUPOLE:
                case MolPropObservable::OCTUPOLE:
                case MolPropObservable::HEXADECAPOLE:
                    {
                       gp = new MolecularMultipole;
                       break;
                    }
                case MolPropObservable::POLARIZABILITY:
                    {
                        gp = new MolecularPolarizability;
                        break;
                    }
                case MolPropObservable::FREQUENCY:
                case MolPropObservable::INTENSITY:
                    {
                        gp = new Harmonics;
                        break;
                    }
                case MolPropObservable::HF:
                case MolPropObservable::DELTAE0:
                case MolPropObservable::INTERACTIONENERGY:
                case MolPropObservable::ELECTROSTATICS:
                case MolPropObservable::EXCHANGE:
                case MolPropObservable::DISPERSION:
                case MolPropObservable::INDUCTION:
                case MolPropObservable::CHARGETRANSFER:
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
                    {
                        gp = new MolecularEnergy;
                        break;
                    }
                case MolPropObservable::POTENTIAL:
                    {
                        gp = new ElectrostaticPotential;
                        break;
                    }
                case MolPropObservable::CHARGE:
                case MolPropObservable::COORDINATES:
                    {
                        gmx_fatal(FARGS, "Don't know what to do...");
                    }
                }
                if (gp)
                {
                    gp->Receive(cr, src);
                    addProperty(mpo, gp);
                }

            }
        } 
        
        //! Receive Atoms
        size_t Natom;
        cr->recv(src, &Natom);
        for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Natom); n++)
        {
            CalcAtom ca;
            cs = ca.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                AddAtom(ca);
            }
        }
    }

    if ((CommunicationStatus::OK != cs) && (nullptr != debug))
    {
        fprintf(debug, "Trying to receive Experiment, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

} // namespace alexandria
