/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include "act/molprop/molprop.h"

#include <cmath>

#include <map>
#include <set>
#include <string>
#include <vector>

#include "act/molprop/composition.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/fatalerror.h"

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

void MolProp::generateComposition()
{
    for (const auto &ex : exper_)
    {
        // Assume that if an experiment has 1 atom, it has all of them
        if (ex.NAtom() > 0)
        {
            composition_.clear();
            bool bHasAtomTypes = true;
            for(const auto &ca : ex.calcAtomConst())
            {
                auto elem = ca.getName().c_str();
                auto cc   = composition_.find(elem);
                if (cc == composition_.end())
                {
                    composition_.insert(std::pair<const char *, int>(elem, 1));
                }
                else
                {
                    cc->second += 1;
                }
                bHasAtomTypes = bHasAtomTypes && !ca.getObtype().empty();
            }
            // Update the MolProp internal variables
            hasAllAtomTypes_ = bHasAtomTypes;
            natom_           = ex.NAtom(); 
        }
    }
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

bool MolProp::BondExists(const Bond &b)
{
    for (auto &bi : bondsConst())
    {
        if (((bi.aI() == b.aI()) && (bi.aJ() == b.aJ())) ||
            ((bi.aI() == b.aJ()) && (bi.aJ() == b.aI())))
        {
            return true;
        }
    }
    return false;
}

int MolProp::totalCharge()
{
    int qtot = 0;
    
    for(const auto &f : fragment_)
    {
        qtot += f.charge();
    }
    return qtot;
}

int MolProp::totalMultiplicity()
{
    int mtot = 1;
    
    for(const auto &f : fragment_)
    {
        int mm = f.multiplicity();
        if (mm % 2 == 0)
        {
            if (mtot == 1)
            {
                mtot = 2;
            }
            else
            {
                mtot = 1;
            }
        }
    }
    return mtot;
}

double MolProp::totalMass() const
{
    double mtot = 0;
    
    for(const auto &f : fragment_)
    {
        mtot += f.mass();
    }
    return mtot;
}

int MolProp::symmetryNumber() const
{
    int symm = 1;
    
    if (fragment_.size() == 1)
    {
        symm = fragment_[0].symmetryNumber();
    }
    return symm;
}

static std::string strNum2aa(const std::string &elem, int n)
{
    if (n > 1)
    {
        return gmx::formatString("%s%d", elem.c_str(), n);
    }
    else
    {
        return elem;
    }
}

static std::string comp2formula(std::map<std::string, int> &comp)
{
    std::string form;
    std::string C("C");
    std::string H("H");
    if (comp.end() != comp.find(C))
    {
        // Treat any compound with C as organic
        form += strNum2aa(C, comp[C]);
        if (comp.end() != comp.find(H))
        {
            form += strNum2aa(H, comp[H]);
        }
        // Other elements
        for(const auto &cc : comp)
        {
            if (cc.first != C && cc.first != H)
            {
                form += strNum2aa(cc.first, cc.second);
            }
        }
    }
    else
    {
        for(const auto &cc : comp)
        {
            form += strNum2aa(cc.first, cc.second);
        }
    }
    return form;
}

static bool fragCompare(const Fragment &a, const Fragment &b)
{
    int amin = *std::min_element(a.atoms().begin(), a.atoms().end());
    int bmin = *std::min_element(b.atoms().begin(), b.atoms().end());
    return amin < bmin;
}

bool MolProp::renumberResidues()
{
    Experiment *myexp = findExperiment(JobType::OPT);
    if (nullptr == myexp)
    {
        myexp = findExperiment(JobType::TOPOLOGY);
    }
    if (nullptr == myexp)
    {
        myexp = LastExperiment();
    }
    if (!myexp)
    {
        return false;
    }
    auto calcAtom = myexp->calcAtom();
    for(size_t j = 0; j < fragment_.size(); j++)
    {
        // TODO: This ignore the possibility that there could be
        // more than one residue in a fragment, e.g. a protein.
        auto atoms_j = fragment_[j].atoms();
        auto res_j   = (*calcAtom)[atoms_j[0]].residueName();
        for(int k : atoms_j)
        {
            (*calcAtom)[k].setResidue(res_j, j);
        }
    }
    return true;
}

void MolProp::generateFragments(const ForceField *pd,
                                double            qtotal)
{
    auto catom = LastExperiment()->calcAtomConst();
    int  natom = catom.size();
    std::vector<int> fragIndex(natom, -1);
    int maxBond = 0;
    for(auto b : bond_)
    {
        int ai = b.aI();
        int aj = b.aJ();
        if (fragIndex[ai] == -1 && fragIndex[aj] == -1)
        {
            fragIndex[ai] = maxBond;
            fragIndex[aj] = maxBond;
            maxBond += 1;
        }
        else if (fragIndex[ai] == -1)
        {
            fragIndex[ai] = fragIndex[aj];
        }
        else if (fragIndex[aj] == -1)
        {
            fragIndex[aj] = fragIndex[ai];
        }
        else if (fragIndex[ai] != fragIndex[aj])
        {
            // Merge them.
            if (fragIndex[ai] < fragIndex[aj])
            {
                int faj = fragIndex[aj];
                for(int i = 0; i < natom; i++)
                {
                    if (faj == fragIndex[i])
                    {
                        fragIndex[i] = fragIndex[ai];
                    }
                }
            }
            else
            {
                int fai = fragIndex[ai];
                for(int i = 0; i < natom; i++)
                {
                    if (fai == fragIndex[i])
                    {
                        fragIndex[i] = fragIndex[aj];
                    }
                }
            }
        }
    }
    for(int i = 0; i < natom; i++)
    {
        if (fragIndex[i] == -1)
        {
            fragIndex[i] = maxBond++;
        }
    }
    clearFragments();
    std::map<int, std::set<int>> fragMols;
    for(int i = 0; i < natom; i++)
    {
        int  fragI = fragIndex[i];
        auto fm    = fragMols.find(fragI);
        if (fragMols.end() == fm)
        {
            // Add empty set
            fragMols.insert({ fragI, { }});
        }
        // Add atom to set
        fragMols[fragI].insert(i);
    }
    if (debug)
    {
        fprintf(debug, "There are %zu distinct compounds in %d atoms and %zu bonds.\n", fragMols.size(),
                natom, bond_.size());
    }
    int  fi    = 0;
    for(auto &fm : fragMols)
    {
        std::vector<int> atomIndices;
        double           mass = 0;
        std::map<std::string, int> comp;
        for(auto &i : fm.second)
        {
            atomIndices.push_back(i);
            // Extract element and mass from force field data
            auto atype = catom[i].getObtype();
            if (pd->hasParticleType(atype))
            {
                auto ptype = pd->findParticleType(atype);
                mass += ptype->mass();
                auto elem  = ptype->element();
                if (comp.find(elem) == comp.end())
                {
                    comp.insert({elem, 0});
                }
                comp[elem] += 1;
            }
            else
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("No particle type %s in force field %s\n",
                                                                   atype.c_str(), pd->filename().c_str()).c_str()));;
            }
        }
        auto formula  = comp2formula(comp);
        auto fragName = gmx::formatString("%d", fi++);
        int  spin     = 1;
        addFragment(Fragment(fragName, mass, qtotal, spin,1, formula, atomIndices));
        // TODO we need proper charges for each fragment
        // Now setting this to zero for subsequent fragments.
        qtotal        = 0;
    }
    std::sort(fragment_.begin(), fragment_.end(), fragCompare);
    
    renumberResidues();
}

int MolProp::Merge(const MolProp *src)
{
    std::string stmp;
    int         nwarn = 0;

    for (auto &si : src->categoryConst())
    {
        AddCategory(si);
    }
    setIndex(src->getIndex());
    for(auto &src_f : src->fragments())
    {
        bool found = false;
        for (auto &dst_f : fragment_)
        {
            // TODO Make this check more rigorous, check for overlaps etc.
            if (//src_f.id() == dst_f.id() &&
                src_f.atoms() == dst_f.atoms() &&
                src_f.multiplicity() == dst_f.multiplicity() &&
                src_f.symmetryNumber() == dst_f.symmetryNumber() &&
                src_f.formula() == dst_f.formula() &&
                src_f.charge() == dst_f.charge())
            {
                found = true;
            }
        }
        if (!found)
        {
            fragment_.push_back(src_f);
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
        for (auto &bi : src->bondsConst())
        {
            alexandria::Bond bb(bi.aI(), bi.aJ(), bi.bondOrder());
            AddBond(bb);
        }
    }
    else
    {
        for (auto &bi : src->bondsConst())
        {
            alexandria::Bond bb(bi.aI(), bi.aJ(), bi.bondOrder());
            if (!BondExists(bb))
            {
                fprintf(stderr, "WARNING bond %d-%d not present in %s.\n",
                        bi.aI(), bi.aJ(), getMolname().c_str());
                for(auto &k : src->bondsConst())
                {
                    fprintf(stderr, "src bond %d %d\n",
                            k.aI(), k.aJ());
                }
                for(auto &k : bondsConst())
                {
                    fprintf(stderr, "dest bond %d %d\n",
                            k.aI(), k.aJ());
                }
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
            Experiment ca(ei.getProgram(), ei.getMethod(),
                          ei.getBasisset(), ei.getReference(),
                          ei.getConformation(), ei.getDatafile(),
                          jtype);
            nwarn += ca.Merge(&ei);
            AddExperiment(ca);
        }
    }
    return nwarn;
}

std::string MolProp::texFormula() const
{
    std::string texform;
    for(size_t i = 0; i < fragment_.size(); i++)
    {
        if (i > 0)
        {
            texform += ".";
        }
        texform += fragment_[i].texFormula();
    }
    return texform;
}

std::string MolProp::formula() const
{
    std::string form;
    for(size_t i = 0; i < fragment_.size(); i++)
    {
        if (i > 0)
        {
            form += ".";
        }
        form += fragment_[i].formula();
    }
    return form;
}

void MolProp::Dump(FILE *fp) const
{
    if (fp)
    {
        fprintf(fp, "molname:      %s\n", getMolname().c_str());
        fprintf(fp, "iupac:        %s\n", getIupac().c_str());
        fprintf(fp, "CAS:          %s\n", getCas().c_str());
        fprintf(fp, "cis:          %s\n", getCid().c_str());
        fprintf(fp, "InChi:        %s\n", getInchi().c_str());
        fprintf(fp, "category:    ");
        for (auto &f : fragments())
        {
            f.dump(fp);
        }
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

bool bCheckTemperature(double Tref, double T)
{
    return (Tref < 0) || (fabs(T - Tref) < 0.05);
}

const Experiment *MolProp::findExperimentConst(JobType job) const
{
    for(auto ei = exper_.begin(); ei < exper_.end(); ++ei)
    {
        if (ei->getJobtype() == job)
        {
            return &(*ei);
        }
    }
    return nullptr;
}

Experiment *MolProp::findExperiment(JobType job)
{
    for(auto ei = exper_.begin(); ei < exper_.end(); ++ei)
    {
        if (ei->getJobtype() == job)
        {
            return &(*ei);
        }
    }
    return nullptr;
}

const GenericProperty *MolProp::qmProperty(MolPropObservable  mpo, 
                                           double             T,
                                           JobType            jt) const
{
    for(auto ei = exper_.begin(); ei < exper_.end(); ++ei)
    {
        if (ei->hasProperty(mpo))
        {
            for (const auto &pp : ei->propertyConst(mpo))
            { 
                if (ei->getJobtype() == jt &&
                    bCheckTemperature(T, pp->getTemperature()))
                {
                    return pp;
                }
            }
        }
    }
    return nullptr;
}

const GenericProperty *MolProp::expProperty(MolPropObservable  mpo, 
                                            double             T) const
{
    for(auto ei = exper_.begin(); ei < exper_.end(); ++ei)
    {
        if (ei->hasProperty(mpo))
        {
            for (const auto &pp : ei->propertyConst(mpo))
            { 
                if (ei->dataSource() == dsExperiment &&
                    bCheckTemperature(T, pp->getTemperature()))
                {
                    return pp;
                }
            }
        }
    }
    return nullptr;
}

CommunicationStatus MolProp::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus                cs = CommunicationStatus::OK;

    /* Generic stuff */
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_int(dest, index_);
        cr->send_str(dest, &molname_);
        cr->send_str(dest, &iupac_);
        cr->send_str(dest, &cas_);
        cr->send_str(dest, &cid_);
        cr->send_str(dest, &inchi_);
        cr->send_int(dest, bond_.size());
        cr->send_int(dest, category_.size());
        cr->send_int(dest, exper_.size());
        cr->send_int(dest, fragment_.size());
        
        /* Send Bonds */
        for (auto &bi : bondsConst())
        {
            cs = bi.Send(cr, dest);
            if (CommunicationStatus::OK != cs)
            {
                break;
            }
        }

        /* Send Categories */
        if (CommunicationStatus::OK == cs)
        {
            for (auto &si : categoryConst())
            {
                if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
                {
                    std::string sii = si.c_str();
                    cr->send_str(dest, &sii);
                    if (nullptr != debug)
                    {
                        fprintf(debug, "Sent category %s\n", si.c_str());
                        fflush(debug);
                    }
                }
            }
        }
        
        /* Send Fragments */
        if (CommunicationStatus::OK == cs)
        {
            for(auto &f : fragments())
            {
                cs = f.Send(cr, dest);
                if (CommunicationStatus::OK != cs)
                {
                    break;
                }
            }
        }
             
        /* Send Experiments */
        if (CommunicationStatus::OK == cs)
        {
            for (auto &ei : experimentConst())
            {
                cs = ei.Send(cr, dest);
                if (CommunicationStatus::OK != cs)
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

CommunicationStatus MolProp::BroadCast(const CommunicationRecord *cr,
                                       int                        root,
                                       MPI_Comm                   comm)
{
    CommunicationStatus cs = cr->bcast_data(comm);

    /* Generic stuff */
    if (CommunicationStatus::OK == cs)
    {
        //! Receive index and more
        cr->bcast(&index_, comm);
        cr->bcast(&molname_, comm);
        cr->bcast(&iupac_, comm);
        cr->bcast(&cas_, comm);
        cr->bcast(&cid_, comm);
        cr->bcast(&inchi_, comm);
        int Nbond = bond_.size();
        cr->bcast(&Nbond, comm);
        if (cr->rank() == root)
        {
            for(auto &b : *bonds())
            {
                b.BroadCast(cr, root, comm);
            }
        }
        else
        {
            //! Receive Bonds
            cs = CommunicationStatus::OK;
            for (int n = 0; (CommunicationStatus::OK == cs) && (n < Nbond); n++)
            {
                Bond b;
                cs = b.BroadCast(cr, root, comm);
                if (CommunicationStatus::OK == cs)
                {
                    AddBond(b);
                }
            }
        }

        int Ncategory = category_.size();
        cr->bcast(&Ncategory, comm);
        //! Receive Categories
        for (int n = 0; (CommunicationStatus::OK == cs) && (n < Ncategory); n++)
        {
            if (CommunicationStatus::OK == cr->bcast_data(comm))
            {
                std::string str;
                cr->bcast(&str, comm);
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

        int Nfrag     = fragment_.size();
        cr->bcast(&Nfrag, comm);
        if (cr->rank() == root)
        {
            for(size_t ii = 0; ii < fragment_.size(); ii++)
            {
                fragment_[ii].BroadCast(cr, root, comm);
            }
        }
        else
        {
            //! Receive Fragments
            for (int n = 0; (CommunicationStatus::OK == cs) && (n < Nfrag); n++)
            {
                Fragment f;
                cs = f.BroadCast(cr, root, comm);
                if (CommunicationStatus::OK == cs)
                {
                    fragment_.push_back(std::move(f));
                }
            }
        }

        int Nexper    = exper_.size();
        cr->bcast(&Nexper, comm);
        //! Receive Experiments
        if (cr->rank() == root)
        {
            for (int n = 0; (CommunicationStatus::OK == cs) && (n < Nexper); n++)
            {
                exper_[n].BroadCast(cr, root, comm);
            }
            if (debug)
            {
                fprintf(debug, "Broadcast %zu experiments for mol %s\n",
                        exper_.size(), getMolname().c_str());
                fflush(debug);
            }
        }
        else
        {
            for (int n = 0; (CommunicationStatus::OK == cs) && (n < Nexper); n++)
            {
                Experiment ex;
                cs = ex.BroadCast(cr, root, comm);
                if (CommunicationStatus::OK == cs)
                {
                    AddExperiment(ex);
                }
            }
            if (nullptr != debug)
            {
                fprintf(debug, "Received %zu experiments from %d for mol %s\n",
                        exper_.size(), root, getMolname().c_str());
                fflush(debug);
            }
        }
    }
    return cs;
}

CommunicationStatus MolProp::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;

    /* Generic stuff */
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        //! Receive index and more
        index_        = cr->recv_int(src);
        cr->recv_str(src, &molname_);
        cr->recv_str(src, &iupac_);
        cr->recv_str(src, &cas_);
        cr->recv_str(src, &cid_);
        cr->recv_str(src, &inchi_);
        int Nbond     = cr->recv_int(src);
        int Ncategory = cr->recv_int(src);
        int Nexper    = cr->recv_int(src);
        int Nfrag     = cr->recv_int(src);

        if (nullptr != debug)
        {
            fprintf(debug, "Got molname %s\n", getMolname().c_str());
        }
        //! Receive Bonds
        cs = CommunicationStatus::OK;
        for (int n = 0; (CommunicationStatus::OK == cs) && (n < Nbond); n++)
        {
            Bond b;
            cs = b.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                AddBond(b);
            }
        }

        //! Receive Categories
        for (int n = 0; (CommunicationStatus::OK == cs) && (n < Ncategory); n++)
        {
            if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
            {
                std::string str;
                cr->recv_str(src, &str);
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

        //! Receive Fragments
        for (int n = 0; (CommunicationStatus::OK == cs) && (n < Nfrag); n++)
        {
            Fragment f;
            cs = f.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                fragment_.push_back(std::move(f));
            }
        }

        //! Receive Experiments
        for (int n = 0; (CommunicationStatus::OK == cs) && (n < Nexper); n++)
        {
            Experiment ex;
            cs = ex.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                AddExperiment(ex);
            }
        }

        if (nullptr != debug)
        {
            fprintf(debug, "Received %zu experiments from %d for mol %s\n",
                    exper_.size(), src, getMolname().c_str());
            fprintf(debug, "Received MolProp %s, status %s\n",
                    getMolname().c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

}
