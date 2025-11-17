/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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

#include "act/basics/msg_handler.h"
#include "act/molprop/composition.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/textwriter.h"

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

void MolProp::generateFragments(MsgHandler       *msg_handler,
                                const ForceField *pd)
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
    if (msg_handler && msg_handler->debug())
    {
        msg_handler->writeDebug(gmx::formatString("There are %zu distinct compounds in %d atoms and %zu bonds.\n",
                                                  fragMols.size(),
                                                  natom, bond_.size()));
    }
    int  fi    = 0;
    for(auto &fm : fragMols)
    {
        std::vector<int>           atomIndices;
        double                     mass = 0;
        std::map<std::string, int> comp;
        for(auto &i : fm.second)
        {
            atomIndices.push_back(i);
            // Extract element and mass from force field data
            auto atype = catom[i].getObtype();
            if (pd->hasParticleType(atype))
            {
                auto ptype = pd->findParticleType(atype);
                mass      += ptype->mass();
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
        // Set charge to 0 and spin to 1 now, fix it later.
        int qtot = 0;
        int spin = 1;
        addFragment(Fragment(fragName, mass, qtot, spin, 1, formula, atomIndices));
    }
    std::sort(fragment_.begin(), fragment_.end(), fragCompare);

    renumberResidues();
}

std::vector<std::string> MolProp::sameCompound(const MolProp *other)
{
    std::vector<std::string> warnings;
    if (other->getMolname() != getMolname())
    {
        warnings.push_back(gmx::formatString("Molnames differ. %s vs %s", getMolname().c_str(),
                                             other->getMolname().c_str()));
    }
    else if (other->fragments().size() != fragment_.size())
    {
        warnings.push_back(gmx::formatString("Fragments size %zu vs %zu", fragment_.size(),
                                             other->fragments().size()));
    }

    for(size_t ff = 0; warnings.empty() && ff < fragment_.size(); ++ff)
    {
        auto &src_f = fragment_[ff];
        auto &dst_f = other->fragments()[ff];

        if (src_f.atoms() != dst_f.atoms())
        {
            warnings.push_back(gmx::formatString("Atoms differ for %s", getMolname().c_str()));
        }
        else if (src_f.formula() != dst_f.formula())
        {
            warnings.push_back(gmx::formatString("Formula: %s vs %s for %s",
                                                 dst_f.formula().c_str(), src_f.formula().c_str(),
                                                 getMolname().c_str()));
        }
        else if (src_f.multiplicity() != dst_f.multiplicity())
        {
            warnings.push_back(gmx::formatString("Multiplicity: %d vs %d for %s",
                                                 dst_f.multiplicity(), src_f.multiplicity(),
                                                 getMolname().c_str()));
        }
        else if (src_f.symmetryNumber() != dst_f.symmetryNumber())
        {
            warnings.push_back(gmx::formatString("Symmetry number: %d vs %d for %s",
                                                 dst_f.symmetryNumber(), src_f.symmetryNumber(),
                                                 getMolname().c_str()));
        }
        else if (src_f.charge() != dst_f.charge())
        {
            warnings.push_back(gmx::formatString("Charge: %d vs %d for %s", dst_f.charge(), src_f.charge(),
                                                 getMolname().c_str()));
        }
    }
    return warnings;
}

std::vector<std::string> MolProp::Merge(MolProp *src)
{
    std::vector<std::string> warnings;
    std::string              stmp;

    if (fragment_.empty())
    {
        for(const auto &f : src->fragments())
        {
            fragment_.push_back(f);
        }
    }
    for (auto &si : src->categoryConst())
    {
        AddCategory(si);
    }
    setIndex(src->getIndex());

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
        std::copy(src->bondsConst().begin(), 
                  src->bondsConst().end(),
                  std::back_inserter(bond_));
    }
    else
    {
        for (auto &bi : src->bondsConst())
        {
            alexandria::Bond bb1(bi.aI(), bi.aJ(), bi.bondOrder());
            alexandria::Bond bb2(bi.aJ(), bi.aI(), bi.bondOrder());
            if (!BondExists(bb1) && !BondExists(bb2))
            {
                warnings.push_back(gmx::formatString("WARNING bond %d-%d not present in %s.",
                                                     bi.aI(), bi.aJ(), getMolname().c_str()));
                for(auto &k : src->bondsConst())
                {
                    warnings.push_back(gmx::formatString("src bond %d %d", k.aI(), k.aJ()));
                }
                for(auto &k : bondsConst())
                {
                    warnings.push_back(gmx::formatString("dest bond %d %d", k.aI(), k.aJ()));
                }
            }
        }
    }
    for(auto ee = src->experiment()->begin(); ee != src->experiment()->end(); ++ee)
    {
        exper_.push_back(std::move(*ee));
    }

    return warnings;
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

void MolProp::Dump(gmx::TextWriter *tw) const
{
    if (tw)
    {
        tw->writeStringFormatted("molname:      %s\n", getMolname().c_str());
        tw->writeStringFormatted("iupac:        %s\n", getIupac().c_str());
        tw->writeStringFormatted("CAS:          %s\n", getCas().c_str());
        tw->writeStringFormatted("cis:          %s\n", getCid().c_str());
        tw->writeStringFormatted("InChi:        %s\n", getInchi().c_str());
        tw->writeStringFormatted("category:    ");
        for (auto &f : fragments())
        {
            f.dump(tw);
        }
        for (auto &si : categoryConst())
        {
            tw->writeStringFormatted(" '%s'", si.c_str());
        }
        tw->writeStringFormatted("\n");
        for (auto &ei : experimentConst())
        {
            ei.Dump(tw);
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

bool MolProp::hasQMProperty(MolPropObservable mpo, 
                            double            T,
                            JobType           jt) const
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
                    return true;
                }
            }
        }
    }
    return false;
}

const std::unique_ptr<GenericProperty> &MolProp::qmProperty(MolPropObservable  mpo, 
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
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such QM property %s", mpo_name(mpo))));
}

const std::unique_ptr<GenericProperty> &MolProp::expProperty(MolPropObservable  mpo, 
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
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such Experimental property %s", mpo_name(mpo))));
}

CommunicationStatus MolProp::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;

    /* Generic stuff */
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send(dest, index_);
        cr->send(dest, molname_);
        cr->send(dest, iupac_);
        cr->send(dest, cas_);
        cr->send(dest, cid_);
        cr->send(dest, inchi_);
        cr->send(dest, bond_.size());
        cr->send(dest, category_.size());
        cr->send(dest, exper_.size());
        cr->send(dest, fragment_.size());
        
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
                    std::string sii(si);
                    cr->send(dest, sii);
                    if (cr->mh() && cr->mh()->debug())
                    {
                        cr->mh()->writeDebug(gmx::formatString("Sent category %s\n", si.c_str()));
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
        if (cr->mh() && cr->mh()->debug())
        {
            cr->mh()->writeDebug(gmx::formatString("Sent MolProp %s, status %s\n",
                                                   getMolname().c_str(), cs_name(cs)));
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
        size_t Nbond = bond_.size();
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
            for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Nbond); n++)
            {
                Bond b;
                cs = b.BroadCast(cr, root, comm);
                if (CommunicationStatus::OK == cs)
                {
                    AddBond(b);
                }
            }
        }

        size_t Ncategory = category_.size();
        cr->bcast(&Ncategory, comm);
        //! Receive Categories
        for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Ncategory); n++)
        {
            if (CommunicationStatus::OK == cr->bcast_data(comm))
            {
                std::string str;
                cr->bcast(&str, comm);
                if (!str.empty())
                {
                    AddCategory(str);
                    if (cr->mh() && cr->mh()->debug())
                    {
                        cr->mh()->writeDebug(gmx::formatString("Received a category %s\n", str.c_str()));
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "A category was promised but I got a nullptr pointer");
                }
            }
        }

        size_t Nfrag     = fragment_.size();
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
            for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Nfrag); n++)
            {
                Fragment f;
                cs = f.BroadCast(cr, root, comm);
                if (CommunicationStatus::OK == cs)
                {
                    fragment_.push_back(std::move(f));
                }
            }
        }

        size_t Nexper    = exper_.size();
        cr->bcast(&Nexper, comm);
        //! Receive Experiments
        if (cr->rank() == root)
        {
            for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Nexper); n++)
            {
                exper_[n].BroadCast(cr, root, comm);
            }
            if (cr->mh() && cr->mh()->debug())
            {
                cr->mh()->writeDebug(gmx::formatString("Broadcast %zu experiments for mol %s\n",
                                                       exper_.size(), getMolname().c_str()));
            }
        }
        else
        {
            for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Nexper); n++)
            {
                Experiment ex;
                cs = ex.BroadCast(cr, root, comm);
                if (CommunicationStatus::OK == cs)
                {
                    AddExperiment(std::move(ex));
                }
            }
            if (cr->mh() && cr->mh()->debug())
            {
                cr->mh()->writeDebug(gmx::formatString("Received %zu experiments from %d for mol %s\n",
                                                       exper_.size(), root, getMolname().c_str()));
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
        cr->recv(src, &index_);
        cr->recv(src, &molname_);
        cr->recv(src, &iupac_);
        cr->recv(src, &cas_);
        cr->recv(src, &cid_);
        cr->recv(src, &inchi_);
        size_t Nbond, Ncategory, Nexper, Nfrag;
        cr->recv(src, &Nbond);
        cr->recv(src, &Ncategory);
        cr->recv(src, &Nexper);
        cr->recv(src, &Nfrag);

        if (cr->mh() && cr->mh()->debug())
        {
            cr->mh()->writeDebug(gmx::formatString("Got molname %s\n", getMolname().c_str()));
        }
        //! Receive Bonds
        cs = CommunicationStatus::OK;
        for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Nbond); n++)
        {
            Bond b;
            cs = b.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                AddBond(b);
            }
        }

        //! Receive Categories
        for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Ncategory); n++)
        {
            if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
            {
                std::string str;
                cr->recv(src, &str);
                if (!str.empty())
                {
                    AddCategory(str);
                    if (cr->mh() && cr->mh()->debug())
                    {
                        cr->mh()->writeDebug(gmx::formatString("Received a category %s\n", str.c_str()));
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "A category was promised but I got a nullptr pointer");
                }
            }
        }

        //! Receive Fragments
        for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Nfrag); n++)
        {
            Fragment f;
            cs = f.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                fragment_.push_back(std::move(f));
            }
        }

        //! Receive Experiments
        for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Nexper); n++)
        {
            Experiment ex;
            cs = ex.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                AddExperiment(std::move(ex));
            }
        }

        if (cr->mh() && cr->mh()->debug())
        {
            cr->mh()->writeDebug(gmx::formatString("Received %zu experiments from %d for mol %s\n",
                                                   exper_.size(), src, getMolname().c_str()));
            cr->mh()->writeDebug(gmx::formatString("Received MolProp %s, status %s\n",
                                                   getMolname().c_str(), cs_name(cs)));
        }
    }
    return cs;
}

}
