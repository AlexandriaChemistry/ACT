/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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

#include <map>
#include <string>
#include <vector>

#include "gromacs/math/utilities.h"
#include "gromacs/utility/fatalerror.h"

#include "communication.h"
#include "composition.h"
#include "gmx_simple_comm.h"
#include "units.h"

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



CommunicationStatus Bond::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, ai_);
        gmx_send_int(cr, dest, aj_);
        gmx_send_double(cr, dest, bondorder_);
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
        bondorder_ = gmx_recv_double(cr, src);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive Bond, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

bool Bond::operator==(const Bond &other) const
{
    return ((ai_ == other.getAi() && aj_ == other.getAj()) ||
            (aj_ == other.getAi() && ai_ == other.getAj()));
}

double MolProp::bondToBondOrder(int ai, int aj) const
{
    Bond   mybond(ai, aj, 1.0);
    auto   bb = std::find(bond_.begin(), bond_.end(), mybond);
    if (bb != bond_.end())
    {
        return bb->getBondOrder();
    }
    else
    {
        Bond   mybond(aj, ai, 1.0);
        auto   bb = std::find(bond_.begin(), bond_.end(), mybond);
        if (bb != bond_.end())
        {
            return bb->getBondOrder();
        }
    }
    printf("Looking for %d-%d\n", ai, aj);
    for(auto &b : bond_)
    {
        printf("Have %d-%d\n", b.getAi(), b.getAj());
    }
    gmx_fatal(FARGS, "Sorry...");
    return 0.0;
}

void MolProp::AddBond(Bond b)
{
    BondConstIterator bi;
    bool              bFound = false;

    for (bi = bondsConst().begin(); bi < bondsConst().end(); ++bi)
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

bool MolProp::BondExists(Bond b)
{
    for (auto &bi : bondsConst())
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
    q = totalCharge();
    if (q == 0)
    {
        SetTotalCharge(src->totalCharge());
    }
    else
    {
        sq = src->totalCharge();
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
        for (auto &bi : src->bondsConst())
        {
            alexandria::Bond bb(bi.getAi(), bi.getAj(), bi.getBondOrder());
            AddBond(bb);
        }
    }
    else
    {
        for (auto &bi : src->bondsConst())
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
        fprintf(fp, "charge:       %d\n", totalCharge());
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

    ncomp.resize(110, 0);
    myform[0]  = '\0';
    texform[0] = '\0';
    for(const auto &cc : composition_)
    { 
        real value;
        if (gmx_atomprop_query(ap, epropElement, "???", cc.first, &value))
        {
            int an = std::lround(value);
            range_check(an, 0, 110);
            if (an > 0)
            {
                ncomp[an] += cc.second;
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

const Experiment *MolProp::findExperimentConst(const std::string &method,
                                               const std::string &basis,
                                               const std::string &conformation) const
{
    for(auto ei = exper_.begin(); ei < exper_.end(); ++ei)
    {
        if (((conformation.size() == 0)   || stringEqual(ei->getConformation(), conformation)) &&
            stringEqual(ei->getMethod(), method) &&
            stringEqual(ei->getBasisset(), basis))
        {
            return &(*ei);
        }
    }
    return nullptr;
}

const GenericProperty *MolProp::findProperty(MolPropObservable  mpo, 
                                             iqmType            iQM,
                                             double             T,
                                             const std::string &method,
                                             const std::string &basis,
                                             const std::string &conf) const
{
    for(auto ei = exper_.begin(); ei < exper_.end(); ++ei)
    {
        if (ei->hasProperty(mpo))
        {
            for (const auto &pp : ei->propertyConst(mpo))
            { 
                if (bCheckTemperature(T, pp->getTemperature()))
                {
                    if ((ei->dataSource() == dsExperiment) &&
                        (iqmType::Exp == iQM || iqmType::Both == iQM))
                    {
                        // For an experiment it is sufficient if the temperature and
                        // the mpo are correct
                        return pp;
                    }
                    else if ((ei->dataSource() == dsTheory) &&
                             (iqmType::QM == iQM || iqmType::Both == iQM))
                    {
                        if (((conf.size() == 0)   || stringEqual(ei->getConformation(), conf)) &&
                            ((method.size() == 0) || stringEqual(ei->getMethod(), method)) &&
                            ((basis.size() == 0)  || stringEqual(ei->getBasisset(), basis)))
                        {
                            return pp;
                        }
                    }
                }
            }
        }
    }
    return nullptr;
}

#ifdef OLD
bool MolProp::getPropRef(MolPropObservable mpo, iqmType iQM,
                         const std::string &method,
                         const std::string &basis,
                         const std::string &conf,
                         double *value, double *error, double *T,
                         std::string *ref, std::string *mylot,
                         std::vector<double> *vec, tensor quad_polar)
{
    bool   done = false;
    double Told = *T;

    if (iQM == iqmType::Both)
    {
        for (auto &ei : experimentConst())
        {
            if ((conf.size() == 0) ||
                stringEqual(ei.getConformation(), conf))
            {
                if (ei.getVal(mpo, value, error, T, vec, quad_polar) &&
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
    if (iQM == iqmType::Exp)
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
    else if (iQM == iqmType::QM)
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
#endif
bool MolProp::getOptHF(double *value)
{
    bool done = false;

    std::string empty;
    
    auto gp = findProperty(MolPropObservable::HF, iqmType::QM, 0.0,
                           empty, empty, empty);
    if (gp)
    {
        *value = gp->getValue();
        done = true;
    }
    return done;
}

int MolProp::NOptSP()
{
    int n = 0;

    for (auto &ei : experimentConst())
    {
        if (ei.getJobtype() == JobType::OPT ||
            ei.getJobtype() == JobType::SP)
        {
            n++;
        }
    }
    return n;
}

#ifdef OLD
bool MolProp::getProp(MolPropObservable mpo, iqmType iQM,
                      const std::string &method,
                      const std::string &basis,
                      const std::string &conf,
                      const std::string &type,
                      double *value, double *error, double *T)
{
    double              myerror;
    std::vector<double> vec;
    tensor              quad;
    bool                bReturn;
    std::string         myref, mylot;

    bReturn = getPropRef(mpo, iQM, method, basis,
                         conf, type, value, &myerror, T,
                         &myref, &mylot, &vec, quad);
    if (nullptr != error)
    {
        *error = myerror;
    }
    return bReturn;
}
#endif

#ifdef OLD
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
                case MolPropObservable::POTENTIAL:
                    done = ci->NPotential() > 0;
                    break;
                case MolPropObservable::DIPOLE:
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
                case MolPropObservable::QUADRUPOLE:
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
                case MolPropObservable::POLARIZABILITY:
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
                case MolPropObservable::DHFORM:
                case MolPropObservable::ENTROPY:
                    for (auto &mdp : ci->molecularEnergyConst())
                    {
                        done = mdp.mpo() == mpo;
                        if (done)
                        {
                            break;
                        }
                    }
                    break;
                case MolPropObservable::COORDINATES:
                    done = NAtom() > 0;
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
#endif
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
        gmx_send_int(cr, dest, category_.size());
        gmx_send_int(cr, dest, exper_.size());

        /* Send Bonds */
        for (auto &bi : bondsConst())
        {
            cs = bi.Send(cr, dest);
            if (CS_OK != cs)
            {
                break;
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
    int                 Nbond, Ncategory, Nexper;

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
