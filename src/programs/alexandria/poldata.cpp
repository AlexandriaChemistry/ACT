/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019 
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

#include "poldata.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"

#include "gmx_simple_comm.h"
#include "stringutil.h"

namespace alexandria
{

void Poldata::setFilename(const std::string &fn2)
{
    GMX_RELEASE_ASSERT((fn2.size() > 0),
                       "Trying to set empty Poldata filename");
    if (0 != filename_.size())
    {
        fprintf(stderr, "Changing Poldata filename_ from %s to %s\n",
                filename_.c_str(), fn2.c_str());

    }
    filename_ = fn2;
}

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Atom STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

void Poldata::addAtype(const std::string &elem,
                       const std::string &desc,
                       const std::string &atype,
                       const std::string &ptype,
                       const std::string &btype,
                       const std::string &ztype,
                       bool               fixVdw,
                       std::string       &vdwparams,
                       const std::string &refEnthalpy)
{

    size_t i;
    for (i = 0; i < alexandria_.size(); i++)
    {
        if (alexandria_[i].getType().compare(atype) == 0)
        {
            break;
        }
    }
    if (i == alexandria_.size())
    {
        Ffatype sp(desc, atype, ptype, btype, ztype,
                   elem, fixVdw, vdwparams, refEnthalpy);

        alexandria_.push_back(sp);
    }
    else
    {
        fprintf(stderr, "Atom type %s was already added to Poldata record\n", atype.c_str());
    }
}

gmx_bool Poldata::strcasestrStart(std::string needle, std::string haystack)
{
    std::string ptr;
    ptr = strcasestr(haystack.c_str(), needle.c_str());
    return (ptr == haystack);
}

bool Poldata::getAtypeRefEnthalpy(const std::string &atype,
                                  double            *Href) const
{
    auto fa = findAtype(atype);
    if (alexandria_.end() != fa)
    {
        *Href = atof(fa->getRefEnthalpy().c_str());
        return true;
    }
    return false;
}

const std::string &Poldata::getDesc(const std::string &atype) const
{
    size_t i;
    if (atype.size() != 0)
    {
        for (i = 0; i < alexandria_.size(); i++)
        {
            if (alexandria_[i].getType().compare(atype) == 0)
            {
                return alexandria_[i].getDesc();
            }
        }
    }
    gmx_fatal(FARGS, "No such atomtype %s", atype.c_str());
}

const std::string &Poldata::getElem(const std::string &atype) const
{
    size_t i;
    if (atype.size() != 0)
    {
        for (i = 0; i < alexandria_.size(); i++)
        {
            if (alexandria_[i].getType().compare(atype) == 0)
            {
                return alexandria_[i].getElem();
            }
        }
    }
    gmx_fatal(FARGS, "No such atomtype %s", atype.c_str());
}

const std::string &Poldata::ztype2elem(const std::string &ztype) const
{
    size_t i;
    if (ztype.size() != 0)
    {
        for (i = 0; i < alexandria_.size(); i++)
        {
            if (alexandria_[i].getZtype().compare(ztype) == 0)
            {
                return alexandria_[i].getElem();
            }
        }
    }
    gmx_fatal(FARGS, "No such zeta type %s", ztype.c_str());
}

std::vector<std::string> Poldata::ztype_names() const
{
    std::vector<std::string> ztype_names;   
    for (auto atpi = alexandria_.begin(); atpi != alexandria_.end(); atpi++)
    { 
        if (std::find(ztype_names.begin(), ztype_names.end(), atpi->getZtype()) == ztype_names.end()) 
        {
            ztype_names.push_back(atpi->getZtype());
        }
    }
    return ztype_names;
}

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Polarizability STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

bool Poldata::atypeToPtype(const std::string &atype,
                           std::string       &ptype) const
{
    if (atype.size() == 0)
    {
        return false;
    }
    auto ai = std::find_if(alexandria_.begin(), alexandria_.end(),
                           [atype](Ffatype const &fa)
                           {
                               return fa.getType().compare(atype) == 0;
                           });
    if (ai != alexandria_.end() && ai->getPtype().size() > 0)
    {
        ptype = ai->getPtype();
        return true;
    }
    return false;
}

bool Poldata::getAtypePol(const std::string &atype,
                          double            *polar,
                          double            *sigPol) const
{
    auto fa = findAtype(atype);
    if (alexandria_.end() != fa)
    {
        return getPtypePol(fa->getPtype(), polar, sigPol);
    }
    return false;
}

bool Poldata::getPtypePol(const std::string &ptype,
                          double            *polar,
                          double            *sigPol) const
{
    size_t j;
    for (j = 0; j < ptype_.size(); j++)
    {
        if (ptype.compare(ptype_[j].getType()) == 0)
        {
            *polar   = ptype_[j].getPolarizability();
            *sigPol  = ptype_[j].getSigPol();
            return true;
        }
    }
    return false;
}

bool Poldata::setPtypePolarizability(const std::string &ptype,
                                     double             polarizability,
                                     double             sigPol)
{
    auto sp = findPtype(ptype);
    if (ptype_.end() != sp)
    {
        sp->setPolarizability(polarizability);
        sp->setSigPol(sigPol);
        return true;
    }
    return false;
}

void Poldata::addPtype(const std::string &ptype,
                       const std::string &miller,
                       const std::string &bosque,
                       double             polarizability,
                       double             sigPol)
{
    size_t i;
    for (i = 0; i < ptype_.size(); i++)
    {
        if (ptype_[i].getType().compare(ptype) == 0)
        {
            break;
        }
    }
    if (i == ptype_.size())
    {
        Ptype sp(ptype, miller, bosque, polarizability, sigPol);
        ptype_.push_back(sp);
    }
    else
    {
        fprintf(stderr, "Polarizability type %s was already added to Poldata record\n", ptype.c_str());
    }
}

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Virtual Site STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

void  Poldata::addVsite(const std::string &atype,
                        const std::string &type,
                        int                number,
                        double             distance,
                        double             angle,
                        int                ncontrolatoms)
{
    size_t i;
    for (i = 0; i < vsite_.size(); i++)
    {
        if (vsite_[i].atype() == atype &&
            vsite_[i].type()  == string2vsiteType(atype.c_str()))
        {
            break;
        }
    }
    if (i == vsite_.size())
    {
        Vsite vs(atype, type, number, distance, angle, ncontrolatoms);
        vsite_.push_back(vs);
    }
    else
    {
        fprintf(stderr, "vsite type %s was already added to Poldata record\n", atype.c_str());
    }
}

/*
 *-+-+-+-+-+-+-+-+
 * Bosque STUFF
 *-+-+-+-+-+-+-+-+
 */

bool Poldata::getBosquePol(const std::string &bosque,
                           double            *polarizability) const
{
    auto bb  = bosque_.begin(), be = bosque_.end();
    auto bos = std::find_if(bb, be, [bosque](Bosque const &b)
                            {
                                return (bosque.compare(b.getBosque()) == 0);
                            });
    if (bosque_.end() != bos)
    {
        *polarizability = bos->getPolarizability();

        return true;
    }
    return false;
}

bool Poldata::ptypeToBosque(const std::string &ptype,
                            std::string       &bosque) const
{
    for (const auto &i : ptype_)
    {
        if (ptype.compare(i.getType()) == 0)
        {
            bosque = i.getBosque();
            return true;
        }
    }
    return false;
}

/*
 *-+-+-+-+-+-+-+-+
 * Miller STUFF
 *-+-+-+-+-+-+-+-+
 */

bool Poldata::ptypeToMiller(const std::string &ptype,
                            std::string       &miller) const
{
    for (const auto &i : ptype_)
    {
        if (ptype.compare(i.getType()) == 0)
        {
            miller = i.getMiller();
            return true;
        }
    }
    return false;
}

void Poldata::addMiller(const std::string &miller,
                        int                atomnumber,
                        double             tauAhc,
                        double             alphaAhp,
                        const std::string &alexandria_equiv)
{
    Miller mil(miller, atomnumber, tauAhc, alphaAhp, alexandria_equiv);
    miller_.push_back(mil);
}


bool Poldata::getMillerPol(const std::string &miller,
                           int               *atomnumber,
                           double            *tauAhc,
                           double            *alphaAhp,
                           std::string       &alexandria_equiv) const
{
    auto mb  = miller_.begin(), me = miller_.end();
    auto mil = std::find_if(mb, me, [miller](Miller const &m)
                            {
                                return (miller.compare(m.getMiller()) == 0);
                            });
    if (miller_.end() != mil)
    {
        *atomnumber      = mil->getAtomnumber();
        *tauAhc          = mil->getTauAhc();
        *alphaAhp        = mil->getAlphaAhp();
        alexandria_equiv = mil->getAlexandriaEquiv();

        return true;
    }
    return false;
}

/*
 *-+-+-+-+-+-+-+-+-+-+
 * LISTED FORCES
 *-+-+-+-+-+-+-+-+-+-+
 */

bool Poldata::atypeToBtype(const std::string &atype,
                           std::string       &btype) const
{
    if (atype.size() == 0)
    {
        return false;
    }
    auto ai = std::find_if(alexandria_.begin(), alexandria_.end(),
                           [atype](Ffatype const &fa)
                           {
                               return fa.getType().compare(atype) == 0;
                           });
    if (ai != alexandria_.end() && ai->getBtype().size() > 0)
    {
        btype = ai->getBtype();
        return true;
    }
    return false;
}

void Poldata::addBtype(const std::string &btype)
{
    size_t i;
    for (i = 0; i < btype_.size(); i++)
    {
        if (btype.compare(btype_[i]) == 0)
        {
            break;
        }
    }
    if (i == btype_.size())
    {
        btype_.push_back(btype);
    }
}


bool Poldata::findForce(std::vector<std::string> &atoms,
                        ListedForceIterator      *force)
{
    for (auto &f : forces_)
    {
        auto tmp = f.findForce(atoms);
        if (f.forceEnd() != tmp)
        {
            *force = tmp;
            return true;
        }
    }
    return false;
}

bool Poldata::findForce(const std::vector<std::string> &atoms,
                        ListedForceConstIterator       *force) const
{
    for (const auto &f : forces_)
    {
        auto tmp = f.findForce(atoms);
        if (f.forceEnd() != tmp)
        {
            *force = tmp;
            return true;
        }
    }
    return false;
}

bool Poldata::findForce(std::vector<std::string> &atoms,
                        ListedForceIterator      *force,
                        size_t                    bondOrder)
{
    for (auto &f : forces_)
    {
        auto tmp = f.findForce(atoms, bondOrder);
        if (f.forceEnd() != tmp)
        {
            *force = tmp;
            return true;
        }
    }
    return false;
}

bool Poldata::findForce(const std::vector<std::string> &atoms,
                        ListedForceConstIterator       *force,
                        size_t                          bondOrder) const
{
    for (const auto &f : forces_)
    {
        auto tmp = f.findForce(atoms, bondOrder);
        if (f.forceEnd() != tmp)
        {
            *force = tmp;
            return true;
        }
    }
    return false;
}

bool Poldata::searchForce(std::vector<std::string> &atoms,
                          std::string              &params,
                          double                   *refValue,
                          double                   *sigma,
                          size_t                   *ntrain) const
{
    for (auto &f : forces_)
    {
        if (f.searchForce(atoms, params, refValue,
                          sigma, ntrain))
        {
            return true;
        }
    }
    return false;
}

bool Poldata::searchForce(std::vector<std::string> &atoms,
                          std::string              &params,
                          double                   *refValue,
                          double                   *sigma,
                          size_t                   *ntrain,
                          InteractionType           iType) const
{
    auto f  = findForces(iType);
    if (forcesEnd() != f)
    {
        if (f->searchForce(atoms, params, refValue,
                           sigma, ntrain))
        {
            return true;
        }
    }
    return false;
}

bool Poldata::searchForce(std::vector<std::string> &atoms,
                          std::string              &params,
                          double                   *refValue,
                          double                   *sigma,
                          size_t                   *ntrain,
                          size_t                    bondOrder) const
{
    for (auto &f : forces_)
    {
        if (f.searchForce(atoms, params, refValue,
                          sigma, ntrain, bondOrder))
        {
            return true;
        }
    }
    return false;
}

bool Poldata::searchForce(std::vector<std::string> &atoms,
                          std::string              &params,
                          double                   *refValue,
                          double                   *sigma,
                          size_t                   *ntrain,
                          size_t                    bondOrder,
                          InteractionType           iType) const
{
    auto f  = findForces(iType);
    if (forcesEnd() != f)
    {
        if (f->searchForce(atoms, params, refValue,
                           sigma, ntrain, bondOrder))
        {
            return true;
        }
    }
    return false;
}


/*
 *-+-+-+-+-+-+-+-+
 * VDW STUFF
 *-+-+-+-+-+-+-+-+
 */

void Poldata::setVdwFunction(const std::string &func)
{
    size_t i;
    for (i = 0; i < F_NRE; i++)
    {
        if (strcasecmp(interaction_function[i].name, func.c_str()) == 0)
        {
            break;
        }
    }
    if (i == F_NRE)
    {
        gmx_fatal(FARGS, "Van der Waals function '%s' does not exist in gromacs", func.c_str());
    }
    gtVdwFtype_    = i;
    gtVdwFunction_ = func;
}

void Poldata::setCombinationRule(const std::string &func)
{
    size_t i;
    for (i = 0; i < eCOMB_NR; i++)
    {
        if (strcasecmp(ecomb_names[i], func.c_str()) == 0)
        {
            break;
        }
    }
    if (i == eCOMB_NR)
    {
        gmx_fatal(FARGS, "Combination rule '%s' does not exist in gromacs", func.c_str());
    }
    gtCombRule_        = i;
    gtCombinationRule_ = func;
}

/*
 *-+-+-+-+-+-+-+-+
 * COULOMB STUFF
 *-+-+-+-+-+-+-+-+
 */

void Poldata::addSymcharges(const std::string &central,
                            const std::string &attached,
                            int                numattach)
{
    size_t        i;
    Symcharges  * sc;
    for (i = 0; i < symcharges_.size(); i++)
    {
        sc = &(symcharges_[i]);
        if ((strcasecmp(sc->getCentral().c_str(), central.c_str()) == 0)   &&
            (strcasecmp(sc->getAttached().c_str(), attached.c_str()) == 0) &&
            (sc->getNumattach() == numattach))
        {
            break;
        }
    }
    if (i == symcharges_.size())
    {
        Symcharges symcharges(central, attached, numattach);
        symcharges_.push_back(symcharges);
    }
}

int Poldata::havePolSupport(const std::string &atype) const
{
    size_t i;
    for (i = 0; i < alexandria_.size(); i++)
    {
        if (atype.compare(alexandria_[i].getType()) == 0)
        {
            return 1;
        }
    }
    return 0;
}

bool Poldata::haveEemSupport(const std::string &atype,
                             gmx_bool           bAllowZeroParameters) const
{
    auto eep = findEem(atype);
    return (eep != EndEemprops() &&
            (bAllowZeroParameters || ((eep->getJ0() > 0) && (eep->getChi0() > 0))));
}

double Poldata::getJ00(const std::string &atype) const
{
    EempropsConstIterator eer;
    if ((eer = findEem(atype)) != EndEemprops())
    {
        return eer->getJ0();
    }
    return -1;
}

const char *Poldata::getQstr(const std::string       &atype) const
{
    EempropsConstIterator eer;
    if ((eer = findEem(atype)) != EndEemprops())
    {
        return eer->getQstr();
    }
    return nullptr;
}

const char *Poldata::getRowstr(const std::string &name) const
{
    EempropsConstIterator eer;
    if ((eer = findEem(name)) != EndEemprops())
    {
        return eer->getRowstr();
    }
    return nullptr;
}

int Poldata::getRow(const std::string       &name,
                    int                      zz) const
{
    EempropsConstIterator eer;
    if ((eer = findEem(name)) != EndEemprops())
    {
        range_check(zz, 0, eer->getNzeta());
        return eer->getRow(zz);
    }
    return -1;
}

double Poldata::getZeta(const std::string       &name,
                        int                      zz) const
{
    EempropsConstIterator eer;
    if ((eer = findEem(name)) != EndEemprops())
    {
        if ((zz < 0) || (zz >= eer->getNzeta()))
        {
            printf("Negative zeta\n");
        }
        range_check(zz, 0, eer->getNzeta());
        return eer->getZeta(zz);
    }
    return -1;
}

int Poldata::getNzeta(const std::string &atype) const
{
    EempropsConstIterator eer;
    if ((eer = findEem(atype)) != EndEemprops())
    {
        return eer->getNzeta();
    }
    return 0;
}

double Poldata::getQ(const std::string       &atype,
                     int                      zz) const
{
    EempropsConstIterator eer;
    if ((eer = findEem(atype)) != EndEemprops())
    {
        range_check(zz, 0, eer->getNzeta());
        return eer->getQ(zz);
    }
    return -1;
}

double Poldata::getChi0(const  std::string &atype) const
{
    EempropsConstIterator eer;
    if ((eer = findEem(atype)) != EndEemprops())
    {
        return eer->getChi0();
    }
    else
    {
        fprintf(stderr, "No chi0 data for atype '%s'", atype.c_str());
        return -1;
    }
}

CommunicationStatus Poldata::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &filename_);
        gmx_send_str(cr, dest, &alexandriaPolarUnit_);
        gmx_send_str(cr, dest, &alexandriaPolarRef_);
        gmx_send_str(cr, dest, &alexandriaForcefield_);
        gmx_send_int(cr, dest, nexcl_);
        gmx_send_double(cr, dest, fudgeQQ_);
        gmx_send_double(cr, dest, fudgeLJ_);
        gmx_send_str(cr, dest, &gtVdwFunction_);
        gmx_send_str(cr, dest, &gtCombinationRule_);
        gmx_send_int(cr, dest, gtVdwFtype_);
        gmx_send_int(cr, dest, gtCombRule_);
        gmx_send_str(cr, dest, &millerTauUnit_);
        gmx_send_str(cr, dest, &millerAhpUnit_);
        gmx_send_str(cr, dest, &millerRef_);
        gmx_send_str(cr, dest, &bosquePolarUnit_);
        gmx_send_str(cr, dest, &bosqueRef_);
        gmx_send_str(cr, dest, &vsite_angle_unit_);
        gmx_send_str(cr, dest, &vsite_length_unit_);
        gmx_send_int(cr, dest, ptype_.size());
        gmx_send_int(cr, dest, vsite_.size());
        gmx_send_int(cr, dest, alexandria_.size());
        gmx_send_int(cr, dest, btype_.size());
        gmx_send_int(cr, dest, forces_.size());
        gmx_send_int(cr, dest, miller_.size());
        gmx_send_int(cr, dest, bosque_.size());
        gmx_send_int(cr, dest, symcharges_.size());
        gmx_send_int(cr, dest, static_cast<int>(eqdModel_));
        gmx_send_str(cr, dest, &eepReference_);
        gmx_send_int(cr, dest, eep_.size());

        /*Send ptype*/
        for (auto &ptype : ptype_)
        {
            cs = ptype.Send(cr, dest);
        }

        /*Send Ffatype*/
        for (auto &alexandria : alexandria_)
        {
            cs = alexandria.Send(cr, dest);
        }

        /*Send Vsite*/
        for (auto &vsite : vsite_)
        {
            cs = vsite.Send(cr, dest);
        }

        /*Send btype*/
        for (auto &btype : btype_)
        {
            gmx_send_str(cr, dest, &btype);
        }

        /*Send Listed Forces*/
        for (auto &force : forces_)
        {
            cs = force.Send(cr, dest);
        }

        /*Send Miller*/
        for (auto &miller : miller_)
        {
            cs = miller.Send(cr, dest);
        }

        /*Send Bosque*/
        for (auto &bosque : bosque_)
        {
            cs = bosque.Send(cr, dest);
        }

        /*Send Symcharge*/
        for (auto &symcharges : symcharges_)
        {
            cs = symcharges.Send(cr, dest);
        }

        /*Send Eemprops*/
        for (auto &eep : eep_)
        {
            cs = eep.Send(cr, dest);
        }
    }
    return cs;
}

CommunicationStatus Poldata::Receive(const t_commrec *cr, int src)
{
    size_t              nptype, nalexandria, nbtype, nforces, nvsite;
    size_t              nmiller, nbosque, nsymcharges, neep;
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &filename_);
        gmx_recv_str(cr, src, &alexandriaPolarUnit_);
        gmx_recv_str(cr, src, &alexandriaPolarRef_);
        gmx_recv_str(cr, src, &alexandriaForcefield_);
        nexcl_                = gmx_recv_int(cr, src);
        fudgeQQ_              = gmx_recv_double(cr, src);
        fudgeLJ_              = gmx_recv_double(cr, src);
        gmx_recv_str(cr, src, &gtVdwFunction_);
        gmx_recv_str(cr, src, &gtCombinationRule_);
        gtVdwFtype_           = gmx_recv_int(cr, src);
        gtCombRule_           = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &millerTauUnit_);
        gmx_recv_str(cr, src, &millerAhpUnit_);
        gmx_recv_str(cr, src, &millerRef_);
        gmx_recv_str(cr, src, &bosquePolarUnit_);
        gmx_recv_str(cr, src, &bosqueRef_);
        gmx_recv_str(cr, src, &vsite_angle_unit_);
        gmx_recv_str(cr, src, &vsite_length_unit_);
        nptype                = gmx_recv_int(cr, src);
        nvsite                = gmx_recv_int(cr, src);
        nalexandria           = gmx_recv_int(cr, src);
        nbtype                = gmx_recv_int(cr, src);
        nforces               = gmx_recv_int(cr, src);
        nmiller               = gmx_recv_int(cr, src);
        nbosque               = gmx_recv_int(cr, src);
        nsymcharges           = gmx_recv_int(cr, src);
        eqdModel_             = static_cast<ChargeDistributionModel>(gmx_recv_int(cr, src));
        gmx_recv_str(cr, src, &eepReference_);
        neep                  = gmx_recv_int(cr, src);


        /*Receive ptype*/
        ptype_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nptype); n++)
        {
            Ptype ptype;
            cs = ptype.Receive(cr, src);
            if (CS_OK == cs)
            {
                ptype_.push_back(ptype);
            }
        }

        /*Receive Ffatype*/
        alexandria_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nalexandria); n++)
        {
            Ffatype alexandria;
            cs = alexandria.Receive(cr, src);
            if (CS_OK == cs)
            {
                alexandria_.push_back(alexandria);
            }
        }

        /*Receive Ffatype*/
        vsite_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nvsite); n++)
        {
            Vsite vsite;
            cs = vsite.Receive(cr, src);
            if (CS_OK == cs)
            {
                vsite_.push_back(vsite);
            }
        }

        /*Receive btype*/
        btype_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nbtype); n++)
        {
            std::string btype;
            gmx_recv_str(cr, src, &btype);
            if (!btype.empty())
            {
                btype_.push_back(btype);
            }
        }

        /* Receive Listed Forces */
        forces_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nforces); n++)
        {
            ListedForces fs;
            cs = fs.Receive(cr, src);
            if (CS_OK == cs)
            {
                forces_.push_back(fs);
            }
        }

        /*Receive Miller*/
        miller_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nmiller); n++)
        {
            Miller miller;
            cs = miller.Receive(cr, src);
            if (CS_OK == cs)
            {
                miller_.push_back(miller);
            }
        }

        /*Receive Bosque*/
        bosque_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nbosque); n++)
        {
            Bosque bosque;
            cs = bosque.Receive(cr, src);
            if (CS_OK == cs)
            {
                bosque_.push_back(bosque);
            }
        }

        /*Receive Symcharges*/
        symcharges_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nsymcharges); n++)
        {
            Symcharges symcharges;
            cs = symcharges.Receive(cr, src);
            if (CS_OK == cs)
            {
                symcharges_.push_back(symcharges);
            }
        }

        /*Receive Eemprops*/
        eep_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < neep); n++)
        {
            Eemprops eep;
            cs = eep.Receive(cr, src);
            if (CS_OK == cs)
            {
                eep_.push_back(eep);
            }
        }
    }
    return cs;
}

void Poldata::broadcast(const t_commrec *cr)
{
    const int src = 0;
    if (MASTER(cr))
    {
        for (int dest = 1; dest < cr->nnodes; dest++)
        {
            auto cs = gmx_send_data(cr, dest);
            if (CS_OK == cs)
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Going to update Poldata on node %d\n", dest);
                }
                Send(cr, dest);
            }
            gmx_send_done(cr, dest);
        }
    }
    else
    {
        auto cs = gmx_recv_data(cr, src);
        if (CS_OK == cs)
        {
            auto cs = Receive(cr, src);
            if (CS_OK == cs)
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Poldata is updated on node %d\n", cr->nodeid);
                }
            }
            else
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Could not update Poldata on node %d\n", cr->nodeid);
                }
            }
        }
        gmx_recv_data(cr, src);
    }
}

EempropsConstIterator Poldata::ztype2Eem(const std::string &ztype) const
{
    return std::find_if(eep_.begin(), eep_.end(),
                        [ztype](Eemprops const &eep)
                        {
                            return (strcasecmp(eep.getName(), ztype.c_str()));
                        });
}

EempropsIterator Poldata::ztype2Eem(const std::string &ztype)
{
    return std::find_if(eep_.begin(), eep_.end(),
                        [ztype](Eemprops const &eep)
                        {
                            return (strcasecmp(eep.getName(), ztype.c_str()) == 0);
                        });
}
                                   
EempropsConstIterator Poldata::findEem(const std::string &atype) const
{
    std::string nn;
    auto        fa = findAtype(atype);
    if (fa != getAtypeEnd())
    {
        auto eqdModel = getEqdModel();
        if (eqdModel == eqdRappe || eqdModel == eqdBultinck ||
            eqdModel == eqdYang)
        {
            nn = fa->getElem();
        }
        else
        {
            nn = fa->getZtype();
        }
    }
    else
    {
        nn = atype;
    }
    return std::find_if(eep_.begin(), eep_.end(),
                        [nn](Eemprops const &eep)
                        {
                            return (strcasecmp(eep.getName(), nn.c_str()) == 0);
                        });
}

EempropsIterator Poldata::findEem(const std::string &atype)
{
    std::string nn;
    auto        fa = findAtype(atype);
    if (fa != getAtypeEnd())
    {
        auto eqdModel = getEqdModel();
        if (eqdModel == eqdRappe || eqdModel == eqdBultinck ||
            eqdModel == eqdYang)
        {
            nn = fa->getElem();
        }
        else
        {
            nn = fa->getZtype();
        }
    }
    else
    {
        nn = atype;
    }
    return std::find_if(eep_.begin(), eep_.end(),
                        [nn](Eemprops const &eep)
                        {
                            return (strcasecmp(eep.getName(), nn.c_str()) == 0);
                        });
}

} // namespace alexandria
