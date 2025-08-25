/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020,2025
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

#include "act/molprop/molprop_util.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "act/alexandria/actmol.h"
#include "act/molprop/molprop_xml.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/textwriter.h"

namespace alexandria
{

void generate_index(std::vector<MolProp> *mp)
{
    int index = 1;
    for (auto mpi = mp->begin(); mpi < mp->end(); ++mpi)
    {
        // We have to set some kind of index here to distinguish
        // the compounds.
        // TODO: implement alexandria ID here.
        mpi->setIndex(index++);
    }
}

/*! \brief
 * Merge multiple molprops for molecules into one, e.g. one molprop for water
 * one for methane etc.
 *
 * \param[inout]  mp        The vector of MolProps
 * \param[in]     doubles   File name for dumping output, can be nullptr
 * \param[in]     bForceMerge If true all molprops for a compound are merged
 * \return the number of remaining molprops
 * \ingroup module_alexandria
 */
static void MergeDoubleMolprops(MsgHandler                       *msg_handler,
                                std::vector<alexandria::MolProp> *mp,
                                std::vector<std::string>         *warnings)
{
    alexandria::MolPropIterator mpi, mmm[2];
    std::string                 molname[2];
    std::string                 form[2];
    int   i, ndouble = 0;
    bool  bDouble;
    int   cur   = 0;
#define prev (1-cur)
    gmx::TextWriter *tw = nullptr;
    if (msg_handler && msg_handler->debug())
    {
        tw = msg_handler->twDebug();
    }
    i = 0;
    for (mpi = mp->begin(); (mpi < mp->end()); )
    {
        bDouble = false;
        mpi->Dump(tw);
        mmm[cur]     = mpi;
        molname[cur] = mpi->getMolname();
        form[cur]    = mpi->formula();
        if (i > 0)
        {
            if (molname[prev] == molname[cur])
            {
                auto warn = mmm[prev]->sameCompound(&(*mmm[cur]));
                if (!warn.empty())
                {
                    for(const auto &w : warn)
                    {
                        warnings->push_back(w);
                    }
                }
                else
                {
                    auto warn = mmm[prev]->Merge(&(*mmm[cur]));
                    if (!warn.empty())
                    {
                        for(const auto &w : warn)
                        {
                            warnings->push_back(w);
                        }
                    }
                    mpi    = mp->erase(mmm[cur]);

                    bDouble = true;
                    ndouble++;
                    if (mpi != mp->end())
                    {
                        mmm[cur]      = mpi;
                        molname[cur]  = mmm[cur]->getMolname();
                        form[cur]     = mmm[cur]->formula();
                    }
                }
            }
        }
        if (!bDouble)
        {
            cur = prev;
            mpi++;
            i++;
        }
    }
    for (mpi = mp->begin(); (mpi < mp->end()); mpi++)
    {
        mpi->Dump(tw);
    }
    if (msg_handler)
    {
        msg_handler->write(gmx::formatString("There were %d double entries, leaving %d after merging.\n",
                                             ndouble, (int)mp->size()));
    }
}

std::vector<std::string> merge_xml(MsgHandler                             *msg_handler,
                                   const gmx::ArrayRef<const std::string> &filens,
                                   std::vector<alexandria::MolProp>       *mpout)
{
    int tmp;

    for (auto &fn : filens)
    {
        std::vector<alexandria::MolProp> mp;
        if (!gmx_fexist(fn.c_str()))
        {
            continue;
        }
        MolPropRead(msg_handler, fn.c_str(), &mp);
        for (auto mpi : mp)
        {
            mpout->push_back(std::move(mpi));
        }
    }
    tmp = mpout->size();
    gmx::TextWriter *tw = nullptr;
    if (msg_handler && msg_handler->debug())
    {
        tw = msg_handler->twDebug();

        tw->writeStringFormatted("mpout->size() = %u mpout->max_size() = %u\n",
                                (unsigned int)mpout->size(), (unsigned int)mpout->max_size());
        for (auto mpi = mpout->begin(); mpi < mpout->end(); ++mpi)
        {
            mpi->Dump(tw);
        }
    }
    MolSelect gms;
    MolPropSort(msg_handler, mpout, MPSA_MOLNAME, nullptr, gms);
    std::vector<std::string> warnings;
    MergeDoubleMolprops(msg_handler, mpout, &warnings);
    printf("There were %d total molecules before merging, %d after.\n",
           tmp, (int)mpout->size());

    return warnings;
}

static bool comp_mp_molname(alexandria::MolProp ma,
                            alexandria::MolProp mb)
{
    std::string mma = ma.getMolname();
    std::string mmb = mb.getMolname();

    return (mma.compare(mmb) < 0);
}

static bool comp_mp_formula(alexandria::MolProp ma,
                            alexandria::MolProp mb)
{
    std::string fma = ma.formula();
    std::string fmb = mb.formula();

    if (fma.compare(fmb) < 0)
    {
        return true;
    }
    else
    {
        return comp_mp_molname(ma, mb);
    }
}

gmx_atomprop_t my_aps;

static int count_elements(std::map<const char *, int> comp, const char *elem)
{
    auto mm = comp.find(elem);
    if (mm == comp.end())
    {
        return 0;
    }
    else
    {
        return mm->second;
    }
}

static bool comp_mp_elem(alexandria::MolProp ma,
                         alexandria::MolProp mb)
{
    int         i;
    auto mcia = ma.composition();
    auto mcib = mb.composition();

    if (!mcia.empty() && !mcib.empty())
    {
        int d = count_elements(mcia, "C") - count_elements(mcib, "C");
        if (d < 0)
        {
            return true;
        }
        else if (d > 0)
        {
            return false;
        }
        else
        {
            for (i = 1; (i <= 109); i++)
            {
                if (i != 6)
                {
                    const char *ee = gmx_atomprop_element(my_aps, i);
                    if (nullptr != ee)
                    {
                        d = count_elements(mcia, ee) - count_elements(mcib, ee);
                        if (d == 0)
                        {
                            continue;
                        }
                        else
                        {
                            return (d < 0);
                        }
                    }
                }
            }
        }
    }
    return comp_mp_molname(ma, mb);
}

static bool comp_mp_index(alexandria::MolProp ma,
                          alexandria::MolProp mb)
{
    return (ma.getIndex() < mb.getIndex());
}

void MolPropSort(MsgHandler                                *msg_handler,
                 std::vector<alexandria::MolProp>          *mp,
                 MolPropSortAlgorithm mpsa, gmx_atomprop_t  apt,
                 const MolSelect                           &gms)
{
    gmx::TextWriter *tw = nullptr;
    if (msg_handler && msg_handler->debug())
    {
        tw = msg_handler->twDebug();
    }
    if (msg_handler)
    {
        msg_handler->write(gmx::formatString("There are %d molprops. Will now sort them.\n", (int)mp->size()));
    }
    if (tw)
    {
        for (auto mpi = mp->begin(); (mpi < mp->end()); mpi++)
        {
            mpi->Dump(tw);
        }
    }
    switch (mpsa)
    {
        case MPSA_MOLNAME:
            std::sort(mp->begin(), mp->end(), comp_mp_molname);
            break;
        case MPSA_FORMULA:
            std::sort(mp->begin(), mp->end(), comp_mp_formula);
            break;
        case MPSA_COMPOSITION:
            if (nullptr != apt)
            {
                my_aps = apt;
                std::sort(mp->begin(), mp->end(), comp_mp_elem);
                my_aps = nullptr;
            }
            else
            {
                gmx_fatal(FARGS, "Requesting a composition sort but atomprop is nullptr");
            }
            break;
        case MPSA_SELECTION:
            if (gms.nMol() > 0)
            {
                for (auto mpi = mp->begin(); mpi < mp->end(); mpi++)
                {
                    int index;
                    if (!gms.index(mpi->getIupac(), &index))
                    {
                        gmx_fatal(FARGS, "Cannot find index for %s", mpi->getIupac().c_str());
                    }
                    mpi->setIndex(index);
                }
                std::sort(mp->begin(), mp->end(), comp_mp_index);
            }
            else
            {
                gmx_fatal(FARGS, "Need molecule selection to sort on");
            }
            break;
        default:
            gmx_incons("Invalid algorithm for MolPropSort");
    }
}

void QmCount::addConf(const std::string &conformation)
{
    auto s = std::find_if(conf_.begin(), conf_.end(),
                          [conformation](const std::string &s)
                          { return (s.compare(conformation) == 0); });
    if (s == conf_.end())
    {
        conf_.push_back(conformation);
    }
}

int QmCount::qmCalcCount(const std::string &method,
                         const std::string &basis,
                         const std::string &type) const
{
    auto a = findCalc(method, basis, type);
    if (qmc_.end() == a)
    {
        return 0;
    }
    else
    {
        return a->count();
    }
}

void QmCount::addCalc(const std::string &method,
                      const std::string &basis,
                      const std::string &type)
{
    auto a = findCalc(method, basis, type);

    if (qmc_.end() == a)
    {
        QmCalc qmc(method, basis, type);
        qmc_.push_back(qmc);
    }
    else
    {
        a->increment();
    }
}

void find_calculations(const std::vector<alexandria::MolProp> &mp,
                       MolPropObservable                       mpo,
                       const char                             *fc_types,
                       QmCount                                *qmc)
{
    std::vector<std::string> types;

    for (auto &mpi : mp)
    {
        for (auto &ei : mpi.experimentConst())
        {
            if (ei.hasProperty(mpo))
            {
                qmc->addConf(ei.getConformation());
            }
        }
    }
    if (nullptr != fc_types)
    {
        std::vector<std::string> qm = split(fc_types, ':');
        for (auto pqm = qm.begin(); (pqm < qm.end()); ++pqm)
        {
            if (pqm->length() > 0)
            {
                std::vector<std::string> ll = split(pqm->c_str(), '/');
                if (ll.size() == 3)
                {
                    std::vector<std::string>::iterator ti;
                    qmc->addCalc(ll[0], ll[1], ll[2]);
                    for (ti = types.begin(); (ti < types.end()); ti++)
                    {
                        if (0 == strcasecmp(ti->c_str(), ll[2].c_str()))
                        {
                            break;
                        }
                    }
                    if (ti == types.end())
                    {
                        types.push_back(ll[2].c_str());
                    }
                }
            }
        }
    }

    for (auto &mpi : mp)
    {
        for (auto &ci : mpi.experimentConst())
        {
            if (dsExperiment ==  ci.dataSource())
            {
                continue;
            }
#ifdef OLD
            for (auto ti = types.begin(); (ti < types.end()); ti++)
            {
                if ((nullptr == fc_str) || (qmc->qmCalcCount(ci.getMethod(),
                                                             ci.getBasisset(),
                                                             *ti) > 0))
                {
                    double              T, value, error;
                    std::vector<double> vec;
                    tensor              quadrupole;
                    if (ci.getVal(ti->c_str(), mpo, &value, &error, &T,
                                  &vec, quadrupole))
                    {
                        qmc->addCalc(ci.getMethod(),
                                     ci.getBasisset(),
                                     *ti);
                    }
                }
            }
#endif
        }
    }
    for (auto q = qmc->beginCalc(); q < qmc->endCalc(); ++q)
    {
        /* Since we initialized these we have counted one extra */
        if (nullptr != fc_types)
        {
            q->decrement();
        }
    }
}

} // namespace alexandria

void splitLot(const char  *lot,
              std::string *method,
              std::string *basis)
{
    method->clear();
    basis->clear();
    if (nullptr != lot && strlen(lot) > 0)
    {
        std::vector<std::string> ll = split(lot, '/');
        if (ll.size() == 2)
        {
            method->assign(ll[0]);
            basis->assign(ll[1]);
        }
    }
}

