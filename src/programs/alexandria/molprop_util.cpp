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

#include "molprop_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"

#include "composition.h"
#include "molprop_xml.h"

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
        mpi->SetIndex(index++);
    }
}

void generate_composition(std::vector<MolProp> &mp,
                          const Poldata        *pd)
{
    int              nOK = 0;
    CompositionSpecs cs;

    for (auto &mpi : mp)
    {
        for (auto csi = cs.beginCS(); (csi < cs.endCS()); ++csi)
        {
            mpi.DeleteComposition(csi->name());
        }
        if (true == mpi.GenerateComposition(pd))
        {
            nOK++;
        }
        else if (debug)
        {
            fprintf(debug, "Failed to make composition for %s\n",
                    mpi.getMolname().c_str());
        }
    }
    if (mp.size() > 1)
    {
        printf("Generated composition for %d out of %d molecules.\n",
               nOK, (int)mp.size());
    }
}

void generate_formula(std::vector<MolProp> &mp,
                      gmx_atomprop_t        ap)
{
    for (auto &mpi : mp)
    {
        mpi.GenerateFormula(ap);
    }
}

int MergeDoubleMolprops(std::vector<alexandria::MolProp> *mp, 
                        char                             *doubles,
                        bool                              bForceMerge)
{
    alexandria::MolPropIterator mpi, mmm[2];
    std::string                 molname[2];
    std::string                 form[2];

    FILE *fp;
    int   i, ndouble = 0;
    bool  bForm, bName, bDouble;
    int   nwarn = 0;
    int   cur   = 0;
#define prev (1-cur)

    if (nullptr != doubles)
    {
        fp = fopen(doubles, "w");
    }
    else
    {
        fp = nullptr;
    }
    i = 0;
    for (mpi = mp->begin(); (mpi < mp->end()); )
    {
        bDouble = false;
        mpi->Dump(debug);
        mmm[cur]     = mpi;
        molname[cur] = mpi->getMolname();
        form[cur]    = mpi->formula();
        if (i > 0)
        {
            bForm = (form[prev] == form[cur]);
            bName = (molname[prev] == molname[cur]);
            if (bName)
            {
                if (!bForm && debug)
                {
                    fprintf(debug, "%s %s with formulae %s - %s\n",
                            bForceMerge ? "Merging molecules" : "Found molecule",
                            molname[prev].c_str(), form[prev].c_str(), form[cur].c_str());
                }
                if (bForceMerge || bForm)
                {
                    if (fp)
                    {
                        fprintf(fp, "%5d  %s\n", ndouble+1, molname[prev].c_str());
                    }
                    nwarn += mmm[prev]->Merge(mmm[cur]);
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
        mpi->Dump(debug);
    }
    if (fp)
    {
        fclose(fp);
    }
    printf("There were %d double entries, leaving %d after merging.\n",
           ndouble, (int)mp->size());
    return nwarn;
}

static void dump_mp(std::vector<alexandria::MolProp> *mp)
{
    alexandria::MolPropIterator mpi;
    FILE *fp;

    fp = fopen("dump_mp.dat", "w");

    for (mpi = mp->begin(); (mpi < mp->end()); mpi++)
    {
        fprintf(fp, "%-20s  %s\n", mpi->formula().c_str(),
                mpi->getMolname().c_str());
    }

    fclose(fp);
}

int merge_xml(gmx::ArrayRef<const std::string> filens,
              std::vector<alexandria::MolProp> *mpout,
              char *outf, char *sorted, char *doubles,
              gmx_atomprop_t ap,
              const Poldata &pd,
              bool bForceMerge)
{
    int npout = 0, tmp;

    for (auto &fn : filens)
    {
        std::vector<alexandria::MolProp> mp;
        if (!gmx_fexist(fn.c_str()))
        {
            continue;
        }
        MolPropRead(fn.c_str(), &mp);
        generate_composition(mp, &pd);
        generate_formula(mp, ap);
        for (auto mpi : mp)
        {
            mpout->push_back(std::move(mpi));
        }
    }
    tmp = mpout->size();
    if (nullptr != debug)
    {
        fprintf(debug, "mpout->size() = %u mpout->max_size() = %u\n",
                (unsigned int)mpout->size(), (unsigned int)mpout->max_size());
        for (auto mpi = mpout->begin(); mpi < mpout->end(); ++mpi)
        {
            mpi->Dump(debug);
        }
    }
    MolSelect gms;
    MolPropSort(mpout, MPSA_MOLNAME, nullptr, gms);
    int nwarn = MergeDoubleMolprops(mpout, doubles, bForceMerge);
    printf("There were %d total molecules before merging, %d after.\n",
           tmp, (int)mpout->size());
    if (outf)
    {
        printf("There are %d entries to store in output file %s\n", npout, outf);
        MolPropWrite(outf, mpout, false);
    }
    if (sorted)
    {
        MolPropSort(mpout, MPSA_FORMULA, nullptr, gms);
        MolPropWrite(sorted, mpout, false);
        dump_mp(mpout);
    }
    return nwarn;
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

static bool comp_mp_elem(alexandria::MolProp ma,
                         alexandria::MolProp mb)
{
    int         i;
    alexandria::MolecularCompositionIterator mcia, mcib;
    std::string bosque("bosque"), C("C");

    mcia = ma.SearchMolecularComposition(bosque);
    mcib = mb.SearchMolecularComposition(bosque);

    if ((mcia != ma.EndMolecularComposition()) &&
        (mcib != mb.EndMolecularComposition()))
    {
        int d = mcia->CountAtoms(C) - mcib->CountAtoms(C);

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
                        std::string elem(ee);
                        d = mcia->CountAtoms(elem) - mcib->CountAtoms(elem);
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

void MolPropSort(std::vector<alexandria::MolProp> *mp,
                 MolPropSortAlgorithm mpsa, gmx_atomprop_t apt,
                 const MolSelect &gms)
{
    printf("There are %d molprops. Will now sort them.\n", (int)mp->size());
    for (auto mpi = mp->begin(); (mpi < mp->end()); mpi++)
    {
        if (nullptr != debug)
        {
            mpi->Dump(debug);
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
                    int index = gms.index(mpi->getIupac());
                    mpi->SetIndex(index);
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

void find_calculations(std::vector<alexandria::MolProp> &mp,
                       MolPropObservable                 mpo,
                       const char                       *fc_str,
                       QmCount                          *qmc)
{
    std::vector<std::string> types;

    for (auto mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        for (auto ei = mpi->BeginExperiment(); (ei < mpi->EndExperiment()); ei++)
        {
            qmc->addConf(ei->getConformation());
        }
    }
    if (nullptr != fc_str)
    {
        std::vector<std::string> qm = split(fc_str, ':');
        int n = 0;
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
            n++;
        }
    }

    for (auto mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        for (auto ci = mpi->BeginExperiment(); (ci < mpi->EndExperiment()); ++ci)
        {
            if (dsExperiment ==  ci->dataSource())
            {
                continue;
            }
            for (auto ti = types.begin(); (ti < types.end()); ti++)
            {
                if ((nullptr == fc_str) || (qmc->qmCalcCount(ci->getMethod(),
                                                             ci->getBasisset(),
                                                             *ti) > 0))
                {
                    double T, value, error;
                    rvec   vec;
                    tensor quadrupole;
                    if (ci->getVal(ti->c_str(), mpo, &value, &error, &T,
                                   vec, quadrupole))
                    {
                        qmc->addCalc(ci->getMethod(),
                                     ci->getBasisset(),
                                     *ti);
                    }
                }
            }
        }
    }
    for (auto q = qmc->beginCalc(); q < qmc->endCalc(); ++q)
    {
        /* Since we initialized these we have counted one extra */
        if (nullptr != fc_str)
        {
            q->decrement();
        }
    }
}

} // namespace alexandria
