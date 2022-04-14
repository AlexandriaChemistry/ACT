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

#include <math.h> 
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <vector>

#ifdef HAVE_LIBSQLITE3
#include <sqlite3.h>
#endif
#include "act/molprop/molpropobservable.h"
#include "act/molprop/molprop_sqlite3.h"

namespace alexandria
{

/*! \brief class to handle synonyms of compound names
 */ 
class Synonym
{
    private:
        //! Molecule name
        std::string molname_;
        //! Internation Union of Pure and Applied Chemistry name
        std::string iupac_;
    public:
        /*! \brief Constructor
         * \param[in] molname The molecule name
         * \param[in] iupac   The iupac name
         */
        Synonym(const std::string molname,
                const std::string iupac) : molname_(molname), iupac_(iupac) {}

        //! \return the iupac name
        const std::string &iupac() const { return iupac_; }
        //! \return the molecule name
        const std::string &molname() const { return molname_; }
};

/*! \brief Class to detail with the compound classes a particular
 * compound is part of. For instance, alkanes, amines.
 */
class Classes
{
    private:
        //! IUPAC name
        std::string              iupac_;
        //! Vector of compound classes.
        std::vector<std::string> classes_;
    public:
        /*! \brief Constructor
         * \param[in] iupac The IUPAC name
         */
        Classes(const std::string iupac) : iupac_(iupac) {}

        /*! \brief Add a class to the list
         * \param[in] klas The new class
         */
        void addClass(const std::string &klas) { classes_.push_back(klas); }

        //! \return the IUPAC name
        const std::string &iupac() const { return iupac_; }

        //! \return the list of classes
        const std::vector<std::string> &classes() { return classes_; }
};

#ifdef HAVE_LIBSQLITE3
static void check_sqlite3(sqlite3 *db, const char *extra, int rc)
{
    const char *msg;

    if (nullptr != db)
    {
        if (SQLITE_OK != sqlite3_errcode(db))
        {
            msg = sqlite3_errmsg(db);
            sqlite3_close(db);
            sqlite3_shutdown();
            gmx_fatal(FARGS, "%s: %s", extra, msg);
        }
    }
    else if (SQLITE_OK != rc)
    {
        gmx_fatal(FARGS, "%s", extra);
    }
}
#endif

static void getSynonyms(sqlite3              *db,
                        std::vector<Synonym> &syn,
                        int                   nMol)
{
    sqlite3_stmt *stmt2 = nullptr;
    char          sql_str[1024];
    int           rc;

    /* Make renaming table */
    snprintf(sql_str, sizeof(sql_str),
             "SELECT syn.name,mol.iupac FROM molecules as mol,synonyms as syn WHERE syn.molid=mol.molid ORDER by syn.name");

    if (nullptr != debug)
    {
        fprintf(debug, "sql_str = '%s'\n", sql_str);
    }

    check_sqlite3(db, "Preparing statement",
                  sqlite3_prepare_v2(db, sql_str, 1+strlen(sql_str), &stmt2, nullptr));
    do
    {
        rc = sqlite3_step(stmt2);
        if (SQLITE_ROW == rc)
        {
            syn.push_back(Synonym((char *)sqlite3_column_text(stmt2, 0),
                                  (char *)sqlite3_column_text(stmt2, 1)));
        }
        else if (SQLITE_DONE != rc)
        {
            check_sqlite3(db, "Stepping", rc);
        }
        else
        {
            printf("There are %d synonyms for %d molecules.\n",
                   static_cast<int>(syn.size()), nMol);
        }
    }
    while (SQLITE_ROW == rc);
    check_sqlite3(db, "Resetting sqlite3 statement",
                  sqlite3_reset(stmt2));
    check_sqlite3(db, "Finalizing sqlite3 statement",
                  sqlite3_finalize(stmt2));
}

static void getClasses(sqlite3              *db,
                       std::vector<Classes> &classes,
                       int                   nMol)
{
    sqlite3_stmt *stmt2 = nullptr;
    char          sql_str[1024];
    int           rc;

    /* Make renaming table */
    snprintf(sql_str, sizeof(sql_str),
             "SELECT mol.iupac,class.class FROM molecules as mol,classification as class,link_mol_class as lmc WHERE (lmc.molid=mol.molid) and (lmc.classid=class.classid) ORDER by mol.iupac");

    if (nullptr != debug)
    {
        fprintf(debug, "sql_str = '%s'\n", sql_str);
    }

    check_sqlite3(db, "Preparing statement",
                  sqlite3_prepare_v2(db, sql_str, 1+strlen(sql_str), &stmt2, nullptr));
    do
    {
        rc = sqlite3_step(stmt2);
        if (SQLITE_ROW == rc)
        {
            const char *iupac = (char *)sqlite3_column_text(stmt2, 0);
            const char *klass = (char *)sqlite3_column_text(stmt2, 1);
            auto        s     = std::find_if(classes.begin(), classes.end(),
                                             [iupac](Classes const &c)
                                             { return c.iupac().compare(iupac) == 0; });
            if (s == classes.end())
            {
                std::string i(iupac);
                classes.push_back(i);
                classes.back().addClass(klass);
            }
            else
            {
                s->addClass(klass);
            }
        }
        else if (SQLITE_DONE != rc)
        {
            check_sqlite3(db, "Stepping", rc);
        }
        else
        {
            printf("There are %d classes for %d molecules.\n",
                   static_cast<int>(classes.size()), nMol);
        }
    }
    while (SQLITE_ROW == rc);
    check_sqlite3(db, "Resetting sqlite3 statement",
                  sqlite3_reset(stmt2));
    check_sqlite3(db, "Finalizing sqlite3 statement",
                  sqlite3_finalize(stmt2));
}

void ReadSqlite3(const char           *sqlite_file,
                 std::vector<MolProp> *mp,
                 double                ref_temperature)
{
#ifdef HAVE_LIBSQLITE3
    std::string                 cas2, csid2;

    sqlite3                    *db   = nullptr;
    sqlite3_stmt               *stmt = nullptr;
    const char                 *cas, *csid, *prop, *unit, *ref;
    double                      value, error, temperature;
    int                         cidx, rc, nbind, nexp_prop, theory, preferred;
    std::vector<Synonym>        synonyms;
    std::vector<Classes>        classes;
    
    if (nullptr == sqlite_file)
    {
        return;
    }   
    check_sqlite3(nullptr, "Initializing sqlite", sqlite3_initialize());
    check_sqlite3(nullptr, "Opening sqlite database in read-only mode",
                  sqlite3_open_v2(sqlite_file, &db, SQLITE_OPEN_READONLY, nullptr));

    /* Now database is open and everything is Hunky Dory */
    printf("Opened SQLite3 database %s\n", sqlite_file);

    // First get the synonyms out.
    getSynonyms(db, synonyms, mp->size());

    // Now get the classes out.
    getClasses(db, classes, mp->size());

    /* Now present a query statement */
    nexp_prop = 0;
    auto sql_str = gmx::formatString("SELECT distinct mol.iupac,mol.cas,mol.csid,pt.prop,pt.unit_text,ref.ref,gp.temperature,gp.value,gp.error, gp.preferred,ds.theory FROM molecules as mol,molproperty as gp,proptypes as pt, datasource as ds,phasetype as ph,reference as ref,link_molprop_ref as lmr WHERE ((gp.phaseid=ph.phaseid) AND (ph.phase='gas') AND (mol.molid = gp.molid) AND (gp.propid = pt.propid) AND (gp.molpropid = lmr.molpropid) AND (lmr.refid = ref.refid) AND (gp.srcid = ds.srcid) AND (upper(?) = upper(mol.iupac)));");
    check_sqlite3(db, "Preparing sqlite3 statement",
                  sqlite3_prepare_v2(db, sql_str.c_str(), 1+sql_str.size(), &stmt, nullptr));

    if (nullptr != debug)
    {
        fprintf(debug, "sql_str = '%s'\nvariable = '%s'\n", sql_str.c_str(),
                sqlite3_bind_parameter_name(stmt, 1));
        nbind = sqlite3_bind_parameter_count(stmt);
        fprintf(debug, "%d binding parameter(s) in the statement\n%s\n", nbind, sql_str.c_str());
    }
    for (auto mpi = mp->begin(); (mpi < mp->end()); mpi++)
    {
        const std::string molname = mpi->getMolname();
        auto              keyptr  = std::find_if(synonyms.begin(), synonyms.end(),
                                                 [molname](Synonym const &s)
                                                 { return molname.compare(s.molname()) == 0; });
        if (synonyms.end() == keyptr)
        {
            fprintf(stderr, "Warning: missing iupac for %s. Will be ignored.\n",
                    molname.c_str());
        }
        else
        {
            if (nullptr != debug)
            {
                fprintf(debug, "Going to query for '%s'\n", keyptr->iupac().c_str());
            }
            check_sqlite3(db, "Binding text",
                          sqlite3_bind_text(stmt, 1, keyptr->iupac().c_str(), -1, SQLITE_STATIC));
            do
            {
                rc = sqlite3_step(stmt);
                if (SQLITE_ROW == rc)
                {
                    cidx   = 0;
                    const char *iupac2 = (char *)sqlite3_column_text(stmt, cidx++);
                    if (strcasecmp(keyptr->iupac().c_str(), iupac2) != 0)
                    {
                        gmx_fatal(FARGS, "Selected '%s' from database but got '%s'. WTF?!",
                                  keyptr->iupac().c_str(), iupac2);
                    }
                    cas            = (char *)sqlite3_column_text(stmt, cidx++);
                    csid           = (char *)sqlite3_column_text(stmt, cidx++);
                    prop           = (char *)sqlite3_column_text(stmt, cidx++);
                    unit           = (char *)sqlite3_column_text(stmt, cidx++);
                    ref            = (char *)sqlite3_column_text(stmt, cidx++);
                    temperature    = sqlite3_column_double(stmt, cidx++);
                    value          = sqlite3_column_double(stmt, cidx++);
                    error          = sqlite3_column_double(stmt, cidx++);
                    preferred      = sqlite3_column_int(stmt, cidx++);
                    theory         = sqlite3_column_int(stmt, cidx++);
                    if (debug)
                    {
                        fprintf(debug, "Found: mol %s prop %s value %g temp %g theory %d preferred %d\n",
                                molname.c_str(), prop, value, temperature, theory, preferred);
                    }
                    if (fabs(ref_temperature-temperature) < 0.1)
                    {
                        iqmType iqm = iqmType::Exp;
                        if (1 == theory)
                        {
                            iqm = iqmType::QM;
                        }
                        MolPropObservable mpo;
                        if (!stringToMolPropObservable(prop, &mpo))
                        {
                            if (debug)
                            {
                                fprintf(debug, "Unknown property %s\n", prop);
                            }
                        }
                        else if ((iqm == iqmType::Exp && 1==preferred) ||
                                 (iqm == iqmType::QM))
                        {
                            if (iqm == iqmType::Exp)
                            {
                                nexp_prop++;
                            }
                            alexandria::Experiment exper(ref, "minimum");
                            std::string exp_type("experiment");
                            GenericProperty *gp;
                            switch (mpo)
                            {
                            case MolPropObservable::POLARIZABILITY:
                                {
                                    gp = new MolecularPolarizability(exp_type, unit,
                                                                     temperature, 0, 0, 0, 0, 0, 0,
                                                                     value, error);
                                    break;
                                }
                            case MolPropObservable::DIPOLE:
                                {
                                    auto mm = new MolecularMultipole(exp_type, unit, temperature,
                                                                MolPropObservable::DIPOLE);
                                    mm->setValue("average", value);
                                    mm->setValue("error", error);
                                    gp = reinterpret_cast<GenericProperty *>(mm);
                                    break;
                                }
                            case MolPropObservable::DGFORM:
                            case MolPropObservable::DHFORM:
                            case MolPropObservable::DSFORM:
                            case MolPropObservable::ENTROPY:
                            case MolPropObservable::STRANS:
                            case MolPropObservable::SROT:
                            case MolPropObservable::SVIB:
                            case MolPropObservable::ZPE:
                            case MolPropObservable::CP:
                            case MolPropObservable::CV:
                                {
                                    gp = new MolecularEnergy(mpo, exp_type, unit, temperature, ePhase::GAS, value, error);
                                    break;
                                }
                            default:
                                gmx_fatal(FARGS, "Unsupported property %s", prop);
                            }
                            exper.addProperty(mpo, gp);
                            mpi->AddExperiment(exper);
                        }
                    }
                    const char *iupac = keyptr->iupac().c_str();
                    auto        cptr  = std::find_if(classes.begin(), classes.end(),
                                                     [iupac](Classes const &c)
                                                     { return c.iupac().compare(iupac) == 0; });
                    if (cptr != classes.end())
                    {
                        for (auto &c : cptr->classes())
                        {
                            mpi->AddCategory(c);
                        }
                    }
                    if (strlen(cas) > 0)
                    {
                        cas2 = mpi->getCas();
                        if ((cas2.length() > 0) &&
                            (strcmp(cas, cas2.c_str()) != 0))
                        {
                            fprintf(stderr, "cas in molprop %s not the same as database %s for %s\n", cas2.c_str(), cas, iupac);
                        }
                        mpi->SetCas(cas);
                    }
                    if (strlen(csid) > 0)
                    {
                        csid2 = mpi->getCid();
                        if ((csid2.length() > 0) &&
                            (strcmp(csid, csid2.c_str()) != 0))
                        {
                            fprintf(stderr, "csid in molprop %s not the same as database %s for %s\n", csid2.c_str(), csid, iupac);
                        }
                        mpi->SetCid(csid);
                    }
                }
                else if (SQLITE_DONE != rc)
                {
                    check_sqlite3(db, "Stepping", rc);
                }
                else if (nullptr != debug)
                {
                    fprintf(debug, "Done finding rows for %s\n", keyptr->iupac().c_str());
                }
            }
            while (SQLITE_ROW == rc);
            sqlite3_clear_bindings(stmt);
            check_sqlite3(db, "Resetting sqlite3 statement", sqlite3_reset(stmt));
        }
    }
    check_sqlite3(db, "Finalizing sqlite3 statement", sqlite3_finalize(stmt));
    check_sqlite3(nullptr, "Closing sqlite database", sqlite3_close(db));
    check_sqlite3(nullptr, "Shutting down sqlite. Sqlite3 code %d.", sqlite3_shutdown());
    printf("Extracted %d data points at %0.2f (K) from sql database\n", nexp_prop, ref_temperature);
    
#else
    fprintf(stderr, "No support for sqlite3 database in this executable.\n");
    fprintf(stderr, "Please rebuild gromacs with cmake flag -DGMX_SQLITE3=ON set.\n");
#endif
}

} // namespace alexandria
