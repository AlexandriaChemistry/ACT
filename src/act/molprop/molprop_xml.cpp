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

#include "act/molprop/molprop_xml.h"

#include <stdlib.h>
#include <string.h>

#include <map>
#include <vector>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "act/molprop/molprop.h"
#include "act/molprop/multipole_names.h"
#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "act/utility/units.h"
#include "act/utility/xml_util.h"

namespace alexandria 
{

std::map<JobType, const char *> job_name =
{
    { JobType::OPT,    "Opt"    },
    { JobType::POP,    "Pop"    },
    { JobType::POLAR,  "POLAR"  },
    { JobType::G2,     "G2"     },
    { JobType::G3,     "G3"     },
    { JobType::G4,     "G4"     },
    { JobType::CBSQB3, "CBSQB3" },
    { JobType::W1U,    "W1U"    },
    { JobType::W1BD,   "W1BD"   },
    { JobType::SP,     "SP"     },
    { JobType::UNKNOWN,"unknown"}
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

static const char *xmltypes[] = {
    nullptr,
    "XML_ELEMENT_NODE",
    "XML_ATTRIBUTE_NODE",
    "XML_TEXT_NODE",
    "XML_CDATA_SECTION_NODE",
    "XML_ENTITY_REF_NODE",
    "XML_ENTITY_NODE",
    "XML_PI_NODE",
    "XML_COMMENT_NODE",
    "XML_DOCUMENT_NODE",
    "XML_DOCUMENT_TYPE_NODE",
    "XML_DOCUMENT_FRAG_NODE",
    "XML_NOTATION_NODE",
    "XML_HTML_DOCUMENT_NODE",
    "XML_DTD_NODE",
    "XML_ELEMENT_DECL",
    "XML_ATTRIBUTE_DECL",
    "XML_ENTITY_DECL",
    "XML_NAMESPACE_DECL",
    "XML_XINCLUDE_START",
    "XML_XINCLUDE_END"
};
#define NXMLTYPES sizeof(xmltypes)/sizeof(xmltypes[0])

enum class MolPropXml {
    MOLECULES,
    MOLECULE,
    FORMULA,
    MOLNAME,
    MASS,
    MOLINFO,
    IUPAC,
    CAS,
    CID,
    INCHI,
    MULTIPLICITY,
    CHARGE,
    CATEGORY,
    CATNAME,
    EXPERIMENT,
    POLARIZABILITY,
    ENERGY,
    DIPOLE,
    QUADRUPOLE,
    OCTUPOLE,
    HEXADECAPOLE,
    POTENTIAL,
    NAME,
    AVERAGE,
    ERROR,
    TEMPERATURE,
    PHASE,
    METHOD,
    REFERENCE,
    TYPE,
    SOURCE,
    BOND,
    AI,
    AJ,
    BONDORDER,
    COMPOSITION,
    COMPNAME,
    CATOM,
    C_NAME,
    C_NUMBER,
    DATASOURCE,
    PROGRAM,
    BASISSET,
    JOBTYPE,
    CONFORMATION,
    DATAFILE,
    UNIT,
    ATOM,
    ATOMID,
    OBTYPE,
    X_UNIT,
    V_UNIT,
    ESPID,
    dX, dY, dZ, dV,
    qXX, qYY, qZZ, qXY, qXZ, qYZ,
    oXXX, oXXY, oXXZ, oXYY, oXYZ, oXZZ, oYYY, oYYZ, oYZZ, oZZZ,
    hXXXX, hXXXY, hXXXZ, hXXYY, hXXYZ, hXXZZ, hXYYY, hXYYZ, hXYZZ, hXZZZ, hYYYY, hYYYZ, hYYZZ, hYZZZ, hZZZZ,
    aQ
};

const std::map<MolPropXml, MolPropObservable> xoMap = {
    { MolPropXml::DIPOLE,       MolPropObservable::DIPOLE       },
    { MolPropXml::QUADRUPOLE,   MolPropObservable::QUADRUPOLE   },
    { MolPropXml::OCTUPOLE,     MolPropObservable::OCTUPOLE     },
    { MolPropXml::HEXADECAPOLE, MolPropObservable::HEXADECAPOLE }
};

std::map<const std::string, MolPropXml> xmlxxx =
{
    { "molecules",        MolPropXml::MOLECULES     },
    { "molecule",         MolPropXml::MOLECULE      },
    { "formula",          MolPropXml::FORMULA       },
    { "molname",          MolPropXml::MOLNAME       },
    { "mass",             MolPropXml::MASS          },
    { "molinfo",          MolPropXml::MOLINFO       },
    { "iupac",            MolPropXml::IUPAC         },
    { "cas",              MolPropXml::CAS           },
    { "cid",              MolPropXml::CID           },
    { "inchi",            MolPropXml::INCHI         },
    { "multiplicity",     MolPropXml::MULTIPLICITY  },
    { "charge",           MolPropXml::CHARGE        },
    { "category",         MolPropXml::CATEGORY      },
    { "catname",          MolPropXml::CATNAME       },
    { "experiment",       MolPropXml::EXPERIMENT    },
    { "polarizability",   MolPropXml::POLARIZABILITY},
    { "energy",           MolPropXml::ENERGY        },
    { "potential",        MolPropXml::POTENTIAL     },
    { "name",             MolPropXml::NAME          },
    { "average",          MolPropXml::AVERAGE       },
    { "error",            MolPropXml::ERROR         },
    { "temperature",      MolPropXml::TEMPERATURE   },
    { "phase",            MolPropXml::PHASE         },
    { "method",           MolPropXml::METHOD        },
    { "reference",        MolPropXml::REFERENCE     },
    { "type",             MolPropXml::TYPE          },
    { "source",           MolPropXml::SOURCE        },
    { "bond",             MolPropXml::BOND          },
    { "ai",               MolPropXml::AI            },
    { "aj",               MolPropXml::AJ            },
    { "bondorder",        MolPropXml::BONDORDER     },
    { "composition",      MolPropXml::COMPOSITION   },
    { "compname",         MolPropXml::COMPNAME      },
    { "catom",            MolPropXml::CATOM         },
    { "cname",            MolPropXml::C_NAME        },
    { "cnumber",          MolPropXml::C_NUMBER      },
    { "datasource",       MolPropXml::DATASOURCE    },
    { "program",          MolPropXml::PROGRAM       },
    { "basisset",         MolPropXml::BASISSET      },
    { "jobtype",          MolPropXml::JOBTYPE       },
    { "conformation",     MolPropXml::CONFORMATION  },
    { "datafile",         MolPropXml::DATAFILE      },
    { "unit",             MolPropXml::UNIT          },
    { "atom",             MolPropXml::ATOM          },
    { "atomid",           MolPropXml::ATOMID        },
    { "obtype",           MolPropXml::OBTYPE        },
    { "coord_unit",       MolPropXml::X_UNIT        },
    { "potential_unit",   MolPropXml::V_UNIT        },
    { "espid",            MolPropXml::ESPID         },
    { "V",                MolPropXml::dV            },
    { "q",                MolPropXml::aQ            }
};

std::map<MolPropXml, const std::string> rmap = {};

static void fillMaps()
{
    rmap.clear();
    for (auto &xo : xoMap)
    {
        xmlxxx.insert({ mpo_name(xo.second), xo.first });
    }
    xmlxxx.insert({ multipoleName({ XX }), MolPropXml::dX });
    xmlxxx.insert({ multipoleName({ YY }), MolPropXml::dY });
    xmlxxx.insert({ multipoleName({ ZZ }), MolPropXml::dZ });
    xmlxxx.insert({ multipoleName({ XX, XX }), MolPropXml::qXX });
    xmlxxx.insert({ multipoleName({ XX, YY }), MolPropXml::qXY });
    xmlxxx.insert({ multipoleName({ XX, ZZ }), MolPropXml::qXZ });
    xmlxxx.insert({ multipoleName({ YY, YY }), MolPropXml::qYY });
    xmlxxx.insert({ multipoleName({ YY, ZZ }), MolPropXml::qYZ });
    xmlxxx.insert({ multipoleName({ ZZ, ZZ }), MolPropXml::qZZ });
    xmlxxx.insert({ multipoleName({ XX, XX, XX }), MolPropXml::oXXX });
    xmlxxx.insert({ multipoleName({ XX, XX, YY }), MolPropXml::oXXY });
    xmlxxx.insert({ multipoleName({ XX, XX, ZZ }), MolPropXml::oXXZ });
    xmlxxx.insert({ multipoleName({ XX, YY, YY }), MolPropXml::oXYY });
    xmlxxx.insert({ multipoleName({ XX, YY, ZZ }), MolPropXml::oXYZ });
    xmlxxx.insert({ multipoleName({ XX, ZZ, ZZ }), MolPropXml::oXZZ });
    xmlxxx.insert({ multipoleName({ YY, YY, YY }), MolPropXml::oYYY });
    xmlxxx.insert({ multipoleName({ YY, YY, ZZ }), MolPropXml::oYYZ });
    xmlxxx.insert({ multipoleName({ YY, ZZ, ZZ }), MolPropXml::oYZZ });
    xmlxxx.insert({ multipoleName({ ZZ, ZZ, ZZ }), MolPropXml::oZZZ });
    xmlxxx.insert({ multipoleName({ XX, XX, XX, XX }), MolPropXml::hXXXX });
    xmlxxx.insert({ multipoleName({ XX, XX, XX, YY }), MolPropXml::hXXXY });
    xmlxxx.insert({ multipoleName({ XX, XX, XX, ZZ }), MolPropXml::hXXXZ });
    xmlxxx.insert({ multipoleName({ XX, XX, YY, YY }), MolPropXml::hXXYY });
    xmlxxx.insert({ multipoleName({ XX, XX, YY, ZZ }), MolPropXml::hXXYZ });
    xmlxxx.insert({ multipoleName({ XX, XX, ZZ, ZZ }), MolPropXml::hXXZZ });
    xmlxxx.insert({ multipoleName({ XX, YY, YY, YY }), MolPropXml::hXYYY });
    xmlxxx.insert({ multipoleName({ XX, YY, YY, ZZ }), MolPropXml::hXYYZ });
    xmlxxx.insert({ multipoleName({ XX, YY, ZZ, ZZ }), MolPropXml::hXYZZ });
    xmlxxx.insert({ multipoleName({ XX, ZZ, ZZ, ZZ }), MolPropXml::hXZZZ });
    xmlxxx.insert({ multipoleName({ YY, YY, YY, YY }), MolPropXml::hYYYY });
    xmlxxx.insert({ multipoleName({ YY, YY, YY, ZZ }), MolPropXml::hYYYZ });
    xmlxxx.insert({ multipoleName({ YY, YY, ZZ, ZZ }), MolPropXml::hYYZZ });
    xmlxxx.insert({ multipoleName({ YY, ZZ, ZZ, ZZ }), MolPropXml::hYZZZ });
    xmlxxx.insert({ multipoleName({ ZZ, ZZ, ZZ, ZZ }), MolPropXml::hZZZZ });

    if (rmap.empty())
    {
        for (auto iter = xmlxxx.begin(); iter != xmlxxx.end(); ++iter)
        {
            rmap.insert({iter->second, iter->first});
        }
    }
}

static bool NN(const std::map<MolPropXml, std::string> *xbuf, MolPropXml index)
{
    auto ptr = xbuf->find(index);
    return (xbuf->end() != ptr) && !ptr->second.empty();
}

static bool xmlFound(const std::map<MolPropXml, std::string> *xbuf, 
                     const std::vector<MolPropXml>           &index)
{
    bool found = true;
    for(auto &ind : index)
    {
        auto ptr = xbuf->find(ind);
        found = found && (xbuf->end() != ptr) && !ptr->second.empty();
    }
    return found;
}

static double xbuf_atof(const std::map<MolPropXml, std::string> *xbuf, MolPropXml index)
{
    auto ptr = xbuf->find(index);
    if (ptr != xbuf->end())
    {
        return my_atof(ptr->second.c_str(), rmap[index].c_str());
    }
    return my_atof("", rmap[index].c_str());
}

static int xbuf_atoi(const std::map<MolPropXml, std::string> *xbuf, MolPropXml index)
{
    auto ptr = xbuf->find(index);
    if (ptr != xbuf->end())
    {
        return my_atoi(ptr->second.c_str(), rmap[index].c_str());
    }
    return my_atoi("", rmap[index].c_str());
}

static void get_attributes(FILE                              *fp, 
                           xmlAttrPtr                         attr,
                           std::map<MolPropXml, std::string> *xbuf)
{
    while (attr != nullptr)
    {
        std::string attrname(reinterpret_cast<const char *>(attr->name));
        std::string attrval(reinterpret_cast<const char *>(attr->children->content));

        auto  iter = xmlxxx.find(attrname);
        if (iter != xmlxxx.end() && !attrval.empty())
        {
            auto ptr = xbuf->find(iter->second);
            if (xbuf->end() == ptr)
            {
                xbuf->insert({iter->second, attrval});
            }
            else
            {
                ptr->second.assign(attrval);
            }
        }
        if (fp)
        {
            fprintf(fp, "Property: '%s' Value: '%s'\n",
                    attrname.c_str(), attrval.c_str());
        }
        attr = attr->next;
    }
}

static void get_molecule_attributes(std::map<MolPropXml, std::string> *xbuf,
                                    MolProp                           *mp)
{
    auto ff = xbuf->find(MolPropXml::FORMULA);
    if (xbuf->end() != ff)
    {
        mp->SetFormula(ff->second);
        xbuf->erase(ff);
    }
    ff = xbuf->find(MolPropXml::MOLNAME);
    if (xbuf->end() != ff)
    {
        mp->SetMolname(ff->second);
        xbuf->erase(ff);
    }
    ff = xbuf->find(MolPropXml::MASS);
    if (xbuf->end() != ff)
    {
        mp->SetMass(my_atof(ff->second.c_str(), rmap[ff->first].c_str()));
        xbuf->erase(ff);
    }
    ff = xbuf->find(MolPropXml::CHARGE);
    if (xbuf->end() != ff)
    {
        mp->SetTotalCharge(atoi(ff->second.c_str()));
        xbuf->erase(ff);
    }
    ff = xbuf->find(MolPropXml::MULTIPLICITY);
    if (xbuf->end() != ff)
    {
        mp->SetMultiplicity(atoi(ff->second.c_str()));
        xbuf->erase(ff);
    }
}

static void get_molinfo_attributes(std::map<MolPropXml, std::string> *xbuf,
                                   MolProp                           *mp)
{
    auto ff = xbuf->find(MolPropXml::IUPAC);
    if (xbuf->end() != ff)
    {
        mp->SetIupac(ff->second);
        xbuf->erase(ff);
    }
    ff = xbuf->find(MolPropXml::CAS);
    if (xbuf->end() != ff)
    {
        mp->SetCas(ff->second);
        xbuf->erase(ff);
    }
    ff = xbuf->find(MolPropXml::CID);
    if (xbuf->end() != ff)
    {
        mp->SetCid(ff->second);
        xbuf->erase(ff);
    }
    ff = xbuf->find(MolPropXml::INCHI);
    if (xbuf->end() != ff)
    {
        mp->SetInchi(ff->second);
        xbuf->erase(ff);
    }
}

static void clean_xbuf(std::map<MolPropXml, std::string> *xbuf,
                       const std::vector<MolPropXml>     &clean)
{
    for (auto &mpx : clean)
    {
        auto ff = xbuf->find(mpx);
        if (xbuf->end() != ff)
        {
            xbuf->erase(ff);
        }
    }
}

static void get_polarizability(std::map<MolPropXml, std::string> *xbuf,
                               Experiment                        *last)
{
    std::string mytype("electronic");
    if (last->dataSource() == dsExperiment)
    {
        mytype.assign("experiment");
    }
    auto mdp = new MolecularPolarizability(mytype,
                                           (*xbuf)[MolPropXml::UNIT],
                                           xbuf_atof(xbuf, MolPropXml::TEMPERATURE),
                                           xbuf_atof(xbuf, MolPropXml::qXX),
                                           xbuf_atof(xbuf, MolPropXml::qYY),
                                           xbuf_atof(xbuf, MolPropXml::qZZ),
                                           xbuf_atof(xbuf, MolPropXml::qXY),
                                           xbuf_atof(xbuf, MolPropXml::qXZ),
                                           xbuf_atof(xbuf, MolPropXml::qYZ),
                                           xbuf_atof(xbuf, MolPropXml::AVERAGE),
                                           xbuf_atof(xbuf, MolPropXml::ERROR));
    last->addProperty(MolPropObservable::POLARIZABILITY, mdp);
    clean_xbuf(xbuf,
               { MolPropXml::TEMPERATURE, MolPropXml::qXX, MolPropXml::qYY, MolPropXml::qZZ,
                 MolPropXml::qXY, MolPropXml::qXZ, MolPropXml::qYZ, MolPropXml::AVERAGE,
                 MolPropXml::UNIT, MolPropXml::ERROR });
}

static void mp_process_tree(FILE                              *fp, 
                            xmlNodePtr                         tree,
                            std::vector<MolProp>              *molprops,
                            std::map<MolPropXml, std::string> *xbuf)
{
    std::string qm_type("electronic");
    std::string exp_type("experiment");
    while (tree != nullptr)
    {
        if (fp)
        {
            if ((tree->type > 0) && ((unsigned)tree->type < NXMLTYPES))
            {
                fprintf(fp, "Node type %s encountered with name %s\n",
                        xmltypes[tree->type], (char *)tree->name);
            }
            else
            {
                fprintf(fp, "Node type %d encountered\n", tree->type);
            }
        }
        if (tree->type == XML_ELEMENT_NODE)
        {
            auto iter = xmlxxx.find((const char *)tree->name);
            if (iter != xmlxxx.end())
            {
                MolPropXml elem = iter->second;
                if (fp)
                {
                    fprintf(fp, "Element node name %s\n", (char *)tree->name);
                }
                if (tree->children && tree->children->content)
                {
                    xbuf->insert({ elem, (const char *)tree->children->content });
                } 
                get_attributes(fp, tree->properties, xbuf);
                
                MolProp *mpt = nullptr;
                if (molprops->size() > 0)
                {
                    mpt = &(molprops->back());
                }
                Experiment *last = nullptr;
                if (nullptr != mpt)
                {
                    last = mpt->LastExperiment();
                }
                switch (elem)
                {
                case MolPropXml::MOLECULES:
                    mp_process_tree(fp, tree->children, molprops, xbuf);
                    clean_xbuf(xbuf, { elem });
                    break;
                case MolPropXml::MOLECULE:
                    {
                        MolProp mp;
                        get_molecule_attributes(xbuf, &mp);
                        molprops->push_back(mp);
                        mpt = &(molprops->back());
                        mp_process_tree(fp, tree->children, molprops, xbuf);
                        clean_xbuf(xbuf, { elem });
                    }
                    break;
                case MolPropXml::MOLINFO:
                    get_molinfo_attributes(xbuf, mpt);
                    clean_xbuf(xbuf, { elem });
                    break;
                case MolPropXml::CATEGORY:
                    {
                        auto ff = xbuf->find(MolPropXml::CATNAME);
                        if (xbuf->end() != ff)
                        {
                            mpt->AddCategory(ff->second);
                            xbuf->erase(ff);
                        }
                    }
                    clean_xbuf(xbuf, { elem });
                    break;
                case MolPropXml::POLARIZABILITY:
                    {
                        std::vector<MolPropXml> clean1 = {
                            MolPropXml::TYPE, MolPropXml::UNIT, MolPropXml::TEMPERATURE
                        };
                        if ((nullptr != last) && xmlFound(xbuf, clean1))
                        {
                            mp_process_tree(fp, tree->children, molprops, xbuf);
                            std::vector<MolPropXml> clean2 = {
                                MolPropXml::AVERAGE, MolPropXml::ERROR
                            };
                            std::vector<MolPropXml> clean3 = {
                                MolPropXml::qXX, MolPropXml::qYY, MolPropXml::qZZ
                            };
                            mp_process_tree(fp, tree->children, molprops, xbuf);
                            if (xmlFound(xbuf, clean2) || xmlFound(xbuf, clean3))
                            {
                                // This routine cleans up as well.
                                get_polarizability(xbuf, last);
                            }
                        }
                    }
                    clean_xbuf(xbuf, { elem });
                    break;
                case MolPropXml::POTENTIAL:
                    {
                        std::vector<MolPropXml> clean1 = {
                            MolPropXml::X_UNIT, MolPropXml::V_UNIT, MolPropXml::ESPID
                        };
                        if ((nullptr != last) && xmlFound(xbuf, clean1))
                        {
                            mp_process_tree(fp, tree->children, molprops, xbuf);
                            std::vector<MolPropXml> clean2 = {
                                MolPropXml::dX, MolPropXml::dY, MolPropXml::dZ, MolPropXml::dV
                            };
                            if (xmlFound(xbuf, clean2))
                            {
                                ElectrostaticPotential ep((*xbuf)[MolPropXml::X_UNIT],
                                                          (*xbuf)[MolPropXml::V_UNIT],
                                                          atoi((*xbuf)[MolPropXml::ESPID].c_str()),
                                                          xbuf_atof(xbuf, MolPropXml::dX),
                                                          xbuf_atof(xbuf, MolPropXml::dY),
                                                          xbuf_atof(xbuf, MolPropXml::dZ),
                                                          xbuf_atof(xbuf, MolPropXml::dV));
                                last->AddPotential(ep);
                                clean_xbuf(xbuf, clean2);
                            }
                            clean_xbuf(xbuf, clean1);
                        }
                    }
                    clean_xbuf(xbuf, { elem });
                    break;
                case MolPropXml::DIPOLE:
                case MolPropXml::QUADRUPOLE:
                case MolPropXml::OCTUPOLE:
                case MolPropXml::HEXADECAPOLE:
                    if ((nullptr != last) &&
                        NN(xbuf, MolPropXml::TYPE) &&
                        NN(xbuf, MolPropXml::UNIT) &&
                        NN(xbuf, MolPropXml::TEMPERATURE))
                    {
                        mp_process_tree(fp, tree->children, molprops, xbuf);
                        MolPropObservable mpo = xoMap.find(elem)->second;
                        auto mq = new MolecularMultipole((*xbuf)[MolPropXml::TYPE],
                                                         (*xbuf)[MolPropXml::UNIT],
                                                         xbuf_atof(xbuf, MolPropXml::TEMPERATURE),
                                                         mpo);
                        std::vector<MolPropXml> myclean;
                        for(auto &x : *xbuf)
                        {
                        
                            if (mq->hasId(rmap[x.first]))
                            {
                                mq->setValue(rmap[x.first], xbuf_atof(xbuf, x.first));
                                myclean.push_back(x.first);
                            }
                        }
                        last->addProperty(mpo, mq);
                        clean_xbuf(xbuf, { MolPropXml::TYPE, MolPropXml::UNIT, MolPropXml::TEMPERATURE });
                        clean_xbuf(xbuf, myclean);
                    }
                    clean_xbuf(xbuf, { elem });
                    break;
                case MolPropXml::BOND:
                    if (NN(xbuf, MolPropXml::AI) && 
                        NN(xbuf, MolPropXml::AJ) &&
                        NN(xbuf, MolPropXml::BONDORDER))
                    {
                        Bond b(xbuf_atoi(xbuf, MolPropXml::AI)-1, 
                               xbuf_atoi(xbuf, MolPropXml::AJ)-1,
                               xbuf_atof(xbuf, MolPropXml::BONDORDER));
                        mpt->AddBond(b);
                        clean_xbuf(xbuf, { MolPropXml::AI, MolPropXml::AJ, MolPropXml::BONDORDER });
                    }
                    clean_xbuf(xbuf, { elem });
                    break;
                case MolPropXml::ENERGY:
                    if (!NN(xbuf, MolPropXml::ENERGY))
                    {
                        mp_process_tree(fp, tree->children, molprops, xbuf);
                    }   
                    if ((nullptr != last) &&
                        NN(xbuf, MolPropXml::TYPE)   && NN(xbuf, MolPropXml::UNIT) &&
                        NN(xbuf, MolPropXml::ENERGY) && NN(xbuf, MolPropXml::TEMPERATURE) &&
                        NN(xbuf, MolPropXml::PHASE))
                    {
                        MolPropObservable mpo;
                        if (stringToMolPropObservable((*xbuf)[MolPropXml::TYPE], &mpo))
                        {
                            std::string mytype(qm_type);
                            if (last->dataSource() == dsExperiment)
                            {
                                mytype = exp_type;
                            }
                            auto me  = new MolecularEnergy(mpo, mytype, (*xbuf)[MolPropXml::UNIT],
                                                           xbuf_atof(xbuf, MolPropXml::TEMPERATURE),
                                                           string2phase((*xbuf)[MolPropXml::PHASE]),
                                                           xbuf_atof(xbuf, MolPropXml::ENERGY),
                                                           xbuf_atof(xbuf, MolPropXml::ERROR));
                            last->addProperty(mpo, me);
                        }
                        else
                        {
                            fprintf(stderr, "Ignoring unknown property %s\n",
                                    (*xbuf)[MolPropXml::TYPE].c_str());
                        }
                        clean_xbuf(xbuf, { MolPropXml::TYPE, MolPropXml::UNIT,
                                          MolPropXml::TEMPERATURE, MolPropXml::PHASE });
                    }
                    clean_xbuf(xbuf, { elem });
                    break;

                case MolPropXml::ATOM:
                    {
                        std::vector<MolPropXml> clean0 = {
                            MolPropXml::NAME, MolPropXml::OBTYPE, MolPropXml::ATOMID
                        };
                        if ((nullptr != last) && xmlFound(xbuf, clean0))
                        {
                            CalcAtom ca((*xbuf)[MolPropXml::NAME],
                                        (*xbuf)[MolPropXml::OBTYPE],
                                        xbuf_atoi(xbuf, MolPropXml::ATOMID));
                            clean_xbuf(xbuf, clean0);
                            mp_process_tree(fp, tree->children, molprops, xbuf);
                        
                            std::vector<MolPropXml> clean1 = {
                                MolPropXml::dX, MolPropXml::dY, 
                                MolPropXml::dZ, MolPropXml::X_UNIT
                            };
                            if (xmlFound(xbuf, clean1))
                            {
                                ca.SetUnit((*xbuf)[MolPropXml::X_UNIT]);
                                ca.SetCoords(xbuf_atof(xbuf, MolPropXml::dX),
                                             xbuf_atof(xbuf, MolPropXml::dY),
                                             xbuf_atof(xbuf, MolPropXml::dZ));
                                clean_xbuf(xbuf, clean1);
                            }
                            std::vector<MolPropXml> clean2 = {
                                MolPropXml::aQ, MolPropXml::TYPE
                            };
                            if (xmlFound(xbuf, clean2))
                            {
                                ca.AddCharge(stringToQtype((*xbuf)[MolPropXml::TYPE]),
                                             xbuf_atof(xbuf, MolPropXml::aQ));
                                clean_xbuf(xbuf, clean2);
                            }
                            /* Now finally add the atom */
                            last->AddAtom(ca);
                        }
                    }
                    clean_xbuf(xbuf, { elem });
                    break;
                    
                case MolPropXml::EXPERIMENT:
                    if (NN(xbuf, MolPropXml::DATASOURCE))
                    {
                        DataSource ds = dataSourceFromName((*xbuf)[MolPropXml::DATASOURCE]);
                        const std::vector<MolPropXml> cleanme = {
                            MolPropXml::PROGRAM,      MolPropXml::REFERENCE,
                            MolPropXml::CONFORMATION, MolPropXml::DATAFILE,
                            MolPropXml::BASISSET,     MolPropXml::METHOD
                        };
                        if (ds == dsTheory && xmlFound(xbuf, cleanme))
                        {
                            Experiment mycalc((*xbuf)[MolPropXml::PROGRAM], (*xbuf)[MolPropXml::METHOD],
                                              (*xbuf)[MolPropXml::BASISSET], (*xbuf)[MolPropXml::REFERENCE],
                                              (*xbuf)[MolPropXml::CONFORMATION], (*xbuf)[MolPropXml::DATAFILE],
                                              string2jobType((*xbuf)[MolPropXml::JOBTYPE]));
                            mpt->AddExperiment(mycalc);
                            clean_xbuf(xbuf, cleanme);
                            clean_xbuf(xbuf, { MolPropXml::JOBTYPE });
                        }
                        else if (ds == dsExperiment)
                        {
                            if (NN(xbuf, MolPropXml::REFERENCE))
                            {
                                std::string conf("unknown");
                                if (NN(xbuf, MolPropXml::CONFORMATION))
                                {
                                    conf = (*xbuf)[MolPropXml::CONFORMATION];
                                }
                                Experiment myexp((*xbuf)[MolPropXml::REFERENCE], conf);
                                mpt->AddExperiment(myexp);
                                clean_xbuf(xbuf, { MolPropXml::CONFORMATION, MolPropXml::REFERENCE  });
                            }
                            else
                            {
                                gmx_fatal(FARGS, "Experimental data without reference");
                            }
                        }
                        mp_process_tree(fp, tree->children, molprops, xbuf);
                        clean_xbuf(xbuf, { MolPropXml::DATASOURCE });
                    }
                    clean_xbuf(xbuf, { elem });
                    break;
                default:
                    break;
                }
            }
        }
        tree = tree->next;
    }
}

void MolPropRead(const char *fn, std::vector<MolProp> *mpt)
{
    xmlDocPtr     doc;
    const char   *db          = "alexandria.ff/molprops.dat";
    std::string   mpfile;

    fillMaps();
    xmlDoValidityCheckingDefaultValue = 0;
    mpfile = gmx::findLibraryFile(fn ? fn : db, true, false);
    if (debug)
    {
        fprintf(debug, "Opening %s\n", mpfile.c_str());
    }
    print_memory_usage(debug);
    if ((doc = xmlParseFile(mpfile.c_str())) == nullptr)
    {
        gmx_fatal(FARGS, "Failed reading XML file %s. Run a syntax checker such as nsgmls.",
                  mpfile.c_str());
    }
    print_memory_usage(debug);
    if (nullptr != debug)
    {
        fprintf(debug, "Reading library file %s\n", fn);
    }
    std::map<MolPropXml, std::string>  xbuf;
    mp_process_tree(nullptr, doc->children, mpt, &xbuf);
    xmlFreeDoc(doc);
    print_memory_usage(debug);
}

static void add_exper_properties(xmlNodePtr                    exp,
                                 const Experiment &exper)
{
    xmlNodePtr child;

    for (auto &props : exper.propertyConst())
    {
        auto mpo = props.first;
        for (auto &prop : props.second)
        {
            std::string outUnit(mpo_unit2(mpo));
            double      fac = convertFromGromacs(1, outUnit);
            switch(mpo)
            {
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
                    double average = prop->getValue();
                    double error   = prop->getError();
                    child = add_xml_child_val(exp, rmap[MolPropXml::ENERGY], gmx_ftoa(average).c_str());
                    add_xml_string(child, rmap[MolPropXml::TYPE], prop->getType());
                    add_xml_string(child, rmap[MolPropXml::UNIT], outUnit);
                    add_xml_double(child, rmap[MolPropXml::TEMPERATURE], prop->getTemperature());
                    add_xml_string(child, rmap[MolPropXml::PHASE], phase2string(prop->getPhase()));
                    add_xml_child_val(child, rmap[MolPropXml::AVERAGE], gmx_ftoa(fac*average).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::ERROR], gmx_ftoa(fac*error).c_str());
                    break;
                }
            case MolPropObservable::DIPOLE:
            case MolPropObservable::QUADRUPOLE:
            case MolPropObservable::OCTUPOLE:
            case MolPropObservable::HEXADECAPOLE:
                {
                    auto mn = multipoleNames(mpo);
                    auto mq = prop->getVector();
                    
                    child = add_xml_child(exp, mpo_name(mpo));
                    add_xml_string(child, rmap[MolPropXml::TYPE], prop->getType());
                    add_xml_string(child, rmap[MolPropXml::UNIT], outUnit);
                    add_xml_double(child, rmap[MolPropXml::TEMPERATURE], prop->getTemperature());
                    for(size_t i = 0; i < mn.size(); i++)
                    {
                        add_xml_child_val(child, mn[i], gmx_ftoa(fac*mq[i]).c_str());
                    }
                    break;
                }
            case MolPropObservable::POLARIZABILITY:
                {
                    auto pprop = static_cast<MolecularPolarizability *>(prop);
                    auto mq    = pprop->getTensor();
                    
                    child = add_xml_child(exp, rmap[MolPropXml::POLARIZABILITY]);
                    add_xml_string(child, rmap[MolPropXml::TYPE], prop->getType());
                    add_xml_string(child, rmap[MolPropXml::UNIT], outUnit);
                    add_xml_double(child, rmap[MolPropXml::TEMPERATURE], prop->getTemperature());
                    add_xml_child_val(child, rmap[MolPropXml::AVERAGE], gmx_ftoa(fac*prop->getValue()).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::ERROR], gmx_ftoa(fac*prop->getError()).c_str());
                    if (mq[XX][XX] > 0 || mq[YY][YY] > 0 || mq[ZZ][ZZ])
                    {
                        add_xml_child_val(child, rmap[MolPropXml::qXX], gmx_ftoa(fac*mq[XX][XX]).c_str());
                        add_xml_child_val(child, rmap[MolPropXml::qYY], gmx_ftoa(fac*mq[YY][YY]).c_str());
                        add_xml_child_val(child, rmap[MolPropXml::qZZ], gmx_ftoa(fac*mq[ZZ][ZZ]).c_str());
                        add_xml_child_val(child, rmap[MolPropXml::qXY], gmx_ftoa(fac*mq[XX][YY]).c_str());
                        add_xml_child_val(child, rmap[MolPropXml::qXZ], gmx_ftoa(fac*mq[XX][YY]).c_str());
                        add_xml_child_val(child, rmap[MolPropXml::qYZ], gmx_ftoa(fac*mq[YY][ZZ]).c_str());
                    }
                    break;
                }
            case MolPropObservable::POTENTIAL:
            case MolPropObservable::CHARGE:
            case MolPropObservable::COORDINATES:
                break;
            }
        }
    }
}

static void add_calc_properties(xmlNodePtr                    exp,
                                const Experiment &calc)
{
    for (auto &ep : calc.electrostaticPotentialConst())
    {
        std::string x_unit, v_unit;
        double      x, y, z, V;
        int         espid;

        ep.get(&x_unit, &v_unit, &espid, &x, &y, &z, &V);
        std::string x_unitOut("pm");
        std::string V_unitOut("Hartree/e");
        double xfac = convertFromGromacs(1.0, x_unitOut);
        double Vfac = convertFromGromacs(1.0, V_unitOut);
        xmlNodePtr child = add_xml_child(exp, rmap[MolPropXml::POTENTIAL]);
        add_xml_char(child, rmap[MolPropXml::X_UNIT], x_unitOut.c_str());
        add_xml_char(child, rmap[MolPropXml::V_UNIT], V_unitOut.c_str());
        add_xml_int(child, rmap[MolPropXml::ESPID], espid);
        if ((x != 0) || (y != 0) || (z != 0) || (V != 0))
        {
            add_xml_child_val(child, rmap[MolPropXml::dX], gmx::formatString("%g", xfac*x).c_str());
            add_xml_child_val(child, rmap[MolPropXml::dY], gmx::formatString("%g", xfac*y).c_str());
            add_xml_child_val(child, rmap[MolPropXml::dZ], gmx::formatString("%g", xfac*z).c_str());
            add_xml_child_val(child, rmap[MolPropXml::dV], gmx::formatString("%g", Vfac*V).c_str());
        }
    }
}

static void add_xml_molprop(xmlNodePtr                 parent,
                            const MolProp &mp)
{
    xmlNodePtr ptr = add_xml_child(parent, rmap[MolPropXml::MOLECULE]);
    add_xml_string(ptr, rmap[MolPropXml::MOLNAME], mp.getMolname());
    add_xml_string(ptr, rmap[MolPropXml::FORMULA], mp.formula());
    add_xml_double(ptr, rmap[MolPropXml::MASS], mp.getMass());
    add_xml_double(ptr, rmap[MolPropXml::CHARGE], mp.totalCharge());
    add_xml_double(ptr, rmap[MolPropXml::MULTIPLICITY], mp.getMultiplicity());

    xmlNodePtr child = add_xml_child(ptr, rmap[MolPropXml::MOLINFO]);
    add_xml_string(child, rmap[MolPropXml::IUPAC], mp.getIupac());
    add_xml_string(child, rmap[MolPropXml::CAS], mp.getCas());
    add_xml_string(child, rmap[MolPropXml::CID], mp.getCid());
    add_xml_string(child, rmap[MolPropXml::INCHI], mp.getInchi());

    for (auto &b : mp.bondsConst())
    {
        xmlNodePtr child = add_xml_child(ptr, rmap[MolPropXml::BOND]);
        add_xml_int(child, rmap[MolPropXml::AI], 1+b.aI());
        add_xml_int(child, rmap[MolPropXml::AJ], 1+b.aJ());
        add_xml_double(child, rmap[MolPropXml::BONDORDER], b.bondOrder());
    }

    for (auto &me : mp.experimentConst())
    {
        xmlNodePtr             child = add_xml_child(ptr, rmap[MolPropXml::EXPERIMENT]);
        DataSource ds    = me.dataSource();
        add_xml_string(child, rmap[MolPropXml::DATASOURCE], dataSourceName(ds));
        add_xml_string(child, rmap[MolPropXml::REFERENCE], me.getReference());
        add_xml_string(child, rmap[MolPropXml::CONFORMATION], me.getConformation());
        if (dsTheory == ds)
        {
            add_xml_string(child, rmap[MolPropXml::PROGRAM], me.getProgram());
            add_xml_string(child, rmap[MolPropXml::METHOD], me.getMethod());
            add_xml_string(child, rmap[MolPropXml::BASISSET], me.getBasisset());
            add_xml_string(child, rmap[MolPropXml::JOBTYPE], jobType2string(me.getJobtype()));
            add_xml_string(child, rmap[MolPropXml::DATAFILE], me.getDatafile());
        }

        add_exper_properties(child, me);
        add_calc_properties(child, me);

        for (auto &ca : me.calcAtomConst())
        {
            xmlNodePtr grandchild = add_xml_child(child, rmap[MolPropXml::ATOM]);
            add_xml_string(grandchild, rmap[MolPropXml::NAME], ca.getName());
            add_xml_string(grandchild, rmap[MolPropXml::OBTYPE], ca.getObtype());
            add_xml_int(grandchild, rmap[MolPropXml::ATOMID], ca.getAtomid());
            add_xml_string(grandchild, rmap[MolPropXml::X_UNIT], ca.getUnit());

            double x, y, z;
            ca.getCoords(&x, &y, &z);

            xmlNodePtr  baby = add_xml_child_val(grandchild, rmap[MolPropXml::dX],
                                                 gmx::formatString("%g", x).c_str());
            baby = add_xml_child_val(grandchild, rmap[MolPropXml::dY], gmx::formatString("%g", y).c_str());
            baby = add_xml_child_val(grandchild, rmap[MolPropXml::dZ], gmx::formatString("%g", z).c_str());

            for (auto &q : ca.chargesConst())
            {
                xmlNodePtr atomptr = add_xml_child_val(grandchild, rmap[MolPropXml::aQ], gmx::formatString("%g", q.second).c_str());
                add_xml_string(atomptr, rmap[MolPropXml::TYPE], qTypeName(q.first));
            }
        }
    }
    for (auto &s : mp.categoryConst())
    {
        xmlNodePtr child = add_xml_child(ptr, rmap[MolPropXml::CATEGORY]);
        add_xml_string(child, rmap[MolPropXml::CATNAME], s);
    }

}

void MolPropWrite(const char                 *fn,
                  const std::vector<MolProp> &mpt,
                  gmx_bool                    bCompress)
{
    xmlDocPtr                   doc;
    xmlDtdPtr                   dtd;
    xmlNodePtr                  myroot;
    xmlChar                    *libdtdname, *dtdname, *gmx;

    gmx        = (xmlChar *) "molecules";
    dtdname    = (xmlChar *) "molprops.dtd";
    libdtdname = dtdname;
    fillMaps();
    if ((doc = xmlNewDoc((xmlChar *)"1.0")) == nullptr)
    {
        gmx_fatal(FARGS, "Creating XML document %s", fn);
    }

    if ((dtd = xmlCreateIntSubset(doc, dtdname, libdtdname, dtdname)) == nullptr)
    {
        gmx_fatal(FARGS, "Creating XML DTD %s", fn);
    }

    if ((myroot = xmlNewDocNode(doc, nullptr, gmx, nullptr)) == nullptr)
    {
        gmx_fatal(FARGS, "Creating root element for %s", fn);
    }
    dtd->next    = myroot;
    myroot->prev = (xmlNodePtr) dtd;

    /* Add molecule definitions */
    for (auto &mp : mpt)
    {
        add_xml_molprop(myroot, mp);
    }
    xmlSetDocCompressMode(doc, (int)bCompress);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fn, doc, "ISO-8859-1", 1) == 0)
    {
        gmx_fatal(FARGS, "Saving file %s", fn);
    }
    xmlFreeDoc(doc);
}

} // namespace alexandria
