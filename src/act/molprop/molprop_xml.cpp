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

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "act/molprop/molprop.h"
#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
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
    dX,
    dY,
    dZ,
    dV,
    qXX,
    qYY,
    qZZ,
    qXY,
    qXZ,
    qYZ,
    aQ
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
    { "dipole",           MolPropXml::DIPOLE        },
    { "quadrupole",       MolPropXml::QUADRUPOLE    },
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
    { "x",                MolPropXml::dX            },
    { "y",                MolPropXml::dY            },
    { "z",                MolPropXml::dZ            },
    { "V",                MolPropXml::dV            },
    { "xx",               MolPropXml::qXX           },
    { "yy",               MolPropXml::qYY           },
    { "zz",               MolPropXml::qZZ           },
    { "xy",               MolPropXml::qXY           },
    { "xz",               MolPropXml::qXZ           },
    { "yz",               MolPropXml::qYZ           },
    { "q",                MolPropXml::aQ            }
};

std::map<MolPropXml, const std::string> rmap = {};

static void add_xml_string(xmlNodePtr ptr, const std::string &name, const std::string &val)
{
    if (xmlSetProp(ptr, xmlCharStrdup(name.c_str()), xmlCharStrdup(val.c_str())) == 0)
    {
        gmx_fatal(FARGS, "Setting %s", name.c_str());
    }
}

static bool NN(const std::map<MolPropXml, std::string> &xbuf, MolPropXml index)
{
    auto ptr = xbuf.find(index);
    return (xbuf.end() != ptr) && !ptr->second.empty();
}

static char *sp(int n, char buf[], int maxindent)
{
    int i;
    if (n >= maxindent)
    {
        n = maxindent-1;
    }

    /* Don't indent more than maxindent characters */
    for (i = 0; (i < n); i++)
    {
        buf[i] = ' ';
    }
    buf[i] = '\0';

    return buf;
}

static double xbuf_atof(const std::map<MolPropXml, std::string> &xbuf, MolPropXml index)
{
    auto ptr = xbuf.find(index);
    if (ptr != xbuf.end())
    {
        return my_atof(ptr->second.c_str(), rmap[index].c_str());
    }
    return my_atof("", rmap[index].c_str());
}

static void get_attributes(FILE                              *fp, 
                           gmx_bool                           bZero,
                           int                                indent,
                           xmlAttrPtr                         attr,
                           std::map<MolPropXml, std::string> *xbuf)
{
    if (bZero)
    {
        xbuf->clear();
    }

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
            char  buf[100];

            fprintf(fp, "%sProperty: '%s' Value: '%s'\n",
                    sp(indent, buf, sizeof(buf)-1),
                    attrname.c_str(), attrval.c_str());
        }
        attr = attr->next;
    }
}

static void process_children(xmlNodePtr tree, std::map<MolPropXml, std::string> *xbuf)
{
    while (nullptr != tree)
    {
        auto iter = xmlxxx.find((const char *)tree->name);
        if (iter != xmlxxx.end() &&
            (nullptr != tree->children) &&
            (nullptr != tree->children->content))
        {
            std::string content(reinterpret_cast<char *>(tree->children->content));
            auto node = xbuf->find(iter->second);
            if (xbuf->end() == node)
            {
                xbuf->insert({ iter->second, content });
            }
            else if (node->second.empty())
            {
                node->second.assign(content);
            }
            else if (debug)
            {
                fprintf(debug, "Multiple children with property '%s' value '%s'\n", tree->name,
                        content.c_str());
            }
        }
        tree = tree->next;
    }
}

static void mp_process_tree(FILE                 *fp, 
                            xmlNodePtr            tree,
                            int                   indent,
                            std::vector<MolProp> *molprops,
                            gmx_bool             *bExperiment)
{
    xmlNodePtr                         tc;
    char                               buf[100];
    MolProp                           *mpt;
    std::map<MolPropXml, std::string>  xbuf;
    MolPropXml                         node;
    std::string                        xxx;

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
        if (molprops->size() > 0)
        {
            mpt = &(molprops->back());
        }
        else
        {
            mpt = nullptr;
        }
        std::string qm_type("electronic");
        std::string exp_type("experiment");
        switch (tree->type)
        {
            case XML_TEXT_NODE:
            {
                //   fprintf(stderr,"Text node %s encountered\n",(char *)tree->name);
                break;
            }
            case XML_ELEMENT_NODE:
            {
                auto iter = xmlxxx.find((const char *)tree->name);
                if (iter != xmlxxx.end())
                {
                    auto elem = iter->second;
                    if (fp)
                    {
                        fprintf(fp, "%sElement node name %s\n", sp(indent, buf, 99),
                                (char *)tree->name);
                    }
                    get_attributes(fp, TRUE, indent, tree->properties, &xbuf);

                    Experiment *last = nullptr;
                    if (nullptr != mpt)
                    {
                        last = mpt->LastExperiment();
                    }
                    switch (elem)
                    {
                        case MolPropXml::MOLECULES:
                            break;
                        case MolPropXml::MOLECULE:
                        {
                            MolProp mp;
                            if (NN(xbuf, MolPropXml::FORMULA))
                            {
                                mp.SetFormula(xbuf[MolPropXml::FORMULA]);
                            }
                            if (NN(xbuf, MolPropXml::MOLNAME))
                            {
                                mp.SetMolname(xbuf[MolPropXml::MOLNAME]);
                            }
                            if (NN(xbuf, MolPropXml::MASS))
                            {
                                mp.SetMass(xbuf_atof(xbuf, MolPropXml::MASS));
                            }
                            if (NN(xbuf, MolPropXml::CHARGE))
                            {
                                mp.SetTotalCharge(atoi(xbuf[MolPropXml::CHARGE].c_str()));
                            }
                            if (NN(xbuf, MolPropXml::MULTIPLICITY))
                            {
                                mp.SetMultiplicity(atoi(xbuf[MolPropXml::MULTIPLICITY].c_str()));
                            }
                            molprops->push_back(mp);
                            mpt = &(molprops->back());
                        }
                        break;
                        /* The items below are handled when treating attributes */
                        case MolPropXml::MOLINFO:
                            if (NN(xbuf, MolPropXml::IUPAC))
                            {
                                mpt->SetIupac(xbuf[MolPropXml::IUPAC]);
                            }
                            if (NN(xbuf, MolPropXml::CAS))
                            {
                                mpt->SetCas(xbuf[MolPropXml::CAS]);
                            }
                            if (NN(xbuf, MolPropXml::CID))
                            {
                                mpt->SetCid(xbuf[MolPropXml::CID]);
                            }
                            if (NN(xbuf, MolPropXml::INCHI))
                            {
                                mpt->SetInchi(xbuf[MolPropXml::INCHI]);
                            }
                            break;
                        case MolPropXml::CATEGORY:
                            if (NN(xbuf, MolPropXml::CATNAME))
                            {
                                mpt->AddCategory(xbuf[MolPropXml::CATNAME]);
                            }
                            break;
                        case MolPropXml::POLARIZABILITY:
                            process_children(tree->children, &xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf, MolPropXml::TYPE)  && NN(xbuf, MolPropXml::UNIT) &&
                                NN(xbuf, MolPropXml::TEMPERATURE) &&
                                ((NN(xbuf, MolPropXml::AVERAGE) && NN(xbuf, MolPropXml::ERROR)) ||
                                 (NN(xbuf, MolPropXml::qXX) && NN(xbuf, MolPropXml::qYY) && NN(xbuf, MolPropXml::qZZ))))
                            {
                                std::string mytype(qm_type);
                                if (last->dataSource() == dsExperiment)
                                {
                                    mytype = exp_type;
                                }
                                auto mdp = new MolecularPolarizability(mytype,
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
                            }
                            break;
                        case MolPropXml::POTENTIAL:
                            process_children(tree->children, &xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf, MolPropXml::X_UNIT) && NN(xbuf, MolPropXml::V_UNIT) &&
                                NN(xbuf, MolPropXml::ESPID) &&
                                NN(xbuf, MolPropXml::dX) && NN(xbuf, MolPropXml::dY) &&
                                NN(xbuf, MolPropXml::dZ) && NN(xbuf, MolPropXml::dV))
                            {
                                ElectrostaticPotential ep(xbuf[MolPropXml::X_UNIT], xbuf[MolPropXml::V_UNIT],
                                                          atoi(xbuf[MolPropXml::ESPID].c_str()),
                                                          xbuf_atof(xbuf, MolPropXml::dX),
                                                          xbuf_atof(xbuf, MolPropXml::dY),
                                                          xbuf_atof(xbuf, MolPropXml::dZ),
                                                          xbuf_atof(xbuf, MolPropXml::dV));
                                last->AddPotential(ep);
                            }
                            break;
                        case MolPropXml::DIPOLE:
                            process_children(tree->children, &xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf, MolPropXml::TYPE) && NN(xbuf, MolPropXml::UNIT) &&
                                NN(xbuf, MolPropXml::AVERAGE) && NN(xbuf, MolPropXml::ERROR) &&
                                NN(xbuf, MolPropXml::TEMPERATURE))
                            {
                                std::string mytype(qm_type);
                                if (last->dataSource() == dsExperiment)
                                {
                                    mytype = exp_type;
                                }
                                auto mdp = new MolecularDipole(mytype,
                                                               xbuf_atof(xbuf, MolPropXml::TEMPERATURE),
                                                               xbuf_atof(xbuf, MolPropXml::dX),
                                                               xbuf_atof(xbuf, MolPropXml::dY),
                                                               xbuf_atof(xbuf, MolPropXml::dZ),
                                                               xbuf_atof(xbuf, MolPropXml::AVERAGE),
                                                               xbuf_atof(xbuf, MolPropXml::ERROR));
                                last->addProperty(MolPropObservable::DIPOLE, mdp);
                            }
                            break;
                        case MolPropXml::QUADRUPOLE:
                            process_children(tree->children, &xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf, MolPropXml::TYPE) && NN(xbuf, MolPropXml::UNIT) &&
                                NN(xbuf, MolPropXml::TEMPERATURE) &&
                                NN(xbuf, MolPropXml::qXX) && NN(xbuf, MolPropXml::qYY) && NN(xbuf, MolPropXml::qZZ) &&
                                NN(xbuf, MolPropXml::qXY) && NN(xbuf, MolPropXml::qXZ) && NN(xbuf, MolPropXml::qYZ))
                            {
                                std::string mytype(qm_type);
                                if (last->dataSource() == dsExperiment)
                                {
                                    mytype = exp_type;
                                }
                                auto mq = new MolecularQuadrupole(mytype,
                                                                  xbuf_atof(xbuf, MolPropXml::TEMPERATURE),
                                                                  xbuf_atof(xbuf, MolPropXml::qXX),
                                                                  xbuf_atof(xbuf, MolPropXml::qYY),
                                                                  xbuf_atof(xbuf, MolPropXml::qZZ),
                                                                  xbuf_atof(xbuf, MolPropXml::qXY),
                                                                  xbuf_atof(xbuf, MolPropXml::qXZ),
                                                                  xbuf_atof(xbuf, MolPropXml::qYZ));
                                last->addProperty(MolPropObservable::QUADRUPOLE, mq);
                            }
                            break;
                        case MolPropXml::BOND:
                            process_children(tree->children, &xbuf);
                            if (NN(xbuf, MolPropXml::AI) && NN(xbuf, MolPropXml::AJ) &&
                                NN(xbuf, MolPropXml::BONDORDER))
                            {
                                Bond b(atoi(xbuf[MolPropXml::AI].c_str())-1, atoi(xbuf[MolPropXml::AJ].c_str())-1,
                                                   xbuf_atof(xbuf, MolPropXml::BONDORDER));
                                mpt->AddBond(b);
                            }
                            break;
                        case MolPropXml::ENERGY:
                            process_children(tree, &xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf, MolPropXml::TYPE) && NN(xbuf, MolPropXml::UNIT) &&
                                NN(xbuf, MolPropXml::ENERGY) && NN(xbuf, MolPropXml::TEMPERATURE) &&
                                NN(xbuf, MolPropXml::PHASE))
                            {
                                MolPropObservable mpo;
                                if (stringToMolPropObservable(xbuf[MolPropXml::TYPE], &mpo))
                                {
                                    std::string mytype(qm_type);
                                    if (last->dataSource() == dsExperiment)
                                    {
                                        mytype = exp_type;
                                    }
                                    auto me  = new MolecularEnergy(mpo, mytype,
                                                                   xbuf_atof(xbuf, MolPropXml::TEMPERATURE),
                                                                   string2phase(xbuf[MolPropXml::PHASE]),
                                                                   xbuf_atof(xbuf, MolPropXml::ENERGY),
                                                                   xbuf_atof(xbuf, MolPropXml::ERROR));
                                    last->addProperty(mpo, me);
                                }
                                else
                                {
                                    fprintf(stderr, "Ignoring unknown property %s\n", xbuf[MolPropXml::TYPE].c_str());
                                }
                                xbuf.erase(xbuf.find(MolPropXml::TYPE));
                                xbuf.erase(xbuf.find(MolPropXml::ENERGY));
                            }
                            break;

                        case MolPropXml::ATOM:
                            if ((nullptr != last) &&
                                NN(xbuf, MolPropXml::NAME) && NN(xbuf, MolPropXml::OBTYPE) && NN(xbuf, MolPropXml::ATOMID))
                            {
                                CalcAtom ca(xbuf[MolPropXml::NAME], xbuf[MolPropXml::OBTYPE],
                                                        atoi(xbuf[MolPropXml::ATOMID].c_str()));
                                xbuf[MolPropXml::NAME].clear();
                                xbuf[MolPropXml::OBTYPE].clear();
                                xbuf[MolPropXml::ATOMID].clear();
                                for (tc = tree->children; (nullptr != tc); tc = tc->next)
                                {
                                    get_attributes(fp, FALSE, indent, tc->properties, &xbuf);
                                    auto iter = xmlxxx.find((char *)tc->name);
                                    if (iter != xmlxxx.end() &&
                                        (nullptr != tc->children) &&
                                        (nullptr != tc->children->content))
                                    {
                                        node       = iter->second;
                                        xbuf[node].assign((char *)tc->children->content);
                                    }

                                    if (NN(xbuf, MolPropXml::dX) && NN(xbuf, MolPropXml::dY) && NN(xbuf, MolPropXml::dZ)
                                        && NN(xbuf, MolPropXml::UNIT))
                                    {
                                        ca.SetUnit(xbuf[MolPropXml::UNIT]);
                                        ca.SetCoords(xbuf_atof(xbuf, MolPropXml::dX),
                                                     xbuf_atof(xbuf, MolPropXml::dY),
                                                     xbuf_atof(xbuf, MolPropXml::dZ));
                                        xbuf[MolPropXml::dX].clear();
                                        xbuf[MolPropXml::dY].clear();
                                        xbuf[MolPropXml::dZ].clear();
                                        xbuf[MolPropXml::UNIT].clear();
                                    }
                                    if (NN(xbuf, MolPropXml::aQ) && NN(xbuf, MolPropXml::TYPE))
                                    {
                                        ca.AddCharge(stringToQtype(xbuf[MolPropXml::TYPE]),
                                                     xbuf_atof(xbuf, MolPropXml::aQ));
                                        xbuf[MolPropXml::aQ].clear();
                                        xbuf[MolPropXml::TYPE].clear();
                                    }
                                }
                                /* Now finally add the atom */
                                last->AddAtom(ca);
                            }
                            break;

                        case MolPropXml::EXPERIMENT:
                            if (NN(xbuf, MolPropXml::DATASOURCE))
                            {
                                DataSource ds = dataSourceFromName(xbuf[MolPropXml::DATASOURCE]);

                                if (ds == dsTheory &&
                                    NN(xbuf, MolPropXml::PROGRAM) && NN(xbuf, MolPropXml::REFERENCE) &&
                                    NN(xbuf, MolPropXml::CONFORMATION) && NN(xbuf, MolPropXml::DATAFILE) &&
                                    NN(xbuf, MolPropXml::BASISSET) && NN(xbuf, MolPropXml::METHOD))
                                {
                                    Experiment mycalc(xbuf[MolPropXml::PROGRAM], xbuf[MolPropXml::METHOD],
                                                      xbuf[MolPropXml::BASISSET], xbuf[MolPropXml::REFERENCE],
                                                      xbuf[MolPropXml::CONFORMATION], xbuf[MolPropXml::DATAFILE],
                                                      string2jobType(xbuf[MolPropXml::JOBTYPE]));
                                    mpt->AddExperiment(mycalc);
                                }
                                else if (ds == dsExperiment)
                                {
                                    if (NN(xbuf, MolPropXml::REFERENCE))
                                    {
                                        std::string conf("unknown");
                                        if (NN(xbuf, MolPropXml::CONFORMATION))
                                        {
                                            conf = xbuf[MolPropXml::CONFORMATION];
                                        }
                                        Experiment myexp(xbuf[MolPropXml::REFERENCE], conf);
                                        mpt->AddExperiment(myexp);
                                    }
                                    else
                                    {
                                        gmx_fatal(FARGS, "Experimental data without reference");
                                    }
                                }
                                break;
                            }
                        default:
                            break;
                    }
                }
                xbuf.clear();
                if (tree->children)
                {
                    auto iter = xmlxxx.find((char *)tree->name);
                    if (iter != xmlxxx.end())
                    {
                        mp_process_tree(fp, tree->children, indent+2,
                                        molprops, bExperiment);
                    }
                }
                break;
            }
            default:
                break;
        }
        tree = tree->next;
    }
}

void MolPropRead(const char *fn, std::vector<MolProp> *mpt)
{
    xmlDocPtr     doc;
    const char   *db          = "alexandria.ff/molprops.dat";
    gmx_bool      bExperiment = false;
    std::string   mpfile;

    rmap.clear();
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
    mp_process_tree(nullptr, 
                    doc->children, 
                    0,
                    mpt, 
                    &bExperiment);
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
            switch(mpo)
            {
            case MolPropObservable::OCTUPOLE:
            case MolPropObservable::HEXADECAPOLE:
                gmx_fatal(FARGS, "Please implement multipoles");
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
                    child = add_xml_child_val(exp, rmap[MolPropXml::ENERGY], gmx_dtoa(average).c_str());
                    add_xml_string(child, rmap[MolPropXml::TYPE], mpo_name(mpo));
                    add_xml_string(child, rmap[MolPropXml::UNIT], mpo_unit(mpo));
                    add_xml_double(child, rmap[MolPropXml::TEMPERATURE], prop->getTemperature());
                    add_xml_string(child, rmap[MolPropXml::PHASE], phase2string(prop->getPhase()));
                    add_xml_child_val(child, rmap[MolPropXml::AVERAGE], gmx_ftoa(average).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::ERROR], gmx_ftoa(error).c_str());
                    break;
                }
            case MolPropObservable::DIPOLE:
                {
                    double average = prop->getValue();
                    double error   = prop->getError();
                    auto dp = prop->getVector();
                    
                    child = add_xml_child(exp, rmap[MolPropXml::DIPOLE]);
                    add_xml_string(child, rmap[MolPropXml::TYPE], prop->getType());
                    add_xml_string(child, rmap[MolPropXml::UNIT], prop->getUnit());
                    add_xml_double(child, rmap[MolPropXml::TEMPERATURE], prop->getTemperature());
                    add_xml_child_val(child, rmap[MolPropXml::AVERAGE], gmx_ftoa(average).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::ERROR], gmx_ftoa(error).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::dX], gmx_ftoa(dp[XX]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::dY], gmx_ftoa(dp[YY]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::dZ], gmx_ftoa(dp[ZZ]).c_str());
                    break;
                }
            case MolPropObservable::QUADRUPOLE:
                {
                    auto mq = prop->getTensor();
                    
                    child = add_xml_child(exp, rmap[MolPropXml::QUADRUPOLE]);
                    add_xml_string(child, rmap[MolPropXml::TYPE], prop->getType());
                    add_xml_string(child, rmap[MolPropXml::UNIT], prop->getUnit());
                    add_xml_double(child, rmap[MolPropXml::TEMPERATURE], prop->getTemperature());
                    add_xml_child_val(child, rmap[MolPropXml::qXX], gmx_ftoa(mq[XX][XX]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qYY], gmx_ftoa(mq[YY][YY]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qZZ], gmx_ftoa(mq[ZZ][ZZ]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qXY], gmx_ftoa(mq[XX][YY]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qXZ], gmx_ftoa(mq[XX][YY]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qYZ], gmx_ftoa(mq[YY][ZZ]).c_str());
                    break;
                }
            case MolPropObservable::POLARIZABILITY:
                {
                    auto mq = prop->getTensor();
                    
                    child = add_xml_child(exp, rmap[MolPropXml::POLARIZABILITY]);
                    add_xml_string(child, rmap[MolPropXml::TYPE], prop->getType());
                    add_xml_string(child, rmap[MolPropXml::UNIT], prop->getUnit());
                    add_xml_double(child, rmap[MolPropXml::TEMPERATURE], prop->getTemperature());
                    add_xml_child_val(child, rmap[MolPropXml::AVERAGE], gmx_ftoa(prop->getValue()).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::ERROR], gmx_ftoa(prop->getError()).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qXX], gmx_ftoa(mq[XX][XX]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qYY], gmx_ftoa(mq[YY][YY]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qZZ], gmx_ftoa(mq[ZZ][ZZ]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qXY], gmx_ftoa(mq[XX][YY]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qXZ], gmx_ftoa(mq[XX][YY]).c_str());
                    add_xml_child_val(child, rmap[MolPropXml::qYZ], gmx_ftoa(mq[YY][ZZ]).c_str());
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

        xmlNodePtr child = add_xml_child(exp, rmap[MolPropXml::POTENTIAL]);
        add_xml_char(child, rmap[MolPropXml::X_UNIT], x_unit.c_str());
        add_xml_char(child, rmap[MolPropXml::V_UNIT], v_unit.c_str());
        add_xml_int(child, rmap[MolPropXml::ESPID], espid);
        if ((x != 0) || (y != 0) || (z != 0) || (V != 0))
        {
            add_xml_child_val(child, rmap[MolPropXml::dX], gmx::formatString("%g", x).c_str());
            add_xml_child_val(child, rmap[MolPropXml::dY], gmx::formatString("%g", y).c_str());
            add_xml_child_val(child, rmap[MolPropXml::dZ], gmx::formatString("%g", z).c_str());
            add_xml_child_val(child, rmap[MolPropXml::dV], gmx::formatString("%g", V).c_str());
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

            double x, y, z;
            ca.getCoords(&x, &y, &z);

            xmlNodePtr  baby = add_xml_child_val(grandchild, rmap[MolPropXml::dX], gmx::formatString("%g", x).c_str());
            add_xml_string(baby, rmap[MolPropXml::UNIT], ca.getUnit());
            baby = add_xml_child_val(grandchild, rmap[MolPropXml::dY], gmx::formatString("%g", y).c_str());
            add_xml_string(baby, rmap[MolPropXml::UNIT], ca.getUnit());
            baby = add_xml_child_val(grandchild, rmap[MolPropXml::dZ], gmx::formatString("%g", z).c_str());
            add_xml_string(baby, rmap[MolPropXml::UNIT], ca.getUnit());

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
    rmap.clear();
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
