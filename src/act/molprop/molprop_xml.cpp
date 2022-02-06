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

#include "gromacs/math/vectypes.h"
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

static bool NN(const std::string &s)
{
    return s.size() > 0;
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

enum {
    exmlMOLECULES      = 0,
    exmlMOLECULE       = 1,
    exmlFORMULA        = 2,
    exmlMOLNAME        = 3,
    exmlMASS           = 4,
    exmlMOLINFO        = 5,
    exmlIUPAC          = 6,
    exmlCAS            = 7,
    exmlCID            = 8,
    exmlINCHI          = 9,
    exmlMULTIPLICITY   = 10,
    exmlCHARGE         = 11,
    exmlCATEGORY       = 12,
    exmlCATNAME        = 13,
    exmlEXPERIMENT     = 14,
    exmlPOLARIZABILITY = 15,
    exmlENERGY         = 16,
    exmlDIPOLE         = 17,
    exmlQUADRUPOLE     = 18,
    exmlPOTENTIAL      = 19,
    exmlNAME           = 20,
    exmlAVERAGE        = 21,
    exmlERROR          = 22,
    exmlTEMPERATURE    = 23,
    exmlPHASE          = 24,
    exmlMETHOD         = 25,
    exmlREFERENCE      = 26,
    exmlTYPE           = 27,
    exmlSOURCE         = 28,
    exmlBOND           = 29,
    exmlAI             = 30,
    exmlAJ             = 31,
    exmlBONDORDER      = 32,
    exmlCOMPOSITION    = 33,
    exmlCOMPNAME       = 34,
    exmlCATOM          = 35,
    exmlC_NAME         = 36,
    exmlC_NUMBER       = 37,
    exmlDATASOURCE     = 38,
    exmlPROGRAM        = 39,
    exmlBASISSET       = 40,
    exmlJOBTYPE        = 41,
    exmlCONFORMATION   = 42,
    exmlDATAFILE       = 43,
    exmlUNIT           = 44,
    exmlATOM           = 45,
    exmlATOMID         = 46,
    exmlOBTYPE         = 47,
    exmlX_UNIT         = 48,
    exmlV_UNIT         = 49,
    exmlESPID          = 50,
    exmlX              = 51,
    exmlY              = 52,
    exmlZ              = 53,
    exmlV              = 54,
    exmlXX             = 55,
    exmlYY             = 56,
    exmlZZ             = 57,
    exmlXY             = 58,
    exmlXZ             = 59,
    exmlYZ             = 60,
    exmlQ              = 61,
    exmlNR             = 62
};

//static const char *exml_names[exmlNR] = {
std::map<const std::string, int> xmlxxx =
{
    { "molecules",        exmlMOLECULES     },
    { "molecule",         exmlMOLECULE      },
    { "formula",          exmlFORMULA       },
    { "molname",          exmlMOLNAME       },
    { "mass",             exmlMASS          },
    { "molinfo",          exmlMOLINFO       },
    { "iupac",            exmlIUPAC         },
    { "cas",              exmlCAS           },
    { "cid",              exmlCID           },
    { "inchi",            exmlINCHI         },
    { "multiplicity",     exmlMULTIPLICITY  },
    { "charge",           exmlCHARGE        },
    { "category",         exmlCATEGORY      },
    { "catname",          exmlCATNAME       },
    { "experiment",       exmlEXPERIMENT    },
    { "polarizability",   exmlPOLARIZABILITY},
    { "energy",           exmlENERGY        },
    { "dipole",           exmlDIPOLE        },
    { "quadrupole",       exmlQUADRUPOLE    },
    { "potential",        exmlPOTENTIAL     },
    { "name",             exmlNAME          },
    { "average",          exmlAVERAGE       },
    { "error",            exmlERROR         },
    { "temperature",      exmlTEMPERATURE   },
    { "phase",            exmlPHASE         },
    { "method",           exmlMETHOD        },
    { "reference",        exmlREFERENCE     },
    { "type",             exmlTYPE          },
    { "source",           exmlSOURCE        },
    { "bond",             exmlBOND          },
    { "ai",               exmlAI            },
    { "aj",               exmlAJ            },
    { "bondorder",        exmlBONDORDER     },
    { "composition",      exmlCOMPOSITION   },
    { "compname",         exmlCOMPNAME      },
    { "catom",            exmlCATOM         },
    { "cname",            exmlC_NAME        },
    { "cnumber",          exmlC_NUMBER      },
    { "datasource",       exmlDATASOURCE    },
    { "program",          exmlPROGRAM       },
    { "basisset",         exmlBASISSET      },
    { "jobtype",          exmlJOBTYPE       },
    { "conformation",     exmlCONFORMATION  },
    { "datafile",         exmlDATAFILE      },
    { "unit",             exmlUNIT          },
    { "atom",             exmlATOM          },
    { "atomid",           exmlATOMID        },
    { "obtype",           exmlOBTYPE        },
    { "coord_unit",       exmlX_UNIT        },
    { "potential_unit",   exmlV_UNIT        },
    { "espid",            exmlESPID         },
    { "x",                exmlX             },
    { "y",                exmlY             },
    { "z",                exmlZ             },
    { "V",                exmlV             },
    { "xx",               exmlXX            },
    { "yy",               exmlYY            },
    { "zz",               exmlZZ            },
    { "xy",               exmlXY            },
    { "xz",               exmlXZ            },
    { "yz",               exmlYZ            },
    { "q",                exmlQ             }
};

std::map<int, const std::string> rmap = {};

static const char *exml_names(int xml)
{
    if (rmap.empty())
    {
        for (auto iter = xmlxxx.begin(); iter != xmlxxx.end(); ++iter)
        {
            rmap.insert({iter->second, iter->first});
        }
    }
    auto iter = rmap.find(xml);
    if (iter != rmap.end())
    {
        return iter->second.c_str();
    }
    return nullptr;
}

static void add_xml_string(xmlNodePtr ptr, const char *name, std::string val)
{
    if (xmlSetProp(ptr, xmlCharStrdup(name), xmlCharStrdup(val.c_str())) == 0)
    {
        gmx_fatal(FARGS, "Setting %s", (char *)name);
    }
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

static double xbuf_atof(std::vector<std::string> xbuf, int xbuf_index)
{
    return my_atof(xbuf[xbuf_index].c_str(), rmap[xbuf_index].c_str());
}

static void get_attributes(FILE *fp, gmx_bool bZero, int indent, xmlAttrPtr attr,
                           std::vector<std::string> &xbuf)
{
    if (bZero)
    {
        for (auto &x : xbuf)
        {
            x.clear();
        }
    }

    while (attr != nullptr)
    {
        char *attrname = (char *)attr->name;
        char *attrval  = (char *)attr->children->content;

        auto  iter = xmlxxx.find(attrname);
        if (iter != xmlxxx.end())
        {
            if (attrval != nullptr)
            {
                xbuf[iter->second].assign(attrval);
            }
        }
        if (fp)
        {
            char  buf[100];

            fprintf(fp, "%sProperty: '%s' Value: '%s'\n",
                    sp(indent, buf, sizeof(buf)-1),
                    attrname, attrval);
        }
        attr = attr->next;
    }
}

static void process_children(xmlNodePtr tree, std::vector<std::string> &xbuf)
{
    int node;

    while (nullptr != tree)
    {
        auto iter = xmlxxx.find((const char *)tree->name);
        if (iter != xmlxxx.end() &&
            (nullptr != tree->children) &&
            (nullptr != tree->children->content))
        {
            node = iter->second;
            if (xbuf[node].size() == 0)
            {
                xbuf[node].assign(reinterpret_cast<char *>(tree->children->content));
            }
        }
        tree = tree->next;
    }
}

static void mp_process_tree(FILE *fp, xmlNodePtr tree,
                            int indent,
                            std::vector<MolProp> *molprops,
                            gmx_bool *bExperiment)
{
    xmlNodePtr                   tc;
    char                         buf[100];
    MolProp         *mpt;
    std::vector<std::string>     xbuf;
    int                          node;
    std::string                  xxx;

    xxx.clear();
    xbuf.resize(exmlNR, xxx);
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
                    get_attributes(fp, TRUE, indent, tree->properties, xbuf);

                    Experiment *last = nullptr;
                    if (nullptr != mpt)
                    {
                        last = mpt->LastExperiment();
                    }
                    switch (elem)
                    {
                        case exmlMOLECULES:
                            break;
                        case exmlMOLECULE:
                        {
                            MolProp mp;
                            if (NN(xbuf[exmlFORMULA]))
                            {
                                mp.SetFormula(xbuf[exmlFORMULA]);
                            }
                            if (NN(xbuf[exmlMOLNAME]))
                            {
                                mp.SetMolname(xbuf[exmlMOLNAME]);
                            }
                            if (NN(xbuf[exmlMASS]))
                            {
                                mp.SetMass(xbuf_atof(xbuf, exmlMASS));
                            }
                            if (NN(xbuf[exmlCHARGE]))
                            {
                                mp.SetTotalCharge(atoi(xbuf[exmlCHARGE].c_str()));
                            }
                            if (NN(xbuf[exmlMULTIPLICITY]))
                            {
                                mp.SetMultiplicity(atoi(xbuf[exmlMULTIPLICITY].c_str()));
                            }
                            molprops->push_back(mp);
                            mpt = &(molprops->back());
                        }
                        break;
                        /* The items below are handled when treating attributes */
                        case exmlMOLINFO:
                            if (NN(xbuf[exmlIUPAC]))
                            {
                                mpt->SetIupac(xbuf[exmlIUPAC]);
                            }
                            if (NN(xbuf[exmlCAS]))
                            {
                                mpt->SetCas(xbuf[exmlCAS]);
                            }
                            if (NN(xbuf[exmlCID]))
                            {
                                mpt->SetCid(xbuf[exmlCID]);
                            }
                            if (NN(xbuf[exmlINCHI]))
                            {
                                mpt->SetInchi(xbuf[exmlINCHI]);
                            }
                            break;
                        case exmlCATEGORY:
                            if (NN(xbuf[exmlCATNAME]))
                            {
                                mpt->AddCategory(xbuf[exmlCATNAME]);
                            }
                            break;
                        case exmlPOLARIZABILITY:
                            process_children(tree->children, xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf[exmlTYPE])  && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlTEMPERATURE]) &&
                                ((NN(xbuf[exmlAVERAGE]) && NN(xbuf[exmlERROR])) ||
                                 (NN(xbuf[exmlXX]) && NN(xbuf[exmlYY]) && NN(xbuf[exmlZZ]))))
                            {
                                std::string mytype(qm_type);
                                if (last->dataSource() == dsExperiment)
                                {
                                    mytype = exp_type;
                                }
                                auto mdp = new MolecularPolarizability(mytype,
                                                                       xbuf_atof(xbuf, exmlTEMPERATURE),
                                                                       xbuf_atof(xbuf, exmlXX),
                                                                       xbuf_atof(xbuf, exmlYY),
                                                                       xbuf_atof(xbuf, exmlZZ),
                                                                       xbuf_atof(xbuf, exmlXY),
                                                                       xbuf_atof(xbuf, exmlXZ),
                                                                       xbuf_atof(xbuf, exmlYZ),
                                                                       xbuf_atof(xbuf, exmlAVERAGE),
                                                                       xbuf_atof(xbuf, exmlERROR));
                                last->addProperty(MolPropObservable::POLARIZABILITY, mdp);
                            }
                            break;
                        case exmlPOTENTIAL:
                            process_children(tree->children, xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf[exmlX_UNIT]) && NN(xbuf[exmlV_UNIT]) &&
                                NN(xbuf[exmlESPID]) &&
                                NN(xbuf[exmlX]) && NN(xbuf[exmlY]) &&
                                NN(xbuf[exmlZ]) && NN(xbuf[exmlV]))
                            {
                                ElectrostaticPotential ep(xbuf[exmlX_UNIT], xbuf[exmlV_UNIT],
                                                          atoi(xbuf[exmlESPID].c_str()),
                                                          xbuf_atof(xbuf, exmlX),
                                                          xbuf_atof(xbuf, exmlY),
                                                          xbuf_atof(xbuf, exmlZ),
                                                          xbuf_atof(xbuf, exmlV));
                                last->AddPotential(ep);
                            }
                            break;
                        case exmlDIPOLE:
                            process_children(tree->children, xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlAVERAGE]) && NN(xbuf[exmlERROR]) &&
                                NN(xbuf[exmlTEMPERATURE]))
                            {
                                std::string mytype(qm_type);
                                if (last->dataSource() == dsExperiment)
                                {
                                    mytype = exp_type;
                                }
                                auto mdp = new MolecularDipole(mytype,
                                                               xbuf_atof(xbuf, exmlTEMPERATURE),
                                                               xbuf_atof(xbuf, exmlX),
                                                               xbuf_atof(xbuf, exmlY),
                                                               xbuf_atof(xbuf, exmlZ),
                                                               xbuf_atof(xbuf, exmlAVERAGE),
                                                               xbuf_atof(xbuf, exmlERROR));
                                last->addProperty(MolPropObservable::DIPOLE, mdp);
                            }
                            break;
                        case exmlQUADRUPOLE:
                            process_children(tree->children, xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlTEMPERATURE]) &&
                                NN(xbuf[exmlXX]) && NN(xbuf[exmlYY]) && NN(xbuf[exmlZZ]) &&
                                NN(xbuf[exmlXY]) && NN(xbuf[exmlXZ]) && NN(xbuf[exmlYZ]))
                            {
                                std::string mytype(qm_type);
                                if (last->dataSource() == dsExperiment)
                                {
                                    mytype = exp_type;
                                }
                                auto mq = new MolecularQuadrupole(mytype,
                                                                  xbuf_atof(xbuf, exmlTEMPERATURE),
                                                                  xbuf_atof(xbuf, exmlXX),
                                                                  xbuf_atof(xbuf, exmlYY),
                                                                  xbuf_atof(xbuf, exmlZZ),
                                                                  xbuf_atof(xbuf, exmlXY),
                                                                  xbuf_atof(xbuf, exmlXZ),
                                                                  xbuf_atof(xbuf, exmlYZ));
                                last->addProperty(MolPropObservable::QUADRUPOLE, mq);
                            }
                            break;
                        case exmlBOND:
                            process_children(tree->children, xbuf);
                            if (NN(xbuf[exmlAI]) && NN(xbuf[exmlAJ]) &&
                                NN(xbuf[exmlBONDORDER]))
                            {
                                Bond b(atoi(xbuf[exmlAI].c_str())-1, atoi(xbuf[exmlAJ].c_str())-1,
                                                   xbuf_atof(xbuf, exmlBONDORDER));
                                mpt->AddBond(b);
                            }
                            break;
                        case exmlENERGY:
                            process_children(tree, xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlENERGY]) && NN(xbuf[exmlTEMPERATURE]) &&
                                NN(xbuf[exmlPHASE]))
                            {
                                MolPropObservable mpo;
                                if (stringToMolPropObservable(xbuf[exmlTYPE], &mpo))
                                {
                                    std::string mytype(qm_type);
                                    if (last->dataSource() == dsExperiment)
                                    {
                                        mytype = exp_type;
                                    }
                                    auto me  = new MolecularEnergy(mpo, mytype,
                                                                   xbuf_atof(xbuf, exmlTEMPERATURE),
                                                                   string2phase(xbuf[exmlPHASE]),
                                                                   xbuf_atof(xbuf, exmlENERGY),
                                                                   xbuf_atof(xbuf, exmlERROR));
                                    last->addProperty(mpo, me);
                                }
                                else
                                {
                                    fprintf(stderr, "Ignoring unknown property %s\n", xbuf[exmlTYPE].c_str());
                                }
                            }
                            break;

                        case exmlATOM:
                            if ((nullptr != last) &&
                                NN(xbuf[exmlNAME]) && NN(xbuf[exmlOBTYPE]) && NN(xbuf[exmlATOMID]))
                            {
                                CalcAtom ca(xbuf[exmlNAME], xbuf[exmlOBTYPE],
                                                        atoi(xbuf[exmlATOMID].c_str()));
                                xbuf[exmlNAME].clear();
                                xbuf[exmlOBTYPE].clear();
                                xbuf[exmlATOMID].clear();
                                for (tc = tree->children; (nullptr != tc); tc = tc->next)
                                {
                                    get_attributes(fp, FALSE, indent, tc->properties, xbuf);
                                    auto iter = xmlxxx.find((char *)tc->name);
                                    if (iter != xmlxxx.end() &&
                                        (nullptr != tc->children) &&
                                        (nullptr != tc->children->content))
                                    {
                                        node       = iter->second;
                                        xbuf[node].assign((char *)tc->children->content);
                                    }

                                    if (NN(xbuf[exmlX]) && NN(xbuf[exmlY]) && NN(xbuf[exmlZ])
                                        && NN(xbuf[exmlUNIT]))
                                    {
                                        ca.SetUnit(xbuf[exmlUNIT]);
                                        ca.SetCoords(xbuf_atof(xbuf, exmlX),
                                                     xbuf_atof(xbuf, exmlY),
                                                     xbuf_atof(xbuf, exmlZ));
                                        xbuf[exmlX].clear();
                                        xbuf[exmlY].clear();
                                        xbuf[exmlZ].clear();
                                        xbuf[exmlUNIT].clear();
                                    }
                                    if (NN(xbuf[exmlQ]) && NN(xbuf[exmlTYPE]))
                                    {
                                        ca.AddCharge(stringToQtype(xbuf[exmlTYPE]),
                                                     xbuf_atof(xbuf, exmlQ));
                                        xbuf[exmlQ].clear();
                                        xbuf[exmlTYPE].clear();
                                    }
                                }
                                /* Now finally add the atom */
                                last->AddAtom(ca);
                            }
                            break;

                        case exmlEXPERIMENT:
                            if (NN(xbuf[exmlDATASOURCE]))
                            {
                                DataSource ds = dataSourceFromName(xbuf[exmlDATASOURCE]);

                                if (ds == dsTheory &&
                                    NN(xbuf[exmlPROGRAM]) && NN(xbuf[exmlREFERENCE]) &&
                                    NN(xbuf[exmlCONFORMATION]) && NN(xbuf[exmlDATAFILE]))
                                {
                                    Experiment mycalc(xbuf[exmlPROGRAM], xbuf[exmlMETHOD],
                                                      xbuf[exmlBASISSET], xbuf[exmlREFERENCE],
                                                      xbuf[exmlCONFORMATION], xbuf[exmlDATAFILE],
                                                      string2jobType(xbuf[exmlJOBTYPE]));
                                    mpt->AddExperiment(mycalc);
                                }
                                else if (ds == dsExperiment)
                                {
                                    if (NN(xbuf[exmlREFERENCE]))
                                    {
                                        const char            *unknown = "unknown";
                                        Experiment myexp(xbuf[exmlREFERENCE],
                                                         NN(xbuf[exmlCONFORMATION]) ? xbuf[exmlCONFORMATION] : unknown);
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
                for (auto &i : xbuf)
                {
                    if (NN(i))
                    {
                        i.clear();
                    }
                }
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
                    child = add_xml_child_val(exp, exml_names(exmlENERGY), gmx_dtoa(average).c_str());
                    add_xml_string(child, exml_names(exmlTYPE), mpo_name(mpo));
                    add_xml_string(child, exml_names(exmlUNIT), mpo_unit(mpo));
                    add_xml_double(child, exml_names(exmlTEMPERATURE), prop->getTemperature());
                    add_xml_string(child, exml_names(exmlPHASE), phase2string(prop->getPhase()));
                    add_xml_child_val(child, exml_names(exmlAVERAGE), gmx_ftoa(average).c_str());
                    add_xml_child_val(child, exml_names(exmlERROR), gmx_ftoa(error).c_str());
                    break;
                }
            case MolPropObservable::DIPOLE:
                {
                    double average = prop->getValue();
                    double error   = prop->getError();
                    auto dp = prop->getVector();
                    
                    child = add_xml_child(exp, exml_names(exmlDIPOLE));
                    add_xml_string(child, exml_names(exmlTYPE), prop->getType());
                    add_xml_string(child, exml_names(exmlUNIT), prop->getUnit());
                    add_xml_double(child, exml_names(exmlTEMPERATURE), prop->getTemperature());
                    add_xml_child_val(child, exml_names(exmlAVERAGE), gmx_ftoa(average).c_str());
                    add_xml_child_val(child, exml_names(exmlERROR), gmx_ftoa(error).c_str());
                    add_xml_child_val(child, exml_names(exmlX), gmx_ftoa(dp[XX]).c_str());
                    add_xml_child_val(child, exml_names(exmlY), gmx_ftoa(dp[YY]).c_str());
                    add_xml_child_val(child, exml_names(exmlZ), gmx_ftoa(dp[ZZ]).c_str());
                    break;
                }
            case MolPropObservable::QUADRUPOLE:
                {
                    auto mq = prop->getTensor();
                    
                    child = add_xml_child(exp, exml_names(exmlQUADRUPOLE));
                    add_xml_string(child, exml_names(exmlTYPE), prop->getType());
                    add_xml_string(child, exml_names(exmlUNIT), prop->getUnit());
                    add_xml_double(child, exml_names(exmlTEMPERATURE), prop->getTemperature());
                    add_xml_child_val(child, exml_names(exmlXX), gmx_ftoa(mq[XX][XX]).c_str());
                    add_xml_child_val(child, exml_names(exmlYY), gmx_ftoa(mq[YY][YY]).c_str());
                    add_xml_child_val(child, exml_names(exmlZZ), gmx_ftoa(mq[ZZ][ZZ]).c_str());
                    add_xml_child_val(child, exml_names(exmlXY), gmx_ftoa(mq[XX][YY]).c_str());
                    add_xml_child_val(child, exml_names(exmlXZ), gmx_ftoa(mq[XX][YY]).c_str());
                    add_xml_child_val(child, exml_names(exmlYZ), gmx_ftoa(mq[YY][ZZ]).c_str());
                    break;
                }
            case MolPropObservable::POLARIZABILITY:
                {
                    auto mq = prop->getTensor();
                    
                    child = add_xml_child(exp, exml_names(exmlPOLARIZABILITY));
                    add_xml_string(child, exml_names(exmlTYPE), prop->getType());
                    add_xml_string(child, exml_names(exmlUNIT), prop->getUnit());
                    add_xml_double(child, exml_names(exmlTEMPERATURE), prop->getTemperature());
                    add_xml_child_val(child, exml_names(exmlAVERAGE), gmx_ftoa(prop->getValue()).c_str());
                    add_xml_child_val(child, exml_names(exmlERROR), gmx_ftoa(prop->getError()).c_str());
                    add_xml_child_val(child, exml_names(exmlXX), gmx_ftoa(mq[XX][XX]).c_str());
                    add_xml_child_val(child, exml_names(exmlYY), gmx_ftoa(mq[YY][YY]).c_str());
                    add_xml_child_val(child, exml_names(exmlZZ), gmx_ftoa(mq[ZZ][ZZ]).c_str());
                    add_xml_child_val(child, exml_names(exmlXY), gmx_ftoa(mq[XX][YY]).c_str());
                    add_xml_child_val(child, exml_names(exmlXZ), gmx_ftoa(mq[XX][YY]).c_str());
                    add_xml_child_val(child, exml_names(exmlYZ), gmx_ftoa(mq[YY][ZZ]).c_str());
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

        xmlNodePtr child = add_xml_child(exp, exml_names(exmlPOTENTIAL));
        add_xml_char(child, exml_names(exmlX_UNIT), x_unit.c_str());
        add_xml_char(child, exml_names(exmlV_UNIT), v_unit.c_str());
        add_xml_int(child, exml_names(exmlESPID), espid);
        if ((x != 0) || (y != 0) || (z != 0) || (V != 0))
        {
            add_xml_child_val(child, exml_names(exmlX), gmx::formatString("%g", x).c_str());
            add_xml_child_val(child, exml_names(exmlY), gmx::formatString("%g", y).c_str());
            add_xml_child_val(child, exml_names(exmlZ), gmx::formatString("%g", z).c_str());
            add_xml_child_val(child, exml_names(exmlV), gmx::formatString("%g", V).c_str());
        }
    }
}

static void add_xml_molprop(xmlNodePtr                 parent,
                            const MolProp &mp)
{
    xmlNodePtr ptr = add_xml_child(parent, exml_names(exmlMOLECULE));
    add_xml_string(ptr, exml_names(exmlMOLNAME), mp.getMolname());
    add_xml_string(ptr, exml_names(exmlFORMULA), mp.formula());
    add_xml_double(ptr, exml_names(exmlMASS), mp.getMass());
    add_xml_double(ptr, exml_names(exmlCHARGE), mp.totalCharge());
    add_xml_double(ptr, exml_names(exmlMULTIPLICITY), mp.getMultiplicity());

    xmlNodePtr child = add_xml_child(ptr, exml_names(exmlMOLINFO));
    add_xml_string(child, exml_names(exmlIUPAC), mp.getIupac());
    add_xml_string(child, exml_names(exmlCAS), mp.getCas());
    add_xml_string(child, exml_names(exmlCID), mp.getCid());
    add_xml_string(child, exml_names(exmlINCHI), mp.getInchi());

    for (auto &b : mp.bondsConst())
    {
        xmlNodePtr child = add_xml_child(ptr, exml_names(exmlBOND));
        add_xml_int(child, exml_names(exmlAI), 1+b.aI());
        add_xml_int(child, exml_names(exmlAJ), 1+b.aJ());
        add_xml_double(child, exml_names(exmlBONDORDER), b.bondOrder());
    }

    for (auto &me : mp.experimentConst())
    {
        xmlNodePtr             child = add_xml_child(ptr, exml_names(exmlEXPERIMENT));
        DataSource ds    = me.dataSource();
        add_xml_string(child, exml_names(exmlDATASOURCE), dataSourceName(ds));
        add_xml_string(child, exml_names(exmlREFERENCE), me.getReference());
        add_xml_string(child, exml_names(exmlCONFORMATION), me.getConformation());
        if (dsTheory == ds)
        {
            add_xml_string(child, exml_names(exmlPROGRAM), me.getProgram());
            add_xml_string(child, exml_names(exmlMETHOD), me.getMethod());
            add_xml_string(child, exml_names(exmlBASISSET), me.getBasisset());
            add_xml_string(child, exml_names(exmlJOBTYPE), jobType2string(me.getJobtype()));
            add_xml_string(child, exml_names(exmlDATAFILE), me.getDatafile());
        }

        add_exper_properties(child, me);
        add_calc_properties(child, me);

        for (auto &ca : me.calcAtomConst())
        {
            xmlNodePtr grandchild = add_xml_child(child, exml_names(exmlATOM));
            add_xml_string(grandchild, exml_names(exmlNAME), ca.getName());
            add_xml_string(grandchild, exml_names(exmlOBTYPE), ca.getObtype());
            add_xml_int(grandchild, exml_names(exmlATOMID), ca.getAtomid());

            double x, y, z;
            ca.getCoords(&x, &y, &z);

            xmlNodePtr  baby = add_xml_child_val(grandchild, exml_names(exmlX), gmx::formatString("%g", x).c_str());
            add_xml_string(baby, exml_names(exmlUNIT), ca.getUnit());
            baby = add_xml_child_val(grandchild, exml_names(exmlY), gmx::formatString("%g", y).c_str());
            add_xml_string(baby, exml_names(exmlUNIT), ca.getUnit());
            baby = add_xml_child_val(grandchild, exml_names(exmlZ), gmx::formatString("%g", z).c_str());
            add_xml_string(baby, exml_names(exmlUNIT), ca.getUnit());

            for (auto &q : ca.chargesConst())
            {
                xmlNodePtr atomptr = add_xml_child_val(grandchild, exml_names(exmlQ), gmx::formatString("%g", q.second).c_str());
                add_xml_string(atomptr, exml_names(exmlTYPE), qTypeName(q.first));
            }
        }
    }
    for (auto &s : mp.categoryConst())
    {
        xmlNodePtr child = add_xml_child(ptr, exml_names(exmlCATEGORY));
        add_xml_string(child, exml_names(exmlCATNAME), s);
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
