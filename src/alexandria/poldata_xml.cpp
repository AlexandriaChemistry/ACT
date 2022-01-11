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

#include "poldata_xml.h"

#include <cstdlib>
#include <cstring>

#include <map>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "gromacs/topology/atoms.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"

#include "forcefieldparameter.h"
#include "forcefieldparameterlist.h"
#include "poldata.h"
#include "poldata_low.h"
#include "xml_util.h"

//extern int xmlDoValidityCheckingDefaultValue;
namespace alexandria
{

const char *xmltypes[] = {
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

//! The different entries to excpect in force field file
enum class xmlEntry {
    GENTOP,
    TIMESTAMP,
    CHECKSUM,
    REFERENCE,
    PARTICLETYPES,
    PARTICLETYPE, 
    COMB_RULE,
    NEXCL,
    VERSION,
    EPSILONR,
    DESC,
    GEOMETRY,
    NUMBONDS,
    VANDERWAALS,
    INTERACTION,
    IDENTIFIER,
    CANSWAP,
    KEY,
    TYPE,
    ATYPE,
    VALUE,
    OPTION,
    PARAMETERLIST,
    PARAMETER,
    UNCERTAINTY,
    MINIMUM,
    MAXIMUM,
    MUTABILITY,
    NONNEGATIVE,
    ATOM1,
    ATOM2,
    ATOM3,
    ATOM4,
    BONDORDER,
    PARAMS,
    REFVALUE,
    UNIT,
    NTRAIN,
    GT_VSITES,
    GT_VSITE,
    SYMMETRIC_CHARGES,
    SYM_CHARGE,
    CENTRAL,
    ATTACHED,
    NUMATTACH,
    MODEL,
    CHARGES,
    REF_CHARGES,
    ANGLE_UNIT,
    LENGTH_UNIT,
    DISTANCE,
    NCONTROLATOMS,
    NUMBER,
    VTYPE,
    ANGLE,
    FUNCTION
};

//! Map from string used in XML file to xmlEntry enum
std::map<const std::string, xmlEntry> xml_pd =
{
    { "gentop",                    xmlEntry::GENTOP           },
    { "timestamp",                 xmlEntry::TIMESTAMP        },
    { "checksum",                  xmlEntry::CHECKSUM         },
    { "reference",                 xmlEntry::REFERENCE        },
    { "particletypes",             xmlEntry::PARTICLETYPES    },
    { "particletype",              xmlEntry::PARTICLETYPE     },
    { "parameterlist",             xmlEntry::PARAMETERLIST    },
    { "parameter",                 xmlEntry::PARAMETER        },
    { "version",                   xmlEntry::VERSION          },
    { "function",                  xmlEntry::FUNCTION         },
    { "nexclusions",               xmlEntry::NEXCL            },
    { "epsilonr",                  xmlEntry::EPSILONR         },
    { "description",               xmlEntry::DESC             },
    //{ "aromatic",                  xmlEntry::AROMATIC         },
    { "geometry",                  xmlEntry::GEOMETRY         },
    { "numbonds",                  xmlEntry::NUMBONDS         },
    //{ "vdwparams",                 xmlEntry::VDWPARAMS        },
    { "vanderwaals",               xmlEntry::VANDERWAALS      },
    { "interaction",               xmlEntry::INTERACTION      },
    { "identifier",                xmlEntry::IDENTIFIER       },
    { "canswap",                   xmlEntry::CANSWAP          },
    { "key",                       xmlEntry::KEY              },
    { "type",                      xmlEntry::TYPE             },
    { "atype",                     xmlEntry::ATYPE            },
    { "value",                     xmlEntry::VALUE            },
    { "option",                    xmlEntry::OPTION           },
    { "uncertainty",               xmlEntry::UNCERTAINTY      },
    { "minimum",                   xmlEntry::MINIMUM          },
    { "maximum",                   xmlEntry::MAXIMUM          },
    { "nonnegative",               xmlEntry::NONNEGATIVE      },
    { "mutability",                xmlEntry::MUTABILITY       },
    { "bondorder",                 xmlEntry::BONDORDER        },
    { "params",                    xmlEntry::PARAMS           },
    { "refValue",                  xmlEntry::REFVALUE         },
    { "unit",                      xmlEntry::UNIT             },
    { "ntrain",                    xmlEntry::NTRAIN           },
    { "gt_vsites",                 xmlEntry::GT_VSITES        },
    { "gt_vsite",                  xmlEntry::GT_VSITE         },
    { "symmetric_charges",         xmlEntry::SYMMETRIC_CHARGES},
    { "sym_charge",                xmlEntry::SYM_CHARGE       },
    { "central",                   xmlEntry::CENTRAL          },
    { "attached",                  xmlEntry::ATTACHED         },
    { "numattach",                 xmlEntry::NUMATTACH        },
    { "model",                     xmlEntry::MODEL            },
    { "angle_unit",                xmlEntry::ANGLE_UNIT       },
    { "length_unit",               xmlEntry::LENGTH_UNIT      },
    { "distance",                  xmlEntry::DISTANCE         },
    { "ncontrolatoms",             xmlEntry::NCONTROLATOMS    },
    { "number",                    xmlEntry::NUMBER           },
    { "vtype",                     xmlEntry::VTYPE            },
    { "angle",                     xmlEntry::ANGLE            }
};

std::map<xmlEntry, const std::string> rmap_pd = {};

typedef std::map<xmlEntry, std::string> xmlBuffer;

static const char *exml_names(xmlEntry xml)
{
    if (rmap_pd.empty())
    {
        for (auto &iter : xml_pd)
        {
            rmap_pd.insert({iter.second, iter.first});
        }
    }
    return rmap_pd[xml].c_str();
}

static bool NNlow(xmlBuffer *xbuf, xmlEntry xml, bool obligatory)
{
    if (xbuf->find(xml) == xbuf->end())
    {
        std::string buf = gmx::formatString("Missing%s variable %d (%s)", obligatory ? " required" : "",
                                            static_cast<int>(xml), exml_names(xml));
        if (obligatory)
        {
            GMX_THROW(gmx::InvalidInputError(buf.c_str()));
        }
        else
        {
            fprintf(stderr, "%s\n", buf.c_str());
            return false;
        }
    }
    return true;
}

static bool NNobligatory(xmlBuffer *xbuf, xmlEntry xml)
{
    return NNlow(xbuf, xml, true);
}

static bool NN(xmlBuffer *xbuf, xmlEntry xml)
{
    return NNlow(xbuf, xml, false);
}

static void sp(int n, char buf[], int maxindent)
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
}

static double xbuf_atof(xmlBuffer *xbuf, xmlEntry  xbuf_index)
{
    auto xb = xbuf->find(xbuf_index);
    if (xb == xbuf->end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No such entry '%d' in xbuf",
                                                       static_cast<int>(xbuf_index)).c_str()));
    }
    auto rm = exml_names(xbuf_index);
    return my_atof(xb->second.c_str(), rm);
}

static void processAttr(FILE       *fp, 
                        xmlAttrPtr  attr,
                        xmlBuffer  *xbuf,
                        xmlEntry    elem,
                        int         indent, 
                        Poldata    *pd)
{
    std::string attrname, attrval;
    char        buf[100];
    //    Identifer   identifier;

    xbuf->clear();
    while (attr != nullptr)
    {
        attrname.assign((char *)attr->name);
        attrval.assign((char *)attr->children->content);

        auto iter = xml_pd.find(attrname);
        if (iter != xml_pd.end())
        {
            xbuf->insert({iter->second, attrval});

            if (nullptr != fp)
            {
                sp(indent, buf, 99);
                fprintf(fp, "%sProperty: '%s' Value: '%s'\n", buf, attrname.c_str(), attrval.c_str());
            }
        }
        else
        {
            fprintf(stderr, "Ignoring invalid attribute %s\n", attrname.c_str());
        }
        attr = attr->next;
    }
    /* Done processing attributes for this element. Let's see if we still need
     *  to interpret them.
     */
#define xbufString(x) xbuf->find(x)->second
    // Some local variables that we need
    static xmlEntry        parentEntry  = xmlEntry::GENTOP;
    static InteractionType currentItype = InteractionType::BONDS;
    static Identifier      myIdentifier;
    switch (elem)
    {
    case xmlEntry::VERSION:
        if (NN(xbuf, xmlEntry::CHECKSUM))
        {
            pd->setCheckSum(xbufString(xmlEntry::CHECKSUM));
        }
        if (NN(xbuf, xmlEntry::TIMESTAMP))
        {
            pd->setTimeStamp(xbufString(xmlEntry::TIMESTAMP));
        }
        break;
    case xmlEntry::PARTICLETYPES:
        if (NNobligatory(xbuf, xmlEntry::NEXCL))
        {
            pd->setNexcl(atoi(xbufString(xmlEntry::NEXCL).c_str()));
        }
        if (NN(xbuf, xmlEntry::EPSILONR))
        {
            pd->setEpsilonR(atof(xbufString(xmlEntry::EPSILONR).c_str()));
        }
        break;
    case xmlEntry::INTERACTION:
        if (NNobligatory(xbuf, xmlEntry::TYPE) &&
            NNobligatory(xbuf, xmlEntry::FUNCTION) && 
            NNobligatory(xbuf, xmlEntry::CANSWAP))
        {
            CanSwap canSwap = stringToCanSwap(xbufString(xmlEntry::CANSWAP));
            ForceFieldParameterList newparam(xbufString(xmlEntry::FUNCTION), canSwap);
            pd->addForces(xbufString(xmlEntry::TYPE), newparam);
            currentItype = stringToInteractionType(xbufString(xmlEntry::TYPE).c_str());
            parentEntry = elem;
        }
        break;
    case xmlEntry::OPTION:
        if (NNobligatory(xbuf, xmlEntry::KEY) &&
            NNobligatory(xbuf, xmlEntry::VALUE))
        {
            if (xmlEntry::INTERACTION == parentEntry)
            {
                auto fpl = pd->findForces(currentItype);
                fpl->addOption(xbufString(xmlEntry::KEY), xbufString(xmlEntry::VALUE));
            }
            else if (xmlEntry::PARTICLETYPE == parentEntry)
            {
                pd->findParticleType(myIdentifier)->setOption(xbufString(xmlEntry::KEY), xbufString(xmlEntry::VALUE));
            }
        }
        break;
    case xmlEntry::PARAMETERLIST:
        if (NNobligatory(xbuf, xmlEntry::IDENTIFIER))
        {
            myIdentifier = Identifier(currentItype,
                                      xbufString(xmlEntry::IDENTIFIER),
                                      pd->findForces(currentItype)->canSwap());
            parentEntry = elem;
        }
        break;
    case xmlEntry::PARAMETER:
        if (NNobligatory(xbuf, xmlEntry::TYPE)    && NNobligatory(xbuf, xmlEntry::UNIT) &&
            NNobligatory(xbuf, xmlEntry::VALUE)   && NNobligatory(xbuf, xmlEntry::UNCERTAINTY) &&
            NNobligatory(xbuf, xmlEntry::MINIMUM) && NNobligatory(xbuf, xmlEntry::MAXIMUM) &&
            NNobligatory(xbuf, xmlEntry::MUTABILITY))
        {
            Mutability mut;
            if (nameToMutability(xbufString(xmlEntry::MUTABILITY), &mut))
            {
                uint64_t ntrain = 0;
                if (NN(xbuf, xmlEntry::NTRAIN))
                {
                    ntrain = atoi(xbufString(xmlEntry::NTRAIN).c_str());
                }
                bool nonNegative = false;
                if (NN(xbuf, xmlEntry::NONNEGATIVE))
                {
                    nonNegative = stringToBoolean(xbufString(xmlEntry::NONNEGATIVE));
                }
                ForceFieldParameter ffp(xbufString(xmlEntry::UNIT),
                                        xbuf_atof(xbuf, xmlEntry::VALUE),
                                        xbuf_atof(xbuf, xmlEntry::UNCERTAINTY),
                                        ntrain,
                                        xbuf_atof(xbuf, xmlEntry::MINIMUM), xbuf_atof(xbuf, xmlEntry::MAXIMUM),
                                        mut, true,
                                        nonNegative,
                                        true);
                if (xmlEntry::PARAMETERLIST == parentEntry)
                {
                    auto fs = pd->findForces(currentItype);
                    fs->addParameter(myIdentifier, xbufString(xmlEntry::TYPE), std::move(ffp));
                }
                else if (xmlEntry::PARTICLETYPE == parentEntry)
                {
                    pd->findParticleType(myIdentifier)->addForceFieldParameter(xbufString(xmlEntry::TYPE), std::move(ffp));
                }
                else
                {
                    GMX_THROW(gmx::InternalError("Don't know what to do with this parameter"));
                }
            }
            else
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Mutability value %s unknown",
                                                                   xbufString(xmlEntry::MUTABILITY).c_str()).c_str()));
            }
        }
        break;
    case xmlEntry::GT_VSITES:
        if (NN(xbuf, xmlEntry::ANGLE_UNIT) && NN(xbuf, xmlEntry::LENGTH_UNIT))
        {
            pd->setVsite_angle_unit(xbufString(xmlEntry::ANGLE_UNIT));
            pd->setVsite_length_unit(xbufString(xmlEntry::LENGTH_UNIT));
        }
        break;
    case xmlEntry::PARTICLETYPE:
        if (NN(xbuf, xmlEntry::TYPE) &&
            NN(xbuf, xmlEntry::IDENTIFIER) &&
            NN(xbuf, xmlEntry::DESC))
        {
            int ept;
            for(ept = 0; ept < eptNR; ept++)
            {
                if (xbufString(xmlEntry::TYPE).compare(ptype_str[ept]) == 0)
                {
                    break;
                }
            }
            if (ept == eptNR)
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such particle type %s", xbufString(xmlEntry::TYPE).c_str()).c_str()));
            }
            myIdentifier = Identifier(xbufString(xmlEntry::IDENTIFIER));
            pd->addParticleType(ParticleType(myIdentifier,
                                             xbufString(xmlEntry::DESC), ept));
            parentEntry = elem;
        }
        break;
    case xmlEntry::GT_VSITE:
        if (NN(xbuf, xmlEntry::ATYPE)  && NN(xbuf, xmlEntry::VTYPE)    &&
            NN(xbuf, xmlEntry::NUMBER) && NN(xbuf, xmlEntry::DISTANCE) &&
            NN(xbuf, xmlEntry::ANGLE)  && NN(xbuf, xmlEntry::NCONTROLATOMS))
            {
                pd->addVsite(xbufString(xmlEntry::ATYPE),
                             xbufString(xmlEntry::VTYPE),
                             atoi(xbufString(xmlEntry::NUMBER).c_str()),
                             atof(xbufString(xmlEntry::DISTANCE).c_str()),
                             atof(xbufString(xmlEntry::ANGLE).c_str()),
                             atoi(xbufString(xmlEntry::NCONTROLATOMS).c_str()));
            }
        break;
    case xmlEntry::SYM_CHARGE:
        if (NN(xbuf, xmlEntry::CENTRAL) && NN(xbuf, xmlEntry::ATTACHED) &&
            NN(xbuf, xmlEntry::NUMATTACH))
        {
            pd->addSymcharges(xbufString(xmlEntry::CENTRAL),
                              xbufString(xmlEntry::ATTACHED),
                              atoi(xbufString(xmlEntry::NUMATTACH).c_str()));
        }
        break;
    default:
        if (nullptr != debug)
        {
            fprintf(debug, "Unknown combination of attributes:\n");
            for (const auto &i : xml_pd)
            {
                xmlEntry ix = i.second;
                if (xbuf->find(ix) != xbuf->end() &&
                    xbuf->find(ix)->second.size() != 0)
                {
                    fprintf(debug, "%s = %s\n", exml_names(ix), xbuf->find(ix)->second.c_str());
                }
            }
        }
    }
#undef xbufString
}

static void processTree(FILE          *fp, 
                        xmlBuffer     *xbuf,
                        xmlNodePtr     tree,
                        int            indent,
                        Poldata       *pd)
{
    char             buf[100];

    while (tree != nullptr)
    {
        if (fp)
        {
            if ((tree->type > 0) && (tree->type < NXMLTYPES))
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
            if (fp)
            {
                sp(indent, buf, 99);
                fprintf(fp, "%sElement node name %s\n", buf, (char *)tree->name);
            }
            auto iter = xml_pd.find((const char *)tree->name);
            if (iter != xml_pd.end())
            {
                auto elem = iter->second;
                if (elem != xmlEntry::GENTOP)
                {
                    processAttr(fp, tree->properties, xbuf, elem, indent+2, pd);
                }

                if (tree->children)
                {
                    processTree(fp, xbuf, tree->children, indent+2, pd);
                }
            }
        }
        tree = tree->next;
    }
}

void readPoldata(const std::string &fileName,
                 Poldata           *pd)
{
    xmlDocPtr   doc;
    std::string fn, fn2;

    rmap_pd.clear();
    fn = fileName;
    if (fn.empty())
    {
        fn.assign("ACM-g_2020.dat");
    }
    fn2 = gmx::findLibraryFile(fn, true, false);
    if (fn2.empty())
    {
        fn  = "alexandria.ff/" + fn;
        fn2 = gmx::findLibraryFile(fn, true, false);
        if (fn2.empty())
        {
            fn  = "top/" + fn;
            fn2 = gmx::findLibraryFile(fn, true, false);
        }
    }
    if (fn2.empty())
    {
        gmx_fatal(FARGS, "Could not find %s\n", fn.c_str());
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Opening library file %s\n", fn2.c_str());
    }
    xmlDoValidityCheckingDefaultValue = 0;
    doc = xmlParseFile(fn2.c_str());
    if (doc == nullptr)
    {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "Error reading XML file '%s'. Run a syntax checker such as nsgmls.",
                 fn2.c_str());
        GMX_THROW(gmx::FileIOError(buf));
    }

    pd->setFilename(fn2);
    xmlBuffer xbuf;
    processTree(debug, &xbuf, doc->children, 0, pd);

    xmlFreeDoc(doc);

    // Generate maps
    pd->checkForPolarizability();
    pd->checkConsistency(debug);
    if (nullptr != debug)
    {
        writePoldata("pdout.dat", pd, false);
    }
}

static void addOption(xmlNodePtr         parent,
                      const std::string &key,
                      const std::string &value)
{
    auto baby = add_xml_child(parent, exml_names(xmlEntry::OPTION));
    add_xml_char(baby, exml_names(xmlEntry::KEY), key.c_str());
    add_xml_char(baby, exml_names(xmlEntry::VALUE), value.c_str());
}

static void addParameter(xmlNodePtr parent, const std::string &type,
                         const ForceFieldParameter &param)
{
    auto baby = add_xml_child(parent, exml_names(xmlEntry::PARAMETER));
    add_xml_char(baby, exml_names(xmlEntry::TYPE), type.c_str());
    add_xml_char(baby, exml_names(xmlEntry::UNIT), param.unit().c_str());
    add_xml_double(baby, exml_names(xmlEntry::VALUE), param.value());
    add_xml_double(baby, exml_names(xmlEntry::UNCERTAINTY), param.uncertainty());
    add_xml_double(baby, exml_names(xmlEntry::MINIMUM), param.minimum());
    add_xml_double(baby, exml_names(xmlEntry::MAXIMUM), param.maximum());
    add_xml_int(baby, exml_names(xmlEntry::NTRAIN), param.ntrain());
    add_xml_char(baby, exml_names(xmlEntry::MUTABILITY), mutabilityName(param.mutability()).c_str());
    add_xml_char(baby, exml_names(xmlEntry::NONNEGATIVE),
                 param.nonNegative() ? "yes" : "no");
}

static void addXmlPoldata(xmlNodePtr parent, const Poldata *pd)
{
    std::string  geometry, name,
        acentral, attached, tau_unit, ahp_unit,
        epref, desc, params, tmp;

    auto child = add_xml_child(parent, exml_names(xmlEntry::VERSION));
    tmp   = pd->checkSum();
    if (!tmp.empty())
    {
        add_xml_char(child, exml_names(xmlEntry::CHECKSUM), tmp.c_str());
    }
    tmp   = pd->timeStamp();
    if (!tmp.empty())
    {
        add_xml_char(child, exml_names(xmlEntry::TIMESTAMP), tmp.c_str());
    }
    
    child = add_xml_child(parent, exml_names(xmlEntry::PARTICLETYPES));
                 
    add_xml_int(child, exml_names(xmlEntry::NEXCL), pd->getNexcl());
    double epsilonr = pd->getEpsilonR();
    add_xml_double(child, exml_names(xmlEntry::EPSILONR), epsilonr);

    for (const auto &aType : pd->particleTypesConst())
    {
        auto grandchild = add_xml_child(child, exml_names(xmlEntry::PARTICLETYPE));
        add_xml_char(grandchild, exml_names(xmlEntry::IDENTIFIER), aType.id().id().c_str());
        add_xml_char(grandchild, exml_names(xmlEntry::TYPE), ptype_str[aType.gmxParticleType()]);
        add_xml_char(grandchild, exml_names(xmlEntry::DESC), aType.description().c_str());
        for(const auto &opt: aType.optionsConst())
        {
            addOption(grandchild, opt.first, opt.second);
        }
        for(const auto &param : aType.parametersConst())
        {
            addParameter(grandchild, param.first, param.second);
        }
    }
    tmp   = pd->getVsite_angle_unit();
    if (0 != tmp.size())
    {
        child = add_xml_child(parent, exml_names(xmlEntry::GT_VSITES));
        add_xml_char(child, exml_names(xmlEntry::ANGLE_UNIT), tmp.c_str());
    }
    tmp   = pd->getVsite_length_unit();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(xmlEntry::LENGTH_UNIT), tmp.c_str());
    }
    for (auto vsite = pd->getVsiteBegin(); vsite != pd->getVsiteEnd(); vsite++)
    {
        auto grandchild = add_xml_child(child, exml_names(xmlEntry::GT_VSITE));
        add_xml_char(grandchild, exml_names(xmlEntry::ATYPE), vsite->atype().c_str());
        add_xml_char(grandchild, exml_names(xmlEntry::VTYPE), vsiteType2string(vsite->type()));
        add_xml_int(grandchild, exml_names(xmlEntry::NUMBER), vsite->nvsite());
        add_xml_double(grandchild, exml_names(xmlEntry::DISTANCE), vsite->distance());
        add_xml_double(grandchild, exml_names(xmlEntry::ANGLE), vsite->angle());
        add_xml_int(grandchild, exml_names(xmlEntry::NCONTROLATOMS), vsite->ncontrolatoms());
    }
    for (auto &fs : pd->forcesConst())
    {
        auto child = add_xml_child(parent,  exml_names(xmlEntry::INTERACTION));
        add_xml_char(child, exml_names(xmlEntry::TYPE), interactionTypeToString(fs.first).c_str());
        add_xml_char(child, exml_names(xmlEntry::FUNCTION), fs.second.function().c_str());
        add_xml_char(child, exml_names(xmlEntry::CANSWAP), 
                     canSwapToString(fs.second.canSwap()).c_str());
        for (auto &option : fs.second.option())
        {
            addOption(child, option.first, option.second);
        }
        for (auto &params : fs.second.parametersConst())
        {
            auto grandchild = add_xml_child(child, exml_names(xmlEntry::PARAMETERLIST));
            add_xml_char(grandchild, exml_names(xmlEntry::IDENTIFIER), params.first.id().c_str());
            for (const auto &param : params.second)
            {
                addParameter(grandchild, param.first, param.second);
            }
        }
    }

    child = add_xml_child(parent, exml_names(xmlEntry::SYMMETRIC_CHARGES));
    for (auto symcharges = pd->getSymchargesBegin();
         symcharges != pd->getSymchargesEnd(); symcharges++)
    {
        auto grandchild = add_xml_child(child, exml_names(xmlEntry::SYM_CHARGE));
        add_xml_char(grandchild, exml_names(xmlEntry::CENTRAL), symcharges->getCentral().c_str());
        add_xml_char(grandchild, exml_names(xmlEntry::ATTACHED), symcharges->getAttached().c_str());
        add_xml_int(grandchild, exml_names(xmlEntry::NUMATTACH), symcharges->getNumattach());
    }
}

void writePoldata(const std::string &fileName,
                  const Poldata     *pd,
                  bool               compress)
{
    xmlDocPtr   doc;
    xmlDtdPtr   dtd;
    xmlNodePtr  myroot;
    xmlChar    *libdtdname, *dtdname, *gmx;

    rmap_pd.clear();
    gmx        = (xmlChar *) "gentop";
    dtdname    = (xmlChar *) "gentop.dtd";
    libdtdname = dtdname;

    if ((doc = xmlNewDoc((xmlChar *)"1.0")) == nullptr)
    {
        gmx_fatal(FARGS, "Creating XML document %s", fileName.c_str());
    }

    if ((dtd = xmlCreateIntSubset(doc, dtdname, libdtdname, dtdname)) == nullptr)
    {
        gmx_fatal(FARGS, "Creating XML DTD in %s", fileName.c_str());
    }

    if ((myroot = xmlNewDocNode(doc, nullptr, gmx, nullptr)) == nullptr)
    {
        gmx_fatal(FARGS, "Creating root element for %s", fileName.c_str());
    }
    dtd->next    = myroot;
    myroot->prev = (xmlNodePtr) dtd;

    /* Add molecule definitions */
    addXmlPoldata(myroot, pd);

    xmlSetDocCompressMode(doc, compress ? 1 : 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fileName.c_str(), doc, "ISO-8859-1", 2) == 0)
    {
        gmx_fatal(FARGS, "Saving file %s", fileName.c_str());
    }
    xmlFreeDoc(doc);
}

}
