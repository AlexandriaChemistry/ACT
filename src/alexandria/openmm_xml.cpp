/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021
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
 * \author Marie-Madeleine Walz <marie-madeleine.walz@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "openmm_xml.h"

#include <cstdlib>
#include <cstring>

#include <map>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"

#include "forcefieldparameter.h"
#include "forcefieldparameterlist.h"
#include "poldata.h"
#include "poldata_low.h"
#include "xml_util.h"

namespace alexandria
{

/*! \brief The different entry types that can be found
 * \TODO Comment each element
 */
enum class xmlEntry {
    GENTOP,
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

std::map<const std::string, xmlEntry> xmlyyy =
{
    { "gentop",                    xmlEntry::GENTOP           },
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
    { "geometry",                  xmlEntry::GEOMETRY         },
    { "numbonds",                  xmlEntry::NUMBONDS         },
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

std::map<xmlEntry, const std::string> rmapyyy = {};

typedef std::map<xmlEntry, std::string> xmlBuffer;

static const char *exml_names(xmlEntry xml)
{
    if (rmapyyy.empty())
    {
        for (auto &iter : xmlyyy)
        {
            rmapyyy.insert({iter.second, iter.first});
        }
    }
    return rmapyyy[xml].c_str();
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

    auto child = add_xml_child(parent, exml_names(xmlEntry::PARTICLETYPES));
    tmp   = pd->getVersion();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(xmlEntry::VERSION), tmp.c_str());
    }
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

void writeOpenMM(const std::string &fileName,
                 const Poldata     *pd,
                 bool               compress)
{
    xmlDocPtr   doc;
    xmlDtdPtr   dtd;
    xmlNodePtr  myroot;
    xmlChar    *libdtdname, *dtdname, *gmx;

    rmapyyy.clear();
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
