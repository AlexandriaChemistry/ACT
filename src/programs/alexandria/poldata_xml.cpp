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

#include "poldata_xml.h"

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

enum xmlEntry {
    exmlFirst,
    exmlGENTOP = exmlFirst,
    exmlREFERENCE,
    exmlATOMTYPES,
    exmlATOMTYPE, 
    exmlCOMB_RULE,
    exmlNEXCL,
    exmlVERSION,
    exmlEPSILONR,
    exmlPOLTYPES,
    exmlPOLTYPE,
    exmlPTYPE,
    exmlELEM,
    exmlMASS,
    exmlATOMNUMBER,
    exmlNAME,
    exmlDESC,
    exmlATYPE,
    exmlVALENCE,
    exmlBTYPE,
    exmlZTYPE,
    exmlACMTYPE,
    exmlROW,
    exmlCHARGE,
    exmlNEIGHBORS,
    exmlAROMATIC,
    exmlGEOMETRY,
    exmlNUMBONDS,
    exmlVDWPARAMS,
    exmlEREF,
    exmlVANDERWAALS,
    exmlINTERACTION,
    exmlIDENTIFIER,
    exmlCANSWAP,
    exmlTYPE,
    exmlVALUE,
    exmlOPTION,
    exmlPARAMETERLIST,
    exmlPARAMETER,
    exmlUNCERTAINTY,
    exmlMINIMUM,
    exmlMAXIMUM,
    exmlMUTABILITY,
    exmlATOM1,
    exmlATOM2,
    exmlATOM3,
    exmlATOM4,
    exmlSIGMA,
    exmlBONDORDER,
    exmlPARAMS,
    exmlREFVALUE,
    exmlUNIT,
    exmlNTRAIN,
    exmlGT_VSITES,
    exmlGT_VSITE,
    exmlTAU_AHC,
    exmlALPHA_AHP,
    exmlSYMMETRIC_CHARGES,
    exmlSYM_CHARGE,
    exmlCENTRAL,
    exmlATTACHED,
    exmlNUMATTACH,
    exmlMODEL,
    exmlCHARGEGENERATIONALGORITHM,
    exmlCHARGES,
    exmlREF_CHARGES,
    exmlANGLE_UNIT,
    exmlLENGTH_UNIT,
    exmlDISTANCE,
    exmlNCONTROLATOMS,
    exmlNUMBER,
    exmlVTYPE,
    exmlANGLE,
    exmlFUNCTION,
    exmlLast
};

std::map<const std::string, xmlEntry> xmlxxx =
{
    { "gentop",                 exmlGENTOP           },
    { "reference",              exmlREFERENCE        },
    { "atomtypes",              exmlATOMTYPES        },
    { "atomtype",               exmlATOMTYPE         },
    { "chargegenerationalgorithm", exmlCHARGEGENERATIONALGORITHM },
    { "parameterlist",          exmlPARAMETERLIST    },
    { "parameter",              exmlPARAMETER        },
    { "version",                exmlVERSION          },
    { "function",               exmlFUNCTION         },
    { "nexclusions",            exmlNEXCL            },
    { "epsilonr",               exmlEPSILONR         },
    { "poltypes",               exmlPOLTYPES         },
    { "poltype",                exmlPOLTYPE          },
    { "ptype",                  exmlPTYPE            },
    { "elem",                   exmlELEM             },
    { "mass",                   exmlMASS             },
    { "charge",                 exmlCHARGE           },
    { "row",                    exmlROW              },
    { "atomnumber",             exmlATOMNUMBER       },
    { "name",                   exmlNAME             },
    { "description",            exmlDESC             },
    { "atype",                  exmlATYPE            },
    { "btype",                  exmlBTYPE            },
    { "ztype",                  exmlZTYPE            },
    { "acmtype",                exmlACMTYPE          },
    { "neighbors",              exmlNEIGHBORS        },
    { "aromatic",               exmlAROMATIC         },
    { "geometry",               exmlGEOMETRY         },
    { "numbonds",               exmlNUMBONDS         },
    { "vdwparams",              exmlVDWPARAMS        },
    { "ref_enthalpy",           exmlEREF             },
    { "vanderwaals",            exmlVANDERWAALS      },
    { "interaction",            exmlINTERACTION      },
    { "identifier",             exmlIDENTIFIER       },
    { "canswap",                exmlCANSWAP          },
    { "type",                   exmlTYPE             },
    { "value",                  exmlVALUE            },
    { "option",                 exmlOPTION           },
    { "uncertainty",            exmlUNCERTAINTY      },
    { "minimum",                exmlMINIMUM          },
    { "maximum",                exmlMAXIMUM          },
    { "mutability",             exmlMUTABILITY       },
    { "sigma",                  exmlSIGMA            },
    { "bondorder",              exmlBONDORDER        },
    { "params",                 exmlPARAMS           },
    { "refValue",               exmlREFVALUE         },
    { "unit",                   exmlUNIT             },
    { "ntrain",                 exmlNTRAIN           },
    { "gt_vsites",              exmlGT_VSITES        },
    { "gt_vsite",               exmlGT_VSITE         },
    { "tau_ahc",                exmlTAU_AHC          },
    { "alpha_ahp",              exmlALPHA_AHP        },
    { "symmetric_charges",      exmlSYMMETRIC_CHARGES},
    { "sym_charge",             exmlSYM_CHARGE       },
    { "central",                exmlCENTRAL          },
    { "attached",               exmlATTACHED         },
    { "numattach",              exmlNUMATTACH        },
    { "model",                  exmlMODEL            },
    { "angle_unit",             exmlANGLE_UNIT       },
    { "length_unit",            exmlLENGTH_UNIT      },
    { "distance",               exmlDISTANCE         },
    { "ncontrolatoms",          exmlNCONTROLATOMS    },
    { "number",                 exmlNUMBER           },
    { "vtype",                  exmlVTYPE            },
    { "angle",                  exmlANGLE            }
};

std::map<xmlEntry, const std::string> rmap = {};

typedef std::map<xmlEntry, std::string> xmlBuffer;

static const char *exml_names(xmlEntry xml)
{
    if (rmap.empty())
    {
        for (auto &iter : xmlxxx)
        {
            rmap.insert({iter.second, iter.first});
        }
    }
    return rmap[xml].c_str();
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
                        int         elem,
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

        auto iter = xmlxxx.find(attrname);
        if (iter != xmlxxx.end())
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
    static InteractionType currentItype = InteractionType::BONDS;
    static Identifier      myIdentifier;
    switch (elem)
    {
    case exmlATOMTYPES:
        if (NNobligatory(xbuf, exmlCHARGEGENERATIONALGORITHM))
        {
            pd->setChargeGenerationAlgorithm(xbufString(exmlCHARGEGENERATIONALGORITHM));
        }
        if (NNobligatory(xbuf, exmlVERSION))
        {
            pd->setVersion(xbufString(exmlVERSION));
        }
        if (NNobligatory(xbuf, exmlNEXCL))
        {
            pd->setNexcl(atoi(xbufString(exmlNEXCL).c_str()));
        }
        if (NN(xbuf, exmlEPSILONR))
        {
            pd->setEpsilonR(atof(xbufString(exmlEPSILONR).c_str()));
        }
        break;
    case exmlINTERACTION:
        if (NNobligatory(xbuf, exmlTYPE) &&
            NNobligatory(xbuf, exmlFUNCTION) && 
            NNobligatory(xbuf, exmlCANSWAP))
        {
            CanSwap canSwap = stringToCanSwap(xbufString(exmlCANSWAP));
            ForceFieldParameterList newparam(xbufString(exmlFUNCTION), canSwap);
            pd->addForces(xbufString(exmlTYPE), newparam);
            currentItype = stringToInteractionType(xbufString(exmlTYPE).c_str());
        }
        break;
    case exmlOPTION:
        if (NNobligatory(xbuf, exmlTYPE) &&
            NNobligatory(xbuf, exmlVALUE))
        {
            auto fpl = pd->findForces(currentItype);
            fpl->addOption(xbufString(exmlTYPE), xbufString(exmlVALUE));
        }
        break;
    case exmlPARAMETERLIST:
        if (NNobligatory(xbuf, exmlIDENTIFIER))
        {
            myIdentifier = Identifier(currentItype,
                                      xbufString(exmlIDENTIFIER),
                                      pd->findForces(currentItype)->canSwap());
        }
        break;
    case exmlPARAMETER:
        if (NNobligatory(xbuf, exmlTYPE)    && NNobligatory(xbuf, exmlUNIT) &&
            NNobligatory(xbuf, exmlVALUE)   && NNobligatory(xbuf, exmlUNCERTAINTY) &&
            NNobligatory(xbuf, exmlMINIMUM) && NNobligatory(xbuf, exmlMAXIMUM) &&
            NNobligatory(xbuf, exmlMUTABILITY))
        {
            Mutability mut;
            if (nameToMutability(xbufString(exmlMUTABILITY), &mut))
            {
                uint64_t ntrain = 0;
                if (NN(xbuf, exmlNTRAIN))
                {
                    ntrain = atoi(xbufString(exmlNTRAIN).c_str());
                }
                ForceFieldParameter ffp(xbufString(exmlUNIT),
                                        xbuf_atof(xbuf, exmlVALUE),
                                        xbuf_atof(xbuf, exmlUNCERTAINTY),
                                        ntrain,
                                        xbuf_atof(xbuf, exmlMINIMUM), xbuf_atof(xbuf, exmlMAXIMUM),
                                        mut, false);
                auto fs = pd->findForces(currentItype);
                fs->addParameter(myIdentifier, xbufString(exmlTYPE), std::move(ffp));
            }
            else
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Mutability value %s unknown",
                                                                   xbufString(exmlMUTABILITY).c_str()).c_str()));
            }
        }
        break;
    case exmlGT_VSITES:
        if (NN(xbuf, exmlANGLE_UNIT) && NN(xbuf, exmlLENGTH_UNIT))
        {
            pd->setVsite_angle_unit(xbufString(exmlANGLE_UNIT));
            pd->setVsite_length_unit(xbufString(exmlLENGTH_UNIT));
        }
        break;
    case exmlATOMTYPE:
        Mutability  mut;
        if (NN(xbuf, exmlELEM) &&
            NN(xbuf, exmlMASS) &&
            NN(xbuf, exmlATOMNUMBER) &&
            NN(xbuf, exmlATYPE) &&
            NN(xbuf, exmlBTYPE) &&
            NN(xbuf, exmlPTYPE) &&
            NN(xbuf, exmlZTYPE) &&
            NN(xbuf, exmlACMTYPE) &&
            NN(xbuf, exmlROW) &&
            NN(xbuf, exmlCHARGE) &&
            NN(xbuf, exmlMUTABILITY) &&
            NN(xbuf, exmlEREF) &&
            nameToMutability(xbufString(exmlMUTABILITY), &mut))
        {
            Ffatype sp(xbufString(exmlDESC),  xbufString(exmlATYPE), xbufString(exmlPTYPE),
                       xbufString(exmlBTYPE), xbufString(exmlZTYPE), 
                       xbufString(exmlACMTYPE), xbufString(exmlELEM),
                       atof(xbufString(exmlMASS).c_str()),
                       atoi(xbufString(exmlATOMNUMBER).c_str()),
                       atof(xbufString(exmlCHARGE).c_str()),
                       atoi(xbufString(exmlROW).c_str()),
                       mut,
                       xbufString(exmlEREF));
                pd->addAtype(std::move(sp));
        }
        break;
    case exmlGT_VSITE:
        if (NN(xbuf, exmlATYPE)  && NN(xbuf, exmlVTYPE)    &&
            NN(xbuf, exmlNUMBER) && NN(xbuf, exmlDISTANCE) &&
            NN(xbuf, exmlANGLE)  && NN(xbuf, exmlNCONTROLATOMS))
            {
                pd->addVsite(xbufString(exmlATYPE),
                             xbufString(exmlVTYPE),
                             atoi(xbufString(exmlNUMBER).c_str()),
                             atof(xbufString(exmlDISTANCE).c_str()),
                             atof(xbufString(exmlANGLE).c_str()),
                             atoi(xbufString(exmlNCONTROLATOMS).c_str()));
            }
        break;
    case exmlSYM_CHARGE:
        if (NN(xbuf, exmlCENTRAL) && NN(xbuf, exmlATTACHED) &&
            NN(xbuf, exmlNUMATTACH))
        {
            pd->addSymcharges(xbufString(exmlCENTRAL),
                              xbufString(exmlATTACHED),
                              atoi(xbufString(exmlNUMATTACH).c_str()));
        }
        break;
    default:
        if (nullptr != debug)
        {
            fprintf(debug, "Unknown combination of attributes:\n");
            for (int i = exmlFirst; i < exmlLast; i++)
            {
                xmlEntry ix = static_cast<xmlEntry>(i);
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
            auto iter = xmlxxx.find((const char *)tree->name);
            if (iter != xmlxxx.end())
            {
                int elem = iter->second;
                if (elem != exmlGENTOP)
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

    rmap.clear();
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
    printf("Reading library file %s\n", fn2.c_str());
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

static void addXmlPoldata(xmlNodePtr parent, const Poldata *pd)
{
    xmlNodePtr   child, grandchild;
    int          nexcl;
    std::string  geometry, name,
                 acentral, attached, tau_unit, ahp_unit,
                 epref, desc, params;
    std::string  neighbors, zeta, qstr, rowstr;
    std::string  tmp, func, blu;

    child = add_xml_child(parent, exml_names(exmlATOMTYPES));
    add_xml_char(child, exml_names(exmlCHARGEGENERATIONALGORITHM),
                 chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str());
    tmp   = pd->getVersion();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(exmlVERSION), tmp.c_str());
    }
    nexcl = pd->getNexcl();
    add_xml_int(child, exml_names(exmlNEXCL), nexcl);
    double epsilonr = pd->getEpsilonR();
    add_xml_double(child, exml_names(exmlEPSILONR), epsilonr);

    for (auto aType = pd->getAtypeBegin();
         aType != pd->getAtypeEnd(); aType++)
    {
        grandchild = add_xml_child(child, exml_names(exmlATOMTYPE));
        add_xml_char(grandchild, exml_names(exmlELEM), aType->getElem().c_str());
        add_xml_int(grandchild, exml_names(exmlATOMNUMBER), aType->atomnumber());
        add_xml_double(grandchild, exml_names(exmlMASS), aType->mass());
        add_xml_double(grandchild, exml_names(exmlCHARGE), aType->charge());
        add_xml_int(grandchild, exml_names(exmlROW), aType->row());
        add_xml_char(grandchild, exml_names(exmlDESC), aType->getDesc().c_str());
        add_xml_char(grandchild, exml_names(exmlATYPE), aType->getType().c_str());
        add_xml_char(grandchild, exml_names(exmlPTYPE), aType->id(InteractionType::POLARIZATION).id().c_str());
        add_xml_char(grandchild, exml_names(exmlBTYPE), aType->id(InteractionType::BONDS).id().c_str());
        add_xml_char(grandchild, exml_names(exmlACMTYPE), aType->id(InteractionType::ELECTRONEGATIVITYEQUALIZATION).id().c_str());
        add_xml_char(grandchild, exml_names(exmlZTYPE), aType->id(InteractionType::CHARGEDISTRIBUTION).id().c_str());
        add_xml_char(grandchild, exml_names(exmlMUTABILITY), 
                     mutabilityName(aType->mutability()).c_str());
        add_xml_char(grandchild, exml_names(exmlEREF), aType->getRefEnthalpy().c_str());
    }
    tmp   = pd->getVsite_angle_unit();
    if (0 != tmp.size())
    {
        child = add_xml_child(parent, exml_names(exmlGT_VSITES));
        add_xml_char(child, exml_names(exmlANGLE_UNIT), tmp.c_str());
    }
    tmp   = pd->getVsite_length_unit();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(exmlLENGTH_UNIT), tmp.c_str());
    }
    for (auto vsite = pd->getVsiteBegin(); vsite != pd->getVsiteEnd(); vsite++)
    {
        grandchild = add_xml_child(child, exml_names(exmlGT_VSITE));
        add_xml_char(grandchild, exml_names(exmlATYPE), vsite->atype().c_str());
        add_xml_char(grandchild, exml_names(exmlVTYPE), vsiteType2string(vsite->type()));
        add_xml_int(grandchild, exml_names(exmlNUMBER), vsite->nvsite());
        add_xml_double(grandchild, exml_names(exmlDISTANCE), vsite->distance());
        add_xml_double(grandchild, exml_names(exmlANGLE), vsite->angle());
        add_xml_int(grandchild, exml_names(exmlNCONTROLATOMS), vsite->ncontrolatoms());
    }
    for (auto &fs : pd->forcesConst())
    {
        child = add_xml_child(parent,  exml_names(exmlINTERACTION));
        add_xml_char(child, exml_names(exmlTYPE), interactionTypeToString(fs.first).c_str());
        add_xml_char(child, exml_names(exmlFUNCTION), fs.second.function().c_str());
        add_xml_char(child, exml_names(exmlCANSWAP), 
                     canSwapToString(fs.second.canSwap()).c_str());
        for (auto &option : fs.second.option())
        {
            auto grandChild = add_xml_child(child, exml_names(exmlOPTION));
            add_xml_char(grandChild, exml_names(exmlTYPE), option.first.c_str());
            add_xml_char(grandChild, exml_names(exmlVALUE), option.second.c_str());
        }
        for (auto &params : fs.second.parametersConst())
        {
            auto grandChild = add_xml_child(child, exml_names(exmlPARAMETERLIST));
            add_xml_char(grandChild, exml_names(exmlIDENTIFIER), params.first.id().c_str());
            for (auto &param : params.second)
            {
                auto baby = add_xml_child(grandChild, exml_names(exmlPARAMETER));
                add_xml_char(baby, exml_names(exmlTYPE), param.first.c_str());
                add_xml_char(baby, exml_names(exmlUNIT), param.second.unit().c_str());
                add_xml_double(baby, exml_names(exmlVALUE), param.second.value());
                add_xml_double(baby, exml_names(exmlMINIMUM), param.second.minimum());
                add_xml_double(baby, exml_names(exmlMAXIMUM), param.second.maximum());
                add_xml_double(baby, exml_names(exmlUNCERTAINTY), param.second.uncertainty());
                add_xml_char(baby, exml_names(exmlMUTABILITY), mutabilityName(param.second.mutability()).c_str());
                add_xml_int(baby, exml_names(exmlNTRAIN), param.second.ntrain());
            }
        }
    }

    child = add_xml_child(parent, exml_names(exmlSYMMETRIC_CHARGES));
    for (auto symcharges = pd->getSymchargesBegin();
         symcharges != pd->getSymchargesEnd(); symcharges++)
    {
        grandchild = add_xml_child(child, exml_names(exmlSYM_CHARGE));
        add_xml_char(grandchild, exml_names(exmlCENTRAL), symcharges->getCentral().c_str());
        add_xml_char(grandchild, exml_names(exmlATTACHED), symcharges->getAttached().c_str());
        add_xml_int(grandchild, exml_names(exmlNUMATTACH), symcharges->getNumattach());
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

    rmap.clear();
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
