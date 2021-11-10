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

enum class xmlEntryOpenMM {
    FORCEFIELD,
    ATOMTYPES,
    TYPE,
    NAME,
    CLASS,
    ELEMENT,
    MASS,
    RESIDUES,
    RESIDUE,
    ATOM_RES,
    TYPE_RES,
    BOND_RES,
    BONDORDER,
    ATOMNAME1_RES,
    ATOMNAME2_RES,
    EXTERNALBOND_RES,
    ATOMNAME_RES,
    VIRTUALSITE,
    SITENAME,
    WEIGHT1,
    WEIGHT2,
    WEIGHT3,
    HARMONICBONDFORCE,
    CLASS1,
    CLASS2,
    CLASS3,
    CLASS4,
    TYPE1,
    TYPE2,
    TYPE3,
    TYPE4,
    LENGTH,
    K,
    HARMONICANGLEFORCE,
    ANGLE_CLASS,
    ANGLE_VALUE,
    CUSTOMBONDFORCE,
    ENERGY,
    GLOBALPARAMETER,
    DEFAULTVALUE,
    PERBONDPARAMETER, 
    CUSTOMANGLEFORCE,
    PERANGLEPARAMETER,
    CUSTOMTORSIONFORCE,
    PERTORSIONPARAMETER,
    PROPER,
    IMPROPER,
    CUSTOMNONBONDEDFORCE,
    PERPARTICLEPARAMETER,
    CUSTOMMANYPARTICLEFORCE,
    DRUDEFORCE,
    PARTICLE
};

std::map<const std::string, xmlEntryOpenMM> xmlyyyOpenMM =
{
    { "ForceField",                xmlEntryOpenMM::FORCEFIELD       },
    { "AtomTypes",                 xmlEntryOpenMM::ATOMTYPES        },
    { "Type",                      xmlEntryOpenMM::TYPE             },
    { "name",                      xmlEntryOpenMM::NAME             },
    { "class",                     xmlEntryOpenMM::CLASS            },
    { "element",                   xmlEntryOpenMM::ELEMENT          },
    { "mass",                      xmlEntryOpenMM::MASS             },
    { "Residues",                  xmlEntryOpenMM::RESIDUES         },
    { "Residue",                   xmlEntryOpenMM::RESIDUE          },
    { "Atom",                      xmlEntryOpenMM::ATOM_RES         },
    { "type",                      xmlEntryOpenMM::TYPE_RES         },
    { "Bond",                      xmlEntryOpenMM::BOND_RES         },
    { "bondorder",                 xmlEntryOpenMM::BONDORDER        },
    { "atomName1",                 xmlEntryOpenMM::ATOMNAME1_RES    },
    { "atomName2",                 xmlEntryOpenMM::ATOMNAME2_RES    },
    { "ExternalBond",              xmlEntryOpenMM::EXTERNALBOND_RES },
    { "atomName",                  xmlEntryOpenMM::ATOMNAME_RES     },
    { "VirtualSite",               xmlEntryOpenMM::VIRTUALSITE      },
    { "siteName",                  xmlEntryOpenMM::SITENAME         },
    { "weight1",                   xmlEntryOpenMM::WEIGHT1            },
    { "weight2",                   xmlEntryOpenMM::WEIGHT2            },
    { "weight3",                   xmlEntryOpenMM::WEIGHT3            },
    { "HarmonicBondForce",         xmlEntryOpenMM::HARMONICBONDFORCE  },
    { "class1",                    xmlEntryOpenMM::CLASS1             },
    { "class2",                    xmlEntryOpenMM::CLASS2             },
    { "class3",                    xmlEntryOpenMM::CLASS3             },
    { "class4",                    xmlEntryOpenMM::CLASS4             },
    { "type1",                     xmlEntryOpenMM::TYPE1             },
    { "type2",                     xmlEntryOpenMM::TYPE2             },
    { "type3",                     xmlEntryOpenMM::TYPE3             },
    { "type4",                     xmlEntryOpenMM::TYPE4             },
    { "length",                    xmlEntryOpenMM::LENGTH             },
    { "k",                         xmlEntryOpenMM::K                  },
    { "HarmonicAngleForce",        xmlEntryOpenMM::HARMONICANGLEFORCE },
    { "Angle",                     xmlEntryOpenMM::ANGLE_CLASS        },
    { "angle",                     xmlEntryOpenMM::ANGLE_VALUE        },
    { "CustomBondForce",           xmlEntryOpenMM::CUSTOMBONDFORCE    },
    { "energy",                    xmlEntryOpenMM::ENERGY             },
    { "GlobalParameter",           xmlEntryOpenMM::GLOBALPARAMETER    },
    { "defaultValue",              xmlEntryOpenMM::DEFAULTVALUE       },
    { "PerBondParameter",          xmlEntryOpenMM::PERBONDPARAMETER   },   
    { "CustomAngleForce",          xmlEntryOpenMM::CUSTOMANGLEFORCE   },
    { "PerAngleParameter",         xmlEntryOpenMM::PERANGLEPARAMETER  },
    { "CustomTorsionForce",        xmlEntryOpenMM::CUSTOMTORSIONFORCE },
    { "PerTorsionParameter",       xmlEntryOpenMM::PERTORSIONPARAMETER },
    { "Proper",                    xmlEntryOpenMM::PROPER             },
    { "Improper",                  xmlEntryOpenMM::IMPROPER           },
    { "CustomNonbondedForce",      xmlEntryOpenMM::CUSTOMNONBONDEDFORCE },
    { "PerParticleParameter",      xmlEntryOpenMM::PERPARTICLEPARAMETER },
    { "CustomManyParticleForce",   xmlEntryOpenMM::CUSTOMMANYPARTICLEFORCE },
    { "DrudeForce",                xmlEntryOpenMM::DRUDEFORCE              },
    { "Particle",                  xmlEntryOpenMM::PARTICLE                }

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

std::map<xmlEntryOpenMM, const std::string> rmapyyyOpenMM = {};

typedef std::map<xmlEntryOpenMM, std::string> xmlBufferOpenMM;

static const char *exml_names(xmlEntryOpenMM xmlOpenMM)
{
    if (rmapyyyOpenMM.empty())
    {
        for (auto &iter : xmlyyyOpenMM)
        {
            rmapyyyOpenMM.insert({iter.second, iter.first});
        }
    }
    return rmapyyyOpenMM[xmlOpenMM].c_str();
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

static void addSpecParameter(xmlNodePtr parent, const std::string &type,
                         const ForceFieldParameter &param, const std::string &specparam)
{
    if (strcmp(type.c_str(), specparam.c_str()) == 0)
    {
    add_xml_double(parent, specparam.c_str(), param.value());
    }
}

static void addSpecOption(xmlNodePtr         parent,
                      const std::string &key,
                      const std::string &value,
                      const std::string &specopt)
{
    if (strcmp(key.c_str(), specopt.c_str()) == 0)
    {    
    add_xml_char(parent, specopt.c_str(), value.c_str());
    }
}

static void addShell(xmlNodePtr         parent,
                      const std::string &key,
                      const std::string &value,
                      const std::string &specopt)
{
    if (strcmp(key.c_str(), specopt.c_str()) == 0)
    {    
    auto baby = add_xml_child(parent, exml_names(xmlEntryOpenMM::ATOM_RES));    
    add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME), value.c_str());
    add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), value.c_str());
    }
}

static void addXmlPoldata(xmlNodePtr parent, const Poldata *pd)
{
    std::string  geometry, name,
        acentral, attached, tau_unit, ahp_unit,
        epref, desc, params, tmp;

    auto child = add_xml_child(parent, exml_names(xmlEntryOpenMM::ATOMTYPES));
    //tmp   = pd->getVersion();
    //if (0 != tmp.size())
    //{
    //    add_xml_char(child, exml_names(xmlEntry::VERSION), tmp.c_str());
    //}
    //add_xml_int(child, exml_names(xmlEntry::NEXCL), pd->getNexcl());
    //double epsilonr = pd->getEpsilonR();
    //add_xml_double(child, exml_names(xmlEntry::EPSILONR), epsilonr);

    for (const auto &aType : pd->particleTypesConst())
    {
        auto grandchild = add_xml_child(child, exml_names(xmlEntryOpenMM::TYPE));
        add_xml_char(grandchild, exml_names(xmlEntryOpenMM::NAME), aType.id().id().c_str());
        add_xml_char(grandchild, exml_names(xmlEntryOpenMM::CLASS), ptype_str[aType.gmxParticleType()]);
        
        //add_xml_char(grandchild, exml_names(xmlEntryOpenMM::NAME), aType.description().c_str());
        
        for(const auto &opt: aType.optionsConst())
        {
            addSpecOption(grandchild, opt.first, opt.second, "element");
        }
        
        for(const auto &param : aType.parametersConst())
        {
            addSpecParameter(grandchild, param.first, param.second, "mass");
        } 

        //for(const auto &opt: aType.optionsConst())
        //{
        //    addOption(grandchild, opt.first, opt.second);
        //}
        //for(const auto &param : aType.parametersConst())
        //{
        //    addParameter(grandchild, param.first, param.second);
        //}
    }

    auto child2 = add_xml_child(parent, exml_names(xmlEntryOpenMM::RESIDUES));

    for (const auto &aType : pd->particleTypesConst())
    {
        if (strcmp(ptype_str[aType.gmxParticleType()], "Atom") == 0)
        {
            auto grandchild = add_xml_child(child2, exml_names(xmlEntryOpenMM::RESIDUE));
            add_xml_char(grandchild, exml_names(xmlEntryOpenMM::NAME), aType.id().id().c_str());
            auto baby = add_xml_child(grandchild, exml_names(xmlEntryOpenMM::ATOM_RES));
            add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME), aType.id().id().c_str());
            add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), aType.id().id().c_str());
        
            for(const auto &opt: aType.optionsConst())
            {
                addShell(grandchild, opt.first, opt.second, "poltype");
            }
        }

    }


    for (auto &fs : pd->forcesConst())
    {
        if (strcmp(interactionTypeToString(fs.first).c_str(), "BONDS") == 0)
        {
            auto child3 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMBONDFORCE));
            
            int i=0;
            for (auto &params : fs.second.parametersConst())
            {            
                for (const auto &param : params.second)
                {
                    auto grandchild1 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
                    add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), param.first.c_str());
                }
                i++;
                if (i == 1) 
                    {
                    break;
                    }
            }
            auto grandchild1 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME),exml_names(xmlEntryOpenMM::BONDORDER));
            
            for (auto &params : fs.second.parametersConst())
            {
                
                auto grandchild2 = add_xml_child(child3, exml_names(xmlEntryOpenMM::BOND_RES));
                std::string s = params.first.id().c_str();
                std::string delimiter = "#";
                std::string type_1 = s.substr(0, s.find(delimiter));
                s.erase(0, s.find(delimiter) + delimiter.length());
                std::string type_2 = s.substr(0, s.find(delimiter));
                s.erase(0, s.find(delimiter) + delimiter.length());
                std::string bondorder = s;
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE1), type_1.c_str());
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), type_2.c_str());
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::BONDORDER), s.c_str());
                for (const auto &param : params.second)
                {
                    add_xml_double(grandchild2, param.first.c_str(), param.second.value());
                }
            }

        }

        if (strcmp(interactionTypeToString(fs.first).c_str(), "ANGLES") == 0)
        {
            auto child4 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMANGLEFORCE));
            
            int i=0;
            for (auto &params : fs.second.parametersConst())
            {            
                for (const auto &param : params.second)
                {
                    auto grandchild1 = add_xml_child(child4, exml_names(xmlEntryOpenMM::PERANGLEPARAMETER));
                    add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), param.first.c_str());
                } 
                i++;
                if (i == 1) 
                    {
                    break;
                    }
            }

            for (auto &params : fs.second.parametersConst())
            {

                auto grandchild2 = add_xml_child(child4, exml_names(xmlEntryOpenMM::ANGLE_CLASS));
                std::string s = params.first.id().c_str();
                std::string delimiter = "#";
                std::string type_1 = s.substr(0, s.find(delimiter));
                s.erase(0, s.find(delimiter) + delimiter.length());
                std::string type_2 = s.substr(0, s.find(delimiter));
                s.erase(0, s.find(delimiter) + delimiter.length());
                std::string type_3 = s;
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE1), type_1.c_str());
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), type_2.c_str());
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE3), type_3.c_str());
                for (const auto &param : params.second)
                {
                    add_xml_double(grandchild2, param.first.c_str(), param.second.value());
                }
            }    
        }
        

        if (strcmp(interactionTypeToString(fs.first).c_str(), "PROPER_DIHEDRALS") == 0)
        {
            auto child6 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMTORSIONFORCE));
            int i=0;
            for (auto &params : fs.second.parametersConst())
            {            
                for (const auto &param : params.second)
                {
                    auto grandchild1 = add_xml_child(child6, exml_names(xmlEntryOpenMM::PERTORSIONPARAMETER));
                    add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), param.first.c_str());
                } 
                i++;
                if (i == 1) 
                    {
                    break;
                    }
            }

            for (auto &params : fs.second.parametersConst())
            {
               
                auto grandchild2 = add_xml_child(child6, exml_names(xmlEntryOpenMM::PROPER));
                std::string s = params.first.id().c_str();
                std::string delimiter = "#";
                std::string type_1 = s.substr(0, s.find(delimiter));
                s.erase(0, s.find(delimiter) + delimiter.length());
                std::string type_2 = s.substr(0, s.find(delimiter));
                s.erase(0, s.find(delimiter) + delimiter.length());
                std::string type_3 = s.substr(0, s.find(delimiter));
                s.erase(0, s.find(delimiter) + delimiter.length());
                std::string type_4 = s;
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE1), type_1.c_str());
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), type_2.c_str());
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE3), type_3.c_str());
                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE4), type_4.c_str());
                for (const auto &param : params.second)
                {
                    add_xml_double(grandchild2, param.first.c_str(), param.second.value());
                }
            }    
        }
        
        if (strcmp(interactionTypeToString(fs.first).c_str(), "VANDERWAALS") == 0)
        {
            auto child5 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMNONBONDEDFORCE));
            int i=0;
            for (auto &params : fs.second.parametersConst())
            {            
                for (const auto &param : params.second)
                {
                    auto grandchild1 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
                    add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), param.first.c_str());
                } 
                i++;
                if (i == 1) 
                    {
                    break;
                    }
            }
            auto grandchild1 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), "charge");
            auto grandchild3 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::NAME), "zeta");

            for (const auto &aType : pd->particleTypesConst())
            {
                    auto grandchild2 = add_xml_child(child5, exml_names(xmlEntryOpenMM::ATOM_RES));
                    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE_RES), aType.id().id().c_str());

 

                    for(const auto &opt: aType.optionsConst())
                    {
                        if (strcmp(opt.first.c_str(), "vdwtype") == 0) // Will there be a vdwtype? If not have to look through interaction type="VANDERWAALS" and search for atype
                        {
                            for (auto &params : fs.second.parametersConst())
                            {
                                for (const auto &param : params.second)
                                {
                                    if (strcmp(opt.second.c_str(), params.first.id().c_str()) == 0)
                                    {    
                                        add_xml_double(grandchild2, param.first.c_str(), param.second.value());
                                    } 
                                }        
                            }
                        }
                        if (strcmp(opt.first.c_str(), "zetatype") == 0)
                        {
                            for (auto &fs : pd->forcesConst())
                            {
                               if (strcmp(interactionTypeToString(fs.first).c_str(), "CHARGEDISTRIBUTION") == 0)
                               {
                                    for (auto &params : fs.second.parametersConst())
                                    {
                                        for (const auto &param : params.second)
                                        {
                                            if (strcmp(opt.second.c_str(), params.first.id().c_str()) == 0)
                                            {
                                                add_xml_double(grandchild2, param.first.c_str(), param.second.value());    
                                            }    
                                        }
                                            
                                    }    
                               } 
                            }
                        }       
                    }

                    for(const auto &param : aType.parametersConst())
                    {
                        addSpecParameter(grandchild2, param.first, param.second, "charge"); 
                    }       
            }
        }      



        if (strcmp(interactionTypeToString(fs.first).c_str(), "POLARIZATION") == 0)
        {
            auto child6 = add_xml_child(parent, exml_names(xmlEntryOpenMM::DRUDEFORCE));
            for (const auto &aType : pd->particleTypesConst())
            {
                if (strcmp(ptype_str[aType.gmxParticleType()], "Atom") == 0)
                {    
                    auto grandchild2 = add_xml_child(child6, exml_names(xmlEntryOpenMM::PARTICLE)); 
                    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE1), aType.id().id().c_str());
                    
                    for(const auto &opt: aType.optionsConst())
                    {
                        if (strcmp(opt.first.c_str(), "poltype") == 0)
                        {    
                            add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), opt.second.c_str());
                        
                            for (auto &params : fs.second.parametersConst())
                            {
                                if (strcmp(params.first.id().c_str(), opt.second.c_str()) == 0)
                                {
                                    for (const auto &param : params.second)
                                    {
                                    add_xml_double(grandchild2, "polarizability", param.second.value());
                                    }     
                                }    
 
                            } 

                            for (const auto &aType : pd->particleTypesConst())
                            {
                                //if (strcmp(ptype_str[aType.gmxParticleType()], "Shell") == 0)
                                //{
                                    if (strcmp(aType.id().id().c_str(), opt.second.c_str()) == 0) 
                                    {
                                        for(const auto &param : aType.parametersConst())
                                        {
                                            addSpecParameter(grandchild2, param.first, param.second, "charge"); 
                                        } 
                                        
                                    }  
                                //}   
                            }
                
                        }
                    }
                    add_xml_double(grandchild2, "thole", 0);
                } 
            }       
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
    gmx        = (xmlChar *) "ForceField";
    dtdname    = (xmlChar *) "ForceField.dtd";
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
