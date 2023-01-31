/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2022
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
 * \author Marie-Madeleine Walz <marie-madeleine.walz@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "openmm_xml.h"

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <map>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/topology/topology.h"


#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_low.h"
#include "act/molprop/molprop_util.h"
#include "actmol.h"
#include "act/utility/xml_util.h"


namespace alexandria
{

/*! \brief The different entry types that can be found
 * TODO: Comment each element
 * xmlEntryOpenMM defines the keys that are in the xmlyyyOpenMM map
 * xmlyyyOpenMM is the map listing value and key used for the OpenMM xml file
 * exml_names returns string value upon calling function and submitting the key 
 * addSpecParameter adds specific parameter for an atomtype upon submitting a string value, here also conversion 
 *       of strings and units are done so that it is suitable for OpenMM, e.g. angles have to be in radians
 * addSpecOption adds specific option for an atomtype 
 * addShell adds shell particle for a specific core atomtype
 * addXmlForceField creates the xml tree ForceField, it consists of Atomtypes, Residues, HarmonicBondForce, HarmonicAngleForce, CustomNonBondedForce, NonBondedForce and DrudeForce
 * !!! The order of the parameters is in certain forces (e.g. CustomNonBondedForce) very important; reordering might break the force field in OpenMM!!!
 * writeOpenMM
 */

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
    CHARGE_RES,
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
    AMOEBA_UB_FORCE,
    UREY_BRADLEY_FORCE,
    CUSTOMBONDFORCE,
    ENERGY,
    GLOBALPARAMETER,
    DEFAULTVALUE,
    PERBONDPARAMETER, 
    CUSTOMANGLEFORCE,
    PERANGLEPARAMETER,
    RBTORSIONFORCE,
    PERIODICTORSIONFORCE,
    CUSTOMTORSIONFORCE,
    NONBONDEDFORCE,
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
    { "charge",                    xmlEntryOpenMM::CHARGE_RES       },
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
    { "AmoebaUreyBradleyForce",    xmlEntryOpenMM::AMOEBA_UB_FORCE    },
    { "UreyBradley",               xmlEntryOpenMM::UREY_BRADLEY_FORCE },
    { "CustomBondForce",           xmlEntryOpenMM::CUSTOMBONDFORCE    },
    { "energy",                    xmlEntryOpenMM::ENERGY             },
    { "GlobalParameter",           xmlEntryOpenMM::GLOBALPARAMETER    },
    { "defaultValue",              xmlEntryOpenMM::DEFAULTVALUE       },
    { "PerBondParameter",          xmlEntryOpenMM::PERBONDPARAMETER   },   
    { "CustomAngleForce",          xmlEntryOpenMM::CUSTOMANGLEFORCE   },
    { "PerAngleParameter",         xmlEntryOpenMM::PERANGLEPARAMETER  },
    { "RBTorsionForce",            xmlEntryOpenMM::RBTORSIONFORCE     },
    { "PeriodicTorsionForce",      xmlEntryOpenMM::PERIODICTORSIONFORCE },
    { "CustomTorsionForce",        xmlEntryOpenMM::CUSTOMTORSIONFORCE },
    { "PerTorsionParameter",       xmlEntryOpenMM::PERTORSIONPARAMETER },
    { "Proper",                    xmlEntryOpenMM::PROPER             },
    { "Improper",                  xmlEntryOpenMM::IMPROPER           },
    { "CustomNonbondedForce",      xmlEntryOpenMM::CUSTOMNONBONDEDFORCE },
    { "NonbondedForce",            xmlEntryOpenMM::NONBONDEDFORCE },
    { "PerParticleParameter",      xmlEntryOpenMM::PERPARTICLEPARAMETER },
    { "CustomManyParticleForce",   xmlEntryOpenMM::CUSTOMMANYPARTICLEFORCE },
    { "DrudeForce",                xmlEntryOpenMM::DRUDEFORCE              },
    { "Particle",                  xmlEntryOpenMM::PARTICLE                } 
};
    
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

static std::string atomTypeOpenMM(const std::string &ffType, size_t index)
{
    return gmx::formatString("%s_%lu", ffType.c_str(), index);
}

std::vector<xmlEntryOpenMM> class_vec = {xmlEntryOpenMM::CLASS1, xmlEntryOpenMM::CLASS2, xmlEntryOpenMM::CLASS3, xmlEntryOpenMM::CLASS4};
std::vector<xmlEntryOpenMM> type_vec  = {xmlEntryOpenMM::TYPE1, xmlEntryOpenMM::TYPE2, xmlEntryOpenMM::TYPE3, xmlEntryOpenMM::TYPE4}; 

static void addSpecParameter(xmlNodePtr                 parent, 
                             const std::string         &type,
                             const ForceFieldParameter &param,
                             const std::string         &specparam)
{
    // TODO remove explicit conversion and use ACT routines.
    double value = param.internalValue();
    std::map<std::string, std::string> gmx2OMM = {
        { morse_name[morseLENGTH], "r0"  },
        { morse_name[morseDE],    "D_e" },
        { morse_name[morseBETA],  "a"   },
        { angle_name[angleKT],    "k"   },
        { angle_name[angleANGLE], "angle" },
        { ub_name[ubKUB],         "k" },
        { ub_name[ubR13],         "d" }
    };
    if (type == specparam)
    {
        // mass for Langevin, subtract shell mass = 0.1
        if (type == "mass")
        {
            add_xml_double(parent, specparam.c_str(), value-0.1); 
        }
        else if (type == "charge")
        {
            add_xml_double(parent, specparam.c_str(), value); 
        }
        else
        {
            if (gmx2OMM.find(type) != gmx2OMM.end())
            {
                add_xml_double(parent, gmx2OMM[type], value);
            }
            else
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Unknown parameter type %s when converting to OpenMM", type.c_str()).c_str()));
            }
        }
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

static void addXmlElemMass(xmlNodePtr parent, const ParticleType &aType)
{
    double mDrude = 0.1;

    if (eptAtom == aType.gmxParticleType())
    {
        add_xml_char(parent, exml_names(xmlEntryOpenMM::ELEMENT),
                     aType.element().c_str());
        double mAtom = aType.mass();
        if (aType.hasOption("poltype"))
        {
            mAtom -= mDrude;
        }
        add_xml_double(parent, exml_names(xmlEntryOpenMM::MASS), mAtom);
    }
    else if (eptShell == aType.gmxParticleType())
    {
        add_xml_double(parent, exml_names(xmlEntryOpenMM::MASS), mDrude);
    }
}

static void addXmlForceField(xmlNodePtr parent, const ForceField *pd, const ACTMol *actmol)
{
    std::string  geometry, name,
        acentral, attached, tau_unit, ahp_unit,
        epref, desc, params, tmp;

    auto child = add_xml_child(parent, exml_names(xmlEntryOpenMM::ATOMTYPES));

    for (const auto &aType : pd->particleTypesConst())
    {
        auto grandchild = add_xml_child(child, exml_names(xmlEntryOpenMM::TYPE));
        add_xml_char(grandchild, exml_names(xmlEntryOpenMM::NAME), aType.id().id().c_str());
        add_xml_char(grandchild, exml_names(xmlEntryOpenMM::CLASS), aType.id().id().c_str());
        addXmlElemMass(grandchild, aType);
    }

    auto myatoms =  actmol -> atomsConst();
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        auto name_ai = atomTypeOpenMM(myatoms[i].ffType(), i);

        auto baby = add_xml_child(child, exml_names(xmlEntryOpenMM::TYPE));
        add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME), name_ai.c_str()); 
        auto ffType = myatoms[i].ffType();
        add_xml_char(baby, exml_names(xmlEntryOpenMM::CLASS), ffType.c_str());

        if (!pd->hasParticleType(ffType))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("No such particle type %s in force field %s", ffType.c_str(), pd->filename().c_str()).c_str()));
        }
        auto aType = pd->findParticleType(ffType);
        addXmlElemMass(baby, *aType);
    }

    auto child2 = add_xml_child(parent, exml_names(xmlEntryOpenMM::RESIDUES));

    std::vector<std::string> mylist{"Li+", "Na+", "K+", "Rb+", "Cs+", "F-", "Cl-", "Br-", "I-"};

    for (const auto &aType : pd->particleTypesConst())
    {
        if (strcmp(ptype_str[aType.gmxParticleType()], "Atom") == 0)
        {
            if (std::find(std::begin(mylist), std::end(mylist), aType.id().id().c_str()) != std::end(mylist))
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
    }

    if (actmol->getMolname().size() > 0)
    {
        auto grandchild = add_xml_child(child2, exml_names(xmlEntryOpenMM::RESIDUE));
        add_xml_char(grandchild, exml_names(xmlEntryOpenMM::NAME), actmol->getMolname().c_str());
        
        auto myatoms =  actmol -> atomsConst();
        for (size_t i = 0; i < myatoms.size(); i++)
        {
            auto name_ai = atomTypeOpenMM(myatoms[i].ffType(), i);

            auto baby = add_xml_child(grandchild, exml_names(xmlEntryOpenMM::ATOM_RES));
            add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME), name_ai.c_str()); 
            //add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), *(myatoms.atomtype[i]));
            add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), name_ai.c_str());  
            add_xml_double(baby, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[i].charge());
        }    

        auto itbonds = InteractionType::BONDS;
        if (pd->interactionPresent(itbonds))
        {
            auto fs = pd->findForcesConst(itbonds);
        
            if (actmol->topology()->hasEntry(itbonds))
            {
                for(const auto topentry : actmol->topology()->entry(itbonds))
                {
                    int ai = topentry->atomIndex(0);
                    int aj = topentry->atomIndex(1);
                    
                    auto name_ai = atomTypeOpenMM(myatoms[ai].ffType(), ai);
                    auto name_aj = atomTypeOpenMM(myatoms[aj].ffType(), aj);
                    
                    auto baby = add_xml_child(grandchild, exml_names(xmlEntryOpenMM::BOND_RES));
                    add_xml_char(baby, exml_names(xmlEntryOpenMM::ATOMNAME1_RES), name_ai.c_str());
                    add_xml_char(baby, exml_names(xmlEntryOpenMM::ATOMNAME2_RES), name_aj.c_str());  
                }   
            }
        }
    } 


    for (auto &fs : pd->forcesConst())
    {
        // This adds the Morse potential;
        // TODO: link the bondorder information to the paramter assignment that happens via atomtypes
        if (strcmp(interactionTypeToString(fs.first).c_str(), "BONDS") == 0)
        {
            auto child3 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMBONDFORCE));
            // The Morse potential could be written as a string here, or it can be added in the openmm python script
            add_xml_double(child3, "energy", 0.0); 
 
            // Specify the per bond parameters
            auto grandchild0 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), "D_e");
            auto grandchild1 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), "a");
            auto grandchild2 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::NAME), "r0");
            
            // Add all bonds
            for (auto &params : fs.second.parametersConst())
            {
                auto grandchild3 = add_xml_child(child3, exml_names(xmlEntryOpenMM::BOND_RES));
                
                int i = 0;
                for (auto &a:params.first.atoms())
                {
                    std::string s = a.c_str();
 
                    if (!s.empty()) 
                    {
                        s.resize(s.size() - 2);
                    }
                    add_xml_char(grandchild3, exml_names(class_vec[i]), s.c_str());
                    i++;
                } 

                for (const auto &param : params.second)
                {
                    addSpecParameter(grandchild3, param.first, param.second, morse_name[morseDE]);
                    addSpecParameter(grandchild3, param.first, param.second, morse_name[morseBETA]);
                    addSpecParameter(grandchild3, param.first, param.second, morse_name[morseLENGTH]);  
                }
            }
        }    
        //    
        //    int i=0;
        //    for (auto &params : fs.second.parametersConst())
        //    {            
        //        for (const auto &param : params.second)
        //        {
        //            auto grandchild1 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
        //            add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), param.first.c_str());
        //        }
        //        i++;
        //        if (i == 1) 
        //            {
        //            break;
        //            }
        //    }
        //    auto grandchild1 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
        //    add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME),exml_names(xmlEntryOpenMM::BONDORDER));

        //    for (auto &params : fs.second.parametersConst())
        //    {
                
        //        auto grandchild2 = add_xml_child(child3, exml_names(xmlEntryOpenMM::BOND_RES));
        //        int i = 0;
        //        for (auto &a:params.first.atoms())
        //        {
        //            add_xml_char(grandchild2, exml_names(class_vec[i]), a.c_str());
        //            i++;
        //        } 

        //        for (auto &b:params.first.bondOrders())
        //        {
        //            add_xml_double(grandchild2, exml_names(xmlEntryOpenMM::BONDORDER), b);
        //        }
        //        
        //        for (const auto &param : params.second)
        //        {
        //            add_xml_double(grandchild2, param.first.c_str(), param.second.value());
        //        }
        //    }

        //}

        // This part is to add a harmonic bond, currently removed as we use Morse potential for bonds
        //if (strcmp(interactionTypeToString(fs.first).c_str(), "BONDS") == 0)
        //{
        //    auto child3 = add_xml_child(parent, exml_names(xmlEntryOpenMM::HARMONICBONDFORCE));  
        //    for (auto &params : fs.second.parametersConst())
        //    {
        //        auto grandchild2 = add_xml_child(child3, exml_names(xmlEntryOpenMM::BOND_RES));
        //        
        //        int i = 0;
        //        for (auto &a:params.first.atoms())
        //        {
        //            std::string s = a.c_str();
        //            if (!s.empty()) 
        //            {
        //                s.resize(s.size() - 2);
        //            }
        //            add_xml_char(grandchild2, exml_names(class_vec[i]), s.c_str());
        //            i++;
        //        } 
        //        for (const auto &param : params.second)
        //        {
        //            addSpecParameter(grandchild2, param.first, param.second, "bondlength");
        //        }
        //        // hardcoded force constant, should be changed if we'd like to use harmonic bond
        //        add_xml_double(grandchild2, "k", 458148 );
        //   }
        //}


        if (strcmp(interactionTypeToString(fs.first).c_str(), "ANGLES") == 0)

        {
            auto child4 = add_xml_child(parent, exml_names(xmlEntryOpenMM::HARMONICANGLEFORCE));

            for (auto &params : fs.second.parametersConst())
            {

                auto grandchild2 = add_xml_child(child4, exml_names(xmlEntryOpenMM::ANGLE_CLASS));
                
                int i=0;
                for (auto &a:params.first.atoms())
                {
                    std::string s = a.c_str();
                    if (!s.empty()) 
                    {
                        s.resize(s.size() - 2);
                    }
                    add_xml_char(grandchild2, exml_names(class_vec[i]), s.c_str());
                    i++;
                } 

                for (const auto &param : params.second)
                {
                    addSpecParameter(grandchild2, param.first, param.second, angle_name[angleANGLE]); 
                    addSpecParameter(grandchild2, param.first, param.second, angle_name[angleKT]); 
                }
            }    
        }

        // Currently removed, as we decided to use harmonic angles
        //if (strcmp(interactionTypeToString(fs.first).c_str(), "ANGLES") == 0)
        //{
        //    auto child4 = add_xml_child(parent, exml_names(xmlEntryOpenMM::AMOEBA_UB_FORCE));
        //
        //    for (auto &params : fs.second.parametersConst())
        //    {
        //        auto grandchild2 = add_xml_child(child4, exml_names(xmlEntryOpenMM::UREY_BRADLEY_FORCE));
        //         
        //        int i = 0;
        //        for (auto &a:params.first.atoms())
        //        {
        //            std::string s = a.c_str();
        //            if (!s.empty()) 
        //            {
        //                s.resize(s.size() - 2);
        //            }
        //            add_xml_char(grandchild2, exml_names(class_vec[i]), s.c_str());
        //            i++;
        //        } 
        //        for (const auto &param : params.second)
        //        {
        //            addSpecParameter(grandchild2, param.first, param.second, "kub"); 
        //            addSpecParameter(grandchild2, param.first, param.second, "r13"); 
        //            //add_xml_double(grandchild2, param.first.c_str(), param.second.value());
        //        }
        //    }    
        //}

        //if (strcmp(interactionTypeToString(fs.first).c_str(), "LINEAR_ANGLES") == 0) 
        //{
        //    auto child4 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMANGLEFORCE));
            
        //    int i=0;
        //    for (auto &params : fs.second.parametersConst())
        //    {            
        //        for (const auto &param : params.second)
        //        {
        //            auto grandchild1 = add_xml_child(child4, exml_names(xmlEntryOpenMM::PERANGLEPARAMETER));
        //            add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), param.first.c_str());
        //        } 
        //        i++;
        //        if (i == 1) 
        //            {
        //            break;
        //            }
        //    }

        //    for (auto &params : fs.second.parametersConst())
        //    {

        //        auto grandchild2 = add_xml_child(child4, exml_names(xmlEntryOpenMM::ANGLE_CLASS));
        //        int i = 0;
        //        for (auto &a:params.first.atoms())
        //        {
                      //std::string s = a.c_str();
                      //if (!s.empty()) 
                      //{
                      //     s.resize(s.size() - 2);
                      //}
        //            add_xml_char(grandchild2, exml_names(class_vec[i]), s.c_str());
        //            i++;
        //        }

        //        for (const auto &param : params.second)
        //        {
        //            add_xml_double(grandchild2, param.first.c_str(), param.second.value());
        //        }
        //    }
        //}
        

        if (strcmp(interactionTypeToString(fs.first).c_str(), "PROPER_DIHEDRALS") == 0)
        {
            auto child6 = add_xml_child(parent, exml_names(xmlEntryOpenMM::RBTORSIONFORCE));

            for (auto &params : fs.second.parametersConst())
            {
               
                auto grandchild2 = add_xml_child(child6, exml_names(xmlEntryOpenMM::PROPER));
                int i=0;
                for (auto &a:params.first.atoms())
                {
                    std::string s = a.c_str();
                    if (!s.empty()) 
                    {
                        s.resize(s.size() - 2);
                    }
                    add_xml_char(grandchild2, exml_names(class_vec[i]), s.c_str());
                    i++;
                } 

                for (const auto &param : params.second)
                {
                    add_xml_double(grandchild2, param.first.c_str(), param.second.value());
                }
                // OpenMM expects 7 parameters for RBTorsion force
                // TODO: Do we also need more parameters for the RBTorsion force in ACT? Currently c0, c1, c2, c3 in ACT.
                add_xml_double(grandchild2, "c4", 0);
                add_xml_double(grandchild2, "c5", 0);
                add_xml_double(grandchild2, "c6", 0);
            }    
        }

        if (strcmp(interactionTypeToString(fs.first).c_str(), "IMPROPER_DIHEDRALS") == 0)
        {
            auto child6 = add_xml_child(parent, exml_names(xmlEntryOpenMM::PERIODICTORSIONFORCE));

            for (auto &params : fs.second.parametersConst())
            {
               
                auto grandchild2 = add_xml_child(child6, exml_names(xmlEntryOpenMM::IMPROPER));
                int i=0;
                for (auto &a:params.first.atoms())
                {
                    std::string s = a.c_str();
                    if (!s.empty()) 
                    {
                        s.resize(s.size() - 2);
                    }
                    add_xml_char(grandchild2, exml_names(class_vec[i]), s.c_str());
                    i++;
                } 

                for (const auto &param : params.second)
                {
                    add_xml_double(grandchild2, param.first.c_str(), param.second.value());
                }
            }    
        }

        
        if (strcmp(interactionTypeToString(fs.first).c_str(), "VANDERWAALS") == 0)
        {
            auto child5 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMNONBONDEDFORCE));
            add_xml_double(child5, "energy", 0.0); 
            add_xml_double(child5, "bondCutoff", 3.0); 

            auto grandchild0 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), "vdW");
            
            //int i=0;
            //for (auto &params : fs.second.parametersConst())
            //{            
            //    for (const auto &param : params.second)
            //    {
            //        auto grandchild1 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            //        add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), param.first.c_str());
            //    } 
            //    i++;
            //    if (i == 1) 
            //        {
            //        break;
            //        }
            //}
 
            // !!!   This order is important !!! Do not change the order of these parameters, otherwise the force field is not working !!!
            auto grandchild2 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::NAME), "sigma");
            auto grandchild3 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::NAME), "epsilon");
            auto grandchild4 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            add_xml_char(grandchild4, exml_names(xmlEntryOpenMM::NAME), "gamma");


            auto grandchild5 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            add_xml_char(grandchild5, exml_names(xmlEntryOpenMM::NAME), "charge");
            auto grandchild6 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            add_xml_char(grandchild6, exml_names(xmlEntryOpenMM::NAME), "beta");



            for (const auto &aType : pd->particleTypesConst())
            {
                    auto grandchild2 = add_xml_child(child5, exml_names(xmlEntryOpenMM::ATOM_RES));
                    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::CLASS), aType.id().id().c_str());
                    if (strcmp( ptype_str[aType.gmxParticleType()], "Atom") == 0) 
                    {
                       add_xml_double(grandchild2, "vdW", 1.0); 
                    }
                    if (strcmp( ptype_str[aType.gmxParticleType()], "Shell") == 0) 
                    {
                        add_xml_double(grandchild2, "vdW", 0.0);
                    }
                    
 
                    for(const auto &opt: aType.optionsConst())
                    {
                        if (strcmp(opt.first.c_str(), "vdwtype") == 0) 
                        {
                            for (auto &params : fs.second.parametersConst())
                            {
                                for (const auto &param : params.second)
                                {
                                    if (strcmp(opt.second.c_str(), params.first.id().c_str()) == 0)
                                    {    
                                        if (strcmp( ptype_str[aType.gmxParticleType()], "Atom") == 0) 
                                        {
                                            add_xml_double(grandchild2, param.first.c_str(), param.second.value());
                                        }
                                        if (strcmp( ptype_str[aType.gmxParticleType()], "Shell") == 0) 
                                        {
                                            if (strcmp(param.first.c_str(), "gamma") == 0)
                                            {
                                               add_xml_double(grandchild2, param.first.c_str(), 7.0); 
                                            }
                                            else 
                                            {
                                                add_xml_double(grandchild2, param.first.c_str(), 1.0); 
                                            }
                                            
                                        }
                                    } 
                                }        
                            }
                        }
                        if (strcmp(opt.first.c_str(), "zetatype") == 0)
                        {
                            for (auto &fs : pd->forcesConst())
                            {
                               if (strcmp(interactionTypeToString(fs.first).c_str(), "COULOMB") == 0)
                               {
                                    for (auto &params : fs.second.parametersConst())
                                    {
                                        for (const auto &param : params.second)
                                        {
                                            if (strcmp(opt.second.c_str(), params.first.id().c_str()) == 0)
                                            {
                                        
                    
                                                for(const auto &param : aType.parametersConst())
                                                {
                                                    addSpecParameter(grandchild2, param.first, param.second, "charge");     
                                                } 
                                                //add_xml_double(grandchild2, param.first.c_str(), param.second.value()); 
                                                add_xml_double(grandchild2, "beta", param.second.value()); 
                                            }    
                                        }
                                            
                                    }    
                               } 
                            }
                        }       
                    }

                    //add_xml_int(grandchild2, "charge", 0);       
            }
            // add customnonbonded (WBH vdW + pg coulomb) for compound
            auto myatoms =  actmol -> atomsConst();
            for (size_t i = 0; i < myatoms.size(); i++)
            {
                auto name_ai = atomTypeOpenMM(myatoms[i].ffType(), i);
                
                auto grandchild3 = add_xml_child(child5, exml_names(xmlEntryOpenMM::ATOM_RES));
                add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::TYPE_RES), name_ai.c_str());  
                
                for (const auto &aType : pd->particleTypesConst())
                {
                    if (aType.id().id() == myatoms[i].ffType())
                    {
                        if (strcmp( ptype_str[aType.gmxParticleType()], "Atom") == 0) 
                        {
                            add_xml_double(grandchild3, "vdW", 1.0); 
                        }
                        if (strcmp( ptype_str[aType.gmxParticleType()], "Shell") == 0)
                        {
                            add_xml_double(grandchild3, "vdW", 0.0);
                        }   
                      
                
                        for(const auto &opt: aType.optionsConst())
                        { 
                            if (strcmp(opt.first.c_str(), "vdwtype") == 0) 
                            {
                                for (auto &params : fs.second.parametersConst())
                                {
                                    for (const auto &param : params.second)
                                    {
                                        if (strcmp(opt.second.c_str(), params.first.id().c_str()) == 0)
                                        {    
                                            if (strcmp( ptype_str[aType.gmxParticleType()], "Atom") == 0) 
                                            {
                                                add_xml_double(grandchild3, param.first.c_str(), param.second.value());
                                            }
                                            if (strcmp( ptype_str[aType.gmxParticleType()], "Shell") == 0) 
                                            {
                                                if (strcmp(param.first.c_str(), "gamma") == 0)
                                                {
                                                   add_xml_double(grandchild3, param.first.c_str(), 7.0); 
                                                }
                                                else 
                                                {
                                                    add_xml_double(grandchild3, param.first.c_str(), 1.0); 
                                                }
                                            
                                            }
                                        } 
                                    }        
                                }
                            }
                            if (strcmp(opt.first.c_str(), "zetatype") == 0)
                            {
                                for (auto &fs : pd->forcesConst())
                                {
                                    if (strcmp(interactionTypeToString(fs.first).c_str(), "COULOMB") == 0)
                                    {
                                        for (auto &params : fs.second.parametersConst())
                                        {
                                            for (const auto &param : params.second)
                                            {
                                                if (strcmp(opt.second.c_str(), params.first.id().c_str()) == 0)
                                                {
                                                    add_xml_double(grandchild3, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[i].charge()); 
                                                    add_xml_double(grandchild3, "beta", param.second.value()); 
                                                }    
                                            }  
                                        }    
                                    } 
                                }
                            }
                        }        
                    } 
                } 
            }        
        }      

        // This part is added to implement PME and LJPME, i.e. long-range interactions are approx. by point charge and 12-6 Lennard-Jones
        // For the different atomtypes it would be good to add optimized sigma and epsilon values 
        // in order to calculate the vdW contribution after the cutoff
        if (strcmp(interactionTypeToString(fs.first).c_str(), "VANDERWAALS") == 0)
        {
            auto child5 = add_xml_child(parent, exml_names(xmlEntryOpenMM::NONBONDEDFORCE));
            add_xml_double(child5, "coulomb14scale", 1.0); 
            add_xml_double(child5, "lj14scale", 1.0); 
            
            for (const auto &aType : pd->particleTypesConst())
            {
                    auto grandchild2 = add_xml_child(child5, exml_names(xmlEntryOpenMM::ATOM_RES));
                    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::CLASS), aType.id().id().c_str());
                    for(const auto &param : aType.parametersConst())
                    {
                        addSpecParameter(grandchild2, param.first, param.second, "charge"); 
                        if (strcmp( ptype_str[aType.gmxParticleType()], "Atom") == 0) 
                        {
                            add_xml_double(grandchild2, "sigma", 0.3);  
                            add_xml_double(grandchild2, "epsilon", 0.05);  
                        }
                        if (strcmp( ptype_str[aType.gmxParticleType()], "Shell") == 0) 
                        {
                            add_xml_double(grandchild2, "sigma", 0.3);  
                            add_xml_double(grandchild2, "epsilon", 0.0);  
                        }

                    }

            } 
            // adding nonbonded (LJ + point coulomb) for compound
            auto myatoms =  actmol -> atomsConst();
            for (size_t i = 0; i < myatoms.size(); i++)
            {
                auto name_ai = atomTypeOpenMM(myatoms[i].ffType(), i);

                auto baby = add_xml_child(child5, exml_names(xmlEntryOpenMM::ATOM_RES));
                add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), name_ai.c_str()); 
                for (const auto &aType : pd->particleTypesConst())
                {
                    if (aType.id().id() == myatoms[i].ffType())
                    {

                        add_xml_double(baby, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[i].charge());
                        if (strcmp( ptype_str[aType.gmxParticleType()], "Atom") == 0) 
                        {
                            add_xml_double(baby, "sigma", 0.3);  
                            add_xml_double(baby, "epsilon", 0.05);  
                        }
                        if (strcmp( ptype_str[aType.gmxParticleType()], "Shell") == 0) 
                        {
                            add_xml_double(baby, "sigma", 0.3);  
                            add_xml_double(baby, "epsilon", 0.0);  
                        }
                           

                    }      
                    
                }        
                
            }       
        }

        // !!! Shell particle has to be type1, core particle has to be type2 !!!
        // 
        if (strcmp(interactionTypeToString(fs.first).c_str(), "POLARIZATION") == 0)
        {  
            if (pd->polarizable())    
            {
                auto child6 = add_xml_child(parent, exml_names(xmlEntryOpenMM::DRUDEFORCE));
                for (const auto &aType : pd->particleTypesConst())
                {
                    if (eptAtom == aType.gmxParticleType())
                    {
                        auto grandchild2 = add_xml_child(child6, exml_names(xmlEntryOpenMM::PARTICLE)); 
                        add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), aType.id().id().c_str());
                        double myalpha = 0;
                        if (aType.hasOption("poltype"))
                        {
                            auto shell_ai = aType.optionValue("poltype");
                            auto alpha = fs.second.findParameterTypeConst(Identifier({shell_ai}),
                                                                          pol_name[polALPHA]);
                            myalpha = alpha.internalValue();
                            add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE1), shell_ai.c_str());
                        }
                        add_xml_double(grandchild2, "polarizability", myalpha);
                        add_xml_double(grandchild2, "thole", 0);
                        const std::string charge("charge");
                        addSpecParameter(grandchild2, exml_names(xmlEntryOpenMM::CHARGE_RES), aType.parameterConst(charge), charge);
                    } 
                } 

                auto myatoms =  actmol -> atomsConst();
                for (size_t i = 0; i < myatoms.size(); i++)
                {
                    auto name_ai = atomTypeOpenMM(myatoms[i].ffType(), i);

                    for (const auto &aType : pd->particleTypesConst())
                    {
                        if (aType.id().id() == myatoms[i].ffType()) 
                        {
                            if (aType.gmxParticleType() == eptAtom &&
                                aType.hasOption("poltype"))
                            {
                                auto grandchild2 = add_xml_child(child6, exml_names(xmlEntryOpenMM::PARTICLE)); 
                                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), name_ai.c_str()); 
                                auto shell_ai = aType.optionValue("poltype");
                                auto alpha = fs.second.findParameterTypeConst(Identifier({shell_ai}),
                                                                              pol_name[polALPHA]);
                                
                                add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE1), shell_ai.c_str());
                                add_xml_double(grandchild2, "polarizability", alpha.internalValue());
                                // TODO: Fix atom number for shell
                                add_xml_double(grandchild2, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[i+1].charge());
                                add_xml_double(grandchild2, "thole", 0);
                            }
                        }
                    }
                }    
            }
        }                                              
      
    }
  
}
   

void writeOpenMM(const std::string &fileName,
                 const ForceField  *pd,
                 const ACTMol      *actmol,
                 bool               compress)
{
    xmlDocPtr   doc;
    xmlDtdPtr   dtd;
    xmlNodePtr  myroot;
    xmlChar    *libdtdname, *dtdname, *gmx;

    rmapyyyOpenMM.clear();
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
    addXmlForceField(myroot, pd, actmol);

    xmlSetDocCompressMode(doc, compress ? 1 : 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fileName.c_str(), doc, "ISO-8859-1", 2) == 0)
    {
        gmx_fatal(FARGS, "Saving file %s", fileName.c_str());
    }
    xmlFreeDoc(doc);
}

}
