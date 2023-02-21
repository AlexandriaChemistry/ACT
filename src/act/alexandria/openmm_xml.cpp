/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2023
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
#include <set>

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

#include <algorithm>
#include <list>
#include <iostream>

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
        { ub_name[ubR13],         "d" },
	{ cubic_name[cubicKB],         "kb" },
	{ cubic_name[cubicRMAX],         "rmax" },
        { cubic_name[cubicDE],         "D_e" },
        { cubic_name[cubicLENGTH],         "r0" },	
	//const char *cubic_name[cubicNR] = { "bondlength", "rmax", "kb", "De" };
	//const char *morse_name[morseNR] = { "beta", "De", "D0", "bondlength" };
    };
    if (type == specparam)
    {
        // TODO make this a parameter
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

//static void addShell(xmlNodePtr         parent,
//                     const std::string &key,
//                     const std::string &value,
//                    const std::string &specopt)
//{
//    if (strcmp(key.c_str(), specopt.c_str()) == 0)
//    {    
//        auto baby = add_xml_child(parent, exml_names(xmlEntryOpenMM::ATOM_RES));    
//        add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME), value.c_str());
//        add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), value.c_str());
//	////std::cout << exml_names(xmlEntryOpenMM::NAME) << "asdf \n";
//    }
//}

static void addXmlElemMass(xmlNodePtr parent, const ParticleType *aType, double mDrude)
{
    if (eptAtom == aType->gmxParticleType())
    {
        add_xml_char(parent, exml_names(xmlEntryOpenMM::ELEMENT),
                     aType->element().c_str());
        double mAtom = aType->mass();
        if (aType->hasOption("poltype"))
        {
            mAtom -= mDrude;
        }
        add_xml_double(parent, exml_names(xmlEntryOpenMM::MASS), mAtom);
    }
    else if (eptShell == aType->gmxParticleType())
    {
        add_xml_double(parent, exml_names(xmlEntryOpenMM::MASS), mDrude);
    }
}

static std::string nameIndex(const std::string &name, int index)
{
    std::string ni = name + gmx::formatString("_%d", index);
    return ni;
}
/////// adding+, real atoms are from 1 to 3 for water for cores and shells alike
static int tellme_RealAtom(int index, const ACTMol *actmol)
{       auto myatoms = actmol->atomsConst();
	int        realAtoms     = 0;
	size_t key_index = index;
	int result = -1;

	for (size_t i = 0; i < myatoms.size(); i++)
        	{
/////		std::cout << i << "aaa\n";	
		
		if (myatoms[i].pType() == eptAtom) 
		{
			realAtoms = realAtoms + 1;
		}
		if (i == key_index)
			{ 
				result = realAtoms;
			}
		}
	return result;
}

static std::set<int> get_unique_residues(const ACTMol *actmol)
{
	auto myatoms =  actmol -> atomsConst();
	std::set<std::string> FfTypeUsed, BondClassUsed;
	std::set<int> Residuelist_for_which_loop_atoms {};


        // First residue will be defined below. 
        xmlNodePtr residuePtr    = nullptr;
        int        residueNumber = -1;
        auto       resnames      = actmol->topology()->residueNames();
        std::set<std::string> ResidueUsed;
        std::set<int>         Atoms_used;
    //    bool       skipAtoms     = false;
        //+
//      int        realAtoms     = 0;
        // Let each residue start with atom number 1 within the residue definition, to do this
        // we store the number of the first atom of each residue.
    //    int        residueStart  = 0;
        for (size_t i = 0; i < myatoms.size(); i++)
        {
            if (myatoms[i].residueNumber() != residueNumber)
            {
                // Time for a new residue, but prevent that we repeat them.
                // We have rely on residue names to mean identical chemical moieties,
                // that is if we for instance have an N-terminal amino acid, it should
                // be a different residue from a mid-chain amino acid.
                residueNumber = myatoms[i].residueNumber();

                if (ResidueUsed.find(resnames[residueNumber]) == ResidueUsed.end())
                {
                    // Check whether we have to terminate the residue by defining bonds
                    if (nullptr != residuePtr)
                    {
                      //  addXmlResidueBonds(residuePtr, pd, actmol, Atoms_used, residueStart);
                        Atoms_used.clear();
                    }
                 //   residuePtr = add_xml_child(child2, exml_names(xmlEntryOpenMM::RESIDUE));
                  //  add_xml_char(residuePtr, exml_names(xmlEntryOpenMM::NAME), resnames[residueNumber].c_str());
                    ResidueUsed.insert(resnames[residueNumber]);
                    Atoms_used.insert(i);
              //      residueStart = i;
                  //  skipAtoms    = false;

                    Residuelist_for_which_loop_atoms.insert(residueNumber);
                  //  Residuelist_for_which_loop_atoms.push_back(residueNumber);

                }
              //  else
              //  {
             //       skipAtoms = true;
              //  }
            }

	}

	return Residuelist_for_which_loop_atoms; 





}
///////
static void addXmlResidueBonds(xmlNodePtr residuePtr, const ForceField *pd, const ACTMol *actmol, 
                               const std::set<int> &Atoms_used, int residueStart)
{
    auto itbonds = InteractionType::BONDS;
    if (pd->interactionPresent(itbonds))
    {
        auto fs = pd->findForcesConst(itbonds);
        
        if (actmol->topology()->hasEntry(itbonds))
        {
            auto myatoms = actmol->atomsConst();
           for(const auto topentry : actmol->topology()->entry(itbonds))
	       {
		int ai = topentry->atomIndex(0);
                int aj = topentry->atomIndex(1);
	////-	std::cout << ai << "AI\n";
/////-		std::cout << tellme_RealAtom(ai, actmol) << "REALi\n";
/////-		std::cout << tellme_RealAtom(aj, actmol) << "REALj\n";

		int reali = tellme_RealAtom(ai, actmol);
		int realj = tellme_RealAtom(aj, actmol);
                if (Atoms_used.find(ai) != Atoms_used.end() &&
                    Atoms_used.find(aj) != Atoms_used.end())
                {
                    //here also myatoms[i].ffType(); but for ai aj
                    auto name_ai = myatoms[ai].name();
                    auto name_aj = myatoms[aj].name();
                    auto baby = add_xml_child(residuePtr, exml_names(xmlEntryOpenMM::BOND_RES));
                    add_xml_char(baby, exml_names(xmlEntryOpenMM::ATOMNAME1_RES), 
                                 nameIndex(name_ai, reali - residueStart).c_str());
                    add_xml_char(baby, exml_names(xmlEntryOpenMM::ATOMNAME2_RES),
                                 nameIndex(name_aj, realj - residueStart).c_str()); 
////		    std::cout << exml_names(xmlEntryOpenMM::ATOMNAME1_RES) << "asdf \n";
                }
            }
        }
    }
}

static bool atomsInClass(const std::vector<std::string> &atoms,
                         const std::set<std::string>    &myClass)
{
    bool inClass = true;
    for (auto &a: atoms)
    {
        if (myClass.find(a) == myClass.end())
        {
            inClass = false;
        }
    }
    return inClass;
}

static void addBondAtoms(xmlNodePtr                      parent,
                         const std::vector<std::string> &atoms)
{
    for (size_t i = 0; i < atoms.size(); i++)
    {
        add_xml_char(parent, exml_names(class_vec[i]), atoms[i].c_str());
    }
}

static void addXmlBonds(xmlNodePtr                     parent,
                        const std::set<std::string>   &BondClassUsed,
                        const ForceFieldParameterList &fs)
{
    switch(fs.gromacsType())
    {
	    ///////////////////////////////////////////////////////////////////////////////////
    case F_MORSE:
        {
            auto child3 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMBONDFORCE));
            // The Morse potential could be written as a string here,
            // or it can be added in the openmm python script
            add_xml_double(child3, "energy", 0.0);
            
            // Specify the per bond parameters
            auto grandchild0 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), "D_e");
            auto grandchild1 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), "a");
            auto grandchild2 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::NAME), "r0");
            
            // Add all bonds
            for (auto &params : fs.parametersConst())
            {
                // To filter out bonds, we should not look for the ffType but for the correspondng bond type
                if (atomsInClass(params.first.atoms(), BondClassUsed))
                {
                    auto grandchild3 = add_xml_child(child3, exml_names(xmlEntryOpenMM::BOND_RES));
                    
                    addBondAtoms(grandchild3, params.first.atoms());
                    
                    for (const auto &param : params.second)
                    {
                        addSpecParameter(grandchild3, param.first, param.second, morse_name[morseDE]);
                        addSpecParameter(grandchild3, param.first, param.second, morse_name[morseBETA]);
                        addSpecParameter(grandchild3, param.first, param.second, morse_name[morseLENGTH]);
                    }
                }
            }
            break;
        }
///////////////////////////////////////////////////////////////////////////////////////

    case F_CUBICBONDS:
        {
            auto child3 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMBONDFORCE));
            // The Morse potential could be written as a string here,
            // or it can be added in the openmm python script
            add_xml_double(child3, "energy", 0.0);

            // Specify the per bond parameters
            auto grandchild0 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), "D_e");
            auto grandchild1 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), "kb");
            auto grandchild2 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::NAME), "r0");
            auto grandchild3 = add_xml_child(child3, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
            add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::NAME), "rmax");

            // Add all bonds
            for (auto &params : fs.parametersConst())
            {
                // To filter out bonds, we should not look for the ffType but for the correspondng bond type
                if (atomsInClass(params.first.atoms(), BondClassUsed))
                {
                    auto grandchild4 = add_xml_child(child3, exml_names(xmlEntryOpenMM::BOND_RES));

                    addBondAtoms(grandchild4, params.first.atoms());

                    for (const auto &param : params.second)
                    {
                        addSpecParameter(grandchild4, param.first, param.second, cubic_name[cubicDE]);
                        addSpecParameter(grandchild4, param.first, param.second, cubic_name[cubicKB]);
                        addSpecParameter(grandchild4, param.first, param.second, cubic_name[cubicLENGTH]);
			addSpecParameter(grandchild4, param.first, param.second, cubic_name[cubicRMAX]);
                    }
                }
            }
            break;
        }




    case F_BONDS:
        {
            // TODO Check for correctness
            auto child3 = add_xml_child(parent, exml_names(xmlEntryOpenMM::HARMONICBONDFORCE));
            // The Morse potential could be written as a string here,
            // or it can be added in the openmm python script
            add_xml_double(child3, "energy", 0.0);
            
            // Add all bonds
            for (auto &params : fs.parametersConst())
            {
                if (atomsInClass(params.first.atoms(), BondClassUsed))
                {
                    auto grandchild2 = add_xml_child(child3, exml_names(xmlEntryOpenMM::BOND_RES));
                    
                    addBondAtoms(grandchild2, params.first.atoms());
                    
                    for (const auto &param : params.second)
                    {
                        // TODO check names and order of parameters. 
                        addSpecParameter(grandchild2, param.first, param.second, "length"); //bond_name[bondLENGTH]);
                        addSpecParameter(grandchild2, param.first, param.second, "k"); // bond_name[bondKB]); 
                    }
                }    
            }
            break;
        }
    default:
        {
            gmx_fatal(FARGS, "Do not know how to treat %s in addXmlBonds", fs.function().c_str());
        }
    }
}

static void addXmlAngles(xmlNodePtr                     parent,
                         const std::set<std::string>   &BondClassUsed,
                         const ForceFieldParameterList &fs)
{
    switch(fs.gromacsType())
    {
    case F_ANGLES:
        {
            auto child4 = add_xml_child(parent, exml_names(xmlEntryOpenMM::HARMONICANGLEFORCE));
            
            for (auto &params : fs.parametersConst())
            {
                if (atomsInClass(params.first.atoms(), BondClassUsed))
                {
                    auto grandchild2 = add_xml_child(child4, exml_names(xmlEntryOpenMM::ANGLE_CLASS));
                    
                    addBondAtoms(grandchild2, params.first.atoms());
                    
                    for (const auto &param : params.second)
                    {
                        addSpecParameter(grandchild2, param.first, param.second, angle_name[angleANGLE]); 
                        addSpecParameter(grandchild2, param.first, param.second, angle_name[angleKT]); 
                    }
                }    
            }
        }
        break;
    case F_UREY_BRADLEY:
        {
            auto child4 = add_xml_child(parent, exml_names(xmlEntryOpenMM::AMOEBA_UB_FORCE));
            for (auto &params : fs.parametersConst())
            {
                auto grandchild2 = add_xml_child(child4, exml_names(xmlEntryOpenMM::UREY_BRADLEY_FORCE));
                addBondAtoms(grandchild2, params.first.atoms());
                            
                for (const auto &param : params.second)
                {
                    addSpecParameter(grandchild2, param.first, param.second, "kub"); 
                    addSpecParameter(grandchild2, param.first, param.second, "r13"); 
                    //add_xml_double(grandchild2, param.first.c_str(), param.second.value());
                }
            }
        }
        break;
    default:
        gmx_fatal(FARGS, "Do not know how to treat %s", fs.function().c_str());
    }
}

static void addXmlLinearAngles(xmlNodePtr                     parent,
                               const std::set<std::string>   &BondClassUsed,
                               const ForceFieldParameterList &fs)
{
    auto child4 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMANGLEFORCE));
    add_xml_double(child4, "energy", 0.0);
    for (auto &params : fs.parametersConst())
    {            
        for (const auto &param : params.second)
        {
            if (atomsInClass(params.first.atoms(), BondClassUsed))
            {
                auto grandchild1 = add_xml_child(child4, exml_names(xmlEntryOpenMM::PERANGLEPARAMETER));
                addBondAtoms(grandchild1, params.first.atoms());
                add_xml_char(grandchild1, exml_names(xmlEntryOpenMM::NAME), param.first.c_str());
            }
        }
    }
#ifdef OLD       
    for (auto &params : fs.parametersConst())
    {
        auto grandchild2 = add_xml_child(child4, exml_names(xmlEntryOpenMM::ANGLE_CLASS));
        addBondAtoms(grandchild2, params.first.atoms());
        
        for (const auto &param : params.second)
        {
            add_xml_double(grandchild2, param.first.c_str(), param.second.value());
        }
    }
#endif
}

static void addXmlPropers(xmlNodePtr                     parent,
                          const std::set<std::string>   &BondClassUsed,
                          const ForceFieldParameterList &fs)
{
    auto child6 = add_xml_child(parent, exml_names(xmlEntryOpenMM::RBTORSIONFORCE));
    
    for (auto &params : fs.parametersConst())
    {
        if (atomsInClass(params.first.atoms(), BondClassUsed))
        {
            auto grandchild2 = add_xml_child(child6, exml_names(xmlEntryOpenMM::PROPER));
            addBondAtoms(grandchild2, params.first.atoms());
            
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
}

static void addXmlImpropers(xmlNodePtr                     parent,
                            const std::set<std::string>   &BondClassUsed,
                            const ForceFieldParameterList &fs)
{
    auto child6 = add_xml_child(parent, exml_names(xmlEntryOpenMM::PERIODICTORSIONFORCE));
    
    for (auto &params : fs.parametersConst())
    {
        if (atomsInClass(params.first.atoms(), BondClassUsed))
        {
            auto grandchild2 = add_xml_child(child6, exml_names(xmlEntryOpenMM::IMPROPER));
            addBondAtoms(grandchild2, params.first.atoms());
            
            for (const auto &param : params.second)
            {
                add_xml_double(grandchild2, param.first.c_str(), param.second.value());
            }
        }
    }    
}

static void addXmlNonbonded(xmlNodePtr                     parent,
                            const ForceField              *pd,
                            const ForceFieldParameterList &fs,
                            const std::vector<ActAtom>    &myatoms,
			    const ACTMol *actmol)
{
    auto child5 = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMNONBONDEDFORCE));
    add_xml_double(child5, "energy", 0.0);
    // TODO Make this a parameter? See https://github.com/dspoel/ACT/issues/26
    // Question whether 0 will work.
    add_xml_int(child5, "bondCutoff", 3);
    
    auto grandchild0 = add_xml_child(child5, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
    add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), "vdW");
    
    // This order is important !!! Do not change the order of these parameters,
    // otherwise the force field is not working !!!
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
    std::set<int> Residuelist_for_which_loop_atoms {};
    Residuelist_for_which_loop_atoms = get_unique_residues(actmol);
    // add customnonbonded (WBH vdW + pg coulomb) for compound
    std::set<std::string> List_used;
    for (size_t i = 0; i < myatoms.size(); i++)
    { int reali = tellme_RealAtom(i, actmol); 
        if (List_used.find(myatoms[i].ffType()) == List_used.end())
        {if (Residuelist_for_which_loop_atoms.find(myatoms[i].residueNumber()) != Residuelist_for_which_loop_atoms.end())
                {
         //   List_used.insert(myatoms[i].ffType());
            auto name_ai = myatoms[i].ffType(); //atomTypeOpenMM(myatoms[i].ffType(), i); DELETE THIS VAR
            auto grandchild3 = add_xml_child(child5, exml_names(xmlEntryOpenMM::ATOM_RES));
      //      auto ffType = myatoms[i].ffType();   // for the class definition
//	    add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::CLASS), ffType.c_str());  //for the class definition instead of type
	    add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::TYPE_RES), nameIndex(myatoms[i].ffType(), reali).c_str());  
//            add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::TYPE_RES), name_ai.c_str());  
            for (const auto &aType : pd->particleTypesConst())
            {
                if (aType.id().id() == myatoms[i].ffType())
                {
                    if (eptAtom == aType.gmxParticleType()) 
                    {
                        add_xml_double(grandchild3, "vdW", 1.0); 
                    }
                    if (eptShell == aType.gmxParticleType())
                    {
                        add_xml_double(grandchild3, "vdW", 0.0);
                    }   
                    
                    for(const auto &opt: aType.optionsConst())
                    { 
                        if (strcmp(opt.first.c_str(), "vdwtype") == 0) 
                        {
                            for (auto &params : fs.parametersConst())
                            {
                                for (const auto &param : params.second)
                                {
                                    if (strcmp(opt.second.c_str(), params.first.id().c_str()) == 0)
                                    {    
                                        if (eptAtom == aType.gmxParticleType()) 
                                        {
                                            add_xml_double(grandchild3, param.first.c_str(), param.second.value());
                                        }
                                        if (eptShell == aType.gmxParticleType()) 
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
                                if (InteractionType::COULOMB == fs.first)
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
	 }//  
        } 
    }
        
    // This part is added to implement PME and LJPME, i.e. long-range interactions are approx. by point charge and 12-6 Lennard-Jones
    // For the different atomtypes it would be good to add optimized sigma and epsilon values
    // in order to calculate the vdW contribution after the cutoff
    child5 = add_xml_child(parent, exml_names(xmlEntryOpenMM::NONBONDEDFORCE));
    add_xml_double(child5, "coulomb14scale", 1.0); 
    add_xml_double(child5, "lj14scale", 1.0); 
    
    // for (const auto &aType : pd->particleTypesConst())
    // {
    //        auto grandchild2 = add_xml_child(child5, exml_names(xmlEntryOpenMM::ATOM_RES));
    //        add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::CLASS), aType.id().id().c_str());
    //  for(const auto &param : aType.parametersConst())
    // {
    //      addSpecParameter(grandchild2, param.first, param.second, "charge"); 
    //      if (strcmp( ptype_str[aType.gmxParticleType()], "Atom") == 0)
    //      {
    //          add_xml_double(grandchild2, "sigma", 0.3);  
    //          add_xml_double(grandchild2, "epsilon", 0.05);  
    //      }
    //      if (strcmp( ptype_str[aType.gmxParticleType()], "Shell") == 0)
    //      {
    //          add_xml_double(grandchild2, "sigma", 0.3);  
    //          add_xml_double(grandchild2, "epsilon", 0.0);  
    //      }
    //                    }
    
    //  }
    // adding nonbonded (LJ + point coulomb) for compound
    // Clear the set
    List_used.clear();	
    for (size_t i = 0; i < myatoms.size(); i++)
    {  int reali = tellme_RealAtom(i, actmol);
        if (List_used.find(myatoms[i].ffType()) == List_used.end())
        {if (Residuelist_for_which_loop_atoms.find(myatoms[i].residueNumber()) != Residuelist_for_which_loop_atoms.end())
                {

	       
		//std::cout << myatoms[i].ffType() << "TYPEEEEEEEEEEEE \n";
         //   List_used.insert(myatoms[i].ffType()); 
            
            auto name_ai = myatoms[i].ffType(); //atomTypeOpenMM(myatoms[i].ffType(), i);
            auto baby = add_xml_child(child5, exml_names(xmlEntryOpenMM::ATOM_RES));
            auto ffType = myatoms[i].ffType();   // for the class definition
//	    add_xml_char(baby, exml_names(xmlEntryOpenMM::CLASS), ffType.c_str());  //for the class definition instead of type
	    add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), nameIndex(myatoms[i].ffType(), reali).c_str()); 
	 //   add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), name_ai.c_str()); 
	  
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
    }
}

static void addXmlPolarization(xmlNodePtr                     parent,
                               const ForceField              *pd,
                               const ForceFieldParameterList &fs,
                               const std::vector<ActAtom>    &myatoms,
			       const ACTMol *actmol)
{
    // !!! Shell particle has to be type1, core particle has to be type2 !!!
    // 
    if (pd->polarizable())    
    {
        auto child6 = add_xml_child(parent, exml_names(xmlEntryOpenMM::DRUDEFORCE));
        std::set<std::string> List_used;
        for (size_t i = 0; i < myatoms.size(); i++)
        { int reali = tellme_RealAtom(i, actmol);

            auto ffType = myatoms[i].ffType();
            if (List_used.find(ffType) == List_used.end())
            {
                List_used.insert(ffType); 
                if (eptAtom == myatoms[i].pType())
                {
                    auto ptp     = pd->findParticleType(ffType);
                    if (ptp->hasOption("poltype"))
                    {
                        auto grandchild2 = add_xml_child(child6, exml_names(xmlEntryOpenMM::PARTICLE)); 

                        add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), nameIndex(myatoms[i].ffType(), reali).c_str()); 
                    //    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), ffType.c_str()); 
                        auto shell_ai = ptp->optionValue("poltype");
                        auto alpha    = fs.findParameterTypeConst(Identifier({shell_ai}),
                                                                  pol_name[polALPHA]);
	//		auto shell = myatoms[i].shells();
		//	auto shelli = 
//			std::cout << myatoms[i] << "ai \n";

                     //  std::cout << myatoms[i].shells() << "ai \n";
		    //    const auto shell = myatoms[i].shells();
                    for(const auto &shell: myatoms[i].shells())

		    	{ 
			//	std::cout << typeid(myatoms[i]).name() << typeid(myatoms[i].shells()).name() << "si \n";
			//	if (myatoms[i] == shell)
			//	{
                        	add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE1), nameIndex(myatoms[shell].ffType(), reali).c_str());

                        	add_xml_double(grandchild2, "polarizability", alpha.internalValue());
                        	auto stp = pd->findParticleType(shell_ai);
                        	add_xml_double(grandchild2, exml_names(xmlEntryOpenMM::CHARGE_RES), stp->charge());
                        	add_xml_double(grandchild2, "thole", 0);
				}
		//	}	
                    }
                }
            }    
        }
    }
}

static void addXmlForceField(xmlNodePtr parent, const ForceField *pd, const ACTMol *actmol, double mDrude)
{
    auto child = add_xml_child(parent, exml_names(xmlEntryOpenMM::ATOMTYPES));
    auto myatoms =  actmol -> atomsConst();
    std::set<std::string> FfTypeUsed, BondClassUsed;
    std::set<int> Residuelist_for_which_loop_atoms {};
    Residuelist_for_which_loop_atoms = get_unique_residues(actmol);
    for (size_t i = 0; i < myatoms.size(); i++)
    {
	for (auto i = Residuelist_for_which_loop_atoms.begin(); i != Residuelist_for_which_loop_atoms.end(); i++)
                                {  // std::cout << *i << "this is the number of residues in list";
                                       }
        auto ffType = myatoms[i].ffType();
        if (FfTypeUsed.find(ffType) == FfTypeUsed.end())
        {
		if (Residuelist_for_which_loop_atoms.find(myatoms[i].residueNumber()) != Residuelist_for_which_loop_atoms.end())
		{ 
				
                

//+        //    FfTypeUsed.insert(ffType); REMOVE THIS CHECK ALTOGETHER
            int reali = tellme_RealAtom(i, actmol);
            auto baby = add_xml_child(child, exml_names(xmlEntryOpenMM::TYPE));
            add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME), nameIndex(myatoms[i].ffType(), reali).c_str());
            add_xml_char(baby, exml_names(xmlEntryOpenMM::CLASS), ffType.c_str());
            
            if (!pd->hasParticleType(ffType))
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("No such particle type %s in force field %s", ffType.c_str(), pd->filename().c_str()).c_str()));
            }
            auto aType = pd->findParticleType(ffType);
            addXmlElemMass(baby, &(*aType), mDrude);
            std::string btype("bondtype");
            if (aType->hasOption(btype))
            {
                BondClassUsed.insert(aType->optionValue(btype));
            }
        }
    }
   //  }	
    }
    auto child2  = add_xml_child(parent, exml_names(xmlEntryOpenMM::RESIDUES));

    if (myatoms.size() > 0)
    {
        // First residue will be defined below. 
        xmlNodePtr residuePtr    = nullptr;
        int        residueNumber = -1;
        auto       resnames      = actmol->topology()->residueNames();
        std::set<std::string> ResidueUsed;
        std::set<int>         Atoms_used;
        bool       skipAtoms     = false;
	//+
//	int        realAtoms     = 0;
        // Let each residue start with atom number 1 within the residue definition, to do this
        // we store the number of the first atom of each residue.
        int        residueStart  = 0;
        for (size_t i = 0; i < myatoms.size(); i++)
        {   
            if (myatoms[i].residueNumber() != residueNumber)
            {    
                // Time for a new residue, but prevent that we repeat them.
                // We have rely on residue names to mean identical chemical moieties,
                // that is if we for instance have an N-terminal amino acid, it should
                // be a different residue from a mid-chain amino acid.
                residueNumber = myatoms[i].residueNumber();
            
                if (ResidueUsed.find(resnames[residueNumber]) == ResidueUsed.end())
                {
                    // Check whether we have to terminate the residue by defining bonds
                    if (nullptr != residuePtr)
                    {
                        addXmlResidueBonds(residuePtr, pd, actmol, Atoms_used, residueStart);
                        Atoms_used.clear();
                    }
                    residuePtr = add_xml_child(child2, exml_names(xmlEntryOpenMM::RESIDUE));
                    add_xml_char(residuePtr, exml_names(xmlEntryOpenMM::NAME), resnames[residueNumber].c_str());
                    ResidueUsed.insert(resnames[residueNumber]);
                    Atoms_used.insert(i);
                    residueStart = i;
                    skipAtoms    = false;
		   // std::cout << "resolves" << " ";

		    Residuelist_for_which_loop_atoms.insert(residueNumber);
	//	    for (auto i = Residuelist_for_which_loop_atoms.begin(); i != Residuelist_for_which_loop_atoms.end(); i++)
	//	    {   std::cout << *i << "this is the number ";
	//			}

                }
                else
                {
                    skipAtoms = true;
                }
            }
            if (!skipAtoms)
            {
		//    std::cout << realAtoms << "core \n";
		//    std::cout << i << "I \n";
                int reali = tellme_RealAtom(i, actmol);

                if (myatoms[i].pType() == eptAtom)
                {
			//realAtoms = realAtoms + 1;   
                    auto baby = add_xml_child(residuePtr, exml_names(xmlEntryOpenMM::ATOM_RES));
                    add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME),
                                 nameIndex(myatoms[i].name(), reali).c_str());

                    add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), nameIndex(myatoms[i].ffType(), reali).c_str());
                   // add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), nameIndex(myatoms[i].ffType(), reali - residueStart).c_str());  // this messes up numbering of different residues
                    add_xml_double(baby, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[i].charge());

                    for(const auto &shell: myatoms[i].shells())
                    {
                        auto baby = add_xml_child(residuePtr, exml_names(xmlEntryOpenMM::ATOM_RES));
                        add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME),
                                     nameIndex(myatoms[shell].ffType(), reali).c_str());
                        add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), nameIndex(myatoms[shell].ffType(), reali).c_str());
		//	std::cout << shell << "myshell \n";
			//add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), myatoms[shell].ffType().c_str());
                        add_xml_double(baby, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[shell].charge());     
                    }
                    Atoms_used.insert(i);
                }
            }
        }
        // Check whether we have to terminate the residue by defining bonds
        if (nullptr != residuePtr)
        {
            addXmlResidueBonds(residuePtr, pd, actmol, Atoms_used, residueStart);
            Atoms_used.clear();
        }
    }
    
    for (auto &fs : pd->forcesConst())
    {
        switch(fs.first)
        {
        case InteractionType::BONDS:
            // This adds the Morse potential
            // TODO: link the bondorder information to the parameter assignment that happens via atomtypes
            addXmlBonds(parent, BondClassUsed, fs.second);
            break;
            
        case InteractionType::ANGLES:
            addXmlAngles(parent, BondClassUsed, fs.second);
            break;
            
        case InteractionType::LINEAR_ANGLES:
            addXmlLinearAngles(parent, BondClassUsed, fs.second);
            break;
            
        case InteractionType::PROPER_DIHEDRALS:
            addXmlPropers(parent, BondClassUsed, fs.second);
            break;

        case InteractionType::IMPROPER_DIHEDRALS:
            addXmlImpropers(parent, BondClassUsed, fs.second);
            break;
        
        case InteractionType::VDW:
            addXmlNonbonded(parent, pd, fs.second, myatoms, actmol);
            break;
        case InteractionType::POLARIZATION:
            addXmlPolarization(parent, pd, fs.second, myatoms, actmol);
            break;
        case InteractionType::COULOMB:
            // Coulomb taken into account with other non-bonded interactions
        case InteractionType::BONDCORRECTIONS:
        case InteractionType::ELECTRONEGATIVITYEQUALIZATION:
            break;
        default:
            fprintf(stderr, "Wanrning: no OpenMM support for %s is present or implemented.\n",
                    interactionTypeToString(fs.first).c_str());
        }
    }
}

void writeOpenMM(const std::string &fileName,
                 const ForceField  *pd,
                 const ACTMol      *actmol,
                 double             mDrude,
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
    addXmlForceField(myroot, pd, actmol,mDrude);

    xmlSetDocCompressMode(doc, compress ? 1 : 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fileName.c_str(), doc, "ISO-8859-1", 2) == 0)
    {
        gmx_fatal(FARGS, "Saving file %s", fileName.c_str());
    }
    xmlFreeDoc(doc);
}

}
