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

static std::string nameIndex(const std::string &name, int index)
{
    std::string ni = name + gmx::formatString("_%d", index);
    return ni;
}

/////// adding+, real atoms are from 1 to 3 for water for cores and shells alike
static int tellme_RealAtom(int                         index, 
                           const std::vector<ActAtom> &myatoms)
{
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

/*! Put code in class for easier structuring
 */
class OpenMMWriter
{
private:
    // Mass of the drude particle
    double mDrude_                = 0.1;
    // Whether or not to add numbers to atomtypes
    bool   addNumbersToAtomTypes_ = true;
    // Number of exclusions
    int    nrexcl_                = 3;     
    // Mapping from InteractionType to xml branch
    std::map<InteractionType, xmlNodePtr> xmlMap_;
    
    void addXmlElemMass(xmlNodePtr parent,
                        const ParticleType &aType);
    
    void addXmlResidueBonds(xmlNodePtr           residuePtr,
                            const ForceField    *pd,
                            const Topology      &topol);

    void addBondAtoms(xmlNodePtr                      parent,
                      const std::vector<std::string> &atoms);
                      
    void addXmlNonbonded(xmlNodePtr                       parent,
                         const ForceField                *pd,
                         const std::map<std::string, int> &ffTypeMap);

    void addXmlBond(xmlNodePtr                      parent,
                    xmlEntryOpenMM                  xmlEntry,
                    const std::vector<std::string> &atoms,
                    const char                     *param_names[],
                    const std::vector<double>      &params);
    
    void addXmlPolarization(xmlNodePtr                       parent,
                            const ForceField                *pd,
                            const std::map<std::string, int> &ffTypeMap);
                            
    void addXmlForceField(xmlNodePtr                 parent,
                          const ForceField          *pd,
                          const std::vector<ACTMol> &actmols);
                          
    void addTopologyEntries(const ForceField                                  *pd,
                            std::map<InteractionType, std::set<std::string> > *BondClassUsed,
                            const Topology                                    *topology);
                          
    void makeXmlMap(xmlNodePtr        parent,
                    const ForceField *pd);
               
public:
    OpenMMWriter(double mDrude,
                 bool   addNumbersToAtoms) : mDrude_(mDrude), addNumbersToAtomTypes_(addNumbersToAtoms)
    {
    }
    
    /*! \brief Do the writing
     */
    void write(const std::string         &fileName,
               const ForceField          *pd,
               const std::vector<ACTMol> &actmols,
               bool                       compress);
};

void OpenMMWriter::addXmlElemMass(xmlNodePtr parent, const ParticleType &aType)
{
    if (eptAtom == aType.gmxParticleType())
    {
        add_xml_char(parent, exml_names(xmlEntryOpenMM::ELEMENT),
                     aType.element().c_str());
        double mAtom = aType.mass();
        if (aType.hasOption("poltype"))
        {
            mAtom -= mDrude_;
        }
        add_xml_double(parent, exml_names(xmlEntryOpenMM::MASS), mAtom);
    }
    else if (eptShell == aType.gmxParticleType())
    {
        add_xml_double(parent, exml_names(xmlEntryOpenMM::MASS), mDrude_);
    }
}


void OpenMMWriter::addXmlResidueBonds(xmlNodePtr           residuePtr, 
                                      const ForceField    *pd, 
                                      const Topology      &topol)
{
    auto itbonds = InteractionType::BONDS;
    if (pd->interactionPresent(itbonds))
    {
        auto fs = pd->findForcesConst(itbonds);
        
        if (topol.hasEntry(itbonds))
        {
            auto myatoms = topol.atoms();
            for(const auto topentry : topol.entry(itbonds))
            {
                int ai = topentry->atomIndex(0);
                int aj = topentry->atomIndex(1);

                int reali = tellme_RealAtom(ai, myatoms);
                int realj = tellme_RealAtom(aj, myatoms);
                //here also myatoms[i].ffType(); but for ai aj
                auto name_ai = myatoms[ai].name();
                auto name_aj = myatoms[aj].name();
                auto baby = add_xml_child(residuePtr, exml_names(xmlEntryOpenMM::BOND_RES));
                add_xml_char(baby, exml_names(xmlEntryOpenMM::ATOMNAME1_RES), 
                             nameIndex(name_ai, reali).c_str());
                add_xml_char(baby, exml_names(xmlEntryOpenMM::ATOMNAME2_RES),
                             nameIndex(name_aj, realj).c_str()); 
            }
        }
    }
}

void OpenMMWriter::addBondAtoms(xmlNodePtr                      parent,
                                const std::vector<std::string> &atoms)
{
    for (size_t i = 0; i < atoms.size(); i++)
    {
        add_xml_char(parent, exml_names(class_vec[i]), atoms[i].c_str());
    }
}

void OpenMMWriter::addXmlBond(xmlNodePtr                      parent,
                              xmlEntryOpenMM                  xmlEntry,
                              const std::vector<std::string> &atoms,
                              const char                     *param_names[],
                              const std::vector<double>      &params)
{
    // Add some bond parameters
    auto grandchild3 = add_xml_child(parent, exml_names(xmlEntry));
                    
    addBondAtoms(grandchild3, atoms);
    for(size_t i = 0; i < params.size(); i++)
    {
        add_xml_double(grandchild3, param_names[i], params[i]);
    }
}

void OpenMMWriter::addXmlNonbonded(xmlNodePtr                       parent,
                                   const ForceField                *pd,
                                   const std::map<std::string, int> &ffTypeMap)
{
    auto fs     = pd->findForcesConst(InteractionType::VDW);
    auto fsCoul = pd->findForcesConst(InteractionType::COULOMB);
    // Custom non-bonded force
    auto fsPtr  = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMNONBONDEDFORCE));
    add_xml_int(fsPtr, "bondCutoff", nrexcl_);
    
    // This part is added to implement PME and LJPME, i.e. long-range interactions are approx. by point charge and 12-6 Lennard-Jones
    // For the different atomtypes it would be good to add optimized sigma and epsilon values
    // in order to calculate the vdW contribution after the cutoff
    auto ljPtr        = add_xml_child(parent, exml_names(xmlEntryOpenMM::NONBONDEDFORCE));
    add_xml_double(ljPtr, "coulomb14scale", 1.0); 
    add_xml_double(ljPtr, "lj14scale", 1.0); 
    
    auto grandchild0 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
    add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), "vdW");
    
    // This order is important !!! Do not change the order of these parameters,
    // otherwise the force field is not working !!!
    auto grandchild2 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::NAME), "sigma");
    auto grandchild3 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
    add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::NAME), "epsilon");
    auto grandchild4 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
    add_xml_char(grandchild4, exml_names(xmlEntryOpenMM::NAME), "gamma");
    
    auto grandchild5 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
    add_xml_char(grandchild5, exml_names(xmlEntryOpenMM::NAME), "charge");
    auto grandchild6 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
    add_xml_char(grandchild6, exml_names(xmlEntryOpenMM::NAME), "beta");

    // add customnonbonded (WBH vdW + pg coulomb) for compound
    for(const auto &fft: ffTypeMap)
    {
        auto aType = pd->findParticleType(fft.first);
        for(int i = 1; i <= fft.second; i++)
        {
            auto param = fs.findParametersConst(Identifier(fft.first));
            
            auto grandchild3 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::ATOM_RES));
            add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::TYPE_RES), nameIndex(fft.first, i).c_str());  
            if (eptAtom == aType->gmxParticleType()) 
            {
                add_xml_double(grandchild3, "vdW", 1.0); 
            }
            else if (eptShell == aType->gmxParticleType())
            {
                add_xml_double(grandchild3, "vdW", 0.0);
            }
            switch (fs.gromacsType())
            {
            case F_BHAM:
                for(size_t j = 0; j < param.size(); j++)
                {
                    if (Mutability::Dependent != param[wbh_name[j]].mutability())
                    {
                        add_xml_double(grandchild3, wbh_name[j], param[wbh_name[j]].internalValue());
                    }
                }
                break;
            case F_GBHAM:
                for(size_t j = 0; j < param.size(); j++)
                {
                    if (Mutability::Dependent != param[gbh_name[j]].mutability())
                    {
                        add_xml_double(grandchild3, gbh_name[j], param[gbh_name[j]].internalValue());
                    }
                }
                break;
            default:
                fprintf(stderr, "Unknown non-bonded force type %d\n", fs.gromacsType());
            }
            auto ztype = "zetatype";
            if (aType->hasOption(ztype))
            {
                auto   ztp  = aType->optionValue(ztype);
                double zeta = 0;
                if (fsCoul.parameterExists(ztp))
                {
                    zeta = fsCoul.findParameterTypeConst(ztp, "zeta").internalValue();
                }
                add_xml_double(grandchild3, exml_names(xmlEntryOpenMM::CHARGE_RES), aType->charge()); 
                add_xml_double(grandchild3, "zeta", zeta); 
            }
            
            // Normal Lennard Jones
            {
                auto grandchild3 = add_xml_child(ljPtr, exml_names(xmlEntryOpenMM::ATOM_RES));
                add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::TYPE_RES), nameIndex(fft.first, i).c_str());  
                std::vector<double> se_param = { 0.3, 0.05 };
                const char *se_name[] = { "sigma", "epsilon" };
                if (eptShell == aType->gmxParticleType())
                {
                    se_param[1] = 0;
                }

                for(size_t j = 0; j < se_param.size(); j++)
                {
                    add_xml_double(grandchild3, se_name[j], se_param[j]);
                }
            }
        }  
    }    
}

void OpenMMWriter::addXmlPolarization(xmlNodePtr                       parent,
                                      const ForceField                *pd,
                                      const std::map<std::string, int> &ffTypeMap)
{
    if (!pd->polarizable())
    {
        return;
    }
    // !!! Shell particle has to be type1, core particle has to be type2 !!!
    auto fs     = pd->findForcesConst(InteractionType::POLARIZATION);
    auto polPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::DRUDEFORCE));

    for(const auto &fft: ffTypeMap)
    {
        auto aType = pd->findParticleType(fft.first);
        for(int i = 1; i <= fft.second; i++)
        {
            if (eptAtom == aType->gmxParticleType())
            {
                if (aType->hasOption("poltype"))
                {
                    auto grandchild2 = add_xml_child(polPtr, exml_names(xmlEntryOpenMM::PARTICLE)); 
                            
                    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), nameIndex(fft.first, i).c_str()); 
                    auto shell_ai = aType->optionValue("poltype");
                    auto param = fs.findParametersConst(Identifier(shell_ai));
            
                    auto alpha    = fs.findParameterTypeConst(Identifier({shell_ai}), pol_name[polALPHA]);
                    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE1), nameIndex(shell_ai, i).c_str());
                        
                    add_xml_double(grandchild2, "polarizability", alpha.internalValue());
                    auto stp = pd->findParticleType(shell_ai);
                    add_xml_double(grandchild2, exml_names(xmlEntryOpenMM::CHARGE_RES), stp->charge());
                    add_xml_double(grandchild2, "thole", 0);
                }
            }
        }
    }
}

void OpenMMWriter::makeXmlMap(xmlNodePtr        parent,
                              const ForceField *pd)
{
    for(const auto &fs: pd->forcesConst())
    {
        xmlNodePtr fsPtr = nullptr;
        switch(fs.second.gromacsType())
        {
        case F_BONDS:
            {
                fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::HARMONICBONDFORCE));
                add_xml_double(fsPtr, "energy", 0.0);
            }
            break;
        case F_MORSE:
            {
                fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMBONDFORCE));
                add_xml_double(fsPtr, "energy", 0.0);
                // Specify the per bond parameters
                for(int i = 0; i < morseNR; i++)
                {
                    auto grandchild0 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
                    add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), morse_name[i]);
                }
            }
            break;
        case F_CUBICBONDS:
            {
                fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMBONDFORCE));
                // The cubic bonds potential could be written as a string here,
                // or it can be added in the openmm python script
                add_xml_double(fsPtr, "energy", 0.0);
                
                // Specify the per bond parameters
                for(int i = 0; i < cubicNR; i++)
                {
                    auto grandchild0 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
                    add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), cubic_name[i]);
                }
            }
            break;
        case F_ANGLES:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::HARMONICANGLEFORCE));
            break;
        case F_LINEAR_ANGLES:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMANGLEFORCE));
            break;
        case F_UREY_BRADLEY:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMANGLEFORCE));
            break;
        case F_PDIHS:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::RBTORSIONFORCE));
            break;
        case F_IDIHS:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::PERIODICTORSIONFORCE));
            break;
        case F_POLARIZATION:
        case F_BHAM:
            break;
        default:
            // Do nothing
            break;
        }
        if (fsPtr)
        {
            xmlMap_.insert({fs.first, fsPtr});
        }
    }
}

void OpenMMWriter::addTopologyEntries(const ForceField                                  *pd,
                                      std::map<InteractionType, std::set<std::string> > *BondClassUsed,
                                      const Topology                                    *topology)
{
    for (auto &fs : pd->forcesConst())
    {
        if (!topology->hasEntry(fs.first))
        {
            continue;
        }
        if (BondClassUsed->end() == BondClassUsed->find(fs.first))
        {
            // Add empty set
            BondClassUsed->insert({fs.first, {}});
        }
        auto &ClassUsed = BondClassUsed->find(fs.first)->second;
        
        for(auto &entry : topology->entry(fs.first))
        {
            const auto &bondId = entry->id().id();
            const auto &atoms  = entry->id().atoms();
            if (ClassUsed.end() == ClassUsed.find(bondId))
            {
                ClassUsed.insert(bondId);
                switch(fs.second.gromacsType())
                {
                case F_BONDS:
                    {
                        const char *omm_bonds[] = { "k", "length", "De" };
                        addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::BOND_RES,
                                   atoms, omm_bonds, entry->params());
                    }
                    break;
                case F_CUBICBONDS:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::BOND_RES,
                               atoms, cubic_name, entry->params());
                    break;
                case F_MORSE:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::BOND_RES,
                               atoms, morse_name, entry->params());
                    break;
                case F_ANGLES:
                    {
                        const char *omm_angles[] = { "k", "angle" };
                        addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::ANGLE_CLASS,
                                   atoms, omm_angles, entry->params());
                    }
                    break;
                case F_UREY_BRADLEY:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::ANGLE_CLASS,
                               atoms, ub_name, entry->params());
                    break;
                    
                case F_LINEAR_ANGLES:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::ANGLE_CLASS,
                               atoms, linang_name, entry->params());
                    break;
                    
                case F_PDIHS:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::PROPER,
                               atoms, fdih_name, entry->params());
                    break;
                    
                case F_IDIHS:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::IMPROPER,
                               atoms, idih_name, entry->params());
                    break;
                    
                case F_COUL_SR:
                case F_LJ:
                case F_BHAM:
                case F_POLARIZATION:
                    break;
                default:
                    fprintf(stderr, "Wanrning: no OpenMM support for %s is present or implemented.\n",
                            interactionTypeToString(fs.first).c_str());
                }
            }
        }
    }
 }

void OpenMMWriter::addXmlForceField(xmlNodePtr                 parent,
                                    const ForceField          *pd,
                                    const std::vector<ACTMol> &actmols)
{
    // We have to make a list of atomtypes, and residue definitions at the same time.
    // Therefore we need to keep track of multiple branches in the xml tree at the same
    // time.
    auto xmlAtomTypePtr     = add_xml_child(parent, exml_names(xmlEntryOpenMM::ATOMTYPES));
    auto xmlResiduePtr      = add_xml_child(parent, exml_names(xmlEntryOpenMM::RESIDUES));
    makeXmlMap(parent, pd);
    // See comment below about fftypeLocalMap.
    std::map<std::string, int> fftypeGlobalMap;
    // List of compounds that we have encountered
    std::set<std::string> InchiUsed;
    // The bond types used.
    std::map<InteractionType, std::set<std::string> > BondClassUsed;
    // Loop over all molecules
    for(auto &actmol : actmols)
    {
        auto fragmentHandler  = actmol.fragmentHandler();
        auto fragIds          = fragmentHandler->ids();
        // First residue will be defined below.
        xmlNodePtr residuePtr = nullptr;
        auto topologies       = fragmentHandler->topologies();
        // The atoms that are part of the present residue
        std::set<int>         AtomsInResidue;
        for(size_t fff = 0; fff < topologies.size(); fff++)
        {
            if (InchiUsed.find(fragIds[fff]) != InchiUsed.end())
            {
                continue;
            }
            // Time for a new residue, but prevent that we repeat them.
            // We have to rely on residue names to indicate identical chemical moieties,
            // that is if we for instance have an N-terminal amino acid, it should
            // be a different residue from a mid-chain amino acid.
            // This is taken care of under the hood when reading through the OpenBabel
            // interface, which replaces the residue names by InChi strings.
            auto myatoms       = topologies[fff].atoms();
                    
            // Check whether we have to terminate the residue by defining bonds
            residuePtr = add_xml_child(xmlResiduePtr, exml_names(xmlEntryOpenMM::RESIDUE));
            add_xml_char(residuePtr, exml_names(xmlEntryOpenMM::NAME), fragIds[fff].c_str());
            InchiUsed.insert(fragIds[fff]);
            // If we need to add numbers to types to distinguish types within a compound,
            // then we need to store what we have used. The map is from force field type
            // to a number index. It restarts from 1 for every molecule, then if we add a
            // new copy of the same force field type we increase the number.
            std::map<std::string, int> fftypeLocalMap;
            
            for (size_t i = 0; i < myatoms.size(); i++)
            {
                // First do atom type stuff
                auto ffType            = myatoms[i].ffType();
                std::string openMMtype = ffType;
                // Do we have copies of this atom type in the local map already?
                int  localIndex = 0;
                auto ffLocalPtr = fftypeLocalMap.find(ffType);
                if (fftypeLocalMap.end() != ffLocalPtr)
                {
                    // Increase localPtr if needed and store it in variable
                    if (addNumbersToAtomTypes_)
                    {
                        ffLocalPtr->second++;
                    }
                    localIndex = ffLocalPtr->second;
                }
                // First check whether this type is in the global map
                auto ffGlobalPtr = fftypeGlobalMap.find(ffType);
                if (fftypeGlobalMap.end() == ffGlobalPtr || 
                    localIndex > ffGlobalPtr->second)
                {
                    // Check whether this type exist in the force field
                    if (!pd->hasParticleType(ffType))
                    {
                        GMX_THROW(gmx::InternalError(gmx::formatString("No such particle type %s in force field %s",
                                                                       ffType.c_str(), pd->filename().c_str()).c_str()));
                    }
                    auto aType = pd->findParticleType(ffType);
                    
                    // New global atomtype?
                    if (fftypeGlobalMap.end() == ffGlobalPtr)
                    {
                        fftypeGlobalMap.insert({ ffType, 0});
                        ffGlobalPtr = fftypeGlobalMap.find(ffType);
                        fftypeLocalMap.insert({ffType, 1});
                        localIndex = 1;
                    }
                    // Make openMM type
                    if (addNumbersToAtomTypes_)
                    {
                        ffGlobalPtr->second++;
                        openMMtype += gmx::formatString("_%d", ffGlobalPtr->second);
                        // printf("openMMtype %s counter %d\n", openMMtype.c_str(), ffGlobalPtr->second);
                    }
                    std::string openMMclass = ffType;
                    std::string btype("bondtype");
                    if (aType->hasOption(btype))
                    {
                        openMMclass = aType->optionValue(btype);
                    }
                    
                    auto newAtype = add_xml_child(xmlAtomTypePtr, exml_names(xmlEntryOpenMM::TYPE));
                    add_xml_char(newAtype, exml_names(xmlEntryOpenMM::NAME), openMMtype.c_str());
                    add_xml_char(newAtype, exml_names(xmlEntryOpenMM::CLASS), openMMclass.c_str());
                    addXmlElemMass(newAtype, *aType);
                }
                else
                {
                    // Re-use an existing atomtype within this residue
                    // Make openMM type
                    if (addNumbersToAtomTypes_)
                    {
                        openMMtype = gmx::formatString("%s_%d", ffType.c_str(), localIndex);
                    }
                }
                int reali = tellme_RealAtom(i, myatoms);

                if (myatoms[i].pType() == eptAtom)
                {
                    auto baby = add_xml_child(residuePtr, exml_names(xmlEntryOpenMM::ATOM_RES));
                    add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME),
                                 nameIndex(myatoms[i].name(), reali).c_str());
                    
                    add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), openMMtype.c_str());
                    add_xml_double(baby, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[i].charge());
                    
                    auto nshells  = myatoms[i].shells().size();
                    size_t ishell = 0;
                    for(const auto &shell: myatoms[i].shells())
                    {
                        auto baby      = add_xml_child(residuePtr, exml_names(xmlEntryOpenMM::ATOM_RES));
                        auto shellName = nameIndex(myatoms[shell].ffType(), reali);
                        std::string shellExt;
                        if (nshells > 1)
                        {
                            shellExt   = gmx::formatString("%c", static_cast<char>('a'+ishell));
                            shellName += shellExt;
                        }
                        auto shellType = myatoms[shell].ffType();
                        if (addNumbersToAtomTypes_)
                        {
                            shellType += gmx::formatString("_%d%s", localIndex, shellExt.c_str());
                        }
                        add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME), shellName.c_str());
                        add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), shellType.c_str());
                        add_xml_double(baby, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[shell].charge());
                        ishell += 1;
                    }
                }
            }
            // Add bonds
            addXmlResidueBonds(residuePtr, pd, topologies[fff]);
            // Add bond types.
            addTopologyEntries(pd, &BondClassUsed, actmol.topology());
        }
    }
    addXmlNonbonded(parent, pd, fftypeGlobalMap);
    addXmlPolarization(parent, pd, fftypeGlobalMap);
}

void OpenMMWriter::write(const std::string         &fileName,
                         const ForceField          *pd,
                         const std::vector<ACTMol> &actmols,
                         bool                       compress)
{
    xmlDocPtr   doc;
    xmlDtdPtr   dtd;
    xmlNodePtr  myroot;
    xmlChar    *libdtdname, *dtdname, *gmx;
    
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
    addXmlForceField(myroot, pd, actmols);

    xmlSetDocCompressMode(doc, compress ? 1 : 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fileName.c_str(), doc, "ISO-8859-1", 2) == 0)
    {
        gmx_fatal(FARGS, "Saving file %s", fileName.c_str());
    }
    xmlFreeDoc(doc);
}


void writeOpenMM(const std::string         &fileName,
                 const ForceField          *pd,
                 const std::vector<ACTMol> &actmols,
                 double                     mDrude,
                 bool                       compress,
                 bool                       addNumbersToAtoms)
{
    rmapyyyOpenMM.clear();
    
    OpenMMWriter writer(mDrude, addNumbersToAtoms);
    writer.write(fileName, pd, actmols, compress);
}

} // namespace alexandria
