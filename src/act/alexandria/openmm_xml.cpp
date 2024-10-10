/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2024
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

#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <set>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "act/alexandria/actmol.h"
#include "act/basics/allmols.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/symcharges.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/potential.h"
#include "act/forces/combinationrules.h"
#include "act/molprop/molprop_util.h"
#include "act/utility/stringutil.h"
#include "act/utility/xml_util.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

/*! \brief The different entry types that can be found
 * TODO: Comment each element
 * xmlEntryOpenMM defines the keys that are in the xmlyyyOpenMM map
 * xmlyyyOpenMM is the map listing value and key used for the OpenMM xml file
 * exml_names returns string value upon calling function and submitting the key 
 * addSpecOption adds specific option for an atomtype 
 * addShell adds shell particle for a specific core atomtype
 * addXmlForceField creates the xml tree ForceField, it consists of Atomtypes, 
 * Residues, HarmonicBondForce, HarmonicAngleForce, CustomNonBondedForce, NonBondedForce, CustomBondForce and DrudeForce.
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
    VSITE_RES,
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
    VSITE2_CLASS,
    VSITE3_CLASS,
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
    USEATTRIBUTEFROMRESIDUE,
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
    { "VirtualSite",               xmlEntryOpenMM::VSITE_RES        },
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
    { "Vsite2",                    xmlEntryOpenMM::VSITE2_CLASS       },
    { "Vsite3",                    xmlEntryOpenMM::VSITE3_CLASS       },
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
    { "UseAttributeFromResidue",   xmlEntryOpenMM::USEATTRIBUTEFROMRESIDUE },
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
		if (myatoms[i].pType() == ActParticle::Atom) 
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
    // Mapping from InteractionType to xml branch
    std::map<InteractionType, xmlNodePtr> xmlMap_;
    
    void addXmlElemMass(xmlNodePtr parent,
                        const ParticleType *aType);
    
    void addXmlResidueBonds(xmlNodePtr           residuePtr,
                            const ForceField    *pd,
                            const Topology      *topol);

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
    
    void addXmlPolarization(xmlNodePtr                        parent,
                            const ForceField                 *pd,
                            const std::map<std::string, int> &ffTypeMap,
                            double                            epsr_fac);

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
    
    /*! \brief Do the writing of the force field file
     * \param[in] fileName The name of the file
     * \param[in] pd       The ACT force field
     * \param[in] actmols  The ACT molecules
     */
    void writeXml(const std::string         &fileName,
                  const ForceField          *pd,
                  const std::vector<ACTMol> &actmols);

    /*! \brief Do the writing of the parameter file
     * \param[in] fileName The name of the file
     * \param[in] pd       The ACT force field
     */
    void writeDat(const std::string         &fileName,
                  const ForceField          *pd);
};

void OpenMMWriter::addXmlElemMass(xmlNodePtr parent, const ParticleType *aType)
{
    switch (aType->apType())
    {
    case ActParticle::Atom:
        {
            add_xml_char(parent, exml_names(xmlEntryOpenMM::ELEMENT),
                         aType->element().c_str());
            double mAtom = aType->mass();
            if (aType->hasOption("poltype"))
            {
                mAtom -= mDrude_;
            }
            add_xml_double(parent, exml_names(xmlEntryOpenMM::MASS), mAtom);
        }
        break;
    case ActParticle::Shell:
        add_xml_double(parent, exml_names(xmlEntryOpenMM::MASS), mDrude_);
        break;
    case ActParticle::Vsite:
        if (0 != aType->mass())
        {
            fprintf(stderr, "Warning: mass should be zero for vsites in OpenMM, not %g. Setting it to zero.\n", aType->mass());
        }
        add_xml_double(parent, exml_names(xmlEntryOpenMM::MASS), 0.0);
        break;
    default:
        GMX_THROW(gmx::InternalError("Don't know what to do. Help!\n"));
    }
}

static void addXmlResidueBond(xmlNodePtr         residuePtr,
                              const std::string &atom1,
                              const std::string &atom2)
{
    auto baby = add_xml_child(residuePtr, exml_names(xmlEntryOpenMM::BOND_RES));
    add_xml_char(baby, exml_names(xmlEntryOpenMM::ATOMNAME1_RES), 
                 atom1.c_str());
    add_xml_char(baby, exml_names(xmlEntryOpenMM::ATOMNAME2_RES),
                 atom2.c_str()); 
}
                              
void OpenMMWriter::addXmlResidueBonds(xmlNodePtr        residuePtr,
                                      const ForceField *pd, 
                                      const Topology   *topol)
{
    auto itbonds = InteractionType::BONDS;
    if (pd->interactionPresent(itbonds))
    {
        auto fs = pd->findForcesConst(itbonds);
        
        if (topol->hasEntry(itbonds))
        {
            auto myatoms = topol->atoms();
            for(const auto &topentry : topol->entry(itbonds))
            {
                int ai = topentry->atomIndex(0);
                int aj = topentry->atomIndex(1);

                int reali = tellme_RealAtom(ai, myatoms);
                int realj = tellme_RealAtom(aj, myatoms);
                //here also myatoms[i].ffType(); but for ai aj
                auto name_ai = myatoms[ai].name();
                auto name_aj = myatoms[aj].name();
                addXmlResidueBond(residuePtr,
                                  nameIndex(name_ai, reali).c_str(),
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
        if (param_names[i])
        {
            add_xml_double(grandchild3, param_names[i], params[i]);
        }
    }
}

void OpenMMWriter::addXmlNonbonded(xmlNodePtr                       parent,
                                   const ForceField                *pd,
                                   const std::map<std::string, int> &ffTypeMap)
{
    auto fs     = pd->findForcesConst(InteractionType::VDW);
    auto fsCoul = pd->findForcesConst(InteractionType::ELECTROSTATICS);
    std::string nnn("nexcl");
    //int nrexcl  = std::max(my_atoi(fs.optionValue(nnn), "nrexclvdw"),
    //                     my_atoi(fsCoul.optionValue(nnn), "nrexclqq"));
    xmlNodePtr fsPtr = nullptr;
    // Custom non-bonded force is needed if we use Buckingham (not LJ)
    if (fs.potential() != Potential::LJ12_6)
    {
        fsPtr  = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMNONBONDEDFORCE));
        add_xml_double(fsPtr, "energy", 0.0);
        add_xml_int(fsPtr, "bondCutoff", 3);
        auto uafr = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::USEATTRIBUTEFROMRESIDUE));
        add_xml_char(uafr, exml_names(xmlEntryOpenMM::NAME), "charge");
        
        // This order is important!
        // Do not change the order of these parameters, otherwise the force field is not working
        auto grandchild2 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
        if (fs.potential() == Potential::WANG_BUCKINGHAM || fs.potential() == Potential::LJ14_7)
        {
            add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::NAME), "sigma");
        }
        else
        {
            add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::NAME), "rmin");
        }
        auto grandchild3 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
        add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::NAME), "epsilon");
        auto grandchild4 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
        add_xml_char(grandchild4, exml_names(xmlEntryOpenMM::NAME), "gamma");
        if (fs.potential() == Potential::GENERALIZED_BUCKINGHAM || fs.potential() == Potential::LJ14_7)
        {
            auto grandchild9 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
            add_xml_char(grandchild9, exml_names(xmlEntryOpenMM::NAME), "delta");
        }
        
        auto grandchild5 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
        add_xml_char(grandchild5, exml_names(xmlEntryOpenMM::NAME), "charge");
        auto grandchild6 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERPARTICLEPARAMETER));
        add_xml_char(grandchild6, exml_names(xmlEntryOpenMM::NAME), "zeta");
    }
    
    // This part is added to implement PME and LJPME, i.e. long-range interactions are approximated
    // by point charge and 12-6 Lennard-Jones (unless LJ_SR is the real potential).
    // For the different atomtypes it would be good to add optimized sigma and epsilon values
    // in order to calculate the vdW contribution after the cutoff
    auto ljPtr        = add_xml_child(parent, exml_names(xmlEntryOpenMM::NONBONDEDFORCE));
    add_xml_double(ljPtr, "coulomb14scale", 1.0); 
    add_xml_double(ljPtr, "lj14scale", 1.0); 
    add_xml_double(ljPtr, "energy", 0.0);
    add_xml_int(ljPtr, "bondCutoff", 3);
    auto uafr = add_xml_child(ljPtr, exml_names(xmlEntryOpenMM::USEATTRIBUTEFROMRESIDUE));
    add_xml_char(uafr, exml_names(xmlEntryOpenMM::NAME), "charge");
    
    // Add custom and/or regular nonbonded parameters.
    for(const auto &fft: ffTypeMap)
    {
        auto aType = pd->findParticleType(fft.first);
        for(int i = 1; i <= fft.second; i++)
        {
            std::string vdwtype("vdwtype");
            if (!aType->hasOption(vdwtype))
            {
                gmx_fatal(FARGS, "Incorrect force field. No option %s for %s.\n",
                          vdwtype.c_str(), fft.first.c_str());
            }
            auto vdwId = Identifier(aType->optionValue(vdwtype));
            auto param = fs.findParametersConst(Identifier(vdwId));
            
            std::string type1 = fft.first;
            if (addNumbersToAtomTypes_)
            {
                type1 = nameIndex(type1, i);
            }
            xmlNodePtr grandchild3 = nullptr;
            if (nullptr != fsPtr)
            {
                grandchild3 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::ATOM_RES));
                add_xml_char(grandchild3, exml_names(xmlEntryOpenMM::TYPE_RES), type1.c_str());  
            }
            double sigma = 0, epsilon = 0;
            switch (fs.potential())
            {
            case Potential::WANG_BUCKINGHAM:
                // TODO: optimize values
                sigma   = param[wbh_name[wbhSIGMA]].internalValue();
                epsilon = param[wbh_name[wbhEPSILON]].internalValue();
                for(size_t j = 0; j < param.size(); j++)
                {
                    if (Mutability::Dependent != param[wbh_name[j]].mutability())
                    {
                        add_xml_double(grandchild3, wbh_name[j], param[wbh_name[j]].internalValue());
                    }
                }
                break;
            case Potential::GENERALIZED_BUCKINGHAM:
                // TODO: optimize values
                sigma   = param[gbh_name[gbhRMIN]].internalValue()/std::pow(2,1.0/6.0);
                epsilon = param[gbh_name[gbhEPSILON]].internalValue();
                for(size_t j = 0; j < param.size(); j++)
                {
                    if (Mutability::Dependent != param[gbh_name[j]].mutability())
                    {
                        add_xml_double(grandchild3, gbh_name[j], param[gbh_name[j]].internalValue());
                    }
                }
                break;
            case Potential::LJ14_7:
                // TODO: optimize values
                sigma   = param[lj14_7_name[lj14_7SIGMA]].internalValue();
                epsilon = param[lj14_7_name[lj14_7EPSILON]].internalValue();
                for(size_t j = 0; j < param.size(); j++)
                {
                    if (Mutability::Dependent != param[lj14_7_name[j]].mutability())
                    {
                        add_xml_double(grandchild3, lj14_7_name[j], param[lj14_7_name[j]].internalValue());
                    }
                }
                break;		
            case Potential::LJ12_6:
                if (nullptr == fsPtr)
                {
                    // If we use "native" Lennard Jones we need to do this here:
                    auto ljchild = add_xml_child(ljPtr, exml_names(xmlEntryOpenMM::ATOM_RES));
                    add_xml_char(ljchild, exml_names(xmlEntryOpenMM::TYPE_RES), type1.c_str());
                    sigma   = param[lj12_6_name[lj12_6SIGMA]].internalValue();
                    epsilon = param[lj12_6_name[lj12_6EPSILON]].internalValue();
                    for(size_t j = 0; j < param.size(); j++)
                    {
                        if (Mutability::Dependent != param[lj12_6_name[j]].mutability())
                        {
                            add_xml_double(ljchild, lj12_6_name[j], param[lj12_6_name[j]].internalValue());
                        }
                    }
                }
                break;
            default:
                fprintf(stderr, "Unknown non-bonded force type %d %s\n", 
                        static_cast<int>(fs.potential()),
                        potentialToString(fs.potential()).c_str());
            }
            if (nullptr != fsPtr)
            {
                auto ztype = "zetatype";
                if (aType->hasOption(ztype))
                {
                    auto   ztp  = aType->optionValue(ztype);
                    double zeta = 0;
                    if (fsCoul.parameterExists(ztp))
                    {
                        zeta = fsCoul.findParameterTypeConst(ztp, "zeta").internalValue();
                    }
                    add_xml_double(grandchild3, "zeta", zeta); 
                }
                // Normal Lennard Jones is always needed
                {
                    auto ljchild = add_xml_child(ljPtr, exml_names(xmlEntryOpenMM::ATOM_RES));
                    add_xml_char(ljchild, exml_names(xmlEntryOpenMM::TYPE_RES), type1.c_str());
                    // Convert the sigma to get the correct position of the minimum.
                    double fac2 = std::pow(2.0, 1.0/6.0);
                    std::vector<double> se_param = { sigma/fac2, epsilon };
                    const char *se_name[] = { "sigma", "epsilon" };
                    
                    for(size_t j = 0; j < se_param.size(); j++)
                    {
                        add_xml_double(ljchild, se_name[j], se_param[j]);
                    }
                }
            }
        }  
    }    
}

void OpenMMWriter::addXmlPolarization(xmlNodePtr                        parent,
                                      const ForceField                 *pd,
                                      const std::map<std::string, int> &ffTypeMap,
                                      double                            epsr_fac)
{
    if (!pd->polarizable())
    {
        return;
    }
    // !!! Shell particle has to be type1, core particle has to be type2 !!!
    auto fs         = pd->findForcesConst(InteractionType::POLARIZATION);
    auto polPtr     = add_xml_child(parent, exml_names(xmlEntryOpenMM::DRUDEFORCE));
    for(const auto &fft: ffTypeMap)
    {
        auto aType = pd->findParticleType(fft.first);
        for(int i = 1; i <= fft.second; i++)
        {
            if (ActParticle::Atom == aType->apType())
            {
                if (aType->hasOption("poltype"))
                {
                    auto grandchild2 = add_xml_child(polPtr, exml_names(xmlEntryOpenMM::PARTICLE)); 

                    std::string type1 = aType->optionValue("poltype");
                    auto param        = fs.findParametersConst(Identifier(type1));
                    auto alpha        = fs.findParameterTypeConst(Identifier({type1}), pol_name[polALPHA]);
                    auto stp          = pd->findParticleType(type1);
                    
                    std::string type2 = fft.first;                    
                    if (addNumbersToAtomTypes_)
                    {
                        type1 = nameIndex(type1, i);
                        type2 = nameIndex(type2, i);
                    }
                    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE2), type2.c_str()); 
                    add_xml_char(grandchild2, exml_names(xmlEntryOpenMM::TYPE1), type1.c_str());
                        
                    add_xml_double(grandchild2, "polarizability", alpha.internalValue());
                    add_xml_double(grandchild2, exml_names(xmlEntryOpenMM::CHARGE_RES), stp->charge()*epsr_fac);
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
        switch(fs.second.potential())
        {
        case Potential::HARMONIC_BONDS:
            {
                fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::HARMONICBONDFORCE));
            }
            break;
        case Potential::MORSE_BONDS:
            {
                fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMBONDFORCE));
                // The Morse bonds potential is written as a string here:
                auto energy = gmx::formatString("(%s*((1 - exp(-%s*(r-%s)))^2-1)+%s);", 
                                                morse_name[morseDE], morse_name[morseBETA],
                                                morse_name[morseLENGTH], morse_name[morseD0]);
                add_xml_char(fsPtr, "energy", energy.c_str());
                // Specify the per bond parameters
                for(int i = 0; i < morseNR; i++)
                {
                    auto grandchild0 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
                    add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), morse_name[i]);
                }
            }
            break;
        case Potential::CUBIC_BONDS:
            {
                fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMBONDFORCE));
                // The cubic bonds potential is written as a string here:
                add_xml_char(fsPtr, "energy", "(kb*(r-bondlength)^2 * (rmax-r) - D_e)*step(2*rmax+bondlength-3*r) + step(3*r-2*rmax-bondlength)*(-D_e - (4*kb/27)*(bondlength-rmax)^2);");
                // Specify the per bond parameters
                for(int i = 0; i < cubicNR; i++)
                {
                    auto grandchild0 = add_xml_child(fsPtr, exml_names(xmlEntryOpenMM::PERBONDPARAMETER));
                    add_xml_char(grandchild0, exml_names(xmlEntryOpenMM::NAME), cubic_name[i]);
                }
            }
            break;
        case Potential::HARMONIC_ANGLES:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::HARMONICANGLEFORCE));
            add_xml_double(fsPtr, "energy", 0);
            break;
        case Potential::LINEAR_ANGLES:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMANGLEFORCE));
            add_xml_double(fsPtr, "energy", 0);
            break;
        case Potential::UREY_BRADLEY_ANGLES:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::CUSTOMANGLEFORCE));
            add_xml_double(fsPtr, "energy", 0);
            break;
        case Potential::FOURIER_DIHEDRALS:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::RBTORSIONFORCE));
            add_xml_double(fsPtr, "energy", 0);
            break;
        case Potential::PROPER_DIHEDRALS:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::PERIODICTORSIONFORCE));
            add_xml_double(fsPtr, "energy", 0);
            break;
        case Potential::HARMONIC_DIHEDRALS:
            fsPtr = add_xml_child(parent, exml_names(xmlEntryOpenMM::PERIODICTORSIONFORCE));
            add_xml_double(fsPtr, "energy", 0);
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
            if (entry->id().atoms().empty())
            {
                if (debug)
                {
                    fprintf(debug, "Skipping interaction %s\n",
                           interactionTypeToString(fs.first).c_str());
                }
                continue;
            }
            const auto &bondId = entry->id().id();
            const auto &atoms  = entry->id().atoms();
            if (ClassUsed.end() == ClassUsed.find(bondId))
            {
                ClassUsed.insert(bondId);
                switch(fs.second.potential())
                {
                case Potential::HARMONIC_BONDS:
                    {
                        const char *omm_bonds[] = { "k", "length", nullptr };
                        addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::BOND_RES,
                                   atoms, omm_bonds, entry->params());
                    }
                    break;
                case Potential::CUBIC_BONDS:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::BOND_RES,
                               atoms, cubic_name, entry->params());
                    break;
                case Potential::MORSE_BONDS:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::BOND_RES,
                               atoms, morse_name, entry->params());
                    break;
                case Potential::HARMONIC_ANGLES:
                    {
                        const char *omm_angles[] = { "k", "angle" };
                        addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::ANGLE_CLASS,
                                   atoms, omm_angles, entry->params());
                    }
                    break;
                case Potential::UREY_BRADLEY_ANGLES:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::ANGLE_CLASS,
                               atoms, ub_name, entry->params());
                    break;
                    
                case Potential::LINEAR_ANGLES:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::ANGLE_CLASS,
                               atoms, linang_name, entry->params());
                    break;
                    
                case Potential::FOURIER_DIHEDRALS:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::RBTORSIONFORCE,
                               atoms, fdih_name, entry->params());
                    break;
                    
                case Potential::PROPER_DIHEDRALS:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::PERIODICTORSIONFORCE,
                               atoms, pdih_name, entry->params());
                    break;
                    
                case Potential::HARMONIC_DIHEDRALS:
                    addXmlBond(xmlMap_[fs.first], xmlEntryOpenMM::IMPROPER,
                               atoms, idih_name, entry->params());
                    break;
                    
                case Potential::COULOMB_POINT:
                case Potential::COULOMB_GAUSSIAN:
                case Potential::COULOMB_SLATER:
                case Potential::LJ12_6:
                case Potential::LJ14_7:
                case Potential::LJ8_6:
                case Potential::WANG_BUCKINGHAM:
                case Potential::GENERALIZED_BUCKINGHAM:
                case Potential::POLARIZATION:
                case Potential::VSITE3OUTS:
                case Potential::VSITE3OUT:
                case Potential::VSITE3S:
                case Potential::VSITE3:
                case Potential::VSITE2:
                case Potential::VSITE1:
                    break;
                default:
                    fprintf(stderr, "Warning: no OpenMM support for %s is present or implemented.\n",
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
    // Compute epsilon r if relevant
    auto coul       = pd->findForcesConst(InteractionType::ELECTROSTATICS);
    double epsr_fac = 1;
    std::string epsr("epsilonr");
    if (coul.optionExists(epsr))
    {
        auto epsval   = coul.optionValue(epsr);
        auto epsilonr = my_atof(epsval.c_str(), "EpsilonR");
        if (epsilonr <= 0)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Incorrect epsilonr '%s'", epsval.c_str()).c_str()));
        }
        epsr_fac = std::sqrt(1.0/epsilonr);
        if (epsilonr != 1)
        {
            printf("Scaling charges to accomodate for epsilonr = %g\n", epsilonr);
        }
    }

    // For looking up InChi codes
    AlexandriaMols amols;    
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
        if (!fragmentHandler)
        {
            printf("No complete information for %s, skipping conversion to OpenMM.\n", actmol.getMolname().c_str());
            continue;
        }
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
            auto myatoms       = topologies[fff]->atoms();
                    
            // Check whether we have to terminate the residue by defining bonds
            residuePtr = add_xml_child(xmlResiduePtr, exml_names(xmlEntryOpenMM::RESIDUE));
            auto amol  = amols.find(fragIds[fff]);
            if (nullptr != amol)
            {
                add_xml_char(residuePtr, exml_names(xmlEntryOpenMM::NAME), amol->iupac.c_str());
            }
            else
            {
                add_xml_char(residuePtr, exml_names(xmlEntryOpenMM::NAME), fragIds[fff].c_str());
            }
            InchiUsed.insert(fragIds[fff]);
            // If we need to add numbers to types to distinguish types within a compound,
            // then we need to store what we have used. The map is from force field type
            // to a number index. It restarts from 1 for every molecule, then if we add a
            // new copy of the same force field type we increase the number.
            std::map<std::string, int> fftypeLocalMap;
            // We need to keep track of the atom names within the molecule as well
            std::vector<std::string> inames;
            
            for (size_t i = 0; i < myatoms.size(); i++)
            {
                int reali = tellme_RealAtom(i, myatoms);
                auto iname = nameIndex(myatoms[i].name(), reali);
                inames.push_back(iname);
            }
            for (size_t i = 0; i < myatoms.size(); i++)
            {
                // First do atom type stuff
                auto ffType            = myatoms[i].ffType();
                // Local index, default 1
                int localIndex         = 1;
                auto ffLocalPtr        = fftypeLocalMap.find(ffType);
                if (fftypeLocalMap.end() != ffLocalPtr)
                {
                    // Increase localPtr if needed and store it in variable
                    if (addNumbersToAtomTypes_)
                    {
                        ffLocalPtr->second++;
                        localIndex = ffLocalPtr->second;
                    }
                }
                else
                {
                    fftypeLocalMap.insert({ffType, localIndex});
                }
                // Make OpenMM type
                std::string openMMtype = ffType;
                if (addNumbersToAtomTypes_)
                {
                    openMMtype += gmx::formatString("_%d", localIndex);
                    if (debug)
                    {
                        fprintf(debug, "Adding openMMtype %s counter %d\n", openMMtype.c_str(), localIndex);
                    }
                }
                // Now check whether this type is in the global map
                bool addGlobalAtomType = false;
                auto ffGlobalPtr = fftypeGlobalMap.find(ffType);
                if (fftypeGlobalMap.end() == ffGlobalPtr)
                {
                    // New global atomtype
                    fftypeGlobalMap.insert({ ffType, 1});
                    ffGlobalPtr = fftypeGlobalMap.find(ffType);

                    addGlobalAtomType = true;
                }
                else if (addNumbersToAtomTypes_ && localIndex > ffGlobalPtr->second)
                {
                    ffGlobalPtr->second = localIndex;
                    addGlobalAtomType = true;
                }
                if (addGlobalAtomType)
                {
                    auto aType = pd->findParticleType(ffType);
                    // Make OpenMM class, that does not need the numbering.
                    std::string openMMclass = ffType;
                    std::string btype("bondtype");
                    if (aType->hasOption(btype))
                    {
                        openMMclass = aType->optionValue(btype);
                    }
                    auto newAtype = add_xml_child(xmlAtomTypePtr, exml_names(xmlEntryOpenMM::TYPE));
                    add_xml_char(newAtype, exml_names(xmlEntryOpenMM::NAME), openMMtype.c_str());
                    add_xml_char(newAtype, exml_names(xmlEntryOpenMM::CLASS), openMMclass.c_str());
                    addXmlElemMass(newAtype, aType);
                }
                int reali = tellme_RealAtom(i, myatoms);

                if (myatoms[i].pType() == ActParticle::Atom || myatoms[i].pType() == ActParticle::Vsite)
                {
                    auto baby = add_xml_child(residuePtr, exml_names(xmlEntryOpenMM::ATOM_RES));
                    auto iname = nameIndex(myatoms[i].name(), reali);
                    add_xml_char(baby, exml_names(xmlEntryOpenMM::NAME), iname.c_str());
                    
                    add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), openMMtype.c_str());
                    add_xml_double(baby, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[i].charge()*epsr_fac);
                    
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
                        add_xml_double(baby, exml_names(xmlEntryOpenMM::CHARGE_RES), myatoms[shell].charge()*epsr_fac);
                        ishell += 1;
                    }
                    if (myatoms[i].pType() == ActParticle::Vsite)
                    {
                        auto baby = add_xml_child(residuePtr, exml_names(xmlEntryOpenMM::VSITE_RES));
                        // Could be different vsite types, remember to implement those.
                        std::map<InteractionType, const char *> itypes = { 
                            { InteractionType::VSITE1,     "average2"   },
                            { InteractionType::VSITE2,     "average2"   },
                            { InteractionType::VSITE3OUT,  "outOfPlane" },
                            { InteractionType::VSITE3OUTS, "outOfPlane" },
                            { InteractionType::VSITE3,     "average3"   },
                            { InteractionType::VSITE3S,    "average3"   },
                            { InteractionType::VSITE3FAD,  "average3"   }
                        };
                        for(const auto &itp : itypes)
                        {
                            if (!topologies[fff]->hasEntry(itp.first))
                            {
                                continue;
                            }
                            for(const auto &ee : topologies[fff]->entry(itp.first))
                            {
                                auto indices = ee->atomIndices();
                                // Check whether the present vsite is the last in the
                                // atom indices. Let's ignore multiple vsite constructors
                                // for the same vsite for now.
                                if (static_cast<int>(i) == indices[indices.size()-1])
                                {
                                    add_xml_char(baby, exml_names(xmlEntryOpenMM::TYPE_RES), itp.second);
                                    add_xml_char(baby, exml_names(xmlEntryOpenMM::SITENAME), iname.c_str());
                                    
                                    // TODO get data from topology instead of making stuff up.
                                    
                                    std::vector<int> cores = myatoms[i].cores();
                                    if (itp.first ==  InteractionType::VSITE3OUTS ||
                                        itp.first ==  InteractionType::VSITE3OUT)
                                    {
                                        // We have to swap cores 0 and 1
                                        std::swap(cores[0], cores[1]);
                                    }
                                    int ppp = 1;
                                    for(const auto &parent : cores)
                                    {
                                        auto an = gmx::formatString("atomName%d", ppp);
                                        if (static_cast<size_t>(parent) < inames.size())
                                        {
                                            add_xml_char(baby, an.c_str(), inames[parent].c_str());
                                            // Add an (artificial) bond such that OpenMM will generate
                                            // exclusions for this. See
                                            // https://github.com/openmm/openmm/issues/811
                                            addXmlResidueBond(residuePtr,
                                                              inames[i].c_str(), inames[parent].c_str());
                                        }
                                        else
                                        {
                                            add_xml_char(baby, an.c_str(), myatoms[parent].name().c_str());
                                        }
                                        ppp += 1;
                                    }
                                    switch (itp.first)
                                    {
                                    case InteractionType::VSITE2:
                                        {
                                            auto p = ee->params()[0];
                                            add_xml_double(baby, "weight1", 1-p);
                                            add_xml_double(baby, "weight2", p);
                                        }
                                        break;
                                    case InteractionType::VSITE3S:
                                        {
                                            auto p = ee->params()[0];
                                            add_xml_double(baby, "weight1", p);
                                            add_xml_double(baby, "weight2", 1-2*p);
                                            add_xml_double(baby, "weight3", p);
                                        }
                                        break;
                                    case InteractionType::VSITE3:
                                        {
                                            auto p = ee->params();
                                            add_xml_double(baby, "weight1", p[0]);
                                            add_xml_double(baby, "weight2", 1-p[0]-p[1]);
                                            add_xml_double(baby, "weight3", p[1]);
                                        }
                                        break;
                                    case InteractionType::VSITE3OUT:
                                        {
                                            auto p = ee->params();
                                            add_xml_double(baby, "weight1", p[0]);
                                            add_xml_double(baby, "weight2", p[1]);
                                            // To make the two vsites different, we assume they are next to each other
                                            // in the row of atoms. In which order the vsites have +/- does not matter.
                                            int sign = 2*(i % 2)-1;
                                            add_xml_double(baby, "weight3", sign*p[2]);
                                        }
                                        break;
                                    case InteractionType::VSITE3OUTS:
                                        {
                                            auto p = ee->params();
                                            add_xml_double(baby, "weight1", p[0]);
                                            add_xml_double(baby, "weight2", p[0]);
                                            // See above
                                            int sign = 2*(i % 2)-1;
                                            add_xml_double(baby, "weight3", sign*p[1]);
                                        }
                                        break;
                                    default:
                                        GMX_THROW(gmx::InternalError(gmx::formatString("InteractionType %s not supported yet.",
                                                                                       interactionTypeToString(itp.first).c_str()).c_str()));
                                    }
                                }
                            }
                        }
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
    addXmlPolarization(parent, pd, fftypeGlobalMap, epsr_fac);
}

void OpenMMWriter::writeXml(const std::string         &fileName,
                            const ForceField          *pd,
                            const std::vector<ACTMol> &actmols)
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

    // Do not compress the OpenMM file.
    xmlSetDocCompressMode(doc, 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fileName.c_str(), doc, "ISO-8859-1", 2) == 0)
    {
        gmx_fatal(FARGS, "Saving file %s", fileName.c_str());
    }
    xmlFreeDoc(doc);
}

static void writeCombinationRules(FILE *fp,
                                  const std::map<const std::string, CombRule> cmap)
{
    fprintf(fp, "combinationRule =");
    for(const auto &cm : cmap)
    {
        fprintf(fp, " %s %s", cm.first.c_str(), combinationRuleName(cm.second).c_str());
    }
    fprintf(fp, "\n");
}

void OpenMMWriter::writeDat(const std::string &fileName,
                            const ForceField  *pd)
{
    std::map<InteractionType, std::vector<std::pair<std::string, std::string>>> act2omm = {
        { InteractionType::ELECTROSTATICS, 
          { { "chargetype", "chargeDistribution" },
            { "epsilonr", "dielectricConstant" },
            { "nexcl", "nexclqq" } } },
        { InteractionType::VDW,
          { {  "nexcl", "nexclvdw" } } }
    };

    FILE *fp = gmx_ffopen(fileName.c_str(), "w");
    fprintf(fp, "rigidWater  = False\n");
    fprintf(fp, "constraints = None\n");
    for(const auto &a2o : act2omm)
    {
        if (pd->interactionPresent(a2o.first))
        {
            auto fs = pd->findForcesConst(a2o.first);
            for(const auto &opt : a2o.second)
            {
                if (fs.optionExists(opt.first))
                {
                    fprintf(fp, "%s = %s\n", opt.second.c_str(),
                            fs.optionValue(opt.first).c_str());
                }
            }
            if (InteractionType::VDW == a2o.first)
            {
                fprintf(fp, "vanderwaals = %s\n", 
                        potentialToString(fs.potential()).c_str());
                writeCombinationRules(fp, getCombinationRule(fs));
            }
        }
    }
    
    gmx_ffclose(fp);
}

void writeOpenMM(const std::string         &fileName,
                 const std::string         &simParams,
                 const ForceField          *pd,
                 const std::vector<ACTMol> &actmols,
                 double                     mDrude,
                 bool                       addNumbersToAtoms)
{
    rmapyyyOpenMM.clear();
    
    OpenMMWriter writer(mDrude, addNumbersToAtoms);
    writer.writeXml(fileName, pd, actmols);
    writer.writeDat(simParams, pd);
}

} // namespace alexandria
