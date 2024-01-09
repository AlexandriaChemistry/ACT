/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include "babel_io.h"

#include "config.h"

#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>

#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "act/alexandria/actmol.h"
#include "act/alexandria/atype_mapping.h"
#include "act/basics/allmols.h"
#include "act/forcefield/forcefield.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/multipole_names.h"
#include "act/utility/stringutil.h"
#include "act/utility/units.h"

// Include Open Babel classes for OBMol and OBConversion
#if HAVE_LIBOPENBABEL3
// Hack to make this compile!
#undef ANGSTROM
#ifdef HAVE_SYS_TIME_H
#define KOKO HAVE_SYS_TIME_H
#undef HAVE_SYS_TIME_H
#endif
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wdeprecated-copy-with-user-provided-dtor"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-dtor"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#include <openbabel/atom.h>
#include <openbabel/babelconfig.h>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/data_utilities.h>
#include <openbabel/elements.h>
#include <openbabel/forcefield.h>
#include <openbabel/generic.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/residue.h>
#include <openbabel/typer.h>
#include <openbabel/math/vector3.h>
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

#ifdef KOKO
#ifndef HAVE_SYS_TIME_H
#define HAVE_SYS_TIME_H KOKO
#endif
#undef KOKO
#endif

using namespace alexandria;

enum einformat{
    einfGaussian    = 0,
    einfNotGaussian = 1,
    einfNR
};

BabelFile::BabelFile(BabelFileType      ftype,
                     const std::string &ext,
                     const std::string &InFormat)
    :
      ftype_(ftype),
      ext_(ext),
      InFormat_(InFormat)
{}

BabelFiles::BabelFiles ()
{
    bfiles_.push_back(BabelFile(ebftPDB,  ".pdb",  "pdb"));
    bfiles_.push_back(BabelFile(ebftXYZ,  ".xyz",  "xyz"));
    bfiles_.push_back(BabelFile(ebftSDF,  ".sdf",  "sdf"));
    bfiles_.push_back(BabelFile(ebftMOL,  ".mol",  "mol"));
    bfiles_.push_back(BabelFile(ebftMOL2, ".mol2", "mol2"));
    bfiles_.push_back(BabelFile(ebftG09,  ".log",  "g03"));
}

BabelFileIterator BabelFiles::findBabelFile(const std::string &fn)
{
    auto              extension = fn.substr(fn.find_last_of("."));
    BabelFileIterator fb        = bfiles_.begin(), fe = bfiles_.end();
    return std::find_if(fb, fe, [extension](const BabelFile &bf)
                        {
                            return (extension == bf.ext());
                        });
}

static bool isGzipFile(const std::string &fileName,
                       std::string       *strippedFileName)
{
    auto gzPosition = fileName.find(".gz");
    if (gzPosition == fileName.size()-3)
    {
        strippedFileName->assign(fileName.substr(0, gzPosition));
        return true;
    }
    return false;
}

static void checkBondOrders(alexandria::MolProp *mpt)
{
    // We cannot be sure that this is the right one
    auto exper = mpt->LastExperiment();
    std::vector<std::string> atomName;
    for(auto &ca : exper->calcAtomConst())
    {
        atomName.push_back(ca.getObtype());
    }
    std::map<std::string, std::string> my_pairs =
        { { "no", "on" }, { "cm", "om" }, { "p5", "om" }, { "s6", "om" }, { "s4", "om" }, { "s3", "om" }, { "py", "om" }, { "cz", "n" } };
    for(const auto &my_pair : my_pairs)
    {
        for(size_t i = 0; i < atomName.size(); i++)
        {
            if (atomName[i] == my_pair.first)
            {
                std::vector<size_t> bondIndex;
                auto allBonds = mpt->bonds();
                for(size_t b = 0; b < (*allBonds).size(); b++)
                {
                    size_t ai = (*allBonds)[b].aI();
                    size_t aj = (*allBonds)[b].aJ();
                    if ((ai == i && atomName[aj] == my_pair.second) ||
                        (aj == i && atomName[ai] == my_pair.second))
                    {
                        bondIndex.push_back(b);
                    }
                }
                // The first atom has to have at least two bonds to the
                // second atom in the map.
                if (bondIndex.size() >= 2)
                {
                    for (size_t j = 0; j < bondIndex.size(); j++)
                    {
                        (*allBonds)[bondIndex[j]].setBondOrder(0, 1.5);
                    }
                }
            }
        }
    }
}

static bool getBondsFromOpenBabel(OpenBabel::OBMol    *mol,
                                  alexandria::MolProp *mpt,
                                  const char          *inputFile,
                                  bool                 changeAromaticBondOrders)
{
    // Bonds
    auto OBbi = mol->BeginBonds();
    if (OBbi != mol->EndBonds())
    {
        for (auto OBb = mol->BeginBond(OBbi); (nullptr != OBb); OBb = mol->NextBond(OBbi))
        {
            double bo = OBb->GetBondOrder();
            if (changeAromaticBondOrders && OBb->IsAromatic())
            {
                bo = 1.5;
            }
            alexandria::Bond ab(OBb->GetBeginAtom()->GetIndex(),
                                OBb->GetEndAtom()->GetIndex(),
                                bo);
            mpt->AddBond(ab);
        }
        if (changeAromaticBondOrders)
        {
            checkBondOrders(mpt);
        }
        return true;
    }
    else
    {
        fprintf(stderr, "No bond is found for %s\n", inputFile);
        return false;
    }
}

static bool addInchiToFragments(const AlexandriaMols    &amols,
                                OpenBabel::OBConversion *conv,
                                OpenBabel::OBMol        *mol,
                                std::vector<Fragment>   *fragptr)
{
    conv->SetOutFormat("inchi");
    std::string alanine("ALA");
    size_t fff = 0;
    for(auto fptr = fragptr->begin(); fptr < fragptr->end(); fptr++)
    {
        // Copy input molecule
        OpenBabel::OBMol      fmol;
        // Atoms in the fragment
        auto                  fatoms = fptr->atoms();
        std::map<int, int>    renumber;
        int                   count = 1;
        FOR_ATOMS_OF_MOL (atom, *mol)
        {
            int idx = atom->GetIdx();
            if (std::find(fatoms.begin(), fatoms.end(), idx-1) != fatoms.end())
            {
                OpenBabel::OBAtom    newatom;
                newatom.Duplicate(&(*atom));
                // Make a copy of the residue information
                OpenBabel::OBResidue residue = *(atom->GetResidue());
                residue.SetNum(fff+1);
                newatom.SetIdx(count++);
                newatom.SetResidue(&residue);
                if (!fmol.AddAtom(newatom))
                {
                    fprintf(stderr, "Could not add atom %d to fmol\n", atom->GetIdx());
                }
                renumber.insert({idx, newatom.GetIdx()});
            }
        }
        fff += 1;
        auto OBbi = mol->BeginBonds();
        if (OBbi != mol->EndBonds())
        {
            for (auto OBb = mol->BeginBond(OBbi); (nullptr != OBb); OBb = mol->NextBond(OBbi))
            {
                int ai = OBb->GetBeginAtom()->GetIdx();
                int aj = OBb->GetEndAtom()->GetIdx();
                int bo = OBb->GetBondOrder();
                if (std::find(fatoms.begin(), fatoms.end(), ai-1) != fatoms.end() &&
                    std::find(fatoms.begin(), fatoms.end(), aj-1) != fatoms.end())
                {
                    if (!fmol.AddBond(renumber[ai], renumber[aj], bo))
                    {
                        fprintf(stderr, "Could not add bond to fmol\n");
                    }
                }
            }
        }
        auto inchi = conv->WriteString(&fmol, true);
        auto amol  = amols.find(inchi);
        if (nullptr != amol)
        {
            fptr->setId(amol->iupac);
            fptr->setCharge(amol->charge);
            fptr->setMass(amol->mass);
            fptr->setFormula(amol->formula);
        }
        else
        {
            fptr->setId(inchi);
        }
    }
    return true;
}

static bool babel2ACT(const ForceField                         *pd,
                      const std::map<std::string, std::string> &g2a,
                      const AlexandriaMols                     &amols,
                      OpenBabel::OBMol                         *mol,
                      alexandria::MolProp                      *mpt,
                      const char                               *molnm,
                      const char                               *iupac,
                      const char                               *conformation,
                      std::string                              *method,
                      std::string                              *basisset,
                      int                                       maxPotential,
                      int                                       nsymm,
                      const char                               *jobType,
                      double                                   *qtot,
                      bool                                      addHydrogen,
                      const char                               *g09,
                      einformat                                 inputformat)
{
    std::string                formula;
    std::string                attr;
    std::string                value;
    const char                *reference   = "Spoel2022a";
    const char                *actmol       = "AMM";
    const char                *myprogram   = "ACT2022";

    /* Variables to read a Gaussian log file */
    char                      *g09ptr;
    alexandria::JobType jobtype = alexandria::string2jobType(jobType);
    auto conv = new OpenBabel::OBConversion(&std::cin, &std::cout);

    // Chemical Categories
    if (conv->SetOutFormat("fpt"))
    {
        std::set<std::string> excluded_categories = {
            ">", "C_ONS_bond", "Rotatable_bond", "Conjugated_double_bond", "Conjugated_triple_bond",
            "Chiral_center_specified", "Cis_double_bond", "Bridged_rings", "Conjugated_tripple_bond",
            "Trans_double_bond"
        };
        
        conv->AddOption("f", OpenBabel::OBConversion::OUTOPTIONS, "FP4");
        conv->AddOption("s");
        conv->Convert();
        
        // We need a copy here because WriteString removes the H.
        OpenBabel::OBMol  mol_copy = *mol;
        std::vector<std::string> categories = gmx::splitString(conv->WriteString(&mol_copy, false));
        for (const auto &category : categories)
        {
            if (excluded_categories.find(category) == excluded_categories.end())
            {
                std::string dup = category;
                std::replace_if(dup.begin(), dup.end(), [](const char c) {return c == '_'; }, ' ');
                mpt->AddCategory(dup);
            }
        }
    }

    // Basis Set
    std::string basis;
    auto OBpd = (OpenBabel::OBPairData *)mol->GetData("basis");
    if (!basisset->empty())
    {
        basis.assign(*basisset);
    }
    else if (nullptr != OBpd)
    {
        basis = OBpd->GetValue();
        size_t p = basis.find(" (5D, 7F)");
        if (p != basis.npos)
        {
            basis.erase(p, basis.npos);
        }
        basisset->assign(basis);
    }
    else
    {
        basis.clear();
    }

    // QM Program
    std::string program(myprogram);
    OBpd = (OpenBabel::OBPairData *)mol->GetData("program");
    if (nullptr != OBpd)
    {
        program.assign(OBpd->GetValue());
    }

    // Method
    OBpd = (OpenBabel::OBPairData *)mol->GetData("method");
    if (nullptr != OBpd)
    {
        method->assign(OBpd->GetValue());
    }
    g09ptr = (char *) strrchr(g09, '/');
    if (nullptr == g09ptr)
    {
        g09ptr = (char *)g09;
    }
    else
    {
        g09ptr++;
        if (strlen(g09ptr) == 0)
        {
            g09ptr = (char *)g09;
        }
    }

    alexandria::Experiment exp(program, *method, basis, reference, conformation, g09ptr, jobtype);
    // We don't just set this here, since the user may override the value
    // However, it seems that OB does not extract this correctly from
    // the input, unless it is a Gaussian log file.
    if (*qtot != 0)
    {
        double qtest = mol->GetTotalCharge();
        if (qtest != *qtot)
        {
            fprintf(stderr,"WARNING: OpenBabel found a total charge of %g, user specified %g. File %s.\n",
                    qtest, *qtot, g09);
        }
        mol->SetTotalCharge(*qtot);
    }
    else
    {
        *qtot = mol->GetTotalCharge();
    }
    mpt->AddExperiment(exp);
    if (nullptr != molnm)
    {
        mpt->SetMolname(molnm);
    }
    else
    {
        mpt->SetMolname(actmol);
    }

    if (nullptr != iupac)
    {
        mpt->SetIupac(iupac);
    }
    else
    {
        mpt->SetIupac(actmol);
    }
    
    // Units for conversion
    std::string energyUnit("kcal/mol");
    std::string entropyUnit("cal/mol K");
    std::string qm_type("electronic");

    // Thermochemistry
    if (inputformat == einfGaussian)
    {
        double              temperature = 0, DeltaHf0 = 0, DeltaHfT = 0;
        double              DeltaGfT    = 0, DeltaSfT = 0, S0T      = 0;
        double              CVT         = 0, CPT      = 0, ZPE      = 0;
        std::vector<double> Scomponents;
        if (extract_thermochemistry(*mol, false, &nsymm,
                                    0, 0.0,
                                    &temperature,
                                    &DeltaHf0,
                                    &DeltaHfT,
                                    &DeltaGfT,
                                    &DeltaSfT,
                                    &S0T,
                                    &CVT,
                                    &CPT,
                                    Scomponents,
                                    &ZPE))
        {
            {
                auto mpo = MolPropObservable::DHFORM;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, energyUnit,
                                                           0, ePhase::GAS,
                                                           alexandria::convertToGromacs(DeltaHf0, energyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::DHFORM;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, energyUnit,
                                                           temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(DeltaHfT, energyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::DGFORM;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, energyUnit,
                                                           temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(DeltaGfT, energyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::DSFORM;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, entropyUnit,
                                                           temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(DeltaSfT, entropyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::ENTROPY;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, entropyUnit,
                                                           temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(S0T, entropyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::CP;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, entropyUnit,
                                                           temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(CPT, entropyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            std::vector<MolPropObservable> mpos = { 
                MolPropObservable::STRANS, 
                MolPropObservable::SROT,
                MolPropObservable::SVIB
            };
            for (size_t i = 0; (i < mpos.size()); i++)
            {
                auto me = new alexandria::MolecularEnergy(mpos[i], qm_type, entropyUnit,
                                                          temperature, ePhase::GAS,
                                                          alexandria::convertToGromacs(Scomponents[i], entropyUnit), 0);
                mpt->LastExperiment()->addProperty(mpos[i], me);
            }
            {
                auto mpo = MolPropObservable::ZPE;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, energyUnit,
                                                           0, ePhase::GAS,
                                                           alexandria::convertToGromacs(ZPE, energyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
        }
    }

    // HF Eenergy
    { 
        auto mpo = MolPropObservable::HF;
        auto me  = new alexandria::MolecularEnergy(mpo, qm_type, energyUnit,
                                                   0, ePhase::GAS, 
                                                   alexandria::convertToGromacs(mol->GetEnergy(), energyUnit), 0);
        mpt->LastExperiment()->addProperty(mpo, me);
    }

    if (addHydrogen)
    {
        mol->AddHydrogens();
    }
    // Frequencies
    auto vibtype = OpenBabel::OBGenericDataType::VibrationData;
    if (mol->HasData(vibtype))
    {
        auto vibdata = static_cast<OpenBabel::OBVibrationData *>(mol->GetData(vibtype));
        auto freq    = vibdata->GetFrequencies();
        if (!freq.empty())
        {
            auto mpo = MolPropObservable::FREQUENCY;
            auto hf  = new Harmonics(mpo_unit2(mpo), 0, mpo);
            for (const auto &f : freq)
            {
                hf->addValue(f);
            }
            mpt->LastExperiment()->addProperty(mpo, hf);
        }
        auto inten   = vibdata->GetIntensities();
        if (!inten.empty())
        {
            auto mpo = MolPropObservable::INTENSITY;
            auto hf  = new Harmonics(mpo_unit2(mpo), 0, mpo);
            for (const auto &f : inten)
            {
                hf->addValue(f);
            }
            mpt->LastExperiment()->addProperty(mpo, hf);
        }
    }

    // Atoms
    const std::string forcefield("alexandria");
    auto *ff = OpenBabel::OBForceField::FindForceField(forcefield);
    std::vector<int> atomIndices;
    int oldresnum = -1;
    if (ff && (ff->Setup(*mol)))
    {
        ff->GetAtomTypes(*mol);
        FOR_ATOMS_OF_MOL (atom, *mol)
        {
            // For our molecular fragments the indices start at 0
            atomIndices.push_back(atom->GetIdx()-1);
            OpenBabel::OBPairData *type = (OpenBabel::OBPairData*) atom->GetData("FFAtomType");
            if (nullptr == type)
            {
                fprintf(stderr, "Error: OpenBabel cannot find the %s atomtype for atom %d\n",
                        forcefield.c_str(), static_cast<int>(atom->GetIdx()));
                return false;
            }
            if (nullptr != debug)
            {
                fprintf(debug, "atom %d gafftype %s OBtype %s\n", atom->GetIdx(), type->GetValue().c_str(), atom->GetType());
            }

            alexandria::CalcAtom ca(OpenBabel::OBElements::GetSymbol(atom->GetAtomicNum()),
                                    type->GetValue(), atom->GetIdx());

            ca.setCoordUnit("Angstrom");
            ca.setCoords(atom->x(), atom->y(), atom->z());
            auto myres = atom->GetResidue();
            // Workaround for incorrect residue numbers coming from babel.
            // For instance for OHH one typically gets 1 0 0 as residue numbers. 
            if (myres->GetNum() > oldresnum)
            {
                oldresnum = myres->GetNum();
            }
            ca.setResidue(myres->GetName(), oldresnum);
            ca.SetChain(myres->GetChainNum(), myres->GetChain());
            if (inputformat == einfGaussian)
            {
                for (const auto &cs : qTypes())
                {
                    if (cs.first == qType::Mulliken ||
                        cs.first == qType::ESP ||
                        cs.first == qType::Hirshfeld ||
                        cs.first == qType::CM5)
                    {
                        std::string qstr = cs.second;
                        qstr.append(" charges");
                        OBpd = (OpenBabel::OBPairData *) mol->GetData(qstr.c_str());
                        if (nullptr != OBpd)
                        {
                            auto OBpc = (OpenBabel::OBPcharge *) mol->GetData(qstr.c_str());
                            if (OBpc && !OBpc->GetPartialCharge().empty())
                            {
                                ca.AddCharge(stringToQtype(cs.second),
                                             OBpc->GetPartialCharge()[atom->GetIdx()-1]);
                            }
                            else
                            {
                                fprintf(stderr, "Inconsistency reading %s from %s", qstr.c_str(), g09);
                                return false;
                                
                            }
                        }
                    }
                }
            }
            mpt->LastExperiment()->AddAtom(ca);
        }
    }
    else
    {
        fprintf(stderr, "Error: OpenBabel cannot read the '%s' force field\n",
                forcefield.c_str());
        return false;
    }
    // Bonds
    getBondsFromOpenBabel(mol, mpt, g09, forcefield.compare("alexandria") == 0);

    // Convert atom types to Alexandria
    if (!g2a.empty())
    {
        if (!renameAtomTypes(mpt, g2a))
        {
            return false;
        }
    }
    // Fragment information will be generated
    mpt->generateFragments(pd, *qtot);
    addInchiToFragments(amols, conv, mol, mpt->fragmentPtr());
    
    // Dipole
    auto my_dipole = mol->GetData("Dipole Moment");
    if (nullptr != my_dipole)
    {
        auto dipole = (OpenBabel::OBVectorData *) my_dipole;
        OpenBabel::vector3 v3  = dipole->GetData();
        auto               mpo = MolPropObservable::DIPOLE;
        auto               dp  = new alexandria::MolecularMultipole(qm_type, "D", 0.0, mpo);
        dp->setValue(multipoleName({ XX }), v3.GetX());
        dp->setValue(multipoleName({ YY }), v3.GetY());
        dp->setValue(multipoleName({ ZZ }), v3.GetZ());
        mpt->LastExperiment()->addProperty(mpo, dp);
    }
    
    // Quadrupole
    auto my_quadrupole = mol->GetData("Traceless Quadrupole Moment");
    if (nullptr != my_quadrupole)
    {
        auto quadrupole = (OpenBabel::OBMatrixData *) my_quadrupole;
        OpenBabel::matrix3x3            m3 = quadrupole->GetData();
        double                          mm[9];
        m3.GetArray(mm);
        auto mpo = MolPropObservable::QUADRUPOLE;
        auto mq  = new alexandria::MolecularMultipole(qm_type, "B", 0.0, mpo);
        mq->setValue(multipoleName({ XX, XX }), mm[0]);
        mq->setValue(multipoleName({ XX, YY }), mm[1]);
        mq->setValue(multipoleName({ XX, ZZ }), mm[2]);
        mq->setValue(multipoleName({ YY, YY }), mm[4]);
        mq->setValue(multipoleName({ YY, ZZ }), mm[5]);
        mq->setValue(multipoleName({ ZZ, ZZ }), mm[8]);
        mpt->LastExperiment()->addProperty(mpo, mq);
    }

    // Polarizability
    auto my_pol_tensor = mol->GetData("Exact polarizability");
    if (nullptr != my_pol_tensor)
    {
        auto pol_tensor = (OpenBabel::OBMatrixData *) my_pol_tensor;
        OpenBabel::matrix3x3 m3 = pol_tensor->GetData();
        double               mm[9], alpha;
        m3.GetArray(mm);
        alpha = (mm[0]+mm[4]+mm[8])/3.0;

        auto mdp = new alexandria::MolecularPolarizability(qm_type, "Bohr3",
                                                           0.0, mm[0], mm[4], mm[8],
                                                           mm[1], mm[2], mm[5], alpha, 0);
        mpt->LastExperiment()->addProperty(MolPropObservable::POLARIZABILITY, mdp);
    }

    // Electrostatic potential
    auto esp = mol->GetData("Electrostatic Potential");
    if (nullptr != esp && maxPotential > 0)
    {
        auto                                            espptr = (OpenBabel::OBFreeGrid *) esp;
        OpenBabel::OBFreeGridPoint                     *fgp;
        OpenBabel::OBFreeGridPointIterator              fgpi;
        int                                             espid  = 0;
        auto espv = new ElectrostaticPotential("Angstrom", "Hartree/e");
        
        fgpi = espptr->BeginPoints();
        for (fgp = espptr->BeginPoint(fgpi); (nullptr != fgp); fgp = espptr->NextPoint(fgpi))
        {
            espv->addPoint(++espid, fgp->GetX(), fgp->GetY(), fgp->GetZ(), fgp->GetV());
        }
        mpt->LastExperiment()->addProperty(MolPropObservable::POTENTIAL, espv);
    }
    return true;
}

static bool readBabel(const std::string               &g09,
                      std::vector<OpenBabel::OBMol *> *mols,
                      einformat                       *inputformat)
{
    std::string baseFileName(g09);
    std::string fileName(g09);
    bool        isGzip = false;
    if (!gmx_fexist(fileName.c_str()))
    {
        fileName += ".gz";
        printf("fileName2: %s\n", fileName.c_str());
        if (!gmx_fexist(fileName.c_str()))
        {
            fprintf(stderr, "Can neither find file %s or %s\n", g09.c_str(), fileName.c_str());
            return false;
        }
        isGzip = true;
    }
    else
    {
        isGzip = isGzipFile(g09, &baseFileName);
    }
    
    std::ifstream g09f;
    
    g09f.open(fileName, std::ios::in);
    if (!g09f.is_open())
    {
        fprintf(stderr, "Cannot open file %s for reading", g09.c_str());
        return false;
    }
    // Read from g09f
    auto        conv        = new OpenBabel::OBConversion(&g09f, &std::cout);
    BabelFiles  babelfiles;
    std::string informat;
    if (isGzip)
    {
        informat = babelfiles.findBabelFile(baseFileName)->informat();
    }
    else
    {
        informat = babelfiles.findBabelFile(fileName)->informat();
    }

    if (informat == "g03" || informat == "g09" || informat == "g16")
    {
        *inputformat = einfGaussian;
    }

    bool read_ok = false;
    if (conv->SetInFormat(informat.c_str(), isGzip))
    {
        try
        {
            do
            {
                read_ok = false;
                if (mols->empty())
                {
                    auto mol = new OpenBabel::OBMol;
                    if (conv->ReadFile(mol, fileName) && !mol->Empty())
                    {
                        mols->push_back(mol);
                        read_ok = true;
                    }
                }
                else
                {
                    auto mol = new OpenBabel::OBMol;
                    if (conv->Read(mol, nullptr) && !mol->Empty())
                    {
                        mols->push_back(mol);
                        read_ok = true;
                    }
                }
                // TODO: check when we need to read just one molecule
                if (*inputformat == einfGaussian)
                {
                    read_ok = false;
                }
            }
            while (read_ok);

            // If we read anything at all we are happy.
            read_ok = mols->size() > 0;
        }
        catch (const std::exception &ex)
        {
            gmx::printFatalErrorMessage(stderr, ex);
        }
    }
    else
    {
        fprintf(stderr, "Input file %s has incomprehensible format.\n", g09.c_str());
    }
    if (read_ok)
    {
        g09f.close();
    }
    return read_ok;
}

static void OBUnitCell2Box(OpenBabel::OBUnitCell *ob, matrix box)
{
    double fa = 0.1*ob->GetA();
    double fb = 0.1*ob->GetB();
    double fc = 0.1*ob->GetC();
    double alpha = ob->GetAlpha();
    double beta  = ob->GetBeta();
    double gamma = ob->GetGamma();
    clear_mat(box);
    box[XX][XX] = fa;
    if ((alpha != 90.0) || (beta != 90.0) || (gamma != 90.0))
    {
        double cosa = 0, cosb = 0, cosg = 0, sing = 0;
        if (alpha != 90.0)
        {
            cosa = std::cos(alpha*DEG2RAD);
        }
        else
        {
            cosa = 0;
        }
        if (beta != 90.0)
        {
            cosb = std::cos(beta*DEG2RAD);
        }
        else
        {
            cosb = 0;
        }
        if (gamma != 90.0)
        {
            cosg = std::cos(gamma*DEG2RAD);
            sing = std::sin(gamma*DEG2RAD);
        }
        else
        {
            cosg = 0;
            sing = 1;
        }
        box[YY][XX] = fb*cosg;
        box[YY][YY] = fb*sing;
        box[ZZ][XX] = fc*cosb;
        box[ZZ][YY] = fc*(cosa - cosb*cosg)/sing;
        box[ZZ][ZZ] = std::sqrt(fc*fc
                                - box[ZZ][XX]*box[ZZ][XX] - box[ZZ][YY]*box[ZZ][YY]);
    }
    else
    {
        box[YY][YY] = fb;
        box[ZZ][ZZ] = fc;
    }
}

bool readBabel(const ForceField                 *pd,
               const char                       *g09,
               std::vector<alexandria::MolProp> *mpt,
               const char                       *molnm,
               const char                       *iupac,
               const char                       *conformation,
               std::string                      *method,
               std::string                      *basisset,
               int                               maxPotential,
               int                               nsymm,
               const char                      *jobType,
               double                          *qtot,
               bool                             addHydrogen,
               matrix                           box,
               bool                             renameAtoms)
{
    std::vector<OpenBabel::OBMol *> mols;
    einformat                       inputformat = einfNotGaussian;
    bool                            read_ok     = readBabel(g09, &mols, &inputformat);
    if (!read_ok)
    {
        fprintf(stderr, "Failed reading %s\n", g09);
        return false;
    }
    std::map<std::string, std::string> g2a;
    if (renameAtoms)
    {
        gaffToAlexandria("", &g2a);
    }
    AlexandriaMols amols;
    
    for(auto &mol : mols)
    {
        alexandria::MolProp mp;
        if (babel2ACT(pd, g2a, amols, mol, &mp, molnm, iupac, conformation, method, basisset, 
                      maxPotential, nsymm, jobType, qtot, addHydrogen, g09,
                      inputformat))
        {
            mpt->push_back(mp);
        }
        //delete mol;
    }
    if (mols[0]->HasData(OpenBabel::OBGenericDataType::UnitCell))
    {
        auto OBUC = (OpenBabel::OBUnitCell *)mols[0]->GetData(OpenBabel::OBGenericDataType::UnitCell);
        OBUnitCell2Box(OBUC, box);
    }
    return true;
}

bool SetMolpropAtomTypesAndBonds(alexandria::MolProp *mmm)
{
    OpenBabel::OBMol mol;
    mol.BeginModify();
    auto frags = mmm->fragments();
    if (!frags.empty())
    {
        mol.SetTotalCharge(frags[0].charge());
        mol.SetTotalSpinMultiplicity(frags[0].multiplicity());
    }
    auto ei  = mmm->experiment()->begin();
    mol.ReserveAtoms(ei->NAtom());
    int  idx = 0;
    for (auto &ca : ei->calcAtomConst())
    {
        int  atomicNum = OpenBabel::OBElements::GetAtomicNum(ca.getName().c_str());
        auto atom      = mol.NewAtom(idx++);
        atom->SetVector(ca.getX(), ca.getY(), ca.getZ());
        atom->SetAtomicNum(atomicNum);
    }
    mol.ConnectTheDots();
    mol.PerceiveBondOrders();
    mol.EndModify();
    if (!frags.empty())
    {
        frags[0].setMass(mol.GetMolWt());
        frags[0].setFormula(mol.GetFormula());
    }
    const char *forcefield = "alexandria";
    auto        pFF        = OpenBabel::OBForceField::FindForceField(forcefield);
    if (!pFF)
    {
        fprintf(stderr, "could not find forcefield %s\n", forcefield);
        return false;
    }
    if (debug)
    {
        pFF->SetLogFile(&std::cout);
        pFF->SetLogLevel(OBFF_LOGLVL_HIGH);
    }
    if (!pFF->Setup(mol))
    {
        fprintf(stderr, "could not setup force field %s.\n", forcefield);
        return false;
    }
    pFF->GetAtomTypes(mol);
    auto ca = ei->calcAtom()->begin();
    for (OpenBabel::OBMolAtomIter a(mol); a; ++a, ++ca)
    {
        auto type = static_cast<OpenBabel::OBPairData*>(a->GetData("FFAtomType"));
        ca->setObtype(type->GetValue());
    }
    // Bonds
    getBondsFromOpenBabel(&mol, mmm, "psi4 input file", true);

    return true;
}

#endif
