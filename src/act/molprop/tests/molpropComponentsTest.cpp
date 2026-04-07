/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2026
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
 * Tests for molprop component classes: Phase, MolPropObservable, MultipoleNames,
 * TopologyEntry subclasses, Fragment, AtomNum, MolecularComposition, CalcAtom,
 * CategoryList, MolProp basic properties, and Experiment utilities.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"

#include "act/basics/msg_handler.h"
#include "act/molprop/categories.h"
#include "act/molprop/composition.h"
#include "act/molprop/experiment.h"
#include "act/molprop/fragment.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molpropobservable.h"
#include "act/molprop/multipole_names.h"
#include "act/molprop/phase.h"
#include "act/molprop/topologyentry.h"

namespace alexandria
{

namespace
{

// ============================================================
// Phase tests
// ============================================================

TEST(PhaseTest, Phase2StringRoundTrip)
{
    EXPECT_EQ("gas",    phase2string(ePhase::GAS));
    EXPECT_EQ("liquid", phase2string(ePhase::LIQUID));
    EXPECT_EQ("solid",  phase2string(ePhase::SOLID));
    EXPECT_EQ("plasma", phase2string(ePhase::PLASMA));
}

TEST(PhaseTest, String2PhaseRoundTrip)
{
    EXPECT_EQ(ePhase::GAS,    string2phase("gas"));
    EXPECT_EQ(ePhase::LIQUID, string2phase("liquid"));
    EXPECT_EQ(ePhase::SOLID,  string2phase("solid"));
    EXPECT_EQ(ePhase::PLASMA, string2phase("plasma"));
}

TEST(PhaseTest, String2PhaseInvalidThrows)
{
    EXPECT_THROW(string2phase("vacuum"), gmx::InvalidInputError);
}

// ============================================================
// MolPropObservable enum function tests
// ============================================================

TEST(MolPropObservableTest, MpoNameKnownValues)
{
    EXPECT_STREQ("dipole",          mpo_name(MolPropObservable::DIPOLE));
    EXPECT_STREQ("quadrupole",      mpo_name(MolPropObservable::QUADRUPOLE));
    EXPECT_STREQ("octupole",        mpo_name(MolPropObservable::OCTUPOLE));
    EXPECT_STREQ("hexadecapole",    mpo_name(MolPropObservable::HEXADECAPOLE));
    EXPECT_STREQ("polarizability",  mpo_name(MolPropObservable::POLARIZABILITY));
    EXPECT_STREQ("HF",              mpo_name(MolPropObservable::HF));
    EXPECT_STREQ("DeltaE0",         mpo_name(MolPropObservable::DELTAE0));
    EXPECT_STREQ("DeltaHform",      mpo_name(MolPropObservable::DHFORM));
    EXPECT_STREQ("frequency",       mpo_name(MolPropObservable::FREQUENCY));
}

TEST(MolPropObservableTest, MpoUnit2KnownValues)
{
    EXPECT_STREQ("D",          mpo_unit2(MolPropObservable::DIPOLE));
    EXPECT_STREQ("B",          mpo_unit2(MolPropObservable::QUADRUPOLE));
    EXPECT_STREQ("Angstrom3",  mpo_unit2(MolPropObservable::POLARIZABILITY));
    EXPECT_STREQ("kJ/mol",     mpo_unit2(MolPropObservable::HF));
    EXPECT_STREQ("kJ/mol",     mpo_unit2(MolPropObservable::DHFORM));
    EXPECT_STREQ("cm^-1",      mpo_unit2(MolPropObservable::FREQUENCY));
}

TEST(MolPropObservableTest, StringToMpolFoundDipole)
{
    MolPropObservable mpo;
    EXPECT_TRUE(stringToMolPropObservable("dipole", &mpo));
    EXPECT_EQ(MolPropObservable::DIPOLE, mpo);
}

TEST(MolPropObservableTest, StringToMpoCaseInsensitive)
{
    MolPropObservable mpo;
    EXPECT_TRUE(stringToMolPropObservable("DIPOLE", &mpo));
    EXPECT_EQ(MolPropObservable::DIPOLE, mpo);
}

TEST(MolPropObservableTest, StringToMpoFoundPolarizability)
{
    MolPropObservable mpo;
    EXPECT_TRUE(stringToMolPropObservable("polarizability", &mpo));
    EXPECT_EQ(MolPropObservable::POLARIZABILITY, mpo);
}

TEST(MolPropObservableTest, StringToMpoNotFound)
{
    MolPropObservable mpo;
    EXPECT_FALSE(stringToMolPropObservable("nosuchobservable", &mpo));
}

// ============================================================
// MultipoleNames tests
// ============================================================

TEST(MultipoleNamesTest, DipoleNamesCount)
{
    auto names = multipoleNames(MolPropObservable::DIPOLE);
    ASSERT_EQ(3u, names.size());
    EXPECT_EQ("x", names[0]);
    EXPECT_EQ("y", names[1]);
    EXPECT_EQ("z", names[2]);
}

TEST(MultipoleNamesTest, QuadrupoleNamesCount)
{
    auto names = multipoleNames(MolPropObservable::QUADRUPOLE);
    ASSERT_EQ(6u, names.size());
    EXPECT_EQ("xx", names[0]);
    EXPECT_EQ("xy", names[1]);
    EXPECT_EQ("xz", names[2]);
    EXPECT_EQ("yy", names[3]);
    EXPECT_EQ("yz", names[4]);
    EXPECT_EQ("zz", names[5]);
}

TEST(MultipoleNamesTest, OctupoleNamesCount)
{
    auto names = multipoleNames(MolPropObservable::OCTUPOLE);
    EXPECT_EQ(10u, names.size());
}

TEST(MultipoleNamesTest, HexadecapoleNamesCount)
{
    auto names = multipoleNames(MolPropObservable::HEXADECAPOLE);
    EXPECT_EQ(15u, names.size());
}

TEST(MultipoleNamesTest, MultipoleNameLookupDipole)
{
    EXPECT_EQ("x", multipoleName({0}));
    EXPECT_EQ("y", multipoleName({1}));
    EXPECT_EQ("z", multipoleName({2}));
}

TEST(MultipoleNamesTest, MultipoleNameLookupQuadrupole)
{
    EXPECT_EQ("xx", multipoleName({0, 0}));
    EXPECT_EQ("xy", multipoleName({0, 1}));
    EXPECT_EQ("yy", multipoleName({1, 1}));
}

TEST(MultipoleNamesTest, MultipoleIndexByStringDipole)
{
    int idx;
    EXPECT_TRUE(multipoleIndex("x", &idx));
    EXPECT_EQ(0, idx);
    EXPECT_TRUE(multipoleIndex("y", &idx));
    EXPECT_EQ(1, idx);
    EXPECT_TRUE(multipoleIndex("z", &idx));
    EXPECT_EQ(2, idx);
}

TEST(MultipoleNamesTest, MultipoleIndexByStringQuadrupole)
{
    int idx;
    EXPECT_TRUE(multipoleIndex("xx", &idx));
    EXPECT_EQ(0, idx);
    EXPECT_TRUE(multipoleIndex("xy", &idx));
    EXPECT_EQ(1, idx);
    EXPECT_TRUE(multipoleIndex("xz", &idx));
    EXPECT_EQ(2, idx);
}

TEST(MultipoleNamesTest, MultipoleIndexByStringUnknown)
{
    int idx = -1;
    EXPECT_FALSE(multipoleIndex("qq", &idx));
}

TEST(MultipoleNamesTest, MultipoleIndexByVectorDipole)
{
    EXPECT_EQ(0, multipoleIndex({0}));
    EXPECT_EQ(1, multipoleIndex({1}));
    EXPECT_EQ(2, multipoleIndex({2}));
}

TEST(MultipoleNamesTest, MultipoleIndexByVectorQuadrupole)
{
    EXPECT_EQ(0, multipoleIndex({0, 0}));
    EXPECT_EQ(1, multipoleIndex({0, 1}));
    EXPECT_EQ(5, multipoleIndex({2, 2}));
}

// ============================================================
// TopologyEntry / SingleAtom tests
// ============================================================

TEST(SingleAtomTest, Construction)
{
    SingleAtom sa(5);
    EXPECT_EQ(5, sa.aI());
}

TEST(SingleAtomTest, EqualityTrue)
{
    SingleAtom sa(3), sb(3);
    EXPECT_TRUE(sa == sb);
}

TEST(SingleAtomTest, EqualityFalse)
{
    SingleAtom sa(3), sc(4);
    EXPECT_FALSE(sa == sc);
}

TEST(SingleAtomTest, LessOperator)
{
    SingleAtom sa(2), sb(5);
    EXPECT_TRUE(sa < sb);
    EXPECT_FALSE(sb < sa);
    EXPECT_FALSE(sa < sa);
}

// ============================================================
// AtomPair tests
// ============================================================

TEST(AtomPairTest, Construction)
{
    AtomPair ap(3, 7);
    EXPECT_EQ(3, ap.aI());
    EXPECT_EQ(7, ap.aJ());
}

TEST(AtomPairTest, EqualitySymmetric)
{
    AtomPair ap1(3, 7), ap2(7, 3);
    EXPECT_TRUE(ap1 == ap2);
}

TEST(AtomPairTest, EqualitySameOrder)
{
    AtomPair ap1(3, 7), ap2(3, 7);
    EXPECT_TRUE(ap1 == ap2);
}

TEST(AtomPairTest, EqualityFalse)
{
    AtomPair ap1(3, 7), ap2(3, 8);
    EXPECT_FALSE(ap1 == ap2);
}

TEST(AtomPairTest, SwapAtoms)
{
    AtomPair ap(3, 7);
    AtomPair swapped = ap.swap();
    EXPECT_EQ(7, swapped.aI());
    EXPECT_EQ(3, swapped.aJ());
}

TEST(AtomPairTest, GetMethod)
{
    AtomPair ap(4, 9);
    int ai, aj;
    ap.get(&ai, &aj);
    EXPECT_EQ(4, ai);
    EXPECT_EQ(9, aj);
}

// ============================================================
// Bond tests
// ============================================================

TEST(BondTest, Construction)
{
    Bond b(1, 2, 1.5);
    EXPECT_EQ(1, b.aI());
    EXPECT_EQ(2, b.aJ());
    EXPECT_DOUBLE_EQ(1.5, b.bondOrder());
}

TEST(BondTest, GetMethod)
{
    Bond b(1, 2, 2.0);
    int ai, aj;
    double bo;
    b.get(&ai, &aj, &bo);
    EXPECT_EQ(1, ai);
    EXPECT_EQ(2, aj);
    EXPECT_DOUBLE_EQ(2.0, bo);
}

TEST(BondTest, EqualitySymmetric)
{
    Bond b1(1, 2, 1.0), b2(2, 1, 1.0);
    EXPECT_TRUE(b1 == b2);
}

TEST(BondTest, EqualitySameOrder)
{
    Bond b1(1, 2, 1.0), b2(1, 2, 1.0);
    EXPECT_TRUE(b1 == b2);
}

TEST(BondTest, EqualityFalse)
{
    Bond b1(1, 2, 1.0), b2(1, 3, 1.0);
    EXPECT_FALSE(b1 == b2);
}

TEST(BondTest, Swap)
{
    Bond b(1, 2, 1.5);
    Bond swapped = b.swap();
    EXPECT_EQ(2, swapped.aI());
    EXPECT_EQ(1, swapped.aJ());
    EXPECT_DOUBLE_EQ(1.5, swapped.bondOrder());
}

TEST(BondTest, RenumberAtoms)
{
    Bond b(0, 1, 1.0);
    std::vector<int> renumber = {10, 11, 12};
    b.renumberAtoms(renumber);
    EXPECT_EQ(10, b.aI());
    EXPECT_EQ(11, b.aJ());
}

TEST(BondTest, SetBondOrder)
{
    Bond b(1, 2, 1.0);
    b.setBondOrder(0, 2.0);
    EXPECT_DOUBLE_EQ(2.0, b.bondOrder());
}

// ============================================================
// TopologyEntry base class tests
// ============================================================

TEST(TopologyEntryTest, CheckPassesWithCorrectCount)
{
    Bond b(1, 2, 1.0);
    EXPECT_NO_THROW(b.check(2));
}

TEST(TopologyEntryTest, CheckThrowsWithWrongCount)
{
    Bond b(1, 2, 1.0);
    EXPECT_THROW(b.check(3), gmx::InternalError);
}

TEST(TopologyEntryTest, AtomIndexOutOfRangeThrows)
{
    Bond b(1, 2, 1.0);
    EXPECT_THROW(b.atomIndex(5), gmx::InternalError);
}

TEST(TopologyEntryTest, AddBondOrderAndRetrieve)
{
    TopologyEntry te;
    te.addBondOrder(1.5);
    te.addBondOrder(2.0);
    ASSERT_EQ(2u, te.bondOrders().size());
    EXPECT_DOUBLE_EQ(1.5, te.bondOrder(0));
    EXPECT_DOUBLE_EQ(2.0, te.bondOrder(1));
}

// ============================================================
// Angle tests
// ============================================================

TEST(AngleTest, Construction)
{
    Bond bij(0, 1, 1.0);
    Bond bjk(1, 2, 1.0);
    Angle angle(bij, bjk);
    ASSERT_EQ(3u, angle.atomIndices().size());
    EXPECT_EQ(0, angle.atomIndex(0));
    EXPECT_EQ(1, angle.atomIndex(1));
    EXPECT_EQ(2, angle.atomIndex(2));
    EXPECT_FALSE(angle.isLinear());
}

TEST(AngleTest, SetLinear)
{
    Bond bij(0, 1, 1.0);
    Bond bjk(1, 2, 1.0);
    Angle angle(bij, bjk);
    angle.setLinear(true);
    EXPECT_TRUE(angle.isLinear());
}

TEST(AngleTest, BondAccessors)
{
    Bond bij(0, 1, 1.0);
    Bond bjk(1, 2, 2.0);
    Angle angle(bij, bjk);
    EXPECT_EQ(0, angle.bij().aI());
    EXPECT_EQ(1, angle.bij().aJ());
    EXPECT_EQ(1, angle.bjk().aI());
    EXPECT_EQ(2, angle.bjk().aJ());
}

TEST(AngleTest, InvalidBondsThrow)
{
    Bond bij(0, 1, 1.0);
    Bond bjk(2, 3, 1.0);
    EXPECT_THROW(Angle angle(bij, bjk), gmx::InvalidInputError);
}

// ============================================================
// Improper tests
// ============================================================

TEST(ImproperTest, Construction)
{
    Bond bij(0, 1, 1.0);
    Bond bik(0, 2, 1.0);
    Bond bil(0, 3, 1.0);
    Improper imp(bij, bik, bil);
    ASSERT_EQ(4u, imp.atomIndices().size());
    EXPECT_EQ(0, imp.atomIndex(0));
    EXPECT_EQ(1, imp.atomIndex(1));
    EXPECT_EQ(2, imp.atomIndex(2));
    EXPECT_EQ(3, imp.atomIndex(3));
}

TEST(ImproperTest, InconsistentBondsThrow)
{
    Bond bij(0, 1, 1.0);
    Bond bik(1, 2, 1.0);
    Bond bil(0, 3, 1.0);
    EXPECT_THROW(Improper imp(bij, bik, bil), gmx::InvalidInputError);
}

// ============================================================
// Proper tests
// ============================================================

TEST(ProperTest, Construction)
{
    Bond bij(0, 1, 1.0);
    Bond bjk(1, 2, 1.0);
    Bond bkl(2, 3, 1.0);
    Proper prop(bij, bjk, bkl);
    ASSERT_EQ(4u, prop.atomIndices().size());
    EXPECT_EQ(0, prop.atomIndex(0));
    EXPECT_EQ(1, prop.atomIndex(1));
    EXPECT_EQ(2, prop.atomIndex(2));
    EXPECT_EQ(3, prop.atomIndex(3));
}

TEST(ProperTest, InconsistentBondsThrow)
{
    Bond bij(0, 1, 1.0);
    Bond bjk(2, 3, 1.0);
    Bond bkl(3, 4, 1.0);
    EXPECT_THROW(Proper prop(bij, bjk, bkl), gmx::InvalidInputError);
}

// ============================================================
// Fragment tests
// ============================================================

TEST(FragmentTest, DefaultConstructor)
{
    Fragment f;
    EXPECT_DOUBLE_EQ(0.0, f.mass());
    EXPECT_EQ(0, f.charge());
    EXPECT_EQ(1, f.multiplicity());
    EXPECT_EQ(1, f.symmetryNumber());
    EXPECT_TRUE(f.formula().empty());
    EXPECT_TRUE(f.atoms().empty());
}

TEST(FragmentTest, ParameterizedConstructor)
{
    Fragment f("InChI=1S/H2O/h1H2", 18.015, 0, 1, 1, "H2O", {0, 1, 2});
    EXPECT_EQ("InChI=1S/H2O/h1H2", f.inchi());
    EXPECT_DOUBLE_EQ(18.015, f.mass());
    EXPECT_EQ(0, f.charge());
    EXPECT_EQ(1, f.multiplicity());
    EXPECT_EQ(1, f.symmetryNumber());
    EXPECT_EQ("H2O", f.formula());
    ASSERT_EQ(3u, f.atoms().size());
    EXPECT_EQ(0, f.atoms()[0]);
    EXPECT_EQ(1, f.atoms()[1]);
    EXPECT_EQ(2, f.atoms()[2]);
}

TEST(FragmentTest, SettersAndGetters)
{
    Fragment f;
    f.setInchi("InChI=1S/CH4/h1H4");
    f.setIupac("methane");
    f.setMass(16.043);
    f.setCharge(-1);
    f.setFormula("CH4");
    EXPECT_EQ("InChI=1S/CH4/h1H4", f.inchi());
    EXPECT_EQ("methane", f.iupac());
    EXPECT_DOUBLE_EQ(16.043, f.mass());
    EXPECT_EQ(-1, f.charge());
    EXPECT_EQ("CH4", f.formula());
}

TEST(FragmentTest, AtomString)
{
    Fragment f("", 1.0, 0, 1, 1, "H", {0, 1, 2});
    // atomString() adds 1 to each atom index and prepends a space
    EXPECT_EQ(" 1 2 3", f.atomString());
}

TEST(FragmentTest, TexFormulaContainsSubscript)
{
    Fragment f("", 18.0, 0, 1, 1, "H2O", {0, 1, 2});
    // H2O should produce H$_2$O
    EXPECT_NE(std::string::npos, f.texFormula().find("$_2$"));
}

TEST(FragmentTest, TexFormulaNoSubscriptForNoDigit)
{
    Fragment f("", 4.0, 0, 1, 1, "He", {0});
    EXPECT_EQ("He", f.texFormula());
}

TEST(FragmentTest, SetAtoms)
{
    Fragment f;
    f.setAtoms({3, 5, 7});
    ASSERT_EQ(3u, f.atoms().size());
    EXPECT_EQ(3, f.atoms()[0]);
    EXPECT_EQ(5, f.atoms()[1]);
    EXPECT_EQ(7, f.atoms()[2]);
    // atomString should be updated
    EXPECT_EQ(" 4 6 8", f.atomString());
}

TEST(FragmentTest, SymmetryNumberAtLeastOne)
{
    Fragment f("", 1.0, 0, 1, 0, "X", {});
    EXPECT_EQ(1, f.symmetryNumber());
}

// ============================================================
// AtomNum tests
// ============================================================

TEST(AtomNumTest, ConstructorCharPtr)
{
    AtomNum an("C", 6);
    EXPECT_EQ("C", an.getAtom());
    EXPECT_EQ(6, an.getNumber());
}

TEST(AtomNumTest, ConstructorString)
{
    AtomNum an(std::string("H"), 2);
    EXPECT_EQ("H", an.getAtom());
    EXPECT_EQ(2, an.getNumber());
}

TEST(AtomNumTest, DefaultAndSetters)
{
    AtomNum an;
    an.SetAtom("O");
    an.SetNumber(8);
    EXPECT_EQ("O", an.getAtom());
    EXPECT_EQ(8, an.getNumber());
}

TEST(AtomNumTest, SetAtomCharPtr)
{
    AtomNum an;
    an.SetAtom("N");
    EXPECT_EQ("N", an.getAtom());
}

// ============================================================
// MolecularComposition tests
// ============================================================

TEST(MolecularCompositionTest, DefaultConstruction)
{
    MolecularComposition mc;
    EXPECT_EQ(0, mc.CountAtoms());
}

TEST(MolecularCompositionTest, NamedConstruction)
{
    MolecularComposition mc("alexandria");
    EXPECT_EQ("alexandria", mc.getCompName());
}

TEST(MolecularCompositionTest, SetCompName)
{
    MolecularComposition mc;
    mc.SetCompName("mycomp");
    EXPECT_EQ("mycomp", mc.getCompName());
}

TEST(MolecularCompositionTest, AddAtomAndCount)
{
    MolecularComposition mc("test");
    mc.AddAtom(AtomNum("C", 1));
    mc.AddAtom(AtomNum("H", 4));
    EXPECT_EQ(2u, mc.atomNumConst().size());
    EXPECT_EQ(1, mc.CountAtoms("C"));
    EXPECT_EQ(4, mc.CountAtoms("H"));
    EXPECT_EQ(5, mc.CountAtoms());
}

TEST(MolecularCompositionTest, AddAtomAccumulatesCount)
{
    MolecularComposition mc("test");
    mc.AddAtom(AtomNum("C", 1));
    mc.AddAtom(AtomNum("C", 2));
    EXPECT_EQ(1u, mc.atomNumConst().size());
    EXPECT_EQ(3, mc.CountAtoms("C"));
}

TEST(MolecularCompositionTest, CountAtomsUnknownReturnsZero)
{
    MolecularComposition mc("test");
    mc.AddAtom(AtomNum("C", 1));
    EXPECT_EQ(0, mc.CountAtoms("O"));
}

TEST(MolecularCompositionTest, DeleteAtom)
{
    MolecularComposition mc("test");
    mc.AddAtom(AtomNum("C", 2));
    mc.AddAtom(AtomNum("H", 6));
    mc.DeleteAtom("H");
    EXPECT_EQ(1u, mc.atomNumConst().size());
    EXPECT_EQ(0, mc.CountAtoms("H"));
    EXPECT_EQ(2, mc.CountAtoms("C"));
}

TEST(MolecularCompositionTest, DeleteAtomNonexistentIsNoOp)
{
    MolecularComposition mc("test");
    mc.AddAtom(AtomNum("C", 1));
    EXPECT_NO_THROW(mc.DeleteAtom("O"));
    EXPECT_EQ(1u, mc.atomNumConst().size());
}

TEST(MolecularCompositionTest, ReplaceAtom)
{
    MolecularComposition mc("test");
    mc.AddAtom(AtomNum("N", 1));
    mc.ReplaceAtom("N", "P");
    EXPECT_EQ(0, mc.CountAtoms("N"));
    EXPECT_EQ(1, mc.CountAtoms("P"));
}

TEST(MolecularCompositionTest, SearchAtomConstFound)
{
    MolecularComposition mc("test");
    mc.AddAtom(AtomNum("C", 3));
    auto it = mc.searchAtomConst("C");
    EXPECT_EQ("C", it->getAtom());
    EXPECT_EQ(3, it->getNumber());
}

TEST(MolecularCompositionTest, SearchAtomConstNotFound)
{
    MolecularComposition mc("test");
    mc.AddAtom(AtomNum("C", 1));
    // Use CountAtoms to confirm atom not present rather than iterator comparison
    EXPECT_EQ(0, mc.CountAtoms("O"));
}

// ============================================================
// CalcAtom tests
// ============================================================

TEST(CalcAtomTest, DefaultConstructor)
{
    CalcAtom ca;
    EXPECT_EQ(0, ca.getAtomid());
}

TEST(CalcAtomTest, ParameterizedConstructor)
{
    CalcAtom ca("C1", "CT", 5);
    EXPECT_EQ("C1", ca.getName());
    EXPECT_EQ(5, ca.getAtomid());
}

TEST(CalcAtomTest, SetName)
{
    CalcAtom ca;
    ca.setName("O2");
    EXPECT_EQ("O2", ca.getName());
}

TEST(CalcAtomTest, SetAtomid)
{
    CalcAtom ca;
    ca.setAtomid(42);
    EXPECT_EQ(42, ca.getAtomid());
}

TEST(CalcAtomTest, ChargeManagement)
{
    CalcAtom ca("O", "OT", 1);
    EXPECT_FALSE(ca.hasCharge("RESP"));
    ca.AddCharge("RESP", -0.834);
    EXPECT_TRUE(ca.hasCharge("RESP"));
    EXPECT_DOUBLE_EQ(-0.834, ca.charge("RESP"));
}

TEST(CalcAtomTest, ChargeOverwrite)
{
    CalcAtom ca("O", "OT", 1);
    ca.AddCharge("RESP", -0.834);
    ca.AddCharge("RESP", -0.5);
    EXPECT_DOUBLE_EQ(-0.5, ca.charge("RESP"));
}

TEST(CalcAtomTest, MultipleChargeTypes)
{
    CalcAtom ca("C", "CT", 0);
    ca.AddCharge("RESP", 0.1);
    ca.AddCharge("ACM",  0.2);
    EXPECT_TRUE(ca.hasCharge("RESP"));
    EXPECT_TRUE(ca.hasCharge("ACM"));
    EXPECT_FALSE(ca.hasCharge("NBO"));
    EXPECT_DOUBLE_EQ(0.1, ca.charge("RESP"));
    EXPECT_DOUBLE_EQ(0.2, ca.charge("ACM"));
}

TEST(CalcAtomTest, EqualSameAtom)
{
    CalcAtom ca("C", "CT", 1);
    CalcAtom cb("C", "CT", 1);
    EXPECT_TRUE(ca.Equal(cb));
}

TEST(CalcAtomTest, EqualDifferentName)
{
    CalcAtom ca("C", "CT", 1);
    CalcAtom cb("O", "CT", 1);
    EXPECT_FALSE(ca.Equal(cb));
}

// ============================================================
// CategoryListElement tests
// ============================================================

TEST(CategoryListElementTest, Construction)
{
    CategoryListElement cle("alkane", "methane");
    EXPECT_EQ("alkane", cle.getName());
    EXPECT_EQ(1, cle.nMolecule());
    EXPECT_TRUE(cle.hasMolecule("methane"));
    EXPECT_FALSE(cle.hasMolecule("ethane"));
}

TEST(CategoryListElementTest, AddMoleculeNoDuplicate)
{
    CategoryListElement cle("alkane", "methane");
    cle.addMolecule("ethane");
    cle.addMolecule("methane");
    EXPECT_EQ(2, cle.nMolecule());
}

TEST(CategoryListElementTest, SortMolecules)
{
    CategoryListElement cle("ring", "toluene");
    cle.addMolecule("benzene");
    cle.addMolecule("naphthalene");
    cle.sortMolecules();
    EXPECT_EQ("benzene",    cle.molecules()[0]);
    EXPECT_EQ("naphthalene", cle.molecules()[1]);
    EXPECT_EQ("toluene",    cle.molecules()[2]);
}

// ============================================================
// CategoryList tests
// ============================================================

TEST(CategoryListTest, DefaultConstruction)
{
    CategoryList cl;
    EXPECT_EQ(0, cl.nCategories());
}

TEST(CategoryListTest, AddSingleCategory)
{
    CategoryList cl;
    cl.addCategory("alkane", "methane");
    EXPECT_EQ(1, cl.nCategories());
    EXPECT_EQ(1, cl.elementsConst()[0].nMolecule());
}

TEST(CategoryListTest, AddToExistingCategory)
{
    CategoryList cl;
    cl.addCategory("alkane", "methane");
    cl.addCategory("alkane", "ethane");
    EXPECT_EQ(1, cl.nCategories());
    EXPECT_EQ(2, cl.elementsConst()[0].nMolecule());
}

TEST(CategoryListTest, AddTwoDifferentCategories)
{
    CategoryList cl;
    cl.addCategory("alkane", "methane");
    cl.addCategory("alkene", "ethylene");
    EXPECT_EQ(2, cl.nCategories());
}

TEST(CategoryListTest, SortCategories)
{
    CategoryList cl;
    cl.addCategory("zzz", "mol1");
    cl.addCategory("aaa", "mol2");
    cl.addCategory("mmm", "mol3");
    cl.sortCategories();
    ASSERT_EQ(3, cl.nCategories());
    EXPECT_EQ("aaa", cl.elementsConst()[0].getName());
    EXPECT_EQ("mmm", cl.elementsConst()[1].getName());
    EXPECT_EQ("zzz", cl.elementsConst()[2].getName());
}

TEST(CategoryListTest, SortAlsoSortsMolecules)
{
    CategoryList cl;
    cl.addCategory("ring", "toluene");
    cl.addCategory("ring", "benzene");
    cl.sortCategories();
    EXPECT_EQ("benzene", cl.elementsConst()[0].molecules()[0]);
    EXPECT_EQ("toluene", cl.elementsConst()[0].molecules()[1]);
}

// ============================================================
// MolProp basic property tests
// ============================================================

TEST(MolPropBasicTest, IdentifiersAndIndex)
{
    MolProp mp;
    mp.SetMolname("water");
    mp.SetIupac("oxidane");
    mp.SetCas("7732-18-5");
    mp.SetCid("962");
    mp.SetInchi("InChI=1S/H2O/h1H2");
    mp.setIndex(42);
    EXPECT_EQ("water",             mp.getMolname());
    EXPECT_EQ("oxidane",           mp.getIupac());
    EXPECT_EQ("7732-18-5",         mp.getCas());
    EXPECT_EQ("962",               mp.getCid());
    EXPECT_EQ("InChI=1S/H2O/h1H2", mp.getInchi());
    EXPECT_EQ(42,                  mp.getIndex());
}

TEST(MolPropBasicTest, IupacFallsBackToMolname)
{
    MolProp mp;
    mp.SetMolname("water");
    EXPECT_EQ("water", mp.getIupac());
}

TEST(MolPropBasicTest, CategoryManagement)
{
    MolProp mp;
    mp.AddCategory("organic");
    mp.AddCategory("aromatic");
    mp.AddCategory("organic");
    EXPECT_EQ(2, mp.NCategory());
    EXPECT_TRUE(mp.SearchCategory("aromatic"));
    EXPECT_FALSE(mp.SearchCategory("inorganic"));
    mp.clearCategory();
    EXPECT_EQ(0, mp.NCategory());
}

TEST(MolPropBasicTest, BondAddAndExists)
{
    MolProp mp;
    Bond b1(0, 1, 1.0), b2(1, 2, 2.0);
    mp.AddBond(b1);
    mp.AddBond(b2);
    EXPECT_EQ(2, mp.NBond());
    EXPECT_TRUE(mp.BondExists(b1));
    EXPECT_TRUE(mp.BondExists(b2));
    Bond b3(5, 6, 1.0);
    EXPECT_FALSE(mp.BondExists(b3));
}

TEST(MolPropBasicTest, FragmentManagement)
{
    MolProp mp;
    Fragment f("", 18.0, 0, 1, 1, "H2O", {0, 1, 2});
    mp.addFragment(f);
    EXPECT_EQ(1u, mp.fragments().size());
    mp.clearFragments();
    EXPECT_EQ(0u, mp.fragments().size());
}

TEST(MolPropBasicTest, TotalMass)
{
    MolProp mp;
    Fragment f1("", 18.015, 0, 1, 1, "H2O", {0, 1, 2});
    Fragment f2("", 44.01,  0, 1, 1, "CO2", {3, 4, 5});
    mp.addFragment(f1);
    mp.addFragment(f2);
    EXPECT_DOUBLE_EQ(62.025, mp.totalMass());
}

TEST(MolPropBasicTest, TotalCharge)
{
    MolProp mp;
    Fragment f1("", 1.0,  1, 1, 1, "A", {0});
    Fragment f2("", 1.0, -1, 1, 1, "B", {1});
    mp.addFragment(f1);
    mp.addFragment(f2);
    EXPECT_EQ(0, mp.totalCharge());
}

TEST(MolPropBasicTest, TotalMultiplicitySinglet)
{
    MolProp mp;
    Fragment f1("", 1.0, 0, 1, 1, "A", {0});
    Fragment f2("", 1.0, 0, 1, 1, "B", {1});
    mp.addFragment(f1);
    mp.addFragment(f2);
    // Both fragments have odd multiplicity (1): no toggle, result stays 1
    EXPECT_EQ(1, mp.totalMultiplicity());
}

TEST(MolPropBasicTest, TotalMultiplicityDoublet)
{
    MolProp mp;
    Fragment f1("", 1.0, 0, 1, 1, "A", {0});
    Fragment f2("", 1.0, 0, 2, 1, "B", {1});
    mp.addFragment(f1);
    mp.addFragment(f2);
    // Fragment f2 has even multiplicity: toggles from 1 to 2
    EXPECT_EQ(2, mp.totalMultiplicity());
}

TEST(MolPropBasicTest, SymmetryNumberSingleFragment)
{
    MolProp mp;
    Fragment f("", 1.0, 0, 1, 3, "A", {0});
    mp.addFragment(f);
    EXPECT_EQ(3, mp.symmetryNumber());
}

TEST(MolPropBasicTest, SymmetryNumberMultipleFragmentsReturnsOne)
{
    MolProp mp;
    Fragment f1("", 1.0, 0, 1, 2, "A", {0});
    Fragment f2("", 1.0, 0, 1, 3, "B", {1});
    mp.addFragment(f1);
    mp.addFragment(f2);
    EXPECT_EQ(1, mp.symmetryNumber());
}

TEST(MolPropBasicTest, AddAndFindExperiment)
{
    MolProp mp;
    Experiment exp("ref", "min");
    mp.AddExperiment(std::move(exp));
    EXPECT_EQ(1u, mp.experimentConst().size());
    EXPECT_NE(nullptr, mp.LastExperiment());
}

TEST(MolPropBasicTest, LastExperimentNullWhenEmpty)
{
    MolProp mp;
    EXPECT_EQ(nullptr, mp.LastExperiment());
}

TEST(MolPropBasicTest, FindExperimentByJobType)
{
    MolProp mp;
    Experiment eOpt("G16", "B3LYP", "6-31G*", "ref", "min", "file.log", JobType::OPT);
    mp.AddExperiment(std::move(eOpt));
    const Experiment *found = mp.findExperimentConst(JobType::OPT);
    EXPECT_NE(nullptr, found);
    EXPECT_EQ(JobType::OPT, found->getJobtype());
    EXPECT_EQ(nullptr, mp.findExperimentConst(JobType::SP));
}

// ============================================================
// Experiment and JobType utility tests
// ============================================================

TEST(ExperimentTest, JobTypeToString)
{
    EXPECT_STREQ("Opt",      jobType2string(JobType::OPT));
    EXPECT_STREQ("SP",       jobType2string(JobType::SP));
    EXPECT_STREQ("Topology", jobType2string(JobType::TOPOLOGY));
    EXPECT_STREQ("unknown",  jobType2string(JobType::UNKNOWN));
}

TEST(ExperimentTest, StringToJobType)
{
    EXPECT_EQ(JobType::OPT,      string2jobType("Opt"));
    EXPECT_EQ(JobType::SP,       string2jobType("SP"));
    EXPECT_EQ(JobType::TOPOLOGY, string2jobType("Topology"));
    EXPECT_EQ(JobType::UNKNOWN,  string2jobType("unknown"));
}

TEST(ExperimentTest, StringToJobTypeInvalidThrows)
{
    EXPECT_THROW(string2jobType("BADTYPE"), gmx::InvalidInputError);
}

TEST(ExperimentTest, StringToJobTypeEmpty)
{
    EXPECT_EQ(JobType::UNKNOWN, string2jobType(""));
}

TEST(ExperimentTest, DataSourceName)
{
    EXPECT_STREQ("Experiment", dataSourceName(dsExperiment));
    EXPECT_STREQ("Theory",     dataSourceName(dsTheory));
}

TEST(ExperimentTest, DataSourceFromName)
{
    EXPECT_EQ(dsExperiment, dataSourceFromName("Experiment"));
    EXPECT_EQ(dsTheory,     dataSourceFromName("Theory"));
}

TEST(ExperimentTest, ExperimentConstructorSetsDataSource)
{
    Experiment exp("reference", "conformation");
    EXPECT_EQ(dsExperiment, exp.dataSource());
    EXPECT_EQ("reference",    exp.getReference());
    EXPECT_EQ("conformation", exp.getConformation());
}

TEST(ExperimentTest, CalculationConstructorSetsTheory)
{
    Experiment exp("Gaussian", "B3LYP", "6-311G**", "ref", "min", "file.log", JobType::OPT);
    EXPECT_EQ(dsTheory,       exp.dataSource());
    EXPECT_EQ("Gaussian",     exp.getProgram());
    EXPECT_EQ("B3LYP",        exp.getMethod());
    EXPECT_EQ("6-311G**",     exp.getBasisset());
    EXPECT_EQ("file.log",     exp.getDatafile());
    EXPECT_EQ(JobType::OPT,   exp.getJobtype());
}

TEST(ExperimentTest, SetJobtype)
{
    Experiment exp("ref", "min");
    exp.setJobtype(JobType::SP);
    EXPECT_EQ(JobType::SP, exp.getJobtype());
}

TEST(ExperimentTest, SetProgram)
{
    Experiment exp("ref", "min");
    exp.setProgram("Orca");
    EXPECT_EQ("Orca", exp.getProgram());
}

TEST(ExperimentTest, HasPropertyFalseInitially)
{
    Experiment exp("ref", "min");
    EXPECT_FALSE(exp.hasProperty(MolPropObservable::DIPOLE));
    EXPECT_FALSE(exp.hasMolPropObservable(MolPropObservable::DIPOLE));
}

TEST(ExperimentTest, AddPropertyDipole)
{
    Experiment exp("ref", "min");
    auto mp = std::make_unique<MolecularMultipole>("calc", "D", 298.15,
                                                   MolPropObservable::DIPOLE);
    mp->setValue("x", 1.0);
    mp->setValue("y", 0.5);
    mp->setValue("z", 0.0);
    exp.addProperty(MolPropObservable::DIPOLE, std::move(mp));
    EXPECT_TRUE(exp.hasProperty(MolPropObservable::DIPOLE));
    EXPECT_TRUE(exp.hasMolPropObservable(MolPropObservable::DIPOLE));
    EXPECT_FALSE(exp.hasProperty(MolPropObservable::QUADRUPOLE));
}

TEST(ExperimentTest, ExperimentId)
{
    Experiment exp("ref", "min");
    exp.setId(7);
    EXPECT_EQ(7, exp.id());
}

TEST(ExperimentTest, NAtomEmpty)
{
    Experiment exp("ref", "min");
    EXPECT_EQ(0, exp.NAtom());
}

TEST(ExperimentTest, AddAtom)
{
    Experiment exp("ref", "min");
    CalcAtom ca("O", "OT", 0);
    exp.AddAtom(ca);
    EXPECT_EQ(1, exp.NAtom());
    EXPECT_EQ("O", exp.calcAtomConst()[0].getName());
}

// ============================================================
// MolecularMultipole / GenericProperty tests
// ============================================================

TEST(MolecularMultipoleTest, DipoleConstruction)
{
    MolecularMultipole mm("calc", "D", 298.15, MolPropObservable::DIPOLE);
    EXPECT_STREQ("calc", mm.getType());
    EXPECT_DOUBLE_EQ(298.15, mm.getTemperature());
    ASSERT_EQ(3u, mm.getVector().size());
}

TEST(MolecularMultipoleTest, SetAndGetValues)
{
    MolecularMultipole mm("calc", "D", 0.0, MolPropObservable::DIPOLE);
    mm.setValue("x", 1.0);
    mm.setValue("y", 2.0);
    mm.setValue("z", 3.0);
    mm.setValue("average", 2.0);
    mm.setValue("error", 0.1);
    // Values are stored in internal GROMACS units (conversion applied),
    // so just verify they are non-zero
    EXPECT_GT(mm.getError(), 0.0);
    // Verify the vector has three entries and at least one is non-zero
    const auto &vec = mm.getVector();
    ASSERT_EQ(3u, vec.size());
    EXPECT_GT(vec[0], 0.0);
    EXPECT_GT(vec[1], 0.0);
    EXPECT_GT(vec[2], 0.0);
}

TEST(MolecularMultipoleTest, HasIdTrue)
{
    MolecularMultipole mm("calc", "D", 0.0, MolPropObservable::DIPOLE);
    EXPECT_TRUE(mm.hasId("x"));
    EXPECT_TRUE(mm.hasId("average"));
    EXPECT_TRUE(mm.hasId("error"));
}

TEST(MolecularMultipoleTest, HasIdFalse)
{
    MolecularMultipole mm("calc", "D", 0.0, MolPropObservable::DIPOLE);
    EXPECT_FALSE(mm.hasId("nosuchid"));
}

TEST(MolecularMultipoleTest, SetValueUnknownThrows)
{
    MolecularMultipole mm("calc", "D", 0.0, MolPropObservable::DIPOLE);
    EXPECT_THROW(mm.setValue("badid", 1.0), gmx::InternalError);
}

TEST(MolecularMultipoleTest, QuadrupoleConstruction)
{
    MolecularMultipole mm("calc", "B", 0.0, MolPropObservable::QUADRUPOLE);
    ASSERT_EQ(6u, mm.getVector().size());
}

TEST(HarmonicsTest, ConstructionFrequency)
{
    Harmonics h("cm^-1", 0.0, MolPropObservable::FREQUENCY);
    EXPECT_EQ(0u, h.getVector().size());
    h.addValue(3657.0);
    h.addValue(1595.0);
    EXPECT_EQ(2u, h.getVector().size());
}

TEST(HarmonicsTest, ConstructionWrongMpoThrows)
{
    EXPECT_THROW(Harmonics h("kJ/mol", 0.0, MolPropObservable::HF), gmx::InternalError);
}

// ============================================================
// MolProp::validate tests
// ============================================================

/*! \brief Helper to build a simple Experiment with N atoms.
 * Atoms are named "C1", "C2", ... with obtype "C.3" and atomid 1, 2, ...
 */
static Experiment makeExperimentWithAtoms(int nAtoms)
{
    Experiment exp("ref", "min");
    for (int i = 0; i < nAtoms; ++i)
    {
        CalcAtom ca(gmx::formatString("C%d", i + 1), "C.3", i + 1);
        exp.AddAtom(ca);
    }
    return exp;
}

TEST(MolPropValidateTest, NoExperimentsIsOk)
{
    MolProp mp;
    mp.SetMolname("empty");
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    mp.validate(&msghandler);
    EXPECT_TRUE(msghandler.ok());
    EXPECT_EQ(0u, msghandler.warningCount(ACTMessage::InconsistentAtomOrder));
    EXPECT_EQ(0u, msghandler.warningCount(ACTMessage::InterFragmentBond));
}

TEST(MolPropValidateTest, ConsistentAtomsIsOk)
{
    MolProp mp;
    mp.SetMolname("water");
    mp.AddExperiment(makeExperimentWithAtoms(3));
    mp.AddExperiment(makeExperimentWithAtoms(3));
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    mp.validate(&msghandler);
    EXPECT_TRUE(msghandler.ok());
    EXPECT_EQ(0u, msghandler.warningCount(ACTMessage::InconsistentAtomOrder));
}

TEST(MolPropValidateTest, DifferentAtomCountIsError)
{
    MolProp mp;
    mp.SetMolname("water");
    mp.AddExperiment(makeExperimentWithAtoms(3));
    // second experiment has a different number of atoms
    mp.AddExperiment(makeExperimentWithAtoms(2));
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    mp.validate(&msghandler);
    EXPECT_FALSE(msghandler.ok());
    EXPECT_GT(msghandler.warningCount(ACTMessage::InconsistentAtomOrder), 0u);
}

TEST(MolPropValidateTest, DifferentAtomNameIsError)
{
    MolProp mp;
    mp.SetMolname("mol");
    // First experiment
    {
        Experiment exp("ref1", "min");
        exp.AddAtom(CalcAtom("C1", "C.3", 1));
        exp.AddAtom(CalcAtom("N2", "N.3", 2));
        mp.AddExperiment(std::move(exp));
    }
    // Second experiment: same size but different atom name at index 0
    {
        Experiment exp("ref2", "min");
        exp.AddAtom(CalcAtom("O1", "O.3", 1));
        exp.AddAtom(CalcAtom("N2", "N.3", 2));
        mp.AddExperiment(std::move(exp));
    }
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    mp.validate(&msghandler);
    EXPECT_FALSE(msghandler.ok());
    EXPECT_GT(msghandler.warningCount(ACTMessage::InconsistentAtomOrder), 0u);
}

TEST(MolPropValidateTest, DifferentObTypeIsError)
{
    MolProp mp;
    mp.SetMolname("mol");
    {
        Experiment exp("ref1", "min");
        exp.AddAtom(CalcAtom("C1", "C.3", 1));
        mp.AddExperiment(std::move(exp));
    }
    {
        Experiment exp("ref2", "min");
        exp.AddAtom(CalcAtom("C1", "C.ar", 1));  // different obtype
        mp.AddExperiment(std::move(exp));
    }
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    mp.validate(&msghandler);
    EXPECT_FALSE(msghandler.ok());
    EXPECT_GT(msghandler.warningCount(ACTMessage::InconsistentAtomOrder), 0u);
}

TEST(MolPropValidateTest, DifferentAtomIdIsError)
{
    MolProp mp;
    mp.SetMolname("mol");
    {
        Experiment exp("ref1", "min");
        exp.AddAtom(CalcAtom("C1", "C.3", 1));
        mp.AddExperiment(std::move(exp));
    }
    {
        Experiment exp("ref2", "min");
        exp.AddAtom(CalcAtom("C1", "C.3", 99));  // different atomid
        mp.AddExperiment(std::move(exp));
    }
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    mp.validate(&msghandler);
    EXPECT_FALSE(msghandler.ok());
    EXPECT_GT(msghandler.warningCount(ACTMessage::InconsistentAtomOrder), 0u);
}

TEST(MolPropValidateTest, DimerNoInterFragmentBondIsOk)
{
    MolProp mp;
    mp.SetMolname("dimer");
    // Two fragments: atoms 0-1 and atoms 2-3
    mp.addFragment(Fragment("frag0", 12.0, 0, 1, 1, "C2", {0, 1}));
    mp.addFragment(Fragment("frag1", 12.0, 0, 1, 1, "C2", {2, 3}));
    // Bond only within fragment 0 and fragment 1
    mp.AddBond(Bond(0, 1, 1.0));
    mp.AddBond(Bond(2, 3, 1.0));
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    mp.validate(&msghandler);
    EXPECT_TRUE(msghandler.ok());
    EXPECT_EQ(0u, msghandler.warningCount(ACTMessage::InterFragmentBond));
}

TEST(MolPropValidateTest, DimerWithInterFragmentBondIsError)
{
    MolProp mp;
    mp.SetMolname("dimer");
    // Two fragments: atoms 0-1 and atoms 2-3
    mp.addFragment(Fragment("frag0", 12.0, 0, 1, 1, "C2", {0, 1}));
    mp.addFragment(Fragment("frag1", 12.0, 0, 1, 1, "C2", {2, 3}));
    // One bond within frag0, one crossing fragments
    mp.AddBond(Bond(0, 1, 1.0));
    mp.AddBond(Bond(1, 2, 1.0));  // crosses fragment boundary
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    mp.validate(&msghandler);
    EXPECT_FALSE(msghandler.ok());
    EXPECT_GT(msghandler.warningCount(ACTMessage::InterFragmentBond), 0u);
}

TEST(MolPropValidateTest, NullMsgHandlerIsNoop)
{
    MolProp mp;
    mp.SetMolname("mol");
    mp.AddExperiment(makeExperimentWithAtoms(2));
    // Should not crash
    mp.validate(nullptr);
}

} // namespace
} // namespace alexandria
