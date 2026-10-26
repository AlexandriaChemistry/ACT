/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
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
 * Tests for the alexandria nma tool and the thermochemistry.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include "act/alexandria/alex_modules.h"
#include "act/basics/atomization_energy.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/molprop/molprop_xml.h"
#include "act/properties/thermochemistry.h"
#include "act/topology/actmol.h"
#include "act/utility/units.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

// Physical constants (CODATA 2018), independent of the ACT/GROMACS ones.
//! Gas constant (J/mol K)
const double cR      = 8.314462618;
//! Second radiation constant h c / k (cm K)
const double cC2     = 1.438776877;
//! Energy of one wavenumber (kJ/mol per cm^-1)
const double cWaveKJ = 0.01196265663;
//! Boltzmann constant (J/K)
const double cKB     = 1.380649e-23;
//! Planck constant (J s)
const double cH      = 6.62607015e-34;
//! Atomic mass unit (kg)
const double cAmu    = 1.66053906660e-27;
//! Room temperature used by the nma tool (K)
const double cRoomT  = 298.15;

//! Atomic masses (Dalton) by element
const std::map<std::string, double> elementMass = {
    { "H", 1.008 }, { "C", 12.011 }, { "N", 14.007 },
    { "O", 15.999 }, { "S", 32.06 }, { "Cl", 35.45 }
};

//! An atom from a pdb file
struct PdbAtom
{
    //! Element
    std::string element;
    //! Coordinates (Angstrom)
    double      x[3];
};

//! Read atoms from a pdb file
std::vector<PdbAtom> readPdb(const std::string &filename)
{
    std::vector<PdbAtom> atoms;
    std::ifstream        in(filename);
    EXPECT_TRUE(in.good()) << "Cannot open " << filename;
    std::string          line;
    while (std::getline(in, line))
    {
        if (line.compare(0, 4, "ATOM") == 0 || line.compare(0, 6, "HETATM") == 0)
        {
            PdbAtom atom;
            for(int m = 0; m < 3; m++)
            {
                atom.x[m] = std::stod(line.substr(30+8*m, 8));
            }
            std::string elem = line.size() >= 78 ? line.substr(76, 2) : "";
            elem.erase(0, elem.find_first_not_of(' '));
            atom.element = elem;
            atoms.push_back(atom);
        }
    }
    return atoms;
}

//! \return the value of a JSON element written by JsonTree::addValueUnit
double jsonValue(const nlohmann::json &obj)
{
    return std::stod(obj.at("value").get<std::string>());
}

//! \return the values of an object with keys "1", "2", ... in numerical order
std::vector<double> jsonSeries(const nlohmann::json &obj)
{
    std::vector<double> vals(obj.size());
    for(const auto &item : obj.items())
    {
        size_t k = std::stoul(item.key());
        EXPECT_TRUE(k >= 1 && k <= vals.size()) << item.key();
        if (k >= 1 && k <= vals.size())
        {
            vals[k-1] = jsonValue(item.value());
        }
    }
    return vals;
}

//! Harmonic oscillator thermochemistry for frequencies in cm^-1
struct Vibration
{
    //! Internal energy (kJ/mol)
    double E  = 0;
    //! Heat capacity (J/mol K)
    double cv = 0;
    //! Entropy (J/mol K)
    double S  = 0;
};

//! \return the vibrational contributions at temperature T
Vibration harmonic(const std::vector<double> &freq, double T)
{
    Vibration vib;
    for(auto nu : freq)
    {
        if (nu <= 0)
        {
            continue;
        }
        double theta = cC2*nu;
        double x     = theta/T;
        double em1   = std::expm1(x);
        vib.E  += cR*theta*(0.5 + 1/em1)/1000;
        vib.cv += cR*x*x*std::exp(x)/(em1*em1);
        vib.S  += cR*(x/em1 - std::log1p(-std::exp(-x)));
    }
    return vib;
}

//! \return the translational entropy (J/mol K) at 1 bar
double sackurTetrode(double massDalton, double T)
{
    double m  = massDalton*cAmu;
    double kT = cKB*T;
    double qt = std::pow(2*M_PI*m*kT/(cH*cH), 1.5)*kT/1e5;
    return cR*(std::log(qt) + 2.5);
}

//! \return the rotational entropy (J/mol K) for a rigid rotor
double rigidRotor(const std::vector<PdbAtom> &atoms, bool linear, double sigma, double T)
{
    // Center of mass
    double mtot   = 0;
    double xcm[3] = { 0, 0, 0 };
    for(const auto &a : atoms)
    {
        double m = elementMass.at(a.element);
        mtot += m;
        for(int d = 0; d < 3; d++)
        {
            xcm[d] += m*a.x[d];
        }
    }
    for(int d = 0; d < 3; d++)
    {
        xcm[d] /= mtot;
    }
    // Inertia tensor (Dalton Angstrom^2)
    double I[3][3] = { { 0 } };
    for(const auto &a : atoms)
    {
        double m = elementMass.at(a.element);
        double r[3], r2 = 0;
        for(int d = 0; d < 3; d++)
        {
            r[d] = a.x[d] - xcm[d];
            r2  += r[d]*r[d];
        }
        for(int d = 0; d < 3; d++)
        {
            for(int e = 0; e < 3; e++)
            {
                I[d][e] += m*((d == e ? r2 : 0) - r[d]*r[e]);
            }
        }
    }
    // Principal moments are the eigenvalues of the symmetric 3x3 tensor
    double p1 = I[0][1]*I[0][1] + I[0][2]*I[0][2] + I[1][2]*I[1][2];
    double q  = (I[0][0] + I[1][1] + I[2][2])/3;
    double p2 = (I[0][0]-q)*(I[0][0]-q) + (I[1][1]-q)*(I[1][1]-q) + (I[2][2]-q)*(I[2][2]-q) + 2*p1;
    double p  = std::sqrt(p2/6);
    double B[3][3];
    for(int d = 0; d < 3; d++)
    {
        for(int e = 0; e < 3; e++)
        {
            B[d][e] = (I[d][e] - (d == e ? q : 0))/p;
        }
    }
    double detB = (B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1]) -
                   B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0]) +
                   B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]));
    double phi  = std::acos(std::max(-1.0, std::min(1.0, detB/2)))/3;
    double eig[3];
    eig[0] = q + 2*p*std::cos(phi);
    eig[2] = q + 2*p*std::cos(phi + 2*M_PI/3);
    eig[1] = 3*q - eig[0] - eig[2];
    // Rotational temperatures (K), inertia in kg m^2
    double factor = cH*cH/(8*M_PI*M_PI*cKB);
    double toSI   = cAmu*1e-20;
    if (linear)
    {
        double Imax = std::max(eig[0], std::max(eig[1], eig[2]))*toSI;
        return cR*(1 + std::log(T*Imax/(sigma*factor)));
    }
    double thetaProd = 1;
    for(int d = 0; d < 3; d++)
    {
        thetaProd *= factor/(eig[d]*toSI);
    }
    return cR*(1.5 + std::log(std::sqrt(M_PI*T*T*T/thetaProd)/sigma));
}

//! Description of a test molecule
struct NmaMolecule
{
    //! File in alexandria/tests/mols
    std::string file;
    //! Number of normal modes, 3N-6 or 3N-5
    size_t      nModes;
    //! Whether the molecule is linear
    bool        linear;
    //! Rotational symmetry number
    double      sigma;
};

//! Changes the working directory and changes it back at the end of the scope
class ScopedWorkingDirectory
{
private:
    std::filesystem::path previous_;
public:
    explicit ScopedWorkingDirectory(const std::string &dir) : previous_(std::filesystem::current_path())
    {
        std::filesystem::current_path(dir);
    }
    ~ScopedWorkingDirectory() { std::filesystem::current_path(previous_); }
};

//! Output of one nma run
struct NmaOutput
{
    //! The JSON output (empty if not produced)
    nlohmann::json json;
    //! Minimized structure
    std::string    confout;
    //! Infrared spectrum
    std::string    spectrum;
};

class NmaTest : public gmx::test::CommandLineTestBase
{
protected:
    //! \return path to an input file in the alexandria tests
    std::string alexandriaFile(const std::string &name)
    {
        return fileManager().getInputFilePath("../../alexandria/tests/" + name);
    }

    /*! \brief Run the nma tool
     * \param[in]  molFile Molecule file in alexandria/tests/mols
     * \param[in]  label   Unique label for the output files
     * \param[in]  extra   Extra command line options
     * \param[out] output  The output
     * \return The exit status of the tool
     */
    int runNma(const std::string              &molFile,
               const std::string              &label,
               const std::vector<std::string> &extra,
               NmaOutput                      *output)
    {
        gmx::test::CommandLine cmd;
        cmd.append("nma");
        cmd.addOption("-ff", alexandriaFile("ACS-g.xml"));
        cmd.addOption("-f", alexandriaFile("mols/" + molFile));
        output->confout  = fileManager().getTemporaryFilePath(label + "-confout.pdb");
        output->spectrum = fileManager().getTemporaryFilePath(label + "-ir.xvg");
        cmd.addOption("-c", output->confout);
        cmd.addOption("-ir", output->spectrum);
        cmd.addOption("-g", fileManager().getTemporaryFilePath(label + ".log"));
        cmd.addOption("-json");
        for(const auto &e : extra)
        {
            cmd.append(e);
        }
        int status;
        {
            // The tool writes nma.json to the working directory
            ScopedWorkingDirectory cwd(fileManager().getOutputTempDirectory());
            std::filesystem::remove("nma.json");
            status = nma(cmd.argc(), cmd.argv());
            std::ifstream in("nma.json");
            if (in.good())
            {
                output->json = nlohmann::json::parse(in);
                in.close();
                std::filesystem::remove("nma.json");
            }
        }
        return status;
    }

    //! Check the thermochemistry against textbook formulas
    void checkThermochemistry(const NmaMolecule          &mol,
                              const nlohmann::json       &thermo,
                              const std::vector<double>  &freq,
                              const std::vector<PdbAtom> &atoms)
    {
        // The JSON has three decimals for values between 1 and 100, otherwise
        // six significant digits (see gmx_ftoa).
        auto   relTol = [](double ref, double tol) { return std::abs(ref)*tol + 6e-4; };
        double mass   = 0;
        for(const auto &a : atoms)
        {
            mass += elementMass.at(a.element);
        }
        // Zero point energy, the same at all temperatures
        double zpe = 0;
        for(auto nu : freq)
        {
            zpe += 0.5*cWaveKJ*nu;
        }
        const auto &rt = thermo.at(gmx::formatString("T=%g", cRoomT));
        const auto &t0 = thermo.at("T=0");
        EXPECT_NEAR(zpe, jsonValue(rt.at("Zero point energy")), relTol(zpe, 1e-5));
        EXPECT_NEAR(zpe, jsonValue(t0.at("Zero point energy")), relTol(zpe, 1e-5));

        // Room temperature
        auto   vib    = harmonic(freq, cRoomT);
        double cvRot  = mol.linear ? cR : 1.5*cR;
        double ERot   = cvRot*cRoomT/1000;
        double STrans = sackurTetrode(mass, cRoomT);
        double SRot   = rigidRotor(atoms, mol.linear, mol.sigma, cRoomT);
        const auto &S  = rt.at("Standard entropy");
        const auto &cv = rt.at("Heat capacity cV");
        const auto &E  = rt.at("Internal energy");
        EXPECT_NEAR(STrans, jsonValue(S.at("Translation")), relTol(STrans, 1e-4));
        EXPECT_NEAR(SRot, jsonValue(S.at("Rotation")), relTol(SRot, 2e-3));
        EXPECT_NEAR(vib.S, jsonValue(S.at("Vibration")), relTol(vib.S, 1e-4));
        EXPECT_NEAR(1.5*cR, jsonValue(cv.at("Translation")), relTol(cR, 1e-5));
        EXPECT_NEAR(cvRot, jsonValue(cv.at("Rotation")), relTol(cR, 1e-5));
        EXPECT_NEAR(vib.cv, jsonValue(cv.at("Vibration")), relTol(vib.cv, 1e-4));
        EXPECT_NEAR(1.5*cR*cRoomT/1000, jsonValue(E.at("Translation")), relTol(cR, 1e-5));
        EXPECT_NEAR(ERot, jsonValue(E.at("Rotation")), relTol(ERot, 1e-5));
        EXPECT_NEAR(vib.E, jsonValue(E.at("Vibration")), relTol(vib.E, 1e-4));
        // Totals are the sum of the components
        for(const auto *prop : { &S, &cv, &E })
        {
            double sum = 0;
            for(const auto *comp : { "Translation", "Rotation", "Vibration" })
            {
                sum += jsonValue(prop->at(comp));
            }
            // Each of the four numbers is rounded in the JSON
            EXPECT_NEAR(sum, jsonValue(prop->at("Total")), relTol(sum, 1e-5) + 1.8e-3);
        }
        // At T = 0 there is no thermal energy and no vibrational entropy
        // or heat capacity.
        for(const auto *comp : { "Translation", "Rotation", "Vibration", "Total" })
        {
            EXPECT_EQ(0, jsonValue(t0.at("Internal energy").at(comp))) << comp;
        }
        EXPECT_EQ(0, jsonValue(t0.at("Standard entropy").at("Vibration")));
        EXPECT_EQ(0, jsonValue(t0.at("Heat capacity cV").at("Vibration")));
    }

    //! Run the nma and check everything for one molecule
    void testNma(const NmaMolecule &mol)
    {
        NmaOutput out;
        ASSERT_EQ(0, runNma(mol.file, "nma", {}, &out));
        ASSERT_FALSE(out.json.empty());
        const auto &sim  = out.json.at("simulate");
        auto        freq = jsonSeries(sim.at("Frequencies").at("Alexandria"));
        auto        inte = jsonSeries(sim.at("Intensities").at("Alexandria"));
        EXPECT_EQ(mol.nModes, freq.size());
        EXPECT_EQ(freq.size(), inte.size());
        for(size_t i = 0; i < freq.size(); i++)
        {
            // A minimized structure has no imaginary modes
            EXPECT_GT(freq[i], 0) << "mode " << i+1;
            EXPECT_GE(inte[i], 0) << "mode " << i+1;
            if (i > 0)
            {
                EXPECT_GE(freq[i], freq[i-1]) << "modes should be sorted";
            }
        }
        auto atoms = readPdb(out.confout);
        ASSERT_FALSE(atoms.empty());
        const auto &thermo = sim.at("Thermochemistry").at("Alexandria");
        EXPECT_EQ(mol.sigma, std::stoi(thermo.at("Symmetry number").get<std::string>()));
        checkThermochemistry(mol, thermo, freq, atoms);

        // Reference data
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        // The JSON contains about six significant digits
        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-5));
        checker.checkSequence(freq.begin(), freq.end(), "Frequencies");
        checker.checkSequence(inte.begin(), inte.end(), "Intensities");
        checker.checkDouble(jsonValue(thermo.at("Epot")), "Epot");
        for(const auto *temp : { "T=0", "T=298.15" })
        {
            const auto &tt = thermo.at(temp);
            auto tchecker  = checker.checkCompound("Thermochemistry", temp);
            tchecker.checkDouble(jsonValue(tt.at("Zero point energy")), "Zero point energy");
            tchecker.checkDouble(jsonValue(tt.at("Delta H formation")), "Delta H formation");
            for(const auto *prop : { "Standard entropy", "Heat capacity cV", "Internal energy" })
            {
                for(const auto *comp : { "Translation", "Rotation", "Vibration", "Total" })
                {
                    tchecker.checkDouble(jsonValue(tt.at(prop).at(comp)),
                                         gmx::formatString("%s %s", prop, comp).c_str());
                }
            }
        }
    }

    //! Read the infrared spectrum
    static std::vector<std::pair<double, double>> readSpectrum(const std::string &filename)
    {
        std::vector<std::pair<double, double>> spec;
        std::ifstream in(filename);
        EXPECT_TRUE(in.good()) << "Cannot open " << filename;
        std::string   line;
        while (std::getline(in, line))
        {
            if (line.empty() || line[0] == '#' || line[0] == '@' || line[0] == '&')
            {
                continue;
            }
            std::istringstream iss(line);
            double f, a;
            if (iss >> f >> a)
            {
                spec.push_back({ f, a });
            }
        }
        return spec;
    }
};

TEST_F(NmaTest, Water)
{
    testNma({ "water-3-oep.log.pdb", 3, false, 2 });
}

TEST_F(NmaTest, CarbonDioxide)
{
    testNma({ "carbon-dioxide.sdf", 4, true, 2 });
}

TEST_F(NmaTest, HydrogenChloride)
{
    testNma({ "hydrogen-chloride.sdf", 1, true, 1 });
}

TEST_F(NmaTest, Methanethiol)
{
    testNma({ "methanethiol.sdf", 12, false, 1 });
}

TEST_F(NmaTest, Acetone)
{
    testNma({ "acetone.sdf", 24, false, 2 });
}

TEST_F(NmaTest, Uracil)
{
    testNma({ "uracil.sdf", 30, false, 1 });
}

TEST_F(NmaTest, WaterDimer)
{
    testNma({ "water_dimer.pdb", 12, false, 1 });
}

// Test: The symmetry number is computed from the structure and the
// rotational entropy includes it.
TEST_F(NmaTest, SymmetryNumber)
{
    for(const auto &mol : { NmaMolecule({ "water-3-oep.log.pdb", 3, false, 2 }),
                            NmaMolecule({ "carbon-dioxide.sdf", 4, true, 2 }),
                            NmaMolecule({ "ammonia.sdf", 6, false, 3 }),
                            NmaMolecule({ "chloroform.sdf", 9, false, 3 }),
                            NmaMolecule({ "dichloromethane.sdf", 9, false, 2 }),
                            NmaMolecule({ "1,3,5-triazine.sdf", 21, false, 6 }) })
    {
        SCOPED_TRACE(mol.file);
        NmaOutput out;
        ASSERT_EQ(0, runNma(mol.file, "sigma", {}, &out));
        auto atoms = readPdb(out.confout);
        double SRot = rigidRotor(atoms, mol.linear, mol.sigma, cRoomT);
        const auto &thermo = out.json.at("simulate").at("Thermochemistry").at("Alexandria");
        EXPECT_EQ(mol.sigma, std::stoi(thermo.at("Symmetry number").get<std::string>()));
        const auto &S = thermo.at(gmx::formatString("T=%g", cRoomT)).at("Standard entropy");
        EXPECT_NEAR(SRot, jsonValue(S.at("Rotation")), 2e-3*SRot + 6e-4);
    }
}

// Test: The infrared spectrum has its peak at the most intense mode
TEST_F(NmaTest, InfraredSpectrum)
{
    for(const auto *molFile : { "water-3-oep.log.pdb", "acetone.sdf" })
    {
        SCOPED_TRACE(molFile);
        std::map<double, double> maxIntensity;
        for(double linewidth : { 24.0, 6.0 })
        {
            NmaOutput out;
            ASSERT_EQ(0, runNma(molFile, gmx::formatString("ir%g", linewidth),
                                { "-linewidth", gmx::formatString("%g", linewidth) }, &out));
            const auto &sim  = out.json.at("simulate");
            auto        freq = jsonSeries(sim.at("Frequencies").at("Alexandria"));
            auto        inte = jsonSeries(sim.at("Intensities").at("Alexandria"));
            auto        spec = readSpectrum(out.spectrum);
            ASSERT_FALSE(spec.empty());
            // The spectrum covers all modes
            EXPECT_LE(spec.front().first, freq.front());
            EXPECT_GE(spec.back().first, freq.back());
            // The peak is close to a mode with a significant intensity
            auto   peak    = std::max_element(spec.begin(), spec.end(),
                                              [](const auto &a, const auto &b) { return a.second < b.second; });
            double imax    = *std::max_element(inte.begin(), inte.end());
            bool   nearby  = false;
            for(size_t i = 0; i < freq.size(); i++)
            {
                if (inte[i] > 0.1*imax && std::abs(peak->first - freq[i]) <= linewidth)
                {
                    nearby = true;
                }
            }
            EXPECT_TRUE(nearby) << "Peak at " << peak->first;
            maxIntensity[linewidth] = peak->second;
        }
        // A narrower line gives a higher peak
        EXPECT_GT(maxIntensity[6.0], maxIntensity[24.0]);
    }
}

// Test: The minimized structure is written and differs from the input
TEST_F(NmaTest, ConfOut)
{
    std::string molFile("water-3-oep.log.pdb");
    NmaOutput   out;
    ASSERT_EQ(0, runNma(molFile, "confout", {}, &out));
    auto xin  = readPdb(alexandriaFile("mols/" + molFile));
    auto xout = readPdb(out.confout);
    ASSERT_EQ(xin.size(), xout.size());
    double maxMove = 0;
    for(size_t i = 0; i < xin.size(); i++)
    {
        for(int m = 0; m < 3; m++)
        {
            maxMove = std::max(maxMove, std::abs(xin[i].x[m] - xout[i].x[m]));
        }
    }
    EXPECT_GT(maxMove, 1e-3);
}

// Test: The normal mode analysis can be run without minimization
TEST_F(NmaTest, NoMinimize)
{
    NmaOutput out;
    EXPECT_EQ(0, runNma("water-3-oep.log.pdb", "nomin", { "-nominimize" }, &out));
    ASSERT_FALSE(out.json.empty());
    auto freq = jsonSeries(out.json.at("simulate").at("Frequencies").at("Alexandria"));
    EXPECT_EQ(3u, freq.size());
}

// Test: A molecule file that does not exist
TEST_F(NmaTest, MissingMolecule)
{
    NmaOutput out;
    EXPECT_EQ(1, runNma("does-not-exist.sdf", "missing", {}, &out));
    EXPECT_TRUE(out.json.empty());
}

// Test fixture for the ThermoChemistry class
class ThermoChemistryTest : public gmx::test::CommandLineTestBase
{
protected:
    MsgHandler        msghandler_;
    ACTMol            water_;
    AtomizationEnergy atomenergy_;

    void SetUp() override
    {
        auto pd = getForceField("../../alexandria/tests/ACS-g2");
        std::vector<MolProp> mps;
        MolPropRead(&msghandler_, fileManager().getInputFilePath("../../molprop/tests/molprop.xml").c_str(), &mps);
        for(auto &mp : mps)
        {
            if (mp.getMolname() == "water")
            {
                water_.Merge(&mp);
            }
        }
        water_.GenerateTopology(&msghandler_, pd, missingParameters::Ignore);
        ASSERT_TRUE(msghandler_.ok());
        atomenergy_.read(&msghandler_);
    }

    //! \return the thermochemistry for frequencies in cm^-1
    ThermoChemistry thermo(const std::vector<double> &freqCm,
                           double                     T,
                           double                     P     = 1,
                           double                     scale = 1)
    {
        std::vector<double> freq;
        for(auto f : freqCm)
        {
            freq.push_back(convertToGromacs(f, "cm^-1"));
        }
        return ThermoChemistry(&msghandler_, &water_, water_.xOriginal(), atomenergy_,
                               freq, 0.0, T, P, scale);
    }
};

// Test: The zero point energy scales with the frequency scale factor
TEST_F(ThermoChemistryTest, ScaleFactor)
{
    std::vector<double> freq = { 1600, 3700, 3800 };
    auto tc1 = thermo(freq, cRoomT);
    auto tc2 = thermo(freq, cRoomT, 1, 0.95);
    EXPECT_NEAR(0.5*cWaveKJ*(1600+3700+3800), tc1.ZPE(), 1e-3);
    EXPECT_NEAR(0.95*tc1.ZPE(), tc2.ZPE(), 1e-6*tc1.ZPE());
    // Scaled frequencies give the vibrational properties of the scaled frequencies
    auto tc3 = thermo({ 0.95*1600, 0.95*3700, 0.95*3800 }, cRoomT);
    EXPECT_NEAR(tc3.S0(TCComponent::Vibration), tc2.S0(TCComponent::Vibration), 1e-6);
    EXPECT_NEAR(tc3.cv(TCComponent::Vibration), tc2.cv(TCComponent::Vibration), 1e-6);
}

// Test: Modes with zero or imaginary (negative) frequencies do not contribute
TEST_F(ThermoChemistryTest, NonPositiveFrequencies)
{
    auto tc1 = thermo({ 1600, 3700, 3800 }, cRoomT);
    auto tc2 = thermo({ -200, 0, 1600, 3700, 3800 }, cRoomT);
    for(const auto tcc : { TCComponent::Vibration, TCComponent::Total })
    {
        EXPECT_DOUBLE_EQ(tc1.S0(tcc), tc2.S0(tcc));
        EXPECT_DOUBLE_EQ(tc1.cv(tcc), tc2.cv(tcc));
        EXPECT_DOUBLE_EQ(tc1.Einternal(tcc), tc2.Einternal(tcc));
    }
}

// Test: Imaginary (negative) frequencies do not contribute to the zero point energy
TEST_F(ThermoChemistryTest, ZeroPointEnergyNonPositiveFrequencies)
{
    auto tc1 = thermo({ 1600, 3700, 3800 }, cRoomT);
    auto tc2 = thermo({ -200, 0, 1600, 3700, 3800 }, cRoomT);
    EXPECT_DOUBLE_EQ(tc1.ZPE(), tc2.ZPE());
}

// Test: At high temperature each mode contributes R to the heat capacity
TEST_F(ThermoChemistryTest, ClassicalLimit)
{
    auto tc = thermo({ 100, 200, 300 }, 1e5);
    EXPECT_NEAR(3*cR, tc.cv(TCComponent::Vibration), 1e-3*cR);
}

// Test: The translational entropy decreases with R ln(P)
TEST_F(ThermoChemistryTest, Pressure)
{
    std::vector<double> freq = { 1600, 3700, 3800 };
    auto tc1  = thermo(freq, cRoomT, 1);
    auto tc10 = thermo(freq, cRoomT, 10);
    EXPECT_NEAR(tc1.S0(TCComponent::Translation) - cR*std::log(10.0),
                tc10.S0(TCComponent::Translation), 1e-6);
    // Other components do not depend on pressure
    EXPECT_DOUBLE_EQ(tc1.S0(TCComponent::Rotation), tc10.S0(TCComponent::Rotation));
    EXPECT_DOUBLE_EQ(tc1.S0(TCComponent::Vibration), tc10.S0(TCComponent::Vibration));
}

// Test: Without frequencies there are no vibrational contributions
TEST_F(ThermoChemistryTest, NoFrequencies)
{
    auto tc = thermo({}, cRoomT);
    EXPECT_EQ(0, tc.ZPE());
    EXPECT_EQ(0, tc.S0(TCComponent::Vibration));
    EXPECT_EQ(0, tc.cv(TCComponent::Vibration));
    EXPECT_EQ(0, tc.Einternal(TCComponent::Vibration));
}

} // namespace

} // namespace alexandria
