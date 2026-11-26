/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/basics/atomization_energy.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/properties/pathintegralmd.h"
#include "act/properties/thermochemistry.h"
#include "act/topology/actmol.h"
#include "act/topology/topology.h"
#include "act/utility/jsontree.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "pathintegralthermo.h"

namespace alexandria
{

namespace
{

//! The longest that a classical simulation at a lower temperature is, in units of the given length
const double maxLonger = 16;

//! The largest number of steps of a simulation
const double maxSteps  = 2e9;

/*! \brief h nu / 2 k_B T
 * \param[in] frequency   The frequency in 1/ps
 * \param[in] temperature The temperature in K
 * \return the ratio
 */
double halfQuantum(double frequency, double temperature)
{
    return 0.5*PLANCK*frequency/(BOLTZ*temperature);
}

} // namespace

void gaussLegendre(int n, std::vector<double> *nodes, std::vector<double> *weights)
{
    if (n < 1)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The number of points %d of the quadrature should be one or more", n)));
    }
    nodes->assign(n, 0.0);
    weights->assign(n, 0.0);
    for(int i = 0; i < n; i++)
    {
        // A root of the Legendre polynomial between minus one and one, by
        // Newton's method from the approximation of Tricomi
        double x  = std::cos(M_PI*(i + 0.75)/(n + 0.5));
        double dp = 1;
        for(int iter = 0; iter < 100; iter++)
        {
            double p0 = 1, p1 = x;
            for(int k = 2; k <= n; k++)
            {
                double p2 = ((2*k - 1)*x*p1 - (k - 1)*p0)/k;
                p0 = p1;
                p1 = p2;
            }
            dp        = n*(x*p1 - p0)/(x*x - 1);
            double dx = p1/dp;
            x        -= dx;
            if (std::abs(dx) < 1e-15)
            {
                break;
            }
        }
        // The derivative at the root
        {
            double p0 = 1, p1 = x;
            for(int k = 2; k <= n; k++)
            {
                double p2 = ((2*k - 1)*x*p1 - (k - 1)*p0)/k;
                p0 = p1;
                p1 = p2;
            }
            dp = n*(x*p1 - p0)/(x*x - 1);
        }
        // The roots come in decreasing order
        (*nodes)[n - 1 - i]   = 0.5*(x + 1);
        (*weights)[n - 1 - i] = 1.0/((1 - x*x)*dp*dp);
    }
}

double harmonicBeadEnergy(const std::vector<double> &frequencies,
                          double                     temperature,
                          int                        nbeads,
                          double                     planckScale,
                          double                    *heatCapacity)
{
    const double kT     = BOLTZ*temperature;
    double       energy = 0, cv = 0;
    for(auto nu : frequencies)
    {
        if (nu <= 0)
        {
            continue;
        }
        const double x = planckScale*halfQuantum(nu, temperature);
        if (nbeads <= 0)
        {
            // The quantum oscillator, with the series for small x
            energy += (x < 1e-4) ? 1 + x*x/3 : x/std::tanh(x);
            cv     += (x < 1e-4) ? 1 - x*x/3 : (x/std::sinh(x))*(x/std::sinh(x));
        }
        else
        {
            for(int k = 0; k < nbeads; k++)
            {
                const double s = nbeads*std::sin(M_PI*k/nbeads);
                const double r = x*x/(x*x + s*s);
                energy += r;
                cv     += r*(2*r - 1);
            }
        }
    }
    if (heatCapacity)
    {
        *heatCapacity = BOLTZ*cv;
    }
    return kT*energy;
}

double harmonicQuantumFreeEnergy(const std::vector<double> &frequencies,
                                 double                     temperature)
{
    double sum = 0;
    for(auto nu : frequencies)
    {
        if (nu <= 0)
        {
            continue;
        }
        double x = halfQuantum(nu, temperature);
        if (x < 0.01)
        {
            sum += x*x*(1.0/6 - x*x/180);
        }
        else
        {
            // ln(sinh x / x) without overflow
            sum += x + std::log1p(-std::exp(-2*x)) - std::log(2*x);
        }
    }
    return BOLTZ*temperature*sum;
}

double HarmonicReference::classicalEnergy(double temperature) const
{
    return (nvibration() + 0.5*nrotation() + 1.5)*BOLTZ*temperature;
}

double HarmonicReference::classicalEntropy(double temperature) const
{
    double vib = 0;
    for(auto nu : frequencies)
    {
        if (nu > 0)
        {
            vib += 1 - std::log(2*halfQuantum(nu, temperature));
        }
    }
    return entropyTranslation + entropyRotation + RGAS*vib;
}

HarmonicReference harmonicReference(MsgHandler                   *msghandler,
                                    const ForceField             *pd,
                                    const ACTMol                 *actmol,
                                    const ForceComputer          *forceComp,
                                    const std::vector<gmx::RVec> &xmin,
                                    double                        temperature)
{
    HarmonicReference ref;
    ref.natom  = actmol->NAtom();
    ref.linear = ref.natom > 1 && actmol->linearMolecule();
    auto coords = xmin;
    // The energy of the minimum
    {
        std::map<InteractionType, double> energies;
        std::vector<gmx::RVec>            forces(coords.size());
        forceComp->compute(nullptr, pd, actmol->topology(), &coords, &forces, &energies);
        ref.epot = energies[InteractionType::EPOT];
    }
    std::vector<double> frequencies;
    if (ref.natom > 1)
    {
        const MolHandler         molhandler;
        std::vector<double>      intensities;
        std::vector<std::string> output;
        molhandler.nma(pd, actmol, forceComp, &coords, &frequencies, &intensities, &output);
        for(const auto &line : output)
        {
            msghandler->msg(ACTStatus::Verbose, line);
        }
        if (static_cast<int>(frequencies.size()) != ref.nvibration())
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("The normal mode analysis gave %zu vibrations where %d are expected for %zu atoms, the thermochemistry cannot be computed",
                                                               frequencies.size(), ref.nvibration(), ref.natom)));
        }
        for(auto nu : frequencies)
        {
            if (!(nu > 0))
            {
                GMX_THROW(gmx::InvalidInputError("The minimized structure is not a minimum of the energy: a vibration has no frequency. The thermochemistry cannot be computed."));
            }
        }
    }
    ref.frequencies = frequencies;
    AtomizationEnergy atomenergy;
    atomenergy.read(msghandler);
    const ThermoChemistry tc(msghandler, actmol, coords, atomenergy, frequencies, ref.epot,
                             temperature, 1.0, 1.0);
    ref.symmetryNumber     = tc.symmetryNumber();
    ref.zeroPointEnergy    = tc.ZPE();
    ref.entropyTranslation = tc.S0(TCComponent::Translation);
    ref.internalEnergy     = tc.Einternal(TCComponent::Translation) + tc.Einternal(TCComponent::Vibration);
    ref.heatCapacity       = tc.cv(TCComponent::Translation) + tc.cv(TCComponent::Vibration);
    ref.entropy            = tc.S0(TCComponent::Translation) + tc.S0(TCComponent::Vibration);
    if (ref.natom > 1)
    {
        // An atom does not rotate
        ref.entropyRotation = tc.S0(TCComponent::Rotation);
        ref.internalEnergy += tc.Einternal(TCComponent::Rotation);
        ref.heatCapacity   += tc.cv(TCComponent::Rotation);
        ref.entropy        += tc.S0(TCComponent::Rotation);
    }
    return ref;
}

PathIntegralThermo computeThermochemistry(const HarmonicReference   &harmonic,
                                          const PathIntegralResults &main,
                                          const ThermoIntegration   *integration)
{
    PathIntegralThermo thermo;
    const double       T  = main.temperature;
    const double       kT = BOLTZ*T;
    thermo.temperature    = T;
    thermo.nbeads         = main.nbeads;
    thermo.harmonic       = harmonic;
    thermo.internalEnergy         = main.etot;
    thermo.internalEnergy.average = main.etot.average - harmonic.epot;
    thermo.enthalpy               = thermo.internalEnergy;
    thermo.enthalpy.average      += kT;
    thermo.heatCapacity           = main.heatCapacity;
    thermo.heatCapacity.average  *= KILO;
    thermo.heatCapacity.error    *= KILO;
    thermo.quantumHarmonic        = harmonicQuantumFreeEnergy(harmonic.frequencies, T);
    const bool quantum            = main.nbeads > 1;
    if (quantum)
    {
        // What harmonic oscillators lack with the beads of the simulation
        double cvAll = 0, cvBeads = 0;
        thermo.beadEnergy = (harmonicBeadEnergy(harmonic.frequencies, T, 0, 1.0, &cvAll) -
                             harmonicBeadEnergy(harmonic.frequencies, T, main.nbeads, 1.0, &cvBeads));
        thermo.beadHeatCapacity        = KILO*(cvAll - cvBeads);
        thermo.internalEnergy.average += thermo.beadEnergy;
        thermo.enthalpy.average       += thermo.beadEnergy;
        thermo.heatCapacity.average   += thermo.beadHeatCapacity;
    }
    if (nullptr == integration)
    {
        return thermo;
    }
    const size_t npoint  = integration->nodes.size();
    if (integration->weights.size() != npoint || integration->classicalEpot.size() != npoint ||
        (quantum && (integration->quantumEkin.size() != npoint || integration->quantumBeads.size() != npoint)))
    {
        GMX_THROW(gmx::InvalidInputError("The simulations for the entropy are not there for all points of the integrals"));
    }
    // The anharmonicity of the classical compound, from the temperature
    // zero where it is harmonic
    double anharmonic = 0, anharmonicVar = 0;
    for(size_t i = 0; i < npoint; i++)
    {
        const double u    = integration->nodes[i];
        const double w    = integration->weights[i]/(u*u);
        const auto  &epot = integration->classicalEpot[i];
        anharmonic       -= w*(epot.average - harmonic.epot - 0.5*harmonic.nvibration()*kT*u);
        anharmonicVar    += w*w*epot.error*epot.error;
    }
    thermo.anharmonic.average = anharmonic;
    thermo.anharmonic.error   = std::sqrt(anharmonicVar);
    thermo.anharmonic.n       = npoint;
    // From the classical to the quantum compound. With one bead the
    // simulation is classical and there is nothing to add.
    double dquantum = 0, quantumVar = 0;
    if (quantum)
    {
        dquantum = thermo.quantumHarmonic;
        const double kclassical = 1.5*harmonic.natom*kT;
        // The kinetic energy of classical oscillators
        const double hclassical = 0.5*harmonicBeadEnergy(harmonic.frequencies, T, 1, 1.0);
        for(size_t i = 0; i < npoint; i++)
        {
            const double g    = integration->nodes[i];
            const double w    = 2*integration->weights[i]/g;
            const auto  &ekin = integration->quantumEkin[i];
            // Half of the energy of a harmonic oscillator is kinetic
            const double hkin = 0.5*harmonicBeadEnergy(harmonic.frequencies, T, integration->quantumBeads[i], g);
            dquantum   += w*(ekin.average - kclassical - (hkin - hclassical));
            quantumVar += w*w*ekin.error*ekin.error;
        }
        thermo.quantum.n = npoint;
    }
    thermo.quantum.average = dquantum;
    thermo.quantum.error   = std::sqrt(quantumVar);
    // The free energy relative to the minimum, without the term of the entropy
    // of the classical rigid rotor and harmonic oscillators
    const double eclassical = harmonic.classicalEnergy(T);
    const double sclassical = harmonic.classicalEntropy(T);
    thermo.haveEntropy      = true;
    thermo.entropy.n        = main.etot.n;
    thermo.entropy.average  = sclassical + KILO*(thermo.internalEnergy.average - eclassical - anharmonic - dquantum)/T;
    thermo.entropy.error    = KILO*std::sqrt(main.etot.error*main.etot.error + anharmonicVar + quantumVar)/T;
    // In the Gibbs energy the internal energy of the simulation cancels
    thermo.gibbsEnergy.n       = main.etot.n;
    thermo.gibbsEnergy.average = kT + eclassical - T*sclassical/KILO + anharmonic + dquantum;
    thermo.gibbsEnergy.error   = std::sqrt(anharmonicVar + quantumVar);
    return thermo;
}

PathIntegralThermo pathIntegralThermochemistry(MsgHandler                         *msghandler,
                                               const ForceField                   *pd,
                                               const ACTMol                       *actmol,
                                               const std::vector<gmx::RVec>       &xmin,
                                               const ForceComputer                *forceComp,
                                               const SimulationConfigHandler      &simConfig,
                                               const HarmonicReference            &harmonic,
                                               const PathIntegralResults          &main,
                                               const std::function<void(int, int)> &progress)
{
    const int nquad = simConfig.nquad();
    if (nquad <= 0 || !main.finished || 0 == main.etot.n)
    {
        return computeThermochemistry(harmonic, main, nullptr);
    }
    ThermoIntegration integration;
    gaussLegendre(nquad, &integration.nodes, &integration.weights);
    const bool                quantum   = simConfig.nbeads() > 1;
    const int                 ntotal    = quantum ? 2*nquad : nquad;
    const std::vector<size_t> atomStart = { 0, xmin.size() };
    const double              T         = simConfig.temperature();
    int                       ndone     = 0;
    const auto seed = [&simConfig](int number)
    {
        // Other random numbers for each simulation, unless they are random anyway
        return simConfig.seed() == 0 ? 0 : simConfig.seed() + 100 + number;
    };
    msghandler->write(gmx::formatString("\nEntropy: %d classical simulation%s at lower temperatures for the anharmonicity%s.",
                                        nquad, nquad == 1 ? "" : "s",
                                        quantum ? gmx::formatString(", and %d with %d beads and a smaller constant of Planck for the quantum effects",
                                                                    nquad, simConfig.nbeads()).c_str() : ""));
    for(int i = 0; i < nquad; i++)
    {
        const double u      = integration.nodes[i];
        const double longer = std::min(maxLonger, 1/u);
        auto         sch    = simConfig;
        sch.setSimulation(static_cast<int>(std::min(maxSteps, longer*simConfig.nsteps())),
                          simConfig.deltat(), u*T, seed(ndone), 0, 0);
        sch.setBeads(1, simConfig.tau(),
                     static_cast<int>(std::min(maxSteps, longer*simConfig.nequil())));
        const auto res = pathIntegralMD(msghandler, pd, actmol->topology(), atomStart, xmin,
                                        forceComp, sch, sch.seed(), PathIntegralFiles());
        integration.classicalEpot.push_back(res.epot);
        ndone++;
        if (progress)
        {
            progress(ndone, ntotal);
        }
    }
    for(int i = 0; quantum && i < nquad; i++)
    {
        // Atoms that are heavier need fewer beads
        const double g      = integration.nodes[i];
        const int    nbeads = std::max(1, static_cast<int>(std::floor(g*simConfig.nbeads())));
        integration.quantumBeads.push_back(nbeads);
        if (1 == nbeads)
        {
            // Classical, of which the kinetic energy is known
            AverageError classical;
            classical.average = 1.5*harmonic.natom*BOLTZ*T;
            integration.quantumEkin.push_back(classical);
        }
        else
        {
            auto sch = simConfig;
            sch.setSimulation(simConfig.nsteps(), simConfig.deltat(), T, seed(ndone), 0, 0);
            sch.setBeads(nbeads, simConfig.tau(), simConfig.nequil());
            msghandler->write(gmx::formatString("\nThe constant of Planck is %.4f times its value.", g));
            const auto res = pathIntegralMD(msghandler, pd, actmol->topology(), atomStart, xmin,
                                            forceComp, sch, sch.seed(), PathIntegralFiles(), nullptr, g);
            integration.quantumEkin.push_back(res.ekin);
        }
        ndone++;
        if (progress)
        {
            progress(ndone, ntotal);
        }
    }
    return computeThermochemistry(harmonic, main, &integration);
}

void reportThermochemistry(MsgHandler               *msghandler,
                           const PathIntegralThermo &thermo,
                           JsonTree                 *jtree)
{
    const auto  &h  = thermo.harmonic;
    const double T  = thermo.temperature;
    const double kT = BOLTZ*T;
    msghandler->write(gmx::formatString("\nThermochemistry at %g K and 1 bar from simulations with %d bead%s, next to the rigid rotor and harmonic oscillators of the normal mode analysis.",
                                        T, thermo.nbeads, thermo.nbeads == 1 ? "" : "s"));
    msghandler->write(gmx::formatString("The energies are relative to the minimum of %.4f kJ/mol. The rotational symmetry number is %d.",
                                        h.epot, h.symmetryNumber));
    if (1 == thermo.nbeads)
    {
        msghandler->write("With one bead the simulation is classical, while the harmonic oscillators are quantum mechanical.");
    }
    JsonTree tree("Thermochemistry");
    tree.addValueUnit("Temperature", T, "K");
    tree.addValueUnit("Pressure", 1.0, "bar");
    tree.addObject("Beads", thermo.nbeads);
    tree.addObject("Symmetry number", h.symmetryNumber);
    tree.addValueUnit("Minimum energy", h.epot, "kJ/mol");
    if (thermo.nbeads > 1)
    {
        msghandler->write(gmx::formatString("The values of the simulation include what harmonic oscillators lack with %d beads: %.3f kJ/mol in the energy and %.3f J/mol K in the heat capacity.",
                                            thermo.nbeads, thermo.beadEnergy, thermo.beadHeatCapacity));
        JsonTree btree("Bead correction");
        btree.addValueUnit("Energy", thermo.beadEnergy, "kJ/mol");
        btree.addValueUnit("Heat capacity", thermo.beadHeatCapacity, "J/mol K");
        tree.addObject(btree);
    }
    msghandler->write(gmx::formatString("  %-28s %12s %12s", "", "Harmonic", "Simulation"));
    const auto add = [&](const std::string &name, const char *unit, double harmonic,
                         const AverageError *simulation)
    {
        const auto label = gmx::formatString("%s (%s)", name.c_str(), unit);
        JsonTree   sub(name);
        sub.addObject("Harmonic", harmonic);
        if (simulation)
        {
            msghandler->write(gmx::formatString("  %-28s %12.3f %12.3f +/- %.3f", label.c_str(), harmonic,
                                                simulation->average, simulation->error));
            sub.addObject("Average", simulation->average);
            sub.addObject("Standard error", simulation->error);
        }
        else
        {
            msghandler->write(gmx::formatString("  %-28s %12.3f", label.c_str(), harmonic));
        }
        sub.addObject("Unit", std::string(unit));
        tree.addObject(sub);
    };
    add("Zero point energy", "kJ/mol", h.zeroPointEnergy, nullptr);
    add("Internal energy", "kJ/mol", h.internalEnergy, &thermo.internalEnergy);
    add("Enthalpy", "kJ/mol", h.internalEnergy + kT, &thermo.enthalpy);
    add("Heat capacity cV", "J/mol K", h.heatCapacity, &thermo.heatCapacity);
    auto cp     = thermo.heatCapacity;
    cp.average += RGAS;
    add("Heat capacity cP", "J/mol K", h.heatCapacity + RGAS, &cp);
    const double gharmonic = h.internalEnergy + kT - T*h.entropy/KILO;
    if (thermo.haveEntropy)
    {
        add("Standard entropy", "J/mol K", h.entropy, &thermo.entropy);
        add("Gibbs energy", "kJ/mol", gharmonic, &thermo.gibbsEnergy);
        msghandler->write("Contributions to the free energy of the simulations:");
        msghandler->write(gmx::formatString("  %-28s %12s %12.4f +/- %.4f", "Anharmonicity (kJ/mol)", "",
                                            thermo.anharmonic.average, thermo.anharmonic.error));
        JsonTree ftree("Free energy");
        JsonTree atree("Anharmonicity");
        atree.addObject("Average", thermo.anharmonic.average);
        atree.addObject("Standard error", thermo.anharmonic.error);
        ftree.addObject(atree);
        if (thermo.nbeads > 1)
        {
            msghandler->write(gmx::formatString("  %-28s %12.4f %12.4f +/- %.4f", "Quantum effects (kJ/mol)",
                                                thermo.quantumHarmonic, thermo.quantum.average, thermo.quantum.error));
            JsonTree qtree("Quantum effects");
            qtree.addObject("Harmonic", thermo.quantumHarmonic);
            qtree.addObject("Average", thermo.quantum.average);
            qtree.addObject("Standard error", thermo.quantum.error);
            ftree.addObject(qtree);
        }
        ftree.addObject("Unit", std::string("kJ/mol"));
        tree.addObject(ftree);
    }
    else
    {
        add("Standard entropy", "J/mol K", h.entropy, nullptr);
        add("Gibbs energy", "kJ/mol", gharmonic, nullptr);
    }
    if (jtree)
    {
        jtree->addObject(tree);
    }
}

} // namespace alexandria
