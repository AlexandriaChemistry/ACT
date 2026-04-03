# forces

This directory implements the computation of energies and forces from a `Topology` and a set of atomic coordinates.

## Contents

- **forcecomputer.h / forcecomputer.cpp** — `ForceComputer`: the public interface for energy and force evaluation. Accepts a `Topology`, coordinate array, and a pre-sized force vector; dispatches to the implementation layer.
- **forcecomputerimpl.h / forcecomputerimpl.cpp** — `ForceComputerImpl`: internal implementation that iterates over all bonded and non-bonded interactions and accumulates contributions to energies and forces.
- **forcecomputerutils.h / forcecomputerutils.cpp** — Shared mathematical helpers used by `ForceComputerImpl` (distance and angle calculations, potential-energy derivatives, etc.).
- **vsitehandler.h / vsitehandler.cpp** — `VsiteHandler`: constructs virtual-site positions from the positions of the real atoms that define them, and distributes virtual-site forces back onto those real atoms.

**Important:** The `forces` vector passed to `ForceComputer::compute` / `computeOnce` must be pre-sized to `topology->atoms().size()` by the caller; the force computer does not resize it.
