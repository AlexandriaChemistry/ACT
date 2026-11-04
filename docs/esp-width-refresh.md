# Keep ESP widths current during training and output

Training only a particle type's `zetaVS` used to change the trial force-field
parameter without refreshing the width cached by the ESP evaluator. The refresh
was conditional on the literal fit name `zeta`. Thus an unrelated selection
flag determined whether the ESP objective evaluated the trial cloud width.

`EspDevComputer` now refreshes ordinary and second-charge widths immediately
before calculating the potential. Its constructor no longer takes a fit-name
flag. This copies particle data; it does not add optimization variables, change
mutability, or make `zeta` an alias for `zetaVS`.

Final ESP statistics and cube/density generation also refresh widths before
copying/evaluating their RESP objects. Previously a SAPT-only fit could print
ESP statistics for cached initial widths rather than its saved force field.
Existing coordinate, charge-generation and dielectric conventions are retained.

The OpenMM wrapper recognizes the Slater metadata but implements distributed
screening only for Gaussian clouds. It now raises `NotImplementedError` for
Slater after detecting the potential and before creating a simulation with
point-charge screening. This is an explicit rejection, not Slater support.

Regression coverage:

- Change an existing molecule's cloud width from 10 to 20 and back. Compare
  pointwise ESP values to a newly initialized RESP object. Cover Gaussian and
  Slater rows 1/2, with retained virtual sites and converted second charges.
- Deliberately restore stale cached widths, write a potential cube, and compare
  its values with fresh evaluation at the written coordinates.
- Detect real Slater OpenMM metadata and reject it before altering the native
  force's direct-space setting. Existing integration tests exercise supported
  potentials.

On the reported water data, external command-line regressions also compare
`-fit "zetaVS"` and `-fit "zeta zetaVS"` in serial, without explicit SIMD, on
four population workers and on two workers with helpers. They select the same
single parameter and give identical paired outcomes. A mutable ordinary zeta
remains excluded when only `zetaVS` is requested.

SAPT-only training of 3P row-1, 3P row-2 and 6P row-1 water models now prints
ESP RMS values 13.2317, 13.4963 and 4.0494 kJ/mol/e respectively, matching fresh
saved-XML evaluation and the pointwise correlation output. Each fit's 1908
Coulomb energies also match the reload at exported precision.

The small residuals in an earlier combined-geometry monomer-swap diagnostic
were investigated independently. The two input waters have slightly different
internal geometries and therefore slightly different geometry-dependent ACM
charges. Evaluating both permutations in one molecule object retained charges
from its first geometry. Separate initialization of the original and swapped
configurations gives identical charges under permutation and zero printed
Coulomb-energy difference for the checked 4P/6P Slater cases and Gaussian
control. No speculative Coulomb-kernel change is included.

Input Slater rows remain physical model choices and are preserved in conversion.
The source fix does not turn row-2 clouds into row 1, rescale widths, or change
the treatment of legacy shared zetatypes.
