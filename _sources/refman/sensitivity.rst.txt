********************
Sensitivity Analysis
********************


The sensitivity analysis is a post-optimization tool that evaluates how sensitive
the total chi-squared (:math:`\chi^2`) deviation is to small changes in each
individual force field parameter. It helps identify which parameters have the most
influence on the fitting quality and whether the optimization has converged to a
well-defined minimum.

The sensitivity analysis runs automatically after parameter optimization because the
*-sensitivity* flag defaults to true. It can be disabled with *-nosensitivity*. The
analysis can also be run on an existing force field without re-training by combining
it with the *-nooptimize* and *-norandom_init* flags.

.. note::

   The sensitivity analysis is only available with the MCMC and HYBRID optimizers.
   It is not supported with the Genetic Algorithm (*-optimizer GA*).

============================
Running Sensitivity Analysis
============================

The sensitivity analysis runs automatically after training with MCMC or HYBRID. To
disable it, pass the *-nosensitivity* flag::

   alexandria train_ff -nosensitivity -optimizer MCMC [other flags...]

To analyze an existing force field without re-training, use::

   alexandria train_ff -nooptimize -norandom_init \
       -fc_inter -fit 'param1 param2 param3' [other flags...]

Here, *-fc_inter* selects intermolecular energy as the fitting target (alternatively
*-fc_epot* for intramolecular energy), and *-fit* specifies the space-separated
names of the parameters to analyze.

=========
Algorithm
=========

For each trainable force field parameter :math:`p_i`, the algorithm evaluates
:math:`\chi^2` at three nearby points:

* :math:`p_\mathrm{low} = \max\!\left(p_i - \delta p,\; p_\mathrm{lower}\right)` — below the current value
* :math:`p_0 = \tfrac{1}{2}\!\left(p_\mathrm{low} + p_\mathrm{high}\right)` — midpoint
* :math:`p_\mathrm{high} = \min\!\left(p_i + \delta p,\; p_\mathrm{upper}\right)` — above the current value

where :math:`\delta p = (p_\mathrm{upper} - p_\mathrm{lower}) / 200` (approximately
0.5% of the allowed range on each side of the current value) and
:math:`p_\mathrm{lower}`, :math:`p_\mathrm{upper}` are the allowed bounds for
the parameter.

The three :math:`(p,\,\chi^2)` pairs are then fit to a parabola

.. math::

   \chi^2(p) = a \cdot p^2 + b \cdot p + c

using least-squares regression. The minimum of this parabola is located at

.. math::

   p_\mathrm{opt} = -\frac{b}{2a}, \qquad
   \chi^2_\mathrm{min} = a\, p_\mathrm{opt}^2 + b\, p_\mathrm{opt} + c.

All other parameters are kept at their trained values while each parameter is
perturbed in turn, so the analysis captures the *local* sensitivity around the
training solution. The analysis is performed on the training set only.

=============
Output Format
=============

The sensitivity analysis writes its results to the log file. The output has the
following structure:

.. code-block:: none

   Starting sensitivity analysis. chi2_0 = 0.482  nParam = 42

   Sensitivity epsilon_OW Fit to parabola: a  1.35e+00 b -1.36e+01 c  3.44e+01
       p[0] 4.900   chi2[0] 0.496
       p[1] 5.000   chi2[1] 0.482
       p[2] 5.100   chi2[2] 0.495
       pmin 5.002  chi2min 0.482 (estimate based on parabola)

   Sensitivity C12_VDW Fit to parabola: a  4.00e-12 b -4.01e-07 c  4.82e-01
       p[0] 4.95e+04  chi2[0] 0.4820
       p[1] 5.00e+04  chi2[1] 0.4820
       p[2] 5.05e+04  chi2[2] 0.4820
       pmin 5.01e+04  chi2min 0.4820 (estimate based on parabola)

   Sensitivity analysis done.

The output fields are:

* **chi2_0** — the :math:`\chi^2` of the current (trained) parameter set before
  any perturbation.
* **nParam** — the number of trainable parameters being analyzed.
* **Parameter name** — the name of the force field parameter (interaction type and
  atom-type label).
* **a, b, c** — coefficients of the fitted parabola
  :math:`\chi^2(p) = a p^2 + b p + c`.
* **p[0], p[1], p[2]** — the three parameter values at which :math:`\chi^2` was
  evaluated.
* **chi2[0], chi2[1], chi2[2]** — the corresponding :math:`\chi^2` values.
* **pmin** — the parameter value at the estimated parabolic minimum.
* **chi2min** — the estimated minimum :math:`\chi^2` from the parabola fit.

==============
Interpretation
==============

The coefficient *a* of the fitted parabola quantifies how sharply :math:`\chi^2`
responds to changes in the parameter:

* **Large** :math:`|a|` — the :math:`\chi^2` surface is steeply curved near the
  current value; the parameter is *sensitive* and must be determined accurately.
  A large *a* also indicates a well-defined minimum.
* **Small** :math:`|a| \approx 0` — the :math:`\chi^2` surface is nearly flat; the
  parameter has little influence on the training targets in the explored range. Such
  parameters may be poorly constrained by the data and are candidates for fixing at
  a reasonable value or removing from the fit.

The estimated minimum (**pmin**, **chi2min**) gives a first-order estimate of the
optimal parameter value and the achievable :math:`\chi^2` from the local curvature.
If *pmin* is far from the trained value, the optimization may not have fully
converged for that parameter, or the training data does not strongly constrain it.
When *a* is very small (insensitive parameter), the estimated *pmin* extrapolates
far outside the sampled range and should be treated with caution.

If the fitted parabola has *negative* curvature (:math:`a < 0`), the current
parameter value sits at a local *maximum* of the chi-squared surface, which suggests
an ill-conditioned or degenerate optimization landscape for that parameter.

The sensitivity analysis is particularly useful for:

* Identifying insensitive (poorly constrained) parameters that can be fixed or
  removed from the fitting.
* Checking whether the optimization has converged near a well-defined minimum.
* Guiding subsequent optimization runs by revealing parameters that may benefit from
  a broader or finer search range.
