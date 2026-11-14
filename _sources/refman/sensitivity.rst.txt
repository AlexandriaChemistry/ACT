********************
Sensitivity Analysis
********************


The sensitivity analysis is a post-optimization tool that evaluates how sensitive
the total chi-squared (:math:`\chi^2`) deviation is to small changes in each
individual force field parameter and, optionally, in each pair of parameters. It
helps identify which parameters have the most influence on the fitting quality,
which parameters are correlated, and whether the optimization has converged to a
well-defined minimum.

The sensitivity analysis runs automatically after parameter optimization. The
*-sensitivity* option selects it: *1D* (the default) changes each parameter on its
own, *2D* also changes every pair of parameters, and *No* disables the analysis. The
analysis can also be run on an existing force field without re-training by combining
it with the *-nooptimize* and *-norandom_init* flags.

.. note::

   The sensitivity analysis is only available with the MCMC and HYBRID optimizers.
   It is not supported with the Genetic Algorithm (*-optimizer GA*).

============================
Running Sensitivity Analysis
============================

The sensitivity analysis runs automatically after training with MCMC or HYBRID. To
disable it, pass *-sensitivity No*::

   alexandria train_ff -sensitivity No -optimizer MCMC [other flags...]

To analyze an existing force field without re-training, use::

   alexandria train_ff -nooptimize -norandom_init \
       -fc_inter -fit 'param1 param2 param3' [other flags...]

Here, *-fc_inter* selects intermolecular energy as the fitting target (alternatively
*-fc_epot* for intramolecular energy), and *-fit* specifies the space-separated
names of the parameters to analyze. Add *-sensitivity 2D* for the two-dimensional
analysis, and *-sensout* to choose the name of the output file.

The analysis is done for the parameters as they are. Parameters that are zero in a
force field that is not trained yet, like the Lennard-Jones parameters of a
template force field, may not change :math:`\chi^2` at all; they are reported as
*flat*, see below.

=========
Algorithm
=========

For each trainable force field parameter :math:`p_i`, the algorithm evaluates
:math:`\chi^2` at three nearby points, one of which is the parameter itself:

* :math:`p_i - \delta p,\; p_i,\; p_i + \delta p` if these are within the bounds
* :math:`p_i,\; p_i + \delta p,\; p_i + 2\delta p` if :math:`p_i - \delta p` is below the lower bound
* :math:`p_i - 2\delta p,\; p_i - \delta p,\; p_i` if :math:`p_i + \delta p` is above the upper bound

where :math:`\delta p = (p_\mathrm{upper} - p_\mathrm{lower}) / 200` (0.5% of
the allowed range) and :math:`p_\mathrm{lower}`, :math:`p_\mathrm{upper}` are the
allowed bounds for the parameter. The :math:`\chi^2` of the parameter itself is
that of the current parameter set, so two evaluations are needed per parameter.

The three :math:`(p,\,\chi^2)` pairs are then fit to a parabola

.. math::

   \chi^2(p) = a \cdot p^2 + b \cdot p + c

using least-squares regression. The fit is done in coordinates centred on
:math:`p_i` and scaled by :math:`\delta p`, such that it stays accurate for
parameters with narrow bounds, and converted to :math:`a, b, c` afterwards. A
curvature smaller than the rounding errors of :math:`\chi^2` is set to zero. If
:math:`\chi^2` is the same at all three points, the parameter is *flat*: it does
not change :math:`\chi^2` here, and has neither curvature nor slope. The minimum
of this parabola is located at

.. math::

   p_\mathrm{opt} = -\frac{b}{2a}, \qquad
   \chi^2_\mathrm{min} = a\, p_\mathrm{opt}^2 + b\, p_\mathrm{opt} + c.

All other parameters are kept at their trained values while each parameter is
perturbed in turn, so the analysis captures the *local* sensitivity around the
training solution. The analysis is performed on the training set only.

A parameter is *converged* if

1. the parabola curves upwards, :math:`a > 0`, such that it has a minimum, and
2. the parameter is within 5% of its allowed range of that minimum,
   :math:`|p_i - p_\mathrm{opt}| \le 0.05\,(p_\mathrm{upper} - p_\mathrm{lower})`.

Here :math:`p_\mathrm{opt}` is limited to the bounds, such that a parameter at a
bound is converged if the minimum of the parabola lies beyond that bound. The
parameter set is in a minimum if all parameters are converged.

=========================
Two-dimensional analysis
=========================

Each parameter on its own can be at the minimum of its parabola while the
parameter set is not in a minimum: two strongly correlated parameters form a
valley in which each one is optimal given the other, and a saddle point looks like
a minimum along each axis. With *-sensitivity 2D* every pair of parameters
:math:`(p_i, p_j)` is analysed as well, on the :math:`3 \times 3` grid of the
three values of each parameter from the one-dimensional analysis. Five points of
the grid are known from that analysis, the four corners, where both parameters
differ from their own value, are computed. The surface is described as the sum of
the two parabolas plus a coupling term

.. math::

   \chi^2(x, y) \approx P_i(x) + P_j(y) - \chi^2_0 + a_{ij}\,(x - p_i)(y - p_j)

where :math:`P_i` and :math:`P_j` are the parabolas of the parameters and
:math:`\chi^2_0` is :math:`\chi^2` of the parameter set. The coupling is fitted to
what the parabolas do not describe at the corners,

.. math::

   a_{ij} = \frac{\sum_k r_k\,\Delta x_k \Delta y_k}{\sum_k (\Delta x_k \Delta y_k)^2},
   \qquad r_k = \chi^2_k - \left[P_i(x_k) + P_j(y_k) - \chi^2_0\right]

with :math:`\Delta x_k = x_k - p_i` and :math:`\Delta y_k = y_k - p_j`. On a symmetric
grid this is the mixed central difference, and it is exact for a quadratic
surface, also at bounds. The root mean square of the corners that the surface does
not describe is reported as *corner_rms*; a large value means that :math:`\chi^2` is
not quadratic on the scale of :math:`\delta p`. The Hessian and the gradient of the
pair at the parameter set are

.. math::

   H = \begin{pmatrix} 2a_i & a_{ij} \\ a_{ij} & 2a_j \end{pmatrix}, \qquad
   g = \begin{pmatrix} 2a_i p_i + b_i \\ 2a_j p_j + b_j \end{pmatrix}, \qquad
   \rho_{ij} = \frac{a_{ij}}{2\sqrt{a_i a_j}}

where the dimensionless coupling :math:`\rho_{ij}` is defined if both curvatures are
positive, and is between -1 and 1 for a minimum. The surface of a pair is

* a *minimum* if both curvatures and the determinant of :math:`H` are positive,
  that is :math:`|\rho_{ij}| < 1`,
* a *saddle* if the determinant is negative,
* a *maximum* if both curvatures are negative and the determinant is positive,
* *degenerate* if a curvature is zero, or the determinant is zero compared to its
  terms, for instance for a flat parameter or a valley without a minimum.

A pair is *converged* if its surface is a minimum, and both parameters are within
5% of their range of the minimum of the surface within the bounds. If the
stationary point :math:`p - H^{-1} g` lies outside the bounds, that minimum is on an
edge of the bounds, where the surface is a parabola in the other parameter. In 2D
the parameter set is in a minimum if all parameters and all pairs are converged.
The analysis takes :math:`2N` evaluations of :math:`\chi^2` in 1D and :math:`2N^2`
in 2D for :math:`N` parameters:

.. table:: Number of evaluations of :math:`\chi^2` of the training set

   ================ ======== ======== ========
   Parameters *N*   1D       2D       Pairs
   ================ ======== ======== ========
   10               20       200      45
   42               84       3,528    861
   100              200      20,000   4,950
   ================ ======== ======== ========

They are divided over the master and the middlemen. The Hessian of all parameters
follows from the pairs as well. It is scaled by the ranges of the parameters,
:math:`H_{ij} (p^\mathrm{upper}_i - p^\mathrm{lower}_i)(p^\mathrm{upper}_j - p^\mathrm{lower}_j)`,
such that its eigenvalues do not depend on the units of the parameters, and the
analysis reports whether it is positive definite, its extreme eigenvalues and the
largest component of the Newton step :math:`-H^{-1} g` relative to the range of the
parameter. All pairs can be minima while the Hessian of all parameters is not
positive definite, so this is a stricter test; it is reported but not used for the
convergence.

============================
Convergence during training
============================

During training with the HYBRID optimizer, the *-max_local_minimum_generations*
option stops the training when the best parameter set has been in a minimum for
that number of generations. The test for a minimum is two-dimensional with
*-sensitivity 2D*, and one-dimensional otherwise, also with *-sensitivity No*. It is
done in each generation, so the two-dimensional test costs :math:`2N^2` evaluations
of :math:`\chi^2` per generation. With *-v 3* the log file reports for each
generation how many parameters, and in 2D pairs, are converged. Adaptive step sizes
(*-adaptive_steps*) use the force constants of the one-dimensional analysis.

=============
Output Format
=============

The results are written to a JSON file, set with *-sensout* (default
*sensitivity.json*). The log file only gets a summary::

   Sensitivity analysis (1D) of 42 parameters on the Train set, 84 chi2 evaluations, chi2_0 = 0.482
   Converged: 41/42 parameters, not in a minimum.
   Details of the sensitivity analysis are in sensitivity.json

or, in 2D::

   Sensitivity analysis (2D) of 42 parameters and 861 pairs on the Train set, 3528 chi2 evaluations, chi2_0 = 0.482
   Converged: 41/42 parameters, 850/861 pairs (858 minimum, 3 saddle, 0 maximum, 0 degenerate), not in a minimum.
   The Hessian of all parameters is positive definite.
   Details of the sensitivity analysis are in sensitivity.json

If there are flat parameters, a line says how many.

and the json output of *train_ff* (*-json*) gets the same summary under the key
*sensitivity*, with the name of the file. The file has the following structure,
shown here for one parameter:

.. code-block:: json

   {
     "sensitivity": {
       "mode": "1D",
       "dataset": "Train",
       "chi2_0": 0.482,
       "n_parameters": 42,
       "n_evaluations": 84,
       "tolerance_fraction": 0.05,
       "minimum": false,
       "n_parameters_converged": 41,
       "n_parameters_flat": 0,
       "parameters": [
         {
           "index": 0,
           "name": "OW epsilon",
           "particle": "OW",
           "parameter": "epsilon",
           "lower": 0.0,
           "upper": 20.0,
           "p_orig": 5.0,
           "chi2_orig": 0.482,
           "points": [ { "p": 4.9, "chi2": 0.496 },
                       { "p": 5.0, "chi2": 0.482 },
                       { "p": 5.1, "chi2": 0.495 } ],
           "fit_ok": true,
           "flat": false,
           "a": 1.35, "b": -13.6, "c": 34.4,
           "p_min": 5.002,
           "chi2_min": 0.482,
           "p_min_bounded": 5.002,
           "converged": true,
           "delta": 0.1,
           "force_constant": 0.054
         }
       ]
     }
   }

The fields are:

* **mode** - the analysis that was done, *1D* or *2D*.
* **chi2_0** - the :math:`\chi^2` of the current (trained) parameter set before
  any perturbation.
* **n_evaluations** - the number of :math:`\chi^2` computations, not counting
  *chi2_0*.
* **minimum** - whether all parameters are converged, see above.
* **name**, **particle**, **parameter** - the force field parameter.
* **lower**, **upper**, **delta** - its bounds and the step
  :math:`\delta p`.
* **points** - the parameter values at which :math:`\chi^2` was computed, in
  increasing order, with their :math:`\chi^2`.
* **a, b, c** - coefficients of the fitted parabola
  :math:`\chi^2(p) = a p^2 + b p + c`. The fit is done in coordinates
  centred on the parameter and scaled by the step, such that it stays
  accurate for parameters with narrow bounds. **fit_ok** is false if the
  parabola could not be fitted.
* **p_min**, **chi2_min** - the extremum of the parabola, *null* if the
  parabola is a straight line. **p_min_bounded** is *p_min* limited to the
  bounds.
* **flat** - whether :math:`\chi^2` is the same at all points, and
  **n_parameters_flat** the number of flat parameters.
* **converged** - whether the parameter is converged, see above.
* **force_constant** - the dimensionless force constant :math:`a/p^2`, used for
  adaptive step sizes.

In 2D the file has these as well:

.. code-block:: json

   {
     "sensitivity": {
       "n_pairs": 861,
       "n_pairs_converged": 850,
       "pair_shapes": { "minimum": 858, "saddle": 3, "maximum": 0, "degenerate": 0 },
       "pairs": [
         {
           "i": 0, "j": 1,
           "names": [ "OW epsilon", "OW sigma" ],
           "grid": { "p_i": [ 4.9, 5.0, 5.1 ],
                     "p_j": [ 0.312, 0.315, 0.318 ],
                     "chi2": [ [ 0.509, 0.496, 0.490 ],
                               [ 0.488, 0.482, 0.484 ],
                               [ 0.488, 0.495, 0.509 ] ] },
           "a_ij": 66.1, "rho": 0.91, "corner_rms": 2.1e-6,
           "hessian": [ [ 2.7, 66.1 ], [ 66.1, 1960.0 ] ],
           "gradient": [ -0.05, 0.42 ],
           "shape": "minimum",
           "p_min": [ 5.03, 0.3138 ], "chi2_min": 0.4818,
           "p_min_bounded": [ 5.03, 0.3138 ],
           "converged": true
         }
       ],
       "hessian_diagnostic": {
         "positive_definite": true,
         "eigenvalue_min": 0.21,
         "eigenvalue_max": 5.4e3,
         "max_newton_step_fraction": 0.012
       }
     }
   }

* **grid** - the values of both parameters and :math:`\chi^2` on the grid,
  *chi2[k][l]* at *p_i[k]* and *p_j[l]*.
* **a_ij**, **rho**, **corner_rms**, **hessian**, **gradient**, **shape** - the
  coupling, the dimensionless coupling, the deviation of the corners, the Hessian,
  the gradient and the shape of the pair, see above.
* **p_min**, **chi2_min** - the stationary point of a minimum and its
  :math:`\chi^2`, **p_min_bounded** the minimum within the bounds, all *null* if
  the pair is not a minimum.
* **hessian_diagnostic** - the Hessian of all parameters, see above.

Numbers are written with full precision, flags as *true* or *false* and undefined
values as *null*.

==============
Interpretation
==============

The coefficient *a* of the fitted parabola quantifies how sharply :math:`\chi^2`
responds to changes in the parameter:

* **Large** :math:`|a|` - the :math:`\chi^2` surface is steeply curved near the
  current value; the parameter is *sensitive* and must be determined accurately.
  A large *a* also indicates a well-defined minimum.
* **Small** :math:`|a| \approx 0` - the :math:`\chi^2` surface is nearly flat; the
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

A *flat* parameter does not influence the training targets at all at its current
value. In a trained force field it is not constrained by the data; in a force field
that is not trained yet it is often zero together with a parameter it multiplies,
like :math:`\sigma` with :math:`\epsilon = 0`, or a virtual site without charge. All
pairs with a flat parameter are degenerate.

In 2D, :math:`\rho_{ij}` close to 1 or -1 means that two parameters are strongly
correlated: the data constrain a combination of them better than each one, and the
trained values of each may not be meaningful on their own. A *saddle* means that the
parameter set is not at a minimum even if both parameters are converged on their
own; the training can still lower :math:`\chi^2` by changing both.

The sensitivity analysis is particularly useful for:

* Identifying insensitive (poorly constrained) parameters that can be fixed or
  removed from the fitting.
* Identifying correlated parameters, with the two-dimensional analysis.
* Checking whether the optimization has converged near a well-defined minimum.
* Guiding subsequent optimization runs by revealing parameters that may benefit from
  a broader or finer search range.
