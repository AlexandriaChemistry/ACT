********************
Sensitivity Analysis
********************


The sensitivity analysis is a post-optimization tool that evaluates how sensitive
the total chi-squared (:math:`\chi^2`) deviation is to small changes in each
individual force field parameter. It helps identify which parameters have the most
influence on the fitting quality and whether the optimization has converged to a
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
names of the parameters to analyze.

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

using least-squares regression. The minimum of this parabola is located at

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
parameter set is in a minimum if all parameters are converged. During training
with the HYBRID optimizer, the *-max_local_minimum_generations* option stops the
training when the best parameter set has been in a minimum for that number of
generations.

=============
Output Format
=============

The results are written to a JSON file, set with *-sensout* (default
*sensitivity.json*). The log file only gets a summary::

   Sensitivity analysis (1D) of 42 parameters on the Train set, 84 chi2 evaluations, chi2_0 = 0.482
   Converged: 41/42 parameters, not in a minimum.
   Details of the sensitivity analysis are in sensitivity.json

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

* **mode** — the analysis that was done, *1D* or *2D*.
* **chi2_0** — the :math:`\chi^2` of the current (trained) parameter set before
  any perturbation.
* **n_evaluations** — the number of :math:`\chi^2` computations, not counting
  *chi2_0*.
* **minimum** — whether all parameters are converged, see above.
* **name**, **particle**, **parameter** — the force field parameter.
* **lower**, **upper**, **delta** — its bounds and the step
  :math:`\delta p`.
* **points** — the parameter values at which :math:`\chi^2` was computed, in
  increasing order, with their :math:`\chi^2`.
* **a, b, c** — coefficients of the fitted parabola
  :math:`\chi^2(p) = a p^2 + b p + c`. The fit is done in coordinates
  centred on the parameter and scaled by the step, such that it stays
  accurate for parameters with narrow bounds. **fit_ok** is false if the
  parabola could not be fitted.
* **p_min**, **chi2_min** — the extremum of the parabola, *null* if the
  parabola is a straight line. **p_min_bounded** is *p_min* limited to the
  bounds.
* **converged** — whether the parameter is converged, see above.
* **force_constant** — the dimensionless force constant :math:`a/p^2`, used for
  adaptive step sizes.

Numbers are written with full precision and flags as *true* or *false*.

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
