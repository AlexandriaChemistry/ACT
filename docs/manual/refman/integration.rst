.. _sec-integration:

********************************
Shell Relaxation and Integration
********************************
In polarizable force fields the positions of the shells (or Drude particles) are not determined by the input structure, they have to be computed for each structure. This is done by minimizing the energy with respect to the shell positions, keeping the atoms fixed, as described in this section. The second part of the section describes the integration of the equations of motion in *alexandria simulate*.

.. _sec-shellrelax:

=====================
Relaxation of shells
=====================
The shells have no mass. In the adiabatic approximation they follow the atoms instantaneously, that is, for each configuration of the atoms they are at the minimum of the potential energy. The shells are moved to this minimum iteratively, using the forces on them, as follows.

#. The shells are put on top of their core particle and the energy and forces are computed.
#. The mean square force on the :math:`N_s` shells is computed

   .. math:: F_{ms} ~=~ \frac{1}{N_s}\sum_{s=1}^{N_s} \left|{\bf F}_s\right|^2.
      :label: fmsshell

   If :math:`F_{ms}` is at most the tolerance, the shells are converged. The tolerance is given in (kJ/mol nm)\ :sup:`2`, see below.
#. Otherwise each shell is displaced along the force on it,

   .. math:: \Delta{\bf r}_s ~=~ \lambda\,\frac{{\bf F}_s}{k_s},\qquad k_s ~=~ \frac{1}{4\pi\varepsilon_0}\frac{q_s^2}{\alpha_c},
      :label: shellstep

   where :math:`k_s` is the force constant of the core-shell spring, see Eqn. :eq:`vpol`. With :math:`\lambda = 1` this step would bring an isolated shell to its minimum, since the force of the spring then cancels the force :math:`{\bf F}_s` from the other particles. With hyperpolarization, the potential of the core-shell pair has an additional term :math:`k_h\,(r_{cs}-r_h)^4` for distances :math:`r_{cs}` larger than :math:`r_h`. For a shell away from its core, :math:`k_s` in Eqn. :eq:`shellstep` is then replaced by the ratio of the core-shell force and the distance

.. math:: k_s + 4k_h\,\frac{(r_{cs}-r_h)^3}{r_{cs}}\qquad\mbox{if}\quad r_{cs} > r_h.
   :label: kshyper

For a shell on top of its core, the quartic term only contributes if :math:`r_h = 0`, and then the displacement :math:`\Delta r` is the solution of :math:`\left|{\bf F}_s\right| = k_s\,\Delta r + 4k_h\,\Delta r^3`. The force constants are computed from the current polarizabilities and shell charges, such that they are correct during training as well.
#. Since the shells interact with each other, a full step (:math:`\lambda = 1`) can overshoot. If :math:`F_{ms}` increases, the step is reduced to :math:`\lambda = 1/2` for the rest of the relaxation. Smaller steps are not used, since they can freeze the shells far from their minimum.
#. A shell with one core is not allowed to move further away from it than a maximum distance :math:`d_{max}`. If the step brings a shell further away, it is moved back to the sphere with radius :math:`d_{max}` around its core, along the core-shell vector. This prevents the polarization catastrophe, where shells of neighbouring molecules attract each other without limit.
#. The forces are computed again and the procedure is repeated from step 2, until the shells are converged or the maximum number of iterations is reached.

A shell held at :math:`d_{max}` by a force pointing away from its core can never reach a zero force. For such a shell only the component of the force perpendicular to the core-shell vector :math:`{\bf r}_{cs}` is used in Eqn. :eq:`fmsshell`

.. math:: {\bf F}_{s,\perp} ~=~ {\bf F}_s - \frac{{\bf F}_s\cdot{\bf r}_{cs}}{r_{cs}^2}\,{\bf r}_{cs}\qquad\mbox{if}\quad r_{cs} = d_{max}\quad\mbox{and}\quad{\bf F}_s\cdot{\bf r}_{cs} > 0,
   :label: fperp

which is zero at the minimum of the energy on the sphere.

The options that control the relaxation are listed in :numref:`tab-shellrelax`.

.. table:: Options for the shell relaxation.
   :name: tab-shellrelax

   +----------------------+---------------------+--------------------+--------------------------------------+
   | **Quantity**         | **train_ff**        | **simulate**       | **Default**                          |
   +======================+=====================+====================+======================================+
   | Tolerance on         | -shellTolerance     | -shelltoler        | 0.001 (train_ff), 1e-6 (simulate),   |
   | :math:`F_{ms}`       |                     |                    | at most a tenth of -toler in simulate|
   +----------------------+---------------------+--------------------+--------------------------------------+
   | Maximum number of    | -shell_max_iter     | -maxiter           | 25 (train_ff), 100 (simulate)        |
   | iterations           |                     |                    |                                      |
   +----------------------+---------------------+--------------------+--------------------------------------+
   | :math:`d_{max}` (nm) | -shell_max_distance |                    | 0.02 (train_ff), 0.04 (simulate),    |
   |                      |                     |                    | 0 means no maximum                   |
   +----------------------+---------------------+--------------------+--------------------------------------+

----------------------------
Shell relaxation in training
----------------------------
During training the energies of all conformations (structures) of a molecule are computed together, see :ref:`sec-training`. The shells of all conformations are relaxed at the same time, and the conformations for which the shells converged are not computed again in the next iterations. For a dimer, the shells of each monomer are relaxed first in the isolated monomer, then all shells in the dimer. The induction energy is the change of the electrostatic energy upon relaxation of the shells. If the induction correction is trained separately (see :ref:`sec-indcorr`), the second order induction is computed by relaxing the shells of one monomer in the field of the other one, of which the shells stay at their positions in the isolated monomer, and the other way around. The rest of the induction is the induction correction.

Without polarization the interaction energy of a dimer is computed from the intermolecular interactions only, unless forces are needed, since the intramolecular energies of the dimer and the monomers cancel.

.. _sec-integrator:

=======================
Integration in simulate
=======================
The command *alexandria simulate* integrates the equations of motion in the gas phase (see :ref:`sec-simulations`) using the velocity Verlet algorithm :cite:p:`Swope1982a`. For each atom :math:`i` with mass :math:`m_i`, coordinates :math:`{\bf r}_i` and velocity :math:`{\bf v}_i`, a step of length :math:`\Delta t` is

.. math:: {\bf r}_i(t+\Delta t) ~=~ {\bf r}_i(t) + \Delta t\,{\bf v}_i(t) + \frac{\Delta t^2}{2m_i}{\bf F}_i(t),
   :label: vvpos

after which the forces :math:`{\bf F}_i(t+\Delta t)` are computed for the new coordinates and the velocities are updated

.. math:: {\bf v}_i(t+\Delta t) ~=~ {\bf v}_i(t) + \frac{\Delta t}{2m_i}\left[{\bf F}_i(t) + {\bf F}_i(t+\Delta t)\right].
   :label: vvvel

Only particles with a mass are integrated. Virtual sites are constructed from the atoms before the forces are computed, and the forces on them are distributed over the atoms they are constructed from. The shells are relaxed for each new configuration of the atoms as described in :ref:`sec-shellrelax`, starting from their cores. Bond lengths and angles can be constrained, as described in :ref:`sec-constraints`.

With *-temp* larger than zero, the initial velocities are generated from a Maxwell distribution at that temperature, using the random number seed *-seed*. The motion of the center of mass is removed, the velocity components along the constraints are removed, and finally the velocities are scaled such that the temperature, Eqn. :eq:`tempmd`, is exactly the requested one. There is no thermostat, the simulation conserves the total energy, which can be used to validate the forces. The energy file contains the potential energy terms, the kinetic energy :math:`E_{kin} = \frac{1}{2}\sum_i m_i{\bf v}_i^2`, the total energy and the temperature

.. math:: T ~=~ \frac{2E_{kin}}{N_{df}\,k_B},
   :label: tempmd

which are written every *-nstener* steps before the step is taken. The number of degrees of freedom is

.. math:: N_{df} ~=~ 3N - 3 - N_c,
   :label: ndf

for :math:`N` atoms with a mass and :math:`N_c` constraints, see :ref:`sec-constraints`. Shells and virtual sites do not contribute, since they have no mass. Three degrees of freedom are subtracted for the center of mass, which does not move, since there are no external forces and the initial velocities are either zero or have the center of mass motion removed. Since the potential energy takes up part of the initial kinetic energy, the temperature during the simulation is lower than the initial one, about half of it if the starting structure is close to the minimum of the energy. The time step is set with *-deltat* (ps) and the number of steps with *-nsteps*.

.. _sec-constraints:

-----------
Constraints
-----------
With the option *-constrain* the fastest motions can be removed from the simulation, using the RATTLE algorithm :cite:p:`Ryckaert1977a,Andersen1983a`. Each constraint keeps the distance :math:`d_c` between two atoms constant. The option has three choices:

* *no*: there are no constraints (the default).
* *bonds*: the length of each bond is kept at the reference bond length :math:`b_0` of the force field.
* *angles*: in addition to the bonds, for each angle :math:`i-j-k` the distance between the outer atoms :math:`i` and :math:`k` is kept constant. Together with the constrained bonds :math:`i-j` and :math:`j-k` this keeps the angle at its reference value :math:`\theta_0`, since

  .. math:: d_{ik} ~=~ \sqrt{b_{ij}^2 + b_{jk}^2 - 2\,b_{ij}\,b_{jk}\cos\theta_0}.
     :label: dik

Not all angles can be constrained independently. The positions of :math:`n` atoms bonded to the same atom, relative to each other, are determined by :math:`2n-3` angles, for instance five of the six angles around a tetrahedral carbon. Around an atom with three neighbours of which the reference angles add up to more than 355\ :math:`^{\circ}`, only two angles are constrained, since the third one follows from the other two in a planar group. Angles larger than 175\ :math:`^{\circ}` are not constrained, since the distance :math:`d_{ik}` does not determine a (nearly) linear angle. Finally, an angle constraint that depends linearly on the other constraints in the starting coordinates is left out as well. This is the case in rings, where the angles also have to close the ring. The log file lists the number of bond and angle constraints and the number of angles that were left out. The angles that are not constrained still contribute to the potential energy. Note that when the reference values of the angles in a ring are not consistent with a planar ring, the ring has to pucker to satisfy the constraints, which may not be possible from a planar starting structure. In that case the simulation stops with an error, and *-constrain bonds* can be used instead.

In the first part of RATTLE, the coordinates :math:`{\bf r}_i'` from Eqn. :eq:`vvpos` are corrected by constraint forces along the distance vectors :math:`{\bf r}_c(t)` of the constrained pairs at time :math:`t`

.. math:: {\bf r}_i(t+\Delta t) ~=~ {\bf r}_i' + \frac{1}{m_i}\sum_c g_c\,\sigma_{ic}\,{\bf r}_c(t),
   :label: shake

where :math:`\sigma_{ic}` is 1 if atom :math:`i` is the first atom of constraint :math:`c`, -1 if it is the second and 0 otherwise. The multipliers :math:`g_c` follow from the conditions :math:`|{\bf r}_c(t+\Delta t)|^2 = d_c^2`. Since angle constraints are strongly coupled to the bond constraints, these non-linear equations are solved for all constraints together by Newton iterations, as in the matrix SHAKE algorithm :cite:p:`Kraeutler2001a`, rather than one constraint at a time as in the original SHAKE algorithm, which converges very slowly for coupled constraints. The step is halved when it does not reduce the deviations from the constraints. The displacements due to the constraints, divided by :math:`\Delta t`, are added to the velocities. In the second part of RATTLE, after Eqn. :eq:`vvvel`, the velocity components along the constraints are removed by solving the linear equations

.. math:: {\bf r}_c\cdot\left({\bf v}_{c_1} - {\bf v}_{c_2}\right) ~=~ 0
   :label: rattle

for all constraints together, where :math:`c_1` and :math:`c_2` are the atoms of constraint :math:`c`. The starting coordinates are constrained as well. The relative tolerance on the constrained distances and the maximum number of iterations are set with *-constr_toler* and *-constr_maxiter*. Each constraint removes one degree of freedom, such that :math:`N_c` in Eqn. :eq:`ndf` is the number of constraints. For acetone, for instance, there are 9 bond constraints and 12 angle constraints, leaving 6 degrees of freedom: the rotation of the molecule, the rotation of the two methyl groups and the out of plane motion of the carbonyl carbon.

The constraints are applied in energy minimization (option *-minimize*) as well, with all three algorithms that can be selected with *-minalg*. The starting coordinates are constrained first. After each computation of the forces, the components of the forces along the constraints are removed in the same way as the velocity components in RATTLE, Eqn. :eq:`rattle`, such that the forces vanish at a minimum of the energy under the constraints. These are the forces used by the minimization algorithms and for the convergence criterion *-toler*. After each step the coordinates are constrained with SHAKE, using the distance vectors before the step. With the Newton algorithm the step follows from the equations

.. math:: {\bf H}\,\Delta{\bf x} + {\bf J}^T\boldsymbol{\mu} ~=~ {\bf F},\qquad {\bf J}\,\Delta{\bf x} ~=~ 0,
   :label: newtoncons

where :math:`{\bf H}` is the second derivative matrix of the energy, :math:`{\bf F}` the forces, the rows of :math:`{\bf J}` are the derivatives of the constraints with respect to the coordinates and :math:`\boldsymbol{\mu}` are Lagrange multipliers. Since the steps close to the minimum are small, the constraints are satisfied to a relative tolerance of :math:`10^3` times the machine precision in minimization, rather than *-constr_toler*. Frozen atoms (option *-freeze*) are not moved by the constraints. A planar ring of which the angle constraints are not consistent with a planar ring, as discussed above, has to pucker during the minimization, which the L-BFGS and Newton algorithms may not achieve from an exactly planar starting structure. The steepest descent algorithm, or *-constrain bonds*, can be used instead.
