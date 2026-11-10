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

   where :math:`k_s` is the force constant of the core-shell spring, see Eqn. :eq:`vpol`. With :math:`\lambda = 1` this step would bring an isolated shell to its minimum, since the force of the spring then cancels the force :math:`{\bf F}_s` from the other particles. With hyperpolarization, the potential of the core-shell pair has an additional term :math:`k_h\,(r_{cs}-r_h)^4` for distances :math:`r_{cs}` larger than :math:`r_h`. The stiffness of the spring then depends on the direction. For a shell away from its core, the component of :math:`{\bf F}_s` perpendicular to the core-shell vector is divided by the ratio of the core-shell force and the distance

.. math:: k_s + 4k_h\,\frac{(r_{cs}-r_h)^3}{r_{cs}}\qquad\mbox{if}\quad r_{cs} > r_h,
   :label: kshyper

and the component along the core-shell vector by the second derivative of the potential

.. math:: k_s + 12k_h\,(r_{cs}-r_h)^2\qquad\mbox{if}\quad r_{cs} > r_h.
   :label: kshyperpar

This is a Newton step for an isolated shell. For a shell on top of its core, the quartic term only contributes if :math:`r_h = 0`, and then the displacement :math:`\Delta r` is the solution of :math:`\left|{\bf F}_s\right| = k_s\,\Delta r + 4k_h\,\Delta r^3`. The force constants are computed from the current polarizabilities and shell charges, such that they are correct during training as well.
#. Since the shells interact with each other, a full step (:math:`\lambda = 1`) can overshoot. If :math:`F_{ms}` increases, the step is reduced to :math:`\lambda = 1/2` for the rest of the relaxation. Smaller steps are not used, since they can freeze the shells far from their minimum.
#. A shell with one core is not allowed to move further away from it than a maximum distance :math:`d_{max}`. If the step brings a shell further away, it is moved back to the sphere with radius :math:`d_{max}` around its core, along the core-shell vector. This prevents the polarization catastrophe, where shells of neighbouring molecules attract each other without limit. With hyperpolarization the quartic term usually prevents it as well, then :math:`d_{max}` can be switched off.
#. The forces are computed again and the procedure is repeated from step 2, until the shells are converged or the maximum number of iterations is reached.

A shell held at :math:`d_{max}` by a force pointing away from its core can never reach a zero force. For such a shell only the component of the force perpendicular to the core-shell vector :math:`{\bf r}_{cs}` is used in Eqn. :eq:`fmsshell`

.. math:: {\bf F}_{s,\perp} ~=~ {\bf F}_s - \frac{{\bf F}_s\cdot{\bf r}_{cs}}{r_{cs}^2}\,{\bf r}_{cs}\qquad\mbox{if}\quad r_{cs} = d_{max}\quad\mbox{and}\quad{\bf F}_s\cdot{\bf r}_{cs} > 0,
   :label: fperp

which is zero at the minimum of the energy on the sphere. When a strong field holds the shell there, which is common with the poor parameters early in training, :math:`{\bf F}_{s,\perp}` is the difference of two large force components. It can then stay well above the tolerance on :math:`F_{ms}`, in a cycle or unchanged, even though the shell is practically at the minimum. Moving the shell over a small angle :math:`\theta` on the sphere changes its energy by about

.. math:: \Delta E ~\approx~ \frac{1}{2}\left|{\bf F}_s\right|\,d_{max}\,\theta^2,
   :label: ewall

since the force away from the core acts as a force constant :math:`\left|{\bf F}_s\right|/d_{max}` across the core-shell vector. Therefore, the force of such a shell counts as zero in Eqn. :eq:`fmsshell` if :math:`\left|{\bf F}_{s,\perp}\right| \le 10^{-4}\left|{\bf F}_s\right|`. For :math:`\left|{\bf F}_s\right| = 10^4` kJ/(mol nm) and :math:`d_{max}` = 0.02 nm the energy can then change by about :math:`10^{-6}` kJ/mol. The perpendicular force is computed as :math:`\left|{\bf F}_s\times{\bf r}_{cs}\right|/r_{cs}`, since the difference in Eqn. :eq:`fperp` loses its precision when the perpendicular force is that small compared to the force, in particular with forces in single precision. Without this, a few conformations with shells held at :math:`d_{max}` keep the relaxation of all the others going until the maximum number of iterations.

The precision of the forces and of the coordinates also limits how small :math:`F_{ms}` can get. The forces are sums of terms that can be much larger than the force itself. When :math:`F_{ms}` no longer decreases, the relaxation therefore stops if :math:`F_{ms}` is within :math:`(32\,\varepsilon\,F_{max})^2`, where :math:`\varepsilon` is the machine precision and :math:`F_{max}` the largest force component of any particle. Terms that cancel can be larger than :math:`F_{max}`, however. In the test force field *PG-pol-ions.xml*, for instance, chloride has a shell charge of -9 e on a core of about +8 e, which are practically on top of each other since its polarizability is at the lower bound of :math:`10^{-6}` Å\ :sup:`3`. A water shell at 0.3 nm then feels forces of about :math:`6\cdot 10^4` kJ/(mol nm) from both, which nearly cancel. In single precision their rounding errors add up to a few hundredths of a kJ/(mol nm), and in a dimer of water and chloride :math:`F_{ms}` went down to :math:`10^{-5}` (kJ/mol nm)\ :sup:`2` and then cycled between :math:`3\cdot 10^{-4}` and :math:`10^{-3}`. A tolerance of :math:`10^{-8}` was then not reached, while in double precision it was. The energies are as precise as the forces allow, but the conformation counts as not converged, takes the maximum number of iterations, and its shells are not kept with *-shell_keep*. With the default tolerance of :math:`10^{-3}`, and with :math:`10^{-6}`, all conformations converged. Therefore, in single precision, and hence on a GPU, do not use a tolerance far below the default.

The coordinates are integers in units of :math:`2^{-56}` nm, about :math:`1.4\cdot 10^{-17}` nm. A shell with a very stiff spring, such as that of the chloride above with :math:`k_s \approx 10^{13}` kJ/(mol nm\ :sup:`2`), does not move when its displacement in Eqn. :eq:`shellstep` is smaller than that unit, which happens at forces of about :math:`10^{-4}` kJ/(mol nm). A conformation of which no shell moved in a step counts as converged, since the next steps would be the same.

The maximum distance :math:`d_{max}` is part of the force field, since its parameters are trained with it. It is the option *shellmaxdistance* (in nm) of the polarization in the force field file,

.. code-block:: xml

   <interaction type="POLARIZATION" function="POLARIZATION" canswap="true">
     <option key="shellmaxdistance" value="0.02"/>

where 0 means that there is no maximum. It can be set with *alexandria gen_ff -shellmaxdistance*, or by editing the file, but not when training or using the force field. Without the option it is 0.02 nm. The other options that control the relaxation are listed in :numref:`tab-shellrelax`.

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
   | Start from the       | -shell_keep         |                    | yes (train_ff), from the cores in    |
   | previous shells      |                     |                    | simulate, see                        |
   |                      |                     |                    | :ref:`sec-shelltraining`             |
   +----------------------+---------------------+--------------------+--------------------------------------+

.. _sec-shelltraining:

----------------------------
Shell relaxation in training
----------------------------
During training the energies of all conformations (structures) of a molecule are computed together, see :ref:`sec-training`. The shells of all conformations are relaxed at the same time, and the conformations for which the shells converged are not computed again in the next iterations. For a dimer, the shells of each monomer are relaxed first in the isolated monomer, then all shells in the dimer. The induction energy is the change of the electrostatic energy upon relaxation of the shells. If the induction correction is trained separately (see :ref:`sec-indcorr`), the second order induction is computed by relaxing the shells of one monomer in the field of the other one, of which the shells stay at their positions in the isolated monomer, and the other way around. The rest of the induction is the induction correction.

A training step changes the parameters only a little, so the relaxed shells of the previous evaluation are much closer to the new minimum than their cores are. Therefore, by default (*train_ff -shell_keep*), the relaxation of each conformation starts from the shells of the previous evaluation by the same process, in the monomers and in the dimer. The shells are still put on their cores first to compute the energies before relaxation, so the induction energy is defined as before. The shells of conformations that did not converge are not kept, they start from their cores in the next evaluation. For the PC+GS force field of the electrostatics benchmark, starting from trained parameters, the number of iterations per conformation decreased from 2.6 to 0.4, which made training about 40% faster on the processor and about 17% faster on a GPU, where the energy calculations take a larger part of the time.

The relaxed positions depend on where the relaxation started, within the tolerance on :math:`F_{ms}`, which has no noticeable effect on :math:`\chi^2` with the default tolerance of 0.001 (kJ/mol nm)\ :sup:`2`. A smaller tolerance is not needed, and in single precision, for instance on a GPU, a tolerance of :math:`10^{-4}` cannot always be reached, such that more conformations fail to converge. With poor parameters, such as those of *-random_init* at the start of training, the shells can have more than one minimum, however. Shells of neighbouring molecules may then attract each other until one of them is held at :math:`d_{max}`, the polarization catastrophe, and which shell ends up at :math:`d_{max}` depends on where the relaxation started. A shell held at :math:`d_{max}` in a strong field can also have minima on different sides of its core. In both cases the energy of the conformation, and hence :math:`\chi^2`, depends on the previous evaluation. In a test with random parameters, 0.08% of the conformations ended in a different minimum than when starting from the cores, while with trained parameters none did. As a consequence, a training started from random parameters depends slightly on the order in which the processes evaluate the parameter sets, such that two runs with the same random number seed need not be identical, and the :math:`\chi^2` of the best parameters, which is recomputed at the end of the training, is checked to agree to within :math:`10^{-3}` relative to its value only, instead of :math:`10^{-3}` absolute. With *-noshell_keep* the shells start from their cores in every evaluation, which gives reproducible results.

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

With *-temp* larger than zero, the initial velocities are generated from a Maxwell distribution at that temperature, using the random number seed *-seed*. The motion of the center of mass is removed (see below), the velocity components along the constraints are removed, and finally the velocities are scaled such that the temperature, Eqn. :eq:`tempmd`, is exactly the requested one. There is no thermostat, the simulation conserves the total energy, which can be used to validate the forces. The energy file contains the potential energy terms, the kinetic energy :math:`E_{kin} = \frac{1}{2}\sum_i m_i{\bf v}_i^2`, the total energy and the temperature

.. math:: T ~=~ \frac{2E_{kin}}{N_{df}\,k_B},
   :label: tempmd

which are written every *-nstener* steps before the step is taken. The number of degrees of freedom is

.. math:: N_{df} ~=~ 3N - N_{com} - N_c,
   :label: ndf

for :math:`N` atoms with a mass, :math:`N_{com}` degrees of freedom of the motion of the center of mass that is removed and :math:`N_c` constraints, see :ref:`sec-constraints`. Shells and virtual sites do not contribute, since they have no mass.

The option *-comm* sets which motion of the center of mass is removed, from the initial velocities and every *-nstcomm* steps (default 100) after the velocities are updated, or only from the initial velocities if *-nstcomm* is 0. The choices are:

* *angular*: the translation of and the rotation around the center of mass (the default). This is needed for stable simulations of clusters in vacuum, since the rotation of the cluster otherwise takes up part of the kinetic energy. The velocity :math:`{\bf v}_{cm}` of the center of mass and the rigid rotation :math:`\boldsymbol\omega\times({\bf r}_i - {\bf r}_{cm})` are subtracted from the velocities, where the angular velocity :math:`\boldsymbol\omega = {\bf I}^{-1}{\bf L}` follows from the inertia tensor :math:`{\bf I}` and the angular momentum :math:`{\bf L}` around the center of mass :math:`{\bf r}_{cm}`. For a molecule that is linear in the starting coordinates, the rotation around the axis, with the smallest moment of inertia, is not removed. This removes :math:`N_{com} = 6` degrees of freedom, or 5 for a linear molecule.
* *linear*: only the translation of the center of mass, :math:`N_{com} = 3`.
* *none*: nothing, :math:`N_{com} = 0`.

Since this is a rigid body motion, the constraints are still satisfied afterwards. Without external forces, the linear and angular momentum are conserved, apart from numerical errors, such as those of the shell relaxation and the constraints, which are removed every *-nstcomm* steps. With the flat-bottom position restraints of *-kfbpr*, which act around the origin, the center of mass does move. Since the potential energy takes up part of the initial kinetic energy, the temperature during the simulation is lower than the initial one, about half of it if the starting structure is close to the minimum of the energy. The time step is set with *-deltat* (ps) and the number of steps with *-nsteps*.

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

The constraints are applied in energy minimization (*alexandria minimize* or *alexandria simulate -minimize*) as well, with all three algorithms that can be selected with *-minalg*. The starting coordinates are constrained first. After each computation of the forces, the components of the forces along the constraints are removed in the same way as the velocity components in RATTLE, Eqn. :eq:`rattle`, such that the forces vanish at a minimum of the energy under the constraints. These are the forces used by the minimization algorithms and for the convergence criterion *-toler*. After each step the coordinates are constrained with SHAKE, using the distance vectors before the step. With the Newton algorithm the step follows from the equations

.. math:: {\bf H}\,\Delta{\bf x} + {\bf J}^T\boldsymbol{\mu} ~=~ {\bf F},\qquad {\bf J}\,\Delta{\bf x} ~=~ 0,
   :label: newtoncons

where :math:`{\bf H}` is the second derivative matrix of the energy, :math:`{\bf F}` the forces, the rows of :math:`{\bf J}` are the derivatives of the constraints with respect to the coordinates and :math:`\boldsymbol{\mu}` are Lagrange multipliers. Since the steps close to the minimum are small, the constraints are satisfied to a relative tolerance of :math:`10^3` times the machine precision in minimization, rather than *-constr_toler*. Frozen atoms (option *-freeze*) are not moved by the constraints. A planar ring of which the angle constraints are not consistent with a planar ring, as discussed above, has to pucker during the minimization, which the L-BFGS and Newton algorithms may not achieve from an exactly planar starting structure. The steepest descent algorithm, or *-constrain bonds*, can be used instead.
