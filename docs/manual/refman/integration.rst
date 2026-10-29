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

   where :math:`k_s` is the force constant of the core-shell spring, see Eqn. :eq:`vpol`. With :math:`\lambda = 1` this step would bring an isolated shell to its minimum, since the force of the spring then cancels the force :math:`{\bf F}_s` from the other particles. With hyperpolarization, the force constant of the quartic term at the current core-shell distance is added to :math:`k_s`. The force constants are computed from the current polarizabilities and shell charges, such that they are correct during training as well.
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

Only particles with a mass are integrated. Virtual sites are constructed from the atoms before the forces are computed, and the forces on them are distributed over the atoms they are constructed from. The shells are relaxed for each new configuration of the atoms as described in :ref:`sec-shellrelax`, starting from their cores. With the option *-constrain* all bond lengths are kept constant at the reference bond lengths of the force field, using SHAKE after Eqn. :eq:`vvpos`, and the velocity components along the bonds are removed using RATTLE after Eqn. :eq:`vvvel` :cite:p:`Ryckaert1977a,Andersen1983a`. The starting coordinates are constrained as well. The relative tolerance and the maximum number of iterations of SHAKE and RATTLE are set with *-constr_toler* and *-constr_maxiter*.

With *-temp* larger than zero, the initial velocities are generated from a Maxwell distribution at that temperature, using the random number seed *-seed*. The motion of the center of mass is removed, the velocity components along constrained bonds are removed, and finally the velocities are scaled such that the temperature, Eqn. :eq:`tempmd`, is exactly the requested one. There is no thermostat, the simulation conserves the total energy, which can be used to validate the forces. The energy file contains the potential energy terms, the kinetic energy :math:`E_{kin} = \frac{1}{2}\sum_i m_i{\bf v}_i^2`, the total energy and the temperature

.. math:: T ~=~ \frac{2E_{kin}}{N_{df}\,k_B},
   :label: tempmd

which are written every *-nstener* steps before the step is taken. The number of degrees of freedom is

.. math:: N_{df} ~=~ 3N - 3 - N_c,
   :label: ndf

for :math:`N` atoms with a mass and :math:`N_c` constraints (the number of bonds when *-constrain* is used, zero otherwise). Shells and virtual sites do not contribute, since they have no mass. Three degrees of freedom are subtracted for the center of mass, which does not move, since there are no external forces and the initial velocities are either zero or have the center of mass motion removed. Since the potential energy takes up part of the initial kinetic energy, the temperature during the simulation is lower than the initial one, about half of it if the starting structure is close to the minimum of the energy. The time step is set with *-deltat* (ps) and the number of steps with *-nsteps*.
