.. _sec-properties:
        
*******************************
Predicting Molecular Properties
*******************************
The ACT can be used to perform MD simulations of clusters in the gas-phase using the *alexandria simulate* command. This module includes the possibility to perform energy minimizations with the flag *-minimize*.
For simulations employing periodic boundaries the OpenMM package~ :cite:p:`Eastman2023a` should be used instead. 
For more details about MD simulations, see Section :ref:`sec-simulations`.

Below we describe some of the properties that can be computed using the ACT.

.. _sec-esp:

=======================
Electrostatic Potential
=======================
The charge distribution :math:`\rho` of a molecule is determined by the nuclear position of the :math:`N` atoms at positions :math:`\mathbf{x}` and the electron density :math:`n(\mathbf{r})`. The molecular electrostatic potential (MEP) at a point in space :math:`r^\prime` is thus given by

.. math:: \Phi(\mathbf{r}^\prime) = \frac{1}{4\pi\varepsilon_0}\left[\sum_{i=1}^N \frac{z_i}{ \|\mathbf{x}_i-\mathbf{r}^\prime\|} - \int \frac{n(\mathbf{r})}{\|\mathbf{r}-\mathbf{r}^\prime\|}{\rm d}r \right]
   :label: phi

where :math:`N` is the number of atoms, :math:`z_i` are the nuclear charges, :math:`\varepsilon_0` is the permittivity of vacuum and integration is over the entire space. The minus sign before the integral is due to the negative charge of electrons. 

Accurate knowledge of the MEP contributes to, for example, the understanding of interactions and function of biological macromolecules in solution~ :cite:p:`Larsson2012b`. For a molecule in the gas phase, Eqn. :eq:`phi` can be evaluated using density functional theory and wave function quantum chemistry, albeit at a significant computational cost. 
Databases of such calculations for small molecules are available to facilitate reuse :cite:p:`Kriz2023a`. For  large condensed-phase systems, however, it is common to apply classical force fields, where electrons are not taken into account explicitly. Instead, effective partial charges on atoms are used. The electronic degrees of freedom, charge polarization, is sometimes taken into account through induced point dipoles and higher electrostatic moments, or by using a core-shell model~ :cite:p:`Dick1958a,Jordan1995a,Maaren2001a`. For additional background we refer to some excellent reviews~ :cite:p:`Dauber-Osguthorpe2019a,Hagler2019a,Jing2019a`.

The MEP can be used as a target in model development in the ACT, however we recommend against that for both fundamental and practical reasons~ :cite:p:`Hosseini2026a`.

=====================
Electrostatic Moments
=====================
If point :math:`{\bf r}^\prime` in Eqn. :eq:`phi` is outside the distribution of electron
density and :math:`{\bf r}^\prime >> {\bf r}`, the electrostatic potential can be evaluated through the  
Taylor expansion of :math:`|{\bf r}-{\bf r}^\prime|^{-1}`  :cite:p:`Berendsen2007a`:

.. math:: \frac{1}{|{\bf r^\prime}-{\bf r}|} \approx \frac{1}{r} + \left({\bf \hat r} \cdot
   {\bf r} \right)\frac{1}{r^2} + \frac{1}{2} \left[3\left({\bf \hat r} \cdot
   {\bf r} \right)^2 - r^{2}{\rm \bf I}\right]\frac{1}{ r^3} + \cdots
   :label: Rexpansion

where :math:`{\bf \hat r} = {\bf r}^\prime/r` and :math:`{\rm \bf I}` is the
identity matrix. By inserting Eqn. :eq:`Rexpansion` into Eqn. :eq:`phi`, we get

.. math:: \Phi({\bf r}^\prime) \approx \frac{1}{4 \pi \epsilon_0} { \frac{Q}{r}
   + \frac{{\bf \mu_0}}{r^2} + \frac{\boldsymbol{\Theta_0} }{r^3} + \cdots}

:math:`Q` is the monopole moment, sometimes called the zeroth
moment of the molecular electron density.  In principle this is the total charge of the molecule. In what follows we combine both the atomic cores and the electron density into :math:`n(r)`. Then, :math:`Q` is given by

.. math:: Q = \int n(\mathbf{r}) d {\bf r}  

:math:`{\bf \mu_0}` is the vector of the permanent dipole moment, which measures 
the polarity of the molecular electron density. :math:`{\bf \mu_0}` is given by

.. math:: {\bf \mu_0} = \int n(\mathbf{r}) \left({\bf \hat r} \cdot {\bf r}\right) d {\bf r} 

:math:`\boldsymbol{\Theta_0}` is the tensor of the permanent quadrupole moment,  which exhibits 
the deviation of the distribution of the molecular electron density
from spherical symmetry. :math:`\boldsymbol{\Theta_0}` is given by

.. math:: \boldsymbol{\Theta_0} = \frac{1}{2}\int n(\mathbf{r}) \left[3\left({\bf \hat r}
   \cdot {\bf r}\right)^2 - r^2 {\rm \bf I} \right] d {\bf r} 

The quadrupole tensor can be written as a traceless :math:`3 \times 3` matrix 
in the Cartesian coordinate if one writes :math:`\left({\bf \hat r} \cdot
{\bf r} \right)^2` as:

.. math:: \left({\bf \hat r} \cdot {\bf r} \right)^2 = {\bf \hat r} \cdot \left({\bf r} {\bf r} \right) \cdot {\bf \hat r}

where :math:`{\bf rr}` is the outer product of vector :math:`{\bf r}` with itself. This
results in a matrix that can be written in terms of the Cartesian 
components of the vector :math:`{\bf r}` as follows:

.. math:: {\bf rr} = 
   \left[\begin{array}{cccc}
   x^2 & xy   & zx \\
   yx  &  y^2 & yz \\
   zx  &  zy   & z^2 \\
   \end{array}\right]

Finally, we get

.. math:: \begin{align}
   \frac{1}{2} \left[3\left({\bf \hat r} \cdot
   {\bf r} \right)^2 - r^{2} {\rm \bf I}\right] = {\bf \hat r} \cdot
   \frac{1}{2}\left[\begin{array}{cccc}
   3x^2- r^2 & 3xy   & 3zx \\
   3yx  &  3y^2- r^2 & 3yz \\
   3zx  &  3zy   & 3z^2- r^2 \\
   \end{array}\right]
   \cdot
   {\bf \hat r} 
   \end{align}
   
==============
Polarizability
==============
The shape of the molecular electron density changes when it interacts with an external electric
field; hence, the total energy of the molecule
changes. The static response of a molecule to a homogeneous external
electric field, (:math:`\bf F`), can be studied by expanding its energy in a
Taylor series  :cite:p:`Jensen2007a,Calaminici1998a`: 

.. math:: E \left( {\bf F} \right) = E\left({\bf 0}\right) + \frac{\partial E}{\partial {\bf F}}\biggm\lvert_{{\bf F}=0}
   {\bf F} + \frac{1}{2}
   \frac{\partial^2E}{\partial {\bf F}^2 }\biggm\lvert_{{\bf F}=0}
   {\bf F}^2 + \frac{1}{6}
   \frac{\partial^3E}{\partial {\bf F}^3 }\biggm\lvert_{{\bf F}=0}
   {\bf F}^3 + \cdots

where

.. math:: -\frac{\partial E}{\partial {\bf F}}\biggm\lvert_{{\bf F}=0}  = {\bf \mu_0}

.. math:: - \frac{\partial^2E}{\partial {\bf F}^2}\biggm\lvert_{{\bf F}=0} = {\bf \alpha}

.. math:: - \frac{\partial^3E}{\partial {\bf F}^3}\biggm\lvert_{{\bf F}=0} = {\bf \beta}

where :math:`{\bf \mu_0}` is the vector of permanent dipole moment, :math:`{\bf \alpha}`
is the tensor of polarizability, which is the 
linear part of the response of the molecular electron density 
with respect to the external electric field, and :math:`\boldsymbol{\beta}` is the
first hyperpolarizability  :cite:p:`Jensen2007a`. 

Instead of expanding the energy, we can expand the dipole moment of a molecule in an
external electric field  :cite:p:`Calaminici1998a`, written as

.. math:: \mu = \mu_0 + {\bf \alpha} {\bf F} +
   \frac{1}{2} {\bf \beta} {\bf F}^2 + \cdots

where :math:`{\bf \alpha} {\bf F}` gives the vector of induced dipole
moment, :math:`{\bf \mu_1}`  :cite:p:`Jensen2007a`: 

.. math:: {\bf \mu_1} = {\bf \alpha}{\bf F}

that can be written in matrix form as

.. math:: \left[\begin{array}{c}
   \mu_x \\
   \mu_y \\
   \mu_z \\
   \end{array}
   \right] =
   \left[\begin{array}{cccc}
   \alpha_{xx} & \alpha_{xy} & \alpha_{xz} \\
   \alpha_{yx} & \alpha_{yy} & \alpha_{yz} \\
   \alpha_{zx} &\alpha_{zy}  & \alpha_{zz} \\
   \end{array}
   \right]
   \left[\begin{array}{c}
   F_x\\
   F_y\\
   F_z\\
   \end{array}
   \right]
   :label: matrix2

From the polarizability tensor the polarizability isotropy
:cite:p:`Stone2013a,Calaminici1998a`,

.. math:: \bar \alpha = \frac{\left(\alpha_{xx} + \alpha_{yy} + \alpha_{zz} \right)}{3}

and the polarizability anisotropy  :cite:p:`Stone2013a`, 

.. math:: \Delta \alpha = \sqrt{[(\alpha_{xx} - \alpha_{yy})^2 + (\alpha_{xx} -
   \alpha_{zz})^2 + (\alpha_{yy} - \alpha_{zz})^2 + 6 (\alpha_{xy}^2 +
   \alpha_{xz}^2 + \alpha_{yz}^2)]/2}

can be calculated.

.. _sec-nma:

====================
Normal Mode Analysis
====================
Please note that the text below is largely taken (with permission) from a paper by Henschel {\em et al.}~ :cite:p:`Henschel2020a`.

The ACT contains the *alexandria nma* tool that performs a normal mode analysis to determine the vibrational frequencies of a compound.
Vibrational frequencies are required to compute the IR spectra and thermochemistry of molecules. The normal modes of molecular vibrations can be obtained by eigenvalue decomposition of the Hessian matrix, whose elements are the second derivatives of the energy with respect to the atomic coordinates :math:`q`.

.. math:: H_{ij} = \frac{\partial^2E}{\partial q_i \partial q_j}

where :math:`i` and :math:`j` run from :math:`0` to :math:`N-1`, where :math:`N` is the number of atoms in the molecule. If virtual sites :math:`v` are used, for example, to model the :math:`\sigma`-hole for halogen atoms, the energy,
E, depends on the positions of both atoms and virtual sites;
that is, :math:`E= E(q_0,...,q_{N-1},v_0,...,v_{M-1})`, where the positions of the :math:`M`
virtual sites, :math:`v`, in the compound are a function of the atomic
coordinates q. 

The Hessian is computed numerically---the
N atoms are moved independently in all three spatial dimensions,
and the forces are computed. From these forces, the second
derivative of the energy is then evaluated numerically. Note that the positions of the virtual sites are
updated before each force calculation, which means that their influence on the Hessian is taken into account explicitly when computing :math:`H`.

.. _sec-irspectrum:

----------------
Infrared Spectra
----------------
For the calculation of a full IR spectrum, in addition to the vibrational frequencies, the intensities and the line shapes are required.

In case of the quantum chemical calculations both the eigenfrequencies and the corresponding IR intensities are produced by default when a frequency calculation is requested in, for instance, the Gaussian software~ :cite:p:`g16`. The frequencies for about 5000 compounds  are available from the Alexandria library~ :cite:p:`Ghahremanpour2018a` at the B3LYP/aug-cc-pvtz level of theory~ :cite:p:`Becke1988a,Kendall1992a,Woon1993a,Woon1993b`.Details of the quantum chemical calculations from which the frequencies were obtained have been presented previously :cite:p:`Ghahremanpour2016a,Ghahremanpour2018a`.

For the force field calculations, the intensities :math:`I_{n}` were derived from the transition dipole derivatives:

.. math:: I_{n}=\sum_{k=1}^{3}\left(\frac{\partial p_{k}}{\partial Q_{n}}\right)^2
   :label: In

where :math:`k` iterates over cartesian dimensions, :math:`p` is the dipole moment of the molecule, and :math:`Q_n` the normal coordinate :math:`n`. In order to take into account virtual sites :math:`v` we note that 

.. math:: p_{k} = p_{k}(q_0,...,q_{N-1},v_0,...,v_{M-1})

and rewrite equation :eq:`In` as:

.. math:: I_{n}=\sum_{k=1}^{3}\left(\frac{\partial p_{k}}{\partial q_s}\frac{\partial q_{s}}{\partial Q_{n}}\right)^2

where :math:`s` iterates over the :math:`N` atomic coordinates. 
To make the calculation of intensities practical, the numerical derivative of the dipole moment with respect to the atomic coordinates :math:`\displaystyle{\frac{\partial p_{k}}{\partial{q_s}}}` is stored in a text file during the normal mode analysis and finally we note that the term :math:`\displaystyle{\frac{\partial q_{s}}{\partial Q_{n}}}` corresponds to one over component :math:`s` of eigenvector :math:`n`.

As an example, an infrared spectrum can be generated using a force field file and a molprop file using::

    alexandria nma -ff OPLS2020 -charges OPLS2020-charges -db ethanol -ir ir-ethanol

yielding the spectrum in Fig. :numref:`fig-ethanol`.

.. figure:: ../images/ethanol-irspectrum-opls2020.pdf
   :name: fig-ethanol
   :width: 90%
   
   Simulated infrared vibrational spectrum for ethanol, based on the OPLS2020 force field~ :cite:p:`Jorgensen2023a`.

---------------
Thermochemistry
---------------
The `Canonical`_ ensemble :math:`Q(N,V,T)` can be used to compute the molecular `Internal`_ energy :math:`E` (:eq:`intE`),
the standard `Entropy`_ :math:`S^o` (:eq:`S0`), and the `Heat-Capacity`_ :math:`C_v` (:eq:`Cv`) at constant volume, using

.. _Canonical: https://en.wikipedia.org/wiki/Canonical_ensemble

.. _Internal: https://en.wikipedia.org/wiki/Internal_energy

.. _Entropy: https://en.wikipedia.org/wiki/Standard_molar_entropy

.. _Heat-Capacity: https://en.wikipedia.org/wiki/Heat_capacity

.. math:: E = RT^2\left(\frac{\partial {\rm ln} Q}{\partial T} \right)_{N,V}
   :label: intE

.. math:: S^o = R {\rm ln} Q + RT\left(\frac{\partial {\rm ln} Q}{\partial T} \right)_{N,V}
   :label: S0

.. math:: C_v = 2RT\left(\frac{\partial {\rm ln} Q}{\partial T} \right)_{N,V} + RT^2 \left(\frac{\partial^2 {\rm ln} Q}{\partial T^2} \right)_{N,V}
   :label: Cv

where :math:`R` is the ideal gas constant and :math:`T` the absolute temperature. For an ideal 
gas, :math:`Q(N,V,T)` can be decomposed into `Partition-Function`_ of different degrees of freedom: electronic (el), translational (tr), rotational (rot) and vibrational (vib) motions. Therefore, for a molecular ideal gas, :math:`Q(N,V,T)`  can be expressed as:

.. _Partition-Function: https://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)

.. math:: Q(N,V,T) = \frac{\left(q_{\rm el} q_{\rm tr} q_{\rm rot} q_{\rm vib}\right)^N}{N!}
   :label: Partition


The `Rigid-Rotor`_ and the `Quantum-Harmonic-Oscillator`_ can be used to approximate the contribution of the rotational and vibrational motions. The partition function of a rigid rotator is defined as 

.. _Rigid-Rotor: https://en.wikipedia.org/wiki/Rigid_rotor

.. _Quantum-Harmonic-Oscillator: https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator

.. math:: q_{\rm rot} = \frac{T}{\sigma \Theta_{\rm rot}}

where :math:`\sigma` is the symmetry number, :math:`\Theta_{\rm rot}` the rotational temperature defined as 

.. math:: \Theta_{\rm rot} = \frac{h^2}{8\pi^2 {\rm I} k_{\beta}}

where :math:`{\rm I}` is the moment of inertia, :math:`k_{\beta}` the Boltzmann constant, and :math:`h` the Planck constant. The partition function of a quantum harmonic oscillator is defined as 

.. math:: q_{\rm vib} = \frac{e^{\frac{-\beta h \nu}{2}}}{1-e^{-\beta h \nu}}

where :math:`\nu` is the vibrational frequency of the oscillator. If we define the vibrational temperature as :math:`\Theta_{\rm vib} = \frac{h \nu}{k_{\beta}}`, then we will have 

.. math:: q_{\rm vib} = \frac{e^{\frac{\Theta_{\rm vib}}{2T}}}{1- e^{\frac{\Theta_{\rm vib}}{T}}}

Applying Eqn. :eq:`Partition` to Eqn. :eq:`intE` followed by the multiplication rule in logarithm yields:

.. math:: E = E_{tr} + E_{rot} + E_{vib}

and similarly, 

.. math:: C_v = C_{tr} + C_{rot} + C_{vib}

.. math:: S^o = S_{tr} + S_{rot} + S_{vib}

Thermochemical properties are computed automatically by the *alexandria nma* command. For more details, please see Van der Spoel *et al.*  :cite:p:`Spoel2018a`.

The *alexandria nma* used above to generate the infrared spectrum in Fig. :numref:`fig-ethanol`  computes the thermochemical variables as well (Table :numref:`tab-thermo`).
Information on calculation of the enthalpy of formation will be pubslished in the near future.

.. table:: Thermochemical values for ethanol at 298.15 K based on the OPLS2020 force field~ :cite:p:`Jorgensen2023a`.
   :name: tab-thermo
   
   +-----------------------+---------------------------------+-----------+
   | Property              | Experiment                      | OPLS2020  |
   +=======================+=================================+===========+
   | :math:`S^o` (J/mol K) | 281 :cite:p:`Handbook2022a`     | 273.0     |
   +-----------------------+---------------------------------+-----------+
   | :math:`C_v` (J/mol K) | 57 :cite:p:`Ghahremanpour2016a` | 62.2      |
   +-----------------------+---------------------------------+-----------+

=========================
Second Virial Coefficient
=========================
The second `Virial-Coefficient`_ is the second term in the `Virial-Expansion`_ that describes the deviation from the ideal gas law for real gases:

.. math:: \displaystyle{\frac{P}{RT\rho}} = A + B_2(T)\rho + C_3(T)\rho^2 + ...

.. _Virial-Coefficient: https://en.wikipedia.org/wiki/Virial_coefficient

.. _Virial-Expansion: https://en.wikipedia.org/wiki/Virial_expansion

with :math:`P` the pressure, :math:`R` the gas constant, :math:`T` the temperature and :math:`\rho` the density.

:math:`B_2(T)` is a useful property gauging interactions in the gas phase because experimental values are available for close to 2000 compounds as a function of temperature. It is computed from an integral weighting the interaction between two molecules over three-dimensional space.

.. math:: B_2^{cl}(T) = -\frac{1}{2}\int_0^{\infty} \left< e^{-\beta u_{12}({\mathbf r})} - 1\right> d{\mathbf r}

where :math:`u_{12}({\mathbf r})` is the interaction energy between two compounds (Eqn. :eq:`vinter`), :math:`\beta = 1/k_BT` and
the integral is over all space and relative orientations of the compounds. If we sample these adequately (including at close, repulsive, distance) we can simplify the integral to a one-dimensional one:

.. math:: B_2^{cl}(T) = -2\pi\int_0^{\infty} r^2 \left< e^{-\beta u_{12}(r)} -1 \right> dr  

The above equation is entirely classical and quantum corrections have to be added according to:

.. math:: B_2^F(T) = \frac{\hbar^2}{24(k_BT)^3} \sum_{j=1}^{2} \left[\frac{\left< {\mathbf F}_j^2 \right>}{m_j}\right]

for the force on the compounds and

.. math:: B_2^{\tau}(T) = \frac{\hbar^2}{24(k_BT)^3} \sum_{j=1}^2\left[\sum_{\alpha=x,y,z} \frac{\left< \tau^2_{j,\alpha}\right>}{I_{j,\alpha}} \right]

for the torque on the compounds,
where :math:`m_j` is the mass of the molecules :math:`j` and :math:`\mathbf{F}^2` is the averaged square force on one molecule given by

.. math:: \left<\mathbf{F}^2\right> = k_B T \int_0^{\infty} \left< e^{\displaystyle{- \beta u_{12}(\mathbf{r})}}\left[\nabla u_{12}(\mathbf{r})\right]^2\right> d \mathbf{r} 

and where :math:`I` is the moment of inertia of the molecule and
:math:`\tau^2` is the average square torque on one molecule defined by

.. math:: \left<\tau^2_{j,\alpha}\right> = k_B T \int_0^{\infty}  \left< e^{{- \beta u_{12}(\mathbf{r})}} \left[\nabla_{\omega} u_{12}(\mathbf{r})\right]_{j,\alpha}^2\right> d \mathbf{r} .

The change in :math:`B_2(T)` as a function of temperature can be used to scrutinize the repulsive and attractive components of the potential energy. :math:`B_2(T)` is negative at low temperatures due to attraction forces, while it becomes positive at higher temperatures as repulsion forces start to dominate, and passes through a maximum and eventually decreases at very high temperatures where repulsion force are fully dominant~ :cite:p:`Amdur1958a`.  

Code to compute the second virial coefficient is available in the *alexandria b2* command. Since the calculation is relatively expensive it has been implemented to use parallel processing using the message passing library. You can run it on a 16-core machine like::

  mpirun -n 16 alexandria b2 -v 3 -g TIP4P -b2 TIP4P -ff TIP4P
  -f water#water.pdb -T1 373.15 -T2 673.15 -dT 25.0
  -maxdimer 32768

where TIP4P corresponds to the well-known water model~ :cite:p:`Jorgensen1983a`, the :math:`T_1` and :math:`T_2` are the temperature limits (note gas-phase for water), :math:`dT` is the temperature interval and maxdimer determines how many relative orientations will be evaluated.
Due to the underlying algorithm for `Quasi-Random`_ numbers, this should be a power of two. The result is plotted in Fig. :numref:`fig-b2`.

.. _Quasi-Random: https://en.wikipedia.org/wiki/Sobol_sequence

.. figure:: ../images/TIP4P-b2.pdf
   :align: center
   :width: 90%
   :name: fig-b2

   Sample second virial coefficient and components for water using the TIP4P model~ :cite:p:`Jorgensen1983a`.

