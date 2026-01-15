****************************
Force Field Training Targets
****************************
A number of different targets can be used for training force fields in the ACT.
Below, we mention the most important ones including,  where needed, the technical details necessary to appreciate the methods.
For information on obtaining or generating training data we refer to Section :ref:`sec-data`.

Computing interaction energies and components of interaction energies is a crucial part of force field development. In the ACT we have hitherto used data from symmetry-adapted perturbation theory (SAPT) calculations :cite:p:`Jeziorski1994a,Parker2014a`. 
It is not entirely trivial to match the energy components from SAPT to force field terms, however.

=============================================
Algorithm to compute energy components in ACT
=============================================
The component of the interaction energy of a dimer can be computed from the difference between dimer and monomers :math:`A` and :math:`B` :cite:p:`chalasinski2000state`:

.. math:: E_x^{inter}(AB) = E_x^{total}(AB) - \left(E_x(A) + E_x(B)\right)
   :label: vinter

where :math:`x` is exchange or dispersion and :math:`total` indicates that the energy includes both the intra- and intermolecular interactions.
To compute the electrostatics or induction energy of a dimer in the gas phase, the ACT first computes the relaxed energy of the two monomers :math:`A` and :math:`B`, that is, the energy of the shell particles is minimized with respect to their positions, yielding :math:`E_x(A)` and :math:`E_x(B)`. The rationale for this is that SAPT computes electrostatic energies between unperturbed monomers based on the response/relaxation of monomer Hartree-Fock (HF) orbitals in the electric field of the interacting partner.
Then, the energies of the dimer :math:`AB` are computed in three steps: 

#. electrostatics is computed with shells located in the relaxed monomer positions, yielding :math:`E_{elec}^{total}(AB)`,
#. the shells of compound :math:`A` are allowed to relax (further) in the electric field from compound :math:`B`, while shells of compound :math:`B` remain at their monomer positions, and vice versa, yielding the second order relaxation :math:`E_{induc}^{inter(2)}` (see ref.  :cite:p:`McDaniel2013a` for details),
#. the shells are allowed to relax completely, yielding the total :math:`E_{induc}` from which the higher order terms, named :math:`E_{induc}^{inter(3)}` here for convenience, can be derived by subtracting  :math:`E_{induc}^{inter(2)}`. According to ref. :cite:p:`McDaniel2013a`, parameters of models corresponding to the higher other terms, including, potentially, charge transfer, can be trained to the :math:`\delta_{HF}` contribution of the SAPT induction energy. Here, we have added the exponential term proposed by McDaniel and Schmidt for this purpose (section :ref:`sec-indcorr`).
 

The terms below can be compared directly to SAPT:

.. math:: E_{elec} ^{inter}(AB) = E_{elec}^{total}(AB) - \left(E_{ei}(A) + E_{ei}(B)\right)
   :label: velec

.. math:: E_{induc}^{inter(2)}(AB) = \left(E^{total}_{ei}\|_A + E^{total}_{ei}\|_B\right)-2E_{elec} ^{total}(AB)
   :label: vind2
   
.. math:: E_{induc}^{inter(3)}(AB) = E_{ei}^{total}(AB) - E_{induc}^{(2)}(AB) - E_{elec}^{total}(AB)
   :label: vind3

where :math:`ei` is short for :math:`elec+induc` and the notation :math:`\|_X` indicates that the shells of compound :math:`X` are kept fixed in the relaxed monomer conformation. 
If we sum Eqns. :eq:`velec`- :eq:`vind3` we recover Eqn. :eq:`vinter` where :math:`x` equals :math:`ei`.
Eqn. :eq:`velec` corresponds to the electrostatics in SAPT.

To summarize, ACT computes the induction term in two parts: a second order term :math:`E_{induc}^{inter(2)}(AB)` (Eqn. :eq:`vind2`) and a third-and-higher-order term :math:`E_{induc}^{inter(3)}(AB)` (Eqn. :eq:`vind3`), the sum of which corresponds to the total induction from SAPT. 

McDaniel and Schmidt :cite:p:`McDaniel2013a` proposed that  Eqn. :eq:`vind3` should be equal to the :math:`\delta_{HF}` + :math:`\delta_{MP2}` terms from SAPT while Eqn. :eq:`vind2` would correspond to the polarization energy.
They then continue to suggest how this can be implemented in a force field:

.. math:: E_{induc} ~=~ E_{shell} + E_{pol} + E_{\delta HF} 

where :math:`E_{shell}` is Eqn. :eq:`vpol` and

.. math:: E_{pol} ~=~ A^{ind}_{ij} {\rm exp}(-b_{ij} r_{ij})

where :math:`r_{ij}` is the interatomic distance and :math:`b_{ij}` a constant, and

.. math:: E_{\delta HF} ~=~ A^{\delta HF}_{ij} {\rm exp}(-b_{ij} r_{ij})

where both :math:`A^{ind}_{ij}` and :math:`A^{\delta_{HF}}_{ij}` are determined from a negative geometric combination rule

.. math:: A_{ij} ~=~ -\sqrt{A_i A_j}

meaning these terms are always attractive. These potentials use the :math:`b_{ij}` that is used in the exchange energy (using a Buckingham potential, Eqn. :eq:`vbh`).
Whether this is the most appropriate way of splitting terms and reproducing SAPT energies remains to be determined.

In the output the ACT training module \actcmd{train\_ff} there are two terms related to induction. The term INDUCTIONCORRECTION refers to Eqn. :eq:`vind3`. If that is present, the term INDUCTION refers to Eqn. :eq:`vind2`, if not it refers to the sum of the two.

===========================
Monomer Energies and Forces
===========================
Typically, a series of single-point quantum calculations are done at a user-defined level of theory. These calculations can then be converted to a molprop file. More information to come.

================
Other Properties
================
In principle, all the molecular properties mentioned in Section :ref:`sec-properties` can be used for training, but it is highly recommended to leave some properties for validation. 
