******************
Force field design
******************

===================
Physical Background
===================
The Alexandria force-fields have a functional form consisting of van der Waals (*vdw*), electrostatics (*coul*), polarization (*pol*)  and bonded terms including a radial (*b*), angular (*a*), out-of-plane dihedral (*i*) and torsion (*d*) terms.

.. math:: E ~=~ V_{vdw}(r_{ij})+V_{coul}(r_{ij})+V_{pol}(r_{cs})+V_b(r_{ij})+V_a(\theta_{ijk})+V_i(\phi_{ijkl}) + V_{d}(\phi_{ijkl})
   :label: ACT_terms

where :math:`r_{ij}` refers to the distance between atoms :math:`i` and :math:`j`, :math:`r_{cs}` to the distance between a core and a shell (Drude) particle, :math:`\theta_{ijk}` to the angle given by three atoms :math:`i`, :math:`j` and :math:`k`, and :math:`\phi_{ijkl}` to a torsional angle between the planes given by atoms :math:`i,j,k` and atoms :math:`j,k,l`.
For each of the terms, multiple functional forms are available, so that within the Alexandria framework, different force-fields, including previously published ones, can be reconstructed and compared to one another in a systematic manner :cite:p:`Spoel2021a`. For details on these potentials we refer to Chapter :ref:`sec-energy`.


In total there are seven *atom* parameter types and 7-14 *bond* parameter types. 
A variety of virtual sites can be used, like those used in water models :cite:p:`Maaren2001a,Lamoureux2003b` or to model anisotropy due to :math:`\sigma`-holes on halogen atoms or water :cite:p:`Kriz2024b`.

.. _sec-partialq:

===========================
Determining partial charges
===========================
The electrostatic potential (ESP, Section :ref:`sec-esp`) has historically been used 
to determine partial charges :cite:p:`Besler1990a` and the ACT supports training models to reproduce the ESP.
However, in a very recent paper we have shown that fitting charges to reproduce the ESP in a limited volume around a compound is fundamentally flawed due to lack of information :cite:p:`Hosseini2026a`. 
If the purpose is to build models that reproduce electrostatic interactions, this can be done directly by training models to reproduce SAPT energy components. Here, the split-charge equilibration (SQE) algorithm :cite:p:`Verstraelen2009a` is used to generate the effective partial charge on each atom in a molecule. SQE, in turn, is based on the electronegativity equalization method (EEM), as developed by Rapp{\'e} and Goddard :cite:p:`Rappe1991a`. In brief, EEM uses the atomic hardness :math:`\eta` and electronegativity :math:`\chi` to determine the atomic charges in a molecule from a Taylor expansion of the molecular energy in terms of charges. 
The SQE algorithm introduces a correction to the atomic electronegativities for bonded atoms :math:`\Delta\chi` as well as a bond hardness :math:`\Delta\eta`. With this addition, charge can ``flow" through bonds only, which overcomes issues with over-polarization in the EEM :cite:p:`Nistor2006a`. 

The ACT code implements the possibility to generate charges for compounds in dimers or clusters where  charge transfer between compounds is disallowed which is a reasonable approximation since charge transfer has been shown to have limited impact on the binding energy of non-covalent complexes :cite:p:`CYin2019a`. For the SQE algorithm two atomic parameters  (:math:`\chi` and :math:`\eta`) as well as two bond parameter types (:math:`\Delta\chi` and :math:`\Delta\eta`) need to be determined and the ACT can train SQE parameters to reproduce electrostatic and induction energies :cite:p:`Hosseini2026a`. 
For background information we refer the reader to an excellent review by Jensen :cite:p:`Jensen2023a`, but below follows a break-down of using SQE with shells or virtual sites.

=================================================
Charge Equilibration with Shells or Virtual Sites
=================================================
Among the approaches to modeling the charge-dependent component of a force field, those rooted in the chemical potential equalization principle are especially notable, as the principle stems directly from density functional theory :cite:p:`itskowitz1997chemical`. The first computational implementation of the chemical potential equalization principle was the electronegativity equalization method (EEM) :cite:p:`Mortier1986a,Rappe1991a`. However, due to limitations of this model, Chelli {\em et al.} proposed the atom-atom charge transfer (AACT) model :cite:p:`chelli1999electrical`. Later, Nistor and co-workers combined the EEM and AACT approaches into a single framework, the split-charge equilibration (SQE) model, which fulfills the essential criteria for a successful charge-transfer potential :cite:p:`Nistor2006a,Verstraelen2009a`.
The ACT implements both the EEM and the SQE as algorithms for determining partial charges. 

In brief, EEM minimizes an empirical model of the intramolecular electrostatic energy (computed from the  atomic electronegativity :math:`\chi_i` and atomic hardness :math:`\eta_i`) with respect to the atomic partial charges :math:`q_i`, where :math:`i` are the atoms. This method comprises a second order expansion of the molecular energy :math:`E_{\mathrm{EEM}}` in terms of the partial charges :math:`q_i`:

.. math:: E_{\mathrm{EEM}}(q_1, q_2, \dots, q_N) = 
   \sum_{i=1}^{N} \Bigg[ \chi_i q_i + \frac{1}{2} \eta_i q_i^2 + \frac{1}{2} \sum_{\substack{l=1\\ l \ne i}}^{N} q_i q_l J_{il} \Bigg],
   :label: eem

where :math:`N` is the number of atoms, :math:`\chi_i` are the atomic electronegativities, :math:`\eta_i` the atomic hardness, and :math:`J_{il}` the Coulomb interaction between atoms. The factor :math:`\frac{1}{2}` before the Coulomb matrix is to avoid double counting.

In this work, the SQE method is used, which addresses a shortcoming of the EEM, namely that molecules tend to get over-polarized :cite:p:`Nistor2006a`. For more background, we refer to the recent review on charge flow models by Jensen :cite:p:`Jensen2023a`.

Verstraelen and co-workers proposed the following variant of the molecular energy:

.. math:: E_{\mathrm{SQE}} = E_{\mathrm{EEM}} +\sum_{i,j}^{M}\left(\frac{1}{2}\Delta\eta_{ij} p_{ij}^2 + \Delta\chi_{ij}(q_i - q_j)\right)
   :label: sqe

where :math:`p_{ij}` corresponds to the (intramolecular) charge transfer  over bonds, :math:`\Delta\eta_{ij}` is the bond hardness and :math:`\Delta\chi_{ij}` is the bond electronegativity correction.
Therefore, the charge variables :math:`q_i` are replaced by charge-transfer variables :math:`p_{ij}` which are related by 

.. math:: q_i = \frac{q_{\mathrm{tot}}}{N} + \sum_{\substack{i,j\\ \text{bonds}}} p_{ij},
   :label: qi

where :math:`q_{\mathrm{tot}}` is the net charge on the compound and :math:`p_{ij} = -p_{ji}`. Although it is trivial to determine the partial charges :math:`q_i` from the charge transfer :math:`p_{ij}`, the reverse is not necessarily true.
As outlined by Chen {\em et al.} :cite:p:`Chen2008b`, the problem can be solved by expressing the energy in terms of the charge transfer variables.
By substituting :math:`J_{ii} = \eta_i` in Eq. :eq:`eem`, inserting Eq. :eq:`eem` into Eq. :eq:`sqe` and introducing :math:`M_x` as the number of bonds for species :math:`x`, we obtain:

.. math:: \begin{align}
   E_{\mathrm{SQE}} &= \sum_{n=1}^{N}\Bigg[\left(\frac{q_{\mathrm{tot}}}{N} + \sum_{m=1}^{M_n} p_{nm}\right)\Bigg(\chi_n + \frac{1}{2}\eta_n \left(\frac{q_{\mathrm{tot}}}{N} + \sum_{m=1}^{M_n} p_{nm}\right) \nonumber\\
   &\quad + \frac{1}{2}\sum_{\substack{l=1\\ l\neq n}}^{N}\left(\frac{q_{\mathrm{tot}}}{N} + \sum_{m=1}^{M_l} p_{lm}\right)J_{nl}\Bigg)\Bigg] \nonumber\\
   &\quad + \sum_{i,j}^{M}\left[\frac{1}{2}\zeta_{ij}p_{ij}^2+\Delta\chi_{ij}\left(\left(\frac{q_{\mathrm{tot}}}{N} + \sum_{m=1}^{M_i} p_{im}\right)-\left(\frac{q_{\mathrm{tot}}}{N} +\sum_{m=1}^{M_j} p_{jm}\right)\right)\right].
   \end{align}

The next step is to determine the :math:`p_{ij}` that minimize :math:`E_{\mathrm{SQE}}`. Since all summations run over atoms :math:`i,j,k,l`, we take the derivative with respect to :math:`p_{ij}` and equate it to zero:

.. math:: \begin{align}
   0 &= \frac{\partial E_{\mathrm{SQE}}}{\partial p_{ij}} \nonumber\\
   &= \left(\chi_i - \chi_j + \frac{q_{\mathrm{tot}}}{N}(\eta_i - \eta_j) + \sum_{k=1}^{M_i}\Delta\chi_{ik} - \sum_{k=1}^{M_j}\Delta\chi_{jk}\right) \nonumber\\
   &\quad + \frac{1}{2}\left(\sum_{l=1}^{N}J_{il}\left(\frac{q_{\mathrm{tot}}}{N}+\sum_{m=1}^{M_l} p_{lm}\right) - \sum_{l=1}^{N}J_{li}\left(\frac{q_{\mathrm{tot}}}{N}+\sum_{k=1}^{M_i} p_{ik}\right)\right) + p_{ij}\zeta_{ij},
   \end{align}
   :label: deriv

using the identity of the Coulomb-matrix elements (:math:`J_{ij} = J_{ji}`), the terms involving the atomic hardness :math:`\eta_x` are incorporated into the diagonals of the :math:`J_{xy}` matrix, excluding the contribution from the total charge :math:`q_{\mathrm{tot}}`. Note that the :math:`q_{\mathrm{tot}}` terms in the two sums cancel. The first term in Eq. :eq:`deriv` is the difference in electronegativity between atoms :math:`i` and :math:`j` sharing the bond plus the correction; the second term represents the difference in electrostatic potentials at the atoms, and the third term accounts for the interaction between :math:`i` and :math:`j` times the charge transfer.
This results in a coupled set of equations, written in matrix form as

.. math:: \mathbf{M}\, \mathbf{P} = \mathbf{R},

where :math:`\mathbf{M}` is a square matrix of dimension equal to the number of bonds, :math:`\mathbf{P}` is the vector of the charge transfers for all bonds, and :math:`\mathbf{R}` is the right-hand side of the equations. The  matrix elements are given by

.. math:: M_{ij,kl} = J_{ik} - J_{il} - J_{jk} + J_{jl} + \delta_{ij,kl}\Delta\eta_{ij}

where :math:`\delta_{ij,kl}` is one if bond :math:`ij` is identical to :math:`kl` and zero otherwise.
The right-hand side is defined by the electronegativity terms according to

.. math:: R_{ij} = \chi_j - \chi_i + \sum_{k=1}^{M_j}\Delta\chi_{jk} - \sum_{k=1}^{M_i}\Delta\chi_{ik} + \frac{q_{\mathrm{tot}}}{N}\left(\sum_{l=1}^{N} J_{jl}-\sum_{l=1}^{N} J_{il}\right).


The charges of the shells (and virtual sites) are treated as constant in this algorithm, meaning that :math:`q_{tot}` becomes the sum of the charges of the shells and virtual sites and the total charge of the compound.
During force field training, all of these charges can be modified alongside the SQE parameters.
As noted in the paper describing the ACT software :cite:p:`Spoel2025b`, the SQE algorithm may not be flexible enough to reproduce optimal charge distributions, and other algorithms :cite:p:`Verstraelen2013a,Jensen2023a` may need to be implemented in future versions of the software.
