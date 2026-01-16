.. _sec-energy:

***************
Energy Function
***************
The ACT supports a range of different potentials for classical force field simulations. These potentials are described briefly in what follows. 

=======================
Non-bonded interactions
=======================

---------------------------------------
Coulomb Interaction using Point charges
---------------------------------------
The Coulomb interaction can be computed with point charges:

.. math:: V_{coul}(r_{ij}) ~=~ \frac{q_i q_j} {4\pi\varepsilon_0\varepsilon_r r_{ij}}.
   :label: vcoul

------------------------------------------
Coulomb Interaction using Gaussian charges
------------------------------------------
Alternatively, the Coulomb interaction can be described using Gaussian shielded charges

.. math:: V_{coul}(r_{ij}) ~=~ \frac{q_i q_j {\rm erf}\left(\zeta_{ij} r_{ij} \right)}{4\pi\varepsilon_0\varepsilon_r r_{ij}}, \quad \zeta_{ij} = \frac{\zeta_i \zeta_j}{\sqrt{\zeta_i^2 + \zeta_j^2}},
   :label: vcoulg

where :math:`\varepsilon_0` is the permittivity of vacuum, :math:`q_{i,j}` are the charges and :math:`\zeta_{i,j}` are the charge distribution widths (screening factors). Note that Eqn. :eq:`vcoulg` includes a relative dielectric constant :math:`\varepsilon_r` that can be used to parameterize a non-polarizable force field using charge-scaling :cite:p:`Kirby2019a`.

.. _sec-slater:

----------------------------------------
Coulomb Interaction using Slater charges
----------------------------------------
In addition, Slater distributed charges can be used as described in ref. :cite:p:`Ghahremanpour2018b`, from which the text below is an adapted version.

The spherical Slater orbital wavefunction is given by

.. math:: \psi_{n}({\bf r}) ~=~ \sqrt{\frac{(2\zeta)^{2n+1}}{4\pi(2n)!}} {\bf r}^{n-1}{\rm e}^{-\zeta {\bf r}}
   :label: SlaterWF

where :math:`n` is the highest  quantum number of the element and :math:`\zeta`
is the  exponent giving the width of the charge density. 
The distribution of  charge  can  be described  
by the charge density, which is the square of the wave function ( :eq:`SlaterWF` ). 
The Coulomb integral
can then be written as :cite:p:`Rick1994a,Hentschke2004a`:

.. math:: J_{ij}({\bf r}) \sim \int \int |\psi_n({\bf r_i})|^2\frac{q_i q_j}{|{\bf r_i}- {\bf r_j}|}|\psi_m({\bf r_j})|^2 d{\bf r_i} d{\bf r_j}
   :label: SlaterIntegral

where :math:`\psi_n` and :math:`\psi_m` are the Slater wave functions of atoms :math:`i` and
:math:`j` with quantum numbers :math:`n` and :math:`m` and with partial charges :math:`q_i` and :math:`q_j`,  respectively.
Both :eq:`vcoulg` and :eq:`SlaterIntegral` have a finite limit as :math:`r \rightarrow 0`.
As a result,  these functions are well behaved at small  :math:`r`. 
A number of methods :cite:p:`Ohrn2016a,Guseinov1970a` have been proposed to evaluate :eq:`SlaterIntegral`. 
Hentschke gives an  analytical solution :cite:p:`Hentschke2004a`:

.. math:: J_{ij}({\bf r}) = \left\{ \begin{aligned}
   & \frac{1}{4 \pi \epsilon_0}\frac{q_iq_j}{\|{\bf r_i} - {\bf r_j}\|} \frac{4\zeta_i^{2n+1}\zeta_j^{2m+1}}{(2n)!(2m)!} \frac{\partial^{2n-2}\partial^{2m-2}}{\partial\zeta_i^{2n-2}\partial\zeta_j^{2m-2}} \left(\frac{1}{\zeta_i^3 \zeta_j^3}\right)\\
   & \times\biggl[ 1 -\frac{(3\zeta_i^2 - \zeta_j^2)\zeta_j^4}{(\zeta_i - \zeta_j)^3(\zeta_i + \zeta_j)^3} e^{-2\zeta_ir_{ij}} -\frac{(\zeta_i^2 - 3\zeta_j^2)\zeta_i^4}{(\zeta_i - \zeta_j)^3(\zeta_i + \zeta_j)^3} e^{-2\zeta_jr_{ij}} \\
   & - \frac{\zeta_i\zeta_j^4}{(\zeta_i-\zeta_j)^2(\zeta_i+\zeta_j)^2} r_{ij} e^{-2\zeta_ir_{ij}} - \frac{\zeta_i^4\zeta_j}{(\zeta_i-\zeta_j)^2(\zeta_i+\zeta_j)^2} r_{ij} e^{-2\zeta_jr_{ij}}
   \biggr] \end{aligned} \right.
   :label: hentschke


:eq:`hentschke` was implemented in a Mathematica :sup`TM` program from which C++ code was generated for the
analytical computation of :math:`J_{ij}`  and its analytical derivatives with respect to :math:`r`, which are necessary for computing forces.
Due to the nature of :eq:`hentschke`, there are many terms with large powers, particularly for :math:`n>3`. Thus,  the equations have to be implemented using the arbitrary precision
arithmetic library Class Library for Numbers ( `CLN`_ ).
to avoid numerical instabilities.  
It should be noted that the  arbitrary  precision library significantly increases the  computational cost to analytically solve :eq:`hentschke`

.. _CLN: https://www.ginac.de/CLN/


-------------------------------
Long range Coulomb interactions
-------------------------------
It is good practice :cite:p:`Spoel2006a` to use the particle-mesh Ewald :cite:p:`Darden1993a,Essmann1995a` method to treat long-range electrostatics interactions in the condensed phase. In short, the method splits the Coulomb potential, either Eqn. :eq:`vcoul`, Eqn. :eq:`vcoulg` or something similar into a short and a long-range part as follows

.. math:: V_{coul}(r_{ij}) ~=~ V_{coul} \,{\rm erfc}(\alpha r_{ij}) +  V_{coul}\, {\rm erf}(\alpha r_{ij})
   :label: vcoullr

using the fact that the complementary error function {\rm erfc} equal 1 - {\rm erf}. The first term here is the short-range part, that is computed in real space, and the second term is the long-range part that is computed in Fourier (reciprocal) space by moving charges on a regular grid :cite:p:`Essmann1995a`.
The constant :math:`\alpha` determines how quickly the short-range potential decays to zero, and it can be computed from a user-specified relative error tolerance related to moving charges on a grid. Given the cut-off distance :math:`rc` and an error tolerance :math:`\epsilon`, typically :math:`10^{-4}`, we have

.. math:: \alpha ~=~ \frac{\sqrt{-{\rm ln}\,2\epsilon}}{rc}
   :label: alphaEwald

where ln is the natural logarithm.
This shows that :math:`\alpha` has the dimension of 1 over distance in the units of rc. Fig. :ref:`coulomb` shows how the division over short and long-range interactions works in practice, and how the long range contribution should be incorporated into simulations using a modified Coulomb function.

.. figure:: ../images/coulomb.pdf
   :name: coulomb
   :width: 90%
   :align: center

   Total Coulomb energy V as a function of distance r for two like unit charges, Gaussian-shielded Coulomb and short and long-range parts of the energy function. Ewald-tolerance :math:`\epsilon` 0.0001, rc 0.9 nm, :math:`\alpha` 3.24269/nm and Gaussian screening :math:`\zeta` 5/nm.

------------
Polarization
------------
Polarization can be treated explicitly as well in the ACT using the  core-shell model :cite:p:`Dick1958a,Jordan1995a` 

.. math:: V_{pol}(r_{cs}) ~=~ \frac{1}{4\pi\varepsilon_0} \frac{q_s^2}{2\alpha_c}r_{cs}^2,
   :label: vpol

where :math:`q_s` is the shell charge, :math:`\alpha_c` is the polarizability and :math:`r_{cs}` the distance between core and shell particles.  

.. _sec-indcorr:

--------------------
Induction correction
--------------------
\label{indcorr}
A term to add extra attraction between atoms was proposed in ref. :cite:p:`McDaniel2013a`

.. math:: V_{ic}(r_{ij}) ~=~ -A_{ij} {\rm e}^{-b_{ij} r_{ij}}
   :label: vic

with :math:`A_{ij}` and :math:`b_{ij}` constants to be trained. 
The geometric combination rule is used for
:math:`A_{ij}` and the arithmetic combination rule for :math:`b_{ij}`.
It is not certain whether the exponential functional form is optimal to describe this interaction and there is no theory behind this.

-------------------------------------
Pauli Repulsion and London Dispersion
-------------------------------------
The repulsion and dispersion interactions acting between atoms :math:`i` and :math:`j` can be described by a number of potentials. In a recent study on noble gases :cite:p:`Kriz2024a` we found that the 14-7 potential due to Halgren :cite:p:`Halgren1992a` was one of the most accurate ones:

.. math:: V_{14-7}(r_{ij}) ~=~ \epsilon_{ij} \left(\frac{1+\delta_{ij}}{\frac{r_{ij}}{\sigma_{ij}} + \delta_{ij}} \right)^7\left(  \frac{1+\gamma_{ij}}{(\frac{r_{ij}}{\sigma_{ij}})^7 +\gamma_{ij}}    -2      \right)       
   :label: 14_7

where  :math:`\epsilon_{ij}` is the well depth at minimum, :math:`\gamma_{ij}` and :math:`\delta_{ij}` are dimensionless numbers that were originally shared for all elements. However, in our previous work :cite:p:`Kriz2024a`, we treated :math:`\gamma` and :math:`\delta` as free atom-specific parameters subject to optimization and combination rules.

Furthermore, the generalized 4-parameter Buckingham potential with an adjustable long-range attraction is implemented as :cite:p:`Werhahn2015a` 

.. math:: V_{GBH}(r_{ij}) = \epsilon_{ij}\frac{\delta_{ij} + 2\gamma_{ij} + 6}{2\gamma_{ij}} \frac{1}{1+ (\frac{r}{\sigma_{ij}})^6 } \left[  \frac{6+\delta_{ij}}{\delta_{ij} +2\gamma_{ij}+6}  e^{\gamma_{ij} (1-\frac{r}{\sigma_{ij}})}  - 1  \right] -\frac{\epsilon_{ij}}{1+(\frac{r_{ij}}{\sigma})^{\delta}}
   :label: GBH

where :math:`\gamma` and :math:`\delta` again are dimensionless constants.


Another potential that is supported is the buffered (Wang) Buckingham potential :cite:p:`Wang2013a`:

.. math:: V_{WBH}(r_{ij}) = \left(\frac{2\epsilon_{ij}}{1-\frac{3}{\gamma_{ij}+3}}\right) \left( \frac{\sigma_{ij}^6}{\sigma_{ij}^6+r_{ij}^6} \right) \left[\frac{3}{\gamma_{ij}+3}e^{\displaystyle{\gamma_{ij}\left(1-\frac{r_{ij}}{\sigma_{ij}}\right)}} - 1 \right] 
   :label: vwbh

where :math:`\epsilon_{ij}` is the well depth at minimum, :math:`\sigma_{ij}` is the Van der Waals radius and :math:`\gamma_{ij}` determines the steepness of the potential and :math:`r_{ij}` is the distance between particles :math:`i` and :math:`j`. 
The buffered Buckingham potential mentioned above have been used to develop an accurate phase-transferable  model for alkali-halides :cite:p:`Walz2018a,Walz2019a,Walz2019b,Walz2019c,Walz2021a`
and the original Buckingham potential :cite:p:`Buckingham1938a` is implemented as

.. math:: V_{BH}(r_{ij}) ~=~ A_{ij}{\rm e}^{\displaystyle -b_{ij}r_{ij}} - \frac{C_{ij}}{r_{ij}^6}
   :label: vbh

where :math:`A_{ij}`, :math:`b_{ij}` and :math:`C_{ij}` are parameters to be trained.
An alternative way of formulating the Buckingham potential :cite:p:`mason_transport_2004` is

.. math:: E_{MBH} = \frac{\epsilon}{1-\frac{6}{\gamma}} \left[  \frac{6}{\gamma}  e^{\gamma (1-\frac{r}{\sigma})}  - \left(\frac{\sigma}{r}\right)^6  \right]
   :label: MBH

where :math:`\epsilon` is the well-depth :math:`\sigma` reflects the position of the minimum and :math:`\gamma` is again a dimensionless constant.
Although mathematically identical to Eqn. :eq:`vbh` we have shown that the potentials behave differently when combination rules are applied :cite:p:`Kriz2024a`.


The well-known Tang-Toennies potential :cite:p:`Tang1984a,Tang2003a` containing five parameters is implemented according to:

.. math:: V_{TT}(x) = Ae^{-b x} - \sum_{n=3}^5\left[1-e^{-b x}\sum_{k=0}^{2n}\frac{(b x)^k}{k!}\right]\frac{C_{2n}}{x^{2n}}
   :label: TT

where repulsion and dispersion terms shared the parameter :math:`b`. Here, all five parameters can be trained by the ACT.
This potential was used in our recent work on noble gases published where we investigated combination rules :cite:p:`Kriz2024a`.

In other works, the Tang-Toennies potential is extended with an additional parameter :cite:p:`Sheng2020a` giving the exponential functions a different decay length :math:`b`:

.. math:: V_{TT2}(x) = Ae^{-b_{exp} x} - \sum_{n=3}^5\left[1-e^{-b_{disp} x}\sum_{k=0}^{2n}\frac{(b_{disp} x)^k}{k!}\right]\frac{C_{2n}}{x^{2n}}
   :label: TT2

and the ACT supports both these functions.
Van Vleet {\em et al.} described a more accurate formula for the Pauli repulsion :cite:p:`Vleet2016a` that was applied to derive a force field later :cite:p:`Vleet2018a`. In this model the exchange repulsion is given by the Slater-ISA formula

.. math:: V_{exch} ~=~ A\left(\frac{1}{3}(br)^2 + br + 1\right){\rm exp}^{-br}\label{eqn:slater-isa}

where it can be noted, that no additional parameters are needed for equation~:eq:`slater-isa` compared to the repulsion part in the Tang-Toennies potential (Eqn. :eq:`TT`). In the ACT the Slater-ISA exchange can be combined with the damped dispersion part of the latter potential.

The Lennard-Jones 12-6 potential :cite:p:`Lennard-Jones1924b` is available as well:

.. math:: V_{LJ}(r_{ij}) ~=~ 4\varepsilon_{ij}\left(\frac{\sigma_{ij}^{12}}{r_{ij}^{12}}-\frac{\sigma_{ij}^{6}}{r_{ij}^6}\right)
   :label: 12_6

and, in addition, the 12-6-4 potential 

.. math:: V_{LJ}(r_{ij}) ~=~ 4\varepsilon_{ij}\left(\frac{\sigma_{ij}^{12}}{r_{ij}^{12}}-\frac{\sigma_{ij}^{6}}{r_{ij}^6}\right) - \frac{\gamma}{r^4}
   :label: 12_6_4

which is meant to model ion-dipole interactions is supported too.

It is worth noting that in all these potentials (except Buckingham, Eqn. :eq:`vbh`) parameters describe both the exchange and dispersion interactions at the same time, which necessitates simultaneous training of these interaction.
To account for anisotropy in exchange, a correction can be implemented using a virtual site particle on the atom that has a :math:`\sigma`-hole interacting with other particles :cite:p:`Kriz2024b`: 

.. math:: V_{EXCH,Corr}(r_{ij}) = -A_{ij}{\rm e}^{\displaystyle -b_{ij}r_{ij}}
   :label: exch_corr

where :math:`A_{ij}` and :math:`b_{ij}` are parameters to be optimized. This correction term (Eqn. :eq:`exch_corr`) does not need to be applied to all possible atom pairs, therefore a new combination rule, dubbed "Kronecker" is introduced

.. math:: A_{ij} ~=~ (1-\delta_{ij})\frac{A_i + A_j}{2}

where :math:`\delta_{ij}` is the Kronecker delta operating on particle types :math:`i` and :math:`j`. This means that only interactions between :math:`\sigma`-holes (represented by a virtual site) and atoms are non-zero. To avoid parameter explosion, the atomic :math:`A_i` are set to zero and only the virtual site :math:`A_j` is non-zero. This is reasonable since the virtual site represents a property of the :math:`\sigma`-hole of the atom it is connected to.
The arithmetic combination rule is used for :math:`b_{ij}`.

Multiple different combination rules can be used for parameters involved in van der Waals interactions, both within ACT and through the interface to OpenMM. For details, see the paper by K\v{r}{\'i}\v{z} {\em et al.} :cite:p:`Kriz2024a`.

----------------------------------
Treatment of long range dispersion
----------------------------------
The different Van der Waals potentials described above can be used with Lennard-Jones PME :cite:p:`Darden1993a,Wennberg2013a` if treated carefully. LJ-PME assumes a simple dispersion interaction given by

.. math:: V_{LJPME}(r_{ij}) ~=~ -\frac{C_{ij}}{r_{ij}^6}
   :label: ljpme

for atoms :math:`i` and :math:`j`, which differs from the potentials mentioned above. For LJ-PME, it is beneficial to use the geometric combination rule, hence

.. math:: V_{LJPME}(r_{ij}) ~=~ -\frac{C_iC_j}{r_{ij}^6}

where :math:`C_i` and :math:`C_j`, in contrast to Eqn. :eq:`ljpme`, have units of distance to the third power.

To minimize the difference between the potential used and the LJ dispersion over the whole volume outside the cut-off, and to specify correct inputs to OpenMM, starting from e.g. Eqn. :eq:`vwbh` we have to solve the following equation:

.. math:: 0 ~=~ \int_{rc}^{\infty} 4\pi r_{ij}^2\left[\left(-\frac{2\epsilon_{ij}}{1-\frac{3}{\gamma_{ij}+3}}\right) \left( \frac{\sigma_{ij}^6}{\sigma_{ij}^6+r_{ij}^6} \right) + \frac{C_iC_j}{r_{ij}^6} \right] {\rm d}r_{ij}
   :label: dispint

where we integrate from the cut-off :math:`rc` to infinity.
In this notation, we have to insert the combination rules for :math:`\epsilon`, :math:`\sigma` and :math:`\gamma`. For the special case that :math:`i` = :math:`j`,
the constant :math:`C_{ii}` can be derived using Mathematica :sup:`TM` as

.. math:: C_{ii} ~=~ \epsilon_{i}\frac{3+\gamma_{i}}{\gamma_{i}}rc^3\sigma_{i}^3\left(\pi-2 {\rm arctan}\left[\left(\frac{rc}{\sigma_{i}}\right)^3\right]\right)

which we can then convert back to a :math:`\sigma'` compatible with standard Lennard-Jones by

.. math:: \sigma'_{i} = \left(\frac{C_{ii}}{4\epsilon_{i}}\right)^{1/6}.

Note however, that we need to get effective :math:`\sigma'` for all atoms, which means the problem turns into a minimization problem where 

.. math:: \varepsilon^2 ~=~ \sum_{i,j}\left[ \int_{rc}^{\infty} 4\pi r_{ij}^2\left[\left(-\frac{2\epsilon_{ij}}{1-\frac{3}{\gamma_{ij}+3}}\right) \left( \frac{\sigma_{ij}^6}{\sigma_{ij}^6+r_{ij}^6} \right) + \frac{C_iC_j}{r_{ij}^6} \right] {\rm d}r_{ij}\right]^2
   :label: dispintsum

has to be minimized with respect to :math:`C_i` independently for each set of combination rules. In the Alexandria alkali-halide model we used the Wang-Buckingham potential (Eqn. :eq:`vwbh`) in conjunction with the combination rules due to Kong :cite:p:`Kong1973a` in which :math:`\sigma_{ij}` depends on all the :math:`\sigma`, :math:`\epsilon` and :math:`\gamma`. This makes it cumbersome to analytically derive integrals like Eqn. :eq:`dispintsum` for many combinations of potentials and combination rules.

Another problem is, that the parameters :math:`C_i` that we need to compute the long-range dispersion interaction now depend on the cut-off, which means they should preferably not be computed beforehand. In summary, the :math:`C_i` should be derived numerically given the potential and combination rules at the start of a simulation.


----------------------------------------------
Combination Rules for Van der Waals Potentials
----------------------------------------------
In a recent paper :cite:p:`Kriz2024a` we studied the effect of combination rules on potentials for noble gases and introduced two new combination rules (Eqn. :eq:`cr_volumetric` and Eqn. :eq:`cr_invsquare` ). Here is a brief summary, copied (with permission) from that paper.

Combination rules reduce the number of parameters required for the pair-wise potentials introduced  above because
the parameters describing the interaction between dissimilar  X--Y atoms are reconstructed from parameters of homodimers X--X and Y--Y. In this way, it is necessary only to fit atomic parameters on data for homodimers.
Different sets of mathematical expression were considered, as detailed below.

The two simplest expressions that have historically been used as combination rules are the geometric :cite:p:`Berthelot1898a` and arithmetic :cite:p:`Lorentz1881a` averages, respectively:

.. math:: X_{12} = \sqrt{x_{1} x_{2}}
   :label: cr_geometric

%was used for all three parameters :math:`\epsilon`, :math:`\gamma`, :math:`\sigma` in "geometric" rules. It is also used for :math:`\epsilon` in "arithmetic" rules and for :math:`\sigma` in "Kong-Mason" rules \cite{}.

.. math:: X_{12} = \frac{x_{1}+ x_{2}}{2}
   :label: cr_arithmetic

where :math:`x_1` and :math:`x_2` are the atomic parameters. Both these rules can in principle be applied to all parameters in the Van der Waals potentials described above and the rules are well-behaved mathematically.

Hogervorst introduced a set of combination rules for 12-6 Lennard-Jones and exp-6, modified Buckingham, potential (Eqn. :eq:`MBH`) :cite:p:`Hogervorst1971a`. He proposed using Eqn. :eq:`cr_hogervorst_epsilon` for :math:`\epsilon` for both potential forms. For the :math:`\sigma` of 12-6 potential (Eqn. :eq:`12_6`), the
arithmetic mean (Eqn. :eq:`cr_arithmetic`) was used. In the case of the modified Buckingham potential (eqn. :eq:`MBH`), expression~:eq:`cr_sigma5`, which depends on the combined parameters :math:`\epsilon_{1,2}` and :math:`\gamma_{1,2}` was advocated for the :math:`\sigma`, along with arithmetic mean for :math:`\gamma`. 

.. math:: X_{12} = \frac{2 x_{1} x_{2}}{x_{1} + x_{2}}
   :label: cr_hogervorst_epsilon

.. math:: \sigma_{12}^6 = \sqrt{\frac{\epsilon_1\gamma_1\sigma_1^6}{\gamma_1 - 6}  \frac{\epsilon_2\gamma_2\sigma_2^6}{\gamma_2 - 6}} \frac{(\gamma_{1,2}-6)}{\gamma_{1,2}\epsilon_{1,2}}
   :label: cr_sigma5

Where the :math:`\gamma`:sup:`1,2` used Eqn. :eq:`cr_arithmetic` and :math:`\epsilon`:sup:`1,2` Eqn. :eq:`cr_hogervorst_epsilon`. 
It should be noted that equation :eq:`cr_hogervorst_epsilon` is ill-behaved if both :math:`x_1` and :math:`x_2` are zero, while Eqn. :eq:`cr_sigma5`
is ill-behaved if either :math:`\gamma_{1,2}` or :math:`\epsilon_{1,2}` is zero.
Yang {\em et al.} :cite:p:`Yang2018_combination` introduced an expression  for the Morse potential :cite:p:`Morse1929a`, using Eqn. :eq:`cr_hogervorst_epsilon` for :math:`\epsilon`, while eqn. :eq:`cr_yang` is used for :math:`\sigma` and  :math:`\gamma`. 

.. math:: X_{12} = \frac{x_{1} x_{2}(x_1 + x_2)}{x_{1}^2 + x_{2}^2}
   :label: cr_yang

The expression for :math:`\gamma`  proposed by Mason :cite:p:`Mason1955a` for the exp-6 potential has the form
%\red{while used other rels for epsilon, sigma in alexandria}

.. math:: \gamma_{12} = \sqrt{\sigma_1 \sigma_2} \left(\frac{\gamma_1}{2\sigma_1}+  \frac{\gamma_2}{2\sigma_2}\right)
   :label: f10

Waldman and Hagler :cite:p:`Waldman1993a` introduced expressions :eq:`cr_whepsilon` for :math:`\epsilon` and  :eq:`cr_whsigma` for :math:`\sigma`
 to reproduce experimental well-depths and interaction distances.

.. math:: \epsilon_{12} = \sqrt{\epsilon_1  \epsilon_2}\frac{2\sigma_1^3\sigma_2^3}{\sigma_1^6+\sigma_2^6}\\
   :label: cr_whepsilon

.. math:: X_{12} = \left(\frac{X_1^6+X_2^6}{2}\right)^{1/6}
   :label: cr_whsigma

Although Eqn. :eq:`cr_whsigma` was devised for :math:`\sigma` we have evaluated it for other parameters here as well, hence the notation with :math:`X`. Qi and coworkers advocated the use of buffered 14-7 Lennard-Jones potential (Eqn. :eq:`14_7`) due to Halgren :cite:p:`Halgren1992a`, alongside  combination expressions :eq:`cr_whepsilon` for :math:`\epsilon` and  :eq:`f7` for :math:`\sigma`: :cite:p:`qi_general_2016` 

.. math:: X_{12} = \frac{x_1^3+x_2^3}{x_1^2+x_2^2}
   :label: f7

A further relation, the harmonic mean rule, was proposed by Halgren :cite:p:`Halgren1992a`:

.. math:: X_{12} =  \frac{4 x_1 x_2}{(x_1^{1/2}+x_2^{1/2})^2}
   :label: cr_harmonic_mean


Finally, we introduce two new combination rules that we have not seen published previously. Since the :math:`\sigma` in most potentials can be interpreted as a Van der Waals radius, we introduced a relation averaging  third powers, corresponding to an atomic volume :cite:p:`Kriz2024a`:

.. math:: \sigma_{12} = \left(\frac{\sigma_1^3+\sigma_2^3}{2}\right)^{1/3}.
   :label: cr_volumetric

In addition, we applied the following rule for in particular :math:`\epsilon` since it yield an :math:`X_{12}` that is smaller than the geometric one (Eqn. :eq:`cr_geometric`):

.. math:: X_{12} = \left(\frac{2}{x_1^{-2}+x_2^{-2}}\right)^{1/2}
   :label: cr_invsquare

The combination relations described above were permuted with each other into new combination rules. In this way, relations that depend on only one parameter type were used for any parameter. Relations depending on multiple parameters  were used only for the specific parameter type combination they depend on (e.g. Eqn. :eq:`cr_whepsilon` was only used for :math:`\epsilon`, using homodimer :math:`\epsilon` and :math:`\sigma`). 
In our previous work on alkali halides :cite:p:`Walz2018a` we used combination rules according to eqn. :eq:`cr_sigma5` for :math:`\sigma`, eqn. :eq:`cr_hogervorst_epsilon` for :math:`\epsilon` and eqn. :eq:`cr_arithmetic` for :math:`\gamma` with the Wang-Buckingham potential (Eqn. :eq:`vwbh`).

A number of these combination rules can be written using the generalized mean equation :cite:p:`Hohm2025a`

.. math:: M_p(x_1,x_2,...,x_N) = \left(\frac{1}{N}\sum_{i=1}^N x_i^p\right)^{1/p}
   :label: cr_genmean

for :math:`p \ne 0`. For :math:`p = 0`, Eqn. :eq:`cr_genmean` turns into the geometric rule. Table :ref:`tab-genmean` lists other well-known combination rules and their respective exponent :math:`p`. Hohm also describes combinations of the generalized mean and other similar expressions, but the most important observation he made was that the exponent :math:`p` can be varied at will :cite:p:`Hohm2025a`. In the ACT it is possible to train this parameter along with the Van der Waals parameters.

.. table:: Correspondence between well-known combination rules and the generalized mean equation :eq:`cr_genmean`.
   :name: tab-genmean
   
   +------+-------------------------------+-----------------------------+
   | p    |  Name                         | Equation                    |
   +======+===============================+=============================+
   |  -2  | Inverse Square                | :eq:`cr_invsquare`          |
   +------+-------------------------------+-----------------------------+
   |  -1  | Hogervorst :math:`\epsilon`   | :eq:`cr_hogervorst_epsilon` |
   +------+-------------------------------+-----------------------------+
   | -1/2 | Harmonic Mean                 | :eq:`cr_harmonic_mean`      |
   +------+-------------------------------+-----------------------------+
   |   0  | Geometric                     | :eq:`cr_geometric`          |
   +------+-------------------------------+-----------------------------+
   |   1  | Arithmetic                    | :eq:`cr_arithmetic`         |
   +------+-------------------------------+-----------------------------+
   |   3  | Volumetric                    | :eq:`cr_volumetric`         |
   +------+-------------------------------+-----------------------------+
   |   6  | Waldman-Hagler :math:`\sigma` | :eq:`cr_whsigma`            |
   +------+-------------------------------+-----------------------------+

===================
Bonded interactions
===================

------------------
Harmonic potential
------------------
Bond vibrations can be described using a harmonic term based on the bond length r:math:`_{ij}`

.. math:: V_b(r_{ij}) ~=~ \frac{k_{ij}^{b}}{2}\left(r_{ij}-r_{ij}^0\right)^2,

where :math:`k_{ij}^{b}` is the force constant, and r:math:`_{ij}^0` is the equilibrium bond length. 

------------------
Morse potential
------------------
A Morse potential :cite:p:`Morse1929a` can be used, with one addition:

.. math:: V_M(r_{ij}) ~=~ D_{ij}^e\left[{\rm e}^{\displaystyle{-2\beta_{ij}(r_{ij}-r_{ij}^0)}}-2{\rm e}^{\displaystyle{-\beta_{ij}(r_{ij}-r_{ij}^0)}}\right] + D_{ij}^0
   :label: morse

The term :math:`D_{ij}^e` roughly corresponds to a dissociation energy, however since there are Coulomb and/or Buckingham interactions between the atoms as well, the total "bond" potential is given by the sum of three terms and a correction term :math:`D_{ij}^0` is needed to get the correct energy minimum. 

------------------
Wei-Hua potential
------------------
In a very recent study we found the potential due to Hua :cite:p:`Hua1990a,Hua1990b` to be the best compromise between accuracy of the vibrational frequencies and computational cost :cite:p:`Maaren2025a`. 
It is given by:

.. math:: U(r) ~=~ D_e\left(\left[ \frac{1-e^{-b(r-r_e)}}{1-ce^{-b(r-r_e)}} \right]^2 - 1 \right)
   :label: hua

where :math:`D_e` is the well-depth, :math:`r_e` the equilibrium bond length, and :math:`b` and :math:`c` are constants with :math:`\|c\| < 1`. 

---------------
Angle potential
---------------
Angle vibrations are described using a harmonic term based on the angle :math:`\theta_{ijk}`

.. math:: V_a(\theta_{ijk}) ~=~ \frac{k_{ijk}^{\theta}}{2}\left(\theta_{ijk}-\theta_{ijk}^0\right)^2,

where :math:`k_{ijk}^{\theta}` is the force constant, and :math:`\theta_{ijk}^0` is the equilibrium angle. 

------------------------------------
Angle potential for linear compounds
------------------------------------
The reference position, corresponding to a minimum energy structure, :math:`\mathbf{x}_j^0` for a central atom :math:`j` in a linear triplet of atoms 
:math:`i,j,k` is given by

.. math:: \mathbf{x}_j^0 ~=~ {\rm a}\, \mathbf{x}_i + (1-{\rm a})\, \mathbf{x}_k 
   :label: xxx

where :math:`a` is a constant defined by the bond-lengths :math:`i-j` and :math:`j-k`.  In a group with bonds :math:`i-j` and :math:`j-k` with lengths :math:`b_{ij}` and :math:`b_{jk}` respectively, the constant is

.. math:: a ~=~ \frac{b_{jk}}{b_{ij}+b_{jk}}.
   :label: eqna

If the order of atoms is flipped  :math:`a` will change to :math:`1-a`. The potential :math:`V_{lin}` is then given by

.. math:: V_{lin} ~=~ \frac{{\rm k}_{lin}}{2} \left(\mathbf{x}_j - \mathbf{x}_j^0\right)^2 
   :label: vlin

with k :sub:`lin` the force constant :cite:p:`Spoel2020a`. 

-----------------------
Out-of-plane vibrations
-----------------------
Finally, out-of-plane vibrations are treated by another harmonic potential

.. math:: V_i(\phi_{ijkl}) ~=~ \frac{k_{ijkl}^{\phi}}{2}\phi_{ijkl}^2,

where :math:`k_{ijkl}^{\phi}` is the force constant and :math:`\phi_{ijkl}` is defined by the angle between the two planes :math:`i,j,k` and :math:`j,k,l`.

------------------
Torsion potential
------------------
A torsion potential is implemented using a Fourier series:

.. math:: V_{d}(\phi_{ijkl}) ~=~ \sum_{n=0}^5 c_n {\rm cos}^n(\pi+\phi_{ijkl})

where :math:`c_n` are constants and the torsion angle is defined as above. The constant :math:`\pi` is added to be compatible with the Ryckaert-Bellemans potential :cite:p:`Ryckaert1975a` that is implemented in simulation codes like GROMACS :cite:p:`Spoel2005a` and OpenMM :cite:p:`Eastman2010a`.

==================
Special potentials
==================
The ACT includes a flat-bottom position restraint potential according to

.. math::
   \begin{align}
    V_{fbpr}(r) &=& 0 &\hspace{2cm} r <= r_0\\
                &  & \frac{k}{2}(r-r_0)^2 &\hspace{2cm} r > r_0
   \end{align}

where :math:`k` is the force constant and :math:`r_0` the radius of the sphere (centered at the origin) in which the potential is zero.
The flat-bottom potential is activated by flags to the *alexandria simulate* command. It is useful
mainly to keep molecules close to the origin and prevent them from flying into outer space.

==================
Virtual Sites
==================
A virtual site is an extra point, located at a defined position in a molecule. A variety of virtual site options is currently implemented within the ACT framework:

  * a virtual site on top of an atom (VSITE1) :cite:p:`Spoel2025a`
  * a virtual sites along the bond (VSITE2) for the description of anisotropic charge distribution and exchange :cite:p:`Kriz2024b` such as encountered in :math:`\sigma` holes
  * a virtual site on the bisector of a angle, like in the TIP4P water model :cite:p:`Jorgensen1983a` (VSITE3S) or in alcohol (VSITE3)
  * off-plane virtual sites for modeling lone-pairs in sp:math:`^3` hybridized compounds, such as water (VSITE3SOUT, symmetric) or asymmetric for compounds like alcohols :cite:p:`Mahoney2000a,Kriz2024b`.
  * Four-particle virtual sites (VSITE4) to model a lone-pair on an amine group.


==================
Total energy
==================
The total energy :math:`E` of a compound then follows from 

.. math::
    E ~=~ V_{vdw}(r_{ij})+V_{coul}(r_{ij})+V_{pol}(r_{cs})+V_b(r_{ij})+V_a(\theta_{ijk})+V_i(\phi_{ijkl}) + V_{d}(\phi_{ijkl}).
   :label: energy
   
Finally, it should be noted that the number of excluded neighbors is user-configurable. That means that atoms that are covalently bonded can interact both through the Buckingham (or Lennard-Jones) and Coulomb potentials, and through the bonded potentials. The main reason for this is that the
short-range Coulomb interactions yield polarization anisotropy that is difficult to reproduce by a non-interacting model.
To make sure that the forces on the atoms in a molecule are zero in the reference  minimum-energy structure from quantum chemistry, both bond lengths :math:`r_{ij}^0` and angles :math:`\theta_{ijk}^0` can be treated as free parameters, that may differ substantially from the reference geometry. 
The number of exclusions can be selected separately for Coulomb and Van der Waals forces.

In total there are up to seven "atom" parameter types (:math:`\epsilon`, :math:`\sigma`, :math:`\gamma`, :math:`\delta`, :math:`\zeta`, :math:`q_s`,  and :math:`\alpha`) and 7-14 "bond" parameter types (:math:`k^b`, :math:`D^e`, :math:`r^0`, :math:`r^{max}`, :math:`k^{\theta}`, :math:`\theta^0`, :math:`k^{lin}` and :math:`c_{n}`) where :math:`n` is the dihedral term index running from 0 to 6 to determine. 
All the atomic parameters are taken to be hybridization state dependent (corresponding to, for instance, sp :sup:`1`, sp :sup:`2` and sp :sup:`3` carbon atoms).
%It is straightforward to add support for other potentials, such as proper dihedrals.

