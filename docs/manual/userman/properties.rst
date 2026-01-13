*******************************
Predicting Molecular Properties
*******************************
\label{properties}
The ACT can be used to perform MD simulations of clusters in the gas-phase using the \actcmd{simulate} command. This module includes the possibility to perform energy minimizations with the \actflag{-minimize}.
For simulations employing periodic boundaries the OpenMM package~\cite{Eastman2023a} should be used instead. 
For more details about MD simulations, see Section~\ref{simulations}.

Below we describe some of the properties that can be computed using the ACT.

=======================
Electrostatic Potential
=======================
\label{sec:esp}
The charge distribution $\rho$ of a molecule is determined by the nuclear position of the $N$ atoms at positions $\mathbf{x}$ and the electron density $n(\mathbf{r})$. The molecular electrostatic potential (MEP) at a point in space $r^\prime$ is thus given by
\begin{equation}
    \Phi(\mathbf{r}^\prime) = \frac{1}{4\pi\varepsilon_0}\left[\sum_{i=1}^N \frac{z_i}{ \|\mathbf{x}_i-\mathbf{r}^\prime\|} - \int \frac{n(\mathbf{r})}{\|\mathbf{r}-\mathbf{r}^\prime\|}{\rm d}r \right]
\label{phi}
\end{equation}
where $N$ is the number of atoms, $z_i$ are the nuclear charges, $\varepsilon_0$ is the permittivity of vacuum and integration is over the entire space. The minus sign before the integral is due to the negative charge of electrons. 

Accurate knowledge of the MEP contributes to, for example, the understanding of interactions and function of biological macromolecules in solution~\cite{Larsson2012b}. For a molecule in the gas phase, Eqn.~\ref{phi} can be evaluated using density functional theory and wave function quantum chemistry, albeit at a significant computational cost. Databases of such calculations for small molecules are available to facilitate reuse~\cite{Kriz2023a}. For  large condensed-phase systems, however, it is common to apply classical force fields, where electrons are not taken into account explicitly. Instead, effective partial charges on atoms are used. The electronic degrees of freedom, charge polarization, is sometimes taken into account through induced point dipoles and higher electrostatic moments, or by using a core-shell model~\cite{Dick1958a,Jordan1995a,Maaren2001a}. For additional background we refer to some excellent reviews~\cite{Dauber-Osguthorpe2019a,Hagler2019a,Jing2019a}.

The MEP can be used as a target in model development in the ACT, however we recommend against that for both fundamental and practical reasons~\cite{Hosseini2025a}.

=====================
Electrostatic Moments
=====================
If point ${\bf r}^\prime$ in \eqnref{phi} is outside the distribution of electron
density and ${\bf r}^\prime >> {\bf r}$, the electrostatic potential can be evaluated through the  
Taylor expansion of $|{\bf r}-{\bf r}^\prime|^{-1}$ \cite{Berendsen2007a}:

\begin{equation}
\frac{1}{|{\bf r^\prime}-{\bf r}|} \approx \frac{1}{r} + \left({\bf \hat r} \cdot
  {\bf r} \right)\frac{1}{r^2} + \frac{1}{2} \left[3\left({\bf \hat r} \cdot
  {\bf r} \right)^2 - r^{2}{\rm \bf I}\right]\frac{1}{ r^3} + \cdots
\label{eq:Rexpansion}
\end{equation}
where ${\bf \hat r} = {\bf r}^\prime/r$ and ${\rm \bf I}$ is the
identity matrix. By inserting \eqnref{eq:Rexpansion} into
\eqnref{phi}, we get

\begin{equation}
  \Phi({\bf r}^\prime) \approx \frac{1}{4 \pi \epsilon_0} { \frac{Q}{r}
    + \frac{{\bf \mu_0}}{r^2} + \frac{\boldsymbol{\Theta_0} }{r^3} + \cdots}
\end{equation}
$Q$ is the monopole moment, sometimes called the zeroth
moment of the molecular electron density.  In principle this is the total charge of the molecule. In what follows we combine both the atomic cores and the electron density into $n(r)$. Then, $Q$ is given by
\begin{equation}
 Q = \int n(\mathbf{r}) d {\bf r}  
\end{equation}
${\bf \mu_0}$ is the vector of the permanent dipole moment, which measures 
the polarity of the molecular electron density. ${\bf \mu_0}$ is given by
\begin{equation}
 {\bf \mu_0} = \int n(\mathbf{r}) \left({\bf \hat r} \cdot {\bf r}\right) d {\bf r} 
\end{equation}
$\boldsymbol{\Theta_0}$ is the tensor of the permanent quadrupole moment,  which exhibits 
the deviation of the distribution of the molecular electron density
from spherical symmetry. $\boldsymbol{\Theta_0}$ is given by
\begin{equation}
 \boldsymbol{\Theta_0} = \frac{1}{2}\int n(\mathbf{r}) \left[3\left({\bf \hat r}
     \cdot {\bf r}\right)^2 - r^2 {\rm \bf I} \right] d {\bf r} 
\end{equation}
The quadrupole tensor can be written as a traceless $3 \times 3$ matrix 
in the Cartesian coordinate if one writes $\left({\bf \hat r} \cdot
  {\bf r} \right)^2$ as:
\begin{eqnarray}
\left({\bf \hat r} \cdot {\bf r} \right)^2 = {\bf \hat r} \cdot \left({\bf r} {\bf r} \right) \cdot {\bf \hat r}
\end{eqnarray}
where ${\bf rr}$ is the outer product of vector ${\bf r}$ with itself. This
results in a matrix that can be written in terms of the Cartesian 
components of the vector ${\bf r}$ as follows:
\begin{equation}
{\bf rr} = 
\left[\begin{array}{cccc}
 x^2 & xy   & zx \\
 yx  &  y^2 & yz \\
 zx  &  zy   & z^2 \\
\end{array}
\right]
\end{equation}
Finally, we get
\begin{align}
\frac{1}{2} \left[3\left({\bf \hat r} \cdot
  {\bf r} \right)^2 - r^{2} {\rm \bf I}\right] = {\bf \hat r} \cdot
\frac{1}{2}\left[\begin{array}{cccc}
3x^2- r^2 & 3xy   & 3zx \\
3yx  &  3y^2- r^2 & 3yz \\
3zx  &  3zy   & 3z^2- r^2 \\
\end{array}
\right]
\cdot
{\bf \hat r} 
\end{align}
\iffalse
Therefore, the components of $\boldsymbol{\Theta_0}$ can be computed as follows:
\begin{equation}
\Theta_{ij} = \sum_{i=x,y,z} \sum_{j=x,y,z}  \frac{1}{2} \left(3 {\bf r}_{i} {\bf r}_{j} - r^2 \delta_{ij}\right)
\end{equation}
where $\rm r^2 = |{\bf r}|^2$ and $\delta$ is the Kronecker delta
function. 
\fi

==============
Polarizability
==============
The shape of the molecular electron density changes when it interacts with an external electric
field; hence, the total energy of the molecule
changes. The static response of a molecule to a homogeneous external
electric field, ($\bf F$), can be studied by expanding its energy in a
Taylor series \cite{Jensen2007a, Calaminici1998a}: 

\begin{eqnarray}
E \left( {\bf F} \right) = E\left({\bf 0}\right) + \frac{\partial E}{\partial {\bf F}}\biggm\lvert_{{\bf F}=0}
                      {\bf F} + \frac{1}{2}
                     \frac{\partial^2E}{\partial {\bf F}^2 }\biggm\lvert_{{\bf F}=0}
                     {\bf F}^2 + \frac{1}{6}
                     \frac{\partial^3E}{\partial {\bf F}^3 }\biggm\lvert_{{\bf F}=0}
                     {\bf F}^3 + \cdots
\end{eqnarray}
where
\begin{equation}
-\frac{\partial E}{\partial {\bf F}}\biggm\lvert_{{\bf F}=0}  = {\bf \mu_0}
\end{equation}
\begin{equation}
- \frac{\partial^2E}{\partial {\bf F}^2}\biggm\lvert_{{\bf F}=0} = {\bf \alpha}
\end{equation}
\begin{equation}
- \frac{\partial^3E}{\partial {\bf F}^3}\biggm\lvert_{{\bf F}=0} = {\bf \beta}
\end{equation}
where ${\bf \mu_0}$ is the vector of permanent dipole moment, ${\bf \alpha}$
is the tensor of polarizability, which is the 
linear part of the response of the molecular electron density 
with respect to the external electric field, and $\boldsymbol{\beta}$ is the
first hyperpolarizability \cite{Jensen2007a}. 

Instead of expanding the energy, we can expand the dipole moment of a molecule in an
external electric field \cite{Calaminici1998a}, written as
\begin{equation}
\mu = \mu_0 + {\bf \alpha} {\bf F} +
\frac{1}{2} {\bf \beta} {\bf F}^2 + \cdots
\end{equation}
where ${\bf \alpha} {\bf F}$ gives the vector of induced dipole
moment, ${\bf \mu_1}$ \cite{Jensen2007a}: 
\begin{equation}
{\bf \mu_1} = {\bf \alpha}{\bf F}
\end{equation}
that can be written in matrix form as
\begin{equation}
\left[\begin{array}{c}
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
\label{eqn:matrix2}
\end{equation}
From the polarizability tensor the polarizability isotropy
\cite{Stone2013a,Calaminici1998a},
\begin{equation}
\bar \alpha = \frac{\left(\alpha_{xx} + \alpha_{yy} + \alpha_{zz} \right)}{3}
\end{equation}
and the polarizability anisotropy \cite{Stone2013a}, 
\begin{equation}
\Delta \alpha = \sqrt{[(\alpha_{xx} - \alpha_{yy})^2 + (\alpha_{xx} -
  \alpha_{zz})^2 + (\alpha_{yy} - \alpha_{zz})^2 + 6 (\alpha_{xy}^2 +
  \alpha_{xz}^2 + \alpha_{yz}^2)]/2}
\end{equation}
can be calculated.

====================
Normal Mode Analysis
====================
\label{sec:nma}
Please note that the text below is largely taken (with permission) from a paper by Henschel {\em et al.}~\cite{Henschel2020a}.

The ACT contains the \actcmd{nma} tool that performs a normal mode analysis to determine the vibrational frequencies of a compound.
Vibrational frequencies are required to compute the IR spectra and thermochemistry of molecules. The normal modes of molecular vibrations can be obtained by eigenvalue decomposition of the Hessian matrix, whose elements are the second derivatives of the energy with respect to the atomic coordinates $q$.

\begin{equation}
H_{ij} = \frac{\partial^2E}{\partial q_i \partial q_j}
\end{equation}

where $i$ and $j$ run from $0$ to $N-1$, where $N$ is the number of atoms in the molecule. If virtual sites $v$ are used, for example, to model the $\sigma$-hole for halogen atoms, the energy,
E, depends on the positions of both atoms and virtual sites;
that is, $E= E(q_0,...,q_{N-1},v_0,...,v_{M-1})$, where the positions of the $M$
virtual sites, $v$, in the compound are a function of the atomic
coordinates q. 

The Hessian is computed numerically---the
N atoms are moved independently in all three spatial dimensions,
and the forces are computed. From these forces, the second
derivative of the energy is then evaluated numerically. Note that the positions of the virtual sites are
updated before each force calculation, which means that their influence on the Hessian is taken into account explicitly when computing $H$.

----------------
Infrared Spectra
----------------
\label{sec:irspectrum}
For the calculation of a full IR spectrum, in addition to the vibrational frequencies, the intensities and the line shapes are required.

In case of the quantum chemical calculations both the eigenfrequencies and the corresponding IR intensities are produced by default when a frequency calculation is requested in, for instance, the Gaussian software~\cite{g16}. The frequencies for about 5000 compounds  are available from the Alexandria library~\cite{Ghahremanpour2018a} at the B3LYP/aug-cc-pvtz level of theory~\cite{Becke1988a,Kendall1992a,Woon1993a,Woon1993b}.Details of the quantum chemical calculations from which the frequencies were obtained have been presented previously (refs.~\citenum{Ghahremanpour2016a,Ghahremanpour2018a}).

For the force field calculations, the intensities $I_{n}$ were derived from the transition dipole derivatives:
\begin{equation}
I_{n}=\sum_{k=1}^{3}\left(\frac{\partial p_{k}}{\partial Q_{n}}\right)^2
\label{In}
\end{equation}
\noindent where $k$ iterates over cartesian dimensions, $p$ is the dipole moment of the molecule, and $Q_n$ the normal coordinate $n$. In order to take into account virtual sites $v$ we note that 
\begin{equation}
p_{k} = p_{k}(q_0,...,q_{N-1},v_0,...,v_{M-1})
\end{equation}
and rewrite equation~\ref{In} as:
\begin{equation}
I_{n}=\sum_{k=1}^{3}\left(\frac{\partial p_{k}}{\partial q_s}\frac{\partial q_{s}}{\partial Q_{n}}\right)^2
\end{equation}
where $s$ iterates over the $N$ atomic coordinates. 
To make the calculation of intensities practical, the numerical derivative of the dipole moment with respect to the atomic coordinates $\displaystyle{\frac{\partial p_{k}}{\partial{q_s}}}$ is stored in a text file during the normal mode analysis and finally we note that the term $\displaystyle{\frac{\partial q_{s}}{\partial Q_{n}}}$ corresponds to one over component $s$ of eigenvector $n$.

As an example, an infrared spectrum can be generated using a force field file and a molprop file using
\begin{verbatim}
    alexandria nma -ff OPLS2020 -charges OPLS2020-charges
    -db ethanol -ir ir-ethanol
\end{verbatim}
yielding the spectrum in Fig.~\ref{fig:ethanol}.

\begin{figure}
    \centering
    \includegraphics[width=0.9\linewidth]{images/ethanol-irspectrum-opls2020.pdf}
        \caption{Simulated infrared vibrational spectrum for ethanol, based on the OPLS2020 force field~\cite{Jorgensen2023a}.}
    \label{fig:ethanol}
\end{figure}

---------------
Thermochemistry
---------------
The \href{(https://en.wikipedia.org/wiki/Canonical_ensemble)}{canonical ensemble} $Q(N,V,T)$ can be used to compute the molecular \href{(https://en.wikipedia.org/wiki/Internal_energy)}{internal energy}, the \href{https://en.wikipedia.org/wiki/Standard_molar_entropy}{standard entropy}, and the \href{https://en.wikipedia.org/wiki/Heat_capacity}{heat capacity} at constant volume. 

\begin{equation}
E = RT^2\left(\frac{\partial {\rm ln} Q}{\partial T} \right)_{N,V}
\label{intE}
\end{equation}

\begin{equation}
C_v = 2RT\left(\frac{\partial {\rm ln} Q}{\partial T} \right)_{N,V} + RT^2 \left(\frac{\partial^2 {\rm ln} Q}{\partial T^2} \right)_{N,V}
\end{equation}

\begin{equation}
S^o = R {\rm ln} Q + RT\left(\frac{\partial {\rm ln} Q}{\partial T} \right)_{N,V}
\end{equation}

where $R$ is the ideal gas constant and $T$ the absolute temperature. For an ideal gas, $Q(N,V,T)$ can be decomposed into \href{https://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)}{partition function} of different degrees of freedom: electronic (el), translational (tr), rotational (rot) and vibrational (vib) motions. Therefore, for a molecular ideal gas, $Q(N,V,T)$  can be expressed as:

\begin{equation}
Q(N,V,T) = \frac{\left(q_{\rm el} q_{\rm tr} q_{\rm rot} q_{\rm vib}\right)^N}{N!}
\label{Partition}
\end{equation}

The \href{https://en.wikipedia.org/wiki/Rigid_rotor}{rigid rotator} and the \href{https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator}{quantum harmonic oscillator} can be used to approximate the contribution of the rotational and vibrational motions. The partition function of a rigid rotator is defined as 

\begin{equation}
q_{\rm rot} = \frac{T}{\sigma \Theta_{\rm rot}}
\end{equation}

where $\sigma$ is the symmetry number, $\Theta_{\rm rot}$ the rotational temperature defined as 

\begin{equation}
\Theta_{\rm rot} = \frac{h^2}{8\pi^2 {\rm I} k_{\beta}}
\end{equation}

where ${\rm I}$ is the moment of inertia, $k_{\beta}$ the Boltzmann constant, and $h$ the Planck constant. The partition function of a quantum harmonic oscillator is defined as 

\begin{equation}
q_{\rm vib} = \frac{e^{\frac{-\beta h \nu}{2}}}{1-e^{-\beta h \nu}}
\end{equation}

where $\nu$ is the vibrational frequency of the oscillator. If we define the vibrational temperature as $\Theta_{\rm vib} = \frac{h \nu}{k_{\beta}}$, then we will have 

\begin{equation}
q_{\rm vib} = \frac{e^{\frac{\Theta_{\rm vib}}{2T}}}{1- e^{\frac{\Theta_{\rm vib}}{T}}}
\end{equation}

Applying Eqn.\ref{Partition} to Eqn.\ref{intE} followed by the multiplication rule in logarithm yields:

\begin{equation}
E = E_{tr} + E_{rot} + E_{vib}
\end{equation}

and similarly, 
\begin{equation}
C_v = C_{tr} + C_{rot} + C_{vib}
\end{equation}

\begin{equation}
S^o = S_{tr} + S_{rot} + S_{vib}
\end{equation}

Thermochemical properties are computed automatically by the \actcmd{nma} command. For more details, please see Van der Spoel \textit{et al} \cite{Spoel2018a}.

The \actcmd{nma} used above to generate the infrared spectrum in Fig.~\ref{fig:ethanol}  computes the thermochemical variables as well (Table~\ref{tab:thermo}).
Information on calculation of the enthalpy of formation will be pubslished in the near future.

\begin{table}[ht]
    \centering
    \caption{Thermochemical values for ethanol at 298.15 K based on the OPLS2020 force field~\cite{Jorgensen2023a}.}
    \label{tab:thermo}
    \begin{tabular}{lcc}
\hline
    Property & Experiment & OPLS2020 \\
\hline
$S^o$ (J/mol K) & 281~\cite{Handbook2022a} & 273.0  \\
$C_v$ (J/mol K) & 57~\cite{Ghahremanpour2016a} & 62.2 \\
\hline
\end{tabular}
\end{table}

=========================
Second Virial Coefficient
=========================
The second \href{https://en.wikipedia.org/wiki/Virial_coefficient}{virial coefficient} is the second term in the \href{https://en.wikipedia.org/wiki/Virial_expansion}{virial expansion} that describes the deviation from the ideal gas law for real gases:
\begin{equation}
    \displaystyle{\frac{P}{RT\rho}} = A + B_2(T)\rho + C_3(T)\rho^2 + ...
\end{equation}
with $P$ the pressure, $R$ the gas constant, $T$ the temperature and $\rho$ the density.

$B_2(T)$ is a useful property gauging interactions in the gas phase because experimental values are available for close to 2000 compounds as a function of temperature. It is computed from an integral weighting the interaction between two molecules over three-dimensional space.
\begin{equation}
    B_2^{cl}(T) = -\frac{1}{2}\int_0^{\infty} \left< e^{-\beta u_{12}({\mathbf r})} - 1\right> d{\mathbf r}
\end{equation}
where u$_{12}$(r) is the interaction energy between two compounds (\textit{See} \ref{Eint}), $\beta = 1/k_BT$ and
the integral is over all space and relative orientations of the compounds. If we sample these adequately (including at close, repulsive, distance) we can simplify the integral to a one-dimensional one:
\begin{equation}
  B_2^{cl}(T) = -2\pi\int_0^{\infty} r^2 \left< e^{-\beta u_{12}(r)} -1 \right> dr  
\end{equation}
The above equation is entirely classical and quantum corrections have to be added according to:
\begin{equation}
     B_2^F(T) = \frac{\hbar^2}{24(k_BT)^3} \sum_{j=1}^{2} \left[\frac{\left< {\mathbf F}_j^2 \right>}{m_j}\right]
\end{equation}
for the force on the compounds and
\begin{equation}
    B_2^{\tau}(T) = \frac{\hbar^2}{24(k_BT)^3} \sum_{j=1}^2\left[\sum_{\alpha=x,y,z} \frac{\left< \tau^2_{j,\alpha}\right>}{I_{j,\alpha}} \right]
\end{equation}
for the torque on the compounds,
where $m_j$ is the mass of the molecules $j$ and $\mathbf{F}^2$ is the averaged square force on one molecule given by
\begin{equation}
   \left<\mathbf{F}^2\right> = k_B T \int_0^{\infty} \left< e^{\displaystyle{- \beta u_{12}(\mathbf{r})}}\left[\nabla u_{12}(\mathbf{r})\right]^2\right> d \mathbf{r} 
\end{equation}
and where $I$ is the moment of inertia of the molecule and
 $\tau^2$ is the average square torque on one molecule defined by
\begin{equation}
     \left<\tau^2_{j,\alpha}\right> = k_B T \int_0^{\infty}  \left< e^{{- \beta u_{12}(\mathbf{r})}} \left[\nabla_{\omega} u_{12}(\mathbf{r})\right]_{j,\alpha}^2\right> d \mathbf{r} .
\end{equation}

The change in $B_2(T)$ as a function of temperature can be used to scrutinize the repulsive and attractive components of the potential energy. $B_2(T)$ is negative at low temperatures due to attraction forces, while it becomes positive at higher temperatures as repulsion forces start to dominate, and passes through a maximum and eventually decreases at very high temperatures where repulsion force are fully dominant~\cite{Amdur1958a}.  

Code to compute the second virial coefficient is available in the \actcmd{b2} command. Since the calculation is relatively expensive it has been implemented to use parallel processing using the message passing library. You can run it on a 16-core machine like
\begin{verbatim}
    mpirun -n 16 alexandria b2 -v 3 -g TIP4P -b2 TIP4P -ff TIP4P
    -f water#water.pdb -T1 373.15 -T2 673.15 -dT 25.0
    -maxdimer 32768
\end{verbatim}
where TIP4P corresponds to the well-known water model~\cite{Jorgensen1983a}, the $T_1$ and $T_2$ are the temperature limits (note gas-phase for water), $dT$ is the temperature interval and maxdimer determines how many relative orientations will be evaluated.
Due to the underlying algorithm for \href{https://en.wikipedia.org/wiki/Sobol_sequence}{quasi-random numbers} this should be a power of two. The result is plotted in Fig.~\ref{fig:b2}.

\begin{figure}[ht]
    \centering
    \includegraphics[width=0.9\linewidth]{images/TIP4P-b2.pdf}
    \caption{Sample second virial coefficient and components for water using the TIP4P model~\cite{Jorgensen1983a}.}
    \label{fig:b2}
\end{figure}

\clearpage
\begin{huge}
\vspace{8cm}
\centerline{Reference Manual}
\end{huge}
\addcontentsline{toc}{section}{Reference Manual}
