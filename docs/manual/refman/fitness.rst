****************
Fitness function
****************
In genetic algorithm terminology, the idea is to maximize the fitness of the population. However, in practice the ACT minimizes a loss function :math:`\chi^2` where zero indicates perfect correspondence to reference data and algorithmic bounds.
The loss function consists of multiple components that can be combined at will by the user. According to the ACT paper, we have

.. math:: F(\Theta) = \Omega \mathbf{X}^p(\Theta) + \Lambda(\Theta)

where :math:`\Theta` is the force field genome, containing all the parameters and :math:`\mathbf{X}^n` is a vector representation of residuals. In case that p equals 2, the loss function corresponds to a least squares form, if p equals 1 the loss function is the the mean absolute error. 
Finally, :math:`\Lambda` adds a penalty to keep (some of the) parameters within bounds. 
The :math:`\chi^2` that is printed in the ACT log file corresponds to :math:`F(\Theta)`.
The first term can be expanded into a sum, since multiple observable can be trained at the same time.

.. math::  F(\Theta) = \sum_{i=1}^{N_{obs}} \Omega_i \sum_{j=1}^{N_i}\left|x_{i,j}^{calc}(\Theta)-x_{i,j}^{ref}\right|^p + \Lambda(\Theta)

where :math:`N_{obs}` is the number of observables, the :math:`\Omega_i` are the corresponding user-supplied weights for each of these and :math:`N_{i}` is the number of individual data points :math:`x` for observable i.

--------------------------
Supported training targets
--------------------------

Table :ref:`tab-fitness` lists the supported observables and the corresponding command line flags.

.. table:: Overview of command-line flags and corresponding observables that can be used in training a force field in the ACT.
   :name: tab-fitness

   +--------------+------------------------------------------------------------+
   | **Flag**     | **Observable**                                             |
   +==============+============================================================+
   | -fc_epot     | Potential energy of a compound.                            |
   +--------------+------------------------------------------------------------+
   | -fc_force    | Force components on each atom in a compound.               |
   +--------------+------------------------------------------------------------+
   | -fc_freq     | Vibrational frequencies from a normal mode analysis (NMA). |
   +--------------+------------------------------------------------------------+
   | -fc_inten    | Intensities corresponding to the frequencies from NMA.     |
   +--------------+------------------------------------------------------------+
   | -fc_mu       | Dipole components of a compound.                           |
   +--------------+------------------------------------------------------------+
   | -fc_quad     | Quadrupole components of a compound.                       |
   +--------------+------------------------------------------------------------+
   | -fc_oct      | Octopule components.                                       |
   +--------------+------------------------------------------------------------+
   | -fc_hexadec  | Hexadecapole components.                                   |
   +--------------+------------------------------------------------------------+
   | -fc_esp      | Electrostatic potential (ESP) on all of the grid points.   |
   |              | This can be modulated by the -watoms flag which introduces |
   |              | a weight for the ESP on the atoms.                         |
   +--------------+------------------------------------------------------------+
   | -fc_polar    | Six independent components of the molecular polarizability |
   |              | tensor.                                                    |
   +--------------+------------------------------------------------------------+
   | -fc_inter    | The total interaction energy from SAPT calculations.       |
   +--------------+------------------------------------------------------------+
   | -fc_elec     | Electrostatic component of the interaction energy.         |
   +--------------+------------------------------------------------------------+
   | -fc_induc    | Induction  component of the interaction energy.            |
   +--------------+------------------------------------------------------------+
   | -fc_allelec  | Sum of the electrostatic and induction components.         |
   +--------------+------------------------------------------------------------+
   | -fc_disp     | Dispersion component of the interaction energy.            |
   +--------------+------------------------------------------------------------+
   | -fc_exch     | Exchange component of the interaction energy.              |
   +--------------+------------------------------------------------------------+
   | -fc_exchind  | Sum of exchange and induction components.                  |
   +--------------+------------------------------------------------------------+
   | -fc_deltahf  | The :math:`\Delta` HF part of the induction can be trained |
   |              | independently of the induction component. If this is used, |
   |              | this part should be subtracted from the induction to       |
   |              | prevent double counting.                                   |
   +--------------+------------------------------------------------------------+
   | -fc_bound    | Penalty for going outside bounds (see below).              |
   +--------------+------------------------------------------------------------+

--------------------
Parameter mutability
--------------------
A force field consists many parameters and typically only a subset is trained at any one time. Therefore, parameters are qualified by a mutability value as indicated in Table :ref:`tab-mutability`.

.. table:: Parameter mutabilities and their interpretation in ACT.
   :name: tab-mutability

   +----------------+------------------------------------------------------------+
   | **Mutability** | **Meaning**                                                |
   +================+============================================================+
   | Fixed          | This parameter will not be modified in training.           |
   +----------------+------------------------------------------------------------+
   | Bounded        | Parameter values will be strictly between the minimum      |
   |                | and maximum value specified in the force field. If the     |
   |                | initial value is outside the bounds, it will be reset to   |
   |                | the nearest bound (either minimum or maximum). If the      |
   |                | minimum is equal to the maximum the mutability will be set |
   |                | to Fixed.                                                  |
   +----------------+------------------------------------------------------------+
   | Dependent      | Parameter value is not trained directly, but may change    |
   |                | due to changes in other parameters. An example would be    |
   |                | a pair parameter such as :math:`\sigma_{ij}` that is       |
   |                | derived from :math:`\sigma_{j}` and :math:`\sigma_{j}`     |
   |                | using combination rules.                                   |
   +----------------+------------------------------------------------------------+
   | ACM            | Alexandria Charge Model is used for charges that are       |
   |                | generated using the EEM or SQE algorithms                  |
   |                | (section :ref:`sec-charges`). If different minimum and     |
   |                | maximum values are given in the force field and the        |
   |                | *-fc_bound* flag is used, a bounds penalty will be applied.|
   +----------------+------------------------------------------------------------+
   | Free           | Parameter can adopt any value.                             |
   +----------------+------------------------------------------------------------+
   


--------------
Bounds Penalty
--------------
A bounds penalty can be applied to charges when the mutability is ACM (Table :ref:`tab-mutability`). In that case the penalty V(q) is

.. math:: V(q) = \begin{split}
   & (q - q_{min})^2 & \hskip2em\mathrm{if}\hskip1em q < q_{min}\\
   & (q - q_{max})^2 & \hskip2em\mathrm{if}\hskip1em q > q_{max}\\
   & 0               & \hskip2em\mathrm{else}
   \end{split}

Importantly, if a point+distributed model is used :cite:p:`Spoel2025a`, then the sum of these charges is used as the target q. So if your force field specifies::

  <particletype identifier="hw" type="Atom" description="hydrogen">
    <parameter type="charge" unit="e" value="0.5" uncertainty="0" minimum="0.3" maximum="0.7" ntrain="0" mutability="ACM" nonnegative="no"/>
  </particletype>

and the hw atom type is equiped with a distributed charge v1hw at the same position, it is the sum of these charges that will be penalized if it is outside the interval 0.3 - 0.7.
