*******************************
Force Field Training Algorithms
*******************************

*This chapter was written by Juli{\'a}n Marrades and is based in part on his M.Sc. thesis~\cite{Marrades2022a}.*

Within the \actcmd{train_ff} module of the Alexandria Chemistry Toolkit you can choose among three algorithms to optimize force field parameters:
\begin{itemize}
\item  Markov Chain Monte Carlo (\actflag{-optimizer MCMC})
\item Genetic Algorithm (\actflag{-optimizer GA})
\item  Hybrid GA/MCMC (\actflag{-optimizer HYBRID})
\end{itemize}
Let us see what they do behind the scenes and how to control them.

====
MCMC
====
Assume we want to set the values for five parameters (nParam = 5) by performing ten MCMC iterations (\actflag{-maxiter 10}). Then, the code "reads"

.. code-block:: c++

   for (int i = 0; i < 10; i++)
     for (int j = 0; j < 5; j++)
       stepMCMC();

That is, we do 10 X 5 = 50 MCMC steps. What is a MCMC step though? In essence,

* `prevDev`: previous deviation from data
* Choose a parameter and alter it
* `newDev`: new deviation from data
* If `newDev` < `prevDev`, we accept the parameter change and continue into the next step. Else, apply Metropolis Criterion. That is, with some probability we accept the change. Otherwise, we restore the parameter to its previous value and proceed into the next step.


Now, let us add some more detail to the algorithm.

#. we already have the deviation from data of the previous step in `prevDev`.
#. we arbitrarily choose a parameter out of the `param` vector, say `param[i]`, which is bounded by its maximum `pmax[i]` and its minimum `pmin[i]`, yielding the range `prange[i] = pmax[i] - pmin[i]`. Here is where the `-step [0, 1]` flag comes into play by specifying a fraction.\\
#. we draw a value of `delta`, which is uniformly distributed in `[- step * prange[i], + step * prange[i]]` and add it to the parameter value `param[i] += delta`. If the parameter has gone beyond its maximum, we set its value to the maximum. Same goes for the minimum.<br>
#. we compute the deviation `newDev` from data of the modified parameter vector. If `newDev < prevDev`, we accept the change and finish the MCMC step. Otherwise, we apply the Metropolis Criterion and finish the step. Then go back to 1.

--------------------
Metropolis Criterion
--------------------
If Mathematics is your thing and you {\em really} want to know the nitty-gritty stuff, you have thorough explanations of this method on \href{https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm}{Wikipedia} and Shuyi Qin's Master thesis~\cite{Qin2021a}.
Here, we shall say that the Metropolis Criterion allows us to take steps that do not get us closer to a minimum, giving the opportunity to explore the parameter space and avoid local minima. How exactly?

The probability of accepting a "bad" parameter change is controlled by the \actflag{-temp T} flag and follows the equation

.. math:: prob = exp(-deltaDev/T)  

where deltaDev = newDev - prevDev. Note that the unit of temperature is the same as that of the deviation.

Given deltaDev, a higher temperature gives a higher probability  of acceptance, as `-deltaDev/T` tends to `0` and `exp(0) = 1`. On the other hand, a lower `T` gives less chances of acceptance, as `-deltaDev/T` tends to $-\infty$ and $exp(-\infty)$ = 0.

ACT gives us the possibility to lower the temperature (simulated annealing) during the MCMC optimization with the \actflag{-anneal [0, 1]} flag. If we use \actflag{-anneal 0.5}, the temperature will be \actflag{-temp T} until 50\% of the iterations have been completed. Then, it will be linearly decreased until it reaches `0` on the last iteration.

There are two remarks to be made here:

*  The temperature is kept constant during the `nParam` MCMC steps that take place for a given iteration.
* Since division by `0` is not defined, we set `T = 1e-6` on the last iteration.

Fig.~\ref{fig:annealing} shows a schematic example of the temperature over time when we use actflag{\-maxiter 10 -temp 5 -anneal 0.5}.

.. figure:: ../images/annealing.png
   :name: fig-annealing
   :width: 60%
   :align: center

   Annealing during a MCMC run.

------------------
Multiple MCMC runs
------------------
We can optimize several candidate solutions in parallel using the \actflag{-pop_size} flag. If \actflag{-random_init} is used (default), each candidate solution will be initialized arbitrarily and respecting parameter ranges. If \actflag{-norandom_init} is employed, the candidate solution(s) will be initialized as specified by the force field file(s) provided by the user.

=================
Genetic Algorithm
=================
Quoting `WikipediaGA`_

  | In computer science and operations research, a genetic algorithm is a metaheuristic 
  | inspired by the process of natural selection that belongs to the larger class of
  | evolutionary algorithms. Genetic algorithms are commonly used to generate 
  | high-quality solutions to optimization and search problems by relying on biologically inspired operators such as mutation, crossover and selection.

.. _WikipediaGA: https://en.wikipedia.org/wiki/Genetic_algorithm

In this section, we describe our implementation of a genetic algorithm for force field parameterization. If this is the first time you hear about genetic algorithms and want to acquaint yourself, we recommend you to read Chapter 2 of Steven F. van Dijk's PhD thesis~\cite{VanDijk2001a} and/or head over to `Youtube`_.

.. _Youtube: https://www.youtube.com/

Our implementation is based upon a class definition and an external function for computing random numbers (Listing :ref:`listing-genome` ).

.. literalinclude:: ../code/genome.c++
   :linenos:
   :language: c++
   :name: listing-genome
   :caption: Class definition.

then the genome evolution 
boils down to the code in Listing :ref:`listing-evolve` 

.. literalinclude:: ../code/evolve.c++
   :linenos:
   :language: c++
   :name: listing-evolve
   :caption: Schematic of the evolution algorithm


\actflag{popSize} and \actflag{nElites} are assumed to be even.
Let us explore the different stages of the process.

--------------
Initialization
--------------
This stage fills the `params` field in the `Genome` class, generating popSize (\actflag{-pop_size}) genomes.

If we are using \actflag{-norandom_init}, the genomes will be initialized as specified by the force field file(s) provided. If we are employing \actflag{-random_init}, each genome will be initialized by setting a random value for each parameter, uniformly distributed over the allowed range.

-------------------
Deviation from data
-------------------
Here we fill the `deviation` field in for each `Genome` in the population by computing the deviation from the dataset.

-------
Sorting
-------
Sorting is not a mandatory step but may be required depending on the GA components selected by the user.
We sort the population in ascending order of `deviation`. Whether we sort or not is controlled by the \actflag{-sort} flag.

---------
Penalties
---------
At this stage, we may alter the population if certain conditions are met, with the main goal of preventing premature convergence and enforcing solution diversity.

  | Covering such a small portion of the space you are. Broaden your search, you should. - Yoda

To that end, we have a function `penalize()` which returns `true` if the population was penalized and `false` otherwise.

For now, there are two components in this function:

* **Volume.** This option enables \actflag{-sort}. If the volume spanned by the population divided by the total volume of the parameter space is smaller than \actflag{-vfp_vol_frac_limit [0, 1]|}, the {\em worst} fraction \actflag{-vfp_pop_frac [0, 1]} of genomes in the population will be randomly reinitialized. If \actflag{-log_volume} is used, volumes will be computed in logarithmic scale. * But wait, then the volume could be negative! Yes, we have to fix that!*
    
  | Death comes equally to us all, and makes us all equal when it comes. - John Donne
  
* **Catastrophe.** Each \actflag{-cp_gen_interval} generations, a fraction \actflag{-cp_pop_frac [0, 1]} of the genomes in the population will be randomly reinitialized. Genomes to reinitialize are arbitrarily chosen.


-----------------------
Selection probabilities
-----------------------
We provide three options for computing selection probabilities:

#. Rank (\actflag{-prob_computer RANK})
#. Fitness (\actflag{-prob_computer FITNESS})
#. Boltzmann (\actflag{-prob_computer BOLTZMANN})

The sum of the selection probabilities of the genomes in the population, is 1.

----
Rank
----
This option enables \actflag{-sort}.
The selection probability depends exclusively on the index (rank) of genome in the population (Listing :ref:`listing-rank`).

.. literalinclude:: ../code/rank.c++
   :language: c++
   :linenos:
   :name: listing-rank
   :caption: Calculation of the probability from the order of probabilities.

That is, the lower the index, the higher the probability of being selected.
The independence of the `deviation` avoids the possible phenomena of a genome with a very high selection probability dominating the population.

-------
Fitness
-------
The selection probability is inversely proportional to the `deviation` (Listing :ref:`listing-fitness`).

.. literalinclude:: ../code/fitness.c++
   :language: c++
   :linenos:
   :name: listing-fitness
   :caption: Calculation of the probability from the deviations from data.

---------
Boltzmann
---------
The `temperature` parameter is specified by the \actflag{-boltz_temp} flag and controls the smoothing of the selection probabilities (Listing :ref:`listing-boltz`). A higher value will avoid polarization in the probability values and vice versa. The \actflag{-boltz_anneal} flag allows us to decrease the temperature over time and has the same logic as \actflag{-anneal}, except that it targets the Boltzmann selection temperature and operates based on the maximum amount of generations.

.. literalinclude:: ../code/boltz.c++
   :language: c++
   :linenos:
   :name: listing-boltz
   :caption: Use of Boltzmann-weighting when calculating the probability

-------
Elitism
-------
In order to avoid losing the best candidate solutions found so far, the GA will move the top nElites \actflag{-n_elites} genomes, {\em unchanged}, into the new population. That means, the genome will not undergo crossover nor mutation.

When \actflag{-n_elites $>$ 0}, \actflag{-sort} will be enabled.

---------
Selection
---------
Once the selection probabilities are computed, the population becomes a \href{https://en.wikipedia.org/wiki/Probability_density_function}{probability density function} from which we can sample genomes based on their `probability`.

As of now, only a vanilla selector is available. It only looks at the probability and can select the same genome to be `parent1` and `parent2`.

---------
Crossover
---------
With certain probability `prCross`, controlled by the \actflag{-pr_cross} flag, two parents will combine their parameters to form two children.
Right now, only an $N$-point crossover is available, where $N$ is defined by the \actflag{-n_crossovers} flag.

If `N = 1`, we arbitrarily select on crossover point (`v`) and::

      v                        v
  |X|X|X|X|X|X|X|X|X|  --> |X|X|Y|Y|Y|Y|Y|Y|Y|
  |Y|Y|Y|Y|Y|Y|Y|Y|Y|  --> |Y|Y|X|X|X|X|X|X|X|

If `N = 2`, we arbitrarily select two crossover points (`v`) and::

    v     v                  v     v
  |X|X|X|X|X|X|X|X|X|  --> |X|Y|Y|Y|X|X|X|X|X|
  |Y|Y|Y|Y|Y|Y|Y|Y|Y|  --> |Y|X|X|X|Y|Y|Y|Y|Y|

If `N = 3`, we arbitrarily select three crossover points (`v`) and::

         v     v     v            v     v     v
  |X|X|X|X|X|X|X|X|X|  --> |X|X|Y|Y|Y|X|X|X|Y|
  |Y|Y|Y|Y|Y|Y|Y|Y|Y|  --> |Y|Y|X|X|X|Y|Y|Y|X| 



If `N = 4`... you get the idea (hopefully).

--------
Mutation
--------
Given the mutation probability `prMut`, controlled by \actflag{-pr_mut}, we iterate through the parameters. If the probability is met, we alter the parameter, otherwise, we leave it unchanged.

.. code-block:: c++
   
   for (int i = 0; i < nParams; i++)
   {
     if (random() <= prMut)
     {
       changeParam(params, i);
     }
   }

The parameter is changed in the same way as in MCMC, by a fraction of its allowed range and not allowing values outside of it. The fraction in this case is controlled by the `-percentage` flag, which has the same meaning as `-step`, but applies to GA instead of MCMC.

-----------
Termination
-----------
This stage decides whether the GA evolution should continue or halt. We allow the user to tweak the termination criteria with several flags:

#. \actflag{-max_generations}: evolution will halt after so many generations.
#. \actflag{-max_test_generations} (disabled by default): evolution will halt if in the last \actflag{-max_test_generations} the best `deviation` of the test set found so far has not improved.

======
HYBRID
======
Even though this optimizer has a very fancy name, it is just a GA with MCMC as its mutator engine.
When MCMC acts as a mutator, it will always alter the genomes independently of \actflag{-pr_mut}. Also, it is important to note that the simulated annealing by default is applied independently in each MCMC run.
For instance, in case we would use \actflag{-max_generations 2 -maxiter 10 -temp 5 -anneal 0.5}, Fig.~\ref{fig:anneal-hybrid} shows the temperature during the MCMC part.

.. figure:: ../images/annealing_hybrid.png
   :name: fig-anneal-hybrid
   :align: center

   Annealing in the hybrid algorithm

However, when using the flag \actflag{-anneal_globally}, the starting temperature of the annealing will be decreased in steps at the beginning of each generation.
