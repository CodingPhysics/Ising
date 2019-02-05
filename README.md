# Ising
This is a P5.js animation of a simple [Ising model](https://en.wikipedia.org/wiki/Ising_model) for a ferromagnet on a two-dimensional square lattice.

This program visualizes magnetic domains and can simulate temperature or magnetic field curves. In the file 'sketch.js', the settings of the simulation are specified (The example setting showcases a cooling curve). If filenames are given, the simulation data and indiviual frames are automatically saved for later evaluation or presentation. In particular, the captured frames can be used to create a video animation.

## Physical Model
The state of spin system with the idividual spins s_i is described by a discrete spin vector:

    S = (s_1, s_2, ..., s_N),	s_i = -1, +1
	
The energy of this system is given in the Ising model by the Hamiltonian:

    H(S) = - J Σ_<i,j> s_i s_j - μ h Σ_i s_i

where J is the exchange coupling, μ is the magnetic moment and h is the magnetic field. 
The first sum is taken over all pairs i,j of adjacent spins.

In thermal equilibrium the probability of a state S is given by a Boltzmann distribution:

    p(S) = const * exp(- H(S)/k T)

where k is the Boltzmann constant and T is the absolute temperature.

## Metropolis Algorithm

The simulation the applies the [Metropolis algorithm](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm).
Thereby, the system is simulated by a Markov chain of different spin configurations as a function of time t:

    S(0), S(1), ... , S(t), ...

The rate at which a state appears in the Markov chain should correspond to its probability in the Boltzmann distribution, so that the ensemble mean can be replaced by a time average. Successive elements of this Markov chain differ only by a single spin flip:

	S  = (s_1, ... , s_i , ... , s_N),
	S' = (s_1, ... , -s_i, ... , s_N)
	
In order to approximate a Boltzmann distribution, the transition S -> S' is assigned to a probability
 
	p(S -> S') = 1                        	if H(S') < H(S)
	             exp((H(S) - H(S'))/k T)    elsewhere

that fulfills the detailed balance p(S) p(S -> S') = p(S') p(S' -> S) in thermal equilibrium.

Thus, the algorithm is summarized by:

1. Choose an initial state S(0)
2. Choose a random i and calculate the probability p(S -> S')
3. Perform a spin flip S(t + 1) = S'(t) (with probability p(S -> S')) or keep the original state S(t + 1) = S(t). 
4. Repeat 2. and 3.
				 
## References

[1] Wolfgang Kinzel, Georg Reents, "Physik per Computer", Spectrum, 1996
