# Ising
This is a P5.js animation of a simple [Ising model](https://en.wikipedia.org/wiki/Ising_model) for a ferromagnet on a two-dimensional square lattice. This program visualizes magnetic domains and can simulate temperature or magnetic field curves. 
If wanted, the simulation data and indiviual frames are automatically saved for later evaluation or presentation. 
(In particular, the captured frames can be used to create a video animation.)

## Physical Model
The state of a Ising ferromagnet is described by a discrete spin vector S 

	S = (s_1, s_2, ..., s_N) ∈ {-1, +1}^N
	
where each spin s_i can assume only two possible values -1 or +1 (spin-down or spin-up). The energy of a certain state is given by the dimensionless Hamilitonian (J = mu = 1):

	H(S) = - ∑_<i,j> s_i s_j - h ∑_i s_i = - ∑_i (h_i + h) s_i
	
where the first sum is taken over all pairs i,j of adjacent spins. h denotes the external magnetic field and h_i is the so-called exchange field. In thermodynamic equilibrium, the systems shows a Boltzmann distribution (k = 1):

	p(S) = const * exp( - H(S) / T) 

where T is the temperature of the spin system. 

## Metropolis Algorithm
The Ising system is simulated by a Markov chain of spin states S(t) as a function of discrete time t, which fulfills the Boltzmann distributions. The algorithm can summarized by the following:

1. Choose an initial spin configuration S(0)
2. Choose a random i and calculate the energy change for a single spin-flip s_i → -s_i:
   ```
	∆E = - 2 s_i (h_i + h)
   ```	
3. Attempt a spin-flip with the transition probability
   ```
	                ┌
	                │ 1          	∆E ≥ 0
	p(s_i → -s_i) = ┤
	                │ exp(∆E / T)	∆E < 0
	                └
   ```
   For this, generate a random number r ∈ [0,1] and flip the spin, if r < p(s_i → -s_i). Otherwise keep the current state.
4. Go to 2. and repeat

## Usage

Ising is called with a settings object as a single argument, that summarizes all parameters of the simulation:
```javascript
var mySettings = { 
	width:            Number                                  // width of the spin lattice 
	height:           Number,                                 // height of the spin lattice
	mcStepsPerFrame:  Number,                                 // number of iterations in a simulation step
	averagingFrames:  Number,                                 // number of simulation steps for averaging 
	temperature:      Number,                                 // initial temperature 
	magneticField:    Number,                                 // initial magnetic field 
	randomizeLattice: Boolean,                                // initial spin configuration 
	loopMode:         'TSWEEP'|'TLOOP'|'HSWEEP'|'HLOOP'|null, // specifies temperature/field curve
	loopTargetValue:  Number,                                 // target value for temperature/field 
	loopIncrement:    Number,                                 // step size for temperature/field 
	imageFile:        String|null,                            // file name for animation frames
	dataFile:         String|null                             // file name for the simulation data
};
var mySimulation = new Ising(mySettings);
```
Only the parameters imageFile, dataFile and loopMode are optional and should be set to null if not used. In order to excute and animate the simulation, the function simulationStep() must be called in the draw function:
```javascript
function draw() {
	mySimulation.simulationStep();
}
```
Image and data files are automatically saved to the local download directory.
