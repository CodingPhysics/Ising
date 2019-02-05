/* Ising simulation example of a temperature sweep on a 100x100 spin lattice */

// Display options:

const CANVAS_WIDTH  = 900;
const CANVAS_HEIGHT = 600;
const FRAME_RATE    = 20;

// Simulation parameters:

const LATTICE_WIDTH                = 100;
const LATTICE_HEIGHT               = 100;
const MARKOV_CHAIN_STEPS_PER_FRAME = 10000;
const FRAMES_FOR_AVERAGING         = 20;
const RANDOMIZE                    = true;
const TEMPERATURE                  = 4;
const MAGNETIC_FIELD               = 0;
const STEP                         = -0.05;
const FINAL_VALUE                  = 0.05;
const LOOP_MODE                    = 'TSWEEP';
const IMAGE_FILE                   = null;
const DATA_FILE                    = 'simulationData';

var ferromagnet;

// Create canvas and Ising-simulation object:

function setup() {
	frameRate(FRAME_RATE);
	createCanvas(CANVAS_WIDTH, CANVAS_HEIGHT);
	background(0);
	
	ferromagnet = new Ising( LATTICE_WIDTH, LATTICE_HEIGHT,                        // size of the lattice
							 MARKOV_CHAIN_STEPS_PER_FRAME, FRAMES_FOR_AVERAGING,   // iteration steps
							 RANDOMIZE, TEMPERATURE, MAGNETIC_FIELD,               // lattice initialization
							 STEP, FINAL_VALUE, LOOP_MODE,                         // loop specifications
							 IMAGE_FILE, DATA_FILE                              ); // image and data export
}

// Draw loop:

function draw() {
	ferromagnet.simulationStep();
}