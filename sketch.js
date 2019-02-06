/* Ising simulation example of a temperature sweep on a 100x100 spin lattice */

// Display options:

const CANVAS_WIDTH  = 900;
const CANVAS_HEIGHT = 600;
const FRAME_RATE    = 20;

// Simulation parameters:

var isingSettings = {
	width:           100,
	height:          100,
	mcStepsPerFrame: 10000,
	averagingFrames: 20,
	randomize:       true,
	temperature:     4,
	magenticField:   0,
	loopMode:        'TSWEEP',
	loopIncrement:   -0.05,
	imageFile:       null,
	dataFile:        'simulationData'
	}
var ferromagnet;

// Create canvas and Ising-simulation object:

function setup() {
	frameRate(FRAME_RATE);
	createCanvas(CANVAS_WIDTH, CANVAS_HEIGHT);
	background(0);
	
	ferromagnet = new Ising(isingSettings);
}

// Draw loop:

function draw() {
	ferromagnet.simulationStep();
}
