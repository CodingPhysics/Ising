/* MONTE CARLO SIMULATION ON AN ISING SPIN SYSTEM ON A TWO-DIMENSIONAL PERIODIC LATTICE */

/*************************************************************************************
From: 	Coding Phyisics
email:	CodingPhysicSimulation@gmail.com
github:	https://github.com/CodingPhysics
date:   2019-02-05
*************************************************************************************/

/*************************************************************************************
PHYSICAL MODEL:

The state of spin system is described by a discrete spin vector:

	S = (s_1, s_2, ..., s_N),	s_i = -1, +1 
	
The energy of this system is given in the Ising model by the Hamiltonian:

	H(S) = - J ∑_{<i,j>} s_i s_j - μ h ∑_i s_i

where J is the exchange coupling, μ is the magnetic moment and h is the magnetic field. 
The first sum is taken over all pairs i,j of adjacent spins.

In thermal equilibrium the probability of a state S is given by a Boltzmann distribution:

	p(S) = const * exp(- H(S)/k T)

where k is the Boltzmann constant and T is the absolute temperature.

ALGORITHM:

The system is simulated by a Markov chain of different spin configurations:

	S(0), S(1), ... , S(t), ...

The ensemble mean of an obserable A is replaced by the average over a time period:

    <A> = ∑_{t = t0}^{t0 + T} A(S(t)) / T

Successive elements of this Markov chain differ only by a spin flip:

	S  = (s_1, ... , s_i , ... , s_N),	 
    S' = (s_1, ... , -s_i, ... , s_N)
	
The probability for a transition S -> S' in the Markov chain is given by 
 
	p(S -> S') = 1                        	if H(S') < H(S)
	             exp((H(S) - H(S'))/k T)    elsewhere

The algorithm can summarized by:

1. Choose an initial state S(0)
2. Choose a random i and calculate the probability p(S -> S')
3. Perform a spin flip S(t + 1) = S'(t) (with probability p(S -> S')) 
   or keep the original state S(t + 1) = S(t) (with probability 1 - p(S -> S')) 
4. Repeat 2. and 3.
				 
REFERENCES:

[1] Wolfgang Kinzel, Georg Reents, "Physik per Computer", Spectrum, 1996

***************************************************************************************/

function Ising(nx, ny, mcSteps, averagingFrames, randomize, temp, hField, step, finalValue, loopMode, imageFile, dataFile){

	/***************
	** VARIABLES: **
	***************/

	this.nx              = nx;              // width of the spin-lattice
	this.ny              = ny;              // height of the spin-lattice
	this.mcSteps         = mcSteps;         // number of markov chains steps in an iteration
	this.averagingFrames = averagingFrames; // number of iteration for averaging
	this.temp            = temp;            // intial temperature (in units of J/k) 
	this.hField          = hField;          // intial magnetic field (in units of J/mu)
	this.step            = step;            // temperature/field step in the loop
	this.finalValue      = finalValue;      // final temperature/field value in the loop
	this.loopMode        = loopMode;        // specifies the loop type ('TSWEEP','HSWEEP','TLOOP', 'HLOOP')
	this.imageFile       = imageFile;       // filename for simulation images
	this.dataFile        = dataFile;        // filename for simulation data
	
	console.log(	'----- SIMULATION PARAMETERS -----' 					+ '\n' +
					'  Lattice size:            ' + this.nx + 'x' + this.ny + '\n' +
					'  MC steps per interation: ' + this.mcSteps 			+ '\n' +
					'  Interations for average: ' + this.averagingFrames 	+ '\n' +
					'  Initial temperature:     ' + this.temp 				+ '\n' +
					'  Initial magnetic field:  ' + this.hField             + '\n' +
					'  Random initialization:   ' + randomize               + '\n' +
					'  Loop type:               ' + this.loopMode           + '\n' +
					'  Final value:             ' + this.finalValue         + '\n' +
					'  Iteration step:          ' + this.step               + '\n' +
					'---------------------------------'									);
	
	this.iter          = 1;                    // iteration variable
	this.timeStamp     = new Date().getTime(); // starting time
	this.dataTable     = new p5.Table();       // table for simulation data
	this.dataTable.addColumn('temperature');
	this.dataTable.addColumn('field');
	this.dataTable.addColumn('magnetization');
	this.dataTable.addColumn('order parameter');
	
	this.backward      = false;                                      // sweep direction in a loop
	this.initialValue  = (loopMode.charAt(0) == 'T')? temp : hField; // initial value in a loop 

	this.mag           = 0;     // magnetization of the spin-lattice 
	this.order         = 0;     // local order parameter
	this.averageMag    = 0;     // average magnetization
	this.averageOrder  = 0;     // average order parameter
	this.spin          = [];    // spin matrix (spin[i][j] = -1, +1)
	this.exchangeField = [];    // exchange field (sum of adjacent spins)
	this.spinBuffer    = [];    // buffer for faster image drawing
	
	/****************************************
	** INITIALIZATION OF THE SPIN LATTICE: **
	****************************************/
	
	for(var i = 0; i < this.nx; i++){
		this.spin.push([]);
		this.spinBuffer.push([]);
		for(var j = 0; j < this.ny; j++){
			if(randomize){
				this.spin[i].push(2*Math.floor(2*Math.random()) - 1);	// random value (-1 or +1)
			} else {
				this.spin[i].push(1);									// high spin state
			}
			this.spinBuffer[i].push(0);
		}
	}
	for(var i = 0; i < this.nx; i++){
		this.exchangeField.push([]);
		for(var j = 0; j < this.ny; j++){
			// spin sum over all next neighbors gives exchange Field in units of J:
			this.exchangeField[i].push(	this.spin[(i + 1)%this.nx][j] + 
										this.spin[(i - 1 + this.nx)%this.nx][j] + 
										this.spin[i][(j + 1)%this.ny] + 
										this.spin[i][(j - 1 + this.ny)%this.ny]		);			
		}
	}
	
	/**********************************
	** SIMULATION STEP IN DRAW LOOP: **
	**********************************/
	
	this.simulationStep = function(){
		console.log('ITERATION ' + this.iter);
		
		this.metropolis();	// execute monte carlo steps
		this.update();		// update physical quantities
		this.show();		// display spin lattice
		this.iterate();		// increment temperature or field and repeat while unfinished
	}
	
	/**************************
	** METROPOLIS ALGORITHM: **
	**************************/
	
	this.metropolis = function(){
		console.log('-> Metropolis algorithm');
		
		for(var k = 0; k < this.mcSteps; k++){
			var i = Math.floor(this.nx*Math.random());
			var j = Math.floor(this.ny*Math.random());
			var deltaE = - 2*this.spin[i][j]*(this.exchangeField[i][j] + this.hField);
			
			if(deltaE >= 0 || Math.random() < Math.exp(deltaE/this.temp)){
				this.exchangeField[(i + 1)%this.nx][j]           -= 2*this.spin[i][j];
				this.exchangeField[(i - 1 + this.nx)%this.nx][j] -= 2*this.spin[i][j];
				this.exchangeField[i][(j + 1)%this.ny]           -= 2*this.spin[i][j];
				this.exchangeField[i][(j - 1 + this.ny)%this.ny] -= 2*this.spin[i][j];
				this.spin[i][j] = -this.spin[i][j];
			}
		}
	}
		
	/********************************
    ** UPDATE PHYSICAL QUANTITIES: **
    ********************************/
	
	this.update = function(){
		console.log('-> Update physical quantities');
		
		// calculate magnetization:
		
		var weight = 1./this.nx/this.ny;
		this.mag   = 0;
		
		for(var i = 0; i < this.nx; i++){
			for(var j = 0; j < this.ny; j++){
				this.mag += this.spin[i][j]*weight;
			}
		}
		
		this.averageMag   += this.mag/this.averagingFrames;
		
		// calculate local order parameter:
		
		weight    /= 4;
		this.order = 0;
		
		for(var i = 0; i < this.nx; i++){
			for(var j = 0; j < this.ny; j++){
				this.order += this.exchangeField[i][j]*this.spin[i][j]*weight;
			}
		}
		
		this.averageOrder += this.order/this.averagingFrames;
		
		// update datatable if averaging is completed:
		
		if(this.iter % this.averagingFrames == 0){	
			var newRow = this.dataTable.addRow();
			newRow.setString('temperature'    , this.temp.toExponential(6)   		); 
			newRow.setString('field'          , this.hField.toExponential(6) 		); 
			newRow.setString('magnetization'  , this.averageMag.toExponential(6)    ); 
			newRow.setString('order parameter', this.averageOrder.toExponential(6)  );
			
			console.log(	'------ UPDATED DATATABLE ------' 		  + '\n' +
							'  Temperature:     ' + this.temp         + '\n' +
							'  Magnetic field:  ' + this.hField 	  + '\n' +
							'  Magnetization:   ' + this.averageMag   + '\n' +
							'  Order parameter: ' + this.averageOrder + '\n' +
							'-------------------------------'					);
			this.averageMag   = 0;
			this.averageOrder = 0; 
		}
	}
		
	/**************************
	** DISPLAY SPIN LATTICE: **
	**************************/
		
	this.show = function(){
		console.log('-> Display spin lattice');
		
		// draw spin lattice:
		
		var boxSize = height/this.ny;		  
		var margin = (width - boxSize*this.nx)/2;    
		
		noStroke();
		for(var i = 0; i < this.nx; i++){
			for(var j = 0; j < this.ny; j++){
				if(this.spin[i][j] !== this.spinBuffer[i][j]){
					fill( (this.spin[i][j] < 0) ? 30 : 225 );
					rect(i*boxSize + margin, j*boxSize, boxSize, boxSize);
					this.spinBuffer[i][j] = this.spin[i][j];
				}
			}
		}
		
		// draw scales for temperature, magnetic field, magnetization and order parameter:
		
		fill(0);
		rect(0,0,margin,height);
		rect(width-margin,0,margin,height);
		
		this.drawScale(this.temp,    0, 4, 'T',   0.3*margin,         Math.floor(margin/10), 5);
		this.drawScale(this.hField, -1, 1, 'H',   0.7*margin,         Math.floor(margin/10), 5);
		this.drawScale(this.mag,    -1, 1, 'M',   width - 0.7*margin, Math.floor(margin/10), 5);
		this.drawScale(this.order,  -1, 1, 'eta', width - 0.3*margin, Math.floor(margin/10), 5);
		
		// save canvas as png image:
		
		if(this.imageFile != null){
			saveCanvas(this.imageFile + this.iter + '.png');
			console.log('   Canvas saved as ' + this.imageFile + this.iter + '.png');
		}
	}
	
	this.drawScale = function(u, u1, u2, label, posx, textHeight, tics){
		var y2      = 2*textHeight;
		var y1      = height-2*textHeight;
		var ticSize = 0.5*textHeight;
		var posy    = y1 + (y2 - y1)*(u - u1)/(u2 -  u1);
		
		noStroke();
		fill(255);
		triangle(posx, posy, posx + 2*ticSize, posy + ticSize, posx + 2*ticSize, posy - ticSize);
		
		textSize(textHeight);
		textAlign(CENTER,TOP);
		text(label, posx, 0);
		textAlign(CENTER,BOTTOM);
		text(u.toFixed(3), posx, height);
		textAlign(RIGHT,CENTER);
		
		for(var tic = 0; tic < tics; tic++){
			var y     = y1 + (y2 - y1)*tic/(tics - 1);
			var value = u1 + (u2 - u1)*tic/(tics - 1);
			
			strokeWeight(0);
			text(value, posx - 1.5*ticSize, y);
			strokeWeight(2);
			stroke(255);
			line(posx, y, posx - ticSize, y);
		}
		line(posx, y2, posx, y1); 
	}
	
	/*******************
	** ITERATION STEP **
	*******************/
	
	this.iterate = function(){
		
		if(this.iter % this.averagingFrames == 0){	// if averaging is completed
			if(this.exitCondition()){
				this.terminate();					// end simulation ..
			} else {								
				this.increment();					// ... or change temperature/field
			}
		}
		this.iter++;
	}
	
	this.increment = function(){
		switch(this.loopMode){
			case 'TLOOP':
				if(!this.backward && Math.abs(this.temp - this.finalValue) < 0.5*Math.abs(this.step)){
					this.backward = true;
					this.step = -this.step;
				}
			case 'TSWEEP':
				this.temp += this.step;
				break;
			case 'HLOOP':
				if(!this.backward && Math.abs(this.hField - this.finalValue) < 0.5*Math.abs(this.step)){
					this.backward = true;
					this.step = -this.step;
				}
			case 'HSWEEP':
				this.hField += this.step;
				break;
		}
	}
	
	this.exitCondition = function(){
		switch(this.loopMode){
			case 'TSWEEP':
				return Math.abs(this.temp - this.finalValue) < 0.5*Math.abs(this.step);
			case 'HSWEEP':
				return Math.abs(this.hField - this.finalValue) < 0.5*Math.abs(this.step);
			case 'TLOOP':
				return this.backward && Math.abs(this.temp - this.initialValue) < 0.5*Math.abs(this.step);
			case 'HLOOP':
				return this.backward && Math.abs(this.hField - this.initialValue) < 0.5*Math.abs(this.step);
			default:
				return this.iter >= this.finalValue;
		}
	}
	
	this.terminate = function(){
		console.log( '------ SIMULATION FINISHED ------\n' +
					 '  Execution time:    ' + (new Date().getTime() - this.timeStamp) + ' ms\n' +
					 '  Monte carlo steps: ' + this.iter * this.mcSteps 						);
					 
		if(this.dataFile != null){
			saveTable(this.dataTable, this.dataFile + '.csv');
			console.log('  Simulation data saved as ' + this.dataFile + '.csv');
		}
		noLoop();
	}
	
}