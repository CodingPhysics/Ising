/********* MONTE CARLO SIMULATION ON AN ISING SPIN SYSTEM ON A TWO-DIMENSIONAL PERIODIC LATTICE **************/

/*************************************************************************************************************
From: 	Coding Phyisics
email:	CodingPhysicSimulation@gmail.com
github:	https://github.com/CodingPhysics
**************************************************************************************************************/

/***************************************** D E S C R I P T I O N *********************************************
This is a P5.js animation of a simple Ising model for a ferromagnet on a two-dimensional square lattice.
This program visualizes magnetic domains and can simulate temperature or magnetic field curves. 
If wanted, the simulation data and indiviual frames are automatically saved for later evaluation or presentation. 
(In particular, the captured frames can be used to create a video animation.)
*************************************************************************************************************/

/************************************** P H Y S I C A L   M O D E L ******************************************
The state of a Ising ferromagnet is described by a dicrete spin vector S 

	S = (s_1, s_2, ..., s_N) ∈ {-1, +1}^N

where each spin s_i can assume only two possible values -1 or +1 (spin-down or spin-up).
The energy of a certain state is given by the dimensionless Hamilitonian (J = mu = 1):

	H(S) = - ∑_<i,j> s_i s_j - h ∑_i s_i = - ∑_i (h_i + h) s_i

where the first sum is taken over all pairs i,j of adjacent spins. 
h denotes the external magnetic field and h_i is the so-called exchange field.
In thermodynamic equilibrium, the systems shows a Boltzmann distribution (k = 1):

	p(S) = const * exp( - H(S) / T) 

where T is the temperature. 
*************************************************************************************************************/

/********************************** M E T R O P O L I S   A L G O R I T M *********************************** 
The Ising system is simulated by a Markov chain of spin states S(t) as a function of discrete time t, 
which fulfills the Boltzmann distributions. The algorithm can summarized by the following:

1. Choose an initial spin configuration S(0)
2. Choose a random i and calculate the energy change for a single spin-flip s_i → -s_i:

	∆E = - 2 s_i (h_i + h)
	
3. Attempt a spin-flip with the transition probability
	                ┌
	                │ 1          	∆E ≥ 0
	p(s_i → -s_i) = ┤
	                │ exp(∆E / T)	∆E < 0
	                └
   For this, generate a random number r ∈ [0,1] and flip the i-th spin, if r < p(s_i → -s_i).
   Otherwise keep the current state.
4. Go to 2. and repeat

************************************************************************************************************/

/********************************************** U S A G E ****************************************************
Ising is called with a settings object as a single argument, that summarizes all parameters of the simulation:

var mySettings = { 
	width:            <Number>                                // width of the spin lattice 
	height:           <Number>,                               // height of the spin lattice
	mcStepsPerFrame:  <Number>,                               // number of iterations in a simulation step
	averagingFrames:  <Number>,                               // number of simulation steps for averaging 
	temperature:      <Number>,                               // initial temperature 
	magneticField:    <Number>,                               // initial magnetic field 
	randomizeLattice: <Boolean>,                              // initial spin configuration 
	loopMode:         'TSWEEP'|'TLOOP'|'HSWEEP'|'HLOOP'|null, // specifies temperature/field curve
	loopTargetValue:  <Number>,                               // target value for temperature/field in the loop
	loopIncrement:    <Number>,                               // step size of temperature/field in the loop
	imageFile:        <String>|null,                          // file name for capturing of the animation frames
	dataFile:         <String>|null                           // file name for the simulation data
};

var mySimulation = new Ising(mySettings);

Only the parameters imageFile, dataFile and loopMode are optional and should be set to null if not used. 
In order to excute and animate the simulation, the function simulationStep() must be called in the draw function:

function draw() {
	mySimulation.simulationStep();
}

Image and data files are automatically saved to the local download directory.
***************************************************************************************************************/


function Ising(settings){
	
	this.width            = settings.width;            
	this.height           = settings.height;           
	this.mcStepsPerFrame  = settings.mcStepsPerFrame; 
	this.averagingFrames  = settings.averagingFrames;  
	this.temperature      = settings.temperature;      
	this.magneticField    = settings.magneticField;    
	this.randomizeLattice = settings.randomizeLattice; 
	this.loopMode         = settings.loopMode;         
	this.loopTargetValue  = settings.loopTargetValue;  
	this.loopIncrement    = settings.loopIncrement;    
	this.imageFile        = settings.imageFile;        
	this.dataFile         = settings.dataFile;         
	
	console.log( '----- SIMULATION PARAMETERS -----'                            + '\n' +
	             '  Lattice size:            ' + this.width + 'x' + this.height + '\n' +
	             '  MC steps per interation: ' + this.mcStepsPerFrame           + '\n' +
	             '  Interations for average: ' + this.averagingFrames           + '\n' +
	             '  Initial temperature:     ' + this.temperature               + '\n' +
	             '  Initial magnetic field:  ' + this.magneticField             + '\n' +
	             '  Random initialization:   ' + this.randomizeLattice          + '\n' +
	             '  Loop type:               ' + this.loopMode                  + '\n' +
	             '  Target value of loop:    ' + this.loopTargetValue           + '\n' +
	             '  Increment:               ' + this.loopIncrement             + '\n' +
		     '---------------------------------'                               );
	
	this.frameCounter  = 1;                    
	this.startingTime  = new Date().getTime(); 
	this.dataTable     = new p5.Table();       
	this.dataTable.addColumn('temperature');
	this.dataTable.addColumn('magnetic field');
	this.dataTable.addColumn('magnetization');
	this.dataTable.addColumn('order parameter');
	
	this.loopBackward   = false; 
	this.loopStartValue = (this.loopMode.charAt(0) == 'T')? this.temperature : this.magneticField; 

	this.magnetization         = 0;     // magnetization of the spin-lattice 
	this.orderParameter        = 0;     // local order parameter
	this.averageMagnetization  = 0;     // average magnetization
	this.averageOrderParameter = 0;     // average order parameter
	this.spin                  = [];    // spin matrix (spin[i][j] = -1, +1)
	this.spinBuffer            = [];    // buffer for faster image drawing
	this.exchangeField         = [];    // exchange field (sum of adjacent spins)
	
	
	this.initialize = function(){
		
		this.spin          = [];
		this.spinBuffer    = [];
		this.exchangeField = [];
		
		for(var i = 0; i < this.width; i++){
			this.spin.push([]);
			this.spinBuffer.push([]);
			
			for(var j = 0; j < this.height; j++){
				var s = (this.randomizeLattice)? 2*Math.floor(2*Math.random()) - 1 : 1;
				this.spin[i].push(s);
				this.spinBuffer[i].push(0);
			}
		}
		
		for(var i = 0; i < this.width; i++){
			this.exchangeField.push([]);
			
			for(var j = 0; j < this.height; j++){
				var right = (i + 1) % this.width;
				var left  = (i - 1 + this.width) % this.width;
				var above = (j + 1) % this.height;
				var below = (j - 1 + this.height) % this.height;
				
				var sumOverNextNeighbors = this.spin[right][j] + this.spin[left][j] + 
				                           this.spin[i][above] + this.spin[i][below] ;
				
				this.exchangeField[i].push(sumOverNextNeighbors);
			}
		}
	}
	this.initialize();
	
	
	this.simulationStep = function(){
		console.log('ITERATION ' + this.frameCounter);
		
		this.metropolis();	// execute Markov chain steps
		this.update();		// update physical quantities
		this.show();		// display spin lattice
		this.iterate();		// increment temperature or field and repeat while unfinished
	}
	
	
	this.metropolis = function(){
		console.log('-> Metropolis algorithm');
		
		for(var mcSteps = 0; mcSteps < this.mcStepsPerFrame; mcSteps++){
			// choose a random cell:
			var i = Math.floor(this.width*Math.random());
			var j = Math.floor(this.height*Math.random());
			
			// calculate energy for spin-flip:
			var deltaE = - 2*this.spin[i][j]*(this.exchangeField[i][j] + this.magneticField);
			
			// attempt spin-flip with certain transtion probability:
			if(deltaE >= 0 || Math.random() < Math.exp(deltaE/this.temperature)){
				// flip spin:
				this.spin[i][j] = -this.spin[i][j];
				
				// update exchange field of adjacent spins:
				
				var right = (i + 1) % this.width;
				var left  = (i - 1 + this.width) % this.width;
				var above = (j + 1) % this.height;
				var below = (j - 1 + this.height) % this.height;
				
				this.exchangeField[right][j] += 2*this.spin[i][j];
				this.exchangeField[left][j]  += 2*this.spin[i][j];
				this.exchangeField[i][above] += 2*this.spin[i][j];
				this.exchangeField[i][below] += 2*this.spin[i][j];
			}
		}
	}
	
	
	this.update = function(){
		console.log('-> Update physical quantities');
		
		// calculate magnetization:
		
		var weight = 1./this.width/this.height;
		this.magnetization = 0;
		
		for(var i = 0; i < this.width; i++){
			for(var j = 0; j < this.height; j++){
				this.magnetization += this.spin[i][j]*weight;
			}
		}
		
		this.averageMagnetization += this.magnetization/this.averagingFrames;
		
		// calculate local order parameter:
		
		weight /= 4;
		this.orderParameter = 0;
		
		for(var i = 0; i < this.width; i++){
			for(var j = 0; j < this.height; j++){
				this.orderParameter += this.exchangeField[i][j]*this.spin[i][j]*weight;
			}
		}
		
		this.averageOrderParamter += this.orderParameter/this.averagingFrames;
		
		// update datatable if averaging is completed:
		
		if(this.frameCounter % this.averagingFrames == 0){	
			this.updateDataTable();
			this.averageMagnetization  = 0;
			this.averageOrderParameter = 0; 
		}
	}
	
	this.updateDataTable = function(){
		var newRow = this.dataTable.addRow();
		newRow.setString('temperature'    , this.temperature.toExponential(6)           ); 
		newRow.setString('magnetic field' , this.magneticField.toExponential(6)         ); 
		newRow.setString('magnetization'  , this.averageMagnetization.toExponential(6)  ); 
		newRow.setString('order parameter', this.averageOrderParameter.toExponential(6) );
			
		console.log( '------ UPDATED DATATABLE ------'                  + '\n' +
		             '  Temperature:     ' + this.temperature           + '\n' +
		             '  Magnetic field:  ' + this.magneticField         + '\n' +
		             '  Magnetization:   ' + this.averageMagnetization  + '\n' +
		             '  Order parameter: ' + this.averageOrderParameter + '\n' +
		             '-------------------------------'                           );
	}
		
	this.show = function(){
		console.log('-> Display spin lattice');
		
		// draw spin lattice:
		
		var boxSize = height/this.height;		  
		var margin = (width - boxSize*this.width)/2;    
		
		noStroke();
		for(var i = 0; i < this.width; i++){
			for(var j = 0; j < this.height; j++){
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
		
		this.drawScale(this.temperature,     0, 4, 'T',   0.3*margin,         Math.floor(margin/10), 5);
		this.drawScale(this.magneticField,  -1, 1, 'H',   0.7*margin,         Math.floor(margin/10), 5);
		this.drawScale(this.magnetization,  -1, 1, 'M',   width - 0.7*margin, Math.floor(margin/10), 5);
		this.drawScale(this.orderParameter, -1, 1, 'eta', width - 0.3*margin, Math.floor(margin/10), 5);
		
		// save frame as png image:
		
		if(this.imageFile !== null){
			saveCanvas(this.imageFile + this.iter + '.png');
			console.log('   Frame saved as ' + this.imageFile + this.iter + '.png');
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
	
	
	this.iterate = function(){
		
		if(this.frameCounter % this.averagingFrames == 0){ // if averaging is completed
			if(this.exitCondition()){
				this.terminate();          // end simulation ..
			} else {								
				this.increment();          // ... or change temperature/field and continue
			}
		}
		this.frameCounter++;
	}
	
	
	this.increment = function(){
		switch(this.loopMode){
			case 'TLOOP':
				if(!this.loopBackward && Math.abs(this.temperature - this.loopTargetValue) < 1E-6){
					this.loopBackward = true;
					this.loopIncrement = -this.loopIncrement;
				}
			case 'TSWEEP':
				this.temperature += this.loopIncrement;
				break;
			case 'HLOOP':
				if(!this.loopBackward && Math.abs(this.magneticField - this.loopTargetValue) < 1E-6){
					this.loopBackward = true;
					this.loopIncrement = -this.loopIncrement;
				}
			case 'HSWEEP':
				this.magneticField += this.loopIncrement;
				break;
		}
	}
	
	
	this.exitCondition = function(){
		switch(this.loopMode){
			case 'TSWEEP':
				return Math.abs(this.temperature - this.loopTargetValue) < 1E-6;
			case 'HSWEEP':
				return Math.abs(this.magneticField - this.loopTargetValue) < 1E-6;
			case 'TLOOP':
				return this.loopBackward && Math.abs(this.temperature - this.loopStartValue) < 1E-6;
			case 'HLOOP':
				return this.loopBackward && Math.abs(this.magneticField - this.loopStartValue) < 1E-6;
			default:
				return this.frameCounter >= this.loopTargetValue;
		}
	}
	
	
	this.terminate = function(){
		console.log( '------ SIMULATION FINISHED ------\n' +
		             '  Execution time:    ' + (new Date().getTime() - this.startingTime) + ' ms\n' +
		             '  Monte carlo steps: ' + this.frameCounter * this.mcSteps 			);
					 
		if(this.dataFile !== null){
			saveTable(this.dataTable, this.dataFile + '.csv');
			console.log('  Simulation data saved as ' + this.dataFile + '.csv');
		}
		noLoop();
	}
	
}
