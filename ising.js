/* MONTE CARLO SIMULATION ON AN ISING SPIN SYSTEM ON A TWO-DIMENSIONAL PERIODIC LATTICE */

/*************************************************************************************
From: 	Coding Phyisics
email:	CodingPhysicSimulation@gmail.com
github:	https://github.com/CodingPhysics
*************************************************************************************/

function Ising(settings){
	
	this.width            = settings.width;            // width of the spin lattice
	this.height           = settings.height;           // height of the spin lattice
	this.mcStepsPerFrame  = settings.mcStepsPerFrame;  // number of Markov chain steps in a simulation step
	this.averagingFrames  = settings.averagingFrames;  // number of simulations step for averaging of measured quantites
	this.temperature      = settings.temperature;      // initial temperature (in units of J/k)
	this.magneticField    = settings.magneticField;    // initial magnetic field (in units of J/mu)
	this.randomizeLattice = settings.randomizeLattice; // initial spin configuration (true: random state, false: high-spin state)
	this.loopMode         = settings.loopMode;         // 'TSWEEP', 'TLOOP', 'HSWEEP', 'HLOOP' for temperature/field curves
	this.loopTargetValue  = settings.loopTargetValue;  // target value for temperature/field in the loop
	this.loopIncrement    = settings.loopIncrement;    // step size of temperature/field in the loop
	this.imageFile        = settings.imageFile;        // file name for capturing of the animation frames
	this.dataFile         = settings.dataFile;         // file name for the simulation data
	
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
	this.loopStartValue = (loopMode.charAt(0) == 'T')? this.temperature : this.magneticField; 

	this.magnetization         = 0;     // magnetization of the spin-lattice 
	this.orderParameter        = 0;     // local order parameter
	this.averageMagnetization  = 0;     // average magnetization
	this.averageOrderParameter = 0;     // average order parameter
	this.spin                  = [];    // spin matrix (spin[i][j] = -1, +1)
	this.spinBuffer            = [];    // buffer for faster image drawing
	this.exchangeField         = [];    // exchange field (sum of adjacent spins)
	
	
	this.initialize = function(){
		
		this.spin       = [];
		this.spinBuffer = [];
		
		for(var i = 0; i < this.width; i++){
			this.spin.push([]);
			this.spinBuffer.push([]);
			
			for(var j = 0; j < this.height; j++){
				var s = (this.randomizeLattice)? 2*Math.floor(2*Math.random()) - 1) : 1;
				this.spin[i].push(s);
				this.spinBuffer[i].push(0);
			}
		}
		
		this.exchangeField = [];
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
		
		for(var i = 0; i < this.nx; i++){
			for(var j = 0; j < this.ny; j++){
				this.orderParameter += this.exchangeField[i][j]*this.spin[i][j]*weight;
			}
		}
		
		this.averageOrderParamter += this.orderParameter/this.averagingFrames;
		
		// update datatable if averaging is completed:
		
		if(this.frameCounter % this.averagingFrames == 0){	
			var newRow = this.dataTable.addRow();
			newRow.setString('temperature'    , this.temperature.toExponential(6)           ); 
			newRow.setString('magnetic field' , this.magneticField.toExponential(6)         ); 
			newRow.setString('magnetization'  , this.averageMagneticField.toExponential(6)  ); 
			newRow.setString('order parameter', this.averageOrderParameter.toExponential(6) );
			
			console.log( '------ UPDATED DATATABLE ------'                  + '\n' +
			             '  Temperature:     ' + this.temperature           + '\n' +
			             '  Magnetic field:  ' + this.magneticField         + '\n' +
			             '  Magnetization:   ' + this.averageMagnetization  + '\n' +
			             '  Order parameter: ' + this.averageOrderParameter + '\n' +
			             '-------------------------------'                           );
			this.averageMagnetization  = 0;
			this.averageOrderParameter = 0; 
		}
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
		
		if(this.imageFile != null){
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
					 
		if(this.dataFile != null){
			saveTable(this.dataTable, this.dataFile + '.csv');
			console.log('  Simulation data saved as ' + this.dataFile + '.csv');
		}
		noLoop();
	}
	
}
