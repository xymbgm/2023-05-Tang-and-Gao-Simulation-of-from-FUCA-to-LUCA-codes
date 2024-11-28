/**
* Name: LUCA(300FUCA)
* Author: Gao Ming
**/


model LUCA

global {
	Pond pond;
	float pContact; 
	list<Protocell> protocellGroup;  
	int minPeptide;
	int maxPeptide;
	int minRNA;
	int maxRNA;
	bool isPause <- false;

	init {
		create Pond number: 1 returns: pondGroup;
		pond <- first(pondGroup);
		pond.currentVolume <- 100;
		
		create V6 number: 300 returns: newV6 {
			do AAMatrixInit;
		}
		protocellGroup <<+ newV6;
		
		list<int> peptideCount <- [];
		list<int> RNACount <- [];
		
		loop i over: V6 {
			peptideCount <+ i.numPeptide;
			RNACount <+ i.numRNA;			
		}
		
		minPeptide <- min(peptideCount);
		maxPeptide <- max(peptideCount);
		minRNA <- min(RNACount);
		maxRNA <- max(RNACount);
		
		save [] to: "/output/LUCA.csv" type:"csv" rewrite: false;		
		file lucaFile <- csv_file("/output/LUCA.csv", ",");
		matrix data <- matrix(lucaFile);
		
		if (length(data) = 0) {
			save ["Simulation", "MinPeptide", "MaxPeptide", "MinRNA", "MaxRNA", "NumPeptideOfLUCA", "NumRNAOfLUCA", "Cycles", "LUCA"] to: "/output/LUCA.csv" type:"csv" rewrite: false;			
		}			
		
		write "Initiation is done!\n";
		write "===========================================";
	}		

	reflex WetOrDry {
		write "\nCycle " + cycle + "\n";
		if (cycle mod 2 = 0) {
			ask pond {
				isWet <- true;
				write "Pond is in Wet Phase.";
			
				if (cycle != 0) {
					do GetWetPhaseVolume;
				}
			}
			write "Current volume:" + pond.currentVolume;
		}
		else {
			ask pond {
				isWet <- false;
				write "\nPond is in Dry Phase";
				
				do GetDryPhaseVolume;
			}
			write "Current volume:" + pond.currentVolume;			
		}
	}
	
	reflex Contact {
		protocellGroup <- shuffle(protocellGroup); 
		
		list<Protocell> tempProtocellGroup <- []; 
		
		loop i from: 0 to: length(protocellGroup) - 1 {
			int j <- i + 1;
			if (j <= (length(protocellGroup) - 1)) {
				loop k from: j to: length(protocellGroup) - 1 {
					if (flip(pContact) and protocellGroup[i].isAlive and protocellGroup[k].isAlive) {
						float pPair <- 0.5; 
						if (flip(pPair)) { 							
							create V6 number: 1 returns: newV6;							
							first(newV6).numPeptide <- protocellGroup[i].numPeptide + protocellGroup[k].numPeptide;
							first(newV6).numRNA <- protocellGroup[i].numRNA + protocellGroup[k].numRNA;							
											 
							loop m from: 0 to: 19 {
								if (protocellGroup[i].AAMatrix[m] = 1) {
									first(newV6).AAMatrix[m] <- protocellGroup[i].AAMatrix[m]; 
									first(newV6).AACodon[m] <- protocellGroup[i].AACodon[m];
								}
								if (protocellGroup[k].AAMatrix[m] = 1) {
									first(newV6).AAMatrix[m] <- protocellGroup[k].AAMatrix[m]; 
									first(newV6).AACodon[m] <- protocellGroup[k].AACodon[m]; 
								}
							}
							
							write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV6);
							
							add first(newV6) to: tempProtocellGroup;							
							protocellGroup[i].isAlive <- false;
							protocellGroup[k].isAlive <- false;		
												
							int numSynchronized <- 0;
							loop m from: 0 to: 19 {
								numSynchronized <- numSynchronized + first(newV6).AACodon[m];
							}
							if (numSynchronized = 20) {
								first(newV6).isSynchronized <- true;
							}									
						}
					}
				}
			}
		}		
	
		if (length(protocellGroup select (each.isAlive = false)) > 0) {
			list<Protocell> readyToKillGroup <- [];  
			readyToKillGroup <<+ protocellGroup select (each.isAlive = false);
			
			protocellGroup >>- protocellGroup select (each.isAlive = false);	
			
			loop p over: readyToKillGroup {
				ask p {
					do Kill;
				}
			}
		}	
		
		protocellGroup <<+ tempProtocellGroup;
	}
	
	reflex GetOtherAAs when: length(V6) > 0 and (cycle mod 2 != 0) {
		write "isAlive:" + V6 count(each.isAlive = true);
		write "isNotAlive:" + V6 count(each.isAlive = false);
		
		loop i over: V6 {
			ask i {			
				do CheckGetAllAAs;
				do SynthesiseLastTenAAs;						
			}
		}
	}
	
	reflex AssignAAWithCodons when: length(V6) > 0 and (cycle mod 2 != 0) {	
		loop i over: V6 {
			ask i {				
				if (isGetAllAAs and !isSynchronized) {
					do SynchronizeAACodon;
				}
			}
		}
	}
	
	reflex checkPeptideRNA {
		list<int> peptideCount <- [];
		list<int> RNACount <- [];
		
		loop i over: V6 {
			peptideCount <+ i.numPeptide;
			RNACount <+ i.numRNA;			
		}
		
		minPeptide <- min(peptideCount);
		maxPeptide <- max(peptideCount);
		minRNA <- min(RNACount);
		maxRNA <- max(RNACount);
	}

	reflex CheckLUCA when: length(V6) > 0 {
		file lucaFile <- csv_file("/output/LUCA.csv", ",");
		matrix data <- matrix(lucaFile);
		
		loop i over: V6 {
			ask i {
				if(isSynchronized) {
					i.isLUCA <- true;
					write "LUCA(" + i + ") appears!"color: #red;
					
					loop j from: 0 to: 19 {
						write "AACodon(" + (j + 1) + "): " + AACodon[j];
					}
					
					isPause <- true;
					
					if (length(data) > 1) {
						save [int(data[0, data.rows -1]) + 1, minPeptide, maxPeptide, minRNA, maxRNA, i.numPeptide, i.numRNA, cycle, i] to: "/output/LUCA.csv" type:"csv" rewrite: false;	
					}
					else {
						save [1, minPeptide, maxPeptide, minRNA, maxRNA, i.numPeptide, i.numRNA, cycle, i] to: "/output/LUCA.csv" type:"csv" rewrite: false;		
					}					
				}
			}
		}
		
		if(isPause) {
			do pause;
		}
	}  
}

species Pond {
	bool isWet <- false;
	int currentVolume;

	action GetDryPhaseVolume {
		currentVolume <- rnd(50, 80);
		pContact <- rnd(10 ^ -5, 10 ^ -6) * ((100.0 / currentVolume) ^ 2);  	
	}
	
	action GetWetPhaseVolume {
		currentVolume <- rnd(80, 100);
		pContact <- rnd(10 ^ -5, 10 ^ -6) * ((100.0 / currentVolume) ^ 2);  
	}
}

species Protocell {
	int type; 
	bool isAlive <- true;
	int numPeptide <- rnd(100, 150);  
	int numRNA <- rnd(300, 450);   
	matrix<int> AAMatrix <- 0 as_matrix({20, 1});
	bool isGetAllAAs <- false;  
	matrix<int> AACodon <- 0 as_matrix({20, 1});  
	bool isSynchronized <- false;  
	bool isLUCA <- false;
	
	action SynthesiseLastTenAAs {
		bool isTrue <- true;  
		 
		if (!isGetAllAAs) {  
			loop while: isTrue {
				int rndIndexAA <- rnd(11, 20);
				 				
				if(AAMatrix[rndIndexAA - 1] = 0) {
					switch rndIndexAA {
						match 11 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}
						match 12 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}
						match 13 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}
						match 14 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}
						match 15 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}
						match 16 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}
						match 17 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}
						match 18 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}
						match 19 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}
						match 20 {
							AAMatrix[rndIndexAA - 1] <- 1;
							isTrue <- false;
						}						
					}			
				}
			} 		
		}	
	}
	
	action CheckGetAllAAs {  
		if (!isGetAllAAs) {
			int tmpSum;
			
			loop i from: 10 to: 19 {
				tmpSum <- tmpSum + AAMatrix[i];
			}
			
			if (tmpSum = 10) {
				isGetAllAAs <- true;
				write "" + self + " gets all AAs.";
			}
			else {
				write "" + self + " doesn't get all AAs.";
			}
		}
	}
	
	action SynchronizeAACodon {
		bool isTrue <- true;  
		int numSynchronized <- 0; 
		float part;
		float pSynchronize;
		loop i from: 0 to: 19 {
			write "AACodon(" + (i + 1) + "): " + AACodon[i];
			numSynchronized <- numSynchronized + AACodon[i];
		}		
		//part <- 20 / (20 - numSynchronized);
		//pSynchronize <- rnd(10 ^ -3, 10 ^ -2) * (part ^ 2);
		pSynchronize <- 0.999;
		write "pSynchronize: " + pSynchronize;
		write "numSynchronized:" + numSynchronized;
		if(flip(pSynchronize) and numSynchronized < 20) {
			loop while: isTrue {
				int rndIndexAA <- rnd(1, 20);
				 				
				if(AACodon[rndIndexAA - 1] = 0) {
					AACodon[rndIndexAA - 1] <- 1;
					write "" + self + ".AA" + rndIndexAA + " has assigned with one set of codons.";
					isTrue <- false;			
				}
			}
			
			numSynchronized <- 0;
			loop i over: AACodon {
				numSynchronized <- numSynchronized + i;
			}
			if (numSynchronized = 20) {
				isSynchronized <- true;
			}
		}		
	}
	
	action Kill {
		do die;
	}
}

species V6 parent: Protocell {
	init {
		type <- 6;	
	}
	
	action AAMatrixInit {
		loop i from: 0 to: 9 {
			AAMatrix[i] <- 1;  
		}
	}
}

experiment LUCA type: gui{
	output {
		monitor "Volumn of Pond" value: pond.currentVolume;
		monitor "Number of FUCA" value: V6 count(each.isAlive = true);
		monitor "Number of LUCA" value: V6 count(each.isLUCA = true);
	}
}




