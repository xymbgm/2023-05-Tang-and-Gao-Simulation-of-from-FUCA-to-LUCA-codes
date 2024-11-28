/***
* Name: FUCA
* Author: Gao Ming
***/

model FUCA

global {
	Pond pond;
	float pContact; 
	list<Protocell> protocellGroup;  
	list<Protocell> tempProtocells;  
	list<Protocell> FUCAGroup;  
	matrix<float> mergeP <- 0.0 as_matrix({6,6});
		
	init {
		mergeP[0, 0] <- 0.5; mergeP[1, 0] <- 0.0; mergeP[2, 0] <- 0.0; mergeP[3, 0] <- 0.0; mergeP[4, 0] <- 0.0; mergeP[5, 0] <- 0.0;  
		mergeP[0, 1] <- 1.0; mergeP[1, 1] <- 0.5; mergeP[2, 1] <- 0.4; mergeP[3, 1] <- 0.3; mergeP[4, 1] <- 0.2; mergeP[5, 1] <- 0.1; 
		mergeP[0, 2] <- 1.0; mergeP[1, 2] <- 0.6; mergeP[2, 2] <- 0.5; mergeP[3, 2] <- 0.4; mergeP[4, 2] <- 0.3; mergeP[5, 2] <- 0.2; 
		mergeP[0, 3] <- 1.0; mergeP[1, 3] <- 0.7; mergeP[2, 3] <- 0.6; mergeP[3, 3] <- 0.5; mergeP[4, 3] <- 0.4; mergeP[5, 3] <- 0.3; 
		mergeP[0, 4] <- 1.0; mergeP[1, 4] <- 0.8; mergeP[2, 4] <- 0.7; mergeP[3, 4] <- 0.6; mergeP[4, 4] <- 0.5; mergeP[5, 4] <- 0.4;
		mergeP[0, 5] <- 1.0; mergeP[1, 5] <- 0.9; mergeP[2, 5] <- 0.8; mergeP[3, 5] <- 0.7; mergeP[4, 5] <- 0.6; mergeP[5, 5] <- 0.5;
		
		create Pond number: 1 returns: pondGroup;
		pond <- first(pondGroup);
		pContact <- rnd(10^-5,10^-6)*((100.0 / pond.currentVolume) ^ 2);  
		write "Volumn of pond is " + pond.currentVolume;
		write "pContact is " + pContact;
		
		create V1 number: pond.v1Num returns: v1Group;
		create V2 number: pond.v2Num returns: v2Group;
		write "Number of Vesicles: " + pond.v1v2Num;
		write "Number of V1: " + length(v1Group);
		write "Number of V2: " + length(v2Group);	
		
		protocellGroup <<+ v1Group;
		protocellGroup <<+ v2Group;
		write "Number of protocellGroup: " + length(protocellGroup);
		loop i over: protocellGroup {
			write i;
		}
		
		protocellGroup <- shuffle(protocellGroup);
		write "After shuffling: ";
		loop i over: protocellGroup {
			write i;
		}	
		
		write "Initiation is done!";
		write "===========================================\n";
	}
	
	reflex WetDry {
		write "\nCycle " + cycle + "\n";
		write "Volumn of pond is " + pond.currentVolume;
		write "pContact is " + pContact;
		write "Number of protocellGroup: " + length(protocellGroup);
		
		tempProtocells <- [];
		
		loop i from: 0 to: length(protocellGroup) - 1 {
			int j <- i + 1;
			if (j <= (length(protocellGroup) - 1)) {
				loop k from: j to: length(protocellGroup) - 1 {
					if(flip(pContact) and protocellGroup[i].alive and protocellGroup[k].alive) {
						float pairP <- mergeP[protocellGroup[i].type - 1, protocellGroup[k].type - 1];  
						if(flip(pairP)) { 
							int peptideSum <- protocellGroup[i].peptide + protocellGroup[k].peptide;
							int RNASum <- protocellGroup[i].RNA + protocellGroup[k].RNA;
							if(peptideSum <= 3) { 
								create V1 number: 1 returns: newV1;
								first(newV1).peptide <- peptideSum;
								first(newV1).RNA <- RNASum;
								add first(newV1) to: tempProtocells;
								protocellGroup[i].alive <- false;
								protocellGroup[k].alive <- false;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV1);
							} 
							if(peptideSum <= 10 and peptideSum >= 4) { 
								create V2 number: 1 returns: newV2;
								first(newV2).peptide <- peptideSum;
								first(newV2).RNA <- RNASum;
								add first(newV2) to: tempProtocells;
								protocellGroup[i].alive <- false;
								protocellGroup[k].alive <- false;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV2);
							}
							if(peptideSum <= 20 and peptideSum >= 11) { 
								create V3 number: 1 returns: newV3;
								first(newV3).peptide <- peptideSum;
								first(newV3).RNA <- RNASum;
								add first(newV3) to: tempProtocells;
								protocellGroup[i].alive <- false;
								protocellGroup[k].alive <- false;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV3);
							}
							if(peptideSum <= 30 and peptideSum >= 21) { 
								create V4 number: 1 returns: newV4;
								first(newV4).peptide <- peptideSum;
								first(newV4).RNA <- RNASum;
								add first(newV4) to: tempProtocells;
								protocellGroup[i].alive <- false;
								protocellGroup[k].alive <- false;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV4);
							}
							if(peptideSum <= 40 and peptideSum >= 31) { 
								create V5 number: 1 returns: newV5;
								first(newV5).peptide <- peptideSum;
								first(newV5).RNA <- RNASum;
								add first(newV5) to: tempProtocells;
								protocellGroup[i].alive <- false;
								protocellGroup[k].alive <- false;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV5);
							}
							if(peptideSum >= 41) { 
								create V6 number: 1 returns: newV6;
								first(newV6).peptide <- peptideSum;
								first(newV6).RNA <- RNASum;
								add first(newV6) to: tempProtocells;
								protocellGroup[i].alive <- false;
								protocellGroup[k].alive <- false;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV6);
								write "FUCA: " + first(newV6) + " appears!" color: #red;
							}
						}
					}
				}
			}
		}
		
		
		if (length(protocellGroup select (each.alive = false)) > 0) {
			list<Protocell> readyToKillGroup <- [];  
			readyToKillGroup <<+ protocellGroup select (each.alive = false);
			
			protocellGroup >>- protocellGroup select (each.alive = false);

			loop p over: readyToKillGroup {
				ask p {
					do Kill;
				}
			}
		}
		
		protocellGroup <<+ tempProtocells;
		
		pond.currentVolume <- rnd(50, 100);
		
		pond.v1v2Num <- rnd(50, 100);  
		pond.v1Num <- int(rnd(0.55, 0.6) * pond.v1v2Num);
		pond.v2Num <- pond.v1v2Num - pond.v1Num; 
		create V1 number: pond.v1Num returns: v1Group;
		create V2 number: pond.v2Num returns: v2Group;
		write "" + pond.v1v2Num + " vesicles are produced."; 
				
		protocellGroup <<+ v1Group;
		protocellGroup <<+ v2Group;
	
		protocellGroup <- shuffle(protocellGroup);
		
		pContact <- rnd(10^-5, 10^-6)*((100.0 / pond.currentVolume) ^ 2);  
	}
	
	reflex WriteToFile {
		if((V6 count(each.alive = true)) >= 300) {
			loop m from: 0 to: length(V6) - 1 {
				if(V6[m].alive = true) {
					add V6[m] to: FUCAGroup;
				}
			}	
			save [1, FUCAGroup[0], FUCAGroup[0].peptide, FUCAGroup[0].RNA] to: "/output/FUCA.csv" type:"csv" rewrite: true;
			loop n from: 1 to: length(FUCAGroup) - 1 {  
				save [n+1, FUCAGroup[n], FUCAGroup[n].peptide, FUCAGroup[n].RNA] to: "/output/FUCA.csv" type:"csv" rewrite: false;
			}				
	
			do pause;
														
		}			
	}
}


species V1 parent: Protocell {
	init {
		type <- 1;
		peptide <- rnd(2, 3);
		RNA <- rnd(5, 10);
	}	
}

species V2 parent: Protocell {
	init {
		type <- 2;
		peptide <- rnd(4, 10);
		RNA <- rnd(10, 30);
	}
}

species V3 parent: Protocell {
	init {
		type <- 3;
	}
}

species V4 parent: Protocell {
	init {
		type <- 4;
	}
}

species V5 parent: Protocell {
	init {
		type <- 5;
	}
}

species V6 parent: Protocell {
	init {
		type <- 6;
	}
}

species Protocell {
	int peptide <- 0;
	int RNA <- 0;
	int type; 
	bool merge <- false;
	bool alive <- true; 
	
	action Kill {
		do die;
	}
}

species Pond {
	int currentVolume;  
	int v1v2Num; 
	int v1Num; 
	int v2Num; 
	 
	init {
		currentVolume <- rnd(50, 100);  
		v1v2Num <- rnd(50, 100);  
		v1Num <- int(rnd(0.55, 0.6) * v1v2Num);
		v2Num <- v1v2Num - v1Num; 
	}
}

experiment FUCA type: gui {
	output {
		display Chart type:java2D {
			chart "FUCA" type: series x_label: "Cycle(s)"  y_label: "Vesicle Number" {
				data "V1" value: V1 count(each.alive = true) color: #violet;
				data "V2" value: V2 count(each.alive = true) color: #green;
				data "V3" value: V3 count(each.alive = true) color: #orange;
				data "V4" value: V4 count(each.alive = true) color: #lightgreen;
				data "V5" value: V5 count(each.alive = true) color: #blue;
				data "FUCA" value: V6 count(each.alive = true) color: #red;
			}
		}
				
		monitor "Volumn of Pond" value: pond.currentVolume;
		monitor "Number of V1" value: V1 count(each.alive = true);
		monitor "Number of V2" value: V2 count(each.alive = true);
		monitor "Number of V3" value: V3 count(each.alive = true);
		monitor "Number of V4" value: V4 count(each.alive = true);
		monitor "Number of V5" value: V5 count(each.alive = true);
		monitor "Number of FUCA" value: V6 count(each.alive = true);
		monitor "Number in total" value: V1 count(each.alive = true)+V2 count(each.alive = true)+V3 count(each.alive = true)+V4 count(each.alive = true)+V5 count(each.alive = true)+V6 count(each.alive = true);
	}	
}



