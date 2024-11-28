/**
* Name: LUCA(3FUCA)
* Author: Gao Ming
**/


model LUCA

global {
	float pSynchronize; 
	list<FUCA> vesicleGroup; 
	bool LUCAAppear <- false;
	
	init {
		create FUCA number: 3 returns: FUCAGroup;
				vesicleGroup <- FUCAGroup;					
		write "Number of VesicleGroup: " + length(vesicleGroup);
		loop i over: vesicleGroup {
			write i;			
			write "Length of FUCA(Proteins--RNAs): " + length(i.peptideComb) + "--" + length(i.RNAComb); 
			write "----------------";	
		}		
		write "Initiation is done!";
		write "===========================================\n";
	}
	
	reflex Process {
		write "\nCycles: " + cycle + "\n";
		write "Number of VesicleGroup: " + length(vesicleGroup);
		int pondVolume <- rnd(50, 100);
		float pContact <- rnd(10^-4,10^-5)*((100.0 / pondVolume) ^ 2);
		
		loop i from: 0 to: length(vesicleGroup) - 1 {
			loop j from: 0 to: length(vesicleGroup[i].peptideComb) - 1 {
				loop k from: 0 to: length(vesicleGroup[i].peptideComb[j].AAComb) - 1 {
					if(vesicleGroup[i].peptideComb[j].AAComb[k].key != vesicleGroup[i].peptideComb[j].codonComb[k].key) {
						if(vesicleGroup[i].peptideComb[j].accelerate = 0) {
							pSynchronize <- rnd(10^-2, 10^-3);
							if(flip(pSynchronize)) {
								vesicleGroup[i].peptideComb[j].AAComb[k].key <- vesicleGroup[i].peptideComb[j].codonComb[k].key;
								vesicleGroup[i].tightlyLinkingAAType[vesicleGroup[i].peptideComb[j].AAComb[k].type - 1] <- 1;	
								write "" + vesicleGroup[i] + "." + vesicleGroup[i].peptideComb[j] + "." + vesicleGroup[i].peptideComb[j].codonComb[k] + ": is synchronized.";	
								write "pSynchronize: " + pSynchronize;
								write "accelerate:" + vesicleGroup[i].peptideComb[j].accelerate;
							}					
						}
						else {
							ask vesicleGroup[i] {
								do GetTightlyLinkingAATypeCount;
								int n <- vesicleGroup[i].tightlyLinkingAATypeCount - 1;
								if(n < 0) {
									n <- 0;
								}
								float part <- 20 / (20 - n);
								pSynchronize <- rnd(10^-2, 10^-3)*(part^2);
								if(flip(pSynchronize)) {
									vesicleGroup[i].peptideComb[j].AAComb[k].key <- vesicleGroup[i].peptideComb[j].codonComb[k].key;
									vesicleGroup[i].tightlyLinkingAAType[vesicleGroup[i].peptideComb[j].AAComb[k].type - 1] <- 1;	
									write "" + vesicleGroup[i] + "." + vesicleGroup[i].peptideComb[j] + "." + vesicleGroup[i].peptideComb[j].codonComb[k] + ": is synchronized.";	
									write "pSynchronize: " + pSynchronize;
									write "accelerate:" + vesicleGroup[i].peptideComb[j].accelerate;						
								}
							}
						}					
					}
			  	}
			}
		}
			
		loop i from: 0 to: length(vesicleGroup) - 1 {
			int j <- i + 1;
			if (j < (length(vesicleGroup) - 1)) {
				loop k from: j to: length(vesicleGroup) - 1 {
					if(flip(pContact) and vesicleGroup[i].alive and vesicleGroup[k].alive) {
						loop m from: 0 to: length(vesicleGroup[k].peptideComb) - 1 {
							add vesicleGroup[k].peptideComb[m] to: vesicleGroup[i].peptideComb;
						}
						
						loop m from: 0 to: length(vesicleGroup[k].RNAComb) - 1 {
							add vesicleGroup[k].RNAComb[m] to: vesicleGroup[i].RNAComb;
						}
						
						write string(vesicleGroup[i]) + " and " + vesicleGroup[k] + " are merged.";
						vesicleGroup[k].alive <- false;
						//accelerate
						list<int> iHaveAAType <- list_with(20, 0);
						list<int> kHaveAAType <- list_with(20, 0);
						loop s from: 0 to: length(vesicleGroup[i].peptideComb) - 1 {
							loop t from: 0 to: length(vesicleGroup[i].peptideComb[s].AAComb) - 1 {
								iHaveAAType[vesicleGroup[i].peptideComb[s].AAComb[t].type - 1] <- 1;
							}
						}
						loop s from: 0 to: length(vesicleGroup[k].peptideComb) - 1 {
							loop t from: 0 to: length(vesicleGroup[k].peptideComb[s].AAComb) - 1 {
								kHaveAAType[vesicleGroup[k].peptideComb[s].AAComb[t].type - 1] <- 1;
							}
						}
						loop tmp from: 0 to: 19 {
							if(iHaveAAType[tmp] != kHaveAAType[tmp]) {
								loop u from: 0 to: length(vesicleGroup[i].peptideComb) - 1 {
									vesicleGroup[i].peptideComb[u].accelerate <- 1;
								}
								write string(vesicleGroup[i]) + ".peptideComb gets acceleration.";
								break;
							}
						}
						
						vesicleGroup[k].peptideComb <- [];
						vesicleGroup[k].RNAComb <- [];
					}						
				}				
			}			
		}

		bool proteinOK;
		bool RNAOK; 
		bool AAOK;
		bool NAOK;
		list<int> AATypeOK;
		int AATypeOKCount;
		bool AACodonMatchOK;
		
		loop i from: 0 to: length(vesicleGroup) - 1 {
			if(vesicleGroup[i].alive){
				LUCAAppear <- false;
				proteinOK <- false;
				RNAOK <- false;
				AAOK <- true;
				NAOK <- true;
				if(length(vesicleGroup[i].peptideComb) >= 100 and length(vesicleGroup[i].RNAComb) >= 300) {
					proteinOK <- true;
					RNAOK <- true;
					loop j from: 0 to: length(vesicleGroup[i].peptideComb) - 1 {
						if(length(vesicleGroup[i].peptideComb[j].AAComb) <= 50) {
							AAOK <- false;
							break;
						}
					}
					loop j from: 0 to: length(vesicleGroup[i].RNAComb) - 1 {
						if(length(vesicleGroup[i].RNAComb[i].NAComb) <= 150) {
							NAOK <- false;
							break;
						}
					}
				}
				
				ask vesicleGroup[i] {
					if(vesicleGroup[i].alive){
						do GetTightlyLinkingAATypeCount;
					}				
				}
	
				AACodonMatchOK <- true;
				int AACodonMatchFalseCount <- 0;
				loop s from: 0 to: length(vesicleGroup[i].peptideComb) - 1 {
					loop t from: 0 to: length(vesicleGroup[i].peptideComb[s].AAComb) - 1 {
						if(vesicleGroup[i].peptideComb[s].AAComb[t].key != vesicleGroup[i].peptideComb[s].codonComb[t].key) {
							AACodonMatchOK <- false;			
							AACodonMatchFalseCount <- AACodonMatchFalseCount + 1;
						}
						vesicleGroup[i].AACodonMatchFalseCount <- AACodonMatchFalseCount;
					}
				}
				if(proteinOK and RNAOK and AAOK and NAOK and vesicleGroup[i].tightlyLinkingAATypeCount = 20 and AACodonMatchOK) {
					LUCAAppear <- true;
					write "LUCA appears!!!" color: #red;
					write "LUCA: " + vesicleGroup[i];
					write "Number of Proteins: " + length(vesicleGroup[i].peptideComb);
					write "Number of RNAs: " + length(vesicleGroup[i].RNAComb);
					do pause;
				}
				else {				
					if(cycle mod 50 = 0) {
						write "" + vesicleGroup[i] + "(" + "ProteinOK: " + proteinOK + ", RNAOK: " + RNAOK + ", AAOK: " + AAOK + ", NAOK: " + NAOK + ", TightlyLinkingAATypeCount: " 
							+ vesicleGroup[i].tightlyLinkingAATypeCount + ", AACodonMatchOK: " + AACodonMatchOK + ").";
					}
					
					write "" + vesicleGroup[i] + " still has " + AACodonMatchFalseCount + " to synchronize.";
				}
			}			
		}
		
		if(LUCAAppear) {
			//
		}
		else {		
			vesicleGroup >>-  vesicleGroup select (each.alive = false);
			vesicleGroup <- shuffle(vesicleGroup);
		}		
	}
}

species FUCA {
	list<Peptide> peptideComb;
	list<RNA> RNAComb;
	bool alive <- true; 
	int tightlyLinkingAATypeCount;
	list<int> tightlyLinkingAAType <- list_with(20, 0); 
	int AACodonMatchFalseCount; 
		
	init {
		int peptideCount <- rnd(100, 150); 
		int RNACount <- rnd(300, 450); 	 
		create Peptide number: peptideCount returns: peptideGroup;
		create RNA number:RNACount returns: RNAGroup;
		peptideComb <- peptideGroup;
		RNAComb <- RNAGroup;				
	}
	
	action DoDie{
		do die;
	}
	
	action GetTightlyLinkingAATypeCount {		 
		loop i from: 0 to: length(peptideComb) - 1 {
			loop j from: 0 to: length(peptideComb[i].AAComb) - 1 {
				tightlyLinkingAAType[peptideComb[i].AAComb[j].type - 1] <- 1;	
			}			
		}		
		tightlyLinkingAATypeCount <- tightlyLinkingAAType count(each = 1);
	}
}

species Peptide {
	list<AA> AAComb;
	list<Codon> codonComb;
	int keySameCount;
	int accelerate <- rnd (0, 1); 
 	
	init {
		int AACount <- rnd(51, 100);	
		loop i from: 0 to: AACount - 1 {
			int rndIndexAA <- rnd(1, 20);
			switch rndIndexAA {
				match 1 {
					create AA1 number: 1 returns: AA1Group;
					add first(AA1Group) to: AAComb;
				}
				match 2 {
					create AA2 number: 1 returns: AA2Group;
					add first(AA2Group) to: AAComb;
				}
				match 3 {
					create AA3 number: 1 returns: AA3Group;
					add first(AA3Group) to: AAComb;
				}
				match 4 {
					create AA4 number: 1 returns: AA4Group;
					add first(AA4Group) to: AAComb;
				}
				match 5 {
					create AA5 number: 1 returns: AA5Group;
					add first(AA5Group) to: AAComb;
				}
				match 6 {
					create AA6 number: 1 returns: AA6Group;
					add first(AA6Group) to: AAComb;
				}
				match 7 {
					create AA7 number: 1 returns: AA7Group;
					add first(AA7Group) to: AAComb;
				}
				match 8 {
					create AA8 number: 1 returns: AA8Group;
					add first(AA8Group) to: AAComb;
				}
				match 9 {
					create AA9 number: 1 returns: AA9Group;
					add first(AA9Group) to: AAComb;
				}
				match 10 {
					create AA10 number: 1 returns: AA10Group;
					add first(AA10Group) to: AAComb;
				}
				match 11 {
					create AA11 number: 1 returns: AA11Group;
					add first(AA11Group) to: AAComb;
				}
				match 12 {
					create AA12 number: 1 returns: AA12Group;
					add first(AA12Group) to: AAComb;
				}
				match 13 {
					create AA13 number: 1 returns: AA13Group;
					add first(AA13Group) to: AAComb;
				}
				match 14 {
					create AA14 number: 1 returns: AA14Group;
					add first(AA14Group) to: AAComb;
				}
				match 15 {
					create AA15 number: 1 returns: AA15Group;
					add first(AA15Group) to: AAComb;
				}
				match 16 {
					create AA16 number: 1 returns: AA16Group;
					add first(AA16Group) to: AAComb;
				}
				match 17 {
					create AA17 number: 1 returns: AA17Group;
					add first(AA17Group) to: AAComb;
				}
				match 18 {
					create AA18 number: 1 returns: AA18Group;
					add first(AA18Group) to: AAComb;
				}
				match 19 {
					create AA19 number: 1 returns: AA19Group;
					add first(AA19Group) to: AAComb;
				}
				match 20 {
					create AA20 number: 1 returns: AA20Group;
					add first(AA20Group) to: AAComb;
				}
			}					
		}
		
		create Codon number: AACount returns: codonGroup;
		loop i over: codonGroup {
			add i to: codonComb;
		}		
	}
	
	action GetKeySameCount {
		keySameCount <- 0;
		loop i from: 0 to: length(AAComb) - 1 {
			if(AAComb[i].key = Codon[i].key) {
				keySameCount <- keySameCount + 1;
			}
		}
	}
}

species RNA {
	list<NA> NAComb;	
	
	init {
		int NACount <- rnd(151, 300); 
		loop i from: 0 to: NACount - 1 {
			int rndIndexNA <- rnd(1, 4);
			switch rndIndexNA {
				match 1 {
					create NA1 number: 1 returns: NA1Group;
					add first(NA1Group) to: NAComb;
				}
				match 2 {
					create NA2 number: 1 returns: NA2Group;
					add first(NA2Group) to: NAComb;
				}
				match 3 {
					create NA3 number: 1 returns: NA3Group;
					add first(NA3Group) to: NAComb;
				}
				match 4 {
					create NA4 number: 1 returns: NA4Group;
					add first(NA4Group) to: NAComb;
				}
			}	
		}		
	}
}

species AA {
	int key; 
	int type;  
}

species AA1 parent: AA {
	init {
		type <- 1;
		key <- rnd(0, 1);
	}
}

species AA2 parent: AA {
	init {
		type <- 2;
		key <- rnd(0, 1);
	}
}

species AA3 parent: AA {
	init {
		type <- 3;
		key <- rnd(0, 1);
	}
}

species AA4 parent: AA {
	init {
		type <- 4;
		key <- rnd(0, 1);
	}
}

species AA5 parent: AA {
	init {
		type <- 5;
		key <- rnd(0, 1);
	}
}

species AA6 parent: AA {
	init {
		type <- 6;
		key <- rnd(0, 1);
	}
}

species AA7 parent: AA {
	init {
		type <- 7;
		key <- rnd(0, 1);
	}
}

species AA8 parent: AA {
	init {
		type <- 8;
		key <- rnd(0, 1);
	}
}

species AA9 parent: AA {
	init {
		type <- 9;
		key <- rnd(0, 1);
	}
}

species AA10 parent: AA {
	init {
		type <- 10;
		key <- rnd(0, 1);
	}
}

species AA11 parent: AA {
	init {
		type <- 11;
		key <- rnd(0, 1);
	}
}

species AA12 parent: AA {
	init {
		type <- 12;
		key <- rnd(0, 1);
	}
}

species AA13 parent: AA {
	init {
		type <- 13;
		key <- rnd(0, 1);
	}
}

species AA14 parent: AA {
	init {
		type <- 14;
		key <- rnd(0, 1);
	}
}

species AA15 parent: AA {
	init {
		type <- 15;
		key <- rnd(0, 1);
	}
}

species AA16 parent: AA {
	init {
		type <- 16;
		key <- rnd(0, 1);
	}
}

species AA17 parent: AA {
	init {
		type <- 17;
		key <- rnd(0, 1);
	}
}

species AA18 parent: AA {
	init {
		type <- 18;
		key <- rnd(0, 1);
	}
}

species AA19 parent: AA {
	init {
		type <- 19;
		key <- rnd(0, 1);
	}
}

species AA20 parent: AA {
	init {
		type <- 20;
		key <- rnd(0, 1);
	}
}

species NA {
	int type;  //种类
}

species NA1 parent: NA {
	init {
		type <- 1;
	}
}

species NA2 parent: NA {
	init {
		type <- 2;
	}
}

species NA3 parent: NA {
	init {
		type <- 3;
	}
}

species NA4 parent: NA {
	init {
		type <- 4;
	}
}

species Codon {
	int key <- rnd (0, 1); 
}

experiment LUCA type: gui{
	output {
		display Chart type:java2D {
			chart "LUCA" type: series x_label: "Cycle(s)" y_label:"Number of AA-Codon to Be Synchronized" {
				data string(FUCA[0]) value: FUCA[0].AACodonMatchFalseCount color: #violet;
      	 		data string(FUCA[1]) value: FUCA[1].AACodonMatchFalseCount color: #green;
      	 		data string(FUCA[2]) value: FUCA[2].AACodonMatchFalseCount color: #red;
			}
		}
	}
}


