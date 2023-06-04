/**
* Name: FUCA_LUCA
* Author: GAO Ming
*/


model FUCA_LUCA

global {
	int tickCount <- 0; 
	Pond pond;
	float pContact; 
	list<Protocell> protocellGroup;  
	list<Protocell> FUCAGroup;  
	matrix<float> mergeP <- 0.0 as_matrix({6, 6});
	matrix<float> conjugateP <- 0.0 as_matrix({4 ,4});
	int numDeadProtocell;
	int numTotalProtocell;
	float pAbsorbMoreAANAOutside <- 5 * 10 ^ -3;
	bool isPause <- false;

	init {
		mergeP[0, 0] <- 0.5; mergeP[1, 0] <- 0.0; mergeP[2, 0] <- 0.0; mergeP[3, 0] <- 0.0; mergeP[4, 0] <- 0.0; mergeP[5, 0] <- 0.0;  
		mergeP[0, 1] <- 1.0; mergeP[1, 1] <- 0.5; mergeP[2, 1] <- 0.4; mergeP[3, 1] <- 0.3; mergeP[4, 1] <- 0.2; mergeP[5, 1] <- 0.1; 
		mergeP[0, 2] <- 1.0; mergeP[1, 2] <- 0.6; mergeP[2, 2] <- 0.5; mergeP[3, 2] <- 0.4; mergeP[4, 2] <- 0.3; mergeP[5, 2] <- 0.2; 
		mergeP[0, 3] <- 1.0; mergeP[1, 3] <- 0.7; mergeP[2, 3] <- 0.6; mergeP[3, 3] <- 0.5; mergeP[4, 3] <- 0.4; mergeP[5, 3] <- 0.3; 
		mergeP[0, 4] <- 1.0; mergeP[1, 4] <- 0.8; mergeP[2, 4] <- 0.7; mergeP[3, 4] <- 0.6; mergeP[4, 4] <- 0.5; mergeP[5, 4] <- 0.4;
		mergeP[0, 5] <- 1.0; mergeP[1, 5] <- 0.9; mergeP[2, 5] <- 0.8; mergeP[3, 5] <- 0.7; mergeP[4, 5] <- 0.6; mergeP[5, 5] <- 0.5;
		
		conjugateP[0, 0] <- 1 * 10 ^ -3;    conjugateP[1, 0] <- 0.6 * 10 ^ -3;   conjugateP[2, 0] <- 0.3 * 10 ^ -3;  conjugateP[3, 0] <- 0.1 * 10 ^ -3;  
		conjugateP[0, 1] <- 0.6 * 10 ^ -3;  conjugateP[1, 1] <- 0.3 * 10 ^ -3;   conjugateP[2, 1] <- 0.1 * 10 ^ -3;  conjugateP[3, 1] <- 0.05 * 10 ^ -3; 
		conjugateP[0, 2] <- 0.3 * 10 ^ -3;  conjugateP[1, 2] <- 0.1 * 10 ^ -3;   conjugateP[2, 2] <- 0.05 * 10 ^ -3; conjugateP[3, 2] <- 0.02 * 10 ^ -3; 
		conjugateP[0, 3] <- 0.1 * 10 ^ -3;  conjugateP[1, 3] <- 0.05 * 10 ^ -3;  conjugateP[2, 3] <- 0.02 * 10 ^ -3; conjugateP[3, 3] <- 0.01 * 10 ^ -3; 
		
		create Pond number: 1 returns: pondGroup;
		pond <- first(pondGroup);
		pond.currentVolume <- 100;
        
		write "Initiation is done!\n";
		write "===========================================";
	}
	
	reflex Ticks {
		if(cycle mod 2 = 0 ) {
			tickCount <- tickCount + 1;
			write "\n[Tick #" + tickCount + "]";
		}
	}
	
	reflex WetOrDry {
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
						float pPair1 <- mergeP[protocellGroup[k].type - 1, protocellGroup[i].type - 1];  
						float pPair2 <- mergeP[protocellGroup[i].type - 1, protocellGroup[k].type - 1];  
						if (flip(pPair1) or flip(pPair2)) { 
							list<Peptide> peptideCombSum <- [];
							list<RNA> RNACombSum <- [];
							peptideCombSum <<+ protocellGroup[i].peptideComb;
							peptideCombSum <<+ protocellGroup[k].peptideComb;
							RNACombSum <<+ protocellGroup[i].RNAComb;
							RNACombSum <<+ protocellGroup[k].RNAComb;
							int numPeptideSum <- length(peptideCombSum);
							int numRNASum <- length(RNACombSum);
							if (numPeptideSum >= 41 and numRNASum >= 120) { 
								create V6 number: 1 returns: newV6;
								numTotalProtocell <- numTotalProtocell + 1;
								first(newV6).peptideComb <- [];
								first(newV6).RNAComb <- [];
								first(newV6).peptideComb <<+ peptideCombSum;
								first(newV6).RNAComb <<+ RNACombSum;
								first(newV6).numPeptide <- numPeptideSum;
								first(newV6).numRNA <- numRNASum;		
								 
								int indexV6Flag <- 0;   
								int indexV6 <- -1;  
								if (protocellGroup[i].type = 6) {
									indexV6 <- i;
									indexV6Flag <- indexV6Flag + 1;
								}				
								if (protocellGroup[k].type = 6) {
									indexV6 <- k;
									indexV6Flag <- indexV6Flag + 1;
								}
								if (indexV6Flag = 1) { 
									loop m from: 0 to: 19 {										
										first(newV6).AAMatrix[m] <- protocellGroup[indexV6].AAMatrix[m]; 
										first(newV6).AACodon[m] <- protocellGroup[indexV6].AACodon[m];
									}
								}
								if (indexV6Flag = 2) {  
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
								}						

								int numSynchronized <- 0;
								loop m from: 0 to: 19 {
									numSynchronized <- numSynchronized + first(newV6).AACodon[m];
								}
								if (numSynchronized = 20) {
									first(newV6).isSynchronized <- true;
								}	
								add first(newV6) to: tempProtocellGroup;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV6);
								write "FUCA(" + first(newV6) + ") appears!" color: #red;
								protocellGroup[i].isAlive <- false;
								protocellGroup[k].isAlive <- false;								
							}
							else if (numPeptideSum >= 31 and numRNASum >= 75) { 
								create V5 number: 1 returns: newV5;
								numTotalProtocell <- numTotalProtocell + 1;
								first(newV5).peptideComb <- [];
								first(newV5).RNAComb <- [];
								first(newV5).peptideComb <<+ peptideCombSum;
								first(newV5).RNAComb <<+ RNACombSum;
								first(newV5).numPeptide <- numPeptideSum;
								first(newV5).numRNA <- numRNASum;								
								add first(newV5) to: tempProtocellGroup;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV5);
								protocellGroup[i].isAlive <- false;
								protocellGroup[k].isAlive <- false;								
							}
							else if (numPeptideSum >= 21 and numRNASum >= 50) { 
								create V4 number: 1 returns: newV4;
								numTotalProtocell <- numTotalProtocell + 1;
								first(newV4).peptideComb <- [];
								first(newV4).RNAComb <- [];
								first(newV4).peptideComb <<+ peptideCombSum;
								first(newV4).RNAComb <<+ RNACombSum;
								first(newV4).numPeptide <- numPeptideSum;
								first(newV4).numRNA <- numRNASum;								
								add first(newV4) to: tempProtocellGroup;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV4);
								protocellGroup[i].isAlive <- false;
								protocellGroup[k].isAlive <- false;								
							}
							else if (numPeptideSum >= 11 and numRNASum >= 25) { 
								create V3 number: 1 returns: newV3;
								numTotalProtocell <- numTotalProtocell + 1;
								first(newV3).peptideComb <- [];
								first(newV3).RNAComb <- [];
								first(newV3).peptideComb <<+ peptideCombSum;
								first(newV3).RNAComb <<+ RNACombSum;
								first(newV3).numPeptide <- numPeptideSum;
								first(newV3).numRNA <- numRNASum;								
								add first(newV3) to: tempProtocellGroup;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV3);
								protocellGroup[i].isAlive <- false;
								protocellGroup[k].isAlive <- false;								
							}
							else if (numPeptideSum >= 5 and numRNASum >= 10) { 
								create V2 number: 1 returns: newV2 {
									do init;
								}
								numTotalProtocell <- numTotalProtocell + 1;
								first(newV2).peptideComb <- [];
								first(newV2).RNAComb <- [];
								first(newV2).peptideComb <<+ peptideCombSum;
								first(newV2).RNAComb <<+ RNACombSum;
								first(newV2).numPeptide <- numPeptideSum;
								first(newV2).numRNA <- numRNASum;								
								add first(newV2) to: tempProtocellGroup;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV2);
								protocellGroup[i].isAlive <- false;
								protocellGroup[k].isAlive <- false;								
							}
							else { 
								create V1 number: 1 returns: newV1 {
									do init;
								}
								numTotalProtocell <- numTotalProtocell + 1;
								first(newV1).peptideComb <- [];
								first(newV1).RNAComb <- [];
								first(newV1).peptideComb <<+ peptideCombSum;
								first(newV1).RNAComb <<+ RNACombSum;
								first(newV1).numPeptide <- numPeptideSum;
								first(newV1).numRNA <- numRNASum;								
								add first(newV1) to: tempProtocellGroup;
								write string(protocellGroup[i]) + " and " + protocellGroup[k] + " are merged as " + first(newV1);
								protocellGroup[i].isAlive <- false;
								protocellGroup[k].isAlive <- false;								
							}
						}
					}
				}
			}
		}
		
		int naturaldeath;
		loop i over: protocellGroup select (each.isAlive = true) {
			ask i {
				do GetFitnessScore;
				
				if (!flip(i.pSurvival)) {
					i.isAlive <- false;
					naturaldeath <- naturaldeath + 1;
				}
			}
		}
		write "Natural death: " + naturaldeath;
		 
		numDeadProtocell <- numDeadProtocell + length(protocellGroup select (each.isAlive = false));	
		
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
		
		list<Protocell> tempMoreProtocellGroup <- [];  
		
		loop i over: tempProtocellGroup { 
			list<Peptide> tempPeptideGroup <- [];  
			
			loop j from: 0 to: length(i.peptideComb) - 1 {
				int k <- j + 1;
				if (k <= length(i.peptideComb) - 1) {
					loop m from: k to: length(i.peptideComb) - 1 {
						if (i.peptideComb[j].isAlive and i.peptideComb[m].isAlive) { 
							int jLength <- length(i.peptideComb[j].AAComb);
							int mLength <- length(i.peptideComb[m].AAComb);
							int jIndex;
							int mIndex;
							
							if (jLength <= 10) {
								jIndex <- 0;
							}
							else if (jLength >= 11 and jLength <= 25) {
								jIndex <- 1;
							}
							else if (jLength >= 26 and jLength <= 50) {
								jIndex <- 2;
							}
							else {
								jIndex <- 3;
							}
							
							if (mLength <= 10) {
								mIndex <- 0;
							}
							else if (mLength >= 11 and mLength <= 25) {
								mIndex <- 1;
							}
							else if (mLength >= 26 and mLength <= 50) {
								mIndex <- 2;
							}
							else {
								mIndex <- 3;
							}
							
							float pConjugate <- conjugateP[mIndex, jIndex];
							
							if(flip(pConjugate)) {
								list<AA> AACombSum <- [];
								AACombSum <<+ i.peptideComb[j].AAComb;
								AACombSum <<+ i.peptideComb[m].AAComb;
								int numAACombSum <- length(AACombSum);								
								create Peptide number: 1 returns: newPeptide;
								first(newPeptide).numAA <- numAACombSum;
								first(newPeptide).AAComb <- [];
								first(newPeptide).AAComb <<+ AACombSum;
								add first(newPeptide) to: tempPeptideGroup;
								write string(i.peptideComb[j]) + " and " + i.peptideComb[m] + " are conjugated as " + first(newPeptide);
								i.peptideComb[j].isAlive <- false;
								i.peptideComb[m].isAlive <- false;
								i.isConjugated <- true;
							}	
						}												
					}
				}				
			}
			
			if (length(i.peptideComb select (each.isAlive = false)) > 0) {	
				list<Peptide> readyToKillPeptide <- [];   
				readyToKillPeptide <<+ i.peptideComb select (each.isAlive = false);
				i.peptideComb >>- i.peptideComb select (each.isAlive = false);	
				loop p over: readyToKillPeptide {
					ask p {
						do Kill;
					}
				}
				i.peptideComb <<+ tempPeptideGroup;	
				i.numPeptide <- length(i.peptideComb);				
				i.isConjugated <- true;
			}			
		 
			list<RNA> tempRNAGroup <- [];  
			
			loop j from: 0 to: length(i.RNAComb) - 1 {
				int k <- j + 1;
				if (k <= length(i.RNAComb) - 1) {
					loop m from: k to: length(i.RNAComb) - 1 {
						if (i.RNAComb[j].isAlive and i.RNAComb[m].isAlive) {   
							int jLength <- i.RNAComb[j].numNBase;
							int mLength <- i.RNAComb[m].numNBase;
							int jIndex;
							int mIndex;

							if (jLength <= 30) {
								jIndex <- 0;
							}
							else if (jLength >= 31 and jLength <= 60) {
								jIndex <- 1;
							}
							else if (jLength >= 61 and jLength <= 100) {
								jIndex <- 2;
							}
							else {
								jIndex <- 3;
							}
							
							if (mLength <= 30) {
								mIndex <- 0;
							}
							else if (mLength >= 31 and mLength <= 60) {
								mIndex <- 1;
							}
							else if (mLength >= 61 and mLength <= 100) {
								mIndex <- 2;
							}
							else {
								mIndex <- 3;
							}

							float pConjugate <- conjugateP[mIndex, jIndex];

							if(flip(pConjugate)) {
								int numNBaseSum <- 0;
								numNBaseSum <- numNBaseSum + i.RNAComb[j].numNBase;
								numNBaseSum <- numNBaseSum + i.RNAComb[m].numNBase;
								create RNA number: 1 returns: newRNA;
								first(newRNA).numNBase <- numNBaseSum;
								add first(newRNA) to: tempRNAGroup;
								write string(i.RNAComb[j]) + " and " + i.RNAComb[m] + " are conjugated as " + first(newRNA);
								i.RNAComb[j].isAlive <- false;
								i.RNAComb[m].isAlive <- false;
								i.isConjugated <- true;
							}
						}
					}
				}
			}

			if (length(i.RNAComb select (each.isAlive = false)) > 0) {
				list<RNA> readyToKillRNA <- [];  
				readyToKillRNA <<+ i.RNAComb select (each.isAlive = false);				
				i.RNAComb >>- i.RNAComb select (each.isAlive = false);
				loop p over: readyToKillRNA {
					ask p {
						do Kill;
					}
				}
				i.RNAComb <<+ tempRNAGroup;
				i.numRNA <- length(i.RNAComb);
				i.isConjugated <- true;
			}			
			 
			if (i.isConjugated) {
				int numPeptide <- length(i.peptideComb);
				int numRNA <- length(i.RNAComb);

				if (((numPeptide >= 31 and numPeptide <= 40) or (numRNA >= 75 and numRNA < 120)) and i.type > 5) { 								
					create V5 number: 1 returns: newV5;
					first(newV5).peptideComb <- [];
					first(newV5).RNAComb <- [];
					first(newV5).peptideComb <<+ i.peptideComb;  
					first(newV5).RNAComb <<+ i.RNAComb;  
					first(newV5).numPeptide <- length(i.peptideComb);
					first(newV5).numRNA <- length(i.RNAComb);
					add first(newV5) to: tempMoreProtocellGroup;
					write string(i) + " has changed into " + first(newV5);
					i.isAlive <- false;								
				}
				else if (((numPeptide >= 21 and numPeptide <= 30) or (numRNA >= 50 and numRNA < 75)) and i.type > 4) { 								
					create V4 number: 1 returns: newV4;
					first(newV4).peptideComb <- [];
					first(newV4).RNAComb <- [];
					first(newV4).peptideComb <<+ i.peptideComb;   
					first(newV4).RNAComb <<+ i.RNAComb;  
					first(newV4).numPeptide <- length(i.peptideComb);
					first(newV4).numRNA <- length(i.RNAComb);	
					add first(newV4) to: tempMoreProtocellGroup;
					write string(i) + " has changed into " + first(newV4);
					i.isAlive <- false;							
				}
				else if (((numPeptide >= 11 and numPeptide <= 20) or (numRNA >= 25 and numRNA < 50)) and i.type > 3) { 								
					create V3 number: 1 returns: newV3;
					first(newV3).peptideComb <- [];
					first(newV3).RNAComb <- [];
					first(newV3).peptideComb <<+ i.peptideComb;   
					first(newV3).RNAComb <<+ i.RNAComb; 
					first(newV3).numPeptide <- length(i.peptideComb);
					first(newV3).numRNA <- length(i.RNAComb);	
					add first(newV3) to: tempMoreProtocellGroup;
					write string(i) + " has changed into " + first(newV3);
					i.isAlive <- false;							
				}
				else if (((numPeptide >= 5 and numPeptide <= 10) or (numRNA >= 10 and numRNA < 25)) and i.type > 2) { 								
					create V2 number: 1 returns: newV2 {
						type <- 2;						
					}
					first(newV2).peptideComb <- [];
					first(newV2).RNAComb <- [];
					first(newV2).peptideComb <<+ i.peptideComb;   
					first(newV2).RNAComb <<+ i.RNAComb;  
					first(newV2).numPeptide <- length(i.peptideComb);
					first(newV2).numRNA <- length(i.RNAComb);	
					add first(newV2) to: tempMoreProtocellGroup;
					write string(i) + " has changed into " + first(newV2);
					i.isAlive <- false;								
				}
				else if ((numPeptide <= 4 or numRNA < 10) and i.type > 1) { 	
					create V1 number: 1 returns: newV1 {
						type <- 1;
					}
					first(newV1).peptideComb <- [];
					first(newV1).RNAComb <- [];
					first(newV1).peptideComb <<+ i.peptideComb;   
					first(newV1).RNAComb <<+ i.RNAComb; 
					first(newV1).numPeptide <- length(i.peptideComb);
					first(newV1).numRNA <- length(i.RNAComb);
					add first(newV1) to: tempMoreProtocellGroup;
					write string(i) + " has changed into " + first(newV1);
					i.isAlive <- false;								
				}
			}
		
		}
		
		if (length(tempProtocellGroup select (each.isAlive = false)) > 0) {			
			list<Protocell> readyToKillGroup2 <- [];   
			readyToKillGroup2 <<+ tempProtocellGroup select (each.isAlive = false);	
			tempProtocellGroup >>- tempProtocellGroup select (each.isAlive = false);
			loop p over: readyToKillGroup2 {
				ask p {
					do Kill;
				}
			}
			protocellGroup <<+ tempProtocellGroup;
			protocellGroup <<+ tempMoreProtocellGroup;		
		}
		else {
			protocellGroup <<+ tempProtocellGroup;
		}
	}
	  
	reflex AbsorbMore when: cycle mod 2 != 0 {
		loop i over: protocellGroup {
			ask i {
				if (flip(pAbsorbMoreAANAOutside)) {
					do AbsorbMoreAANAOutside;
				}
			}
		}
	}
	
	reflex GetOtherAAs when: length(V6) > 0 and (cycle mod 2 != 0) {
		loop i over: V6 {
			ask i {
				do CheckGetFirstTenAAs;
				
				if (isGetFirstTenAAs) {
					do CheckGetAllAAs;
					do SynthesiseLastTenAAs;						
				}
				
				if (isGetAllAAs) {  
					do SynthesiseMoreLastTenAAs;
				}							
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

	reflex CheckLUCA when: length(V6) > 0 {
		loop i over: V6 {
			int numPeptide;
			int numRNA;			
			bool OKAA;
			bool OKNBase;
			
			ask i {
				if (isGetAllAAs and isSynchronized) {
					numPeptide <- length(i.peptideComb);
					numRNA <- length(i.RNAComb);
					write "Check LUCA";
					write "Peptides: " + numPeptide;
					write "RNAs: " + numRNA;
							
					if(numPeptide >= 100 and numRNA >= 300) {
						int numOKAAPeptide <- 0;
						int numTotalAAPeptide <- 0;
						int numOKNBaseRNA <- 0;
						int numTotalNBaseRNA <- 0;
			
						loop j over: i.peptideComb {
							write "" + j + ".AA:" + length(j.AAComb);
							
							if (length(j.AAComb) >= 50) {
								numOKAAPeptide <- numOKAAPeptide + 1;
							}
							
							numTotalAAPeptide <- numTotalAAPeptide + 1;
						}
						
						if (numOKAAPeptide / numTotalAAPeptide >= 0.8) {
							OKAA <- true;
						}
						else {
							OKAA <- false;
						}
						
						loop k over: i.RNAComb {
							write "" + k + ".NBase:" + k.numNBase;
							
							if (k.numNBase >= 150) {
								numOKNBaseRNA <- numOKNBaseRNA + 1;
							}
							
							numTotalNBaseRNA <- numTotalNBaseRNA + 1;
						}
						
						if (numOKNBaseRNA / numTotalNBaseRNA >= 0.8) {
							OKNBase <- true;
						}
						else {
							OKNBase <- false;
						}
						
						if(OKAA and OKNBase) {							
							i.isLUCA <- true;
							isPause <- true;
							write "LUCA(" + i + ") appears!" color: #red;
							write "Peptides: " + numPeptide;
							write "RNAs: " + numRNA;
							loop j over: i.peptideComb {
								write "" + j + ".AA:" + length(j.AAComb);
							}
							loop k over: i.RNAComb {
								write "" + k + ".NBase:" + k.numNBase;
							}													
						}
					}
				}
			}
			
			if(isPause) {
				do pause;
			}
		}
	}  
}

species Pond {
	bool isWet <- false;
	int currentVolume;
	int V1InitNum <- 3000;
	int V2InitNum <- 2000;
	float V1Ratio <- V1InitNum / (V1InitNum + V2InitNum);
	
	init {
		create V1 number: V1InitNum returns: V1Group {
			do init;
		}
		create V2 number: V2InitNum returns: V2Group {
			do init;	
		}
		numTotalProtocell <- V1InitNum + V2InitNum;
		protocellGroup <<+ V1Group;
		protocellGroup <<+ V2Group;		
	}
	
	action GetDryPhaseVolume {
		currentVolume <- rnd(50, 80);
		pContact <- rnd(10 ^ -5, 10 ^ -6) * ((100.0 / currentVolume) ^ 2); 
		do GetNewV1V2_Dry;			
	}
	
	action GetWetPhaseVolume {
		currentVolume <- rnd(80, 100);
		pContact <- rnd(10 ^ -5, 10 ^ -6) * ((100.0 / currentVolume) ^ 2);  
		do GetNewV1V2_Wet;
	}
	
	action GetNewV1V2_Dry {
		int totalNum <- rnd(2000, 3000);
		create V1 number: int(totalNum * V1Ratio) returns: newV1Group {
			do init;	
		}
		create V2 number: totalNum - int(totalNum * V1Ratio) returns: newV2Group {
			do init;
		}
		numTotalProtocell <- numTotalProtocell + totalNum;
		protocellGroup <<+ newV1Group;
		protocellGroup <<+ newV2Group;	
	}
	
	action GetNewV1V2_Wet {
		int totalNum <- rnd(4000, 5000);			
		create V1 number: int(totalNum * V1Ratio) returns: newV1Group {
			do init;
		}
		create V2 number: totalNum - int(totalNum * V1Ratio) returns: newV2Group {
			do init;
		}
		numTotalProtocell <- numTotalProtocell + totalNum;
		protocellGroup <<+ newV1Group;
		protocellGroup <<+ newV2Group;	
	}
}

species AA {
	int type;  
}

species AA1 parent: AA {
	init {
		type <- 1;
	}
}

species AA2 parent: AA {
	init {
		type <- 2;
	}
}

species AA3 parent: AA {
	init {
		type <- 3;
	}
}

species AA4 parent: AA {
	init {
		type <- 4;
	}
}

species AA5 parent: AA {
	init {
		type <- 5;
	}
}

species AA6 parent: AA {
	init {
		type <- 6;
	}
}

species AA7 parent: AA {
	init {
		type <- 7;
	}
}

species AA8 parent: AA {
	init {
		type <- 8;
	}
}

species AA9 parent: AA {
	init {
		type <- 9;
	}
}

species AA10 parent: AA {
	init {
		type <- 10;
	}
}

species AA11 parent: AA {
	init {
		type <- 11;
	}
}

species AA12 parent: AA {
	init {
		type <- 12;
	}
}

species AA13 parent: AA {
	init {
		type <- 13;
	}
}

species AA14 parent: AA {
	init {
		type <- 14;
	}
}

species AA15 parent: AA {
	init {
		type <- 15;
	}
}

species AA16 parent: AA {
	init {
		type <- 16;
	}
}

species AA17 parent: AA {
	init {
		type <- 17;
	}
}

species AA18 parent: AA {
	init {
		type <- 18;
	}
}

species AA19 parent: AA {
	init {
		type <- 19;
	}
}

species AA20 parent: AA {
	init {
		type <- 20;
	}
}

species Peptide {
	int numAA;
	list<AA> AAComb;
	bool isAlive <- true;
 	
	action init {
		loop i from: 0 to: numAA - 1 {
			int rndIndexAA <- rnd(1, 10);  
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
			}					
		}
	}

	action Kill {
		do die;
	}
}

species RNA {
	int numNBase;
	bool isAlive <- true;

	action Kill {
		do die;
	}
}

species Protocell {
	int type; 
	int numPeptide;
	int numRNA;
	list<Peptide> peptideComb;
	list<RNA> RNAComb;
	float pSurvival;	
	bool isAlive <- true;
	matrix<int> AAMatrix <- 0 as_matrix({20, 1});
	bool isGetFirstTenAAs <- false;   
	int numFirstTenAAs <- 0; 
	bool isGetAllAAs <- false;  
	bool isConjugated <- false;  
	matrix<int> AACodon <- 0 as_matrix({20, 1});  
	bool isSynchronized <- false;  
	bool isLUCA <- false; 
	
	action GetFitnessScore {
		int AACount <- 0;
		int fitnessScore <- 0;
		
		loop i over: peptideComb {
			AACount <- AACount + length(i.AAComb);
		}
		
		fitnessScore <- AACount * length(peptideComb);
		pSurvival <- ln(fitnessScore) / 10;
	}
	
	action CheckGetFirstTenAAs {
		loop i over: peptideComb {
			loop j over: i.AAComb {
				if (j.type <= 10) {
					numFirstTenAAs <- numFirstTenAAs + 1;
					AAMatrix[j.type - 1] <- 1;
				}				
			}
		}
		
		if (!isGetFirstTenAAs) {
			int tmpSum;
			
			loop i from: 0 to: 9 {
				tmpSum <- tmpSum + AAMatrix[i];
			}
			
			if (tmpSum = 10) {
				isGetFirstTenAAs <- true;
				write "" + self + " gets first ten AAs.";
			}
			else {
				isGetFirstTenAAs <- false;
				write "" + self + " doesn't get first ten AAs.";
			}
		}
	}	
	
	action SynthesiseLastTenAAs {
		list<AA> lastTenGroup;
		bool isTrue <- true;  
		 
		if (!isGetAllAAs) {  
			loop while: isTrue {
				int rndIndexAA <- rnd(11, 20);
				 				
				if(AAMatrix[rndIndexAA - 1] = 0) {
					switch rndIndexAA {
						match 11 {
							create AA11 number: rnd(20, numFirstTenAAs / 20) returns: AA11Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA11Group;
							isTrue <- false;
						}
						match 12 {
							create AA12 number: rnd(20, numFirstTenAAs / 20) returns: AA12Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA12Group;
							isTrue <- false;
						}
						match 13 {
							create AA13 number: rnd(20, numFirstTenAAs / 20) returns: AA13Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA13Group;
							isTrue <- false;
						}
						match 14 {
							create AA14 number: rnd(20, numFirstTenAAs / 20) returns: AA14Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA14Group;
							isTrue <- false;
						}
						match 15 {
							create AA15 number: rnd(20, numFirstTenAAs / 20) returns: AA15Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA15Group;
							isTrue <- false;
						}
						match 16 {
							create AA16 number: rnd(20, numFirstTenAAs / 20) returns: AA16Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA16Group;
							isTrue <- false;
						}
						match 17 {
							create AA17 number: rnd(20, numFirstTenAAs / 20) returns: AA17Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA17Group;
							isTrue <- false;
						}
						match 18 {
							create AA18 number: rnd(20, numFirstTenAAs / 20) returns: AA18Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA18Group;
							isTrue <- false;
						}
						match 19 {
							create AA19 number: rnd(20, numFirstTenAAs / 20) returns: AA19Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA19Group;
							isTrue <- false;
						}
						match 20 {
							create AA20 number: rnd(20, numFirstTenAAs / 20) returns: AA20Group;
							AAMatrix[rndIndexAA - 1] <- 1;
							lastTenGroup <<+ AA20Group;
							isTrue <- false;
						}						
					}			
				}
			} 
			
			int numMoreNBase <- 3 * length(lastTenGroup);  
			
			loop while: length(lastTenGroup) > 0 { 
				loop i over: peptideComb {
					if (flip(rnd(0.1))) {   
						if (length(lastTenGroup) > 0) {
							int numAA <- rnd(1, int(length(lastTenGroup) / 20));    
 							loop j from: 1 to: numAA {
								i.AAComb <+ first(lastTenGroup);
								lastTenGroup >- first(lastTenGroup);							
							}						
							write "" + self + "." + i + " synthesises " + numAA + " AA" + last(i.AAComb).type + "(s).";
						}
						else {
							break;
						}	
					}											
				}
			}
			
			loop while: numMoreNBase > 0 {  
				loop i over: RNAComb {
					if (flip(rnd(0.03))) {  
						if (numMoreNBase > 0) {
							int numNN <- rnd(3, int(numMoreNBase / 30));   
							loop while: (numNN mod 3 != 0) {  
								numNN <- rnd(3, numMoreNBase);
							}
							i.numNBase <- i.numNBase + numNN;
							numMoreNBase <- numMoreNBase - numNN;
							write "" + self + "." + i + " gets " + numNN + " NN(s).";
						}
						else {
							break;
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
	
	action SynthesiseMoreLastTenAAs {
		list<AA> lastTenGroup;
		 
		int rndIndexAA <- rnd(11, 20);
			 				
		switch rndIndexAA {
			match 11 {
				create AA11 number: rnd(30, numFirstTenAAs / 20) returns: AA11Group;
				lastTenGroup <<+ AA11Group;					
			}
			match 12 {
				create AA12 number: rnd(30, numFirstTenAAs / 20) returns: AA12Group;
				lastTenGroup <<+ AA12Group;
			}
			match 13 {
				create AA13 number: rnd(30, numFirstTenAAs / 20) returns: AA13Group;
				lastTenGroup <<+ AA13Group;
			}
			match 14 {
				create AA14 number: rnd(30, numFirstTenAAs / 20) returns: AA14Group;
				lastTenGroup <<+ AA14Group;
			}
			match 15 {
				create AA15 number: rnd(30, numFirstTenAAs / 20) returns: AA15Group;
				lastTenGroup <<+ AA15Group;
			}
			match 16 {
				create AA16 number: rnd(30, numFirstTenAAs / 20) returns: AA16Group;
				lastTenGroup <<+ AA16Group;
			}
			match 17 {
				create AA17 number: rnd(30, numFirstTenAAs / 20) returns: AA17Group;
				lastTenGroup <<+ AA17Group;
			}
			match 18 {
				create AA18 number: rnd(30, numFirstTenAAs / 20) returns: AA18Group;
				lastTenGroup <<+ AA18Group;
			}
			match 19 {
				create AA19 number: rnd(30, numFirstTenAAs / 20) returns: AA19Group;
				lastTenGroup <<+ AA19Group;
			}
			match 20 {
				create AA20 number: rnd(30, numFirstTenAAs / 20) returns: AA20Group;
				lastTenGroup <<+ AA20Group;
			}						
		}
		
		int numMoreNBase <- 3 * length(lastTenGroup); 
		
		loop while: length(lastTenGroup) > 0 { 
			loop i over: peptideComb {
				if (flip(rnd(0.1))) {  
					if (length(lastTenGroup) > 0) {
						int numAA <- rnd(1, int(length(lastTenGroup) / 20));   
						loop j from: 1 to: numAA {
							i.AAComb <+ first(lastTenGroup);
							lastTenGroup >- first(lastTenGroup);							
						}						
						write "" + self + "." + i + " synthesises " + numAA + " AA" + last(i.AAComb).type + "(s).";
					}
					else {
						break;
					}	
				}											
			}
		}
		
		loop while: numMoreNBase > 0 {  
			loop i over: RNAComb {
				if (flip(rnd(0.03))) {  
					if (numMoreNBase > 0) {
						int numNN <- rnd(3, int(numMoreNBase / 30));   
						loop while: (numNN mod 3 != 0) {  
							numNN <- rnd(3, numMoreNBase);
						}
						i.numNBase <- i.numNBase + numNN;
						numMoreNBase <- numMoreNBase - numNN;
						write "" + self + "." + i + " gets " + numNN + " NN(s).";
					}
					else {
						break;
					}	
				}					
			}				
		}	
	}
	
	action AbsorbMoreAANAOutside {
		list<AA> firstTenGroup;
		 
		int rndIndexAA <- rnd(1, 10);
			 				
		switch rndIndexAA {
			match 1 {
				create AA1 number: rnd(30, 60) returns: AA1Group;
				firstTenGroup <<+ AA1Group;					
			}
			match 2 {
				create AA2 number: rnd(30, 60) returns: AA2Group;
				firstTenGroup <<+ AA2Group;	
			}
			match 3 {
				create AA3 number: rnd(30, 60) returns: AA3Group;
				firstTenGroup <<+ AA3Group;	
			}
			match 4 {
				create AA4 number: rnd(30, 60) returns: AA4Group;
				firstTenGroup <<+ AA4Group;	
			}
			match 5 {
				create AA5 number: rnd(30, 60) returns: AA5Group;
				firstTenGroup <<+ AA5Group;	
			}
			match 6 {
				create AA6 number: rnd(30, 60) returns: AA6Group;
				firstTenGroup <<+ AA6Group;	
			}
			match 7 {
				create AA7 number: rnd(30, 60) returns: AA7Group;
				firstTenGroup <<+ AA7Group;	
			}
			match 8 {
				create AA8 number: rnd(30, 60) returns: AA8Group;
				firstTenGroup <<+ AA8Group;	
			}
			match 9 {
				create AA9 number: rnd(30, 60) returns: AA9Group;
				firstTenGroup <<+ AA9Group;	
			}
			match 10 {
				create AA10 number: rnd(30, 60) returns: AA10Group;
				firstTenGroup <<+ AA10Group;	
			}						
		}
		
		int numMoreNBase <- 3 * length(firstTenGroup);  
		
		loop while: length(firstTenGroup) > 0 { 
			loop i over: peptideComb {
				if (flip(rnd(0.1))) {  
					if (length(firstTenGroup) > 0) {
						int numAA <- rnd(1, int(length(firstTenGroup) / 20)); 
						loop j from: 1 to: numAA {
							i.AAComb <+ first(firstTenGroup);
							firstTenGroup >- first(firstTenGroup);							
						}						
						write "" + self + "." + i + " absorbs " + numAA + " AA" + last(i.AAComb).type + "(s).";
					}
					else {
						break;
					}	
				}											
			}
		}
		
		loop while: numMoreNBase > 0 {  
			loop i over: RNAComb {
				if (flip(rnd(0.03))) {   
					if (numMoreNBase > 0) {
						int numNN <- rnd(3, int(numMoreNBase / 30));   
						loop while: (numNN mod 3 != 0) {  
							numNN <- rnd(3, numMoreNBase);
						}
						i.numNBase <- i.numNBase + numNN;
						numMoreNBase <- numMoreNBase - numNN;
						write "" + self + "." + i + " absorbs " + numNN + " NN(s).";
					}
					else {
						break;
					}	
				}					
			}				
		}
	}
	
	action SynchronizeAACodon {
		bool isTrue <- true;  
		int numSynchronized <- 0; 
		float part;
		float pSynchronize;
		loop i over: AACodon {
			write "AACodon:" + i;
			numSynchronized <- numSynchronized + i;
		}		
		//part <- 20 / (20 - numSynchronized);
		//pSynchronize <- (10 ^ -3) * (part ^ 2);
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

species V1 parent: Protocell {
	action init {
		type <- 1;
		numPeptide <- rnd(2, 4);
		numRNA <- rnd(5, 10);
		
		create Peptide number: numPeptide returns: peptideGroup {
			numAA <- rnd(5, 6);
			do init;
		}
		peptideComb <<+ peptideGroup;
		
		create RNA number: numRNA returns: RNAGroup {
			numNBase <- 3 * rnd(5, 6);
		}
		RNAComb <<+ RNAGroup;		
	}
}

species V2 parent: Protocell {
	action init {
		type <- 2;
		numPeptide <- rnd(5, 10);
		numRNA <- rnd(10, 30);
		
		create Peptide number: numPeptide returns: peptideGroup {
			numAA <- rnd(6, 8);
			do init;
		}
		peptideComb <<+ peptideGroup;
		
		create RNA number: numRNA returns: RNAGroup {
			numNBase <- 3 * rnd(6, 8);
		}
		RNAComb <<+ RNAGroup;		
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

experiment LUCA_FUCA type: gui{
	output {
		display Chart type:java2D {
			chart "FUCA_LUCA" type: series x_label: "Cycle(s)"  y_label: "Vesicle Number" {
 				data "V1" value: length(V1) color: #violet;
 				data "V2" value: length(V2) color: #green;
 				data "V3" value: length(V3) color: #orange;
 				data "V4" value: length(V4) color: #lightgreen;
 				data "V5" value: length(V5) color: #blue;
 				data "V6(FUCA)" value: length(V6 select (each.isLUCA = false)) color: #red;
			}
		}
				
		monitor "Volumn of Pond" value: pond.currentVolume;
		monitor "Number of V1" value: length(V1);
		monitor "Number of V2" value: length(V2);
		monitor "Number of V3" value: length(V3);
		monitor "Number of V4" value: length(V4);
		monitor "Number of V5" value: length(V5);
		monitor "Number of V6(FUCA)" value: length(V6 select (each.isLUCA = false));
		monitor "Number of LUCA" value: length(V6 select (each.isLUCA = true));
		monitor "Number of Alive Protocells" value: length(V1) + length(V2) + length(V3) + length(V4) + length(V5) + length(V6);
		monitor "Number of Dead Protocells" value: numDeadProtocell;
		monitor "Number of Total Protocells" value: numTotalProtocell;
	}
}

