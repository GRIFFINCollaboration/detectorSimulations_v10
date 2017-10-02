#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <cctype> //isdigit etc

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "BeamDistribution.hh"//single-file now

BeamDistribution::BeamDistribution() {
  fDigitBool = false;
}
BeamDistribution::~BeamDistribution() {}

void BeamDistribution::ReadIn(G4String filename){//Reads in the data from the file
    
    std::ifstream distrodata;//create input object
    distrodata.open(filename.c_str());//open data file
    
    if(distrodata.is_open()) {//opens as gets in here
	G4cout<<"\n---> " "Reading the distribution from "<<filename<<" ... "<<std::endl;
	
	char c; int counter;
	while(fDigitBool == false){
	  c = distrodata.get();
	  if ( !std::isdigit(c) ) {
	      std::getline(distrodata,fInfo);
	      counter++;
	  }else fDigitBool = true;
	}
	
	for(fCol = 0; std::getline(distrodata, fInfo); ++fCol){}
	std::cout << "Number of data points in text file: " << fCol << std::endl;
	
	distrodata.clear();
	distrodata.seekg(0,std::ios::beg);//reset istream to beginning of file
	

	for(fIter = 0; fIter<counter; fIter++){//skip text lines
	  std::getline(distrodata,fInfo);
	}
	
	std::cout << "Reading values ... "<< std::endl;
	fSpiceTargetThickness *= 1000.; fUsedLines = 0;
	for(int i=0; i<fCol; i++){
// 	  G4cout << i << std::endl;
	  distrodata >> fDist >> fRelProb ;
	  if(fDist <= fSpiceTargetThickness){//only for SPICE at the moment, could use booleans if want to apply more generally
	    fDistVec.push_back(fDist);
	    fRelProbVec.push_back(fRelProb);
	    fUsedLines++; //Used for iterator limits elsewhere
	  }
	}
	distrodata.close();
    }else std::cout << "File not opened" << std::endl; 
}

void BeamDistribution::ReadOut(){//reads out the data, mainly used for debugging. Called in PGA if need to turn off
    std::cout <<  "CALCIUM DATA: " <<  std::endl;
    fIter = 0;
    while(fIter < fUsedLines){
      std::cout << "Stopping Distance (um): " << fDistVec[fIter] << "	Relative Probability: " << fRelProbVec[fIter] << std::endl;
      fIter++;
    }
}

void BeamDistribution::SumProb(){
    G4double sum = 0.;
    for(fIter = 0; fIter < fUsedLines; fIter++){
      sum += fRelProbVec[fIter];//feed sum value (if sum=1 is relative, else divide by sum to normalise)
    }
    std::cout << "Total inputted \"Probability\": " << sum << std::endl;
    if(sum != 1.){
      std::cout << "Normalising probability" << std::endl << std::endl;
      AdjustProb(sum);
    }
}

void BeamDistribution::AdjustProb(double sum){//runs through the probability vector and normalises
    for(fIter = 0; fIter < fUsedLines; fIter++){
	fRelProbVec[fIter] /= sum;
    }
}

double BeamDistribution::SelectDist(G4double RelProb){//the stopping distance for the beam is chosen here on an event by event basis using a random number generator to specify the probabilty
    G4double counter = 0;				//using a random number generator to specify the probabilty value chosen
    for(fIter = 0; fIter < fUsedLines; fIter++){
      counter += fRelProbVec[fIter];
      if(counter>RelProb) {
	return fDistVec[fIter]*(micrometer);//sends into PGA
	break;
      }
    }
    std::cout << "SOMETHING WENT WRONG" << std::endl;
    return 0.;//EXIT criteria if failed
}

//Bismuth PGA
/*
	quick bismuth energy peak func
	CE K 481.7 1.537 %
	  CE L   553.8 0.442 %
	  CE M  565.9 0.111 %
	  CE K 975.7   7.08 %
	  CE L    1047.8  1.84 %
	  CE M 1059.8  0.44 %
	G4double BiEn[6] =  {481.7,553.8,565.9,975.7,1047.8,1059.8};
	G4double EnSplit[6] = {0.1342358079,0.1728384279,0.1825327511,0.8008733624,0.9615720524,1.};
	G4double BiRand = G4UniformRand();
	G4cout << "Rando" << BiRand << G4endl;
	G4int BiIter;
	for(BiIter = 0; BiIter <6; BiIter ++){
	 // G4cout << "Iterator" << BiIter << G4endl;
	  if(BiRand < EnSplit[BiIter] ) {
	    fEffEnergy = BiEn[BiIter]*keV;
	    break;
	  } 
	}
	
	*/