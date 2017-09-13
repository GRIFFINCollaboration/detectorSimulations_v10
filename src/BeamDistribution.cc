#include <iostream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "BeamDistribution.hh"//single-file now

BeamDistribution::BeamDistribution() {}
BeamDistribution::~BeamDistribution() {}

void BeamDistribution::ReadIn(G4String filename){
    
    std::ifstream distrodata;//create input object
    distrodata.open(filename.c_str());//open data file
    
    if(distrodata.is_open()) {//opens as gets in here
	G4cout<<"\n---> " "Reading the distribution from "<<filename<<" ... "<<std::endl;
	
	std::cout << "Reading number of lines: "<< std::endl;
	for (fCol = 0; std::getline(distrodata, rubbish); ++fCol){}
	fCol -= 3;//expected 3 lines of text
	std::cout << "Number of usable lines in text file: " << fCol << std::endl;
	
	distrodata.clear();
	distrodata.seekg(0,std::ios::beg);//reset istream to beginning of file
	  
	for(fIter = 0; fIter<3; fIter++){//three lines of stuff in file moved
	  std::getline(distrodata,rubbish);
	}
	std::cout << "Reading values: "<< std::endl;
	
	for(int i=0; i<fCol; i++){//hard coded for now
// 	  G4cout << i << std::endl;
	  distrodata >> fDist >> fRelProb ;
	  fDistVec.push_back(fDist);
	  fRelProbVec.push_back(fRelProb);
	}
    }else std::cout << "File not opened" << std::endl;
    distrodata.close();
}

void BeamDistribution::ReadOut(){
    std::cout <<  "CALCIUM DATA: " <<  std::endl;
    fIter = 0;
    while(fIter < fCol){
      std::cout << "STOPPING POWER: " << fDistVec[fIter] << "	RELPROB: " << fRelProbVec[fIter] << std::endl;
      fIter++;
    }
}

void BeamDistribution::SumProb(){
    G4double sum = 0.;
    for(fIter = 0; fIter < fCol; fIter++){
      sum += fRelProbVec[fIter];//feed sum value (if sum=1 is relative, else divide by sum to normalise)
    }
    std::cout << "REL PROB: " << sum << std::endl;
    if(sum != 1.) std::cout << " ADJUSTING " << std::endl;
    if(sum != 1.) AdjustProb(sum);
}

void BeamDistribution::AdjustProb(double sum){
    for(fIter = 0; fIter < fCol; fIter++){
	fRelProbVec[fIter] /= sum;
    }
}

double BeamDistribution::SelectDist(G4double RelProb){
    G4double counter = 0;
    for(fIter = 0; fIter < fCol; fIter++){
      counter += fRelProbVec[fIter];
      if(counter>RelProb) {
	//std::cout << fPowerVec[fIter] << std::endl;
	return fDistVec[fIter]*(mm);
	break;
      }
    }
    std::cout << "SOMETHING WENT WRONG" << std::endl;
    return 0.;
}