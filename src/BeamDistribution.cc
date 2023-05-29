#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <cctype> //isdigit etc

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "BeamDistribution.hh"

BeamDistribution::BeamDistribution() { 
  fLoadGood = false;
}

BeamDistribution::~BeamDistribution() {}

BeamDistribution::BeamDistribution(G4String filename){//Reads in the data from the file
	fLoadGood = false;
	LoadDistribution(filename);
}


void BeamDistribution::LoadDistribution(G4String filename){//Reads in the data from the file
	fLoadGood = false;
	std::ifstream distrodata;//create input object
	distrodata.open(filename.c_str());//open data file
	if(!distrodata.is_open()){
		G4cout<<"File "<<filename<<" not opened."<<G4endl; 
		return;
	}

	double x,prob;
	std::vector<double> xv,pv;
	double probsum = 0.;

	while(distrodata>>x>>prob){
		xv.push_back(x);
		probsum+=prob;
		pv.push_back(prob);
	}
	distrodata.close();

	if(xv.size()<2){
		G4cout<<"Insufficient data points."<<G4endl; 
		return;
	}
		
	fDistStep=(xv[1]-xv[0])*um;

	double probstep=probsum;
	for(unsigned int i=0;i<xv.size();i++){
		fDistVec.push_back(xv[i]*um-fDistStep*0.5);
		probstep-=pv[i];
		fRelProbVec.push_back(probstep/probsum);
	}
	fLoadGood = true;
}

G4double BeamDistribution::GetRandom(){
	if(!fLoadGood)return 0;
	G4double prob=G4UniformRand();
	for(unsigned int i=0;i<fRelProbVec.size();i++){
		if(prob>fRelProbVec[i]){
			return fDistVec[i]+fDistStep*G4UniformRand();
		}
	}
	return 0;
}
