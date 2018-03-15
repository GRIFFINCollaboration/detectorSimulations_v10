#ifndef BEAMDISTRIBUTION_HH
#define BEAMDISTRIBUTION_HH

#include "globals.hh"
#include "G4ios.hh"

#include <fstream>
#include <vector>
#include <cmath>
#include <string>

class BeamDistribution
{
  public:
      BeamDistribution();
      BeamDistribution(G4String);
      ~BeamDistribution();

  public:
     void LoadDistribution(G4String filename);
     G4double GetRandom();

  private:
      bool fLoadGood;
      std::vector<G4double> fDistVec, fRelProbVec; //these vectors will hld all of the data that is read in
      G4double  fDistStep;

  public:
	  bool Good(){return fLoadGood;}
	  
}; 
#endif
