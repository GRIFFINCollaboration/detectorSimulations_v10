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
      ~BeamDistribution();

  public:
      void ReadIn(G4String);//gets 2 columns of values from standardised file
      void ReadOut();//reads names (for now)
      void SumProb();//checks if normalised
      void AdjustProb(G4double);//normalises probabilities
      G4double SelectDist(G4double);

  private:
      std::vector<G4double> fDistVec, fRelProbVec;
      std::string rubbish;
      G4double  fDist, fRelProb;
      G4int fCol, fIter;//fIter is iterator used throughout
};
#endif