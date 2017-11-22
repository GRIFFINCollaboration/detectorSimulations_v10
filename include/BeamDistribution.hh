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
      //Functions used, all called externally
      void ReadIn(G4String);//gets 2 columns of values from standardised file
      void ReadOut();//reads names to termial for debugging (for now)
      void SumProb();//checks if normalised
      void AdjustProb(G4double);//normalises probabilities
      G4double SelectDist(G4double); //Called by PGA and selects the stopping distance for input to the particle gun
      
      //Variable called externally
      void GetThickness(G4double val){ fSpiceTargetThickness = val; };

  private:
      std::vector<G4double> fDistVec, fRelProbVec; //these vectors will hld all of the data that is read in
      std::string fInfo; //this is information in the data file that is not used to create the vectors
      G4double  fDist, fRelProb; //These will hold the data as the valeus are read in
      G4int fCol, fIter, fUsedLines;//fIter is iterator used throughout
      G4bool fDigitBool; //When the data file starts (and extra info is finished) the bool will flip
 			  // telling the program to start collecting data from the file
      G4double fSpiceTargetThickness;
};
#endif