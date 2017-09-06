//
// Code developed by:
//  S.Larsson and J. Generowicz.
//
//    *************************************
//    *                                   *
//    *    PurgMagTabulatedField3D.hh     *
//    *                                   *
//    *************************************
//
// $Id: PurgMagTabulatedField3D.hh,v 1.3 2006-06-29 16:06:05 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-01 $
//

#ifndef TABULATEDMAGNETICFIELD_HH
#define TABULATEDMAGNETICFIELD_HH

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ios.hh"

#include <fstream>
#include <vector>
#include <cmath>

class TabulatedMagneticField
#ifndef STANDALONE
: public G4MagneticField
#endif
{
public:
  TabulatedMagneticField(const char* filename, G4double zOffset, G4double zRotation);
  void  GetFieldValue(const G4double Point[4], G4double* Bfield) const;
  
private:
  // Storage space for the table
  std::vector<std::vector<std::vector<G4double> > > fXField;
  std::vector<std::vector<std::vector<G4double> > > fYField;
  std::vector<std::vector<std::vector<G4double> > > fZField;
  // The dimensions of the table
  G4int fNx;
  G4int fNy;
  G4int fNz; 
  // The physical limits of the defined region
  G4double fMinx;
  G4double fMaxx;
  G4double fMiny;
  G4double fMaxy;
  G4double fMinz;
  G4double fMaxz;
  //The physical limits of the defined region
  G4double fMaxbx;
  G4double fMaxby;
  G4double fMaxbz;
  // The physical extent of the defined region
  G4double fDx;
  G4double fDy;
  G4double fDz;
  G4double fZoffset;
  G4double fZrotation;
  G4bool   fInvertX;
  G4bool   fInvertY;
  G4bool   fInvertZ;
};
#endif
