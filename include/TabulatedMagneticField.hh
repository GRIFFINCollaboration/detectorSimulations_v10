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

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ios.hh"

#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class TabulatedMagneticField
#ifndef STANDALONE
: public G4MagneticField
#endif
{
  
  // Storage space for the table
  vector< vector< vector< G4double > > > xField;
  vector< vector< vector< G4double > > > yField;
  vector< vector< vector< G4double > > > zField;
  // The dimensions of the table
  G4int nx,ny,nz; 
  // The physical limits of the defined region
  G4double minx, maxx, miny, maxy, minz, maxz;
  //The physical limits of the defined region
  G4double maxbx=0, maxby=0, maxbz=0;
  // The physical extent of the defined region
  G4double dx, dy, dz;
  G4double fZoffset;
  G4double fZrotation;
  G4bool invertX, invertY, invertZ;

public:
  TabulatedMagneticField(const char* filename, G4double zOffset, G4double zRotation  );
  void  GetFieldValue( const  G4double Point[4],
		       G4double *Bfield          ) const;
};

