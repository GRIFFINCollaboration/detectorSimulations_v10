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
  vector< vector< vector< double > > > xField;
  vector< vector< vector< double > > > yField;
  vector< vector< vector< double > > > zField;
  // The dimensions of the table
  int nx,ny,nz; 
  // The physical limits of the defined region
  double minx, maxx, miny, maxy, minz, maxz;
  // The physical extent of the defined region
  double dx, dy, dz;
  double fZoffset;
  double fZrotation;
  bool invertX, invertY, invertZ;

public:
  TabulatedMagneticField(const char* filename, G4double zOffset, G4double zRotation  );
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;
};

