#include "TabulatedMagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"

#include <algorithm> //for max_element

TabulatedMagneticField::TabulatedMagneticField( const char* filename, G4double zOffset, G4double zRotation ) //called in 'nonuniformfield'
  :fZoffset(zOffset),fZrotation(zRotation),invertX(false),invertY(false),invertZ(false)
{    
 
  G4double lenUnit= mm;
  G4double fieldUnit= tesla; 
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      Magnetic field"
	 << "\n-----------------------------------------------------------";
    
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl; 

  //G4cout << "\n ----> Read ? "<< endl;  G4cin.get(); //what does this do?
 
  ifstream file( filename ); // Open the file for reading.
  
  // Ignore first blank line
  char buffer[256];
  file.getline(buffer,256);
  
  if (!file.is_open()) {
	  G4cout <<"\n\ncannot open file : " << filename <<"\n\n" ;
	  G4cin.get();
	  }


  // Read table dimensions 
  file >> nx >> ny >> nz; // Note dodgy order //20/7 what? xyz looks fine!

  G4cout << "  [ Number of values x,y,z: " 
	 << nx << " " << ny << " " << nz << " ] "//111*111*106 gives 1306026 rows
	 << endl;

  // Set up storage space for table using values from above- does not initialise values
  xField.resize( nx );
  yField.resize( nx );
  zField.resize( nx );
  G4int ix, iy, iz; //iterators for below
  for (ix=0; ix<nx; ix++) {
    xField[ix].resize(ny);
    yField[ix].resize(ny);
    zField[ix].resize(ny);
    for (iy=0; iy<ny; iy++) {
      xField[ix][iy].resize(nz);
      yField[ix][iy].resize(nz);
      zField[ix][iy].resize(nz);
    }
  }//three field values per point, hence three arrays needed
  
  // Ignores other header information    
  // The first line whose second character is '0' is considered to
  // be the last line of the header.
  do {
    file.getline(buffer,256);
  } while ( buffer[1]!='0');
  
  // Read in the data and fill arrays
  G4double xval,yval,zval,bx,by,bz;
  G4double permeability; // Not used in this example. (is 0 for spice)
  for (iz=0; iz<nz; iz++) {
    for (iy=0; iy<ny; iy++) {
      for (ix=0; ix<nx; ix++) {
        file >> xval >> yval >> zval >> bx >> by >> bz >> permeability;
	bx =bx/1000.;
	by =by/1000.;
	bz =bz/1000.;
        if ( ix==0 && iy==0 && iz==0 ) {//min is first value?? theres no sort/search 
          minx = xval * lenUnit;
          miny = yval * lenUnit;
          minz = zval * lenUnit;
        }
        xField[ix][iy][iz] = bx * fieldUnit;
        yField[ix][iy][iz] = by * fieldUnit;
        zField[ix][iy][iz] = bz * fieldUnit;
	
	if(bx > maxbx) maxbx = bx; // no unit as read-out in terminal
	if(by > maxby) maxby = by;
	if(bz > maxbz) maxbz = bz;
      }
    }
  }
  file.close(); //internally stored values, so can close file
  
  //my attempts - max field value in each column
  G4cout << "\t\t\tMax Bx value = "<< maxbx << " Tesla." << G4endl;
  G4cout << "\t\t\tMax By value = "<< maxby << " Tesla." << G4endl;
  G4cout << "\t\t\tMax Bz value = "<< maxbz << " Tesla." << G4endl;

  maxx = xval * lenUnit; //now max dimension values are the last values - i.e. post-loop values
  maxy = yval * lenUnit;
  maxz = zval * lenUnit;

  G4cout << "\n ---> ... done reading " << endl;

  // G4cout << " Read values of field from file " << filename << endl; 
  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << "\n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm "
	 << "\n ---> The field will be offset by " << zOffset/cm << " cm " << endl;

  // Should really check that the limits are not the wrong way around.
  if (maxx < minx) {swap(maxx,minx); invertX = true;} 
  if (maxy < miny) {swap(maxy,miny); invertY = true;} 
  if (maxz < minz) {swap(maxz,minz); invertZ = true;} 
  G4cout << "\nAfter reordering if neccesary"  
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << " \n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm ";

  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
  G4cout << "\n ---> Dif values x,y,z (range): " 
	 << dx/cm << " " << dy/cm << " " << dz/cm << " cm in z "
	 << "\n-----------------------------------------------------------" << endl;
}

void TabulatedMagneticField::GetFieldValue(const G4double point[4],
				      G4double *Bfield ) const
{

 //G4cout << "---> values x,y,z, t (GetFieldValue): " << point[0]/mm << " " << point[1]/mm << " " << point[2]/mm << " " << point[3]/nanosecond << " " << G4endl ; 
	 
  G4double x = point[0];
  G4double y = point[1];
  G4double z = point[2] + fZoffset;
  
  //Rotation treatment Mhd : 25 Mar 2015
	  /*
	  x,y,z   ---Rot (angle) --->  x',y',z'
	  			            		 |
	  					   			 v
	 Bx,By,Bz <---Rot (-angle)--- Bx',By',Bz'	 
	  */
  
  // Rotate the position here : Mhd : 25 Mar 2015 
  	G4ThreeVector R( x,  y,  z); 
	R.rotateZ(-fZrotation*deg); // rotation made in the opposite direction of the lens
	x = R.getX();
	y = R.getY();
	z = R.getZ();
         

  // Check that the point is within the defined region 
  if ( x>=minx && x<=maxx &&
       y>=miny && y<=maxy && 
       z>=minz && z<=maxz ) {
        
    // Position of given point within region, normalized to the range
    // [0,1]
    G4double xfraction = (x - minx) / dx;
    G4double yfraction = (y - miny) / dy; 
    G4double zfraction = (z - minz) / dz;

    if (invertX) { xfraction = 1 - xfraction;}
    if (invertY) { yfraction = 1 - yfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    G4double xdindex, ydindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    G4double xlocal = ( std::modf(xfraction*(nx-1), &xdindex));
    G4double ylocal = ( std::modf(yfraction*(ny-1), &ydindex));
    G4double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));
    
    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);
    

#ifdef DEBUG_INTERPOLATING_FIELD
    G4cout << "Local x,y,z: " << xlocal << " " << ylocal << " " << zlocal << endl;
    G4cout << "Index x,y,z: " << xindex << " " << yindex << " " << zindex << endl;
    G4double valx0z0, mulx0z0, valx1z0, mulx1z0;
    G4double valx0z1, mulx0z1, valx1z1, mulx1z1;
    valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
    valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
    valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
    valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
#endif

        // Full 3-dimensional version
    Bfield[0] =
      xField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      xField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      xField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      xField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      xField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      xField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      xField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      xField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[1] =
      yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      yField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      yField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      yField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[2] =
      zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      zField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      zField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
      
      // Rotate the BField here : Mhd : 25 Mar 2015  
	G4ThreeVector B(Bfield[0],  Bfield[1],  Bfield[2]);
	B.rotateZ(fZrotation*deg); // rotation made in the same direction of the lens
	Bfield[0] = B.getX();
	Bfield[1] = B.getY();
	Bfield[2] = B.getZ();

  } else {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
  }
  
}
