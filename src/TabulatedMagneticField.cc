#include "TabulatedMagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"

#include <algorithm> //for max_element

TabulatedMagneticField::TabulatedMagneticField(const char* filename, G4double zOffset, G4double zRotation) //called in 'nonuniformfield'
: fZoffset(zOffset), fZrotation(zRotation), fInvertX(false), fInvertY(false), fInvertZ(false)
{    
	G4double lenUnit= mm;
	G4double fieldUnit= tesla; 
	G4cout<<"\n-----------------------------------------------------------"
		<<"\n      Magnetic field"
		<<"\n-----------------------------------------------------------";

	G4cout<<"\n ---> " "Reading the field grid from "<<filename<<" ... "<<std::endl; 

	std::ifstream file(filename); // Open the file for reading.

	// Ignore first blank line
	char buffer[256];
	file.getline(buffer,256);

	if(!file.is_open()) {
		G4cout<<"\n\ncannot open file : "<<filename<<"\n\n" ;
		G4cin.get();
	}


	// Read table dimensions 
	file>>fNx>>fNy>>fNz; // Note dodgy order //20/7 what? xyz looks fine!

	G4cout<<"  [ Number of values x,y,z: " 
		<<fNx<<" "<<fNy<<" "<<fNz<<" ] "//111*111*106 gives 1306026 rows
		<<std::endl;

	// Set up storage space for table using values from above- does not initialise values
	fXField.resize(fNx);
	fYField.resize(fNx);
	fZField.resize(fNx);
	for(G4int ix=0; ix<fNx; ix++) {
		fXField[ix].resize(fNy);
		fYField[ix].resize(fNy);
		fZField[ix].resize(fNy);
		for(G4int iy=0; iy<fNy; iy++) {
			fXField[ix][iy].resize(fNz);
			fYField[ix][iy].resize(fNz);
			fZField[ix][iy].resize(fNz);
		}
	}//three field values per point, hence three arrays needed

	// Ignores other header information    
	// The first line whose second character is '0' is considered to
	// be the last line of the header.
	do {
		file.getline(buffer,256);
	} while (buffer[1] != '0');

	// Read in the data and fill arrays
	G4double xval,yval,zval,bx,by,bz;
	G4double permeability; // Not used in this example. (is 0 for spice)
	for(G4int iz=0; iz<fNz; iz++) {
		for(G4int iy=0; iy<fNy; iy++) {
			for(G4int ix=0; ix<fNx; ix++) {
				file >> xval >> yval >> zval >> bx >> by >> bz >> permeability;
				bx =bx/1000.;
				by =by/1000.;
				bz =bz/1000.;
				if(ix==0 && iy==0 && iz==0) {//min is first value?? theres no sort/search 
					fMinx = xval * lenUnit;
					fMiny = yval * lenUnit;
					fMinz = zval * lenUnit;
				}
				fXField[ix][iy][iz] = bx * fieldUnit;
				fYField[ix][iy][iz] = by * fieldUnit;
				fZField[ix][iy][iz] = bz * fieldUnit;

				if(bx > fMaxbx) fMaxbx = bx; // no unit as read-out in terminal
				if(by > fMaxby) fMaxby = by;
				if(bz > fMaxbz) fMaxbz = bz;
			}
		}
	}
	file.close(); //internally stored values, so can close file

	//my attempts - max field value in each column
	G4cout<<"\t\t\tMax Bx value = "<< fMaxbx<<" Tesla."<<G4endl;
	G4cout<<"\t\t\tMax By value = "<< fMaxby<<" Tesla."<<G4endl;
	G4cout<<"\t\t\tMax Bz value = "<< fMaxbz<<" Tesla."<<G4endl;

	fMaxx = xval * lenUnit; //now max dimension values are the last values - i.e. post-loop values
	fMaxy = yval * lenUnit;
	fMaxz = zval * lenUnit;

	G4cout<<"\n ---> ... done reading "<<std::endl;

	// G4cout<<" Read values of field from file "<<filename<<std::endl; 
	G4cout<<" ---> assumed the order:  x, y, z, Bx, By, Bz "
		<< "\n ---> Min values x,y,z: " 
		<< fMinx/cm<<" "<<fMiny/cm<<" "<<fMinz/cm<<" cm "
		<< "\n ---> Max values x,y,z: " 
		<< fMaxx/cm<<" "<<fMaxy/cm<<" "<<fMaxz/cm<<" cm "
		<< "\n ---> The field will be offset by "<<zOffset/cm<<" cm "<<std::endl;

	// Should really check that the limits are not the wrong way around.
	if(fMaxx < fMinx) {
		std::swap(fMaxx,fMinx);
		fInvertX = true;
	} 
	if(fMaxy < fMiny) {
		std::swap(fMaxy,fMiny);
		fInvertY = true;
	} 
	if(fMaxz < fMinz) {
		std::swap(fMaxz,fMinz); 
		fInvertZ = true;
	} 
	G4cout<<"\nAfter reordering if neccesary"  
		<<"\n ---> Min values x,y,z: " 
		<<fMinx/cm<<" "<<fMiny/cm<<" "<<fMinz/cm<<" cm "
		<<" \n ---> Max values x,y,z: " 
		<<fMaxx/cm<<" "<<fMaxy/cm<<" "<<fMaxz/cm<<" cm ";

	fDx = fMaxx - fMinx;
	fDy = fMaxy - fMiny;
	fDz = fMaxz - fMinz;
	G4cout<<"\n ---> Dif values x,y,z (range): " 
		<<fDx/cm<<" "<<fDy/cm<<" "<<fDz/cm<<" cm in z "
		<<"\n-----------------------------------------------------------"<<std::endl;
}

void TabulatedMagneticField::GetFieldValue(const G4double point[4], G4double* Bfield) const
{
	G4double x = point[0];
	G4double y = point[1];
	G4double z = point[2] + fZoffset;

	//Rotation treatment Mhd : 25 Mar 2015
	//
	//  x,y,z   ---Rot (angle) --->  x',y',z'
	//  |
	//  v
	//  Bx,By,Bz <---Rot (-angle)--- Bx',By',Bz'	 
	//

	// Rotate the position here : Mhd : 25 Mar 2015 
	G4ThreeVector R(x,  y,  z); 
	R.rotateZ(-fZrotation*deg); // rotation made in the opposite direction of the lens
	x = R.getX();
	y = R.getY();
	z = R.getZ();

	// Check that the point is within the defined region 
	if(x>=fMinx && x<=fMaxx &&
			y>=fMiny && y<=fMaxy && 
			z>=fMinz && z<=fMaxz) {

		// Position of given point within region, normalized to the range
		// [0,1]
		G4double xfraction = (x - fMinx)/fDx;
		G4double yfraction = (y - fMiny)/fDy; 
		G4double zfraction = (z - fMinz)/fDz;

		if(fInvertX) { 
			xfraction = 1 - xfraction;
		}
		if(fInvertY) { 
			yfraction = 1 - yfraction;
		}
		if(fInvertZ) { 
			zfraction = 1 - zfraction;
		}

		// Need addresses of these to pass to modf below.
		// modf uses its second argument as an OUTPUT argument.
		G4double xdindex, ydindex, zdindex;

		// Position of the point within the cuboid defined by the
		// nearest surrounding tabulated points
		G4double xlocal = std::modf(xfraction*(fNx-1), &xdindex);
		G4double ylocal = std::modf(yfraction*(fNy-1), &ydindex);
		G4double zlocal = std::modf(zfraction*(fNz-1), &zdindex);

		// The indices of the nearest tabulated point whose coordinates
		// are all less than those of the given point
		int xindex = static_cast<int>(xdindex);
		int yindex = static_cast<int>(ydindex);
		int zindex = static_cast<int>(zdindex);


#ifdef DEBUG_INTERPOLATING_FIELD
		G4cout<<"Local x,y,z: "<<xlocal<<" "<<ylocal<<" "<<zlocal<<std::endl;
		G4cout<<"Index x,y,z: "<<xindex<<" "<<yindex<<" "<<zindex<<std::endl;
		G4double valx0z0, mulx0z0, valx1z0, mulx1z0;
		G4double valx0z1, mulx0z1, valx1z1, mulx1z1;
		valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
		valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
		valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
		valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
#endif

		// Full 3-dimensional version
		Bfield[0] =
			fXField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
			fXField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
			fXField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
			fXField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
			fXField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
			fXField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
			fXField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
			fXField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
		Bfield[1] =
			fYField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
			fYField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
			fYField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
			fYField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
			fYField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
			fYField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
			fYField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
			fYField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
		Bfield[2] =
			fZField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
			fZField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
			fZField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
			fZField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
			fZField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
			fZField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
			fZField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
			fZField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

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

