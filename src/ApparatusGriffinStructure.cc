#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "ApparatusGriffinStructure.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

ApparatusGriffinStructure::ApparatusGriffinStructure() :
	// LogicalVolumes
	fSquarePieceLog(0),
	fTrianglePieceLog(0),
	fRingPieceLog(0),
	fRodPieceLog(0)
{
	fStructureMaterial = "G4_Al";
	fRodMaterial = "G4_STAINLESS-STEEL";
	fInnerDistanceToSquareFace = 470.0*mm; // inner radius to the centre of square plate - GOOD
	fSquareAndTriangleFaceThickness = 95.0*mm; // thickness of the square and triangle plates. In reality these thicknesses are different, but we will approximate them to be the same. - GOOD.
	fSquareHoleLength = 314.0*mm; // total legth of square hole, from inside structure. - GOOD.
	fSquareHoleRoundedEdgeRadius = 20.0*mm; // the roundness of the square holes. This is approximated, won't have a large effect. - OKAY
	fTriangleOuterHoleRadius = 120.65*mm; // a hole radius in the triangle face. Radius to the rounded edge, this is the larger radius. - GOOD.
	fTriangleInnerHoleRadius = 101.7*mm; // a hole radius in the triangle face. Radius to the square edge, this is the smaller radius. - GOOD.
	fSquareRodHoleRadius = (31.75*mm)/2.0; // the radius of the holes for the steel rods. - GOOD.
	fSquareRodHoleInnerRadius = 568.14*mm; // the radius from the origin to the steel rod. - GOOD.
	fSquareRodPositionFromCenterOfSquareFace = (350.0*mm)/2.0; // the position the holes are from the origin. This is NOT a radius, this is along x or y. - GOOD.
	fSquareRodTotalLength = 570.0*mm; // the total length of the rods. - GOOD.
	fSquareBgoCutLength = (118.5*mm)*2.0; // the length of the bgo slot. - GOOD.
	fSquareBgoCutWidth = 30.0*mm; // the width of the bgo slot. - GOOD.
	// ref http://en.wikipedia.org/wiki/Rhombicuboctahedron
	fSquareSquareAngle = 144.74*deg; // the angle between square faces in a rhombicuboctahedron. - GOOD.
	fSquareTriangleAngle = 135.0*deg; // the angle between a square face and a triangle face in a rhombicuboctahedron. - GOOD.
	fLargeRingFrameInnerRadius = 850.0*mm; // the inner radius of one of the the large rings that holds up the structure. - GOOD.
	fLargeRingFrameOuterRadius = 1100.0*mm; // the outer radius of one of the the large rings that holds up the structure. - GOOD.
	fLargeRingFrameThickness = 50.8*mm; // the thickness of one of the the large rings that holds up the structure. - GOOD.
	fLargeRingFrameZPosition = 240.0*mm; // the positive z distance the ring is from the origin. - GOOD.
	//G4double triangleThetaAngle = (180/M_PI)*(atan((1/sqrt(3))/sqrt((11/12) + (1/sqrt(2))))+atan((sqrt(2))/(1+sqrt(2))))*deg;
	fTriangleThetaAngle = 54.735610317245360*deg;

	// Check surfaces to determine any problematic overlaps. Turn this on to have Geant4 check the surfaces.
	// Do not leave this on, it will slow the DetectorConstruction process!
	// This was last check on May 29th, 2015. - GOOD!
	fSurfCheck = false;

	/////////////////////////////////////////////////////////////////////
	// griffinCoords for GRIFFIN
	// Note that the GRIFFIN lampshade angles are rotated by 45 degrees with respect to those of TIGRESS.
	// Modified griffinCoords for TIGRESS are below!
	/////////////////////////////////////////////////////////////////////
	// theta
	fGriffinCoords[0][0] 	= 45.0*deg;
	fGriffinCoords[1][0] 	= 45.0*deg;
	fGriffinCoords[2][0] 	= 45.0*deg;
	fGriffinCoords[3][0] 	= 45.0*deg;
	fGriffinCoords[4][0] 	= 90.0*deg;
	fGriffinCoords[5][0] 	= 90.0*deg;
	fGriffinCoords[6][0] 	= 90.0*deg;
	fGriffinCoords[7][0] 	= 90.0*deg;
	fGriffinCoords[8][0] 	= 90.0*deg;
	fGriffinCoords[9][0] 	= 90.0*deg;
	fGriffinCoords[10][0] 	= 90.0*deg;
	fGriffinCoords[11][0] 	= 90.0*deg;
	fGriffinCoords[12][0] 	= 135.0*deg;
	fGriffinCoords[13][0] 	= 135.0*deg;
	fGriffinCoords[14][0] 	= 135.0*deg;
	fGriffinCoords[15][0] 	= 135.0*deg;
	// phi
	fGriffinCoords[0][1] 	= 67.5*deg;
	fGriffinCoords[1][1] 	= 157.5*deg;
	fGriffinCoords[2][1] 	= 247.5*deg;
	fGriffinCoords[3][1] 	= 337.5*deg;
	fGriffinCoords[4][1] 	= 22.5*deg;
	fGriffinCoords[5][1] 	= 67.5*deg;
	fGriffinCoords[6][1] 	= 112.5*deg;
	fGriffinCoords[7][1] 	= 157.5*deg;
	fGriffinCoords[8][1] 	= 202.5*deg;
	fGriffinCoords[9][1] 	= 247.5*deg;
	fGriffinCoords[10][1] 	= 292.5*deg;
	fGriffinCoords[11][1] 	= 337.5*deg;
	fGriffinCoords[12][1] 	= 67.5*deg;
	fGriffinCoords[13][1] 	= 157.5*deg;
	fGriffinCoords[14][1] 	= 247.5*deg;
	fGriffinCoords[15][1] 	= 337.5*deg;
	// yaw (alpha)
	fGriffinCoords[0][2] 	= 0.0*deg;
	fGriffinCoords[1][2] 	= 0.0*deg;
	fGriffinCoords[2][2] 	= 0.0*deg;
	fGriffinCoords[3][2] 	= 0.0*deg;
	fGriffinCoords[4][2] 	= 0.0*deg;
	fGriffinCoords[5][2] 	= 0.0*deg;
	fGriffinCoords[6][2] 	= 0.0*deg;
	fGriffinCoords[7][2] 	= 0.0*deg;
	fGriffinCoords[8][2] 	= 0.0*deg;
	fGriffinCoords[9][2] 	= 0.0*deg;
	fGriffinCoords[10][2] 	= 0.0*deg;
	fGriffinCoords[11][2] 	= 0.0*deg;
	fGriffinCoords[12][2] 	= 0.0*deg;
	fGriffinCoords[13][2] 	= 0.0*deg;
	fGriffinCoords[14][2] 	= 0.0*deg;
	fGriffinCoords[15][2] 	= 0.0*deg;
	// pitch (beta)
	fGriffinCoords[0][3] 	= -45.0*deg;
	fGriffinCoords[1][3] 	= -45.0*deg;
	fGriffinCoords[2][3] 	= -45.0*deg;
	fGriffinCoords[3][3] 	= -45.0*deg;
	fGriffinCoords[4][3] 	= 0.0*deg;
	fGriffinCoords[5][3] 	= 0.0*deg;
	fGriffinCoords[6][3] 	= 0.0*deg;
	fGriffinCoords[7][3] 	= 0.0*deg;
	fGriffinCoords[8][3] 	= 0.0*deg;
	fGriffinCoords[9][3] 	= 0.0*deg;
	fGriffinCoords[10][3] 	= 0.0*deg;
	fGriffinCoords[11][3] 	= 0.0*deg;
	fGriffinCoords[12][3] 	= 45.0*deg;
	fGriffinCoords[13][3] 	= 45.0*deg;
	fGriffinCoords[14][3] 	= 45.0*deg;
	fGriffinCoords[15][3] 	= 45.0*deg;
	// roll (gamma)
	fGriffinCoords[0][4] 	= 67.5*deg;
	fGriffinCoords[1][4] 	= 157.5*deg;
	fGriffinCoords[2][4] 	= 247.5*deg;
	fGriffinCoords[3][4] 	= 337.5*deg;
	fGriffinCoords[4][4] 	= 22.5*deg;
	fGriffinCoords[5][4] 	= 67.5*deg;
	fGriffinCoords[6][4] 	= 112.5*deg;
	fGriffinCoords[7][4] 	= 157.5*deg;
	fGriffinCoords[8][4] 	= 202.5*deg;
	fGriffinCoords[9][4] 	= 247.5*deg;
	fGriffinCoords[10][4] 	= 292.5*deg;
	fGriffinCoords[11][4] 	= 337.5*deg;
	fGriffinCoords[12][4] 	= 67.5*deg;
	fGriffinCoords[13][4] 	= 157.5*deg;
	fGriffinCoords[14][4] 	= 247.5*deg;
	fGriffinCoords[15][4] 	= 337.5*deg;

	// Triangle sections for ancillary detectors coords
	// theta
	fAncillaryCoords[0][0] 	= fTriangleThetaAngle;
	fAncillaryCoords[1][0] 	= fTriangleThetaAngle;
	fAncillaryCoords[2][0] 	= fTriangleThetaAngle;
	fAncillaryCoords[3][0] 	= fTriangleThetaAngle;
	fAncillaryCoords[4][0] 	= 180.0*deg - fTriangleThetaAngle;
	fAncillaryCoords[5][0] 	= 180.0*deg - fTriangleThetaAngle;
	fAncillaryCoords[6][0] 	= 180.0*deg - fTriangleThetaAngle;
	fAncillaryCoords[7][0] 	= 180.0*deg - fTriangleThetaAngle;
	// phi
	fAncillaryCoords[0][1] 	= 22.5*deg;
	fAncillaryCoords[1][1] 	= 112.5*deg;
	fAncillaryCoords[2][1] 	= 202.5*deg;
	fAncillaryCoords[3][1] 	= 292.5*deg;
	fAncillaryCoords[4][1] 	= 22.5*deg;
	fAncillaryCoords[5][1] 	= 112.5*deg;
	fAncillaryCoords[6][1] 	= 202.5*deg;
	fAncillaryCoords[7][1] 	= 292.5*deg;
	// yaw (alpha)
	fAncillaryCoords[0][2] 	= 0.0*deg;
	fAncillaryCoords[1][2] 	= 0.0*deg;
	fAncillaryCoords[2][2] 	= 0.0*deg;
	fAncillaryCoords[3][2] 	= 0.0*deg;
	fAncillaryCoords[4][2] 	= 0.0*deg;
	fAncillaryCoords[5][2] 	= 0.0*deg;
	fAncillaryCoords[6][2] 	= 0.0*deg;
	fAncillaryCoords[7][2] 	= 0.0*deg;
	// pitch (beta)
	fAncillaryCoords[0][3] 	= fTriangleThetaAngle;
	fAncillaryCoords[1][3] 	= fTriangleThetaAngle;
	fAncillaryCoords[2][3] 	= fTriangleThetaAngle;
	fAncillaryCoords[3][3] 	= fTriangleThetaAngle;
	fAncillaryCoords[4][3] 	= 180.0*deg - fTriangleThetaAngle;
	fAncillaryCoords[5][3] 	= 180.0*deg - fTriangleThetaAngle;
	fAncillaryCoords[6][3] 	= 180.0*deg - fTriangleThetaAngle;
	fAncillaryCoords[7][3] 	= 180.0*deg - fTriangleThetaAngle;
	// roll (gamma)
	fAncillaryCoords[0][4] 	= 22.5*deg;
	fAncillaryCoords[1][4] 	= 112.5*deg;
	fAncillaryCoords[2][4] 	= 202.5*deg;
	fAncillaryCoords[3][4] 	= 292.5*deg;
	fAncillaryCoords[4][4] 	= 22.5*deg;
	fAncillaryCoords[5][4] 	= 112.5*deg;
	fAncillaryCoords[6][4] 	= 202.5*deg;
	fAncillaryCoords[7][4] 	= 292.5*deg;
}

ApparatusGriffinStructure::~ApparatusGriffinStructure() {
	delete fSquarePieceLog;
	delete fTrianglePieceLog;
	delete fRingPieceLog;
	delete fRodPieceLog;
}

G4int ApparatusGriffinStructure::Build() { 
	// Setup assembly volumes
	fAssemblySquare = new G4AssemblyVolume();
	fAssemblyTriangle = new G4AssemblyVolume();
	fAssemblyRing = new G4AssemblyVolume();
	fAssemblyRod = new G4AssemblyVolume();

	//G4cout<<"BuildSquarePiece"<<G4endl;
	BuildSquarePiece();

	//G4cout<<"BuildTrianglePiece"<<G4endl;
	BuildTrianglePiece();

	//G4cout<<"BuildLargeRingFrame"<<G4endl;
	BuildLargeRingFrame();

	//G4cout<<"BuildRodPiece"<<G4endl;
	BuildRodPiece();

	return 1;
}//end ::Build


G4int ApparatusGriffinStructure::Place(G4LogicalVolume* expHallLog, G4int selector) {
	G4ThreeVector move = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
	G4RotationMatrix* rotate = new G4RotationMatrix;

	G4double alpha;
	G4double beta;
	G4double gamma;

	G4int iSquare1 = 0;
	G4int iSquare2 = 15;
	G4int iTriangle1 =0;
	G4int iTriangle2 =7;

	if(selector == 0) { // Full structure
		iSquare1 = 0;
		iSquare2 = 15;
		iTriangle1 =0;
		iTriangle2 =7;
	} else if(selector == 1) { // upstream-half of structure: corona plus upstream lampshade
		iSquare1 = 4;
		iSquare2 = 15;
		iTriangle1 =4;
		iTriangle2 =7;
	} else if(selector == 2) { // downstream-half of structure: corona plus downstream lampshade
		iSquare1 = 0;
		iSquare2 = 11;
		iTriangle1 =0;
		iTriangle2 =3;
	} else if(selector == 3) { // just the corona, no lampshades
		iSquare1 = 4;
		iSquare2 = 11;
		iTriangle1 =4;
		iTriangle2 =3;
	}

	for(G4int i=iSquare1; i<=iSquare2; i++) { // loop over all square pieces
		alpha  = fGriffinCoords[i][2]; // yaw
		beta   = fGriffinCoords[i][3]; // pitch
		gamma  = fGriffinCoords[i][4]; // roll

		rotate = new G4RotationMatrix;
		rotate->rotateY(M_PI/2.0);
		rotate->rotateX(M_PI/2.0);
		rotate->rotateX(alpha);
		rotate->rotateY(beta);
		rotate->rotateZ(gamma);

		fAssemblySquare->MakeImprint(expHallLog, move, rotate, 0, fSurfCheck);
	}

	for(G4int i=iSquare1; i<=iSquare2; i++) { // loop over all rods
		alpha  = fGriffinCoords[i][2]; // yaw
		beta   = fGriffinCoords[i][3]; // pitch
		gamma  = fGriffinCoords[i][4]; // roll

		rotate = new G4RotationMatrix;
		rotate->rotateY(M_PI/2.0);
		rotate->rotateX(M_PI/2.0);
		rotate->rotateX(alpha);
		rotate->rotateY(beta);
		rotate->rotateZ(gamma);

		fAssemblyRod->MakeImprint(expHallLog, move, rotate, 0, fSurfCheck);
	}

	for(G4int i=iTriangle1; i<=iTriangle2; i++){ // loop over all triangle pieces
		alpha  = fAncillaryCoords[i][2]; // yaw
		beta   = fAncillaryCoords[i][3]; // pitch
		gamma  = fAncillaryCoords[i][4]; // roll

		rotate = new G4RotationMatrix;
		if(i>=4) rotate->rotateZ(60.0*deg);
		rotate->rotateY(M_PI);
		rotate->rotateY(M_PI+beta);
		rotate->rotateZ(gamma);

		fAssemblyTriangle->MakeImprint(expHallLog, move, rotate, 0, fSurfCheck);
	}

	rotate = new G4RotationMatrix;
	rotate->rotateY(0.0*deg);
	fAssemblyRing->MakeImprint(expHallLog, move, rotate, 0, fSurfCheck);

	rotate = new G4RotationMatrix;
	rotate->rotateY(180.0*deg);
	fAssemblyRing->MakeImprint(expHallLog, move, rotate, 0, fSurfCheck);

	return 1;
}

G4int ApparatusGriffinStructure::BuildSquarePiece() {
	G4Material* material = G4Material::GetMaterial(fStructureMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fStructureMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
	visAtt->SetVisibility(true);

	G4SubtractionSolid* squarePiece = SquarePiece();

	// Define rotation and movement objects
	G4ThreeVector direction = G4ThreeVector(0,0,1);
	G4double zPosition		= fInnerDistanceToSquareFace + (fSquareAndTriangleFaceThickness)/2.0;
	G4ThreeVector move 		= zPosition * direction;

	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fSquarePieceLog == nullptr) {
		fSquarePieceLog = new G4LogicalVolume(squarePiece, material, "squarePieceLog", 0, 0, 0);
		fSquarePieceLog->SetVisAttributes(visAtt);
	}

	fAssemblySquare->AddPlacedVolume(fSquarePieceLog, move, rotate);

	return 1;
}

G4int ApparatusGriffinStructure::BuildTrianglePiece() {
	G4Material* material = G4Material::GetMaterial(fStructureMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fStructureMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
	visAtt->SetVisibility(true);

	G4SubtractionSolid* trianglePiece = TrianglePiece();

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);

	G4double r          = fInnerDistanceToSquareFace +  fSquareAndTriangleFaceThickness/2.0; // to centre
	G4double x          = r*tan(22.5*deg); // half length of square side
	G4double h          = (x)*tan(30.0*deg); // height from bottom of the triangle to origin.
	G4double c          = x + h*sin(fTriangleThetaAngle);
	G4double zPosition	= ((c)/(sin(90.0*deg - fTriangleThetaAngle)))  ;

	G4ThreeVector move 		= zPosition * direction;

	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fTrianglePieceLog == nullptr) {
		fTrianglePieceLog = new G4LogicalVolume(trianglePiece, material, "trianglePieceLog", 0, 0, 0);
		fTrianglePieceLog->SetVisAttributes(visAtt);
	}

	fAssemblyTriangle->AddPlacedVolume(fTrianglePieceLog, move, rotate);

	return 1;
}

G4int ApparatusGriffinStructure::BuildLargeRingFrame()
{
	G4Material* material = G4Material::GetMaterial(fStructureMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fStructureMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
	visAtt->SetVisibility(true);

	G4Tubs* ringPiece = LargeRingFrame();

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition         = fLargeRingFrameZPosition + fLargeRingFrameThickness/2.0;
	G4ThreeVector move          = zPosition * direction;

	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fRingPieceLog == nullptr) {
		fRingPieceLog = new G4LogicalVolume(ringPiece, material, "ringPieceLog", 0, 0, 0);
		fRingPieceLog->SetVisAttributes(visAtt);
	}

	fAssemblyRing->AddPlacedVolume(fRingPieceLog, move, rotate);

	return 1;
}

G4int ApparatusGriffinStructure::BuildRodPiece() {
	G4Material* material = G4Material::GetMaterial(fRodMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fRodMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.30,0.30,0.30));
	visAtt->SetVisibility(true);

	G4Tubs* rodPiece = RodPiece();

	// Define rotation and movement objects
	G4ThreeVector move 	= G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	G4double zPosition = sqrt((fSquareRodHoleInnerRadius * fSquareRodHoleInnerRadius) - (fSquareRodPositionFromCenterOfSquareFace * fSquareRodPositionFromCenterOfSquareFace)) + fSquareRodTotalLength/2.0 ;

	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fRodPieceLog == nullptr) {
		fRodPieceLog = new G4LogicalVolume(rodPiece, material, "rodPieceLog", 0, 0, 0);
		fRodPieceLog->SetVisAttributes(visAtt);
	}

	move = G4ThreeVector(1.0*fSquareRodPositionFromCenterOfSquareFace,1.0*fSquareRodPositionFromCenterOfSquareFace,zPosition);
	fAssemblyRod->AddPlacedVolume(fRodPieceLog, move, rotate);

	move = G4ThreeVector(-1.0*fSquareRodPositionFromCenterOfSquareFace,1.0*fSquareRodPositionFromCenterOfSquareFace,zPosition);
	fAssemblyRod->AddPlacedVolume(fRodPieceLog, move, rotate);

	move = G4ThreeVector(1.0*fSquareRodPositionFromCenterOfSquareFace,-1.0*fSquareRodPositionFromCenterOfSquareFace,zPosition);
	fAssemblyRod->AddPlacedVolume(fRodPieceLog, move, rotate);

	move = G4ThreeVector(-1.0*fSquareRodPositionFromCenterOfSquareFace,-1.0*fSquareRodPositionFromCenterOfSquareFace,zPosition);
	fAssemblyRod->AddPlacedVolume(fRodPieceLog, move, rotate);

	return 1;
}

G4SubtractionSolid* ApparatusGriffinStructure::SquarePiece() {
	G4double innerRadius;
	G4double outerRadius;
	G4double startPhi = 0.0*deg;
	G4double endPhi = 360.0*deg;

	G4double r              = fInnerDistanceToSquareFace + fSquareAndTriangleFaceThickness/2.0; // to centre of square solid
	G4double x              = r*tan(22.5*deg); // half length of square side
	G4double xt             = (fSquareAndTriangleFaceThickness/2.0)*tan(22.5*deg); // add extra thickess for the edges
	G4double halfLengthX  = x + xt;
	G4double halfLengthY  = halfLengthX;
	G4double halfLengthZ  = fSquareAndTriangleFaceThickness/2.0;

	G4Box* griffinStructureSquarePlate = new G4Box("griffinStructureSquarePlate", halfLengthX, halfLengthY, halfLengthZ);

	G4double cutHalfLengthX;
	G4double cutHalfLengthY;
	G4double cutHalfLengthZ;
	G4double extraCut;

	G4ThreeVector moveCut = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	G4RotationMatrix* rotateCut;
	rotateCut = new G4RotationMatrix;

	//Now cut bgo slots
	extraCut           = 50.0*mm;
	cutHalfLengthX   = ((fSquareBgoCutWidth/2.0 + extraCut)*cos(22.5*deg));
	cutHalfLengthY   = fSquareBgoCutLength/2.0;
	cutHalfLengthZ   = fSquareAndTriangleFaceThickness;

	G4Box* cutBgo = new G4Box("cutBgo", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	// snip first bgo slot
	moveCut = G4ThreeVector(-1.0*(x-fSquareBgoCutWidth/2.0 + extraCut),0.0*mm,0.0*mm); // halfLengthX, the outer length of the original square
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateY(22.5*deg);

	G4SubtractionSolid* griffinStructureSquarePlateCut1 = new G4SubtractionSolid("griffinStructureSquarePlateCut1", griffinStructureSquarePlate, cutBgo, rotateCut, moveCut);

	// snip first bgo slot
	moveCut = G4ThreeVector(1.0*(x-fSquareBgoCutWidth/2.0 + extraCut),0.0*mm,0.0*mm); // halfLengthX, the outer length of the original square
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(180.0*deg);
	rotateCut->rotateY(22.5*deg);

	G4SubtractionSolid* griffinStructureSquarePlateCut2 = new G4SubtractionSolid("griffinStructureSquarePlateCut2", griffinStructureSquarePlateCut1, cutBgo, rotateCut, moveCut);

	// snip next bgo slot
	moveCut = G4ThreeVector(0.0*mm,1.0*(x-fSquareBgoCutWidth/2.0 + extraCut),0.0*mm); // halfLengthX, the outer length of the original square
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(90.0*deg);
	rotateCut->rotateY(22.5*deg);

	G4SubtractionSolid* griffinStructureSquarePlateCut3 = new G4SubtractionSolid("griffinStructureSquarePlateCut3", griffinStructureSquarePlateCut2, cutBgo, rotateCut, moveCut);

	// snip last bgo slot
	moveCut = G4ThreeVector(0.0*mm,-1.0*(x-fSquareBgoCutWidth/2.0 + extraCut),0.0*mm); // halfLengthX, the outer length of the original square
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(270.0*deg);
	rotateCut->rotateY(22.5*deg);

	G4SubtractionSolid* griffinStructureSquarePlateCut4 = new G4SubtractionSolid("griffinStructureSquarePlateCut4", griffinStructureSquarePlateCut3, cutBgo, rotateCut, moveCut);


	// Now cut the edges so that they fit nicely together
	cutHalfLengthX = halfLengthZ;
	cutHalfLengthY = 2.0*r;
	cutHalfLengthZ = 2.0*r;

	G4Box* cutSide = new G4Box("cutSide", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	moveCut = G4ThreeVector(-1.0*cutHalfLengthX/cos(22.5*deg),0.0*mm,-1.0*r);
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateY(22.5*deg);

	// snip 1st edge
	G4SubtractionSolid* griffinStructureSquarePlateCut5 = new G4SubtractionSolid("griffinStructureSquarePlateCut5", griffinStructureSquarePlateCut4, cutSide, rotateCut, moveCut);

	moveCut = G4ThreeVector(1.0*cutHalfLengthX/cos(22.5*deg),0.0*mm,-1.0*r);
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(180.0*deg);
	rotateCut->rotateY(22.5*deg);

	// snip next edge
	G4SubtractionSolid* griffinStructureSquarePlateCut6 = new G4SubtractionSolid("griffinStructureSquarePlateCut6", griffinStructureSquarePlateCut5, cutSide, rotateCut, moveCut);

	moveCut = G4ThreeVector(0.0*mm,1.0*cutHalfLengthX/cos(22.5*deg),-1.0*r);
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(90.0*deg);
	rotateCut->rotateY(22.5*deg);

	// snip next edge
	G4SubtractionSolid* griffinStructureSquarePlateCut7 = new G4SubtractionSolid("griffinStructureSquarePlateCut7", griffinStructureSquarePlateCut6, cutSide, rotateCut, moveCut);

	moveCut = G4ThreeVector(0.0*mm,-1.0*cutHalfLengthX/cos(22.5*deg),-1.0*r);
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(270.0*deg);
	rotateCut->rotateY(22.5*deg);

	// snip last edge
	G4SubtractionSolid* griffinStructureSquarePlateCut8 = new G4SubtractionSolid("griffinStructureSquarePlateCut8", griffinStructureSquarePlateCut7, cutSide, rotateCut, moveCut);


	// Now cut the holes for the steel rods
	innerRadius = 0.0*mm;
	outerRadius = fSquareRodHoleRadius;
	halfLengthZ = fSquareRodTotalLength/2.0;

	G4Tubs* rodHole = new G4Tubs("rodHole", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);


	G4double zPositionRod     = sqrt((fSquareRodHoleInnerRadius * fSquareRodHoleInnerRadius) - (fSquareRodPositionFromCenterOfSquareFace * fSquareRodPositionFromCenterOfSquareFace)) + fSquareRodTotalLength/2.0 ;
	G4double zPositionSquare  = fInnerDistanceToSquareFace + (fSquareAndTriangleFaceThickness)/2.0;

	// snip first hole
	moveCut = G4ThreeVector(1.0*fSquareRodPositionFromCenterOfSquareFace,1.0*fSquareRodPositionFromCenterOfSquareFace,zPositionRod-zPositionSquare);

	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid* griffinStructureSquarePlateCut9 = new G4SubtractionSolid("griffinStructureSquarePlateCut9", griffinStructureSquarePlateCut8, rodHole, rotateCut, moveCut);

	// snip next hole
	moveCut = G4ThreeVector(-1.0*fSquareRodPositionFromCenterOfSquareFace,1.0*fSquareRodPositionFromCenterOfSquareFace,zPositionRod-zPositionSquare);
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid* griffinStructureSquarePlateCut10 = new G4SubtractionSolid("griffinStructureSquarePlateCut10", griffinStructureSquarePlateCut9, rodHole, rotateCut, moveCut);

	// snip next hole
	moveCut = G4ThreeVector(1.0*fSquareRodPositionFromCenterOfSquareFace,-1.0*fSquareRodPositionFromCenterOfSquareFace,zPositionRod-zPositionSquare);
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid* griffinStructureSquarePlateCut11 = new G4SubtractionSolid("griffinStructureSquarePlateCut11", griffinStructureSquarePlateCut10, rodHole, rotateCut, moveCut);

	// snip last hole
	moveCut = G4ThreeVector(-1.0*fSquareRodPositionFromCenterOfSquareFace,-1.0*fSquareRodPositionFromCenterOfSquareFace,zPositionRod-zPositionSquare);
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid* griffinStructureSquarePlateCut12 = new G4SubtractionSolid("griffinStructureSquarePlateCut12", griffinStructureSquarePlateCut11, rodHole, rotateCut, moveCut);


	// cut out the square hole with rounded corners
	moveCut = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid* squareWithRoundedCorners = SquareWithRoundedCorners();

	G4SubtractionSolid* griffinStructureSquarePlateCut13 = new G4SubtractionSolid("griffinStructureSquarePlateCut13", griffinStructureSquarePlateCut12, squareWithRoundedCorners, rotateCut, moveCut);

	return griffinStructureSquarePlateCut13;
}


G4SubtractionSolid* ApparatusGriffinStructure::TrianglePiece()
{
	// grab the scissors, time to cut.
	G4double startPhi = 0.0*deg;
	G4double endPhi = 360.0*deg;
	G4double phi = 0.0*deg;

	G4double innerRadius;
	G4double outerRadius;
	G4double cutHalfLengthX;
	G4double cutHalfLengthY;
	G4double cutHalfLengthZ;

	G4ThreeVector move, moveCut;
	G4RotationMatrix* rotateCut;

	// to centre
	G4double r              = fInnerDistanceToSquareFace + fSquareAndTriangleFaceThickness/2.0;
	G4double x              = r*tan(22.5*deg);
	G4double xt             = (fSquareAndTriangleFaceThickness/2.0)*tan(22.5*deg); // add thickness for edges
	G4double halfLengthX  = x + xt;
	G4double halfLengthZ  = fSquareAndTriangleFaceThickness/2.0;

	// equilateral triangle is inscribed in a circle of radius a.
	innerRadius = 0.0*mm;
	outerRadius = (2.0*halfLengthX)/(sqrt(3));
	G4Tubs* griffinStructureTrianglePlate = new G4Tubs("griffinStructureTrianglePlate", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	G4SubtractionSolid* truncatedThreeSidedCylinder = TruncatedThreeSidedCylinder();

	move = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	rotateCut = new G4RotationMatrix;

	// cut out inside truncatedThreeSidedCylinder hole
	G4SubtractionSolid* griffinStructureTrianglePlateCut1 = new G4SubtractionSolid("griffinStructureTrianglePlateCut1", griffinStructureTrianglePlate, truncatedThreeSidedCylinder, rotateCut, move);

	r = fInnerDistanceToSquareFace + fSquareAndTriangleFaceThickness/2.0;
	x = r*tan(22.5*deg);

	// Now cut the edges so that the pieces fit nicely together
	cutHalfLengthX = 1.0*r;
	cutHalfLengthY = 2.0*r;
	cutHalfLengthZ = 2.0*r;
	G4Box* cutSide = new G4Box("cutSide", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	G4double theta = ((90.0*deg - fTriangleThetaAngle) - 22.5*deg);

	move = G4ThreeVector(cutHalfLengthX/cos(theta),0.0*mm,-1.0*(x*tan(30.0*deg))/(tan(theta))); // phi = 0;

	// 1st edge
	phi = 0.0*deg;

	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(-1.0*phi);
	rotateCut->rotateY(-1.0*theta);

	moveCut 	= G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

	// snip snip
	G4SubtractionSolid* griffinStructureTrianglePlateCut2 = new G4SubtractionSolid("griffinStructureTrianglePlateCut2", griffinStructureTrianglePlateCut1, cutSide, rotateCut, moveCut);

	// 2nd edge
	phi = 120.0*deg;
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(-1.0*phi);
	rotateCut->rotateY(-1.0*theta);

	moveCut 	= G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

	// snip snip
	G4SubtractionSolid* griffinStructureTrianglePlateCut3 = new G4SubtractionSolid("griffinStructureTrianglePlateCut3", griffinStructureTrianglePlateCut2, cutSide, rotateCut, moveCut);

	// last edge
	phi = 240.0*deg;

	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(-1.0*phi);
	rotateCut->rotateY(-1.0*theta);

	moveCut 	= G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

	// snip snip
	G4SubtractionSolid* griffinStructureTrianglePlateCut4 = new G4SubtractionSolid("griffinStructureTrianglePlateCut4", griffinStructureTrianglePlateCut3, cutSide, rotateCut, moveCut);

	return griffinStructureTrianglePlateCut4;
}

G4Tubs* ApparatusGriffinStructure::LargeRingFrame()
{
	G4double startPhi = 0.0*deg;
	G4double endPhi = 360.0*deg;

	G4double innerRadius = fLargeRingFrameInnerRadius;
	G4double outerRadius = fLargeRingFrameOuterRadius;
	G4double halfLengthZ = fLargeRingFrameThickness/2.0;

	G4Tubs* ringCylinder = new G4Tubs("ringCylinder", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return ringCylinder;
}

G4Tubs* ApparatusGriffinStructure::RodPiece()
{
	G4double startPhi = 0.0*deg;
	G4double endPhi = 360.0*deg;

	G4double innerRadius = 0.0*mm;
	G4double outerRadius = fSquareRodHoleRadius;
	G4double halfLengthZ = fSquareRodTotalLength/2.0;

	G4Tubs* rod = new G4Tubs("rod", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return rod;
}

G4SubtractionSolid* ApparatusGriffinStructure::SquareWithRoundedCorners()
{
	G4double halfLengthX = fSquareHoleLength/2.0;
	G4double halfLengthY = fSquareHoleLength/2.0;
	G4double halfLengthZ = fSquareAndTriangleFaceThickness;

	G4ThreeVector move = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	G4RotationMatrix* rotate;

	G4Box* griffinStructureSquarePlate1 = new G4Box("griffinStructureSquarePlate1", halfLengthX, halfLengthY, halfLengthX);

	// cut rounded corners
	G4double startPhi = 0.0*deg;
	G4double endPhi = 90.0*deg;

	G4double innerRadius = fSquareHoleRoundedEdgeRadius;
	G4double outerRadius = fSquareHoleRoundedEdgeRadius + 2.0*halfLengthX;

	G4Tubs* roundedCorner = new G4Tubs("roundedCorner", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// 1st corner
	move = G4ThreeVector(1.0*(halfLengthX - fSquareHoleRoundedEdgeRadius),1.0*(halfLengthX - fSquareHoleRoundedEdgeRadius),0.0*mm);
	rotate = new G4RotationMatrix;

	G4SubtractionSolid* griffinStructureSquarePlate2 = new G4SubtractionSolid("griffinStructureSquarePlate2", griffinStructureSquarePlate1, roundedCorner, rotate, move);

	// next corner
	move = G4ThreeVector(1.0*(halfLengthX - fSquareHoleRoundedEdgeRadius),-1.0*(halfLengthX - fSquareHoleRoundedEdgeRadius),0.0*mm);
	rotate = new G4RotationMatrix;
	rotate->rotateZ(90.0*deg);

	G4SubtractionSolid* griffinStructureSquarePlate3 = new G4SubtractionSolid("griffinStructureSquarePlate3", griffinStructureSquarePlate2, roundedCorner, rotate, move);

	// next corner
	move = G4ThreeVector(-1.0*(halfLengthX - fSquareHoleRoundedEdgeRadius),-1.0*(halfLengthX - fSquareHoleRoundedEdgeRadius),0.0*mm);
	rotate = new G4RotationMatrix;
	rotate->rotateZ(180.0*deg);

	G4SubtractionSolid* griffinStructureSquarePlate4 = new G4SubtractionSolid("griffinStructureSquarePlate4", griffinStructureSquarePlate3, roundedCorner, rotate, move);

	// last corner
	move = G4ThreeVector(-1.0*(halfLengthX - fSquareHoleRoundedEdgeRadius),1.0*(halfLengthX - fSquareHoleRoundedEdgeRadius),0.0*mm);
	rotate = new G4RotationMatrix;
	rotate->rotateZ(270.0*deg);

	G4SubtractionSolid* griffinStructureSquarePlate5 = new G4SubtractionSolid("griffinStructureSquarePlate5", griffinStructureSquarePlate4, roundedCorner, rotate, move);

	return griffinStructureSquarePlate5;
}

G4SubtractionSolid* ApparatusGriffinStructure::TruncatedThreeSidedCylinder()
{
	G4double startPhi = 0.0*deg;
	G4double endPhi = 360.0*deg;

	G4double innerRadius = 0.0*mm;
	G4double outerRadius = fTriangleOuterHoleRadius;
	G4double halfLengthZ = fSquareAndTriangleFaceThickness;

	G4double phi;
	G4ThreeVector move, moveCut;
	G4RotationMatrix* rotateCut;

	G4Tubs* holeCylinder = new G4Tubs("holeCylinder", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// cut flat sides
	G4double cutHalfLengthX = 2.0*outerRadius;
	G4double cutHalfLengthY = 2.0*outerRadius;
	G4double cutHalfLengthZ = 2.0*halfLengthZ;

	G4Box* cutEdge = new G4Box("cutEdge", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	move = G4ThreeVector(fTriangleInnerHoleRadius + cutHalfLengthX, 0.0*mm,0.0*mm);

	// 1st cut
	phi = 0.0*deg;
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(-1.0*phi);

	moveCut = G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

	G4SubtractionSolid* holeCylinder1 = new G4SubtractionSolid("holeCylinder1", holeCylinder, cutEdge, rotateCut, moveCut);

	// 2nd cut
	phi = 120.0*deg;
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(-1.0*phi);

	moveCut = G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

	G4SubtractionSolid* holeCylinder2 = new G4SubtractionSolid("holeCylinder2", holeCylinder1, cutEdge, rotateCut, moveCut);

	// 3rd cut
	phi = 240.0*deg;
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(-1.0*phi);

	moveCut = G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

	G4SubtractionSolid* holeCylinder3 = new G4SubtractionSolid("holeCylinder3", holeCylinder2, cutEdge, rotateCut, moveCut);

	return holeCylinder3;
}
