#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "Global.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemAncillaryBGO.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


DetectionSystemAncillaryBGO::DetectionSystemAncillaryBGO() :
	// LogicalVolumes
	fDetectorVolumeLog(0),
	fBgoBlockLog(0),
	fVacuumBlockLog(0),
	fCanCylinderLog(0),
	fCanCylinderBackCoverLog(0),
	fHevimetBlockLog(0)
{
	// can //
	fCanLength              = 190.0*mm;
	fCanThickness           = 0.5*mm;
	fCanThicknessFront     = 3.0*mm;
	fCanInnerRadius        = 32.5*mm;
	fCanMaterial            = "G4_Al";
	// BGO //
	fBgoLength              = 106.5*mm;
	fBgoThickness           = 20.0*mm;
	fBgoMaterial            = "G4_BGO";
	// Gap (Vacuum) //
	fGapThickness           = 0.5*mm;
	fGapThicknessOuter     = 3.5*mm;
	fGapMaterial            = "Vacuum";

	// distance used in chopping process //
	fCanFaceToBuildOrigin= 165.0*mm;
	// Chopping angles //
	fChoppingOuterAngle    = 12.76*deg;
	fChoppingInnerAngle    = 8.0*deg;
	// Hevimet //
	fHevimetThickness       = 10.0*mm;
	fHevimetMaterial        = "Hevimetal";

	//G4double triangleThetaAngle = (180/M_PI)*(atan((1/sqrt(3))/sqrt((11/12) + (1/sqrt(2))))+atan((sqrt(2))/(1+sqrt(2))))*deg;
	G4double triangleThetaAngle = 54.735610317245360*deg;

	// theta
	fDetectorAngles[0][0] 	= triangleThetaAngle;
	fDetectorAngles[1][0] 	= triangleThetaAngle;
	fDetectorAngles[2][0] 	= triangleThetaAngle;
	fDetectorAngles[3][0] 	= triangleThetaAngle;
	fDetectorAngles[4][0] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[5][0] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[6][0] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[7][0] 	= 180.0*deg - triangleThetaAngle;
	// phi
	fDetectorAngles[0][1] 	= 22.5*deg;
	fDetectorAngles[1][1] 	= 112.5*deg;
	fDetectorAngles[2][1] 	= 202.5*deg;
	fDetectorAngles[3][1] 	= 292.5*deg;
	fDetectorAngles[4][1] 	= 22.5*deg;
	fDetectorAngles[5][1] 	= 112.5*deg;
	fDetectorAngles[6][1] 	= 202.5*deg;
	fDetectorAngles[7][1] 	= 292.5*deg;
	// yaw (alpha)
	fDetectorAngles[0][2] 	= 0.0*deg;
	fDetectorAngles[1][2] 	= 0.0*deg;
	fDetectorAngles[2][2] 	= 0.0*deg;
	fDetectorAngles[3][2] 	= 0.0*deg;
	fDetectorAngles[4][2] 	= 0.0*deg;
	fDetectorAngles[5][2] 	= 0.0*deg;
	fDetectorAngles[6][2] 	= 0.0*deg;
	fDetectorAngles[7][2] 	= 0.0*deg;
	// pitch (beta)
	fDetectorAngles[0][3] 	= triangleThetaAngle;
	fDetectorAngles[1][3] 	= triangleThetaAngle;
	fDetectorAngles[2][3] 	= triangleThetaAngle;
	fDetectorAngles[3][3] 	= triangleThetaAngle;
	fDetectorAngles[4][3] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[5][3] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[6][3] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[7][3] 	= 180.0*deg - triangleThetaAngle;
	// roll (gamma)
	fDetectorAngles[0][4] 	= 22.5*deg;
	fDetectorAngles[1][4] 	= 112.5*deg;
	fDetectorAngles[2][4] 	= 202.5*deg;
	fDetectorAngles[3][4] 	= 292.5*deg;
	fDetectorAngles[4][4] 	= 22.5*deg;
	fDetectorAngles[5][4] 	= 112.5*deg;
	fDetectorAngles[6][4] 	= 202.5*deg;
	fDetectorAngles[7][4] 	= 292.5*deg;
}

DetectionSystemAncillaryBGO::~DetectionSystemAncillaryBGO() {
	// LogicalVolumes
	delete fDetectorVolumeLog;
	delete fBgoBlockLog;
	delete fVacuumBlockLog;
	delete fCanCylinderLog;
	delete fCanCylinderBackCoverLog;
	delete fHevimetBlockLog;
}


G4int DetectionSystemAncillaryBGO::Build() {
	// Build assembly volume
	fAssembly = new G4AssemblyVolume();
	fAssemblyHevimet = new G4AssemblyVolume();

	G4cout<<"BuildBGOVolume"<<G4endl;
	BuildBGOVolume();
	G4cout<<"BuildVacuumVolume"<<G4endl;
	BuildVacuumVolume();
	G4cout<<"BuildAluminumCanVolume"<<G4endl;
	BuildAluminumCanVolume();
	G4cout<<"BuildAluminumCanBackCoverVolume"<<G4endl;
	BuildAluminumCanBackCoverVolume();
	G4cout<<"BuildHevimetPiece"<<G4endl;
	BuildHevimetPiece();

	//  G4cout<<"BuildAluminumCanVolume"<<G4endl;
	//  BuildAluminumCanVolume();
	//  G4cout<<"BuildPackingVolume"<<G4endl;
	//  BuildPackingVolume();
	//  G4cout<<"BuildDiscVolume"<<G4endl;
	//  BuildDiscVolume();
	//  G4cout<<"BuildSealVolume"<<G4endl;
	//  BuildSealVolume();

	return 1;
}

G4int DetectionSystemAncillaryBGO::PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber, G4double radialpos, G4int hevimetopt) {
	G4int detectorCopyID = 0;

	G4cout<<"AncillaryBGO Detector Number = "<<detectorNumber<<G4endl;

	G4int copyNumber = detectorCopyID + detectorNumber;

	G4double position = radialpos + (fCanLength / 2.0) ;

	G4double theta  = fDetectorAngles[detectorNumber][0];
	G4double phi    = fDetectorAngles[detectorNumber][1];
	//G4double alpha  = fDetectorAngles[detectorNumber][2]; // yaw
	G4double beta   = fDetectorAngles[detectorNumber][3]; // pitch
	G4double gamma  = fDetectorAngles[detectorNumber][4]; // roll

	G4double x = 0;
	G4double y = 0;
	G4double z = position;

	G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
	G4ThreeVector move;
	for(G4int i=0; i<3; i++) {
		rotate->set(0,0,0);
		rotate->rotateZ((120*deg)*i);
		if(detectorNumber > 3) {
			rotate->rotateZ(180*deg);
		}
		rotate->rotateY(M_PI);
		rotate->rotateY(M_PI+beta);
		rotate->rotateZ(gamma);
		move = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,y,z,theta));
		fAssembly->MakeImprint(expHallLog, move, rotate, copyNumber);
	}

	if(hevimetopt == 1) {
		for(G4int i=0; i<3; i++) {
			rotate->set(0,0,0);
			rotate->rotateZ((120*deg)*i);
			if(detectorNumber > 3) {
				rotate->rotateZ(180*deg);
			}
			rotate->rotateY(M_PI);
			rotate->rotateY(M_PI+beta);
			rotate->rotateZ(gamma);
			move = G4ThreeVector(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,y,z,theta));
			fAssemblyHevimet->MakeImprint(expHallLog, move, rotate, copyNumber);
		}
	}

	return 1;
}

G4int DetectionSystemAncillaryBGO::BuildBGOVolume() {
	G4Material* material = G4Material::GetMaterial(fBgoMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fBgoMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.2,1.0,0.3));
	visAtt->SetVisibility(true);

	G4SubtractionSolid* bgoBlock = BGOPiece();

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= fCanThicknessFront + fGapThickness + (fBgoLength - fCanLength)/2.0;
	G4ThreeVector move 		= zPosition * direction;

	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fBgoBlockLog == nullptr) {
		fBgoBlockLog = new G4LogicalVolume(bgoBlock, material, "ancillaryBgoBlockLog", 0, 0, 0);
		fBgoBlockLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fBgoBlockLog, move, rotate);

	return 1;
}


G4int DetectionSystemAncillaryBGO::BuildVacuumVolume() {
	G4Material* material = G4Material::GetMaterial(fGapMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fGapMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	visAtt->SetVisibility(true);

	G4SubtractionSolid* vacuumBlock = VacuumPiece();

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= fCanThicknessFront + (fBgoLength + fGapThickness - fCanLength)/2.0;
	G4ThreeVector move 		= zPosition * direction;

	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fVacuumBlockLog == nullptr) {
		fVacuumBlockLog = new G4LogicalVolume(vacuumBlock, material, "ancillaryVacuumBlockLog", 0, 0, 0);
		fVacuumBlockLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fVacuumBlockLog, move, rotate);

	return 1;
}

G4int DetectionSystemAncillaryBGO::BuildAluminumCanVolume() {
	G4Material* material = G4Material::GetMaterial(fCanMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fCanMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	visAtt->SetVisibility(true);

	G4SubtractionSolid* canCylinder = AluminumCanPiece();

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= 0.0*mm;
	G4ThreeVector move 		= zPosition * direction;

	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fCanCylinderLog == nullptr) {
		fCanCylinderLog = new G4LogicalVolume(canCylinder, material, "ancillaryCanCylinderLog", 0, 0, 0);
		fCanCylinderLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fCanCylinderLog, move, rotate);

	return 1;
}

G4int DetectionSystemAncillaryBGO::BuildAluminumCanBackCoverVolume() {
	G4Material* material = G4Material::GetMaterial(fCanMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fCanMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	visAtt->SetVisibility(true);

	G4Tubs* canCylinder = AluminumCanBackCover();

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= (fCanLength - fCanThickness)/2.0;
	G4ThreeVector move 		= zPosition * direction;

	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fCanCylinderBackCoverLog == nullptr) {
		fCanCylinderBackCoverLog = new G4LogicalVolume(canCylinder, material, "ancillaryCanCylinderBackCoverLog", 0, 0, 0);
		fCanCylinderBackCoverLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fCanCylinderBackCoverLog, move, rotate);

	return 1;
}

G4int DetectionSystemAncillaryBGO::BuildHevimetPiece() {
	G4Material* material = G4Material::GetMaterial(fHevimetMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fHevimetMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	visAtt->SetVisibility(true);

	G4SubtractionSolid* canCylinder = HevimetPiece();

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= -1.0*(fHevimetThickness + fCanLength)/2.0;
	G4ThreeVector move 		= zPosition * direction;

	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fHevimetBlockLog == nullptr) {
		fHevimetBlockLog = new G4LogicalVolume(canCylinder, material, "ancillaryHevimetBlockLog", 0, 0, 0);
		fHevimetBlockLog->SetVisAttributes(visAtt);
	}

	fAssemblyHevimet->AddPlacedVolume(fHevimetBlockLog, move, rotate);

	return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4SubtractionSolid* DetectionSystemAncillaryBGO::BGOPiece() {
	G4double startPhi = 0.0*deg;
	G4double endPhi = 120.0*deg;

	G4double innerRadius = fCanInnerRadius+fCanThickness+fGapThickness;
	G4double outerRadius = innerRadius+fBgoThickness;
	G4double halfLengthZ = fBgoLength/2.0;

	G4Tubs* bgoBlock = new G4Tubs("bgoBlock", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// turn the block to be symmetric about 0 degress so we can cut it.
	bgoBlock->SetStartPhiAngle(-60.0*deg,true);
	bgoBlock->SetDeltaPhiAngle(120.0*deg);

	// grab the scissors, time to cut.
	// make the cut surface unnecessarily large.
	G4double cutHalfLengthZ = 2.0*fCanFaceToBuildOrigin;
	G4double cutHalfLengthX = 2.0*fCanFaceToBuildOrigin;
	G4double cutHalfLengthY = 2.0*fCanFaceToBuildOrigin;

	G4Box* cutPlate = new G4Box("cutPlate", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	// move cutting box to the "build origin" (see technical drawings)
	G4double moveCutZ = -1.0*(fCanFaceToBuildOrigin + fGapThickness + fCanThicknessFront + fBgoLength/2.0);
	G4ThreeVector moveCut((cutHalfLengthX-fGapThickness-fCanThickness)/(cos(fChoppingOuterAngle)),0,moveCutZ);

	// rotate our cutting block to the "chopping" angle
	G4RotationMatrix* rotateCut = new G4RotationMatrix;
	rotateCut->rotateY(-1.0*fChoppingOuterAngle);

	// snip snip
	G4SubtractionSolid* bgoBlockCut = new G4SubtractionSolid("bgoBlockCut", bgoBlock, cutPlate, rotateCut, moveCut);

	return bgoBlockCut;
}

G4SubtractionSolid* DetectionSystemAncillaryBGO::VacuumPiece() {
	G4double startPhi = 0.0*deg;
	G4double endPhi = 120.0*deg;

	G4double innerRadius = fCanInnerRadius+fCanThickness;
	G4double outerRadius = innerRadius+fGapThickness+fBgoThickness+fGapThicknessOuter;
	G4double halfLengthZ = (fBgoLength + fGapThickness)/2.0;

	G4Tubs* vacuumBlock = new G4Tubs("vacuumBlock", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// turn the block to be symmetric about 0 degress so we can cut it.
	vacuumBlock->SetStartPhiAngle(-60.0*deg,true);
	vacuumBlock->SetDeltaPhiAngle(120.0*deg);

	// grab the scissors, time to cut.
	// make the cut surface unnecessarily large.
	G4double cutHalfLengthZ = 2.0*fCanFaceToBuildOrigin;
	G4double cutHalfLengthX = 2.0*fCanFaceToBuildOrigin;
	G4double cutHalfLengthY = 2.0*fCanFaceToBuildOrigin;

	G4Box* cutPlate = new G4Box("cutPlate", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	// move cutting box to the "build origin" (see technical drawings)
	G4double moveCutZ = -1.0*(fCanFaceToBuildOrigin + fGapThickness + fCanThicknessFront + fBgoLength/2.0);
	G4ThreeVector moveCut((cutHalfLengthX-fCanThickness)/(cos(fChoppingOuterAngle)),0,moveCutZ);

	// rotate our cutting block to the "chopping" angle
	G4RotationMatrix* rotateCut = new G4RotationMatrix;
	rotateCut->rotateY(-1.0*fChoppingOuterAngle);

	// snip snip
	G4SubtractionSolid* vacuumBlockCut = new G4SubtractionSolid("vacuumBlockCut", vacuumBlock, cutPlate, rotateCut, moveCut);

	// Now we have to do another cut on the vacuum.
	// We need to remove the area where the BGO sits
	// let's take this piece larger in phi to make sure we cut everything
	startPhi = 0.0*deg;
	endPhi = 140.0*deg;

	innerRadius = fCanInnerRadius+fCanThickness+fGapThickness;
	outerRadius = innerRadius+fBgoThickness;
	halfLengthZ = (fBgoLength + fGapThickness)/2.0;

	G4Tubs* vacuumBlockCutBgo = new G4Tubs("vacuumBlockCutBgo", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// turn the block to be symmetric about 0 degress so we can cut it.
	vacuumBlockCutBgo->SetStartPhiAngle(-70.0*deg,true);
	vacuumBlockCutBgo->SetDeltaPhiAngle(140.0*deg);

	// grab the scissors, time to cut.
	// make the cut surface unnecessarily large.
	cutHalfLengthZ = 2.0*fCanFaceToBuildOrigin;
	cutHalfLengthX = 2.0*fCanFaceToBuildOrigin;
	cutHalfLengthY = 2.0*fCanFaceToBuildOrigin;

	G4Box* cutPlateForBgo = new G4Box("cutPlateForBgo", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	// move cutting box to the "build origin" (see technical drawings)
	moveCutZ = -1.0*(fCanFaceToBuildOrigin + fGapThickness + fCanThicknessFront + (fBgoLength + fGapThickness)/2.0);
	G4ThreeVector moveCutForBgo((cutHalfLengthX-fCanThickness-fGapThickness)/(cos(fChoppingOuterAngle)),0,moveCutZ);

	// rotate our cutting block to the "chopping" angle
	G4RotationMatrix* rotateCutForBgo = new G4RotationMatrix;
	rotateCutForBgo->rotateY(-1.0*fChoppingOuterAngle);

	// snip snip
	G4SubtractionSolid* bgoBlockToBeCutOut = new G4SubtractionSolid("bgoBlockToBeCutOut", vacuumBlockCutBgo, cutPlateForBgo, rotateCutForBgo, moveCutForBgo);

	// now we can remove the BGO area from inside the vacuum
	// move cut just along z-axis by the amount of the gapThickness (both vacuum and bgo blacks are same size)
	G4ThreeVector moveCutMe(0, 0, fGapThickness);
	// no rotation
	G4RotationMatrix* rotateCutMe = new G4RotationMatrix;
	// final snip snip
	G4SubtractionSolid* vacuumBlockWithBgoCutOut = new G4SubtractionSolid("vacuumBlockWithBgoCutOut", vacuumBlockCut, bgoBlockToBeCutOut, rotateCutMe, moveCutMe);

	return vacuumBlockWithBgoCutOut;
}

G4SubtractionSolid* DetectionSystemAncillaryBGO::AluminumCanPiece() {
	G4double startPhi = 0.0*deg;
	G4double endPhi = 120.0*deg;

	G4double innerRadius = fCanInnerRadius;
	G4double outerRadius = innerRadius+fCanThickness+fGapThickness+fBgoThickness+fGapThicknessOuter+fCanThickness; //2*fCanThickness ???
	G4double halfLengthZ = (fCanLength)/2.0;

	G4Tubs* canBlock = new G4Tubs("canBlock", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// turn the block to be symmetric about 0 degress so we can cut it.
	canBlock->SetStartPhiAngle(-60.0*deg,true);
	canBlock->SetDeltaPhiAngle(120.0*deg);

	// grab the scissors, time to cut.
	// make the cut surface unnecessarily large.
	G4double cutHalfLengthZ = 2.0*fCanFaceToBuildOrigin;
	G4double cutHalfLengthX = 2.0*fCanFaceToBuildOrigin;
	G4double cutHalfLengthY = 2.0*fCanFaceToBuildOrigin;

	G4Box* cutPlate = new G4Box("cutPlate", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	// move cutting box to the "build origin" (see technical drawings)
	G4double moveCutZ = -1.0*(fCanFaceToBuildOrigin + fCanLength/2.0);
	G4ThreeVector moveCut((cutHalfLengthX)/(cos(fChoppingOuterAngle)),0,moveCutZ);

	// rotate our cutting block to the "chopping" angle
	G4RotationMatrix* rotateCut = new G4RotationMatrix;
	rotateCut->rotateY(-1.0*fChoppingOuterAngle);

	// snip snip
	G4SubtractionSolid* canBlockCut = new G4SubtractionSolid("canBlockCut", canBlock, cutPlate, rotateCut, moveCut);

	// Now we have to do another cut on the can.
	// We need to remove the area where the vacuum and BGO sits
	// let's take this piece larger in phi to make sure we cut everything
	startPhi = 0.0*deg;
	endPhi = 140.0*deg;

	innerRadius = fCanInnerRadius+fCanThickness;
	outerRadius = innerRadius+fGapThickness+fBgoThickness+fGapThicknessOuter;
	halfLengthZ = (fCanLength)/2.0;

	G4Tubs* canBlockCutVacuum = new G4Tubs("canBlockCutVacuum", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// turn the block to be symmetric about 0 degress so we can cut it.
	canBlockCutVacuum->SetStartPhiAngle(-70.0*deg,true);
	canBlockCutVacuum->SetDeltaPhiAngle(140.0*deg);

	// grab the scissors, time to cut.
	// make the cut surface unnecessarily large.
	cutHalfLengthZ = 2.0*fCanFaceToBuildOrigin;
	cutHalfLengthX = 2.0*fCanFaceToBuildOrigin;
	cutHalfLengthY = 2.0*fCanFaceToBuildOrigin;

	G4Box* cutPlate3 = new G4Box("cutPlate3", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	// move cutting box to the "build origin" (see technical drawings)
	moveCutZ = -1.0*(fCanFaceToBuildOrigin + fCanThicknessFront + fCanLength/2.0);
	G4ThreeVector moveCut3((cutHalfLengthX-fCanThickness)/(cos(fChoppingOuterAngle)),0,moveCutZ);

	// rotate our cutting block to the "chopping" angle
	G4RotationMatrix* rotateCut3 = new G4RotationMatrix;
	rotateCut3->rotateY(-1.0*fChoppingOuterAngle);

	// snip snip
	G4SubtractionSolid* vacuumBlockToBeCutOut = new G4SubtractionSolid("vacuumBlockToBeCutOut", canBlockCutVacuum, cutPlate3, rotateCut3, moveCut3);

	// now we can remove the vacuum/BGO area from inside the can
	// move cut just along z-axis by the amount of the canThicknessFront (both vacuum and can blacks are same size)
	G4ThreeVector moveCutMe(0,0,fCanThicknessFront);
	// no rotation
	G4RotationMatrix* rotateCutMe = new G4RotationMatrix;
	// final snip snip
	G4SubtractionSolid* canBlockWithVacuumCutOut = new G4SubtractionSolid("canBlockCut", canBlockCut, vacuumBlockToBeCutOut, rotateCutMe, moveCutMe);

	return canBlockWithVacuumCutOut;
}


G4Tubs* DetectionSystemAncillaryBGO::AluminumCanBackCover() {
	G4double startPhi = 0.0*deg;
	G4double endPhi = 120.0*deg;

	G4double innerRadius = fCanInnerRadius+fCanThickness;
	G4double outerRadius = innerRadius+fCanThickness+fGapThickness+fBgoThickness+fGapThicknessOuter;
	G4double halfLengthZ = (fCanThickness)/2.0;

	G4Tubs* canBlock = new G4Tubs("canBlock", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// turn the block to be symmetric about 0 degress
	canBlock->SetStartPhiAngle(-60.0*deg,true);
	canBlock->SetDeltaPhiAngle(120.0*deg);

	return canBlock;
}

G4SubtractionSolid* DetectionSystemAncillaryBGO::HevimetPiece() {// this has been fix 
	G4double startPhi = 0.0*deg;
	G4double endPhi = 120.0*deg;

	G4double innerRadius = 0.0*mm;
	G4double outerRadius = fCanInnerRadius+fCanThickness+fGapThickness+fBgoThickness+fGapThicknessOuter+fCanThickness; //2*fCanThickness ???
	G4double halfLengthZ = (fHevimetThickness)/2.0;

	G4Tubs* hevimetBlock = new G4Tubs("hevimetBlock", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// turn the block to be symmetric about 0 degress so we can cut it.
	hevimetBlock->SetStartPhiAngle(-60.0*deg,true);
	hevimetBlock->SetDeltaPhiAngle(120.0*deg);

	// grab the scissors, time to cut.
	// make the cut surface unnecessarily large.
	G4double cutHalfLengthZ = 2.0*fCanFaceToBuildOrigin;
	G4double cutHalfLengthX = 2.0*fCanFaceToBuildOrigin;
	G4double cutHalfLengthY = 2.0*fCanFaceToBuildOrigin;

	G4Box* cutPlate = new G4Box("cutPlate", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	// move cutting box to the "build origin" (see technical drawings)
	G4double moveCutZ = -1.0*(fCanFaceToBuildOrigin)+fHevimetThickness/2.0;
	G4ThreeVector moveCut((cutHalfLengthX)/(cos(fChoppingOuterAngle)),0,moveCutZ);

	// rotate our cutting block to the "chopping" angle
	G4RotationMatrix* rotateCut = new G4RotationMatrix;
	rotateCut->rotateY(-1.0*fChoppingOuterAngle);

	// snip snip
	G4SubtractionSolid* hevimetBlockCut = new G4SubtractionSolid("hevimetBlockCut", hevimetBlock, cutPlate, rotateCut, moveCut);


	// Now we need to cone in the middle.
	G4double coneInnerRadius = 0.0*mm;
	G4double coneLowerOuterRadius = 0.0*mm;
	G4double coneUpperOuterRadius = tan(fChoppingInnerAngle)*(fCanFaceToBuildOrigin+fHevimetThickness);
	G4double coneHalfLengthZ = (fCanFaceToBuildOrigin+fHevimetThickness)/2.0;

	G4Cons* hevimetCutCone = new G4Cons("hevimetCutCone", coneInnerRadius, coneLowerOuterRadius, coneInnerRadius, coneUpperOuterRadius, coneHalfLengthZ, 0.0*deg, 360.0*deg);

	// now we can remove the inner radius from the Hevimet
	// move cut just along z-axis by the amount of the canThicknessFront (both vacuum and can blacks are same size)
	G4ThreeVector moveCutMe(0,0,-1.0*(fCanFaceToBuildOrigin)/2.0+fHevimetThickness/2.0);
	// no rotation
	G4RotationMatrix* rotateCutMe = new G4RotationMatrix;
	// final snip snip
	G4SubtractionSolid* hevimetBlockCutWithConeCut = new G4SubtractionSolid("hevimetBlockCutWithConeCut", hevimetBlockCut, hevimetCutCone, rotateCutMe, moveCutMe);


	return hevimetBlockCutWithConeCut;
}


//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemAncillaryBGO::GetDirectionXYZ(G4double theta, G4double phi) {
	G4double x,y,z;
	x = sin(theta) * cos(phi);
	y = sin(theta) * sin(phi);
	z = cos(theta);

	G4ThreeVector direction = G4ThreeVector(x,y,z);

	return direction;
}//end ::GetDirection
