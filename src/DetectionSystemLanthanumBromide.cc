#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "Global.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemLanthanumBromide.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


DetectionSystemLanthanumBromide::DetectionSystemLanthanumBromide() :
	// LogicalVolumes
	fDetectorVolumeLog(0),
	fCrystalBlockLog(0),
	fCanCylinderLog(0),
	fCanFrontLidLog(0),
	fCanBackLidLog(0),
	fSealFrontLidLog(0),
	fDiscFrontLidLog(0),
	fPackingCylinderLog(0),
	fPackingFrontLidLog(0)
{

	// NaI Crystal and Can Physical Properties
	// From: Scintillation Spectrometry gamma-ray catalogue - R. L. Heath, 1964

	/////////////////////////////////////////////////////////////////////
	//  0.019"
	//  "can" - Al Container
	//
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% 0.040"
	//% "seal" - Compressed Neoprene Rubber
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//-------------------------------------------------------------------
	//- 0.006"
	//- "disc" - Polyethylene Disc
	//-------------------------------------------------------------------
	/////////////////////////////////////////////////////////////////////
	//- 0.0625"
	//- "packing" - Packed Aluminum Oxide
	/////////////////////////////////////////////////////////////////////
	//
	//
	//  NaI Crystal - Used as first-order estimate of Lanthanum Bromide Crystal
	//
	//
	//
	//
	//


	fDetailViewEndAngle	      = 360.0*deg;
	fCrystalMaterial            = "CeriumDopedLanthanumBromide";

	fCanMaterial                = "G4_Al";
	fSealMaterial               = "G4_RUBBER_NEOPRENE";
	fDiscMaterial               = "G4_POLYETHYLENE";
	fPackingMaterial            = "G4_ALUMINUM_OXIDE";

	fCrystalLengthZ            = 2.0*inchtocm*cm;
	fCrystalInnerRadius 		  = 0.0*cm;
	fCrystalOuterRadius 		  = 1.0*inchtocm*cm;

	fPackingLengthZ            = fCrystalLengthZ;
	fPackingInnerRadius 		  = fCrystalOuterRadius;
	fPackingOuterRadius 		  = 0.0625*inchtocm*cm + fCrystalOuterRadius;

	fPackingLidInnerRadius    = 0.0*cm;
	fPackingLidOuterRadius 	  = fPackingOuterRadius;
	fPackingFrontLidThickness = 0.0625*inchtocm*cm;

	fDiscLidInnerRadius       = 0.0*cm;
	fDiscLidOuterRadius       = fPackingLidOuterRadius;
	fDiscFrontLidThickness	= 0.006*inchtocm*cm;

	fSealLidInnerRadius       = 0.0*cm;
	fSealLidOuterRadius       = fPackingLidOuterRadius;
	fSealFrontLidThickness	  = 0.04*inchtocm*cm;

	fCanLengthZ                = fCrystalLengthZ + fPackingFrontLidThickness + fDiscFrontLidThickness + fSealFrontLidThickness;
	fCanInnerRadius            = fPackingOuterRadius;
	fCanOuterRadius            = 0.019*inchtocm*cm + fPackingOuterRadius;

	fCanLidInnerRadius        = 0.0*cm;
	fCanLidOuterRadius        = fCanOuterRadius;
	fCanFrontLidThickness     = 0.019*inchtocm*cm;
	fCanBackLidThickness      = 0.019*inchtocm*cm;


	fDetectorLengthZ           = fCrystalLengthZ +
		fPackingFrontLidThickness +
		fDiscFrontLidThickness +
		fSealFrontLidThickness +
		fCanFrontLidThickness +
		fCanBackLidThickness;
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

DetectionSystemLanthanumBromide::~DetectionSystemLanthanumBromide() {
	// LogicalVolumes
	delete fDetectorVolumeLog;
	delete fCrystalBlockLog;
	delete fCanCylinderLog;
	delete fCanFrontLidLog;
	delete fCanBackLidLog;
	delete fSealFrontLidLog;
	delete fDiscFrontLidLog;
	delete fPackingCylinderLog;
	delete fPackingFrontLidLog;
}


G4int DetectionSystemLanthanumBromide::Build() {

	// Build assembly volume
	fAssembly = new G4AssemblyVolume();

	G4cout<<"BuildCrystalVolume"<<G4endl;
	BuildCrystalVolume();
	G4cout<<"BuildAluminumCanVolume"<<G4endl;
	BuildAluminumCanVolume();
	G4cout<<"BuildPackingVolume"<<G4endl;
	BuildPackingVolume();
	G4cout<<"BuildDiscVolume"<<G4endl;
	BuildDiscVolume();
	G4cout<<"BuildSealVolume"<<G4endl;
	BuildSealVolume();

	return 1;
}

G4double DetectionSystemLanthanumBromide::GetR() {
	// to crystal face

	G4double position = fSetRadialPos - (fPackingFrontLidThickness + fDiscFrontLidThickness + fCanFrontLidThickness - fCanBackLidThickness);
	return position;
}

G4double DetectionSystemLanthanumBromide::GetTheta(G4int i) {
	// to crystal face

	return fDetectorAngles[i][0]; //theta
}

G4double DetectionSystemLanthanumBromide::GetPhi(G4int i) {
	// to crystal face

	return fDetectorAngles[i][1]; //phi
}

G4double DetectionSystemLanthanumBromide::GetYaw(G4int i) {
	// to crystal face

	return fDetectorAngles[i][2]; //yaw
}

G4double DetectionSystemLanthanumBromide::GetPitch(G4int i) {
	// to crystal face

	return fDetectorAngles[i][3]; //pitch
}

G4double DetectionSystemLanthanumBromide::GetRoll(G4int i) {
	// to crystal face

	return fDetectorAngles[i][4]; //roll
}

G4int DetectionSystemLanthanumBromide::PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber, G4double radialpos) {
	G4int detectorCopyID = 0;

	G4cout<<"LanthanumBromide Detector Number = "<<detectorNumber<<G4endl;

	fCopyNumber = detectorCopyID + detectorNumber;

	G4double position = radialpos + fCanFrontLidThickness +(fCanLengthZ / 2.0) ;
	fSetRadialPos = radialpos;

	G4double theta  = fDetectorAngles[detectorNumber][0];
	G4double phi    = fDetectorAngles[detectorNumber][1];
	//G4double alpha  = fDetectorAngles[detectorNumber][2]; // yaw
	G4double beta   = fDetectorAngles[detectorNumber][3]; // pitch
	G4double gamma  = fDetectorAngles[detectorNumber][4]; // roll

	G4double x = 0;
	G4double y = 0;
	G4double z = position;

	G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
	rotate->rotateY(M_PI);
	rotate->rotateY(M_PI+beta);
	rotate->rotateZ(gamma);

	G4ThreeVector move(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,y,z,theta));

	fAssembly->MakeImprint(expHallLog, move, rotate, fCopyNumber);

	return 1;
}

G4int DetectionSystemLanthanumBromide::BuildCrystalVolume() {
	G4Material* material = G4Material::GetMaterial(fCrystalMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fCrystalMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.2,1.0,0.3));
	visAtt->SetVisibility(true);

	G4Tubs* crystalBlock = BuildCrystal();

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= ((fCanFrontLidThickness + fSealFrontLidThickness + fDiscFrontLidThickness + fPackingFrontLidThickness) - (fCanBackLidThickness))/2.0;
	G4ThreeVector move 		= zPosition * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;

	//logical volume
	if(fCrystalBlockLog == nullptr) {
		fCrystalBlockLog = new G4LogicalVolume(crystalBlock, material, "lanthanumBromideCrystalBlockLog", 0, 0, 0);
		fCrystalBlockLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fCrystalBlockLog, move, rotate);

	return 1;
}

G4int DetectionSystemLanthanumBromide::BuildAluminumCanVolume() {
	G4Material* material = G4Material::GetMaterial(fCanMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fCanMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.55, 0.38, 1.0));
	visAtt->SetVisibility(true);

	G4ThreeVector direction =  G4ThreeVector(0,0,1);
	G4double zPosition;
	G4ThreeVector move;
	G4RotationMatrix* rotate = new G4RotationMatrix;

	/////////////////////////////////////////////////////////////////////
	// Build and Place Aluminum Can
	/////////////////////////////////////////////////////////////////////
	G4Tubs* canCylinder = BuildAluminumCan();

	//logical volume
	if(fCanCylinderLog == nullptr) {
		fCanCylinderLog = new G4LogicalVolume(canCylinder, material, "canCylinderLog", 0, 0, 0);
		fCanCylinderLog->SetVisAttributes(visAtt);
	}

	// place front canLid
	zPosition 	= 0;
	move 		= zPosition * direction;

	//add physical cylinder
	fAssembly->AddPlacedVolume(fCanCylinderLog, move, rotate);

	/////////////////////////////////////////////////////////////////////
	// Build and Place Aluminum Front Lid
	/////////////////////////////////////////////////////////////////////
	G4Tubs* canFrontLid = BuildAluminumCanFrontLid();

	// logical volume
	if(fCanFrontLidLog == nullptr) {
		fCanFrontLidLog = new G4LogicalVolume(canFrontLid, material, "canFrontLidLog", 0, 0, 0);
		fCanFrontLidLog->SetVisAttributes(visAtt);
	}

	// place front canLid

	zPosition 	= -(fDetectorLengthZ/2.0) + (fCanFrontLidThickness/2.0);
	move 		= zPosition * direction;

	//add physical front canLid
	fAssembly->AddPlacedVolume(fCanFrontLidLog, move, rotate);

	/////////////////////////////////////////////////////////////////////
	// Build and Place Aluminum Back Lid
	/////////////////////////////////////////////////////////////////////
	G4Tubs* canBackLid = BuildAluminumCanBackLid();

	// logical volume
	if(fCanBackLidLog == nullptr) {
		fCanBackLidLog = new G4LogicalVolume(canBackLid, material, "canBackLidLog", 0, 0, 0);
		fCanBackLidLog->SetVisAttributes(visAtt);
	}

	// place back canLid
	zPosition 	= (fDetectorLengthZ/2.0) - (fCanBackLidThickness/2.0);
	move 		= zPosition * direction;

	// add physical back canLid
	fAssembly->AddPlacedVolume(fCanBackLidLog, move, rotate);

	return 1;
}

G4int DetectionSystemLanthanumBromide::BuildPackingVolume() {
	G4Material* material = G4Material::GetMaterial(fPackingMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fPackingMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	visAtt->SetVisibility(true);

	G4ThreeVector direction =  G4ThreeVector(0,0,1);
	G4double zPosition;
	G4ThreeVector move;
	G4RotationMatrix* rotate = new G4RotationMatrix;

	/////////////////////////////////////////////////////////////////////
	// Build and Place Packing
	/////////////////////////////////////////////////////////////////////
	G4Tubs* packingCylinder = BuildPacking();

	//logical volume
	if(fPackingCylinderLog == nullptr) {
		fPackingCylinderLog = new G4LogicalVolume(packingCylinder, material, "packingCylinderLog", 0, 0, 0);
		fPackingCylinderLog->SetVisAttributes(visAtt);
	}

	// place cylinder
	zPosition 	= ((fCanFrontLidThickness + fSealFrontLidThickness + fDiscFrontLidThickness + fPackingFrontLidThickness) - (fCanBackLidThickness))/2.0;
	move          = zPosition * direction;

	//add physical cylinder
	fAssembly->AddPlacedVolume(fPackingCylinderLog, move, rotate);

	/////////////////////////////////////////////////////////////////////
	// Build and Place Front Lid
	/////////////////////////////////////////////////////////////////////
	G4Tubs* packingFrontLid = BuildPackingFrontLid();

	// logical volume
	if(fPackingFrontLidLog == nullptr) {
		fPackingFrontLidLog = new G4LogicalVolume(packingFrontLid, material, "packingFrontLidLog", 0, 0, 0);
		fPackingFrontLidLog->SetVisAttributes(visAtt);
	}

	// place front packingLid
	zPosition 	= -(fDetectorLengthZ/2.0) + (fCanFrontLidThickness + fSealFrontLidThickness + fDiscFrontLidThickness + (fPackingFrontLidThickness/2.0));
	move 		= zPosition * direction;

	//add physical front packingLid
	fAssembly->AddPlacedVolume(fPackingFrontLidLog, move, rotate);

	return 1;
}

G4int DetectionSystemLanthanumBromide::BuildDiscVolume() {
	G4Material* material = G4Material::GetMaterial(fDiscMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fDiscMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	visAtt->SetVisibility(true);

	G4ThreeVector direction =  G4ThreeVector(0,0,1);
	G4double zPosition;
	G4ThreeVector move;
	G4RotationMatrix* rotate = new G4RotationMatrix;

	/////////////////////////////////////////////////////////////////////
	// Build and Place Front Lid
	/////////////////////////////////////////////////////////////////////
	G4Tubs* discFrontLid = BuildDiscFrontLid();

	// logical volume
	if(fDiscFrontLidLog == nullptr) {
		fDiscFrontLidLog = new G4LogicalVolume(discFrontLid, material, "discFrontLidLog", 0, 0, 0);
		fDiscFrontLidLog->SetVisAttributes(visAtt);
	}

	// place front discLid
	zPosition 	= -(fDetectorLengthZ/2.0) + (fCanFrontLidThickness + fSealFrontLidThickness + (fDiscFrontLidThickness/2.0));
	move 		= zPosition * direction;

	//add physical front discLid
	fAssembly->AddPlacedVolume(fDiscFrontLidLog, move, rotate);

	return 1;
}

G4int DetectionSystemLanthanumBromide::BuildSealVolume() {
	G4Material* material = G4Material::GetMaterial(fSealMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fSealMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
	visAtt->SetVisibility(true);

	G4ThreeVector direction =  G4ThreeVector(0,0,1);
	G4double zPosition;
	G4ThreeVector move;
	G4RotationMatrix* rotate = new G4RotationMatrix;

	/////////////////////////////////////////////////////////////////////
	// Build and Place Front Lid
	/////////////////////////////////////////////////////////////////////
	G4Tubs* sealFrontLid = BuildSealFrontLid();

	// logical volume
	if(fSealFrontLidLog == nullptr) {
		fSealFrontLidLog = new G4LogicalVolume(sealFrontLid, material, "sealFrontLidLog", 0, 0, 0);
		fSealFrontLidLog->SetVisAttributes(visAtt);
	}

	// place front sealLid
	zPosition 	= -(fDetectorLengthZ/2.0) + (fCanFrontLidThickness + (fSealFrontLidThickness/2.0));
	move 		= zPosition * direction;

	//add physical front sealLid
	fAssembly->AddPlacedVolume(fSealFrontLidLog, move, rotate);

	return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Tubs* DetectionSystemLanthanumBromide::BuildCrystal() {
	G4double startPhi = 0.0;
	G4double endPhi = fDetailViewEndAngle;

	G4double innerRadius = fCrystalInnerRadius;
	G4double outerRadius = fCrystalOuterRadius;
	G4double halfLengthZ = (fCrystalLengthZ)/2.0;

	G4Tubs* crystalBlock = new G4Tubs("crystalBlock", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return crystalBlock;
}//end ::BuildCrystal


G4Tubs* DetectionSystemLanthanumBromide::BuildAluminumCan() {
	G4double startPhi = 0.0;
	G4double endPhi = fDetailViewEndAngle;

	G4double innerRadius 	= fCanInnerRadius;
	G4double outerRadius 	= fCanOuterRadius;
	G4double halfLengthZ 	= fCanLengthZ/2.0;

	G4Tubs* canCylinder = new G4Tubs("canCylinder", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return canCylinder;
}//end ::BuildAluminumCan

G4Tubs* DetectionSystemLanthanumBromide::BuildAluminumCanFrontLid() {
	G4double startPhi = 0.0;
	G4double endPhi = fDetailViewEndAngle;

	G4double innerRadius = fCanLidInnerRadius;
	G4double outerRadius = fCanLidOuterRadius;
	G4double halfLengthZ = fCanFrontLidThickness/2.0;

	G4Tubs* canLid = new G4Tubs("canLid", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return canLid;
}//end ::BuildAluminumFrontLid

G4Tubs* DetectionSystemLanthanumBromide::BuildAluminumCanBackLid() {
	G4double startPhi = 0.0;
	G4double endPhi = fDetailViewEndAngle;

	G4double innerRadius = fCanLidInnerRadius;
	G4double outerRadius = fCanLidOuterRadius;
	G4double halfLengthZ = fCanBackLidThickness/2.0;

	G4Tubs* canLid = new G4Tubs("canLid", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return canLid;
}//end ::BuildAluminumBackLid

G4Tubs* DetectionSystemLanthanumBromide::BuildPacking() {
	G4double startPhi = 0.0;
	G4double endPhi = fDetailViewEndAngle;

	G4double innerRadius 	= fPackingInnerRadius;
	G4double outerRadius 	= fPackingOuterRadius;
	G4double halfLengthZ 	= fPackingLengthZ/2.0;

	G4Tubs* packingCylinder = new G4Tubs("packingCylinder", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return packingCylinder;
}//end ::BuildPacking

G4Tubs* DetectionSystemLanthanumBromide::BuildPackingFrontLid() {
	G4double startPhi = 0.0;
	G4double endPhi = fDetailViewEndAngle;

	G4double innerRadius = fPackingLidInnerRadius;
	G4double outerRadius = fPackingLidOuterRadius;
	G4double halfLengthZ = fPackingFrontLidThickness/2.0;

	G4Tubs* packingLid = new G4Tubs("packingLid", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return packingLid;
}//end ::BuildPackingFrontLid

G4Tubs* DetectionSystemLanthanumBromide::BuildDiscFrontLid() {
	G4double startPhi = 0.0;
	G4double endPhi = fDetailViewEndAngle;

	G4double innerRadius = fDiscLidInnerRadius;
	G4double outerRadius = fDiscLidOuterRadius;
	G4double halfLengthZ = fDiscFrontLidThickness/2.0;

	G4Tubs* discLid = new G4Tubs("discLid", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return discLid;
}//end ::BuildDiscFrontLid

G4Tubs* DetectionSystemLanthanumBromide::BuildSealFrontLid() {
	G4double startPhi = 0.0;
	G4double endPhi = fDetailViewEndAngle;

	G4double innerRadius = fSealLidInnerRadius;
	G4double outerRadius = fSealLidOuterRadius;
	G4double halfLengthZ = fSealFrontLidThickness/2.0;

	G4Tubs* sealLid = new G4Tubs("sealLid", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	return sealLid;
}//end ::BuildSealFrontLid

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemLanthanumBromide::GetDirectionXYZ(G4double theta, G4double phi) {
	G4double x,y,z;
	x = sin(theta) * cos(phi);
	y = sin(theta) * sin(phi);
	z = cos(theta);

	G4ThreeVector direction = G4ThreeVector(x,y,z);

	return direction;
}//end ::GetDirection
