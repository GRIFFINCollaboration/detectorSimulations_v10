#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
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

#include "DetectionSystemS3.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

DetectionSystemS3::DetectionSystemS3() :
	// LogicalVolumes
	fS3InnerGuardRingLog(0),
	fS3OuterGuardRingLog(0)
{
	/////////////////////////////////////////////////////////////////////
	// SPICE Physical Properties
	/////////////////////////////////////////////////////////////////////

	fWaferMaterial             = "Silicon";

	//-----------------------------//
	// parameters for the annular  //
	// planar detector crystal     //
	//-----------------------------//
	fS3DetCrystalOuterDiameter = 70.*mm;
	fS3DetCrystalInnerDiameter = 22.*mm;
	fS3DetCrystalThickness = .15*mm;
	fS3DetRadialSegments = 24.;
	fS3DetPhiSegments = 32.;

	//-----------------------------//
	// parameters for guard ring   //
	//-----------------------------//
	fS3DetGuardRingOuterDiameter = 76*mm;
	fS3DetGuardRingInnerDiameter = 20*mm;
}

DetectionSystemS3::~DetectionSystemS3() {
	for(int i = 0; i < 24; ++i) {
		if(fSiDetS3RingLog[i] != nullptr) {
			delete fSiDetS3RingLog[i];
		}
	}
	delete fS3InnerGuardRingLog;
	delete fS3OuterGuardRingLog;

}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemS3::Build() {
	fAssembly =  new G4AssemblyVolume();

	for(int ringID=0; ringID<24; ringID++) {	// Loop through each ring ...
		fAssemblyS3Ring[ringID] = new G4AssemblyVolume(); 		// Build assembly volumes
		BuildSiliconWafer(ringID);		// Build Silicon Ring
	} // end for(int ringID)

	// Build Guard Rings
	BuildInnerGuardRing();
	BuildOuterGuardRing();

	return 1;
}

//---------------------------------------------------------//
// "place" function called in DetectorMessenger            //
// if detector is added                                    //
//---------------------------------------------------------//
G4int DetectionSystemS3::PlaceDetector(G4LogicalVolume* expHallLog, G4ThreeVector move, G4int ringNumber, G4int Seg, G4int detectorNumber) {
	G4RotationMatrix* rotate = new G4RotationMatrix;
	G4int NumberSeg = (G4int) fS3DetPhiSegments;
	G4double angle = (360./NumberSeg)*(Seg-0.5)*deg;
	rotate->rotateZ(angle);
	fAssemblyS3Ring[ringNumber]->MakeImprint(expHallLog, move, rotate, detectorNumber);

	return 1;
}

G4int DetectionSystemS3::PlaceGuardRing(G4LogicalVolume* expHallLog, G4ThreeVector move) {
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);
	fAssembly->MakeImprint(expHallLog, move, rotate, 0);

	return 1;
}

//---------------------------------------------------------//
// build functions for different parts                     //
// called in main build function                           //
//---------------------------------------------------------//
G4int DetectionSystemS3::BuildSiliconWafer(G4int RingID) {
	// Define the material, return error if not found
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fWaferMaterial
			<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	visAtt->SetVisibility(true);

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= -(fS3DetCrystalThickness/2.);
	G4ThreeVector move 		= zPosition * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);

	// construct solid
	G4Tubs* siDetS3RingSec = BuildCrystal(RingID);

	// construct logical volume if it doesn't already exist
	if(!fSiDetS3RingLog[RingID]) {
		G4String s3name = "siDetS3Ring_";
		s3name += G4UIcommand::ConvertToString(RingID);
		s3name += "_Log";

		fSiDetS3RingLog[RingID] = new G4LogicalVolume(siDetS3RingSec, material, s3name, 0, 0, 0);
		fSiDetS3RingLog[RingID]->SetVisAttributes(visAtt);
	}
	fAssemblyS3Ring[RingID]->AddPlacedVolume(fSiDetS3RingLog[RingID], move, rotate);

	return 1;
}

G4int DetectionSystemS3::BuildInnerGuardRing() {
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fWaferMaterial<<" not found, cannot build the inner guard ring of the S3 detector! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	visAtt->SetVisibility(true);

	G4Tubs* innerGuardRing = new G4Tubs("innerGuardRing",
			fS3DetGuardRingInnerDiameter/2.,
			fS3DetCrystalInnerDiameter/2.,
			fS3DetCrystalThickness/2.,0,360);

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= -(fS3DetCrystalThickness/2.);
	G4ThreeVector move 		= zPosition * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);

	//logical volume
	if(fS3InnerGuardRingLog == nullptr) {
		fS3InnerGuardRingLog = new G4LogicalVolume(innerGuardRing, material, "innerGuardRing", 0,0,0);
		fS3InnerGuardRingLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fS3InnerGuardRingLog, move, rotate);

	return 1;
}

G4int DetectionSystemS3::BuildOuterGuardRing() {
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fWaferMaterial<<" not found, cannot build the outer guard ring of the S3 detector! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	visAtt->SetVisibility(true);

	G4Tubs* outerGuardRing = new G4Tubs("outerGuardRing",
			fS3DetCrystalOuterDiameter/2.,
			fS3DetGuardRingOuterDiameter/2.,
			fS3DetCrystalThickness/2.,0,360);

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= -(fS3DetCrystalThickness/2.);
	G4ThreeVector move 		= zPosition * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);

	//logical volume
	if(fS3OuterGuardRingLog == nullptr) {
		fS3OuterGuardRingLog = new G4LogicalVolume(outerGuardRing, material, "outerGuardRing", 0,0,0);
		fS3OuterGuardRingLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fS3OuterGuardRingLog, move, rotate);

	return 1;
}

///////////////////////////////////////////////////////
// Build one segment of S3
// the geometry depends on the distance from the center
///////////////////////////////////////////////////////
G4Tubs* DetectionSystemS3::BuildCrystal(G4int RingID) {
	// define angle, length, thickness, and inner and outer diameter
	// of silicon detector segment
	G4double tubeElementLength = (fS3DetCrystalOuterDiameter - fS3DetCrystalInnerDiameter)/(2*(fS3DetRadialSegments));
	G4double tubeElementAngularWidth = (360./fS3DetPhiSegments)*deg;
	G4double tubeElementInnerRadius = (fS3DetCrystalInnerDiameter)/2.0 + tubeElementLength*(RingID);
	G4double tubeElementOuterRadius = ((G4double)fS3DetCrystalInnerDiameter)/2.0 + tubeElementLength*(RingID+1);
	G4double tubeElementHalfThickness = (fS3DetCrystalThickness)/2.0;

	// establish solid
	G4Tubs* crystalBlock = new G4Tubs("crystalBlock",tubeElementInnerRadius,tubeElementOuterRadius,tubeElementHalfThickness,0,tubeElementAngularWidth);

	return crystalBlock;
}//end ::BuildCrystal

