#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

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

#include "DetectionSystemTestcan.hh"

#include "G4SystemOfUnits.hh"

#include <string>

DetectionSystemTestcan::DetectionSystemTestcan(G4double length, G4double radius) :
	// LogicalVolumes
	fTestcanAlumCasingLog(0),
	fTestcanScintillatorLog(0),
	fTestcanQuartzWindowLog(0)
{
	// can properties
	fScintillatorLength        = length;
	fAlumCanThickness          = 1.0*mm;
	fScintillatorInnerRadius   = 0.0*mm;
	fScintillatorOuterRadius   = radius;
	fQuartzThickness           = 6.35*mm;
	fQuartzRadius              = radius + fAlumCanThickness;
	fCanMaterial               = "G4_Al";
	fLiquidMaterial            = "BC537";
	fQuartzMaterial            = "G4_SILICON_DIOXIDE";

	fStartPhi               = 0.0*deg;
	fEndPhi                 = 360.0*deg;

	fLiquidColour           = G4Colour(0.0/255.0,255.0/255.0,225.0/255.0);
	fGreyColour             = G4Colour(0.5, 0.5, 0.5); 
	fQuartzColour           = G4Colour(1.0, 0.0, 1.0);    

}


DetectionSystemTestcan::~DetectionSystemTestcan() {
	// LogicalVolumes
	delete fTestcanAlumCasingLog;
	delete fTestcanScintillatorLog;
	delete fTestcanQuartzWindowLog;
}


G4int DetectionSystemTestcan::Build() {
	fAssemblyTestcan = new G4AssemblyVolume(); 

	BuildTestcan();

	return 1;
}

G4int DetectionSystemTestcan::PlaceDetector(G4LogicalVolume* expHallLog) {
	G4RotationMatrix * rotate = new G4RotationMatrix;
	G4ThreeVector move = G4ThreeVector(0., 0., 0.);

	fAssemblyTestcan->MakeImprint(expHallLog, move, rotate);

	return 1;
}

G4int DetectionSystemTestcan::BuildTestcan() {
	G4ThreeVector move, direction;
	G4RotationMatrix* rotate;

	G4Material* canG4material = G4Material::GetMaterial(fCanMaterial);
	if(!canG4material) {
		G4cout<<" ----> Material "<<fCanMaterial<<" not found, cannot build! "<<G4endl;
		return 0;
	}

	G4Material* liquidG4material = G4Material::GetMaterial(fLiquidMaterial);
	if(!liquidG4material) {
		G4cout<<" ----> Material "<<fLiquidMaterial<<" not found, cannot build! "<<G4endl;
		return 0;
	}

	G4Material* quartzG4material = G4Material::GetMaterial(fQuartzMaterial);
	if(!quartzG4material) {
		G4cout<<" ----> Material "<<fQuartzMaterial<<" not found, cannot build! "<<G4endl;
		return 0;
	}

	// building the aluminum can
	G4double cutExtra = 10.0*cm;
	G4Tubs* cut = new G4Tubs("cutting volume", fScintillatorInnerRadius, fScintillatorOuterRadius, fScintillatorLength/2.0 + cutExtra, fStartPhi, fEndPhi);
	G4Tubs* can = new G4Tubs("can volume before cut", fScintillatorInnerRadius, fScintillatorOuterRadius + fAlumCanThickness, fScintillatorLength/2.0 + fAlumCanThickness/2.0, fStartPhi, fEndPhi);

	move = G4ThreeVector(0., 0., -0.5*fAlumCanThickness - cutExtra);
	rotate = new G4RotationMatrix;
	G4SubtractionSolid* testcanAlumCasing = new G4SubtractionSolid("testcanAlumCan", can, cut, rotate, move);

	// building the scintillator
	G4Tubs* testcanScintillator = new G4Tubs("testcanScintillator", fScintillatorInnerRadius, fScintillatorOuterRadius, fScintillatorLength/2.0, fStartPhi, fEndPhi);

	// building the quartz window
	G4Tubs* testcanQuartzWindow = new G4Tubs("testcanQuartzWindow", fScintillatorInnerRadius, fQuartzRadius, fQuartzThickness/2.0, fStartPhi, fEndPhi);

	// Set visualization attributes
	G4VisAttributes* canVisAtt = new G4VisAttributes(fGreyColour);
	canVisAtt->SetVisibility(true);
	G4VisAttributes* liquidVisAtt = new G4VisAttributes(fLiquidColour);
	liquidVisAtt->SetVisibility(true);
	G4VisAttributes* quartzVisAtt = new G4VisAttributes(fQuartzColour);
	quartzVisAtt->SetVisibility(true);
	//quartzVisAtt->SetForceWireframe(true);

	// Define rotation and movement objects for aluminum can
	direction 	  = G4ThreeVector(0., 0., 1.);
	move          = G4ThreeVector(0., 0., -1.0*fScintillatorLength/2.0 + fAlumCanThickness/2.0);
	rotate = new G4RotationMatrix;

	//logical volume for aluminum can
	if(fTestcanAlumCasingLog == nullptr) {
		fTestcanAlumCasingLog = new G4LogicalVolume(testcanAlumCasing, canG4material, "testcanAlumCasingLog", 0, 0, 0);
		fTestcanAlumCasingLog->SetVisAttributes(canVisAtt);
	}
	fAssemblyTestcan->AddPlacedVolume(fTestcanAlumCasingLog, move, rotate);

	// Define rotation and movement objects for scintillator
	direction 	  = G4ThreeVector(0., 0., 1.);
	move          = G4ThreeVector(0., 0., -1.0*fScintillatorLength/2.0);
	rotate = new G4RotationMatrix;

	// logical volume for scintillator
	if(fTestcanScintillatorLog == nullptr) {
		fTestcanScintillatorLog = new G4LogicalVolume(testcanScintillator, liquidG4material, "testcanScintillatorLog", 0, 0, 0);
		fTestcanScintillatorLog->SetVisAttributes(liquidVisAtt);
	}
	fAssemblyTestcan->AddPlacedVolume(fTestcanScintillatorLog, move, rotate);

	// Define rotation and movement objects for quartzWindow
	direction 	  = G4ThreeVector(0., 0., 1.);
	move         = G4ThreeVector(0., 0., -1.0*fQuartzThickness/2.0 - 1.0*fScintillatorLength);
	rotate = new G4RotationMatrix;

	// logical volume for quartz window
	if(fTestcanQuartzWindowLog == nullptr) {
		fTestcanQuartzWindowLog = new G4LogicalVolume(testcanQuartzWindow, quartzG4material, "testcanQuartzWindowLog", 0, 0, 0);
		fTestcanQuartzWindowLog->SetVisAttributes(quartzVisAtt);
	}
	fAssemblyTestcan->AddPlacedVolume(fTestcanQuartzWindowLog, move, rotate);

	return 1;
}

