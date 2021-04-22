#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemZDS.hh"

#include "G4SystemOfUnits.hh"

#include <string>

DetectionSystemZDS::DetectionSystemZDS()
{
	// detector dimensions  properties
	fScintillatorWidth         = 1.*mm;
	fPMTWidth         = 1.*mm;
	fDiameter         = 25.*mm;
	fRadialDistance = 3.*mm; //Says in GRIFFIN nim can move from a few mm to 5 cm.  At closest covers roughly 25% of 4pi
	fStartPhi = 0.;
	fDeltaPhi = 2*M_PI;
	fWrapThickness = 0.5 * mm; //Factor of 10 should be  applied later for visualization purposes.  0.05 for simulations, 0.5 for visualization
	fAddWrap = true;
	fWrapMaterial = "Teflon";	
	fZDSMaterial = "BC422";
	//fZDSMaterial = "Dense"; // for Solid Angle test
	fPMTMaterial = "G4_SILICON_DIOXIDE";

	bronze=G4Color(0.8,0.5,0.2);
	silver=G4Color(0.75,0.75,0.75);
	black=G4Color(0.,0.,0.);

	G4cout << "Calling ZDS Constructor" << G4endl;
}
/////
///////
DetectionSystemZDS::~DetectionSystemZDS() {
	// LogicalVolumes
		delete fZDSLog;
		delete fWrapLog;
		delete fPMTLog;

}
////////
/////////
G4int DetectionSystemZDS::Build() {
	fAssemblyZDS = new G4AssemblyVolume(); 
	G4cout << "Calling ZDS Build function" << G4endl;
	BuildZDS();

	return 1;
}
////////
///////
G4int DetectionSystemZDS::PlaceDetector(G4LogicalVolume* expHallLog) {
	G4RotationMatrix * rotate = new G4RotationMatrix;
	G4ThreeVector move = G4ThreeVector(0., 0., 0.);

	//To check overlaps
	fAssemblyZDS->MakeImprint(expHallLog, move, rotate, 0, true);
	G4cout << "Calling ZDS place detector" << G4endl;
	return 1;
}
///////////
/////////
G4int DetectionSystemZDS::BuildZDS() {

	G4cout << "Calling Build ZDS" << G4endl;


	G4ThreeVector move, direction;
	G4RotationMatrix* rotate;

	G4Material* zdsG4material = G4Material::GetMaterial(fZDSMaterial);
	if( !zdsG4material ) {
		G4cout << " ----> Material " << fZDSMaterial << " not found, cannot build! " << G4endl;
		return 0;
	}
	else {
		G4cout << zdsG4material->GetName() << " is the name of the detector material" << G4endl;
	}
	G4Material* pmtG4material = G4Material::GetMaterial(fPMTMaterial);
	if( !pmtG4material ) {
		G4cout << " ----> Material " << fPMTMaterial << " not found, cannot build! " << G4endl;
		return 0;
	}
	else {
		G4cout << pmtG4material->GetName() << " is the name of the pmt material" << G4endl;
	}
	G4Material* wrapG4material = G4Material::GetMaterial(fWrapMaterial);
	if( !wrapG4material ) {
		G4cout << " ----> Material " << fWrapMaterial << " not found, cannot build! " << G4endl;
		return 0;
	}
	else {
		G4cout << wrapG4material->GetName() << " is the name of the wrapping material" << G4endl;
	}
	///// Building the ZDS Geometry /////
	G4Tubs * zds = new G4Tubs("zds", 0., fDiameter/2., fScintillatorWidth/2., fStartPhi, fDeltaPhi);
	G4Tubs * pmt = new G4Tubs("pmt", 0., fDiameter/2., fPMTWidth/2., fStartPhi, fDeltaPhi);

	//For placing volume
	rotate = new G4RotationMatrix;
	
	//Adding wrapper
	G4Tubs * wrapFull = new G4Tubs("wrapFull", 0., fDiameter/2.+fWrapThickness, (fPMTWidth+fScintillatorWidth)/2.+fWrapThickness, fStartPhi, fDeltaPhi);
	G4Tubs * wrapSub = new G4Tubs("wrapSub", 0., fDiameter/2., (fPMTWidth+fScintillatorWidth)/2., fStartPhi, fDeltaPhi);
	G4SubtractionSolid * wrap = new G4SubtractionSolid("wrap", wrapFull, wrapSub, rotate, G4ThreeVector(0,0,0));
	


	//Set visual attributes
	G4VisAttributes * zds_vis = new G4VisAttributes(silver);
	zds_vis->SetVisibility(true);
	G4VisAttributes * pmt_vis = new G4VisAttributes(bronze);
	pmt_vis->SetVisibility(true);
	G4VisAttributes * wrap_vis = new G4VisAttributes(black);
	wrap_vis->SetVisibility(true);

	//Names
	G4String nameLogZDS = "ZDS";
	G4String nameLogPMT = "zdsWindow";
	G4String nameLogWrap = "zdsWrap";
	//Assign Logical Volume for detectors and wrapping affected by beamline
	fZDSLog = new G4LogicalVolume(zds, zdsG4material, nameLogZDS,0,0,0);
	fPMTLog = new G4LogicalVolume(pmt, pmtG4material, nameLogPMT,0,0,0);
	fWrapLog = new G4LogicalVolume(wrap, wrapG4material, nameLogWrap,0,0,0);

	//Give everything colour
	fZDSLog->SetVisAttributes(zds_vis);
	fPMTLog->SetVisAttributes(pmt_vis);
	fWrapLog->SetVisAttributes(wrap_vis);
	move = G4ThreeVector(0., 0., fRadialDistance);
	fAssemblyZDS->AddPlacedVolume(fZDSLog, move, rotate);
	move = G4ThreeVector(0., 0., fRadialDistance+fScintillatorWidth);
	fAssemblyZDS->AddPlacedVolume(fPMTLog, move, rotate);
	move = G4ThreeVector(0., 0., fRadialDistance+(fScintillatorWidth+fPMTWidth)/2. - fWrapThickness);
	if(fAddWrap == true){
	fAssemblyZDS->AddPlacedVolume(fWrapLog, move, rotate);
	}




	return 1;
}

