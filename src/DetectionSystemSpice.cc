#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "ApparatusSpiceTargetChamber.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"

#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemSpice.hh"
#include "HistoManager.hh" // for spice histogram array when placing detector

DetectionSystemSpice::DetectionSystemSpice() 
{
	/////////////////////////////////////////////////////////////////////
	// SPICE Physical Properties
	/////////////////////////////////////////////////////////////////////
	//-----------------------------//
	// Materials     			   //
	//-----------------------------//
	fWaferMaterial          = "Silicon";
	fDetectorMountMaterial = "Aluminum"; // CHECK ?
	fAnnularClampMaterial  = "Peek"; // CHECK ?

	// ----------------------------
	// Dimensions of Detector Mount
	// ----------------------------
	fDetectorMountLength = 148*mm;
	fDetectorMountWidth = 130*mm;
	fDetectorMountThickness = 8*mm;
	fDetectorMountInnerRadius = 52.8*mm;
	fDetectorMountLipRadius = 48.3*mm;
	fDetectorMountLipThickness = 1*mm;
	fDetectorMountAngularOffset = 0*deg;


	fDetectorToTargetDistance = ApparatusSpiceTargetChamber::fTargetChamberDistanceFromTarget +117.6*mm;
	fDetectorThickness = 6*mm;

	// ---------------------------
	// Dimensions of Annular Clamp
	// ---------------------------
	fAnnularClampThickness = 2*mm;
	fAnnularClampLength = 19.2*mm;
	fAnnularClampWidth = 20*mm;
	fAnnularClampPlaneOffset = 47.8*mm; //from beam line to first edge  

	//-----------------------------//
	// parameters for the annular  //
	// planar detector crystal     //
	//-----------------------------//
	fSiDetCrystalOuterDiameter = 102.*mm;
	fSiDetCrystalInnerDiameter = 10.*mm;
	fSiDetCrystalActiveInnerDiameter = 16.0*mm;
	fSiDetCrystalActiveOuterDiameter = 94.0*mm;
	fSiDetRingRadii = 3.9*mm;
	fSiDetCrystalThickness = 6.15*mm;
	fSiDetRadialSegments = 10;
	fSiDetPhiSegments = 12;

	//-----------------------------//
	// parameters for guard ring   //
	//-----------------------------//
	fSiDetGuardRingOuterDiameter = 102*mm;
	fSiDetGuardRingInnerDiameter = 10*mm;

}

DetectionSystemSpice::~DetectionSystemSpice()
{   	
	delete fDetectorMountLog;
	delete fAnnularClampLog;
	for(unsigned int i=0;i<fSiliLog.size();i++)delete fSiliLog[i];
}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemSpice::BuildPlace(G4LogicalVolume* ExpHallLog)
{
	BuildPlaceSiliconWafer(ExpHallLog);  

	BuildDetectorMount();
	PlaceDetectorMount(ExpHallLog);

	BuildAnnularClamps();
	PlaceAnnularClamps(ExpHallLog);

	return 1;
} // end Build



void DetectionSystemSpice::PlaceDetectorMount(G4LogicalVolume* exp_hall_log)
{
	G4double annularDetectorDistance = fDetectorToTargetDistance;
	G4ThreeVector pos(0,0,-annularDetectorDistance); 
	G4double detector_mount_gap = fDetectorMountThickness 
		- fDetectorMountLipThickness - fDetectorThickness;

	G4double z_offset = - fDetectorMountThickness/2. + detector_mount_gap ;
	G4ThreeVector offset(0, 0, z_offset);

	pos = pos + offset ; 
	G4RotationMatrix* rotate = new G4RotationMatrix(fDetectorMountAngularOffset, 0, 0);
	fDetectorMountPhys = new G4PVPlacement(rotate, pos, fDetectorMountLog,
			"detector_mount", exp_hall_log, false, 0);

} // end::PlaceDetectorMount()

void DetectionSystemSpice::PlaceAnnularClamps(G4LogicalVolume* exp_hall_log) {
	G4double annularDetectorDistance = fDetectorToTargetDistance;
	G4ThreeVector pos(0,0,-annularDetectorDistance); 
	G4double z_offset = fAnnularClampThickness/2. ;
	G4double x_offset = (fAnnularClampPlaneOffset
			+ fAnnularClampLength/2.) 
		* cos(fDetectorMountAngularOffset + 45*deg);
	G4double y_offset = (fAnnularClampPlaneOffset
			+ fAnnularClampLength/2.)
		* sin(fDetectorMountAngularOffset + 45*deg);
	G4ThreeVector offset(-x_offset, -y_offset, z_offset);

	pos = pos + offset ; 	 
	G4RotationMatrix* rotate = new G4RotationMatrix(fDetectorMountAngularOffset + 45*deg, 0, 0);
	fAnnularClampPhys = new G4PVPlacement(rotate, pos, fAnnularClampLog,"annular_clamp", exp_hall_log,
			false,0);

} // end::PlaceAnnularClamps()

//---------------------------------------------------------//
// build functions for different parts                     //
// called in main build function                           //
//---------------------------------------------------------//

G4int DetectionSystemSpice::BuildPlaceSiliconWafer(G4LogicalVolume* ExpHallLog) 
{
	// Define the material, return error if not found
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fWaferMaterial<< " not found, cannot build the SPICE detecto! "<<G4endl;
		return 0;
	}
	G4Material* matWorld = G4Material::GetMaterial("Vacuum");
	if(!matWorld) {
		G4cout<<" ----> Material Vacuum not found! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	vis_att->SetVisibility(true);  
	G4VisAttributes* vis_att_hid = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	vis_att_hid->SetVisibility(false);

	// Define rotation and movement objects
	G4RotationMatrix* rotate = new G4RotationMatrix;

	G4double tube_dPhi = 2.* M_PI * rad;

	// Detector crystal mother
	G4Tubs* siDetMother = new G4Tubs("siDetMother",(fSiDetGuardRingInnerDiameter)/2.0,(fSiDetGuardRingOuterDiameter)/2.0 ,(fSiDetCrystalThickness)/2.0,0, tube_dPhi);
	G4LogicalVolume* siDetMotherLog = new G4LogicalVolume(siDetMother, matWorld, "fDetActiveLog", 0,0,0);  
	siDetMotherLog->SetVisAttributes(vis_att_hid);
	fSiliLog.push_back(siDetMotherLog);

	//Guard Rings
	//Solids
	G4Tubs* siDetInnerGuardRing = new G4Tubs("siDetInnerGuardRing",(fSiDetGuardRingInnerDiameter)/2.0,(fSiDetCrystalActiveInnerDiameter)/2.0 ,(fSiDetCrystalThickness)/2.0,0, tube_dPhi);
	G4Tubs* siDetOuterGuardRing = new G4Tubs("siDetOuterGuardRing",(fSiDetCrystalActiveOuterDiameter)/2.0,(fSiDetGuardRingOuterDiameter)/2.0 ,(fSiDetCrystalThickness)/2.0,0, tube_dPhi); 
	//Logical
	G4LogicalVolume* fInnerGuardRingLog = new G4LogicalVolume(siDetInnerGuardRing, material, "fInnerGuardRingLog", 0,0,0);  
	fInnerGuardRingLog->SetVisAttributes(vis_att);
	fSiliLog.push_back(fInnerGuardRingLog);
	G4LogicalVolume* fOuterGuardRingLog = new G4LogicalVolume(siDetOuterGuardRing, material, "fOuterGuardRingLog", 0,0,0);  
	fOuterGuardRingLog->SetVisAttributes(vis_att);
	fSiliLog.push_back(fOuterGuardRingLog);
	//Place physical
	new G4PVPlacement(rotate,G4ThreeVector(0,0,0), fInnerGuardRingLog, "SiLiInnerGuard",  siDetMotherLog, false, 0);
	new G4PVPlacement(rotate,G4ThreeVector(0,0,0), fOuterGuardRingLog, "SiLiOuterGuard",  siDetMotherLog, false, 0);

	G4double divided_tube_dPhi = tube_dPhi/fSiDetPhiSegments;

	//Prepare rotation matrices 
	std::vector<G4RotationMatrix*> rotsec;
	G4double phioffset=180*deg;
	for(int sec=0;sec<fSiDetPhiSegments;sec++){
		G4RotationMatrix* rot = new G4RotationMatrix;
		rot->rotateZ(phioffset - 30*deg*sec);
		rotsec.push_back(rot);
	}

	for(int r=0;r<fSiDetRadialSegments;r++){
		G4Tubs* siDetSegment = new G4Tubs("siDetSegment",(fSiDetCrystalActiveInnerDiameter/2.0)+(fSiDetRingRadii*(r)),(fSiDetCrystalActiveInnerDiameter/2.0)+(fSiDetRingRadii*(r+1)),(fSiDetCrystalThickness)/2.0,0,divided_tube_dPhi);
		G4LogicalVolume* fSiSegmentLog =new G4LogicalVolume(siDetSegment, material, "fSiSegmentLog", 0,0,0);
		fSiSegmentLog->SetVisAttributes(vis_att);
		fSiliLog.push_back(fSiSegmentLog);

		for(int sec=0;sec<fSiDetPhiSegments;sec++){
			std::stringstream temp;
			temp<<"SiSegmentPhys_"<<r<<"_"<<sec;
			new G4PVPlacement(rotsec[sec],G4ThreeVector(0,0,0), fSiSegmentLog, temp.str(),  siDetMotherLog, false, 0);
		}	
	}

	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double z_position		= -(fSiDetCrystalThickness/2.)-fDetectorToTargetDistance;
	G4ThreeVector move 		= z_position * direction;

	new G4PVPlacement(rotate,move, siDetMotherLog, "SiLiDetActive",  ExpHallLog, false, 0);

	//   //// Alternate construction using recommended G4PVReplica class should actually be more efficient but causes particles volumes that arent observed with the simple above method.
	//   G4Tubs* siDetActiveArea = new G4Tubs("siDetActiveArea",(fSiDetCrystalActiveInnerDiameter)/2.0,(fSiDetCrystalActiveOuterDiameter)/2.0 ,(fSiDetCrystalThickness)/2.0,0, tube_dPhi);
	//   G4Tubs* siDetSector = new G4Tubs("siDetSector",(fSiDetCrystalActiveInnerDiameter)/2.0,(fSiDetCrystalActiveOuterDiameter)/2.0 ,(fSiDetCrystalThickness)/2.0,0,divided_tube_dPhi);
	//   G4Tubs* siDetSegment = new G4Tubs("siDetSegment",(fSiDetCrystalActiveInnerDiameter)/2.0,((fSiDetCrystalActiveInnerDiameter)/2.0)+(fSiDetRingRadii),(fSiDetCrystalThickness)/2.0,0,divided_tube_dPhi);
	//   G4LogicalVolume* fDetActiveLog = new G4LogicalVolume(siDetActiveArea, matWorld, "fDetActiveLog", 0,0,0);  
	//   fDetActiveLog->SetVisAttributes(vis_att_hid);
	//   G4LogicalVolume* fSiSectorLog = new G4LogicalVolume(siDetSector, matWorld, "fSiSectorLog", 0,0,0);  
	//   fSiSectorLog->SetVisAttributes(vis_att_hid);
	//   G4LogicalVolume* fSiSegmentLog = new G4LogicalVolume(siDetSegment, material, "fSiSegmentLog", 0,0,0);  
	//   fSiSegmentLog->SetVisAttributes(vis_att);
	//   new G4PVReplica("SiSectorPhys", fSiSectorLog, fDetActiveLog, kPhi, fSiDetPhiSegments, divided_tube_dPhi);
	//   new G4PVReplica("SiSegmentPhys", fSiSegmentLog, fSiSectorLog, kRho, fSiDetRadialSegments, fSiDetRingRadii,fSiDetCrystalActiveInnerDiameter/2.0);
	//   ////A quirk of G4PVReplica, cannot sub-divide a kRho so the order is important
	//   new G4PVPlacement(rotate,move, fDetActiveLog, "SiLiDetActive",  ExpHallLog, false, 0);

	//   //// Other alternate which incorporates guard rings and should be even better, but causes more geometry warnings
	//   G4Tubs* siDetActiveArea = new G4Tubs("siDetActiveArea",(fSiDetGuardRingInnerDiameter)/2.0,(fSiDetGuardRingOuterDiameter)/2.0 ,(fSiDetCrystalThickness)/2.0,0, tube_dPhi);
	//   G4Tubs* siDetSector = new G4Tubs("siDetSector",(fSiDetGuardRingInnerDiameter)/2.0,(fSiDetGuardRingOuterDiameter)/2.0 ,(fSiDetCrystalThickness)/2.0,0,divided_tube_dPhi);
	//   G4Tubs* siDetSegment = new G4Tubs("siDetSegment",(fSiDetCrystalActiveInnerDiameter)/2.0,((fSiDetCrystalActiveInnerDiameter)/2.0)+(fSiDetRingRadii),(fSiDetCrystalThickness)/2.0,0,divided_tube_dPhi);
	//   G4LogicalVolume* fDetActiveLog = new G4LogicalVolume(siDetActiveArea, material, "fDetActiveLog", 0,0,0);  
	//   fDetActiveLog->SetVisAttributes(vis_att);
	//   G4LogicalVolume* fSiSectorLog = new G4LogicalVolume(siDetSector, material, "fSiSectorLog", 0,0,0);  
	//   fSiSectorLog->SetVisAttributes(vis_att_hid);
	//   G4LogicalVolume* fSiSegmentLog = new G4LogicalVolume(siDetSegment, material, "fSiSegmentLog", 0,0,0);  
	//   fSiSegmentLog->SetVisAttributes(vis_att);
	//   new G4PVReplica("SiSectorPhys", fSiSectorLog, fDetActiveLog, kPhi, 12, divided_tube_dPhi);
	//   new G4PVDivision("SiSegmentPhys",fSiSegmentLog,fSiSectorLog,kRho,10,fSiDetRingRadii,(fSiDetCrystalActiveInnerDiameter-fSiDetGuardRingInnerDiameter)/2.);
	//   G4VPhysicalVolume* PlaceActice = new G4PVPlacement(rotate,move, fDetActiveLog, "SiLiDetActive",  ExpHallLog, false, 0);


	return 1;
}

void DetectionSystemSpice::BuildDetectorMount() {

	// ** Visualisation
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
	vis_att->SetVisibility(true);

	// ** Dimensions
	// Box
	G4double box_half_width = fDetectorMountWidth/2.;
	G4double box_half_length = fDetectorMountLength/2.;
	G4double box_half_thickness = fDetectorMountThickness/2.;
	// Inner Radius
	G4double lip_radius = fDetectorMountLipRadius;
	//G4double lip_half_thickness = fDetectorMountLipThickness/2.;
	G4double box_cut_radius = fDetectorMountInnerRadius;
	// Annular Clamp
	G4double clamp_half_thickness = fAnnularClampThickness;
	G4double clamp_half_width = fAnnularClampWidth/2.;
	G4double clamp_half_length = fAnnularClampLength/2.;

	// ** Shapes
	G4Box* mount_box = new G4Box("mount_box", box_half_width, box_half_length, box_half_thickness);
	G4Tubs* inner_radius_cut = new G4Tubs("inner_radius_cut", 0, lip_radius, 2*box_half_thickness, 0, 360*deg);
	G4Tubs* lip_cut = new G4Tubs("lip_cut", 0, box_cut_radius, box_half_thickness, 0, 360*deg);
	G4Box* annular_clamp = new G4Box("annular_clamp", clamp_half_width, clamp_half_length, clamp_half_thickness);

	G4SubtractionSolid* detector_mount_pre = new G4SubtractionSolid("detector_mount_pre", mount_box, inner_radius_cut);
	G4ThreeVector trans(0, 0, fDetectorMountLipThickness);
	G4SubtractionSolid* detector_mount = new G4SubtractionSolid("detector_mount", detector_mount_pre, lip_cut, 0, trans);

	G4double plane_offset = (fAnnularClampPlaneOffset + clamp_half_length) / sqrt(2.);
	G4double z_offset = box_half_thickness;
	G4ThreeVector move(plane_offset, plane_offset, z_offset);
	G4RotationMatrix* rotate = new G4RotationMatrix(45*deg, 0, 0);
	G4SubtractionSolid* detector_mount2 = new G4SubtractionSolid("detector_mount2", detector_mount, annular_clamp, rotate, move);
	move.setX(-plane_offset);
	rotate->rotateZ(90*deg);
	G4SubtractionSolid* detector_mount3 = new G4SubtractionSolid("detector_mount3", detector_mount2, annular_clamp, rotate, move);
	move.setY(-plane_offset);
	rotate->rotateZ(90*deg);
	G4SubtractionSolid* detector_mount4 = new G4SubtractionSolid("detector_mount4", detector_mount3, annular_clamp, rotate, move);
	move.setX(plane_offset);
	rotate->rotateZ(90*deg);
	G4SubtractionSolid* detector_mount5 = new G4SubtractionSolid("detector_mount5", detector_mount4, annular_clamp, rotate, move);

	// ** Logical
	G4Material* DetectorMountMaterial = G4Material::GetMaterial(fDetectorMountMaterial);
	fDetectorMountLog = new G4LogicalVolume(detector_mount5, DetectorMountMaterial, "fDetectorMountLog", 0, 0, 0);
	fDetectorMountLog->SetVisAttributes(vis_att);

	//fAssembly->AddPlacedVolume(fDetectorMountLog, move, rotate); CHECK!!

} // end::BuildDetectorMount()

void DetectionSystemSpice::BuildAnnularClamps() {

	// ** Visualisation
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(PEEK_COL));
	vis_att->SetVisibility(true);

	// ** Dimensions
	G4double clamp_half_length = fAnnularClampLength/2.;
	G4double clamp_half_width = fAnnularClampWidth/2.;
	G4double clamp_half_thickness = fAnnularClampThickness/2.;
	// Distance
	G4double beam_clamp_distance = fAnnularClampPlaneOffset + clamp_half_length;

	// ** Shapes
	G4Box* annular_clamp = new G4Box("annular_clamp", clamp_half_width, clamp_half_length, clamp_half_thickness);

	G4ThreeVector move(2*beam_clamp_distance, 0, 0);
	G4UnionSolid* double_clamps = new G4UnionSolid("double_clamps", annular_clamp, annular_clamp, 0, move);

	G4Box* annular_clamp2 = new G4Box("annular_clamp2", clamp_half_length, clamp_half_width, clamp_half_thickness);
	G4ThreeVector trans(0, 2*beam_clamp_distance, 0);
	G4UnionSolid* double_clamps2 = new G4UnionSolid("double_clamps2", annular_clamp2, annular_clamp2, 0, trans);

	G4ThreeVector trans2(beam_clamp_distance, -beam_clamp_distance, 0);
	G4UnionSolid* four_clamps = new G4UnionSolid("four_clamps", double_clamps, double_clamps2, 0, trans2);

	// ** Logical
	G4Material* AnnularClampMaterial = G4Material::GetMaterial(fAnnularClampMaterial);
	fAnnularClampLog = new G4LogicalVolume(four_clamps, AnnularClampMaterial, "fAnnularClampLog", 0, 0, 0);
	fAnnularClampLog->SetVisAttributes(vis_att);

} // end::BuildAnnularClamps()  

