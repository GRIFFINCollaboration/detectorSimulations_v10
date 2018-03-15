#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

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

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemSpice.hh"
#include "HistoManager.hh" // for spice histogram array when placing detector

DetectionSystemSpice::DetectionSystemSpice() :
	// LogicalVolumes
	fSiInnerGuardRingLog(0),
	fSiOuterGuardRingLog(0)
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
	fDetectorToTargetDistance = 115*mm;
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
	fSiDetCrystalOuterDiameter = 94.*mm;
	fSiDetCrystalInnerDiameter = 16.*mm;
	fSiDetCrystalThickness = 6.15*mm;
	fSiDetRadialSegments = 10.;
	fSiDetPhiSegments = 12.;

	//-----------------------------//
	// parameters for guard ring   //
	//-----------------------------//
	fSiDetGuardRingOuterDiameter = 102*mm;
	fSiDetGuardRingInnerDiameter = 10*mm;

}

DetectionSystemSpice::~DetectionSystemSpice()
{   	
	// LogicalVolumes in ConstructSPICEDetectionSystem
	delete fDetectorMountLog;
	delete fAnnularClampLog;
	for(int i = 0; i < 10; ++i) {
		delete fSiDetSpiceRingLog[i];
	}
	delete fSiInnerGuardRingLog;
	delete fSiOuterGuardRingLog;

	//Physical volumes 
	delete fDetectorMountPhys;
	delete fAnnularClampPhys;
}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemSpice::Build()
{

	fAssembly = new G4AssemblyVolume();

	// Loop through each ring ...
	for(int ringID=0; ringID<10; ringID++) {

		// Build fAssembly volumes
		fAssemblySiRing[ringID] = new G4AssemblyVolume();
		// Build the logical segment of the Silicon Ring
		BuildSiliconWafer(ringID);  

	} // end for(int ringID)

	BuildInnerGuardRing();
	BuildOuterGuardRing();
	BuildDetectorMount();
	BuildAnnularClamps();


	return 1;
} // end Build

//---------------------------------------------------------//
// "place" function called in DetectorMessenger            //
// if detector is added                                    //
//---------------------------------------------------------//
G4int DetectionSystemSpice::PlaceDetector(G4LogicalVolume* exp_hall_log, G4int nRings)
{
	PlaceDetectorMount(exp_hall_log);
	PlaceAnnularClamps(exp_hall_log);   
	PlaceGuardRing(exp_hall_log);//clears the code in det construction


	//To see how the segments fill put 1 detector ring in (via macro) and 1,2,3 segments (changed here) to see how they add.
	//segments add following JanalysisToolkit convention
	G4int NumberSeg = 12; // Segments in Phi for loop below
	G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
	G4ThreeVector pos(0,0,-annularDetectorDistance); 
	for(int ring=0; ring<nRings; ring++)
	{
		for(int Seg=0; Seg<NumberSeg; Seg++)
		{
			G4int TotalNumberSeg = (G4int)fSiDetPhiSegments; // total number of segments = 12
			G4double angle = (360.*deg/TotalNumberSeg)*(-Seg); // Seg = {0, ...,11} //changed to -seg and fills in correctly
			G4RotationMatrix* rotate = new G4RotationMatrix;
			rotate->rotateZ(-angle); // the axis are inverted, this operation will correct for it  [MHD : 03 April 2014]
			//affecting value here (originally 210) may help mapping
			//180 deg gave correct start position for 1st ring, needs to now fill opposite direction
			//could bring in loops to here as to iterate and allow for staggered start position (add 30 deg for each ring etc)
			//as filling in reverse, can establish the histos in reverse then fill in normally 
			// 		  G4double anglerad = (-165*deg-angle)*180/pi;//180 to 165 for true mapping
			fAssemblySiRing[ring]->MakeImprint(exp_hall_log, pos, rotate, Seg);

			// 		  G4cout << "DetectionSystemSpice::Ring number:"<< ring <<  " SegNumber in ring:"<< Seg << " SegmentNumber in detector:" <<segmentID<< " at phi angle"<<anglerad<< G4endl;
			//above for mapping output

		} // end for(int Seg=0; Seg<NumberSeg; Seg++)
	} // end for(int ring = 0; ring<nRings; ring++)

	return 1;
}

G4int DetectionSystemSpice::PlaceGuardRing(G4LogicalVolume* exp_hall_log)
{
	G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
	G4ThreeVector pos(0,0,-annularDetectorDistance); 
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);
	fAssembly->MakeImprint(exp_hall_log, pos, rotate, 0);

	return 1;
}

void DetectionSystemSpice::PlaceDetectorMount(G4LogicalVolume* exp_hall_log)
{
	G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
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
	G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
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
G4int DetectionSystemSpice::BuildSiliconWafer(G4int RingID)  // RingID = { 0, 9 } //Splits each in 2???
{
	// Define the material, return error if not found
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if( !material ) {
		G4cout << " ----> Material " << fWaferMaterial 
			<< " not found, cannot build the detector shell! " << G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	vis_att->SetVisibility(true);

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double z_position		= -(fSiDetCrystalThickness/2.);
	G4ThreeVector move 		= z_position * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);

	// construct solid
	G4Tubs* siDetSpiceRingSec = BuildCrystal(RingID);
	// construct logical volume
	if( !fSiDetSpiceRingLog[RingID] ) {
		G4String ringName = "siDetSpiceRing_";
		G4String ringID = G4UIcommand::ConvertToString(RingID);
		ringName += ringID;
		ringName += "_Log";

		//G4cout << "DetectionSystemSpice : ringName "<< ringName <<" RingID" <<  RingID<< G4endl ; 
		//G4cin.get();

		fSiDetSpiceRingLog[RingID] = new G4LogicalVolume(siDetSpiceRingSec, material, ringName, 0, 0, 0);
		fSiDetSpiceRingLog[RingID]->SetVisAttributes(vis_att);
	}

	fAssemblySiRing[RingID]->AddPlacedVolume(fSiDetSpiceRingLog[RingID], move, rotate);

	return 1;
}

G4int DetectionSystemSpice::BuildInnerGuardRing()
{
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if( !material ) {
		G4cout << " ----> Material " << fWaferMaterial << " not found, cannot build the inner guard ring of Spice! " << G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	vis_att->SetVisibility(true);

	G4Tubs* innerGuardRing = new G4Tubs("innerGuardRing",
			fSiDetGuardRingInnerDiameter/2.,
			fSiDetCrystalInnerDiameter/2.,
			fSiDetCrystalThickness/2.,0,360);

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double z_position		= -(fSiDetCrystalThickness/2.);
	G4ThreeVector move 		= z_position * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);

	//logical volume
	if( fSiInnerGuardRingLog == NULL )
	{
		fSiInnerGuardRingLog = new G4LogicalVolume(innerGuardRing, material, "innerGuardRing", 0,0,0);
		fSiInnerGuardRingLog->SetVisAttributes(vis_att);
	}

	fAssembly->AddPlacedVolume(fSiInnerGuardRingLog, move, rotate);

	return 1;
}

G4int DetectionSystemSpice::BuildOuterGuardRing()
{
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if( !material ) {
		G4cout << " ----> Material " << fWaferMaterial << " not found, cannot build the outer guard ring of Spice! " << G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	vis_att->SetVisibility(true);

	G4Tubs* outerGuardRing = new G4Tubs("outerGuardRing",
			fSiDetCrystalOuterDiameter/2.,
			fSiDetGuardRingOuterDiameter/2.,
			fSiDetCrystalThickness/2.,0,360);

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double z_position		= -(fSiDetCrystalThickness/2.);
	G4ThreeVector move 		= z_position * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);

	//logical volume
	if( fSiOuterGuardRingLog == NULL )
	{
		fSiOuterGuardRingLog = new G4LogicalVolume(outerGuardRing, material, "outerGuardRing", 0,0,0);
		fSiOuterGuardRingLog->SetVisAttributes(vis_att);
	}

	fAssembly->AddPlacedVolume(fSiOuterGuardRingLog, move, rotate);

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


/////////////////////////////////////////////////////////
// Build one segment of Spice, 			       //
// the geometry depends on the distance from the center//
/////////////////////////////////////////////////////////
G4Tubs* DetectionSystemSpice::BuildCrystal(G4int RingID)//called for wafer construction 
{
	// define angle, length, thickness, and inner and outer diameter
	// of silicon detector segment
	G4double tube_element_length = (fSiDetCrystalOuterDiameter - fSiDetCrystalInnerDiameter)/(2*(fSiDetRadialSegments));
	G4double tube_element_angular_width = (360./fSiDetPhiSegments)*deg;
	G4double tube_element_inner_radius = (fSiDetCrystalInnerDiameter)/2.0 + tube_element_length*(RingID);
	G4double tube_element_outer_radius = ((G4double)fSiDetCrystalInnerDiameter)/2.0 + tube_element_length*(RingID+1);
	G4double tube_element_half_thickness = (fSiDetCrystalThickness)/2.0;

	// establish solid
	G4Tubs* crystal_block = new G4Tubs("crystal_block",tube_element_inner_radius,tube_element_outer_radius,tube_element_half_thickness,0,tube_element_angular_width);

	return crystal_block;
}//end ::BuildCrystal
