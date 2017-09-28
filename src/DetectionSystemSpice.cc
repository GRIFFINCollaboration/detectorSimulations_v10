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

using namespace CLHEP;

DetectionSystemSpice::DetectionSystemSpice() :
  // LogicalVolumes
  siInnerGuardRing_log(0),
  siOuterGuardRing_log(0)
{
    /////////////////////////////////////////////////////////////////////
    // SPICE Physical Properties
    /////////////////////////////////////////////////////////////////////
	//-----------------------------//
    // Materials     			   //
    //-----------------------------//
    this->wafer_material          = "Silicon";
  	this->detector_mount_material = "Aluminum"; // CHECK ?
  	this->annular_clamp_material  = "Peek"; // CHECK ?
  
	// ----------------------------
	// Dimensions of Detector Mount
	// ----------------------------
	this->detector_mount_length = 148*mm;
	this->detector_mount_width = 130*mm;
	this->detector_mount_thickness = 8*mm;
	this->detector_mount_inner_radius = 52.8*mm;
	this->detector_mount_lip_radius = 48.3*mm;
	this->detector_mount_lip_thickness = 1*mm;
	this->detector_mount_angular_offset = 0*deg;
	this->detector_to_target_distance = 115*mm;
	this->detector_thickness = 6*mm;

	// ---------------------------
	// Dimensions of Annular Clamp
	// ---------------------------
	this->annular_clamp_thickness = 2*mm;
	this->annular_clamp_length = 19.2*mm;
	this->annular_clamp_width = 20*mm;
	this->annular_clamp_plane_offset = 47.8*mm; //from beam line to first edge  

    //-----------------------------//
    // parameters for the annular  //
    // planar detector crystal     //
    //-----------------------------//
    this->siDetCrystalOuterDiameter = 94.*mm;
    this->siDetCrystalInnerDiameter = 16.*mm;
    this->siDetCrystalThickness = 6.15*mm;
    this->siDetRadialSegments = 10.;
    this->siDetPhiSegments = 12.;

    //-----------------------------//
    // parameters for guard ring   //
    //-----------------------------//
    this->siDetGuardRingOuterDiameter = 102*mm;
    this->siDetGuardRingInnerDiameter = 10*mm;
    
}

DetectionSystemSpice::~DetectionSystemSpice()
{   
	// LogicalVolumes in ConstructSPICEDetectionSystem
	delete detector_mount_log;
	delete annular_clamp_log;
	delete siDetSpiceRing_log;//[] siDetSpiceRing_log; 12/7
	delete siInnerGuardRing_log;
	delete siOuterGuardRing_log;
	
	//Physical volumes 
	delete detector_mount_phys;
	delete annular_clamp_phys;
}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemSpice::Build()
{
	
    this->assembly = new G4AssemblyVolume();
   
	// Loop through each ring ...
	for(int ringID=0; ringID<10; ringID++) {
				
		// Build assembly volumes
	  	this->assemblySiRing[ringID] = new G4AssemblyVolume();
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
  G4int NumberSeg = 12; // Segments in Phi for loop below - originally 12
  G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
  G4ThreeVector pos(0,0,-annularDetectorDistance); 
    for(int ring=0; ring<nRings; ring++)
    {
      for(int Seg=0; Seg<NumberSeg; Seg++)
		{
		  G4int TotalNumberSeg = (G4int)this->siDetPhiSegments; // total number of segments = 12
		  G4double angle = (360.*deg/TotalNumberSeg)*(-Seg); // Seg = {0, ...,11} //changed to -seg and fills in correctly
		  G4RotationMatrix* rotate = new G4RotationMatrix;
		  rotate->rotateZ(-180*deg-angle); // the axis are inverted, this operation will correct for it  [MHD : 03 April 2014]
		      //affecting value here (originally 210) may help mapping
		      //180 deg gave correct start position for 1st ring, needs to now fill opposite direction
		      //could bring in loops to here as to iterate and allow for staggered start position (add 30 deg for each ring etc)
		      //as filling in reverse, can establish the histos in reverse then fill in normally 
// 		  G4double anglerad = (-165*deg-angle)*180/pi;//180 to 165 for true mapping
		  assemblySiRing[ring]->MakeImprint(exp_hall_log, pos, rotate, Seg);
		  
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
  assembly->MakeImprint(exp_hall_log, pos, rotate, 0);

  return 1;
}

void DetectionSystemSpice::PlaceDetectorMount(G4LogicalVolume* exp_hall_log)
{
  G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
  G4ThreeVector pos(0,0,-annularDetectorDistance); 
  G4double detector_mount_gap = this->detector_mount_thickness 
    - this->detector_mount_lip_thickness - this->detector_thickness;

  G4double z_offset = - this->detector_mount_thickness/2. + detector_mount_gap ;
  G4ThreeVector offset(0, 0, z_offset);

  pos = pos + offset ; 
  G4RotationMatrix* rotate = new G4RotationMatrix(this->detector_mount_angular_offset, 0, 0);
  detector_mount_phys = new G4PVPlacement(rotate, pos, detector_mount_log,
	"detector_mount", exp_hall_log, false, 0);

} // end::PlaceDetectorMount()

void DetectionSystemSpice::PlaceAnnularClamps(G4LogicalVolume* exp_hall_log) {
    G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
  G4ThreeVector pos(0,0,-annularDetectorDistance); 
  G4double z_offset = this->annular_clamp_thickness/2. ;
  G4double x_offset = (this->annular_clamp_plane_offset
		       + this->annular_clamp_length/2.) 
    * cos(this->detector_mount_angular_offset + 45*deg);
  G4double y_offset = (this->annular_clamp_plane_offset
		       + this->annular_clamp_length/2.)
    * sin(this->detector_mount_angular_offset + 45*deg);
  G4ThreeVector offset(-x_offset, -y_offset, z_offset);
 
  pos = pos + offset ; 	 
  G4RotationMatrix* rotate = new G4RotationMatrix(this->detector_mount_angular_offset + 45*deg, 0, 0);
  annular_clamp_phys = new G4PVPlacement(rotate, pos, annular_clamp_log,"annular_clamp", exp_hall_log,
					 false,0);
  
} // end::PlaceAnnularClamps()

//---------------------------------------------------------//
// build functions for different parts                     //
// called in main build function                           //
//---------------------------------------------------------//
G4int DetectionSystemSpice::BuildSiliconWafer(G4int RingID)  // RingID = { 0, 9 } //Splits each in 2???
{
	// Define the material, return error if not found
  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
  	G4cout << " ----> Material " << this->wafer_material 
  				 << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }

  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  vis_att->SetVisibility(true);

  // Define rotation and movement objects
  G4ThreeVector direction 	= G4ThreeVector(0,0,1);
  G4double z_position		= -(this->siDetCrystalThickness/2.);
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  // construct solid
  G4Tubs* siDetSpiceRingSec = BuildCrystal(RingID);
  // construct logical volume
  if( !siDetSpiceRing_log[RingID] ) {
		G4String ringName = "siDetSpiceRing_";
		G4String ringID = G4UIcommand::ConvertToString(RingID);
		ringName += ringID;
		ringName += "_Log";
		
		//G4cout << "DetectionSystemSpice : ringName "<< ringName <<" RingID" <<  RingID<< G4endl ; 
		//G4cin.get();
		
		siDetSpiceRing_log[RingID] = new G4LogicalVolume(siDetSpiceRingSec, material, ringName, 0, 0, 0);
		siDetSpiceRing_log[RingID]->SetVisAttributes(vis_att);
	}
	
	this->assemblySiRing[RingID]->AddPlacedVolume(siDetSpiceRing_log[RingID], move, rotate);

  return 1;
}

G4int DetectionSystemSpice::BuildInnerGuardRing()
{
  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
    G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the inner guard ring of Spice! " << G4endl;
    return 0;
  }

  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  vis_att->SetVisibility(true);

  G4Tubs* innerGuardRing = new G4Tubs("innerGuardRing",
					 this->siDetGuardRingInnerDiameter/2.,
					 this->siDetCrystalInnerDiameter/2.,
					 this->siDetCrystalThickness/2.,0,360);

  // Define rotation and movement objects
  G4ThreeVector direction 	= G4ThreeVector(0,0,1);
  G4double z_position		= -(this->siDetCrystalThickness/2.);
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  //logical volume
  if( siInnerGuardRing_log == NULL )
  {
    siInnerGuardRing_log = new G4LogicalVolume(innerGuardRing, material, "innerGuardRing", 0,0,0);
    siInnerGuardRing_log->SetVisAttributes(vis_att);
  }

  this->assembly->AddPlacedVolume(siInnerGuardRing_log, move, rotate);

  return 1;
}

G4int DetectionSystemSpice::BuildOuterGuardRing()
{
  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
    G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the outer guard ring of Spice! " << G4endl;
    return 0;
  }

  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  vis_att->SetVisibility(true);

  G4Tubs* outerGuardRing = new G4Tubs("outerGuardRing",
					 this->siDetCrystalOuterDiameter/2.,
					 this->siDetGuardRingOuterDiameter/2.,
					 this->siDetCrystalThickness/2.,0,360);

  // Define rotation and movement objects
  G4ThreeVector direction 	= G4ThreeVector(0,0,1);
  G4double z_position		= -(this->siDetCrystalThickness/2.);
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  //logical volume
  if( siOuterGuardRing_log == NULL )
  {
    siOuterGuardRing_log = new G4LogicalVolume(outerGuardRing, material, "outerGuardRing", 0,0,0);
    siOuterGuardRing_log->SetVisAttributes(vis_att);
  }

  this->assembly->AddPlacedVolume(siOuterGuardRing_log, move, rotate);

  return 1;
}


void DetectionSystemSpice::BuildDetectorMount() {

  // ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
  vis_att->SetVisibility(true);
  
  // ** Dimensions
  // Box
  G4double box_half_width = this->detector_mount_width/2.;
  G4double box_half_length = this->detector_mount_length/2.;
  G4double box_half_thickness = this->detector_mount_thickness/2.;
  // Inner Radius
  G4double lip_radius = this->detector_mount_lip_radius;
  //G4double lip_half_thickness = this->detector_mount_lip_thickness/2.;
  G4double box_cut_radius = this->detector_mount_inner_radius;
  // Annular Clamp
  G4double clamp_half_thickness = this->annular_clamp_thickness;
  G4double clamp_half_width = this->annular_clamp_width/2.;
  G4double clamp_half_length = this->annular_clamp_length/2.;
  
  // ** Shapes
  G4Box* mount_box = new G4Box("mount_box", box_half_width, box_half_length, box_half_thickness);
  G4Tubs* inner_radius_cut = new G4Tubs("inner_radius_cut", 0, lip_radius, 2*box_half_thickness, 0, 360*deg);
  G4Tubs* lip_cut = new G4Tubs("lip_cut", 0, box_cut_radius, box_half_thickness, 0, 360*deg);
  G4Box* annular_clamp = new G4Box("annular_clamp", clamp_half_width, clamp_half_length, clamp_half_thickness);
  
  G4SubtractionSolid* detector_mount_pre = new G4SubtractionSolid("detector_mount_pre", mount_box, inner_radius_cut);
  G4ThreeVector trans(0, 0, this->detector_mount_lip_thickness);
  G4SubtractionSolid* detector_mount = new G4SubtractionSolid("detector_mount", detector_mount_pre, lip_cut, 0, trans);
  
  G4double plane_offset = (this->annular_clamp_plane_offset + clamp_half_length) / sqrt(2.);
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
  G4Material* detector_mount_material = G4Material::GetMaterial(this->detector_mount_material);
  detector_mount_log = new G4LogicalVolume(detector_mount5, detector_mount_material, "detector_mount_log", 0, 0, 0);
  detector_mount_log->SetVisAttributes(vis_att);
  
    //this->assembly->AddPlacedVolume(detector_mount_log, move, rotate); CHECK!!
    
} // end::BuildDetectorMount()

void DetectionSystemSpice::BuildAnnularClamps() {

  // ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(PEEK_COL));
  vis_att->SetVisibility(true);
  
  // ** Dimensions
  G4double clamp_half_length = this->annular_clamp_length/2.;
  G4double clamp_half_width = this->annular_clamp_width/2.;
  G4double clamp_half_thickness = this->annular_clamp_thickness/2.;
  // Distance
  G4double beam_clamp_distance = this->annular_clamp_plane_offset + clamp_half_length;
  
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
  G4Material* annular_clamp_material = G4Material::GetMaterial(this->annular_clamp_material);
  annular_clamp_log = new G4LogicalVolume(four_clamps, annular_clamp_material, "annular_clamp_log", 0, 0, 0);
  annular_clamp_log->SetVisAttributes(vis_att);
  
} // end::BuildAnnularClamps()  


/////////////////////////////////////////////////////////
// Build one segment of Spice, 			       //
// the geometry depends on the distance from the center//
/////////////////////////////////////////////////////////
G4Tubs* DetectionSystemSpice::BuildCrystal(G4int RingID)//called for wafer construction, splits here?? 
{
  // define angle, length, thickness, and inner and outer diameter
  // of silicon detector segment
  G4double tube_element_length = (this->siDetCrystalOuterDiameter - this->siDetCrystalInnerDiameter)/(2*(this->siDetRadialSegments));
  G4double tube_element_angular_width = (360./this->siDetPhiSegments)*deg;
  G4double tube_element_inner_radius = (this->siDetCrystalInnerDiameter)/2.0 + tube_element_length*(RingID);
  G4double tube_element_outer_radius = ((G4double)this->siDetCrystalInnerDiameter)/2.0 + tube_element_length*(RingID+1);
  G4double tube_element_half_thickness = (this->siDetCrystalThickness)/2.0;

  // establish solid
  G4Tubs* crystal_block = new G4Tubs("crystal_block",tube_element_inner_radius,tube_element_outer_radius,tube_element_half_thickness,0,tube_element_angular_width);

  return crystal_block;
}//end ::BuildCrystal
