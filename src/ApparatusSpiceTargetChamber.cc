#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4VSolid.hh"
#include "G4Trap.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "ApparatusSpiceTargetChamber.hh"
#include <math.h>

using namespace CLHEP;

////////////////////////////////////////////////
// building and placing the different parts   //
// of the SPICE target chamber                //
// including chamber itself, magnets and      //
// mounts for magnets and detectors           //
// not including detectors                    //
////////////////////////////////////////////////


//////////////////////////////////////////////////////
// define physical properties of ApparatusSpiceTargetChamber //
//////////////////////////////////////////////////////

ApparatusSpiceTargetChamber::ApparatusSpiceTargetChamber(G4String MedLo)//parameter chooses which lens is in place.
{

  this->NUMBER_OF_MAGNETS = 4;
  this->fNumberOfFrames = 3;

  // Materials
  this->magnet_material = "NdFeB"; 
  this->magnet_clamp_material = "Aluminum";
  this->target_chamber_material = "Delrin"; 
  this->photon_shield_layer_one_material = "WTa"; 
  this->photon_shield_layer_two_material = "Tin"; 
  this->photon_shield_layer_three_material = "Copper";
  this->downstream_cone_material = "Aluminum"; 
  this->ps_clamp_material = "Aluminum";
  this->target_wheel_material = "Peek"; 
  this->target_wheel_gear_material_a = "Delrin";
  this->target_wheel_gear_material_b = "Aluminum";
  this->gear_plate_material = "Aluminum";
  this->target_mount_plate_material = "Peek"; 
  this->bias_plate_material = "Aluminum"; 
  this->gear_stick_material = "Peek"; 
  this->electro_box_material = "Aluminum"; 
  this->small_bolt_material = "Brass"; 
  this->large_bolt_material = "Titanium";
  this->shield_cover_material = "Kapton";
  this->magnet_cover_material = "Peek";
  this->cold_finger_material = "Copper";
  this->s3_cable_case_material = "Delrin";
  this->beam_pipe_material = "Copper";
   
  //-----------------------------
  // Dimensions of Target Chamber
  //-----------------------------
  this->target_chamber_cylinder_inner_radius = 94.6*mm;
  this->target_chamber_cylinder_outer_radius = 102*mm;
  this->target_chamber_cylinder_length = 89*mm;
  this->target_chamber_distance_from_target = 1*mm;
  this->target_chamber_cylinder_lip_inner_radius = 48*mm;
  this->target_chamber_cylinder_lip_thickness = 6*mm;
  this->target_chamber_cylinder_det_outer_radius = 128*mm;
  this->target_chamber_cylinder_det_thickness = 10*mm;
  this->target_chamber_cylinder_hole_length = 32*mm;
  this->target_chamber_cylinder_hole_breadth = 17.7*mm;
  this->target_chamber_cylinder_hole_distance_y = 31.39*mm;
  this->target_chamber_cylinder_hole_distance_x = 74.64*mm;
  
  // Sphere
  this->target_chamber_sphere_inner_radius = 97*mm;
  this->target_chamber_sphere_outer_radius = 102*mm;
  this->target_chamber_sphere_centre = 6*mm;
  this->target_chamber_sphere_cutting_outer_radius = 21*mm;
  this->target_chamber_sphere_cut_off_point = 92.8*mm;
  // Front Ring
  this->front_ring_inner_radius = 15*mm;
  this->front_ring_outer_radius = 34*mm;
  this->front_ring_length = 11.4*mm;
  this->front_ring_cylinder_outer_radius = 19*mm;
  this->front_ring_cylinder_length = 52.607*mm;
  
  this->cone_front_inner_rad = 48.801*mm;
  this->cone_front_outer_rad = 52.801*mm;
  this->cone_length = 92.868*mm;
  
  // --------------------------
  // Dimensions of Target Wheel
  // --------------------------
  this->target_wheel_radius = 40*mm;
  this->target_wheel_thickness = 1.5*mm; // 3.0 for old
  this->target_wheel_offset = 15*mm;
  this->target_radius = 9*mm; //6*mm; // 4mm for old
  this->target_offset = 15*mm;
  this->collimator_radius = 9*mm;
  
  // Bias Plate
  this->bias_plate_outer_radius = 61*mm;
  this->bias_plate_thickness = 0.8*mm;
  this->bias_plate_offset_z = -0.2*mm;
  this->bias_plate_cut = 100*mm;
  
  // Gears
  this->first_gear_radius = 12.7*mm;
  this->first_gear_plane_offset = 80.25*mm;
  this->second_gear_radius = 15.875*mm;
  this->second_gear_plane_offset = 51.675*mm;
  this->third_gear_radius = 50.8*mm;
  this->third_gear_plane_offset = 15*mm;
  this->all_gear_thickness = 3.175*mm;
  this->all_gear_z_offset = 6.459*mm;
  
  // Gear Plates
  this->gear_plate_one_radius = 5*mm;
  this->gear_plate_two_radius = 6*mm;
  this->gear_plate_thickness = 3.866*mm;
   
  // Mount Plate
  this->target_mount_plate_radius = 87*mm;
  this->target_mount_plate_thickness = 5*mm;
  this->target_mount_plate_z_offset = 13.5*mm;
  
  // Gear Stick
  this->gear_stick_length = 250*mm;
  this->gear_stick_radius = 3.18*mm;
  this->gear_stick_z_offset = 5.072*mm;
  
  // Target Frame
  this->fDimTargetFrameOutZ = 0.5*mm;
  this->fDimTargetFrameOutY = 12.*mm;
  this->fDimTargetFrameCutY = 8.*mm;

  //----------------------------
  // Dimensions of Photon Shield
  //----------------------------
  this->photon_shield_front_radius = 11.65*mm; //12.64*mm;
  this->photon_shield_back_radius = 23.19*mm; //24.108*mm;
  this->photon_shield_length = 30.*mm;
  this->photon_shield_inner_radius = 5.*mm;
  this->photon_shield_back_face_pos = -58.*mm;
 
  this->photon_shield_layer_one_thickness = 25*mm;
  this->photon_shield_layer_two_thickness = 4*mm;
  this->photon_shield_layer_three_thickness = 1*mm;
  
  // ----------------------------------
  // Dimensions of Photon Shield Clamps
  // ----------------------------------
  this->ps_target_clamp_radius = 11*mm;
  this->ps_target_clamp_thickness = 2*mm;
  this->ps_target_clamp_extension = 8*mm;
  this->ps_det_clamp_radius = 21.5*mm;
  this->ps_det_clamp_thickness = 3*mm;
  
  // ---------------------------------
  // Dimensions of Photon Shield Bolts
  // ---------------------------------
  this->ps_target_bolt_length = 9.5*mm;
  this->ps_target_bolt_radius = 1.5*mm;
  this->ps_target_bolt_plane_offset = 9*mm;
  this->ps_target_bolt_head_radius = 3.2*mm;
  this->ps_target_bolt_head_length = 1.2*mm;
  this->ps_det_bolt_length = 18.05*mm;
  this->ps_det_bolt_radius = 2.413*mm;
  this->ps_det_bolt_plane_offset = 11*mm;
  this->ps_det_bolt_head_radius = 4.854*mm;
  this->ps_det_bolt_head_length = 3.2*mm;

  //---------------------
  // Dimensions of Magnet
  //---------------------
  this->plate_one_thickness = 3*mm;
  this->plate_one_length = 75*mm; 
  this->plate_one_height = 50*mm; // z-component
  this->plate_one_lower_height = 20*mm;
  if (MedLo = "Med")  {//this works - seen output on terminal via simpel G4cout command
	this->plate_two_thickness = 3.4*mm; // LEL is 5.0, MEL is 3.4
    	this->plate_two_length = 55.0*mm; // LEL is 30.0, MEL is 55.0
  }
  else if (MedLo = "Lo") {
	this->plate_two_thickness = 5.0*mm; // LEL is 5.0, MEL is 3.4
    	this->plate_two_length = 30.0*mm; // LEL is 30.0, MEL is 55.0
  }
  this->plate_two_height = 50*mm;
  
  this->cutting_box_angle = 60*deg;
  this->distance_from_target = 4*mm; // in z-direction
 
  this->plate_one_edge_x = (92*mm - this->plate_one_length); // 92mm is radius of chamber
  
  // ------------------------------------
  // Dimensions of Magnet Clamp (Chamber)
  // ------------------------------------
  this->mclamp_chamber_thickness = this->plate_one_thickness 
    + (2*this->plate_two_thickness) + (2*3.25*mm);
  this->mclamp_chamber_length = 6.5*mm;
  this->mclamp_chamber_height = 83*mm;
  
  // ------------------------------------------
  // Dimensions of Magnet Clamp (Photon Shield)
  // ------------------------------------------
  this->psclamp_psbuffer = 1*mm;
  this->psclamp_chamber_buffer = 5*mm;
  this->psclamp_thickness_buffer = 1.25*mm; // on each side
  this->psclamp_length = this->plate_one_length 
    - this->plate_two_length + this->psclamp_psbuffer 
    + this->psclamp_chamber_buffer;
  
  // -------------------------
  // Dimensions of ElectroBox
  // -------------------------
  this->electrobox_outer_box_length = 269*mm;
  this->electrobox_inner_box_length = 243.5*mm;
  this->electrobox_mid_length = 660*mm;
  this->electrobox_lip_length = 21.5*mm;
  this->electrobox_lip_inner_radius = 102*mm;
  this->electrobox_z_offset = -90*mm;
 
  // ------------------------------
  // Dimensions of Plastic Coatings
  // ------------------------------
  this->magnet_coating_thickness = 1.25*mm;
  this->shield_coating_thickness = 0.2*mm;
    
  // -------------------------------------
  // individual offsets for visualisation
  // -------------------------------------
  this->frontDomeOffset = 0*mm;
  this->targetWheelOffset = 0*mm;
  this->middleRingOffset = 0*mm;
  this->backDetAndAlBoxOffset = 0*mm;

  // -------------------------
  // Dimensions of ColdFinger
  // -------------------------
  this->coldfinger_thickness = 5*mm ; // CHECK
  this->coldfinger_length = 150*mm ; // CHECK
  this->coldfinger_width = 200*mm  ; // CHECK
  this->coldfinger_z_offset = -122*mm ; // CHECK
  this->coldfinger_hole_radius = 10*mm ; // CHECK
  
  // -----------------------
  // Dimensions of Beam Pipe
  // -----------------------
  this->pipe_inner_radius = 5.0*mm; // 5
  this->pipe_outer_radius = 16.0*mm;
  this->pipe_z_length = 30.0*mm;

  
} // end ApparatusSpiceTargetChamber

//////////////////////////////////////////////////
// delete ApparatusSpiceTargetChamber
/////////////////////////////////////////////////

ApparatusSpiceTargetChamber::~ApparatusSpiceTargetChamber()
{

  // logical volumes in ApparatusSpiceTargetChamber
  delete target_chamber_front_ring_log;
  delete target_chamber_front_cone_log;
  delete target_chamber_sphere_log;
  delete target_chamber_cylinder_down_log;
  delete target_wheel_log;
  delete fFrameLogical;
  delete fCollimLogical;
  delete first_gear_log;
  delete second_gear_log;
  delete third_gear_log;
  delete gear_plate_one_log;
  delete gear_plate_two_log;
  delete gear_stick_log;
  delete target_mount_plate_log;
  delete bias_plate_log;
  delete photon_shield_layer_one_log;
  delete photon_shield_layer_two_log;
  delete photon_shield_layer_three_log;
  delete ps_target_clamp_log;
  delete ps_detector_clamp_log;
  delete target_bolt_log;
  delete detector_bolt_log;
  delete magnet_log;
  delete mclamp_chamber_log;
  delete mclamp_shield_log;
  delete electro_box_log;
  delete shield_cover_log;
  delete magnet_cover_log;
  delete cold_finger_log;
  delete beam_pipe_log;
 
  // physical volumes in ApparatusSpiceTargetChamber
  delete target_chamber_front_ring_phys;
  delete target_chamber_front_cone_phys;
  delete target_chamber_sphere_phys;
  delete target_chamber_cylinder_down_phys;
  delete target_wheel_phys;
  delete fFramePhysical;
  delete fCollimPhysical;
  delete first_gear_phys;
  delete second_gear_phys;
  delete third_gear_phys;
  delete gear_plate_one_phys;
  delete gear_plate_two_phys;
  delete gear_stick_phys;
  delete target_mount_plate_phys;
  delete bias_plate_phys;
  delete photon_shield_layer_one_phys;
  delete photon_shield_layer_two_phys;
  delete photon_shield_layer_three_phys;
  delete ps_target_clamp_phys;
  delete ps_detector_clamp_phys;
  delete target_bolt_phys;
  delete detector_bolt_phys;
  delete magnet_phys;
  delete mclamp_chamber_phys;
  delete mclamp_shield_phys;
  delete electro_box_phys;
  delete shield_cover_phys;
  delete magnet_cover_phys;
  delete cold_finger_phys;
  delete beam_pipe_phys;

} // end ~ApparatusSpiceTargetChamber

///////////////////////////////////////
// build and place components        //
///////////////////////////////////////

void ApparatusSpiceTargetChamber::Build(G4LogicalVolume* exp_hall_log)
{
  this->expHallLog = exp_hall_log; // why? should be defined at start?
				    //segmentation fault otherwise - mike
				    
  BuildTargetChamberFrontRing();
  BuildTargetChamberSphere();
  BuildTargetChamberCylinderDownstream();
  BuildTargetWheel();
  BuildTargetWheelExtension();
  BuildTargetWheelRods();
  BuildTargetFrame();
  BuildCollimator();
  BuildTargetWheelGears();
  BuildTargetWheelGearPlates();
  BuildGearStick();
  BuildTargetMountPlate();
  BuildBiasPlate();
  BuildPhotonShield();
  BuildPhotonShieldClamps();
  BuildPhotonShieldClampBolts();
  BuildCollectorMagnet();
  BuildMagnetClampChamber();
  BuildMagnetClampPhotonShield();
  BuildElectroBox();
  BuildShieldCovering();
  BuildMagnetCovering();
  BuildColdFinger();  
  BuildS3CableHolder();
  BuildBeamPipe();

  PlaceTargetChamberFrontRing();
  PlaceTargetChamberSphere();

  PlaceTargetChamberCylinderDownstream(); 

  PlaceTargetWheel(); 
  PlaceTargetWheelExtension(); 
  PlaceCollimator();  
  
  PlaceTargetWheelGears(); 
  PlaceTargetWheelGearPlates(); 
  PlaceGearStick(); 
  PlaceTargetMountPlate(); 
  PlaceBiasPlate(); //commented out originally


  PlacePhotonShield(); 
  PlaceShieldCovering();
  PlacePhotonShieldClamps();
  PlaceElectroBox();
  PlaceColdFinger(); 

  PlaceS3CableHolder();

  PlaceBeamPipe(); //commented out originally (blocked beam) as inner radius was 0

  G4double fFrameDegrees[2] = {0*deg, -110*deg};
  for(G4int i=0; i<2; i++)
    PlaceTargetFrame(fFrameDegrees[i]);
    
  G4double fRodDegrees[7] = {40*deg, -40*deg, -70*deg, -150*deg, 91.2*deg, 125.*deg, 158.8*deg};
  G4double fRodRadius[7] = {10.624, 10.624, 10.624, 10.624, 18.166, 8.4, 18.166 };
  for(G4int i=0; i<7; i++)
    PlaceTargetWheelRods(fRodDegrees[i], fRodRadius[i]);
  
  for(G4int copyID=0; copyID<this->NUMBER_OF_MAGNETS; copyID++)    {
      PlaceCollectorMagnet(copyID);
      PlaceMagnetClampChamber(copyID);
      PlaceMagnetClampPhotonShield(copyID);
      PlacePhotonShieldClampBolts(copyID);
      PlaceMagnetCovering(copyID);
      
    }
  
} // end Build

////////////////////////////////////////////////////
// methods used to build and place the components:
////////////////////////////////////////////////////

void ApparatusSpiceTargetChamber::BuildTargetWheelRods() {

  // Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(AL_COL));
	sVisAtt->SetVisibility(true);
	
	// Dimensions
	G4double sDimRadOuter = 1.0*mm;
	G4double sDimHalfThick = 3.75*mm;
		
	// Shapes
  G4Tubs *sRodTubs = new G4Tubs("sRodTubs", 0, sDimRadOuter, sDimHalfThick, 0, 360*deg);
	
  // Logical
	G4Material* sMaterial = G4Material::GetMaterial(this->target_wheel_material);
	fRodLogical = new G4LogicalVolume(sRodTubs, sMaterial, "target_rods_log", 0, 0, 0);
	fRodLogical->SetVisAttributes(sVisAtt);


}

void ApparatusSpiceTargetChamber::BuildBeamPipe() {

  // Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(AL_COL));
	sVisAtt->SetVisibility(true);
	
	// Dimensions
	G4double sDimRadOuter = this->pipe_outer_radius;
	G4double sDimRadInner = this->pipe_inner_radius;
	G4double sDimHalfThick = this->pipe_z_length / 2.;

	// Shapes
  G4Tubs *sPipe = new G4Tubs("sPipe", sDimRadInner, sDimRadOuter, sDimHalfThick, 0, 360*deg);
	
  // Logical
	G4Material* sMaterial = G4Material::GetMaterial(this->beam_pipe_material);
	beam_pipe_log = new G4LogicalVolume(sPipe, sMaterial, "beam_pipe_log", 0, 0, 0);
	beam_pipe_log->SetVisAttributes(sVisAtt);


}

void ApparatusSpiceTargetChamber::BuildTargetFrame() {

  // Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(0.0, 0.9, 0.1));
	sVisAtt->SetVisibility(true);	

  // Dimensions
	G4double sDimOuterRadius = 7.00*mm;
	G4double sDimInnerRadius = 4.00*mm;
	G4double sDimHalfThick = 0.25*mm;
		
	// Shapes
  G4Tubs* sTubs = new G4Tubs("sTubs", sDimInnerRadius, sDimOuterRadius, sDimHalfThick, 0, 360.*deg);
  G4Tubs* sLimit = new G4Tubs("sLimit", 20.648, 25.0, 2*sDimHalfThick, 0, 360.*deg);
   
  G4ThreeVector sTrans(15.*mm, 0., 0.);
  G4SubtractionSolid *sSub = new G4SubtractionSolid("sSub", sTubs, sLimit, 0, sTrans);
  
  G4Box* sBox = new G4Box("sBox", 2.0*mm, 1.0*mm, sDimHalfThick);
  sTrans.setY(9.0*sin(-45*deg));
  sTrans.setX(9.0*cos(-45*deg));
  
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  sRotate->rotateZ(45.*deg);
  G4UnionSolid *sUnion = new G4UnionSolid("sUnion", sSub, sBox, sRotate, sTrans);
  
  sTrans.setY(9.0*sin(45*deg));
  sTrans.setX(9.0*cos(45*deg));
  
  sRotate->rotateZ(-90.*deg);
  G4UnionSolid *sUnion2 = new G4UnionSolid("sUnion", sUnion, sBox, sRotate, sTrans);
  
  	
	// Logical
	G4Material *sMaterial = G4Material::GetMaterial("Aluminum");
	fFrameLogical = new G4LogicalVolume(sUnion2, sMaterial, "target_frame", 0, 0, 0);
	fFrameLogical->SetVisAttributes(sVisAtt);

}

void ApparatusSpiceTargetChamber::BuildCollimator() {

  // Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	sVisAtt->SetVisibility(true);	
	
	// Dimensions
	
	G4Box* sRectangle = new G4Box("sRectangle", 0.57*mm, 12.*mm, 0.5*mm);
	G4Trd* sTrd = new G4Trd("sTrd", 0.5*mm, 0.5*mm, 5.45*mm, 12.*mm, 3.7*mm);

  G4ThreeVector sTrans(-4.27, 0., 0.);	
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  sRotate->rotateY(-90.*deg);
	G4UnionSolid *sUnion = new G4UnionSolid("sUnion", sRectangle, sTrd, sRotate, sTrans);
	
	G4Tubs *sTubs = new G4Tubs("sTubs", 15., 19.8, 0.5*mm, -37.3*deg, 74.6*deg);
	sTrans.setX(-15.2);
	G4UnionSolid *sUnion2 = new G4UnionSolid("sUnion2", sUnion, sTubs, 0, sTrans);
	
	
	G4Tubs *sTubs1 = new G4Tubs("sTubs1", 0., 1.0, 1.0*mm, 0., 360.*deg);
	sTrans.setX(-1);
	sTrans.setY(-5.13);
  G4SubtractionSolid *sSub = new G4SubtractionSolid("sSub", sUnion2, sTubs1, 0, sTrans);
	
  G4Tubs *sTubs5 = new G4Tubs("sTubs5", 0., 2.5, 1.0*mm, 0., 360.*deg);
	sTrans.setY(5.13);
  G4SubtractionSolid *sSub2 = new G4SubtractionSolid("sSub2", sSub, sTubs5, 0, sTrans);
	
	// Logical
	G4Material *CollimatorMaterial = G4Material::GetMaterial("Titanium"); 
	fCollimLogical = new G4LogicalVolume(sSub2, CollimatorMaterial, "collimator", 0, 0, 0);
	fCollimLogical->SetVisAttributes(sVisAtt);

}

void ApparatusSpiceTargetChamber::BuildTargetChamberFrontRing()
{

  // Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
  vis_att->SetVisibility(true);

  // ** Dimensions
  // Upstream Lip
  G4double tube_inner_radius = this->front_ring_inner_radius;
  G4double tube_outer_radius = this->front_ring_outer_radius;
  G4double tube_length = this->front_ring_length/2.;
  // Middle Cylinder
  G4double cylinder_outer_radius = this->front_ring_cylinder_outer_radius;
  G4double cylinder_half_length = this->front_ring_cylinder_length/2.;
  // Front Cone
  G4double cone_front_inner = this->cone_front_inner_rad;
  G4double cone_front_outer = this->cone_front_outer_rad;
  G4double cone_half_length = this->cone_length/2.;
  
  // ** ShapesF
  // Upstream Lip
  G4Tubs* upstream_lip = new G4Tubs("upstream_lip", tube_inner_radius, 
				    tube_outer_radius, tube_length,
				    0 , 360*deg);
  // Middle Cylinder
  G4Tubs* middle_cylinder = new G4Tubs("middle_cylinder", tube_inner_radius, 
				       cylinder_outer_radius, cylinder_half_length, 
				       0, 360*deg);
  // Cone
  G4Cons* front_cone = new G4Cons("front_cone", tube_inner_radius, 
				  cylinder_outer_radius, cone_front_inner, 
				  cone_front_outer, cone_half_length, 0., 360*deg);
  
	G4ThreeVector trans(0, 0, tube_length + cylinder_half_length);
	G4UnionSolid* front_ring = new G4UnionSolid("front_ring", upstream_lip, middle_cylinder, 0, trans);

  // ** Logical Volume
  G4Material* front_ring_material = G4Material::GetMaterial(this->downstream_cone_material);
  target_chamber_front_ring_log = new G4LogicalVolume(front_ring,front_ring_material,"target_chamber_front_ring_log",0,0,0);
  target_chamber_front_ring_log->SetVisAttributes(vis_att);
  target_chamber_front_cone_log = new G4LogicalVolume(front_cone,front_ring_material,"target_chamber_front_cone_log",0,0,0);
  target_chamber_front_cone_log->SetVisAttributes(vis_att);

}//end BuildTargetChamberFrontRing

void ApparatusSpiceTargetChamber::BuildTargetChamberSphere()
{

  // ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(TARGET_CHAMBER_COL));
  vis_att->SetVisibility(true);
  
  // ** Dimensions
  // Sphere
  G4double sphere_outer_radius = this->target_chamber_sphere_outer_radius;
  G4double sphere_inner_radius = this->target_chamber_sphere_inner_radius;
  // Upstream Addon Cylinder
  G4double addon_cylinder_length = (this->target_chamber_sphere_centre + this->target_chamber_distance_from_target)/2.;
  // Downstream Cutting Cylinder
  G4double cutting_outer_radius = this->target_chamber_sphere_cutting_outer_radius;
  // Downstream Cutting Box
  G4double cutting_box_half_length = 50*mm;
  G4double cutting_box_from_sphere_centre = this->target_chamber_sphere_cut_off_point;

	// ** Shapes
	// Sphere
	G4Sphere* sphere = new G4Sphere("sphere", sphere_inner_radius, sphere_outer_radius, 0*deg, 360*deg, 0*deg, 90*deg);
	// Upstream Addon Cylinder
	G4Tubs* addon_cylinder = new G4Tubs("addon_cylinder", sphere_inner_radius, sphere_outer_radius, addon_cylinder_length, 0, 360*deg);	
	G4ThreeVector trans(0, 0, -addon_cylinder_length);
	G4UnionSolid* target_chamber_sphere_pre = new G4UnionSolid("target_chamber_sphere_pre", sphere, addon_cylinder, 0, trans);
	// Downstream Cutting Cylinder
	G4Tubs* cutting_cylinder = new G4Tubs("cutting_cylinder", 0, cutting_outer_radius, 200*mm, 0, 360*deg);
	G4SubtractionSolid* target_chamber_sphere_pre2 = new G4SubtractionSolid("target_chamber_sphere_pre2", target_chamber_sphere_pre, cutting_cylinder);
	// Downstream Cutting Box
	G4Box* cutting_box = new G4Box("cutting_box", cutting_box_half_length, cutting_box_half_length, cutting_box_half_length);
	trans.setZ(cutting_box_from_sphere_centre + cutting_box_half_length);
	G4SubtractionSolid* target_chamber_sphere = new G4SubtractionSolid("target_chamber_sphere", target_chamber_sphere_pre2, cutting_box, 0, trans);

  G4Material* sphere_material = G4Material::GetMaterial(this->target_chamber_material);
  target_chamber_sphere_log = new G4LogicalVolume(target_chamber_sphere,sphere_material,"target_chamber_sphere_log",0,0,0);
  target_chamber_sphere_log->SetVisAttributes(vis_att);

}

void ApparatusSpiceTargetChamber::BuildTargetChamberCylinderDownstream()
{

  // ** Visualization
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(TARGET_CHAMBER_COL));
  vis_att->SetVisibility(true);

	// ** Dimensions
  // Main Cylinder
  G4double inner_radius = this->target_chamber_cylinder_inner_radius;
  G4double outer_radius = this->target_chamber_cylinder_outer_radius;
  G4double length = this->target_chamber_cylinder_length/2.;
  // Cylinder Target Lip
  G4double target_lip_inner_radius = this->target_chamber_cylinder_lip_inner_radius;
  G4double target_lip_length = this->target_chamber_cylinder_lip_thickness/2.;
  // Cylinder Detector Lip
  G4double detector_lip_outer_radius = this->target_chamber_cylinder_det_outer_radius;
  G4double detector_lip_length = this->target_chamber_cylinder_det_thickness/2.;  
  // Target Lip Cutting Box
  G4double lip_cutting_half_x = this->target_chamber_cylinder_hole_length/2.;
  G4double lip_cutting_half_y = this->target_chamber_cylinder_hole_breadth/2.;
  G4double lip_cutting_half_z = this->target_chamber_cylinder_lip_thickness;

  // ** Shapes
  G4Tubs* main_cylinder = new G4Tubs("main_cylinder",inner_radius,outer_radius,length,0,360*deg);
  G4Tubs* target_lip_cylinder = new G4Tubs("target_lip_cylinder", target_lip_inner_radius, outer_radius, target_lip_length, 0, 360*deg);
  G4Tubs* detector_lip_cylinder = new G4Tubs("detector_lip_cylinder", outer_radius, detector_lip_outer_radius, detector_lip_length, 0, 360*deg);
  G4Box* target_lip_cutting_box = new G4Box("target_lip_cutting_box", lip_cutting_half_x, lip_cutting_half_y, lip_cutting_half_z);
 
  // Union
  G4ThreeVector trans(0,0, length-target_lip_length);
  G4UnionSolid* target_det_cylinder_pre = new G4UnionSolid("target_det_cylinder_pre", main_cylinder, target_lip_cylinder, 0, trans);
  trans.setZ(-(length-detector_lip_length));
  G4UnionSolid* target_det_cylinder = new G4UnionSolid("target_det_cylinder", target_det_cylinder_pre, detector_lip_cylinder, 0, trans);
  trans.setX(-this->target_chamber_cylinder_hole_distance_x);
  trans.setY(-this->target_chamber_cylinder_hole_distance_y);
  trans.setZ(length-target_lip_length);
  G4RotationMatrix* rotate = new G4RotationMatrix(111.8*deg , 0, 0);
	G4SubtractionSolid* target_det_chamber = new G4SubtractionSolid("target_det_chamber", target_det_cylinder, target_lip_cutting_box, rotate, trans);
  
  // ** Logical Volume
  G4Material* cylinder_material = G4Material::GetMaterial(this->target_chamber_material);
  target_chamber_cylinder_down_log = new G4LogicalVolume(target_det_chamber,cylinder_material,"target_chamber_cylinder_down_log",0,0,0);
  target_chamber_cylinder_down_log->SetVisAttributes(vis_att);

}

void ApparatusSpiceTargetChamber::BuildTargetWheel(){
	
	// ** Visualisation
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
	vis_att->SetVisibility(true);
	
	// ** Dimensions
	G4double wheel_radius = this->target_wheel_radius;
	G4double wheel_half_thickness = this->target_wheel_thickness/2.;
	// Target Cut-Out
	G4double target_radius = this->target_radius;
	G4double target_thickness = this->target_wheel_thickness;
	// Collimator Cut-Out
	G4double collimator_radius = this->collimator_radius;
	
	// ** Shapes
	G4Tubs* target_wheel_pre = new G4Tubs("target_wheel_pre", 0, wheel_radius, wheel_half_thickness, 0, 360*deg);
	G4Tubs* target = new G4Tubs("target", 0, target_radius, target_thickness, 0, 360*deg);
	G4Box* collimator = new G4Box("collimator", 4.5*mm, 8.5*mm, this->target_wheel_thickness);
	
	G4ThreeVector trans(0, 0, 0);
	
	trans.setX(this->target_offset*cos(70*deg));
	trans.setY(this->target_offset*sin(70*deg));
	G4SubtractionSolid* target_wheel0 = new G4SubtractionSolid("target_wheel0", target_wheel_pre, target, 0, trans);
	
	trans.setX(this->target_offset*cos(180*deg));
	trans.setY(this->target_offset*sin(180*deg));
	G4SubtractionSolid* target_wheel1 = new G4SubtractionSolid("target_wheel1", target_wheel0, target, 0, trans);
	
	trans.setX(this->target_offset*cos(-55*deg));
	trans.setY(this->target_offset*sin(-55*deg));
	G4RotationMatrix* sRotate = new G4RotationMatrix(125*deg , 0, 0);
	G4SubtractionSolid* target_wheel = new G4SubtractionSolid("target_wheel", target_wheel1, collimator, sRotate, trans);

	
	// ** Logical
	G4Material* target_wheel_material = G4Material::GetMaterial(this->target_wheel_material); //problem
	target_wheel_log = new G4LogicalVolume(target_wheel, target_wheel_material, "target_wheel_log", 0, 0, 0);
	target_wheel_log->SetVisAttributes(vis_att);

} // end BuildTargetWheel()

void ApparatusSpiceTargetChamber::BuildTargetWheelExtension() {

  // Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(AL_COL));
	sVisAtt->SetVisibility(true);	
	
	// Dimensions
	G4double wheel_radius = this->target_wheel_radius;
	G4double inner_radius = wheel_radius - 5.0*mm;
	G4double wheel_half_thickness = 5.*mm;
	
	// Shapes
	G4Tubs* extension = new G4Tubs("extension", inner_radius, wheel_radius, wheel_half_thickness, 0, 360*deg);

  // Logical
  G4Material* target_wheel_material = G4Material::GetMaterial(this->target_wheel_material);
  fExtLogical = new G4LogicalVolume(extension, target_wheel_material, "fExtLogical", 0, 0, 0);
  fExtLogical->SetVisAttributes(sVisAtt);

}

void ApparatusSpiceTargetChamber::BuildTargetWheelGears(){

	// ** Visualisation
	G4VisAttributes* vis_att_b = new G4VisAttributes(G4Colour(AL_COL));
	vis_att_b->SetVisibility(true);
	G4VisAttributes* vis_att_a = new G4VisAttributes(G4Colour(DELRIN_COL));
	vis_att_a->SetVisibility(true);
	
	// ** Dimensions
	G4double first_gear_radius = this->first_gear_radius;
	G4double second_gear_radius = this->second_gear_radius;
	G4double third_gear_outer_radius = this->third_gear_radius;
	G4double third_gear_inner_radius = this->target_wheel_radius;
	G4double gear_half_thickness = this->all_gear_thickness/2.;
	
	// ** Shapes
	G4Tubs* first_gear = new G4Tubs("first_gear", 0, first_gear_radius, gear_half_thickness, 0, 360*deg);
	G4Tubs* second_gear = new G4Tubs("second_gear", 0, second_gear_radius, gear_half_thickness, 0, 360*deg);
	G4Tubs* third_gear = new G4Tubs("third_gear", third_gear_inner_radius, third_gear_outer_radius, gear_half_thickness, 0, 360*deg);
	
	// ** Logical
	G4Material* gear_material_b = G4Material::GetMaterial(this->target_wheel_gear_material_b);
	first_gear_log = new G4LogicalVolume(first_gear, gear_material_b, "first_gear_log", 0, 0, 0);
	G4Material* gear_material_a = G4Material::GetMaterial(this->target_wheel_gear_material_a);
	second_gear_log = new G4LogicalVolume(second_gear, gear_material_a, "second_gear_log", 0, 0, 0);
	third_gear_log = new G4LogicalVolume(third_gear, gear_material_a, "third_gear_log", 0, 0, 0);
	
	first_gear_log->SetVisAttributes(vis_att_b);
	second_gear_log->SetVisAttributes(vis_att_a);
	third_gear_log->SetVisAttributes(vis_att_a);

} // end BuildTargetWheelGears()

void ApparatusSpiceTargetChamber::BuildTargetWheelGearPlates() {

	// ** Visualisation
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
	vis_att->SetVisibility(true);
	
	// ** Dimensions
	G4double gear_one_radius = this->gear_plate_one_radius;
	G4double gear_two_radius = this->gear_plate_two_radius;
	G4double gear_half_thickness = this->gear_plate_thickness/2.;
	
	// ** Shapes
	G4Tubs* plate_one = new G4Tubs("plate_one", 0, gear_one_radius, gear_half_thickness, 0, 360*deg);
	G4Tubs* plate_two = new G4Tubs("plate_two", 0, gear_two_radius, gear_half_thickness, 0, 360*deg);
	
	// ** Logical
	G4Material* gear_plate_material = G4Material::GetMaterial(this->gear_plate_material);
	gear_plate_one_log = new G4LogicalVolume(plate_one, gear_plate_material, "gear_plate_one_log", 0, 0, 0);
	gear_plate_one_log->SetVisAttributes(vis_att);
	gear_plate_two_log = new G4LogicalVolume(plate_two, gear_plate_material, "gear_plate_two_log", 0, 0, 0);
	gear_plate_two_log->SetVisAttributes(vis_att);
	
} // end::BuildTargetWheelGearPlates()

void ApparatusSpiceTargetChamber::BuildGearStick(){

	// ** Visualisation
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(PEEK_COL));
	vis_att->SetVisibility(true);
	
	// ** Dimensions
	G4double gear_stick_radius = this->gear_stick_radius;
	G4double gear_stick_half_length = this->gear_stick_length/2.;
	
	// ** Shapes
	G4Tubs* gear_stick = new G4Tubs("gear_stick", 0, gear_stick_radius, gear_stick_half_length, 0, 360*deg);
	
	// ** Logical
	G4Material* gear_stick_material = G4Material::GetMaterial(this->gear_stick_material);
	gear_stick_log = new G4LogicalVolume(gear_stick, gear_stick_material, "gear_stick_log", 0, 0, 0);
	gear_stick_log->SetVisAttributes(vis_att);
	
}
	
void ApparatusSpiceTargetChamber::BuildTargetMountPlate(){

	// ** Visualisation
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(PEEK_COL));
	vis_att->SetVisibility(true);
	
	// ** Dimensions
	G4double mount_plate_radius = this->target_mount_plate_radius;
	G4double mount_plate_half_thickness = this->target_mount_plate_thickness/2.;
	// Target Wheel Cut
	G4double wheel_radius = this->target_wheel_radius;
	G4double wheel_half_thickness = this->target_wheel_thickness;
	
	// ** Shapes
	G4Tubs* mount_plate_pre = new G4Tubs("mount_plate_pre", 0, mount_plate_radius, mount_plate_half_thickness, 0, 360*deg);
	
	G4Tubs* target_wheel = new G4Tubs("target_wheel", 0, wheel_radius, 2.*mount_plate_half_thickness, 0, 360*deg);
	G4double offset = this->target_wheel_offset / sqrt(2.);
	G4ThreeVector move(offset, offset, 0);
	G4SubtractionSolid* mount_plate = new G4SubtractionSolid("mount_plate", mount_plate_pre, target_wheel, 0, move);
	
	// ** Logical
	G4Material* target_mount_plate_material = G4Material::GetMaterial(this->target_mount_plate_material);
	target_mount_plate_log = new G4LogicalVolume(mount_plate, target_mount_plate_material, "target_mount_plate_log", 0, 0, 0);
	target_mount_plate_log->SetVisAttributes(vis_att);

} // end BuildTargetMountPlate
	

void ApparatusSpiceTargetChamber::BuildBiasPlate(){

	// ** Visualisation
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
	vis_att->SetVisibility(true);
	
	// ** Dimensions
	G4double plate_outer_radius = this->bias_plate_outer_radius;
	G4double bias_plate_half_thickness = this->bias_plate_thickness/2.;
	// Target Wheel Cut
	G4double wheel_radius = this->target_wheel_radius;
	G4double wheel_half_thickness = this->target_wheel_thickness;
	// Notch
	G4double box_length = this->bias_plate_cut/2.;
	
	// ** Shapes
	G4Tubs* bias_plate_pre = new G4Tubs("bias_plate_pre", 0, plate_outer_radius, bias_plate_half_thickness, 0, 360*deg);
	G4Tubs* target_wheel = new G4Tubs("target_wheel", 0, wheel_radius, wheel_half_thickness, 0, 360*deg);
	G4Box* notch_box = new G4Box("notch_box", box_length, box_length, box_length);	
	
	G4double offset = this->target_wheel_offset / sqrt(2.);
	G4ThreeVector move(offset, offset, 0);
	G4SubtractionSolid* bias_plate_pre2 = new G4SubtractionSolid("bias_plate_pre2", bias_plate_pre, target_wheel, 0, move);
	G4ThreeVector move2(plate_outer_radius, plate_outer_radius, 0);
	G4SubtractionSolid* bias_plate = new G4SubtractionSolid("bias_plate", bias_plate_pre2, notch_box, 0, move2);
	
	// ** Logical
	G4Material* bias_plate_material = G4Material::GetMaterial(this->bias_plate_material);
	bias_plate_log = new G4LogicalVolume(bias_plate, bias_plate_material, "bias_plate_log", 0, 0, 0);
}

void ApparatusSpiceTargetChamber::BuildPhotonShield()
{
	// ** Visualization
  G4VisAttributes* one_vis_att = new G4VisAttributes(G4Colour(PB_COL));
  G4VisAttributes* two_vis_att = new G4VisAttributes(G4Colour(SN_COL));
  G4VisAttributes* three_vis_att = new G4VisAttributes(G4Colour(CU_COL));
  one_vis_att->SetVisibility(true);
  two_vis_att->SetVisibility(true);
  three_vis_att->SetVisibility(true);
  
  // ** Build Photon Shield Solid (First Layer)
 // Define Variables
  G4double inner_radius = this->photon_shield_inner_radius;
  G4double total_length = this->photon_shield_layer_one_thickness + this->photon_shield_layer_two_thickness + this->photon_shield_layer_three_thickness;
  G4double half_length = total_length/2.;

  G4double back_outer_radius = this->photon_shield_back_radius;
  G4double front_outer_radius = this->photon_shield_front_radius;

  G4Cons* photon_shield = new G4Cons("photon_shield", inner_radius, back_outer_radius, 
				     inner_radius, front_outer_radius,	half_length, 0, 2*pi);
  
  // ** Build Magnet Clamp
  G4double clamp_half_thickness = this->plate_one_thickness/2. + this->psclamp_thickness_buffer;
  G4Box* magnet_clamp = new G4Box("magnet_clamp", this->psclamp_length/2., clamp_half_thickness, this->photon_shield_length);
  
  
  // ** Cut Plate One from Photon Shield
  G4RotationMatrix* cutRot = new G4RotationMatrix; 
  G4double transRadial = this->plate_one_edge_x + this->psclamp_length/2. - this->psclamp_psbuffer;
  G4double transX, transY;
  G4double transZ = 0;
  G4ThreeVector trans;

  cutRot->rotateZ(-(45.+0*90.)*deg);
  transX = transRadial*cos((0*90.+45.)*deg);
  transY = transRadial*sin((0*90.+45.)*deg);
  trans.set(transX,transY,transZ);
  G4SubtractionSolid* photon_shield_temp_1 = new G4SubtractionSolid("photon_shield",photon_shield,magnet_clamp,cutRot,trans);      
  cutRot->rotateZ(-90.*deg);
  transX = transRadial*cos((1*90.+45.)*deg);
  transY = transRadial*sin((1*90.+45.)*deg);
  trans.set(transX,transY,transZ);
  G4SubtractionSolid* photon_shield_temp_2 = new G4SubtractionSolid("photon_shield",photon_shield_temp_1,magnet_clamp,cutRot,trans);      
  cutRot->rotateZ(-90.*deg);
  transX = transRadial*cos((2*90.+45.)*deg);
  transY = transRadial*sin((2*90.+45.)*deg);
  trans.set(transX,transY,transZ);
  G4SubtractionSolid* photon_shield_temp_3 = new G4SubtractionSolid("photon_shield",photon_shield_temp_2,magnet_clamp,cutRot,trans);      
  cutRot->rotateZ(-90.*deg);
  transX = transRadial*cos((3*90.+45.)*deg);
  transY = transRadial*sin((3*90.+45.)*deg);
  trans.set(transX,transY,transZ);
  G4SubtractionSolid* photon_shield_temp_4 = new G4SubtractionSolid("photon_shield",photon_shield_temp_3,magnet_clamp,cutRot,trans); 
  
  // ** Bolt Holes
  G4double target_bolt_half_length = this->ps_target_bolt_length/2.;
  G4double target_bolt_radius = this->ps_target_bolt_radius;
  G4double target_bolt_offset = this->ps_target_bolt_plane_offset / sqrt(2.);
  G4double target_bolt_z_offset = this->photon_shield_length/2. - target_bolt_half_length;
  G4double detector_bolt_half_length = this->ps_det_bolt_length/2.;
  G4double detector_bolt_radius = this->ps_det_bolt_radius;
  G4double detector_bolt_offset = this->ps_det_bolt_plane_offset / sqrt(2.);
  G4double detector_bolt_z_offset = -this->photon_shield_length/2. + detector_bolt_half_length;
  
  G4Tubs* target_bolt = new G4Tubs("target_bolt", 0, target_bolt_radius, target_bolt_half_length, 0, 360*deg);
  G4Tubs* detector_bolt = new G4Tubs("detector_bolt", 0, detector_bolt_radius, detector_bolt_half_length, 0, 360*deg);
  
  // the four front (target end) bolts are created using a union solid then deducted from the previous volume
  G4double distance_between_target_bolts = this->ps_target_bolt_plane_offset * sqrt(2.);
  G4ThreeVector bolt_distance(distance_between_target_bolts, 0, 0);
  G4UnionSolid* target_bolts_two = new G4UnionSolid("target_bolts_two", target_bolt, target_bolt, 0, bolt_distance);
  G4ThreeVector bolt_dist2(0, distance_between_target_bolts, 0);
  G4UnionSolid* target_bolts = new G4UnionSolid("target_bolts", target_bolts_two, target_bolts_two, 0, bolt_dist2);
  G4ThreeVector bolt_move(-target_bolt_offset, -target_bolt_offset, target_bolt_z_offset);
  G4SubtractionSolid* photon_shield_temp_5 = new G4SubtractionSolid("photon_shield", photon_shield_temp_4, target_bolts, 0, bolt_move);
  
  G4double distance_between_det_bolts = this->ps_det_bolt_plane_offset * sqrt(2.);
  G4ThreeVector bolt_dist3(distance_between_det_bolts, 0, 0);
  G4UnionSolid* det_bolts_two = new G4UnionSolid("det_bolts_two", detector_bolt, detector_bolt, 0, bolt_dist3); 
  G4ThreeVector bolt_dist4(0, distance_between_det_bolts, 0);
  G4UnionSolid* det_bolts = new G4UnionSolid("det_bolts", det_bolts_two, det_bolts_two, 0, bolt_dist4);
  G4ThreeVector bolt_move2(-detector_bolt_offset, -detector_bolt_offset, detector_bolt_z_offset);
  G4SubtractionSolid* photon_shield_temp_6 = new G4SubtractionSolid("photon_shield", photon_shield_temp_5, det_bolts, 0, bolt_move2);
  

  // photon_shield_temp_12 can be considered a template.  To segment the shield into layers I will create a cutting box
  G4Box* photon_shield_slicer = new G4Box("photon_shield_slicer",50*mm, 50*mm, 50*mm);

  // Cut out the first layar
  G4ThreeVector slicer(0, 0, -(50+half_length-(this->photon_shield_layer_three_thickness+this->photon_shield_layer_two_thickness)));
  G4SubtractionSolid* first_layer = new G4SubtractionSolid("first_layer", photon_shield_temp_6, photon_shield_slicer, 0, slicer);

  // Cut out a second layer
  slicer.set(0,0, (50+half_length-this->photon_shield_layer_one_thickness));
  G4SubtractionSolid* second_pre = new G4SubtractionSolid("second_pre", photon_shield_temp_6, photon_shield_slicer, 0, slicer);
  slicer.set(0, 0, -(50+half_length-this->photon_shield_layer_three_thickness));
  G4SubtractionSolid* second_layer = new G4SubtractionSolid("second_layer", second_pre, photon_shield_slicer, 0, slicer);

  // Cut out a third layer
  slicer.set(0, 0, (50+half_length-(this->photon_shield_layer_one_thickness+this->photon_shield_layer_two_thickness)));
  G4SubtractionSolid* third_layer = new G4SubtractionSolid("third_layer", photon_shield_temp_6, photon_shield_slicer, 0, slicer);
 
  // Create the three separate logs 
  G4Material* layer_one_material = G4Material::GetMaterial(this->photon_shield_layer_one_material);
  photon_shield_layer_one_log = new G4LogicalVolume(first_layer, layer_one_material,"photon_shield_layer_one_log",0,0,0);
  photon_shield_layer_one_log->SetVisAttributes(one_vis_att);
  G4Material* layer_two_material = G4Material::GetMaterial(this->photon_shield_layer_two_material);
  photon_shield_layer_two_log = new G4LogicalVolume(second_layer, layer_two_material, "photon_shield_layer_two_log",0,0,0);
  photon_shield_layer_two_log->SetVisAttributes(two_vis_att);
  G4Material* layer_three_material = G4Material::GetMaterial(this->photon_shield_layer_three_material);
  photon_shield_layer_three_log = new G4LogicalVolume(third_layer, layer_three_material, "photon_shield_layer_three_log",0,0,0);
  photon_shield_layer_three_log->SetVisAttributes(three_vis_att);
  
} // end BuildPhotonShield()

void ApparatusSpiceTargetChamber::BuildShieldCovering()
{
	// ** Visualization
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(KAPTON_COL));
  vis_att->SetVisibility(true);
  
  // ** Build Photon Shield Solid
  //G4double inner_radius = this->photon_shield_inner_radius; 
  G4double front_outer_radius = this->photon_shield_front_radius;
  G4double back_outer_radius = this->photon_shield_back_radius;
  G4double half_length = this->photon_shield_length/2.;  
  G4Cons* photon_shield = new G4Cons("photon_shield", 0, back_outer_radius, 0, front_outer_radius,	half_length, 0, 2*pi);

	// ------ 																							
	// Photon Shield Clamps
	// -- Dimensions  																							
  G4double detector_end_radius = this->ps_det_clamp_radius;
  G4double detector_end_half = this->ps_det_clamp_thickness/2.;
  G4double target_end_radius = this->ps_target_clamp_radius;
  G4double target_end_half = this->ps_target_clamp_thickness/2.;				
  // -- Shapes																			
	G4Tubs* target_clamp = new G4Tubs("target_clamp", 0, target_end_radius, target_end_half, 0, 360*deg);
	G4Tubs* det_clamp = new G4Tubs("det_clamp", 0, detector_end_radius, detector_end_half, 0, 360*deg);
	// -- Union
	G4ThreeVector move(0, 0, half_length + target_end_half - 0.01*mm);
	G4UnionSolid* photon_shield_2 = new G4UnionSolid("photon_shield_2", photon_shield, target_clamp, 0, move);
	move.setZ(-half_length);
	G4UnionSolid* photon_shield_3 = new G4UnionSolid("photon_shield_3", photon_shield_2, det_clamp, 0, move);
	// ------
	
																								
  // ** Build Shield Covering
  G4double theta = atan(2 * half_length / (back_outer_radius - front_outer_radius));
  G4double alpha = 90*deg - theta;
  G4double cover_front_radius = front_outer_radius - (this->shield_coating_thickness / tan(theta)) +
  																							(this->shield_coating_thickness * cos(alpha));
  G4double cover_back_radius = back_outer_radius + (this->shield_coating_thickness / tan(theta)) +
  																							(this->shield_coating_thickness * cos(alpha));

  G4Cons* shield_cover_pre = new G4Cons("shield_cover_pre", 0, cover_back_radius, 0, cover_front_radius,	
  																							half_length + this->shield_coating_thickness, 0, 2*pi);
  																							
	// ** Deduct Shield from Cover
	G4SubtractionSolid* shield_cover = new G4SubtractionSolid("shield_cover", shield_cover_pre, photon_shield_3);
	
	// ------
	// Magnet Clamps
	// -- Dimensions
	G4double clamp_half_thickness = this->plate_one_thickness/2. + this->psclamp_thickness_buffer;
	G4Box* magnet_clamp = new G4Box("magnet_clamp", this->psclamp_length/2., clamp_half_thickness, this->photon_shield_length);
	// -- Rotation & Translation
	G4RotationMatrix* cutRot = new G4RotationMatrix; 
  G4double transRadial = this->plate_one_edge_x + this->psclamp_length/2. - this->psclamp_psbuffer;
  G4double transX, transY;
  G4double transZ = 0;
  G4ThreeVector trans;
	// -- Union
  cutRot->rotateZ(-(45.+0*90.)*deg);
  transX = transRadial*cos((0*90.+45.)*deg);
  transY = transRadial*sin((0*90.+45.)*deg);
  trans.set(transX,transY,transZ);
  G4SubtractionSolid* shield_cover_2 = new G4SubtractionSolid("photon_shield", shield_cover,magnet_clamp,cutRot,trans);      
  cutRot->rotateZ(-90.*deg);
  transX = transRadial*cos((1*90.+45.)*deg);
  transY = transRadial*sin((1*90.+45.)*deg);
  trans.set(transX,transY,transZ);
  G4SubtractionSolid* shield_cover_3 = new G4SubtractionSolid("photon_shield",shield_cover_2,magnet_clamp,cutRot,trans);      
  cutRot->rotateZ(-90.*deg);
  transX = transRadial*cos((2*90.+45.)*deg);
  transY = transRadial*sin((2*90.+45.)*deg);
  trans.set(transX,transY,transZ);
  G4SubtractionSolid* shield_cover_4 = new G4SubtractionSolid("photon_shield",shield_cover_3,magnet_clamp,cutRot,trans);      
  cutRot->rotateZ(-90.*deg);
  transX = transRadial*cos((3*90.+45.)*deg);
  transY = transRadial*sin((3*90.+45.)*deg);
  trans.set(transX,transY,transZ);
  G4SubtractionSolid* shield_cover_5 = new G4SubtractionSolid("photon_shield",shield_cover_4,magnet_clamp,cutRot,trans); 
  // -------
	
	G4Material* shield_cover_material = G4Material::GetMaterial(this->shield_cover_material);
  shield_cover_log = new G4LogicalVolume(shield_cover_5, shield_cover_material, "shield_cover_log", 0, 0, 0);
  shield_cover_log->SetVisAttributes(vis_att);

}

void ApparatusSpiceTargetChamber::BuildPhotonShieldClamps() {

	// ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
  vis_att->SetVisibility(true);
  
  // ** Dimensions
  G4double detector_end_radius = this->ps_det_clamp_radius;
  G4double detector_end_half = this->ps_det_clamp_thickness/2.;
  G4double target_end_radius = this->ps_target_clamp_radius;
  G4double target_end_half = this->ps_target_clamp_thickness/2.;
  G4double inner_radius = this->photon_shield_inner_radius;
  // Extensions
  G4double ext_half_length = this->ps_target_clamp_radius + this->ps_target_clamp_extension;
  G4double ext_half_width = this->plate_one_thickness/2. + this->psclamp_thickness_buffer;
  // Bolts
  G4double target_bolt_radius = this->ps_target_bolt_radius;
  G4double target_bolt_offset = this->ps_target_bolt_plane_offset;
  G4double detector_bolt_radius = this->ps_det_bolt_radius;
  G4double detector_bolt_offset = this->ps_det_bolt_plane_offset;
  
  
  // ** Shapes
  G4Tubs* target_end_clamp_pre = new G4Tubs("target_end_clamp_pre", inner_radius, target_end_radius, target_end_half, 0, 360*deg);
  G4Box* target_clamp_ext = new G4Box("target_clamp_ext", ext_half_length, ext_half_width, target_end_half);
  G4Box* target_clamp_ext2 = new G4Box("target_clamp_ext2", ext_half_width, ext_half_length, target_end_half);
  G4UnionSolid* target_extensions = new G4UnionSolid("target_extensions", target_clamp_ext, target_clamp_ext2);
  G4Tubs* extension_cut = new G4Tubs("extension_cut", 0, inner_radius + 2*mm, 2*target_end_half, 0, 360*deg);
  G4SubtractionSolid* target_ext_hole = new G4SubtractionSolid("target_ext_hole", target_extensions, extension_cut);  
  G4UnionSolid* target_clamp_pre_bolts = new G4UnionSolid("target_clamp_pre_bolts", target_end_clamp_pre, target_ext_hole);  
  
  G4Tubs* target_bolt = new G4Tubs("target_bolt", 0, target_bolt_radius, 2.*target_end_half, 0, 360*deg);
  G4ThreeVector move_bolt(target_bolt_offset, 0, 0);
  G4SubtractionSolid* target_end_clamp0 = new G4SubtractionSolid("target_end_clamp0", target_clamp_pre_bolts, target_bolt, 0, move_bolt);
  move_bolt.setX(-target_bolt_offset);
  G4SubtractionSolid* target_end_clamp1 = new G4SubtractionSolid("target_end_clamp1", target_end_clamp0, target_bolt, 0, move_bolt);
  move_bolt.setX(0);
  move_bolt.setY(target_bolt_offset);
  G4SubtractionSolid* target_end_clamp2 = new G4SubtractionSolid("target_end_clamp2", target_end_clamp1, target_bolt, 0, move_bolt);
  move_bolt.setY(-target_bolt_offset);
  G4SubtractionSolid* target_end_clamp3 = new G4SubtractionSolid("target_end_clamp3", target_end_clamp2, target_bolt, 0, move_bolt);
  
  G4Tubs* det_end_clamp_pre = new G4Tubs("det_end_clamp_pre", inner_radius, detector_end_radius, detector_end_half, 0, 360*deg);
  G4Tubs* detector_bolt = new G4Tubs("detector_bolt", 0, detector_bolt_radius, 2.* detector_end_half, 0, 360*deg);
	move_bolt.setY(detector_bolt_offset);
	G4SubtractionSolid* det_end_clamp0 = new G4SubtractionSolid("det_end_clamp0", det_end_clamp_pre, detector_bolt, 0, move_bolt);
	move_bolt.setY(-detector_bolt_offset);
	G4SubtractionSolid* det_end_clamp1 = new G4SubtractionSolid("det_end_clamp1", det_end_clamp0, detector_bolt, 0, move_bolt);
	move_bolt.setY(0);
	move_bolt.setX(detector_bolt_offset);
	G4SubtractionSolid* det_end_clamp2 = new G4SubtractionSolid("det_end_clamp2", det_end_clamp1, detector_bolt, 0, move_bolt);
	move_bolt.setX(-detector_bolt_offset);
	G4SubtractionSolid* det_end_clamp3 = new G4SubtractionSolid("det_end_clamp3", det_end_clamp2, detector_bolt, 0, move_bolt);
	
  
  // ** Logical
  G4Material* ps_clamp_material = G4Material::GetMaterial(this->ps_clamp_material);
  ps_target_clamp_log = new G4LogicalVolume(target_end_clamp3, ps_clamp_material, "ps_target_clamp_log", 0, 0, 0);
  ps_target_clamp_log->SetVisAttributes(vis_att);
  ps_detector_clamp_log = new G4LogicalVolume(det_end_clamp3, ps_clamp_material, "ps_detector_clamp_log", 0, 0, 0);
  ps_detector_clamp_log->SetVisAttributes(vis_att);
  
} // end::BuildPhotonShieldClamps

void ApparatusSpiceTargetChamber::BuildPhotonShieldClampBolts() {
	
	// ** Visualisation
  G4VisAttributes* small_vis_att = new G4VisAttributes(G4Colour(CU_COL));
  small_vis_att->SetVisibility(true);
  G4VisAttributes* large_vis_att = new G4VisAttributes(G4Colour(SN_COL));
  large_vis_att->SetVisibility(true);
  
  // ** Dimensions
  G4double target_bolt_half_length = this->ps_target_bolt_length/2.;
  G4double target_bolt_radius = this->ps_target_bolt_radius;
  G4double target_bolt_head_half = this->ps_target_bolt_head_length/2.;
  G4double target_bolt_head_rad = this->ps_target_bolt_head_radius;
  G4double detector_bolt_half_length = this->ps_det_bolt_length/2. + //this->pipe_z_length/2. + this->ps_det_clamp_thickness/2.;
                                                     this->ps_det_clamp_thickness/2.;
  G4double detector_bolt_radius = this->ps_det_bolt_radius;
  G4double detector_head_half = this->ps_det_bolt_head_length/2.;
  G4double detector_head_radius = this->ps_det_bolt_head_radius;
  
  // ** Shapes
  G4Tubs* target_bolt_shaft = new G4Tubs("target_bolt_shaft", 0, target_bolt_radius, target_bolt_half_length, 0, 360*deg);
  G4Tubs* target_bolt_head = new G4Tubs("target_bolt_head", 0, target_bolt_head_rad, target_bolt_head_half, 0, 360*deg);
	G4ThreeVector trans(0, 0, target_bolt_half_length + target_bolt_head_half);
  G4UnionSolid* target_bolt = new G4UnionSolid("target_bolt", target_bolt_shaft, target_bolt_head, 0, trans);
  G4Tubs* detector_bolt_shaft = new G4Tubs("detector_bolt_shaft", 0, detector_bolt_radius, detector_bolt_half_length, 0, 360*deg);
  G4Tubs* detector_bolt_head = new G4Tubs("detector_bolt_head", 0, detector_head_radius, detector_head_half, 0, 360*deg);
  trans.setZ(-detector_bolt_half_length - detector_head_half);
  G4UnionSolid* detector_bolt = new G4UnionSolid("detector_bolt", detector_bolt_shaft, detector_bolt_head, 0, trans);
  
  // ** Logical
  G4Material* small_bolt_material = G4Material::GetMaterial(this->small_bolt_material);
  target_bolt_log = new G4LogicalVolume(target_bolt, small_bolt_material, "target_bolt_log", 0, 0, 0);
  target_bolt_log->SetVisAttributes(small_vis_att);
  G4Material* large_bolt_material = G4Material::GetMaterial(this->large_bolt_material);
  detector_bolt_log = new G4LogicalVolume(detector_bolt, large_bolt_material, "detector_bolt_log", 0, 0, 0);
  detector_bolt_log->SetVisAttributes(large_vis_att);
  
} // end::BuildPhotonShieldClampBolts()
  


void ApparatusSpiceTargetChamber::BuildCollectorMagnet()
{

  // ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(NDFEB_COL));
  vis_att->SetVisibility(true);
  
	// ** Build Plate One
	G4Box* plate_one_pre_cut = new G4Box("plate_one_pre_cut", this->plate_one_length/2., 
																						this->plate_one_thickness/2., this->plate_one_height/2.);
	G4Box* plate_one_cut_box = new G4Box("plate_one_cut_box", this->plate_one_length/2., 
																						this->plate_one_thickness/2. + 1., this->plate_one_height/2.);
	
	// ( Cut box is rotated and moved, exact movement calculated below )
	G4double angle_to_corner = atan(this->plate_one_length/this->plate_one_height);
	G4double hypotenuse_to_corner = sqrt(pow(this->plate_one_height/2., 2) + pow(this->plate_one_length/2., 2));
	G4double angle_difference, translate_cut_x, translate_cut_z;
	if (angle_to_corner < this->cutting_box_angle)
	{
		angle_difference = this->cutting_box_angle - angle_to_corner;
		translate_cut_x = (hypotenuse_to_corner * sin(angle_difference)) + this->plate_one_length/2.;
		translate_cut_z = (hypotenuse_to_corner * cos(angle_difference)) - this->plate_one_height/2.
																						+ this->plate_one_lower_height;
	}
	else if (angle_to_corner > this->cutting_box_angle)
	{
		angle_difference = angle_to_corner - this->cutting_box_angle;
		translate_cut_x = this->plate_one_length/2. - (hypotenuse_to_corner * sin(angle_difference));
		translate_cut_z = (hypotenuse_to_corner * cos(angle_difference)) - this->plate_one_height/2.
																						+ this->plate_one_lower_height;
	}
	G4ThreeVector translate_cut_box(-translate_cut_x, 0, translate_cut_z);
	G4RotationMatrix* rotate_cut_box = new G4RotationMatrix;
  rotate_cut_box->rotateY(this->cutting_box_angle);
	G4SubtractionSolid* plate_one = new G4SubtractionSolid("plate_one", plate_one_pre_cut, plate_one_cut_box,
																						rotate_cut_box, translate_cut_box);
																						
	// ** Build Plate Two
	G4Box* plate_two = new G4Box("plate_two", this->plate_two_length/2., 
																						(this->plate_two_thickness + this->plate_one_thickness/2.),
																						this->plate_two_height/2.);
	
	// ** Combine Plates
	G4ThreeVector translate_plate_two(this->plate_one_length/2. - this->plate_two_length/2., 0, 0);
	G4UnionSolid* plate_combo = new G4UnionSolid("plate_combo", plate_one, plate_two, 0, translate_plate_two);
	  
  // ** Logical Volume
  G4Material* magnet_material = G4Material::GetMaterial(this->magnet_material);
  magnet_log = new G4LogicalVolume(plate_combo, magnet_material, "magnet_log", 0, 0, 0);
  magnet_log->SetVisAttributes(vis_att);

} // end:BuildCollectorMagnet()

void ApparatusSpiceTargetChamber::BuildMagnetCovering()
{
	// ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(PEEK_COL));
  vis_att->SetVisibility(true);
  
	// ** Build Plate One
	G4Box* plate_one_pre_cut = new G4Box("plate_one_pre_cut", this->plate_one_length/2., 
																						this->plate_one_thickness/2., this->plate_one_height/2.);
	G4Box* plate_one_cut_box = new G4Box("plate_one_cut_box", this->plate_one_length/2., 
																						this->plate_one_thickness/2. + 1., this->plate_one_height/2.);
	
	// ( Cut box is rotated and moved, exact movement calculated below )
	G4double angle_to_corner = atan(this->plate_one_length/this->plate_one_height);
	G4double hypotenuse_to_corner = sqrt(pow(this->plate_one_height/2., 2) + pow(this->plate_one_length/2., 2));
	G4double angle_difference, translate_cut_x, translate_cut_z;
	if (angle_to_corner < this->cutting_box_angle)
	{
		angle_difference = this->cutting_box_angle - angle_to_corner;
		translate_cut_x = (hypotenuse_to_corner * sin(angle_difference)) + this->plate_one_length/2.;
		translate_cut_z = (hypotenuse_to_corner * cos(angle_difference)) - this->plate_one_height/2.
																						+ this->plate_one_lower_height;
	}
	else if (angle_to_corner > this->cutting_box_angle)
	{
		angle_difference = angle_to_corner - this->cutting_box_angle;
		translate_cut_x = this->plate_one_length/2. - (hypotenuse_to_corner * sin(angle_difference));
		translate_cut_z = (hypotenuse_to_corner * cos(angle_difference)) - this->plate_one_height/2.
																						+ this->plate_one_lower_height;
	}
	G4ThreeVector translate_cut_box(-translate_cut_x, 0, translate_cut_z);
	G4RotationMatrix* rotate_cut_box = new G4RotationMatrix;
  rotate_cut_box->rotateY(this->cutting_box_angle);
	G4SubtractionSolid* plate_one = new G4SubtractionSolid("plate_one", plate_one_pre_cut, plate_one_cut_box,
																						rotate_cut_box, translate_cut_box);
																						
	// ** Build Plate Two
	G4Box* plate_two = new G4Box("plate_two", this->plate_two_length/2., 
																						(this->plate_two_thickness + this->plate_one_thickness/2.),
																						this->plate_two_height/2.);
	
	// ** Combine Plates
	G4ThreeVector translate_plate_two(this->plate_one_length/2. - this->plate_two_length/2., 0, 0);
	G4UnionSolid* plate_combo = new G4UnionSolid("plate_combo", plate_one, plate_two, 0, translate_plate_two);
	
	// ------------
	// Photon Shield Clamp
	// -- Dimensions
  G4double clamp_half_thickness = this->plate_one_thickness/2. + this->psclamp_thickness_buffer;
  G4double clamp_half_length = (this->plate_one_length - this->plate_two_length + this->psclamp_psbuffer + this->psclamp_chamber_buffer)/2.;
  G4double clamp_half_height = this->photon_shield_length/2.;
  // -- Shapes
  G4Box* clamp_main = new G4Box("clamp_main", clamp_half_length, clamp_half_thickness, clamp_half_height);
  // -- Union
  G4double xtrans = (this->plate_one_length + this->psclamp_length)/2. - this->psclamp_length + this->psclamp_psbuffer;
	G4double ztrans = -this->photon_shield_back_face_pos - this->plate_one_height/2. - this->distance_from_target - this->photon_shield_length/2.;
  G4ThreeVector move_clamp(-xtrans, 0, -ztrans);
  G4UnionSolid* plate_and_clamp = new G4UnionSolid("plate_and_clamp", plate_combo, clamp_main, 0, move_clamp);
  // --------------
  
  
  // --------------
  // Chamber Clamps
  // -- Dimensions
  clamp_half_thickness = this->mclamp_chamber_thickness/2.;
  clamp_half_length = this->mclamp_chamber_length/2.;
  clamp_half_height = this->mclamp_chamber_height/2.;
  // -- Shape
  G4Box* clamp_main2 = new G4Box("clamp_main", clamp_half_length, clamp_half_thickness, clamp_half_height);
  xtrans = -(clamp_half_length + this->plate_one_length/2.) + 
  					(this->mclamp_chamber_length - (this->target_chamber_cylinder_inner_radius - (this->plate_one_length + this->plate_one_edge_x)));
  G4ThreeVector move_again(-xtrans, 0, -clamp_half_height + this->plate_one_height/2.);
  G4UnionSolid* total_deduction = new G4UnionSolid("total_deduction", plate_and_clamp, clamp_main2, 0, move_again);
  // ---------------
  
	
	// ------------------
	// Build Magnet Cover
	
	// ** Cover of Plate One
	G4Box* plate_one_cover_pre_cut = new G4Box("plate_one_cover_pre_cut", this->plate_one_length/2. + this->magnet_coating_thickness, 
																						this->plate_one_thickness/2. + this->magnet_coating_thickness, 
																						this->plate_one_height/2. + this->magnet_coating_thickness);
	G4Box* plate_one_cover_cut_box = new G4Box("plate_one__cover_cut_box", 
																						this->plate_one_length/2., 
																						this->plate_one_thickness + this->magnet_coating_thickness, 
																						this->plate_one_height/2. - this->magnet_coating_thickness);
	G4SubtractionSolid* plate_one_cover = new G4SubtractionSolid("plate_one_cover", plate_one_cover_pre_cut, 
																						plate_one_cover_cut_box, rotate_cut_box, translate_cut_box); 
																						
	// ** Build Plate Two Cover
	G4Box* plate_two_cover = new G4Box("plate_two_cover", this->plate_two_length/2. + this->magnet_coating_thickness, 
																						this->plate_two_thickness + this->plate_one_thickness/2. + this->magnet_coating_thickness,
																						this->plate_two_height/2. + this->magnet_coating_thickness);
																						
	// ** Combine Plate Covers
	G4UnionSolid* cover_combo = new G4UnionSolid("cover_combo", plate_one_cover, plate_two_cover, 0, translate_plate_two);
	
	// ---------------------
	// Deduct Magnet from Cover
	G4SubtractionSolid* magnet_cover = new G4SubtractionSolid("magnet_cover", cover_combo, total_deduction);
	
  // ** Logical Volume
  G4Material* magnet_cover_material = G4Material::GetMaterial(this->magnet_cover_material);
  magnet_cover_log = new G4LogicalVolume(magnet_cover, magnet_cover_material, "magnet_cover_log", 0, 0, 0);
  magnet_cover_log->SetVisAttributes(vis_att);

}

void ApparatusSpiceTargetChamber::BuildMagnetClampChamber(){

	// ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
  vis_att->SetVisibility(true);
  
  // ** Dimensions
  G4double clamp_half_thickness = this->mclamp_chamber_thickness/2.;
  G4double clamp_half_length = this->mclamp_chamber_length/2.;
  G4double clamp_half_height = this->mclamp_chamber_height/2.;
  
  G4double magnet_half_thickness = (this->plate_one_thickness + 2.*this->plate_two_thickness)/2.;
  G4double magnet_half_length = this->plate_two_length/2.;
  G4double magnet_half_height = this->plate_two_height/2.;
  
  // ** Shapes
  G4Box* clamp_main = new G4Box("clamp_main", clamp_half_length, clamp_half_thickness, clamp_half_height);
  G4Box* magnet_cut = new G4Box("magnet_cut", magnet_half_length, magnet_half_thickness, magnet_half_height);
  G4double xtrans = -(clamp_half_length + magnet_half_length) + 
  					(this->mclamp_chamber_length - (this->target_chamber_cylinder_inner_radius - (this->plate_one_length + this->plate_one_edge_x)));
  G4ThreeVector trans(xtrans, 0, clamp_half_height - magnet_half_height);
  G4SubtractionSolid* magnet_clamp = new G4SubtractionSolid("magnet_clamp", clamp_main, magnet_cut, 0, trans);
  
  // ** Logical
  G4Material* clamp_material = G4Material::GetMaterial(this->magnet_clamp_material);
  mclamp_chamber_log = new G4LogicalVolume(magnet_clamp, clamp_material, "magnet_clamp_chamber_log", 0, 0, 0);
  mclamp_chamber_log->SetVisAttributes(vis_att);


} // end::BuildMagnetClampChamber()

void ApparatusSpiceTargetChamber::BuildMagnetClampPhotonShield(){

	// ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
  vis_att->SetVisibility(true);
  
  // ** Dimensions
  G4double clamp_half_thickness = this->plate_one_thickness/2. + this->psclamp_thickness_buffer;
  G4double clamp_half_length = (this->plate_one_length - this->plate_two_length + this->psclamp_psbuffer + this->psclamp_chamber_buffer)/2.;
  G4double clamp_half_height = this->photon_shield_length/2.;
  
  // ** Shapes
  G4Box* clamp_main = new G4Box("clamp_main", clamp_half_length, clamp_half_thickness, clamp_half_height);
  
  // ------------------------
  // Rebuild Magnets
  // ** Build Plate One
  G4Box* plate_one_pre_cut = new G4Box("plate_one_pre_cut", this->plate_one_length/2., 
				       this->plate_one_thickness/2., this->plate_one_height/2.);
  G4Box* plate_one_cut_box = new G4Box("plate_one_cut_box", this->plate_one_length/2., 
				       this->plate_one_thickness/2. + 1., this->plate_one_height/2.);
  
  // ( Cut box is rotated and moved, exact movement calculated below )
  G4double angle_to_corner = atan(this->plate_one_length/this->plate_one_height);
  G4double hypotenuse_to_corner = sqrt(pow(this->plate_one_height/2., 2) + pow(this->plate_one_length/2., 2));
  G4double angle_difference, translate_cut_x, translate_cut_z;
  if (angle_to_corner < this->cutting_box_angle)
    {
      angle_difference = this->cutting_box_angle - angle_to_corner;
      translate_cut_x = (hypotenuse_to_corner * sin(angle_difference)) + this->plate_one_length/2.;
      translate_cut_z = (hypotenuse_to_corner * cos(angle_difference)) - this->plate_one_height/2.
	+ this->plate_one_lower_height;
    }
  else if (angle_to_corner > this->cutting_box_angle)
    {
      angle_difference = angle_to_corner - this->cutting_box_angle;
      translate_cut_x = this->plate_one_length/2. - (hypotenuse_to_corner * sin(angle_difference));
      translate_cut_z = (hypotenuse_to_corner * cos(angle_difference)) - this->plate_one_height/2.
	+ this->plate_one_lower_height;
    }
  G4ThreeVector translate_cut_box(-translate_cut_x, 0, translate_cut_z);
  G4RotationMatrix* rotate_cut_box = new G4RotationMatrix;
  rotate_cut_box->rotateY(this->cutting_box_angle);
  G4SubtractionSolid* plate_one = new G4SubtractionSolid("plate_one", plate_one_pre_cut, plate_one_cut_box,
							 rotate_cut_box, translate_cut_box);
  
  // ** Build Plate Two
  G4Box* plate_two = new G4Box("plate_two", this->plate_two_length/2., 
			       (this->plate_two_thickness + this->plate_one_thickness/2.),
			       this->plate_two_height/2.);
  
  // ** Combine Plates
  G4ThreeVector translate_plate_two(this->plate_one_length/2. - this->plate_two_length/2., 0, 0);
  G4UnionSolid* plate_combo = new G4UnionSolid("plate_combo", plate_one, plate_two, 0, translate_plate_two);
  // ---------------------------------
  
  G4double xtrans = (this->plate_one_length + this->psclamp_length)/2. - this->psclamp_length + this->psclamp_psbuffer;
  G4double ztrans = -this->photon_shield_back_face_pos - this->plate_one_height/2. - this->distance_from_target - this->photon_shield_length/2.;
  G4ThreeVector trans(xtrans, 0, ztrans);
  G4SubtractionSolid* magnet_clamp = new G4SubtractionSolid("magnet_clamp", clamp_main, plate_combo, 0, trans);
  
  
  // ** Logical
  G4Material* clamp_material = G4Material::GetMaterial(this->magnet_clamp_material);
  mclamp_shield_log = new G4LogicalVolume(magnet_clamp, clamp_material, "magnet_clamp_shield_log", 0, 0, 0);
  mclamp_shield_log->SetVisAttributes(vis_att);
  
} // end::BuildMagnetClampPhotonShield()

void ApparatusSpiceTargetChamber::BuildElectroBox(){

	// ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(ELECTROBOX_COL));
  vis_att->SetVisibility(true);
  
  // ** Dimensions
  G4double outer_box_half_width = this->electrobox_outer_box_length/2.;
  G4double outer_box_half_length = this->electrobox_mid_length/2.;
  G4double inner_box_half_width = this->electrobox_inner_box_length/2.;
  G4double inner_box_half_length = this->electrobox_mid_length/2.;
  G4double inner_cyl_radius = this->electrobox_lip_inner_radius;
  
  
  // ** Shapes
  G4Box* outer_box = new G4Box("outer_box", outer_box_half_width, outer_box_half_width, outer_box_half_length);
  G4Box* inner_box = new G4Box("inner_box", inner_box_half_width, inner_box_half_width, inner_box_half_length);
  G4Tubs* inner_cyl = new G4Tubs("inner_cyl", 0, inner_cyl_radius, 2*inner_box_half_length, 0, 360*deg);
  G4UnionSolid* total_inner = new G4UnionSolid("total_inner", inner_box, inner_cyl);
  G4ThreeVector trans(0, 0, -this->electrobox_lip_length);
  G4SubtractionSolid* box_main = new G4SubtractionSolid("box_main", outer_box, total_inner, 0, trans);
  
  // ** Logical
  G4Material* electro_box_material = G4Material::GetMaterial(this->electro_box_material);
  electro_box_log = new G4LogicalVolume(box_main, electro_box_material, "electro_box_log", 0, 0, 0);
  electro_box_log->SetVisAttributes(vis_att);
  
} // end:BuildElectroBox()


void ApparatusSpiceTargetChamber::BuildColdFinger(){

	// ** Visualisation
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(CU_COL));
  vis_att->SetVisibility(true);
  
  // ** Dimensions
  G4double coldfinger_half_thickness = this->coldfinger_thickness/2.;
  G4double coldfinger_half_length = this->coldfinger_length/2.;
  G4double coldfinger_half_width = this->coldfinger_width/2.;
  G4double inner_hole_radius = this->coldfinger_hole_radius;

  // ** Shapes
  G4Box* plate_box = new G4Box("plate_box", coldfinger_half_length, coldfinger_half_width, coldfinger_half_thickness);
  G4Tubs* inner_hole = new G4Tubs("inner_hole", 0, inner_hole_radius, 2*coldfinger_half_thickness, 0, 360*deg);
  G4SubtractionSolid* plate_sub_hole = new G4SubtractionSolid("plate-hole", plate_box , inner_hole, 0, G4ThreeVector(0,0,0));
  
  // ** Logical
  G4Material* cold_finger_material = G4Material::GetMaterial(this->cold_finger_material);
  cold_finger_log = new G4LogicalVolume(plate_sub_hole, cold_finger_material, " cold_finger_log", 0, 0, 0);
  cold_finger_log->SetVisAttributes(vis_att);
  
} // end:BuildColdFinger()

void ApparatusSpiceTargetChamber::BuildS3CableHolder() {

  // Visualisation
  G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(CU_COL));
  sVisAtt->SetVisibility(true);
  
  // Dimensions
  
  
  
  // Shapes
  G4Box *sBox = new G4Box("sBox", 8.1*mm, 12.*mm, 40.5*mm);
  
  // Logical
  G4Material* sMaterial = G4Material::GetMaterial(this->s3_cable_case_material);
  fS3CaseLogical = new G4LogicalVolume(sBox, sMaterial, "s3_case_log", 0, 0, 0);
  fS3CaseLogical->SetVisAttributes(sVisAtt);



}




// **************************************************************************************************
// **************************************************************************************************
// *************************************PLACEMENT****************************************************
// **************************************************************************************************
// **************************************************************************************************

void ApparatusSpiceTargetChamber::PlaceTargetFrame(G4double d) {

 
  G4double sOffset = this->target_wheel_offset;
  G4double sOffsetAdj = sOffset - 15.00;

  G4double sPosZ = -7.75*mm; 
  G4double sPosY = (sOffset-sOffsetAdj)*sin(180*deg + d);
  G4double sPosX = sOffset + (sOffset-sOffsetAdj)*cos(180*deg + d);

  G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(-d);   
	
  
  G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);
  fFramePhysical = new G4PVPlacement(sRotate, sTranslate, fFrameLogical, "target_frame", expHallLog, false, 0);
  
}

void ApparatusSpiceTargetChamber::PlaceBeamPipe() {

  G4double sOffset = this->photon_shield_back_face_pos - (this->pipe_z_length / 2.) - (this->ps_det_clamp_thickness);
  G4ThreeVector move(0,0,sOffset);
  
  G4RotationMatrix* rotate = new G4RotationMatrix;
  
  beam_pipe_phys = new G4PVPlacement(rotate,move,beam_pipe_log, "beam_pipe",expHallLog,false,0);
  
}

void ApparatusSpiceTargetChamber::PlaceTargetWheelRods(G4double a, G4double r) {

 
  G4double sOffset = this->target_wheel_offset; 
  G4double sOffsetAdj = sOffset - r;

  G4double sPosZ = -3.75*mm; 
  G4double sPosY = (sOffset-sOffsetAdj)*sin(180*deg + a);
  G4double sPosX = sOffset + (sOffset-sOffsetAdj)*cos(180*deg + a);

  G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(-a);   
	//sRotate->rotateX(180.*deg); 
	//sRotate->rotateY(90.*deg);
  
  G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);
  fRodPhysical = new G4PVPlacement(sRotate, sTranslate, fRodLogical, "target_rods", expHallLog, false, 0);
  
}

void ApparatusSpiceTargetChamber::PlaceCollimator() {

  G4double sOffset = this->target_wheel_offset; 
  G4double sOffsetAdj = sOffset - 15.2;

  G4double sPosZ = -7.75*mm; 
  G4double sPosY = (sOffset-sOffsetAdj)*sin(-55*deg);
  G4double sPosX = sOffset + (sOffset-sOffsetAdj)*cos(-55*deg);

  G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(55*deg);   
  
  G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);
  fCollimPhysical = new G4PVPlacement(sRotate, sTranslate, fCollimLogical, "collimator", expHallLog, false, 0);

}

void ApparatusSpiceTargetChamber::PlaceTargetWheelExtension() {

  G4double sPosZ = 8.0*mm; 
  G4double sPosY = 0.*mm; 
  G4double sPosX = this->target_wheel_offset;

  G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);
  fExtPhysical = new G4PVPlacement(0, sTranslate, fExtLogical, "wheel_extension", expHallLog, false, 0);

}


void ApparatusSpiceTargetChamber::PlaceTargetChamberFrontRing()
{

  G4double z_position = this->target_chamber_sphere_centre
    + this->target_chamber_sphere_cut_off_point + this->front_ring_length/2. + this->frontDomeOffset;

  G4ThreeVector move(0,0,z_position);
  target_chamber_front_ring_phys = new G4PVPlacement(0,move,target_chamber_front_ring_log,
						     "target_chamber_front_ring",expHallLog,
						     false,0);
						     
  move.setZ(z_position +  this->front_ring_length/2. + this->front_ring_cylinder_length + this->cone_length/2.);
  target_chamber_front_cone_phys = new G4PVPlacement(0, move, target_chamber_front_cone_log,
						     "target_chamber_front_cone", expHallLog,
						     false,0);

}// end:PlaceTargetChamberFront()

void ApparatusSpiceTargetChamber::PlaceTargetChamberSphere()
{
  G4double z_position = this->target_chamber_sphere_centre + this->frontDomeOffset;
  G4ThreeVector move(0,0,z_position);
  target_chamber_sphere_phys = new G4PVPlacement(0,move,target_chamber_sphere_log,
						 "target_chamber_sphere",expHallLog,
						 false,0);

}// end:PlaceTargetChamberSphere()

void ApparatusSpiceTargetChamber::PlaceTargetChamberCylinderDownstream()
{

  G4double z_position = -(this->target_chamber_distance_from_target + this->target_chamber_cylinder_length/2. + this->middleRingOffset);
  G4ThreeVector move(0,0,z_position);
  target_chamber_cylinder_down_phys = new G4PVPlacement(0,move,target_chamber_cylinder_down_log,
						   "target_chamber_cylinder_downstream",expHallLog,
						   false,0);

}// end:PlaceTargetChamberCylinderDownstream()

void ApparatusSpiceTargetChamber::PlaceTargetWheel()
{
    G4RotationMatrix* rotate = new G4RotationMatrix;  
    rotate->rotateZ((0)*deg); 
    
  G4double offset = this->target_wheel_offset;
  G4double z_position =  this->target_wheel_thickness/2. + this->targetWheelOffset; 
  //G4double z_position = this->target_mount_plate_z_offset + this->target_mount_plate_thickness/2. + this->targetWheelOffset;
  
  G4ThreeVector move(offset, 0, z_position);
  target_wheel_phys = new G4PVPlacement(rotate, move, target_wheel_log,
					"target_wheel", expHallLog,
					false,0);
  
}// end:PlaceTargetWheel()

void ApparatusSpiceTargetChamber::PlaceTargetWheelGears()
{
    
  G4double z_offset = this->all_gear_z_offset + this->all_gear_thickness/2. + this->targetWheelOffset;
  G4double first_gear_offset = this->first_gear_plane_offset;
  G4double second_gear_offset = this->second_gear_plane_offset;
  G4double third_gear_offset = this->third_gear_plane_offset;
  
  G4ThreeVector move(-first_gear_offset, 0, z_offset);
  first_gear_phys = new G4PVPlacement(0, move, first_gear_log,
				      "first_gear", expHallLog,
				      false,0);
				      
  move.setX(-second_gear_offset);
  second_gear_phys = new G4PVPlacement(0, move, second_gear_log,
				       "second_gear", expHallLog,
				       false,0);
				       
  move.setX(third_gear_offset);
  third_gear_phys = new G4PVPlacement(0, move, third_gear_log,
				      "third_gear", expHallLog,
				      false,0);
  
} // end PlaceTargetWheelGears()

void ApparatusSpiceTargetChamber::PlaceTargetWheelGearPlates()
{
    
  G4double z_offset = this->all_gear_z_offset 
    + this->all_gear_thickness 
    + this->gear_plate_thickness/2.
    + this->targetWheelOffset;
  
  G4double first_gear_offset = this->first_gear_plane_offset;
  G4double second_gear_offset = this->second_gear_plane_offset;
  
  G4ThreeVector move(-first_gear_offset, 0, z_offset);
  gear_plate_one_phys = new G4PVPlacement(0, move, gear_plate_one_log,
					  "gear_plate_one", expHallLog,
					  false,0);
  move.setX(-second_gear_offset);
  gear_plate_one_phys = new G4PVPlacement(0, move, gear_plate_two_log,
					  "gear_plate_two", expHallLog,
					  false,0);	
  
}	// end::PlaceTargetWheelGearPlates()

void ApparatusSpiceTargetChamber::PlaceGearStick()
{
         
	G4double z_offset = -this->gear_stick_length/2. + this->gear_stick_z_offset + this->targetWheelOffset;
	G4double plane_offset = this->first_gear_plane_offset;
	
	G4ThreeVector move(-plane_offset, 0 , z_offset);
	gear_stick_phys = new G4PVPlacement(0, move, gear_stick_log,
					    "gear_stick", expHallLog,
					    false,0);

} // end PlaceGearStick();

void ApparatusSpiceTargetChamber::PlaceTargetMountPlate() 
{

    G4RotationMatrix* rotate = new G4RotationMatrix;  
    rotate->rotateZ(45*deg); 
    
	G4double z_offset = this->target_mount_plate_z_offset + this->target_mount_plate_thickness/2. + this->targetWheelOffset;
	G4ThreeVector move(0, 0, z_offset);
	target_mount_plate_phys = new G4PVPlacement(rotate, move, target_mount_plate_log,
						    "target_mount_plate", expHallLog,
						    false,0);

} // end PlaceTargetMountPlate()

void ApparatusSpiceTargetChamber::PlaceBiasPlate()
{

    G4RotationMatrix* rotate = new G4RotationMatrix;  
    rotate->rotateZ(45*deg); 
      
  G4double z_offset = this->bias_plate_offset_z - bias_plate_thickness/2. + this->targetWheelOffset;
  G4ThreeVector move(0,0,z_offset);
  bias_plate_phys = new G4PVPlacement(rotate, move, bias_plate_log,
				      "bias_plate", expHallLog, 
				      false,0);
  
} // end:PlaceBiasPlate()

void ApparatusSpiceTargetChamber::PlacePhotonShield()
{

  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  G4double z_position =  this->photon_shield_back_face_pos
    + this->photon_shield_layer_three_thickness/2.
    + this->photon_shield_layer_two_thickness/2. 
    + this->photon_shield_layer_one_thickness/2.
    + this->middleRingOffset;
  
  G4ThreeVector move(0,0,z_position);
  photon_shield_layer_one_phys = new G4PVPlacement(rotate,move,photon_shield_layer_one_log, "photon_shield_layer_one",expHallLog,
						   false,0);
  
  move.set(0,0, this->photon_shield_back_face_pos
	   + this->photon_shield_layer_three_thickness/2.
	   + this->photon_shield_layer_two_thickness/2.
	   + this->photon_shield_layer_one_thickness/2.);
    photon_shield_layer_two_phys = new G4PVPlacement(rotate,move,photon_shield_layer_two_log,
  						   "photon_shield_layer_two", expHallLog,
  						   false,0);
  
  move.set(0,0, this->photon_shield_back_face_pos
	   + this->photon_shield_layer_three_thickness/2. 
	   + this->photon_shield_layer_two_thickness/2. 
	   + this->photon_shield_layer_one_thickness/2.);
    photon_shield_layer_three_phys = new G4PVPlacement(rotate,move,photon_shield_layer_three_log,
  						     "photon_shield_layer_three", expHallLog,
  						     false,0);
  
}// end:PlacePhotonShield()

void ApparatusSpiceTargetChamber::PlaceShieldCovering()
{

  G4double z_position = this->photon_shield_back_face_pos + this->photon_shield_length/2. + this->middleRingOffset;

  G4ThreeVector move(0,0,z_position);

  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateX(0);
  rotate->rotateY(0);
  rotate->rotateZ(0*deg);
 
  shield_cover_phys = new G4PVPlacement(rotate,move,shield_cover_log,
					"shield_cover",expHallLog,
					false,0

);

}// end:PlaceShieldCovering()

void ApparatusSpiceTargetChamber::PlacePhotonShieldClamps() 
{
  
   G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(45*deg);
  
  G4double z_position = this->photon_shield_back_face_pos + this->middleRingOffset;
  
  G4ThreeVector move(0, 0, z_position - this->ps_det_clamp_thickness/2.);
  ps_detector_clamp_phys = new G4PVPlacement(rotate, move, ps_detector_clamp_log,
					     "photon_shield_detector_clamp",expHallLog,
					     false,0);
  
  move.setZ(z_position + this->photon_shield_length + this->ps_target_clamp_thickness/2.);
  ps_target_clamp_phys = new G4PVPlacement(rotate, move, ps_target_clamp_log,
					   "photon_shield_target_clamp",expHallLog,
					   false,0);
  
} // end::PlacePhotonShieldClamps()

void ApparatusSpiceTargetChamber::PlacePhotonShieldClampBolts(G4int copyID)
{
 
   G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(45*deg);
   
  G4double target_bolt_z_position = this->photon_shield_back_face_pos 
    + this->photon_shield_length + this->ps_target_clamp_thickness 
    - this->ps_target_bolt_length/2. + this->middleRingOffset;
  G4double distance_from_beam_line = this->ps_target_bolt_plane_offset;
  
  G4ThreeVector move = TranslateMagnets(copyID, distance_from_beam_line, target_bolt_z_position);
  target_bolt_phys = new G4PVPlacement(rotate, move, target_bolt_log,
				       "target_bolt", expHallLog, 
				       false,0);
  
  G4double detector_bolt_z_position = this->photon_shield_back_face_pos  + this->ps_det_bolt_length/2. + this->pipe_z_length/2. + //this->middleRingOffset - this->pipe_z_length - this->ps_det_clamp_thickness/2.;
  this->middleRingOffset - this->ps_det_clamp_thickness/2.;
  distance_from_beam_line = this->ps_det_bolt_plane_offset;
  G4ThreeVector move2 = TranslateMagnets(copyID, distance_from_beam_line, detector_bolt_z_position);

  detector_bolt_phys = new G4PVPlacement(rotate, move2, detector_bolt_log,
					 "detector_bolt", expHallLog,
					 false,0);
  
} // end::PlacePhotonShieldClampBoots()

void ApparatusSpiceTargetChamber::PlaceCollectorMagnet(G4int copyID)
{

  // ** Position Co-ordinates
  G4double magnetPosX = this->plate_one_edge_x + this->plate_one_length/2.;
  G4double magnetPosZ = -(this->distance_from_target + this->plate_one_height/2.) + this->middleRingOffset; // - (4mm + 50/2mm)  + 0 

  G4RotationMatrix* rotation;
  rotation = RotateMagnets(copyID);
  G4double radial_position = magnetPosX; 
 
  G4ThreeVector move = TranslateMagnets(copyID, radial_position, magnetPosZ);

  // ** Physical Volume
  magnet_phys = new G4PVPlacement(rotation, move, magnet_log,
				  "collector_magnet",expHallLog,
				  false,0);
				  
} // end:PlaceCollectorMagnet()

void ApparatusSpiceTargetChamber::PlaceMagnetCovering(G4int copyID)
{
  
  // ** Position Co-ordinates
  G4double magnetPosX = this->plate_one_edge_x + this->plate_one_length/2.;
  G4double magnetPosZ = -(this->distance_from_target + this->plate_one_height/2. ) + this->middleRingOffset; 

  G4RotationMatrix* rotation;
  rotation = RotateMagnets(copyID);
  G4double radial_position = magnetPosX; 
 
  G4ThreeVector move = TranslateMagnets(copyID, radial_position, magnetPosZ);

  // ** Physical Volume
  magnet_cover_phys = new G4PVPlacement(rotation, move, magnet_cover_log,
					"magnet_covering", expHallLog,
					false,0);

}

void ApparatusSpiceTargetChamber::PlaceMagnetClampChamber(G4int copyID)
{
  
  G4double z_position = -(this->distance_from_target + this->mclamp_chamber_height/2.) + this->middleRingOffset;
  G4double x_position = this->target_chamber_cylinder_inner_radius - this->mclamp_chamber_length/2.;
  
  G4RotationMatrix* rotation;
  rotation = RotateMagnets(copyID);
  
  G4ThreeVector move = TranslateMagnets(copyID, x_position, z_position);
  
  mclamp_chamber_phys = new G4PVPlacement(rotation, move, mclamp_chamber_log,
					  "magnet_clamp_chamber", expHallLog,
					  false,0);
  
} // end::PlaceMagnetClampChamber()

void ApparatusSpiceTargetChamber::PlaceMagnetClampPhotonShield(G4int copyID)
{
  
  G4double x_position = this->plate_one_edge_x 
    + (this->plate_one_length - this->plate_two_length
       + this->psclamp_psbuffer + this->psclamp_chamber_buffer)/2.
    - this->psclamp_psbuffer;
  G4double z_position = this->photon_shield_back_face_pos + this->photon_shield_length/2. + this->middleRingOffset;
  
  G4RotationMatrix* rotation;
  rotation = RotateMagnets(copyID);
  
  G4ThreeVector move = TranslateMagnets(copyID, x_position, z_position);
  
  mclamp_shield_phys = new G4PVPlacement(rotation, move, mclamp_shield_log,
					 "magnet_clamp_shield", expHallLog, 
					 false, 0);
  
} // end::PlaceMagnetClampPhotonShield()

void ApparatusSpiceTargetChamber::PlaceElectroBox()
{
  
  G4double z_offset = electrobox_z_offset - electrobox_mid_length/2. + this->backDetAndAlBoxOffset;
  
  G4ThreeVector move(0, 0, z_offset);
  
  electro_box_phys = new G4PVPlacement(0, move, electro_box_log,
				       "electro_box", expHallLog, 
				       false, 0);
  
} // end:PlaceElectroBox()

void ApparatusSpiceTargetChamber::PlaceColdFinger()
{
  
  G4double z_offset = coldfinger_z_offset  - coldfinger_thickness/2.  ;
  
  G4ThreeVector move(0, 0, z_offset);
  
  cold_finger_phys = new G4PVPlacement(0, move, cold_finger_log,
				       "cold_finger", expHallLog, 
				       false, 0);
  
} // end:PlaceColdFinger()

void ApparatusSpiceTargetChamber::PlaceS3CableHolder() {

  G4double sPosZ = -(40.5 + 6.)*mm;   
  G4double sPosY = 87.1*sin(-22.5*deg);
  G4double sPosX = -87.1*cos(-22.5*deg);

  G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);
  
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  sRotate->rotateZ(-22.5*deg);
  
  fS3CasePhysical = new G4PVPlacement(sRotate, sTranslate, fS3CaseLogical, "s3_cable_case", expHallLog, false, 0);

}


////////////////////////////////////////////////////////////////////
//Translations and Rotations of Magnets (same for all parts)      //
////////////////////////////////////////////////////////////////////

G4RotationMatrix* ApparatusSpiceTargetChamber::RotateMagnets(G4int copyID)
{

  G4RotationMatrix* rotate = new G4RotationMatrix;
  if(this->NUMBER_OF_MAGNETS == 8) 
  	rotate->rotateZ(-(copyID+0.5)*45.*deg);
  else 
  	rotate->rotateZ(-(copyID+0.5)*90.*deg);

  return rotate;
}

G4ThreeVector ApparatusSpiceTargetChamber::TranslateMagnets(G4int copyID, G4double radial_position, G4double z_position)
{
  G4double x_position(0);
  G4double y_position(0);
  if(this->NUMBER_OF_MAGNETS == 8){
    x_position = radial_position*cos((copyID+0.5)*45.*deg);
    y_position = radial_position*sin((copyID+0.5)*45.*deg);
  }
  else{
    x_position = radial_position*cos((copyID+0.5)*90.*deg);
    y_position = radial_position*sin((copyID+0.5)*90.*deg);
  }
  return G4ThreeVector(x_position, y_position, z_position);
}


