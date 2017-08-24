//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ApparatusSpiceTargetChamber_h
#define ApparatusSpiceTargetChamber_h 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;


#define SN_COL 1.0, 1.0, 1.0
#define CU_COL 1.0, 1.0, 0.0
//#define PEEK_COL 0.0, 1.0, 0.0 //was commented out - error as couldn't find
#define PEEK_COL 0.5, 0.5, 0.0
#define KAPTON_COL 0.2, 0.7, 0.1
#define PB_COL 0.0, 0.0, 0.0
#define NDFEB_COL 0.7,0.3,0.3
#define DELRIN_COL 0.0, 0.0, 1.0
#define AL_COL 0.5, 0.5, 0.5

// define custom colours for visualisation
#define TARGET_CHAMBER_COL 0.15,0.15,0.15
#define ELECTROBOX_COL 0.15,0.15,0.15


class ApparatusSpiceTargetChamber
{
public:
  ApparatusSpiceTargetChamber(G4String, G4double);
  ~ApparatusSpiceTargetChamber();

public:
  void Build(G4LogicalVolume*);
  
private:
  G4LogicalVolume* expHallLog;

private:
  ////////////////////////////////////////////
  // Logical Volumes used in ApparatusSpiceTargetChamber
  ////////////////////////////////////////////
  G4LogicalVolume* target_chamber_front_ring_log;
  G4LogicalVolume* target_chamber_front_cone_log;
  G4LogicalVolume* target_chamber_sphere_log;
  G4LogicalVolume* target_chamber_cylinder_down_log;
  G4LogicalVolume* target_wheel_log;
  G4LogicalVolume* fFrameLogical;
  G4LogicalVolume* fCollimLogical;
  G4LogicalVolume* fExtLogical;
  G4LogicalVolume* fRodLogical;
  G4LogicalVolume* first_gear_log;
  G4LogicalVolume* second_gear_log;
  G4LogicalVolume* third_gear_log;
  G4LogicalVolume* gear_plate_one_log;
  G4LogicalVolume* gear_plate_two_log;
  G4LogicalVolume* gear_stick_log;
  G4LogicalVolume* target_mount_plate_log;
  G4LogicalVolume* bias_plate_log;
  G4LogicalVolume* photon_shield_layer_one_log;
  G4LogicalVolume* photon_shield_layer_two_log;
  G4LogicalVolume* photon_shield_layer_three_log;
  G4LogicalVolume* ps_target_clamp_log;
  G4LogicalVolume* ps_detector_clamp_log;
  G4LogicalVolume* target_bolt_log;
  G4LogicalVolume* detector_bolt_log;
  G4LogicalVolume* magnet_log;
  G4LogicalVolume* magnet_cover_log;
  G4LogicalVolume* mclamp_chamber_log;
  G4LogicalVolume* mclamp_shield_log;
  G4LogicalVolume* electro_box_log;
  G4LogicalVolume* shield_cover_log;
  G4LogicalVolume* cold_finger_log;
  G4LogicalVolume* fS3CaseLogical;
  G4LogicalVolume* beam_pipe_log;
  G4LogicalVolume* conical_collimator_log;//11/8
  G4LogicalVolume* xray_insert_log;
  
  
private:
  ////////////////////////////////////////////
  // Physical Volumes used in ApparatusSpiceTargetChamber
  ////////////////////////////////////////////
  G4VPhysicalVolume* target_chamber_front_ring_phys;
  G4VPhysicalVolume* target_chamber_front_cone_phys;
  G4VPhysicalVolume* target_chamber_sphere_phys;
  G4VPhysicalVolume* target_chamber_cylinder_down_phys;
  G4VPhysicalVolume* target_wheel_phys;
  G4VPhysicalVolume* fFramePhysical;
  G4VPhysicalVolume* fCollimPhysical;
  G4VPhysicalVolume* fExtPhysical;
  G4VPhysicalVolume* fRodPhysical;
  G4VPhysicalVolume* first_gear_phys;
  G4VPhysicalVolume* second_gear_phys;
  G4VPhysicalVolume* third_gear_phys;
  G4VPhysicalVolume* gear_plate_one_phys;
  G4VPhysicalVolume* gear_plate_two_phys;
  G4VPhysicalVolume* gear_stick_phys;
  G4VPhysicalVolume* target_mount_plate_phys;
  G4VPhysicalVolume* bias_plate_phys;
  G4VPhysicalVolume* photon_shield_layer_one_phys;
  G4VPhysicalVolume* photon_shield_layer_two_phys;
  G4VPhysicalVolume* photon_shield_layer_three_phys;
  G4VPhysicalVolume* ps_target_clamp_phys;
  G4VPhysicalVolume* ps_detector_clamp_phys;
  G4VPhysicalVolume* target_bolt_phys;
  G4VPhysicalVolume* detector_bolt_phys;
  G4VPhysicalVolume* magnet_phys;
  G4VPhysicalVolume* magnet_cover_phys;
  G4VPhysicalVolume* mclamp_chamber_phys;
  G4VPhysicalVolume* mclamp_shield_phys;
  G4VPhysicalVolume* electro_box_phys;
  G4VPhysicalVolume* shield_cover_phys;
  G4VPhysicalVolume* cold_finger_phys;
  G4VPhysicalVolume* fS3CasePhysical;
  G4VPhysicalVolume* beam_pipe_phys;
  G4VPhysicalVolume* conical_collimator_phys;//11/8
  G4VPhysicalVolume* xray_insert_phys;
  
private:
  ////////////////////////////////////////////
  // Properties used in ApparatusSpiceTargetChamber
  ////////////////////////////////////////////

  //-------------------------
  // Magnet properties:
  //-------------------------
  G4int NUMBER_OF_MAGNETS;
  G4int fNumberOfFrames;
  G4double zOffset;

  //-------------------------
  // Materials:
  //-------------------------
  G4String magnet_material; //
  G4String target_chamber_material;//
  G4String downstream_cone_material;//
  G4String mount_lead; //doesn't exist
  G4String photon_shield_layer_one_material;//
  G4String photon_shield_layer_two_material;//
  G4String photon_shield_layer_three_material;//
  G4String ps_clamp_material;//
  G4String magnet_clamp_material;//
  G4String target_wheel_material;//
  G4String gear_stick_material;//
  G4String bias_plate_material;//
  G4String target_wheel_gear_material_a;//
  G4String target_wheel_gear_material_b;//
  G4String gear_plate_material;//
  G4String target_mount_plate_material;//
  G4String electro_box_material;//
  G4String small_bolt_material;//
  G4String large_bolt_material;//
  G4String shield_cover_material;//
  G4String magnet_cover_material;//
  G4String cold_finger_material;//
  G4String s3_cable_case_material;//
  G4String beam_pipe_material;//
  G4String conical_collimator_material;//11/8
  G4String xray_insert_material;//
  
  //-------------------------
  // Dimensions:
  //-------------------------

  //-----------------------------
  // Dimensions of Target Chamber
  //-----------------------------
  G4double target_chamber_cylinder_inner_radius;
  G4double target_chamber_cylinder_outer_radius;
  G4double target_chamber_cylinder_length;
  G4double target_chamber_distance_from_target;
  G4double target_chamber_cylinder_lip_inner_radius;
  G4double target_chamber_cylinder_lip_thickness;
  G4double target_chamber_cylinder_det_outer_radius;
  G4double target_chamber_cylinder_det_thickness;
  G4double target_chamber_cylinder_hole_length;
  G4double target_chamber_cylinder_hole_breadth;
  G4double target_chamber_cylinder_hole_distance_x;
  G4double target_chamber_cylinder_hole_distance_y;
  
  G4double target_chamber_sphere_inner_radius;
  G4double target_chamber_sphere_outer_radius;
  G4double target_chamber_sphere_centre;
  G4double target_chamber_sphere_cutting_outer_radius;
  G4double target_chamber_sphere_cut_off_point;
  G4double front_ring_inner_radius;
  G4double front_ring_outer_radius;
  G4double front_ring_length;
  G4double front_ring_cylinder_outer_radius;
  G4double front_ring_cylinder_length;
  
  G4double cone_front_inner_rad;
  G4double cone_front_outer_rad;
  G4double cone_length;

  
  // --------------------------
  // Dimensions of Target Wheel
  // --------------------------
  G4double target_wheel_radius;
  G4double target_wheel_thickness;
  G4double target_wheel_offset;
  G4double target_radius;
  G4double target_offset;
  G4double collimator_radius;
  // Bias Plate
  G4double bias_plate_outer_radius;
  G4double bias_plate_thickness;
  G4double bias_plate_offset_z;
  G4double bias_plate_cut;
  // Gears
  G4double first_gear_radius;
  G4double first_gear_plane_offset;
  G4double second_gear_radius;
  G4double second_gear_plane_offset;
  G4double third_gear_radius;
  G4double third_gear_plane_offset;
  G4double all_gear_thickness;
  G4double all_gear_z_offset;
  // Gear Plates
  G4double gear_plate_one_radius;
  G4double gear_plate_two_radius;
  G4double gear_plate_thickness;
  // Mount Plate
  G4double target_mount_plate_radius;
  G4double target_mount_plate_thickness;
  G4double target_mount_plate_z_offset;
  // Gear Stick
  G4double gear_stick_length;
  G4double gear_stick_radius;
  G4double gear_stick_z_offset;
  // Target Frame
  G4double fDimTargetFrameOutZ;
	G4double fDimTargetFrameOutY;
	G4double fDimTargetFrameCutY;

  //----------------------------
  // Dimensions of Photon Shield
  //----------------------------
  G4double photon_shield_front_radius;
  G4double photon_shield_back_radius;
  G4double photon_shield_inner_radius;
  G4double photon_shield_length;
  G4double photon_shield_back_face_pos;

  G4double photon_shield_layer_one_thickness;
  G4double photon_shield_layer_two_thickness;
  G4double photon_shield_layer_three_thickness;
  
  // ----------------------------------
  // Dimensions of Photon Shield Clamps
  // ----------------------------------
  G4double ps_target_clamp_radius;
  G4double ps_target_clamp_thickness;
  G4double ps_target_clamp_extension;
  G4double ps_det_clamp_radius;
  G4double ps_det_clamp_thickness;
  
  // ---------------------------------
  // Dimensions of Photon Shield Bolts
  // ---------------------------------
  G4double ps_target_bolt_length;
  G4double ps_target_bolt_radius;
  G4double ps_target_bolt_plane_offset;
  G4double ps_target_bolt_head_radius;
  G4double ps_target_bolt_head_length;
  G4double ps_det_bolt_length;
  G4double ps_det_bolt_radius;
  G4double ps_det_bolt_plane_offset;
  G4double ps_det_bolt_head_radius;
  G4double ps_det_bolt_head_length;
  
  //---------------------
  // Dimensions of Magnet
  //---------------------
  G4double plate_one_thickness;
  G4double plate_one_length;
  G4double plate_one_height;
  G4double plate_one_lower_height;
  
  G4double plate_one_edge_x;
  G4double plate_two_thickness;
  G4double plate_two_length;
  G4double plate_two_height;
  
  G4double distance_from_target;
  G4double cutting_box_angle;
  
  // ------------------------------------
  // Dimensions of Magnet Clamp (Chamber)
  // ------------------------------------
  G4double mclamp_chamber_thickness;
  G4double mclamp_chamber_length;
  G4double mclamp_chamber_height;
  
  // ------------------------------------------
  // Dimensions of Magnet Clamp (Photon Shield)
  // ------------------------------------------
  G4double psclamp_psbuffer;
  G4double psclamp_chamber_buffer;
  G4double psclamp_thickness_buffer;
  G4double psclamp_length;
	
	// -------------------------
  // Dimensions of ElectroBox
  // -------------------------
  G4double electrobox_outer_box_length;
  G4double electrobox_inner_box_length;
  G4double electrobox_mid_length;
  G4double electrobox_lip_length;
  G4double electrobox_lip_inner_radius;
  G4double electrobox_z_offset;

  // -------------------------
  // Dimensions of ColdFinger
  // -------------------------
  G4double coldfinger_thickness;
  G4double coldfinger_length;
  G4double coldfinger_width;
  G4double coldfinger_z_offset;
  G4double coldfinger_hole_radius;
    
  // ------------------------------
  // Dimensions of Plastic Coatings
  // ------------------------------
  G4double magnet_coating_thickness;
  G4double shield_coating_thickness;
  
  // ---------------------------
  // individual offsets for visualisation
  // ---------------------------
  G4double frontDomeOffset;
  G4double targetWheelOffset;
  G4double middleRingOffset;
  G4double backDetAndAlBoxOffset;
  
  // -----------------------
  // Dimensions of Beam Pipe
  // -----------------------
  G4double pipe_inner_radius;
  G4double pipe_outer_radius;
  G4double pipe_z_length;
  G4double pipe_z_offset;
  
  // --------------------------------
  // Dimensions of Conical Collimator
  // --------------------------------
  G4double mid_inner_radius;
  G4double mid_outer_radius;
  G4double mid_z_length;
  G4double outercyl_outer_radius;
  G4double outercyl_inner_radius;
  G4double outercyl_z_length;
  G4double innercyl_outer_radius;
  G4double innercyl_inner_radius;
  G4double innercyl_z_length;
  G4double edge_inner_radius;
  G4double edge_outer_radius;
  G4double edge_z_length;
  
  // --------------------------
  // Dimensions of X-ray Insert
  // --------------------------
  G4double insert_cyl_outer_radius;
  G4double insert_cyl_inner_radius;
  G4double insert_cyl_z_length;
  G4double insert_edge_inner_radius;
  G4double insert_edge_outer_radius;
  G4double insert_edge_z_length;
  G4double insert_hole_radius;
  G4double insert_hole_length;
  
  //-----------------------------
  // copy numbers
  //-----------------------------
  G4int targetChamberDownstreamCopyNumber;
  G4int photonShieldCopyNumber;
  G4int photonShieldClampCopyNumber;
  G4int magnetsCopyNumber;
  G4int magnetClampCopyNumber;
  G4int targetWheelCopyNumber;
  G4int biasPlateCopyNumber;
  G4int gearCopyNumber;
  G4int gearStickCopyNumber;
  G4int electroBoxCopyNumber;
  G4int photonShieldClampBoltCopyNumber;
  G4int shieldCoveringCopyNumber;
  G4int magnetCoveringCopyNumber;
  G4int coldFingerCopyNumber;
  
private:
  //////////////////////////////////////////////////////
  // internal methods and functions in ApparatusSpiceTargetChamber::Build()
  //////////////////////////////////////////////////////

  // methods
  void BuildTargetChamberFrontRing();
  void BuildTargetChamberSphere();
  void BuildTargetChamberCylinderDownstream();
  void BuildTargetWheel();
  void BuildTargetWheelSimple();
  void BuildTargetWheelRods();
  void BuildTargetWheelGears();
  void BuildTargetWheelGearPlates();
  void BuildTargetFrame();
  void BuildTargetWheelExtension();
  void BuildCollimator();
  void BuildGearStick();
  void BuildBiasPlate();
  void BuildTargetMountPlate();
  void BuildPhotonShield();
  void BuildPhotonShieldClamps();
  void BuildPhotonShieldClampBolts();
  void BuildCollectorMagnet();
  void BuildMagnetCovering();
  void BuildMagnetClampChamber();
  void BuildMagnetClampPhotonShield();
  void BuildElectroBox();
  void BuildShieldCovering();
  void BuildColdFinger();
  void BuildS3CableHolder();
  void BuildBeamPipe();
  void BuildConicalCollimator();
  void BuildXrayInsert();
  
  void PlaceTargetChamberFrontRing();
  void PlaceTargetChamberSphere();
  void PlaceTargetChamberCylinderDownstream();
  void PlaceTargetWheel();
  void PlaceTargetWheelRods(G4double, G4double);
  void PlaceTargetWheelGears();
  void PlaceTargetWheelGearPlates();
  void PlaceTargetFrame(G4double);
  void PlaceCollimator();
  void PlaceGearStick();
  void PlaceTargetMountPlate();
  void PlaceTargetWheelExtension();
  void PlaceBiasPlate();
  void PlacePhotonShield();
  void PlacePhotonShieldClamps();
  void PlacePhotonShieldClampBolts(G4int);
  void PlaceCollectorMagnet(G4int);
  void PlaceMagnetCovering(G4int);
  void PlaceMagnetClampChamber(G4int);
  void PlaceMagnetClampPhotonShield(G4int);
  void PlaceElectroBox();
  void PlaceShieldCovering();
  void PlaceColdFinger();
  void PlaceS3CableHolder();
  void PlaceBeamPipe();
  void PlaceConicalCollimator();
  void PlaceXrayInsert();
  
  G4double targetz;
  // functions
  G4RotationMatrix* RotateMagnets(G4int);
  G4ThreeVector TranslateMagnets(G4int,G4double,G4double);
};
#endif
