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

#ifndef DetectionSystemGriffin_h
#define DetectionSystemGriffin_h 1

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectionSystemGriffin
{
public:
    DetectionSystemGriffin(G4int sel, G4int suppSwitch, G4double detRad, G4int hevimetSel );
    ~DetectionSystemGriffin();

    void Build() ;
    // For detector specific dead layers
    void BuildDeadLayerSpecificCrystal(G4int det);
    void BuildEverythingButCrystals();
    G4double GetCrystalDistanceFromOrigin() {return crystal_dist_from_origin;}

    G4double transX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi);

    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector moveBAH, G4RotationMatrix* rotateBAH, G4int detector_number);
    // For detector specific dead layers
    G4int PlaceDeadLayerSpecificCrystal(G4LogicalVolume* exp_hall_log, G4int detector_number, G4int position_number, G4bool posTigress);
    G4int PlaceEverythingButCrystals(G4LogicalVolume* exp_hall_log, G4int detector_number, G4int position_number, G4bool posTigress);

private:
    G4String sdName0;
    G4String sdName1;
    G4String sdName2;
    G4String sdName3;
    G4String sdName4;
    G4String sdName5;
    G4String colNameGe;
    G4String colNameLeftCasing;
    G4String colNameRightCasing;
    G4String colNameLeftExtension;
    G4String colNameRightExtension;
    G4String colNameBackPlug;

    G4bool include_extension_suppressors;
    G4bool include_side_suppressors;
    G4bool include_back_suppressors;

    G4String back_suppressor_material;
    G4String BGO_material;

    G4RotationMatrix* rotate_null;
    G4ThreeVector move_null;

    G4double cut_clearance;
    G4double extra_cut_length;

    double coords[20][5];
    
    //Jun 21, 2005, Epapr7.80: modification to the prototype suppressor shield,
    //since the hevimet collimators they provided were too big, and therefore the
    //extensions do not come far enough forward
    G4double extension_accidental_back_shift;

    G4String hevimet_choice;

    //Settings for the suppressor shield
    G4double thickness_double;
    G4double side_thickness_double;
    G4double extension_thickness_double;

    G4double radial_distance;
    G4int number_of_detectors;
    G4double origin_to_crystal_distance;
    G4bool dead_layer_include_flag;
    G4double inner_dead_layer_thickness;
    G4double outer_dead_layer_thickness;

    //Cold Finger
    G4double cold_finger_outer_shell_radius;
    G4double cold_finger_shell_thickness;
    G4double cold_finger_shell_length;
    G4double rear_plate_thickness;
    G4double cold_finger_end_plate_thickness;
    G4double cold_finger_radius;
    G4double cold_finger_space;
    G4double cold_finger_length;
    G4double coolant_length;
    G4double coolant_radius;
    
    //New cooling structures (Jan 2005)
    G4double extra_block_thickness;
    G4double extra_block_distance_from_back_plate;
    G4double extra_block_inner_diameter;
    G4double triangle_posts_distance_from_crystals;
    G4double triangle_post_starting_depth;
    G4double fet_air_hole_radius;
    G4double cooling_side_block_thickness;
    G4double cooling_side_block_width;
    G4double cooling_side_block_horizontal_depth;
    G4double structureMat_cold_finger_thickness;
    G4double structureMat_cold_finger_radius;
    G4double cooling_bar_thickness;
    G4double cooling_bar_width;

    G4int suppressor_design;
    G4int suppressor_position_selector;
    G4int hevimet_selector;

    G4double forward_inner_radius;
    G4double back_inner_radius;

    G4double suppressor_forward_radius;
    G4double suppressor_back_radius;


    // For the optimization of the depth segmentation
    G4double depth_segmentation_adjustment;

    //the germanium detector's variables (for one "leaf" of the clover)
    G4double germanium_outer_radius;
    G4double germanium_hole_radius;
    G4double germanium_width;
    G4double germanium_length;
    G4double germanium_hole_dist_from_face;
    G4double germanium_dist_from_can_face;
    G4double germanium_bent_length;
    G4double germanium_shift;    		//the amount by which the two sides adjacent to the other
    //germanium crystals are cut closer to the center

    G4double germanium_separation;  		//the space between the quarter detectors
    G4double inter_crystal_electrodeMat_thickness;
    G4double electrodeMat_starting_depth;
    G4double germanium_corner_cone_end_length; 	//the ending length of the cones
    //at the corners of each leaf
    
    //basic germanium detector's values
    G4double detector_block_length;
    G4double detector_block_height;
    G4double detector_block_trapezoidal_height;

    G4double detector_total_width;
    G4double detector_total_length;
    G4double bent_end_angle;
    G4double bent_end_length;
    G4double can_face_thickness;
    G4double can_side_thickness;

    G4double hevimet_tip_thickness;
    G4double hevimet_tip_angle;

    //Values for the BGO
    G4double suppressor_cut_extension;
    G4double suppressor_shell_thickness;
    
    G4double back_BGO_thickness;
    G4double BGO_chopped_tip;
    G4double BGO_movable_space;
    G4double side_suppressor_back_offset;
    G4double side_BGO_thickness;
    G4double BGO_can_seperation;
    G4double side_BGO_length;
    G4double side_suppressor_length;
    G4double BGO_trap_length;
    G4double suppressor_extension_thickness;
    G4double suppressor_extension_length;
    G4double suppressor_extension_angle;
    
    G4double suppressor_extension_length_det ;

    //Values for the HeavyMet
    G4double HeavyMet_thickness;
    G4double HeavyMet_inside_angle;

    G4double air_box_front_width;
    G4double air_box_front_length;
    G4double air_box_back_length;
    
    G4double air_box_back_length_det ;
    G4double air_box_front_length_det ;
    G4double air_box_front_width_det ;

    G4double shift;
    G4double suppShift ;

    G4int copy_number;
    G4int copy_number_two;
    G4int germanium_copy_number;
    G4int left_suppressor_side_copy_number;
    G4int right_suppressor_side_copy_number;
    G4int left_suppressor_extension_copy_number;
    G4int right_suppressor_extension_copy_number;
    G4int back_suppressor_copy_number;
    
    G4double rhombi_diameter;
    G4double new_rhombi_radius;

    G4double new_rhombi_radius_det;

    G4double detector_position_shift;
    G4double applied_back_shift;

    G4int germanium_selector;
    G4int can_selector;
    G4int BGO_selector;
    G4int cold_finger_selector;

    //Redacted parameters/////////////
    G4double detectorPlacementCxn;
    G4double trianglePostDim;
    G4double suppressorExtRightX;
    G4double suppressorExtRightY;
    G4double suppressorExtRightZ;
    G4double suppressorExtLeftX;
    G4double suppressorExtLeftY;
    G4double suppressorExtLeftZ;
    G4double wedgeDim;
    G4double quarterDetectorCxn;
    G4double quarterDetectorCxnB;
    G4String electrodeMaterial;
    G4String structureMaterial;

    // Assembly volumes
    G4AssemblyVolume* assembly;
    G4AssemblyVolume* germaniumAssembly;
    G4AssemblyVolume* leftSuppressorCasingAssembly;
    G4AssemblyVolume* rightSuppressorCasingAssembly;
    G4AssemblyVolume* leftSuppressorExtensionAssembly;
    G4AssemblyVolume* rightSuppressorExtensionAssembly;
    G4AssemblyVolume* suppressorBackAssembly;
    // For detector specific dead layers
    G4AssemblyVolume* assemblyCry[4];
    G4AssemblyVolume* germaniumAssemblyCry[4];
    G4AssemblyVolume* leftSuppressorCasingAssemblyCry[4];
    G4AssemblyVolume* rightSuppressorCasingAssemblyCry[4];
    G4AssemblyVolume* leftSuppressorExtensionAssemblyCry[4];
    G4AssemblyVolume* rightSuppressorExtensionAssemblyCry[4];
    G4AssemblyVolume* suppressorBackAssemblyCry[4];
    G4AssemblyVolume* suppressorShellAssembly;
    G4AssemblyVolume* backAndSideSuppressorShellAssembly ;
    G4AssemblyVolume* extensionSuppressorShellAssembly ;
    G4AssemblyVolume* hevimetAssembly ;
    

    // Logical volumes
    G4LogicalVolume* air_box_log;

    // methods to construct all of the components of the detector
    void ConstructNewSuppressorCasingWithShells();
    void BuildelectrodeMatElectrodes();
    void ConstructComplexDetectorBlockWithDeadLayer();

    // For detector specific dead layers
    void ConstructComplexDetectorBlockWithDetectorSpecificDeadLayer(G4int det, G4int cry);
    void ConstructNewSuppressorCasingDetectorSpecificDeadLayer(G4int det, G4int cry);
    void ConstructDetector();
    void ConstructComplexDetectorBlock();
    void ConstructColdFinger();
    void ConstructNewHeavyMet();

    //LogicalVolumes used in ConstructBasicDetectorBlock
    G4LogicalVolume* germanium_block_log;

    //Logical Volumes used in ConstructComplexDetectorBlock:
    G4LogicalVolume* germanium_block1_log;
    G4LogicalVolume* germanium_hole_log;
    G4LogicalVolume* inner_dead_layer_log;
    G4LogicalVolume* inner_dead_layer_cap_log;
    G4LogicalVolume* outer_dead_layer_log;
    
    G4LogicalVolume* inter_crystal_electrodeMat_back_log;
    G4LogicalVolume* inter_crystal_electrodeMat_front_log;
    
    //Logical Volumes used in ConstructBGOCasing:
    G4LogicalVolume* back_BGO_log;
    G4LogicalVolume* BGO_casing_log;

    //Logical Volumes used in ConstructNewSuppressorCasing:
    G4LogicalVolume* back_quarter_suppressor_shell_log;
    G4LogicalVolume* right_suppressor_shell_log;
    G4LogicalVolume* left_suppressor_shell_log;
    G4LogicalVolume* right_suppressor_shell_extension_log;
    G4LogicalVolume* left_suppressor_shell_extension_log;

    G4LogicalVolume* cap_for_right_suppressor_log;

    G4LogicalVolume* back_quarter_suppressor_log;
    G4LogicalVolume* right_suppressor_log;
    G4LogicalVolume* left_suppressor_log;
    G4LogicalVolume* right_suppressor_extension_log;
    G4LogicalVolume* left_suppressor_extension_log;
    
    //Logical Volumes used in ConstructDetector:
    G4LogicalVolume* front_face_log;
    G4LogicalVolume* right_bent_piece_log;
    G4LogicalVolume* left_bent_piece_log;
    G4LogicalVolume* top_bent_piece_log;
    G4LogicalVolume* bottom_bent_piece_log;
    G4LogicalVolume* right_wedge_log;
    G4LogicalVolume* left_wedge_log;
    G4LogicalVolume* top_wedge_log;
    G4LogicalVolume* bottom_wedge_log;
    G4LogicalVolume* upper_right_cone_log;
    G4LogicalVolume* lower_right_cone_log;
    G4LogicalVolume* upper_left_cone_log;
    G4LogicalVolume* lower_left_cone_log;
    G4LogicalVolume* upper_right_tube_log;
    G4LogicalVolume* lower_right_tube_log;
    G4LogicalVolume* upper_left_tube_log;
    G4LogicalVolume* lower_left_tube_log;
    G4LogicalVolume* right_side_panel_log;
    G4LogicalVolume* left_side_panel_log;
    G4LogicalVolume* top_side_panel_log;
    G4LogicalVolume* bottom_side_panel_log;
    G4LogicalVolume* rear_plate_log;
    G4LogicalVolume* finger_shell_log;
    G4LogicalVolume* tank_log;

    //Logical Volumes used in ConstructColdFinger:
    G4LogicalVolume* end_plate_log;
    G4LogicalVolume* finger_log;
    G4LogicalVolume* extra_cold_block_log;
    G4LogicalVolume* triangle_post_log;
    G4LogicalVolume* fet_air_hole_log;
    G4LogicalVolume* cooling_bar_log;
    G4LogicalVolume* cooling_side_block_log;
    G4LogicalVolume* structureMat_cold_finger_log;
    
    //Logical Volumes used in ConstructNewHeavyMet:
    G4LogicalVolume* hevimet_log;

    //internal methods for ConstructCan()
    G4Box* squareFrontFace();
    G4Trap* cornerWedge();
    G4Para* bentSidePiece();
    G4Box* otherBentSidePiece();
    G4Cons* roundedEndEdge();
    G4Tubs* cornerTube();
    G4Box* sidePanel();
    G4SubtractionSolid* rearPlate();
    G4Tubs* coldFingerShell();
    G4Tubs* liquidNitrogenTank();

    //internal methods for ConstructBasicDetectorBlock()
    G4Box* rectangularSegment();
    G4Trd* trapezoidalSegment();

    //internal methods for ConstructComplexDetectorBlock()
    G4SubtractionSolid* quarterDetector();
    // For detector specific dead layers
    G4SubtractionSolid* quarterSpecificDeadLayerDetector(G4int det, G4int cry);

    //internal methods for ConstructComplexDetectorBlockWithPlastic()
    G4UnionSolid* interCrystalelectrodeMatBack();
    G4UnionSolid* interCrystalelectrodeMatFront();
    
    //internal methods for ConstructColdFinger()
    G4Tubs* airHole();
    G4Tubs* airHoleCut();
    G4Box* endPlate();
    G4Tubs* finger();
    G4SubtractionSolid* extraColdBlock();
    G4Trd* trianglePost();
    G4Box* coolingBar();
    G4Box* coolingSideBlock();
    G4Tubs* structureMatColdFinger();
    
    //internal methods for ConstructBGOCasing()
    G4SubtractionSolid* backBGO();
    G4SubtractionSolid* BGOCasing();
    G4SubtractionSolid* frontBGO();
    G4Trd* sideBGO();

    //internal methods for ConstructNewSuppressorCasing()
    G4SubtractionSolid* backSuppressorQuarter();
    G4SubtractionSolid* frontSlantSuppressor(G4String sidePosition, G4bool choppingSuppressor) ;
    G4SubtractionSolid* sideSuppressorExtension(G4String sidePosition, G4bool choppingSuppressor) ;

    //internal methods for New SuppressorCasingWithShells
    G4SubtractionSolid* shellForBackSuppressorQuarter();
    G4SubtractionSolid* shellForFrontSlantSuppressor(G4String sidePosition) ;
    G4SubtractionSolid* shellForSuppressorExtension(G4String sidePosition);

    //internal methods for ConstructNewHeavyMet()
    G4SubtractionSolid* newHeavyMet();

    G4String crystal_material;
    G4String can_material;
    G4String vacuum_material;
    G4double crystal_length_x;
    G4double crystal_length_y;
    G4double crystal_length_z;
    G4double crystal_inner_radius;
    G4double crystal_outer_radius;
    G4double can_thickness;
    G4double can_inner_radius;
    G4double can_lid_inner_radius;
    G4double can_lid_outer_radius;
    G4double can_front_lid_thickness;
    G4double can_back_lid_thickness;
    G4double can_face_dist_from_origin;
    G4double crystal_dist_from_can_face;
    G4double crystal_dist_from_can_back;
    G4double can_length_z;
    G4double crystal_dist_from_origin;

    // For detector specific dead layers
    G4double griffinDeadLayers[16][4];
    G4Colour griffinCrystalColours[4];
    G4Colour griffinDeadLayerColours[4];

    // internal methods
    void BuildOneDetector();
    //    void PlaceDetector(G4int detector_number);


};

#endif

