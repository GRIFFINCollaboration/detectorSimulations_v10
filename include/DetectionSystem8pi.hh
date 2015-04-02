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
// $Id: DetectionSystem8pi.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectionSystem8pi_h
#define DetectionSystem8pi_h 1

#include "DetectionSystem8pi.hh"
#include "globals.hh"

#include "G4AssemblyVolume.hh"

class DetectionSystem8pi
{
public:
    DetectionSystem8pi();
    ~DetectionSystem8pi();

    G4int Build() ;//G4SDManager* mySDman);
    G4int BuildOneDetector();
    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate, G4int detector_number);
    //    G4double const GetDetectorLengthOfUnitsCM() {return this->can_length_z;};

    // Assembly volumes
    G4AssemblyVolume* assembly;
    G4AssemblyVolume* assemblyGe;
    G4AssemblyVolume* assemblyInnerBGO;
    G4AssemblyVolume* assemblyOuterLowerBGO;
    G4AssemblyVolume* assemblyOuterUpperBGO;

    G4double cut_clearance;

    G4double detail_view_end_angle;
    G4double dist_from_origin;
    G4double crystal_dist_from_origin;
    G4double crystal_outer_radius;
    G4double crystal_inner_radius;
    G4double hole_starting_depth;
    G4double crystal_length;
    G4double dead_layer_thickness;
    G4double electrode_radius;
    G4double inner_can_thickness;
    G4double inner_can_extends_past_crystal;
    G4double inner_can_lid_thickness;
    G4double inner_can_lid_separation;

    G4double structureMat_cooling_rod_radius;
    G4double electrodeMat_cooling_rod_radius;
    G4double electrodeMat_cooling_rod_length;
    G4double electrodeMat_cooling_rod_dist_from_inner_lid;
    G4double cooling_rod_cover_clearance;
    G4double outer_can_length;
    G4double outer_can_dist_from_inner_can;
    G4double outer_can_thickness;
    G4double outer_can_extends_past_crystal;
    G4double outer_can_innerRadius;
    G4double beryllium_dist_from_crystal;
    G4double beryllium_thickness;
    G4double beryllium_radius;
    G4double inner_BGO_annulus_length;
    G4double inner_BGO_clearance;
    G4double inner_BGO_innerRadius;
    G4double inner_BGO_outerRadius;
    G4double outer_BGO_bottom_thickness;
    G4double outer_BGO_taper_height;
    G4double outer_BGO_top_outer_radius;
    G4double outer_BGO_total_length;
    G4double outer_BGO_clearance;
    G4double outer_BGO_displacement;
    G4double liquid_N2_length;
    G4double liquid_N2_radius;

    G4double hevimetal_thickness;
    G4double hevimetal_front_side_length;
    G4double hevimetal_rear_side_length;
    G4double auxMat_plug_neg_radius;
    G4double auxMat_plug_pos_radius;
    G4double auxMat_layer_thickness;
    G4double auxMat_layer_front_side_length;

    //half z-lengths for various parts
    //used to pass information from logical to physical volumes
    G4double dead_layer_half_length_z;
    G4double germanium_vacuum_core_ZHalfLength;
    G4double lower_electrodeMat_electrode_ZHalfLength;
    G4double upper_electrodeMat_electrode_ZHalfLength;
    G4double inner_cage_ZHalfLength;
    G4double inner_cage_bottom_ZHalfLength;
    G4double inner_cage_lid_ZHalfLength;
    G4double structureMat_cooling_rod_ZHalfLength;
    G4double electrodeMat_cooling_rod_ZHalfLength;
    G4double cooling_rod_cover_ZHalfLength;
    G4double cooling_rod_cover_lid_ZHalfLength;
    G4double outer_can_side_ZHalfLength;
    G4double outer_can_lid_ZHalfLength;
    G4double outer_can_bottom_ZHalfLength;
    G4double beryllium_window_ZHalfLength;
    G4double inner_BGO_annulus_ZHalfLength;
    G4double structureMat_sheath_ZHalfLength;
    G4double outer_lower_BGO_annulus_ZHalfLength;
    G4double outer_upper_BGO_annulus_ZHalfLength;
    G4double liquid_N2_ZHalfLength;
    G4double liquid_N2_side_ZHalfLength;
    G4double liquid_N2_lid_ZHalfLength;
    G4double liquid_N2_bottom_ZHalfLength;

    //suppressed:
    G4double extraClearance;
    G4String electrodeMat;
    G4String structureMat;
    G4String auxMat;

private:
    // Logical volumes
    //
    G4LogicalVolume* exp_hall_log;
    G4LogicalVolume* germanium_block_log;
    G4LogicalVolume* germanium_dead_layer_log;
    G4LogicalVolume* germanium_vacuum_core_log;
    G4LogicalVolume* lower_electrodeMat_electrode_log;
    G4LogicalVolume* upper_electrodeMat_electrode_log;
    G4LogicalVolume* inner_cage_1_log;
    G4LogicalVolume* inner_cage_2_log;
    G4LogicalVolume* inner_cage_3_log;
    G4LogicalVolume* inner_cage_4_log;
    G4LogicalVolume* inner_cage_bottom_log;
    G4LogicalVolume* inner_cage_lid_log;
    G4LogicalVolume* structureMat_cooling_rod_log;
    G4LogicalVolume* electrodeMat_cooling_rod_log;
    G4LogicalVolume* cooling_rod_cover_log;
    G4LogicalVolume* cooling_rod_cover_lid_log;
    G4LogicalVolume* outer_can_side_log;
    G4LogicalVolume* outer_can_lid_log;
    G4LogicalVolume* outer_can_bottom_log;
    G4LogicalVolume* beryllium_window_log;
    G4LogicalVolume* inner_BGO_annulus_log;
    G4LogicalVolume* structureMat_sheath_log;
    G4LogicalVolume* outer_lower_BGO_annulus_log;
    G4LogicalVolume* outer_upper_BGO_annulus_log;
    G4LogicalVolume* liquid_N2_log;
    G4LogicalVolume* liquid_N2_side_log;
    G4LogicalVolume* liquid_N2_lid_log;
    G4LogicalVolume* liquid_N2_bottom_log;
    G4LogicalVolume* hevimetal_log;
    G4LogicalVolume* auxMat_plug_log;
    G4LogicalVolume* auxMat_layer_log;

    //    SensitiveDetector* germanium_block_SD;
    //    SensitiveDetector* inner_BGO_annulus_SD;
    //    SensitiveDetector* outer_lower_BGO_SD;
    //    SensitiveDetector* outer_upper_BGO_SD;

    G4ThreeVector GetDirectionXYZ(G4double, G4double);

    G4int AddGermanium();
    G4int AddGermanium_Logical();
    G4int AddGermanium_Physical();
    G4int AddGermaniumDeadLayer();
    G4int AddGermaniumCore();
    G4int AddElectrodeMatElectrode();
    G4int AddStructureMatCage();
    G4int AddInnerStructureMatLids();
    G4int AddStructureMatCoolingRod();
    G4int AddElectrodeMatCoolingRod();
    G4int AddBerylliumWindow();
    G4int AddOuterStructureMatCan();
    G4int AddCoolingRodCover();
    G4int AddInnerBGOAnnulus();
    G4int AddStructureMatBGOSheath();
    G4int AddOuterBGOAnnulus();
    G4int AddLiquidN2Container();
    G4int AddHevimetalCollimator();
    G4int AddAuxMatPlug();
    G4int AddThinAuxMatLayer();
};

#endif
