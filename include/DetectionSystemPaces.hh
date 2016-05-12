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
// $Id: DetectionSystemPaces.hh,v 1.1 2012-11-14 14:00:00 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectionSystemPaces_h
#define DetectionSystemPaces_h 1

#include "DetectionSystemPaces.hh"
#include "globals.hh"

class DetectionSystemPaces
{
public:
    DetectionSystemPaces();
    ~DetectionSystemPaces();

    // Assembly volumes
    G4AssemblyVolume* assembly;
    G4AssemblyVolume* assemblyDetector;
    G4AssemblyVolume* assemblySilicon;

private:
    // Logical volumes
    G4LogicalVolume* aluminum_hemisphere_log;
    G4LogicalVolume* aluminum_annulus_top_log;
    G4LogicalVolume* aluminum_annulus_bot_log;
    G4LogicalVolume* canister_log;
    G4LogicalVolume* silicon_block_log;
    G4LogicalVolume* silicon_dead_layer_log;
    //G4LogicalVolume* screw_log;
    G4LogicalVolume* teflon_annulus_top_log;
    G4LogicalVolume* teflon_annulus_bot_log;
    G4LogicalVolume* delrin_hemisphere_log;

    //SensitiveDetector* silicon_block_SD;

public:
    G4double cut_clearance;

    // Solid parameters
    G4double aluminum_hemisphere_inner_radius;
    G4double aluminum_hemisphere_outer_radius;
    G4double aluminum_hemisphere_beam_hole_radius;
    G4double aluminum_hemisphere_beam_hole_rim_height;
    G4double aluminum_annulus_top_inner_radius;
    G4double aluminum_annulus_top_outer_radius;
    G4double aluminum_annulus_top_thickness;
    G4double aluminum_annulus_bot_inner_radius;
    G4double aluminum_annulus_bot_outer_radius;
    G4double aluminum_annulus_bot_thickness;
    G4double canister_inner_radius;
    G4double canister_outer_radius;
    G4double canister_thickness;
    G4double silicon_block_radius;
    G4double silicon_block_thickness;
    G4double silicon_dead_layer_thickness;
    G4double screw_radius;
    G4double screw_placement_radius;
    //G4double screw_length;
    G4double teflon_annulus_top_inner_radius;
    G4double teflon_annulus_top_outer_radius;
    G4double teflon_annulus_top_thickness;
    G4double teflon_annulus_bot_inner_radius;
    G4double teflon_annulus_bot_outer_radius;
    G4double teflon_annulus_bot_thickness;
    G4double delrin_hemisphere_inner_radius;
    G4double delrin_hemisphere_outer_radius;
    G4double delrin_hemisphere_beam_hole_radius;

    G4double aluminum_hemisphere_dist;
    G4double delrin_hemisphere_dist;

    // distances from FRONT of detector facing the source, for assembly
    G4double aluminum_annulus_top_dist; // will be half the annulus thickness
    G4double aluminum_annulus_bot_dist;
    G4double canister_dist;
    G4double silicon_block_dist;
    G4double silicon_dead_layer_front_dist;
    G4double silicon_dead_layer_back_dist;
    G4double teflon_annulus_top_dist;
    G4double teflon_annulus_bot_dist;

    // detector placement and orientation scalars
    G4double paces_placement_distance[5];
    G4double paces_placement_theta[5];
    G4double paces_placement_phi[5];
    G4double paces_orientation_theta[5];
    G4double paces_orientation_phi[5];

    // assembly volume sizes
    G4double detector_assembly_radius;
    G4double detector_assembly_thickness;

public:
    G4int Build();//G4SDManager* mySDman);
    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4int ndet);

private:
    // Construction methods
    G4int AddAluminumHemisphere();
    G4int AddAluminumAnnulusTop();
    G4int AddAluminumAnnulusBot();
    G4int AddCanister();
    G4int AddSiliconBlock();
    G4int AddSiliconDeadLayer();
    G4int AddScrews();
    G4int AddTeflonAnnulusTop();
    G4int AddTeflonAnnulusBot();
    G4int AddDelrinHemisphere();
    G4int CombineAssemblySilicon();
    G4int CombineAssemblyDetector();
    G4int CombineAssembly();
};

#endif
