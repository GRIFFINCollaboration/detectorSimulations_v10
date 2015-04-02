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

#ifndef DetectionSystemLanthanumBromide_h
#define DetectionSystemLanthanumBromide_h 1

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

const G4double inchtocm = 2.54;

class DetectionSystemLanthanumBromide
{
public:
    DetectionSystemLanthanumBromide();
    ~DetectionSystemLanthanumBromide();
    
    G4int Build() ; //G4SDManager* mySDman);
    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector_number, G4double radialpos);
    G4double const GetDetectorLengthOfUnitsCM() {return this->detector_length_z;};
    G4double GetCrystalRadius() {return this->crystal_outer_radius;};
    G4double GetCrystalLength() {return this->crystal_length_z;};
    G4double GetR();
    G4double GetTheta(G4int i);
    G4double GetPhi(G4int i);
    G4double GetYaw(G4int i);
    G4double GetPitch(G4int i);
    G4double GetRoll(G4int i);
    G4double GetCrystalRadialPosition() {return this->packing_front_lid_thickness+this->disc_front_lid_thickness+this->seal_front_lid_thickness+this->can_front_lid_thickness+12.5*cm ;};

private:
    // Logical volumes
    G4LogicalVolume* detector_volume_log;
    G4LogicalVolume* crystal_block_log;

    G4LogicalVolume* can_cylinder_log;
    G4LogicalVolume* can_front_lid_log;
    G4LogicalVolume* can_back_lid_log;

    G4LogicalVolume* seal_front_lid_log;
    G4LogicalVolume* disc_front_lid_log;

    G4LogicalVolume* packing_cylinder_log;
    G4LogicalVolume* packing_front_lid_log;

    // Assembly volumes
    G4AssemblyVolume* assembly;

    //    SensitiveDetector* crystal_block_SD;

    G4int copy_number;
    G4int number_of_detectors;
    G4int number_of_segments;

    G4String crystal_material;
    G4String can_material;
    G4String seal_material;
    G4String disc_material;
    G4String packing_material;

    G4double detail_view_end_angle;
    G4double crystal_length_z;
    G4double crystal_inner_radius;
    G4double crystal_outer_radius;
    G4double packing_length_z;
    G4double packing_inner_radius;
    G4double packing_outer_radius;
    G4double packing_lid_inner_radius;
    G4double packing_lid_outer_radius;
    G4double packing_front_lid_thickness;
    G4double disc_lid_inner_radius;
    G4double disc_lid_outer_radius;
    G4double disc_front_lid_thickness;
    G4double seal_lid_inner_radius;
    G4double seal_lid_outer_radius ;
    G4double seal_front_lid_thickness;
    G4double can_length_z;
    G4double can_inner_radius;
    G4double can_outer_radius;
    G4double can_lid_inner_radius;
    G4double can_lid_outer_radius;
    G4double can_front_lid_thickness;
    G4double can_back_lid_thickness;

    G4double detector_length_z;
    G4double detectorAngles[8][5];

    G4double set_radial_pos;

    G4Tubs* BuildCrystal();
    G4Tubs* BuildAluminumCan();
    G4Tubs* BuildAluminumCanFrontLid();
    G4Tubs* BuildAluminumCanBackLid();
    G4Tubs* BuildPacking();
    G4Tubs* BuildPackingFrontLid();
    G4Tubs* BuildDiscFrontLid();
    G4Tubs* BuildSealFrontLid();
    
    G4int BuildSealVolume();
    G4int BuildDiscVolume();
    G4int BuildPackingVolume();
    G4int BuildAluminumCanVolume();
    G4int BuildCrystalVolume();
    G4int BuildOneDetector();

    G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

    G4double transX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi);
};

#endif

