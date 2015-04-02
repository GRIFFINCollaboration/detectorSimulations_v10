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

#ifndef DetectionSystemAncillaryBGO_h
#define DetectionSystemAncillaryBGO_h 1

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units
#include "G4IntersectionSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

class DetectionSystemAncillaryBGO
{
public:
    DetectionSystemAncillaryBGO();
    ~DetectionSystemAncillaryBGO();
    
    G4int Build() ; //G4SDManager* mySDman);
    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector_number, G4double radialpos, G4int hevimetopt);


private:
    // Logical volumes
    G4LogicalVolume* detector_volume_log;
    G4LogicalVolume* bgo_block_log;
    G4LogicalVolume* vacuum_block_log;
    G4LogicalVolume* can_cylinder_log;
    G4LogicalVolume* can_cylinder_back_cover_log;
    G4LogicalVolume* hevimet_block_log;

    // Assembly volumes
    G4AssemblyVolume* assembly;
    G4AssemblyVolume* assemblyHevimet;

    G4double can_length;
    G4double can_thickness;
    G4double can_thickness_front;
    G4double can_inner_radius;
    G4String can_material;

    G4double bgo_length;
    G4double bgo_thickness;
    G4String bgo_material;

    G4double gap_thickness;
    G4double gap_thickness_outer;
    G4String gap_material;

    G4double can_face_to_build_origin;

    G4double chopping_outer_angle;
    G4double chopping_inner_angle;

    G4double hevimet_thickness;
    G4String hevimet_material;

    G4double detectorAngles[8][5];

    G4SubtractionSolid* BGOPiece();
    G4SubtractionSolid* VacuumPiece();
    G4SubtractionSolid* AluminumCanPiece();
    G4Tubs*             AluminumCanBackCover();
    G4SubtractionSolid* HevimetPiece();

    G4int BuildBGOVolume();
    G4int BuildVacuumVolume();
    G4int BuildAluminumCanVolume();
    G4int BuildAluminumCanBackCoverVolume();
    G4int BuildHevimetPiece();



    G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

    G4double transX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi);
};

#endif

