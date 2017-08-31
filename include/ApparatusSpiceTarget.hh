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

#ifndef ApparatusSpiceTarget_h
#define ApparatusSpiceTarget_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
//forward declare for placement
class G4VPhysicalVolume;

class ApparatusSpiceTarget
{
  public:
    ApparatusSpiceTarget(G4double);
    ~ApparatusSpiceTarget();

  private:
    // Logical volumes
    G4LogicalVolume* expHallLog;
    G4LogicalVolume* target_spice_log;
    G4LogicalVolume* target_Backer_spice_log;
    G4LogicalVolume* target_Protector_spice_log;
    
    //Physical Volumes
    G4VPhysicalVolume* target_phys;
    G4VPhysicalVolume* target_Backer_phys;
    G4VPhysicalVolume* target_Protector_phys;
  private://variables for build, set in the class as read in
    G4String target_material;
    G4double target_material_density;
    G4double target_radius = 9.0*mm;//this value is equal to the target frame hole radius //used for backer, and protector also
    G4double target_thickness;
    G4double target_surface_density;
    G4String target_Backer_material;
    G4double target_Backer_material_density;
    G4double target_Backer_thickness;
    G4double target_Backer_surface_density;
    G4String target_Protector_material;
    G4double target_Protector_material_density;
    G4double target_Protector_thickness;
    G4double target_Protector_surface_density;
    
    
  public: 
    G4int BuildTarget(G4String target_material, G4double target_thickness, G4double target_material_density);
    G4int BuildBacker(G4String target_Backer_material, G4double target_Backer_thickness, G4double target_Backer_material_density);
    G4int BuildProtector(G4String target_Protector_material, G4double target_Protector_thickness, G4double target_Protector_material_density);
    void PlaceTarget(G4LogicalVolume* expHallLog);
    void PlaceTargetBacker(G4LogicalVolume* expHallLog);
    void PlaceTargetProtector(G4LogicalVolume* expHallLog);
    
    G4double beampos;
             

};

#endif