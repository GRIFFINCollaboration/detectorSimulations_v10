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

#ifndef ApparatusGenericTarget_h
#define ApparatusGenericTarget_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class ApparatusGenericTarget
{
  public:
    ApparatusGenericTarget();
    ~ApparatusGenericTarget();

  private:
    // Logical volumes        
    G4LogicalVolume* target_log;    

    // Assembly volumes
    G4AssemblyVolume* assembly;

  private: 
    G4String target_material;
    G4double target_length_x;   
    G4double target_length_y;
    G4double target_length_z;

  public: 
    G4int Build(G4String target_material, G4double target_length_x, G4double target_length_y, G4double target_length_z); 
    G4int PlaceApparatus(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate); 

  private: 
    G4Box* BuildTarget();
    
  private: 
    G4int BuildTargetVolume();           
        
  private:
    G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);
};

#endif

