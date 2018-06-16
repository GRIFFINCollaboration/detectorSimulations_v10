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

    void Build(G4String, G4double, G4double, G4double); 
    void PlaceApparatus(G4LogicalVolume*, G4ThreeVector, G4RotationMatrix*); 

  private:
    G4Box* BuildTarget();
    void BuildTargetVolume();           
    G4ThreeVector GetDirectionXYZ(G4double, G4double);

    // Logical volumes        
    G4LogicalVolume* fTargetLog;    

    // Assembly volumes
    G4AssemblyVolume* fAssembly;

    G4String fTargetMaterial;
    G4double fTargetLengthX;   
    G4double fTargetLengthY;
    G4double fTargetLengthZ;
};

#endif

