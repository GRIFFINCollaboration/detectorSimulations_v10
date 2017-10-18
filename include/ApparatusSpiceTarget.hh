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
    //G4LogicalVolume* expHallLog;
    G4LogicalVolume* fTargetSpiceLog;
    G4LogicalVolume* fTargetBackerSpiceLog;
    G4LogicalVolume* fTargetProtectorSpiceLog;
    G4LogicalVolume* fTargetBracketLog;
    G4LogicalVolume* fTargetHolderLog;
    
    //Physical Volumes
    G4VPhysicalVolume* fTargetPhys;
    G4VPhysicalVolume* fTargetBackerPhys;
    G4VPhysicalVolume* fTargetProtectorPhys;
    G4VPhysicalVolume* fTargetBracketPhys;
    G4VPhysicalVolume* fTargetHolderPhys;
    
    //Material manager
//  G4NistManager* man = G4NistManager::Instance();
    
  private://variables for build, set in the class as read in
    
    //Target (middle)
    G4String fTargetMaterial;
    G4double fTargetMaterialDensity;
    G4double fTargetRadius;//this value is equal to the target frame hole radius //used for backer, and protector also
    G4double fTargetThickness;
    G4double fTargetSurfaceDensity;
    
    //Backer
    G4String fTargetBackerMaterial;
    G4double fTargetBackerMaterialDensity;
    G4double fTargetBackerThickness;
    G4double fTargetBackerSurfaceDensity;
    
    //Protector:
    G4String fTargetProtectorMaterial;
    G4double fTargetProtectorMaterialDensity;
    G4double fTargetProtectorThickness;
    G4double fTargetProtectorSurfaceDensity;
    
    //Bracket
    G4String fBracketMaterial;
    G4double fBracketSideLength;
    G4double fBracketBackLength;
    G4double fBracketBackTopThickness;//x/y-axis
    G4double fBracketSideTopThickness;
    G4double fBracketDepth;//z-axis
    G4double fScrewRadius;
    G4String fBracketScrewMaterial;
    
    //Holder
    G4String fHolderMaterial;
    G4double fHolderOuterRadii;
    G4double fHolderInnerRadii;
    G4double fHolderThickness;
    
  public: 
    G4int BuildTarget(G4String fTargetMaterial, G4double fTargetThickness, G4double fTargetMaterialDensity);
    G4int BuildBacker(G4String fTargetBackerMaterial, G4double fTargetBackerThickness, G4double fTargetBackerMaterialDensity);
    G4int BuildProtector(G4String fTargetProtectorMaterial, G4double fTargetProtectorThickness, G4double fTargetProtectorMaterialDensity);

    void  PlaceTarget(G4LogicalVolume*);
    void  PlaceTargetBacker(G4LogicalVolume*);
    void  PlaceTargetProtector(G4LogicalVolume*);
    
    G4int BuildBracket();
    void  PlaceBracket(G4LogicalVolume*);
    G4int BuildHolder();
    void  PlaceHolder(G4LogicalVolume*);
    G4double fBeamPos;//uses beam position to place targets
             

};

#endif