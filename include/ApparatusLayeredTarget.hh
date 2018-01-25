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

#ifndef ApparatusLayeredTarget_HH
#define ApparatusLayeredTarget_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

//forward declare for placement
class G4VPhysicalVolume;
class G4UserLimits;
class ApparatusLayeredTarget
{
  public:
	  
    ApparatusLayeredTarget(G4double);
    ~ApparatusLayeredTarget();

  private:
   std::vector< G4LogicalVolume* > fTargetLayerLog;
   std::vector< G4VPhysicalVolume* > fTargetLayerPhys;
   std::vector< G4double > fTargetLayerThick;
    
  private://variables for build, set in the class as read in
    
    G4double fTargetRadius;
    
    //Independent stepping
    G4UserLimits* fStepLimit;
    
  public: 
    G4int BuildTargetLayer(G4String LayerMaterial, G4double LayerArealDensity);
    void  PlaceTarget(G4LogicalVolume*);
    G4double LayerStart(int);
    
    G4double fBeamPos;//uses beam position to place targets
             

};

#endif