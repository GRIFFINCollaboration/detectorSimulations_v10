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

#ifndef DetectionSystemSodiumIodide_h
#define DetectionSystemSodiumIodide_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

const G4double inch2cm = 2.54;

class DetectionSystemSodiumIodide
{
public:
    DetectionSystemSodiumIodide();
    ~DetectionSystemSodiumIodide();
    
    G4int Build() ; //G4SDManager* mySDman);
    G4int PlaceDetector(G4LogicalVolume* expHallLog, G4ThreeVector move, G4RotationMatrix* rotate, G4int detectorNumber);
    G4double GetDetectorLengthOfUnitsCM() { return fDetectorLengthZ; }

private:
    // Logical volumes
    G4LogicalVolume* fDetectorVolumeLog;
    G4LogicalVolume* fCrystalBlockLog;

    G4LogicalVolume* fCanCylinderLog;
    G4LogicalVolume* fCanFrontLidLog;
    G4LogicalVolume* fCanBackLidLog;

    G4LogicalVolume* fSealFrontLidLog;
    G4LogicalVolume* fDiscFrontLidLog;

    G4LogicalVolume* fPackingCylinderLog;
    G4LogicalVolume* fPackingFrontLidLog;

    // Assembly volumes
    G4AssemblyVolume* fAssembly;

    //    SensitiveDetector* crystalBlockSD;

    G4int fCopyNumber;
    G4int fNumberOfDetectors;
    G4int fNumberOfSegments;

    G4String fCrystalMaterial;
    G4String fCanMaterial;
    G4String fSealMaterial;
    G4String fDiscMaterial;
    G4String fPackingMaterial;

    G4double fDetailViewEndAngle;
    G4double fCrystalLengthZ;
    G4double fCrystalInnerRadius;
    G4double fCrystalOuterRadius;
    G4double fPackingLengthZ;
    G4double fPackingInnerRadius;
    G4double fPackingOuterRadius;
    G4double fPackingLidInnerRadius;
    G4double fPackingLidOuterRadius;
    G4double fPackingFrontLidThickness;
    G4double fDiscLidInnerRadius;
    G4double fDiscLidOuterRadius;
    G4double fDiscFrontLidThickness;
    G4double fSealLidInnerRadius;
    G4double fSealLidOuterRadius ;
    G4double fSealFrontLidThickness;
    G4double fCanLengthZ;
    G4double fCanInnerRadius;
    G4double fCanOuterRadius;
    G4double fCanLidInnerRadius;
    G4double fCanLidOuterRadius;
    G4double fCanFrontLidThickness;
    G4double fCanBackLidThickness;

    G4double fDetectorLengthZ;

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
};

#endif

