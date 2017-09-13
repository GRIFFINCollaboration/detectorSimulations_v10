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
    G4int PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber, G4double radialpos, G4int hevimetopt);


private:
    // Logical volumes
    G4LogicalVolume* fDetectorVolumeLog;
    G4LogicalVolume* fBgoBlockLog;
    G4LogicalVolume* fVacuumBlockLog;
    G4LogicalVolume* fCanCylinderLog;
    G4LogicalVolume* fCanCylinderBackCoverLog;
    G4LogicalVolume* fHevimetBlockLog;

    // Assembly volumes
    G4AssemblyVolume* fAssembly;
    G4AssemblyVolume* fAssemblyHevimet;

    G4double fCanLength;
    G4double fCanThickness;
    G4double fCanThicknessFront;
    G4double fCanInnerRadius;
    G4String fCanMaterial;

    G4double fBgoLength;
    G4double fBgoThickness;
    G4String fBgoMaterial;

    G4double fGapThickness;
    G4double fGapThicknessOuter;
    G4String fGapMaterial;

    G4double fCanFaceToBuildOrigin;

    G4double fChoppingOuterAngle;
    G4double fChoppingInnerAngle;

    G4double fHevimetThickness;
    G4String fHevimetMaterial;

    G4double fDetectorAngles[8][5];

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
};

#endif

