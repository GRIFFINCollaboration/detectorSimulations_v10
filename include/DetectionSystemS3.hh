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

#ifndef DetectionSystemS3_h
#define DetectionSystemS3_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#define ALCOL 0.5,0.5,0.5

class DetectionSystemS3
{
public:
    DetectionSystemS3();
    ~DetectionSystemS3();

    //------------------------------------------------//
    // logical and physical volumes
    //------------------------------------------------//
private:
    G4AssemblyVolume* fAssembly;
    G4AssemblyVolume* fAssemblyS3Ring[24];

public:
    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* expHallLog, G4ThreeVector move,
                        G4int ringNumber, G4int nRadSeg, G4int detectorNumber);
    G4int PlaceGuardRing(G4LogicalVolume* expHallLog, G4ThreeVector move);

private:
    G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

    G4LogicalVolume* fS3InnerGuardRingLog;
    G4LogicalVolume* fS3OuterGuardRingLog;
    G4LogicalVolume* fSiDetS3RingLog[24];

    //--------------------------------------------------------//
    // SPICE physical properties
    // OBS: crystal properties are public, others are private
    //--------------------------------------------------------//
private:
    G4String fWaferMaterial;

    //-----------------------------//
    // parameters for the annular  //
    // planar detector crystal     //
    //-----------------------------//
public:
    G4double fS3DetCrystalThickness;
    G4double fS3DetCrystalOuterDiameter;
    G4double fS3DetCrystalInnerDiameter;
    G4double fS3DetRadialSegments;
    G4double fS3DetPhiSegments;

    //-------------------------------//
    // parameters for the guard ring //
    //-------------------------------//
private:
    G4double fS3DetGuardRingInnerDiameter;
    G4double fS3DetGuardRingOuterDiameter;

    //------------------------------------------------//
    // internal methods in Build()
    //------------------------------------------------//
private:
    G4int BuildSiliconWafer(G4int ringID);
    G4int BuildInnerGuardRing();
    G4int BuildOuterGuardRing();

    G4Tubs*             BuildCrystal(G4int myRingID);
};

#endif
