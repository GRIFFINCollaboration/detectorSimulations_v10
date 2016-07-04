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

#ifndef DetectionSystemSpiceV02_h
#define DetectionSystemSpiceV02_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#define ALCOL 0.5,0.5,0.5

class DetectionSystemSpiceV02
{
public:
    DetectionSystemSpiceV02();
    ~DetectionSystemSpiceV02();
    
    G4int Build() ; //G4SDManager* mySDman);
    G4int PlaceDetector(G4LogicalVolume* expHallLog, G4ThreeVector move, G4RotationMatrix* rotate, G4int detectorNumber);
    //G4double const GetDetectorLengthOfUnitsCM() {return this->canLengthZ;};

    //-----------------------------//
    // parameters for the square   //
    // planar detector crystal     //
    //-----------------------------//
    G4double fSquareDetCrystalLength;
    G4double fSquareDet90DegCrystalRadius;
    G4double fSquareDet45DegCrystalRadius;
    G4double fSquareDetCrystalThickness;
    G4double fDetectorFaceTargetDistance;

    //------------------------------------------------//
    // logical and physical volumes
    //------------------------------------------------//
private:
    G4AssemblyVolume* fAssembly;
    G4AssemblyVolume* fAssemblySi;
    //    SensitiveDetector* crystalBlockSD;

    G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

    G4LogicalVolume* fDetectorCasingSideLog;
    G4LogicalVolume* fDetectorCasingBackLog;
    G4LogicalVolume* fCrystalBlockLog;


    //--------------------------------------------------------//
    // SPICE physical properties
    // OBS: crystal properties are public, others are private
    //--------------------------------------------------------//

    G4String fCasingMaterial;
    G4String fWaferMaterial;

    //-----------------------------//
    // parameters for the square   //
    // planar detector casing      //
    //-----------------------------//

    G4double fSquareDetCasingLength;
    G4double fSquareDetCasingThickness;

    G4int fSiDetCopyNumber;
    G4int fSiDetCasingSideCopyNumber;
    G4int fSiDetCasingBackCopyNumber;

    //------------------------------------------------//
    // internal methods in Build()
    //------------------------------------------------//

    G4int BuildSiliconWafer();
    G4int BuildAluminiumCasingSide();
    G4int BuildAluminiumCasingBack();

    //    void PlaceSiliconWafer();
    //    void PlaceAluminiumCasingSide();
    //    void PlaceAluminiumCasingBack();

    //    void BuildSPICETargetChamber();

    //    G4ThreeVector Translate45DegDetector(G4int);
    //    G4ThreeVector Translate90DegDetector(G4int);
    //    G4RotationMatrix* RotateDetector(G4int);
    //------------------------------------------------//
    // public methods
    //------------------------------------------------//
    G4Box*              BuildCrystal();
    G4SubtractionSolid* BuildCasingSide();
    G4Box*              BuildCasingBack();
};

#endif

