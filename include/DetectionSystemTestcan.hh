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

#ifndef DetectionSystemTestcan_h
#define DetectionSystemTestcan_h 1

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

class DetectionSystemTestcan
{
public:
    DetectionSystemTestcan(G4double length, G4double radius);
    ~DetectionSystemTestcan();

    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* expHallLog);

private:
    // Logical volumes
    G4LogicalVolume* fTestcanAlumCasingLog;
    G4LogicalVolume* fTestcanScintillatorLog;   
    G4LogicalVolume* fTestcanQuartzWindowLog;

    // Assembly volumes
    G4AssemblyVolume* fAssemblyTestcan;                 // Contains all non-sensitive materials

    G4double fScintillatorLength;
    G4double fScintillatorInnerRadius;
    G4double fScintillatorOuterRadius;

    G4double fAlumCanThickness;
    G4double fQuartzThickness;
    G4double fQuartzRadius;

    G4double fStartPhi;
    G4double fEndPhi;

    G4String fCanMaterial;
    G4String fLiquidMaterial;
    G4String fQuartzMaterial;

    // The colours of the detectors
    G4Colour fLiquidColour; // Scintillator colour
    G4Colour fGreyColour;   // can colour
    G4Colour fQuartzColour;

    G4int BuildTestcan();

};

#endif

