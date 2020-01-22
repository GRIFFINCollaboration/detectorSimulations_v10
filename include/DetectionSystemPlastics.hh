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

#ifndef DETECTIONSYSTEMPLASTICS_HH
#define DETECTIONSYSTEMPLASTICS_HH

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

class DetectionSystemPlastics
{
public:
    DetectionSystemPlastics(G4double thickness, G4int material, G4double numDet);
    ~DetectionSystemPlastics();

    void SetWrapping(G4bool wrap){fAddWrap = wrap;};
    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* expHallLog);

private:

    // Assembly volumes
    G4AssemblyVolume* fAssemblyPlastics;                 // Contains all non-sensitive materials

    G4double fScintillatorLength;
    G4double fScintillatorHeight;
    G4double fScintillatorWidth;

    G4double fRadialDistance;
    G4double fLeadShieldThickness;
    G4double fSpacing;
    G4double fWrapThickness;
    G4double fAirGap;
   
    G4bool fAddWrap;

    G4double fNumDet;
    G4double fNumDetBot;

    G4String fPMTMaterial;
    G4String fWrapMaterial;
    G4String fPlasticMaterial;
   
    std::vector<G4LogicalVolume*>  fPlasticLogArray;
    std::vector<G4LogicalVolume*>  fWrapLogArray;
    std::vector<G4LogicalVolume*>  fPMT1LogArray;
    std::vector<G4LogicalVolume*>  fPMT2LogArray;
   
    G4Color bronze;
    G4Color black;
    G4Color silver;

    G4int BuildPlastics();

};

#endif

