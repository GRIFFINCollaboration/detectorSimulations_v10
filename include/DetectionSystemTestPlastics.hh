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

#ifndef DETECTIONSYSTEMTESTPLASTICS_HH
#define DETECTIONSYSTEMTESTPLASTICS_HH

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

class DetectionSystemTestPlastics
{
public:
    DetectionSystemTestPlastics(G4double thickness, G4int material, G4double numDet);
    ~DetectionSystemTestPlastics();

    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* expHallLog);

private:

    // Assembly volumes
    G4AssemblyVolume* fAssemblyTestPlastics;                


    G4double fWrapThickness;
    G4double fAirGap;

    G4String fAlMaterial;
    G4String fPMTMaterial;
    G4String fWrapMaterial;
    G4String fPlasticMaterial;
    G4String fZDSMaterial;
   
    G4LogicalVolume*  fPlasticLog;
    G4LogicalVolume*  fWrapLog;
    G4LogicalVolume*  fPMTLog;
    G4LogicalVolume*  fPMT1Log;
    G4LogicalVolume*  fPMT2Log;
    G4LogicalVolume*  fPMT3Log;
    G4LogicalVolume*  fPMT4Log;
    G4LogicalVolume*  fZDSLog;
    G4LogicalVolume*  fZDSPMTLog;
    G4LogicalVolume*  fSourceLog;
    G4LogicalVolume*  fAlLog;
    
    G4double fStartPhi;
    G4double fDeltaPhi;
    G4double fDiameter;
    G4double fRadialDistance;
    G4double fScintillatorWidth;
    G4double fPMTWidth;
   
   
    G4Color bronze;
    G4Color black;
    G4Color silver;

    G4int BuildTestPlastics();

};

#endif

