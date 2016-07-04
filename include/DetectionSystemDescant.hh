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

#ifndef DetectionSystemDescant_h
#define DetectionSystemDescant_h 1

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SubtractionSolid.hh"

class G4AssemblyVolume;

class DetectionSystemDescant
{
public:
    DetectionSystemDescant(G4bool leadShield);
    ~DetectionSystemDescant();

    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber);
    G4int PlaceDetectorAuxPorts(G4LogicalVolume* expHallLog, G4int detectorNumber, G4double radialpos);

    G4int PlaceDetector(G4LogicalVolume* expHallLog, G4String color, G4ThreeVector pos, G4ThreeVector rot);
  
private:
    // Logical volumes
    G4LogicalVolume* fBlueVolumeLog;
    G4LogicalVolume* fGreenVolumeLog;
    G4LogicalVolume* fRedVolumeLog;
    G4LogicalVolume* fWhiteVolumeLog;
    G4LogicalVolume* fYellowVolumeLog;

    G4LogicalVolume* fBlueScintillatorVolumeLog;
    G4LogicalVolume* fGreenScintillatorVolumeLog;
    G4LogicalVolume* fRedScintillatorVolumeLog;
    G4LogicalVolume* fWhiteScintillatorVolumeLog;
    G4LogicalVolume* fYellowScintillatorVolumeLog;

    G4LogicalVolume* fBlueLeadVolumeLog;
    G4LogicalVolume* fGreenLeadVolumeLog;
    G4LogicalVolume* fRedLeadVolumeLog;
    G4LogicalVolume* fWhiteLeadVolumeLog;
    G4LogicalVolume* fYellowLeadVolumeLog;
    
    G4LogicalVolume* fQuartzWindow5inchLog;
    G4LogicalVolume* fQuartzWindow3inchLog;   

    // Assembly volumes
    G4AssemblyVolume* fAssemblyBlue;                 // Contains all non-sensitive materials
    G4AssemblyVolume* fAssemblyBlueScintillator;     // Contains all only sensitive materials, eg. the scintillation material
    G4AssemblyVolume* fAssemblyGreen;
    G4AssemblyVolume* fAssemblyGreenScintillator;
    G4AssemblyVolume* fAssemblyRed;
    G4AssemblyVolume* fAssemblyRedScintillator;
    G4AssemblyVolume* fAssemblyWhite;
    G4AssemblyVolume* fAssemblyWhiteScintillator;
    G4AssemblyVolume* fAssemblyYellow;
    G4AssemblyVolume* fAssemblyYellowScintillator;

    G4double fCanLength;
    G4double fCanThickness;
    G4double fCanThicknessFront;
    G4double fCanInnerRadius;
    G4double fLeadShieldThickness;
    G4double fRadialDistance;
    G4String fCanMaterial;
    G4String fLiquidMaterial;
    G4String fLeadMaterial;
    G4bool   fIncludeLead;



    G4double fCanBackThickness;



    // for making PMT cuts
    G4double fStartPhi;
    G4double fEndPhi;
    G4double fInnerRadiusWindow;
    G4double fOuterRadiusWindow5inch;
    G4double fOuterRadiusWindow3inch;
    G4double fHalfLengthZWindowCut;
    G4double fWhitePMTOffsetY;
    G4double fBluePMTOffsetX;
    G4double fRedPMTOffsetX;
    G4double fYellowGreenPMTOffsetX;
    G4double fYellowGreenPMTOffsetY;
    G4String fQuartzMaterial;
    G4double fOpticalWindowHalfThickness;

    // Saint Gobain data files, 6 points around the front face of the can, and 6 on the back
    // from Dan Brennan
    G4double fBlueDetector[12][3];
    G4double fGreenDetector[12][3];
    G4double fRedDetector[12][3];
    G4double fWhiteDetector[12][3];
    G4double fYellowDetector[12][3];

    G4double fTrimGreenY;
    G4double fTrimYellowY;
    G4double fTrimBlueX;

    // These are the angles of the cuts for each of the 6 sides of the detectors
    G4double fBluePhi[6];
    G4double fGreenPhi[6];
    G4double fRedPhi[6];
    G4double fWhitePhi[6];
    G4double fYellowPhi[6];

    // The Euler angles from James' MSc thesis which gives us the detector positions
    // Some of the angles for the green and yellow detectors are wrong in James' thesis,
    // note the +180 on a few angles.
    G4double fBlueAlphaBetaGamma[15][3];
    G4double fGreenAlphaBetaGamma[10][3];
    G4double fRedAlphaBetaGamma[15][3];
    G4double fWhiteAlphaBetaGamma[20][3];
    G4double fYellowAlphaBetaGamma[10][3];

    // for LaBr3 detector locations
    G4double fDetectorAngles[8][5];
    G4double fSetRadialPos;

    // The colours of the detectors
    G4Colour fBlueColour;
    G4Colour fGreenColour;
    G4Colour fRedColour;
    G4Colour fWhiteColour;
    G4Colour fYellowColour;
    G4Colour fLiquidColour; // Scintillator colour
    G4Colour fGreyColour;   // stainless steel pmt tube colour
    G4Colour fMagentaColour;  // optical window colour
    G4Colour fBlackColour;

    G4bool fSurfCheck;

    G4int BuildCanVolume();
    G4int BuildDetectorVolume();

    G4SubtractionSolid* CanVolume(G4bool insideVol, G4double volumeLength, G4double detector[12][3], G4double detectorPhi[6]);
    G4SubtractionSolid* CutVolumeOnFourPoints(G4int idx, G4bool insideVol, G4double volumeLength, G4double detectorPhi[6], G4SubtractionSolid* volume, G4ThreeVector frontP1, G4ThreeVector frontP2);

    G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

    G4ThreeVector SolveLineEquationX(G4ThreeVector p1, G4ThreeVector p2, G4double x);
    G4ThreeVector SolveLineEquationY(G4ThreeVector p1, G4ThreeVector p2, G4double y);
    G4ThreeVector SolveLineEquationZ(G4ThreeVector p1, G4ThreeVector p2, G4double z);
};

#endif

