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
//#include "G4IntersectionSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SubtractionSolid.hh"

class G4AssemblyVolume;
class DetectionSystemDescant
{
public:
    DetectionSystemDescant();
    ~DetectionSystemDescant();

    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector_number);

private:
    // Logical volumes
    G4LogicalVolume* blue_volume_log;
    G4LogicalVolume* green_volume_log;
    G4LogicalVolume* red_volume_log;
    G4LogicalVolume* white_volume_log;
    G4LogicalVolume* yellow_volume_log;

    G4LogicalVolume* blue_scintillator_volume_log;
    G4LogicalVolume* green_scintillator_volume_log;
    G4LogicalVolume* red_scintillator_volume_log;
    G4LogicalVolume* white_scintillator_volume_log;
    G4LogicalVolume* yellow_scintillator_volume_log;

    G4LogicalVolume* blue_lead_volume_log;
    G4LogicalVolume* green_lead_volume_log;
    G4LogicalVolume* red_lead_volume_log;
    G4LogicalVolume* white_lead_volume_log;
    G4LogicalVolume* yellow_lead_volume_log;

    // Assembly volumes
    G4AssemblyVolume* assemblyBlue;                 // Contains all non-sensitive materials
    G4AssemblyVolume* assemblyBlueScintillator;     // Contains all only sensitive materials, eg. the scintillation material
    G4AssemblyVolume* assemblyGreen;
    G4AssemblyVolume* assemblyGreenScintillator;
    G4AssemblyVolume* assemblyRed;
    G4AssemblyVolume* assemblyRedScintillator;
    G4AssemblyVolume* assemblyWhite;
    G4AssemblyVolume* assemblyWhiteScintillator;
    G4AssemblyVolume* assemblyYellow;
    G4AssemblyVolume* assemblyYellowScintillator;

    G4double can_length;
    G4double can_thickness;
    G4double can_thickness_front;
    G4double can_inner_radius;
    G4double lead_shield_thickness;
    G4double radial_distance;
    G4String can_material;
    G4String liquid_material;
    G4String lead_material;

    // Saint Gobain data files, 6 points around the front face of the can, and 6 on the back
    // from Dan Brennan
    G4double blue_detector[12][3];
    G4double green_detector[12][3];
    G4double red_detector[12][3];
    G4double white_detector[12][3];
    G4double yellow_detector[12][3];

    // These are the angles of the cuts for each of the 6 sides of the detectors
    G4double blue_phi[6];
    G4double green_phi[6];
    G4double red_phi[6];
    G4double white_phi[6];
    G4double yellow_phi[6];

    // The Euler angles from James' MSc thesis which gives us the detector positions
    // Some of the angles for the green and yellow detectors are wrong in James' thesis,
    // note the +180 on a few angles.
    G4double blue_alpha_beta_gamma[15][3];
    G4double green_alpha_beta_gamma[10][3];
    G4double red_alpha_beta_gamma[15][3];
    G4double white_alpha_beta_gamma[20][3];
    G4double yellow_alpha_beta_gamma[10][3];

    // The colours of the detectors
    G4Colour blue_colour;
    G4Colour green_colour;
    G4Colour red_colour;
    G4Colour white_colour;
    G4Colour yellow_colour;
    G4Colour liquid_colour; // Scintillator colour

    G4int BuildCanVolume();
    G4int BuildDetectorVolume();

    G4bool includeLead;


    G4SubtractionSolid* CanVolume(G4bool insideVol, G4double volume_length, G4double detector[12][3], G4double detector_phi[6]);
    G4SubtractionSolid* CutVolumeOnFourPoints(G4int idx, G4bool insideVol, G4double volume_length, G4double detector_phi[6], G4SubtractionSolid* volume, G4ThreeVector front_p1, G4ThreeVector front_p2, G4ThreeVector back_p1, G4ThreeVector back_p2);

    G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

    G4double transX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi);

    G4ThreeVector SolveLineEquationX(G4ThreeVector p1, G4ThreeVector p2, G4double x);
    G4ThreeVector SolveLineEquationY(G4ThreeVector p1, G4ThreeVector p2, G4double y);
    G4ThreeVector SolveLineEquationZ(G4ThreeVector p1, G4ThreeVector p2, G4double z);
};

#endif

