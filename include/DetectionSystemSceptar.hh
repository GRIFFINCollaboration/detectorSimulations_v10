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

#ifndef DetectionSystemSceptar_h
#define DetectionSystemSceptar_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectionSystemSceptar	
{
public:
    DetectionSystemSceptar();
    ~DetectionSystemSceptar();

    G4int Build() ; //G4SDManager* mySDman);
    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detectorNumber);

private:
    // Assembly volumes
    G4AssemblyVolume* assembly;
    G4AssemblyVolume* assemblySquare;
    G4AssemblyVolume* assemblyAngled;
    G4AssemblyVolume* assemblySquareSD;
    G4AssemblyVolume* assemblyAngledSD;

    //    SensitiveDetector* square_scint_SD;
    //    SensitiveDetector* angled_scint_SD;

    G4double convert;
    G4double square_scintillator_length;
    G4double square_scintillator_width;
    G4double square_scintillator_thickness;
    G4double angled_scintillator_length;
    G4double angled_scintillator_long_width;
    G4double angled_scintillator_short_width;
    G4double angled_scintillator_thickness;
    G4double mylar_thickness;
    G4double scint_gap;
    G4double scint_angle1;
    G4double scint_angle2;
    G4double scint_angle3;
    G4double scint_angle4;
    G4double scint_angle_move;
    G4double square_scint_radial_distance;
    G4double angled_scint_radial_distance;
    G4double angled_scint_move_back;
    G4double Delrin_inner_radius;
    G4double Delrin_outer_radius;
    G4double Delrin2_inner_radius;
    G4double Delrin2_outer_radius;
    G4double Hevimet_inner_radius;
    G4double Hevimet_outer_radius;
    G4double Delrin_hole_radius;
    G4double separate_hemispheres;
    
    //##################################################################
    //### methods to construct all of the components of the detector ###
    //##################################################################
    G4int ConstructScintillator();
    G4int ConstructDelrinShell();
    G4int Construct2ndDelrinShell();
    G4int ConstructHevimetShell();

    //Logical Volumes
    //
    G4LogicalVolume* square_mylar_log;
    G4LogicalVolume* angled_mylar_log;
    G4LogicalVolume* square_scintillator_log;
    G4LogicalVolume* angled_scintillator_log;
    G4LogicalVolume* Delrin_shell_log;
    G4LogicalVolume* Delrin_shell2_log;
    G4LogicalVolume* Hevimet_shell_log;

    G4Trd* squareMylar();
    G4SubtractionSolid* squareMylarWithCut();
    G4SubtractionSolid* angledMylar();
    G4SubtractionSolid* angledMylarWithCut();
    G4Trd* squareScintillator();
    G4SubtractionSolid* angledScintillator();
    G4SubtractionSolid* DelrinShell();
    G4SubtractionSolid* DelrinShell2();
    G4SubtractionSolid* HevimetShell();



};

#endif

