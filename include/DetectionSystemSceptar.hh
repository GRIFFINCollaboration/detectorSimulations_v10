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
    G4int PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber);

private:
    // Assembly volumes
    G4AssemblyVolume* fAssembly;
    G4AssemblyVolume* fAssemblySquare;
    G4AssemblyVolume* fAssemblyAngled;
    G4AssemblyVolume* fAssemblySquareSD;
    G4AssemblyVolume* fAssemblyAngledSD;

    //    SensitiveDetector* squareScintSD;
    //    SensitiveDetector* angledScintSD;

    G4double fConvert;
    G4double fSquareScintillatorLength;
    G4double fSquareScintillatorWidth;
    G4double fSquareScintillatorThickness;
    G4double fAngledScintillatorLength;
    G4double fAngledScintillatorLongWidth;
    G4double fAngledScintillatorShortWidth;
    G4double fAngledScintillatorThickness;
    G4double fMylarThickness;
    G4double fScintGap;
    G4double fScintAngle1;
    G4double fScintAngle2;
    G4double fScintAngle3;
    G4double fScintAngle4;
    G4double fScintAngleMove;
    G4double fSquareScintRadialDistance;
    G4double fAngledScintRadialDistance;
    G4double fAngledScintMoveBack;
    G4double fDelrinInnerRadius;
    G4double fDelrinOuterRadius;
    G4double fDelrin2InnerRadius;
    G4double fDelrin2OuterRadius;
    G4double fHevimetInnerRadius;
    G4double fHevimetOuterRadius;
    G4double fDelrinHoleRadius;
    G4double fSeparateHemispheres;
    
    //##################################################################
    //### methods to construct all of the components of the detector ###
    //##################################################################
    G4int ConstructScintillator();
    G4int ConstructDelrinShell();
    G4int Construct2ndDelrinShell();
    G4int ConstructHevimetShell();

    //Logical Volumes
    //
    G4LogicalVolume* fSquareMylarLog;
    G4LogicalVolume* fAngledMylarLog;
    G4LogicalVolume* fSquareScintillatorLog;
    G4LogicalVolume* fAngledScintillatorLog;
    G4LogicalVolume* fDelrinShellLog;
    G4LogicalVolume* fDelrinShell2Log;
    G4LogicalVolume* fHevimetShellLog;

    G4Trd* SquareMylar();
    G4SubtractionSolid* SquareMylarWithCut();
    G4SubtractionSolid* AngledMylar();
    G4SubtractionSolid* AngledMylarWithCut();
    G4Trd* SquareScintillator();
    G4SubtractionSolid* AngledScintillator();
    G4SubtractionSolid* DelrinShell();
    G4SubtractionSolid* DelrinShell2();
    G4SubtractionSolid* HevimetShell();

};

#endif

