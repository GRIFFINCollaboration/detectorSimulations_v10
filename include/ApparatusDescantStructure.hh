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

#ifndef ApparatusDescantStructure_h
#define ApparatusDescantStructure_h 1

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units
//#include "G4IntersectionSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SubtractionSolid.hh"
//#include "G4Sphere.hh"

class G4AssemblyVolume;
class ApparatusDescantStructure
{
public:
    ApparatusDescantStructure();
    ~ApparatusDescantStructure();

    G4int Build();

    G4int PlaceDescantStructure(G4LogicalVolume* exp_hall_log);

private:
    // Logical volumes
    G4LogicalVolume* fDescantStructureLog;

    // Assembly volumes
    G4AssemblyVolume* fAssemblyDescantStructure;

    G4String fStructureMaterial;
    G4double fStructureInnerRadius;
    G4double fStructureOuterRadius;
    G4double fStartPhi;
    G4double fDeltaPhi;
    G4double fStartTheta;
    G4double fDeltaTheta;
    G4double fStructureRadialDistance;

    //for subtraction polyhedra
    G4int fNumSides;
    G4int fNumZplanes;
    G4double fZPlanes[2];

    G4double fPolyRadiusRed;
    G4double fInnerRadiusRed[2];
    G4double fOuterRadiusRed[2];

    G4double fPolyRadiusGreenYellow;
    G4double fInnerRadiusGreenYellow[2];
    G4double fOuterRadiusGreenYellow[2];

    G4double fPolyRadiusWhite;
    G4double fInnerRadiusWhite[2];
    G4double fOuterRadiusWhite[2];

    //for the subtraction box
    G4double fSubBoxX, fSubBoxY, fSubBoxZ;
    G4double fMoveCut;

    //for back polyhedra
    G4double fBackRadius;
    G4double fInnerRadiusBack[2];
    G4double fOuterRadiusBack[2];
    G4double fZPlanesBack[2];

    // for subtraction cylinder
    G4double fMinRadius;
    G4double fMaxRadius;
    G4double fHalfLength;
    G4double fSubtractionStartPhi;
    G4double fSubtractionDeltaPhi;

    // The Euler angles from James' MSc thesis which gives us the detector positions
    // Some of the angles for the green and yellow detectors are wrong in James' thesis,
    // note the +180 on a few angles.
    G4double fBlueAlphaBetaGamma[15][3];
    G4double fGreenAlphaBetaGamma[10][3];
    G4double fRedAlphaBetaGamma[15][3];
    G4double fWhiteAlphaBetaGamma[20][3];
    G4double fYellowAlphaBetaGamma[10][3];

    // Descant Structure Colour
    G4Colour fGreyColour;

    G4int BuildDescantShell();
};

#endif

