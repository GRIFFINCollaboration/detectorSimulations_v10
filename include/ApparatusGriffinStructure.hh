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

#ifndef ApparatusGriffinStructure_h
#define ApparatusGriffinStructure_h 1

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units
#include "G4IntersectionSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

///////////////////////////////////////////////////////////////////////
// ApparatusGriffinStructure
///////////////////////////////////////////////////////////////////////
class ApparatusGriffinStructure
{
public:
    ApparatusGriffinStructure();
    ~ApparatusGriffinStructure();

    G4int Build();
    G4int Place(G4LogicalVolume* expHallLog, G4int selector);

private:
    G4LogicalVolume* fSquarePieceLog;
    G4LogicalVolume* fTrianglePieceLog;
    G4LogicalVolume* fRingPieceLog;
    G4LogicalVolume* fRodPieceLog;


    G4AssemblyVolume* fAssemblySquare;
    G4AssemblyVolume* fAssemblyTriangle;
    G4AssemblyVolume* fAssemblyRing;
    G4AssemblyVolume* fAssemblyRod;

private:
    // Materials
    G4String fStructureMaterial;
    G4String fRodMaterial;

    // Dimensions
    G4double fInnerDistanceToSquareFace; // inner radius to the centre of square plate
    G4double fSquareAndTriangleFaceThickness; // thickness of the square and triangle plates
    G4double fSquareHoleLength; // total legth of square hole, from inside structure.
    G4double fSquareHoleRoundedEdgeRadius; // the roundness of the square holes.
    G4double fSquareRodHoleRadius; // the radius of the holes for the steel rods
    G4double fSquareRodHoleInnerRadius; // the radius from the origin to the steel rod
    G4double fSquareRodPositionFromCenterOfSquareFace; // the position the holes are from the origin. This is NOT a radius, this is along x or y
    G4double fSquareRodTotalLength; // the total length of the rods
    G4double fSquareBgoCutLength; // the length of the bgo slot
    G4double fSquareBgoCutWidth; // the width of the bgo slot
    G4double fSquareSquareAngle; // the angle between square faces in a rhombicuboctahedron
    G4double fSquareTriangleAngle; // the angle between a square face and a triangle face in a rhombicuboctahedron
    G4double fTriangleOuterHoleRadius; // a hole radius in the triangle face. Radius to the rounded edge, this is the larger radius.
    G4double fTriangleInnerHoleRadius; // a hole radius in the triangle face. Radius to the square edge, this is the smaller radius.
    G4double fLargeRingFrameInnerRadius; // the inner radius of one of the the large rings that holds up the structure.
    G4double fLargeRingFrameOuterRadius; // the outer radius of one of the the large rings that holds up the structure.
    G4double fLargeRingFrameThickness; // the thickness of one of the the large rings that holds up the structure.
    G4double fLargeRingFrameZPosition; // the positive z distance the ring is from the origin.

    G4double fTriangleThetaAngle;

    G4double fGriffinCoords[16][5];
    G4double fAncillaryCoords[8][5];

    G4bool fSurfCheck;

    // Methods
    G4int BuildSquarePiece();
    G4int BuildTrianglePiece();
    G4int BuildLargeRingFrame();
    G4int BuildRodPiece();


    // Shapes
    G4SubtractionSolid* SquarePiece();
    G4SubtractionSolid* TrianglePiece();
    G4Tubs*             LargeRingFrame();
    G4Tubs*             RodPiece();


    // Hole Shapes
    G4SubtractionSolid* SquareWithRoundedCorners();
    G4SubtractionSolid* TruncatedThreeSidedCylinder();
};

#endif
