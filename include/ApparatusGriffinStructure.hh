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
    G4int Place(G4LogicalVolume* exp_hall_log, G4int selector);

private:
    G4LogicalVolume* square_piece_log;
    G4LogicalVolume* triangle_piece_log;
    G4LogicalVolume* ring_piece_log;
    G4LogicalVolume* rod_piece_log;


    G4AssemblyVolume* assemblySquare;
    G4AssemblyVolume* assemblyTriangle;
    G4AssemblyVolume* assemblyRing;
    G4AssemblyVolume* assemblyRod;

private:
    // Materials
    G4String structure_material;
    G4String rod_material;

    // Dimensions
    G4double inner_distance_to_square_face; // inner radius to the centre of square plate
    G4double square_and_triangle_face_thickness; // thickness of the square and triangle plates
    G4double square_hole_length; // total legth of square hole, from inside structure.
    G4double square_hole_rounded_edge_radius; // the roundness of the square holes.
    G4double square_rod_hole_radius; // the radius of the holes for the steel rods
    G4double square_rod_hole_inner_radius; // the radius from the origin to the steel rod
    G4double square_rod_position_from_center_of_square_face; // the position the holes are from the origin. This is NOT a radius, this is along x or y
    G4double square_rod_total_length; // the total length of the rods
    G4double square_bgo_cut_length; // the length of the bgo slot
    G4double square_bgo_cut_width; // the width of the bgo slot
    G4double square_square_angle; // the angle between square faces in a rhombicuboctahedron
    G4double square_triangle_angle; // the angle between a square face and a triangle face in a rhombicuboctahedron
    G4double triangle_outer_hole_radius; // a hole radius in the triangle face. Radius to the rounded edge, this is the larger radius.
    G4double triangle_inner_hole_radius; // a hole radius in the triangle face. Radius to the square edge, this is the smaller radius.
    G4double large_ring_frame_inner_radius; // the inner radius of one of the the large rings that holds up the structure.
    G4double large_ring_frame_outer_radius; // the outer radius of one of the the large rings that holds up the structure.
    G4double large_ring_frame_thickness; // the thickness of one of the the large rings that holds up the structure.
    G4double large_ring_frame_z_position; // the positive z distance the ring is from the origin.

    G4double triangleThetaAngle;

    G4double griffinCoords[16][5];
    G4double ancillaryCoords[8][5];

    G4bool surfCheck;

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

    // General Spherical Translations
    G4double transX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi);
};

#endif
