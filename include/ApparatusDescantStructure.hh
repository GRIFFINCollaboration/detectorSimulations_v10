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
    G4LogicalVolume* descant_structure_log;

    // Assembly volumes
    G4AssemblyVolume* assemblyDescantStructure;

    G4String structure_material;
    G4double structure_inner_radius;
    G4double structure_outer_radius;
    G4double start_Phi;
    G4double delta_Phi;
    G4double start_Theta;
    G4double delta_Theta;
    G4double structure_radial_distance;

    //for subtraction polyhedra
    G4int numSides;
    G4int numZplanes;
    G4double cutExtra;
    G4double zPlanes[2];

    G4double polyRadiusRed;
    G4double rInnerRed[2];
    G4double rOuterRed[2];

    G4double polyRadiusGreenYellow;
    G4double rInnerGreenYellow[2];
    G4double rOuterGreenYellow[2];

    G4double polyRadiusWhite;
    G4double rInnerWhite[2];
    G4double rOuterWhite[2];

    //for the subtraction box
    G4double subBoxX, subBoxY, subBoxZ;
    G4double moveCut;

    //for back polyhedra
    G4double backRadius;
    G4double rInnerBack[2];
    G4double rOuterBack[2];
    G4double zPlanesBack[2];

    // for subtraction cylinder
    G4double rMin;
    G4double rMax;
    G4double halfLength;
    G4double sPhi;
    G4double dPhi;

    // The Euler angles from James' MSc thesis which gives us the detector positions
    // Some of the angles for the green and yellow detectors are wrong in James' thesis,
    // note the +180 on a few angles.
    G4double blue_alpha_beta_gamma[15][3];
    G4double green_alpha_beta_gamma[10][3];
    G4double red_alpha_beta_gamma[15][3];
    G4double white_alpha_beta_gamma[20][3];
    G4double yellow_alpha_beta_gamma[10][3];



    // Descant Structure Colour
    G4Colour grey_colour;

    G4int BuildDescantShell();


};

#endif

