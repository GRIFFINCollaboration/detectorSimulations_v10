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

#ifndef DetectionSystemGrid_h
#define DetectionSystemGrid_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;
class FieldSetup;

class DetectionSystemGrid
{
public:
    DetectionSystemGrid(G4double xLengthIn, G4double yLengthIn, G4double zLengthIn, G4double gridSizeIn, G4String gridMat, G4ThreeVector gridColour, G4ThreeVector posOffset);
    ~DetectionSystemGrid();

    //    G4int Build(G4SDManager* mySDman);
    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* expHallLog);

private:
    // Logical volumes
    G4LogicalVolume* fGridcellLog;

    // Assembly volumes
    G4AssemblyVolume* fAssembly;

    //    SensitiveDetector* crystalBlockSD;

    FieldSetup* fEmFieldSetup;

    G4double fXLength;
    G4double fYLength;
    G4double fZLength;

    G4double fXPosOffset;
    G4double fYPosOffset;
    G4double fZPosOffset;

    G4double fCubeSize;

    G4int fCellsXRow;
    G4int fCellsYRow;
    G4int fCellsZRow;


    G4int fCopyNumber;
    G4int fNumberOfCells;
    G4int fCellsPerRow;
    G4double fCellWidth;

    G4String fCellMaterial;

    G4ThreeVector fCellColour;

    G4Box* BuildCell();

    G4int BuildCellVolume();
};

#endif

