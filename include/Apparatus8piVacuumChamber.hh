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

#ifndef Apparatus8piVacuumChamber_h
#define Apparatus8piVacuumChamber_h 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

///////////////////////////////////////////////////////////////////////
// Apparatus8piVacuumChamber
///////////////////////////////////////////////////////////////////////
class Apparatus8piVacuumChamber
{
public:
    Apparatus8piVacuumChamber();
    ~Apparatus8piVacuumChamber();

    void Build(G4LogicalVolume*);

private:
    G4LogicalVolume* expHallLog;

    // LogicalVolumes used in Apparatus8piVacuumChamber()
    G4LogicalVolume* vacuum_chamber_sphere_log;
    //    G4LogicalVolume* vacuum_chamber_sphere_vacuum_log;
    
    // Physical Volumes used in Apparatus8piVacuumChamber()

    G4VPhysicalVolume* vacuum_chamber_sphere_phys;
    //    G4VPhysicalVolume* vacuum_chamber_sphere_vacuum_phys;

    ///////////////////////////////////////////////////////////////////
    // Apparatus8piVacuumChamberCylinder Properties
    ///////////////////////////////////////////////////////////////////

    // Materials
    G4String vacuum_material;
    G4String vacuum_chamber_sphere_material;

    // Dimensions
    G4double vacuum_chamber_inner_radius;
    G4double vacuum_chamber_outer_radius;

    // internal methods for ConstructApparatus8piVacuumChamber() for Building
    void BuildApparatus8piVacuumChamberSphere();

    // internal methods for ConstructApparatus8piVacuumChamber() for Placing
    void PlaceApparatus8piVacuumChamberSphere();
    
};

#endif
