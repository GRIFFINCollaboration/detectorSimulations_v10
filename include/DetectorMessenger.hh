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
// $Id: DetectorMessenger.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
public:
    DetectorMessenger( DetectorConstruction* ) ;
    ~DetectorMessenger();

    void SetNewValue(G4UIcommand*, G4String);

private:
    DetectorConstruction* Detector;

    G4UIdirectory*             detDir;
    G4UIdirectory*             appDir;
    G4UIdirectory*             worldDir;
    G4UIdirectory*             DetSysDir;
    G4UIcmdWithAString*        WorldMaterialCmd;
    G4UIcmdWith3VectorAndUnit* WorldDimensionsCmd;
    G4UIcmdWithABool*          WorldVisCmd;
    G4UIcmdWith3VectorAndUnit* WorldMagneticFieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;

    // Generic Target Apparatus
    G4UIcmdWithAString*        GenericTargetCmd;
    G4UIcmdWith3VectorAndUnit* GenericTargetDimensionsCmd;
    G4UIcmdWith3VectorAndUnit* GenericTargetPositionCmd;

    G4UIcmdWithAString*        FieldBoxMaterialCmd;
    G4UIcmdWith3VectorAndUnit* FieldBoxDimensionsCmd;
    G4UIcmdWith3VectorAndUnit* FieldBoxPositionCmd;
    G4UIcmdWith3VectorAndUnit* FieldBoxMagneticFieldCmd;

    G4UIcmdWithoutParameter*   AddApparatusSpiceTargetChamberCmd;
    G4UIcmdWithoutParameter*   AddApparatus8piVacuumChamberCmd;
    G4UIcmdWithADoubleAndUnit*      AddApparatus8piVacuumChamberAuxMatShellCmd;
    G4UIcmdWithAnInteger*   AddApparatusGriffinStructureCmd;

    G4UIcmdWithAString*         addBoxMatCmd;
    G4UIcmdWithADoubleAndUnit*  addBoxThicknessCmd;
    G4UIcmdWith3VectorAndUnit*  addBoxInnerDimensionsCmd;
    G4UIcmdWith3Vector*         addBoxColourCmd;
    G4UIcmdWithoutParameter*    addBoxCmd;

    G4UIcmdWithAString*         addGridMatCmd;
    G4UIcmdWithADoubleAndUnit*  addGridSizeCmd;
    G4UIcmdWith3VectorAndUnit*  addGridDimensionsCmd;
    G4UIcmdWith3Vector*         addGridColourCmd;
    G4UIcmdWith3VectorAndUnit*  addGridPosOffsetCmd;
    G4UIcmdWithoutParameter*    addGridCmd;


    // Detection Systems
    G4UIcmdWithAnInteger*       AddDetectionSystemGammaTrackingCmd;
    G4UIcmdWithAnInteger*       AddDetectionSystemSodiumIodideCmd;
    G4UIcmdWith3Vector*         AddDetectionSystemLanthanumBromideCmd;
    G4UIcmdWith3Vector*         AddDetectionSystemAncillaryBGOCmd;

    G4UIcmdWithAnInteger*       AddDetectionSystem8piCmd;
    G4UIcmdWithAnInteger*       AddDetectionSystem8piDetectorCmd;

    G4UIcmdWithAnInteger*       AddDetectionSystemDescantCmd;
    G4UIcmdWith3Vector*         AddDetectionSystemDescantAuxPortsCmd;
    
    G4UIcmdWithoutParameter*	  AddApparatusDescantStructureCmd;
    G4UIcmdWithAString*         SetDetectionSystemDescantColorCmd;
    G4UIcmdWith3Vector*         SetDetectionSystemDescantRotationCmd;
    G4UIcmdWith3VectorAndUnit*  AddDetectionSystemDescantCartCmd;
    G4UIcmdWith3VectorAndUnit*  AddDetectionSystemDescantSpherCmd;

    G4UIcmdWithAnInteger*       AddDetectionSystemSceptarCmd;
    G4UIcmdWithAnInteger*       AddDetectionSystemGriffinForwardCmd;
    G4UIcmdWithAnInteger*       AddDetectionSystemGriffinForwardDetectorCmd;
    G4UIcmdWithAnInteger*       AddDetectionSystemGriffinBackCmd;
    G4UIcmdWithAnInteger*       AddDetectionSystemGriffinBackDetectorCmd;
    //G4UIcmdWith3Vector*       AddDetectionSystemGriffinPositionConfigCmd;
    G4UIcmdWithAnInteger*       AddDetectionSystemSpiceCmd;
    G4UIcmdWithAnInteger*       AddDetectionSystemSpiceV02Cmd;
    G4UIcmdWithAnInteger*       AddDetectionSystemPacesCmd;
    G4UIcmdWithAnInteger*       AddDetectionSystemGriffinHevimetCmd ;

    G4UIcmdWithAnInteger*       AddDetectionSystemGriffinCustomDetectorCmd ;
    G4UIcmdWithAnInteger*	     AddDetectionSystemGriffinCustomCmd ;
    G4UIcmdWithAnInteger*       AddDetectionSystemGriffinShieldSelectCmd ;
    G4UIcmdWithADoubleAndUnit*  AddDetectionSystemGriffinSetRadialDistanceCmd ;
    G4UIcmdWithAnInteger*       AddDetectionSystemGriffinSetExtensionSuppLocationCmd ;
    G4UIcmdWith3Vector*         AddDetectionSystemGriffinSetDeadLayerCmd ;
    G4UIcmdWithABool*           UseTIGRESSPositionsCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

