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
    DetectorConstruction* fDetector;

    G4UIdirectory*             fDetDir;
    G4UIdirectory*             fAppDir;
    G4UIdirectory*             fWorldDir;
    G4UIdirectory*             fDetSysDir;
    G4UIcmdWithAString*        fWorldMaterialCmd;
    G4UIcmdWith3VectorAndUnit* fWorldDimensionsCmd;
    G4UIcmdWithABool*          fWorldVisCmd;
    G4UIcmdWith3VectorAndUnit* fWorldMagneticFieldCmd;
    G4UIcmdWithoutParameter*   fUpdateCmd;

    // Generic Target Apparatus
    G4UIcmdWithAString*        fGenericTargetCmd;
    G4UIcmdWith3VectorAndUnit* fGenericTargetDimensionsCmd;
    G4UIcmdWith3VectorAndUnit* fGenericTargetPositionCmd;

    G4UIcmdWithAString*        fFieldBoxMaterialCmd;
    G4UIcmdWith3VectorAndUnit* fFieldBoxDimensionsCmd;
    G4UIcmdWith3VectorAndUnit* fFieldBoxPositionCmd;
    G4UIcmdWith3VectorAndUnit* fFieldBoxMagneticFieldCmd;

    G4UIcmdWithoutParameter*   fAddApparatusSpiceTargetChamberCmd;
    G4UIcmdWithoutParameter*   fAddApparatus8piVacuumChamberCmd;
    G4UIcmdWithADoubleAndUnit* fAddApparatus8piVacuumChamberAuxMatShellCmd;
    G4UIcmdWithAnInteger*      fAddApparatusGriffinStructureCmd;

    G4UIcmdWithAString*         fAddBoxMatCmd;
    G4UIcmdWithADoubleAndUnit*  fAddBoxThicknessCmd;
    G4UIcmdWith3VectorAndUnit*  fAddBoxInnerDimensionsCmd;
    G4UIcmdWith3Vector*         fAddBoxColourCmd;
    G4UIcmdWithoutParameter*    fAddBoxCmd;

    G4UIcmdWithAString*         fAddGridMatCmd;
    G4UIcmdWithADoubleAndUnit*  fAddGridSizeCmd;
    G4UIcmdWith3VectorAndUnit*  fAddGridDimensionsCmd;
    G4UIcmdWith3Vector*         fAddGridColourCmd;
    G4UIcmdWith3VectorAndUnit*  fAddGridPosOffsetCmd;
    G4UIcmdWithoutParameter*    fAddGridCmd;


    // Detection Systems
    G4UIcmdWithAnInteger*       fAddDetectionSystemGammaTrackingCmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystemSodiumIodideCmd;
    G4UIcmdWith3Vector*         fAddDetectionSystemLanthanumBromideCmd;
    G4UIcmdWith3Vector*         fAddDetectionSystemAncillaryBGOCmd;

    G4UIcmdWithAnInteger*       fAddDetectionSystem8piCmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystem8piDetectorCmd;

    G4UIcmdWithAnInteger*       fAddDetectionSystemDescantCmd;
    G4UIcmdWith3Vector*         fAddDetectionSystemDescantAuxPortsCmd;
    
    G4UIcmdWithoutParameter*	  fAddApparatusDescantStructureCmd;
    G4UIcmdWithAString*         fSetDetectionSystemDescantColorCmd;
    G4UIcmdWith3Vector*         fSetDetectionSystemDescantRotationCmd;
    G4UIcmdWith3VectorAndUnit*  fAddDetectionSystemDescantCartCmd;
    G4UIcmdWith3VectorAndUnit*  fAddDetectionSystemDescantSpherCmd;

    G4UIcmdWith3Vector*         fAddDetectionSystemTestcanCmd;

    G4UIcmdWithAnInteger*       fAddDetectionSystemSceptarCmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystemGriffinForwardCmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystemGriffinForwardDetectorCmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystemGriffinBackCmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystemGriffinBackDetectorCmd;
    //G4UIcmdWith3Vector*       fAddDetectionSystemGriffinPositionConfigCmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystemSpiceCmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystemSpiceV02Cmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystemPacesCmd;
    G4UIcmdWithAnInteger*       fAddDetectionSystemGriffinHevimetCmd ;

    G4UIcmdWithAnInteger*       fAddDetectionSystemGriffinCustomDetectorCmd ;
    G4UIcmdWithAnInteger*	     fAddDetectionSystemGriffinCustomCmd ;
    G4UIcmdWithAnInteger*       fAddDetectionSystemGriffinShieldSelectCmd ;
    G4UIcmdWithADoubleAndUnit*  fAddDetectionSystemGriffinSetRadialDistanceCmd ;
    G4UIcmdWithAnInteger*       fAddDetectionSystemGriffinSetExtensionSuppLocationCmd ;
    G4UIcmdWith3Vector*         fAddDetectionSystemGriffinSetDeadLayerCmd ;
    G4UIcmdWithABool*           fUseTIGRESSPositionsCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

