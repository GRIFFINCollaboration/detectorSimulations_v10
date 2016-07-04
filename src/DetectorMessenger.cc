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
// $Id: DetectorMessenger.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
    :fDetector(Det)
{
    fDetSysDir = new G4UIdirectory("/DetSys/");
    fDetSysDir->SetGuidance("UI commands of this example");

    fDetDir = new G4UIdirectory("/DetSys/det/");
    fDetDir->SetGuidance("detector control");

    fAppDir = new G4UIdirectory("/DetSys/app/");
    fAppDir->SetGuidance("apparatus control");

    fWorldDir = new G4UIdirectory("/DetSys/world/");
    fWorldDir->SetGuidance("world control");

    fWorldMaterialCmd = new G4UIcmdWithAString("/DetSys/world/material",this);
    fWorldMaterialCmd->SetGuidance("Select material for the world.");
    fWorldMaterialCmd->SetParameterName("choice",false);
    fWorldMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fWorldDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/world/dimensions",this);
    fWorldDimensionsCmd->SetGuidance("Set world dimensions - x y z unit.");
    fWorldDimensionsCmd->SetUnitCategory("Length");
    fWorldDimensionsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fWorldVisCmd = new G4UIcmdWithABool("/DetSys/world/vis",this);
    fWorldVisCmd->SetGuidance("Set the visulization of the world");
    fWorldVisCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fWorldMagneticFieldCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/world/magneticField",this);
    fWorldMagneticFieldCmd->SetGuidance("Set world magnetic field - x y z unit.");
    fWorldMagneticFieldCmd->SetUnitCategory("Magnetic flux density");
    fWorldMagneticFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fGenericTargetCmd = new G4UIcmdWithAString("/DetSys/app/genericTarget",this);
    fGenericTargetCmd->SetGuidance("Select material of the target.");
    fGenericTargetCmd->SetParameterName("choice",false);
    fGenericTargetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fGenericTargetDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/genericTargetDimensions",this);
    fGenericTargetDimensionsCmd->SetGuidance("Set target dimensions - x y z unit.");
    fGenericTargetDimensionsCmd->SetUnitCategory("Length");
    fGenericTargetDimensionsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fGenericTargetPositionCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/genericTargetPosition",this);
    fGenericTargetPositionCmd->SetGuidance("Set target position - x y z unit.");
    fGenericTargetPositionCmd->SetUnitCategory("Length");
    fGenericTargetPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fFieldBoxMaterialCmd = new G4UIcmdWithAString("/DetSys/app/fieldBoxMaterial",this);
    fFieldBoxMaterialCmd->SetGuidance("Select material of the target.");
    fFieldBoxMaterialCmd->SetParameterName("choice",false);
    fFieldBoxMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fFieldBoxDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/fieldBoxDimensions",this);
    fFieldBoxDimensionsCmd->SetGuidance("Set target dimensions - x y z unit.");
    fFieldBoxDimensionsCmd->SetUnitCategory("Length");
    fFieldBoxDimensionsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fFieldBoxPositionCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/fieldBoxPosition",this);
    fFieldBoxPositionCmd->SetGuidance("Set target position - x y z unit.");
    fFieldBoxPositionCmd->SetUnitCategory("Length");
    fFieldBoxPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fFieldBoxMagneticFieldCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/fieldBoxMagneticField",this);
    fFieldBoxMagneticFieldCmd->SetGuidance("Set magnetic field - x y z unit.");
    fFieldBoxMagneticFieldCmd->SetUnitCategory("Magnetic flux density");
    fFieldBoxMagneticFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Box Stuff
    fAddBoxMatCmd = new G4UIcmdWithAString("/DetSys/det/boxMat",this);
    fAddBoxMatCmd->SetGuidance("Set box material.");
    fAddBoxMatCmd->SetParameterName("choice",false);
    fAddBoxMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddBoxThicknessCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/det/boxThickness",this);
    fAddBoxThicknessCmd->SetGuidance("Set box thickness.");
    fAddBoxThicknessCmd->SetUnitCategory("Length");
    fAddBoxThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddBoxInnerDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/det/boxInnerDimensions",this);
    fAddBoxInnerDimensionsCmd->SetGuidance("Set box inner dimensions.");
    fAddBoxInnerDimensionsCmd->SetUnitCategory("Length");
    fAddBoxInnerDimensionsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddBoxColourCmd = new G4UIcmdWith3Vector("/DetSys/det/boxColour",this);
    fAddBoxColourCmd->SetGuidance("Set box colour.");
    fAddBoxColourCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddBoxCmd = new G4UIcmdWithoutParameter("/DetSys/det/addBox",this);
    fAddBoxCmd->SetGuidance("Add a box.");
    fAddBoxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Grid Stuff
    fAddGridMatCmd = new G4UIcmdWithAString("/DetSys/det/gridMat",this);
    fAddGridMatCmd->SetGuidance("Set grid material.");
    fAddGridMatCmd->SetParameterName("choice",false);
    fAddGridMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddGridSizeCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/det/gridSize",this);
    fAddGridSizeCmd->SetGuidance("Set grid size.");
    fAddGridSizeCmd->SetUnitCategory("Length");
    fAddGridSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddGridDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/det/gridDimensions",this);
    fAddGridDimensionsCmd->SetGuidance("Set grid dimensions.");
    fAddGridDimensionsCmd->SetUnitCategory("Length");
    fAddGridDimensionsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddGridColourCmd = new G4UIcmdWith3Vector("/DetSys/det/gridColour",this);
    fAddGridColourCmd->SetGuidance("Set grid colour.");
    fAddGridColourCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddGridCmd = new G4UIcmdWithoutParameter("/DetSys/det/addGrid",this);
    fAddGridCmd->SetGuidance("Add a grid.");
    fAddGridCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddGridPosOffsetCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/det/gridPosOffset",this);
    fAddGridPosOffsetCmd->SetGuidance("Set grid offset.");
    fAddGridPosOffsetCmd->SetUnitCategory("Length");
    fAddGridPosOffsetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddApparatusSpiceTargetChamberCmd = new G4UIcmdWithoutParameter("/DetSys/app/addSpiceTargetChamber",this);
    fAddApparatusSpiceTargetChamberCmd->SetGuidance("Add SPICE target chamber.");
    fAddApparatusSpiceTargetChamberCmd->AvailableForStates(G4State_Idle);

    fAddApparatus8piVacuumChamberCmd = new G4UIcmdWithoutParameter("/DetSys/app/add8piVacuumChamber",this);
    fAddApparatus8piVacuumChamberCmd->SetGuidance("Add 8pi vacuum chamber.");
    fAddApparatus8piVacuumChamberCmd->AvailableForStates(G4State_Idle);

    fAddApparatus8piVacuumChamberAuxMatShellCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/app/add8piVacuumChamberAuxMatShell",this);
    fAddApparatus8piVacuumChamberAuxMatShellCmd->SetGuidance("Add AuxMat shell around 8pi vacuum chamber");
    fAddApparatus8piVacuumChamberAuxMatShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddApparatusGriffinStructureCmd = new G4UIcmdWithAnInteger("/DetSys/app/addGriffinStructure",this);
    fAddApparatusGriffinStructureCmd->SetGuidance("Add Griffin Structure");
    fAddApparatusGriffinStructureCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fUpdateCmd = new G4UIcmdWithoutParameter("/DetSys/det/update",this);
    fUpdateCmd->SetGuidance("Update geometry.");
    fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
    fUpdateCmd->AvailableForStates(G4State_Idle);

    fAddDetectionSystemGammaTrackingCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGammaTracking",this);
    fAddDetectionSystemGammaTrackingCmd->SetGuidance("Add Detection System GammaTracking");
    fAddDetectionSystemGammaTrackingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystemSodiumIodideCmd = new G4UIcmdWithAnInteger("/DetSys/det/addSodiumIodide",this);
    fAddDetectionSystemSodiumIodideCmd->SetGuidance("Add Detection System SodiumIodide");
    fAddDetectionSystemSodiumIodideCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystemLanthanumBromideCmd = new G4UIcmdWith3Vector("/DetSys/det/addLanthanumBromide",this);
    fAddDetectionSystemLanthanumBromideCmd->SetGuidance("Add Detection System LanthanumBromide - number of dets, radius in cm, empty");
    fAddDetectionSystemLanthanumBromideCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystemAncillaryBGOCmd = new G4UIcmdWith3Vector("/DetSys/det/addAncillaryBGO",this);
    fAddDetectionSystemAncillaryBGOCmd->SetGuidance("Add Detection System AncillaryBGO - number of dets, radius in cm, hevimet yes(1) no(0)");
    fAddDetectionSystemAncillaryBGOCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystem8piCmd = new G4UIcmdWithAnInteger("/DetSys/det/add8pi",this);
    fAddDetectionSystem8piCmd->SetGuidance("Add Detection System 8pi");
    fAddDetectionSystem8piCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystem8piDetectorCmd = new G4UIcmdWithAnInteger("/DetSys/det/add8piDetector",this);
    fAddDetectionSystem8piDetectorCmd->SetGuidance("Add 8pi Detector");
    fAddDetectionSystem8piDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fAddDetectionSystemDescantCmd = new G4UIcmdWithAnInteger("/DetSys/det/addDescant",this);
    fAddDetectionSystemDescantCmd->SetGuidance("Add Detection System DESCANT");
    fAddDetectionSystemDescantCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    fAddDetectionSystemDescantAuxPortsCmd = new G4UIcmdWith3Vector("/DetSys/det/addDescantAuxPorts",this);
    fAddDetectionSystemDescantAuxPortsCmd->SetGuidance("Add 8 DESCANT detectors in the auxillary LaBr3 detector locations");
    fAddDetectionSystemDescantAuxPortsCmd->SetGuidance("/DetSys/det/addDescantAuxPorts _nDet_ _radialPos_ _leadShield_");
    fAddDetectionSystemDescantAuxPortsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddApparatusDescantStructureCmd = new G4UIcmdWithoutParameter("/DetSys/det/addDescantStructure",this);
    fAddApparatusDescantStructureCmd->SetGuidance("Add DESCANT structure");
    fAddApparatusDescantStructureCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddDetectionSystemTestcanCmd = new G4UIcmdWith3Vector("/DetSys/det/addTestcan",this);
    fAddDetectionSystemTestcanCmd->SetGuidance("Add Testcan Detection System");
    fAddDetectionSystemTestcanCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
	 fSetDetectionSystemDescantColorCmd = new G4UIcmdWithAString("/DetSys/det/setDescantColor", this);
	 fSetDetectionSystemDescantColorCmd->SetGuidance("Set color of next descant detector to be added via addDescantCart or addDescantSpher");
	 fSetDetectionSystemDescantColorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 fSetDetectionSystemDescantColorCmd->SetCandidates("white blue red green yellow");
	 
	 fSetDetectionSystemDescantRotationCmd = new G4UIcmdWith3Vector("/DetSys/det/setDescantRotation", this);
	 fSetDetectionSystemDescantRotationCmd->SetGuidance("Set rotation of next descant detector to be added via addDescantCart or addDescantSpher (alhpa beta gamma)");
	 fSetDetectionSystemDescantRotationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
	 fAddDetectionSystemDescantCartCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/det/addDescantCart", this);
	 fAddDetectionSystemDescantCartCmd->SetGuidance("Add single DESCANT detector at provided cartesian coordinates");
	 fAddDetectionSystemDescantCartCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	 
	 fAddDetectionSystemDescantSpherCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/det/addDescantSpher", this);
	 fAddDetectionSystemDescantSpherCmd->SetGuidance("Add single DESCANT detector at provided spherical coordinates");
	 fAddDetectionSystemDescantSpherCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	 
    fAddDetectionSystemGriffinForwardCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGriffinForward",this);
    fAddDetectionSystemGriffinForwardCmd->SetGuidance("Add Detection System GriffinForward");
    fAddDetectionSystemGriffinForwardCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystemGriffinForwardDetectorCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGriffinForwardDetector",this);
    fAddDetectionSystemGriffinForwardDetectorCmd->SetGuidance("Add GriffinForward Detector");
    fAddDetectionSystemGriffinForwardDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystemGriffinBackCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGriffinBack",this);
    fAddDetectionSystemGriffinBackCmd->SetGuidance("Add Detection System GriffinBack");
    fAddDetectionSystemGriffinBackCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystemGriffinBackDetectorCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGriffinBackDetector",this);
    fAddDetectionSystemGriffinBackDetectorCmd->SetGuidance("Add GriffinBack Detector");
    fAddDetectionSystemGriffinBackDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAddDetectionSystemGriffinCustomDetectorCmd = new G4UIcmdWithAnInteger( "/DetSys/det/addGriffinCustomDetector", this ) ;
    fAddDetectionSystemGriffinCustomDetectorCmd->SetGuidance( "Adds a detector using the paramaters specified") ;
    fAddDetectionSystemGriffinCustomDetectorCmd->AvailableForStates( G4State_PreInit, G4State_Idle ) ;
	 
    fAddDetectionSystemGriffinCustomCmd = new G4UIcmdWithAnInteger( "/DetSys/det/addGriffinCustom", this ) ;
    fAddDetectionSystemGriffinCustomCmd->SetGuidance( "Adds a detection system using the paramaters specified" ) ;
    fAddDetectionSystemGriffinCustomCmd->AvailableForStates( G4State_PreInit, G4State_Idle ) ;
	 
    //////// Commands that are required for addGriffinCustom
    fAddDetectionSystemGriffinShieldSelectCmd = new G4UIcmdWithAnInteger( "/DetSys/det/SetCustomShieldsPresent", this );
    fAddDetectionSystemGriffinShieldSelectCmd->SetGuidance( "Selects whether or not the detector suppressors are included" ) ;
    fAddDetectionSystemGriffinShieldSelectCmd->AvailableForStates( G4State_PreInit, G4State_Idle ) ;
	 
    fAddDetectionSystemGriffinSetRadialDistanceCmd = new G4UIcmdWithADoubleAndUnit( "/DetSys/det/SetCustomRadialDistance", this ) ;
    fAddDetectionSystemGriffinSetRadialDistanceCmd->SetGuidance( "Selects the radial distance for the detector from the origin") ;
    fAddDetectionSystemGriffinSetRadialDistanceCmd->AvailableForStates( G4State_PreInit, G4State_Idle ) ;
	 
    fAddDetectionSystemGriffinSetExtensionSuppLocationCmd = new G4UIcmdWithAnInteger( "/DetSys/det/SetCustomExtensionSuppressorLocation", this ) ;
    fAddDetectionSystemGriffinSetExtensionSuppLocationCmd->SetGuidance( "Selects a position for the extension suppressors. Either forward (0) or back (1)." ) ;
    fAddDetectionSystemGriffinSetExtensionSuppLocationCmd->AvailableForStates( G4State_PreInit, G4State_Idle ) ;
	 
    fAddDetectionSystemGriffinSetDeadLayerCmd = new G4UIcmdWith3Vector( "/DetSys/det/SetCustomDeadLayer", this ) ;
    fAddDetectionSystemGriffinSetDeadLayerCmd->SetGuidance( "Sets the dead layer and position for the detector placed in the next call to addGriffinCustom." ) ;
    fAddDetectionSystemGriffinSetDeadLayerCmd->AvailableForStates( G4State_PreInit, G4State_Idle ) ;
    ////////
	 
    fAddDetectionSystemGriffinHevimetCmd = new G4UIcmdWithAnInteger( "/DetSys/det/includeGriffinHevimet", this ) ;
    fAddDetectionSystemGriffinHevimetCmd->SetGuidance( "Includes the Hevimet for a Griffin detector." ) ;
    fAddDetectionSystemGriffinHevimetCmd->AvailableForStates( G4State_PreInit, G4State_Idle ) ;
	 
    fAddDetectionSystemSceptarCmd = new G4UIcmdWithAnInteger("/DetSys/det/addSceptar",this);
    fAddDetectionSystemSceptarCmd->SetGuidance("Add Detection System Sceptar");
    fAddDetectionSystemSceptarCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystemSpiceCmd = new G4UIcmdWithAnInteger("/DetSys/det/addSpice",this);
    fAddDetectionSystemSpiceCmd->SetGuidance("Add Detection System Spice");
    fAddDetectionSystemSpiceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystemSpiceV02Cmd = new G4UIcmdWithAnInteger("/DetSys/det/addSpiceV02",this);
    fAddDetectionSystemSpiceV02Cmd->SetGuidance("Add Detection System SpiceV02");
    fAddDetectionSystemSpiceV02Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fAddDetectionSystemPacesCmd = new G4UIcmdWithAnInteger("/DetSys/det/addPaces",this);
    fAddDetectionSystemPacesCmd->SetGuidance("Add Detection System Paces");
    fAddDetectionSystemPacesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	 
    fUseTIGRESSPositionsCmd = new G4UIcmdWithABool("/DetSys/det/UseTIGRESSPositions",this);
    fUseTIGRESSPositionsCmd->SetGuidance("Use TIGRESS detector positions rather than GRIFFIN");
    fUseTIGRESSPositionsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
    delete fDetDir;
    delete fAppDir;
    delete fWorldDir;
    delete fWorldMaterialCmd;
    delete fWorldDimensionsCmd;
    delete fWorldVisCmd;
    delete fWorldMagneticFieldCmd;
    delete fDetSysDir;
    delete fUpdateCmd;
    delete fGenericTargetCmd;
    delete fGenericTargetDimensionsCmd;
    delete fGenericTargetPositionCmd;
    delete fFieldBoxMaterialCmd;
    delete fFieldBoxDimensionsCmd;
    delete fFieldBoxPositionCmd;
    delete fFieldBoxMagneticFieldCmd;
    delete fAddBoxMatCmd;
    delete fAddBoxThicknessCmd;
    delete fAddBoxInnerDimensionsCmd;
    delete fAddBoxColourCmd;
    delete fAddBoxCmd;
    delete fAddGridMatCmd;
    delete fAddGridSizeCmd;
    delete fAddGridDimensionsCmd;
    delete fAddGridColourCmd;
    delete fAddGridPosOffsetCmd;
    delete fAddGridCmd;
    delete fAddApparatusSpiceTargetChamberCmd;
    delete fAddDetectionSystemGammaTrackingCmd;
    delete fAddApparatus8piVacuumChamberCmd;
    delete fAddApparatus8piVacuumChamberAuxMatShellCmd;
    delete fAddApparatusGriffinStructureCmd;


    delete fAddDetectionSystemSodiumIodideCmd;
    delete fAddDetectionSystemLanthanumBromideCmd;
    delete fAddDetectionSystemAncillaryBGOCmd;

    delete fAddDetectionSystem8piCmd;
    delete fAddDetectionSystem8piDetectorCmd;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    delete fAddDetectionSystemDescantCmd;
    delete fAddDetectionSystemDescantAuxPortsCmd;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    delete fAddApparatusDescantStructureCmd;

    delete fAddDetectionSystemTestcanCmd;

	 delete fSetDetectionSystemDescantColorCmd;
	 delete fSetDetectionSystemDescantRotationCmd;	 
	 delete fAddDetectionSystemDescantCartCmd;
	 delete fAddDetectionSystemDescantSpherCmd;

    delete fAddDetectionSystemSceptarCmd;
    delete fAddDetectionSystemSpiceCmd;
    delete fAddDetectionSystemSpiceV02Cmd;
    delete fAddDetectionSystemPacesCmd;

    delete fAddDetectionSystemGriffinForwardCmd;
    delete fAddDetectionSystemGriffinForwardDetectorCmd;
    delete fAddDetectionSystemGriffinBackCmd;
    delete fAddDetectionSystemGriffinBackDetectorCmd;
    delete fAddDetectionSystemGriffinCustomDetectorCmd ;
    delete fAddDetectionSystemGriffinCustomCmd ;
    delete fAddDetectionSystemGriffinHevimetCmd ;

    delete fAddDetectionSystemGriffinShieldSelectCmd ;
    delete fAddDetectionSystemGriffinSetRadialDistanceCmd ;
    delete fAddDetectionSystemGriffinSetExtensionSuppLocationCmd ;
    delete fAddDetectionSystemGriffinSetDeadLayerCmd ;

    delete fUseTIGRESSPositionsCmd;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
    if(command == fWorldMaterialCmd ) {
        fDetector->SetWorldMaterial(newValue);
		  return;
    }
    if(command == fWorldDimensionsCmd ) {
        fDetector->SetWorldDimensions(fWorldDimensionsCmd->GetNew3VectorValue(newValue));
    }
    if(command == fWorldVisCmd ) {
        fDetector->SetWorldVis(fWorldVisCmd->GetNewBoolValue(newValue));
    }
    //  if(command == fWorldMagneticFieldCmd ) {
    //    fDetector->SetWorldMagneticField(fWorldMagneticFieldCmd->GetNew3VectorValue(newValue));
    //  }
    if(command == fUpdateCmd ) {
        fDetector->UpdateGeometry();
    }
    //  if(command == fGenericTargetCmd ) {
    //    fDetector->SetGenericTargetMaterial(newValue);
    //  }
    //  if(command == fGenericTargetDimensionsCmd ) {
    //    fDetector->SetGenericTargetDimensions(fGenericTargetDimensionsCmd->GetNew3VectorValue(newValue));
    //  }
    //  if(command == fGenericTargetPositionCmd ) {
    //    fDetector->SetGenericTargetPosition(fGenericTargetPositionCmd->GetNew3VectorValue(newValue));
    //  }
    //  if(command == fFieldBoxMaterialCmd ) {
    //    fDetector->SetFieldBoxMaterial(newValue);
    //  }
    //  if(command == fFieldBoxDimensionsCmd ) {
    //    fDetector->SetFieldBoxDimensions(fFieldBoxDimensionsCmd->GetNew3VectorValue(newValue));
    //  }
    //  if(command == fFieldBoxPositionCmd ) {
    //    fDetector->SetFieldBoxPosition(fFieldBoxPositionCmd->GetNew3VectorValue(newValue));
    //  }
    //  if(command == fFieldBoxMagneticFieldCmd ) {
    //    fDetector->SetFieldBoxMagneticField(fFieldBoxMagneticFieldCmd->GetNew3VectorValue(newValue));
    //  }
    //  if(command == fAddBoxMatCmd ) {
    //    fDetector->SetBoxMat(newValue);
    //  }
    //  if(command == fAddBoxThicknessCmd ) {
    //    fDetector->SetBoxThickness(fAddBoxThicknessCmd->GetNewDoubleValue(newValue));
    //  }
    //  if(command == fAddBoxInnerDimensionsCmd ) {
    //    fDetector->SetBoxInnerDimensions(fAddBoxInnerDimensionsCmd->GetNew3VectorValue(newValue));
    //  }
    //  if(command == fAddBoxColourCmd ) {
    //    fDetector->SetBoxColour(fAddBoxColourCmd->GetNew3VectorValue(newValue));
    //  }
    //  if(command == fAddBoxCmd ) {
    //    fDetector->AddBox();
    //  }
    if(command == fAddGridMatCmd ) {
        fDetector->SetGridMat(newValue);
    }
    if(command == fAddGridSizeCmd ) {
        fDetector->SetGridSize(fAddGridSizeCmd->GetNewDoubleValue(newValue));
    }
    if(command == fAddGridDimensionsCmd ) {
        fDetector->SetGridDimensions(fAddGridDimensionsCmd->GetNew3VectorValue(newValue));
    }
    if(command == fAddGridColourCmd ) {
        fDetector->SetGridColour(fAddGridColourCmd->GetNew3VectorValue(newValue));
    }
    if(command == fAddGridPosOffsetCmd ) {
        fDetector->SetGridPosOffset(fAddGridPosOffsetCmd->GetNew3VectorValue(newValue));
    }
    if(command == fAddGridCmd ) {
        fDetector->AddGrid();
    }
    //  if(command == fAddApparatusSpiceTargetChamberCmd ) {
    //    fDetector->fAddApparatusSpiceTargetChamber();
    //  }
    if(command == fAddApparatus8piVacuumChamberCmd ) {
        fDetector->AddApparatus8piVacuumChamber();
    }
    if(command == fAddApparatus8piVacuumChamberAuxMatShellCmd ) {
        fDetector->AddApparatus8piVacuumChamberAuxMatShell(fAddApparatus8piVacuumChamberAuxMatShellCmd->GetNewDoubleValue(newValue));
    }
    if(command == fAddApparatusGriffinStructureCmd ) {
        fDetector->AddApparatusGriffinStructure(fAddApparatusGriffinStructureCmd->GetNewIntValue(newValue));
    }
    //  if(command == fAddDetectionSystemGammaTrackingCmd ) {
    //    fDetector->AddDetectionSystemGammaTracking(fAddDetectionSystemGammaTrackingCmd->GetNewIntValue(newValue));
    //  }
    //  if(command == fAddDetectionSystemSodiumIodideCmd ) {
    //    fDetector->AddDetectionSystemSodiumIodide(fAddDetectionSystemSodiumIodideCmd->GetNewIntValue(newValue));
    //  }
    if(command == fAddDetectionSystemLanthanumBromideCmd ) {
        fDetector->AddDetectionSystemLanthanumBromide(fAddDetectionSystemLanthanumBromideCmd->GetNew3VectorValue(newValue));
    }
    if(command == fAddDetectionSystemAncillaryBGOCmd ) {
        fDetector->AddDetectionSystemAncillaryBGO(fAddDetectionSystemAncillaryBGOCmd->GetNew3VectorValue(newValue));
    }
    if(command == fAddDetectionSystem8piCmd ) {
        fDetector->AddDetectionSystem8pi(fAddDetectionSystem8piCmd->GetNewIntValue(newValue));
    }
    if(command == fAddDetectionSystem8piDetectorCmd ) {
        fDetector->AddDetectionSystem8piDetector(fAddDetectionSystem8piDetectorCmd->GetNewIntValue(newValue));
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(command == fAddDetectionSystemDescantCmd ) {
        fDetector->AddDetectionSystemDescant(fAddDetectionSystemDescantCmd->GetNewIntValue(newValue));
    }
    if(command == fAddDetectionSystemDescantAuxPortsCmd )  {
        fDetector->AddDetectionSystemDescantAuxPorts(fAddDetectionSystemDescantAuxPortsCmd->GetNew3VectorValue(newValue));
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(command == fAddApparatusDescantStructureCmd ) {
        fDetector->AddApparatusDescantStructure();
    }

	 if(command == fSetDetectionSystemDescantColorCmd ) {
		fDetector->SetDetectionSystemDescantColor(newValue);
	 }
	 if(command == fSetDetectionSystemDescantRotationCmd ) {
		fDetector->SetDetectionSystemDescantRotation(fSetDetectionSystemDescantRotationCmd->GetNew3VectorValue(newValue));
	 }
	 if(command == fAddDetectionSystemDescantCartCmd ) {
		fDetector->AddDetectionSystemDescantCart(fAddDetectionSystemDescantCartCmd->GetNew3VectorValue(newValue));
	 }
	 if(command == fAddDetectionSystemDescantSpherCmd ) {
		fDetector->AddDetectionSystemDescantSpher(fAddDetectionSystemDescantSpherCmd->GetNew3VectorRawValue(newValue), fAddDetectionSystemDescantSpherCmd->GetNewUnitValue(newValue));
	 }

    if(command == fAddDetectionSystemTestcanCmd ) { 
        fDetector->AddDetectionSystemTestcan(fAddDetectionSystemTestcanCmd->GetNew3VectorValue(newValue));
    }

    if(command == fAddDetectionSystemSceptarCmd ) {
        fDetector->AddDetectionSystemSceptar(fAddDetectionSystemSceptarCmd->GetNewIntValue(newValue));
    }
    if(command == fAddDetectionSystemGriffinForwardCmd ) {
        fDetector->AddDetectionSystemGriffinForward(fAddDetectionSystemGriffinForwardCmd->GetNewIntValue(newValue));
    }
    if(command == fAddDetectionSystemGriffinForwardDetectorCmd ) {
        fDetector->AddDetectionSystemGriffinForwardDetector(fAddDetectionSystemGriffinForwardDetectorCmd->GetNewIntValue(newValue));
    }
    if(command == fAddDetectionSystemGriffinBackCmd ) {
        fDetector->AddDetectionSystemGriffinBack(fAddDetectionSystemGriffinBackCmd->GetNewIntValue(newValue));
    }
    if(command == fAddDetectionSystemGriffinBackDetectorCmd ) {
        fDetector->AddDetectionSystemGriffinBackDetector(fAddDetectionSystemGriffinBackDetectorCmd->GetNewIntValue(newValue));
    }
    //  if(command == fAddDetectionSystemGriffinPositionConfigCmd ) {
    //    fDetector->AddDetectionSystemGriffinPositionConfig(fAddDetectionSystemGriffinPositionConfigCmd->GetNew3VectorValue(newValue));
    //	}
    if(command == fAddDetectionSystemGriffinCustomDetectorCmd ) {
        fDetector->AddDetectionSystemGriffinCustomDetector(fAddDetectionSystemGriffinCustomDetectorCmd->GetNewIntValue(newValue));
    }
    if(command == fAddDetectionSystemGriffinCustomCmd ) {
        fDetector->AddDetectionSystemGriffinCustom(fAddDetectionSystemGriffinCustomCmd->GetNewIntValue(newValue));
    }
    if(command == fAddDetectionSystemGriffinHevimetCmd ) {
        fDetector->AddDetectionSystemGriffinHevimet(fAddDetectionSystemGriffinHevimetCmd->GetNewIntValue(newValue));
    }
    if(command == fAddDetectionSystemGriffinShieldSelectCmd ) {
        fDetector->AddDetectionSystemGriffinShieldSelect(fAddDetectionSystemGriffinShieldSelectCmd->GetNewIntValue(newValue)) ;
    }
    if(command == fAddDetectionSystemGriffinSetRadialDistanceCmd ) {
        fDetector->AddDetectionSystemGriffinSetRadialDistance(fAddDetectionSystemGriffinSetRadialDistanceCmd->GetNewDoubleValue(newValue)) ;
    }
    if(command == fAddDetectionSystemGriffinSetExtensionSuppLocationCmd ) {
        fDetector->AddDetectionSystemGriffinSetExtensionSuppLocation(fAddDetectionSystemGriffinSetExtensionSuppLocationCmd->GetNewIntValue(newValue)) ;
    }
    if(command == fAddDetectionSystemGriffinSetDeadLayerCmd ) {
        fDetector->AddDetectionSystemGriffinSetDeadLayer(fAddDetectionSystemGriffinSetDeadLayerCmd->GetNew3VectorValue(newValue ) ) ;
    }
    
    //  if(command == fAddDetectionSystemSpiceCmd ) {
    //    fDetector->AddDetectionSystemSpice(fAddDetectionSystemSpiceCmd->GetNewIntValue(newValue));
    //  }
    //  if(command == fAddDetectionSystemSpiceV02Cmd ) {
    //    fDetector->AddDetectionSystemSpiceV02(fAddDetectionSystemSpiceV02Cmd->GetNewIntValue(newValue));
    //  }
    
    if(command == fAddDetectionSystemPacesCmd ) {
        fDetector->AddDetectionSystemPaces(fAddDetectionSystemPacesCmd->GetNewIntValue(newValue));
    }
    if(command == fUseTIGRESSPositionsCmd ) {
        fDetector->UseTIGRESSPositions(fUseTIGRESSPositionsCmd->GetNewBoolValue(newValue));
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
