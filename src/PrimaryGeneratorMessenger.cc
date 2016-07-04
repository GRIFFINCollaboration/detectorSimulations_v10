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
// $Id: PrimaryGeneratorMessenger.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun)
    :fAction(Gun)
{
    fNumberOfDecayingLaBrDetectorsCmd = new G4UIcmdWithAnInteger("/DetSys/gun/numberOfDecayingLaBrDetectors",this);
    fNumberOfDecayingLaBrDetectorsCmd->SetGuidance("Set the number of radioactive LaBr detectors");
    fNumberOfDecayingLaBrDetectorsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fEfficiencyEnergyCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/efficiencyEnergy",this);
    fEfficiencyEnergyCmd->SetGuidance("Set gamma efficiency energy");
    fEfficiencyEnergyCmd->SetUnitCategory("Energy");
    fEfficiencyEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fEfficiencyDirectionCmd = new G4UIcmdWith3Vector("/DetSys/gun/direction",this);
    fEfficiencyDirectionCmd->SetGuidance("Set momentum direction.");
    fEfficiencyDirectionCmd->SetGuidance("Direction needs not to be a unit vector.");
    fEfficiencyDirectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fEfficiencyPositionCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/gun/position",this);
    fEfficiencyPositionCmd->SetGuidance("Set position.");
    fEfficiencyPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fEfficiencyParticleCmd = new G4UIcmdWithAString("/DetSys/gun/particle",this);
    fEfficiencyParticleCmd->SetGuidance("Set particle.");
    fEfficiencyParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fEfficiencyPolarizationCmd = new G4UIcmdWith3Vector("/DetSys/gun/polarization",this);
    fEfficiencyPolarizationCmd->SetGuidance("Set gamma polarization direction.");
    fEfficiencyPolarizationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fEfficiencyBeamRadiusCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/beamRadius",this);
    fEfficiencyBeamRadiusCmd->SetGuidance("Set beam radius");
    fEfficiencyBeamRadiusCmd->SetUnitCategory("Length");
    fEfficiencyBeamRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {
    delete fNumberOfDecayingLaBrDetectorsCmd;
    delete fEfficiencyEnergyCmd;
    delete fEfficiencyDirectionCmd;
    delete fEfficiencyPolarizationCmd;
    delete fEfficiencyBeamRadiusCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if(command == fNumberOfDecayingLaBrDetectorsCmd) {
        fAction->SetNumberOfDecayingLaBrDetectors(fNumberOfDecayingLaBrDetectorsCmd->GetNewIntValue(newValue));
		  return;
    }
    if(command == fEfficiencyEnergyCmd ) {
        fAction->SetEfficiencyEnergy(fEfficiencyEnergyCmd->GetNewDoubleValue(newValue));
		  return;
    }
    if( command == fEfficiencyDirectionCmd ) {
        fAction->SetEfficiencyDirection(fEfficiencyDirectionCmd->GetNew3VectorValue(newValue));
		  return;
    }
    if( command == fEfficiencyPositionCmd ) {
        fAction->SetEfficiencyPosition(fEfficiencyPositionCmd->GetNew3VectorValue(newValue));
		  return;
    }
    if( command == fEfficiencyParticleCmd ) {
        fAction->SetEfficiencyParticle(newValue);
		  return;
    }
    if( command == fEfficiencyPolarizationCmd ) {
        fAction->SetEfficiencyPolarization(fEfficiencyPolarizationCmd->GetNew3VectorValue(newValue));
		  return;
    }
    if( command == fEfficiencyBeamRadiusCmd ) {
        fAction->SetEfficiencyBeamRadius(fEfficiencyBeamRadiusCmd->GetNewDoubleValue(newValue));
		  return;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

