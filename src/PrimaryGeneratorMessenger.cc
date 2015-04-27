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
    :Action(Gun)
{
    numberOfDecayingLaBrDetectorsCmd = new G4UIcmdWithAnInteger("/DetSys/gun/numberOfDecayingLaBrDetectors",this);
    numberOfDecayingLaBrDetectorsCmd->SetGuidance("Set the number of radioactive LaBr detectors");
    numberOfDecayingLaBrDetectorsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    efficiencyEnergyCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/efficiencyEnergy",this);
    efficiencyEnergyCmd->SetGuidance("Set gamma efficiency energy");
    efficiencyEnergyCmd->SetUnitCategory("Energy");
    efficiencyEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    efficiencyDirectionCmd = new G4UIcmdWith3Vector("/DetSys/gun/direction",this);
    efficiencyDirectionCmd->SetGuidance("Set momentum direction.");
    efficiencyDirectionCmd->SetGuidance("Direction needs not to be a unit vector.");
    efficiencyDirectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    efficiencyPositionCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/gun/position",this);
    efficiencyPositionCmd->SetGuidance("Set position.");
    efficiencyPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    efficiencyParticleCmd = new G4UIcmdWithAString("/DetSys/gun/particle",this);
    efficiencyParticleCmd->SetGuidance("Set particle.");
    efficiencyParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    efficiencyPolarizationCmd = new G4UIcmdWith3Vector("/DetSys/gun/polarization",this);
    efficiencyPolarizationCmd->SetGuidance("Set gamma polarization direction.");
    efficiencyPolarizationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    efficiencyBeamRadiusCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/beamRadius",this);
    efficiencyBeamRadiusCmd->SetGuidance("Set beam radius");
    efficiencyBeamRadiusCmd->SetUnitCategory("Length");
    efficiencyBeamRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
    delete numberOfDecayingLaBrDetectorsCmd;
    delete efficiencyEnergyCmd;
    delete efficiencyDirectionCmd;
    delete efficiencyPolarizationCmd;
    delete efficiencyBeamRadiusCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue( G4UIcommand* command, G4String newValue )
{
    if( command == numberOfDecayingLaBrDetectorsCmd ) {
        Action->SetNumberOfDecayingLaBrDetectors(numberOfDecayingLaBrDetectorsCmd->GetNewIntValue(newValue));
    }
    if( command == efficiencyEnergyCmd ) {
        Action->SetEfficiencyEnergy(efficiencyEnergyCmd->GetNewDoubleValue(newValue));
    }
    if( command == efficiencyDirectionCmd ) {
        Action->SetEfficiencyDirection(efficiencyDirectionCmd->GetNew3VectorValue(newValue));
    }
    if( command == efficiencyPositionCmd ) {
        Action->SetEfficiencyPosition(efficiencyPositionCmd->GetNew3VectorValue(newValue));
    }
    if( command == efficiencyParticleCmd ) {
        Action->SetEfficiencyParticle(newValue);
    }
    if( command == efficiencyPolarizationCmd ) {
        Action->SetEfficiencyPolarization(efficiencyPolarizationCmd->GetNew3VectorValue(newValue));
    }
    if( command == efficiencyBeamRadiusCmd ) {
        Action->SetEfficiencyBeamRadius(efficiencyBeamRadiusCmd->GetNewDoubleValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

