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
#include "DetectorConstruction.hh"//for SPICE target pedestal tunnelling

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

	fConeAngleCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/coneMaxAngle",this);//SPICE cone angle value
	fConeAngleCmd->SetGuidance("Set cone value for outer theta - use deg (0-90)");
	fConeAngleCmd->SetUnitCategory("Angle");
	fConeAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fConeMinAngleCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/coneMinAngle",this);//SPICE cone angle value
	fConeMinAngleCmd->SetGuidance("Set cone value for inner theta - use deg (0-90) - default is 0 if none specified");
	fConeMinAngleCmd->SetUnitCategory("Angle");
	fConeMinAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fBeamSpotSigmaCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/BeamSpot",this);//Beam spot sigma
	fBeamSpotSigmaCmd->SetGuidance("Set sigma for a realistic beamspot");
	fBeamSpotSigmaCmd->SetUnitCategory("Length");
	fBeamSpotSigmaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);   

	fBeamDistroCmd = new G4UIcmdWithAnInteger("/DetSys/gun/TargetLayer",this);//with target, can apply a distribution
	fBeamDistroCmd->SetGuidance("Set beam distribution within a target layer, zero indexed");
	fBeamDistroCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fBeamFileCmd = new G4UIcmdWithAString("/DetSys/gun/FileDistro",this);//with target, can apply a distribution
	fBeamFileCmd->SetGuidance("Set beam distribution within a target using definitions in a data file");
	fBeamFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fKentuckyEnergyCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/KentuckyEnergy",this);
	fKentuckyEnergyCmd->SetGuidance("Set neutron energy at zero-degree for Kentucky experiment.");
	fKentuckyEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fKentuckyReactionCmd = new G4UIcmdWithAString("/DetSys/gun/KentuckyReaction",this);
	fKentuckyReactionCmd->SetGuidance("Set reaction for Kentucky experiment (can be 'p,t','d,t', or 'd,d').");
	fKentuckyReactionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fMinimumPhiCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/MinimumPhi",this);
	fMinimumPhiCmd->SetGuidance("Set minimum phi for primary particle (has to be used before first call to KentuckyEnergy/KentuckyReaction!!!)");
	fMinimumPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fMaximumPhiCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/MaximumPhi",this);
	fMaximumPhiCmd->SetGuidance("Set maximum phi for primary particle (has to be used before first call to KentuckyEnergy/KentuckyReaction!!!)");
	fMaximumPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fMinimumThetaCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/MinimumTheta",this);
	fMinimumThetaCmd->SetGuidance("Set minimum theta for primary particle (has to be used before first call to KentuckyEnergy/KentuckyReaction!!!)");
	fMinimumThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fMaximumThetaCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/MaximumTheta",this);
	fMaximumThetaCmd->SetGuidance("Set maximum theta for primary particle (has to be used before first call to KentuckyEnergy/KentuckyReaction!!!)");
	fMaximumThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fVerbosityCmd = new G4UIcmdWithAnInteger("/DetSys/gun/verbose", this);
	fVerbosityCmd->SetGuidance("Set verbosity of PrimaryGeneratorAction. Call before setting Kentucky properties if this should be applied to that class as well.");
	fVerbosityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {
	delete fNumberOfDecayingLaBrDetectorsCmd;
	delete fEfficiencyEnergyCmd;
	delete fEfficiencyDirectionCmd;
	delete fEfficiencyPolarizationCmd;
	delete fConeAngleCmd;
	delete fConeMinAngleCmd;
	delete fBeamSpotSigmaCmd;
	delete fBeamDistroCmd;
	delete fBeamFileCmd;
	delete fKentuckyEnergyCmd;
	delete fKentuckyReactionCmd;
	delete fMinimumPhiCmd;
	delete fMaximumPhiCmd;
	delete fMinimumThetaCmd;
	delete fMaximumThetaCmd;
	delete fVerbosityCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
	if(command == fNumberOfDecayingLaBrDetectorsCmd) {
		fAction->SetNumberOfDecayingLaBrDetectors(fNumberOfDecayingLaBrDetectorsCmd->GetNewIntValue(newValue));
		return;
	}
	if(command == fEfficiencyEnergyCmd) {
		fAction->SetEfficiencyEnergy(fEfficiencyEnergyCmd->GetNewDoubleValue(newValue));
		return;
	}
	if(command == fEfficiencyDirectionCmd) {
		fAction->SetEfficiencyDirection(fEfficiencyDirectionCmd->GetNew3VectorValue(newValue));
		return;
	}
	if(command == fEfficiencyPositionCmd) {
		fAction->SetEfficiencyPosition(fEfficiencyPositionCmd->GetNew3VectorValue(newValue));
		return;
	}
	if(command == fEfficiencyParticleCmd) {
		fAction->SetEfficiencyParticle(newValue);
		return;
	}
	if(command == fEfficiencyPolarizationCmd) {
		fAction->SetEfficiencyPolarization(fEfficiencyPolarizationCmd->GetNew3VectorValue(newValue));
		return;
	}
	if(command == fConeAngleCmd) {
		fAction->SetConeMaxAngle(fConeAngleCmd->GetNewDoubleValue(newValue));
		G4cout<<"Cone Beam via Angle selected"<<G4endl;
		return;
	}
	if(command == fConeMinAngleCmd) {
		fAction->SetConeMinAngle(fConeMinAngleCmd->GetNewDoubleValue(newValue));
		G4cout<<"Cone Beam minimum angle supplied"<<G4endl;
		return;
	}
	if(command == fBeamSpotSigmaCmd) {
		fAction->SetBeamSpotSigma(fBeamSpotSigmaCmd->GetNewDoubleValue(newValue));
		G4cout<<"Beam Spot sigma supplied"<<G4endl;
		return;
	}  
	if(command == fBeamDistroCmd) {
		fAction->SetLayeredTargetBeamDistro(fBeamDistroCmd->GetNewIntValue(newValue));
		return;
	}
	if(command == fBeamFileCmd) {
		G4cout<<"Beam Distribution from file "<<newValue<<" selected "<< G4endl;
		fAction->PrepareBeamFile(newValue);
		return;
	}
	if(command == fKentuckyEnergyCmd) {
		fAction->SetKentuckyEnergy(fKentuckyEnergyCmd->GetNewDoubleValue(newValue));
		return;
	}
	if(command == fKentuckyReactionCmd) {
		fAction->SetKentuckyReaction(newValue);
		return;
	}
	if(command == fMinimumPhiCmd) {
		fAction->SetMinimumPhi(fMinimumPhiCmd->GetNewDoubleValue(newValue));
		return;
	}
	if(command == fMaximumPhiCmd) {
		fAction->SetMaximumPhi(fMaximumPhiCmd->GetNewDoubleValue(newValue));
		return;
	}
	if(command == fMinimumThetaCmd) {
		fAction->SetMinimumTheta(fMinimumThetaCmd->GetNewDoubleValue(newValue));
		return;
	}
	if(command == fMaximumThetaCmd) {
		fAction->SetMaximumTheta(fMaximumThetaCmd->GetNewDoubleValue(newValue));
		return;
	}
	if(command == fVerbosityCmd) {
		fAction->SetVerbosityLevel(fVerbosityCmd->GetNewIntValue(newValue));
		return;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
