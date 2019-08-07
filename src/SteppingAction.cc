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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcess.hh"
#include "G4ProcessType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* detcon,
		EventAction* evt)
: G4UserSteppingAction(),
	fDetector(detcon), fEventAction(evt)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep) {
	G4bool trackSteps   = false;
	G4int processType   = -1;
	G4int evntNb        = 0;

	// Get volume of the current step
	G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	G4String volname = volume->GetName();
	// collect energy and track length step by step
	// As it's called more than once, get the Track and assign to variable
	G4double edep = aStep->GetTotalEnergyDeposit();
	G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();

	G4Track* theTrack = aStep->GetTrack();
	G4int stepNumber = theTrack->GetCurrentStepNumber();

	// Track particle type in EVERY step
	G4int particleType  = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

	const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
	G4int targetZ = -1;
	if(process != nullptr && process->GetProcessType() == fHadronic) {
		G4HadronicProcess* hadrProcess = (G4HadronicProcess*) process;
		const G4Isotope* target = nullptr;
		target = hadrProcess->GetTargetIsotope();
		if(target != nullptr) {
			targetZ = target->GetZ();
		}
	}

	// this can be modified to add more processes
	if(theTrack->GetCreatorProcess() != nullptr) {
		G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
		if(processName == "RadioactiveDecay")      processType = 1;
		else if(processName == "eIoni")            processType = 2;
		else if(processName == "msc")              processType = 3;
		else if(processName == "Scintillation")    processType = 4;
		else if(processName == "Cerenkov")         processType = 5;
		else if(processName == "hadElastic")         processType = 6;
		else if(processName == "neutronInelastic")         processType = 7;
		else if(processName == "nCapture")         processType = 8;
		else if(processName == "hIoni")         processType = 9;
		else if(processName == "ionIoni")         processType = 10;
		else if(processName == "compt")         processType = 11;
		else if(processName == "phot")         processType = 12;
		else if(processName == "conv")         processType = 13;
		else if(processName == "eBrem")         processType = 14;
		else if(processName == "annihil")         processType = 15;
		else if(processName == "dInelastic")         processType = 16;
		else if(processName == "CoulombScat")         processType = 17;
		else if(processName == "protonInelastic")         processType = 18;
		else if(processName == "alphaInelastic")         processType = 19;
		else if(processName == "tInelastic")         processType = 20;
		else if(processName == "photonNuclear")         processType = 21;
		else if(processName == "He3Inelastic")         processType = 22;
		else {                                      processType = 0;
		G4cout << "unknown process -> " << processName << G4endl;
		}
	}
	else {
		processType = -1;
	}

	evntNb =  fEventAction->GetEventNumber();

	// Get initial momentum direction & energy of particle
	G4int trackID = theTrack->GetTrackID();
	G4int parentID = theTrack->GetParentID();

	G4StepPoint* prePoint = aStep->GetPreStepPoint();
	G4StepPoint* postPoint = aStep->GetPostStepPoint();
	G4ThreeVector postPos = postPoint->GetPosition();
	G4double postTime = postPoint->GetGlobalTime();

	size_t found;

	//Set Lab Angle
	G4int nSecondaries = aStep->GetSecondary()->size();
	G4double lab_angle = -1;
	found = volname.find("PlasticDet");
<<<<<<< HEAD
	//G4cout << "Found " << found << G4endl;
	if(postPoint->GetProcessDefinedStep()->GetProcessName() == "hadElastic" && fEventAction->GetLabAngle() == -1 && found != G4String::npos) {
=======
	G4cout << "Found " << found << G4endl;
	if(postPoint->GetProcessDefinedStep()->GetProcessName() == "hadElastic" && nSecondaries == 1 && fEventAction->GetLabAngle() == -1 && found != G4String::npos) {
>>>>>>> 6828e1eec50f49505f6682d4e9b4096511be47f0
	G4ThreeVector momentum_1 = prePoint->GetMomentum();
	G4ThreeVector momentum_2 = postPoint->GetMomentum();
	lab_angle = momentum_2.angle(momentum_1);
	fEventAction->SetLabAngle(lab_angle);
	G4cout << "lab_angle: " << lab_angle << G4endl;
	G4cout << "GetLabAngle(): " << fEventAction->GetLabAngle() << G4endl;
	}

<<<<<<< HEAD
	//Get angle when leaving Detector
	G4double final_angle = -1;
	found = volname.find("PlasticDet");
	//G4cout << "Found " << found << G4endl;
	if(aStep->GetTrack()->GetParentID() == 0 && fEventAction->GetFinalAngle() == -1 && found != G4String::npos && aStep->IsFirstStepInVolume() == true) { 
	G4ThreeVector momentum_3 = prePoint->GetMomentum();
	//Set initial momenturm 
	fEventAction->SetInitialMomentum(momentum_3);
	G4cout<< "Setting Initial Momentum as : " << fEventAction->GetInitialMomentum() << G4endl;
	}
	if(aStep->GetTrack()->GetParentID() == 0 && fEventAction->GetFinalAngle() == -1 && found != G4String::npos && aStep->IsLastStepInVolume() == true) { 
	G4ThreeVector momentum_4 = postPoint->GetMomentum();
	//Set final momentum
	fEventAction->SetFinalMomentum(momentum_4);
	G4cout<< "Setting Final Momentum as : " << fEventAction->GetFinalMomentum() << G4endl;
	}
	G4ThreeVector check = G4ThreeVector(0., 0., 0.);
	if(aStep->GetTrack()->GetParentID() == 0 && fEventAction->GetFinalAngle() == -1 && found != G4String::npos && fEventAction->GetInitialMomentum() != check && fEventAction->GetFinalMomentum() != check) {
	G4ThreeVector initialM = fEventAction->GetInitialMomentum();
	G4ThreeVector finalM = fEventAction->GetFinalMomentum();
	final_angle = finalM.angle(initialM);
	fEventAction->SetFinalAngle(final_angle);
	G4cout<<"Calling Get Final Angle: " << fEventAction->GetFinalAngle() << G4endl;
	}
	final_angle= fEventAction->GetFinalAngle();
	G4cout<<"Calling Get Final Angle: " << fEventAction->GetFinalAngle() << G4endl;


	//Counting hits for efficiencies
	//By not initilalizing counters in event action to zero, they keep counting for whole run, which is good
	found = volname.find("PlasticDet");
	//G4cout << "Found " << found << G4endl;
	if(found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {
	if(postPoint->GetProcessDefinedStep()->GetProcessName() == "hadElastic") {
	fEventAction->totalCounter();
	fEventAction->elasticCounter();
	G4cout << "GetTotalCounter(): " << fEventAction->GetTotalCounter() << G4endl;
	G4cout << "GetElasticCounter(): " << fEventAction->GetElasticCounter() << G4endl;
	}
	if(postPoint->GetProcessDefinedStep()->GetProcessName() == "neutronInelastic" || postPoint->GetProcessDefinedStep()->GetProcessName() == "nCapture" || postPoint->GetProcessDefinedStep()->GetProcessName() == "nFission") {
	fEventAction->totalCounter();
	fEventAction->inelasticCounter();
	G4cout << "GetTotalCounter(): " << fEventAction->GetTotalCounter() << G4endl;
	G4cout << "GetInelasticCounter(): " << fEventAction->GetInelasticCounter() << G4endl;
	}
}

G4int total = fEventAction->GetTotalCounter();
G4int elastic = fEventAction->GetElasticCounter();
G4int inelastic = fEventAction->GetInelasticCounter();
=======
	//Counting hits for efficiencies
	//By not initilalizing counters in event action to zero, they keep counting for whole run, which is good
	found = volname.find("PlasticDet");
	G4cout << "Found " << found << G4endl;
	if(found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {
	if(postPoint->GetProcessDefinedStep()->GetProcessName() == "hadElastic") {
	fEventAction->totalCounter();
	fEventAction->elasticCounter();
	G4cout << "GetTotalCounter(): " << fEventAction->GetTotalCounter() << G4endl;
	G4cout << "GetElasticCounter(): " << fEventAction->GetElasticCounter() << G4endl;
	}
	if(postPoint->GetProcessDefinedStep()->GetProcessName() == "neutronInelastic" || postPoint->GetProcessDefinedStep()->GetProcessName() == "nCapture" || postPoint->GetProcessDefinedStep()->GetProcessName() == "nFission") {
	fEventAction->totalCounter();
	fEventAction->inelasticCounter();
	G4cout << "GetTotalCounter(): " << fEventAction->GetTotalCounter() << G4endl;
	G4cout << "GetInelasticCounter(): " << fEventAction->GetInelasticCounter() << G4endl;
	}
}

G4int total = fEventAction->GetTotalCounter();
G4int elastic = fEventAction->GetElasticCounter();
G4int inelastic = fEventAction->GetInelasticCounter();


	//Counting number of scintillating photons -> Setting to zero at beginning of event ---- Not sure if it does though...
	G4double numScintPhotons;
	G4double numCollectedPhotons;
	found = volname.find("PlasticDet");
	G4cout << "Found " << found << G4endl;
	const std::vector<const G4Track*> *secondaries = aStep->GetSecondaryInCurrentStep();
	if (secondaries->size()>0) {
		for(unsigned int i=0; i<secondaries->size(); ++i) {
			if(secondaries->at(i)->GetParentID()>0) {
				if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())  {
					if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Scintillation") {
						fEventAction->CountOneScintPhoton();
						G4cout<< "GetTotScintPhoton() "<< fEventAction->GetTotScintPhoton() <<G4endl;
					}
				}
			}
		}
	}
numScintPhotons = fEventAction->GetTotScintPhoton();

		//G4cout << "edep out of loop " << edep << G4endl;
		//G4cout << "Has Properties out of loop " << fDetector->HasProperties(volume) << G4endl;
		//G4cout << "Get Properties out of loop " << G4endl;
		//fDetector->GetProperties(volume);
<<<<<<< HEAD
	
	
		//Instead of looking at energy deposited, look at kinetic energy of neutrons before they hit descant
		//Comment OUT! ONLY for DELRIN TEST
		/*
		G4String volume_name = volume->GetName();
		found  = volume_name.find("blueScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
	//	if (aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		//theTrack->SetTrackStatus(fKillTrackAndSecondaries);
		}
		found  = volume_name.find("greenScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		}
		found  = volume_name.find("redScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		}
		found  = volume_name.find("whiteScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		}
		found  = volume_name.find("yellowScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		}
		if(aStep->GetTrack()->GetParentID() > 0 ){ //Should maybe be >0 but oh well
		theTrack->SetTrackStatus(fStopAndKill);
		G4cout<< " Calling kill" << G4endl;
		edep =0;
		}
		*/
>>>>>>> 6828e1eec50f49505f6682d4e9b4096511be47f0

	// check if this volume has its properties set, i.e. it's an active detector
	if((edep > 0 || (fDetector->GridCell() && ekin > 0)) && fDetector->HasProperties(volume)) {
		G4cout << "edep in loop " << edep << G4endl;
=======

	//Counting number of scintillating photons -> Setting to zero at beginning of event ---- Not sure if it does though...
	G4double numScintPhotons;
	G4double numCollectedPhotons;
	found = volname.find("PlasticDet");
	//G4cout << "Found " << found << G4endl;
	const std::vector<const G4Track*> *secondaries = aStep->GetSecondaryInCurrentStep();
	if (secondaries->size()>0) {
		for(unsigned int i=0; i<secondaries->size(); ++i) {
			if(secondaries->at(i)->GetParentID()>0) {
				if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())  {
					if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Scintillation") {
						fEventAction->CountOneScintPhoton();
						G4cout<< "GetTotScintPhoton() "<< fEventAction->GetTotScintPhoton() <<G4endl;
					}
				}
			}
		}
	}
numScintPhotons = fEventAction->GetTotScintPhoton();

		//G4cout << "edep out of loop " << edep << G4endl;
		//G4cout << "Has Properties out of loop " << fDetector->HasProperties(volume) << G4endl;
		//G4cout << "Get Properties out of loop " << G4endl;
		//fDetector->GetProperties(volume);
	
	
		//Instead of looking at energy deposited, look at kinetic energy of neutrons before they hit descant
		//Comment OUT! ONLY for DELRIN TEST
		/*
		G4String volume_name = volume->GetName();
		found  = volume_name.find("blueScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
	//	if (aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		//theTrack->SetTrackStatus(fKillTrackAndSecondaries);
		}
		found  = volume_name.find("greenScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		}
		found  = volume_name.find("redScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		}
		found  = volume_name.find("whiteScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		}
		found  = volume_name.find("yellowScintillatorVolumeLog");
		if (found != G4String::npos && aStep->GetTrack()->GetParentID() == 0 && aStep->IsFirstStepInVolume() == true) {	
		G4cout << "edep before assigning ekin " << edep << G4endl;
		edep = ekin;
		G4cout << "edep after assigning ekin " << edep << G4endl;		
		theTrack->SetTrackStatus(fStopAndKill);
		}
		if(aStep->GetTrack()->GetParentID() > 0 ){ //Should maybe be >0 but oh well
		theTrack->SetTrackStatus(fStopAndKill);
		G4cout<< " Calling kill" << G4endl;
		edep =0;
		}
		*/

	// check if this volume has its properties set, i.e. it's an active detector
	if((edep > 0 || (fDetector->GridCell() && ekin > 0)) && fDetector->HasProperties(volume)) {
		//G4cout << "edep in loop " << edep << G4endl;
<<<<<<< HEAD
	G4cout<<"Calling Get Loop Final Angle: " << fEventAction->GetFinalAngle() << G4endl;
=======
>>>>>>> 2c0eeb1b5df7aa4f21d8f8eff1ea115c68f24a66
>>>>>>> 6828e1eec50f49505f6682d4e9b4096511be47f0
	
		DetectorProperties prop = fDetector->GetProperties(volume);

		if(fDetector->GridCell()) {
			G4String volumeName = volume->GetName();
			if(volumeName.find("gridcellLog") != G4String::npos) {
				// use ekin as edep
				edep = ekin;
				// Now kill the track!
				theTrack->SetTrackStatus(fStopAndKill);
			}
		}
		if(fDetector->Spice()) {
			G4double stepl = 0.;
			if(theTrack->GetDefinition()->GetPDGCharge() != 0.) {
				stepl = aStep->GetStepLength();
			}
			fEventAction->SpiceDet(edep, stepl, prop.detectorNumber, prop.crystalNumber);
		}




		// check edep again in case we use the grid cell but haven't hit it
		//G4cout << "edep " << edep << G4endl; //Testing PLastic fillling ntuple
		if(edep <= 0) return;

<<<<<<< HEAD
		fEventAction->AddHitTracker(prop, evntNb, trackID, parentID, stepNumber, particleType, processType, edep, postPos, postTime, targetZ, total, elastic, inelastic, numScintPhotons, lab_angle, final_angle);

		if(trackSteps) {
			fEventAction->AddStepTracker(prop, evntNb, trackID, parentID, stepNumber, particleType, processType, edep, postPos, postTime, targetZ, total, elastic, inelastic, numScintPhotons, lab_angle, final_angle);
=======
		fEventAction->AddHitTracker(prop, evntNb, trackID, parentID, stepNumber, particleType, processType, edep, postPos, postTime, targetZ, total, elastic, inelastic, numScintPhotons, lab_angle);

		if(trackSteps) {
			fEventAction->AddStepTracker(prop, evntNb, trackID, parentID, stepNumber, particleType, processType, edep, postPos, postTime, targetZ, total, elastic, inelastic, numScintPhotons, lab_angle);
>>>>>>> 6828e1eec50f49505f6682d4e9b4096511be47f0
		}
	}// if(fDetector->HasProperties(volume))
}

