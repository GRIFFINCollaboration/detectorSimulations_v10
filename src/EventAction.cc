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
/// \file analysis/shared/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
// $Id: EventAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "Randomize.hh"
#include "RunAction.hh"

#include <iostream>
#include <fstream>

#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* run)
    :G4UserEventAction(),
      fRunAct(run),
      fPrintModulo(1000)
{
    fNumberOfHits = 0;
    fNumberOfSteps = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt) {  
	fEvtNb = evt->GetEventID();
	if(fEvtNb%fPrintModulo == 0) {
		G4cout<<"---> Begin of event: "<<fEvtNb<<G4endl;
	}

	ClearVariables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*) {
	for (G4int i = 0 ; i < fNumberOfHits; i++) {
		fRunAct->GetHistoManager()->FillHitNtuple(fHitTrackerI[0][i], fHitTrackerI[1][i], fHitTrackerI[2][i], fHitTrackerI[3][i],  fHitTrackerI[4][i], fHitTrackerI[5][i], fHitTrackerI[6][i], fHitTrackerI[7][i], fHitTrackerI[8][i], fHitTrackerD[0][i]/keV, fHitTrackerD[1][i]/mm, fHitTrackerD[2][i]/mm, fHitTrackerD[3][i]/mm, fHitTrackerD[4][i]/second, fHitTrackerI[9][i]);
	}
	for (G4int i = 0 ; i < fNumberOfSteps; i++) {
		fRunAct->GetHistoManager()->FillStepNtuple(fStepTrackerI[0][i], fStepTrackerI[1][i], fStepTrackerI[2][i], fStepTrackerI[3][i],  fStepTrackerI[4][i], fStepTrackerI[5][i], fStepTrackerI[6][i], fStepTrackerI[7][i], fStepTrackerI[8][i], fStepTrackerD[0][i]/keV, fStepTrackerD[1][i]/mm, fStepTrackerD[2][i]/mm, fStepTrackerD[3][i]/mm, fStepTrackerD[4][i]/second, fStepTrackerI[9][i]);
	}

	ClearVariables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::ClearVariables() {
	if(fRunAct->GetHistoManager()->GetStepTrackerBool()) {
		fStepIndex = 0;
		fNumberOfSteps = 0;
		for (G4int i = 0 ; i < MAXSTEPS; i++) {
			for (G4int j = 0 ; j < NUMSTEPVARS; j++) {
				fStepTrackerI[j][i] = 0;
				fStepTrackerD[j][i] = 0.0;
			}
		}
	}

	if(fRunAct->GetHistoManager()->GetHitTrackerBool()) {
		fHitIndex = 0;
		fNumberOfHits = 0;
		fPTrackID = -1;
		fPParentID = -1;

		for (G4int i = 0 ; i < MAXHITS; i++) {
			fPHitMnemonic[i] = "XXX00XX00X";
			for (G4int j = 0 ; j < NUMSTEPVARS; j++) {
				fHitTrackerI[j][i] = 0;
				fHitTrackerD[j][i] = 0.0;
			}
		}
	}
}

void EventAction::AddHitTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ) {
	G4bool newhit = true;
	for (G4int i = 0 ; i < fNumberOfHits; i++) {
		if(fPHitMnemonic[i] == mnemonic) {
			// sum the new enery
			fHitTrackerD[0][i] = fHitTrackerD[0][i] + depEnergy;
			newhit = false;
			break;
		}
	}
	if (newhit) { // new hit
		fPHitMnemonic[fHitIndex] = mnemonic;
		fPTrackID = trackID;
		fPParentID = parentID;
		fHitTrackerI[0][fHitIndex] = eventNumber;
		fHitTrackerI[1][fHitIndex] = trackID;
		fHitTrackerI[2][fHitIndex] = parentID;
		fHitTrackerI[3][fHitIndex] = stepNumber;
		fHitTrackerI[4][fHitIndex] = particleType;
		fHitTrackerI[5][fHitIndex] = processType;
		fHitTrackerI[6][fHitIndex] = systemID;
		fHitTrackerI[7][fHitIndex] = cryNumber;
		fHitTrackerI[8][fHitIndex] = detNumber;
		fHitTrackerI[9][fHitIndex] = targetZ;
		fHitTrackerD[0][fHitIndex] = depEnergy;
		fHitTrackerD[1][fHitIndex] = posx;
		fHitTrackerD[2][fHitIndex] = posy;
		fHitTrackerD[3][fHitIndex] = posz;
		fHitTrackerD[4][fHitIndex] = time;

		fHitIndex++;
		fNumberOfHits = fHitIndex;

		if(fNumberOfHits >= MAXHITS) {
			G4cout << "ERROR! Too many hits!" << G4endl;
		}
	}
}

void EventAction::AddStepTracker(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ) {
	G4bool newstep = true;
	if(newstep) { // new step
		fStepTrackerI[0][fStepIndex] = eventNumber;
		fStepTrackerI[1][fStepIndex] = trackID;
		fStepTrackerI[2][fStepIndex] = parentID;
		fStepTrackerI[3][fStepIndex] = stepNumber;
		fStepTrackerI[4][fStepIndex] = particleType;
		fStepTrackerI[5][fStepIndex] = processType;
		fStepTrackerI[6][fStepIndex] = systemID;
		fStepTrackerI[7][fStepIndex] = cryNumber;
		fStepTrackerI[8][fStepIndex] = detNumber;
		fStepTrackerI[9][fStepIndex] = targetZ;
		fStepTrackerD[0][fStepIndex] = depEnergy;
		fStepTrackerD[1][fStepIndex] = posx;
		fStepTrackerD[2][fStepIndex] = posy;
		fStepTrackerD[3][fStepIndex] = posz;
		fStepTrackerD[4][fStepIndex] = time;

		fStepIndex++;

		fNumberOfSteps = fStepIndex;
		if(fNumberOfSteps >= MAXSTEPS) {
			G4cout << "ERROR! Too many steps!" << G4endl;
		}
	}
}
