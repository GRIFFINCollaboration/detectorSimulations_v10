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

EventAction::EventAction(RunAction* run, HistoManager* hist)
:G4UserEventAction(),
	fRunAction(run),
	fHistoManager(hist),
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

	if(fHistoManager != nullptr) ClearVariables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt) {
	if(fHistoManager != nullptr) {
		for(G4int i = 0; i < fNumberOfHits; i++) {
			fHistoManager->FillHitNtuple(fHitTrackerI[0][i], fHitTrackerI[1][i], fHitTrackerI[2][i], fHitTrackerI[3][i],  fHitTrackerI[4][i], fHitTrackerI[5][i], fHitTrackerI[6][i], fHitTrackerI[7][i], fHitTrackerI[8][i], fHitTrackerD[0][i]/keV, fHitTrackerD[1][i]/mm, fHitTrackerD[2][i]/mm, fHitTrackerD[3][i]/mm, fHitTrackerD[4][i]/second, fHitTrackerI[9][i]);
		}
		for(G4int i = 0; i < fNumberOfSteps; i++) {
			fHistoManager->FillStepNtuple(fStepTrackerI[0][i], fStepTrackerI[1][i], fStepTrackerI[2][i], fStepTrackerI[3][i],  fStepTrackerI[4][i], fStepTrackerI[5][i], fStepTrackerI[6][i], fStepTrackerI[7][i], fStepTrackerI[8][i], fStepTrackerD[0][i]/keV, fStepTrackerD[1][i]/mm, fStepTrackerD[2][i]/mm, fStepTrackerD[3][i]/mm, fStepTrackerD[4][i]/second, fStepTrackerI[9][i]);
		}

		if(fNumberOfHits == 0 && fHistoManager->RecordAll()) fHistoManager->FillHitNtuple(evt->GetEventID());
		ClearVariables();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddHitTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& targetZ, const G4double& kinEnergy) {
	for(G4int i = 0; i < fNumberOfHits; i++) {
		if(fProperties[i] == properties) {
			// sum the new enery
			fHitTrackerD[0][i] = fHitTrackerD[0][i] + depEnergy;
			return;
		}
	}
	// new hit
	fProperties[fNumberOfHits] = properties;
	fHitTrackerI[0][fNumberOfHits] = eventNumber;
	fHitTrackerI[1][fNumberOfHits] = trackID;
	fHitTrackerI[2][fNumberOfHits] = parentID;
	fHitTrackerI[3][fNumberOfHits] = stepNumber;
	fHitTrackerI[4][fNumberOfHits] = particleType;
	fHitTrackerI[5][fNumberOfHits] = processType;
	fHitTrackerI[6][fNumberOfHits] = properties.systemID;
	fHitTrackerI[7][fNumberOfHits] = properties.crystalNumber;
	fHitTrackerI[8][fNumberOfHits] = properties.detectorNumber;
	fHitTrackerI[9][fNumberOfHits] = targetZ;
	fHitTrackerD[0][fNumberOfHits] = depEnergy;
	fHitTrackerD[1][fNumberOfHits] = pos.x();
	fHitTrackerD[2][fNumberOfHits] = pos.y();
	fHitTrackerD[3][fNumberOfHits] = pos.z();
	fHitTrackerD[4][fNumberOfHits] = time;

	++fNumberOfHits;

	if(fNumberOfHits >= MAXHITS) {
		G4cout<<"ERROR! Too many hits!"<<G4endl;
		throw;
	}

	if(fHistoManager->GetDetectorConstruction()->Descant() || fHistoManager->GetDetectorConstruction()->Testcan()) {
		fHistoManager->PushBack(depEnergy, kinEnergy, particleType); //, time, processType);
	}
}

void EventAction::AddStepTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& targetZ) {
	// new step
	fStepTrackerI[0][fNumberOfSteps] = eventNumber;
	fStepTrackerI[1][fNumberOfSteps] = trackID;
	fStepTrackerI[2][fNumberOfSteps] = parentID;
	fStepTrackerI[3][fNumberOfSteps] = stepNumber;
	fStepTrackerI[4][fNumberOfSteps] = particleType;
	fStepTrackerI[5][fNumberOfSteps] = processType;
	fStepTrackerI[6][fNumberOfSteps] = properties.systemID;
	fStepTrackerI[7][fNumberOfSteps] = properties.crystalNumber;
	fStepTrackerI[8][fNumberOfSteps] = properties.detectorNumber;
	fStepTrackerI[9][fNumberOfSteps] = targetZ;
	fStepTrackerD[0][fNumberOfSteps] = depEnergy;
	fStepTrackerD[1][fNumberOfSteps] = pos.x();
	fStepTrackerD[2][fNumberOfSteps] = pos.y();
	fStepTrackerD[3][fNumberOfSteps] = pos.z();
	fStepTrackerD[4][fNumberOfSteps] = time;

	++fNumberOfSteps;

	if(fNumberOfSteps >= MAXSTEPS) {
		G4cout<<"ERROR! Too many steps!"<<G4endl;
		throw;
	}
}

void EventAction::ClearVariables() {
	if(fHistoManager->GetStepTrackerBool()) {
		fNumberOfSteps = 0;
		for(G4int i = 0; i < MAXSTEPS; i++) {
			for(G4int j = 0; j < NUMSTEPVARS; j++) {
				fStepTrackerI[j][i] = 0;
				fStepTrackerD[j][i] = 0.0;
			}
		}
	}

	if(fHistoManager->GetHitTrackerBool()) {
		fNumberOfHits = 0;

		for(G4int i = 0; i < MAXHITS; i++) {
			fProperties[i].Clear();
			for(G4int j = 0; j < NUMSTEPVARS; j++) {
				fHitTrackerI[j][i] = 0;
				fHitTrackerD[j][i] = 0.0;
			}
		}
	}

	if(fHistoManager->GetDetectorConstruction()->Descant() || fHistoManager->GetDetectorConstruction()->Testcan()) {
		fHistoManager->ClearVariables();
	}

	if(fHistoManager->GetDetectorConstruction()->Spice()) {
		for(G4int i = 0; i < 10; i++){
			for(G4int j=0; j<12; j++) {
				fSpiceEnergyDet[i][j] = 0;
				fSpiceTrackDet[i][j]  = 0;
			}
		}
	}
}


//void EventAction::FillSpice() {
//	G4double energySumDet = 0;
//	G4int fSpiceMultiplicity = 0;
//	G4double SpiceEnergy,SpiceEnergyRaw;
//	for(G4int ring=0; ring < MAXNUMDETSPICE; ring++) {
//		for(G4int seg=0; seg < 12; seg++) {
//			if(fSpiceEnergyDet[ring][seg] > 20.*CLHEP::keV&&WRITEEDEPHISTOS) {
//				SpiceEnergyRaw = fSpiceEnergyDet[ring][seg];// SpiceRawEnergy. No resolution applied
//				SpiceEnergy=SpiceEnergyRaw;
//				if(fHistoManager->GetDetectorConstruction()->SpiceRes()){
//					SpiceEnergy = ApplySpiceRes(SpiceEnergyRaw);
//					fHistoManager->FillHistogram(fHistoManager->AngleDistro(8),SpiceEnergy*1000);
//				}
//
//				//fill energies across all segments
//				fHistoManager->FillHistogram(fHistoManager->SpiceHistNumbers(0),(SpiceEnergy)*1000.);//+4keV for Bismuth
//				//fill standard energy spectra - no add-back
//
//				//2D all seg energies
//				fHistoManager->Fill2DHistogram(fHistoManager->AngleDistro(1), (G4double) MAXNUMSEGSPICE*ring+seg,SpiceEnergy*1000., 1.0);
//
//
//				fSpiceMultiplicity++;//iterates every deposition per fill to track multiplicity per event
//				energySumDet += SpiceEnergy;
//
//				//Gated angular histos
//				if(abs(SpiceEnergyRaw - fHistoManager->BeamEnergy()) < 10.*keV){//gating kinematics on (raw) detected energy
//
//					//Mapping phi to overlay mirrored segments
//					G4int remainder = seg%3;
//					G4double PhiMap = ((int)seg/3)*CLHEP::pi/2;
//					PhiMap = fHistoManager->BeamPhi()-PhiMap;
//					PhiMap = asin(sin(PhiMap));//if we went over 360
//
//					//THETA vs. PHI
//					fHistoManager->Fill2DHistogram(fHistoManager->SpiceAngleHists(MAXNUMSEGSPICE*ring+remainder), fHistoManager->BeamTheta(), PhiMap, 1.);//PHI M<AP!!!!!!!!!!!!!
//				}
//
//
//			}//end of filling histos with SPICE detections per individual seg
//		}//seg loop
//	}//ring loop
//
//	//wiritng add-backed energies
//	if(energySumDet > MINENERGYTHRES) {//after exiting loops for all rings/segs, will input summed (add-back) energy if above threshold
//		if(WRITEEDEPHISTOS)     fHistoManager->FillHistogram(fHistoManager->SpiceHistNumbers(1), energySumDet*1000.);
//	}
//}


G4bool EventAction::SpiceTest(){//is SPICE inputted
	return fHistoManager->GetDetectorConstruction()->Spice();
}
