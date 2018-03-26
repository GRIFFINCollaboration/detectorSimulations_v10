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

	 if(fHistoManager->GetDetectorConstruction()->Spice()) {
		 SetupSpiceErfc();
	 }
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

void EventAction::EndOfEventAction(const G4Event*) {
	if(fHistoManager != nullptr) {
		if(fHistoManager->GetDetectorConstruction()->Spice()) {
			FillSpice();
		} else {
			for(G4int i = 0; i < fNumberOfHits; i++) {
				fHistoManager->FillHitNtuple(fHitTrackerI[0][i], fHitTrackerI[1][i], fHitTrackerI[2][i], fHitTrackerI[3][i],  fHitTrackerI[4][i], fHitTrackerI[5][i], fHitTrackerI[6][i], fHitTrackerI[7][i], fHitTrackerI[8][i], fHitTrackerD[0][i]/keV, fHitTrackerD[1][i]/mm, fHitTrackerD[2][i]/mm, fHitTrackerD[3][i]/mm, fHitTrackerD[4][i]/second, fHitTrackerI[9][i]);
			}
			for(G4int i = 0; i < fNumberOfSteps; i++) {
				fHistoManager->FillStepNtuple(fStepTrackerI[0][i], fStepTrackerI[1][i], fStepTrackerI[2][i], fStepTrackerI[3][i],  fStepTrackerI[4][i], fStepTrackerI[5][i], fStepTrackerI[6][i], fStepTrackerI[7][i], fStepTrackerI[8][i], fStepTrackerD[0][i]/keV, fStepTrackerD[1][i]/mm, fStepTrackerD[2][i]/mm, fStepTrackerD[3][i]/mm, fStepTrackerD[4][i]/second, fStepTrackerI[9][i]);
			}
		}

		ClearVariables();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddHitTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ) {
	G4bool newhit = true;
	for(G4int i = 0; i < fNumberOfHits; i++) {
		if(fPHitMnemonic[i] == mnemonic) {
			// sum the new enery
			fHitTrackerD[0][i] = fHitTrackerD[0][i] + depEnergy;
			newhit = false;
			break;
		}
	}
	if(newhit) { // new hit
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

void EventAction::ClearVariables() {
	if(fHistoManager->GetStepTrackerBool()) {
		fStepIndex = 0;
		fNumberOfSteps = 0;
		for(G4int i = 0; i < MAXSTEPS; i++) {
			for(G4int j = 0; j < NUMSTEPVARS; j++) {
				fStepTrackerI[j][i] = 0;
				fStepTrackerD[j][i] = 0.0;
			}
		}
	}

	if(fHistoManager->GetHitTrackerBool()) {
		fHitIndex = 0;
		fNumberOfHits = 0;
		fPTrackID = -1;
		fPParentID = -1;

		for(G4int i = 0; i < MAXHITS; i++) {
			fPHitMnemonic[i] = "XXX00XX00X";
			for(G4int j = 0; j < NUMSTEPVARS; j++) {
				fHitTrackerI[j][i] = 0;
				fHitTrackerD[j][i] = 0.0;
			}
		}
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


void EventAction::FillSpice() {
	G4double energySumDet = 0;
	G4int fSpiceMultiplicity = 0;
	G4double SpiceEnergy,SpiceEnergyRaw;
	for(G4int ring=0; ring < MAXNUMDETSPICE; ring++) {
		for(G4int seg=0; seg < 12; seg++) {
			if(fSpiceEnergyDet[ring][seg] > 20.*CLHEP::keV&&WRITEEDEPHISTOS) {
				SpiceEnergyRaw = fSpiceEnergyDet[ring][seg];// SpiceRawEnergy. No resolution applied
				SpiceEnergy=SpiceEnergyRaw;
				if(fHistoManager->GetDetectorConstruction()->SpiceRes()){
					SpiceEnergy = ApplySpiceRes(SpiceEnergyRaw);
					fHistoManager->FillHistogram(fHistoManager->AngleDistro(8),SpiceEnergy*1000);
				}

				//fill energies across all segments
				fHistoManager->FillHistogram(fHistoManager->SpiceHistNumbers(0),(SpiceEnergy)*1000.);//+4keV for Bismuth
				//fill standard energy spectra - no add-back

				//2D all seg energies
				fHistoManager->Fill2DHistogram(fHistoManager->AngleDistro(1), (G4double) MAXNUMSEGSPICE*ring+seg,SpiceEnergy*1000., 1.0);


				fSpiceMultiplicity++;//iterates every deposition per fill to track multiplicity per event
				energySumDet += SpiceEnergy;

				//Gated angular histos
				if(abs(SpiceEnergyRaw - fHistoManager->BeamEnergy()) < 10.*keV){//gating kinematics on (raw) detected energy

					//Mapping phi to overlay mirrored segments
					G4int remainder = seg%3;
					G4double PhiMap = ((int)seg/3)*CLHEP::pi/2;
					PhiMap = fHistoManager->BeamPhi()-PhiMap;
					PhiMap = asin(sin(PhiMap));//if we went over 360

					//THETA vs. PHI
					fHistoManager->Fill2DHistogram(fHistoManager->SpiceAngleHists(MAXNUMSEGSPICE*ring+remainder), fHistoManager->BeamTheta(), PhiMap, 1.);//PHI M<AP!!!!!!!!!!!!!
				}


			}//end of filling histos with SPICE detections per individual seg
		}//seg loop
	}//ring loop

	//wiritng add-backed energies
	if(energySumDet > MINENERGYTHRES) {//after exiting loops for all rings/segs, will input summed (add-back) energy if above threshold
		if(WRITEEDEPHISTOS)     fHistoManager->FillHistogram(fHistoManager->SpiceHistNumbers(1), energySumDet*1000.);
	}
}


G4bool EventAction::SpiceTest(){//is SPICE inputted
	return fHistoManager->GetDetectorConstruction()->Spice();
}


G4double EventAction::ApplySpiceRes(G4double energy) {
	//pre-calculated functions that scale with energy to give accurate resolution split
	//depending on result, will be part of one of two Gaussians, or a ERFC-bsaed decay
	G4double gaussPick = G4UniformRand(); 
	G4double result = 0.;
	if(gaussPick < 37./100.) {
		G4double Sigma1 = (0.7E-6*energy/CLHEP::keV + 0.0009);//include ->keV conversion
		result = G4RandGauss::shoot(energy, Sigma1);
	} else if(gaussPick > 77./100.){
		result = SpiceErfc()*(0.5E-6*energy + 0.0012)+energy/CLHEP::keV/1000.;//ERFC from (from http://radware.phy.ornl.gov/gf3/)
	} else {
		G4double Sigma2 = (1.8E-6*energy/CLHEP::keV + 0.00500);
		result = G4RandGauss::shoot(energy, Sigma2);
	}
	return result;
}

void EventAction::SetupSpiceErfc() {
	//10000 bins between -0.5 and 0.5 containing decay function information for SPICE pre-amps
	fAmpTot = 0.;

	//making fixed ERFC func - accessed via a getter (local to this class)
	for(int i = 0; i < 10000; ++i) {
		//Full function from -80 to +20
		double channel = static_cast<double>(i);
		channel /=100.;
		channel -= 80.;
		fAmp[i] = exp(channel/30.)*erfc(channel/sqrt(2.)/sqrt(2.) + 1./(sqrt(2.)*30.));
		fAmpx[i] = channel; //have a corresponding x-array to the Ampitudes
		fAmpTot += fAmp[i];
	}
}

G4double EventAction::SpiceErfc() {
	G4double erfcRand = G4UniformRand()*fAmpTot;
	G4double c = 0.;
	for(int i = 0; i < 10000; ++i) {
		c += fAmp[i];
		if(c > erfcRand) {
			return fAmpx[i];
		}
	}

	return 0.;
}
