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
	fFinalAngle = -1;
	pLabAngle = -1;
	fTOF = 0.;
	SetTOFPos(G4ThreeVector(0.,0.,0.));
	fTOFMulti = 0.;
	SetTOFPosMulti(G4ThreeVector(0.,0.,0.));
	SetFinalMomentum(G4ThreeVector(0.,0.,0.));
	SetInitialMomentum(G4ThreeVector(0.,0.,0.));
	fPEdep =-1.;
	fPEkin=-1.;
	SetTotScintPhoton(0);
	fNumberOfHits = 0;
	fNumberOfSteps = 0;
	SetOldTrackID(-1);
	foldTrackID = -1;

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
	
	//Set Lab angle
	pLabAngle = -1;
	fFinalAngle = -1;
	//Set momentums
	SetFinalMomentum(G4ThreeVector(0.,0.,0.));
	SetInitialMomentum(G4ThreeVector(0.,0.,0.));
	
	//Reset Optical Photon Counter
	SetTotScintPhoton(0);
	//Reset Scatter Counters
	SetElasticCounter(0);
	SetInelasticCounter(0);
	SetTotalCounter(0);
	//Reset TOF
	SetTOF(0);
	SetTOFPos(G4ThreeVector(0.,0.,0.));
	SetTOFMulti(0);
	SetTOFPosMulti(G4ThreeVector(0.,0.,0.));
	//Reset PEdep
	SetPEdep(-1);
	SetPEkin(-1);
	SetOldTrackID(-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*) {
		//G4cout << "End of event activated "  << G4endl;

	if(fHistoManager != nullptr) {
	//G4cout << "testing if statement nullptr " << fNumberOfHits << " " << fNumberOfSteps << G4endl;
		/*if(fHistoManager->GetDetectorConstruction()->Spice()) {
			FillSpice();
		} else*/ {  
			for(G4int i = 0; i < fNumberOfHits; i++) {
		//G4cout << "Filling HitNTuple, testing printing fHitTrackerD[4][i] time " << fHitTrackerD[4][i] << ", plastic hits "<<fPlasticHits<<G4endl;
				G4int j;
				for(j = 0; j < fPlasticHits; ++j) {
//G4cout<<i<<", "<<j<<": "<<fProperties[i].systemID<<", "<<fProperties[i].detectorNumber<<" == "<<fPlasticNumber[j]<<G4endl;
					if(fProperties[i].systemID >= 8700 && fProperties[i].systemID <= 8800 && fProperties[i].detectorNumber == fPlasticNumber[j]) {	
						fHistoManager->FillHitNtuple(fHitTrackerI[0][i], fHitTrackerI[1][i], fHitTrackerI[2][i], fHitTrackerI[3][i],  fHitTrackerI[4][i], fHitTrackerI[5][i], fHitTrackerI[6][i], fHitTrackerI[7][i], fHitTrackerI[8][i], fHitTrackerD[0][i]/keV, fHitTrackerD[1][i]/mm, fHitTrackerD[2][i]/mm, fHitTrackerD[3][i]/mm, fHitTrackerD[4][i]/second, fHitTrackerI[9][i], GetTotalCounter(), GetElasticCounter(), GetInelasticCounter(), GetTotScintPhoton(), fHitTrackerD[9][i]/degree, fHitTrackerD[10][i]/degree, fHitTrackerD[11][i]/nanosecond, fHitTrackerD[12][i]/cm, fHitTrackerD[13][i]/cm, fHitTrackerD[14][i]/cm, fHitTrackerD[15][i]/nanosecond, fHitTrackerD[16][i]/cm, fHitTrackerD[17][i]/cm, fHitTrackerD[18][i]/cm, GetPEkin()/keV, GetPEdep()/keV, fTop1Counter[j], fTop2Counter[j], fTop3Counter[j], fBottom1Counter[j], fBottom2Counter[j], fBottom3Counter[j], fFrontTop1Counter[j], fFrontTop2Counter[j], fFrontMid1Counter[j], fFrontMid2Counter[j], fFrontBottom1Counter[j], fFrontBottom2Counter[j], fCollectionTimeTop1Vector[j], fCollectionTimeTop2Vector[j], fCollectionTimeTop3Vector[j], fCollectionTimeBottom1Vector[j], fCollectionTimeBottom2Vector[j], fCollectionTimeBottom3Vector[j], fCollectionTimeFrontTop1Vector[j], fCollectionTimeFrontTop2Vector[j], fCollectionTimeFrontMid1Vector[j], fCollectionTimeFrontMid2Vector[j], fCollectionTimeFrontBottom1Vector[j], fCollectionTimeFrontBottom2Vector[j], fOpTimeVector[j], fOpEnergyVector[j]);
						break;
					}					
				}			
				if(j == fPlasticHits) {
					fHistoManager->FillHitNtuple(fHitTrackerI[0][i], fHitTrackerI[1][i], fHitTrackerI[2][i], fHitTrackerI[3][i],  fHitTrackerI[4][i], fHitTrackerI[5][i], fHitTrackerI[6][i], fHitTrackerI[7][i], fHitTrackerI[8][i], fHitTrackerD[0][i]/keV, fHitTrackerD[1][i]/mm, fHitTrackerD[2][i]/mm, fHitTrackerD[3][i]/mm, fHitTrackerD[4][i]/second, fHitTrackerI[9][i], GetTotalCounter(), GetElasticCounter(), GetInelasticCounter(), GetTotScintPhoton(), fHitTrackerD[9][i]/degree, fHitTrackerD[10][i]/degree, fHitTrackerD[11][i]/nanosecond, fHitTrackerD[12][i]/cm, fHitTrackerD[13][i]/cm, fHitTrackerD[14][i]/cm, fHitTrackerD[15][i]/nanosecond, fHitTrackerD[16][i]/cm, fHitTrackerD[17][i]/cm, fHitTrackerD[18][i]/cm, GetPEkin()/keV, GetPEdep()/keV);
				}
			}		
			for(G4int i = 0; i < fNumberOfSteps; i++) {
	//	G4cout << "Filling StepNTuple, testing printing fStepTrackerD[4][i] time " << fStepTrackerD[4][i] << G4endl;
				G4int j;
				for(j = 0; j < fPlasticHits; ++j) {
					if(fProperties[i].systemID >= 8700 && fProperties[i].systemID < 8800 && fProperties[i].detectorNumber == fPlasticNumber[j]) {	
				fHistoManager->FillStepNtuple(fStepTrackerI[0][i], fStepTrackerI[1][i], fStepTrackerI[2][i], fStepTrackerI[3][i],  fStepTrackerI[4][i], fStepTrackerI[5][i], fStepTrackerI[6][i], fStepTrackerI[7][i], fStepTrackerI[8][i], fStepTrackerD[0][i]/keV, fStepTrackerD[1][i]/mm, fStepTrackerD[2][i]/mm, fStepTrackerD[3][i]/mm, fStepTrackerD[4][i]/second, fStepTrackerI[9][i], fStepTrackerD[5][i], fStepTrackerD[6][i], fStepTrackerD[7][i], fStepTrackerD[8][i], fStepTrackerD[9][i], fStepTrackerD[10][i], fStepTrackerD[11][i]/nanosecond, fStepTrackerD[12][i]/cm, fStepTrackerD[13][i]/cm, fStepTrackerD[14][i]/cm, fStepTrackerD[15][i]/nanosecond, fStepTrackerD[16][i]/cm, fStepTrackerD[17][i]/cm, fStepTrackerD[18][i]/cm, fStepTrackerD[19][i]/keV, fStepTrackerD[20][i]/keV, fTop1Counter[j], fTop2Counter[j], fTop3Counter[j], fBottom1Counter[j], fBottom2Counter[j], fBottom3Counter[j], fFrontTop1Counter[j], fFrontTop2Counter[j], fFrontMid1Counter[j], fFrontMid2Counter[j], fFrontBottom1Counter[j], fFrontBottom2Counter[j], fCollectionTimeTop1Vector[j], fCollectionTimeTop2Vector[j], fCollectionTimeTop3Vector[j], fCollectionTimeBottom1Vector[j], fCollectionTimeBottom2Vector[j], fCollectionTimeBottom3Vector[j], fCollectionTimeFrontTop1Vector[j], fCollectionTimeFrontTop2Vector[j], fCollectionTimeFrontMid1Vector[j], fCollectionTimeFrontMid2Vector[j], fCollectionTimeFrontBottom1Vector[j], fCollectionTimeFrontBottom2Vector[j], fOpTimeVector[j], fOpEnergyVector[j]);
						break;
					}
			}

			if(j==fPlasticHits){
				fHistoManager->FillStepNtuple(fStepTrackerI[0][i], fStepTrackerI[1][i], fStepTrackerI[2][i], fStepTrackerI[3][i],  fStepTrackerI[4][i], fStepTrackerI[5][i], fStepTrackerI[6][i], fStepTrackerI[7][i], fStepTrackerI[8][i], fStepTrackerD[0][i]/keV, fStepTrackerD[1][i]/mm, fStepTrackerD[2][i]/mm, fStepTrackerD[3][i]/mm, fStepTrackerD[4][i]/second, fStepTrackerI[9][i], fStepTrackerD[5][i], fStepTrackerD[6][i], fStepTrackerD[7][i], fStepTrackerD[8][i], fStepTrackerD[9][i], fStepTrackerD[10][i], fStepTrackerD[11][i]/nanosecond, fStepTrackerD[12][i]/cm, fStepTrackerD[13][i]/cm, fStepTrackerD[14][i]/cm, fStepTrackerD[15][i]/nanosecond, fStepTrackerD[16][i]/cm, fStepTrackerD[17][i]/cm, fStepTrackerD[18][i]/cm, fStepTrackerD[19][i]/keV, fStepTrackerD[20][i]/keV);
				}
		}
		}
		ClearVariables();
	}


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddHitTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& targetZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos, G4double TOFMulti, G4ThreeVector TOFPosMulti, G4double PEkin, G4double PEdep) {

//      G4cout << "ParentID,  edep, time " << parentID << "  " << depEnergy << "  " << time << G4endl;

	for(G4int i = 0; i < fNumberOfHits; i++) {
		if(fProperties[i] == properties) {
			// sum the new enery
			// Optical photons do not obey energy conservation and their energy must not be tallied
			//if (particleType == 8) continue;
			fHitTrackerD[0][i] = fHitTrackerD[0][i] + depEnergy;
			//fHitTrackerD[20][i] = fHitTrackerD[20][i] + PEdep;

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
	fHitTrackerD[5][fNumberOfHits] = total;
	fHitTrackerD[6][fNumberOfHits] = elastic;
	fHitTrackerD[7][fNumberOfHits] = inelastic;
	fHitTrackerD[8][fNumberOfHits] = numScintPhotons;
	fHitTrackerD[9][fNumberOfHits] = lab_angle; //pLabAngle?
	fHitTrackerD[10][fNumberOfHits] = final_angle; //fFinalAngle?
	fHitTrackerD[11][fNumberOfHits] = TOF; //based on first scatter
	fHitTrackerD[12][fNumberOfHits] = TOFPos.x(); //based on first scatter
	fHitTrackerD[13][fNumberOfHits] = TOFPos.y(); //based on first scatter
	fHitTrackerD[14][fNumberOfHits] = TOFPos.z(); //based on first scatter
	fHitTrackerD[15][fNumberOfHits] = TOFMulti;  //based on >2 scatter
	fHitTrackerD[16][fNumberOfHits] = TOFPosMulti.x(); //based on >2 scatter
	fHitTrackerD[17][fNumberOfHits] = TOFPosMulti.y(); //based on >2 scatter
	fHitTrackerD[18][fNumberOfHits] = TOFPosMulti.z(); //based on >2 scatter
	//fHitTrackerD[19][fNumberOfHits] = PEkin; //kinetic energy of neutron in plastic
	fHitTrackerD[19][fNumberOfHits] = GetPEkin(); //kinetic energy of neutron in plastic
	fHitTrackerD[20][fNumberOfHits] = GetPEdep(); //energy is plastic

	++fNumberOfHits;

	if(fNumberOfHits >= MAXHITS) {
		G4cout<<"ERROR! Too many hits!"<<G4endl;
		throw;
	}
}

void EventAction::AddStepTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& targetZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos, G4double TOFMulti, G4ThreeVector TOFPosMulti, G4double PEkin, G4double PEdep) {
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
	fStepTrackerD[5][fNumberOfSteps] = total;
	fStepTrackerD[6][fNumberOfSteps] = elastic;
	fStepTrackerD[7][fNumberOfSteps] = inelastic;
	fStepTrackerD[8][fNumberOfSteps] = numScintPhotons;
	fStepTrackerD[9][fNumberOfSteps] = lab_angle;
	fStepTrackerD[10][fNumberOfSteps] = final_angle;
	fStepTrackerD[11][fNumberOfSteps] = TOF;
	fStepTrackerD[12][fNumberOfSteps] = TOFPos.x();
	fStepTrackerD[13][fNumberOfSteps] = TOFPos.y();
	fStepTrackerD[14][fNumberOfSteps] = TOFPos.z();
	fStepTrackerD[15][fNumberOfSteps] = TOFMulti;
	fStepTrackerD[16][fNumberOfSteps] = TOFPosMulti.x();
	fStepTrackerD[17][fNumberOfSteps] = TOFPosMulti.y();
	fStepTrackerD[18][fNumberOfSteps] = TOFPosMulti.z();
	fStepTrackerD[19][fNumberOfSteps] = PEkin;
	fStepTrackerD[20][fNumberOfSteps] = PEdep;

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

	if(fHistoManager->GetDetectorConstruction()->Spice()) {
		for(G4int i = 0; i < 10; i++){
			for(G4int j=0; j<12; j++) {
				fSpiceEnergyDet[i][j] = 0;
				fSpiceTrackDet[i][j]  = 0;
			}
		}
	}

	for(int i=0; i<fPlasticHits; ++i){
	fCollectionTimeTop1Vector[i].clear();
	fCollectionTimeTop2Vector[i].clear();
	fCollectionTimeTop3Vector[i].clear();
	fCollectionTimeBottom1Vector[i].clear();
	fCollectionTimeBottom2Vector[i].clear();
	fCollectionTimeBottom3Vector[i].clear();
	fCollectionTimeFrontTop1Vector[i].clear();
	fCollectionTimeFrontTop2Vector[i].clear();
	fCollectionTimeFrontMid1Vector[i].clear();
	fCollectionTimeFrontMid2Vector[i].clear();
	fCollectionTimeFrontBottom1Vector[i].clear();
	fCollectionTimeFrontBottom2Vector[i].clear();
	//fPlasticNumber[i].Clear();
	fPlasticNumber[i]=-1;
	fTop1Counter[i] = 0;
	fTop2Counter[i] = 0;
	fTop3Counter[i] = 0;
	fBottom1Counter[i] = 0;
	fBottom2Counter[i] = 0;
	fBottom3Counter[i] = 0;
	fFrontTop1Counter[i] = 0;
	fFrontTop2Counter[i] = 0;
	fFrontMid1Counter[i] = 0;
	fFrontMid2Counter[i] = 0;
	fFrontBottom1Counter[i] = 0;
	fFrontBottom2Counter[i] = 0;
	//fOpTimeVector[i].clear();
	//fOpEnergyVector[i].clear();
	}
	fPlasticHits = 0;
	for(int i=0; i<fPlasticHits2; ++i){
	fOpTimeVector[i].clear();
	fOpEnergyVector[i].clear();
	}
	fPlasticHits2 = 0;



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

bool EventAction::PhotonDetectionEfficiency(G4double ekin) {
	G4double uni = G4UniformRand();
	G4double wave = 1239.8/(ekin*1.e6); //nm

	//G4cout << "Uni: " << uni << G4endl;
	//G4cout << "ekin: " << ekin << G4endl;
	//G4cout << "wavelength: " << wave << G4endl;
        G4double prob = 0;
	G4bool pde = false;
	if (wave < 250) prob = 0.05;
	else if (wave >= 250 && wave < 350) prob = 0.2;
	else if (wave >= 325 && wave < 350) prob = 0.3;
	else if (wave >= 350 && wave < 375) prob = 0.35;
	else if (wave >= 375 && wave < 390) prob = 0.38;
	else if (wave >= 390 && wave < 460) prob = 0.4;
	else if (wave >= 460 && wave < 475) prob = 0.38;
	else if (wave >= 475 && wave < 500) prob = 0.32;
	else if (wave >= 500 && wave < 525) prob = 0.28;
	else if (wave >= 525 && wave < 550) prob = 0.2;
	else if (wave >= 550 && wave < 625) prob = 0.15;
	else if (wave >= 625 && wave < 675) prob = 0.1;
	else if (wave > 675) prob = 0.05;
	
	//prob = 0.4;
	
	if (uni <= prob) pde = true;
	//pde=true;//Using the scaling instead to limit computation time
	return pde;
}
void EventAction::SetScintPhotonTimeTop1(G4double top, G4int detNum) {
	//	G4cout << "Top hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i]== detNum) {

			fCollectionTimeTop1Vector[i].push_back(top);
			fTop1Counter[i]++;
//	G4cout << "fTop1Counter function: " << fTop1Counter[i] << G4endl;
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeTop1Vector[fPlasticHits].push_back(top);
	fTop1Counter[fPlasticHits]++;
	fPlasticHits++;
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}
void EventAction::SetScintPhotonTimeTop2(G4double top, G4int detNum) {
//	G4cout << "Top hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i]== detNum) {

			fCollectionTimeTop2Vector[i].push_back(top);
			fTop2Counter[i]++;
//	G4cout << "fTop2Counter function: " << fTop2Counter[i] << G4endl;
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeTop2Vector[fPlasticHits].push_back(top);
	fTop2Counter[fPlasticHits]++;
	fPlasticHits++;
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}
void EventAction::SetScintPhotonTimeTop3(G4double top, G4int detNum) {
//	G4cout << "Top hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i]== detNum) {

			fCollectionTimeTop3Vector[i].push_back(top);
			fTop3Counter[i]++;
//	G4cout << "fTop3Counter function: " << fTop3Counter[i] << G4endl;
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeTop3Vector[fPlasticHits].push_back(top);
	fTop3Counter[fPlasticHits]++;
	fPlasticHits++;
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}
void EventAction::SetScintPhotonTimeBottom1(G4double bottom, G4int detNum) {
//	G4cout << "Bottom hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i] == detNum) {

			fCollectionTimeBottom1Vector[i].push_back(bottom);
			fBottom1Counter[i]++;
//	G4cout << "fBottom1Counter function: " << fBottom1Counter[i] << G4endl;
			
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeBottom1Vector[fPlasticHits].push_back(bottom);
	fBottom1Counter[fPlasticHits]++;
	fPlasticHits++;
	
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}
void EventAction::SetScintPhotonTimeBottom2(G4double bottom, G4int detNum) {
//	G4cout << "Bottom hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i] == detNum) {

			fCollectionTimeBottom2Vector[i].push_back(bottom);
			fBottom2Counter[i]++;
//	G4cout << "fBottom2Counter function: " << fBottom2Counter[i] << G4endl;
			
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeBottom2Vector[fPlasticHits].push_back(bottom);
	fBottom2Counter[fPlasticHits]++;
	fPlasticHits++;
	
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}
void EventAction::SetScintPhotonTimeBottom3(G4double bottom, G4int detNum) {
//	G4cout << "Bottom hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i] == detNum) {

			fCollectionTimeBottom3Vector[i].push_back(bottom);
			fBottom3Counter[i]++;
//	G4cout << "fBottom3Counter function: " << fBottom3Counter[i] << G4endl;
			
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeBottom3Vector[fPlasticHits].push_back(bottom);
	fBottom3Counter[fPlasticHits]++;
	fPlasticHits++;
	
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}


void EventAction::SetScintPhotonTimeFrontTop1(G4double frontTop, G4int detNum) {
//	G4cout << "Front Top hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i]== detNum) {

			fCollectionTimeFrontTop1Vector[i].push_back(frontTop);
			fFrontTop1Counter[i]++;
//	G4cout << "fFrontTop1Counter function: " << fFrontTop1Counter[i] << G4endl;
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeFrontTop1Vector[fPlasticHits].push_back(frontTop);
	fFrontTop1Counter[fPlasticHits]++;
	fPlasticHits++;
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}
void EventAction::SetScintPhotonTimeFrontTop2(G4double frontTop, G4int detNum) {
//	G4cout << "Front Top hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i]== detNum) {

			fCollectionTimeFrontTop2Vector[i].push_back(frontTop);
			fFrontTop2Counter[i]++;
//	G4cout << "fFrontTop2Counter function: " << fFrontTop2Counter[i] << G4endl;
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeFrontTop2Vector[fPlasticHits].push_back(frontTop);
	fFrontTop2Counter[fPlasticHits]++;
	fPlasticHits++;
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}

void EventAction::SetScintPhotonTimeFrontMid1(G4double frontMid, G4int detNum) {
//	G4cout << "Front Mid hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i]== detNum) {

			fCollectionTimeFrontMid1Vector[i].push_back(frontMid);
			fFrontMid1Counter[i]++;
//	G4cout << "Front Mid1 Counter function: " << fFrontMid1Counter[i] << G4endl;
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeFrontMid1Vector[fPlasticHits].push_back(frontMid);
	fFrontMid1Counter[fPlasticHits]++;
	fPlasticHits++;
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}
void EventAction::SetScintPhotonTimeFrontMid2(G4double frontMid, G4int detNum) {
//	G4cout << "Front Mid hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i]== detNum) {

			fCollectionTimeFrontMid2Vector[i].push_back(frontMid);
			fFrontMid2Counter[i]++;
//	G4cout << "Front Mid2 Counter function: " << fFrontMid2Counter[i] << G4endl;
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeFrontMid2Vector[fPlasticHits].push_back(frontMid);
	fFrontMid2Counter[fPlasticHits]++;
	fPlasticHits++;
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}

void EventAction::SetScintPhotonTimeFrontBottom1(G4double frontBot, G4int detNum) {
//	G4cout << "Front Bot hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i]== detNum) {

			fCollectionTimeFrontBottom1Vector[i].push_back(frontBot);
			fFrontBottom1Counter[i]++;
//	G4cout << "fFrontBottom1Counter function: " << fFrontBottom1Counter[i] << G4endl;
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeFrontBottom1Vector[fPlasticHits].push_back(frontBot);
	fFrontBottom1Counter[fPlasticHits]++;
	fPlasticHits++;
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}
void EventAction::SetScintPhotonTimeFrontBottom2(G4double frontBot, G4int detNum) {
//	G4cout << "Front Bot hit: " << fPlasticHits << G4endl;
	for(G4int i = 0; i < fPlasticHits; i++) {
		if(fPlasticNumber[i]== detNum) {

			fCollectionTimeFrontBottom2Vector[i].push_back(frontBot);
			fFrontBottom2Counter[i]++;
//	G4cout << "fFrontBottom2Counter function: " << fFrontBottom2Counter[i] << G4endl;
			return;
		}
	}

	fPlasticNumber[fPlasticHits] = detNum;
	fCollectionTimeFrontBottom2Vector[fPlasticHits].push_back(frontBot);
	fFrontBottom2Counter[fPlasticHits]++;
	fPlasticHits++;
	if(fPlasticHits >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics hits!"<<G4endl;
		throw;
	}
}


void EventAction::SetScintPhotonEnergyTime(G4double OpTime, G4double OpEnergy, G4int detNum) {
	for(G4int i = 0; i < fPlasticHits2; i++) {
		if(fPlasticNumber2[i] == detNum) {

			fOpTimeVector[i].push_back(OpTime);
			fOpEnergyVector[i].push_back(OpEnergy);
			
			return;
		}
	}

	fPlasticNumber2[fPlasticHits2] = detNum;
	fOpTimeVector[fPlasticHits2].push_back(OpTime);
	fOpEnergyVector[fPlasticHits2].push_back(OpEnergy);
	fPlasticHits2++;
	
	if(fPlasticHits2 >= MAXHITS) {
		G4cout<<"ERROR! Too many plastics2 hits!"<<G4endl;
		throw;
	}
}




