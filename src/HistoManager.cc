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
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id: HistoManager.cc 74272 2013-10-02 14:48:50Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager(DetectorConstruction* detectorConstruction) {
	fFileName = "g4out";
	fFactoryOn = false;

	// Only fill one NTuple at a time. If fStepTrackerBool is true, then fHitTrackerBool should be false, or vise-versa.
	// There is no need to have the hit NTuple and the step NTuple.
	fHitTrackerBool = true;
	fStepTrackerBool = false;
	// ntuple
	for(G4int k=0; k<MAXNTCOL; k++) {
		//        fNtColId[k] = 0;
		fNtColIdHit[k] = 0;
		fNtColIdStep[k] = 0;
	}

	fDetectorConstruction = detectorConstruction;
	fMessenger = new HistoMessenger(this);
	fFirstRecordingId = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
	delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book() {
	// Create or get analysis manager
	// The choice of analysis technology is done via selection of a namespace
	// in HistoManager.hh
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	// this merges the root files at the end, but puts all branches inside a branch called ntuple ...
	//analysisManager->SetNtupleMerging(true);
	//analysisManager->SetVerboseLevel(4);

	// Create directories
	// Open an output file
	G4bool fileOpen = analysisManager->OpenFile(fFileName); 
	if(!fileOpen) {
		G4String extension = analysisManager->GetFileType();
		G4cout<<"---> HistoManager::book(): cannot open "<<fFileName<<"."<<extension<<G4endl;
		return;
	}

	///////////////////////////////////////////////////////////////////
	// Create 1 ntuple
	if(fHitTrackerBool) {
		analysisManager->CreateNtuple("ntuple", "HitTracker");
		int colId = 0;
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("eventNumber");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("trackID");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("parentID");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("stepNumber");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("particleType");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("processType");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("systemID");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("cryNumber");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("detNumber");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("depEnergy");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("posx");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("posy");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("posz");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("time");
		fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("targetZ"); // 14 here
		if(fDetectorConstruction->Descant() || fDetectorConstruction->Testcan()) {
			G4cout<<"Filling descant settings after "<<colId<<G4endl;
			fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("eDepVector", fEdepVector);
			fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("eKinVector", fEkinVector);
			fNtColIdHit[colId++] = analysisManager->CreateNtupleIColumn("particleTypeVector", fParticleTypeVector); // 17 here
		}
		if(fRecordGun) {
			fFirstRecordingId = colId; // this is actually one less than the first recording ID!
			G4cout<<"Creating gun settings after "<<fFirstRecordingId<<G4endl;
			fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("primaryE");
			fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("primaryTheta");
			fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("primaryPhi");
			fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("originX");
			fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("originY");
			fNtColIdHit[colId++] = analysisManager->CreateNtupleDColumn("originZ"); // 20 or 23 here
		}
		analysisManager->FinishNtuple();
	}


	if(fStepTrackerBool) {
		analysisManager->CreateNtuple("ntuple", "StepTracker");
		fNtColIdStep[0] = analysisManager->CreateNtupleIColumn("eventNumber");
		fNtColIdStep[1] = analysisManager->CreateNtupleIColumn("trackID");
		fNtColIdStep[2] = analysisManager->CreateNtupleIColumn("parentID");
		fNtColIdStep[3] = analysisManager->CreateNtupleIColumn("stepNumber");
		fNtColIdStep[4] = analysisManager->CreateNtupleIColumn("particleType");
		fNtColIdStep[5] = analysisManager->CreateNtupleIColumn("processType");
		fNtColIdStep[6] = analysisManager->CreateNtupleIColumn("systemID");
		fNtColIdStep[7] = analysisManager->CreateNtupleIColumn("cryNumber");
		fNtColIdStep[8] = analysisManager->CreateNtupleIColumn("detNumber");
		fNtColIdStep[9] = analysisManager->CreateNtupleDColumn("depEnergy");
		fNtColIdStep[10] = analysisManager->CreateNtupleDColumn("posx");
		fNtColIdStep[11] = analysisManager->CreateNtupleDColumn("posy");
		fNtColIdStep[12] = analysisManager->CreateNtupleDColumn("posz");
		fNtColIdStep[13] = analysisManager->CreateNtupleDColumn("time");
		fNtColIdStep[14] = analysisManager->CreateNtupleIColumn("targetZ");
		analysisManager->FinishNtuple();
	}

	fFactoryOn = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save() {
	if(fFactoryOn) {
		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		analysisManager->Write();
		analysisManager->CloseFile();

		delete G4AnalysisManager::Instance();
		fFactoryOn = false;
	}
}

void HistoManager::FillHitNtuple(G4int eventNumber) {
	FillHitNtuple(eventNumber, -1, -1, -1, -1, -1, -1, -1, -1, -1., -1., -1., -1., -1., -1);
}

void HistoManager::FillHitNtuple(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ) {
	if(fHitTrackerBool) {
		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		analysisManager->FillNtupleIColumn(fNtColIdHit[0], eventNumber);
		analysisManager->FillNtupleIColumn(fNtColIdHit[1], trackID);
		analysisManager->FillNtupleIColumn(fNtColIdHit[2], parentID);
		analysisManager->FillNtupleIColumn(fNtColIdHit[3], stepNumber);
		analysisManager->FillNtupleIColumn(fNtColIdHit[4], particleType);
		analysisManager->FillNtupleIColumn(fNtColIdHit[5], processType);
		analysisManager->FillNtupleIColumn(fNtColIdHit[6], systemID);
		analysisManager->FillNtupleIColumn(fNtColIdHit[7], cryNumber);
		analysisManager->FillNtupleIColumn(fNtColIdHit[8], detNumber);
		analysisManager->FillNtupleDColumn(fNtColIdHit[9], depEnergy);
		analysisManager->FillNtupleDColumn(fNtColIdHit[10], posx);
		analysisManager->FillNtupleDColumn(fNtColIdHit[11], posy);
		analysisManager->FillNtupleDColumn(fNtColIdHit[12], posz);
		analysisManager->FillNtupleDColumn(fNtColIdHit[13], time);
		analysisManager->FillNtupleIColumn(fNtColIdHit[14], targetZ);

		if(fRecordGun) {
			analysisManager->FillNtupleDColumn(fNtColIdHit[fFirstRecordingId+0], fBeamEnergy/CLHEP::keV);
			analysisManager->FillNtupleDColumn(fNtColIdHit[fFirstRecordingId+1], fBeamTheta);
			analysisManager->FillNtupleDColumn(fNtColIdHit[fFirstRecordingId+2], fBeamPhi);
			analysisManager->FillNtupleDColumn(fNtColIdHit[fFirstRecordingId+3], fBeamPos.x());
			analysisManager->FillNtupleDColumn(fNtColIdHit[fFirstRecordingId+4], fBeamPos.y());
			analysisManager->FillNtupleDColumn(fNtColIdHit[fFirstRecordingId+5], fBeamPos.z());
		}

		analysisManager->AddNtupleRow();
	}
}

void HistoManager::FillStepNtuple(G4int eventNumber) {
	FillStepNtuple(eventNumber, -1, -1, -1, -1, -1, -1, -1, -1, -1., -1., -1., -1., -1., -1);
}

void HistoManager::FillStepNtuple(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ) {
	if(fStepTrackerBool) {
		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		analysisManager->FillNtupleIColumn(fNtColIdStep[0], eventNumber);
		analysisManager->FillNtupleIColumn(fNtColIdStep[1], trackID);
		analysisManager->FillNtupleIColumn(fNtColIdStep[2], parentID);
		analysisManager->FillNtupleIColumn(fNtColIdStep[3], stepNumber);
		analysisManager->FillNtupleIColumn(fNtColIdStep[4], particleType);
		analysisManager->FillNtupleIColumn(fNtColIdStep[5], processType);
		analysisManager->FillNtupleIColumn(fNtColIdStep[6], systemID);
		analysisManager->FillNtupleIColumn(fNtColIdStep[7], cryNumber);
		analysisManager->FillNtupleIColumn(fNtColIdStep[8], detNumber);
		analysisManager->FillNtupleDColumn(fNtColIdStep[9], depEnergy);
		analysisManager->FillNtupleDColumn(fNtColIdStep[10], posx);
		analysisManager->FillNtupleDColumn(fNtColIdStep[11], posy);
		analysisManager->FillNtupleDColumn(fNtColIdStep[12], posz);
		analysisManager->FillNtupleDColumn(fNtColIdStep[13], time);
		analysisManager->FillNtupleIColumn(fNtColIdStep[14], targetZ);
		analysisManager->AddNtupleRow();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String HistoManager::G4intToG4String(G4int value) {
	G4String theString;
	std::stringstream out;
	out<<value;
	theString = out.str();
	return theString;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::MakeHistogram(G4AnalysisManager* analysisManager, G4String name,  G4String title, G4double xmin, G4double xmax, G4int nbins) {
	fMakeHistogramIndex++;
	if(fMakeHistogramIndex >= MAXHISTO) {
		G4cout<<"---> Exceeded maximum number of histograms. Increase MAXHISTO in HistoManager.hh"<<G4endl;
		exit(1);
	}

	fHistId[fMakeHistogramIndex] = analysisManager->CreateH1(name, title, nbins, xmin, xmax);
	fHistPt[fMakeHistogramIndex] = analysisManager->GetH1(fHistId[fMakeHistogramIndex]);
}

void HistoManager::Make2DHistogram(G4AnalysisManager* analysisManager, const G4String& name, const G4String& title,
		G4int nxbins, G4double xmin, G4double xmax, 
		G4int nybins, G4double ymin, G4double ymax) {
	fMakeHistogramIndex++; //global in Histomanager so allowed in this scope unaltered
	if(fMakeHistogramIndex >= MAXHISTO) {
		G4cout<<"---> Exceeded maximum number of histograms. Increase MAXHISTO in HistoManager.hh"<<G4endl;
		exit(1);
	}

	fHistId[fMakeHistogramIndex] = analysisManager->CreateH2(name, title, nxbins, xmin, xmax, nybins, ymin, ymax);
	fHistPt2[fMakeHistogramIndex] = analysisManager->GetH2(fHistId[fMakeHistogramIndex]);
}

void HistoManager::FillHistogram(G4int ih, G4double xbin, G4double weight) {
	if(ih >= MAXHISTO) {
		G4cout<<"---> warning from HistoManager::FillHistogram() : histo "<<ih
			<<"does not exist; xbin= "<<xbin<<" weight= "<<weight<<G4endl;
		return;
	}
	if(fHistPt[ih]) fHistPt[ih]->fill(xbin, weight);
}

void HistoManager::Fill2DHistogram(G4int ih, G4double xbin, G4double ybin, G4double weight) {
	if(ih >= MAXHISTO) {
		G4cout<<"---> warning from HistoManager::Fill2DHistogram() : histo "<<ih
			<< "does not exist; xbin= "<<xbin<<" ybin= "<< ybin<<" weight= "<<weight<<G4endl;
		return;
	}
	if(fHistPt2[ih]) fHistPt2[ih]->fill(xbin, ybin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

