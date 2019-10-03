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
	fFileName[0] = "g4out";
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book() {
	// Create or get analysis manager
	// The choice of analysis technology is done via selection of a namespace
	// in HistoManager.hh
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	//analysisManager->SetVerboseLevel(2);
	G4String extension = analysisManager->GetFileType();
	fFileName[1] = fFileName[0] + "." + extension; // creating root output file in build folder

	// Create directories
	// Open an output file
	G4bool fileOpen = analysisManager->OpenFile(fFileName[0]); 
	if(!fileOpen) {
		G4cout<<"---> HistoManager::book(): cannot open "<<fFileName[1]<<G4endl;
		return;
	}

	///////////////////////////////////////////////////////////////////
	// Create 1 ntuple
	if(fHitTrackerBool) {
		analysisManager->CreateNtuple("ntuple", "HitTracker");
		fNtColIdHit[0] = analysisManager->CreateNtupleIColumn("eventNumber");
		fNtColIdHit[1] = analysisManager->CreateNtupleIColumn("trackID");
		fNtColIdHit[2] = analysisManager->CreateNtupleIColumn("parentID");
		fNtColIdHit[3] = analysisManager->CreateNtupleIColumn("stepNumber");
		fNtColIdHit[4] = analysisManager->CreateNtupleIColumn("particleType");
		fNtColIdHit[5] = analysisManager->CreateNtupleIColumn("processType");
		fNtColIdHit[6] = analysisManager->CreateNtupleIColumn("systemID");
		fNtColIdHit[7] = analysisManager->CreateNtupleIColumn("cryNumber");
		fNtColIdHit[8] = analysisManager->CreateNtupleIColumn("detNumber");
		fNtColIdHit[9] = analysisManager->CreateNtupleDColumn("depEnergy");
		fNtColIdHit[10] = analysisManager->CreateNtupleDColumn("posx");
		fNtColIdHit[11] = analysisManager->CreateNtupleDColumn("posy");
		fNtColIdHit[12] = analysisManager->CreateNtupleDColumn("posz");
		fNtColIdHit[13] = analysisManager->CreateNtupleDColumn("time");
		fNtColIdHit[14] = analysisManager->CreateNtupleIColumn("targetZ");
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

	if(fDetectorConstruction->Spice()) {
		BookSpiceHistograms();
	}

	fFactoryOn = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save() {
	if(fFactoryOn) {
		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		analysisManager->Write();
		analysisManager->CloseFile();
		//G4cout<<"----> Histogram Tree is saved in "<<fFileName[1]<<G4endl;

		delete G4AnalysisManager::Instance();
		fFactoryOn = false;
	}
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
		analysisManager->AddNtupleRow();
	}
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

void HistoManager::BookSpiceHistograms() {
	// create selected histograms
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	G4String name;
	G4String title;
	G4String detString;
	G4String cryString;
	G4String segString;//for spice segment for-loop
	G4double xmin;
	G4double xmax;
	G4int    nbins;

	analysisManager->SetFirstHistoId(1);
	name  = "astats_particle_type_in_each_step";
	title     = "Particle Creation";
	nbins     = 20;
	xmin      = 0.;
	xmax      = 20.;
	MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

	name  = "astats_particle_type_in_each_event";
	title     = "Number of Particle Types in Event";
	nbins     = 100;
	xmin      = 0.;
	xmax      = 100.;
	MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

	if(fDetectorConstruction->GridCell() && WRITEEKINHISTOS) {
		for(G4int i=0; i < MAXNUMDET; i++) {
			detString = G4intToG4String(i);
			name  = "gridcell_electron_ekin_det" + detString;
			title     = "EKin in cell (keV)";
			nbins     = EKINNBINS;
			xmin      = EKINXMIN;
			xmax      = EKINXMAX;
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
		}
	}

	if(fDetectorConstruction->GridCell() && WRITETRACKLHISTOS) {
		for(G4int i=0; i < MAXNUMDET; i++) {
			detString = G4intToG4String(i);
			name  = "gridcell_electron_trackl_det" + detString;
			title     = "Trackl in cell (keV)";
			nbins     = TRACKLNBINS;
			xmin      = TRACKLXMIN;
			xmax      = TRACKLXMAX;
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
		}
	}

	if(fDetectorConstruction->GridCell() && WRITEEKINHISTOS) {
		for(G4int i=0; i < MAXNUMDET; i++) {
			detString = G4intToG4String(i);
			name  = "gridcell_gamma_ekin_det" + detString;
			title     = "EKin in cell (keV)";
			nbins     = EKINNBINS;
			xmin      = EKINXMIN;
			xmax      = EKINXMAX;
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
		}
	}

	if(fDetectorConstruction->GridCell() && WRITETRACKLHISTOS) {
		for(G4int i=0; i < MAXNUMDET; i++) {
			detString = G4intToG4String(i);
			name  = "gridcell_gamma_trackl_det" + detString;
			title     = "Trackl in cell (keV)";
			nbins     = TRACKLNBINS;
			xmin      = TRACKLXMIN;
			xmax      = TRACKLXMAX;
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
		}
	}

	if(fDetectorConstruction->GridCell() && WRITEEKINHISTOS) {
		for(G4int i=0; i < MAXNUMDET; i++) {
			detString = G4intToG4String(i);
			name  = "gridcell_neutron_ekin_det" + detString;
			title     = "EKin in cell (keV)";
			nbins     = EKINNBINS;
			xmin      = EKINXMIN;
			xmax      = EKINXMAX;
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
		}
	}

	if(fDetectorConstruction->GridCell() && WRITETRACKLHISTOS) {
		for(G4int i=0; i < MAXNUMDET; i++) {
			detString = G4intToG4String(i);
			name  = "gridcell_neutron_trackl_det" + detString;
			title     = "Trackl in cell (keV)";
			nbins     = TRACKLNBINS;
			xmin      = TRACKLXMIN;
			xmax      = TRACKLXMAX;
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
		}
	}

	if(WRITEEDEPHISTOS) {
		// Variables and title used for all detectors
		nbins     = EDEPNBINS;
		xmin      = EDEPXMIN;
		xmax      = EDEPXMAX;
		title     = "Edep in crystal (keV)";
		if(fDetectorConstruction->Griffin()) {
			// Griffin Suppressors
			name  = "griffin_crystal_sup_edep";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			name  = "griffin_crystal_sup_edep_cry";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			name  = "griffin_crystal_sup_edep_sum";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			for(G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
				detString = G4intToG4String(i);

				name  = "griffin_crystal_sup_edep_det" + detString ;
				MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
			}

			for(G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
				for(G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
					detString = G4intToG4String(i);
					cryString = G4intToG4String(j);

					name  = "griffin_crystal_sup_edep_det" + detString + "_cry" + cryString;
					MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
				}
			}

			// Griffin Crystal Unsuppressed
			name  = "griffin_crystal_unsup_edep";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			name  = "griffin_crystal_unsup_edep_cry";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			name  = "griffin_crystal_unsup_edep_sum";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			for(G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
				detString = G4intToG4String(i);

				name  = "griffin_crystal_unsup_edep_det" + detString;
				MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
			}

			for(G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
				for(G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
					detString = G4intToG4String(i);
					cryString = G4intToG4String(j);

					name  = "griffin_crystal_unsup_edep_det" + detString + "_cry" + cryString;
					MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
				}
			}
		}//if(fGriffin)

		if(fDetectorConstruction->LaBr()) {
			// Brilliance Detector
			name  = "labr_crystal_edep";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			name  = "labr_crystal_edep_sum";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			for(G4int i=0; i < MAXNUMDET; i++) {
				detString = G4intToG4String(i);

				name  = "labr_crystal_edep_det" + detString;
				MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
			}
		}//if(LaBr)

		if(fDetectorConstruction->AncBgo()) {
			// Ancillary BGO Detector
			name  = "ancillary_bgo_crystal_edep";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			name  = "ancillary_bgo_crystal_edep_sum";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			for(G4int i=0; i < MAXNUMDET; i++) {
				detString = G4intToG4String(i);

				name  = "ancillary_bgo_crystal_edep_det" + detString;
				MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
			}
		}//if(fAncBgo)

		if(fDetectorConstruction->NaI()) {
			// Sodium Iodide detector
			name  = "sodiumIodide_crystal_edep";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			name  = "sodiumIodide_crystal_edep_sum";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			for(G4int i=0; i < MAXNUMDET; i++) {
				detString = G4intToG4String(i);

				name  = "sodiumIodide_crystal_edep_det" + detString;
				MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
			}
		}//if(fNaI)

		if(fDetectorConstruction->Sceptar()) {
			// Sceptar detector
			name  = "sceptar_edep";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			name  = "sceptar_edep_sum";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			for(G4int i=0; i < MAXNUMDET; i++) {
				detString = G4intToG4String(i);

				name  = "sceptar_edep_det" + detString;
				MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
			}
		}//if(fSceptar)

		if(fDetectorConstruction->EightPi()) {
			// 8pi detector
			name  = "Eightpi_crystal_edep";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			name  = "Eightpi_crystal_edep_sum";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

			for(G4int i=0; i < MAXNUMDET; i++) {
				detString = G4intToG4String(i);

				name  = "Eightpi_crystal_edep_det" + detString;
				MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
			}
		}//if(fEightPi)

		if(fDetectorConstruction->Spice()) {// spice detector
			nbins     = EDEPNBINSSPICE;
			xmin      = EDEPXMINSPICE;
			xmax      = EDEPXMAXSPICE;

			title     = "Edep in crystal (keV)";
			name  = "FullEdep"; //edep = energy deposition
			MakeHistogram(analysisManager, name,  title, xmin*1000., xmax*1000., nbins);
			fSpiceHistNumbers[0] = fMakeHistogramIndex; //assigns array item the histo index, for reference by FillHisto()

			name  = "Edep_Addback";// summed
			MakeHistogram(analysisManager, name,  title, xmin*1000., xmax*1000., nbins);
			fSpiceHistNumbers[1] = fMakeHistogramIndex;

			name = "BeamEnergy";
			title ="Energy of Input Beam";
			MakeHistogram(analysisManager, name, title, xmin, xmax, nbins);
			fAngleDistro[0] = fMakeHistogramIndex;

			name = "AllSegEnergies";
			title = "Energy vs. Segment";
			Make2DHistogram(analysisManager, name, title,
					120, 0., 120., nbins, 0., xmax*1000.);
			fAngleDistro[1] = fMakeHistogramIndex;

			name  = "x-y";
			title     = "X-Y 2D distribution";
			nbins = 100;
			xmin      = -6.;
			xmax      = 6.;
			Make2DHistogram(analysisManager, name,  title, nbins, xmin, xmax, nbins, xmin, xmax);
			fAngleDistro[2] = fMakeHistogramIndex;

			title  = "z-distribution";
			name     = "z-distro";
			MakeHistogram(analysisManager, name,  title, -0.1, 0.1, 1000);
			fAngleDistro[3] = fMakeHistogramIndex;

			for(G4int ring=0; ring < MAXNUMDETSPICE; ring++){
				for(G4int seg=0; seg < 3; seg++) {

					title  = "AngR" + std::to_string(ring) + "S" +std::to_string(seg);
					name     = "AngR" + std::to_string(ring) + "S" +std::to_string(seg);
					nbins     = 200;
					xmin      = 0.;//explicitly defined for clarity
					xmax      = CLHEP::pi;
					Make2DHistogram(analysisManager, name, title,
							nbins/2, xmin, xmax, nbins, -xmax, xmax);

					fSpiceAngleHists[ring*MAXNUMSEGSPICE+seg] = fMakeHistogramIndex;
				}
			}
		}
		if(fDetectorConstruction->Paces()) {// paces detector
			name  = "paces_crystal_edep";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
			fPacesHistNumbers[0] = fMakeHistogramIndex;

			name  = "paces_crystal_edep_sum";
			MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);
			fPacesHistNumbers[1] = fMakeHistogramIndex;

			for(G4int i=0; i < MAXNUMDETPACES; i++) {//[MAXNUMDET];
				detString = G4intToG4String(i);

				name  = "paces_crystal_edep_det" + detString;
				MakeHistogram(analysisManager, name,  title, xmin, xmax, nbins);

				fPacesHistNumbers[i+2] = fMakeHistogramIndex;
			}
		}//if(fPaces)
	}//if(WRITEEDEPHISTOS)
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

