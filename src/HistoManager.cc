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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
    fileName[0] = "g4out";
    factoryOn = false;

    // Only fill one NTuple at a time. If stepTrackerBool is true, then hitTrackerBool should be false, or vise-versa.
    // There is no need to have the hit NTuple and the step NTuple.
    hitTrackerBool = true;
    stepTrackerBool = false;

    makeHistoIndex = 0;

    // histograms
    for (G4int k=0; k<MAXHISTO; k++) {
        fHistId[k] = 0;
        fHistPt[k] = 0;
    }
    // ntuple
    for (G4int k=0; k<MAXNTCOL; k++) {
//        fNtColId[k] = 0;
        fNtColIdHit[k] = 0;
        fNtColIdStep[k] = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{
    G4String filename;
    G4String title;
    G4String detString;
    G4String cryString;
    G4double xmin;
    G4double xmax;
    G4int    nbins;
    //  G4int	 	 detStringNum ;

    // determine the maximum number of file names we will need
    // this corresponds to the maximum number of detectors.
    //	if( ( MAXNUMDET >= MAXNUMDETGRIFFIN ) && MAXNUMDET != 0 )
    //		detStringNum = MAXNUMDET ;
    //	else
    //		detStringNum = MAXNUMDETGRIFFIN ;
    //
    //	G4String detString[detStringNum] ;
    //
    //	// This should reduce calls to G4intToG4String in the future.
    //	for( G4int i = 0 ; i < detStringNum ; i++ )
    //		detString[i] = G4intToG4String(i);



    // Create or get analysis manager
    // The choice of analysis technology is done via selection of a namespace
    // in HistoManager.hh
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel(2);
    G4String extension = analysisManager->GetFileType();
    fileName[1] = fileName[0] + "." + extension;

    // Create directories
    analysisManager->SetHistoDirectoryName("histo");
    analysisManager->SetNtupleDirectoryName("ntuple");

    // Open an output file
    G4bool fileOpen = analysisManager->OpenFile(fileName[0]);
    if (!fileOpen) {
        G4cout << "\n---> HistoManager::book(): cannot open " << fileName[1]
               << G4endl;
        return;
    }

    // create selected histograms
    analysisManager->SetFirstHistoId(1);
    filename  = "astats_particle_type_in_each_step";
    title     = "Particle Creation";
    nbins     = 20;
    xmin      = 0.;
    xmax      = 20.;
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "astats_particle_type_in_each_event";
    title     = "Number of Particle Types in Event";
    nbins     = 100;
    xmin      = 0.;
    xmax      = 100.;
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    if(WRITEEKINHISTOS) {
        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);
            filename  = "gridcell_electron_ekin_det" + detString;
            title     = "EKin in cell (keV)";
            nbins     = EKINNBINS;
            xmin      = EKINXMIN;
            xmax      = EKINXMAX;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }
    }

    if(WRITETRACKLHISTOS) {
        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);
            filename  = "gridcell_electron_trackl_det" + detString;
            title     = "Trackl in cell (keV)";
            nbins     = TRACKLNBINS;
            xmin      = TRACKLXMIN;
            xmax      = TRACKLXMAX;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }
    }

    if(WRITEEKINHISTOS) {
        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);
            filename  = "gridcell_gamma_ekin_det" + detString;
            title     = "EKin in cell (keV)";
            nbins     = EKINNBINS;
            xmin      = EKINXMIN;
            xmax      = EKINXMAX;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }
    }

    if(WRITETRACKLHISTOS) {
        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);
            filename  = "gridcell_gamma_trackl_det" + detString;
            title     = "Trackl in cell (keV)";
            nbins     = TRACKLNBINS;
            xmin      = TRACKLXMIN;
            xmax      = TRACKLXMAX;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }
    }

    if(WRITEEKINHISTOS) {
        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);
            filename  = "gridcell_neutron_ekin_det" + detString;
            title     = "EKin in cell (keV)";
            nbins     = EKINNBINS;
            xmin      = EKINXMIN;
            xmax      = EKINXMAX;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }
    }

    if(WRITETRACKLHISTOS) {
        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);
            filename  = "gridcell_neutron_trackl_det" + detString;
            title     = "Trackl in cell (keV)";
            nbins     = TRACKLNBINS;
            xmin      = TRACKLXMIN;
            xmax      = TRACKLXMAX;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }
    }


    if(WRITEEDEPHISTOS)
    {
        // Variables and title used for all detectors
        nbins     = EDEPNBINS;
        xmin      = EDEPXMIN;
        xmax      = EDEPXMAX;
        title     = "Edep in crystal (keV)";


        // Griffin Suppressors
        filename  = "griffin_crystal_sup_edep";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "griffin_crystal_sup_edep_cry";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "griffin_crystal_sup_edep_sum";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
            detString = G4intToG4String(i);

            filename  = "griffin_crystal_sup_edep_det" + detString ;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }

        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
                detString = G4intToG4String(i);
                cryString = G4intToG4String(j);

                filename  = "griffin_crystal_sup_edep_det" + detString + "_cry" + cryString;
                MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
            }
        }

        // Griffin Crystal Unsuppressed
        filename  = "griffin_crystal_unsup_edep";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "griffin_crystal_unsup_edep_cry";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "griffin_crystal_unsup_edep_sum";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
            detString = G4intToG4String(i);

            filename  = "griffin_crystal_unsup_edep_det" + detString;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }

        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
                detString = G4intToG4String(i);
                cryString = G4intToG4String(j);

                filename  = "griffin_crystal_unsup_edep_det" + detString + "_cry" + cryString;
                MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
            }
        }

        // Brilliance Detector
        filename  = "labr_crystal_edep";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "labr_crystal_edep_sum";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);

            filename  = "labr_crystal_edep_det" + detString;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }

        // Ancillary BGO Detector
        filename  = "ancillary_bgo_crystal_edep";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "ancillary_bgo_crystal_edep_sum";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);

            filename  = "ancillary_bgo_crystal_edep_det" + detString;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }


        // Sodium Iodide detector
        filename  = "sodiumIodide_crystal_edep";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "sodiumIodide_crystal_edep_sum";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);

            filename  = "sodiumIodide_crystal_edep_det" + detString;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }

        // Sceptar detector
        filename  = "sceptar_edep";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "sceptar_edep_sum";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);

            filename  = "sceptar_edep_det" + detString;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }

        // 8pi detector
        filename  = "Eightpi_crystal_edep";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "Eightpi_crystal_edep_sum";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);

            filename  = "Eightpi_crystal_edep_det" + detString;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }

        // spice detector
        filename  = "spice_crystal_edep";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "spice_crystal_edep_sum";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);

            filename  = "spice_crystal_edep_det" + detString;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }

        // paces detector
        filename  = "paces_crystal_edep";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        filename  = "paces_crystal_edep_sum";
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

        for (G4int i=0; i < MAXNUMDET; i++) {
            detString = G4intToG4String(i);

            filename  = "paces_crystal_edep_det" + detString;
            MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
        }
    }


    ///////////////////////////////////////////////////////////////////
    // Create 1 ntuple
    if(hitTrackerBool) {
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
        analysisManager->FinishNtuple();
    }


    if(stepTrackerBool) {
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
        analysisManager->FinishNtuple();
    }

    factoryOn = true;
    G4cout << "\n----> Histogram Tree is opened in " << fileName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{
    if (factoryOn) {
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        analysisManager->Write();
        analysisManager->CloseFile();
        G4cout << "\n----> Histogram Tree is saved in " << fileName[1] << G4endl;

        delete G4AnalysisManager::Instance();
        factoryOn = false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::MakeHisto(G4AnalysisManager* analysisManager, G4String filename,  G4String title, G4double xmin, G4double xmax, G4int nbins)
{
    makeHistoIndex++;
    if (makeHistoIndex >= MAXHISTO) {
        G4cout << "---> Exceeded maximum number of histograms. Increase MAXHISTO in HistoManager.hh" << G4endl;
        exit(1);
    }
    fHistId[makeHistoIndex] = analysisManager->CreateH1(filename, title, nbins, xmin, xmax);
    fHistPt[makeHistoIndex] = analysisManager->GetH1(fHistId[makeHistoIndex]);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
    if (ih > MAXHISTO) {
        G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
               << "does note xist; xbin= " << xbin << " w= " << weight << G4endl;
        return;
    }
    if (fHistPt[ih]) fHistPt[ih]->fill(xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
    if (ih >= MAXHISTO) {
        G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
               << "  fac= " << fac << G4endl;
        return;
    }
    if (fHistPt[ih]) fHistPt[ih]->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void HistoManager::FillNtuple(G4double eventNumber, G4double stepNumber, G4double cryNumber, G4double detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time)
//{
//    if(stepTrackerBool) {
//        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//        analysisManager->FillNtupleDColumn(fNtColId[0], eventNumber);
//        analysisManager->FillNtupleDColumn(fNtColId[1], stepNumber);
//        analysisManager->FillNtupleDColumn(fNtColId[2], cryNumber);
//        analysisManager->FillNtupleDColumn(fNtColId[3], detNumber);
//        analysisManager->FillNtupleDColumn(fNtColId[4], depEnergy);
//        analysisManager->FillNtupleDColumn(fNtColId[5], posx);
//        analysisManager->FillNtupleDColumn(fNtColId[6], posy);
//        analysisManager->FillNtupleDColumn(fNtColId[7], posz);
//        analysisManager->FillNtupleDColumn(fNtColId[8], time);
//        analysisManager->FillNtupleDColumn(fNtColId[9], 0);
//        analysisManager->AddNtupleRow();
//    }
//}


void HistoManager::FillHitNtuple(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time)
{
    if(hitTrackerBool) {
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
        analysisManager->AddNtupleRow();
    }
}

void HistoManager::FillStepNtuple(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time)
{
    if(stepTrackerBool) {
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
        analysisManager->AddNtupleRow();
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String HistoManager::G4intToG4String(G4int value)
{
    G4String theString;
    std::stringstream out;
    out << value;
    theString = out.str();
    return theString;
}
