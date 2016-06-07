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
/// \file analysis/AnaEx01/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
//
// $Id: HistoManager.hh 74272 2013-10-02 14:48:50Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"

#include "g4root.hh"
////#include "g4xml.hh"
////#include "g4hbook.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

const G4bool WRITEEKINHISTOS    = true;
const G4bool WRITEEDEPHISTOS    = true;
const G4bool WRITETRACKLHISTOS  = true;

const G4int MAXHISTO            = 500;
const G4int MAXNTCOL            = 15;
const G4int MAXNUMDET           = 20;
const G4int MAXNUMDETGRIFFIN    = 16;
const G4int MAXNUMCRYGRIFFIN    = 4;
const G4int NUMPARTICLETYPES    = 20;

// ekin histo properties    ///////////////////////
const G4int     EKINNBINS  = 10000;
const G4double  EKINXMIN   = 0.5*keV;
const G4double  EKINXMAX   = 10000.5*keV;

// edep histo properties    ///////////////////////
const G4int     EDEPNBINS  = 10000;
const G4double  EDEPXMIN   = 0.5*keV;
const G4double  EDEPXMAX   = 10000.5*keV;

// trackl histo properties  ///////////////////////
const G4int     TRACKLNBINS = 5000;
const G4double  TRACKLXMIN  = 0.5*mm;
const G4double  TRACKLXMAX  = 5000.5*mm;

///////////////////////////////////////////////////

const G4double MINENERGYTHRES   = 0.001*keV;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager {
public:
    static HistoManager& Instance() {
		static HistoManager histManager;
		return histManager;
	 }

    void Book();
    void Save();

    void MakeHisto(G4AnalysisManager* analysisManager, G4String filename,  G4String title, G4double xmin, G4double xmax, G4int nbins);
    void FillHisto(G4int ih, G4double e, G4double weight = 1.0);
    void Normalize(G4int id, G4double fac);

//    void FillNtuple(G4double eventNumber, G4double stepNumber, G4double cryNumber, G4double detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time);
    void FillHitNtuple(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ);
    void FillStepNtuple(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ);

    void PrintStatistic();

    G4bool GetStepTrackerBool() { return fStepTrackerBool;};
    G4bool GetHitTrackerBool()  { return fHitTrackerBool;};

	 // setter
	 void GridCell(G4bool val) { fGridCell = val; }
 	 void Griffin (G4bool val) { fGriffin  = val; }
	 void LaBr 	  (G4bool val) { fLaBr     = val; }
	 void AncBgo  (G4bool val) { fAncBgo   = val; }
	 void NaI 	  (G4bool val) { fNaI      = val; }
	 void Sceptar (G4bool val) { fSceptar  = val; }
	 void EightPi (G4bool val) { fEightPi  = val; }
	 void Spice   (G4bool val) { fSpice    = val; }
	 void Paces   (G4bool val) { fPaces    = val; }  
	 void Descant (G4bool val) { fDescant  = val; }  
	 void Testcan (G4bool val) { fTestcan  = val; }  
    // getter
	 G4bool GridCell() { return fGridCell; }
 	 G4bool Griffin () { return fGriffin;  }
	 G4bool LaBr 	 () { return fLaBr;     }
	 G4bool AncBgo  () { return fAncBgo;   }
	 G4bool NaI 	 () { return fNaI;      }
	 G4bool Sceptar () { return fSceptar;  }
	 G4bool EightPi () { return fEightPi;  }
	 G4bool Spice 	 () { return fSpice;    }
	 G4bool Paces   () { return fPaces;    }  
	 G4bool Descant () { return fDescant;  }  
	 G4bool Testcan () { return fTestcan;  }  

private:
    HistoManager();
    ~HistoManager();
    // for C++11 we could instead use public delete of these methods
    HistoManager(HistoManager const&);   // don't implement
    void operator=(HistoManager const&); // don't implement

    static HistoManager* fHistoManager;
    G4String G4intToG4String(G4int value);

    G4bool        fFactoryOn;
    G4int         fMakeHistoIndex;

    G4String      fFileName[2];

    G4int         fHistId[MAXHISTO];
    G4AnaH1*      fHistPt[MAXHISTO];

    G4int         fNtColId[MAXNTCOL];
    G4int         fNtColIdHit[MAXNTCOL];
    G4int         fNtColIdStep[MAXNTCOL];

    G4bool fStepTrackerBool;
    G4bool fHitTrackerBool;

	 //booleans which control which histograms are created (these are set by the detector construction)
	 G4bool fGridCell;
 	 G4bool fGriffin;
	 G4bool fLaBr;
	 G4bool fAncBgo;
	 G4bool fNaI;
	 G4bool fSceptar;
	 G4bool fEightPi;
	 G4bool fSpice;
	 G4bool fPaces;
	 G4bool fDescant;
	 G4bool fTestcan;
};

enum HISTONAME
{
    null_name = 0,
    astats_particle_type_in_each_step,
    astats_particle_type_in_each_event,
    gridcell_electron_ekin_det0,
    gridcell_electron_ekin_det1,
    gridcell_electron_ekin_det2,
    gridcell_electron_ekin_det3,
    gridcell_electron_ekin_det4,
    gridcell_electron_ekin_det5,
    gridcell_electron_ekin_det6,
    gridcell_electron_ekin_det7,
    gridcell_electron_ekin_det8,
    gridcell_electron_ekin_det9,
    gridcell_electron_ekin_det10,
    gridcell_electron_ekin_det11,
    gridcell_electron_ekin_det12,
    gridcell_electron_ekin_det13,
    gridcell_electron_ekin_det14,
    gridcell_electron_ekin_det15,
    gridcell_electron_ekin_det16,
    gridcell_electron_ekin_det17,
    gridcell_electron_ekin_det18,
    gridcell_electron_ekin_det19,
    gridcell_electron_trackl_det0,
    gridcell_electron_trackl_det1,
    gridcell_electron_trackl_det2,
    gridcell_electron_trackl_det3,
    gridcell_electron_trackl_det4,
    gridcell_electron_trackl_det5,
    gridcell_electron_trackl_det6,
    gridcell_electron_trackl_det7,
    gridcell_electron_trackl_det8,
    gridcell_electron_trackl_det9,
    gridcell_electron_trackl_det10,
    gridcell_electron_trackl_det11,
    gridcell_electron_trackl_det12,
    gridcell_electron_trackl_det13,
    gridcell_electron_trackl_det14,
    gridcell_electron_trackl_det15,
    gridcell_electron_trackl_det16,
    gridcell_electron_trackl_det17,
    gridcell_electron_trackl_det18,
    gridcell_electron_trackl_det19,
    gridcell_gamma_ekin_det0,
    gridcell_gamma_ekin_det1,
    gridcell_gamma_ekin_det2,
    gridcell_gamma_ekin_det3,
    gridcell_gamma_ekin_det4,
    gridcell_gamma_ekin_det5,
    gridcell_gamma_ekin_det6,
    gridcell_gamma_ekin_det7,
    gridcell_gamma_ekin_det8,
    gridcell_gamma_ekin_det9,
    gridcell_gamma_ekin_det10,
    gridcell_gamma_ekin_det11,
    gridcell_gamma_ekin_det12,
    gridcell_gamma_ekin_det13,
    gridcell_gamma_ekin_det14,
    gridcell_gamma_ekin_det15,
    gridcell_gamma_ekin_det16,
    gridcell_gamma_ekin_det17,
    gridcell_gamma_ekin_det18,
    gridcell_gamma_ekin_det19,
    gridcell_gamma_trackl_det0,
    gridcell_gamma_trackl_det1,
    gridcell_gamma_trackl_det2,
    gridcell_gamma_trackl_det3,
    gridcell_gamma_trackl_det4,
    gridcell_gamma_trackl_det5,
    gridcell_gamma_trackl_det6,
    gridcell_gamma_trackl_det7,
    gridcell_gamma_trackl_det8,
    gridcell_gamma_trackl_det9,
    gridcell_gamma_trackl_det10,
    gridcell_gamma_trackl_det11,
    gridcell_gamma_trackl_det12,
    gridcell_gamma_trackl_det13,
    gridcell_gamma_trackl_det14,
    gridcell_gamma_trackl_det15,
    gridcell_gamma_trackl_det16,
    gridcell_gamma_trackl_det17,
    gridcell_gamma_trackl_det18,
    gridcell_gamma_trackl_det19,
    gridcell_neutron_ekin_det0,
    gridcell_neutron_ekin_det1,
    gridcell_neutron_ekin_det2,
    gridcell_neutron_ekin_det3,
    gridcell_neutron_ekin_det4,
    gridcell_neutron_ekin_det5,
    gridcell_neutron_ekin_det6,
    gridcell_neutron_ekin_det7,
    gridcell_neutron_ekin_det8,
    gridcell_neutron_ekin_det9,
    gridcell_neutron_ekin_det10,
    gridcell_neutron_ekin_det11,
    gridcell_neutron_ekin_det12,
    gridcell_neutron_ekin_det13,
    gridcell_neutron_ekin_det14,
    gridcell_neutron_ekin_det15,
    gridcell_neutron_ekin_det16,
    gridcell_neutron_ekin_det17,
    gridcell_neutron_ekin_det18,
    gridcell_neutron_ekin_det19,
    gridcell_neutron_trackl_det0,
    gridcell_neutron_trackl_det1,
    gridcell_neutron_trackl_det2,
    gridcell_neutron_trackl_det3,
    gridcell_neutron_trackl_det4,
    gridcell_neutron_trackl_det5,
    gridcell_neutron_trackl_det6,
    gridcell_neutron_trackl_det7,
    gridcell_neutron_trackl_det8,
    gridcell_neutron_trackl_det9,
    gridcell_neutron_trackl_det10,
    gridcell_neutron_trackl_det11,
    gridcell_neutron_trackl_det12,
    gridcell_neutron_trackl_det13,
    gridcell_neutron_trackl_det14,
    gridcell_neutron_trackl_det15,
    gridcell_neutron_trackl_det16,
    gridcell_neutron_trackl_det17,
    gridcell_neutron_trackl_det18,
    gridcell_neutron_trackl_det19,
    griffin_crystal_sup_edep,
    griffin_crystal_sup_edep_cry,
    griffin_crystal_sup_edep_sum,
    griffin_crystal_sup_edep_det0,
    griffin_crystal_sup_edep_det1,
    griffin_crystal_sup_edep_det2,
    griffin_crystal_sup_edep_det3,
    griffin_crystal_sup_edep_det4,
    griffin_crystal_sup_edep_det5,
    griffin_crystal_sup_edep_det6,
    griffin_crystal_sup_edep_det7,
    griffin_crystal_sup_edep_det8,
    griffin_crystal_sup_edep_det9,
    griffin_crystal_sup_edep_det10,
    griffin_crystal_sup_edep_det11,
    griffin_crystal_sup_edep_det12,
    griffin_crystal_sup_edep_det13,
    griffin_crystal_sup_edep_det14,
    griffin_crystal_sup_edep_det15,
    griffin_crystal_sup_edep_det0_cry0,
    griffin_crystal_sup_edep_det1_cry0,
    griffin_crystal_sup_edep_det2_cry0,
    griffin_crystal_sup_edep_det3_cry0,
    griffin_crystal_sup_edep_det4_cry0,
    griffin_crystal_sup_edep_det5_cry0,
    griffin_crystal_sup_edep_det6_cry0,
    griffin_crystal_sup_edep_det7_cry0,
    griffin_crystal_sup_edep_det8_cry0,
    griffin_crystal_sup_edep_det9_cry0,
    griffin_crystal_sup_edep_det10_cry0,
    griffin_crystal_sup_edep_det11_cry0,
    griffin_crystal_sup_edep_det12_cry0,
    griffin_crystal_sup_edep_det13_cry0,
    griffin_crystal_sup_edep_det14_cry0,
    griffin_crystal_sup_edep_det15_cry0,
    griffin_crystal_sup_edep_det0_cry1,
    griffin_crystal_sup_edep_det1_cry1,
    griffin_crystal_sup_edep_det2_cry1,
    griffin_crystal_sup_edep_det3_cry1,
    griffin_crystal_sup_edep_det4_cry1,
    griffin_crystal_sup_edep_det5_cry1,
    griffin_crystal_sup_edep_det6_cry1,
    griffin_crystal_sup_edep_det7_cry1,
    griffin_crystal_sup_edep_det8_cry1,
    griffin_crystal_sup_edep_det9_cry1,
    griffin_crystal_sup_edep_det10_cry1,
    griffin_crystal_sup_edep_det11_cry1,
    griffin_crystal_sup_edep_det12_cry1,
    griffin_crystal_sup_edep_det13_cry1,
    griffin_crystal_sup_edep_det14_cry1,
    griffin_crystal_sup_edep_det15_cry1,
    griffin_crystal_sup_edep_det0_cry2,
    griffin_crystal_sup_edep_det1_cry2,
    griffin_crystal_sup_edep_det2_cry2,
    griffin_crystal_sup_edep_det3_cry2,
    griffin_crystal_sup_edep_det4_cry2,
    griffin_crystal_sup_edep_det5_cry2,
    griffin_crystal_sup_edep_det6_cry2,
    griffin_crystal_sup_edep_det7_cry2,
    griffin_crystal_sup_edep_det8_cry2,
    griffin_crystal_sup_edep_det9_cry2,
    griffin_crystal_sup_edep_det10_cry2,
    griffin_crystal_sup_edep_det11_cry2,
    griffin_crystal_sup_edep_det12_cry2,
    griffin_crystal_sup_edep_det13_cry2,
    griffin_crystal_sup_edep_det14_cry2,
    griffin_crystal_sup_edep_det15_cry2,
    griffin_crystal_sup_edep_det0_cry3,
    griffin_crystal_sup_edep_det1_cry3,
    griffin_crystal_sup_edep_det2_cry3,
    griffin_crystal_sup_edep_det3_cry3,
    griffin_crystal_sup_edep_det4_cry3,
    griffin_crystal_sup_edep_det5_cry3,
    griffin_crystal_sup_edep_det6_cry3,
    griffin_crystal_sup_edep_det7_cry3,
    griffin_crystal_sup_edep_det8_cry3,
    griffin_crystal_sup_edep_det9_cry3,
    griffin_crystal_sup_edep_det10_cry3,
    griffin_crystal_sup_edep_det11_cry3,
    griffin_crystal_sup_edep_det12_cry3,
    griffin_crystal_sup_edep_det13_cry3,
    griffin_crystal_sup_edep_det14_cry3,
    griffin_crystal_sup_edep_det15_cry3,
    griffin_crystal_unsup_edep,
    griffin_crystal_unsup_edep_cry,
    griffin_crystal_unsup_edep_sum,
    griffin_crystal_unsup_edep_det0,
    griffin_crystal_unsup_edep_det1,
    griffin_crystal_unsup_edep_det2,
    griffin_crystal_unsup_edep_det3,
    griffin_crystal_unsup_edep_det4,
    griffin_crystal_unsup_edep_det5,
    griffin_crystal_unsup_edep_det6,
    griffin_crystal_unsup_edep_det7,
    griffin_crystal_unsup_edep_det8,
    griffin_crystal_unsup_edep_det9,
    griffin_crystal_unsup_edep_det10,
    griffin_crystal_unsup_edep_det11,
    griffin_crystal_unsup_edep_det12,
    griffin_crystal_unsup_edep_det13,
    griffin_crystal_unsup_edep_det14,
    griffin_crystal_unsup_edep_det15,
    griffin_crystal_unsup_edep_det0_cry0,
    griffin_crystal_unsup_edep_det1_cry0,
    griffin_crystal_unsup_edep_det2_cry0,
    griffin_crystal_unsup_edep_det3_cry0,
    griffin_crystal_unsup_edep_det4_cry0,
    griffin_crystal_unsup_edep_det5_cry0,
    griffin_crystal_unsup_edep_det6_cry0,
    griffin_crystal_unsup_edep_det7_cry0,
    griffin_crystal_unsup_edep_det8_cry0,
    griffin_crystal_unsup_edep_det9_cry0,
    griffin_crystal_unsup_edep_det10_cry0,
    griffin_crystal_unsup_edep_det11_cry0,
    griffin_crystal_unsup_edep_det12_cry0,
    griffin_crystal_unsup_edep_det13_cry0,
    griffin_crystal_unsup_edep_det14_cry0,
    griffin_crystal_unsup_edep_det15_cry0,
    griffin_crystal_unsup_edep_det0_cry1,
    griffin_crystal_unsup_edep_det1_cry1,
    griffin_crystal_unsup_edep_det2_cry1,
    griffin_crystal_unsup_edep_det3_cry1,
    griffin_crystal_unsup_edep_det4_cry1,
    griffin_crystal_unsup_edep_det5_cry1,
    griffin_crystal_unsup_edep_det6_cry1,
    griffin_crystal_unsup_edep_det7_cry1,
    griffin_crystal_unsup_edep_det8_cry1,
    griffin_crystal_unsup_edep_det9_cry1,
    griffin_crystal_unsup_edep_det10_cry1,
    griffin_crystal_unsup_edep_det11_cry1,
    griffin_crystal_unsup_edep_det12_cry1,
    griffin_crystal_unsup_edep_det13_cry1,
    griffin_crystal_unsup_edep_det14_cry1,
    griffin_crystal_unsup_edep_det15_cry1,
    griffin_crystal_unsup_edep_det0_cry2,
    griffin_crystal_unsup_edep_det1_cry2,
    griffin_crystal_unsup_edep_det2_cry2,
    griffin_crystal_unsup_edep_det3_cry2,
    griffin_crystal_unsup_edep_det4_cry2,
    griffin_crystal_unsup_edep_det5_cry2,
    griffin_crystal_unsup_edep_det6_cry2,
    griffin_crystal_unsup_edep_det7_cry2,
    griffin_crystal_unsup_edep_det8_cry2,
    griffin_crystal_unsup_edep_det9_cry2,
    griffin_crystal_unsup_edep_det10_cry2,
    griffin_crystal_unsup_edep_det11_cry2,
    griffin_crystal_unsup_edep_det12_cry2,
    griffin_crystal_unsup_edep_det13_cry2,
    griffin_crystal_unsup_edep_det14_cry2,
    griffin_crystal_unsup_edep_det15_cry2,
    griffin_crystal_unsup_edep_det0_cry3,
    griffin_crystal_unsup_edep_det1_cry3,
    griffin_crystal_unsup_edep_det2_cry3,
    griffin_crystal_unsup_edep_det3_cry3,
    griffin_crystal_unsup_edep_det4_cry3,
    griffin_crystal_unsup_edep_det5_cry3,
    griffin_crystal_unsup_edep_det6_cry3,
    griffin_crystal_unsup_edep_det7_cry3,
    griffin_crystal_unsup_edep_det8_cry3,
    griffin_crystal_unsup_edep_det9_cry3,
    griffin_crystal_unsup_edep_det10_cry3,
    griffin_crystal_unsup_edep_det11_cry3,
    griffin_crystal_unsup_edep_det12_cry3,
    griffin_crystal_unsup_edep_det13_cry3,
    griffin_crystal_unsup_edep_det14_cry3,
    griffin_crystal_unsup_edep_det15_cry3,
    labr_crystal_edep,
    labr_crystal_edep_sum,
    labr_crystal_edep_det0,
    labr_crystal_edep_det1,
    labr_crystal_edep_det2,
    labr_crystal_edep_det3,
    labr_crystal_edep_det4,
    labr_crystal_edep_det5,
    labr_crystal_edep_det6,
    labr_crystal_edep_det7,
    labr_crystal_edep_det8,
    labr_crystal_edep_det9,
    labr_crystal_edep_det10,
    labr_crystal_edep_det11,
    labr_crystal_edep_det12,
    labr_crystal_edep_det13,
    labr_crystal_edep_det14,
    labr_crystal_edep_det15,
    labr_crystal_edep_det16,
    labr_crystal_edep_det17,
    labr_crystal_edep_det18,
    labr_crystal_edep_det19,
    ancillary_bgo_crystal_edep,
    ancillary_bgo_crystal_edep_sum,
    ancillary_bgo_crystal_edep_det0,
    ancillary_bgo_crystal_edep_det1,
    ancillary_bgo_crystal_edep_det2,
    ancillary_bgo_crystal_edep_det3,
    ancillary_bgo_crystal_edep_det4,
    ancillary_bgo_crystal_edep_det5,
    ancillary_bgo_crystal_edep_det6,
    ancillary_bgo_crystal_edep_det7,
    ancillary_bgo_crystal_edep_det8,
    ancillary_bgo_crystal_edep_det9,
    ancillary_bgo_crystal_edep_det10,
    ancillary_bgo_crystal_edep_det11,
    ancillary_bgo_crystal_edep_det12,
    ancillary_bgo_crystal_edep_det13,
    ancillary_bgo_crystal_edep_det14,
    ancillary_bgo_crystal_edep_det15,
    ancillary_bgo_crystal_edep_det16,
    ancillary_bgo_crystal_edep_det17,
    ancillary_bgo_crystal_edep_det18,
    ancillary_bgo_crystal_edep_det19,
    sodiumIodide_crystal_edep,
    sodiumIodide_crystal_edep_sum,
    sodiumIodide_crystal_edep_det0,
    sodiumIodide_crystal_edep_det1,
    sodiumIodide_crystal_edep_det2,
    sodiumIodide_crystal_edep_det3,
    sodiumIodide_crystal_edep_det4,
    sodiumIodide_crystal_edep_det5,
    sodiumIodide_crystal_edep_det6,
    sodiumIodide_crystal_edep_det7,
    sodiumIodide_crystal_edep_det8,
    sodiumIodide_crystal_edep_det9,
    sodiumIodide_crystal_edep_det10,
    sodiumIodide_crystal_edep_det11,
    sodiumIodide_crystal_edep_det12,
    sodiumIodide_crystal_edep_det13,
    sodiumIodide_crystal_edep_det14,
    sodiumIodide_crystal_edep_det15,
    sodiumIodide_crystal_edep_det16,
    sodiumIodide_crystal_edep_det17,
    sodiumIodide_crystal_edep_det18,
    sodiumIodide_crystal_edep_det19,
    sceptar_edep,
    sceptar_edep_sum,
    sceptar_edep_det0,
    sceptar_edep_det1,
    sceptar_edep_det2,
    sceptar_edep_det3,
    sceptar_edep_det4,
    sceptar_edep_det5,
    sceptar_edep_det6,
    sceptar_edep_det7,
    sceptar_edep_det8,
    sceptar_edep_det9,
    sceptar_edep_det10,
    sceptar_edep_det11,
    sceptar_edep_det12,
    sceptar_edep_det13,
    sceptar_edep_det14,
    sceptar_edep_det15,
    sceptar_edep_det16,
    sceptar_edep_det17,
    sceptar_edep_det18,
    sceptar_edep_det19,
    Eightpi_crystal_edep,
    Eightpi_crystal_edep_sum,
    Eightpi_crystal_edep_det0,
    Eightpi_crystal_edep_det1,
    Eightpi_crystal_edep_det2,
    Eightpi_crystal_edep_det3,
    Eightpi_crystal_edep_det4,
    Eightpi_crystal_edep_det5,
    Eightpi_crystal_edep_det6,
    Eightpi_crystal_edep_det7,
    Eightpi_crystal_edep_det8,
    Eightpi_crystal_edep_det9,
    Eightpi_crystal_edep_det10,
    Eightpi_crystal_edep_det11,
    Eightpi_crystal_edep_det12,
    Eightpi_crystal_edep_det13,
    Eightpi_crystal_edep_det14,
    Eightpi_crystal_edep_det15,
    Eightpi_crystal_edep_det16,
    Eightpi_crystal_edep_det17,
    Eightpi_crystal_edep_det18,
    Eightpi_crystal_edep_det19,
    spice_crystal_edep,
    spice_crystal_edep_sum,
    spice_crystal_edep_det0,
    spice_crystal_edep_det1,
    spice_crystal_edep_det2,
    spice_crystal_edep_det3,
    spice_crystal_edep_det4,
    spice_crystal_edep_det5,
    spice_crystal_edep_det6,
    spice_crystal_edep_det7,
    spice_crystal_edep_det8,
    spice_crystal_edep_det9,
    spice_crystal_edep_det10,
    spice_crystal_edep_det11,
    spice_crystal_edep_det12,
    spice_crystal_edep_det13,
    spice_crystal_edep_det14,
    spice_crystal_edep_det15,
    spice_crystal_edep_det16,
    spice_crystal_edep_det17,
    spice_crystal_edep_det18,
    spice_crystal_edep_det19,
    paces_crystal_edep,		// This is probably the wrong number of detectors for paces.
    paces_crystal_edep_sum,
    paces_crystal_edep_det0,
    paces_crystal_edep_det1,
    paces_crystal_edep_det2,
    paces_crystal_edep_det3,
    paces_crystal_edep_det4,
    paces_crystal_edep_det5,
    paces_crystal_edep_det6,
    paces_crystal_edep_det7,
    paces_crystal_edep_det8,
    paces_crystal_edep_det9,
    paces_crystal_edep_det10,
    paces_crystal_edep_det11,
    paces_crystal_edep_det12,
    paces_crystal_edep_det13,
    paces_crystal_edep_det14,
    paces_crystal_edep_det15,
    paces_crystal_edep_det16,
    paces_crystal_edep_det17,
    paces_crystal_edep_det18,
    paces_crystal_edep_det19

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

