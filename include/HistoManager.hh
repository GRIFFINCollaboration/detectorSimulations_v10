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

#ifndef HISTOMANAGER_HH
#define HISTOMANAGER_HH

#include "globals.hh"
#include "g4root.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

const G4bool WRITEEKINHISTOS    = true;//bools needed to write histos
const G4bool WRITEEDEPHISTOS    = true;
const G4bool WRITETRACKLHISTOS  = true;

const G4int MAXHISTO            = 5000;//max number of histos in root file
const G4int MAXNTCOL            = 15;
const G4int MAXNUMDET           = 20;
const G4int MAXNUMDETSPICE	= 10;
const G4int MAXNUMSEGSPICE	= 12;
const G4int MAXNUMDETPACES	= 5;
const G4int MAXNUMDETGRIFFIN    = 16;
const G4int MAXNUMCRYGRIFFIN    = 4;
const G4int NUMPARTICLETYPES    = 20;

// Anita & Rishita
const G4int MAXNUMANG_GG = 11;
const G4int MAXNUMANG_GE = 18;
const G4int MAXNUMANG = std::max(MAXNUMANG_GG,MAXNUMANG_GE);
const G4int MAXANGCORRHISTO = 3;

// ekin histo properties    ///////////////////////
const G4int     EKINNBINS  = 10000;
const G4double  EKINXMIN   = 0.5*keV;
const G4double  EKINXMAX   = 10000.5*keV;
const G4int     EKINNBINSSPICE  = 10000;//spice histos range different
const G4double  EKINXMINSPICE   = 0.5*keV;
const G4double  EKINXMAXSPICE   = 10000.5*keV;

// edep histo properties    ///////////////////////
const G4int     EDEPNBINS  = 10000;//was 10000
const G4double  EDEPXMIN   = 0.5*keV;
const G4double  EDEPXMAX   = 10000.5*keV;//was 10000.5
const G4int     EDEPNBINSSPICE  = 10000;//was 10000	//spice histos range different
const G4double  EDEPXMINSPICE   = 0.*keV;
const G4double  EDEPXMAXSPICE   = 2100.0*keV;//was 10000.5
const G4int 	EDEPNBINSGRIFFIN = 1500;
const G4double	EDEPXMINGRIFFIN = 0.0;
const G4double 	EDEPXMAXGRIFFIN	= 1500.0;

// trackl histo properties  ///////////////////////
const G4int     TRACKLNBINS = 5000;
const G4double  TRACKLXMIN  = 0.5*mm;
const G4double  TRACKLXMAX  = 5000.5*mm;
const G4int     TRACKLNBINSSPICE = 5000;//spice histos range different
const G4double  TRACKLXMINSPICE  = 0.5*mm;
const G4double  TRACKLXMAXSPICE  = 5000.5*mm;

///////////////////////////////////////////////////
const G4double MINENERGYTHRES   = 0.001*keV;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager {
public:
    static HistoManager& Instance() { // static so only ever one instance
		static HistoManager histManager;
		return histManager;
	 }

    void Book();
    void Save();

// Anita & Rishita ----------------------------
// replace MAXNUMANG_GE with the highest valued MAXNUMANG variable
short AngCorrNumbers[MAXNUMANG_GE*MAXANGCORRHISTO+MAXANGCORRHISTO];
short fAngCorrAngles[1];
G4bool   gg = false;
G4bool   ge = true;
G4bool   ee = false;
//---------------------------------------------
    
    //arrays allow histos to be made independently of GRIFFIN
    //	sizeof() arrays appear double the elements, due to 2 bits per element
    short fPacesHistNumbers[MAXNUMDETPACES+2]; //+2 for edep and sum histos 
    short fSpiceHistNumbers[MAXNUMDETSPICE*MAXNUMSEGSPICE+2]; //+2 for edep and sum histos 
    short fSegmentHisto[120];//this array will hold segment IDs to be transferred to reference for histos
    short fAngleDistro[10]; //this variable will hold the histogran ID for the angular distributions from the cone
    
    // BEGIN GRIFFIN arrays
    short GriffinSup[MAXNUMDETGRIFFIN+2];
	// edep is stored at index 0
	// edep_sum is stored at index 1
    short GriffinSupCry[MAXNUMDETGRIFFIN*MAXNUMCRYGRIFFIN+1];
	// edep_cry is stored at index 0
    short GriffinUnSup[MAXNUMDETGRIFFIN+2];
    short GriffinUnSupCry[MAXNUMDETGRIFFIN*MAXNUMCRYGRIFFIN+1];
    // END GRIFFIN arrays

    void MakeHisto(G4AnalysisManager* analysisManager, G4String filename,  G4String title, G4double xmin, G4double xmax, G4int nbins);
    void MakeHistoWithAxisTitles(G4AnalysisManager* analysisManager, G4String name, 
					   G4String title, G4double xmin, G4double xmax, 
					   G4int nbins, const G4String& unitName, const G4String& fcnName);
    void Make2DHistoWithAxisTitles(G4AnalysisManager* analysisManager, const G4String& name, const G4String& title,
                 G4int nxbins, G4double xmin, G4double xmax, 
                 G4int nybins, G4double ymin, G4double ymax);
    void FillHisto(G4int ih, G4double e, G4double weight = 1.0);
    void Fill2DHisto(G4int ih, G4double xbin, G4double ybin, G4double weight = 1.0);
    void Normalize(G4int id, G4double fac);

    void FillHitNtuple(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ);
    void FillStepNtuple(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ);

    void PrintStatistic();

    G4bool GetStepTrackerBool() { return fStepTrackerBool;};
    G4bool GetHitTrackerBool()  { return fHitTrackerBool;};

	 // setter
	 void GridCell  (G4bool val)   { fGridCell = val; }
 	 void Griffin   (G4bool val)   { fGriffin  = val; }
	 void LaBr 	    (G4bool val)   { fLaBr     = val; }
	 void AncBgo    (G4bool val)   { fAncBgo   = val; }
	 void NaI 	    (G4bool val)   { fNaI      = val; }
	 void Sceptar   (G4bool val)   { fSceptar  = val; }
	 void EightPi   (G4bool val)   { fEightPi  = val; }
	 void Spice     (G4bool val)   { fSpice    = val; }
	 void Paces     (G4bool val)   { fPaces    = val; }  
	 void Descant   (G4bool val)   { fDescant  = val; }  
	 void Testcan   (G4bool val)   { fTestcan  = val; }  
	 void BeamEnergy(G4double val) { fBeamEnergy = val; }

	 // getter
	 G4bool   GridCell()   { return fGridCell;   }
	 G4bool   Griffin()    { return fGriffin;    }
	 G4bool   LaBr()       { return fLaBr;       }
	 G4bool   AncBgo()     { return fAncBgo;     }
	 G4bool   NaI()        { return fNaI;        }
	 G4bool   Sceptar()    { return fSceptar;    }
	 G4bool   EightPi()    { return fEightPi;    }
	 G4bool   Spice()      { return fSpice;      }
	 G4bool   Paces()      { return fPaces;      }  
	 G4bool   Descant()    { return fDescant;    }  
	 G4bool   Testcan()    { return fTestcan;    } 
	 G4double BeamEnergy() { return fBeamEnergy; }
	 
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
    G4AnaH2*      fHistPt2[MAXHISTO];
    
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
    G4bool fDescant;
    G4bool fTestcan;	 
    G4bool fSpice;
    G4bool fPaces;
    G4bool fAngCorr; //Anita
	 G4double fBeamEnergy;
};

enum HISTONAME
{
    kNullName = 0,
    kAstatsParticleTypeInEachStep,
    kAstatsParticleTypeInEachEvent,
   /* kSpiceEdep,
    kSpiceEdepSum,
    kSpiceEdepDet0,
    kSpiceEdepDet01,
    kSpiceEdepDet02,
    kSpiceEdepDet03,
    kSpiceEdepDet04,
    kSpiceEdepDet05,
    kSpiceEdepDet06,
    kSpiceEdepDet07,
    kSpiceEdepDet08,
    kSpiceEdepDet09, //max number to be reached
    kSpiceEdepDet10,
    kSpiceEdepDet11,
    kSpiceEdepDet12,
    kSpiceEdepDet13,
    kSpiceEdepDet14,//these are no longer used
    kSpiceEdepDet15,
    kSpiceEdepDet16,
    kSpiceEdepDet17,
    kSpiceEdepDet18,
    kSpiceEdepDet19,*/
    kPacesCrystalEdep,		// This is probably the wrong number of detectors for paces.
    kPacesCrystalEdepSum,
    kPacesCrystalEdepDet0,
    kPacesCrystalEdepDet1,
    kPacesCrystalEdepDet2,
    kPacesCrystalEdepDet3,
    kPacesCrystalEdepDet4, // reduced to 5 total (zero-indexed) as number of detectors 22/6
  /* kPacesCrystalEdepDet5,
    kPacesCrystalEdepDet6,
    kPacesCrystalEdepDet7,
    kPacesCrystalEdepDet8,
    kPacesCrystalEdepDet9,
    kPacesCrystalEdepDet10,
    kPacesCrystalEdepDet11,
    kPacesCrystalEdepDet12,
    kPacesCrystalEdepDet13,
    kPacesCrystalEdepDet14,
    kPacesCrystalEdepDet15,
    kPacesCrystalEdepDet16,
    kPacesCrystalEdepDet17,
    kPacesCrystalEdepDet18,
    kPacesCrystalEdepDet19, */
    kGridcellElectronEkinDet0,
    kGridcellElectronEkinDet1,
    kGridcellElectronEkinDet2,
    kGridcellElectronEkinDet3,
    kGridcellElectronEkinDet4,
    kGridcellElectronEkinDet5,
    kGridcellElectronEkinDet6,
    kGridcellElectronEkinDet7,
    kGridcellElectronEkinDet8,
    kGridcellElectronEkinDet9,
    kGridcellElectronEkinDet10,
    kGridcellElectronEkinDet11,
    kGridcellElectronEkinDet12,
    kGridcellElectronEkinDet13,
    kGridcellElectronEkinDet14,
    kGridcellElectronEkinDet15,
    kGridcellElectronEkinDet16,
    kGridcellElectronEkinDet17,
    kGridcellElectronEkinDet18,
    kGridcellElectronEkinDet19,
    kGridcellElectronTracklDet0,
    kGridcellElectronTracklDet1,
    kGridcellElectronTracklDet2,
    kGridcellElectronTracklDet3,
    kGridcellElectronTracklDet4,
    kGridcellElectronTracklDet5,
    kGridcellElectronTracklDet6,
    kGridcellElectronTracklDet7,
    kGridcellElectronTracklDet8,
    kGridcellElectronTracklDet9,
    kGridcellElectronTracklDet10,
    kGridcellElectronTracklDet11,
    kGridcellElectronTracklDet12,
    kGridcellElectronTracklDet13,
    kGridcellElectronTracklDet14,
    kGridcellElectronTracklDet15,
    kGridcellElectronTracklDet16,
    kGridcellElectronTracklDet17,
    kGridcellElectronTracklDet18,
    kGridcellElectronTracklDet19,
    kGridcellGammaEkinDet0,
    kGridcellGammaEkinDet1,
    kGridcellGammaEkinDet2,
    kGridcellGammaEkinDet3,
    kGridcellGammaEkinDet4,
    kGridcellGammaEkinDet5,
    kGridcellGammaEkinDet6,
    kGridcellGammaEkinDet7,
    kGridcellGammaEkinDet8,
    kGridcellGammaEkinDet9,
    kGridcellGammaEkinDet10,
    kGridcellGammaEkinDet11,
    kGridcellGammaEkinDet12,
    kGridcellGammaEkinDet13,
    kGridcellGammaEkinDet14,
    kGridcellGammaEkinDet15,
    kGridcellGammaEkinDet16,
    kGridcellGammaEkinDet17,
    kGridcellGammaEkinDet18,
    kGridcellGammaEkinDet19,
    kGridcellGammaTracklDet0,
    kGridcellGammaTracklDet1,
    kGridcellGammaTracklDet2,
    kGridcellGammaTracklDet3,
    kGridcellGammaTracklDet4,
    kGridcellGammaTracklDet5,
    kGridcellGammaTracklDet6,
    kGridcellGammaTracklDet7,
    kGridcellGammaTracklDet8,
    kGridcellGammaTracklDet9,
    kGridcellGammaTracklDet10,
    kGridcellGammaTracklDet11,
    kGridcellGammaTracklDet12,
    kGridcellGammaTracklDet13,
    kGridcellGammaTracklDet14,
    kGridcellGammaTracklDet15,
    kGridcellGammaTracklDet16,
    kGridcellGammaTracklDet17,
    kGridcellGammaTracklDet18,
    kGridcellGammaTracklDet19,
    kGridcellNeutronEkinDet0,
    kGridcellNeutronEkinDet1,
    kGridcellNeutronEkinDet2,
    kGridcellNeutronEkinDet3,
    kGridcellNeutronEkinDet4,
    kGridcellNeutronEkinDet5,
    kGridcellNeutronEkinDet6,
    kGridcellNeutronEkinDet7,
    kGridcellNeutronEkinDet8,
    kGridcellNeutronEkinDet9,
    kGridcellNeutronEkinDet10,
    kGridcellNeutronEkinDet11,
    kGridcellNeutronEkinDet12,
    kGridcellNeutronEkinDet13,
    kGridcellNeutronEkinDet14,
    kGridcellNeutronEkinDet15,
    kGridcellNeutronEkinDet16,
    kGridcellNeutronEkinDet17,
    kGridcellNeutronEkinDet18,
    kGridcellNeutronEkinDet19,
    kGridcellNeutronTracklDet0,
    kGridcellNeutronTracklDet1,
    kGridcellNeutronTracklDet2,
    kGridcellNeutronTracklDet3,
    kGridcellNeutronTracklDet4,
    kGridcellNeutronTracklDet5,
    kGridcellNeutronTracklDet6,
    kGridcellNeutronTracklDet7,
    kGridcellNeutronTracklDet8,
    kGridcellNeutronTracklDet9,
    kGridcellNeutronTracklDet10,
    kGridcellNeutronTracklDet11,
    kGridcellNeutronTracklDet12,
    kGridcellNeutronTracklDet13,
    kGridcellNeutronTracklDet14,
    kGridcellNeutronTracklDet15,
    kGridcellNeutronTracklDet16,
    kGridcellNeutronTracklDet17,
    kGridcellNeutronTracklDet18,
    kGridcellNeutronTracklDet19,
    kGriffinCrystalSupEdep,
    kGriffinCrystalSupEdepCry,
    kGriffinCrystalSupEdepSum,
    kGriffinCrystalSupEdepDet0,
    kGriffinCrystalSupEdepDet1,
    kGriffinCrystalSupEdepDet2,
    kGriffinCrystalSupEdepDet3,
    kGriffinCrystalSupEdepDet4,
    kGriffinCrystalSupEdepDet5,
    kGriffinCrystalSupEdepDet6,
    kGriffinCrystalSupEdepDet7,
    kGriffinCrystalSupEdepDet8,
    kGriffinCrystalSupEdepDet9,
    kGriffinCrystalSupEdepDet10,
    kGriffinCrystalSupEdepDet11,
    kGriffinCrystalSupEdepDet12,
    kGriffinCrystalSupEdepDet13,
    kGriffinCrystalSupEdepDet14,
    kGriffinCrystalSupEdepDet15,
    kGriffinCrystalSupEdepDet0Cry0,
    kGriffinCrystalSupEdepDet1Cry0,
    kGriffinCrystalSupEdepDet2Cry0,
    kGriffinCrystalSupEdepDet3Cry0,
    kGriffinCrystalSupEdepDet4Cry0,
    kGriffinCrystalSupEdepDet5Cry0,
    kGriffinCrystalSupEdepDet6Cry0,
    kGriffinCrystalSupEdepDet7Cry0,
    kGriffinCrystalSupEdepDet8Cry0,
    kGriffinCrystalSupEdepDet9Cry0,
    kGriffinCrystalSupEdepDet10Cry0,
    kGriffinCrystalSupEdepDet11Cry0,
    kGriffinCrystalSupEdepDet12Cry0,
    kGriffinCrystalSupEdepDet13Cry0,
    kGriffinCrystalSupEdepDet14Cry0,
    kGriffinCrystalSupEdepDet15Cry0,
    kGriffinCrystalSupEdepDet0Cry1,
    kGriffinCrystalSupEdepDet1Cry1,
    kGriffinCrystalSupEdepDet2Cry1,
    kGriffinCrystalSupEdepDet3Cry1,
    kGriffinCrystalSupEdepDet4Cry1,
    kGriffinCrystalSupEdepDet5Cry1,
    kGriffinCrystalSupEdepDet6Cry1,
    kGriffinCrystalSupEdepDet7Cry1,
    kGriffinCrystalSupEdepDet8Cry1,
    kGriffinCrystalSupEdepDet9Cry1,
    kGriffinCrystalSupEdepDet10Cry1,
    kGriffinCrystalSupEdepDet11Cry1,
    kGriffinCrystalSupEdepDet12Cry1,
    kGriffinCrystalSupEdepDet13Cry1,
    kGriffinCrystalSupEdepDet14Cry1,
    kGriffinCrystalSupEdepDet15Cry1,
    kGriffinCrystalSupEdepDet0Cry2,
    kGriffinCrystalSupEdepDet1Cry2,
    kGriffinCrystalSupEdepDet2Cry2,
    kGriffinCrystalSupEdepDet3Cry2,
    kGriffinCrystalSupEdepDet4Cry2,
    kGriffinCrystalSupEdepDet5Cry2,
    kGriffinCrystalSupEdepDet6Cry2,
    kGriffinCrystalSupEdepDet7Cry2,
    kGriffinCrystalSupEdepDet8Cry2,
    kGriffinCrystalSupEdepDet9Cry2,
    kGriffinCrystalSupEdepDet10Cry2,
    kGriffinCrystalSupEdepDet11Cry2,
    kGriffinCrystalSupEdepDet12Cry2,
    kGriffinCrystalSupEdepDet13Cry2,
    kGriffinCrystalSupEdepDet14Cry2,
    kGriffinCrystalSupEdepDet15Cry2,
    kGriffinCrystalSupEdepDet0Cry3,
    kGriffinCrystalSupEdepDet1Cry3,
    kGriffinCrystalSupEdepDet2Cry3,
    kGriffinCrystalSupEdepDet3Cry3,
    kGriffinCrystalSupEdepDet4Cry3,
    kGriffinCrystalSupEdepDet5Cry3,
    kGriffinCrystalSupEdepDet6Cry3,
    kGriffinCrystalSupEdepDet7Cry3,
    kGriffinCrystalSupEdepDet8Cry3,
    kGriffinCrystalSupEdepDet9Cry3,
    kGriffinCrystalSupEdepDet10Cry3,
    kGriffinCrystalSupEdepDet11Cry3,
    kGriffinCrystalSupEdepDet12Cry3,
    kGriffinCrystalSupEdepDet13Cry3,
    kGriffinCrystalSupEdepDet14Cry3,
    kGriffinCrystalSupEdepDet15Cry3,
    kGriffinCrystalUnsupEdep,
    kGriffinCrystalUnsupEdepCry,
    kGriffinCrystalUnsupEdepSum,
    kGriffinCrystalUnsupEdepDet0,
    kGriffinCrystalUnsupEdepDet1,
    kGriffinCrystalUnsupEdepDet2,
    kGriffinCrystalUnsupEdepDet3,
    kGriffinCrystalUnsupEdepDet4,
    kGriffinCrystalUnsupEdepDet5,
    kGriffinCrystalUnsupEdepDet6,
    kGriffinCrystalUnsupEdepDet7,
    kGriffinCrystalUnsupEdepDet8,
    kGriffinCrystalUnsupEdepDet9,
    kGriffinCrystalUnsupEdepDet10,
    kGriffinCrystalUnsupEdepDet11,
    kGriffinCrystalUnsupEdepDet12,
    kGriffinCrystalUnsupEdepDet13,
    kGriffinCrystalUnsupEdepDet14,
    kGriffinCrystalUnsupEdepDet15,
    kGriffinCrystalUnsupEdepDet0Cry0,
    kGriffinCrystalUnsupEdepDet1Cry0,
    kGriffinCrystalUnsupEdepDet2Cry0,
    kGriffinCrystalUnsupEdepDet3Cry0,
    kGriffinCrystalUnsupEdepDet4Cry0,
    kGriffinCrystalUnsupEdepDet5Cry0,
    kGriffinCrystalUnsupEdepDet6Cry0,
    kGriffinCrystalUnsupEdepDet7Cry0,
    kGriffinCrystalUnsupEdepDet8Cry0,
    kGriffinCrystalUnsupEdepDet9Cry0,
    kGriffinCrystalUnsupEdepDet10Cry0,
    kGriffinCrystalUnsupEdepDet11Cry0,
    kGriffinCrystalUnsupEdepDet12Cry0,
    kGriffinCrystalUnsupEdepDet13Cry0,
    kGriffinCrystalUnsupEdepDet14Cry0,
    kGriffinCrystalUnsupEdepDet15Cry0,
    kGriffinCrystalUnsupEdepDet0Cry1,
    kGriffinCrystalUnsupEdepDet1Cry1,
    kGriffinCrystalUnsupEdepDet2Cry1,
    kGriffinCrystalUnsupEdepDet3Cry1,
    kGriffinCrystalUnsupEdepDet4Cry1,
    kGriffinCrystalUnsupEdepDet5Cry1,
    kGriffinCrystalUnsupEdepDet6Cry1,
    kGriffinCrystalUnsupEdepDet7Cry1,
    kGriffinCrystalUnsupEdepDet8Cry1,
    kGriffinCrystalUnsupEdepDet9Cry1,
    kGriffinCrystalUnsupEdepDet10Cry1,
    kGriffinCrystalUnsupEdepDet11Cry1,
    kGriffinCrystalUnsupEdepDet12Cry1,
    kGriffinCrystalUnsupEdepDet13Cry1,
    kGriffinCrystalUnsupEdepDet14Cry1,
    kGriffinCrystalUnsupEdepDet15Cry1,
    kGriffinCrystalUnsupEdepDet0Cry2,
    kGriffinCrystalUnsupEdepDet1Cry2,
    kGriffinCrystalUnsupEdepDet2Cry2,
    kGriffinCrystalUnsupEdepDet3Cry2,
    kGriffinCrystalUnsupEdepDet4Cry2,
    kGriffinCrystalUnsupEdepDet5Cry2,
    kGriffinCrystalUnsupEdepDet6Cry2,
    kGriffinCrystalUnsupEdepDet7Cry2,
    kGriffinCrystalUnsupEdepDet8Cry2,
    kGriffinCrystalUnsupEdepDet9Cry2,
    kGriffinCrystalUnsupEdepDet10Cry2,
    kGriffinCrystalUnsupEdepDet11Cry2,
    kGriffinCrystalUnsupEdepDet12Cry2,
    kGriffinCrystalUnsupEdepDet13Cry2,
    kGriffinCrystalUnsupEdepDet14Cry2,
    kGriffinCrystalUnsupEdepDet15Cry2,
    kGriffinCrystalUnsupEdepDet0Cry3,
    kGriffinCrystalUnsupEdepDet1Cry3,
    kGriffinCrystalUnsupEdepDet2Cry3,
    kGriffinCrystalUnsupEdepDet3Cry3,
    kGriffinCrystalUnsupEdepDet4Cry3,
    kGriffinCrystalUnsupEdepDet5Cry3,
    kGriffinCrystalUnsupEdepDet6Cry3,
    kGriffinCrystalUnsupEdepDet7Cry3,
    kGriffinCrystalUnsupEdepDet8Cry3,
    kGriffinCrystalUnsupEdepDet9Cry3,
    kGriffinCrystalUnsupEdepDet10Cry3,
    kGriffinCrystalUnsupEdepDet11Cry3,
    kGriffinCrystalUnsupEdepDet12Cry3,
    kGriffinCrystalUnsupEdepDet13Cry3,
    kGriffinCrystalUnsupEdepDet14Cry3,
    kGriffinCrystalUnsupEdepDet15Cry3,
    kLabrCrystalEdep,
    kLabrCrystalEdepSum,
    kLabrCrystalEdepDet0,
    kLabrCrystalEdepDet1,
    kLabrCrystalEdepDet2,
    kLabrCrystalEdepDet3,
    kLabrCrystalEdepDet4,
    kLabrCrystalEdepDet5,
    kLabrCrystalEdepDet6,
    kLabrCrystalEdepDet7,
    kLabrCrystalEdepDet8,
    kLabrCrystalEdepDet9,
    kLabrCrystalEdepDet10,
    kLabrCrystalEdepDet11,
    kLabrCrystalEdepDet12,
    kLabrCrystalEdepDet13,
    kLabrCrystalEdepDet14,
    kLabrCrystalEdepDet15,
    kLabrCrystalEdepDet16,
    kLabrCrystalEdepDet17,
    kLabrCrystalEdepDet18,
    kLabrCrystalEdepDet19,
    kAncillaryBgoCrystalEdep,
    kAncillaryBgoCrystalEdepSum,
    kAncillaryBgoCrystalEdepDet0,
    kAncillaryBgoCrystalEdepDet1,
    kAncillaryBgoCrystalEdepDet2,
    kAncillaryBgoCrystalEdepDet3,
    kAncillaryBgoCrystalEdepDet4,
    kAncillaryBgoCrystalEdepDet5,
    kAncillaryBgoCrystalEdepDet6,
    kAncillaryBgoCrystalEdepDet7,
    kAncillaryBgoCrystalEdepDet8,
    kAncillaryBgoCrystalEdepDet9,
    kAncillaryBgoCrystalEdepDet10,
    kAncillaryBgoCrystalEdepDet11,
    kAncillaryBgoCrystalEdepDet12,
    kAncillaryBgoCrystalEdepDet13,
    kAncillaryBgoCrystalEdepDet14,
    kAncillaryBgoCrystalEdepDet15,
    kAncillaryBgoCrystalEdepDet16,
    kAncillaryBgoCrystalEdepDet17,
    kAncillaryBgoCrystalEdepDet18,
    kAncillaryBgoCrystalEdepDet19,
    kSodiumIodideCrystalEdep,
    kSodiumIodideCrystalEdepSum,
    kSodiumIodideCrystalEdepDet0,
    kSodiumIodideCrystalEdepDet1,
    kSodiumIodideCrystalEdepDet2,
    kSodiumIodideCrystalEdepDet3,
    kSodiumIodideCrystalEdepDet4,
    kSodiumIodideCrystalEdepDet5,
    kSodiumIodideCrystalEdepDet6,
    kSodiumIodideCrystalEdepDet7,
    kSodiumIodideCrystalEdepDet8,
    kSodiumIodideCrystalEdepDet9,
    kSodiumIodideCrystalEdepDet10,
    kSodiumIodideCrystalEdepDet11,
    kSodiumIodideCrystalEdepDet12,
    kSodiumIodideCrystalEdepDet13,
    kSodiumIodideCrystalEdepDet14,
    kSodiumIodideCrystalEdepDet15,
    kSodiumIodideCrystalEdepDet16,
    kSodiumIodideCrystalEdepDet17,
    kSodiumIodideCrystalEdepDet18,
    kSodiumIodideCrystalEdepDet19,
    kSceptarEdep,
    kSceptarEdepSum,
    kSceptarEdepDet0,
    kSceptarEdepDet1,
    kSceptarEdepDet2,
    kSceptarEdepDet3,
    kSceptarEdepDet4,
    kSceptarEdepDet5,
    kSceptarEdepDet6,
    kSceptarEdepDet7,
    kSceptarEdepDet8,
    kSceptarEdepDet9,
    kSceptarEdepDet10,
    kSceptarEdepDet11,
    kSceptarEdepDet12,
    kSceptarEdepDet13,
    kSceptarEdepDet14,
    kSceptarEdepDet15,
    kSceptarEdepDet16,
    kSceptarEdepDet17,
    kSceptarEdepDet18,
    kSceptarEdepDet19,
    kEightpiCrystalEdep,
    kEightpiCrystalEdepSum,
    kEightpiCrystalEdepDet0,
    kEightpiCrystalEdepDet1,
    kEightpiCrystalEdepDet2,
    kEightpiCrystalEdepDet3,
    kEightpiCrystalEdepDet4,
    kEightpiCrystalEdepDet5,
    kEightpiCrystalEdepDet6,
    kEightpiCrystalEdepDet7,
    kEightpiCrystalEdepDet8,
    kEightpiCrystalEdepDet9,
    kEightpiCrystalEdepDet10,
    kEightpiCrystalEdepDet11,
    kEightpiCrystalEdepDet12,
    kEightpiCrystalEdepDet13,
    kEightpiCrystalEdepDet14,
    kEightpiCrystalEdepDet15,
    kEightpiCrystalEdepDet16,
    kEightpiCrystalEdepDet17,
    kEightpiCrystalEdepDet18,
    kEightpiCrystalEdepDet19
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

