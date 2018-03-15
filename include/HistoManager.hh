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

const G4int MAXNTCOL            = 15;

///////////////////////////////////////////////////
const G4double MINENERGYTHRES   = 0.001*keV;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager {
public:
    HistoManager();
    ~HistoManager();

    void Book();
    void Save();
    
    
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
	 G4bool LaBr    () { return fLaBr;     }
	 G4bool AncBgo  () { return fAncBgo;   }
	 G4bool NaI 	() { return fNaI;      }
	 G4bool Sceptar () { return fSceptar;  }
	 G4bool EightPi () { return fEightPi;  }
	 G4bool Spice 	() { return fSpice;    }
	 G4bool Paces   () { return fPaces;    }  
	 G4bool Descant () { return fDescant;  }  
	 G4bool Testcan () { return fTestcan;  } 

private:
    G4String G4intToG4String(G4int value);

    G4bool        fFactoryOn;
    G4String      fFileName[2];

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
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

