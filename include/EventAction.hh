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
/// \file analysis/shared/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
// $Id: EventAction.hh 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "HistoManager.hh"

class RunAction;
class HistoManager;

static const int MAXSTEPS       = 1000;
static const int MAXHITS        = 100;
static const int NUMSTEPVARS    = 15;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction*);
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event*);
    virtual void    EndOfEventAction(const G4Event*);

    G4int GetEventNumber() { return fEvtNb;};

    // particle types
    void AddParticleType(G4int index) { fParticleTypes[index] += 1;};

    // Energy deposit in detection systems
    void AddGriffinCrystDet(G4double de, G4double dl, G4int det, G4int cry) { fGriffinCrystEnergyDet[det][cry] += de; fGriffinCrystTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorBackDet(G4double de, G4double dl, G4int det, G4int cry) { fGriffinSuppressorBackEnergyDet[det][cry] += de; fGriffinSuppressorBackTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorLeftExtensionDet(G4double de, G4double dl, G4int det, G4int cry) { fGriffinSuppressorLeftExtensionEnergyDet[det][cry] += de; fGriffinSuppressorLeftExtensionTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorLeftSideDet(G4double de, G4double dl, G4int det, G4int cry) { fGriffinSuppressorLeftSideEnergyDet[det][cry] += de; fGriffinSuppressorLeftSideTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorRightExtensionDet(G4double de, G4double dl, G4int det, G4int cry) { fGriffinSuppressorRightExtensionEnergyDet[det][cry] += de; fGriffinSuppressorRightExtensionTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorRightSideDet(G4double de, G4double dl, G4int det, G4int cry) { fGriffinSuppressorRightSideEnergyDet[det][cry] += de; fGriffinSuppressorRightSideTrackDet[det][cry] += dl;};

    void Add8piCrystDet(G4double de, G4double dl, G4int det) { fEightPiCrystEnergyDet[det] += de; fEightPiCrystTrackDet[det] += dl;} ;

    void AddLaBrCrystDet(G4double de, G4double dl, G4int det) { fLaBrCrystEnergyDet[det] += de; fLaBrCrystTrackDet[det] += dl;} ;

    void AncillaryBgoDet(G4double de, G4double dl, G4int det) { fAncillaryBgoEnergyDet[det] += de; fAncillaryBgoTrackDet[det] += dl;} ;

    void SceptarDet(G4double de, G4double dl, G4int det) { fSceptarEnergyDet[det] += de; fSceptarTrackDet[det] += dl;} ;

    void AddGridCellElectron(G4double de, G4double dl, G4int det) { fGridCellElectronEKinDet[det] += de; fGridCellElectronTrackDet[det] += dl;} ;
    void AddGridCellGamma(G4double de, G4double dl, G4int det) { fGridCellGammaEKinDet[det] += de; fGridCellGammaTrackDet[det] += dl;} ;
    void AddGridCellNeutron(G4double de, G4double dl, G4int det) { fGridCellNeutronEKinDet[det] += de; fGridCellNeutronTrackDet[det] += dl;} ;

    void AddHitTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int trackerZ);
    void AddStepTracker(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int trackerZ);


private:
    RunAction*    fRunAct;

    G4int     fPrintModulo;
    G4int     fEvtNb;

    void ClearVariables();
    void FillParticleType();
    void FillGriffinCryst();
    void Fill8piCryst() ;
    void FillLaBrCryst() ;
    void FillAncillaryBgo() ;
    void FillSceptar() ;
    void FillGridCell() ;

    // Tracking info
    G4int    fHitTrackerI[NUMSTEPVARS][MAXHITS];
    G4double fHitTrackerD[NUMSTEPVARS][MAXHITS];
    G4int    fHitIndex;
    G4int    fNumberOfHits;
    G4String fPHitMnemonic[MAXHITS];
				 
    G4int    fStepTrackerI[NUMSTEPVARS][MAXSTEPS];
    G4double fStepTrackerD[NUMSTEPVARS][MAXSTEPS];
    G4int    fStepIndex;
    G4int    fNumberOfSteps;
				 
    G4int    fPTrackID;
    G4int    fPParentID;


    // Particle types in simulation
    G4int fParticleTypes[NUMPARTICLETYPES];

    // Energy deposit in detection systems
    G4double fGriffinCrystEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double fGriffinCrystTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double fGriffinSuppressorBackEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double fGriffinSuppressorBackTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
				 
    G4double fGriffinSuppressorLeftExtensionEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double fGriffinSuppressorLeftExtensionTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double fGriffinSuppressorLeftSideEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double fGriffinSuppressorLeftSideTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
				 
    G4double fGriffinSuppressorRightExtensionEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double fGriffinSuppressorRightExtensionTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double fGriffinSuppressorRightSideEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double fGriffinSuppressorRightSideTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
				 
    G4double fEightPiCrystEnergyDet[MAXNUMDET] ;
    G4double fEightPiCrystTrackDet[MAXNUMDET] ;
				 
    G4double fLaBrCrystEnergyDet[MAXNUMDET] ;
    G4double fLaBrCrystTrackDet[MAXNUMDET] ;
				 
    G4double fAncillaryBgoEnergyDet[MAXNUMDET] ;
    G4double fAncillaryBgoTrackDet[MAXNUMDET] ;
				 
    G4double fSceptarEnergyDet[MAXNUMDET] ;
    G4double fSceptarTrackDet[MAXNUMDET] ;
				 
    G4double fGridCellElectronEKinDet[MAXNUMDET] ;
    G4double fGridCellElectronTrackDet[MAXNUMDET] ;
    G4double fGridCellGammaEKinDet[MAXNUMDET] ;
    G4double fGridCellGammaTrackDet[MAXNUMDET] ;
    G4double fGridCellNeutronEKinDet[MAXNUMDET] ;
    G4double fGridCellNeutronTrackDet[MAXNUMDET] ;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


