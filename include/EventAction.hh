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

static const int MAXSTEPS = 1000;
static const int MAXHITS = 100;
static const int NUMSTEPVARS = 14;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction*, HistoManager*);
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event*);
    virtual void    EndOfEventAction(const G4Event*);
    
    void AddAbs(G4double de, G4double dl) {fEnergyAbs += de; fTrackLAbs += dl;};
    void AddGap(G4double de, G4double dl) {fEnergyGap += de; fTrackLGap += dl;};

    G4int GetEventNumber(){return evtNb;};

    // particle types
    void AddParticleType(G4int index) {particleTypes[index] += 1;};

    // Energy deposit in detection systems
    void AddGriffinCrystDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinCrystEnergyDet[det][cry] += de; GriffinCrystTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorBackDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorBackEnergyDet[det][cry] += de; GriffinSuppressorBackTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorLeftExtensionDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorLeftExtensionEnergyDet[det][cry] += de; GriffinSuppressorLeftExtensionTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorLeftSideDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorLeftSideEnergyDet[det][cry] += de; GriffinSuppressorLeftSideTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorRightExtensionDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorRightExtensionEnergyDet[det][cry] += de; GriffinSuppressorRightExtensionTrackDet[det][cry] += dl;};
    void AddGriffinSuppressorRightSideDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorRightSideEnergyDet[det][cry] += de; GriffinSuppressorRightSideTrackDet[det][cry] += dl;};

    void Add8piCrystDet(G4double de, G4double dl, G4int det) {EightPiCrystEnergyDet[det] += de; EightPiCrystTrackDet[det] += dl;} ;

    void AddLaBrCrystDet(G4double de, G4double dl, G4int det) {LaBrCrystEnergyDet[det] += de; LaBrCrystTrackDet[det] += dl;} ;

    void AncillaryBgoDet(G4double de, G4double dl, G4int det) {AncillaryBgoEnergyDet[det] += de; AncillaryBgoTrackDet[det] += dl;} ;

    void SceptarDet(G4double de, G4double dl, G4int det) {SceptarEnergyDet[det] += de; SceptarTrackDet[det] += dl;} ;

    void AddGridCellElectron(G4double de, G4double dl, G4int det) {GridCellElectronEKinDet[det] += de; GridCellElectronTrackDet[det] += dl;} ;
    void AddGridCellGamma(G4double de, G4double dl, G4int det) {GridCellGammaEKinDet[det] += de; GridCellGammaTrackDet[det] += dl;} ;
    void AddGridCellNeutron(G4double de, G4double dl, G4int det) {GridCellNeutronEKinDet[det] += de; GridCellNeutronTrackDet[det] += dl;} ;

    void AddHitTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time);


private:
    RunAction*    fRunAct;
    HistoManager* fHistoManager;


    G4double  fEnergyAbs, fEnergyGap;
    G4double  fTrackLAbs, fTrackLGap;

    G4int     fPrintModulo;
    G4int     evtNb;

    void ClearVariables();
    void FillParticleType();
    void FillGriffinCryst();
    void Fill8piCryst() ;
    void FillLaBrCryst() ;
    void FillAncillaryBgo() ;
    void FillSceptar() ;
    void FillGridCell() ;

    // tracking info
    G4double stepTracker[NUMSTEPVARS][MAXSTEPS];

    G4String stepVolume[MAXSTEPS]; // volume at each step
    G4int    stepIndex;

    G4int    hitTrackerI[NUMSTEPVARS][MAXHITS];
    G4double hitTrackerD[NUMSTEPVARS][MAXHITS];
    G4int    hitIndex;
    G4int    pTrackID;
    G4int    pParentID;

    G4int    numberOfHits;
    G4String pHitMnemonic[MAXHITS];


    // Particle types in simulation
    G4int particleTypes[NUMPARTICLETYPES];

    // Energy deposit in detection systems
    G4double GriffinCrystEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double GriffinCrystTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double GriffinSuppressorBackEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double GriffinSuppressorBackTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];

    G4double GriffinSuppressorLeftExtensionEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double GriffinSuppressorLeftExtensionTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double GriffinSuppressorLeftSideEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double GriffinSuppressorLeftSideTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];

    G4double GriffinSuppressorRightExtensionEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double GriffinSuppressorRightExtensionTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double GriffinSuppressorRightSideEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
    G4double GriffinSuppressorRightSideTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];

    G4double EightPiCrystEnergyDet[MAXNUMDET] ;
    G4double EightPiCrystTrackDet[MAXNUMDET] ;

    G4double LaBrCrystEnergyDet[MAXNUMDET] ;
    G4double LaBrCrystTrackDet[MAXNUMDET] ;

    G4double AncillaryBgoEnergyDet[MAXNUMDET] ;
    G4double AncillaryBgoTrackDet[MAXNUMDET] ;

    G4double SceptarEnergyDet[MAXNUMDET] ;
    G4double SceptarTrackDet[MAXNUMDET] ;

    G4double GridCellElectronEKinDet[MAXNUMDET] ;
    G4double GridCellElectronTrackDet[MAXNUMDET] ;
    G4double GridCellGammaEKinDet[MAXNUMDET] ;
    G4double GridCellGammaTrackDet[MAXNUMDET] ;
    G4double GridCellNeutronEKinDet[MAXNUMDET] ;
    G4double GridCellNeutronTrackDet[MAXNUMDET] ;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


