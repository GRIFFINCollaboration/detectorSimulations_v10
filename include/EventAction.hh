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

#ifndef EVENTACTION_HH
#define EVENTACTION_HH

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

	void AddHitTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int trackerZ);
	void AddStepTracker(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int trackerZ);

    G4bool SpiceTest();
private:
	RunAction*    fRunAct;
	HistoManager* fHistoManager;

	G4int     fPrintModulo;
	G4int     fEvtNb;
	G4double BeamInputEnergy;

	void ClearVariables();

	G4int SPICEPhiRemainder(G4int);
	G4int SPICEPhiMap(G4int);

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
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


