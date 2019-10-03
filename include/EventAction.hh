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

#include "DetectorConstruction.hh" // for DetectorProperties

class RunAction;
class HistoManager;

static const int MAXSTEPS       = 1000;
static const int MAXHITS        = 100;
static const int NUMSTEPVARS    = 21;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
	EventAction(RunAction*, HistoManager*);
	virtual ~EventAction();

	virtual void  BeginOfEventAction(const G4Event*);
	virtual void    EndOfEventAction(const G4Event*);

	G4int GetEventNumber() { return fEvtNb;};

<<<<<<< HEAD
	void AddHitTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos, G4double TOFMulti, G4ThreeVector TOFPosMulti, G4double PEkin, G4double PEdep);
	void AddStepTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos, G4double TOFMulti, G4ThreeVector TOFPosMulti, G4double PEkin, G4double PEdep);
=======
<<<<<<< HEAD
	void AddHitTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos);
	void AddStepTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos);
=======
<<<<<<< HEAD
	void AddHitTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos);
	void AddStepTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos);
=======
<<<<<<< HEAD
	void AddHitTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos);
	void AddStepTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle, G4double TOF, G4ThreeVector TOFPos);
=======
<<<<<<< HEAD
	void AddHitTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle);
	void AddStepTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle, G4double final_angle);
=======
	void AddHitTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle);
	void AddStepTracker(const DetectorProperties& properties, const G4int& eventNumber, const G4int& trackID, const G4int& parentID, const G4int& stepNumber, const G4int& particleType, const G4int& processType, const G4double& depEnergy, const G4ThreeVector& pos, const G4double& time, const G4int& trackerZ, G4int total, G4int elastic, G4int inelastic, G4int numScintPhotons, G4double lab_angle);
>>>>>>> 6828e1eec50f49505f6682d4e9b4096511be47f0
>>>>>>> 47c502e8e4318a844071ec8793822e1e8af02bf9
>>>>>>> 3064246b7a7818cc5e2250e0a2139fa8a7c21204
>>>>>>> 6487c29a178d340690563b4cc1ec43aceac6d944
>>>>>>> 77f9ed17803d5c660639c899c05ad4a9b1bdcd29


	// Energy deposit in detection systems
	void SpiceDet(G4double de, G4double dl, G4int det, G4int seg) { fSpiceEnergyDet[det][seg] += de; fSpiceTrackDet[det][seg] += dl;};///19/7

	G4bool SpiceTest();

	//Lab Angle
	G4double GetLabAngle() {return pLabAngle; };
	void SetLabAngle(G4double angle) {pLabAngle = angle; };
	
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 47c502e8e4318a844071ec8793822e1e8af02bf9
>>>>>>> 3064246b7a7818cc5e2250e0a2139fa8a7c21204
>>>>>>> 6487c29a178d340690563b4cc1ec43aceac6d944
>>>>>>> 77f9ed17803d5c660639c899c05ad4a9b1bdcd29
	//Final Angle
	G4double GetFinalAngle() {return fFinalAngle; };
	void SetFinalAngle(G4double angle2) {fFinalAngle = angle2; };
	void SetFinalMomentum(G4ThreeVector finalMomentum) {fFinalMomentum = finalMomentum; };
	void SetInitialMomentum(G4ThreeVector initialMomentum) {fInitialMomentum = initialMomentum; };
	G4ThreeVector GetFinalMomentum() {return fFinalMomentum;};
	G4ThreeVector GetInitialMomentum() {return fInitialMomentum;};
	

<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
=======
>>>>>>> 6828e1eec50f49505f6682d4e9b4096511be47f0
>>>>>>> 47c502e8e4318a844071ec8793822e1e8af02bf9
>>>>>>> 3064246b7a7818cc5e2250e0a2139fa8a7c21204
>>>>>>> 6487c29a178d340690563b4cc1ec43aceac6d944
>>>>>>> 77f9ed17803d5c660639c899c05ad4a9b1bdcd29
	//Counting for Efficiencies
	void totalCounter() {++TotalCounter;};
	G4int GetTotalCounter() {return TotalCounter;};
	void SetTotalCounter(int count) {TotalCounter = count;};
	void elasticCounter() {++ElasticCounter;};
	G4int GetElasticCounter() {return ElasticCounter;};
	void SetElasticCounter(int count1) {ElasticCounter=count1;};
	void inelasticCounter() {++InelasticCounter;};
	G4int GetInelasticCounter() {return InelasticCounter;};
	void SetInelasticCounter(int count2) {InelasticCounter=count2;};

<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 3064246b7a7818cc5e2250e0a2139fa8a7c21204
>>>>>>> 6487c29a178d340690563b4cc1ec43aceac6d944
>>>>>>> 77f9ed17803d5c660639c899c05ad4a9b1bdcd29
	//Getting TOF
	void SetTOF(G4double scatterTime) {fTOF = scatterTime;};
	G4double GetTOF() {return fTOF;};
	void SetTOFPos(G4ThreeVector scatterPos) {fTOFPos = scatterPos;};
	G4ThreeVector GetTOFPos() {return fTOFPos;};
<<<<<<< HEAD
	
	void SetTOFMulti(G4double scatterTimeMulti) {fTOFMulti = scatterTimeMulti;};
	G4double GetTOFMulti() {return fTOFMulti;};
	void SetTOFPosMulti(G4ThreeVector scatterPosMulti) {fTOFPosMulti = scatterPosMulti;};
	G4ThreeVector GetTOFPosMulti() {return fTOFPosMulti;};

=======

<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
=======
>>>>>>> 47c502e8e4318a844071ec8793822e1e8af02bf9
>>>>>>> 3064246b7a7818cc5e2250e0a2139fa8a7c21204
>>>>>>> 6487c29a178d340690563b4cc1ec43aceac6d944
>>>>>>> 77f9ed17803d5c660639c899c05ad4a9b1bdcd29
	//Counting Scintillation photons
	void CountOneScintPhoton() {++totScintPhotons;};
	G4int GetTotScintPhoton() {return totScintPhotons;};
	void SetTotScintPhoton(G4int total) {totScintPhotons = total;};
<<<<<<< HEAD
	
	//Edep in Plastics
	void AddPEdep(G4double pedep) {fPEdep = fPEdep + pedep;};
	G4double GetPEdep() {return fPEdep;};
	void SetPEdep(G4double Pedep) {fPEdep = Pedep;};
	G4double GetPEkin() {return fPEkin;};
	void SetPEkin(G4double Pekin) {fPEkin = Pekin;};
=======
>>>>>>> 77f9ed17803d5c660639c899c05ad4a9b1bdcd29

private:
	RunAction*    fRunAction;
	HistoManager* fHistoManager;

	G4int     fPrintModulo;
	G4int     fEvtNb;

	void ClearVariables();
	void FillSpice();

	G4double ApplySpiceRes(G4double);
	void SetupSpiceErfc();
	G4double SpiceErfc();

	//Applying a resolution to SPCIE energies if desired
	G4double fAmp[10000];//Bin amp (effective y)
	G4double fAmpx[10000];
	G4double fAmpTot;

	// Tracking info
	G4int    fHitTrackerI[NUMSTEPVARS][MAXHITS];
	G4double fHitTrackerD[NUMSTEPVARS][MAXHITS];
	G4int    fNumberOfHits;
	DetectorProperties fProperties[MAXHITS];

	G4int    fStepTrackerI[NUMSTEPVARS][MAXSTEPS];
	G4double fStepTrackerD[NUMSTEPVARS][MAXSTEPS];
	G4int    fNumberOfSteps;

	// Energy deposit in detection systems
	G4double fSpiceEnergyDet[10][12];
	G4double fSpiceTrackDet[10][12];

	//Lab Angle
	G4double pLabAngle;
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 3064246b7a7818cc5e2250e0a2139fa8a7c21204
>>>>>>> 6487c29a178d340690563b4cc1ec43aceac6d944
>>>>>>> 77f9ed17803d5c660639c899c05ad4a9b1bdcd29
	G4double fFinalAngle;
	G4ThreeVector fFinalMomentum;
	G4ThreeVector fInitialMomentum;
	//TOF
	G4double fTOF;
	G4ThreeVector fTOFPos;
<<<<<<< HEAD
	G4double fTOFMulti;
	G4ThreeVector fTOFPosMulti;
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
=======
<<<<<<< HEAD
	G4double fFinalAngle;
	G4ThreeVector fFinalMomentum;
	G4ThreeVector fInitialMomentum;
	
=======
>>>>>>> 6828e1eec50f49505f6682d4e9b4096511be47f0
>>>>>>> 47c502e8e4318a844071ec8793822e1e8af02bf9
>>>>>>> 3064246b7a7818cc5e2250e0a2139fa8a7c21204
>>>>>>> 6487c29a178d340690563b4cc1ec43aceac6d944
>>>>>>> 77f9ed17803d5c660639c899c05ad4a9b1bdcd29
	//Counters
	G4int TotalCounter;
	G4int InelasticCounter;
	G4int ElasticCounter;
	//Optical Scintillation photon counter
	G4int totScintPhotons;
<<<<<<< HEAD
	//Plastic edep
	G4double fPEdep;
	G4double fPEkin;
=======
>>>>>>> 77f9ed17803d5c660639c899c05ad4a9b1bdcd29
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
