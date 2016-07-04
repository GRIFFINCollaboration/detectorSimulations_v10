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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcess.hh"
#include "G4ProcessType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* detcon,
                               EventAction* evt)
    : G4UserSteppingAction(),
      fDetector(detcon), fEventAction(evt)
{
    fGriffinDetectorMapSet = false;
    fNumberOfAssemblyVols = 13;
    fStepNumber = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep) {
    G4bool trackSteps   = false;
    G4int particleType  = 0;
    G4int processType   = 0;
    G4int systemID      = 9999;
    G4int evntNb;

    fDet = 0;
    fCry = 0;

    G4String particleName;
    G4String mnemonic = "XXX00XX00X";

    // Get volume of the current step
    G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    G4String volname = volume->GetName();

    // collect energy and track length step by step
    // As it's called more than once, get the Track and assign to variable
    G4double edep = aStep->GetTotalEnergyDeposit();
    G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();

    G4Track* theTrack = aStep->GetTrack();
    G4double stepl = 0.;
    if (theTrack->GetDefinition()->GetPDGCharge() != 0.)
        stepl = aStep->GetStepLength();

    fStepNumber = theTrack->GetCurrentStepNumber();

    // Track particle type in EVERY step
    //G4cout << "Particle name = " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << G4endl;
    particleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    if (particleName == "gamma")         particleType = 1;
    else if (particleName == "e-")       particleType = 2;
    else if (particleName == "e+")       particleType = 3;
    else if (particleName == "proton")   particleType = 4;
    else if (particleName == "neutron")  particleType = 5;
    else if (particleName == "deuteron") particleType = 6;
    else if (particleName == "C12")      particleType = 7;
    else particleType = 0;

	 const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
	 G4int targetZ = -1;
	 if(process != NULL && process->GetProcessType() == fHadronic) {
		G4HadronicProcess* hadrProcess = (G4HadronicProcess*) process;
		const G4Isotope* target = NULL;
		target = hadrProcess->GetTargetIsotope();
		if(target != NULL) {
		  //	G4cout<<particleName<<", "<<process->GetProcessName()<<" on "<<target->GetName()<<G4endl;
		  targetZ = target->GetZ();
		}
	 }

    // this can be modified to add more processes
	 if(theTrack->GetCreatorProcess() != NULL) {
		G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
		//G4cout<<"found secondary, particle "<<particleName<<", creation process "<<process<<G4endl;
      if (processName == "RadioactiveDecay")      processType = 1;
      else if (processName == "eIoni")            processType = 2;
      else if (processName == "msc")              processType = 3;
      else if (processName == "Scintillation")    processType = 4;
      else if (processName == "Cerenkov")         processType = 5;
      else  processType = 0;
	 } else {
		processType = -1;
	 }

    fEventAction->AddParticleType(particleType);
    evntNb =  fEventAction->GetEventNumber();

    //G4cout << "Found Edep = " << edep/keV << " keV in " << volname << G4endl;
    // example volname
    //volname = av1_impr6_sodiumIodideCrystalBlockLogPv0

    // Get initial momentum direction & energy of particle
    G4int trackID = theTrack->GetTrackID();
    G4int parentID = theTrack->GetParentID();

    G4StepPoint* point1 = aStep->GetPreStepPoint();
    G4StepPoint* point2 = aStep->GetPostStepPoint();

    G4ThreeVector pos1 = point1->GetPosition();
    G4ThreeVector pos2 = point2->GetPosition();

    //G4double time1 = point1->GetGlobalTime();
    G4double time2 = point2->GetGlobalTime();

    size_t found;

    // Griffin energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("germaniumBlock1");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinCrystDet(edep,stepl,fDet-1,fCry-1);
        mnemonic.replace(0,3,"GRG");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"G");
        systemID = 1000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("backQuarterSuppressor");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorBackDet(edep,stepl,fDet-1,fCry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"E");
        systemID = 1050;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("leftSuppressorExtension");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorLeftExtensionDet(edep,stepl,fDet-1,fCry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"A");
        systemID = 1010;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("rightSuppressorExtension");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorRightExtensionDet(edep,stepl,fDet-1,fCry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"B");
        systemID = 1020;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("leftSuppressorCasing");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorLeftSideDet(edep,stepl,fDet-1,fCry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"C");
        systemID = 1030;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("rightSuppressorCasing");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorRightSideDet(edep,stepl,fDet-1,fCry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"D");
        systemID = 1040;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // Dead layer specific code
    found = volname.find("germaniumDlsBlock1");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForDeadLayerSpecificGriffinCrystal(volname);
        fEventAction->AddGriffinCrystDet(edep,stepl,fDet-1,fCry-1);
        mnemonic.replace(0,3,"GRG");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"G");
        systemID = 1000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // LaBr detector energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("lanthanumBromideCrystalBlock");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        fEventAction->AddLaBrCrystDet(edep,stepl,fDet-1);
        mnemonic.replace(0,3,"LAB");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        systemID = 2000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // Ancillary BGO detector energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("ancillaryBgoBlock");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForAncillaryBGODetector(volname);
        fEventAction->AncillaryBgoDet(edep,stepl,fDet-1);
        mnemonic.replace(0,3,"ABG");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        systemID = 3000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // NaI energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("sodiumIodideCrystalBlock");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"NAI");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        systemID = 4000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // Sceptar energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("sceptarSquareScintillatorLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        //if(fDet >= 6) fDet = fDet + 5; // to number SCPETAR paddles correctly

        if(fDet == 1) fDet = 6;
        else if(fDet == 2) fDet = 10;
        else if(fDet == 3) fDet = 9;
        else if(fDet == 4) fDet = 8;
        else if(fDet == 5) fDet = 7;
        else if(fDet == 6) fDet = 14;
        else if(fDet == 7) fDet = 13;
        else if(fDet == 8) fDet = 12;
        else if(fDet == 9) fDet = 11;
        else if(fDet == 10) fDet = 15;
        fEventAction->SceptarDet(edep,stepl,fDet-1);
        mnemonic.replace(0,3,"SCP");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        systemID = 5000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }
    found = volname.find("sceptarAngledScintillatorLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        if(fDet == 1) fDet = 1; // to number SCPETAR paddles correctly
        else if(fDet == 2) fDet = 5;
        else if(fDet == 3) fDet = 4;
        else if(fDet == 4) fDet = 3;
        else if(fDet == 5) fDet = 2;
        else if(fDet == 6) fDet = 19;
        else if(fDet == 7) fDet = 18;
        else if(fDet == 8) fDet = 17;
        else if(fDet == 9) fDet = 16;
        else if(fDet == 10) fDet = 20;
        fEventAction->SceptarDet(edep,stepl,fDet-1);
        mnemonic.replace(0,3,"SCP");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        systemID = 5000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // 8PI energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("8piGermaniumBlockLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        fEventAction->Add8piCrystDet(edep,stepl,fDet-1);
        mnemonic.replace(0,3,"8PI");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        systemID = 6000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("8piInnerBGOAnnulus");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"8PI");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"A");
        systemID = 6010;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("8piOuterLowerBGOAnnulus");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"8PI");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"B");
        systemID = 6020;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("8piOuterUpperBGOAnnulus");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"8PI");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        mnemonic.replace(6,1,"C");
        systemID = 6030;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("gridcellLog");
    if (ekin != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        if(particleType == 1)fEventAction->AddGridCellGamma(ekin,stepl,fDet-1);
        if(particleType == 2)fEventAction->AddGridCellElectron(ekin,stepl,fDet-1);
        if(particleType == 5)fEventAction->AddGridCellNeutron(ekin,stepl,fDet-1);
        mnemonic.replace(0,3,"GRD");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        systemID = 7000;
        // edep is ekin in this case. It would be useful if we also had the momentum positiin of the particles...
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, ekin, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
        // Now kill the track!
        theTrack->SetTrackStatus(fStopAndKill);
    }

    // DESCANT Detectors ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("blueScintillatorVolumeLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8010;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("greenScintillatorVolumeLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8020;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("redScintillatorVolumeLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8030;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("whiteScintillatorVolumeLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8040;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("yellowScintillatorVolumeLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8050;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("pacesSiliconBlockLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        fEventAction->Add8piCrystDet(edep,stepl,fDet-1);
        mnemonic.replace(0,3,"PCS");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        systemID = 9000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("testcanScintillatorLog");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"XXX");
        mnemonic.replace(3,2,G4intToG4String(fDet));
        mnemonic.replace(5,1,GetCrystalColour(fCry));
        systemID = 8500;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }


    if(trackSteps) {
        fEventAction->AddStepTracker(evntNb, trackID, parentID, fStepNumber, particleType, processType, systemID, fCry-1, fDet-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
	 }
}

void SteppingAction::SetDetAndCryNumberForGriffinComponent(G4String volname) {
    const char *cstr = volname.c_str();
    G4int av;
    G4int impr;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if( avOver9 == 47 ) { // under 10
        av = cstr[3]-'0';
        impr = cstr[10]-'0';
    }
    else if( avOver99 == 47 ) { // under 100
        av = (cstr[3]-'0')*10+(cstr[4]-'0');
        impr = cstr[11]-'0';
    }
    else { // OVER 100
        av = (cstr[3]-'0')*100+(cstr[4]-'0')*10+(cstr[5]-'0');  // This was fixed
        impr = cstr[12]-'0';
    }

    fDet = (G4int)(ceil(((G4double)(av)-5.0)/(G4double)(fNumberOfAssemblyVols))); // This was fixed
    fCry = impr;

    fDet = FindTrueGriffinDetector(fDet);

    //G4cout << "Found Edep in " << volname <<  " fCry = " << fCry << " fDet = " << fDet << " av = " << av << G4endl;
}

void SteppingAction::SetDetAndCryNumberForDeadLayerSpecificGriffinCrystal(G4String volname) {
    const char *cstr = volname.c_str();
    G4int av;
    G4int impr;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if(avOver9 == 47) { // under 10
        av = cstr[3]-'0';
        impr = cstr[10]-'0';
    }
    else if(avOver99 == 47) { // under 100
        av = (cstr[3]-'0')*10+(cstr[4]-'0');
        impr = cstr[11]-'0';
    }
    else { // OVER 100
        av = (cstr[3]-'0')*100+(cstr[4]-'0')*10+(cstr[5]-'0');
        impr = cstr[12]-'0';
    }

    fDet = (G4int)(ceil((G4double)(av)/(G4double)(fNumberOfAssemblyVols)));
    fCry = av - fNumberOfAssemblyVols*(fDet-1);

    fDet = FindTrueGriffinDetector(fDet);

    //G4cout << "Found Edep in " << volname <<  " fCry = " << fCry << " fDet = " << fDet << " av = " << av << G4endl;
}

void SteppingAction::SetDetNumberForGenericDetector(G4String volname) {
    const char *cstr = volname.c_str();
    G4int volNameOver9;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if(avOver9 == 47) { // under 10
        volNameOver9 = cstr[11]-'0';
        if(volNameOver9 == 47) {
            fDet = cstr[10]-'0';
        }
        else {
            fDet = ((cstr[10]-'0')*10)+volNameOver9 ;
        }
    }
    else if(avOver99 == 47) { // under 100
        volNameOver9 = cstr[12]-'0';
        if(volNameOver9 == 47) {
            fDet = cstr[11]-'0';
        }
        else {
            fDet = ((cstr[11]-'0')*10)+volNameOver9 ;
        }
    }
    else { // OVER 100
        volNameOver9 = cstr[13]-'0';
        if(volNameOver9 == 47) {
            fDet = cstr[12]-'0';
        }
        else {
            fDet = ((cstr[12]-'0')*10)+volNameOver9 ;
        }
    }
    //G4cout << "Stepping Action :: Found electron ekin in " << volname << " fDet = " << fDet << G4endl;
}


void SteppingAction::SetDetNumberForAncillaryBGODetector(G4String volname) {
    const char *cstr = volname.c_str();
    G4int av;
    G4int impr;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if( avOver9 == 47 ) { // under 10
        av = cstr[3]-'0';
        if((cstr[12]-'0') == 47) { // impr > 10 and < 100
            impr = (cstr[10]-'0')*10+(cstr[11]-'0');
        }
        else { // assume < 10
            impr = cstr[10]-'0';
        }
    }
    else if( avOver99 == 47 ) { // under 100
        av = (cstr[3]-'0')*10+(cstr[4]-'0');
        if((cstr[13]-'0') == 47) { // impr > 10 and < 100
            impr = (cstr[11]-'0')*10+(cstr[12]-'0');
        }
        else { // assume < 10
            impr = cstr[11]-'0';
        }
    }
    else { // OVER 100
        av = (cstr[3]-'0')*100+(cstr[4]-'0')*10+(cstr[5]-'0');  // This was fixed
        if((cstr[14]-'0') == 47) { // impr > 10 and < 100
            impr = (cstr[12]-'0')*10+(cstr[13]-'0');
        }
        else { // assume < 10
            impr = cstr[12]-'0';
        }
    }
    fDet = (G4int)((ceil)(G4double(impr)/3.0));
    fCry = impr-((fDet-1)*3);

    //    G4cout << "Found Edep in " << volname <<  " fCry = " << fCry << " fDet = " << fDet << " av = " << av << G4endl;
}

G4int SteppingAction::FindTrueGriffinDetector(G4int detval) {
    G4int trueDet;
    trueDet = fDetector->fGriffinDetectorsMap[detval-1];

    return trueDet;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String SteppingAction::G4intToG4String(G4int value) {
    G4String theString;
    G4String output = "00";
    std::stringstream out;
    out << value;
    theString = out.str();

    if(value < 10) {
        output.replace(1,1,theString);
    }
    else {
        output = theString;
    }

    return output;
}

G4String SteppingAction::GetCrystalColour(G4int value) {
    G4String output = "X";
    if(value == 1) {
        output = "B";
    }
    else if(value == 2) {
        output = "G";
    }
    else if(value == 3) {
        output = "R";
    }
    else if(value == 4) {
        output = "W";
    }
    return output;
}
