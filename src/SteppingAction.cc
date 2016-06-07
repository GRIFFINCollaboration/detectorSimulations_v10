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
    griffinDetectorMapSet = false;
    numberOfAssemblyVols = 13;
    stepNumber = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    G4bool trackSteps   = false;
    G4int particleType  = 0;
    G4int processType   = 0;
    G4int systemID      = 9999;
    G4int evntNb;

    det = 0;
    cry = 0;

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

    stepNumber = theTrack->GetCurrentStepNumber();

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
		G4HadronicProcess* hadr_process = (G4HadronicProcess*) process;
		const G4Isotope* target = NULL;
		target = hadr_process->GetTargetIsotope();
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
    //volname = av_1_impr_6_sodium_iodide_crystal_block_log_pv_0

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
    found = volname.find("germanium_block1");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinCrystDet(edep,stepl,det-1,cry-1);
        mnemonic.replace(0,3,"GRG");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"G");
        systemID = 1000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("back_quarter_suppressor");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorBackDet(edep,stepl,det-1,cry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"E");
        systemID = 1050;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("left_suppressor_extension");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorLeftExtensionDet(edep,stepl,det-1,cry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"A");
        systemID = 1010;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("right_suppressor_extension");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorRightExtensionDet(edep,stepl,det-1,cry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"B");
        systemID = 1020;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("left_suppressor_casing");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorLeftSideDet(edep,stepl,det-1,cry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"C");
        systemID = 1030;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("right_suppressor_casing");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForGriffinComponent(volname);
        fEventAction->AddGriffinSuppressorRightSideDet(edep,stepl,det-1,cry-1);
        mnemonic.replace(0,3,"GRS");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"D");
        systemID = 1040;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // Dead layer specific code
    found = volname.find("germanium_dls_block1");
    if (edep != 0 && found!=G4String::npos) {
        SetDetAndCryNumberForDeadLayerSpecificGriffinCrystal(volname);
        fEventAction->AddGriffinCrystDet(edep,stepl,det-1,cry-1);
        mnemonic.replace(0,3,"GRG");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"G");
        systemID = 1000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // LaBr detector energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("lanthanum_bromide_crystal_block");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        fEventAction->AddLaBrCrystDet(edep,stepl,det-1);
        mnemonic.replace(0,3,"LAB");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        systemID = 2000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // Ancillary BGO detector energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("ancillary_bgo_block");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForAncillaryBGODetector(volname);
        fEventAction->AncillaryBgoDet(edep,stepl,det-1);
        mnemonic.replace(0,3,"ABG");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        systemID = 3000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // NaI energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("sodium_iodide_crystal_block");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"NAI");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        systemID = 4000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // Sceptar energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("sceptar_square_scintillator_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        //if(det >= 6) det = det + 5; // to number SCPETAR paddles correctly

        if(det == 1) det = 6;
        else if(det == 2) det = 10;
        else if(det == 3) det = 9;
        else if(det == 4) det = 8;
        else if(det == 5) det = 7;
        else if(det == 6) det = 14;
        else if(det == 7) det = 13;
        else if(det == 8) det = 12;
        else if(det == 9) det = 11;
        else if(det == 10) det = 15;
        fEventAction->SceptarDet(edep,stepl,det-1);
        mnemonic.replace(0,3,"SCP");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        systemID = 5000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }
    found = volname.find("sceptar_angled_scintillator_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        if(det == 1) det = 1; // to number SCPETAR paddles correctly
        else if(det == 2) det = 5;
        else if(det == 3) det = 4;
        else if(det == 4) det = 3;
        else if(det == 5) det = 2;
        else if(det == 6) det = 19;
        else if(det == 7) det = 18;
        else if(det == 8) det = 17;
        else if(det == 9) det = 16;
        else if(det == 10) det = 20;
        fEventAction->SceptarDet(edep,stepl,det-1);
        mnemonic.replace(0,3,"SCP");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        systemID = 5000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    // 8PI energy deposits ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("8pi_germanium_block_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        fEventAction->Add8piCrystDet(edep,stepl,det-1);
        mnemonic.replace(0,3,"8PI");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        systemID = 6000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("8pi_inner_BGO_annulus");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"8PI");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"A");
        systemID = 6010;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("8pi_outer_lower_BGO_annulus");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"8PI");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"B");
        systemID = 6020;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("8pi_outer_upper_BGO_annulus");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"8PI");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        mnemonic.replace(6,1,"C");
        systemID = 6030;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("gridcell_log");
    if (ekin != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        if(particleType == 1)fEventAction->AddGridCellGamma(ekin,stepl,det-1);
        if(particleType == 2)fEventAction->AddGridCellElectron(ekin,stepl,det-1);
        if(particleType == 5)fEventAction->AddGridCellNeutron(ekin,stepl,det-1);
        mnemonic.replace(0,3,"GRD");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        systemID = 7000;
        // edep is ekin in this case. It would be useful if we also had the momentum positiin of the particles...
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, ekin, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
        // Now kill the track!
        theTrack->SetTrackStatus(fStopAndKill);
    }

    // DESCANT Detectors ////////////////////////////////////////////////////////////////////////////////
    found = volname.find("blue_scintillator_volume_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8010;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("green_scintillator_volume_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8020;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("red_scintillator_volume_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8030;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("white_scintillator_volume_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8040;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("yellow_scintillator_volume_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"DSC");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
		  if(targetZ < 10) mnemonic.replace(6,1,G4intToG4String(targetZ));
        systemID = 8050;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("paces_silicon_block_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        fEventAction->Add8piCrystDet(edep,stepl,det-1);
        mnemonic.replace(0,3,"PCS");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        systemID = 9000;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }

    found = volname.find("testcan_scintillator_log");
    if (edep != 0 && found!=G4String::npos) {
        SetDetNumberForGenericDetector(volname);
        mnemonic.replace(0,3,"XXX");
        mnemonic.replace(3,2,G4intToG4String(det));
        mnemonic.replace(5,1,GetCrystalColour(cry));
        systemID = 8500;
        fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
    }


    //  // gamma angular correlations in world
    //  found = volname.find("World");
    //  if (ekin != 0 && found!=G4String::npos && particleType == 1) {
    //      if(ekin <= (100.+2.)*keV && ekin >= (100.-2.)*keV) {
    //          det = 1;
    //      }
    //      else if(ekin <= (200.+2.)*keV && ekin >= (200.-2.)*keV) {
    //          det = 2;
    //      }
    //      else if(ekin <= (300.+2.)*keV && ekin >= (300.-2.)*keV) {
    //          det = 3;
    //      }
    //      else if(ekin <= (400.+2.)*keV && ekin >= (400.-2.)*keV) {
    //          det = 4;
    //      }
    //      else if(ekin <= (500.+2.)*keV && ekin >= (500.-2.)*keV) {
    //          det = 5;
    //      }
    //      else if(ekin <= (600.+2.)*keV && ekin >= (600.-2.)*keV) {
    //          det = 6;
    //      }
    //      else if(ekin <= (700.+2.)*keV && ekin >= (700.-2.)*keV) {
    //          det = 7;
    //      }
    //      else if(ekin <= (800.+2.)*keV && ekin >= (800.-2.)*keV) {
    //          det = 8;
    //      }
    //      else if(ekin <= (900.+2.)*keV && ekin >= (900.-2.)*keV) {
    //          det = 9;
    //      }
    //      else if(ekin <= (1173.2+2.)*keV && ekin >= (1173.2-2.)*keV) { //60Co decay
    //          det = 10;
    //      }
    //      else if(ekin <= (1332.5+2.)*keV && ekin >= (1332.5-2.)*keV) {
    //          det = 11;
    //      }
    //      else if(ekin <= (1808.660+2.)*keV && ekin >= (1808.660-2.)*keV) { //26Na decay
    //          det = 12;
    //      }
    //      else if(ekin <= (1896.720+2.)*keV && ekin >= (1896.720-2.)*keV) {
    //          det = 13;
    //      }
    //      else if(ekin <= (2541.220+2.)*keV && ekin >= (2541.220-2.)*keV) {
    //          det = 14;
    //      }
    //      else if(ekin <= (1129.580+2.)*keV && ekin >= (1129.580-2.)*keV) {
    //          det = 15;
    //      }
    //      else if(ekin <= (953.8+2.)*keV && ekin >= (953.8-2.)*keV) { //62Ga decay
    //          det = 16;
    //      }
    //      else if(ekin <= (850.8+2.)*keV && ekin >= (850.8-2.)*keV) {
    //          det = 17;
    //      }
    //      else if(ekin <= (1388.300+2.)*keV && ekin >= (1388.300-2.)*keV) {
    //          det = 18;
    //      }
    //      else {
    //          det = 0;
    //      }
    //      cry = 1;
    //      systemID = 9999;
    //      mnemonic.replace(0,3,"GAC");
    //      mnemonic.replace(3,2,G4intToG4String(det));
    //      mnemonic.replace(5,1,GetCrystalColour(cry));
    //      fEventAction->AddHitTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, ekin, pos2.x(), pos2.y(), pos2.z(), time2);
    //  }

    if(trackSteps) {
        fEventAction->AddStepTracker(mnemonic, evntNb, trackID, parentID, stepNumber, particleType, processType, systemID, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, targetZ);
	 }
}

void SteppingAction::SetDetAndCryNumberForGriffinComponent(G4String volname)
{
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

    det = (G4int)(ceil(((G4double)(av)-5.0)/(G4double)(numberOfAssemblyVols))); // This was fixed
    cry = impr;

    det = FindTrueGriffinDetector(det);

    //G4cout << "Found Edep in " << volname <<  " cry = " << cry << " det = " << det << " av = " << av << G4endl;
}

void SteppingAction::SetDetAndCryNumberForDeadLayerSpecificGriffinCrystal(G4String volname)
{
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

    det = (G4int)(ceil((G4double)(av)/(G4double)(numberOfAssemblyVols)));
    cry = av - numberOfAssemblyVols*(det-1);

    det = FindTrueGriffinDetector(det);

    //G4cout << "Found Edep in " << volname <<  " cry = " << cry << " det = " << det << " av = " << av << G4endl;
}

void SteppingAction::SetDetNumberForGenericDetector(G4String volname)
{
    const char *cstr = volname.c_str();
    G4int volNameOver9;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if(avOver9 == 47) { // under 10
        volNameOver9 = cstr[11]-'0';
        if(volNameOver9 == 47) {
            det = cstr[10]-'0';
        }
        else {
            det = ((cstr[10]-'0')*10)+volNameOver9 ;
        }
    }
    else if(avOver99 == 47) { // under 100
        volNameOver9 = cstr[12]-'0';
        if(volNameOver9 == 47) {
            det = cstr[11]-'0';
        }
        else {
            det = ((cstr[11]-'0')*10)+volNameOver9 ;
        }
    }
    else { // OVER 100
        volNameOver9 = cstr[13]-'0';
        if(volNameOver9 == 47) {
            det = cstr[12]-'0';
        }
        else {
            det = ((cstr[12]-'0')*10)+volNameOver9 ;
        }
    }
    //G4cout << "Stepping Action :: Found electron ekin in " << volname << " det = " << det << G4endl;
}


void SteppingAction::SetDetNumberForAncillaryBGODetector(G4String volname)
{
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
    det = (G4int)((ceil)(G4double(impr)/3.0));
    cry = impr-((det-1)*3);

    //    G4cout << "Found Edep in " << volname <<  " cry = " << cry << " det = " << det << " av = " << av << G4endl;
}

G4int SteppingAction::FindTrueGriffinDetector(G4int detval)
{
    G4int trueDet;
    trueDet = fDetector->griffinDetectorsMap[detval-1];

    return trueDet;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String SteppingAction::G4intToG4String(G4int value)
{
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

G4String SteppingAction::GetCrystalColour(G4int value)
{
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
