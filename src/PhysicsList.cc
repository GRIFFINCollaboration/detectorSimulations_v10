
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
// $Id: PhysicsList.cc 73290 2013-08-23 10:04:20Z gcosmo $
//
/// \file radioactivedecay/rdecay02/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "PhysListParticles.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "PhysListHadron.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4StepLimiter.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"

#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4DecayPhysics.hh"

#include "G4SystemOfUnits.hh"

#include "G4Threading.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4Scintillation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() :
	G4VModularPhysicsList(),
	fCutForGamma(1.*mm), fCutForElectron(0.01*mm),//e- cut to follow old code
	fCutForPositron(1.*mm),//just sets defaults, can be altered through commands in a macro
	fEmPhysicsList(0),
	fRaddecayList(0),
	fParticleList(0),
	fHadPhysicsList(0),
	fNhadcomp(0),
	fPMessenger(0),fDetectorCuts(0), fTargetCuts(0),
	fScintProcess(0)
{
	G4LossTableManager::Instance();
	defaultCutValue =0.1*mm;

	fPMessenger = new PhysicsListMessenger(this);

	SetVerboseLevel(1);

	//default physics
	fParticleList = new G4DecayPhysics();

	//default radioactive physics
	fRaddecayList = new G4RadioactiveDecayPhysics();

	// EM physics
	fEmPhysicsList = new G4EmStandardPhysics();

	// Scintillation phyics
	fScintProcess = new G4Scintillation();
	fScintProcess->SetScintillationYieldFactor(1.);
	fScintProcess->SetTrackSecondariesFirst(true);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
	delete fPMessenger;
	delete fRaddecayList;
	delete fEmPhysicsList;
	if(fHadPhysicsList) delete fHadPhysicsList;
	if(fNhadcomp > 0) {
		for(G4int i=0; i<fNhadcomp; i++) {
			delete fHadronPhys[i];
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
	fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
	AddTransportation();
	// em
	fEmPhysicsList->ConstructProcess();
	// decays
	fParticleList->ConstructProcess();
	fRaddecayList->ConstructProcess();
	// had
	if(fNhadcomp > 0) {
		for(G4int i=0; i<fNhadcomp; i++) {
			(fHadronPhys[i])->ConstructProcess();
		}
	}
	if(fHadPhysicsList) fHadPhysicsList->ConstructProcess();
	G4cout<<"### PhysicsList::ConstructProcess is done"<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SelectPhysicsList(const G4String& name)
{
	if(verboseLevel>1) {
		G4cout<<"PhysicsList::SelectPhysicsList: <"<<name<<">"<<G4endl;
	}
	// default  Had physics
	if(name == "Hadron" && !fHadPhysicsList) {
		fHadPhysicsList = new PhysListHadron("hadron");
	} else if(name == "QGSP_BERT") {
		AddExtraBuilders(false);
		fHadPhysicsList = new G4HadronPhysicsQGSP_BERT(verboseLevel);
	} else if(name == "QGSP_BIC" && !fHadPhysicsList) {
		AddExtraBuilders(false);
		fHadPhysicsList = new G4HadronPhysicsQGSP_BIC(verboseLevel);
	} else if(name == "QGSP_BERT_HP"  && !fHadPhysicsList) {
		AddExtraBuilders(true);
		fHadPhysicsList = new G4HadronPhysicsQGSP_BERT_HP(verboseLevel);
	} else if(name == "QGSP_BIC_HP"  && !fHadPhysicsList) {
		AddExtraBuilders(true);
		fHadPhysicsList = new G4HadronPhysicsQGSP_BIC_HP(verboseLevel);
	} else if (name == "FTFP_BERT_HP" && !fHadPhysicsList) {
		AddExtraBuilders(true);
		fHadPhysicsList = new G4HadronPhysicsFTFP_BERT_HP(verboseLevel);
	} else if(name == "emlivermore") {
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmLivermorePhysics(verboseLevel);
	} else if(name == "emlivermorepolarized") {
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmLivermorePolarizedPhysics(verboseLevel);
	} else if(name == "empenelope") {
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmPenelopePhysics(verboseLevel);
	} else if(name == "emstandard") {
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmStandardPhysics(verboseLevel);
	} else if(name == "emstandard_opt4") {
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmStandardPhysics_option4(verboseLevel);
	} else {
		G4cout<<"PhysicsList WARNING wrong or unknown <"
			<< name<<"> Physics "<<G4endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PhysicsList::ConstructOp(G4bool constructOp)
{
	if(constructOp == false) { G4cout<<"NOT Building Optical Phyiscs"<<G4endl; return; }
	else if(constructOp == true) 
	{
		G4cout<<"Building Optical Phyiscs... "<<G4endl;
		// G4Cerenkov* cerenkovProcess = new G4Cerenkov("Cerenkov");

		//   turning on particle-specific scintillation process
		G4Scintillation* scintillationProcess = new G4Scintillation("Scintillation");
		scintillationProcess->SetScintillationByParticleType(true);

		G4OpAbsorption* absorptionProcess = new G4OpAbsorption();
		G4OpRayleigh* rayleighScatteringProcess = new G4OpRayleigh();
		//G4OpMieHG* mieHGScatteringProcess = new G4OpMieHG();
		//G4OpBoundaryProcess* boundaryProcess = new G4OpBoundaryProcess();

		// Use Birks Correction in the Scintillation process
		if(!G4Threading::IsWorkerThread())
		{
			G4EmSaturation* emSaturation =
				G4LossTableManager::Instance()->EmSaturation();
			fScintProcess->AddSaturation(emSaturation);
		}

		auto aParticleIterator = GetParticleIterator();
		aParticleIterator->reset();
		while((*aParticleIterator)())
		{
			G4ParticleDefinition* particle = aParticleIterator->value();
			G4ProcessManager* pmanager = particle->GetProcessManager();
			G4String particleName = particle->GetParticleName();

			if(scintillationProcess->IsApplicable(*particle)) {
				pmanager->AddProcess(scintillationProcess);
				pmanager->SetProcessOrderingToLast(scintillationProcess, idxAtRest);
				pmanager->SetProcessOrderingToLast(scintillationProcess, idxPostStep);
			}

			if(particleName == "opticalphoton") {
				G4cout<<"AddDiscreteProcess to OpticalPhoton "<<G4endl;
				pmanager->AddDiscreteProcess(absorptionProcess);
				pmanager->AddDiscreteProcess(rayleighScatteringProcess);
				// pmanager->AddDiscreteProcess(mieHGScatteringProcess);
				//  pmanager->AddDiscreteProcess(boundaryProcess);
			}
		}
		G4cout<<"Done Building Optical Physics"<<G4endl;
	}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SpiceStepper(G4bool step)
{
	if(step == false) { G4cout<<"No SPICE step"<<G4endl; return; }
	else if(step == true) 
	{
		G4ProcessManager* pmanager = G4Electron::Electron()->GetProcessManager();
		pmanager->AddDiscreteProcess(new G4StepLimiter);
		//Electrons will now use user defined maxstep set by LogicVolume->SetUserLimits(new G4UserLimits(1.*mm));

		G4cout<<"SPICE stepper"<<G4endl;
	}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddExtraBuilders(G4bool flagHP)
{
	fNhadcomp = 5;
	fHadronPhys.push_back(new G4EmExtraPhysics(verboseLevel));
	if(flagHP) {
		fHadronPhys.push_back(new G4HadronElasticPhysicsHP(verboseLevel));
	} else {
		fHadronPhys.push_back(new G4HadronElasticPhysics(verboseLevel));
	}
	fHadronPhys.push_back(new G4StoppingPhysics(verboseLevel));
	fHadronPhys.push_back(new G4IonBinaryCascadePhysics(verboseLevel));
	fHadronPhys.push_back(new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
	SetCutValue(fCutForGamma, "gamma");
	SetCutValue(fCutForElectron, "e-");
	SetCutValue(fCutForPositron, "e+");

	SetCutValue(1.0*nm, "proton");
	SetCutValue(1.0*nm, "deuteron");
	SetCutValue(1.0*nm, "C12");
	SetCutValue(1.0*um, "C13");
	SetCutValue(1.0*um, "Be9");
	SetCutValue(1.0*um, "B10");

	G4cout<<"world cuts are set"<<G4endl;

	//  if(!fTargetCuts) SetTargetCut(fCutForElectron);
	//  G4Region* region = (G4RegionStore::GetInstance())->GetRegion("Target");
	//  region->SetProductionCuts(fTargetCuts);
	//  G4cout<<"Target cuts are set"<<G4endl;

	//  if(!fDetectorCuts) SetDetectorCut(fCutForElectron);
	//  region = (G4RegionStore::GetInstance())->GetRegion("Detector");
	//  region->SetProductionCuts(fDetectorCuts);
	//  G4cout<<"Detector cuts are set"<<G4endl;

	if(verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
	fCutForGamma = cut;
	SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
	fCutForElectron = cut;
	SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
	fCutForPositron = cut;
	SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetTargetCut(G4double cut)
{
	if(!fTargetCuts) fTargetCuts = new G4ProductionCuts();

	fTargetCuts->SetProductionCut(cut, idxG4GammaCut);
	fTargetCuts->SetProductionCut(cut, idxG4ElectronCut);
	fTargetCuts->SetProductionCut(cut, idxG4PositronCut);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetDetectorCut(G4double cut)
{
	if(!fDetectorCuts) fDetectorCuts = new G4ProductionCuts();

	fDetectorCuts->SetProductionCut(cut, idxG4GammaCut);
	fDetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
	fDetectorCuts->SetProductionCut(cut, idxG4PositronCut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
