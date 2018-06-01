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
// $Id: PhysListHadron.cc 82660 2014-07-02 08:43:20Z gcosmo $
//
/// \file radioactivedecay/rdecay02/src/PhysListHadron.cc
/// \brief Implementation of the PhysListHadron class
//

#include "PhysListHadron.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
//#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4HadronElastic.hh"
#include "G4LFission.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4CascadeInterface.hh"

#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4NuclNuclDiffuseElastic.hh"

//HPNeutron

#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
//c-s
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
//-------------------------------------
//#include "G4GGNuclNuclCrossSection.hh"
#include "G4ComponentGGNuclNuclXsc.hh" 
//-------------------------------------
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

// RadioactiveDecay
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"

#include "G4SystemOfUnits.hh"


PhysListHadron::PhysListHadron(const G4String& name)
: G4VPhysicsConstructor(name),
	fTheNeutronElasticProcess(0), fTheFissionProcess(0),
	fTheCaptureProcess(0),fTheDeuteronInelasticProcess(0),
	fTheTritonInelasticProcess(0), fTheAlphaInelasticProcess(0),
	fTheIonInelasticProcess(0)
{}

PhysListHadron::~PhysListHadron()
{}


void PhysListHadron::ConstructProcess()
{
	G4ProcessManager* pManager = 0;

	G4TheoFSGenerator* theTheoModel = new G4TheoFSGenerator;
	// all models for treatment of thermal nucleus
	G4Evaporation* theEvaporation = new G4Evaporation;
	//G4FermiBreakUp* theFermiBreakUp = new G4FermiBreakUp;
	G4StatMF* theMF = new G4StatMF;
	// Evaporation logic
	G4ExcitationHandler* theHandler = new G4ExcitationHandler;
	theHandler->SetEvaporation(theEvaporation);
	//theHandler->SetFermiModel(theFermiBreakUp);
	theHandler->SetMultiFragmentation(theMF);
	theHandler->SetMaxAandZForFermiBreakUp(12, 6);
	theHandler->SetMinEForMultiFrag(5*MeV);
	// Pre equilibrium stage
	G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(theHandler);

	// a no-cascade generator-precompound interaface
	G4GeneratorPrecompoundInterface* theCascade =
		new G4GeneratorPrecompoundInterface;
	theCascade->SetDeExcitation(thePreEquilib);

	// High energy parts
	// String model; still not quite according to design
	// Explicit use of the forseen interfaces
	G4VPartonStringModel* theStringModel;
	theStringModel = new G4QGSModel<G4QGSParticipants>;
	theTheoModel->SetTransport(theCascade);
	theTheoModel->SetHighEnergyGenerator(theStringModel);
	theTheoModel->SetMinEnergy(10*GeV);  // 15 GeV may be the right limit
	theTheoModel->SetMaxEnergy(100*TeV);

	G4VLongitudinalStringDecay * theFragmentation = new G4QGSMFragmentation;
	G4ExcitedStringDecay * theStringDecay =
		new G4ExcitedStringDecay(theFragmentation);
	theStringModel->SetFragmentationModel(theStringDecay);

	// Elastic Process

	fTheElasticProcess.RegisterMe(new G4HadronElastic());

	// Hadron elastic process for all particles except neutrons and generic ions
	auto aParticleIterator = GetParticleIterator();
	aParticleIterator->reset();
	while ((*aParticleIterator)()) {
		G4ParticleDefinition* particle = aParticleIterator->value();
		G4String particleName = particle->GetParticleName();
		if(particleName != "neutron" && particleName != "GenericIon" &&
				particleName != "He3") {
			pManager = particle->GetProcessManager();
			if(particle->GetPDGMass() > 110.*MeV &&
					fTheElasticProcess.IsApplicable(*particle) &&
					!particle->IsShortLived()) {
				pManager->AddDiscreteProcess(&fTheElasticProcess);
				//        G4cout<<"### Elastic model are registered for "
				//      <<particle->GetParticleName()
				//      <<G4endl;
			}
		}
	}

	// Proton
	pManager = G4Proton::Proton()->GetProcessManager();
	// add inelastic process
	// Binary Cascade
	G4BinaryCascade * theBC = new G4BinaryCascade;
	theBC->SetMaxEnergy(10.5*GeV);
	fTheProtonInelastic.RegisterMe(theBC);
	// Higher energy
	fTheProtonInelastic.RegisterMe(theTheoModel);
	// now the cross-sections.
	G4ProtonInelasticCrossSection* theProtonData =
		new G4ProtonInelasticCrossSection;
	fTheProtonInelastic.AddDataSet(theProtonData);
	pManager->AddDiscreteProcess(&fTheProtonInelastic);

	// Neutron
	pManager = G4Neutron::Neutron()->GetProcessManager();
	// add process
	// elastic scattering
	fTheNeutronElasticProcess = new G4HadronElasticProcess();
	G4HadronElastic* theElasticModel1 = new G4HadronElastic;
	G4NeutronHPElastic * theElasticNeutron = new G4NeutronHPElastic;
	fTheNeutronElasticProcess->RegisterMe(theElasticModel1);
	theElasticModel1->SetMinEnergy(19.*MeV);
	fTheNeutronElasticProcess->RegisterMe(theElasticNeutron);
	theElasticNeutron->SetMaxEnergy(20.*MeV);

	G4NeutronHPElasticData * theNeutronData = new G4NeutronHPElasticData;
	fTheNeutronElasticProcess->AddDataSet(theNeutronData);
	pManager->AddDiscreteProcess(fTheNeutronElasticProcess);

	// inelastic
	G4NeutronHPInelastic * theHPNeutronInelasticModel =
		new G4NeutronHPInelastic;
	theHPNeutronInelasticModel->SetMaxEnergy(20.*MeV);
	fTheNeutronInelastic.RegisterMe(theHPNeutronInelasticModel);
	G4NeutronHPInelasticData * theNeutronData1 = new G4NeutronHPInelasticData;
	fTheNeutronInelastic.AddDataSet(theNeutronData1);
	// binary
	G4BinaryCascade * neutronBC = new G4BinaryCascade;
	neutronBC->SetMinEnergy(19.*MeV);
	neutronBC->SetMaxEnergy(10.5*GeV);
	fTheNeutronInelastic.RegisterMe(neutronBC);
	// higher energy
	fTheNeutronInelastic.RegisterMe(theTheoModel);
	// now the cross-sections.
	G4NeutronInelasticCrossSection * theNeutronData2 =
		new G4NeutronInelasticCrossSection;
	fTheNeutronInelastic.AddDataSet(theNeutronData2);
	pManager->AddDiscreteProcess(&fTheNeutronInelastic);

	// fission
	fTheFissionProcess = new G4HadronFissionProcess;
	G4LFission* theFissionModel = new G4LFission;
	fTheFissionProcess->RegisterMe(theFissionModel);
	pManager->AddDiscreteProcess(fTheFissionProcess);

	//capture
	fTheCaptureProcess = new G4HadronCaptureProcess;
	G4NeutronRadCapture* theCaptureModel = new G4NeutronRadCapture;
	theCaptureModel->SetMinEnergy(19.*MeV);
	fTheCaptureProcess->RegisterMe(theCaptureModel);
	fTheCaptureProcess->AddDataSet(new G4NeutronCaptureXS);

	G4NeutronHPCapture * theHPNeutronCaptureModel = new G4NeutronHPCapture;
	fTheCaptureProcess->RegisterMe(theHPNeutronCaptureModel);
	G4NeutronHPCaptureData * theNeutronData3 = new G4NeutronHPCaptureData;
	fTheCaptureProcess->AddDataSet(theNeutronData3);
	pManager->AddDiscreteProcess(fTheCaptureProcess);

	// now light ions
	// light Ion BC
	G4BinaryLightIonReaction* theIonBC= new G4BinaryLightIonReaction;
	theIonBC->SetMinEnergy(1*MeV);
	theIonBC->SetMaxEnergy(20*GeV);
	G4TripathiCrossSection * TripathiCrossSection= new G4TripathiCrossSection;
	G4IonsShenCrossSection * aShen = new G4IonsShenCrossSection;

	// deuteron
	pManager = G4Deuteron::Deuteron()->GetProcessManager();
	fTheDeuteronInelasticProcess = new G4DeuteronInelasticProcess("inelastic");
	fTheDeuteronInelasticProcess->AddDataSet(TripathiCrossSection);
	fTheDeuteronInelasticProcess->AddDataSet(aShen);
	fTheDeuteronInelasticProcess->RegisterMe(theIonBC);
	fTheDeuteronInelasticProcess->RegisterMe(theTheoModel);
	pManager->AddDiscreteProcess(fTheDeuteronInelasticProcess);

	// triton
	pManager = G4Triton::Triton()->GetProcessManager();
	fTheTritonInelasticProcess = new G4TritonInelasticProcess("inelastic");
	fTheTritonInelasticProcess->AddDataSet(TripathiCrossSection);
	fTheTritonInelasticProcess->AddDataSet(aShen);
	fTheTritonInelasticProcess->RegisterMe(theIonBC);
	fTheTritonInelasticProcess->RegisterMe(theTheoModel);
	pManager->AddDiscreteProcess(fTheTritonInelasticProcess);

	// alpha
	pManager = G4Alpha::Alpha()->GetProcessManager();
	fTheAlphaInelasticProcess = new G4AlphaInelasticProcess("inelastic");
	fTheAlphaInelasticProcess->AddDataSet(TripathiCrossSection);
	fTheAlphaInelasticProcess->AddDataSet(aShen);
	fTheAlphaInelasticProcess->RegisterMe(theIonBC);
	fTheAlphaInelasticProcess->RegisterMe(theTheoModel);
	pManager->AddDiscreteProcess(fTheAlphaInelasticProcess);

	// Generic ion elastic
	G4NuclNuclDiffuseElastic* ionElasticModel = new G4NuclNuclDiffuseElastic;

	// ----------------------
	//G4ComponentGGNuclNuclXsc * ionElasticXS = new G4ComponentGGNuclNuclXsc;
	G4VCrossSectionDataSet* ionElasticXS = new G4VCrossSectionDataSet;
	//G4GGNuclNuclCrossSection* ionElasticXS = new G4GGNuclNuclCrossSection;
	//--------------------------

	fIonElasticProcess.RegisterMe(ionElasticModel);
	fIonElasticProcess.AddDataSet(ionElasticXS);
	pManager = G4GenericIon::GenericIon()->GetProcessManager();
	pManager->AddDiscreteProcess(&fIonElasticProcess);

	// Generic ion inelastic
	fTheIonInelasticProcess = new G4IonInelasticProcess();
	fTheIonInelasticProcess->AddDataSet(TripathiCrossSection);
	fTheIonInelasticProcess->AddDataSet(aShen);
	fTheIonInelasticProcess->RegisterMe(theIonBC);
	fTheIonInelasticProcess->RegisterMe(theTheoModel);
	pManager->AddDiscreteProcess(fTheIonInelasticProcess);
}
