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
/// \file analysis/shared/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh" // for command-based data input

#include "Global.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Geantino.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh" //for detector based information
#include "BeamDistribution.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	PrimaryGeneratorAction::PrimaryGeneratorAction(HistoManager* histoManager)
: G4VUserPrimaryGeneratorAction(),
	fParticleGun(nullptr),
	fDetector(histoManager->GetDetectorConstruction()),
	fHistoManager(histoManager)
{
	G4int nParticle = 1;
	fParticleGun  = new G4ParticleGun(nParticle); //In our code, the gun is called fParticleGun
	//create a messenger for this class
	fGunMessenger = new PrimaryGeneratorMessenger(this);

	//these 3 lines initialise the Gun, basic values
	fParticleGun->SetParticleEnergy(0*eV);
	fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

	//defaults for gun properties
	fNumberOfDecayingLaBrDetectors = 0;
	fEffEnergy = 0.0;
	fEffDirection = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	fEffDirectionBool = false;//initialises bools, if command entered will go to true and loops below entered
	fEffPositionBool = false;
	fEffParticleBool = false;
	fEffPolarization = false;
	LaBrinit(); //sets up default variables - messy having them all declared here

	fBeamSpotSigma = 0.*mm;

	fTargetDistro=false;
	fNeedFileDistro = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fGunMessenger;
	delete fBeamDistribution;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//G4cout<<G4endl<<fParticleGun->GetParticleDefinition()->GetParticleName()<<G4endl;
	if(fNumberOfDecayingLaBrDetectors != 0) {
		G4double crystalRadius    = 2.54*cm;
		G4double crystalLength    = 5.08*cm;
		// The detector material material is LaBr3:Ce. 95% is LaBr3 and 5% is Ce.
		// Of the 95% that is LaBr3, there is one atom of La for every 3 of Br.
		// La has a molar mass of 138.9055*g/mole,
		// Br has a molar mass of 79.904*g/mole,
		// Ce has a molar mass of 140.116*g/mole.

		// We need to figure out how atoms of La are in our simulation.
		// Of the total number of La atoms, 0.08881(71)% are 139La and are radioative.
		// Calculated activity = 146.82908820569 Bq
		G4double prob, sumProb;
		G4int detnumber = 0;

		// Choose random detector
		// this is done in a weird way, why not just use numberOfDecayingLaBrDetectors*G4UniformRand() and cast it to a integer?
		prob = 1.0/((G4double)(fNumberOfDecayingLaBrDetectors));
		sumProb = 0.0;
		G4double randomDet = G4UniformRand();
		for(G4int j = 0 ; j < fNumberOfDecayingLaBrDetectors ; j++) {  // get the number of particles in decay and loop over them
			sumProb = sumProb + prob;
			if(randomDet <= sumProb) {
				detnumber = j;
				break;
			}
		}

		// Choose energy
		// 33.6% Beta Decay 788.742 keV Gamma-ray and 66.4% EC Decay 1435.795 keV Gamma-Ray
		prob = 0.336;
		sumProb = 0;
		G4double thisEnergy;
		G4double randomEnergy = G4UniformRand();

		if(randomEnergy <= prob) {
			thisEnergy = 788.742*keV;
		} else {
			thisEnergy = 1435.795*keV;
		}

		G4double detTheta       = fDetectorAnglesLaBr3[detnumber][0];
		G4double detPhi         = fDetectorAnglesLaBr3[detnumber][1];
		G4double detRadialPos   = 12.82385*cm;

		G4double randomZ   = crystalLength*G4UniformRand() + detRadialPos;
		G4double randomPhi = 2.0*(360.*deg)*G4UniformRand();
		G4double randomR   = crystalRadius*pow(G4UniformRand(),0.5);

		G4double x2=randomR*sin(randomPhi);
		G4double y2=randomR*cos(randomPhi);
		G4double z2=randomZ;

		G4ThreeVector pos2 = G4ThreeVector(x2,y2,z2);
		pos2.rotateY(M_PI);
		pos2.rotateY(M_PI+fDetectorAnglesLaBr3[detnumber][3]);
		pos2.rotateZ(fDetectorAnglesLaBr3[detnumber][4]);

		//cart
		G4double x = pos2.x();
		G4double y = pos2.y();
		G4double z = pos2.z();

		G4ThreeVector thisPosition(TransX(x,y,z,detTheta,detPhi), TransY(x,y,z,detTheta,detPhi), TransZ(x,y,z,detTheta));

		G4ParticleDefinition* agamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");

		// random direction
		G4double randcostheta = 2.*G4UniformRand()-1.0;
		G4double randsintheta = sqrt(1. - randcostheta*randcostheta);
		G4double randphi      = (360.*deg)*G4UniformRand();
		G4ThreeVector thisDirection = G4ThreeVector(randsintheta*cos(randphi), randsintheta*sin(randphi), randcostheta);

		fParticleGun->SetParticleDefinition(agamma);
		fParticleGun->SetParticlePosition(thisPosition);
		fParticleGun->SetParticleMomentumDirection(thisDirection);
		fParticleGun->SetParticleEnergy(thisEnergy);
	} else  {
		// Changed so that most grsi "/Detsys/gun/" commands still effect gun when using
		// Underlying geant4 commands such as '/gun/particle ion" & "/gun/ion"
		if(fEffParticleBool) {
			G4ParticleDefinition* effPart;
			if(fEffParticle == "electron" || fEffParticle == "e-") {
				effPart = G4ParticleTable::GetParticleTable()->FindParticle("e-");
			} else if(fEffParticle == "positron" || fEffParticle == "e+") {
				effPart = G4ParticleTable::GetParticleTable()->FindParticle("e+");
			} else if(fEffParticle == "neutron"){
				effPart = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
			} else {
				effPart = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
			}
			fParticleGun->SetParticleDefinition(effPart);
		}


		//////////
		//ORIGIN CONTROLS 
		//////////        

		G4double x = 0.*mm;
		G4double y = 0.*mm;
		G4double z = 0.*mm;
		if(fEffPositionBool){
			x = fEffPosition.x();
			y = fEffPosition.y();
			z = fEffPosition.z();	
		}
		G4double zZero = z;

		// If we want to simulate a realistic beam spot, instead of perfect pencil beam.
		if(fBeamSpotSigma>0){
			x = G4RandGauss::shoot(x,fBeamSpotSigma)*mm;
			y = G4RandGauss::shoot(y,fBeamSpotSigma)*mm;
		}

		if(fTargetDistro){
			z = fLayerStart+G4UniformRand()*fLayerLength;
		} else if(fNeedFileDistro){
			z += fBeamDistribution->GetRandom();
		}

		G4ThreeVector thisEffPosition = G4ThreeVector(x,y,z);//in constructor

		//////////
		//DIRECTION CONTROLS (directly forced, beam, or cone)
		//////////
		G4double effRandCosTheta, effRandSinTheta, effRandPhi;
		G4ThreeVector effdirection;
		if(fEffDirectionBool) { 
			effdirection = fEffDirection;

			if(fConeAngleBool){
				//min max input order doesnt actually matter
				G4double cmin =cos(fAngleMinInit);
				G4double cmax =cos(fAngleInit);
				G4double CosTheta = G4UniformRand()*abs(cmax-cmin);
				if(cmin<cmax)CosTheta+=cmin;
				else CosTheta+=cmax;

				// 	      G4cout<<asin(SinTheta)<<G4endl;
				G4double SinTheta = sqrt(1. - pow(CosTheta, 2.0));
				G4double Phi      = (2.0*CLHEP::pi)*G4UniformRand();

				effdirection = G4ThreeVector(SinTheta*cos(Phi), SinTheta*sin(Phi), CosTheta);
				// G4cout<<effdirection<<G4endl;

			}	  
		} else {
			//G4cout<<"Random "<< G4endl; //may offer the solution, an altered 2pi rando. Using 4pi for efficiency
			// random direction if no preference provided
			effRandCosTheta = 2.*G4UniformRand()-1.0; //cos(theta) = 2cos^2(0.5theta)-1 ??
			effRandSinTheta = sqrt(1. - effRandCosTheta*effRandCosTheta); //from sin^2(theta)+cos^2(theta)=1
			effRandPhi      = (360.*deg)*G4UniformRand();
			effdirection = G4ThreeVector(effRandSinTheta*cos(effRandPhi), effRandSinTheta*sin(effRandPhi), effRandCosTheta);
			//converts from Spherical polar(physics def.) to cartesian via (rsin(theta)cos(phi),rsin(theta)cos(phi),rcos(theta)) r=1,unit length
		}

		//after running through if-statements above we now have particle type definition, position, mom. direction, and the energy (or their initialised values)
		fParticleGun->SetParticlePosition(thisEffPosition);
		fParticleGun->SetParticleMomentumDirection(effdirection);
		fParticleGun->SetParticleEnergy(fEffEnergy);

		if(fDetector->Spice()) {
			fHistoManager->Fill2DHistogram(fHistoManager->AngleDistro(2), thisEffPosition.x(), thisEffPosition.y(), 1.0);
			fHistoManager->FillHistogram(fHistoManager->AngleDistro(3), thisEffPosition.z() - zZero);
			fHistoManager->FillHistogram(fHistoManager->AngleDistro(0), fEffEnergy);
			fHistoManager->BeamEnergy(fEffEnergy);
			fHistoManager->BeamTheta(acos(effdirection.z()/effdirection.mag()));
			fHistoManager->BeamPhi(atan2(effdirection.y(),effdirection.x()));
		}
	}


	// Set Optional Polarization
	if(fEffPolarization) {
		fParticleGun->SetParticlePolarization(fEffPolarizationVector);
	}

	//create vertex
	//
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::PassEfficiencyPosition(G4ThreeVector  num){
	fDetector->PassEfficiencyPosition(num);
}

void PrimaryGeneratorAction::PrepareBeamFile(G4String filename){
	fBeamDistribution = new BeamDistribution(filename);
	if(fBeamDistribution->Good())fNeedFileDistro = true;
}


void PrimaryGeneratorAction::SetLayeredTargetBeamDistro(G4int layer){
	fTargetDistro=true;
	fLayerStart=fDetector->LayeredTargetLayerStart(layer);
	fLayerLength=fDetector->LayeredTargetLayerStart(layer+1);
	fLayerLength-=fLayerStart;
}


void PrimaryGeneratorAction::LaBrinit() {
	//default LaBr properties
	G4double triangleThetaAngle = 54.735610317245360*deg;
	// theta
	fDetectorAnglesLaBr3[0][0] 	= triangleThetaAngle;
	fDetectorAnglesLaBr3[1][0] 	= triangleThetaAngle;
	fDetectorAnglesLaBr3[2][0] 	= triangleThetaAngle;
	fDetectorAnglesLaBr3[3][0] 	= triangleThetaAngle;
	fDetectorAnglesLaBr3[4][0] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAnglesLaBr3[5][0] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAnglesLaBr3[6][0] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAnglesLaBr3[7][0] 	= 180.0*deg - triangleThetaAngle;
	// phi
	fDetectorAnglesLaBr3[0][1] 	= 22.5*deg;
	fDetectorAnglesLaBr3[1][1] 	= 112.5*deg;
	fDetectorAnglesLaBr3[2][1] 	= 202.5*deg;
	fDetectorAnglesLaBr3[3][1] 	= 292.5*deg;
	fDetectorAnglesLaBr3[4][1] 	= 22.5*deg;
	fDetectorAnglesLaBr3[5][1] 	= 112.5*deg;
	fDetectorAnglesLaBr3[6][1] 	= 202.5*deg;
	fDetectorAnglesLaBr3[7][1] 	= 292.5*deg;
	// yaw (alpha)
	fDetectorAnglesLaBr3[0][2] 	= 0.0*deg;
	fDetectorAnglesLaBr3[1][2] 	= 0.0*deg;
	fDetectorAnglesLaBr3[2][2] 	= 0.0*deg;
	fDetectorAnglesLaBr3[3][2] 	= 0.0*deg;
	fDetectorAnglesLaBr3[4][2] 	= 0.0*deg;
	fDetectorAnglesLaBr3[5][2] 	= 0.0*deg;
	fDetectorAnglesLaBr3[6][2] 	= 0.0*deg;
	fDetectorAnglesLaBr3[7][2] 	= 0.0*deg;
	// pitch (beta)
	fDetectorAnglesLaBr3[0][3] 	= triangleThetaAngle;
	fDetectorAnglesLaBr3[1][3] 	= triangleThetaAngle;
	fDetectorAnglesLaBr3[2][3] 	= triangleThetaAngle;
	fDetectorAnglesLaBr3[3][3] 	= triangleThetaAngle;
	fDetectorAnglesLaBr3[4][3] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAnglesLaBr3[5][3] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAnglesLaBr3[6][3] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAnglesLaBr3[7][3] 	= 180.0*deg - triangleThetaAngle;
	// roll (gamma)
	fDetectorAnglesLaBr3[0][4] 	= 22.5*deg;
	fDetectorAnglesLaBr3[1][4] 	= 112.5*deg;
	fDetectorAnglesLaBr3[2][4] 	= 202.5*deg;
	fDetectorAnglesLaBr3[3][4] 	= 292.5*deg;
	fDetectorAnglesLaBr3[4][4] 	= 22.5*deg;
	fDetectorAnglesLaBr3[5][4] 	= 112.5*deg;
	fDetectorAnglesLaBr3[6][4] 	= 202.5*deg;
	fDetectorAnglesLaBr3[7][4] 	= 292.5*deg;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
