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
#include "DetectorConstruction.hh" //for detector based information

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Geantino.hh"
#include "Randomize.hh"

#include "HistoManager.hh"
using namespace CLHEP; //for pi

//#include <stdlib.h> //for abs

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(0),
      fDetector(DC)

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
    fEffBeam = false;
    LaBrinit(); //sets up default variables - messy having them all declared here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
    delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
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
        for( G4int j = 0 ; j < fNumberOfDecayingLaBrDetectors ; j++ ) {  // get the number of particles in decay and loop over them
            sumProb = sumProb + prob;
            if(randomDet <= sumProb ) {
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

        if(randomEnergy <= prob ) {
            thisEnergy = 788.742*keV;
        }
        else {
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
        G4double randsintheta = sqrt( 1. - randcostheta*randcostheta );
        G4double randphi      = (360.*deg)*G4UniformRand();
        G4ThreeVector thisDirection = G4ThreeVector(randsintheta*cos(randphi), randsintheta*sin(randphi), randcostheta);

        fParticleGun->SetParticleDefinition(agamma);
        fParticleGun->SetParticlePosition(thisPosition);
        fParticleGun->SetParticleMomentumDirection(thisDirection);
        fParticleGun->SetParticleEnergy(thisEnergy);
    }
	
    else if(fEffEnergy != 0.0) {
        G4ParticleDefinition* effPart;
        if(fEffParticleBool) {

            if(fEffParticle == "electron" || fEffParticle == "e-") {
                effPart = G4ParticleTable::GetParticleTable()->FindParticle("e-");
            }
            else if(fEffParticle == "positron" || fEffParticle == "e+") {
                effPart = G4ParticleTable::GetParticleTable()->FindParticle("e+");
            }
            else if (fEffParticle == "gamma" || fEffParticle == "photon"){
                effPart = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
            }
            else if (fEffParticle == "neutron"){
                effPart = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
            }
        }
        else {

            effPart = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        }
        
        G4ThreeVector thisEffPosition = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);//in constructor?? //19/7
        if(fEffPositionBool) thisEffPosition = fEffPosition;

	//DIRECTION CONTROLS (directly forced, beam, or cone)
        G4double effRandCosTheta, effRandSinTheta, effRandPhi;
        G4ThreeVector effdirection;
        if(fEffDirectionBool) { 
            effdirection = fEffDirection;
            // If we want to simulate a realistic beam spot, instead of perfect pencil beam.
            if(fEffBeam) {// following SetEfficiencyBeamRadius command
	      G4cout << "Beam " << G4endl;
                G4double xMonte = 10000.0*m; //monte for monte-carlo, effective errors to a central beam, creates a distribution
                G4double yMonte = 10000.0*m;// why these values?
                G4double zMonte = 0.0*m;//SPICE should have z/y errors?? as x along beam axis
                G4ThreeVector vecMonte;

                G4double directionTheta = fEffDirection.theta();
                G4double directionPhi = fEffDirection.phi();

                while( pow(xMonte,2) + pow(yMonte,2) > pow(fEffBeamRadius,2) ) {//xmonte^2 + ymonte^2 > BeamRadius^2
                    xMonte = (2.*G4UniformRand()-1.0)*fEffBeamRadius;
                    yMonte = (2.*G4UniformRand()-1.0)*fEffBeamRadius;
                }

                vecMonte = G4ThreeVector(xMonte,yMonte,zMonte); //rotate# funcs are CLHEP
                vecMonte.rotateY(directionTheta);
                vecMonte.rotateZ(directionPhi);

                thisEffPosition = thisEffPosition + vecMonte;
            }

	    if(fConeRadiusBool){// following SetConeRadius command - emits here in 2pi (Not conal as no x-y limits placed
	   // G4cout << "Conal " << G4endl;
	      effdirection = SetCone(fConeRadius);//unit direction in -Z due to SPICE's downstream detector
	  //  G4cout << "Direction " << effdirection << G4endl << " Cone radius " << fConeRadius << "mm"<< G4endl;
	    }
	    if(fConeValueBool){
	      effdirection = SetCone(fConeRValue, fConeZValue);
	    }
	    if(fConeAngleBool){
	      
	      //808's and magic //limit theta to input, phi unrestricted
	      G4double theta = G4UniformRand()*fAngleInit;
	      G4double CosTheta = cos(theta);
	      G4double SinTheta = sqrt(1. - pow(CosTheta, 2.0));
	      G4double Phi      = (2.0*pi)*G4UniformRand();
	      
	      effdirection = G4ThreeVector(SinTheta*cos(Phi), SinTheta*sin(Phi), -CosTheta);
	      
	   if (fConeAngleBool) HistoManager::Instance().FillHisto(HistoManager::Instance().angledistro[1], theta);
  if (fConeAngleBool) HistoManager::Instance().FillHisto(HistoManager::Instance().angledistro[2], sin(theta));      
  if (fConeAngleBool) HistoManager::Instance().FillHisto(HistoManager::Instance().angledistro[0], SinTheta);//tan(phi) if yCone/xCone 	      
	  if (fConeAngleBool) HistoManager::Instance().FillHisto(HistoManager::Instance().angledistro[3], 2*pi*(1.-CosTheta));//tan(phi) if yCone/xCone 	      
	     // effdirection = SetCone(tan(fAngleInit),1);//read in as degree, GEANT auto-convert to rads
	    }
	} else {
	    //G4cout << "Random " << G4endl; //may offer the solution, an altered 2pi rando. Using 4pi for efficiency
            // random direction if no preference provided
            effRandCosTheta = 2.*G4UniformRand()-1.0; //cos(theta) = 2cos^2(0.5theta)-1 ??
            effRandSinTheta = sqrt( 1. - effRandCosTheta*effRandCosTheta ); //from sin^2(theta)+cos^2(theta)=1
            effRandPhi      = (360.*deg)*G4UniformRand();
            effdirection = G4ThreeVector(effRandSinTheta*cos(effRandPhi), effRandSinTheta*sin(effRandPhi), effRandCosTheta);
	    //converts from Spherical polar(physics def.) to cartesian via (rsin(theta)cos(phi),rsin(theta)cos(phi),rcos(theta)) r=1,unit length
        }
	//std::cout<<thisEffPosition<<effPart<<effdirection<<fEffEnergy<<std::endl; //gives beam starting pos as expected
	//effpart is mem location, but all the same indicating just the electron //effdirection is random as expected
		
	//after running through if-statements above we now have particle type definition, position, mom. direction, and the energy (or their initialised values)
        fParticleGun->SetParticleDefinition(effPart);
        fParticleGun->SetParticlePosition(thisEffPosition);
        fParticleGun->SetParticleMomentumDirection(effdirection);
        fParticleGun->SetParticleEnergy(fEffEnergy);
    }


    // Set Optional Polarization
    if(fEffPolarization) {
        fParticleGun->SetParticlePolarization(fEffPolarizationVector);
    } 

    //create vertex
    //
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
G4ThreeVector PrimaryGeneratorAction::SetCone(G4double ConeRadius, G4double zVal) {
  G4double xCone=0.0, yCone=0.0;//107.5mm from particle creation to SiLi, 47mm SiLi radius ->23 deg cone, 23/90 = 0.26237 for angleratio
  G4ThreeVector coneDirection;
  
  yCone = G4UniformRand()*(ConeRadius);//Random x up to ConeRadius
  xCone = G4UniformRand()*(sqrt(pow(ConeRadius,2)-pow(yCone,2)));//Random y from whatever is left from the maximum of r
 
  xCone=atan(xCone/zVal);//a length at this point, need an angle ratio
  yCone=atan(yCone/zVal);
  
  G4double xSign = G4UniformRand(), ySign= G4UniformRand(); 
  if (xSign>0.5) xCone=-xCone;//50/50 chance of + or - (sign set here) direction for cone
  if (ySign>0.5) yCone=-yCone;

  //G4cout << xCone << "\t" << yCone << "\t" << atan(yCone/xCone) << G4endl;
  coneDirection = G4ThreeVector(xCone,yCone,-1.0);//(SinTheta*cos(Phi), SinTheta*sin(Phi), CosTheta)
  //sin(theta)=(SinTheta*cos(Phi)^2 + SinTheta*sin(Phi)^2)^1/2
  G4double sinthetaoutput=sqrt(pow(xCone,2.0)+pow(yCone,2.0));

  return coneDirection;
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
