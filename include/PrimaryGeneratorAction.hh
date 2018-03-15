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
/// \file analysis/shared/include/PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.hh 68015 2013-03-13 13:27:27Z gcosmo $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"


class G4ParticleGun;
class G4Event;

class PrimaryGeneratorMessenger;
class DetectorConstruction;
class HistoManager;
class BeamDistribution;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
private:    
    G4ParticleGun*                fParticleGun;  //pointer a to G4 class
    DetectorConstruction*         fDetector;     //pointer to the geometry
    PrimaryGeneratorMessenger*    fGunMessenger; //messenger of this class
    BeamDistribution*		  fBeamDistribution;//pointer to the BeamDistribution class

public:
    PrimaryGeneratorAction(DetectorConstruction*);
    virtual ~PrimaryGeneratorAction();
    virtual void GeneratePrimaries(G4Event*);
    

    void SetNumberOfDecayingLaBrDetectors( G4int num ) {fNumberOfDecayingLaBrDetectors = num;} ;
    void SetEfficiencyEnergy( G4double num ) {fEffEnergy = num;} ;
    void SetEfficiencyDirection( G4ThreeVector num ) {fEffDirection = num; fEffDirectionBool = true;} ;
    
    void PassEfficiencyPosition( G4ThreeVector num );
    void SetEfficiencyPosition( G4ThreeVector num ) {fEffPosition = num; fEffPositionBool = true;PassEfficiencyPosition(fEffPosition);}
    
    void SetEfficiencyParticle( G4String val ) {fEffParticle = val; fEffParticleBool = true;} ;
    void SetEfficiencyPolarization( G4ThreeVector num ) {fEffPolarizationVector = num; fEffPolarization = true;} ;
    void SetConeMaxAngle( G4double num1 ) {fAngleInit = num1; fConeAngleBool = true; fEffDirectionBool = true;};
    void SetConeMinAngle( G4double num1 ) {fAngleMinInit = num1;};
    void SetBeamSpotSigma( G4double num1 ) {fBeamSpotSigma = num1;};
    void PassTarget(G4double);
    void PrepareBeamFile(G4String);
    void SetLayeredTargetBeamDistro(G4int layer);
 
    
    void SetSourceNeeded(G4bool needed){fSourceNeeded = needed;};
    void SetSourceName(G4String input){fSourceName = input;};
    
private:
    //variables
    G4int fNumberOfDecayingLaBrDetectors;
    G4double fEffEnergy;
    G4ThreeVector fEffDirection;
    G4ThreeVector fEffPosition;
    G4bool fEffDirectionBool;
    G4bool fEffPositionBool;
    G4String fEffParticle;
    G4bool fEffParticleBool;
    G4double fDetectorAnglesLaBr3[8][5];
    G4bool fEffPolarization;
    G4ThreeVector fEffPolarizationVector;
    
    G4double fAngleInit;
    G4bool fConeAngleBool;
    G4double fAngleMinInit;
    G4double fBeamSpotSigma;

    G4bool fTargetDistro;
    G4double fLayerStart;
    G4double fLayerLength;
    G4bool fNeedFileDistro;
    
    G4bool fSourceNeeded;
    G4String fSourceName;
    
    //functions
    void LaBrinit();
        
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


