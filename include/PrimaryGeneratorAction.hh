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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(DetectorConstruction*);
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    void SetNumberOfDecayingLaBrDetectors( G4int num ) {numberOfDecayingLaBrDetectors = num;} ;
    void SetEfficiencyEnergy( G4double num ) {effEnergy = num;} ;
    void SetEfficiencyDirection( G4ThreeVector num ) {effDirection = num; effDirectionBool = true;} ;
    void SetEfficiencyPosition( G4ThreeVector num ) {effPosition = num; effPositionBool = true;} ;
    void SetEfficiencyParticle( G4String val ) {effParticle = val; effParticleBool = true;} ;


private:
    G4ParticleGun*                fParticleGun;  //pointer a to G4 class
    DetectorConstruction*         fDetector;     //pointer to the geometry
    PrimaryGeneratorMessenger*    fGunMessenger; //messenger of this class

    G4int numberOfDecayingLaBrDetectors;
    G4double effEnergy;
    G4ThreeVector effDirection;
    G4ThreeVector effPosition;
    G4bool effDirectionBool;
    G4bool effPositionBool;
    G4String effParticle;
    G4bool effParticleBool;
    G4double detectorAnglesLaBr3[8][5];

    G4double transX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


