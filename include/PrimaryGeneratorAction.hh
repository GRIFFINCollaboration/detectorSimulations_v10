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
class DetectorConstruction;
class PrimaryGeneratorMessenger;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(DetectorConstruction*);
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    void SetNumberOfDecayingLaBrDetectors( G4int num ) {fNumberOfDecayingLaBrDetectors = num;} ;
    void SetEfficiencyEnergy( G4double num ) {fEffEnergy = num;} ;
    void SetEfficiencyDirection( G4ThreeVector num ) {fEffDirection = num; fEffDirectionBool = true;} ;
    void SetEfficiencyPosition( G4ThreeVector num ) {fEffPosition = num; fEffPositionBool = true;} ;
    void SetEfficiencyParticle( G4String val ) {fEffParticle = val; fEffParticleBool = true;} ;
    void SetEfficiencyPolarization( G4ThreeVector num ) {fEffPolarizationVector = num; fEffPolarization = true;} ;
    void SetEfficiencyBeamRadius( G4double num ) {fEffBeamRadius = num; fEffBeam = true; } ;
    void SetConeRadius( G4double num ) {fConeRadius = num; fConeRadiusBool = true; fEffDirectionBool = true;} ;//Direction needed, should not require explicit initialisation command
    void SetConeZValue( G4double num ) {fConeZValue = num; fConeValueBool = true; fEffDirectionBool = true;};
    void SetConeRValue( G4double num ) {fConeRValue = num; fConeValueBool = true; fEffDirectionBool = true;};
    void SetConeMaxAngle( G4double num1 ) {fAngleInit = num1; fConeAngleBool = true; fEffDirectionBool = true;};
    void SetConeMinAngle( G4double num1 ) {fAngleMinInit = num1;};
    //booleans (initially false), above, true if a command has been entered for the loops in source file to be entered
    void SendBeamEnergyToHist(G4double);
    void PassTarget(G4double);
    G4bool fNeedBeamDistro;
    

private:
    G4ParticleGun*                fParticleGun;  //pointer a to G4 class
    DetectorConstruction*         fDetector;     //pointer to the geometry
    PrimaryGeneratorMessenger*    fGunMessenger; //messenger of this class

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
    G4bool fEffBeam;
    G4double fEffBeamRadius;
    G4double fConeRadius;
    G4bool fConeRadiusBool;
    G4double fConeZValue;
    G4double fConeRValue;
    G4bool fConeValueBool;
    G4double fAngleInit;
    G4bool fConeAngleBool;
    G4double fAngleMinInit;

    G4ThreeVector SetCone(G4double fConeRadius, G4double zVal=107.50685);//value is disttoSiliFromParticleCreation
    void LaBrInit();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


