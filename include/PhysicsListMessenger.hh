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
// $Id: PhysicsListMessenger.hh 68007 2013-03-13 11:28:03Z gcosmo $
//
/// \file radioactivedecay/rdecay02/include/PhysicsListMessenger.hh
/// \brief Definition of the PhysicsListMessenger class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PHYSICSLISTMESSENGER_HH
#define PHYSICSLISTMESSENGER_HH

#include "globals.hh"
#include "G4UImessenger.hh"

class PhysicsList;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsListMessenger: public G4UImessenger
{
public:

    PhysicsListMessenger(PhysicsList* );
    virtual ~PhysicsListMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

private:

    PhysicsList* fPPhysicsList;

    G4UIdirectory*             fPhysDir;
    G4UIcmdWithADoubleAndUnit* fGammaCutCmd;
    G4UIcmdWithADoubleAndUnit* fElectCutCmd;
    G4UIcmdWithADoubleAndUnit* fProtoCutCmd;
    G4UIcmdWithADoubleAndUnit* fAllCutCmd;
    G4UIcmdWithADoubleAndUnit* fMCutCmd;
    G4UIcmdWithADoubleAndUnit* fECutCmd;
    G4UIcmdWithAString*        fPListCmd;
    G4UIcmdWithABool*          fConstructOpCmd;
    G4UIcmdWithABool*          fSpiceStepperCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
