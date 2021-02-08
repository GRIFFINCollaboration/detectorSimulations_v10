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
/// \file analysis/AnaEx01/AnaEx01.cc
/// \brief Main program of the analysis/AnaEx01 example
//
//
// $Id: AnaEx01.cc 73919 2013-09-17 07:38:47Z gcosmo $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <ctime>
#include <functional>
#include <string>
#include <unistd.h>

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#include "Randomize.hh"
#include "G4ParticleHPManager.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include "G4VisExecutive.hh"

#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	// Calculate seed from the product of PID, hostname hash, and time
	pid_t pid = getpid();

	char* hostname = new char[1024];
	gethostname(hostname, 1024);
	hostname[1023] = '\0';
	size_t hostnamehash = std::hash<std::string>{}(hostname);

	size_t time = std::time(nullptr);

	size_t seed = time*pid*hostnamehash;

	G4Random::setTheSeed(seed);

	// Construct the default run manager
#ifdef G4MULTITHREADED
	G4int nThreads = 2;
	if(argc == 3) {
		nThreads = strtol(argv[2], nullptr, 10);
	}
	G4cout<<"RUNNING MULTITHREADED WITH "<<nThreads<<" THREADS"<<G4endl;
	G4MTRunManager* runManager = new G4MTRunManager;
	runManager->SetNumberOfThreads(nThreads);
#else
	G4cout<<"NOT RUNNING MULTITHREADED"<<G4endl;
	G4RunManager* runManager = new G4RunManager;
#endif

	// turn off messages from particle HP manager (/process/had/particle_hp/verbose command does not work?)
	G4ParticleHPManager::GetInstance()->SetVerboseLevel(0);

	// Set mandatory initialization classes
	DetectorConstruction* detector = new DetectorConstruction;
	runManager->SetUserInitialization(detector);
	runManager->SetUserInitialization(new PhysicsList);
	runManager->SetUserInitialization(new ActionInitialization(detector));

	// We don't initialize the G4 kernel at run time so the physics list can be changed!

	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();

	if(argc != 1) { // batch mode
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);
	} else { // interactive mode : define visualization and UI terminal
		G4UIExecutive* ui = new G4UIExecutive(argc, argv);
		UImanager->ApplyCommand("/control/execute vis.mac");

		ui->SessionStart();

		delete ui;
	}

	delete visManager;

	// Job termination
	delete runManager;

	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
