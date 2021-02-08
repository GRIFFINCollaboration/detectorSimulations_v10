#include "HistoMessenger.hh"
#include "HistoManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithADouble.hh"

#include <string>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoMessenger::HistoMessenger(HistoManager* hist)
:fHist(hist)
{
	fDir = new G4UIdirectory("/Histo/");
	fDir->SetGuidance("UI commands for histo manager");
    
	fRecordGunCmd = new G4UIcmdWithABool("/Histo/RecordGun",this);
	fRecordGunCmd->SetGuidance("Record the particle for each event in the tree");
	fRecordGunCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
    
	fRecordAllCmd = new G4UIcmdWithABool("/Histo/RecordAll",this);
	fRecordAllCmd->SetGuidance("Record each generated event in the tree (and not only those hitting detectors)");
	fRecordAllCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

	fFileNameCmd = new G4UIcmdWithAString("/Histo/FileName", this);
	fFileNameCmd->SetGuidance("Set output file name, default is g4out. Warning, .root or _tX.root will automatically be added!");
	fFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

HistoMessenger::~HistoMessenger()
{
	delete fDir;
	delete fRecordGunCmd;
	delete fRecordAllCmd;
	delete fFileNameCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fRecordGunCmd) {
		fHist->RecordGun(fRecordGunCmd->GetNewBoolValue(newValue));
		return;
	}
	if(command == fRecordAllCmd) {
		fHist->RecordAll(fRecordAllCmd->GetNewBoolValue(newValue));
		return;
	}
	if(command == fFileNameCmd) {
		fHist->FileName(newValue);
		return;
	}
}

