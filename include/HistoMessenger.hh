//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HISTOMESSENGER_HH
#define HISTOMESSENGER_HH

#include "globals.hh"
#include "G4UImessenger.hh"

class HistoManager;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoMessenger: public G4UImessenger
{
public:
	HistoMessenger(HistoManager*);
	~HistoMessenger();

	void SetNewValue(G4UIcommand*, G4String);

private:
	HistoManager*               fHist;

	G4UIdirectory*              fDir;
	
	G4UIcmdWithABool*           fRecordGunCmd;
	G4UIcmdWithABool*           fRecordAllCmd;
	G4UIcmdWithAString*			 fFileNameCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

