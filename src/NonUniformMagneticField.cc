#include "NonUniformMagneticField.hh"

#include "TabulatedMagneticField.hh" //called in constructor

#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

NonUniformMagneticField::NonUniformMagneticField(const char* fieldName="./", G4double zOffset=0., G4double zRotation=0.) //hard code here? as no changes made??
: fChordFinder(0), fStepper(0)
{
	fMagneticField = new TabulatedMagneticField(fieldName, zOffset, zRotation);//sets field to that read in from the table
	GetGlobalFieldManager()->CreateChordFinder(fMagneticField);

	fEquation = new G4Mag_UsualEqRhs(fMagneticField); 

	fMinStep     = 0.001*CLHEP::mm; // minimal step of 1 mm is default //0.25 potentially but seems too high //chosen is neither
	fStepperType = 4;      // ClassicalRK4 is default stepper

	fFieldManager = GetGlobalFieldManager();

	UpdateField();
}


NonUniformMagneticField::~NonUniformMagneticField()
{
	delete fMagneticField;
	delete fChordFinder;
	delete fStepper;
}

void NonUniformMagneticField::UpdateField()
{
	SetStepper();
	G4cout<<"The minimal step is equal to "<<fMinStep/CLHEP::mm<<" mm"<<G4endl;

	fFieldManager->SetDetectorField(fMagneticField);

	delete fChordFinder; //why delete here and remake below? - allows for new steppers/minimum steps/or field, but this is never invoked may allow for greater optimisation

	fChordFinder = new G4ChordFinder(fMagneticField, fMinStep, fStepper);

	fFieldManager->SetChordFinder(fChordFinder);
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void NonUniformMagneticField::SetStepper()
{
	delete fStepper;

	switch(fStepperType) 
	{
		case 0:  
			fStepper = new G4ExplicitEuler(fEquation); 
			G4cout<<"G4ExplicitEuler is called"<<G4endl;     
			break;
		case 1:  
			fStepper = new G4ImplicitEuler(fEquation);      
			G4cout<<"G4ImplicitEuler is called"<<G4endl;     
			break;
		case 2:  
			fStepper = new G4SimpleRunge(fEquation);        
			G4cout<<"G4SimpleRunge is called"<<G4endl;     
			break;
		case 3:  
			fStepper = new G4SimpleHeum(fEquation);         
			G4cout<<"G4SimpleHeum is called"<<G4endl;     
			break;
		case 4:  
			fStepper = new G4ClassicalRK4(fEquation);       
			G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;     
			break;
		case 5:  
			fStepper = new G4HelixExplicitEuler(fEquation); 
			G4cout<<"G4HelixExplicitEuler is called"<<G4endl;     
			break;
		case 6:  
			fStepper = new G4HelixImplicitEuler(fEquation); 
			G4cout<<"G4HelixImplicitEuler is called"<<G4endl;     
			break;
		case 7:  
			fStepper = new G4HelixSimpleRunge(fEquation);   
			G4cout<<"G4HelixSimpleRunge is called"<<G4endl;     
			break;
		case 8:  
			fStepper = new G4CashKarpRKF45(fEquation);      
			G4cout<<"G4CashKarpRKF45 is called"<<G4endl;     
			break;
		case 9:  
			fStepper = new G4RKG3_Stepper(fEquation);       
			G4cout<<"G4RKG3_Stepper is called"<<G4endl;     
			break;
		default: 
			fStepper = nullptr;
	}
}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  NonUniformMagneticField::GetGlobalFieldManager()
{
	return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}
