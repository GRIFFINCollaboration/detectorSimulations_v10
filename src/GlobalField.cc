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
//
//

#include <time.h>

#include "Randomize.hh"
#include "G4TransportationManager.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"

#include "GlobalField.hh"

using namespace CLHEP;

GlobalField* GlobalField::object = 0;

GlobalField::GlobalField() : G4ElectroMagneticField(),
                             minStep(0.01*mm), deltaChord(3.0*mm),
                             deltaOneStep(0.01*mm), deltaIntersection(0.1*mm),
                             epsMin(2.5e-7*mm), epsMax(0.05*mm),
                             fEquation(0), fFieldManager(0), 
                             fFieldPropagator(0), fStepper(0), fChordFinder(0)
//GlobalField::GlobalField() : G4MagneticField(),
//                             minStep(0.01*mm), deltaChord(3.0*mm),
//                             deltaOneStep(0.01*mm), deltaIntersection(0.1*mm),
//                             epsMin(2.5e-7*mm), epsMax(0.05*mm),
//                             fEquation(0), fFieldManager(0),
//                             fFieldPropagator(0), fStepper(0), fChordFinder(0)
{
  fFieldMessenger = new FieldMessenger(this);

  fields = new FieldList();

  fStepperType = 4 ;       // ClassicalRK4 is default stepper

  //  set object

  object = this;

  updateField();
}

GlobalField::~GlobalField()
{
  clear();

  delete fFieldMessenger;

  if (fEquation)        delete fEquation;
  if (fFieldManager)    delete fFieldManager;
  if (fFieldPropagator) delete fFieldPropagator;
  if (fStepper)         delete fStepper;
  if (fChordFinder)     delete fChordFinder;
}

void GlobalField::updateField()
{
  first = true;

  nfp = 0;
  fp = 0;

  clear();

  //  Construct equ. of motion of particles through B fields
//  fEquation = new G4Mag_EqRhs(this);
  //  Construct equ. of motion of particles through e.m. fields
//  fEquation = new G4EqMagElectricField(this);
  //  Construct equ. of motion of particles including spin through B fields
//  fEquation = new G4Mag_SpinEqRhs(this);
  //  Construct equ. of motion of particles including spin through e.m. fields
  fEquation = new G4EqEMFieldWithSpin(this);

  //  Get transportation, field, and propagator managers
  G4TransportationManager* fTransportManager =
         G4TransportationManager::GetTransportationManager();

  fFieldManager = GetGlobalFieldManager();

  fFieldPropagator = fTransportManager->GetPropagatorInField();

  //  Need to SetFieldChangesEnergy to account for a time varying electric
  //  field (r.f. fields)
  fFieldManager->SetFieldChangesEnergy(true);

  //  Set the field
  fFieldManager->SetDetectorField(this);

  //  Choose a stepper for integration of the equation of motion
  SetStepper();

  //  Create a cord finder providing the (global field, min step length,
  //  a pointer to the stepper)
  fChordFinder = new G4ChordFinder((G4MagneticField*)this,minStep,fStepper);

  // Set accuracy parameters
  fChordFinder->SetDeltaChord( deltaChord );

  fFieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

  fFieldManager->SetDeltaIntersection(deltaIntersection);

  fFieldPropagator->SetMinimumEpsilonStep(epsMin);
  fFieldPropagator->SetMaximumEpsilonStep(epsMax);

  G4cout << "Accuracy Parameters:" <<
            " MinStep=" << minStep <<
            " DeltaChord=" << deltaChord <<
            " DeltaOneStep=" << deltaOneStep << G4endl;
  G4cout << "                    " <<
            " DeltaIntersection=" << deltaIntersection <<
            " EpsMin=" << epsMin <<
            " EpsMax=" << epsMax <<  G4endl;

  fFieldManager->SetChordFinder(fChordFinder);

}

GlobalField* GlobalField::getObject()
{
  if (!object) new GlobalField();
  return object;
}

void GlobalField::SetStepper()
{
  if(fStepper) delete fStepper;

  switch ( fStepperType )
  {
    case 0:
//      fStepper = new G4ExplicitEuler( fEquation, 8 ); // no spin tracking
      fStepper = new G4ExplicitEuler( fEquation, 12 ); // with spin tracking
      G4cout << "G4ExplicitEuler is called" << G4endl;
      break;
    case 1:
//      fStepper = new G4ImplicitEuler( fEquation, 8 ); // no spin tracking
      fStepper = new G4ImplicitEuler( fEquation, 12 ); // with spin tracking
      G4cout << "G4ImplicitEuler is called" << G4endl;
      break;
    case 2:
//      fStepper = new G4SimpleRunge( fEquation, 8 ); // no spin tracking
      fStepper = new G4SimpleRunge( fEquation, 12 ); // with spin tracking
      G4cout << "G4SimpleRunge is called" << G4endl;
      break;
    case 3:
//      fStepper = new G4SimpleHeum( fEquation, 8 ); // no spin tracking
      fStepper = new G4SimpleHeum( fEquation, 12 ); // with spin tracking
      G4cout << "G4SimpleHeum is called" << G4endl;
      break;
    case 4:
//      fStepper = new G4ClassicalRK4( fEquation, 8 ); // no spin tracking
      fStepper = new G4ClassicalRK4( fEquation, 12 ); // with spin tracking
      G4cout << "G4ClassicalRK4 (default) is called" << G4endl;
      break;
    case 5:
//      fStepper = new G4CashKarpRKF45( fEquation, 8 ); // no spin tracking
      fStepper = new G4CashKarpRKF45( fEquation, 12 ); // with spin tracking
      G4cout << "G4CashKarpRKF45 is called" << G4endl;
      break;
    default: fStepper = 0;
  }
}

G4FieldManager* GlobalField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
                                ->GetFieldManager();
}

void GlobalField::GetFieldValue(const G4double* point, G4double* field) const
{
  //G4cout << "GetFieldValue : xpos = " << point[0]/mm << " mm ypos = " << point[1]/mm << " mm zpos = " << point[2]/mm << " mm" << G4endl;

  // NOTE: this routine dominates the CPU time for tracking.
  // Using the simple array fp[] instead of fields[] 
  // directly sped it up

	// TB: Not sure I understand the point of all these if statements. In every case all of the field[] 
	// arrays will be set to zero. 

  G4double x_length = 0.15*m;
  G4double y_length = 0.3*m;
  G4double z_length = 0.01*m;

  if((point[0] < (x_length/2.0)) && (point[0] > (-1.0*x_length/2.0))) 
  {
  	if((point[1] < (y_length/2.0)) && (point[1] > (-1.0*y_length/2.0))) 
  	{
 			if((point[2] < (z_length/2.0)) && (point[2] > (-1.0*z_length/2.0))) 
			{
      	field[0] = field[1] = field[2] = field[3] = field[4] = field[5] = 0.0;
			}
  		else 
			{
	  	  field[0] = field[1] = field[2] = field[3] = field[4] = field[5] = 0.0;
			}
  	}
  	else 
  	{
    	field[0] = field[1] = field[2] = field[3] = field[4] = field[5] = 0.0;
  	}
  }
  else 
  {
    field[0] = field[1] = field[2] = field[3] = field[4] = field[5] = 0.0;
  }
//  field[0] = field[1] = field[2] = field[3] = field[4] = field[5] = 0.0;

  // protect against Geant4 bug that calls us with point[] NaN.
  if(point[0] != point[0]) return;

  // (can't use nfp or fp, as they may change)
  if (first) ((GlobalField*)this)->setupArray();   // (cast away const)

  for (int i=0; i<nfp; ++i) {
      const ElementField* p = fp[i];
      if (p->isInBoundingBox(point)) {
         p->addFieldValue(point,field);
      }
  }

//  G4double rad_unitless = point[0]*point[1]/(point[2]*point[2]+1.0);

//  field[0] = 0.0*tesla;
//  field[1] = 0.0*tesla;
//  field[2] = 0.025*tesla*rad_unitless;

}

void GlobalField::clear()
{
  if (fields) {
     if (fields->size()>0) {
        FieldList::iterator i;
        for (i=fields->begin(); i!=fields->end(); ++i) delete *i;
        fields->clear();
     }
  }

  if (fp) delete[] fp;

  first = true;

  nfp = 0;
  fp = NULL;
}

void GlobalField::setupArray()
{
  first = false;
  nfp = fields->size();
  fp = new const ElementField* [nfp+1]; // add 1 so it's never 0
  for (int i=0; i<nfp; ++i) fp[i] = (*fields)[i];
}
