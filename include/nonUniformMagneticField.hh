#ifndef nonUniformMagneticField_H
#define nonUniformMagneticField_H 1

#include "G4MagneticField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
//class F03FieldMessenger;

class nonUniformMagneticField
{

public:

  nonUniformMagneticField(const char* fieldName, double zOffset, double zRotation); // field read from file
  ~nonUniformMagneticField() ;
      
  void SetStepperType( G4int i) { fStepperType = i ; }

  void SetStepper();

  void SetMinStep(G4double s) { fMinStep = s ; }

  void UpdateField();

protected:

  G4FieldManager*         GetGlobalFieldManager() ;
    // Returns the global Field Manager

  G4FieldManager*         fFieldManager ;
 
  G4ChordFinder*          fChordFinder ;
 
  G4Mag_UsualEqRhs*       fEquation ; 
 
  G4MagneticField*        fMagneticField ; 
 
  G4MagIntegratorStepper* fStepper ;
  G4int                   fStepperType ;

  G4double                fMinStep ;
 
//  F03FieldMessenger*      fFieldMessenger;

};

#endif
