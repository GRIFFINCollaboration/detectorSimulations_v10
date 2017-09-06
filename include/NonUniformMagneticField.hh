#ifndef NONUNIFORMMAGNETICFIELD_HH
#define NONUNIFORMMAGNETICFIELD_HH

#include "G4MagneticField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
//class F03FieldMessenger;

class NonUniformMagneticField
{

public:

  NonUniformMagneticField(const char* fieldName, double zOffset, double zRotation); // field read from file
  ~NonUniformMagneticField() ;
      
  void SetStepperType(G4int i) { fStepperType = i; }
  void SetStepper();
  void SetMinStep(G4double step) { fMinStep = step; }
  void UpdateField();

protected:

  G4FieldManager*         GetGlobalFieldManager();    // Returns the global Field Manager

  G4FieldManager*         fFieldManager;
  G4ChordFinder*          fChordFinder;
  G4Mag_UsualEqRhs*       fEquation; 
  G4MagneticField*        fMagneticField; 
  G4MagIntegratorStepper* fStepper;
  G4int                   fStepperType;
  G4double                fMinStep;
};

#endif
