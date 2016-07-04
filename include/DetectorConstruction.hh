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
/// \file analysis/shared/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.hh 77256 2013-11-22 10:10:23Z gcosmo $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
//class DetectionSystemGammaTracking;
class DetectionSystemGriffin;
class DetectionSystem8pi;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DetectionSystemDescant;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ApparatusDescantStructure;

class DetectionSystemDescant;

class DetectionSystemSceptar;
class DetectionSystemSpice;
class DetectionSystemSpiceV02;
class DetectionSystemPaces;
class DetectionSystemSodiumIodide;
class DetectionSystemLanthanumBromide;
class DetectionSystemBox;
class DetectionSystemAncillaryBGO;

//class MagneticField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction();

    G4int fGriffinDetectorsMapIndex;
    G4int fGriffinDetectorsMap[16];

    void SetWorldMaterial( G4String );
    void SetWorldDimensions( G4ThreeVector );
    void SetWorldVis( G4bool );
    //    void SetWorldMagneticField( G4ThreeVector );

    //    void SetGenericTargetMaterial( G4String );
    //    void SetGenericTargetDimensions( G4ThreeVector );
    //    void SetGenericTargetPosition( G4ThreeVector );
    //    void SetGenericTarget( );
    //    void SetFieldBoxMaterial( G4String );
    //    void SetFieldBoxDimensions( G4ThreeVector );
    //    void SetFieldBoxPosition( G4ThreeVector );
    //    void SetFieldBoxMagneticField( G4ThreeVector );
    //    void SetFieldBox( );

    //    void SetBoxMat( G4String input )                   {boxMat = input;};
    //        void SetBoxThickness( G4double input )             {boxThickness = input;};
    //        void SetBoxInnerDimensions( G4ThreeVector input )  {boxInnerDimensions = input;};
    //        void SetBoxColour( G4ThreeVector input )           {boxColour = input;};
    //        void AddBox();

    // Grid Functions
    void SetGridMat( G4String input )                  {fGridMat = input;};
    void SetGridSize( G4double input )                 {fGridSize = input;};
    void SetGridDimensions( G4ThreeVector input )      {fGridDimensions = input;};
    void SetGridColour( G4ThreeVector input )          {fGridColour = input;};
    void SetGridPosOffset( G4ThreeVector input )          {fGridOffset = input;};
    void AddGrid();

    //    void AddApparatusSpiceTargetChamber();
    void AddApparatus8piVacuumChamber();
    void AddApparatus8piVacuumChamberAuxMatShell(G4double thickness);
    void AddApparatusGriffinStructure(G4int selector);

    G4double GetWorldSizeX()           {return fWorldSizeX;};
    G4double GetWorldSizeY()           {return fWorldSizeY;};
    G4double GetWorldSizeZ()           {return fWorldSizeZ;};

    const G4VPhysicalVolume* GetphysiWorld() {return fPhysiWorld;};

    G4VPhysicalVolume* Construct();

    void UpdateGeometry();

    //    void AddDetectionSystemGammaTracking(G4int ndet);
    void AddDetectionSystemSodiumIodide(G4int ndet);
    void AddDetectionSystemLanthanumBromide(G4ThreeVector input);
    void AddDetectionSystemAncillaryBGO(G4ThreeVector input);

    void AddDetectionSystem8pi(G4int ndet);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void AddDetectionSystemDescant(G4int ndet);
    void AddDetectionSystemDescantAuxPorts(G4ThreeVector input);
 
    void SetDetectionSystemDescantRotation(G4ThreeVector input);
    void SetDetectionSystemDescantColor(G4String input);
    void AddDetectionSystemDescantCart(G4ThreeVector input);
    void AddDetectionSystemDescantSpher(G4ThreeVector input, G4double unit);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void AddApparatusDescantStructure();

    void AddDetectionSystemTestcan(G4ThreeVector input);

    void AddDetectionSystem8piDetector(G4int ndet);
    void AddDetectionSystemGriffinForward(G4int ndet);
    void AddDetectionSystemGriffinForwardDetector(G4int ndet);
    void AddDetectionSystemGriffinBack(G4int ndet);
    void AddDetectionSystemGriffinBackDetector(G4int ndet);
    //void AddDetectionSystemGriffinPositionConfig(G4ThreeVector input);
    void AddDetectionSystemGriffinHevimet( G4int input ) ;
    void AddDetectionSystemGriffinCustom( G4int ndet ) ;
    void AddDetectionSystemGriffinCustomDetector( G4int ndet ) ;
    void AddDetectionSystemGriffinShieldSelect( G4int ShieldSelect ) ;
    void AddDetectionSystemGriffinSetRadialDistance( G4double detectorDist ) ;
    void AddDetectionSystemGriffinSetExtensionSuppLocation( G4int detectorPos ) ;
    void AddDetectionSystemGriffinSetDeadLayer( G4ThreeVector params ) ;

    void AddDetectionSystemSceptar(G4int ndet);
    void AddDetectionSystemPaces(G4int ndet);
    void AddDetectionSystemSpice(G4int ndet);
    void AddDetectionSystemSpiceV02(G4int ndet);

    G4double GetLanthanumBromideCrystalRadius();
    G4double GetLanthanumBromideCrystalLength();
    G4double GetLanthanumBromideR();
    G4double GetLanthanumBromideTheta(G4int i);
    G4double GetLanthanumBromidePhi(G4int i);
    G4double GetLanthanumBromideYaw(G4int i);
    G4double GetLanthanumBromidePitch(G4int i);
    G4double GetLanthanumBromideRoll(G4int i);
    G4double GetLanthanumBromideCrystalRadialPosition();


    void UseTIGRESSPositions( G4bool input )                  {fUseTigressPositions = input;};
private:

    //    MagneticField* worldMagField;

    G4double  fWorldSizeX;
    G4double  fWorldSizeY;
    G4double  fWorldSizeZ;
    G4bool    fWorldVis;
    G4bool    fBuiltDetectors;
    G4double  fGriffinFwdBackPosition;
    G4int     fDetectorShieldSelect ;
    G4double  fDetectorRadialDistance ;
    G4int     fExtensionSuppressorLocation ;
    G4int     fCustomDetectorNumber ;
    G4int     fCustomDetectorPosition ;
    G4int     fCustomDetectorVal ;
    G4int     fHevimetSelector ;
    G4bool    fUseTigressPositions;

    // Box
    G4String           fBoxMat;
    G4double           fBoxThickness;
    G4ThreeVector      fBoxInnerDimensions;
    G4ThreeVector      fBoxColour;

    G4Box*             fSolidWorld;    //pointer to the solid World
    G4LogicalVolume*   fLogicWorld;    //pointer to the logical World
    G4VPhysicalVolume* fPhysiWorld;    //pointer to the physical World

    // Grid
    G4String           fGridMat;
    G4double           fGridSize;
    G4ThreeVector      fGridDimensions;
    G4ThreeVector      fGridColour;
    G4ThreeVector      fGridOffset;

    void DefineSuppressedParameters();
    void DefineMaterials();

    G4double fCoords[20][5];
    G4bool        fSetGenericTargetMaterial;
    G4bool        fSetGenericTargetDimensions;
    G4bool        fSetGenericTargetPosition;
    G4String      fGenericTargetMaterial;
    G4ThreeVector fGenericTargetDimensions;
    G4ThreeVector fGenericTargetPosition;

    G4bool        fSetFieldBoxMaterial;
    G4bool        fSetFieldBoxDimensions;
    G4bool        fSetFieldBoxPosition;
    G4bool        fSetFieldBoxMagneticField;
    G4String      fFieldBoxMaterial;
    G4ThreeVector fFieldBoxDimensions;
    G4ThreeVector fFieldBoxPosition;
    G4ThreeVector fFieldBoxMagneticField;

    G4String fMatWorldName;

    DetectorMessenger* fDetectorMessenger;

    G4ThreeVector fDescantRotation;
    G4String fDescantColor;

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


