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
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectionSystemSpiceV02_h
#define DetectionSystemSpiceV02_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#define AL_COL 0.5,0.5,0.5

class DetectionSystemSpiceV02
{
public:
    DetectionSystemSpiceV02();
    ~DetectionSystemSpiceV02();
    
    G4int Build() ; //G4SDManager* mySDman);
    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate, G4int detector_number);
    //G4double const GetDetectorLengthOfUnitsCM() {return this->can_length_z;};

    //-----------------------------//
    // parameters for the square   //
    // planar detector crystal     //
    //-----------------------------//
    G4double squareDetCrystalLength;
    G4double squareDet90DegCrystalRadius;
    G4double squareDet45DegCrystalRadius;
    G4double squareDetCrystalThickness;
    G4double detector_face_target_distance;

    //------------------------------------------------//
    // logical and physical volumes
    //------------------------------------------------//
private:
    G4AssemblyVolume* assembly;
    G4AssemblyVolume* assemblySi;
    //    SensitiveDetector* crystal_block_SD;

    G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

    G4LogicalVolume* detector_casing_side_log;
    G4LogicalVolume* detector_casing_back_log;
    G4LogicalVolume* crystal_block_log;


    //--------------------------------------------------------//
    // SPICE physical properties
    // OBS: crystal properties are public, others are private
    //--------------------------------------------------------//

    G4String casing_material;
    G4String wafer_material;

    //-----------------------------//
    // parameters for the square   //
    // planar detector casing      //
    //-----------------------------//

    G4double squareDetCasingLength;
    G4double squareDetCasingThickness;

    G4int siDetCopyNumber;
    G4int siDetCasingSideCopyNumber;
    G4int siDetCasingBackCopyNumber;

    //------------------------------------------------//
    // internal methods in Build()
    //------------------------------------------------//

    G4int BuildSiliconWafer();
    G4int BuildAluminiumCasingSide();
    G4int BuildAluminiumCasingBack();

    //    void PlaceSiliconWafer();
    //    void PlaceAluminiumCasingSide();
    //    void PlaceAluminiumCasingBack();

    //    void BuildSPICETargetChamber();

    //    G4ThreeVector Translate45DegDetector(G4int);
    //    G4ThreeVector Translate90DegDetector(G4int);
    //    G4RotationMatrix* RotateDetector(G4int);
    //------------------------------------------------//
    // public methods
    //------------------------------------------------//
    G4Box*              BuildCrystal();
    G4SubtractionSolid* BuildCasingSide();
    G4Box*              BuildCasingBack();
};

#endif










//#ifndef SPICEDetectionSystem_h
//#define SPICEDetectionSystem_h 1

//class DetectableSD;
//class SPICETargetChamber;

//#include "G4SDManager.hh"
//#include "DetectableSD.hh"

//#include "DetectorConstruction.hh"

//#include "G4LogicalVolume.hh"
//#include "G4VPhysicalVolume.hh"
//#include "SPICETargetChamber.hh"
//#include "G4AssemblyVolume.hh"

//#include "SPICE_globals.hh"

//#define AL_COL 0.5,0.5,0.5

/////////////////////////////////////////////////////////////////////////
//// SPICEDetectionSystem
/////////////////////////////////////////////////////////////////////////
//class SPICEDetectionSystem
//{

//public:
//  DetectableSD* myCrystal_block_SD;

//public:
//  SPICEDetectionSystem(G4SDManager*, EventDataOutput*);
//  ~SPICEDetectionSystem();

//public:
//  void Build(G4LogicalVolume* exp_hall_log);

//  G4double GetDetector2Origin(SPICEDetectionSystem*);

//private:
//  Materials* set_of_materials;
//  G4SDManager* mySDman;
//  G4LogicalVolume* expHallLog;

//private:
//  DetectableSD* crystal_block_SD;
//  EventDataOutput* my_data_output;

//  //------------------------------------------------//
//  // logical and physical volumes
//  //------------------------------------------------//
//private:
//  G4LogicalVolume* detector_casing_side_log;
//  G4LogicalVolume* detector_casing_back_log;
//  G4LogicalVolume* crystal_block_log;
//  G4VPhysicalVolume* crystal_45Deg_phys;
//  G4VPhysicalVolume* crystal_90Deg_phys;
//  G4VPhysicalVolume* casing_side_45Deg_phys;
//  G4VPhysicalVolume* casing_side_90Deg_phys;
//  G4VPhysicalVolume* casing_back_45Deg_phys;
//  G4VPhysicalVolume* casing_back_90Deg_phys;

//  //--------------------------------------------------------//
//  // SPICE physical properties
//  // OBS: crystal properties are public, others are private
//  //--------------------------------------------------------//
//private:
//  G4String casing_material;
//  G4String wafer_material;

//  //-----------------------------//
//  // parameters for the square   //
//  // planar detector crystal     //
//  //-----------------------------//
//public:
//  G4double squareDetCrystalLength;
//  G4double squareDet90DegCrystalRadius;
//  G4double squareDet45DegCrystalRadius;
//  G4double squareDetCrystalThickness;
//  G4double detector_face_target_distance;

//  //-----------------------------//
//  // parameters for the square   //
//  // planar detector casing      //
//  //-----------------------------//
//private:
//  G4double squareDetCasingLength;
//  G4double squareDetCasingThickness;

//  G4int siDetCopyNumber;
//  G4int siDetCasingSideCopyNumber;
//  G4int siDetCasingBackCopyNumber;

//  //------------------------------------------------//
//  // internal methods in Build()
//  //------------------------------------------------//
//private:
//  void BuildSiliconWafer();
//  void BuildAluminiumCasingSide();
//  void BuildAluminiumCasingBack();
//  void PlaceSiliconWafer();
//  void PlaceAluminiumCasingSide();
//  void PlaceAluminiumCasingBack();
//  void BuildSPICETargetChamber();

//  G4ThreeVector Translate45DegDetector(G4int);
//  G4ThreeVector Translate90DegDetector(G4int);
//  G4RotationMatrix* RotateDetector(G4int);
//  //------------------------------------------------//
//  // public methods
//  //------------------------------------------------//
//public:
//  DetectableSD* getDetectableSD();

//};

//#endif
