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

#ifndef DETECTIONSYSTEMGRIFFIN_HH
#define DETECTIONSYSTEMGRIFFIN_HH

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectionSystemGriffin
{
public:
    DetectionSystemGriffin(G4int sel, G4int suppSwitch, G4double detRad, G4int hevimetSel );
    ~DetectionSystemGriffin();

    void Build() ;
    // For detector specific dead layers
    void BuildDeadLayerSpecificCrystal(G4int det);

    void BuildEverythingButCrystals(G4int det);
    G4double GetCrystalDistanceFromOrigin() {return fCrystalDistFromOrigin;}

    G4double TransX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double TransY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double TransZ(G4double x, G4double z, G4double theta);

    // For detector specific dead layers
    G4int PlaceDeadLayerSpecificCrystal(G4LogicalVolume* exp_hall_log, G4int detector_number, G4int position_number, G4bool posTigress);
    G4int PlaceEverythingButCrystals(G4LogicalVolume* exp_hall_log, G4int detector_number, G4int position_number, G4bool posTigress);

	 void SetDeadLayer(G4int detNum, G4int cryNum, G4double deadLayer) { fGriffinDeadLayers[detNum][cryNum] = deadLayer; }

private:
    G4String fSdName0;
    G4String fSdName1;
    G4String fSdName2;
    G4String fSdName3;
    G4String fSdName4;
    G4String fSdName5;
    G4String fColNameGe;
    G4String fColNameLeftCasing;
    G4String fColNameRightCasing;
    G4String fColNameLeftExtension;
    G4String fColNameRightExtension;
    G4String fColNameBackPlug;

    G4bool fIncludeExtensionSuppressors;
    G4bool fIncludeSideSuppressors;
    G4bool fIncludeBackSuppressors;

    G4String fBackSuppressorMaterial;
    G4String fBGOMaterial;

    G4RotationMatrix* fRotateNull;
    G4ThreeVector fMoveNull;

    G4double fCutClearance;
    G4double fExtraCutLength;

    G4bool fSurfCheck;

    double fCoords[20][5];
    
    //Jun 21, 2005, Epapr7.80: modification to the prototype suppressor shield,
    //since the hevimet collimators they provided were too big, and therefore the
    //extensions do not come far enough forward
    G4double fExtensionAccidentalBackShift;

    G4String fHevimetChoice;

    //Settings for the suppressor shield
    G4double fThicknessDouble;
    G4double fSideThicknessDouble;
    G4double fExtensionThicknessDouble;

    G4double fRadialDistance;
    G4int fNumberOfDetectors;
    G4double fOriginToCrystalDistance;
    G4bool fDeadLayerIncludeFlag;
    G4double fInnerDeadLayerThickness;
    G4double fOuterDeadLayerThickness;

    //Cold Finger
    G4double fColdFingerOuterShellRadius;
    G4double fColdFingerShellThickness;
    G4double fColdFingerShellLength;
    G4double fRearPlateThickness;
    G4double fColdFingerEndPlateThickness;
    G4double fColdFingerRadius;
    G4double fColdFingerSpace;
    G4double fColdFingerLength;
    G4double fCoolantLength;
    G4double fCoolantRadius;
    G4double fCoolantThickness;

    //New cooling structures (Jan 2005)
    G4double fExtraBlockThickness;
    G4double fExtraBlockDistanceFromBackPlate;
    G4double fExtraBlockInnerDiameter;
    G4double fTrianglePostsDistanceFromCrystals;
    G4double fTrianglePostStartingDepth;
    G4double fFetAirHoleRadius;
    G4double fCoolingSideBlockThickness;
    G4double fCoolingSideBlockWidth;
    G4double fCoolingSideBlockHorizontalDepth;
    G4double fStructureMatColdFingerThickness;
    G4double fStructureMatColdFingerRadius;
    G4double fCoolingBarThickness;
    G4double fCoolingBarWidth;

    G4int fSuppressorDesign;
    G4int fSuppressorPositionSelector;
    G4int fHevimetSelector;

    G4double fForwardInnerRadius;
    G4double fBackInnerRadius;

    G4double fSuppressorForwardRadius;
    G4double fSuppressorBackRadius;


    // For the optimization of the depth segmentation
    G4double fDepthSegmentationAdjustment;

    //the germanium detector's variables (for one "leaf" of the clover)
    G4double fGermaniumOuterRadius;
    G4double fGermaniumHoleRadius;
    G4double fGermaniumWidth;
    G4double fGermaniumLength;
    G4double fGermaniumHoleDistFromFace;
    G4double fGermaniumDistFromCanFace;
    G4double fGermaniumBentLength;
    G4double fGermaniumShift;    		//the amount by which the two sides adjacent to the other
    //germanium crystals are cut closer to the center

    G4double fGermaniumSeparation;  		//the space between the quarter detectors
    G4double fInterCrystalElectrodeMatThickness;
    G4double fElectrodeMatStartingDepth;
    G4double fGermaniumCornerConeEndLength; 	//the ending length of the cones
    //at the corners of each leaf
    
    //basic germanium detector's values
    G4double fDetectorBlockLength;
    G4double fDetectorBlockHeight;
    G4double fDetectorBlockTrapezoidalHeight;

    G4double fDetectorTotalWidth;
    G4double fDetectorTotalLength;
    G4double fBentEndAngle;
    G4double fBentEndLength;
    G4double fCanFaceThickness;
    G4double fCanSideThickness;

    G4double fHevimetTipThickness;
    G4double fHevimetTipAngle;

    //Values for the BGO
    G4double fSuppressorCutExtension;
    G4double fSuppressorShellThickness;
    
    G4double fBackBGOThickness;
    G4double fBGOChoppedTip;
    G4double fBGOMovableSpace;
    G4double fSideSuppressorBackOffset;
    G4double fSideBGOThickness;
    G4double fBGOCanSeperation;
    G4double fSideBGOLength;
    G4double fSideSuppressorLength;
    G4double fBGOTrapLength;
    G4double fSuppressorExtensionThickness;
    G4double fSuppressorExtensionLength;
    G4double fSuppressorExtensionAngle;
    
    G4double fSuppressorExtensionLengthDet ;

    //Values for the HeavyMet
    G4double fHeavyMetThickness;
    G4double fHeavyMetInsideAngle;

    G4double fAirBoxFrontWidth;
    G4double fAirBoxFrontLength;
    G4double fAirBoxBackLength;
    
    G4double fAirBoxBackLengthDet ;
    G4double fAirBoxFrontLengthDet ;
    G4double fAirBoxFrontWidthDet ;

    G4double fShift;
    G4double fSuppShift ;

    G4int fCopyNumber;
    G4int fCopyNumberTwo;
    G4int fGermaniumCopyNumber;
    G4int fLeftSuppressorSideCopyNumber;
    G4int fRightSuppressorSideCopyNumber;
    G4int fLeftSuppressorExtensionCopyNumber;
    G4int fRightSuppressorExtensionCopyNumber;
    G4int fBackSuppressorCopyNumber;
    
    G4double fRhombiDiameter;
    G4double fNewRhombiRadius;

    G4double fNewRhombiRadiusDet;

    G4double fDetectorPositionShift;
    G4double fAppliedBackShift;

    G4int fGermaniumSelector;
    G4int fCanSelector;
    G4int fBGOSelector;
    G4int fColdFingerSelector;

    //Redacted parameters/////////////
    G4double fDetectorPlacementCxn;
    G4double fTrianglePostDim;
    G4double fSuppressorExtRightX;
    G4double fSuppressorExtRightY;
    G4double fSuppressorExtRightZ;
    G4double fSuppressorExtLeftX;
    G4double fSuppressorExtLeftY;
    G4double fSuppressorExtLeftZ;
    G4double fWedgeDim;
    G4double fQuarterDetectorCxn;
    G4double fQuarterDetectorCxnB;
    G4String fElectrodeMaterial;
    G4String fStructureMaterial;

    // Assembly volumes
    G4AssemblyVolume* fAssembly;
    G4AssemblyVolume* fGermaniumAssembly;
    G4AssemblyVolume* fLeftSuppressorCasingAssembly;
    G4AssemblyVolume* fRightSuppressorCasingAssembly;
    G4AssemblyVolume* fLeftSuppressorExtensionAssembly;
    G4AssemblyVolume* fRightSuppressorExtensionAssembly;
    G4AssemblyVolume* fSuppressorBackAssembly;
    // For detector specific dead layers
    G4AssemblyVolume* fAssemblyCry[4];
    G4AssemblyVolume* fGermaniumAssemblyCry[4];
    G4AssemblyVolume* fLeftSuppressorCasingAssemblyCry[4];
    G4AssemblyVolume* fRightSuppressorCasingAssemblyCry[4];
    G4AssemblyVolume* fLeftSuppressorExtensionAssemblyCry[4];
    G4AssemblyVolume* fRightSuppressorExtensionAssemblyCry[4];
    G4AssemblyVolume* fSuppressorBackAssemblyCry[4];
    G4AssemblyVolume* fSuppressorShellAssembly;
    G4AssemblyVolume* fBackAndSideSuppressorShellAssembly ;
    G4AssemblyVolume* fExtensionSuppressorShellAssembly ;
    G4AssemblyVolume* fHevimetAssembly ;
    

    // Logical volumes
    G4LogicalVolume* fAirBoxLog;

    // methods to construct all of the components of the detector
    void ConstructNewSuppressorCasingWithShells(G4int det);
    void BuildelectrodeMatElectrodes();
    void ConstructComplexDetectorBlockWithDeadLayer();

    // For detector specific dead layers
    void ConstructComplexDetectorBlockWithDetectorSpecificDeadLayer(G4int det, G4int cry);
    void ConstructNewSuppressorCasingDetectorSpecificDeadLayer(G4int det, G4int cry);
    void ConstructDetector();
    void ConstructComplexDetectorBlock();
    void ConstructColdFinger();
    void ConstructNewHeavyMet();

    //LogicalVolumes used in ConstructBasicDetectorBlock
    G4LogicalVolume* fGermaniumBlockLog;

    //Logical Volumes used in ConstructComplexDetectorBlock:
    G4LogicalVolume* fGermaniumBlock1Log;
    G4LogicalVolume* fGermaniumHoleLog;
    G4LogicalVolume* fInnerDeadLayerLog;
    G4LogicalVolume* fInnerDeadLayerCapLog;
    G4LogicalVolume* fOuterDeadLayerLog;
    
    G4LogicalVolume* fInterCrystalElectrodeMatBackLog;
    G4LogicalVolume* fInterCrystalElectrodeMatFrontLog;
    
    //Logical Volumes used in ConstructBGOCasing:
    G4LogicalVolume* fBackBGOLog;
    G4LogicalVolume* fBGOCasingLog;

    //Logical Volumes used in ConstructNewSuppressorCasing:
    G4LogicalVolume* fBackQuarterSuppressorShellLog;
    G4LogicalVolume* fRightSuppressorShellLog;
    G4LogicalVolume* fLeftSuppressorShellLog;
    G4LogicalVolume* fRightSuppressorShellExtensionLog;
    G4LogicalVolume* fLeftSuppressorShellExtensionLog;

    G4LogicalVolume* fCapForRightSuppressorLog;

    G4LogicalVolume* fBackQuarterSuppressorLog;
    G4LogicalVolume* fRightSuppressorLog;
    G4LogicalVolume* fLeftSuppressorLog;
    G4LogicalVolume* fRightSuppressorExtensionLog;
    G4LogicalVolume* fLeftSuppressorExtensionLog;
    
    //Logical Volumes used in ConstructDetector:
    G4LogicalVolume* fFrontFaceLog;
    G4LogicalVolume* fRightBentPieceLog;
    G4LogicalVolume* fLeftBentPieceLog;
    G4LogicalVolume* fTopBentPieceLog;
    G4LogicalVolume* fBottomBentPieceLog;
    G4LogicalVolume* fRightWedgeLog;
    G4LogicalVolume* fLeftWedgeLog;
    G4LogicalVolume* fTopWedgeLog;
    G4LogicalVolume* fBottomWedgeLog;
    G4LogicalVolume* fUpperRightConeLog;
    G4LogicalVolume* fLowerRightConeLog;
    G4LogicalVolume* fUpperLeftConeLog;
    G4LogicalVolume* fLowerLeftConeLog;
    G4LogicalVolume* fUpperRightTubeLog;
    G4LogicalVolume* fLowerRightTubeLog;
    G4LogicalVolume* fUpperLeftTubeLog;
    G4LogicalVolume* fLowerLeftTubeLog;
    G4LogicalVolume* fRightSidePanelLog;
    G4LogicalVolume* fLeftSidePanelLog;
    G4LogicalVolume* fTopSidePanelLog;
    G4LogicalVolume* fBottomSidePanelLog;
    G4LogicalVolume* fRearPlateLog;
    G4LogicalVolume* fFingerShellLog;
    G4LogicalVolume* fTankLog;
    G4LogicalVolume* fTankLid1Log;
    G4LogicalVolume* fTankLid2Log;
    G4LogicalVolume* fTankLiquidLog;

    //Logical Volumes used in ConstructColdFinger:
    G4LogicalVolume* fEndPlateLog;
    G4LogicalVolume* fFingerLog;
    G4LogicalVolume* fExtraColdBlockLog;
    G4LogicalVolume* fTrianglePostLog;
    G4LogicalVolume* fFetAirHoleLog;
    G4LogicalVolume* fCoolingBarLog;
    G4LogicalVolume* fCoolingSideBlockLog;
    G4LogicalVolume* fStructureMatColdFingerLog;
    
    //Logical Volumes used in ConstructNewHeavyMet:
    G4LogicalVolume* fHevimetLog;

    //internal methods for ConstructCan()
    G4Box* SquareFrontFace();
    G4Trap* CornerWedge();
    G4Para* BentSidePiece();
    G4Box* OtherBentSidePiece();
    G4Cons* RoundedEndEdge();
    G4Tubs* CornerTube();
    G4Box* SidePanel();
    G4SubtractionSolid* RearPlate();
    G4Tubs* ColdFingerShell();
    G4Tubs* LiquidNitrogenTank();
    G4Tubs* LiquidNitrogenTankLid();
    G4Tubs* LiquidNitrogen();


    //internal methods for ConstructBasicDetectorBlock()
    G4Box* RectangularSegment();
    G4Trd* TrapezoidalSegment();

    //internal methods for ConstructComplexDetectorBlock()
    G4SubtractionSolid* QuarterDetector();
    // For detector specific dead layers
    G4SubtractionSolid* QuarterSpecificDeadLayerDetector(G4int det, G4int cry);

    //internal methods for ConstructComplexDetectorBlockWithPlastic()
    G4UnionSolid* InterCrystalelectrodeMatBack();
    G4UnionSolid* InterCrystalelectrodeMatFront();
    
    //internal methods for ConstructColdFinger()
    G4Tubs* AirHole();
    G4Tubs* AirHoleCut();
    G4Box* EndPlate();
    G4Tubs* Finger();
    G4SubtractionSolid* ExtraColdBlock();
    G4Trd* TrianglePost();
    G4Box* CoolingBar();
    G4Box* CoolingSideBlock();
    G4Tubs* StructureMatColdFinger();
    
    //internal methods for ConstructBGOCasing()
    G4SubtractionSolid* BackBGO();
    G4SubtractionSolid* BGOCasing();
    G4SubtractionSolid* FrontBGO();
    G4Trd* SideBGO();

    //internal methods for ConstructNewSuppressorCasing()
    G4SubtractionSolid* BackSuppressorQuarter();
    G4SubtractionSolid* FrontSlantSuppressor(G4String sidePosition, G4bool choppingSuppressor) ;
    G4SubtractionSolid* SideSuppressorExtension(G4String sidePosition, G4bool choppingSuppressor) ;
    G4Trap*             SideSuppressorExtensionUncut() ;

    //internal methods for New SuppressorCasingWithShells
    G4SubtractionSolid* ShellForBackSuppressorQuarter();
    G4SubtractionSolid* ShellForFrontSlantSuppressor(G4String sidePosition) ;
    G4SubtractionSolid* ShellForSuppressorExtension(G4String sidePosition);

    //internal methods for ConstructNewHeavyMet()
    G4SubtractionSolid* NewHeavyMet();

    G4String fCrystalMaterial;
    G4String fCanMaterial;
    G4String fVacuumMaterial;
    G4double fCrystalLengthX;
    G4double fCrystalLengthY;
    G4double fCrystalLengthZ;
    G4double fCrystalInnerRadius;
    G4double fCrystalOuterRadius;
    G4double fCanThickness;
    G4double fCanInnerRadius;
    G4double fCanLidInnerRadius;
    G4double fCanLidOuterRadius;
    G4double fCanFrontLidThickness;
    G4double fCanBackLidThickness;
    G4double fCanFaceDistFromOrigin;
    G4double fCrystalDistFromCanFace;
    G4double fCrystalDistFromCanBack;
    G4double fCanLengthZ;
    G4double fCrystalDistFromOrigin;

    // For detector specific dead layers
    G4double fGriffinDeadLayers[16][4];
    G4Colour fGriffinCrystalColours[4];
    G4Colour fGriffinDeadLayerColours[4];

    // internal methods
    void BuildOneDetector(G4int det);
    //    void PlaceDetector(G4int detectorNumber);
};

#endif

