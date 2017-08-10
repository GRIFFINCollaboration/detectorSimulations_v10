// Includes Physical Constants and System of Units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemGriffin.hh"

#include "DetectionSystem8pi.hh"
#include "Apparatus8piVacuumChamber.hh"
#include "Apparatus8piVacuumChamberAuxMatShell.hh"
#include "ApparatusGriffinStructure.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


///////////////////////////////////////////////////////////////////////
// The ::DetectionSystemGriffin constructor instatiates all the
// Logical and Physical Volumes used in the detector geometery, and the
// ::~DetectionSystemGriffin destructor deletes them from the stack
// when they go out of scope
///////////////////////////////////////////////////////////////////////
DetectionSystemGriffin::DetectionSystemGriffin( G4int sel, G4int suppSwitch, G4double detRad, G4int hevimetSel ):
    ///////////////////////////////////////////////////////////////////
    // LogicalVolumes
    ///////////////////////////////////////////////////////////////////

    //LogicalVolumes used in ConstructBasicDetectorBlock
    fGermaniumBlockLog(NULL),

    //Logical Volumes used in ConstructComplexDetectorBlock:
    fGermaniumBlock1Log(NULL), fGermaniumHoleLog(NULL),
    fInnerDeadLayerLog(NULL), fInnerDeadLayerCapLog(NULL),
    fOuterDeadLayerLog(NULL),
    fInterCrystalElectrodeMatBackLog(NULL), fInterCrystalElectrodeMatFrontLog(NULL),
    
    //Logical Volumes used in ConstructBGOCasing:
    fBackBGOLog(NULL), fBGOCasingLog(NULL),

    //Logical Volumes used in ConstructNewSuppressorCasing:
    fBackQuarterSuppressorShellLog(NULL), fRightSuppressorShellLog(NULL),
    fLeftSuppressorShellLog(NULL),
    fRightSuppressorShellExtensionLog(NULL), fLeftSuppressorShellExtensionLog(NULL),
    fCapForRightSuppressorLog(NULL), fBackQuarterSuppressorLog(NULL),
    fRightSuppressorLog(NULL), fLeftSuppressorLog(NULL),
    fRightSuppressorExtensionLog(NULL), fLeftSuppressorExtensionLog(NULL),
    
    //Logical Volumes used in ConstructDetector:
    fFrontFaceLog(NULL), fRightBentPieceLog(NULL),
    fLeftBentPieceLog(NULL), fTopBentPieceLog(NULL),
    fBottomBentPieceLog(NULL), fRightWedgeLog(NULL),
    fLeftWedgeLog(NULL), fTopWedgeLog(NULL),
    fBottomWedgeLog(NULL), fUpperRightConeLog(NULL),
    fLowerRightConeLog(NULL), fUpperLeftConeLog(NULL),
    fLowerLeftConeLog(NULL), fUpperRightTubeLog(NULL),
    fLowerRightTubeLog(NULL), fUpperLeftTubeLog(NULL),
    fLowerLeftTubeLog(NULL), fRightSidePanelLog(NULL),
    fLeftSidePanelLog(NULL), fTopSidePanelLog(NULL),
    fBottomSidePanelLog(NULL), fRearPlateLog(NULL),
    fFingerShellLog(NULL), fTankLog(NULL),
    fTankLid1Log(NULL), fTankLid2Log(NULL),
    fTankLiquidLog(NULL),

    //Logical Volumes used in ConstructColdFinger:
    fEndPlateLog(NULL), fFingerLog(NULL),
    fExtraColdBlockLog(NULL), fTrianglePostLog(NULL),
    fFetAirHoleLog(NULL), fCoolingBarLog(NULL),
    fCoolingSideBlockLog(NULL), fStructureMatColdFingerLog(NULL),
    
    //Logical Volumes used in ConstructNewHeavyMet:
    fHevimetLog(NULL)

{
    G4ThreeVector myMoveNull(0, 0, 1);
    fMoveNull = myMoveNull * 0;
    G4RotationMatrix* myRotateNull = new G4RotationMatrix;
    fRotateNull = myRotateNull;

    /////////////////////////////////////////////////////////////////////
    // Coords for GRIFFIN
    // Note that the GRIFFIN lampshade angles are rotated by 45 degrees with respect to those of TIGRESS.
    // Modified coords for TIGRESS are below!
    /////////////////////////////////////////////////////////////////////
    // theta
    fCoords[0][0] 	= 45.0;
    fCoords[1][0] 	= 45.0;
    fCoords[2][0] 	= 45.0;
    fCoords[3][0] 	= 45.0;
    fCoords[4][0] 	= 90.0;
    fCoords[5][0] 	= 90.0;
    fCoords[6][0] 	= 90.0;
    fCoords[7][0] 	= 90.0;
    fCoords[8][0] 	= 90.0;
    fCoords[9][0] 	= 90.0;
    fCoords[10][0] 	= 90.0;
    fCoords[11][0] 	= 90.0;
    fCoords[12][0] 	= 135.0;
    fCoords[13][0] 	= 135.0;
    fCoords[14][0] 	= 135.0;
    fCoords[15][0] 	= 135.0;
    // phi
    fCoords[0][1] 	= 67.5;
    fCoords[1][1] 	= 157.5;
    fCoords[2][1] 	= 247.5;
    fCoords[3][1] 	= 337.5;
    fCoords[4][1] 	= 22.5;
    fCoords[5][1] 	= 67.5;
    fCoords[6][1] 	= 112.5;
    fCoords[7][1] 	= 157.5;
    fCoords[8][1] 	= 202.5;
    fCoords[9][1] 	= 247.5;
    fCoords[10][1] 	= 292.5;
    fCoords[11][1] 	= 337.5;
    fCoords[12][1] 	= 67.5;
    fCoords[13][1] 	= 157.5;
    fCoords[14][1] 	= 247.5;
    fCoords[15][1] 	= 337.5;
    // yaw (alpha)
    fCoords[0][2] 	= 0.0;
    fCoords[1][2] 	= 0.0;
    fCoords[2][2] 	= 0.0;
    fCoords[3][2] 	= 0.0;
    fCoords[4][2] 	= 0.0;
    fCoords[5][2] 	= 0.0;
    fCoords[6][2] 	= 0.0;
    fCoords[7][2] 	= 0.0;
    fCoords[8][2] 	= 0.0;
    fCoords[9][2] 	= 0.0;
    fCoords[10][2] 	= 0.0;
    fCoords[11][2] 	= 0.0;
    fCoords[12][2] 	= 0.0;
    fCoords[13][2] 	= 0.0;
    fCoords[14][2] 	= 0.0;
    fCoords[15][2] 	= 0.0;
    // pitch (beta)
    fCoords[0][3] 	= -45.0;
    fCoords[1][3] 	= -45.0;
    fCoords[2][3] 	= -45.0;
    fCoords[3][3] 	= -45.0;
    fCoords[4][3] 	= 0.0;
    fCoords[5][3] 	= 0.0;
    fCoords[6][3] 	= 0.0;
    fCoords[7][3] 	= 0.0;
    fCoords[8][3] 	= 0.0;
    fCoords[9][3] 	= 0.0;
    fCoords[10][3] 	= 0.0;
    fCoords[11][3] 	= 0.0;
    fCoords[12][3] 	= 45.0;
    fCoords[13][3] 	= 45.0;
    fCoords[14][3] 	= 45.0;
    fCoords[15][3] 	= 45.0;
    // roll (gamma)
    fCoords[0][4] 	= 67.5;
    fCoords[1][4] 	= 157.5;
    fCoords[2][4] 	= 247.5;
    fCoords[3][4] 	= 337.5;
    fCoords[4][4] 	= 22.5;
    fCoords[5][4] 	= 67.5;
    fCoords[6][4] 	= 112.5;
    fCoords[7][4] 	= 157.5;
    fCoords[8][4] 	= 202.5;
    fCoords[9][4] 	= 247.5;
    fCoords[10][4] 	= 292.5;
    fCoords[11][4] 	= 337.5;
    fCoords[12][4] 	= 67.5;
    fCoords[13][4] 	= 157.5;
    fCoords[14][4] 	= 247.5;
    fCoords[15][4] 	= 337.5;

    // Note that the GRIFFIN lampshade angles are rotated by 45 degrees with respect to those of TIGRESS.
    // Uncomment this for TIGRESS!
    // phi
    //  fCoords[0][1] 	= fCoords[0][1] - 45.0;
    //  fCoords[1][1] 	= fCoords[1][1] - 45.0;
    //  fCoords[2][1] 	= fCoords[2][1] - 45.0;
    //  fCoords[3][1] 	= fCoords[3][1] - 45.0;
    //  fCoords[12][1] 	= fCoords[12][1] - 45.0;
    //  fCoords[13][1] 	= fCoords[13][1] - 45.0;
    //  fCoords[14][1] 	= fCoords[14][1] - 45.0;
    //  fCoords[15][1] 	= fCoords[15][1] - 45.0;
    //  // roll (gamma)
    //  fCoords[0][4] 	= fCoords[0][4] - 45.0;
    //  fCoords[1][4] 	= fCoords[1][4] - 45.0;
    //  fCoords[2][4] 	= fCoords[2][4] - 45.0;
    //  fCoords[3][4] 	= fCoords[3][4] - 45.0;
    //  fCoords[12][4] 	= fCoords[12][4] - 45.0;
    //  fCoords[13][4] 	= fCoords[13][4] - 45.0;
    //  fCoords[14][4] 	= fCoords[14][4] - 45.0;
    //  fCoords[15][4] 	= fCoords[15][4] - 45.0;

    // Surfaces were check on June 1st 2015, no geometry overlaps!
    fSurfCheck = false;

    /////////////////////////////////////////////////////////////////////
    // GRIFFIN/TIGRESS Detector Properties
    /////////////////////////////////////////////////////////////////////
    fCutClearance = 0.01*mm;
    fExtraCutLength = 10.0*mm;

    // Killing Suppressors!
    fIncludeExtensionSuppressors = true;
    fIncludeSideSuppressors      = true;
    fIncludeBackSuppressors      = true;		//BM: right place to turn on suppressors?

    //////
    G4int canPos                               = 1;
    G4int canPositionSpecifier                 = canPos;  	// 1=fully forward, 2=flush with BGO, 3=fully back
    G4int detectorFormat                       = 2; 		   // 1=simple, 2=complex, 3=none
    G4String includeDetectorCan                = "y"; 		// "y" or "n"
    G4String includeBGO                        = "y";
    G4String includeColdFinger                 = "y";
    ///// Included in new system

    if(suppSwitch  == 1) { // Include Suppressors
        fIncludeExtensionSuppressors = true;
        fIncludeSideSuppressors      = true;
        fIncludeBackSuppressors      = true;
    } else if(suppSwitch  == -1) { // just side and back
        fIncludeExtensionSuppressors = false;
        fIncludeSideSuppressors      = true;
        fIncludeBackSuppressors      = true;
    } else { // Dont include Suppressors
        fIncludeExtensionSuppressors = false;
        fIncludeSideSuppressors      = false;
        fIncludeBackSuppressors      = false;
        includeBGO						 = "n" ;
    }

    // GRIFFIN Dead Layers
    // Generated Oct 10th, 2013
    //  fGriffinDeadLayers[0][0] = 1.38*mm;
    //  fGriffinDeadLayers[0][1] = 1.30*mm;
    //  fGriffinDeadLayers[0][2] = 1.61*mm;
    //  fGriffinDeadLayers[0][3] = 1.84*mm;
    //  fGriffinDeadLayers[1][0] = 0.324*mm;
    //  fGriffinDeadLayers[1][1] = 0.506*mm;
    //  fGriffinDeadLayers[1][2] = 0.43*mm;
    //  fGriffinDeadLayers[1][3] = 0.40*mm;
    //  fGriffinDeadLayers[2][0] = 0.203*mm;
    //  fGriffinDeadLayers[2][1] = 0.457*mm;
    //  fGriffinDeadLayers[2][2] = 0.228*mm;
    //  fGriffinDeadLayers[2][3] = 0.614*mm;
    //  fGriffinDeadLayers[3][0] = 1.33*mm;
    //  fGriffinDeadLayers[3][1] = 1.17*mm;
    //  fGriffinDeadLayers[3][2] = 1.83*mm;
    //  fGriffinDeadLayers[3][3] = 1.38*mm;
    //  fGriffinDeadLayers[4][0] = 0.485*mm;
    //  fGriffinDeadLayers[4][1] = 0.488*mm;
    //  fGriffinDeadLayers[4][2] = 0.495*mm;
    //  fGriffinDeadLayers[4][3] = 0.508*mm;
    //  fGriffinDeadLayers[5][0] = 1.75*mm;
    //  fGriffinDeadLayers[5][1] = 1.48*mm;
    //  fGriffinDeadLayers[5][2] = 1.66*mm;
    //  fGriffinDeadLayers[5][3] = 1.47*mm;
    //  fGriffinDeadLayers[6][0] = 1.25*mm;
    //  fGriffinDeadLayers[6][1] = 1.31*mm;
    //  fGriffinDeadLayers[6][2] = 1.21*mm;
    //  fGriffinDeadLayers[6][3] = 1.35*mm;
    //  fGriffinDeadLayers[7][0] = 1.19*mm;
    //  fGriffinDeadLayers[7][1] = 2.07*mm;
    //  fGriffinDeadLayers[7][2] = 1.22*mm;
    //  fGriffinDeadLayers[7][3] = 1.39*mm;


    // Using average dead layers
    fGriffinDeadLayers[0][0] = 1.07*mm;
    fGriffinDeadLayers[0][1] = 1.07*mm;
    fGriffinDeadLayers[0][2] = 1.07*mm;
    fGriffinDeadLayers[0][3] = 1.07*mm;
    fGriffinDeadLayers[1][0] = 1.07*mm;
    fGriffinDeadLayers[1][1] = 1.07*mm;
    fGriffinDeadLayers[1][2] = 1.07*mm;
    fGriffinDeadLayers[1][3] = 1.07*mm;
    fGriffinDeadLayers[2][0] = 1.07*mm;
    fGriffinDeadLayers[2][1] = 1.07*mm;
    fGriffinDeadLayers[2][2] = 1.07*mm;
    fGriffinDeadLayers[2][3] = 1.07*mm;
    fGriffinDeadLayers[3][0] = 1.07*mm;
    fGriffinDeadLayers[3][1] = 1.07*mm;
    fGriffinDeadLayers[3][2] = 1.07*mm;
    fGriffinDeadLayers[3][3] = 1.07*mm;
    fGriffinDeadLayers[4][0] = 1.07*mm;
    fGriffinDeadLayers[4][1] = 1.07*mm;
    fGriffinDeadLayers[4][2] = 1.07*mm;
    fGriffinDeadLayers[4][3] = 1.07*mm;
    fGriffinDeadLayers[5][0] = 1.07*mm;
    fGriffinDeadLayers[5][1] = 1.07*mm;
    fGriffinDeadLayers[5][2] = 1.07*mm;
    fGriffinDeadLayers[5][3] = 1.07*mm;
    fGriffinDeadLayers[6][0] = 1.07*mm;
    fGriffinDeadLayers[6][1] = 1.07*mm;
    fGriffinDeadLayers[6][2] = 1.07*mm;
    fGriffinDeadLayers[6][3] = 1.07*mm;
    fGriffinDeadLayers[7][0] = 1.07*mm;
    fGriffinDeadLayers[7][1] = 1.07*mm;
    fGriffinDeadLayers[7][2] = 1.07*mm;
    fGriffinDeadLayers[7][3] = 1.07*mm;
    fGriffinDeadLayers[8][0] = 1.07*mm;
    fGriffinDeadLayers[8][1] = 1.07*mm;
    fGriffinDeadLayers[8][2] = 1.07*mm;
    fGriffinDeadLayers[8][3] = 1.07*mm;
    fGriffinDeadLayers[9][0] = 1.07*mm;
    fGriffinDeadLayers[9][1] = 1.07*mm;
    fGriffinDeadLayers[9][2] = 1.07*mm;
    fGriffinDeadLayers[9][3] = 1.07*mm;
    fGriffinDeadLayers[10][0] = 1.07*mm;
    fGriffinDeadLayers[10][1] = 1.07*mm;
    fGriffinDeadLayers[10][2] = 1.07*mm;
    fGriffinDeadLayers[10][3] = 1.07*mm;
    fGriffinDeadLayers[11][0] = 1.07*mm;
    fGriffinDeadLayers[11][1] = 1.07*mm;
    fGriffinDeadLayers[11][2] = 1.07*mm;
    fGriffinDeadLayers[11][3] = 1.07*mm;
    fGriffinDeadLayers[12][0] = 1.07*mm;
    fGriffinDeadLayers[12][1] = 1.07*mm;
    fGriffinDeadLayers[12][2] = 1.07*mm;
    fGriffinDeadLayers[12][3] = 1.07*mm;
    fGriffinDeadLayers[13][0] = 1.07*mm;
    fGriffinDeadLayers[13][1] = 1.07*mm;
    fGriffinDeadLayers[13][2] = 1.07*mm;
    fGriffinDeadLayers[13][3] = 1.07*mm;
    fGriffinDeadLayers[14][0] = 1.07*mm;
    fGriffinDeadLayers[14][1] = 1.07*mm;
    fGriffinDeadLayers[14][2] = 1.07*mm;
    fGriffinDeadLayers[14][3] = 1.07*mm;
    fGriffinDeadLayers[15][0] = 1.07*mm;
    fGriffinDeadLayers[15][1] = 1.07*mm;
    fGriffinDeadLayers[15][2] = 1.07*mm;
    fGriffinDeadLayers[15][3] = 1.07*mm;

    fGriffinCrystalColours[0] = G4Colour(0.0,0.0,1.0);
    fGriffinCrystalColours[1] = G4Colour(0.0,1.0,0.0);
    fGriffinCrystalColours[2] = G4Colour(1.0,0.0,0.0);
    fGriffinCrystalColours[3] = G4Colour(1.0,1.0,1.0);

    fGriffinDeadLayerColours[0] = G4Colour(0.0,0.0,0.50);
    fGriffinDeadLayerColours[1] = G4Colour(0.0,0.50,0.0);
    fGriffinDeadLayerColours[2] = G4Colour(0.50,0.0,0.0);
    fGriffinDeadLayerColours[3] = G4Colour(0.3,0.3,0.3);

    fBGOMaterial                  = "G4_BGO";
    fBackSuppressorMaterial 	    = "G4_CESIUM_IODIDE";

    if(sel == 0) {
        fSdName0 =               "/sd/allGriffinForward0";
        fSdName1 =               "/sd/allGriffinForward1";
        fSdName2 =               "/sd/allGriffinForward2";
        fSdName3 =               "/sd/allGriffinForward3";
        fSdName4 =               "/sd/allGriffinForward4";
        fSdName5 =               "/sd/allGriffinForward5";

        fColNameGe =             "CollectionGriffinForwardGe";
        fColNameLeftCasing =     "CollectionGriffinForwardLeftCasing";
        fColNameRightCasing =    "CollectionGriffinForwardRightCasing";
        fColNameLeftExtension =  "CollectionGriffinForwardLeftExtension";
        fColNameRightExtension = "CollectionGriffinForwardRightExtension";
        fColNameBackPlug =       "CollectionGriffinForwardBackPlug";
    } else if(sel == 1) {
        fSdName0 =               "/sd/allGriffinBack0";
        fSdName1 =               "/sd/allGriffinBack1";
        fSdName2 =               "/sd/allGriffinBack2";
        fSdName3 =               "/sd/allGriffinBack3";
        fSdName4 =               "/sd/allGriffinBack4";
        fSdName5 =               "/sd/allGriffinBack5";

        fColNameGe =             "CollectionGriffinBackGe";
        fColNameLeftCasing =     "CollectionGriffinBackLeftCasing";
        fColNameRightCasing =    "CollectionGriffinBackRightCasing";
        fColNameLeftExtension =  "CollectionGriffinBackLeftExtension";
        fColNameRightExtension = "CollectionGriffinBackRightExtension";
        fColNameBackPlug =       "CollectionGriffinBackBackPlug";
    } else {
        G4cout << "Error: 123109812 " << G4endl;
        exit(1);
    }

    // These keep track of the number of sensitive detectors. Initialized to zero
    fGermaniumCopyNumber                = 4000;
    fLeftSuppressorSideCopyNumber       = 4100;
    fRightSuppressorSideCopyNumber      = 4200;
    fLeftSuppressorExtensionCopyNumber  = 4300;
    fRightSuppressorExtensionCopyNumber = 4400;
    fBackSuppressorCopyNumber           = 4500;

    fSuppressorPositionSelector         = sel;  	// 0 = detector forward, 1 = detector back

    fHevimetSelector                    = hevimetSel ;

    fInnerDeadLayerThickness 	          = 1.0*mm;
    fDepthSegmentationAdjustment 	    = 0.0*cm;


    // Determine whether the structureMat shells around the suppressors will be included
    fForwardInnerRadius 		          = detRad ;
    fBackInnerRadius 		             = fForwardInnerRadius + 3.5*cm;
    fSuppressorForwardRadius            = 11.0*cm ;
    fSuppressorBackRadius               = fSuppressorForwardRadius + 3.5*cm ;

    // "Outer" dead layer description: thickness determined from R.G.Helmer paper
    fOuterDeadLayerThickness	          = 2.5*micrometer;

    //2.00 new: the thickness of the structureMat shell around the suppression material
    fSuppressorShellThickness 	       = 0.5*mm;

    //Jun 21, 2005: Epapr7.80: the suppressor extensions were accidentally made too large,
    //so they do not come as far forward as they should.  Shift for back, with hevimet was 1.0*cm
    //Jun 28, 2005: Just! for the back, no hevimet case, change this to 0.5*cm
    fExtensionAccidentalBackShift       = 0.5*cm;
    //Set (in this) whether the structureMat shells around the suppressors will be included
    fGermaniumOuterRadius 		          = 3.0*cm;	// Radius of the circular part of the germanium
    fGermaniumHoleRadius 		          = 0.5*cm;
    fGermaniumLength 		             = 9.0*cm;
    fGermaniumHoleDistFromFace 	       = 1.5*cm;
    fGermaniumDistFromCanFace 	       = 0.55*cm;
    fGermaniumSeparation 		          = 0.6*mm;  	// Separation between quarter detectors: changed Jan 2005 from 0.7mm

    // This has been disabled Apr 11, 2005 for version 7.69: setting this
    fInterCrystalElectrodeMatThickness  = 0.6*mm;
    fGermaniumWidth 		                = 56.5*mm;  	// Width of one quarter detector
    fGermaniumBentLength 		          = 3.62*cm;
    fGermaniumShift                     = 1.05*mm;  	// this can't be more than 2.75mm. It is the amount by which
    // one side is cut closer to the center than the other
    fGermaniumCornerConeEndLength       = 3.0*cm; 	// the ending length of the cones

    fElectrodeMatStartingDepth          = fGermaniumCornerConeEndLength-1*micrometer;
    fDetectorBlockLength                = 2.0*fGermaniumWidth; 	// Obsolete
    fDetectorBlockHeight                = fGermaniumLength - fGermaniumBentLength;
    fDetectorTotalWidth                 = 12.4*cm;             // Width of the detector can
    fDetectorTotalLength                = 16.0*cm;             // Length of the detector can
    fBentEndAngle                       = 22.5*M_PI/180;
    fBentEndLength                      = 4.325*cm;
    fCanFaceThickness                   = 0.15*cm;
    fCanSideThickness                   = 0.2*cm;
    fDetectorBlockTrapezoidalHeight     = fGermaniumBentLength;

    fColdFingerOuterShellRadius         = 2.5*cm;
    fColdFingerShellThickness           = 0.2*cm;
    fColdFingerShellLength              = 14.5*cm;

    // Values for the cold finger, and cooling structure (added Jan 2005)
    fRearPlateThickness                 = 0.2*cm;
    fColdFingerEndPlateThickness        = 10.0*mm; 	// Changed Jan 2005 from 9.0mm
    fColdFingerRadius                   = 9.0*mm; 	// Changed Jan 2005 from 8.25mm
    fColdFingerSpace                    = 1.0*cm;  	// Space between the cold plate and the germanium crystals
    fColdFingerLength                   = fColdFingerShellLength + fDetectorTotalLength
		- fCanFaceThickness - fGermaniumDistFromCanFace
		- fGermaniumLength - fColdFingerSpace
		- fColdFingerEndPlateThickness;
    // New stuff
    fExtraBlockThickness                = 9.0*mm;
    fExtraBlockDistanceFromBackPlate    = 10.0*mm;
    fExtraBlockInnerDiameter            = 34.0*mm;

    fTrianglePostsDistanceFromCrystals  = 2.0*mm;
    fTrianglePostStartingDepth          = 38.0*mm;
    fFetAirHoleRadius                   = 26.0*mm/2.0;
    fCoolingSideBlockThickness          = 15.0*mm;
    fCoolingSideBlockWidth              = 20.0*mm;
    fCoolingSideBlockHorizontalDepth    = 10.0*mm;
    fStructureMatColdFingerThickness    = 26.0*mm;
    fStructureMatColdFingerRadius       = 22.0*mm/2.0;
    fCoolingBarThickness                = 10.0*mm;
    fCoolingBarWidth                    = 8.0*mm ;

    // Liquid nitrogen tank stuff
    fCoolantLength                      = 39.3*cm;   // Length of the Liquid Nitrogen tank
    fCoolantRadius                      = 12.5*cm;
    fCoolantThickness                   = 0.5*cm;


    fHevimetTipThickness                = 1.27*cm;
    fHevimetTipAngle                    = fBentEndAngle - atan((fGermaniumWidth + 0.5*fGermaniumSeparation
																					 - fGermaniumBentLength*tan(fBentEndAngle))/
																					(fSuppressorBackRadius + fCanFaceThickness + fGermaniumDistFromCanFace));

    // the difference between the 22.5 angle and the angle that the heavymet makes to the germanium
    fBackBGOThickness                   = 3.7*cm;    // 6.0*cm;
    fBGOChoppedTip                      = 5.0*mm;
    fBGOMovableSpace                    = 3.0*cm;  	// amount the can face can move back into the BGO
    fSideSuppressorBackOffset           = 1.0*cm;
    fSideBGOThickness                   = 2.0*cm;
    fBGOCanSeperation                   = 0.5*mm;  	// Distance seperating the BGO and the detector can

    fSideBGOLength                      = fBGOMovableSpace + fDetectorTotalLength + fBackBGOThickness + fBGOCanSeperation;

    fSideSuppressorLength               = fDetectorTotalLength - fBentEndLength + fBGOCanSeperation + fBackBGOThickness
		- fSideSuppressorBackOffset - (fBGOCanSeperation + fBGOChoppedTip)/tan(fBentEndAngle);

    fBGOTrapLength                      = (fSideBGOThickness/tan(fBentEndAngle)) - (fBGOChoppedTip/tan(fBentEndAngle));

    fSuppressorExtensionThickness       = 1.0*cm;

    fSuppressorExtensionLength          = (fSuppressorBackRadius + fBentEndLength + (fBGOCanSeperation + fSideBGOThickness) / tan(fBentEndAngle)
														 - fSuppressorExtensionThickness * sin(fBentEndAngle)) / cos(fBentEndAngle) - fHevimetTipThickness - fSuppressorForwardRadius;

    fSuppressorExtensionLengthDet       = (fBackInnerRadius + fBentEndLength
														 + (fBGOCanSeperation + fSideBGOThickness) / tan(fBentEndAngle)
														 - fSuppressorExtensionThickness * sin(fBentEndAngle)) / cos(fBentEndAngle)
		- fHevimetTipThickness - fForwardInnerRadius;

    fSuppressorExtensionAngle 	       = atan(((fSuppressorBackRadius + fBentEndLength + (fBGOCanSeperation
																														 + fSideBGOThickness) / tan(fBentEndAngle)
																 - fSuppressorExtensionThickness *sin(fBentEndAngle))
																* tan(fBentEndAngle) - (fSuppressorForwardRadius + fHevimetTipThickness)
																* sin(fBentEndAngle)) / (fSuppressorExtensionLength));

    fHeavyMetThickness                  = 3.0*cm;   	// This can't be more than fBentEndLength + (choppedTip+seperation)/tan(fBentAngle)

    fHeavyMetInsideAngle                = atan((fDetectorTotalWidth/2.0)/((fDetectorTotalWidth/2.0 +fBGOChoppedTip
																									+ fBGOCanSeperation)/tan(fBentEndAngle)
																								  + fBGOMovableSpace +fBentEndLength));

    // used only in new suppressor design, but don't affect old suppressor design
    // As of September 2013, fSuppressorShellsIncludeFlag will always be true, so any dependence on it should be removed.

    fAirBoxFrontWidth 		             = 2.0*(fDetectorTotalWidth/2.0 +fBGOCanSeperation
															+ (fSideBGOThickness + fSuppressorShellThickness*2.0)
															+ ( fSuppressorExtensionLength
																 +  (fSuppressorShellThickness*2.0
																	  * (1.0/tan(fBentEndAngle)-tan(fBentEndAngle))) )
															* sin(fBentEndAngle) +(fSuppressorExtensionThickness
																						  + fSuppressorShellThickness*2.0)/cos(fBentEndAngle)) ;
    // the longer width of the box's trapezoid
    fAirBoxFrontWidthDet                = 2.0*(fDetectorTotalWidth/2.0 +fBGOCanSeperation
															  + (fSideBGOThickness + fSuppressorShellThickness*2.0)
															  + ( fSuppressorExtensionLengthDet
																	+  (fSuppressorShellThickness*2.0
																		 * (1.0/tan(fBentEndAngle)-tan(fBentEndAngle))) )
															  * sin(fBentEndAngle) +(fSuppressorExtensionThickness
																							 + fSuppressorShellThickness*2.0)/cos(fBentEndAngle)) ;

    fAirBoxFrontLength                  = (fAirBoxFrontWidth/2.0 - fSuppressorForwardRadius
														  * sin(fBentEndAngle)) / tan(fBentEndAngle);
	 
    fAirBoxFrontLengthDet               =  (fAirBoxFrontWidthDet/2.0 - fForwardInnerRadius
														  * sin(fBentEndAngle)) / tan(fBentEndAngle) ;

    fAirBoxBackLength                   = fDetectorTotalLength + fColdFingerShellLength
		+ fCoolantLength + ( fSuppressorBackRadius
									- fSuppressorForwardRadius * cos(fBentEndAngle))
		- fAirBoxFrontLength;
	 
    fAirBoxBackLengthDet                = fDetectorTotalLength + fColdFingerShellLength
		+ fCoolantLength + (fBackInnerRadius
								  - fForwardInnerRadius * cos(fBentEndAngle))
		- fAirBoxFrontLengthDet;

    // shift" is applied to everything to place them relative to the front of the airBox
    // the negative is because the shift is along the -ive X-xis

    fSuppShift                          = - (fAirBoxBackLength/2.0 + fAirBoxFrontLength
															- (fSuppressorForwardRadius * (1 - cos(fBentEndAngle)))
															- fCanFaceThickness/2.0);


    fShift                              = - (fAirBoxBackLengthDet/2.0 + fAirBoxFrontLengthDet
															- (fForwardInnerRadius - fForwardInnerRadius*cos(fBentEndAngle))
															- fCanFaceThickness/2.0); // Original


    // Gives the proper diameter of the Tigress array based on can face length
    fRhombiDiameter                     = fDetectorTotalWidth
		- 2.0*fBentEndLength
		* tan(fBentEndAngle)
		+ 2.0*(fDetectorTotalWidth
				 - 2.0*fBentEndLength
				 * tan(fBentEndAngle))
		* cos(M_PI/4.0);

    // This ones used by the new suppressor design
    fNewRhombiRadiusDet                 = fForwardInnerRadius*cos(fBentEndAngle); // Original
    fNewRhombiRadius                    = fSuppressorForwardRadius*cos(fBentEndAngle);


    if(canPositionSpecifier == 1) {
		fDetectorPositionShift              =  fBentEndLength +(fBGOChoppedTip + fBGOCanSeperation)/tan(fBentEndAngle);
    } else if (canPositionSpecifier == 2) {
		fDetectorPositionShift              = 0.0*cm;
    } else if (canPositionSpecifier == 3) {
		fDetectorPositionShift              = -(fBGOMovableSpace);
    } else {
		fDetectorPositionShift              = 0.0*cm;
    }

    // this is used to actually apply the necessary detector shift to all the pieces involved
    if(fSuppressorPositionSelector == 0) {
        fAppliedBackShift                   = 0.0*cm;
    } else if(fSuppressorPositionSelector == 1) {
        fAppliedBackShift                   = fBackInnerRadius - fForwardInnerRadius; // 3.5cm
    }

    if(detectorFormat == 1) {
		fGermaniumSelector                   = 0;
    } else if (detectorFormat == 2) {
		fGermaniumSelector                   = 1;
    } else if(detectorFormat == 3) {
		fGermaniumSelector                   = 2;
    } else {
		fGermaniumSelector                   = 2;
    }

    if (includeDetectorCan == "y") {
		fCanSelector                         = 1;
    } else if(includeDetectorCan == "n") {
		fCanSelector                         = 0;
    } else {
		fCanSelector                         = 0;
    }

    if(includeBGO == "y") {
		fBGOSelector                         = 1;
    } else if(includeBGO == "n") {
		fBGOSelector                         = 0;
    } else {
		fBGOSelector                         = 0;
    }

    if(includeColdFinger == "y") {
		fColdFingerSelector                 = 1;
    } else if(includeColdFinger == "n") {
		fColdFingerSelector                 = 0;
    } else {
		fColdFingerSelector                 = 0;
    }

    if (fSuppressorPositionSelector == 0) {
		fCrystalDistFromOrigin              = fForwardInnerRadius + fGermaniumDistFromCanFace + fCanFaceThickness;
    } else if (fSuppressorPositionSelector == 1) {
		fCrystalDistFromOrigin              = fBackInnerRadius    + fGermaniumDistFromCanFace + fCanFaceThickness;
    }

    fRadialDistance                     = fForwardInnerRadius*cos(fBentEndAngle) + fAirBoxBackLength/2.0 +fAirBoxFrontLength;

    //Redacted parameters/////////////////////////////////////////////////////////////
    fDetectorPlacementCxn = 0.4*mm;
    fTrianglePostDim = 1.0*micrometer;

    fSuppressorExtRightX = -2.55*mm;
    fSuppressorExtRightY = 0.24*mm;
    fSuppressorExtRightZ = -0.1*mm;

    fSuppressorExtLeftX = -2.55*mm;
    fSuppressorExtLeftY = 0.24*mm;
    fSuppressorExtLeftZ = 0.1*mm;

    fWedgeDim = 0.001*mm;

    fQuarterDetectorCxn = 0.01*mm;
    fQuarterDetectorCxnB = 0.1*mm;
	 
	 fElectrodeMaterial = "G4_Cu";
    fStructureMaterial = "Aluminum";

}// end ::DetectionSystemGriffin



///// Legacy 8Pi content /////////////////////////////////////////////////////////////////////
DetectionSystem8pi::DetectionSystem8pi() :
    // Logical Volumes
    fGermaniumBlockLog(0),
    fGermaniumDeadLayerLog(0),
    fGermaniumVacuumCoreLog(0),
    fLowerElectrodeMatElectrodeLog(0),
    fUpperElectrodeMatElectrodeLog(0),
    fInnerCage1Log(0),
    fInnerCage2Log(0),
    fInnerCage3Log(0),
    fInnerCage4Log(0),
    fInnerCageBottomLog(0),
    fInnerCageLidLog(0),
    fStructureMatCoolingRodLog(0),
    fElectrodeMatCoolingRodLog(0),
    fCoolingRodCoverLog(0),
    fCoolingRodCoverLidLog(0),
    fOuterCanSideLog(0),
    fOuterCanLidLog(0),
    fBerylliumWindowLog(0),
    fInnerBGOAnnulusLog(0),
    fStructureMatSheathLog(0),
    fOuterLowerBGOAnnulusLog(0),
    fOuterUpperBGOAnnulusLog(0),
    fLiquidN2Log(0),
    fLiquidN2SideLog(0),
    fLiquidN2LidLog(0),
    fLiquidN2BottomLog(0),
    fHevimetalLog(0),
    fAuxMatPlugLog(0),
    fAuxMatLayerLog(0)
{
    /*Set detail view angle (360 for full, 180 for half view)*/
    //myDetector->detailViewEndAngle = 180.0*deg; //cross-section view
    fDetailViewEndAngle = 360.0*deg; //fully enclosed

    fCutClearance = 0.1*mm;

    /*Set detector measurements for myDetector*/
    //~   denotes estimated from blueprints
    //*   denotes obtained/converted from blueprints
    //?   denotes unknown/estimated
    //OK  denotes confirmed measurements

    //distance from origin to front face of Hevimetal
    //  distFromOrigin = 9.906*cm; //OK
    fDistFromOrigin = 0.0*cm; //OK

    //Germanium crystal, core, dead layer & electrode
    fCrystalOuterRadius = 2.50*cm; //OK    *** Changed from spec (2.585*cm) to give realistic efficiencies *** -Evan
    fCrystalInnerRadius = 0.50*cm; //OK
    fHoleStartingDepth = 2.0*cm; //OK
    fCrystalLength = 5.0*cm; //OK     *** Changed from spec (5.620*cm) to give realistic efficiencies *** -Evan
    fDeadLayerThickness = 1.0*mm; //OK
    fElectrodeRadius = 1.0*mm; //OK

    //StructureMat cage surrounding Ge crystal
    fInnerCanThickness = 0.5*mm; //OK
    fInnerCanExtendsPastCrystal = 2.0*cm; //~2.0cm
    fInnerCanLidThickness = 1.78*mm; //*
    fInnerCanLidSeparation = 4.0*mm; //OK

    //Beryllium window
    fBerylliumDistFromCrystal = 2.0*mm; //~   //lower cage ring thickness omitted
    fBerylliumThickness = 0.5*mm; //OK

    //StructureMat can surrounding Ge crystal & cage
    fOuterCanLength = 11.5*cm; //OK
    fOuterCanThickness = 0.75*mm; //OK

    //StructureMat & ElectrodeMat Cooling Rods & Cover
    fStructureMatCoolingRodRadius = 3.0*mm; //~

    fElectrodeMatCoolingRodRadius = 4.76*mm; //OK   //must be > structureMat rod radius
    fElectrodeMatCoolingRodLength = 36.75*cm; //OK    //length from where structureMat rod ENDS to LN2 container

    //Inner BGO surrounding cooling rod
    fInnerBGOAnnulusLength = 8.1*cm; //OK from blueprint
    fInnerBGOClearance = 0.5*mm; //~ dist. b/w BGO & surrounding structureMat cover
    fInnerBGOInnerRadius = 9.9*mm; //OK    from blueprint
    fInnerBGOOuterRadius = 30.7*mm; //OK   from blueprint

    //Outer BGO surrounding entire detector
    fOuterBGOBottomThickness = 1.35*cm; //OK
    fOuterBGOTaperHeight = 6.2*cm; //*
    fOuterBGOTopOuterRadius = 6.55*cm;
    fOuterBGOTotalLength = 19.0*cm;//*
    fOuterBGOClearance = 0.5*mm;//~
    fOuterBGODisplacement = 3.25*cm; //OK (distance that detector is recessed from BGO) *** Changed from spec (1.0*cm) to give realistic efficiencies *** -Evan

    //Liquid N2 container measurements
    fLiquidN2Length = 22.8*cm; //OK
    fLiquidN2Radius = 7.64*cm;  //OK

    //Hevimetal collimator measurements
    fHevimetalThickness = 2.54*cm; //OK
    fHevimetalFrontSideLength = 4.369*cm; //OK
    fHevimetalRearSideLength = 5.489*cm; //OK

    //AuxMat
    fAuxMatPlugNegRadius = 1.64*cm; //OK
    fAuxMatPlugPosRadius = 2.06*cm; //OK
    fAuxMatLayerThickness = 9.6*mm; //OK
    //fAuxMatLayerThickness = 1.50*cm; //OK

    //suppressed:
    fExtraClearance = 0.5*mm;
    fElectrodeMat = "Copper";
    fStructureMat = "Aluminum";
    fAuxMat = "Delrin";

    /*****************************************************************************/
    /*** definitions for commonly used parameters - DO NOT EDIT (unless you really have to) ***/
    /**/  //distance to center of Ge crystal from world origin
    /**/  fCrystalDistFromOrigin = fCrystalLength/2.0
            /**/                 + fBerylliumDistFromCrystal
            /**/                 + fBerylliumThickness
            /**/                 + fOuterBGODisplacement
            /**/                 + fHevimetalThickness
            /**/                 + fDistFromOrigin;
    /**/  //distance from top of crystal to top of outer can
    /**/  fOuterCanExtendsPastCrystal =  fOuterCanLength
            /**/                    - fCrystalLength
            /**/                    - fBerylliumDistFromCrystal
            /**/                    - fBerylliumThickness;
    /**/  //inner radius of outer structureMat can
    /**/  fOuterCanInnerRadius =   fInnerBGOOuterRadius
            /**/                + fInnerBGOClearance
            /**/                - fOuterCanThickness;
    /**/  //auxMat front side length calculated as extension from hevimetal
    /**/  fAuxMatLayerFrontSideLength =  fHevimetalFrontSideLength
            /**/                    - (fAuxMatLayerThickness
                                       /**/                    * (fHevimetalRearSideLength - fHevimetalFrontSideLength)
                                       /**/                    / fHevimetalThickness);
    /**************************************************************/
}

///////////////////////////////////////////////////////////////////////
// The ::Apparatus8piVacuumChamber constructor initiates all the Logical
// and Physical Volumes used in the vacuum Chamber geometery, and the
//  ::Apparatus8piVacuumChamber  destructor deletes them from the stack
// when they go out of scope
///////////////////////////////////////////////////////////////////////
Apparatus8piVacuumChamber::Apparatus8piVacuumChamber() :
    // LogicalVolumes
    fVacuumChamberSphereLog(0)
{
    /////////////////////////////////////////////////////////////////////
    // Apparatus8piVacuumChamber Physical Properties
    /////////////////////////////////////////////////////////////////////

    fVacuumMaterial                               = "Vacuum";
    fVacuumChamberSphereMaterial                  = "Delrin";
    fVacuumChamberOuterRadius                     = 89.4*mm;

}// end ::Apparatus8piVacuumChamber

///////////////////////////////////////////////////////////////////////
// The ::Apparatus8piVacuumChamberAuxMatShell constructor initiates all the Logical
// and Physical Volumes used in the vacuum Chamber geometery, and the
//  ::Apparatus8piVacuumChamberAuxMatShell  destructor deletes them from the stack
// when they go out of scope
///////////////////////////////////////////////////////////////////////
Apparatus8piVacuumChamberAuxMatShell::Apparatus8piVacuumChamberAuxMatShell() :
    // LogicalVolumes
    fVacuumChamberAuxSphereLog(0)
{
    /////////////////////////////////////////////////////////////////////
    // Apparatus8piVacuumChamberAuxMatShell Physical Properties
    /////////////////////////////////////////////////////////////////////

    fVacuumChamberSphereMaterial                  = "Delrin";
    fVacuumChamberOuterRadius                     = 89.4*mm;

}// end ::Apparatus8piVacuumChamberAuxMatShell



void DetectorConstruction::AddDetectionSystem8pi(G4int ndet)
{
    // Describe Placement
    G4double theta,phi,position;
    G4ThreeVector move,direction;

    G4double hexagonAlign[20][5] = {
        {-0.982247, 0.000000,-0.187593,  100.812,  180.000},
        {-0.303531, 0.934172,-0.187592,  100.812,  108.000},
        { 0.794654, 0.577351,-0.187592,  100.812,   36.000},
        { 0.794655,-0.577350,-0.187592,  100.812,  324.000},
        {-0.303531,-0.934172,-0.187592,  100.812,  252.000},
        { 0.982247, 0.000000, 0.187593,   79.188,    0.000},
        { 0.303531,-0.934172, 0.187592,   79.188,  288.000},
        {-0.794654,-0.577351, 0.187592,   79.188,  216.000},
        {-0.794655, 0.577350, 0.187592,   79.188,  144.000},
        { 0.303531, 0.934172, 0.187592,   79.188,   72.000},
        {-0.607062, 0.000000,-0.794655,  142.623,  180.000},
        {-0.187592, 0.577350,-0.794655,  142.623,  108.000},
        { 0.491124, 0.356822,-0.794654,  142.623,   36.000},
        { 0.491124,-0.356822,-0.794654,  142.623,  324.000},
        {-0.187592,-0.577350,-0.794655,  142.623,  252.000},
        { 0.607062, 0.000000, 0.794655,   37.377,    0.000},
        { 0.187592,-0.577350, 0.794655,   37.377,  288.000},
        {-0.491124,-0.356822, 0.794654,   37.377,  216.000},
        {-0.491124, 0.356822, 0.794654,   37.377,  144.000},
        { 0.187592, 0.577350, 0.794655,   37.377,   72.000}};

    DetectionSystem8pi* pDetectionSystem8pi = new DetectionSystem8pi() ;
    pDetectionSystem8pi->Build() ;

    for(G4int detectorNumber = 0; detectorNumber < ndet; detectorNumber++)
    {
        theta = hexagonAlign[detectorNumber][3]*deg;
        phi = hexagonAlign[detectorNumber][4]*deg;

        G4RotationMatrix* rotate = new G4RotationMatrix; // rotation matrix corresponding to direction vector
        rotate->rotateZ(30.0*deg);
        rotate->rotateY(theta);
        rotate->rotateZ(phi);

        direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        position = 9.906*cm;
        move = position * direction;

        pDetectionSystem8pi->PlaceDetector(fLogicWorld, move, rotate, detectorNumber);
    }

}

void DetectorConstruction::AddDetectionSystem8piDetector(G4int ndet)
{
    // Describe Placement
    G4double theta,phi,position;
    G4ThreeVector move,direction;

    G4double hexagonAlign[20][5] = {
        {-0.982247, 0.000000,-0.187593,  100.812,  180.000},
        {-0.303531, 0.934172,-0.187592,  100.812,  108.000},
        { 0.794654, 0.577351,-0.187592,  100.812,   36.000},
        { 0.794655,-0.577350,-0.187592,  100.812,  324.000},
        {-0.303531,-0.934172,-0.187592,  100.812,  252.000},
        { 0.982247, 0.000000, 0.187593,   79.188,    0.000},
        { 0.303531,-0.934172, 0.187592,   79.188,  288.000},
        {-0.794654,-0.577351, 0.187592,   79.188,  216.000},
        {-0.794655, 0.577350, 0.187592,   79.188,  144.000},
        { 0.303531, 0.934172, 0.187592,   79.188,   72.000},
        {-0.607062, 0.000000,-0.794655,  142.623,  180.000},
        {-0.187592, 0.577350,-0.794655,  142.623,  108.000},
        { 0.491124, 0.356822,-0.794654,  142.623,   36.000},
        { 0.491124,-0.356822,-0.794654,  142.623,  324.000},
        {-0.187592,-0.577350,-0.794655,  142.623,  252.000},
        { 0.607062, 0.000000, 0.794655,   37.377,    0.000},
        { 0.187592,-0.577350, 0.794655,   37.377,  288.000},
        {-0.491124,-0.356822, 0.794654,   37.377,  216.000},
        {-0.491124, 0.356822, 0.794654,   37.377,  144.000},
        { 0.187592, 0.577350, 0.794655,   37.377,   72.000}};

    DetectionSystem8pi* pDetectionSystem8pi = new DetectionSystem8pi() ;
    pDetectionSystem8pi->Build() ;

    theta = hexagonAlign[ndet-1][3]*deg;
    phi = hexagonAlign[ndet-1][4]*deg;

    G4RotationMatrix* rotate = new G4RotationMatrix; //rotation matrix corresponding to direction vector
    rotate->rotateZ(30.0*deg);
    rotate->rotateY(theta);
    rotate->rotateZ(phi);

    direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    position = 9.906*cm;
    move = position * direction;

    pDetectionSystem8pi->PlaceDetector(fLogicWorld, move, rotate, ndet-1);
}
