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

#ifndef APPARATUSSPICETARGETCHAMBER_HH
#define APPARATUSSPICETARGETCHAMBER_HH

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

// define custom colours for visualisation
#define SN_COL 1.0, 1.0, 1.0
#define CU_COL 1.0, 1.0, 0.0
#define PEEK_COL 0.5, 0.5, 0.0
#define KAPTON_COL 0.2, 0.7, 0.1
#define PB_COL 0.6, 0.1, 0.1
#define NDFEB_COL 0.7,0.3,0.3
#define DELRIN_COL 0.0, 0.0, 1.0
#ifndef AL_COL
#define AL_COL 0.5, 0.5, 0.5
#endif

#define TARGET_CHAMBER_COL 0.15,0.15,0.15
#define ELECTROBOX_COL 0.15,0.15,0.15

class ApparatusSpiceTargetChamber
{
public:
	ApparatusSpiceTargetChamber(G4String, G4double);
	~ApparatusSpiceTargetChamber();

public:
	void Build(G4LogicalVolume*);

private:
	////////////////////////////////////////////
	// Logical Volumes used in ApparatusSpiceTargetChamber
	////////////////////////////////////////////
	G4LogicalVolume* fTargetChamberFrontRingLog;
	G4LogicalVolume* fTargetChamberFrontConeLog;
	G4LogicalVolume* fTargetChamberSphereLog;
	G4LogicalVolume* fTargetChamberCylinderDownLog;
	G4LogicalVolume* fTargetWheelLog;
	G4LogicalVolume* fFrameLogical;
	G4LogicalVolume* fCollimLogical;
	G4LogicalVolume* fExtLogical;
	G4LogicalVolume* fRodLogical;
	G4LogicalVolume* fFirstGearLog;
	G4LogicalVolume* fSecondGearLog;
	G4LogicalVolume* fThirdGearLog;
	G4LogicalVolume* fGearPlateOneLog;
	G4LogicalVolume* fGearPlateTwoLog;
	G4LogicalVolume* fGearStickLog;
	G4LogicalVolume* fTargetMountPlateLog;
	G4LogicalVolume* fBiasPlateLog;
	G4LogicalVolume* fPhotonShieldLayerOneLog;
	G4LogicalVolume* fPhotonShieldLayerTwoLog;
	G4LogicalVolume* fPhotonShieldLayerThreeLog;
	G4LogicalVolume* fPsTargetClampLog;
	G4LogicalVolume* fPsDetectorClampLog;
	G4LogicalVolume* fTargetBoltLog;
	G4LogicalVolume* fDetectorBoltLog;
	G4LogicalVolume* fMagnetLog;
	G4LogicalVolume* fMagnetCoverLog;
	G4LogicalVolume* fMclampChamberLog;
	G4LogicalVolume* fElectroBoxLog;
	G4LogicalVolume* fShieldCoverLog;
	G4LogicalVolume* fColdFingerLog;
	G4LogicalVolume* fS3CaseLogical;
	G4LogicalVolume* fConicalCollimatorLog;//11/8
	G4LogicalVolume* fXRayInsertLog;
	
private:
	////////////////////////////////////////////
	// Physical Volumes used in ApparatusSpiceTargetChamber
	////////////////////////////////////////////
	G4VPhysicalVolume* fTargetChamberFrontRingPhys;
	G4VPhysicalVolume* fTargetChamberFrontConePhys;
	G4VPhysicalVolume* fTargetChamberSpherePhys;
	G4VPhysicalVolume* fTargetChamberCylinderDownPhys;
	G4VPhysicalVolume* fTargetWheelPhys;
	G4VPhysicalVolume* fFramePhysical;
	G4VPhysicalVolume* fCollimPhysical;
	G4VPhysicalVolume* fExtPhysical;
	G4VPhysicalVolume* fRodPhysical;
	G4VPhysicalVolume* fFirstGearPhys;
	G4VPhysicalVolume* fSecondGearPhys;
	G4VPhysicalVolume* fThirdGearPhys;
	G4VPhysicalVolume* fGearPlateOnePhys;
	G4VPhysicalVolume* fGearPlateTwoPhys;
	G4VPhysicalVolume* fGearStickPhys;
	G4VPhysicalVolume* fTargetMountPlatePhys;
	G4VPhysicalVolume* fBiasPlatePhys;
	G4VPhysicalVolume* fPhotonShieldLayerOnePhys;
	G4VPhysicalVolume* fPhotonShieldLayerTwoPhys;
	G4VPhysicalVolume* fPhotonShieldLayerThreePhys;
	G4VPhysicalVolume* fPsTargetClampPhys;
	G4VPhysicalVolume* fPsDetectorClampPhys;
	G4VPhysicalVolume* fTargetBoltPhys;
	G4VPhysicalVolume* fDetectorBoltPhys;
	G4VPhysicalVolume* fMagnetPhys;
	G4VPhysicalVolume* fMagnetCoverPhys;
	G4VPhysicalVolume* fMclampChamberPhys;
	G4VPhysicalVolume* fMclampShieldPhys;
	G4VPhysicalVolume* fElectroBoxPhys;
	G4VPhysicalVolume* fShieldCoverPhys;
	G4VPhysicalVolume* fColdFingerPhys;
	G4VPhysicalVolume* fS3CasePhysical;
	G4VPhysicalVolume* fConicalCollimatorPhys;//11/8
	G4VPhysicalVolume* fXRayInsertPhys;

private:
	////////////////////////////////////////////
	// Properties used in ApparatusSpiceTargetChamber
	////////////////////////////////////////////

	//-------------------------
	// Magnet properties:
	//-------------------------
	G4int fNumberOfMagnets;
	G4int fNumberOfFrames;
	G4double fZOffset;

	//-------------------------
	// Materials:
	//-------------------------
	G4String fMagnetMaterial; //
	G4String fTargetChamberMaterial;//
	G4String fDownstreamConeMaterial;//
	G4String fMountLead; //doesn't exist
	G4String fPhotonShieldLayerOneMaterial;//
	G4String fPhotonShieldLayerTwoMaterial;//
	G4String fPhotonShieldLayerThreeMaterial;//
	G4String fPsClampMaterial;//
	G4String fMagnetClampMaterial;//
	G4String fTargetWheelMaterial;//
	G4String fGearStickMaterial;//
	G4String fBiasPlateMaterial;//
	G4String fTargetWheelGearMaterialA;//
	G4String fTargetWheelGearMaterialB;//
	G4String fGearPlateMaterial;//
	G4String fTargetMountPlateMaterial;//
	G4String fElectroBoxMaterial;//
	G4String fSmallBoltMaterial;//
	G4String fLargeBoltMaterial;//
	G4String fShieldCoverMaterial;//
	G4String fMagnetCoverMaterial;//
	G4String fColdFingerMaterial;//
	G4String fS3CableCaseMaterial;//
	G4String fConicalCollimatorMaterial;//11/8
        G4String fXRayInsertMaterial;//
	//-------------------------
	// Dimensions:
	//-------------------------

	//-----------------------------
	// Dimensions of Target Chamber
	//-----------------------------
	G4double fTargetChamberCylinderInnerRadius;
	G4double fTargetChamberCylinderOuterRadius;
	G4double fTargetChamberCylinderLength;
	G4double fTargetChamberDistanceFromTarget;
	G4double fTargetChamberCylinderLipInnerRadius;
	G4double fTargetChamberCylinderLipThickness;
	G4double fTargetChamberCylinderDetOuterRadius;
	G4double fTargetChamberCylinderDetThickness;
	G4double fTargetChamberCylinderHoleLength;
	G4double fTargetChamberCylinderHoleBreadth;
	G4double fTargetChamberCylinderHoleDistanceX;
	G4double fTargetChamberCylinderHoleDistanceY;

	G4double fTargetChamberSphereInnerRadius;
	G4double fTargetChamberSphereOuterRadius;
	G4double fTargetChamberSphereCentre;
	G4double fTargetChamberSphereCuttingOuterRadius;
	G4double fTargetChamberSphereCutOffPoint;
	G4double fFrontRingInnerRadius;
	G4double fFrontRingOuterRadius;
	G4double fFrontRingLength;
	G4double fFrontRingCylinderOuterRadius;
	G4double fFrontRingCylinderLength;

	G4double fConeFrontInnerRad;
	G4double fConeFrontOuterRad;
	G4double fConeLength;


	// --------------------------
	// Dimensions of Target Wheel
	// --------------------------
	G4double fTargetWheelRadius;
	G4double fTargetWheelThickness;
	G4double fTargetWheelOffset;
	G4double fTargetRadius;
	G4double fTargetOffset;
	G4double fCollimatorRadius;
	// Bias Plate
	G4double fBiasPlateOuterRadius;
	G4double fBiasPlateThickness;
	G4double fBiasPlateOffsetZ;
	G4double fBiasPlateCut;
	// Gears
	G4double fFirstGearRadius;
	G4double fFirstGearPlaneOffset;
	G4double fSecondGearRadius;
	G4double fSecondGearPlaneOffset;
	G4double fThirdGearRadius;
	G4double fThirdGearPlaneOffset;
	G4double fAllGearThickness;
	G4double fAllGearZOffset;
	// Gear Plates
	G4double fGearPlateOneRadius;
	G4double fGearPlateTwoRadius;
	G4double fGearPlateThickness;
	// Mount Plate
	G4double fTargetMountPlateRadius;
	G4double fTargetMountPlateThickness;
	G4double fTargetMountPlateZOffset;
	// Gear Stick
	G4double fGearStickLength;
	G4double fGearStickRadius;
	G4double fGearStickZOffset;
	// Target Frame
	G4double fDimTargetFrameOutZ;
	G4double fDimTargetFrameOutY;
	G4double fDimTargetFrameCutY;

	//----------------------------
	// Dimensions of Photon Shield
	//----------------------------
	G4double fPhotonShieldFrontRadius;
	G4double fPhotonShieldBackRadius;
	G4double fPhotonShieldInnerRadius;
	G4double fPhotonShieldLength;
	G4double fPhotonShieldBackFacePos;

	G4double fPhotonShieldLayerOneThickness;
	G4double fPhotonShieldLayerTwoThickness;
	G4double fPhotonShieldLayerThreeThickness;

	// ----------------------------------
	// Dimensions of Photon Shield Clamps
	// ----------------------------------
	G4double fPsTargetClampRadius;
	G4double fPsTargetClampThickness;
	G4double fPsTargetClampExtension;
	G4double fPsDetClampRadius;
	G4double fPsDetClampThickness;

	// ---------------------------------
	// Dimensions of Photon Shield Bolts
	// ---------------------------------
	G4double fPsTargetBoltLength;
	G4double fPsTargetBoltRadius;
	G4double fPsTargetBoltPlaneOffset;
	G4double fPsTargetBoltHeadRadius;
	G4double fPsTargetBoltHeadLength;
	G4double fPsDetBoltLength;
	G4double fPsDetBoltRadius;
	G4double fPsDetBoltPlaneOffset;
	G4double fPsDetBoltHeadRadius;
	G4double fPsDetBoltHeadLength;

	//---------------------
	// Dimensions of Magnet
	//---------------------
	G4double fPlateOneThickness;
	G4double fPlateOneLength;
	G4double fPlateOneHeight;
	G4double fPlateOneLowerHeight;

	G4double fPlateOneEdgeX;
	G4double fPlateTwoThickness;
	G4double fPlateTwoLength;
	G4double fPlateTwoHeight;

	G4double fDistanceFromTarget;
	G4double fCuttingBoxAngle;

	// ------------------------------------
	// Dimensions of Magnet Clamp (Chamber)
	// ------------------------------------
	G4double fMclampChamberThickness;
	G4double fMclampChamberLength;
	G4double fMclampChamberHeight;

	// ------------------------------------------
	// Dimensions of Magnet Clamp (Photon Shield)
	// ------------------------------------------
	G4double fPsClampPsbuffer;
	G4double fPsClampChamberBuffer;
	G4double fPsClampThicknessBuffer;
	G4double fPsClampLength;

	// -------------------------
	// Dimensions of ElectroBox
	// -------------------------
	G4double fElectroboxOuterBoxLength;
	G4double fElectroboxInnerBoxLength;
	G4double fElectroboxMidLength;
	G4double fElectroboxLipLength;
	G4double fElectroboxLipInnerRadius;
	G4double fElectroboxZOffset;

	// -------------------------
	// Dimensions of ColdFinger
	// -------------------------
	G4double fColdfingerThickness;
	G4double fColdfingerLength;
	G4double fColdfingerWidth;
	G4double fColdfingerZOffset;
	G4double fColdfingerHoleRadius;

	// ------------------------------
	// Dimensions of Plastic Coatings
	// ------------------------------
	G4double fMagnetCoatingThickness;
	G4double fShieldCoatingThickness;

	// ---------------------------
	// individual offsets for visualisation
	// ---------------------------
	G4double fFrontDomeOffset;
	//G4double fTargetWheelOffset;
	G4double fMiddleRingOffset;
	G4double fBackDetAndAlBoxOffset;

	// -----------------------
	// Dimensions of Beam Pipe
	// -----------------------
	G4double fPipeInnerRadius;
	G4double fPipeOuterRadius;
	G4double fPipeZLength;
	G4double fPipeZOffset;
	
	// --------------------------------
	// Dimensions of Conical Collimator
	// --------------------------------
	G4double fMidInnerRadius;
	G4double fMidOuterRadius;
	G4double fMidZLength;
	G4double fOuterCylOuterRadius;
	G4double fOuterCylInnerRadius;
	G4double fOuterCylZLength;
	G4double fInnerCylOuterRadius;
	G4double fInnerCylInnerRadius;
	G4double fInnerCylZLength;
	G4double fEdgeInnerRadius;
	G4double fEdgeOuterRadius;
	G4double fEdgeZLength;
	
	// --------------------------
	// Dimensions of X-ray Insert
	// --------------------------
	G4double fInsertCylOuterRadius;
	G4double fInsertCylInnerRadius;
	G4double fInsertCylZLength;
	G4double fInsertEdgeInnerRadius;
	G4double fInsertEdgeOuterRadius;
	G4double fInsertEdgeZLength;
	G4double fInsertHoleRadius;
	G4double fInsertHoleLength;
	
	//-----------------------------
	// copy numbers
	//-----------------------------
	G4int fTargetChamberDownstreamCopyNumber;
	G4int fPhotonShieldCopyNumber;
	G4int fPhotonShieldClampCopyNumber;
	G4int fMagnetsCopyNumber;
	G4int fMagnetClampCopyNumber;
	G4int fTargetWheelCopyNumber;
	G4int fBiasPlateCopyNumber;
	G4int fGearCopyNumber;
	G4int fGearStickCopyNumber;
	G4int fElectroBoxCopyNumber;
	G4int fPhotonShieldClampBoltCopyNumber;
	G4int fShieldCoveringCopyNumber;
	G4int fMagnetCoveringCopyNumber;
	G4int fColdFingerCopyNumber;

private:
	//////////////////////////////////////////////////////
	// internal methods and functions in ApparatusSpiceTargetChamber::Build()
	//////////////////////////////////////////////////////

	// methods
	void BuildTargetChamberFrontRing();
	void BuildTargetChamberSphere();
	void BuildTargetChamberCylinderDownstream();
	void BuildTargetWheel();
	void BuildTargetWheelSimple();
	void BuildTargetWheelRods();
	void BuildTargetWheelGears();
	void BuildTargetWheelGearPlates();
	void BuildTargetFrame();
	void BuildTargetWheelExtension();
	void BuildCollimator();
	void BuildGearStick();
	void BuildBiasPlate();
	void BuildTargetMountPlate();
	void BuildPhotonShield();
	void BuildPhotonShieldClamps();
	void BuildPhotonShieldClampBolts();
	void BuildCollectorMagnet();
	void BuildMagnetCovering();
	void BuildMagnetClampChamber();
	void BuildElectroBox();
	void BuildShieldCovering();
	void BuildColdFinger();
	void BuildS3CableHolder();
	void BuildConicalCollimator();
	void BuildXrayInsert();
	
	void PlaceTargetChamberFrontRing(G4LogicalVolume*);
	void PlaceTargetChamberSphere(G4LogicalVolume*);
	void PlaceTargetChamberCylinderDownstream(G4LogicalVolume*);
	void PlaceTargetWheel(G4LogicalVolume*);
	void PlaceTargetWheelRods(G4double, G4double, G4LogicalVolume*);
	void PlaceTargetWheelGears(G4LogicalVolume*);
	void PlaceTargetWheelGearPlates(G4LogicalVolume*);
	void PlaceTargetFrame(G4double, G4LogicalVolume*);
	void PlaceCollimator(G4LogicalVolume*);
	void PlaceGearStick(G4LogicalVolume*);
	void PlaceTargetMountPlate(G4LogicalVolume*);
	void PlaceTargetWheelExtension(G4LogicalVolume*);
	void PlaceBiasPlate(G4LogicalVolume*);
	void PlacePhotonShield(G4LogicalVolume*);
	void PlacePhotonShieldClamps(G4LogicalVolume*);
	void PlacePhotonShieldClampBolts(G4int, G4LogicalVolume*);
	void PlaceCollectorMagnet(G4int, G4LogicalVolume*);
	void PlaceMagnetCovering(G4int, G4LogicalVolume*);
	void PlaceMagnetClampChamber(G4int, G4LogicalVolume*);
	void PlaceElectroBox(G4LogicalVolume*);
	void PlaceShieldCovering(G4LogicalVolume*);
	void PlaceColdFinger(G4LogicalVolume*);
	void PlaceS3CableHolder(G4LogicalVolume*);
        void PlaceConicalCollimator(G4LogicalVolume*);
	void PlaceXrayInsert(G4LogicalVolume*);
	
	G4double fTargetZ;
	
	// functions
	G4RotationMatrix* RotateMagnets(G4int);
	G4ThreeVector TranslateMagnets(G4int, G4double, G4double);
};
#endif