#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "PrimaryGeneratorMessenger.hh"//PGM
#include "PrimaryGeneratorAction.hh"//PGA, for PGM object

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4VSolid.hh"
#include "G4Trap.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <math.h>

#include "HistoManager.hh"
#include "ApparatusSpiceTargetChamber.hh"

////////////////////////////////////////////////
// building and placing the different parts   //
// of the SPICE target chamber                //
// including chamber itself, magnets and      //
// mounts for magnets and detectors           //
// not including detectors                    //
////////////////////////////////////////////////


//////////////////////////////////////////////////////
// define physical properties of ApparatusSpiceTargetChamber //
//////////////////////////////////////////////////////

ApparatusSpiceTargetChamber::ApparatusSpiceTargetChamber(G4String MedLo, G4double TargetPedestal)//parameter chooses which lens is in place.
{
	fTargetZ = TargetPedestal;
	G4cout << TargetPedestal << " <- Beam Position in constructor" << G4endl;

	fNumberOfMagnets = 4;
	fNumberOfFrames = 3;

	// Materials
	fMagnetMaterial = "NdFeB"; 
	fMagnetClampMaterial = "Aluminum";
	fTargetChamberMaterial = "Delrin"; 
	fPhotonShieldLayerOneMaterial = "WTa"; 
	fPhotonShieldLayerTwoMaterial = "Tin"; 
	fPhotonShieldLayerThreeMaterial = "Copper";
	fDownstreamConeMaterial = "Aluminum"; 
	fPsClampMaterial = "Aluminum";
	fTargetWheelMaterial = "Peek"; 
	fTargetWheelGearMaterialA = "Delrin";
	fTargetWheelGearMaterialB = "Aluminum";
	fGearPlateMaterial = "Aluminum";
	fTargetMountPlateMaterial = "Peek"; 
	fGearStickMaterial = "Peek"; 
	fElectroBoxMaterial = "Aluminum"; 
	fSmallBoltMaterial = "Brass"; 
	fLargeBoltMaterial = "Titanium";
	fShieldCoverMaterial = "Kapton";
	fMagnetCoverMaterial = "Peek";
	fColdFingerMaterial = "Copper";
	fS3CableCaseMaterial = "Delrin";
	fConicalCollimatorMaterial = "WTa";//11/8
	fXRayInsertMaterial="Copper";

	//-----------------------------
	// Dimensions of Target Chamber
	//-----------------------------
	fTargetChamberCylinderInnerRadius = 94.6*CLHEP::mm;
	fTargetChamberCylinderOuterRadius = 102*CLHEP::mm;
	fTargetChamberCylinderLength = 89*CLHEP::mm;
	fTargetChamberDistanceFromTarget = 1*CLHEP::mm;
	fTargetChamberCylinderLipInnerRadius = 48*CLHEP::mm;
	fTargetChamberCylinderLipThickness = 6*CLHEP::mm;
	fTargetChamberCylinderDetOuterRadius = 128*CLHEP::mm;
	fTargetChamberCylinderDetThickness = 10*CLHEP::mm;
	fTargetChamberCylinderHoleLength = 32*CLHEP::mm;
	fTargetChamberCylinderHoleBreadth = 17.7*CLHEP::mm;
	fTargetChamberCylinderHoleDistanceY = 31.39*CLHEP::mm;
	fTargetChamberCylinderHoleDistanceX = 74.64*CLHEP::mm;

	// Sphere
	fTargetChamberSphereInnerRadius = 97*CLHEP::mm;
	fTargetChamberSphereOuterRadius = 102*CLHEP::mm;
	fTargetChamberSphereCentre = 6*CLHEP::mm;
	fTargetChamberSphereCuttingOuterRadius = 21*CLHEP::mm;
	fTargetChamberSphereCutOffPoint = 92.8*CLHEP::mm;
	// Front Ring
	fFrontRingInnerRadius = 15*CLHEP::mm;
	fFrontRingOuterRadius = 34*CLHEP::mm;
	fFrontRingLength = 11.4*CLHEP::mm;
	fFrontRingCylinderOuterRadius = 19*CLHEP::mm;
	fFrontRingCylinderLength = 52.607*CLHEP::mm;

	fConeFrontInnerRad = 48.801*CLHEP::mm;
	fConeFrontOuterRad = 52.801*CLHEP::mm;
	fConeLength = 92.868*CLHEP::mm;

	// --------------------------
	// Dimensions of Target Wheel
	// --------------------------
	fTargetWheelRadius = 40*CLHEP::mm;
	fTargetWheelThickness = 1.5*CLHEP::mm; // 3.0 for old
	fTargetWheelOffset = 15*CLHEP::mm;
	fTargetRadius = 9*CLHEP::mm; //6*CLHEP::mm; // 4mm for old
	fTargetOffset = 15*CLHEP::mm;
	fCollimatorRadius = 9*CLHEP::mm;

	// Gears
	fFirstGearRadius = 12.7*CLHEP::mm;
	fFirstGearPlaneOffset = 80.25*CLHEP::mm;
	fSecondGearRadius = 15.875*CLHEP::mm;
	fSecondGearPlaneOffset = 51.675*CLHEP::mm;
	fThirdGearRadius = 50.8*CLHEP::mm;
	fThirdGearPlaneOffset = 15*CLHEP::mm;
	fAllGearThickness = 3.175*CLHEP::mm;
	fAllGearZOffset = 6.459*CLHEP::mm;

	// Gear Plates
	fGearPlateOneRadius = 5*CLHEP::mm;
	fGearPlateTwoRadius = 6*CLHEP::mm;
	fGearPlateThickness = 3.866*CLHEP::mm;

	// Mount Plate
	fTargetMountPlateRadius = 87*CLHEP::mm;
	fTargetMountPlateThickness = 5*CLHEP::mm;
	fTargetMountPlateZOffset = 13.5*CLHEP::mm;

	// Gear Stick
	fGearStickLength = 250*CLHEP::mm;
	fGearStickRadius = 3.18*CLHEP::mm;
	fGearStickZOffset = 5.072*CLHEP::mm;

	// Target Frame
	fDimTargetFrameOutZ = 0.5*CLHEP::mm;
	fDimTargetFrameOutY = 12.*CLHEP::mm;
	fDimTargetFrameCutY = 8.*CLHEP::mm;

	//----------------------------
	// Dimensions of Photon Shield
	//----------------------------
	fPhotonShieldFrontRadius = 11.65*CLHEP::mm; //12.64*CLHEP::mm;
	fPhotonShieldBackRadius = 23.19*CLHEP::mm; //24.108*CLHEP::mm;
	fPhotonShieldLength = 30.*CLHEP::mm;
	fPhotonShieldInnerRadius = 5.*CLHEP::mm;
	fPhotonShieldBackFacePos = -58.*CLHEP::mm;

	fPhotonShieldLayerOneThickness = 25*CLHEP::mm;
	fPhotonShieldLayerTwoThickness = 4*CLHEP::mm;
	fPhotonShieldLayerThreeThickness = 1*CLHEP::mm;

	// ----------------------------------
	// Dimensions of Photon Shield Clamps
	// ----------------------------------
	fPsTargetClampRadius = 11*CLHEP::mm;
	fPsTargetClampThickness = 2*CLHEP::mm;
	fPsTargetClampExtension = 8*CLHEP::mm;
	fPsDetClampRadius = 21.5*CLHEP::mm;
	fPsDetClampThickness = 3*CLHEP::mm;

	// ---------------------------------
	// Dimensions of Photon Shield Bolts
	// ---------------------------------
	fPsTargetBoltLength = 9.5*CLHEP::mm;
	fPsTargetBoltRadius = 1.5*CLHEP::mm;
	fPsTargetBoltPlaneOffset = 9*CLHEP::mm;
	fPsTargetBoltHeadRadius = 3.2*CLHEP::mm;
	fPsTargetBoltHeadLength = 1.2*CLHEP::mm;
	fPsDetBoltLength = 18.05*CLHEP::mm;
	fPsDetBoltRadius = 2.413*CLHEP::mm;
	fPsDetBoltPlaneOffset = 11*CLHEP::mm;
	fPsDetBoltHeadRadius = 4.854*CLHEP::mm;
	fPsDetBoltHeadLength = 3.2*CLHEP::mm;

	//---------------------
	// Dimensions of Magnet
	//---------------------
	fPlateOneThickness = 3*CLHEP::mm;
	fPlateOneLength = 75*CLHEP::mm; 
	fPlateOneHeight = 50*CLHEP::mm; // z-component
	fPlateOneLowerHeight = 20*CLHEP::mm;

	if(MedLo=="Med"||MedLo=="med"||MedLo=="MED"||MedLo=="MEL"||MedLo=="Medium"||MedLo=="medium"){
		fPlateTwoThickness = 3.4*CLHEP::mm; // LEL is 5.0, MEL is 3.4
		fPlateTwoLength = 55.0*CLHEP::mm; // LEL is 30.0, MEL is 55.0
	} else {
		fPlateTwoThickness = 5.0*CLHEP::mm; // LEL is 5.0, MEL is 3.4
		fPlateTwoLength = 30.0*CLHEP::mm; // LEL is 30.0, MEL is 55.0
	}
	fPlateTwoHeight = 50*CLHEP::mm;

	fCuttingBoxAngle = 60*CLHEP::deg;
	fDistanceFromTarget = 4*CLHEP::mm; // in z-direction

	fPlateOneEdgeX = (92*CLHEP::mm - fPlateOneLength); // 92mm is radius of chamber

	// ------------------------------------
	// Dimensions of Magnet Clamp (Chamber)
	// ------------------------------------
	fMclampChamberThickness = fPlateOneThickness 
		+ (2*fPlateTwoThickness) + (2*3.25*CLHEP::mm);
	fMclampChamberLength = 6.5*CLHEP::mm;
	fMclampChamberHeight = 83*CLHEP::mm;

	// ------------------------------------------
	// Dimensions of Magnet Clamp (Photon Shield)
	// ------------------------------------------
	fPsClampPsbuffer = 1*CLHEP::mm;
	fPsClampChamberBuffer = 5*CLHEP::mm;
	fPsClampThicknessBuffer = 1.25*CLHEP::mm; // on each side
	fPsClampLength = fPlateOneLength 
		- fPlateTwoLength + fPsClampPsbuffer 
		+ fPsClampChamberBuffer;

	// -------------------------
	// Dimensions of ElectroBox
	// -------------------------
	fElectroboxOuterBoxLength = 269*CLHEP::mm;
	fElectroboxInnerBoxLength = 243.5*CLHEP::mm;
	fElectroboxMidLength = 660*CLHEP::mm;
	fElectroboxLipLength = 21.5*CLHEP::mm;
	fElectroboxLipInnerRadius = 102*CLHEP::mm;
	fElectroboxZOffset = -90*CLHEP::mm;

	// ------------------------------
	// Dimensions of Plastic Coatings
	// ------------------------------
	fMagnetCoatingThickness = 1.25*CLHEP::mm;
	fShieldCoatingThickness = 0.2*CLHEP::mm;

	// -------------------------------------
	// individual offsets for visualisation
	// -------------------------------------
	fFrontDomeOffset = 0*CLHEP::mm;
	fMiddleRingOffset = 0*CLHEP::mm;
	fBackDetAndAlBoxOffset = 0*CLHEP::mm;

	// -------------------------
	// Dimensions of ColdFinger
	// -------------------------
	fColdfingerThickness = 5*CLHEP::mm ; // CHECK
	fColdfingerLength = 150*CLHEP::mm ; // CHECK
	fColdfingerWidth = 200*CLHEP::mm  ; // CHECK
	fColdfingerZOffset = -122*CLHEP::mm ; // CHECK
	fColdfingerHoleRadius = 10*CLHEP::mm ; // CHECK

	// -----------------------
	// Dimensions of Beam Pipe
	// -----------------------
	fPipeInnerRadius = 5.0*CLHEP::mm; // 5
	fPipeOuterRadius = 16.0*CLHEP::mm;
	fPipeZLength = 30.0*CLHEP::mm;
} // end ApparatusSpiceTargetChamber

//////////////////////////////////////////////////
// delete ApparatusSpiceTargetChamber
/////////////////////////////////////////////////

ApparatusSpiceTargetChamber::~ApparatusSpiceTargetChamber()
{
	// logical volumes in ApparatusSpiceTargetChamber
	delete fTargetChamberFrontRingLog;
	delete fTargetChamberFrontConeLog;
	delete fTargetChamberSphereLog;
	delete fTargetChamberCylinderDownLog;
	delete fTargetWheelLog;
	delete fFrameLogical;
	delete fCollimLogical;
	delete fFirstGearLog;
	delete fSecondGearLog;
	delete fThirdGearLog;
	delete fGearPlateOneLog;
	delete fGearPlateTwoLog;
	delete fGearStickLog;
	delete fTargetMountPlateLog;
	delete fPhotonShieldLayerOneLog;
	delete fPhotonShieldLayerTwoLog;
	delete fPhotonShieldLayerThreeLog;
	delete fPsTargetClampLog;
	delete fPsDetectorClampLog;
	delete fTargetBoltLog;
	delete fDetectorBoltLog;
	delete fMagnetLog;
	delete fMclampChamberLog;
	delete fElectroBoxLog;
	delete fShieldCoverLog;
	delete fMagnetCoverLog;
	delete fColdFingerLog;
	// physical volumes in ApparatusSpiceTargetChamber
	delete fTargetChamberFrontRingPhys;
	delete fTargetChamberFrontConePhys;
	delete fTargetChamberSpherePhys;
	delete fTargetChamberCylinderDownPhys;
	delete fTargetWheelPhys;
	delete fFramePhysical;
	delete fCollimPhysical;
	delete fFirstGearPhys;
	delete fSecondGearPhys;
	delete fThirdGearPhys;
	delete fGearPlateOnePhys;
	delete fGearPlateTwoPhys;
	delete fGearStickPhys;
	delete fTargetMountPlatePhys;
	delete fPhotonShieldLayerOnePhys;
	delete fPhotonShieldLayerTwoPhys;
	delete fPhotonShieldLayerThreePhys;
	delete fPsTargetClampPhys;
	delete fPsDetectorClampPhys;
	delete fTargetBoltPhys;
	delete fDetectorBoltPhys;
	delete fMagnetPhys;
	delete fMclampChamberPhys;
	delete fMclampShieldPhys;
	delete fElectroBoxPhys;
	delete fShieldCoverPhys;
	delete fMagnetCoverPhys;
	delete fColdFingerPhys;

} // end ~ApparatusSpiceTargetChamber

///////////////////////////////////////
// build and place components        //
///////////////////////////////////////

void ApparatusSpiceTargetChamber::Build(G4LogicalVolume* expHallLog)
{


	// 	BuildTargetChamberFrontRing();
	BuildTargetChamberSphere();
	BuildTargetChamberCylinderDownstream();
	BuildTargetWheel();
	BuildTargetWheelExtension();
	BuildTargetWheelRods();
	BuildTargetWheelGears();
	BuildTargetWheelGearPlates();
	BuildGearStick();
	BuildTargetMountPlate();
	BuildPhotonShield();
	BuildPhotonShieldClamps();
	BuildPhotonShieldClampBolts();
	BuildCollectorMagnet();
	BuildMagnetClampChamber();
	BuildElectroBox();
	BuildShieldCovering();
	BuildMagnetCovering();
	BuildColdFinger();  
	BuildS3CableHolder();
	BuildConicalCollimator();
	BuildXrayInsert();
	BuildCollimator();

	// 	PlaceTargetChamberFrontRing(expHallLog);
	PlaceTargetChamberSphere(expHallLog);
	PlaceTargetChamberCylinderDownstream(expHallLog); 
	PlaceTargetWheelGears(expHallLog); 
	PlaceTargetWheelGearPlates(expHallLog); 
	PlaceGearStick(expHallLog); 
	PlaceTargetMountPlate(expHallLog); 
	PlacePhotonShield(expHallLog); 
	PlaceShieldCovering(expHallLog);
	PlacePhotonShieldClamps(expHallLog);
	PlaceElectroBox(expHallLog);
	PlaceColdFinger(expHallLog); 
	PlaceS3CableHolder(expHallLog);
	PlaceConicalCollimator(expHallLog);
	PlaceXrayInsert(expHallLog);


	if(ApparatusSpiceTargetChamber::fTargetZ < -4.00){
		G4cout << "Building target pedestals" << G4endl;

		PlaceTargetWheel(expHallLog); 
		PlaceTargetWheelExtension(expHallLog);
		PlaceCollimator(expHallLog); 

		BuildTargetFrame();//First draft of target bool


		G4double fFrameDegrees[2] = {0*deg, -110*deg};
		for(G4int i=0; i<2; i++)
			PlaceTargetFrame(fFrameDegrees[i], expHallLog);

		G4double fRodDegrees[7] = {40*deg, -40*deg, -70*deg, -150*deg, 91.2*deg, 125.*deg, 158.8*deg};
		G4double fRodRadius[7] = {10.624, 10.624, 10.624, 10.624, 18.166, 8.4, 18.166 };
		for(G4int i=0; i<4; i++){
			PlaceTargetWheelRods(fRodDegrees[i], fRodRadius[i], expHallLog);
		}
	}else{
		//Th
		PlaceTargetWheel(expHallLog); 
		PlaceTargetWheelExtension(expHallLog);
		PlaceCollimator(expHallLog); 

		G4cout << "Target Pedestals not built due to beam position" << G4endl;
	}


	for(G4int copyID=0; copyID<fNumberOfMagnets; copyID++)    {
		PlaceCollectorMagnet(copyID, expHallLog);
		PlaceMagnetClampChamber(copyID, expHallLog);
		PlacePhotonShieldClampBolts(copyID, expHallLog);
		PlaceMagnetCovering(copyID, expHallLog);
	}

} // end Build

////////////////////////////////////////////////////
// methods used to build and place the components:
////////////////////////////////////////////////////

void ApparatusSpiceTargetChamber::BuildTargetWheelRods() {
	// Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(AL_COL));
	sVisAtt->SetVisibility(true);

	// Dimensions
	G4double sDimRadOuter = 1.0*CLHEP::mm;
	G4double sDimHalfThick = 3.75*CLHEP::mm;

	// Shapes
	G4Tubs* sRodTubs = new G4Tubs("sRodTubs", 0, sDimRadOuter, sDimHalfThick, 0, 360*CLHEP::deg);

	// Logical
	G4Material* sMaterial = G4Material::GetMaterial(fTargetWheelMaterial);
	fRodLogical = new G4LogicalVolume(sRodTubs, sMaterial, "TargetRodsLog", 0, 0, 0);
	fRodLogical->SetVisAttributes(sVisAtt);
}

void ApparatusSpiceTargetChamber::BuildTargetFrame() {
	// Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(0.0, 0.9, 0.1));
	sVisAtt->SetVisibility(true);	

	// Dimensions
	G4double sDimOuterRadius = 7.00*CLHEP::mm;
	G4double sDimInnerRadius = 4.00*CLHEP::mm;
	G4double sDimHalfThick = 0.25*CLHEP::mm;

	// Shapes
	G4Tubs* sTubs = new G4Tubs("sTubs", sDimInnerRadius, sDimOuterRadius, sDimHalfThick, 0, 360.*CLHEP::deg);
	G4Tubs* sLimit = new G4Tubs("sLimit", 20.648, 25.0, 2*sDimHalfThick, 0, 360.*CLHEP::deg);

	G4ThreeVector sTrans(15.*CLHEP::mm, 0., 0.);
	G4SubtractionSolid* sSub = new G4SubtractionSolid("sSub", sTubs, sLimit, 0, sTrans);

	G4Box* sBox = new G4Box("sBox", 2.0*CLHEP::mm, 1.0*CLHEP::mm, sDimHalfThick);
	sTrans.setY(9.0*sin(-45*CLHEP::deg));
	sTrans.setX(9.0*cos(-45*CLHEP::deg));

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(45.*CLHEP::deg);
	G4UnionSolid* sUnion = new G4UnionSolid("sUnion", sSub, sBox, sRotate, sTrans);

	sTrans.setY(9.0*sin(45*CLHEP::deg));
	sTrans.setX(9.0*cos(45*CLHEP::deg));

	sRotate->rotateZ(-90.*CLHEP::deg);
	G4UnionSolid* sUnion2 = new G4UnionSolid("sUnion", sUnion, sBox, sRotate, sTrans);


	// Logical
	G4Material* sMaterial = G4Material::GetMaterial("Aluminum");
	fFrameLogical = new G4LogicalVolume(sUnion2, sMaterial, "TargetFrame", 0, 0, 0);
	fFrameLogical->SetVisAttributes(sVisAtt);
}

void ApparatusSpiceTargetChamber::BuildCollimator() {
	// Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	sVisAtt->SetVisibility(true);	

	// Dimensions

	G4Box* sRectangle = new G4Box("sRectangle", 0.57*CLHEP::mm, 12.*CLHEP::mm, 0.5*CLHEP::mm);
	G4Trd* sTrd = new G4Trd("sTrd", 0.5*CLHEP::mm, 0.5*CLHEP::mm, 5.45*CLHEP::mm, 12.*CLHEP::mm, 3.7*CLHEP::mm);

	G4ThreeVector sTrans(-4.27, 0., 0.);	
	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateY(-90.*CLHEP::deg);
	G4UnionSolid* sUnion = new G4UnionSolid("sUnion", sRectangle, sTrd, sRotate, sTrans);

	G4Tubs* sTubs = new G4Tubs("sTubs", 15., 19.8, 0.5*CLHEP::mm, -37.3*CLHEP::deg, 74.6*CLHEP::deg);
	sTrans.setX(-15.2);
	G4UnionSolid* sUnion2 = new G4UnionSolid("sUnion2", sUnion, sTubs, 0, sTrans);


	G4Tubs* sTubs1 = new G4Tubs("sTubs1", 0., 1.0, 1.0*CLHEP::mm, 0., 360.*CLHEP::deg);
	sTrans.setX(-1);
	sTrans.setY(-5.13);
	G4SubtractionSolid* sSub = new G4SubtractionSolid("sSub", sUnion2, sTubs1, 0, sTrans);

	G4Tubs* sTubs5 = new G4Tubs("sTubs5", 0., 2.5, 1.0*CLHEP::mm, 0., 360.*CLHEP::deg);
	sTrans.setY(5.13);
	G4SubtractionSolid* sSub2 = new G4SubtractionSolid("sSub2", sSub, sTubs5, 0, sTrans);

	// Logical
	G4Material* CollimatorMaterial = G4Material::GetMaterial("Titanium"); 
	fCollimLogical = new G4LogicalVolume(sSub2, CollimatorMaterial, "collimator", 0, 0, 0);
	fCollimLogical->SetVisAttributes(sVisAtt);
}

void ApparatusSpiceTargetChamber::BuildTargetChamberFrontRing()
{
	// Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(AL_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	// Upstream Lip
	G4double TubeInnerRadius = fFrontRingInnerRadius;
	G4double TubeOuterRadius = fFrontRingOuterRadius;
	G4double TubeLength = fFrontRingLength/2.;
	// Middle Cylinder
	G4double CylinderOuterRadius = fFrontRingCylinderOuterRadius;
	G4double CylinderHalfLength = fFrontRingCylinderLength/2.;
	// Front Cone
	G4double ConeFrontInner = fConeFrontInnerRad;
	G4double ConeFrontOuter = fConeFrontOuterRad;
	G4double ConeHalfLength = fConeLength/2.;

	// ** ShapesF
	// Upstream Lip
	G4Tubs* UpstreamLip = new G4Tubs("UpstreamLip", TubeInnerRadius, 
			TubeOuterRadius, TubeLength,
			0 , 360*CLHEP::deg);
	// Middle Cylinder
	G4Tubs* MiddleCylinder = new G4Tubs("MiddleCylinder", TubeInnerRadius, 
			CylinderOuterRadius, CylinderHalfLength, 
			0, 360*CLHEP::deg);
	// Cone
	G4Cons* FrontCone = new G4Cons("FrontCone", TubeInnerRadius, 
			CylinderOuterRadius, ConeFrontInner, 
			ConeFrontOuter, ConeHalfLength, 0., 360*CLHEP::deg);

	G4ThreeVector trans(0, 0, TubeLength + CylinderHalfLength);
	G4UnionSolid* FrontRing = new G4UnionSolid("FrontRing", UpstreamLip, MiddleCylinder, 0, trans);

	// ** Logical Volume
	G4Material* frontRingMaterial = G4Material::GetMaterial(fDownstreamConeMaterial);
	fTargetChamberFrontRingLog = new G4LogicalVolume(FrontRing, frontRingMaterial,"TargetChamberFrontRingLog",0,0,0);
	fTargetChamberFrontRingLog->SetVisAttributes(VisAtt);
	fTargetChamberFrontConeLog = new G4LogicalVolume(FrontCone, frontRingMaterial,"TargetChamberFrontConeLog",0,0,0);
	fTargetChamberFrontConeLog->SetVisAttributes(VisAtt);
}//end BuildTargetChamberFrontRing

void ApparatusSpiceTargetChamber::BuildTargetChamberSphere()
{
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(TARGET_CHAMBER_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	// Sphere
	G4double SphereOuterRadius = fTargetChamberSphereOuterRadius;
	G4double SphereInnerRadius = fTargetChamberSphereInnerRadius;
	// Upstream Addon Cylinder
	G4double AddonCylinderLength = (fTargetChamberSphereCentre + fTargetChamberDistanceFromTarget)/2.;
	// Downstream Cutting Cylinder
	G4double CuttingOuterRadius = fTargetChamberSphereCuttingOuterRadius;
	// Downstream Cutting Box
	G4double CuttingBoxHalfLength = 50*CLHEP::mm;
	G4double CuttingBoxFromSphereCentre = fTargetChamberSphereCutOffPoint;

	// ** Shapes
	// Sphere
	G4Sphere* sphere = new G4Sphere("sphere", SphereInnerRadius, SphereOuterRadius, 0*CLHEP::deg, 360*CLHEP::deg, 0*CLHEP::deg, 90*CLHEP::deg);
	// Upstream Addon Cylinder
	G4Tubs* AddonCylinder = new G4Tubs("AddonCylinder", SphereInnerRadius, SphereOuterRadius, AddonCylinderLength, 0, 360*CLHEP::deg);	
	G4ThreeVector trans(0, 0, -AddonCylinderLength);
	G4UnionSolid* TargetChamberSpherePre = new G4UnionSolid("TargetChamberSpherePre", sphere, AddonCylinder, 0, trans);
	// Downstream Cutting Cylinder
	G4Tubs* CuttingCylinder = new G4Tubs("CuttingCylinder", 0, CuttingOuterRadius, 200*CLHEP::mm, 0, 360*CLHEP::deg);
	G4SubtractionSolid* TargetChamberSpherePre2 = new G4SubtractionSolid("TargetChamberSpherePre2", TargetChamberSpherePre, CuttingCylinder);
	// Downstream Cutting Box
	G4Box* CuttingBox = new G4Box("CuttingBox", CuttingBoxHalfLength, CuttingBoxHalfLength, CuttingBoxHalfLength);
	trans.setZ(CuttingBoxFromSphereCentre + CuttingBoxHalfLength);
	G4SubtractionSolid* TargetChamberSphere = new G4SubtractionSolid("TargetChamberSphere", TargetChamberSpherePre2, CuttingBox, 0, trans);

	G4Material* sphereMaterial = G4Material::GetMaterial(fTargetChamberMaterial);
	fTargetChamberSphereLog = new G4LogicalVolume(TargetChamberSphere, sphereMaterial,"TargetChamberSphereLog",0,0,0);
	fTargetChamberSphereLog->SetVisAttributes(VisAtt);
}

void ApparatusSpiceTargetChamber::BuildTargetChamberCylinderDownstream()
{
	// ** Visualization
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(TARGET_CHAMBER_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	// Main Cylinder
	G4double InnerRadius = fTargetChamberCylinderInnerRadius;
	G4double OuterRadius = fTargetChamberCylinderOuterRadius;
	G4double length = fTargetChamberCylinderLength/2.;
	// Cylinder Target Lip
	G4double TargetLipInnerRadius = fTargetChamberCylinderLipInnerRadius;
	G4double TargetLipLength = fTargetChamberCylinderLipThickness/2.;
	// Cylinder Detector Lip
	G4double DetectorLipOuterRadius = fTargetChamberCylinderDetOuterRadius;
	G4double DetectorLipLength = fTargetChamberCylinderDetThickness/2.;  
	// Target Lip Cutting Box
	G4double LipCuttingHalfX = fTargetChamberCylinderHoleLength/2.;
	G4double LipCuttingHalfY = fTargetChamberCylinderHoleBreadth/2.;
	G4double LipCuttingHalfZ = fTargetChamberCylinderLipThickness;

	// ** Shapes
	G4Tubs* MainCylinder = new G4Tubs("MainCylinder",InnerRadius,OuterRadius,length,0,360*CLHEP::deg);
	G4Tubs* TargetLipCylinder = new G4Tubs("TargetLipCylinder", TargetLipInnerRadius, OuterRadius, TargetLipLength, 0, 360*CLHEP::deg);
	G4Tubs* DetectorLipCylinder = new G4Tubs("DetectorLipCylinder", OuterRadius, DetectorLipOuterRadius, DetectorLipLength, 0, 360*CLHEP::deg);
	G4Box* TargetLipCuttingBox = new G4Box("TargetLipCuttingBox", LipCuttingHalfX, LipCuttingHalfY, LipCuttingHalfZ);

	// Union
	G4ThreeVector trans(0,0, length-TargetLipLength);
	G4UnionSolid* TargetDetCylinderPre = new G4UnionSolid("TargetDetCylinderPre", MainCylinder, TargetLipCylinder, 0, trans);
	trans.setZ(-(length-DetectorLipLength));
	G4UnionSolid* TargetDetCylinder = new G4UnionSolid("TargetDetCylinder", TargetDetCylinderPre, DetectorLipCylinder, 0, trans);
	trans.setX(-fTargetChamberCylinderHoleDistanceX);
	trans.setY(-fTargetChamberCylinderHoleDistanceY);
	trans.setZ(length-TargetLipLength);
	G4RotationMatrix* rotate = new G4RotationMatrix(111.8*CLHEP::deg , 0, 0);
	G4SubtractionSolid* TargetDetChamber = new G4SubtractionSolid("TargetDetChamber", TargetDetCylinder, TargetLipCuttingBox, rotate, trans);

	// ** Logical Volume
	G4Material* cylinderMaterial = G4Material::GetMaterial(fTargetChamberMaterial);
	fTargetChamberCylinderDownLog = new G4LogicalVolume(TargetDetChamber, cylinderMaterial,"TargetChamberCylinderDownLog",0,0,0);
	fTargetChamberCylinderDownLog->SetVisAttributes(VisAtt);
}

void ApparatusSpiceTargetChamber::BuildTargetWheel(){
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(AL_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	G4double WheelHalfThickness = fTargetWheelThickness/2.;

	// ** Shapes
	G4Tubs* TargetWheelPre = new G4Tubs("TargetWheelPre", 0, fTargetWheelRadius, WheelHalfThickness, 0, 360*CLHEP::deg);
	G4Tubs* target = new G4Tubs("target", 0, fTargetRadius, fTargetWheelThickness, 0, 360*CLHEP::deg);
	G4Box* collimator = new G4Box("collimator", 4.5*CLHEP::mm, 8.5*CLHEP::mm, fTargetWheelThickness);

	G4ThreeVector trans(0, 0, 0);

	trans.setX(fTargetOffset*cos(70*CLHEP::deg));
	trans.setY(fTargetOffset*sin(70*CLHEP::deg));
	G4SubtractionSolid* TargetWheel0 = new G4SubtractionSolid("TargetWheel0", TargetWheelPre, target, 0, trans);

	trans.setX(fTargetOffset*cos(180*CLHEP::deg));
	trans.setY(fTargetOffset*sin(180*CLHEP::deg));
	G4SubtractionSolid* TargetWheel1 = new G4SubtractionSolid("TargetWheel1", TargetWheel0, target, 0, trans);

	trans.setX(fTargetOffset*cos(-55*CLHEP::deg));
	trans.setY(fTargetOffset*sin(-55*CLHEP::deg));
	G4RotationMatrix* sRotate = new G4RotationMatrix(125*CLHEP::deg , 0, 0);
	G4SubtractionSolid* TargetWheel = new G4SubtractionSolid("TargetWheel", TargetWheel1, collimator, sRotate, trans);


	// ** Logical
	G4Material* targetWheelMaterial = G4Material::GetMaterial(fTargetWheelMaterial); //problem
	fTargetWheelLog = new G4LogicalVolume(TargetWheel, targetWheelMaterial, "TargetWheelLog", 0, 0, 0);
	fTargetWheelLog->SetVisAttributes(VisAtt);
} // end BuildTargetWheel()

void ApparatusSpiceTargetChamber::BuildTargetWheelExtension() {
	// Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(AL_COL));
	sVisAtt->SetVisibility(true);	

	// Dimensions
	G4double WheelRadius = fTargetWheelRadius;
	G4double InnerRadius = WheelRadius - 5.0*CLHEP::mm;
	G4double WheelHalfThickness = 5.*CLHEP::mm;

	// Shapes
	G4Tubs* extension = new G4Tubs("extension", InnerRadius, WheelRadius, WheelHalfThickness, 0, 360*CLHEP::deg);

	// Logical
	G4Material* targetWheelMaterial = G4Material::GetMaterial(fTargetWheelMaterial);
	fExtLogical = new G4LogicalVolume(extension, targetWheelMaterial, "fExtLogical", 0, 0, 0);
	fExtLogical->SetVisAttributes(sVisAtt);
}

void ApparatusSpiceTargetChamber::BuildTargetWheelGears(){
	// ** Visualisation
	G4VisAttributes* VisAttB = new G4VisAttributes(G4Colour(AL_COL));
	VisAttB->SetVisibility(true);
	G4VisAttributes* VisAttA = new G4VisAttributes(G4Colour(DELRIN_COL));
	VisAttA->SetVisibility(true);

	// ** Dimensions
	G4double FirstGearRadius = fFirstGearRadius;
	G4double SecondGearRadius = fSecondGearRadius;
	G4double ThirdGearOuterRadius = fThirdGearRadius;
	G4double ThirdGearInnerRadius = fTargetWheelRadius;
	G4double GearHalfThickness = fAllGearThickness/2.;

	// ** Shapes
	G4Tubs* FirstGear = new G4Tubs("FirstGear", 0, FirstGearRadius, GearHalfThickness, 0, 360*CLHEP::deg);
	G4Tubs* SecondGear = new G4Tubs("SecondGear", 0, SecondGearRadius, GearHalfThickness, 0, 360*CLHEP::deg);
	G4Tubs* ThirdGear = new G4Tubs("ThirdGear", ThirdGearInnerRadius, ThirdGearOuterRadius, GearHalfThickness, 0, 360*CLHEP::deg);

	// ** Logical
	G4Material* gearMaterialB = G4Material::GetMaterial(fTargetWheelGearMaterialB);
	fFirstGearLog = new G4LogicalVolume(FirstGear, gearMaterialB, "FirstGearLog", 0, 0, 0);
	G4Material* gearMaterialA = G4Material::GetMaterial(fTargetWheelGearMaterialA);
	fSecondGearLog = new G4LogicalVolume(SecondGear, gearMaterialA, "SecondGearLog", 0, 0, 0);
	fThirdGearLog = new G4LogicalVolume(ThirdGear, gearMaterialA, "ThirdGearLog", 0, 0, 0);

	fFirstGearLog->SetVisAttributes(VisAttB);
	fSecondGearLog->SetVisAttributes(VisAttA);
	fThirdGearLog->SetVisAttributes(VisAttA);
} // end BuildTargetWheelGears()

void ApparatusSpiceTargetChamber::BuildTargetWheelGearPlates() {
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(AL_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	G4double GearOneRadius = fGearPlateOneRadius;
	G4double GearTwoRadius = fGearPlateTwoRadius;
	G4double GearHalfThickness = fGearPlateThickness/2.;

	// ** Shapes
	G4Tubs* PlateOne = new G4Tubs("PlateOne", 0, GearOneRadius, GearHalfThickness, 0, 360*CLHEP::deg);
	G4Tubs* PlateTwo = new G4Tubs("PlateTwo", 0, GearTwoRadius, GearHalfThickness, 0, 360*CLHEP::deg);

	// ** Logical
	G4Material* gearPlateMaterial = G4Material::GetMaterial(fGearPlateMaterial);
	fGearPlateOneLog = new G4LogicalVolume(PlateOne, gearPlateMaterial, "GearPlateOneLog", 0, 0, 0);
	fGearPlateOneLog->SetVisAttributes(VisAtt);
	fGearPlateTwoLog = new G4LogicalVolume(PlateTwo, gearPlateMaterial, "GearPlateTwoLog", 0, 0, 0);
	fGearPlateTwoLog->SetVisAttributes(VisAtt);
} // end::BuildTargetWheelGearPlates()

void ApparatusSpiceTargetChamber::BuildGearStick(){
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(PEEK_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	G4double GearStickHalfLength = fGearStickLength/2.;

	// ** Shapes
	G4Tubs* GearStick = new G4Tubs("GearStick", 0, fGearStickRadius, GearStickHalfLength, 0, 360*CLHEP::deg);

	// ** Logical
	G4Material* gearStickMaterial = G4Material::GetMaterial(fGearStickMaterial);
	fGearStickLog = new G4LogicalVolume(GearStick, gearStickMaterial, "GearStickLog", 0, 0, 0);
	fGearStickLog->SetVisAttributes(VisAtt);
}

void ApparatusSpiceTargetChamber::BuildTargetMountPlate(){
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(PEEK_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	G4double MountPlateRadius = fTargetMountPlateRadius;
	G4double MountPlateHalfThickness = fTargetMountPlateThickness/2.;
	// Target Wheel Cut
	G4double WheelRadius = fTargetWheelRadius;

	// ** Shapes
	G4Tubs* MountPlatePre = new G4Tubs("MountPlatePre", 0, MountPlateRadius, MountPlateHalfThickness, 0, 360*CLHEP::deg);

	G4Tubs* TargetWheel = new G4Tubs("TargetWheel", 0, WheelRadius, 2.*MountPlateHalfThickness, 0, 360*CLHEP::deg);
	G4double offset = fTargetWheelOffset / sqrt(2.);
	G4ThreeVector move(offset, offset, 0);
	G4SubtractionSolid* MountPlate = new G4SubtractionSolid("MountPlate", MountPlatePre, TargetWheel, 0, move);

	// ** Logical
	G4Material* targetMountPlateMaterial = G4Material::GetMaterial(fTargetMountPlateMaterial);
	fTargetMountPlateLog = new G4LogicalVolume(MountPlate, targetMountPlateMaterial, "TargetMountPlateLog", 0, 0, 0);
	fTargetMountPlateLog->SetVisAttributes(VisAtt);
} // end BuildTargetMountPlate

void ApparatusSpiceTargetChamber::BuildPhotonShield()
{
	// ** Visualization
	G4VisAttributes* OneVisAtt = new G4VisAttributes(G4Colour(PB_COL));
	G4VisAttributes* TwoVisAtt = new G4VisAttributes(G4Colour(SN_COL));
	G4VisAttributes* ThreeVisAtt = new G4VisAttributes(G4Colour(CU_COL));
	OneVisAtt->SetVisibility(true);
	TwoVisAtt->SetVisibility(true);
	ThreeVisAtt->SetVisibility(true);

	// ** Build Photon Shield Solid (First Layer)
	// Define Variables
	G4double InnerRadius = fPhotonShieldInnerRadius;
	G4double TotalLength = fPhotonShieldLayerOneThickness + fPhotonShieldLayerTwoThickness + fPhotonShieldLayerThreeThickness;
	G4double HalfLength = TotalLength/2.;

	G4double BackOuterRadius = fPhotonShieldBackRadius;
	G4double FrontOuterRadius = fPhotonShieldFrontRadius;

	G4Cons* PhotonShield = new G4Cons("PhotonShield", InnerRadius, BackOuterRadius, 
			InnerRadius, FrontOuterRadius,	HalfLength, 0, 2*CLHEP::pi);

	// ** Build Magnet Clamp
	G4double ClampHalfThickness = fPlateOneThickness/2. + fPsClampThicknessBuffer;
	G4Box* MagnetClamp = new G4Box("MagnetClamp", fPsClampLength/2., ClampHalfThickness, fPhotonShieldLength);


	// ** Cut Plate One from Photon Shield
	G4RotationMatrix* cutRot = new G4RotationMatrix; 
	G4double transRadial = fPlateOneEdgeX + fPsClampLength/2. - fPsClampPsbuffer;
	G4double transX, transY;
	G4double transZ = 0;
	G4ThreeVector trans;

	cutRot->rotateZ(-(45.+0*90.)*CLHEP::deg);
	transX = transRadial*cos((0*90.+45.)*CLHEP::deg);
	transY = transRadial*sin((0*90.+45.)*CLHEP::deg);
	trans.set(transX,transY,transZ);
	G4SubtractionSolid* PhotonShieldTemp_1 = new G4SubtractionSolid("PhotonShield",PhotonShield,MagnetClamp,cutRot,trans);      
	cutRot->rotateZ(-90.*CLHEP::deg);
	transX = transRadial*cos((1*90.+45.)*CLHEP::deg);
	transY = transRadial*sin((1*90.+45.)*CLHEP::deg);
	trans.set(transX,transY,transZ);
	G4SubtractionSolid* PhotonShieldTemp_2 = new G4SubtractionSolid("PhotonShield",PhotonShieldTemp_1,MagnetClamp,cutRot,trans);      
	cutRot->rotateZ(-90.*CLHEP::deg);
	transX = transRadial*cos((2*90.+45.)*CLHEP::deg);
	transY = transRadial*sin((2*90.+45.)*CLHEP::deg);
	trans.set(transX,transY,transZ);
	G4SubtractionSolid* PhotonShieldTemp_3 = new G4SubtractionSolid("PhotonShield",PhotonShieldTemp_2,MagnetClamp,cutRot,trans);      
	cutRot->rotateZ(-90.*CLHEP::deg);
	transX = transRadial*cos((3*90.+45.)*CLHEP::deg);
	transY = transRadial*sin((3*90.+45.)*CLHEP::deg);
	trans.set(transX,transY,transZ);
	G4SubtractionSolid* PhotonShieldTemp_4 = new G4SubtractionSolid("PhotonShield",PhotonShieldTemp_3,MagnetClamp,cutRot,trans); 

	// ** Bolt Holes
	G4double TargetBoltHalfLength = fPsTargetBoltLength/2.;
	G4double TargetBoltRadius = fPsTargetBoltRadius;
	G4double TargetBoltOffset = fPsTargetBoltPlaneOffset / sqrt(2.);
	G4double TargetBoltZOffset = fPhotonShieldLength/2. - TargetBoltHalfLength;
	G4double DetectorBoltHalfLength = fPsDetBoltLength/2.;
	G4double DetectorBoltRadius = fPsDetBoltRadius;
	G4double DetectorBoltOffset = fPsDetBoltPlaneOffset / sqrt(2.);
	G4double DetectorBoltZOffset = -fPhotonShieldLength/2. + DetectorBoltHalfLength;

	G4Tubs* TargetBolt = new G4Tubs("TargetBolt", 0, TargetBoltRadius, TargetBoltHalfLength, 0, 360*CLHEP::deg);
	G4Tubs* DetectorBolt = new G4Tubs("DetectorBolt", 0, DetectorBoltRadius, DetectorBoltHalfLength, 0, 360*CLHEP::deg);

	// the four front (target end) bolts are created using a union solid then deducted from the previous volume
	G4double DistanceBetweenTargetBolts = fPsTargetBoltPlaneOffset * sqrt(2.);
	G4ThreeVector BoltDistance(DistanceBetweenTargetBolts, 0, 0);
	G4UnionSolid* TargetBoltsTwo = new G4UnionSolid("TargetBoltsTwo", TargetBolt, TargetBolt, 0, BoltDistance);
	G4ThreeVector BoltDist2(0, DistanceBetweenTargetBolts, 0);
	G4UnionSolid* TargetBolts = new G4UnionSolid("TargetBolts", TargetBoltsTwo, TargetBoltsTwo, 0, BoltDist2);
	G4ThreeVector BoltMove(-TargetBoltOffset, -TargetBoltOffset, TargetBoltZOffset);
	G4SubtractionSolid* PhotonShieldTemp_5 = new G4SubtractionSolid("PhotonShield", PhotonShieldTemp_4, TargetBolts, 0, BoltMove);

	G4double DistanceBetweenDetBolts = fPsDetBoltPlaneOffset * sqrt(2.);
	G4ThreeVector BoltDist3(DistanceBetweenDetBolts, 0, 0);
	G4UnionSolid* DetBoltsTwo = new G4UnionSolid("DetBoltsTwo", DetectorBolt, DetectorBolt, 0, BoltDist3); 
	G4ThreeVector BoltDist4(0, DistanceBetweenDetBolts, 0);
	G4UnionSolid* DetBolts = new G4UnionSolid("DetBolts", DetBoltsTwo, DetBoltsTwo, 0, BoltDist4);
	G4ThreeVector BoltMove2(-DetectorBoltOffset, -DetectorBoltOffset, DetectorBoltZOffset);
	G4SubtractionSolid* PhotonShieldTemp_6 = new G4SubtractionSolid("PhotonShield", PhotonShieldTemp_5, DetBolts, 0, BoltMove2);


	// PhotonShieldTemp_12 can be considered a template.  To segment the shield into layers I will create a cutting box
	G4Box* PhotonShieldSlicer = new G4Box("PhotonShieldSlicer",50*CLHEP::mm, 50*CLHEP::mm, 50*CLHEP::mm);

	// Cut out the first layer
	G4ThreeVector slicer(0, 0, -(50+HalfLength-(fPhotonShieldLayerThreeThickness+fPhotonShieldLayerTwoThickness)));
	G4SubtractionSolid* FirstLayer = new G4SubtractionSolid("FirstLayer", PhotonShieldTemp_6, PhotonShieldSlicer, 0, slicer);

	// Cut out a second layer
	slicer.set(0,0, (50+HalfLength-fPhotonShieldLayerOneThickness));
	G4SubtractionSolid* SecondPre = new G4SubtractionSolid("SecondPre", PhotonShieldTemp_6, PhotonShieldSlicer, 0, slicer);
	slicer.set(0, 0, -(50+HalfLength-fPhotonShieldLayerThreeThickness));
	G4SubtractionSolid* SecondLayer = new G4SubtractionSolid("SecondLayer", SecondPre, PhotonShieldSlicer, 0, slicer);

	// Cut out a third layer
	slicer.set(0, 0, (50+HalfLength-(fPhotonShieldLayerOneThickness+fPhotonShieldLayerTwoThickness)));
	G4SubtractionSolid* ThirdLayer = new G4SubtractionSolid("ThirdLayer", PhotonShieldTemp_6, PhotonShieldSlicer, 0, slicer);

	// Create the three separate logs 
	G4Material* layerOneMaterial = G4Material::GetMaterial(fPhotonShieldLayerOneMaterial);
	fPhotonShieldLayerOneLog = new G4LogicalVolume(FirstLayer, layerOneMaterial,"PhotonShieldLayerOneLog",0,0,0);
	fPhotonShieldLayerOneLog->SetVisAttributes(OneVisAtt);
	G4Material* layerTwoMaterial = G4Material::GetMaterial(fPhotonShieldLayerTwoMaterial);
	fPhotonShieldLayerTwoLog = new G4LogicalVolume(SecondLayer, layerTwoMaterial, "PhotonShieldLayerTwoLog",0,0,0);
	fPhotonShieldLayerTwoLog->SetVisAttributes(TwoVisAtt);
	G4Material* layerThreeMaterial = G4Material::GetMaterial(fPhotonShieldLayerThreeMaterial);
	fPhotonShieldLayerThreeLog = new G4LogicalVolume(ThirdLayer, layerThreeMaterial, "PhotonShieldLayerThreeLog",0,0,0);
	fPhotonShieldLayerThreeLog->SetVisAttributes(ThreeVisAtt);
} // end BuildPhotonShield()

void ApparatusSpiceTargetChamber::BuildShieldCovering()
{
	// ** Visualization
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0.0,0.2,0.8));
	VisAtt->SetVisibility(true);

	// ** Build Photon Shield Solid
	//G4double InnerRadius = fPhotonShieldInnerRadius; 
	G4double FrontOuterRadius = fPhotonShieldFrontRadius;
	G4double BackOuterRadius = fPhotonShieldBackRadius;
	G4double HalfLength = fPhotonShieldLength/2.;  
	G4Cons* PhotonShield = new G4Cons("PhotonShield", (BackOuterRadius + 0.1*CLHEP::mm), (BackOuterRadius + 0.2*CLHEP::mm),(FrontOuterRadius + 0.1*CLHEP::mm),
			(FrontOuterRadius + 0.2*CLHEP::mm), HalfLength, 0, 2*CLHEP::pi);//made smaller

	/*// ------ 																							
	// Photon Shield Clamps
	// -- Dimensions  																							
	G4double DetectorEndRadius = fPsDetClampRadius;
	G4double DetectorEndHalf = fPsDetClampThickness/2.;
	G4double TargetEndRadius = fPsTargetClampRadius;
	G4double TargetEndHalf = fPsTargetClampThickness/2.;				
	// -- Shapes																			
	G4Tubs* TargetClamp = new G4Tubs("TargetClamp", 0, TargetEndRadius, TargetEndHalf, 0, 360*CLHEP::deg);
	G4Tubs* DetClamp = new G4Tubs("DetClamp", 0, DetectorEndRadius, DetectorEndHalf, 0, 360*CLHEP::deg);
	// -- Union
	G4ThreeVector move(0, 0, HalfLength + TargetEndHalf - 0.01*CLHEP::mm);
	G4UnionSolid* PhotonShield_2 = new G4UnionSolid("PhotonShield_2", PhotonShield, TargetClamp, 0, move);

	G4UnionSolid* PhotonShield_3 = new G4UnionSolid("PhotonShield_3", PhotonShield_2, DetClamp, 0, move);	
	move.setZ(-HalfLength);
	// ------*/


	/*	// ** Build Shield Covering
		G4double theta = atan(2 * HalfLength / (BackOuterRadius - FrontOuterRadius));
		G4double alpha = 90*CLHEP::deg - theta;
		G4double CoverFrontRadius = FrontOuterRadius - (fShieldCoatingThickness / tan(theta)) +
		(fShieldCoatingThickness * cos(alpha));
		G4double CoverBackRadius = BackOuterRadius + (fShieldCoatingThickness / tan(theta)) +
		(fShieldCoatingThickness * cos(alpha));

		G4Cons* ShieldCoverPre = new G4Cons("ShieldCoverPre", 0, CoverBackRadius, 0, CoverFrontRadius,	
		HalfLength + fShieldCoatingThickness, 0, 2*CLHEP::pi);

	// ** Deduct Shield from Cover
	G4SubtractionSolid* ShieldCover = new G4SubtractionSolid("ShieldCover", ShieldCoverPre, PhotonShield_3);//LOOP*/


	// ------
	// Magnet Clamps
	// -- Dimensions
	G4double ClampHalfThickness = fPlateOneThickness/2. + fPsClampThicknessBuffer;
	G4Box* MagnetClamp = new G4Box("MagnetClamp", fPsClampLength/2., ClampHalfThickness, fPhotonShieldLength);

	// -- Rotation & Translation
	G4RotationMatrix* cutRot = new G4RotationMatrix; 
	G4double transRadial = fPlateOneEdgeX + fPsClampLength/2. - fPsClampPsbuffer;
	G4double transX, transY;
	G4double transZ = 0;
	G4ThreeVector trans;
	// -- Union
	cutRot->rotateZ(-(45.+0*90.)*CLHEP::deg);
	transX = transRadial*cos((0*90.+45.)*CLHEP::deg);
	transY = transRadial*sin((0*90.+45.)*CLHEP::deg);
	trans.set(transX,transY,transZ);
	G4SubtractionSolid* ShieldCover_2 = new G4SubtractionSolid("PhotonShield", PhotonShield,MagnetClamp,cutRot,trans);      
	cutRot->rotateZ(-90.*CLHEP::deg);
	transX = transRadial*cos((1*90.+45.)*CLHEP::deg);
	transY = transRadial*sin((1*90.+45.)*CLHEP::deg);
	trans.set(transX,transY,transZ);
	G4SubtractionSolid* ShieldCover_3 = new G4SubtractionSolid("PhotonShield",ShieldCover_2,MagnetClamp,cutRot,trans);      
	cutRot->rotateZ(-90.*CLHEP::deg);
	transX = transRadial*cos((2*90.+45.)*CLHEP::deg);
	transY = transRadial*sin((2*90.+45.)*CLHEP::deg);
	trans.set(transX,transY,transZ);
	G4SubtractionSolid* ShieldCover_4 = new G4SubtractionSolid("PhotonShield",ShieldCover_3,MagnetClamp,cutRot,trans);      
	cutRot->rotateZ(-90.*CLHEP::deg);
	transX = transRadial*cos((3*90.+45.)*CLHEP::deg);
	transY = transRadial*sin((3*90.+45.)*CLHEP::deg);
	trans.set(transX,transY,transZ);
	G4SubtractionSolid* ShieldCover_5 = new G4SubtractionSolid("PhotonShield",ShieldCover_4,MagnetClamp,cutRot,trans); 
	// -------

	G4Material* shieldCoverMaterial = G4Material::GetMaterial(fShieldCoverMaterial);
	fShieldCoverLog = new G4LogicalVolume(ShieldCover_5, shieldCoverMaterial, "ShieldCoverLog", 0, 0, 0);
	fShieldCoverLog->SetVisAttributes(VisAtt);
}

void ApparatusSpiceTargetChamber::BuildPhotonShieldClamps() {

	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(AL_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	G4double DetectorEndRadius = fPsDetClampRadius;
	G4double DetectorEndHalf = fPsDetClampThickness/2.;
	G4double TargetEndRadius = fPsTargetClampRadius;
	G4double TargetEndHalf = fPsTargetClampThickness/2.;
	G4double InnerRadius = fPhotonShieldInnerRadius;
	// Extensions
	G4double ExtHalfLength = fPsTargetClampRadius + fPsTargetClampExtension;
	G4double ExtHalfWidth = fPlateOneThickness/2. + fPsClampThicknessBuffer;
	// Bolts
	G4double TargetBoltRadius = fPsTargetBoltRadius;
	G4double TargetBoltOffset = fPsTargetBoltPlaneOffset;
	G4double DetectorBoltRadius = fPsDetBoltRadius;
	G4double DetectorBoltOffset = fPsDetBoltPlaneOffset;


	// ** Shapes
	G4Tubs* TargetEndClampPre = new G4Tubs("TargetEndClampPre", InnerRadius, TargetEndRadius, TargetEndHalf, 0, 360*CLHEP::deg);
	G4Box* TargetClampExt = new G4Box("TargetClampExt", ExtHalfLength, ExtHalfWidth, TargetEndHalf);
	G4Box* TargetClampExt2 = new G4Box("TargetClampExt2", ExtHalfWidth, ExtHalfLength, TargetEndHalf);
	G4UnionSolid* TargetExtensions = new G4UnionSolid("TargetExtensions", TargetClampExt, TargetClampExt2);
	G4Tubs* ExtensionCut = new G4Tubs("ExtensionCut", 0, InnerRadius + 2*CLHEP::mm, 2*TargetEndHalf, 0, 360*CLHEP::deg);
	G4SubtractionSolid* TargetExtHole = new G4SubtractionSolid("TargetExtHole", TargetExtensions, ExtensionCut);  
	G4UnionSolid* TargetClampPreBolts = new G4UnionSolid("TargetClampPreBolts", TargetEndClampPre, TargetExtHole);  

	G4Tubs* TargetBolt = new G4Tubs("TargetBolt", 0, TargetBoltRadius, 2.*TargetEndHalf, 0, 360*CLHEP::deg);
	G4ThreeVector MoveBolt(TargetBoltOffset, 0, 0);
	G4SubtractionSolid* TargetEndClamp0 = new G4SubtractionSolid("TargetEndClamp0", TargetClampPreBolts, TargetBolt, 0, MoveBolt);
	MoveBolt.setX(-TargetBoltOffset);
	G4SubtractionSolid* TargetEndClamp1 = new G4SubtractionSolid("TargetEndClamp1", TargetEndClamp0, TargetBolt, 0, MoveBolt);
	MoveBolt.setX(0);
	MoveBolt.setY(TargetBoltOffset);
	G4SubtractionSolid* TargetEndClamp2 = new G4SubtractionSolid("TargetEndClamp2", TargetEndClamp1, TargetBolt, 0, MoveBolt);
	MoveBolt.setY(-TargetBoltOffset);
	G4SubtractionSolid* TargetEndClamp3 = new G4SubtractionSolid("TargetEndClamp3", TargetEndClamp2, TargetBolt, 0, MoveBolt);

	G4Tubs* DetEndClampPre = new G4Tubs("DetEndClampPre", InnerRadius, DetectorEndRadius, DetectorEndHalf, 0, 360*CLHEP::deg);
	G4Tubs* DetectorBolt = new G4Tubs("DetectorBolt", 0, DetectorBoltRadius, 2.* DetectorEndHalf, 0, 360*CLHEP::deg);
	MoveBolt.setY(DetectorBoltOffset);
	G4SubtractionSolid* DetEndClamp0 = new G4SubtractionSolid("DetEndClamp0", DetEndClampPre, DetectorBolt, 0, MoveBolt);
	MoveBolt.setY(-DetectorBoltOffset);
	G4SubtractionSolid* DetEndClamp1 = new G4SubtractionSolid("DetEndClamp1", DetEndClamp0, DetectorBolt, 0, MoveBolt);
	MoveBolt.setY(0);
	MoveBolt.setX(DetectorBoltOffset);
	G4SubtractionSolid* DetEndClamp2 = new G4SubtractionSolid("DetEndClamp2", DetEndClamp1, DetectorBolt, 0, MoveBolt);
	MoveBolt.setX(-DetectorBoltOffset);
	G4SubtractionSolid* DetEndClamp3 = new G4SubtractionSolid("DetEndClamp3", DetEndClamp2, DetectorBolt, 0, MoveBolt);


	// ** Logical
	G4Material* psClampMaterial = G4Material::GetMaterial(fPsClampMaterial);
	fPsTargetClampLog = new G4LogicalVolume(TargetEndClamp3, psClampMaterial, "PsTargetClampLog", 0, 0, 0);
	fPsTargetClampLog->SetVisAttributes(VisAtt);
	fPsDetectorClampLog = new G4LogicalVolume(DetEndClamp3, psClampMaterial, "PsDetectorClampLog", 0, 0, 0);
	fPsDetectorClampLog->SetVisAttributes(VisAtt);
} // end::BuildPhotonShieldClamps

void ApparatusSpiceTargetChamber::BuildPhotonShieldClampBolts() {
	// ** Visualisation
	G4VisAttributes* SmallVisAtt = new G4VisAttributes(G4Colour(CU_COL));
	SmallVisAtt->SetVisibility(true);
	G4VisAttributes* LargeVisAtt = new G4VisAttributes(G4Colour(SN_COL));
	LargeVisAtt->SetVisibility(true);

	// ** Dimensions
	G4double TargetBoltHalfLength = fPsTargetBoltLength/2.;
	G4double TargetBoltRadius = fPsTargetBoltRadius;
	G4double TargetBoltHeadHalf = fPsTargetBoltHeadLength/2.;
	G4double TargetBoltHeadRad = fPsTargetBoltHeadRadius;
	G4double DetectorBoltHalfLength = fPsDetBoltLength/2. + //PipeZLength/2. + PsDetClampThickness/2.;
	fPsDetClampThickness/2.;
	G4double DetectorBoltRadius = fPsDetBoltRadius;
	G4double DetectorHeadHalf = fPsDetBoltHeadLength/2.;
	G4double DetectorHeadRadius = fPsDetBoltHeadRadius;

	// ** Shapes
	G4Tubs* TargetBoltShaft = new G4Tubs("TargetBoltShaft", 0, TargetBoltRadius, TargetBoltHalfLength, 0, 360*CLHEP::deg);
	G4Tubs* TargetBoltHead = new G4Tubs("TargetBoltHead", 0, TargetBoltHeadRad, TargetBoltHeadHalf, 0, 360*CLHEP::deg);
	G4ThreeVector trans(0, 0, TargetBoltHalfLength + TargetBoltHeadHalf);
	G4UnionSolid* TargetBolt = new G4UnionSolid("TargetBolt", TargetBoltShaft, TargetBoltHead, 0, trans);
	G4Tubs* DetectorBoltShaft = new G4Tubs("DetectorBoltShaft", 0, DetectorBoltRadius, DetectorBoltHalfLength, 0, 360*CLHEP::deg);
	G4Tubs* DetectorBoltHead = new G4Tubs("DetectorBoltHead", 0, DetectorHeadRadius, DetectorHeadHalf, 0, 360*CLHEP::deg);
	trans.setZ(-DetectorBoltHalfLength - DetectorHeadHalf);
	G4UnionSolid* DetectorBolt = new G4UnionSolid("DetectorBolt", DetectorBoltShaft, DetectorBoltHead, 0, trans);

	// ** Logical
	G4Material* smallBoltMaterial = G4Material::GetMaterial(fSmallBoltMaterial);
	fTargetBoltLog = new G4LogicalVolume(TargetBolt, smallBoltMaterial, "TargetBoltLog", 0, 0, 0);
	fTargetBoltLog->SetVisAttributes(SmallVisAtt);
	G4Material* largeBoltMaterial = G4Material::GetMaterial(fLargeBoltMaterial);
	fDetectorBoltLog = new G4LogicalVolume(DetectorBolt, largeBoltMaterial, "DetectorBoltLog", 0, 0, 0);
	fDetectorBoltLog->SetVisAttributes(LargeVisAtt);
} // end::BuildPhotonShieldClampBolts()

void ApparatusSpiceTargetChamber::BuildCollectorMagnet()
{

	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(NDFEB_COL));
	VisAtt->SetVisibility(true);

	// ** Build Plate One
	G4Box* PlateOnePreCut = new G4Box("PlateOnePreCut", fPlateOneLength/2., 
			fPlateOneThickness/2., fPlateOneHeight/2.);
	G4Box* PlateOneCutBox = new G4Box("PlateOneCutBox", fPlateOneLength/2., 
			fPlateOneThickness/2. + 1., fPlateOneHeight/2.);

	// ( Cut box is rotated and moved, exact movement calculated below )
	G4double AngleToCorner = atan(fPlateOneLength/fPlateOneHeight);
	G4double HypotenuseToCorner = sqrt(pow(fPlateOneHeight/2., 2) + pow(fPlateOneLength/2., 2));
	G4double AngleDifference = 0.;
	G4double TranslateCutX = 0.;
	G4double TranslateCutZ = 0.;
	if(AngleToCorner < fCuttingBoxAngle)	{
		AngleDifference = fCuttingBoxAngle - AngleToCorner;
		TranslateCutX = (HypotenuseToCorner * sin(AngleDifference)) + fPlateOneLength/2.;
		TranslateCutZ = (HypotenuseToCorner * cos(AngleDifference)) - fPlateOneHeight/2.
			+ fPlateOneLowerHeight;
	}
	else if(AngleToCorner > fCuttingBoxAngle)
	{
		AngleDifference = AngleToCorner - fCuttingBoxAngle;
		TranslateCutX = fPlateOneLength/2. - (HypotenuseToCorner * sin(AngleDifference));
		TranslateCutZ = (HypotenuseToCorner * cos(AngleDifference)) - fPlateOneHeight/2.
			+ fPlateOneLowerHeight;
	}
	G4ThreeVector TranslateCutBox(-TranslateCutX, 0, TranslateCutZ);
	G4RotationMatrix* RotateCutBox = new G4RotationMatrix;
	RotateCutBox->rotateY(fCuttingBoxAngle);
	G4SubtractionSolid* PlateOne = new G4SubtractionSolid("PlateOne", PlateOnePreCut, PlateOneCutBox,
			RotateCutBox, TranslateCutBox);

	// ** Build Plate Two
	G4Box* PlateTwo = new G4Box("PlateTwo", fPlateTwoLength/2., 
			(fPlateTwoThickness + fPlateOneThickness/2.),
			fPlateTwoHeight/2.);

	// ** Combine Plates
	G4ThreeVector TranslatePlateTwo(fPlateOneLength/2. - fPlateTwoLength/2., 0, 0);
	G4UnionSolid* PlateCombo = new G4UnionSolid("PlateCombo", PlateOne, PlateTwo, 0, TranslatePlateTwo);

	// ** Logical Volume
	G4Material* magnetMaterial = G4Material::GetMaterial(fMagnetMaterial);
	fMagnetLog = new G4LogicalVolume(PlateCombo, magnetMaterial, "MagnetLog", 0, 0, 0);
	fMagnetLog->SetVisAttributes(VisAtt);

} // end:BuildCollectorMagnet()

void ApparatusSpiceTargetChamber::BuildMagnetCovering()
{
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(PEEK_COL));
	VisAtt->SetVisibility(true);

	// ( Cut box is rotated and moved, exact movement calculated below )
	G4double AngleToCorner = atan(fPlateOneLength/fPlateOneHeight);
	G4double HypotenuseToCorner = sqrt(pow(fPlateOneHeight/2., 2) + pow(fPlateOneLength/2., 2));
	G4double AngleDifference = 0.;
	G4double TranslateCutX = 0.;
	G4double TranslateCutZ = 0.;
	if(AngleToCorner < fCuttingBoxAngle)
	{
		AngleDifference = fCuttingBoxAngle - AngleToCorner;
		TranslateCutX = (HypotenuseToCorner * sin(AngleDifference)) + fPlateOneLength/2.;
		TranslateCutZ = (HypotenuseToCorner * cos(AngleDifference)) - fPlateOneHeight/2.
			+ fPlateOneLowerHeight;
	}
	else if(AngleToCorner > fCuttingBoxAngle)
	{
		AngleDifference = AngleToCorner - fCuttingBoxAngle;
		TranslateCutX = fPlateOneLength/2. - (HypotenuseToCorner * sin(AngleDifference));
		TranslateCutZ = (HypotenuseToCorner * cos(AngleDifference)) - fPlateOneHeight/2.
			+ fPlateOneLowerHeight;
	}
	G4ThreeVector TranslateCutBox(-TranslateCutX, 0, TranslateCutZ);
	G4RotationMatrix* RotateCutBox = new G4RotationMatrix;
	RotateCutBox->rotateY(fCuttingBoxAngle);

	// ** Combine Plates
	G4ThreeVector TranslatePlateTwo(fPlateOneLength/2. - fPlateTwoLength/2., 0, 0);

	// ------------
	// Photon Shield Clamp
	// -- Dimensions
	G4double ClampHalfLength = (fPlateOneLength - fPlateTwoLength + fPsClampPsbuffer + fPsClampChamberBuffer)/2.;
	G4double ClampHalfHeight = fPhotonShieldLength/2.;
	// -- Union
	G4double xtrans = (fPlateOneLength + fPsClampLength)/2. - fPsClampLength + fPsClampPsbuffer;
	G4double ztrans = -fPhotonShieldBackFacePos - fPlateOneHeight/2. - fDistanceFromTarget - fPhotonShieldLength/2.;
	G4ThreeVector MoveClamp(-xtrans, 0, -ztrans);
	// --------------


	// --------------
	// Chamber Clamps
	// -- Dimensions
	ClampHalfLength = fMclampChamberLength/2.;
	ClampHalfHeight = fMclampChamberHeight/2.;
	// -- Shape
	xtrans = -(ClampHalfLength + fPlateOneLength/2.) + 
		(fMclampChamberLength - (fTargetChamberCylinderInnerRadius - (fPlateOneLength + fPlateOneEdgeX)));
	G4ThreeVector MoveAgain(-xtrans, 0, -ClampHalfHeight + fPlateOneHeight/2.);
	// ---------------


	// ------------------
	// Build Magnet Cover

	// ** Cover of Plate One
	G4Box* PlateOneCoverPreCut = new G4Box("PlateOneCoverPreCut", fPlateOneLength/2. + fMagnetCoatingThickness, 
			fPlateOneThickness/2. + fMagnetCoatingThickness, 
			fPlateOneHeight/2. + fMagnetCoatingThickness);
	G4Box* PlateOneCoverCutBox = new G4Box("PlateOne_CoverCutBox", 
			fPlateOneLength/2., 
			fPlateOneThickness + fMagnetCoatingThickness, 
			fPlateOneHeight/2. - fMagnetCoatingThickness);
	G4SubtractionSolid* PlateOneCover = new G4SubtractionSolid("PlateOneCover", PlateOneCoverPreCut, 
			PlateOneCoverCutBox, RotateCutBox, TranslateCutBox); 

	// ** Build Plate Two Cover
	G4Box* PlateTwoCover = new G4Box("PlateTwoCover", fPlateTwoLength/2. + fMagnetCoatingThickness, 
			fPlateTwoThickness + fPlateOneThickness/2. + fMagnetCoatingThickness,
			fPlateTwoHeight/2. + fMagnetCoatingThickness);

	// ** Combine Plate Covers
	G4UnionSolid* CoverCombo = new G4UnionSolid("CoverCombo", PlateOneCover, PlateTwoCover, 0, TranslatePlateTwo);

	// ---------------------
	// ** Logical Volume
	G4Material* magnetCoverMaterial = G4Material::GetMaterial(fMagnetCoverMaterial);
	fMagnetCoverLog = new G4LogicalVolume(CoverCombo, magnetCoverMaterial, "MagnetCoverLog", 0, 0, 0);
	fMagnetCoverLog->SetVisAttributes(VisAtt);
}

void ApparatusSpiceTargetChamber::BuildMagnetClampChamber(){

	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(AL_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	G4double ClampHalfThickness = fMclampChamberThickness/2.;
	G4double ClampHalfLength = fMclampChamberLength/2.;
	G4double ClampHalfHeight = fMclampChamberHeight/2.;

	G4double MagnetHalfThickness = (fPlateOneThickness + 2.*fPlateTwoThickness)/2.;
	G4double MagnetHalfLength = fPlateTwoLength/2.;
	G4double MagnetHalfHeight = fPlateTwoHeight/2.;

	// ** Shapes
	G4Box* ClampMain = new G4Box("ClampMain", ClampHalfLength, ClampHalfThickness, ClampHalfHeight);
	G4Box* MagnetCut = new G4Box("MagnetCut", MagnetHalfLength, MagnetHalfThickness, MagnetHalfHeight);
	G4double xtrans = -(ClampHalfLength + MagnetHalfLength) + 
		(fMclampChamberLength - (fTargetChamberCylinderInnerRadius - (fPlateOneLength + fPlateOneEdgeX)));
	G4ThreeVector trans(xtrans, 0, ClampHalfHeight - MagnetHalfHeight);
	G4SubtractionSolid* MagnetClamp = new G4SubtractionSolid("MagnetClamp", ClampMain, MagnetCut, 0, trans);

	// ** Logical
	G4Material* clampMaterial = G4Material::GetMaterial(fMagnetClampMaterial);
	fMclampChamberLog = new G4LogicalVolume(MagnetClamp, clampMaterial, "MagnetClampChamberLog", 0, 0, 0);
	fMclampChamberLog->SetVisAttributes(VisAtt);
} // end::BuildMagnetClampChamber()


void ApparatusSpiceTargetChamber::BuildElectroBox()
{
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(ELECTROBOX_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	G4double OuterBoxHalfWidth = fElectroboxOuterBoxLength/2.;
	G4double OuterBoxHalfLength = fElectroboxMidLength/2.;
	G4double InnerBoxHalfWidth = fElectroboxInnerBoxLength/2.;
	G4double InnerBoxHalfLength = fElectroboxMidLength/2.;
	G4double InnerCylRadius = fElectroboxLipInnerRadius;


	// ** Shapes
	G4Box* OuterBox = new G4Box("OuterBox", OuterBoxHalfWidth, OuterBoxHalfWidth, OuterBoxHalfLength);
	G4Box* InnerBox = new G4Box("InnerBox", InnerBoxHalfWidth, InnerBoxHalfWidth, InnerBoxHalfLength);
	G4Tubs* InnerCyl = new G4Tubs("InnerCyl", 0, InnerCylRadius, 2*InnerBoxHalfLength, 0, 360*CLHEP::deg);
	G4UnionSolid* TotalInner = new G4UnionSolid("TotalInner", InnerBox, InnerCyl);
	G4ThreeVector trans(0, 0, -fElectroboxLipLength);
	G4SubtractionSolid* BoxMain = new G4SubtractionSolid("BoxMain", OuterBox, TotalInner, 0, trans);

	// ** Logical
	G4Material* electroBoxMaterial = G4Material::GetMaterial(fElectroBoxMaterial);
	fElectroBoxLog = new G4LogicalVolume(BoxMain, electroBoxMaterial, "ElectroBoxLog", 0, 0, 0);
	fElectroBoxLog->SetVisAttributes(VisAtt);
} // end:BuildElectroBox()


void ApparatusSpiceTargetChamber::BuildColdFinger()
{
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(CU_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	G4double ColdfingerHalfThickness = fColdfingerThickness/2.;
	G4double ColdfingerHalfLength = fColdfingerLength/2.;
	G4double ColdfingerHalfWidth = fColdfingerWidth/2.;
	G4double InnerHoleRadius = fColdfingerHoleRadius;

	// ** Shapes
	G4Box* PlateBox = new G4Box("PlateBox", ColdfingerHalfLength, ColdfingerHalfWidth, ColdfingerHalfThickness);
	G4Tubs* InnerHole = new G4Tubs("InnerHole", 0, InnerHoleRadius, 2*ColdfingerHalfThickness, 0, 360*CLHEP::deg);
	G4SubtractionSolid* PlateSubHole = new G4SubtractionSolid("plate-hole", PlateBox , InnerHole, 0, G4ThreeVector(0,0,0));

	// ** Logical
	G4Material* coldFingerMaterial = G4Material::GetMaterial(fColdFingerMaterial);
	fColdFingerLog = new G4LogicalVolume(PlateSubHole, coldFingerMaterial, " ColdFingerLog", 0, 0, 0);
	fColdFingerLog->SetVisAttributes(VisAtt);
} // end:BuildColdFinger()

void ApparatusSpiceTargetChamber::BuildS3CableHolder()
{
	// Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(CU_COL));
	sVisAtt->SetVisibility(true);

	// Dimensions
	// Shapes
	G4Box* sBox = new G4Box("sBox", 8.1*CLHEP::mm, 12.*CLHEP::mm, 40.5*CLHEP::mm);

	// Logical
	G4Material* sMaterial = G4Material::GetMaterial(fS3CableCaseMaterial);
	fS3CaseLogical = new G4LogicalVolume(sBox, sMaterial, "s3CaseLog", 0, 0, 0);
	fS3CaseLogical->SetVisAttributes(sVisAtt);
}
void ApparatusSpiceTargetChamber::BuildConicalCollimator() {

	// ** Dimensions


	// ** shapes
	G4Tubs* solid_mid_cyl = new G4Tubs("mid_cyl", 3.15*mm, 4.5*mm, 20.25*mm, 0, 360*CLHEP::deg);
	G4Tubs* solid_inner_cyl = new G4Tubs("inner_cyl", 2.9*mm, 3.15*mm, 7.75*mm, 0, 360*CLHEP::deg);
	G4Tubs* solid_outer_cyl = new G4Tubs("outer_cyl", 4.5*mm, 4.9*mm, 14.25*mm, 0, 360*CLHEP::deg);
	G4Tubs* edge1 = new G4Tubs("edges", 4.9*mm, 6.9*mm, 1.0*mm, -20.*CLHEP::deg, 40.*CLHEP::deg);
	G4Tubs* edge2 = new G4Tubs("edges", 4.9*mm, 6.9*mm, 1.0*mm, 70.*CLHEP::deg, 40.*CLHEP::deg);
	G4Tubs* edge3 = new G4Tubs("edges", 4.9*mm, 6.9*mm, 1.0*mm, 160.*CLHEP::deg, 40.*CLHEP::deg);
	G4Tubs* edge4 = new G4Tubs("edges", 4.9*mm, 6.9*mm, 1.0*mm, 250.*CLHEP::deg, 40.*CLHEP::deg);

	G4ThreeVector sTrans1(0., 0., -12.5*mm);
	G4ThreeVector sTrans2(0., 0., -6.*mm);
	G4ThreeVector sTrans3(0, 0, -19.25*mm);

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	G4UnionSolid* midandinner = new G4UnionSolid("Mid+Outer", solid_mid_cyl, solid_inner_cyl, sRotate, sTrans1);
	G4UnionSolid* threecyl = new G4UnionSolid("ThreeCyl", midandinner, solid_outer_cyl, sRotate, sTrans2);
	G4UnionSolid* threecyl1edge = new G4UnionSolid("threecyl1edge", threecyl, edge1, sRotate, sTrans3);
	G4UnionSolid* threecyl2 = new G4UnionSolid("Threecyl2", threecyl1edge, edge2, sRotate, sTrans3);
	G4UnionSolid* threecyl3 = new G4UnionSolid("ThreeCyl3", threecyl2, edge3, sRotate, sTrans3);
	G4UnionSolid* threecyl4edge = new G4UnionSolid("ThreeCyl4edge", threecyl3, edge4, sRotate, sTrans3);

	// ** Logical
	G4Material* sMaterial = G4Material::GetMaterial(this->fConicalCollimatorMaterial);
	fConicalCollimatorLog = new G4LogicalVolume(threecyl4edge,         //its solid
			sMaterial,          //its material
			"threecyl4edge");           //its name

}

void ApparatusSpiceTargetChamber::BuildXrayInsert(){
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(CU_COL));
	sVisAtt->SetVisibility(true);

	// ** shapes
	G4Tubs* solid_insert_cyl = new G4Tubs("insert_cyl", 3.68*mm, 8.5*mm, 1.5*mm, 0, 360*CLHEP::deg);
	G4Tubs* solid_insert_edges = new G4Tubs("insert_edge", 4.55*mm, 8.5*mm, 1.8*mm, 25.*CLHEP::deg, 40.*CLHEP::deg);
	G4Tubs* solid_insert_holes = new G4Tubs("insert_hole", 0.*mm, 2.*mm, 5.*mm, 0, 360*CLHEP::deg);

	G4RotationMatrix* sRotate1 = new G4RotationMatrix;
	G4ThreeVector sTrans4(0., 0., -3.3*mm);
	G4UnionSolid* cyland1 = new G4UnionSolid("cyl+edge1", solid_insert_cyl, solid_insert_edges, sRotate1, sTrans4);
	sRotate1->rotateZ(90.*CLHEP::deg);
	G4UnionSolid* cyland2 = new G4UnionSolid("cyl+edge2", cyland1, solid_insert_edges, sRotate1, sTrans4);
	sRotate1->rotateZ(90.*CLHEP::deg);
	G4UnionSolid* cyland3 = new G4UnionSolid("cyl+edge3", cyland2, solid_insert_edges, sRotate1, sTrans4);
	sRotate1->rotateZ(90.*CLHEP::deg);
	G4UnionSolid* cyland4 = new G4UnionSolid("cyl+edge4", cyland3, solid_insert_edges, sRotate1, sTrans4);

	sRotate1->rotateZ(45.*CLHEP::deg);
	G4ThreeVector sTrans5(0., 8.5*mm, 3.3*mm);
	G4SubtractionSolid* c4hole4edge1 = new G4SubtractionSolid("c4hole4", cyland4, solid_insert_holes, sRotate1, sTrans5);
	sTrans5 = G4ThreeVector(0., -8.5*mm, 3.3*mm);
	G4SubtractionSolid* c4hole4edge2 = new G4SubtractionSolid("c4hole4", c4hole4edge1, solid_insert_holes, sRotate1, sTrans5);
	sTrans5 = G4ThreeVector(8.5*mm, 0., 3.3*mm);
	G4SubtractionSolid* c4hole4edge3 = new G4SubtractionSolid("c4hole4", c4hole4edge2, solid_insert_holes, sRotate1, sTrans5);
	sTrans5 = G4ThreeVector(-8.5*mm, 0., 3.3*mm);
	G4SubtractionSolid* c4hole4complete = new G4SubtractionSolid("c4hole4", c4hole4edge3, solid_insert_holes, sRotate1, sTrans5);

	// ** Logical  
	G4Material* sMaterial = G4Material::GetMaterial(this->fXRayInsertMaterial);
	fXRayInsertLog = new G4LogicalVolume(c4hole4complete,         //its solid
			sMaterial,          //its material
			"c4hole4");           //its name
	fXRayInsertLog->SetVisAttributes(sVisAtt);
}

// **************************************************************************************************
// **************************************************************************************************
// *************************************PLACEMENT****************************************************
// **************************************************************************************************
// **************************************************************************************************

void ApparatusSpiceTargetChamber::PlaceTargetFrame(G4double d, G4LogicalVolume* expHallLog) {
	G4double sOffset = fTargetWheelOffset;
	G4double sOffsetAdj = sOffset - 15.00;

	G4double sPosZ = -7.75*CLHEP::mm; 
	G4double sPosY = (sOffset-sOffsetAdj)*sin(180*CLHEP::deg + d);
	G4double sPosX = sOffset + (sOffset-sOffsetAdj)*cos(180*CLHEP::deg + d);

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(-d);   

	G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);
	fFramePhysical = new G4PVPlacement(sRotate, sTranslate, fFrameLogical, "TargetFrame", expHallLog, false, 0);
}

void ApparatusSpiceTargetChamber::PlaceTargetWheelRods(G4double a, G4double r, G4LogicalVolume* expHallLog) {
	G4double sOffset = fTargetWheelOffset; 
	G4double sOffsetAdj = sOffset - r;

	G4double sPosZ = -3.75*CLHEP::mm; 
	G4double sPosY = (sOffset-sOffsetAdj)*sin(180*CLHEP::deg + a);
	G4double sPosX = sOffset + (sOffset-sOffsetAdj)*cos(180*CLHEP::deg + a);

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(-a);   

	G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);
	fRodPhysical = new G4PVPlacement(sRotate, sTranslate, fRodLogical, "TargetRods", expHallLog, false, 0);
}

void ApparatusSpiceTargetChamber::PlaceCollimator(G4LogicalVolume* expHallLog) {
	G4double sOffset = fTargetWheelOffset; 
	G4double sOffsetAdj = sOffset - 15.2;

	// 	G4double sPosZ = -7.75*CLHEP::mm; 
	G4double sPosZ = -0.5*CLHEP::mm; 
	G4double sPosY = (sOffset-sOffsetAdj)*sin(-55*CLHEP::deg);
	G4double sPosX = sOffset + (sOffset-sOffsetAdj)*cos(-55*CLHEP::deg);

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(55*CLHEP::deg);   

	G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);
	fCollimPhysical = new G4PVPlacement(sRotate, sTranslate, fCollimLogical, "collimator", expHallLog, false, 0);
}

void ApparatusSpiceTargetChamber::PlaceTargetWheelExtension(G4LogicalVolume* expHallLog) {
	G4double sPosZ = 8.0*CLHEP::mm; 
	G4double sPosY = 0.*CLHEP::mm; 
	G4double sPosX = fTargetWheelOffset;

	G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);
	fExtPhysical = new G4PVPlacement(0, sTranslate, fExtLogical, "WheelExtension", expHallLog, false, 0);
}


void ApparatusSpiceTargetChamber::PlaceTargetChamberFrontRing(G4LogicalVolume* expHallLog)
{
	G4double ZPosition = fTargetChamberSphereCentre
		+ fTargetChamberSphereCutOffPoint + fFrontRingLength/2. + fFrontDomeOffset;

	G4ThreeVector move(0,0,ZPosition);
	fTargetChamberFrontRingPhys = new G4PVPlacement(0,move,fTargetChamberFrontRingLog,
			"TargetChamberFrontRing",expHallLog,
			false,0);

	move.setZ(ZPosition +  fFrontRingLength/2. + fFrontRingCylinderLength + fConeLength/2.);
	fTargetChamberFrontConePhys = new G4PVPlacement(0, move, fTargetChamberFrontConeLog,
			"TargetChamberFrontCone", expHallLog,
			false,0);
}// end:PlaceTargetChamberFront()

void ApparatusSpiceTargetChamber::PlaceTargetChamberSphere(G4LogicalVolume* expHallLog)
{
	G4double ZPosition = fTargetChamberSphereCentre + fFrontDomeOffset;
	G4ThreeVector move(0,0,ZPosition);
	fTargetChamberSpherePhys = new G4PVPlacement(0,move, fTargetChamberSphereLog,
			"TargetChamberSphere",expHallLog,
			false,0);
}// end:PlaceTargetChamberSphere()

void ApparatusSpiceTargetChamber::PlaceTargetChamberCylinderDownstream(G4LogicalVolume* expHallLog)
{
	G4double ZPosition = -(fTargetChamberDistanceFromTarget + fTargetChamberCylinderLength/2. + fMiddleRingOffset);
	G4ThreeVector move(0,0,ZPosition);
	fTargetChamberCylinderDownPhys = new G4PVPlacement(0,move,fTargetChamberCylinderDownLog,
			"TargetChamberCylinderDownstream",expHallLog,
			false,0);
}// end:PlaceTargetChamberCylinderDownstream()

void ApparatusSpiceTargetChamber::PlaceTargetWheel(G4LogicalVolume* expHallLog)
{
	G4RotationMatrix* rotate = new G4RotationMatrix;  
	rotate->rotateZ((0)*CLHEP::deg); 

	G4double offset = fTargetWheelOffset;
	G4double ZPosition =  fTargetWheelThickness/2.; 

	G4ThreeVector move(offset, 0, ZPosition);
	fTargetWheelPhys = new G4PVPlacement(rotate, move, fTargetWheelLog,
			"TargetWheel", expHallLog,
			false,0);
}// end:PlaceTargetWheel()

void ApparatusSpiceTargetChamber::PlaceTargetWheelGears(G4LogicalVolume* expHallLog)
{
	G4double ZOffset = fAllGearZOffset + fAllGearThickness/2. + fTargetWheelOffset;
	G4double FirstGearOffset = fFirstGearPlaneOffset;
	G4double SecondGearOffset = fSecondGearPlaneOffset;
	G4double ThirdGearOffset = fThirdGearPlaneOffset;

	G4ThreeVector move(-FirstGearOffset, 0, ZOffset);
	fFirstGearPhys = new G4PVPlacement(0, move, fFirstGearLog,
			"FirstGear", expHallLog,
			false,0);

	move.setX(-SecondGearOffset);
	fSecondGearPhys = new G4PVPlacement(0, move, fSecondGearLog,
			"SecondGear", expHallLog,
			false,0);

	move.setX(ThirdGearOffset);
	fThirdGearPhys = new G4PVPlacement(0, move, fThirdGearLog,
			"ThirdGear", expHallLog,
			false,0);
} // end PlaceTargetWheelGears()

void ApparatusSpiceTargetChamber::PlaceTargetWheelGearPlates(G4LogicalVolume* expHallLog)
{
	G4double ZOffset = fAllGearZOffset 
		+ fAllGearThickness 
		+ fGearPlateThickness/2.
		+ fTargetWheelOffset;

	G4double FirstGearOffset = fFirstGearPlaneOffset;
	G4double SecondGearOffset = fSecondGearPlaneOffset;

	G4ThreeVector move(-FirstGearOffset, 0, ZOffset);
	fGearPlateOnePhys = new G4PVPlacement(0, move, fGearPlateOneLog,
			"GearPlateOne", expHallLog,
			false,0);
	move.setX(-SecondGearOffset);
	fGearPlateOnePhys = new G4PVPlacement(0, move, fGearPlateTwoLog,
			"GearPlateTwo", expHallLog,
			false,0);	
}	// end::PlaceTargetWheelGearPlates()

void ApparatusSpiceTargetChamber::PlaceGearStick(G4LogicalVolume* expHallLog)
{
	G4double ZOffset = -fGearStickLength/2. + fGearStickZOffset + fTargetWheelOffset;
	G4double PlaneOffset = fFirstGearPlaneOffset;

	G4ThreeVector move(-PlaneOffset, 0 , ZOffset);
	fGearStickPhys = new G4PVPlacement(0, move, fGearStickLog,
			"GearStick", expHallLog,
			false,0);
} // end PlaceGearStick();

void ApparatusSpiceTargetChamber::PlaceTargetMountPlate(G4LogicalVolume* expHallLog) 
{
	G4RotationMatrix* rotate = new G4RotationMatrix;  
	rotate->rotateZ(45*CLHEP::deg); 

	G4double ZOffset = fTargetMountPlateZOffset + fTargetMountPlateThickness/2. + fTargetWheelOffset;
	G4ThreeVector move(0, 0, ZOffset);
	fTargetMountPlatePhys = new G4PVPlacement(rotate, move, fTargetMountPlateLog,
			"TargetMountPlate", expHallLog,
			false,0);
} // end PlaceTargetMountPlate()

void ApparatusSpiceTargetChamber::PlacePhotonShield(G4LogicalVolume* expHallLog)
{
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*CLHEP::deg);

	G4double ZPosition =  fPhotonShieldBackFacePos
		+ fPhotonShieldLayerThreeThickness/2.
		+ fPhotonShieldLayerTwoThickness/2. 
		+ fPhotonShieldLayerOneThickness/2.
		+ fMiddleRingOffset;

	G4ThreeVector move(0,0,ZPosition);
	fPhotonShieldLayerOnePhys = new G4PVPlacement(rotate,move, fPhotonShieldLayerOneLog, "PhotonShieldLayerOne",expHallLog,
			false,0);

	move.set(0,0, fPhotonShieldBackFacePos
			+ fPhotonShieldLayerThreeThickness/2.
			+ fPhotonShieldLayerTwoThickness/2.
			+ fPhotonShieldLayerOneThickness/2.);
	fPhotonShieldLayerTwoPhys = new G4PVPlacement(rotate,move, fPhotonShieldLayerTwoLog,
			"PhotonShieldLayerTwo", expHallLog,
			false,0);

	move.set(0,0, fPhotonShieldBackFacePos
			+ fPhotonShieldLayerThreeThickness/2. 
			+ fPhotonShieldLayerTwoThickness/2. 
			+ fPhotonShieldLayerOneThickness/2.);
	fPhotonShieldLayerThreePhys = new G4PVPlacement(rotate,move, fPhotonShieldLayerThreeLog,
			"PhotonShieldLayerThree", expHallLog,
			false,0);

}// end:PlacePhotonShield()

void ApparatusSpiceTargetChamber::PlaceShieldCovering(G4LogicalVolume* expHallLog)
{

	G4double ZPosition = fPhotonShieldBackFacePos + fPhotonShieldLength/2. + fMiddleRingOffset;

	G4ThreeVector move(0,0,ZPosition);

	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateX(0);
	rotate->rotateY(0);
	rotate->rotateZ(0*CLHEP::deg);

	fShieldCoverPhys = new G4PVPlacement(rotate,move,fShieldCoverLog,
			"ShieldCover",expHallLog,
			false, 0);
}// end:PlaceShieldCovering()

void ApparatusSpiceTargetChamber::PlacePhotonShieldClamps(G4LogicalVolume* expHallLog)
{
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(45*CLHEP::deg);

	G4double ZPosition = fPhotonShieldBackFacePos + fMiddleRingOffset;

	G4ThreeVector move(0, 0, ZPosition - fPsDetClampThickness/2.);
	fPsDetectorClampPhys = new G4PVPlacement(rotate, move, fPsDetectorClampLog,
			"PhotonShieldDetectorClamp",expHallLog,
			false,0);

	move.setZ(ZPosition + fPhotonShieldLength + fPsTargetClampThickness/2.);
	fPsTargetClampPhys = new G4PVPlacement(rotate, move, fPsTargetClampLog,
			"PhotonShieldTargetClamp",expHallLog,
			false,0);
} // end::PlacePhotonShieldClamps()

void ApparatusSpiceTargetChamber::PlacePhotonShieldClampBolts(G4int copyID, G4LogicalVolume* expHallLog)
{
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(45*CLHEP::deg);

	G4double TargetBoltZPosition = fPhotonShieldBackFacePos 
		+ fPhotonShieldLength + fPsTargetClampThickness 
		- fPsTargetBoltLength/2. + fMiddleRingOffset;
	G4double DistanceFromBeamLine = fPsTargetBoltPlaneOffset;

	G4ThreeVector move = TranslateMagnets(copyID, DistanceFromBeamLine, TargetBoltZPosition);
	fTargetBoltPhys = new G4PVPlacement(rotate, move, fTargetBoltLog,
			"TargetBolt", expHallLog, 
			false,0);

	G4double DetectorBoltZPosition = fPhotonShieldBackFacePos  + fPsDetBoltLength/2. + fPipeZLength/2. + //fMiddleRingOffset - PipeZLength - PsDetClampThickness/2.;
	fMiddleRingOffset - fPsDetClampThickness/2.;
	DistanceFromBeamLine = fPsDetBoltPlaneOffset;
	G4ThreeVector move2 = TranslateMagnets(copyID, DistanceFromBeamLine, DetectorBoltZPosition);

	fDetectorBoltPhys = new G4PVPlacement(rotate, move2, fDetectorBoltLog,
			"DetectorBolt", expHallLog,
			false,0);
} // end::PlacePhotonShieldClampBoots()

void ApparatusSpiceTargetChamber::PlaceCollectorMagnet(G4int copyID, G4LogicalVolume* expHallLog)
{
	// ** Position Co-ordinates
	G4double magnetPosX = fPlateOneEdgeX + fPlateOneLength/2.;
	G4double magnetPosZ = -(fDistanceFromTarget + fPlateOneHeight/2.) + fMiddleRingOffset; // - (4mm + 50/2mm)  + 0 

	G4RotationMatrix* rotation;
	rotation = RotateMagnets(copyID);
	G4double RadialPosition = magnetPosX; 

	G4ThreeVector move = TranslateMagnets(copyID, RadialPosition, magnetPosZ);

	// ** Physical Volume
	fMagnetPhys = new G4PVPlacement(rotation, move, fMagnetLog,
			"CollectorMagnet",expHallLog,
			false,0);

} // end:PlaceCollectorMagnet()

void ApparatusSpiceTargetChamber::PlaceMagnetCovering(G4int copyID, G4LogicalVolume* expHallLog)
{
	// ** Position Co-ordinates
	G4double magnetPosX = fPlateOneEdgeX + fPlateOneLength/2.;
	G4double magnetPosZ = -(fDistanceFromTarget + fPlateOneHeight/2. ) + fMiddleRingOffset; 

	G4RotationMatrix* rotation;
	rotation = RotateMagnets(copyID);
	G4double RadialPosition = magnetPosX; 

	G4ThreeVector move = TranslateMagnets(copyID, RadialPosition, magnetPosZ);

	// ** Physical Volume
	fMagnetCoverPhys = new G4PVPlacement(rotation, move, fMagnetCoverLog,
			"MagnetCovering", expHallLog,
			false,0);
}

void ApparatusSpiceTargetChamber::PlaceMagnetClampChamber(G4int copyID, G4LogicalVolume* expHallLog)
{
	G4double ZPosition = -(fDistanceFromTarget + fMclampChamberHeight/2.) + fMiddleRingOffset;
	G4double XPosition = fTargetChamberCylinderInnerRadius - fMclampChamberLength/2.;

	G4RotationMatrix* rotation;
	rotation = RotateMagnets(copyID);

	G4ThreeVector move = TranslateMagnets(copyID, XPosition, ZPosition);

	fMclampChamberPhys = new G4PVPlacement(rotation, move, fMclampChamberLog,
			"MagnetClampChamber", expHallLog,
			false,0);
} // end::PlaceMagnetClampChamber()


void ApparatusSpiceTargetChamber::PlaceElectroBox(G4LogicalVolume* expHallLog)
{
	G4double ZOffset = fElectroboxZOffset - fElectroboxMidLength/2. + fBackDetAndAlBoxOffset;

	G4ThreeVector move(0, 0, ZOffset);

	fElectroBoxPhys = new G4PVPlacement(0, move, fElectroBoxLog,
			"ElectroBox", expHallLog, 
			false, 0);
} // end:PlaceElectroBox()

void ApparatusSpiceTargetChamber::PlaceColdFinger(G4LogicalVolume* expHallLog)
{
	G4double ZOffset = fColdfingerZOffset  - fColdfingerThickness/2.  ;

	G4ThreeVector move(0, 0, ZOffset);

	fColdFingerPhys = new G4PVPlacement(0, move, fColdFingerLog,
			"ColdFinger", expHallLog, 
			false, 0);
} // end:PlaceColdFinger()

void ApparatusSpiceTargetChamber::PlaceS3CableHolder(G4LogicalVolume* expHallLog)
{
	G4double sPosZ = -(40.5 + 6.)*CLHEP::mm;   
	G4double sPosY = 87.1*sin(-22.5*CLHEP::deg);
	G4double sPosX = -87.1*cos(-22.5*CLHEP::deg);

	G4ThreeVector sTranslate(sPosX, sPosY, sPosZ);

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(-22.5*CLHEP::deg);

	fS3CasePhysical = new G4PVPlacement(sRotate, sTranslate, fS3CaseLogical, "s3CableCase", expHallLog, false, 0);
}
void ApparatusSpiceTargetChamber::PlaceConicalCollimator(G4LogicalVolume* expHallLog) {
	G4ThreeVector pos1 = G4ThreeVector(0, 0, -44.*mm);//addition to the photon shield
	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateY(180*deg);
	fConicalCollimatorPhys = new G4PVPlacement(sRotate,                       //no rotation
			pos1,                    //at position
			fConicalCollimatorLog,             //its logical volume
			"ConicalCollimator",                //its name
			expHallLog,                //its mother  volume
			false,                   //no boolean operation
			false,                       //copy number
			0);          //overlaps checking

}

void ApparatusSpiceTargetChamber::PlaceXrayInsert(G4LogicalVolume* expHallLog){

	G4double ZPosition = fPhotonShieldBackFacePos + fMiddleRingOffset -fPsDetClampThickness -5.1*mm;
	G4ThreeVector move(0, 0, ZPosition);

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateY(180*deg);
	sRotate->rotateZ(45*deg);
	fXRayInsertPhys = new G4PVPlacement(sRotate,                       //rotation
			move,                    //at position
			fXRayInsertLog,             //its logical volume
			"XrayClamp",                //its name
			expHallLog,                //its mother  volume
			false,                   //no boolean operation
			false,                       //copy number
			0);          //overlaps checking
}
////////////////////////////////////////////////////////////////////
//Translations and Rotations of Magnets (same for all parts)      //
////////////////////////////////////////////////////////////////////

G4RotationMatrix* ApparatusSpiceTargetChamber::RotateMagnets(G4int copyID)
{
	G4RotationMatrix* rotate = new G4RotationMatrix;
	if(fNumberOfMagnets == 8) rotate->rotateZ(-(copyID+0.5)*45.*CLHEP::deg);
	else                      rotate->rotateZ(-(copyID+0.5)*90.*CLHEP::deg);

	return rotate;
}

G4ThreeVector ApparatusSpiceTargetChamber::TranslateMagnets(G4int copyID, G4double RadialPosition, G4double ZPosition)
{
	G4double XPosition = 0.;
	G4double YPosition = 0.;
	if(fNumberOfMagnets == 8){
		XPosition = RadialPosition*cos((copyID+0.5)*45.*CLHEP::deg);
		YPosition = RadialPosition*sin((copyID+0.5)*45.*CLHEP::deg);
	} else {
		XPosition = RadialPosition*cos((copyID+0.5)*90.*CLHEP::deg);
		YPosition = RadialPosition*sin((copyID+0.5)*90.*CLHEP::deg);
	}
	return G4ThreeVector(XPosition, YPosition, ZPosition);
}
