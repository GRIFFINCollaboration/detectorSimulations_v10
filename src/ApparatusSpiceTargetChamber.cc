#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "PrimaryGeneratorMessenger.hh"//PGM

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

G4int ApparatusSpiceTargetChamber::fTargetWheelTargetPos[6]={0,0,1,0,1,0};

//Z distance to Z=0 mm (actual target may be elsewhere)
G4double ApparatusSpiceTargetChamber::fTargetChamberDistanceFromTarget = 1*CLHEP::mm;	
G4double ApparatusSpiceTargetChamber::fTargetWheelFaceZ=0*CLHEP::mm;
// The whole magnet & chamber apparatus are placed in the world relative to the Z fTargetChamberDistanceFromTarget
// The target assembly is placed in the world relative to the Z fTargetWheelFaceZ
//
// The rods connecting the two have been set such that they match and an error statement prints if distance changed.

ApparatusSpiceTargetChamber::ApparatusSpiceTargetChamber(G4String Options)//parameter chooses which lens is in place.
{
	BuildMEL=false;
	Empty=false;
	Light=false;
	Bunny=false;
	Eff=false;
	Source=false;
	Col=false;

	G4cout<<G4endl<<"SPICE Command String : "<<Options;
	if(Options.contains("Empty")||Options.contains("EMPTY")||Options.contains("empty")){
		Empty=true;
		G4cout<<G4endl<<" Using Empty SPICE Lens Chamber";
	}else{
		if(Options.contains("Med")||Options.contains("med")||Options.contains("MED")||Options.contains("MEL")||Options.contains("Medium")||Options.contains("medium")){
			BuildMEL=true;
			G4cout<<G4endl<<" Building SPICE Medium Energy Lens";	
		}else{
			G4cout<<G4endl<<" Building SPICE Low Energy Lens";
		}

		if(Options.contains("Col")||Options.contains("COL")||Options.contains("col")||Options.contains("collimator")||Options.contains("Collimator")||Options.contains("COLLIMATOR")){
			Col=true;
			G4cout<<G4endl<<" ADDING NEW COLLIMATOR INSTALLED IN AUGUST 2017.";
		}
	}

	if(Options.contains("Light")||Options.contains("LIGHT")||Options.contains("light")||Options.contains("lite")||Options.contains("LITE")||Options.contains("Lite")){
		Light=true;
		G4cout<<G4endl<<" Building Only Crucial SPICE Components";
	}

	if(Options.contains("Eff")||Options.contains("EFF")||Options.contains("eff")||Options.contains("efficiency")||Options.contains("efficiency")||Options.contains("EFFICIENCY")){
		Eff=true;
		G4cout<<G4endl<<" Building the SPICE On-Wheel Efficiency Source Holder.";
	}else if(Options.contains("Source")||Options.contains("source")||Options.contains("SOURCE")){
		Source=true;
		G4cout<<G4endl<<" Building Z=0 NON-Efficiency Source Holder.";
	}else if(Options.contains("Bunny")||Options.contains("BUNNY")||Options.contains("bunny")||Options.contains("bun")||Options.contains("Bun")||Options.contains("BUN")){
		Bunny=true;
		G4cout<<G4endl<<" Building Bunny-ear Target Pedestals.";
	}

	G4cout<<G4endl;

	// Materials
	fMagnetChamberClampMaterial = "Peek";
	fDownstreamConeMaterial = "Aluminum"; 
	fTargetChamberMaterial = "Delrin"; 
	fElectroBoxMaterial = "Aluminum";
	fColdFingerMaterial = "Copper"; 
	fS3CableCaseMaterial = "Delrin";

	fMagnetMaterial = "NdFeB"; 
	fPhotonShieldLayerOneMaterial = "WTa"; 
	fPhotonShieldLayerTwoMaterial = "Tin"; 
	fPhotonShieldLayerThreeMaterial = "Copper";
	fPsClampMaterial = "Aluminum";
	fShieldCoverMaterial = "Kapton";
	fMagnetCoverMaterial = "Peek";
	fPSBoltMaterial = "Brass"; 

	fTargetWheelMaterial = "Aluminum"; 
	fTargetWheelGearMaterialA = "Delrin";
	fTargetWheelGearMaterialB = "Aluminum";
	fGearPlateMaterial = "Aluminum";
	fTargetMountPlateMaterial = "Peek";
	fGearStickMaterial = "Peek"; 

	fConicalCollimatorMaterial = "WTa";//11/8
	fXRayInsertMaterial="Copper";

	if(Source)fTargetMountPlateMaterial = "Aluminum";

	//-----------------------------
	// Dimensions of Target Chamber
	//-----------------------------


	fTargetChamberCylinderInnerRadius = 97*CLHEP::mm;
	fTargetChamberCylinderOuterRadius = 102*CLHEP::mm;
	fTargetChamberCylinderLength = 89*CLHEP::mm;
	fTargetChamberCylinderFlatRadius = 94.6*CLHEP::mm;
	fTargetChamberMELTabRad= 34.5*CLHEP::mm;
	fTargetChamberMELTabWidth= 11.8*CLHEP::mm;
	fTargetChamberMagnetIndent= 3*CLHEP::mm;

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

	// -------------------------
	// Dimensions of ElectroBox
	// -------------------------
	fElectroboxOuterBoxLength = 269*CLHEP::mm;
	fElectroboxInnerBoxLength = 243.5*CLHEP::mm;
	fElectroboxMidLength = 660*CLHEP::mm;
	fElectroboxLipLength = 21.5*CLHEP::mm;
	fElectroboxLipInnerRadius = 102*CLHEP::mm;
	fElectroboxZOffset = -90*CLHEP::mm;

	// -------------------------
	// Dimensions of ColdFinger
	// -------------------------
	fColdfingerThickness = 8*CLHEP::mm ; // CHECK
	fColdfingerLength = 150*CLHEP::mm ; // CHECK
	fColdfingerWidth = 200*CLHEP::mm  ; // CHECK
	fColdfingerZOffset = -130*CLHEP::mm ; // CHECK
	fColdfingerHoleRadius = 10*CLHEP::mm ; // CHECK

	// ----------------------
	// Dimensions of Beam Pipe
	// -----------------------
	fPipeInnerRadius = 5.0*CLHEP::mm; // 5
	fPipeOuterRadius = 16.0*CLHEP::mm;
	fPipeZLength = 30.0*CLHEP::mm;


	// --------------------------
	// Dimensions of Target Wheel
	// --------------------------

	fTargetWheelRotation=30*CLHEP::deg;
	fTargetWheelHoleSpaceAng=60*CLHEP::deg;
	fTargetWheelRadius = 40*CLHEP::mm;
	fTargetWheelThickness = 3.1*CLHEP::mm;
	fTargetWheelIndentRadius= 28.9*CLHEP::mm;

	fTargetWheelExtLength= 13.5*CLHEP::mm;
	fTargetWheelExtInnerRad1= 32*CLHEP::mm;
	fTargetWheelExtInnerRad2= 38*CLHEP::mm;
	fTargetWheelExtBaseRad= 48.6*CLHEP::mm;
	fTargetWheelExtTargLip= 3*CLHEP::mm;
	fTargetWheelExtBaseLip= 3.2*CLHEP::mm;

	fTargetOffsetRadius = 15*CLHEP::mm;

	// Mount Plate
	fTargetMountPlateRadius = 87*CLHEP::mm;
	fTargetMountPlateThickness = 5*CLHEP::mm;
	fTargetMountPlateCutDepth = 3.1*CLHEP::mm;
	fTargetMountPlateHoleRad = 40*CLHEP::mm;

	fMountPlateRodRadius=3.125*CLHEP::mm;
	fMountPlateRodPlaceR=79*CLHEP::mm;
	fMountPlateRodLength=14.5*CLHEP::mm;

	fTargetMountPlateZOffset = fTargetWheelFaceZ+fTargetWheelThickness+fTargetWheelExtLength-fTargetMountPlateCutDepth+fTargetMountPlateThickness/2.;

	fMountPlateRodZoffset=fTargetMountPlateZOffset-fTargetMountPlateThickness/2.-fMountPlateRodLength/2;

	G4double AlternatefMountPlateRodZoffset=-fTargetChamberDistanceFromTarget+fMountPlateRodLength/2;	
	G4double delta=	abs((fMountPlateRodZoffset-AlternatefMountPlateRodZoffset)/CLHEP::mm);
	if(delta>0.1)G4cout<<G4endl<<"BEWARE DIFFERENCE IN Z=0 BETWEEN LENS AND TARGET WHEEL CONSTRUCTIONS "<<delta<<" mm"<<G4endl;

	if(Eff||Bunny){
		//Bunny ears baseplate
		fTargetWheelIndentDepth= 1.4*CLHEP::mm;
		fTargetWheelTargetHoleRad=9*CLHEP::mm;
	}else{
		// 		Evaporator targets Z=0
		fTargetWheelIndentDepth= 0.4*CLHEP::mm;
		fTargetWheelTargetHoleRad=5*CLHEP::mm;
	}
	//Old 5 hole "superman" baseplate
	// 	fTargetWheelTargetPos[1]=1;
	// 	fTargetWheelTargetPos[3]=1;
	// 	fTargetWheelTargetPos[5]=1;
	// 	fTargetWheelIndentDepth= 1.4*CLHEP::mm;
	// 	fTargetWheelTargetHoleRad=6*CLHEP::mm;


	//Collimator
	fTargColRadius1=1*CLHEP::mm;
	fTargColRadius2=2.5*CLHEP::mm;
	fTargColHoleHalfSep=5.13*CLHEP::mm;
	fTargColPlateAng=76*CLHEP::deg;
	fTargColPlateThick=1*CLHEP::mm;
	fTargColPlateFlat=10.8*CLHEP::mm;
	fTargColPlateRad=19.8*CLHEP::mm;

	fTargetWheelOffset = fTargetOffsetRadius;

	// Gears
	fFirstGearRadius = 12.7*CLHEP::mm;
	fFirstGearPlaneOffset = 80.25*CLHEP::mm;
	fSecondGearRadius = 15.875*CLHEP::mm;
	fSecondGearPlaneOffset = 51.675*CLHEP::mm;
	fThirdGearRadius = 50.8*CLHEP::mm;
	fThirdGearPlaneOffset = fTargetWheelOffset;
	fAllGearThickness = 3.175*CLHEP::mm;
	fAllGearZOffset =fTargetWheelFaceZ + fAllGearThickness/2.+ 6.5*CLHEP::mm;

	// Gear Plates
	fGearPlateOneRadius = 5*CLHEP::mm;
	fGearPlateTwoRadius = 6*CLHEP::mm;
	fGearPlateThickness = 3.866*CLHEP::mm;


	// Gear Stick
	fGearStickLength = 250*CLHEP::mm;
	fGearStickRadius = 3.18*CLHEP::mm;

	// Targets
	fBunnyTargetThickness=0.5*CLHEP::mm;
	fBunnyRodLength=7.5*CLHEP::mm;

	fEvapTargetThickness=0.5*CLHEP::mm;
	fEvapTargetRadius=9.0*CLHEP::mm;
	fEvapTargetInner=4.0*CLHEP::mm;
	fEvapTargetTabWidth=14.0*CLHEP::mm;
	fEvapTargetTabExtra=0.5*CLHEP::mm;
	fEvapTargetTabIndent=1.6*CLHEP::mm;
	fEvapTargetScrewRad=2.0*CLHEP::mm;
	fEvapTargetScrewHeight=1.6*CLHEP::mm;

	// ------------------------------------------
	// Dimensions of Magnet Cover
	// ------------------------------------------
	fPsClampThicknessTop =1.5*CLHEP::mm;
	fPsClampThicknessBottom =4*mm;
	fPsClampThicknessCut= 1*CLHEP::mm;
	fPsClampThicknessInner = 3*CLHEP::mm;
	fPsClampThicknessOuter = -4.9*CLHEP::mm;
	fPsClampThicknessFace = 1*CLHEP::mm;

	// Note both Magnets & Magnet Cover overlap with Target Chamber
	// This was much easier than making the relevant cut-outs of Target Chamber & Magnet Cover
	// The exposed not overlapped geometries are correct as a result.

	//---------------------
	// Dimensions of Magnet
	//---------------------

	fPlateOneThickness = 3*CLHEP::mm;
	fPlateOneLength = 75*CLHEP::mm; 
	fPlateOneHeight = 50*CLHEP::mm; // z-component
	fPlateOneLowerHeight = 20*CLHEP::mm;

	if(BuildMEL){
		fPlateTwoThickness = 3.4*CLHEP::mm; // LEL is 5.0, MEL is 3.4
		fPlateTwoLength = 55.0*CLHEP::mm; // LEL is 30.0, MEL is 55.0
	} else {
		fPlateTwoThickness = 5.0*CLHEP::mm; // LEL is 5.0, MEL is 3.4
		fPlateTwoLength = 30.0*CLHEP::mm; // LEL is 30.0, MEL is 55.0
	}
	fPlateTwoHeight = fPlateOneHeight;

	fCuttingBoxAngle = 60*CLHEP::deg;

	//Z direction of magnet nearest edge to Z=0 mm (actual target may be elsewhere)
	fMagDistanceFromTarget =fTargetChamberDistanceFromTarget+fTargetChamberCylinderLipThickness-fTargetChamberMagnetIndent;

	fMagnetRadiusClampOffset = 2.4*CLHEP::mm;

	//Radial distance to the outer edge of the magnet
	fPlateOneEdgeX = (fTargetChamberCylinderFlatRadius - fPlateOneLength -fMagnetRadiusClampOffset); 

	// ------------------------------------
	// Dimensions of Magnet Clamp (Chamber)
	// ------------------------------------
	fMclampChamberThickness = fPlateOneThickness+fPlateTwoThickness*2+7.1*CLHEP::mm;
	fMclampChamberLength = 6.5*CLHEP::mm;
	fMclampChamberHeight = fTargetChamberCylinderLength-fTargetChamberCylinderLipThickness;

	//----------------------------
	// Dimensions of Photon Shield
	//----------------------------
	fPhotonShieldFrontRadius = 13.889*CLHEP::mm; 
	fPhotonShieldBackRadius = 26.031*CLHEP::mm;
	fPhotonShieldLength = 30.*CLHEP::mm;
	fPhotonShieldInnerRadius = 5.*CLHEP::mm;
	fPhotonShieldBackFacePos = -(fMagDistanceFromTarget+fPlateOneHeight+fPsClampThicknessBottom);

	fPhotonShieldLayerOneThickness = 25*CLHEP::mm;
	fPhotonShieldLayerTwoThickness = 4*CLHEP::mm;
	fPhotonShieldLayerThreeThickness = fPhotonShieldLength-fPhotonShieldLayerOneThickness-fPhotonShieldLayerTwoThickness;

	fPsClampTabHeight = fPhotonShieldLength;

	// ----------------------------------
	// Dimensions of Photon Shield Clamps
	// ----------------------------------
	fPsTargetClampRadius = 11.8*CLHEP::mm;
	fPsTargetClampThickness = 2*CLHEP::mm;
	fPsTargetClampExtension = 17.*CLHEP::mm;
	fPsTargetClampExtWidth = 4.8*CLHEP::mm;
	fPsDetClampRadius = 21.5*CLHEP::mm;
	fPsDetClampThickness = 3*CLHEP::mm;

	// ---------------------------------
	// Dimensions of Photon Shield Bolts
	// ---------------------------------
	fPsBoltHoleRadius = 1.7*CLHEP::mm;
	fPsBoltPlaceRadius = 9*CLHEP::mm;

	fPsBoltLength=50*CLHEP::mm;
	fPsBoltRadius=1.5*CLHEP::mm;
	fPsBoltBackShield=10.5*CLHEP::mm;
	fPsNutLength=2.4*CLHEP::mm;
	fPsNutRad=3.15*CLHEP::mm;

	// ------------------------------
	// Dimensions of Kapton Coatings
	// ------------------------------
	fShieldCoatingThickness = 0.2*CLHEP::mm;

	// --------------------------
	// Dimensions of Source 
	// --------------------------

	fEffSourceHolderHeight=8.7*CLHEP::mm;
	fEffSourceHolderRad=15.875*CLHEP::mm;
	fEffSourceHolderInner=9.475*CLHEP::mm;
	fEffSourceHoldPocket=3.4*CLHEP::mm;
	fEffSourceHoldPocketRad=12.75*CLHEP::mm;
	fEffSourceScrewRad=1.9*CLHEP::mm;
	fEffSourceScrewHeight=1.9*CLHEP::mm;

	fMEezagSourceHeight=3.18*CLHEP::mm;
	fMEezagSourceRad=12.2*CLHEP::mm;
	fMEezagSourceWall=0.5*CLHEP::mm;
	fMEezagSourceRim=3.*CLHEP::mm;


	G4Material* MatAc = G4Material::GetMaterial("Acrylic");
	G4double fMEezagSourceAcrylic=0.2*CLHEP::mg / CLHEP::cm2;
	fMEezagSourceAcrylicLength=fMEezagSourceAcrylic/MatAc->GetDensity();

	fBracketSideLength = 19.5*mm;
	fBracketBackLength = 30.7*mm;
	fBracketDepth = 3.2*mm;
	fBracketBackWidth = 3.5*mm;
	fBracketSideWidth = 4.0*mm;
	fBracketoffset = 2.95*mm;
	fBracketAngle = 180.*CLHEP::deg;
	fBracketScrewRad = 2.9*mm;
	fBracketScrewHeight = 2.5*mm;

	if(Eff){
		fMEezagSourceZ=fTargetWheelFaceZ-fEffSourceHolderHeight+fMEezagSourceHeight/2.;
		// 	fMEezagSourceZ=fTargetWheelFaceZ-fEffSourceHolderHeight+fEffSourceHoldPocket-fMEezagSourceHeight/2.;
	}else{		
		fMEezagSourceZ=-fTargetChamberDistanceFromTarget+fMountPlateRodLength-11*CLHEP::mm-fMEezagSourceHeight/2.;
	}
	fMEezagSourceActiveZ=fMEezagSourceZ-fMEezagSourceHeight/2.+fMEezagSourceWall;

	// --------------------------------
	// Dimensions of Conical Collimator
	// --------------------------------
	fConicalInnerRad1=2.9*CLHEP::mm;
	fConicalInnerRad2=3.15*CLHEP::mm;
	fConicalOuterRad1=4.9*CLHEP::mm;
	fConicalOuterRad2=4.5*CLHEP::mm;
	fConicalLength1=15.5*CLHEP::mm;
	fConicalLength2=13*CLHEP::mm;
	fConicalLength3=12*CLHEP::mm;
	fConicalTabRad=6.9*CLHEP::mm;
	fConicalTabHeight=2*CLHEP::mm;
	fConicalTabWidth=4*CLHEP::mm;


	// --------------------------
	// Dimensions of X-ray Insert
	// --------------------------

	fInsertCylOuterRadius=8.5*CLHEP::mm;
	fInsertCylInnerRadius=4.55*CLHEP::mm;
	fInsertCylZLength=3.2*CLHEP::mm;

	fInsertConeZLength=3.0*CLHEP::mm;
	fInsertConeR1=3.15*CLHEP::mm;
	fInsertConeR2=3.68*CLHEP::mm;

	fInsertEdgeAngle=37.275*CLHEP::deg;
	fInsertEdgeZLength=3.4*CLHEP::mm;

	fInsertHoleRadius=2*CLHEP::mm;

	// -------------------------------------
	// individual offsets for visualisation
	// -------------------------------------
	fFrontDomeOffset = 0*CLHEP::mm;
	fMiddleRingOffset = 0*CLHEP::mm;
	fBackDetAndAlBoxOffset = 0*CLHEP::mm;

	if(Eff||Source){
		G4cout<<G4endl<<"Using Source. GUN MUST BE AT Z="<<fMEezagSourceActiveZ/mm<<" mm";
	}else  if(Bunny){
		G4cout<<G4endl<<"Bunny Frames. Targets & Gun should be from -ve from Z="<<(fTargetWheelFaceZ-fBunnyRodLength-fBunnyTargetThickness)/mm<<" mm (or +ve from Z="<<(fTargetWheelFaceZ-fBunnyRodLength)/mm<<" mm)";
	}else{
		G4cout<<G4endl<<"Z=0 Target Frames. Targets & Gun should be from -ve from Z="<<(fTargetWheelFaceZ-fEvapTargetThickness)/mm<<" mm";
	}
	G4cout<<G4endl;


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
	delete fShieldBoltLog;
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

	if(!Light){
		BuildTargetChamberFrontRing();
		PlaceTargetChamberFrontRing(expHallLog);
		BuildElectroBox();
		PlaceElectroBox(expHallLog);	
		BuildColdFinger(); 
		PlaceColdFinger(expHallLog); 
		BuildTargetChamberSphere();
		PlaceTargetChamberSphere(expHallLog);		
	}

	BuildS3CableHolder();
	PlaceS3CableHolder(expHallLog);

	BuildTargetChamberCylinderDownstream(BuildMEL);
	PlaceTargetChamberCylinderDownstream(expHallLog); 

	if(!Empty){
		BuildMagnetClampChamber();
		BuildCollectorMagnetandCover();
		BuildPhotonShieldClampBolts();
		for(G4int copyID=0; copyID<4; copyID++)    {
			PlaceMagnetAndCover(copyID, expHallLog);
			PlacePhotonShieldClampBolts(copyID, expHallLog);
		}

		BuildPhotonShield();
		PlacePhotonShield(expHallLog); 

		BuildShieldCovering();
		PlaceShieldCovering(expHallLog);

		BuildPhotonShieldClamps();
		PlacePhotonShieldClamps(expHallLog);

		if(Col){
			BuildConicalCollimator();
			PlaceConicalCollimator(expHallLog);
			BuildXrayInsert();	
			PlaceXrayInsert(expHallLog);
		}
	}

	BuildTargetWheel();
	PlaceTargetWheel(expHallLog); 	
	if(!Light){
		BuildGearStick();
		PlaceGearStick(expHallLog); 

		BuildTargetMountPlate();
		PlaceTargetMountPlate(expHallLog); 		

		if(!Source){
			BuildTargetWheelGears();
			PlaceTargetWheelGears(expHallLog); 

			BuildTargetWheelGearPlates();
			PlaceTargetWheelGearPlates(expHallLog); 

			BuildCollimator();		
			PlaceCollimator(expHallLog); 
		}
	}

	if(Bunny){
		BuildBunnyTarget();
		PlaceBunnyTargets(expHallLog);
	}else if(Eff){
		BuildPlaceEfficiencySource(expHallLog);
	}else if(Source){
		BuildPlaceSource(expHallLog);
	}else{
		//Normal Z=0 target
		BuildEvapTarget();
		PlaceEvapTargets(expHallLog);
	}
} // end Build

////////////////////////////////////////////////////
// methods used to build and place the components:
////////////////////////////////////////////////////


void ApparatusSpiceTargetChamber::BuildBunnyTarget() {
	// Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(0.0, 0.9, 0.1));
	sVisAtt->SetVisibility(true);	

	// Dimensions
	G4double sDimOuterRadius = 7.00*CLHEP::mm;
	G4double sDimInnerRadius = 4.00*CLHEP::mm;
	G4double sDimHalfThick = fBunnyTargetThickness/2.;

	// Shapes
	G4Tubs* sTubs = new G4Tubs("sTubs", sDimInnerRadius, sDimOuterRadius, sDimHalfThick, 0, 360.*CLHEP::deg);
	G4Tubs* sLimit = new G4Tubs("sLimit", 20.648, 25.0, 2*sDimHalfThick, 0, 360.*CLHEP::deg);

	G4ThreeVector sTrans(15.*CLHEP::mm, 0., 0.);
	G4SubtractionSolid* sSub = new G4SubtractionSolid("sSub", sTubs, sLimit, 0, sTrans);

	G4Box* sBox = new G4Box("sBox", 2.0*CLHEP::mm, 1.0*CLHEP::mm, sDimHalfThick);
	G4Tubs* sRound = new G4Tubs("sTubs", 0, 1.0*CLHEP::mm, sDimHalfThick, -90.*CLHEP::deg, 180.*CLHEP::deg);


	G4UnionSolid* sRoundBox = new G4UnionSolid("sRoundBox", sBox, sRound, 0, G4ThreeVector(2.0*CLHEP::mm,0,0));

	G4ThreeVector RoundP(8.0*CLHEP::mm, 0., 0.);
	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(45.*CLHEP::deg);
	RoundP.rotateZ(-45.*CLHEP::deg);

	G4UnionSolid* sUnion = new G4UnionSolid("sUnion", sSub, sRoundBox, sRotate, RoundP);

	sRotate->rotateZ(-90.*CLHEP::deg);
	RoundP.rotateZ(90.*CLHEP::deg);
	G4UnionSolid* sUnion2 = new G4UnionSolid("sUnion2", sUnion, sRoundBox, sRotate, RoundP);

	// Logical
	G4Material* sMaterial = G4Material::GetMaterial("Aluminum");
	fFrameLogical = new G4LogicalVolume(sUnion2, sMaterial, "TargetFrame", 0, 0, 0);
	fFrameLogical->SetVisAttributes(sVisAtt);

	//Build Rods

	// Visualisation
	sVisAtt = new G4VisAttributes(G4Colour(AL_COL));
	sVisAtt->SetVisibility(true);

	// Dimensions
	sDimHalfThick = fBunnyRodLength/2;

	// Shapes
	G4Tubs* sRodTubs = new G4Tubs("sRodTubs", 0, 1.0*CLHEP::mm, sDimHalfThick, 0, 360*CLHEP::deg);

	// Logical
	sMaterial = G4Material::GetMaterial(fTargetWheelMaterial);
	fRodLogical = new G4LogicalVolume(sRodTubs, sMaterial, "TargetRodsLog", 0, 0, 0);
	fRodLogical->SetVisAttributes(sVisAtt);
}

void ApparatusSpiceTargetChamber::BuildEvapTarget(){
	// Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(0.0, 0.9, 0.1));
	sVisAtt->SetVisibility(true);	

	// Dimensions
	G4double sDimHalfThick = fEvapTargetThickness/2.;
	G4double sTabHeightHalf = (fEvapTargetTabExtra+fEvapTargetRadius-fEvapTargetInner)/2.;

	// Shapes
	G4Tubs* sTubs = new G4Tubs("sTubs", fEvapTargetInner, fEvapTargetRadius, sDimHalfThick, 0, 360.*CLHEP::deg);
	G4Box* sTab = new G4Box("sTab",fEvapTargetTabWidth/2.,sTabHeightHalf, sDimHalfThick);

	G4ThreeVector sTrans(0., -(sTabHeightHalf+fEvapTargetInner), 0.);
	G4UnionSolid* sSub = new G4UnionSolid("sSub", sTubs, sTab, 0, sTrans);


	G4Box* sBox = new G4Box("sBox", fEvapTargetTabIndent,fEvapTargetTabExtra, fEvapTargetThickness);
	G4Tubs* sRound = new G4Tubs("sTubs", 0, fEvapTargetTabIndent, fEvapTargetThickness, 0, 180.*CLHEP::deg);

	G4UnionSolid* sRoundBox = new G4UnionSolid("sRoundBox", sRound, sBox, 0, G4ThreeVector(0,-fEvapTargetTabExtra,0));

	G4SubtractionSolid* sTarg = new G4SubtractionSolid("sTarg", sSub, sRoundBox, 0, G4ThreeVector(0,-fEvapTargetRadius,0));

	// Logical
	G4Material* sMaterial = G4Material::GetMaterial("Aluminum");
	fFrameLogical = new G4LogicalVolume(sTarg, sMaterial, "TargetFrame", 0, 0, 0);
	fFrameLogical->SetVisAttributes(sVisAtt);

	//Build Rods

	// Visualisation
	sVisAtt = new G4VisAttributes(G4Colour(CU_COL));
	sVisAtt->SetVisibility(true);
	// Shapes
	G4Tubs* sRodTubs = new G4Tubs("sRodTubs", 0, fEvapTargetScrewRad, fEvapTargetScrewHeight/2., 0, 360*CLHEP::deg);

	// Logical
	sMaterial = G4Material::GetMaterial(fPSBoltMaterial);
	fRodLogical = new G4LogicalVolume(sRodTubs, sMaterial, "TargetRodsLog", 0, 0, 0);
	fRodLogical->SetVisAttributes(sVisAtt);
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

void ApparatusSpiceTargetChamber::BuildTargetChamberCylinderDownstream(bool MEL)
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
	// Target Lip Cutting Box
	G4double FlatHalfZ = length;
	G4double FlatHalfY = sqrt(pow(fTargetChamberCylinderInnerRadius,2)-pow(fTargetChamberCylinderFlatRadius,2));
	G4double FlatHalfX = (fTargetChamberCylinderInnerRadius-fTargetChamberCylinderFlatRadius)/2.;	

	G4double MELTabHalfZ = TargetLipLength;
	G4double MELTabHalfY = fTargetChamberMELTabWidth/2.;
	G4double MELTabHalfX = fTargetChamberMELTabWidth;


	// ** Shapes
	G4Tubs* MainCylinder = new G4Tubs("MainCylinder",InnerRadius,OuterRadius,length,0,360*CLHEP::deg);
	G4Tubs* TargetLipCylinder = new G4Tubs("TargetLipCylinder", TargetLipInnerRadius, OuterRadius, TargetLipLength, 0, 360*CLHEP::deg);
	G4Tubs* DetectorLipCylinder = new G4Tubs("DetectorLipCylinder", OuterRadius, DetectorLipOuterRadius, DetectorLipLength, 0, 360*CLHEP::deg);
	G4Box* TargetLipCuttingBox = new G4Box("TargetLipCuttingBox", LipCuttingHalfX, LipCuttingHalfY, LipCuttingHalfZ);
	G4Box* CylinderFlatBox = new G4Box("CylinderFlatBox", FlatHalfX, FlatHalfY, FlatHalfZ);
	G4Box* MELTabBox = new G4Box("MELTabBox", MELTabHalfX, MELTabHalfY, MELTabHalfZ);

	// Union


	G4VSolid* Lip=TargetLipCylinder;
	G4VSolid* MainCylinderFlats=MainCylinder;
	for(int i=0;i<4;i++){
		G4double ang=(45.0+90*i)*CLHEP::deg ;
		G4RotationMatrix* rotate = new G4RotationMatrix(0, 0,ang);
		G4ThreeVector trans(fTargetChamberCylinderFlatRadius+FlatHalfX,0,0);
		trans.rotateZ(ang);
		MainCylinderFlats = new G4UnionSolid("MainCylinderFlats", MainCylinderFlats, CylinderFlatBox, rotate, trans);

		if(MEL){
			G4ThreeVector transs(fTargetChamberMELTabRad+MELTabHalfX,0,0);
			transs.rotateZ(ang);
			Lip = new G4UnionSolid("TargetDetChamber", Lip, MELTabBox, rotate, transs);
		}
	}

	G4ThreeVector trans(0,0, length-TargetLipLength);
	G4UnionSolid* TargetDetCylinderPre = new G4UnionSolid("TargetDetCylinderPre", MainCylinderFlats, Lip, 0, trans);
	trans.setZ(-(length-DetectorLipLength));
	G4UnionSolid* TargetDetCylinder = new G4UnionSolid("TargetDetCylinder", TargetDetCylinderPre, DetectorLipCylinder, 0, trans);
	trans.setX(fTargetChamberCylinderHoleDistanceX);
	trans.setY(-fTargetChamberCylinderHoleDistanceY);
	trans.setZ(length-TargetLipLength);
	G4RotationMatrix* rotate = new G4RotationMatrix(-111.8*CLHEP::deg , 0, 0);
	G4VSolid* TargetDetChamber = new G4SubtractionSolid("TargetDetChamber", TargetDetCylinder, TargetLipCuttingBox, rotate, trans);

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
	G4double WheelHalfThickPlate = (fTargetWheelThickness-fTargetWheelIndentDepth)/2.;
	G4double WheelHalfThickRing1 = (fTargetWheelIndentDepth)/2.;
	G4double WheelHalfThickRing2 = (fTargetWheelExtTargLip)/2.;
	G4double WheelHalfThickRing3 = (fTargetWheelExtLength-fTargetWheelExtTargLip-fTargetWheelExtBaseLip)/2.;
	G4double WheelHalfThickRing4 = (fTargetWheelExtBaseLip)/2.;

	// ** Shapes
	G4Tubs* TargetWheelBase = new G4Tubs("TargetWheelBase", 0, fTargetWheelRadius, WheelHalfThickPlate, 0, 360*CLHEP::deg);
	G4Tubs* targethole = new G4Tubs("targethole", 0, fTargetWheelTargetHoleRad, fTargetWheelThickness, 0, 360*CLHEP::deg);

	G4Tubs* ring1 = new G4Tubs("ring1", fTargetWheelIndentRadius, fTargetWheelRadius, WheelHalfThickRing1, 0, 360*CLHEP::deg);
	G4Tubs* ring2 = new G4Tubs("ring2", fTargetWheelExtInnerRad1, fTargetWheelRadius, WheelHalfThickRing2, 0, 360*CLHEP::deg);
	G4Tubs* ring3 = new G4Tubs("ring3", fTargetWheelExtInnerRad2, fTargetWheelRadius, WheelHalfThickRing3, 0, 360*CLHEP::deg);
	G4Tubs* ring4 = new G4Tubs("ring4", fTargetWheelExtInnerRad2, fTargetWheelExtBaseRad, WheelHalfThickRing4, 0, 360*CLHEP::deg);

	G4VSolid* TargetWheel=TargetWheelBase;
	for(int i=0;i<6;i++){
		if(!fTargetWheelTargetPos[i])continue;
		G4ThreeVector trans(0, -fTargetOffsetRadius, 0);
		trans.rotateZ(fTargetWheelHoleSpaceAng*i);
		TargetWheel= new G4SubtractionSolid("TargetWheel", TargetWheel, targethole, 0, trans);
	}

	G4double colholerad=fTargColRadius2+1.0*mm;
	G4Box* colrect = new G4Box("colrect",fTargColHoleHalfSep,colholerad, fTargetWheelThickness);
	G4Tubs* colcurv = new G4Tubs("colcurv", 0, colholerad, fTargetWheelThickness, 0, 360*CLHEP::deg);
	G4UnionSolid* collimator = new G4UnionSolid("collimator", colrect, colcurv, 0, G4ThreeVector(fTargColHoleHalfSep+0.1*mm,0,0));
	collimator = new G4UnionSolid("collimator", collimator, colcurv, 0, G4ThreeVector(-(fTargColHoleHalfSep+0.1*mm),0,0));

	G4double collimatoroffser=sqrt(pow(fTargetOffsetRadius,2)-pow(fTargColHoleHalfSep,2));
	TargetWheel= new G4SubtractionSolid("TargetWheel", TargetWheel, collimator, 0, G4ThreeVector(0,-collimatoroffser,0));


	G4double Zoffset=WheelHalfThickPlate+WheelHalfThickRing1;
	TargetWheel = new G4UnionSolid("BoltsTwo", TargetWheel, ring1, 0, G4ThreeVector(0,0,Zoffset));
	Zoffset+=WheelHalfThickRing1+WheelHalfThickRing2;
	TargetWheel = new G4UnionSolid("BoltsTwo", TargetWheel, ring2, 0, G4ThreeVector(0,0,Zoffset));
	Zoffset+=WheelHalfThickRing2+WheelHalfThickRing3;
	TargetWheel = new G4UnionSolid("BoltsTwo", TargetWheel, ring3, 0, G4ThreeVector(0,0,Zoffset));
	Zoffset+=WheelHalfThickRing3+WheelHalfThickRing4;
	TargetWheel = new G4UnionSolid("BoltsTwo", TargetWheel, ring4, 0, G4ThreeVector(0,0,Zoffset));

	// ** Logical
	G4Material* targetWheelMaterial = G4Material::GetMaterial(fTargetWheelMaterial); //problem
	fTargetWheelLog = new G4LogicalVolume(TargetWheel, targetWheelMaterial, "TargetWheelLog", 0, 0, 0);
	fTargetWheelLog->SetVisAttributes(VisAtt);
} // end BuildTargetWheel()


void ApparatusSpiceTargetChamber::BuildCollimator() {
	// Visualisation
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	sVisAtt->SetVisibility(true);	

	// Dimensions
	G4Tubs* platearc = new G4Tubs("platearc", 0, fTargColPlateRad, fTargColPlateThick/2., (270*CLHEP::deg)-(fTargColPlateAng/2.),fTargColPlateAng);


	G4double flatdist=fTargColPlateFlat/(2*tan(fTargColPlateAng/2.));
	G4Box* sRectangle = new G4Box("sRectangle", flatdist*2, flatdist,fTargColPlateThick);

	G4SubtractionSolid* plate = new G4SubtractionSolid("plate", platearc, sRectangle, 0, G4ThreeVector(0,0,0));

	G4double collimatoroffser=sqrt(pow(fTargetOffsetRadius,2)-pow(fTargColHoleHalfSep,2));
	G4Tubs* sTubs1 = new G4Tubs("sTubs1", 0., fTargColRadius1, fTargColPlateThick, 0., 360.*CLHEP::deg);
	G4SubtractionSolid* sSub = new G4SubtractionSolid("sSub", plate, sTubs1, 0, G4ThreeVector(-fTargColHoleHalfSep,-collimatoroffser,0));

	G4Tubs* sTubs5 = new G4Tubs("sTubs5", 0., fTargColRadius2,fTargColPlateThick, 0., 360.*CLHEP::deg);
	G4SubtractionSolid* sSub2 = new G4SubtractionSolid("sSub2", sSub, sTubs5, 0, G4ThreeVector(fTargColHoleHalfSep,-collimatoroffser,0));	


	// Logical
	G4Material* CollimatorMaterial = G4Material::GetMaterial("Titanium"); 
	fCollimLogical = new G4LogicalVolume(sSub2, CollimatorMaterial, "collimator", 0, 0, 0);
	fCollimLogical->SetVisAttributes(sVisAtt);
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

	// ** Shapes
	G4Tubs* MountPlatePre = new G4Tubs("MountPlatePre", 0, fTargetMountPlateRadius, fTargetMountPlateThickness/2., 0, 360*CLHEP::deg);
	G4Tubs* TargetWheel = new G4Tubs("TargetWheel", 0, fTargetMountPlateHoleRad, fTargetMountPlateThickness, 0, 360*CLHEP::deg);
	G4Tubs* Lip = new G4Tubs("TargetWheel", 0, fTargetWheelExtBaseRad,fTargetMountPlateCutDepth, 0, 360*CLHEP::deg);

	G4ThreeVector move(-fTargetWheelOffset, 0, 0);
	G4SubtractionSolid* MountPlate = new G4SubtractionSolid("MountPlate", MountPlatePre, TargetWheel, 0, move);
	G4ThreeVector moveB(-fTargetWheelOffset, 0, -fTargetMountPlateThickness/2.);
	MountPlate = new G4SubtractionSolid("MountPlateB", MountPlate, Lip, 0, moveB);

	// ** Logical
	G4Material* targetMountPlateMaterial = G4Material::GetMaterial(fTargetMountPlateMaterial);
	fTargetMountPlateLog = new G4LogicalVolume(MountPlate, targetMountPlateMaterial, "TargetMountPlateLog", 0, 0, 0);
	fTargetMountPlateLog->SetVisAttributes(VisAtt);

	G4Tubs* MountRod = new G4Tubs("MountRod", 0, fMountPlateRodRadius, fMountPlateRodLength/2., 0, 360*CLHEP::deg);
	fTargetMountRodLog = new G4LogicalVolume(MountRod, targetMountPlateMaterial, "MountRodLog", 0, 0, 0);
	fTargetMountRodLog->SetVisAttributes(VisAtt);	

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
	G4double TotalLength = fPhotonShieldLength;
	G4double HalfLength = TotalLength/2.;

	G4double BackOuterRadius = fPhotonShieldBackRadius;
	G4double FrontOuterRadius = fPhotonShieldFrontRadius;

	G4Cons* PhotonShield = new G4Cons("PhotonShield", InnerRadius, BackOuterRadius, 
			InnerRadius, FrontOuterRadius,	HalfLength, 0, 2*CLHEP::pi);

	// ** Build Magnet Clamp
	G4double ClampHalfThickness = fPlateOneThickness/2. + fPsClampThicknessFace;
	G4Box* MagnetClamp = new G4Box("MagnetClamp", fPsClampTabHeight/2., ClampHalfThickness, fPhotonShieldLength);


	// ** Cut Plate One from Photon Shield
	G4RotationMatrix* cutRot = new G4RotationMatrix; 
	G4double transRadial = fPlateOneEdgeX + fPsClampTabHeight/2. - 0.1*mm;
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
	G4double BoltHalfLength = fPhotonShieldLength/2.;
	G4double BoltRadius = fPsBoltHoleRadius;
	G4double BoltOffset = fPsBoltPlaceRadius / sqrt(2.);

	G4Tubs* ShieldBolt = new G4Tubs("ShieldBolt", 0, BoltRadius, BoltHalfLength, 0, 360*CLHEP::deg);

	// the four front (target end) bolts are created using a union solid then deducted from the previous volume
	G4double DistanceBetweenBolts = fPsBoltPlaceRadius * sqrt(2.);
	G4ThreeVector BoltDistance(DistanceBetweenBolts, 0, 0);
	G4UnionSolid* BoltsTwo = new G4UnionSolid("BoltsTwo", ShieldBolt, ShieldBolt, 0, BoltDistance);
	G4ThreeVector BoltDist2(0, DistanceBetweenBolts, 0);
	G4UnionSolid* Bolts = new G4UnionSolid("Bolts", BoltsTwo, BoltsTwo, 0, BoltDist2);
	G4ThreeVector BoltMove(-BoltOffset, -BoltOffset, 0);
	G4SubtractionSolid* PhotonShieldTemp_6 = new G4SubtractionSolid("PhotonShield", PhotonShieldTemp_4, Bolts, 0, BoltMove);


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
	G4double ClampHalfThickness = fPlateOneThickness/2. + fPsClampThicknessFace;
	G4Box* MagnetClamp = new G4Box("MagnetClamp", fPsClampTabHeight/2., ClampHalfThickness, fPhotonShieldLength);

	// -- Rotation & Translation
	G4RotationMatrix* cutRot = new G4RotationMatrix; 
	G4double transRadial = fPlateOneEdgeX + fPsClampTabHeight/2. - 0.1*mm;
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
	G4double ExtHalfWidth = fPsTargetClampExtWidth/2.;
	G4double ExtHalfLength = fPsTargetClampExtension-TargetEndRadius*cos(asin(ExtHalfWidth/TargetEndRadius));

	// ** Shapes
	G4VSolid* TargetEndClamp = new G4Tubs("TargetEndClampPre", InnerRadius, TargetEndRadius, TargetEndHalf, 0, 360*CLHEP::deg);
	G4Box* TargetClampExt = new G4Box("TargetClampExt", ExtHalfLength, ExtHalfWidth, TargetEndHalf);
	G4Tubs* TargetBolt = new G4Tubs("TargetBolt", 0, fPsBoltHoleRadius,TargetEndHalf+DetectorEndHalf, 0, 360*CLHEP::deg);
	G4VSolid* DetEndClampPre = new G4Tubs("DetEndClampPre", InnerRadius, DetectorEndRadius, DetectorEndHalf, 0, 360*CLHEP::deg);

	for(int i=0;i<4;i++){
		G4double Ang=(i*90)*deg;
		G4RotationMatrix* rotation=new G4RotationMatrix();
		rotation->rotateZ(Ang);
		G4ThreeVector move(fPsTargetClampExtension-ExtHalfLength, 0, 0);
		move.rotateZ(-Ang);
		TargetEndClamp = new G4UnionSolid("TargetClamp", TargetEndClamp, TargetClampExt,rotation,move);

		G4ThreeVector boltplace(fPsBoltPlaceRadius, 0, 0);
		boltplace.rotateZ(-Ang);
		TargetEndClamp = new G4SubtractionSolid("TargetClamp", TargetEndClamp, TargetBolt,0,boltplace);
		DetEndClampPre = new G4SubtractionSolid("DetEndClamp", DetEndClampPre, TargetBolt,0,boltplace);	
	}

	// ** Logical
	G4Material* psClampMaterial = G4Material::GetMaterial(fPsClampMaterial);
	fPsTargetClampLog = new G4LogicalVolume(TargetEndClamp, psClampMaterial, "PsTargetClampLog", 0, 0, 0);
	fPsTargetClampLog->SetVisAttributes(VisAtt);
	fPsDetectorClampLog = new G4LogicalVolume(DetEndClampPre, psClampMaterial, "PsDetectorClampLog", 0, 0, 0);
	fPsDetectorClampLog->SetVisAttributes(VisAtt);
} // end::BuildPhotonShieldClamps

void ApparatusSpiceTargetChamber::BuildPhotonShieldClampBolts() {

	// ** Visualisation
	G4VisAttributes* SmallVisAtt = new G4VisAttributes(G4Colour(CU_COL));
	SmallVisAtt->SetVisibility(true);

	// 	// ** Dimensions

	G4double BoltHalfLength = fPsBoltLength/2.;
	G4double NutHalfLength = fPsNutLength/2.;

	G4Tubs* BoltShaft = new G4Tubs("BoltShaft", 0, fPsBoltRadius, BoltHalfLength, 0, 360*CLHEP::deg);
	G4Tubs* Nut = new G4Tubs("Nut", 0, fPsNutRad, NutHalfLength, 0, 360*CLHEP::deg);
	G4ThreeVector BoltVecOne(0,0,-BoltHalfLength+fPsBoltBackShield-fPsDetClampThickness-NutHalfLength);
	G4ThreeVector BoltVecTwo(0,0,-BoltHalfLength+fPsBoltBackShield+fPhotonShieldLength+fPsTargetClampThickness+NutHalfLength);

	G4UnionSolid* TargetBolt1 = new G4UnionSolid("TargetBolt1", BoltShaft, Nut, 0, BoltVecOne);
	G4UnionSolid* TargetBolt2 = new G4UnionSolid("TargetBolt2", TargetBolt1, Nut, 0, BoltVecTwo);

	// ** Logical
	G4Material* smallBoltMaterial = G4Material::GetMaterial(fPSBoltMaterial);
	fShieldBoltLog = new G4LogicalVolume(TargetBolt2, smallBoltMaterial, "TargetBoltLog", 0, 0, 0);
	fShieldBoltLog->SetVisAttributes(SmallVisAtt);

} // end::BuildPhotonShieldClampBolts()

void ApparatusSpiceTargetChamber::BuildCollectorMagnetandCover()
{
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(NDFEB_COL));
	VisAtt->SetVisibility(true);

	// ** Build Plate One
	G4Box* PlateOnePreCut = new G4Box("PlateOnePreCut", fPlateOneLength/2., 
			fPlateOneThickness/2., fPlateOneHeight/2.);
	G4Box* PlateOneCutBox = new G4Box("PlateOneCutBox", fPlateOneLength/2., 
			fPlateOneThickness/2. + 1., fPlateOneHeight/2.);

	// (Cut box is rotated and moved, exact movement calculated below)
	G4double AngleToCorner = atan(fPlateOneLength/fPlateOneHeight);
	G4double HypotenuseToCorner = sqrt(pow(fPlateOneHeight/2., 2) + pow(fPlateOneLength/2., 2));
	G4double AngleDifference, TranslateCutX, TranslateCutZ;
	if(AngleToCorner < fCuttingBoxAngle)	{
		AngleDifference = fCuttingBoxAngle - AngleToCorner;
		TranslateCutX = (HypotenuseToCorner * sin(AngleDifference)) + fPlateOneLength/2.;
		TranslateCutZ = (HypotenuseToCorner * cos(AngleDifference)) - fPlateOneHeight/2.
			+ fPlateOneLowerHeight;
	}
	// 	else if(AngleToCorner > fCuttingBoxAngle)
	else
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

	BuildMagnetCovering(PlateCombo);
} // end:BuildCollectorMagnet()

void ApparatusSpiceTargetChamber::BuildMagnetCovering(G4VSolid* Magnet)
{
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(PEEK_COL));
	VisAtt->SetVisibility(true);

	G4double PlateOneHalfLength = (fPlateOneLength+fPsClampThicknessInner+fPsClampThicknessOuter)/2;
	G4double PlateOneHalfThickness  = (fPlateOneThickness/2.)+fPsClampThicknessFace;
	G4double PlateOneHalfHeight = (fPlateOneHeight+fPsClampThicknessBottom+fPsClampThicknessTop)/2; 

	// ** Build Plate One
	G4Box* PlateOnePreCut = new G4Box("PlateOnePreCut", PlateOneHalfLength, 
			PlateOneHalfThickness, PlateOneHalfHeight);
	G4Box* PlateOneCutBox =  new G4Box("PlateOnePreCut", PlateOneHalfLength, 
			PlateOneHalfThickness +1, PlateOneHalfHeight);
	G4Box* PlateOneTabBox =  new G4Box("PlateOnePreCut", fPsClampTabHeight/2., 
			PlateOneHalfThickness, fPsClampTabHeight/2.);

	G4double CutBottomCornerX=PlateOneHalfLength;
	G4double CutBottomCornerZ=PlateOneHalfHeight;
	CutBottomCornerZ-=(fPlateOneLowerHeight+cos(fCuttingBoxAngle)*fPsClampThicknessCut);

	G4ThreeVector CutBottomCorner(-CutBottomCornerX, 0,-CutBottomCornerZ);
	G4ThreeVector Offset(-PlateOneHalfHeight*sin(fCuttingBoxAngle),0,PlateOneHalfHeight*cos(fCuttingBoxAngle));
	G4ThreeVector OffsetSlide(PlateOneHalfHeight*0.5*cos(fCuttingBoxAngle),0,PlateOneHalfHeight*0.5*sin(fCuttingBoxAngle));	

	G4ThreeVector TranslateCutBox=CutBottomCorner+Offset+OffsetSlide;
	G4RotationMatrix* RotateCutBox = new G4RotationMatrix;
	RotateCutBox->rotateY(fCuttingBoxAngle);
	G4VSolid* PlateOne = new G4SubtractionSolid("PlateOne", PlateOnePreCut, PlateOneCutBox,
			RotateCutBox, TranslateCutBox);	

	G4ThreeVector TranslateTabBox(-PlateOneHalfLength+fPsClampTabHeight/2.,0,-PlateOneHalfHeight+fPsClampTabHeight/2.);
	G4VSolid* PlateOneTab = new G4UnionSolid("PlateOne", PlateOne, PlateOneTabBox,0, TranslateTabBox);		

	G4double PlateTwoHalfLength = (fPlateTwoLength+fPsClampThicknessFace+fPsClampThicknessOuter)/2;
	G4double PlateTwoHalfThickness  = fPlateTwoThickness+(fPlateOneThickness/2.)+fPsClampThicknessFace;
	G4double PlateTwoHalfHeight = PlateOneHalfHeight; 

	// ** Build Plate Two
	G4Box* PlateTwo = new G4Box("PlateTwo", PlateTwoHalfLength,PlateTwoHalfThickness,PlateTwoHalfHeight);

	// ** Combine Plates
	G4ThreeVector TranslatePlateTwo(PlateOneHalfLength - PlateTwoHalfLength, 0, 0);
	G4UnionSolid* PlateCombo = new G4UnionSolid("PlateCombo", PlateOneTab, PlateTwo, 0, TranslatePlateTwo);

	G4Material* magnetCoverMaterial = G4Material::GetMaterial(fMagnetCoverMaterial);

	// ** Logical Volume
	fMagnetCoverLog = new G4LogicalVolume(PlateCombo, magnetCoverMaterial, "MagnetCoverLog", 0, 0, 0);

	//ONLY for viewing internal thickness is correct
	// 	G4Box* xray = new G4Box("xray", 5*cm,5*cm,5*cm);
	// 	G4VSolid* XRAY = new G4SubtractionSolid("MagnetCover", PlateCombo, xray,0,G4ThreeVector(0,-5*cm,0));
	// 	fMagnetCoverLog = new G4LogicalVolume(XRAY, magnetCoverMaterial, "MagnetCoverLog", 0, 0, 0);

	// This is the way it should be done but it has some glitches, so have left overlap with magnet rather than cutout
	// 	G4double MagCoverOffsetX = -(fPsClampThicknessInner-fPsClampThicknessOuter)/2;
	// 	G4double MagCoverOffsetZ = -(fPsClampThicknessBottom-fPsClampThicknessTop)/2; 
	// 	G4ThreeVector MagOffset(MagCoverOffsetX,0,-MagCoverOffsetZ);
	// 	G4VSolid* MagnetCover = new G4SubtractionSolid("MagnetCover", PlateCombo, Magnet,0,MagOffset);
	// 	fMagnetCoverLog = new G4LogicalVolume(MagnetCover, magnetCoverMaterial, "MagnetCoverLog", 0, 0, 0);
	Magnet->GetName();//Just here to suppress warning;

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
		(fMclampChamberLength - (fTargetChamberCylinderFlatRadius - (fPlateOneLength + fPlateOneEdgeX)));
	G4ThreeVector trans(xtrans, 0, ClampHalfHeight - MagnetHalfHeight);
	G4SubtractionSolid* MagnetClamp = new G4SubtractionSolid("MagnetClamp", ClampMain, MagnetCut, 0, trans);

	// ** Logical
	G4Material* clampMaterial = G4Material::GetMaterial(fMagnetChamberClampMaterial);
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
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(DELRIN_COL));
	sVisAtt->SetVisibility(true);

	// Dimensions
	// Shapes
	// 	G4Box* sBox = new G4Box("sBox", 8.1*CLHEP::mm, 12.*CLHEP::mm, 40.5*CLHEP::mm);
	G4Box* sBox = new G4Box("sBox", 6.3*CLHEP::mm, 12.*CLHEP::mm, fTargetChamberCylinderLength/2.);

	// Logical
	G4Material* sMaterial = G4Material::GetMaterial(fS3CableCaseMaterial);
	fS3CaseLogical = new G4LogicalVolume(sBox, sMaterial, "s3CaseLog", 0, 0, 0);
	fS3CaseLogical->SetVisAttributes(sVisAtt);
}
void ApparatusSpiceTargetChamber::BuildConicalCollimator() {

	// ** shapes
	G4Tubs* solid_cyl1 = new G4Tubs("solid_cyl1", fConicalInnerRad1, fConicalOuterRad1, fConicalLength1/2., 0, 360*CLHEP::deg);
	G4Tubs* solid_cyl2 = new G4Tubs("solid_cyl2", fConicalInnerRad2, fConicalOuterRad1, fConicalLength2/2., 0, 360*CLHEP::deg);
	G4Tubs* solid_cyl3 = new G4Tubs("solid_cyl3", fConicalInnerRad2, fConicalOuterRad2, fConicalLength3/2., 0, 360*CLHEP::deg);

	G4Box* sBox = new G4Box("sBox",(fConicalTabRad-fConicalInnerRad1)/2.,fConicalTabWidth/2, fConicalTabHeight/2.);
	G4VSolid* tabson =solid_cyl1;
	for(int i=0;i<4;i++){
		G4double ang=i*90*CLHEP::deg;
		G4ThreeVector moveR((fConicalTabRad+fConicalInnerRad1)/2, 0, (fConicalLength1-fConicalTabHeight)/2.);
		moveR.rotateZ(ang);
		G4RotationMatrix* sRotate = new G4RotationMatrix;
		sRotate->rotateZ(ang);
		tabson = new G4UnionSolid("Mid+Outer", tabson, sBox, sRotate, moveR);
	}

	G4Tubs* cut_cyl1 = new G4Tubs("solid_cyl1", fConicalTabRad, fConicalTabRad*1.5, fConicalLength1, 0, 360*CLHEP::deg);
	G4SubtractionSolid* Fullfirstcyl = new G4SubtractionSolid("Fullfirstcyl", tabson, cut_cyl1, 0, G4ThreeVector());

	G4double Z2=-(fConicalLength1+fConicalLength2)/2.;
	G4double Z3=Z2-(fConicalLength2+fConicalLength3)/2.;

	G4UnionSolid* TwoCyl = new G4UnionSolid("Mid+Outer", Fullfirstcyl, solid_cyl2, 0, G4ThreeVector(0,0,Z2));
	G4UnionSolid* ThreeCyl = new G4UnionSolid("Mid+Outer", TwoCyl, solid_cyl3, 0, G4ThreeVector(0,0,Z3));

	// ** Logical
	G4Material* sMaterial = G4Material::GetMaterial(this->fConicalCollimatorMaterial);
	fConicalCollimatorLog = new G4LogicalVolume(ThreeCyl,sMaterial,"fConicalCollimatorLog"); 

}

void ApparatusSpiceTargetChamber::BuildXrayInsert(){
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(CU_COL));
	sVisAtt->SetVisibility(true);

	G4Cons* InnerCone = new G4Cons("InnerCone",fInsertConeR1, fInsertCylInnerRadius, fInsertConeR2, fInsertCylInnerRadius,fInsertConeZLength/2., 0, 2*CLHEP::pi);

	G4double wideto=fInsertCylOuterRadius-1*mm+fInsertCylZLength;
	G4Cons* OuterCut = new G4Cons("OuterCut",wideto-fInsertCylZLength, wideto, wideto, wideto,fInsertCylZLength/2.+0.01*mm, 0, 360*CLHEP::deg);

	// ** shapes
	G4Tubs* solid_insert_cyl = new G4Tubs("insert_cyl", fInsertCylInnerRadius, fInsertCylOuterRadius, fInsertCylZLength/2, 0, 360*CLHEP::deg);
	G4Tubs* solid_insert_edges = new G4Tubs("insert_edge", fInsertCylInnerRadius,fInsertCylOuterRadius,fInsertEdgeZLength/2., -fInsertEdgeAngle/2.,fInsertEdgeAngle);
	G4Tubs* solid_insert_holes = new G4Tubs("insert_hole", 0.,fInsertHoleRadius,fInsertConeZLength, 0, 360*CLHEP::deg);

	G4VSolid* topcyl =solid_insert_cyl;
	G4ThreeVector moveH(fInsertCylOuterRadius+0.5*mm, 0,0);
	moveH.rotateZ(45*CLHEP::deg);
	for(int i=0;i<4;i++){
		moveH.rotateZ(90*CLHEP::deg);
		topcyl = new G4SubtractionSolid("topcyl", topcyl, solid_insert_holes, 0, moveH);
	}

	G4SubtractionSolid* edgetrim = new G4SubtractionSolid("edgetrim", topcyl, OuterCut, 0, G4ThreeVector());
	G4UnionSolid* incone = new G4UnionSolid("incone", edgetrim, InnerCone, 0, G4ThreeVector(0,0,(fInsertConeZLength-fInsertCylZLength)/2.));


	G4VSolid* fullcyl =incone;
	G4ThreeVector moveR(0, 0, (fInsertCylZLength+fInsertEdgeZLength)/2.);
	G4RotationMatrix* sRotate1 = new G4RotationMatrix();
	for(int i=0;i<4;i++){
		sRotate1->rotateZ(90*CLHEP::deg);
		fullcyl = new G4UnionSolid("fullcyl", fullcyl, solid_insert_edges, sRotate1, moveR);
	}

	// ** Logical  
	G4Material* sMaterial = G4Material::GetMaterial(this->fXRayInsertMaterial);
	fXRayInsertLog = new G4LogicalVolume(fullcyl,sMaterial,"fXRayInsertLog");
	fXRayInsertLog->SetVisAttributes(sVisAtt);
}




void ApparatusSpiceTargetChamber::BuildPlaceEfficiencySource(G4LogicalVolume* expHallLog){

	//Build shapes for the source holder

	G4double base=fEffSourceHolderHeight-fEffSourceHoldPocket;
	G4Tubs* Base = new G4Tubs("Base", fEffSourceHolderInner,fEffSourceHolderRad, base/2., 0, 360*CLHEP::deg);
	G4Tubs* Pocket = new G4Tubs("insert_cyl", fEffSourceHoldPocketRad, fEffSourceHolderRad, fEffSourceHoldPocket/2., 0, 360*CLHEP::deg);
	G4UnionSolid* Teflon = new G4UnionSolid("Teflon", Base, Pocket, 0, G4ThreeVector(0,0,-fEffSourceHolderHeight/2.));
	G4Tubs* Screw = new G4Tubs("SourcePlate", 0, fEffSourceScrewRad, fEffSourceScrewHeight/2., 0, 360*CLHEP::deg);

	//Build logical and place the source holder
	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(PEEK_COL));
	sVisAtt->SetVisibility(true);
	G4Material* sMaterialT = G4Material::GetMaterial("Teflon");
	G4LogicalVolume* TLog = new G4LogicalVolume(Teflon,sMaterialT,"TLog"); 
	TLog->SetVisAttributes(sVisAtt);
	new G4PVPlacement(new G4RotationMatrix,G4ThreeVector(0,0,fTargetWheelFaceZ-base/2.), TLog, "SourceHolderEff", expHallLog, false, 0);

	//Build logical and place the source holder screws
	sVisAtt = new G4VisAttributes(G4Colour(CU_COL));
	G4Material* sMaterialSR = G4Material::GetMaterial("Brass");
	G4LogicalVolume* SScrewLog = new G4LogicalVolume(Screw,sMaterialSR,"SScrewLog"); 
	SScrewLog->SetVisAttributes(sVisAtt);
	G4double screwradplacement=((fEffSourceHoldPocketRad-fEffSourceHolderRad)/2.) + fEffSourceHolderRad;
	new G4PVPlacement(new G4RotationMatrix,G4ThreeVector(0,screwradplacement,fTargetWheelFaceZ-fEffSourceHolderHeight-fEffSourceScrewHeight/2.), SScrewLog, "SourceHolderScrew1", expHallLog, false, 0);
	new G4PVPlacement(new G4RotationMatrix,G4ThreeVector(0,-screwradplacement,fTargetWheelFaceZ-fEffSourceHolderHeight-fEffSourceScrewHeight/2.), SScrewLog, "SourceHolderScrew2", expHallLog, false, 0);


	BuildMEezagSource(expHallLog);
}

void  ApparatusSpiceTargetChamber::BuildMEezagSource(G4LogicalVolume* expHallLog){
	//Build shapes for the source frame

	G4Tubs* SourceWall = new G4Tubs("SourceWall", fMEezagSourceRad-fMEezagSourceWall,fMEezagSourceRad, fMEezagSourceHeight/2., 0, 360*CLHEP::deg);
	G4Tubs* SourceRim = new G4Tubs("SourceRim", fMEezagSourceRad-fMEezagSourceRim, fMEezagSourceRad, fMEezagSourceWall/2., 0, 360*CLHEP::deg);
	G4Tubs* SourcePlate = new G4Tubs("SourcePlate", 0, fMEezagSourceRad-0.1*mm, fMEezagSourceWall/2., 0, 360*CLHEP::deg);


	G4UnionSolid* WallRim = new G4UnionSolid("WallRim", SourceWall, SourceRim, 0, G4ThreeVector(0,0,fMEezagSourceWall/2.-fMEezagSourceHeight/2.));
	G4UnionSolid* WallRimPlate = new G4UnionSolid("WallRimPlate", WallRim, SourcePlate, 0, G4ThreeVector(0,0,fMEezagSourceWall*1.5-fMEezagSourceHeight/2.));

	G4VisAttributes* sVisAtt = new G4VisAttributes(G4Colour(AL_COL));
	sVisAtt->SetVisibility(true);
	G4Material* sMaterialsf = G4Material::GetMaterial("Aluminum");
	G4LogicalVolume* SFlog = new G4LogicalVolume(WallRimPlate,sMaterialsf,"SFlog");
	SFlog->SetVisAttributes(sVisAtt);	
	new G4PVPlacement(new G4RotationMatrix,G4ThreeVector(0,0,fMEezagSourceZ), SFlog, "SourceFrame", expHallLog, false, 0);

	G4Tubs* Acrylic = new G4Tubs("Acrylic", 0, fMEezagSourceRad-fMEezagSourceRim-0.5*mm, fMEezagSourceAcrylicLength/2., 0, 360*CLHEP::deg);

	G4double AcrylicZ=fMEezagSourceActiveZ-fMEezagSourceAcrylicLength/2.;
	G4Material* sMaterialAc = G4Material::GetMaterial("Acrylic");
	G4LogicalVolume* Acrlog = new G4LogicalVolume(Acrylic,sMaterialAc,"Acrlog");
	new G4PVPlacement(new G4RotationMatrix,G4ThreeVector(0,0,AcrylicZ), Acrlog, "SourceAcrylic", expHallLog, false, 0);


}
void  ApparatusSpiceTargetChamber::BuildPlaceSource(G4LogicalVolume* expHallLog){
	//Build shapes for the source frame

	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
	vis_att->SetVisibility(true); 


	G4double inner=fBracketBackLength/2.-fBracketSideWidth;

	G4Box* OverBracket = new G4Box("Bracket", fBracketSideLength/2., fBracketBackLength/2., fBracketDepth/2.);//Larger rectangle
	G4Box* InnerBracket = new G4Box("TakeAway", fBracketSideLength/2.,inner, fBracketDepth);//Smaller rectangle

	G4ThreeVector sTrans(-fBracketBackWidth, 0., 0.);

	G4SubtractionSolid* Bracket = new G4SubtractionSolid("Bracket", OverBracket, InnerBracket, 0, sTrans);

	G4Tubs* screw = new G4Tubs("screw", 0, fBracketScrewRad, fBracketScrewHeight/2., 0, 360*CLHEP::deg);

	G4Material* material = G4Material::GetMaterial("Aluminum");
	G4LogicalVolume* fSourceBracketLog = new G4LogicalVolume(Bracket, material, "Source_Bracket", 0, 0, 0);
	fSourceBracketLog->SetVisAttributes(vis_att);

	vis_att = new G4VisAttributes(G4Colour(CU_COL));
	vis_att->SetVisibility(true); 
	material = G4Material::GetMaterial("Brass");
	G4LogicalVolume* fScrewLog = new G4LogicalVolume(screw, material, "screw", 0, 0, 0);
	fScrewLog->SetVisAttributes(vis_att);	


	G4double placeZ = fMEezagSourceZ-fMEezagSourceHeight/2.-fBracketDepth/2.;
	G4ThreeVector move = G4ThreeVector(fBracketoffset,0.,placeZ);
	move.rotateZ(fBracketAngle);  
	G4RotationMatrix* rot=new G4RotationMatrix();
	rot->rotateZ(-fBracketAngle);  
	new G4PVPlacement(rot,move,fSourceBracketLog, "Source_Bracket", expHallLog,false, 0); 

	placeZ -= (fBracketScrewHeight+fBracketDepth)/2.;
	G4double placeX = fBracketoffset+(fBracketSideLength)/2.-fBracketScrewRad;
	G4double placeY =fBracketBackLength/2.-fBracketScrewRad;
	G4ThreeVector moveB = G4ThreeVector(placeX,placeY,placeZ);
	G4ThreeVector moveC = G4ThreeVector(placeX,-placeY,placeZ);
	moveB.rotateZ(fBracketAngle);  
	moveC.rotateZ(fBracketAngle); 

	new G4PVPlacement(new G4RotationMatrix,moveB,fScrewLog, "bracket_screw", expHallLog,false, 0); 
	new G4PVPlacement(new G4RotationMatrix,moveC,fScrewLog, "bracket_screw", expHallLog,false, 0); 

	BuildMEezagSource(expHallLog);

}

// **************************************************************************************************
// **************************************************************************************************
// *************************************PLACEMENT****************************************************
// **************************************************************************************************
// **************************************************************************************************


void ApparatusSpiceTargetChamber::PlaceBunnyTargets(G4LogicalVolume* expHallLog){
	PlaceBunnyTarget(2,expHallLog);
	PlaceBunnyTarget(4,expHallLog);
}

void ApparatusSpiceTargetChamber::PlaceBunnyTarget(G4int i, G4LogicalVolume* expHallLog) {


	G4ThreeVector WheelOffest(-fTargetOffsetRadius,0,0);

	G4ThreeVector TargetOffset(0,-fTargetOffsetRadius,fTargetWheelFaceZ-fBunnyRodLength-fBunnyTargetThickness/2.);

	double ang=fTargetWheelRotation+fTargetWheelHoleSpaceAng*i;
	TargetOffset.rotateZ(-ang);

	TargetOffset+=WheelOffest;

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(ang-90*CLHEP::deg);   

	fFramePhysical = new G4PVPlacement(sRotate, TargetOffset, fFrameLogical, "BunnyTargetFrame", expHallLog, false, 0);



	sRotate = new G4RotationMatrix;
	G4ThreeVector RodOffset(0,10*CLHEP::mm,(fBunnyRodLength+fBunnyTargetThickness)/2.);
	RodOffset.rotateZ(-45*CLHEP::deg);
	RodOffset.rotateZ(-ang);

	fRodPhysical = new G4PVPlacement(sRotate, TargetOffset+RodOffset, fRodLogical, "BunnyTargetRods", expHallLog, false, 0);
	RodOffset.rotateZ(90*CLHEP::deg);
	fRodPhysical = new G4PVPlacement(sRotate, TargetOffset+RodOffset, fRodLogical, "BunnyTargetRods", expHallLog, false, 0);


}


void ApparatusSpiceTargetChamber::PlaceEvapTargets(G4LogicalVolume* expHallLog){
	PlaceEvapTarget(2,expHallLog);
	PlaceEvapTarget(4,expHallLog);
}

void ApparatusSpiceTargetChamber::PlaceEvapTarget(G4int i, G4LogicalVolume* expHallLog) {

	G4ThreeVector WheelOffest(-fTargetOffsetRadius,0,0);

	G4ThreeVector TargetOffset(0,-fTargetOffsetRadius,fTargetWheelFaceZ-fEvapTargetThickness/2.);

	double ang=fTargetWheelRotation+fTargetWheelHoleSpaceAng*i;
	TargetOffset.rotateZ(-ang);

	TargetOffset+=WheelOffest;

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(fTargetWheelRotation);   

	fFramePhysical = new G4PVPlacement(sRotate, TargetOffset, fFrameLogical, "EvapTargetFrame", expHallLog, false, 0);

	sRotate = new G4RotationMatrix;
	G4ThreeVector RodOffset(0,-fEvapTargetRadius,-(fEvapTargetThickness+fEvapTargetScrewHeight)/2.);
	RodOffset.rotateZ(-fTargetWheelRotation);

	fRodPhysical = new G4PVPlacement(sRotate, TargetOffset+RodOffset, fRodLogical, "EvapTargetScrew", expHallLog, false, 0);

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
	rotate->rotateZ(fTargetWheelRotation); 

	G4double offset = -fTargetWheelOffset;
	G4double ZPosition =fTargetWheelFaceZ+(fTargetWheelThickness-fTargetWheelIndentDepth)/2.; 

	G4ThreeVector move(offset, 0, ZPosition);
	fTargetWheelPhys = new G4PVPlacement(rotate, move, fTargetWheelLog,
			"TargetWheel", expHallLog,
			false,0);
}// end:PlaceTargetWheel()

void ApparatusSpiceTargetChamber::PlaceCollimator(G4LogicalVolume* expHallLog) {
	G4RotationMatrix* rotate = new G4RotationMatrix;  
	rotate->rotateZ(fTargetWheelRotation); 

	G4double offset = -fTargetWheelOffset;
	G4double ZPosition =fTargetWheelFaceZ-(fTargColPlateThick/2.); 
	G4ThreeVector move(offset, 0, ZPosition);

	fCollimPhysical = new G4PVPlacement(rotate, move, fCollimLogical, "collimator", expHallLog, false, 0);
}

void ApparatusSpiceTargetChamber::PlaceTargetWheelGears(G4LogicalVolume* expHallLog)
{
	G4double ZOffset = fAllGearZOffset;
	G4double FirstGearOffset = fFirstGearPlaneOffset;
	G4double SecondGearOffset = fSecondGearPlaneOffset;
	G4double ThirdGearOffset = fThirdGearPlaneOffset;

	G4ThreeVector move(FirstGearOffset, 0, ZOffset);
	fFirstGearPhys = new G4PVPlacement(0, move, fFirstGearLog,
			"FirstGear", expHallLog,
			false,0);

	move.setX(SecondGearOffset);
	fSecondGearPhys = new G4PVPlacement(0, move, fSecondGearLog,
			"SecondGear", expHallLog,
			false,0);

	move.setX(-ThirdGearOffset);
	fThirdGearPhys = new G4PVPlacement(0, move, fThirdGearLog,
			"ThirdGear", expHallLog,
			false,0);
} // end PlaceTargetWheelGears()

void ApparatusSpiceTargetChamber::PlaceTargetWheelGearPlates(G4LogicalVolume* expHallLog)
{
	G4double ZOffset = fAllGearZOffset + fAllGearThickness/2. + fGearPlateThickness/2.;

	G4double FirstGearOffset = fFirstGearPlaneOffset;
	G4double SecondGearOffset = fSecondGearPlaneOffset;

	G4ThreeVector move(FirstGearOffset, 0, ZOffset);
	fGearPlateOnePhys = new G4PVPlacement(0, move, fGearPlateOneLog,
			"GearPlateOne", expHallLog,
			false,0);
	move.setX(SecondGearOffset);
	fGearPlateOnePhys = new G4PVPlacement(0, move, fGearPlateTwoLog,
			"GearPlateTwo", expHallLog,
			false,0);	
}	// end::PlaceTargetWheelGearPlates()

void ApparatusSpiceTargetChamber::PlaceGearStick(G4LogicalVolume* expHallLog)
{
	G4double ZOffset =fAllGearZOffset -fGearStickLength/2. -fAllGearThickness/2.;
	G4double PlaneOffset = fFirstGearPlaneOffset;

	G4ThreeVector move(PlaneOffset, 0 , ZOffset);
	fGearStickPhys = new G4PVPlacement(0, move, fGearStickLog,
			"GearStick", expHallLog,
			false,0);
} // end PlaceGearStick();

void ApparatusSpiceTargetChamber::PlaceTargetMountPlate(G4LogicalVolume* expHallLog) 
{
	G4RotationMatrix* rotate = new G4RotationMatrix;  
	G4ThreeVector move(0, 0, fTargetMountPlateZOffset);
	fTargetMountPlatePhys = new G4PVPlacement(rotate, move, fTargetMountPlateLog,
			"TargetMountPlate", expHallLog,
			false,0);

	for(int i=0;i<7;i++){
		G4double ang=(22.5+i*45)*CLHEP::deg;
		G4ThreeVector moveR(fMountPlateRodPlaceR, 0, fMountPlateRodZoffset);
		moveR.rotateZ(ang);
		fTargetMountRodPhys = new G4PVPlacement(rotate, moveR, fTargetMountRodLog,"TargetMountRod", expHallLog,false,0);
	}

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

	G4double Ang=(45+copyID*90)*deg;
	G4ThreeVector move(fPsBoltPlaceRadius, 0, fPhotonShieldBackFacePos+fPsBoltLength/2.-fPsBoltBackShield);
	move.rotateZ(-Ang);

	fTargetBoltPhys = new G4PVPlacement(0, move, fShieldBoltLog,
			"ShieldBolt", expHallLog, 
			false,0);

} // end::PlacePhotonShieldClampBoots()

void ApparatusSpiceTargetChamber::PlaceMagnetAndCover(G4int copyID, G4LogicalVolume* expHallLog)
{
	// ** Position Co-ordinates
	G4double magnetPosX = fPlateOneEdgeX + fPlateOneLength/2.;
	G4double magnetPosZ = -(fMagDistanceFromTarget + fPlateOneHeight/2.) + fMiddleRingOffset; // - (4mm + 50/2mm)  + 0 
	G4double Ang=(45+copyID*90)*deg;
	G4RotationMatrix* rotation=new G4RotationMatrix();
	rotation->rotateZ(Ang);
	G4ThreeVector move(magnetPosX, 0, magnetPosZ);
	move.rotateZ(-Ang);
	fMagnetPhys = new G4PVPlacement(rotation, move, fMagnetLog,"CollectorMagnet",expHallLog,false,0);

	G4double MagCoverOffsetX = -(fPsClampThicknessInner-fPsClampThicknessOuter)/2;
	G4double MagCoverOffsetZ = -(fPsClampThicknessBottom-fPsClampThicknessTop)/2; 
	G4ThreeVector movecov(magnetPosX+MagCoverOffsetX, 0, magnetPosZ+MagCoverOffsetZ);
	movecov.rotateZ(-Ang);
	fMagnetCoverPhys = new G4PVPlacement(rotation, movecov, fMagnetCoverLog,"MagnetCovering",expHallLog,false,0);

	G4double ZPosition = -(fMagDistanceFromTarget + fMclampChamberHeight/2.) + fMiddleRingOffset;
	G4double XPosition = fTargetChamberCylinderFlatRadius - fMclampChamberLength/2.;
	G4ThreeVector movecmp(XPosition, 0, ZPosition);
	movecmp.rotateZ(-Ang);
	fMclampChamberPhys = new G4PVPlacement(rotation, movecmp, fMclampChamberLog,"MagnetClampChamber", expHallLog,false,0);
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
	G4ThreeVector sTranslate(fTargetChamberCylinderFlatRadius-6.3*CLHEP::mm,0,-fTargetChamberCylinderLength/2.-fTargetChamberDistanceFromTarget);
	sTranslate.rotateZ(-22.5*CLHEP::deg);

	G4RotationMatrix* sRotate = new G4RotationMatrix;
	sRotate->rotateZ(22.5*CLHEP::deg);

	fS3CasePhysical = new G4PVPlacement(sRotate, sTranslate, fS3CaseLogical, "s3CableCase", expHallLog, false, 0);
}
void ApparatusSpiceTargetChamber::PlaceConicalCollimator(G4LogicalVolume* expHallLog) {

	G4double ZPosition = fPhotonShieldBackFacePos + fMiddleRingOffset+ fPhotonShieldLength + fPsTargetClampThickness;
	ZPosition+=fConicalTabHeight-fConicalLength1/2.;
	fConicalCollimatorPhys = new G4PVPlacement(0,G4ThreeVector(0, 0, ZPosition),fConicalCollimatorLog,"ConicalCollimator",expHallLog,false,0);

}

void ApparatusSpiceTargetChamber::PlaceXrayInsert(G4LogicalVolume* expHallLog){

	G4double CollimatorEnd = fPhotonShieldBackFacePos + fMiddleRingOffset+ fPhotonShieldLength + fPsTargetClampThickness;
	CollimatorEnd+=fConicalTabHeight-fConicalLength1-fConicalLength2-fConicalLength3;
	CollimatorEnd-=fInsertCylZLength/2.;
	CollimatorEnd+=fInsertCylZLength-fInsertConeZLength;

	G4double ShieldBack = fPhotonShieldBackFacePos + fMiddleRingOffset- fPsDetClampThickness;
	ShieldBack-=(fInsertEdgeZLength+fInsertCylZLength/2);

	G4double ZPosition = CollimatorEnd;
	if(ShieldBack<ZPosition)ZPosition=ShieldBack;
	G4ThreeVector move(0, 0, ZPosition);

	fXRayInsertPhys = new G4PVPlacement(new G4RotationMatrix(),                       //rotation
			move,                    //at position
			fXRayInsertLog,             //its logical volume
			"XrayClamp",                //its name
			expHallLog,                //its mother  volume
			false,                   //no boolean operation
			false,                       //copy number
			0);          //overlaps checking


}

