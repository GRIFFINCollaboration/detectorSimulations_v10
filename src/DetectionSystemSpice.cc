#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemSpice.hh"
#include "HistoManager.hh" // for spice histogram array when placing detector

DetectionSystemSpice::DetectionSystemSpice() :
	// LogicalVolumes
	fSiInnerGuardRingLog(0),
	fSiOuterGuardRingLog(0)
{
	/////////////////////////////////////////////////////////////////////
	// SPICE Physical Properties
	/////////////////////////////////////////////////////////////////////
	//-----------------------------//
	// Materials     			   //
	//-----------------------------//
	fWaferMaterial          = "Silicon";
	fDetectorMountMaterial = "Aluminum"; // CHECK ?
	fAnnularClampMaterial  = "Peek"; // CHECK ?

	// ----------------------------
	// Dimensions of Detector Mount
	// ----------------------------
	fDetectorMountLength = 148*mm;
	fDetectorMountWidth = 130*mm;
	fDetectorMountThickness = 8*mm;
	fDetectorMountInnerRadius = 52.8*mm;
	fDetectorMountLipRadius = 48.3*mm;
	fDetectorMountLipThickness = 1*mm;
	fDetectorMountAngularOffset = 0*deg;
	fDetectorToTargetDistance = 115*mm;
	fDetectorThickness = 6*mm;

	// ---------------------------
	// Dimensions of Annular Clamp
	// ---------------------------
	fAnnularClampThickness = 2*mm;
	fAnnularClampLength = 19.2*mm;
	fAnnularClampWidth = 20*mm;
	fAnnularClampPlaneOffset = 47.8*mm; //from beam line to first edge  

	//-----------------------------//
	// parameters for the annular  //
	// planar detector crystal     //
	//-----------------------------//
	fSiDetCrystalOuterDiameter = 94.*mm;
	fSiDetCrystalInnerDiameter = 16.*mm;
	fSiDetCrystalThickness = 6.15*mm;
	fSiDetRadialSegments = 10.;
	fSiDetPhiSegments = 12.;

	//-----------------------------//
	// parameters for guard ring   //
	//-----------------------------//
	fSiDetGuardRingOuterDiameter = 102*mm;
	fSiDetGuardRingInnerDiameter = 10*mm;
}

DetectionSystemSpice::~DetectionSystemSpice()
{   
	// LogicalVolumes in ConstructSPICEDetectionSystem
	delete fDetectorMountLog;
	delete fAnnularClampLog;
	for(int i = 0; i < 10; ++i) {
		delete fSiDetSpiceRingLog[i];
	}
	delete fSiInnerGuardRingLog;
	delete fSiOuterGuardRingLog;

	//Physical volumes 
	delete fDetectorMountPhys;
	delete fAnnularClampPhys;
}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemSpice::Build()
{
	fAssembly = new G4AssemblyVolume();

	// Loop through each ring ...
	for(int ringID=0; ringID<10; ringID++) {
		// Build fAssembly volumes
		fAssemblySiRing[ringID] = new G4AssemblyVolume();
		// Build the logical segment of the Silicon Ring
		BuildSiliconWafer(ringID);  
	} // end for(int ringID)

	BuildInnerGuardRing();
	BuildOuterGuardRing();
	BuildDetectorMount();
	BuildAnnularClamps();

	return 1;
} // end Build

//---------------------------------------------------------//
// "place" function called in DetectorMessenger            //
// if detector is added                                    //
//---------------------------------------------------------//
G4int DetectionSystemSpice::PlaceDetector(G4LogicalVolume* ExpHallLog, G4int nRings)
{
	PlaceDetectorMount(ExpHallLog);
	PlaceAnnularClamps(ExpHallLog);   
	PlaceGuardRing(ExpHallLog);//clears the code in det construction
	//To see how the segments fill put 1 detector ring in (via macro) and 1,2,3 segments (changed here) to see how they add.
	//segments add following JanalysisToolkit convention
	G4int NumberSeg = 12; // Segments in Phi for loop below - originally 12
	G4int segmentID=0;
	G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
	G4ThreeVector pos(0,0,-annularDetectorDistance); 
	for(int ring=0; ring<nRings; ring++) {
		for(int Seg=0; Seg<NumberSeg; Seg++) {
			G4int TotalNumberSeg = (G4int)fSiDetPhiSegments; // total number of segments = 12
			G4double angle = (360.*deg/TotalNumberSeg)*(-Seg); // Seg = {0, ...,11} //changed to -seg and fills in correctly
			G4RotationMatrix* rotate = new G4RotationMatrix;
			rotate->rotateZ(-180*deg-angle); // the axis are inverted, this operation will correct for it  [MHD : 03 April 2014]
			//affecting value here (originally 210) may help mapping
			//180 deg gave correct start position for 1st ring, needs to now fill opposite direction
			//could bring in loops to here as to iterate and allow for staggered start position (add 30 deg for each ring etc)
			//as filling in reverse, can establish the histos in reverse then fill in normally 
			fAssemblySiRing[ring]->MakeImprint(ExpHallLog, pos, rotate, Seg);
			segmentID++;
		} // end for(int Seg=0; Seg<NumberSeg; Seg++)
	} // end for(int ring = 0; ring<nRings; ring++)

	return 1;
}

G4int DetectionSystemSpice::PlaceGuardRing(G4LogicalVolume* ExpHallLog)
{
	G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
	G4ThreeVector pos(0,0,-annularDetectorDistance); 
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);
	fAssembly->MakeImprint(ExpHallLog, pos, rotate, 0);

	return 1;
}

void DetectionSystemSpice::PlaceDetectorMount(G4LogicalVolume* ExpHallLog)
{
	G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
	G4ThreeVector pos(0,0,-annularDetectorDistance); 
	G4double DetectorMountGap = fDetectorMountThickness 
		- fDetectorMountLipThickness - fDetectorThickness;

	G4double ZOffset = - fDetectorMountThickness/2. + DetectorMountGap ;
	G4ThreeVector offset(0, 0, ZOffset);

	pos = pos + offset ; 
	G4RotationMatrix* rotate = new G4RotationMatrix(fDetectorMountAngularOffset, 0, 0);
	fDetectorMountPhys = new G4PVPlacement(rotate, pos, fDetectorMountLog,
			"DetectorMount", ExpHallLog, false, 0);
} // end::PlaceDetectorMount()

void DetectionSystemSpice::PlaceAnnularClamps(G4LogicalVolume* ExpHallLog) {
	G4double annularDetectorDistance = 115*mm /*+ 150*mm*/;
	G4ThreeVector pos(0,0,-annularDetectorDistance); 
	G4double ZOffset = fAnnularClampThickness/2. ;
	G4double XOffset = (fAnnularClampPlaneOffset
			+ fAnnularClampLength/2.) 
		* cos(fDetectorMountAngularOffset + 45*deg);
	G4double YOffset = (fAnnularClampPlaneOffset
			+ fAnnularClampLength/2.)
		* sin(fDetectorMountAngularOffset + 45*deg);
	G4ThreeVector offset(-XOffset, -YOffset, ZOffset);

	pos = pos + offset ; 	 
	G4RotationMatrix* rotate = new G4RotationMatrix(fDetectorMountAngularOffset + 45*deg, 0, 0);
	fAnnularClampPhys = new G4PVPlacement(rotate, pos, fAnnularClampLog,"fAnnularClamp", ExpHallLog,
			false,0);
} // end::PlacefAnnularClamps()

//---------------------------------------------------------//
// build functions for different parts                     //
// called in main build function                           //
//---------------------------------------------------------//
G4int DetectionSystemSpice::BuildSiliconWafer(G4int RingID)  // RingID = { 0, 9 }
{
	// Define the material, return error if not found
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if( !material ) {
		G4cout << " ----> Material " << fWaferMaterial 
			<< " not found, cannot build the detector shell! " << G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	VisAtt->SetVisibility(true);

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double ZPosition		= -(fSiDetCrystalThickness/2.);
	G4ThreeVector move 		= ZPosition * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);

	// construct solid
	G4Tubs* fSiDetSpiceRingSec = BuildCrystal(RingID);
	// construct logical volume
	if(fSiDetSpiceRingLog[RingID] == NULL) {
		G4String ringName = "fSiDetSpiceRing_";
		G4String ringID = G4UIcommand::ConvertToString(RingID);
		ringName += ringID;
		ringName += "_Log";

		fSiDetSpiceRingLog[RingID] = new G4LogicalVolume(fSiDetSpiceRingSec, material, ringName, 0, 0, 0);
		fSiDetSpiceRingLog[RingID]->SetVisAttributes(VisAtt);
	}

	fAssemblySiRing[RingID]->AddPlacedVolume(fSiDetSpiceRingLog[RingID], move, rotate);

	return 1;
}

G4int DetectionSystemSpice::BuildInnerGuardRing()
{
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if(material == NULL) {
		G4cout << " ----> Material " << fWaferMaterial << " not found, cannot build the inner guard ring of Spice! " << G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	VisAtt->SetVisibility(true);

	G4Tubs* innerGuardRing = new G4Tubs("innerGuardRing",
			fSiDetGuardRingInnerDiameter/2.,
			fSiDetCrystalInnerDiameter/2.,
			fSiDetCrystalThickness/2.,0,360);

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double ZPosition		= -(fSiDetCrystalThickness/2.);
	G4ThreeVector move 		= ZPosition * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);

	//logical volume
	if(fSiInnerGuardRingLog == NULL) {
		fSiInnerGuardRingLog = new G4LogicalVolume(innerGuardRing, material, "innerGuardRing", 0,0,0);
		fSiInnerGuardRingLog->SetVisAttributes(VisAtt);
	}

	fAssembly->AddPlacedVolume(fSiInnerGuardRingLog, move, rotate);

	return 1;
}

G4int DetectionSystemSpice::BuildOuterGuardRing()
{
	G4Material* material = G4Material::GetMaterial(fWaferMaterial);
	if(material == NULL) {
		G4cout << " ----> Material " << fWaferMaterial << " not found, cannot build the outer guard ring of Spice! " << G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	VisAtt->SetVisibility(true);

	G4Tubs* outerGuardRing = new G4Tubs("outerGuardRing",
			fSiDetCrystalOuterDiameter/2.,
			fSiDetGuardRingOuterDiameter/2.,
			fSiDetCrystalThickness/2.,0,360);

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double ZPosition		= -(fSiDetCrystalThickness/2.);
	G4ThreeVector move 		= ZPosition * direction;
	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateZ(0*deg);

	//logical volume
	if(fSiOuterGuardRingLog == NULL)	{
		fSiOuterGuardRingLog = new G4LogicalVolume(outerGuardRing, material, "outerGuardRing", 0,0,0);
		fSiOuterGuardRingLog->SetVisAttributes(VisAtt);
	}

	fAssembly->AddPlacedVolume(fSiOuterGuardRingLog, move, rotate);

	return 1;
}


void DetectionSystemSpice::BuildDetectorMount() {

	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(AL_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	// Box
	G4double BoxHalfWidth = fDetectorMountWidth/2.;
	G4double BoxHalfLength = fDetectorMountLength/2.;
	G4double BoxHalfThickness = fDetectorMountThickness/2.;
	// Inner Radius
	G4double LipRadius = fDetectorMountLipRadius;
	//G4double LipHalfThickness = fDetectorMountLipThickness/2.;
	G4double BoxCutRadius = fDetectorMountInnerRadius;
	// Annular Clamp
	G4double ClampHalfThickness = fAnnularClampThickness;
	G4double ClampHalfWidth = fAnnularClampWidth/2.;
	G4double ClampHalfLength = fAnnularClampLength/2.;

	// ** Shapes
	G4Box* MountBox = new G4Box("MountBox", BoxHalfWidth, BoxHalfLength, BoxHalfThickness);
	G4Tubs* InnerRadiusCut = new G4Tubs("InnerRadiusCut", 0, LipRadius, 2*BoxHalfThickness, 0, 360*deg);
	G4Tubs* LipCut = new G4Tubs("LipCut", 0, BoxCutRadius, BoxHalfThickness, 0, 360*deg);
	G4Box* fAnnularClamp = new G4Box("fAnnularClamp", ClampHalfWidth, ClampHalfLength, ClampHalfThickness);

	G4SubtractionSolid* DetectorMountPre = new G4SubtractionSolid("DetectorMountPre", MountBox, InnerRadiusCut);
	G4ThreeVector trans(0, 0, fDetectorMountLipThickness);
	G4SubtractionSolid* DetectorMount = new G4SubtractionSolid("DetectorMount", DetectorMountPre, LipCut, 0, trans);

	G4double PlaneOffset = (fAnnularClampPlaneOffset + ClampHalfLength) / sqrt(2.);
	G4double ZOffset = BoxHalfThickness;
	G4ThreeVector move(PlaneOffset, PlaneOffset, ZOffset);
	G4RotationMatrix* rotate = new G4RotationMatrix(45*deg, 0, 0);
	G4SubtractionSolid* DetectorMount2 = new G4SubtractionSolid("DetectorMount2", DetectorMount, fAnnularClamp, rotate, move);
	move.setX(-PlaneOffset);
	rotate->rotateZ(90*deg);
	G4SubtractionSolid* DetectorMount3 = new G4SubtractionSolid("DetectorMount3", DetectorMount2, fAnnularClamp, rotate, move);
	move.setY(-PlaneOffset);
	rotate->rotateZ(90*deg);
	G4SubtractionSolid* DetectorMount4 = new G4SubtractionSolid("DetectorMount4", DetectorMount3, fAnnularClamp, rotate, move);
	move.setX(PlaneOffset);
	rotate->rotateZ(90*deg);
	G4SubtractionSolid* DetectorMount5 = new G4SubtractionSolid("DetectorMount5", DetectorMount4, fAnnularClamp, rotate, move);

	// ** Logical
	G4Material* detectorMountMaterial = G4Material::GetMaterial(fDetectorMountMaterial);
	fDetectorMountLog = new G4LogicalVolume(DetectorMount5, detectorMountMaterial, "DetectorMountLog", 0, 0, 0);
	fDetectorMountLog->SetVisAttributes(VisAtt);
} // end::BuildDetectorMount()

void DetectionSystemSpice::BuildAnnularClamps() {
	// ** Visualisation
	G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(PEEK_COL));
	VisAtt->SetVisibility(true);

	// ** Dimensions
	G4double ClampHalfLength = fAnnularClampLength/2.;
	G4double ClampHalfWidth = fAnnularClampWidth/2.;
	G4double ClampHalfThickness = fAnnularClampThickness/2.;
	// Distance
	G4double BeamClampDistance = fAnnularClampPlaneOffset + ClampHalfLength;

	// ** Shapes
	G4Box* fAnnularClamp = new G4Box("fAnnularClamp", ClampHalfWidth, ClampHalfLength, ClampHalfThickness);

	G4ThreeVector move(2*BeamClampDistance, 0, 0);
	G4UnionSolid* DoubleClamps = new G4UnionSolid("DoubleClamps", fAnnularClamp, fAnnularClamp, 0, move);

	G4Box* fAnnularClamp2 = new G4Box("fAnnularClamp2", ClampHalfLength, ClampHalfWidth, ClampHalfThickness);
	G4ThreeVector trans(0, 2*BeamClampDistance, 0);
	G4UnionSolid* DoubleClamps2 = new G4UnionSolid("DoubleClamps2", fAnnularClamp2, fAnnularClamp2, 0, trans);

	G4ThreeVector trans2(BeamClampDistance, -BeamClampDistance, 0);
	G4UnionSolid* FourClamps = new G4UnionSolid("FourClamps", DoubleClamps, DoubleClamps2, 0, trans2);

	// ** Logical
	G4Material* annularClampMaterial = G4Material::GetMaterial(fAnnularClampMaterial);
	fAnnularClampLog = new G4LogicalVolume(FourClamps, annularClampMaterial, "fAnnularClampLog", 0, 0, 0);
	fAnnularClampLog->SetVisAttributes(VisAtt);
} // end::BuildfAnnularClamps()  

/////////////////////////////////////////////////////////
// Build one segment of Spice, 			       //
// the geometry depends on the distance from the center//
/////////////////////////////////////////////////////////
G4Tubs* DetectionSystemSpice::BuildCrystal(G4int RingID)
{
	// define angle, length, thickness, and inner and outer diameter
	// of silicon detector segment
	G4double TubeElementLength = (fSiDetCrystalOuterDiameter - fSiDetCrystalInnerDiameter)/(2*(fSiDetRadialSegments));
	G4double TubeElementAngularWidth = (360./fSiDetPhiSegments)*deg;
	G4double TubeElementInnerRadius = (fSiDetCrystalInnerDiameter)/2.0 + TubeElementLength*(RingID);
	G4double TubeElementOuterRadius = ((G4double)fSiDetCrystalInnerDiameter)/2.0 + TubeElementLength*(RingID+1);
	G4double TubeElementHalfThickness = (fSiDetCrystalThickness)/2.0;

	// establish solid
	G4Tubs* CrystalBlock = new G4Tubs("CrystalBlock",TubeElementInnerRadius,TubeElementOuterRadius,TubeElementHalfThickness,0,TubeElementAngularWidth);

	return CrystalBlock;
}//end ::BuildCrystal
