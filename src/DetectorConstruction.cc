
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
/// \file analysis/shared/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.cc 77256 2013-11-22 10:10:23Z gcosmo $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include <iostream>
#include <sstream>

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4UserLimits.hh"

#include "G4Threading.hh"

//#include "DetectionSystemGammaTracking.hh"
#include "DetectionSystem8pi.hh"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "DetectionSystemDescant.hh"
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "ApparatusDescantStructure.hh"

#include "DetectionSystemTestcan.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "NonUniformMagneticField.hh"

#include "DetectionSystemGriffin.hh"
#include "DetectionSystemSceptar.hh"
#include "DetectionSystemSpice.hh"
#include "DetectionSystemTrific.hh"
#include "DetectionSystemPaces.hh"
#include "DetectionSystemSodiumIodide.hh"
#include "DetectionSystemLanthanumBromide.hh"
#include "ApparatusGenericTarget.hh"
#include "ApparatusLayeredTarget.hh"
#include "ApparatusSpiceTargetChamber.hh"
#include "Apparatus8piVacuumChamber.hh"
#include "Apparatus8piVacuumChamberAuxMatShell.hh"
#include "ApparatusGriffinStructure.hh"

#include "DetectionSystemBox.hh" // New file
#include "DetectionSystemGrid.hh"

#include "DetectionSystemAncillaryBGO.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool operator==(const DetectorProperties& lhs, const DetectorProperties& rhs)
{
	return (lhs.systemID == rhs.systemID && lhs.detectorNumber == rhs.detectorNumber && lhs.crystalNumber == rhs.crystalNumber);
}

DetectorConstruction::DetectorConstruction() :
	fSolidWorld(nullptr),
	fLogicWorld(nullptr),
	fPhysiWorld(nullptr)
{
	fWorldSizeX  = fWorldSizeY = fWorldSizeZ = 10.0*m;

	fBoxMat = "G4_WATER";
	fBoxThickness = 0.0*mm;
	fBoxInnerDimensions = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	fBoxColour = G4ThreeVector(0.0,0.0,1.0);

	fGridMat = "G4_WATER";
	fGridSize = 0.0*mm;
	fGridDimensions = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	fGridColour = G4ThreeVector(1.0,0.0,0.0);

	// materials
	DefineMaterials();

	//  builtDetectors = false;


	fMatWorldName = "G4_AIR";

	// Generic Target Apparatus
	fSetGenericTargetMaterial   = false;
	fSetGenericTargetDimensions = false;
	fSetGenericTargetPosition   = false;

	// Field Box
	fSetFieldBoxMaterial= false;//think this has been removed 17/8
	fSetFieldBoxDimensions= false;
	fSetFieldBoxPosition= false;
	fSetFieldBoxMagneticField= false;

	// parameters to suppress:

	DefineSuppressedParameters();

	// Shield Selection Default

	fUseTigressPositions = false;

	fDetectorShieldSelect = 1 ; // Include suppressors by default.
	fExtensionSuppressorLocation = 0 ; // Back by default (Detector Forward)
	fHevimetSelector = 0 ; // Chooses whether or not to include a hevimet

	fCustomDetectorNumber 		= 1 ; // detNum
	fCustomDetectorPosition  = 1 ; // posNum

	// create commands for interactive definition

	fDetectorMessenger = new DetectorMessenger(this);

	// ensure the global field is initialized
	//(void)GlobalField::getObject();

	//expHallMagField = new MagneticField(); // Global field is set to zero

	fGriffinDetectorsMapIndex = 0;
	for(G4int i = 0; i < 16; ++i) {
		fGriffinDetectorsMap[i] = 0;
		for(G4int j = 0; j < 4; ++j) {
			fGriffinDeadLayer[i][j] = -1;
		}
	}

	fDescantColor = "white";
	fDescantRotation.setX(M_PI);
	fDescantRotation.setY(0.);
	fDescantRotation.setZ(0.);

	fApparatusLayeredTarget=0;

	fGridCell = false;
	fGriffin  = false;
	fLaBr     = false;
	fAncBgo   = false;
	fNaI      = false;
	fSceptar  = false;
	fEightPi  = false;
	fSpice    = false;
	fPaces    = false;
	fDescant  = false;
	fTestcan  = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {	
	delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {

	// Replaced by ConstructDetectionSystems

	// Experimental hall (world volume)
	// search the world material by its name

	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();

	G4Material* matWorld = G4Material::GetMaterial(fMatWorldName);

	if(!matWorld) {
		G4cout<<" ----> Material "<<fMatWorldName<<" not found, cannot build world! "<<G4endl;
		return 0;
	}

	fSolidWorld = new G4Box("World", fWorldSizeX/2,fWorldSizeY/2,fWorldSizeZ/2);

	fLogicWorld = new G4LogicalVolume(fSolidWorld,		//its solid
			matWorld,	//its material
			"World");		//its name

	fPhysiWorld = new G4PVPlacement(0,                  //no rotation
			G4ThreeVector(),	//at (0,0,0)
			fLogicWorld,         //its logical volume
			"World",            //its name
			0,                  //its mother  volume
			false,              //no boolean operation
			0);                 //copy number

	// Visualization Attributes

	fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible); // The following block of code works too.

	//  G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	//  worldVisAtt->SetForceWireframe(true);
	//  worldVisAtt->SetVisibility(worldVis);
	//  fLogicWorld->SetVisAttributes(worldVisAtt);
	//  fLogicWorld = fLogicWorld;

	return fPhysiWorld;
}

void DetectorConstruction::SetWorldMaterial(G4String name) {
	fMatWorldName = name;
	UpdateGeometry(); // auto update
}

void DetectorConstruction::SetWorldDimensions(G4ThreeVector vec) {
	fWorldSizeX = vec.x() ;
	fWorldSizeY = vec.y() ;
	fWorldSizeZ = vec.z() ;
	UpdateGeometry(); // auto update
}

void DetectorConstruction::SetWorldVis(G4bool vis) {
	fWorldVis = vis;
	UpdateGeometry(); // auto update
}

void DetectorConstruction::SetWorldStepLimit(G4double step) {
	if(fLogicWorld == nullptr) {
		Construct();
	}
	fLogicWorld->SetUserLimits(new G4UserLimits(step));
}

//void DetectorConstruction::SetWorldMagneticField(G4ThreeVector vec)
//{
//    //expHallMagField->SetFieldValue(G4ThreeVector(vec.x(),vec.y(),vec.z()));
//}
void DetectorConstruction::LayeredTargetAdd(G4String Material, G4double Areal)
{
	if(!fApparatusLayeredTarget)fApparatusLayeredTarget=new ApparatusLayeredTarget(fDetEffPosition.z());
	if(fApparatusLayeredTarget->BuildTargetLayer(Material,Areal)){
		if(fLogicWorld == nullptr) {
			Construct();
		}  
		fApparatusLayeredTarget->PlaceTarget(fLogicWorld);
	}
}


void DetectorConstruction::SetTabMagneticField(G4String PathAndTableName, G4double z_offset, G4double z_rotation)
{
	//const char * c = PathAndTableName.c_str();///conversion attempt using .c_str() (a string member function)
	new NonUniformMagneticField(PathAndTableName,z_offset,z_rotation);///addition of field for SPICE
}

void DetectorConstruction::UpdateGeometry() {
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

void DetectorConstruction::SetGenericTargetMaterial(G4String name)
{
	DetectorConstruction::fSetGenericTargetMaterial = true;
	DetectorConstruction::fGenericTargetMaterial = name;
	DetectorConstruction::SetGenericTarget();
}

void DetectorConstruction::SetGenericTargetDimensions(G4ThreeVector vec)
{
	DetectorConstruction::fSetGenericTargetDimensions = true;
	DetectorConstruction::fGenericTargetDimensions = vec;
	DetectorConstruction::SetGenericTarget();
}

void DetectorConstruction::SetGenericTargetPosition(G4ThreeVector vec)
{
	DetectorConstruction::fSetGenericTargetPosition = true;
	DetectorConstruction::fGenericTargetPosition = vec;
	DetectorConstruction::SetGenericTarget();
}

void DetectorConstruction::SetGenericTarget()
{
	if(fSetGenericTargetMaterial) {
		if(fSetGenericTargetDimensions) {
			if(fSetGenericTargetPosition) {
				G4String name = fGenericTargetMaterial;
				G4double vecX = fGenericTargetDimensions.x()/mm;
				G4double vecY = fGenericTargetDimensions.y()/mm;
				G4double vecZ = fGenericTargetDimensions.z()/mm;

				ApparatusGenericTarget* pApparatusGenericTarget = new ApparatusGenericTarget();
				pApparatusGenericTarget->Build(name, vecX, vecY, vecZ);
				G4RotationMatrix* rotate = new G4RotationMatrix;
				pApparatusGenericTarget->PlaceApparatus(fLogicWorld, fGenericTargetPosition, rotate);
			}
		}
	}
}


G4double DetectorConstruction::LayeredTargetLayerStart(int layer){
	if(fApparatusLayeredTarget){
		return fApparatusLayeredTarget->LayerStart(layer);
	}
	return fDetEffPosition.z();
}


void DetectorConstruction::AddGrid() {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	if(fGridSize != 0.0*mm)
	{
		DetectionSystemGrid* pGrid = new DetectionSystemGrid(	fGridDimensions.x(),
				fGridDimensions.y(),
				fGridDimensions.z(),
				fGridSize,
				fGridMat,
				fGridColour,
				fGridOffset) ;
		pGrid->Build();
		pGrid->PlaceDetector(fLogicWorld);
	}
}

void DetectorConstruction::AddApparatusSpiceTargetChamber(G4String Options)//parameter sets lens for SPICE - should be matched with field
{
	//Create Target Chamber
	ApparatusSpiceTargetChamber* fApparatusSpiceTargetChamber = new ApparatusSpiceTargetChamber(Options);
	fApparatusSpiceTargetChamber->Build(fLogicWorld);
}

void DetectorConstruction::AddApparatus8piVacuumChamber() {
	//Create Vacuum Chamber
	if(fLogicWorld == nullptr) {
		Construct();
	}

	Apparatus8piVacuumChamber* pApparatus8piVacuumChamber = new Apparatus8piVacuumChamber();
	pApparatus8piVacuumChamber->Build(fLogicWorld);
}

void DetectorConstruction::AddApparatus8piVacuumChamberAuxMatShell(G4double thickness) {
	//Create Shell Around Vacuum Chamber
	if(fLogicWorld == nullptr) {
		Construct();
	}

	Apparatus8piVacuumChamberAuxMatShell* pApparatus8piVacuumChamberAuxMatShell = new Apparatus8piVacuumChamberAuxMatShell();
	pApparatus8piVacuumChamberAuxMatShell->Build(fLogicWorld, thickness);
}

void DetectorConstruction::AddApparatusGriffinStructure(G4int selector) {
	//Create Shell Around Vacuum Chamber
	if(fLogicWorld == nullptr) {
		Construct();
	}

	ApparatusGriffinStructure* pApparatusGriffinStructure = new ApparatusGriffinStructure();
	pApparatusGriffinStructure->Build();

	pApparatusGriffinStructure->Place(fLogicWorld, selector);
}


void DetectorConstruction::AddDetectionSystemSodiumIodide(G4int ndet) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	// Describe Placement
	G4double detectorAngles[8][2] = {{}};
	G4double theta,phi,position;
	G4ThreeVector move,direction;

	detectorAngles[0][0] 	= 0.0;
	detectorAngles[1][0] 	= 45.0;
	detectorAngles[2][0] 	= 90.0;
	detectorAngles[3][0] 	= 135.0;
	detectorAngles[4][0] 	= 180.0;
	detectorAngles[5][0] 	= 225.0;
	detectorAngles[6][0] 	= 270.0;
	detectorAngles[7][0] 	= 315.0;
	detectorAngles[0][1] 	= 90.0;
	detectorAngles[1][1] 	= 90.0;
	detectorAngles[2][1] 	= 90.0;
	detectorAngles[3][1] 	= 90.0;
	detectorAngles[4][1] 	= 90.0;
	detectorAngles[5][1] 	= 90.0;
	detectorAngles[6][1] 	= 90.0;
	detectorAngles[7][1] 	= 90.0;

	DetectionSystemSodiumIodide* pSodiumIodide = new DetectionSystemSodiumIodide() ;
	pSodiumIodide->Build() ;

	for(G4int detectorNumber = 0; detectorNumber < ndet; detectorNumber++)
	{
		phi = detectorAngles[detectorNumber][0]*deg; // Creates a ring in phi plane
		theta = detectorAngles[detectorNumber][1]*deg;

		direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
		position = 25.0*cm + (pSodiumIodide->GetDetectorLengthOfUnitsCM()/2.0);
		move = position * direction;

		G4RotationMatrix* rotate = new G4RotationMatrix; 		//rotation matrix corresponding to direction vector
		rotate->rotateX(theta);
		rotate->rotateY(0);
		rotate->rotateZ(phi+0.5*M_PI);

		pSodiumIodide->PlaceDetector(fLogicWorld, move, rotate, detectorNumber) ;
	}
	fNaI = true;
}

void DetectorConstruction::AddDetectionSystemLanthanumBromide(G4ThreeVector input) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	G4int ndet = G4int(input.x());
	G4double radialpos = input.y()*cm;

	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	pDetectionSystemLanthanumBromide->Build();

	for(G4int detectorNumber = 0; detectorNumber < ndet; detectorNumber++)
	{
		pDetectionSystemLanthanumBromide->PlaceDetector(fLogicWorld, detectorNumber, radialpos);
	}
	fLaBr = true;
}

void DetectorConstruction::AddDetectionSystemAncillaryBGO(G4ThreeVector input) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	G4int ndet = G4int(input.x());
	G4double radialpos = input.y()*cm;
	G4int hevimetopt = G4int(input.z());

	DetectionSystemAncillaryBGO* pDetectionSystemAncillaryBGO = new DetectionSystemAncillaryBGO();
	pDetectionSystemAncillaryBGO->Build();

	for(G4int detectorNumber = 0; detectorNumber < ndet; detectorNumber++) {
		pDetectionSystemAncillaryBGO->PlaceDetector(fLogicWorld, detectorNumber, radialpos, hevimetopt);
	}
	fAncBgo = true;
}



G4double DetectorConstruction::GetLanthanumBromideCrystalRadius() {
	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	G4double out = pDetectionSystemLanthanumBromide->GetCrystalRadius();
	return out;
}

G4double DetectorConstruction::GetLanthanumBromideCrystalLength() {
	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	G4double out = pDetectionSystemLanthanumBromide->GetCrystalLength();
	return out;
}

G4double DetectorConstruction::GetLanthanumBromideR() {
	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	G4double out = pDetectionSystemLanthanumBromide->GetR();
	return out;
}

G4double DetectorConstruction::GetLanthanumBromideTheta(G4int i) {
	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	G4double out = pDetectionSystemLanthanumBromide->GetTheta(i);
	return out;
}

G4double DetectorConstruction::GetLanthanumBromidePhi(G4int i) {
	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	G4double out = pDetectionSystemLanthanumBromide->GetPhi(i);
	return out;
}
G4double DetectorConstruction::GetLanthanumBromideYaw(G4int i) {
	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	G4double out = pDetectionSystemLanthanumBromide->GetYaw(i);
	return out;
}

G4double DetectorConstruction::GetLanthanumBromidePitch(G4int i) {
	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	G4double out = pDetectionSystemLanthanumBromide->GetPitch(i);
	return out;
}

G4double DetectorConstruction::GetLanthanumBromideRoll(G4int i) {
	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	G4double out = pDetectionSystemLanthanumBromide->GetRoll(i);
	return out;
}

G4double DetectorConstruction::GetLanthanumBromideCrystalRadialPosition() {
	DetectionSystemLanthanumBromide* pDetectionSystemLanthanumBromide = new DetectionSystemLanthanumBromide();
	G4double out = pDetectionSystemLanthanumBromide->GetCrystalRadialPosition();
	return out;
}


// Temporary Function for testing purposes
void DetectorConstruction::AddDetectionSystemGriffinCustomDetector(G4int) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	fGriffinDetectorsMap[fGriffinDetectorsMapIndex] = fCustomDetectorNumber ;
	fGriffinDetectorsMapIndex++;

	// NOTE: ndet served no purpose in this case but I left it in just in case this needs to be modified later. The position of a detector placed using this function must be set using
	// SetDeadLayer.

	DetectionSystemGriffin* pGriffinCustom = new DetectionSystemGriffin(fExtensionSuppressorLocation, fDetectorShieldSelect, fDetectorRadialDistance, fHevimetSelector); // Select Forward (0) or Back (1)

	// check for each crystal if a custom dead layer has been specified
	for(G4int crystal = 0; crystal < 4; ++crystal) {
		if(fGriffinDeadLayer[fCustomDetectorNumber-1][crystal] >= 0.) {
			pGriffinCustom->SetDeadLayer(fCustomDetectorNumber-1, crystal, fGriffinDeadLayer[fCustomDetectorNumber-1][crystal]);
		}
	}

	pGriffinCustom->BuildDeadLayerSpecificCrystal(fCustomDetectorNumber-1);

	pGriffinCustom->PlaceDeadLayerSpecificCrystal(fLogicWorld, fCustomDetectorNumber-1, fCustomDetectorPosition-1, fUseTigressPositions);

	pGriffinCustom->BuildEverythingButCrystals();

	pGriffinCustom->PlaceEverythingButCrystals(fLogicWorld, fCustomDetectorNumber-1, fCustomDetectorPosition-1, fUseTigressPositions);

	fGriffin = true;
}

void DetectorConstruction::AddDetectionSystemGriffinCustom(G4int ndet) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	G4int detNum;
	G4int posNum;

	for(detNum = 1; detNum <= ndet; detNum++) {
		posNum = detNum;

		fGriffinDetectorsMap[fGriffinDetectorsMapIndex] = detNum;
		fGriffinDetectorsMapIndex++;

		DetectionSystemGriffin* pGriffinCustom = new DetectionSystemGriffin(fExtensionSuppressorLocation, fDetectorShieldSelect, fDetectorRadialDistance, fHevimetSelector) ; // Select Forward (0) or Back (1)

		pGriffinCustom->BuildDeadLayerSpecificCrystal(detNum-1);
		pGriffinCustom->PlaceDeadLayerSpecificCrystal(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;
		pGriffinCustom->BuildEverythingButCrystals();
		pGriffinCustom->PlaceEverythingButCrystals(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;

	}

	fGriffin = true;
}

void DetectorConstruction::AddDetectionSystemGriffinShieldSelect(G4int ShieldSelect){
	fDetectorShieldSelect = ShieldSelect ;
}

void DetectorConstruction::AddDetectionSystemGriffinSetRadialDistance(G4double detectorDist){
	fDetectorRadialDistance = detectorDist ;
}

void DetectorConstruction::AddDetectionSystemGriffinSetExtensionSuppLocation(G4int detectorPos){
	fExtensionSuppressorLocation = detectorPos ;
}

void DetectorConstruction::AddDetectionSystemGriffinSetPosition(G4ThreeVector params) {
	G4int detNum = static_cast<G4int>(params.x());
	G4int detPos = static_cast<G4int>(params.y());
	if(detNum < 1 || detNum > 16 || detPos < 1 || detPos > 16) {
		G4cout<<"Warning, can't set detector "<<detNum<<" to position "<<detPos<<", indices out of range!"<<G4endl;
		return;
	}
	fCustomDetectorNumber    = detNum;
	fCustomDetectorPosition  = detPos;
}

void DetectorConstruction::AddDetectionSystemGriffinSetDeadLayer(G4ThreeVector params) {
	/// x = detector number, y = crystal number, z = dead layer thickness
	G4int detNum = static_cast<G4int>(params.x());
	G4int cryNum = static_cast<G4int>(params.y());
	if(detNum < 0 || detNum > 15 || cryNum < 0 || cryNum > 3) {
		G4cout<<"Warning, can't set dead layer for detector "<<detNum<<", crystal "<<cryNum<<", indices out of range!"<<G4endl;
		return;
	}
	fGriffinDeadLayer[detNum][cryNum] = params.z();
}

void DetectorConstruction::AddDetectionSystemGriffinForward(G4int ndet) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	G4int detNum;
	G4int posNum;
	G4int config  = 0;
	for(detNum = 1; detNum <= ndet; detNum++) {
		posNum = detNum;

		fGriffinDetectorsMap[fGriffinDetectorsMapIndex] = detNum;
		fGriffinDetectorsMapIndex++;

		DetectionSystemGriffin* pGriffinDLS = new DetectionSystemGriffin(config, 1, fGriffinFwdBackPosition, fHevimetSelector); // Select Forward (0) or Back (1)

		pGriffinDLS->BuildDeadLayerSpecificCrystal(detNum-1);
		pGriffinDLS->PlaceDeadLayerSpecificCrystal(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;
		pGriffinDLS->BuildEverythingButCrystals();
		pGriffinDLS->PlaceEverythingButCrystals(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;
	}

	fGriffin = true;
}

void DetectorConstruction::AddDetectionSystemGriffinForwardDetector(G4int ndet) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	G4int detNum = ndet;
	G4int posNum = ndet;
	G4int config  = 0;
	fGriffinDetectorsMap[fGriffinDetectorsMapIndex] = detNum;
	fGriffinDetectorsMapIndex++;

	DetectionSystemGriffin* pGriffinDLS = new DetectionSystemGriffin(config, 1, fGriffinFwdBackPosition, fHevimetSelector); // Select Forward (0) or Back (1)

	pGriffinDLS->BuildDeadLayerSpecificCrystal(detNum-1);

	pGriffinDLS->PlaceDeadLayerSpecificCrystal(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;

	pGriffinDLS->BuildEverythingButCrystals();

	pGriffinDLS->PlaceEverythingButCrystals(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;

	fGriffin = true;
}

void DetectorConstruction::AddDetectionSystemGriffinBack(G4int ndet) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	G4int detNum;
	G4int posNum;
	G4int config  = 1;

	for(detNum = 1; detNum <= ndet; detNum++) {
		posNum = detNum;

		fGriffinDetectorsMap[fGriffinDetectorsMapIndex] = detNum;
		fGriffinDetectorsMapIndex++;

		DetectionSystemGriffin* pGriffinDLS = new DetectionSystemGriffin(config, 1, fGriffinFwdBackPosition, fHevimetSelector); // Select Forward (0) or Back (1)

		pGriffinDLS->BuildDeadLayerSpecificCrystal(detNum-1);
		pGriffinDLS->PlaceDeadLayerSpecificCrystal(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;
		pGriffinDLS->BuildEverythingButCrystals();
		pGriffinDLS->PlaceEverythingButCrystals(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;
	}

	fGriffin = true;
}

void DetectorConstruction::AddDetectionSystemGriffinBackDetector(G4int ndet) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	G4int detNum = ndet;
	G4int posNum = ndet;
	G4int config  = 1;

	fGriffinDetectorsMap[fGriffinDetectorsMapIndex] = detNum;
	fGriffinDetectorsMapIndex++;

	DetectionSystemGriffin* pGriffinDLS = new DetectionSystemGriffin(config, 1, fGriffinFwdBackPosition, fHevimetSelector); // Select Forward (0) or Back (1)

	pGriffinDLS->BuildDeadLayerSpecificCrystal(detNum-1);
	pGriffinDLS->PlaceDeadLayerSpecificCrystal(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;
	pGriffinDLS->BuildEverythingButCrystals();
	pGriffinDLS->PlaceEverythingButCrystals(fLogicWorld, detNum-1, posNum-1, fUseTigressPositions) ;

	fGriffin = true;
}

void DetectorConstruction::AddDetectionSystemGriffinHevimet(G4int input) {
	// Includes hevimet.
	fHevimetSelector = input ;

}

void DetectorConstruction::AddDetectionSystemSceptar(G4int ndet) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	DetectionSystemSceptar* pSceptar = new DetectionSystemSceptar() ;
	pSceptar->Build() ;
	pSceptar->PlaceDetector(fLogicWorld, ndet) ;

	fSceptar = true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DetectorConstruction::AddDetectionSystemDescant(G4int ndet) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	DetectionSystemDescant* pDetectionSystemDescant = new DetectionSystemDescant(true) ;
	pDetectionSystemDescant->Build() ;
	pDetectionSystemDescant->PlaceDetector(fLogicWorld, ndet) ;

	fDescant = true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DetectorConstruction::AddDetectionSystemDescantAuxPorts(G4ThreeVector input) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	G4int ndet = G4int(input.x());
	G4double radialpos = input.y()*cm;
	G4int leadShield = G4bool(input.z());

	if(leadShield) G4cout<<"Building DESCANT detectors WITH lead shielding"<<G4endl;
	if(!leadShield)G4cout<<"Building DESCANT detectors WITHOUT lead shielding"<<G4endl;

	DetectionSystemDescant *  pDetectionSystemDescant = new DetectionSystemDescant(leadShield) ;
	pDetectionSystemDescant->Build() ;

	for(G4int detectorNumber = 0; detectorNumber < ndet; detectorNumber++)
	{
		pDetectionSystemDescant->PlaceDetectorAuxPorts(fLogicWorld, detectorNumber, radialpos);
	}

	fDescant = true;
}

void DetectorConstruction::SetDetectionSystemDescantRotation(G4ThreeVector input) {
	//convert from degree to rad
	fDescantRotation.setX(input.x()/180.*M_PI);
	fDescantRotation.setY(input.y()/180.*M_PI);
	fDescantRotation.setZ(input.z()/180.*M_PI);
}

void DetectorConstruction::SetDetectionSystemDescantColor(G4String input) {
	fDescantColor = input;
}

void DetectorConstruction::AddDetectionSystemDescantCart(G4ThreeVector input) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	DetectionSystemDescant* pDetectionSystemDescant = new DetectionSystemDescant(true) ;
	pDetectionSystemDescant->Build() ;
	pDetectionSystemDescant->PlaceDetector(fLogicWorld, fDescantColor, input, fDescantRotation) ;

	fDescant = true;
}

void DetectorConstruction::AddDetectionSystemDescantSpher(G4ThreeVector input, G4double unit) {
	G4ThreeVector sphericalVector;
	//convert from degree to rad
	sphericalVector.setRThetaPhi(input.x()*unit, input.y()/180.*M_PI, input.z()/180.*M_PI);
	AddDetectionSystemDescantCart(sphericalVector);

	fDescant = true;
}

void DetectorConstruction::AddApparatusDescantStructure() {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	ApparatusDescantStructure* pApparatusDescantStructure = new ApparatusDescantStructure() ;
	pApparatusDescantStructure->Build() ;
	pApparatusDescantStructure->PlaceDescantStructure(fLogicWorld);

	//pApparatusDescantStructure->PlaceDetector(fLogicWorld, ndet) ;
}

void DetectorConstruction::AddDetectionSystemTestcan(G4ThreeVector input) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	G4double length = G4double(input.x())*cm;
	G4double radius = G4double(input.y())*cm;

	DetectionSystemTestcan* pDetectionSystemTestcan = new DetectionSystemTestcan(length, radius);
	pDetectionSystemTestcan->Build();
	pDetectionSystemTestcan->PlaceDetector(fLogicWorld);

	fTestcan = true;
}

void DetectorConstruction::AddDetectionSystemSpice() {
	if(fLogicWorld == nullptr) {
		Construct();
	}
	DetectionSystemSpice* pSpice = new DetectionSystemSpice() ;
	pSpice->BuildPlace(fLogicWorld); 

	fSpice = true;
	//HistoManager::Instance().Spice(true);//boolean needed to make SPICE histograms
}

void DetectorConstruction::AddDetectionSystemTrific(G4double Torr) {
	if(fLogicWorld == nullptr) {
		Construct();
	}
	DetectionSystemTrific* pTrific = new DetectionSystemTrific(Torr) ;
	pTrific->BuildPlace(fLogicWorld); 

}

void DetectorConstruction::AddDetectionSystemPaces(G4int ndet) {
	if(fLogicWorld == nullptr) {
		Construct();
	}

	DetectionSystemPaces* pPaces = new DetectionSystemPaces() ;
	pPaces->Build() ;

	pPaces->PlaceDetector(fLogicWorld, ndet) ;

	fPaces = true;
}

void DetectorConstruction::SetProperties() {
	// loop over all existing daughters of the world volume
	// check if their properties are set and if not, set them
	if(fLogicWorld == nullptr) return;
	// thread id is -1 for master, -2 in sequential mode
	// so this only outputs the number of volumes once
	if(G4Threading::G4GetThreadId() < 0) {
		G4cout<<fLogicWorld->GetNoDaughters()<<" daughter volumes"<<std::endl;
	}
	for(int i = 0; i < fLogicWorld->GetNoDaughters(); ++i) {
		if(!HasProperties(fLogicWorld->GetDaughter(i)) && CheckVolumeName(fLogicWorld->GetDaughter(i)->GetName())) {
			fPropertiesMap[fLogicWorld->GetDaughter(i)] = ParseVolumeName(fLogicWorld->GetDaughter(i)->GetName());
		}
	}
}

void DetectorConstruction::Print() {
	if(fLogicWorld == nullptr) {
		std::cout<<"World volume does not exist yet!"<<std::endl;
		return;
	}

	SetProperties();

	for(int i = 0; i < fLogicWorld->GetNoDaughters(); ++i) {
		std::cout<<i<<": "<<fLogicWorld->GetDaughter(i)<<" - "<<fLogicWorld->GetDaughter(i)->GetName();
		if(HasProperties(fLogicWorld->GetDaughter(i))) {
			auto prop = GetProperties(fLogicWorld->GetDaughter(i));
			std::cout<<" - "<<prop.detectorNumber<<", "<<prop.crystalNumber<<", "<<prop.systemID;
		}
		std::cout<<std::endl;
	}
}

bool DetectorConstruction::CheckVolumeName(G4String volumeName) {
	if(volumeName.find("germaniumBlock1") != G4String::npos) return true;
	if(volumeName.find("backQuarterSuppressor") != G4String::npos) return true;
	if(volumeName.find("leftSuppressorExtension") != G4String::npos) return true;
	if(volumeName.find("rightSuppressorExtension") != G4String::npos) return true;
	if(volumeName.find("leftSuppressorCasing") != G4String::npos) return true;
	if(volumeName.find("rightSuppressorCasing") != G4String::npos) return true;
	if(volumeName.find("germaniumDlsBlock1") != G4String::npos) return true;
	if(volumeName.find("lanthanumBromideCrystalBlock") != G4String::npos) return true;
	if(volumeName.find("ancillaryBgoBlock") != G4String::npos) return true;
	if(volumeName.find("sodiumIodideCrystalBlock") != G4String::npos) return true;
	if(volumeName.find("SiSegmentPhys") != G4String::npos) return true;
	if(volumeName.find("TrificGasCell") != G4String::npos) return true;
	if(volumeName.find("pacesSiliconBlockLog") != G4String::npos) return true;
	if(volumeName.find("sceptarSquareScintillatorLog") != G4String::npos) return true;
	if(volumeName.find("sceptarAngledScintillatorLog") != G4String::npos) return true;
	if(volumeName.find("8piGermaniumBlockLog") != G4String::npos) return true;
	if(volumeName.find("8piInnerBGOAnnulus") != G4String::npos) return true;
	if(volumeName.find("8piOuterLowerBGOAnnulus") != G4String::npos) return true;
	if(volumeName.find("8piOuterUpperBGOAnnulus") != G4String::npos) return true;
	if(volumeName.find("gridcellLog") != G4String::npos) return true;
	if(volumeName.find("blueScintillatorVolumeLog") != G4String::npos) return true;
	if(volumeName.find("greenScintillatorVolumeLog") != G4String::npos) return true;
	if(volumeName.find("redScintillatorVolumeLog") != G4String::npos) return true;
	if(volumeName.find("whiteScintillatorVolumeLog") != G4String::npos) return true;
	if(volumeName.find("yellowScintillatorVolumeLog") != G4String::npos) return true;
	if(volumeName.find("testcanScintillatorLog") != G4String::npos) return true;
	return false;
}

DetectorProperties DetectorConstruction::ParseVolumeName(G4String volumeName) {
	DetectorProperties result;
	// GRIFFIN detectors have the detector and crystal number in their names
	if(volumeName.find("germaniumBlock1") != G4String::npos) {
		// strip "germaniumBlock1_" (16 characters) and everything before from the string
		std::string tmpString = volumeName.substr(volumeName.find("germaniumBlock1")+16);
		// replace all '_' with spaces so we can just use istringstream::operator>>
		std::replace(tmpString.begin(), tmpString.end(), '_', ' ');
		// create istringstream from the stripped and converted stream, and read detector and crystal number
		std::istringstream is(tmpString);
		is>>result.detectorNumber>>result.crystalNumber;
		// converting this number to a "true" detector number isn't necessary anymore since we use the real number and not the assembly/imprint number
		result.systemID = 1000;
		return result;
	}

	if(volumeName.find("germaniumDlsBlock1") != G4String::npos) {
		// strip "germaniumDlsBlock1_" (16 characters) and everything before from the string
		std::string tmpString = volumeName.substr(volumeName.find("germaniumDlsBlock1")+19);
		// replace all '_' with spaces so we can just use istringstream::operator>>
		std::replace(tmpString.begin(), tmpString.end(), '_', ' ');
		// create istringstream from the stripped and converted stream, and read detector and crystal number
		std::istringstream is(tmpString);
		is>>result.detectorNumber>>result.crystalNumber;
		// converting this number to a "true" detector number isn't necessary anymore since we use the real number and not the assembly/imprint number
		result.systemID = 1000;
		return result;
	}

	// Spice detectors also have the number in their name
	if(volumeName.find("SiSegmentPhys") != G4String::npos) {
		// strip "SiSegmentPhys_" (13 characters) and everything before from the string
		std::string tmpString = volumeName.substr(volumeName.find("SiSegmentPhys")+13);
		// replace all '_' with spaces so we can just use istringstream::operator>>
		std::replace(tmpString.begin(), tmpString.end(), '_', ' ');
		// create istringstream from the stripped and converted stream, and read detector and crystal number
		std::istringstream is(tmpString);
		is>>result.detectorNumber>>result.crystalNumber;
		result.systemID = 10;
		return result;
	}

	// Trific detectors have the number in their name too
	if(volumeName.find("TrificGasCell") != G4String::npos) {
		// strip "TrificGasCell_" (13 characters) and everything before from the string
		std::string tmpString = volumeName.substr(volumeName.find("TrificGasCell")+13);
		// replace all '_' with spaces so we can just use istringstream::operator>>
		std::replace(tmpString.begin(), tmpString.end(), '_', ' ');
		// create istringstream from the stripped and converted stream, and read detector and crystal number
		std::istringstream is(tmpString);
		is>>result.detectorNumber;
		result.systemID = 10;//?? why the same as SiSegmentPhys??
		return result;
	}

	// for all other detectors we need to get the assembly and imprint number
	// the name has the form av_<assembly volume number>_impr_<imprint number>_<volume name>_pv_<??? number>
	if(volumeName.find("av_") != 0 || volumeName.find("_impr_") == G4String::npos) {
		std::cerr<<"wrongly formed volume name '"<<volumeName<<"', should start with 'av_' and have '_impr_' in it?"<<std::endl;
		throw std::invalid_argument("unexpected volume name");
	}
	// create new string that starts with the assembly number and use stringstream to read it
	std::string tmpString = volumeName.substr(3);
	// replace all '_' with spaces so we can just use istringstream::operator>>
	std::replace(tmpString.begin(), tmpString.end(), '_', ' ');
	std::istringstream is(tmpString);
	G4int assemblyNumber;
	is>>assemblyNumber;

	// create new string that starts with the imprint number ('_' have been replace by ' ') and use stringstream to read it
	tmpString = tmpString.substr(tmpString.find(" impr ")+6);
	is.str(tmpString);
	G4int imprintNumber;
	is>>imprintNumber;

	// for griffin "components" (aka suppressors) we get the detector number from the assembly volume number (av-5)/fNumberOfAssemblyVols
	// and the crystal number is the imprint number (fNumberOfAssemblyVols was hard-coded to be 13 in the constructor)
	result.detectorNumber = static_cast<G4int>(ceil((assemblyNumber-5.)/13.));
	result.crystalNumber = imprintNumber;
	if(volumeName.find("backQuarterSuppressor") != G4String::npos) {
		result.systemID = 1050;
		return result;
	}

	if(volumeName.find("leftSuppressorExtension") != G4String::npos) {
		result.systemID = 1010;
		return result;
	}

	if(volumeName.find("rightSuppressorExtension") != G4String::npos) {
		result.systemID = 1020;
		return result;
	}

	if(volumeName.find("leftSuppressorCasing") != G4String::npos) {
		result.systemID = 1030;
		return result;
	}

	if(volumeName.find("rightSuppressorCasing") != G4String::npos) {
		result.systemID = 1040;
		return result;
	}

	// for ancillary BGOs the detector number is the (ceiling of the) imprint number divided by 3
	// and the crystal number is the imprint number minus 3 times the detector number minus one
	result.detectorNumber = static_cast<G4int>(ceil(imprintNumber/3.0));
	result.crystalNumber = imprintNumber - (result.detectorNumber-1)*3;
	if(volumeName.find("ancillaryBgoBlock") != G4String::npos) {
		result.systemID = 3000;
		return result;
	}

	// for "generic" detectors the detector number is the imprint number
	result.detectorNumber = imprintNumber;
	result.crystalNumber = 0;
	if(volumeName.find("lanthanumBromideCrystalBlock") != G4String::npos) {
		result.systemID = 2000;
		return result;
	}

	if(volumeName.find("sodiumIodideCrystalBlock") != G4String::npos) {
		result.systemID = 4000;
		return result;
	}

	if(volumeName.find("pacesSiliconBlockLog") != G4String::npos) {
		result.systemID = 50;
		return result;
	}

	if(volumeName.find("sceptarSquareScintillatorLog") != G4String::npos) {
		// to number SCEPTAR paddles correctly:
		switch(result.detectorNumber) {
			case 0:
				result.detectorNumber = 5;
				break;
			case 1:
				result.detectorNumber = 9;
				break;
			case 2:
				result.detectorNumber = 8;
				break;
			case 3:
				result.detectorNumber = 7;
				break;
			case 4:
				result.detectorNumber = 6;
				break;
			case 5:
				result.detectorNumber = 13;
				break;
			case 6:
				result.detectorNumber = 12;
				break;
			case 7:
				result.detectorNumber = 11;
				break;
			case 8:
				result.detectorNumber = 10;
				break;
			case 9:
				result.detectorNumber = 14;
				break;
			default:
				std::cerr<<"Unknown detector number "<<result.detectorNumber<<" for square SCEPTAR!"<<std::endl;
				break;
		}
		result.systemID = 5000;
		return result;
	}

	if(volumeName.find("sceptarAngledScintillatorLog") != G4String::npos) {
		// to number SCEPTAR paddles correctly (1 stays 1):
		switch(result.detectorNumber) {
			case 0:
				result.detectorNumber = 0;
				break;
			case 1:
				result.detectorNumber = 4;
				break;
			case 2:
				result.detectorNumber = 3;
				break;
			case 3:
				result.detectorNumber = 2;
				break;
			case 4:
				result.detectorNumber = 1;
				break;
			case 5:
				result.detectorNumber = 18;
				break;
			case 6:
				result.detectorNumber = 17;
				break;
			case 7:
				result.detectorNumber = 16;
				break;
			case 8:
				result.detectorNumber = 15;
				break;
			case 9:
				result.detectorNumber = 19;
				break;
			default:
				std::cerr<<"Unknown detector number "<<result.detectorNumber<<" for angled SCEPTAR!"<<std::endl;
				break;
		}
		result.systemID = 5000;
		return result;
	}

	if(volumeName.find("8piGermaniumBlockLog") != G4String::npos) {
		result.systemID = 6000;
		return result;
	}

	if(volumeName.find("8piInnerBGOAnnulus") != G4String::npos) {
		result.systemID = 6010;
		return result;
	}

	if(volumeName.find("8piOuterLowerBGOAnnulus") != G4String::npos) {
		result.systemID = 6020;
		return result;
	}

	if(volumeName.find("8piOuterUpperBGOAnnulus") != G4String::npos) {
		result.systemID = 6030;
		return result;
	}

	if(volumeName.find("gridcellLog") != G4String::npos) {
		result.systemID = 7000;
		return result;
	}

	if(volumeName.find("blueScintillatorVolumeLog") != G4String::npos) {
		result.systemID = 8010;
		return result;
	}

	if(volumeName.find("greenScintillatorVolumeLog") != G4String::npos) {
		result.systemID = 8020;
		return result;
	}

	if(volumeName.find("redScintillatorVolumeLog") != G4String::npos) {
		result.systemID = 8030;
		return result;
	}

	if(volumeName.find("whiteScintillatorVolumeLog") != G4String::npos) {
		result.systemID = 8040;
		return result;
	}

	if(volumeName.find("yellowScintillatorVolumeLog") != G4String::npos) {
		result.systemID = 8050;
		return result;
	}

	if(volumeName.find("testcanScintillatorLog") != G4String::npos) {
		result.systemID = 8500;
		return result;
	}

	return result;
}

