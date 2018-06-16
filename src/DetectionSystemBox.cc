#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemBox.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

DetectionSystemBox::DetectionSystemBox(G4double xLengthIn, G4double yLengthIn, G4double zLengthIn, G4double thicknessIn, G4String boxMat, G4ThreeVector boxColour) :
	// LogicalVolumes
	fCellLog(0),
	fSideNormXLog(0),
	fSideNormYLog(0),
	fSideNormZLog(0)
{ 

	fXLength = xLengthIn;
	fYLength = yLengthIn;
	fZLength = zLengthIn;

	fThickness = thicknessIn;

	fCellMaterial = boxMat;

	fCellColour = boxColour;
}

DetectionSystemBox::~DetectionSystemBox()
{
	// LogicalVolumes
	delete fCellLog;
	delete fSideNormXLog;
	delete fSideNormYLog;
	delete fSideNormZLog;

	//    delete crystalBlockSD;
}

//G4int DetectionSystemBox::Build(G4SDManager* mySDman)
G4int DetectionSystemBox::Build()
{ 
	// Build assembly volume
	fAssembly = new G4AssemblyVolume();

	G4cout<<"BuildXNormalVolumes"<<G4endl;
	BuildXNormalVolumes();
	BuildYNormalVolumes();
	BuildZNormalVolumes();

	return 1;
}

G4int DetectionSystemBox::PlaceDetector(G4LogicalVolume* expHallLog)
{
	G4int cellNumber = 1;

	G4ThreeVector position = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	G4RotationMatrix* rotate  = new G4RotationMatrix;

	fAssembly->MakeImprint(expHallLog, position, rotate, cellNumber);

	return 1;
}

G4int DetectionSystemBox::BuildXNormalVolumes()
{
	G4Material* material = G4Material::GetMaterial(fCellMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fCellMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(fCellColour));
	visAtt->SetVisibility(true);


	G4Box* sideNormX = BuildXNormalSide();

	// place first side
	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(1,0,0);
	G4double zPosition		= (fXLength+fThickness)/2.0;
	G4ThreeVector move 		= zPosition * direction;
	G4RotationMatrix* rotate  = new G4RotationMatrix;

	//logical volume
	if(fSideNormXLog == nullptr)
	{
		fSideNormXLog = new G4LogicalVolume(sideNormX, material, "sideNormXLog", 0, 0, 0);
		fSideNormXLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fSideNormXLog, move, rotate);

	// place second side
	// Define rotation and movement objects
	direction     = G4ThreeVector(-1,0,0);
	zPosition    = (fXLength+fThickness)/2.0;
	move          = zPosition * direction;

	fAssembly->AddPlacedVolume(fSideNormXLog, move, rotate);

	return 1;
}

G4int DetectionSystemBox::BuildYNormalVolumes()
{
	G4Material* material = G4Material::GetMaterial(fCellMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fCellMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(fCellColour));
	visAtt->SetVisibility(true);


	G4Box* sideNormY = BuildYNormalSide();

	// place first side
	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,1,0);
	G4double zPosition		= (fYLength+fThickness)/2.0;
	G4ThreeVector move 		= zPosition * direction;
	G4RotationMatrix* rotate  = new G4RotationMatrix;

	//logical volume
	if(fSideNormYLog == nullptr)
	{
		fSideNormYLog = new G4LogicalVolume(sideNormY, material, "sideNormYLog", 0, 0, 0);
		fSideNormYLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fSideNormYLog, move, rotate);

	// place second side
	// Define rotation and movement objects
	direction     = G4ThreeVector(0,-1,0);
	zPosition    = (fYLength+fThickness)/2.0;
	move          = zPosition * direction;

	fAssembly->AddPlacedVolume(fSideNormYLog, move, rotate);

	return 1;
}

G4int DetectionSystemBox::BuildZNormalVolumes()
{
	G4Material* material = G4Material::GetMaterial(fCellMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fCellMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(fCellColour));
	visAtt->SetVisibility(true);


	G4Box* sideNormZ = BuildZNormalSide();

	// place first side
	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,1);
	G4double zPosition		= (fZLength+fThickness)/2.0;
	G4ThreeVector move 		= zPosition * direction;
	G4RotationMatrix* rotate  = new G4RotationMatrix;

	//logical volume
	if(fSideNormZLog == nullptr) {
		fSideNormZLog = new G4LogicalVolume(sideNormZ, material, "sideNormZLog", 0, 0, 0);
		fSideNormZLog->SetVisAttributes(visAtt);
	}

	fAssembly->AddPlacedVolume(fSideNormZLog, move, rotate);

	// place second side
	// Define rotation and movement objects
	direction    = G4ThreeVector(0,0,-1);
	zPosition    = (fZLength+fThickness)/2.0;
	move         = zPosition * direction;

	fAssembly->AddPlacedVolume(fSideNormZLog, move, rotate);

	return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemBox::BuildXNormalSide()
{
	G4double halfLengthX = fThickness/2.0;
	G4double halfLengthY = (fYLength+(2.0*fThickness))/2.0;
	G4double halfLengthZ = (fZLength+(2.0*fThickness))/2.0;

	G4Box* cellX = new G4Box("cellX", halfLengthX, halfLengthY, halfLengthZ);

	return cellX;
}//end ::BuildXNormalSide

G4Box* DetectionSystemBox::BuildYNormalSide()
{
	G4double halfLengthX = (fXLength)/2.0;
	G4double halfLengthY = fThickness/2.0;
	G4double halfLengthZ = (fZLength+(2.0*fThickness))/2.0;

	G4Box* cellY = new G4Box("cellY", halfLengthX, halfLengthY, halfLengthZ);

	return cellY;
}//end ::BuildXNormalSide

G4Box* DetectionSystemBox::BuildZNormalSide()
{
	G4double halfLengthX = (fXLength)/2.0;
	G4double halfLengthY = (fYLength)/2.0;
	G4double halfLengthZ = fThickness/2.0;

	G4Box* cellZ = new G4Box("cellZ", halfLengthX, halfLengthY, halfLengthZ);

	return cellZ;
}//end ::BuildXNormalSide
