#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

//#include "FieldSetup.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

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

#include "DetectionSystemGrid.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


DetectionSystemGrid::DetectionSystemGrid(G4double xLengthIn, G4double yLengthIn, G4double zLengthIn, G4double cubeSizeIn, G4String gridMat, G4ThreeVector gridColour, G4ThreeVector posOffset) :
	// LogicalVolumes
	fGridcellLog(0)
{ 
	fXLength = xLengthIn;
	fYLength = yLengthIn;
	fZLength = zLengthIn;

	fXPosOffset = posOffset.x();
	fYPosOffset = posOffset.y();
	fZPosOffset = posOffset.z();

	fCubeSize = cubeSizeIn;

	fCellsXRow = (G4int)(fXLength/cubeSizeIn);
	fCellsYRow = (G4int)(fYLength/cubeSizeIn);
	fCellsZRow = (G4int)(fZLength/cubeSizeIn);

	fNumberOfCells = fCellsXRow*fCellsYRow*fCellsZRow;
	fCellsPerRow   = 10;
	fCellWidth      = fCubeSize;

	fCellMaterial = gridMat;

	fCellColour = gridColour;

	//fEmFieldSetup = theFieldSetup;

	G4cout<<"xLength = "<<fXLength/cm<<G4endl;
	G4cout<<"yLength = "<<fYLength/cm<<G4endl;
	G4cout<<"zLength = "<<fZLength/cm<<G4endl;

	G4cout<<"cellsXRow = "<<fCellsXRow<<G4endl;
}

DetectionSystemGrid::~DetectionSystemGrid() {
	// LogicalVolumes
	delete fGridcellLog;
}


G4int DetectionSystemGrid::Build() { 

	// Build assembly volume
	fAssembly = new G4AssemblyVolume(); 

	G4cout<<"BuildCellVolume"<<G4endl;
	BuildCellVolume();

	return 1;
}

G4int DetectionSystemGrid::PlaceDetector(G4LogicalVolume* expHallLog) {
	G4double startX = (-1.0*fXLength/2.0);
	G4double startY = (-1.0*fYLength/2.0);
	G4double startZ = (-1.0*fZLength/2.0);
	G4double posX;
	G4double posY;
	G4double posZ;

	G4int cellNumber = 0;

	G4ThreeVector position = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
	G4RotationMatrix* rotate  = new G4RotationMatrix;

	for(G4int rowX = 0; rowX < fCellsXRow; rowX++) {
		posX = startX+(fCellWidth*rowX)+(fCellWidth/2.0);
		position.setX(posX + fXPosOffset);

		for(G4int rowY = 0; rowY < fCellsYRow; rowY++) {
			posY = startY+(fCellWidth*rowY)+(fCellWidth/2.0);
			position.setY(posY + fYPosOffset);

			for(G4int rowZ = 0; rowZ < fCellsZRow; rowZ++) {
				posZ = startZ+(fCellWidth*rowZ)+(fCellWidth/2.0);
				position.setZ(posZ + fZPosOffset);

				cellNumber++;
				//G4cout<<"--------------------------- Grid Detector Number = "<<cellNumber<<G4endl;

				fAssembly->MakeImprint(expHallLog, position, rotate, cellNumber);
			}
		}
	}

	return 1;
}

G4int DetectionSystemGrid::BuildCellVolume() {
	G4Material* material = G4Material::GetMaterial(fCellMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fCellMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// Set visualization attributes
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(fCellColour));
	visAtt->SetVisibility(true);


	G4Box* gridcell = BuildCell();

	// Define rotation and movement objects
	G4ThreeVector direction 	= G4ThreeVector(0,0,0);
	G4double zPosition		= 0.0*mm;
	G4ThreeVector move 		= zPosition * direction;
	G4RotationMatrix* rotate  = new G4RotationMatrix;

	//logical volume
	if(fGridcellLog == nullptr) {
		fGridcellLog = new G4LogicalVolume(gridcell, material, "gridcellLog", 0, 0, 0);
		fGridcellLog->SetVisAttributes(visAtt);

		//// Set local field manager and local field in radiator and its daughters:
		//G4bool allLocal = true ;
		//gridcellLog->SetFieldManager(fEmFieldSetup->GetLocalFieldManager(), allLocal) ;
	}


	fAssembly->AddPlacedVolume(fGridcellLog, move, rotate);

	return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemGrid::BuildCell() {
	G4double halfLengthX = fCellWidth/2.0;
	G4double halfLengthY = fCellWidth/2.0;
	G4double halfLengthZ = fCellWidth/2.0;

	G4Box* cell = new G4Box("cell", halfLengthX, halfLengthY, halfLengthZ);

	return cell;
}//end ::BuildCell
