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

DetectionSystemBox::DetectionSystemBox(G4double x_length_in, G4double y_length_in, G4double z_length_in, G4double thickness_in, G4String box_mat, G4ThreeVector box_colour) :
    // LogicalVolumes
    cell_log(0),
    side_norm_x_log(0),
    side_norm_y_log(0),
    side_norm_z_log(0)
{ 

    x_length = x_length_in;
    y_length = y_length_in;
    z_length = z_length_in;

    thickness = thickness_in;

    cell_material = box_mat;

    cell_colour = box_colour;
}

DetectionSystemBox::~DetectionSystemBox()
{
    // LogicalVolumes
    delete cell_log;
    delete side_norm_x_log;
    delete side_norm_y_log;
    delete side_norm_z_log;

    //    delete crystal_block_SD;
}

//G4int DetectionSystemBox::Build(G4SDManager* mySDman)
G4int DetectionSystemBox::Build()
{ 
    // Build assembly volume
    G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    this->assembly = myAssembly;

    G4cout << "BuildXNormalVolumes" << G4endl;
    BuildXNormalVolumes();
    BuildYNormalVolumes();
    BuildZNormalVolumes();

    return 1;
}

G4int DetectionSystemBox::PlaceDetector(G4LogicalVolume* exp_hall_log)
{
    G4int cell_number = 1;

    G4ThreeVector position = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
    G4RotationMatrix* rotate  = new G4RotationMatrix;

    assembly->MakeImprint(exp_hall_log, position, rotate, cell_number);

    return 1;
}

G4int DetectionSystemBox::BuildXNormalVolumes()
{
    G4Material* material = G4Material::GetMaterial(this->cell_material);
    if( !material ) {
        G4cout << " ----> Material " << this->cell_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(cell_colour));
    vis_att->SetVisibility(true);


    G4Box* side_norm_x = BuildXNormalSide();

    // place first side
    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(1,0,0);
    G4double z_position		= (x_length+thickness)/2.0;
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate  = new G4RotationMatrix;

    //logical volume
    if( side_norm_x_log == NULL )
    {
        side_norm_x_log = new G4LogicalVolume(side_norm_x, material, "side_norm_x_log", 0, 0, 0);
        side_norm_x_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(side_norm_x_log, move, rotate);

    // place second side
    // Define rotation and movement objects
    direction     = G4ThreeVector(-1,0,0);
    z_position    = (x_length+thickness)/2.0;
    move          = z_position * direction;

    this->assembly->AddPlacedVolume(side_norm_x_log, move, rotate);

    return 1;
}

G4int DetectionSystemBox::BuildYNormalVolumes()
{
    G4Material* material = G4Material::GetMaterial(this->cell_material);
    if( !material ) {
        G4cout << " ----> Material " << this->cell_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(cell_colour));
    vis_att->SetVisibility(true);


    G4Box* side_norm_y = BuildYNormalSide();

    // place first side
    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,1,0);
    G4double z_position		= (y_length+thickness)/2.0;
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate  = new G4RotationMatrix;

    //logical volume
    if( side_norm_y_log == NULL )
    {
        side_norm_y_log = new G4LogicalVolume(side_norm_y, material, "side_norm_y_log", 0, 0, 0);
        side_norm_y_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(side_norm_y_log, move, rotate);

    // place second side
    // Define rotation and movement objects
    direction     = G4ThreeVector(0,-1,0);
    z_position    = (y_length+thickness)/2.0;
    move          = z_position * direction;

    this->assembly->AddPlacedVolume(side_norm_y_log, move, rotate);

    return 1;
}

G4int DetectionSystemBox::BuildZNormalVolumes()
{
    G4Material* material = G4Material::GetMaterial(this->cell_material);
    if( !material ) {
        G4cout << " ----> Material " << this->cell_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(cell_colour));
    vis_att->SetVisibility(true);


    G4Box* side_norm_z = BuildZNormalSide();

    // place first side
    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= (z_length+thickness)/2.0;
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate  = new G4RotationMatrix;

    //logical volume
    if( side_norm_z_log == NULL )
    {
        side_norm_z_log = new G4LogicalVolume(side_norm_z, material, "side_norm_z_log", 0, 0, 0);
        side_norm_z_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(side_norm_z_log, move, rotate);

    // place second side
    // Define rotation and movement objects
    direction     = G4ThreeVector(0,0,-1);
    z_position    = (z_length+thickness)/2.0;
    move          = z_position * direction;

    this->assembly->AddPlacedVolume(side_norm_z_log, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemBox::BuildXNormalSide()
{
    G4double half_length_x = thickness/2.0;
    G4double half_length_y = (y_length+(2.0*thickness))/2.0;
    G4double half_length_z = (z_length+(2.0*thickness))/2.0;

    G4Box* cellX = new G4Box("cellX", half_length_x, half_length_y, half_length_z);

    return cellX;
}//end ::BuildXNormalSide

G4Box* DetectionSystemBox::BuildYNormalSide()
{
    G4double half_length_x = (x_length)/2.0;
    G4double half_length_y = thickness/2.0;
    G4double half_length_z = (z_length+(2.0*thickness))/2.0;

    G4Box* cellY = new G4Box("cellY", half_length_x, half_length_y, half_length_z);

    return cellY;
}//end ::BuildXNormalSide

G4Box* DetectionSystemBox::BuildZNormalSide()
{
    G4double half_length_x = (x_length)/2.0;
    G4double half_length_y = (y_length)/2.0;
    G4double half_length_z = thickness/2.0;

    G4Box* cellZ = new G4Box("cellZ", half_length_x, half_length_y, half_length_z);

    return cellZ;
}//end ::BuildXNormalSide
