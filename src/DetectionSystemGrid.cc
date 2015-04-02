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


DetectionSystemGrid::DetectionSystemGrid(G4double x_length_in, G4double y_length_in, G4double z_length_in, G4double cube_size_in, G4String grid_mat, G4ThreeVector grid_colour, G4ThreeVector pos_offset) :
    // LogicalVolumes
    gridcell_log(0)
{ 
    x_length = x_length_in;
    y_length = y_length_in;
    z_length = z_length_in;

    x_pos_offset = pos_offset.x();
    y_pos_offset = pos_offset.y();
    z_pos_offset = pos_offset.z();

    cube_size = cube_size_in;

    cells_x_row = (G4int)(x_length/cube_size_in);
    cells_y_row = (G4int)(y_length/cube_size_in);
    cells_z_row = (G4int)(z_length/cube_size_in);

    number_of_cells = cells_x_row*cells_y_row*cells_z_row;
    cells_per_row   = 10;
    cell_width      = cube_size;

    cell_material = grid_mat;

    cell_colour = grid_colour;

    //fEmFieldSetup = theFieldSetup;

    G4cout << "x_length = " << x_length/cm << G4endl;
    G4cout << "y_length = " << y_length/cm << G4endl;
    G4cout << "z_length = " << z_length/cm << G4endl;

    G4cout << "cells_x_row = " << cells_x_row << G4endl;
}

DetectionSystemGrid::~DetectionSystemGrid()
{
    // LogicalVolumes
    delete gridcell_log;

}


G4int DetectionSystemGrid::Build()
{ 

    // Build assembly volume
    G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    this->assembly = myAssembly;

    G4cout << "BuildCellVolume" << G4endl;
    BuildCellVolume();

    return 1;
}

G4int DetectionSystemGrid::PlaceDetector(G4LogicalVolume* exp_hall_log)
{
    G4double start_x = (-1.0*x_length/2.0);
    G4double start_y = (-1.0*y_length/2.0);
    G4double start_z = (-1.0*z_length/2.0);
    G4double pos_x;
    G4double pos_y;
    G4double pos_z;

    G4int cell_number = 0;

    G4ThreeVector position = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
    G4RotationMatrix* rotate  = new G4RotationMatrix;

    for(G4int row_x = 0; row_x < cells_x_row; row_x++)
    {
        pos_x = start_x+(cell_width*row_x)+(cell_width/2.0);
        position.setX(pos_x + this->x_pos_offset);

        for(G4int row_y = 0; row_y < cells_y_row; row_y++)
        {
            pos_y = start_y+(cell_width*row_y)+(cell_width/2.0);
            position.setY(pos_y + this->y_pos_offset);

            for(G4int row_z = 0; row_z < cells_z_row; row_z++)
            {
                pos_z = start_z+(cell_width*row_z)+(cell_width/2.0);
                position.setZ(pos_z + this->z_pos_offset);

                cell_number++;
                //G4cout << "--------------------------- Grid Detector Number = " << cell_number << G4endl;

                assembly->MakeImprint(exp_hall_log, position, rotate, cell_number);
            }
        }
    }

    return 1;
}

G4int DetectionSystemGrid::BuildCellVolume()
{
    G4Material* material = G4Material::GetMaterial(this->cell_material);
    if( !material ) {
        G4cout << " ----> Material " << this->cell_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(cell_colour));
    vis_att->SetVisibility(true);


    G4Box* gridcell = BuildCell();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,0);
    G4double z_position		= 0.0*mm;
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate  = new G4RotationMatrix;

    //logical volume
    if( gridcell_log == NULL )
    {
        gridcell_log = new G4LogicalVolume(gridcell, material, "gridcell_log", 0, 0, 0);
        gridcell_log->SetVisAttributes(vis_att);

        //// Set local field manager and local field in radiator and its daughters:
        //G4bool allLocal = true ;
        //gridcell_log->SetFieldManager( fEmFieldSetup->GetLocalFieldManager(), allLocal ) ;
    }


    this->assembly->AddPlacedVolume(gridcell_log, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemGrid::BuildCell()
{
    G4double half_length_x = cell_width/2.0;
    G4double half_length_y = cell_width/2.0;
    G4double half_length_z = cell_width/2.0;

    G4Box* cell = new G4Box("cell", half_length_x, half_length_y, half_length_z);

    return cell;
}//end ::BuildCell
