#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemSodiumIodide.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

DetectionSystemSodiumIodide::DetectionSystemSodiumIodide() :
    // LogicalVolumes
    detector_volume_log(0),
    crystal_block_log(0),
    can_cylinder_log(0),
    can_front_lid_log(0),
    can_back_lid_log(0),
    seal_front_lid_log(0),
    disc_front_lid_log(0),
    packing_cylinder_log(0),
    packing_front_lid_log(0)
{ 

    // NaI Crystal and Can Physical Properties
    // From: Scintillation Spectrometry gamma-ray catalogue - R. L. Heath, 1964



    /////////////////////////////////////////////////////////////////////
    //  0.019"
    //  "can" - Al Container
    //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //% 0.040"
    //% "seal" - Compressed Neoprene Rubber
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //-------------------------------------------------------------------
    //- 0.006"
    //- "disc" - Polyethylene Disc
    //-------------------------------------------------------------------
    /////////////////////////////////////////////////////////////////////
    //- 0.0625"
    //- "packing" - Packed Aluminum Oxide
    /////////////////////////////////////////////////////////////////////
    //
    //
    //  NaI Crystal
    //
    //
    //
    //
    //


    this->detail_view_end_angle	 	= 360.0*deg;
    this->crystal_material            = "G4_SODIUM_IODIDE";
    this->can_material                = "G4_Al";
    this->seal_material               = "G4_RUBBER_NEOPRENE";
    this->disc_material               = "G4_POLYETHYLENE";
    this->packing_material            = "G4_ALUMINUM_OXIDE";

    this->crystal_length_z            = 3.0*inch2cm*cm;
    this->crystal_inner_radius 		= 0.0*cm;
    this->crystal_outer_radius 		= 1.5*inch2cm*cm;

    this->packing_length_z            = this->crystal_length_z;
    this->packing_inner_radius 		= this->crystal_outer_radius;
    this->packing_outer_radius 		= 0.0625*inch2cm*cm + this->crystal_outer_radius;

    this->packing_lid_inner_radius    = 0.0*cm;
    this->packing_lid_outer_radius 	= this->packing_outer_radius;
    this->packing_front_lid_thickness	= 0.0625*inch2cm*cm;

    this->disc_lid_inner_radius       = 0.0*cm;
    this->disc_lid_outer_radius       = this->packing_lid_outer_radius;
    this->disc_front_lid_thickness	= 0.006*inch2cm*cm;

    this->seal_lid_inner_radius       = 0.0*cm;
    this->seal_lid_outer_radius       = this->packing_lid_outer_radius;
    this->seal_front_lid_thickness	= 0.04*inch2cm*cm;

    this->can_length_z                = this->crystal_length_z +
            this->packing_front_lid_thickness +
            this->disc_front_lid_thickness +
            this->seal_front_lid_thickness;
    this->can_inner_radius            = this->packing_outer_radius;
    this->can_outer_radius            = 0.019*inch2cm*cm + this->packing_outer_radius;

    this->can_lid_inner_radius        = 0.0*cm;
    this->can_lid_outer_radius        = this->can_outer_radius;
    this->can_front_lid_thickness     = 0.019*inch2cm*cm;
    this->can_back_lid_thickness      = 0.019*inch2cm*cm;


    this->detector_length_z           = this->crystal_length_z +
            this->packing_front_lid_thickness +
            this->disc_front_lid_thickness +
            this->seal_front_lid_thickness +
            this->can_front_lid_thickness +
            this->can_back_lid_thickness;
}

DetectionSystemSodiumIodide::~DetectionSystemSodiumIodide()
{
    // LogicalVolumes
    delete detector_volume_log;
    delete crystal_block_log;
    delete can_cylinder_log;
    delete can_front_lid_log;
    delete can_back_lid_log;
    delete seal_front_lid_log;
    delete disc_front_lid_log;
    delete packing_cylinder_log;
    delete packing_front_lid_log;

}


G4int DetectionSystemSodiumIodide::Build()//G4SDManager* mySDman)
{

    // Build assembly volume
    G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    this->assembly = myAssembly;

    G4cout << "BuildCrystalVolume" << G4endl;
    BuildCrystalVolume();
    G4cout << "BuildAluminumCanVolume" << G4endl;
    BuildAluminumCanVolume();
    G4cout << "BuildPackingVolume" << G4endl;
    BuildPackingVolume();
    G4cout << "BuildDiscVolume" << G4endl;
    BuildDiscVolume();
    G4cout << "BuildSealVolume" << G4endl;
    BuildSealVolume();

    return 1;
}

G4int DetectionSystemSodiumIodide::PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate, G4int detector_number)
{
    G4int detector_copy_ID = 0;

    G4cout << "SodiumIodide Detector Number = " << detector_number << G4endl;

    copy_number = detector_copy_ID + detector_number;

    assembly->MakeImprint(exp_hall_log, move, rotate, copy_number);

    return 1;
}

G4int DetectionSystemSodiumIodide::BuildCrystalVolume()
{
    G4Material* material = G4Material::GetMaterial(this->crystal_material);
    if( !material ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.2,1.0,0.3));
    vis_att->SetVisibility(true);

    G4Tubs* crystal_block = BuildCrystal();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= ( (can_front_lid_thickness + seal_front_lid_thickness + disc_front_lid_thickness + packing_front_lid_thickness) - (can_back_lid_thickness) )/2.0;
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( crystal_block_log == NULL )
    {
        crystal_block_log = new G4LogicalVolume(crystal_block, material, "sodium_iodide_crystal_block_log", 0, 0, 0);
        crystal_block_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(crystal_block_log, move, rotate);

    return 1;
}

G4int DetectionSystemSodiumIodide::BuildAluminumCanVolume()
{
    G4Material* material = G4Material::GetMaterial(this->can_material);
    if( !material ) {
        G4cout << " ----> Material " << this->can_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.6,0.6,0.6));
    vis_att->SetVisibility(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double z_position;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Place Aluminum Can
    /////////////////////////////////////////////////////////////////////
    G4Tubs* can_cylinder = BuildAluminumCan();

    //logical volume
    if( can_cylinder_log == NULL )
    {
        can_cylinder_log = new G4LogicalVolume(can_cylinder, material, "can_cylinder_log", 0, 0, 0);
        can_cylinder_log->SetVisAttributes(vis_att);
    }

    // place front can_lid
    z_position 	= 0;
    move 		= z_position * direction;

    //add physical cylinder
    this->assembly->AddPlacedVolume(can_cylinder_log, move, rotate);

    /////////////////////////////////////////////////////////////////////
    // Build and Place Aluminum Front Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* can_front_lid = BuildAluminumCanFrontLid();

    // logical volume
    if( can_front_lid_log == NULL )
    {
        can_front_lid_log = new G4LogicalVolume(can_front_lid, material, "can_front_lid_log", 0, 0, 0);
        can_front_lid_log->SetVisAttributes(vis_att);
    }

    // place front can_lid

    z_position 	= -(detector_length_z/2.0) + (can_front_lid_thickness/2.0);
    move 		= z_position * direction;

    //add physical front can_lid
    this->assembly->AddPlacedVolume(can_front_lid_log, move, rotate);
    
    /////////////////////////////////////////////////////////////////////
    // Build and Place Aluminum Back Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* can_back_lid = BuildAluminumCanBackLid();
    
    // logical volume
    if( can_back_lid_log == NULL )
    {
        can_back_lid_log = new G4LogicalVolume(can_back_lid, material, "can_back_lid_log", 0, 0, 0);
        can_back_lid_log->SetVisAttributes(vis_att);
    }

    // place back can_lid
    z_position 	= (detector_length_z/2.0) - (can_back_lid_thickness/2.0);
    move 		= z_position * direction;

    // add physical back can_lid
    this->assembly->AddPlacedVolume(can_back_lid_log, move, rotate);

    return 1;
}

G4int DetectionSystemSodiumIodide::BuildPackingVolume()
{
    G4Material* material = G4Material::GetMaterial(this->packing_material);
    if( !material ) {
        G4cout << " ----> Material " << this->packing_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    vis_att->SetVisibility(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double z_position;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Place Packing
    /////////////////////////////////////////////////////////////////////
    G4Tubs* packing_cylinder = BuildPacking();

    //logical volume
    if( packing_cylinder_log == NULL )
    {
        packing_cylinder_log = new G4LogicalVolume(packing_cylinder, material, "packing_cylinder_log", 0, 0, 0);
        packing_cylinder_log->SetVisAttributes(vis_att);
    }

    // place cylinder
    z_position 	= ( (can_front_lid_thickness + seal_front_lid_thickness + disc_front_lid_thickness + packing_front_lid_thickness) - (can_back_lid_thickness) )/2.0;
    move          = z_position * direction;

    //add physical cylinder
    this->assembly->AddPlacedVolume(packing_cylinder_log, move, rotate);

    /////////////////////////////////////////////////////////////////////
    // Build and Place Front Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* packing_front_lid = BuildPackingFrontLid();

    // logical volume
    if( packing_front_lid_log == NULL )
    {
        packing_front_lid_log = new G4LogicalVolume(packing_front_lid, material, "packing_front_lid_log", 0, 0, 0);
        packing_front_lid_log->SetVisAttributes(vis_att);
    }

    // place front packing_lid
    z_position 	= -(detector_length_z/2.0) + (can_front_lid_thickness + seal_front_lid_thickness + disc_front_lid_thickness + (packing_front_lid_thickness/2.0));
    move 		= z_position * direction;

    //add physical front packing_lid
    this->assembly->AddPlacedVolume(packing_front_lid_log, move, rotate);

    return 1;
}

G4int DetectionSystemSodiumIodide::BuildDiscVolume()
{
    G4Material* material = G4Material::GetMaterial(this->disc_material);
    if( !material ) {
        G4cout << " ----> Material " << this->disc_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    vis_att->SetVisibility(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double z_position;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Place Front Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* disc_front_lid = BuildDiscFrontLid();

    // logical volume
    if( disc_front_lid_log == NULL )
    {
        disc_front_lid_log = new G4LogicalVolume(disc_front_lid, material, "disc_front_lid_log", 0, 0, 0);
        disc_front_lid_log->SetVisAttributes(vis_att);
    }

    // place front disc_lid
    z_position 	= -(detector_length_z/2.0) + (can_front_lid_thickness + seal_front_lid_thickness + (disc_front_lid_thickness/2.0));
    move 		= z_position * direction;

    //add physical front disc_lid
    this->assembly->AddPlacedVolume(disc_front_lid_log, move, rotate);

    return 1;
}

G4int DetectionSystemSodiumIodide::BuildSealVolume()
{
    G4Material* material = G4Material::GetMaterial(this->seal_material);
    if( !material ) {
        G4cout << " ----> Material " << this->seal_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
    vis_att->SetVisibility(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double z_position;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Place Front Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* seal_front_lid = BuildSealFrontLid();

    // logical volume
    if( seal_front_lid_log == NULL )
    {
        seal_front_lid_log = new G4LogicalVolume(seal_front_lid, material, "seal_front_lid_log", 0, 0, 0);
        seal_front_lid_log->SetVisAttributes(vis_att);
    }

    // place front seal_lid
    z_position 	= -(detector_length_z/2.0) + (can_front_lid_thickness + (seal_front_lid_thickness/2.0));
    move 		= z_position * direction;

    //add physical front seal_lid
    this->assembly->AddPlacedVolume(seal_front_lid_log, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Tubs* DetectionSystemSodiumIodide::BuildCrystal()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = crystal_inner_radius;
    G4double outer_radius = crystal_outer_radius;
    G4double half_length_z = (crystal_length_z)/2.0;

    G4Tubs* crystal_block = new G4Tubs("crystal_block", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return crystal_block;
}//end ::BuildCrystal


G4Tubs* DetectionSystemSodiumIodide::BuildAluminumCan()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius 	= can_inner_radius;
    G4double outer_radius 	= can_outer_radius;
    G4double half_length_z 	= can_length_z/2.0;

    G4Tubs* can_cylinder = new G4Tubs("can_cylinder", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return can_cylinder;
}//end ::BuildAluminumCan

G4Tubs* DetectionSystemSodiumIodide::BuildAluminumCanFrontLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = can_lid_inner_radius;
    G4double outer_radius = can_lid_outer_radius;
    G4double half_length_z = can_front_lid_thickness/2.0;

    G4Tubs* can_lid = new G4Tubs("can_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return can_lid;
}//end ::BuildAluminumFrontLid

G4Tubs* DetectionSystemSodiumIodide::BuildAluminumCanBackLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = can_lid_inner_radius;
    G4double outer_radius = can_lid_outer_radius;
    G4double half_length_z = can_back_lid_thickness/2.0;

    G4Tubs* can_lid = new G4Tubs("can_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return can_lid;
}//end ::BuildAluminumBackLid

G4Tubs* DetectionSystemSodiumIodide::BuildPacking()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius 	= packing_inner_radius;
    G4double outer_radius 	= packing_outer_radius;
    G4double half_length_z 	= packing_length_z/2.0;

    G4Tubs* packing_cylinder = new G4Tubs("packing_cylinder", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return packing_cylinder;
}//end ::BuildPacking

G4Tubs* DetectionSystemSodiumIodide::BuildPackingFrontLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = packing_lid_inner_radius;
    G4double outer_radius = packing_lid_outer_radius;
    G4double half_length_z = packing_front_lid_thickness/2.0;

    G4Tubs* packing_lid = new G4Tubs("packing_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return packing_lid;
}//end ::BuildPackingFrontLid

G4Tubs* DetectionSystemSodiumIodide::BuildDiscFrontLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = disc_lid_inner_radius;
    G4double outer_radius = disc_lid_outer_radius;
    G4double half_length_z = disc_front_lid_thickness/2.0;

    G4Tubs* disc_lid = new G4Tubs("disc_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return disc_lid;
}//end ::BuildDiscFrontLid

G4Tubs* DetectionSystemSodiumIodide::BuildSealFrontLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = seal_lid_inner_radius;
    G4double outer_radius = seal_lid_outer_radius;
    G4double half_length_z = seal_front_lid_thickness/2.0;

    G4Tubs* seal_lid = new G4Tubs("seal_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return seal_lid;
}//end ::BuildSealFrontLid

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemSodiumIodide::GetDirectionXYZ(G4double theta, G4double phi)
{
    G4double x,y,z;
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);

    G4ThreeVector direction = G4ThreeVector(x,y,z);

    return direction;
}//end ::GetDirection

