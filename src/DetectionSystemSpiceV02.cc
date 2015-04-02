#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemSpiceV02.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

DetectionSystemSpiceV02::DetectionSystemSpiceV02() :
    // LogicalVolumes
    detector_casing_side_log(0),
    detector_casing_back_log(0),
    crystal_block_log(0)
{
    /////////////////////////////////////////////////////////////////////
    // SPICE Physical Properties
    /////////////////////////////////////////////////////////////////////

    this->casing_material            = "Aluminum";
    this->wafer_material             = "Silicon";

    //-----------------------------//
    // parameters for the square   //
    // planar detector crystal     //
    //-----------------------------//
    this->squareDetCrystalLength = 38*mm;
    this->squareDet90DegCrystalRadius = 63.64*mm;
    this->squareDet45DegCrystalRadius = 45*mm;
    this->squareDetCrystalThickness = 5*mm;

    //-----------------------------//
    // parameters for the square   //
    // planar detector casing      //
    //-----------------------------//
    this->squareDetCasingLength = 42*mm;
    this->squareDetCasingThickness = 14*mm;

    // minimial distance 122.2mm, assuming 2mm gap between casings
    // calculated by distance of back plate - detector casing thickness (z)
    // - 150.65 (back face)
    // + 20 (2 times mount cylinders)
    // + 14 (casing size z)
    this->detector_face_target_distance              = -116.65*mm;

    //this->my_data_output = data_output;

    // These keep track of the number of sensitive detectors. Initialized to zero
    //this->siDetCopyNumber=SI_DET_COPY_NUMBER;
    //this->siDetCasingSideCopyNumber=SI_DET_CASING_SIDE_COPY_NUMBER;
    //this->siDetCasingBackCopyNumber=SI_DET_CASING_BACK_COPY_NUMBER;

#define NUMBER_OF_CELLS 4

}

DetectionSystemSpiceV02::~DetectionSystemSpiceV02()
{
    // LogicalVolumes in ConstructSPICEDetectionSystem
    delete crystal_block_log;
    delete detector_casing_side_log;
    delete detector_casing_back_log;

}

G4int DetectionSystemSpiceV02::Build() 
{

    // Build assembly volume
    G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    this->assembly = myAssembly;
    G4AssemblyVolume* myAssemblySi = new G4AssemblyVolume();
    this->assemblySi = myAssemblySi;

    G4cout << "BuildSiliconWafer" << G4endl;
    BuildSiliconWafer();
    G4cout << "BuildAluminiumCasingSide" << G4endl;
    BuildAluminiumCasingSide();
    G4cout << "BuildAluminiumCasingBack" << G4endl;
    BuildAluminiumCasingBack();

    return 1;
}

G4int DetectionSystemSpiceV02::PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate, G4int detector_number)
{
    G4int detector_copy_ID = 0;

    G4cout << "SpiceV02 Detector Number = " << detector_number << G4endl;

    G4int copy_number = detector_copy_ID + detector_number;

    assemblySi->MakeImprint(exp_hall_log, move, rotate, copy_number);
    assembly->MakeImprint(exp_hall_log, move, rotate, copy_number);

    return 1;
}

G4int DetectionSystemSpiceV02::BuildSiliconWafer()
{
    G4Material* material = G4Material::GetMaterial(this->wafer_material);
    if( !material ) {
        G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    vis_att->SetVisibility(true);

    G4Box* crystal_block = BuildCrystal();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= 0.0*mm;
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( crystal_block_log == NULL )
    {
        crystal_block_log = new G4LogicalVolume(crystal_block, material, "crystal_block_log", 0, 0, 0);
        crystal_block_log->SetVisAttributes(vis_att);
    }

    this->assemblySi->AddPlacedVolume(crystal_block_log, move, rotate);

    return 1;
}

G4int DetectionSystemSpiceV02::BuildAluminiumCasingSide()
{
    G4Material* material = G4Material::GetMaterial(this->casing_material);
    if( !material ) {
        G4cout << " ----> Material " << this->casing_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
    vis_att->SetVisibility(true);

    G4SubtractionSolid* detector_casing_side = BuildCasingSide();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= 0.0*mm;
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( detector_casing_side_log == NULL )
    {
        detector_casing_side_log = new G4LogicalVolume(detector_casing_side, material, "detector_casing_log", 0,0,0);
        detector_casing_side_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(detector_casing_side_log, move, rotate);

    return 1;
}

G4int DetectionSystemSpiceV02::BuildAluminiumCasingBack()
{
    G4Material* material = G4Material::GetMaterial(this->casing_material);
    if( !material ) {
        G4cout << " ----> Material " << this->casing_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
    vis_att->SetVisibility(true);

    G4Box* detector_casing_back = BuildCasingBack();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= 0.0*mm;
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( detector_casing_back_log == NULL )
    {
        detector_casing_back_log = new G4LogicalVolume(detector_casing_back, material,"detector_casing_back_log", 0,0,0);
        detector_casing_back_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(detector_casing_back_log, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemSpiceV02::BuildCrystal()
{
    // silicon detector dimensions
    G4double half_length_x = (this->squareDetCrystalLength)/2.0;
    G4double half_length_y = (this->squareDetCrystalLength)/2.0;
    G4double half_length_z = (this->squareDetCrystalThickness)/2.0;

    // establish solid
    G4Box* crystal_block = new G4Box("crystal_block", half_length_x,half_length_y,half_length_z);

    return crystal_block;
}//end ::BuildCrystal

G4SubtractionSolid* DetectionSystemSpiceV02::BuildCasingSide()
{
    // dimensions of casing
    G4double half_length_x = (this->squareDetCasingLength)/2.0;
    G4double half_length_y = (this->squareDetCasingLength)/2.0;
    G4double half_length_z = (this->squareDetCrystalThickness)/2.0;
    // dimensions of silicon detector
    G4double crys_half_length_z = (this->squareDetCrystalThickness)/2.0;

    // relative distance
    G4double z_transl = (half_length_z - crys_half_length_z);
    G4ThreeVector transl(0,0,z_transl);

    // establish solid
    G4Box* filled_box = new G4Box("filled_box",half_length_x,half_length_y,half_length_z);
    G4Box* crystal_block = BuildCrystal();
    G4SubtractionSolid* detector_casing_side = new G4SubtractionSolid("detector_casing_side",filled_box,crystal_block,0,transl);

    return detector_casing_side;
}//end ::BuildCasingSide

G4Box* DetectionSystemSpiceV02::BuildCasingBack()
{
    // dimensions of casing
    G4double half_length_x = (this->squareDetCasingLength)/2.0;
    G4double half_length_y = (this->squareDetCasingLength)/2.0;
    G4double half_length_z = (this->squareDetCasingThickness-this->squareDetCrystalThickness)/2.0;

    // establish solid
    G4Box* detector_casing_back = new G4Box("detector_casing_back", half_length_x,half_length_y,half_length_z);

    return detector_casing_back;
}//end ::BuildCasingBack

