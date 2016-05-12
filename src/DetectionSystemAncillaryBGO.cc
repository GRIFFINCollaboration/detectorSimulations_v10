#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
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

#include "DetectionSystemAncillaryBGO.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


DetectionSystemAncillaryBGO::DetectionSystemAncillaryBGO() :
    // LogicalVolumes
    detector_volume_log(0),
    bgo_block_log(0),
    vacuum_block_log(0),
    can_cylinder_log(0),
    can_cylinder_back_cover_log(0),
    hevimet_block_log(0)
{
    // can //
    can_length              = 190.0*mm;
    can_thickness           = 0.5*mm;
    can_thickness_front     = 3.0*mm;
    can_inner_radius        = 32.5*mm;
    can_material            = "G4_Al";
    // BGO //
    bgo_length              = 106.5*mm;
    bgo_thickness           = 20.0*mm;
    bgo_material            = "G4_BGO";
    // Gap (Vacuum) //
    gap_thickness           = 0.5*mm;
    gap_thickness_outer     = 3.5*mm;
    gap_material            = "Vacuum";

    // distance used in chopping process //
    can_face_to_build_origin= 165.0*mm;
    // Chopping angles //
    chopping_outer_angle    = 12.76*deg;
    chopping_inner_angle    = 8.0*deg;
    // Hevimet //
    hevimet_thickness       = 10.0*mm;
    hevimet_material        = "Hevimetal";

    //G4double triangleThetaAngle = (180/M_PI)*(atan((1/sqrt(3))/sqrt((11/12) + (1/sqrt(2))) )+atan((sqrt(2))/(1+sqrt(2))))*deg;
    G4double triangleThetaAngle = 54.735610317245360*deg;

    // theta
    this->detectorAngles[0][0] 	= triangleThetaAngle;
    this->detectorAngles[1][0] 	= triangleThetaAngle;
    this->detectorAngles[2][0] 	= triangleThetaAngle;
    this->detectorAngles[3][0] 	= triangleThetaAngle;
    this->detectorAngles[4][0] 	= 180.0*deg - triangleThetaAngle;
    this->detectorAngles[5][0] 	= 180.0*deg - triangleThetaAngle;
    this->detectorAngles[6][0] 	= 180.0*deg - triangleThetaAngle;
    this->detectorAngles[7][0] 	= 180.0*deg - triangleThetaAngle;
    // phi
    this->detectorAngles[0][1] 	= 22.5*deg;
    this->detectorAngles[1][1] 	= 112.5*deg;
    this->detectorAngles[2][1] 	= 202.5*deg;
    this->detectorAngles[3][1] 	= 292.5*deg;
    this->detectorAngles[4][1] 	= 22.5*deg;
    this->detectorAngles[5][1] 	= 112.5*deg;
    this->detectorAngles[6][1] 	= 202.5*deg;
    this->detectorAngles[7][1] 	= 292.5*deg;
    // yaw (alpha)
    this->detectorAngles[0][2] 	= 0.0*deg;
    this->detectorAngles[1][2] 	= 0.0*deg;
    this->detectorAngles[2][2] 	= 0.0*deg;
    this->detectorAngles[3][2] 	= 0.0*deg;
    this->detectorAngles[4][2] 	= 0.0*deg;
    this->detectorAngles[5][2] 	= 0.0*deg;
    this->detectorAngles[6][2] 	= 0.0*deg;
    this->detectorAngles[7][2] 	= 0.0*deg;
    // pitch (beta)
    this->detectorAngles[0][3] 	= triangleThetaAngle;
    this->detectorAngles[1][3] 	= triangleThetaAngle;
    this->detectorAngles[2][3] 	= triangleThetaAngle;
    this->detectorAngles[3][3] 	= triangleThetaAngle;
    this->detectorAngles[4][3] 	= 180.0*deg - triangleThetaAngle;
    this->detectorAngles[5][3] 	= 180.0*deg - triangleThetaAngle;
    this->detectorAngles[6][3] 	= 180.0*deg - triangleThetaAngle;
    this->detectorAngles[7][3] 	= 180.0*deg - triangleThetaAngle;
    // roll (gamma)
    this->detectorAngles[0][4] 	= 22.5*deg;
    this->detectorAngles[1][4] 	= 112.5*deg;
    this->detectorAngles[2][4] 	= 202.5*deg;
    this->detectorAngles[3][4] 	= 292.5*deg;
    this->detectorAngles[4][4] 	= 22.5*deg;
    this->detectorAngles[5][4] 	= 112.5*deg;
    this->detectorAngles[6][4] 	= 202.5*deg;
    this->detectorAngles[7][4] 	= 292.5*deg;
}

DetectionSystemAncillaryBGO::~DetectionSystemAncillaryBGO()
{
    // LogicalVolumes
    delete detector_volume_log;
    delete bgo_block_log;
    delete vacuum_block_log;
    delete can_cylinder_log;
    delete can_cylinder_back_cover_log;
    delete hevimet_block_log;
}


G4int DetectionSystemAncillaryBGO::Build()//G4SDManager* mySDman)
{
    // Build assembly volume
    G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    this->assembly = myAssembly;

    G4AssemblyVolume* myAssemblyHevimet = new G4AssemblyVolume();
    this->assemblyHevimet = myAssemblyHevimet;

    G4cout << "BuildBGOVolume" << G4endl;
    BuildBGOVolume();
    G4cout << "BuildVacuumVolume" << G4endl;
    BuildVacuumVolume();
    G4cout << "BuildAluminumCanVolume" << G4endl;
    BuildAluminumCanVolume();
    G4cout << "BuildAluminumCanBackCoverVolume" << G4endl;
    BuildAluminumCanBackCoverVolume();
    G4cout << "BuildHevimetPiece" << G4endl;
    BuildHevimetPiece();

    //  G4cout << "BuildAluminumCanVolume" << G4endl;
    //  BuildAluminumCanVolume();
    //  G4cout << "BuildPackingVolume" << G4endl;
    //  BuildPackingVolume();
    //  G4cout << "BuildDiscVolume" << G4endl;
    //  BuildDiscVolume();
    //  G4cout << "BuildSealVolume" << G4endl;
    //  BuildSealVolume();

    return 1;
}

G4int DetectionSystemAncillaryBGO::PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector_number, G4double radialpos, G4int hevimetopt)
{
    G4int detector_copy_ID = 0;

    G4cout << "AncillaryBGO Detector Number = " << detector_number << G4endl;

    G4int copy_number = detector_copy_ID + detector_number;

    G4double position = radialpos + ( this->can_length / 2.0 ) ;

    G4double theta  = this->detectorAngles[detector_number][0];
    G4double phi    = this->detectorAngles[detector_number][1];
    //G4double alpha  = this->detectorAngles[detector_number][2]; // yaw
    G4double beta   = this->detectorAngles[detector_number][3]; // pitch
    G4double gamma  = this->detectorAngles[detector_number][4]; // roll

    G4double x = 0;
    G4double y = 0;
    G4double z = position;

    G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
    G4ThreeVector move;
    for (G4int i=0; i<3; i++) {
        rotate->set(0,0,0);
        rotate->rotateZ((120*deg)*i);
        if(detector_number > 3) {
            rotate->rotateZ(180*deg);
        }
        rotate->rotateY(M_PI);
        rotate->rotateY(M_PI+beta);
        rotate->rotateZ(gamma);
        move = G4ThreeVector(transX(x,y,z,theta,phi), transY(x,y,z,theta,phi), transZ(x,y,z,theta,phi));
        assembly->MakeImprint(exp_hall_log, move, rotate, copy_number);
    }

    if(hevimetopt == 1) {
        for (G4int i=0; i<3; i++) {
            rotate->set(0,0,0);
            rotate->rotateZ((120*deg)*i);
            if(detector_number > 3) {
                rotate->rotateZ(180*deg);
            }
            rotate->rotateY(M_PI);
            rotate->rotateY(M_PI+beta);
            rotate->rotateZ(gamma);
            move = G4ThreeVector(transX(x,y,z,theta,phi), transY(x,y,z,theta,phi), transZ(x,y,z,theta,phi));
            assemblyHevimet->MakeImprint(exp_hall_log, move, rotate, copy_number);
        }
    }

    return 1;
}

G4int DetectionSystemAncillaryBGO::BuildBGOVolume()
{
    G4Material* material = G4Material::GetMaterial(this->bgo_material);
    if( !material ) {
        G4cout << " ----> Material " << this->bgo_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.2,1.0,0.3));
    vis_att->SetVisibility(true);

    G4SubtractionSolid* bgo_block = BGOPiece();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= can_thickness_front + gap_thickness + (bgo_length - can_length)/2.0;
    G4ThreeVector move 		= z_position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( bgo_block_log == NULL )
    {
        bgo_block_log = new G4LogicalVolume(bgo_block, material, "ancillary_bgo_block_log", 0, 0, 0);
        bgo_block_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(bgo_block_log, move, rotate);

    return 1;
}


G4int DetectionSystemAncillaryBGO::BuildVacuumVolume()
{
    G4Material* material = G4Material::GetMaterial(this->gap_material);
    if( !material ) {
        G4cout << " ----> Material " << this->gap_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    vis_att->SetVisibility(true);

    G4SubtractionSolid* vacuum_block = VacuumPiece();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= can_thickness_front + (bgo_length + gap_thickness - can_length)/2.0;
    G4ThreeVector move 		= z_position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( vacuum_block_log == NULL )
    {
        vacuum_block_log = new G4LogicalVolume(vacuum_block, material, "ancillary_vacuum_block_log", 0, 0, 0);
        vacuum_block_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(vacuum_block_log, move, rotate);

    return 1;
}

G4int DetectionSystemAncillaryBGO::BuildAluminumCanVolume()
{
    G4Material* material = G4Material::GetMaterial(this->can_material);
    if( !material ) {
        G4cout << " ----> Material " << this->can_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    vis_att->SetVisibility(true);

    G4SubtractionSolid* can_cylinder = AluminumCanPiece();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= 0.0*mm;
    G4ThreeVector move 		= z_position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( can_cylinder_log == NULL )
    {
        can_cylinder_log = new G4LogicalVolume(can_cylinder, material, "ancillary_can_cylinder_log", 0, 0, 0);
        can_cylinder_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(can_cylinder_log, move, rotate);

    return 1;
}

G4int DetectionSystemAncillaryBGO::BuildAluminumCanBackCoverVolume()
{
    G4Material* material = G4Material::GetMaterial(this->can_material);
    if( !material ) {
        G4cout << " ----> Material " << this->can_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    vis_att->SetVisibility(true);

    G4Tubs* can_cylinder = AluminumCanBackCover();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= (can_length - can_thickness)/2.0;
    G4ThreeVector move 		= z_position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( can_cylinder_back_cover_log == NULL )
    {
        can_cylinder_back_cover_log = new G4LogicalVolume(can_cylinder, material, "ancillary_can_cylinder_back_cover_log", 0, 0, 0);
        can_cylinder_back_cover_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(can_cylinder_back_cover_log, move, rotate);

    return 1;
}

G4int DetectionSystemAncillaryBGO::BuildHevimetPiece()
{
    G4Material* material = G4Material::GetMaterial(this->hevimet_material);
    if( !material ) {
        G4cout << " ----> Material " << this->can_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    vis_att->SetVisibility(true);

    G4SubtractionSolid* can_cylinder = HevimetPiece();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= -1.0*(hevimet_thickness + can_length)/2.0;
    G4ThreeVector move 		= z_position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( hevimet_block_log == NULL )
    {
        hevimet_block_log = new G4LogicalVolume(can_cylinder, material, "ancillary_hevimet_block_log", 0, 0, 0);
        hevimet_block_log->SetVisAttributes(vis_att);
    }

    this->assemblyHevimet->AddPlacedVolume(hevimet_block_log, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4SubtractionSolid* DetectionSystemAncillaryBGO::BGOPiece()
{
    G4double start_phi = 0.0*deg;
    G4double end_phi = 120.0*deg;

    G4double inner_radius = can_inner_radius+can_thickness+gap_thickness;
    G4double outer_radius = inner_radius+bgo_thickness;
    G4double half_length_z = bgo_length/2.0;

    G4Tubs* bgo_block = new G4Tubs("bgo_block", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // turn the block to be symmetric about 0 degress so we can cut it.
    bgo_block->SetStartPhiAngle(-60.0*deg,true);
    bgo_block->SetDeltaPhiAngle(120.0*deg);

    // grab the scissors, time to cut.
    // make the cut surface unnecessarily large.
    G4double cut_half_length_z = 2.0*can_face_to_build_origin;
    G4double cut_half_length_x = 2.0*can_face_to_build_origin;
    G4double cut_half_length_y = 2.0*can_face_to_build_origin;

    G4Box* cut_plate = new G4Box("cut_plate", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    // move cutting box to the "build origin" (see technical drawings)
    G4double move_cut_z = -1.0*(can_face_to_build_origin + gap_thickness + can_thickness_front + bgo_length/2.0);
    G4ThreeVector move_cut((cut_half_length_x-gap_thickness-can_thickness)/(cos(chopping_outer_angle)),0,move_cut_z);

    // rotate our cutting block to the "chopping" angle
    G4RotationMatrix* rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateY(-1.0*chopping_outer_angle);

    // snip snip
    G4SubtractionSolid* bgo_block_cut = new G4SubtractionSolid("bgo_block_cut", bgo_block, cut_plate, rotate_cut, move_cut);

    return bgo_block_cut;
}

G4SubtractionSolid* DetectionSystemAncillaryBGO::VacuumPiece()
{
    G4double start_phi = 0.0*deg;
    G4double end_phi = 120.0*deg;

    G4double inner_radius = can_inner_radius+can_thickness;
    G4double outer_radius = inner_radius+gap_thickness+bgo_thickness+gap_thickness_outer;
    G4double half_length_z = (bgo_length + gap_thickness)/2.0;

    G4Tubs* vacuum_block = new G4Tubs("vacuum_block", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // turn the block to be symmetric about 0 degress so we can cut it.
    vacuum_block->SetStartPhiAngle(-60.0*deg,true);
    vacuum_block->SetDeltaPhiAngle(120.0*deg);

    // grab the scissors, time to cut.
    // make the cut surface unnecessarily large.
    G4double cut_half_length_z = 2.0*can_face_to_build_origin;
    G4double cut_half_length_x = 2.0*can_face_to_build_origin;
    G4double cut_half_length_y = 2.0*can_face_to_build_origin;

    G4Box* cut_plate = new G4Box("cut_plate", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    // move cutting box to the "build origin" (see technical drawings)
    G4double move_cut_z = -1.0*(can_face_to_build_origin + gap_thickness + can_thickness_front + bgo_length/2.0);
    G4ThreeVector move_cut((cut_half_length_x-can_thickness)/(cos(chopping_outer_angle)),0,move_cut_z);

    // rotate our cutting block to the "chopping" angle
    G4RotationMatrix* rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateY(-1.0*chopping_outer_angle);

    // snip snip
    G4SubtractionSolid* vacuum_block_cut = new G4SubtractionSolid("vacuum_block_cut", vacuum_block, cut_plate, rotate_cut, move_cut);

    // Now we have to do another cut on the vacuum.
    // We need to remove the area where the BGO sits
    // let's take this piece larger in phi to make sure we cut everything
    start_phi = 0.0*deg;
    end_phi = 140.0*deg;

    inner_radius = can_inner_radius+can_thickness+gap_thickness;
    outer_radius = inner_radius+bgo_thickness;
    half_length_z = (bgo_length + gap_thickness)/2.0;

    G4Tubs* vacuum_block_cut_bgo = new G4Tubs("vacuum_block_cut_bgo", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // turn the block to be symmetric about 0 degress so we can cut it.
    vacuum_block_cut_bgo->SetStartPhiAngle(-70.0*deg,true);
    vacuum_block_cut_bgo->SetDeltaPhiAngle(140.0*deg);

    // grab the scissors, time to cut.
    // make the cut surface unnecessarily large.
    cut_half_length_z = 2.0*can_face_to_build_origin;
    cut_half_length_x = 2.0*can_face_to_build_origin;
    cut_half_length_y = 2.0*can_face_to_build_origin;

    G4Box* cut_plate_for_bgo = new G4Box("cut_plate_for_bgo", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    // move cutting box to the "build origin" (see technical drawings)
    move_cut_z = -1.0*(can_face_to_build_origin + gap_thickness + can_thickness_front + (bgo_length + gap_thickness)/2.0);
    G4ThreeVector move_cut_for_bgo((cut_half_length_x-can_thickness-gap_thickness)/(cos(chopping_outer_angle)),0,move_cut_z);

    // rotate our cutting block to the "chopping" angle
    G4RotationMatrix* rotate_cut_for_bgo = new G4RotationMatrix;
    rotate_cut_for_bgo->rotateY(-1.0*chopping_outer_angle);

    // snip snip
    G4SubtractionSolid* bgo_block_to_be_cut_out = new G4SubtractionSolid("bgo_block_to_be_cut_out", vacuum_block_cut_bgo, cut_plate_for_bgo, rotate_cut_for_bgo, move_cut_for_bgo);

    // now we can remove the BGO area from inside the vacuum
    // move cut just along z-axis by the amount of the gap_thickness (both vacuum and bgo blacks are same size)
    G4ThreeVector move_cut_me(0,0,gap_thickness);
    // no rotation
    G4RotationMatrix* rotate_cut_me = new G4RotationMatrix;
    // final snip snip
    G4SubtractionSolid* vacuum_block_with_bgo_cut_out = new G4SubtractionSolid("vacuum_block_with_bgo_cut_out", vacuum_block_cut, bgo_block_to_be_cut_out, rotate_cut_me, move_cut_me);

    return vacuum_block_with_bgo_cut_out;
}

G4SubtractionSolid* DetectionSystemAncillaryBGO::AluminumCanPiece()
{
    G4double start_phi = 0.0*deg;
    G4double end_phi = 120.0*deg;

    G4double inner_radius = can_inner_radius;
    G4double outer_radius = inner_radius+can_thickness+gap_thickness+bgo_thickness+gap_thickness_outer+can_thickness;
    G4double half_length_z = (can_length)/2.0;

    G4Tubs* can_block = new G4Tubs("can_block", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // turn the block to be symmetric about 0 degress so we can cut it.
    can_block->SetStartPhiAngle(-60.0*deg,true);
    can_block->SetDeltaPhiAngle(120.0*deg);

    // grab the scissors, time to cut.
    // make the cut surface unnecessarily large.
    G4double cut_half_length_z = 2.0*can_face_to_build_origin;
    G4double cut_half_length_x = 2.0*can_face_to_build_origin;
    G4double cut_half_length_y = 2.0*can_face_to_build_origin;

    G4Box* cut_plate = new G4Box("cut_plate", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    // move cutting box to the "build origin" (see technical drawings)
    G4double move_cut_z = -1.0*(can_face_to_build_origin + can_length/2.0);
    G4ThreeVector move_cut((cut_half_length_x)/(cos(chopping_outer_angle)),0,move_cut_z);

    // rotate our cutting block to the "chopping" angle
    G4RotationMatrix* rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateY(-1.0*chopping_outer_angle);

    // snip snip
    G4SubtractionSolid* can_block_cut = new G4SubtractionSolid("can_block_cut", can_block, cut_plate, rotate_cut, move_cut);

    // Now we have to do another cut on the can.
    // We need to remove the area where the vacuum and BGO sits
    // let's take this piece larger in phi to make sure we cut everything
    start_phi = 0.0*deg;
    end_phi = 140.0*deg;

    inner_radius = can_inner_radius+can_thickness;
    outer_radius = inner_radius+gap_thickness+bgo_thickness+gap_thickness_outer;
    half_length_z = (can_length)/2.0;

    G4Tubs* can_block_cut_vacuum = new G4Tubs("can_block_cut_vacuum", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // turn the block to be symmetric about 0 degress so we can cut it.
    can_block_cut_vacuum->SetStartPhiAngle(-70.0*deg,true);
    can_block_cut_vacuum->SetDeltaPhiAngle(140.0*deg);

    // grab the scissors, time to cut.
    // make the cut surface unnecessarily large.
    cut_half_length_z = 2.0*can_face_to_build_origin;
    cut_half_length_x = 2.0*can_face_to_build_origin;
    cut_half_length_y = 2.0*can_face_to_build_origin;

    G4Box* cut_plate3 = new G4Box("cut_plate3", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    // move cutting box to the "build origin" (see technical drawings)
    move_cut_z = -1.0*(can_face_to_build_origin + can_thickness_front + can_length/2.0);
    G4ThreeVector move_cut3((cut_half_length_x-can_thickness)/(cos(chopping_outer_angle)),0,move_cut_z);

    // rotate our cutting block to the "chopping" angle
    G4RotationMatrix* rotate_cut3 = new G4RotationMatrix;
    rotate_cut3->rotateY(-1.0*chopping_outer_angle);

    // snip snip
    G4SubtractionSolid* vacuum_block_to_be_cut_out = new G4SubtractionSolid("vacuum_block_to_be_cut_out", can_block_cut_vacuum, cut_plate3, rotate_cut3, move_cut3);

    // now we can remove the vacuum/BGO area from inside the can
    // move cut just along z-axis by the amount of the can_thickness_front (both vacuum and can blacks are same size)
    G4ThreeVector move_cut_me(0,0,can_thickness_front);
    // no rotation
    G4RotationMatrix* rotate_cut_me = new G4RotationMatrix;
    // final snip snip
    G4SubtractionSolid* can_block_with_vacuum_cut_out = new G4SubtractionSolid("can_block_cut", can_block_cut, vacuum_block_to_be_cut_out, rotate_cut_me, move_cut_me);

    return can_block_with_vacuum_cut_out;
}


G4Tubs* DetectionSystemAncillaryBGO::AluminumCanBackCover()
{
    G4double start_phi = 0.0*deg;
    G4double end_phi = 120.0*deg;

    G4double inner_radius = can_inner_radius+can_thickness;
    G4double outer_radius = inner_radius+can_thickness+gap_thickness+bgo_thickness+gap_thickness_outer;
    G4double half_length_z = (can_thickness)/2.0;

    G4Tubs* can_block = new G4Tubs("can_block", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // turn the block to be symmetric about 0 degress
    can_block->SetStartPhiAngle(-60.0*deg,true);
    can_block->SetDeltaPhiAngle(120.0*deg);

    return can_block;
}

G4SubtractionSolid* DetectionSystemAncillaryBGO::HevimetPiece() // this has been fix
{
    G4double start_phi = 0.0*deg;
    G4double end_phi = 120.0*deg;

    G4double inner_radius = 0.0*mm;
    G4double outer_radius = can_inner_radius+can_thickness+gap_thickness+bgo_thickness+gap_thickness_outer+can_thickness;
    G4double half_length_z = (hevimet_thickness)/2.0;

    G4Tubs* hevimet_block = new G4Tubs("hevimet_block", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // turn the block to be symmetric about 0 degress so we can cut it.
    hevimet_block->SetStartPhiAngle(-60.0*deg,true);
    hevimet_block->SetDeltaPhiAngle(120.0*deg);

    // grab the scissors, time to cut.
    // make the cut surface unnecessarily large.
    G4double cut_half_length_z = 2.0*can_face_to_build_origin;
    G4double cut_half_length_x = 2.0*can_face_to_build_origin;
    G4double cut_half_length_y = 2.0*can_face_to_build_origin;

    G4Box* cut_plate = new G4Box("cut_plate", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    // move cutting box to the "build origin" (see technical drawings)
    G4double move_cut_z = -1.0*(can_face_to_build_origin)+hevimet_thickness/2.0;
    G4ThreeVector move_cut((cut_half_length_x)/(cos(chopping_outer_angle)),0,move_cut_z);

    // rotate our cutting block to the "chopping" angle
    G4RotationMatrix* rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateY(-1.0*chopping_outer_angle);

    // snip snip
    G4SubtractionSolid* hevimet_block_cut = new G4SubtractionSolid("hevimet_block_cut", hevimet_block, cut_plate, rotate_cut, move_cut);


    // Now we need to cone in the middle.
    G4double cone_inner_radius = 0.0*mm;
    G4double cone_lower_outer_radius = 0.0*mm;
    G4double cone_upper_outer_radius = tan(chopping_inner_angle)*(can_face_to_build_origin+hevimet_thickness);
    G4double cone_half_length_z = (can_face_to_build_origin+hevimet_thickness)/2.0;

    G4Cons* hevimet_cut_cone = new G4Cons("hevimet_cut_cone", cone_inner_radius, cone_lower_outer_radius, cone_inner_radius, cone_upper_outer_radius, cone_half_length_z, 0.0*deg, 360.0*deg);

    // now we can remove the inner radius from the Hevimet
    // move cut just along z-axis by the amount of the can_thickness_front (both vacuum and can blacks are same size)
    G4ThreeVector move_cut_me(0,0,-1.0*(can_face_to_build_origin)/2.0+hevimet_thickness/2.0);
    // no rotation
    G4RotationMatrix* rotate_cut_me = new G4RotationMatrix;
    // final snip snip
    G4SubtractionSolid* hevimet_block_cut_with_cone_cut = new G4SubtractionSolid("hevimet_block_cut_with_cone_cut", hevimet_block_cut, hevimet_cut_cone, rotate_cut_me, move_cut_me);


    return hevimet_block_cut_with_cone_cut;
}


//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemAncillaryBGO::GetDirectionXYZ(G4double theta, G4double phi)
{
    G4double x,y,z;
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);

    G4ThreeVector direction = G4ThreeVector(x,y,z);

    return direction;
}//end ::GetDirection

G4double DetectionSystemAncillaryBGO::transX(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*cos(phi) );
}

G4double DetectionSystemAncillaryBGO::transY(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*sin(phi) );
}

G4double DetectionSystemAncillaryBGO::transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*cos(theta) );
}
