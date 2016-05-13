#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polyhedra.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystem8pi.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

//constructor suppressed

DetectionSystem8pi::~DetectionSystem8pi()
{
    delete germanium_block_log;
    delete germanium_dead_layer_log;
    delete germanium_vacuum_core_log;
    delete lower_electrodeMat_electrode_log;
    delete upper_electrodeMat_electrode_log;
    delete inner_cage_1_log;
    delete inner_cage_2_log;
    delete inner_cage_3_log;
    delete inner_cage_4_log;
    delete inner_cage_bottom_log;
    delete inner_cage_lid_log;
    delete structureMat_cooling_rod_log;
    delete electrodeMat_cooling_rod_log;
    delete cooling_rod_cover_log;
    delete cooling_rod_cover_lid_log;
    delete outer_can_side_log;
    delete outer_can_lid_log;
    delete outer_can_bottom_log;
    delete beryllium_window_log;
    delete inner_BGO_annulus_log;
    delete structureMat_sheath_log;
    delete outer_lower_BGO_annulus_log;
    delete outer_upper_BGO_annulus_log;
    delete liquid_N2_log;
    delete liquid_N2_side_log;
    delete liquid_N2_lid_log;
    delete liquid_N2_bottom_log;
    delete hevimetal_log;
    delete auxMat_plug_log;
    delete auxMat_layer_log;

}

G4int DetectionSystem8pi::Build()//G4SDManager* mySDman)
{

    // Build assembly volumes
    G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    this->assembly = myAssembly;
    G4AssemblyVolume* myAssemblyGe = new G4AssemblyVolume();
    this->assemblyGe = myAssemblyGe;
    G4AssemblyVolume* myAssemblyInnerBGO = new G4AssemblyVolume();
    this->assemblyInnerBGO = myAssemblyInnerBGO;
    G4AssemblyVolume* myAssemblyOuterLowerBGO = new G4AssemblyVolume();
    this->assemblyOuterLowerBGO = myAssemblyOuterLowerBGO;
    G4AssemblyVolume* myAssemblyOuterUpperBGO = new G4AssemblyVolume();
    this->assemblyOuterUpperBGO = myAssemblyOuterUpperBGO;

    AddGermanium();

    AddGermaniumDeadLayer();

    AddGermaniumCore();
    //Add electrodeMat electrode inside germanium core

    AddElectrodeMatElectrode();
    //Add the inner structureMat cage around Ge crystal

    AddStructureMatCage();
    //Add inner structureMat can lids

    AddInnerStructureMatLids();
    //Add beryllium window

    AddStructureMatCoolingRod();
    //Add electrodeMat cooling rod section

    AddElectrodeMatCoolingRod();
    //Add the outer can

    AddBerylliumWindow();
    //Add structureMat cooling rod section

    AddOuterStructureMatCan();
    //add structureMat cover around cooling rod

    AddCoolingRodCover();
    //add structureMat sheath around inner BGO annulus

    AddStructureMatBGOSheath();
    //add inner BGO annulus around cooling rod

    AddInnerBGOAnnulus();
    //add outer BGO annulus

    AddOuterBGOAnnulus();
    //add liquid N2 cooling container

    AddLiquidN2Container();
    //add hevimetal collimator & core

    AddHevimetalCollimator();
    //Add AuxMat plug in Hevimetal collimator

    AddAuxMatPlug();
    //Add thin auxMat layer

    AddThinAuxMatLayer();
    //Add hevimetal and auxMat plug

    return 1;
}//end ::Build 

G4int DetectionSystem8pi::PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate, G4int detector_number)
{
    //G4int detector_copy_ID = 0;

    //G4int copy_number = detector_copy_ID + detector_number;

    assemblyGe->MakeImprint(exp_hall_log, move, rotate, detector_number);
    assemblyInnerBGO->MakeImprint(exp_hall_log, move, rotate, detector_number+(20*1));
    assemblyOuterLowerBGO->MakeImprint(exp_hall_log, move, rotate, detector_number+(20*2));
    assemblyOuterUpperBGO->MakeImprint(exp_hall_log, move, rotate, detector_number+(20*3));
    assembly->MakeImprint(exp_hall_log, move, rotate, detector_number+(20*4));

    return 1;
}

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystem8pi::GetDirectionXYZ(G4double theta, G4double phi)
{
    G4double x,y,z;
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);
    G4ThreeVector direction = G4ThreeVector(x,y,z);
    return direction;
}//end ::end GetDirectionXYZ


//Add the germanium crystal
G4int DetectionSystem8pi::AddGermanium()
{
    //material
    G4Material* material = G4Material::GetMaterial("Germanium");
    if( !material ) {
        G4cout << " ----> Material " << "Germanium" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    //vis attributes
    G4VisAttributes* germanium_block_vis_att = new G4VisAttributes(G4Colour(0.0,0.5,0.0));
    germanium_block_vis_att->SetVisibility(true);

    // germanium detector measurements
    G4double inner_radius = 0.0*cm;
    G4double outer_radius = crystal_outer_radius;
    G4double half_length_z = crystal_length/2.0;

    //primtive volume
    G4Tubs* germanium_block = new G4Tubs("germanium_block", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    // Cut Out Dead Layer! ///////////////////////////////
    //measurements
    G4double cut_inner_radius = 0.0*cm;
    G4double cut_outer_radius =  crystal_inner_radius + dead_layer_thickness + this->cut_clearance;
    G4double cut_half_length_z =  crystal_length/2.0 - hole_starting_depth/2.0 + dead_layer_thickness/2.0 + this->cut_clearance/2.0;
    //primitive volume
    G4Tubs* cut_germanium_dead_layer = new G4Tubs("cut_germanium_dead_layer", cut_inner_radius, cut_outer_radius, cut_half_length_z, 0.0*deg, detail_view_end_angle);

    //positioning
    G4double cut_z_position = hole_starting_depth/2.0 - dead_layer_thickness/2.0;
    G4ThreeVector move_cut = G4ThreeVector(0, 0, cut_z_position);

    G4SubtractionSolid* germanium_block_with_cavity = new G4SubtractionSolid("germanium_block_with_cavity", germanium_block, cut_germanium_dead_layer, 0, move_cut);
    ////////////////////////////////////////////////////////

    //positioning
    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( germanium_block_log == NULL )
    {
        germanium_block_log = new G4LogicalVolume(germanium_block_with_cavity, material, "8pi_germanium_block_log", 0, 0, 0);
        germanium_block_log->SetVisAttributes(germanium_block_vis_att);
    }

    //physical volume
    //  germanium_block_phys = new G4PVPlacement(rotate, move,
    //       "germanium_block_phys", germanium_block_log, exp_hall_phys,
    //       false, detector_number);

    this->assemblyGe->AddPlacedVolume(germanium_block_log, move, rotate);

    return 1;
}//end ::end AddGermanium

//Add the germanium crystal dead layer
G4int DetectionSystem8pi::AddGermaniumDeadLayer()
{
    //material
    G4Material* material = G4Material::GetMaterial("Germanium");
    if( !material ) {
        G4cout << " ----> Material " << "Germanium" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    //vis attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.5,0.0));
    vis_att->SetVisibility(true);

    //measurements
    G4double inner_radius = 0.0*cm;
    G4double outer_radius =  crystal_inner_radius + dead_layer_thickness;
    G4double half_length_z =  crystal_length/2.0 - hole_starting_depth/2.0 + dead_layer_thickness/2.0;
    //primitive volume
    G4Tubs* germanium_dead_layer = new G4Tubs("germanium_dead_layer", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    // Cut Out Vacuum Layer! ///////////////////////////////
    //measurements
    G4double cut_inner_radius = 0.0*cm;
    G4double cut_outer_radius = crystal_inner_radius + this->cut_clearance;
    G4double cut_half_length_z = ((crystal_length - hole_starting_depth) + this->cut_clearance)/2.0;

    //primitive volume
    G4Tubs* cut_germanium_vacuum_core = new G4Tubs("cut_germanium_vacuum_core", cut_inner_radius, cut_outer_radius, cut_half_length_z, 0.0*deg, detail_view_end_angle);

    //positioning
    G4double cut_z_position = dead_layer_thickness/2.0;
    G4ThreeVector move_cut = G4ThreeVector(0, 0, cut_z_position);

    G4SubtractionSolid* germanium_dead_layer_with_cavity = new G4SubtractionSolid("germanium_dead_layer_with_cavity", germanium_dead_layer, cut_germanium_vacuum_core, 0, move_cut);
    ////////////////////////////////////////////////////////

    //positioning
    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_length/2.0 - half_length_z + crystal_dist_from_origin; // adjust for placement
    G4ThreeVector move = G4ThreeVector(0, 0, z_position);

    //logical volume
    if( germanium_dead_layer_log == NULL)
    {
        germanium_dead_layer_log = new G4LogicalVolume(germanium_dead_layer_with_cavity, material, "germanium_dead_layer_log", 0, 0, 0);
        germanium_dead_layer_log->SetVisAttributes(vis_att);
    }

    //physical volume
    //  germanium_dead_layer_phys = new G4PVPlacement(0, move,
    //       germanium_dead_layer_log, "germanium_dead_layer_phys",  germanium_block_log,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(germanium_dead_layer_log, move, rotate);

    return 1;
}//end ::end AddGermaniumDeadLayer


//Add the logical germanium crystal vacuum core
G4int DetectionSystem8pi::AddGermaniumCore()
{
    //material
    G4Material* material = G4Material::GetMaterial("Vacuum");
    if( !material ) {
        G4cout << " ----> Material " << "Vacuum" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    //vis attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    vis_att->SetVisibility(true);

    //measurements
    G4double inner_radius = 0.0*cm;
    G4double outer_radius = crystal_inner_radius;
    G4double half_length_z = (crystal_length - hole_starting_depth)/2.0;

    //primitive volume
    G4Tubs* germanium_vacuum_core = new G4Tubs("germanium_vacuum_core", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    // Cut Out ElectrodeMat Electrode! ///////////////////////////////
    //measurements
    G4double cut_inner_radius = 0.0*cm;
    G4double cut_outer_radius = electrode_radius + this->cut_clearance/2.0;
    G4double cut_half_length_z =  crystal_length/2.0 - hole_starting_depth/2.0 + this->cut_clearance/2.0;
    //primitive volume
    G4Tubs* cut_germanium_vacuum_core = new G4Tubs("cut_germanium_vacuum_core", cut_inner_radius, cut_outer_radius, cut_half_length_z, 0.0*deg, detail_view_end_angle);

    //positioning
    G4double cut_z_position = 0.0;
    G4ThreeVector move_cut = G4ThreeVector(0, 0, cut_z_position);

    G4SubtractionSolid* germanium_vacuum_core_with_cavity = new G4SubtractionSolid("germanium_vacuum_core_with_cavity", germanium_vacuum_core, cut_germanium_vacuum_core, 0, move_cut);
    ////////////////////////////////////////////////////////

    //positioning
    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    //G4double z_position = dead_layer_thickness/2.0 + crystal_length/2.0 - half_length_z + crystal_dist_from_origin; // adjust for placement
    G4double z_position = crystal_length/2.0 - half_length_z + crystal_dist_from_origin; // adjust for placement
    G4ThreeVector move = G4ThreeVector(0, 0, z_position);

    //logical volume
    if( germanium_vacuum_core_log == NULL)
    {
        germanium_vacuum_core_log = new G4LogicalVolume(germanium_vacuum_core_with_cavity, material, "germanium_vacuum_core_log", 0, 0, 0);
        germanium_vacuum_core_log->SetVisAttributes(vis_att);
    }

    //physical volume
    //  germanium_vacuum_core_phys = new G4PVPlacement(0, move,
    //       germanium_vacuum_core_log, "germanium_vacuum_core_phys", germanium_dead_layer_log,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(germanium_vacuum_core_log, move, rotate);

    return 1;
}//end ::end AddGermaniumCore


//Add the electrodeMat electrode inside the germanium crystal core
G4int DetectionSystem8pi::AddElectrodeMatElectrode()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->electrodeMat);
    if( !material ) {
        G4cout << " ----> Material " << "ElectrodeMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    //vis attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.75, 0.45, 0.2));
    vis_att->SetVisibility(true);

    //lower portion of electrodeMat electrode
    //measurements
    G4double inner_radius = 0.0*cm;
    G4double outer_radius = electrode_radius;
    G4double half_length_z =  crystal_length/2.0
            - hole_starting_depth/2.0;
    //primitive volume
    G4Tubs* lower_electrodeMat_cylinder = new G4Tubs("lower_electrodeMat_cylinder", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    //positioning
    G4double lower_z_position = crystal_dist_from_origin + hole_starting_depth/2.0; // check this! // adjust position

    //logical volume
    if( lower_electrodeMat_electrode_log == NULL)
    {
        lower_electrodeMat_electrode_log = new G4LogicalVolume(lower_electrodeMat_cylinder, material, "lower_electrodeMat_electrode_log", 0, 0, 0);
        lower_electrodeMat_electrode_log->SetVisAttributes(vis_att);
    }

    //upper portion of electrodeMat electrode
    //measurements
    half_length_z =  inner_can_extends_past_crystal/2.0;

    //primitive volume
    G4Tubs* upper_electrodeMat_cylinder = new G4Tubs("upper_electrodeMat_cylinder", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    //positioning
    G4double z_position = crystal_dist_from_origin + crystal_length/2.0 + half_length_z;

    //logical volume
    if( upper_electrodeMat_electrode_log == NULL )
    {
        upper_electrodeMat_electrode_log = new G4LogicalVolume(upper_electrodeMat_cylinder, material, "upper_electrodeMat_electrode_log", 0, 0, 0);
        upper_electrodeMat_electrode_log->SetVisAttributes(vis_att);
    }

    //position vector
    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4ThreeVector move_lower = lower_z_position * direction;
    G4ThreeVector move_upper = z_position * direction;
    G4ThreeVector move_null = G4ThreeVector(0, 0, 0);


    //physical volumes
    //  lower_electrodeMat_electrode_phys = new G4PVPlacement(/*rotate*/0, /*move*/move_null,
    //       lower_electrodeMat_electrode_log,"lower_electrodeMat_electrode_phys", germanium_vacuum_core_log,
    //       false, detector_number);
    //  upper_electrodeMat_electrode_phys = new G4PVPlacement(rotate, move_upper,
    //       "upper_electrodeMat_electrode_phys", upper_electrodeMat_electrode_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(lower_electrodeMat_electrode_log, move_lower, rotate);
    this->assembly->AddPlacedVolume(upper_electrodeMat_electrode_log, move_upper, rotate);

    return 1;
}//end ::end AddElectrodeMatElectrode


//Add the cage directly around the germanium crystal
G4int DetectionSystem8pi::AddStructureMatCage()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->structureMat);
    if( !material ) {
        G4cout << " ----> Material " << "StructureMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    vis_att->SetVisibility(true);

    // inner_cage measurements
    G4double inner_radius = crystal_outer_radius + extraClearance; //clearance
    G4double outer_radius = inner_radius + inner_can_thickness;
    G4double half_length_z = (crystal_length + inner_can_extends_past_crystal)/2.0;

    G4Tubs* inner_cage_1 = new G4Tubs("inner_cage_1", inner_radius, outer_radius, half_length_z, 0.0*deg, /*detail_view_end_angle*/15.0*deg);
    G4Tubs* inner_cage_2 = new G4Tubs("inner_cage_2", inner_radius, outer_radius, half_length_z, 90.0*deg, /*detail_view_end_angle*/15.0*deg);
    G4Tubs* inner_cage_3 = new G4Tubs("inner_cage_3", inner_radius, outer_radius, half_length_z, 180.0*deg, /*detail_view_end_angle*/15.0*deg);
    G4Tubs* inner_cage_4 = new G4Tubs("inner_cage_4", inner_radius, outer_radius, half_length_z, 270.0*deg, /*detail_view_end_angle*/15.0*deg);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin + inner_can_extends_past_crystal/2.0;
    G4ThreeVector move = z_position * direction;

    //4 logical sides (quadrants)
    if( inner_cage_1_log == NULL && inner_cage_2_log == NULL && inner_cage_3_log == NULL && inner_cage_4_log == NULL)
    {
        inner_cage_1_log = new G4LogicalVolume(inner_cage_1, material, "inner_cage_1_log", 0, 0, 0);
        inner_cage_1_log->SetVisAttributes(vis_att);
        inner_cage_2_log = new G4LogicalVolume(inner_cage_2, material, "inner_cage_2_log", 0, 0, 0);
        inner_cage_2_log->SetVisAttributes(vis_att);
        inner_cage_3_log = new G4LogicalVolume(inner_cage_3, material, "inner_cage_3_log", 0, 0, 0);
        inner_cage_3_log->SetVisAttributes(vis_att);
        inner_cage_4_log = new G4LogicalVolume(inner_cage_4, material, "inner_cage_4_log", 0, 0, 0);
        inner_cage_4_log->SetVisAttributes(vis_att);
    }

    //4 physical sides
    //  inner_cage_1_phys = new G4PVPlacement(rotate, move,
    //       "inner_cage_1_phys", inner_cage_1_log, exp_hall_phys,
    //       false, detector_number);

    //  inner_cage_2_phys = new G4PVPlacement(rotate, move,
    //       "inner_cage_2_phys", inner_cage_2_log, exp_hall_phys,
    //       false, detector_number);

    //  inner_cage_3_phys = new G4PVPlacement(rotate, move,
    //       "inner_cage_3_phys", inner_cage_3_log, exp_hall_phys,
    //       false, detector_number);

    //  inner_cage_4_phys = new G4PVPlacement(rotate, move,
    //       "inner_cage_4_phys", inner_cage_4_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(inner_cage_1_log, move, rotate);
    this->assembly->AddPlacedVolume(inner_cage_2_log, move, rotate);
    this->assembly->AddPlacedVolume(inner_cage_3_log, move, rotate);
    this->assembly->AddPlacedVolume(inner_cage_4_log, move, rotate);

    //add bottom to inner_can
    inner_radius = /*0.0*mm*/ outer_radius - 2.0*mm; //2.0mm 'lip'
    half_length_z = inner_can_thickness/2.0;

    G4Tubs* inner_cage_bottom = new G4Tubs("inner_cage_bottom", inner_radius,
                                           outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    z_position = crystal_dist_from_origin - crystal_length/2.0 - half_length_z;
    move =  z_position * direction;

    //logical volume
    if( inner_cage_bottom_log == NULL )
    {
        inner_cage_bottom_log = new G4LogicalVolume(inner_cage_bottom, material, "inner_cage_bottom_log", 0, 0, 0);
        inner_cage_bottom_log->SetVisAttributes(vis_att);
    }

    //  inner_cage_bottom_phys = new G4PVPlacement(rotate, move,
    //       "inner_cage_bottom_phys", inner_cage_bottom_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(inner_cage_bottom_log, move, rotate);

    return 1;
}//end ::end AddInnerStructureMatCage

G4int DetectionSystem8pi::AddInnerStructureMatLids()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->structureMat);
    if( !material ) {
        G4cout << " ----> Material " << "StructureMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    vis_att->SetVisibility(true);

    //ring measurements
    G4double inner_radius = structureMat_cooling_rod_radius;
    G4double outer_radius = crystal_outer_radius + inner_can_thickness + extraClearance; //extra clearance
    G4double half_length_z = inner_can_lid_thickness/2.0;

    G4Tubs* lid = new G4Tubs("lid", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin + crystal_length/2.0 + inner_can_extends_past_crystal + half_length_z;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( inner_cage_lid_log == NULL )
    {
        inner_cage_lid_log = new G4LogicalVolume(lid, material, "inner_cage_lid_log", 0, 0, 0);
        inner_cage_lid_log->SetVisAttributes(vis_att);
    }

    //place 1st lid
    //  inner_cage_lid_1_phys = new G4PVPlacement(rotate, move,
    //       "inner_cage_lid_1_phys", inner_cage_lid_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(inner_cage_lid_log, move, rotate);

    //separate 2 lids
    z_position =  z_position + inner_can_lid_separation + inner_can_lid_thickness;
    move =  z_position * direction;

    //place 2nd lid
    //  inner_cage_lid_2_phys = new G4PVPlacement(rotate, move,
    //       "inner_cage_lid_2_phys", inner_cage_lid_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(inner_cage_lid_log, move, rotate);

    return 1;
}//end ::end AddInnerStructureMatLids

G4int DetectionSystem8pi::AddStructureMatCoolingRod()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->structureMat);
    if( !material ) {
        G4cout << " ----> Material " << "StructureMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    vis_att->SetVisibility(true);

    //ring measurements
    G4double inner_radius = 0.0*mm;
    G4double outer_radius = structureMat_cooling_rod_radius;
    G4double half_length_z = (outer_can_extends_past_crystal - inner_can_extends_past_crystal )/2.0;

    G4Tubs* cooling_rod = new G4Tubs("cooling_rod", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4double z_position = crystal_dist_from_origin + crystal_length/2.0 + inner_can_extends_past_crystal + half_length_z;

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( structureMat_cooling_rod_log == NULL )
    {
        structureMat_cooling_rod_log = new G4LogicalVolume(cooling_rod, material, "structureMat_cooling_rod_log", 0, 0, 0);
        structureMat_cooling_rod_log->SetVisAttributes(vis_att);
    }

    //  structureMat_cooling_rod_phys = new G4PVPlacement(rotate, move,
    //       "structureMat_cooling_rod_phys", structureMat_cooling_rod_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(structureMat_cooling_rod_log, move, rotate);

    return 1;
}//end ::end AddStructureMatCoolingRod


G4int DetectionSystem8pi::AddElectrodeMatCoolingRod()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->electrodeMat);
    if( !material ) {
        G4cout << " ----> Material " << "ElectrodeMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.75, 0.45, 0.2));
    vis_att->SetVisibility(true);

    //upper cooling rod measurements
    G4double inner_radius = 0.0*mm;
    G4double outer_radius = electrodeMat_cooling_rod_radius;
    G4double half_length_z = electrodeMat_cooling_rod_length/2.0;

    G4Tubs* cooling_rod = new G4Tubs("cooling_rod", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4double z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal + half_length_z;

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( electrodeMat_cooling_rod_log == NULL)
    {
        electrodeMat_cooling_rod_log = new G4LogicalVolume(cooling_rod, material, "electrodeMat_cooling_rod_log", 0, 0, 0);
        electrodeMat_cooling_rod_log->SetVisAttributes(vis_att);
    }

    //still need to add 'ring' section where structureMat rod penetrates electrodeMat rod

    //  electrodeMat_cooling_rod_phys = new G4PVPlacement(rotate, move,
    //       "electrodeMat_cooling_rod_phys", electrodeMat_cooling_rod_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(electrodeMat_cooling_rod_log, move, rotate);

    return 1;
}//end ::end AddElectrodeMatCoolingRod

//Add the beryllium window in front of the crystal
G4int DetectionSystem8pi::AddBerylliumWindow()
{
    //material
    G4Material* material = G4Material::GetMaterial("Beryllium");
    if( !material ) {
        G4cout << " ----> Material " << "Beryllium" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.0, 0.0));
    vis_att->SetVisibility(true);

    // beryllium window measurements
    G4double inner_radius = 0.0*cm;
    G4double outer_radius = outer_can_innerRadius; /*crystal_outer_radius + inner_can_thickness + outer_can_dist_from_inner_can;*/
    G4double half_length_z = beryllium_thickness/2.0;

    G4Tubs* beryllium_window = new G4Tubs("beryllium_window", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin -(crystal_length/2.0 + beryllium_dist_from_crystal + half_length_z);
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( beryllium_window_log == NULL )
    {
        beryllium_window_log = new G4LogicalVolume(beryllium_window, material, "beryllium_window_log", 0, 0, 0);
        beryllium_window_log->SetVisAttributes(vis_att);
    }

    //  beryllium_window_phys = new G4PVPlacement(rotate, move,
    //       "beryllium_window_phys", beryllium_window_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(beryllium_window_log, move, rotate);

    return 1;
}//end ::end AddBerylliumWindow

//Add the outer structureMat can around the germanium crystal
G4int DetectionSystem8pi::AddOuterStructureMatCan()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->structureMat);
    if( !material ) {
        G4cout << " ----> Material " << "StructureMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    vis_att->SetVisibility(true);

    // measurements
    G4double inner_radius = inner_BGO_outerRadius + inner_BGO_clearance - outer_can_thickness;
    G4double outer_radius = inner_radius + outer_can_thickness;
    G4double half_length_z = outer_can_length/2.0;

    G4Tubs* outer_can_side = new G4Tubs("outer_can_side", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin + half_length_z - crystal_length/2.0 - beryllium_dist_from_crystal - beryllium_thickness;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( outer_can_side_log == NULL )
    {
        outer_can_side_log = new G4LogicalVolume(outer_can_side, material, "outer_can_side_log", 0, 0, 0);
        outer_can_side_log->SetVisAttributes(vis_att);
    }

    //add physical outer_can_side
    //  outer_can_side_phys = new G4PVPlacement(rotate, move,
    //       "outer_can_side_phys", outer_can_side_log, exp_hall_phys,
    //       false, detector_number);
    this->assembly->AddPlacedVolume(outer_can_side_log, move, rotate);

    //add outer can lid
    inner_radius = inner_BGO_innerRadius - inner_BGO_clearance - outer_can_thickness;
    outer_radius = outer_can_innerRadius;
    half_length_z = outer_can_thickness/2.0;

    G4Tubs* outer_can_lid = new G4Tubs("outer_can_lid", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    //logical volume
    if( outer_can_lid_log == NULL )
    {
        outer_can_lid_log = new G4LogicalVolume(outer_can_lid, material, "outer_can_lid_log", 0, 0, 0);
        outer_can_lid_log->SetVisAttributes(vis_att);
    }

    //place outer_can_lid
    z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal - outer_can_thickness + half_length_z;
    move = z_position * direction;

    //  outer_can_lid_phys = new G4PVPlacement(rotate, move,
    //       "outer_can_lid_phys", outer_can_lid_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(outer_can_lid_log, move, rotate);

    return 1;
}//end ::end AddOuterStructureMatCan

G4int DetectionSystem8pi::AddCoolingRodCover()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->structureMat);
    if( !material ) {
        G4cout << " ----> Material " << "StructureMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    vis_att->SetVisibility(true);

    //cooling rod cover measurements
    G4double inner_radius = inner_BGO_innerRadius - inner_BGO_clearance - outer_can_thickness;
    G4double outer_radius = inner_radius + outer_can_thickness;
    G4double half_length_z = electrodeMat_cooling_rod_length/2.0;

    G4Tubs* cooling_rod_cover = new G4Tubs("cooling_rod_cover", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal + half_length_z;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( cooling_rod_cover_log == NULL )
    {
        cooling_rod_cover_log = new G4LogicalVolume(cooling_rod_cover, material, "cooling_rod_cover_log", 0, 0, 0);
        cooling_rod_cover_log->SetVisAttributes(vis_att);
    }

    //place physical cooling_rod_cover
    //  cooling_rod_cover_phys = new G4PVPlacement(rotate, move,
    //       "cooling_rod_cover_phys", cooling_rod_cover_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(cooling_rod_cover_log, move, rotate);

    //cooling rod lid measurements
    inner_radius = inner_radius + outer_can_thickness;
    outer_radius = inner_BGO_outerRadius + inner_BGO_clearance;
    /*crystal_outer_radius + inner_can_thickness + outer_can_dist_from_inner_can;*/

    half_length_z = outer_can_thickness/2.0;

    G4Tubs* cooling_rod_cover_lid = new G4Tubs("cooling_rod_lid", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    //logical volume
    if( cooling_rod_cover_lid_log == NULL )
    {
        cooling_rod_cover_lid_log = new G4LogicalVolume(cooling_rod_cover_lid, material, "cooling_rod_cover_lid_log", 0, 0, 0);
        cooling_rod_cover_lid_log->SetVisAttributes(vis_att);
    }

    //place upper physical cooling_rod_cover_lid
    z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal + electrodeMat_cooling_rod_length - half_length_z;
    move = z_position * direction;

    //  cooling_rod_cover_lid_phys = new G4PVPlacement(rotate, move,
    //       "cooling_rod_cover_lid_phys", cooling_rod_cover_lid_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(cooling_rod_cover_lid_log, move, rotate);

    return 1;
}//end ::end AddCoolingRodCover


//Add the inner BGO annulus around the cooling rod
G4int DetectionSystem8pi::AddInnerBGOAnnulus()
{
    //material
    G4Material* material = G4Material::GetMaterial("BGO");
    if( !material ) {
        G4cout << " ----> Material " << "BGO" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0, 0.0, 0.5));
    vis_att->SetVisibility(true);

    // measurements
    G4double inner_radius = inner_BGO_innerRadius;
    G4double outer_radius = inner_BGO_outerRadius;

    G4double half_length_z = inner_BGO_annulus_length/2.0;

    G4Tubs* inner_BGO_annulus = new G4Tubs("inner_BGO_annulus", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal + inner_BGO_clearance + half_length_z;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( inner_BGO_annulus_log == NULL )
    {
        inner_BGO_annulus_log = new G4LogicalVolume(inner_BGO_annulus, material, "8pi_inner_BGO_annulus", 0, 0, 0);
        inner_BGO_annulus_log->SetVisAttributes(vis_att);
    }

    //  inner_BGO_annulus_phys = new G4PVPlacement(rotate, move,
    //       "inner_BGO_annulus_phys", inner_BGO_annulus_log, exp_hall_phys,
    //       false, detector_number);

    this->assemblyInnerBGO->AddPlacedVolume(inner_BGO_annulus_log, move, rotate);

    return 1;
}//end ::end AddInnerBGOAnnulus

//Add the outer structureMat sheath around the inner BGO annulus
G4int DetectionSystem8pi::AddStructureMatBGOSheath()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->structureMat);
    if( !material ) {
        G4cout << " ----> Material " << "StructureMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    vis_att->SetVisibility(true);

    // measurements
    G4double inner_radius = inner_BGO_outerRadius + inner_BGO_clearance;
    G4double outer_radius = inner_radius + outer_can_thickness; //structureMat_sheath_thickness?
    G4double half_length_z = electrodeMat_cooling_rod_length/2.0; //- outer_can_thickness;

    G4Tubs* structureMat_sheath = new G4Tubs("structureMat_sheath", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal + half_length_z;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( structureMat_sheath_log == NULL )
    {
        structureMat_sheath_log = new G4LogicalVolume(structureMat_sheath, material, "structureMat_sheath_log", 0, 0, 0);
        structureMat_sheath_log->SetVisAttributes(vis_att);
    }

    //add physical outer_can_side
    //  structureMat_sheath_phys = new G4PVPlacement(rotate, move,
    //       "structureMat_sheath_phys", structureMat_sheath_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(structureMat_sheath_log, move, rotate);

    return 1;
}//end ::end AddStructureMatBGOSheath

//Add the outer BGO annulus around the detector
G4int DetectionSystem8pi::AddOuterBGOAnnulus()
{
    //material
    G4Material* material = G4Material::GetMaterial("BGO");
    if( !material ) {
        G4cout << " ----> Material " << "BGO" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0, 0.0, 0.5));
    vis_att->SetVisibility(true);

    //BOTTOM (cone) portion of outer BGO
    // measurements
    G4double inner_radius = inner_BGO_outerRadius + inner_BGO_clearance + outer_can_thickness + outer_BGO_clearance;
    G4double lower_outer_radius = inner_radius + outer_BGO_bottom_thickness;
    G4double upper_outer_radius = outer_BGO_top_outer_radius;
    G4double half_length_z = outer_BGO_taper_height/2.0;

    G4Cons* outer_lower_BGO_annulus = new G4Cons("outer_lower_BGO_annulus", inner_radius, lower_outer_radius, inner_radius, upper_outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin + half_length_z - crystal_length/2.0 - beryllium_dist_from_crystal - beryllium_thickness - outer_BGO_displacement;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( outer_lower_BGO_annulus_log == NULL )
    {
        outer_lower_BGO_annulus_log = new G4LogicalVolume(outer_lower_BGO_annulus, material, "8pi_outer_lower_BGO_annulus", 0, 0, 0);
        outer_lower_BGO_annulus_log->SetVisAttributes(vis_att);
    }

    //  outer_lower_BGO_annulus_phys = new G4PVPlacement(rotate, move,
    //       "outer_lower_BGO_annulus_phys", outer_lower_BGO_annulus_log, exp_hall_phys,
    //       false, detector_number);

    this->assemblyOuterLowerBGO->AddPlacedVolume(outer_lower_BGO_annulus_log, move, rotate);

    //UPPER (annulus) portion of outer BGO
    half_length_z = (outer_BGO_total_length - outer_BGO_taper_height)/2.0;

    G4Tubs* outer_upper_BGO_annulus = new G4Tubs("outer_upper_BGO_annulus", inner_radius, upper_outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    z_position = crystal_dist_from_origin + half_length_z + outer_BGO_taper_height - crystal_length/2.0 - beryllium_dist_from_crystal - beryllium_thickness - outer_BGO_displacement;
    move = z_position * direction;

    //logical volume
    if( outer_upper_BGO_annulus_log == NULL )
    {
        outer_upper_BGO_annulus_log = new G4LogicalVolume(outer_upper_BGO_annulus, material, "8pi_outer_upper_BGO_annulus", 0, 0, 0);
        outer_upper_BGO_annulus_log->SetVisAttributes(vis_att);
    }

    //  outer_upper_BGO_annulus_phys = new G4PVPlacement(rotate, move,
    //       "outer_upper_BGO_annulus_phys", outer_upper_BGO_annulus_log, exp_hall_phys,
    //       false, detector_number);

    this->assemblyOuterUpperBGO->AddPlacedVolume(outer_upper_BGO_annulus_log, move, rotate);

    return 1;
}//end ::end AddOuterBGOAnnulus

//Add the liquid N2 container at top of detector
G4int DetectionSystem8pi::AddLiquidN2Container()
{
    //material
    G4Material* material = G4Material::GetMaterial("Liquid_N2");
    if( !material ) {
        G4cout << " ----> Material " << "Liquid_N2" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* N2_vis_att = new G4VisAttributes(G4Colour(0.0, 0.5, 0.5));
    N2_vis_att->SetVisibility(true);
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    vis_att->SetVisibility(true);

    // measurements
    G4double inner_radius = 0;
    G4double outer_radius = liquid_N2_radius - outer_can_thickness;
    G4double half_length_z = (liquid_N2_length - outer_can_thickness*2.0)/2.0;

    G4Tubs* liquid_N2 = new G4Tubs("liquid_N2", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal + electrodeMat_cooling_rod_length + half_length_z;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( liquid_N2_log == NULL )
    {
        liquid_N2_log = new G4LogicalVolume(liquid_N2, material, "liquid_N2_log", 0, 0, 0);
        liquid_N2_log->SetVisAttributes(N2_vis_att);
    }

    //add physical liquid N2
    //  liquid_N2_phys = new G4PVPlacement(rotate, move,
    //       "liquid_N2_phys", liquid_N2_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(liquid_N2_log, move, rotate);

    //material
    material = G4Material::GetMaterial(this->structureMat);
    if( !material ) {
        G4cout << " ----> Material " << "StructureMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    //place structureMat container surrounding N2
    //bottom section
    inner_radius = inner_BGO_outerRadius + inner_BGO_clearance + outer_can_thickness;
    outer_radius = liquid_N2_radius - outer_can_thickness;
    half_length_z = outer_can_thickness/2.0;

    G4Tubs* liquid_N2_bottom = new G4Tubs("liquid_N2_bottom", inner_radius, outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);


    z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal + electrodeMat_cooling_rod_length + half_length_z - outer_can_thickness;

    move = z_position * direction;

    //logical volume
    if( liquid_N2_bottom_log == NULL)
    {
        liquid_N2_bottom_log = new G4LogicalVolume(liquid_N2_bottom, material, "liquid_N2_bottom_log", 0, 0, 0);
        liquid_N2_bottom_log->SetVisAttributes(vis_att);
    }

    //  liquid_N2_bottom_phys = new G4PVPlacement(rotate, move,
    //       "liquid_N2_bottom_phys", liquid_N2_bottom_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(liquid_N2_bottom_log, move, rotate);

    //side section
    inner_radius = liquid_N2_radius - outer_can_thickness;
    outer_radius = liquid_N2_radius;
    half_length_z = liquid_N2_length/2.0;

    G4Tubs* liquid_N2_side = new G4Tubs("liquid_N2_side", inner_radius,
                                        outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);


    z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal + electrodeMat_cooling_rod_length + half_length_z - outer_can_thickness;
    move = z_position * direction;

    //logical volume
    if( liquid_N2_side_log == NULL)
    {
        liquid_N2_side_log = new G4LogicalVolume(liquid_N2_side, material, "liquid_N2_side_log", 0, 0, 0);
        liquid_N2_side_log->SetVisAttributes(vis_att);
    }

    //  liquid_N2_side_phys = new G4PVPlacement(rotate, move,
    //       "liquid_N2_side_phys", liquid_N2_side_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(liquid_N2_side_log, move, rotate);

    //top (lid) section
    inner_radius = 0;
    outer_radius = liquid_N2_radius - outer_can_thickness;
    half_length_z = outer_can_thickness/2.0;

    G4Tubs* liquid_N2_lid = new G4Tubs("liquid_N2_lid", inner_radius,
                                       outer_radius, half_length_z, 0.0*deg, detail_view_end_angle);


    z_position = crystal_dist_from_origin + crystal_length/2.0 + outer_can_extends_past_crystal + electrodeMat_cooling_rod_length + liquid_N2_length + half_length_z - outer_can_thickness*2.0;

    move = z_position * direction;

    //logical volume
    if( liquid_N2_lid_log == NULL)
    {
        liquid_N2_lid_log = new G4LogicalVolume(liquid_N2_lid, material, "liquid_N2_lid_log", 0, 0, 0);
        liquid_N2_lid_log->SetVisAttributes(vis_att);
    }

    //  liquid_N2_lid_phys = new G4PVPlacement(rotate, move,
    //       "liquid_N2_lid_phys", liquid_N2_lid_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(liquid_N2_lid_log, move, rotate);

    return 1;
}//end ::end AddLiquidN2Container


G4int DetectionSystem8pi::AddHevimetalCollimator()
{
    //material
    G4Material* material = G4Material::GetMaterial("Hevimetal");
    if( !material ) {
        G4cout << " ----> Material " << "Hevimetal" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* hevimetal_vis_att = new G4VisAttributes(G4Colour(0.25, 0.25, 0.25));
    hevimetal_vis_att->SetVisibility(true);

    // hevimetal measurements
    G4double front_tangent_dist = hevimetal_front_side_length * tan(60.0*deg) / 2.0; //tangent distance to inner surface
    G4double rear_tangent_dist =  hevimetal_rear_side_length * tan(60.0*deg) / 2.0;
    G4double half_length_z = hevimetal_thickness/2.0;

    G4double phiStart = 0.0*deg;		//start angle
    //G4double phiEnd = 360.0*deg;		//end angle
    G4double phiEnd = detail_view_end_angle;
    G4int		 numSide = 6;						//number of sides => 6: normal hexagon
    G4int		 numZPlanes = 2; 				//number of Z planes
    G4double zPlane[] = {-half_length_z, half_length_z}; //distance between Z planes
    G4double rInner[] = {0.0, 0.0}; //solid inside
    G4double rOuter[] = {front_tangent_dist, rear_tangent_dist}; //front face and rear face tangent distances

    G4Polyhedra* hevimetal = new G4Polyhedra("hevimetal", phiStart, phiEnd, numSide, numZPlanes, zPlane, rInner, rOuter);

    // Cut Out AuxMat Plug! ///////////////////////////////
    //Place a conical slice of auxMat in the hevimetal
    G4double inner_radius_neg = 0.0*cm;
    G4double outer_radius_neg = auxMat_plug_neg_radius + this->cut_clearance; // 0.1*mm cut clearance
    G4double inner_radius_pos = 0.0*cm;
    G4double outer_radius_pos = auxMat_plug_pos_radius + this->cut_clearance; // 0.1*mm cut clearance

    G4double new_half_length_z = (hevimetal_thickness + this->cut_clearance)/2.0; // 0.1*mm cut clearance

    G4Cons* auxMat_plug = new G4Cons("auxMat_plug", inner_radius_neg, outer_radius_neg, inner_radius_pos, outer_radius_pos, new_half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector move_cut(0,0,0);

    G4SubtractionSolid* hevimetal_with_cavity = new G4SubtractionSolid("hevimetal_with_cavity", hevimetal, auxMat_plug, 0, move_cut);
    ////////////////////////////////////////////////////////

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin - crystal_length/2.0 - beryllium_dist_from_crystal - beryllium_thickness - outer_BGO_displacement - half_length_z;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( hevimetal_log == NULL)
    {
        hevimetal_log = new G4LogicalVolume(hevimetal_with_cavity, material, "hevimetal_log", 0, 0, 0);
        hevimetal_log->SetVisAttributes(hevimetal_vis_att);
    }

    //     hevimetal_phys = new G4PVPlacement(rotate, move,
    //       "hevimetal_phys", hevimetal_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(hevimetal_log, move, rotate);

    return 1;
}//end ::end AddHevimetalCollimator


G4int DetectionSystem8pi::AddAuxMatPlug()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->auxMat);
    if( !material ) {
        G4cout << " ----> Material " << "AuxMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
    vis_att->SetVisibility(true);

    //Place a conical slice of auxMat in the hevimetal
    G4double inner_radius_neg = 0.0*cm;
    G4double outer_radius_neg = auxMat_plug_neg_radius;
    G4double inner_radius_pos = 0.0*cm;
    G4double outer_radius_pos = auxMat_plug_pos_radius;

    G4double half_length_z = hevimetal_thickness/2.0;

    G4Cons* auxMat_plug = new G4Cons("auxMat_plug", inner_radius_neg, outer_radius_neg, inner_radius_pos, outer_radius_pos, half_length_z, 0.0*deg, detail_view_end_angle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin - crystal_length/2.0 - beryllium_dist_from_crystal - beryllium_thickness - outer_BGO_displacement - half_length_z;

    //Nest the auxMat plug in the hevimetal
    G4ThreeVector move = z_position * direction ; // check this!

    //logical volume
    if( auxMat_plug_log == NULL)
    {
        auxMat_plug_log = new G4LogicalVolume(auxMat_plug, material, "auxMat_plug_log", 0, 0, 0);
        auxMat_plug_log->SetVisAttributes(vis_att);
    }

    G4ThreeVector move_null = G4ThreeVector(0, 0, 0);

    //  auxMat_plug_phys = new G4PVPlacement(/*no rotate*/ 0, /*no move*/ move_null,
    //       auxMat_plug_log, "auxMat_plug_phys", hevimetal_log,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(auxMat_plug_log, move, rotate);

    return 1;
}//end ::end AddAuxMatPlug

G4int DetectionSystem8pi::AddThinAuxMatLayer()
{
    //material
    G4Material* material = G4Material::GetMaterial(this->auxMat);
    if( !material ) {
        G4cout << " ----> Material " << "AuxMat" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
    vis_att->SetVisibility(true);

    // auxMat layer measurements
    G4double front_tangent_dist = auxMat_layer_front_side_length * tan(60.0*deg) / 2.0;//front tangent distance to edge
    G4double rear_tangent_dist =  hevimetal_front_side_length * tan(60.0*deg) / 2.0; //rear tangent distance to edge
    G4double half_length_z = auxMat_layer_thickness/2.0;

    G4double phiStart = 0.0*deg;		//start angle
    //G4double phiEnd = 360.0*deg;		//end angle
    G4double phiEnd = detail_view_end_angle;
    G4int		 numSide = 6;						//number of sides => 6: normal hexagon
    G4int		 numZPlanes = 2; 				//number of Z planes
    G4double zPlane[] = {-half_length_z, half_length_z}; //distance between Z planes
    G4double rInner[] = {0.0, 0.0}; //solid inside
    G4double rOuter[] = {front_tangent_dist, rear_tangent_dist}; //front face and rear face tangent distances


    G4Polyhedra* auxMat_layer = new G4Polyhedra("auxMat_layer", phiStart, phiEnd, numSide, numZPlanes, zPlane, rInner, rOuter);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double z_position = crystal_dist_from_origin - crystal_length/2.0 - beryllium_dist_from_crystal - beryllium_thickness - outer_BGO_displacement - hevimetal_thickness - half_length_z;
    G4ThreeVector move = z_position * direction;

    //logical volume
    if( auxMat_layer_log == NULL)
    {
        auxMat_layer_log = new G4LogicalVolume(auxMat_layer, material, "auxMat_layer_log", 0, 0, 0);
        auxMat_layer_log->SetVisAttributes(vis_att);
    }

    //  auxMat_layer_phys = new G4PVPlacement(rotate, move,
    //       "auxMat_layer_phys", auxMat_layer_log, exp_hall_phys,
    //       false, detector_number);

    this->assembly->AddPlacedVolume(auxMat_layer_log, move, rotate);

    return 1;
}//end ::end AddThinAuxMatLayer
