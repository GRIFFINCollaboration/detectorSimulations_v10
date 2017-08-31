#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

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

#include "ApparatusSpiceTarget.hh"

ApparatusSpiceTarget::ApparatusSpiceTarget(G4double beaminput) { 
  
  beampos = beaminput*mm;
  //G4cout << "\n\n\n\n\n\ntarget z in spice target: " << beampos << "\n\n\n\n\n\n"<< G4endl; 
  this->target_material = "Gold";
  this->target_Backer_material = "Gold";
  this->target_Protector_material = "Gold";//initialising - may remove
}

ApparatusSpiceTarget::~ApparatusSpiceTarget()
{ 
  // LogicalVolumes 
    delete target_spice_log;
    delete target_Backer_spice_log;
    delete target_Protector_spice_log;
    delete target_phys;
    delete target_Backer_phys;
    delete target_Protector_phys;
}

G4int ApparatusSpiceTarget::BuildTarget(G4String input_material, G4double input_surface_density, G4double input_density)
{
    //G4cout << "\n\n\n\n\n\ntarget z in spice target build: " << beampos << "\n\n\n\n\n\n"<< G4endl; 
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.8,0.0,0.2));
    vis_att->SetVisibility(true);  

    //Get Materials 
    this->target_material = input_material;
    this->target_material_density = input_density*(g/cm3);//need units
    G4cout << "Post-units: " << target_material_density << G4endl;
    this->target_surface_density = input_surface_density*(mg/cm2);//need units
    G4cout << "Post-units: " << input_surface_density << G4endl;
    this->target_thickness = (target_surface_density/target_material_density);
    G4cout << "Thicc: " << target_thickness << G4endl;

    // Build assembly volume
    //G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    //this->assembly = myAssembly;
   
    G4Material* material = G4Material::GetMaterial(this->target_material);
      G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    if( !material ) {
      G4cout << " ----> Target Material " << this->target_material << " not found, cannot build!" << G4endl;
      return 0;
    } 
    
    G4Tubs* target = new G4Tubs("Spice_Target", 0.0, target_radius, target_thickness/2.,
                0.0, 360*deg);
    
    //logical volume
      target_spice_log = new G4LogicalVolume(target, material, "target_spice_log", 0, 0, 0);
      target_spice_log->SetVisAttributes(vis_att);
      return 1;
}

G4int ApparatusSpiceTarget::BuildBacker(G4String input_material, G4double input_surface_density, G4double input_density)
{
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.2,0.8));
    vis_att->SetVisibility(true);  
    
    this->target_Backer_material = input_material;
    this->target_Backer_material_density = input_density*(g/cm3);//need units
    this->target_Backer_surface_density = input_surface_density*(mg/cm2);//need units
    this->target_Backer_thickness = (target_Backer_surface_density/target_Backer_material_density);
    
    G4Material* material = G4Material::GetMaterial(this->target_Backer_material);
    if( !material ) {
      G4cout << " ----> Target Backer Material " << this->target_Backer_material << " not found, cannot build!" << G4endl;
      return 0;
    }
    
    G4Tubs* target_Backer = new G4Tubs("Spice_Target_Backer", 0.0, target_radius, target_Backer_thickness/2.,
                0.0, 360*deg);
    
    //logical volume
     target_Backer_spice_log = new G4LogicalVolume(target_Backer, material, "target_Backer_spice_log", 0, 0, 0);
     target_Backer_spice_log->SetVisAttributes(vis_att);
     return 1;
}

G4int ApparatusSpiceTarget::BuildProtector(G4String input_material, G4double input_surface_density, G4double input_density)
{
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.8,0.2));
    vis_att->SetVisibility(true);  
    
    this->target_Protector_material = input_material;
    this->target_Protector_material_density = input_density*(g/cm3);//need units  
    this->target_Protector_surface_density = input_surface_density*(mg/cm2);//need units
    this->target_Protector_thickness = (target_Protector_surface_density/target_Protector_material_density);
    
    G4Material* material = G4Material::GetMaterial(this->target_Protector_material);
    if( !material ) {
      G4cout << " ----> Target Protector Material " << this->target_Protector_material << " not found, cannot build!" << G4endl;
      return 0;
    }  

    G4Tubs* target_Protector = new G4Tubs("Spice_Target_Protector", 0.0, target_radius, target_Protector_thickness/2.,
                0.0, 360*deg);
    
    //logical volume
    target_Protector_spice_log = new G4LogicalVolume(target_Protector, material, "target_Protector_spice_log", 0, 0, 0);
    target_Protector_spice_log->SetVisAttributes(vis_att);
    return 1;
}

void ApparatusSpiceTarget::PlaceTarget(G4LogicalVolume* oexpHallLog)
{
  G4cout << "Spice Target implemented " << G4endl;
  G4cout << "targ-thick: " << target_thickness << " BACKER IN TARG: " <<this->target_Backer_thickness << G4endl;
  G4double placement = beampos*mm - target_Backer_thickness - target_thickness/2.;//0.5mm=half thickness of target wheel where 0,0,0 is placed
  if(beampos < -3.9*mm) placement = placement - 0.5*mm;
  G4cout << " Placement target : " << placement << " " << placement*mm << G4endl;
  //NEED TO LINK BEAM POS TO HERE for move = three-vector for mosition
  //beampos*mm-(target_Backer_thickness)*cm
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
   
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  target_phys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    target_spice_log,             //its logical volume
                    "Target_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
}

void ApparatusSpiceTarget::PlaceTargetBacker(G4LogicalVolume* oexpHallLog)
{
  G4cout << "Spice Target backer implemented " << G4endl;
  G4cout << "Back-Thick: " << target_Backer_thickness << G4endl;
  //NEED TO LINK BEAM POS TO HERE for move = three-vector for mosition
  G4double placement = beampos*mm - (this->target_Backer_thickness/2.);
  if(beampos < -3.9*mm) placement = placement - 0.5*mm;
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  target_Backer_phys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    target_Backer_spice_log,             //its logical volume
                    "Target_backer_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
}

void ApparatusSpiceTarget::PlaceTargetProtector(G4LogicalVolume* oexpHallLog)
{
  G4cout << "Spice Target protector implemented " << G4endl;
//   beampos-(target_Backer_thickness/target_Backer_material_density)*2.0-(target_thickness/target_material_density)*2.0
  //NEED TO LINK BEAM POS TO HERE for move = three-vector for mosition
  G4cout << " Pro-thick: " << target_Protector_thickness << G4endl;
  G4double placement = beampos*mm - (this->target_Backer_thickness) - (this->target_thickness) - (this->target_Protector_thickness/2.);
  if(beampos < -3.9*mm) placement = placement - 0.5*mm;
  G4cout << " Placement pro : " << placement << " " << placement*mm << " targetthick " << (this->target_thickness)<<" "<<G4endl;
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  target_Protector_phys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    target_Protector_spice_log,             //its logical volume
                    "Target_Protector_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
}