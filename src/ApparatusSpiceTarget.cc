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
  
  fBeamPos = beaminput*mm;
  //G4cout << "\n\n\n\n\n\ntarget z in spice target: " << fBeamPos << "\n\n\n\n\n\n"<< G4endl; 
  this->fTargetMaterial = "Gold";
  this->fTargetBackerMaterial = "Gold";
  this->fTargetProtectorMaterial = "Gold";//initialising - may remove
}

ApparatusSpiceTarget::~ApparatusSpiceTarget()
{ 
  // LogicalVolumes 
    delete fTargetSpiceLog;
    delete fTargetBackerSpiceLog;
    delete fTargetProtectorSpiceLog;
    delete fTargetPhys;
    delete fTargetBackerPhys;
    delete fTargetProtectorPhys;
}

G4int ApparatusSpiceTarget::BuildTarget(G4String input_material, G4double input_surface_density, G4double input_density)
{
    //G4cout << "\n\n\n\n\n\ntarget z in spice target build: " << fBeamPos << "\n\n\n\n\n\n"<< G4endl; 
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.8,0.0,0.2));
    vis_att->SetVisibility(true);  

    //Get Materials 
    this->fTargetMaterial = input_material;
    this->fTargetMaterialDensity = input_density*(g/cm3);//need units
    G4cout << "Post-units: " << fTargetMaterialDensity << G4endl;
    this->fTargetSurfaceDensity = input_surface_density*(mg/cm2);//need units
    G4cout << "Post-units: " << fTargetSurfaceDensity << G4endl;
    this->fTargetThickness = (fTargetSurfaceDensity/fTargetMaterialDensity);
    G4cout << "Thicc: " << fTargetThickness << G4endl;

    // Build assembly volume
    //G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    //this->assembly = myAssembly;
   
    G4Material* material = G4Material::GetMaterial(this->fTargetMaterial);
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    if( !material ) {
      G4cout << " ----> Target Material " << this->fTargetMaterial << " not found, cannot build!" << G4endl;
      return 0;
    } 
    
    G4Tubs* target = new G4Tubs("Spice_Target", 0.0, fTargetRadius, fTargetThickness/2.,
                0.0, 360*deg);
    
    //logical volume
      fTargetSpiceLog = new G4LogicalVolume(target, material, "target_spice_log", 0, 0, 0);
      fTargetSpiceLog->SetVisAttributes(vis_att);
      return 1;
}

G4int ApparatusSpiceTarget::BuildBacker(G4String input_material, G4double input_surface_density, G4double input_density)
{
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.2,0.8));
    vis_att->SetVisibility(true);  
    
    this->fTargetBackerMaterial = input_material;
    this->fTargetBackerMaterialDensity = input_density*(g/cm3);//need units
    this->fTargetBackerSurfaceDensity = input_surface_density*(mg/cm2);//need units
    this->fTargetBackerThickness = (fTargetBackerSurfaceDensity/fTargetBackerMaterialDensity);
    
    G4Material* material = G4Material::GetMaterial(this->fTargetBackerMaterial);
    if( !material ) {
      G4cout << " ----> Target Backer Material " << this->fTargetBackerMaterial << " not found, cannot build!" << G4endl;
      return 0;
    }
    
    G4Tubs* target_Backer = new G4Tubs("Spice_Target_Backer", 0.0, fTargetRadius, fTargetBackerThickness/2.,
                0.0, 360*deg);
    
    //logical volume
     fTargetBackerSpiceLog = new G4LogicalVolume(target_Backer, material, "target_Backer_spice_log", 0, 0, 0);
     fTargetBackerSpiceLog->SetVisAttributes(vis_att);
     return 1;
}

G4int ApparatusSpiceTarget::BuildProtector(G4String input_material, G4double input_surface_density, G4double input_density)
{
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.8,0.2));
    vis_att->SetVisibility(true);  
    
    this->fTargetProtectorMaterial = input_material;
    this->fTargetProtectorMaterialDensity = input_density*(g/cm3);//need units  
    this->fTargetProtectorSurfaceDensity = input_surface_density*(mg/cm2);//need units
    this->fTargetProtectorThickness = (fTargetProtectorSurfaceDensity/fTargetProtectorMaterialDensity);
    
    G4Material* material = G4Material::GetMaterial(this->fTargetProtectorMaterial);
    if( !material ) {
      G4cout << " ----> Target Protector Material " << this-> fTargetProtectorMaterial<< " not found, cannot build!" << G4endl;
      return 0;
    }  

    G4Tubs* target_Protector = new G4Tubs("Spice_Target_Protector", 0.0, fTargetRadius, fTargetProtectorThickness/2.,
                0.0, 360*deg);
    
    //logical volume
    fTargetProtectorSpiceLog = new G4LogicalVolume(target_Protector, material, "target_Protector_spice_log", 0, 0, 0);
    fTargetProtectorSpiceLog->SetVisAttributes(vis_att);
    return 1;
}

void ApparatusSpiceTarget::PlaceTarget(G4LogicalVolume* oexpHallLog)
{
  G4cout << "Spice Target implemented " << G4endl;
  G4cout << "targ-thick: " << fTargetThickness << " BACKER IN TARG: " <<this->fTargetBackerThickness << G4endl;
  G4double placement = fBeamPos*mm - fTargetBackerThickness - fTargetThickness/2.;//0.5mm=half thickness of target wheel where 0,0,0 is placed
  if(fBeamPos < -3.9*mm) placement = placement - 0.5*mm;
  G4cout << " Placement target : " << placement << " " << placement*mm << G4endl;
  //NEED TO LINK BEAM POS TO HERE for move = three-vector for mosition
  //fBeamPos*mm-(target_Backer_thickness)*cm
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
   
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  fTargetPhys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    fTargetSpiceLog,             //its logical volume
                    "Target_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
}

void ApparatusSpiceTarget::PlaceTargetBacker(G4LogicalVolume* oexpHallLog)
{
  G4cout << "Spice Target backer implemented " << G4endl;
  G4cout << "Back-Thick: " << fTargetBackerThickness<< G4endl;
  //NEED TO LINK BEAM POS TO HERE for move = three-vector for mosition
  G4double placement = fBeamPos*mm - (this->fTargetBackerThickness/2.);
  if(fBeamPos < -3.9*mm) placement = placement - 0.5*mm;
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  fTargetBackerPhys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    fTargetBackerSpiceLog,             //its logical volume
                    "Target_backer_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
}

void ApparatusSpiceTarget::PlaceTargetProtector(G4LogicalVolume* oexpHallLog)
{
  G4cout << "Spice Target protector implemented " << G4endl;
//   fBeamPos-(target_Backer_thickness/target_Backer_material_density)*2.0-(target_thickness/target_material_density)*2.0
  //NEED TO LINK BEAM POS TO HERE for move = three-vector for mosition
  G4cout << " Pro-thick: " << fTargetProtectorThickness << G4endl;
  G4double placement = fBeamPos*mm - (this->fTargetBackerThickness) - (this->fTargetThickness) - (this->fTargetProtectorThickness/2.);
  if(fBeamPos < -3.9*mm) placement = placement - 0.5*mm;
  G4cout << " Placement pro : " << placement << " " << placement*mm << " targetthick " << (this->fTargetThickness)<<" "<<G4endl;
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  fTargetProtectorPhys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    fTargetProtectorSpiceLog,             //its logical volume
                    "Target_Protector_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
}