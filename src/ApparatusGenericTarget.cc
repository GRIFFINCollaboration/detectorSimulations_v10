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

#include "ApparatusGenericTarget.hh"

ApparatusGenericTarget::ApparatusGenericTarget() :
    // LogicalVolumes 
    target_log(0)
{ 
  this->target_material = "Gold";
  this->target_length_x = 0.0*cm; 
  this->target_length_y = 0.0*cm;
  this->target_length_z = 0.0*cm;
}

ApparatusGenericTarget::~ApparatusGenericTarget()
{
    // LogicalVolumes 
    delete target_log;    
}

G4int ApparatusGenericTarget::Build(G4String target_material, G4double target_length_x, G4double target_length_y, G4double target_length_z)
{
  this->target_material = target_material;
  this->target_length_x = target_length_x*mm;
  this->target_length_y = target_length_y*mm;
  this->target_length_z = target_length_z*mm;

  // Build assembly volume
  G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
  this->assembly = myAssembly;

  G4cout << "BuildTargetVolume" << G4endl;
  BuildTargetVolume();      

  return 1;
}

G4int ApparatusGenericTarget::PlaceApparatus(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate)
{
  G4int copy_ID = 0;

  this->assembly->MakeImprint(exp_hall_log, move, rotate, copy_ID);

  return 1;
}

G4int ApparatusGenericTarget::BuildTargetVolume()
{
  G4Material* material = G4Material::GetMaterial(this->target_material);
  if( !material ) {
    G4cout << " ----> Material " << this->target_material << " not found, cannot build!" << G4endl;
    return 0;
  }
  
  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.8,0.0,0.2));
  vis_att->SetVisibility(true);  

  G4ThreeVector move; 
  G4RotationMatrix* rotate = new G4RotationMatrix;

  G4Box* target = BuildTarget(); 

  //logical volume
  if( target_log == NULL )
  {
    target_log = new G4LogicalVolume(target, material, "target_log", 0, 0, 0);
    target_log->SetVisAttributes(vis_att);
  }

  this->assembly->AddPlacedVolume(target_log, move, rotate);

  return 1;
}

G4Box* ApparatusGenericTarget::BuildTarget()
{
  G4double half_length_x = this->target_length_x/2.0;
  G4double half_length_y = this->target_length_y/2.0;
  G4double half_length_z = this->target_length_z/2.0;

  G4Box* target = new G4Box("target", half_length_x, half_length_y, half_length_z);
  return target;
}

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector ApparatusGenericTarget::GetDirectionXYZ(G4double theta, G4double phi)
{
  G4double x,y,z;
  x = sin(theta) * cos(phi);
  y = sin(theta) * sin(phi);
  z = cos(theta);
	
  G4ThreeVector direction = G4ThreeVector(x,y,z);
	
  return direction;
}
