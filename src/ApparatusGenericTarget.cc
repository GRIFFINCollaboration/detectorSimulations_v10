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
    fTargetLog(0)
{ 
  fTargetMaterial = "Gold";
  fTargetLengthX = 0.0*cm; 
  fTargetLengthY = 0.0*cm;
  fTargetLengthZ = 0.0*cm;
}

ApparatusGenericTarget::~ApparatusGenericTarget()
{
    // LogicalVolumes 
    delete fTargetLog;    
}

void ApparatusGenericTarget::Build(G4String target_material_in, G4double target_length_x_in, G4double target_length_y_in, G4double target_length_z_in)
{
  fTargetMaterial = target_material_in;
  fTargetLengthX = target_length_x_in*mm;
  fTargetLengthY = target_length_y_in*mm;
  fTargetLengthZ = target_length_z_in*mm;

  // Build assembly volume
  fAssembly = new G4AssemblyVolume();

  G4cout<<"BuildTargetVolume"<<G4endl;
  BuildTargetVolume();      
}

void ApparatusGenericTarget::PlaceApparatus(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate)
{
  G4int copy_ID = 0;

  fAssembly->MakeImprint(exp_hall_log, move, rotate, copy_ID);
}

void ApparatusGenericTarget::BuildTargetVolume()
{
  G4Material* material = G4Material::GetMaterial(fTargetMaterial);
  if(!material) {
    G4cout<<" ----> Material "<<fTargetMaterial<<" not found, cannot build!"<<G4endl;
    throw std::invalid_argument("Material not found, cannot build!");
  }
  
  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.8,0.0,0.2));
  vis_att->SetVisibility(true);  

  G4ThreeVector move; 
  G4RotationMatrix* rotate = new G4RotationMatrix;

  G4Box* target = BuildTarget(); 

  //logical volume
  if(fTargetLog == nullptr)
  {
    fTargetLog = new G4LogicalVolume(target, material, "target_log", 0, 0, 0);
    fTargetLog->SetVisAttributes(vis_att);
  }

  fAssembly->AddPlacedVolume(fTargetLog, move, rotate);
}

G4Box* ApparatusGenericTarget::BuildTarget()
{
  G4double half_length_x = fTargetLengthX/2.0;
  G4double half_length_y = fTargetLengthY/2.0;
  G4double half_length_z = fTargetLengthZ/2.0;

  return new G4Box("target", half_length_x, half_length_y, half_length_z);
}

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector ApparatusGenericTarget::GetDirectionXYZ(G4double theta, G4double phi)
{
  G4double x,y,z;
  x = sin(theta) * cos(phi);
  y = sin(theta) * sin(phi);
  z = cos(theta);
	
  return G4ThreeVector(x,y,z);
}
