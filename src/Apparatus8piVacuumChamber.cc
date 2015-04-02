#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "Apparatus8piVacuumChamber.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

//constructor suppressed

Apparatus8piVacuumChamber::~Apparatus8piVacuumChamber()
{
    // LogicalVolumes in ConstructApparatus8piVacuumChamber
    delete vacuum_chamber_sphere_log;
    //   delete vacuum_chamber_sphere_vacuum_log;

}// end ::~Apparatus8piVacuumChamber

///////////////////////////////////////////////////////////////////////
//ConstructApparatus8piVacuumChamber builds the vacuum Chamber at the origin
///////////////////////////////////////////////////////////////////////
void Apparatus8piVacuumChamber::Build(G4LogicalVolume* exp_hall_log)
{ 
    this->expHallLog = exp_hall_log;

    BuildApparatus8piVacuumChamberSphere();               //Includes a vacuum

    PlaceApparatus8piVacuumChamberSphere();

}//end ::Build

void Apparatus8piVacuumChamber::BuildApparatus8piVacuumChamberSphere()
{
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.4,0.2,1.0));
    vis_att->SetVisibility(true);

    G4double startPhi = 0;
    G4double endPhi = 2*M_PI;

    G4double startTheta = 0;
    G4double endTheta = M_PI;

    G4double inner_radius = this->vacuum_chamber_inner_radius;
    G4double outer_radius = this->vacuum_chamber_outer_radius;

    G4Sphere* vacuum_chamber_sphere = new G4Sphere("vacuum_chamber_sphere",  inner_radius, outer_radius, startPhi, endPhi, startTheta, endTheta);
    G4Material* vacuum_chamber_sphere_material = G4Material::GetMaterial(this->vacuum_chamber_sphere_material);
    vacuum_chamber_sphere_log = new G4LogicalVolume(vacuum_chamber_sphere, vacuum_chamber_sphere_material, "vacuum_chamber_sphere_log", 0, 0, 0);
    vacuum_chamber_sphere_log->SetVisAttributes(vis_att);

    //   // Fill with Vacuum
    //   inner_radius = 0.0*mm;
    //   outer_radius = this->vacuum_chamber_inner_radius;

    //   G4Sphere* vacuum_chamber_sphere_vacuum = new G4Sphere("vacuum_chamber_sphere_vacuum",  inner_radius, outer_radius, startPhi, endPhi, startTheta, endTheta);
    //   G4Material* vacuum_chamber_sphere_vacuum_material = G4Material::GetMaterial(this->vacuum_material);
    //   vacuum_chamber_sphere_vacuum_log = new G4LogicalVolume(vacuum_chamber_sphere_vacuum, vacuum_chamber_sphere_vacuum_material, "vacuum_chamber_sphere_vacuum_log", 0, 0, 0);
    //   vacuum_chamber_sphere_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);

}//end ::Apparatus8piVacuumChamberCylinder


///////////////////////////////////////////////////////////////////////
//methods used in Build()
///////////////////////////////////////////////////////////////////////

void Apparatus8piVacuumChamber::PlaceApparatus8piVacuumChamberSphere()
{
    G4double z_position;
    z_position = 0.0;
    G4ThreeVector move(0, 0, z_position);
    // Establish physical volumes
    vacuum_chamber_sphere_phys = new G4PVPlacement(0, move, vacuum_chamber_sphere_log, "vacuum_chamber_sphere_phys", expHallLog, false, 0);
    //   vacuum_chamber_sphere_vacuum_phys = new G4PVPlacement(0, move, vacuum_chamber_sphere_vacuum_log, "vacuum_chamber_sphere_vacuum_phys", expHallLog, false, 0);
}//end ::PlaceApparatus8piVacuumChamberCylinder()
