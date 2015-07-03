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

#include "Apparatus8piVacuumChamberAuxMatShell.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

//constructor suppressed

Apparatus8piVacuumChamberAuxMatShell::~Apparatus8piVacuumChamberAuxMatShell()
{
    // LogicalVolumes in ConstructApparatus8piVacuumChamberAuxMatShell
    delete vacuum_chamber_aux_sphere_log;

}// end ::~Apparatus8piVacuumChamberAuxMatShell

///////////////////////////////////////////////////////////////////////
//ConstructApparatus8piVacuumChamberAuxMatShell builds the vacuum Chamber at the origin
///////////////////////////////////////////////////////////////////////
void Apparatus8piVacuumChamberAuxMatShell::Build(G4LogicalVolume* exp_hall_log, G4double thickness)
{ 
    this->expHallLog = exp_hall_log;

    BuildApparatus8piVacuumChamberAuxMatShellSphere(thickness);

    PlaceApparatus8piVacuumChamberAuxMatShellSphere();

}//end ::Build

void Apparatus8piVacuumChamberAuxMatShell::BuildApparatus8piVacuumChamberAuxMatShellSphere(G4double thickness)
{
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.4,0.4,1.0));
    vis_att->SetVisibility(true);

    G4double startPhi = 0;
    G4double endPhi = 2*M_PI;

    G4double startTheta = 0;
    G4double endTheta = M_PI;

    G4double inner_radius = this->vacuum_chamber_outer_radius;
    G4double outer_radius = this->vacuum_chamber_outer_radius + thickness;

    G4Sphere* vacuum_chamber_sphere = new G4Sphere("vacuum_chamber_sphere",  inner_radius, outer_radius, startPhi, endPhi, startTheta, endTheta);
    G4Material* vacuum_chamber_sphere_material = G4Material::GetMaterial(this->vacuum_chamber_sphere_material);

    vacuum_chamber_aux_sphere_log = new G4LogicalVolume(vacuum_chamber_sphere, vacuum_chamber_sphere_material, "vacuum_chamber_aux_sphere_log", 0, 0, 0);
    vacuum_chamber_aux_sphere_log->SetVisAttributes(vis_att);

}//end ::Apparatus8piVacuumChamberAuxMatShellCylinder


///////////////////////////////////////////////////////////////////////
//methods used in Build()
///////////////////////////////////////////////////////////////////////

void Apparatus8piVacuumChamberAuxMatShell::PlaceApparatus8piVacuumChamberAuxMatShellSphere()
{
    G4double z_position;
    z_position = 0.0;
    G4ThreeVector move(0, 0, z_position);
    // Establish physical volumes
    vacuum_chamber_sphere_phys = new G4PVPlacement(0, move, vacuum_chamber_aux_sphere_log, "vacuum_chamber_sphere_phys", expHallLog, false, 0);
}//end ::PlaceApparatus8piVacuumChamberAuxMatShellCylinder()
