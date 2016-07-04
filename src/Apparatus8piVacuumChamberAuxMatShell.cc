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
    delete fVacuumChamberAuxSphereLog;

}// end ::~Apparatus8piVacuumChamberAuxMatShell

///////////////////////////////////////////////////////////////////////
//ConstructApparatus8piVacuumChamberAuxMatShell builds the vacuum Chamber at the origin
///////////////////////////////////////////////////////////////////////
void Apparatus8piVacuumChamberAuxMatShell::Build(G4LogicalVolume* expHallLog, G4double thickness)
{ 
    fExpHallLog = expHallLog;

    BuildApparatus8piVacuumChamberAuxMatShellSphere(thickness);

    PlaceApparatus8piVacuumChamberAuxMatShellSphere();

}//end ::Build

void Apparatus8piVacuumChamberAuxMatShell::BuildApparatus8piVacuumChamberAuxMatShellSphere(G4double thickness)
{
    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.4,0.4,1.0));
    visAtt->SetVisibility(true);

    G4double startPhi = 0;
    G4double endPhi = 2*M_PI;

    G4double startTheta = 0;
    G4double endTheta = M_PI;

    G4double innerRadius = fVacuumChamberOuterRadius;
    G4double outerRadius = fVacuumChamberOuterRadius + thickness;

    G4Sphere* vacuumChamberSphere = new G4Sphere("VacuumChamberSphere", innerRadius, outerRadius, startPhi, endPhi, startTheta, endTheta);
    G4Material* vacuumChamberSphereMat = G4Material::GetMaterial(fVacuumChamberSphereMaterial);

    fVacuumChamberAuxSphereLog = new G4LogicalVolume(vacuumChamberSphere, vacuumChamberSphereMat, "VacuumChamberAuxSphereLog", 0, 0, 0);
    fVacuumChamberAuxSphereLog->SetVisAttributes(visAtt);

}//end ::Apparatus8piVacuumChamberAuxMatShellCylinder


///////////////////////////////////////////////////////////////////////
//methods used in Build()
///////////////////////////////////////////////////////////////////////

void Apparatus8piVacuumChamberAuxMatShell::PlaceApparatus8piVacuumChamberAuxMatShellSphere()
{
    G4double zPosition;
    zPosition = 0.0;
    G4ThreeVector move(0, 0, zPosition);
    // Establish physical volumes
    fVacuumChamberSpherePhys = new G4PVPlacement(0, move, fVacuumChamberAuxSphereLog, "VacuumChamberSpherePhys", fExpHallLog, false, 0);
}//end ::PlaceApparatus8piVacuumChamberAuxMatShellCylinder()
