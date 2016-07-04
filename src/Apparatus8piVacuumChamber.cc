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
    delete fVacuumChamberSphereLog;

}// end ::~Apparatus8piVacuumChamber

///////////////////////////////////////////////////////////////////////
//ConstructApparatus8piVacuumChamber builds the vacuum Chamber at the origin
///////////////////////////////////////////////////////////////////////
void Apparatus8piVacuumChamber::Build(G4LogicalVolume* expHallLog)
{ 
    fExpHallLog = expHallLog;

    BuildApparatus8piVacuumChamberSphere();               //Includes a vacuum

    PlaceApparatus8piVacuumChamberSphere();

}//end ::Build

void Apparatus8piVacuumChamber::BuildApparatus8piVacuumChamberSphere()
{
    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.4,0.2,1.0));
    visAtt->SetVisibility(true);

    G4double startPhi = 0;
    G4double endPhi = 2*M_PI;

    G4double startTheta = 0;
    G4double endTheta = M_PI;

    G4double innerRadius = fVacuumChamberInnerRadius;
    G4double outerRadius = fVacuumChamberOuterRadius;

    G4Sphere* fVacuumChamberSphere = new G4Sphere("VacuumChamberSphere",  innerRadius, outerRadius, startPhi, endPhi, startTheta, endTheta);
    G4Material* vacuumChamberSphereMat = G4Material::GetMaterial(fVacuumChamberSphereMaterial);
    fVacuumChamberSphereLog = new G4LogicalVolume(fVacuumChamberSphere, vacuumChamberSphereMat, "VacuumChamberSphereLog", 0, 0, 0);
    fVacuumChamberSphereLog->SetVisAttributes(visAtt);

}//end ::Apparatus8piVacuumChamberCylinder


///////////////////////////////////////////////////////////////////////
//methods used in Build()
///////////////////////////////////////////////////////////////////////

void Apparatus8piVacuumChamber::PlaceApparatus8piVacuumChamberSphere()
{
    G4double zPosition;
    zPosition = 0.0;
    G4ThreeVector move(0, 0, zPosition);
    // Establish physical volumes
    fVacuumChamberSpherePhys = new G4PVPlacement(0, move, fVacuumChamberSphereLog, "VacuumChamberSpherePhys", fExpHallLog, false, 0);
}//end ::PlaceApparatus8piVacuumChamberCylinder()
