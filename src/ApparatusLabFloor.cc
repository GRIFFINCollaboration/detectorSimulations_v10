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

#include "ApparatusLabFloor.hh"

#include "G4SystemOfUnits.hh"

#include <string>

ApparatusLabFloor::ApparatusLabFloor()
{
    // LogicalVolumes
    fFloorLog = NULL;
    // can properties
    fFloorLength        = 10.*m;
    fFloorHeight        = 1.*m;
    fFloorWidth         = 10.*m;
   fFloorMaterial = "G4_CONCRETE";

}
/////
///////
ApparatusLabFloor::~ApparatusLabFloor() {
    // LogicalVolumes
delete fFloorLog;

}
////////
/////////
G4int ApparatusLabFloor::Build() {
    fAssemblyFloor = new G4AssemblyVolume(); 
    BuildFloor();

    return 1;
}
////////
///////
G4int ApparatusLabFloor::PlaceLabFloor(G4LogicalVolume* expHallLog) {
    G4RotationMatrix * rotate = new G4RotationMatrix;
    G4ThreeVector move = G4ThreeVector(0., 0., 0.);

    fAssemblyFloor->MakeImprint(expHallLog, move, rotate);
    return 1;
}
///////////
/////////
G4int ApparatusLabFloor::BuildFloor() {

    G4ThreeVector move, direction;
    G4RotationMatrix* rotate;

    G4Material* floorG4material = G4Material::GetMaterial(fFloorMaterial);
    if( !floorG4material ) {
        G4cout << " ----> Material " << fFloorMaterial << " not found, cannot build! " << G4endl;
        return 0;
    }
      else {
G4cout << floorG4material->GetName() << " is the name of the floor material" << G4endl;
}

G4Box * box = new G4Box("LabFloor", fFloorLength, fFloorHeight, fFloorWidth);
   move = G4ThreeVector(0., -2500., 0.);
    rotate = new G4RotationMatrix;
    direction 	  = G4ThreeVector(0., 0., 0.);
   
    

    if(fFloorLog == NULL ) {
        fFloorLog = new G4LogicalVolume(box, floorG4material, "LabFloor", 0, 0, 0);
        //fTestcanAlumCasingLog->SetVisAttributes(canVisAtt);
    }
    fAssemblyFloor->AddPlacedVolume(fFloorLog, move, rotate);


 
    return 1;
}

