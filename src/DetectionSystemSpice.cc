#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
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

#include "Randomize.hh"

#include "DetectionSystemSpice.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

double DetectionSystemSpice::fSpiceResolution[2];

DetectionSystemSpice::DetectionSystemSpice() :
    // LogicalVolumes
    fSiInnerGuardRingLog(0),
    fSiOuterGuardRingLog(0)
{
    /////////////////////////////////////////////////////////////////////
    // SPICE Physical Properties
    /////////////////////////////////////////////////////////////////////

    fWaferMaterial             = "Silicon";

    //-----------------------------//
    // parameters for the annular  //
    // planar detector crystal     //
    //-----------------------------//
    fSiDetCrystalOuterDiameter = 94.*mm;
    fSiDetCrystalInnerDiameter = 16.*mm;
    fSiDetCrystalThickness = 6.15*mm;
    fSiDetRadialSegments = 10.;
    fSiDetPhiSegments = 12.;

    //-----------------------------//
    // parameters for guard ring   //
    //-----------------------------//
    fSiDetGuardRingOuterDiameter = 102*mm;
    fSiDetGuardRingInnerDiameter = 10*mm;
    
}

DetectionSystemSpice::~DetectionSystemSpice() {    // LogicalVolumes in ConstructSPICEDetectionSystem
  for(int i = 0; i < 10; ++i) {
    if(fSiDetSpiceRingLog[i] != NULL) {
      delete fSiDetSpiceRingLog[i];
    }
  }

    delete fSiInnerGuardRingLog;
    delete fSiOuterGuardRingLog;

}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemSpice::Build() {

    fAssembly = new G4AssemblyVolume();

    // Loop through each ring ...
    for(int ringID=0; ringID<10; ringID++) {

        // Build assembly volumes
        fAssemblySiRing[ringID] = new G4AssemblyVolume();
        // Build the logical segment of the Silicon Ring
        BuildSiliconWafer(ringID);

    } // end for(int ringID)

    BuildInnerGuardRing();
    BuildOuterGuardRing();

    return 1;
} // end Build

//---------------------------------------------------------//
// "place" function called in DetectorMessenger            //
// if detector is added                                    //
//---------------------------------------------------------//
G4int DetectionSystemSpice::PlaceDetector(G4LogicalVolume* expHallLog, G4ThreeVector move, G4int ringNumber, G4int Seg, G4int SegmentNumber) {

    //G4cout << "DetectionSystemSpice :: Ring number : "<< ringNumber <<  " SegNumber in ring : "<< Seg << " SegmentNumber in detector : " <<SegmentNumber << G4endl;
    //G4cin.get();
    G4int NumberSeg = (G4int)fSiDetPhiSegments; // total number of segments = 12
    G4double angle = (360.*deg/NumberSeg)*(Seg); // Seg = {0, ...,11}
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(-210*deg-angle); // the axis are inverted, this operation will correct for it  [MHD : 03 April 2014]

    fAssemblySiRing[ringNumber]->MakeImprint(expHallLog, move, rotate, SegmentNumber);

    return 1;
}

G4int DetectionSystemSpice::PlaceGuardRing(G4LogicalVolume* expHallLog, G4ThreeVector move) {
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);
    fAssembly->MakeImprint(expHallLog, move, rotate, 0);

    return 1;
}

//---------------------------------------------------------//
// build functions for different parts                     //
// called in main build function                           //
//---------------------------------------------------------//
G4int DetectionSystemSpice::BuildSiliconWafer(G4int RingID)  {// RingID = { 0, 9 }
    // Define the material, return error if not found
    G4Material* material = G4Material::GetMaterial(fWaferMaterial);
    if( !material ) {
        G4cout << " ----> Material " << fWaferMaterial
               << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    visAtt->SetVisibility(true);

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double zPosition		= -(fSiDetCrystalThickness/2.);
    G4ThreeVector move 		= zPosition * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    // construct solid
    G4Tubs* siDetSpiceRingSec = BuildCrystal(RingID);
    // construct logical volume
    if(!fSiDetSpiceRingLog[RingID]) {
        G4String ringName = "siDetSpiceRing_";
        G4String ringID = G4UIcommand::ConvertToString(RingID);
        ringName += ringID;
        ringName += "_Log";

        //G4cout << "DetectionSystemSpice : ringName "<< ringName <<" RingID" <<  RingID<< G4endl ;
        //G4cin.get();

        fSiDetSpiceRingLog[RingID] = new G4LogicalVolume(siDetSpiceRingSec, material, ringName, 0, 0, 0);
        fSiDetSpiceRingLog[RingID]->SetVisAttributes(visAtt);
    }

    fAssemblySiRing[RingID]->AddPlacedVolume(fSiDetSpiceRingLog[RingID], move, rotate);

    return 1;
}

G4int DetectionSystemSpice::BuildInnerGuardRing() {
    G4Material* material = G4Material::GetMaterial(fWaferMaterial);
    if( !material ) {
        G4cout << " ----> Material " << fWaferMaterial << " not found, cannot build the inner guard ring of Spice! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    visAtt->SetVisibility(true);

    G4Tubs* innerGuardRing = new G4Tubs("innerGuardRing",
                                        fSiDetGuardRingInnerDiameter/2.,
                                        fSiDetCrystalInnerDiameter/2.,
                                        fSiDetCrystalThickness/2.,0,360);

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double zPosition		= -(fSiDetCrystalThickness/2.);
    G4ThreeVector move 		= zPosition * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    //logical volume
    if(fSiInnerGuardRingLog == NULL) {
        fSiInnerGuardRingLog = new G4LogicalVolume(innerGuardRing, material, "innerGuardRing", 0,0,0);
        fSiInnerGuardRingLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fSiInnerGuardRingLog, move, rotate);

    return 1;
}

G4int DetectionSystemSpice::BuildOuterGuardRing() {
    G4Material* material = G4Material::GetMaterial(fWaferMaterial);
    if( !material ) {
        G4cout << " ----> Material " << fWaferMaterial << " not found, cannot build the outer guard ring of Spice! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    visAtt->SetVisibility(true);

    G4Tubs* outerGuardRing = new G4Tubs("outerGuardRing",
                                        fSiDetCrystalOuterDiameter/2.,
                                        fSiDetGuardRingOuterDiameter/2.,
                                        fSiDetCrystalThickness/2.,0,360);

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double zPosition		= -(fSiDetCrystalThickness/2.);
    G4ThreeVector move 		= zPosition * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    //logical volume
    if(fSiOuterGuardRingLog == NULL) {
        fSiOuterGuardRingLog = new G4LogicalVolume(outerGuardRing, material, "outerGuardRing", 0,0,0);
        fSiOuterGuardRingLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fSiOuterGuardRingLog, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////
// Build one segment of Spice, 
// the geometry depends on the distance from the center
///////////////////////////////////////////////////////
G4Tubs* DetectionSystemSpice::BuildCrystal(G4int RingID) {
    // define angle, length, thickness, and inner and outer diameter
    // of silicon detector segment
    G4double tubeElementLength = (fSiDetCrystalOuterDiameter - fSiDetCrystalInnerDiameter)/(2*(fSiDetRadialSegments));
    G4double tubeElementAngularWidth = (360./fSiDetPhiSegments)*deg;
    G4double tubeElementInnerRadius = (fSiDetCrystalInnerDiameter)/2.0 + tubeElementLength*(RingID);
    G4double tubeElementOuterRadius = ((G4double)fSiDetCrystalInnerDiameter)/2.0 + tubeElementLength*(RingID+1);
    G4double tubeElementHalfThickness = (fSiDetCrystalThickness)/2.0;

    // establish solid
    G4Tubs* crystalBlock = new G4Tubs("crystalBlock",tubeElementInnerRadius,tubeElementOuterRadius,tubeElementHalfThickness,0,tubeElementAngularWidth);

    return crystalBlock;
}//end ::BuildCrystal

G4int DetectionSystemSpice::AssignSpiceResolution(G4double intercept, G4double gradient) {
    fSpiceResolution[0] = intercept;
    fSpiceResolution[1] = gradient;

    return 1;
} //end ::SetResolution

G4double DetectionSystemSpice::ApplySpiceResolution(G4double energy) {
    G4double rand01 = G4UniformRand();
    G4double rand02 = G4UniformRand();
    G4double randStdNormal = sqrt(-2.0 * log(rand01)) * sin(2.0 * M_PI * rand02); //random normal(0,1)
    energy += ((fSpiceResolution[1] * energy) + fSpiceResolution[0]) * randStdNormal;

    return energy;
} //end ::ApplySpiceResolution

