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

#include "DetectionSystemSpiceV02.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

#define NUMBEROFCELLS 4

DetectionSystemSpiceV02::DetectionSystemSpiceV02() :
    // LogicalVolumes
    fDetectorCasingSideLog(0),
    fDetectorCasingBackLog(0),
    fCrystalBlockLog(0)
{
    /////////////////////////////////////////////////////////////////////
    // SPICE Physical Properties
    /////////////////////////////////////////////////////////////////////

    fCasingMaterial            = "Aluminum";
    fWaferMaterial             = "Silicon";

    //-----------------------------//
    // parameters for the square   //
    // planar detector crystal     //
    //-----------------------------//
    fSquareDetCrystalLength = 38*mm;
    fSquareDet90DegCrystalRadius = 63.64*mm;
    fSquareDet45DegCrystalRadius = 45*mm;
    fSquareDetCrystalThickness = 5*mm;

    //-----------------------------//
    // parameters for the square   //
    // planar detector casing      //
    //-----------------------------//
    fSquareDetCasingLength = 42*mm;
    fSquareDetCasingThickness = 14*mm;

    // minimial distance 122.2mm, assuming 2mm gap between casings
    // calculated by distance of back plate - detector casing thickness (z)
    // - 150.65 (back face)
    // + 20 (2 times mount cylinders)
    // + 14 (casing size z)
    fDetectorFaceTargetDistance              = -116.65*mm;

    //myDataOutput = dataOutput;

    // These keep track of the number of sensitive detectors. Initialized to zero
    //siDetCopyNumber=SIDETCOPYNUMBER;
    //siDetCasingSideCopyNumber=SIDETCASINGSIDECOPYNUMBER;
    //siDetCasingBackCopyNumber=SIDETCASINGBACKCOPYNUMBER;

}

DetectionSystemSpiceV02::~DetectionSystemSpiceV02() {
    // LogicalVolumes in ConstructSPICEDetectionSystem
    delete fCrystalBlockLog;
    delete fDetectorCasingSideLog;
    delete fDetectorCasingBackLog;

}

G4int DetectionSystemSpiceV02::Build()  {

    // Build assembly volume
    fAssembly = new G4AssemblyVolume();
    fAssemblySi = new G4AssemblyVolume();

    G4cout << "BuildSiliconWafer" << G4endl;
    BuildSiliconWafer();
    G4cout << "BuildAluminiumCasingSide" << G4endl;
    BuildAluminiumCasingSide();
    G4cout << "BuildAluminiumCasingBack" << G4endl;
    BuildAluminiumCasingBack();

    return 1;
}

G4int DetectionSystemSpiceV02::PlaceDetector(G4LogicalVolume* expHallLog, G4ThreeVector move, G4RotationMatrix* rotate, G4int detectorNumber) {
    G4int detectorCopyID = 0;

    G4cout << "SpiceV02 Detector Number = " << detectorNumber << G4endl;

    G4int copyNumber = detectorCopyID + detectorNumber;

    fAssemblySi->MakeImprint(expHallLog, move, rotate, copyNumber);
    fAssembly->MakeImprint(expHallLog, move, rotate, copyNumber);

    return 1;
}

G4int DetectionSystemSpiceV02::BuildSiliconWafer() {
    G4Material* material = G4Material::GetMaterial(fWaferMaterial);
    if( !material ) {
        G4cout << " ----> Material " << fWaferMaterial << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    visAtt->SetVisibility(true);

    G4Box* crystalBlock = BuildCrystal();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double zPosition		= 0.0*mm;
    G4ThreeVector move 		= zPosition * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if(fCrystalBlockLog == NULL) {
        fCrystalBlockLog = new G4LogicalVolume(crystalBlock, material, "crystalBlockLog", 0, 0, 0);
        fCrystalBlockLog->SetVisAttributes(visAtt);
    }

    fAssemblySi->AddPlacedVolume(fCrystalBlockLog, move, rotate);

    return 1;
}

G4int DetectionSystemSpiceV02::BuildAluminiumCasingSide() {
    G4Material* material = G4Material::GetMaterial(fCasingMaterial);
    if( !material ) {
        G4cout << " ----> Material " << fCasingMaterial << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(ALCOL));
    visAtt->SetVisibility(true);

    G4SubtractionSolid* detectorCasingSide = BuildCasingSide();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double zPosition		= 0.0*mm;
    G4ThreeVector move 		= zPosition * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if(fDetectorCasingSideLog == NULL) {
        fDetectorCasingSideLog = new G4LogicalVolume(detectorCasingSide, material, "detectorCasingLog", 0,0,0);
        fDetectorCasingSideLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fDetectorCasingSideLog, move, rotate);

    return 1;
}

G4int DetectionSystemSpiceV02::BuildAluminiumCasingBack() {
    G4Material* material = G4Material::GetMaterial(fCasingMaterial);
    if(!material) {
        G4cout << " ----> Material " << fCasingMaterial << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(ALCOL));
    visAtt->SetVisibility(true);

    G4Box* detectorCasingBack = BuildCasingBack();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double zPosition		= 0.0*mm;
    G4ThreeVector move 		= zPosition * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if(fDetectorCasingBackLog == NULL) {
        fDetectorCasingBackLog = new G4LogicalVolume(detectorCasingBack, material,"detectorCasingBackLog", 0,0,0);
        fDetectorCasingBackLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fDetectorCasingBackLog, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemSpiceV02::BuildCrystal() {
    // silicon detector dimensions
    G4double halfLengthX = (fSquareDetCrystalLength)/2.0;
    G4double halfLengthY = (fSquareDetCrystalLength)/2.0;
    G4double halfLengthZ = (fSquareDetCrystalThickness)/2.0;

    // establish solid
    G4Box* crystalBlock = new G4Box("crystalBlock", halfLengthX,halfLengthY,halfLengthZ);

    return crystalBlock;
}//end ::BuildCrystal

G4SubtractionSolid* DetectionSystemSpiceV02::BuildCasingSide() {
    // dimensions of casing
    G4double halfLengthX = (fSquareDetCasingLength)/2.0;
    G4double halfLengthY = (fSquareDetCasingLength)/2.0;
    G4double halfLengthZ = (fSquareDetCrystalThickness)/2.0;
    // dimensions of silicon detector
    G4double crysHalfLengthZ = (fSquareDetCrystalThickness)/2.0;

    // relative distance
    G4double zTransl = (halfLengthZ - crysHalfLengthZ);
    G4ThreeVector transl(0,0,zTransl);

    // establish solid
    G4Box* filledBox = new G4Box("filledBox",halfLengthX,halfLengthY,halfLengthZ);
    G4Box* crystalBlock = BuildCrystal();
    G4SubtractionSolid* detectorCasingSide = new G4SubtractionSolid("detectorCasingSide",filledBox,crystalBlock,0,transl);

    return detectorCasingSide;
}//end ::BuildCasingSide

G4Box* DetectionSystemSpiceV02::BuildCasingBack() {
    // dimensions of casing
    G4double halfLengthX = (fSquareDetCasingLength)/2.0;
    G4double halfLengthY = (fSquareDetCasingLength)/2.0;
    G4double halfLengthZ = (fSquareDetCasingThickness-fSquareDetCrystalThickness)/2.0;

    // establish solid
    G4Box* detectorCasingBack = new G4Box("detectorCasingBack", halfLengthX,halfLengthY,halfLengthZ);

    return detectorCasingBack;
}//end ::BuildCasingBack

