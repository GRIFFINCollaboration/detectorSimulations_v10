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

#include "DetectionSystemS3.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

DetectionSystemS3::DetectionSystemS3() :
    // LogicalVolumes
    S3InnerGuardRing_log(0),
    S3OuterGuardRing_log(0)
{
    /////////////////////////////////////////////////////////////////////
    // SPICE Physical Properties
    /////////////////////////////////////////////////////////////////////

    this->wafer_material             = "Silicon";

    //-----------------------------//
    // parameters for the annular  //
    // planar detector crystal     //
    //-----------------------------//
    this->S3DetCrystalOuterDiameter = 70.*mm;
    this->S3DetCrystalInnerDiameter = 22.*mm;
    this->S3DetCrystalThickness = .15*mm;
    this->S3DetRadialSegments = 24.;
    this->S3DetPhiSegments = 32.;

    //-----------------------------//
    // parameters for guard ring   //
    //-----------------------------//
    this->S3DetGuardRingOuterDiameter = 76*mm;
    this->S3DetGuardRingInnerDiameter = 20*mm;
}

DetectionSystemS3::~DetectionSystemS3()
{   
    delete [] siDetS3Ring_log;
    delete S3InnerGuardRing_log;
    delete S3OuterGuardRing_log;

}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemS3::Build()
{
    this->assembly =  new G4AssemblyVolume();

    for (int ringID=0; ringID<24; ringID++) {	// Loop through each ring ...

        this->assemblyS3Ring[ringID] = new G4AssemblyVolume(); 		// Build assembly volumes
        BuildSiliconWafer(ringID);		// Build Silicon Ring
    } // end for(int ringID)

    // Build Guard Rings
    BuildInnerGuardRing();
    BuildOuterGuardRing();

    return 1;
}

//---------------------------------------------------------//
// "place" function called in DetectorMessenger            //
// if detector is added                                    //
//---------------------------------------------------------//
G4int DetectionSystemS3::PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4int ringNumber, G4int Seg, G4int detectorNumber)
{
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4int NumberSeg = (G4int)this->S3DetPhiSegments;
    G4double angle = (360./NumberSeg)*(Seg-0.5)*deg;
    rotate->rotateZ(angle);
    assemblyS3Ring[ringNumber]->MakeImprint(exp_hall_log, move, rotate, detectorNumber);

    return 1;
}

G4int DetectionSystemS3::PlaceGuardRing(G4LogicalVolume* exp_hall_log, G4ThreeVector move)
{
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);
    assembly->MakeImprint(exp_hall_log, move, rotate, 0);

    return 1;
}

//---------------------------------------------------------//
// build functions for different parts                     //
// called in main build function                           //
//---------------------------------------------------------//
G4int DetectionSystemS3::BuildSiliconWafer(G4int RingID)
{
    // Define the material, return error if not found
    G4Material* material = G4Material::GetMaterial(this->wafer_material);
    if( !material ) {
        G4cout << " ----> Material " << this->wafer_material
               << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    vis_att->SetVisibility(true);

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= -(this->S3DetCrystalThickness/2.);
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    // construct solid
    G4Tubs* siDetS3RingSec = BuildCrystal(RingID);

    // construct logical volume if it doesn't already exist
    if( !siDetS3Ring_log[RingID] )
    {
        G4String s3name = "siDetS3Ring_";
        s3name += G4UIcommand::ConvertToString(RingID);
        s3name += "_Log";

        siDetS3Ring_log[RingID] = new G4LogicalVolume(siDetS3RingSec, material, s3name, 0, 0, 0);
        siDetS3Ring_log[RingID]->SetVisAttributes(vis_att);
    }
    this->assemblyS3Ring[RingID]->AddPlacedVolume(siDetS3Ring_log[RingID], move, rotate);
    
    return 1;
}

G4int DetectionSystemS3::BuildInnerGuardRing()
{
    G4Material* material = G4Material::GetMaterial(this->wafer_material);
    if( !material ) {
        G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the inner guard ring of the S3 detector! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    vis_att->SetVisibility(true);

    G4Tubs* innerGuardRing = new G4Tubs("innerGuardRing",
                                        this->S3DetGuardRingInnerDiameter/2.,
                                        this->S3DetCrystalInnerDiameter/2.,
                                        this->S3DetCrystalThickness/2.,0,360);

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= -(this->S3DetCrystalThickness/2.);
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    //logical volume
    if( S3InnerGuardRing_log == NULL )
    {
        S3InnerGuardRing_log = new G4LogicalVolume(innerGuardRing, material, "innerGuardRing", 0,0,0);
        S3InnerGuardRing_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(S3InnerGuardRing_log, move, rotate);

    return 1;
}

G4int DetectionSystemS3::BuildOuterGuardRing()
{
    G4Material* material = G4Material::GetMaterial(this->wafer_material);
    if( !material ) {
        G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the outer guard ring of the S3 detector! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    vis_att->SetVisibility(true);

    G4Tubs* outerGuardRing = new G4Tubs("outerGuardRing",
                                        this->S3DetCrystalOuterDiameter/2.,
                                        this->S3DetGuardRingOuterDiameter/2.,
                                        this->S3DetCrystalThickness/2.,0,360);

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= -(this->S3DetCrystalThickness/2.);
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    //logical volume
    if( S3OuterGuardRing_log == NULL )
    {
        S3OuterGuardRing_log = new G4LogicalVolume(outerGuardRing, material, "outerGuardRing", 0,0,0);
        S3OuterGuardRing_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(S3OuterGuardRing_log, move, rotate);

    return 1;
}


///////////////////////////////////////////////////////
// Build one segment of S3
// the geometry depends on the distance from the center
///////////////////////////////////////////////////////
G4Tubs* DetectionSystemS3::BuildCrystal(G4int RingID)
{
    // define angle, length, thickness, and inner and outer diameter
    // of silicon detector segment
    G4double tube_element_length = (this->S3DetCrystalOuterDiameter - this->S3DetCrystalInnerDiameter)/(2*(this->S3DetRadialSegments));
    G4double tube_element_angular_width = (360./this->S3DetPhiSegments)*deg;
    G4double tube_element_inner_radius = (this->S3DetCrystalInnerDiameter)/2.0 + tube_element_length*(RingID);
    G4double tube_element_outer_radius = ((G4double)this->S3DetCrystalInnerDiameter)/2.0 + tube_element_length*(RingID+1);
    G4double tube_element_half_thickness = (this->S3DetCrystalThickness)/2.0;

    // establish solid
    G4Tubs* crystal_block = new G4Tubs("crystal_block",tube_element_inner_radius,tube_element_outer_radius,tube_element_half_thickness,0,tube_element_angular_width);

    return crystal_block;
}//end ::BuildCrystal

