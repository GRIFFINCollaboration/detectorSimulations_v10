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

double DetectionSystemSpice::SpiceResolution[2];

DetectionSystemSpice::DetectionSystemSpice() :
    // LogicalVolumes
    siInnerGuardRing_log(0),
    siOuterGuardRing_log(0)
{
    /////////////////////////////////////////////////////////////////////
    // SPICE Physical Properties
    /////////////////////////////////////////////////////////////////////

    this->wafer_material             = "Silicon";

    //-----------------------------//
    // parameters for the annular  //
    // planar detector crystal     //
    //-----------------------------//
    this->siDetCrystalOuterDiameter = 94.*mm;
    this->siDetCrystalInnerDiameter = 16.*mm;
    this->siDetCrystalThickness = 6.15*mm;
    this->siDetRadialSegments = 10.;
    this->siDetPhiSegments = 12.;

    //-----------------------------//
    // parameters for guard ring   //
    //-----------------------------//
    this->siDetGuardRingOuterDiameter = 102*mm;
    this->siDetGuardRingInnerDiameter = 10*mm;
    
}

DetectionSystemSpice::~DetectionSystemSpice()
{    // LogicalVolumes in ConstructSPICEDetectionSystem
    delete [] siDetSpiceRing_log;

    delete siInnerGuardRing_log;
    delete siOuterGuardRing_log;

}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemSpice::Build()
{

    this->assembly = new G4AssemblyVolume();

    // Loop through each ring ...
    for(int ringID=0; ringID<10; ringID++) {

        // Build assembly volumes
        this->assemblySiRing[ringID] = new G4AssemblyVolume();
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
G4int DetectionSystemSpice::PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4int ringNumber, G4int Seg, G4int SegmentNumber)
{

    //G4cout << "DetectionSystemSpice :: Ring number : "<< ringNumber <<  " SegNumber in ring : "<< Seg << " SegmentNumber in detector : " <<SegmentNumber << G4endl;
    //G4cin.get();
    G4int NumberSeg = (G4int)this->siDetPhiSegments; // total number of segments = 12
    G4double angle = (360.*deg/NumberSeg)*(Seg); // Seg = {0, ...,11}
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(-210*deg-angle); // the axis are inverted, this operation will correct for it  [MHD : 03 April 2014]

    assemblySiRing[ringNumber]->MakeImprint(exp_hall_log, move, rotate, SegmentNumber);

    return 1;
}

G4int DetectionSystemSpice::PlaceGuardRing(G4LogicalVolume* exp_hall_log, G4ThreeVector move)
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
G4int DetectionSystemSpice::BuildSiliconWafer(G4int RingID)  // RingID = { 0, 9 }
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
    G4double z_position		= -(this->siDetCrystalThickness/2.);
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    // construct solid
    G4Tubs* siDetSpiceRingSec = BuildCrystal(RingID);
    // construct logical volume
    if( !siDetSpiceRing_log[RingID] ) {
        G4String ringName = "siDetSpiceRing_";
        G4String ringID = G4UIcommand::ConvertToString(RingID);
        ringName += ringID;
        ringName += "_Log";

        //G4cout << "DetectionSystemSpice : ringName "<< ringName <<" RingID" <<  RingID<< G4endl ;
        //G4cin.get();

        siDetSpiceRing_log[RingID] = new G4LogicalVolume(siDetSpiceRingSec, material, ringName, 0, 0, 0);
        siDetSpiceRing_log[RingID]->SetVisAttributes(vis_att);
    }

    this->assemblySiRing[RingID]->AddPlacedVolume(siDetSpiceRing_log[RingID], move, rotate);


    return 1;
}

G4int DetectionSystemSpice::BuildInnerGuardRing()
{
    G4Material* material = G4Material::GetMaterial(this->wafer_material);
    if( !material ) {
        G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the inner guard ring of Spice! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    vis_att->SetVisibility(true);

    G4Tubs* innerGuardRing = new G4Tubs("innerGuardRing",
                                        this->siDetGuardRingInnerDiameter/2.,
                                        this->siDetCrystalInnerDiameter/2.,
                                        this->siDetCrystalThickness/2.,0,360);

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= -(this->siDetCrystalThickness/2.);
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    //logical volume
    if( siInnerGuardRing_log == NULL )
    {
        siInnerGuardRing_log = new G4LogicalVolume(innerGuardRing, material, "innerGuardRing", 0,0,0);
        siInnerGuardRing_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(siInnerGuardRing_log, move, rotate);

    return 1;
}

G4int DetectionSystemSpice::BuildOuterGuardRing()
{
    G4Material* material = G4Material::GetMaterial(this->wafer_material);
    if( !material ) {
        G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the outer guard ring of Spice! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    vis_att->SetVisibility(true);

    G4Tubs* outerGuardRing = new G4Tubs("outerGuardRing",
                                        this->siDetCrystalOuterDiameter/2.,
                                        this->siDetGuardRingOuterDiameter/2.,
                                        this->siDetCrystalThickness/2.,0,360);

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= -(this->siDetCrystalThickness/2.);
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    //logical volume
    if( siOuterGuardRing_log == NULL )
    {
        siOuterGuardRing_log = new G4LogicalVolume(outerGuardRing, material, "outerGuardRing", 0,0,0);
        siOuterGuardRing_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(siOuterGuardRing_log, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////
// Build one segment of Spice, 
// the geometry depends on the distance from the center
///////////////////////////////////////////////////////
G4Tubs* DetectionSystemSpice::BuildCrystal(G4int RingID)
{
    // define angle, length, thickness, and inner and outer diameter
    // of silicon detector segment
    G4double tube_element_length = (this->siDetCrystalOuterDiameter - this->siDetCrystalInnerDiameter)/(2*(this->siDetRadialSegments));
    G4double tube_element_angular_width = (360./this->siDetPhiSegments)*deg;
    G4double tube_element_inner_radius = (this->siDetCrystalInnerDiameter)/2.0 + tube_element_length*(RingID);
    G4double tube_element_outer_radius = ((G4double)this->siDetCrystalInnerDiameter)/2.0 + tube_element_length*(RingID+1);
    G4double tube_element_half_thickness = (this->siDetCrystalThickness)/2.0;

    // establish solid
    G4Tubs* crystal_block = new G4Tubs("crystal_block",tube_element_inner_radius,tube_element_outer_radius,tube_element_half_thickness,0,tube_element_angular_width);

    return crystal_block;
}//end ::BuildCrystal

G4int DetectionSystemSpice::AssignSpiceResolution(G4double intercept, G4double gradient)
{
    DetectionSystemSpice::SpiceResolution[0] = intercept;
    DetectionSystemSpice::SpiceResolution[1] = gradient;

    return 1;
} //end ::SetResolution

G4double DetectionSystemSpice::ApplySpiceResolution(G4double energy)
{
    G4double rand01 = G4UniformRand();
    G4double rand02 = G4UniformRand();
    G4double randStdNormal = sqrt(-2.0 * log(rand01)) * sin(2.0 * M_PI * rand02); //random normal(0,1)
    energy += ((DetectionSystemSpice::SpiceResolution[1] * energy) + DetectionSystemSpice::SpiceResolution[0]) * randStdNormal;

    return energy;
} //end ::ApplySpiceResolution

