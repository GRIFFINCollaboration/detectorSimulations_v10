#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

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

#include "DetectionSystemLanthanumBromide.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


DetectionSystemLanthanumBromide::DetectionSystemLanthanumBromide() :
    // LogicalVolumes
    detector_volume_log(0),
    crystal_block_log(0),
    can_cylinder_log(0),
    can_front_lid_log(0),
    can_back_lid_log(0),
    seal_front_lid_log(0),
    disc_front_lid_log(0),
    packing_cylinder_log(0),
    packing_front_lid_log(0)
{

    // NaI Crystal and Can Physical Properties
    // From: Scintillation Spectrometry gamma-ray catalogue - R. L. Heath, 1964

    /////////////////////////////////////////////////////////////////////
    //  0.019"
    //  "can" - Al Container
    //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //% 0.040"
    //% "seal" - Compressed Neoprene Rubber
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //-------------------------------------------------------------------
    //- 0.006"
    //- "disc" - Polyethylene Disc
    //-------------------------------------------------------------------
    /////////////////////////////////////////////////////////////////////
    //- 0.0625"
    //- "packing" - Packed Aluminum Oxide
    /////////////////////////////////////////////////////////////////////
    //
    //
    //  NaI Crystal - Used as first-order estimate of Lanthanum Bromide Crystal
    //
    //
    //
    //
    //


    this->detail_view_end_angle	      = 360.0*deg;
    this->crystal_material            = "Cerium_Doped_Lanthanum_Bromide";

    this->can_material                = "G4_Al";
    this->seal_material               = "G4_RUBBER_NEOPRENE";
    this->disc_material               = "G4_POLYETHYLENE";
    this->packing_material            = "G4_ALUMINUM_OXIDE";

    this->crystal_length_z            = 2.0*inchtocm*cm;
    this->crystal_inner_radius 		  = 0.0*cm;
    this->crystal_outer_radius 		  = 1.0*inchtocm*cm;

    this->packing_length_z            = this->crystal_length_z;
    this->packing_inner_radius 		  = this->crystal_outer_radius;
    this->packing_outer_radius 		  = 0.0625*inchtocm*cm + this->crystal_outer_radius;

    this->packing_lid_inner_radius    = 0.0*cm;
    this->packing_lid_outer_radius 	  = this->packing_outer_radius;
    this->packing_front_lid_thickness = 0.0625*inchtocm*cm;

    this->disc_lid_inner_radius       = 0.0*cm;
    this->disc_lid_outer_radius       = this->packing_lid_outer_radius;
    this->disc_front_lid_thickness	= 0.006*inchtocm*cm;

    this->seal_lid_inner_radius       = 0.0*cm;
    this->seal_lid_outer_radius       = this->packing_lid_outer_radius;
    this->seal_front_lid_thickness	  = 0.04*inchtocm*cm;

    this->can_length_z                = this->crystal_length_z + this->packing_front_lid_thickness + this->disc_front_lid_thickness + this->seal_front_lid_thickness;
    this->can_inner_radius            = this->packing_outer_radius;
    this->can_outer_radius            = 0.019*inchtocm*cm + this->packing_outer_radius;

    this->can_lid_inner_radius        = 0.0*cm;
    this->can_lid_outer_radius        = this->can_outer_radius;
    this->can_front_lid_thickness     = 0.019*inchtocm*cm;
    this->can_back_lid_thickness      = 0.019*inchtocm*cm;


    this->detector_length_z           = this->crystal_length_z +
            this->packing_front_lid_thickness +
            this->disc_front_lid_thickness +
            this->seal_front_lid_thickness +
            this->can_front_lid_thickness +
            this->can_back_lid_thickness;
    {
        //G4double triangleThetaAngle = (180/M_PI)*(atan((1/sqrt(3))/sqrt((11/12) + (1/sqrt(2))) )+atan((sqrt(2))/(1+sqrt(2))))*deg;
        G4double triangleThetaAngle = 54.735610317245360*deg;

        // theta
        this->detectorAngles[0][0] 	= triangleThetaAngle;
        this->detectorAngles[1][0] 	= triangleThetaAngle;
        this->detectorAngles[2][0] 	= triangleThetaAngle;
        this->detectorAngles[3][0] 	= triangleThetaAngle;
        this->detectorAngles[4][0] 	= 180.0*deg - triangleThetaAngle;
        this->detectorAngles[5][0] 	= 180.0*deg - triangleThetaAngle;
        this->detectorAngles[6][0] 	= 180.0*deg - triangleThetaAngle;
        this->detectorAngles[7][0] 	= 180.0*deg - triangleThetaAngle;
        // phi
        this->detectorAngles[0][1] 	= 22.5*deg;
        this->detectorAngles[1][1] 	= 112.5*deg;
        this->detectorAngles[2][1] 	= 202.5*deg;
        this->detectorAngles[3][1] 	= 292.5*deg;
        this->detectorAngles[4][1] 	= 22.5*deg;
        this->detectorAngles[5][1] 	= 112.5*deg;
        this->detectorAngles[6][1] 	= 202.5*deg;
        this->detectorAngles[7][1] 	= 292.5*deg;
        // yaw (alpha)
        this->detectorAngles[0][2] 	= 0.0*deg;
        this->detectorAngles[1][2] 	= 0.0*deg;
        this->detectorAngles[2][2] 	= 0.0*deg;
        this->detectorAngles[3][2] 	= 0.0*deg;
        this->detectorAngles[4][2] 	= 0.0*deg;
        this->detectorAngles[5][2] 	= 0.0*deg;
        this->detectorAngles[6][2] 	= 0.0*deg;
        this->detectorAngles[7][2] 	= 0.0*deg;
        // pitch (beta)
        this->detectorAngles[0][3] 	= triangleThetaAngle;
        this->detectorAngles[1][3] 	= triangleThetaAngle;
        this->detectorAngles[2][3] 	= triangleThetaAngle;
        this->detectorAngles[3][3] 	= triangleThetaAngle;
        this->detectorAngles[4][3] 	= 180.0*deg - triangleThetaAngle;
        this->detectorAngles[5][3] 	= 180.0*deg - triangleThetaAngle;
        this->detectorAngles[6][3] 	= 180.0*deg - triangleThetaAngle;
        this->detectorAngles[7][3] 	= 180.0*deg - triangleThetaAngle;
        // roll (gamma)
        this->detectorAngles[0][4] 	= 22.5*deg;
        this->detectorAngles[1][4] 	= 112.5*deg;
        this->detectorAngles[2][4] 	= 202.5*deg;
        this->detectorAngles[3][4] 	= 292.5*deg;
        this->detectorAngles[4][4] 	= 22.5*deg;
        this->detectorAngles[5][4] 	= 112.5*deg;
        this->detectorAngles[6][4] 	= 202.5*deg;
        this->detectorAngles[7][4] 	= 292.5*deg;
    }
}

DetectionSystemLanthanumBromide::~DetectionSystemLanthanumBromide()
{
    // LogicalVolumes
    delete detector_volume_log;
    delete crystal_block_log;
    delete can_cylinder_log;
    delete can_front_lid_log;
    delete can_back_lid_log;
    delete seal_front_lid_log;
    delete disc_front_lid_log;
    delete packing_cylinder_log;
    delete packing_front_lid_log;
}


G4int DetectionSystemLanthanumBromide::Build()//G4SDManager* mySDman)
{

    // Build assembly volume
    G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    this->assembly = myAssembly;

    G4cout << "BuildCrystalVolume" << G4endl;
    BuildCrystalVolume();
    G4cout << "BuildAluminumCanVolume" << G4endl;
    BuildAluminumCanVolume();
    G4cout << "BuildPackingVolume" << G4endl;
    BuildPackingVolume();
    G4cout << "BuildDiscVolume" << G4endl;
    BuildDiscVolume();
    G4cout << "BuildSealVolume" << G4endl;
    BuildSealVolume();

    return 1;
}

G4double DetectionSystemLanthanumBromide::GetR()
{
    // to crystal face

    G4double position = set_radial_pos - (packing_front_lid_thickness + disc_front_lid_thickness + can_front_lid_thickness - can_back_lid_thickness);
    return position;
}

G4double DetectionSystemLanthanumBromide::GetTheta(G4int i)
{
    // to crystal face

    G4double theta  = this->detectorAngles[i][0];

    return theta;
}

G4double DetectionSystemLanthanumBromide::GetPhi(G4int i)
{
    // to crystal face

    G4double phi    = this->detectorAngles[i][1];
    return phi;
}
G4double DetectionSystemLanthanumBromide::GetYaw(G4int i)
{
    // to crystal face

    G4double yaw    = this->detectorAngles[i][2];
    return yaw;
}
G4double DetectionSystemLanthanumBromide::GetPitch(G4int i)
{
    // to crystal face

    G4double pitch    = this->detectorAngles[i][3];
    return pitch;
}
G4double DetectionSystemLanthanumBromide::GetRoll(G4int i)
{
    // to crystal face

    G4double roll    = this->detectorAngles[i][4];
    return roll;
}

G4int DetectionSystemLanthanumBromide::PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector_number, G4double radialpos)
{
    G4int detector_copy_ID = 0;

    G4cout << "LanthanumBromide Detector Number = " << detector_number << G4endl;

    G4int copy_number = detector_copy_ID + detector_number;

    G4double position = radialpos + can_front_lid_thickness +( this->can_length_z / 2.0 ) ;
    set_radial_pos = radialpos;

    G4double theta  = this->detectorAngles[detector_number][0];
    G4double phi    = this->detectorAngles[detector_number][1];
    G4double alpha  = this->detectorAngles[detector_number][2]; // yaw
    G4double beta   = this->detectorAngles[detector_number][3]; // pitch
    G4double gamma  = this->detectorAngles[detector_number][4]; // roll

    G4double x = 0;
    G4double y = 0;
    G4double z = position;

    G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
    rotate->rotateY(M_PI);
    rotate->rotateY(M_PI+beta);
    rotate->rotateZ(gamma);

    G4ThreeVector move(transX(x,y,z,theta,phi), transY(x,y,z,theta,phi), transZ(x,y,z,theta,phi));

    assembly->MakeImprint(exp_hall_log, move, rotate, copy_number);

    return 1;
}

G4int DetectionSystemLanthanumBromide::BuildCrystalVolume()
{
    G4Material* material = G4Material::GetMaterial(this->crystal_material);
    if( !material ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.2,1.0,0.3));
    vis_att->SetVisibility(true);

    G4Tubs* crystal_block = BuildCrystal();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position		= ( (can_front_lid_thickness + seal_front_lid_thickness + disc_front_lid_thickness + packing_front_lid_thickness) - (can_back_lid_thickness) )/2.0;
    G4ThreeVector move 		= z_position * direction;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( crystal_block_log == NULL )
    {
        crystal_block_log = new G4LogicalVolume(crystal_block, material, "lanthanum_bromide_crystal_block_log", 0, 0, 0);
        crystal_block_log->SetVisAttributes(vis_att);
    }

    this->assembly->AddPlacedVolume(crystal_block_log, move, rotate);

    return 1;
}

G4int DetectionSystemLanthanumBromide::BuildAluminumCanVolume()
{
    G4Material* material = G4Material::GetMaterial(this->can_material);
    if( !material ) {
        G4cout << " ----> Material " << this->can_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.55, 0.38, 1.0));
    vis_att->SetVisibility(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double z_position;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Place Aluminum Can
    /////////////////////////////////////////////////////////////////////
    G4Tubs* can_cylinder = BuildAluminumCan();

    //logical volume
    if( can_cylinder_log == NULL )
    {
        can_cylinder_log = new G4LogicalVolume(can_cylinder, material, "can_cylinder_log", 0, 0, 0);
        can_cylinder_log->SetVisAttributes(vis_att);
    }

    // place front can_lid
    z_position 	= 0;
    move 		= z_position * direction;

    //add physical cylinder
    this->assembly->AddPlacedVolume(can_cylinder_log, move, rotate);

    /////////////////////////////////////////////////////////////////////
    // Build and Place Aluminum Front Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* can_front_lid = BuildAluminumCanFrontLid();

    // logical volume
    if( can_front_lid_log == NULL )
    {
        can_front_lid_log = new G4LogicalVolume(can_front_lid, material, "can_front_lid_log", 0, 0, 0);
        can_front_lid_log->SetVisAttributes(vis_att);
    }

    // place front can_lid

    z_position 	= -(detector_length_z/2.0) + (can_front_lid_thickness/2.0);
    move 		= z_position * direction;

    //add physical front can_lid
    this->assembly->AddPlacedVolume(can_front_lid_log, move, rotate);

    /////////////////////////////////////////////////////////////////////
    // Build and Place Aluminum Back Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* can_back_lid = BuildAluminumCanBackLid();

    // logical volume
    if( can_back_lid_log == NULL )
    {
        can_back_lid_log = new G4LogicalVolume(can_back_lid, material, "can_back_lid_log", 0, 0, 0);
        can_back_lid_log->SetVisAttributes(vis_att);
    }

    // place back can_lid
    z_position 	= (detector_length_z/2.0) - (can_back_lid_thickness/2.0);
    move 		= z_position * direction;

    // add physical back can_lid
    this->assembly->AddPlacedVolume(can_back_lid_log, move, rotate);

    return 1;
}

G4int DetectionSystemLanthanumBromide::BuildPackingVolume()
{
    G4Material* material = G4Material::GetMaterial(this->packing_material);
    if( !material ) {
        G4cout << " ----> Material " << this->packing_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    vis_att->SetVisibility(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double z_position;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Place Packing
    /////////////////////////////////////////////////////////////////////
    G4Tubs* packing_cylinder = BuildPacking();

    //logical volume
    if( packing_cylinder_log == NULL )
    {
        packing_cylinder_log = new G4LogicalVolume(packing_cylinder, material, "packing_cylinder_log", 0, 0, 0);
        packing_cylinder_log->SetVisAttributes(vis_att);
    }

    // place cylinder
    z_position 	= ( (can_front_lid_thickness + seal_front_lid_thickness + disc_front_lid_thickness + packing_front_lid_thickness) - (can_back_lid_thickness) )/2.0;
    move          = z_position * direction;

    //add physical cylinder
    this->assembly->AddPlacedVolume(packing_cylinder_log, move, rotate);

    /////////////////////////////////////////////////////////////////////
    // Build and Place Front Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* packing_front_lid = BuildPackingFrontLid();

    // logical volume
    if( packing_front_lid_log == NULL )
    {
        packing_front_lid_log = new G4LogicalVolume(packing_front_lid, material, "packing_front_lid_log", 0, 0, 0);
        packing_front_lid_log->SetVisAttributes(vis_att);
    }

    // place front packing_lid
    z_position 	= -(detector_length_z/2.0) + (can_front_lid_thickness + seal_front_lid_thickness + disc_front_lid_thickness + (packing_front_lid_thickness/2.0));
    move 		= z_position * direction;

    //add physical front packing_lid
    this->assembly->AddPlacedVolume(packing_front_lid_log, move, rotate);

    return 1;
}

G4int DetectionSystemLanthanumBromide::BuildDiscVolume()
{
    G4Material* material = G4Material::GetMaterial(this->disc_material);
    if( !material ) {
        G4cout << " ----> Material " << this->disc_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    vis_att->SetVisibility(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double z_position;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Place Front Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* disc_front_lid = BuildDiscFrontLid();

    // logical volume
    if( disc_front_lid_log == NULL )
    {
        disc_front_lid_log = new G4LogicalVolume(disc_front_lid, material, "disc_front_lid_log", 0, 0, 0);
        disc_front_lid_log->SetVisAttributes(vis_att);
    }

    // place front disc_lid
    z_position 	= -(detector_length_z/2.0) + (can_front_lid_thickness + seal_front_lid_thickness + (disc_front_lid_thickness/2.0));
    move 		= z_position * direction;

    //add physical front disc_lid
    this->assembly->AddPlacedVolume(disc_front_lid_log, move, rotate);

    return 1;
}

G4int DetectionSystemLanthanumBromide::BuildSealVolume()
{
    G4Material* material = G4Material::GetMaterial(this->seal_material);
    if( !material ) {
        G4cout << " ----> Material " << this->seal_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
    vis_att->SetVisibility(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double z_position;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Place Front Lid
    /////////////////////////////////////////////////////////////////////
    G4Tubs* seal_front_lid = BuildSealFrontLid();

    // logical volume
    if( seal_front_lid_log == NULL )
    {
        seal_front_lid_log = new G4LogicalVolume(seal_front_lid, material, "seal_front_lid_log", 0, 0, 0);
        seal_front_lid_log->SetVisAttributes(vis_att);
    }

    // place front seal_lid
    z_position 	= -(detector_length_z/2.0) + (can_front_lid_thickness + (seal_front_lid_thickness/2.0));
    move 		= z_position * direction;

    //add physical front seal_lid
    this->assembly->AddPlacedVolume(seal_front_lid_log, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Tubs* DetectionSystemLanthanumBromide::BuildCrystal()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = crystal_inner_radius;
    G4double outer_radius = crystal_outer_radius;
    G4double half_length_z = (crystal_length_z)/2.0;

    G4Tubs* crystal_block = new G4Tubs("crystal_block", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return crystal_block;
}//end ::BuildCrystal


G4Tubs* DetectionSystemLanthanumBromide::BuildAluminumCan()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius 	= can_inner_radius;
    G4double outer_radius 	= can_outer_radius;
    G4double half_length_z 	= can_length_z/2.0;

    G4Tubs* can_cylinder = new G4Tubs("can_cylinder", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return can_cylinder;
}//end ::BuildAluminumCan

G4Tubs* DetectionSystemLanthanumBromide::BuildAluminumCanFrontLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = can_lid_inner_radius;
    G4double outer_radius = can_lid_outer_radius;
    G4double half_length_z = can_front_lid_thickness/2.0;

    G4Tubs* can_lid = new G4Tubs("can_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return can_lid;
}//end ::BuildAluminumFrontLid

G4Tubs* DetectionSystemLanthanumBromide::BuildAluminumCanBackLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = can_lid_inner_radius;
    G4double outer_radius = can_lid_outer_radius;
    G4double half_length_z = can_back_lid_thickness/2.0;

    G4Tubs* can_lid = new G4Tubs("can_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return can_lid;
}//end ::BuildAluminumBackLid

G4Tubs* DetectionSystemLanthanumBromide::BuildPacking()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius 	= packing_inner_radius;
    G4double outer_radius 	= packing_outer_radius;
    G4double half_length_z 	= packing_length_z/2.0;

    G4Tubs* packing_cylinder = new G4Tubs("packing_cylinder", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return packing_cylinder;
}//end ::BuildPacking

G4Tubs* DetectionSystemLanthanumBromide::BuildPackingFrontLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = packing_lid_inner_radius;
    G4double outer_radius = packing_lid_outer_radius;
    G4double half_length_z = packing_front_lid_thickness/2.0;

    G4Tubs* packing_lid = new G4Tubs("packing_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return packing_lid;
}//end ::BuildPackingFrontLid

G4Tubs* DetectionSystemLanthanumBromide::BuildDiscFrontLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = disc_lid_inner_radius;
    G4double outer_radius = disc_lid_outer_radius;
    G4double half_length_z = disc_front_lid_thickness/2.0;

    G4Tubs* disc_lid = new G4Tubs("disc_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return disc_lid;
}//end ::BuildDiscFrontLid

G4Tubs* DetectionSystemLanthanumBromide::BuildSealFrontLid()
{
    G4double start_phi = 0.0;
    G4double end_phi = this->detail_view_end_angle;

    G4double inner_radius = seal_lid_inner_radius;
    G4double outer_radius = seal_lid_outer_radius;
    G4double half_length_z = seal_front_lid_thickness/2.0;

    G4Tubs* seal_lid = new G4Tubs("seal_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return seal_lid;
}//end ::BuildSealFrontLid

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemLanthanumBromide::GetDirectionXYZ(G4double theta, G4double phi)
{
    G4double x,y,z;
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);

    G4ThreeVector direction = G4ThreeVector(x,y,z);

    return direction;
}//end ::GetDirection

G4double DetectionSystemLanthanumBromide::transX(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*cos(phi) );
}

G4double DetectionSystemLanthanumBromide::transY(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*sin(phi) );
}

G4double DetectionSystemLanthanumBromide::transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*cos(theta) );
}
