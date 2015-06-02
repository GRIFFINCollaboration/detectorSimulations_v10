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

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "ApparatusGriffinStructure.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

ApparatusGriffinStructure::ApparatusGriffinStructure() :
    // LogicalVolumes
    square_piece_log(0),
    triangle_piece_log(0),
    ring_piece_log(0),
    rod_piece_log(0)
{
    this->structure_material = "G4_Al";
    this->rod_material = "G4_STAINLESS-STEEL";
    this->inner_distance_to_square_face = 470.0*mm; // inner radius to the centre of square plate - GOOD
    this->square_and_triangle_face_thickness = 95.0*mm; // thickness of the square and triangle plates. In reality these thicknesses are different, but we will approximate them to be the same. - GOOD.
    this->square_hole_length = 314.0*mm; // total legth of square hole, from inside structure. - GOOD.
    this->square_hole_rounded_edge_radius = 20.0*mm; // the roundness of the square holes. This is approximated, won't have a large effect. - OKAY
    this->triangle_outer_hole_radius = 120.65*mm; // a hole radius in the triangle face. Radius to the rounded edge, this is the larger radius. - GOOD.
    this->triangle_inner_hole_radius = 101.7*mm; // a hole radius in the triangle face. Radius to the square edge, this is the smaller radius. - GOOD.
    this->square_rod_hole_radius = (31.75*mm)/2.0; // the radius of the holes for the steel rods. - GOOD.
    this->square_rod_hole_inner_radius = 568.14*mm; // the radius from the origin to the steel rod. - GOOD.
    this->square_rod_position_from_center_of_square_face = (350.0*mm)/2.0; // the position the holes are from the origin. This is NOT a radius, this is along x or y. - GOOD.
    this->square_rod_total_length = 570.0*mm; // the total length of the rods. - GOOD.
    this->square_bgo_cut_length = (118.5*mm)*2.0; // the length of the bgo slot. - GOOD.
    this->square_bgo_cut_width = 30.0*mm; // the width of the bgo slot. - GOOD.
    // ref http://en.wikipedia.org/wiki/Rhombicuboctahedron
    this->square_square_angle = 144.74*deg; // the angle between square faces in a rhombicuboctahedron. - GOOD.
    this->square_triangle_angle = 135.0*deg; // the angle between a square face and a triangle face in a rhombicuboctahedron. - GOOD.
    this->large_ring_frame_inner_radius = 850.0*mm; // the inner radius of one of the the large rings that holds up the structure. - GOOD.
    this->large_ring_frame_outer_radius = 1100.0*mm; // the outer radius of one of the the large rings that holds up the structure. - GOOD.
    this->large_ring_frame_thickness = 50.8*mm; // the thickness of one of the the large rings that holds up the structure. - GOOD.
    this->large_ring_frame_z_position = 240.0*mm; // the positive z distance the ring is from the origin. - GOOD.
    //G4double triangleThetaAngle = (180/M_PI)*(atan((1/sqrt(3))/sqrt((11/12) + (1/sqrt(2))) )+atan((sqrt(2))/(1+sqrt(2))))*deg;
    this->triangleThetaAngle = 54.735610317245360*deg;

    // Check surfaces to determine any problematic overlaps. Turn this on to have Geant4 check the surfaces.
    // Do not leave this on, it will slow the DetectorConstruction process!
    // This was last check on May 29th, 2015. - GOOD!
    surfCheck = false;

    /////////////////////////////////////////////////////////////////////
    // griffinCoords for GRIFFIN
    // Note that the GRIFFIN lampshade angles are rotated by 45 degrees with respect to those of TIGRESS.
    // Modified griffinCoords for TIGRESS are below!
    /////////////////////////////////////////////////////////////////////
    // theta
    this->griffinCoords[0][0] 	= 45.0*deg;
    this->griffinCoords[1][0] 	= 45.0*deg;
    this->griffinCoords[2][0] 	= 45.0*deg;
    this->griffinCoords[3][0] 	= 45.0*deg;
    this->griffinCoords[4][0] 	= 90.0*deg;
    this->griffinCoords[5][0] 	= 90.0*deg;
    this->griffinCoords[6][0] 	= 90.0*deg;
    this->griffinCoords[7][0] 	= 90.0*deg;
    this->griffinCoords[8][0] 	= 90.0*deg;
    this->griffinCoords[9][0] 	= 90.0*deg;
    this->griffinCoords[10][0] 	= 90.0*deg;
    this->griffinCoords[11][0] 	= 90.0*deg;
    this->griffinCoords[12][0] 	= 135.0*deg;
    this->griffinCoords[13][0] 	= 135.0*deg;
    this->griffinCoords[14][0] 	= 135.0*deg;
    this->griffinCoords[15][0] 	= 135.0*deg;
    // phi
    this->griffinCoords[0][1] 	= 67.5*deg;
    this->griffinCoords[1][1] 	= 157.5*deg;
    this->griffinCoords[2][1] 	= 247.5*deg;
    this->griffinCoords[3][1] 	= 337.5*deg;
    this->griffinCoords[4][1] 	= 22.5*deg;
    this->griffinCoords[5][1] 	= 67.5*deg;
    this->griffinCoords[6][1] 	= 112.5*deg;
    this->griffinCoords[7][1] 	= 157.5*deg;
    this->griffinCoords[8][1] 	= 202.5*deg;
    this->griffinCoords[9][1] 	= 247.5*deg;
    this->griffinCoords[10][1] 	= 292.5*deg;
    this->griffinCoords[11][1] 	= 337.5*deg;
    this->griffinCoords[12][1] 	= 67.5*deg;
    this->griffinCoords[13][1] 	= 157.5*deg;
    this->griffinCoords[14][1] 	= 247.5*deg;
    this->griffinCoords[15][1] 	= 337.5*deg;
    // yaw (alpha)
    this->griffinCoords[0][2] 	= 0.0*deg;
    this->griffinCoords[1][2] 	= 0.0*deg;
    this->griffinCoords[2][2] 	= 0.0*deg;
    this->griffinCoords[3][2] 	= 0.0*deg;
    this->griffinCoords[4][2] 	= 0.0*deg;
    this->griffinCoords[5][2] 	= 0.0*deg;
    this->griffinCoords[6][2] 	= 0.0*deg;
    this->griffinCoords[7][2] 	= 0.0*deg;
    this->griffinCoords[8][2] 	= 0.0*deg;
    this->griffinCoords[9][2] 	= 0.0*deg;
    this->griffinCoords[10][2] 	= 0.0*deg;
    this->griffinCoords[11][2] 	= 0.0*deg;
    this->griffinCoords[12][2] 	= 0.0*deg;
    this->griffinCoords[13][2] 	= 0.0*deg;
    this->griffinCoords[14][2] 	= 0.0*deg;
    this->griffinCoords[15][2] 	= 0.0*deg;
    // pitch (beta)
    this->griffinCoords[0][3] 	= -45.0*deg;
    this->griffinCoords[1][3] 	= -45.0*deg;
    this->griffinCoords[2][3] 	= -45.0*deg;
    this->griffinCoords[3][3] 	= -45.0*deg;
    this->griffinCoords[4][3] 	= 0.0*deg;
    this->griffinCoords[5][3] 	= 0.0*deg;
    this->griffinCoords[6][3] 	= 0.0*deg;
    this->griffinCoords[7][3] 	= 0.0*deg;
    this->griffinCoords[8][3] 	= 0.0*deg;
    this->griffinCoords[9][3] 	= 0.0*deg;
    this->griffinCoords[10][3] 	= 0.0*deg;
    this->griffinCoords[11][3] 	= 0.0*deg;
    this->griffinCoords[12][3] 	= 45.0*deg;
    this->griffinCoords[13][3] 	= 45.0*deg;
    this->griffinCoords[14][3] 	= 45.0*deg;
    this->griffinCoords[15][3] 	= 45.0*deg;
    // roll (gamma)
    this->griffinCoords[0][4] 	= 67.5*deg;
    this->griffinCoords[1][4] 	= 157.5*deg;
    this->griffinCoords[2][4] 	= 247.5*deg;
    this->griffinCoords[3][4] 	= 337.5*deg;
    this->griffinCoords[4][4] 	= 22.5*deg;
    this->griffinCoords[5][4] 	= 67.5*deg;
    this->griffinCoords[6][4] 	= 112.5*deg;
    this->griffinCoords[7][4] 	= 157.5*deg;
    this->griffinCoords[8][4] 	= 202.5*deg;
    this->griffinCoords[9][4] 	= 247.5*deg;
    this->griffinCoords[10][4] 	= 292.5*deg;
    this->griffinCoords[11][4] 	= 337.5*deg;
    this->griffinCoords[12][4] 	= 67.5*deg;
    this->griffinCoords[13][4] 	= 157.5*deg;
    this->griffinCoords[14][4] 	= 247.5*deg;
    this->griffinCoords[15][4] 	= 337.5*deg;

    // Triangle sections for ancillary detectors coords
    // theta
    this->ancillaryCoords[0][0] 	= triangleThetaAngle;
    this->ancillaryCoords[1][0] 	= triangleThetaAngle;
    this->ancillaryCoords[2][0] 	= triangleThetaAngle;
    this->ancillaryCoords[3][0] 	= triangleThetaAngle;
    this->ancillaryCoords[4][0] 	= 180.0*deg - triangleThetaAngle;
    this->ancillaryCoords[5][0] 	= 180.0*deg - triangleThetaAngle;
    this->ancillaryCoords[6][0] 	= 180.0*deg - triangleThetaAngle;
    this->ancillaryCoords[7][0] 	= 180.0*deg - triangleThetaAngle;
    // phi
    this->ancillaryCoords[0][1] 	= 22.5*deg;
    this->ancillaryCoords[1][1] 	= 112.5*deg;
    this->ancillaryCoords[2][1] 	= 202.5*deg;
    this->ancillaryCoords[3][1] 	= 292.5*deg;
    this->ancillaryCoords[4][1] 	= 22.5*deg;
    this->ancillaryCoords[5][1] 	= 112.5*deg;
    this->ancillaryCoords[6][1] 	= 202.5*deg;
    this->ancillaryCoords[7][1] 	= 292.5*deg;
    // yaw (alpha)
    this->ancillaryCoords[0][2] 	= 0.0*deg;
    this->ancillaryCoords[1][2] 	= 0.0*deg;
    this->ancillaryCoords[2][2] 	= 0.0*deg;
    this->ancillaryCoords[3][2] 	= 0.0*deg;
    this->ancillaryCoords[4][2] 	= 0.0*deg;
    this->ancillaryCoords[5][2] 	= 0.0*deg;
    this->ancillaryCoords[6][2] 	= 0.0*deg;
    this->ancillaryCoords[7][2] 	= 0.0*deg;
    // pitch (beta)
    this->ancillaryCoords[0][3] 	= triangleThetaAngle;
    this->ancillaryCoords[1][3] 	= triangleThetaAngle;
    this->ancillaryCoords[2][3] 	= triangleThetaAngle;
    this->ancillaryCoords[3][3] 	= triangleThetaAngle;
    this->ancillaryCoords[4][3] 	= 180.0*deg - triangleThetaAngle;
    this->ancillaryCoords[5][3] 	= 180.0*deg - triangleThetaAngle;
    this->ancillaryCoords[6][3] 	= 180.0*deg - triangleThetaAngle;
    this->ancillaryCoords[7][3] 	= 180.0*deg - triangleThetaAngle;
    // roll (gamma)
    this->ancillaryCoords[0][4] 	= 22.5*deg;
    this->ancillaryCoords[1][4] 	= 112.5*deg;
    this->ancillaryCoords[2][4] 	= 202.5*deg;
    this->ancillaryCoords[3][4] 	= 292.5*deg;
    this->ancillaryCoords[4][4] 	= 22.5*deg;
    this->ancillaryCoords[5][4] 	= 112.5*deg;
    this->ancillaryCoords[6][4] 	= 202.5*deg;
    this->ancillaryCoords[7][4] 	= 292.5*deg;
}

ApparatusGriffinStructure::~ApparatusGriffinStructure()
{
    delete square_piece_log;
    delete triangle_piece_log;
    delete ring_piece_log;
    delete rod_piece_log;
}

G4int ApparatusGriffinStructure::Build()
{ 
    // Setup assembly volumes
    G4AssemblyVolume* myassemblySquare = new G4AssemblyVolume();
    this->assemblySquare = myassemblySquare;

    G4AssemblyVolume* myassemblyTriangle = new G4AssemblyVolume();
    this->assemblyTriangle = myassemblyTriangle;

    G4AssemblyVolume* myassemblyRing = new G4AssemblyVolume();
    this->assemblyRing = myassemblyRing;

    G4AssemblyVolume* myassemblyRod = new G4AssemblyVolume();
    this->assemblyRod = myassemblyRod;

    G4cout << "BuildSquarePiece" << G4endl;
    BuildSquarePiece();

    G4cout << "BuildTrianglePiece" << G4endl;
    BuildTrianglePiece();

    G4cout << "BuildLargeRingFrame" << G4endl;
    BuildLargeRingFrame();

    G4cout << "BuildRodPiece" << G4endl;
    BuildRodPiece();

    return 1;
}//end ::Build


G4int ApparatusGriffinStructure::Place(G4LogicalVolume* exp_hall_log, G4int selector)
{
    G4ThreeVector move = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
    G4RotationMatrix* rotate = new G4RotationMatrix;

    G4double theta;
    G4double phi;
    G4double alpha;
    G4double beta;
    G4double gamma;

    G4int iSquare1 = 0;
    G4int iSquare2 = 15;
    G4int iTriangle1 =0;
    G4int iTriangle2 =7;

    if(selector == 0) { // Full structure
        iSquare1 = 0;
        iSquare2 = 15;
        iTriangle1 =0;
        iTriangle2 =7;
    }
    else if(selector == 1) { // upstream-half of structure: corona plus upstream lampshade
        iSquare1 = 4;
        iSquare2 = 15;
        iTriangle1 =4;
        iTriangle2 =7;
    }
    else if(selector == 2) { // downstream-half of structure: corona plus downstream lampshade
        iSquare1 = 0;
        iSquare2 = 11;
        iTriangle1 =0;
        iTriangle2 =3;
    }
    else if(selector == 3) { // just the corona, no lampshades
        iSquare1 = 4;
        iSquare2 = 11;
        iTriangle1 =4;
        iTriangle2 =3;
    }

    for(G4int i=iSquare1; i<=iSquare2; i++){ // loop over all square pieces
        theta  = this->griffinCoords[i][0];
        phi    = this->griffinCoords[i][1];
        alpha  = this->griffinCoords[i][2]; // yaw
        beta   = this->griffinCoords[i][3]; // pitch
        gamma  = this->griffinCoords[i][4]; // roll

        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI/2.0);
        rotate->rotateX(M_PI/2.0);
        rotate->rotateX(alpha);
        rotate->rotateY(beta);
        rotate->rotateZ(gamma);

        this->assemblySquare->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);
    }

    for(G4int i=iSquare1; i<=iSquare2; i++){ // loop over all rods
        theta  = this->griffinCoords[i][0];
        phi    = this->griffinCoords[i][1];
        alpha  = this->griffinCoords[i][2]; // yaw
        beta   = this->griffinCoords[i][3]; // pitch
        gamma  = this->griffinCoords[i][4]; // roll

        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI/2.0);
        rotate->rotateX(M_PI/2.0);
        rotate->rotateX(alpha);
        rotate->rotateY(beta);
        rotate->rotateZ(gamma);

        this->assemblyRod->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);
    }

    for(G4int i=iTriangle1; i<=iTriangle2; i++){ // loop over all triangle pieces
        theta  = this->ancillaryCoords[i][0];
        phi    = this->ancillaryCoords[i][1];
        alpha  = this->ancillaryCoords[i][2]; // yaw
        beta   = this->ancillaryCoords[i][3]; // pitch
        gamma  = this->ancillaryCoords[i][4]; // roll

        rotate = new G4RotationMatrix;
        if(i>=4)
            rotate->rotateZ(60.0*deg);
        rotate->rotateY(M_PI);
        rotate->rotateY(M_PI+beta);
        rotate->rotateZ(gamma);

        this->assemblyTriangle->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);
    }

    rotate = new G4RotationMatrix;
    rotate->rotateY(0.0*deg);
    this->assemblyRing->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);

    rotate = new G4RotationMatrix;
    rotate->rotateY(180.0*deg);
    this->assemblyRing->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);

    return 1;
}

G4int ApparatusGriffinStructure::BuildSquarePiece()
{
    G4Material* material = G4Material::GetMaterial(this->structure_material);
    if( !material ) {
        G4cout << " ----> Material " << this->structure_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
    vis_att->SetVisibility(true);

    G4SubtractionSolid* square_piece = SquarePiece();

    // Define rotation and movement objects
    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4double z_position		= this->inner_distance_to_square_face + (this->square_and_triangle_face_thickness)/2.0;
    G4ThreeVector move 		= z_position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( square_piece_log == NULL )
    {
        square_piece_log = new G4LogicalVolume(square_piece, material, "square_piece_log", 0, 0, 0);
        square_piece_log->SetVisAttributes(vis_att);
    }

    this->assemblySquare->AddPlacedVolume(square_piece_log, move, rotate);

    return 1;
}

G4int ApparatusGriffinStructure::BuildTrianglePiece()
{
    G4Material* material = G4Material::GetMaterial(this->structure_material);
    if( !material ) {
        G4cout << " ----> Material " << this->structure_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
    vis_att->SetVisibility(true);

    G4SubtractionSolid* triangle_piece = TrianglePiece();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);

    G4double r          = this->inner_distance_to_square_face +  this->square_and_triangle_face_thickness/2.0; // to centre
    G4double x          = r*tan(22.5*deg); // half length of square side
    G4double h          = (x)*tan(30.0*deg); // height from bottom of the triangle to origin.
    G4double c          = x + h*sin(triangleThetaAngle);
    G4double z_position	= ((c)/(sin(90.0*deg - triangleThetaAngle)))  ;

    G4ThreeVector move 		= z_position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( triangle_piece_log == NULL )
    {
        triangle_piece_log = new G4LogicalVolume(triangle_piece, material, "triangle_piece_log", 0, 0, 0);
        triangle_piece_log->SetVisAttributes(vis_att);
    }

    this->assemblyTriangle->AddPlacedVolume(triangle_piece_log, move, rotate);

    return 1;
}

G4int ApparatusGriffinStructure::BuildLargeRingFrame()
{
    G4Material* material = G4Material::GetMaterial(this->structure_material);
    if( !material ) {
        G4cout << " ----> Material " << this->structure_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
    vis_att->SetVisibility(true);

    G4Tubs* ring_piece = LargeRingFrame();

    // Define rotation and movement objects
    G4ThreeVector direction 	= G4ThreeVector(0,0,1);
    G4double z_position         = this->large_ring_frame_z_position + this->large_ring_frame_thickness/2.0;
    G4ThreeVector move          = z_position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( ring_piece_log == NULL )
    {
        ring_piece_log = new G4LogicalVolume(ring_piece, material, "ring_piece_log", 0, 0, 0);
        ring_piece_log->SetVisAttributes(vis_att);
    }

    this->assemblyRing->AddPlacedVolume(ring_piece_log, move, rotate);

    return 1;
}

G4int ApparatusGriffinStructure::BuildRodPiece()
{
    G4Material* material = G4Material::GetMaterial(this->rod_material);
    if( !material ) {
        G4cout << " ----> Material " << this->rod_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.30,0.30,0.30));
    vis_att->SetVisibility(true);

    G4Tubs* rod_piece = RodPiece();

    // Define rotation and movement objects
    G4ThreeVector move 	= G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
    G4double z_position = sqrt( (this->square_rod_hole_inner_radius * this->square_rod_hole_inner_radius ) - (this->square_rod_position_from_center_of_square_face * this->square_rod_position_from_center_of_square_face) ) + this->square_rod_total_length/2.0 ;

    G4RotationMatrix* rotate = new G4RotationMatrix;

    //logical volume
    if( rod_piece_log == NULL )
    {
        rod_piece_log = new G4LogicalVolume(rod_piece, material, "rod_piece_log", 0, 0, 0);
        rod_piece_log->SetVisAttributes(vis_att);
    }

    move = G4ThreeVector(1.0*this->square_rod_position_from_center_of_square_face,1.0*this->square_rod_position_from_center_of_square_face,z_position);
    this->assemblyRod->AddPlacedVolume(rod_piece_log, move, rotate);

    move = G4ThreeVector(-1.0*this->square_rod_position_from_center_of_square_face,1.0*this->square_rod_position_from_center_of_square_face,z_position);
    this->assemblyRod->AddPlacedVolume(rod_piece_log, move, rotate);

    move = G4ThreeVector(1.0*this->square_rod_position_from_center_of_square_face,-1.0*this->square_rod_position_from_center_of_square_face,z_position);
    this->assemblyRod->AddPlacedVolume(rod_piece_log, move, rotate);

    move = G4ThreeVector(-1.0*this->square_rod_position_from_center_of_square_face,-1.0*this->square_rod_position_from_center_of_square_face,z_position);
    this->assemblyRod->AddPlacedVolume(rod_piece_log, move, rotate);

    return 1;
}

G4SubtractionSolid* ApparatusGriffinStructure::SquarePiece()
{
    G4double inner_radius;
    G4double outer_radius;
    G4double start_phi = 0.0*deg;
    G4double end_phi = 360.0*deg;

    G4double r              = this->inner_distance_to_square_face + this->square_and_triangle_face_thickness/2.0; // to centre of square solid
    G4double x              = r*tan(22.5*deg); // half length of square side
    G4double xt             = (this->square_and_triangle_face_thickness/2.0)*tan(22.5*deg); // add extra thickess for the edges
    G4double half_length_x  = x + xt;
    G4double half_length_y  = half_length_x;
    G4double half_length_z  = this->square_and_triangle_face_thickness/2.0;

    G4Box* griffin_structure_square_plate = new G4Box("griffin_structure_square_plate", half_length_x, half_length_y, half_length_z);

    G4double cut_half_length_x;
    G4double cut_half_length_y;
    G4double cut_half_length_z;
    G4double extra_cut;

    G4ThreeVector move_cut = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
    G4RotationMatrix* rotate_cut;
    rotate_cut = new G4RotationMatrix;

    //Now cut bgo slots
    extra_cut           = 50.0*mm;
    cut_half_length_x   = ((this->square_bgo_cut_width/2.0 + extra_cut)*cos(22.5*deg));
    cut_half_length_y   = this->square_bgo_cut_length/2.0;
    cut_half_length_z   = this->square_and_triangle_face_thickness;

    G4Box* cut_bgo = new G4Box("cut_bgo", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    // snip first bgo slot
    move_cut = G4ThreeVector(-1.0*(x-this->square_bgo_cut_width/2.0 + extra_cut),0.0*mm,0.0*mm); // half_length_x, the outer length of the original square
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateY(22.5*deg);

    G4SubtractionSolid* griffin_structure_square_plate_cut1 = new G4SubtractionSolid("griffin_structure_square_plate_cut1", griffin_structure_square_plate, cut_bgo, rotate_cut, move_cut);

    // snip first bgo slot
    move_cut = G4ThreeVector(1.0*(x-this->square_bgo_cut_width/2.0 + extra_cut),0.0*mm,0.0*mm); // half_length_x, the outer length of the original square
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(180.0*deg);
    rotate_cut->rotateY(22.5*deg);

    G4SubtractionSolid* griffin_structure_square_plate_cut2 = new G4SubtractionSolid("griffin_structure_square_plate_cut2", griffin_structure_square_plate_cut1, cut_bgo, rotate_cut, move_cut);

    // snip next bgo slot
    move_cut = G4ThreeVector(0.0*mm,1.0*(x-this->square_bgo_cut_width/2.0 + extra_cut),0.0*mm); // half_length_x, the outer length of the original square
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(90.0*deg);
    rotate_cut->rotateY(22.5*deg);

    G4SubtractionSolid* griffin_structure_square_plate_cut3 = new G4SubtractionSolid("griffin_structure_square_plate_cut3", griffin_structure_square_plate_cut2, cut_bgo, rotate_cut, move_cut);

    // snip last bgo slot
    move_cut = G4ThreeVector(0.0*mm,-1.0*(x-this->square_bgo_cut_width/2.0 + extra_cut),0.0*mm); // half_length_x, the outer length of the original square
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(270.0*deg);
    rotate_cut->rotateY(22.5*deg);

    G4SubtractionSolid* griffin_structure_square_plate_cut4 = new G4SubtractionSolid("griffin_structure_square_plate_cut4", griffin_structure_square_plate_cut3, cut_bgo, rotate_cut, move_cut);


    // Now cut the edges so that they fit nicely together
    cut_half_length_x = half_length_z;
    cut_half_length_y = 2.0*r;
    cut_half_length_z = 2.0*r;

    G4Box* cut_side = new G4Box("cut_side", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    move_cut = G4ThreeVector(-1.0*cut_half_length_x/cos(22.5*deg),0.0*mm,-1.0*r);
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateY(22.5*deg);

    // snip 1st edge
    G4SubtractionSolid* griffin_structure_square_plate_cut5 = new G4SubtractionSolid("griffin_structure_square_plate_cut5", griffin_structure_square_plate_cut4, cut_side, rotate_cut, move_cut);

    move_cut = G4ThreeVector(1.0*cut_half_length_x/cos(22.5*deg),0.0*mm,-1.0*r);
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(180.0*deg);
    rotate_cut->rotateY(22.5*deg);

    // snip next edge
    G4SubtractionSolid* griffin_structure_square_plate_cut6 = new G4SubtractionSolid("griffin_structure_square_plate_cut6", griffin_structure_square_plate_cut5, cut_side, rotate_cut, move_cut);

    move_cut = G4ThreeVector(0.0*mm,1.0*cut_half_length_x/cos(22.5*deg),-1.0*r);
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(90.0*deg);
    rotate_cut->rotateY(22.5*deg);

    // snip next edge
    G4SubtractionSolid* griffin_structure_square_plate_cut7 = new G4SubtractionSolid("griffin_structure_square_plate_cut7", griffin_structure_square_plate_cut6, cut_side, rotate_cut, move_cut);

    move_cut = G4ThreeVector(0.0*mm,-1.0*cut_half_length_x/cos(22.5*deg),-1.0*r);
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(270.0*deg);
    rotate_cut->rotateY(22.5*deg);

    // snip last edge
    G4SubtractionSolid* griffin_structure_square_plate_cut8 = new G4SubtractionSolid("griffin_structure_square_plate_cut8", griffin_structure_square_plate_cut7, cut_side, rotate_cut, move_cut);


    // Now cut the holes for the steel rods
    inner_radius = 0.0*mm;
    outer_radius = this->square_rod_hole_radius;
    half_length_z = this->square_rod_total_length/2.0;

    G4Tubs* rod_hole = new G4Tubs("rod_hole", inner_radius, outer_radius, half_length_z, start_phi, end_phi);


    G4double z_position_rod     = sqrt( (this->square_rod_hole_inner_radius * this->square_rod_hole_inner_radius ) - (this->square_rod_position_from_center_of_square_face * this->square_rod_position_from_center_of_square_face) ) + this->square_rod_total_length/2.0 ;
    G4double z_position_square  = this->inner_distance_to_square_face + (this->square_and_triangle_face_thickness)/2.0;

    // snip first hole
    move_cut = G4ThreeVector(1.0*this->square_rod_position_from_center_of_square_face,1.0*this->square_rod_position_from_center_of_square_face,z_position_rod-z_position_square);

    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* griffin_structure_square_plate_cut9 = new G4SubtractionSolid("griffin_structure_square_plate_cut9", griffin_structure_square_plate_cut8, rod_hole, rotate_cut, move_cut);

    // snip next hole
    move_cut = G4ThreeVector(-1.0*this->square_rod_position_from_center_of_square_face,1.0*this->square_rod_position_from_center_of_square_face,z_position_rod-z_position_square);
    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* griffin_structure_square_plate_cut10 = new G4SubtractionSolid("griffin_structure_square_plate_cut10", griffin_structure_square_plate_cut9, rod_hole, rotate_cut, move_cut);

    // snip next hole
    move_cut = G4ThreeVector(1.0*this->square_rod_position_from_center_of_square_face,-1.0*this->square_rod_position_from_center_of_square_face,z_position_rod-z_position_square);
    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* griffin_structure_square_plate_cut11 = new G4SubtractionSolid("griffin_structure_square_plate_cut11", griffin_structure_square_plate_cut10, rod_hole, rotate_cut, move_cut);

    // snip last hole
    move_cut = G4ThreeVector(-1.0*this->square_rod_position_from_center_of_square_face,-1.0*this->square_rod_position_from_center_of_square_face,z_position_rod-z_position_square);
    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* griffin_structure_square_plate_cut12 = new G4SubtractionSolid("griffin_structure_square_plate_cut12", griffin_structure_square_plate_cut11, rod_hole, rotate_cut, move_cut);


    // cut out the square hole with rounded corners
    move_cut = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* squareWithRoundedCorners = SquareWithRoundedCorners();

    G4SubtractionSolid* griffin_structure_square_plate_cut13 = new G4SubtractionSolid("griffin_structure_square_plate_cut13", griffin_structure_square_plate_cut12, squareWithRoundedCorners, rotate_cut, move_cut);

    return griffin_structure_square_plate_cut13;
}


G4SubtractionSolid* ApparatusGriffinStructure::TrianglePiece()
{
    // grab the scissors, time to cut.
    G4double start_phi = 0.0*deg;
    G4double end_phi = 360.0*deg;
    G4double phi = 0.0*deg;

    G4double inner_radius;
    G4double outer_radius;
    G4double cut_half_length_x;
    G4double cut_half_length_y;
    G4double cut_half_length_z;

    G4ThreeVector move, move_cut;
    G4RotationMatrix* rotate_cut;

    // to centre
    G4double r              = this->inner_distance_to_square_face + this->square_and_triangle_face_thickness/2.0;
    G4double x              = r*tan(22.5*deg);
    G4double xt             = (this->square_and_triangle_face_thickness/2.0)*tan(22.5*deg); // add thickness for edges
    G4double half_length_x  = x + xt;
    G4double half_length_z  = this->square_and_triangle_face_thickness/2.0;

    // equilateral triangle is inscribed in a circle of radius a.
    inner_radius = 0.0*mm;
    outer_radius = (2.0*half_length_x)/(sqrt(3));
    G4Tubs* griffin_structure_triangle_plate = new G4Tubs("griffin_structure_triangle_plate", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    G4SubtractionSolid* truncatedThreeSidedCylinder = TruncatedThreeSidedCylinder();

    move = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
    rotate_cut = new G4RotationMatrix;

    // cut out inside truncatedThreeSidedCylinder hole
    G4SubtractionSolid* griffin_structure_triangle_plate_cut1 = new G4SubtractionSolid("griffin_structure_triangle_plate_cut1", griffin_structure_triangle_plate, truncatedThreeSidedCylinder, rotate_cut, move);

    r = this->inner_distance_to_square_face + this->square_and_triangle_face_thickness/2.0;
    x = r*tan(22.5*deg);

    // Now cut the edges so that the pieces fit nicely together
    cut_half_length_x = 1.0*r;
    cut_half_length_y = 2.0*r;
    cut_half_length_z = 2.0*r;
    G4Box* cut_side = new G4Box("cut_side", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    G4double theta = ((90.0*deg - triangleThetaAngle) - 22.5*deg);

    move = G4ThreeVector(cut_half_length_x/cos(theta),0.0*mm,-1.0*(x*tan(30.0*deg))/(tan(theta))); // phi = 0;

    // 1st edge
    phi = 0.0*deg;

    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(-1.0*phi);
    rotate_cut->rotateY(-1.0*theta);

    move_cut 	= G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

    // snip snip
    G4SubtractionSolid* griffin_structure_triangle_plate_cut2 = new G4SubtractionSolid("griffin_structure_triangle_plate_cut2", griffin_structure_triangle_plate_cut1, cut_side, rotate_cut, move_cut);

    // 2nd edge
    phi = 120.0*deg;
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(-1.0*phi);
    rotate_cut->rotateY(-1.0*theta);

    move_cut 	= G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

    // snip snip
    G4SubtractionSolid* griffin_structure_triangle_plate_cut3 = new G4SubtractionSolid("griffin_structure_triangle_plate_cut3", griffin_structure_triangle_plate_cut2, cut_side, rotate_cut, move_cut);

    // last edge
    phi = 240.0*deg;

    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(-1.0*phi);
    rotate_cut->rotateY(-1.0*theta);

    move_cut 	= G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

    // snip snip
    G4SubtractionSolid* griffin_structure_triangle_plate_cut4 = new G4SubtractionSolid("griffin_structure_triangle_plate_cut4", griffin_structure_triangle_plate_cut3, cut_side, rotate_cut, move_cut);

    return griffin_structure_triangle_plate_cut4;
}

G4Tubs* ApparatusGriffinStructure::LargeRingFrame()
{
    G4double start_phi = 0.0*deg;
    G4double end_phi = 360.0*deg;

    G4double inner_radius = this->large_ring_frame_inner_radius;
    G4double outer_radius = this->large_ring_frame_outer_radius;
    G4double half_length_z = this->large_ring_frame_thickness/2.0;

    G4Tubs* ring_cylinder = new G4Tubs("ring_cylinder", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return ring_cylinder;
}

G4Tubs* ApparatusGriffinStructure::RodPiece()
{
    G4double start_phi = 0.0*deg;
    G4double end_phi = 360.0*deg;

    G4double inner_radius = 0.0*mm;
    G4double outer_radius = this->square_rod_hole_radius;
    G4double half_length_z = this->square_rod_total_length/2.0;

    G4Tubs* rod = new G4Tubs("rod", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    return rod;
}

G4SubtractionSolid* ApparatusGriffinStructure::SquareWithRoundedCorners()
{
    G4double half_length_x = this->square_hole_length/2.0;
    G4double half_length_y = this->square_hole_length/2.0;
    G4double half_length_z = this->square_and_triangle_face_thickness;

    G4ThreeVector move = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
    G4RotationMatrix* rotate;

    G4Box* griffin_structure_square_plate_1 = new G4Box("griffin_structure_square_plate_1", half_length_x, half_length_y, half_length_x);

    // cut rounded corners
    G4double start_phi = 0.0*deg;
    G4double end_phi = 90.0*deg;

    G4double inner_radius = this->square_hole_rounded_edge_radius;
    G4double outer_radius = this->square_hole_rounded_edge_radius + 2.0*half_length_x;

    G4Tubs* rounded_corner = new G4Tubs("rounded_corner", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // 1st corner
    move = G4ThreeVector(1.0*(half_length_x - this->square_hole_rounded_edge_radius),1.0*(half_length_x - this->square_hole_rounded_edge_radius),0.0*mm);
    rotate = new G4RotationMatrix;

    G4SubtractionSolid* griffin_structure_square_plate_2 = new G4SubtractionSolid("griffin_structure_square_plate_2", griffin_structure_square_plate_1, rounded_corner, rotate, move);

    // next corner
    move = G4ThreeVector(1.0*(half_length_x - this->square_hole_rounded_edge_radius),-1.0*(half_length_x - this->square_hole_rounded_edge_radius),0.0*mm);
    rotate = new G4RotationMatrix;
    rotate->rotateZ(90.0*deg);

    G4SubtractionSolid* griffin_structure_square_plate_3 = new G4SubtractionSolid("griffin_structure_square_plate_3", griffin_structure_square_plate_2, rounded_corner, rotate, move);

    // next corner
    move = G4ThreeVector(-1.0*(half_length_x - this->square_hole_rounded_edge_radius),-1.0*(half_length_x - this->square_hole_rounded_edge_radius),0.0*mm);
    rotate = new G4RotationMatrix;
    rotate->rotateZ(180.0*deg);

    G4SubtractionSolid* griffin_structure_square_plate_4 = new G4SubtractionSolid("griffin_structure_square_plate_4", griffin_structure_square_plate_3, rounded_corner, rotate, move);

    // last corner
    move = G4ThreeVector(-1.0*(half_length_x - this->square_hole_rounded_edge_radius),1.0*(half_length_x - this->square_hole_rounded_edge_radius),0.0*mm);
    rotate = new G4RotationMatrix;
    rotate->rotateZ(270.0*deg);

    G4SubtractionSolid* griffin_structure_square_plate_5 = new G4SubtractionSolid("griffin_structure_square_plate_5", griffin_structure_square_plate_4, rounded_corner, rotate, move);

    return griffin_structure_square_plate_5;
}

G4SubtractionSolid* ApparatusGriffinStructure::TruncatedThreeSidedCylinder()
{
    G4double start_phi = 0.0*deg;
    G4double end_phi = 360.0*deg;

    G4double inner_radius = 0.0*mm;
    G4double outer_radius = this->triangle_outer_hole_radius;
    G4double half_length_z = this->square_and_triangle_face_thickness;

    G4double phi;
    G4ThreeVector move, move_cut;
    G4RotationMatrix* rotate_cut;

    G4Tubs* hole_cylinder = new G4Tubs("hole_cylinder", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // cut flat sides
    G4double cut_half_length_x = 2.0*outer_radius;
    G4double cut_half_length_y = 2.0*outer_radius;
    G4double cut_half_length_z = 2.0*half_length_z;

    G4Box* cut_edge = new G4Box("cut_edge", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    move = G4ThreeVector(this->triangle_inner_hole_radius + cut_half_length_x, 0.0*mm,0.0*mm);

    // 1st cut
    phi = 0.0*deg;
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(-1.0*phi);

    move_cut = G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

    G4SubtractionSolid* hole_cylinder_1 = new G4SubtractionSolid("hole_cylinder_1", hole_cylinder, cut_edge, rotate_cut, move_cut);

    // 2nd cut
    phi = 120.0*deg;
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(-1.0*phi);

    move_cut = G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

    G4SubtractionSolid* hole_cylinder_2 = new G4SubtractionSolid("hole_cylinder_2", hole_cylinder_1, cut_edge, rotate_cut, move_cut);

    // 3rd cut
    phi = 240.0*deg;
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(-1.0*phi);

    move_cut = G4ThreeVector(move.x()*cos(phi),move.x()*sin(phi),move.z());

    G4SubtractionSolid* hole_cylinder_3 = new G4SubtractionSolid("hole_cylinder_3", hole_cylinder_2, cut_edge, rotate_cut, move_cut);

    return hole_cylinder_3;
}

// General Spherical Translations
G4double ApparatusGriffinStructure::transX(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*cos(phi) );
}

G4double ApparatusGriffinStructure::transY(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*sin(phi) );
}

G4double ApparatusGriffinStructure::transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*cos(theta) );
}
