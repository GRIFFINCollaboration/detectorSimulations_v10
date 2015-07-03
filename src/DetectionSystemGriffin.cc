#include <math.h>

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemGriffin.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


///////////////////////////////////////////////////////////////////////
// The ::DetectionSystemGriffin constructor instatiates all the 
// Logical and Physical Volumes used in the detector geometery (found
// in our local files as it contains NDA-protected info), and the
// ::~DetectionSystemGriffin destructor deletes them from the stack 
// when they go out of scope
///////////////////////////////////////////////////////////////////////

DetectionSystemGriffin::~DetectionSystemGriffin()
{

    // LogicalVolumes in ConstructNewHeavyMet
    delete hevimet_log;

    // LogicalVolumes in ConstructNewSuppressorCasing
    delete left_suppressor_extension_log;
    delete right_suppressor_extension_log;
    delete left_suppressor_log;
    delete right_suppressor_log;
    delete back_quarter_suppressor_log;

    delete back_quarter_suppressor_shell_log;
    delete right_suppressor_shell_log;
    delete left_suppressor_shell_log;
    delete right_suppressor_shell_extension_log;
    delete left_suppressor_shell_extension_log;
    // LogicalVolumes in ConstructBGOCasing

    delete BGO_casing_log;
    delete back_BGO_log;

    // LogicalVolumes in ConstructColdFinger
    delete cooling_side_block_log;
    delete structureMat_cold_finger_log;
    delete cooling_bar_log;
    delete fet_air_hole_log;
    delete triangle_post_log;
    delete extra_cold_block_log;
    delete finger_log;
    delete end_plate_log;

    // LogicalVolumes in ConstructComplexDetectorBlock
    delete germanium_hole_log;
    delete germanium_block1_log;
    delete inner_dead_layer_log;
    delete inter_crystal_electrodeMat_back_log;
    delete inter_crystal_electrodeMat_front_log;

    // LogicalVolumes in ConstructBasicDetectorBlock
    delete germanium_block_log;

    // LogicalVolumes in ConstructDetector
    delete tank_log;
    delete tank_lid1_log;
    delete tank_lid2_log;
    delete finger_shell_log;
    delete rear_plate_log;
    delete bottom_side_panel_log;
    delete top_side_panel_log;
    delete left_side_panel_log;
    delete right_side_panel_log;
    delete lower_left_tube_log;
    delete upper_left_tube_log;
    delete lower_right_tube_log;
    delete upper_right_tube_log;
    delete lower_left_cone_log;
    delete upper_left_cone_log;
    delete lower_right_cone_log;
    delete upper_right_cone_log;
    delete bottom_wedge_log;
    delete top_wedge_log;
    delete left_wedge_log;
    delete right_wedge_log;
    delete bottom_bent_piece_log;
    delete top_bent_piece_log;
    delete left_bent_piece_log;
    delete right_bent_piece_log;
    delete front_face_log;
}// end ::~DetectionSystemGriffin

///////////////////////////////////////////////////////////////////////
// ConstructDetectionSystemGriffin builds the DetectionSystemGriffin 
// at the origin
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::Build()
{ 

    surfCheck = true;

    BuildOneDetector();

}//end ::Build

G4double DetectionSystemGriffin::transX(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return (x*cos(theta)+z*sin(theta))*cos(phi)-y*sin(phi);
}

G4double DetectionSystemGriffin::transY(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return (x*cos(theta)+z*sin(theta))*sin(phi)+y*cos(phi);
}

G4double DetectionSystemGriffin::transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return -x*sin(theta)+z*cos(theta);
}

// PlaceEverythingButCrystals and PlaceDeadLayerSpecificCrystal repalces the need for this method
G4int DetectionSystemGriffin::PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector moveBAH, G4RotationMatrix* rotateBAH, G4int detector_number)
{
//    G4double theta  = this->coords[detector_number][0]*deg;
//    G4double phi    = this->coords[detector_number][1]*deg;
//    G4double alpha  = this->coords[detector_number][2]*deg; // yaw
//    G4double beta   = this->coords[detector_number][3]*deg; // pitch
//    G4double gamma  = this->coords[detector_number][4]*deg; // roll

//    G4double x;
//    G4double y;
//    G4double z;

//    G4double x0, y0, z0;

//    G4int i;

//    G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
//    rotate->rotateX(M_PI/2.0);
//    rotate->rotateX(alpha);
//    rotate->rotateY(beta);
//    rotate->rotateZ(gamma);

//    // positioning
//    G4double dist_from_origin = this->air_box_back_length/2.0 +this->air_box_front_length + this->new_rhombi_radius ;
//    G4double dist_from_origin_det = this->air_box_back_length_det/2.0 +this->air_box_front_length_det + this->new_rhombi_radius_det ;

//    x = 0;
//    y = 0;
//    z = dist_from_origin_det;

//    G4ThreeVector move(DetectionSystemGriffin::transX(x,y,z,theta,phi), DetectionSystemGriffin::transY(x,y,z,theta,phi), DetectionSystemGriffin::transZ(x,y,z,theta,phi));

//    this->assembly->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);
//    this->backAndSideSuppressorShellAssembly->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);

//    z = dist_from_origin ; // This isolates the motion of the extension suppressors from the rest of the suppressors.

//    move = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi), DetectionSystemGriffin::transY(x,y,z,theta,phi), DetectionSystemGriffin::transZ(x,y,z,theta,phi));

//    this->extensionSuppressorShellAssembly->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);


//    x0 = (this->germanium_width + this->germanium_separation)/2.0;
//    y0 = (this->germanium_width + this->germanium_separation)/2.0;
//    z0 = this->germanium_length / 2.0 + this->can_face_thickness / 2.0 + this->germanium_dist_from_can_face + this->shift + this->applied_back_shift + dist_from_origin_det;

//    /////////////////////////////////////////////////////////////////////
//    // now we place all 4 of the 1/4 detectors using the LogicalVolume
//    // of the 1st, ensuring that they too will have the hole
//    /////////////////////////////////////////////////////////////////////
//    copy_number = germanium_copy_number + detector_number*4;

//    G4RotationMatrix* rotate_germanium[4];
//    G4ThreeVector move_germanium[4];

//    for(i=0; i<4; i++){
//        rotate_germanium[i] = new G4RotationMatrix;

//        rotate_germanium[i]->rotateY( -M_PI/2.0 ) ;
//        rotate_germanium[i]->rotateX( -M_PI/2.0 * i ) ;
//        rotate_germanium[i]->rotateX( alpha ) ;
//        rotate_germanium[i]->rotateY( beta ) ;
//        rotate_germanium[i]->rotateZ( gamma ) ;

//        x = x0*pow(-1, floor(i/2) );
//        y = y0*pow(-1, floor((i+1)/2) );

//        z = z0;

//        move_germanium[i] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi), DetectionSystemGriffin::transY(x,y,z,theta,phi), DetectionSystemGriffin::transZ(x,y,z,theta,phi));
//        this->germaniumAssembly->MakeImprint(exp_hall_log, move_germanium[i], rotate_germanium[i], copy_number++, surfCheck);

//    }
    
//    /////////////////////////////////////////////////////////////////////
//    // end germanium_block1_log
//    /////////////////////////////////////////////////////////////////////
//    copy_number = back_suppressor_copy_number + detector_number*4;

//    // Back Suppressors
//    if(include_back_suppressors) {
//        x0 = (this->detector_total_width/4.0);
//        y0 = (this->detector_total_width/4.0);
//        z0 = (this->back_BGO_thickness - this->can_face_thickness)/2.0 + this->suppressor_shell_thickness
//                + this->detector_total_length +this->BGO_can_seperation + this->shift + this->applied_back_shift + dist_from_origin_det;

//        G4RotationMatrix* rotate_back_quarter_suppressor[4];
//        G4ThreeVector move_back_quarter_suppressor[4];

//        for(i=0; i<4; i++){
//            rotate_back_quarter_suppressor[i] = new G4RotationMatrix;
//            rotate_back_quarter_suppressor[i]->rotateX(-M_PI/2.0*i);
//            rotate_back_quarter_suppressor[i]->rotateX(alpha);
//            rotate_back_quarter_suppressor[i]->rotateY(beta);
//            rotate_back_quarter_suppressor[i]->rotateZ(gamma);

//            x = x0*pow(-1, floor(i/2) );
//            y = y0*pow(-1, floor((i+1)/2) );
//            z = z0;

//            move_back_quarter_suppressor[i] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));

//            this->suppressorBackAssembly->MakeImprint(exp_hall_log, move_back_quarter_suppressor[i], rotate_back_quarter_suppressor[i], copy_number++, surfCheck);
//        }
//    }

//    // Now Side Suppressors
//    /////////////////////////////////////////////////////////////////////
//    // Note : Left and Right are read from the BACK of the detector
//    // Suppressors 1 and 2 cover germanium 1
//    // Suppressors 3 and 4 cover germanium 2
//    // Suppressors 5 and 6 cover germanium 3
//    // Suppressors 7 and 8 cover germanium 4
//    /////////////////////////////////////////////////////////////////////
//    copy_number = right_suppressor_side_copy_number + detector_number*4;
//    copy_number_two = left_suppressor_side_copy_number + detector_number*4;

//    // Replacement for this->side_suppressor_length
//    G4double shell_side_suppressor_length = this->side_suppressor_length
//            + (this->suppressor_shell_thickness*2.0);

//    // Replacement for this->suppressor_extension_length
//    G4double shell_suppressor_extension_length = this->suppressor_extension_length
//            + (this->suppressor_shell_thickness*2.0)*(1.0/tan(this->bent_end_angle)
//                                                      - tan(this->bent_end_angle));

//    // Replacement for this->suppressor_extension_angle: must totally recalculate
//    G4double shell_suppressor_extension_angle = atan((( this->suppressor_back_radius
//                                                        + this->bent_end_length +(this->BGO_can_seperation
//                                                                                  + this->side_BGO_thickness + this->suppressor_shell_thickness*2.0)
//                                                        / tan(this->bent_end_angle)
//                                                        - ( this->suppressor_extension_thickness + this->suppressor_shell_thickness * 2.0 )
//                                                        * sin(this->bent_end_angle))
//                                                      * tan(this->bent_end_angle) -( this->suppressor_forward_radius + this->hevimet_tip_thickness )
//                                                      * sin(this->bent_end_angle))/(shell_suppressor_extension_length));

//    // these two parameters are for shifting the extensions back and out when in their BACK position
//    G4double extension_back_shift = this->air_box_front_length
//            - (this->hevimet_tip_thickness + shell_suppressor_extension_length
//               + (this->suppressor_extension_thickness + this->suppressor_shell_thickness*2.0)
//               * tan(this->bent_end_angle))
//            * cos(this->bent_end_angle);

//    G4double extension_radial_shift = extension_back_shift*tan(this->bent_end_angle);

//    if(include_side_suppressors) {

//        G4RotationMatrix* rotateSideSuppressor[8];
//        G4ThreeVector moveInnerSuppressor[8];

//        x0 = this->side_BGO_thickness/2.0 + this->suppressor_shell_thickness + this->detector_total_width/2.0 +this->BGO_can_seperation;
//        y0 = (this->detector_total_width/2.0 +this->BGO_can_seperation + this->side_BGO_thickness/2.0)/2.0;
//        z0 = (shell_side_suppressor_length/2.0 -this->can_face_thickness/2.0 + this->bent_end_length
//              +(this->BGO_can_seperation + this->BGO_chopped_tip)/tan(this->bent_end_angle) +this->shift
//              + this->applied_back_shift - this->suppressor_shell_thickness/2.0 + dist_from_origin_det);

//        for(i=0; i<4; i++){

//            rotateSideSuppressor[2*i] = new G4RotationMatrix;
//            rotateSideSuppressor[2*i]->rotateZ(M_PI/2.0);
//            rotateSideSuppressor[2*i]->rotateY(M_PI/2.0);
//            rotateSideSuppressor[2*i]->rotateX(-M_PI/2.0*i);
//            rotateSideSuppressor[2*i]->rotateX(alpha);
//            rotateSideSuppressor[2*i]->rotateY(beta);
//            rotateSideSuppressor[2*i]->rotateZ(gamma);

//            x = x0*cos(i*M_PI/2.0) + y0*sin(i*M_PI/2);
//            y = -x0*sin(i*M_PI/2.0) + y0*cos(i*M_PI/2);
//            z = z0;

//            moveInnerSuppressor[i*2] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));

//            this->rightSuppressorCasingAssembly->MakeImprint(exp_hall_log, moveInnerSuppressor[2*i], rotateSideSuppressor[2*i], copy_number++, surfCheck);


//            rotateSideSuppressor[2*i+1] = new G4RotationMatrix;
//            rotateSideSuppressor[2*i+1]->rotateY(-M_PI/2.0);
//            rotateSideSuppressor[2*i+1]->rotateX(-M_PI/2.0*(i+2));
//            rotateSideSuppressor[2*i+1]->rotateX(alpha);
//            rotateSideSuppressor[2*i+1]->rotateY(beta);
//            rotateSideSuppressor[2*i+1]->rotateZ(gamma);

//            x = x0*sin(i*M_PI/2.0) + y0*cos(i*M_PI/2);
//            y = x0*cos(i*M_PI/2.0) - y0*sin(i*M_PI/2);
//            z = z0;

//            moveInnerSuppressor[i*2+1] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));

//            this->leftSuppressorCasingAssembly->MakeImprint(exp_hall_log, moveInnerSuppressor[2*i+1], rotateSideSuppressor[2*i+1], copy_number_two++, surfCheck);
//        }

//    }

//    // now we add the side pieces of suppressor that extend out in front of the can when it's in the back position
//    copy_number = right_suppressor_extension_copy_number + detector_number*4;
//    copy_number_two = left_suppressor_extension_copy_number + detector_number*4;

//    G4RotationMatrix* rotateExtension[8];
//    G4ThreeVector moveInnerExtension[8];

//    // this->back_inner_radius changed to this->suppressor_back_radius
//    // this->forward_inner_radius changed to this->suppressor_forward_radius
//    x0 =  - ((this->suppressor_extension_thickness/2.0 + this->suppressor_shell_thickness)
//             / cos(this->bent_end_angle)
//             + (shell_suppressor_extension_length/2.0 - (this->suppressor_extension_thickness
//                                                         + this->suppressor_shell_thickness*2.0)
//                * tan(this->bent_end_angle)/2.0) * sin(this->bent_end_angle) - (this->suppressor_back_radius
//                                                                                + this->bent_end_length + (this->BGO_can_seperation + this->side_BGO_thickness
//                                                                                                           + this->suppressor_shell_thickness*2.0)
//                                                                                / tan(this->bent_end_angle) - (this->suppressor_extension_thickness
//                                                                                                               + this->suppressor_shell_thickness*2.0)
//                                                                                * sin(this->bent_end_angle))
//             * tan(this->bent_end_angle));


//    // This is the exact same as the y0 when the suppressors are forward.
//    y0 =  (shell_suppressor_extension_length*tan(shell_suppressor_extension_angle)/2.0
//           + (this->suppressor_forward_radius
//              + this->hevimet_tip_thickness)*sin(this->bent_end_angle) - this->BGO_can_seperation*2.0 - this->detectorPlacementCxn)/2.0;

//    // only difference between this z0 and the detectors forward z0 is that this one does
//    // not have extension_back_shift added.
//    z0 =  - this->can_face_thickness/2.0 -(shell_suppressor_extension_length/2.0
//                                           - (this->suppressor_extension_thickness
//                                              + this->suppressor_shell_thickness*2.0)
//                                           * tan(this->bent_end_angle)/2.0)
//            * cos(this->bent_end_angle) +this->bent_end_length +(this->BGO_can_seperation
//                                                                 + this->side_BGO_thickness
//                                                                 + this->suppressor_shell_thickness*2.0)/tan(this->bent_end_angle)
//            - (this->suppressor_extension_thickness + this->suppressor_shell_thickness*2.0)
//            * sin(this->bent_end_angle) +this->suppShift + this->suppressor_back_radius
//            - this->suppressor_forward_radius
//            - this->suppressor_shell_thickness/2.0 + dist_from_origin;



//    if( !(this->suppressor_position_selector) && include_extension_suppressors )
//    {

//        x0 += extension_radial_shift ;

//        z0 += extension_back_shift ;

//        for( i = 0 ; i < 4 ; i++ )
//        {
//            rotateExtension[2*i] = new G4RotationMatrix;
//            rotateExtension[2*i]->rotateZ(M_PI/2.0);
//            rotateExtension[2*i]->rotateY(this->bent_end_angle);
//            rotateExtension[2*i]->rotateX(-M_PI/2.0*i);
//            rotateExtension[2*i]->rotateX(alpha);
//            rotateExtension[2*i]->rotateY(beta);
//            rotateExtension[2*i]->rotateZ(gamma);

//            x = x0*cos(i*M_PI/2.0) + y0*sin(i*M_PI/2);
//            y = -x0*sin(i*M_PI/2.0) + y0*cos(i*M_PI/2);
//            z = z0;

//            moveInnerExtension[2*i] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));
//            this->rightSuppressorExtensionAssembly->MakeImprint(exp_hall_log, moveInnerExtension[2*i], rotateExtension[2*i], copy_number++, surfCheck);

//            rotateExtension[2*i+1] = new G4RotationMatrix;
//            rotateExtension[2*i+1]->rotateY(M_PI/2.0);
//            rotateExtension[2*i+1]->rotateZ(M_PI/2.0 + this->bent_end_angle);
//            rotateExtension[2*i+1]->rotateX(-M_PI/2.0*i);
//            rotateExtension[2*i+1]->rotateX(alpha);
//            rotateExtension[2*i+1]->rotateY(beta);
//            rotateExtension[2*i+1]->rotateZ(gamma);

//            x = x0*sin(i*M_PI/2.0) + y0*cos(i*M_PI/2);
//            y = x0*cos(i*M_PI/2.0) - y0*sin(i*M_PI/2);
//            z = z0;

//            moveInnerExtension[2*i+1] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));
//            this->leftSuppressorExtensionAssembly->MakeImprint(exp_hall_log, moveInnerExtension[2*i+1], rotateExtension[2*i+1], copy_number_two++, surfCheck);
//        }

//    }//end if(detectors forward) statement

//    // Otherwise, put them forward

//    else if( this->suppressor_position_selector && include_extension_suppressors)
//    {

//        for( i = 0 ; i < 4 ; i++ )
//        {
//            rotateExtension[2*i] = new G4RotationMatrix;
//            rotateExtension[2*i]->rotateZ(M_PI/2.0);
//            rotateExtension[2*i]->rotateY(this->bent_end_angle);
//            rotateExtension[2*i]->rotateX(-M_PI/2.0*i);
//            rotateExtension[2*i]->rotateX(alpha);
//            rotateExtension[2*i]->rotateY(beta);
//            rotateExtension[2*i]->rotateZ(gamma);

//            x = x0*cos(i*M_PI/2.0) + y0*sin(i*M_PI/2);
//            y = -x0*sin(i*M_PI/2.0) + y0*cos(i*M_PI/2);
//            z = z0;

//            moveInnerExtension[2*i] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));
//            this->rightSuppressorExtensionAssembly->MakeImprint(exp_hall_log, moveInnerExtension[2*i], rotateExtension[2*i], copy_number++, surfCheck);

//            rotateExtension[2*i+1] = new G4RotationMatrix;
//            rotateExtension[2*i+1]->rotateY(M_PI/2.0);
//            rotateExtension[2*i+1]->rotateZ(M_PI/2.0 + this->bent_end_angle);
//            rotateExtension[2*i+1]->rotateX(-M_PI/2.0*i);
//            rotateExtension[2*i+1]->rotateX(alpha);
//            rotateExtension[2*i+1]->rotateY(beta);
//            rotateExtension[2*i+1]->rotateZ(gamma);

//            x = x0*sin(i*M_PI/2.0) + y0*cos(i*M_PI/2);
//            y = x0*cos(i*M_PI/2.0) - y0*sin(i*M_PI/2);
//            z = z0;


//            moveInnerExtension[2*i+1] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));
//            this->leftSuppressorExtensionAssembly->MakeImprint(exp_hall_log, moveInnerExtension[2*i+1], rotateExtension[2*i+1], copy_number_two++, surfCheck);
//        }

//    }//end if(detectors back) statement
//    return 1;
}
// end PlaceDetector    

G4int DetectionSystemGriffin::PlaceEverythingButCrystals(G4LogicalVolume* exp_hall_log, G4int detector_number, G4int position_number, G4bool posTigress)
{
    G4double theta  = this->coords[position_number][0]*deg;
    G4double phi    = this->coords[position_number][1]*deg;
    G4double alpha  = this->coords[position_number][2]*deg; // yaw
    G4double beta   = this->coords[position_number][3]*deg; // pitch
    G4double gamma  = this->coords[position_number][4]*deg; // roll

    //In GRIFFIN upstream is now downstream relative to TIGRESS. However, we want to maintain the same numbering scheme as TIGRESS. Net result is that the lampshade angles change by 45 degrees.
    if(posTigress && (position_number < 4 || position_number > 11) ) {
        phi       = this->coords[position_number][1]*deg - 45.0*deg;
        gamma     = this->coords[position_number][4]*deg - 45.0*deg;
    }

    G4double x;
    G4double y;
    G4double z;

    G4double x0, y0, z0;

    G4int i;

    G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
    rotate->rotateX(M_PI/2.0);
    rotate->rotateX(alpha);
    rotate->rotateY(beta);
    rotate->rotateZ(gamma);

    // positioning
    G4double dist_from_origin = this->air_box_back_length/2.0 +this->air_box_front_length + this->new_rhombi_radius ;
    G4double dist_from_origin_det = this->air_box_back_length_det/2.0 +this->air_box_front_length_det + this->new_rhombi_radius_det ;

    x = 0;
    y = 0;
    z = dist_from_origin_det;

    G4ThreeVector move(DetectionSystemGriffin::transX(x,y,z,theta,phi), DetectionSystemGriffin::transY(x,y,z,theta,phi), DetectionSystemGriffin::transZ(x,y,z,theta,phi));

    this->assembly->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);
    this->backAndSideSuppressorShellAssembly->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);

    z = dist_from_origin ; // This isolates the motion of the extension suppressors from the rest of the suppressors.

    move = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi), DetectionSystemGriffin::transY(x,y,z,theta,phi), DetectionSystemGriffin::transZ(x,y,z,theta,phi));

    this->extensionSuppressorShellAssembly->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck);
    this->hevimetAssembly->MakeImprint(exp_hall_log, move, rotate, 0, surfCheck) ;

    copy_number = back_suppressor_copy_number + detector_number*4;

    // Back Suppressors
    if(include_back_suppressors)
    {
        x0 = (this->detector_total_width/4.0);
        y0 = (this->detector_total_width/4.0);
        z0 = (this->back_BGO_thickness - this->can_face_thickness)/2.0 + this->suppressor_shell_thickness
                + this->detector_total_length +this->BGO_can_seperation + this->shift + this->applied_back_shift + dist_from_origin_det;

        G4RotationMatrix* rotate_back_quarter_suppressor[4];
        G4ThreeVector move_back_quarter_suppressor[4];

        for(i=0; i<4; i++){

            rotate_back_quarter_suppressor[i] = new G4RotationMatrix;
            rotate_back_quarter_suppressor[i]->rotateX(M_PI/2.0-M_PI/2.0*i);
            rotate_back_quarter_suppressor[i]->rotateX(alpha);
            rotate_back_quarter_suppressor[i]->rotateY(beta);
            rotate_back_quarter_suppressor[i]->rotateZ(gamma);

            x = -x0*pow(-1, floor((i+1)/2) );
            y = -y0*pow(-1, floor((i+2)/2) );
            z = z0;

            move_back_quarter_suppressor[i] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));

            this->suppressorBackAssembly->MakeImprint(exp_hall_log, move_back_quarter_suppressor[i], rotate_back_quarter_suppressor[i], copy_number++, surfCheck);
        }
    }

    // Now Side Suppressors
    /////////////////////////////////////////////////////////////////////
    // Note : Left and Right are read from the BACK of the detector
    // Suppressors 1 and 2 cover germanium 1
    // Suppressors 3 and 4 cover germanium 2
    // Suppressors 5 and 6 cover germanium 3
    // Suppressors 7 and 8 cover germanium 4
    /////////////////////////////////////////////////////////////////////
    copy_number = right_suppressor_side_copy_number + detector_number*4;
    copy_number_two = left_suppressor_side_copy_number + detector_number*4;

    // Replacement for this->side_suppressor_length
    G4double shell_side_suppressor_length = this->side_suppressor_length
            + (this->suppressor_shell_thickness*2.0);

    // Replacement for this->suppressor_extension_length
    G4double shell_suppressor_extension_length = this->suppressor_extension_length
            + (this->suppressor_shell_thickness*2.0)*(1.0/tan(this->bent_end_angle)
                                                      - tan(this->bent_end_angle));

    // Replacement for this->suppressor_extension_angle: must totally recalculate
    G4double shell_suppressor_extension_angle = atan((( this->suppressor_back_radius
                                                        + this->bent_end_length +(this->BGO_can_seperation
                                                                                  + this->side_BGO_thickness + this->suppressor_shell_thickness*2.0)
                                                        / tan(this->bent_end_angle)
                                                        - ( this->suppressor_extension_thickness + this->suppressor_shell_thickness * 2.0 )
                                                        * sin(this->bent_end_angle))
                                                      * tan(this->bent_end_angle) -( this->suppressor_forward_radius + this->hevimet_tip_thickness )
                                                      * sin(this->bent_end_angle))/(shell_suppressor_extension_length));

    // these two parameters are for shifting the extensions back and out when in their BACK position
    G4double extension_back_shift = this->air_box_front_length
            - (this->hevimet_tip_thickness + shell_suppressor_extension_length
               + (this->suppressor_extension_thickness + this->suppressor_shell_thickness*2.0)
               * tan(this->bent_end_angle))
            * cos(this->bent_end_angle);

    G4double extension_radial_shift = extension_back_shift*tan(this->bent_end_angle);

    if(include_side_suppressors)
    {

        G4RotationMatrix* rotateSideSuppressor[8];
        G4ThreeVector moveInnerSuppressor[8];

        x0 = this->side_BGO_thickness/2.0 + this->suppressor_shell_thickness + this->detector_total_width/2.0 +this->BGO_can_seperation;
        y0 = (this->detector_total_width/2.0 +this->BGO_can_seperation + this->side_BGO_thickness/2.0)/2.0;
        z0 = (shell_side_suppressor_length/2.0 -this->can_face_thickness/2.0 + this->bent_end_length
              +(this->BGO_can_seperation + this->BGO_chopped_tip)/tan(this->bent_end_angle) +this->shift
              + this->applied_back_shift - this->suppressor_shell_thickness/2.0 + dist_from_origin_det);

        for(i=0; i<4; i++)
        {

            rotateSideSuppressor[2*i] = new G4RotationMatrix;
            rotateSideSuppressor[2*i]->rotateZ(M_PI/2.0);
            rotateSideSuppressor[2*i]->rotateY(M_PI/2.0);
            rotateSideSuppressor[2*i]->rotateX(M_PI/2.0-M_PI/2.0*i);
            rotateSideSuppressor[2*i]->rotateX(alpha);
            rotateSideSuppressor[2*i]->rotateY(beta);
            rotateSideSuppressor[2*i]->rotateZ(gamma);

            x = x0*cos((i-1)*M_PI/2.0) + y0*sin((i-1)*M_PI/2);
            y = -x0*sin((i-1)*M_PI/2.0) + y0*cos((i-1)*M_PI/2);
            z = z0;

            moveInnerSuppressor[i*2] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));

            this->rightSuppressorCasingAssembly->MakeImprint(exp_hall_log, moveInnerSuppressor[2*i], rotateSideSuppressor[2*i], copy_number++, surfCheck);

            rotateSideSuppressor[2*i+1] = new G4RotationMatrix;
            rotateSideSuppressor[2*i+1]->rotateY(-M_PI/2.0);
            rotateSideSuppressor[2*i+1]->rotateX(-M_PI/2.0*(i+2));
            rotateSideSuppressor[2*i+1]->rotateX(alpha);
            rotateSideSuppressor[2*i+1]->rotateY(beta);
            rotateSideSuppressor[2*i+1]->rotateZ(gamma);

            x = x0*sin(i*M_PI/2.0) + y0*cos(i*M_PI/2);
            y = x0*cos(i*M_PI/2.0) - y0*sin(i*M_PI/2);
            z = z0;

            moveInnerSuppressor[i*2+1] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));

            this->leftSuppressorCasingAssembly->MakeImprint(exp_hall_log, moveInnerSuppressor[2*i+1], rotateSideSuppressor[2*i+1], copy_number_two++, surfCheck);
        }

    }

    // now we add the side pieces of suppressor that extend out in front of the can when it's in the back position
    copy_number = right_suppressor_extension_copy_number + detector_number*4;
    copy_number_two = left_suppressor_extension_copy_number + detector_number*4;

    G4RotationMatrix* rotateExtension[8];
    G4ThreeVector moveInnerExtension[8];

    x0 =  - ((this->suppressor_extension_thickness/2.0 + this->suppressor_shell_thickness)
             / cos(this->bent_end_angle)
             + (shell_suppressor_extension_length/2.0 - (this->suppressor_extension_thickness
                                                         + this->suppressor_shell_thickness*2.0)
                * tan(this->bent_end_angle)/2.0) * sin(this->bent_end_angle) - (this->suppressor_back_radius
                                                                                + this->bent_end_length + (this->BGO_can_seperation + this->side_BGO_thickness
                                                                                                           + this->suppressor_shell_thickness*2.0)
                                                                                / tan(this->bent_end_angle) - (this->suppressor_extension_thickness
                                                                                                               + this->suppressor_shell_thickness*2.0)
                                                                                * sin(this->bent_end_angle))
             * tan(this->bent_end_angle));

    y0 =  (shell_suppressor_extension_length*tan(shell_suppressor_extension_angle)/2.0
           + (this->suppressor_forward_radius
              + this->hevimet_tip_thickness)*sin(this->bent_end_angle) - this->BGO_can_seperation*2.0 - this->detectorPlacementCxn)/2.0;

    z0 =  - this->can_face_thickness/2.0 -(shell_suppressor_extension_length/2.0
                                           - (this->suppressor_extension_thickness
                                              + this->suppressor_shell_thickness*2.0)
                                           * tan(this->bent_end_angle)/2.0)
            * cos(this->bent_end_angle) +this->bent_end_length +(this->BGO_can_seperation
                                                                 + this->side_BGO_thickness
                                                                 + this->suppressor_shell_thickness*2.0)/tan(this->bent_end_angle)
            - (this->suppressor_extension_thickness + this->suppressor_shell_thickness*2.0)
            * sin(this->bent_end_angle) +this->suppShift + this->suppressor_back_radius
            - this->suppressor_forward_radius
            - this->suppressor_shell_thickness/2.0 + dist_from_origin;


    if( !(this->suppressor_position_selector) && include_extension_suppressors )
    {
        // If the detectors are forward, put the extensions in the back position
        // the placement of the extensions matches the placement of the side_suppressor pieces

        x0 += extension_radial_shift ;

        z0 += extension_back_shift ;


        for(i=0; i<4; i++)
        {
            rotateExtension[2*i] = new G4RotationMatrix;
            rotateExtension[2*i]->rotateZ(M_PI/2.0);
            rotateExtension[2*i]->rotateY(this->bent_end_angle);
            rotateExtension[2*i]->rotateX(M_PI/2.0 - M_PI/2.0*i);
            rotateExtension[2*i]->rotateX(alpha);
            rotateExtension[2*i]->rotateY(beta);
            rotateExtension[2*i]->rotateZ(gamma);

            x = x0*cos((i-1)*M_PI/2.0) + y0*sin((i-1)*M_PI/2);
            y = -x0*sin((i-1)*M_PI/2.0) + y0*cos((i-1)*M_PI/2);
            z = z0;

            moveInnerExtension[2*i] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));
            this->rightSuppressorExtensionAssembly->MakeImprint(exp_hall_log, moveInnerExtension[2*i], rotateExtension[2*i], copy_number++, surfCheck);

            rotateExtension[2*i+1] = new G4RotationMatrix;
            rotateExtension[2*i+1]->rotateY(M_PI/2.0);
            rotateExtension[2*i+1]->rotateZ(M_PI/2.0 + this->bent_end_angle);
            rotateExtension[2*i+1]->rotateX(-M_PI/2.0*i);
            rotateExtension[2*i+1]->rotateX(alpha);
            rotateExtension[2*i+1]->rotateY(beta);
            rotateExtension[2*i+1]->rotateZ(gamma);

            x = x0*sin(i*M_PI/2.0) + y0*cos(i*M_PI/2);
            y = x0*cos(i*M_PI/2.0) - y0*sin(i*M_PI/2);
            z = z0;

            moveInnerExtension[2*i+1] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));
            this->leftSuppressorExtensionAssembly->MakeImprint(exp_hall_log, moveInnerExtension[2*i+1], rotateExtension[2*i+1], copy_number_two++, surfCheck);
        }

    }//end if(detectors forward) statement

    // Otherwise, put them forward
    else if( ( this->suppressor_position_selector == 1 ) && include_extension_suppressors)
    {

        for(i=0; i<4; i++)
        {

            rotateExtension[2*i] = new G4RotationMatrix;
            rotateExtension[2*i]->rotateZ(M_PI/2.0);
            rotateExtension[2*i]->rotateY(this->bent_end_angle);
            rotateExtension[2*i]->rotateX(M_PI/2.0-M_PI/2.0*i);
            rotateExtension[2*i]->rotateX(alpha);
            rotateExtension[2*i]->rotateY(beta);
            rotateExtension[2*i]->rotateZ(gamma);

            x = x0*cos((i-1)*M_PI/2.0) + y0*sin((i-1)*M_PI/2);
            y = -x0*sin((i-1)*M_PI/2.0) + y0*cos((i-1)*M_PI/2);
            z = z0;

            moveInnerExtension[2*i] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));
            this->rightSuppressorExtensionAssembly->MakeImprint(exp_hall_log, moveInnerExtension[2*i], rotateExtension[2*i], copy_number++, surfCheck);

            rotateExtension[2*i+1] = new G4RotationMatrix;
            rotateExtension[2*i+1]->rotateY(M_PI/2.0);
            rotateExtension[2*i+1]->rotateZ(M_PI/2.0 + this->bent_end_angle);
            rotateExtension[2*i+1]->rotateX(-M_PI/2.0*i);
            rotateExtension[2*i+1]->rotateX(alpha);
            rotateExtension[2*i+1]->rotateY(beta);
            rotateExtension[2*i+1]->rotateZ(gamma);

            x = x0*sin(i*M_PI/2.0) + y0*cos(i*M_PI/2);
            y = x0*cos(i*M_PI/2.0) - y0*sin(i*M_PI/2);
            z = z0;

            moveInnerExtension[2*i+1] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi),DetectionSystemGriffin::transY(x,y,z,theta,phi),DetectionSystemGriffin::transZ(x,y,z,theta,phi));
            this->leftSuppressorExtensionAssembly->MakeImprint(exp_hall_log, moveInnerExtension[2*i+1], rotateExtension[2*i+1], copy_number_two++, surfCheck);
        }

    }//end if(detectors back) statement
    return 1;
} // End PlaceEverythingButCrystals

G4int DetectionSystemGriffin::PlaceDeadLayerSpecificCrystal(G4LogicalVolume* exp_hall_log, G4int detector_number, G4int position_number, G4bool posTigress)
{

    G4double theta  = this->coords[position_number][0]*deg;
    G4double phi    = this->coords[position_number][1]*deg;
    G4double alpha  = this->coords[position_number][2]*deg; // yaw
    G4double beta   = this->coords[position_number][3]*deg; // pitch
    G4double gamma  = this->coords[position_number][4]*deg; // roll

    //In GRIFFIN upstream is now downstream relative to TIGRESS. However, we want to maintain the same numbering scheme as TIGRESS. Net result is that the lampshade angles change by 45 degrees.
    if(posTigress && (position_number < 4 || position_number > 11) ) {
        phi       = this->coords[position_number][1]*deg - 45.0*deg;
        gamma     = this->coords[position_number][4]*deg - 45.0*deg;
    }

    G4double x;
    G4double y;
    G4double z;

    G4double x0, y0, z0;

    G4int i;

    // positioning
    G4double dist_from_origin_det = this->air_box_back_length_det/2.0 +this->air_box_front_length_det + this->new_rhombi_radius_det;

    x0 = (this->germanium_width + this->germanium_separation)/2.0;
    y0 = (this->germanium_width + this->germanium_separation)/2.0;
    z0 = this->germanium_length/2.0 + this->can_face_thickness/2.0 + this->germanium_dist_from_can_face + this->shift + this->applied_back_shift + dist_from_origin_det;

    /////////////////////////////////////////////////////////////////////
    // now we place all 4 of the 1/4 detectors using the LogicalVolume
    // of the 1st, ensuring that they too will have the hole
    /////////////////////////////////////////////////////////////////////
    copy_number = germanium_copy_number + detector_number*4;

    G4RotationMatrix* rotate_germanium[4];
    G4ThreeVector move_germanium[4];

    for(i=0; i<4; i++){
        rotate_germanium[i] = new G4RotationMatrix;
        rotate_germanium[i]->rotateY(-M_PI/2.0);
        rotate_germanium[i]->rotateX(M_PI/2.0-M_PI/2.0*i);
        rotate_germanium[i]->rotateX(alpha);
        rotate_germanium[i]->rotateY(beta);
        rotate_germanium[i]->rotateZ(gamma);

        x = -x0*pow(-1, floor((i+1)/2) );
        y = -y0*pow(-1, floor((i+2)/2) );
        z = z0;

        move_germanium[i] = G4ThreeVector(DetectionSystemGriffin::transX(x,y,z,theta,phi), DetectionSystemGriffin::transY(x,y,z,theta,phi), DetectionSystemGriffin::transZ(x,y,z,theta,phi));

        this->germaniumAssemblyCry[i]->MakeImprint(exp_hall_log, move_germanium[i], rotate_germanium[i], copy_number++, surfCheck);
    }

    /////////////////////////////////////////////////////////////////////
    // end germanium_block1_log
    /////////////////////////////////////////////////////////////////////

    return 1;
} // end PlaceDeadLayerSpecificCrystal()


///////////////////////////////////////////////////////////////////////
// BuildOneDetector()
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::BuildOneDetector()
{
    // Build assembly volumes
    // Holds all pieces that are not a detector (ie. the can, nitrogen tank, cold finger, electrodes, etc.)
    this->assembly                           = new G4AssemblyVolume();

    // Holds germanium cores
    this->germaniumAssembly                  = new G4AssemblyVolume();

    // Holds left suppressors
    this->leftSuppressorCasingAssembly       = new G4AssemblyVolume();

    // Holds right suppressors
    this->rightSuppressorCasingAssembly      = new G4AssemblyVolume();

    // holds left suppressor extensions
    this->leftSuppressorExtensionAssembly    = new G4AssemblyVolume();

    // holds right suppressor extensions
    this->rightSuppressorExtensionAssembly   = new G4AssemblyVolume();

    // Holds back suppressors
    this->suppressorBackAssembly             = new G4AssemblyVolume();

    // Holds the extension suppressor shells
    this->extensionSuppressorShellAssembly   = new G4AssemblyVolume() ;

    // Holds the back and side suppressor shells
    this->backAndSideSuppressorShellAssembly = new G4AssemblyVolume() ;

    // Holds the Hevimets
    this->hevimetAssembly                    = new G4AssemblyVolume() ;

    ConstructComplexDetectorBlockWithDeadLayer();
    BuildelectrodeMatElectrodes();

    ConstructDetector();

    // Include BGOs?
    if (this->BGO_selector == 1)
    {
        ConstructNewSuppressorCasingWithShells() ;
    }
    else if (this->BGO_selector == 0)
    {
        G4cout << "Not building BGO " << G4endl ;
    }
    else
    {
        G4cout << "Error 234235" << G4endl ;
        exit(1);
    }

    ConstructColdFinger();

    if( this->suppressor_position_selector && this->hevimet_selector )
        ConstructNewHeavyMet();

} // end BuildOneDetector()

void DetectionSystemGriffin::BuildEverythingButCrystals()
{
    // Build assembly volumes
    this->assembly                           = new G4AssemblyVolume();
    this->leftSuppressorCasingAssembly       = new G4AssemblyVolume();
    this->rightSuppressorCasingAssembly      = new G4AssemblyVolume();
    this->leftSuppressorExtensionAssembly    = new G4AssemblyVolume();
    this->rightSuppressorExtensionAssembly   = new G4AssemblyVolume();
    this->suppressorBackAssembly             = new G4AssemblyVolume();

    // Holds the extension suppressor shells
    this->extensionSuppressorShellAssembly   = new G4AssemblyVolume() ;

    // Holds the back and side suppressor shells
    this->backAndSideSuppressorShellAssembly = new G4AssemblyVolume() ;

    // Holds the hevimets
    this->hevimetAssembly                    = new G4AssemblyVolume() ;

    BuildelectrodeMatElectrodes();

    ConstructDetector();

    // Include BGOs?
    if (this->BGO_selector == 1) {
        ConstructNewSuppressorCasingWithShells();
    }
    else if (this->BGO_selector == 0) {
        G4cout << "Not building BGO " << G4endl;
    }
    else {
        G4cout << "Error 234235" << G4endl;
        exit(1);
    }

    ConstructColdFinger();

    if( this->suppressor_position_selector &&  this->hevimet_selector )
        ConstructNewHeavyMet();


} // end BuildEverythingButCrystals()


void DetectionSystemGriffin::BuildDeadLayerSpecificCrystal(G4int det)
{
    this->germaniumAssemblyCry[0] = new G4AssemblyVolume();
    this->germaniumAssemblyCry[1] = new G4AssemblyVolume();
    this->germaniumAssemblyCry[2] = new G4AssemblyVolume();
    this->germaniumAssemblyCry[3] = new G4AssemblyVolume();

    for(G4int i=0; i<4; i++)
        ConstructComplexDetectorBlockWithDetectorSpecificDeadLayer(det,i);

}


///////////////////////////////////////////////////////////////////////
// ConstructComplexDetectorBlock builds four quarters of germanium
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructComplexDetectorBlock()
{
    G4Material* materialGe = G4Material::GetMaterial("G4_Ge");
    if( !materialGe ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }
    G4Material* materialVacuum = G4Material::GetMaterial("Vacuum");
    if( !materialVacuum ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }

    G4VisAttributes* germanium_block1_vis_att = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    germanium_block1_vis_att->SetVisibility(true);

    G4VisAttributes* dead_layer_vis_att = new G4VisAttributes(G4Colour(0.80, 0.0, 0.0));
    dead_layer_vis_att->SetVisibility(true);

    G4VisAttributes* air_vis_att = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    air_vis_att->SetVisibility(true);

    G4VisAttributes* electrodeMat_vis_att = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    electrodeMat_vis_att->SetVisibility(true);

    G4SubtractionSolid* detector1 = this->quarterDetector();

    germanium_block1_log = new G4LogicalVolume(detector1, materialGe, "germanium_block1_log", 0, 0, 0);
    germanium_block1_log->SetVisAttributes(germanium_block1_vis_att);

    this->germaniumAssembly->AddPlacedVolume(germanium_block1_log, this->move_null, this->rotate_null);

    // now we make the hole that will go in the back of each quarter detector,
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;
    G4double hole_radius = this->germanium_hole_radius;
    G4double hole_half_length_z = (this->germanium_length -this->germanium_hole_dist_from_face)/2.0;


    G4Tubs* hole_tubs = new G4Tubs("hole_tubs", 0.0, hole_radius, hole_half_length_z, start_angle, final_angle);

    G4ThreeVector move_hole(this->germanium_shift, -(this->germanium_shift),
                            -((this->germanium_hole_dist_from_face)));

    germanium_hole_log = new G4LogicalVolume(hole_tubs, materialVacuum, "germanium_hole_log", 0, 0, 0);
    germanium_hole_log->SetVisAttributes(air_vis_att);

    this->germaniumAssembly->AddPlacedVolume(germanium_hole_log, move_hole, this->rotate_null);

}//end ::ConstructComplexDetectorBlock

///////////////////////////////////////////////////////////////////////
// ConstructComplexDetectorBlockWithDeadLayer builds four quarters of 
// germanium, with dead layers
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructComplexDetectorBlockWithDeadLayer()
{
    G4Material* materialGe = G4Material::GetMaterial("G4_Ge");
    if( !materialGe ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }
    G4Material* materialVacuum = G4Material::GetMaterial("Vacuum");
    if( !materialVacuum ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }

    G4VisAttributes* germanium_block1_vis_att = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    germanium_block1_vis_att->SetVisibility(true);

    G4VisAttributes* dead_layer_vis_att = new G4VisAttributes(G4Colour(0.80, 0.0, 0.0));
    dead_layer_vis_att->SetVisibility(true);

    G4VisAttributes* air_vis_att = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    air_vis_att->SetVisibility(true);

    G4VisAttributes* electrodeMat_vis_att = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    electrodeMat_vis_att->SetVisibility(true);

    G4SubtractionSolid* detector1 = this->quarterDetector();

    germanium_block1_log = new G4LogicalVolume(detector1, materialGe, "germanium_block1_log", 0, 0, 0);
    germanium_block1_log->SetVisAttributes(germanium_block1_vis_att);

    this->germaniumAssembly->AddPlacedVolume(germanium_block1_log, this->move_null, this->rotate_null);

    /////////////////////////////////////////////////////////////////////
    // Since if we are using this method, the inner dead layer must be
    // simulated, make a cylinder of germanium, place the air cylinder
    // inside it, and then place the whole thing in the crystal
    /////////////////////////////////////////////////////////////////////
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;
    G4double dead_layer_radius = this->germanium_hole_radius + this->inner_dead_layer_thickness;
    G4double dead_layer_half_length_z = (this->germanium_length
                                         - this->germanium_hole_dist_from_face
                                         + this->inner_dead_layer_thickness)/2.0;

    G4Tubs* dead_layer_tubs = new G4Tubs("dead_layer_tubs", this->germanium_hole_radius,
                                         dead_layer_radius, dead_layer_half_length_z, start_angle,final_angle);

    G4ThreeVector move_dead_layer(this->germanium_shift, -(this->germanium_shift),
                                  -((this->germanium_hole_dist_from_face
                                     - this->inner_dead_layer_thickness)/2.0));

    inner_dead_layer_log = new G4LogicalVolume(dead_layer_tubs, materialGe, "inner_dead_layer_log", 0, 0, 0);
    inner_dead_layer_log->SetVisAttributes(dead_layer_vis_att);

    this->germaniumAssembly->AddPlacedVolume(inner_dead_layer_log, move_dead_layer, this->rotate_null);

    // dead layer cap
    dead_layer_radius = this->germanium_hole_radius;
    G4double dead_layer_cap_half_length_z = (this->inner_dead_layer_thickness)/2.0;

    G4Tubs* dead_layer_cap = new G4Tubs("dead_layer_cap", 0.0,
                                        dead_layer_radius, dead_layer_cap_half_length_z, start_angle,final_angle);

    G4ThreeVector move_dead_layer_cap(this->germanium_shift, -(this->germanium_shift),
                                      -((this->germanium_hole_dist_from_face - this->inner_dead_layer_thickness)/2.0 - dead_layer_half_length_z + this->inner_dead_layer_thickness/2.0));

    inner_dead_layer_cap_log = new G4LogicalVolume(dead_layer_cap, materialGe, "inner_dead_layer_cap_log", 0, 0, 0);
    inner_dead_layer_cap_log->SetVisAttributes(dead_layer_vis_att);

    this->germaniumAssembly->AddPlacedVolume(inner_dead_layer_cap_log, move_dead_layer_cap, this->rotate_null);


    // now we make the hole that will go in the back of each quarter detector,
    // and place it inside the dead layer cylinder
    start_angle = 0.0*M_PI;
    final_angle = 2.0*M_PI;
    G4double hole_radius = this->germanium_hole_radius;
    G4double hole_half_length_z = (this->germanium_length -this->germanium_hole_dist_from_face)/2.0;

    G4Tubs* hole_tubs = new G4Tubs("hole_tubs", 0.0, hole_radius, hole_half_length_z, start_angle, final_angle);

    G4ThreeVector move_hole(this->germanium_shift, -(this->germanium_shift),
                            -((this->germanium_hole_dist_from_face
                               - this->inner_dead_layer_thickness)/2.0)
                            -(this->inner_dead_layer_thickness/2.0));

    germanium_hole_log = new G4LogicalVolume(hole_tubs, materialVacuum, "germanium_hole_log", 0, 0, 0);
    germanium_hole_log->SetVisAttributes(air_vis_att);

    this->germaniumAssembly->AddPlacedVolume(germanium_hole_log, move_hole, this->rotate_null);

}//end ::ConstructComplexDetectorBlockWithDeadLayer

void DetectionSystemGriffin::ConstructComplexDetectorBlockWithDetectorSpecificDeadLayer(G4int det, G4int cry)
{
    G4String strdet = G4UIcommand::ConvertToString(det);
    G4String strcry = G4UIcommand::ConvertToString(cry);

    G4Material* materialGe = G4Material::GetMaterial("G4_Ge");
    if( !materialGe ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }
    G4Material* materialVacuum = G4Material::GetMaterial("Vacuum");
    if( !materialVacuum ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }

    G4VisAttributes* germanium_block1_vis_att = new G4VisAttributes(this->griffinCrystalColours[cry]);
    germanium_block1_vis_att->SetVisibility(true);

    G4VisAttributes* dead_layer_vis_att = new G4VisAttributes(this->griffinDeadLayerColours[cry]);
    dead_layer_vis_att->SetVisibility(true);

    G4VisAttributes* air_vis_att = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    air_vis_att->SetVisibility(true);

    G4VisAttributes* electrodeMat_vis_att = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    electrodeMat_vis_att->SetVisibility(true);

    G4SubtractionSolid* detector1 = this->quarterSpecificDeadLayerDetector(det, cry);

    G4String germanium_block1_name = "germanium_dls_block1_" + strdet + "_" + strcry + "_log" ;

    germanium_block1_log = new G4LogicalVolume(detector1, materialGe, germanium_block1_name, 0, 0, 0);
    germanium_block1_log->SetVisAttributes(germanium_block1_vis_att);

    this->germaniumAssemblyCry[cry]->AddPlacedVolume(germanium_block1_log, this->move_null, this->rotate_null);

    /////////////////////////////////////////////////////////////////////
    // Since if we are using this method, the inner dead layer must be
    // simulated, make a cylinder of germanium, place the air cylinder
    // inside it, and then place the whole thing in the crystal
    /////////////////////////////////////////////////////////////////////
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;
    G4double dead_layer_radius = this->germanium_hole_radius + this->griffinDeadLayers[det][cry];
    G4double dead_layer_half_length_z = (this->germanium_length
                                         - this->germanium_hole_dist_from_face
                                         + this->griffinDeadLayers[det][cry])/2.0;

    G4Tubs* dead_layer_tubs = new G4Tubs("dead_layer_tubs", this->germanium_hole_radius,
                                         dead_layer_radius, dead_layer_half_length_z, start_angle,final_angle);

    G4ThreeVector move_dead_layer(this->germanium_shift, -(this->germanium_shift),
                                  -((this->germanium_hole_dist_from_face
                                     - this->griffinDeadLayers[det][cry])/2.0));

    inner_dead_layer_log = new G4LogicalVolume(dead_layer_tubs, materialGe, "inner_dead_layer_log", 0, 0, 0);
    inner_dead_layer_log->SetVisAttributes(dead_layer_vis_att);

    this->germaniumAssemblyCry[cry]->AddPlacedVolume(inner_dead_layer_log, move_dead_layer, this->rotate_null);

    // dead layer cap
    dead_layer_radius = this->germanium_hole_radius;
    G4double dead_layer_cap_half_length_z = (this->griffinDeadLayers[det][cry])/2.0;

    G4Tubs* dead_layer_cap = new G4Tubs("dead_layer_cap", 0.0,
                                        dead_layer_radius, dead_layer_cap_half_length_z, start_angle,final_angle);

    G4ThreeVector move_dead_layer_cap(this->germanium_shift, -(this->germanium_shift),
                                      -((this->germanium_hole_dist_from_face - this->griffinDeadLayers[det][cry])/2.0 - dead_layer_half_length_z + this->griffinDeadLayers[det][cry]/2.0));

    inner_dead_layer_cap_log = new G4LogicalVolume(dead_layer_cap, materialGe, "inner_dead_layer_cap_log", 0, 0, 0);
    inner_dead_layer_cap_log->SetVisAttributes(dead_layer_vis_att);

    this->germaniumAssemblyCry[cry]->AddPlacedVolume(inner_dead_layer_cap_log, move_dead_layer_cap, this->rotate_null);


    // now we make the hole that will go in the back of each quarter detector,
    // and place it inside the dead layer cylinder
    start_angle = 0.0*M_PI;
    final_angle = 2.0*M_PI;
    G4double hole_radius = this->germanium_hole_radius;
    G4double hole_half_length_z = (this->germanium_length -this->germanium_hole_dist_from_face)/2.0;

    G4Tubs* hole_tubs = new G4Tubs("hole_tubs", 0.0, hole_radius, hole_half_length_z, start_angle, final_angle);

    G4ThreeVector move_hole(this->germanium_shift, -(this->germanium_shift),
                            -((this->germanium_hole_dist_from_face
                               - this->griffinDeadLayers[det][cry])/2.0)
                            -(this->griffinDeadLayers[det][cry]/2.0));

    germanium_hole_log = new G4LogicalVolume(hole_tubs, materialVacuum, "germanium_hole_log", 0, 0, 0);
    germanium_hole_log->SetVisAttributes(air_vis_att);

    this->germaniumAssemblyCry[cry]->AddPlacedVolume(germanium_hole_log, move_hole, this->rotate_null);

}//end ::ConstructComplexDetectorBlockWithDetectorSpecificDeadLayer

///////////////////////////////////////////////////////////////////////
// Builds a layer of electrodeMat between germanium crystals to
// approximate electrodes, etc. that Eurisys won't tell us 
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::BuildelectrodeMatElectrodes()
{

    G4Material* electrodeMaterial = G4Material::GetMaterial(this->electrodeMaterial);
    if( !electrodeMaterial ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }

    G4VisAttributes* electrodeMat_vis_att = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    electrodeMat_vis_att->SetVisibility(true);

    // electrodeMat layers between crystals - back part
    G4RotationMatrix* rotate_electrodeMat_back = new G4RotationMatrix;
    rotate_electrodeMat_back->rotateZ(M_PI/2.0);
    G4ThreeVector move_inter_crystal_electrodeMat_back(this->germanium_length/2.0
                                                       + this->germanium_bent_length/2.0
                                                       + this->can_face_thickness/2.0
                                                       + this->germanium_dist_from_can_face + this->shift
                                                       + this->applied_back_shift, 0,0);

    G4UnionSolid* inter_crystal_electrodeMat_back = this->interCrystalelectrodeMatBack();

    inter_crystal_electrodeMat_back_log = new G4LogicalVolume(inter_crystal_electrodeMat_back,
                                                              electrodeMaterial, "inter_crystal_electrodeMat_back_log", 0, 0, 0);
    inter_crystal_electrodeMat_back_log->SetVisAttributes(electrodeMat_vis_att);

    this->assembly->AddPlacedVolume(inter_crystal_electrodeMat_back_log, move_inter_crystal_electrodeMat_back, rotate_electrodeMat_back);

    // electrodeMat between crystals - front part
    G4RotationMatrix* rotate_electrodeMat_front = new G4RotationMatrix;
    rotate_electrodeMat_front->rotateY(M_PI/2.0);
    G4UnionSolid* inter_crystal_electrodeMat_front = this->interCrystalelectrodeMatFront();

    G4ThreeVector move_inter_crystal_electrodeMat_front(this->germanium_bent_length/2.0
                                                        + this->electrodeMat_starting_depth/2.0
                                                        + this->can_face_thickness/2.0
                                                        + this->germanium_dist_from_can_face + this->shift
                                                        + this->applied_back_shift, 0,0);

    inter_crystal_electrodeMat_front_log = new G4LogicalVolume(inter_crystal_electrodeMat_front,
                                                               electrodeMaterial, "inter_crystal_electrodeMat_front_log", 0, 0, 0);
    inter_crystal_electrodeMat_front_log->SetVisAttributes(electrodeMat_vis_att);

    this->assembly->AddPlacedVolume(inter_crystal_electrodeMat_front_log, move_inter_crystal_electrodeMat_front, rotate_electrodeMat_front);
}//end ::BuildelectrodeMatElectrodes

///////////////////////////////////////////////////////////////////////
// ConstructDetector builds the structureMat can, cold finger shell, 
// and liquid nitrogen tank. Since the can's face is built first, 
// everything else is placed relative to it  
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructDetector() 
{

    G4double x0, y0, z0, yMov, zMov;
    G4int i;
    G4RotationMatrix* rotate_piece[4] ;
    G4ThreeVector move_piece[4] ;

    G4Material* structureMaterial = G4Material::GetMaterial(this->structureMaterial);
    if( !structureMaterial ) {
        G4cout << " ----> Material " << this->structureMaterial << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }
    G4Material* liquidN2Material = G4Material::GetMaterial("Liquid_N2");
    if( !liquidN2Material ) {
        G4cout << " ----> Material " << "Liquid_N2" << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }
    // first we make the can's front face
    G4Box* front_face = this->squareFrontFace();
    front_face_log = new G4LogicalVolume(front_face, structureMaterial, "front_face_log", 0, 0, 0);

    G4ThreeVector move_front_face(this->shift +this->applied_back_shift, 0, 0);

    this->assembly->AddPlacedVolume(front_face_log, move_front_face, this->rotate_null);

    // now we put on the four angled side pieces

    x0 = ((this->bent_end_length)/2.0) + this->shift + this->applied_back_shift ;

    z0 = ((this->can_face_thickness - this->bent_end_length) * tan(this->bent_end_angle)
          + this->detector_total_width - this->can_face_thickness)/2.0 ;

    G4Para* side_piece[4] ;
    
    for( i = 0 ; i < 4 ; i++ )
    {
        // Top, Right, Bottom, Left
        side_piece[i] = this->bentSidePiece() ;
        rotate_piece[i] = new G4RotationMatrix ;

        // The left side is slightly different than the others, so this if statement is required. The order of rotation matters.
        if(i == 3)
        {
            rotate_piece[i]->rotateX( -M_PI/2.0 + i*M_PI/2.0 ) ;
            rotate_piece[i]->rotateY(this->bent_end_angle) ;
        }
        else
        {
            rotate_piece[i]->rotateY( -this->bent_end_angle) ;
            rotate_piece[i]->rotateX( -M_PI/2.0 + i*M_PI/2.0 ) ;
        }

        yMov = z0 * cos( i * M_PI/2.0 ) ;
        zMov = z0 * sin( i * M_PI/2.0 ) ;

        move_piece[i] = G4ThreeVector(x0, yMov, zMov) ;
    }

    top_bent_piece_log = new G4LogicalVolume(side_piece[0], structureMaterial, "top_bent_piece", 0, 0, 0);
    this->assembly->AddPlacedVolume(top_bent_piece_log, move_piece[0], rotate_piece[0]);

    right_bent_piece_log = new G4LogicalVolume(side_piece[1], structureMaterial, "right_bent_piece", 0, 0, 0);
    this->assembly->AddPlacedVolume(right_bent_piece_log, move_piece[1], rotate_piece[1]);

    bottom_bent_piece_log = new G4LogicalVolume(side_piece[2], structureMaterial, "right_bent_piece", 0, 0, 0);
    this->assembly->AddPlacedVolume(right_bent_piece_log, move_piece[2], rotate_piece[2]);

    left_bent_piece_log = new G4LogicalVolume(side_piece[3], structureMaterial, "right_bent_piece", 0, 0, 0);
    this->assembly->AddPlacedVolume(right_bent_piece_log, move_piece[3], rotate_piece[3]);


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // now we add the wedges to the edges of the face. These complete the angled side pieces
    G4Trap* side_wedge[4] ;

    x0 = this->shift +this->applied_back_shift ;

    y0 = (this->detector_total_width/2.0) - (this->bent_end_length * tan(this->bent_end_angle))
            + (this->can_face_thickness/2.0) *tan(this->bent_end_angle)/2.0 ;

    for( i = 0 ; i < 4 ; i++ )
    {
        // Top, Right, Bottom, Left
        side_wedge[i] = this->cornerWedge() ;
        rotate_piece[i] = new G4RotationMatrix ;

        rotate_piece[i]->rotateX( M_PI/2.0 ) ;
        rotate_piece[i]->rotateY( -M_PI/2.0 ) ;
        rotate_piece[i]->rotateX( -M_PI/2.0 + i*M_PI/2.0 ) ;

        yMov = y0 * cos( i*M_PI/2.0 ) ;
        zMov = y0 * sin( i*M_PI/2.0 ) ;

        move_piece[i] = G4ThreeVector(x0, yMov, zMov) ;
    }

    top_wedge_log = new G4LogicalVolume(side_wedge[0], structureMaterial, "top_wedge_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(top_wedge_log, move_piece[0], rotate_piece[0]);

    right_wedge_log = new G4LogicalVolume(side_wedge[1], structureMaterial, "right_wedge_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(right_wedge_log, move_piece[1], rotate_piece[1]);

    bottom_wedge_log = new G4LogicalVolume(side_wedge[2], structureMaterial, "bottom_wedge_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(bottom_wedge_log, move_piece[2], rotate_piece[2]);

    left_wedge_log = new G4LogicalVolume(side_wedge[3], structureMaterial, "left_wedge_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(left_wedge_log, move_piece[3], rotate_piece[3]);


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // now we add the rounded corners, made from a quarter of a cone

    // Correction factor used to move the corner cones so as to avoid
    // holes left in the can caused by the bent side pieces

    G4double hole_eliminator = this->can_face_thickness *(1.0 -tan(this->bent_end_angle));
    G4Cons* cone_location[4] ;

    x0 = this->bent_end_length/2.0 + this->shift + this->applied_back_shift ;

    y0 = (this->detector_total_width/2.0) - (this->bent_end_length * tan(this->bent_end_angle)) - hole_eliminator ;

    for( i = 0 ; i < 4 ; i++ )
    {
        // Lower Left, Upper Left, Upper Right, Lower Right
        cone_location[i] = this->roundedEndEdge() ;
        rotate_piece[i] = new G4RotationMatrix ;

        rotate_piece[i]->rotateY( M_PI/2.0 ) ;
        rotate_piece[i]->rotateX( -M_PI/2.0 + i*M_PI/2.0 ) ;

        yMov = y0 * ( -cos(i*M_PI/2.0) + sin(i*M_PI/2.0) ) ;
        zMov = y0 * ( -cos(i*M_PI/2.0) - sin(i*M_PI/2.0) ) ;

        move_piece[i] = G4ThreeVector(x0, yMov, zMov) ;
    }

    lower_left_cone_log = new G4LogicalVolume(cone_location[0], structureMaterial, "lower_left_cone_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(lower_left_cone_log, move_piece[0],  rotate_piece[0]);

    upper_left_cone_log = new G4LogicalVolume(cone_location[1], structureMaterial, "upper_left_cone_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(upper_left_cone_log, move_piece[1], rotate_piece[1]);

    upper_right_cone_log = new G4LogicalVolume(cone_location[2], structureMaterial, "upper_right_cone_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(upper_right_cone_log, move_piece[2], rotate_piece[2]);

    lower_right_cone_log = new G4LogicalVolume(cone_location[3], structureMaterial, "lower_right_cone_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(lower_right_cone_log, move_piece[3], rotate_piece[3]);


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /* now we add the corner tubes which extend the rounded corners to the back of the can
  *
  * This is extremely similar to the previous loop but they are not exactly the same. While it would be
  * possible to consolidate them into one, they would only share the rotation variables, so it
  * might not be a huge benefit.
  */

    G4Tubs* tube_location[4] ;

    x0 = (this->detector_total_length +this->bent_end_length
          - this->rear_plate_thickness -this->can_face_thickness)/2.0
            + this->shift +this->applied_back_shift ;

    y0 = (this->detector_total_width/2.0) - (this->bent_end_length * tan(this->bent_end_angle)) ;

    for( i = 0 ; i < 4 ; i++ )
    {
        // Lower Left, Upper Left, Upper Right, Lower Right
        tube_location[i] = this->cornerTube() ;
        rotate_piece[i] = new G4RotationMatrix ;
        rotate_piece[i]->rotateY(M_PI/2.0) ;
        rotate_piece[i]->rotateX( -M_PI/2.0 + i*M_PI/2.0 ) ;

        yMov = y0 * ( -cos(i*M_PI/2.0) + sin(i*M_PI/2.0) ) ;
        zMov = y0 * ( -cos(i*M_PI/2.0) - sin(i*M_PI/2.0) ) ;

        move_piece[i] = G4ThreeVector(x0, yMov, zMov) ;

    }

    lower_left_tube_log = new G4LogicalVolume(tube_location[0], structureMaterial, "lower_left_tube_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(lower_left_tube_log, move_piece[0], rotate_piece[0]);

    upper_left_tube_log = new G4LogicalVolume(tube_location[1], structureMaterial, "upper_left_tube_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(upper_left_tube_log, move_piece[1], rotate_piece[1]);

    upper_right_tube_log = new G4LogicalVolume(tube_location[2], structureMaterial, "upper_right_tube_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(upper_right_tube_log, move_piece[2], rotate_piece[2]);

    lower_right_tube_log = new G4LogicalVolume(tube_location[3], structureMaterial, "lower_right_tube_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(lower_right_tube_log, move_piece[3], rotate_piece[3]);


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // now we add the side panels that extend from the bent pieces to the back of the can


    G4Box* panel_location[4] ;

    x0 = (this->detector_total_length + this->bent_end_length - this->rear_plate_thickness - this->can_face_thickness)/2.0
            + this->shift + this->applied_back_shift ;

    y0 = 0 ;

    z0 = (this->detector_total_width - this->can_side_thickness)/2.0 ;

    for( i = 0 ; i < 4 ; i++ )
    {
        // Right, Top, Left, Bottom
        panel_location[i] = this->sidePanel() ;
        rotate_piece[i] = new G4RotationMatrix ;
        rotate_piece[i]->rotateX( M_PI/2.0 - i*M_PI/2.0 ) ;

        yMov = z0*sin( i*M_PI/2.0 ) ;
        zMov = z0*cos( i*M_PI/2.0 ) ;

        move_piece[i] = G4ThreeVector( x0, yMov, zMov ) ;

    }

    right_side_panel_log = new G4LogicalVolume(panel_location[0], structureMaterial, "right_side_panel_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(right_side_panel_log, move_piece[0], rotate_piece[0]);

    top_side_panel_log = new G4LogicalVolume(panel_location[1], structureMaterial, "top_side_panel_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(top_side_panel_log, move_piece[1], this->rotate_null);

    left_side_panel_log = new G4LogicalVolume(panel_location[2], structureMaterial, "left_side_panel_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(left_side_panel_log, move_piece[2], rotate_piece[2]);

    bottom_side_panel_log = new G4LogicalVolume(panel_location[3], structureMaterial, "bottom_side_panel_log", 0, 0, 0);
    this->assembly->AddPlacedVolume(bottom_side_panel_log, move_piece[3], rotate_piece[3]);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // now we add the rear plate, which has a hole for the cold finger to pass through
    G4SubtractionSolid* rear_plate = this->rearPlate();

    G4RotationMatrix* rotate_rear_plate = new G4RotationMatrix;
    rotate_rear_plate->rotateY(M_PI/2.0);

    x0 = this->detector_total_length -this->can_face_thickness/2.0
            - this->rear_plate_thickness/2.0 +this->shift
            + this->applied_back_shift ;

    G4ThreeVector move_rear_plate(this->detector_total_length -this->can_face_thickness/2.0
                                  - this->rear_plate_thickness/2.0 +this->shift
                                  + this->applied_back_shift, 0, 0);

    rear_plate_log = new G4LogicalVolume(rear_plate, structureMaterial, "rear_plate_log", 0, 0, 0);

    this->assembly->AddPlacedVolume(rear_plate_log, move_rear_plate, rotate_rear_plate);

    // we know add the cold finger shell, which extends out the back of the can
    G4Tubs* finger_shell = this->coldFingerShell();

    G4RotationMatrix* rotate_finger_shell = new G4RotationMatrix;
    rotate_finger_shell->rotateY(M_PI/2.0);

    G4ThreeVector move_finger_shell(this->cold_finger_shell_length/2.0 -this->can_face_thickness/2.0
                                    + this->detector_total_length +this->shift
                                    + this->applied_back_shift, 0, 0);

    finger_shell_log = new G4LogicalVolume(finger_shell, structureMaterial, "finger_shell_log", 0, 0, 0);

    this->assembly->AddPlacedVolume(finger_shell_log, move_finger_shell, rotate_finger_shell);

    // lastly we add the liquid nitrogen tank at the back, past the cold finger shell
    G4Tubs* tank = this->liquidNitrogenTank();

    G4RotationMatrix* rotate_tank = new G4RotationMatrix;
    rotate_tank->rotateY(M_PI/2.0);

    G4ThreeVector move_tank((this->coolant_length -this->can_face_thickness)/2.0 + this->detector_total_length
                            + this->cold_finger_shell_length + this->shift +this->applied_back_shift, 0, 0);

    tank_log = new G4LogicalVolume(tank, structureMaterial, "tank_log", 0, 0, 0);

    this->assembly->AddPlacedVolume(tank_log, move_tank, rotate_tank);

    // lids
    G4Tubs* lid = this->liquidNitrogenTankLid();

    G4RotationMatrix* rotate_lid = new G4RotationMatrix;
    rotate_lid->rotateY(M_PI/2.0);

    G4ThreeVector move_lid1( this->coolant_length/2.0 - this->coolant_thickness/2.0 + (this->coolant_length -this->can_face_thickness)/2.0 + this->detector_total_length
                            + this->cold_finger_shell_length + this->shift +this->applied_back_shift, 0, 0);
    G4ThreeVector move_lid2( -1.0*this->coolant_length/2.0 + this->coolant_thickness/2.0 + (this->coolant_length -this->can_face_thickness)/2.0 + this->detector_total_length
                            + this->cold_finger_shell_length + this->shift +this->applied_back_shift, 0, 0);

    tank_lid1_log = new G4LogicalVolume(lid, structureMaterial, "tank_lid1_log", 0, 0, 0);
    tank_lid2_log = new G4LogicalVolume(lid, structureMaterial, "tank_lid2_log", 0, 0, 0);

    this->assembly->AddPlacedVolume(tank_lid1_log, move_lid1, rotate_lid);
    this->assembly->AddPlacedVolume(tank_lid2_log, move_lid2, rotate_lid);

    // LN2
    G4Tubs* liquidn2 = this->liquidNitrogen();

    G4RotationMatrix* rotate_liquidn2 = new G4RotationMatrix;
    rotate_liquidn2->rotateY(M_PI/2.0);

    G4ThreeVector move_liquidn2((this->coolant_length -this->can_face_thickness)/2.0 + this->detector_total_length
                            + this->cold_finger_shell_length + this->shift +this->applied_back_shift, 0, 0);

    tank_liquid_log = new G4LogicalVolume(liquidn2, liquidN2Material, "tank_lid1_log", 0, 0, 0);

    this->assembly->AddPlacedVolume(tank_liquid_log, move_liquidn2, rotate_liquidn2);
}// end::ConstructDetector

///////////////////////////////////////////////////////////////////////
// ConstructColdFinger has changed in this version!  It now 
// incorporates all the internal details released in Nov 2004 by 
// Eurisys. ConstructColdFinger builds the cold finger as well as the 
// cold plate. The finger extends to the Liquid Nitrogen tank, while 
// the plate is always the same distance from the back of the germanium  
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructColdFinger()
{

    G4Material* materialAir = G4Material::GetMaterial("Air");
    if( !materialAir ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }
    G4Material* structureMaterial = G4Material::GetMaterial(this->structureMaterial);
    if( !structureMaterial ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }
    G4Material* electrodeMaterial = G4Material::GetMaterial(this->electrodeMaterial);
    if( !electrodeMaterial ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }

    G4VisAttributes* cold_finger_vis_att = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    cold_finger_vis_att->SetVisibility(true);

    G4VisAttributes* structureMat_vis_att = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    structureMat_vis_att->SetVisibility(true);

    G4VisAttributes* air_vis_att = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    air_vis_att->SetVisibility(true);

    G4Box* end_plate = this->endPlate();

    // cut out holes
    G4double air_hole_distance = this->germanium_separation/2.0
            + (this->germanium_width/2.0 - this->germanium_shift); //centre of the middle hole

    G4RotationMatrix* rotate_fet_air_hole = new G4RotationMatrix;
    rotate_fet_air_hole->rotateY(M_PI/2.0);

    G4Tubs* fet_air_hole_cut = this->airHoleCut();

    G4ThreeVector move_piece[8] ;
    G4SubtractionSolid* end_plate_cut[4] ;
    G4int i ;
    G4double x, y, z, y0, z0, y01, y02;

    y0 = air_hole_distance ;
    z0 = air_hole_distance ;

    for( i = 0 ; i < 4 ; i++ )
    {
        y = y0 * ( -cos( i*M_PI/2.0 ) + sin( i*M_PI/2.0 ) ) ;
        z = z0 * ( cos( i*M_PI/2.0 ) + sin( i*M_PI/2.0 ) ) ;
        move_piece[i] = G4ThreeVector( 0 , y, z ) ;

        if( i == 0 )
            // Slightly different than the rest. Sorry I couldnt think of a better workaround
            end_plate_cut[i] = new G4SubtractionSolid("end_plate_cut", end_plate, fet_air_hole_cut, rotate_fet_air_hole, move_piece[i] ) ;
        else
            end_plate_cut[i] = new G4SubtractionSolid("end_plate_cut", end_plate_cut[i-1], fet_air_hole_cut, rotate_fet_air_hole, move_piece[i] ) ;
    }

    x = this->cold_finger_end_plate_thickness/2.0 +this->can_face_thickness/2.0
            + this->germanium_dist_from_can_face +this->germanium_length
            + this->cold_finger_space +this->shift +this->applied_back_shift ;
    
    G4ThreeVector move_end_plate( x, 0, 0 ) ;



    // Fix overlap.

    G4double cut_tube_x0, cut_tube_y0, cut_yMov, cut_zMov;
    G4RotationMatrix* cut_rotate_piece[4] ;
    G4ThreeVector cut_move_piece[4] ;
    G4Tubs* cut_tube_location[4] ;

    cut_tube_x0 = (this->detector_total_length +this->bent_end_length
          - this->rear_plate_thickness -this->can_face_thickness)/2.0
            + this->shift +this->applied_back_shift ;

    cut_tube_y0 = (this->detector_total_width/2.0) - (this->bent_end_length * tan(this->bent_end_angle)) ;

    for( i = 0 ; i < 4 ; i++ )
    {
        // Lower Left, Upper Left, Upper Right, Lower Right
        cut_tube_location[i] = this->cornerTube() ;
        cut_rotate_piece[i] = new G4RotationMatrix ;
        cut_rotate_piece[i]->rotateY(M_PI/2.0) ;
        cut_rotate_piece[i]->rotateZ( -M_PI/2.0 + (i-1)*M_PI/2.0 ) ;

        //cut_rotate_piece[i]->rotateX( -M_PI/2.0 + i*M_PI/2.0 ) ;

        cut_yMov = cut_tube_y0 * ( -cos(i*M_PI/2.0) + sin(i*M_PI/2.0) ) ;
        cut_zMov = cut_tube_y0 * ( -cos(i*M_PI/2.0) - sin(i*M_PI/2.0) ) ;

        cut_move_piece[i] = G4ThreeVector(0, cut_yMov, cut_zMov) ;


        end_plate_cut[i+4] = new G4SubtractionSolid("end_plate_cut", end_plate_cut[i+3], cut_tube_location[i], cut_rotate_piece[i], cut_move_piece[i] ) ;

    }












    end_plate_log = new G4LogicalVolume(end_plate_cut[7], structureMaterial, "end_plate_log", 0, 0, 0);

    end_plate_log->SetVisAttributes(structureMat_vis_att);

    this->assembly->AddPlacedVolume(end_plate_log, move_end_plate, 0);

    // Place holes in the old cold plate
    G4Tubs* fet_air_hole = this->airHole();

    fet_air_hole_log = new G4LogicalVolume(fet_air_hole, materialAir, "fet_air_hole_log", 0, 0, 0);
    fet_air_hole_log->SetVisAttributes(air_vis_att);

    for( i = 0 ; i < 4 ; i++ )
    {
        y = y0 * ( -cos( i*M_PI/2.0 ) + sin( i*M_PI/2.0 ) ) ;
        z = z0 * ( cos( i*M_PI/2.0 ) + sin( i*M_PI/2.0 ) ) ;
        move_piece[i] = G4ThreeVector( x, y, z ) ;
        this->assembly->AddPlacedVolume( fet_air_hole_log, move_piece[i], rotate_fet_air_hole ) ;
    }

    // Cold Finger
    G4Tubs* finger = this->finger();

    G4RotationMatrix* rotate_finger = new G4RotationMatrix;
    rotate_finger->rotateY(M_PI/2.0);

    G4ThreeVector move_finger((this->cold_finger_length +this->can_face_thickness)/2.0
                              + this->germanium_dist_from_can_face +this->germanium_length
                              + this->cold_finger_space +this->cold_finger_end_plate_thickness
                              + this->shift +this->applied_back_shift
                              + this->structureMat_cold_finger_thickness/2.0 , 0, 0); //changed Jan 2005

    finger_log = new G4LogicalVolume(finger, electrodeMaterial, "finger_log", 0, 0, 0);
    finger_log->SetVisAttributes(cold_finger_vis_att);

    this->assembly->AddPlacedVolume(finger_log, move_finger, rotate_finger);

    // The extra block of structureMat that's in the diagram "Cut C"
    G4SubtractionSolid* extra_cold_block = this->extraColdBlock();
    G4RotationMatrix* rotate_extra_cold_block = new G4RotationMatrix;
    rotate_extra_cold_block->rotateY(M_PI/2.0);

    G4ThreeVector move_extra_cold_block(this->detector_total_length -this->can_face_thickness/2.0
                                        - this->rear_plate_thickness +this->shift
                                        + this->applied_back_shift -this->extra_block_distance_from_back_plate
                                        - this->extra_block_thickness/2.0, 0, 0);

    extra_cold_block_log = new G4LogicalVolume(extra_cold_block, structureMaterial, "extra_cold_block_log", 0, 0, 0);

    this->assembly->AddPlacedVolume(extra_cold_block_log, move_extra_cold_block, rotate_extra_cold_block);

    // The structures above the cooling end plate
    G4Tubs* structureMat_cold_finger = this->structureMatColdFinger();
    structureMat_cold_finger_log = new G4LogicalVolume(structureMat_cold_finger, structureMaterial, "structureMat_cold_finger_log", 0, 0, 0);
    structureMat_cold_finger_log->SetVisAttributes(structureMat_vis_att);

    G4RotationMatrix* rotate_structureMat_cold_finger = new G4RotationMatrix;
    rotate_structureMat_cold_finger->rotateY(M_PI/2.0);
    G4ThreeVector move_structureMat_cold_finger(this->structureMat_cold_finger_thickness/2.0
                                                + this->cold_finger_end_plate_thickness
                                                + this->can_face_thickness/2.0
                                                + this->germanium_dist_from_can_face +this->germanium_length
                                                + this->cold_finger_space +this->shift +this->applied_back_shift, 0, 0);

    this->assembly->AddPlacedVolume(structureMat_cold_finger_log, move_structureMat_cold_finger, rotate_structureMat_cold_finger);

    // Cooling Side Block
    G4Box* cooling_side_block = this->coolingSideBlock();
    cooling_side_block_log = new G4LogicalVolume(cooling_side_block, structureMaterial, "cooling_side_block_log", 0, 0, 0);
    cooling_side_block_log->SetVisAttributes(structureMat_vis_att);

    G4Box* cooling_bar = this->coolingBar();
    cooling_bar_log = new G4LogicalVolume(cooling_bar, structureMaterial, "cooling_bar_log", 0, 0, 0);
    cooling_bar_log->SetVisAttributes(structureMat_vis_att);

    G4RotationMatrix* rotate_piece[8] ;

    x = this->cooling_side_block_thickness/2.0
            + this->cold_finger_end_plate_thickness
            + this->can_face_thickness/2.0
            + this->germanium_dist_from_can_face +this->germanium_length
            + this->cold_finger_space +this->shift +this->applied_back_shift ;

    y01 = this->germanium_width - this->cooling_side_block_horizontal_depth/2.0 ;

    y02 = this->germanium_width
            - this->cooling_side_block_horizontal_depth
            - (this->germanium_width
               - this->cooling_side_block_horizontal_depth
               - this->structureMat_cold_finger_radius)/2.0 ;

    for( i = 0 ; i < 4 ; i++ )
    {
        // Add cooling side blocks and cooling bars at the same time.
        rotate_piece[2*i] = new G4RotationMatrix ; // cooling side blocks
        rotate_piece[2*i]->rotateZ( M_PI/2.0 ) ;
        rotate_piece[2*i]->rotateX( M_PI/2.0 - i*M_PI/2.0 ) ;

        y = -y01 * cos( i*M_PI/2.0 ) ;
        z = -y01 * sin( i*M_PI/2.0 ) ;

        move_piece[2*i] = G4ThreeVector( x, y, z ) ;

        rotate_piece[2*i+1] = new G4RotationMatrix ; // cooling bars
        rotate_piece[2*i+1]->rotateZ( M_PI/2.0 ) ;
        rotate_piece[2*i+1]->rotateX( M_PI/2.0 - i*M_PI/2.0 ) ;

        y = -y02 * cos( i*M_PI/2.0 ) ;
        z = -y02 * sin( i*M_PI/2.0 ) ;

        move_piece[2*i+1] = G4ThreeVector( x, y, z ) ;

        this->assembly->AddPlacedVolume( cooling_side_block_log, move_piece[2*i], rotate_piece[2*i] );
        this->assembly->AddPlacedVolume( cooling_bar_log, move_piece[2*i+1], rotate_piece[2*i+1] );
    }

    // Triangle posts from "Cut A"
    // First, find how far from the centre to place the tips of the triangles

    G4double distance_of_the_tip = this->germanium_separation/2.0
            + (this->germanium_width/2.0 - this->germanium_shift) //centre of the middle hole
            + sqrt( pow((this->germanium_outer_radius
                         + this->triangle_posts_distance_from_crystals), 2.0)
                    - pow(this->germanium_width/2.0 - this->germanium_shift, 2.0) );

    // The distance of the base from the detector centre
    G4double distance_of_the_base = this->germanium_separation/2.0 + this->germanium_width;

    G4double triangle_post_length = this->germanium_length - this->triangle_post_starting_depth
            + this->cold_finger_space;

    G4Trd* triangle_post = this->trianglePost();
    triangle_post_log = new G4LogicalVolume(triangle_post, structureMaterial, "triangle_post_log", 0, 0, 0);

    x = this->can_face_thickness/2.0
            + this->germanium_dist_from_can_face +this->germanium_length
            + this->cold_finger_space +this->shift +this->applied_back_shift
            - triangle_post_length/2.0 ;

    y0 = (distance_of_the_base + distance_of_the_tip)/2.0 ;

    for( i = 0 ; i < 4 ; i++ )
    {
        rotate_piece[i] = new G4RotationMatrix ;
        rotate_piece[i]->rotateY(0) ; // Initially included...seems a little redundant.
        rotate_piece[i]->rotateX( M_PI/2.0 - i*M_PI/2.0 ) ;  // 4, 1, 2, 3

        y = -y0 * cos( i*M_PI/2.0 ) ;
        z =  y0 * sin( i*M_PI/2.0 ) ;

        move_piece[i] = G4ThreeVector( x, y, z ) ;

        this->assembly->AddPlacedVolume(triangle_post_log, move_piece[i], rotate_piece[i]);
    }
}//end ::ConstructColdFinger

///////////////////////////////////////////////////////////////////////
// methods used in ConstructColdFinger()
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemGriffin::endPlate()
{
    G4double half_thickness_x = this->cold_finger_end_plate_thickness/2.0;
    G4double half_length_y = this->germanium_width;  //Since it must cover the back of the detector
    G4double half_length_z = half_length_y;  //Since it is symmetric

    G4Box* end_plate = new G4Box("end_plate", half_thickness_x, half_length_y, half_length_z);

    return end_plate;
}//end ::endPlate


G4Tubs* DetectionSystemGriffin::finger()
{
    G4double inner_radius = 0.0*cm;
    G4double outer_radius = this->cold_finger_radius;
    G4double half_length_z = (this->cold_finger_length
                              - this->structureMat_cold_finger_thickness)/2.0;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* finger = new G4Tubs("finger", inner_radius, outer_radius, half_length_z, start_angle, final_angle);

    return finger;
}//end ::finger


G4SubtractionSolid* DetectionSystemGriffin::extraColdBlock()
{
    G4double half_width_x = this->germanium_width;
    G4double half_width_y = half_width_x;
    G4double half_thickness_z = this->extra_block_thickness/2.0;

    G4Box* plate = new G4Box("plate", half_width_x, half_width_y, half_thickness_z);

    G4double inner_radius = 0.0;
    G4double outer_radius = this->extra_block_inner_diameter/2.0;

    G4double half_height_z = half_thickness_z +1.0*cm;  //+1.0*cm just to make sure the hole goes completely through the plate
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* hole = new G4Tubs("hole", inner_radius, outer_radius, half_height_z, start_angle, final_angle);

    G4SubtractionSolid* extra_cold_block = new G4SubtractionSolid("extra_cold_block", plate, hole);



    // Fix overlap
    G4SubtractionSolid* fix_extra_cold_block[4];
    G4double cut_tube_x0, cut_tube_y0, cut_yMov, cut_zMov;
    G4RotationMatrix* cut_rotate_piece[4] ;
    G4ThreeVector cut_move_piece[4] ;
    G4Tubs* cut_tube_location[4] ;

    cut_tube_x0 = (this->detector_total_length +this->bent_end_length
          - this->rear_plate_thickness -this->can_face_thickness)/2.0
            + this->shift +this->applied_back_shift ;

    cut_tube_y0 = (this->detector_total_width/2.0) - (this->bent_end_length * tan(this->bent_end_angle)) ;

    for( G4int i = 0 ; i < 4 ; i++ )
    {
        // Lower Left, Upper Left, Upper Right, Lower Right
        cut_tube_location[i] = this->cornerTube() ;
        cut_rotate_piece[i] = new G4RotationMatrix ;
        cut_rotate_piece[i]->rotateZ( -M_PI/2.0 + (i-1)*M_PI/2.0 ) ;

        cut_yMov = cut_tube_y0 * ( -cos(i*M_PI/2.0) + sin(i*M_PI/2.0) ) ;
        cut_zMov = cut_tube_y0  * ( -cos(i*M_PI/2.0) - sin(i*M_PI/2.0) ) ;

        cut_move_piece[i] = G4ThreeVector(cut_zMov, cut_yMov, 0) ;

        if(i == 0) {
            fix_extra_cold_block[i] = new G4SubtractionSolid("extra_cold_block", extra_cold_block, cut_tube_location[i], cut_rotate_piece[i], cut_move_piece[i] ) ;
        }
        else {
            fix_extra_cold_block[i] = new G4SubtractionSolid("extra_cold_block", fix_extra_cold_block[i-1], cut_tube_location[i], cut_rotate_piece[i], cut_move_piece[i] ) ;
        }
    }






    return fix_extra_cold_block[3];
}//end ::extraColdBlock


G4Trd* DetectionSystemGriffin::trianglePost()
{
    // Calculations that are also done in ConstructColdFinger for positioning
    // First, find how far from the centre to place the tips of the triangles
    G4double distance_of_the_tip	= this->germanium_separation/2.0
            + (this->germanium_width/2.0 - this->germanium_shift) //centre of the middle hole
            + sqrt( pow((this->germanium_outer_radius
                         + this->triangle_posts_distance_from_crystals), 2.0)
                    - pow(this->germanium_width/2.0 - this->germanium_shift, 2.0) );

    // The distance of the base from the detector centre
    G4double distance_of_the_base = this->germanium_separation/2.0
            + this->germanium_width;

    // The distance away from the boundary between crystals of the side points
    G4double distance_of_the_side_points 	= sqrt( pow((this->germanium_outer_radius
                                                         + this->triangle_posts_distance_from_crystals), 2.0)
                                                    - pow(this->germanium_width/2.0 +/*notice*/ this->germanium_shift, 2.0) );

    // Measurements to make the posts with
    G4double length = this->germanium_length - this->triangle_post_starting_depth
            + this->cold_finger_space;

    G4double base_to_tip_height = (distance_of_the_base-distance_of_the_tip);

    G4double half_width_of_base = distance_of_the_side_points;

    G4double half_width_of_top = this->trianglePostDim; //the easiest way to make a triangle
    //G4Trd( const G4String& pName,G4double  dx1, G4double dx2, G4double  dy1, G4double dy2,G4double  dz )
    G4Trd* triangle_post = new G4Trd(	"triangle_post", length / 2.0, length / 2.0,
                                        half_width_of_top, half_width_of_base, base_to_tip_height / 2.0 );

    return triangle_post;
}//end ::trianglePost


G4Tubs* DetectionSystemGriffin::airHole()
{
    G4double inner_radius = 0.0*cm;
    G4double outer_radius = this->fet_air_hole_radius;
    G4double half_length_z = this->cold_finger_end_plate_thickness/2.0;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* fet_air_hole = new G4Tubs("fet_air_hole", inner_radius,
                                      outer_radius, half_length_z, start_angle, final_angle);

    return fet_air_hole;
}//end ::airHole

G4Tubs* DetectionSystemGriffin::airHoleCut()
{
    G4double inner_radius = 0.0*cm;
    G4double outer_radius = this->fet_air_hole_radius;
    G4double half_length_z = this->cold_finger_end_plate_thickness/2.0 + 1.0*cm;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* fet_air_hole = new G4Tubs("fet_air_hole", inner_radius,
                                      outer_radius, half_length_z, start_angle, final_angle);

    return fet_air_hole;
}//end ::airHoleCut

G4Box* DetectionSystemGriffin::coolingBar()
{
    G4double half_thickness_x = this->cooling_bar_thickness/2.0;
    G4double half_length_y = this->cooling_bar_width/2.0;
    G4double half_length_z = (this->germanium_width
                              - this->cooling_side_block_horizontal_depth
                              - this->structureMat_cold_finger_radius)/2.0;

    G4Box* cooling_bar = new G4Box("cooling_bar", half_thickness_x, half_length_y, half_length_z);

    return cooling_bar;
}//end ::coolingBar


G4Box* DetectionSystemGriffin::coolingSideBlock()
{
    G4double half_thickness_x = this->cooling_side_block_width/2.0;
    G4double half_length_y = this->cooling_side_block_thickness/2.0;
    G4double half_length_z = this->cooling_side_block_horizontal_depth/2.0;

    G4Box* cooling_side_block = new G4Box("cooling_side_block",
                                          half_thickness_x, half_length_y, half_length_z);

    return cooling_side_block;
}//end ::coolingSideBlock


G4Tubs* DetectionSystemGriffin::structureMatColdFinger()
{
    G4double inner_radius = 0.0*cm;
    G4double outer_radius = this->structureMat_cold_finger_radius;
    G4double half_length_z = this->structureMat_cold_finger_thickness/2.0;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* structureMat_cold_finger= new G4Tubs("structureMat_cold_finger",
                                                 inner_radius, outer_radius, half_length_z, start_angle, final_angle);

    return structureMat_cold_finger;
}//end ::structureMatColdFinger


///////////////////////////////////////////////////////////////////////
// ConstructNewSuppressorCasingWithShells builds the suppressor that 
// surrounds the can. It tapers off at the front, and there is a thick 
// piece covering the back of the can, which is divided in the middle,
// as per the design in the specifications. Also, there is a layer
// of structureMat that surrounds the physical pieces.
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructNewSuppressorCasingWithShells()
{
    G4int i;
    G4double x0, y0, z0, x, y, z;

    G4Material* structureMaterial = G4Material::GetMaterial(this->structureMaterial);
    if( !structureMaterial ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }
    G4Material* materialBGO = G4Material::GetMaterial(this->BGO_material);
    if( !materialBGO ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }

    // Change some values to accomodate the shells
    // Replacement for this->side_suppressor_length
    G4double shell_side_suppressor_length = this->side_suppressor_length
            + (this->suppressor_shell_thickness*2.0);

    // Replacement for this->suppressor_extension_length
    G4double shell_suppressor_extension_length = this->suppressor_extension_length
            + (this->suppressor_shell_thickness*2.0)*(1.0/tan(this->bent_end_angle)
                                                      - tan(this->bent_end_angle));

    // Replacement for this->suppressor_extension_angle: must totally recalculate
    G4double shell_suppressor_extension_angle = atan(((this->suppressor_back_radius
                                                       + this->bent_end_length +(this->BGO_can_seperation
                                                                                 + this->side_BGO_thickness + this->suppressor_shell_thickness*2.0)
                                                       / tan(this->bent_end_angle)
                                                       - (this->suppressor_extension_thickness + this->suppressor_shell_thickness*2.0)
                                                       * sin(this->bent_end_angle))
                                                      * tan(this->bent_end_angle) -(this->suppressor_forward_radius +this->hevimet_tip_thickness)
                                                      * sin(this->bent_end_angle))/(shell_suppressor_extension_length));

    G4VisAttributes* Suppressor_vis_att = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
    Suppressor_vis_att->SetVisibility(true);
    G4VisAttributes* innards_vis_att = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
    innards_vis_att->SetVisibility(true);
    G4VisAttributes* back_innards_vis_att = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    back_innards_vis_att->SetVisibility(true);

    // first we add the four pieces of back suppressor, and their shells

    G4SubtractionSolid* back_quarter_suppressor = this->backSuppressorQuarter();

    G4SubtractionSolid* back_quarter_suppressor_shell = this->shellForBackSuppressorQuarter();

    // Insert the structureMat shell first.  The shells must be given numbers, rather than
    // the suppressor pieces themselves.  As an error checking method, the suppressor
    // pieces are given a copy number value far out of range of any useful copy number.
    if(include_back_suppressors)
    {
        back_quarter_suppressor_shell_log = new G4LogicalVolume(
                    back_quarter_suppressor_shell, structureMaterial,
                    "back_quarter_suppressor_log", 0, 0, 0);

        back_quarter_suppressor_shell_log->SetVisAttributes( Suppressor_vis_att );

        G4Material* back_material = G4Material::GetMaterial(this->back_suppressor_material);
        if( !back_material )
        {
            G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
            exit(1);
        }

        back_quarter_suppressor_log = new G4LogicalVolume(back_quarter_suppressor, back_material,
                                                          "back_quarter_suppressor_log", 0, 0, 0);
        back_quarter_suppressor_log->SetVisAttributes(back_innards_vis_att);

        this->suppressorBackAssembly->AddPlacedVolume(back_quarter_suppressor_log, this->move_null, this->rotate_null);

        x0 = 	( this->back_BGO_thickness -this->can_face_thickness)/2.0 + this->suppressor_shell_thickness
                + this->detector_total_length +this->BGO_can_seperation
                + this->shift + this->applied_back_shift ;

        y0 = this->detector_total_width/4.0 ;

        z0 = this->detector_total_width/4.0 ;

        G4RotationMatrix* rotate_back_suppressor_shells[4] ;
        G4ThreeVector move_back_quarter_suppressor[4] ;

        for( G4int i = 0 ; i < 4 ; i++ )
        {
            rotate_back_suppressor_shells[i] = new G4RotationMatrix ;
            rotate_back_suppressor_shells[i]->rotateX( -M_PI / 2.0 + i * M_PI / 2.0 ) ;

            x = x0 ;

            y = y0 * ( sin( i * M_PI/2.0 ) - cos( i * M_PI/2.0 ) ) ;
            z = -z0 * ( sin( i * M_PI/2.0 ) + cos( i * M_PI/2.0 ) ) ;

            move_back_quarter_suppressor[i] = G4ThreeVector( x, y, z ) ;

            this->backAndSideSuppressorShellAssembly->AddPlacedVolume( back_quarter_suppressor_shell_log, move_back_quarter_suppressor[i], rotate_back_suppressor_shells[i] ) ;
        }
    }

    // now we add the side pieces of suppressor that taper off towards the front of the can

    // Define the structureMat shell logical volume
    G4cout << "Calling: shellForFrontRightSlantSuppressor" << G4endl ;
    G4SubtractionSolid* right_suppressor_shell = this->shellForFrontSlantSuppressor("right");


    right_suppressor_shell_log = new G4LogicalVolume(right_suppressor_shell, structureMaterial,
                                                     "right_suppressor_shell_log", 0,0,0);
    right_suppressor_shell_log->SetVisAttributes(Suppressor_vis_att);

    G4cout << "Calling: shellForFrontLeftSlantSuppressor" << G4endl ;
    G4SubtractionSolid* left_suppressor_shell = this->shellForFrontSlantSuppressor("left");

    left_suppressor_shell_log = new G4LogicalVolume(left_suppressor_shell, structureMaterial,
                                                    "left_suppressor_shell_log", 0,0,0);
    left_suppressor_shell_log->SetVisAttributes(Suppressor_vis_att);

    G4SubtractionSolid* right_suppressor = this->frontSlantSuppressor("right", false); // Right, non-chopping.

    right_suppressor_log = new G4LogicalVolume(right_suppressor, materialBGO, "right_suppressor_casing_log", 0, 0, 0);
    right_suppressor_log->SetVisAttributes(innards_vis_att);

    G4SubtractionSolid* left_suppressor = this->frontSlantSuppressor("left", false); // Left, non-chopping.

    left_suppressor_log = new G4LogicalVolume(left_suppressor, materialBGO, "left_suppressor_casing_log", 0, 0, 0);
    left_suppressor_log->SetVisAttributes(innards_vis_att);

    this->rightSuppressorCasingAssembly->AddPlacedVolume(right_suppressor_log, this->move_null, this->rotate_null);
    this->leftSuppressorCasingAssembly->AddPlacedVolume(left_suppressor_log, this->move_null, this->rotate_null);


    /////////////////////////////////////////////////////////////////////
    // Note : Left and Right are read from the BACK of the detector
    // Suppressors 1 and 2 cover germanium 1
    // Suppressors 3 and 4 cover germanium 2
    // Suppressors 5 and 6 cover germanium 3
    // Suppressors 7 and 8 cover germanium 4
    /////////////////////////////////////////////////////////////////////

    x0 =    shell_side_suppressor_length/2.0 -this->can_face_thickness/2.0
            + this->bent_end_length + (this->BGO_can_seperation
                                       + this->BGO_chopped_tip) / tan(this->bent_end_angle)
            + this->shift + this->applied_back_shift ;

    y0 =    this->side_BGO_thickness/2.0 + this->suppressor_shell_thickness
            + this->detector_total_width/2.0 + this->BGO_can_seperation ;

    z0 =    this->side_BGO_thickness/2.0 + this->suppressor_shell_thickness
            + this->detector_total_width/2.0 + this->BGO_can_seperation ;

    G4RotationMatrix* rotate_suppressor_extension_shell[8] ;
    G4ThreeVector move_suppressor_extension_shell[8] ;

    for( i = 0 ; i < 4 ; i++ )
    {
        //********************** RIGHT **********************//
        rotate_suppressor_extension_shell[2*i] = new G4RotationMatrix ;
        rotate_suppressor_extension_shell[2*i]->rotateZ( M_PI/2.0 ) ;
        rotate_suppressor_extension_shell[2*i]->rotateY( M_PI/2.0 ) ;
        rotate_suppressor_extension_shell[2*i]->rotateX( M_PI * ( 1 - i/2.0 ) ) ;

        x = x0 ;

        y = y0 * ( -cos( i * M_PI/2.0 ) + sin( i * M_PI/2.0 ) ) / ( 1 + ( i + 1 ) % 2 ) ; // -y/2, y, y/2, -y

        z = z0 * ( cos( i * M_PI/2.0 ) + sin( i * M_PI/2.0 ) ) / ( 1 + i % 2 ) ;  // z, z/2, -z, -z/2

        move_suppressor_extension_shell[2*i] = G4ThreeVector( x, y, z ) ;

        this->backAndSideSuppressorShellAssembly->AddPlacedVolume(right_suppressor_shell_log, move_suppressor_extension_shell[2*i], rotate_suppressor_extension_shell[2*i] ) ;

        //*********************** LEFT ***********************//
        rotate_suppressor_extension_shell[2*i+1] = new G4RotationMatrix ;
        rotate_suppressor_extension_shell[2*i+1]->rotateY( -M_PI/2.0 ) ;
        rotate_suppressor_extension_shell[2*i+1]->rotateX( M_PI * ( 1 - i/2.0 ) ) ;

        x = x0 ;

        y = y0 * ( cos( i * M_PI/2.0 ) - sin( i * M_PI/2.0 ) ) / ( 1 + i % 2 ) ; // y, -y/2, -y, y/2

        z = -z0 * ( cos( i * M_PI/2.0 ) + sin( i * M_PI/2.0 ) ) / ( 1 + ( i + 1 ) % 2 ) ; // -z/2, -z, z/2, z

        move_suppressor_extension_shell[2*i+1] = G4ThreeVector( x, y, z ) ;

        this->backAndSideSuppressorShellAssembly->AddPlacedVolume(left_suppressor_shell_log, move_suppressor_extension_shell[2*i+1], rotate_suppressor_extension_shell[2*i+1] ) ;

    }


    // now we add the side pieces of suppressor that extend out in front of the can when it's in the back position

    // Define the shell right logical volume
    G4cout << "Calling: shellForRightSuppressorExtension" << G4endl ;
    // G4SubtractionSolid* right_suppressor_shell_extension = this->shellForRightSuppressorExtension();
    G4SubtractionSolid* right_suppressor_shell_extension = this->shellForSuppressorExtension("right");

    right_suppressor_shell_extension_log = new G4LogicalVolume(right_suppressor_shell_extension,
                                                               materialBGO, "right_suppressor_shell_extension_log", 0, 0, 0);
    right_suppressor_shell_extension_log->SetVisAttributes(Suppressor_vis_att);

    G4SubtractionSolid* right_suppressor_extension = this->sideSuppressorExtension( "right", false ); // Right, non-chopping // CALLED

    right_suppressor_extension_log = new G4LogicalVolume(right_suppressor_extension, materialBGO,
                                                         "right_suppressor_extension_log", 0, 0, 0);
    right_suppressor_extension_log->SetVisAttributes(innards_vis_att);

    // Define the left shell logical volume
    G4cout << "Calling: shellForLeftSuppressorExtension" << G4endl ;
    // G4SubtractionSolid* left_suppressor_shell_extension = this->shellForLeftSuppressorExtension();
    G4SubtractionSolid* left_suppressor_shell_extension = this->shellForSuppressorExtension("left");


    left_suppressor_shell_extension_log = new G4LogicalVolume(left_suppressor_shell_extension,
                                                              materialBGO, "left_suppressor_shell_extension_log", 0, 0, 0);
    left_suppressor_shell_extension_log->SetVisAttributes(Suppressor_vis_att);

    G4SubtractionSolid* left_suppressor_extension = this->sideSuppressorExtension( "left", false ); // Left, Non-chopping // CALLED

    left_suppressor_extension_log = new G4LogicalVolume(left_suppressor_extension, materialBGO,
                                                        "left_suppressor_extension_log", 0, 0, 0);
    left_suppressor_extension_log->SetVisAttributes(innards_vis_att);

    this->rightSuppressorExtensionAssembly->AddPlacedVolume(right_suppressor_extension_log, this->move_null, this->rotate_null);
    this->leftSuppressorExtensionAssembly->AddPlacedVolume(left_suppressor_extension_log, this->move_null, this->rotate_null);

    //geometry objects for the following:
    G4RotationMatrix* rotateExtension[8];
    G4ThreeVector moveExtension[8];

    x0 =  - this->can_face_thickness/2.0 -(shell_suppressor_extension_length/2.0
                                           - (this->suppressor_extension_thickness + this->suppressor_shell_thickness*2.0)
                                           * tan(this->bent_end_angle)/2.0)
            * cos(this->bent_end_angle) +this->bent_end_length +(this->BGO_can_seperation
                                                                 + this->side_BGO_thickness
                                                                 + this->suppressor_shell_thickness*2.0)/tan(this->bent_end_angle)
            - (this->suppressor_extension_thickness + this->suppressor_shell_thickness*2.0)
            * sin(this->bent_end_angle) +this->suppShift +this->suppressor_back_radius
            - this->suppressor_forward_radius ;

    y0 =  - ( shell_suppressor_extension_length * tan( shell_suppressor_extension_angle ) / 2.0
              + ( this->suppressor_forward_radius + this->hevimet_tip_thickness )
              * sin(this->bent_end_angle))/2.0;

    z0 =  - ((this->suppressor_extension_thickness/2.0 + this->suppressor_shell_thickness)
             / cos(this->bent_end_angle)
             + (shell_suppressor_extension_length/2.0 -(this->suppressor_extension_thickness
                                                        + this->suppressor_shell_thickness*2.0)
                * tan(this->bent_end_angle)/2.0)*sin(this->bent_end_angle) -(this->suppressor_back_radius
                                                                             + this->bent_end_length +(this->BGO_can_seperation
                                                                                                       + this->side_BGO_thickness + this->suppressor_shell_thickness*2.0)
                                                                             / tan(this->bent_end_angle) -(this->suppressor_extension_thickness
                                                                                                           + this->suppressor_shell_thickness*2.0)
                                                                             * sin(this->bent_end_angle)) * tan(this->bent_end_angle)) ;

    if( this->suppressor_position_selector == 0 && this->include_extension_suppressors)
    {

        // these two parameters are for shifting the extensions back and out when in their BACK position
        G4double extension_back_shift =   this->air_box_front_length
                - (this->hevimet_tip_thickness +shell_suppressor_extension_length
                   + (this->suppressor_extension_thickness + this->suppressor_shell_thickness*2.0)
                   * tan(this->bent_end_angle)) * cos(this->bent_end_angle) ;

        G4double extension_radial_shift = extension_back_shift * tan(this->bent_end_angle) ;

        // The suppressors are put into the back position

        x0 += extension_back_shift ;

        z0 += extension_radial_shift ;

        for(i=0; i<4; i++)
        {
            rotateExtension[i*2] = new G4RotationMatrix;
            rotateExtension[i*2]->rotateZ(M_PI/2.0);
            rotateExtension[i*2]->rotateY(this->bent_end_angle);
            rotateExtension[i*2]->rotateX(M_PI - M_PI/2.0*i);

            x = x0;
            y = y0*cos(i*M_PI/2) + z0*sin(i*M_PI/2);
            z = z0*cos(i*M_PI/2) - y0*sin(i*M_PI/2);

            moveExtension[i*2] = G4ThreeVector(x, y, z);

            this->extensionSuppressorShellAssembly->AddPlacedVolume(right_suppressor_shell_extension_log, moveExtension[i*2], rotateExtension[i*2]);

            rotateExtension[i*2+1] = new G4RotationMatrix;
            rotateExtension[i*2+1]->rotateY(M_PI/2.0);
            rotateExtension[i*2+1]->rotateZ(M_PI/2.0 + this->bent_end_angle);
            rotateExtension[i*2+1]->rotateX(M_PI - M_PI/2.0*i);

            x = x0;
            y = -z0*cos(i*M_PI/2) - y0*sin(i*M_PI/2);
            z = -y0*cos(i*M_PI/2) + z0*sin(i*M_PI/2);

            moveExtension[i*2+1] = G4ThreeVector(x, y, z);

            this->extensionSuppressorShellAssembly->AddPlacedVolume(left_suppressor_shell_extension_log, moveExtension[i*2+1], rotateExtension[i*2+1]);
        }


    }//end if(detectors forward) statement

    // Otherwise, put them forward
    else if( this->suppressor_position_selector == 1 && this->include_extension_suppressors)
    {

        for( i = 0 ; i < 4 ; i++){
            rotateExtension[2*i] = new G4RotationMatrix;
            rotateExtension[2*i]->rotateZ(M_PI/2.0);
            rotateExtension[2*i]->rotateY(this->bent_end_angle);
            rotateExtension[2*i]->rotateX(M_PI - i*M_PI/2.0);

            x = x0;
            y = y0*cos(i*M_PI/2) + z0*sin(i*M_PI/2);
            z = z0*cos(i*M_PI/2) - y0*sin(i*M_PI/2);

            moveExtension[i*2] = G4ThreeVector(x, y, z);

            this->extensionSuppressorShellAssembly->AddPlacedVolume(right_suppressor_shell_extension_log, moveExtension[2*i], rotateExtension[2*i]);

            rotateExtension[2*i+1] = new G4RotationMatrix;
            rotateExtension[2*i+1]->rotateY(M_PI/2.0);
            rotateExtension[2*i+1]->rotateZ(M_PI/2.0 + this->bent_end_angle);
            rotateExtension[2*i+1]->rotateX(M_PI - i*M_PI/2);

            x =  x0;
            y =  -z0 * cos(i*M_PI/2) - y0 * sin(i*M_PI/2);
            z =  -y0 * cos(i*M_PI/2) + z0 * sin(i*M_PI/2);

            moveExtension[i*2+1] = G4ThreeVector(x, y, z);

            this->extensionSuppressorShellAssembly->AddPlacedVolume(left_suppressor_shell_extension_log, moveExtension[2*i+1], rotateExtension[2*i+1]);
        }

    }//end if(detectors back) statement

}//end ::ConstructNewSuppressorCasingWithShells


void DetectionSystemGriffin::ConstructNewSuppressorCasingDetectorSpecificDeadLayer(G4int det, G4int cry)
{
    G4String strdet = G4UIcommand::ConvertToString(det);
    G4String strcry = G4UIcommand::ConvertToString(cry);

    G4Material* structureMaterial = G4Material::GetMaterial(this->structureMaterial);
    if( !structureMaterial ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }
    G4Material* materialBGO = G4Material::GetMaterial(this->BGO_material);
    if( !materialBGO ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }

    G4VisAttributes* Suppressor_vis_att = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
    Suppressor_vis_att->SetVisibility(true);
    G4VisAttributes* innards_vis_att = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
    innards_vis_att->SetVisibility(true);
    G4VisAttributes* back_innards_vis_att = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    back_innards_vis_att->SetVisibility(true);

    // first we add the four pieces of back suppressor, and their shells

    G4SubtractionSolid* back_quarter_suppressor = this->backSuppressorQuarter();

    // Insert the structureMat shell first.  The shells must be given numbers, rather than
    // the suppressor pieces themselves.  As an error checking method, the suppressor
    // pieces are given a copy number value far out of range of any useful copy number.
    if(include_back_suppressors)
    {
        G4Material* back_material = G4Material::GetMaterial(this->back_suppressor_material);

        if( !back_material )
        {
            G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
            exit(1);
        }

        G4String back_quarter_suppressor_name = "back_dls_quarter_suppressor_" + strdet + "_" + strcry + "_log" ;

        back_quarter_suppressor_log = new G4LogicalVolume(back_quarter_suppressor, back_material,
                                                          back_quarter_suppressor_name, 0, 0, 0);

        back_quarter_suppressor_log->SetVisAttributes(back_innards_vis_att);

        this->suppressorBackAssemblyCry[cry]->AddPlacedVolume(back_quarter_suppressor_log, this->move_null, this->rotate_null);

    }

    // now we add the side pieces of suppressor that taper off towards the front of the can

    G4String right_suppressor_casing_name = "right_dls_suppressor_casing_" + strdet + "_" + strcry + "_log" ;
    G4String left_suppressor_casing_name = "left_dls_suppressor_casing_" + strdet + "_" + strcry + "_log" ;

    G4SubtractionSolid* right_suppressor = this->frontSlantSuppressor("right", false); // Right, non-chopping. // CALLED

    right_suppressor_log = new G4LogicalVolume(right_suppressor, materialBGO, right_suppressor_casing_name, 0, 0, 0);
    right_suppressor_log->SetVisAttributes(innards_vis_att);

    G4SubtractionSolid* left_suppressor = this->frontSlantSuppressor("left", false); // Left, non-chopping. // CALLED

    left_suppressor_log = new G4LogicalVolume(left_suppressor, materialBGO, left_suppressor_casing_name, 0, 0, 0);
    left_suppressor_log->SetVisAttributes(innards_vis_att);

    this->rightSuppressorCasingAssemblyCry[cry]->AddPlacedVolume(right_suppressor_log, this->move_null, this->rotate_null);
    this->leftSuppressorCasingAssemblyCry[cry]->AddPlacedVolume(left_suppressor_log, this->move_null, this->rotate_null);

    G4String right_suppressor_extension_name = "right_dls_suppressor_extension_" + strdet + "_" + strcry + "_log" ;
    G4String left_suppressor_extension_name = "left_dls_suppressor_extension_" + strdet + "_" + strcry + "_log" ;

    // now we add the side pieces of suppressor that extend out in front of the can when it's in the back position

    G4SubtractionSolid* right_suppressor_extension = this->sideSuppressorExtension( "right", false ); // Right, non-chopping // CALLED

    right_suppressor_extension_log = new G4LogicalVolume(	right_suppressor_extension, materialBGO,
                                                            right_suppressor_extension_name, 0, 0, 0);
    right_suppressor_extension_log->SetVisAttributes(innards_vis_att);

    G4SubtractionSolid* left_suppressor_extension = this->sideSuppressorExtension( "left", false ); // Left, Non-chopping // CALLED

    left_suppressor_extension_log = new G4LogicalVolume(left_suppressor_extension, materialBGO,
                                                        left_suppressor_extension_name, 0, 0, 0);
    left_suppressor_extension_log->SetVisAttributes(innards_vis_att);

    this->rightSuppressorExtensionAssemblyCry[cry]->AddPlacedVolume(right_suppressor_extension_log, this->move_null, this->rotate_null);
    this->leftSuppressorExtensionAssemblyCry[cry]->AddPlacedVolume(left_suppressor_extension_log, this->move_null, this->rotate_null);


}//end ::ConstructNewSuppressorCasingDetectorSpecificDeadLayer


///////////////////////////////////////////////////////////////////////
// This is a messy place to put it, but these are the methods to make
// the electrodeMat layers between the crystals
///////////////////////////////////////////////////////////////////////
G4UnionSolid* DetectionSystemGriffin::interCrystalelectrodeMatBack()
{
    G4double distance_of_the_triangle_tips = this->germanium_separation/2.0
            + (this->germanium_width/2.0 - this->germanium_shift) //centre of the middle hole
            + sqrt( pow((this->germanium_outer_radius
                         + this->triangle_posts_distance_from_crystals), 2.0)
                    - pow(this->germanium_width/2.0 - this->germanium_shift, 2.0) );

    G4double extent_of_the_electrodeMat_pieces = distance_of_the_triangle_tips
            - this->triangle_posts_distance_from_crystals;

    G4Box* electrodeMat_piece1 = new G4Box( "electrodeMat_piece1", extent_of_the_electrodeMat_pieces,
                                            (this->germanium_length-this->germanium_bent_length)/2.0,
                                            this->inter_crystal_electrodeMat_thickness/2.0);

    G4Box* electrodeMat_piece2 = new G4Box( "electrodeMat_piece2", extent_of_the_electrodeMat_pieces,
                                            (this->germanium_length-this->germanium_bent_length)/2.0,
                                            this->inter_crystal_electrodeMat_thickness/2.0);

    G4RotationMatrix* rotate_piece2 = new G4RotationMatrix;
    rotate_piece2->rotateY(M_PI/2.0);

    G4ThreeVector move_zero(0,0,0);

    G4UnionSolid* inter_crystal_electrodeMat_back = new G4UnionSolid(
                "inter_crystal_electrodeMat_back", electrodeMat_piece1, electrodeMat_piece2,
                rotate_piece2, move_zero);

    return inter_crystal_electrodeMat_back;
} //end ::interCrystalelectrodeMatBack

G4UnionSolid* DetectionSystemGriffin::interCrystalelectrodeMatFront()
{

    G4double distance_of_the_triangle_tips = this->germanium_separation/2.0
            + (this->germanium_width/2.0 - this->germanium_shift) //centre of the middle hole
            + sqrt( pow((this->germanium_outer_radius
                         + this->triangle_posts_distance_from_crystals), 2.0)
                    - pow(this->germanium_width/2.0 - this->germanium_shift, 2.0) );

    G4Trd* electrodeMat_piece1 = new G4Trd("electrodeMat_piece1", distance_of_the_triangle_tips
                                           - this->triangle_posts_distance_from_crystals,
                                           this->germanium_width - this->germanium_bent_length
                                           *tan(this->bent_end_angle),
                                           this->germanium_separation/2.0,
                                           this->germanium_separation/2.0, (this->germanium_bent_length
                                                                            - this->electrodeMat_starting_depth)/2.0);

    G4Trd* electrodeMat_piece2 = new G4Trd("electrodeMat_piece2", distance_of_the_triangle_tips
                                           - this->triangle_posts_distance_from_crystals,
                                           this->germanium_width - this->germanium_bent_length
                                           *tan(this->bent_end_angle),
                                           this->germanium_separation/2.0,
                                           this->germanium_separation/2.0, (this->germanium_bent_length
                                                                            - this->electrodeMat_starting_depth)/2.0);

    G4RotationMatrix* rotate_piece2 = new G4RotationMatrix;
    rotate_piece2->rotateZ(M_PI/2.0);

    G4ThreeVector move_zero(0,0,0);

    G4UnionSolid* inter_crystal_electrodeMat_front = new G4UnionSolid(
                "inter_crystal_electrodeMat_front", electrodeMat_piece1, electrodeMat_piece2,
                rotate_piece2, move_zero);

    return inter_crystal_electrodeMat_front;
} //end ::interCrystalelectrodeMatFront

///////////////////////////////////////////////////////////////////////
// ConstructNewHeavyMet builds the heavy metal that goes on the front 
// of the Suppressor. It only goes on when the extensions are in their 
// forward position   
///////////////////////////////////////////////////////////////////////
void DetectionSystemGriffin::ConstructNewHeavyMet()
{
    G4Material* materialHevimetal = G4Material::GetMaterial("Hevimetal");
    if( !materialHevimetal ) {
        G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
        exit(1);
    }

    G4VisAttributes* HeviMet_vis_att = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    HeviMet_vis_att->SetVisibility(true);

    G4SubtractionSolid* hevimet = this->newHeavyMet();
    
    G4RotationMatrix* rotate_hevimet = new G4RotationMatrix;
    rotate_hevimet->rotateY(-1.0*M_PI/2.0);

    G4ThreeVector move_hevimet( (( this->hevimet_tip_thickness / cos( this->hevimet_tip_angle ) ) * cos( this->bent_end_angle -this->hevimet_tip_angle )
                                 + this->suppressor_forward_radius* tan(this->hevimet_tip_angle)*sin(this->bent_end_angle))/2.0
                                - this->air_box_back_length/2.0 -this->air_box_front_length, 0.0, 0.0);

    // NOTE** The hevimet does not require the "shift" parameter, as the air_box has been designed
    // so that the hevimet sits right at the front, and has been moved accordingly. Also, the
    // "applied_back_shift" parameter is missing, as the hevimet only goes on in the "detector back" position

    hevimet_log = new G4LogicalVolume(hevimet, materialHevimetal, "hevimet_log", 0, 0, 0);

    hevimet_log->SetVisAttributes(HeviMet_vis_att);

    // this->assembly->AddPlacedVolume(hevimet_log, move_hevimet, rotate_hevimet);
    this->hevimetAssembly->AddPlacedVolume(hevimet_log, move_hevimet, rotate_hevimet);


}//end ::ConstructHeavyMet


///////////////////////////////////////////////////////////////////////
// methods used in ConstructDetector()
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemGriffin::squareFrontFace()
{

    G4double fix_overlap = 250*nm;

    G4double half_thickness_x = this->can_face_thickness/2.0;
    
    G4double half_length_y = ((this->detector_total_width/2.0) -(this->bent_end_length
                                                                *tan(this->bent_end_angle)))-fix_overlap;
    
    G4double half_length_z = half_length_y;

    G4Box* front_face = new G4Box("front_face", half_thickness_x, half_length_y, half_length_z);

    return front_face;

}//end ::squareFrontFace


G4Trap* DetectionSystemGriffin::cornerWedge()
{
    G4double height_y = this->can_face_thickness;
    G4double length_z = this->detector_total_width -2.0*(this->bent_end_length *tan(this->bent_end_angle));
    G4double extension_length_x = height_y *tan(this->bent_end_angle);

    G4Trap* wedge = new G4Trap("wedge", length_z, height_y, extension_length_x, this->wedgeDim);

    return wedge;
}//end ::cornerWedge


///////////////////////////////////////////////////////////////////////
// The bent side pieces attach to the front faceface they angle
// outwards, so the face is smaller than the main part of 
// the can four are needed for a square-faced detector
///////////////////////////////////////////////////////////////////////
G4Para* DetectionSystemGriffin::bentSidePiece()
{

    G4double half_length_x = ((this->bent_end_length -this->can_face_thickness)
                              /cos(this->bent_end_angle))/2.0;

    G4double half_length_y = (this->detector_total_width/2.0) -(this->bent_end_length
                                                                *tan(this->bent_end_angle)) -this->can_face_thickness
            *(1.0 -tan(this->bent_end_angle));

    // Last two operations in "half_length_y" are correction
    // factor to keep the bent pieces from overlapping
    G4double half_length_z = this->can_face_thickness *cos(this->bent_end_angle)/2.0;
    G4double alpha_angle = 0.0*M_PI;
    G4double polar_angle = this->bent_end_angle;
    G4double azimuthal_angle = 0.0*M_PI;

    G4Para* bent_side_piece = new G4Para("bent_side_piece", half_length_x,
                                         half_length_y, half_length_z, alpha_angle, polar_angle, azimuthal_angle);

    return bent_side_piece;

}//end ::bentSidePiece


// this is a conic piece that attaches two bent end pieces
G4Cons* DetectionSystemGriffin::roundedEndEdge()
{
    G4double hole_eliminator = this->can_face_thickness *(1.0 -tan(this->bent_end_angle));

    // this is correction factor to avoid holes caused by bent pieces
    G4double half_length_z = (this->bent_end_length -this->can_face_thickness)/2.0;
    G4double inside_radius_at_apex = 0.0;
    G4double outside_radius_at_apex = this->can_face_thickness *tan(this->bent_end_angle) +hole_eliminator;
    G4double outside_radius_at_base = this->bent_end_length *tan(this->bent_end_angle) +hole_eliminator;
    G4double inside_radius_at_base = outside_radius_at_base -this->can_face_thickness;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 0.5*M_PI;

    G4Cons* rounded_end_edge = new G4Cons("rounded_end_edge", inside_radius_at_apex,
                                          outside_radius_at_apex, inside_radius_at_base, outside_radius_at_base,
                                          half_length_z, start_angle, final_angle);

    return rounded_end_edge;

}//end ::roundedEndEdge


G4Tubs* DetectionSystemGriffin::cornerTube()
{
    G4double outer_radius 	= this->bent_end_length *tan(this->bent_end_angle);
    G4double inner_radius 	= outer_radius - this->can_side_thickness;
    G4double half_height_z 	= (this->detector_total_length -this->bent_end_length
                               -this->rear_plate_thickness)/2.0;

    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 0.5*M_PI;

    G4Tubs* corner_tube = new G4Tubs("corner_tube", inner_radius, outer_radius, half_height_z, start_angle, final_angle);

    return corner_tube;

}//end ::cornerTube


G4Box* DetectionSystemGriffin::sidePanel()
{
    G4double half_length_x = (this->detector_total_length -this->bent_end_length -this->rear_plate_thickness)/2.0;
    G4double half_thickness_y = (this->can_side_thickness)/2.0;
    G4double half_width_z = (this->detector_total_width)/2.0 -(this->bent_end_length *tan(this->bent_end_angle));

    G4Box* side_panel = new G4Box("side_panel", half_length_x, half_thickness_y, half_width_z);

    return side_panel;

}//end ::sidePanel


G4SubtractionSolid* DetectionSystemGriffin::rearPlate()
{

    G4double half_width_x = this->detector_total_width/2.0;
    G4double half_width_y = half_width_x;
    G4double half_thickness_z = this->rear_plate_thickness/2.0;

    G4Box* plate = new G4Box("plate", half_width_x, half_width_y, half_thickness_z);

    G4double inner_radius = 0.0;
    G4double outer_radius = this->cold_finger_outer_shell_radius;
    G4double half_height_z = half_thickness_z +1.0*cm;  // +1.0*cm just to make sure the hole goes completely through the plate
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* hole = new G4Tubs("hole", inner_radius, outer_radius, half_height_z, start_angle, final_angle);

    G4SubtractionSolid* rear_plate = new G4SubtractionSolid("rear_plate", plate, hole);

    return rear_plate;

}//end ::rearPlate


G4Tubs* DetectionSystemGriffin::coldFingerShell()
{
    G4double outer_radius = this->cold_finger_outer_shell_radius;
    G4double inner_radius = outer_radius - this->cold_finger_shell_thickness;
    G4double half_length_z = this->cold_finger_shell_length/2.0;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* cold_finger_shell = new G4Tubs("cold_finger_shell", inner_radius,
                                           outer_radius, half_length_z, start_angle, final_angle);

    return cold_finger_shell;

}//end ::coldFingerShell


G4Tubs* DetectionSystemGriffin::liquidNitrogenTank()
{
    G4double inner_radius = this->coolant_radius - this->coolant_thickness;
    G4double outer_radius = this->coolant_radius;
    G4double half_length_z = (this->coolant_length - 2.0*this->coolant_thickness)/2.0;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* tank = new G4Tubs("tank", inner_radius, outer_radius, half_length_z, start_angle, final_angle);

    return tank;

}//end ::liquidNitrogenTank

G4Tubs* DetectionSystemGriffin::liquidNitrogenTankLid()
{
    G4double inner_radius = 0.0*mm;
    G4double outer_radius = this->coolant_radius;
    G4double half_length_z = (this->coolant_thickness)/2.0;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* tank = new G4Tubs("tank", inner_radius, outer_radius, half_length_z, start_angle, final_angle);

    return tank;

}//end ::liquidNitrogenTankLid

G4Tubs* DetectionSystemGriffin::liquidNitrogen()
{
    G4double inner_radius = 0.0*mm;
    G4double outer_radius = this->coolant_radius - this->coolant_thickness;
    G4double half_length_z = (this->coolant_length - 2.0*this->coolant_thickness)/2.0;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4Tubs* tank = new G4Tubs("tank", inner_radius, outer_radius, half_length_z, start_angle, final_angle);

    return tank;

}//end ::liquidNitrogen

///////////////////////////////////////////////////////////////////////
// methods used in ConstructBasicDetectorBlock()
// Both the rectangularSegment and the trapezoidalSegment join
// to form a UnionSolid
///////////////////////////////////////////////////////////////////////
G4Box* DetectionSystemGriffin::rectangularSegment()
{
    G4double half_length_x = this->detector_block_height/2.0;
    G4double half_length_y = this->detector_block_length/2.0;
    G4double half_length_z = half_length_y; 	// Since it is symmetric

    G4Box* rectangular_segment = new G4Box("rectangular_segment", half_length_x, half_length_y, half_length_z);

    return rectangular_segment;

}// end ::rectangularSegment


G4Trd* DetectionSystemGriffin::trapezoidalSegment()
{
    G4double half_base_x 		= this->detector_block_length/2.0;
    G4double half_top_x 		= half_base_x - (this->detector_block_trapezoidal_height *tan(this->bent_end_angle));
    G4double half_base_y 		= half_base_x; 	// Since it is symmetric
    G4double half_top_y 		= half_top_x;
    G4double half_height_z 	= this->detector_block_trapezoidal_height/2.0;

    G4Trd* trapezoidal_segment = new G4Trd("trapezoidal_segment",
                                           half_base_x, half_top_x, half_base_y, half_top_y, half_height_z);

    return trapezoidal_segment;

}//end ::trapezoidalSegment

///////////////////////////////////////////////////////////////////////
// methods used in ConstructComplexDetectorBlock()
// Starting with a rectangle, it gets chopped by an off-centered
// cylinder, and then the two edges are chopped diagonaly
///////////////////////////////////////////////////////////////////////
G4SubtractionSolid* DetectionSystemGriffin::quarterDetector()
{

    G4double half_width_x 	= this->germanium_width/2.0;
    G4double half_width_y 	= half_width_x;
    G4double half_length_z 	= this->germanium_length/2.0;

    G4Box* rectangular_germanium 	= new G4Box("rectangular_germanium", half_width_x, half_width_y, half_length_z);

    G4double outer_radius = ((this->germanium_width +this->germanium_shift)*1.5)/2.0;

    // 1.5 is almost root 2, so this makes sure the corners are chopped off

    G4double inner_radius = this->germanium_outer_radius;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4ThreeVector move_chopping_cylinder(this->germanium_shift, -(this->germanium_shift), 0);
    G4Tubs* chopping_cylinder = new G4Tubs("chopping_cylinder", inner_radius,
                                           outer_radius, half_length_z+1.0*cm, start_angle, final_angle);

    G4SubtractionSolid* germanium_rounded_corners = new G4SubtractionSolid("germanium_rounded_corners",
                                                                           rectangular_germanium, chopping_cylinder, 0, move_chopping_cylinder);

    G4double base_inner_radius = this->germanium_outer_radius;
    G4double base_outer_radius = 2.0*base_inner_radius;

    // G4double tip_inner_radius = base_inner_radius - this->germanium_bent_length*tan(this->bent_end_angle);
    G4double tip_inner_radius = base_inner_radius
            - this->germanium_corner_cone_end_length*tan(this->bent_end_angle);
    G4double tip_outer_radius = base_outer_radius;

    // G4double con_half_length_z = this->germanium_bent_length/2.0;
    G4double con_half_length_z = this->germanium_corner_cone_end_length/2.0 + this->quarterDetectorCxn;

    G4double initial_angle = acos((this->germanium_width/2.0 +this->germanium_shift)
                                  /this->germanium_outer_radius);
    G4double total_angle = M_PI/2.0 -2.0*initial_angle;

    G4Cons* rounded_edge = new G4Cons("rounded_edge", base_inner_radius,
                                      base_outer_radius, tip_inner_radius, tip_outer_radius,
                                      con_half_length_z, initial_angle, total_angle);


    G4RotationMatrix* rotate_rounded_edge = new G4RotationMatrix;
    rotate_rounded_edge->rotateZ(-M_PI/2.0);

    G4ThreeVector move_rounded_edge(this->germanium_shift, -this->germanium_shift,
                                    this->germanium_length/2.0 -this->germanium_corner_cone_end_length/2.0
                                    + this->quarterDetectorCxnB);

    G4SubtractionSolid* germanium_rounded_edge = new G4SubtractionSolid("germanium_rounded_edge",
                                                                        germanium_rounded_corners, rounded_edge, rotate_rounded_edge,
                                                                        move_rounded_edge);

    // now we make the diagonal slices
    G4double half_chop_piece_width_x 	= this->germanium_width/2.0;
    G4double half_chop_piece_width_y 	= this->germanium_width/2.0;
    G4double half_chop_piece_length_z 	= (this->germanium_bent_length/cos(this->bent_end_angle))/2.0;
    G4Box* chop_piece = new G4Box("chop_piece", half_chop_piece_width_x,
                                  half_chop_piece_width_y, half_chop_piece_length_z);

    G4RotationMatrix* rotate_chop_piece1 = new G4RotationMatrix;
    rotate_chop_piece1->rotateX(-(this->bent_end_angle));

    G4ThreeVector move_chop_piece1(0, this->germanium_width/2.0 +(this->germanium_bent_length
                                                                  /tan(this->bent_end_angle) +this->germanium_width
                                                                  *cos(this->bent_end_angle))/2.0 -((this->germanium_bent_length
                                                                                                     /cos(this->bent_end_angle))/2.0) /sin(this->bent_end_angle),
                                   this->germanium_length/2.0 -this->germanium_bent_length
                                   +(this->germanium_bent_length +this->germanium_width
                                     *sin(this->bent_end_angle))/2.0);

    G4SubtractionSolid* chopped_germanium1 = new G4SubtractionSolid("chopped_germanium1",
                                                                    germanium_rounded_edge, chop_piece, rotate_chop_piece1, move_chop_piece1);


    G4RotationMatrix* rotate_chop_piece2 = new G4RotationMatrix;
    rotate_chop_piece2->rotateY(-(this->bent_end_angle));

    G4ThreeVector move_chop_piece2(-(this->germanium_width/2.0 +(this->germanium_bent_length
                                                                 /tan(this->bent_end_angle) +this->germanium_width
                                                                 *cos(this->bent_end_angle))/2.0 -((this->germanium_bent_length
                                                                                                    /cos(this->bent_end_angle))/2.0) /sin(this->bent_end_angle)), 0,
                                   this->germanium_length/2.0 -this->germanium_bent_length
                                   +(this->germanium_bent_length +this->germanium_width
                                     *sin(this->bent_end_angle))/2.0);

    G4SubtractionSolid* chopped_germanium2 = new G4SubtractionSolid("chopped_germanium2",
                                                                    chopped_germanium1, chop_piece, rotate_chop_piece2, move_chop_piece2);


    // now we make the hole that will go in the back of each quarter detector, /////////////////////////////////////////////////////////
    start_angle = 0.0*M_PI;
    final_angle = 2.0*M_PI;
    G4double hole_radius = this->germanium_hole_radius;
    G4double hole_half_length_z = (this->germanium_length -this->germanium_hole_dist_from_face)/2.0;

    G4Tubs* hole_tubs = new G4Tubs("hole_tubs", 0.0, hole_radius, hole_half_length_z, start_angle, final_angle);
    G4ThreeVector move_hole(this->germanium_shift, -(this->germanium_shift),-((this->germanium_hole_dist_from_face)));

    G4SubtractionSolid* chopped_germanium3 = new G4SubtractionSolid("chopped_germanium3",
                                                                    chopped_germanium2, hole_tubs, 0, move_hole);

    start_angle = 0.0*M_PI;
    final_angle = 2.0*M_PI;
    G4double dead_layer_radius = this->germanium_hole_radius + this->inner_dead_layer_thickness;
    G4double dead_layer_half_length_z = (this->germanium_length
                                         - this->germanium_hole_dist_from_face
                                         + this->inner_dead_layer_thickness)/2.0;

    // now dead layer
    G4Tubs* dead_layer_tubs = new G4Tubs("dead_layer_tubs", 0.0,
                                         dead_layer_radius, dead_layer_half_length_z, start_angle,final_angle);

    G4ThreeVector move_dead_layer(this->germanium_shift, -(this->germanium_shift),
                                  -((this->germanium_hole_dist_from_face
                                     - this->inner_dead_layer_thickness)/2.0));

    G4SubtractionSolid* chopped_germanium4 = new G4SubtractionSolid("chopped_germanium4",
                                                                    chopped_germanium3, dead_layer_tubs, 0, move_dead_layer);

    return chopped_germanium4;

}//end ::quarterDetector


G4SubtractionSolid* DetectionSystemGriffin::quarterSpecificDeadLayerDetector(G4int det, G4int cry)
{

    G4double half_width_x 	= this->germanium_width/2.0;
    G4double half_width_y 	= half_width_x;
    G4double half_length_z 	= this->germanium_length/2.0;

    G4Box* rectangular_germanium 	= new G4Box("rectangular_germanium", half_width_x, half_width_y, half_length_z);

    G4double outer_radius = ((this->germanium_width +this->germanium_shift)*1.5)/2.0;

    // 1.5 is almost root 2, so this makes sure the corners are chopped off

    G4double inner_radius = this->germanium_outer_radius;
    G4double start_angle = 0.0*M_PI;
    G4double final_angle = 2.0*M_PI;

    G4ThreeVector move_chopping_cylinder(this->germanium_shift, -(this->germanium_shift), 0);
    G4Tubs* chopping_cylinder = new G4Tubs("chopping_cylinder", inner_radius,
                                           outer_radius, half_length_z+1.0*cm, start_angle, final_angle);

    G4SubtractionSolid* germanium_rounded_corners = new G4SubtractionSolid("germanium_rounded_corners",
                                                                           rectangular_germanium, chopping_cylinder, 0, move_chopping_cylinder);

    G4double base_inner_radius = this->germanium_outer_radius;
    G4double base_outer_radius = 2.0*base_inner_radius;

    G4double tip_inner_radius = base_inner_radius - this->germanium_corner_cone_end_length*tan(this->bent_end_angle);
    G4double tip_outer_radius = base_outer_radius;

    G4double con_half_length_z = this->germanium_corner_cone_end_length/2.0 + this->quarterDetectorCxn;

    G4double initial_angle = acos((this->germanium_width/2.0 +this->germanium_shift)
                                  /this->germanium_outer_radius);
    G4double total_angle = M_PI/2.0 -2.0*initial_angle;

    G4Cons* rounded_edge = new G4Cons("rounded_edge", base_inner_radius,
                                      base_outer_radius, tip_inner_radius, tip_outer_radius,
                                      con_half_length_z, initial_angle, total_angle);

    G4RotationMatrix* rotate_rounded_edge = new G4RotationMatrix;
    rotate_rounded_edge->rotateZ(-M_PI/2.0);

    G4ThreeVector move_rounded_edge(this->germanium_shift, -this->germanium_shift,
                                    this->germanium_length/2.0 -this->germanium_corner_cone_end_length/2.0
                                    + this->quarterDetectorCxnB);

    G4SubtractionSolid* germanium_rounded_edge = new G4SubtractionSolid("germanium_rounded_edge",
                                                                        germanium_rounded_corners, rounded_edge, rotate_rounded_edge,
                                                                        move_rounded_edge);

    // now we make the diagonal slices
    G4double half_chop_piece_width_x 	= this->germanium_width/2.0;
    G4double half_chop_piece_width_y 	= this->germanium_width/2.0;
    G4double half_chop_piece_length_z 	= (this->germanium_bent_length/cos(this->bent_end_angle))/2.0;
    G4Box* chop_piece = new G4Box("chop_piece", half_chop_piece_width_x,
                                  half_chop_piece_width_y, half_chop_piece_length_z);

    G4RotationMatrix* rotate_chop_piece1 = new G4RotationMatrix;
    rotate_chop_piece1->rotateX(-(this->bent_end_angle));

    G4ThreeVector move_chop_piece1(0, this->germanium_width/2.0 +(this->germanium_bent_length
                                                                  /tan(this->bent_end_angle) +this->germanium_width
                                                                  *cos(this->bent_end_angle))/2.0 -((this->germanium_bent_length
                                                                                                     /cos(this->bent_end_angle))/2.0) /sin(this->bent_end_angle),
                                   this->germanium_length/2.0 -this->germanium_bent_length
                                   +(this->germanium_bent_length +this->germanium_width
                                     *sin(this->bent_end_angle))/2.0);

    G4SubtractionSolid* chopped_germanium1 = new G4SubtractionSolid("chopped_germanium1",
                                                                    germanium_rounded_edge, chop_piece, rotate_chop_piece1, move_chop_piece1);


    G4RotationMatrix* rotate_chop_piece2 = new G4RotationMatrix;
    rotate_chop_piece2->rotateY(-(this->bent_end_angle));

    G4ThreeVector move_chop_piece2(-(this->germanium_width/2.0 +(this->germanium_bent_length
                                                                 /tan(this->bent_end_angle) +this->germanium_width
                                                                 *cos(this->bent_end_angle))/2.0 -((this->germanium_bent_length
                                                                                                    /cos(this->bent_end_angle))/2.0) /sin(this->bent_end_angle)), 0,
                                   this->germanium_length/2.0 -this->germanium_bent_length
                                   +(this->germanium_bent_length +this->germanium_width
                                     *sin(this->bent_end_angle))/2.0);

    G4SubtractionSolid* chopped_germanium2 = new G4SubtractionSolid("chopped_germanium2",
                                                                    chopped_germanium1, chop_piece, rotate_chop_piece2, move_chop_piece2);


    // now we make the hole that will go in the back of each quarter detector, /////////////////////////////////////////////////////////
    start_angle = 0.0*M_PI;
    final_angle = 2.0*M_PI;
    G4double hole_radius = this->germanium_hole_radius;
    G4double hole_half_length_z = (this->germanium_length -this->germanium_hole_dist_from_face)/2.0;

    G4Tubs* hole_tubs = new G4Tubs("hole_tubs", 0.0, hole_radius, hole_half_length_z, start_angle, final_angle);
    G4ThreeVector move_hole(this->germanium_shift, -(this->germanium_shift),-((this->germanium_hole_dist_from_face)));

    G4SubtractionSolid* chopped_germanium3 = new G4SubtractionSolid("chopped_germanium3",
                                                                    chopped_germanium2, hole_tubs, 0, move_hole);

    start_angle = 0.0*M_PI;
    final_angle = 2.0*M_PI;
    G4double dead_layer_radius = this->germanium_hole_radius + this->griffinDeadLayers[det][cry];
    G4double dead_layer_half_length_z = (this->germanium_length
                                         - this->germanium_hole_dist_from_face
                                         + this->griffinDeadLayers[det][cry])/2.0;

    // now dead layer
    G4Tubs* dead_layer_tubs = new G4Tubs("dead_layer_tubs", 0.0,
                                         dead_layer_radius, dead_layer_half_length_z, start_angle,final_angle);

    G4ThreeVector move_dead_layer(this->germanium_shift, -(this->germanium_shift),
                                  -((this->germanium_hole_dist_from_face
                                     - this->griffinDeadLayers[det][cry])/2.0));

    G4SubtractionSolid* chopped_germanium4 = new G4SubtractionSolid("chopped_germanium4",
                                                                    chopped_germanium3, dead_layer_tubs, 0, move_dead_layer);

    return chopped_germanium4;

}//end ::quarterSpecificDeadLayerDetector

///////////////////////////////////////////////////////////////////////
// methods used in ConstructNewSuppressorCasingWithShells
///////////////////////////////////////////////////////////////////////

G4SubtractionSolid* DetectionSystemGriffin::shellForBackSuppressorQuarter()
{
    G4double half_thickness_x = this->back_BGO_thickness/2.0 + this->suppressor_shell_thickness;
    G4double half_length_y = this->detector_total_width/4.0;
    G4double half_length_z = half_length_y;
    G4Box* quarter_suppressor_shell = new G4Box("quarter_suppressor_shell", half_thickness_x, half_length_y, half_length_z);

    G4double shell_hole_radius = this->cold_finger_outer_shell_radius;
    G4Tubs* shell_hole = new G4Tubs("shell_hole", 0, shell_hole_radius, half_thickness_x +1.0*cm, 0.0*M_PI, 2.0*M_PI);

    G4RotationMatrix* rotate_shell_hole = new G4RotationMatrix;
    rotate_shell_hole->rotateY(M_PI/2.0);
    G4ThreeVector move_shell_hole(0, -half_length_y, half_length_z);

    G4SubtractionSolid* quarter_suppressor_shell_with_hole = new G4SubtractionSolid("quarter_suppressor_shell_with_hole",
                                                                                    quarter_suppressor_shell, shell_hole, rotate_shell_hole, move_shell_hole);

    //now we need to cut out inner cavity, first we define the cavity
    half_thickness_x = this->back_BGO_thickness/2.0;
    half_length_y = this->detector_total_width/4.0 - this->suppressor_shell_thickness;
    half_length_z = half_length_y;
    G4Box* quarter_suppressor = new G4Box("quarter_suppressor", half_thickness_x, half_length_y, half_length_z);

    G4ThreeVector move_cut(0,0,0);

    // cut
    G4SubtractionSolid* quarter_suppressor_shell_with_hole_and_cavity = new G4SubtractionSolid("quarter_suppressor_shell_with_hole_and_cavity",
                                                                                               quarter_suppressor_shell_with_hole, quarter_suppressor, 0, move_cut);

    return quarter_suppressor_shell_with_hole_and_cavity;

}//end ::shellForBackSuppressorQuarter


G4SubtractionSolid* DetectionSystemGriffin::shellForFrontSlantSuppressor(G4String sidePosition)
{
    // Change some values to accomodate the shells
    // Replacement for this->side_suppressor_length
    G4double shell_side_suppressor_shell_length = this->side_suppressor_length + (this->suppressor_shell_thickness*2.0);

    G4double length_z         = shell_side_suppressor_shell_length;
    G4double length_y         = this->side_BGO_thickness + this->suppressor_shell_thickness*2.0;
    G4double length_longer_x  = this->detector_total_width/2.0 +this->BGO_can_seperation
            + this->side_BGO_thickness + this->suppressor_shell_thickness*2.0;
    G4double length_shorter_x = length_longer_x -this->side_BGO_thickness
            - this->suppressor_shell_thickness*2.0;

    G4Trap* suppressor_shell = new G4Trap("suppressor_shell", length_z, length_y, length_longer_x, length_shorter_x);

    G4double half_length_x      = length_longer_x/2.0 +1.0*cm;
    G4double half_thickness_y   = this->side_BGO_thickness/2.0 + this->suppressor_shell_thickness;
    G4double half_length_z      = ((this->side_BGO_thickness + this->suppressor_shell_thickness*2.0
                                    - this->BGO_chopped_tip)/sin(this->bent_end_angle))/2.0;

    G4Box* chopping_shell_box = new G4Box("chopping_shell_box", half_length_x, half_thickness_y, half_length_z);

    G4RotationMatrix* rotate_chopping_shell_box = new G4RotationMatrix;

    G4double y0 = this->side_BGO_thickness/2.0 - this->suppressor_shell_thickness * sin(this->bent_end_angle)
            - this->BGO_chopped_tip - 0.5 * sqrt(pow((2.0 * half_length_z), 2.0)
                                                 + pow((2.0 * half_thickness_y), 2.0)) * cos(M_PI/2.0 - this->bent_end_angle - atan(half_thickness_y/half_length_z)) ;

    G4double z0 = length_z/2.0 - 0.5 * sqrt(pow((2.0 * half_length_z), 2.0) + pow((2.0 * half_thickness_y), 2.0)) * sin(M_PI/2.0
                                                                                                                        - this->bent_end_angle - atan(half_thickness_y/half_length_z)) ;
    G4ThreeVector move_chopping_shell_box ;

    if(sidePosition == "left")
    {
        rotate_chopping_shell_box->rotateX(this->bent_end_angle);
        move_chopping_shell_box = G4ThreeVector( 0, y0, z0 ) ;
    }
    else if(sidePosition == "right")
    {
        rotate_chopping_shell_box->rotateX(-this->bent_end_angle);
        move_chopping_shell_box = G4ThreeVector( 0, y0, -z0 ) ;
    }

    G4SubtractionSolid* side_suppressor_shell = new G4SubtractionSolid("side_suppressor_shell",
                                                                       suppressor_shell, chopping_shell_box, rotate_chopping_shell_box, move_chopping_shell_box);

    G4SubtractionSolid* side_suppressor_shell_with_cavity ;

    // cut out cavity
    if(sidePosition == "left")
    {
        G4ThreeVector move_cut(-(this->suppressor_shell_thickness + ( this->extra_cut_length/2.0 - this->suppressor_shell_thickness)/2.0 ), 0, (this->suppressor_shell_thickness/2.0) );

        side_suppressor_shell_with_cavity = new G4SubtractionSolid( "side_suppressor_shell_with_cavity", side_suppressor_shell,
                                                                    this->frontSlantSuppressor("left", true), 0, move_cut); // chopping, left from frontSlantSuppressor.
    }
    else if(sidePosition == "right")
    {
        G4ThreeVector move_cut(-(this->suppressor_shell_thickness + ( this->extra_cut_length/2.0 - this->suppressor_shell_thickness)/2.0 ), 0, -(this->suppressor_shell_thickness/2.0) );

        side_suppressor_shell_with_cavity = new G4SubtractionSolid( "side_suppressor_shell_with_cavity", side_suppressor_shell,
                                                                    this->frontSlantSuppressor("right", true), 0, move_cut); // chopping, right from frontSlantSuppressor.
    }

    return side_suppressor_shell_with_cavity;
} // end ::shellForFrontSlantSuppressor


G4SubtractionSolid* DetectionSystemGriffin::shellForSuppressorExtension(G4String sidePosition)
{

    // Replacement for this->suppressor_extension_length
    G4double shell_suppressor_shell_extension_length = this->suppressor_extension_length
            + (this->suppressor_shell_thickness*2.0)*(1.0/tan(this->bent_end_angle)
                                                      - tan(this->bent_end_angle));

    G4double thickness_z = this->suppressor_extension_thickness + this->suppressor_shell_thickness*2.0;
    G4double length_y = shell_suppressor_shell_extension_length;

    G4double longer_length_x =  ( this->suppressor_back_radius +this->bent_end_length +(this->BGO_can_seperation
                                                                                        + this->side_BGO_thickness
                                                                                        + this->suppressor_shell_thickness*2.0)
                                  / tan(this->bent_end_angle)
                                  - (this->suppressor_extension_thickness
                                     + this->suppressor_shell_thickness*2.0)
                                  * sin(this->bent_end_angle))*tan(this->bent_end_angle);

    G4double shorter_length_x = (this->suppressor_forward_radius + this->hevimet_tip_thickness) * sin(this->bent_end_angle) ;

    G4Trap* uncut_extension_shell = new G4Trap("uncut_extension_shell", thickness_z, length_y, longer_length_x, shorter_length_x);

    // because these pieces are rotated in two planes, there are multiple angles that need to be calculated to make sure
    // all of the extensions join up

    G4double beta = atan((longer_length_x -shorter_length_x)/(length_y));
    G4double phi  = atan(1/cos(this->bent_end_angle));

    G4double chopping_half_length_x = ( thickness_z / sin(phi) ) / 2.0;
    G4double chopping_half_length_y = length_y / ( 2.0 * cos(beta) );
    G4double chopping_half_length_z = chopping_half_length_x;
    G4double y_angle = -beta;
    G4double x_angle = 0.0;
    G4double z_angle ;

    if( sidePosition == "left" )
        z_angle = M_PI/2.0 - phi ;
    else if( sidePosition == "right" )
        z_angle = phi - M_PI/2.0 ;


    G4Para* chopping_para_shell = new G4Para("chopping_para_shell", chopping_half_length_x,
                                             chopping_half_length_y, chopping_half_length_z, y_angle, z_angle, x_angle);

    G4ThreeVector move_para_shell(((longer_length_x -shorter_length_x)/2.0 +shorter_length_x)/2.0
                                  + chopping_half_length_x -chopping_half_length_x*cos(phi), 0.0, 0.0);

    G4SubtractionSolid* extension_shell = new G4SubtractionSolid("extension_shell",
                                                                 uncut_extension_shell, chopping_para_shell, 0, move_para_shell);

    G4ThreeVector move_cut ;
    G4SubtractionSolid* extension_suppressor_shell_with_cavity ;
    G4SubtractionSolid* extension_suppressor_shell_with_cavity_pre;

    // cut out cavity ----
    // We make two cuts, one "crude" one from a G4Trap, and a more refined cut using a G4SubtractionSolid. There are two reasons why we make two cuts.
    // Firstly, the shape of these suppressors are strange. To make the inside cut large enough to have an opening such that the left and right suppressors touch would
    // require some clever math, and that's a pain. The second (less lazy) reason is Geant4 doesn't like to make a G4SubtractionSolid out of two G4SubtractionSolids.
    // My experience is that the visulization will fail if two sides of the G4SubtractionSolids are close, that is do not have a LARGE (>~10mm) overlap. To get around this,
    // I simply make a "crude" first cut using a G4Trap.
    // The cuts are slightly larger than suppressor to allow for small gaps between the shell and the suppressor.
    if( sidePosition == "left" )
    {
        move_cut = G4ThreeVector(-30.0*mm,0.40*mm,0.25*mm);
        extension_suppressor_shell_with_cavity_pre = new G4SubtractionSolid("extension_suppressor_shell_with_cavity_pre", extension_shell, this->sideSuppressorExtensionUncut(), 0, (move_cut)/2.0 ); // Left, Chopping sideSuppressorExtension

        move_cut = G4ThreeVector(-1.15*mm,0.40*mm,0.25*mm);
        extension_suppressor_shell_with_cavity = new G4SubtractionSolid("extension_suppressor_shell_with_cavity", extension_suppressor_shell_with_cavity_pre, this->sideSuppressorExtension( "left", true ), 0, (move_cut)/2.0 );
    }
    else if( sidePosition == "right" )
    {
        move_cut = G4ThreeVector(-30.0*mm,0.40*mm,-0.25*mm);
        extension_suppressor_shell_with_cavity_pre = new G4SubtractionSolid("extension_suppressor_shell_with_cavity_pre", extension_shell, this->sideSuppressorExtensionUncut(), 0, (move_cut)/2.0 ); // Left, Chopping sideSuppressorExtension

        move_cut = G4ThreeVector(-1.15*mm,0.40*mm,-0.25*mm);
        extension_suppressor_shell_with_cavity = new G4SubtractionSolid("extension_suppressor_shell_with_cavity", extension_suppressor_shell_with_cavity_pre, this->sideSuppressorExtension( "right", true ), 0, (move_cut)/2.0 );
    }

    return extension_suppressor_shell_with_cavity;
} // end ::shellForSuppressorExtension


// Back to making the suppressor volumes themselves
G4SubtractionSolid* DetectionSystemGriffin::backSuppressorQuarter()
{
    G4double half_thickness_x = this->back_BGO_thickness/2.0;
    G4double half_length_y = this->detector_total_width/4.0 - this->suppressor_shell_thickness;
    G4double half_length_z = half_length_y;
    G4Box* quarter_suppressor = new G4Box("quarter_suppressor", half_thickness_x, half_length_y, half_length_z);

    G4double hole_radius = this->cold_finger_outer_shell_radius;
    G4Tubs* hole = new G4Tubs("hole", 0, hole_radius, half_thickness_x +1.0*cm, 0.0*M_PI, 2.0*M_PI);

    G4RotationMatrix* rotate_hole = new G4RotationMatrix;
    rotate_hole->rotateY(M_PI/2.0);
    G4ThreeVector move_hole(0, -half_length_y, half_length_z);

    G4SubtractionSolid* quarter_suppressor_with_hole = new G4SubtractionSolid("quarter_suppressor_with_hole",
                                                                              quarter_suppressor, hole, rotate_hole, move_hole);

    return quarter_suppressor_with_hole;

}//end ::backSuppressorQuarter


G4SubtractionSolid* DetectionSystemGriffin::frontSlantSuppressor( G4String sidePosition, G4bool choppingSuppressor )
{
    // If chop is true, it will be a chopping suppressor.

    G4double length_z     = this->side_suppressor_length ;
    G4double length_y     = this->side_BGO_thickness ;
    G4double length_longer_x  = 0 ;

    if( choppingSuppressor )
        length_longer_x  = this->detector_total_width / 2.0 + this->BGO_can_seperation + this->side_BGO_thickness + this->extra_cut_length / 2.0 ;
    else
        length_longer_x  = this->detector_total_width / 2.0 + this->BGO_can_seperation + this->side_BGO_thickness;
    
    G4double length_shorter_x   = length_longer_x - this->side_BGO_thickness;

    G4Trap* suppressor = new G4Trap( "suppressor", length_z, length_y, length_longer_x, length_shorter_x ) ;

    G4double half_length_x  = length_longer_x / 2.0 + 1.0*cm ;
    G4double half_thickness_y   = this->side_BGO_thickness / 2.0;
    G4double half_length_z  = (( this->side_BGO_thickness - this->BGO_chopped_tip ) / sin( this->bent_end_angle ) ) / 2.0;

    G4Box* chopping_box = new G4Box("chopping_box", half_length_x, half_thickness_y, half_length_z);

    G4RotationMatrix* rotate_chopping_box = new G4RotationMatrix;
    G4ThreeVector move_chopping_box ;
    G4double y0, z0 ;

    y0 =  this->side_BGO_thickness / 2.0 - this->BGO_chopped_tip - 0.5 * sqrt( pow( (2.0 * half_length_z), 2.0 )
                                                                               + pow( (2.0 * half_thickness_y), 2.0 ) ) * cos( M_PI/2.0 - this->bent_end_angle - atan( half_thickness_y / half_length_z ) ) ;

    z0 =  length_z / 2.0 - 0.5 * sqrt( pow( (2.0 * half_length_z), 2.0 ) + pow( (2.0 * half_thickness_y), 2.0) )
            * sin(M_PI/2.0 - this->bent_end_angle - atan( half_thickness_y / half_length_z ) ) ;

    if( sidePosition == "left" )
    {
        rotate_chopping_box->rotateX( this->bent_end_angle ) ;
        move_chopping_box = G4ThreeVector( 0, y0, z0 ) ;
    }
    else if( sidePosition == "right" )
    {
        rotate_chopping_box->rotateX( -this->bent_end_angle ) ;
        move_chopping_box = G4ThreeVector( 0, y0, -z0 ) ;
    }
    
    G4SubtractionSolid* side_suppressor = new G4SubtractionSolid("side_suppressor",
                                                                 suppressor, chopping_box, rotate_chopping_box, move_chopping_box);

    return side_suppressor;

}// end ::frontSlantSuppressor


G4SubtractionSolid* DetectionSystemGriffin::sideSuppressorExtension(G4String sidePosition, G4bool choppingSuppressor) 
{

    G4double thickness_z      =   this->suppressor_extension_thickness ;
    G4double length_y         =   this->suppressor_extension_length ;

    G4double longer_length_x  =   (( this->suppressor_back_radius + this->bent_end_length
                                    + (this->BGO_can_seperation
                                       + this->side_BGO_thickness) / tan(this->bent_end_angle)
                                    - this->suppressor_extension_thickness
                                    * sin(this->bent_end_angle)) * tan(this->bent_end_angle) - this->BGO_can_seperation * 2.0) ;  // - this->BGO_can_seperation

    G4double shorter_length_x =   (( this->suppressor_forward_radius + this->hevimet_tip_thickness)
            * sin(this->bent_end_angle) - this->BGO_can_seperation * 2.0) ; // - this->BGO_can_seperation

    if( choppingSuppressor ) {
        thickness_z += 0.1*mm;
        length_y += 0.1*mm;
    }

    G4Trap* uncut_extension   = new G4Trap("uncut_extension", thickness_z, length_y, longer_length_x, shorter_length_x);

    // because these pieces are rotated in two planes, there are multiple angles that need to be calculated to make sure
    // all of the extensions join up
    G4double beta = atan( ( longer_length_x - shorter_length_x ) / ( length_y ) ) ;
    G4double phi  = atan( 1 / cos( this->bent_end_angle ) ) ;

    G4double chopping_half_length_x = thickness_z / ( 2.0 * sin( phi ) ) ;
    G4double chopping_half_length_y = length_y / ( 2.0 * cos( beta ) ) ;
    G4double chopping_half_length_z = chopping_half_length_x ;

    G4double x_angle = 0.0 ;
    G4double y_angle = -beta ;
    G4double z_angle = phi - M_PI/2.0 ;

    if( sidePosition == "left" )
        z_angle *= -1 ;

    G4Para* chopping_para = new G4Para("chopping_para", chopping_half_length_x,
                                       chopping_half_length_y, chopping_half_length_z,
                                       y_angle, z_angle, x_angle);

    G4ThreeVector move_para(((longer_length_x - shorter_length_x ) / 2.0
                             + shorter_length_x ) / 2.0 + chopping_half_length_x
                            - chopping_half_length_x * cos( phi ), 0.0, 0.0);

    G4SubtractionSolid* right_extension = new G4SubtractionSolid( "right_extension", uncut_extension, chopping_para, 0, move_para );

    return right_extension;

}// end ::sideSuppressorExtension()


G4Trap* DetectionSystemGriffin::sideSuppressorExtensionUncut()
{

    G4double thickness_z      =   this->suppressor_extension_thickness ;
    G4double length_y         =   this->suppressor_extension_length ;

    G4double longer_length_x  =   (( this->suppressor_back_radius + this->bent_end_length
                                    + (this->BGO_can_seperation
                                       + this->side_BGO_thickness) / tan(this->bent_end_angle)
                                    - this->suppressor_extension_thickness
                                    * sin(this->bent_end_angle)) * tan(this->bent_end_angle) - this->BGO_can_seperation * 2.0) ;  // - this->BGO_can_seperation

    G4double shorter_length_x =   (( this->suppressor_forward_radius + this->hevimet_tip_thickness)
            * sin(this->bent_end_angle) - this->BGO_can_seperation * 2.0) ; // - this->BGO_can_seperation

    // Similar to the "true" choppingSuppressor, add a bit on the z and y lengths to give a bit of a gap.
    thickness_z += 0.1*mm;
    length_y += 0.1*mm;

    G4Trap* uncut_extension   = new G4Trap("uncut_extension", thickness_z, length_y, longer_length_x, shorter_length_x);

    return uncut_extension;

}// end ::sideSuppressorExtensionUncut()


///////////////////////////////////////////////////////////////////////
// methods used in ConstructNewHeavyMet
///////////////////////////////////////////////////////////////////////
G4SubtractionSolid* DetectionSystemGriffin::newHeavyMet()
{

    G4double half_thickness_z = ((this->hevimet_tip_thickness/cos(this->hevimet_tip_angle))
                                 * cos(this->bent_end_angle -this->hevimet_tip_angle)
                                 + this->suppressor_forward_radius * tan(this->hevimet_tip_angle)
                                 * sin(this->bent_end_angle))/2.0;


    G4double half_length_shorter_x  = this->suppressor_forward_radius * sin(this->bent_end_angle);
    G4double half_length_longer_x 	= half_length_shorter_x +2.0*half_thickness_z*tan(this->bent_end_angle);
    G4double  half_length_shorter_y = half_length_shorter_x;
    G4double half_length_longer_y 	= half_length_longer_x;

    G4Trd* uncut_hevimet = new G4Trd("uncut_hevimet", half_length_longer_x,
                                     half_length_shorter_x, half_length_longer_y, half_length_shorter_y, half_thickness_z);

    G4double half_shorter_length_x = half_length_shorter_x
            - this->suppressor_forward_radius * tan( this->hevimet_tip_angle ) * cos(this->bent_end_angle);

    G4double half_height_z =  ( 2.0 * half_thickness_z + ( half_length_longer_x
                                                           - ( ( this->suppressor_forward_radius + this->hevimet_tip_thickness )
                                                               * tan( this->hevimet_tip_angle ) ) / cos( this->bent_end_angle )
                                                           - half_shorter_length_x)*tan( this->bent_end_angle )) / 2.0;

    G4double half_longer_length_x 	= half_shorter_length_x + 2.0 * half_height_z/tan( this->bent_end_angle ) ;
    G4double half_shorter_length_y 	= half_shorter_length_x;
    G4double half_longer_length_y 	= half_longer_length_x;

    G4Trd* intersector = new G4Trd( "intersector", half_shorter_length_x,
                                    half_longer_length_x, half_shorter_length_y, half_longer_length_y, half_height_z );

    G4Trd* chopper = new G4Trd( "chopper", half_shorter_length_x, half_longer_length_x,
                                half_shorter_length_y, half_longer_length_y, half_height_z );

    G4ThreeVector move_chopper(0.0, 0.0, half_height_z + half_thickness_z
                               - this->suppressor_forward_radius * tan(this->hevimet_tip_angle)
                               * sin(this->bent_end_angle) ) ;

    G4SubtractionSolid* chopped_hevimet = new G4SubtractionSolid("chopped_hevimet",
                                                                 uncut_hevimet, chopper, 0, move_chopper);

    G4ThreeVector move_intersector(0.0, 0.0, -half_height_z + half_thickness_z);

    G4IntersectionSolid* intersected_hevimet = new G4IntersectionSolid("intersected_hevimet",
                                                                       chopped_hevimet, intersector, 0, move_intersector);

    G4double longer_half_length = half_length_longer_x -((this->suppressor_forward_radius
                                                          + this->hevimet_tip_thickness)*tan(this->hevimet_tip_angle)) / cos(this->bent_end_angle) ;

    G4double shorter_half_length = longer_half_length -2.0*(half_thickness_z +0.1*mm)
            * tan(this->bent_end_angle -this->hevimet_tip_angle);

    G4Trd* middle_chopper = new G4Trd("middle_chopper", longer_half_length, shorter_half_length,
                                      longer_half_length, shorter_half_length, half_thickness_z +0.1*mm);

    G4SubtractionSolid* hevimet = new G4SubtractionSolid("hevimet", intersected_hevimet, middle_chopper);

    return hevimet;

}//end ::newHeavyMet


