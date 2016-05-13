#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemSceptar.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

DetectionSystemSceptar::DetectionSystemSceptar() :
    // LogicalVolumes
    square_mylar_log(0),
    angled_mylar_log(0),
    square_scintillator_log(0),
    angled_scintillator_log(0),
    Delrin_shell_log(0),
    Delrin_shell2_log(0),
    Hevimet_shell_log(0)
{

    /////////////////////////////////////////////////////////////////////
    // Detector Properties
    /////////////////////////////////////////////////////////////////////
    this->separate_hemispheres = 1.0*mm;
    this->convert = 0.0254*m; //Needed to this->convert inches to meters

    this->square_scintillator_length = 1.194*this->convert;
    this->square_scintillator_width = 0.449*this->convert;
    this->square_scintillator_thickness = this->convert/16.0;
    this->angled_scintillator_length  = 1.003*this->convert;
    this->angled_scintillator_long_width = 1.045*this->convert;
    this->angled_scintillator_short_width = 0.690*this->convert;
    this->angled_scintillator_thickness = this->convert/16.0;
    this->mylar_thickness = 0.01*mm;
    this->scint_gap = 0.75*mm;
    this->scint_angle1 = 36.0*deg;
    this->scint_angle_move = 18.0*deg;

    //*********This angle is the tapering along the length (10.05)
    this->scint_angle2 = atan((this->angled_scintillator_long_width -this->angled_scintillator_short_width)/(2.0*this->angled_scintillator_length));

    //*********This angle is the slight tapering of the sides (34.765))
    this->scint_angle3 = acos(cos(this->scint_angle1)/cos(this->scint_angle2));

    //%%%%% This is the angle that defines the packed
    //%%%%% configuration of the angled_scints
    //%%%%%
    this->scint_angle4 = asin((this->angled_scintillator_long_width -this->angled_scintillator_short_width)/(2.0*this->angled_scintillator_length*tan(this->scint_angle1)));

    //%%%%% This is the radial distance from the center of the
    //%%%%% pentagon formed by the square_detectors to the
    //%%%%% center of any one of the square_detectors in that pentagon
    //%%%%% The gap between scintillators and the mylar coating are taken into account
    //%%%%%
    this->square_scint_radial_distance = ((this->square_scintillator_length +2.0*this->mylar_thickness +this->scint_gap/cos(this->scint_angle1))*tan((M_PI/2.0)-this->scint_angle1)-(this->square_scintillator_thickness+2.0*this->mylar_thickness))/2.0;

    //%%%%% This is the radial distance from the axis of the
    //%%%%% pentagon formed by the angled_detectors to the
    //%%%%% center of any one of the angled_detectors in that pentagon
    //%%%%% The gap between the scintillators and the mylar coating are taken into account
    //%%%%%
    this->angled_scint_radial_distance = 0.5*(((this->angled_scintillator_short_width +2.0*this->mylar_thickness +this->scint_gap/cos(this->scint_angle1))/tan(this->scint_angle1)) +(this->angled_scintillator_length+2.0*this->mylar_thickness)*sin(this->scint_angle4)-(this->angled_scintillator_thickness +2.0*this->mylar_thickness)*cos(this->scint_angle4));

    //%%%%% This is the distance from the origin, along-y that
    //%%%%% the angled_scints need to be moved to give
    //%%%%% proper alignment with the square_scints
    //%%%%% The mylar coating is taken into account
    //%%%%%
    this->angled_scint_move_back = this->square_scintillator_width -(this->square_scint_radial_distance -0.5*this->square_scintillator_thickness -0.5*this->angled_scintillator_long_width/tan(this->scint_angle1))/tan(atan((this->square_scint_radial_distance -0.5*this->square_scintillator_thickness)/this->square_scintillator_width))+0.5*this->angled_scintillator_length*cos(this->scint_angle4);

    this->Delrin_inner_radius = 3.275*this->convert;
    this->Delrin_outer_radius = 3.475*this->convert;
    this->Delrin2_inner_radius = 8.90*cm;
    this->Delrin2_outer_radius = 9.90*cm;
    this->Hevimet_inner_radius = 9.906*cm;
    this->Hevimet_outer_radius = 12.446*cm;
    this->Delrin_hole_radius = 1.25*this->convert;
}

DetectionSystemSceptar::~DetectionSystemSceptar()
{
    // LogicalVolumes
    delete square_mylar_log;
    delete angled_mylar_log;
    delete square_scintillator_log;
    delete angled_scintillator_log;
    delete Delrin_shell_log;
    delete Delrin_shell2_log;
    delete Hevimet_shell_log;
}

G4int DetectionSystemSceptar::Build()
{
    // Build assembly volume
    G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    this->assembly = myAssembly;
    G4AssemblyVolume* myAssemblySquare = new G4AssemblyVolume();
    this->assemblySquare = myAssemblySquare;
    G4AssemblyVolume* myAssemblyAngled = new G4AssemblyVolume();
    this->assemblyAngled = myAssemblyAngled;
    G4AssemblyVolume* myAssemblySquareSD = new G4AssemblyVolume();
    this->assemblySquareSD = myAssemblySquareSD;
    G4AssemblyVolume* myAssemblyAngledSD = new G4AssemblyVolume();
    this->assemblyAngledSD = myAssemblyAngledSD;

    ConstructScintillator();
    //  ConstructDelrinShell();
    //  Construct2ndDelrinShell();
    //  ConstructHevimetShell();

    return 1;
}

G4int DetectionSystemSceptar::PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detectorNumber)
{

    //G4RotationMatrix* rotate_null = new G4RotationMatrix;
    G4ThreeVector move_null(0.0,0.0,0.0);

    G4double extra = 2.0*this->mylar_thickness;

    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4ThreeVector move;

    G4RotationMatrix* rotate_square_scint1 = new G4RotationMatrix;
    rotate_square_scint1->rotateX(M_PI/2.0);
    rotate_square_scint1->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint1(	-this->square_scint_radial_distance ,
                                        0,
                                        -(this->square_scintillator_width +extra + this->separate_hemispheres)/2.0);

    G4RotationMatrix* rotate_square_scint2 = new G4RotationMatrix;
    rotate_square_scint2->rotateY(-1.0*(2.0*this->scint_angle1));
    rotate_square_scint2->rotateX(M_PI/2.0);
    rotate_square_scint2->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint2(	-this->square_scint_radial_distance*cos(2.0*this->scint_angle1),
                                        this->square_scint_radial_distance*sin(2.0*this->scint_angle1),
                                        -(this->square_scintillator_width +extra + this->separate_hemispheres) / 2.0 ) ;

    G4RotationMatrix* rotate_square_scint3 = new G4RotationMatrix;
    rotate_square_scint3->rotateY(-1.0*(M_PI -this->scint_angle1));
    rotate_square_scint3->rotateX(M_PI/2.0);
    rotate_square_scint3->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint3(	this->square_scint_radial_distance*cos(this->scint_angle1),
                                        this->square_scint_radial_distance*sin(this->scint_angle1),
                                        -(this->square_scintillator_width +extra + this->separate_hemispheres) / 2.0 ) ;

    G4RotationMatrix* rotate_square_scint4 = new G4RotationMatrix;
    rotate_square_scint4->rotateY(-1.0*(this->scint_angle1 -M_PI));
    rotate_square_scint4->rotateX(M_PI/2.0);
    rotate_square_scint4->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint4(	this->square_scint_radial_distance*cos(this->scint_angle1),
                                        -this->square_scint_radial_distance*sin(this->scint_angle1),
                                        -(this->square_scintillator_width +extra + this->separate_hemispheres) / 2.0 ) ;

    G4RotationMatrix* rotate_square_scint5 = new G4RotationMatrix;
    rotate_square_scint5->rotateY(-1.0*(-2.0*this->scint_angle1));
    rotate_square_scint5->rotateX(M_PI/2.0);
    rotate_square_scint5->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint5(	-this->square_scint_radial_distance*cos(2.0*this->scint_angle1),
                                        -this->square_scint_radial_distance*sin(2.0*this->scint_angle1),
                                        -(this->square_scintillator_width +extra + this->separate_hemispheres) / 2.0 ) ;

    G4RotationMatrix* rotate_square_scint6 = new G4RotationMatrix;
    rotate_square_scint6->rotateY(M_PI);
    rotate_square_scint6->rotateX(M_PI/2.0);
    rotate_square_scint6->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint6;
    move_square_scint6 = - 1.0 * move_square_scint1;

    G4RotationMatrix* rotate_square_scint7 = new G4RotationMatrix;
    rotate_square_scint7->rotateY(-1.0*(2.0*this->scint_angle1 -M_PI));
    rotate_square_scint7->rotateX(M_PI/2.0);
    rotate_square_scint7->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint7;
    move_square_scint7 = -1.0*move_square_scint2;

    G4RotationMatrix* rotate_square_scint8 = new G4RotationMatrix;
    rotate_square_scint8->rotateY(-1.0*(-this->scint_angle1));
    rotate_square_scint8->rotateX(M_PI/2.0);
    rotate_square_scint8->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint8;
    move_square_scint8 = -1.0*move_square_scint3;

    G4RotationMatrix* rotate_square_scint9 = new G4RotationMatrix;
    rotate_square_scint9->rotateY(-1.0*(this->scint_angle1));
    rotate_square_scint9->rotateX(M_PI/2.0);
    rotate_square_scint9->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint9;
    move_square_scint9 = -1.0*move_square_scint4;

    G4RotationMatrix* rotate_square_scint10 = new G4RotationMatrix;
    rotate_square_scint10->rotateY(-1.0*(-2.0*this->scint_angle1 -M_PI));
    rotate_square_scint10->rotateX(M_PI/2.0);
    rotate_square_scint10->rotateZ(M_PI/2.0);
    G4ThreeVector move_square_scint10;
    move_square_scint10 = -1.0*move_square_scint5;

    G4RotationMatrix* rotate_angled_scint1 = new G4RotationMatrix;
    rotate_angled_scint1->rotateX(-1.0*(this->scint_angle4 -M_PI/2.0));
    rotate_angled_scint1->rotateX(M_PI/2.0);
    rotate_angled_scint1->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint1(	-this->angled_scint_radial_distance,
                                        0,
                                        -this->angled_scint_move_back);

    G4RotationMatrix* rotate_angled_scint2 = new G4RotationMatrix;
    rotate_angled_scint2->rotateX(-1.0*(this->scint_angle4));
    rotate_angled_scint2->rotateZ(-1.0*(-2.0*this->scint_angle1));
    rotate_angled_scint2->rotateX(M_PI);
    rotate_angled_scint2->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint2(	-this->angled_scint_radial_distance*cos(2.0*this->scint_angle1),
                                        this->angled_scint_radial_distance*sin(2.0*this->scint_angle1),
                                        -this->angled_scint_move_back);

    G4RotationMatrix* rotate_angled_scint3 = new G4RotationMatrix;
    rotate_angled_scint3->rotateX(-1.0*(this->scint_angle4));
    rotate_angled_scint3->rotateZ(-1.0*(this->scint_angle1 -M_PI));
    rotate_angled_scint3->rotateX(M_PI);
    rotate_angled_scint3->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint3(	this->angled_scint_radial_distance*cos(this->scint_angle1),
                                        this->angled_scint_radial_distance*sin(this->scint_angle1),
                                        -this->angled_scint_move_back);

    G4RotationMatrix* rotate_angled_scint4 = new G4RotationMatrix;
    rotate_angled_scint4->rotateX(-1.0*(this->scint_angle4));
    rotate_angled_scint4->rotateZ(-1.0*(M_PI -this->scint_angle1));
    rotate_angled_scint4->rotateX(M_PI);
    rotate_angled_scint4->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint4(	this->angled_scint_radial_distance*cos(this->scint_angle1),
                                        -this->angled_scint_radial_distance*sin(this->scint_angle1),
                                        -this->angled_scint_move_back);

    G4RotationMatrix* rotate_angled_scint5 = new G4RotationMatrix;
    rotate_angled_scint5->rotateX(-1.0*(this->scint_angle4));
    rotate_angled_scint5->rotateZ(-1.0*(2.0*this->scint_angle1));
    rotate_angled_scint5->rotateX(M_PI);
    rotate_angled_scint5->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint5(	-this->angled_scint_radial_distance*cos(2.0*this->scint_angle1),
                                        -this->angled_scint_radial_distance*sin(2.0*this->scint_angle1),
                                        -this->angled_scint_move_back);

    G4RotationMatrix* rotate_angled_scint6 = new G4RotationMatrix;
    rotate_angled_scint6->rotateX(-1.0*(this->scint_angle4 +M_PI/2.0));
    rotate_angled_scint6->rotateX(M_PI/2.0);
    rotate_angled_scint6->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint6;
    move_angled_scint6 = -1.0*move_angled_scint1;

    G4RotationMatrix* rotate_angled_scint7 = new G4RotationMatrix;
    rotate_angled_scint7->rotateX(-1.0*(this->scint_angle4));
    rotate_angled_scint7->rotateZ(-1.0*(2.0*this->scint_angle1));
    rotate_angled_scint7->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint7;
    move_angled_scint7 = -1.0*move_angled_scint2;

    G4RotationMatrix* rotate_angled_scint8 = new G4RotationMatrix;
    rotate_angled_scint8->rotateX(-1.0*(this->scint_angle4));
    rotate_angled_scint8->rotateZ(-1.0*(M_PI -this->scint_angle1));
    rotate_angled_scint8->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint8;
    move_angled_scint8 = -1.0*move_angled_scint3;

    G4RotationMatrix* rotate_angled_scint9 = new G4RotationMatrix;
    rotate_angled_scint9->rotateX(-1.0*(this->scint_angle4));
    rotate_angled_scint9->rotateZ(-1.0*(this->scint_angle1 -M_PI));
    rotate_angled_scint9->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint9;
    move_angled_scint9 = -1.0*move_angled_scint4;

    G4RotationMatrix* rotate_angled_scint10 = new G4RotationMatrix;
    rotate_angled_scint10->rotateX(-1.0*(this->scint_angle4));
    rotate_angled_scint10->rotateZ(-1.0*(-2.0*this->scint_angle1));
    rotate_angled_scint10->rotateZ(M_PI/2.0);
    G4ThreeVector move_angled_scint10;
    move_angled_scint10 = -1.0*move_angled_scint5;

    for(G4int detector_number = 0; detector_number < detectorNumber; detector_number++)
    {
        if(detector_number == 0)
        {
            rotate = rotate_square_scint6;
            move = move_square_scint6;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }
        else if(detector_number == 1)
        {
            rotate = rotate_square_scint7;
            move = move_square_scint7;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }
        else if(detector_number == 2)
        {
            rotate = rotate_square_scint8;
            move = move_square_scint8;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }
        else if(detector_number == 3)
        {
            rotate = rotate_square_scint9;
            move = move_square_scint9;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }
        else if(detector_number == 4)
        {
            rotate = rotate_square_scint10;
            move = move_square_scint10;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }
        else if(detector_number == 5)
        {
            rotate = rotate_angled_scint6;
            move = move_angled_scint6;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }
        else if(detector_number == 6)
        {
            rotate = rotate_angled_scint7;
            move = move_angled_scint7;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }
        else if(detector_number == 7)
        {
            rotate = rotate_angled_scint8;
            move = move_angled_scint8;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }
        else if(detector_number == 8)
        {
            rotate = rotate_angled_scint9;
            move = move_angled_scint9;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }
        else if(detector_number == 9)
        {
            rotate = rotate_angled_scint10;
            move = move_angled_scint10;
            move.rotateZ(this->scint_angle_move);
            rotate->rotateZ(this->scint_angle_move);
        }

        else if(detector_number == 10)
        {
            rotate = rotate_square_scint1;
            move = move_square_scint1;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }
        else if(detector_number == 11)
        {
            rotate = rotate_square_scint2;
            move = move_square_scint2;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }
        else if(detector_number == 12)
        {
            rotate = rotate_square_scint3;
            move = move_square_scint3;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }
        else if(detector_number == 13)
        {
            rotate = rotate_square_scint4;
            move = move_square_scint4;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }
        else if(detector_number == 14)
        {
            rotate = rotate_square_scint5;
            move = move_square_scint5;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }
        else if(detector_number == 15)
        {
            rotate = rotate_angled_scint1;
            move = move_angled_scint1;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }
        else if(detector_number == 16)
        {
            rotate = rotate_angled_scint2;
            move = move_angled_scint2;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }
        else if(detector_number == 17)
        {
            rotate = rotate_angled_scint3;
            move = move_angled_scint3;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }
        else if(detector_number == 18)
        {
            rotate = rotate_angled_scint4;
            move = move_angled_scint4;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }
        else if(detector_number == 19)
        {
            rotate = rotate_angled_scint5;
            move = move_angled_scint5;
            move.rotateZ(this->scint_angle_move+this->scint_angle1);
            rotate->rotateZ(this->scint_angle_move+this->scint_angle1);
        }



        if( (detector_number < 5) || (detector_number >= 10 && detector_number < 15) )
        {
            assemblySquareSD->MakeImprint(exp_hall_log, move, rotate, 0);
            assemblySquare->MakeImprint(exp_hall_log, move, rotate, 0);
        }
        if( (detector_number >= 5 && detector_number < 10) || (detector_number >= 15) )
        {
            assemblyAngledSD->MakeImprint(exp_hall_log, move, rotate, 0);
            assemblyAngled->MakeImprint(exp_hall_log, move, rotate, 0);
        }

    }

    return 1;
}

G4int DetectionSystemSceptar::ConstructScintillator()
{
    G4Material* material_mylar = G4Material::GetMaterial("Mylar");
    if( !material_mylar ) {
        G4cout << " ----> Material " << "Mylar" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }
    G4Material* material_bc404 = G4Material::GetMaterial("BC404");
    if( !material_bc404 ) {
        G4cout << " ----> Material " << "BC404" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4ThreeVector move_new;

    // Set visualization attributes
    G4VisAttributes* scintillator_vis_att = new G4VisAttributes(G4Colour(0.0,0.3,0.7));
    scintillator_vis_att->SetVisibility(true);
    G4VisAttributes* mylar_vis_att = new G4VisAttributes(G4Colour(0.6,0.6,0.6));
    mylar_vis_att->SetVisibility(true);

    //G4double extra = 2.0*this->mylar_thickness;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&&&&&&& Making the square-ish shaped scintillators 1-5 are &&&&&&&
    //&&&&&&& neg-y of origin, while 6-10 are in pos-y of origin &&&&&&&
    //&&&&&&& NOTE: this numbering scheme does not reflect the   &&&&&&&
    //&&&&&&&    actual way the detectors will be numbered       &&&&&&&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    //******* Making the Mylar coating

    //  G4Trd* square_mylar = this->squareMylar();

    G4SubtractionSolid* square_mylar = this->squareMylarWithCut();

    square_mylar_log = new G4LogicalVolume(square_mylar, material_mylar, "square_scintillator_log", 0, 0, 0);
    square_mylar_log->SetVisAttributes(mylar_vis_att);

    //******* Making the square scint and putting it in the Mylar

    G4Trd* square_scintillator = this->squareScintillator();

    square_scintillator_log = new G4LogicalVolume(square_scintillator, material_bc404, "sceptar_square_scintillator_log", 0, 0, 0);
    square_scintillator_log->SetVisAttributes(scintillator_vis_att);

    G4RotationMatrix* rotate_null = new G4RotationMatrix;
    G4ThreeVector move_null(0.0,0.0,0.0);

    this->assemblySquare->AddPlacedVolume(square_mylar_log, move_null, rotate_null);
    this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_null, rotate_null);




    ////  this->assembly->AddPlacedVolume(square_scintillator_log, move_null, rotate_null);

    ////  ******* Placing all the square scints, which involves placing the mylar

    //  G4ThreeVector move_square_scint1(0, -(this->square_scintillator_width +extra)/2.0, -this->square_scint_radial_distance);

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_square_scint1, rotate_null);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_square_scint1, rotate_null);

    //  G4RotationMatrix* rotate_square_scint2 = new G4RotationMatrix;
    //  rotate_square_scint2->rotateY(-1.0*(2.0*this->scint_angle1));

    //  G4ThreeVector move_square_scint2(this->square_scint_radial_distance*sin(2.0*this->scint_angle1), -(this->square_scintillator_width +extra)/2.0, -this->square_scint_radial_distance*cos(2.0*this->scint_angle1));

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_square_scint2, rotate_square_scint2);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_square_scint2, rotate_square_scint2);


    //  G4RotationMatrix* rotate_square_scint3 = new G4RotationMatrix;
    //  rotate_square_scint3->rotateY(-1.0*(M_PI -this->scint_angle1));

    //  G4ThreeVector move_square_scint3(this->square_scint_radial_distance*sin(this->scint_angle1), -(this->square_scintillator_width +extra)/2.0, this->square_scint_radial_distance*cos(this->scint_angle1));

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_square_scint3, rotate_square_scint3);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_square_scint3, rotate_square_scint3);

    //  G4RotationMatrix* rotate_square_scint4 = new G4RotationMatrix;
    //  rotate_square_scint4->rotateY(-1.0*(this->scint_angle1 -M_PI));

    //  G4ThreeVector move_square_scint4(-this->square_scint_radial_distance*sin(this->scint_angle1), -(this->square_scintillator_width +extra)/2.0, this->square_scint_radial_distance*cos(this->scint_angle1));

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_square_scint4, rotate_square_scint4);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_square_scint4, rotate_square_scint4);


    //  G4RotationMatrix* rotate_square_scint5 = new G4RotationMatrix;
    //  rotate_square_scint5->rotateY(-1.0*(-2.0*this->scint_angle1));

    //  G4ThreeVector move_square_scint5(-this->square_scint_radial_distance*sin(2.0*this->scint_angle1), -(this->square_scintillator_width +extra)/2.0, -this->square_scint_radial_distance*cos(2.0*this->scint_angle1));

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_square_scint5, rotate_square_scint5);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_square_scint5, rotate_square_scint5);

    //  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //  //&&&&&&& Since 6-10 are just 1-5, rotated through 180*deg, use          &&&&&&&
    //  //&&&&&&& neg-move_scint for corresponding detectors, with new rotation  &&&&&&&
    //  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    //  G4RotationMatrix* rotate_square_scint6 = new G4RotationMatrix;
    //  rotate_square_scint6->rotateY(M_PI);

    //  move_new = -1.0*move_square_scint1;

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_new, rotate_square_scint6);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_new, rotate_square_scint6);

    //  G4RotationMatrix* rotate_square_scint7 = new G4RotationMatrix;
    //  rotate_square_scint7->rotateY(-1.0*(2.0*this->scint_angle1 -M_PI));

    //  move_new = -1.0*move_square_scint2;

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_new, rotate_square_scint7);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_new, rotate_square_scint7);

    //  G4RotationMatrix* rotate_square_scint8 = new G4RotationMatrix;
    //  rotate_square_scint8->rotateY(-1.0*(-this->scint_angle1));

    //  move_new = -1.0*move_square_scint3;

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_new, rotate_square_scint8);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_new, rotate_square_scint8);

    //  G4RotationMatrix* rotate_square_scint9 = new G4RotationMatrix;
    //  rotate_square_scint9->rotateY(-1.0*(this->scint_angle1));

    //  move_new = -1.0*move_square_scint4;

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_new, rotate_square_scint9);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_new, rotate_square_scint9);

    //  G4RotationMatrix* rotate_square_scint10 = new G4RotationMatrix;
    //  rotate_square_scint10->rotateY(-1.0*(-2.0*this->scint_angle1 -M_PI));

    //  move_new = -1.0*move_square_scint5;

    //  this->assemblySquare->AddPlacedVolume(square_mylar_log, move_new, rotate_square_scint10);
    ////  this->assemblySquareSD->AddPlacedVolume(square_scintillator_log, move_new, rotate_square_scint10);


    //~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0
    //~o0o~~o0o~ Now we make the angled_scintillators, which are subtraction solids ~o0o~~o0o~
    //~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0

    //******* Making the Mylar coating

    //  G4SubtractionSolid* angled_mylar = this->angledMylar();

    G4SubtractionSolid* angled_mylar = this->angledMylarWithCut();

    angled_mylar_log = new G4LogicalVolume(angled_mylar, material_mylar, "angled_mylar_log", 0, 0, 0);
    angled_mylar_log->SetVisAttributes(mylar_vis_att);

    //******* Making the angled scint and putting it in the Mylar

    G4SubtractionSolid* angled_scintillator = this->angledScintillator();

    angled_scintillator_log = new G4LogicalVolume(angled_scintillator, material_bc404, "sceptar_angled_scintillator_log", 0, 0, 0);
    angled_scintillator_log->SetVisAttributes(scintillator_vis_att);

    this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_null, rotate_null);
    this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_null, rotate_null);






    ////  this->assembly->AddPlacedVolume(angled_mylar_log, move_null, rotate_null);

    //  //******* Placing all the angled scints, which involves placing the mylar

    //  G4RotationMatrix* rotate_angled_scint1 = new G4RotationMatrix;
    //  rotate_angled_scint1->rotateX(-1.0*(this->scint_angle4 -M_PI/2.0));

    //  G4ThreeVector move_angled_scint1(0, -this->angled_scint_move_back, -this->angled_scint_radial_distance);

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_angled_scint1, rotate_angled_scint1);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_angled_scint1, rotate_angled_scint1);

    //  G4RotationMatrix* rotate_angled_scint2 = new G4RotationMatrix;
    //  rotate_angled_scint2->rotateX(-1.0*(this->scint_angle4));
    //  rotate_angled_scint2->rotateZ(-1.0*(-2.0*this->scint_angle1));
    //  rotate_angled_scint2->rotateX(-1.0*(-M_PI/2.0));

    //  G4ThreeVector move_angled_scint2(this->angled_scint_radial_distance*sin(2.0*this->scint_angle1), -this->angled_scint_move_back, -this->angled_scint_radial_distance*cos(2.0*this->scint_angle1));

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_angled_scint2, rotate_angled_scint2);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_angled_scint2, rotate_angled_scint2);

    //  G4RotationMatrix* rotate_angled_scint3 = new G4RotationMatrix;
    //  rotate_angled_scint3->rotateX(-1.0*(this->scint_angle4));
    //  rotate_angled_scint3->rotateZ(-1.0*(this->scint_angle1 -M_PI));
    //  rotate_angled_scint3->rotateX(-1.0*(-M_PI/2.0));

    //  G4ThreeVector move_angled_scint3(this->angled_scint_radial_distance*sin(this->scint_angle1), -this->angled_scint_move_back, this->angled_scint_radial_distance*cos(this->scint_angle1));

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_angled_scint3, rotate_angled_scint3);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_angled_scint3, rotate_angled_scint3);

    //  G4RotationMatrix* rotate_angled_scint4 = new G4RotationMatrix;
    //  rotate_angled_scint4->rotateX(-1.0*(this->scint_angle4));
    //  rotate_angled_scint4->rotateZ(-1.0*(M_PI -this->scint_angle1));
    //  rotate_angled_scint4->rotateX(-1.0*(-M_PI/2.0));

    //  G4ThreeVector move_angled_scint4(-this->angled_scint_radial_distance*sin(this->scint_angle1), -this->angled_scint_move_back, this->angled_scint_radial_distance*cos(this->scint_angle1));

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_angled_scint4, rotate_angled_scint4);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_angled_scint4, rotate_angled_scint4);

    //  G4RotationMatrix* rotate_angled_scint5 = new G4RotationMatrix;
    //  rotate_angled_scint5->rotateX(-1.0*(this->scint_angle4));
    //  rotate_angled_scint5->rotateZ(-1.0*(2.0*this->scint_angle1));
    //  rotate_angled_scint5->rotateX(-1.0*(-M_PI/2.0));

    //  G4ThreeVector move_angled_scint5(-this->angled_scint_radial_distance*sin(2.0*this->scint_angle1), -this->angled_scint_move_back, -this->angled_scint_radial_distance*cos(2.0*this->scint_angle1));

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_angled_scint5, rotate_angled_scint5);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_angled_scint5, rotate_angled_scint5);

    //  //~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~
    //  //~o0o~~o0o~ Now we make the angled_scints 6-10, which are just 1-5 rotated throught   ~o0o~~o0o~
    //  //~o0o~~o0o~ 180*deg, so use neg-move from 1-5 to move 6-10, and apply proper rotation ~o0o~~o0o~
    //  //~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~~o0o~

    //  G4RotationMatrix* rotate_angled_scint6 = new G4RotationMatrix;
    //  rotate_angled_scint6->rotateX(-1.0*(this->scint_angle4 +M_PI/2.0));

    //  move_new = -1.0*move_angled_scint1;

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_new, rotate_angled_scint6);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_new, rotate_angled_scint6);

    //  G4RotationMatrix* rotate_angled_scint7 = new G4RotationMatrix;
    //  rotate_angled_scint7->rotateX(-1.0*(this->scint_angle4));
    //  rotate_angled_scint7->rotateZ(-1.0*(2.0*this->scint_angle1));
    //  rotate_angled_scint7->rotateX(-1.0*(M_PI/2.0));

    //  move_new = -1.0*move_angled_scint2;

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_new, rotate_angled_scint7);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_new, rotate_angled_scint7);

    //  G4RotationMatrix* rotate_angled_scint8 = new G4RotationMatrix;
    //  rotate_angled_scint8->rotateX(-1.0*(this->scint_angle4));
    //  rotate_angled_scint8->rotateZ(-1.0*(M_PI -this->scint_angle1));
    //  rotate_angled_scint8->rotateX(-1.0*(M_PI/2.0));

    //  move_new = -1.0*move_angled_scint3;

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_new, rotate_angled_scint8);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_new, rotate_angled_scint8);

    //  G4RotationMatrix* rotate_angled_scint9 = new G4RotationMatrix;
    //  rotate_angled_scint9->rotateX(-1.0*(this->scint_angle4));
    //  rotate_angled_scint9->rotateZ(-1.0*(this->scint_angle1 -M_PI));
    //  rotate_angled_scint9->rotateX(-1.0*(M_PI/2.0));

    //  move_new = -1.0*move_angled_scint4;

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_new, rotate_angled_scint9);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_new, rotate_angled_scint9);

    //  G4RotationMatrix* rotate_angled_scint10 = new G4RotationMatrix;
    //  rotate_angled_scint10->rotateX(-1.0*(this->scint_angle4));
    //  rotate_angled_scint10->rotateZ(-1.0*(-2.0*this->scint_angle1));
    //  rotate_angled_scint10->rotateX(-1.0*(M_PI/2.0));

    //  move_new = -1.0*move_angled_scint5;

    //  this->assemblyAngled->AddPlacedVolume(angled_mylar_log, move_new, rotate_angled_scint10);
    ////  this->assemblyAngledSD->AddPlacedVolume(angled_scintillator_log, move_new, rotate_angled_scint10);

    return 1;
}


//*****************************************************************
//ConstructDelrinShell builds a spherical shell of Delrin around
//the detectors, with two holes in it for the beam line
//*****************************************************************

G4int DetectionSystemSceptar::ConstructDelrinShell()
{
    G4Material* material_delrin = G4Material::GetMaterial("Delrin");
    if( !material_delrin ) {
        G4cout << " ----> Material " << "Delrin" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* Delrin_vis_att = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
    Delrin_vis_att->SetVisibility(true);

    G4SubtractionSolid* shell = this->DelrinShell();

    Delrin_shell_log = new G4LogicalVolume(shell, material_delrin, "Delrin_shell_log", 0, 0, 0);
    Delrin_shell_log->SetVisAttributes(Delrin_vis_att);

    G4RotationMatrix* rotate_null = new G4RotationMatrix;
    G4ThreeVector move_null(0.0,0.0,0.0);

    this->assembly->AddPlacedVolume(Delrin_shell_log, move_null, rotate_null);

    return 1;
} //end ::ConstructDelrinShell


//*******************************************************************
//Construct2ndDelrinShell builds a spherical shell of Delrin around
//the first sphere of Delrin, with two holes in it for the beam line.
//This is to simulate the Delrin in front of the germanium detectors
//*******************************************************************

G4int DetectionSystemSceptar::Construct2ndDelrinShell()
{
    G4Material* material_delrin = G4Material::GetMaterial("Delrin");
    if( !material_delrin ) {
        G4cout << " ----> Material " << "Delrin" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* Delrin2_vis_att = new G4VisAttributes(G4Colour(1,1,1));
    Delrin2_vis_att->SetVisibility(true);

    G4SubtractionSolid* shell2 = this->DelrinShell2();

    Delrin_shell2_log = new G4LogicalVolume(shell2, material_delrin, "Delrin_shell2_log", 0, 0, 0);
    Delrin_shell2_log->SetVisAttributes(Delrin2_vis_att);

    G4RotationMatrix* rotate_null = new G4RotationMatrix;
    G4ThreeVector move_null(0.0,0.0,0.0);

    this->assembly->AddPlacedVolume(Delrin_shell2_log, move_null, rotate_null);

    return 1;

} //end ::Construct2ndDelrinShell


//*******************************************************************
//ConstructHevimetShell builds a spherical shell of Hevimet around
//the two spheres of Delrin, with two holes in it for the beam line.
//This is to simulate the Hevimet in front of the germanium detectors
//*******************************************************************

G4int DetectionSystemSceptar::ConstructHevimetShell()
{
    G4Material* material_hevimet = G4Material::GetMaterial("Hevimet");
    if( !material_hevimet ) {
        G4cout << " ----> Material " << "Hevimet" << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    G4VisAttributes* Hevimet_vis_att = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
    Hevimet_vis_att->SetVisibility(true);

    G4SubtractionSolid* shell3 = this->HevimetShell();

    Hevimet_shell_log = new G4LogicalVolume(shell3, material_hevimet, "Hevimet_shell_log", 0, 0, 0);
    Hevimet_shell_log->SetVisAttributes(Hevimet_vis_att);

    G4RotationMatrix* rotate_null = new G4RotationMatrix;
    G4ThreeVector move_null(0.0,0.0,0.0);

    this->assembly->AddPlacedVolume(Hevimet_shell_log, move_null, rotate_null);

    return 1;

} //end ::Construct2ndDelrinShell

//G4int DetectionSystemSceptar::BuildCanVacuumVolume()
//{
//  G4Material* material = G4Material::GetMaterial(this->vacuum_material);
//  if( !material ) {
//    G4cout << " ----> Material " << this->vacuum_material << " not found, cannot build the detector shell! " << G4endl;
//    return 0;
//  }
//
//  // Set visualization attributes
//  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.1,0.0,0.9));
//  vis_att->SetVisibility(true);

//  G4ThreeVector direction =  G4ThreeVector(0,0,1);
//  G4double z_position;
//  G4ThreeVector move;
//  G4RotationMatrix* rotate = new G4RotationMatrix;

//  /////////////////////////////////////////////////////////////////////
//  // Build and Place Vacuum Can
//  /////////////////////////////////////////////////////////////////////
//  G4Tubs* can_vacuum_cylinder = BuildCanVacuum();

//  // Define rotation and movement objects
//  z_position	= ( (can_back_lid_thickness + crystal_dist_from_can_back) - (can_front_lid_thickness + crystal_dist_from_can_face) )/2.0;
//  move 		= z_position * direction;

//  // logical volume
//  if( can_vacuum_cylinder_log == NULL )
//  {
//    can_vacuum_cylinder_log = new G4LogicalVolume(can_vacuum_cylinder, material, "can_vacuum_cylinder_log", 0, 0, 0);
//    can_vacuum_cylinder_log->SetVisAttributes(vis_att);
//  }
//
//  // add physical can cylinder
//  this->assembly->AddPlacedVolume(can_vacuum_cylinder_log, move, rotate);
//
//  /////////////////////////////////////////////////////////////////////
//  // Build and Place Vacuum Front Lid
//  /////////////////////////////////////////////////////////////////////
//  G4Tubs* can_vacuum_front_lid = BuildCanVacuumFrontLid();

//  // logical volume
//  if( can_vacuum_front_lid_log == NULL )
//  {
//    can_vacuum_front_lid_log = new G4LogicalVolume(can_vacuum_front_lid, material, "can_vacuum_front_lid_log", 0, 0, 0);
//    can_vacuum_front_lid_log->SetVisAttributes(vis_att);
//  }

//  // place can vacuum front lid
//  z_position 	= (can_length_z/2.0) - (can_front_lid_thickness) - (crystal_dist_from_can_face/2.0);
//  move 		= z_position * direction;
//
//  // add physical BACK outer_can_lid
//  this->assembly->AddPlacedVolume(can_vacuum_front_lid_log, move, rotate);
//
//  /////////////////////////////////////////////////////////////////////
//  // Build and Place Vacuum Back Lid
//  /////////////////////////////////////////////////////////////////////
//  G4Tubs* can_vacuum_back_lid = BuildCanVacuumBackLid();

//  // logical volume
//  if( can_vacuum_back_lid_log == NULL )
//  {
//    can_vacuum_back_lid_log = new G4LogicalVolume(can_vacuum_back_lid, material, "can_vacuum_back_lid_log", 0, 0, 0);
//    can_vacuum_back_lid_log->SetVisAttributes(vis_att);
//  }

//  // place can vacuum front lid
//  z_position 	= - (can_length_z/2.0) + (can_back_lid_thickness) + (crystal_dist_from_can_back/2.0);
//  move 		= z_position * direction;
//
//  // add physical BACK outer_can_lid
//  this->assembly->AddPlacedVolume(can_vacuum_back_lid_log, move, rotate);
//
//  return 1;
//}//end ::BuildCanVacuumVolume

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Trd* DetectionSystemSceptar::squareMylar()
{
    G4double extra = 2.0*this->mylar_thickness;

    G4double half_length_longer_x = (this->square_scintillator_length +extra)/2.0;
    G4double half_length_shorter_x = half_length_longer_x - tan(this->scint_angle1)*(this->square_scintillator_thickness +extra);
    G4double half_length_y = (this->square_scintillator_width +extra)/2.0;
    G4double half_length_z = (this->square_scintillator_thickness +extra)/2.0;

    G4Trd* square_mylar = new G4Trd("square_mylar", half_length_longer_x, half_length_shorter_x,half_length_y, half_length_y, half_length_z);

    return square_mylar;
} //end ::squareMylar

G4SubtractionSolid* DetectionSystemSceptar::squareMylarWithCut()
{
    G4double extra = 2.0*this->mylar_thickness;

    G4double half_length_longer_x = (this->square_scintillator_length +extra)/2.0;
    G4double half_length_shorter_x = half_length_longer_x - tan(this->scint_angle1)*(this->square_scintillator_thickness +extra);
    G4double half_length_y = (this->square_scintillator_width +extra)/2.0;
    G4double half_length_z = (this->square_scintillator_thickness +extra)/2.0;

    G4Trd* square_mylar = new G4Trd("square_mylar", half_length_longer_x, half_length_shorter_x,half_length_y, half_length_y, half_length_z);

    G4Trd* scintCut = this->squareScintillator();

    G4RotationMatrix* rotateCut = new G4RotationMatrix;
    G4ThreeVector moveCut(0, 0, 0);

    G4SubtractionSolid* square_mylar_with_cut = new G4SubtractionSolid("square_mylar_with_cut", square_mylar, scintCut, rotateCut, moveCut);

    return square_mylar_with_cut;
} //end ::squareMylar

G4SubtractionSolid* DetectionSystemSceptar::angledMylar()
{
    G4double extra = 2.0*this->mylar_thickness;

    G4double half_length_neg_x = (this->angled_scintillator_long_width +extra)/2.0;
    G4double half_length_pos_x = (this->angled_scintillator_short_width +extra)/2.0;
    G4double half_length_z = (this->angled_scintillator_length +extra)/2.0;
    G4double half_length_y = (this->angled_scintillator_thickness +extra)/2.0;

    G4Trd* mylar = new G4Trd("mylar", half_length_neg_x, half_length_pos_x, half_length_y, half_length_y, half_length_z);

    G4double half_thickness = (this->angled_scintillator_thickness +extra)/(2.0*cos(this->scint_angle3));
    G4double half_length = this->angled_scintillator_length +extra;

    G4Box* box = new G4Box("box", half_thickness, half_thickness, half_length);

    G4RotationMatrix* rotate_box1 = new G4RotationMatrix;
    rotate_box1->rotateY(this->scint_angle2);
    rotate_box1->rotateZ(-this->scint_angle3);

    G4ThreeVector move_box1((this->angled_scintillator_short_width +extra +tan(this->scint_angle2)*(this->angled_scintillator_length +extra) -tan(this->scint_angle3)*(this->angled_scintillator_thickness +extra)/cos(this->scint_angle2))/2.0 +cos(this->scint_angle3)*half_thickness/cos(this->scint_angle2), half_thickness*sin(this->scint_angle3), 0);

    G4SubtractionSolid* angled_mylar1 = new G4SubtractionSolid("angled_mylar1", mylar, box, rotate_box1, move_box1);

    G4RotationMatrix* rotate_box2 = new G4RotationMatrix;
    rotate_box2->rotateY(-this->scint_angle2);
    rotate_box2->rotateZ(this->scint_angle3);

    G4ThreeVector move_box2(-(this->angled_scintillator_short_width +extra +tan(this->scint_angle2)*(this->angled_scintillator_length +extra) -tan(this->scint_angle3)*(this->angled_scintillator_thickness +extra)/cos(this->scint_angle2))/2.0 -cos(this->scint_angle3)*half_thickness/cos(this->scint_angle2), half_thickness*sin(this->scint_angle3), 0);

    G4SubtractionSolid* angled_mylar = new G4SubtractionSolid("angled_mylar", angled_mylar1, box, rotate_box2, move_box2);

    return angled_mylar;

} //end ::angledMylar


G4SubtractionSolid* DetectionSystemSceptar::angledMylarWithCut()
{
    G4double extra = 2.0*this->mylar_thickness;

    G4double half_length_neg_x = (this->angled_scintillator_long_width +extra)/2.0;
    G4double half_length_pos_x = (this->angled_scintillator_short_width +extra)/2.0;
    G4double half_length_z = (this->angled_scintillator_length +extra)/2.0;
    G4double half_length_y = (this->angled_scintillator_thickness +extra)/2.0;

    G4Trd* mylar = new G4Trd("mylar", half_length_neg_x, half_length_pos_x, half_length_y, half_length_y, half_length_z);

    G4double half_thickness = (this->angled_scintillator_thickness +extra)/(2.0*cos(this->scint_angle3));
    G4double half_length = this->angled_scintillator_length +extra;

    G4Box* box = new G4Box("box", half_thickness, half_thickness, half_length);

    G4RotationMatrix* rotate_box1 = new G4RotationMatrix;
    rotate_box1->rotateY(this->scint_angle2);
    rotate_box1->rotateZ(-this->scint_angle3);

    G4ThreeVector move_box1((this->angled_scintillator_short_width +extra +tan(this->scint_angle2)*(this->angled_scintillator_length +extra) -tan(this->scint_angle3)*(this->angled_scintillator_thickness +extra)/cos(this->scint_angle2))/2.0 +cos(this->scint_angle3)*half_thickness/cos(this->scint_angle2), half_thickness*sin(this->scint_angle3), 0);

    G4SubtractionSolid* angled_mylar1 = new G4SubtractionSolid("angled_mylar1", mylar, box, rotate_box1, move_box1);

    G4RotationMatrix* rotate_box2 = new G4RotationMatrix;
    rotate_box2->rotateY(-this->scint_angle2);
    rotate_box2->rotateZ(this->scint_angle3);

    G4ThreeVector move_box2(-(this->angled_scintillator_short_width +extra +tan(this->scint_angle2)*(this->angled_scintillator_length +extra) -tan(this->scint_angle3)*(this->angled_scintillator_thickness +extra)/cos(this->scint_angle2))/2.0 -cos(this->scint_angle3)*half_thickness/cos(this->scint_angle2), half_thickness*sin(this->scint_angle3), 0);

    G4SubtractionSolid* angled_mylar = new G4SubtractionSolid("angled_mylar", angled_mylar1, box, rotate_box2, move_box2);

    G4SubtractionSolid* scintCut = this->angledScintillator();

    G4RotationMatrix* rotateCut = new G4RotationMatrix;
    G4ThreeVector moveCut(0, 0, 0);

    G4SubtractionSolid* angled_mylar_with_cut = new G4SubtractionSolid("angled_mylar_with_cut", angled_mylar, scintCut, rotateCut, moveCut);

    return angled_mylar_with_cut;

} //end ::angledMylar


G4Trd* DetectionSystemSceptar::squareScintillator()
{
    G4double half_length_longer_x = this->square_scintillator_length/2.0;
    G4double half_length_shorter_x = half_length_longer_x - tan(this->scint_angle1)*this->square_scintillator_thickness;
    G4double half_length_y = this->square_scintillator_width/2.0;
    G4double half_length_z = this->square_scintillator_thickness/2.0;

    G4Trd* square_scintillator = new G4Trd("square_scintillator", half_length_longer_x, half_length_shorter_x,half_length_y, half_length_y, half_length_z);

    return square_scintillator;
} //end ::squareScintillator


G4SubtractionSolid* DetectionSystemSceptar::angledScintillator()
{
    G4double half_length_neg_x = this->angled_scintillator_long_width/2.0;
    G4double half_length_pos_x = this->angled_scintillator_short_width/2.0;
    G4double half_length_z = this->angled_scintillator_length/2.0;
    G4double half_length_y = this->angled_scintillator_thickness/2.0;

    G4Trd* scint = new G4Trd("scint", half_length_neg_x, half_length_pos_x, half_length_y, half_length_y, half_length_z);

    G4double half_thickness = this->angled_scintillator_thickness/(2.0*cos(this->scint_angle3));
    G4double half_length = this->angled_scintillator_length;

    G4Box* box = new G4Box("box", half_thickness, half_thickness, half_length);

    G4RotationMatrix* rotate_box1 = new G4RotationMatrix;
    rotate_box1->rotateY(this->scint_angle2);
    rotate_box1->rotateZ(-this->scint_angle3);

    G4ThreeVector move_box1((this->angled_scintillator_short_width +tan(this->scint_angle2)*this->angled_scintillator_length -tan(this->scint_angle3)*this->angled_scintillator_thickness/cos(this->scint_angle2))/2.0 +cos(this->scint_angle3)*half_thickness/cos(this->scint_angle2), half_thickness*sin(this->scint_angle3), 0);

    G4SubtractionSolid* angled_scint1 = new G4SubtractionSolid("angled_scint1", scint, box, rotate_box1, move_box1);

    G4RotationMatrix* rotate_box2 = new G4RotationMatrix;
    rotate_box2->rotateY(-this->scint_angle2);
    rotate_box2->rotateZ(this->scint_angle3);

    G4ThreeVector move_box2(-(this->angled_scintillator_short_width +tan(this->scint_angle2)*this->angled_scintillator_length -tan(this->scint_angle3)*this->angled_scintillator_thickness/cos(this->scint_angle2))/2.0 -cos(this->scint_angle3)*half_thickness/cos(this->scint_angle2), half_thickness*sin(this->scint_angle3), 0);

    G4SubtractionSolid* angled_scint = new G4SubtractionSolid("angled_scint", angled_scint1, box, rotate_box2, move_box2);

    return angled_scint;

} //end ::angledScintillator


G4SubtractionSolid* DetectionSystemSceptar::DelrinShell()
{
    G4Sphere* sphere = new G4Sphere("sphere", this->Delrin_inner_radius, this->Delrin_outer_radius, 0, 2.0*M_PI, 0, M_PI);

    G4Tubs* chop_tub = new G4Tubs("chop_tub", 0, this->Delrin_hole_radius, this->Delrin_outer_radius +1.0, 0, 2.0*M_PI);

    G4RotationMatrix* rotate_chop_tub = new G4RotationMatrix;
    rotate_chop_tub->rotateX(-M_PI/2.0);

    G4ThreeVector move_null(0.0,0.0,0.0);

    G4SubtractionSolid* Delrin_shell = new G4SubtractionSolid("Delrin_shell", sphere, chop_tub, rotate_chop_tub, move_null);

    return Delrin_shell;
} //end ::DelrinShell


G4SubtractionSolid* DetectionSystemSceptar::DelrinShell2()
{
    G4Sphere* sphere = new G4Sphere("sphere", this->Delrin2_inner_radius, this->Delrin2_outer_radius, 0, 2.0*M_PI, 0, M_PI);

    G4Tubs* chop_tub = new G4Tubs("chop_tub", 0, this->Delrin_hole_radius, this->Delrin2_outer_radius +1.0, 0, 2.0*M_PI);

    G4RotationMatrix* rotate_chop_tub = new G4RotationMatrix;
    rotate_chop_tub->rotateX(-M_PI/2.0);

    G4ThreeVector move_null(0.0,0.0,0.0);

    G4SubtractionSolid* Delrin_shell2 = new G4SubtractionSolid("Delrin_shell2", sphere, chop_tub, rotate_chop_tub, move_null);

    return Delrin_shell2;
} //end ::DelrinShell2


G4SubtractionSolid* DetectionSystemSceptar::HevimetShell()
{
    G4Sphere* sphere = new G4Sphere("sphere", this->Hevimet_inner_radius, this->Hevimet_outer_radius,0, 2.0*M_PI, 0, M_PI);

    G4Tubs* chop_tub = new G4Tubs("chop_tub", 0, this->Delrin_hole_radius, this->Hevimet_outer_radius +1.0,0, 2.0*M_PI);

    G4RotationMatrix* rotate_chop_tub = new G4RotationMatrix;
    rotate_chop_tub->rotateX(-M_PI/2.0);

    G4ThreeVector move_null(0.0,0.0,0.0);

    G4SubtractionSolid* Hevimet_shell = new G4SubtractionSolid("Hevimet_shell", sphere, chop_tub, rotate_chop_tub, move_null);

    return Hevimet_shell;
} //end ::DelrinShell2

