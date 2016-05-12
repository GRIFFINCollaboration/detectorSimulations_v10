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

#include "DetectionSystemDescant.hh"
#include "ApparatusDescantStructure.hh"

#include "G4SystemOfUnits.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4VSolid.hh"
#include "G4Polyhedra.hh"

#include <vector>


ApparatusDescantStructure::ApparatusDescantStructure() :
    descant_structure_log(0)
{

    // DESCANT structure properties
    structure_material 	       = "G4_Al";
    structure_inner_radius     = 915.82*mm;// using freeCAD to calculate thickness of the shell from STEP file
    structure_outer_radius     = 1007.9566*mm;
    start_Phi                  = 0.0*deg;
    delta_Phi                  = (360.0/5)*deg;
    start_Theta		           = 13.0*deg;
    delta_Theta                = 67.8743593*deg - start_Theta;
    structure_radial_distance  = 85.75*cm;

    // for red polyhedra subtraction
    numSides = 6;
    numZplanes = 2;
    cutExtra = 4.0*mm;
    polyRadiusRed = (171.57*mm)/2. + cutExtra; // taken from furthest distance between opposide sides on a red detector

    G4double rInnerRedArray[2] = {0.0*mm, 0.0*mm};
    memcpy(rInnerRed, rInnerRedArray, sizeof(rInnerRed));
    G4double rOuterRedArray[2] = {polyRadiusRed, polyRadiusRed};
    memcpy(rOuterRed, rOuterRedArray, sizeof(rOuterRed));
    G4double zPlanesArray[2] = {-500.0*mm, 500.0*mm};
    memcpy(zPlanes, zPlanesArray, sizeof(zPlanes));

    // for green/yellow polyhedra subtraction
    polyRadiusGreenYellow = (150.96*mm)/2. + cutExtra; // taken from furthest distance between opposide sides on a red detector

    G4double rInnerGreenYellowArray[2] = {0.0*mm, 0.0*mm};
    memcpy(rInnerGreenYellow, rInnerGreenYellowArray, sizeof(rInnerGreenYellow));
    G4double rOuterGreenYellowArray[2] = {polyRadiusGreenYellow, polyRadiusGreenYellow};
    memcpy(rOuterGreenYellow, rOuterGreenYellowArray, sizeof(rOuterGreenYellow));
    //G4double zPlanesArray[2] = {-500.0*mm, 500.0*mm}; // already done for red subtraction
    memcpy(zPlanes, zPlanesArray, sizeof(zPlanes));

    // for white polyhedra subtraction
    cutExtra = 10.0*mm;
    polyRadiusWhite = (150.86*mm)/2. + cutExtra;

    G4double rInnerWhiteArray[2] = {0.0*mm, 0.0*mm};
    memcpy(rInnerWhite, rInnerWhiteArray, sizeof(rInnerWhite));
    G4double rOuterWhiteArray[2] = {polyRadiusWhite, polyRadiusWhite};
    memcpy(rOuterWhite, rOuterWhiteArray, sizeof(rOuterWhite));
    //G4double zPlanesArray[2] = {-500.0*mm, 500.0*mm}; // already done for red subtraction
    memcpy(zPlanes, zPlanesArray, sizeof(zPlanes));

    // for polyhedra backing
    backRadius = 105.0*mm;

    G4double rInnerBackArray[2] = {0.0*mm, 0.0*mm};
    memcpy(rInnerBack, rInnerBackArray, sizeof(rInnerBack));
    G4double rOuterBackArray[2] = {backRadius, backRadius};
    memcpy(rOuterBack, rOuterBackArray, sizeof(rOuterBack));
    G4double zPlanesBackArray[2] = {-9.0*mm, 9.0*mm};
    memcpy(zPlanesBack, zPlanesBackArray, sizeof(zPlanesBack));

    // for cylinder subtraction from back plates
    rMin = 0.0*mm;
    rMax = 50.0*mm;
    halfLength = 2000.0*mm;
    sPhi = 0.0*deg;
    dPhi = 360.0*deg;

    // for subtraction box
    subBoxX = subBoxY = subBoxZ = 1000.0*mm;
    moveCut = 100.0*mm;

    // The Euler angles from James' MSc thesis which gives us the detector positions
    // Some of the angles for the green and yellow detectors are wrong in James' thesis,
    // note the +180 on a few angles.

    G4double blue_alpha_beta_gamma_array[15][3] = {
        {	-22.39*deg,	-33.76*deg,	-71.19*deg},
        {	-49.61*deg,	-33.76*deg,	71.19*deg},
        {	-36.00*deg,	-46.06*deg,	180.00*deg},
        {	85.61*deg,	33.76*deg,	108.81*deg},
        {	58.39*deg,	33.76*deg,	251.19*deg},
        {	72.00*deg,	46.06*deg,	360.00*deg},
        {	13.61*deg,	33.76*deg,	108.81*deg},
        {	-13.61*deg,	33.76*deg,	251.19*deg},
        {	0.00*deg,	46.06*deg,	360.00*deg},
        {	-58.39*deg,	33.76*deg,	108.81*deg},
        {	-85.61*deg,	33.76*deg,	251.19*deg},
        {	-72.00*deg,	46.06*deg,	360.00*deg},
        {	49.61*deg,	-33.76*deg,	-71.19*deg},
        {	22.39*deg,	-33.76*deg,	71.19*deg},
        {	36.00*deg,	-46.06*deg,	180.00*deg}
    };

    memcpy(blue_alpha_beta_gamma, blue_alpha_beta_gamma_array, sizeof(blue_alpha_beta_gamma));

    G4double green_alpha_beta_gamma_array[10][3] = {
        {	-12.45*deg,	-60.47*deg,	-12.45*deg},
        {	-27.73*deg,	-58.55*deg,	-4.54*deg},
        {	-84.45*deg,	-60.47*deg,	-12.45*deg},
        {	80.27*deg,	58.55*deg,	(-4.54+180)*deg},
        {	23.55*deg,	60.47*deg,	167.55*deg},
        {	8.27*deg,	58.55*deg,	175.46*deg},
        {	-48.45*deg,	60.47*deg,	167.55*deg},
        {	-63.73*deg,	58.55*deg,	175.46*deg},
        {	59.55*deg,	-60.47*deg,	-12.45*deg},
        {	44.27*deg,	-58.55*deg,	-4.54*deg}
    };

    memcpy(green_alpha_beta_gamma, green_alpha_beta_gamma_array, sizeof(green_alpha_beta_gamma));

    G4double red_alpha_beta_gamma_array[15][3] = {
        {	-36.00*deg,	-20.39*deg,	180.00*deg},
        {	-16.04*deg,	-47.84*deg,	45.10*deg},
        {	-55.96*deg,	-47.84*deg,	314.90*deg},
        {	72.00*deg,	20.39*deg,	0.00*deg},
        {	-88.04*deg,	-47.84*deg,	45.10*deg},
        {	52.04*deg,	47.84*deg,	134.90*deg},
        {	0.00*deg,	20.39*deg,	0.00*deg},
        {	19.96*deg,	47.84*deg,	225.10*deg},
        {	-19.96*deg,	47.84*deg,	134.90*deg},
        {	-72.00*deg,	20.39*deg,	0.00*deg},
        {	-52.04*deg,	47.84*deg,	225.10*deg},
        {	88.04*deg,	-47.84*deg,	314.90*deg},
        {	36.00*deg,	-20.39*deg,	180.00*deg},
        {	55.96*deg,	-47.84*deg,	45.10*deg},
        {	16.04*deg,	-47.84*deg,	314.90*deg}
    };

    memcpy(red_alpha_beta_gamma, red_alpha_beta_gamma_array, sizeof(red_alpha_beta_gamma));

    G4double white_alpha_beta_gamma_array[20][3] = {
        {	0.00*deg,	-11.37*deg,	90.00*deg},
        {	0.00*deg,	-24.67*deg,	90.00*deg},
        {	0.00*deg,	-38.76*deg,	-90.00*deg},
        {	0.00*deg,	-52.06*deg,	-90.00*deg},
        {	-72.00*deg,	-11.37*deg,	90.00*deg},
        {	-72.00*deg,	-24.67*deg,	90.00*deg},
        {	-72.00*deg,	-38.76*deg,	-90.00*deg},
        {	-72.00*deg,	-52.06*deg,	-90.00*deg},
        {	36.00*deg,	11.37*deg,	270.00*deg},
        {	36.00*deg,	24.67*deg,	270.00*deg},
        {	36.00*deg,	38.76*deg,	90.00*deg},
        {	36.00*deg,	52.06*deg,	90.00*deg},
        {	-36.00*deg,	11.37*deg,	270.00*deg},
        {	-36.00*deg,	24.67*deg,	270.00*deg},
        {	-36.00*deg,	38.76*deg,	90.00*deg},
        {	-36.00*deg,	52.06*deg,	90.00*deg},
        {	72.00*deg,	-11.37*deg,	90.00*deg},
        {	72.00*deg,	-24.67*deg,	90.00*deg},
        {	72.00*deg,	-38.76*deg,	-90.00*deg},
        {	72.00*deg,	-52.06*deg,	-90.00*deg} //59
    };

    memcpy(white_alpha_beta_gamma, white_alpha_beta_gamma_array, sizeof(white_alpha_beta_gamma));

    G4double yellow_alpha_beta_gamma_array[10][3] = {
        {	-44.27*deg,	-58.55*deg,	(4.54+180)*deg},
        {	-59.55*deg,	-60.47*deg,	192.45*deg},
        {	63.73*deg,	58.55*deg,	4.54*deg},
        {	48.45*deg,	60.47*deg,	(192.45+180)*deg},
        {	-8.27*deg,	58.55*deg,	(184.54+180)*deg},
        {	-23.55*deg,	60.47*deg,	372.45*deg},
        {	-80.27*deg,	58.55*deg,	(184.54+180)*deg},
        {	84.45*deg,	-60.47*deg,	(372.45+180)*deg},
        {	27.73*deg,	-58.55*deg,	(4.54+180)*deg},
        {	12.45*deg,	-60.47*deg,	192.45*deg}
    };

    memcpy(yellow_alpha_beta_gamma, yellow_alpha_beta_gamma_array, sizeof(yellow_alpha_beta_gamma));

    grey_colour     = G4Colour(0.9,0.9,0.9);
}


ApparatusDescantStructure::~ApparatusDescantStructure()
{
    // LogicalVolumes
    delete descant_structure_log;
}



G4int ApparatusDescantStructure::Build()
{
// Why is this needed?
    G4AssemblyVolume* myDescantStructureAssembly = new G4AssemblyVolume();
    this->assemblyDescantStructure = myDescantStructureAssembly;

    BuildDescantShell();

    return 1;
}



G4int ApparatusDescantStructure::PlaceDescantStructure(G4LogicalVolume* exp_hall_log)
{
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4ThreeVector move = G4ThreeVector(0.0, 0.0, 0.0);

    for(G4int i = 0; i < 5; i++){
        rotate = new G4RotationMatrix;
        rotate->rotateZ((i*72.0)*deg);
        assemblyDescantStructure->MakeImprint(exp_hall_log, move, rotate);
    }
        //assemblyDescantStructure->MakeImprint(exp_hall_log, move, rotate);

	return 0;
}


G4int ApparatusDescantStructure::BuildDescantShell()
{
    G4ThreeVector move = G4ThreeVector(0.0, 0.0, 0.0);
    G4RotationMatrix * rotate = new G4RotationMatrix;
    G4int detector_number = 70;
    G4int idx;
    structure_radial_distance = structure_outer_radius + 8.0*mm;

    // for cutting the detector shapes
    G4Box * subtractionBox = new G4Box("subtraction_box", subBoxX, subBoxY, subBoxZ);

    // for construction of main DESCANT shell structure
	G4Sphere * alum_sphere = new G4Sphere("descant sphere",
                                       structure_inner_radius, structure_outer_radius, start_Phi, delta_Phi, start_Theta, delta_Theta);

    // for the subtraction of each detector shape from the solid shell
        //red - cut
    G4Polyhedra * polyRed = new G4Polyhedra("poly_red", sPhi, dPhi, numSides, numZplanes, zPlanes, rInnerRed, rOuterRed);
    rotate = new G4RotationMatrix;
    rotate->rotateZ(20.0*deg);
    moveCut = 120.0*mm;
    move = G4ThreeVector(1.0*subBoxX + moveCut, 0.0, 0.0);
    //G4SubtractionSolid * polyRedCut1 = new G4SubtractionSolid("poly_red_cut_1", polyRed, subtractionBox, rotate, move);
    rotate = new G4RotationMatrix;
    rotate->rotateZ(-20.0*deg);
    moveCut = 120.0*mm;
    move = G4ThreeVector(1.0*subBoxX + moveCut, 0.0, 0.0);
    //G4SubtractionSolid * polyRedCut2 = new G4SubtractionSolid("poly_red_cut_2", polyRedCut1, subtractionBox, rotate, move);
        //white
    G4Polyhedra * polyWhite = new G4Polyhedra("poly_white", sPhi, dPhi, numSides, numZplanes, zPlanes, rInnerWhite, rOuterWhite);
        //green or yellow before cut
    G4Polyhedra * polyGreenYellow = new G4Polyhedra("poly_green_yellow", sPhi, dPhi, numSides, numZplanes, zPlanes, rInnerGreenYellow, rOuterGreenYellow);
        //green after cut
    rotate = new G4RotationMatrix;
    moveCut = 65.0*mm;
    move = G4ThreeVector(1.0*subBoxX + moveCut, 0.0, 0.0);
    G4SubtractionSolid * polyGreenCut = new G4SubtractionSolid("poly_green_cut", polyGreenYellow, subtractionBox, rotate, move);
        //yellow after cut
    rotate = new G4RotationMatrix;
    moveCut = 65.0*mm;
    move = G4ThreeVector(-1.0*subBoxX - moveCut, 0.0, 0.0);
    G4SubtractionSolid * polyYellowCut = new G4SubtractionSolid("poly_yellow_cut", polyGreenYellow, subtractionBox, rotate, move);

    // for construction of the back plates of each DESCANT detector
        // general
    //rotate = new G4RotationMatrix;
    //move = G4ThreeVector(0.0, 0.0, 0.0);
    //G4Polyhedra * shellBack = new G4Polyhedra("shell_back", sPhi, dPhi, numSides, numZplanes, zPlanesBack, rInnerBack, rOuterBack);
    //G4Tubs * cylinder = new G4Tubs("cylinder", rMin, rMax, halfLength, sPhi, dPhi);
    //G4SubtractionSolid * shellBackHole = new G4SubtractionSolid("shell_back_with_holes", shellBack, cylinder, rotate, move);
        // yellow
    //moveCut = 90.0*mm;
    //rotate = new G4RotationMatrix;
    //move = G4ThreeVector(-1.0*subBoxX - moveCut, 0.0, 0.0);
    //G4SubtractionSolid * shellBackYellow = new G4SubtractionSolid("shell_back_with_holes_yellow", shellBackHole, subtractionBox, rotate, move);
        // green
    //moveCut = 90.0*mm;
    //rotate = new G4RotationMatrix;
    //move = G4ThreeVector(1.0*subBoxX + moveCut, 0.0, 0.0);
    //G4SubtractionSolid * shellBackGreen = new G4SubtractionSolid("shell_back_with_holes_green", shellBackHole, subtractionBox, rotate, move);

    // Colour
    G4VisAttributes* grey_vis_att = new G4VisAttributes(grey_colour);
    grey_vis_att->SetVisibility(true);

	 // Material
	 G4Material* structure_g4material = G4Material::GetMaterial(this->structure_material);
	 if( !structure_g4material ) {
		G4cout << " ----> Material " << this->structure_material << " not found, cannot build! " << G4endl;
		return 0;
	 }
	 
	 
	 //G4LogicalVolume * logicBack = new G4LogicalVolume(shellBackHole, structure_g4material, "logic_back", 0, 0, 0);
	 //logicBack->SetVisAttributes(grey_colour);
	 
	 //G4LogicalVolume * logicBackYellow = new G4LogicalVolume(shellBackYellow, structure_g4material, "logic_back_yellow", 0, 0, 0);
	 //logicBackYellow->SetVisAttributes(grey_colour);
	 
    //G4LogicalVolume * logicBackGreen = new G4LogicalVolume(shellBackGreen, structure_g4material, "logic_back_green", 0, 0, 0);
    //logicBackYellow->SetVisAttributes(grey_colour);




    // for making the DESCANT shell
    G4SubtractionSolid * pre_subtraction;
    G4SubtractionSolid * post_subtraction;
    std::vector<G4SubtractionSolid*> descant_shell;

    	// subtract Blue Detector
    for(G4int i=0; i<(detector_number-55); i++){
        idx = i;
        //
        move = G4ThreeVector(0.0, 0.0, structure_radial_distance);
        move.rotateZ(blue_alpha_beta_gamma[idx][2]);
        move.rotateY(blue_alpha_beta_gamma[idx][1]);
        move.rotateZ(blue_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(blue_alpha_beta_gamma[idx][2]);
        rotate->rotateY(blue_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(blue_alpha_beta_gamma[idx][0]);
        //if(i==4||i==6||i==8){assemblyDescantStructure->AddPlacedVolume(logicBack, move, rotate);}
        rotate->invert();
        std::stringstream temp;
        temp << "blue_subtraction_" << i;
        G4String name = temp.str();
        if(i==4||i==5||i==6||i==8){
            if(i == 4){ post_subtraction = new G4SubtractionSolid(name, alum_sphere, polyWhite, rotate, move);
                        descant_shell.push_back(post_subtraction);}
            else{       pre_subtraction = descant_shell.back();
                        post_subtraction = new G4SubtractionSolid(name, pre_subtraction, polyWhite, rotate, move);
                        descant_shell.push_back(post_subtraction);}
        }
    }

	 // subtract Green Detector
    for(G4int i=15; i<(detector_number-45); i++){
        idx = i-15;
        move = G4ThreeVector(0.0, 0.0, structure_radial_distance);
        move.rotateZ(green_alpha_beta_gamma[idx][2]);
        move.rotateY(green_alpha_beta_gamma[idx][1]);
        move.rotateZ(green_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(green_alpha_beta_gamma[idx][2]);
        rotate->rotateY(green_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(green_alpha_beta_gamma[idx][0]);
        //if(i==19||i==20){assemblyDescantStructure->AddPlacedVolume(logicBackGreen, move, rotate);}
        rotate->invert();
        std::stringstream temp;
        temp << "green_subtraction_" << i;
        G4String name = temp.str();
        if(i==19||i==20){
        pre_subtraction = descant_shell.back();
        post_subtraction = new G4SubtractionSolid(name, pre_subtraction, polyGreenCut, rotate, move);
        descant_shell.push_back(post_subtraction);
        }
    }

	// subtract Red Detector
    for(G4int i=25; i<(detector_number-30); i++){
        idx = i-25;
        move = G4ThreeVector(0.0, 0.0, structure_radial_distance);
        move.rotateZ(red_alpha_beta_gamma[idx][2]);
        move.rotateY(red_alpha_beta_gamma[idx][1]);
        move.rotateZ(red_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(red_alpha_beta_gamma[idx][2]);
        rotate->rotateY(red_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(red_alpha_beta_gamma[idx][0]);
        //if(i==30||i==31||i==32){assemblyDescantStructure->AddPlacedVolume(logicBack, move, rotate);}
        rotate->invert();
        std::stringstream temp;
        temp << "red_subtraction_" << i;
        G4String name = temp.str();
        if(i==28||i==30||i==31||i==32){
        pre_subtraction = descant_shell.back();
        post_subtraction = new G4SubtractionSolid(name, pre_subtraction, polyRed, rotate, move);
        //post_subtraction = new G4SubtractionSolid(name, pre_subtraction, polyRedCut2, rotate, move);
        descant_shell.push_back(post_subtraction);
        }
    }

    // subtract White Detector
    for(G4int i=40; i<(detector_number-10); i++){
        idx = i-40;
        move = G4ThreeVector(0.0, 0.0, structure_radial_distance);
        move.rotateZ(white_alpha_beta_gamma[idx][2]);
        move.rotateY(white_alpha_beta_gamma[idx][1]);
        move.rotateZ(white_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(white_alpha_beta_gamma[idx][2]);
        rotate->rotateY(white_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(white_alpha_beta_gamma[idx][0]);
        //if(i==48||i==49||i==50||i==51){assemblyDescantStructure->AddPlacedVolume(logicBack, move, rotate);}
        rotate->invert();
        std::stringstream temp;
        temp << "white_subtraction_" << i;
        G4String name = temp.str();
        if(i==48||i==49||i==50||i==51){
        pre_subtraction = descant_shell.back();
        post_subtraction = new G4SubtractionSolid(name, pre_subtraction, polyWhite, rotate, move);
        descant_shell.push_back(post_subtraction);
        }
    }
    // subtract Yellow Detector
    for(G4int i=60; i<(detector_number); i++){
        idx = i-60;
        move = G4ThreeVector(0.0, 0.0, structure_radial_distance);
        move.rotateZ(yellow_alpha_beta_gamma[idx][2]);
        move.rotateY(yellow_alpha_beta_gamma[idx][1]);
        move.rotateZ(yellow_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(yellow_alpha_beta_gamma[idx][2]);
        rotate->rotateY(yellow_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(yellow_alpha_beta_gamma[idx][0]);
        //if(i==62||i==63){assemblyDescantStructure->AddPlacedVolume(logicBackYellow, move, rotate);}
        rotate->invert();
        std::stringstream temp;
        temp << "yellow_subtraction_" << i;
        G4String name = temp.str();
        if(i==62||i==63){
        pre_subtraction = descant_shell.back();
        post_subtraction = new G4SubtractionSolid(name, pre_subtraction, polyYellowCut, rotate, move);
        descant_shell.push_back(post_subtraction);
        }
    }

	// Logical Volume
	if( descant_structure_log == NULL ){
        descant_structure_log = new G4LogicalVolume(descant_shell.back(), structure_g4material, "descant_structure_log", 0, 0, 0);
		descant_structure_log->SetVisAttributes(grey_vis_att);
    	}

        rotate = new G4RotationMatrix;
        move = G4ThreeVector(0.0,0.0,0.0);
        this->assemblyDescantStructure->AddPlacedVolume(descant_structure_log, move, rotate);


////////////////////////////////////////////////////////////


        //G4LogicalVolume * red_cut = new G4LogicalVolume(polyRedCut2, structure_g4material, "red_test", 0, 0, 0);
        //red_cut->SetVisAttributes(grey_vis_att);
        //this->assemblyDescantStructure->AddPlacedVolume(red_cut, move, rotate);


/////////////////////////////////////////////////////////////////

    return 1;
}


