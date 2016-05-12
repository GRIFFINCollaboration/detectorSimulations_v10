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

#include "G4SystemOfUnits.hh"

//#include "G4String.hh"
#include <string>

DetectionSystemDescant::DetectionSystemDescant(G4bool leadShield) :
    // LogicalVolumes
    blue_volume_log(0),
    green_volume_log(0),
    red_volume_log(0),
    white_volume_log(0),
    yellow_volume_log(0),
    blue_scintillator_volume_log(0),
    green_scintillator_volume_log(0),
    red_scintillator_volume_log(0),
    white_scintillator_volume_log(0),
    yellow_scintillator_volume_log(0),
    blue_lead_volume_log(0),
    green_lead_volume_log(0),
    red_lead_volume_log(0),
    white_lead_volume_log(0),
    yellow_lead_volume_log(0)
{
    // can properties
    can_length              = 150.0*mm;
    can_thickness           = 1.5*mm;
    can_material            = "G4_Al";
    liquid_material         = "Deuterated Scintillator";
    lead_material           = "G4_Pb";
    lead_shield_thickness   = 6.35*mm;
    //includeLead             = true;
    includeLead = leadShield;

    // The face of the detector is 50 cm from the origin, this does NOT include the lead shield.
    radial_distance         = 50*cm;
	// the radial distance of 84 cm projects the detectors to approximately where they should be cut out of the partial sphere for the descant shell
    //radial_distance 	    = 75.75*cm;
    // Saint Gobain data files, 6 points around the front face of the can, and 6 on the back
    // from Dan Brennan

    G4double blue[12][3] = {
        { 40.60*mm,  61.10*mm,    0.00*mm},
        {-34.70*mm,  61.10*mm,    0.00*mm},
        {-75.10*mm,   0.00*mm,    0.00*mm},
        {-34.70*mm, -61.10*mm,    0.00*mm},
        { 40.60*mm, -61.10*mm,    0.00*mm},
        { 76.30*mm,   0.00*mm,    0.00*mm},
        { 52.80*mm,  79.40*mm, -150.00*mm},
        {-45.10*mm,  79.40*mm, -150.00*mm},
        {-97.60*mm,   0.00*mm, -150.00*mm},
        {-45.10*mm, -79.40*mm, -150.00*mm},
        { 52.80*mm, -79.40*mm, -150.00*mm},
        { 99.20*mm,   0.00*mm, -150.00*mm}
    };

    G4double green[12][3] = {
        { 31.90*mm,  61.20*mm,    0.00*mm},
        {-31.90*mm,  61.20*mm,    0.00*mm},
        {-71.50*mm,   0.00*mm,    0.00*mm},
        {-31.90*mm, -55.30*mm,    0.00*mm},
        { 31.90*mm, -55.30*mm,    0.00*mm},
        { 47.90*mm,  36.60*mm,    0.00*mm},
        { 41.50*mm,  79.60*mm, -150.00*mm},
        {-41.50*mm,  79.60*mm, -150.00*mm},
        {-93.00*mm,   0.00*mm, -150.00*mm},
        {-41.50*mm, -71.90*mm, -150.00*mm},
        { 41.50*mm, -71.90*mm, -150.00*mm},
        { 62.30*mm,  47.60*mm, -150.00*mm}
    };

    G4double red[12][3] = {
        { 28.80*mm,  67.00*mm,    0.00*mm},
        {-39.70*mm,  64.10*mm,    0.00*mm},
        {-78.30*mm,   0.00*mm,    0.00*mm},
        {-39.70*mm, -64.10*mm,    0.00*mm},
        { 28.80*mm, -67.00*mm,    0.00*mm},
        { 56.20*mm,   0.00*mm,    0.00*mm},
        { 37.40*mm,  87.10*mm, -150.00*mm},
        {-51.60*mm,  83.30*mm, -150.00*mm},
        {-101.80*mm,  0.00*mm, -150.00*mm},
        {-51.60*mm, -83.30*mm, -150.00*mm},
        { 37.40*mm, -87.10*mm, -150.00*mm},
        { 73.10*mm,   0.00*mm, -150.00*mm}
    };

    G4double white[12][3] = {
        { 31.90*mm,  61.20*mm,    0.00*mm},
        {-31.90*mm,  61.20*mm,    0.00*mm},
        {-71.50*mm,   0.00*mm,    0.00*mm},
        {-31.90*mm, -55.30*mm,    0.00*mm},
        { 31.90*mm, -55.30*mm,    0.00*mm},
        { 71.50*mm,   0.00*mm,    0.00*mm},
        { 41.50*mm,  79.60*mm, -150.00*mm},
        {-41.50*mm,  79.60*mm, -150.00*mm},
        {-93.00*mm,   0.00*mm, -150.00*mm},
        {-41.50*mm, -71.90*mm, -150.00*mm},
        { 41.50*mm, -71.90*mm, -150.00*mm},
        { 93.00*mm,   0.00*mm, -150.00*mm}
    };

    G4double yellow[12][3] = {
        { 31.90*mm,   61.20*mm,    0.00*mm},
        {-31.90*mm,   61.20*mm,    0.00*mm},
        {-47.90*mm,   36.60*mm,    0.00*mm},
        {-31.90*mm,  -55.30*mm,    0.00*mm},
        { 31.90*mm,  -55.30*mm,    0.00*mm},
        { 71.50*mm,    0.00*mm,    0.00*mm},
        { 41.50*mm,   79.60*mm, -150.00*mm},
        {-41.50*mm,   79.60*mm, -150.00*mm},
        {-62.30*mm,   47.60*mm, -150.00*mm},
        {-41.50*mm,  -71.90*mm, -150.00*mm},
        { 41.50*mm,  -71.90*mm, -150.00*mm},
        { 93.00*mm,    0.00*mm, -150.00*mm}
    };

/////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(blue_detector,   blue,   sizeof(blue_detector));
    memcpy(green_detector,  green,  sizeof(green_detector));
    memcpy(red_detector,    red,    sizeof(red_detector));
    memcpy(white_detector,  white,  sizeof(white_detector));
    memcpy(yellow_detector, yellow, sizeof(yellow_detector));

    // These are the angles of the cuts for each of the 6 sides of the detectors
    // These were worked out by hand.
    G4double blue_phi_array[6] = {
        (90.0*deg),
        (180.0*deg-(33.4731*deg)),
        -1.0*(180.0*deg-(33.4731*deg)),
        -1.0*(90.0*deg),
        -1.0*(30.2972*deg),
        (30.2972*deg)
    };

    G4double green_phi_array[6] = {
        (90.0*deg),
        (180.0*deg-(32.9052*deg)),
        -1.0*(180.0*deg-(35.6062*deg)),
        -1.0*(90.0*deg),
        -1.0*(9.8763*deg),
        (33.0402*deg)
    };

    G4double red_phi_array[6] = {
        (90.0*deg + 2.4242*deg),
        (180.0*deg-(31.0557*deg)),
        -1.0*(180.0*deg-(31.0557*deg)),
        -1.0*(90.0*deg + 2.4242*deg),
        -1.0*(22.2424*deg),
        (22.2424*deg)
    };

    G4double white_phi_array[6] = {
        (90.0*deg),
        (180.0*deg-(32.9052*deg)),
        -1.0*(180.0*deg-(35.6062*deg)),
        -1.0*(90.0*deg),
        -1.0*(35.6062*deg),
        (32.9052*deg)
    };

    G4double yellow_phi_array[6] = {
        (90.0*deg),
        (180.0*deg-(33.0402*deg)),
        -1.0*(180.0*deg-(9.8763*deg)),
        -1.0*(90.0*deg),
        -1.0*(35.6062*deg),
        1.0*(32.9052*deg)
    };

    memcpy(blue_phi,    blue_phi_array,     sizeof(blue_phi));
    memcpy(green_phi,   green_phi_array,    sizeof(green_phi));
    memcpy(red_phi,     red_phi_array,      sizeof(red_phi));
    memcpy(white_phi,   white_phi_array,   sizeof(white_phi));
    memcpy(yellow_phi,  yellow_phi_array,   sizeof(yellow_phi));

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
        {	72.00*deg,	-52.06*deg,	-90.00*deg}
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

    blue_colour     = G4Colour(0.0/255.0,0.0/255.0,255.0/255.0);
    green_colour    = G4Colour(0.0/255.0,255.0/255.0,0.0/255.0);
    red_colour      = G4Colour(255.0/255.0,0.0/255.0,0.0/255.0);
    white_colour    = G4Colour(255.0/255.0,255.0/255.0,255.0/255.0);
    yellow_colour   = G4Colour(255.0/255.0,255.0/255.0,0.0/255.0);
    liquid_colour   = G4Colour(0.0/255.0,255.0/255.0,225.0/255.0);


    // for lanthinum bromide locations

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


DetectionSystemDescant::~DetectionSystemDescant()
{
    // LogicalVolumes
    delete blue_volume_log;
    delete green_volume_log;
    delete red_volume_log;
    delete white_volume_log;
    delete yellow_volume_log;
    delete blue_scintillator_volume_log;
    delete green_scintillator_volume_log;
    delete red_scintillator_volume_log;
    delete white_scintillator_volume_log;
    delete yellow_scintillator_volume_log;
    delete blue_lead_volume_log;
    delete green_lead_volume_log;
    delete red_lead_volume_log;
    delete white_lead_volume_log;
    delete yellow_lead_volume_log;
}


G4int DetectionSystemDescant::Build()
{
    // Build assembly volumes
    // assemblyBlue contains non-sensitive volumes of the blue detector
    // assemblyBlueScintillator contains only sensitive volumes of the blue detector
    G4AssemblyVolume* myAssemblyBlue = new G4AssemblyVolume();
    this->assemblyBlue = myAssemblyBlue;
    G4AssemblyVolume* myassemblyBlueScintillator = new G4AssemblyVolume();
    this->assemblyBlueScintillator = myassemblyBlueScintillator;

    G4AssemblyVolume* myAssemblyGreen = new G4AssemblyVolume();
    this->assemblyGreen = myAssemblyGreen;
    G4AssemblyVolume* myassemblyGreenScintillator = new G4AssemblyVolume();
    this->assemblyGreenScintillator = myassemblyGreenScintillator;

    G4AssemblyVolume* myAssemblyRed = new G4AssemblyVolume();
    this->assemblyRed = myAssemblyRed;
    G4AssemblyVolume* myassemblyRedScintillator = new G4AssemblyVolume();
    this->assemblyRedScintillator = myassemblyRedScintillator;

    G4AssemblyVolume* myAssemblyWhite = new G4AssemblyVolume();
    this->assemblyWhite = myAssemblyWhite;
    G4AssemblyVolume* myassemblyWhiteScintillator = new G4AssemblyVolume();
    this->assemblyWhiteScintillator = myassemblyWhiteScintillator;

    G4AssemblyVolume* myAssemblyYellow = new G4AssemblyVolume();
    this->assemblyYellow = myAssemblyYellow;
    G4AssemblyVolume* myassemblyYellowScintillator = new G4AssemblyVolume();
    this->assemblyYellowScintillator = myassemblyYellowScintillator;

    G4cout << "Building DESCANT..." << G4endl;

    G4cout << "BuildCanVolume" << G4endl;
    BuildCanVolume();

    G4cout << "BuildDetectorVolume" << G4endl;
    BuildDetectorVolume();

    return 1;
}

G4int DetectionSystemDescant::PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector_number)
{
    G4double position = radial_distance;

    G4int idx;
    G4double x = 0;
    G4double y = 0;
    G4double z = position;

    G4ThreeVector move;
    G4RotationMatrix* rotate;

    // Place Blue Detector
    for(G4int i=0; i<(detector_number-55); i++){
        move = G4ThreeVector(x,y,z);
        idx = i;
        move.rotateZ(blue_alpha_beta_gamma[idx][2]);
        move.rotateY(blue_alpha_beta_gamma[idx][1]);
        move.rotateZ(blue_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(blue_alpha_beta_gamma[idx][2]);
        rotate->rotateY(blue_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(blue_alpha_beta_gamma[idx][0]);
        assemblyBlue->MakeImprint(exp_hall_log, move, rotate, i);
        assemblyBlueScintillator->MakeImprint(exp_hall_log, move, rotate, i);

    }
    // Place Green Detector
    for(G4int i=15; i<(detector_number-45); i++){
        idx = i-15;
        move = G4ThreeVector(x,y,z);
        move.rotateZ(green_alpha_beta_gamma[idx][2]);
        move.rotateY(green_alpha_beta_gamma[idx][1]);
        move.rotateZ(green_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(green_alpha_beta_gamma[idx][2]);
        rotate->rotateY(green_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(green_alpha_beta_gamma[idx][0]);
        assemblyGreen->MakeImprint(exp_hall_log, move, rotate, i);
        assemblyGreenScintillator->MakeImprint(exp_hall_log, move, rotate, i);
    }
    // Place Red Detector
    for(G4int i=25; i<(detector_number-30); i++){
        idx = i-25;
        move = G4ThreeVector(x,y,z);
        move.rotateZ(red_alpha_beta_gamma[idx][2]);
        move.rotateY(red_alpha_beta_gamma[idx][1]);
        move.rotateZ(red_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(red_alpha_beta_gamma[idx][2]);
        rotate->rotateY(red_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(red_alpha_beta_gamma[idx][0]);
        assemblyRed->MakeImprint(exp_hall_log, move, rotate, i);
        assemblyRedScintillator->MakeImprint(exp_hall_log, move, rotate, i);
    }
    // Place White Detector
    for(G4int i=40; i<(detector_number-10); i++){
        idx = i-40;
        move = G4ThreeVector(x,y,z);
        move.rotateZ(white_alpha_beta_gamma[idx][2]);
        move.rotateY(white_alpha_beta_gamma[idx][1]);
        move.rotateZ(white_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(white_alpha_beta_gamma[idx][2]);
        rotate->rotateY(white_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(white_alpha_beta_gamma[idx][0]);
        assemblyWhite->MakeImprint(exp_hall_log, move, rotate, i);
        assemblyWhiteScintillator->MakeImprint(exp_hall_log, move, rotate, i);
    }
    // Place Yellow Detector
    for(G4int i=60; i<(detector_number); i++){
        idx = i-60;
        move = G4ThreeVector(x,y,z);
        move.rotateZ(yellow_alpha_beta_gamma[idx][2]);
        move.rotateY(yellow_alpha_beta_gamma[idx][1]);
        move.rotateZ(yellow_alpha_beta_gamma[idx][0]);
        rotate = new G4RotationMatrix;
        rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
        rotate->rotateZ(yellow_alpha_beta_gamma[idx][2]);
        rotate->rotateY(yellow_alpha_beta_gamma[idx][1]);
        rotate->rotateZ(yellow_alpha_beta_gamma[idx][0]);
        assemblyYellow->MakeImprint(exp_hall_log, move, rotate, i);
        assemblyYellowScintillator->MakeImprint(exp_hall_log, move, rotate, i);
    }

    return 1;
}

G4int DetectionSystemDescant::PlaceDetector(G4LogicalVolume* exp_hall_log, G4String color, G4ThreeVector pos, G4ThreeVector rot)
{
    G4RotationMatrix* rotate;

	 G4double alpha = rot.x();
	 G4double beta  = rot.y();
	 G4double gamma = rot.z();

    // Place Detector
	 rotate = new G4RotationMatrix;
	 rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
	 rotate->rotateZ(gamma);
	 rotate->rotateY(beta);
	 rotate->rotateZ(alpha);

	 if(color.compareTo("blue", G4String::ignoreCase) == 0) {
		assemblyBlue->MakeImprint(exp_hall_log, pos, rotate, 1);
		assemblyBlueScintillator->MakeImprint(exp_hall_log, pos, rotate, 1);
	 } else if(color.compareTo("red", G4String::ignoreCase) == 0) {
		assemblyRed->MakeImprint(exp_hall_log, pos, rotate, 1);
		assemblyRedScintillator->MakeImprint(exp_hall_log, pos, rotate, 1);
	 } else if(color.compareTo("white", G4String::ignoreCase) == 0) {
		assemblyWhite->MakeImprint(exp_hall_log, pos, rotate, 1);
		assemblyWhiteScintillator->MakeImprint(exp_hall_log, pos, rotate, 1);
	 } else if(color.compareTo("green", G4String::ignoreCase) == 0) {
		assemblyGreen->MakeImprint(exp_hall_log, pos, rotate, 1);
		assemblyGreenScintillator->MakeImprint(exp_hall_log, pos, rotate, 1);
	 } else if(color.compareTo("yellow", G4String::ignoreCase) == 0) {
		assemblyYellow->MakeImprint(exp_hall_log, pos, rotate, 1);
		assemblyYellowScintillator->MakeImprint(exp_hall_log, pos, rotate, 1);
	 }

    return 1;
}

G4int DetectionSystemDescant::PlaceDetectorAuxPorts(G4LogicalVolume* exp_hall_log, G4int detector_number, G4double radialpos)
{
    G4int detector_copy_ID = 0;

    G4cout << "Descant Detector Number = " << detector_number << G4endl;

    G4int copy_number = detector_copy_ID + detector_number;

    //G4double position = radialpos + can_front_lid_thickness +( this->can_length_z / 2.0 ) ;
    G4double position = radialpos + can_thickness_front +( this->can_length / 2.0 ) ;
    set_radial_pos = radialpos;

    G4double theta  = this->detectorAngles[detector_number][0];
    G4double phi    = this->detectorAngles[detector_number][1];
    //G4double alpha  = this->detectorAngles[detector_number][2]; // yaw
    G4double beta   = this->detectorAngles[detector_number][3]; // pitch
    G4double gamma  = this->detectorAngles[detector_number][4]; // roll

    G4double x = 0;
    G4double y = 0;
    G4double z = position;

    G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
    //rotate->rotateY(M_PI);
    rotate->rotateZ(90.0*deg);
    rotate->rotateY(M_PI+beta);
    rotate->rotateZ(gamma);

    G4ThreeVector move(transX(x,y,z,theta,phi), transY(x,y,z,theta,phi), transZ(x,y,z,theta,phi));

    assemblyWhite->MakeImprint(exp_hall_log, move, rotate, copy_number);
    assemblyWhiteScintillator->MakeImprint(exp_hall_log, move, rotate, copy_number);

    return 1;



}


G4int DetectionSystemDescant::BuildCanVolume()
{
    G4ThreeVector move_cut, move, direction;
    G4RotationMatrix* rotate_cut;
    G4RotationMatrix* rotate;
    G4double z_position;
    G4double extra_cut_length = 100.0*mm;

    G4Material* can_g4material = G4Material::GetMaterial(this->can_material);
    if( !can_g4material ) {
        return 0;
    }

    G4Material* lead_g4material = G4Material::GetMaterial(this->lead_material);
    if( !lead_g4material ) {
        G4cout << " ----> Material " << this->lead_material << " not found, cannot build! " << G4endl;
        return 0;
    }

    // BLUE
    // Set visualization attributes
    G4VisAttributes* blue_vis_att = new G4VisAttributes(blue_colour);
    blue_vis_att->SetVisibility(true);

    G4SubtractionSolid* blue_volume_1      = CanVolume(false, this->can_length, blue_detector, blue_phi);
    G4SubtractionSolid* blue_volume_1_cut  = CanVolume(true, this->can_length+extra_cut_length, blue_detector, blue_phi);

    move_cut = G4ThreeVector(0.0*mm,0.0*mm,-1.0*(this->can_thickness + (extra_cut_length/2.0)));
    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* blue_volume_2 = new G4SubtractionSolid("blue_volume_2", blue_volume_1, blue_volume_1_cut, rotate_cut, move_cut);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( blue_volume_log == NULL )
    {
        blue_volume_log = new G4LogicalVolume(blue_volume_2, can_g4material, "descant_blue_volume_log", 0, 0, 0);
        blue_volume_log->SetVisAttributes(blue_vis_att);
    }

    this->assemblyBlue->AddPlacedVolume(blue_volume_log, move, rotate);

    // BLUE LEAD
    // Set visualization attributes
    G4VisAttributes* blue_lead_vis_att = new G4VisAttributes(blue_colour);
    blue_lead_vis_att->SetVisibility(true);

    G4SubtractionSolid* blue_lead_volume_1 = CanVolume(false, this->lead_shield_thickness, blue_detector, blue_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = this->lead_shield_thickness/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( blue_lead_volume_log == NULL )
    {
        blue_lead_volume_log = new G4LogicalVolume(blue_lead_volume_1, lead_g4material, "descant_blue_lead_volume_log", 0, 0, 0);
        blue_lead_volume_log->SetVisAttributes(blue_lead_vis_att);
    }

    if(includeLead)
        this->assemblyBlue->AddPlacedVolume(blue_lead_volume_log, move, rotate);

    // GREEN
    // Set visualization attributes
    G4VisAttributes* green_vis_att = new G4VisAttributes(green_colour);
    green_vis_att->SetVisibility(true);

    G4SubtractionSolid* green_volume_1      = CanVolume(false, this->can_length, green_detector, green_phi);
    G4SubtractionSolid* green_volume_1_cut  = CanVolume(true, this->can_length+extra_cut_length, green_detector, green_phi);

    move_cut = G4ThreeVector(0.0*mm,0.0*mm,-1.0*(this->can_thickness + (extra_cut_length/2.0)));
    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* green_volume_2 = new G4SubtractionSolid("green_volume_2", green_volume_1, green_volume_1_cut, rotate_cut, move_cut);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( green_volume_log == NULL )
    {
        green_volume_log = new G4LogicalVolume(green_volume_2, can_g4material, "descant_green_volume_log", 0, 0, 0);
        green_volume_log->SetVisAttributes(green_vis_att);
    }

    this->assemblyGreen->AddPlacedVolume(green_volume_log, move, rotate);

    // GREEN LEAD
    // Set visualization attributes
    G4VisAttributes* green_lead_vis_att = new G4VisAttributes(green_colour);
    green_lead_vis_att->SetVisibility(true);

    G4SubtractionSolid* green_lead_volume_1 = CanVolume(false, this->lead_shield_thickness, green_detector, green_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = this->lead_shield_thickness/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( green_lead_volume_log == NULL )
    {
        green_lead_volume_log = new G4LogicalVolume(green_lead_volume_1, lead_g4material, "descant_green_lead_volume_log", 0, 0, 0);
        green_lead_volume_log->SetVisAttributes(green_lead_vis_att);
    }

    if(includeLead)
        this->assemblyGreen->AddPlacedVolume(green_lead_volume_log, move, rotate);

    // RED
    // Set visualization attributes
    G4VisAttributes* red_vis_att = new G4VisAttributes(red_colour);
    red_vis_att->SetVisibility(true);

    G4SubtractionSolid* red_volume_1      = CanVolume(false, this->can_length, red_detector, red_phi);
    G4SubtractionSolid* red_volume_1_cut  = CanVolume(true, this->can_length+extra_cut_length, red_detector, red_phi);

    move_cut = G4ThreeVector(0.0*mm,0.0*mm,-1.0*(this->can_thickness + (extra_cut_length/2.0)));
    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* red_volume_2 = new G4SubtractionSolid("red_volume_2", red_volume_1, red_volume_1_cut, rotate_cut, move_cut);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( red_volume_log == NULL )
    {
        red_volume_log = new G4LogicalVolume(red_volume_2, can_g4material, "descant_red_volume_log", 0, 0, 0);
        red_volume_log->SetVisAttributes(red_vis_att);
    }

    this->assemblyRed->AddPlacedVolume(red_volume_log, move, rotate);

    // RED LEAD
    // Set visualization attributes
    G4VisAttributes* red_lead_vis_att = new G4VisAttributes(red_colour);
    red_lead_vis_att->SetVisibility(true);

    G4SubtractionSolid* red_lead_volume_1 = CanVolume(false, this->lead_shield_thickness, red_detector, red_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = this->lead_shield_thickness/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( red_lead_volume_log == NULL )
    {
        red_lead_volume_log = new G4LogicalVolume(red_lead_volume_1, lead_g4material, "descant_red_lead_volume_log", 0, 0, 0);
        red_lead_volume_log->SetVisAttributes(red_lead_vis_att);
    }

    if(includeLead)
        this->assemblyRed->AddPlacedVolume(red_lead_volume_log, move, rotate);

    // WHITE
    // Set visualization attributes
    G4VisAttributes* white_vis_att = new G4VisAttributes(white_colour);
    //G4VisAttributes* white_vis_att = new G4VisAttributes(red_colour);
    white_vis_att->SetVisibility(true);

    G4SubtractionSolid* white_volume_1      = CanVolume(false, this->can_length, white_detector, white_phi);
    G4SubtractionSolid* white_volume_1_cut  = CanVolume(true, this->can_length+extra_cut_length, white_detector, white_phi);

    move_cut = G4ThreeVector(0.0*mm,0.0*mm,-1.0*(this->can_thickness + (extra_cut_length/2.0)));
    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* white_volume_2 = new G4SubtractionSolid("white_volume_2", white_volume_1, white_volume_1_cut, rotate_cut, move_cut);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( white_volume_log == NULL )
    {
        white_volume_log = new G4LogicalVolume(white_volume_2, can_g4material, "descant_white_volume_log", 0, 0, 0);
        white_volume_log->SetVisAttributes(white_vis_att);
    }

    this->assemblyWhite->AddPlacedVolume(white_volume_log, move, rotate);

    // WHITE LEAD
    // Set visualization attributes
    G4VisAttributes* white_lead_vis_att = new G4VisAttributes(white_colour);
    //G4VisAttributes* white_lead_vis_att = new G4VisAttributes(red_colour);
    white_lead_vis_att->SetVisibility(true);

    G4SubtractionSolid* white_lead_volume_1 = CanVolume(false, this->lead_shield_thickness, white_detector, white_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = this->lead_shield_thickness/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( white_lead_volume_log == NULL )
    {
        white_lead_volume_log = new G4LogicalVolume(white_lead_volume_1, lead_g4material, "descant_white_lead_volume_log", 0, 0, 0);
        white_lead_volume_log->SetVisAttributes(white_lead_vis_att);
    }

    if(includeLead)
        this->assemblyWhite->AddPlacedVolume(white_lead_volume_log, move, rotate);

    // YELLOW
    // Set visualization attributes
    G4VisAttributes* yellow_vis_att = new G4VisAttributes(yellow_colour);
    yellow_vis_att->SetVisibility(true);

    G4SubtractionSolid* yellow_volume_1      = CanVolume(false, this->can_length, yellow_detector, yellow_phi);
    G4SubtractionSolid* yellow_volume_1_cut  = CanVolume(true, this->can_length+extra_cut_length, yellow_detector, yellow_phi);

    move_cut = G4ThreeVector(0.0*mm,0.0*mm,-1.0*(this->can_thickness + (extra_cut_length/2.0)));
    rotate_cut = new G4RotationMatrix;

    G4SubtractionSolid* yellow_volume_2 = new G4SubtractionSolid("yellow_volume_2", yellow_volume_1, yellow_volume_1_cut, rotate_cut, move_cut);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( yellow_volume_log == NULL )
    {
        yellow_volume_log = new G4LogicalVolume(yellow_volume_2, can_g4material, "descant_yellow_volume_log", 0, 0, 0);
        yellow_volume_log->SetVisAttributes(yellow_vis_att);
    }

    this->assemblyYellow->AddPlacedVolume(yellow_volume_log, move, rotate);

    // YELLOW LEAD
    // Set visualization attributes
    G4VisAttributes* yellow_lead_vis_att = new G4VisAttributes(yellow_colour);
    yellow_lead_vis_att->SetVisibility(true);

    G4SubtractionSolid* yellow_lead_volume_1 = CanVolume(false, this->lead_shield_thickness, yellow_detector, yellow_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = this->lead_shield_thickness/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( yellow_lead_volume_log == NULL )
    {
        yellow_lead_volume_log = new G4LogicalVolume(yellow_lead_volume_1, lead_g4material, "descant_yellow_lead_volume_log", 0, 0, 0);
        yellow_lead_volume_log->SetVisAttributes(yellow_lead_vis_att);
    }

    if(includeLead)
        this->assemblyYellow->AddPlacedVolume(yellow_lead_volume_log, move, rotate);

    return 1;
}

G4int DetectionSystemDescant::BuildDetectorVolume()
{
    G4ThreeVector move, direction;
    G4RotationMatrix* rotate;
    G4double z_position;

    G4Material* material = G4Material::GetMaterial(this->liquid_material);
    if( !material ) {
        G4cout << " ----> Material " << this->liquid_material << " not found, cannot build the detector shell! " << G4endl;
        return 0;
    }

    // BLUE
    // Set visualization attributes
    G4VisAttributes* blue_vis_att = new G4VisAttributes(blue_colour);
    blue_vis_att->SetVisibility(true);

    G4SubtractionSolid* blue_volume_1_cut  = CanVolume(true, this->can_length-this->can_thickness, blue_detector, blue_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length+this->can_thickness)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( blue_scintillator_volume_log == NULL )
    {
        blue_scintillator_volume_log = new G4LogicalVolume(blue_volume_1_cut, material, "blue_scintillator_volume_log", 0, 0, 0);
        blue_scintillator_volume_log->SetVisAttributes(blue_vis_att);
    }

    this->assemblyBlueScintillator->AddPlacedVolume(blue_scintillator_volume_log, move, rotate);

    //GREEN
    // Set visualization attributes
    G4VisAttributes* green_vis_att = new G4VisAttributes(green_colour);
    green_vis_att->SetVisibility(true);

    G4SubtractionSolid* green_volume_1_cut  = CanVolume(true, this->can_length-this->can_thickness, green_detector, green_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length+this->can_thickness)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( green_scintillator_volume_log == NULL )
    {
        green_scintillator_volume_log = new G4LogicalVolume(green_volume_1_cut, material, "green_scintillator_volume_log", 0, 0, 0);
        green_scintillator_volume_log->SetVisAttributes(green_vis_att);
    }

    this->assemblyGreenScintillator->AddPlacedVolume(green_scintillator_volume_log, move, rotate);

    // RED
    // Set visualization attributes
    G4VisAttributes* red_vis_att = new G4VisAttributes(red_colour);
    red_vis_att->SetVisibility(true);

    G4SubtractionSolid* red_volume_1_cut  = CanVolume(true, this->can_length-this->can_thickness, red_detector, red_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length+this->can_thickness)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( red_scintillator_volume_log == NULL )
    {
        red_scintillator_volume_log = new G4LogicalVolume(red_volume_1_cut, material, "red_scintillator_volume_log", 0, 0, 0);
        red_scintillator_volume_log->SetVisAttributes(red_vis_att);
    }

    this->assemblyRedScintillator->AddPlacedVolume(red_scintillator_volume_log, move, rotate);

    // WHITE
    // Set visualization attributes
    G4VisAttributes* white_vis_att = new G4VisAttributes(white_colour);
    white_vis_att->SetVisibility(true);

    G4SubtractionSolid* white_volume_1_cut  = CanVolume(true, this->can_length-this->can_thickness, white_detector, white_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length+this->can_thickness)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( white_scintillator_volume_log == NULL )
    {
        white_scintillator_volume_log = new G4LogicalVolume(white_volume_1_cut, material, "white_scintillator_volume_log", 0, 0, 0);
        white_scintillator_volume_log->SetVisAttributes(white_vis_att);
    }

    this->assemblyWhiteScintillator->AddPlacedVolume(white_scintillator_volume_log, move, rotate);

    // YELLOW
    // Set visualization attributes
    G4VisAttributes* yellow_vis_att = new G4VisAttributes(yellow_colour);
    yellow_vis_att->SetVisibility(true);

    G4SubtractionSolid* yellow_volume_1_cut  = CanVolume(true, this->can_length-this->can_thickness, yellow_detector, yellow_phi);

    // Define rotation and movement objects
    direction 	= G4ThreeVector(0,0,1);
    z_position    = -1.0*(this->can_length+this->can_thickness)/2.0 ;
    move          = z_position * direction;

    rotate = new G4RotationMatrix;

    //logical volume
    if( yellow_scintillator_volume_log == NULL )
    {
        yellow_scintillator_volume_log = new G4LogicalVolume(yellow_volume_1_cut, material, "yellow_scintillator_volume_log", 0, 0, 0);
        yellow_scintillator_volume_log->SetVisAttributes(yellow_vis_att);
    }

    this->assemblyYellowScintillator->AddPlacedVolume(yellow_scintillator_volume_log, move, rotate);

    return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
// This function constructs a generic 6 sided shape defined by 12 input points, and the cut phi angles.
G4SubtractionSolid* DetectionSystemDescant::CanVolume(G4bool insideVol, G4double volume_length, G4double detector[12][3], G4double detector_phi[6])
{
    G4int idx1;
    G4int idx2;

    G4ThreeVector front_p1;
    G4ThreeVector front_p2;
    G4ThreeVector back_p1;
    G4ThreeVector back_p2;

    G4double start_phi      = 0.0*deg;
    G4double end_phi        = 360.0*deg;
    G4double inner_radius;
    G4double outer_radius;
    G4double half_length_z;

    inner_radius = 0.0*mm;
    outer_radius = 200.0*mm;
    half_length_z = volume_length/2.0;

    G4Tubs* can_volume = new G4Tubs("can_volume", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

    // Cut out a small cube out of the front of the detector to get a G4SubtractionSolid,
    // needed for CutVolumeOnFourPoints methods
    G4double cut_half_length_x = 5.0*mm;
    G4double cut_half_length_y = 5.0*mm;
    G4double cut_half_length_z = 5.0*mm;

    G4Box* cut_plate = new G4Box("cut_plate", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    G4ThreeVector move_cut = G4ThreeVector(0,outer_radius,1.0*volume_length/2.0);
    G4RotationMatrix* rotate_cut = new G4RotationMatrix;

    // snip snip
    G4SubtractionSolid* can_volume_cut_1 = new G4SubtractionSolid("can_volume_cut_1", can_volume, cut_plate, rotate_cut, move_cut);

    // first set of points
    idx1 = 0;
    idx2 = idx1+1;

    front_p1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
    front_p2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);
    back_p1     = G4ThreeVector(detector[6+idx1][0],detector[6+idx1][1],detector[6+idx1][2]);
    back_p2     = G4ThreeVector(detector[6+idx2][0],detector[6+idx2][1],detector[6+idx2][2]);

    G4SubtractionSolid* can_volume_cut_2 = CutVolumeOnFourPoints(idx1, insideVol, volume_length, detector_phi, can_volume_cut_1, front_p1, front_p2, back_p1, back_p2);

    // next set of points
    idx1 = 1;
    idx2 = idx1+1;

    front_p1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
    front_p2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);
    back_p1     = G4ThreeVector(detector[6+idx1][0],detector[6+idx1][1],detector[6+idx1][2]);
    back_p2     = G4ThreeVector(detector[6+idx2][0],detector[6+idx2][1],detector[6+idx2][2]);

    G4SubtractionSolid* can_volume_cut_3 = CutVolumeOnFourPoints(idx1, insideVol, volume_length, detector_phi, can_volume_cut_2, front_p1, front_p2, back_p1, back_p2);

    idx1 = 2;
    idx2 = idx1+1;

    front_p1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
    front_p2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);
    back_p1     = G4ThreeVector(detector[6+idx1][0],detector[6+idx1][1],detector[6+idx1][2]);
    back_p2     = G4ThreeVector(detector[6+idx2][0],detector[6+idx2][1],detector[6+idx2][2]);

    G4SubtractionSolid* can_volume_cut_4 = CutVolumeOnFourPoints(idx1, insideVol, volume_length, detector_phi, can_volume_cut_3, front_p1, front_p2, back_p1, back_p2);

    idx1 = 3;
    idx2 = idx1+1;

    front_p1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
    front_p2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);
    back_p1     = G4ThreeVector(detector[6+idx1][0],detector[6+idx1][1],detector[6+idx1][2]);
    back_p2     = G4ThreeVector(detector[6+idx2][0],detector[6+idx2][1],detector[6+idx2][2]);

    G4SubtractionSolid* can_volume_cut_5 = CutVolumeOnFourPoints(idx1, insideVol, volume_length, detector_phi, can_volume_cut_4, front_p1, front_p2, back_p1, back_p2);

    idx1 = 4;
    idx2 = idx1+1;

    front_p1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
    front_p2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);
    back_p1     = G4ThreeVector(detector[6+idx1][0],detector[6+idx1][1],detector[6+idx1][2]);
    back_p2     = G4ThreeVector(detector[6+idx2][0],detector[6+idx2][1],detector[6+idx2][2]);

    G4SubtractionSolid* can_volume_cut_6 = CutVolumeOnFourPoints(idx1, insideVol, volume_length, detector_phi, can_volume_cut_5, front_p1, front_p2, back_p1, back_p2);

    idx1 = 5;
    idx2 = 0;

    front_p1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
    front_p2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);
    back_p1     = G4ThreeVector(detector[6+idx1][0],detector[6+idx1][1],detector[6+idx1][2]);
    back_p2     = G4ThreeVector(detector[6+idx2][0],detector[6+idx2][1],detector[6+idx2][2]);

    G4SubtractionSolid* can_volume_cut_7 = CutVolumeOnFourPoints(idx1, insideVol, volume_length, detector_phi, can_volume_cut_6, front_p1, front_p2, back_p1, back_p2);

    return can_volume_cut_7;
}

// This method does the actual cutting of the surface.
G4SubtractionSolid* DetectionSystemDescant::CutVolumeOnFourPoints(G4int idx, G4bool insideVol, G4double volume_length, G4double detector_phi[6], G4SubtractionSolid* volume, G4ThreeVector front_p1, G4ThreeVector front_p2, G4ThreeVector back_p1, G4ThreeVector back_p2)
{
    G4double maxZ = 500.0*mm;
    G4double theta = 0;
    G4double phi = 0;
    G4double move_z_cut = 0;

    // make the cuts excessively large.
    G4double cut_half_length_x = 1.0*maxZ;
    G4double cut_half_length_y = 1.0*maxZ;
    G4double cut_half_length_z = 2.0*maxZ;

    G4ThreeVector move_cut;
    G4RotationMatrix* rotate_cut;

    //    G4double front_mid_x = (front_p2.x()+front_p1.x())/2.0;
    //    G4double front_mid_y = (front_p2.y()+front_p1.y())/2.0;
    //    G4double back_mid_x = (back_p2.x()+back_p1.x())/2.0;
    //    G4double back_mid_y = (back_p2.y()+back_p1.y())/2.0;

    // find the equation of the line defined by the two front face points.
    // then using a circle of radius r, solve for the minimum r to get the minimum theta angle to that plane!
    G4double slope  = (front_p2.y() - front_p1.y())/(front_p2.x() - front_p1.x());
    G4double intcp  = (front_p1.y()) - (slope*front_p1.x());
    G4double minx   = (-1.0*slope*intcp)/(1+slope*slope);
    G4double miny   = slope*minx + intcp;

    // make vectors, and get angle between vectors
    G4ThreeVector new_v1 = G4ThreeVector(minx,miny,500.0*mm);
    G4ThreeVector new_v2 = G4ThreeVector(0.0*mm,0.0*mm,500.0*mm);

    theta   = acos((new_v1.dot(new_v2))/(new_v1.mag() * new_v2.mag()));
    phi     = detector_phi[idx];

    if(insideVol)
        move_z_cut = 1.0*(volume_length/2.0 + maxZ - (((this->can_thickness)/(sin(theta)) ) - this->can_thickness));
    else
        move_z_cut = 1.0*(volume_length/2.0 + maxZ);

    if(volume_length == lead_shield_thickness) // this is the lead shield, push it a bit closer so the tapers match up
        move_z_cut = move_z_cut - volume_length;

    G4Box* cut_plate = new G4Box("cut_plate", cut_half_length_x, cut_half_length_y, cut_half_length_z);

    move_cut = G4ThreeVector(cut_half_length_x/cos(theta),0.0*mm,move_z_cut);
    move_cut.rotateZ(1.0*phi);

    // rotate our cutting block to the "chopping" angle
    rotate_cut = new G4RotationMatrix;
    rotate_cut->rotateZ(-1.0*phi);
    rotate_cut->rotateY(1.0*theta);

    // snip snip
    G4SubtractionSolid* can_volume_cut = new G4SubtractionSolid("can_volume_cut", volume, cut_plate, rotate_cut, move_cut);

    return can_volume_cut;
}

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemDescant::GetDirectionXYZ(G4double theta, G4double phi)
{
    G4double x,y,z;
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);

    G4ThreeVector direction = G4ThreeVector(x,y,z);

    return direction;
}

G4double DetectionSystemDescant::transX(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*cos(phi) );
}

G4double DetectionSystemDescant::transY(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*sin(phi) );
}

G4double DetectionSystemDescant::transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi){
    return ( pow(x*x+y*y+z*z,0.5)*cos(theta) );
}

G4ThreeVector DetectionSystemDescant::SolveLineEquationX(G4ThreeVector p1, G4ThreeVector p2, G4double x){
    G4ThreeVector v, result;
    v = p2 - p1;
    G4double t = (x - p1.x())/(v.x());
    result = G4ThreeVector(x,p1.y()+v.y()*t,p1.z()+v.z()*t);
    return result;
}

G4ThreeVector DetectionSystemDescant::SolveLineEquationY(G4ThreeVector p1, G4ThreeVector p2, G4double y){
    G4ThreeVector v, result;
    v = p2 - p1;
    G4double t = (y - p1.y())/(v.y());
    result = G4ThreeVector(p1.x()+v.x()*t,y,p1.z()+v.z()*t);
    return result;
}

G4ThreeVector DetectionSystemDescant::SolveLineEquationZ(G4ThreeVector p1, G4ThreeVector p2, G4double z){
    G4ThreeVector v, result;
    v = p2 - p1;
    G4double t = (z - p1.z())/(v.z());
    result = G4ThreeVector(p1.x()+v.x()*t,p1.y()+v.y()*t,z);
    return result;
}
