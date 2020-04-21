#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "Global.hh"

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

#include <string>

DetectionSystemDescant::DetectionSystemDescant(G4bool leadShield) :
	// LogicalVolumes
	fBlueVolumeLog(0),
	fGreenVolumeLog(0),
	fRedVolumeLog(0),
	fWhiteVolumeLog(0),
	fYellowVolumeLog(0),
	fBlueScintillatorVolumeLog(0),
	fGreenScintillatorVolumeLog(0),
	fRedScintillatorVolumeLog(0),
	fWhiteScintillatorVolumeLog(0),
	fYellowScintillatorVolumeLog(0),
	fBlueLeadVolumeLog(0),
	fGreenLeadVolumeLog(0),
	fRedLeadVolumeLog(0),
	fWhiteLeadVolumeLog(0),
	fYellowLeadVolumeLog(0),
	fQuartzWindow5inchLog(0),
	fQuartzWindow3inchLog(0)
{
	// can properties
	fCanLength              = 150.0*mm;
	fCanThickness           = 1.5*mm;
	fCanMaterial            = "G4_Al";
	fLiquidMaterial         = "BC537";
	fLeadMaterial           = "G4_Pb";
	fLeadShieldThickness   = 6.35*mm;
	fCanBackThickness      = 12.7*mm;   

	// set DESCANT lead shield on/off (set by run macro)
	fIncludeLead = leadShield;

	// for cutting the PMT port on the detectors
	fInnerRadiusWindow       = 0.0*mm;     
	fOuterRadiusWindow5inch  = (127.0*mm)/2.0; 
	fOuterRadiusWindow3inch  = (76.0*mm)/2.0;  
	fHalfLengthZWindowCut    = 100.0*mm;
	fStartPhi                = 0.0*deg;
	fEndPhi                  = 360.0*deg;

	fWhitePMTOffsetY         = 4.0*mm;
	fBluePMTOffsetX          = 3.0*mm;
	fRedPMTOffsetX           = 12.7*mm;
	fYellowGreenPMTOffsetX   = 10.0*mm; 
	fYellowGreenPMTOffsetY   = 5.0*mm;

	fQuartzMaterial                 = "G4_SILICON_DIOXIDE";
	fOpticalWindowHalfThickness   = (9.0*mm)/2.0;

	// The face of the detector is 50 cm from the origin, this does NOT include the lead shield.
	fRadialDistance         = 50*cm;

	// Check surfaces to determine any problematic overlaps. Turn this on to have Geant4 check the surfaces.
	// Do not leave this on, it will slow the DetectorConstruction process!
	// This was last check on May 29th, 2015. - GOOD!, but...
	// The green, yellow and blue detectors had small overlaps, this is probably due to taking the pure mathematical model
	// of DESCANT, and then Saint Gobain rounding the result to generate the data files (below).
	fSurfCheck               = false;

	// Geometry overlaps were found for green, yellow and blue detectors.               
	// We need to squeeze the green and yellow detectors along the y direction,                
	// and we need to squeeze the blue detectors on one x edge of the blue detector.           
	fTrimGreenY    = 10.0*um;         
	fTrimYellowY   = 10.0*um;         
	fTrimBlueX     = 2.2*um;

	// Saint Gobain data files, 6 points around the front face of the can, and 6 on the back
	// from Dan Brennan

	G4double blue[12][3] = {
		{ 40.60*mm,  61.10*mm,    0.00*mm},
		{-34.70*mm,  61.10*mm,    0.00*mm},
		{-75.10*mm,   0.00*mm,    0.00*mm},
		{-34.70*mm, -61.10*mm,    0.00*mm},
		{ 40.60*mm, -61.10*mm,    0.00*mm},
		{ 76.30*mm-fTrimBlueX,   0.00*mm,    0.00*mm},
		{ 52.80*mm,  79.40*mm, -150.00*mm},
		{-45.10*mm,  79.40*mm, -150.00*mm},
		{-97.60*mm,   0.00*mm, -150.00*mm},
		{-45.10*mm, -79.40*mm, -150.00*mm},
		{ 52.80*mm, -79.40*mm, -150.00*mm},
		{ 99.20*mm-fTrimBlueX,   0.00*mm, -150.00*mm}
	};

	G4double green[12][3] = {
		{ 31.90*mm,  61.20*mm-fTrimGreenY,    0.00*mm},
		{-31.90*mm,  61.20*mm-fTrimGreenY,    0.00*mm},
		{-71.50*mm,   0.00*mm,    0.00*mm},
		{-31.90*mm, -55.30*mm+fTrimGreenY,    0.00*mm},
		{ 31.90*mm, -55.30*mm+fTrimGreenY,    0.00*mm},
		{ 47.90*mm,  36.60*mm,    0.00*mm},
		{ 41.50*mm,  79.60*mm-fTrimGreenY, -150.00*mm},
		{-41.50*mm,  79.60*mm-fTrimGreenY, -150.00*mm},
		{-93.00*mm,   0.00*mm, -150.00*mm},
		{-41.50*mm, -71.90*mm+fTrimGreenY, -150.00*mm},
		{ 41.50*mm, -71.90*mm+fTrimGreenY, -150.00*mm},
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
		{ 31.90*mm,   61.20*mm-fTrimYellowY,    0.00*mm},
		{-31.90*mm,   61.20*mm-fTrimYellowY,    0.00*mm},
		{-47.90*mm,   36.60*mm,    0.00*mm},
		{-31.90*mm,  -55.30*mm+fTrimYellowY,    0.00*mm},
		{ 31.90*mm,  -55.30*mm+fTrimYellowY,    0.00*mm},
		{ 71.50*mm,    0.00*mm,    0.00*mm},
		{ 41.50*mm,   79.60*mm-fTrimYellowY, -150.00*mm},
		{-41.50*mm,   79.60*mm-fTrimYellowY, -150.00*mm},
		{-62.30*mm,   47.60*mm, -150.00*mm},
		{-41.50*mm,  -71.90*mm+fTrimYellowY, -150.00*mm},
		{ 41.50*mm,  -71.90*mm+fTrimYellowY, -150.00*mm},
		{ 93.00*mm,    0.00*mm, -150.00*mm}
	};


	memcpy(fBlueDetector,   blue,   sizeof(fBlueDetector));
	memcpy(fGreenDetector,  green,  sizeof(fGreenDetector));
	memcpy(fRedDetector,    red,    sizeof(fRedDetector));
	memcpy(fWhiteDetector,  white,  sizeof(fWhiteDetector));
	memcpy(fYellowDetector, yellow, sizeof(fYellowDetector));

	// These are the angles of the cuts for each of the 6 sides of the detectors
	// These were worked out by hand.
	G4double bluePhiArray[6] = {
		(90.0*deg),
		(180.0*deg-(33.4731*deg)),
		-1.0*(180.0*deg-(33.4731*deg)),
		-1.0*(90.0*deg),
		-1.0*(30.2972*deg),
		(30.2972*deg)
	};

	G4double greenPhiArray[6] = {
		(90.0*deg),
		(180.0*deg-(32.9052*deg)),
		-1.0*(180.0*deg-(35.6062*deg)),
		-1.0*(90.0*deg),
		-1.0*(9.8763*deg),
		(33.0402*deg)
	};

	G4double redPhiArray[6] = {
		(90.0*deg + 2.4242*deg),
		(180.0*deg-(31.0557*deg)),
		-1.0*(180.0*deg-(31.0557*deg)),
		-1.0*(90.0*deg + 2.4242*deg),
		-1.0*(22.2424*deg),
		(22.2424*deg)
	};

	G4double whitePhiArray[6] = {
		(90.0*deg),
		(180.0*deg-(32.9052*deg)),
		-1.0*(180.0*deg-(35.6062*deg)),
		-1.0*(90.0*deg),
		-1.0*(35.6062*deg),
		(32.9052*deg)
	};

	G4double yellowPhiArray[6] = {
		(90.0*deg),
		(180.0*deg-(33.0402*deg)),
		-1.0*(180.0*deg-(9.8763*deg)),
		-1.0*(90.0*deg),
		-1.0*(35.6062*deg),
		1.0*(32.9052*deg)
	};

	memcpy(fBluePhi,    bluePhiArray,     sizeof(fBluePhi));
	memcpy(fGreenPhi,   greenPhiArray,    sizeof(fGreenPhi));
	memcpy(fRedPhi,     redPhiArray,      sizeof(fRedPhi));
	memcpy(fWhitePhi,   whitePhiArray,    sizeof(fWhitePhi));
	memcpy(fYellowPhi,  yellowPhiArray,   sizeof(fYellowPhi));

	// The Euler angles from James' MSc thesis which gives us the detector positions
	// Some of the angles for the green and yellow detectors are wrong in James' thesis,
	// note the +180 on a few angles.
	G4double blueAlphaBetaGammaArray[15][3] = {
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

	G4double greenAlphaBetaGammaArray[10][3] = {
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

	G4double redAlphaBetaGammaArray[15][3] = {
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

	G4double whiteAlphaBetaGammaArray[20][3] = {
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

	G4double yellowAlphaBetaGammaArray[10][3] = {
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

	memcpy(fBlueAlphaBetaGamma,   blueAlphaBetaGammaArray,   sizeof(fBlueAlphaBetaGamma));
	memcpy(fGreenAlphaBetaGamma,  greenAlphaBetaGammaArray,  sizeof(fGreenAlphaBetaGamma));
	memcpy(fRedAlphaBetaGamma,    redAlphaBetaGammaArray,    sizeof(fRedAlphaBetaGamma));
	memcpy(fWhiteAlphaBetaGamma,  whiteAlphaBetaGammaArray,  sizeof(fWhiteAlphaBetaGamma));
	memcpy(fYellowAlphaBetaGamma, yellowAlphaBetaGammaArray, sizeof(fYellowAlphaBetaGamma));

	fBlueColour     = G4Colour(0.0/255.0,0.0/255.0,255.0/255.0);
	fGreenColour    = G4Colour(0.0/255.0,255.0/255.0,0.0/255.0);
	fRedColour      = G4Colour(255.0/255.0,0.0/255.0,0.0/255.0);
	fWhiteColour    = G4Colour(255.0/255.0,255.0/255.0,255.0/255.0);
	fYellowColour   = G4Colour(255.0/255.0,255.0/255.0,0.0/255.0);
	fLiquidColour   = G4Colour(0.0/255.0,255.0/255.0,225.0/255.0);
	fGreyColour     = G4Colour(0.5, 0.5, 0.5); 
	fMagentaColour  = G4Colour(1.0, 0.0, 1.0); 
	fBlackColour    = G4Colour(0.0, 0.0, 0.0);

	// for lanthanum bromide locations
	//G4double triangleThetaAngle = (180/M_PI)*(atan((1/sqrt(3))/sqrt((11/12) + (1/sqrt(2))))+atan((sqrt(2))/(1+sqrt(2))))*deg;
	G4double triangleThetaAngle = 54.735610317245360*deg;

	// theta
	fDetectorAngles[0][0] 	= triangleThetaAngle;
	fDetectorAngles[1][0] 	= triangleThetaAngle;
	fDetectorAngles[2][0] 	= triangleThetaAngle;
	fDetectorAngles[3][0] 	= triangleThetaAngle;
	fDetectorAngles[4][0] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[5][0] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[6][0] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[7][0] 	= 180.0*deg - triangleThetaAngle;
	// phi
	fDetectorAngles[0][1] 	= 22.5*deg;
	fDetectorAngles[1][1] 	= 112.5*deg;
	fDetectorAngles[2][1] 	= 202.5*deg;
	fDetectorAngles[3][1] 	= 292.5*deg;
	fDetectorAngles[4][1] 	= 22.5*deg;
	fDetectorAngles[5][1] 	= 112.5*deg;
	fDetectorAngles[6][1] 	= 202.5*deg;
	fDetectorAngles[7][1] 	= 292.5*deg;
	// yaw (alpha)
	fDetectorAngles[0][2] 	= 0.0*deg;
	fDetectorAngles[1][2] 	= 0.0*deg;
	fDetectorAngles[2][2] 	= 0.0*deg;
	fDetectorAngles[3][2] 	= 0.0*deg;
	fDetectorAngles[4][2] 	= 0.0*deg;
	fDetectorAngles[5][2] 	= 0.0*deg;
	fDetectorAngles[6][2] 	= 0.0*deg;
	fDetectorAngles[7][2] 	= 0.0*deg;
	// pitch (beta)
	fDetectorAngles[0][3] 	= triangleThetaAngle;
	fDetectorAngles[1][3] 	= triangleThetaAngle;
	fDetectorAngles[2][3] 	= triangleThetaAngle;
	fDetectorAngles[3][3] 	= triangleThetaAngle;
	fDetectorAngles[4][3] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[5][3] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[6][3] 	= 180.0*deg - triangleThetaAngle;
	fDetectorAngles[7][3] 	= 180.0*deg - triangleThetaAngle;
	// roll (gamma)
	fDetectorAngles[0][4] 	= 22.5*deg;
	fDetectorAngles[1][4] 	= 112.5*deg;
	fDetectorAngles[2][4] 	= 202.5*deg;
	fDetectorAngles[3][4] 	= 292.5*deg;
	fDetectorAngles[4][4] 	= 22.5*deg;
	fDetectorAngles[5][4] 	= 112.5*deg;
	fDetectorAngles[6][4] 	= 202.5*deg;
	fDetectorAngles[7][4] 	= 292.5*deg;
}


DetectionSystemDescant::~DetectionSystemDescant() {
	// LogicalVolumes
	delete fBlueVolumeLog;
	delete fGreenVolumeLog;
	delete fRedVolumeLog;
	delete fWhiteVolumeLog;
	delete fYellowVolumeLog;
	delete fBlueScintillatorVolumeLog;
	delete fGreenScintillatorVolumeLog;
	delete fRedScintillatorVolumeLog;
	delete fWhiteScintillatorVolumeLog;
	delete fYellowScintillatorVolumeLog;
	delete fBlueLeadVolumeLog;
	delete fGreenLeadVolumeLog;
	delete fRedLeadVolumeLog;
	delete fWhiteLeadVolumeLog;
	delete fYellowLeadVolumeLog;
	delete fQuartzWindow5inchLog;
	delete fQuartzWindow3inchLog;
}


G4int DetectionSystemDescant::Build() {
	// Build assembly volumes
	// assemblyBlue contains non-sensitive volumes of the blue detector
	// assemblyBlueScintillator contains only sensitive volumes of the blue detector
	fAssemblyBlue = new G4AssemblyVolume();
	fAssemblyBlueScintillator = new G4AssemblyVolume();

	fAssemblyGreen = new G4AssemblyVolume();
	fAssemblyGreenScintillator = new G4AssemblyVolume();

	fAssemblyRed = new G4AssemblyVolume();
	fAssemblyRedScintillator = new G4AssemblyVolume();

	fAssemblyWhite = new G4AssemblyVolume();
	fAssemblyWhiteScintillator = new G4AssemblyVolume();

	fAssemblyYellow = new G4AssemblyVolume();
	fAssemblyYellowScintillator = new G4AssemblyVolume();

	G4cout<<"Building DESCANT..."<<G4endl;

	G4cout<<"BuildCanVolume"<<G4endl;
	BuildCanVolume();

	G4cout<<"BuildDetectorVolume"<<G4endl;
	BuildDetectorVolume();

	return 1;
}

G4int DetectionSystemDescant::PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber) {
	G4double position = fRadialDistance;

	G4double x = 0;
	G4double y = 0;
	G4double z = position;

	G4ThreeVector move = G4ThreeVector(x,y,z);
	G4RotationMatrix rotate;
	rotate.rotateY(M_PI); // flip the detector so that the face is pointing upstream.

	for(int i = 0; i < detectorNumber; ++i) {
		// brute force map to DESCANT positions, could be improved by creating a single array of angles?
		switch(i) {
			case 0:
				move.rotateZ(fWhiteAlphaBetaGamma[8][2]);
				move.rotateY(fWhiteAlphaBetaGamma[8][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[8][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[8][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[8][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[8][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 1:
				move.rotateZ(fWhiteAlphaBetaGamma[4][2]);
				move.rotateY(fWhiteAlphaBetaGamma[4][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[4][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[4][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[4][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[4][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 2:
				move.rotateZ(fWhiteAlphaBetaGamma[0][2]);
				move.rotateY(fWhiteAlphaBetaGamma[0][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[0][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[0][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[0][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[0][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 3:
				move.rotateZ(fWhiteAlphaBetaGamma[16][2]);
				move.rotateY(fWhiteAlphaBetaGamma[16][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[16][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[16][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[16][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[16][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 4:
				move.rotateZ(fWhiteAlphaBetaGamma[12][2]);
				move.rotateY(fWhiteAlphaBetaGamma[12][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[12][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[12][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[12][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[12][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 5:
				move.rotateZ(fWhiteAlphaBetaGamma[9][2]);
				move.rotateY(fWhiteAlphaBetaGamma[9][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[9][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[9][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[9][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[9][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 6:
				move.rotateZ(fRedAlphaBetaGamma[3][2]);
				move.rotateY(fRedAlphaBetaGamma[3][1]);
				move.rotateZ(fRedAlphaBetaGamma[3][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[3][2]);
				rotate.rotateY(fRedAlphaBetaGamma[3][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[3][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 7:
				move.rotateZ(fWhiteAlphaBetaGamma[5][2]);
				move.rotateY(fWhiteAlphaBetaGamma[5][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[5][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[5][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[5][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[5][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 8:
				move.rotateZ(fRedAlphaBetaGamma[0][2]);
				move.rotateY(fRedAlphaBetaGamma[0][1]);
				move.rotateZ(fRedAlphaBetaGamma[0][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[0][2]);
				rotate.rotateY(fRedAlphaBetaGamma[0][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[0][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 9:
				move.rotateZ(fWhiteAlphaBetaGamma[1][2]);
				move.rotateY(fWhiteAlphaBetaGamma[1][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[1][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[1][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[1][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[1][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 10:
				move.rotateZ(fRedAlphaBetaGamma[12][2]);
				move.rotateY(fRedAlphaBetaGamma[12][1]);
				move.rotateZ(fRedAlphaBetaGamma[12][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[12][2]);
				rotate.rotateY(fRedAlphaBetaGamma[12][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[12][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 11:
				move.rotateZ(fWhiteAlphaBetaGamma[17][2]);
				move.rotateY(fWhiteAlphaBetaGamma[17][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[17][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[17][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[17][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[17][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 12:
				move.rotateZ(fRedAlphaBetaGamma[9][2]);
				move.rotateY(fRedAlphaBetaGamma[9][1]);
				move.rotateZ(fRedAlphaBetaGamma[9][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[9][2]);
				rotate.rotateY(fRedAlphaBetaGamma[9][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[9][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 13:
				move.rotateZ(fWhiteAlphaBetaGamma[13][2]);
				move.rotateY(fWhiteAlphaBetaGamma[13][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[13][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[13][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[13][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[13][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 14:
				move.rotateZ(fRedAlphaBetaGamma[6][2]);
				move.rotateY(fRedAlphaBetaGamma[6][1]);
				move.rotateZ(fRedAlphaBetaGamma[6][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[6][2]);
				rotate.rotateY(fRedAlphaBetaGamma[6][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[6][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 15:
				move.rotateZ(fWhiteAlphaBetaGamma[10][2]);
				move.rotateY(fWhiteAlphaBetaGamma[10][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[10][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[10][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[10][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[10][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 16:
				move.rotateZ(fBlueAlphaBetaGamma[4][2]);
				move.rotateY(fBlueAlphaBetaGamma[4][1]);
				move.rotateZ(fBlueAlphaBetaGamma[4][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[4][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[4][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[4][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 17:
				move.rotateZ(fBlueAlphaBetaGamma[3][2]);
				move.rotateY(fBlueAlphaBetaGamma[3][1]);
				move.rotateZ(fBlueAlphaBetaGamma[3][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[3][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[3][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[3][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 18:
				move.rotateZ(fWhiteAlphaBetaGamma[6][2]);
				move.rotateY(fWhiteAlphaBetaGamma[6][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[6][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[6][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[6][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[6][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 19:
				move.rotateZ(fBlueAlphaBetaGamma[1][2]);
				move.rotateY(fBlueAlphaBetaGamma[1][1]);
				move.rotateZ(fBlueAlphaBetaGamma[1][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[1][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[1][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[1][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 20:
				move.rotateZ(fBlueAlphaBetaGamma[0][2]);
				move.rotateY(fBlueAlphaBetaGamma[0][1]);
				move.rotateZ(fBlueAlphaBetaGamma[0][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[0][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[0][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[0][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 21:
				move.rotateZ(fWhiteAlphaBetaGamma[2][2]);
				move.rotateY(fWhiteAlphaBetaGamma[2][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[2][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[2][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[2][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[2][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 22:
				move.rotateZ(fBlueAlphaBetaGamma[13][2]);
				move.rotateY(fBlueAlphaBetaGamma[13][1]);
				move.rotateZ(fBlueAlphaBetaGamma[13][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[13][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[13][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[13][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 23:
				move.rotateZ(fBlueAlphaBetaGamma[12][2]);
				move.rotateY(fBlueAlphaBetaGamma[12][1]);
				move.rotateZ(fBlueAlphaBetaGamma[12][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[12][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[12][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[12][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 24:
				move.rotateZ(fWhiteAlphaBetaGamma[18][2]);
				move.rotateY(fWhiteAlphaBetaGamma[18][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[18][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[18][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[18][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[18][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 25:
				move.rotateZ(fBlueAlphaBetaGamma[10][2]);
				move.rotateY(fBlueAlphaBetaGamma[10][1]);
				move.rotateZ(fBlueAlphaBetaGamma[10][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[10][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[10][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[10][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 26:
				move.rotateZ(fBlueAlphaBetaGamma[9][2]);
				move.rotateY(fBlueAlphaBetaGamma[9][1]);
				move.rotateZ(fBlueAlphaBetaGamma[9][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[9][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[9][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[9][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 27:
				move.rotateZ(fWhiteAlphaBetaGamma[14][2]);
				move.rotateY(fWhiteAlphaBetaGamma[14][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[14][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[14][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[14][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[14][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 28:
				move.rotateZ(fBlueAlphaBetaGamma[7][2]);
				move.rotateY(fBlueAlphaBetaGamma[7][1]);
				move.rotateZ(fBlueAlphaBetaGamma[7][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[7][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[7][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[7][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 29:
				move.rotateZ(fBlueAlphaBetaGamma[6][2]);
				move.rotateY(fBlueAlphaBetaGamma[6][1]);
				move.rotateZ(fBlueAlphaBetaGamma[6][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[6][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[6][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[6][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 30:
				move.rotateZ(fWhiteAlphaBetaGamma[11][2]);
				move.rotateY(fWhiteAlphaBetaGamma[11][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[11][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[11][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[11][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[11][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 31:
				move.rotateZ(fRedAlphaBetaGamma[5][2]);
				move.rotateY(fRedAlphaBetaGamma[5][1]);
				move.rotateZ(fRedAlphaBetaGamma[5][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[5][2]);
				rotate.rotateY(fRedAlphaBetaGamma[5][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[5][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 32:
				move.rotateZ(fBlueAlphaBetaGamma[5][2]);
				move.rotateY(fBlueAlphaBetaGamma[5][1]);
				move.rotateZ(fBlueAlphaBetaGamma[5][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[5][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[5][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[5][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 33:
				move.rotateZ(fRedAlphaBetaGamma[4][2]);
				move.rotateY(fRedAlphaBetaGamma[4][1]);
				move.rotateZ(fRedAlphaBetaGamma[4][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[4][2]);
				rotate.rotateY(fRedAlphaBetaGamma[4][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[4][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 34:
				move.rotateZ(fWhiteAlphaBetaGamma[7][2]);
				move.rotateY(fWhiteAlphaBetaGamma[7][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[7][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[7][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[7][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[7][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 35:
				move.rotateZ(fRedAlphaBetaGamma[2][2]);
				move.rotateY(fRedAlphaBetaGamma[2][1]);
				move.rotateZ(fRedAlphaBetaGamma[2][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[2][2]);
				rotate.rotateY(fRedAlphaBetaGamma[2][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[2][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 36:
				move.rotateZ(fBlueAlphaBetaGamma[2][2]);
				move.rotateY(fBlueAlphaBetaGamma[2][1]);
				move.rotateZ(fBlueAlphaBetaGamma[2][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[2][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[2][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[2][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 37:
				move.rotateZ(fRedAlphaBetaGamma[1][2]);
				move.rotateY(fRedAlphaBetaGamma[1][1]);
				move.rotateZ(fRedAlphaBetaGamma[1][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[1][2]);
				rotate.rotateY(fRedAlphaBetaGamma[1][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[1][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 38:
				move.rotateZ(fWhiteAlphaBetaGamma[3][2]);
				move.rotateY(fWhiteAlphaBetaGamma[3][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[3][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[3][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[3][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[3][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 39:
				move.rotateZ(fRedAlphaBetaGamma[14][2]);
				move.rotateY(fRedAlphaBetaGamma[14][1]);
				move.rotateZ(fRedAlphaBetaGamma[14][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[14][2]);
				rotate.rotateY(fRedAlphaBetaGamma[14][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[14][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 40:
				move.rotateZ(fBlueAlphaBetaGamma[14][2]);
				move.rotateY(fBlueAlphaBetaGamma[14][1]);
				move.rotateZ(fBlueAlphaBetaGamma[14][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[14][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[14][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[14][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 41:
				move.rotateZ(fRedAlphaBetaGamma[13][2]);
				move.rotateY(fRedAlphaBetaGamma[13][1]);
				move.rotateZ(fRedAlphaBetaGamma[13][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[13][2]);
				rotate.rotateY(fRedAlphaBetaGamma[13][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[13][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 42:
				move.rotateZ(fWhiteAlphaBetaGamma[19][2]);
				move.rotateY(fWhiteAlphaBetaGamma[19][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[19][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[19][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[19][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[19][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 43:
				move.rotateZ(fRedAlphaBetaGamma[10][2]);
				move.rotateY(fRedAlphaBetaGamma[10][1]);
				move.rotateZ(fRedAlphaBetaGamma[10][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[10][2]);
				rotate.rotateY(fRedAlphaBetaGamma[10][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[10][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 44:
				move.rotateZ(fBlueAlphaBetaGamma[11][2]);
				move.rotateY(fBlueAlphaBetaGamma[11][1]);
				move.rotateZ(fBlueAlphaBetaGamma[11][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[11][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[11][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[11][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 45:
				move.rotateZ(fRedAlphaBetaGamma[11][2]);
				move.rotateY(fRedAlphaBetaGamma[11][1]);
				move.rotateZ(fRedAlphaBetaGamma[11][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[11][2]);
				rotate.rotateY(fRedAlphaBetaGamma[11][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[11][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 46:
				move.rotateZ(fWhiteAlphaBetaGamma[15][2]);
				move.rotateY(fWhiteAlphaBetaGamma[15][1]);
				move.rotateZ(fWhiteAlphaBetaGamma[15][0]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[15][2]);
				rotate.rotateY(fWhiteAlphaBetaGamma[15][1]);
				rotate.rotateZ(fWhiteAlphaBetaGamma[15][0]);
				fAssemblyWhite->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 47:
				move.rotateZ(fRedAlphaBetaGamma[8][2]);
				move.rotateY(fRedAlphaBetaGamma[8][1]);
				move.rotateZ(fRedAlphaBetaGamma[8][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[8][2]);
				rotate.rotateY(fRedAlphaBetaGamma[8][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[8][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 48:
				move.rotateZ(fBlueAlphaBetaGamma[8][2]);
				move.rotateY(fBlueAlphaBetaGamma[8][1]);
				move.rotateZ(fBlueAlphaBetaGamma[8][0]);
				rotate.rotateZ(fBlueAlphaBetaGamma[8][2]);
				rotate.rotateY(fBlueAlphaBetaGamma[8][1]);
				rotate.rotateZ(fBlueAlphaBetaGamma[8][0]);
				fAssemblyBlue->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyBlueScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 49:
				move.rotateZ(fRedAlphaBetaGamma[7][2]);
				move.rotateY(fRedAlphaBetaGamma[7][1]);
				move.rotateZ(fRedAlphaBetaGamma[7][0]);
				rotate.rotateZ(fRedAlphaBetaGamma[7][2]);
				rotate.rotateY(fRedAlphaBetaGamma[7][1]);
				rotate.rotateZ(fRedAlphaBetaGamma[7][0]);
				fAssemblyRed->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyRedScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 50:
				move.rotateZ(fGreenAlphaBetaGamma[4][2]);
				move.rotateY(fGreenAlphaBetaGamma[4][1]);
				move.rotateZ(fGreenAlphaBetaGamma[4][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[4][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[4][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[4][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 51:
				move.rotateZ(fYellowAlphaBetaGamma[3][2]);
				move.rotateY(fYellowAlphaBetaGamma[3][1]);
				move.rotateZ(fYellowAlphaBetaGamma[3][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[3][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[3][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[3][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 52:
				move.rotateZ(fYellowAlphaBetaGamma[2][2]);
				move.rotateY(fYellowAlphaBetaGamma[2][1]);
				move.rotateZ(fYellowAlphaBetaGamma[2][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[2][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[2][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[2][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 53:
				move.rotateZ(fGreenAlphaBetaGamma[3][2]);
				move.rotateY(fGreenAlphaBetaGamma[3][1]);
				move.rotateZ(fGreenAlphaBetaGamma[3][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[3][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[3][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[3][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 54:
				move.rotateZ(fGreenAlphaBetaGamma[2][2]);
				move.rotateY(fGreenAlphaBetaGamma[2][1]);
				move.rotateZ(fGreenAlphaBetaGamma[2][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[2][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[2][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[2][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 55:
				move.rotateZ(fYellowAlphaBetaGamma[1][2]);
				move.rotateY(fYellowAlphaBetaGamma[1][1]);
				move.rotateZ(fYellowAlphaBetaGamma[1][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[1][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[1][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[1][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 56:
				move.rotateZ(fYellowAlphaBetaGamma[0][2]);
				move.rotateY(fYellowAlphaBetaGamma[0][1]);
				move.rotateZ(fYellowAlphaBetaGamma[0][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[0][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[0][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[0][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 57:
				move.rotateZ(fGreenAlphaBetaGamma[1][2]);
				move.rotateY(fGreenAlphaBetaGamma[1][1]);
				move.rotateZ(fGreenAlphaBetaGamma[1][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[1][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[1][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[1][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 58:
				move.rotateZ(fGreenAlphaBetaGamma[0][2]);
				move.rotateY(fGreenAlphaBetaGamma[0][1]);
				move.rotateZ(fGreenAlphaBetaGamma[0][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[0][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[0][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[0][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 59:
				move.rotateZ(fYellowAlphaBetaGamma[9][2]);
				move.rotateY(fYellowAlphaBetaGamma[9][1]);
				move.rotateZ(fYellowAlphaBetaGamma[9][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[9][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[9][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[9][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 60:
				move.rotateZ(fYellowAlphaBetaGamma[8][2]);
				move.rotateY(fYellowAlphaBetaGamma[8][1]);
				move.rotateZ(fYellowAlphaBetaGamma[8][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[8][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[8][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[8][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 61:
				move.rotateZ(fGreenAlphaBetaGamma[9][2]);
				move.rotateY(fGreenAlphaBetaGamma[9][1]);
				move.rotateZ(fGreenAlphaBetaGamma[9][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[9][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[9][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[9][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 62:
				move.rotateZ(fGreenAlphaBetaGamma[8][2]);
				move.rotateY(fGreenAlphaBetaGamma[8][1]);
				move.rotateZ(fGreenAlphaBetaGamma[8][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[8][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[8][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[8][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 63:
				move.rotateZ(fYellowAlphaBetaGamma[7][2]);
				move.rotateY(fYellowAlphaBetaGamma[7][1]);
				move.rotateZ(fYellowAlphaBetaGamma[7][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[7][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[7][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[7][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 64:
				move.rotateZ(fYellowAlphaBetaGamma[6][2]);
				move.rotateY(fYellowAlphaBetaGamma[6][1]);
				move.rotateZ(fYellowAlphaBetaGamma[6][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[6][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[6][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[6][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 65:
				move.rotateZ(fGreenAlphaBetaGamma[7][2]);
				move.rotateY(fGreenAlphaBetaGamma[7][1]);
				move.rotateZ(fGreenAlphaBetaGamma[7][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[7][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[7][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[7][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 66:
				move.rotateZ(fGreenAlphaBetaGamma[6][2]);
				move.rotateY(fGreenAlphaBetaGamma[6][1]);
				move.rotateZ(fGreenAlphaBetaGamma[6][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[6][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[6][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[6][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 67:
				move.rotateZ(fYellowAlphaBetaGamma[5][2]);
				move.rotateY(fYellowAlphaBetaGamma[5][1]);
				move.rotateZ(fYellowAlphaBetaGamma[5][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[5][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[5][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[5][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 68:
				move.rotateZ(fYellowAlphaBetaGamma[4][2]);
				move.rotateY(fYellowAlphaBetaGamma[4][1]);
				move.rotateZ(fYellowAlphaBetaGamma[4][0]);
				rotate.rotateZ(fYellowAlphaBetaGamma[4][2]);
				rotate.rotateY(fYellowAlphaBetaGamma[4][1]);
				rotate.rotateZ(fYellowAlphaBetaGamma[4][0]);
				fAssemblyYellow->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyYellowScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			case 69:
				move.rotateZ(fGreenAlphaBetaGamma[5][2]);
				move.rotateY(fGreenAlphaBetaGamma[5][1]);
				move.rotateZ(fGreenAlphaBetaGamma[5][0]);
				rotate.rotateZ(fGreenAlphaBetaGamma[5][2]);
				rotate.rotateY(fGreenAlphaBetaGamma[5][1]);
				rotate.rotateZ(fGreenAlphaBetaGamma[5][0]);
				fAssemblyGreen->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				fAssemblyGreenScintillator->MakeImprint(expHallLog, move, &rotate, i, fSurfCheck);
				break;
			default:
				G4cout<<"Error, unknown DESCANT detector number "<<i<<"/"<<detectorNumber<<std::endl;
				break;
		}
		// reset move vector and rotation matrix

		move.set(x,y,z);
		rotate = G4RotationMatrix::IDENTITY;
		rotate.rotateY(M_PI); // flip the detector so that the face is pointing upstream.
	}

	return 1;
}


G4int DetectionSystemDescant::PlaceDetector(G4LogicalVolume* expHallLog, G4String color, G4ThreeVector pos, G4ThreeVector rot) {
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
		fAssemblyBlue->MakeImprint(expHallLog, pos, rotate, 1);
		fAssemblyBlueScintillator->MakeImprint(expHallLog, pos, rotate, 1);
	} else if(color.compareTo("red", G4String::ignoreCase) == 0) {
		fAssemblyRed->MakeImprint(expHallLog, pos, rotate, 1);
		fAssemblyRedScintillator->MakeImprint(expHallLog, pos, rotate, 1);
	} else if(color.compareTo("white", G4String::ignoreCase) == 0) {
		fAssemblyWhite->MakeImprint(expHallLog, pos, rotate, 1);
		fAssemblyWhiteScintillator->MakeImprint(expHallLog, pos, rotate, 1);
	} else if(color.compareTo("green", G4String::ignoreCase) == 0) {
		fAssemblyGreen->MakeImprint(expHallLog, pos, rotate, 1);
		fAssemblyGreenScintillator->MakeImprint(expHallLog, pos, rotate, 1);
	} else if(color.compareTo("yellow", G4String::ignoreCase) == 0) {
		fAssemblyYellow->MakeImprint(expHallLog, pos, rotate, 1);
		fAssemblyYellowScintillator->MakeImprint(expHallLog, pos, rotate, 1);
	}

	return 1;
}

G4int DetectionSystemDescant::PlaceDetectorAuxPorts(G4LogicalVolume* expHallLog, G4int detectorNumber, G4double radialpos) {
	G4int detectorCopyID = 0;

	G4cout<<"Descant Detector Number = "<<detectorNumber<<G4endl;

	G4int copyNumber = detectorCopyID + detectorNumber;

	G4double position = radialpos + fCanThicknessFront +(fCanLength / 2.0) ;
	fSetRadialPos = radialpos;

	G4double theta  = fDetectorAngles[detectorNumber][0];
	G4double phi    = fDetectorAngles[detectorNumber][1];
	//G4double alpha  = fDetectorAngles[detectorNumber][2]; // yaw
	G4double beta   = fDetectorAngles[detectorNumber][3]; // pitch
	G4double gamma  = fDetectorAngles[detectorNumber][4]; // roll

	G4double x = 0;
	G4double y = 0;
	G4double z = position;

	G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
	//rotate->rotateY(M_PI);
	rotate->rotateZ(90.0*deg);
	rotate->rotateY(M_PI+beta);
	rotate->rotateZ(gamma);

	G4ThreeVector move(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,y,z,theta));

	fAssemblyWhite->MakeImprint(expHallLog, move, rotate, copyNumber);
	fAssemblyWhiteScintillator->MakeImprint(expHallLog, move, rotate, copyNumber);

	return 1;
}


G4int DetectionSystemDescant::BuildCanVolume() {
	G4ThreeVector moveCut, move, direction;
	G4RotationMatrix* rotateCut;
	G4RotationMatrix* rotate;
	G4double zPosition;

	G4Material* canG4material = G4Material::GetMaterial(fCanMaterial);
	if(!canG4material) {
		G4cout<<" ----> Material "<<fCanMaterial<<" not found, cannot build! "<<G4endl;
		return 0;
	}

	G4Material* leadG4material = G4Material::GetMaterial(fLeadMaterial);
	if(!leadG4material) {
		G4cout<<" ----> Material "<<fLeadMaterial<<" not found, cannot build! "<<G4endl;
		return 0;
	}

	G4Material* quartzG4material = G4Material::GetMaterial(fQuartzMaterial);
	if(!quartzG4material) {
		G4cout<<" ----> Material "<<fQuartzMaterial<<" not found, cannot build! "<<G4endl;
		return 0;
	}

	// for cutting PMT window
	G4Tubs * windowCut5inch = new G4Tubs("windowCut5inch", fInnerRadiusWindow, fOuterRadiusWindow5inch, fHalfLengthZWindowCut, fStartPhi, fEndPhi);
	G4Tubs * windowCut3inch = new G4Tubs("windowCut3inch", fInnerRadiusWindow, fOuterRadiusWindow3inch, fHalfLengthZWindowCut, fStartPhi, fEndPhi);

	// quartz windows
	//G4VisAttributes* quartzVisAtt = new G4VisAttributes(fMagentaColour);
	G4VisAttributes* quartzVisAtt = new G4VisAttributes(fBlackColour);
	quartzVisAtt->SetForceWireframe(true);
	quartzVisAtt->SetVisibility(true);

	G4Tubs * quartzWindow5inch = new G4Tubs("window5inch", fInnerRadiusWindow, fOuterRadiusWindow5inch, fOpticalWindowHalfThickness, fStartPhi, fEndPhi);
	G4Tubs * quartzWindow3inch = new G4Tubs("window3inch", fInnerRadiusWindow, fOuterRadiusWindow3inch, fOpticalWindowHalfThickness, fStartPhi, fEndPhi);

	if(fQuartzWindow5inchLog == nullptr) {
		fQuartzWindow5inchLog = new G4LogicalVolume(quartzWindow5inch, quartzG4material, "quartzWindow5inchLog", 0, 0, 0);
		fQuartzWindow5inchLog->SetVisAttributes(quartzVisAtt);
	}
	if(fQuartzWindow3inchLog == nullptr) {
		fQuartzWindow3inchLog = new G4LogicalVolume(quartzWindow3inch, quartzG4material, "quartzWindow3inchLog", 0, 0, 0);
		fQuartzWindow3inchLog->SetVisAttributes(quartzVisAtt);
	}


	// BLUE
	// Set visualization attributes
	G4VisAttributes* blueVisAtt = new G4VisAttributes(fBlueColour);
	blueVisAtt->SetVisibility(true);

	G4SubtractionSolid* blueVolume1      = CanVolume(false, fCanLength + fCanBackThickness, fBlueDetector, fBluePhi);
	G4SubtractionSolid* blueVolume1_cut  = CanVolume(true, fCanLength - fCanThickness, fBlueDetector, fBluePhi);

	moveCut = G4ThreeVector(0.0*mm,0.0*mm, (fCanBackThickness)/2.0 - (fCanThickness)/2.0);
	rotateCut = new G4RotationMatrix;
	G4SubtractionSolid* blueVolume2 = new G4SubtractionSolid("blueVolume2", blueVolume1, blueVolume1_cut, rotateCut, moveCut);

	moveCut = G4ThreeVector(fBluePMTOffsetX, 0.0 , -1.0*(fCanLength));
	rotateCut = new G4RotationMatrix;  
	G4SubtractionSolid * blueVolume3 = new G4SubtractionSolid("blueVolume3", blueVolume2, windowCut5inch, rotateCut, moveCut);

	// Define rotation and movement objects
	direction   = G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength + fCanBackThickness)/2.0 ;
	move          = zPosition * direction;
	rotate = new G4RotationMatrix;

	// logical volume blue
	if(fBlueVolumeLog == nullptr) {
		fBlueVolumeLog = new G4LogicalVolume(blueVolume3, canG4material, "descantBlueVolumeLog", 0, 0, 0);
		fBlueVolumeLog->SetVisAttributes(blueVisAtt);
	}
	fAssemblyBlue->AddPlacedVolume(fBlueVolumeLog, move, rotate);

	// quartz window blue
	move = G4ThreeVector(fBluePMTOffsetX, 0.0, -1.0*(fOpticalWindowHalfThickness+fCanLength));
	rotate = new G4RotationMatrix;
	fAssemblyBlue->AddPlacedVolume(fQuartzWindow5inchLog, move, rotate);

	// BLUE LEAD
	// Set visualization attributes
	G4VisAttributes* blueLeadVisAtt = new G4VisAttributes(fBlueColour);
	blueLeadVisAtt->SetVisibility(true);

	G4SubtractionSolid* blueLeadVolume1 = CanVolume(false, fLeadShieldThickness, fBlueDetector, fBluePhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = fLeadShieldThickness/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fBlueLeadVolumeLog == nullptr) {
		fBlueLeadVolumeLog = new G4LogicalVolume(blueLeadVolume1, leadG4material, "descantBlueLeadVolumeLog", 0, 0, 0);
		fBlueLeadVolumeLog->SetVisAttributes(blueLeadVisAtt);
	}

	if(fIncludeLead)
		fAssemblyBlue->AddPlacedVolume(fBlueLeadVolumeLog, move, rotate);

	// GREEN
	// Set visualization attributes
	G4VisAttributes* greenVisAtt = new G4VisAttributes(fGreenColour);
	greenVisAtt->SetVisibility(true);

	G4SubtractionSolid* greenVolume1      = CanVolume(false, fCanLength + fCanBackThickness, fGreenDetector, fGreenPhi);
	G4SubtractionSolid* greenVolume1_cut  = CanVolume(true, fCanLength -  fCanThickness, fGreenDetector, fGreenPhi);

	moveCut = G4ThreeVector(0.0*mm,0.0*mm, (fCanBackThickness)/2.0 - (fCanThickness)/2.0);
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid* greenVolume2 = new G4SubtractionSolid("greenVolume2", greenVolume1, greenVolume1_cut, rotateCut, moveCut);

	moveCut = G4ThreeVector(-1.0*fYellowGreenPMTOffsetX, fYellowGreenPMTOffsetY , -1.0*(fCanLength));
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid * greenVolume3 = new G4SubtractionSolid("greenVolume3", greenVolume2, windowCut3inch, rotateCut, moveCut);

	// Define rotation and movement objects
	direction   = G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength + fCanBackThickness)/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fGreenVolumeLog == nullptr) {
		fGreenVolumeLog = new G4LogicalVolume(greenVolume3, canG4material, "descantGreenVolumeLog", 0, 0, 0);
		fGreenVolumeLog->SetVisAttributes(greenVisAtt);
	}
	fAssemblyGreen->AddPlacedVolume(fGreenVolumeLog, move, rotate);

	// quartz window green
	move = G4ThreeVector(-1.0*fYellowGreenPMTOffsetX, fYellowGreenPMTOffsetY, -1.0*(fOpticalWindowHalfThickness+fCanLength));
	rotate = new G4RotationMatrix;
	fAssemblyGreen->AddPlacedVolume(fQuartzWindow3inchLog, move, rotate);

	// GREEN LEAD
	// Set visualization attributes
	G4VisAttributes* greenLeadVisAtt = new G4VisAttributes(fGreenColour);
	greenLeadVisAtt->SetVisibility(true);

	G4SubtractionSolid* greenLeadVolume1 = CanVolume(false, fLeadShieldThickness, fGreenDetector, fGreenPhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = fLeadShieldThickness/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fGreenLeadVolumeLog == nullptr) {
		fGreenLeadVolumeLog = new G4LogicalVolume(greenLeadVolume1, leadG4material, "descantGreenLeadVolumeLog", 0, 0, 0);
		fGreenLeadVolumeLog->SetVisAttributes(greenLeadVisAtt);
	}

	if(fIncludeLead)
		fAssemblyGreen->AddPlacedVolume(fGreenLeadVolumeLog, move, rotate);

	// RED
	// Set visualization attributes
	G4VisAttributes* redVisAtt = new G4VisAttributes(fRedColour);
	redVisAtt->SetVisibility(true);

	G4SubtractionSolid* redVolume1      = CanVolume(false, fCanLength + fCanBackThickness, fRedDetector, fRedPhi);
	G4SubtractionSolid* redVolume1_cut  = CanVolume(true, fCanLength -  fCanThickness, fRedDetector, fRedPhi);

	moveCut = G4ThreeVector(0.0*mm,0.0*mm, (fCanBackThickness)/2.0 - (fCanThickness)/2.0);
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid* redVolume2 = new G4SubtractionSolid("redVolume2", redVolume1, redVolume1_cut, rotateCut, moveCut);

	moveCut = G4ThreeVector(-1.0*fRedPMTOffsetX , 0.0 , -1.0*(fCanLength));
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid * redVolume3 = new G4SubtractionSolid("redVolume3", redVolume2, windowCut5inch, rotateCut, moveCut);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength + fCanBackThickness)/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fRedVolumeLog == nullptr) {
		fRedVolumeLog = new G4LogicalVolume(redVolume3, canG4material, "descantRedVolumeLog", 0, 0, 0);
		fRedVolumeLog->SetVisAttributes(redVisAtt);
	}
	fAssemblyRed->AddPlacedVolume(fRedVolumeLog, move, rotate);

	move = G4ThreeVector(-1.0*fRedPMTOffsetX, 0.0, -1.0*(fOpticalWindowHalfThickness+fCanLength));
	rotate = new G4RotationMatrix;
	fAssemblyRed->AddPlacedVolume(fQuartzWindow5inchLog, move, rotate);

	// RED LEAD
	// Set visualization attributes
	G4VisAttributes* redLeadVisAtt = new G4VisAttributes(fRedColour);
	redLeadVisAtt->SetVisibility(true);

	G4SubtractionSolid* redLeadVolume1 = CanVolume(false, fLeadShieldThickness, fRedDetector, fRedPhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = fLeadShieldThickness/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fRedLeadVolumeLog == nullptr) {
		fRedLeadVolumeLog = new G4LogicalVolume(redLeadVolume1, leadG4material, "descantRedLeadVolumeLog", 0, 0, 0);
		fRedLeadVolumeLog->SetVisAttributes(redLeadVisAtt);
	}

	if(fIncludeLead)
		fAssemblyRed->AddPlacedVolume(fRedLeadVolumeLog, move, rotate);

	// WHITE
	// Set visualization attributes
	G4VisAttributes* whiteVisAtt = new G4VisAttributes(fWhiteColour);
	whiteVisAtt->SetVisibility(true);

	G4SubtractionSolid* whiteVolume1      = CanVolume(false, fCanLength + fCanBackThickness, fWhiteDetector, fWhitePhi);
	G4SubtractionSolid* whiteVolume1_cut  = CanVolume(true, fCanLength - fCanThickness, fWhiteDetector, fWhitePhi);

	moveCut = G4ThreeVector(0.0*mm,0.0*mm, (fCanBackThickness)/2.0 - (fCanThickness)/2.0);
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid* whiteVolume2 = new G4SubtractionSolid("whiteVolume2", whiteVolume1, whiteVolume1_cut, rotateCut, moveCut);

	moveCut = G4ThreeVector(0.0, fWhitePMTOffsetY, -1.0*(fCanLength));
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid * whiteVolume3 = new G4SubtractionSolid("whiteVolume3", whiteVolume2, windowCut5inch, rotateCut, moveCut);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength + fCanBackThickness)/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fWhiteVolumeLog == nullptr) {
		fWhiteVolumeLog = new G4LogicalVolume(whiteVolume3, canG4material, "descantWhiteVolumeLog", 0, 0, 0);
		fWhiteVolumeLog->SetVisAttributes(whiteVisAtt);

	}
	fAssemblyWhite->AddPlacedVolume(fWhiteVolumeLog, move, rotate);

	// QUARTZ WINDOW
	move = G4ThreeVector(0.0, fWhitePMTOffsetY, -1.0*(fOpticalWindowHalfThickness+fCanLength));
	rotate = new G4RotationMatrix;
	fAssemblyWhite->AddPlacedVolume(fQuartzWindow5inchLog, move, rotate);

	// WHITE LEAD
	// Set visualization attributes
	G4VisAttributes* whiteLeadVisAtt = new G4VisAttributes(fWhiteColour);
	whiteLeadVisAtt->SetVisibility(true);

	G4SubtractionSolid* whiteLeadVolume1 = CanVolume(false, fLeadShieldThickness, fWhiteDetector, fWhitePhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = fLeadShieldThickness/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fWhiteLeadVolumeLog == nullptr) {
		fWhiteLeadVolumeLog = new G4LogicalVolume(whiteLeadVolume1, leadG4material, "descantWhiteLeadVolumeLog", 0, 0, 0);
		fWhiteLeadVolumeLog->SetVisAttributes(whiteLeadVisAtt);
	}

	if(fIncludeLead)
		fAssemblyWhite->AddPlacedVolume(fWhiteLeadVolumeLog, move, rotate);

	// YELLOW
	// Set visualization attributes
	G4VisAttributes* yellowVisAtt = new G4VisAttributes(fYellowColour);
	yellowVisAtt->SetVisibility(true);

	G4SubtractionSolid* yellowVolume1      = CanVolume(false, fCanLength + fCanBackThickness, fYellowDetector, fYellowPhi);
	G4SubtractionSolid* yellowVolume1_cut  = CanVolume(true, fCanLength -  fCanThickness, fYellowDetector, fYellowPhi);

	moveCut = G4ThreeVector(0.0*mm,0.0*mm, (fCanBackThickness)/2.0 - (fCanThickness)/2.0);
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid* yellowVolume2 = new G4SubtractionSolid("yellowVolume2", yellowVolume1, yellowVolume1_cut, rotateCut, moveCut);

	moveCut = G4ThreeVector(1.0*fYellowGreenPMTOffsetX, fYellowGreenPMTOffsetY , -1.0*(fCanLength));
	rotateCut = new G4RotationMatrix;

	G4SubtractionSolid * yellowVolume3 = new G4SubtractionSolid("yellowVolume3", yellowVolume2, windowCut3inch, rotateCut, moveCut);

	// Define rotation and movement objects
	direction   = G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength + fCanBackThickness)/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fYellowVolumeLog == nullptr) {
		fYellowVolumeLog = new G4LogicalVolume(yellowVolume3, canG4material, "descantYellowVolumeLog", 0, 0, 0);
		fYellowVolumeLog->SetVisAttributes(yellowVisAtt);
	}
	fAssemblyYellow->AddPlacedVolume(fYellowVolumeLog, move, rotate);

	// QUARTZ WINDOW
	move = G4ThreeVector(1.0*fYellowGreenPMTOffsetX, fYellowGreenPMTOffsetY, -1.0*(fOpticalWindowHalfThickness+fCanLength));
	rotate = new G4RotationMatrix;
	fAssemblyYellow->AddPlacedVolume(fQuartzWindow3inchLog, move, rotate);

	// YELLOW LEAD
	// Set visualization attributes
	G4VisAttributes* yellowLeadVisAtt = new G4VisAttributes(fYellowColour);
	yellowLeadVisAtt->SetVisibility(true);

	G4SubtractionSolid* yellowLeadVolume1 = CanVolume(false, fLeadShieldThickness, fYellowDetector, fYellowPhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = fLeadShieldThickness/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fYellowLeadVolumeLog == nullptr) {
		fYellowLeadVolumeLog = new G4LogicalVolume(yellowLeadVolume1, leadG4material, "descantYellowLeadVolumeLog", 0, 0, 0);
		fYellowLeadVolumeLog->SetVisAttributes(yellowLeadVisAtt);
	}

	if(fIncludeLead)
		fAssemblyYellow->AddPlacedVolume(fYellowLeadVolumeLog, move, rotate);

	return 1;
}

G4int DetectionSystemDescant::BuildDetectorVolume()
{
	G4ThreeVector move, direction;
	G4RotationMatrix* rotate;
	G4double zPosition;

	G4Material* material = G4Material::GetMaterial(fLiquidMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fLiquidMaterial<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	// BLUE
	// Set visualization attributes
	G4VisAttributes* blueVisAtt = new G4VisAttributes(fLiquidColour);
	blueVisAtt->SetVisibility(true);

	G4SubtractionSolid* blueVolume1_cut  = CanVolume(true, fCanLength-fCanThickness, fBlueDetector, fBluePhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength+fCanThickness)/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fBlueScintillatorVolumeLog == nullptr)
	{
		fBlueScintillatorVolumeLog = new G4LogicalVolume(blueVolume1_cut, material, "blueScintillatorVolumeLog", 0, 0, 0);
		fBlueScintillatorVolumeLog->SetVisAttributes(blueVisAtt);
	}

	fAssemblyBlueScintillator->AddPlacedVolume(fBlueScintillatorVolumeLog, move, rotate);

	//GREEN
	// Set visualization attributes
	G4VisAttributes* greenVisAtt = new G4VisAttributes(fLiquidColour);
	greenVisAtt->SetVisibility(true);

	G4SubtractionSolid* greenVolume1_cut  = CanVolume(true, fCanLength-fCanThickness, fGreenDetector, fGreenPhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength+fCanThickness)/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fGreenScintillatorVolumeLog == nullptr)
	{
		fGreenScintillatorVolumeLog = new G4LogicalVolume(greenVolume1_cut, material, "greenScintillatorVolumeLog", 0, 0, 0);
		fGreenScintillatorVolumeLog->SetVisAttributes(greenVisAtt);
	}

	fAssemblyGreenScintillator->AddPlacedVolume(fGreenScintillatorVolumeLog, move, rotate);

	// RED
	// Set visualization attributes
	G4VisAttributes* redVisAtt = new G4VisAttributes(fLiquidColour);
	redVisAtt->SetVisibility(true);

	G4SubtractionSolid* redVolume1_cut  = CanVolume(true, fCanLength-fCanThickness, fRedDetector, fRedPhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength+fCanThickness)/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fRedScintillatorVolumeLog == nullptr)
	{
		fRedScintillatorVolumeLog = new G4LogicalVolume(redVolume1_cut, material, "redScintillatorVolumeLog", 0, 0, 0);
		fRedScintillatorVolumeLog->SetVisAttributes(redVisAtt);
	}

	fAssemblyRedScintillator->AddPlacedVolume(fRedScintillatorVolumeLog, move, rotate);

	// WHITE
	// Set visualization attributes
	G4VisAttributes* whiteVisAtt = new G4VisAttributes(fLiquidColour);
	whiteVisAtt->SetVisibility(true);

	G4SubtractionSolid* whiteVolume1_cut  = CanVolume(true, fCanLength-fCanThickness, fWhiteDetector, fWhitePhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength+fCanThickness)/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fWhiteScintillatorVolumeLog == nullptr)
	{
		fWhiteScintillatorVolumeLog = new G4LogicalVolume(whiteVolume1_cut, material, "whiteScintillatorVolumeLog", 0, 0, 0);
		fWhiteScintillatorVolumeLog->SetVisAttributes(whiteVisAtt);
	}

	fAssemblyWhiteScintillator->AddPlacedVolume(fWhiteScintillatorVolumeLog, move, rotate);

	// YELLOW
	// Set visualization attributes
	G4VisAttributes* yellowVisAtt = new G4VisAttributes(fLiquidColour);
	yellowVisAtt->SetVisibility(true);

	G4SubtractionSolid* yellowVolume1_cut  = CanVolume(true, fCanLength-fCanThickness, fYellowDetector, fYellowPhi);

	// Define rotation and movement objects
	direction 	= G4ThreeVector(0,0,1);
	zPosition    = -1.0*(fCanLength+fCanThickness)/2.0 ;
	move          = zPosition * direction;

	rotate = new G4RotationMatrix;

	//logical volume
	if(fYellowScintillatorVolumeLog == nullptr)
	{
		fYellowScintillatorVolumeLog = new G4LogicalVolume(yellowVolume1_cut, material, "yellowScintillatorVolumeLog", 0, 0, 0);
		fYellowScintillatorVolumeLog->SetVisAttributes(yellowVisAtt);
	}

	fAssemblyYellowScintillator->AddPlacedVolume(fYellowScintillatorVolumeLog, move, rotate);

	return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
// This function constructs a generic 6 sided shape defined by 12 input points, and the cut phi angles.
G4SubtractionSolid* DetectionSystemDescant::CanVolume(G4bool insideVol, G4double volumeLength, G4double detector[12][3], G4double detectorPhi[6])
{
	G4int idx1;
	G4int idx2;

	G4ThreeVector frontP1;
	G4ThreeVector frontP2;

	G4double startPhi      = 0.0*deg;
	G4double endPhi        = 360.0*deg;
	G4double innerRadius;
	G4double outerRadius;
	G4double halfLengthZ;

	innerRadius = 0.0*mm;
	outerRadius = 200.0*mm;
	halfLengthZ = volumeLength/2.0;

	G4Tubs* canVolume = new G4Tubs("canVolume", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

	// Cut out a small cube out of the front of the detector to get a G4SubtractionSolid,
	// needed for CutVolumeOnFourPoints methods
	G4double cutHalfLengthX = 5.0*mm;
	G4double cutHalfLengthY = 5.0*mm;
	G4double cutHalfLengthZ = 5.0*mm;

	G4Box* cutPlate = new G4Box("cutPlate", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	G4ThreeVector moveCut = G4ThreeVector(0,outerRadius,1.0*volumeLength/2.0);
	G4RotationMatrix* rotateCut = new G4RotationMatrix;

	// snip snip
	G4SubtractionSolid* canVolumeCut1 = new G4SubtractionSolid("canVolumeCut1", canVolume, cutPlate, rotateCut, moveCut);

	// first set of points
	idx1 = 0;
	idx2 = idx1+1;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut2 = CutVolumeOnFourPoints(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut1, frontP1, frontP2);

	// next set of points
	idx1 = 1;
	idx2 = idx1+1;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut3 = CutVolumeOnFourPoints(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut2, frontP1, frontP2);

	idx1 = 2;
	idx2 = idx1+1;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut4 = CutVolumeOnFourPoints(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut3, frontP1, frontP2);

	idx1 = 3;
	idx2 = idx1+1;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut5 = CutVolumeOnFourPoints(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut4, frontP1, frontP2);

	idx1 = 4;
	idx2 = idx1+1;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut6 = CutVolumeOnFourPoints(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut5, frontP1, frontP2);

	idx1 = 5;
	idx2 = 0;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut7 = CutVolumeOnFourPoints(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut6, frontP1, frontP2);

	return canVolumeCut7;
}

// This method does the actual cutting of the surface.
G4SubtractionSolid* DetectionSystemDescant::CutVolumeOnFourPoints(G4int idx, G4bool insideVol, G4double volumeLength, G4double detectorPhi[6], G4SubtractionSolid* volume, G4ThreeVector frontP1, G4ThreeVector frontP2)
{
	G4double maxZ = 500.0*mm;
	G4double theta = 0;
	G4double phi = 0;
	G4double moveZCut = 0;

	// make the cuts excessively large.
	G4double cutHalfLengthX = 1.0*maxZ;
	G4double cutHalfLengthY = 1.0*maxZ;
	G4double cutHalfLengthZ = 2.0*maxZ;

	G4ThreeVector moveCut;
	G4RotationMatrix* rotateCut;

	// find the equation of the line defined by the two front face points.
	// then using a circle of radius r, solve for the minimum r to get the minimum theta angle to that plane!
	G4double slope  = (frontP2.y() - frontP1.y())/(frontP2.x() - frontP1.x());
	G4double intcp  = (frontP1.y()) - (slope*frontP1.x());
	G4double minx   = (-1.0*slope*intcp)/(1+slope*slope);
	G4double miny   = slope*minx + intcp;

	// make vectors, and get angle between vectors
	G4ThreeVector newV1 = G4ThreeVector(minx,miny,500.0*mm);
	G4ThreeVector newV2 = G4ThreeVector(0.0*mm,0.0*mm,500.0*mm);

	theta   = acos((newV1.dot(newV2))/(newV1.mag() * newV2.mag()));
	phi     = detectorPhi[idx];

	if(insideVol)
		moveZCut = 1.0*(volumeLength/2.0 + maxZ - (((fCanThickness)/(sin(theta))) - fCanThickness));
	else
		moveZCut = 1.0*(volumeLength/2.0 + maxZ);

	if(volumeLength == fLeadShieldThickness) // this is the lead shield, push it a bit closer so the tapers match up
		moveZCut = moveZCut - volumeLength;

	G4Box* cutPlate = new G4Box("cutPlate", cutHalfLengthX, cutHalfLengthY, cutHalfLengthZ);

	moveCut = G4ThreeVector(cutHalfLengthX/cos(theta),0.0*mm,moveZCut);
	moveCut.rotateZ(1.0*phi);

	// rotate our cutting block to the "chopping" angle
	rotateCut = new G4RotationMatrix;
	rotateCut->rotateZ(-1.0*phi);
	rotateCut->rotateY(1.0*theta);

	// snip snip
	G4SubtractionSolid* canVolumeCut = new G4SubtractionSolid("canVolumeCut", volume, cutPlate, rotateCut, moveCut);

	return canVolumeCut;
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
