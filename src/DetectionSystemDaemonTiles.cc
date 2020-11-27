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

#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4OpticalSurface.hh"
#include "G4OpticalPhysics.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemDaemonTiles.hh"

#include "G4SystemOfUnits.hh"

#include <string>

//DetectionSystemDaemonTiles::DetectionSystemDaemonTiles(G4bool leadShield) :
DetectionSystemDaemonTiles::DetectionSystemDaemonTiles(G4double thickness, G4int material) :
	// LogicalVolumes
	fBluePlasticLog(0),
	fBluePMTLog(0),
	fBlueWrapLog(0),
	fWhitePlasticLog(0),
	fWhitePMTLog(0),
	fWhiteWrapLog(0),
	fRedPlasticLog(0),
	fRedPMTLog(0),
	fRedWrapLog(0),
	fGreenPlasticLog(0),
	fGreenPMTLog(0),
	fGreenWrapLog(0),
	fYellowPlasticLog(0),
	fYellowPMTLog(0),
	fYellowWrapLog(0)
{
	//Logical Volumes
	fBluePlasticLog = NULL;
	fBluePMTLog = NULL;
	fBlueWrapLog = NULL;
	fWhitePlasticLog = NULL;
	fWhitePMTLog = NULL;
	fWhiteWrapLog = NULL;
	fRedPlasticLog = NULL;
	fRedPMTLog = NULL;
	fRedWrapLog = NULL;
	fGreenPlasticLog = NULL;
	fGreenPMTLog = NULL;
	fGreenWrapLog = NULL;
	fYellowPlasticLog = NULL;
	fYellowPMTLog = NULL;
	fYellowWrapLog = NULL;
	//Segment Arrays
	fBluePlasticLogArray.resize(4, NULL);
	fBluePMTLogArray.resize(4, NULL);
	fWhitePlasticLogArray.resize(4, NULL);
	fWhitePMTLogArray.resize(4, NULL);
	fRedPlasticLogArray.resize(4, NULL);
	fRedPMTLogArray.resize(4, NULL);
	fGreenPlasticLogArray.resize(4, NULL);
	fGreenPMTLogArray.resize(4, NULL);
	fYellowPlasticLogArray.resize(4, NULL);
	fYellowPMTLogArray.resize(4, NULL);

	// can properties
	fWrapThickness = 0.5 * mm; 
	fPlasticThickness   = thickness; //cm
	fAirGap = 0.001*mm;

	//Materials
	fWrapMaterial = "Teflon";
	fPMTMaterial = "G4_SILICON_DIOXIDE";

	if(material == 1)  fPlasticMaterial = "BC404";
	else if (material == 2) fPlasticMaterial = "BC408";
	else if (material == 3) fPlasticMaterial = "deuterium";
	else if (material == 4) fPlasticMaterial = "Hydrogen";
	else if (material == 5) fPlasticMaterial = "Carbon";
	else if (material == 6) fPlasticMaterial = "Deuterated Scintillator";
	else if (material == 7) fPlasticMaterial = "Dense";
	else if (material == 8) fPlasticMaterial = "BC537";
	else G4cout<< "Material Unknown" << G4endl;

	// Check surfaces to determine any problematic overlaps. Turn this on to have Geant4 check the surfaces.
	// Do not leave this on, it will slow the DetectorConstruction process!
	// This was last check on May 29th, 2015. - GOOD!, but...
	// The green, yellow and blue detectors had small overlaps, this is probably due to taking the pure mathematical model
	// of DESCANT, and then Saint Gobain rounding the result to generate the data files (below).
	fSurfCheck = false;
	// set DAEMON Wrap and segment on/off
	fAddSegment = true;
	fAddWrap = true;




	// The face of the detector is 50 cm from the origin, this does NOT include the lead shield.
	fRadialDistance         = 50*cm;
	//fRadialDistance         = 49.5*cm;


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
	fBronzeColour   = G4Colour(0.8,0.5,0.2);
	fSilverColour   = G4Colour(0.75,0.75,0.75);

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


DetectionSystemDaemonTiles::~DetectionSystemDaemonTiles() {
	// LogicalVolumes
	delete fBluePlasticLog;
	delete fBluePMTLog;
	delete fBlueWrapLog;
	delete fWhitePlasticLog;
	delete fWhitePMTLog;
	delete fWhiteWrapLog;
	delete fRedPlasticLog;
	delete fRedPMTLog;
	delete fRedWrapLog;
	delete fGreenPlasticLog;
	delete fGreenPMTLog;
	delete fGreenWrapLog;
	delete fYellowPlasticLog;
	delete fYellowPMTLog;
	delete fYellowWrapLog;

	for ( int i = 0 ; i<4; i++){
		delete fBluePlasticLogArray[i];
		delete fBluePMTLogArray[i];
		delete fWhitePlasticLogArray[i];
		delete fWhitePMTLogArray[i];
		delete fRedPlasticLogArray[i];
		delete fRedPMTLogArray[i];
		delete fGreenPlasticLogArray[i];
		delete fGreenPMTLogArray[i];
		delete fYellowPlasticLogArray[i];
		delete fYellowPMTLogArray[i];
	}

}


G4int DetectionSystemDaemonTiles::Build() {
	// Build assembly volumes
	// assemblyBlue contains non-sensitive volumes of the blue detector
	// assemblyBlueScintillator contains only sensitive volumes of the blue detector
	G4AssemblyVolume* myAssemblyBlue = new G4AssemblyVolume();
	fAssemblyBlue = myAssemblyBlue;

	G4AssemblyVolume* myAssemblyGreen = new G4AssemblyVolume();
	fAssemblyGreen = myAssemblyGreen;

	G4AssemblyVolume* myAssemblyRed = new G4AssemblyVolume();
	fAssemblyRed = myAssemblyRed;

	G4AssemblyVolume* myAssemblyWhite = new G4AssemblyVolume();
	fAssemblyWhite = myAssemblyWhite;

	G4AssemblyVolume* myAssemblyYellow = new G4AssemblyVolume();
	fAssemblyYellow = myAssemblyYellow;

	G4cout<<"Building DAEMON Tiles..."<<G4endl;

	G4cout<<"BuildCanVolume"<<G4endl;
	BuildCanVolume();

	G4cout<<"BuildDetectorVolume"<<G4endl;
	//BuildDetectorVolume();

	return 1;
}

G4int DetectionSystemDaemonTiles::PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber) {
	G4double position = fRadialDistance;

	G4int idx;
	G4double x = 0;
	G4double y = 0;
	G4double z = position;
	//For location of pmt's
	G4double pmtSizeZ = 0.5*cm;
	G4double zpmt  = fPlasticThickness/2. + fWrapThickness;
	G4double xypmt = 2.5*cm;
	G4double x2=0, y2=0, z2=0, R=0;
	G4ThreeVector move2;

	G4ThreeVector move;
	G4RotationMatrix* rotate;

	G4cout << "Blue detectors" << G4endl;	
	G4cout << "DetNum (i), x, y, z" << G4endl;
	// Place Blue Detector
	for(G4int i=0; i<(detectorNumber-55); i++) {
		move = G4ThreeVector(x,y,z);
		idx = i;
		move.rotateZ(fBlueAlphaBetaGamma[idx][2]);
		move.rotateY(fBlueAlphaBetaGamma[idx][1]);
		move.rotateZ(fBlueAlphaBetaGamma[idx][0]);
		rotate = new G4RotationMatrix;
		rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
		rotate->rotateZ(fBlueAlphaBetaGamma[idx][2]);
		rotate->rotateY(fBlueAlphaBetaGamma[idx][1]);
		rotate->rotateZ(fBlueAlphaBetaGamma[idx][0]);
		fAssemblyBlue->MakeImprint(expHallLog, move, rotate, i, fSurfCheck);

		//To get coordinates of pmt
		//Unsegmented
		if (fAddSegment == false) {
			x2 = 0;
			y2 = 0;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][2]);
			move2.rotateY(fBlueAlphaBetaGamma[idx][1]);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			//G4cout << "{" << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
		}
		//Segmented
		if (fAddSegment == true) {
			x2 = xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][2]);
			move2.rotateY(fBlueAlphaBetaGamma[idx][1]);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
	/*		x2 = -xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][2]);
			move2.rotateY(fBlueAlphaBetaGamma[idx][1]);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][2]);
			move2.rotateY(fBlueAlphaBetaGamma[idx][1]);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = -xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][2]);
			move2.rotateY(fBlueAlphaBetaGamma[idx][1]);
			move2.rotateZ(fBlueAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
		*/}

	}
	
	G4cout << "White detectors" << G4endl;	
	G4cout << "DetNum (i), x, y, z" << G4endl;
	// Place White Detector
	for(G4int i=40; i<(detectorNumber-10); i++) {
		idx = i-40;
		move = G4ThreeVector(x,y,z);
		move.rotateZ(fWhiteAlphaBetaGamma[idx][2]);
		move.rotateY(fWhiteAlphaBetaGamma[idx][1]);
		move.rotateZ(fWhiteAlphaBetaGamma[idx][0]);
		rotate = new G4RotationMatrix;
		rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
		rotate->rotateZ(fWhiteAlphaBetaGamma[idx][2]);
		rotate->rotateY(fWhiteAlphaBetaGamma[idx][1]);
		rotate->rotateZ(fWhiteAlphaBetaGamma[idx][0]);
		fAssemblyWhite->MakeImprint(expHallLog, move, rotate, i, fSurfCheck);
		//To get coordinates of pmt
		//Unsegmented
		if (fAddSegment == false) {
			x2 = 0;
			y2 = 0;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][2]);
			move2.rotateY(fWhiteAlphaBetaGamma[idx][1]);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			//G4cout << "{" << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
		}
		//Segmented
		if (fAddSegment == true) {
			x2 = xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][2]);
			move2.rotateY(fWhiteAlphaBetaGamma[idx][1]);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
	/*		x2 = -xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][2]);
			move2.rotateY(fWhiteAlphaBetaGamma[idx][1]);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][2]);
			move2.rotateY(fWhiteAlphaBetaGamma[idx][1]);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = -xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][2]);
			move2.rotateY(fWhiteAlphaBetaGamma[idx][1]);
			move2.rotateZ(fWhiteAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
		*/}
	}
	G4cout << "Red detectors" << G4endl;	
	G4cout << "DetNum (i), x, y, z" << G4endl;
	// Place Red Detector
	for(G4int i=25; i<(detectorNumber-30); i++) {
		idx = i-25;
		move = G4ThreeVector(x,y,z);
		move.rotateZ(fRedAlphaBetaGamma[idx][2]);
		move.rotateY(fRedAlphaBetaGamma[idx][1]);
		move.rotateZ(fRedAlphaBetaGamma[idx][0]);
		rotate = new G4RotationMatrix;
		rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
		rotate->rotateZ(fRedAlphaBetaGamma[idx][2]);
		rotate->rotateY(fRedAlphaBetaGamma[idx][1]);
		rotate->rotateZ(fRedAlphaBetaGamma[idx][0]);
		fAssemblyRed->MakeImprint(expHallLog, move, rotate, i, fSurfCheck);
		//To get coordinates of pmt
		//Unsegmented
		if (fAddSegment == false) {
			x2 = 0;
			y2 = 0;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fRedAlphaBetaGamma[idx][2]);
			move2.rotateY(fRedAlphaBetaGamma[idx][1]);
			move2.rotateZ(fRedAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			//G4cout << "{" << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
		}
		//Segmented
		if (fAddSegment == true) {
			x2 = xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fRedAlphaBetaGamma[idx][2]);
			move2.rotateY(fRedAlphaBetaGamma[idx][1]);
			move2.rotateZ(fRedAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
	/*		x2 = -xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fRedAlphaBetaGamma[idx][2]);
			move2.rotateY(fRedAlphaBetaGamma[idx][1]);
			move2.rotateZ(fRedAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fRedAlphaBetaGamma[idx][2]);
			move2.rotateY(fRedAlphaBetaGamma[idx][1]);
			move2.rotateZ(fRedAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = -xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fRedAlphaBetaGamma[idx][2]);
			move2.rotateY(fRedAlphaBetaGamma[idx][1]);
			move2.rotateZ(fRedAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
		*/}
	}

	G4cout << "Green detectors" << G4endl;	
	G4cout << "DetNum (i), x, y, z" << G4endl;
	// Place Green Detector
	for(G4int i=15; i<(detectorNumber-45); i++) {
		idx = i-15;
		move = G4ThreeVector(x,y,z);
		move.rotateZ(fGreenAlphaBetaGamma[idx][2]);
		move.rotateY(fGreenAlphaBetaGamma[idx][1]);
		move.rotateZ(fGreenAlphaBetaGamma[idx][0]);
		rotate = new G4RotationMatrix;
		rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
		rotate->rotateZ(fGreenAlphaBetaGamma[idx][2]);
		rotate->rotateY(fGreenAlphaBetaGamma[idx][1]);
		rotate->rotateZ(fGreenAlphaBetaGamma[idx][0]);
		fAssemblyGreen->MakeImprint(expHallLog, move, rotate, i, fSurfCheck);
		//To get coordinates of pmt
		//Unsegmented
		if (fAddSegment == false) {
			x2 = 0;
			y2 = 0;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][2]);
			move2.rotateY(fGreenAlphaBetaGamma[idx][1]);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			//G4cout << "{" << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
		}
		//Segmented
		if (fAddSegment == true) {
			x2 = xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][2]);
			move2.rotateY(fGreenAlphaBetaGamma[idx][1]);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
	/*		x2 = -xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][2]);
			move2.rotateY(fGreenAlphaBetaGamma[idx][1]);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][2]);
			move2.rotateY(fGreenAlphaBetaGamma[idx][1]);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = -xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][2]);
			move2.rotateY(fGreenAlphaBetaGamma[idx][1]);
			move2.rotateZ(fGreenAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
		*/}
	}

	G4cout << "Yellow detectors" << G4endl;	
	G4cout << "DetNum (i), x, y, z" << G4endl;
	// Place Yellow Detector
	for(G4int i=60; i<(detectorNumber); i++) {
		idx = i-60;
		move = G4ThreeVector(x,y,z);
		move.rotateZ(fYellowAlphaBetaGamma[idx][2]);
		move.rotateY(fYellowAlphaBetaGamma[idx][1]);
		move.rotateZ(fYellowAlphaBetaGamma[idx][0]);
		rotate = new G4RotationMatrix;
		rotate->rotateY(M_PI); // flip the detector so that the face is pointing upstream.
		rotate->rotateZ(fYellowAlphaBetaGamma[idx][2]);
		rotate->rotateY(fYellowAlphaBetaGamma[idx][1]);
		rotate->rotateZ(fYellowAlphaBetaGamma[idx][0]);
		fAssemblyYellow->MakeImprint(expHallLog, move, rotate, i, fSurfCheck);
		//To get coordinates of pmt
		//Unsegmented
		if (fAddSegment == false) {
			x2 = 0;
			y2 = 0;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][2]);
			move2.rotateY(fYellowAlphaBetaGamma[idx][1]);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			//G4cout << "{" << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;

		}
		//Segmented
		if (fAddSegment == true) {
			x2 = xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][2]);
			move2.rotateY(fYellowAlphaBetaGamma[idx][1]);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
	/*		x2 = -xypmt;
			y2 = xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][2]);
			move2.rotateY(fYellowAlphaBetaGamma[idx][1]);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Top PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][2]);
			move2.rotateY(fYellowAlphaBetaGamma[idx][1]);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
			x2 = -xypmt;
			y2 = -xypmt;
			z2 = position - zpmt;
			move2 = G4ThreeVector(x2,y2,z2);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][2]);
			move2.rotateY(fYellowAlphaBetaGamma[idx][1]);
			move2.rotateZ(fYellowAlphaBetaGamma[idx][0]);
			x2 = move2.getX();
			y2 = move2.getY();
			z2 = move2.getZ();
			R = move2.mag();
			G4cout << "Front Bottom PMT: DetNum (i), x, y, z" << G4endl;
			G4cout << "{" << R << " , " << x2 << " , " <<  y2 << " , " << z2 << "}," << G4endl;
		*/}
	}

	return 1;
}

G4int DetectionSystemDaemonTiles::BuildCanVolume() {
	G4ThreeVector moveCut, move, direction;
	G4RotationMatrix* rotate;
	G4double zPosition;


	G4Material* plasticG4material = G4Material::GetMaterial(fPlasticMaterial);
	if( !plasticG4material ) {
		G4cout << " ----> Material " << fPlasticMaterial << " not found, cannot build! " << G4endl;
		return 0;
	}
	else {
		G4cout << plasticG4material->GetName() << " is the name of the detector material" << G4endl;
	}
	G4Material* wrapG4material = G4Material::GetMaterial(fWrapMaterial);
	if( !wrapG4material ) {
		G4cout << " ----> Material " << fWrapMaterial << " not found, cannot build! " << G4endl;
		return 0;
	}
	else {
		G4cout << wrapG4material->GetName() << " is the name of the wrapping material" << G4endl;
	}
	G4Material* PMTG4material = G4Material::GetMaterial(fPMTMaterial);
	if( !PMTG4material ) {
		G4cout << " ----> Material " << fPMTMaterial << " not found, cannot build! " << G4endl;
		return 0;
	}
	else {
		G4cout << PMTG4material->GetName() << " is the name of the pmt material" << G4endl;
	}



	////////Scintillation Properties ////////  --------- Might have to be put before the material is constructed
	//Based on BC408 data
	G4MaterialPropertiesTable * scintillatorMPT = new G4MaterialPropertiesTable();

	//If no scintillation yield for p,d,t,a,C then they all default to the electron yield.
	//Have to uncomment line in ConstructOp that allows for this to work with boolean (true)
	//The following data is for BC400, very similar in properties and composition then BC408.
	const G4int num2 = 4;
	G4double e_range[num2] = {1.*keV, 0.1*MeV, 1.*MeV, 10.*MeV};
	G4double yield_e[num2] = {10., 1000., 10000., 100000.};//More realistic
	//G4double yield_e[num2] = {1000., 10000., 10000., 100000.}; //testing
	G4double yield_p[num2] = {1., 65., 1500., 45000.};
	G4double yield_d[num2] = {1., 65., 1500., 45000.};//no data provided, assume same order of magnitude as proton
	G4double yield_t[num2] = {1., 65., 1500., 45000.};//no data provided, assume same order of magnitude as proton
	G4double yield_a[num2] = {1., 20., 200., 14000.};
	G4double yield_C[num2] = {1., 10., 70., 600.};
	assert(sizeof(e_test) == sizeof(num_test));
	assert(sizeof(p_test) == sizeof(num_test));
	assert(sizeof(d_test) == sizeof(num_test));
	assert(sizeof(t_test) == sizeof(num_test));
	assert(sizeof(a_test) == sizeof(num_test));
	assert(sizeof(C_test) == sizeof(num_test));

	scintillatorMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", e_range, yield_e, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("PROTONSCINTILLATIONYIELD", e_range, yield_p, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("DEUTERONSCINTILLATIONYIELD", e_range, yield_d, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("TRITONSCINTILLATIONYIELD", e_range, yield_t, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("ALPHASCINTILLATIONYIELD", e_range, yield_a, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("IONSCINTILLATIONYIELD", e_range, yield_C, num2)->SetSpline(true);

	//scintillatorMPT->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV); //Scintillation Efficiency - characteristic light yield //10000./MeV
	///////
	//

	if (fPlasticMaterial == "BC408"){
		const G4int num = 12; //BC408
		G4cout << "BC408 and num = " << num << G4endl;
		G4double photonEnergy[num] = {1.7*eV, 2.38*eV, 2.48*eV, 2.58*eV, 2.70*eV, 2.76*eV, 2.82*eV, 2.91*eV, 2.95*eV, 3.1*eV, 3.26*eV, 3.44*eV}; //BC408 emission spectra & corresponding energies
		G4double RIndex1[num] = {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58}; //BC408
		G4double absorption[num] = {380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm}; ///light attenuation BC408
		G4double scint[num] = {3., 3., 8., 18., 43., 55., 80., 100., 80., 20., 7., 3. }; ///// Based off emission spectra for BC408

		assert(sizeof(RIndex1) == sizeof(photonEnergy));
		const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);
		assert(sizeof(absorption) == sizeof(photonEnergy));
		assert(sizeof(scint) == sizeof(photonEnergy));
		G4cout << "nEntries = " << nEntries << G4endl;

		scintillatorMPT->AddProperty("FASTCOMPONENT", photonEnergy, scint, nEntries)->SetSpline(true); // BC408 emission spectra
		scintillatorMPT->AddProperty("SLOWCOMPONENT", photonEnergy, scint, nEntries)->SetSpline(true); // BC408 emission spectra
		scintillatorMPT->AddProperty("RINDEX", photonEnergy, RIndex1, nEntries);  //refractive index can change with energy
		//note if photon is created outside of energy range it will have no index of refraction
		scintillatorMPT->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)->SetSpline(true); //absorption length doesnt change with energy - examples showing it can...
		//scintillatorMPT->AddConstProperty("ABSLENGTH", 380.*cm); //Bulk light attenuation 

		scintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns); //only one decay constant given - BC408
		scintillatorMPT->AddConstProperty("SLOWTIMECONSTANT", 2.1*ns); //only one decay constant given - BC408
		//scintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 0.00000000000000000*ns); // for testing the effective speed of light
		//scintillatorMPT->AddConstProperty("SLOWTIMECONSTANT", 0.00000000000000000*ns); // for testing the effective speed of light

		//Should these be in the physics list?
		//G4OpticalPhysics * opticalPhysics = new G4OpticalPhysics();
		//opticalPhysics->SetFiniteRiseTime(true);
		scintillatorMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.9*ns); //default rise time is 0ns, have to set manually BC408
		scintillatorMPT->AddConstProperty("SLOWSCINTILLATIONRISETIME", 0.9*ns); //default rise time is 0ns, have to set manually BC408

		//scintillatorMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.*ns); // For testing speed of light
		//scintillatorMPT->AddConstProperty("SLOWSCINTILLATIONRISETIME", 0.*ns); // For testing speed of light
	}

	if(fPlasticMaterial == "BC404"){
		const G4int num = 20; //BC404
		G4cout << "BC404 and num = " << num << G4endl;
		G4double photonEnergy[num] = {1.7*eV, 2.38*eV, 2.48*eV, 2.58*eV, 2.70*eV, 2.76*eV, 2.82*eV, 2.91*eV, 2.95*eV, 2.97*eV, 3.0*eV, 3.02*eV, 3.04*eV,  3.06*eV, 3.1*eV, 3.14*eV, 3.18*eV, 3.21*eV, 3.26*eV, 3.44*eV}; //BC404 emission spectra & corresponding energies
		G4double RIndex1[num] = {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58,  1.58, 1.58, 1.58}; //BC404
		G4double absorption[num] = {160.*cm, 160*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm}; ///light attenuation BC404
		G4double scint[num] = {0., 1., 2., 5., 13., 20., 35., 50., 55., 60., 85., 93., 100., 96., 87., 70., 38., 18., 5., 1. }; ///// Based off emission spectra for BC404

		assert(sizeof(RIndex1) == sizeof(photonEnergy));
		const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);
		assert(sizeof(absorption) == sizeof(photonEnergy));
		assert(sizeof(scint) == sizeof(photonEnergy));
		G4cout << "nEntries = " << nEntries << G4endl;

		scintillatorMPT->AddProperty("FASTCOMPONENT", photonEnergy, scint, nEntries)->SetSpline(true); // BC408 emission spectra
		scintillatorMPT->AddProperty("SLOWCOMPONENT", photonEnergy, scint, nEntries)->SetSpline(true); // BC408 emission spectra
		scintillatorMPT->AddProperty("RINDEX", photonEnergy, RIndex1, nEntries);  //refractive index can change with energy

		//note if photon is created outside of energy range it will have no index of refraction
		scintillatorMPT->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)->SetSpline(true); //absorption length doesnt change with energy - examples showing it can...
		//scintillatorMPT->AddConstProperty("ABSLENGTH", 160.*cm); //Bulk light attenuation 

		scintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 1.8*ns); //only one decay constant given - BC404
		scintillatorMPT->AddConstProperty("SLOWTIMECONSTANT", 1.8*ns); //only one decay constant given - BC404
		//scintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 0.00000000000000000*ns); // for testing the effective speed of light
		//scintillatorMPT->AddConstProperty("SLOWTIMECONSTANT", 0.00000000000000000*ns); // for testing the effective speed of light

		//Should these be in the physics list?
		//G4OpticalPhysics * opticalPhysics = new G4OpticalPhysics();
		//opticalPhysics->SetFiniteRiseTime(true);
		scintillatorMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.7*ns); //default rise time is 0ns, have to set manually BC404
		scintillatorMPT->AddConstProperty("SLOWSCINTILLATIONRISETIME", 0.7*ns); //default rise time is 0ns, have to set manually BC404
		//scintillatorMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.*ns); // For testing speed of light
		//scintillatorMPT->AddConstProperty("SLOWSCINTILLATIONRISETIME", 0.*ns); // For testing speed of light
	}

	// The number of photons produced per interaction is sampled from a Gaussian distribution with a full-width at half-maximum set to 20% of the number of produced photons. From Joeys Thesis
	scintillatorMPT->AddConstProperty("RESOLUTIONSCALE", 1.2); // broadens the statistical distribution of generated photons, sqrt(num generated)* resScale, gaussian based on SCINTILLATIONYIELD, >1 broadens, 0 no distribution. 20%
	scintillatorMPT->AddConstProperty("YIELDRATIO", 1.0); //The relative strength of the fast component as a fraction of total scintillation yield is given by the YIELDRATIO.

	//properties I may be missing: scintillation, rayleigh
	plasticG4material->SetMaterialPropertiesTable(scintillatorMPT);

	const G4int numShort =3;
	G4double photonEnergyShort[numShort] = {1.7*eV,   2.82*eV, 3.44*eV}; //BC408 emission spectra & corresponding energies
	const G4int nEntriesShort = sizeof(photonEnergyShort)/sizeof(G4double);
	//////Optical Surface - Teflon wrapping //////
	G4OpticalSurface * ScintWrapper = new G4OpticalSurface("wrapper");
	G4MaterialPropertiesTable * ScintWrapperMPT = new G4MaterialPropertiesTable();
	ScintWrapper->SetModel(unified);  // unified or glisur
	ScintWrapper->SetType(dielectric_dielectric);  // dielectric and dielectric or metal?
	// teflon wrapping on polished surface->front/back painted // Teflon should be Lambertian in air, specular in optical grease
	//polished front painted is more simplified. Only specular spike reflection
	ScintWrapper->SetFinish(polishedfrontpainted);  

	/*	//poished back painted is maybe more realistic, need to then include sigma alpha (angle of the micro facet to the average normal surface) 
		ScintWrapper->SetFinish(polishedbackpainted);	
		ScintWrapper->SetSigmaAlpha(0.1); // 0 for smooth, 1 for max roughness
		const G4int NUM =3;
		G4double pp[NUM] = {2.038*eV, 4.144*eV};
		G4double specularlobe[NUM] = {0.033, 0.033};
		G4double specularspike[NUM] = {0.9, 0.9};
		G4double backscatter[NUM] = {0.033, 0.033};
	//Diffuse lobe constant is implicit, but spec lobe, spec spike, and backscatter all need to add up to 1. //diffuse lobe constant is the probability of internal lambertian refection
	ScintWrapperMPT->AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM); //reflection probability about the normal of the micro facet
	ScintWrapperMPT->AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM); //reflection probability about average surface normal
	ScintWrapperMPT->AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM); //probability of exact back scatter based on mutiple reflections within the deep groove
	//end of polished back painted
	*/
	G4double rIndex_Teflon[numShort] = {1.35, 1.35, 1.35}; //Taken from wikipedia
	ScintWrapperMPT->AddProperty("RINDEX", photonEnergyShort, rIndex_Teflon, nEntriesShort)->SetSpline(true);  //refractive index can change with energy
	G4double reflectivity[numShort] = {0.95, 0.95, 0.95};
	//G4double reflectivity[numShort] = {0.9999, 0.9999, 0.9999};
	ScintWrapperMPT->AddProperty("REFLECTIVITY", photonEnergyShort, reflectivity, nEntriesShort)->SetSpline(true);  // light reflected / incident light.
	//G4double efficiency[numShort] = {0.95, 0.95, 0.95};
	//ScintWrapperMPT->AddProperty("EFFICIENCY",photonEnergyShort,efficiency,nEntriesShort); //This is the Quantum Efficiency of the photocathode = # electrons / # of incident photons
	ScintWrapper->SetMaterialPropertiesTable(ScintWrapperMPT);
	ScintWrapper->DumpInfo();

	//////// Quartz  ////////////
	G4MaterialPropertiesTable * QuartzMPT = new G4MaterialPropertiesTable();
	G4double rIndex_Quartz[numShort] = {1.474, 1.474, 1.474}; //Taken from Joey github
	QuartzMPT->AddProperty("RINDEX", photonEnergyShort, rIndex_Quartz, nEntriesShort)->SetSpline(true);  //refractive index can change with energy
	QuartzMPT->AddConstProperty("ABSLENGTH", 40.*cm); //from Joeys github
	PMTG4material->SetMaterialPropertiesTable(QuartzMPT);

	//Visualization properties
	//Uniform
	G4VisAttributes* PlasticVisAtt = new G4VisAttributes(fSilverColour);
	PlasticVisAtt->SetVisibility(true);
	G4VisAttributes* PMTVisAtt = new G4VisAttributes(fBronzeColour);
	PMTVisAtt->SetVisibility(true);
	G4VisAttributes* WrapVisAtt = new G4VisAttributes(fBlackColour);
	WrapVisAtt->SetVisibility(true);
	//Blue
	G4VisAttributes* bluePlasticVisAtt = new G4VisAttributes(fSilverColour);
	bluePlasticVisAtt->SetVisibility(true);
	G4VisAttributes* blueWrapVisAtt = new G4VisAttributes(fBlueColour);
	blueWrapVisAtt->SetVisibility(true);
	G4VisAttributes* bluePMTVisAtt = new G4VisAttributes(fBronzeColour);
	bluePMTVisAtt->SetVisibility(true);
	//White
	G4VisAttributes* whitePlasticVisAtt = new G4VisAttributes(fSilverColour);
	whitePlasticVisAtt->SetVisibility(true);
	G4VisAttributes* whiteWrapVisAtt = new G4VisAttributes(fWhiteColour);
	whiteWrapVisAtt->SetVisibility(true);
	G4VisAttributes* whitePMTVisAtt = new G4VisAttributes(fBronzeColour);
	whitePMTVisAtt->SetVisibility(true);
	//Red
	G4VisAttributes* redPlasticVisAtt = new G4VisAttributes(fSilverColour);
	redPlasticVisAtt->SetVisibility(true);
	G4VisAttributes* redWrapVisAtt = new G4VisAttributes(fRedColour);
	redWrapVisAtt->SetVisibility(true);
	G4VisAttributes* redPMTVisAtt = new G4VisAttributes(fBronzeColour);
	redPMTVisAtt->SetVisibility(true);
	//Green
	G4VisAttributes* greenPlasticVisAtt = new G4VisAttributes(fSilverColour);
	greenPlasticVisAtt->SetVisibility(true);
	G4VisAttributes* greenWrapVisAtt = new G4VisAttributes(fGreenColour);
	greenWrapVisAtt->SetVisibility(true);
	G4VisAttributes* greenPMTVisAtt = new G4VisAttributes(fBronzeColour);
	greenPMTVisAtt->SetVisibility(true);
	//Yellow
	G4VisAttributes* yellowPlasticVisAtt = new G4VisAttributes(fSilverColour);
	yellowPlasticVisAtt->SetVisibility(true);
	G4VisAttributes* yellowWrapVisAtt = new G4VisAttributes(fYellowColour);
	yellowWrapVisAtt->SetVisibility(true);
	G4VisAttributes* yellowPMTVisAtt = new G4VisAttributes(fBronzeColour);
	yellowPMTVisAtt->SetVisibility(true);

	//Names
	G4String nameWrapperB = "wrapperBlue";
	G4String nameWrapperW = "wrapperWhite";
	G4String nameWrapperR = "wrapperRed";
	G4String nameWrapperG = "wrapperGreen";
	G4String nameWrapperY = "wrapperYellow";

	//PMT sizes 
	G4double pmtSize, pmtSizeZ, pmtx, pmty;
	//PMT sizes Unsegmented
	if(!fAddSegment){
		pmtSize = 1.2*cm;
		pmtSizeZ = 0.5*cm;
		pmtx = 0.*cm;
		pmty = 0.*cm;
	}
	//PMT sizes Segmented
	if(fAddSegment){
		pmtSize = 0.6*cm;
		pmtSizeZ = 0.5*cm;
		pmtx = 2.5*cm;
		pmty = 2.5*cm;
	}	
	//For Segment cutting
	G4double thick = 5.*cm;
	G4double SegPosition = thick + fWrapThickness + fAirGap;
	//G4double SegPositionWrap = thick + fWrapThickness;


	// BLUE Detectors

	if(!fAddSegment){

		G4SubtractionSolid* bluePlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fBlueDetector, fBluePhi, false);
		G4SubtractionSolid* blueWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fBlueDetector, fBluePhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;
		G4SubtractionSolid* blueWrapVolume2 = new G4SubtractionSolid("blueWrapVolume2", blueWrapVolume1, bluePlasticVolume1, rotate, move);
		G4Box * bluePMT = new G4Box("bluePMT", pmtSize, pmtSize, pmtSizeZ);
		move = G4ThreeVector(0.0, 0.0, fPlasticThickness/2.); //+pmtSize?
		G4SubtractionSolid* blueWrapVolume3 = new G4SubtractionSolid("blueWrapVolume3", blueWrapVolume2, bluePMT, rotate, move);

		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;

		if(fBlueWrapLog == nullptr) {
			fBlueWrapLog = new G4LogicalVolume(blueWrapVolume3, wrapG4material, nameWrapperB, 0, 0, 0);
			fBlueWrapLog->SetVisAttributes(blueWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperB, fBlueWrapLog, ScintWrapper);
		}
		if(fBluePMTLog == nullptr) {
			fBluePMTLog = new G4LogicalVolume(bluePMT, PMTG4material, "BluePMT_1", 0, 0, 0);
			fBluePMTLog->SetVisAttributes(PMTVisAtt);
		}
		if(fBluePlasticLog == nullptr) {
			fBluePlasticLog = new G4LogicalVolume(bluePlasticVolume1, plasticG4material, "BlueTilePlasticDet_1", 0, 0, 0);
			fBluePlasticLog->SetVisAttributes(PlasticVisAtt);
		}

		fAssemblyBlue->AddPlacedVolume(fBluePlasticLog, move, rotate);
		if(fAddWrap){
			fAssemblyBlue->AddPlacedVolume(fBlueWrapLog, move, rotate);
		}
		direction 	= G4ThreeVector(0,0,1);
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		move          = zPosition * direction;
		fAssemblyBlue->AddPlacedVolume(fBluePMTLog, move, rotate);
	}

	if(fAddSegment){

		G4SubtractionSolid* bluePlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fBlueDetector, fBluePhi, false);
		G4SubtractionSolid* blueWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fBlueDetector, fBluePhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;

		G4Box * blueCut = new G4Box("blueCut", thick, thick, thick);
		//For Plastics
		G4IntersectionSolid * blueTopLeft= new G4IntersectionSolid("blueTopLeft", bluePlasticVolume1, blueCut, rotate, G4ThreeVector(SegPosition, SegPosition,0));
		G4IntersectionSolid * blueTopRight= new G4IntersectionSolid("blueTopRight", bluePlasticVolume1, blueCut, rotate, G4ThreeVector(-SegPosition, SegPosition,0));
		G4IntersectionSolid * blueBottomLeft= new G4IntersectionSolid("blueBottomLeft", bluePlasticVolume1, blueCut, rotate, G4ThreeVector(SegPosition, -SegPosition,0));
		G4IntersectionSolid * blueBottomRight= new G4IntersectionSolid("blueBottomRight", bluePlasticVolume1, blueCut, rotate, G4ThreeVector(-SegPosition,-SegPosition,0));

		G4Box * bluePMT = new G4Box("bluePMT", pmtSize, pmtSize, pmtSizeZ);

		//For Wrap Subtraction With Larger plastic volumes
		/*		G4SubtractionSolid* bluePlasticVolumeWrap1 = CanVolumeDaemon(true, fPlasticThickness, fBlueDetector, fBluePhi, true);
				G4IntersectionSolid * blueTopLeftWrap = new G4IntersectionSolid("blueTopLeftWrap", bluePlasticVolumeWrap1, blueCut, rotate, G4ThreeVector(SegPositionWrap, SegPositionWrap,0));
				G4IntersectionSolid * blueTopRightWrap = new G4IntersectionSolid("blueTopRightWrap", bluePlasticVolumeWrap1, blueCut, rotate, G4ThreeVector(-SegPositionWrap, SegPositionWrap,0));
				G4IntersectionSolid * blueBottomLeftWrap = new G4IntersectionSolid("blueBottomLeftWrap", bluePlasticVolumeWrap1, blueCut, rotate, G4ThreeVector(SegPositionWrap, -SegPositionWrap,0));
				G4IntersectionSolid * blueBottomRightWrap = new G4IntersectionSolid("blueBottomRightWrap", bluePlasticVolumeWrap1, blueCut, rotate, G4ThreeVector(-SegPositionWrap,-SegPositionWrap,0));		
				G4SubtractionSolid* blueWrapVolumeF = CutWrapperVolumeDaemon(blueWrapVolume1, blueTopRightWrap, blueTopLeftWrap, blueBottomRightWrap, blueBottomLeftWrap, bluePMT);
				*/		
		//For Wrap Subtraction With Normal plastic volumes
		G4SubtractionSolid* blueWrapVolumeF = CutWrapperVolumeDaemon(blueWrapVolume1, blueTopRight, blueTopLeft, blueBottomRight, blueBottomLeft, bluePMT);

		//Subtraction without CutWrapperVolume function
		/*		G4SubtractionSolid* blueWrapVolume2 = new G4SubtractionSolid("blueWrapVolume2", blueWrapVolume1, blueTopLeft, rotate, move);
				G4SubtractionSolid* blueWrapVolume3 = new G4SubtractionSolid("blueWrapVolume3", blueWrapVolume2, blueTopRight, rotate, move);
				G4SubtractionSolid* blueWrapVolume4 = new G4SubtractionSolid("blueWrapVolume4", blueWrapVolume3, blueBottomLeft, rotate, move);
				G4SubtractionSolid* blueWrapVolume5 = new G4SubtractionSolid("blueWrapVolume5", blueWrapVolume4, blueBottomRight, rotate, move);


				move = G4ThreeVector(pmtx, pmty, pmtz);
				G4SubtractionSolid* blueWrapVolume6 = new G4SubtractionSolid("blueWrapVolume6", blueWrapVolume5, bluePMT, rotate, move);
				move = G4ThreeVector(-pmtx, pmty, pmtz);
				G4SubtractionSolid* blueWrapVolume7 = new G4SubtractionSolid("blueWrapVolume7", blueWrapVolume6, bluePMT, rotate, move);
				move = G4ThreeVector(pmtx, -pmty, pmtz);
				G4SubtractionSolid* blueWrapVolume8 = new G4SubtractionSolid("blueWrapVolume8", blueWrapVolume7, bluePMT, rotate, move);
				move = G4ThreeVector(-pmtx, -pmty, pmtz);
				G4SubtractionSolid* blueWrapVolume9 = new G4SubtractionSolid("blueWrapVolume9", blueWrapVolume8, bluePMT, rotate, move);
				*/
		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;


		if(fBlueWrapLog == nullptr) {
			fBlueWrapLog = new G4LogicalVolume(blueWrapVolumeF, wrapG4material, nameWrapperB, 0, 0, 0);
			fBlueWrapLog->SetVisAttributes(blueWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperB, fBlueWrapLog, ScintWrapper);
		}
		if(fBluePMTLogArray[0] == nullptr) {
			fBluePMTLogArray[0] = new G4LogicalVolume(bluePMT, PMTG4material, "BluePMT_1", 0, 0, 0);
			fBluePMTLogArray[0]->SetVisAttributes(PMTVisAtt);
		}
		if(fBluePMTLogArray[1] == nullptr) {
			fBluePMTLogArray[1] = new G4LogicalVolume(bluePMT, PMTG4material, "BluePMT_2", 0, 0, 0);
			fBluePMTLogArray[1]->SetVisAttributes(PMTVisAtt);
		}
		if(fBluePMTLogArray[2] == nullptr) {
			fBluePMTLogArray[2] = new G4LogicalVolume(bluePMT, PMTG4material, "BluePMT_3", 0, 0, 0);
			fBluePMTLogArray[2]->SetVisAttributes(PMTVisAtt);
		}
		if(fBluePMTLogArray[3] == nullptr) {
			fBluePMTLogArray[3] = new G4LogicalVolume(bluePMT, PMTG4material, "BluePMT_4", 0, 0, 0);
			fBluePMTLogArray[3]->SetVisAttributes(PMTVisAtt);
		}
		if(fBluePlasticLogArray[0] == nullptr) {
			fBluePlasticLogArray[0] = new G4LogicalVolume(blueTopLeft, plasticG4material, "BlueTilePlasticDet_1", 0, 0, 0);
			fBluePlasticLogArray[0]->SetVisAttributes(PlasticVisAtt);
		}
		if(fBluePlasticLogArray[1] == nullptr) {
			fBluePlasticLogArray[1] = new G4LogicalVolume(blueTopRight, plasticG4material, "BlueTilePlasticDet_2", 0, 0, 0);
			fBluePlasticLogArray[1]->SetVisAttributes(PlasticVisAtt);
		}
		if(fBluePlasticLogArray[2] == nullptr) {
			fBluePlasticLogArray[2] = new G4LogicalVolume(blueBottomLeft, plasticG4material, "BlueTilePlasticDet_3", 0, 0, 0);
			fBluePlasticLogArray[2]->SetVisAttributes(PlasticVisAtt);
		}
		if(fBluePlasticLogArray[3] == nullptr) {
			fBluePlasticLogArray[3] = new G4LogicalVolume(blueBottomRight, plasticG4material, "BlueTilePlasticDet_4", 0, 0, 0);
			fBluePlasticLogArray[3]->SetVisAttributes(PlasticVisAtt);
		}
		fAssemblyBlue->AddPlacedVolume(fBluePlasticLogArray[0], move, rotate);
		fAssemblyBlue->AddPlacedVolume(fBluePlasticLogArray[1], move, rotate);
		fAssemblyBlue->AddPlacedVolume(fBluePlasticLogArray[2], move, rotate);
		fAssemblyBlue->AddPlacedVolume(fBluePlasticLogArray[3], move, rotate);

		if(fAddWrap){
			fAssemblyBlue->AddPlacedVolume(fBlueWrapLog, move, rotate);
		}
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		//direction 	= G4ThreeVector(0,0,1);
		//move          = zPosition * direction;
		move = G4ThreeVector(pmtx, pmty, zPosition);
		fAssemblyBlue->AddPlacedVolume(fBluePMTLogArray[0], move, rotate);
		move = G4ThreeVector(-pmtx, pmty, zPosition);
		fAssemblyBlue->AddPlacedVolume(fBluePMTLogArray[1], move, rotate);
		move = G4ThreeVector(pmtx, -pmty, zPosition);
		fAssemblyBlue->AddPlacedVolume(fBluePMTLogArray[2], move, rotate);
		move = G4ThreeVector(-pmtx, -pmty, zPosition);
		fAssemblyBlue->AddPlacedVolume(fBluePMTLogArray[3], move, rotate);
	}

	// WHITE Detectors

	if(!fAddSegment){

		G4SubtractionSolid* whitePlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fWhiteDetector, fWhitePhi, false);
		G4SubtractionSolid* whiteWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fWhiteDetector, fWhitePhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;
		G4SubtractionSolid* whiteWrapVolume2 = new G4SubtractionSolid("whiteWrapVolume2", whiteWrapVolume1, whitePlasticVolume1, rotate, move);
		G4Box * whitePMT = new G4Box("whitePMT", pmtSize, pmtSize, pmtSizeZ);
		move = G4ThreeVector(0.0, 0.0, fPlasticThickness/2.); //+pmtSize?
		G4SubtractionSolid* whiteWrapVolume3 = new G4SubtractionSolid("whiteWrapVolume3", whiteWrapVolume2, whitePMT, rotate, move);

		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;

		if(fWhiteWrapLog == nullptr) {
			fWhiteWrapLog = new G4LogicalVolume(whiteWrapVolume3, wrapG4material, nameWrapperW, 0, 0, 0);
			fWhiteWrapLog->SetVisAttributes(whiteWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperW, fWhiteWrapLog, ScintWrapper);
		}
		if(fWhitePMTLog == nullptr) {
			fWhitePMTLog = new G4LogicalVolume(whitePMT, PMTG4material, "WhitePMT_1", 0, 0, 0);
			fWhitePMTLog->SetVisAttributes(PMTVisAtt);
		}
		if(fWhitePlasticLog == nullptr) {
			fWhitePlasticLog = new G4LogicalVolume(whitePlasticVolume1, plasticG4material, "WhiteTilePlasticDet_1", 0, 0, 0);
			fWhitePlasticLog->SetVisAttributes(PlasticVisAtt);
		}

		fAssemblyWhite->AddPlacedVolume(fWhitePlasticLog, move, rotate);
		if(fAddWrap){
			fAssemblyWhite->AddPlacedVolume(fWhiteWrapLog, move, rotate);
		}
		direction 	= G4ThreeVector(0,0,1);
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		move          = zPosition * direction;
		fAssemblyWhite->AddPlacedVolume(fWhitePMTLog, move, rotate);
	}

	if(fAddSegment){

		G4SubtractionSolid* whitePlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fWhiteDetector, fWhitePhi, false);
		G4SubtractionSolid* whiteWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fWhiteDetector, fWhitePhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;

		G4Box * whiteCut = new G4Box("whiteCut", thick, thick, thick);
		//For Plastics
		G4IntersectionSolid * whiteTopLeft= new G4IntersectionSolid("whiteTopLeft", whitePlasticVolume1, whiteCut, rotate, G4ThreeVector(SegPosition, SegPosition,0));
		G4IntersectionSolid * whiteTopRight= new G4IntersectionSolid("whiteTopRight", whitePlasticVolume1, whiteCut, rotate, G4ThreeVector(-SegPosition, SegPosition,0));
		G4IntersectionSolid * whiteBottomLeft= new G4IntersectionSolid("whiteBottomLeft", whitePlasticVolume1, whiteCut, rotate, G4ThreeVector(SegPosition, -SegPosition,0));
		G4IntersectionSolid * whiteBottomRight= new G4IntersectionSolid("whiteBottomRight", whitePlasticVolume1, whiteCut, rotate, G4ThreeVector(-SegPosition,-SegPosition,0));

		G4Box * whitePMT = new G4Box("whitePMT", pmtSize, pmtSize, pmtSizeZ);

		//For Wrap Subtraction With Normal plastic volumes
		G4SubtractionSolid* whiteWrapVolumeF = CutWrapperVolumeDaemon(whiteWrapVolume1, whiteTopRight, whiteTopLeft, whiteBottomRight, whiteBottomLeft, whitePMT);

		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;


		if(fWhiteWrapLog == nullptr) {
			fWhiteWrapLog = new G4LogicalVolume(whiteWrapVolumeF, wrapG4material, nameWrapperW, 0, 0, 0);
			fWhiteWrapLog->SetVisAttributes(whiteWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperW, fWhiteWrapLog, ScintWrapper);
		}
		if(fWhitePMTLogArray[0] == nullptr) {
			fWhitePMTLogArray[0] = new G4LogicalVolume(whitePMT, PMTG4material, "WhitePMT_1", 0, 0, 0);
			fWhitePMTLogArray[0]->SetVisAttributes(PMTVisAtt);
		}
		if(fWhitePMTLogArray[1] == nullptr) {
			fWhitePMTLogArray[1] = new G4LogicalVolume(whitePMT, PMTG4material, "WhitePMT_2", 0, 0, 0);
			fWhitePMTLogArray[1]->SetVisAttributes(PMTVisAtt);
		}
		if(fWhitePMTLogArray[2] == nullptr) {
			fWhitePMTLogArray[2] = new G4LogicalVolume(whitePMT, PMTG4material, "WhitePMT_3", 0, 0, 0);
			fWhitePMTLogArray[2]->SetVisAttributes(PMTVisAtt);
		}
		if(fWhitePMTLogArray[3] == nullptr) {
			fWhitePMTLogArray[3] = new G4LogicalVolume(whitePMT, PMTG4material, "WhitePMT_4", 0, 0, 0);
			fWhitePMTLogArray[3]->SetVisAttributes(PMTVisAtt);
		}
		if(fWhitePlasticLogArray[0] == nullptr) {
			fWhitePlasticLogArray[0] = new G4LogicalVolume(whiteTopLeft, plasticG4material, "WhiteTilePlasticDet_1", 0, 0, 0);
			fWhitePlasticLogArray[0]->SetVisAttributes(PlasticVisAtt);
		}
		if(fWhitePlasticLogArray[1] == nullptr) {
			fWhitePlasticLogArray[1] = new G4LogicalVolume(whiteTopRight, plasticG4material, "WhiteTilePlasticDet_2", 0, 0, 0);
			fWhitePlasticLogArray[1]->SetVisAttributes(PlasticVisAtt);
		}
		if(fWhitePlasticLogArray[2] == nullptr) {
			fWhitePlasticLogArray[2] = new G4LogicalVolume(whiteBottomLeft, plasticG4material, "WhiteTilePlasticDet_3", 0, 0, 0);
			fWhitePlasticLogArray[2]->SetVisAttributes(PlasticVisAtt);
		}
		if(fWhitePlasticLogArray[3] == nullptr) {
			fWhitePlasticLogArray[3] = new G4LogicalVolume(whiteBottomRight, plasticG4material, "WhiteTilePlasticDet_4", 0, 0, 0);
			fWhitePlasticLogArray[3]->SetVisAttributes(PlasticVisAtt);
		}
		fAssemblyWhite->AddPlacedVolume(fWhitePlasticLogArray[0], move, rotate);
		fAssemblyWhite->AddPlacedVolume(fWhitePlasticLogArray[1], move, rotate);
		fAssemblyWhite->AddPlacedVolume(fWhitePlasticLogArray[2], move, rotate);
		fAssemblyWhite->AddPlacedVolume(fWhitePlasticLogArray[3], move, rotate);

		if(fAddWrap){
			fAssemblyWhite->AddPlacedVolume(fWhiteWrapLog, move, rotate);
		}
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		//direction 	= G4ThreeVector(0,0,1);
		//move          = zPosition * direction;
		move = G4ThreeVector(pmtx, pmty, zPosition);
		fAssemblyWhite->AddPlacedVolume(fWhitePMTLogArray[0], move, rotate);
		move = G4ThreeVector(-pmtx, pmty, zPosition);
		fAssemblyWhite->AddPlacedVolume(fWhitePMTLogArray[1], move, rotate);
		move = G4ThreeVector(pmtx, -pmty, zPosition);
		fAssemblyWhite->AddPlacedVolume(fWhitePMTLogArray[2], move, rotate);
		move = G4ThreeVector(-pmtx, -pmty, zPosition);
		fAssemblyWhite->AddPlacedVolume(fWhitePMTLogArray[3], move, rotate);
	}
	// RED Detectors

	if(!fAddSegment){

		G4SubtractionSolid* redPlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fRedDetector, fRedPhi, false);
		G4SubtractionSolid* redWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fRedDetector, fRedPhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;
		G4SubtractionSolid* redWrapVolume2 = new G4SubtractionSolid("redWrapVolume2", redWrapVolume1, redPlasticVolume1, rotate, move);
		G4Box * redPMT = new G4Box("redPMT", pmtSize, pmtSize, pmtSizeZ);
		move = G4ThreeVector(0.0, 0.0, fPlasticThickness/2.); //+pmtSize?
		G4SubtractionSolid* redWrapVolume3 = new G4SubtractionSolid("redWrapVolume3", redWrapVolume2, redPMT, rotate, move);

		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;

		if(fRedWrapLog == nullptr) {
			fRedWrapLog = new G4LogicalVolume(redWrapVolume3, wrapG4material, nameWrapperR, 0, 0, 0);
			fRedWrapLog->SetVisAttributes(redWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperR, fRedWrapLog, ScintWrapper);
		}
		if(fRedPMTLog == nullptr) {
			fRedPMTLog = new G4LogicalVolume(redPMT, PMTG4material, "RedPMT_1", 0, 0, 0);
			fRedPMTLog->SetVisAttributes(PMTVisAtt);
		}
		if(fRedPlasticLog == nullptr) {
			fRedPlasticLog = new G4LogicalVolume(redPlasticVolume1, plasticG4material, "RedTilePlasticDet_1", 0, 0, 0);
			fRedPlasticLog->SetVisAttributes(PlasticVisAtt);
		}

		fAssemblyRed->AddPlacedVolume(fRedPlasticLog, move, rotate);
		if(fAddWrap){
			fAssemblyRed->AddPlacedVolume(fRedWrapLog, move, rotate);
		}
		direction 	= G4ThreeVector(0,0,1);
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		move          = zPosition * direction;
		fAssemblyRed->AddPlacedVolume(fRedPMTLog, move, rotate);
	}

	if(fAddSegment){

		G4SubtractionSolid* redPlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fRedDetector, fRedPhi, false);
		G4SubtractionSolid* redWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fRedDetector, fRedPhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;

		G4Box * redCut = new G4Box("redCut", thick, thick, thick);
		//For Plastics
		G4IntersectionSolid * redTopLeft= new G4IntersectionSolid("redTopLeft", redPlasticVolume1, redCut, rotate, G4ThreeVector(SegPosition, SegPosition,0));
		G4IntersectionSolid * redTopRight= new G4IntersectionSolid("redTopRight", redPlasticVolume1, redCut, rotate, G4ThreeVector(-SegPosition, SegPosition,0));
		G4IntersectionSolid * redBottomLeft= new G4IntersectionSolid("redBottomLeft", redPlasticVolume1, redCut, rotate, G4ThreeVector(SegPosition, -SegPosition,0));
		G4IntersectionSolid * redBottomRight= new G4IntersectionSolid("redBottomRight", redPlasticVolume1, redCut, rotate, G4ThreeVector(-SegPosition,-SegPosition,0));

		G4Box * redPMT = new G4Box("redPMT", pmtSize, pmtSize, pmtSizeZ);

		//For Wrap Subtraction With Normal plastic volumes
		G4SubtractionSolid* redWrapVolumeF = CutWrapperVolumeDaemon(redWrapVolume1, redTopRight, redTopLeft, redBottomRight, redBottomLeft, redPMT);

		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;


		if(fRedWrapLog == nullptr) {
			fRedWrapLog = new G4LogicalVolume(redWrapVolumeF, wrapG4material, nameWrapperR, 0, 0, 0);
			fRedWrapLog->SetVisAttributes(redWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperR, fRedWrapLog, ScintWrapper);
		}
		if(fRedPMTLogArray[0] == nullptr) {
			fRedPMTLogArray[0] = new G4LogicalVolume(redPMT, PMTG4material, "RedPMT_1", 0, 0, 0);
			fRedPMTLogArray[0]->SetVisAttributes(PMTVisAtt);
		}
		if(fRedPMTLogArray[1] == nullptr) {
			fRedPMTLogArray[1] = new G4LogicalVolume(redPMT, PMTG4material, "RedPMT_2", 0, 0, 0);
			fRedPMTLogArray[1]->SetVisAttributes(PMTVisAtt);
		}
		if(fRedPMTLogArray[2] == nullptr) {
			fRedPMTLogArray[2] = new G4LogicalVolume(redPMT, PMTG4material, "RedPMT_3", 0, 0, 0);
			fRedPMTLogArray[2]->SetVisAttributes(PMTVisAtt);
		}
		if(fRedPMTLogArray[3] == nullptr) {
			fRedPMTLogArray[3] = new G4LogicalVolume(redPMT, PMTG4material, "RedPMT_4", 0, 0, 0);
			fRedPMTLogArray[3]->SetVisAttributes(PMTVisAtt);
		}
		if(fRedPlasticLogArray[0] == nullptr) {
			fRedPlasticLogArray[0] = new G4LogicalVolume(redTopLeft, plasticG4material, "RedTilePlasticDet_1", 0, 0, 0);
			fRedPlasticLogArray[0]->SetVisAttributes(PlasticVisAtt);
		}
		if(fRedPlasticLogArray[1] == nullptr) {
			fRedPlasticLogArray[1] = new G4LogicalVolume(redTopRight, plasticG4material, "RedTilePlasticDet_2", 0, 0, 0);
			fRedPlasticLogArray[1]->SetVisAttributes(PlasticVisAtt);
		}
		if(fRedPlasticLogArray[2] == nullptr) {
			fRedPlasticLogArray[2] = new G4LogicalVolume(redBottomLeft, plasticG4material, "RedTilePlasticDet_3", 0, 0, 0);
			fRedPlasticLogArray[2]->SetVisAttributes(PlasticVisAtt);
		}
		if(fRedPlasticLogArray[3] == nullptr) {
			fRedPlasticLogArray[3] = new G4LogicalVolume(redBottomRight, plasticG4material, "RedTilePlasticDet_4", 0, 0, 0);
			fRedPlasticLogArray[3]->SetVisAttributes(PlasticVisAtt);
		}
		fAssemblyRed->AddPlacedVolume(fRedPlasticLogArray[0], move, rotate);
		fAssemblyRed->AddPlacedVolume(fRedPlasticLogArray[1], move, rotate);
		fAssemblyRed->AddPlacedVolume(fRedPlasticLogArray[2], move, rotate);
		fAssemblyRed->AddPlacedVolume(fRedPlasticLogArray[3], move, rotate);

		if(fAddWrap){
			fAssemblyRed->AddPlacedVolume(fRedWrapLog, move, rotate);
		}
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		//direction 	= G4ThreeVector(0,0,1);
		//move          = zPosition * direction;
		move = G4ThreeVector(pmtx, pmty, zPosition);
		fAssemblyRed->AddPlacedVolume(fRedPMTLogArray[0], move, rotate);
		move = G4ThreeVector(-pmtx, pmty, zPosition);
		fAssemblyRed->AddPlacedVolume(fRedPMTLogArray[1], move, rotate);
		move = G4ThreeVector(pmtx, -pmty, zPosition);
		fAssemblyRed->AddPlacedVolume(fRedPMTLogArray[2], move, rotate);
		move = G4ThreeVector(-pmtx, -pmty, zPosition);
		fAssemblyRed->AddPlacedVolume(fRedPMTLogArray[3], move, rotate);
	}
	// Green Detectors

	if(!fAddSegment){

		G4SubtractionSolid* greenPlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fGreenDetector, fGreenPhi, false);
		G4SubtractionSolid* greenWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fGreenDetector, fGreenPhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;
		G4SubtractionSolid* greenWrapVolume2 = new G4SubtractionSolid("greenWrapVolume2", greenWrapVolume1, greenPlasticVolume1, rotate, move);
		G4Box * greenPMT = new G4Box("greenPMT", pmtSize, pmtSize, pmtSizeZ);
		move = G4ThreeVector(0.0, 0.0, fPlasticThickness/2.); //+pmtSize?
		G4SubtractionSolid* greenWrapVolume3 = new G4SubtractionSolid("greenWrapVolume3", greenWrapVolume2, greenPMT, rotate, move);

		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;

		if(fGreenWrapLog == nullptr) {
			fGreenWrapLog = new G4LogicalVolume(greenWrapVolume3, wrapG4material, nameWrapperG, 0, 0, 0);
			fGreenWrapLog->SetVisAttributes(greenWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperG, fGreenWrapLog, ScintWrapper);
		}
		if(fGreenPMTLog == nullptr) {
			fGreenPMTLog = new G4LogicalVolume(greenPMT, PMTG4material, "GreenPMT_1", 0, 0, 0);
			fGreenPMTLog->SetVisAttributes(PMTVisAtt);
		}
		if(fGreenPlasticLog == nullptr) {
			fGreenPlasticLog = new G4LogicalVolume(greenPlasticVolume1, plasticG4material, "GreenTilePlasticDet_1", 0, 0, 0);
			fGreenPlasticLog->SetVisAttributes(PlasticVisAtt);
		}

		fAssemblyGreen->AddPlacedVolume(fGreenPlasticLog, move, rotate);
		if(fAddWrap){
			fAssemblyGreen->AddPlacedVolume(fGreenWrapLog, move, rotate);
		}
		direction 	= G4ThreeVector(0,0,1);
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		move          = zPosition * direction;
		fAssemblyGreen->AddPlacedVolume(fGreenPMTLog, move, rotate);
	}

	if(fAddSegment){

		G4SubtractionSolid* greenPlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fGreenDetector, fGreenPhi, false);
		G4SubtractionSolid* greenWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fGreenDetector, fGreenPhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;

		G4Box * greenCut = new G4Box("greenCut", thick, thick, thick);
		//For Plastics
		G4IntersectionSolid * greenTopLeft= new G4IntersectionSolid("greenTopLeft", greenPlasticVolume1, greenCut, rotate, G4ThreeVector(SegPosition, SegPosition,0));
		G4IntersectionSolid * greenTopRight= new G4IntersectionSolid("greenTopRight", greenPlasticVolume1, greenCut, rotate, G4ThreeVector(-SegPosition, SegPosition,0));
		G4IntersectionSolid * greenBottomLeft= new G4IntersectionSolid("greenBottomLeft", greenPlasticVolume1, greenCut, rotate, G4ThreeVector(SegPosition, -SegPosition,0));
		G4IntersectionSolid * greenBottomRight= new G4IntersectionSolid("greenBottomRight", greenPlasticVolume1, greenCut, rotate, G4ThreeVector(-SegPosition,-SegPosition,0));

		G4Box * greenPMT = new G4Box("greenPMT", pmtSize, pmtSize, pmtSizeZ);

		//For Wrap Subtraction With Normal plastic volumes
		G4SubtractionSolid* greenWrapVolumeF = CutWrapperVolumeDaemon(greenWrapVolume1, greenTopRight, greenTopLeft, greenBottomRight, greenBottomLeft, greenPMT);

		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;


		if(fGreenWrapLog == nullptr) {
			fGreenWrapLog = new G4LogicalVolume(greenWrapVolumeF, wrapG4material, nameWrapperG, 0, 0, 0);
			fGreenWrapLog->SetVisAttributes(greenWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperG, fGreenWrapLog, ScintWrapper);
		}
		if(fGreenPMTLogArray[0] == nullptr) {
			fGreenPMTLogArray[0] = new G4LogicalVolume(greenPMT, PMTG4material, "GreenPMT_1", 0, 0, 0);
			fGreenPMTLogArray[0]->SetVisAttributes(PMTVisAtt);
		}
		if(fGreenPMTLogArray[1] == nullptr) {
			fGreenPMTLogArray[1] = new G4LogicalVolume(greenPMT, PMTG4material, "GreenPMT_2", 0, 0, 0);
			fGreenPMTLogArray[1]->SetVisAttributes(PMTVisAtt);
		}
		if(fGreenPMTLogArray[2] == nullptr) {
			fGreenPMTLogArray[2] = new G4LogicalVolume(greenPMT, PMTG4material, "GreenPMT_3", 0, 0, 0);
			fGreenPMTLogArray[2]->SetVisAttributes(PMTVisAtt);
		}
		if(fGreenPMTLogArray[3] == nullptr) {
			fGreenPMTLogArray[3] = new G4LogicalVolume(greenPMT, PMTG4material, "GreenPMT_4", 0, 0, 0);
			fGreenPMTLogArray[3]->SetVisAttributes(PMTVisAtt);
		}
		if(fGreenPlasticLogArray[0] == nullptr) {
			fGreenPlasticLogArray[0] = new G4LogicalVolume(greenTopLeft, plasticG4material, "GreenTilePlasticDet_1", 0, 0, 0);
			fGreenPlasticLogArray[0]->SetVisAttributes(PlasticVisAtt);
		}
		if(fGreenPlasticLogArray[1] == nullptr) {
			fGreenPlasticLogArray[1] = new G4LogicalVolume(greenTopRight, plasticG4material, "GreenTilePlasticDet_2", 0, 0, 0);
			fGreenPlasticLogArray[1]->SetVisAttributes(PlasticVisAtt);
		}
		if(fGreenPlasticLogArray[2] == nullptr) {
			fGreenPlasticLogArray[2] = new G4LogicalVolume(greenBottomLeft, plasticG4material, "GreenTilePlasticDet_3", 0, 0, 0);
			fGreenPlasticLogArray[2]->SetVisAttributes(PlasticVisAtt);
		}
		if(fGreenPlasticLogArray[3] == nullptr) {
			fGreenPlasticLogArray[3] = new G4LogicalVolume(greenBottomRight, plasticG4material, "GreenTilePlasticDet_4", 0, 0, 0);
			fGreenPlasticLogArray[3]->SetVisAttributes(PlasticVisAtt);
		}
		fAssemblyGreen->AddPlacedVolume(fGreenPlasticLogArray[0], move, rotate);
		fAssemblyGreen->AddPlacedVolume(fGreenPlasticLogArray[1], move, rotate);
		fAssemblyGreen->AddPlacedVolume(fGreenPlasticLogArray[2], move, rotate);
		fAssemblyGreen->AddPlacedVolume(fGreenPlasticLogArray[3], move, rotate);

		if(fAddWrap){
			fAssemblyGreen->AddPlacedVolume(fGreenWrapLog, move, rotate);
		}
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		//direction 	= G4ThreeVector(0,0,1);
		//move          = zPosition * direction;
		move = G4ThreeVector(pmtx, pmty, zPosition);
		fAssemblyGreen->AddPlacedVolume(fGreenPMTLogArray[0], move, rotate);
		move = G4ThreeVector(-pmtx, pmty, zPosition);
		fAssemblyGreen->AddPlacedVolume(fGreenPMTLogArray[1], move, rotate);
		move = G4ThreeVector(pmtx, -pmty, zPosition);
		fAssemblyGreen->AddPlacedVolume(fGreenPMTLogArray[2], move, rotate);
		move = G4ThreeVector(-pmtx, -pmty, zPosition);
		fAssemblyGreen->AddPlacedVolume(fGreenPMTLogArray[3], move, rotate);
	}

	// Yellow Detectors

	if(!fAddSegment){

		G4SubtractionSolid* yellowPlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fYellowDetector, fYellowPhi, false);
		G4SubtractionSolid* yellowWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fYellowDetector, fYellowPhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;
		G4SubtractionSolid* yellowWrapVolume2 = new G4SubtractionSolid("yellowWrapVolume2", yellowWrapVolume1, yellowPlasticVolume1, rotate, move);
		G4Box * yellowPMT = new G4Box("yellowPMT", pmtSize, pmtSize, pmtSizeZ);
		move = G4ThreeVector(0.0, 0.0, fPlasticThickness/2.); //+pmtSize?
		G4SubtractionSolid* yellowWrapVolume3 = new G4SubtractionSolid("yellowWrapVolume3", yellowWrapVolume2, yellowPMT, rotate, move);

		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;

		if(fYellowWrapLog == nullptr) {
			fYellowWrapLog = new G4LogicalVolume(yellowWrapVolume3, wrapG4material, nameWrapperY, 0, 0, 0);
			fYellowWrapLog->SetVisAttributes(yellowWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperY, fYellowWrapLog, ScintWrapper);
		}
		if(fYellowPMTLog == nullptr) {
			fYellowPMTLog = new G4LogicalVolume(yellowPMT, PMTG4material, "YellowPMT_1", 0, 0, 0);
			fYellowPMTLog->SetVisAttributes(PMTVisAtt);
		}
		if(fYellowPlasticLog == nullptr) {
			fYellowPlasticLog = new G4LogicalVolume(yellowPlasticVolume1, plasticG4material, "YellowTilePlasticDet_1", 0, 0, 0);
			fYellowPlasticLog->SetVisAttributes(PlasticVisAtt);
		}

		fAssemblyYellow->AddPlacedVolume(fYellowPlasticLog, move, rotate);
		if(fAddWrap){
			fAssemblyYellow->AddPlacedVolume(fYellowWrapLog, move, rotate);
		}
		direction 	= G4ThreeVector(0,0,1);
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		move          = zPosition * direction;
		fAssemblyYellow->AddPlacedVolume(fYellowPMTLog, move, rotate);
	}

	if(fAddSegment){

		G4SubtractionSolid* yellowPlasticVolume1 = CanVolumeDaemon(true, fPlasticThickness, fYellowDetector, fYellowPhi, false);
		G4SubtractionSolid* yellowWrapVolume1 = CanVolumeDaemon(false, fPlasticThickness+fWrapThickness, fYellowDetector, fYellowPhi, false);
		move = G4ThreeVector(0.0, 0.0, 0.0);
		rotate = new G4RotationMatrix;

		G4Box * yellowCut = new G4Box("yellowCut", thick, thick, thick);
		//For Plastics
		G4IntersectionSolid * yellowTopLeft= new G4IntersectionSolid("yellowTopLeft", yellowPlasticVolume1, yellowCut, rotate, G4ThreeVector(SegPosition, SegPosition,0));
		G4IntersectionSolid * yellowTopRight= new G4IntersectionSolid("yellowTopRight", yellowPlasticVolume1, yellowCut, rotate, G4ThreeVector(-SegPosition, SegPosition,0));
		G4IntersectionSolid * yellowBottomLeft= new G4IntersectionSolid("yellowBottomLeft", yellowPlasticVolume1, yellowCut, rotate, G4ThreeVector(SegPosition, -SegPosition,0));
		G4IntersectionSolid * yellowBottomRight= new G4IntersectionSolid("yellowBottomRight", yellowPlasticVolume1, yellowCut, rotate, G4ThreeVector(-SegPosition,-SegPosition,0));

		G4Box * yellowPMT = new G4Box("yellowPMT", pmtSize, pmtSize, pmtSizeZ);

		//For Wrap Subtraction With Normal plastic volumes
		G4SubtractionSolid* yellowWrapVolumeF = CutWrapperVolumeDaemon(yellowWrapVolume1, yellowTopRight, yellowTopLeft, yellowBottomRight, yellowBottomLeft, yellowPMT);

		// Define rotation and movement objects
		direction 	= G4ThreeVector(0,0,1);
		// fWrapThickness/2.  stops overlap with descant
		zPosition    = fPlasticThickness/2.0+fWrapThickness/2.;
		move          = zPosition * direction;

		rotate = new G4RotationMatrix;


		if(fYellowWrapLog == nullptr) {
			fYellowWrapLog = new G4LogicalVolume(yellowWrapVolumeF, wrapG4material, nameWrapperY, 0, 0, 0);
			fYellowWrapLog->SetVisAttributes(yellowWrapVisAtt);
			//Set Logical Skin for optical photons on wrapping
			G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapperY, fYellowWrapLog, ScintWrapper);
		}
		if(fYellowPMTLogArray[0] == nullptr) {
			fYellowPMTLogArray[0] = new G4LogicalVolume(yellowPMT, PMTG4material, "YellowPMT_1", 0, 0, 0);
			fYellowPMTLogArray[0]->SetVisAttributes(PMTVisAtt);
		}
		if(fYellowPMTLogArray[1] == nullptr) {
			fYellowPMTLogArray[1] = new G4LogicalVolume(yellowPMT, PMTG4material, "YellowPMT_2", 0, 0, 0);
			fYellowPMTLogArray[1]->SetVisAttributes(PMTVisAtt);
		}
		if(fYellowPMTLogArray[2] == nullptr) {
			fYellowPMTLogArray[2] = new G4LogicalVolume(yellowPMT, PMTG4material, "YellowPMT_3", 0, 0, 0);
			fYellowPMTLogArray[2]->SetVisAttributes(PMTVisAtt);
		}
		if(fYellowPMTLogArray[3] == nullptr) {
			fYellowPMTLogArray[3] = new G4LogicalVolume(yellowPMT, PMTG4material, "YellowPMT_4", 0, 0, 0);
			fYellowPMTLogArray[3]->SetVisAttributes(PMTVisAtt);
		}
		if(fYellowPlasticLogArray[0] == nullptr) {
			fYellowPlasticLogArray[0] = new G4LogicalVolume(yellowTopLeft, plasticG4material, "YellowTilePlasticDet_1", 0, 0, 0);
			fYellowPlasticLogArray[0]->SetVisAttributes(PlasticVisAtt);
		}
		if(fYellowPlasticLogArray[1] == nullptr) {
			fYellowPlasticLogArray[1] = new G4LogicalVolume(yellowTopRight, plasticG4material, "YellowTilePlasticDet_2", 0, 0, 0);
			fYellowPlasticLogArray[1]->SetVisAttributes(PlasticVisAtt);
		}
		if(fYellowPlasticLogArray[2] == nullptr) {
			fYellowPlasticLogArray[2] = new G4LogicalVolume(yellowBottomLeft, plasticG4material, "YellowTilePlasticDet_3", 0, 0, 0);
			fYellowPlasticLogArray[2]->SetVisAttributes(PlasticVisAtt);
		}
		if(fYellowPlasticLogArray[3] == nullptr) {
			fYellowPlasticLogArray[3] = new G4LogicalVolume(yellowBottomRight, plasticG4material, "YellowTilePlasticDet_4", 0, 0, 0);
			fYellowPlasticLogArray[3]->SetVisAttributes(PlasticVisAtt);
		}
		fAssemblyYellow->AddPlacedVolume(fYellowPlasticLogArray[0], move, rotate);
		fAssemblyYellow->AddPlacedVolume(fYellowPlasticLogArray[1], move, rotate);
		fAssemblyYellow->AddPlacedVolume(fYellowPlasticLogArray[2], move, rotate);
		fAssemblyYellow->AddPlacedVolume(fYellowPlasticLogArray[3], move, rotate);

		if(fAddWrap){
			fAssemblyYellow->AddPlacedVolume(fYellowWrapLog, move, rotate);
		}
		zPosition    = fPlasticThickness + pmtSizeZ + fWrapThickness/2.;
		//direction 	= G4ThreeVector(0,0,1);
		//move          = zPosition * direction;
		move = G4ThreeVector(pmtx, pmty, zPosition);
		fAssemblyYellow->AddPlacedVolume(fYellowPMTLogArray[0], move, rotate);
		move = G4ThreeVector(-pmtx, pmty, zPosition);
		fAssemblyYellow->AddPlacedVolume(fYellowPMTLogArray[1], move, rotate);
		move = G4ThreeVector(pmtx, -pmty, zPosition);
		fAssemblyYellow->AddPlacedVolume(fYellowPMTLogArray[2], move, rotate);
		move = G4ThreeVector(-pmtx, -pmty, zPosition);
		fAssemblyYellow->AddPlacedVolume(fYellowPMTLogArray[3], move, rotate);
	}

	return 1;
}
///////////////////////////////////////////////////////////////////////
// Methods used to build shapes Daemon
///////////////////////////////////////////////////////////////////////
// This function constructs a generic 6 sided shape defined by 12 input points, and the cut phi angles.
G4SubtractionSolid* DetectionSystemDaemonTiles::CanVolumeDaemon(G4bool insideVol, G4double volumeLength, G4double detector[12][3], G4double detectorPhi[6], G4bool wrap)
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

	G4SubtractionSolid* canVolumeCut2 = CutVolumeOnFourPointsDaemon(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut1, frontP1, frontP2, wrap);

	// next set of points
	idx1 = 1;
	idx2 = idx1+1;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut3 = CutVolumeOnFourPointsDaemon(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut2, frontP1, frontP2, wrap);

	idx1 = 2;
	idx2 = idx1+1;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut4 = CutVolumeOnFourPointsDaemon(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut3, frontP1, frontP2, wrap);

	idx1 = 3;
	idx2 = idx1+1;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut5 = CutVolumeOnFourPointsDaemon(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut4, frontP1, frontP2, wrap);

	idx1 = 4;
	idx2 = idx1+1;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut6 = CutVolumeOnFourPointsDaemon(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut5, frontP1, frontP2, wrap);

	idx1 = 5;
	idx2 = 0;

	frontP1    = G4ThreeVector(detector[0+idx1][0],detector[0+idx1][1],detector[0+idx1][2]);
	frontP2    = G4ThreeVector(detector[0+idx2][0],detector[0+idx2][1],detector[0+idx2][2]);

	G4SubtractionSolid* canVolumeCut7 = CutVolumeOnFourPointsDaemon(idx1, insideVol, volumeLength, detectorPhi, canVolumeCut6, frontP1, frontP2, wrap);

	return canVolumeCut7;
}

// This method does the actual cutting of the surface.
G4SubtractionSolid* DetectionSystemDaemonTiles::CutVolumeOnFourPointsDaemon(G4int idx, G4bool insideVol, G4double volumeLength, G4double detectorPhi[6], G4SubtractionSolid* volume, G4ThreeVector frontP1, G4ThreeVector frontP2, G4bool wrap)
{
	//G4double maxZ = 500.0*mm;
	G4double maxZ = 495.0*mm;
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
	//have y = m*x + b.  Now wnt to find the min x and y to go from origin to the line. therefore setup y = 1/m*x+c is perpendicular line from origin. 
	//Sub in point on new line (0, 0) since comes from origin, get c = 0.
	//Want point where 2 lines meet: set y = 1/m *x = m*x+b.  Get formula below.
	G4double minx   = (-1.0*slope*intcp)/(1+slope*slope);
	G4double miny   = slope*minx + intcp;

	// make vectors, and get angle between vectors
	G4ThreeVector newV1 = G4ThreeVector(minx,miny, maxZ);
	G4ThreeVector newV2 = G4ThreeVector(0.0*mm,0.0*mm,maxZ);
	//Theta is definition of angle between vectors
	theta   = acos((newV1.dot(newV2))/(newV1.mag() * newV2.mag()));
	phi     = detectorPhi[idx];

	if(insideVol && !wrap)
		moveZCut = 1.0*(volumeLength/2.0 + maxZ  - 3.*fWrapThickness); //Adjusts for cutting inside the can
	else if(insideVol && wrap)
		moveZCut = 1.0*(volumeLength/2.0 + maxZ  - 3.*fWrapThickness+fAirGap); //Adjusts for cutting inside the can
	//moveZCut = 1.0*(volumeLength/2.0 + maxZ  - fWrapThickness); //Adjusts for cutting inside the can
	//moveZCut = 1.0*(volumeLength/2.0 + maxZ - 15.*(((fWrapThickness)/(sin(theta))) - fWrapThickness)); //Adjusts for cutting inside the can
	else
		moveZCut = 1.0*(volumeLength/2.0 + maxZ - 2.*fWrapThickness);
	//moveZCut = 1.0*(volumeLength/2.0 + maxZ - 10.*(((fWrapThickness)/(sin(theta))) - fWrapThickness)); //Adjusts for cutting inside the can

	//if(volumeLength == (fPlasticThickness+fWrapThickness) ) // this is the lead shield, push it a bit closer so the tapers match up
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


G4SubtractionSolid* DetectionSystemDaemonTiles::CutWrapperVolumeDaemon(G4SubtractionSolid * initialW, G4IntersectionSolid * topRight, G4IntersectionSolid * topLeft, G4IntersectionSolid * bottomRight, G4IntersectionSolid * bottomLeft, G4Box * pmt){

	G4ThreeVector move;
	G4RotationMatrix* rotate;
	move = G4ThreeVector(0.0, 0.0, 0.0);
	rotate = new G4RotationMatrix;

	//Cut out the wholes for the plastic
	G4SubtractionSolid* WrapVolume1 = new G4SubtractionSolid("WrapVolume1", initialW, topLeft, rotate, move);
	G4SubtractionSolid* WrapVolume2 = new G4SubtractionSolid("WrapVolume2", WrapVolume1, topRight, rotate, move);
	G4SubtractionSolid* WrapVolume3 = new G4SubtractionSolid("WrapVolume3", WrapVolume2, bottomRight, rotate, move);
	G4SubtractionSolid* WrapVolume4 = new G4SubtractionSolid("WrapVolume4", WrapVolume3, bottomLeft, rotate, move);

	G4double pmtx = 2.5*cm;
	G4double pmty = 2.5*cm;
	G4double pmtz = fPlasticThickness/2.;

	//Cut out the wholes for the pmts
	move = G4ThreeVector(pmtx, pmty, pmtz);
	G4SubtractionSolid* WrapVolume5 = new G4SubtractionSolid("WrapVolume5", WrapVolume4, pmt, rotate, move);
	move = G4ThreeVector(-pmtx, pmty, pmtz);
	G4SubtractionSolid* WrapVolume6 = new G4SubtractionSolid("WrapVolume6", WrapVolume5, pmt, rotate, move);
	move = G4ThreeVector(pmtx, -pmty, pmtz);
	G4SubtractionSolid* WrapVolume7 = new G4SubtractionSolid("WrapVolume7", WrapVolume6, pmt, rotate, move);
	move = G4ThreeVector(-pmtx, -pmty, pmtz);
	G4SubtractionSolid* WrapVolume8 = new G4SubtractionSolid("WrapVolume8", WrapVolume7, pmt, rotate, move);

	return WrapVolume8;

}





//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemDaemonTiles::GetDirectionXYZ(G4double theta, G4double phi)
{
	G4double x,y,z;
	x = sin(theta) * cos(phi);
	y = sin(theta) * sin(phi);
	z = cos(theta);

	G4ThreeVector direction = G4ThreeVector(x,y,z);

	return direction;
}

G4ThreeVector DetectionSystemDaemonTiles::SolveLineEquationX(G4ThreeVector p1, G4ThreeVector p2, G4double x){
	G4ThreeVector v, result;
	v = p2 - p1;
	G4double t = (x - p1.x())/(v.x());
	result = G4ThreeVector(x,p1.y()+v.y()*t,p1.z()+v.z()*t);
	return result;
}

G4ThreeVector DetectionSystemDaemonTiles::SolveLineEquationY(G4ThreeVector p1, G4ThreeVector p2, G4double y){
	G4ThreeVector v, result;
	v = p2 - p1;
	G4double t = (y - p1.y())/(v.y());
	result = G4ThreeVector(p1.x()+v.x()*t,y,p1.z()+v.z()*t);
	return result;
}

G4ThreeVector DetectionSystemDaemonTiles::SolveLineEquationZ(G4ThreeVector p1, G4ThreeVector p2, G4double z){
	G4ThreeVector v, result;
	v = p2 - p1;
	G4double t = (z - p1.z())/(v.z());
	result = G4ThreeVector(p1.x()+v.x()*t,p1.y()+v.y()*t,z);
	return result;
}
