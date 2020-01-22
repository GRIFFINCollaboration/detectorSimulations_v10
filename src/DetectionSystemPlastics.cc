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
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//// Added ///  G4MaterialTable.hh is alreadu elsewhere
#include "G4OpticalSurface.hh"
#include "G4OpticalPhysics.hh"
#include "G4LogicalSkinSurface.hh"
////////

#include "DetectionSystemPlastics.hh"

#include "G4SystemOfUnits.hh"

#include <string>

DetectionSystemPlastics::DetectionSystemPlastics(G4double thickness, G4int material, G4double numDet)
{
	fAddWrap = true;
	// detector dimensions  properties
	fScintillatorLength        = 6.*cm;
	fScintillatorHeight        = 6.*cm;
	fScintillatorWidth         = thickness;
	fRadialDistance = 50*cm;
	fLeadShieldThickness = 6.35*mm;
	fNumDet = numDet ;
	fPlasticLogArray.resize(numDet, NULL);
	fWrapLogArray.resize(numDet, NULL);
	fPMT1LogArray.resize(numDet, NULL);
	fPMT2LogArray.resize(numDet, NULL);

	fWrapThickness = 0.5 * mm; //Factor of 10 should be  applied later for visualization purposes.  0.05 for simulations, 0.5 for visualization
	fSpacing = fWrapThickness; //assuming no lead on DESCANT but taking into account the optical wrapping
	//fWrapThickness = 0.1 * cm;
	//fAirGap = 0.1 * fWrapThickness;
	fAirGap = 0.;
	fWrapMaterial = "Teflon";
	fPMTMaterial = "G4_SILICON_DIOXIDE";

	//blue=G4Color(0.,0.,1.);
	bronze=G4Color(0.8,0.5,0.2);
	//cyan=G4Color(0.,1.,1.);
	silver=G4Color(0.75,0.75,0.75);
	black=G4Color(0.,0.,0.);

	if(material == 1)  fPlasticMaterial = "BC408";
	else if (material == 2) fPlasticMaterial = "deuterium";
	else if (material == 3) fPlasticMaterial = "Hydrogen";
	else if (material == 4) fPlasticMaterial = "Carbon";
	else if (material == 5) fPlasticMaterial = "Deuterated Scintillator";
	else if (material == 6) fPlasticMaterial = "BC537";
	else G4cout<< "Material Unknown" << G4endl;


	G4cout << "Calling Constructor" << G4endl;
}
/////
///////
DetectionSystemPlastics::~DetectionSystemPlastics() {
	// LogicalVolumes
	for (int i = 0; i<(fNumDet+fNumDetBot); i++) {
		delete fPlasticLogArray[i];
		delete fWrapLogArray[i];
		delete fPMT1LogArray[i];
		delete fPMT2LogArray[i];
	}
	G4cout << "Calling Destructor" << G4endl;

}
////////
/////////
G4int DetectionSystemPlastics::Build() {
	fAssemblyPlastics = new G4AssemblyVolume(); 
	G4cout << "Calling Build function" << G4endl;
	BuildPlastics();

	return 1;
}
////////
///////
G4int DetectionSystemPlastics::PlaceDetector(G4LogicalVolume* expHallLog) {
	G4RotationMatrix * rotate = new G4RotationMatrix;
	G4ThreeVector move = G4ThreeVector(0., 0., 0.);

	//To not check overlaps
	//    fAssemblyPlastics->MakeImprint(expHallLog, move, rotate);
	//To check overlaps
	fAssemblyPlastics->MakeImprint(expHallLog, move, rotate, 0, true);
	G4cout << "Calling place detector" << G4endl;
	return 1;
}
///////////
/////////
G4int DetectionSystemPlastics::BuildPlastics() {

	G4cout << "Calling Build PLastics" << G4endl;


	G4ThreeVector move, direction;
	G4RotationMatrix* rotate;

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
		G4cout << PMTG4material->GetName() << " is the name of the wrapping material" << G4endl;
	}

	////////Scintillation Properties ////////  --------- Might have to be put before the material is constructed
	//Based on BC408 data
	G4MaterialPropertiesTable * scintillatorMPT = new G4MaterialPropertiesTable();

	////Can potentially comment this all out because i dont know what e yield should be.
	//If I dont specify ny own scintillation yield for p,d,t,a,C then they all default to the electron.. Might be a good test...
	//This is SCINTILLATION YIELD as a function of energy and particle type,
	//Have to uncomment line in ConstructOp that allows for this to work with boolean (true)
	/*
		const G4int num2 = 4;
		G4double e_test[num2] = {1.*keV, 0.1*MeV, 1.*MeV, 10.*MeV};
		G4double p_test[num2] = {1.*keV, 0.1*MeV, 1.*MeV, 10.*MeV};
		G4double d_test[num2] = {1.*keV, 0.1*MeV, 1.*MeV, 10.*MeV};
		G4double t_test[num2] = {1.*keV, 0.1*MeV, 1.*MeV, 10.*MeV};
		G4double a_test[num2] = {1.*keV, 0.1*MeV, 1.*MeV, 10.*MeV};
		G4double C_test[num2] = {1.*keV, 0.1*MeV, 1.*MeV, 10.*MeV};
		G4double num_test[num2] = {10., 1000., 10000., 100000.};
		assert(sizeof(e_test) == sizeof(num_test));
		assert(sizeof(p_test) == sizeof(num_test));
		assert(sizeof(d_test) == sizeof(num_test));
		assert(sizeof(t_test) == sizeof(num_test));
		assert(sizeof(a_test) == sizeof(num_test));
		assert(sizeof(C_test) == sizeof(num_test));

		scintillatorMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", e_test, num_test, num2);
		scintillatorMPT->AddProperty("PROTONSCINTILLATIONYIELD", p_test, num_test, num2);
		scintillatorMPT->AddProperty("DEUTERONSCINTILLATIONYIELD", d_test, num_test, num2);
		scintillatorMPT->AddProperty("TRITONSCINTILLATIONYIELD", t_test, num_test, num2);
		scintillatorMPT->AddProperty("ALPHASCINTILLATIONYIELD", a_test, num_test, num2);
		scintillatorMPT->AddProperty("IONSCINTILLATIONYIELD", C_test, num_test, num2);
	///////
	//
	*/
	const G4int num = 12;

	G4double photonEnergy[num] = {1.7*eV, 2.38*eV, 2.48*eV, 2.58*eV, 2.70*eV, 2.76*eV, 2.82*eV, 2.91*eV, 2.95*eV, 3.1*eV, 3.26*eV, 3.44*eV}; //BC408 emission spectra & corresponding energies

	G4double RIndex1[num] = {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58};
	assert(sizeof(RIndex1) == sizeof(photonEnergy));
	const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

	//G4double absorption[num] = {380.*cm, 380*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm, 380.*cm}; ///light attenuation
	//assert(sizeof(absorption) == sizeof(photonEnergy));

	G4double scint[num] = {3., 3., 8., 18., 43., 55., 80., 100., 80., 20., 7., 3. }; ///// Based off emission spectra for BC408
	assert(sizeof(scint) == sizeof(photonEnergy));


	scintillatorMPT->AddProperty("FASTCOMPONENT", photonEnergy, scint, nEntries)->SetSpline(true); // BC408 emission spectra
	scintillatorMPT->AddProperty("SLOWCOMPONENT", photonEnergy, scint, nEntries)->SetSpline(true); // BC408 emission spectra
	scintillatorMPT->AddProperty("RINDEX", photonEnergy, RIndex1, nEntries);  //refractive index can change with energy
	//note if photon is created outside of energy range it will have no index of refraction
	//scintillatorMPT->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries); //absorption length doesnt change with energy - examples showing it can...
	scintillatorMPT->AddConstProperty("ABSLENGTH", 380.*cm); //Scintillation Efficiency - characteristic light yield
	scintillatorMPT->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV); //Scintillation Efficiency - characteristic light yield //10000./MeV
	scintillatorMPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // broadens the statistical distribution of generated photons, gaussian based on SCINTILLATIONYIELD, >1 broadens, 0 no distribution
	scintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns); //only one decay constant given
	scintillatorMPT->AddConstProperty("SLOWTIMECONSTANT", 2.1*ns); //only one decay constant given - triplet-triplet annihilation 
	scintillatorMPT->AddConstProperty("YIELDRATIO", 1.0); //The relative strength of the fast component as a fraction of total scintillation yield is given by the YIELDRATIO.
	//Should these be in the physics list?
	G4OpticalPhysics * opticalPhysics = new G4OpticalPhysics();
	opticalPhysics->SetFiniteRiseTime(true);
	scintillatorMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.7*ns); //default rise time is 0ns, have to set manually
	scintillatorMPT->AddConstProperty("SLOWSCINTILLATIONRISETIME", 0.7*ns); //default rise time is 0ns, have to set manually
	//properties I may be missing: scintillation, rayleigh
	plasticG4material->SetMaterialPropertiesTable(scintillatorMPT);

	const G4int numShort =3;
	G4double photonEnergyShort[numShort] = {1.7*eV,   2.82*eV, 3.44*eV}; //BC408 emission spectra & corresponding energies
	const G4int nEntriesShort = sizeof(photonEnergyShort)/sizeof(G4double);
	//////Optical Surface - Teflon wrapping //////
	G4OpticalSurface * ScintWrapper = new G4OpticalSurface("wrapper");
	///Test 1
	ScintWrapper->SetModel(unified);  // unified or glisur
	ScintWrapper->SetType(dielectric_dielectric);  // dielectric and dielectric or metal?
	ScintWrapper->SetFinish(polishedfrontpainted);  // teflon wrapping on polished surface->front/back painted // Teflon should be Lambertian in air, specular in optical grease
	//ScintWrapper->SetPolish(0.9);  // specfic to the glisur model

	G4MaterialPropertiesTable * ScintWrapperMPT = new G4MaterialPropertiesTable();
	G4double rIndex_Teflon[numShort] = {1.35, 1.35, 1.35}; //Taken from wikipedia
	ScintWrapperMPT->AddProperty("RINDEX", photonEnergyShort, rIndex_Teflon, nEntriesShort)->SetSpline(true);  //refractive index can change with energy
	G4double reflectivity[numShort] = {0.95, 0.95, 0.95};
	ScintWrapperMPT->AddProperty("REFLECTIVITY", photonEnergyShort, reflectivity, nEntriesShort)->SetSpline(true);  //
	ScintWrapper->SetMaterialPropertiesTable(ScintWrapperMPT);
	ScintWrapper->DumpInfo();


	//////// Quartz  ////////////
	G4MaterialPropertiesTable * QuartzMPT = new G4MaterialPropertiesTable();
	G4double rIndex_Quartz[numShort] = {1.474, 1.474, 1.474}; //Taken from Joey github
	QuartzMPT->AddProperty("RINDEX", photonEnergyShort, rIndex_Quartz, nEntriesShort)->SetSpline(true);  //refractive index can change with energy
	QuartzMPT->AddConstProperty("ABSLENGTH", 40.*cm); //from Joeys github
	PMTG4material->SetMaterialPropertiesTable(QuartzMPT);


	///// Building the Plastic Geometry /////

	//Opening angle of DESCANT is 65.5 degrees, approx 1.143 radians. Limiting to 1.13, but 1.23 helps for making horizontal cuts for PMTs
	G4double startTheta = 0.;
	//G4double deltaTheta = 1.13;
	G4double deltaTheta = 1.23;
	G4double deltaThetaWrap = 1.23;
	G4double startPhi = 0.;
	G4double deltaPhi = 2*M_PI;
	G4double startThetaPMT = 0.;
	G4double deltaThetaPMT = 1.23;
	G4double startPhiPMT1 = 0.;
	G4double deltaPhiPMT1 = M_PI;
	G4double startPhiPMT2 = M_PI;
	G4double deltaPhiPMT2 = M_PI;

	//Detector Placement
	//place outer radius of plastics at position of DESCANT detectors, taking into account lead shield
	//G4double outerRadius = fRadialDistance - fLeadShieldThickness - fSpacing;
	//place outer radius of plastics at position of DESCANT detectors, without lead shield
	G4double outerRadius = fRadialDistance - fSpacing;
	G4double innerRadius = outerRadius - fScintillatorWidth;

	//Optical Wrap dimensions
	G4double outerWrap = outerRadius + fWrapThickness+fAirGap;
	G4double innerWrap = innerRadius - fWrapThickness-fAirGap;
	//G4double outerWrap = outerRadius + 10.*(fWrapThickness)+fAirGap; //factor of 10 for visualization purposes
	//G4double innerWrap = innerRadius - 10.*(fWrapThickness)-fAirGap;//factor of 10 for visualization purposes 

	//To determine where the detectors start (left side ie +ve x-axis)
	G4double length = 80.*cm; //at maximum should be 2*outerRadius*sin(deltaTheta)?
	G4double detWidth = (length-2.*(fWrapThickness+fAirGap)*fNumDet)/fNumDet;
	G4double startPos = length/2.-detWidth/2.-fWrapThickness-fAirGap;
	G4double wrapBoxThick = (detWidth+2.*fWrapThickness+2.*fAirGap)/2.;
	G4cout << "detWidth: " << detWidth <<G4endl;

	//For placing volume
	move = G4ThreeVector(0., 0., 0.);
	rotate = new G4RotationMatrix;
	//For intersection solid
	G4RotationMatrix *rot1 = new G4RotationMatrix(0,0,0);

	//G4Solids
	G4double pmtSize = 4.5*cm;
	//G4double pmtSize = 3.5*cm;
	G4double YPosPMT1 = length/2.;
	G4VSolid * boxBars = new G4Box("boxBars", detWidth/2., 1.*m , 1.*m);
	G4VSolid * boxBarsPMTSmaller = new G4Box("boxBarsPMTSmaller", detWidth/2.-0.1*mm, 1.*m , 1.*m);
	G4VSolid * boxBarsPMT = new G4Box("boxBarsPMT", detWidth/2.+fAirGap+0.01*cm, 2.*pmtSize , 1.*m);
	G4VSolid * boxBarsLarger = new G4Box("boxBarsLarger", detWidth/2.+fAirGap, 1.*m , 1.*m);
	G4VSolid * SphereBars = new G4Sphere("SphereBars", innerRadius, outerRadius, startPhi, deltaPhi, startTheta, deltaTheta);
	G4IntersectionSolid * interSolidBars;
	G4IntersectionSolid * interSolidBarsLarger;
	G4VSolid * SphereWrap = new G4Sphere("SphereWrap", innerWrap, outerWrap, startPhi, deltaPhi, startTheta, deltaThetaWrap);
	G4VSolid * boxWrap = new G4Box("boxWrap", wrapBoxThick, 1.*m , 1.*m);
	G4IntersectionSolid * interSolidWrap;
	G4SubtractionSolid * subtractSolidWrap;
	G4VSolid * SphereInnerWrap = new G4Sphere("SphereInnerWrap", innerRadius-fAirGap, outerRadius+fAirGap, startPhi, deltaPhi, startTheta, deltaThetaWrap); //changed from innerRadius->Wrap
	G4VSolid * SpherePMT1 = new G4Sphere("SpherePMT1", innerRadius+0.1*mm, outerRadius-0.1*mm, startPhiPMT1, deltaPhiPMT1, startThetaPMT, deltaThetaPMT); //changed from innerRadius->Wrap
	G4VSolid * SpherePMT2 = new G4Sphere("SpherePMT2", innerRadius+0.1*mm, outerRadius-0.1*mm, startPhiPMT2, deltaPhiPMT2, startThetaPMT, deltaThetaPMT);//changed from innerRadius->Wrap
	G4double BeamLineXY = 6.5*cm;
	G4SubtractionSolid * subtractSolidBeamLine_temp;
	G4SubtractionSolid * subtractSolidBeamLine_Bar;
	G4SubtractionSolid * subtractSolidBeamLine_Wrap;
	G4SubtractionSolid * subtractSolidBeamLine_PMT;
	G4SubtractionSolid * bars_PMT1_Corrected;
	G4SubtractionSolid * bars_PMT_Corrected;
	G4IntersectionSolid * interSolidBarsPMT1;
	G4IntersectionSolid * interSolidBarsPMT2;
	G4SubtractionSolid * PMT1;
	G4SubtractionSolid * PMT2;

	//Set visual attributes
	G4VisAttributes * plastic_vis = new G4VisAttributes(silver);
	plastic_vis->SetVisibility(true);
	G4VisAttributes * wrap_vis = new G4VisAttributes(black);
	wrap_vis->SetVisibility(true);
	G4VisAttributes * pmt_vis = new G4VisAttributes(bronze);
	pmt_vis->SetVisibility(true);

	//Counter for bottom detectors 
	G4int bCounter=0;
	G4double bStartPos=0.;
	G4double startPosY = -1.*m+BeamLineXY;

	//Build array of logical volumes including top detectors above beam line
	for (int i= 0 ; i < fNumDet ; i++) {
		//Determine where to make horizontal cuts on detectors to allow for PMTs.  This changes with x-position
		if (startPos>0.){
			YPosPMT1 = sqrt(pow(innerRadius, 2.0) - pow(startPos+detWidth/2.,2.0));
		}
		else if (startPos<0.){
			YPosPMT1 = sqrt(pow(innerRadius, 2.0) - pow(startPos-detWidth/2.,2.0));
		}
		else if (startPos == 0) {
			YPosPMT1 = sqrt(pow(innerRadius, 2.0) - pow(startPos-detWidth/2.,2.0));
		}
		//For Plastics
		interSolidBars = new G4IntersectionSolid("interSolidBars", SphereBars, boxBars, rot1, G4ThreeVector(startPos, 0, 50*cm));
		bars_PMT1_Corrected = new G4SubtractionSolid("bars_PMT1_Corrected", interSolidBars, boxBarsPMT, rot1, G4ThreeVector(startPos, YPosPMT1 , 50*cm));
		bars_PMT_Corrected = new G4SubtractionSolid("bars_PMT_Corrected", bars_PMT1_Corrected, boxBarsPMT, rot1, G4ThreeVector(startPos, -YPosPMT1 , 50*cm));
		//For PMT's
		interSolidBarsPMT1 = new G4IntersectionSolid("interSolidBarsPMT1", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, 0, 50*cm));
		interSolidBarsPMT2 = new G4IntersectionSolid("interSolidBarsPMT2", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, 0, 50*cm));
		PMT1 = new G4SubtractionSolid("PMT1", interSolidBarsPMT1, bars_PMT_Corrected, rot1, G4ThreeVector(0., 0., 0.));
		PMT2 = new G4SubtractionSolid("PMT2", interSolidBarsPMT2, bars_PMT_Corrected, rot1, G4ThreeVector(0., 0., 0.));

		//For wrapping
		interSolidBarsLarger = new G4IntersectionSolid("interSolidBarsLarger", SphereInnerWrap, boxBarsLarger, rot1, G4ThreeVector(startPos, 0, 50*cm));
		interSolidWrap = new G4IntersectionSolid("interSolidWrap", SphereWrap, boxWrap, rot1, G4ThreeVector(startPos, 0, 50*cm));
		subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, interSolidBarsLarger, rot1, G4ThreeVector(0,0,0));

		//Names
		G4String name0 = "PlasticDet_";
		G4String nameLog=name0+std::to_string(i);

		G4String name1 = "wrapper_";
		G4String nameWrapper = name1+std::to_string(i);

		G4String name2 = "PMT1_top_";
		G4String namePMT1 = name2+std::to_string(i);

		G4String name3 = "PMT2_bottom_";
		G4String namePMT2 = name3+std::to_string(i);

		//G4cout << "name: " << name << G4endl;
		G4cout << "startPos: " << startPos << G4endl;

		//Putting in subtraction Beam Line 
		if( abs(startPos-detWidth/2.) < BeamLineXY || abs(startPos+detWidth/2.) < BeamLineXY || abs(startPos) < BeamLineXY) {

			//increment bottom counter and get the final position for use with below beamline detectors
			bCounter++;
			bStartPos = startPos;
			G4cout << "startPos In loop: " << startPos << G4endl;
			//Implemented in a way where no partial subtractions occur.  change wrapBoxThick<->BeamLineXY in x position to revert to partial subtractions.
			//G4VSolid * boxBeamLine = new G4Box("boxBeamLine", wrapBoxThick, BeamLineXY, 1*m); 
			//Only make top detector
			G4double boxheight = 1.*m;
			//		G4double pmtSize = 4.5*cm;
			G4double subtractHeight = 2*boxheight + startPosY + pmtSize;
			G4double subtractHeight2 = startPosY + pmtSize;
			G4VSolid * boxBeamLine = new G4Box("boxBeamLine", wrapBoxThick+0.5*mm, boxheight, 1*m); 
			//Implemented in a way where no partial subtractions occur.  change startPos<->0. in x position to revert to partial subtractions.
			subtractSolidBeamLine_temp = new G4SubtractionSolid("subtractSolidBeamLine_temp", bars_PMT_Corrected, boxBeamLine, rot1, G4ThreeVector(startPos,subtractHeight,50*cm));
			subtractSolidBeamLine_Bar = new G4SubtractionSolid("subtractSolidBeamLine_Bar", bars_PMT_Corrected, boxBeamLine, rot1, G4ThreeVector(startPos,subtractHeight2,50*cm));
			subtractSolidBeamLine_Wrap = new G4SubtractionSolid("subtractSolidBeamLine_Wrap", subtractSolidWrap, boxBeamLine, rot1, G4ThreeVector(startPos,startPosY,50*cm));
			//Make top inner PMT
			subtractSolidBeamLine_PMT = new G4SubtractionSolid("subtractSolidBeamLine_PMT", subtractSolidBeamLine_temp, boxBeamLine, rot1, G4ThreeVector(startPos,startPosY,50*cm));
			//Assign Logical Volume for detectors and wrapping affected by beamline
			fPlasticLogArray[i] = new G4LogicalVolume(subtractSolidBeamLine_Bar, plasticG4material, nameLog,0,0,0);
			fWrapLogArray[i] = new G4LogicalVolume(subtractSolidBeamLine_Wrap, wrapG4material, nameWrapper,0,0,0); 
			fPMT1LogArray[i] = new G4LogicalVolume(PMT1, PMTG4material, namePMT1,0,0,0);
			fPMT2LogArray[i] = new G4LogicalVolume(subtractSolidBeamLine_PMT, PMTG4material, namePMT2,0,0,0);
		}

		else 
		{
			//Assign Logical Volume for detectors and wrapping unaffected by beamline
			fPlasticLogArray[i] = new G4LogicalVolume(bars_PMT_Corrected, plasticG4material, nameLog,0,0,0);
			fWrapLogArray[i] = new G4LogicalVolume(subtractSolidWrap, wrapG4material, nameWrapper,0,0,0);
			fPMT1LogArray[i] = new G4LogicalVolume(PMT1, PMTG4material, namePMT1,0,0,0);
			fPMT2LogArray[i] = new G4LogicalVolume(PMT2, PMTG4material, namePMT2,0,0,0);
		}

		//Set Logical Skin for optical photons on wrapping
		G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapper, fWrapLogArray[i], ScintWrapper);
		//Give everything colour
		fPlasticLogArray[i]->SetVisAttributes(plastic_vis);
		fWrapLogArray[i]->SetVisAttributes(wrap_vis);
		fPMT1LogArray[i]->SetVisAttributes(pmt_vis);
		fPMT2LogArray[i]->SetVisAttributes(pmt_vis);

		//build every second detector
		//if(i%2==0) {}
		//Place detectors
		fAssemblyPlastics->AddPlacedVolume(fPlasticLogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT1LogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT2LogArray[i], move, rotate);
		if(fAddWrap ==true){
			fAssemblyPlastics->AddPlacedVolume(fWrapLogArray[i], move, rotate);
		}
		//Shift to the right ie all detectors except bottom detectors are built left to right (+ve x to -ve)
		startPos = startPos-detWidth-2.*fWrapThickness- 2.*fAirGap;
		//G4cout << "startPos after : " << startPos << G4endl;

	}
	// Assign variables for bottom detectors
	fNumDetBot = bCounter;
	startPos = bStartPos;
	startPosY = 1.*m-BeamLineXY;
	//Resize Logical Volume Arrays
	fPlasticLogArray.resize(fNumDetBot+fNumDet, NULL);
	fWrapLogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMT1LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMT2LogArray.resize(fNumDetBot+fNumDet, NULL);


	//Loop for building bottom detectors
	for(G4int k=0; k<fNumDetBot; ++k) {
		G4cout << "startPos In Bottom loop: " << startPos << G4endl;
		//Naming
		G4int detNumBottom = fNumDet + k;
		G4String name0 = "PlasticDet_";
		G4String nameLog=name0+std::to_string(detNumBottom);
		G4String name1 = "wrapper_";
		G4String nameWrapper = name1+std::to_string(detNumBottom);
		G4String name2 = "PMT1_top_";
		G4String namePMT1 = name2+std::to_string(detNumBottom);
		G4String name3 = "PMT2_bottom_";
		G4String namePMT2 = name3+std::to_string(detNumBottom);

		if (startPos>0.){
			YPosPMT1 = sqrt(pow(innerRadius, 2.0) - pow(startPos+detWidth/2.,2.0));
		}
		else if (startPos<0.){
			YPosPMT1 = sqrt(pow(innerRadius, 2.0) - pow(startPos-detWidth/2.,2.0));
		}
		else if (startPos == 0) {
			YPosPMT1 = sqrt(pow(innerRadius, 2.0) - pow(startPos-detWidth/2.,2.0));
		}
		//For Plastics
		interSolidBars = new G4IntersectionSolid("interSolidBars", SphereBars, boxBars, rot1, G4ThreeVector(startPos, 0, 50*cm));
		bars_PMT1_Corrected = new G4SubtractionSolid("bars_PMT1_Corrected", interSolidBars, boxBarsPMT, rot1, G4ThreeVector(startPos, YPosPMT1 , 50*cm));
		bars_PMT_Corrected = new G4SubtractionSolid("bars_PMT_Corrected", bars_PMT1_Corrected, boxBarsPMT, rot1, G4ThreeVector(startPos, -YPosPMT1 , 50*cm));
		//For PMT's
		interSolidBarsPMT1 = new G4IntersectionSolid("interSolidBarsPMT1", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, 0, 50*cm));
		interSolidBarsPMT2 = new G4IntersectionSolid("interSolidBarsPMT2", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, 0, 50*cm));
		PMT1 = new G4SubtractionSolid("PMT1", interSolidBarsPMT1, bars_PMT_Corrected, rot1, G4ThreeVector(0., 0., 0.));
		PMT2 = new G4SubtractionSolid("PMT2", interSolidBarsPMT2, bars_PMT_Corrected, rot1, G4ThreeVector(0., 0., 0.));
		//For wrapping
		interSolidBarsLarger = new G4IntersectionSolid("interSolidBarsLarger", SphereInnerWrap, boxBarsLarger, rot1, G4ThreeVector(startPos, 0, 50*cm));
		interSolidWrap = new G4IntersectionSolid("interSolidWrap", SphereWrap, boxWrap, rot1, G4ThreeVector(startPos, 0, 50*cm));
		subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, interSolidBarsLarger, rot1, G4ThreeVector(0,0,0));

		G4double boxheight = 1.*m;
		//G4double pmtSize = 4.5*cm;
		G4double subtractHeight = -2*boxheight + startPosY - pmtSize;
		G4double subtractHeight2 = startPosY - pmtSize;
		G4VSolid * boxBeamLine = new G4Box("boxBeamLine", wrapBoxThick+0.5*mm, boxheight, 1*m); 
		//Implemented in a way where no partial subtractions occur.  change startPos<->0. in x position to revert to partial subtractions.
		subtractSolidBeamLine_temp = new G4SubtractionSolid("subtractSolidBeamLine_temp", bars_PMT_Corrected, boxBeamLine, rot1, G4ThreeVector(startPos,subtractHeight,50*cm));
		subtractSolidBeamLine_Bar = new G4SubtractionSolid("subtractSolidBeamLine_Bar", bars_PMT_Corrected, boxBeamLine, rot1, G4ThreeVector(startPos,subtractHeight2,50*cm));
		subtractSolidBeamLine_Wrap = new G4SubtractionSolid("subtractSolidBeamLine_Wrap", subtractSolidWrap, boxBeamLine, rot1, G4ThreeVector(startPos,startPosY,50*cm));
		//Make top inner PMT
		subtractSolidBeamLine_PMT = new G4SubtractionSolid("subtractSolidBeamLine_PMT", subtractSolidBeamLine_temp, boxBeamLine, rot1, G4ThreeVector(startPos,startPosY,50*cm));
		//Assign Logical Volume for detectors and wrapping affected by beamline
		fPlasticLogArray[detNumBottom] = new G4LogicalVolume(subtractSolidBeamLine_Bar, plasticG4material, nameLog,0,0,0);
		fWrapLogArray[detNumBottom] = new G4LogicalVolume(subtractSolidBeamLine_Wrap, wrapG4material, nameWrapper,0,0,0); 
		fPMT1LogArray[detNumBottom] = new G4LogicalVolume(subtractSolidBeamLine_PMT, PMTG4material, namePMT1,0,0,0);
		fPMT2LogArray[detNumBottom] = new G4LogicalVolume(PMT2, PMTG4material, namePMT2,0,0,0);
		//Set Logical Skin for optical photons on wrapping
		G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapper, fWrapLogArray[detNumBottom], ScintWrapper);
		//Give everything colour
		fPlasticLogArray[detNumBottom]->SetVisAttributes(plastic_vis);
		fWrapLogArray[detNumBottom]->SetVisAttributes(wrap_vis);
		fPMT1LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMT2LogArray[detNumBottom]->SetVisAttributes(pmt_vis);

		//Place detectors
		fAssemblyPlastics->AddPlacedVolume(fPlasticLogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT1LogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT2LogArray[detNumBottom], move, rotate);
		if(fAddWrap ==true){
			fAssemblyPlastics->AddPlacedVolume(fWrapLogArray[detNumBottom], move, rotate);
		}
		startPos = startPos+detWidth+2.*fWrapThickness+ 2.*fAirGap;


	}




	return 1;
}

