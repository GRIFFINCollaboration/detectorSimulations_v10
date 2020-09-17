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
	fAddFrontPMTs = false;
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
	fPMTFace1LogArray.resize(numDet, NULL);
	fPMTFace2LogArray.resize(numDet, NULL);
	fPMTFaceMidLogArray.resize(numDet, NULL);

	fWrapThickness = 5. * mm; //Factor of 10 should be  applied later for visualization purposes.  0.05 for simulations, 0.5 for visualization
	fSpacing = fWrapThickness; //assuming no lead on DESCANT but taking into account the optical wrapping
	//fWrapThickness = 0.1 * cm;
	//fAirGap = 0.1 * fWrapThickness;
	fAirGap = 0.000001*m;
	fWrapMaterial = "Teflon";
	fPMTMaterial = "G4_SILICON_DIOXIDE";

	//blue=G4Color(0.,0.,1.);
	bronze=G4Color(0.8,0.5,0.2);
	//cyan=G4Color(0.,1.,1.);
	silver=G4Color(0.75,0.75,0.75);
	black=G4Color(0.,0.,0.);

	if(material == 1)  fPlasticMaterial = "BC404";
	else if (material == 2) fPlasticMaterial = "BC408";
	else if (material == 3) fPlasticMaterial = "deuterium";
	else if (material == 4) fPlasticMaterial = "Hydrogen";
	else if (material == 5) fPlasticMaterial = "Carbon";
	else if (material == 6) fPlasticMaterial = "Deuterated Scintillator";
	else if (material == 7) fPlasticMaterial = "Dense";
	else if (material == 8) fPlasticMaterial = "BC537";
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
		delete fPMTFace1LogArray[i];
		delete fPMTFace2LogArray[i];
		delete fPMTFaceMidLogArray[i];
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
	///////
	//

	const G4int num = 20;

	//G4double photonEnergy[num] = {1.7*eV, 2.38*eV, 2.48*eV, 2.58*eV, 2.70*eV, 2.76*eV, 2.82*eV, 2.91*eV, 2.95*eV, 3.1*eV, 3.26*eV, 3.44*eV}; //BC408 emission spectra & corresponding energies
	G4double photonEnergy[num] = {1.7*eV, 2.38*eV, 2.48*eV, 2.58*eV, 2.70*eV, 2.76*eV, 2.82*eV, 2.91*eV, 2.95*eV, 2.97*eV, 3.0*eV, 3.02*eV, 3.04*eV,  3.06*eV, 3.1*eV, 3.14*eV, 3.18*eV, 3.21*eV, 3.26*eV, 3.44*eV}; //BC404 emission spectra & corresponding energies

	G4double RIndex1[num] = {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58,  1.58, 1.58, 1.58};
	assert(sizeof(RIndex1) == sizeof(photonEnergy));
	const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

	G4double absorption[num] = {160.*cm, 160*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm, 160.*cm}; ///light attenuation
	//	G4double absorption[num] = {4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm, 4000.*cm}; ///light attenuation
	assert(sizeof(absorption) == sizeof(photonEnergy));

	//G4double scint[num] = {3., 3., 8., 18., 43., 55., 80., 100., 80., 20., 7., 3. }; ///// Based off emission spectra for BC408
	//G4double scint[num] = {1., 1., 2., 5., 13., 20., 35., 50., 55., 60., 85., 93., 100., 96., 87., 70., 38., 18., 5., 1. }; ///// Based off emission spectra for BC404
	G4double scint[num] = {0., 1., 2., 5., 13., 20., 35., 50., 55., 60., 85., 93., 100., 96., 87., 70., 38., 18., 5., 1. }; ///// Based off emission spectra for BC404
	//G4double scint[num] = {0.0011976048, 0.0011976048, 0.0023952096, 0.0059880240, 0.015568862, 0.023952096, 0.041916168, 0.059880240, 0.065868263, 0.071856287, 0.10179641, 0.11137725, 0.11976048, 0.11497006, 0.10419162, 0.083832335, 0.045508982, 0.021556886, 0.0059880240, 0.0011976048 }; ///// Based off emission spectra for BC404 . Normalized (total 835)
	assert(sizeof(scint) == sizeof(photonEnergy));


	scintillatorMPT->AddProperty("FASTCOMPONENT", photonEnergy, scint, nEntries)->SetSpline(true); // BC408 emission spectra
	scintillatorMPT->AddProperty("SLOWCOMPONENT", photonEnergy, scint, nEntries)->SetSpline(true); // BC408 emission spectra
	scintillatorMPT->AddProperty("RINDEX", photonEnergy, RIndex1, nEntries);  //refractive index can change with energy
	//note if photon is created outside of energy range it will have no index of refraction
	scintillatorMPT->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)->SetSpline(true); //absorption length doesnt change with energy - examples showing it can...
	//scintillatorMPT->AddConstProperty("ABSLENGTH", 160.*cm); //Bulk light attenuation - 380 for 408
	//scintillatorMPT->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV); //Scintillation Efficiency - characteristic light yield //10000./MeV

	// The number of photons produced per interaction is sampled from a Gaussian distribution with a full-width at half-maximum set to 20% of the number of produced photons. From Joeys Thesis
	scintillatorMPT->AddConstProperty("RESOLUTIONSCALE", 1.2); // broadens the statistical distribution of generated photons, sqrt(num generated)* resScale, gaussian based on SCINTILLATIONYIELD, >1 broadens, 0 no distribution. 20%
	scintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 1.8*ns); //only one decay constant given - 2.1ns for 408, 1.8ns for 404
	scintillatorMPT->AddConstProperty("SLOWTIMECONSTANT", 1.8*ns); //only one decay constant given - triplet-triplet annihilation 
	//	scintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 0.00000000000000000*ns); // for testing the effective speed of light
	//	scintillatorMPT->AddConstProperty("SLOWTIMECONSTANT", 0.00000000000000000*ns); // for testing the effective speed of light
	scintillatorMPT->AddConstProperty("YIELDRATIO", 1.0); //The relative strength of the fast component as a fraction of total scintillation yield is given by the YIELDRATIO.
	//Should these be in the physics list?
	G4OpticalPhysics * opticalPhysics = new G4OpticalPhysics();
	opticalPhysics->SetFiniteRiseTime(true);
	scintillatorMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.7*ns); //default rise time is 0ns, have to set manually
	scintillatorMPT->AddConstProperty("SLOWSCINTILLATIONRISETIME", 0.7*ns); //default rise time is 0ns, have to set manually
	//	scintillatorMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.*ns); // For testing speed of light
	//	scintillatorMPT->AddConstProperty("SLOWSCINTILLATIONRISETIME", 0.*ns); // For testing speed of light
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


	///// Building the Plastic Geometry /////

	//Opening angle of DESCANT is 65.5 degrees, approx 1.143 radians. Limiting to 1.13, but 1.23 helps for making horizontal cuts for PMTs
	G4double realAngle = 1.13;
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
	G4double zStartPos = innerRadius * cos (realAngle);
	G4cout << "zStartPos: " << zStartPos <<G4endl;

	//For placing volume
	move = G4ThreeVector(0., 0., 0.);
	rotate = new G4RotationMatrix;
	//For intersection solid
	G4RotationMatrix *rot1 = new G4RotationMatrix(0,0,0);

	//G4Solids
	G4double pmtSize = 4.5*cm;
	//G4double pmtSize = 3.5*cm;
	G4double pmtThick = 0.6*cm;
	G4double calcRadii = (outerRadius - innerRadius - pmtThick)/2.;
	//G4double pmt_width = 0.5*cm;
	G4double pmt_width = pmtThick;
	G4double BeamLineXY = 6.5*cm;
	//G4VSolid * boxBars = new G4Box("boxBars", detWidth/2., 1.*m , 1.*m); //old geometry, new is inside loop currently
	//G4VSolid * boxBarsPMTSmaller = new G4Box("boxBarsPMTSmaller", detWidth/2.-0.1*mm, 1.*m , 1.*m); //old geometry, new is inside loop
	//G4VSolid * boxBarsPMT = new G4Box("boxBarsPMT", detWidth/2.+fAirGap+0.01*cm, 2.*pmtSize , 1.*m);
	//G4VSolid * boxBarsPMT = new G4Box("boxBarsPMT", detWidth/2.+fAirGap+0.01*cm, pmtSize , 1.*m);
	//G4VSolid * boxBarsLarger = new G4Box("boxBarsLarger", detWidth/2.+fAirGap, 1.*m , 1.*m);
	G4VSolid * SphereBars = new G4Sphere("SphereBars", innerRadius, outerRadius, startPhi, deltaPhi, startTheta, deltaTheta);
	G4IntersectionSolid * interSolidBars;
	G4IntersectionSolid * interSolidBarsLarger;
	G4VSolid * SphereWrap = new G4Sphere("SphereWrap", innerWrap, outerWrap, startPhi, deltaPhi, startTheta, deltaThetaWrap);
	//G4VSolid * boxWrap = new G4Box("boxWrap", wrapBoxThick, 1.*m , 1.*m); //old geomtery, new is inside loop
	G4IntersectionSolid * interSolidWrap;
	G4SubtractionSolid * subtractSolidWrap;
	G4VSolid * SphereInnerWrap = new G4Sphere("SphereInnerWrap", innerRadius-fAirGap, outerRadius+fAirGap, startPhi, deltaPhi, startTheta, deltaThetaWrap); //changed from innerRadius->Wrap
	G4VSolid * SpherePMT_Face1 = new G4Sphere("SpherePMT_Face1", innerRadius-pmt_width, innerRadius, startPhi, deltaPhi, startTheta, deltaTheta);
	G4VSolid * SpherePMT_FaceLarger = new G4Sphere("SpherePMT_FaceLarger", innerRadius-pmt_width, innerRadius+fWrapThickness, startPhi, deltaPhi, startTheta, deltaTheta);
	G4VSolid * SpherePMT1 = new G4Sphere("SpherePMT1", innerRadius+calcRadii, outerRadius-calcRadii, startPhiPMT1, deltaPhiPMT1, startThetaPMT, deltaThetaPMT); //changed from innerRadius->Wrap
	G4VSolid * SpherePMT2 = new G4Sphere("SpherePMT2", innerRadius+calcRadii, outerRadius-calcRadii, startPhiPMT2, deltaPhiPMT2, startThetaPMT, deltaThetaPMT);//changed from innerRadius->Wrap
	//G4VSolid * SpherePMT1 = new G4Sphere("SpherePMT1", innerRadius+0.1*mm, outerRadius-0.1*mm, startPhiPMT1, deltaPhiPMT1, startThetaPMT, deltaThetaPMT); //changed from innerRadius->Wrap
	//G4VSolid * SpherePMT2 = new G4Sphere("SpherePMT2", innerRadius+0.1*mm, outerRadius-0.1*mm, startPhiPMT2, deltaPhiPMT2, startThetaPMT, deltaThetaPMT);//changed from innerRadius->Wrap
	//G4IntersectionSolid * interSolidBarsPMT1;
	//G4IntersectionSolid * interSolidBarsPMT2;
	G4IntersectionSolid * PMT1;
	G4IntersectionSolid * PMT2;
	G4IntersectionSolid * PMT2BeamLine;
	G4IntersectionSolid * PMT1BeamLine;
	//G4SubtractionSolid * PMT1; //old
	//G4SubtractionSolid * PMT2; //old
	G4UnionSolid * Bar_PMT1;
	G4UnionSolid * outsidePMT_bar;	
	G4UnionSolid * outsidePMT_bar_Face1;	
	G4UnionSolid * outsidePMT_bar_AllFace;	
	G4IntersectionSolid * PMT_Front1;
	G4IntersectionSolid * PMT_Front2;
	G4IntersectionSolid * PMT_FrontMid;
	G4IntersectionSolid * PMT_FrontLarger1;
	G4IntersectionSolid * PMT_FrontLarger2;
	G4IntersectionSolid * PMT_FrontLargerMid;
	G4IntersectionSolid * PMT1Larger;
	G4IntersectionSolid * PMT2Larger;


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

	//Build array of logical volumes including top detectors above beam line
	for (int i= 0 ; i < fNumDet ; i++) {
		//Names
		G4String name0 = "PlasticDet_";
		G4String nameLog=name0+std::to_string(i);
		G4String name1 = "wrapper_";
		G4String nameWrapper = name1+std::to_string(i);
		G4String name2 = "PMT1_top_";
		G4String namePMT1 = name2+std::to_string(i);
		G4String name3 = "PMT2_bottom_";
		G4String namePMT2 = name3+std::to_string(i);
		G4String name4 = "PMTFront1_top_";
		G4String namePMTFront1 = name4+std::to_string(i);
		G4String name5 = "PMTFront2_bottom_";
		G4String namePMTFront2 = name5+std::to_string(i);
		G4String name6 = "PMTFront_mid_";
		G4String namePMTFrontMid = name6+std::to_string(i);

		//For new geometry placement of the outside PMTs ie the vertical bar cutoff
		G4double YFacePMT1 = sqrt(pow(innerRadius, 2.0) - pow(abs(startPos)+detWidth/2.,2.0) - pow(zStartPos, 2.));

		//For Plastics New
		G4VSolid * boxBars = new G4Box("boxBars", detWidth/2., YFacePMT1 , 1.*m);
		interSolidBars = new G4IntersectionSolid("interSolidBars", SphereBars, boxBars, rot1, G4ThreeVector(startPos, 0, 50*cm));

		//For PMT's Top & Bottom New
		G4VSolid * boxBarsPMTSmaller = new G4Box("boxBarsPMTSmaller", detWidth/2.-0.1*mm, pmt_width/2. , 1.*m);
		PMT1 = new G4IntersectionSolid("PMT1", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, YFacePMT1+pmt_width/2., 50*cm));
		PMT2 = new G4IntersectionSolid("PMT2", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));

		//Calc beamline middle face pmt placement above/below beamline
		G4double y1=YFacePMT1, y2=(BeamLineXY+pmtSize) , z1 = zStartPos, z2, y3, z3;
		G4double r1 = sqrt(y1*y1 + z1*z1);
		z2 = sqrt(r1*r1 - y2*y2);
		y3 = (y1+y2)/2.;
		z3 = (z1+z2)/2.;
		G4double r3 = sqrt(y3*y3 + z3*z3);
		G4double BeamlineY = y3*r1/r3;

		//Face PMTs, placing them along the vertical at -1/3 and +1/3
		G4VSolid * boxBarsPMTFace = new G4Box("boxBarsPMTFace", detWidth/2.-0.1*mm, pmtThick/2. , 1.*m);
		PMT_Front1 = new G4IntersectionSolid("PMT_Front1", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, YFacePMT1/3., 50*cm));
		PMT_Front2 = new G4IntersectionSolid("PMT_Front2", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -YFacePMT1/3., 50.*cm));

		//Face PMTs, placing the mid pmts  the vertical at 0 for first and last, or the middle wrt to the beamline if above or below
		if (i == 0 || i == (fNumDet-1) ) {
			PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, 0., 0.));
		} else{
			//PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 0.)); //Y centered
			PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, BeamlineY, 0.)); //Center of arc
		}


		//For wrapping Made slightly larger versions of all pmts and bar to be subtracted from the optical wrapping
		//Slightly larger Bar
		G4VSolid * boxBarsLarger = new G4Box("boxBarsLarger", detWidth/2.+fAirGap, YFacePMT1+2.*fAirGap , 1.*m);
		interSolidBarsLarger = new G4IntersectionSolid("interSolidBarsLarger", SphereInnerWrap, boxBarsLarger, rot1, G4ThreeVector(startPos, 0, 50*cm));
		//The actual optical wrap that will be subtracted from
		G4VSolid * boxWrap = new G4Box("boxWrap", wrapBoxThick, YFacePMT1+2.*fAirGap+fWrapThickness , 1.*m);
		interSolidWrap = new G4IntersectionSolid("interSolidWrap", SphereWrap, boxWrap, rot1, G4ThreeVector(startPos, 0, 50*cm));
		//Slightly larger PMTs
		G4VSolid * boxBarsPMTSmaller2 = new G4Box("boxBarsPMTSmaller2", detWidth/2.-0.1*mm, pmt_width/2.+fWrapThickness , 1.*m);
		PMT1Larger = new G4IntersectionSolid("PMT1Larger", SpherePMT1, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos, YFacePMT1+pmt_width/2., 50*cm));
		PMT2Larger = new G4IntersectionSolid("PMT2Larger", SpherePMT2, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));

		PMT_FrontLarger1 = new G4IntersectionSolid("PMT_FrontLarger1", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, YFacePMT1/3., 50*cm));
		PMT_FrontLarger2 = new G4IntersectionSolid("PMT_FrontLarger2", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -YFacePMT1/3., 50.*cm));
		if (i == 0 || i == (fNumDet-1) ) {
			PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, 0., 0.));
		} else{
			//PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 0.)); //YCentered
			PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, BeamlineY, 0.)); //Center of Arc
		}


		//Adding all Larger volumes together to be subtracted later
		Bar_PMT1 = new G4UnionSolid("Bar_PMT1", PMT1Larger, interSolidBarsLarger, rot1, G4ThreeVector(0.,0.,0.));
		outsidePMT_bar = new G4UnionSolid("outsidePMT_bar", Bar_PMT1, PMT2Larger, rot1, G4ThreeVector(0.,0.,0.));

		if( abs(startPos-detWidth/2.) < BeamLineXY || abs(startPos+detWidth/2.) < BeamLineXY || abs(startPos) < BeamLineXY || i == 0 || i == (fNumDet-1) ) {
			outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar, PMT_FrontLargerMid, rot1, G4ThreeVector(0.,0.,0.));
		} else{
			outsidePMT_bar_Face1 = new G4UnionSolid("outsidePMT_bar_Face1", outsidePMT_bar, PMT_FrontLarger1, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar_Face1, PMT_FrontLarger2, rot1, G4ThreeVector(0.,0.,0.));
		}

		//Subtraction to create actual optical wrapping shape
		if(fAddFrontPMTs == true) {
			subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar_AllFace, rot1, G4ThreeVector(0,0,0));
		} else subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar, rot1, G4ThreeVector(0,0,0));

		//G4cout << "name: " << name << G4endl;
		G4cout << "startPos: " << startPos << G4endl;
		G4cout << "YPMT1: " << YFacePMT1 << G4endl;
		G4cout << "ZPos: " <<  sqrt(pow(innerRadius, 2.0) - pow(startPos,2.0) - pow(YFacePMT1, 2.)) << G4endl;

		//Putting in subtraction Beam Line 
		if( abs(startPos-detWidth/2.) < BeamLineXY || abs(startPos+detWidth/2.) < BeamLineXY || abs(startPos) < BeamLineXY) {

			//increment bottom counter and get the final position for use with below beamline detectors
			bCounter++;
			bStartPos = startPos;
			G4cout << "startPos In loop: " << startPos << G4endl;
			//Implemented in a way where no partial subtractions occur.  change wrapBoxThick<->BeamLineXY in x position to revert to partial subtractions.
			//For Plastics New
			G4VSolid * boxBars2 = new G4Box("boxBars2", detWidth/2., (YFacePMT1-pmtSize-BeamLineXY)/2. , 1.*m);
			interSolidBars = new G4IntersectionSolid("interSolidBars", SphereBars, boxBars2, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 50*cm));

			//Only Need to only make bottom pmt here, note SpherePMT1 because above y=0
			PMT2BeamLine = new G4IntersectionSolid("PMT2BeamLine", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, pmtSize+BeamLineXY-pmt_width/2., 50*cm));

			//For wrapping Made slightly larger versions of all pmts and bar to be subtracted from the optical wrapping
			//Slightly larger Bar
			G4VSolid * boxBarsLarger2 = new G4Box("boxBarsLarger2", detWidth/2.+fAirGap, (YFacePMT1-pmtSize-BeamLineXY)/2.+2.*fAirGap , 1.*m);
			interSolidBarsLarger = new G4IntersectionSolid("interSolidBarsLarger", SphereInnerWrap, boxBarsLarger2, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 50*cm));

			//The actual optical wrap that will be subtracted from
			G4VSolid * boxWrap2 = new G4Box("boxWrap2", wrapBoxThick, (YFacePMT1-pmtSize-BeamLineXY)/2.+2.*fAirGap+fWrapThickness , 1.*m);
			interSolidWrap = new G4IntersectionSolid("interSolidWrap", SphereWrap, boxWrap2, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 50*cm));

			//Slightly larger PMTs
			G4VSolid * boxBarsPMTSmaller3 = new G4Box("boxBarsPMTSmaller3", detWidth/2.-0.1*mm, pmt_width/2.+fWrapThickness , 1.*m);
			PMT1Larger = new G4IntersectionSolid("PMT1Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, YFacePMT1+pmt_width/2., 50*cm));
			PMT2Larger = new G4IntersectionSolid("PMT2Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, pmtSize+BeamLineXY-pmt_width/2., 50*cm));

			//Adding all Larger volumes together to be subtracted later
			Bar_PMT1 = new G4UnionSolid("Bar_PMT1", PMT1Larger, interSolidBarsLarger, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar = new G4UnionSolid("outsidePMT_bar", Bar_PMT1, PMT2Larger, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar, PMT_FrontLargerMid, rot1, G4ThreeVector(0.,0.,0.));

			//Subtraction to create actual optical wrapping shape
			if(fAddFrontPMTs == true) {
				subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar_AllFace, rot1, G4ThreeVector(0,0,0));
			} else subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar, rot1, G4ThreeVector(0,0,0));

			//Assign Logical Volume for detectors and wrapping affected by beamline
			fPlasticLogArray[i] = new G4LogicalVolume(interSolidBars, plasticG4material, nameLog,0,0,0);
			fWrapLogArray[i] = new G4LogicalVolume(subtractSolidWrap, wrapG4material, nameWrapper,0,0,0); 
			fPMT1LogArray[i] = new G4LogicalVolume(PMT1, PMTG4material, namePMT1,0,0,0);
			fPMT2LogArray[i] = new G4LogicalVolume(PMT2BeamLine, PMTG4material, namePMT2,0,0,0);
			fPMTFace1LogArray[i] = new G4LogicalVolume(PMT_Front1, PMTG4material, namePMTFront1,0,0,0);
			fPMTFace2LogArray[i] = new G4LogicalVolume(PMT_Front2, PMTG4material, namePMTFront2,0,0,0);
			fPMTFaceMidLogArray[i] = new G4LogicalVolume(PMT_FrontMid, PMTG4material, namePMTFrontMid,0,0,0);
		}

		else 
		{
			//Assign Logical Volume for detectors and wrapping unaffected by beamline
			fPlasticLogArray[i] = new G4LogicalVolume(interSolidBars, plasticG4material, nameLog,0,0,0);
			fWrapLogArray[i] = new G4LogicalVolume(subtractSolidWrap, wrapG4material, nameWrapper,0,0,0);
			fPMT1LogArray[i] = new G4LogicalVolume(PMT1, PMTG4material, namePMT1,0,0,0);
			fPMT2LogArray[i] = new G4LogicalVolume(PMT2, PMTG4material, namePMT2,0,0,0);
			fPMTFace1LogArray[i] = new G4LogicalVolume(PMT_Front1, PMTG4material, namePMTFront1,0,0,0);
			fPMTFace2LogArray[i] = new G4LogicalVolume(PMT_Front2, PMTG4material, namePMTFront2,0,0,0);
			fPMTFaceMidLogArray[i] = new G4LogicalVolume(PMT_FrontMid, PMTG4material, namePMTFrontMid,0,0,0);
		}

		//Set Logical Skin for optical photons on wrapping
		G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapper, fWrapLogArray[i], ScintWrapper);

		//Give everything colour
		fPlasticLogArray[i]->SetVisAttributes(plastic_vis);
		fWrapLogArray[i]->SetVisAttributes(wrap_vis);
		fPMT1LogArray[i]->SetVisAttributes(pmt_vis);
		fPMT2LogArray[i]->SetVisAttributes(pmt_vis);
		fPMTFaceMidLogArray[i]->SetVisAttributes(pmt_vis);
		fPMTFace1LogArray[i]->SetVisAttributes(pmt_vis);
		fPMTFace2LogArray[i]->SetVisAttributes(pmt_vis);

		//build every second detector
		//if(i%2==0) {}
		//Place detectors
		fAssemblyPlastics->AddPlacedVolume(fPlasticLogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT1LogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT2LogArray[i], move, rotate);
		if(fAddFrontPMTs == true) {
			if( abs(startPos-detWidth/2.) < BeamLineXY || abs(startPos+detWidth/2.) < BeamLineXY || abs(startPos) < BeamLineXY || i == 0 || i == (fNumDet-1) ) {
				fAssemblyPlastics->AddPlacedVolume(fPMTFaceMidLogArray[i], move, rotate);
			} else {
				fAssemblyPlastics->AddPlacedVolume(fPMTFace1LogArray[i], move, rotate);
				fAssemblyPlastics->AddPlacedVolume(fPMTFace2LogArray[i], move, rotate);
			}
		}
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
	//Resize Logical Volume Arrays
	fPlasticLogArray.resize(fNumDetBot+fNumDet, NULL);
	fWrapLogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMT1LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMT2LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMTFace1LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMTFace2LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMTFaceMidLogArray.resize(fNumDetBot+fNumDet, NULL);


	//Loop for building bottom detectors
	for(G4int k=0; k<fNumDetBot; ++k) {
		G4cout << "startPos In Bottom loop: " << startPos << G4endl;
		//Names
		G4int detNumBottom = fNumDet + k;
		G4String name0 = "PlasticDet_";
		G4String nameLog=name0+std::to_string(detNumBottom);
		G4String name1 = "wrapper_";
		G4String nameWrapper = name1+std::to_string(detNumBottom);
		G4String name2 = "PMT1_top_";
		G4String namePMT1 = name2+std::to_string(detNumBottom);
		G4String name3 = "PMT2_bottom_";
		G4String namePMT2 = name3+std::to_string(detNumBottom);
		G4String name4 = "PMTFront1_top_";
		G4String namePMTFront1 = name4+std::to_string(detNumBottom);
		G4String name5 = "PMTFront2_bottom_";
		G4String namePMTFront2 = name5+std::to_string(detNumBottom);
		G4String name6 = "PMTFront_mid_";
		G4String namePMTFrontMid = name6+std::to_string(detNumBottom);

		//For new geometry placement of the outside PMTs ie the vertical bar cutoff
		G4double YFacePMT1 = sqrt(pow(innerRadius, 2.0) - pow(abs(startPos)+detWidth/2.,2.0) - pow(zStartPos, 2.));

		//For Plastics New
		G4VSolid * boxBars2 = new G4Box("boxBars2", detWidth/2., (YFacePMT1-pmtSize-BeamLineXY)/2. , 1.*m);
		interSolidBars = new G4IntersectionSolid("interSolidBars", SphereBars, boxBars2, rot1, G4ThreeVector(startPos, (-YFacePMT1-pmtSize-BeamLineXY)/2., 50*cm));

		//For PMT's Top & Bottom New
		G4VSolid * boxBarsPMTSmaller = new G4Box("boxBarsPMTSmaller", detWidth/2.-0.1*mm, pmt_width/2. , 1.*m);
		PMT2 = new G4IntersectionSolid("PMT2", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		//Only Need to only make top pmt here, note SpherePMT2 because below y=0
		PMT1BeamLine = new G4IntersectionSolid("PMT1BeamLine", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));

		//Calc beamline middle face pmt placement above/below beamline
		G4double y1=YFacePMT1, y2=(BeamLineXY+pmtSize) , z1 = zStartPos, z2, y3, z3;
		G4double r1 = sqrt(y1*y1 + z1*z1);
		z2 = sqrt(r1*r1 - y2*y2);
		y3 = (y1+y2)/2.;
		z3 = (z1+z2)/2.;
		G4double r3 = sqrt(y3*y3 + z3*z3);
		G4double BeamlineY = y3*r1/r3;

		//Face PMTs, placing them along the vertical at -1/3 and +1/3
		G4VSolid * boxBarsPMTFace = new G4Box("boxBarsPMTFace", detWidth/2.-0.1*mm, pmtThick/2. , 1.*m);
		PMT_Front1 = new G4IntersectionSolid("PMT_Front1", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, YFacePMT1/3., 50*cm));
		PMT_Front2 = new G4IntersectionSolid("PMT_Front2", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -YFacePMT1/3., 50.*cm));

		//Face PMTs, placing the mid pmts  the vertical at 0 for first and last, or the middle wrt to the beamline if above or below
		//PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 0.)); //Y centered
		PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -BeamlineY, 0.)); //Center of arc

		//For wrapping Made slightly larger versions of all pmts and bar to be subtracted from the optical wrapping
		//Slightly larger Bar
		G4VSolid * boxBarsLarger2 = new G4Box("boxBarsLarger2", detWidth/2.+fAirGap, (YFacePMT1-pmtSize-BeamLineXY)/2.+2.*fAirGap , 1.*m);
		interSolidBarsLarger = new G4IntersectionSolid("interSolidBarsLarger", SphereInnerWrap, boxBarsLarger2, rot1, G4ThreeVector(startPos, (-YFacePMT1-pmtSize-BeamLineXY)/2., 50*cm));

		//The actual optical wrap that will be subtracted from
		G4VSolid * boxWrap2 = new G4Box("boxWrap2", wrapBoxThick, (YFacePMT1-pmtSize-BeamLineXY)/2.+2.*fAirGap+fWrapThickness , 1.*m);
		interSolidWrap = new G4IntersectionSolid("interSolidWrap", SphereWrap, boxWrap2, rot1, G4ThreeVector(startPos, (-YFacePMT1-pmtSize-BeamLineXY)/2., 50*cm));

		//Slightly larger PMTs
		G4VSolid * boxBarsPMTSmaller3 = new G4Box("boxBarsPMTSmaller3", detWidth/2.-0.1*mm, pmt_width/2.+fWrapThickness , 1.*m);
		PMT2Larger = new G4IntersectionSolid("PMT2Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		PMT1Larger = new G4IntersectionSolid("PMT1Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));
		//PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 0.)); //YCentered
		PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -BeamlineY, 0.)); //Center of Arc

		//Adding all Larger volumes together to be subtracted later
		Bar_PMT1 = new G4UnionSolid("Bar_PMT1", PMT1Larger, interSolidBarsLarger, rot1, G4ThreeVector(0.,0.,0.));
		outsidePMT_bar = new G4UnionSolid("outsidePMT_bar", Bar_PMT1, PMT2Larger, rot1, G4ThreeVector(0.,0.,0.));
		outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar, PMT_FrontLargerMid, rot1, G4ThreeVector(0.,0.,0.));

		//Subtraction to create actual optical wrapping shape
		if(fAddFrontPMTs == true) {
			subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar_AllFace, rot1, G4ThreeVector(0,0,0));
		} else subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar, rot1, G4ThreeVector(0,0,0));

		//Assign Logical Volume for detectors and wrapping affected by beamline
		fPlasticLogArray[detNumBottom] = new G4LogicalVolume(interSolidBars, plasticG4material, nameLog,0,0,0);
		fWrapLogArray[detNumBottom] = new G4LogicalVolume(subtractSolidWrap, wrapG4material, nameWrapper,0,0,0); 
		fPMT1LogArray[detNumBottom] = new G4LogicalVolume(PMT1BeamLine, PMTG4material, namePMT1,0,0,0);
		fPMT2LogArray[detNumBottom] = new G4LogicalVolume(PMT2, PMTG4material, namePMT2,0,0,0);
		fPMTFace1LogArray[detNumBottom] = new G4LogicalVolume(PMT_Front1, PMTG4material, namePMTFront1,0,0,0);
		fPMTFace2LogArray[detNumBottom] = new G4LogicalVolume(PMT_Front2, PMTG4material, namePMTFront2,0,0,0);
		fPMTFaceMidLogArray[detNumBottom] = new G4LogicalVolume(PMT_FrontMid, PMTG4material, namePMTFrontMid,0,0,0);

		//Set Logical Skin for optical photons on wrapping
		G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapper, fWrapLogArray[detNumBottom], ScintWrapper);

		//Give everything colour
		fPlasticLogArray[detNumBottom]->SetVisAttributes(plastic_vis);
		fWrapLogArray[detNumBottom]->SetVisAttributes(wrap_vis);
		fPMT1LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMT2LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMTFaceMidLogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMTFace1LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMTFace2LogArray[detNumBottom]->SetVisAttributes(pmt_vis);

		//build every second detector
		//if(i%2==0) {}
		//Place detectors
		fAssemblyPlastics->AddPlacedVolume(fPlasticLogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT1LogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT2LogArray[detNumBottom], move, rotate);
		if(fAddFrontPMTs == true) {
			fAssemblyPlastics->AddPlacedVolume(fPMTFaceMidLogArray[detNumBottom], move, rotate);
		}
		if(fAddWrap ==true){
			fAssemblyPlastics->AddPlacedVolume(fWrapLogArray[detNumBottom], move, rotate);
		}


		startPos = startPos+detWidth+2.*fWrapThickness+ 2.*fAirGap;


	}




	return 1;
}

