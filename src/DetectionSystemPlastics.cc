#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Polyhedra.hh"
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
	fAddFrontPMTs = true;
	fAddWrap = true;
	// detector dimensions  properties
	fScintillatorLength        = 6.*cm;
	fScintillatorHeight        = 6.*cm;
	fScintillatorWidth         = thickness;
	fRadialDistance = 50*cm;
	fLeadShieldThickness = 6.35*mm;
	fNumDet = numDet ;
	
//	fPlasticTileLogArray.resize(4*numDet, NULL);
//	fPMTTileLogArray.resize(4*numDet, NULL);
	
	fPlasticLogArray.resize(numDet, NULL);
	fWrapLogArray.resize(numDet, NULL);
	fPMT11LogArray.resize(numDet, NULL);
	fPMT12LogArray.resize(numDet, NULL);
	fPMT13LogArray.resize(numDet, NULL);
	fPMT21LogArray.resize(numDet, NULL);
	fPMT22LogArray.resize(numDet, NULL);
	fPMT23LogArray.resize(numDet, NULL);
	fPMTFace11LogArray.resize(numDet, NULL);
	fPMTFace12LogArray.resize(numDet, NULL);
	fPMTFace21LogArray.resize(numDet, NULL);
	fPMTFace22LogArray.resize(numDet, NULL);
	fPMTFaceMid1LogArray.resize(numDet, NULL);
	fPMTFaceMid2LogArray.resize(numDet, NULL);

	fWrapThickness = 0.5 * mm; //Factor of 10 should be  applied later for visualization purposes.  0.05 for simulations, 0.5 for visualization
	fSpacing = fWrapThickness; //assuming no lead on DESCANT but taking into account the optical wrapping
	//fWrapThickness = 0.1 * cm;
	//fAirGap = 0.1 * fWrapThickness;
	//fAirGap = 0.000001*m;
	fAirGap = 0.;
	fWrapMaterial = "Teflon";
	fPMTMaterial = "G4_SILICON_DIOXIDE";

	//blue=G4Color(0.,0.,1.);
	bronze=G4Color(0.8,0.5,0.2);
	//cyan=G4Color(0.,1.,1.);
	silver=G4Color(0.75,0.75,0.75);
	black=G4Color(0.,0.,0.);

	if(material == 1)  fPlasticMaterial = "BC408";
	else if (material == 2) fPlasticMaterial = "BC404";
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
		delete fPMT11LogArray[i];
		delete fPMT12LogArray[i];
		delete fPMT13LogArray[i];
		delete fPMT21LogArray[i];
		delete fPMT22LogArray[i];
		delete fPMT23LogArray[i];
		delete fPMTFace11LogArray[i];
		delete fPMTFace12LogArray[i];
		delete fPMTFace21LogArray[i];
		delete fPMTFace22LogArray[i];
		delete fPMTFaceMid1LogArray[i];
		delete fPMTFaceMid2LogArray[i];
	}
/*	for (int i = 0; i<4*(fNumDet+fNumDetBot); i++) {
		delete fPlasticTileLogArray[i];
		delete fPMTTileLogArray[i];
	}
*/	G4cout << "Calling Destructor" << G4endl;

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
/*	const G4int num2 = 4;
	G4double e_range[num2] = {1.*keV, 0.1*MeV, 1.*MeV, 10.*MeV};
	G4double yield_e[num2] = {10., 1000., 10000., 100000.};//More realistic
	//G4double yield_e[num2] = {1000., 10000., 10000., 100000.}; //testing
	G4double yield_p[num2] = {1., 65., 1500., 45000.};
	G4double yield_d[num2] = {1., 65., 1500., 45000.};//no data provided, assume same order of magnitude as proton
	G4double yield_t[num2] = {1., 65., 1500., 45000.};//no data provided, assume same order of magnitude as proton
	G4double yield_a[num2] = {1., 20., 200., 14000.};
	G4double yield_C[num2] = {1., 10., 70., 600.};

	scintillatorMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", e_range, yield_e, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("PROTONSCINTILLATIONYIELD", e_range, yield_p, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("DEUTERONSCINTILLATIONYIELD", e_range, yield_d, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("TRITONSCINTILLATIONYIELD", e_range, yield_t, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("ALPHASCINTILLATIONYIELD", e_range, yield_a, num2)->SetSpline(true);
	scintillatorMPT->AddProperty("IONSCINTILLATIONYIELD", e_range, yield_C, num2)->SetSpline(true);
*/	
	int energyPoints = 26;
	G4double pEF = 0.4; //SiPM efficiency (TODO - discuss it)
	//G4double protonScalingFact = 1.35;
	G4double protonScalingFact = 1;
	G4double psF = protonScalingFact;


	//light yield - data taken form V.V. Verbinski et al, Nucl. Instrum. & Meth. 65 (1968) 8-25
	G4double particleEnergy[] = { 0.001*MeV, 0.1*MeV, 0.13*MeV, 0.17*MeV, 
		0.2*MeV, 0.24*MeV, 0.3*MeV, 0.34*MeV, 
		0.4*MeV, 0.48*MeV, 0.6*MeV, 0.72*MeV, 0.84*MeV,
		1.*MeV, 1.3*MeV, 1.7*MeV, 2.*MeV, 
		2.4*MeV, 3.*MeV, 3.4*MeV, 4.*MeV, 4.8*MeV, 
		6.*MeV, 7.2*MeV, 8.4*MeV, 10.*MeV };
	G4double electronYield[] = {1*pEF, 1000*pEF, 1300*pEF, 1700*pEF,
		2000*pEF, 2400*pEF, 3000*pEF, 3400*pEF,
		4000*pEF, 4800*pEF, 6000*pEF, 7200*pEF,
		8400*pEF,10000*pEF, 13000*pEF, 17000*pEF,
		20000*pEF, 24000*pEF, 30000*pEF, 34000*pEF,
		40000*pEF, 48000*pEF, 60000*pEF, 72000*pEF,
		84000*pEF, 100000*pEF };
	G4double protonYield[] = { 0.6*pEF*psF, 67.1*pEF*psF, 88.6*pEF*psF,
		120.7*pEF*psF,
		146.5*pEF*psF, 183.8*pEF*psF, 246*pEF*psF, 290*pEF*psF,
		365*pEF*psF, 483*pEF*psF, 678*pEF*psF, 910*pEF*psF,
		1175*pEF*psF, 1562*pEF*psF, 2385*pEF*psF, 3660*pEF*psF,
		4725*pEF*psF,6250*pEF*psF, 8660*pEF*psF, 10420*pEF*psF,
		13270*pEF*psF,17180*pEF*psF, 23100*pEF*psF,
		29500*pEF*psF, 36200*pEF*psF, 45500*pEF*psF};
	G4double alphaYield[] = { 0.2*pEF, 16.4*pEF, 20.9*pEF, 27.2*pEF,
		32*pEF, 38.6*pEF, 49*pEF, 56.4*pEF,
		67.5*pEF, 83*pEF, 108*pEF, 135*pEF,
		165.6*pEF, 210*pEF, 302*pEF, 441*pEF,
		562*pEF, 750*pEF, 1100*pEF, 1365*pEF,
		1815*pEF, 2555*pEF, 4070*pEF, 6070*pEF,
		8700*pEF, 13200*pEF };

	G4double ionYield[] = { 0.2*pEF, 10.4*pEF, 12.7*pEF, 15.7*pEF,
		17.9*pEF, 20.8*pEF, 25.1*pEF, 27.9*pEF,
		31.9*pEF, 36.8*pEF, 43.6*pEF, 50.2*pEF,
		56.9*pEF, 65.7*pEF, 81.3*pEF, 101.6*pEF,
		116.5*pEF, 136.3*pEF, 166.15*pEF, 187.1*pEF,
		218.6*pEF, 260.54*pEF, 323.5*pEF, 387.5*pEF,
		451.54*pEF, 539.9*pEF };
	assert(sizeof(electronYield) == sizeof(particleEnergy));
	assert(sizeof(protonYield) == sizeof(particleEnergy));
	assert(sizeof(alphaYield) == sizeof(particleEnergy));
	assert(sizeof(ionYield) == sizeof(particleEnergy));
	scintillatorMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield, energyPoints)->SetSpline(true);
	scintillatorMPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield, energyPoints)->SetSpline(true);
	scintillatorMPT->AddProperty("DEUTERONSCINTILLATIONYIELD", particleEnergy, protonYield, energyPoints)->SetSpline(true);
	scintillatorMPT->AddProperty("TRITONSCINTILLATIONYIELD", particleEnergy, protonYield, energyPoints)->SetSpline(true);
	scintillatorMPT->AddProperty("ALPHASCINTILLATIONYIELD", particleEnergy, alphaYield, energyPoints)->SetSpline(true);
	scintillatorMPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield, energyPoints)->SetSpline(true);
	
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
	//ScintWrapper->SetType(dielectric_metal);  // dielectric and dielectric or metal?
	// teflon wrapping on polished surface->front/back painted // Teflon should be Lambertian in air, specular in optical grease
	//polished front painted is more simplified. Only specular spike reflection
	ScintWrapper->SetFinish(polishedfrontpainted);  
	//ScintWrapper->SetFinish(polished);  
	//ground front painted is more simplified. Only lambertain reflection
	//ScintWrapper->SetFinish(groundfrontpainted);  

/*		//poished back painted is maybe more realistic, need to then include sigma alpha (angle of the micro facet to the average normal surface) 
		ScintWrapper->SetFinish(polishedbackpainted);	
		//ScintWrapper->SetFinish(groundbackpainted);	
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
/*	ScintWrapper->SetModel(DAVIS);  // unified or glisur
	ScintWrapper->SetType(dielectric_LUTDAVIS);  // dielectric and dielectric or metal?
	ScintWrapper->SetFinish(PolishedESR_LUT);  
*/	//ScintWrapper->SetFinish(polishedtyvekair);  
	
	
	G4double rIndex_Teflon[numShort] = {1.35, 1.35, 1.35}; //Taken from wikipedia
	ScintWrapperMPT->AddProperty("RINDEX", photonEnergyShort, rIndex_Teflon, nEntriesShort)->SetSpline(true);  //refractive index can change with energy
	//G4double reflectivity[numShort] = {0.95, 0.95, 0.95};
	G4double reflectivity[numShort] = {0.99, 0.99, 0.99};
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
	//G4double pmtSize = 4.5*cm;
	//G4double pmtThick = 0.6*cm;//for 1 by 8 array pmtThick = 6 mm
	G4double pmtThick = 1.2*cm;
	pmtThick = fScintillatorWidth*pmtThick/1.5/10.; //Normalize to 1.5cm thickness
	G4double detWidth14 = 56.1429; //mm
	G4double calcRadii = (outerRadius - innerRadius - pmtThick)/2.; 
	//G4double pmt_width = 0.5*cm;
	G4double pmt_width = 1.2*cm;
	pmt_width = pmt_width*detWidth/detWidth14;
	G4double pmtSize = pmtThick+1.*cm;
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
	G4VSolid * SpherePMT_Face1 = new G4Sphere("SpherePMT_Face1", innerRadius-pmtThick, innerRadius, startPhi, deltaPhi, startTheta, deltaTheta);
	G4VSolid * SpherePMT_FaceLarger = new G4Sphere("SpherePMT_FaceLarger", innerRadius-pmtThick, innerRadius+fWrapThickness, startPhi, deltaPhi, startTheta, deltaTheta);
	G4VSolid * SpherePMT1 = new G4Sphere("SpherePMT1", innerRadius+calcRadii, outerRadius-calcRadii, startPhiPMT1, deltaPhiPMT1, startThetaPMT, deltaThetaPMT); //changed from innerRadius->Wrap
	G4VSolid * SpherePMT2 = new G4Sphere("SpherePMT2", innerRadius+calcRadii, outerRadius-calcRadii, startPhiPMT2, deltaPhiPMT2, startThetaPMT, deltaThetaPMT);//changed from innerRadius->Wrap
	//G4VSolid * SpherePMT1 = new G4Sphere("SpherePMT1", innerRadius+0.1*mm, outerRadius-0.1*mm, startPhiPMT1, deltaPhiPMT1, startThetaPMT, deltaThetaPMT); //changed from innerRadius->Wrap
	//G4VSolid * SpherePMT2 = new G4Sphere("SpherePMT2", innerRadius+0.1*mm, outerRadius-0.1*mm, startPhiPMT2, deltaPhiPMT2, startThetaPMT, deltaThetaPMT);//changed from innerRadius->Wrap
	//G4IntersectionSolid * interSolidBarsPMT1;
	//G4IntersectionSolid * interSolidBarsPMT2;
	//G4IntersectionSolid * PMT1;
	//G4IntersectionSolid * PMT2;
	G4IntersectionSolid * PMT11;
	G4IntersectionSolid * PMT12;
	G4IntersectionSolid * PMT13;
	G4IntersectionSolid * PMT21;
	G4IntersectionSolid * PMT22;
	G4IntersectionSolid * PMT23;
	//G4IntersectionSolid * PMT2BeamLine;
	G4IntersectionSolid * PMT21BeamLine;
	G4IntersectionSolid * PMT22BeamLine;
	G4IntersectionSolid * PMT23BeamLine;
	//G4IntersectionSolid * PMT1BeamLine;
	G4IntersectionSolid * PMT11BeamLine;
	G4IntersectionSolid * PMT12BeamLine;
	G4IntersectionSolid * PMT13BeamLine;
	//G4SubtractionSolid * PMT1; //old
	//G4SubtractionSolid * PMT2; //old
	//G4UnionSolid * Bar_PMT1;
	G4UnionSolid * Bar_PMT11;
	G4UnionSolid * Bar_PMT12;
	G4UnionSolid * Bar_PMT13;
	G4UnionSolid * Bar_PMT21;
	G4UnionSolid * Bar_PMT22;
	G4UnionSolid * outsidePMT_bar;	
	G4UnionSolid * outsidePMT_bar_Face1;	
	G4UnionSolid * outsidePMT_bar_Face2;	
	G4UnionSolid * outsidePMT_bar_Face3;	
	G4UnionSolid * outsidePMT_bar_AllFace;	
	G4UnionSolid * outsidePMT_bar_AllFace1;	
	G4IntersectionSolid * PMT_Front11;
	G4IntersectionSolid * PMT_Front21;
	G4IntersectionSolid * PMT_Front12;
	G4IntersectionSolid * PMT_Front22;
	G4IntersectionSolid * PMT_FrontMid1;
	G4IntersectionSolid * PMT_FrontMid2;
	G4IntersectionSolid * PMT_FrontLarger11;
	G4IntersectionSolid * PMT_FrontLarger12;
	G4IntersectionSolid * PMT_FrontLarger21;
	G4IntersectionSolid * PMT_FrontLarger22;
	G4IntersectionSolid * PMT_FrontLargerMid1;
	G4IntersectionSolid * PMT_FrontLargerMid2;
	//G4IntersectionSolid * PMT1Larger;
	//G4IntersectionSolid * PMT2Larger;
	G4IntersectionSolid * PMT11Larger;
	G4IntersectionSolid * PMT12Larger;
	G4IntersectionSolid * PMT13Larger;
	G4IntersectionSolid * PMT21Larger;
	G4IntersectionSolid * PMT22Larger;
	G4IntersectionSolid * PMT23Larger;


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
		G4String name0 = "BarsPlasticDet_";
		G4String nameLog=name0+std::to_string(i);
		G4String name1 = "wrapper_";
		G4String nameWrapper = name1+std::to_string(i);
		G4String name21 = "PMT_top1_";
		G4String namePMT11 = name21+std::to_string(i);
		G4String name22 = "PMT_top2_";
		G4String namePMT12 = name22+std::to_string(i);
		G4String name23 = "PMT_top3_";
		G4String namePMT13 = name23+std::to_string(i);
		G4String name31 = "PMT_bottom1_";
		G4String namePMT21 = name31+std::to_string(i);
		G4String name32 = "PMT_bottom2_";
		G4String namePMT22 = name32+std::to_string(i);
		G4String name33 = "PMT_bottom3_";
		G4String namePMT23 = name33+std::to_string(i);
		G4String name41 = "PMTFront_top1_";
		G4String namePMTFront11 = name41+std::to_string(i);
		G4String name42 = "PMTFront_top2_";
		G4String namePMTFront12 = name42+std::to_string(i);
		G4String name51 = "PMTFront_bottom1_";
		G4String namePMTFront21 = name51+std::to_string(i);
		G4String name52 = "PMTFront_bottom2_";
		G4String namePMTFront22 = name52+std::to_string(i);
		G4String name6 = "PMTFront_mid1_";
		G4String namePMTFrontMid1 = name6+std::to_string(i);
		G4String name7 = "PMTFront_mid2_";
		G4String namePMTFrontMid2 = name7+std::to_string(i);

		//For new geometry placement of the outside PMTs ie the vertical bar cutoff
		G4double YFacePMT1 = sqrt(pow(innerRadius, 2.0) - pow(abs(startPos)+detWidth/2.,2.0) - pow(zStartPos, 2.));

		//For Plastics New
		G4VSolid * boxBars = new G4Box("boxBars", detWidth/2., YFacePMT1 , 1.*m);
		interSolidBars = new G4IntersectionSolid("interSolidBars", SphereBars, boxBars, rot1, G4ThreeVector(startPos, 0, 50*cm));

		//For PMT's Top & Bottom New
		//G4VSolid * boxBarsPMTSmaller = new G4Box("boxBarsPMTSmaller", detWidth/2.-0.1*mm, pmt_width/2. , 1.*m);
		//PMT1 = new G4IntersectionSolid("PMT1", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, YFacePMT1+pmt_width/2., 50*cm));
		//PMT2 = new G4IntersectionSolid("PMT2", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		G4VSolid * boxBarsPMTSmaller = new G4Box("boxBarsPMTSmaller", pmt_width/2., pmt_width/2. , 1.*m);
		PMT11 = new G4IntersectionSolid("PMT11", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos - pmt_width -1.*mm, YFacePMT1+pmt_width/2., 50*cm));
		PMT12 = new G4IntersectionSolid("PMT12", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, YFacePMT1+pmt_width/2., 50*cm));
		PMT13 = new G4IntersectionSolid("PMT13", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos + pmt_width+1.*mm, YFacePMT1+pmt_width/2., 50*cm));
		PMT21 = new G4IntersectionSolid("PMT21", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos - pmt_width -1.*mm, -YFacePMT1-pmt_width/2., 50*cm));
		PMT22 = new G4IntersectionSolid("PMT22", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		PMT23 = new G4IntersectionSolid("PMT23", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos + pmt_width +1.*mm, -YFacePMT1-pmt_width/2., 50*cm));

		//Calc beamline middle face pmt placement above/below beamline
		G4double y1=YFacePMT1, y2=(BeamLineXY+pmtSize) , z1, z2, y3, z3;
		z1 = sqrt(innerRadius*innerRadius - startPos*startPos - y1*y1); //is this really different from zStartPos?  -> startPos vs startPos +detWidth/2
		G4double r1 = sqrt(y1*y1 + z1*z1);
		z2 = sqrt(r1*r1 - y2*y2);
		y3 = (y1+y2)/2.;
		z3 = (z1+z2)/2.;
		G4double r3 = sqrt(y3*y3 + z3*z3);
		G4double BeamlineY = y3*r1/r3;
		//Calculate Placefront of front top and front bot pmts
		G4ThreeVector ZDir_V = G4ThreeVector( 0., 0., 1.);
		G4ThreeVector Top_V2D = G4ThreeVector(0., YFacePMT1, z1);
		G4double theta = ZDir_V.angle(Top_V2D)/3.;
		G4double YFrontPlacement = r1*sin(theta);
		G4cout << "YFrontPlacement: " << YFrontPlacement << G4endl;
		//G4cout << "YFacePMT1/3: " << YFacePMT1/3. << G4endl;
		G4cout << "BeamlineY: " << BeamlineY << G4endl;

		//Face PMTs, placing them along the vertical at -1/3 and +1/3
		//G4VSolid * boxBarsPMTFace = new G4Box("boxBarsPMTFace", detWidth/2.-0.1*mm, pmtThick/2. , 1.*m); //for the 1 by 8 configuration
		//G4VSolid * boxBarsPMTFace = new G4Box("boxBarsPMTFace", detWidth/2.-0.1*mm, pmtThick , 1.*m); // for the 2, 4 by 4 configuration
		G4VSolid * boxBarsPMTFace1 = new G4Box("boxBarsPMTFace1", pmt_width, pmt_width , 1.*m); // for the 2, 4 by 4 configuration, independant 1
		//PMT_Front1 = new G4IntersectionSolid("PMT_Front1", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, YFrontPlacement, 50*cm)); //Y used to be YFacePMT1/3.
		//PMT_Front2 = new G4IntersectionSolid("PMT_Front2", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -YFrontPlacement, 50.*cm)); //Y used to be -YFacePMT1/3.
		PMT_Front11 = new G4IntersectionSolid("PMT_Front11", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, YFrontPlacement, 50*cm)); //Y used to be YFacePMT1/3. independant 1
		PMT_Front12 = new G4IntersectionSolid("PMT_Front12", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, YFrontPlacement, 50*cm)); //Y used to be YFacePMT1/3. independant 2
		PMT_Front21 = new G4IntersectionSolid("PMT_Front21", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, -YFrontPlacement, 50.*cm)); //Y used to be -YFacePMT1/3. independant 1
		PMT_Front22 = new G4IntersectionSolid("PMT_Front22", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, -YFrontPlacement, 50.*cm)); //Y used to be -YFacePMT1/3. independant 2

		//Face PMTs, placing the mid pmts  the vertical at 0 for first and last, or the middle wrt to the beamline if above or below
		if (i == 0 || i == (fNumDet-1) ) {
			//PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, 0., 0.));
			PMT_FrontMid1 = new G4IntersectionSolid("PMT_FrontMid1", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, 0., 0.)); //independant 1
			PMT_FrontMid2 = new G4IntersectionSolid("PMT_FrontMid2", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, 0., 0.)); //independant 2
		} else{
			//PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 0.)); //Y centered
			//PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, BeamlineY, 0.)); //Center of arc
			PMT_FrontMid1 = new G4IntersectionSolid("PMT_FrontMid1", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, BeamlineY, 0.)); //Center of arc, independant 1
			PMT_FrontMid2 = new G4IntersectionSolid("PMT_FrontMid2", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, BeamlineY, 0.)); //Center of arc, independant 2
		}


		//For wrapping Made slightly larger versions of all pmts and bar to be subtracted from the optical wrapping
		//Slightly larger Bar
		G4VSolid * boxBarsLarger = new G4Box("boxBarsLarger", detWidth/2.+fAirGap, YFacePMT1+2.*fAirGap , 1.*m);
		interSolidBarsLarger = new G4IntersectionSolid("interSolidBarsLarger", SphereInnerWrap, boxBarsLarger, rot1, G4ThreeVector(startPos, 0, 50*cm));
		//The actual optical wrap that will be subtracted from
		G4VSolid * boxWrap = new G4Box("boxWrap", wrapBoxThick, YFacePMT1+2.*fAirGap+fWrapThickness , 1.*m);
		interSolidWrap = new G4IntersectionSolid("interSolidWrap", SphereWrap, boxWrap, rot1, G4ThreeVector(startPos, 0, 50*cm));
		//Slightly larger PMTs
		//G4VSolid * boxBarsPMTSmaller2 = new G4Box("boxBarsPMTSmaller2", detWidth/2.-0.1*mm, pmt_width/2.+fWrapThickness , 1.*m);
		//PMT1Larger = new G4IntersectionSolid("PMT1Larger", SpherePMT1, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos, YFacePMT1+pmt_width/2., 50*cm));
		//PMT2Larger = new G4IntersectionSolid("PMT2Larger", SpherePMT2, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		G4VSolid * boxBarsPMTSmaller2 = new G4Box("boxBarsPMTSmaller2", pmt_width/2., pmt_width/2.+fWrapThickness , 1.*m);
		PMT11Larger = new G4IntersectionSolid("PMT11Larger", SpherePMT1, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos- pmt_width -1.*mm, YFacePMT1+pmt_width/2., 50*cm));
		PMT12Larger = new G4IntersectionSolid("PMT12Larger", SpherePMT1, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos, YFacePMT1+pmt_width/2., 50*cm));
		PMT13Larger = new G4IntersectionSolid("PMT13Larger", SpherePMT1, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos+ pmt_width +1.*mm, YFacePMT1+pmt_width/2., 50*cm));
		PMT21Larger = new G4IntersectionSolid("PMT21Larger", SpherePMT2, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos- pmt_width -1.*mm, -YFacePMT1-pmt_width/2., 50*cm));
		PMT22Larger = new G4IntersectionSolid("PMT22Larger", SpherePMT2, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		PMT23Larger = new G4IntersectionSolid("PMT23Larger", SpherePMT2, boxBarsPMTSmaller2, rot1, G4ThreeVector(startPos+ pmt_width +1.*mm, -YFacePMT1-pmt_width/2., 50*cm));

		//PMT_FrontLarger1 = new G4IntersectionSolid("PMT_FrontLarger1", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, YFrontPlacement, 50*cm)); //Y used to be YFacePMT1/3.
		//PMT_FrontLarger2 = new G4IntersectionSolid("PMT_FrontLarger2", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -YFrontPlacement, 50.*cm)); //Y used to be -YFacePMT1/3.
		PMT_FrontLarger11 = new G4IntersectionSolid("PMT_FrontLarger11", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, YFrontPlacement, 50*cm)); //Y used to be YFacePMT1/3., independant 1
		PMT_FrontLarger12 = new G4IntersectionSolid("PMT_FrontLarger12", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, YFrontPlacement, 50*cm)); //Y used to be YFacePMT1/3., independant 2
		PMT_FrontLarger21 = new G4IntersectionSolid("PMT_FrontLarger21", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, -YFrontPlacement, 50.*cm)); //Y used to be -YFacePMT1/3., independant 1
		PMT_FrontLarger22 = new G4IntersectionSolid("PMT_FrontLarger22", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, -YFrontPlacement, 50.*cm)); //Y used to be -YFacePMT1/3., independant 2
		if (i == 0 || i == (fNumDet-1) ) {
			//PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, 0., 0.));
			PMT_FrontLargerMid1 = new G4IntersectionSolid("PMT_FrontLargerMid1", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, 0., 0.)); //independant 1
			PMT_FrontLargerMid2 = new G4IntersectionSolid("PMT_FrontLargerMid2", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, 0., 0.)); //independant 2
		} else{
			//PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 0.)); //YCentered
			//PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, BeamlineY, 0.)); //Center of Arc
			PMT_FrontLargerMid1 = new G4IntersectionSolid("PMT_FrontLargerMid1", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, BeamlineY, 0.)); //Center of Arc, independant 1
			PMT_FrontLargerMid2 = new G4IntersectionSolid("PMT_FrontLargerMid2", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, BeamlineY, 0.)); //Center of Arc, independant 2
		}


		//Adding all Larger volumes together to be subtracted later
		//Bar_PMT1 = new G4UnionSolid("Bar_PMT1", PMT1Larger, interSolidBarsLarger, rot1, G4ThreeVector(0.,0.,0.));
		//outsidePMT_bar = new G4UnionSolid("outsidePMT_bar", Bar_PMT1, PMT2Larger, rot1, G4ThreeVector(0.,0.,0.));
		Bar_PMT11 = new G4UnionSolid("Bar_PMT11", PMT11Larger, interSolidBarsLarger, rot1, G4ThreeVector(0.,0.,0.));
		Bar_PMT12 = new G4UnionSolid("Bar_PMT12", Bar_PMT11, PMT12Larger, rot1, G4ThreeVector(0.,0.,0.));
		Bar_PMT13 = new G4UnionSolid("Bar_PMT13", Bar_PMT12, PMT13Larger, rot1, G4ThreeVector(0.,0.,0.));
		Bar_PMT21 = new G4UnionSolid("Bar_PMT21", Bar_PMT13, PMT21Larger, rot1, G4ThreeVector(0.,0.,0.));
		Bar_PMT22 = new G4UnionSolid("Bar_PMT22", Bar_PMT21, PMT22Larger, rot1, G4ThreeVector(0.,0.,0.));
		outsidePMT_bar = new G4UnionSolid("outsidePMT_bar", Bar_PMT22, PMT23Larger, rot1, G4ThreeVector(0.,0.,0.));

		if( abs(startPos-detWidth/2.) < BeamLineXY || abs(startPos+detWidth/2.) < BeamLineXY || abs(startPos) < BeamLineXY || i == 0 || i == (fNumDet-1) ) {
			//outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar, PMT_FrontLargerMid, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_AllFace1 = new G4UnionSolid("outsidePMT_bar_AllFace1", outsidePMT_bar, PMT_FrontLargerMid1, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar_AllFace1, PMT_FrontLargerMid2, rot1, G4ThreeVector(0.,0.,0.));
		} else{
			//outsidePMT_bar_Face1 = new G4UnionSolid("outsidePMT_bar_Face1", outsidePMT_bar, PMT_FrontLarger1, rot1, G4ThreeVector(0.,0.,0.));
			//outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar_Face1, PMT_FrontLarger2, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_Face1 = new G4UnionSolid("outsidePMT_bar_Face1", outsidePMT_bar, PMT_FrontLarger11, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_Face2 = new G4UnionSolid("outsidePMT_bar_Face2", outsidePMT_bar_Face1, PMT_FrontLarger12, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_Face3 = new G4UnionSolid("outsidePMT_bar_Face3", outsidePMT_bar_Face2, PMT_FrontLarger21, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar_Face3, PMT_FrontLarger22, rot1, G4ThreeVector(0.,0.,0.));
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
			//PMT2BeamLine = new G4IntersectionSolid("PMT2BeamLine", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, pmtSize+BeamLineXY-pmt_width/2., 50*cm));
			PMT21BeamLine = new G4IntersectionSolid("PMT21BeamLine1", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos- pmt_width -1.*mm, pmtSize+BeamLineXY-pmt_width/2., 50*cm));
			PMT22BeamLine = new G4IntersectionSolid("PMT22BeamLine2", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, pmtSize+BeamLineXY-pmt_width/2., 50*cm));
			PMT23BeamLine = new G4IntersectionSolid("PMT23BeamLine3", SpherePMT1, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos+ pmt_width +1.*mm, pmtSize+BeamLineXY-pmt_width/2., 50*cm));

			//For wrapping Made slightly larger versions of all pmts and bar to be subtracted from the optical wrapping
			//Slightly larger Bar
			G4VSolid * boxBarsLarger2 = new G4Box("boxBarsLarger2", detWidth/2.+fAirGap, (YFacePMT1-pmtSize-BeamLineXY)/2.+2.*fAirGap , 1.*m);
			interSolidBarsLarger = new G4IntersectionSolid("interSolidBarsLarger", SphereInnerWrap, boxBarsLarger2, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 50*cm));

			//The actual optical wrap that will be subtracted from
			G4VSolid * boxWrap2 = new G4Box("boxWrap2", wrapBoxThick, (YFacePMT1-pmtSize-BeamLineXY)/2.+2.*fAirGap+fWrapThickness , 1.*m);
			interSolidWrap = new G4IntersectionSolid("interSolidWrap", SphereWrap, boxWrap2, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 50*cm));

			//Slightly larger PMTs
			//G4VSolid * boxBarsPMTSmaller3 = new G4Box("boxBarsPMTSmaller3", detWidth/2.-0.1*mm, pmt_width/2.+fWrapThickness , 1.*m);
			//PMT1Larger = new G4IntersectionSolid("PMT1Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, YFacePMT1+pmt_width/2., 50*cm));
			//PMT2Larger = new G4IntersectionSolid("PMT2Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, pmtSize+BeamLineXY-pmt_width/2., 50*cm));
			G4VSolid * boxBarsPMTSmaller3 = new G4Box("boxBarsPMTSmaller3", pmt_width/2., pmt_width/2.+fWrapThickness , 1.*m);
			PMT11Larger = new G4IntersectionSolid("PMT11Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos- pmt_width -1.*mm, YFacePMT1+pmt_width/2., 50*cm));
			PMT12Larger = new G4IntersectionSolid("PMT12Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, YFacePMT1+pmt_width/2., 50*cm));
			PMT13Larger = new G4IntersectionSolid("PMT13Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos+ pmt_width +1.*mm, YFacePMT1+pmt_width/2., 50*cm));
			PMT21Larger = new G4IntersectionSolid("PMT21Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos- pmt_width -1.*mm, pmtSize+BeamLineXY-pmt_width/2., 50*cm));
			PMT22Larger = new G4IntersectionSolid("PMT22Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, pmtSize+BeamLineXY-pmt_width/2., 50*cm));
			PMT23Larger = new G4IntersectionSolid("PMT23Larger", SpherePMT1, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos+ pmt_width +1.*mm, pmtSize+BeamLineXY-pmt_width/2., 50*cm));

			//Adding all Larger volumes together to be subtracted later
			//Bar_PMT1 = new G4UnionSolid("Bar_PMT1", PMT1Larger, interSolidBarsLarger, rot1, G4ThreeVector(0.,0.,0.));
			//outsidePMT_bar = new G4UnionSolid("outsidePMT_bar", Bar_PMT1, PMT2Larger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT11 = new G4UnionSolid("Bar_PMT11", PMT11Larger, interSolidBarsLarger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT12 = new G4UnionSolid("Bar_PMT12", Bar_PMT11, PMT12Larger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT13 = new G4UnionSolid("Bar_PMT13", Bar_PMT12, PMT13Larger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT21 = new G4UnionSolid("Bar_PMT21", Bar_PMT13, PMT21Larger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT22 = new G4UnionSolid("Bar_PMT22", Bar_PMT21, PMT22Larger, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar = new G4UnionSolid("outsidePMT_bar", Bar_PMT22, PMT23Larger, rot1, G4ThreeVector(0.,0.,0.));
			//outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar, PMT_FrontLargerMid, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_AllFace1 = new G4UnionSolid("outsidePMT_bar_AllFace1", outsidePMT_bar, PMT_FrontLargerMid1, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar_AllFace1, PMT_FrontLargerMid2, rot1, G4ThreeVector(0.,0.,0.));

			//Subtraction to create actual optical wrapping shape
			if(fAddFrontPMTs == true) {
				subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar_AllFace, rot1, G4ThreeVector(0,0,0));
			} else subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar, rot1, G4ThreeVector(0,0,0));

			//Assign Logical Volume for detectors and wrapping affected by beamline
			fPlasticLogArray[i] = new G4LogicalVolume(interSolidBars, plasticG4material, nameLog,0,0,0);
			fWrapLogArray[i] = new G4LogicalVolume(subtractSolidWrap, wrapG4material, nameWrapper,0,0,0); 
			fPMT11LogArray[i] = new G4LogicalVolume(PMT11, PMTG4material, namePMT11,0,0,0);
			fPMT12LogArray[i] = new G4LogicalVolume(PMT12, PMTG4material, namePMT12,0,0,0);
			fPMT13LogArray[i] = new G4LogicalVolume(PMT13, PMTG4material, namePMT13,0,0,0);
			fPMT21LogArray[i] = new G4LogicalVolume(PMT21BeamLine, PMTG4material, namePMT21,0,0,0);
			fPMT22LogArray[i] = new G4LogicalVolume(PMT22BeamLine, PMTG4material, namePMT22,0,0,0);
			fPMT23LogArray[i] = new G4LogicalVolume(PMT23BeamLine, PMTG4material, namePMT23,0,0,0);
			fPMTFace11LogArray[i] = new G4LogicalVolume(PMT_Front11, PMTG4material, namePMTFront11,0,0,0);
			fPMTFace12LogArray[i] = new G4LogicalVolume(PMT_Front12, PMTG4material, namePMTFront12,0,0,0);
			fPMTFace21LogArray[i] = new G4LogicalVolume(PMT_Front21, PMTG4material, namePMTFront21,0,0,0);
			fPMTFace22LogArray[i] = new G4LogicalVolume(PMT_Front22, PMTG4material, namePMTFront22,0,0,0);
			fPMTFaceMid1LogArray[i] = new G4LogicalVolume(PMT_FrontMid1, PMTG4material, namePMTFrontMid1,0,0,0);
			fPMTFaceMid2LogArray[i] = new G4LogicalVolume(PMT_FrontMid2, PMTG4material, namePMTFrontMid2,0,0,0);
		}

		else 
		{
			//Assign Logical Volume for detectors and wrapping unaffected by beamline
			fPlasticLogArray[i] = new G4LogicalVolume(interSolidBars, plasticG4material, nameLog,0,0,0);
			fWrapLogArray[i] = new G4LogicalVolume(subtractSolidWrap, wrapG4material, nameWrapper,0,0,0);
			fPMT11LogArray[i] = new G4LogicalVolume(PMT11, PMTG4material, namePMT11,0,0,0);
			fPMT12LogArray[i] = new G4LogicalVolume(PMT12, PMTG4material, namePMT12,0,0,0);
			fPMT13LogArray[i] = new G4LogicalVolume(PMT13, PMTG4material, namePMT13,0,0,0);
			fPMT21LogArray[i] = new G4LogicalVolume(PMT21, PMTG4material, namePMT21,0,0,0);
			fPMT22LogArray[i] = new G4LogicalVolume(PMT22, PMTG4material, namePMT22,0,0,0);
			fPMT23LogArray[i] = new G4LogicalVolume(PMT23, PMTG4material, namePMT23,0,0,0);
			fPMTFace11LogArray[i] = new G4LogicalVolume(PMT_Front11, PMTG4material, namePMTFront11,0,0,0);
			fPMTFace12LogArray[i] = new G4LogicalVolume(PMT_Front12, PMTG4material, namePMTFront12,0,0,0);
			fPMTFace21LogArray[i] = new G4LogicalVolume(PMT_Front21, PMTG4material, namePMTFront21,0,0,0);
			fPMTFace22LogArray[i] = new G4LogicalVolume(PMT_Front22, PMTG4material, namePMTFront22,0,0,0);
			fPMTFaceMid1LogArray[i] = new G4LogicalVolume(PMT_FrontMid1, PMTG4material, namePMTFrontMid1,0,0,0);
			fPMTFaceMid2LogArray[i] = new G4LogicalVolume(PMT_FrontMid2, PMTG4material, namePMTFrontMid2,0,0,0);
		}

		//Set Logical Skin for optical photons on wrapping
		G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapper, fWrapLogArray[i], ScintWrapper);

		//Give everything colour
		fPlasticLogArray[i]->SetVisAttributes(plastic_vis);
		fWrapLogArray[i]->SetVisAttributes(wrap_vis);
		fPMT11LogArray[i]->SetVisAttributes(pmt_vis);
		fPMT12LogArray[i]->SetVisAttributes(pmt_vis);
		fPMT13LogArray[i]->SetVisAttributes(pmt_vis);
		fPMT21LogArray[i]->SetVisAttributes(pmt_vis);
		fPMT22LogArray[i]->SetVisAttributes(pmt_vis);
		fPMT23LogArray[i]->SetVisAttributes(pmt_vis);
		fPMTFaceMid1LogArray[i]->SetVisAttributes(pmt_vis);
		fPMTFaceMid2LogArray[i]->SetVisAttributes(pmt_vis);
		fPMTFace11LogArray[i]->SetVisAttributes(pmt_vis);
		fPMTFace12LogArray[i]->SetVisAttributes(pmt_vis);
		fPMTFace21LogArray[i]->SetVisAttributes(pmt_vis);
		fPMTFace22LogArray[i]->SetVisAttributes(pmt_vis);

		//build every second detector
		//if(i%2==0) {}
		//Place detectors
		fAssemblyPlastics->AddPlacedVolume(fPlasticLogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT11LogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT12LogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT13LogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT21LogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT22LogArray[i], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT23LogArray[i], move, rotate);
		if(fAddFrontPMTs == true) {
			if( abs(startPos-detWidth/2.) < BeamLineXY || abs(startPos+detWidth/2.) < BeamLineXY || abs(startPos) < BeamLineXY || i == 0 || i == (fNumDet-1) ) {
				fAssemblyPlastics->AddPlacedVolume(fPMTFaceMid1LogArray[i], move, rotate);
				fAssemblyPlastics->AddPlacedVolume(fPMTFaceMid2LogArray[i], move, rotate);
			} else {
				fAssemblyPlastics->AddPlacedVolume(fPMTFace11LogArray[i], move, rotate);
				fAssemblyPlastics->AddPlacedVolume(fPMTFace12LogArray[i], move, rotate);
				fAssemblyPlastics->AddPlacedVolume(fPMTFace21LogArray[i], move, rotate);
				fAssemblyPlastics->AddPlacedVolume(fPMTFace22LogArray[i], move, rotate);
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
	fPMT11LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMT12LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMT13LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMT21LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMT22LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMT23LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMTFace11LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMTFace12LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMTFace21LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMTFace22LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMTFaceMid1LogArray.resize(fNumDetBot+fNumDet, NULL);
	fPMTFaceMid2LogArray.resize(fNumDetBot+fNumDet, NULL);


	//Loop for building bottom detectors
	for(G4int k=0; k<fNumDetBot; ++k) {
		G4cout << "startPos In Bottom loop: " << startPos << G4endl;
		//Names
		G4int detNumBottom = fNumDet + k;
		G4String name0 = "BarsPlasticDet_";
		G4String nameLog=name0+std::to_string(detNumBottom);
		G4String name1 = "wrapper_";
		G4String nameWrapper = name1+std::to_string(detNumBottom);
		G4String name21 = "PMT_top1_";
		G4String namePMT11 = name21+std::to_string(detNumBottom);
		G4String name22 = "PMT_top2_";
		G4String namePMT12 = name22+std::to_string(detNumBottom);
		G4String name23 = "PMT_top3_";
		G4String namePMT13 = name23+std::to_string(detNumBottom);
		G4String name31 = "PMT_bottom1_";
		G4String namePMT21 = name31+std::to_string(detNumBottom);
		G4String name32 = "PMT_bottom2_";
		G4String namePMT22 = name32+std::to_string(detNumBottom);
		G4String name33 = "PMT_bottom3_";
		G4String namePMT23 = name33+std::to_string(detNumBottom);
		G4String name41 = "PMTFront_top1_";
		G4String namePMTFront11 = name41+std::to_string(detNumBottom);
		G4String name42 = "PMTFront_top2_";
		G4String namePMTFront12 = name42+std::to_string(detNumBottom);
		G4String name51 = "PMTFront_bottom1_";
		G4String namePMTFront21 = name51+std::to_string(detNumBottom);
		G4String name52 = "PMTFront_bottom2_";
		G4String namePMTFront22 = name52+std::to_string(detNumBottom);
		G4String name6 = "PMTFront_mid1_";
		G4String namePMTFrontMid1 = name6+std::to_string(detNumBottom);
		G4String name7 = "PMTFront_mid2_";
		G4String namePMTFrontMid2 = name7+std::to_string(detNumBottom);

		//For new geometry placement of the outside PMTs ie the vertical bar cutoff
		G4double YFacePMT1 = sqrt(pow(innerRadius, 2.0) - pow(abs(startPos)+detWidth/2.,2.0) - pow(zStartPos, 2.));

		//For Plastics New
		G4VSolid * boxBars2 = new G4Box("boxBars2", detWidth/2., (YFacePMT1-pmtSize-BeamLineXY)/2. , 1.*m);
		interSolidBars = new G4IntersectionSolid("interSolidBars", SphereBars, boxBars2, rot1, G4ThreeVector(startPos, (-YFacePMT1-pmtSize-BeamLineXY)/2., 50*cm));

		//For PMT's Top & Bottom New
		//G4VSolid * boxBarsPMTSmaller = new G4Box("boxBarsPMTSmaller", detWidth/2.-0.1*mm, pmt_width/2. , 1.*m);
		//PMT2 = new G4IntersectionSolid("PMT2", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		G4VSolid * boxBarsPMTSmaller = new G4Box("boxBarsPMTSmaller", pmt_width/2., pmt_width/2. , 1.*m);
		PMT21 = new G4IntersectionSolid("PMT21", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos - pmt_width -1.*mm, -YFacePMT1-pmt_width/2., 50*cm));
		PMT22 = new G4IntersectionSolid("PMT22", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		PMT23 = new G4IntersectionSolid("PMT23", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos + pmt_width +1.*mm, -YFacePMT1-pmt_width/2., 50*cm));
		//Only Need to only make top pmt here, note SpherePMT2 because below y=0
		//PMT1BeamLine = new G4IntersectionSolid("PMT1BeamLine", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));
		PMT11BeamLine = new G4IntersectionSolid("PMT11BeamLine1", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos- pmt_width -1.*mm, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));
		PMT12BeamLine = new G4IntersectionSolid("PMT12BeamLine2", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));
		PMT13BeamLine = new G4IntersectionSolid("PMT13BeamLine3", SpherePMT2, boxBarsPMTSmaller, rot1, G4ThreeVector(startPos+ pmt_width +1.*mm, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));

		//Calc beamline middle face pmt placement above/below beamline
		G4double y1=YFacePMT1, y2=(BeamLineXY+pmtSize) , z1, z2, y3, z3;
		z1 = sqrt(innerRadius*innerRadius - startPos*startPos - y1*y1); //is this really different from zStartPos?  -> startPos vs startPos +detWidth/2
		G4double r1 = sqrt(y1*y1 + z1*z1);
		z2 = sqrt(r1*r1 - y2*y2);
		y3 = (y1+y2)/2.;
		z3 = (z1+z2)/2.;
		G4double r3 = sqrt(y3*y3 + z3*z3);
		G4double BeamlineY = y3*r1/r3;
		//Calculate Placefront of front top and front bot pmts
		G4ThreeVector ZDir_V = G4ThreeVector( 0., 0., 1.);
		G4ThreeVector Top_V2D = G4ThreeVector(0., YFacePMT1, z1);
		G4double theta = ZDir_V.angle(Top_V2D)/3.;
		G4double YFrontPlacement = r1*sin(theta);
		G4cout << "YFrontPlacement: " << YFrontPlacement << G4endl;
		G4cout << "YFacePMT1/3: " << YFacePMT1/3. << G4endl;

		//Face PMTs, placing them along the vertical at -1/3 and +1/3
		//G4VSolid * boxBarsPMTFace = new G4Box("boxBarsPMTFace", detWidth/2.-0.1*mm, pmtThick/2. , 1.*m); // for the 1 by 8 configuration
		//G4VSolid * boxBarsPMTFace = new G4Box("boxBarsPMTFace", detWidth/2.-0.1*mm, pmtThick , 1.*m); // for the 2, 4 by 4 configuration
		G4VSolid * boxBarsPMTFace1 = new G4Box("boxBarsPMTFace1", pmt_width, pmt_width , 1.*m); // for the 2, 4 by 4 configuration, independant 1
		//PMT_Front1 = new G4IntersectionSolid("PMT_Front1", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, YFrontPlacement, 50*cm));
		//PMT_Front2 = new G4IntersectionSolid("PMT_Front2", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -YFrontPlacement, 50.*cm));
		PMT_Front11 = new G4IntersectionSolid("PMT_Front11", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, YFrontPlacement, 50*cm)); //Y used to be YFacePMT1/3. independant 1
		PMT_Front12 = new G4IntersectionSolid("PMT_Front12", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, YFrontPlacement, 50*cm)); //Y used to be YFacePMT1/3. independant 2
		PMT_Front21 = new G4IntersectionSolid("PMT_Front21", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, -YFrontPlacement, 50.*cm)); //Y used to be -YFacePMT1/3. independant 1
		PMT_Front22 = new G4IntersectionSolid("PMT_Front22", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, -YFrontPlacement, 50.*cm)); //Y used to be -YFacePMT1/3. independant 2

		//Face PMTs, placing the mid pmts  the vertical at 0 for first and last, or the middle wrt to the beamline if above or below
		//PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 0.)); //Y centered
		//PMT_FrontMid = new G4IntersectionSolid("PMT_FrontMid", SpherePMT_Face1, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -BeamlineY, 0.)); //Center of arc
		PMT_FrontMid1 = new G4IntersectionSolid("PMT_FrontMid1", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, -BeamlineY, 0.)); //Center of arc, independant 1
		PMT_FrontMid2 = new G4IntersectionSolid("PMT_FrontMid2", SpherePMT_Face1, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, -BeamlineY, 0.)); //Center of arc, independant 2

		//For wrapping Made slightly larger versions of all pmts and bar to be subtracted from the optical wrapping
		//Slightly larger Bar
		G4VSolid * boxBarsLarger2 = new G4Box("boxBarsLarger2", detWidth/2.+fAirGap, (YFacePMT1-pmtSize-BeamLineXY)/2.+2.*fAirGap , 1.*m);
		interSolidBarsLarger = new G4IntersectionSolid("interSolidBarsLarger", SphereInnerWrap, boxBarsLarger2, rot1, G4ThreeVector(startPos, (-YFacePMT1-pmtSize-BeamLineXY)/2., 50*cm));

		//The actual optical wrap that will be subtracted from
		G4VSolid * boxWrap2 = new G4Box("boxWrap2", wrapBoxThick, (YFacePMT1-pmtSize-BeamLineXY)/2.+2.*fAirGap+fWrapThickness , 1.*m);
		interSolidWrap = new G4IntersectionSolid("interSolidWrap", SphereWrap, boxWrap2, rot1, G4ThreeVector(startPos, (-YFacePMT1-pmtSize-BeamLineXY)/2., 50*cm));

		//Slightly larger PMTs
		//G4VSolid * boxBarsPMTSmaller3 = new G4Box("boxBarsPMTSmaller3", detWidth/2.-0.1*mm, pmt_width/2.+fWrapThickness , 1.*m);
		//PMT2Larger = new G4IntersectionSolid("PMT2Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		//PMT1Larger = new G4IntersectionSolid("PMT1Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));
		G4VSolid * boxBarsPMTSmaller3 = new G4Box("boxBarsPMTSmaller3", pmt_width/2., pmt_width/2.+fWrapThickness , 1.*m);
		PMT21Larger = new G4IntersectionSolid("PMT21Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos- pmt_width -1.*mm, -YFacePMT1-pmt_width/2., 50*cm));
		PMT22Larger = new G4IntersectionSolid("PMT22Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, -YFacePMT1-pmt_width/2., 50*cm));
		PMT23Larger = new G4IntersectionSolid("PMT23Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos+ pmt_width +1.*mm, -YFacePMT1-pmt_width/2., 50*cm));
		PMT11Larger = new G4IntersectionSolid("PMT11Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos- pmt_width -1.*mm, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));
		PMT12Larger = new G4IntersectionSolid("PMT12Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));
		PMT13Larger = new G4IntersectionSolid("PMT13Larger", SpherePMT2, boxBarsPMTSmaller3, rot1, G4ThreeVector(startPos+ pmt_width +1.*mm, -pmtSize-BeamLineXY+pmt_width/2., 50*cm));
		//PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, (YFacePMT1+pmtSize+BeamLineXY)/2., 0.)); //YCentered
		//PMT_FrontLargerMid = new G4IntersectionSolid("PMT_FrontLargerMid", SpherePMT_FaceLarger, boxBarsPMTFace, rot1, G4ThreeVector(startPos, -BeamlineY, 0.)); //Center of Arc
		PMT_FrontLargerMid1 = new G4IntersectionSolid("PMT_FrontLargerMid1", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos-pmt_width-1.*mm, -BeamlineY, 0.)); //Center of Arc, independant 1
		PMT_FrontLargerMid2 = new G4IntersectionSolid("PMT_FrontLargerMid2", SpherePMT_FaceLarger, boxBarsPMTFace1, rot1, G4ThreeVector(startPos+pmt_width+1.*mm, -BeamlineY, 0.)); //Center of Arc, independant 2

		//Adding all Larger volumes together to be subtracted later
		//Bar_PMT1 = new G4UnionSolid("Bar_PMT1", PMT1Larger, interSolidBarsLarger, rot1, G4ThreeVector(0.,0.,0.));
		//outsidePMT_bar = new G4UnionSolid("outsidePMT_bar", Bar_PMT1, PMT2Larger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT11 = new G4UnionSolid("Bar_PMT11", PMT11Larger, interSolidBarsLarger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT12 = new G4UnionSolid("Bar_PMT12", Bar_PMT11, PMT12Larger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT13 = new G4UnionSolid("Bar_PMT13", Bar_PMT12, PMT13Larger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT21 = new G4UnionSolid("Bar_PMT21", Bar_PMT13, PMT21Larger, rot1, G4ThreeVector(0.,0.,0.));
			Bar_PMT22 = new G4UnionSolid("Bar_PMT22", Bar_PMT21, PMT22Larger, rot1, G4ThreeVector(0.,0.,0.));
			outsidePMT_bar = new G4UnionSolid("outsidePMT_bar", Bar_PMT22, PMT23Larger, rot1, G4ThreeVector(0.,0.,0.));
		//outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar, PMT_FrontLargerMid, rot1, G4ThreeVector(0.,0.,0.));
		outsidePMT_bar_AllFace1 = new G4UnionSolid("outsidePMT_bar_AllFace1", outsidePMT_bar, PMT_FrontLargerMid1, rot1, G4ThreeVector(0.,0.,0.));
		outsidePMT_bar_AllFace = new G4UnionSolid("outsidePMT_bar_AllFace", outsidePMT_bar_AllFace1, PMT_FrontLargerMid2, rot1, G4ThreeVector(0.,0.,0.));

		//Subtraction to create actual optical wrapping shape
		if(fAddFrontPMTs == true) {
			subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar_AllFace, rot1, G4ThreeVector(0,0,0));
		} else subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, outsidePMT_bar, rot1, G4ThreeVector(0,0,0));

		//Assign Logical Volume for detectors and wrapping affected by beamline
		fPlasticLogArray[detNumBottom] = new G4LogicalVolume(interSolidBars, plasticG4material, nameLog,0,0,0);
		fWrapLogArray[detNumBottom] = new G4LogicalVolume(subtractSolidWrap, wrapG4material, nameWrapper,0,0,0); 
		fPMT11LogArray[detNumBottom] = new G4LogicalVolume(PMT11BeamLine, PMTG4material, namePMT11,0,0,0);
		fPMT12LogArray[detNumBottom] = new G4LogicalVolume(PMT12BeamLine, PMTG4material, namePMT12,0,0,0);
		fPMT13LogArray[detNumBottom] = new G4LogicalVolume(PMT13BeamLine, PMTG4material, namePMT13,0,0,0);
		fPMT21LogArray[detNumBottom] = new G4LogicalVolume(PMT21, PMTG4material, namePMT21,0,0,0);
		fPMT22LogArray[detNumBottom] = new G4LogicalVolume(PMT22, PMTG4material, namePMT22,0,0,0);
		fPMT23LogArray[detNumBottom] = new G4LogicalVolume(PMT23, PMTG4material, namePMT23,0,0,0);
		fPMTFace11LogArray[detNumBottom] = new G4LogicalVolume(PMT_Front11, PMTG4material, namePMTFront11,0,0,0);
		fPMTFace12LogArray[detNumBottom] = new G4LogicalVolume(PMT_Front12, PMTG4material, namePMTFront12,0,0,0);
		fPMTFace21LogArray[detNumBottom] = new G4LogicalVolume(PMT_Front21, PMTG4material, namePMTFront21,0,0,0);
		fPMTFace22LogArray[detNumBottom] = new G4LogicalVolume(PMT_Front22, PMTG4material, namePMTFront22,0,0,0);
		fPMTFaceMid1LogArray[detNumBottom] = new G4LogicalVolume(PMT_FrontMid1, PMTG4material, namePMTFrontMid1,0,0,0);
		fPMTFaceMid2LogArray[detNumBottom] = new G4LogicalVolume(PMT_FrontMid2, PMTG4material, namePMTFrontMid2,0,0,0);

		//Set Logical Skin for optical photons on wrapping
		G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapper, fWrapLogArray[detNumBottom], ScintWrapper);

		//Give everything colour
		fPlasticLogArray[detNumBottom]->SetVisAttributes(plastic_vis);
		fWrapLogArray[detNumBottom]->SetVisAttributes(wrap_vis);
		fPMT11LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMT12LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMT13LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMT21LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMT22LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMT23LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMTFaceMid1LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMTFaceMid2LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMTFace11LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMTFace12LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMTFace21LogArray[detNumBottom]->SetVisAttributes(pmt_vis);
		fPMTFace22LogArray[detNumBottom]->SetVisAttributes(pmt_vis);

		//build every second detector
		//if(i%2==0) {}
		//Place detectors
		fAssemblyPlastics->AddPlacedVolume(fPlasticLogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT11LogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT12LogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT13LogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT21LogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT22LogArray[detNumBottom], move, rotate);
		fAssemblyPlastics->AddPlacedVolume(fPMT23LogArray[detNumBottom], move, rotate);
		if(fAddFrontPMTs == true) {
			fAssemblyPlastics->AddPlacedVolume(fPMTFaceMid1LogArray[detNumBottom], move, rotate);
			fAssemblyPlastics->AddPlacedVolume(fPMTFaceMid2LogArray[detNumBottom], move, rotate);
		}
		if(fAddWrap ==true){
			fAssemblyPlastics->AddPlacedVolume(fWrapLogArray[detNumBottom], move, rotate);
		}


		startPos = startPos+detWidth+2.*fWrapThickness+ 2.*fAirGap;


	}




	return 1;
}

