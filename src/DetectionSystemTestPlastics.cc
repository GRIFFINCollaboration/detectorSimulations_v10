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

#include "DetectionSystemTestPlastics.hh"

#include "G4SystemOfUnits.hh"

#include <string>

DetectionSystemTestPlastics::DetectionSystemTestPlastics(G4double thickness, G4int material, G4double numDet)
{

	fPMTWidth         = 1.*mm;
	fScintillatorWidth         = 1.*mm;
	fDiameter         = 25.*mm;
	fRadialDistance = 5.*mm; //Says in GRIFFIN nim can move from a few mm to 5 cm.  At closest covers roughly 25% of 4pi
	fStartPhi = 0.;
	fDeltaPhi = 2*M_PI;

	fWrapThickness = 0.5 * mm; //Factor of 10 should be  applied later for visualization purposes.  0.05 for simulations, 0.5 for visualization
	fAirGap = 0.;
	fWrapMaterial = "Teflon";
	fPMTMaterial = "G4_SILICON_DIOXIDE";
	fZDSMaterial = "BC422";
	//fZDSMaterial = "Dense"; // for Solid Angle test

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
DetectionSystemTestPlastics::~DetectionSystemTestPlastics() {
	// LogicalVolumes
	delete fPlasticLog;
	delete fWrapLog;
	delete fPMTLog;
	delete fZDSLog;
	delete fZDSPMTLog;
	G4cout << "Calling Destructor" << G4endl;

}
////////
/////////
G4int DetectionSystemTestPlastics::Build() {
	fAssemblyTestPlastics = new G4AssemblyVolume(); 
	G4cout << "Calling Build function" << G4endl;
	BuildTestPlastics();

	return 1;
}
////////
///////
G4int DetectionSystemTestPlastics::PlaceDetector(G4LogicalVolume* expHallLog) {
	G4RotationMatrix * rotate = new G4RotationMatrix;
	G4ThreeVector move = G4ThreeVector(0., 0., 0.);

	//To not check overlaps
	//    fAssemblyPlastics->MakeImprint(expHallLog, move, rotate);
	//To check overlaps
	fAssemblyTestPlastics->MakeImprint(expHallLog, move, rotate, 0, true);
	G4cout << "Calling place detector" << G4endl;
	return 1;
}
///////////
/////////
G4int DetectionSystemTestPlastics::BuildTestPlastics() {

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
	G4Material* zdsG4material = G4Material::GetMaterial(fZDSMaterial);
	if( !zdsG4material ) {
		G4cout << " ----> Material " << fZDSMaterial << " not found, cannot build! " << G4endl;
		return 0;
	}
	else {
		G4cout << zdsG4material->GetName() << " is the name of the detector material" << G4endl;
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
	// teflon wrapping on polished surface->front/back painted // Teflon should be Lambertian in air, specular in optical grease
	//polished front painted is more simplified. Only specular spike reflection
	ScintWrapper->SetFinish(polishedfrontpainted);  

	//ground front painted is more simplified. Only lambertain reflection
	//ScintWrapper->SetFinish(groundfrontpainted);  


	/*		//poished back painted is maybe more realistic, need to then include sigma alpha (angle of the micro facet to the average normal surface) 
	//ScintWrapper->SetFinish(polishedbackpainted);	
	ScintWrapper->SetFinish(groundbackpainted);	
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


	G4RotationMatrix *rot1 = new G4RotationMatrix(0,0,0);

	
	///////// ZDS and 1x1x1 with 1 SiPM ///////////////////

	double length = 5.*cm;
	G4VSolid * Scint = new G4Box("Scint", 1.*cm/2., length/2. , 1.*cm/2.);
	//G4VSolid * Scint = new G4Box("Scint", 1.*cm/2., 1.*cm/2. , 1.*cm/2.);
	G4VSolid * Wrap_Bigger = new G4Box("Wrap_Bigger", 1.*cm/2.+fWrapThickness, length/2.+fWrapThickness, 1.*cm/2.+fWrapThickness);
	G4SubtractionSolid * subtractWrap = new G4SubtractionSolid("subtractWrap", Wrap_Bigger, Scint, rot1, G4ThreeVector(0,0,0));
	//G4VSolid * Wrap_Bigger = new G4Box("Wrap_Bigger", 1.*cm/2.+fWrapThickness, 1.*cm/2.+fWrapThickness, 1.*cm/2.+fWrapThickness);
	//G4SubtractionSolid * subtractWrap = new G4SubtractionSolid("subtractWrap", Wrap_Bigger, Scint, rot1, G4ThreeVector(0,0,0));
	G4VSolid * PMT = new G4Box("PMT", 0.4*cm/2., 0.4*cm/2. , 0.4*cm/2.);
	G4SubtractionSolid * subtractWrap2 = new G4SubtractionSolid("subtractWrap2", subtractWrap, PMT, rot1, G4ThreeVector(0,-0.5*length-fWrapThickness,0));
	//G4SubtractionSolid * subtractWrap2 = new G4SubtractionSolid("subtractWrap2", subtractWrap, PMT, rot1, G4ThreeVector(0,-0.5*cm-fWrapThickness,0));

	///// Building the ZDS Geometry /////
	G4Tubs * zds = new G4Tubs("zds", 0., fDiameter/2., fScintillatorWidth/2., fStartPhi, fDeltaPhi);
	G4Tubs * pmt = new G4Tubs("pmt", 0., fDiameter/2., fPMTWidth/2., fStartPhi, fDeltaPhi);

	//For placing volume
	rotate = new G4RotationMatrix;

	//Set visual attributes
	G4VisAttributes * plastic_vis = new G4VisAttributes(silver);
	plastic_vis->SetVisibility(true);
	G4VisAttributes * wrapper_vis = new G4VisAttributes(black);
	wrapper_vis->SetVisibility(true);
	G4VisAttributes * zds_vis = new G4VisAttributes(silver);
	zds_vis->SetVisibility(true);
	G4VisAttributes * pmt_vis = new G4VisAttributes(bronze);
	pmt_vis->SetVisibility(true);

	//Names
	G4String nameLog = "TestPlastic";
	G4String nameWrapper = "wrapper";
	G4String namePMT = "TestPMT1";
	G4String nameZDS = "ZDS";
	G4String nameZDSPMT = "zdsWindow";
	//Assign Logical Volume for detectors and wrapping affected by beamline
	fPlasticLog = new G4LogicalVolume(Scint, plasticG4material, nameLog,0,0,0);
	fWrapLog = new G4LogicalVolume(subtractWrap2, wrapG4material, nameWrapper,0,0,0); 
	fPMTLog = new G4LogicalVolume(PMT, PMTG4material, namePMT,0,0,0);
	fZDSLog = new G4LogicalVolume(zds, zdsG4material, nameZDS,0,0,0);
	fZDSPMTLog = new G4LogicalVolume(pmt, PMTG4material, nameZDSPMT,0,0,0);

	//Set Logical Skin for optical photons on wrapping
	G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapper, fWrapLog, ScintWrapper);

	//Give everything colour
	fZDSLog->SetVisAttributes(zds_vis);
	fPMTLog->SetVisAttributes(pmt_vis);
	fPlasticLog->SetVisAttributes(plastic_vis);
	fWrapLog->SetVisAttributes(wrapper_vis);
	fZDSPMTLog->SetVisAttributes(pmt_vis);

	//Place Detectors
	move = G4ThreeVector(0., -length/2. - fRadialDistance, 0.);
	fAssemblyTestPlastics->AddPlacedVolume(fPlasticLog, move, rotate);
	move = G4ThreeVector(0., -length/2. - fRadialDistance, 0.);
	fAssemblyTestPlastics->AddPlacedVolume(fWrapLog, move, rotate);
	move = G4ThreeVector(0., -length-0.2*cm-fRadialDistance, 0.);
	fAssemblyTestPlastics->AddPlacedVolume(fPMTLog, move, rotate);
	move = G4ThreeVector(0., 0., fRadialDistance);
	move.rotateX(-M_PI/2.);
	rotate = new G4RotationMatrix;
	rotate->rotateX(M_PI/2.); // flip the detector so that the face is pointing upstream.
	fAssemblyTestPlastics->AddPlacedVolume(fZDSLog, move, rotate);
	move = G4ThreeVector(0., 0., fRadialDistance+fScintillatorWidth);
	move.rotateX(-M_PI/2.);
	rotate = new G4RotationMatrix;
	rotate->rotateX(M_PI/2.); // flip the detector so that the face is pointing upstream.
	fAssemblyTestPlastics->AddPlacedVolume(fZDSPMTLog, move, rotate);
	
	
	
	///////// 1x1xlength and 2 SiPm ///////////////////
/*
	double length = 3.*cm;
	G4VSolid * Scint = new G4Box("Scint", 1.*cm/2., length/2. , 1.*cm/2.);
	G4VSolid * Wrap_Bigger = new G4Box("Wrap_Bigger", 1.*cm/2.+fWrapThickness, length/2.+fWrapThickness, 1.*cm/2.+fWrapThickness);
	G4SubtractionSolid * subtractWrap = new G4SubtractionSolid("subtractWrap", Wrap_Bigger, Scint, rot1, G4ThreeVector(0,0,0));
	G4VSolid * PMT1 = new G4Box("PMT1", 0.4*cm/2., 0.4*cm/2. , 0.4*cm/2.);
	G4VSolid * PMT2 = new G4Box("PMT2", 0.4*cm/2., 0.4*cm/2. , 0.4*cm/2.);
	// For SiPM on either end ie on square face
	//G4SubtractionSolid * subtractWrap2 = new G4SubtractionSolid("subtractWrap2", subtractWrap, PMT1, rot1, G4ThreeVector(0,-0.5*length-fWrapThickness,0));
	//G4SubtractionSolid * subtractWrap3 = new G4SubtractionSolid("subtractWrap3", subtractWrap2, PMT2, rot1, G4ThreeVector(0,0.5*length+fWrapThickness,0));
	// For SiPM on either end but on same long rectangular face
	G4SubtractionSolid * subtractWrap2 = new G4SubtractionSolid("subtractWrap2", subtractWrap, PMT1, rot1, G4ThreeVector(0, -0.5*length + 0.2*cm, 0.5*cm + fWrapThickness));
	G4SubtractionSolid * subtractWrap3 = new G4SubtractionSolid("subtractWrap3", subtractWrap2, PMT2, rot1, G4ThreeVector(0, 0.5*length - 0.2*cm, 0.5*cm + fWrapThickness));


	//For placing volume
	rotate = new G4RotationMatrix;

	//Set visual attributes
	G4VisAttributes * plastic_vis = new G4VisAttributes(silver);
	plastic_vis->SetVisibility(true);
	G4VisAttributes * wrapper_vis = new G4VisAttributes(black);
	wrapper_vis->SetVisibility(true);
	G4VisAttributes * pmt_vis = new G4VisAttributes(bronze);
	pmt_vis->SetVisibility(true);

	//Names
	G4String nameLog = "TestPlastic";
	G4String nameWrapper = "wrapper";
	G4String namePMT1 = "TestPMT1";
	G4String namePMT2 = "TestPMT2";
	//Assign Logical Volume for detectors and wrapping affected by beamline
	fPlasticLog = new G4LogicalVolume(Scint, plasticG4material, nameLog,0,0,0);
	fWrapLog = new G4LogicalVolume(subtractWrap3, wrapG4material, nameWrapper,0,0,0); 
	fPMT1Log = new G4LogicalVolume(PMT1, PMTG4material, namePMT1,0,0,0);
	fPMT2Log = new G4LogicalVolume(PMT2, PMTG4material, namePMT2,0,0,0);

	//Set Logical Skin for optical photons on wrapping
	G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapper, fWrapLog, ScintWrapper);

	//Give everything colour
	fPMT1Log->SetVisAttributes(pmt_vis);
	fPMT2Log->SetVisAttributes(pmt_vis);
	fPlasticLog->SetVisAttributes(plastic_vis);
	fWrapLog->SetVisAttributes(wrapper_vis);

	//Place Detectors
	//double yOffset = length/2. - 0.5*cm;
	double yOffset = 0.;
	double zOffset = 1.*cm;
	move = G4ThreeVector(0., yOffset, zOffset);
	fAssemblyTestPlastics->AddPlacedVolume(fPlasticLog, move, rotate);
	move = G4ThreeVector(0., yOffset, zOffset);
	fAssemblyTestPlastics->AddPlacedVolume(fWrapLog, move, rotate);

		// For SiPM on either end ie on square face
	//	move = G4ThreeVector(0., yOffset + 0.5*length+0.2*cm, zOffset);
	//	fAssemblyTestPlastics->AddPlacedVolume(fPMT1Log, move, rotate);
	//	move = G4ThreeVector(0., yOffset - 0.5*length-0.2*cm, zOffset);
	//	fAssemblyTestPlastics->AddPlacedVolume(fPMT2Log, move, rotate);
			
	// For SiPM on either end but on same long rectangular face
	move = G4ThreeVector(0., yOffset + 0.5*length - 0.2*cm, zOffset + 0.2*cm + 0.5*cm);
	fAssemblyTestPlastics->AddPlacedVolume(fPMT1Log, move, rotate);
	move = G4ThreeVector(0., yOffset - 0.5*length + 0.2*cm, zOffset + 0.2*cm + 0.5*cm);
	fAssemblyTestPlastics->AddPlacedVolume(fPMT2Log, move, rotate);
*/



	return 1;
}

