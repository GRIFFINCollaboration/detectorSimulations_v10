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
    // detector dimensions  properties
    fScintillatorLength        = 6.*cm;
    fScintillatorHeight        = 6.*cm;
    fScintillatorWidth         = thickness;
    fRadialDistance = 50*cm;
    fLeadShieldThickness = 6.35*mm;
    fSpacing = 0.5*mm; //assuming no lead on DESCANT but taking into account the optical wrapping
    fNumDet = numDet ;
    fPlasticLogArray.resize(numDet, NULL);
    fWrapLogArray.resize(numDet, NULL);
    fPMT1LogArray.resize(numDet, NULL);
    fPMT2LogArray.resize(numDet, NULL);

    fWrapThickness = 0.5 * mm;
 
    fWrapMaterial = "Teflon";
    fPMTMaterial = "G4_SILICON_DIOXIDE";
    
 //   blue=G4Color(0.,0.,1.);
    blue=G4Color(0.,1.,1.);
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
for (int i = 0; i<fNumDet; i++) {
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
ScintWrapper->SetModel(unified);  // unified or glisur
ScintWrapper->SetType(dielectric_dielectric);  // dielectric and dielectric or metal?
ScintWrapper->SetFinish(groundfrontpainted);  // teflon wrapping on polished surface->front/back painted // Teflon should be Lambertian in air, specular in optical grease
//ScintWrapper->SetPolish(0.9);  // specfic to the glisur model
G4MaterialPropertiesTable * ScintWrapperMPT = new G4MaterialPropertiesTable();
G4double rIndex_Teflon[numShort] = {1.35, 1.35, 1.35}; //Taken from wikipedia
ScintWrapperMPT->AddProperty("RINDEX", photonEnergyShort, rIndex_Teflon, nEntriesShort)->SetSpline(true);  //refractive index can change with energy
G4double reflectivity[numShort] = {0.95, 0.95, 0.95};
ScintWrapperMPT->AddProperty("REFLECTIVITY", photonEnergyShort, reflectivity, nEntriesShort);  //
ScintWrapper->SetMaterialPropertiesTable(ScintWrapperMPT);
ScintWrapper->DumpInfo();


//////// Quartz  ////////////
G4MaterialPropertiesTable * QuartzMPT = new G4MaterialPropertiesTable();
G4double rIndex_Quartz[numShort] = {1.474, 1.474, 1.474}; //Taken from Joey github
QuartzMPT->AddProperty("RINDEX", photonEnergyShort, rIndex_Quartz, nEntriesShort)->SetSpline(true);  //refractive index can change with energy
QuartzMPT->AddConstProperty("ABSLENGTH", 40.*cm); //from Joeys github
PMTG4material->SetMaterialPropertiesTable(QuartzMPT);

//Building the Plastic Geometry

//Detector Placement
G4double detWidth = 100.*cm/(fNumDet);
G4double startPos = 50.*cm-detWidth/2.+2*fWrapThickness;

//place outer radius of plastics at position of DESCANT detectors, taking into account lead shield
//G4double outerRadius = fRadialDistance - fLeadShieldThickness - fSpacing;
//place outer radius of plastics at position of DESCANT detectors, without lead shield
G4double outerRadius = fRadialDistance - fSpacing;
G4double innerRadius = outerRadius - fScintillatorWidth;
G4double outerWrap = outerRadius + fWrapThickness;
G4double innerWrap = innerRadius - fWrapThickness;
G4double wrapBoxThick = (detWidth+2.*fWrapThickness)/2.;
//Opening angle of DESCANT is 65.5 degrees, approx 1.143 radians. Limiting to 1.13
G4double startTheta = 0.;
G4double deltaTheta = 1.13;
G4double startPhi = 0. ;
G4double deltaPhi = 2*M_PI;
G4double startThetaPMT = 1.13;
G4double deltaThetaPMT = 0.1;
G4double startPhiPMT1 = 0.;
G4double deltaPhiPMT1 = M_PI;
G4double startPhiPMT2 = M_PI;
G4double deltaPhiPMT2 = M_PI;

//for Plastic Hollow Sphere
//G4Sphere * plasticSphere = new G4Sphere("Plastic Detector", innerRadius, outerRadius, startPhi, deltaPhi, startTheta, deltaTheta);

//For placing volume
move = G4ThreeVector(0., 0., 0.);
rotate = new G4RotationMatrix;
//For intersection solid
G4RotationMatrix *rot1 = new G4RotationMatrix(0,0,0);
//G4cout << "startPos: " << startPos << G4endl;
//G4cout << "detWidth: " << detWidth << G4endl;

//G4Solids
G4VSolid * boxBars = new G4Box("boxBars", detWidth/2., 1.*m , 1.*m);
G4VSolid * SphereBars = new G4Sphere("SphereBars", innerRadius, outerRadius, startPhi, deltaPhi, startTheta, deltaTheta);
G4IntersectionSolid * interSolidBars;
G4VSolid * SphereWrap = new G4Sphere("SphereWrap", innerWrap, outerWrap, startPhi, deltaPhi, startTheta, deltaTheta);
G4VSolid * boxWrap = new G4Box("boxWrap", wrapBoxThick, 1.*m , 1.*m);
G4IntersectionSolid * interSolidWrap;
G4SubtractionSolid * subtractSolidWrap;
G4VSolid * SpherePMT1 = new G4Sphere("SpherePMT1", innerRadius, outerRadius, startPhiPMT1, deltaPhiPMT1, startThetaPMT, deltaThetaPMT);
G4IntersectionSolid * interSolidPMT1;
G4VSolid * SpherePMT2 = new G4Sphere("SpherePMT2", innerRadius, outerRadius, startPhiPMT2, deltaPhiPMT2, startThetaPMT, deltaThetaPMT);
G4IntersectionSolid * interSolidPMT2;
G4double BeamLineXY = 6.5*cm;
G4SubtractionSolid * subtractSolidBeamLine_Bar;
G4SubtractionSolid * subtractSolidBeamLine_Wrap;


//Set visual attributes
G4VisAttributes * plastic_vis = new G4VisAttributes(blue);
plastic_vis->SetVisibility(true);
G4VisAttributes * wrap_vis = new G4VisAttributes(black);
wrap_vis->SetVisibility(true);


//Build array of logical volumes
for (int i= 0 ; i < fNumDet ; i++) {
interSolidBars = new G4IntersectionSolid("interSolidBars", SphereBars, boxBars, rot1, G4ThreeVector(startPos, 0, 50*cm));
G4String name0 = "PlasticDet_";
G4String nameLog=name0+std::to_string(i);

interSolidWrap = new G4IntersectionSolid("interSolidWrap", SphereWrap, boxWrap, rot1, G4ThreeVector(startPos, 0, 50*cm));
subtractSolidWrap = new G4SubtractionSolid("subtractSolidWrap", interSolidWrap, interSolidBars, rot1, G4ThreeVector(0,0,0));

interSolidPMT1 = new G4IntersectionSolid("interSolidPMT1", SpherePMT1, boxWrap, rot1, G4ThreeVector(startPos, 0, 50*cm));
interSolidPMT2 = new G4IntersectionSolid("interSolidPMT2", SpherePMT2, boxWrap, rot1, G4ThreeVector(startPos, 0, 50*cm));

G4String name1 = "wrapper_";
G4String nameWrapper = name1+std::to_string(i);

//Should confirm these are actually top
G4String name2 = "PMT1_top_";
G4String namePMT1 = name2+std::to_string(i);

G4String name3 = "PMT2_bottom_";
G4String namePMT2 = name3+std::to_string(i);

//G4cout << "name: " << name << G4endl;
G4cout << "startPos: " << startPos << G4endl;

//Counter for bottom detectors 
int bCounter;
//Putting in subtraction Beam Line 
if( abs(startPos-detWidth/2.) < BeamLineXY || abs(startPos+detWidth/2.) < BeamLineXY || abs(startPos) < BeamLineXY) {
G4cout << "startPos In loop: " << startPos << G4endl;
//Implemented in a way where no partial subtractions occur.  change wrapBoxThick<->BeamLineXY in x position to revert to partial subtractions.
G4VSolid * boxBeamLine = new G4Box("boxBeamLine", wrapBoxThick, BeamLineXY, 1*m); 
//Implemented in a way where no partial subtractions occur.  change startPos<->0. in x position to revert to partial subtractions.
subtractSolidBeamLine_Bar = new G4SubtractionSolid("subtractSolidBeamLine_Bar", interSolidBars, boxBeamLine, rot1, G4ThreeVector(startPos,0,50*cm));
subtractSolidBeamLine_Wrap = new G4SubtractionSolid("subtractSolidBeamLine_Bar", subtractSolidWrap, boxBeamLine, rot1, G4ThreeVector(startPos,0,50*cm));
fPlasticLogArray[i] = new G4LogicalVolume(subtractSolidBeamLine_Bar, plasticG4material, nameLog,0,0,0);
fWrapLogArray[i] = new G4LogicalVolume(subtractSolidBeamLine_Wrap, wrapG4material, nameWrapper,0,0,0); 
}
else {
fPlasticLogArray[i] = new G4LogicalVolume(interSolidBars, plasticG4material, nameLog,0,0,0);
fWrapLogArray[i] = new G4LogicalVolume(subtractSolidWrap, wrapG4material, nameWrapper,0,0,0);
}
fPMT1LogArray[i] = new G4LogicalVolume(interSolidPMT1, PMTG4material, namePMT1,0,0,0);
fPMT2LogArray[i] = new G4LogicalVolume(interSolidPMT2, PMTG4material, namePMT2,0,0,0);

G4LogicalSkinSurface * Surface = new G4LogicalSkinSurface(nameWrapper, fWrapLogArray[i], ScintWrapper);

fPlasticLogArray[i]->SetVisAttributes(plastic_vis);
fWrapLogArray[i]->SetVisAttributes(wrap_vis);


//build every second detector
//if(i%2==0) {}
fAssemblyPlastics->AddPlacedVolume(fPlasticLogArray[i], move, rotate);
fAssemblyPlastics->AddPlacedVolume(fWrapLogArray[i], move, rotate);
fAssemblyPlastics->AddPlacedVolume(fPMT1LogArray[i], move, rotate);
fAssemblyPlastics->AddPlacedVolume(fPMT2LogArray[i], move, rotate);

startPos = startPos-detWidth-2*fWrapThickness;
//G4cout << "startPos after : " << startPos << G4endl;

}

 
    return 1;
}

