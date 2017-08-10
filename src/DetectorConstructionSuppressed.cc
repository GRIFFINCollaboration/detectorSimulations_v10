// Includes Physical Constants and System of Units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//#include "G4FieldManager.hh"
//#include "G4UniformMagField.hh"
//#include "MagneticField.hh"
//#include "G4TransportationManager.hh"
//#include "Field.hh"
//#include "GlobalField.hh"

#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"

#include "DetectionSystem8pi.hh"
#include "DetectionSystemGriffin.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


void DetectorConstruction::DefineSuppressedParameters()
{
    fGriffinFwdBackPosition = 11.0*cm;
    fDetectorRadialDistance = 11.0*cm ;
}

void DetectorConstruction::DefineMaterials()
{ 
    // use G4-NIST materials data base
    //
    G4NistManager* man = G4NistManager::Instance();
    man->FindOrBuildMaterial("G4_Galactic");
    man->FindOrBuildMaterial("G4_Pb");
    man->FindOrBuildMaterial("G4_lAr");
    man->FindOrBuildMaterial("G4_STAINLESS-STEEL");


    man->FindOrBuildMaterial("G4_Al");
    man->FindOrBuildMaterial("G4_POLYETHYLENE");
    man->FindOrBuildMaterial("G4_RUBBER_NEOPRENE");
    man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    man->FindOrBuildMaterial("G4_BGO");
    man->FindOrBuildMaterial("G4_CESIUM_IODIDE");
    man->FindOrBuildMaterial("G4_Ge");
    man->FindOrBuildMaterial("G4_Cu");
    man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    man->FindOrBuildMaterial("G4_AIR");
    G4Material* SODIUM_IODIDE   = man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    G4Material* Tl              = man->FindOrBuildMaterial("G4_Tl");

    G4double a, z, density, temperature, pressure;
    G4String name, symbol;
    G4int    nelements, ncomponents, natoms;
    G4double fractionmass;

    std::vector<G4Element*>  myElements;    // save pointers here to avoid
    std::vector<G4Material*> myMaterials;   // warnings of unused variables

    // Elements
    G4Element* elH = new G4Element(name="H", symbol="H", z=1., a = 1.00*g/mole);
    myElements.push_back(elH);

    G4Element* elC = new G4Element(name="C", symbol="C", z=6., a = 12.01*g/mole);
    myElements.push_back(elC);

    G4Element* elN  = new G4Element(name="N", symbol="N",  z=7.,  a= 14.00674*g/mole);
    myElements.push_back(elN);

    G4Element* elO  = new G4Element(name="O",   symbol="O",  z=8.,  a= 15.9994 *g/mole);
    myElements.push_back(elO);

    G4Element* elAr  = new G4Element(name="Ar",   symbol="Ar",  z=18.,  a= 39.948 *g/mole);
    myElements.push_back(elAr);

    G4Element* elLa = new G4Element(name="La", symbol="La", z=57., 138.9055*g/mole);
    myElements.push_back(elLa);

    G4Element* elBr = new G4Element(name="Br", symbol="Br", z=35., 79.904*g/mole);
    myElements.push_back(elBr);

    G4Element* elGe = new G4Element(name="Ge", symbol="Ge", z=32., 72.64*g/mole);
    myElements.push_back(elGe);

    G4Element* elI = new G4Element(name="I", symbol="I", z=53., 126.90*g/mole);
    myElements.push_back(elI);

    G4Element* elCs = new G4Element(name="Cs", symbol="Cs", z=55., 132.91*g/mole);
    myElements.push_back(elCs);

    G4Element* elTa = new G4Element(name="Ta", symbol="Ta", z=73., 180.95*g/mole);
    myElements.push_back(elTa);

    G4Element* elW = new G4Element(name="W", symbol="W", z=74., 183.84*g/mole);
    myElements.push_back(elW);

    G4Element* elBi = new G4Element(name="Bi", symbol="Bi", z=83., 208.98*g/mole);
    myElements.push_back(elBi);

    G4Element* elCe = new G4Element(name="Ce", symbol="Ce", z=58., 140.116*g/mole);
    myElements.push_back(elCe);

    G4Element* elNi = new G4Element(name="Ni", symbol="Ni", z=28., 58.69*g/mole);
    myElements.push_back(elNi);

    G4Element* elCu = new G4Element(name="Cu", symbol="Cu", z=29., 63.546*g/mole);
    myElements.push_back(elCu);

    G4Element* elNd = new G4Element(name="Nd", symbol="Nd", z=60., 144.242*g/mole);
    myElements.push_back(elNd);

    G4Element* elFe = new G4Element(name="Fe", symbol="Fe", z=26., 55.845*g/mole);
    myElements.push_back(elFe);

    G4Element* elB = new G4Element(name="B", symbol="B", z=5., 10.811*g/mole);
    myElements.push_back(elB);

    G4Element* elNa = new G4Element(name="Na", symbol="Na", z=11., 22.99*g/mole);
    myElements.push_back(elNa);

    G4Element* elLi = new G4Element(name="Li", symbol="Li", z=3. , a=  6.94   *g/mole);
    myElements.push_back(elLi);

    G4Element* elF  = new G4Element(name="F" , symbol="F" , z=9. , a= 18.9984 *g/mole);
    myElements.push_back(elF);

    G4Element* elSi = new G4Element(name="Si", symbol="Si", z=14., a= 28.0855 *g/mole);
    myElements.push_back(elSi);

    G4Element* elCr = new G4Element(name="Cr", symbol="Cr", z=24., a= 51.9961 *g/mole);
    myElements.push_back(elCr);

    G4Element* elPb = new G4Element(name="Pb", symbol="Pb", z=82., a=207.2    *g/mole);
    myElements.push_back(elPb);

    // Materials
    //  density     = universe_mean_density; //from PhysicalConstants.h
    //  pressure    = 1.e-19*pascal;
    //  temperature = 0.1*kelvin;
    //  G4Material* Galactic = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density, kStateGas,temperature,pressure);
    //  myMaterials.push_back(Galactic);

    G4Material* Air = new G4Material(name="Air", density=1.29*mg/cm3, nelements=2);
    Air->AddElement(elN, .7);
    Air->AddElement(elO, .3);
    myMaterials.push_back(Air);

    //  G4Material* TechVacuum = new G4Material( "TechVacuum", density=1.e-5*g/cm3, 1, kStateGas, temperature=STP_Temperature, pressure=2.e-2*bar );
    //  TechVacuum->AddMaterial( Air, 1. );
    //  myMaterials.push_back(TechVacuum);

    //  G4Material* Vacuum = new G4Material(name="Vacuum", z=1., a= 1.01*g/mole, density= universe_mean_density, kStateGas, temperature = 0.1*kelvin, pressure=1.0e-19*pascal);
    //  myMaterials.push_back(Vacuum); // "Galatic" Vacuum

    //universe_mean_density = 1.e-25*g/cm3;

    density     = 1.e-24*g/cm3;
    G4Material* VacuumDensityYoctogramPerCm3 = new G4Material(name="VacuumDensityYoctogramPerCm3", z=1., a=1.01*g/mole, density);
    myMaterials.push_back(VacuumDensityYoctogramPerCm3);

    density     = 1.e-9*g/cm3;
    G4Material* VacuumDensityNanogramPerCm3 = new G4Material(name="VacuumDensityNanogramPerCm3", z=1., a=1.01*g/mole, density);
    myMaterials.push_back(VacuumDensityNanogramPerCm3);

    density     = 10.e-9*g/cm3;
    G4Material* VacuumDensity10nanogramPerCm3 = new G4Material(name="VacuumDensity10nanogramPerCm3", z=1., a=1.01*g/mole, density);
    myMaterials.push_back(VacuumDensity10nanogramPerCm3);

    density     = 100.e-9*g/cm3;
    G4Material* VacuumDensity100nanogramPerCm3 = new G4Material(name="VacuumDensity100nanogramPerCm3", z=1., a=1.01*g/mole, density);
    myMaterials.push_back(VacuumDensity100nanogramPerCm3);

    density     = 1.e-6*g/cm3;
    G4Material* VacuumDensityMicrogramPerCm3 = new G4Material(name="VacuumDensityMicrogramPerCm3", z=1., a=1.01*g/mole, density);
    myMaterials.push_back(VacuumDensityMicrogramPerCm3);

    density     = 1.e-3*g/cm3;
    G4Material* VacuumDensityMilligramPerCm3 = new G4Material(name="VacuumDensityMilligramPerCm3", z=1., a=1.01*g/mole, density);
    myMaterials.push_back(VacuumDensityMilligramPerCm3);

    density     = 1.e-2*g/cm3;
    G4Material* VacuumDensityCentigramPerCm3 = new G4Material(name="VacuumDensityCentigramPerCm3", z=1., a=1.01*g/mole, density);
    myMaterials.push_back(VacuumDensityCentigramPerCm3);

    density     = 1.e-1*g/cm3;
    G4Material* VacuumDensityDecigramPerCm3 = new G4Material(name="VacuumDensityDecigramPerCm3", z=1., a=1.01*g/mole, density);
    myMaterials.push_back(VacuumDensityDecigramPerCm3);

    density     = 1.0*g/cm3;
    G4Material* VacuumDensityGramPerCm3 = new G4Material(name="VacuumDensityGramPerCm3", z=1., a=1.01*g/mole, density);
    myMaterials.push_back(VacuumDensityGramPerCm3);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // From http://geant4.cern.ch/support/source/geant4/examples/advanced/purgingMagnet/src/PurgMagDetectorConstruction.cc
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Laboratory vacuum: Dry air (average composition)
    density = 1.7836*mg/cm3 ;       // STP
    G4Material* Argon = new G4Material(name="Argon", density, ncomponents=1);
    Argon->AddElement(elAr, 1);

    density = 1.25053*mg/cm3 ;       // STP
    G4Material* Nitrogen = new G4Material(name="N2", density, ncomponents=1);
    Nitrogen->AddElement(elN, 2);

    density = 1.4289*mg/cm3 ;       // STP
    G4Material* Oxygen = new G4Material(name="O2", density, ncomponents=1);
    Oxygen->AddElement(elO, 2);

    // LaboratoryVacuum
    density  = 1.2928*mg/cm3 ;       // STP
    density *= 1.0e-8 ;              // pumped vacuum
    temperature = STP_Temperature;
    pressure = 1.0e-8*STP_Pressure;
    G4Material* Vacuum = new G4Material(name="Vacuum",
                                        density,ncomponents=3,
                                        kStateGas,temperature,pressure);
    Vacuum->AddMaterial( Nitrogen, fractionmass = 0.7557 ) ;
    Vacuum->AddMaterial( Oxygen, fractionmass = 0.2315 ) ;
    Vacuum->AddMaterial( Argon,fractionmass = 0.0128 ) ;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    G4Material* Water = new G4Material(name="Water", density=1000*kg/m3, nelements=2);
    Water->AddElement(elH, natoms=2);
    Water->AddElement(elO, natoms=1);
    myMaterials.push_back(Water);

    G4Material* Al = new G4Material(name="Aluminum", z=13., a= 26.98154*g/mole, density= 2.70  *g/cm3);
    myMaterials.push_back(Al);

    G4Material* W = new G4Material(name="Tungsten", z=74., a= 183.84*g/mole, density= 19.25*g/cm3);
    myMaterials.push_back(W);

//	G4Material* Cu = new G4Material("Copper" , 29., 63.550*g/mole, 8.96*g/cm3);
//	mymaterials.push_back(Cu);



    G4Material* Si = new G4Material(name="Silicon", z=14., a= 28.0855*g/mole, density= 2.330  *g/cm3);
    myMaterials.push_back(Si);

    G4Material* Ti = new G4Material(name="Titanium", z=22., 47.867*g/mole, 4.54*g/cm3); //
    myMaterials.push_back(Ti);

    G4Material* Sn = new G4Material(name="Tin", z=50., 118.71*g/mole, 6.99*g/cm3);//
    myMaterials.push_back(Sn);

    G4Material* Au = new G4Material(name="Gold", z=79., a= 196.9666*g/mole, density= 19.30  *g/cm3);
    myMaterials.push_back(Au);

    G4Material* Pb = new G4Material(name="Lead", z=82., a= 207.19*g/mole, density= 11.35  *g/cm3);
    myMaterials.push_back(Pb);

    G4Material* Ge = new G4Material(name="Germanium", z=32., a= 72.64*g/mole, density=5.323 *g/cm3);
    myMaterials.push_back(Ge);

    G4Material* LaBr = new G4Material(name="LanthanumBromide", density=5.08*g/cm3, nelements=2);
    LaBr->AddElement(elLa, natoms=1);
    LaBr->AddElement(elBr, natoms=3);
    myMaterials.push_back(LaBr);

    G4Material* LaBrCe = new G4Material(name="CeriumDopedLanthanumBromide", density=5.08*g/cm3, nelements=2);
    LaBrCe->AddMaterial(LaBr, 95*perCent);
    LaBrCe->AddElement(elCe, 5*perCent);
    myMaterials.push_back(LaBrCe);

    G4Material* Hevimetal = new G4Material("Hevimetal", density=19.0*g/cm3, nelements=3);
    Hevimetal->AddElement(elTa, 80*perCent);
    Hevimetal->AddElement(elNi, 13*perCent);
    Hevimetal->AddElement(elCu, 7*perCent);
    myMaterials.push_back(Hevimetal);

    G4Material* LN2 = new G4Material("LiquidN2", density=804*kg/m3, ncomponents=1);
    LN2->AddElement(elN, natoms=2);
    myMaterials.push_back(LN2);

    G4Material* Delrin = new G4Material("Delrin", density = 1.4*g/cm3, ncomponents=3);
    Delrin->AddElement(elC, natoms=1);
    Delrin->AddElement(elO, natoms=1);
    Delrin->AddElement(elH, natoms=2);
    myMaterials.push_back(Delrin);

    G4Material* Be = new G4Material("Beryllium", z=4., a=9.012*g/mole, density=1848*kg/m3);
    myMaterials.push_back(Be);

    G4Material* Cu = new G4Material("Copper", z=29., a = 63.546*g/mole, density = 8960*kg/m3);
    myMaterials.push_back(Cu);

    G4Material* Zn = new G4Material("Zinc"   , 30., 65.409*g/mole, 7.14*g/cm3);
    myMaterials.push_back(Zn);

    G4Material* Brass= new G4Material("Brass", 8.5*g/cm3, 2);
    Brass->AddMaterial(Cu  , 70*perCent);
    Brass->AddMaterial(Zn , 30*perCent);
    myMaterials.push_back(Brass);

    G4Material* BGO = new G4Material("BGO", density = 7.3*g/cm3, ncomponents=3);
    BGO->AddElement(elBi, natoms=4);
    BGO->AddElement(elGe, natoms=3);
    BGO->AddElement(elO, natoms=12);
    myMaterials.push_back(BGO);

    G4Material* CsI = new G4Material("CesiumIodide", density = 4.51*g/cm3, ncomponents=2);
    CsI->AddElement(elCs, natoms=1);
    CsI->AddElement(elI, natoms=1);
    myMaterials.push_back(CsI);

    G4Material* BC404 = new G4Material("BC404", density = 1.032*g/cm3, ncomponents=2);
    BC404->AddElement(elH, 52.4*perCent);
    BC404->AddElement(elC, 47.6*perCent);
    myMaterials.push_back(BC404);

    G4Material* Mylar = new G4Material("Mylar", density = 1.397*g/cm3, ncomponents=3);
    Mylar->AddElement(elC, natoms=10);
    Mylar->AddElement(elH, natoms=8);
    Mylar->AddElement(elO, natoms=4);
    myMaterials.push_back(Mylar);

    G4Material* NdFeB = new G4Material("NdFeB", density = 7.45*g/cm3, ncomponents=3);  //From Wikipedia
    NdFeB->AddElement(elNd, natoms=2);
    NdFeB->AddElement(elFe, natoms=14);
    NdFeB->AddElement(elB, natoms=1);
    myMaterials.push_back(NdFeB);

    G4Material* Teflon = new G4Material("Teflon", density = 2.2*g/cm3, ncomponents=2);
    Teflon->AddElement(elC, natoms=2);
    Teflon->AddElement(elF, natoms=4);
    myMaterials.push_back(Teflon);

    G4Material* SiLi = new G4Material("Si(Li)", density = 2.330*g/cm3, ncomponents=2);
    SiLi->AddElement(elSi, 99.999*perCent);
    SiLi->AddElement(elLi,  0.001*perCent);
    myMaterials.push_back(SiLi);

    G4Material* Plate = new G4Material("Plate", density = 55.84*g/cm3, ncomponents=3);
    Plate->AddElement(elFe, 90.*perCent);
    Plate->AddElement(elCr, 9. *perCent);
    Plate->AddElement(elPb, 1. *perCent);
    myMaterials.push_back(Plate);

    G4Material* NaI = new G4Material("NaI", density = 3.67*g/cm3, ncomponents=2);
    NaI->AddElement(elNa, natoms=1);
    NaI->AddElement(elI,  natoms=1);
    NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);
    myMaterials.push_back(NaI);

    G4Material* G4SODIUMIODIDEPLUSTl = new G4Material( "G4SODIUMIODIDEPLUSTl", density=3.667*g/cm3, ncomponents = 2);
    G4SODIUMIODIDEPLUSTl->AddMaterial( SODIUM_IODIDE, 99.9*perCent );
    G4SODIUMIODIDEPLUSTl->AddMaterial( Tl, 0.1*perCent );
    myMaterials.push_back(G4SODIUMIODIDEPLUSTl);

    // Deuterated Scintillator from Joey
    //deuterium
    G4Isotope* De = new G4Isotope("De", 1, 2, a = 2.02*g/mole);
    G4Element* D = new G4Element(name="Deuterium",symbol="D",ncomponents=1);
    D->AddIsotope(De,fractionmass=100.*perCent);

    //DeuScin
    G4Material* DeuScin = new G4Material("Deuterated Scintillator", density = 0.954*g/cm3, 3);
    DeuScin->AddElement(elH,fractionmass=0.0625*perCent);
    DeuScin->AddElement(elC,fractionmass=85.7326*perCent);
    DeuScin->AddElement(D,fractionmass=14.2049*perCent);
    myMaterials.push_back(DeuScin);

    G4Material* Peek = new G4Material("Peek",density = 1.26*g/cm3, ncomponents=3);
  Peek->AddElement(elC, natoms=19);
  Peek->AddElement(elH, natoms=12);
  Peek->AddElement(elO, natoms=3);
  myMaterials.push_back(Peek);

  G4Material* WTa = new G4Material("WTa",density = 18677*kg/m3, ncomponents=2 );
  WTa->AddElement(elW, 80.*perCent);
  WTa->AddElement(elTa, 20.*perCent);
  myMaterials.push_back(WTa);

  G4Material* Kapton = new G4Material("Kapton",density = 1.43*g/cm3, ncomponents=4);
  Kapton->AddElement(elC,natoms=22);
  Kapton->AddElement(elH,natoms=10);
  Kapton->AddElement(elN,natoms=2);
  Kapton->AddElement(elO,natoms=5);
  myMaterials.push_back(Kapton);

    //G4cout << *(G4Material::GetMaterialTable()) << G4endl; //outputs material table to terminal
}
