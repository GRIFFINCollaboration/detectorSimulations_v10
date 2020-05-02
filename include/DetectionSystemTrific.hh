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
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectionSystemTrific_HH
#define DetectionSystemTrific_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"



class DetectionSystemTrific
{
public:
    
	DetectionSystemTrific(G4double setpressure=760,G4double setwindow=6*CLHEP::um,G4bool setal=true,G4bool setflat=true,G4double setdegrader=0,G4String degradermat="");
	~DetectionSystemTrific();

	//------------------------------------------------//
	// logical and physical volumes
	//------------------------------------------------//
public:
	G4int BuildPlace(G4LogicalVolume* ExpHallLog);


private:

	G4LogicalVolume* fGasCell;
// 	G4LogicalVolume* fAnnularClampLog;
// 	std::vector<G4LogicalVolume*> fSiliLog;
// 
// 	G4VPhysicalVolume* fDetectorMountPhys;
// 	G4VPhysicalVolume* fAnnularClampPhys;
	G4AssemblyVolume* fGridPCB;

	//--------------------------------------------------------//
	// Trific physical properties
	//--------------------------------------------------------//
	G4String fWireMaterial;
	G4String fPCBMaterial;
	G4String fChamberMaterial;
	G4String fWindowFrameMaterial;
	G4String fWindowMaterial;
	G4String fWindowSurfaceMaterial;
	G4String fGasMaterial;
	G4String fDegraderMaterial;


	// ----------------------------
	// Dimensions of Detector
	// ----------------------------
	G4double fChamberLength;
	G4double fPCBThickness;
	G4double fGridAngle;
	G4double fSpacerZThickness;
	G4double fGridSpacing;
	G4double fWireD;
	G4double fWirePitch;
	G4double fChamberInnerD;
	G4double fChamberOuterD;
	G4double fPCBInnerD;
	G4double fPCBOuterD;
	G4double fPCBCutDepth;
	G4double fPCBCutHalfHeight;
	G4double fRodD;
	G4double fRodPlaceD;

	G4int fYgrid;
	std::vector<G4int> fYwires;
	G4int fXgrid;
	std::vector<G4int> fXwires;
    
    
	G4double fTargetChamberZOffset;
	G4double fChamberWindowZ;
	G4double fWindowGridZ;
	G4double fWindowWasherZ;
	G4double fWindowThickness;
	G4double fWindowCoatingThickness;
	G4double fWindowInnerD;
	G4double fWindowOuterD;
	G4double fWindowChamberLength;
	G4double fWindowPipeInnerD;
	
	G4double fDegraderThickness;
    
    G4bool fFlatWindow;
    G4bool fAluminised;
    
	
	//------------------------------------------------//
	// internal methods in Build()
	//------------------------------------------------//
private:

	void 	BuildPCB();
	void 	PlacePCBs(G4LogicalVolume*);
    G4VSolid* GasCell(int N=1);
	void 	BuildPlaceYSense(G4LogicalVolume*);
	void 	BuildPlaceXSense(G4LogicalVolume*);
	void 	BuildPlaceWindow(G4LogicalVolume*);
	void 	BuildPlaceFlatWindow(G4LogicalVolume*);
	void 	BuildPlacePipe(G4LogicalVolume*);
	void 	BuildGasVolume(G4LogicalVolume*);
// 	void BuildPCB(G4LogicalVolume* ExpHallLog); 
};

#endif
