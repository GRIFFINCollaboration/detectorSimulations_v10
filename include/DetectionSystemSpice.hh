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

#ifndef DETECTIONSYSTEMSPICE_HH
#define DETECTIONSYSTEMSPICE_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#ifndef AL_COL
#define AL_COL 0.5, 0.5, 0.5
#endif
#define PEEK_COL 0.5, 0.5, 0.0

class DetectionSystemSpice
{
public:
	DetectionSystemSpice();
	~DetectionSystemSpice();

	//------------------------------------------------//
	// logical and physical volumes
	//------------------------------------------------//
public:
	G4int Build();
	G4int PlaceDetector(G4LogicalVolume* ExpHallLog, G4int nRings);
	G4int PlaceGuardRing(G4LogicalVolume* ExpHallLog);
	void PlaceDetectorMount(G4LogicalVolume* ExpHallLog); 
	void PlaceAnnularClamps(G4LogicalVolume* ExpHallLog);

private:
	G4AssemblyVolume* fAssembly;
	G4AssemblyVolume* fAssemblySiRing[10];

	G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

	G4LogicalVolume* fDetectorMountLog;
	G4LogicalVolume* fAnnularClampLog;  
	G4LogicalVolume* fSiInnerGuardRingLog;
	G4LogicalVolume* fSiOuterGuardRingLog;
	G4LogicalVolume* fSiDetSpiceRingLog[10];

	G4VPhysicalVolume* fDetectorMountPhys;
	G4VPhysicalVolume* fAnnularClampPhys;

	//--------------------------------------------------------//
	// SPICE physical properties
	// OBS: crystal properties are public, others are private
	//--------------------------------------------------------//
	G4String fWaferMaterial;
	G4String fDetectorMountMaterial;
	G4String fAnnularClampMaterial;

	// ----------------------------
	// Dimensions of Detector Mount
	// ----------------------------
	G4double fDetectorMountLength;
	G4double fDetectorMountWidth;
	G4double fDetectorMountThickness;
	G4double fDetectorMountInnerRadius;
	G4double fDetectorMountLipRadius;
	G4double fDetectorMountLipThickness;
	G4double fDetectorMountAngularOffset;
	G4double fDetectorToTargetDistance;
	G4double fDetectorThickness;

	// ---------------------------
	// Dimensions of Annular Clamp
	// ---------------------------
	G4double fAnnularClampThickness;
	G4double fAnnularClampLength;
	G4double fAnnularClampWidth;
	G4double fAnnularClampPlaneOffset;

	//-----------------------------
	// copy numbers
	//-----------------------------
	G4int fDetectorMountCopyNumber;
	G4int fAnnularClampCopyNumber;

	//-----------------------------//
	// parameters for the annular  //
	// planar detector crystal     //
	//-----------------------------//
	G4double 	fSiDetCrystalThickness;
	G4double 	fSiDetCrystalOuterDiameter;
	G4double 	fSiDetCrystalInnerDiameter;
	G4double 	fSiDetRadialSegments;
	G4double 	fSiDetPhiSegments;

	//-------------------------------//
	// parameters for the guard ring //
	//-------------------------------//
	G4double 	fSiDetGuardRingInnerDiameter;
	G4double 	fSiDetGuardRingOuterDiameter;

	//------------------------------------------------//
	// internal methods in Build()
	//------------------------------------------------//
private:
	G4int 	BuildSiliconWafer(G4int ringID);
	G4int 	BuildInnerGuardRing();
	G4int 	BuildOuterGuardRing();
	void 	BuildDetectorMount();
	void 	BuildAnnularClamps();   

	G4Tubs*	BuildCrystal(G4int myRingID);
};

#endif
