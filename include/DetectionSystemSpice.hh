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

#ifndef DetectionSystemSpice_h
#define DetectionSystemSpice_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#define AL_COL 0.5,0.5,0.5
#define PEEK_COL 0.5, 0.5, 0.0

class DetectionSystemSpice
{
public:
  DetectionSystemSpice();
  ~DetectionSystemSpice();
  
  //------------------------------------------------//
  // logical and physical volumes
  //------------------------------------------------//
private:
  G4AssemblyVolume* assembly;
  G4AssemblyVolume* assemblySiRing[10];

  
public:
  G4int Build();
  G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4int nRings);
  G4int PlaceGuardRing(G4LogicalVolume* exp_hall_log);
  void PlaceDetectorMount(G4LogicalVolume* exp_hall_log); 
  void PlaceAnnularClamps(G4LogicalVolume* exp_hall_log);
  
private:
  G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);

  G4LogicalVolume* detector_mount_log;
  G4LogicalVolume* annular_clamp_log;  
  G4LogicalVolume* siInnerGuardRing_log;
  G4LogicalVolume* siOuterGuardRing_log;
  G4LogicalVolume* siDetSpiceRing_log[10];

   G4VPhysicalVolume* detector_mount_phys;
   G4VPhysicalVolume* annular_clamp_phys;
      
  //--------------------------------------------------------//
  // SPICE physical properties
  // OBS: crystal properties are public, others are private
  //--------------------------------------------------------//
private:
  G4String 	wafer_material;
  G4String detector_mount_material;
  G4String annular_clamp_material;
  
  // ----------------------------
  // Dimensions of Detector Mount
  // ----------------------------
  G4double detector_mount_length;
  G4double detector_mount_width;
  G4double detector_mount_thickness;
  G4double detector_mount_inner_radius;
  G4double detector_mount_lip_radius;
  G4double detector_mount_lip_thickness;
  G4double detector_mount_angular_offset;
  G4double detector_to_target_distance;
  G4double detector_thickness;
  
  // ---------------------------
  // Dimensions of Annular Clamp
  // ---------------------------
  G4double annular_clamp_thickness;
  G4double annular_clamp_length;
  G4double annular_clamp_width;
  G4double annular_clamp_plane_offset;

  //-----------------------------
  // copy numbers
  //-----------------------------
  G4int detectorMountCopyNumber;
  G4int annularClampCopyNumber;
      
  //-----------------------------//
  // parameters for the annular  //
  // planar detector crystal     //
  //-----------------------------//
public:
  G4double 	siDetCrystalThickness;
  G4double 	siDetCrystalOuterDiameter;
  G4double 	siDetCrystalInnerDiameter;
  G4double 	siDetRadialSegments;
  G4double 	siDetPhiSegments;
  
  //-------------------------------//
  // parameters for the guard ring //
  //-------------------------------//
private:
  G4double 	siDetGuardRingInnerDiameter;
  G4double 	siDetGuardRingOuterDiameter;
   
    //------------------------------------------------//
    // internal methods in Build()
    //------------------------------------------------//
private:
  G4int 	BuildSiliconWafer(G4int ringID);
  G4int 	BuildInnerGuardRing();
  G4int 	BuildOuterGuardRing();
  void 		BuildDetectorMount();
  void 		BuildAnnularClamps();   
      
  G4Tubs*	BuildCrystal(G4int myRingID);
 	
};

#endif
