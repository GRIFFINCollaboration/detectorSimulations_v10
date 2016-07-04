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
// $Id: DetectionSystem8pi.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectionSystem8pi_h
#define DetectionSystem8pi_h 1

#include "DetectionSystem8pi.hh"
#include "globals.hh"

#include "G4AssemblyVolume.hh"

class DetectionSystem8pi
{
public:
    DetectionSystem8pi();
    ~DetectionSystem8pi();

    G4int Build() ;//G4SDManager* mySDman);
    G4int BuildOneDetector();
    G4int PlaceDetector(G4LogicalVolume* expHallLog, G4ThreeVector move, G4RotationMatrix* rotate, G4int detectorNumber);
    //    G4double const GetDetectorLengthOfUnitsCM() {return this->canLengthZ;};

    // Assembly volumes
    G4AssemblyVolume* fAssembly;
    G4AssemblyVolume* fAssemblyGe;
    G4AssemblyVolume* fAssemblyInnerBGO;
    G4AssemblyVolume* fAssemblyOuterLowerBGO;
    G4AssemblyVolume* fAssemblyOuterUpperBGO;

    G4double fCutClearance;

    G4double fDetailViewEndAngle;
    G4double fDistFromOrigin;
    G4double fCrystalDistFromOrigin;
    G4double fCrystalOuterRadius;
    G4double fCrystalInnerRadius;
    G4double fHoleStartingDepth;
    G4double fCrystalLength;
    G4double fDeadLayerThickness;
    G4double fElectrodeRadius;
    G4double fInnerCanThickness;
    G4double fInnerCanExtendsPastCrystal;
    G4double fInnerCanLidThickness;
    G4double fInnerCanLidSeparation;

    G4double fStructureMatCoolingRodRadius;
    G4double fElectrodeMatCoolingRodRadius;
    G4double fElectrodeMatCoolingRodLength;
    G4double fElectrodeMatCoolingRodDistFromInnerLid;
    G4double fCoolingRodCoverClearance;
    G4double fOuterCanLength;
    G4double fOuterCanDistFromInnerCan;
    G4double fOuterCanThickness;
    G4double fOuterCanExtendsPastCrystal;
    G4double fOuterCanInnerRadius;
    G4double fBerylliumDistFromCrystal;
    G4double fBerylliumThickness;
    G4double fBerylliumRadius;
    G4double fInnerBGOAnnulusLength;
    G4double fInnerBGOClearance;
    G4double fInnerBGOInnerRadius;
    G4double fInnerBGOOuterRadius;
    G4double fOuterBGOBottomThickness;
    G4double fOuterBGOTaperHeight;
    G4double fOuterBGOTopOuterRadius;
    G4double fOuterBGOTotalLength;
    G4double fOuterBGOClearance;
    G4double fOuterBGODisplacement;
    G4double fLiquidN2Length;
    G4double fLiquidN2Radius;

    G4double fHevimetalThickness;
    G4double fHevimetalFrontSideLength;
    G4double fHevimetalRearSideLength;
    G4double fAuxMatPlugNegRadius;
    G4double fAuxMatPlugPosRadius;
    G4double fAuxMatLayerThickness;
    G4double fAuxMatLayerFrontSideLength;

    //half z-lengths for various parts
    //used to pass information from logical to physical volumes
    G4double fDeadLayerHalfLengthZ;
    G4double fGermaniumVacuumCoreZHalfLength;
    G4double fLowerElectrodeMatElectrodeZHalfLength;
    G4double fUpperElectrodeMatElectrodeZHalfLength;
    G4double fInnerCageZHalfLength;
    G4double fInnerCageBottomZHalfLength;
    G4double fInnerCageLidZHalfLength;
    G4double fStructureMatCoolingRodZHalfLength;
    G4double fElectrodeMatCoolingRodZHalfLength;
    G4double fCoolingRodCoverZHalfLength;
    G4double fCoolingRodCoverLidZHalfLength;
    G4double fOuterCanSideZHalfLength;
    G4double fOuterCanLidZHalfLength;
    G4double fOuterCanBottomZHalfLength;
    G4double fBerylliumWindowZHalfLength;
    G4double fInnerBGOAnnulusZHalfLength;
    G4double fStructureMatSheathZHalfLength;
    G4double fOuterLowerBGOAnnulusZHalfLength;
    G4double fOuterUpperBGOAnnulusZHalfLength;
    G4double fLiquidN2ZHalfLength;
    G4double fLiquidN2SideZHalfLength;
    G4double fLiquidN2LidZHalfLength;
    G4double fLiquidN2BottomZHalfLength;

    //suppressed:
    G4double fExtraClearance;
    G4String fElectrodeMat;
    G4String fStructureMat;
    G4String fAuxMat;

private:
    // Logical volumes
    //
    G4LogicalVolume* fGermaniumBlockLog;
    G4LogicalVolume* fGermaniumDeadLayerLog;
    G4LogicalVolume* fGermaniumVacuumCoreLog;
    G4LogicalVolume* fLowerElectrodeMatElectrodeLog;
    G4LogicalVolume* fUpperElectrodeMatElectrodeLog;
    G4LogicalVolume* fInnerCage1Log;
    G4LogicalVolume* fInnerCage2Log;
    G4LogicalVolume* fInnerCage3Log;
    G4LogicalVolume* fInnerCage4Log;
    G4LogicalVolume* fInnerCageBottomLog;
    G4LogicalVolume* fInnerCageLidLog;
    G4LogicalVolume* fStructureMatCoolingRodLog;
    G4LogicalVolume* fElectrodeMatCoolingRodLog;
    G4LogicalVolume* fCoolingRodCoverLog;
    G4LogicalVolume* fCoolingRodCoverLidLog;
    G4LogicalVolume* fOuterCanSideLog;
    G4LogicalVolume* fOuterCanLidLog;
    G4LogicalVolume* fOuterCanBottomLog;
    G4LogicalVolume* fBerylliumWindowLog;
    G4LogicalVolume* fInnerBGOAnnulusLog;
    G4LogicalVolume* fStructureMatSheathLog;
    G4LogicalVolume* fOuterLowerBGOAnnulusLog;
    G4LogicalVolume* fOuterUpperBGOAnnulusLog;
    G4LogicalVolume* fLiquidN2Log;
    G4LogicalVolume* fLiquidN2SideLog;
    G4LogicalVolume* fLiquidN2LidLog;
    G4LogicalVolume* fLiquidN2BottomLog;
    G4LogicalVolume* fHevimetalLog;
    G4LogicalVolume* fAuxMatPlugLog;
    G4LogicalVolume* fAuxMatLayerLog;

    //    SensitiveDetector* germaniumBlockSD;
    //    SensitiveDetector* innerBGOAnnulusSD;
    //    SensitiveDetector* outerLowerBGOSD;
    //    SensitiveDetector* outerUpperBGOSD;

    G4ThreeVector GetDirectionXYZ(G4double, G4double);

    G4int AddGermanium();
    G4int AddGermaniumLogical();
    G4int AddGermaniumPhysical();
    G4int AddGermaniumDeadLayer();
    G4int AddGermaniumCore();
    G4int AddElectrodeMatElectrode();
    G4int AddStructureMatCage();
    G4int AddInnerStructureMatLids();
    G4int AddStructureMatCoolingRod();
    G4int AddElectrodeMatCoolingRod();
    G4int AddBerylliumWindow();
    G4int AddOuterStructureMatCan();
    G4int AddCoolingRodCover();
    G4int AddInnerBGOAnnulus();
    G4int AddStructureMatBGOSheath();
    G4int AddOuterBGOAnnulus();
    G4int AddLiquidN2Container();
    G4int AddHevimetalCollimator();
    G4int AddAuxMatPlug();
    G4int AddThinAuxMatLayer();
};

#endif
