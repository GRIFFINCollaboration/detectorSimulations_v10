#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polyhedra.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystem8pi.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

//constructor suppressed

DetectionSystem8pi::~DetectionSystem8pi() {
    delete fGermaniumBlockLog;
    delete fGermaniumDeadLayerLog;
    delete fGermaniumVacuumCoreLog;
    delete fLowerElectrodeMatElectrodeLog;
    delete fUpperElectrodeMatElectrodeLog;
    delete fInnerCage1Log;
    delete fInnerCage2Log;
    delete fInnerCage3Log;
    delete fInnerCage4Log;
    delete fInnerCageBottomLog;
    delete fInnerCageLidLog;
    delete fStructureMatCoolingRodLog;
    delete fElectrodeMatCoolingRodLog;
    delete fCoolingRodCoverLog;
    delete fCoolingRodCoverLidLog;
    delete fOuterCanSideLog;
    delete fOuterCanLidLog;
    delete fOuterCanBottomLog;
    delete fBerylliumWindowLog;
    delete fInnerBGOAnnulusLog;
    delete fStructureMatSheathLog;
    delete fOuterLowerBGOAnnulusLog;
    delete fOuterUpperBGOAnnulusLog;
    delete fLiquidN2Log;
    delete fLiquidN2SideLog;
    delete fLiquidN2LidLog;
    delete fLiquidN2BottomLog;
    delete fHevimetalLog;
    delete fAuxMatPlugLog;
    delete fAuxMatLayerLog;
}

G4int DetectionSystem8pi::Build() {

    // Build assembly volumes
    fAssembly = new G4AssemblyVolume();
    fAssemblyGe = new G4AssemblyVolume();
    fAssemblyInnerBGO = new G4AssemblyVolume();
    fAssemblyOuterLowerBGO = new G4AssemblyVolume();
    fAssemblyOuterUpperBGO = new G4AssemblyVolume();

    AddGermanium();

    AddGermaniumDeadLayer();

    AddGermaniumCore();
    //Add electrodeMat electrode inside germanium core

    AddElectrodeMatElectrode();
    //Add the inner structureMat cage around Ge crystal

    AddStructureMatCage();
    //Add inner structureMat can lids

    AddInnerStructureMatLids();
    //Add beryllium window

    AddStructureMatCoolingRod();
    //Add electrodeMat cooling rod section

    AddElectrodeMatCoolingRod();
    //Add the outer can

    AddBerylliumWindow();
    //Add structureMat cooling rod section

    AddOuterStructureMatCan();
    //add structureMat cover around cooling rod

    AddCoolingRodCover();
    //add structureMat sheath around inner BGO annulus

    AddStructureMatBGOSheath();
    //add inner BGO annulus around cooling rod

    AddInnerBGOAnnulus();
    //add outer BGO annulus

    AddOuterBGOAnnulus();
    //add liquid N2 cooling container

    AddLiquidN2Container();
    //add hevimetal collimator & core

    AddHevimetalCollimator();
    //Add AuxMat plug in Hevimetal collimator

    AddAuxMatPlug();
    //Add thin auxMat layer

    AddThinAuxMatLayer();
    //Add hevimetal and auxMat plug

    return 1;
}//end ::Build 

G4int DetectionSystem8pi::PlaceDetector(G4LogicalVolume* expHallLog, G4ThreeVector move, G4RotationMatrix* rotate, G4int detectorNumber) {
    //G4int detectorCopyID = 0;

    //G4int copyNumber = detectorCopyID + detectorNumber;

    fAssemblyGe->MakeImprint(expHallLog, move, rotate, detectorNumber);
    fAssemblyInnerBGO->MakeImprint(expHallLog, move, rotate, detectorNumber+(20*1));
    fAssemblyOuterLowerBGO->MakeImprint(expHallLog, move, rotate, detectorNumber+(20*2));
    fAssemblyOuterUpperBGO->MakeImprint(expHallLog, move, rotate, detectorNumber+(20*3));
    fAssembly->MakeImprint(expHallLog, move, rotate, detectorNumber+(20*4));

    return 1;
}

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystem8pi::GetDirectionXYZ(G4double theta, G4double phi) {
    G4double x,y,z;
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);
    G4ThreeVector direction = G4ThreeVector(x,y,z);
    return direction;
}//end ::end GetDirectionXYZ


//Add the germanium crystal
G4int DetectionSystem8pi::AddGermanium() {
    //material
    G4Material* material = G4Material::GetMaterial("Germanium");
    if(!material) {
        G4cout<<" ----> Material "<<"Germanium"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    //vis attributes
    G4VisAttributes* germaniumBlockVisAtt = new G4VisAttributes(G4Colour(0.0,0.5,0.0));
    germaniumBlockVisAtt->SetVisibility(true);

    // germanium detector measurements
    G4double innerRadius = 0.0*cm;
    G4double outerRadius = fCrystalOuterRadius;
    G4double halfLengthZ = fCrystalLength/2.0;

    //primtive volume
    G4Tubs* germaniumBlock = new G4Tubs("germaniumBlock", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    // Cut Out Dead Layer! ///////////////////////////////
    //measurements
    G4double cutInnerRadius = 0.0*cm;
    G4double cutOuterRadius =  fCrystalInnerRadius + fDeadLayerThickness + fCutClearance;
    G4double cutHalfLengthZ =  fCrystalLength/2.0 - fHoleStartingDepth/2.0 + fDeadLayerThickness/2.0 + fCutClearance/2.0;
    //primitive volume
    G4Tubs* cutGermaniumDeadLayer = new G4Tubs("cutGermaniumDeadLayer", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0.0*deg, fDetailViewEndAngle);

    //positioning
    G4double cutZPosition = fHoleStartingDepth/2.0 - fDeadLayerThickness/2.0;
    G4ThreeVector moveCut = G4ThreeVector(0, 0, cutZPosition);

    G4SubtractionSolid* germaniumBlockWithCavity = new G4SubtractionSolid("germaniumBlockWithCavity", germaniumBlock, cutGermaniumDeadLayer, 0, moveCut);
    ////////////////////////////////////////////////////////

    //positioning
    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fGermaniumBlockLog == nullptr) {
        fGermaniumBlockLog = new G4LogicalVolume(germaniumBlockWithCavity, material, "8piGermaniumBlockLog", 0, 0, 0);
        fGermaniumBlockLog->SetVisAttributes(germaniumBlockVisAtt);
    }

    //physical volume
    //  germaniumBlockPhys = new G4PVPlacement(rotate, move,
    //       "germaniumBlockPhys", germaniumBlockLog, expHallPhys,
    //       false, detectorNumber);

    fAssemblyGe->AddPlacedVolume(fGermaniumBlockLog, move, rotate);

    return 1;
}//end ::end AddGermanium

//Add the germanium crystal dead layer
G4int DetectionSystem8pi::AddGermaniumDeadLayer() {
    //material
    G4Material* material = G4Material::GetMaterial("Germanium");
    if(!material) {
        G4cout<<" ----> Material "<<"Germanium"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    //vis attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,0.5,0.0));
    visAtt->SetVisibility(true);

    //measurements
    G4double innerRadius = 0.0*cm;
    G4double outerRadius =  fCrystalInnerRadius + fDeadLayerThickness;
    G4double halfLengthZ =  fCrystalLength/2.0 - fHoleStartingDepth/2.0 + fDeadLayerThickness/2.0;
    //primitive volume
    G4Tubs* germaniumDeadLayer = new G4Tubs("germaniumDeadLayer", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    // Cut Out Vacuum Layer! ///////////////////////////////
    //measurements
    G4double cutInnerRadius = 0.0*cm;
    G4double cutOuterRadius = fCrystalInnerRadius + fCutClearance;
    G4double cutHalfLengthZ = ((fCrystalLength - fHoleStartingDepth) + fCutClearance)/2.0;

    //primitive volume
    G4Tubs* cutGermaniumVacuumCore = new G4Tubs("cutGermaniumVacuumCore", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0.0*deg, fDetailViewEndAngle);

    //positioning
    G4double cutZPosition = fDeadLayerThickness/2.0;
    G4ThreeVector moveCut = G4ThreeVector(0, 0, cutZPosition);

    G4SubtractionSolid* germaniumDeadLayerWithCavity = new G4SubtractionSolid("germaniumDeadLayerWithCavity", germaniumDeadLayer, cutGermaniumVacuumCore, 0, moveCut);
    ////////////////////////////////////////////////////////

    //positioning
    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalLength/2.0 - halfLengthZ + fCrystalDistFromOrigin; // adjust for placement
    G4ThreeVector move = G4ThreeVector(0, 0, zPosition);

    //logical volume
    if(fGermaniumDeadLayerLog == nullptr) {
		fGermaniumDeadLayerLog = new G4LogicalVolume(germaniumDeadLayerWithCavity, material, "germaniumDeadLayerLog", 0, 0, 0);
		fGermaniumDeadLayerLog->SetVisAttributes(visAtt);
    }

    //physical volume
    //  germaniumDeadLayerPhys = new G4PVPlacement(0, move,
    //       germaniumDeadLayerLog, "germaniumDeadLayerPhys",  germaniumBlockLog,
    //       false, detectorNumber);

    fAssembly->AddPlacedVolume(fGermaniumDeadLayerLog, move, rotate);

    return 1;
}//end ::end AddGermaniumDeadLayer


//Add the logical germanium crystal vacuum core
G4int DetectionSystem8pi::AddGermaniumCore() {
    //material
    G4Material* material = G4Material::GetMaterial("Vacuum");
    if(!material) {
        G4cout<<" ----> Material "<<"Vacuum"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    //vis attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    visAtt->SetVisibility(true);

    //measurements
    G4double innerRadius = 0.0*cm;
    G4double outerRadius = fCrystalInnerRadius;
    G4double halfLengthZ = (fCrystalLength - fHoleStartingDepth)/2.0;

    //primitive volume
    G4Tubs* germaniumVacuumCore = new G4Tubs("germaniumVacuumCore", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    // Cut Out ElectrodeMat Electrode! ///////////////////////////////
    //measurements
    G4double cutInnerRadius = 0.0*cm;
    G4double cutOuterRadius = fElectrodeRadius + fCutClearance/2.0;
    G4double cutHalfLengthZ =  fCrystalLength/2.0 - fHoleStartingDepth/2.0 + fCutClearance/2.0;
    //primitive volume
    G4Tubs* cutGermaniumVacuumCore = new G4Tubs("cutGermaniumVacuumCore", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0.0*deg, fDetailViewEndAngle);

    //positioning
    G4double cutZPosition = 0.0;
    G4ThreeVector moveCut = G4ThreeVector(0, 0, cutZPosition);

    G4SubtractionSolid* germaniumVacuumCoreWithCavity = new G4SubtractionSolid("germaniumVacuumCoreWithCavity", germaniumVacuumCore, cutGermaniumVacuumCore, 0, moveCut);
    ////////////////////////////////////////////////////////

    //positioning
    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    //G4double zPosition = fDeadLayerThickness/2.0 + fCrystalLength/2.0 - halfLengthZ + fCrystalDistFromOrigin; // adjust for placement
    G4double zPosition = fCrystalLength/2.0 - halfLengthZ + fCrystalDistFromOrigin; // adjust for placement
    G4ThreeVector move = G4ThreeVector(0, 0, zPosition);

    //logical volume
    if(fGermaniumVacuumCoreLog == nullptr) {
		fGermaniumVacuumCoreLog = new G4LogicalVolume(germaniumVacuumCoreWithCavity, material, "germaniumVacuumCoreLog", 0, 0, 0);
		fGermaniumVacuumCoreLog->SetVisAttributes(visAtt);
    }

    //physical volume
    //  germaniumVacuumCorePhys = new G4PVPlacement(0, move,
    //       germaniumVacuumCoreLog, "germaniumVacuumCorePhys", germaniumDeadLayerLog,
    //       false, detectorNumber);

    fAssembly->AddPlacedVolume(fGermaniumVacuumCoreLog, move, rotate);

    return 1;
}//end ::end AddGermaniumCore


//Add the electrodeMat electrode inside the germanium crystal core
G4int DetectionSystem8pi::AddElectrodeMatElectrode() {
    //material
    G4Material* material = G4Material::GetMaterial(fElectrodeMat);
    if(!material) {
        G4cout<<" ----> Material "<<fElectrodeMat<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    //vis attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.75, 0.45, 0.2));
    visAtt->SetVisibility(true);

    //lower portion of electrodeMat electrode
    //measurements
    G4double innerRadius = 0.0*cm;
    G4double outerRadius = fElectrodeRadius;
    G4double halfLengthZ =  fCrystalLength/2.0
            - fHoleStartingDepth/2.0;
    //primitive volume
    G4Tubs* lowerElectrodeMatCylinder = new G4Tubs("lowerElectrodeMatCylinder", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    //positioning
    G4double lowerZPosition = fCrystalDistFromOrigin + fHoleStartingDepth/2.0; // check this! // adjust position

    //logical volume
    if(fLowerElectrodeMatElectrodeLog == nullptr) {
        fLowerElectrodeMatElectrodeLog = new G4LogicalVolume(lowerElectrodeMatCylinder, material, "lowerElectrodeMatElectrodeLog", 0, 0, 0);
        fLowerElectrodeMatElectrodeLog->SetVisAttributes(visAtt);
    }

    //upper portion of electrodeMat electrode
    //measurements
    halfLengthZ =  fInnerCanExtendsPastCrystal/2.0;

    //primitive volume
    G4Tubs* upperElectrodeMatCylinder = new G4Tubs("upperElectrodeMatCylinder", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    //positioning
    G4double zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + halfLengthZ;

    //logical volume
    if(fUpperElectrodeMatElectrodeLog == nullptr) {
        fUpperElectrodeMatElectrodeLog = new G4LogicalVolume(upperElectrodeMatCylinder, material, "upperElectrodeMatElectrodeLog", 0, 0, 0);
        fUpperElectrodeMatElectrodeLog->SetVisAttributes(visAtt);
    }

    //position vector
    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4ThreeVector moveLower = lowerZPosition * direction;
    G4ThreeVector moveUpper = zPosition * direction;
    G4ThreeVector moveNull = G4ThreeVector(0, 0, 0);


    //physical volumes
    fAssembly->AddPlacedVolume(fLowerElectrodeMatElectrodeLog, moveLower, rotate);
    fAssembly->AddPlacedVolume(fUpperElectrodeMatElectrodeLog, moveUpper, rotate);

    return 1;
}//end ::end AddElectrodeMatElectrode


//Add the cage directly around the germanium crystal
G4int DetectionSystem8pi::AddStructureMatCage() {
    //material
    G4Material* material = G4Material::GetMaterial(fStructureMat);
    if(!material) {
        G4cout<<" ----> Material "<<fStructureMat<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    visAtt->SetVisibility(true);

    // innerCage measurements
    G4double innerRadius = fCrystalOuterRadius + fExtraClearance; //clearance
    G4double outerRadius = innerRadius + fInnerCanThickness;
    G4double halfLengthZ = (fCrystalLength + fInnerCanExtendsPastCrystal)/2.0;

    G4Tubs* innerCage1 = new G4Tubs("innerCage1", innerRadius, outerRadius, halfLengthZ, 0.0*deg, /*fDetailViewEndAngle*/15.0*deg);
    G4Tubs* innerCage2 = new G4Tubs("innerCage2", innerRadius, outerRadius, halfLengthZ, 90.0*deg, /*fDetailViewEndAngle*/15.0*deg);
    G4Tubs* innerCage3 = new G4Tubs("innerCage3", innerRadius, outerRadius, halfLengthZ, 180.0*deg, /*fDetailViewEndAngle*/15.0*deg);
    G4Tubs* innerCage4 = new G4Tubs("innerCage4", innerRadius, outerRadius, halfLengthZ, 270.0*deg, /*fDetailViewEndAngle*/15.0*deg);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin + fInnerCanExtendsPastCrystal/2.0;
    G4ThreeVector move = zPosition * direction;

    //4 logical sides (quadrants)
    if(fInnerCage1Log == nullptr && fInnerCage2Log == nullptr && fInnerCage3Log == nullptr && fInnerCage4Log == nullptr) {
        fInnerCage1Log = new G4LogicalVolume(innerCage1, material, "innerCage1Log", 0, 0, 0);
        fInnerCage1Log->SetVisAttributes(visAtt);
        fInnerCage2Log = new G4LogicalVolume(innerCage2, material, "innerCage2Log", 0, 0, 0);
        fInnerCage2Log->SetVisAttributes(visAtt);
        fInnerCage3Log = new G4LogicalVolume(innerCage3, material, "innerCage3Log", 0, 0, 0);
        fInnerCage3Log->SetVisAttributes(visAtt);
        fInnerCage4Log = new G4LogicalVolume(innerCage4, material, "innerCage4Log", 0, 0, 0);
        fInnerCage4Log->SetVisAttributes(visAtt);
    }

    //4 physical sides
    fAssembly->AddPlacedVolume(fInnerCage1Log, move, rotate);
    fAssembly->AddPlacedVolume(fInnerCage2Log, move, rotate);
    fAssembly->AddPlacedVolume(fInnerCage3Log, move, rotate);
    fAssembly->AddPlacedVolume(fInnerCage4Log, move, rotate);

    //add bottom to innerCan
    innerRadius = /*0.0*mm*/ outerRadius - 2.0*mm; //2.0mm 'lip'
    halfLengthZ = fInnerCanThickness/2.0;

    G4Tubs* innerCageBottom = new G4Tubs("innerCageBottom", innerRadius,
                                           outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle
												);

    zPosition = fCrystalDistFromOrigin - fCrystalLength/2.0 - halfLengthZ;
    move =  zPosition * direction;

    //logical volume
    if(fInnerCageBottomLog == nullptr) {
        fInnerCageBottomLog = new G4LogicalVolume(innerCageBottom, material, "innerCageBottomLog", 0, 0, 0);
        fInnerCageBottomLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fInnerCageBottomLog, move, rotate);

    return 1;
}//end ::end AddInnerStructureMatCage

G4int DetectionSystem8pi::AddInnerStructureMatLids() {
    //material
    G4Material* material = G4Material::GetMaterial(fStructureMat);
    if(!material) {
        G4cout<<" ----> Material "<<"StructureMat"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    visAtt->SetVisibility(true);

    //ring measurements
    G4double innerRadius = fStructureMatCoolingRodRadius;
    G4double outerRadius = fCrystalOuterRadius + fInnerCanThickness + fExtraClearance; //extra clearance
    G4double halfLengthZ = fInnerCanLidThickness/2.0;

    G4Tubs* lid = new G4Tubs("lid", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fInnerCanExtendsPastCrystal + halfLengthZ;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fInnerCageLidLog == nullptr) {
        fInnerCageLidLog = new G4LogicalVolume(lid, material, "innerCageLidLog", 0, 0, 0);
        fInnerCageLidLog->SetVisAttributes(visAtt);
    }

    //place 1st lid
    fAssembly->AddPlacedVolume(fInnerCageLidLog, move, rotate);

    //separate 2 lids
    zPosition =  zPosition + fInnerCanLidSeparation + fInnerCanLidThickness;
    move =  zPosition * direction;

    //place 2nd lid
    fAssembly->AddPlacedVolume(fInnerCageLidLog, move, rotate);

    return 1;
}//end ::end AddInnerStructureMatLids

G4int DetectionSystem8pi::AddStructureMatCoolingRod() {
    //material
    G4Material* material = G4Material::GetMaterial(fStructureMat);
    if(!material) {
        G4cout<<" ----> Material "<<fStructureMat<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    visAtt->SetVisibility(true);

    //ring measurements
    G4double innerRadius = 0.0*mm;
    G4double outerRadius = fStructureMatCoolingRodRadius;
    G4double halfLengthZ = (fOuterCanExtendsPastCrystal - fInnerCanExtendsPastCrystal)/2.0;

    G4Tubs* coolingRod = new G4Tubs("coolingRod", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4double zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fInnerCanExtendsPastCrystal + halfLengthZ;

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fStructureMatCoolingRodLog == nullptr) {
        fStructureMatCoolingRodLog = new G4LogicalVolume(coolingRod, material, "structureMatCoolingRodLog", 0, 0, 0);
        fStructureMatCoolingRodLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fStructureMatCoolingRodLog, move, rotate);

    return 1;
}//end ::end AddStructureMatCoolingRod


G4int DetectionSystem8pi::AddElectrodeMatCoolingRod()
{
    //material
    G4Material* material = G4Material::GetMaterial(fElectrodeMat);
    if(!material) {
        G4cout<<" ----> Material "<<fElectrodeMat<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.75, 0.45, 0.2));
    visAtt->SetVisibility(true);

    //upper cooling rod measurements
    G4double innerRadius = 0.0*mm;
    G4double outerRadius = fElectrodeMatCoolingRodRadius;
    G4double halfLengthZ = fElectrodeMatCoolingRodLength/2.0;

    G4Tubs* coolingRod = new G4Tubs("coolingRod", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4double zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal + halfLengthZ;

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fElectrodeMatCoolingRodLog == nullptr) {
        fElectrodeMatCoolingRodLog = new G4LogicalVolume(coolingRod, material, "electrodeMatCoolingRodLog", 0, 0, 0);
        fElectrodeMatCoolingRodLog->SetVisAttributes(visAtt);
    }

    //still need to add 'ring' section where structureMat rod penetrates electrodeMat rod

    fAssembly->AddPlacedVolume(fElectrodeMatCoolingRodLog, move, rotate);

    return 1;
}//end ::end AddElectrodeMatCoolingRod

//Add the beryllium window in front of the crystal
G4int DetectionSystem8pi::AddBerylliumWindow() {
    //material
    G4Material* material = G4Material::GetMaterial("Beryllium");
    if(!material) {
        G4cout<<" ----> Material "<<"Beryllium"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5, 0.0, 0.0));
    visAtt->SetVisibility(true);

    // beryllium window measurements
    G4double innerRadius = 0.0*cm;
    G4double outerRadius = fOuterCanInnerRadius; /*crystalOuterRadius + fInnerCanThickness + fOuterCanDistFromInnerCan;*/
    G4double halfLengthZ = fBerylliumThickness/2.0;

    G4Tubs* berylliumWindow = new G4Tubs("berylliumWindow", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin -(fCrystalLength/2.0 + fBerylliumDistFromCrystal + halfLengthZ);
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fBerylliumWindowLog == nullptr) {
        fBerylliumWindowLog = new G4LogicalVolume(berylliumWindow, material, "berylliumWindowLog", 0, 0, 0);
        fBerylliumWindowLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fBerylliumWindowLog, move, rotate);

    return 1;
}//end ::end AddBerylliumWindow

//Add the outer structureMat can around the germanium crystal
G4int DetectionSystem8pi::AddOuterStructureMatCan() {
    //material
    G4Material* material = G4Material::GetMaterial(fStructureMat);
    if(!material) {
        G4cout<<" ----> Material "<<fStructureMat<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    visAtt->SetVisibility(true);

    // measurements
    G4double innerRadius = fInnerBGOOuterRadius + fInnerBGOClearance - fOuterCanThickness;
    G4double outerRadius = innerRadius + fOuterCanThickness;
    G4double halfLengthZ = fOuterCanLength/2.0;

    G4Tubs* outerCanSide = new G4Tubs("outerCanSide", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin + halfLengthZ - fCrystalLength/2.0 - fBerylliumDistFromCrystal - fBerylliumThickness;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fOuterCanSideLog == nullptr) {
        fOuterCanSideLog = new G4LogicalVolume(outerCanSide, material, "outerCanSideLog", 0, 0, 0);
        fOuterCanSideLog->SetVisAttributes(visAtt);
    }

    //add physical outerCanSide
    fAssembly->AddPlacedVolume(fOuterCanSideLog, move, rotate);

    //add outer can lid
    innerRadius = fInnerBGOInnerRadius - fInnerBGOClearance - fOuterCanThickness;
    outerRadius = fOuterCanInnerRadius;
    halfLengthZ = fOuterCanThickness/2.0;

    G4Tubs* outerCanLid = new G4Tubs("outerCanLid", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    //logical volume
    if(fOuterCanLidLog == nullptr) {
        fOuterCanLidLog = new G4LogicalVolume(outerCanLid, material, "outerCanLidLog", 0, 0, 0);
        fOuterCanLidLog->SetVisAttributes(visAtt);
    }

    //place outerCanLid
    zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal - fOuterCanThickness + halfLengthZ;
    move = zPosition * direction;

    fAssembly->AddPlacedVolume(fOuterCanLidLog, move, rotate);

    return 1;
}//end ::end AddOuterStructureMatCan

G4int DetectionSystem8pi::AddCoolingRodCover()
{
    //material
    G4Material* material = G4Material::GetMaterial(fStructureMat);
    if(!material) {
        G4cout<<" ----> Material "<<"StructureMat"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    visAtt->SetVisibility(true);

    //cooling rod cover measurements
    G4double innerRadius = fInnerBGOInnerRadius - fInnerBGOClearance - fOuterCanThickness;
    G4double outerRadius = innerRadius + fOuterCanThickness;
    G4double halfLengthZ = fElectrodeMatCoolingRodLength/2.0;

    G4Tubs* coolingRodCover = new G4Tubs("coolingRodCover", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal + halfLengthZ;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fCoolingRodCoverLog == nullptr)
    {
        fCoolingRodCoverLog = new G4LogicalVolume(coolingRodCover, material, "coolingRodCoverLog", 0, 0, 0);
        fCoolingRodCoverLog->SetVisAttributes(visAtt);
    }

    //place physical coolingRodCover
    fAssembly->AddPlacedVolume(fCoolingRodCoverLog, move, rotate);

    //cooling rod lid measurements
    innerRadius = innerRadius + fOuterCanThickness;
    outerRadius = fInnerBGOOuterRadius + fInnerBGOClearance;
    /*crystalOuterRadius + fInnerCanThickness + outerCanDistFromInnerCan;*/

    halfLengthZ = fOuterCanThickness/2.0;

    G4Tubs* coolingRodCoverLid = new G4Tubs("coolingRodLid", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    //logical volume
    if(fCoolingRodCoverLidLog == nullptr)
    {
        fCoolingRodCoverLidLog = new G4LogicalVolume(coolingRodCoverLid, material, "coolingRodCoverLidLog", 0, 0, 0);
        fCoolingRodCoverLidLog->SetVisAttributes(visAtt);
    }

    //place upper physical coolingRodCoverLid
    zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal + fElectrodeMatCoolingRodLength - halfLengthZ;
    move = zPosition * direction;

    fAssembly->AddPlacedVolume(fCoolingRodCoverLidLog, move, rotate);

    return 1;
}//end ::end AddCoolingRodCover


//Add the inner BGO annulus around the cooling rod
G4int DetectionSystem8pi::AddInnerBGOAnnulus() {
    //material
    G4Material* material = G4Material::GetMaterial("BGO");
    if(!material) {
        G4cout<<" ----> Material "<<"BGO"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.5));
    visAtt->SetVisibility(true);

    // measurements
    G4double innerRadius = fInnerBGOInnerRadius;
    G4double outerRadius = fInnerBGOOuterRadius;

    G4double halfLengthZ = fInnerBGOAnnulusLength/2.0;

    G4Tubs* innerBGOAnnulus = new G4Tubs("innerBGOAnnulus", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal + fInnerBGOClearance + halfLengthZ;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fInnerBGOAnnulusLog == nullptr) {
        fInnerBGOAnnulusLog = new G4LogicalVolume(innerBGOAnnulus, material, "8piInnerBGOAnnulus", 0, 0, 0);
        fInnerBGOAnnulusLog->SetVisAttributes(visAtt);
    }

    fAssemblyInnerBGO->AddPlacedVolume(fInnerBGOAnnulusLog, move, rotate);

    return 1;
}//end ::end AddInnerBGOAnnulus

//Add the outer structureMat sheath around the inner BGO annulus
G4int DetectionSystem8pi::AddStructureMatBGOSheath() {
    //material
    G4Material* material = G4Material::GetMaterial(fStructureMat);
    if(!material) {
        G4cout<<" ----> Material "<<fStructureMat<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    visAtt->SetVisibility(true);

    // measurements
    G4double innerRadius = fInnerBGOOuterRadius + fInnerBGOClearance;
    G4double outerRadius = innerRadius + fOuterCanThickness; //structureMatSheathThickness?
    G4double halfLengthZ = fElectrodeMatCoolingRodLength/2.0; //- fOuterCanThickness;

    G4Tubs* structureMatSheath = new G4Tubs("structureMatSheath", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal + halfLengthZ;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fStructureMatSheathLog == nullptr) {
        fStructureMatSheathLog = new G4LogicalVolume(structureMatSheath, material, "structureMatSheathLog", 0, 0, 0);
        fStructureMatSheathLog->SetVisAttributes(visAtt);
    }

    //add physical outerCanSide
    fAssembly->AddPlacedVolume(fStructureMatSheathLog, move, rotate);

    return 1;
}//end ::end AddStructureMatBGOSheath

//Add the outer BGO annulus around the detector
G4int DetectionSystem8pi::AddOuterBGOAnnulus() {
    //material
    G4Material* material = G4Material::GetMaterial("BGO");
    if(!material) {
        G4cout<<" ----> Material "<<"BGO"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.5));
    visAtt->SetVisibility(true);

    //BOTTOM (cone) portion of outer BGO
    // measurements
    G4double innerRadius = fInnerBGOOuterRadius + fInnerBGOClearance + fOuterCanThickness + fOuterBGOClearance;
    G4double lowerOuterRadius = innerRadius + fOuterBGOBottomThickness;
    G4double upperOuterRadius = fOuterBGOTopOuterRadius;
    G4double halfLengthZ = fOuterBGOTaperHeight/2.0;

    G4Cons* outerLowerBGOAnnulus = new G4Cons("outerLowerBGOAnnulus", innerRadius, lowerOuterRadius, innerRadius, upperOuterRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin + halfLengthZ - fCrystalLength/2.0 - fBerylliumDistFromCrystal - fBerylliumThickness - fOuterBGODisplacement;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fOuterLowerBGOAnnulusLog == nullptr) {
		  fOuterLowerBGOAnnulusLog = new G4LogicalVolume(outerLowerBGOAnnulus, material, "8piOuterLowerBGOAnnulus", 0, 0, 0);
        fOuterLowerBGOAnnulusLog->SetVisAttributes(visAtt);
    }

    fAssemblyOuterLowerBGO->AddPlacedVolume(fOuterLowerBGOAnnulusLog, move, rotate);

    //UPPER (annulus) portion of outer BGO
    halfLengthZ = (fOuterBGOTotalLength - fOuterBGOTaperHeight)/2.0;

    G4Tubs* outerUpperBGOAnnulus = new G4Tubs("outerUpperBGOAnnulus", innerRadius, upperOuterRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    zPosition = fCrystalDistFromOrigin + halfLengthZ + fOuterBGOTaperHeight - fCrystalLength/2.0 - fBerylliumDistFromCrystal - fBerylliumThickness - fOuterBGODisplacement;
    move = zPosition * direction;

    //logical volume
    if(fOuterUpperBGOAnnulusLog == nullptr) {
        fOuterUpperBGOAnnulusLog = new G4LogicalVolume(outerUpperBGOAnnulus, material, "8piOuterUpperBGOAnnulus", 0, 0, 0);
        fOuterUpperBGOAnnulusLog->SetVisAttributes(visAtt);
    }

    fAssemblyOuterUpperBGO->AddPlacedVolume(fOuterUpperBGOAnnulusLog, move, rotate);

    return 1;
}//end ::end AddOuterBGOAnnulus

//Add the liquid N2 container at top of detector
G4int DetectionSystem8pi::AddLiquidN2Container() {
    //material
    G4Material* material = G4Material::GetMaterial("LiquidN2");
    if(!material) {
        G4cout<<" ----> Material "<<"LiquidN2"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* N2VisAtt = new G4VisAttributes(G4Colour(0.0, 0.5, 0.5));
    N2VisAtt->SetVisibility(true);
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    visAtt->SetVisibility(true);

    // measurements
    G4double innerRadius = 0;
    G4double outerRadius = fLiquidN2Radius - fOuterCanThickness;
    G4double halfLengthZ = (fLiquidN2Length - fOuterCanThickness*2.0)/2.0;

    G4Tubs* liquidN2 = new G4Tubs("liquidN2", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal + fElectrodeMatCoolingRodLength + halfLengthZ;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fLiquidN2Log == nullptr) {
        fLiquidN2Log = new G4LogicalVolume(liquidN2, material, "liquidN2Log", 0, 0, 0);
        fLiquidN2Log->SetVisAttributes(N2VisAtt);
    }

    //add physical liquid N2
    fAssembly->AddPlacedVolume(fLiquidN2Log, move, rotate);

    //material
    material = G4Material::GetMaterial(fStructureMat);
    if(!material) {
        G4cout<<" ----> Material "<<fStructureMat<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    //place structureMat container surrounding N2
    //bottom section
    innerRadius = fInnerBGOOuterRadius + fInnerBGOClearance + fOuterCanThickness;
    outerRadius = fLiquidN2Radius - fOuterCanThickness;
    halfLengthZ = fOuterCanThickness/2.0;

    G4Tubs* liquidN2Bottom = new G4Tubs("liquidN2Bottom", innerRadius, outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle);


    zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal + fElectrodeMatCoolingRodLength + halfLengthZ - fOuterCanThickness;

    move = zPosition * direction;

    //logical volume
    if(fLiquidN2BottomLog == nullptr) {
        fLiquidN2BottomLog = new G4LogicalVolume(liquidN2Bottom, material, "liquidN2BottomLog", 0, 0, 0);
        fLiquidN2BottomLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fLiquidN2BottomLog, move, rotate);

    //side section
    innerRadius = fLiquidN2Radius - fOuterCanThickness;
    outerRadius = fLiquidN2Radius;
    halfLengthZ = fLiquidN2Length/2.0;

    G4Tubs* liquidN2Side = new G4Tubs("liquidN2Side", innerRadius,
                                        outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle
											);


    zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal + fElectrodeMatCoolingRodLength + halfLengthZ - fOuterCanThickness;
    move = zPosition * direction;

    //logical volume
    if(fLiquidN2SideLog == nullptr) {
        fLiquidN2SideLog = new G4LogicalVolume(liquidN2Side, material, "liquidN2SideLog", 0, 0, 0);
        fLiquidN2SideLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fLiquidN2SideLog, move, rotate);

    //top (lid) section
    innerRadius = 0;
    outerRadius = fLiquidN2Radius - fOuterCanThickness;
    halfLengthZ = fOuterCanThickness/2.0;

    G4Tubs* liquidN2Lid = new G4Tubs("liquidN2Lid", innerRadius,
                                       outerRadius, halfLengthZ, 0.0*deg, fDetailViewEndAngle
										);


    zPosition = fCrystalDistFromOrigin + fCrystalLength/2.0 + fOuterCanExtendsPastCrystal + fElectrodeMatCoolingRodLength + fLiquidN2Length + halfLengthZ - fOuterCanThickness*2.0;

    move = zPosition * direction;

    //logical volume
    if(fLiquidN2LidLog == nullptr) {
        fLiquidN2LidLog = new G4LogicalVolume(liquidN2Lid, material, "liquidN2LidLog", 0, 0, 0);
        fLiquidN2LidLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fLiquidN2LidLog, move, rotate);

    return 1;
}//end ::end AddLiquidN2Container


G4int DetectionSystem8pi::AddHevimetalCollimator() {
    //material
    G4Material* material = G4Material::GetMaterial("Hevimetal");
    if(!material) {
        G4cout<<" ----> Material "<<"Hevimetal"<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* hevimetalVisAtt = new G4VisAttributes(G4Colour(0.25, 0.25, 0.25));
    hevimetalVisAtt->SetVisibility(true);

    // hevimetal measurements
    G4double frontTangentDist = fHevimetalFrontSideLength * tan(60.0*deg) / 2.0; //tangent distance to inner surface
    G4double rearTangentDist =  fHevimetalRearSideLength * tan(60.0*deg) / 2.0;
    G4double halfLengthZ = fHevimetalThickness/2.0;

    G4double phiStart = 0.0*deg;		//start angle
    //G4double phiEnd = 360.0*deg;		//end angle
    G4double phiEnd = fDetailViewEndAngle;
    G4int		 numSide = 6;						//number of sides => 6: normal hexagon
    G4int		 numZPlanes = 2; 				//number of Z planes
    G4double zPlane[] = {-halfLengthZ, halfLengthZ}; //distance between Z planes
    G4double rInner[] = {0.0, 0.0}; //solid inside
    G4double rOuter[] = {frontTangentDist, rearTangentDist}; //front face and rear face tangent distances

    G4Polyhedra* hevimetal = new G4Polyhedra("hevimetal", phiStart, phiEnd, numSide, numZPlanes, zPlane, rInner, rOuter);

    // Cut Out AuxMat Plug! ///////////////////////////////
    //Place a conical slice of auxMat in the hevimetal
    G4double innerRadiusNeg = 0.0*cm;
    G4double outerRadiusNeg = fAuxMatPlugNegRadius + fCutClearance; // 0.1*mm cut clearance
    G4double innerRadiusPos = 0.0*cm;
    G4double outerRadiusPos = fAuxMatPlugPosRadius + fCutClearance; // 0.1*mm cut clearance

    G4double newHalfLengthZ = (fHevimetalThickness + fCutClearance)/2.0; // 0.1*mm cut clearance

    G4Cons* auxMatPlug = new G4Cons("auxMatPlug", innerRadiusNeg, outerRadiusNeg, innerRadiusPos, outerRadiusPos, newHalfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector moveCut(0,0,0);

    G4SubtractionSolid* hevimetalWithCavity = new G4SubtractionSolid("hevimetalWithCavity", hevimetal, auxMatPlug, 0, moveCut);
    ////////////////////////////////////////////////////////

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin - fCrystalLength/2.0 - fBerylliumDistFromCrystal - fBerylliumThickness - fOuterBGODisplacement - halfLengthZ;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fHevimetalLog == nullptr) {
        fHevimetalLog = new G4LogicalVolume(hevimetalWithCavity, material, "hevimetalLog", 0, 0, 0);
        fHevimetalLog->SetVisAttributes(hevimetalVisAtt);
    }

    fAssembly->AddPlacedVolume(fHevimetalLog, move, rotate);

    return 1;
}//end ::end AddHevimetalCollimator


G4int DetectionSystem8pi::AddAuxMatPlug() {
    //material
    G4Material* material = G4Material::GetMaterial(fAuxMat);
    if(!material) {
        G4cout<<" ----> Material "<<fAuxMat<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
    visAtt->SetVisibility(true);

    //Place a conical slice of auxMat in the hevimetal
    G4double innerRadiusNeg = 0.0*cm;
    G4double outerRadiusNeg = fAuxMatPlugNegRadius;
    G4double innerRadiusPos = 0.0*cm;
    G4double outerRadiusPos = fAuxMatPlugPosRadius;

    G4double halfLengthZ = fHevimetalThickness/2.0;

    G4Cons* auxMatPlug = new G4Cons("auxMatPlug", innerRadiusNeg, outerRadiusNeg, innerRadiusPos, outerRadiusPos, halfLengthZ, 0.0*deg, fDetailViewEndAngle);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin - fCrystalLength/2.0 - fBerylliumDistFromCrystal - fBerylliumThickness - fOuterBGODisplacement - halfLengthZ;

    //Nest the auxMat plug in the hevimetal
    G4ThreeVector move = zPosition * direction ; // check this!

    //logical volume
    if(fAuxMatPlugLog == nullptr) {
        fAuxMatPlugLog = new G4LogicalVolume(auxMatPlug, material, "auxMatPlugLog", 0, 0, 0);
        fAuxMatPlugLog->SetVisAttributes(visAtt);
    }

    G4ThreeVector moveNull = G4ThreeVector(0, 0, 0);

    fAssembly->AddPlacedVolume(fAuxMatPlugLog, move, rotate);

    return 1;
}//end ::end AddAuxMatPlug

G4int DetectionSystem8pi::AddThinAuxMatLayer() {
    //material
    G4Material* material = G4Material::GetMaterial(fAuxMat);
    if(!material) {
        G4cout<<" ----> Material "<<fAuxMat<<" not found, cannot build the detector shell! "<<G4endl;
        return 0;
    }

    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
    visAtt->SetVisibility(true);

    // auxMat layer measurements
    G4double frontTangentDist = fAuxMatLayerFrontSideLength * tan(60.0*deg) / 2.0;//front tangent distance to edge
    G4double rearTangentDist =  fHevimetalFrontSideLength * tan(60.0*deg) / 2.0; //rear tangent distance to edge
    G4double halfLengthZ = fAuxMatLayerThickness/2.0;

    G4double phiStart = 0.0*deg;		//start angle
    //G4double phiEnd = 360.0*deg;		//end angle
    G4double phiEnd = fDetailViewEndAngle;
    G4int		 numSide = 6;						//number of sides => 6: normal hexagon
    G4int		 numZPlanes = 2; 				//number of Z planes
    G4double zPlane[] = {-halfLengthZ, halfLengthZ}; //distance between Z planes
    G4double rInner[] = {0.0, 0.0}; //solid inside
    G4double rOuter[] = {frontTangentDist, rearTangentDist}; //front face and rear face tangent distances


    G4Polyhedra* auxMatLayer = new G4Polyhedra("auxMatLayer", phiStart, phiEnd, numSide, numZPlanes, zPlane, rInner, rOuter);

    G4ThreeVector direction = G4ThreeVector(0,0,1);
    G4RotationMatrix* rotate = new G4RotationMatrix;
    G4double zPosition = fCrystalDistFromOrigin - fCrystalLength/2.0 - fBerylliumDistFromCrystal - fBerylliumThickness - fOuterBGODisplacement - fHevimetalThickness - halfLengthZ;
    G4ThreeVector move = zPosition * direction;

    //logical volume
    if(fAuxMatLayerLog == nullptr) {
        fAuxMatLayerLog = new G4LogicalVolume(auxMatLayer, material, "auxMatLayerLog", 0, 0, 0);
        fAuxMatLayerLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fAuxMatLayerLog, move, rotate);

    return 1;
}//end ::end AddThinAuxMatLayer
