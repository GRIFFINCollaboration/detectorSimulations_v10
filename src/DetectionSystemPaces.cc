#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemPaces.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


DetectionSystemPaces::DetectionSystemPaces() :
	// Logical Volumes
	fAluminumHemisphereLog(0),
	fAluminumAnnulusTopLog(0),
	fAluminumAnnulusBotLog(0),
	fCanisterLog(0),
	fSiliconBlockLog(0),
	fSiliconDeadLayerLog(0),

	fTeflonAnnulusTopLog(0),
	fTeflonAnnulusBotLog(0),
	fDelrinHemisphereLog(0)
{
	/*Measurement status*/
	//~   estimate
	//OK  confirmed measurement

	// Cut clearance
	fCutClearance = 0.01*mm;

	// Aluminum hemisphere
	fAluminumHemisphereInnerRadius =         (82.50 / 2)*mm; //OK
	fAluminumHemisphereOuterRadius =         (88.86 / 2)*mm; //OK
	fAluminumHemisphereBeamHoleRadius =     (34.50 / 2)*mm; //~
	fAluminumHemisphereBeamHoleRimHeight = (3.0)*mm; //~

	// Annulus that sits in front of silicon crystal
	fAluminumAnnulusTopInnerRadius =        (15.86 / 2)*mm; //OK
	fAluminumAnnulusTopOuterRadius =        (31.71 / 2)*mm; //OK
	fAluminumAnnulusTopThickness =           (2.37)*mm; //OK

	// Annulus that sits in back of silicon crystal
	fAluminumAnnulusBotInnerRadius =        (15.86 / 2)*mm; //~
	fAluminumAnnulusBotOuterRadius =        (31.71 / 2)*mm; //~
	fAluminumAnnulusBotThickness =           (4.50)*mm; //~

	// Canister holding crystal
	fCanisterInnerRadius =                    (18.50 / 2)*mm; //~
	fCanisterOuterRadius =                    (21.00 / 2)*mm; //~
	fCanisterThickness =                       (5.00)*mm; //~

	// Silicon crystal detector
	fSiliconBlockRadius =                     (15.86 / 2)*mm; //~
	fSiliconBlockThickness =                  (4.90)*mm; //~
	fSiliconDeadLayerThickness =             (0.004)*mm; //~ Nominal Structure Stopping Power of Window - Ortec, 2u Si equivalent

	// Front and back Au (0.2 and 2um) contacts.
	//frontContactThickness =  		   (0.0002)*mm; //~ Common front contact thickness for Si(Li) Zhimin
	//backContactThickness =  		   (0.002)*mm; //~ Common back contact thickness for Si(Li) Zhimin

	// Screws
	fScrewRadius =                             (1.00 / 2)*mm; //~
	fScrewPlacementRadius =                   (22.50 / 2)*mm; //~

	// Teflon bases
	fTeflonAnnulusTopInnerRadius =          (21.00 / 2)*mm; //~
	fTeflonAnnulusTopOuterRadius =          (32.00 / 2)*mm; //~
	fTeflonAnnulusTopThickness =             (2.22)*mm; //OK
	fTeflonAnnulusBotInnerRadius =          (31.65 / 2)*mm; //~
	fTeflonAnnulusBotOuterRadius =          (34.50 / 2)*mm; //~
	fTeflonAnnulusBotThickness =             (4.00)*mm; //~

	// delrin hemisphere
	fDelrinHemisphereInnerRadius =           (84.4)*mm; //OK
	fDelrinHemisphereOuterRadius =           (89.4)*mm; //OK
	fDelrinHemisphereBeamHoleRadius =       (20.00 / 2)*mm; //OK

	// placement distances
	fAluminumHemisphereDist =     (0.)*mm;
	fDelrinHemisphereDist =       (0.)*mm;
	fAluminumAnnulusTopDist =    (fAluminumAnnulusTopThickness / 2);
	fCanisterDist =              (fAluminumAnnulusTopThickness) + (fCanisterThickness / 2) ;
	fSiliconBlockDist =          (fAluminumAnnulusTopThickness) + (fSiliconBlockThickness / 2); // Takes dead layer into account
	fSiliconDeadLayerFrontDist = (fAluminumAnnulusTopThickness) + (fSiliconDeadLayerThickness / 2) ;
	fSiliconDeadLayerBackDist  = (fAluminumAnnulusTopThickness) + (fSiliconBlockThickness) - (fSiliconDeadLayerThickness / 2);
	fTeflonAnnulusTopDist =      (fAluminumAnnulusTopThickness) + (fSiliconBlockThickness) - (fTeflonAnnulusTopThickness / 2);
	fTeflonAnnulusBotDist =      (fTeflonAnnulusTopDist) + (fTeflonAnnulusTopThickness / 2) + (fTeflonAnnulusBotThickness / 2);
	fAluminumAnnulusBotDist =    (fTeflonAnnulusTopDist) + (fTeflonAnnulusTopThickness / 2) + (fAluminumAnnulusBotThickness / 2);

	// Placement and Orientation scalars
	/*pacesPlacementDistance[5] = {  30.00*mm ,  30.00*mm ,  30.00*mm ,  30.00*mm ,  30.00*mm };
	  pacesPlacementPhi[5] =    {   0.0*deg ,  72.0*deg , 144.0*deg , 216.0*deg , 288.0*deg };
	  pacesPlacementTheta[5] =      {  60.0*deg ,  60.0*deg ,  60.0*deg ,  60.0*deg ,  60.0*deg };
	  pacesOrientationPhi[5] =  {   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg };
	  pacesOrientationTheta[5] =    {   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg };*/

	// distance is deduced from Alumium hemisphere assembly
	fPacesPlacementDistance[0] = 29.50*mm; // 29.5*mm Deduced from Alumium hemisphere
	fPacesPlacementDistance[1] = 29.50*mm;
	fPacesPlacementDistance[2] = 29.50*mm;
	fPacesPlacementDistance[3] = 29.50*mm;
	fPacesPlacementDistance[4] = 29.50*mm;

	// phi angle from measurement
	fPacesPlacementPhi[0] =   0.0000*deg; //OK
	fPacesPlacementPhi[1] =  -68.7575*deg; //OK
	fPacesPlacementPhi[2] = -138.4082*deg; //OK
	fPacesPlacementPhi[3] = -212.1280*deg; //OK
	fPacesPlacementPhi[4] = -287.4275*deg; //OK

	// theta angle is guess estimate
	fPacesPlacementTheta[0] = (180 - 60.0)*deg; //~
	fPacesPlacementTheta[1] = (180 - 60.0)*deg; //~
	fPacesPlacementTheta[2] = (180 - 60.0)*deg; //~
	fPacesPlacementTheta[3] = (180 - 60.0)*deg; //~
	fPacesPlacementTheta[4] = (180 - 60.0)*deg; //~

	fPacesOrientationPhi[0] = 0.0*deg;
	fPacesOrientationPhi[1] = 0.0*deg;
	fPacesOrientationPhi[2] = 0.0*deg;
	fPacesOrientationPhi[3] = 0.0*deg;
	fPacesOrientationPhi[4] = 0.0*deg;

	fPacesOrientationTheta[0] = 0.1780*deg; //OK
	fPacesOrientationTheta[1] = 0.8275*deg; //OK
	fPacesOrientationTheta[2] =-0.2579*deg; //OK
	fPacesOrientationTheta[3] = 0.2992*deg; //OK
	fPacesOrientationTheta[4] = 0.1934*deg; //OK

	// assembly volume sizes
	fDetectorAssemblyRadius =    (34.50 / 2)*mm; //?
	fDetectorAssemblyThickness = (15.0)*mm; //?
}

DetectionSystemPaces::~DetectionSystemPaces() {
	// logical volumes
	delete fAluminumHemisphereLog;
	delete fAluminumAnnulusTopLog;
	delete fAluminumAnnulusBotLog;
	delete fCanisterLog;
	delete fSiliconBlockLog;
	delete fSiliconDeadLayerLog;
	delete fTeflonAnnulusTopLog;
	delete fTeflonAnnulusBotLog;
	delete fDelrinHemisphereLog;

}

G4int DetectionSystemPaces::Build() {

	// Build assembly volumes
	fAssembly = new G4AssemblyVolume();
	fAssemblyDetector = new G4AssemblyVolume();
	fAssemblySilicon = new G4AssemblyVolume();

	// Add silicon assembly
	AddSiliconBlock();
	AddSiliconDeadLayer();

	// Add detector assembly
	AddCanister();
	AddAluminumAnnulusTop();
	AddTeflonAnnulusTop();
	AddTeflonAnnulusBot();
	AddAluminumAnnulusBot();
	AddScrews();

	// Add total assembly
	AddAluminumHemisphere();

	// Assemble the ensemble
	CombineAssemblySilicon();
	CombineAssemblyDetector();
	CombineAssembly();

	return 1;
}//end ::Build

G4int DetectionSystemPaces::PlaceDetector(G4LogicalVolume* expHallLog, G4int ndet) {
	//place detectors
	G4int dI;
	G4double* ptrPd = fPacesPlacementDistance;
	G4double* ptrPt = fPacesPlacementPhi;//t=phi, but p=theta???????
	G4double* ptrPp = fPacesPlacementTheta;
	G4double* ptrOt = fPacesOrientationPhi;
	G4double* ptrOp = fPacesOrientationTheta;
	/*G4double* ptrPp = fPacesPlacementPhi;//changed code, rotated by 90 degrees
	  G4double* ptrPt = fPacesPlacementTheta;
	  G4double* ptrOp = fPacesOrientationPhi;
	  G4double* ptrOt = fPacesOrientationTheta;*/
	if(ndet > 5 || ndet < 0) ndet = 5;
	G4double dDist, dPhi, dTheta, oriPhi, oriTheta;
	G4RotationMatrix* Ra;
	G4ThreeVector Ta, yprimeaxis;
	for(G4int i=0; i<ndet; i++) {
		dI = i;
		dDist = ptrPd[dI];
		dPhi = ptrPt[dI];
		dTheta = ptrPp[dI];
		Ra = new G4RotationMatrix;
		Ta.setX(dDist * cos(dPhi) * sin(dTheta));
		Ta.setY(dDist * sin(dPhi) * sin(dTheta));
		Ta.setZ(dDist *      1.0     * cos(dTheta));

		oriPhi = dPhi + ptrOt[dI] + M_PI/2; //plus 90 deg
		oriTheta = dTheta + ptrOp[dI];
		yprimeaxis = G4ThreeVector(cos(oriPhi), sin(oriPhi), 0);
		Ra->set(yprimeaxis, oriTheta);

		//G4cout<<"----------- dI = "<<dI<<G4endl;
		fAssemblySilicon->MakeImprint(expHallLog, Ta, Ra, dI+1);
		fAssemblyDetector->MakeImprint(expHallLog, Ta, Ra, dI*20);
	}

	//place hemisphere
	G4RotationMatrix* R0 = new G4RotationMatrix; G4ThreeVector T0;
	fAssembly->MakeImprint(expHallLog, T0, R0, 5*40);

	return 1;
}//end ::PlaceDetector

G4int DetectionSystemPaces::CombineAssemblySilicon() {
	G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
	Ta.setZ(fSiliconBlockDist);
	fAssemblySilicon->AddPlacedVolume(fSiliconBlockLog, Ta, Ra);

	return 1;
}//end ::CombineAssemblySilicon


G4int DetectionSystemPaces::CombineAssemblyDetector() {
	G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
	Ta.setZ(fSiliconDeadLayerFrontDist);
	fAssemblyDetector->AddPlacedVolume(fSiliconDeadLayerLog, Ta, Ra);
	Ta.setZ(fSiliconDeadLayerBackDist);
	fAssemblyDetector->AddPlacedVolume(fSiliconDeadLayerLog, Ta, Ra);
	Ta.setZ(fCanisterDist);
	fAssemblyDetector->AddPlacedVolume(fCanisterLog, Ta, Ra);
	Ta.setZ(fAluminumAnnulusTopDist);
	fAssemblyDetector->AddPlacedVolume(fAluminumAnnulusTopLog, Ta, Ra);
	Ta.setZ(fAluminumAnnulusBotDist);
	fAssemblyDetector->AddPlacedVolume(fAluminumAnnulusBotLog, Ta, Ra);
	Ta.setZ(fTeflonAnnulusTopDist);
	fAssemblyDetector->AddPlacedVolume(fTeflonAnnulusTopLog, Ta, Ra);
	Ta.setZ(fTeflonAnnulusBotDist);
	fAssemblyDetector->AddPlacedVolume(fTeflonAnnulusBotLog, Ta, Ra);

	return 1;
}//end ::CombineAssemblyDetector

G4int DetectionSystemPaces::CombineAssembly() {
	//place aluminum hemisphere
	G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
	Ta.setZ(fAluminumHemisphereDist);
	fAssembly->AddPlacedVolume(fAluminumHemisphereLog, Ta, Ra);

	return 1;
}//end ::CombineAssembly

//Add silicon detector
G4int DetectionSystemPaces::AddSiliconBlock() {
	//material
	G4Material* material = G4Material::GetMaterial("Silicon");
	if(!material) {
		G4cout<<" ----> Material "<<"Silicon"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* siliconBlockVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	siliconBlockVisAtt->SetVisibility(true);

	//silicon block (without dead layer, smaller than actual)
	G4double deadLayerCut = fSiliconDeadLayerThickness + fCutClearance;
	G4double innerRadius = 0.*mm;
	G4double outerRadius = fSiliconBlockRadius;
	G4double halfLengthZ = fSiliconBlockThickness/2.0 - deadLayerCut;

	//primitive volume
	G4Tubs* siliconBlock = new G4Tubs("siliconBlock", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);

	//logical volume
	if(fSiliconBlockLog == nullptr) {
		fSiliconBlockLog = new G4LogicalVolume(siliconBlock, material, "pacesSiliconBlockLog", 0, 0, 0); // Renamed from "siliconBlockLog" to "pacesSiliconBlockLog"
		fSiliconBlockLog->SetVisAttributes(siliconBlockVisAtt);
	}

	return 1;
}//end ::end AddSiliconBlock

//Add detector dead layer
G4int DetectionSystemPaces::AddSiliconDeadLayer() {
	//material
	G4Material* material = G4Material::GetMaterial("Silicon"); // "Nominal Structure Stopping Power of Window-Ortec, 2u Si"
	if(!material) {
		G4cout<<" ----> Material "<<"Silicon"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* siliconDeadLayerVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
	siliconDeadLayerVisAtt->SetVisibility(true);

	//dead layer
	G4double innerRadius = 0.*mm;
	G4double outerRadius = fSiliconBlockRadius;
	G4double halfLengthZ = fSiliconDeadLayerThickness / 2.0;

	//primitive volume
	G4Tubs* siliconDeadLayer = new G4Tubs("siliconDeadLayer", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);

	//logical volume
	if(fSiliconDeadLayerLog == nullptr) {
		fSiliconDeadLayerLog = new G4LogicalVolume(siliconDeadLayer, material, "siliconDeadLayerLog", 0, 0, 0);
		fSiliconDeadLayerLog->SetVisAttributes(siliconDeadLayerVisAtt);
	}

	return 1;
}//end ::end AddSiliconDeadLayer

//Add detector canister
G4int DetectionSystemPaces::AddCanister() {
	//material
	G4Material* material = G4Material::GetMaterial("Silicon");  // would be Silicon, Zhimin
	if(!material) {
		G4cout<<" ----> Material "<<"Silicon"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* canisterVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	canisterVisAtt->SetVisibility(true);

	//main annulus
	G4double innerRadius = fCanisterInnerRadius + fCutClearance;
	G4double outerRadius = fCanisterOuterRadius;
	G4double halfLengthZ = fCanisterThickness/2.0;

	//primitive volume
	G4Tubs* canister = new G4Tubs("canister", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);

	//logical volume
	if(fCanisterLog == nullptr) {
		fCanisterLog = new G4LogicalVolume(canister, material, "canisterLog", 0, 0, 0);
		fCanisterLog->SetVisAttributes(canisterVisAtt);
	}

	return 1;
}//end ::end AddCanister

//Add annulus that is on top of silicon detector
G4int DetectionSystemPaces::AddAluminumAnnulusTop() {
	//material
	G4Material* material = G4Material::GetMaterial("Aluminum");
	if(!material) {
		G4cout<<" ----> Material "<<"Aluminum"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* aluminumAnnulusTopVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.0));
	aluminumAnnulusTopVisAtt->SetVisibility(true);

	//main annulus
	G4double innerRadius = fAluminumAnnulusTopInnerRadius;
	G4double outerRadius = fAluminumAnnulusTopOuterRadius;
	G4double halfLengthZ = fAluminumAnnulusTopThickness/2.0;
	//screw cuts
	G4double cutInnerRadius = 0.*mm;
	G4double cutOuterRadius = fScrewRadius;
	G4double cutHalfLengthZ = 20.*mm;

	//primitive volume
	G4Tubs* alAnToWithoutCuts = new G4Tubs("alAnToWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);
	G4Tubs* cutScrewHole = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0*M_PI, 2*M_PI);

	//cut out screw holes
	G4double cI;
	G4double cR = fScrewPlacementRadius;
	G4double cDt = 45.*deg;
	G4ThreeVector moveCut;
	cI = 0; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnTo0 = new G4SubtractionSolid("alAnTo0", alAnToWithoutCuts, cutScrewHole, 0, moveCut);
	cI = 1; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnTo1 = new G4SubtractionSolid("alAnTo1", alAnTo0, cutScrewHole, 0, moveCut);
	cI = 2; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnTo2 = new G4SubtractionSolid("alAnTo2", alAnTo1, cutScrewHole, 0, moveCut);
	cI = 3; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnTo3 = new G4SubtractionSolid("alAnTo3", alAnTo2, cutScrewHole, 0, moveCut);
	cI = 4; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnTo4 = new G4SubtractionSolid("alAnTo4", alAnTo3, cutScrewHole, 0, moveCut);
	cI = 5; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnTo5 = new G4SubtractionSolid("alAnTo5", alAnTo4, cutScrewHole, 0, moveCut);
	cI = 6; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnTo6 = new G4SubtractionSolid("alAnTo6", alAnTo5, cutScrewHole, 0, moveCut);
	//final cut
	cI = 7; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* aluminumAnnulusTop = new G4SubtractionSolid("aluminumAnnulusTop", alAnTo6, cutScrewHole, 0, moveCut);

	//logical volume
	if(fAluminumAnnulusTopLog == nullptr) {
		fAluminumAnnulusTopLog = new G4LogicalVolume(aluminumAnnulusTop, material, "aluminumAnnulusTopLog", 0, 0, 0);
		fAluminumAnnulusTopLog->SetVisAttributes(aluminumAnnulusTopVisAtt);
	}

	return 1;
}//end ::end AddAluminumAnnulusTop

//Add annulus that is on bottom of silicon detector
G4int DetectionSystemPaces::AddAluminumAnnulusBot() {
	//material
	G4Material* material = G4Material::GetMaterial("Aluminum");
	if(!material) {
		G4cout<<" ----> Material "<<"Aluminum"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* aluminumAnnulusBotVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.0));
	aluminumAnnulusBotVisAtt->SetVisibility(true);

	//main annulus
	G4double innerRadius = fAluminumAnnulusBotInnerRadius;
	G4double outerRadius = fAluminumAnnulusBotOuterRadius;
	G4double halfLengthZ = fAluminumAnnulusBotThickness/2.0;
	//screw cuts
	G4double cutInnerRadius = 0.*mm;
	G4double cutOuterRadius = fScrewRadius;
	G4double cutHalfLengthZ = 20.*mm;

	//primitive volume
	G4Tubs* alAnBoWithoutCuts = new G4Tubs("alAnBoWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);
	G4Tubs* cutScrewHole = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0*M_PI, 2*M_PI);

	//cut out screw holes
	G4double cI;
	G4double cR = fScrewPlacementRadius;
	G4double cDt = 45.*deg;
	G4ThreeVector moveCut;
	cI = 0; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnBo0 = new G4SubtractionSolid("alAnBo0", alAnBoWithoutCuts, cutScrewHole, 0, moveCut);
	cI = 1; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnBo1 = new G4SubtractionSolid("alAnBo1", alAnBo0, cutScrewHole, 0, moveCut);
	cI = 2; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnBo2 = new G4SubtractionSolid("alAnBo2", alAnBo1, cutScrewHole, 0, moveCut);
	cI = 3; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnBo3 = new G4SubtractionSolid("alAnBo3", alAnBo2, cutScrewHole, 0, moveCut);
	cI = 4; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnBo4 = new G4SubtractionSolid("alAnBo4", alAnBo3, cutScrewHole, 0, moveCut);
	cI = 5; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnBo5 = new G4SubtractionSolid("alAnBo5", alAnBo4, cutScrewHole, 0, moveCut);
	cI = 6; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* alAnBo6 = new G4SubtractionSolid("alAnBo6", alAnBo5, cutScrewHole, 0, moveCut);
	//final cut
	cI = 7; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* aluminumAnnulusBot = new G4SubtractionSolid("aluminumAnnulusBot", alAnBo6, cutScrewHole, 0, moveCut);

	//logical volume
	if(fAluminumAnnulusBotLog == nullptr) {
		fAluminumAnnulusBotLog = new G4LogicalVolume(aluminumAnnulusBot, material, "aluminumAnnulusBotLog", 0, 0, 0);
		fAluminumAnnulusBotLog->SetVisAttributes(aluminumAnnulusBotVisAtt);
	}

	return 1;
}//end ::end AddAluminumAnnulusBot

//Add annulus that is on side of silicon detector
G4int DetectionSystemPaces::AddTeflonAnnulusTop() {
	//material
	G4Material* material = G4Material::GetMaterial("Teflon");
	if(!material) {
		G4cout<<" ----> Material "<<"Teflon"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* teflonAnnulusTopVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	teflonAnnulusTopVisAtt->SetVisibility(true);

	//main annulus
	G4double innerRadius = fTeflonAnnulusTopInnerRadius;
	G4double outerRadius = fTeflonAnnulusTopOuterRadius;
	G4double halfLengthZ = fTeflonAnnulusTopThickness/2.0;
	//screw cuts
	G4double cutInnerRadius = 0.*mm;
	G4double cutOuterRadius = fScrewRadius;
	G4double cutHalfLengthZ = 20.*mm;

	//primitive volume
	G4Tubs* teAnToWithoutCuts = new G4Tubs("teAnToWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);
	G4Tubs* cutScrewHole = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0*M_PI, 2*M_PI);

	//cut out screw holes
	G4double cI;
	G4double cR = fScrewPlacementRadius;
	G4double cDt = 45.*deg;
	G4ThreeVector moveCut;
	cI = 0; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnTo0 = new G4SubtractionSolid("teAnTo0", teAnToWithoutCuts, cutScrewHole, 0, moveCut);
	cI = 1; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnTo1 = new G4SubtractionSolid("teAnTo1", teAnTo0, cutScrewHole, 0, moveCut);
	cI = 2; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnTo2 = new G4SubtractionSolid("teAnTo2", teAnTo1, cutScrewHole, 0, moveCut);
	cI = 3; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnTo3 = new G4SubtractionSolid("teAnTo3", teAnTo2, cutScrewHole, 0, moveCut);
	cI = 4; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnTo4 = new G4SubtractionSolid("teAnTo4", teAnTo3, cutScrewHole, 0, moveCut);
	cI = 5; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnTo5 = new G4SubtractionSolid("teAnTo5", teAnTo4, cutScrewHole, 0, moveCut);
	cI = 6; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnTo6 = new G4SubtractionSolid("teAnTo6", teAnTo5, cutScrewHole, 0, moveCut);
	//final cut
	cI = 7; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teflonAnnulusTop = new G4SubtractionSolid("teflonAnnulusTop", teAnTo6, cutScrewHole, 0, moveCut);

	//logical volume
	if(fTeflonAnnulusTopLog == nullptr) {
		fTeflonAnnulusTopLog = new G4LogicalVolume(teflonAnnulusTop, material, "teflonAnnulusTopLog", 0, 0, 0);
		fTeflonAnnulusTopLog->SetVisAttributes(teflonAnnulusTopVisAtt);
	}

	return 1;
}//end ::end AddTeflonAnnulusTop

//Add annulus that is on bottom of silicon detector
G4int DetectionSystemPaces::AddTeflonAnnulusBot() {
	//material
	G4Material* material = G4Material::GetMaterial("Teflon");
	if(!material) {
		G4cout<<" ----> Material "<<"Teflon"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* teflonAnnulusBotVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	teflonAnnulusBotVisAtt->SetVisibility(true);

	//main annulus
	G4double innerRadius = fTeflonAnnulusBotInnerRadius;
	G4double outerRadius = fTeflonAnnulusBotOuterRadius;
	G4double halfLengthZ = fTeflonAnnulusBotThickness/2.0;
	//screw cuts
	G4double cutInnerRadius = 0.*mm;
	G4double cutOuterRadius = fScrewRadius;
	G4double cutHalfLengthZ = 20.*mm;

	//primitive volume
	G4Tubs* teAnBoWithoutCuts = new G4Tubs("teAnBoWithoutCuts", innerRadius, outerRadius, halfLengthZ, 0*M_PI, 2*M_PI);
	G4Tubs* cutScrewHole = new G4Tubs("cutScrewHole", cutInnerRadius, cutOuterRadius, cutHalfLengthZ, 0*M_PI, 2*M_PI);

	//cut out screw holes
	G4double cI;
	G4double cR = fScrewPlacementRadius;
	G4double cDt = 45.*deg;
	G4ThreeVector moveCut;
	cI = 0; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnBo0 = new G4SubtractionSolid("teAnBo0", teAnBoWithoutCuts, cutScrewHole, 0, moveCut);
	cI = 1; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnBo1 = new G4SubtractionSolid("teAnBo1", teAnBo0, cutScrewHole, 0, moveCut);
	cI = 2; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnBo2 = new G4SubtractionSolid("teAnBo2", teAnBo1, cutScrewHole, 0, moveCut);
	cI = 3; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnBo3 = new G4SubtractionSolid("teAnBo3", teAnBo2, cutScrewHole, 0, moveCut);
	cI = 4; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnBo4 = new G4SubtractionSolid("teAnBo4", teAnBo3, cutScrewHole, 0, moveCut);
	cI = 5; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnBo5 = new G4SubtractionSolid("teAnBo5", teAnBo4, cutScrewHole, 0, moveCut);
	cI = 6; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teAnBo6 = new G4SubtractionSolid("teAnBo6", teAnBo5, cutScrewHole, 0, moveCut);
	//final cut
	cI = 7; moveCut.setX(cR*cos(cI*cDt)); moveCut.setY(cR*sin(cI*cDt));
	G4SubtractionSolid* teflonAnnulusBot = new G4SubtractionSolid("teflonAnnulusBot", teAnBo6, cutScrewHole, 0, moveCut);

	//logical volume
	if(fTeflonAnnulusBotLog == nullptr) {
		fTeflonAnnulusBotLog = new G4LogicalVolume(teflonAnnulusBot, material, "teflonAnnulusBotLog", 0, 0, 0);
		fTeflonAnnulusBotLog->SetVisAttributes(teflonAnnulusBotVisAtt);
	}

	return 1;
}//end ::end AddTeflonAnnulusBot

//Add the aluminum hemisphere, housing for PACES
G4int DetectionSystemPaces::AddAluminumHemisphere() {
	//material
	G4Material* material = G4Material::GetMaterial("Aluminum");
	if(!material) {
		G4cout<<" ----> Material "<<"Aluminum"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	//vis attributes
	G4VisAttributes* aluminumHemisphereVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	aluminumHemisphereVisAtt->SetVisibility(true);

	//main hemisphere shell
	G4double innerRadius = fAluminumHemisphereInnerRadius;
	G4double outerRadius = fAluminumHemisphereOuterRadius;
	G4double thickness = std::abs(outerRadius - innerRadius);
	G4double startPhi = M_PI/2;
	G4double endPhi = M_PI - startPhi;
	//beam hole rim
	G4double rimInnerRadius = std::abs(innerRadius - fAluminumHemisphereBeamHoleRimHeight);
	G4double rimOuterRadius = innerRadius;
	G4double rimEndPhi = asin((fAluminumHemisphereBeamHoleRadius + thickness) / innerRadius);
	G4double rimStartPhi = M_PI - rimEndPhi;
	//beam hole cut
	G4double beamHoleInnerRadius = 0.*mm;
	G4double beamHoleOuterRadius = fAluminumHemisphereBeamHoleRadius;
	G4double beamHoleHalfLengthZ = outerRadius;
	//cuts where detectors go
	G4double detectorInnerRadius = 0.*mm;
	G4double detectorOuterRadius = fDetectorAssemblyRadius + fCutClearance;
	G4double detectorHalfLengthZ = fDetectorAssemblyThickness/2.0;

	//primitive volumes
	G4Sphere* alHeShell = new G4Sphere("alHeShell", innerRadius, outerRadius, 0*M_PI, 2*M_PI, startPhi, endPhi);
	G4Sphere* alHeRim = new G4Sphere("alHeRim", rimInnerRadius, rimOuterRadius, 0*M_PI, 2*M_PI, rimStartPhi, rimEndPhi);
	G4Tubs* cutBeamHole = new G4Tubs("cutBeamHole", beamHoleInnerRadius, beamHoleOuterRadius, beamHoleHalfLengthZ, 0*M_PI, 2*M_PI);
	G4Tubs* cutDetector = new G4Tubs("cutDetector", detectorInnerRadius, detectorOuterRadius, detectorHalfLengthZ, 0*M_PI, 2*M_PI);

	//add rim and beam hole
	G4RotationMatrix* R0 = new G4RotationMatrix; G4ThreeVector T0;
	G4UnionSolid* alHeShellRim = new G4UnionSolid("alHeShellRim", alHeShell, alHeRim, R0, T0);
	G4SubtractionSolid* alHeShellRimHole = new G4SubtractionSolid("alHeShellRimHole", alHeShellRim, cutBeamHole, R0, T0);
	//cut out detectors

	G4double* ptrPd = fPacesPlacementDistance;
	G4double* ptrPt = fPacesPlacementPhi;
	G4double* ptrPp = fPacesPlacementTheta;
	G4double* ptrOt = fPacesOrientationPhi;
	G4double* ptrOp = fPacesOrientationTheta;

	//  G4RotationMatrix* rotateCut[5];
	//  G4ThreeVector moveCut[5], yprimeaxis;
	//  G4double dDist, dPhi, dTheta, oriPhi, oriTheta;
	//
	//  for(int i=0; i<5; i++)
	//  {
	//    dDist = ptrPd[i] + detectorHalfLengthZ;
	//    dPhi = ptrPt[i];
	//    dTheta = ptrPp[i];
	//    rotateCut[i] = new G4RotationMatrix;
	//    moveCut[i].setX(dDist * cos(dPhi) * sin(dTheta));
	//    moveCut[i].setY(dDist * sin(dPhi) * sin(dTheta));
	//    moveCut[i].setZ(dDist *    1.0     * cos(dTheta));
	//    oriPhi = dPhi + ptrOt[i] - M_PI/2; //minus 90 deg
	//    oriTheta = dTheta + ptrOp[i];
	////    yprimeaxis = G4ThreeVector(cos(oriPhi), sin(oriPhi), 0);
	//	  yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	//    rotateCut[i]->set(yprimeaxis, oriTheta);
	//  }
	//  G4SubtractionSolid* alHe0 = new G4SubtractionSolid("alHe0", alHeShellRimHole, cutDetector, rotateCut[0], moveCut[0]);
	//  G4SubtractionSolid* alHe1 = new G4SubtractionSolid("alHe1", alHe0, cutDetector, rotateCut[1], moveCut[1]);
	//  G4SubtractionSolid* alHe2 = new G4SubtractionSolid("alHe2", alHe1, cutDetector, rotateCut[2], moveCut[2]);
	//  G4SubtractionSolid* alHe3 = new G4SubtractionSolid("alHe3", alHe2, cutDetector, rotateCut[3], moveCut[3]);
	//  G4SubtractionSolid* aluminumHemisphere = new G4SubtractionSolid("aluminumHemisphere", alHe3, cutDetector, rotateCut[4], moveCut[4]);
	//  for(int i=0; i<5; i++) delete rotateCut[i]; //safety

	G4ThreeVector moveCut, yprimeaxis;
	G4RotationMatrix* rotateCut = new G4RotationMatrix;
	G4double dDist, dPhi, dTheta, oriPhi, oriTheta;
	G4int dI;
	//one
	dI = 0;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* alHe0 = new G4SubtractionSolid("alHe0", alHeShellRimHole, cutDetector, rotateCut, moveCut);
	//two
	dI = 1;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* alHe1 = new G4SubtractionSolid("alHe1", alHe0, cutDetector, rotateCut, moveCut);
	//three
	dI = 2;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* alHe2 = new G4SubtractionSolid("alHe2", alHe1, cutDetector, rotateCut, moveCut);
	//four
	dI = 3;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* alHe3 = new G4SubtractionSolid("alHe3", alHe2, cutDetector, rotateCut, moveCut);
	//five
	dI = 4;
	dDist = ptrPd[dI] + detectorHalfLengthZ;
	dPhi = ptrPt[dI];
	dTheta = ptrPp[dI];
	moveCut.setX(dDist * cos(dPhi) * sin(dTheta));
	moveCut.setY(dDist * sin(dPhi) * sin(dTheta));
	moveCut.setZ(dDist *      1.0     * cos(dTheta));
	oriPhi = dPhi + ptrOt[dI] - M_PI/2;
	oriTheta = dTheta + ptrOp[dI];
	yprimeaxis.set(cos(oriPhi), sin(oriPhi), 0);
	rotateCut->set(yprimeaxis, oriTheta);
	G4SubtractionSolid* aluminumHemisphere = new G4SubtractionSolid("aluminumHemisphere", alHe3, cutDetector, rotateCut, moveCut);
	//  //end

	//logical volume
	if(fAluminumHemisphereLog == nullptr) {
		fAluminumHemisphereLog = new G4LogicalVolume(aluminumHemisphere, material, "aluminumHemisphereLog", 0, 0, 0);
		fAluminumHemisphereLog->SetVisAttributes(aluminumHemisphereVisAtt);
	}

	return 1;
}//end ::end AddAluminumHemisphere


//Add screws
G4int DetectionSystemPaces::AddScrews() {
	return 1;
}//end ::end AddScrews
