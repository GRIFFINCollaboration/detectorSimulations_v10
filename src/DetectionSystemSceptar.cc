#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemSceptar.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

DetectionSystemSceptar::DetectionSystemSceptar() :
	// LogicalVolumes
	fSquareMylarLog(0),
	fAngledMylarLog(0),
	fSquareScintillatorLog(0),
	fAngledScintillatorLog(0),
	fDelrinShellLog(0),
	fDelrinShell2Log(0),
	fHevimetShellLog(0)
{

	/////////////////////////////////////////////////////////////////////
	// Detector Properties
	/////////////////////////////////////////////////////////////////////
	fSeparateHemispheres = 1.0*mm;
	fConvert = 0.0254*m; //Needed to convert inches to meters

	fSquareScintillatorLength = 1.194*fConvert;
	fSquareScintillatorWidth = 0.449*fConvert;
	fSquareScintillatorThickness = fConvert/16.0;
	fAngledScintillatorLength  = 1.003*fConvert;
	fAngledScintillatorLongWidth = 1.045*fConvert;
	fAngledScintillatorShortWidth = 0.690*fConvert;
	fAngledScintillatorThickness = fConvert/16.0;
	fMylarThickness = 0.01*mm;
	fScintGap = 0.75*mm;
	fScintAngle1 = 36.0*deg;
	fScintAngleMove = 18.0*deg;

	//*********This angle is the tapering along the length (10.05)
	fScintAngle2 = atan((fAngledScintillatorLongWidth -fAngledScintillatorShortWidth)/(2.0*fAngledScintillatorLength));

	//*********This angle is the slight tapering of the sides (34.765))
	fScintAngle3 = acos(cos(fScintAngle1)/cos(fScintAngle2));

	//%%%%% This is the angle that defines the packed
	//%%%%% configuration of the angledScints
	//%%%%%
	fScintAngle4 = asin((fAngledScintillatorLongWidth -fAngledScintillatorShortWidth)/(2.0*fAngledScintillatorLength*tan(fScintAngle1)));

	//%%%%% This is the radial distance from the center of the
	//%%%%% pentagon formed by the squareDetectors to the
	//%%%%% center of any one of the squareDetectors in that pentagon
	//%%%%% The gap between scintillators and the mylar coating are taken into account
	//%%%%%
	fSquareScintRadialDistance = ((fSquareScintillatorLength +2.0*fMylarThickness +fScintGap/cos(fScintAngle1))*tan((M_PI/2.0)-fScintAngle1)-(fSquareScintillatorThickness+2.0*fMylarThickness))/2.0;

	//%%%%% This is the radial distance from the axis of the
	//%%%%% pentagon formed by the angledDetectors to the
	//%%%%% center of any one of the angledDetectors in that pentagon
	//%%%%% The gap between the scintillators and the mylar coating are taken into account
	//%%%%%
	fAngledScintRadialDistance = 0.5*(((fAngledScintillatorShortWidth +2.0*fMylarThickness +fScintGap/cos(fScintAngle1))/tan(fScintAngle1)) +(fAngledScintillatorLength+2.0*fMylarThickness)*sin(fScintAngle4)-(fAngledScintillatorThickness +2.0*fMylarThickness)*cos(fScintAngle4));

	//%%%%% This is the distance from the origin, along-y that
	//%%%%% the angledScints need to be moved to give
	//%%%%% proper alignment with the squareScints
	//%%%%% The mylar coating is taken into account
	//%%%%%
	fAngledScintMoveBack = fSquareScintillatorWidth -(fSquareScintRadialDistance -0.5*fSquareScintillatorThickness -0.5*fAngledScintillatorLongWidth/tan(fScintAngle1))/tan(atan((fSquareScintRadialDistance -0.5*fSquareScintillatorThickness)/fSquareScintillatorWidth))+0.5*fAngledScintillatorLength*cos(fScintAngle4);

	fDelrinInnerRadius = 3.275*fConvert;
	fDelrinOuterRadius = 3.475*fConvert;
	fDelrin2InnerRadius = 8.90*cm;
	fDelrin2OuterRadius = 9.90*cm;
	fHevimetInnerRadius = 9.906*cm;
	fHevimetOuterRadius = 12.446*cm;
	fDelrinHoleRadius = 1.25*fConvert;
}

DetectionSystemSceptar::~DetectionSystemSceptar() {
	// LogicalVolumes
	delete fSquareMylarLog;
	delete fAngledMylarLog;
	delete fSquareScintillatorLog;
	delete fAngledScintillatorLog;
	delete fDelrinShellLog;
	delete fDelrinShell2Log;
	delete fHevimetShellLog;
}

G4int DetectionSystemSceptar::Build() {
	// Build assembly volume
	fAssembly = new G4AssemblyVolume();
	fAssemblySquare = new G4AssemblyVolume();
	fAssemblyAngled = new G4AssemblyVolume();
	fAssemblySquareSD = new G4AssemblyVolume();
	fAssemblyAngledSD = new G4AssemblyVolume();

	ConstructScintillator();
	//  ConstructDelrinShell();
	//  Construct2ndDelrinShell();
	//  ConstructHevimetShell();

	return 1;
}

G4int DetectionSystemSceptar::PlaceDetector(G4LogicalVolume* expHallLog, G4int detectorNumber) {

	//G4RotationMatrix* rotateNull = new G4RotationMatrix;
	G4ThreeVector moveNull(0.0,0.0,0.0);

	G4double extra = 2.0*fMylarThickness;

	G4RotationMatrix* rotate = new G4RotationMatrix;
	G4ThreeVector move;

	G4RotationMatrix* rotateSquareScint1 = new G4RotationMatrix;
	rotateSquareScint1->rotateX(M_PI/2.0);
	rotateSquareScint1->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint1(	-fSquareScintRadialDistance ,
			0,
			-(fSquareScintillatorWidth +extra + fSeparateHemispheres)/2.0);

	G4RotationMatrix* rotateSquareScint2 = new G4RotationMatrix;
	rotateSquareScint2->rotateY(-1.0*(2.0*fScintAngle1));
	rotateSquareScint2->rotateX(M_PI/2.0);
	rotateSquareScint2->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint2(	-fSquareScintRadialDistance*cos(2.0*fScintAngle1),
			fSquareScintRadialDistance*sin(2.0*fScintAngle1),
			-(fSquareScintillatorWidth +extra + fSeparateHemispheres) / 2.0) ;

	G4RotationMatrix* rotateSquareScint3 = new G4RotationMatrix;
	rotateSquareScint3->rotateY(-1.0*(M_PI -fScintAngle1));
	rotateSquareScint3->rotateX(M_PI/2.0);
	rotateSquareScint3->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint3(	fSquareScintRadialDistance*cos(fScintAngle1),
			fSquareScintRadialDistance*sin(fScintAngle1),
			-(fSquareScintillatorWidth +extra + fSeparateHemispheres) / 2.0) ;

	G4RotationMatrix* rotateSquareScint4 = new G4RotationMatrix;
	rotateSquareScint4->rotateY(-1.0*(fScintAngle1 -M_PI));
	rotateSquareScint4->rotateX(M_PI/2.0);
	rotateSquareScint4->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint4(	fSquareScintRadialDistance*cos(fScintAngle1),
			-fSquareScintRadialDistance*sin(fScintAngle1),
			-(fSquareScintillatorWidth +extra + fSeparateHemispheres) / 2.0) ;

	G4RotationMatrix* rotateSquareScint5 = new G4RotationMatrix;
	rotateSquareScint5->rotateY(-1.0*(-2.0*fScintAngle1));
	rotateSquareScint5->rotateX(M_PI/2.0);
	rotateSquareScint5->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint5(	-fSquareScintRadialDistance*cos(2.0*fScintAngle1),
			-fSquareScintRadialDistance*sin(2.0*fScintAngle1),
			-(fSquareScintillatorWidth +extra + fSeparateHemispheres) / 2.0) ;

	G4RotationMatrix* rotateSquareScint6 = new G4RotationMatrix;
	rotateSquareScint6->rotateY(M_PI);
	rotateSquareScint6->rotateX(M_PI/2.0);
	rotateSquareScint6->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint6;
	moveSquareScint6 = - 1.0 * moveSquareScint1;

	G4RotationMatrix* rotateSquareScint7 = new G4RotationMatrix;
	rotateSquareScint7->rotateY(-1.0*(2.0*fScintAngle1 -M_PI));
	rotateSquareScint7->rotateX(M_PI/2.0);
	rotateSquareScint7->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint7;
	moveSquareScint7 = -1.0*moveSquareScint2;

	G4RotationMatrix* rotateSquareScint8 = new G4RotationMatrix;
	rotateSquareScint8->rotateY(-1.0*(-fScintAngle1));
	rotateSquareScint8->rotateX(M_PI/2.0);
	rotateSquareScint8->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint8;
	moveSquareScint8 = -1.0*moveSquareScint3;

	G4RotationMatrix* rotateSquareScint9 = new G4RotationMatrix;
	rotateSquareScint9->rotateY(-1.0*(fScintAngle1));
	rotateSquareScint9->rotateX(M_PI/2.0);
	rotateSquareScint9->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint9;
	moveSquareScint9 = -1.0*moveSquareScint4;

	G4RotationMatrix* rotateSquareScint10 = new G4RotationMatrix;
	rotateSquareScint10->rotateY(-1.0*(-2.0*fScintAngle1 -M_PI));
	rotateSquareScint10->rotateX(M_PI/2.0);
	rotateSquareScint10->rotateZ(M_PI/2.0);
	G4ThreeVector moveSquareScint10;
	moveSquareScint10 = -1.0*moveSquareScint5;

	G4RotationMatrix* rotateAngledScint1 = new G4RotationMatrix;
	rotateAngledScint1->rotateX(-1.0*(fScintAngle4 -M_PI/2.0));
	rotateAngledScint1->rotateX(M_PI/2.0);
	rotateAngledScint1->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint1(	-fAngledScintRadialDistance,
			0,
			-fAngledScintMoveBack);

	G4RotationMatrix* rotateAngledScint2 = new G4RotationMatrix;
	rotateAngledScint2->rotateX(-1.0*(fScintAngle4));
	rotateAngledScint2->rotateZ(-1.0*(-2.0*fScintAngle1));
	rotateAngledScint2->rotateX(M_PI);
	rotateAngledScint2->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint2(	-fAngledScintRadialDistance*cos(2.0*fScintAngle1),
			fAngledScintRadialDistance*sin(2.0*fScintAngle1),
			-fAngledScintMoveBack);

	G4RotationMatrix* rotateAngledScint3 = new G4RotationMatrix;
	rotateAngledScint3->rotateX(-1.0*(fScintAngle4));
	rotateAngledScint3->rotateZ(-1.0*(fScintAngle1 -M_PI));
	rotateAngledScint3->rotateX(M_PI);
	rotateAngledScint3->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint3(	fAngledScintRadialDistance*cos(fScintAngle1),
			fAngledScintRadialDistance*sin(fScintAngle1),
			-fAngledScintMoveBack);

	G4RotationMatrix* rotateAngledScint4 = new G4RotationMatrix;
	rotateAngledScint4->rotateX(-1.0*(fScintAngle4));
	rotateAngledScint4->rotateZ(-1.0*(M_PI -fScintAngle1));
	rotateAngledScint4->rotateX(M_PI);
	rotateAngledScint4->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint4(	fAngledScintRadialDistance*cos(fScintAngle1),
			-fAngledScintRadialDistance*sin(fScintAngle1),
			-fAngledScintMoveBack);

	G4RotationMatrix* rotateAngledScint5 = new G4RotationMatrix;
	rotateAngledScint5->rotateX(-1.0*(fScintAngle4));
	rotateAngledScint5->rotateZ(-1.0*(2.0*fScintAngle1));
	rotateAngledScint5->rotateX(M_PI);
	rotateAngledScint5->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint5(	-fAngledScintRadialDistance*cos(2.0*fScintAngle1),
			-fAngledScintRadialDistance*sin(2.0*fScintAngle1),
			-fAngledScintMoveBack);

	G4RotationMatrix* rotateAngledScint6 = new G4RotationMatrix;
	rotateAngledScint6->rotateX(-1.0*(fScintAngle4 +M_PI/2.0));
	rotateAngledScint6->rotateX(M_PI/2.0);
	rotateAngledScint6->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint6;
	moveAngledScint6 = -1.0*moveAngledScint1;

	G4RotationMatrix* rotateAngledScint7 = new G4RotationMatrix;
	rotateAngledScint7->rotateX(-1.0*(fScintAngle4));
	rotateAngledScint7->rotateZ(-1.0*(2.0*fScintAngle1));
	rotateAngledScint7->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint7;
	moveAngledScint7 = -1.0*moveAngledScint2;

	G4RotationMatrix* rotateAngledScint8 = new G4RotationMatrix;
	rotateAngledScint8->rotateX(-1.0*(fScintAngle4));
	rotateAngledScint8->rotateZ(-1.0*(M_PI -fScintAngle1));
	rotateAngledScint8->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint8;
	moveAngledScint8 = -1.0*moveAngledScint3;

	G4RotationMatrix* rotateAngledScint9 = new G4RotationMatrix;
	rotateAngledScint9->rotateX(-1.0*(fScintAngle4));
	rotateAngledScint9->rotateZ(-1.0*(fScintAngle1 -M_PI));
	rotateAngledScint9->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint9;
	moveAngledScint9 = -1.0*moveAngledScint4;

	G4RotationMatrix* rotateAngledScint10 = new G4RotationMatrix;
	rotateAngledScint10->rotateX(-1.0*(fScintAngle4));
	rotateAngledScint10->rotateZ(-1.0*(-2.0*fScintAngle1));
	rotateAngledScint10->rotateZ(M_PI/2.0);
	G4ThreeVector moveAngledScint10;
	moveAngledScint10 = -1.0*moveAngledScint5;

	for(G4int detNum = 0; detNum < detectorNumber; detNum++) {
		if(detNum == 0) {
			rotate = rotateSquareScint6;
			move = moveSquareScint6;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 1) {
			rotate = rotateSquareScint7;
			move = moveSquareScint7;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 2) {
			rotate = rotateSquareScint8;
			move = moveSquareScint8;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 3) {
			rotate = rotateSquareScint9;
			move = moveSquareScint9;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 4) {
			rotate = rotateSquareScint10;
			move = moveSquareScint10;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 5) {
			rotate = rotateAngledScint6;
			move = moveAngledScint6;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 6) {
			rotate = rotateAngledScint7;
			move = moveAngledScint7;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 7) {
			rotate = rotateAngledScint8;
			move = moveAngledScint8;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 8) {
			rotate = rotateAngledScint9;
			move = moveAngledScint9;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 9) {
			rotate = rotateAngledScint10;
			move = moveAngledScint10;
			move.rotateZ(fScintAngleMove);
			rotate->rotateZ(fScintAngleMove);
		} else if(detNum == 10) {
			rotate = rotateSquareScint1;
			move = moveSquareScint1;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		} else if(detNum == 11) {
			rotate = rotateSquareScint2;
			move = moveSquareScint2;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		} else if(detNum == 12) {
			rotate = rotateSquareScint3;
			move = moveSquareScint3;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		} else if(detNum == 13) {
			rotate = rotateSquareScint4;
			move = moveSquareScint4;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		} else if(detNum == 14) {
			rotate = rotateSquareScint5;
			move = moveSquareScint5;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		} else if(detNum == 15) {
			rotate = rotateAngledScint1;
			move = moveAngledScint1;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		} else if(detNum == 16) {
			rotate = rotateAngledScint2;
			move = moveAngledScint2;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		} else if(detNum == 17) {
			rotate = rotateAngledScint3;
			move = moveAngledScint3;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		} else if(detNum == 18) {
			rotate = rotateAngledScint4;
			move = moveAngledScint4;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		} else if(detNum == 19) {
			rotate = rotateAngledScint5;
			move = moveAngledScint5;
			move.rotateZ(fScintAngleMove+fScintAngle1);
			rotate->rotateZ(fScintAngleMove+fScintAngle1);
		}



		if((detNum < 5) || (detNum >= 10 && detNum < 15)) {
			fAssemblySquareSD->MakeImprint(expHallLog, move, rotate, 0);
			fAssemblySquare->MakeImprint(expHallLog, move, rotate, 0);
		}
		if((detNum >= 5 && detNum < 10) || (detNum >= 15)) {
			fAssemblyAngledSD->MakeImprint(expHallLog, move, rotate, 0);
			fAssemblyAngled->MakeImprint(expHallLog, move, rotate, 0);
		}

	}

	return 1;
}

G4int DetectionSystemSceptar::ConstructScintillator() {
	G4Material* materialMylar = G4Material::GetMaterial("Mylar");
	if(!materialMylar) {
		G4cout<<" ----> Material "<<"Mylar"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}
	G4Material* materialBc404 = G4Material::GetMaterial("BC404");
	if(!materialBc404) {
		G4cout<<" ----> Material "<<"BC404"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	G4ThreeVector moveNew;

	// Set visualization attributes
	G4VisAttributes* scintillatorVisAtt = new G4VisAttributes(G4Colour(0.0,0.3,0.7));
	scintillatorVisAtt->SetVisibility(true);
	G4VisAttributes* mylarVisAtt = new G4VisAttributes(G4Colour(0.6,0.6,0.6));
	mylarVisAtt->SetVisibility(true);

	//G4double extra = 2.0*fMylarThickness;

	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//&&&&&&& Making the square-ish shaped scintillators 1-5 are &&&&&&&
	//&&&&&&& neg-y of origin, while 6-10 are in pos-y of origin &&&&&&&
	//&&&&&&& NOTE: this numbering scheme does not reflect the   &&&&&&&
	//&&&&&&&    actual way the detectors will be numbered       &&&&&&&
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	//******* Making the Mylar coating

	//  G4Trd* squareMylar = squareMylar();

	G4SubtractionSolid* squareMylar = SquareMylarWithCut();

	fSquareMylarLog = new G4LogicalVolume(squareMylar, materialMylar, "squareScintillatorLog", 0, 0, 0);
	fSquareMylarLog->SetVisAttributes(mylarVisAtt);

	//******* Making the square scint and putting it in the Mylar

	G4Trd* squareScintillator = SquareScintillator();

	fSquareScintillatorLog = new G4LogicalVolume(squareScintillator, materialBc404, "sceptarSquareScintillatorLog", 0, 0, 0);
	fSquareScintillatorLog->SetVisAttributes(scintillatorVisAtt);

	G4RotationMatrix* rotateNull = new G4RotationMatrix;
	G4ThreeVector moveNull(0.0,0.0,0.0);

	fAssemblySquare->AddPlacedVolume(fSquareMylarLog, moveNull, rotateNull);
	fAssemblySquareSD->AddPlacedVolume(fSquareScintillatorLog, moveNull, rotateNull);
	//******* Making the Mylar coating

	G4SubtractionSolid* angledMylar = AngledMylarWithCut();

	fAngledMylarLog = new G4LogicalVolume(angledMylar, materialMylar, "angledMylarLog", 0, 0, 0);
	fAngledMylarLog->SetVisAttributes(mylarVisAtt);

	//******* Making the angled scint and putting it in the Mylar

	G4SubtractionSolid* angledScintillator = AngledScintillator();

	fAngledScintillatorLog = new G4LogicalVolume(angledScintillator, materialBc404, "sceptarAngledScintillatorLog", 0, 0, 0);
	fAngledScintillatorLog->SetVisAttributes(scintillatorVisAtt);

	fAssemblyAngled->AddPlacedVolume(fAngledMylarLog, moveNull, rotateNull);
	fAssemblyAngledSD->AddPlacedVolume(fAngledScintillatorLog, moveNull, rotateNull);

	return 1;
}


//*****************************************************************
//ConstructDelrinShell builds a spherical shell of Delrin around
//the detectors, with two holes in it for the beam line
//*****************************************************************

G4int DetectionSystemSceptar::ConstructDelrinShell() {
	G4Material* materialDelrin = G4Material::GetMaterial("Delrin");
	if(!materialDelrin) {
		G4cout<<" ----> Material "<<"Delrin"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	G4VisAttributes* DelrinVisAtt = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
	DelrinVisAtt->SetVisibility(true);

	G4SubtractionSolid* shell = DelrinShell();

	fDelrinShellLog = new G4LogicalVolume(shell, materialDelrin, "DelrinShellLog", 0, 0, 0);
	fDelrinShellLog->SetVisAttributes(DelrinVisAtt);

	G4RotationMatrix* rotateNull = new G4RotationMatrix;
	G4ThreeVector moveNull(0.0,0.0,0.0);

	fAssembly->AddPlacedVolume(fDelrinShellLog, moveNull, rotateNull);

	return 1;
} //end ::ConstructDelrinShell


//*******************************************************************
//Construct2ndDelrinShell builds a spherical shell of Delrin around
//the first sphere of Delrin, with two holes in it for the beam line.
//This is to simulate the Delrin in front of the germanium detectors
//*******************************************************************

G4int DetectionSystemSceptar::Construct2ndDelrinShell() {
	G4Material* materialDelrin = G4Material::GetMaterial("Delrin");
	if(!materialDelrin) {
		G4cout<<" ----> Material "<<"Delrin"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	G4VisAttributes* Delrin2VisAtt = new G4VisAttributes(G4Colour(1,1,1));
	Delrin2VisAtt->SetVisibility(true);

	G4SubtractionSolid* shell2 = DelrinShell2();

	fDelrinShell2Log = new G4LogicalVolume(shell2, materialDelrin, "DelrinShell2Log", 0, 0, 0);
	fDelrinShell2Log->SetVisAttributes(Delrin2VisAtt);

	G4RotationMatrix* rotateNull = new G4RotationMatrix;
	G4ThreeVector moveNull(0.0,0.0,0.0);

	fAssembly->AddPlacedVolume(fDelrinShell2Log, moveNull, rotateNull);

	return 1;

} //end ::Construct2ndDelrinShell


//*******************************************************************
//ConstructHevimetShell builds a spherical shell of Hevimet around
//the two spheres of Delrin, with two holes in it for the beam line.
//This is to simulate the Hevimet in front of the germanium detectors
//*******************************************************************

G4int DetectionSystemSceptar::ConstructHevimetShell() {
	G4Material* materialHevimet = G4Material::GetMaterial("Hevimet");
	if(!materialHevimet) {
		G4cout<<" ----> Material "<<"Hevimet"<<" not found, cannot build the detector shell! "<<G4endl;
		return 0;
	}

	G4VisAttributes* HevimetVisAtt = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
	HevimetVisAtt->SetVisibility(true);

	G4SubtractionSolid* shell3 = HevimetShell();

	fHevimetShellLog = new G4LogicalVolume(shell3, materialHevimet, "HevimetShellLog", 0, 0, 0);
	fHevimetShellLog->SetVisAttributes(HevimetVisAtt);

	G4RotationMatrix* rotateNull = new G4RotationMatrix;
	G4ThreeVector moveNull(0.0,0.0,0.0);

	fAssembly->AddPlacedVolume(fHevimetShellLog, moveNull, rotateNull);

	return 1;

} //end ::Construct2ndDelrinShell


///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Trd* DetectionSystemSceptar::SquareMylar() {
	G4double extra = 2.0*fMylarThickness;

	G4double halfLengthLongerX = (fSquareScintillatorLength +extra)/2.0;
	G4double halfLengthShorterX = halfLengthLongerX - tan(fScintAngle1)*(fSquareScintillatorThickness +extra);
	G4double halfLengthY = (fSquareScintillatorWidth +extra)/2.0;
	G4double halfLengthZ = (fSquareScintillatorThickness +extra)/2.0;

	G4Trd* squareMylar = new G4Trd("squareMylar", halfLengthLongerX, halfLengthShorterX,halfLengthY, halfLengthY, halfLengthZ);

	return squareMylar;
} //end ::squareMylar

G4SubtractionSolid* DetectionSystemSceptar::SquareMylarWithCut() {
	G4double extra = 2.0*fMylarThickness;

	G4double halfLengthLongerX = (fSquareScintillatorLength +extra)/2.0;
	G4double halfLengthShorterX = halfLengthLongerX - tan(fScintAngle1)*(fSquareScintillatorThickness +extra);
	G4double halfLengthY = (fSquareScintillatorWidth +extra)/2.0;
	G4double halfLengthZ = (fSquareScintillatorThickness +extra)/2.0;

	G4Trd* squareMylar = new G4Trd("squareMylar", halfLengthLongerX, halfLengthShorterX,halfLengthY, halfLengthY, halfLengthZ);

	G4Trd* scintCut = SquareScintillator();

	G4RotationMatrix* rotateCut = new G4RotationMatrix;
	G4ThreeVector moveCut(0, 0, 0);

	G4SubtractionSolid* squareMylarWithCut = new G4SubtractionSolid("squareMylarWithCut", squareMylar, scintCut, rotateCut, moveCut);

	return squareMylarWithCut;
} //end ::squareMylar

G4SubtractionSolid* DetectionSystemSceptar::AngledMylar() {
	G4double extra = 2.0*fMylarThickness;

	G4double halfLengthNegX = (fAngledScintillatorLongWidth +extra)/2.0;
	G4double halfLengthPosX = (fAngledScintillatorShortWidth +extra)/2.0;
	G4double halfLengthZ = (fAngledScintillatorLength +extra)/2.0;
	G4double halfLengthY = (fAngledScintillatorThickness +extra)/2.0;

	G4Trd* mylar = new G4Trd("mylar", halfLengthNegX, halfLengthPosX, halfLengthY, halfLengthY, halfLengthZ);

	G4double halfThickness = (fAngledScintillatorThickness +extra)/(2.0*cos(fScintAngle3));
	G4double halfLength = fAngledScintillatorLength +extra;

	G4Box* box = new G4Box("box", halfThickness, halfThickness, halfLength);

	G4RotationMatrix* rotateBox1 = new G4RotationMatrix;
	rotateBox1->rotateY(fScintAngle2);
	rotateBox1->rotateZ(-fScintAngle3);

	G4ThreeVector moveBox1((fAngledScintillatorShortWidth +extra +tan(fScintAngle2)*(fAngledScintillatorLength +extra) -tan(fScintAngle3)*(fAngledScintillatorThickness +extra)/cos(fScintAngle2))/2.0 +cos(fScintAngle3)*halfThickness/cos(fScintAngle2), halfThickness*sin(fScintAngle3), 0);

	G4SubtractionSolid* angledMylar1 = new G4SubtractionSolid("angledMylar1", mylar, box, rotateBox1, moveBox1);

	G4RotationMatrix* rotateBox2 = new G4RotationMatrix;
	rotateBox2->rotateY(-fScintAngle2);
	rotateBox2->rotateZ(fScintAngle3);

	G4ThreeVector moveBox2(-(fAngledScintillatorShortWidth +extra +tan(fScintAngle2)*(fAngledScintillatorLength +extra) -tan(fScintAngle3)*(fAngledScintillatorThickness +extra)/cos(fScintAngle2))/2.0 -cos(fScintAngle3)*halfThickness/cos(fScintAngle2), halfThickness*sin(fScintAngle3), 0);

	G4SubtractionSolid* angledMylar = new G4SubtractionSolid("angledMylar", angledMylar1, box, rotateBox2, moveBox2);

	return angledMylar;

} //end ::angledMylar


G4SubtractionSolid* DetectionSystemSceptar::AngledMylarWithCut() {
	G4double extra = 2.0*fMylarThickness;

	G4double halfLengthNegX = (fAngledScintillatorLongWidth +extra)/2.0;
	G4double halfLengthPosX = (fAngledScintillatorShortWidth +extra)/2.0;
	G4double halfLengthZ = (fAngledScintillatorLength +extra)/2.0;
	G4double halfLengthY = (fAngledScintillatorThickness +extra)/2.0;

	G4Trd* mylar = new G4Trd("mylar", halfLengthNegX, halfLengthPosX, halfLengthY, halfLengthY, halfLengthZ);

	G4double halfThickness = (fAngledScintillatorThickness +extra)/(2.0*cos(fScintAngle3));
	G4double halfLength = fAngledScintillatorLength +extra;

	G4Box* box = new G4Box("box", halfThickness, halfThickness, halfLength);

	G4RotationMatrix* rotateBox1 = new G4RotationMatrix;
	rotateBox1->rotateY(fScintAngle2);
	rotateBox1->rotateZ(-fScintAngle3);

	G4ThreeVector moveBox1((fAngledScintillatorShortWidth +extra +tan(fScintAngle2)*(fAngledScintillatorLength +extra) -tan(fScintAngle3)*(fAngledScintillatorThickness +extra)/cos(fScintAngle2))/2.0 +cos(fScintAngle3)*halfThickness/cos(fScintAngle2), halfThickness*sin(fScintAngle3), 0);

	G4SubtractionSolid* angledMylar1 = new G4SubtractionSolid("angledMylar1", mylar, box, rotateBox1, moveBox1);

	G4RotationMatrix* rotateBox2 = new G4RotationMatrix;
	rotateBox2->rotateY(-fScintAngle2);
	rotateBox2->rotateZ(fScintAngle3);

	G4ThreeVector moveBox2(-(fAngledScintillatorShortWidth +extra +tan(fScintAngle2)*(fAngledScintillatorLength +extra) -tan(fScintAngle3)*(fAngledScintillatorThickness +extra)/cos(fScintAngle2))/2.0 -cos(fScintAngle3)*halfThickness/cos(fScintAngle2), halfThickness*sin(fScintAngle3), 0);

	G4SubtractionSolid* angledMylar = new G4SubtractionSolid("angledMylar", angledMylar1, box, rotateBox2, moveBox2);

	G4SubtractionSolid* scintCut = AngledScintillator();

	G4RotationMatrix* rotateCut = new G4RotationMatrix;
	G4ThreeVector moveCut(0, 0, 0);

	G4SubtractionSolid* angledMylarWithCut = new G4SubtractionSolid("angledMylarWithCut", angledMylar, scintCut, rotateCut, moveCut);

	return angledMylarWithCut;

} //end ::angledMylar


G4Trd* DetectionSystemSceptar::SquareScintillator() {
	G4double halfLengthLongerX = fSquareScintillatorLength/2.0;
	G4double halfLengthShorterX = halfLengthLongerX - tan(fScintAngle1)*fSquareScintillatorThickness;
	G4double halfLengthY = fSquareScintillatorWidth/2.0;
	G4double halfLengthZ = fSquareScintillatorThickness/2.0;

	G4Trd* squareScintillator = new G4Trd("squareScintillator", halfLengthLongerX, halfLengthShorterX,halfLengthY, halfLengthY, halfLengthZ);

	return squareScintillator;
} //end ::squareScintillator


G4SubtractionSolid* DetectionSystemSceptar::AngledScintillator() {
	G4double halfLengthNegX = fAngledScintillatorLongWidth/2.0;
	G4double halfLengthPosX = fAngledScintillatorShortWidth/2.0;
	G4double halfLengthZ = fAngledScintillatorLength/2.0;
	G4double halfLengthY = fAngledScintillatorThickness/2.0;

	G4Trd* scint = new G4Trd("scint", halfLengthNegX, halfLengthPosX, halfLengthY, halfLengthY, halfLengthZ);

	G4double halfThickness = fAngledScintillatorThickness/(2.0*cos(fScintAngle3));
	G4double halfLength = fAngledScintillatorLength;

	G4Box* box = new G4Box("box", halfThickness, halfThickness, halfLength);

	G4RotationMatrix* rotateBox1 = new G4RotationMatrix;
	rotateBox1->rotateY(fScintAngle2);
	rotateBox1->rotateZ(-fScintAngle3);

	G4ThreeVector moveBox1((fAngledScintillatorShortWidth +tan(fScintAngle2)*fAngledScintillatorLength -tan(fScintAngle3)*fAngledScintillatorThickness/cos(fScintAngle2))/2.0 +cos(fScintAngle3)*halfThickness/cos(fScintAngle2), halfThickness*sin(fScintAngle3), 0);

	G4SubtractionSolid* angledScint1 = new G4SubtractionSolid("angledScint1", scint, box, rotateBox1, moveBox1);

	G4RotationMatrix* rotateBox2 = new G4RotationMatrix;
	rotateBox2->rotateY(-fScintAngle2);
	rotateBox2->rotateZ(fScintAngle3);

	G4ThreeVector moveBox2(-(fAngledScintillatorShortWidth +tan(fScintAngle2)*fAngledScintillatorLength -tan(fScintAngle3)*fAngledScintillatorThickness/cos(fScintAngle2))/2.0 -cos(fScintAngle3)*halfThickness/cos(fScintAngle2), halfThickness*sin(fScintAngle3), 0);

	G4SubtractionSolid* angledScint = new G4SubtractionSolid("angledScint", angledScint1, box, rotateBox2, moveBox2);

	return angledScint;

} //end ::angledScintillator


G4SubtractionSolid* DetectionSystemSceptar::DelrinShell() {
	G4Sphere* sphere = new G4Sphere("sphere", fDelrinInnerRadius, fDelrinOuterRadius, 0, 2.0*M_PI, 0, M_PI);

	G4Tubs* chopTub = new G4Tubs("chopTub", 0, fDelrinHoleRadius, fDelrinOuterRadius +1.0, 0, 2.0*M_PI);

	G4RotationMatrix* rotateChopTub = new G4RotationMatrix;
	rotateChopTub->rotateX(-M_PI/2.0);

	G4ThreeVector moveNull(0.0,0.0,0.0);

	G4SubtractionSolid* delrinShell = new G4SubtractionSolid("DelrinShell", sphere, chopTub, rotateChopTub, moveNull);

	return delrinShell;
} //end ::DelrinShell


G4SubtractionSolid* DetectionSystemSceptar::DelrinShell2() {
	G4Sphere* sphere = new G4Sphere("sphere", fDelrin2InnerRadius, fDelrin2OuterRadius, 0, 2.0*M_PI, 0, M_PI);

	G4Tubs* chopTub = new G4Tubs("chopTub", 0, fDelrinHoleRadius, fDelrin2OuterRadius +1.0, 0, 2.0*M_PI);

	G4RotationMatrix* rotateChopTub = new G4RotationMatrix;
	rotateChopTub->rotateX(-M_PI/2.0);

	G4ThreeVector moveNull(0.0,0.0,0.0);

	G4SubtractionSolid* delrinShell2 = new G4SubtractionSolid("DelrinShell2", sphere, chopTub, rotateChopTub, moveNull);

	return delrinShell2;
} //end ::DelrinShell2


G4SubtractionSolid* DetectionSystemSceptar::HevimetShell() {
	G4Sphere* sphere = new G4Sphere("sphere", fHevimetInnerRadius, fHevimetOuterRadius,0, 2.0*M_PI, 0, M_PI);

	G4Tubs* chopTub = new G4Tubs("chopTub", 0, fDelrinHoleRadius, fHevimetOuterRadius +1.0,0, 2.0*M_PI);

	G4RotationMatrix* rotateChopTub = new G4RotationMatrix;
	rotateChopTub->rotateX(-M_PI/2.0);

	G4ThreeVector moveNull(0.0,0.0,0.0);

	G4SubtractionSolid* hevimetShell = new G4SubtractionSolid("HevimetShell", sphere, chopTub, rotateChopTub, moveNull);

	return hevimetShell;
} //end ::DelrinShell2

