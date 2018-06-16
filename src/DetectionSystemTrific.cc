#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"

#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemTrific.hh"
#include "HistoManager.hh" // for Trific histogram array when placing detector

DetectionSystemTrific::DetectionSystemTrific(G4double setpressure) 
{
	//-----------------------------//
	// Materials     			   //
	//-----------------------------//
	fWireMaterial    = "Tungsten";
	fPCBMaterial     = "Acrylic";//Need to define G10
	fChamberMaterial = "G4_STAINLESS-STEEL"; 
	fWindowFrameMaterial = "G4_Al"; 
	fWindowMaterial = "Mylar";
	fWindowSurfaceMaterial = "G4_Al";
	fGasMaterial = "TrificCF4";


	G4double a, z, density, temperature, pressure;
	G4String name, symbol;
	G4int    ncomponents, natoms;
	// 	G4double fractionmass;

	G4Element* ElC = new G4Element(name="C", symbol="C", z=6., a = 12.01*g/mole);
	G4Element* ElF  = new G4Element(name="F" , symbol="F" , z=9. , a= 18.9984 *g/mole);

	temperature = CLHEP::STP_Temperature;
	pressure = (setpressure/760.)*CLHEP::STP_Pressure;
	density = (setpressure/760.)*3.72*mg/cm3;

	G4Material* TrificCF4 = new G4Material(name="TrificCF4",  density,ncomponents=2,kStateGas,temperature,pressure);
	// 	G4Material* TrificCF4= new G4Material(name="TrificCF4",density, ncomponents=2);
	TrificCF4->AddElement(ElC,natoms=1);
	TrificCF4->AddElement(ElF,natoms=4);

	G4cout<<G4endl<<"Created CF4 gas for pressure "<<setpressure<< " Torr ("<<1000.*pressure/bar<<" mbar)";
	G4cout<<G4endl<<"Density "<<density/(mg/cm3)<< " mg/cm3";

	// ----------------------------
	// Dimensions 
	// ----------------------------
	fPCBThickness = 2.3*mm;
	fGridAngle = 30*deg;
	fSpacerZThickness = 1.8*cm;
	fGridSpacing = fPCBThickness/cos(fGridAngle)+fSpacerZThickness;
	fWireD = 20*um;
	// 	fWireD = 4*mm;//easier to see for development
	fWirePitch = 2*mm;
	// 	fWirePitch = 10*mm;//simpler for development
	fChamberLength = 60*cm;
	fChamberInnerD = 12*cm;
	fChamberOuterD = 12.8*cm;
	fPCBInnerD = 6*cm;
	fPCBOuterD = 11.5*cm;
	fPCBCutDepth = 2*cm;
	fPCBCutHalfHeight = 2*cm;
	fRodD = 0.5*cm;
	fRodPlaceD = 9*cm;


	fTargetChamberZOffset = 0*cm;
	fChamberWindowZ = 4*cm;
	fWindowGridZ = 1.5*cm;
	fWindowWasherZ = 1*cm;
	fWindowTickness = 5*um;
	fWindowCoatingThickness = 0.9*um;
	fWindowInnerD = 4*cm;
	fWindowOuterD = 8*cm;
	fWindowChamberLength = 2*cm;

}

DetectionSystemTrific::~DetectionSystemTrific()
{   	

}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemTrific::BuildPlace(G4LogicalVolume* ExpHallLog)
{

	BuildPlacePipe(ExpHallLog);
	BuildGasVolume(ExpHallLog);

	return 1;
} // end Build


void DetectionSystemTrific::BuildPCB(){
	fGridPCB = new G4AssemblyVolume();

	///
	/// First we defined the volume of gas that is one "cell" i.e. one grid + spacer
	///

	// Define the material, return error if not found
	G4Material* material = G4Material::GetMaterial(fGasMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fGasMaterial<< " not found, cannot build the Trific detector! "<<G4endl;
		return;
	}

	G4double cell_rad =(fChamberInnerD/2)-0.1*mm;
	G4double tubelength =cell_rad*tan(fGridAngle)+fGridSpacing;

	G4double cutlength =tubelength/cos(fGridAngle);
	G4double cutwidth =0.1*mm+cell_rad;
	G4double cutheight =0.1*mm+cell_rad/cos(fGridAngle)+cutlength*tan(fGridAngle);

	G4double cutfrontZ =cutlength/cos(fGridAngle)+fGridSpacing;
	G4double cutbackZ =-cutlength/cos(fGridAngle);

	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateX(fGridAngle);

	G4Tubs* CellTube = new G4Tubs("CellTube",0,cell_rad, tubelength, 0, 360.*deg);
	G4Box* CutBox = new G4Box("CutBox", cutwidth, cutheight, cutlength);

	G4VSolid* Cell = new G4SubtractionSolid("Cell", CellTube, CutBox,rotate, G4ThreeVector(0,0,cutfrontZ));
	G4VSolid* CellB = new G4SubtractionSolid("CellB", Cell, CutBox,rotate, G4ThreeVector(0,0,cutbackZ));

	fGasCell = new G4LogicalVolume(CellB, material, "TrificGasCell", 0, 0, 0);
	fGasCell->SetVisAttributes (G4VisAttributes::Invisible); 


	//Original version was box with tube cut (so defined verticle and rotated later)
	//but for some reason that didnt work here, though that method is used for the PCB next


	///
	/// Next we construct the PCB and place it in the cell
	///

	// Define the material, return error if not found
	material = G4Material::GetMaterial(fPCBMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fPCBMaterial<< " not found, cannot build the Trific detector! "<<G4endl;
		return;
	}

	// PCB 
	G4double box_half_width = fPCBOuterD/2.;
	G4double box_half_length = fPCBOuterD/(2*cos(fGridAngle))-3*mm;
	G4double box_half_thickness = fPCBThickness/2.;

	G4Box* PCBbase = new G4Box("PCBbase", box_half_width, box_half_length, box_half_thickness);
	//   G4Box* PCBbase = new G4Box("PCBbase", 10*cm, 10*cm, 10*cm);
	G4Box* PCBCut = new G4Box("PCBcutout", fPCBCutDepth, fPCBCutHalfHeight, fPCBThickness);
	G4VSolid* CutPCBa = new G4SubtractionSolid("CutPCBa", PCBbase, PCBCut,0, G4ThreeVector(box_half_width,0.,0.));
	G4VSolid* CutPCBb = new G4SubtractionSolid("CutPCBb", CutPCBa, PCBCut,0, G4ThreeVector(-box_half_width,0.,0.));


	G4Tubs* RodHole = new G4Tubs("RodHole", 0, fRodD/2+0.2*mm, fPCBOuterD, 0, 360.*CLHEP::deg);
	G4Tubs* CentreHole = new G4Tubs("CentreHole", 0, fPCBInnerD/2, fPCBOuterD, 0, 360.*CLHEP::deg);
	G4Tubs* OuterCut = new G4Tubs("OuterCut", fPCBOuterD/2, fPCBOuterD, fPCBOuterD, 0, 360.*CLHEP::deg);


	//These 2 holes are actually cut perpendicular to the PCB but defining that geometry is much harder
	G4VSolid* CutPCBc = new G4SubtractionSolid("CutPCBc", CutPCBb, CentreHole,rotate, G4ThreeVector());
	G4VSolid* CutPCB = new G4SubtractionSolid("CutPCB", CutPCBc, OuterCut,rotate, G4ThreeVector());


	for(int i=0;i<4;i++){
		G4ThreeVector rodplace(fRodPlaceD/2,0,0);
		rodplace.rotateZ(45.*deg);
		rodplace.rotateZ(90.*i*deg);
		// 	rodplace.rotateX(fGridAngle);
		CutPCB = new G4SubtractionSolid("CutPCB", CutPCB, RodHole,rotate,rodplace);
	}


	G4LogicalVolume* fFullPCB = new G4LogicalVolume(CutPCB, material, "fFullPCB", 0, 0, 0);

	// Set visualization attributes
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	vis_att->SetVisibility(true);  
	fFullPCB->SetVisAttributes(vis_att);

	G4double pcbplacedistance = (box_half_thickness+fWireD)/cos(fGridAngle);
	G4ThreeVector moveP(0,0,pcbplacedistance);
	new G4PVPlacement(rotate, moveP, fFullPCB, "TrificPCB", fGasCell,false, 0);

	// 	fGridPCB->AddPlacedVolume(fFullPCB, moveP, rotate); //G4AssemblyVolume version

	///
	/// Next we add the wires next to/on the PCB in the cell
	///

	//Add wires
	material = G4Material::GetMaterial(fWireMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fPCBMaterial<< " not found, cannot build the Trific detector! "<<G4endl;
		return;
	}

	vis_att = new G4VisAttributes(G4Colour(0.0,0.5,1.0));
	vis_att->SetVisibility(true);  
	fFullPCB->SetVisAttributes(vis_att);

	double wireX=fWirePitch/2;
	rotate = new G4RotationMatrix;
	rotate->rotateX(90*deg);
	rotate->rotateX(fGridAngle);

	G4double wireplacedistance = (fWireD/2)/cos(fGridAngle);
	while(wireX<fPCBInnerD/2){
		G4double WireHalfLength=sqrt((fPCBInnerD*fPCBInnerD/4)-wireX*wireX)/cos(fGridAngle)+1*cm;
		G4Tubs* WireSolid = new G4Tubs("WireSolid",0, fWireD/2, WireHalfLength, 0, 360.*CLHEP::deg);
		G4LogicalVolume* WireLogical = new G4LogicalVolume(WireSolid, material, "WireLogical", 0, 0, 0);
		WireLogical->SetVisAttributes(vis_att);

		G4ThreeVector moveL(wireX,0,wireplacedistance);
		new G4PVPlacement(rotate, moveL, WireLogical, "TrificWire", fGasCell,false, 0);
		G4ThreeVector moveR(-wireX,0,wireplacedistance);
		new G4PVPlacement(rotate, moveR, WireLogical, "TrificWire", fGasCell,false, 0);

		wireX+=fWirePitch;
	}

	return;
}


void DetectionSystemTrific::PlacePCBs(G4LogicalVolume*  TrificGasVol){


	G4double Zpos=-fChamberLength/2+fChamberWindowZ+fWindowGridZ;
	for(int i=0;i<24;i++){
		G4ThreeVector move(0,0,Zpos);
		std::stringstream temp;
		temp<<"TrificGasCell_" <<i;
		new G4PVPlacement(0, move, fGasCell, temp.str(), TrificGasVol, false, 0);
		// 		fGridPCB->MakeImprint(TrificGasVol, move, rotate, i);
		Zpos+=fGridSpacing;
	}


}


void DetectionSystemTrific::BuildPlacePipe(G4LogicalVolume*  ExpHallLog){
	// Define the material, return error if not found
	G4Material* material = G4Material::GetMaterial(fChamberMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fPCBMaterial<< " not found, cannot build the Trific detector! "<<G4endl;
		return;
	}

	G4Tubs* Pipe = new G4Tubs("Pipe",fChamberInnerD/2, fChamberOuterD/2, fChamberLength/2, 0, 360.*CLHEP::deg);

	G4LogicalVolume* TrificPipe = new G4LogicalVolume(Pipe, material, "TrificPipe", 0, 0, 0);
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1,1));
	vis_att->SetVisibility(true);  
	TrificPipe->SetVisAttributes(vis_att);

	G4double sPosZ = fTargetChamberZOffset+fChamberLength/2; 
	G4ThreeVector sTranslate(0, 0, sPosZ);
	new G4PVPlacement(0, sTranslate, TrificPipe, "fTrificPipePhysical", ExpHallLog, false, 0);


	material = G4Material::GetMaterial(fWindowFrameMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fWindowFrameMaterial<< " not found, cannot build the Trific detector! "<<G4endl;
		return;
	}


	G4Tubs* Windowipe = new G4Tubs("Pipe",fWindowInnerD/2, fChamberOuterD/2, fWindowChamberLength/2, 0, 360.*CLHEP::deg);

	G4LogicalVolume* WindowPipe = new G4LogicalVolume(Windowipe, material, "WindowPipe", 0, 0, 0);

	sPosZ = fTargetChamberZOffset-fWindowChamberLength/2; 
	G4ThreeVector WPPosZ(0, 0, sPosZ);
	new G4PVPlacement(0, WPPosZ, WindowPipe, "fTrificWindowPipe", ExpHallLog, false, 0);	

}

void DetectionSystemTrific::BuildGasVolume(G4LogicalVolume* ExpHallLog){

	// Define the material, return error if not found
	G4Material* material = G4Material::GetMaterial(fGasMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fGasMaterial<< " not found, cannot build the Trific detector! "<<G4endl;
		return;
	}

	G4Tubs* Pipe = new G4Tubs("Pipe",0, (fChamberInnerD/2)-0.05*mm, fChamberLength/2, 0, 360.*CLHEP::deg);
	G4LogicalVolume* TrificGasVol = new G4LogicalVolume(Pipe, material, "TrificGasVol", 0, 0, 0);
	TrificGasVol->SetVisAttributes(G4VisAttributes::Invisible); 

	BuildPCB();
	PlacePCBs(TrificGasVol);	
	BuildPlaceWindow(TrificGasVol);


	G4double sPosZ = fTargetChamberZOffset+fChamberLength/2; 
	G4ThreeVector sTranslate(0, 0, sPosZ);	
	new G4PVPlacement(0, sTranslate, TrificGasVol, "fTrificGasVol", ExpHallLog, false, 0);
}



void DetectionSystemTrific::BuildPlaceWindow(G4LogicalVolume* TrificGasVol){

	//
	// First we make the part of the window frame that is inside the gas volume
	//

	// Define the material, return error if not found
	G4Material* material = G4Material::GetMaterial(fWindowFrameMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fWindowFrameMaterial<< " not found, cannot build the Trific detector! "<<G4endl;
		return;
	}	

	G4double winframhalfrise =(fWindowOuterD*tan(fGridAngle))/2;
	G4double winframhalflength =(winframhalfrise+fChamberWindowZ)/2;
	G4double cutheight =(fWindowOuterD/2)/cos(fGridAngle)+winframhalfrise*tan(fGridAngle);

	G4Tubs* CellTube = new G4Tubs("CellTube",fWindowInnerD/2,fWindowOuterD/2, winframhalflength, 0, 360.*deg);
	G4Box* CutBox = new G4Box("CutBox", 0.1*mm+fWindowOuterD/2, cutheight, winframhalfrise);

	G4double cutZ =winframhalfrise/cos(fGridAngle)+fChamberWindowZ-winframhalflength;

	G4RotationMatrix* rotate = new G4RotationMatrix;
	rotate->rotateX(fGridAngle);

	G4VSolid* WindowFrame = new G4SubtractionSolid("WindowFrame", CellTube, CutBox,rotate, G4ThreeVector(0,0,cutZ));

	G4LogicalVolume* fWindowFrame = new G4LogicalVolume(WindowFrame, material, "fWindowFrame", 0, 0, 0);

	G4double sPosZ = winframhalflength-fChamberLength/2; 
	G4ThreeVector sTranslate(0, 0, sPosZ);	
	new G4PVPlacement(0, sTranslate, fWindowFrame, "fTrificGasVol", TrificGasVol, false, 0);

	//
	// Next we add vacuum to the small bit of pipe to displace the gas in that region of the mother volume TrificGasVol
	//	

	material = G4Material::GetMaterial("Vacuum");
	if(!material) {
		G4cout<<" ----> Material Vacuum not found, cannot build the Trific detector! "<<G4endl;
		return;
	}	

	G4Tubs* VacTube = new G4Tubs("VacTube",0,(fWindowInnerD/2)-0.1*mm, winframhalflength, 0, 360.*deg);
	G4VSolid* VacSolid = new G4SubtractionSolid("VacSolid", VacTube, CutBox,rotate, G4ThreeVector(0,0,cutZ));
	G4LogicalVolume* VacLog = new G4LogicalVolume(VacSolid, material, "VacLog", 0, 0, 0);
	VacLog->SetVisAttributes(G4VisAttributes::Invisible); 
	new G4PVPlacement(0, sTranslate, VacLog, "fVacLog", TrificGasVol, false, 0);

	//
	// Next we add the window itself
	//	

	material = G4Material::GetMaterial(fWindowMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fWindowMaterial<< " not found, cannot build the Trific detector! "<<G4endl;
		return;
	}	

	G4Tubs* WindowTube = new G4Tubs("CellTube",0,fWindowOuterD/2, winframhalfrise, 0, 360.*deg);
	cutZ =winframhalfrise/cos(fGridAngle)+fWindowTickness/2;
	G4VSolid* WindowA = new G4SubtractionSolid("WindowA", WindowTube, CutBox,rotate, G4ThreeVector(0,0,cutZ));
	G4VSolid* WindowB = new G4SubtractionSolid("WindowB", WindowA, CutBox,rotate, G4ThreeVector(0,0,-cutZ));
	G4LogicalVolume* Window = new G4LogicalVolume(WindowB, material, "Window", 0, 0, 0);
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(1.0,0.5,0.0));
	vis_att->SetVisibility(true);  
	Window->SetVisAttributes(vis_att);

	sPosZ = fChamberWindowZ+fWindowTickness/2-fChamberLength/2; 
	G4ThreeVector sWindTrans(0, 0, sPosZ);	
	new G4PVPlacement(0, sWindTrans, Window, "fWindow", TrificGasVol, false, 0);

	//
	// Next we add the window aluminised layer
	//	

	material = G4Material::GetMaterial(fWindowSurfaceMaterial);
	if(!material) {
		G4cout<<" ----> Material "<<fWindowMaterial<< " not found, cannot build the Trific detector! "<<G4endl;
		return;
	}	

	cutZ =winframhalfrise/cos(fGridAngle)+fWindowCoatingThickness/2;
	G4VSolid* WindowCoverA = new G4SubtractionSolid("WindowCoverA", WindowTube, CutBox,rotate, G4ThreeVector(0,0,cutZ));
	G4VSolid* WindowCoverB = new G4SubtractionSolid("WindowCoverB", WindowCoverA, CutBox,rotate, G4ThreeVector(0,0,-cutZ));
	G4LogicalVolume* WindowCo = new G4LogicalVolume(WindowCoverB, material, "WindowCo", 0, 0, 0);
	vis_att = new G4VisAttributes(G4Colour(1.0,0.5,1.0));
	vis_att->SetVisibility(true);  
	WindowCo->SetVisAttributes(vis_att);

	sPosZ = fChamberWindowZ+fWindowTickness+fWindowCoatingThickness/2-fChamberLength/2; 
	G4ThreeVector sWindCoTrans(0, 0, sPosZ);	
	new G4PVPlacement(0, sWindCoTrans, WindowCo, "fWindowcov", TrificGasVol, false, 0);


	//
	// Finally we add the window washer
	//	

	material = G4Material::GetMaterial(fWindowFrameMaterial);

	cutZ =winframhalfrise/cos(fGridAngle)+fWindowWasherZ/2;
	G4VSolid* WindowWasherA = new G4SubtractionSolid("WindowWasherA", CellTube, CutBox,rotate, G4ThreeVector(0,0,cutZ));
	G4VSolid* WindowWasherB = new G4SubtractionSolid("WindowWasherB", WindowWasherA, CutBox,rotate, G4ThreeVector(0,0,-cutZ));
	G4LogicalVolume* WindowWasher = new G4LogicalVolume(WindowWasherB, material, "WindowWasher", 0, 0, 0);

	sPosZ = fChamberWindowZ+fWindowTickness+fWindowCoatingThickness+fWindowWasherZ/2-fChamberLength/2; 
	G4ThreeVector sWindWashTrans(0, 0, sPosZ);	
	new G4PVPlacement(0, sWindWashTrans, WindowWasher, "fWindowwash", TrificGasVol, false, 0);		


}
