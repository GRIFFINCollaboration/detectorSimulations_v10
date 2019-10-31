#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

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

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh" //individual step-size per logical volume/region

#include "ApparatusLayeredTarget.hh"

ApparatusLayeredTarget::ApparatusLayeredTarget(G4double beaminput) { 
	fBeamPos = beaminput*mm;
	fTargetRadius = 6.0*mm;

	G4double fMaxStep=10.*CLHEP::um; 
	fStepLimit = new G4UserLimits(fMaxStep);
}

ApparatusLayeredTarget::~ApparatusLayeredTarget()
{ 
	// LogicalVolumes 
	for(unsigned int i=0;i<fTargetLayerLog.size();i++){
		delete fTargetLayerLog[i];
	}
	delete fStepLimit;
}


G4int ApparatusLayeredTarget::BuildTargetLayer(G4String LayerMaterial, G4double LayerArealDensity)
{
	//Get Materials 
	G4Material* material = G4Material::GetMaterial(LayerMaterial);
	if(!material) {
		G4cout<<" ----> Target Material "<<LayerMaterial<<" not found, cannot build!"<<G4endl;
		return 0;
	} 

	G4double LayerMaterialDensity=material->GetDensity();
	G4cout<<"Target LayerMaterialDensity: "<<LayerMaterialDensity/(g/cm3) <<" g/cm3"<< G4endl;
	G4double LayerThickness = ((LayerArealDensity*(mg/cm2))/LayerMaterialDensity);
	G4cout<<"Target Thickness (um): "<<LayerThickness/(um)<<G4endl;

	// Build assembly volume
	G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.8+fTargetLayerLog.size()*0.1,0.0,0.2));
	vis_att->SetVisibility(true);  

	std::stringstream ss;
	ss<<"TargetLayer"<<fTargetLayerLog.size();
	G4Tubs* target = new G4Tubs(ss.str(), 0.0, fTargetRadius, std::abs(LayerThickness)/2.,0.0, 360*deg);

	//logical volume
	std::stringstream SS;
	SS<<"TargetLayerLog"<<fTargetLayerLog.size();
	G4LogicalVolume* LayerLog = new G4LogicalVolume(target, material,SS.str(), 0, 0, 0);
	LayerLog->SetVisAttributes(vis_att);
	fTargetLayerLog.push_back(LayerLog);
	fTargetLayerThick.push_back(LayerThickness);

	LayerLog->SetUserLimits(fStepLimit);
	return 1;
}

void ApparatusLayeredTarget::PlaceTarget(G4LogicalVolume* oexpHallLog)
{
	if(fTargetLayerLog.size()<=fTargetLayerPhys.size())return;
	G4cout<<"Adding Target Layer "<<G4endl;

	G4double placement =LayerStart(fTargetLayerPhys.size()) - fTargetLayerThick[fTargetLayerPhys.size()]/2.;
	G4cout<<"Placement target : "<<placement<<" "<<placement/mm<<G4endl;

	G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
	G4RotationMatrix* sRotate = new G4RotationMatrix;

	std::stringstream ss;
	ss<<"TargetLayer"<<fTargetLayerLog.size();
	G4PVPlacement* TargetPhys = new G4PVPlacement(sRotate,                       //no rotation
			move,                    //at position
			fTargetLayerLog[fTargetLayerPhys.size()],             //its logical volume
			ss.str(),                //its name
			oexpHallLog,                //its mother  volume
			false,                   //no boolean operation
			false,                       //copy number
			0);          //overlaps checking

	fTargetLayerPhys.push_back(TargetPhys);
}

G4double ApparatusLayeredTarget::LayerStart(int layeri){
	if(!fTargetLayerThick.size())return 0;
	G4double ret = fBeamPos;

	for(int i=0;i<layeri;i++){
		if((unsigned)i<fTargetLayerThick.size())ret-=fTargetLayerThick[i];
	}
	return ret;
}
