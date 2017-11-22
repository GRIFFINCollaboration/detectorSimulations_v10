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

#include "ApparatusSpiceTarget.hh"

ApparatusSpiceTarget::ApparatusSpiceTarget(G4double beaminput) { 
  
  fBeamPos = beaminput*mm;
  //G4cout << "\n\n\n\n\n\ntarget z in spice target: " << fBeamPos << "\n\n\n\n\n\n"<< G4endl; 
  this->fTargetMaterial = "Gold";
  this->fTargetBackerMaterial = "Gold";
  this->fTargetProtectorMaterial = "Gold";//initialising - may remove
  this->fTargetRadius = 12.5*mm;
  this->fTargetBackerSurfaceDensity=0.;
  
  fBracketMaterial = "G4_Al";
  fBracketSideLength = 19.77*mm;
  fBracketBackLength = 30.77*mm;
  fBracketBackTopThickness = 3.63*mm;//x/y-axis
  fBracketSideTopThickness = 4.34*mm;
  fBracketDepth = 4.08*mm;//z-axis
  
  fHolderMaterial = "G4_Al";
  fHolderOuterRadii = 26.*mm;
  fHolderInnerRadii = 25*mm;
  fHolderThickness = 0.5*mm;
  
  fMaxStep= 10.*CLHEP::um;
  
  fStepLimit = new G4UserLimits(fMaxStep);
}

ApparatusSpiceTarget::~ApparatusSpiceTarget()
{ 
  // LogicalVolumes 
    delete fTargetSpiceLog;
    delete fTargetBackerSpiceLog;
    delete fTargetProtectorSpiceLog;
    delete fTargetBracketLog;
    delete fTargetHolderLog;
    delete fTargetPhys;
    delete fTargetBackerPhys;
    delete fTargetProtectorPhys;
    delete fTargetBracketPhys;
    delete fTargetHolderPhys;
    delete fSourceLog;//May be temp
    delete fSourcePhys;
    
    delete fStepLimit;
}

G4int ApparatusSpiceTarget::BuildTarget(G4String input_material, G4double input_surface_density, G4double input_density)
{
    //G4cout << "\n\n\n\n\n\ntarget z in spice target build: " << fBeamPos << "\n\n\n\n\n\n"<< G4endl; 
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.8,0.0,0.2));
    vis_att->SetVisibility(true);  

    //Get Materials 
    this->fTargetMaterial = input_material;
    this->fTargetMaterialDensity = input_density*(g/cm3);//need units
    G4cout << "Post-units: " << fTargetMaterialDensity << G4endl;
    this->fTargetSurfaceDensity = input_surface_density*(mg/cm2);//need units
    G4cout << "Post-units: " << fTargetSurfaceDensity << G4endl;
    this->fTargetThickness = (fTargetSurfaceDensity/fTargetMaterialDensity);
    G4cout << "Target Thickness: " << fTargetThickness << G4endl;

    // Build assembly volume
    //G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
    //this->assembly = myAssembly;
   
    G4Material* material = G4Material::GetMaterial(this->fTargetMaterial);
   // G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    if( !material ) {
      G4cout << " ----> Target Material " << this->fTargetMaterial << " not found, cannot build!" << G4endl;
      return 0;
    } 
    
    G4Tubs* target = new G4Tubs("Spice_Target", 0.0, fTargetRadius, fTargetThickness/2.,
                0.0, 360*deg);
    
    //logical volume
      fTargetSpiceLog = new G4LogicalVolume(target, material, "target_spice_log", 0, 0, 0);
      fTargetSpiceLog->SetVisAttributes(vis_att);
      return 1;
}

G4int ApparatusSpiceTarget::BuildBacker(G4String input_material, G4double input_surface_density, G4double input_density)
{
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.2,0.8));
    vis_att->SetVisibility(true);  
    
    this->fTargetBackerMaterial = input_material;
    this->fTargetBackerMaterialDensity = input_density*(g/cm3);//need units
    this->fTargetBackerSurfaceDensity = input_surface_density*(mg/cm2);//need units
    this->fTargetBackerThickness = (fTargetBackerSurfaceDensity/fTargetBackerMaterialDensity);
    
    G4Material* material = G4Material::GetMaterial(this->fTargetBackerMaterial);
    if( !material ) {
      G4cout << " ----> Target Backer Material " << this->fTargetBackerMaterial << " not found, cannot build!" << G4endl;
      return 0;
    }
    
    G4Tubs* target_Backer = new G4Tubs("Spice_Target_Backer", 0.0, fTargetRadius, fTargetBackerThickness/2.,
                0.0, 360*deg);
    
    //logical volume
     fTargetBackerSpiceLog = new G4LogicalVolume(target_Backer, material, "target_Backer_spice_log", 0, 0, 0);
     fTargetBackerSpiceLog->SetVisAttributes(vis_att);
     return 1;
}

G4int ApparatusSpiceTarget::BuildProtector(G4String input_material, G4double input_surface_density, G4double input_density)
{
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,0.8,0.2));
    vis_att->SetVisibility(true);  
    
    this->fTargetProtectorMaterial = input_material;
    this->fTargetProtectorMaterialDensity = input_density*(g/cm3);//need units  
    this->fTargetProtectorSurfaceDensity = input_surface_density*(mg/cm2);//need units
    this->fTargetProtectorThickness = (fTargetProtectorSurfaceDensity/fTargetProtectorMaterialDensity);
    
    G4Material* material = G4Material::GetMaterial(this->fTargetProtectorMaterial);
    if( !material ) {
      G4cout << " ----> Target Protector Material " << this-> fTargetProtectorMaterial<< " not found, cannot build!" << G4endl;
      return 0;
    }  

    G4Tubs* target_Protector = new G4Tubs("Spice_Target_Protector", 0.0, fTargetRadius, fTargetProtectorThickness/2.,
                0.0, 360*deg);
    
    //logical volume
    fTargetProtectorSpiceLog = new G4LogicalVolume(target_Protector, material, "target_Protector_spice_log", 0, 0, 0);
    fTargetProtectorSpiceLog->SetVisAttributes(vis_att);
    return 1;
}


G4int ApparatusSpiceTarget::BuildBracket(){
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(132./(404.),135./404.,137./404.));//bracket
    vis_att->SetVisibility(true); 

    
    G4Material* material = G4Material::GetMaterial(this->fBracketMaterial);
    if( !material ) {
      G4cout << " ----> Target Bracket Material " << this-> fBracketMaterial<< " not found, cannot build!" << G4endl;
      return 0;
    }  
    
    G4double BracketBackLengthRemove = 22.09*mm;
    G4double BracketSideLengthRemove = 16.14*mm;
    
    G4Box* OverBracket = new G4Box("Bracket", fBracketBackLength/2., fBracketSideLength/2., fBracketDepth/2.);//Larger rectangle
    G4Box* InnerBracket = new G4Box("TakeAway", (BracketBackLengthRemove/2.)*mm,(BracketSideLengthRemove/2.)*mm, fBracketDepth/1.2);//Smaller rectangle
    
    G4RotationMatrix* sRotate = new G4RotationMatrix;//no rotation but need null-matrix
    G4ThreeVector sTrans(0., -3.64*mm, 0.);//Moves in Y to create 3-sided bracket
    
    G4SubtractionSolid* Bracket = new G4SubtractionSolid("Bracket", OverBracket, InnerBracket, sRotate, sTrans);
    
    fTargetBracketLog = new G4LogicalVolume(Bracket, material, "target_Bracket_spice_log", 0, 0, 0);
    fTargetBracketLog->SetVisAttributes(vis_att);
    
  return 1;
}

G4int ApparatusSpiceTarget::BuildHolder()
{
  
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.4,0.4,0.2));
    vis_att->SetVisibility(true);

    G4Material* material = G4Material::GetMaterial(fHolderMaterial);
   // G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    if( !material ) {
      G4cout << " ----> Target Material " << fHolderMaterial << " not found, cannot build!" << G4endl;
      return 0;
    } 
    fHolderThickness =  fTargetProtectorThickness + fTargetBackerThickness +fTargetThickness+1.*mm; //all data in by this point - as inputs before builds (called upon material choice)
    G4Tubs* holder = new G4Tubs("Spice_Target_Holder", fHolderInnerRadii/2., fHolderOuterRadii/2., fHolderThickness/2.,
                0.0, 360*deg);
    
    //logical volume
      fTargetHolderLog = new G4LogicalVolume(holder, material, "target_spice_holder_log", 0, 0, 0);
      fTargetHolderLog->SetVisAttributes(vis_att);
      return 1;
}

G4int ApparatusSpiceTarget::BuildSource()//may be used in the future
{
   /*G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.338,0.331,0.331));
   vis_att->SetVisibility(true);
   G4Material* material = G4Material::GetMaterial("Barium");
   if( !material ) {
      G4cout << " ----> Source Material " << "Barium" << " not found, cannot build!" << G4endl;
      return 0;
   }
   G4double SourceThickness = 0.27624309392*CLHEP::um*2.;
   G4Tubs* source = new G4Tubs("Spice_Source_Block", 0., 4.*CLHEP::mm, SourceThickness,
                0.0, 360*deg);
   //logical volume
   fSourceLog = new G4LogicalVolume(source, material, "target_spice_source_log", 0, 0, 0);
   fSourceLog->SetVisAttributes(vis_att);*/
   return 1;
}

void ApparatusSpiceTarget::PlaceTarget(G4LogicalVolume* oexpHallLog)
{
  G4cout << "Spice Target implemented " << G4endl;
  G4cout << "targ-thick: " << fTargetThickness << " BACKER IN TARG: " <<this->fTargetBackerThickness << G4endl;
  
  G4double placement = fBeamPos*mm - fTargetBackerThickness - fTargetThickness/2.;
  if(fBeamPos < -3.9*mm) placement = placement - 0.5*mm;//0.5mm=half thickness of target wheel where 0,0,0 is placed
  G4cout << "Placement target : " << placement << " " << placement*mm << G4endl;
  //NEED TO LINK BEAM POS TO HERE for move = three-vector for mosition
  //fBeamPos*mm-(target_Backer_thickness)*cm
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
   
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  fTargetPhys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    fTargetSpiceLog,             //its logical volume
                    "Target_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
  
  fTargetSpiceLog->SetUserLimits(fStepLimit);
 /* BuildSource();
  PlaceSource(oexpHallLog);//(0.27624309392*2.)*CLHEP::um */
}

void ApparatusSpiceTarget::PlaceTargetBacker(G4LogicalVolume* oexpHallLog)
{
  G4cout << "Spice Target backer implemented " << G4endl;
  G4cout << "Back-Thick: " << fTargetBackerThickness<< G4endl;

  G4double placement = fBeamPos*mm - (this->fTargetBackerThickness/2.);
  if(fBeamPos < -3.9*mm) placement = placement - 0.5*mm;//0.5mm=half thickness of target wheel where 0,0,0 is placed
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  fTargetBackerPhys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    fTargetBackerSpiceLog,             //its logical volume
                    "Target_backer_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
  
    fTargetBackerSpiceLog->SetUserLimits(fStepLimit);
}

void ApparatusSpiceTarget::PlaceTargetProtector(G4LogicalVolume* oexpHallLog)
{
  G4cout << "Spice Target protector implemented " << G4endl;
  G4cout << " Pro-thick: " << fTargetProtectorThickness << G4endl;
  
  G4double placement = fBeamPos*mm - (this->fTargetBackerThickness) - (this->fTargetThickness) 
		    - (this->fTargetProtectorThickness/2.)-1.*mm;//1mm for bismuth tests
  if(fBeamPos < -3.9*mm) placement = placement - 0.5*mm;//0.5mm=half thickness of target wheel where 0,0,0 is placed
  G4cout << "Placement pro : " << placement << " " << placement*mm << " targetthick " << (this->fTargetThickness)<<" "<<G4endl;
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  fTargetProtectorPhys = new G4PVPlacement(sRotate,                       //rotation
                    move,                    //at position
                    fTargetProtectorSpiceLog,             //its logical volume
                    "Target_Protector_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
  
  fTargetProtectorSpiceLog->SetUserLimits(fStepLimit);
}

void ApparatusSpiceTarget::PlaceBracket(G4LogicalVolume* oexpHallLog){
  G4cout << "Spice Target bracket implemented " << G4endl;
  
  G4double placement = fBeamPos*mm - (this->fTargetBackerThickness) - (this->fTargetThickness)
	      - (this->fTargetProtectorThickness) - fBracketBackTopThickness/2. - 1.5*mm;//1.5mm for bismuth tests
   G4cout << "Spice Target bracket implemented " << G4endl;
  //if parts of target not used, values are initialised to 0 thickness
   
  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//position of bracket
  G4RotationMatrix* sRotate = new G4RotationMatrix;//rotation of bracket
  sRotate->rotateZ(277.*CLHEP::deg);//Measured rotation
  
  fTargetBracketPhys = new G4PVPlacement(sRotate,                       //rotation
                    move,                    //at position
                    fTargetBracketLog,             //its logical volume
                    "Target_Bracket_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
  
  fTargetBracketLog->SetUserLimits(fStepLimit);
  
  G4cout << "Spice Target bracket implemented " << G4endl;
  
}
void ApparatusSpiceTarget::PlaceHolder(G4LogicalVolume* oexpHallLog){
  
  G4cout << "Spice Target holder implemented " << G4endl;
  G4double placement = fBeamPos*mm - fHolderThickness/2.;//0.5mm=half thickness of target wheel where 0,0,0 is placed
  if(fBeamPos < -3.9*mm) placement -= 0.5*mm;
  G4cout << "Placement holder: " << placement << " " << placement*mm << G4endl;

  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
   
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  fTargetHolderPhys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    fTargetHolderLog,             //its logical volume
                    "Target_Holder_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
  
  fTargetHolderLog->SetUserLimits(fStepLimit);
}

void ApparatusSpiceTarget::PlaceSource(){
  /*G4cout << "Spice source implemented " << G4endl;
  
  G4double placement = fBeamPos*mm - fTargetBackerThickness -(0.27624309392)*CLHEP::um ;//0.5mm=half thickness of target wheel where 0,0,0 is placed
  if(fBeamPos < -3.9*mm) placement -= 0.5*mm;
  G4cout << "Placement source: " << placement << " " << placement*mm << G4endl;

  G4ThreeVector move = G4ThreeVector(0.,0.,placement);//0 for now until read in 
   
  G4RotationMatrix* sRotate = new G4RotationMatrix;
  sRotate->rotateZ(90.*CLHEP::deg);
  fSourcePhys = new G4PVPlacement(sRotate,                       //no rotation
                    move,                    //at position
                    fSourceLog,             //its logical volume
                    "Target_source_SPICE",                //its name
                    oexpHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    false,                       //copy number
                    0);          //overlaps checking
  
  fSourceLog->SetUserLimits(fStepLimit);*/
}