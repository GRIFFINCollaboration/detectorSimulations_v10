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
/// \file analysis/shared/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
// $Id: EventAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "Randomize.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* run)
    :G4UserEventAction(),
      fRunAct(run),
      fPrintModulo(0)
{
    fNumberOfHits = 0;
    fNumberOfSteps = 0;
    fPrintModulo = 1000;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() { 
  G4cout << "1 dep: " << MultiplicityArray[0] << "| 2 dep: " << MultiplicityArray[1] << "| 3 dep: " << MultiplicityArray[2] << 
  "| 4 dep: " << MultiplicityArray[3] << "| 5 dep: " << MultiplicityArray[4] << G4endl;
  
  G4double ME = MultiplicityArray[1] + MultiplicityArray[2] + MultiplicityArray[3] + MultiplicityArray[4]; //multiple events
  
  G4double OS = ME/(MultiplicityArray[0] + ME); //multiple/all events
  
  G4cout << "Old style multiplicity: " << OS*100.0 << G4endl;//% figure for multiplicity
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt) {  
    fEvtNb = evt->GetEventID();
    if (fEvtNb%fPrintModulo == 0)
        //    G4cout << "\n---> Begin of event: " << fEvtNb << G4endl;
        printf( " ---> Ev.# %5d\r", fEvtNb); //lists up-to-date event number //% offers formatted output following printf
    G4cout.flush();

    ClearVariables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*) {
    BeamInputEnergy = HistoManager::Instance().BeamEnergy;
    FillParticleType() ; // was uncommented - otherwise not filled
    FillGriffinCryst() ;
    Fill8piCryst() ;
    FillLaBrCryst() ;
    FillAncillaryBgo() ;
    FillSceptar() ;
    FillGridCell() ;
    FillSpice() ; ///19/6
    FillPacesCryst() ; ///20/6	
    //G4cout << "fNumberOfHits = " << fNumberOfHits << G4endl;
    for (G4int i = 0 ; i < fNumberOfHits; i++) {
		HistoManager::Instance().FillHitNtuple(fHitTrackerI[0][i], fHitTrackerI[1][i], fHitTrackerI[2][i], fHitTrackerI[3][i],  fHitTrackerI[4][i], fHitTrackerI[5][i], fHitTrackerI[6][i], fHitTrackerI[7][i], fHitTrackerI[8][i], fHitTrackerD[0][i]/keV, fHitTrackerD[1][i]/mm, fHitTrackerD[2][i]/mm, fHitTrackerD[3][i]/mm, fHitTrackerD[4][i]/second, fHitTrackerI[9][i]);
    }
    for (G4int i = 0 ; i < fNumberOfSteps; i++) {
		HistoManager::Instance().FillStepNtuple(fStepTrackerI[0][i], fStepTrackerI[1][i], fStepTrackerI[2][i], fStepTrackerI[3][i],  fStepTrackerI[4][i], fStepTrackerI[5][i], fStepTrackerI[6][i], fStepTrackerI[7][i], fStepTrackerI[8][i], fStepTrackerD[0][i]/keV, fStepTrackerD[1][i]/mm, fStepTrackerD[2][i]/mm, fStepTrackerD[3][i]/mm, fStepTrackerD[4][i]/second, fStepTrackerI[9][i]);
    }

    ClearVariables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::ClearVariables() {
    if(HistoManager::Instance().GetStepTrackerBool()) {
        fStepIndex = 0;
        fNumberOfSteps = 0;
        for (G4int i = 0 ; i < MAXSTEPS; i++) {
            for (G4int j = 0 ; j < NUMSTEPVARS; j++) {
                fStepTrackerI[j][i] = 0;
                fStepTrackerD[j][i] = 0.0;
            }
        }
    }

    if(HistoManager::Instance().GetHitTrackerBool()) {
        fHitIndex = 0;
        fNumberOfHits = 0;
        fPTrackID = -1;
        fPParentID = -1;

        for (G4int i = 0 ; i < MAXHITS; i++) {
            fPHitMnemonic[i] = "XXX00XX00X";
            for (G4int j = 0 ; j < NUMSTEPVARS; j++) {
                fHitTrackerI[j][i] = 0;
                fHitTrackerD[j][i] = 0.0;
            }
        }
    }

    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++) {
        fParticleTypes[i]   = 0;
    }
    for (G4int iii = 0; iii < MAXNUMDETPACES; iii++){		//4/7[MAXNUMDETPACES]
	fPacesCrystEnergyDet[iii]	= 0 ;///20/6
	fPacesCrystTrackDet[iii]	= 0 ;//20/6
    }
    for (G4int i = 0; i < MAXNUMDETSPICE; i++){//30/6 //above number needed, doubly sure of cleared variables
      for (G4int j=0; j<12; j++){
	fSpiceEnergyDet[i][j]	= 0 ; ///19/6
  	fSpiceTrackDet[i][j]	= 0 ; ///19/6
    }}
    for (G4int i = 0 ; i < MAXNUMDET; i++) {
        fEightPiCrystEnergyDet[i]  = 0 ;
        fEightPiCrystTrackDet[i]   = 0 ;
		  
        fLaBrCrystEnergyDet[i]  = 0 ;
        fLaBrCrystTrackDet[i]   = 0 ;
		  
        fAncillaryBgoEnergyDet[i]  = 0 ;
        fAncillaryBgoTrackDet[i]   = 0 ;
	
        fSceptarEnergyDet[i]  = 0 ;
        fSceptarTrackDet[i]   = 0 ;
		  
        fGridCellElectronEKinDet[i]  = 0 ;
        fGridCellElectronTrackDet[i]   = 0 ;
        fGridCellGammaEKinDet[i]  = 0 ;
        fGridCellGammaTrackDet[i]   = 0 ;
        fGridCellNeutronEKinDet[i]  = 0 ;
        fGridCellNeutronTrackDet[i]   = 0 ;
    }
	
    for (G4int i = 0 ; i < MAXNUMDETGRIFFIN; i++) {
        for (G4int j = 0 ; j < MAXNUMCRYGRIFFIN; j++) {
            fGriffinCrystEnergyDet[i][j]                 = 0;
            fGriffinCrystTrackDet[i][j]                  = 0;
            fGriffinSuppressorBackEnergyDet[i][j]        = 0;
            fGriffinSuppressorBackTrackDet[i][j]         = 0;
            fGriffinSuppressorLeftExtensionEnergyDet[i][j]   = 0;
            fGriffinSuppressorLeftExtensionTrackDet[i][j]    = 0;
            fGriffinSuppressorLeftSideEnergyDet[i][j]        = 0;
            fGriffinSuppressorLeftSideTrackDet[i][j]         = 0;
            fGriffinSuppressorRightExtensionEnergyDet[i][j]   = 0;
            fGriffinSuppressorRightExtensionTrackDet[i][j]    = 0;
            fGriffinSuppressorRightSideEnergyDet[i][j]        = 0;
            fGriffinSuppressorRightSideTrackDet[i][j]         = 0;
        }
    }  // NOTE: Clear the variables from the new Fill___Cryst functions as last EndofEvent action
}


void EventAction::FillParticleType() {///this works
    G4int numParticleTypes = 0;
    for(G4int i = 0 ; i < NUMPARTICLETYPES; i++) {
        if(fParticleTypes[i] != 0) { // if particle type 'i' has non-zero counts
            for (G4int j = 0 ; j< fParticleTypes[i]; j++) { // loop over the number of time we saw it
                //G4cout << "fParticleTypes[" << i << "] = " << fParticleTypes[i] << G4endl; //super annoying, not optimised
                HistoManager::Instance().FillHisto(kAstatsParticleTypeInEachStep, i);
            }
        }
    }

    // Fill the number of particle types in the event
    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++)
    {
        if (fParticleTypes[i] != 0)
            numParticleTypes++;
    }
    HistoManager::Instance().FillHisto(kAstatsParticleTypeInEachEvent, numParticleTypes);
}

void EventAction::FillGriffinCryst()
{
    G4double  energySum = 0;
    G4double  energySumDet = 0;
    G4bool suppressorBackFired[MAXNUMDETGRIFFIN] = {0};
    G4bool suppressorExtensionFired[MAXNUMDETGRIFFIN] = {0};
    G4bool suppressorSideFired[MAXNUMDETGRIFFIN] = {0};
    G4bool suppressorFired = false;

    // Fill Griffin Histos
    for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
        energySumDet = 0;
        // Find if any suppressors were fired
        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            if(fGriffinSuppressorBackEnergyDet[i][j] > MINENERGYTHRES) {
                suppressorBackFired[i] = true;
            }
            if(fGriffinSuppressorLeftExtensionEnergyDet[i][j] > MINENERGYTHRES) {
                suppressorExtensionFired[i] = true;
            }
            if(fGriffinSuppressorRightExtensionEnergyDet[i][j] > MINENERGYTHRES) {
                suppressorExtensionFired[i] = true;
            }
            if(fGriffinSuppressorLeftSideEnergyDet[i][j] > MINENERGYTHRES) {
                suppressorSideFired[i] = true;
            }
            if(fGriffinSuppressorRightSideEnergyDet[i][j] > MINENERGYTHRES) {
                suppressorSideFired[i] = true;
            }
            if(!suppressorFired && ( suppressorBackFired[i] || suppressorExtensionFired[i] || suppressorSideFired[i] ) ) {
                suppressorFired = true;
            }
        }

        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            if(fGriffinCrystEnergyDet[i][j] > MINENERGYTHRES) {
                // fill energies in each crystal
                if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto((kGriffinCrystalUnsupEdepDet0Cry0+(MAXNUMDETGRIFFIN*j))+i, fGriffinCrystEnergyDet[i][j]);
                if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kGriffinCrystalUnsupEdepCry, fGriffinCrystEnergyDet[i][j]);
                if(!suppressorBackFired[i] && !suppressorExtensionFired[i] && !suppressorSideFired[i]) { // Suppressor fired?
                    if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto((kGriffinCrystalSupEdepDet0Cry0+(MAXNUMDETGRIFFIN*j))+i, fGriffinCrystEnergyDet[i][j]);
                    if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto(kGriffinCrystalSupEdepCry, fGriffinCrystEnergyDet[i][j]);
                }
                energySumDet += fGriffinCrystEnergyDet[i][j];
            }
        }
        if(energySumDet > MINENERGYTHRES) {
            // fill energies in each detector
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kGriffinCrystalUnsupEdepDet0+i, energySumDet);
            // fill standard energy and track spectra
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kGriffinCrystalUnsupEdep, energySumDet);
            if(!suppressorBackFired[i] && !suppressorExtensionFired[i] && !suppressorSideFired[i]) {
                // fill energies in each detector
                if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto(kGriffinCrystalSupEdepDet0+i, energySumDet);
                // fill standard energy and track spectra
                if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto(kGriffinCrystalSupEdep, energySumDet);
            }
        }
        energySum += energySumDet;
    }

    if(energySum > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kGriffinCrystalUnsupEdepSum, energySum);
        if(!suppressorFired) {
            if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto(kGriffinCrystalSupEdepSum, energySum);
        }
    }
}

void EventAction::Fill8piCryst() {
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(fEightPiCrystEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kEightpiCrystalEdep, fEightPiCrystEnergyDet[j]);
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kEightpiCrystalEdepDet0+j, fEightPiCrystEnergyDet[j]);
            energySumDet += fEightPiCrystEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kEightpiCrystalEdepSum, energySumDet);
    }
}


void EventAction::FillLaBrCryst() {
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(fLaBrCrystEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kLabrCrystalEdep, fLaBrCrystEnergyDet[j]);
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kLabrCrystalEdepDet0+j, fLaBrCrystEnergyDet[j]);
            energySumDet += fLaBrCrystEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kLabrCrystalEdepSum, energySumDet);
    }
}

void EventAction::FillAncillaryBgo() {
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(fAncillaryBgoEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kAncillaryBgoCrystalEdep, fAncillaryBgoEnergyDet[j]);
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kAncillaryBgoCrystalEdepDet0+j, fAncillaryBgoEnergyDet[j]);
            energySumDet += fAncillaryBgoEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kAncillaryBgoCrystalEdepSum, energySumDet);
    }
}

void EventAction::FillSceptar() {
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(fSceptarEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kSceptarEdep, fSceptarEnergyDet[j]);
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kSceptarEdepDet0+j, fSceptarEnergyDet[j]);
            energySumDet += fSceptarEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(kSceptarEdepSum, energySumDet);
    }

}

void EventAction::FillGridCell() {
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(WRITEEKINHISTOS && fGridCellElectronEKinDet[j] > MINENERGYTHRES)  HistoManager::Instance().FillHisto(kGridcellElectronEkinDet0+j, fGridCellElectronEKinDet[j]);
        if(WRITEEKINHISTOS && fGridCellGammaEKinDet[j] > MINENERGYTHRES)     HistoManager::Instance().FillHisto(kGridcellGammaEkinDet0+j, fGridCellGammaEKinDet[j]);
        if(WRITEEKINHISTOS && fGridCellNeutronEKinDet[j] > MINENERGYTHRES)   HistoManager::Instance().FillHisto(kGridcellNeutronEkinDet0+j, fGridCellNeutronEKinDet[j]);
    }
}

void EventAction::FillSpice() {
	
	//G4cout << "size of histo array: " << sizeof(HistoManager::Instance().SpiceHistNumbers) << G4endl;
	G4double energySumDet = 0;
	fSpiceMultiplicity = 0;
	G4double Multiplicityenergy = 0.0;
	for (G4int ring=0; ring < MAXNUMDETSPICE; ring++) {
	  for (G4int seg=0; seg < 12; seg++) {
	    if(fSpiceEnergyDet[ring][seg] > MINENERGYTHRES) {
	      //fill energies in each detector
	      if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().SpiceHistNumbers[0], G4RandGauss::shoot(fSpiceEnergyDet[ring][seg],0.002));
	      //fill standard energy spectra
	      if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().SpiceHistNumbers[MAXNUMDETSPICE*ring+seg+2], G4RandGauss::shoot(fSpiceEnergyDet[ring][seg],0.002));
	      fSpiceMultiplicity += 1;
	      Multiplicityenergy += fSpiceEnergyDet[ring][seg];//dont need smeared energies
	      //add sum energies
	      /*Resolution of parameters of a HPGe detector is implemented here, only for the sum as of 8/8
	  G4double A1 = 0.95446 ;
	  G4double B1 = 0.00335 ;
	  G4double C1 = -7.4117e-7 ;
	  G4double TX1 = ((energySumDet)/keV);
	  G4double per_res1 = (sqrt(A1 + (B1*TX1) + (C1*TX1*TX1)))/2.363;
	  G4double sigma1 = (per_res1/1000);
	  energySumDet += fSpiceEnergyDet[ring][seg]; //no smearing
	  G4cout << "Sigma1 = "<<sigma1<<G4endl;*/
	      energySumDet += G4RandGauss::shoot(fSpiceEnergyDet[ring][seg],0.002);//2keV resolution for spice so far - not energy dependent
	      //summing enery deposits per event for add-back method
	    }
	  }
	}
	if(energySumDet > MINENERGYTHRES) {//after exiting loops for all rings/segs, will input energy if > threshold
	  //fill sum energies
	  if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().SpiceHistNumbers[1], energySumDet);
	  }
	  
	  if(energySumDet > (BeamInputEnergy-0.015)){//0.6 = 3 sigma 
	  if(fSpiceMultiplicity>0) HistoManager::Instance().FillHisto(HistoManager::Instance().angledistro[0], fSpiceMultiplicity);
	  switch(fSpiceMultiplicity) {
	  case 1 : MultiplicityArray[0] += 1;
             break;       
	  case 2 : MultiplicityArray[1] += 1;
	     break;
	  case 3 : MultiplicityArray[2] += 1; 
             break;       
	  case 4 : MultiplicityArray[3] += 1;
             break;
	  case 5 : MultiplicityArray[4] += 1;
             break;       
  
	  }}
	  
	//  G4cout << "1 dep: " << MultiplicityArray[0] << "| 2 dep: " << MultiplicityArray[1] << "| 3 dep: " << MultiplicityArray[2] << 
 // "| 4 dep: " << MultiplicityArray[3] << "| 5 dep: " << MultiplicityArray[4] << G4endl;
	//  G4cout << "Total edep: " << MultiplicityArray[0] + MultiplicityArray[1]*2 + MultiplicityArray[2]*3 +MultiplicityArray[3]*4 +MultiplicityArray[4]*5 <<G4endl;
  
}
//if(WRITEEDEPHISTOS && (energySumDet > MINENERGYTHRES))     G4cout << "energysumDet " <<  energySumDet << G4endl;} ///////////////////////////////////////19/6

void EventAction::FillPacesCryst() {
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDETPACES; j++) {
        if(fPacesCrystEnergyDet[j] > MINENERGYTHRES) {
	  
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().PacesHistNumbers[0], fPacesCrystEnergyDet[j]);
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().PacesHistNumbers[j+2], fPacesCrystEnergyDet[j]);
            energySumDet += fPacesCrystEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS && (energySumDet > MINENERGYTHRES))     HistoManager::Instance().FillHisto(HistoManager::Instance().PacesHistNumbers[1], energySumDet);
	
    }
}////////////////////////////////////////////20/6

void EventAction::AddHitTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ) {
    G4bool newhit = true;
    for (G4int i = 0 ; i < fNumberOfHits; i++) {
        if(fPHitMnemonic[i] == mnemonic) {
            // sum the new enery
            fHitTrackerD[0][i] = fHitTrackerD[0][i] + depEnergy;
            newhit = false;
            break;
        }
    }
    if (newhit) { // new hit
        fPHitMnemonic[fHitIndex] = mnemonic;
        fPTrackID = trackID;
        fPParentID = parentID;
        fHitTrackerI[0][fHitIndex] = eventNumber;
        fHitTrackerI[1][fHitIndex] = trackID;
        fHitTrackerI[2][fHitIndex] = parentID;
        fHitTrackerI[3][fHitIndex] = stepNumber;
        fHitTrackerI[4][fHitIndex] = particleType;
        fHitTrackerI[5][fHitIndex] = processType;
        fHitTrackerI[6][fHitIndex] = systemID;
        fHitTrackerI[7][fHitIndex] = cryNumber;
        fHitTrackerI[8][fHitIndex] = detNumber;
        fHitTrackerI[9][fHitIndex] = targetZ;
        fHitTrackerD[0][fHitIndex] = depEnergy;
        fHitTrackerD[1][fHitIndex] = posx;
        fHitTrackerD[2][fHitIndex] = posy;
        fHitTrackerD[3][fHitIndex] = posz;
        fHitTrackerD[4][fHitIndex] = time;

        fHitIndex++;
        fNumberOfHits = fHitIndex;

        if(fNumberOfHits >= MAXHITS) {
            G4cout << "ERROR! Too many hits!" << G4endl;
        }
    }
}

void EventAction::AddStepTracker(G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ) {
    G4bool newstep = true;
    if(newstep) { // new step
        fStepTrackerI[0][fStepIndex] = eventNumber;
        fStepTrackerI[1][fStepIndex] = trackID;
        fStepTrackerI[2][fStepIndex] = parentID;
        fStepTrackerI[3][fStepIndex] = stepNumber;
        fStepTrackerI[4][fStepIndex] = particleType;
        fStepTrackerI[5][fStepIndex] = processType;
        fStepTrackerI[6][fStepIndex] = systemID;
        fStepTrackerI[7][fStepIndex] = cryNumber;
        fStepTrackerI[8][fStepIndex] = detNumber;
        fStepTrackerI[9][fStepIndex] = targetZ;
        fStepTrackerD[0][fStepIndex] = depEnergy;
        fStepTrackerD[1][fStepIndex] = posx;
        fStepTrackerD[2][fStepIndex] = posy;
        fStepTrackerD[3][fStepIndex] = posz;
        fStepTrackerD[4][fStepIndex] = time;

        fStepIndex++;

        fNumberOfSteps = fStepIndex;
        if(fNumberOfSteps >= MAXSTEPS) {
            G4cout << "ERROR! Too many steps!" << G4endl;
        }
    }
}
