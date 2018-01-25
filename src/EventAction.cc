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
#include <iostream>
#include <fstream>

#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* run)
    :G4UserEventAction(),
      fRunAct(run),
      fPrintModulo(0)
{
    G4cout << "EVENT ACTION CONSTRUCTOR" << G4endl;
    fNumberOfHits = 0;
    fNumberOfSteps = 0;
    fPrintModulo = 1000;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {	
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
    BeamInputEnergy = HistoManager::Instance().fBeamEnergy;
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
    
    if(!HistoManager::Instance().Spice()){//Not used in SPICE, saves memory and speeds up simulation //'!SPICE' may be too heavy-handed
      for (G4int i = 0 ; i < fNumberOfHits; i++) {
		HistoManager::Instance().FillHitNtuple(fHitTrackerI[0][i], fHitTrackerI[1][i], fHitTrackerI[2][i], fHitTrackerI[3][i],  fHitTrackerI[4][i], fHitTrackerI[5][i], fHitTrackerI[6][i], fHitTrackerI[7][i], fHitTrackerI[8][i], fHitTrackerD[0][i]/keV, fHitTrackerD[1][i]/mm, fHitTrackerD[2][i]/mm, fHitTrackerD[3][i]/mm, fHitTrackerD[4][i]/second, fHitTrackerI[9][i]);
      }
      for (G4int i = 0 ; i < fNumberOfSteps; i++) {
		HistoManager::Instance().FillStepNtuple(fStepTrackerI[0][i], fStepTrackerI[1][i], fStepTrackerI[2][i], fStepTrackerI[3][i],  fStepTrackerI[4][i], fStepTrackerI[5][i], fStepTrackerI[6][i], fStepTrackerI[7][i], fStepTrackerI[8][i], fStepTrackerD[0][i]/keV, fStepTrackerD[1][i]/mm, fStepTrackerD[2][i]/mm, fStepTrackerD[3][i]/mm, fStepTrackerD[4][i]/second, fStepTrackerI[9][i]);
      }

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
//     G4cout << "FillSpice entered " << G4endl;
    G4double energySumDet = 0;
    G4Int fSpiceMultiplicity = 0;
    G4double SpiceEnergy,SpiceEnergyRaw;
    for (G4int ring=0; ring < MAXNUMDETSPICE; ring++) {
      for (G4int seg=0; seg < 12; seg++) {
	if(fSpiceEnergyDet[ring][seg] > 70.*CLHEP::keV&&WRITEEDEPHISTOS) {
		
		SpiceEnergyRaw = fSpiceEnergyDet[ring][seg];// SpiceRawEnergy. No resolution applied
		SpiceEnergy=SpiceEnergyRaw;
		if(HistoManager::Instance().SpiceResolution()){
			SpiceEnergy = ApplySpiceRes(SpiceEnergyRaw);// + 9.0*CLHEP::keV ;// + 2.1*CLHEP::keV;	  //fill decay histo
			HistoManager::Instance().FillHisto(HistoManager::Instance().fAngleDistro[8],SpiceEnergy*1000);
		}
		
	  //fill energies across all segments
	    HistoManager::Instance().FillHisto(HistoManager::Instance().fSpiceHistNumbers[0],(SpiceEnergy)*1000.);//+4keV for Bismuth
	  //fill standard energy spectra - no add-back

	  //Segment based energies
// 	  HistoManager::Instance().FillHisto(HistoManager::Instance().fSpiceHistNumbers[MAXNUMSEGSPICE*ring+seg+2],SpiceEnergy);
	  
	  //2D all seg energies
	  HistoManager::Instance().Fill2DHisto(HistoManager::Instance().fAngleDistro[1], (G4double) MAXNUMSEGSPICE*ring+seg,SpiceEnergy*1000., 1.0);
	  

	  fSpiceMultiplicity++;//iterates every deposition per fill to track multiplicity per event
	  energySumDet += SpiceEnergy;
	  
	  //Gated angular histos
	  if(abs(SpiceEnergyRaw - HistoManager::Instance().fBeamEnergy) < 10.*keV){//gating kinematics on (raw) detected energy
	    
		//Mapping phi to overlay mirrored segments
		G4int remainder = seg%3;
		G4double PhiMap = ((int)seg/3)*CLHEP::pi/2;
		PhiMap = HistoManager::Instance().fBeamPhi-PhiMap;
		PhiMap = asin(sin(PhiMap));//if we went over 360
	
	      //THETA vs. PHI
	      HistoManager::Instance().Fill2DHisto(HistoManager::Instance().fSpiceAngleHists[MAXNUMSEGSPICE*ring+remainder], HistoManager::Instance().fBeamTheta,PhiMap, 1.);//PHI M<AP!!!!!!!!!!!!!
	  }


	}//end of filling histos with SPICE detections per individual seg
      }//seg loop
    }//ring loop
    
    //wiritng add-backed energies
    if(energySumDet > MINENERGYTHRES) {//after exiting loops for all rings/segs, will input summed (add-back) energy if above threshold
      if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().fSpiceHistNumbers[1], energySumDet*1000.);
    }
}

void EventAction::FillPacesCryst() {
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDETPACES; j++) {
        if(fPacesCrystEnergyDet[j] > MINENERGYTHRES) {
	  
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().fPacesHistNumbers[0], fPacesCrystEnergyDet[j]);
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().fPacesHistNumbers[j+2], fPacesCrystEnergyDet[j]);
            energySumDet += fPacesCrystEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS && (energySumDet > MINENERGYTHRES))     HistoManager::Instance().FillHisto(HistoManager::Instance().fPacesHistNumbers[1], energySumDet);
	
    }
}

G4bool EventAction::SpiceTest(){//is SPICE inputted
  return HistoManager::Instance().Spice();
}

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


G4double EventAction::ApplySpiceRes(G4double EnergyIn){//pre-calculated functions that scale with energy to give accurate resolution split
  G4double GaussPick = G4UniformRand();//depending on result, will be part of one of two Gaussians, or a ERFC-bsaed decay
//   std::cout << "Rand: " << GaussPick <<std::endl; 
  G4double EnergyOut, ERFChold;
  if(GaussPick < 37./100.){
    
      G4double Sigma1 = (0.7E-6*EnergyIn/CLHEP::keV + 0.0009);//include ->keV conversion
  //     std::cout << "Sigma1: " << Sigma1 << " Energy: "<< EnergyIn/CLHEP::keV << std::endl;
      EnergyOut = G4RandGauss::shoot(EnergyIn, Sigma1);
  } else if(GaussPick > 77./100.){
      ERFChold = HistoManager::Instance().SPICE_ERFC();
      EnergyOut = ERFChold*(0.5E-6*EnergyIn + 0.0012)+EnergyIn/CLHEP::keV/1000.;//ERFC from (from http://radware.phy.ornl.gov/gf3/)
      
//        std::cout << "ERFC energy out: " << EnergyOut << std::endl;
//        std::cout << "ERFC changed out: " << ERFChold*(0.5E-6*EnergyIn + 0.0012)+EnergyIn/CLHEP::keV/1000. << std::endl;
  //     EnergyOut /= CLHEP::keV;
  } else {
      G4double Sigma2 = (1.8E-6*EnergyIn/CLHEP::keV + 0.00500);
  //    std::cout << "Sigma2: " << Sigma2 << std::endl;
      EnergyOut = G4RandGauss::shoot(EnergyIn, Sigma2);
//             std::cout << "non - ERFC energy out: " << EnergyOut << std::endl;
  }
  return EnergyOut;
}
