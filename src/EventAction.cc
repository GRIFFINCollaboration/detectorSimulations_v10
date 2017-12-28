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

#include "RunAction.hh"
#include "HistoManager.hh"

#include "Randomize.hh"

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
    BeamInputEnergy = HistoManager::Instance().BeamEnergy();
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
    
// Rishita Gudapati & Anita Mathews ------------------------------------------------------------------------------------------------------------------------------------------------------------------

//	incorrect angles?
//	G4double arr_gg[] = {33.6541, 53.8336, 66.4608, 76.3814, 88.4736, 101.331, 112.5438, 119.8489, 135.6357, 161.2132, 180.00};

	// used for gamma-gamma correlations
	G4double arr_gg[MAXNUMANG_GG] = {33.166, 53.690, 66.891, 76.694, 88.418, 101.57, 112.95, 119.84, 135.66, 160.87, 180.00};
	
	// used for all other correlations
	G4double arr[MAXNUMANG_GE] = {5.237, 16.203, 24.515, 35.016, 44.340, 55.517, 65.511, 74.412, 84.969, 94.932, 105.660, 114.784, 124.771, 135.490, 145.176, 155.214, 164.536, 174.763};

// All 52 angles for germanium-germanium hits, along with their weights, are below.
// to access the angle, use arr_gg[k][0]
// only 11 angles are used. To change this, reset MAXNUMANG_GG in the HistoManager header.

/*	G4double arr_gg[52][2] = {
        {0.0000, 64},
        {19.131, 128},
        {25.235, 64},
        {27.184, 64},
        {31.860, 64},
        {33.166, 48},
        {44.341, 128},
        {46.607, 96},
        {48.703, 128},
        {49.631, 96},
        {53.690, 48},
        {60.157, 96},
        {62.720, 48},
        {63.403, 64},
        {65.195, 96},
        {66.891, 64},
        {67.049, 64},
        {69.473, 64},
        {71.054, 96},
        {72.817, 64},
        {76.694, 96},
        {78.429, 64},
        {82.965, 64},
        {86.164, 64},
        {86.721, 48},
        {88.418, 128},
        {91.582, 128},
        {93.279, 48},
        {93.836, 64},
        {97.035, 64},
        {101.57, 64},
        {103.31, 96},
        {107.18, 64},
        {108.95, 96},
        {110.53, 64},
        {112.95, 64},
        {113.11, 64},
        {114.81, 96},
        {116.60, 64},
        {117.28, 48},
        {119.84, 96},
        {126.31, 48},
        {130.37, 96},
        {131.30, 128},
        {133.39, 96},
        {135.66, 128},
        {146.83, 48},
        {148.14, 64},
        {152.82, 64},
        {154.77, 64},
        {160.87, 128},
        {180.00, 64}
    };

// All 18 angles for paces-germanium hits, along with their weights, are below.
// angle == mean value, not bin centre.

	G4double arr[18][2] = 
	{{5.237, 2},
	{16.203, 3},
	{24.515, 15},
	{35.016, 13},
	{44.340, 23},
	{55.517, 18},
	{65.511, 21},
	{74.412, 21},
	{84.969, 29},
	{94.932, 29},
	{105.660, 24},
	{114.784, 20},
	{124.771, 19},
	{135.490, 24},
	{145.176, 16},
	{155.214, 17},
	{164.536, 4},
	{174.763, 2}};
*/

	// satisfied only when the hit was either in a germanium (ID 1000) or in paces (ID 50)
	if((fHitTrackerI[6][i] == 1000) || (fHitTrackerI[6][i] == 50)) { 
			
		edep1 = fHitTrackerD[0][i]/keV;
		
		G4int det1 = fHitTrackerI[8][i];
	 	G4int cry1 = fHitTrackerI[7][i];
		G4int cry1_index = -1;

		if (fHitTrackerI[6][i] == 1000) {
			cry1_index = det1*MAXNUMCRYGRIFFIN + cry1;
		} else {
			cry1_index = det1;
		}
		
		// start at i + 1 to prevent double counting		
		for(G4int j = i + 1; j < fNumberOfHits; j++) {
			if ((fHitTrackerI[6][j] == 1000) || (fHitTrackerI[6][j] == 50)) {

				edep2 = fHitTrackerD[0][j]/keV;

				G4int det2 = fHitTrackerI[8][j];
				G4int cry2 = fHitTrackerI[7][j];
				G4int cry2_index = -1;

				if (fHitTrackerI[6][j] == 1000) {
					cry2_index = det2*MAXNUMCRYGRIFFIN + cry2;
				} else {
					cry2_index = det2;
				}

				G4double value = -1;
			
				if((fHitTrackerI[6][i] == 1000) && (fHitTrackerI[6][j] == 1000) && gg) { // germanium-germanium
					value = thisGriffinCryMap[cry1_index][cry2_index]; 
	        			for (G4int k = 0; k < MAXNUMANG_GG; k++) {
						if (arr_gg[k] > (value - GGBIN) && arr_gg[k] < (value + GGBIN)) {
						   HistoManager::Instance().Fill2DHisto(HistoManager::Instance().AngCorrNumbers[k*MAXANGCORRHISTO], edep1, edep2, 1.0);					
							break; // once k is found, there's no need to keep looking through the array for it
						}
					}
				} else if ((fHitTrackerI[6][i] == 1000) && (fHitTrackerI[6][j] == 50) && ge) { // germanium-paces
					value = thisPacesGriffinCryMap[cry2_index][cry1_index]; // first index is for paces
	        			for (G4int k = 0; k < MAXNUMANG_GE; k++) {
						if (arr[k] > (value - GEBIN) && arr[k] < (value + GEBIN)) {
							HistoManager::Instance().Fill2DHisto(HistoManager::Instance().AngCorrNumbers[k*MAXANGCORRHISTO+1], edep1, edep2, 1.0);
							break;
						}
					}
				} else if ((fHitTrackerI[6][i] == 50) && (fHitTrackerI[6][j] == 1000) && ge) { // paces-germanium
					value = thisPacesGriffinCryMap[cry1_index][cry2_index]; 
					for (G4int k = 0; k < MAXNUMANG_GE; k++) {
						if (arr[k] > (value - GEBIN) && arr[k] < (value + GEBIN)) {
							// is filled into the same histogram as a germanium-paces
							// edep1 and edep2 switched to keep germanium on x-axis
							HistoManager::Instance().Fill2DHisto(HistoManager::Instance().AngCorrNumbers[k*MAXANGCORRHISTO+1], edep2, edep1, 1.0);
							break;
						}
					}

				// a Paces-Paces angle array is necessary.
			/*	} else if ((fHitTrackerI[6][i] == 50) && (fHitTrackerI[6][j] == 50) && ee) { //paces-paces
					G4double value = thisPacesPacesCryMap[cry1_index][cry2_index]; 
					for (G4int k = 0; k < MAXNUMANG_EE; k++) {
						if (arr[k][0] > (value - 0.005) && arr[k][0] < (value + 0.005)) {
							HistoManager::Instance().Fill2DHisto(HistoManager::Instance().AngCorrNumbers[k*MAXANGCORRHISTO+2], edep1, edep2, 1.0);
							break;
						}
					}
			*/
				}
			
	        		HistoManager::Instance().FillHisto(HistoManager::Instance().fAngCorrAngles[0], value); // keeps track of all values produced in simulation
			} // type of hit on second detector	
		}// for(j = i + 1 to number of hits)

	} // type of hit on first detector
    }// for(i == 0 to numberofhits)    
// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    	for (G4int i = 0 ; i < fNumberOfSteps; i++) {
		HistoManager::Instance().FillStepNtuple(fStepTrackerI[0][i], fStepTrackerI[1][i], fStepTrackerI[2][i], fStepTrackerI[3][i],  fStepTrackerI[4][i], 
			fStepTrackerI[5][i], fStepTrackerI[6][i], fStepTrackerI[7][i], fStepTrackerI[8][i], fStepTrackerD[0][i]/keV, fStepTrackerD[1][i]/mm, 
			fStepTrackerD[2][i]/mm, fStepTrackerD[3][i]/mm, fStepTrackerD[4][i]/second, fStepTrackerI[9][i]);
    	}
    
    ClearVariables();
}// EventAction::EndOfEventAction


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

                if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinUnSupCry[MAXNUMDETGRIFFIN*j+1+i], fGriffinCrystEnergyDet[i][j]/keV);
                if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinUnSupCry[0], fGriffinCrystEnergyDet[i][j]/keV);
                
		if(!suppressorBackFired[i] && !suppressorExtensionFired[i] && !suppressorSideFired[i]) { // Suppressor fired?
                  
		    if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinSupCry[MAXNUMDETGRIFFIN*j+1+i], fGriffinCrystEnergyDet[i][j]/keV);
                    if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinSupCry[0], fGriffinCrystEnergyDet[i][j]/keV);

                }
                energySumDet += fGriffinCrystEnergyDet[i][j];
            }
        }
        if(energySumDet > MINENERGYTHRES) {
            
	    // fill energies in each detector
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinUnSup[i+2], energySumDet/keV);
            
	    // fill standard energy and track spectra
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinUnSup[0], energySumDet/keV);
            
	    if(!suppressorBackFired[i] && !suppressorExtensionFired[i] && !suppressorSideFired[i]) {
                
		// fill energies in each detector
                if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinSup[i+2], energySumDet/keV);
                
		// fill standard energy and track spectra
                if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinSup[0], energySumDet/keV);

            }
        }
        energySum += energySumDet;
    }

    if(energySum > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinUnSup[1], energySum/keV);
        if(!suppressorFired) {
            if(WRITEEDEPHISTOS) HistoManager::Instance().FillHisto(HistoManager::Instance().GriffinSup[1], energySum/keV);
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
	      if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().fSpiceHistNumbers[0], G4RandGauss::shoot(fSpiceEnergyDet[ring][seg],0.002));
	      //fill standard energy spectra
	      if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().fSpiceHistNumbers[MAXNUMDETSPICE*ring+seg+2], G4RandGauss::shoot(fSpiceEnergyDet[ring][seg],0.002));
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
	  if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().fSpiceHistNumbers[1], energySumDet);
	  }
	  
	  if(energySumDet > (BeamInputEnergy-0.015)){//0.6 = 3 sigma 
	  if(fSpiceMultiplicity>0) HistoManager::Instance().FillHisto(HistoManager::Instance().fAngleDistro[0], fSpiceMultiplicity);
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
	  
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().fPacesHistNumbers[0], fPacesCrystEnergyDet[j]/keV);
            if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().fPacesHistNumbers[j+2], fPacesCrystEnergyDet[j]/keV);
            energySumDet += fPacesCrystEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     HistoManager::Instance().FillHisto(HistoManager::Instance().fPacesHistNumbers[1], energySumDet/keV);
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
