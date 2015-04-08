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

#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* run, HistoManager* histo)
    :G4UserEventAction(),
      fRunAct(run),fHistoManager(histo),
      fPrintModulo(0)
{
    numberOfHits = 0;
    numberOfSteps = 0;
    fPrintModulo = 100;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
    evtNb = evt->GetEventID();
    if (evtNb%fPrintModulo == 0)
        //    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
        printf( " ---> Ev.# %5d\r", evtNb);
    G4cout.flush();

    ClearVariables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
    //    FillParticleType() ;
    FillGriffinCryst() ;
    Fill8piCryst() ;
    FillLaBrCryst() ;
    FillAncillaryBgo() ;
    FillSceptar() ;
    FillGridCell() ;

    //G4cout << "numberOfHits = " << numberOfHits << G4endl;
    for (G4int i = 0 ; i < numberOfHits; i++) {
        fHistoManager->FillHitNtuple(hitTrackerI[0][i], hitTrackerI[1][i], hitTrackerI[2][i], hitTrackerI[3][i],  hitTrackerI[4][i], hitTrackerI[5][i], hitTrackerI[6][i], hitTrackerI[7][i], hitTrackerI[8][i], hitTrackerD[0][i]/keV, hitTrackerD[1][i]/mm, hitTrackerD[2][i]/mm, hitTrackerD[3][i]/mm, hitTrackerD[4][i]/second);
    }
    for (G4int i = 0 ; i < numberOfSteps; i++) {
        fHistoManager->FillStepNtuple(stepTrackerI[0][i], stepTrackerI[1][i], stepTrackerI[2][i], stepTrackerI[3][i],  stepTrackerI[4][i], stepTrackerI[5][i], stepTrackerI[6][i], stepTrackerI[7][i], stepTrackerI[8][i], stepTrackerD[0][i]/keV, stepTrackerD[1][i]/mm, stepTrackerD[2][i]/mm, stepTrackerD[3][i]/mm, stepTrackerD[4][i]/second);
    }

    ClearVariables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




void EventAction::ClearVariables()
{
    if(fHistoManager->GetStepTrackerBool()) {
        stepIndex = 0;
        numberOfSteps = 0;
        for (G4int i = 0 ; i < MAXSTEPS; i++) {
            for (G4int j = 0 ; j < NUMSTEPVARS; j++) {
                stepTrackerI[j][i] = 0;
                stepTrackerD[j][i] = 0.0;
            }
        }
    }

    if(fHistoManager->GetHitTrackerBool()) {
        hitIndex = 0;
        numberOfHits = 0;
        pTrackID = -1;
        pParentID = -1;

        for (G4int i = 0 ; i < MAXHITS; i++) {
            pHitMnemonic[i] = "XXX00XX00X";
            for (G4int j = 0 ; j < NUMSTEPVARS; j++) {
                hitTrackerI[j][i] = 0;
                hitTrackerD[j][i] = 0.0;
            }
        }
    }

    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++) {
        particleTypes[i]                  = 0;
    }
    for (G4int i = 0 ; i < MAXNUMDET; i++) {
        EightPiCrystEnergyDet[i]  = 0 ;
        EightPiCrystTrackDet[i]   = 0 ;

        LaBrCrystEnergyDet[i]  = 0 ;
        LaBrCrystTrackDet[i]   = 0 ;

        AncillaryBgoEnergyDet[i]  = 0 ;
        AncillaryBgoTrackDet[i]   = 0 ;

        SceptarEnergyDet[i]  = 0 ;
        SceptarTrackDet[i]   = 0 ;

        GridCellElectronEKinDet[i]  = 0 ;
        GridCellElectronTrackDet[i]   = 0 ;
        GridCellGammaEKinDet[i]  = 0 ;
        GridCellGammaTrackDet[i]   = 0 ;
        GridCellNeutronEKinDet[i]  = 0 ;
        GridCellNeutronTrackDet[i]   = 0 ;

    }
    for (G4int i = 0 ; i < MAXNUMDETGRIFFIN; i++) {
        for (G4int j = 0 ; j < MAXNUMCRYGRIFFIN; j++) {

            GriffinCrystEnergyDet[i][j]                 = 0;
            GriffinCrystTrackDet[i][j]                  = 0;
            GriffinSuppressorBackEnergyDet[i][j]        = 0;
            GriffinSuppressorBackTrackDet[i][j]         = 0;
            GriffinSuppressorLeftExtensionEnergyDet[i][j]   = 0;
            GriffinSuppressorLeftExtensionTrackDet[i][j]    = 0;
            GriffinSuppressorLeftSideEnergyDet[i][j]        = 0;
            GriffinSuppressorLeftSideTrackDet[i][j]         = 0;
            GriffinSuppressorRightExtensionEnergyDet[i][j]   = 0;
            GriffinSuppressorRightExtensionTrackDet[i][j]    = 0;
            GriffinSuppressorRightSideEnergyDet[i][j]        = 0;
            GriffinSuppressorRightSideTrackDet[i][j]         = 0;
        }
    }
    // NOTE: Clear the variables from the new Fill___Cryst functions.
}


void EventAction::FillParticleType()
{
    G4int numParticleTypes = 0;
    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++)
    {
        if (particleTypes[i] != 0)
        { // if particle type 'i' has non-zero counts
            for (G4int j = 0 ; j< particleTypes[i]; j++)
            { // loop over the number of time we saw it
                G4cout << "particleTypes[" << i << "] = " << particleTypes[i] << G4endl;
                fHistoManager->FillHisto(astats_particle_type_in_each_step, i);
            }
        }
    }

    // Fill the number of particle types in the event
    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++)
    {
        if (particleTypes[i] != 0)
            numParticleTypes++;
    }
    fHistoManager->FillHisto(astats_particle_type_in_each_event, numParticleTypes);
}

void EventAction::FillGriffinCryst()
{
    G4double  energySum = 0;
    G4double  energySumDet = 0;
    G4bool SuppressorBackFired[MAXNUMDETGRIFFIN] = {0};
    G4bool SuppressorExtensionFired[MAXNUMDETGRIFFIN] = {0};
    G4bool SuppressorSideFired[MAXNUMDETGRIFFIN] = {0};
    G4bool SuppressorFired = false;

    // Fill Griffin Histos
    for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
        energySumDet = 0;
        // Find if any suppressors were fired
        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            if (GriffinSuppressorBackEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorBackFired[i] = true;
            }
            if (GriffinSuppressorLeftExtensionEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorExtensionFired[i] = true;
            }
            if (GriffinSuppressorRightExtensionEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorExtensionFired[i] = true;
            }
            if (GriffinSuppressorLeftSideEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorSideFired[i] = true;
            }
            if (GriffinSuppressorRightSideEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorSideFired[i] = true;
            }
            if ( !SuppressorFired && ( SuppressorBackFired[i] || SuppressorExtensionFired[i] || SuppressorSideFired[i] ) ) {
                SuppressorFired = true;
            }
        }

        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            if(GriffinCrystEnergyDet[i][j] > MINENERGYTHRES) {
                // fill energies in each crystal
                if(WRITEEDEPHISTOS)     fHistoManager->FillHisto((griffin_crystal_unsup_edep_det0_cry0+(MAXNUMDETGRIFFIN*j))+i, GriffinCrystEnergyDet[i][j]);
                if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(griffin_crystal_unsup_edep_cry, GriffinCrystEnergyDet[i][j]);
                if(!SuppressorBackFired[i] && !SuppressorExtensionFired[i] && !SuppressorSideFired[i]) { // Suppressor fired?
                    if(WRITEEDEPHISTOS) fHistoManager->FillHisto((griffin_crystal_sup_edep_det0_cry0+(MAXNUMDETGRIFFIN*j))+i, GriffinCrystEnergyDet[i][j]);
                    if(WRITEEDEPHISTOS) fHistoManager->FillHisto(griffin_crystal_sup_edep_cry, GriffinCrystEnergyDet[i][j]);
                }
                energySumDet += GriffinCrystEnergyDet[i][j];
            }
        }
        if(energySumDet > MINENERGYTHRES) {
            // fill energies in each detector
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(griffin_crystal_unsup_edep_det0+i, energySumDet);
            // fill standard energy and track spectra
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(griffin_crystal_unsup_edep, energySumDet);
            if(!SuppressorBackFired[i] && !SuppressorExtensionFired[i] && !SuppressorSideFired[i]) {
                // fill energies in each detector
                if(WRITEEDEPHISTOS) fHistoManager->FillHisto(griffin_crystal_sup_edep_det0+i, energySumDet);
                // fill standard energy and track spectra
                if(WRITEEDEPHISTOS) fHistoManager->FillHisto(griffin_crystal_sup_edep, energySumDet);
            }
        }
        energySum += energySumDet;
    }

    if(energySum > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(griffin_crystal_unsup_edep_sum, energySum);
        if(!SuppressorFired) {
            if(WRITEEDEPHISTOS) fHistoManager->FillHisto(griffin_crystal_sup_edep_sum, energySum);
        }
    }
}

void EventAction::Fill8piCryst()
{
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(EightPiCrystEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(Eightpi_crystal_edep, EightPiCrystEnergyDet[j]);
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(Eightpi_crystal_edep_det0+j, EightPiCrystEnergyDet[j]);
            energySumDet += EightPiCrystEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(Eightpi_crystal_edep_sum, energySumDet);
    }

}


void EventAction::FillLaBrCryst()
{
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(LaBrCrystEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(labr_crystal_edep, LaBrCrystEnergyDet[j]);
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(labr_crystal_edep_det0+j, LaBrCrystEnergyDet[j]);
            energySumDet += LaBrCrystEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(labr_crystal_edep_sum, energySumDet);
    }
}

void EventAction::FillAncillaryBgo()
{
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(AncillaryBgoEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(ancillary_bgo_crystal_edep, AncillaryBgoEnergyDet[j]);
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(ancillary_bgo_crystal_edep_det0+j, AncillaryBgoEnergyDet[j]);
            energySumDet += AncillaryBgoEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(ancillary_bgo_crystal_edep_sum, energySumDet);
    }
}

void EventAction::FillSceptar()
{
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(SceptarEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(sceptar_edep, SceptarEnergyDet[j]);
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(sceptar_edep_det0+j, SceptarEnergyDet[j]);
            energySumDet += SceptarEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(sceptar_edep_sum, energySumDet);
    }

}

void EventAction::FillGridCell()
{
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(WRITEEKINHISTOS && GridCellElectronEKinDet[j] > MINENERGYTHRES)     fHistoManager->FillHisto(gridcell_electron_ekin_det0+j, GridCellElectronEKinDet[j]);
        if(WRITEEKINHISTOS && GridCellGammaEKinDet[j] > MINENERGYTHRES)     fHistoManager->FillHisto(gridcell_gamma_ekin_det0+j, GridCellGammaEKinDet[j]);
        if(WRITEEKINHISTOS && GridCellNeutronEKinDet[j] > MINENERGYTHRES)     fHistoManager->FillHisto(gridcell_neutron_ekin_det0+j, GridCellNeutronEKinDet[j]);
    }
}

void EventAction::AddHitTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time)
{
    G4bool newhit = true;
    for (G4int i = 0 ; i < numberOfHits; i++) {
        if(pHitMnemonic[i] == mnemonic) {
            // sum the new enery
            hitTrackerD[0][i] = hitTrackerD[0][i] + depEnergy;
            newhit = false;
            break;
        }
    }
    if (newhit) { // new hit
        pHitMnemonic[hitIndex] = mnemonic;
        pTrackID = trackID;
        pParentID = parentID;
        hitTrackerI[0][hitIndex] = eventNumber;
        hitTrackerI[1][hitIndex] = trackID;
        hitTrackerI[2][hitIndex] = parentID;
        hitTrackerI[3][hitIndex] = stepNumber;
        hitTrackerI[4][hitIndex] = particleType;
        hitTrackerI[5][hitIndex] = processType;
        hitTrackerI[6][hitIndex] = systemID;
        hitTrackerI[7][hitIndex] = cryNumber;
        hitTrackerI[8][hitIndex] = detNumber;
        hitTrackerD[0][hitIndex] = depEnergy;
        hitTrackerD[1][hitIndex] = posx;
        hitTrackerD[2][hitIndex] = posy;
        hitTrackerD[3][hitIndex] = posz;
        hitTrackerD[4][hitIndex] = time;

        hitIndex++;
        numberOfHits = hitIndex;

        if(numberOfHits >= MAXHITS) {
            G4cout << "ERROR! Too many hits!" << G4endl;
        }
    }
}


void EventAction::AddStepTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time)
{
    G4bool newstep = true;
    if (newstep) { // new step
        stepTrackerI[0][stepIndex] = eventNumber;
        stepTrackerI[1][stepIndex] = trackID;
        stepTrackerI[2][stepIndex] = parentID;
        stepTrackerI[3][stepIndex] = stepNumber;
        stepTrackerI[4][stepIndex] = particleType;
        stepTrackerI[5][stepIndex] = processType;
        stepTrackerI[6][stepIndex] = systemID;
        stepTrackerI[7][stepIndex] = cryNumber;
        stepTrackerI[8][stepIndex] = detNumber;
        stepTrackerD[0][stepIndex] = depEnergy;
        stepTrackerD[1][stepIndex] = posx;
        stepTrackerD[2][stepIndex] = posy;
        stepTrackerD[3][stepIndex] = posz;
        stepTrackerD[4][stepIndex] = time;

        stepIndex++;

        numberOfSteps = stepIndex;
        if(numberOfSteps >= MAXSTEPS) {
            G4cout << "ERROR! Too many steps!" << G4endl;
        }
    }
}



