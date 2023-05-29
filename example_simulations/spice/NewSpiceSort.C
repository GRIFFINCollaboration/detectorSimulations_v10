TRandom r1;

double SpiceResolution(double energy) {
    cout<<endl<<energy;
    double sigma=0.7E-3*energy + 0.9;
    if(r1.Uniform()<0.4)sigma=1.8E-3*energy + 2.0;
    //Generally observe 2ish part resolution, possibly channel based rather than random
    
    if(r1.Uniform()<0.33){
		energy -= r1.Exp(1.5E-3*energy + 1.2);
        //Charge trapping tail
    }
    energy=r1.Gaus(energy,sigma);
    cout<<" "<<energy;
    return energy;
}

void NewSpiceSort(const char *root1 = 0, const char *rootout = 0) {

    bool USERESOLUTION=true;
    r1.SetSeed();
    
    // get the merged root file
    TFile *newfile = new TFile("g4out.root");
    TTree *newtree = (TTree*)newfile->Get("ntuple");
    
	int EventID,Ring,Sector;
	double Energy;
	double posx,posy,posz;
	double beamE,beamTheta,beamPhi;
    
	newtree->SetBranchAddress("eventNumber",&EventID);
	newtree->SetBranchAddress("detNumber",&Ring);
	newtree->SetBranchAddress("cryNumber",&Sector);
	newtree->SetBranchAddress("depEnergy",&Energy);
	newtree->SetBranchAddress("originX",&posx);
	newtree->SetBranchAddress("originY",&posy);
	newtree->SetBranchAddress("originZ",&posz);
	newtree->SetBranchAddress("primaryE",&beamE);
	newtree->SetBranchAddress("primaryTheta",&beamTheta);
	newtree->SetBranchAddress("primaryPhi",&beamPhi);
    
	TFile* out= new TFile("SPICE.root","RECREATE");
	out->cd();	
        TH1D* Edep=new TH1D("Edep","Energy Dep Singles;Energy [keV];Counts/keV",2000,0,2000);
        TH1D* Edepadd=new TH1D("Edepadd","Energy Dep Addback;Energy [keV];Counts/keV",2000,0,2000);

		TH2F* SegEnergy=new TH2F("SegEnergy","Segment vs Energy;Segment;Energy [keV];Counts/keV",120,0,120,2000,0,2000);
    
        out->mkdir("BeamDist");
        out->cd("BeamDist");
        		TH2F* SourceXY=new TH2F("SourceXY","Soource XY Position;X [mm];Y [mm]",200,-5,5,200,-5,5);
        		TH1D* SourceZ=new TH1D("SourceZ","Z Position;Z [mm];",100000,-1.0,10.0);
        		TH1D* SourceE=new TH1D("SourceE","Source Energy;Source Energy [keV];Counts/keV",2000,0,2000);
        out->cd();	
        
        out->mkdir("AngDist");
        out->cd("AngDist");
            TH1* AngDist[20];
            for(int r=0;r<10;r++){
            for(int s=0;s<2;s++){
                stringstream ss;
                ss<<"AngR"<<r<<"S"<<s;
        		AngDist[r+s*10]=new TH2F(ss.str().c_str(),(ss.str()+";Theta [rad];Phi [rad]").c_str(),100,0,TMath::Pi(),200,-TMath::Pi(),TMath::Pi());
            }}
        out->cd();	
        
	gROOT->cd();

	long nentries = newtree->GetEntries();
	int count=0;
	long CurrentEvent=0;
	
	vector<double> SPICE_E(120,0.0);//Holder for the event
	double Esum=0;
    
	bool EndOfEvent=true;
    
    double bE,bX,bY,bZ;
    
    TVector3 beam;

	for(long jentry=0;jentry<nentries;jentry++){
		newtree->GetEntry(jentry); 

		// Accumulate one event, detector data spread over several tree entries
        
		if(EventID==CurrentEvent){
			EndOfEvent=false;
			if(Ring>=0&&Ring<10){
                if(Sector>=0&&Sector<12){
                    int seg=Ring*12+Sector;
                    if(SPICE_E[seg]>0){
                        cout<<endl<<"Error double fill of segment "<<seg<<flush;
                    }else{
                        count++;
                        if(USERESOLUTION){Energy=SpiceResolution(Energy);}
                        SPICE_E[seg]=Energy;
                        Esum+=Energy;
                    }
                }
            }
			
            bX=posx;
            bY=posy;
            bZ=posz;
            bE=beamE;
            beam.SetMagThetaPhi(1,beamTheta,beamPhi);
		}else{
			CurrentEvent=EventID;
			EndOfEvent=true;
			jentry--;
		}
		
		if(jentry==(nentries-1))EndOfEvent=true;
		
		if(!EndOfEvent) continue;
		
        // We have all the detector data from a single event so process it bellow.        
        // First build the actual observable signals from the G4 simplifications

        SourceXY->Fill(bX,bY);
        SourceZ->Fill(bZ);
        SourceE->Fill(bE);
        
        
        for(int i=0;i<120;i++){
            double E=SPICE_E[i];
            if(E>0){
                Edep->Fill(E);
                SegEnergy->Fill(i,E);
                
                // Angular stuff
                
                // Only fill these for non-scattered events
                if(count>1)continue;
                if(abs(E-bE)/bE>0.05)continue;
                
                // Various mapings rotations and reflections
                // There are only really 20 unique segments by rotational/reflection symetry
                // So this effectively increases the stats of the distributions by x5

                int r=i/12;
                int s=i%12;
                int q=s/3;
                int ss=s%3;
                
                beam.RotateZ(-(q+2)*TMath::Pi()/2.0);
        
                if(ss==2){
                    ss=0;
                    beam.SetPhi(-beam.Phi());
                    beam.RotateZ(TMath::Pi()/2.0);
                }
                
                AngDist[r+ss*10]->Fill(beam.Theta(),beam.Phi());
            }
        }
        Edepadd->Fill(Esum);
        
        // Reset for next event
        count=0;
        std::fill(SPICE_E.begin(), SPICE_E.end(), 0);
        Esum=0;
        
	}
    
    out->Write();
    out->Close();
    newfile->Close();
    new TBrowser;
}
