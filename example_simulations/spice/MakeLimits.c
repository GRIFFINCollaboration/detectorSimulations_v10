void MakeLimits(){
	
	bool FitPeaks=true;
	int SIMEVENTS=1000000;
	
        TFile output("SPICELimits.root", "RECREATE");
	output.mkdir("Ehist");
	gROOT->cd();
	TGraph EffCurve;
	TGraph FitEffCurve;
	
	
	TGraph shad[2][10];
	
	// All input files should have the format "SPICE_750.root"
	for (int k = 1; k<=2000; k++) {
		
		// See if file for this energy exists
		stringstream ss;
		ss<<"SPICE_"<<k;
		TFile Einp((ss.str()+".root").c_str(),"READ");
		if(!Einp.IsOpen())continue;
		gROOT->cd();
		
		//Load the two histograms
		TH2 *echan = (TH2*)Einp.Get("histo/AllSegEnergies");
		if(!echan)continue;
		
		// Gets the number of simm events from the file
		// This means you could decrease statistics for certain energy ranges 
		TH1 *be = (TH1*)Einp.Get("histo/BeamEnergy");
		if(be)SIMEVENTS=be->GetEntries();
		
		TH1* shadhist[2][10];
		TH1* nadhist[10];
		
		for(unsigned int i=0;i<10;i++){
			for(unsigned int j=0;j<2;j++){
				stringstream nn;
				nn<<"tmp"<<j<<i;
				shadhist[j][i]=echan->ProjectionY(nn.str().c_str(),1,1);
				shadhist[j][i]->Reset();
			}
		}
		
		//Subtract out any disabled segments of the detector
		for(unsigned int i=0;i<120;i++){
			TH1* temph=echan->ProjectionY("tmp",i+1,i+1);
			int ring=i/12;
			if(i%3==1)shadhist[0][ring]->Add(temph);
			else shadhist[1][ring]->Add(temph);
			delete temph;
		}
		
		double rangeup=k+10;
		double rangedown=k-40-(4000./k);//Widen range for low energy where target effects large
		
		int binup=shadhist[0][0]->GetXaxis()->FindBin(rangeup);
		int bindown=shadhist[0][0]->GetXaxis()->FindBin(rangedown);
		
		for(unsigned int i=0;i<10;i++){
			for(unsigned int j=0;j<2;j++){
				double sum=shadhist[j][i]->Integral(bindown,binup);
				double efficiency=(sum/SIMEVENTS)/(4.*(j+1));
				shad[j][i].SetPoint(shad[j][i].GetN(),k,efficiency);
				
				delete shadhist[j][i];
			}
		}
		
		Einp.Close();
	}
	
	double limit=0.0002;
	
	TGraph shadfix[2][10];
	
	output.cd();
	for(unsigned int i=0;i<10;i++){
		for(unsigned int j=0;j<2;j++){
			for(unsigned int p=0;p<shad[j][i].GetN();p++){
				double x,y;
				shad[j][i].GetPoint(p,x,y);
				if(y<limit){
					shadfix[j][i].SetPoint(shadfix[j][i].GetN(),x,y/limit);
				}
				
				if(p+1<shad[j][i].GetN()){
					double xx,yy;
					shad[j][i].GetPoint(p+1,xx,yy);
					
					double grad=(yy-y)/(xx-x);
					double limpos=(((limit-y)/grad)+x);
					
					if((limpos>x)&&(limpos<xx)){
						shadfix[j][i].SetPoint(shadfix[j][i].GetN(),limpos,1);
					}
				}else{	
					if(y>limit){
						shadfix[j][i].SetPoint(shadfix[j][i].GetN(),x,1);
					}
				}
			}
			stringstream name;
			if(!j)name<<"shadow";
			else name<<"noshad";
			name<<i;
			shad[j][i].Write(name.str().c_str());
			name<<"n";
			shadfix[j][i].Write(name.str().c_str());
		}
	}
	
	output.Close();
}



