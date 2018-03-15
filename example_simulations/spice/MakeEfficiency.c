void MakeEfficiency(){
	
	std::vector<int> badchan;

	//Set the numbers of the SPICE channels that are disable or comment out
// 	badchan = {0,99,12};
	
	
	int SIMEVENTS=1000000;
	
        TFile output("SPICEefficiency.root", "RECREATE");
	output.mkdir("Ehist");
	gROOT->cd();
	TGraph EffCurve;
	
	// All input files should have the format "SPICE_750.root"
	for (int k = 1; k<=2000; k++) {
		
		// See if file for this energy exists
		stringstream ss;
		ss<<"SPICE_"<<k;
		TFile Einp((ss.str()+".root").c_str(),"READ");
		if(!Einp.IsOpen())continue;
		gROOT->cd();
		
		//Load the two histograms
		TH1 *full = (TH1*)Einp.Get("histo/FullEdep");
		if(!full)continue;
		TH2 *echan = (TH2*)Einp.Get("histo/AllSegEnergies");
		if(!echan)continue;
		
		// Gets the number of simm events from the file
		// This means you could decrease statistics for certain energy ranges 
		TH1 *be = (TH1*)Einp.Get("histo/BeamEnergy");
		if(be)SIMEVENTS=be->GetEntries();
		
		//Subtract out any disabled segments of the detector
		for(unsigned int i=0;i<badchan.size();i++){
			TH1* temph=echan->ProjectionY("tmp",i+1,i+1);
			full->Add(temph,-1);
			delete temph;
		}
		full->Sumw2(kFALSE);
		
		double rangeup=k+10;
		double rangedown=k-40-(4000./k);//Widen range for low energy where target effects large
		
		int binup=full->GetXaxis()->FindBin(rangeup);
		int bindown=full->GetXaxis()->FindBin(rangedown);
		
		double sum=full->Integral(bindown,binup);
		
		cout<<endl<<"COUNT "<<k<<" "<<sum<<endl;
		
		double efficiency=sum/SIMEVENTS;
		
		EffCurve.SetPoint(EffCurve.GetN(),k,efficiency);
		
		output.cd("Ehist");
			full->SetTitle(ss.str().c_str());
			full->Write(ss.str().c_str());
		gROOT->cd();
		Einp.Close();
	}
	
	output.cd();
		EffCurve.Write("EffCurve");
	gROOT->cd();
	output.Close();
}



