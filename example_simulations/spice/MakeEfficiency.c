void MakeEfficiency(){
    
	ifstream data("EffPoints.txt");
    if(!data.is_open()){
        cout<<endl<<"NO DATA"<<endl;
        return;
    }
    int N;
    data>>N;
    
    TGraphErrors Eff;
    double E,i;
    while(data>>E>>i){
        Eff.SetPoint(Eff.GetN(),E,i/N);
        if(i>=0)Eff.SetPointError(Eff.GetN()-1,0,sqrt(i)/N);
    }
    
    TFile output("SPICEefficiency.root", "RECREATE");
        Eff.Write("EffCurve");
	gROOT->cd();
	output.Close();
}



