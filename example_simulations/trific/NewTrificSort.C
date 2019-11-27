// Grid data
int fYgrid=2;
int fXgrid=4;
double pitchmm=2;
double A=(60./180.)*TMath::Pi();

double targetgrid=53.4;
double windowgrid=2.8;

int XYN(int det){
    if(det==fYgrid||det==fYgrid-1)return 1;
    if(det==fXgrid||det==fXgrid-1)return 2;
    return 0;
}

double xp(int seg){
    int s=seg/2;
    double ret=3*pitchmm*(s+0.5);
    if(seg%2)ret=-ret;
    return ret;
}

double yp(int seg){
    int s=seg/2;
    double ret=-4*pitchmm*(s+0.5);
    if(seg%2)ret=-ret;
    return ret;
}   

void NewTrificSort(bool flatwindow=false, const char *rootout = "TRIFIC.root", const char *rootin = "g4out.root") {
    
    if(flatwindow)windowgrid=6.16;
    double targetwindow=targetgrid-windowgrid;
    
    
    TF1 tail("tail","pol1",0,20);
    tail.SetParameters(14000,1000);
    TF1 peak("peak","[1]-(x<[0])*abs([2])*pow(x-[0],2)-(x>=[0])*abs([3])*pow(x-[0],2)");
    peak.SetParameters(32,20000,1000,10000);
    

    // get the merged root file
    TFile *newfile = new TFile(rootin);
    if(!newfile->IsOpen())return;
    TTree *newtree = (TTree*)newfile->Get("ntuple");
    if(newtree==nullptr)return;
  
    
    double grid0mm=534;
    double dgmm=13;//In the Z axis, 10 mm along normal
    double zY=grid0mm+dgmm*fYgrid;
    double zX=grid0mm+dgmm*fXgrid;
    TVector3 gridnormal(0,-cos(A),sin(A));

    double EMAXGRID=40000;// Set based on ion
    
	
	int EventID;
	int Grid,Segment;
	double Energy;
    
	newtree->SetBranchAddress("eventNumber",&EventID);
	newtree->SetBranchAddress("detNumber",&Grid);
	newtree->SetBranchAddress("cryNumber",&Segment);
	newtree->SetBranchAddress("depEnergy",&Energy);

	TFile out(rootout,"RECREATE");
	out.cd();	
		TH2F* XYPos=new TH2F("XYPos","XY Position;X [mm];Y [mm]",200,-80,80,200,-80,80);

		TH1F* SumEnergy=new TH1F("SumEnergy","Sum Energy;Energy [arb]",10000,0,EMAXGRID*20);
		TH2F* SegmentEnergy=new TH2F("SegmentEnergy","Segment Energy;Segment;Energy [keV]",24,0,24,500,0,EMAXGRID);
		TH2F* ZdedxA=new TH2F("ZdedxA","Zdedx Y Corrected ;Z [cm];de/dx [keV/cm]",2000,0,45,500,0,EMAXGRID);
        
        TH3F* IDnew=new TH3F("IDnew","IDnew;Range [cm];de/dx peak;de/dx zero",100,20,40,200,0,EMAXGRID,200,0,EMAXGRID);
        
        TH2F* rangeZ=new TH2F("rangeZ","Range v max;Range [cm];de/dx peak",200,20,40,500,0,EMAXGRID);
        TH2F* range0=new TH2F("range0","Range v zero;Range [cm];de/dx zero",200,20,40,500,0,EMAXGRID);
        TH2F* Z0=new TH2F("Z0","max v zero;de/dx peak;de/dx zero",500,0,EMAXGRID,500,0,EMAXGRID);
        TH2F* red0Z=new TH2F("red0Z","max v reduced;de/dx peak;de/dx zero / E sum",500,0,EMAXGRID,500,0,0.4);
        TH2F* rangered0=new TH2F("rangered0","Range v reduced;Range [cm];de/dx zero / E sum",200,20,40,500,0,0.4);
        
	gROOT->cd();

	long nentries = newtree->GetEntries();
	int count=0;
	long CurrentEvent=0;
	
	vector<double> TrificE(24,0.0);//Holder for the event
	vector<double> TrificY(12,0.0);//Holder for the event
	vector<double> TrificX(12,0.0);//Holder for the event
    
// 	vector<double> fYmm;
// 	vector<double> fXmm;
// 	vector<double> fYOmm;
// 	vector<double> fXOmm;
    
//     double pos=0;
//     for(unsigned int i=0;i<fYwires.size();i++){
//         double step=pitchmm*fYwires[i];
//         fYmm.push_back(pos+step*0.5);
//         fYmm.push_back(-(pos+step*0.5));
//         pos+=step;
//         fYOmm.push_back(pos);
//         fYOmm.push_back(-pos);
//     }
    
//     pos=0;
//     for(unsigned int i=0;i<fXwires.size();i++){
//         double step=pitchmm*fXwires[i];
//         fXmm.push_back(pos+step*0.5);
//         fXmm.push_back(-(pos+step*0.5));
//         pos+=step;
//         fXOmm.push_back(pos);
//         fXOmm.push_back(-pos);
//     }  
//     double tPhi=0,tTheta=0;
    
    
    TRandom r;
	bool EndOfEvent=true;
    int Mult=0;
    double sum=0;
	for(long jentry=0;jentry<nentries;jentry++){
		newtree->GetEntry(jentry); 

		// Accumulate one event, detector data spread over several tree entries
        
		if(EventID==CurrentEvent){
			EndOfEvent=false;
			if(Grid>=0&&Grid<24){
				if(TrificE[Grid]>0&&!XYN(Grid)){
					cout<<endl<<"Error double fill of segment "<<Grid<<flush;
				}else{
					count++;
					TrificE[Grid]+=Energy;
                    sum+=Energy;
                    if(XYN(Grid)==1)TrificY[Segment]+=Energy;
                    if(XYN(Grid)==2)TrificX[Segment]+=Energy;
				}
			}
		}else{
			CurrentEvent=EventID;
			EndOfEvent=true;
			jentry--;
		}
		
		if(jentry==(nentries-1))EndOfEvent=true;
		if(!EndOfEvent) continue;
		
        SumEnergy->Fill(sum);
        
        /////////////////////////////////////////////////////////////////////////////////
        // In the G4 we have divided into grid cells, but real grids sense both neighbour cells
        // This section reconstructs all that.
        /////////////////////////////////////////////////////////////////////////////////
        
        int last=0;
        for(int i=23;i>0;i--){
            if(last==0&&TrificE[i]>0)last=i;
            TrificE[i]+=TrificE[i-1];
        }

        ///////////////////////////
        // We could add some additional random width and energy thresholds to the signals here to futher replicate real data.
        // However as we are currently conserned with the *relative* improvement comparing different setups this isnt needed   
        //////////////////////////
        
        /////////////////////////////////////
        // Angular correction calculations //
        /////////////////////////////////////
        
        if(last>fXgrid+2){
            double Ymean=0;
            double Xmean=0;
            
            for(int i=0;i<12;i++){
                Ymean+=TrificY[i]*yp(i);
                Xmean+=TrificX[i]*xp(i);
            }
            
            Xmean/=TrificE[fXgrid];
            Ymean/=TrificE[fYgrid];
            
            XYPos->Fill(Xmean,Ymean);
                
            double Yreal=sin(A)*Ymean;
            double Yzreal=zY+(cos(A)*Ymean);
            double theY=Yreal/Yzreal;
            double theX=Xmean/zX;
            TVector3 particle(theX,theY,1);
            double ZratioX=1./gridnormal.Dot(particle.Unit());
            
            
            TGraph tmpbragg;
            double maxE=0;
            double maxX=0;
            for(int i=1;i<24;i++){
                if(TrificE[i]>0){
                    SegmentEnergy->Fill(i,TrificE[i]);
                    
                    double basicZ=(dgmm*i*0.1)+windowgrid;
                    
                    if(flatwindow){
                        double dy=(particle.Unit().Y()/particle.Unit().Z())*targetwindow;
                        double dz=dy/tan(A);
                        basicZ+=dz;
                    }
                    
                    double dE=TrificE[i]/ZratioX;
                    double x=basicZ*ZratioX;
                    ZdedxA->Fill(x,dE);
                    tmpbragg.SetPoint(tmpbragg.GetN(),x,dE);
                
                    if(dE>maxE){
                        maxE=dE;
                        maxX=x;
                    }
                }
            }
            
            tail.SetParameters(14000,1000);
            tmpbragg.Fit(&tail,"QNR");
            
            peak.SetParameters(maxX,maxE,1000,10000);
            peak.SetRange(maxX-5,40);
            tmpbragg.Fit(&peak,"QNR");          
            double dE0=tail.Eval(5);
            double dedxpeak=peak.GetParameter(1);
            double range=peak.GetParameter(0)+sqrt(abs(dedxpeak/peak.GetParameter(3))); 
            IDnew->Fill(range,dedxpeak,dE0);
            
            rangeZ->Fill(range,dedxpeak);
            range0->Fill(range,dE0);
            Z0->Fill(dedxpeak,dE0);
            red0Z->Fill(dedxpeak,dE0/sum);
            rangered0->Fill(range,dE0/sum);
            
//             TCanvas C1;
//             gPad->Update();
//             tmpbragg.Fit(&tail,"+R");
//             tmpbragg.Fit(&peak,"+R");
//             cout<<endl<<"dE0 = "<<dE0;
//             cout<<endl<<"Z = "<<dedxpeak;
//             cout<<endl<<"range = "<<range<<endl;
//             tmpbragg.Draw("al");
//             C1.Modified();
//             C1.Update();
//             C1.WaitPrimitive();
                
                
        }
            
        //Reset things ready for next event
        count=0;
        std::fill(TrificE.begin(), TrificE.end(), 0);
        std::fill(TrificX.begin(), TrificX.end(), 0);
        std::fill(TrificY.begin(), TrificY.end(), 0);
        sum=0;
	}
    

    out.Write();
    delete newfile;
}
