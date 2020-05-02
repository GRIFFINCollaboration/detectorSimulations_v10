
// Sort control bool
bool UseFullBraggFit=true;
bool DoTestDraw=false;
bool SimulateElectronic=true;
	
// These factors should be tuned to data when possible
double noisefactor=300;
double threshold=0+noisefactor;

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

void LatestTrificSort(bool flatwindow=false, const char *rootout = "TRIFIC.root", const char *rootin = "g4out.root") {
    
    if(flatwindow)windowgrid=6.16;
    double targetwindow=targetgrid-windowgrid;
    
    
    TF1 tail("tail","pol1",0,20);
    tail.SetParameters(14000,1000);
    TF1 peak("peak","[1]-(x<[0])*abs([2])*pow(x-[0],2)-(x>=[0])*abs([3])*pow(x-[0],2)");
    peak.SetParameters(32,20000,1000,10000);
    
    TF1 newbragg("newbragg","[0]*exp(([1]*[1]*[2]*[2]*0.5)+([1]*(x-[3])))*TMath::Erfc(([1]*[2]/sqrt(2))-(([3]-x)/([2]*sqrt(2))))",0,50);
	newbragg.SetLineColor(3);

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
	r.SetSeed();
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
        // Now that we have converted from the geant4 events to the signals we would expect out of the detector
        // We add electronic effects of noise and thresholds
        //////////////////////////
   
if(SimulateElectronic){
		for(int i=23;i>0;i--){
			if(TrificE[i]<threshold){
				TrificE[i]=0;
				if(last==i)last--;
			}else{
				TrificE[i]+=r.Gaus(0,noisefactor);
			}
        }
        
		TrificE[fXgrid]=0;
		TrificE[fYgrid]=0;
		for(int i=0;i<12;i++){		
			if(TrificY[i]<threshold){
				TrificY[i]=0;
			}else{
				TrificY[i]+=r.Gaus(0,noisefactor);
			}
			TrificE[fYgrid]+=TrificY[i];
        
			if(TrificX[i]<threshold){
				TrificX[i]=0;
			}else{
				TrificX[i]+=r.Gaus(0,noisefactor);
			}
			TrificE[fXgrid]+=TrificX[i];
		}
}
        
        /////////////////////////////////////
        // Angular correction calculations //
        /////////////////////////////////////
        
        if(last>fXgrid+2&&TrificE[fXgrid]>0&&TrificE[fYgrid]>0){
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
            

            
            peak.SetParameters(maxX,maxE,1000,10000);
            peak.SetRange(maxX-5,40);
            tmpbragg.Fit(&peak,"QNR");          
            double dE0=TrificE[0];
            double dedxpeak=peak.GetParameter(1);
            double range=peak.GetParameter(0)+sqrt(abs(dedxpeak/peak.GetParameter(3))); 
		
if(UseFullBraggFit){
	
			double xmax=peak.GetParameter(0);
			double sig=range-xmax;
			if(sig<1)sig=1;
			newbragg.SetParameter(0,dedxpeak*0.5);
			newbragg.SetParLimits(0,dedxpeak*0.4,dedxpeak*0.6);
			newbragg.SetParameter(1,0.02);
			newbragg.SetParameter(2,sig);
			newbragg.SetParLimits(2,sig*0.1,sig*2);
			newbragg.SetParameter(3,xmax+sig*0.5);
			newbragg.SetParLimits(3,xmax-sig,xmax+sig*2);
            tmpbragg.Fit(&newbragg,"QNR");
         
            dE0=newbragg.Eval(5);
			range=newbragg.GetParameter(3)+newbragg.GetParameter(2)*2;// Approx
			dedxpeak=newbragg.Eval(newbragg.GetParameter(3)-newbragg.GetParameter(2)*2.3);// Approx
}else{
            tail.SetParameters(14000,1000);
            tmpbragg.Fit(&tail,"QNR");      
            dE0=tail.Eval(5);
}
			
            IDnew->Fill(range,dedxpeak,dE0);
            rangeZ->Fill(range,dedxpeak);
            range0->Fill(range,dE0);
            Z0->Fill(dedxpeak,dE0);
            red0Z->Fill(dedxpeak,dE0/sum);
            rangered0->Fill(range,dE0/sum);
            
			
///////////// The fit may need tuning for different regions ///////////
///////////// Use these lines to look at it while testing   ///////////
if(DoTestDraw){
            TCanvas C1;
            gPad->Update();
            tmpbragg.Fit(&tail,"+R");
            tmpbragg.Fit(&peak,"+R");
			
if(UseFullBraggFit){tmpbragg.Fit(&newbragg,"+R");}

            cout<<endl<<"dE0 = "<<dE0;
            cout<<endl<<"Z = "<<dedxpeak;
            cout<<endl<<"range = "<<range<<endl;
			tmpbragg.SetMarkerStyle(20);
            tmpbragg.Draw("ap");
			
            C1.Modified();
            C1.Update();
            C1.WaitPrimitive();
}

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
