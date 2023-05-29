//g++ SpiceMiniSort.C `root-config --cflags --libs` -o sMini
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TH1.h"
#include "TROOT.h"
#include <iostream>

using namespace std;

TH1* DrawSingles(TFile* file){
    gROOT->cd();
    TTree *newtree = (TTree*)file->Get("ntuple");
    if(!newtree)return 0;
    TH1* h1=new TH1D("h1","Esingles;Energy [keV];Counts/keV",2000,0,2000);
    newtree->Draw("depEnergy>>h1","","goff");
    return h1;
}

int SpiceMiniSort(double e=1000,double f=0.05) {
    // get the merged root file
    TFile *newfile = new TFile("g4out.root");
    gROOT->cd();
    TH1* E=DrawSingles(newfile);
    newfile->Close();
    if(E){
        TAxis *a=E->GetXaxis();
        int ret=E->Integral(a->FindBin(e*(1-f)),a->FindBin(e*(1+f)));       
        TFile* out= new TFile("SPICEm.root","RECREATE");
        E->Write();
        out->Close();
        delete E;
        return ret;
    }
    return 0;
}


int main(int argc, char *argv[]){
    double e=1000;
    double f=0.05;
    
    if(argc>1){
        stringstream ss;
        ss<<argv[1];
        ss>>e;
    }
    if(argc>2){
        stringstream ss;
        ss<<argv[2];
        ss>>f;
    }
    
    std::cout<<SpiceMiniSort(e,f);

    return 0;
}
