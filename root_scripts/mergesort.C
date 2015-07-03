// This ROOT macro appends 1 (up to 10) ntuple files together and then sorts them according to the event number

void mergesort(const char *root1 = 0, const char *rootout = 0) {

    // define strings
    TString string1 = root1;
    TString stringout = rootout;

    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");

    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Merge(".mergedroot.root");

    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();

    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);

    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();

    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }

    newtree->Print();
    newtree->AutoSave();

    delete oldfile;
    delete newfile;
}

void mergesort(const char *root1 = 0, const char *root2 = 0, const char *rootout = 0) {
    
    // define strings
    TString string1 = root1;
    TString string2 = root2;
    TString stringout = rootout;
    
    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");
    string2.Append("/ntuple/ntuple");
    
    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Add(string2);
    fChain->Merge(".mergedroot.root");
    
    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();
    
    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();
    
    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }
    
    newtree->Print();
    newtree->AutoSave();
    
    delete oldfile;
    delete newfile;
}

void mergesort(const char *root1 = 0, const char *root2 = 0, const char *root3 = 0, const char *rootout = 0) {
    
    // define strings
    TString string1 = root1;
    TString string2 = root2;
    TString string3 = root3;
    TString stringout = rootout;
    
    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");
    string2.Append("/ntuple/ntuple");
    string3.Append("/ntuple/ntuple");
    
    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Add(string2);
    fChain->Add(string3);
    fChain->Merge(".mergedroot.root");
    
    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();
    
    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();
    
    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }
    
    newtree->Print();
    newtree->AutoSave();
    
    delete oldfile;
    delete newfile;
}

void mergesort(const char *root1 = 0, const char *root2 = 0, const char *root3 = 0, const char *root4 = 0, const char *rootout = 0) {
    
    // define strings
    TString string1 = root1;
    TString string2 = root2;
    TString string3 = root3;
    TString string4 = root4;
    TString stringout = rootout;
    
    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");
    string2.Append("/ntuple/ntuple");
    string3.Append("/ntuple/ntuple");
    string4.Append("/ntuple/ntuple");
    
    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Add(string2);
    fChain->Add(string3);
    fChain->Add(string4);
    fChain->Merge(".mergedroot.root");
    
    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();
    
    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();
    
    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }
    
    newtree->Print();
    newtree->AutoSave();
    
    delete oldfile;
    delete newfile;
}


void mergesort(const char *root1 = 0, const char *root2 = 0, const char *root3 = 0, const char *root4 = 0, const char *root5 = 0, const char *rootout = 0) {
    
    // define strings
    TString string1 = root1;
    TString string2 = root2;
    TString string3 = root3;
    TString string4 = root4;
    TString string5 = root5;
    TString stringout = rootout;
    
    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");
    string2.Append("/ntuple/ntuple");
    string3.Append("/ntuple/ntuple");
    string4.Append("/ntuple/ntuple");
    string5.Append("/ntuple/ntuple");
    
    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Add(string2);
    fChain->Add(string3);
    fChain->Add(string4);
    fChain->Add(string5);
    fChain->Merge(".mergedroot.root");
    
    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();
    
    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();
    
    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }
    
    newtree->Print();
    newtree->AutoSave();
    
    delete oldfile;
    delete newfile;
}


void mergesort(const char *root1 = 0, const char *root2 = 0, const char *root3 = 0, const char *root4 = 0, const char *root5 = 0, const char *root6 = 0, const char *rootout = 0) {
    
    // define strings
    TString string1 = root1;
    TString string2 = root2;
    TString string3 = root3;
    TString string4 = root4;
    TString string5 = root5;
    TString string6 = root6;
    TString stringout = rootout;
    
    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");
    string2.Append("/ntuple/ntuple");
    string3.Append("/ntuple/ntuple");
    string4.Append("/ntuple/ntuple");
    string5.Append("/ntuple/ntuple");
    string6.Append("/ntuple/ntuple");
    
    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Add(string2);
    fChain->Add(string3);
    fChain->Add(string4);
    fChain->Add(string5);
    fChain->Add(string6);
    fChain->Merge(".mergedroot.root");
    
    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();
    
    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();
    
    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }
    
    newtree->Print();
    newtree->AutoSave();
    
    delete oldfile;
    delete newfile;
}

void mergesort(const char *root1 = 0, const char *root2 = 0, const char *root3 = 0, const char *root4 = 0, const char *root5 = 0, const char *root6 = 0, const char *root7 = 0, const char *rootout = 0) {
    
    // define strings
    TString string1 = root1;
    TString string2 = root2;
    TString string3 = root3;
    TString string4 = root4;
    TString string5 = root5;
    TString string6 = root6;
    TString string7 = root7;
    TString stringout = rootout;
    
    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");
    string2.Append("/ntuple/ntuple");
    string3.Append("/ntuple/ntuple");
    string4.Append("/ntuple/ntuple");
    string5.Append("/ntuple/ntuple");
    string6.Append("/ntuple/ntuple");
    string7.Append("/ntuple/ntuple");
    
    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Add(string2);
    fChain->Add(string3);
    fChain->Add(string4);
    fChain->Add(string5);
    fChain->Add(string6);
    fChain->Add(string7);
    fChain->Merge(".mergedroot.root");
    
    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();
    
    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();
    
    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }
    
    newtree->Print();
    newtree->AutoSave();
    
    delete oldfile;
    delete newfile;
}

void mergesort(const char *root1 = 0, const char *root2 = 0, const char *root3 = 0, const char *root4 = 0, const char *root5 = 0, const char *root6 = 0, const char *root7 = 0, const char *root8 = 0, const char *rootout = 0) {
    
    // define strings
    TString string1 = root1;
    TString string2 = root2;
    TString string3 = root3;
    TString string4 = root4;
    TString string5 = root5;
    TString string6 = root6;
    TString string7 = root7;
    TString string8 = root8;
    TString stringout = rootout;
    
    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");
    string2.Append("/ntuple/ntuple");
    string3.Append("/ntuple/ntuple");
    string4.Append("/ntuple/ntuple");
    string5.Append("/ntuple/ntuple");
    string6.Append("/ntuple/ntuple");
    string7.Append("/ntuple/ntuple");
    string8.Append("/ntuple/ntuple");
    
    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Add(string2);
    fChain->Add(string3);
    fChain->Add(string4);
    fChain->Add(string5);
    fChain->Add(string6);
    fChain->Add(string7);
    fChain->Add(string8);
    fChain->Merge(".mergedroot.root");
    
    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();
    
    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();
    
    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }
    
    newtree->Print();
    newtree->AutoSave();
    
    delete oldfile;
    delete newfile;
}


void mergesort(const char *root1 = 0, const char *root2 = 0, const char *root3 = 0, const char *root4 = 0, const char *root5 = 0, const char *root6 = 0, const char *root7 = 0, const char *root8 = 0, const char *root9 = 0, const char *rootout = 0) {
    
    // define strings
    TString string1 = root1;
    TString string2 = root2;
    TString string3 = root3;
    TString string4 = root4;
    TString string5 = root5;
    TString string6 = root6;
    TString string7 = root7;
    TString string8 = root8;
    TString string9 = root9;
    TString stringout = rootout;
    
    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");
    string2.Append("/ntuple/ntuple");
    string3.Append("/ntuple/ntuple");
    string4.Append("/ntuple/ntuple");
    string5.Append("/ntuple/ntuple");
    string6.Append("/ntuple/ntuple");
    string7.Append("/ntuple/ntuple");
    string8.Append("/ntuple/ntuple");
    string9.Append("/ntuple/ntuple");
    
    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Add(string2);
    fChain->Add(string3);
    fChain->Add(string4);
    fChain->Add(string5);
    fChain->Add(string6);
    fChain->Add(string7);
    fChain->Add(string8);
    fChain->Add(string9);
    fChain->Merge(".mergedroot.root");
    
    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();
    
    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();
    
    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }
    
    newtree->Print();
    newtree->AutoSave();
    
    delete oldfile;
    delete newfile;
}


void mergesort(const char *root1 = 0, const char *root2 = 0, const char *root3 = 0, const char *root4 = 0, const char *root5 = 0, const char *root6 = 0, const char *root7 = 0, const char *root8 = 0, const char *root9 = 0, const char *root10 = 0, const char *rootout = 0) {
    
    // define strings
    TString string1 = root1;
    TString string2 = root2;
    TString string3 = root3;
    TString string4 = root4;
    TString string5 = root5;
    TString string6 = root6;
    TString string7 = root7;
    TString string8 = root8;
    TString string9 = root9;
    TString string10 = root10;
    TString stringout = rootout;
    
    // the ntuple in the g4out.root Geant4 files are in a directory "ntuple"
    string1.Append("/ntuple/ntuple");
    string2.Append("/ntuple/ntuple");
    string3.Append("/ntuple/ntuple");
    string4.Append("/ntuple/ntuple");
    string5.Append("/ntuple/ntuple");
    string6.Append("/ntuple/ntuple");
    string7.Append("/ntuple/ntuple");
    string8.Append("/ntuple/ntuple");
    string9.Append("/ntuple/ntuple");
    string10.Append("/ntuple/ntuple");
    
    // we call the new merged ntuple mntuple, and add all ntuples.
    TChain *fChain=new TChain("mntuple");
    fChain->Add(string1);
    fChain->Add(string2);
    fChain->Add(string3);
    fChain->Add(string4);
    fChain->Add(string5);
    fChain->Add(string6);
    fChain->Add(string7);
    fChain->Add(string8);
    fChain->Add(string9);
    fChain->Add(string10);
    fChain->Merge(".mergedroot.root");
    
    // get the merged root file
    TFile *oldfile = new TFile(".mergedroot.root");
    TTree *oldtree = (TTree*)oldfile->Get("mntuple");
    Long64_t nentries = oldtree->GetEntries();
    
    // create a new root file
    TFile *newfile = new TFile(stringout,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    // build the index from the eventNumber
    oldtree->BuildIndex("eventNumber");
    TTreeIndex *index = (TTreeIndex*)oldtree->GetTreeIndex();
    
    // loop over this index and fill the new tree
    for( int i = 0; i < index->GetN() ; ++i ) {
        oldtree->GetEntry(index->GetIndex()[i]);
        newtree->Fill();
    }
    
    newtree->Print();
    newtree->AutoSave();
    
    delete oldfile;
    delete newfile;
}
