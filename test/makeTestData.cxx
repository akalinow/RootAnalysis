///This is a simplae script for making a test data file.
///This script should be executem form ROOT command line:
/// > root
/// [0] .x makeTestData.cxx
///
void makeTestData(){

  ////Data members
  Double_t x,y;

  ///Random number generator
  TRandom3 myRndm;
  
  // Define the Tree
  TFile *file = new TFile("Data.root","RECREATE");
  TTree *tree = new TTree("Data","Fakse events");

  TBranch *branchX = tree->Branch("x",&x);
  TBranch *branchY = tree->Branch("y",&y);

  Int_t nEvents = 1E8;

  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent) {

    x = myRndm.Gaus();
    y = myRndm.Gaus();
    
    //Fill the TTree with current data
    tree->Fill();
  }

  file->Write();
  delete file;
  
}
