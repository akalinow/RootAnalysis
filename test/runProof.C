{

TProof::AddEnvVar("PATH1",gSystem->Getenv("PATH"));

TProof * p = TProof::Open("workers=7"); 
gSystem->Load("libRootAnalysis");

//This makes sure the TSelector library and dictionary are properly
// installed in the remote PROOF servers
p->Exec( ".x proof_remote.C" );

//This creates the 'data set' which defines what files we need to process
// NOTE: the files given must be accessible by the remote systems
TDSet c( "TTree", "Events");

TPython::LoadMacro("tmpConfig.py");
int size = TPython::Eval("len(process.source.fileNames)");
int nEventsToRead = TPython::Eval("process.maxEvents.input.value()");

char text[200]; 
for(int i=0;i<size;++i){
  sprintf(text,"process.source.fileNames[%d]",i);
  std::string fileName = TPython::Eval(text);
  std::string fileNameTmp = fileName.substr(fileName.find("file:")+5,fileName.size());
  c.Add(fileNameTmp.c_str());
}


TString configFile(gSystem->pwd());
configFile+="/tmpConfig.py";

//This makes the actual processing happen
p->Process(&c,"VBFTauTau::FWLiteTSelector",configFile, nEventsToRead);

p->Close();
}
