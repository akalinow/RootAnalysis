{
  //This make the TSelectors in the library available to the remote proof session
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  gSystem->Load("libPFAnalysesVBFHTauTau");

gSystem->Setenv("PATH",gSystem->Getenv("PATH1"));

//gSystem->Exec("which python");
//gSystem->Exec("python -V");
}
