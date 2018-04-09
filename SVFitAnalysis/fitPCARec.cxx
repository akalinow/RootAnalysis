void fitPCARec(){

  TFile f("RootAnalysis_AnalysisMuTau.root");

  TH1F *h1DFlightPathPCARecggHTT125 = (TH1F*)f.Get("svfitAnalyzer/h1DFlightPathPCARecggHTT125");
  TH1F *h1DFlightPathPCARecDYAllJetsMatchT = (TH1F*)f.Get("svfitAnalyzer/h1DFlightPathPCARecDYAllJetsMatchT");

  h1DFlightPathPCARecDYAllJetsMatchT->SetLineColor(2);
/*
  h1DFlightPathPCARecggHTT125->DrawCopy();
  h1DFlightPathPCARecDYAllJetsMatchT->DrawCopy("same");
  return;
*/
  TH1F *h = (TH1F*)h1DFlightPathPCARecggHTT125->Clone("h");
  h->Add(h1DFlightPathPCARecDYAllJetsMatchT);
  h->Scale(1.0/2);

  TF1 * fitFunc = new TF1("fitFunc", "expo(0)",0,0.01);

  Double_t par[6];
   TF1 *f1    = new TF1("g1","expo(0)",0,0.006);
   TF1 *f2    = new TF1("g2","expo(0)",0.006,0.03);
   TF1 *f3    = new TF1("g3","expo(0)",0.04,0.10);
   TF1 *total = new TF1("total","expo(0)+expo(2)+expo(4)",0,0.1);
   total->SetLineColor(3);
   h->Fit(f1,"R");
   h->Fit(f2,"R+");
   h->Fit(f3,"R+");
   f1->GetParameters(&par[0]);
   f2->GetParameters(&par[2]);
   f3->GetParameters(&par[4]);
   total->SetParameters(par);
   h->Fit(total,"R+");
   h->DrawCopy();
   total->SavePrimitive(std::cout);

   TH1F *h1DFlightPathPCARecLeg1ggHTT125 = (TH1F*)f.Get("svfitAnalyzer/h1DFlightPathPCARecLeg1ggHTT125");
   TH1F *h1DFlightPathPCARecLeg1DYAllJetsMatchT= (TH1F*)f.Get("svfitAnalyzer/h1DFlightPathPCARecLeg1DYAllJetsMatchT");
   h1DFlightPathPCARecLeg1ggHTT125->Add(h1DFlightPathPCARecLeg1DYAllJetsMatchT);
/*
   h1DFlightPathPCARecLeg1DYAllJetsMatchT->SetLineColor(2);
   h1DFlightPathPCARecLeg1DYAllJetsMatchT->DrawCopy("same");
   return;
*/

   h->SetLineColor(1);

   h1DFlightPathPCARecLeg1ggHTT125->Scale(1.0/2);
   h1DFlightPathPCARecLeg1ggHTT125->DrawCopy();
   h->DrawCopy("same");
   total->Draw("same");
   h->Print();
   h1DFlightPathPCARecLeg1ggHTT125->Print();
}
