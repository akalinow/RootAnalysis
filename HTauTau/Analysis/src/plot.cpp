void plot(){

TFile * file = new TFile("RootAnalysis_AnalysisMuTau.root");

TH2F *h2DH = (TH2F*)file->Get("HTTAnalyzer/h2DTestHistoggHTT125_mu_pi_RefitPV");
TH2F *h2DA = (TH2F*)file->Get("HTTAnalyzer/h2DTestHistoATT_mu_pi_RefitPV");
TH2F *h2DZ = (TH2F*)file->Get("HTTAnalyzer/h2DTestHistoDY0JetsMatchT_mu_pi_RefitPV");
//TH2F *h2DTestHistoggHTT125_66 = (TH2F*)file->Get("HTTAnalyzer/h2DTestHistoggHTT125_mu_pi_GenPV");


//TH1F *h1DPhiTestHistoggHTT125_66 = (TH1F*)file->Get("HTTAnalyzer/h1DPhiTestHistoggHTT125_mu_pi_RefitPV");
//TH1F *h1DPhiTestHistoggHTT125_66 = (TH1F*)file->Get("HTTAnalyzer/h1DPhi-nVectorsggHTT125_mu_pi_RefitPV");
//TH1F *h1DPhiTestHistoggHTT125_66 = (TH1F*)file->Get("HTTAnalyzer/h1DPhi-nVectorsATT_mu_pi_RefitPV");


TProfile *hProf = (TProfile*)file->Get("HTTAnalyzer/hProfTest1ggHTT125_mu_pi_RefitPV");
std::cout<<hProf<<std::endl;

TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
                                 600,500);

h2DZ->SetXTitle("#gamma_{leg1} + #gamma_{leg2} from SVFit");
h2DZ->SetYTitle("SVfit mass");
h2DZ->SetStats(kFALSE);

h2DH->SetMarkerColor(1);
h2DH->SetLineColor(1);

h2DH->SetStats(kFALSE);
h2DH->Draw("colz text");
h2DZ->Draw("candle");
h2DH->Draw("candle");
TF1 *f1 = new TF1("f1","125.09",0,140);
f1->SetLineColor(2);
f1->Draw("same");
//h2DTestHistoggHTT125_66->Draw("violin");
c->Print("plotH.png");

h2DZ->Draw("candle");
f1->Draw("same");
c->Print("plotZ.png");

h2DH->Draw("candle same");

c->Print("plotH_Z.png");
//h2DTestHistoggHTT125_66->Print("all");
//return;

hProf->Draw();
TGraph *aGr = new TGraph(hProf->GetNbinsX());
for(unsigned int iBin=1;iBin<hProf->GetNbinsX();++iBin){
    float x = hProf->GetXaxis()->GetBinCenter(iBin);
    float y = hProf->GetBinError(iBin);
    aGr->SetPoint(iBin,x,y);
}
aGr->Draw("AC");
//return;

TH1D *h1DH = h2DH->ProjectionX();
TH1D *h1DA = h2DA->ProjectionX();
TH1D *h1DZ = h2DZ->ProjectionX();

h1DA->SetLineColor(2);
h1DZ->SetLineColor(3);

h1DZ->SetStats(kFALSE);
h1DZ->DrawNormalized();
//h1DA->DrawNormalized("same");
h1DH->DrawNormalized("same");
c->Print("plot_projX.png");
//return;

int iBinLow = h2DH->GetXaxis()->FindBin(100.0);
int iBinMax = h2DH->GetXaxis()->FindBin(700.0);
h1DH = h2DH->ProjectionY("h1DH",iBinLow, iBinMax);
h1DA = h2DA->ProjectionY("h1DA",iBinLow, iBinMax);
h1DZ = h2DZ->ProjectionY("h1DZ",iBinLow, iBinMax);

h1DA->SetLineColor(2);
h1DZ->SetLineColor(3);

h1DZ->SetStats(kFALSE);
h1DZ->DrawNormalized();
//h1DA->DrawNormalized("same");
h1DH->DrawNormalized("same");
c->Print("plot_projY_70.png");

}
