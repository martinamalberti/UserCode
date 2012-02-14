//****** simple macro to compute PU weights ******
{
  //*** mc file - Fall 2011
  TChain *ntu = new TChain("ntu");
  ntu->Add("/data2/calibrator/NTUPLES/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11_All.root");
  TH1F *hmc = new TH1F("hmc","hmc",60,0,60);
  ntu->Draw("PUit_NumInteractions >> hmc","","goff");
 
  //*** data file 
  TFile *fda = TFile::Open("/afs/cern.ch/user/a/adavidzh/public/json/111105/pileup/2011_0100_73500.pileup.root");

  TH1F *pileup = (TH1F*)fda->Get("pileup");

  TH1F *hdata = new TH1F("hdata","hdata",60,0,60);
  for (int ibin = 1; ibin < 61; ibin++){
    hdata->SetBinContent(ibin, pileup->GetBinContent(ibin));
  }
 

  //*** compute weights
  TH1F *hweights = (TH1F*)hdata->Clone("hweights");
  hweights->Divide(hdata,hmc,1./hdata->GetSumOfWeights(),1./hmc->GetSumOfWeights());
  TFile *fout = new TFile("./PUweights_2011_0100_73500_WJetsToLL_Fall11_S6_new.root","create");
  hweights->Write("hweights");

}
