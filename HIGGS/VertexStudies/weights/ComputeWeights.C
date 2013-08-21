//****** simple macro to compute PU weights ******

void ComputeWeights(){
  
  //*** mc pileup file
  //  TFile *fmc = TFile::Open("./PU_DYJetsToLL_Summer12_DR53X-PU_S10.root"); // 2012
  //TFile *fmc = TFile::Open("./PU_DYJetsToLL_Summer12_DR53X-PU_RD1.root"); // 2012 RD ReReco
  TFile *fmc = TFile::Open("./PU_DYJetsToLL_Summer11dr53X-PU_S13.root"); // 2011 53X ReReco
  TH1F *hmc  = (TH1F*)fmc->Get("hmc");
  TH1F *hmctrue  = (TH1F*)fmc->Get("hmctrue");
  cout << "Number of bins in MC histogram: " << hmc->GetNbinsX()<< endl;

  //*** data file -- observed PU
  //  TFile *fda   = TFile::Open("../pileup/pileup_2012D_minBiasXsec69400_corr_observed.root");
  //TFile *fda   = TFile::Open("../pileup/pileup_2012ABCD_22Jan2013ReReco_corr_observed.root");
  TFile *fda   = TFile::Open("../pileup/pileup_2011_minBiasXsec68000_pixelcorr_observed.root");
  TH1F  *hdata = (TH1F*)fda->Get("pileup");
  cout << "Number of bins in DATA histogram: " << hdata->GetNbinsX()<< endl;
  
  //*** compute weights
  TH1F *hweights = (TH1F*)hdata->Clone("hweights");
  hweights->Divide(hdata,hmc,1./hdata->GetSumOfWeights(),1./hmc->GetSumOfWeights());

  //  TFile *fout = new TFile("./PUweights_DYJetsToLL_Summer12_DR53X-PU_RD1_minBiasXsec69400_corr_observed_2012ABCD.root","create");
  TFile *fout = new TFile("./PUweights_DYJetsToLL_Summer11dr53X-PU_S13_minBiasXsec68000_pixelcorr_observed_2011.root","create");
  hweights->Write("hweights");
  hdata->Write("hdata");
  hmc->Write("hmc");
  fout->Close();


//   //*** data file -- true PU
//   TFile *fdatrue = TFile::Open("../pileup/pileup_190456-208686_minBiasXsec69400_corr_observed.root");
//   TH1F *hdatatrue = (TH1F*)fdatrue->Get("pileup");
//   cout << hdatatrue->GetNbinsX()<< endl;

//   //*** compute weights
//   TH1F *hweightstrue = (TH1F*)hdatatrue->Clone("hweightstrue");
//   hweightstrue->Divide(hdatatrue,hmctrue,1./hdatatrue->GetSumOfWeights(),1./hmctrue->GetSumOfWeights());
 
//   TFile *fout2 = new TFile("./PUweights_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_true_Run2012ABCD.root","create");
//   hweightstrue->Write("hweights");
//   hdatatrue->Write("hdata");
//   hmctrue->Write("hmc");
//   fout2->Close();
}
