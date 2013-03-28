//****** simple macro to compute PU weights ******

void ComputeWeights(){
  
  //*** mc pileup file
  TFile *fmc = TFile::Open("./PU_DYJetsToLL_Summer12_DR53X-PU_S10.root"); 
  TH1F *hmc  = (TH1F*)fmc->Get("hmc");
  TH1F *hmctrue  = (TH1F*)fmc->Get("hmctrue");
  //cout << hmc->GetNbinsX()<< endl;

  //*** mc ntuple 
  //   TChain* chain = new TChain("MiBiCommonNTTwoPhotons/SimpleNtuple");
  //   chain->Add("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VERTEX/2012/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/MiBiCommonNT_*.root");

//   TH1F *hmc = new TH1F("hmc","hmc",60,0,60);
//   cout << hmc->GetNbinsX()<< endl;
//   chain->Draw("mc_PUit_NumInteractions>>hmc");

//   TH1F *hmctrue = new TH1F("hmctrue","hmctrue",60,0,60);
//   cout << hmctrue->GetNbinsX()<< endl;
//   chain->Draw("mc_PUit_TrueNumInteractions>>hmctrue");

  //*** data file -- observed PU
  TFile *fda   = TFile::Open("../pileup/pileup_2012D_minBiasXsec69400_corr_observed.root");
  TH1F  *hdata = (TH1F*)fda->Get("pileup");
  cout << hdata->GetNbinsX()<< endl;
  
  //*** compute weights
  TH1F *hweights = (TH1F*)hdata->Clone("hweights");
  hweights->Divide(hdata,hmc,1./hdata->GetSumOfWeights(),1./hmc->GetSumOfWeights());

  TFile *fout = new TFile("./PUweights_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012D.root","create");
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
