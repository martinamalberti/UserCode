{

  //*** mc ntuple 
  TChain* chain = new TChain("MiBiCommonNTTwoPhotons/SimpleNtuple");
  //chain->Add("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VERTEX/2012/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/MiBiCommonNT_*.root");
  //chain->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/vertextrees_martina/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/MiBiCommonNT_*.root");
  chain->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/vertextrees_martina/MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11dr53X-PU_S13_START53_LV3-v1_AODSIM/MiBiCommonNT_*.root");

 
  cout << "Observed PU ... "<< endl;
  TH1F *hmc = new TH1F("hmc","hmc",60,0,60);
  cout << hmc->GetNbinsX()<< endl;
  chain->Draw("mc_PUit_NumInteractions>>hmc");


  cout << "True PU ... "<< endl;
  TH1F *hmctrue = new TH1F("hmctrue","hmctrue",60,0,60);
  cout << hmctrue->GetNbinsX()<< endl;
  chain->Draw("mc_PUit_TrueNumInteractions>>hmctrue");

  //TFile *fout = new TFile("./PU_DYJetsToLL_Summer12_DR53X-PU_RD1.root","create");
  TFile *fout = new TFile("./PU_DYJetsToLL_Summer11dr53X-PU_S13.root","create");
  hmc->Write("hmc");
  hmctrue->Write("hmctrue");
  fout->Close();
  
}

