

{
  
  TFile *fz = TFile::Open("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");

  TH1F *hz = new TH1F("hz","hz",100, -2.5, 2.5);  
  TH1F *hzHighR9 = new TH1F("hzHighR9","hzHighR9",100, -2.5, 2.5);  
  TH1F *hzLowR9 = new TH1F("hzLowR9","hzLowR9",100, -2.5, 2.5);
  
  fz->cd();
  ntu->Draw("ele1_scEta >> hz","ele1_e3x3/ele1_scE > 0","goff");
  ntu->Draw("ele1_scEta >> hzHighR9","ele1_e3x3/ele1_scE > 0.90","goff");
  ntu->Draw("ele1_scEta >> hzLowR9","ele1_e3x3/ele1_scE < 0.90","goff");
  hzHighR9 -> Divide(hzHighR9,hz);
  hzLowR9  -> Divide(hzLowR9,hz);


  TFile *fw = TFile::Open("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root"); 

  TH1F *hw = new TH1F("hw","hw",100, -2.5, 2.5);  
  TH1F *hwHighR9 = new TH1F("hwHighR9","hwHighR9",100, -2.5, 2.5);  
  TH1F *hwLowR9 = new TH1F("hwLowR9","hwLowR9",100, -2.5, 2.5);  

  fw->cd();
  ntu->Draw("ele1_scEta >> hw","ele1_e3x3/ele1_scE > 0","goff");
  ntu->Draw("ele1_scEta >> hwHighR9","ele1_e3x3/ele1_scE > 0.90","goff");
  ntu->Draw("ele1_scEta >> hwLowR9","ele1_e3x3/ele1_scE < 0.90","goff");
  hwHighR9->Divide(hwHighR9,hw);
  hwLowR9->Divide(hwLowR9,hw);

  TCanvas *c = new TCanvas("c","c");
  hzHighR9 ->SetLineColor(2);
  hzHighR9 ->GetYaxis()->SetRangeUser(0,1.1);
  hzHighR9 -> Draw("");
  hwHighR9 -> Draw("same");

  TCanvas *cc = new TCanvas("cc","cc");
  hzLowR9 ->SetLineColor(2);
  hzLowR9 ->GetYaxis()->SetRangeUser(0,1.1);
  hzLowR9 -> Draw("");
  hwLowR9 -> Draw("same");

  
    
  

  
}
