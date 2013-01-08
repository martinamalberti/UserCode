{
  gROOT->SetStyle("Plain");
  _file0->cd();

  TH1F *hs0 = new TH1F("hs0","hs0",10,-1,1);
  TH1F *hs2 = new TH1F("hs2","hs2",20,-1,1);
  TH1F *hb = new TH1F("hb","hb",20,-1,1);
  TH1F *hbmc = new TH1F("hbmc","hbmc",3,0,3.15);

//   TCut sig_tight = "(leadPtOverM>0.5 && subleadPtOverM>0.3 && leadJPt>30 && subleadJPt>30 && TMath::Abs(deltaEtaJJ)>3 && TMath::Abs(Zep)<2.5 && MJJ>500 && TMath::Abs(deltaPhiJJGamGam)>2.6)";
//   TCut bkg_tight = "(abs(mass-125.)>5 && leadPtOverM>0.5 && subleadPtOverM>0.3 && leadJPt>30 && subleadJPt>30 && TMath::Abs(deltaEtaJJ)>3 && TMath::Abs(Zep)<2.5 && MJJ>500 && TMath::Abs(deltaPhiJJGamGam)>2.6)";
//   TCut sig_loose = "leadPtOverM>0.5 && subleadPtOverM>0.3 && leadJPt>30 && subleadJPt>20 && TMath::Abs(deltaEtaJJ)>3 && TMath::Abs(Zep)<2.5 && MJJ>250 && TMath::Abs(deltaPhiJJGamGam)>2.6 && !(leadPtOverM>0.5 && subleadPtOverM>0.3 && leadJPt>30 && subleadJPt>30 && TMath::Abs(deltaEtaJJ)>3 && TMath::Abs(Zep)<2.5 && MJJ>500 && TMath::Abs(deltaPhiJJGamGam)>2.6)";
//   TCut bkg_loose = "abs(mass-125.)>5 && leadPtOverM>0.5 && subleadPtOverM>0.3 && leadJPt>30 && subleadJPt>20 && TMath::Abs(deltaEtaJJ)>3 && TMath::Abs(Zep)<2.5 && MJJ>250 && TMath::Abs(deltaPhiJJGamGam)>2.6 && !(leadPtOverM>0.5 && subleadPtOverM>0.3 && leadJPt>30 && subleadJPt>30 && TMath::Abs(deltaEtaJJ)>3 && TMath::Abs(Zep)<2.5 && MJJ>500 && TMath::Abs(deltaPhiJJGamGam)>2.6)";

 
  TCut sig_tight = "diphotonMVA>-0.05 && category==4 && leadJPt>30 && subleadJPt>20 && MJJ>250";
  TCut bkg_tight = "(abs(MGamGam-125.)>5 && diphotonMVA>-0.05 && category==4 && leadJPt>30 && subleadJPt>20 && MJJ>250)";
 
  vbf_m125_8TeV->Draw("VBFSpin_Discriminant>>hs0", sig_tight,"GOFF");
  vbf_spin2_m125_8TeV->Draw("VBFSpin_Discriminant>>hs2", sig_tight,"GOFF");
  Data->Draw("VBFSpin_Discriminant>>hb", bkg_tight,"GOFF");
  
  hs0->SetLineColor(2);
  hs2->SetLineColor(4);
  hs0->Sumw2();
  hs2->Sumw2();
  hb->Sumw2();


  new TCanvas();
  //  hs0->GetXaxis()->SetTitle("#Delta#Phi(jj)");
  hs0->GetXaxis()->SetTitle("LD");
  hs0->GetYaxis()->SetTitle("a.u.");
  hs0->GetYaxis()->SetRangeUser(0.,2* hs0->GetMaximum());
  hs0->DrawNormalized();
  hs2->DrawNormalized("same");
  hb->DrawNormalized("same");
 
  TLegend *l = new TLegend(0.15,0.6,0.45,0.85);
  l->AddEntry(hs0,"MC VBF spin0","PL");
  l->AddEntry(hs2,"MC VBF spin2","PL");
  l->AddEntry(hb,"DATA |m(#gamma#gamma)-125.|>5","PL");
  l->Draw("same");
  
  
}
