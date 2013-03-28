

{

  TFile *f0 = TFile::Open("vtxIdScaleFactorFromZmumu_dz1.root");
  TFile *f1 = TFile::Open("vtxIdScaleFactorFromZmumu_dz1_BSrwNew.root");
  TFile *f2 = TFile::Open("vtxIdScaleFactorFromZmumu_dz01.root");
  TFile *f3 = TFile::Open("vtxIdScaleFactorFromZmumu_dz01_BSrwNew.root");

  TLegend *l = new TLegend(0.6,0.6,0.89,0.89);

  f0->cd();
  hscaleFactor->GetXaxis()->SetRangeUser(0,200);
  hscaleFactor->Draw();
  l->AddEntry(hscaleFactor,"DZ < 1 - NO BS re-weighting");

  f1->cd();
  hscaleFactor->SetLineColor(2);
  hscaleFactor->SetMarkerColor(2);
  hscaleFactor->Draw("same");
  l->AddEntry(hscaleFactor,"DZ < 1 - WITH BS re-weighting");

  f2->cd();
  hscaleFactor->SetLineColor(1);
  hscaleFactor->SetMarkerColor(1);
  hscaleFactor->SetMarkerStyle(24);
  hscaleFactor->Draw("same");
  l->AddEntry(hscaleFactor,"DZ < 0.1 - NO BS re-weighting");

  f3->cd();
  hscaleFactor->SetLineColor(2);
  hscaleFactor->SetMarkerColor(2);
  hscaleFactor->SetMarkerStyle(24);
  hscaleFactor->Draw("same");
  l->AddEntry(hscaleFactor,"DZ < 0.1 - WITH BS re-weighting");

  l->Draw("same") ;




}
