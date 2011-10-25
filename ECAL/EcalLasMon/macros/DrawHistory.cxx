{
 
  gStyle->SetOptStat(0);
  gStyle->SetNdivisions(1020,"X");
  
  TFile *f = TFile::Open("EElaserAnalysis_fed605.root");

  // time interval
  TTimeStamp dateMin(2011, 2, 22, 0, kTRUE, 0); 
  TTimeStamp dateMax(2011, 10,20, 0, kTRUE, 0);   

  // plots
  TProfile *laser1 = f->Get("p_vptpn_las_fed605_harness6");
  TProfile *laser2 = f->Get("p_vptpn_led_fed605_harness6");
  TProfile *ratio  = f->Get("p_ratio_fed605_harness6");

  TProfile *tmax1  = f->Get("p_tmax_las_fed605_harness6");
  TProfile *tmax2  = f->Get("p_tmax_led_fed605_harness6");

  TProfile *matacqampli1     = f->Get("p_matacqampli_las_fed605_harness6");
  TProfile *matacqfwhm1      = f->Get("p_matacqfwhm_las_fed605_harness6");
  TProfile *matacqrisetime1  = f->Get("p_matacqrisetime_las_fed605_harness6");
  
  float ymin, ymax;

  TCanvas *clas = new TCanvas("clas","clas",1200,400);
  ymin = laser1->GetMinimum(0.8);
  ymax = laser1->GetMaximum(1.1);
  clas->SetGridx();
  clas->SetGridy();
  laser1->GetXaxis()->SetTitle("date");
  laser1->GetYaxis()->SetTitle("APD/PN");
  laser1->GetXaxis()->SetTimeFormat("%d/%m%F1970-01-01 00:00:00");
  laser1->GetXaxis()->SetTimeDisplay(1);
  laser1->GetXaxis()->SetRangeUser(dateMin, dateMax);
  laser1->GetXaxis()->SetNdivisions(1020);
  laser1->GetYaxis()->SetRangeUser(ymin,ymax);  
  laser1->Draw();
  laser2->Draw("same");
  ratio->Draw("same");
  TLegend *leg1 = new TLegend(0.12,0.7,0.3,0.89);
  leg1->SetFillColor(0);
  leg1->AddEntry(laser1, "blue laser", "LP");
  leg1->AddEntry(laser2, "blue led", "LP");
  leg1->AddEntry(ratio,  "ratio", "LP");
  leg1->Draw("same");

  TCanvas *ctmax = new TCanvas("ctmax","ctmax",1200,400);
  ymin = tmax1->GetMinimum(0.5);
  ymax = tmax1->GetMaximum(1.5);
  ctmax->SetGridx();
  ctmax->SetGridy();
  tmax1->GetXaxis()->SetNdivisions(1020);
  tmax1->GetXaxis()->SetTitle("date");
  tmax1->GetYaxis()->SetTitle("Tmax");
  tmax1->GetXaxis()->SetTimeDisplay(1);
  tmax1->GetXaxis()->SetRangeUser(dateMin, dateMax);
  tmax1->GetYaxis()->SetRangeUser(ymin,ymax);
  tmax1->Draw();
  tmax2->Draw("same");

  TCanvas *cparam = new TCanvas("cparam","cparam",1200,400);
  ymin = matacqampli1->GetMinimum(0);
  ymax = matacqampli1->GetMaximum(3);
  cparam->SetGridx();
  cparam->SetGridy();
  matacqampli1->GetXaxis()->SetNdivisions(1020);
  matacqampli1->GetXaxis()->SetTitle("date");
  matacqampli1->GetYaxis()->SetTitle("A/A_{0}");
  matacqampli1->GetXaxis()->SetTimeFormat("%d/%m%F1970-01-01 00:00:00");
  matacqampli1->GetXaxis()->SetTimeDisplay(1);
  matacqampli1->GetXaxis()->SetRangeUser(dateMin, dateMax);
  matacqampli1->GetYaxis()->SetRangeUser(ymin,ymax);
  matacqampli1->SetMarkerColor(kGreen+2);
  matacqampli1->SetMarkerStyle(24);
  matacqampli1->Draw();
  matacqfwhm1->SetMarkerColor(kRed);
  matacqfwhm1->SetMarkerStyle(24);
  matacqfwhm1->Draw("same");
  matacqrisetime1->SetMarkerColor(kMagenta+2);
  matacqrisetime1->SetMarkerStyle(24);
  matacqrisetime1->Draw("same");
  tmax1->Draw("same");

  TLegend *leg = new TLegend(0.12,0.7,0.3,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(matacqampli1, "matacq amplitude", "P");
  leg->AddEntry(matacqfwhm1, "matacq fwhm", "P");
  leg->AddEntry(matacqrisetime1, "matacq risetime", "P");
  leg->AddEntry(tmax1, "tmax", "P");
  leg->Draw("same");


  // 


}
