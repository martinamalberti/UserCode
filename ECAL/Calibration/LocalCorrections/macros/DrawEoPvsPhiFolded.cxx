{
  gROOT->SetStyle("Plain");
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.05);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.05);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.05);
  gROOT->ForceStyle();


  TFile *f1 = TFile::Open("outputFiles/Graphs_EoP_vs_LocalPhi_Folded20_crackCorr.root");
  TGraphErrors* g_EoP_MC_crack = (TGraphErrors*)f1->Get("gEoP_MC");
  TGraphErrors* g_EoP_Data_crack = (TGraphErrors*)f1->Get("gEoP_Data");
  g_EoP_MC_crack -> SetMarkerStyle(20);
  g_EoP_MC_crack -> SetMarkerSize(1);
  g_EoP_MC_crack -> SetMarkerColor(kRed); 
  g_EoP_MC_crack -> SetLineColor(kRed); 
  g_EoP_MC_crack -> SetLineStyle(2); 
  g_EoP_MC_crack -> SetLineWidth(2); 
  g_EoP_Data_crack -> SetMarkerStyle(20);
  g_EoP_Data_crack -> SetMarkerSize(1);
  g_EoP_Data_crack -> SetMarkerColor(kGreen); 
  g_EoP_Data_crack -> SetLineColor(kGreen); 
  g_EoP_Data_crack -> SetLineStyle(2); 
  g_EoP_Data_crack -> SetLineWidth(2); 

  bool drackCrackUnfolded = false;
  
  // TFile *f = TFile::Open("outputFiles/Graphs_EoP_vs_LocalPhi_Folded20_regression.root");
  //TFile *f = TFile::Open("outputFiles/Graphs_EoP_vs_LocalPhi_Folded4_regression.root");
  TFile *f = TFile::Open("test_fold4.root");
  
  int cryFold = 4;  

  TGraphErrors* g_EoP_MC   = (TGraphErrors*)f->Get("gEoP_MC");
  TGraphErrors* g_EoC_MC   = (TGraphErrors*)f->Get("gEoC_MC");  
  TGraphErrors* g_EoP_Data = (TGraphErrors*)f->Get("gEoP_Data");
  TGraphErrors* g_EoC_Data = (TGraphErrors*)f->Get("gEoC_Data");

  TGraphErrors* g_Ratio_MC   = (TGraphErrors*)f->Get("gCorr_Uncorr_MC");
  g_Ratio_MC -> SetMarkerStyle(20);
  g_Ratio_MC -> SetMarkerSize(0.7);
  g_Ratio_MC -> SetMarkerColor(kBlue); 
  g_Ratio_MC -> SetLineColor(kBlue); 

  TGraphErrors* g_Ratio_Data = (TGraphErrors*)f->Get("gCorr_Uncorr_Data");
  g_Ratio_Data -> SetMarkerStyle(20);
  g_Ratio_Data -> SetMarkerSize(0.7);
  g_Ratio_Data -> SetMarkerColor(kBlue); 
  g_Ratio_Data -> SetLineColor(kBlue); 

  // --- MC
  TCanvas* c_g_fit_MC = new TCanvas("g_fit_MC", "g_fit_MC",100,100,700,600);
  c_g_fit_MC->Divide(1,2);
  c_g_fit_MC->cd(1);

  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(0,0.96,cryFold,1.02);
  hPad->GetXaxis()->SetTitle("#phi_{SC} (deg)");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.045);

  g_EoP_MC -> SetMarkerStyle(20);
  g_EoP_MC -> SetMarkerSize(1);
  g_EoP_MC -> SetMarkerColor(kRed); 
  g_EoP_MC -> SetLineColor(kRed); 

  g_EoC_MC -> SetMarkerStyle(20);
  g_EoC_MC -> SetMarkerSize(1);
  g_EoC_MC -> SetMarkerColor(kRed+2); 
  g_EoC_MC -> SetLineColor(kRed+2); 

  if (drackCrackUnfolded) g_EoP_MC_crack -> Draw("L");
  g_EoP_MC -> Draw("PL");
  g_EoC_MC -> Draw("PL");
  
  TLegend *tl_MC = new TLegend(0.60,0.15,0.89,0.35);
  tl_MC -> SetFillColor(0);
  if (drackCrackUnfolded)
    tl_MC -> AddEntry(g_EoP_MC_crack,"MC - unfolded crack corr","L");  
  tl_MC -> AddEntry(g_EoP_MC,"MC - default","P");
  tl_MC -> AddEntry(g_EoC_MC,"MC - w/ regression","P");
  tl_MC -> Draw();


  // ratio MC
  c_g_fit_MC->cd(2);
  gPad->SetGrid();
  TH1F *hPad2 = (TH1F*)gPad->DrawFrame(0,0.98,cryFold,1.02);
  hPad2->GetXaxis()->SetTitle("#phi_{SC} (deg)");
  hPad2->GetYaxis()->SetTitle("ratio"); 
  hPad2->GetXaxis()->SetTitleOffset(0.9);
  hPad2->GetYaxis()->SetTitleSize(0.05);
  hPad2->GetYaxis()->SetLabelSize(0.045);

  g_Ratio_MC->Draw("PL");  






  //--- data
  TCanvas* c_g_fit_Data = new TCanvas("g_fit_Data", "g_fit_Data",100,100,700,600);
  c_g_fit_Data->Divide(1,2);
  c_g_fit_Data->cd(1);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(0,0.96,cryFold,1.02);
  hPad->GetXaxis()->SetTitle("#phi_{SC} (deg)");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.045);

  g_EoP_Data -> SetMarkerStyle(20);
  g_EoP_Data -> SetMarkerSize(1);
  g_EoP_Data -> SetMarkerColor(kGreen+1); 
  g_EoP_Data -> SetLineColor(kGreen+1); 

  g_EoC_Data -> SetMarkerStyle(20);
  g_EoC_Data -> SetMarkerSize(1);
  g_EoC_Data -> SetMarkerColor(kGreen+3); 
  g_EoC_Data -> SetLineColor(kGreen+3); 

  if (drackCrackUnfolded) g_EoP_Data_crack -> Draw("L");
  g_EoP_Data -> Draw("PL");
  g_EoC_Data -> Draw("PL");
  
  TLegend *tl_Data = new TLegend(0.60,0.15,0.89,0.35);
  tl_Data -> SetFillColor(0);
  if (drackCrackUnfolded)
    tl_Data -> AddEntry(g_EoP_Data_crack,"DATA - unfolded crack corr","L");  
  tl_Data -> AddEntry(g_EoP_Data,"DATA - default","P");
  tl_Data -> AddEntry(g_EoC_Data,"DATA - w/regression","P");
  tl_Data ->Draw();

  // ratio MC
  c_g_fit_Data->cd(2);
  gPad->SetGrid();
  TH1F *hPad2 = (TH1F*)gPad->DrawFrame(0,0.98,cryFold,1.02);
  hPad2->GetXaxis()->SetTitle("#phi_{SC} (deg)");
  hPad2->GetYaxis()->SetTitle("ratio"); 
  hPad2->GetXaxis()->SetTitleOffset(0.9);
  hPad2->GetYaxis()->SetTitleSize(0.05);
  hPad2->GetYaxis()->SetLabelSize(0.045);

  g_Ratio_Data->Draw("PL");  


  // residual spread data wrt MC
  TH1F* h_spreadRatio = new TH1F("h_spreadRatio", "spreadRatio", 200, 0.95, 1.05);
  for (int i = 1; i < g_EoC_Data -> GetN()+1; i++){
    double xda, yda;
    double  xmc, ymc;
    g_EoC_MC   -> GetPoint(i,xmc,ymc);
    g_EoC_Data -> GetPoint(i,xda,yda);
    cout << yda << "  " << ymc << endl;
    h_spreadRatio->Fill(yda/ymc);
  }

  TCanvas *cspread = new TCanvas("cspread","cspread",500,500);
  h_spreadRatio->Draw();
  

}
