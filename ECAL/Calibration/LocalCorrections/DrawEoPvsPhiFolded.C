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
  
  TFile *f = TFile::Open("Graphs_EoP_vs_LocalPhi_Folded20.root");
  
  int cryFold = 20;  

  TGraphErrors* g_EoP_MC   = (TGraphErrors*)f->Get("gEoP_MC");
  TGraphErrors* g_EoC_MC   = (TGraphErrors*)f->Get("gEoC_MC");  
  TGraphErrors* g_EoP_Data = (TGraphErrors*)f->Get("gEoP_Data");
  TGraphErrors* g_EoC_Data = (TGraphErrors*)f->Get("gEoC_Data");


  // --- MC
  TCanvas* c_g_fit_MC = new TCanvas("g_fit_MC", "g_fit_MC",100,100,800,500);
  c_g_fit_MC->cd(1);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(0,0.93,cryFold,1.03);
  hPad->GetXaxis()->SetTitle("#phi_{SC} (deg)");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.05);

  g_EoP_MC -> SetMarkerStyle(20);
  g_EoP_MC -> SetMarkerSize(1);
  g_EoP_MC -> SetMarkerColor(kRed); 
  g_EoP_MC -> Draw("PL");
  g_EoC_MC -> SetMarkerStyle(20);
  g_EoC_MC -> SetMarkerSize(1);
  g_EoC_MC -> SetMarkerColor(kRed+2); 
  g_EoC_MC -> Draw("PL");

//   g_EoC_Data -> SetMarkerStyle(20);
//   g_EoC_Data -> SetMarkerSize(1);
//   g_EoC_Data -> SetMarkerColor(kGreen+2); 
//   g_EoC_Data -> Draw("PL");

  TLegend *tl_MC = new TLegend(0.60,0.15,0.89,0.35);
  tl_MC -> SetFillColor(0);
  tl_MC -> AddEntry(g_EoP_MC,"MC - uncorrected","P");
  tl_MC -> AddEntry(g_EoC_MC,"MC - corrected","P");
  tl_MC -> Draw();



  //--- data
  TCanvas* c_g_fit_Data = new TCanvas("g_fit_Data", "g_fit_Data",100,100,800,500);
  c_g_fit_Data->cd(1);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(0,0.93,cryFold,1.03);
  hPad->GetXaxis()->SetTitle("#phi_{SC} (deg)");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.05);

  g_EoP_Data -> SetMarkerStyle(20);
  g_EoP_Data -> SetMarkerSize(1);
  g_EoP_Data -> SetMarkerColor(kGreen); 
  g_EoP_Data -> Draw("PL");
  g_EoC_Data -> SetMarkerStyle(20);
  g_EoC_Data -> SetMarkerSize(1);
  g_EoC_Data -> SetMarkerColor(kGreen+2); 
  g_EoC_Data -> Draw("PL");

  TLegend *tl_Data = new TLegend(0.60,0.15,0.89,0.35);
  tl_Data -> SetFillColor(0);
  tl_Data -> AddEntry(g_EoP_Data,"DATA - uncorrected","P");
  tl_Data -> AddEntry(g_EoC_Data,"DATA - corrected","P");
  tl_Data ->Draw();



}
