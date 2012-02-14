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
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(10); 
  gStyle->SetOptFit(1); 
  gROOT->ForceStyle();

  TFile *f = TFile::Open("outputFiles/momentumScale_scM_highR9.root");

  g_EoP -> SetMarkerStyle(20);
  g_EoP -> SetMarkerSize(1);
  g_EoP -> SetMarkerColor(kRed+1); 
  g_EoP -> SetLineColor(kRed+1); 
  g_EoC -> SetMarkerStyle(20);
  g_EoC -> SetMarkerSize(1);
  g_EoC -> SetMarkerColor(kGreen+1);
  g_EoC -> SetLineColor(kGreen+1);
  
  g_Rat -> SetMarkerStyle(20);
  g_Rat -> SetMarkerSize(1);
  g_Rat -> SetMarkerColor(kBlue+2); 
  g_Rat -> SetLineColor(kBlue+2); 

  float etaLim = 115;

  TCanvas* c_g_fit = new TCanvas("g_fit", "g_fit",100,100,800,600);
  TPad *cLower = new TPad("pad_0","pad_0",0.00,0.00,1.00,0.30);
  TPad *cUpper = new TPad("pad_1","pad_1",0.00,0.30,1.00,1.00);

  cLower->SetBottomMargin(0.25); cUpper->SetTopMargin(0.01); 
  cUpper->SetBottomMargin(0.01); 

  cLower->Draw();
  cUpper->Draw();

  float FontSCF = cUpper->GetHNDC()/cLower->GetHNDC(); 
  float tYoffset = 0.8; 
  float labSize = 0.06;

  cUpper -> cd();
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(-etaLim,0.955,etaLim,1.025);
  hPad->GetXaxis()->SetLabelSize(labSize);
  hPad->GetXaxis()->SetTitleSize(labSize);
  hPad->GetYaxis()->SetLabelSize(labSize);
  hPad->GetYaxis()->SetTitleSize(labSize);
  hPad->GetXaxis()->SetTitleOffset(tYoffset);
  hPad->GetYaxis()->SetTitleOffset(tYoffset);
  hPad->GetXaxis()->SetTitle("#eta_{SC}");
  hPad->GetYaxis()->SetTitle("M_{Z}^{2}/M_{ee}^{2} #propto 1/p"); 
  
  g_EoP -> Draw("PL");
  g_EoC -> Draw("PL");

  TLegend *tl = new TLegend(0.80,0.80,0.90,0.95);
  tl -> SetTextFont(40);
  tl -> SetFillColor(0);
  //  tl -> SetBorderSize(0); 
  tl -> AddEntry(g_EoP,"MC","PL");
  tl -> AddEntry(g_EoC,"Data","PL");
  tl -> Draw();

  cLower -> cd();
  gPad->SetGrid();

  TH1F *hRat = (TH1F*)gPad->DrawFrame(-etaLim,0.98,etaLim,1.02);
  hRat->GetYaxis()->SetNdivisions(505);
  hRat->GetXaxis()->SetLabelSize(labSize*FontSCF);
  hRat->GetXaxis()->SetTitleSize(labSize*FontSCF);
  hRat->GetYaxis()->SetLabelSize(labSize*FontSCF);
  hRat->GetYaxis()->SetTitleSize(labSize*FontSCF);
  hRat->GetYaxis()->SetTitleOffset(tYoffset/FontSCF);
  hRat->GetXaxis()->SetTitleOffset(0.6);
  hRat->GetXaxis()->SetTitle("#eta_{SC}"); 
  hRat->GetYaxis()->SetTitle("Data / MC"); 
  g_Rat->Draw("PL"); 




}
