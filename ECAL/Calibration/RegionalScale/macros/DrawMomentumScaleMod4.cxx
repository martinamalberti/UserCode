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

  //TFile *f = TFile::Open("outputFiles/momentumScale_mod4_EBM.root");
  TFile *f = TFile::Open("outputFiles/RelativeScaleVSVertexPosition_mod4_EBM_noPcalibration.root");
  //TFile *f = TFile::Open("outputFiles/zeeScale_mod4_EBP.root");

  g_EoP -> SetMarkerStyle(20);
  g_EoP -> SetMarkerSize(1);
  g_EoP -> SetMarkerColor(kRed+1); 
  g_EoP -> SetLineColor(kRed+1); 
  g_EoC -> SetMarkerStyle(20);
  g_EoC -> SetMarkerSize(1);
  g_EoC -> SetMarkerColor(kGreen+2);
  g_EoC -> SetLineColor(kGreen+2);
  
  g_Rat -> SetMarkerStyle(20);
  g_Rat -> SetMarkerSize(1);
  g_Rat -> SetMarkerColor(kBlue+2); 
  g_Rat -> SetLineColor(kBlue+2); 

  float zLim = 20;

  TCanvas* c_g_fit = new TCanvas("g_fit", "g_fit",100,100,800,600);
  TPad *cLower = new TPad("pad_0","pad_0",0.00,0.00,1.00,0.30);
  TPad *cUpper = new TPad("pad_1","pad_1",0.00,0.30,1.00,1.00);

  cLower->SetBottomMargin(0.2); 
  cUpper->SetTopMargin(0.03); 
  cUpper->SetBottomMargin(0.1); 

  cLower->Draw();
  cUpper->Draw();

  float FontSCF = cUpper->GetHNDC()/cLower->GetHNDC(); 
  float tYoffset = 0.9; 
  float labSize  = 0.05;

  cUpper -> cd();
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(-zLim,0.95,zLim,1.05);
  hPad->GetXaxis()->SetLabelSize(labSize);
  hPad->GetXaxis()->SetTitleSize(labSize);
  hPad->GetYaxis()->SetLabelSize(labSize);
  hPad->GetYaxis()->SetTitleSize(labSize);
  hPad->GetXaxis()->SetTitleOffset(tYoffset);
  hPad->GetYaxis()->SetTitleOffset(tYoffset);
  hPad->GetXaxis()->SetTitle("z_{vtx} (cm)");
  hPad->GetYaxis()->SetTitle("M_{Z}^{2}/M_{ee}^{2} #propto 1/p"); 
  
  g_EoP -> Draw("PL");
  g_EoC -> Draw("PL");

  TLegend *tl = new TLegend(0.70,0.80,0.89,0.95);
  tl -> SetTextFont(40);
  tl -> SetFillColor(0);
  //  tl -> SetBorderSize(0); 
  tl -> AddEntry(g_EoP,"MC","PL");
  tl -> AddEntry(g_EoC,"Data","PL");
  tl -> Draw();

  cLower -> cd();
  gPad->SetGrid();

  TH1F *hRat = (TH1F*)gPad->DrawFrame(-zLim,0.97,zLim,1.03);
  hRat->GetYaxis()->SetNdivisions(505);
  hRat->GetXaxis()->SetLabelSize(labSize*FontSCF);
  hRat->GetXaxis()->SetTitleSize(labSize*FontSCF);
  hRat->GetYaxis()->SetLabelSize(labSize*FontSCF);
  hRat->GetYaxis()->SetTitleSize(labSize*FontSCF);
  hRat->GetYaxis()->SetTitleOffset(tYoffset/FontSCF);
  hRat->GetXaxis()->SetTitleOffset(0.6);
  hRat->GetXaxis()->SetTitle("z_{vtx} (cm)"); 
  hRat->GetYaxis()->SetTitle("Data / MC"); 
  g_Rat->Draw("PL"); 





}
