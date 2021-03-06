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

  float etaMin = -115;
  float etaMax = 115;
  float ymin   = 0.9;
  float ymax   = 1.1;  
  float yrmin   = 0.95;
  float yrmax   = 1.05;  
  
  //TFile *f = TFile::Open("outputFiles/GraphsRegionalScaleEta_highR9_diffTemplates.root");
  //TFile *f = TFile::Open("outputFiles/GraphsRegionalScaleEta_highR9_DYee.root");
  //TFile *f = TFile::Open("outputFiles/GraphsRegionalScaleEta_lowR9_diffTemplates.root");
  //TFile *f = TFile::Open("GraphsRegionalScaleEta_allR9_Pcalib_templatesFromRegressionFromMC.root");
  TFile *f = TFile::Open("outputFiles/GraphsRegionalScaleEta_allR9_Pcalib_templatesFromMC_42XReReco_FT_R_42_V21B.root");
  //TFile *f = TFile::Open("GraphsRegionalScaleEta_allR9_Pcalib_templatesFromMC_42XReReco_FT_R_42_V21B.root");
  
  TGraphErrors* g_EoP_MC   = (TGraphErrors*)f->Get("gEoP_MC");
  g_EoP_MC -> SetMarkerStyle(21);
  g_EoP_MC -> SetMarkerSize(0.7);
  g_EoP_MC -> SetMarkerColor(kRed); 

  TGraphErrors* g_EoC_MC   = (TGraphErrors*)f->Get("gEoC_MC");  
  g_EoC_MC -> SetMarkerStyle(21);
  g_EoC_MC -> SetMarkerSize(0.7);
  g_EoC_MC -> SetMarkerColor(kRed+2); 

  TGraphErrors* g_EoP_Data = (TGraphErrors*)f->Get("gEoP_Data");
  g_EoP_Data -> SetMarkerStyle(20);
  g_EoP_Data -> SetMarkerSize(0.7);
  g_EoP_Data -> SetMarkerColor(kGreen); 
  
  TGraphErrors* g_EoC_Data = (TGraphErrors*)f->Get("gEoC_Data");
  g_EoC_Data -> SetMarkerStyle(20);
  g_EoC_Data -> SetMarkerSize(0.7);
  g_EoC_Data -> SetMarkerColor(kGreen+2); 

  TGraphErrors* g_Ratio_MC   = (TGraphErrors*)f->Get("gCorr_Uncorr_MC");
  g_Ratio_MC -> SetMarkerStyle(21);
  g_Ratio_MC -> SetMarkerSize(0.7);
  g_Ratio_MC -> SetMarkerColor(kBlue); 

  TGraphErrors* g_Ratio_Data = (TGraphErrors*)f->Get("gCorr_Uncorr_Data");
  g_Ratio_Data -> SetMarkerStyle(20);
  g_Ratio_Data -> SetMarkerSize(0.7);
  g_Ratio_Data -> SetMarkerColor(kBlue); 

  TGraphErrors* gP_ratio = new TGraphErrors();
  for(int i =0 ; i < g_EoP_MC ->GetN();i++){
    double v1, v2, v3, xx1, xx2, xx3;
    g_EoP_Data->GetPoint(i,xx1,v1);
    g_EoP_MC->GetPoint(i,xx2,v2);
    if( xx1 != xx2 ){ cout<<"Error 2 !!"<<endl; return;}
    double rat = v1/v2; 
    double ev1 = g_EoP_Data->GetErrorY(i);
    double ev2 = g_EoP_MC->GetErrorY(i);
    
    if (ev1 > 5) {
      g_EoP_Data->SetPointError(i, 0, 0);
      ev1 = 0;
    }
    if (ev2 > 5) {
      g_EoP_MC->SetPointError(i, 0, 0);
      ev2 = 0;
    }

    double erat = rat*sqrt(  (ev1*ev1)/(v1*v1) + (ev2*ev2)/(v2*v2));
    gP_ratio->SetPoint(i, xx1, rat);
    gP_ratio->SetPointError(i, 0, erat);
  }
  
  TGraphErrors* gC_ratio = new TGraphErrors();
  for(int i =0 ; i < g_EoC_MC ->GetN();i++){
    double v1, v2, v3, xx1, xx2, xx3;
    g_EoC_Data->GetPoint(i,xx1,v1);
    g_EoC_MC->GetPoint(i,xx2,v2);
    if( xx1 != xx2 ){ cout<<"Error 2 !!"<<endl; return;}
    double rat = v1/v2; 
    double ev1 = g_EoC_Data->GetErrorY(i);
    double ev2 = g_EoC_MC->GetErrorY(i);
    
    if (ev1 > 5) {
      g_EoC_Data->SetPointError(i, 0, 0);
      ev1 = 0;
    }
    if (ev2 > 5) {
      g_EoC_MC->SetPointError(i, 0, 0);
      ev2 = 0;
    }

    double erat = rat*sqrt(  (ev1*ev1)/(v1*v1) + (ev2*ev2)/(v2*v2));
    gC_ratio->SetPoint(i, xx1, rat);
    gC_ratio->SetPointError(i, 0, erat);
  }


  // --- UNCORRECTED
  TCanvas* c_g_fit = new TCanvas("g_fit", "g_fit",100,100,1000,600);
  c_g_fit->Divide(1,2);
  c_g_fit->cd(1);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,ymin,etaMax,ymax);
  hPad->GetXaxis()->SetTitle("#eta_{SC}");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.05);

  g_EoP_MC -> Draw("PL");
  g_EoP_Data -> Draw("PL");

  TLegend *tl = new TLegend(0.60,0.15,0.89,0.35);
  tl -> SetFillColor(0);
  tl -> AddEntry(g_EoP_MC,"MC - default","PL");
  tl -> AddEntry(g_EoP_Data,"DATA - default","PL");
  tl -> Draw();

  // ratio
  c_g_fit->cd(2);
  gPad->SetGrid();
  
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,ymin,etaMax,ymax);
  hPad->GetXaxis()->SetTitle("#eta_{SC}");
  hPad->GetYaxis()->SetTitle("data/MC"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.05);

  gP_ratio -> SetMarkerStyle(20);
  gP_ratio -> SetMarkerSize(0.7);
  gP_ratio -> SetMarkerColor(kBlue); 
  gP_ratio -> Draw("PL");

  //****************************************************************

  // CORRECTED
  TCanvas* c_g_fitC = new TCanvas("g_fitC", "g_fitC",100,100,1000,600);
  c_g_fitC->Divide(1,2);
  c_g_fitC->cd(1);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,ymin,etaMax,ymax);
  hPad->GetXaxis()->SetTitle("i#eta_{SC}");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.05);

  g_EoC_MC -> Draw("PL");
  g_EoC_Data -> Draw("PL");

  TLegend *tlc = new TLegend(0.60,0.15,0.89,0.35);
  tlc -> SetFillColor(0);
  tlc -> AddEntry(g_EoC_MC,"MC - w/ regression","P");
  tlc -> AddEntry(g_EoC_Data,"DATA -w/ regression ","P");
  tlc -> Draw();

  // ratio
  c_g_fitC->cd(2);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,ymin,etaMax,ymax);
  hPad->GetXaxis()->SetTitle("i#eta_{SC}");
  hPad->GetYaxis()->SetTitle("data/MC"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.05);
  
  gC_ratio -> SetMarkerStyle(20);
  gC_ratio -> SetMarkerSize(0.7);
  gC_ratio -> SetMarkerColor(kBlue); 
  gC_ratio -> Draw("PL");



  //****************************************************************
  
  // MC vs MC corr
  TCanvas* c_g_fitMC = new TCanvas("g_fitMC", "g_fitMC",100,100,1000,600);
  c_g_fitMC->Divide(1,2);
  c_g_fitMC->cd(1);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,ymin,etaMax,ymax);
  hPad->GetXaxis()->SetTitle("i#eta_{SC}");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  hPad->GetXaxis()->SetTitleOffset(0.55);
  hPad->GetXaxis()->SetTitleSize(0.07);
  hPad->GetXaxis()->SetLabelSize(0.07); 
  hPad->GetYaxis()->SetTitleOffset(0.6);
  hPad->GetYaxis()->SetTitleSize(0.07);
  hPad->GetYaxis()->SetLabelSize(0.07);

  g_EoP_MC -> Draw("PL");
  g_EoC_MC -> Draw("PL");

  TLegend *tlmc = new TLegend(0.65,0.18,0.89,0.35);
  tlmc -> SetFillColor(0);
  tlmc -> AddEntry(g_EoP_MC,"MC - default","P");
  tlmc -> AddEntry(g_EoC_MC,"MC -w/ regression ","P");
  tlmc -> Draw();

  // ratio
  c_g_fitMC->cd(2);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,yrmin,etaMax,yrmax);
  hPad->GetXaxis()->SetTitle("i#eta_{SC}");
  hPad->GetYaxis()->SetTitle("regression/default"); 
  hPad->GetXaxis()->SetTitleOffset(0.55);
  hPad->GetXaxis()->SetTitleSize(0.07);
  hPad->GetXaxis()->SetLabelSize(0.07); 
  hPad->GetYaxis()->SetTitleOffset(0.6);
  hPad->GetYaxis()->SetTitleSize(0.07);
  hPad->GetYaxis()->SetLabelSize(0.07);
  g_Ratio_MC -> Draw("PL");

  //****************************************************************

  // DATA vs DATA corr
  TCanvas* c_g_fitDATA = new TCanvas("g_fitDATA", "g_fitDATA",100,100,1000,600);
  c_g_fitDATA->Divide(1,2);
  c_g_fitDATA->cd(1);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,ymin,etaMax,ymax);
  hPad->GetXaxis()->SetTitle("i#eta_{SC}");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  hPad->GetXaxis()->SetTitleOffset(0.55);
  hPad->GetXaxis()->SetTitleSize(0.07);
  hPad->GetXaxis()->SetLabelSize(0.07); 
  hPad->GetYaxis()->SetTitleOffset(0.6);
  hPad->GetYaxis()->SetTitleSize(0.07);
  hPad->GetYaxis()->SetLabelSize(0.07);


  g_EoP_Data -> Draw("PL");
  g_EoC_Data -> Draw("PL");

  TLegend *tlda = new TLegend(0.65,0.18,0.89,0.35);
  tlda -> SetFillColor(0);
  tlda -> AddEntry(g_EoP_Data,"DATA - default","P");
  tlda -> AddEntry(g_EoC_Data,"DATA -w/ regression ","P");
  tlda -> Draw();

  // ratio
  c_g_fitDATA->cd(2);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,yrmin,etaMax,yrmax);
  hPad->GetXaxis()->SetTitle("i#eta_{SC}");
  hPad->GetYaxis()->SetTitle("regression/default"); 
  hPad->GetXaxis()->SetTitleOffset(0.55);
  hPad->GetXaxis()->SetTitleSize(0.07);
  hPad->GetXaxis()->SetLabelSize(0.07); 
  hPad->GetYaxis()->SetTitleOffset(0.6);
  hPad->GetYaxis()->SetTitleSize(0.07);
  hPad->GetYaxis()->SetLabelSize(0.07);
  g_Ratio_Data -> Draw("PL");



  // data/MC ratios
  TCanvas* cratio = new TCanvas("cratio", "cratio",100,100,1000,400);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,ymin, etaMax,ymax);
  hPad->GetXaxis()->SetTitle("i#eta_{SC}");
  hPad->GetYaxis()->SetTitle("data/MC"); 
  hPad->GetXaxis()->SetTitleOffset(0.55);
  hPad->GetXaxis()->SetTitleSize(0.06);
  hPad->GetXaxis()->SetLabelSize(0.06); 
  hPad->GetYaxis()->SetTitleOffset(0.6);
  hPad->GetYaxis()->SetTitleSize(0.07);
  hPad->GetYaxis()->SetLabelSize(0.06);

  gP_ratio ->SetMarkerSize(0.7);
  gP_ratio ->SetMarkerColor(kAzure+2);
  gP_ratio ->SetLineColor(kAzure+2);
  gC_ratio ->SetMarkerSize(0.7);
  gC_ratio ->SetMarkerColor(kBlue+2);
  gC_ratio ->SetLineColor(kBlue+2);

  gP_ratio -> Draw("PL");
  gC_ratio -> Draw("PL");

  TLegend *tlda = new TLegend(0.65,0.18,0.89,0.35);
  tlda -> SetFillColor(0);
  tlda -> AddEntry(gP_ratio,"default","P");
  tlda -> AddEntry(gC_ratio,"w/ regression ","P");
  tlda -> Draw();



  FILE *fcorr = fopen("scalFactorsTommaso.txt","r");
  float scaleFactor;
  float ieta;
  TGraphErrors *gg = new TGraphErrors();
  int nn = 0;
  while(fscanf(fcorr,"%f %f",&ieta,&scaleFactor) !=EOF ) {
    gg->SetPoint(nn,ieta, 1./scaleFactor);
    nn++;
  }

  gg->SetMarkerStyle(20);
  gg->SetMarkerSize(0.7);
  gg->Draw("psame");

  TFile *_file0 = TFile::Open("outputFiles/RegionalScale_Zee.root");
  g_Rat->SetMarkerColor(kOrange);
  g_Rat->SetMarkerStyle(20);
  g_Rat->SetMarkerSize(0.7);
  g_Rat->Draw("psame");


//   // residual spread data wrt MC
//   TH1F* h_spreadRatio = new TH1F("h_spreadRatio", "spreadRatio", 200, 0.95, 1.05);
//   for (int i = 1; i < g_EoC_Data -> GetN()+1; i++){
  
//     double xda, yda;
//     double  xmc, ymc;
//     gC_ratio -> GetPoint(i,xda,yda);
//     if (fabs(xda)<85) h_spreadRatio->Fill(yda);

//   }

//   TCanvas *cspread = new TCanvas("cspread","cspread",500,500);
//   h_spreadRatio->Draw();
// }
}
