

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
  
  TFile *fw = TFile::Open("outputFiles/GraphsRegionalScaleEta_highR9_Pcalib_templatesFromMC_42XReReco_FT_R_42_V21B.root");
  TFile *fz = TFile::Open("outputFiles/GraphsRegionalScaleEta_highR9_Pcalib_templatesFromMC_42XReReco_FT_R_42_V24.root");
  
  char gname[100];
  sprintf(gname,"gEoP_Data");

  TGraphErrors* gw  = (TGraphErrors*)fw->Get(gname);
  gw -> SetMarkerStyle(20);
  gw -> SetMarkerSize(0.7);
  gw -> SetMarkerColor(kRed); 

  TGraphErrors* gz   = (TGraphErrors*)fz->Get(gname);
  gz -> SetMarkerStyle(20);
  gz -> SetMarkerSize(0.7);
  gz -> SetMarkerColor(kBlue); 
  
  // -- build ratio
  TGraphErrors* gP_ratio = new TGraphErrors();
  gP_ratio->SetMarkerColor(kAzure+2);
  gP_ratio->SetMarkerStyle(20);
  gP_ratio->SetMarkerSize(0.7);

  for(int i =0 ; i < gw ->GetN();i++){
    double v1, v2, v3, xx1, xx2, xx3;
    gw->GetPoint(i,xx1,v1);
    gz->GetPoint(i,xx2,v2);
    if( xx1 != xx2 ){ cout<<"Error 2 !!"<<endl; return;}
    double rat = v1/v2; 
    double ev1 = gw->GetErrorY(i);
    double ev2 = gz->GetErrorY(i);

    if (ev1 > 0.5) {
      gw->SetPointError(i, 0, 0);
      ev1 = 0;
    }
    if (ev2 > 0.5) {
      gz->SetPointError(i, 0, 0);
      ev2 = 0;
    }
    double erat = rat*sqrt(  (ev1*ev1)/(v1*v1) + (ev2*ev2)/(v2*v2));
    gP_ratio->SetPoint(i, xx1, rat);
    gP_ratio->SetPointError(i, 0, erat);
  }

  TCanvas* c = new TCanvas("c", "c",100,100,1000,600);
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGrid();gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,ymin, etaMax,ymax);
  hPad->GetXaxis()->SetTitle("i#eta_{SC}");
  hPad->GetYaxis()->SetTitle("relative E/p"); 
  hPad->GetXaxis()->SetTitleOffset(0.55);
  hPad->GetXaxis()->SetTitleSize(0.06);
  hPad->GetXaxis()->SetLabelSize(0.06); 
  hPad->GetYaxis()->SetTitleOffset(0.6);
  hPad->GetYaxis()->SetTitleSize(0.07);
  hPad->GetYaxis()->SetLabelSize(0.06);

  gw -> Draw("PL");
  gz -> Draw("PL");

  TLegend *tlda = new TLegend(0.65,0.18,0.89,0.35);
  tlda -> SetFillColor(0);
  tlda -> AddEntry(gw,"FT_R_42_V21B","P");
  tlda -> AddEntry(gz,"FT_R_42_V24","P");
  tlda -> Draw();
  
  c->cd(2);
  gPad->SetGrid();gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(etaMin,ymin, etaMax,ymax);
  hPad->GetXaxis()->SetTitle("i#eta_{SC}");
  hPad->GetYaxis()->SetTitle("ratio"); 
  hPad->GetXaxis()->SetTitleOffset(0.55);
  hPad->GetXaxis()->SetTitleSize(0.07);
  hPad->GetXaxis()->SetLabelSize(0.07); 
  hPad->GetYaxis()->SetTitleOffset(0.6);
  hPad->GetYaxis()->SetTitleSize(0.07);
  hPad->GetYaxis()->SetLabelSize(0.07);
  gP_ratio-> Draw("PL");

}
