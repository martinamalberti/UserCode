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
  
  TFile *f = TFile::Open("regionalScale.root");
  

  TGraphErrors* g_EoP_MC   = (TGraphErrors*)f->Get("gEoP_MC");
  TGraphErrors* g_EoC_MC   = (TGraphErrors*)f->Get("gEoC_MC");  
  TGraphErrors* g_EoP_Data = (TGraphErrors*)f->Get("gEoP_Data");
  TGraphErrors* g_EoC_Data = (TGraphErrors*)f->Get("gEoC_Data");

  TGraphErrors* g_ratio = new TGraphErrors();
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
    
    g_ratio->SetPoint(i, xx1, rat);
    g_ratio->SetPointError(i, 0, erat);

  }


  // --- MC
  TCanvas* c_g_fit_MC = new TCanvas("g_fit_MC", "g_fit_MC",100,100,1000,600);
  c_g_fit_MC->Divide(1,2);
  c_g_fit_MC->cd(1);
  
  gPad->SetGrid();
  
  TH1F *hPad = (TH1F*)gPad->DrawFrame(-1.45,0.95,1.45,1.03);
  hPad->GetXaxis()->SetTitle("#eta_{SC}");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.05);

  g_EoC_MC -> SetMarkerStyle(20);
  g_EoC_MC -> SetMarkerSize(0.7);
  g_EoC_MC -> SetMarkerColor(kRed+2); 
  g_EoC_MC -> Draw("PL");

  g_EoC_Data -> SetMarkerStyle(20);
  g_EoC_Data -> SetMarkerSize(0.7);
  g_EoC_Data -> SetMarkerColor(kGreen+2); 
  g_EoC_Data -> Draw("PL");

  TLegend *tl = new TLegend(0.60,0.15,0.89,0.35);
  tl -> SetFillColor(0);
  tl -> AddEntry(g_EoC_MC,"MC","P");
  tl -> AddEntry(g_EoC_Data,"DATA","P");
  tl -> Draw();

 
  // ratio
  c_g_fit_MC->cd(2);
  gPad->SetGrid();
  
  TH1F *hPad = (TH1F*)gPad->DrawFrame(-1.45,0.95,1.45,1.03);
  hPad->GetXaxis()->SetTitle("#eta_{SC}");
  hPad->GetYaxis()->SetTitle("data/MC"); 
  hPad->GetXaxis()->SetTitleOffset(0.9);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetLabelSize(0.05);

  g_ratio -> SetMarkerStyle(20);
  g_ratio -> SetMarkerSize(0.7);
  g_ratio -> SetMarkerColor(kBlue); 
  g_ratio -> Draw("PL");
  
}
