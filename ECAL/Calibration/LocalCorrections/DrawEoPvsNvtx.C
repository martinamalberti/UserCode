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
  
  
  const int Ntempl = 1;

  TFile *f = TFile::Open("GraphsEoP_vs_nvtx_EB.root");
  
  TGraphErrors* g_EoP_MC[Ntempl];
  TGraphErrors* g_EoC_MC[Ntempl];
  TGraphErrors* g_EoR_MC[Ntempl];
    
  TGraphErrors* g_EoP_Data[Ntempl];
  TGraphErrors* g_EoC_Data[Ntempl];
  TGraphErrors* g_EoR_Data[Ntempl];
  
  for (int mod=0; mod<Ntempl; mod++){
    char histoName[100];
    
    sprintf(histoName, "gEoP_MC_mod%d", mod+1);
    g_EoP_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoC_MC_mod%d", mod+1);
    g_EoC_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoR_MC_mod%d", mod+1);
    g_EoR_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoP_Data_mod%d", mod+1);
    g_EoP_Data[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoC_Data_mod%d", mod+1);
    g_EoC_Data[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoR_Data_mod%d", mod+1);
    g_EoR_Data[mod]   = (TGraphErrors*)f->Get(histoName); 
    
  }

  TCanvas* c_g_fit[Ntempl];
  TLegend* tl[Ntempl];
  TLegend* tlr[Ntempl];

  for(int mod=0;mod<Ntempl;mod++) {
    char padName[100];
    sprintf(padName, "g_fit_mod%d", mod+1);
    c_g_fit[mod] = new TCanvas(padName,padName,100,100,700,600);

    c_g_fit[mod]->cd(1);

    gPad->SetGrid();
    TH1F *hPad = (TH1F*)gPad->DrawFrame(0,0.985,30,1.015);
    hPad->GetXaxis()->SetTitle("Number of vertices");
    hPad->GetYaxis()->SetTitle("Relative E/p scale");
    hPad->GetYaxis()->SetTitleOffset(1.3);
    hPad->GetYaxis()->SetTitleSize(0.04);
    hPad->GetYaxis()->SetLabelSize(0.04);
    hPad->GetXaxis()->SetTitleSize(0.04);
    hPad->GetXaxis()->SetLabelSize(0.04);

    sprintf(padName, "mod%d", mod+1);
    hPad -> SetTitle(padName);
   
    //----------------- first file 
    //-- mc uncorrected
    g_EoP_MC[mod] -> SetMarkerStyle(20);
    g_EoP_MC[mod] -> SetMarkerSize(1.);
    g_EoP_MC[mod] -> SetMarkerColor(kRed); 
    g_EoP_MC[mod] -> SetLineColor(kRed); 
    
    //-- mc corrected
    g_EoC_MC[mod] -> SetMarkerStyle(20);
    g_EoC_MC[mod] -> SetMarkerSize(1.);
    g_EoC_MC[mod] -> SetMarkerColor(kGreen+1); 
    g_EoC_MC[mod] -> SetLineColor(kGreen+1); 

    //-- mc corrected - regression
    g_EoR_MC[mod] -> SetMarkerStyle(20);
    g_EoR_MC[mod] -> SetMarkerSize(1.);
    g_EoR_MC[mod] -> SetMarkerColor(kBlue+1); 
    g_EoR_MC[mod] -> SetLineColor(kBlue+1); 
      
    // data - uncorrected
    g_EoP_Data[mod] -> SetMarkerStyle(24);
    g_EoP_Data[mod] -> SetMarkerSize(1.);
    g_EoP_Data[mod] -> SetMarkerColor(kRed); 
    g_EoP_Data[mod] -> SetLineColor(kRed); 

    // data - corrected
    g_EoC_Data[mod] -> SetMarkerStyle(24);
    g_EoC_Data[mod] -> SetMarkerSize(1.);
    g_EoC_Data[mod] -> SetMarkerColor(kGreen+1); 
    g_EoC_Data[mod] -> SetLineColor(kGreen+1); 

    // data - corrected regression
    g_EoR_Data[mod] -> SetMarkerStyle(24);
    g_EoR_Data[mod] -> SetMarkerSize(1.);
    g_EoR_Data[mod] -> SetMarkerColor(kBlue+1); 
    g_EoR_Data[mod] -> SetLineColor(kBlue+1); 

  

    // ----- Drawing--------------------
    g_EoP_MC[mod] -> Draw("P");
    g_EoC_MC[mod] -> Draw("P");
    g_EoR_MC[mod] -> Draw("P");
    g_EoP_Data[mod] -> Draw("P");
    g_EoC_Data[mod] -> Draw("P");
    g_EoR_Data[mod] -> Draw("P");

    //----- legend
    tl[mod] = new TLegend(0.60,0.15,0.89,0.45);
    tl[mod] -> SetFillColor(0);

    tl[mod] -> AddEntry(g_EoP_MC[mod],"MC uncorrected","PL");
    tl[mod] -> AddEntry(g_EoC_MC[mod],"MC - dyn. seed clust.","PL");
    tl[mod] -> AddEntry(g_EoR_MC[mod],"MC - regression","PL");
    tl[mod] -> AddEntry(g_EoP_Data[mod],"DATA uncorrected","PL");
    tl[mod] -> AddEntry(g_EoC_Data[mod],"DATA - dyn. seed clust.","PL");
    tl[mod] -> AddEntry(g_EoR_Data[mod],"DATA - regression","PL");
    tl[mod] -> Draw();

  }





}
