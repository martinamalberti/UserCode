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
  
  
  const int Ntempl = 4;

  TFile *f = TFile::Open("GraphsLocalEta_singleBCcorrectionPol2_Zee_nBCgt1.root");
  
  TGraphErrors* g_EoP_MC[Ntempl];
  TGraphErrors* g_EoC_MC[Ntempl];
   
  TGraphErrors* g_EoP_Data[Ntempl];
  TGraphErrors* g_EoC_Data[Ntempl];

  TGraphErrors* g_ratio_MC[Ntempl];
  TGraphErrors* g_ratio_Data[Ntempl];

  TGraphErrors* g_ratio_uncorr[Ntempl];
  TGraphErrors* g_ratio_corr[Ntempl];

  for (int mod=0; mod<Ntempl; mod++){
    char histoName[100];
    
    sprintf(histoName, "gEoP_MC_mod%d", mod+1);
    g_EoP_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoC_MC_mod%d", mod+1);
    g_EoC_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoP_Data_mod%d", mod+1);
    g_EoP_Data[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoC_Data_mod%d", mod+1);
    g_EoC_Data[mod]   = (TGraphErrors*)f->Get(histoName); 

    sprintf(histoName, "gRatio_MC_mod%d", mod+1);
    g_ratio_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
 
    sprintf(histoName, "gRatio_Data_mod%d", mod+1);
    g_ratio_Data[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gRatio_uncorr_mod%d", mod+1);
    g_ratio_uncorr[mod]  = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gRatio_corr_mod%d", mod+1);
    g_ratio_corr[mod]  = (TGraphErrors*)f->Get(histoName);  


    // set colors for plotting
    //-- mc uncorrected
    g_EoP_MC[mod] -> SetMarkerStyle(21);
    g_EoP_MC[mod] -> SetMarkerSize(1.);
    g_EoP_MC[mod] -> SetMarkerColor(kRed); 
    g_EoP_MC[mod] -> SetLineColor(kRed); 
    
    //-- mc corrected
    g_EoC_MC[mod] -> SetMarkerStyle(21);
    g_EoC_MC[mod] -> SetMarkerSize(1.);
    g_EoC_MC[mod] -> SetMarkerColor(kRed+2); 
    g_EoC_MC[mod] -> SetLineColor(kRed+2); 
  
    //-- data uncorrected
    g_EoP_Data[mod] -> SetMarkerStyle(20);
    g_EoP_Data[mod] -> SetMarkerSize(1.);
    g_EoP_Data[mod] -> SetMarkerColor(kGreen); 
    g_EoP_Data[mod] -> SetLineColor(kGreen); 

    //-- data corrected
    g_EoC_Data[mod] -> SetMarkerStyle(20);
    g_EoC_Data[mod] -> SetMarkerSize(1.);
    g_EoC_Data[mod] -> SetMarkerColor(kGreen+2); 
    g_EoC_Data[mod] -> SetLineColor(kGreen+2); 
 

    // ratio
    g_ratio_MC[mod]->SetLineColor(kMagenta+1);
    g_ratio_MC[mod]->SetMarkerColor(kMagenta+1);
    g_ratio_MC[mod]->SetMarkerStyle(21);
    g_ratio_MC[mod]->SetMarkerSize(1);
    
    g_ratio_Data[mod]->SetLineColor(kMagenta+3);
    g_ratio_Data[mod]->SetMarkerColor(kMagenta+3);
    g_ratio_Data[mod]->SetMarkerStyle(20);
    g_ratio_Data[mod]->SetMarkerSize(1);

    g_ratio_uncorr[mod]->SetLineColor(kAzure+2);
    g_ratio_uncorr[mod]->SetMarkerColor(kAzure+2);
    g_ratio_uncorr[mod]->SetMarkerStyle(20);
    g_ratio_uncorr[mod]->SetMarkerSize(1);
    
    g_ratio_corr[mod]->SetLineColor(kBlue+2);
    g_ratio_corr[mod]->SetMarkerColor(kBlue+2);
    g_ratio_corr[mod]->SetMarkerStyle(20);
    g_ratio_corr[mod]->SetMarkerSize(1);

  }
  

  TCanvas* c_g_fit[Ntempl];
  TLegend* tl[Ntempl];
  TLegend* tlr[Ntempl];
  TLegend* tlrr[Ntempl];

  for(int mod=0;mod<4;mod++) {
    char padName[100];
    sprintf(padName, "g_fit_mod%d", mod+1);
    c_g_fit[mod] = new TCanvas(padName,padName,600,300);
    
    float tYoffset = 0.8; 
    float labSize = 0.05;
    gPad->SetGrid();
    TH1F *hPad = (TH1F*)gPad->DrawFrame(-0.55,0.97,0.55,1.03);
    hPad->GetXaxis()->SetTitle("#eta_{SC} (deg)");
    hPad->GetYaxis()->SetTitle("Relative E/p scale");
    hPad->GetXaxis()->SetLabelSize(labSize);
    hPad->GetXaxis()->SetTitleSize(labSize);
    hPad->GetYaxis()->SetLabelSize(labSize);
    hPad->GetYaxis()->SetTitleSize(labSize);
    hPad->GetXaxis()->SetTitleOffset(tYoffset);
    hPad->GetYaxis()->SetTitleOffset(tYoffset);
    
    // ----- Drawing--------------------
    //    g_EoP_MC[mod]   -> Draw("PL");
    ///g_EoC_MC[mod]   -> Draw("PL");
    g_EoP_Data[mod] -> Draw("PL");
    g_EoC_Data[mod] -> Draw("PL");

    //----- legend
    tl[mod] = new TLegend(0.50,0.12,0.89,0.35);
    tl[mod] -> SetFillColor(0);
    
    //tl[mod] -> AddEntry(g_EoP_MC[mod],"MC - default","PL");
    //tl[mod] -> AddEntry(g_EoC_MC[mod],"MC - w/ regression","PL");
    tl[mod] -> AddEntry(g_EoP_Data[mod],"DATA - default","PL");
    tl[mod] -> AddEntry(g_EoC_Data[mod],"DATA - w/ BC correction","PL");
    tl[mod] -> Draw();

    
  }





}
