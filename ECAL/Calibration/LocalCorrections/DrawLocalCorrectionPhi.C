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

  //TFile *f = TFile::Open("GraphsLocalPhi_regression_highR9_etacharge_pos.root");
  TFile *f = TFile::Open("GraphsLocalPhi_regression_DYee.root");
  
  
  TGraphErrors* g_EoP_MC[Ntempl];
  TGraphErrors* g_EoC_MC[Ntempl];
  TGraphErrors* g_ratio_MC[Ntempl];
  
  TGraphErrors* g_EoP_Data[Ntempl];
  TGraphErrors* g_EoC_Data[Ntempl];
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
    c_g_fit[mod] = new TCanvas(padName,padName,100,100,700,700);
    
    //    c_g_fit[mod]->Divide(1,2);
    //c_g_fit[mod]->cd(1);

    TPad *cLower  = new TPad("pad_0","pad_0",0.00,0.00,1.00,0.25);
    TPad *cCenter = new TPad("pad_1","pad_1",0.00,0.25,1.00,0.50);
    TPad *cUpper  = new TPad("pad_2","pad_2",0.00,0.50,1.00,1.00);
    
    cLower->SetBottomMargin(0.21); 
    cLower->SetTopMargin(0.05);

    cCenter->SetBottomMargin(0.1); 
    cCenter->SetTopMargin(0.05);

    cUpper->SetBottomMargin(0.1); 
    cUpper->SetTopMargin(0.05);

    cLower->Draw();
    cCenter->Draw();
    cUpper->Draw();
    
    float FontSCF = cUpper->GetHNDC()/cLower->GetHNDC(); 
    float tYoffset = 0.9; 
    float labSizeY = 0.055;
    float labSizeX = 0.045;

    cUpper-> cd();
    gPad->SetGrid();
    TH1F *hPad = (TH1F*)gPad->DrawFrame(-0.55,0.985,0.55,1.01);
    hPad->GetXaxis()->SetTitle("#phi_{SC} (deg)");
    hPad->GetYaxis()->SetTitle("Relative E/p scale");
    hPad->GetXaxis()->SetLabelSize(labSizeX);
    hPad->GetXaxis()->SetTitleSize(labSizeX);
    hPad->GetYaxis()->SetLabelSize(labSizeY);
    hPad->GetYaxis()->SetTitleSize(labSizeY);
    hPad->GetXaxis()->SetTitleOffset(tYoffset);
    hPad->GetYaxis()->SetTitleOffset(tYoffset);

    // ----- Drawing--------------------
    g_EoP_MC[mod]   -> Draw("PL");
    g_EoC_MC[mod]   -> Draw("PL");
    g_EoP_Data[mod] -> Draw("PL");
    g_EoC_Data[mod] -> Draw("PL");

    //----- legend
    tl[mod] = new TLegend(0.60,0.12,0.89,0.35);
    tl[mod] -> SetFillColor(0);
    
    tl[mod] -> AddEntry(g_EoP_MC[mod],"MC - default","PL");
    tl[mod] -> AddEntry(g_EoC_MC[mod],"MC - w/ regression","PL");
    tl[mod] -> AddEntry(g_EoP_Data[mod],"DATA - default","PL");
    tl[mod] -> AddEntry(g_EoC_Data[mod],"DATA - w/ regression","PL");
    tl[mod] -> Draw();


    //--------- RATIO PLOTS ------------------------------------------

    //    c_g_fit[mod]->cd(2);
    
    cCenter-> cd();
    gPad->SetGrid();
    TH1F *hPad2 = (TH1F*)gPad->DrawFrame(-0.55,0.99,0.55,1.01);
    hPad2->GetXaxis()->SetTitle("#phi_{SC} (deg)");
    hPad2->GetYaxis()->SetTitle("regr./default ratio");
    hPad2->GetXaxis()->SetLabelSize(labSizeX*FontSCF);
    hPad2->GetXaxis()->SetTitleSize(labSizeX*FontSCF);
    hPad2->GetYaxis()->SetLabelSize(labSizeY*FontSCF);
    hPad2->GetYaxis()->SetTitleSize(labSizeY*FontSCF);
    hPad2->GetYaxis()->SetTitleOffset(tYoffset/FontSCF);
    hPad2->GetYaxis()->SetNdivisions(505);
    g_ratio_MC[mod]->Draw("PL");
    g_ratio_Data[mod]->Draw("PL");
    
    tlr[mod] = new TLegend(0.60,0.15,0.89,0.35);
    tlr[mod] -> SetFillColor(0);
    tlr[mod] -> AddEntry(g_ratio_MC[mod],"MC","PL");
    tlr[mod] -> AddEntry(g_ratio_Data[mod],"DATA","PL");
    tlr[mod] -> Draw();


    cLower-> cd();
    gPad->SetGrid();
    TH1F *hPad3 = (TH1F*)gPad->DrawFrame(-0.55,0.993,0.55,1.007);
    hPad3->GetXaxis()->SetTitle("#phi_{SC} (deg)");
    hPad3->GetYaxis()->SetTitle("data/MC ratio");
    hPad3->GetXaxis()->SetLabelSize(labSizeX*FontSCF);
    hPad3->GetXaxis()->SetTitleSize(labSizeX*FontSCF);
    hPad3->GetYaxis()->SetLabelSize(labSizeY*FontSCF);
    hPad3->GetYaxis()->SetTitleSize(labSizeY*FontSCF);
    hPad3->GetYaxis()->SetTitleOffset(tYoffset/FontSCF);
    hPad3->GetXaxis()->SetTitleOffset(0.9);
    hPad3->GetYaxis()->SetNdivisions(505);
    g_ratio_uncorr[mod]->Draw("PL");
    g_ratio_corr[mod]->Draw("PL");
    

    tlrr[mod] = new TLegend(0.60,0.25,0.89,0.45);
    tlrr[mod] -> SetFillColor(0);
    tlrr[mod] -> AddEntry(g_ratio_uncorr[mod],"default","PL");
    tlrr[mod] -> AddEntry(g_ratio_corr[mod],"w/ regression","PL");
    tlrr[mod] -> Draw();

    
  }





}
