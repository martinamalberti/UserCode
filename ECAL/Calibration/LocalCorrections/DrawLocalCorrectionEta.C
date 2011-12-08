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

  TFile *f = TFile::Open("GraphsLocalEta.root");
  
  
  TGraphErrors* g_EoP_MC[Ntempl];
  TGraphErrors* g_EoC_MC[Ntempl];
  TGraphErrors* g_ratio_MC[Ntempl];
  
  TGraphErrors* g_EoP_Data[Ntempl];
  TGraphErrors* g_EoC_Data[Ntempl];
  TGraphErrors* g_ratio_Data[Ntempl];

  TGraphErrors* g_ratio_uncorr[Ntempl];
  TGraphErrors* g_ratio_corr[Ntempl];
 
  for (int mod=0; mod<4; mod++){
    char histoName[100];
    
    sprintf(histoName, "gEoP_MC_mod%d", mod+1);
    g_EoP_MC[mod]   = (TGraphErrors*)f->Get(histoName); 
    
    sprintf(histoName, "gEoC_MC_mod%d", mod+1);
    g_EoC_MC[mod]   = (TGraphErrors*)f->Get(histoName); 

    sprintf(histoName, "gEoP_Data_mod%d", mod+1);
    g_EoP_Data[mod]   = (TGraphErrors*)f->Get(histoName); 

    sprintf(histoName, "gEoC_Data_mod%d", mod+1);
    g_EoC_Data[mod]   = (TGraphErrors*)f->Get(histoName); 

    sprintf(histoName, "gRatio_uncorr_mod%d", mod+1);
    g_ratio_uncorr[mod]  = (TGraphErrors*)f->Get(histoName); 

    sprintf(histoName, "gRatio_corr_mod%d", mod+1);
    g_ratio_corr[mod]  = (TGraphErrors*)f->Get(histoName);  
  }
  

  TCanvas* c_g_fit[Ntempl];
  TLegend* tl[Ntempl];
  TLegend* tlr[Ntempl];

  for(int mod=0;mod<4;mod++) {
    char padName[100];
    sprintf(padName, "g_fit_mod%d", mod+1);
    c_g_fit[mod] = new TCanvas(padName,padName,100,100,700,600);
    c_g_fit[mod]->Divide(1,2);

    c_g_fit[mod]->cd(1);

    gPad->SetGrid();
    TH1F *hPad = (TH1F*)gPad->DrawFrame(-0.55,0.97,0.55,1.02);
    hPad->GetXaxis()->SetTitle("#eta_{SC} (deg)");
    hPad->GetYaxis()->SetTitle("Relative E/p scale");
    hPad->GetXaxis()->SetTitleOffset(0.8);
    hPad->GetYaxis()->SetTitleSize(0.05);
    hPad->GetYaxis()->SetLabelSize(0.05);

    sprintf(padName, "mod%d", mod+1);
    hPad -> SetTitle(padName);
    

    g_EoP_MC[mod] -> SetMarkerStyle(20);
    g_EoP_MC[mod] -> SetMarkerSize(1.);
    g_EoP_MC[mod] -> SetMarkerColor(kRed); 
    g_EoP_MC[mod] -> SetLineColor(kRed); 
    //    g_EoP_MC[mod] -> Draw("PL");
    g_EoC_MC[mod] -> SetMarkerStyle(20);
    g_EoC_MC[mod] -> SetMarkerSize(1.);
    g_EoC_MC[mod] -> SetMarkerColor(kRed+2); 
    g_EoC_MC[mod] -> SetLineColor(kRed+2); 
    g_EoC_MC[mod] -> Draw("PL");
  
    g_EoP_Data[mod] -> SetMarkerStyle(20);
    g_EoP_Data[mod] -> SetMarkerSize(1.);
    g_EoP_Data[mod] -> SetMarkerColor(kGreen); 
    g_EoP_Data[mod] -> SetLineColor(kGreen); 
    //g_EoP_Data[mod] -> Draw("PL");
    g_EoC_Data[mod] -> SetMarkerStyle(20);
    g_EoC_Data[mod] -> SetMarkerSize(1.);
    g_EoC_Data[mod] -> SetMarkerColor(kGreen+2); 
    g_EoC_Data[mod] -> SetLineColor(kGreen+2); 
    g_EoC_Data[mod] -> Draw("PL");
  

    tl[mod] = new TLegend(0.60,0.15,0.89,0.35);
    tl[mod] -> SetFillColor(0);
    tl[mod] -> AddEntry(g_EoC_MC[mod],"MC","PL");
    tl[mod] -> AddEntry(g_EoC_Data[mod],"DATA","PL");
//     tl[mod] -> AddEntry(g_EoP_MC[mod],"MC uncorrected","PL");
//     tl[mod] -> AddEntry(g_EoC_MC[mod],"MC corrected","PL");
//     tl[mod] -> AddEntry(g_EoP_Data[mod],"Data uncorrected","PL");
//     tl[mod] -> AddEntry(g_EoC_Data[mod],"Data corrected","PL");
    tl[mod] -> Draw();

    c_g_fit[mod]->cd(2);
    gPad->SetGrid();
    TH1F *hPad2 = (TH1F*)gPad->DrawFrame(-0.55,0.97,0.55,1.03);
    hPad2->GetXaxis()->SetTitle("#eta_{SC} (deg)");
    hPad2->GetYaxis()->SetTitle("data/MC ratio");
    hPad2->GetXaxis()->SetTitleOffset(0.8);
    hPad2->GetYaxis()->SetTitleSize(0.05);
    hPad2->GetYaxis()->SetLabelSize(0.05);
  
    g_ratio_uncorr[mod]->SetLineColor(kBlue);
    g_ratio_uncorr[mod]->SetMarkerColor(kBlue);
    g_ratio_uncorr[mod]->SetMarkerStyle(20);
    g_ratio_uncorr[mod]->SetMarkerSize(0.7);
    //    g_ratio_uncorr[mod]->Draw("PL");
    
    g_ratio_corr[mod]->SetLineColor(kBlue+3);
    g_ratio_corr[mod]->SetMarkerColor(kBlue+3);
    g_ratio_corr[mod]->SetMarkerStyle(20);
    g_ratio_corr[mod]->SetMarkerSize(0.7);
    g_ratio_corr[mod]->Draw("PL");
   
    tlr[mod] = new TLegend(0.60,0.15,0.89,0.35);
    tlr[mod] -> SetFillColor(0);
    tlr[mod] -> AddEntry(g_ratio_uncorr[mod],"uncorrected","PL");
    tlr[mod] -> AddEntry(g_ratio_corr[mod],"corrected","PL");
    //tlr[mod] -> Draw();

  }





}
