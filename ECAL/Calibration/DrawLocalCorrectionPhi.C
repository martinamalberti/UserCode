{
  // Elettroni in EB+ sono come positroni in EB -

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

  TFile *f[2];
  //f[0] = TFile::Open("GraphsLocalPhi.root");
  //f[1] = TFile::Open("GraphsLocalPhi.root");

  f[0] = TFile::Open("GraphsLocalPhi_etacharge_pos_highR9.root");
  f[1] = TFile::Open("GraphsLocalPhi_etacharge_neg_highR9.root");

  TGraphErrors* g_EoP_MC[Ntempl][2];
  TGraphErrors* g_EoC_MC[Ntempl][2];
  TGraphErrors* g_ratio_MC[Ntempl][2];
  
  TGraphErrors* g_EoP_Data[Ntempl][2];
  TGraphErrors* g_EoC_Data[Ntempl][2];
  TGraphErrors* g_ratio_Data[Ntempl][2];

  TGraphErrors* g_ratio_uncorr[Ntempl][2];
  TGraphErrors* g_ratio_corr[Ntempl][2];
 
  for (int ifile = 0; ifile < 2 ; ifile++){
    for (int mod=0; mod<4; mod++){
      char histoName[100];
      
      sprintf(histoName, "gEoP_MC_mod%d", mod+1);
      g_EoP_MC[mod][ifile]   = (TGraphErrors*)f[ifile]->Get(histoName); 
      
      sprintf(histoName, "gEoC_MC_mod%d", mod+1);
      g_EoC_MC[mod][ifile]   = (TGraphErrors*)f[ifile]->Get(histoName); 
      
      sprintf(histoName, "gEoP_Data_mod%d", mod+1);
      g_EoP_Data[mod][ifile]   = (TGraphErrors*)f[ifile]->Get(histoName); 
      
      sprintf(histoName, "gEoC_Data_mod%d", mod+1);
      g_EoC_Data[mod][ifile]   = (TGraphErrors*)f[ifile]->Get(histoName); 
      
      sprintf(histoName, "gRatio_uncorr_mod%d", mod+1);
      g_ratio_uncorr[mod][ifile]  = (TGraphErrors*)f[ifile]->Get(histoName); 
      
      sprintf(histoName, "gRatio_corr_mod%d", mod+1);
      g_ratio_corr[mod][ifile]  = (TGraphErrors*)f[ifile]->Get(histoName);  
    }
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
    TH1F *hPad = (TH1F*)gPad->DrawFrame(-0.55,0.985,0.55,1.01);
    hPad->GetXaxis()->SetTitle("#phi_{SC} (deg)");
    hPad->GetYaxis()->SetTitle("Relative E/p scale");
    hPad->GetXaxis()->SetTitleOffset(0.8);
    hPad->GetYaxis()->SetTitleSize(0.05);
    hPad->GetYaxis()->SetLabelSize(0.05);

    sprintf(padName, "mod%d", mod+1);
    hPad -> SetTitle(padName);
   
    g_EoP_MC[mod][0] -> SetMarkerStyle(20);
    g_EoP_MC[mod][0] -> SetMarkerSize(1.);
    g_EoP_MC[mod][0] -> SetMarkerColor(kRed); 
    g_EoP_MC[mod][0] -> SetLineColor(kRed); 
    //    g_EoP_MC[mod][0] -> Draw("PL");
    g_EoC_MC[mod][0] -> SetMarkerStyle(20);
    g_EoC_MC[mod][0] -> SetMarkerSize(1.);
    g_EoC_MC[mod][0] -> SetMarkerColor(kRed+2); 
    g_EoC_MC[mod][0] -> SetLineColor(kRed+2); 
    g_EoC_MC[mod][0] -> Draw("PL");
  
    g_EoP_Data[mod][0] -> SetMarkerStyle(20);
    g_EoP_Data[mod][0] -> SetMarkerSize(1.);
    g_EoP_Data[mod][0] -> SetMarkerColor(kGreen); 
    g_EoP_Data[mod][0] -> SetLineColor(kGreen); 
    //g_EoP_Data[mod][0] -> Draw("PL");
    g_EoC_Data[mod][0] -> SetMarkerStyle(20);
    g_EoC_Data[mod][0] -> SetMarkerSize(1.);
    g_EoC_Data[mod][0] -> SetMarkerColor(kGreen+2); 
    g_EoC_Data[mod][0] -> SetLineColor(kGreen+2); 
    g_EoC_Data[mod][0] -> Draw("PL");
  
    g_EoP_MC[mod][1] -> SetMarkerStyle(24);
    g_EoP_MC[mod][1] -> SetMarkerSize(1.);
    g_EoP_MC[mod][1] -> SetMarkerColor(kRed); 
    g_EoP_MC[mod][1] -> SetLineColor(kRed); 
    //    g_EoP_MC[mod][1] -> Draw("PL");
    g_EoC_MC[mod][1] -> SetMarkerStyle(24);
    g_EoC_MC[mod][1] -> SetMarkerSize(1.);
    g_EoC_MC[mod][1] -> SetMarkerColor(kRed+2); 
    g_EoC_MC[mod][1] -> SetLineColor(kRed+2); 
    g_EoC_MC[mod][1] -> Draw("PL");
  
    g_EoP_Data[mod][1] -> SetMarkerStyle(24);
    g_EoP_Data[mod][1] -> SetMarkerSize(1.);
    g_EoP_Data[mod][1] -> SetMarkerColor(kGreen); 
    g_EoP_Data[mod][1] -> SetLineColor(kGreen); 
    //g_EoP_Data[mod][1] -> Draw("PL");
    g_EoC_Data[mod][1] -> SetMarkerStyle(24);
    g_EoC_Data[mod][1] -> SetMarkerSize(1.);
    g_EoC_Data[mod][1] -> SetMarkerColor(kGreen+2); 
    g_EoC_Data[mod][1] -> SetLineColor(kGreen+2); 
    g_EoC_Data[mod][1] -> Draw("PL");

    tl[mod] = new TLegend(0.60,0.15,0.89,0.45);
    tl[mod] -> SetFillColor(0);
//     tl[mod] -> AddEntry(g_EoC_MC[mod][0],"MC","PL");
//     tl[mod] -> AddEntry(g_EoC_Data[mod][0],"DATA","PL");
   
    tl[mod] -> AddEntry(g_EoC_MC[mod][0],"MC #eta*charge > 0","PL");
    tl[mod] -> AddEntry(g_EoC_Data[mod][0],"DATA #eta*charge > 0","PL");
    tl[mod] -> AddEntry(g_EoC_MC[mod][1],"MC #eta*charge < 0","PL");
    tl[mod] -> AddEntry(g_EoC_Data[mod][1],"DATA #eta*charge < 0","PL");
//     tl[mod] -> AddEntry(g_EoP_MC[mod],"MC uncorrected","PL");
//     tl[mod] -> AddEntry(g_EoC_MC[mod],"MC corrected","PL");
//     tl[mod] -> AddEntry(g_EoP_Data[mod],"Data uncorrected","PL");
//     tl[mod] -> AddEntry(g_EoC_Data[mod],"Data corrected","PL");
    tl[mod] -> Draw();

    c_g_fit[mod]->cd(2);
    gPad->SetGrid();
    TH1F *hPad2 = (TH1F*)gPad->DrawFrame(-0.55,0.97,0.55,1.03);
    hPad2->GetXaxis()->SetTitle("#phi_{SC} (deg)");
    hPad2->GetYaxis()->SetTitle("data/MC ratio");
    hPad2->GetXaxis()->SetTitleOffset(0.8);
    hPad2->GetYaxis()->SetTitleSize(0.05);
    hPad2->GetYaxis()->SetLabelSize(0.05);
  
    g_ratio_uncorr[mod][0]->SetLineColor(kBlue);
    g_ratio_uncorr[mod][0]->SetMarkerColor(kBlue);
    g_ratio_uncorr[mod][0]->SetMarkerStyle(20);
    g_ratio_uncorr[mod][0]->SetMarkerSize(0.7);
    //g_ratio_uncorr[mod][0]->Draw("PL");
    
    g_ratio_corr[mod][0]->SetLineColor(kBlue+3);
    g_ratio_corr[mod][0]->SetMarkerColor(kBlue+3);
    g_ratio_corr[mod][0]->SetMarkerStyle(20);
    g_ratio_corr[mod][0]->SetMarkerSize(0.7);
    g_ratio_corr[mod][0]->Draw("PL");

    g_ratio_uncorr[mod][1]->SetLineColor(kBlue);
    g_ratio_uncorr[mod][1]->SetMarkerColor(kBlue);
    g_ratio_uncorr[mod][1]->SetMarkerStyle(24);
    g_ratio_uncorr[mod][1]->SetMarkerSize(0.7);
    //    g_ratio_uncorr[mod][1]->Draw("PL");
    
    g_ratio_corr[mod][1]->SetLineColor(kBlue+3);
    g_ratio_corr[mod][1]->SetMarkerColor(kBlue+3);
    g_ratio_corr[mod][1]->SetMarkerStyle(24);
    g_ratio_corr[mod][1]->SetMarkerSize(0.7);
    g_ratio_corr[mod][1]->Draw("PL");

    tlr[mod] = new TLegend(0.60,0.15,0.89,0.35);
    tlr[mod] -> SetFillColor(0);
    tlr[mod] -> AddEntry(g_ratio_uncorr[mod][0],"uncorrected","PL");
    tlr[mod] -> AddEntry(g_ratio_corr[mod][0],"corrected","PL");
    //tlr[mod] -> Draw();

  }





}
