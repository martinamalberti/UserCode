
{
  gROOT->SetStyle("Plain");
  // gStyle->SetPadTickX(1);
  // gStyle->SetPadTickY(1);
  // gStyle->SetOptTitle(0); 
  // gStyle->SetOptStat(1110);
  // gStyle->SetOptFit(1);
  
  // gStyle->SetTextFont(42);
  // gStyle->SetTextSize(0.03);
  // gStyle->SetTitleFont(42,"xyz");
  // gStyle->SetTitleSize(0.03);
  // gStyle->SetLabelFont(42,"xyz");
  // gStyle->SetLabelSize(0.03);
  // gStyle->SetStatFont(42);
  // gStyle->SetStatFontSize(0.03);
  // gROOT->ForceStyle();
  

  TFile *f[3];
  f[0] = TFile::Open("outputFiles/mcZpeakVsNvtx_nChange0_PU.root");
  f[1] = TFile::Open("outputFiles/mcZpeakVsNvtx_nChange1_PU.root");
  f[2] = TFile::Open("outputFiles/mcZpeakVsNvtx_nChange2_PU.root");
  // f[0] = TFile::Open("outputFiles/dataZpeakVsNvtx_allR9_PU.root");
  // f[1] = TFile::Open("outputFiles/dataZpeakVsNvtx_highR9_PU.root");
  // f[2] = TFile::Open("outputFiles/dataZpeakVsNvtx_lowR9_PU.root");

  TGraphErrors *gEBsigma[3];
  TGraphErrors *gEBsigmac[3];
  TGraphErrors *gEBsigmar[3];
  TGraphErrors *gEBdeltam[3];
  TGraphErrors *gEBdeltamc[3];
  TGraphErrors *gEBdeltamr[3];
  RooPlot *plotEBEB_all[3];
  RooPlot *plotEBEBc_all[3];
  RooPlot *plotEBEBr_all[3];

  TH1F *h_mZ_EBEB_all[3];
  TH1F *h_mZc_EBEB_all[3];
  TH1F *h_mZr_EBEB_all[3];
 

  for (int i = 0; i < 3; i++){

    gEBsigma[i] =(TGraphErrors*)f[i]->Get("gEBsigma");
    gEBsigmac[i]=(TGraphErrors*)f[i]->Get("gEBsigmac");
    gEBsigmar[i]=(TGraphErrors*)f[i]->Get("gEBsigmar");


    gEBdeltam[i] =(TGraphErrors*)f[i]->Get("gEBdeltam");
    gEBdeltamc[i]=(TGraphErrors*)f[i]->Get("gEBdeltamc");
    gEBdeltamr[i]=(TGraphErrors*)f[i]->Get("gEBdeltamr");

    plotEBEB_all[i]  = (RooPlot*)f[i]->Get("plotEBEB_all");
    plotEBEBc_all[i] = (RooPlot*)f[i]->Get("plotEBEBc_all");
    plotEBEBr_all[i] = (RooPlot*)f[i]->Get("plotEBEBr_all");
    
    h_mZ_EBEB_all[i]  = (TH1F*)f[i]->Get("h_mZ_EBEB_all");
    h_mZc_EBEB_all[i] = (TH1F*)f[i]->Get("h_mZc_EBEB_all");
    h_mZr_EBEB_all[i] = (TH1F*)f[i]->Get("h_mZr_EBEB_all");


    gEBsigmac[i]->SetMarkerColor(kGreen+3);
    gEBsigmac[i]->SetLineColor(kGreen+3);
    gEBdeltamc[i]->SetMarkerColor(kGreen+3);
    gEBdeltamc[i]->SetLineColor(kGreen+3);
  }
  
  // -- fraction of high r9
  cout << "@@@ Default: N(both highR9)/N(tot) = " << h_mZ_EBEB_all[1]->Integral()/h_mZ_EBEB_all[0]->Integral() << endl;
  cout << "@@@ Dyn. Clust.: N(both highR9)/N(tot) = " << h_mZc_EBEB_all[1]->Integral()/h_mZc_EBEB_all[0]->Integral() << endl;
  cout << "@@@ Default: N(both lowR9)/N(tot) = " << h_mZ_EBEB_all[2]->Integral()/h_mZ_EBEB_all[0]->Integral() << endl;
  cout << "@@@ Dyn. Clust.: N(both lowR9)/N(tot) = " << h_mZc_EBEB_all[2]->Integral()/h_mZc_EBEB_all[0]->Integral() << endl;

  
  float MZ = 91.188;
  for (int ii= 0; ii < 3 ; ii++){
    for (int i = 0; i < gEBsigma[ii]->GetN(); i++){
      double s , es, nvtx, dm;
      gEBsigma[ii]->GetPoint(i,nvtx,s);
      gEBdeltam[ii]->GetPoint(i,nvtx,dm);
      es = gEBsigma[ii]->GetErrorY(i);
      gEBsigma[ii]->SetPoint(i, nvtx, s/(dm+MZ));
      gEBsigma[ii]->SetPointError(i, 0.5, es/(dm+MZ));

      gEBsigmac[ii]->GetPoint(i,nvtx,s);
      gEBdeltamc[ii]->GetPoint(i,nvtx,dm);
      es = gEBsigmac[ii]->GetErrorY(i);
      gEBsigmac[ii]->SetPoint(i, nvtx, s/(dm+MZ));
      gEBsigmac[ii]->SetPointError(i, 0.5, es/(dm+MZ));

      gEBsigmar[ii]->GetPoint(i,nvtx,s);
      gEBdeltamr[ii]->GetPoint(i,nvtx,dm);
      es = gEBsigmar[ii]->GetErrorY(i);
      gEBsigmar[ii]->SetPoint(i, nvtx, s/(dm+MZ));
      gEBsigmar[ii]->SetPointError(i, 0.5, es/(dm+MZ));
    }
  }

  TLegend * leg = new TLegend(0.15,0.75, 0.35, 0.89);
  leg->SetFillColor(0);
  leg->AddEntry(gEBsigma[0],"default clustering", "LP");
  leg->AddEntry(gEBsigmac[0],"dynamic seed clustering", "LP");


  //-- sigmaCB
  TCanvas *csigma[3];
  char cname[100];
  for (int i = 0; i<3; i++){
    if (i==0) sprintf(cname, "sigmaCB");
    if (i==1) sprintf(cname, "sigmaCB_n1");
    if (i==2) sprintf(cname, "sigmaCB_n2");
    csigma[i] = new TCanvas(cname,cname,500,500);
    csigma[i] ->SetLeftMargin(0.15);
    csigma[i] ->SetTitle(cname);
    csigma[i] ->SetGridx();
    csigma[i] ->SetGridy();
    gEBsigma[i]->GetXaxis()->SetTitle("N_{VTX}");
    //gEBsigma[i]->GetYaxis()->SetTitle("#sigma_{CB} (GeV)");
    gEBsigma[i]->GetYaxis()->SetTitleOffset(1.8);
    gEBsigma[i]->GetYaxis()->SetTitle("#sigma_{CB} (%)");
    gEBsigma[i]->GetXaxis()->SetRangeUser(1.5,20);
    gEBsigma[i]->GetYaxis()->SetRangeUser(0,0.04);
    gEBsigma[i]->Draw("ap");
    gEBsigmac[i]->Draw("psame");
    // gEBsigmar[i]->Draw("psame");
    leg->Draw("same");
  }


  //-- deltaM
  TCanvas *cdeltam[3];
  for (int i = 0; i<3; i++){
    if (i==0) sprintf(cname, "deltaM");
    if (i==1) sprintf(cname, "deltaM_n1");
    if (i==2) sprintf(cname, "deltaM_n2");
    cdeltam[i] = new TCanvas(cname,cname,500,500);
    cdeltam[i] ->SetTitle(cname);
    cdeltam[i] ->SetGridx();
    cdeltam[i] ->SetGridy();
    gEBdeltam[i]->GetXaxis()->SetTitle("N_{VTX}");
    gEBdeltam[i]->GetYaxis()->SetTitle("#Delta M (GeV)");
    gEBdeltam[i]->GetXaxis()->SetRangeUser(1.5,20);
    gEBdeltam[i]->GetYaxis()->SetRangeUser(-4,4);
    gEBdeltam[i]->Draw("ap");
    gEBdeltamc[i]->Draw("psame");
    //gEBdeltamr[i]->Draw("psame");
    leg->Draw("same");
  }


  //-- deltaM
  TCanvas *cZpeak[3];
  for (int i = 0; i<3; i++){
    if (i==0) sprintf(cname, "cZpeak");
    if (i==1) sprintf(cname, "cZpeak_n1");
    if (i==2) sprintf(cname, "cZpeak_n2");
    cZpeak[i] = new TCanvas(cname,cname,500,500);
    plotEBEB_all[i]->Draw("");
  }

  TCanvas *cZpeakPUcleaned[3];
  for (int i = 0; i<3; i++){
    if (i==0) sprintf(cname, "cZpeakPUcleaned");
    if (i==1) sprintf(cname, "cZpeakPUcleaned_n1");
    if (i==2) sprintf(cname, "cZpeakPUcleaned_n2");
    cZpeakPUcleaned[i] = new TCanvas(cname,cname,400,400);
    plotEBEBc_all[i]->Draw("");
  }


}
