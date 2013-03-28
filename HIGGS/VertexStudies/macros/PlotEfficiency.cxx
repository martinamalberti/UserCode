TH1F* MakeRatio(TH1F* h1, TH1F*h2, string hname){

  TH1F *hratio = (TH1F*)h1->Clone(hname.c_str());
  hratio->Divide(h1,h2,1./h1->GetSumOfWeights(),1./h2->GetSumOfWeights());
  return(hratio);
}



void PlotEfficiency(string outdir, string lumi, bool plotVariables=false) 
{

  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();
  gStyle->SetErrorX(0.5);

  int saveScaleFactors    = 0;
  bool useVariableBinning = false;

  // 0 : data
  // 1 : mc

  TFile *f[2];

  //2012ABCD
  f[0] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DoubleMu_Run2012ABCD/testEfficiency.root");
  f[1] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012ABCD_BSrw/testEfficiency.root");
  //f[0] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DoubleMu_Run2012D/testEfficiency.root");
  //f[1] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012D/testEfficiency.root");
    
  std::string text = "#splitline{CMS preliminary}{#sqrt{s} = 8 TeV L = "+lumi+" fb^{-1}}";

  TLatex *latex = new TLatex(0.55,0.85,text.c_str());
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.04);
  
  TLatex *latex2 = new TLatex(0.15,0.85,"Z#rightarrow#mu#mu");
  latex2->SetNDC();
  latex2->SetTextFont(42);
  latex2->SetTextSize(0.06);


  // legend
  string legtitle1    = "DATA Z#rightarrow#mu#mu (RANK)";
  string legtitle1bdt = "DATA Z#rightarrow#mu#mu ";
  string legtitle2    = "MC Z#rightarrow#mu#mu (RANK)";
  string legtitle2bdt = "MC Z#rightarrow#mu#mu ";

  // rebinning
  int nRePt = 2;

  // if variable size bins
  double xbins[10] = {0.,10.,20.,30.,40.,50.,70.,110.,250.,400.};

  // histograms
  TH1F* NvtxAll[2];
  TH1F* NvtxBaseline[2];
  TH1F* NvtxBDT[2];

  TH1F* PtAll[2];
  TH1F* PtBaseline[2];
  TH1F* PtBDT[2];

  TH1* PtAllRebinned[2];
  TH1* PtBaselineRebinned[2];
  TH1* PtBDTRebinned[2];

  TH1F *effVsNvtx_BDT[2];
  TH1F *effVsNvtx_Baseline[2];

  TH1F *effVsPt_BDT[2];
  TH1F *effVsPt_Baseline[2];

  TH1F *effVsEta_BDT[2];
  TH1F *effVsEta_Baseline[2];

  TH1F *BDToutput[2];
  TH1F *BDToutput_sig[2];
  TH1F *BDToutput_bkg[2];

  TH1F *evtBDToutput[2];
  TH1F *evtBDToutput_sig[2];
  TH1F *evtBDToutput_bkg[2];

  TH1F *ptbal_sig[2];
  TH1F *ptbal_bkg[2];

  TH1F *ptasym_sig[2];
  TH1F *ptasym_bkg[2];

  TH1F *sumpt2_sig[2];
  TH1F *sumpt2_bkg[2];

  int mycolor, mystyle;
  char hname[100];
  std::string evtType[2] = {"data","mc"};
  
  for (int i = 0; i < 2; i++){
    if (plotVariables){
      sumpt2_sig[i] = (TH1F*)f[i]->Get("sumpt2_sig");
      sumpt2_bkg[i] = (TH1F*)f[i]->Get("sumpt2_bkg");
      ptbal_sig[i] = (TH1F*)f[i]->Get("ptbal_sig");
      ptbal_bkg[i] = (TH1F*)f[i]->Get("ptbal_bkg");
      ptasym_sig[i] = (TH1F*)f[i]->Get("ptasym_sig");
      ptasym_bkg[i] = (TH1F*)f[i]->Get("ptasym_bkg");
      
      sumpt2_sig[i]->Sumw2();
      sumpt2_bkg[i]->Sumw2();
      ptbal_sig[i]->Sumw2();
      ptbal_bkg[i]->Sumw2();
      ptasym_sig[i]->Sumw2();
      ptasym_bkg[i]->Sumw2();

    }


    BDToutput_sig[i] = (TH1F*) f[i]->Get("BDToutput_sig");
    BDToutput_bkg[i] = (TH1F*) f[i]->Get("BDToutput_bkg");
    evtBDToutput[i]     = (TH1F*) f[i]->Get("perEventBDToutput");
    evtBDToutput_sig[i] = (TH1F*) f[i]->Get("perEventBDToutput_sig");
    evtBDToutput_bkg[i] = (TH1F*) f[i]->Get("perEventBDToutput_bkg");

    BDToutput_sig[i]->Sumw2();
    BDToutput_bkg[i]->Sumw2();
    evtBDToutput[i]->Sumw2();
    evtBDToutput_sig[i]->Sumw2();
    evtBDToutput_bkg[i]->Sumw2();
    
    BDToutput_sig[i]->Rebin(10);
    BDToutput_bkg[i]->Rebin(10);
    evtBDToutput[i]->Rebin(10);
    evtBDToutput_sig[i]->Rebin(10);
    evtBDToutput_bkg[i]->Rebin(10);

    NvtxAll[i]      = (TH1F*) f[i]->Get("NvtAll");
    NvtxBaseline[i] = (TH1F*) f[i]->Get("NvtGood_RANK");
    NvtxBDT[i]      = (TH1F*) f[i]->Get("NvtGood_BDT");

    PtAll[i]      = (TH1F*) f[i]->Get("PtAll");
    PtBaseline[i] = (TH1F*) f[i]->Get("PtGood_RANK");
    PtBDT[i]      = (TH1F*) f[i]->Get("PtGood_BDT");

    NvtxAll[i]     -> Sumw2();
    NvtxBaseline[i]-> Sumw2();
    NvtxBDT[i]     -> Sumw2();

    PtAll[i]       -> Sumw2();
    PtBaseline[i]  -> Sumw2();
    PtBDT[i]       -> Sumw2();

    sprintf(hname,"PtAllRebinned_%d",i);
    PtAllRebinned[i] = PtAll[i]   ->Rebin(9,hname,xbins);
    sprintf(hname,"PtBaselineRebinned_%d",i);
    PtBaselineRebinned[i] = PtBaseline[i]->Rebin(9,hname,xbins);
    sprintf(hname,"PtBDTRebinned_%d",i);
    PtBDTRebinned[i] = PtBDT[i]   ->Rebin(9,"PtBDTRebinned",xbins);

    sprintf(hname, "effVsNvtx_Baseline_%s",evtType[i].c_str());
    effVsNvtx_Baseline[i] = (TH1F*)NvtxBaseline[i]->Clone(hname);
    sprintf(hname, "effVsNvtx_BDT_%s",evtType[i].c_str());
    effVsNvtx_BDT[i]      = (TH1F*)NvtxBDT[i]->Clone(hname);

    if (useVariableBinning){
      sprintf(hname, "effVsPt_Baseline_%s",evtType[i].c_str());
      effVsPt_Baseline[i]   = (TH1F*)PtBaselineRebinned[i]->Clone(hname);
      sprintf(hname, "effVsPt_BDT_%s",evtType[i].c_str());
      effVsPt_BDT[i]        = (TH1F*)PtBDTRebinned[i]->Clone(hname); 
    }
    else {
      sprintf(hname, "effVsPt_Baseline_%s",evtType[i].c_str());
      effVsPt_Baseline[i]   = (TH1F*)PtBaseline[i]->Clone(hname);
      sprintf(hname, "effVsPt_BDT_%s",evtType[i].c_str());
      effVsPt_BDT[i]        = (TH1F*)PtBDT[i]->Clone(hname);
    }

    if (i==0) {mycolor = kBlue+1;}
    else {mycolor = kRed;}
    
    effVsNvtx_BDT[i]->SetLineColor(mycolor);
    effVsNvtx_BDT[i]->SetMarkerColor(mycolor);
    effVsNvtx_BDT[i]->SetMarkerStyle(20);

    effVsPt_BDT[i]->SetLineColor(mycolor);
    effVsPt_BDT[i]->SetMarkerColor(mycolor);
    effVsPt_BDT[i]->SetMarkerStyle(20);

    effVsNvtx_Baseline[i]->SetLineColor(mycolor);
    effVsNvtx_Baseline[i]->SetMarkerColor(mycolor);
    effVsNvtx_Baseline[i]->SetMarkerStyle(24);

    effVsPt_Baseline[i]->SetLineColor(mycolor);
    effVsPt_Baseline[i]->SetMarkerColor(mycolor);
    effVsPt_Baseline[i]->SetMarkerStyle(24);

  }

  
  for (int i = 0; i<2; i++){
    effVsNvtx_Baseline[i] -> Divide(NvtxBaseline[i],NvtxAll[i],1,1 ,"B"); 
    effVsNvtx_BDT[i]      -> Divide(NvtxBDT[i],NvtxAll[i],1,1, "B");
   
    if (!useVariableBinning){ 
      effVsPt_Baseline[i]   -> Divide(PtBaseline[i],PtAll[i],1,1, "B"); 
      effVsPt_BDT[i]        -> Divide(PtBDT[i],PtAll[i],1,1, "B");  
    }    
    else {
      effVsPt_Baseline[i]   -> Divide(PtBaselineRebinned[i],PtAllRebinned[i],1,1, "B"); 
      effVsPt_BDT[i]        -> Divide(PtBDTRebinned[i],PtAllRebinned[i],1,1, "B");  
    }
  }

  int nMBMax = 40;
  
  TLegend legend1(0.68, 0.18, 0.92, 0.38);
  legend1.SetFillColor(kWhite);
  legend1.SetBorderSize(1);
  //legend1.AddEntry(effVsNvtx_Baseline[0],legtitle1.c_str(),"LP");
  //legend1.AddEntry(effVsNvtx_Baseline[1],legtitle2.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_BDT[0],legtitle1bdt.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_BDT[1],legtitle2bdt.c_str(),"LP");

  //*** EFF vs NVTX
  TCanvas c1;
  c1.SetGridx();
  c1.SetGridy();
  TH2F cc("cc","",nMBMax+1,0,nMBMax+1,1000,0.,1.1);
  cc.SetStats(0); 
  cc.GetXaxis()->SetTitle("number of reconstructed vertices"); 
  cc.GetYaxis()->SetTitle("fraction of events");
  cc.Draw();
  for (int i = 0; i< 2; i++){
    //    effVsNvtx_Baseline[i]->Draw("e1,same");
    effVsNvtx_BDT[i]->Draw("e1,same");
  }
  legend1.Draw("same");
  latex->Draw("same");

  //*** EFF vs BOSON PT  
  TCanvas c2;
  c2.SetGridx();
  c2.SetGridy();
  TH2F dd("dd","",250,0,250,1000,0.0,1.1);
  dd.SetStats(0); 
  dd.GetXaxis()->SetTitle("p_{T}(Z) (GeV/c)"); 
  dd.GetYaxis()->SetTitle("fraction of events");
  dd.Draw();
  for (int i = 0; i< 2; i++){
    //effVsPt_Baseline[i]->Draw("e1,same");
    effVsPt_BDT[i]->Draw("e1,same");
  }
  legend1.Draw("same");
  latex->DrawLatex(0.18,0.22,text.c_str());
      
  float ptlow = 0.;
  float pthigh = 250.;

  int bin1 = PtBDT[0]->FindBin(ptlow);
  int bin2 = PtBDT[0]->FindBin(pthigh);
   
  double eff1 = PtBDT[0]->Integral(bin1,bin2)/ PtAll[0]->Integral(bin1,bin2);
  double eff2 = PtBDT[1]->Integral(bin1,bin2)/ PtAll[1]->Integral(bin1,bin2);
    
  cout << "Efficiency integrated in boson pt : [" << ptlow << ","<< pthigh<<"] GeV"<< endl;
  cout << legtitle1bdt.c_str() << " --> eff = " << eff1 << " +/- " 
       << sqrt(eff1*(1-eff1)/PtAll[0]->Integral(bin1,bin2))<< endl;  
  cout << legtitle2bdt.c_str() << " --> eff = " << eff2 << " +/- " 
       << sqrt(eff2*(1-eff2)/PtAll[1]->Integral(bin1,bin2))<<endl;  
    

  
  //EFFICIENCY RATIOs

  //------ RATIO EFF vs PT -------------------------------------------
  TH1F *ratioEffVsPt_Baseline = (TH1F*)effVsPt_Baseline[0]->Clone("ratioEffVsPt_Baseline");
  ratioEffVsPt_Baseline ->Divide(effVsPt_Baseline[1]);

  TH1F *ratioEffVsPt_BDT = (TH1F*)effVsPt_BDT[0]->Clone("ratioEffVsPt_BDT");
  ratioEffVsPt_BDT ->Divide(effVsPt_BDT[1]);
  
  TCanvas cRatioPt("cRatioPt","cRatioPt",600,300);
  cRatioPt.SetGridx();
  cRatioPt.SetGridy();

  ratioEffVsPt_BDT->SetMarkerColor(1);
  ratioEffVsPt_BDT->SetLineColor(1);
  ratioEffVsPt_BDT->GetXaxis()->SetTitle("p_{T}(Z) (GeV/c)");
  ratioEffVsPt_BDT->GetYaxis()->SetTitle("#epsilon(data)/#epsilon(MC)");
  ratioEffVsPt_BDT->GetYaxis()->SetTitleSize(0.07);
  ratioEffVsPt_BDT->GetYaxis()->SetTitleOffset(0.8);
  ratioEffVsPt_BDT->GetXaxis()->SetRangeUser(0,250);
  ratioEffVsPt_BDT->GetYaxis()->SetRangeUser(0.7,1.3);
  ratioEffVsPt_BDT->Draw("e1");
  //ratioEffVsPt_Baseline->Draw("e1same");
  latex2->Draw("same");
  latex->SetTextSize(0.06);
  latex->DrawLatex(0.67,0.82,text.c_str()); 
  //-- save also in TGraphErrors format
  TGraphErrors *gratioEffVsPt_BDT = new TGraphErrors();
  for (int ibin = 0; ibin < ratioEffVsPt_BDT->GetNbinsX();  ibin++){
    float x  = ratioEffVsPt_BDT->GetBinCenter(ibin+1);
    float ex = ratioEffVsPt_BDT->GetBinWidth(ibin+1)/2;
    float y  = ratioEffVsPt_BDT->GetBinContent(ibin+1);
    float ey = ratioEffVsPt_BDT->GetBinError(ibin+1);
    gratioEffVsPt_BDT->SetPoint(ibin,x,y);
    gratioEffVsPt_BDT->SetPointError(ibin,ex,ey);
  }

  TLegend legend4(0.15, 0.2, 0.45, 0.4);
  legend4.SetFillColor(kWhite);
  legend4.SetBorderSize(1);
  legend4.AddEntry(ratioEffVsPt_Baseline,"RANK","LP");
  legend4.AddEntry(ratioEffVsPt_BDT,"BDT","LP");
  //legend4.Draw("same");



  //------ RATIO EFF vs NVTX -------------------------------------------
  TH1F *ratioEffVsNvtx_Baseline = (TH1F*)effVsNvtx_Baseline[0]->Clone("ratioEffVsNvtx_Baseline");
  ratioEffVsNvtx_Baseline ->Divide(effVsNvtx_Baseline[1]);

  TH1F *ratioEffVsNvtx_BDT = (TH1F*)effVsNvtx_BDT[0]->Clone("ratioEffVsNvtx_BDT");
  ratioEffVsNvtx_BDT ->Divide(effVsNvtx_BDT[1]);

  TCanvas cRatioNvtx("cRatioNvtx","cRatioNvtx",600,300);
  cRatioNvtx.SetGridx();
  cRatioNvtx.SetGridy();

  ratioEffVsNvtx_BDT->SetMarkerColor(1);
  ratioEffVsNvtx_BDT->SetLineColor(1);
  ratioEffVsNvtx_BDT->GetXaxis()->SetTitle("number of reconstructed vertices");
  ratioEffVsNvtx_BDT->GetYaxis()->SetTitle("#epsilon(data)/#epsilon(MC)");
  ratioEffVsNvtx_BDT->GetYaxis()->SetTitleSize(0.07);
  ratioEffVsNvtx_BDT->GetYaxis()->SetTitleOffset(0.8);
  ratioEffVsNvtx_BDT->GetXaxis()->SetRangeUser(0,35);
  ratioEffVsNvtx_BDT->GetYaxis()->SetRangeUser(0.7,1.3);
  ratioEffVsNvtx_BDT->Draw("e1");
  //ratioEffVsNvtx_Baseline->Draw("e1");
  latex2->Draw("same");
  latex->SetTextSize(0.06);
  latex->DrawLatex(0.67,0.82,text.c_str());

  TLegend legend5(0.15, 0.2, 0.45, 0.4);
  legend5.SetFillColor(kWhite);
  legend5.SetBorderSize(1);
  legend5.AddEntry(ratioEffVsNvtx_Baseline,"RANK","LP");
  legend5.AddEntry(ratioEffVsNvtx_BDT,"BDT","LP");
  //legend5.Draw("same");



  // NVTX control plot
  TCanvas cNvtx("cNvtx","cNvtx",500,500);
  NvtxAll[1]->SetFillColor(kRed);
  NvtxAll[1]->SetFillStyle(3004);
  NvtxAll[1]->GetXaxis()->SetTitle("Number of vertices");
  NvtxAll[1]->DrawNormalized("histo");
  NvtxAll[0]->DrawNormalized("esame");
  latex->Draw("same");
  float p = NvtxAll[0]->KolmogorovTest(NvtxAll[1],"");
  cout << "DATA < N_vtx > = " << NvtxAll[0]->GetMean() << endl;
  cout << "MC   < N_vtx > = " << NvtxAll[1]->GetMean() << endl;
  cout << "NVTX Kolmogorov : p = " << p << endl;
  TLegend leg (0.6, 0.6,0.89,0.75);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry(NvtxAll[1],"MC Z#rightarrow#mu#mu","F");
  leg.AddEntry(NvtxAll[0],"Data Z#rightarrow#mu#mu","LP");
  leg.Draw("same");



  //----------------- MVA control plots
  // per vertex mva
  BDToutput_sig[0] -> SetMarkerStyle(20);
  BDToutput_sig[0] -> SetMarkerSize(0.8);
  BDToutput_sig[0] -> SetMarkerColor(kGreen+2);
  BDToutput_bkg[0] -> SetMarkerStyle(20);
  BDToutput_bkg[0] -> SetMarkerSize(0.8);
  BDToutput_bkg[0] -> SetMarkerColor(kRed+2);
  BDToutput_sig[1] -> SetFillColor(kGreen+1);
  BDToutput_sig[1] -> SetFillStyle(3002);
  BDToutput_bkg[1] -> SetFillColor(kRed+1);
  BDToutput_bkg[1] -> SetFillStyle(3005);

  TH1F *ratioBDToutput_sig = MakeRatio(BDToutput_sig[0], BDToutput_sig[1], "ratioBDToutput_sig");
  TH1F *ratioBDToutput_bkg = MakeRatio(BDToutput_bkg[0], BDToutput_bkg[1], "ratioBDToutput_bkg");
  
  TLegend legVtxMva (0.6, 0.6,0.89,0.75);
  legVtxMva.SetFillColor(0);
  legVtxMva.SetBorderSize(0);
  legVtxMva.AddEntry(evtBDToutput_sig[1],"right vertex MC","F");
  legVtxMva.AddEntry(evtBDToutput_bkg[1],"wrong vertex MC","F");
  legVtxMva.AddEntry(evtBDToutput_sig[0],"right vertex DATA","LP");
  legVtxMva.AddEntry(evtBDToutput_bkg[0],"wrong vertex DATA","LP");

  TCanvas cVertexMva("cVertexMva","cVertexMva",500,600);
  cVertexMva.Divide(1,2);
  cVertexMva.cd(1);
  BDToutput_bkg[1] -> GetXaxis() -> SetTitle("MVA_{vtx}");
  BDToutput_bkg[1] -> DrawNormalized("histo");
  BDToutput_sig[1] -> DrawNormalized("histo same");
  BDToutput_bkg[0] -> DrawNormalized("esame");
  BDToutput_sig[0] -> DrawNormalized("esame");
  legVtxMva.Draw("same");
  latex->Draw("same");
  cVertexMva.cd(2);
  cVertexMva.cd(2)->SetGridy();
  ratioBDToutput_sig -> GetYaxis()->SetTitle("data/MC");
  ratioBDToutput_sig -> GetYaxis()->SetRangeUser(0.,2.0);
  ratioBDToutput_sig ->Draw();
  ratioBDToutput_bkg ->Draw("same");

 
  // per event mva
  evtBDToutput_sig[0] -> SetMarkerStyle(20);
  evtBDToutput_sig[0] -> SetMarkerSize(0.8);
  evtBDToutput_sig[0] -> SetMarkerColor(kGreen+2);
  evtBDToutput_bkg[0] -> SetMarkerStyle(20);
  evtBDToutput_bkg[0] -> SetMarkerSize(0.8);
  evtBDToutput_bkg[0] -> SetMarkerColor(kRed+2);
  evtBDToutput_sig[1] -> SetFillColor(kGreen+1);
  evtBDToutput_sig[1] -> SetFillStyle(3002);
  evtBDToutput_bkg[1] -> SetFillColor(kRed+1);
  evtBDToutput_bkg[1] -> SetFillStyle(3005);

  TH1F *ratioevtBDToutput_sig = MakeRatio(evtBDToutput_sig[0], evtBDToutput_sig[1], "ratioevtBDToutput_sig");
  TH1F *ratioevtBDToutput_bkg = MakeRatio(evtBDToutput_bkg[0], evtBDToutput_bkg[1], "ratioevtBDToutput_bkg");

  TLegend legEvtMva (0.6, 0.6,0.89,0.75);
  legEvtMva.SetFillColor(0);
  legEvtMva.SetBorderSize(0);
  legEvtMva.AddEntry(evtBDToutput_sig[1],"right vertex MC","F");
  legEvtMva.AddEntry(evtBDToutput_bkg[1],"wrong vertex MC","F");
  legEvtMva.AddEntry(evtBDToutput_sig[0],"right vertex DATA","LP");
  legEvtMva.AddEntry(evtBDToutput_bkg[0],"wrong vertex DATA","LP");

  TCanvas cEventMva("cEventMva","cEventMva",500,600);
  cEventMva.Divide(1,2);
  cEventMva.cd(1);
  evtBDToutput_sig[1] -> GetXaxis() -> SetTitle("MVA_{event}");
  evtBDToutput_sig[1] -> DrawNormalized("histo");
  evtBDToutput_bkg[1] -> DrawNormalized("histo same");
  evtBDToutput_sig[0] -> DrawNormalized("esame");
  evtBDToutput_bkg[0] -> DrawNormalized("esame");
  legEvtMva.Draw("same");
  latex->Draw("same");
  cEventMva.cd(2);
  cEventMva.cd(2)->SetGridy();
  ratioevtBDToutput_sig -> GetYaxis()->SetTitle("data/MC");
  ratioevtBDToutput_sig -> GetYaxis()->SetRangeUser(0.,2.0);
  ratioevtBDToutput_sig ->Draw();
  ratioevtBDToutput_bkg ->Draw("same");

  // vertex probability
  TH1F *hVertexProbability[2];
  hVertexProbability[0]= (TH1F*)evtBDToutput_sig[0]->Clone("hVertexProbability");
  hVertexProbability[1]= (TH1F*)evtBDToutput_sig[1]->Clone("hVertexProbability");
  hVertexProbability[0]->Divide(evtBDToutput[0]);
  hVertexProbability[1]->Divide(evtBDToutput[1]);
  hVertexProbability[0]->SetMarkerColor(kBlue+1);
  hVertexProbability[1]->SetMarkerColor(kRed+1);

  TH1F *hRatioVertexProbability = (TH1F*)hVertexProbability[0]->Clone("hRatioVertexProbability");
  hRatioVertexProbability->Divide(hVertexProbability[1]);


  TF1 *fprob = new TF1("fprob","1+[0]*(x+1)",-1,1);
  fprob->SetParameter(0,-0.49);
  fprob->SetLineColor(kGray+1);
  TCanvas cProbability("cProbability","cProbability",600,600);
  cProbability.Divide(1,2);
  cProbability.cd(1);
  cProbability.cd(1)->SetGridx();
  cProbability.cd(1)->SetGridy();
  hVertexProbability[0]->GetXaxis()->SetTitle("MVA");
  hVertexProbability[0]->GetYaxis()->SetTitle("probability");
  hVertexProbability[0]->Draw();
  hVertexProbability[1]->Draw("same");
  fprob->Draw("same");
  legend1.Draw("same");

  cProbability.cd(2);
  cProbability.cd(2)->SetGridx();
  cProbability.cd(2)->SetGridy();
  hRatioVertexProbability->GetYaxis()->SetRangeUser(0.8,1.2);
  hRatioVertexProbability->GetYaxis()->SetTitle("data/MC");
  hRatioVertexProbability->GetXaxis()->SetTitle("MVA");
  hRatioVertexProbability->Draw();

 
  //input vtx id variables
  if (plotVariables){ 
    sumpt2_sig[0] -> SetMarkerStyle(20);
    sumpt2_sig[0] -> SetMarkerSize(0.8);
    sumpt2_sig[0] -> SetMarkerColor(kGreen+2);
    sumpt2_bkg[0] -> SetMarkerStyle(20);
    sumpt2_bkg[0] -> SetMarkerSize(0.8);
    sumpt2_bkg[0] -> SetMarkerColor(kRed+2);
    sumpt2_sig[1] -> SetFillColor(kGreen+1);
    sumpt2_sig[1] -> SetFillStyle(3002);
    sumpt2_bkg[1] -> SetFillColor(kRed+1);
    sumpt2_bkg[1] -> SetFillStyle(3005);

    TH1F *ratiosumpt2_sig = MakeRatio(sumpt2_sig[0], sumpt2_sig[1], "ratiosumpt2_sig");
    TH1F *ratiosumpt2_bkg = MakeRatio(sumpt2_bkg[0], sumpt2_bkg[1], "ratiosumpt2_bkg");
 
    TCanvas cSumpt2("cSumpt2","cSumpt2",500,600);
    cSumpt2.Divide(1,2);
    cSumpt2.cd(1);
    sumpt2_bkg[1] -> GetXaxis() -> SetTitle("log(sumpt2)");
    sumpt2_bkg[1] -> DrawNormalized("histo");
    sumpt2_sig[1] -> DrawNormalized("histo same");
    sumpt2_bkg[0] -> DrawNormalized("esame");
    sumpt2_sig[0] -> DrawNormalized("esame");
    legVtxMva.Draw("same");
    latex->Draw("same");
    cSumpt2.cd(2);
    cSumpt2.cd(2)->SetGridy();
    ratiosumpt2_sig -> GetYaxis()->SetTitle("data/MC");
    ratiosumpt2_sig -> GetYaxis()->SetRangeUser(0.,2.0);
    ratiosumpt2_sig ->Draw();
    ratiosumpt2_bkg ->Draw("same");

    

    ptbal_sig[0] -> SetMarkerStyle(20);
    ptbal_sig[0] -> SetMarkerSize(0.8);
    ptbal_sig[0] -> SetMarkerColor(kGreen+2);
    ptbal_bkg[0] -> SetMarkerStyle(20);
    ptbal_bkg[0] -> SetMarkerSize(0.8);
    ptbal_bkg[0] -> SetMarkerColor(kRed+2);
    ptbal_sig[1] -> SetFillColor(kGreen+1);
    ptbal_sig[1] -> SetFillStyle(3002);
    ptbal_bkg[1] -> SetFillColor(kRed+1);
    ptbal_bkg[1] -> SetFillStyle(3005);

    TH1F *ratioptbal_sig = MakeRatio(ptbal_sig[0], ptbal_sig[1], "ratioptbal_sig");
    TH1F *ratioptbal_bkg = MakeRatio(ptbal_bkg[0], ptbal_bkg[1], "ratioptbal_bkg");

    TCanvas cPtbal("cPtbal","cPtbal",500,600);
    cPtbal.Divide(1,2);
    cPtbal.cd(1);
    ptbal_bkg[1] -> GetXaxis() -> SetTitle("ptbal");
    ptbal_bkg[1] -> DrawNormalized("histo");
    ptbal_sig[1] -> DrawNormalized("histo same");
    ptbal_bkg[0] -> DrawNormalized("esame");
    ptbal_sig[0] -> DrawNormalized("esame");
    legVtxMva.Draw("same");
    latex->Draw("same");
    cPtbal.cd(2);
    cPtbal.cd(2)->SetGridy();
    ratioptbal_sig -> GetYaxis()->SetTitle("data/MC");
    ratioptbal_sig -> GetYaxis()->SetRangeUser(0.,2.0);
    ratioptbal_sig ->Draw();
    ratioptbal_bkg ->Draw("same");

    ptasym_sig[0] -> SetMarkerStyle(20);
    ptasym_sig[0] -> SetMarkerSize(0.8);
    ptasym_sig[0] -> SetMarkerColor(kGreen+2);
    ptasym_bkg[0] -> SetMarkerStyle(20);
    ptasym_bkg[0] -> SetMarkerSize(0.8);
    ptasym_bkg[0] -> SetMarkerColor(kRed+2);
    ptasym_sig[1] -> SetFillColor(kGreen+1);
    ptasym_sig[1] -> SetFillStyle(3002);
    ptasym_bkg[1] -> SetFillColor(kRed+1);
    ptasym_bkg[1] -> SetFillStyle(3005);

    TH1F *ratioptasym_sig = MakeRatio(ptasym_sig[0], ptasym_sig[1], "ratioptasym_sig");
    TH1F *ratioptasym_bkg = MakeRatio(ptasym_bkg[0], ptasym_bkg[1], "ratioptasym_bkg");

    TCanvas cPtasym("cPtasym","cPtasym",500,600);
    cPtasym.Divide(1,2);
    cPtasym.cd(1);
    ptasym_bkg[1] -> GetXaxis() -> SetTitle("ptasym");
    ptasym_bkg[1] -> DrawNormalized("histo");
    ptasym_sig[1] -> DrawNormalized("histo same");
    ptasym_bkg[0] -> DrawNormalized("esame");
    ptasym_sig[0] -> DrawNormalized("esame");
    legVtxMva.Draw("same");
    latex->Draw("same");
    cPtasym.cd(2);
    cPtasym.cd(2)->SetGridy();
    ratioptasym_sig -> GetYaxis()->SetTitle("data/MC");
    ratioptasym_sig -> GetYaxis()->SetRangeUser(0.,2.0);
    ratioptasym_sig ->Draw();
    ratioptasym_bkg ->Draw("same");
  }


  if (saveScaleFactors){
    TFile *fileout = new TFile("vtxIdScaleFactorFromZmumu_PUweights_minBiasXsec69400_observed_Run2012ABCD_BSreweight.root","recreate");
    ratioEffVsPt_BDT->SetTitle("hscaleFactor");
    ratioEffVsPt_BDT->Write("hscaleFactor");
    gratioEffVsPt_BDT->SetTitle("scaleFactor");
    gratioEffVsPt_BDT->Write("scaleFactor");
    fileout->Close();
    
    TFile *fileout = new TFile("vtxProbRatioFromZmumu_Run2012ABCD.root","recreate");
    hRatioVertexProbability->SetTitle("hRatioVertexProbability");
    hRatioVertexProbability->Write();
  }


  
  //--- save plots
  gSystem->mkdir(outdir.c_str(),true);
  gSystem->cd(outdir.c_str());
  cNvtx.SaveAs("nvtx_zmumu.png");
  c1.SaveAs("efficiency_vs_nvtx_zmumu.png");
  c2.SaveAs("efficiency_vs_pt_zmumu.png");
  cRatioNvtx.SaveAs("efficiencyRatio_vs_nvtx_zmumu.png");
  cRatioPt.SaveAs("efficiencyRatio_vs_pt_zmumu.png");
  cVertexMva.SaveAs("vertex_mva_zmumu.png");
  cEventMva.SaveAs("event_mva_zmumu.png");
  cProbability.SaveAs("vertex_probability_zmumu.png");

  cNvtx.SaveAs("nvtx_zmumu.pdf");
  c1.SaveAs("efficiency_vs_nvtx_zmumu.pdf");
  c2.SaveAs("efficiency_vs_pt_zmumu.pdf");
  cRatioNvtx.SaveAs("efficiencyRatio_vs_nvtx_zmumu.pdf");
  cRatioPt.SaveAs("efficiencyRatio_vs_pt_zmumu.pdf");
  cVertexMva.SaveAs("vertex_mva_zmumu.pdf");
  cEventMva.SaveAs("event_mva_zmumu.pdf");
  cProbability.SaveAs("vertex_probability_zmumu.pdf");

  gApplication->Run();


  string done;
  cout << "Done? " << endl;
  cin  >>  done ;

}
