TH1F* MakeRatio(TH1F* h1, TH1F*h2, string hname){

  TH1F *hratio = (TH1F*)h1->Clone(hname.c_str());
  hratio->Divide(h1,h2,1./h1->GetSumOfWeights(),1./h2->GetSumOfWeights());
  return(hratio);
}



void PlotVertexMvaInVertexBins(int vtxbin, string outdir) 
{

  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();
  gStyle->SetErrorX(0.5);

  // 0 : data
  // 1 : mc

  TFile *f[2];

  //2012ABCD
  f[0] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DoubleMu_Run2012ABCD_vtxbins/testEfficiency.root");
  f[1] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012ABCD_vtxbins/testEfficiency.root");
  
  //std::string text = "#splitline{CMS preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}";
  std::string text = "CMS preliminary #sqrt{s} = 8 TeV L = 19.6 fb^{-1}";
  
  TLatex *latex = new TLatex(0.15,0.96,text.c_str());
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.05);

  // std::string text2 = "10 < N_{vtx} < 20 ";
  std::string text2 = "N_{vtx} > 20 ";
  TLatex *latex2 = new TLatex(0.20,0.85,text2.c_str());
  latex2->SetNDC();
  latex2->SetTextFont(42);
  latex2->SetTextSize(0.06);

  TH1F *vtxBDToutput_sig[2];
  TH1F *vtxBDToutput_bkg[2];

  TH1F *evtBDToutput[2];
  TH1F *evtBDToutput_sig[2];
  TH1F *evtBDToutput_bkg[2];

  int mycolor, mystyle;
  char hname[100];
  std::string evtType[2] = {"data","mc"};
  
  for (int i = 0; i < 2; i++){
    
    vtxBDToutput_sig[i] = (TH1F*) f[i]->Get(Form("perVertexBDToutput_sig_vtxbin%d",vtxbin));
    vtxBDToutput_bkg[i] = (TH1F*) f[i]->Get(Form("perVertexBDToutput_bkg_vtxbin%d",vtxbin));
    evtBDToutput_sig[i] = (TH1F*) f[i]->Get(Form("perEventBDToutput_sig_vtxbin%d",vtxbin));
    evtBDToutput_bkg[i] = (TH1F*) f[i]->Get(Form("perEventBDToutput_bkg_vtxbin%d",vtxbin));
    evtBDToutput[i] = (TH1F*) f[i]->Get(Form("perEventBDToutput_vtxbin%d",vtxbin));

    vtxBDToutput_sig[i]->Sumw2();
    vtxBDToutput_bkg[i]->Sumw2();
    evtBDToutput_sig[i]->Sumw2();
    evtBDToutput_bkg[i]->Sumw2();
    evtBDToutput[i]->Sumw2();

    vtxBDToutput_sig[i]->Rebin(10);
    vtxBDToutput_bkg[i]->Rebin(10);
    evtBDToutput_sig[i]->Rebin(10);
    evtBDToutput_bkg[i]->Rebin(10);
    evtBDToutput[i]->Rebin(10);
  }

  //----------------- MVA control plots
  //---- per vertex mva
  vtxBDToutput_sig[0] -> SetMarkerStyle(20);
  vtxBDToutput_sig[0] -> SetMarkerSize(0.8);
  vtxBDToutput_sig[0] -> SetMarkerColor(kGreen+2);
  vtxBDToutput_bkg[0] -> SetMarkerStyle(20);
  vtxBDToutput_bkg[0] -> SetMarkerSize(0.8);
  vtxBDToutput_bkg[0] -> SetMarkerColor(kRed+2);
  vtxBDToutput_sig[1] -> SetFillColor(kGreen+1);
  vtxBDToutput_sig[1] -> SetFillStyle(3002);
  vtxBDToutput_bkg[1] -> SetFillColor(kRed+1);
  vtxBDToutput_bkg[1] -> SetFillStyle(3005);

  TH1F *ratioBDToutput_sig = MakeRatio(vtxBDToutput_sig[0], vtxBDToutput_sig[1], "ratioBDToutput_sig");
  TH1F *ratioBDToutput_bkg = MakeRatio(vtxBDToutput_bkg[0], vtxBDToutput_bkg[1], "ratioBDToutput_bkg");
  
  TLegend legVtxMva (0.6, 0.7,0.89,0.89);
  legVtxMva.SetFillColor(0);
  legVtxMva.SetBorderSize(0);
  legVtxMva.AddEntry(vtxBDToutput_sig[1],"right vertex MC","F");
  legVtxMva.AddEntry(vtxBDToutput_bkg[1],"wrong vertex MC","F");
  legVtxMva.AddEntry(vtxBDToutput_sig[0],"right vertex DATA","LP");
  legVtxMva.AddEntry(vtxBDToutput_bkg[0],"wrong vertex DATA","LP");

  TCanvas cVertexMva("cVertexMva","cVertexMva",500,600);
  cVertexMva.Divide(1,2);
  cVertexMva.cd(1);
  vtxBDToutput_bkg[1] -> GetXaxis() -> SetTitle("MVA_{vtx}");
  vtxBDToutput_bkg[1] -> GetYaxis() -> SetTitle("a.u.");
  vtxBDToutput_bkg[1] -> DrawNormalized("histo");
  vtxBDToutput_sig[1] -> DrawNormalized("histo same");
  vtxBDToutput_bkg[0] -> DrawNormalized("esame");
  vtxBDToutput_sig[0] -> DrawNormalized("esame");
  legVtxMva.Draw("same");
  latex->Draw("same");
  latex2->Draw("same");
  cVertexMva.cd(2);
  cVertexMva.cd(2)->SetGridy();
  ratioBDToutput_sig -> GetYaxis()->SetTitle("data/MC");
  ratioBDToutput_sig -> GetYaxis()->SetRangeUser(0.,2.0);
  ratioBDToutput_sig ->Draw();
  ratioBDToutput_bkg ->Draw("same");

 
  //---- per event mva
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
  evtBDToutput_sig[1] -> GetYaxis() -> SetTitle("a.u.");
  evtBDToutput_sig[1] -> DrawNormalized("histo");
  evtBDToutput_bkg[1] -> DrawNormalized("histo same");
  evtBDToutput_sig[0] -> DrawNormalized("esame");
  evtBDToutput_bkg[0] -> DrawNormalized("esame");
  legEvtMva.Draw("same");
  latex->Draw("same");
  latex2->Draw("same");
  cEventMva.cd(2);
  cEventMva.cd(2)->SetGridy();
  ratioevtBDToutput_sig -> GetXaxis()->SetTitle("MVA_{event}");
  ratioevtBDToutput_sig -> GetYaxis()->SetTitle("data/MC");
  ratioevtBDToutput_sig -> GetYaxis()->SetRangeUser(0.,2.0);
  ratioevtBDToutput_sig ->Draw();
  ratioevtBDToutput_bkg ->Draw("same");


  //--- vertex probability
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
  hVertexProbability[0]->GetYaxis()->SetRangeUser(0.0,1.2);
  hVertexProbability[0]->GetXaxis()->SetTitle("MVA");
  hVertexProbability[0]->GetYaxis()->SetTitle("probability");
  hVertexProbability[0]->Draw();
  hVertexProbability[1]->Draw("same");
  //fprob->Draw("same");
  latex2->Draw("same");
  TLegend legend1(0.68, 0.68, 0.92, 0.88);
  legend1.SetFillColor(kWhite);
  legend1.SetBorderSize(1);
  legend1.AddEntry(hVertexProbability[0],"data","LP");
  legend1.AddEntry(hVertexProbability[1],"MC","LP");
  legend1.Draw("same");



  cProbability.cd(2);
  cProbability.cd(2)->SetGridx();
  cProbability.cd(2)->SetGridy();
  hRatioVertexProbability->GetYaxis()->SetRangeUser(0.8,1.2);
  hRatioVertexProbability->GetYaxis()->SetTitle("data/MC");
  hRatioVertexProbability->GetXaxis()->SetTitle("MVA");
  hRatioVertexProbability->Draw();


  //--- save plots
  gSystem->mkdir(outdir.c_str(),true);
  gSystem->cd(outdir.c_str());
  cVertexMva.SaveAs("vertex_mva_zmumu.png");
  cEventMva.SaveAs("event_mva_zmumu.png");
  cProbability.SaveAs("vertex_probability_zmumu.png");
  cVertexMva.SaveAs("vertex_mva_zmumu.pdf");
  cEventMva.SaveAs("event_mva_zmumu.pdf");
  cProbability.SaveAs("vertex_probability_zmumu.pdf");
  gApplication->Run();


  string done;
  cout << "Done? " << endl;
  cin  >>  done ;

}
