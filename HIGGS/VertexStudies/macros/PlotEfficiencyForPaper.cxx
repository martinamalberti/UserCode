TH1F* MakeRatio(TH1F* h1, TH1F*h2, string hname){

  TH1F *hratio = (TH1F*)h1->Clone(hname.c_str());
  hratio->Divide(h1,h2,1./h1->GetSumOfWeights(),1./h2->GetSumOfWeights());
  return(hratio);
}



void PlotEfficiencyForPaper(string outdir, string sqrts, string lumi) 
{

  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();
  //gStyle->SetErrorX(0.5);
  gStyle->SetLegendFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.04);

  // 0 : data
  // 1 : mc

  TFile *f[2];

  //2012 rereco
  f[0] = TFile::Open("../lxbatch_scripts/Efficiency_DoubleMu_22JanReReco_2012/testEfficiency.root");
  f[1] = TFile::Open("../lxbatch_scripts/Efficiency_DYJetsToLL_2012_pixelcorr/testEfficiency.root");
  //2011 rereco
  //f[0] = TFile::Open("../lxbatch_scripts/Efficiency_DoubleMu_2011/testEfficiency.root");
  //f[1] = TFile::Open("../lxbatch_scripts/Efficiency_DYJetsToLL_2011/testEfficiency.root");

  std::string text = "CMS #sqrt{s} = "+sqrts+" TeV L = "+lumi+" fb^{-1}";

  TLatex *latex = new TLatex(0.35,0.96,text.c_str());
  latex->SetNDC();

  TH1F *BDToutput[2];
  TH1F *BDToutput_sig[2];
  TH1F *BDToutput_bkg[2];

  TH1F *evtBDToutput[2];
  TH1F *evtBDToutput_sig[2];
  TH1F *evtBDToutput_bkg[2];

  int mycolor, mystyle;
  char hname[100];
  std::string evtType[2] = {"data","mc"};
  
  for (int i = 0; i < 2; i++){
    BDToutput_sig[i]    = (TH1F*) f[i]->Get("BDToutput_sig");
    BDToutput_bkg[i]    = (TH1F*) f[i]->Get("BDToutput_bkg");
    evtBDToutput[i]     = (TH1F*) f[i]->Get("perEventBDToutput");
    evtBDToutput_sig[i] = (TH1F*) f[i]->Get("perEventBDToutput_sig");
    evtBDToutput_bkg[i] = (TH1F*) f[i]->Get("perEventBDToutput_bkg");

    BDToutput_sig[i]   ->Sumw2();
    BDToutput_bkg[i]   ->Sumw2();
    evtBDToutput[i]    ->Sumw2();
    evtBDToutput_sig[i]->Sumw2();
    evtBDToutput_bkg[i]->Sumw2();
    
    BDToutput_sig[i]   ->Rebin(10);
    BDToutput_bkg[i]   ->Rebin(10);
    evtBDToutput[i]    ->Rebin(10);
    evtBDToutput_sig[i]->Rebin(10);
    evtBDToutput_bkg[i]->Rebin(10);
  }

  //----------------- MVA control plots
  // per event mva
  evtBDToutput_sig[0] -> SetMarkerStyle(20);
  evtBDToutput_sig[0] -> SetMarkerSize(1);
  evtBDToutput_bkg[0] -> SetMarkerStyle(21);
  evtBDToutput_bkg[0] -> SetMarkerSize(0.8);
  evtBDToutput_sig[1] -> SetFillColor(kYellow);
  //evtBDToutput_sig[1] -> SetFillStyle(3001);
  evtBDToutput_bkg[1] -> SetFillColor(kGreen);
  evtBDToutput_bkg[1] -> SetFillStyle(3005);

  TLegend legVtxMva (0.5, 0.6,0.89,0.92);
  legVtxMva.SetFillColor(0);
  legVtxMva.SetBorderSize(0);
  legVtxMva.AddEntry(evtBDToutput_sig[1],"Simulation: Right vertex","F");
  legVtxMva.AddEntry(evtBDToutput_bkg[1],"Simulation: Wrong vertex","F");
  legVtxMva.AddEntry(evtBDToutput_sig[0],"Data: Right vertex ","LP");
  legVtxMva.AddEntry(evtBDToutput_bkg[0],"Data: Wrong vertex ","LP");

  TCanvas cEventMva("cEventMva","cEventMva",500,500);
  cEventMva.cd(1);
  evtBDToutput_sig[1] -> GetXaxis() -> SetTitle("BDT score");
  evtBDToutput_sig[1] -> GetYaxis() -> SetTitle("Events/0.04");
  evtBDToutput_sig[1] -> GetXaxis() -> SetLabelSize(0.04);
  evtBDToutput_sig[1] -> GetYaxis() -> SetLabelSize(0.04);
  evtBDToutput_sig[1] -> GetXaxis() -> SetTitleSize(0.04);
  evtBDToutput_sig[1] -> GetYaxis() -> SetTitleSize(0.04);
  evtBDToutput_sig[1] -> GetXaxis() -> SetTitleOffset(1.3);
  evtBDToutput_sig[1] -> GetYaxis() -> SetTitleOffset(1.4);

  float totData = evtBDToutput_sig[0]->GetSumOfWeights() + evtBDToutput_bkg[0]->GetSumOfWeights();
  float totMC   = evtBDToutput_sig[1]->GetSumOfWeights() + evtBDToutput_bkg[1]->GetSumOfWeights();
  evtBDToutput_sig[1] -> DrawNormalized("histo", evtBDToutput_sig[1] ->GetSumOfWeights()*totData/totMC); //normalize to data entries
  evtBDToutput_sig[0] -> Draw("esame");
  evtBDToutput_bkg[1] -> DrawNormalized("histo same", evtBDToutput_bkg[1] ->GetSumOfWeights()*totData/totMC);
  evtBDToutput_bkg[0] -> Draw("esame");
  evtBDToutput_sig[0] -> Draw("esame");
  legVtxMva.Draw("same");
  latex->Draw("same");
  
  //--- save plots
  gSystem->mkdir(outdir.c_str(),true);
  gSystem->cd(outdir.c_str());
  cEventMva.SaveAs("event_mva_zmumu.png");
  cEventMva.SaveAs("event_mva_zmumu.pdf");
  cEventMva.SaveAs("event_mva_zmumu.C");
  cEventMva.SaveAs("event_mva_zmumu.root");

  gApplication->Run();


  string done;
  cout << "Done? " << endl;
  cin  >>  done ;

}
