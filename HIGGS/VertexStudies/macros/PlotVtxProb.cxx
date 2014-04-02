void PlotVtxProb(string outdir) 
{

  //gROOT->LoadMacro("~/setTDRStyle.C");
  //setTDRStyle();
  //gStyle->SetErrorX(0.5);
  //gStyle->SetLegendFont(42);
  //gStyle->SetTextFont(42);
  //gStyle->SetTextSize(0.04);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);



  // 0 : data
  // 1 : mc

  TFile *f[2];

  //2012 rereco
  f[0] = TFile::Open("./Efficiency_DoubleMu_22JanReReco_2012_v1/testEfficiency.root");
  f[1] = TFile::Open("./Efficiency_DYJetsToLL_2012_pixelcorr_v1/testEfficiency.root");
  //2011 rereco
  //f[0] = TFile::Open("../lxbatch_scripts/Efficiency_DoubleMu_2011/testEfficiency.root");
  //f[1] = TFile::Open("../lxbatch_scripts/Efficiency_DYJetsToLL_2011/testEfficiency.root");

  std::string text = "CMS  #sqrt{s} = 8 TeV  L = 19.7 fb^{-1}";

  TLatex *latex = new TLatex(0.10,0.92,text.c_str());
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextFont(42);

  TH1F *evtBDToutput_sig[2];
  TH1F *evtBDToutput_bkg[2];
  TH1F *vtxProbability_sig[2];
  TH1F *vtxProbability_bkg[2];

  int mycolor, mystyle;
  char hname[100];
  std::string evtType[2] = {"data","mc"};
  
  for (int i = 0; i < 2; i++){
    evtBDToutput_sig[i] = (TH1F*) f[i]->Get("perEventBDToutput_sig");
    evtBDToutput_bkg[i] = (TH1F*) f[i]->Get("perEventBDToutput_bkg");
    evtBDToutput_sig[i]->Sumw2();
    evtBDToutput_bkg[i]->Sumw2();
    evtBDToutput_sig[i]->Rebin(10);
    evtBDToutput_bkg[i]->Rebin(10);

    vtxProbability_sig[i] = (TH1F*) f[i]->Get("vtxProbability_sig");
    vtxProbability_bkg[i] = (TH1F*) f[i]->Get("vtxProbability_bkg");
    vtxProbability_sig[i]->Sumw2();
    vtxProbability_bkg[i]->Sumw2();
    vtxProbability_sig[i]->Rebin(5);
    vtxProbability_bkg[i]->Rebin(5);
  }

  cout << vtxProbability_sig[0]->GetBinWidth(1) <<endl;

  //----------------- MVA control plots
  // per event mva
  evtBDToutput_sig[0] -> SetMarkerStyle(21);
  evtBDToutput_sig[0] -> SetMarkerSize(1);
  evtBDToutput_bkg[0] -> SetMarkerStyle(20);
  evtBDToutput_bkg[0] -> SetMarkerSize(1);
  evtBDToutput_sig[1] -> SetFillColor(kYellow);
  TColor *color = new TColor(1756, 0.0, 1.0, 0.0, "", 0.4); 
  evtBDToutput_bkg[1] -> SetFillColor(1756);

  TLegend legmva (0.5, 0.6,0.89,0.92);
  legmva.SetFillColor(0);
  legmva.SetBorderSize(0);
  legmva.AddEntry(evtBDToutput_sig[1],"Simulation: Right vertex","F");
  legmva.AddEntry(evtBDToutput_bkg[1],"Simulation: Wrong vertex","F");
  legmva.AddEntry(evtBDToutput_sig[0],"Data: Right vertex ","LP");
  legmva.AddEntry(evtBDToutput_bkg[0],"Data: Wrong vertex ","LP");

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
  legmva.Draw("same");
  latex->Draw("same");
  

  // vertex probability
  vtxProbability_sig[0] -> SetMarkerStyle(21);
  vtxProbability_sig[0] -> SetMarkerSize(1);
  vtxProbability_bkg[0] -> SetMarkerStyle(20);
  vtxProbability_bkg[0] -> SetMarkerSize(1);
  vtxProbability_sig[1] -> SetFillColor(kYellow);
  vtxProbability_bkg[1] -> SetFillColor(1756);

  //mio
  //TLegend legprob (0.15, 0.6,0.54,0.92);
  //legprob.SetFillColor(0);
  //legprob.SetBorderSize(0);
  //nancy
  TLegend legprob (0.15,0.5,0.6,0.85,"","brNDC");
  legprob.SetBorderSize(0);
  legprob.SetFillColor(0);
  legprob.SetTextFont(42);
  legprob.SetTextSize(0.035);
  legprob.SetTextAlign(12);
  legprob.AddEntry(evtBDToutput_sig[1],"Simulation: right vertex","F");
  legprob.AddEntry(evtBDToutput_bkg[1],"Simulation: wrong vertex","F");
  legprob.AddEntry(evtBDToutput_sig[0],"Data: right vertex ","P");
  legprob.AddEntry(evtBDToutput_bkg[0],"Data: wrong vertex ","P");

  //  TCanvas cVtxProb("cVtxProb","cVtxProb",500,500);
  TCanvas cVtxProb("cVtxProb","cVtxProb",430, 10, 700,700);
  cVtxProb.cd(1);
  
  vtxProbability_sig[1] -> GetXaxis() -> SetTitle("Vertex probability estimate");
  vtxProbability_sig[1] -> GetYaxis() -> SetTitle("Events #times 10^{4}/0.02");
  vtxProbability_sig[1] -> GetYaxis()->SetTitleSize(0.04);
  vtxProbability_sig[1] -> GetXaxis()->SetTitleSize(0.04);
  vtxProbability_sig[1] -> GetXaxis()->SetTitleOffset(1.25);//
  vtxProbability_sig[1] -> GetYaxis()->SetTitleOffset(1.28);
  vtxProbability_sig[1] -> GetYaxis()->SetTitleFont(42);
  vtxProbability_sig[1] -> GetXaxis()->SetTitleFont(42);
  vtxProbability_sig[1] -> GetXaxis()->SetNdivisions(405);
  vtxProbability_sig[1] -> GetXaxis()->SetLabelSize(0.05);
  vtxProbability_sig[1] -> GetYaxis()->SetLabelSize(0.05);
  vtxProbability_sig[1] -> GetXaxis()->SetLabelFont(42);
  vtxProbability_sig[1] -> GetYaxis()->SetLabelFont(42);
  //vtxProbability_sig[1] -> GetXaxis() -> SetLabelSize(0.04);
  //vtxProbability_sig[1] -> GetYaxis() -> SetLabelSize(0.04);
  //vtxProbability_sig[1] -> GetXaxis() -> SetTitleSize(0.04);
  //vtxProbability_sig[1] -> GetYaxis() -> SetTitleSize(0.04);
  //vtxProbability_sig[1] -> GetXaxis() -> SetTitleOffset(1.3);
  //vtxProbability_sig[1] -> GetYaxis() -> SetTitleOffset(1.4);
  //vtxProbability_sig[1] -> GetXaxis()->SetNdivisions(505);

     
  vtxProbability_sig[0] -> Scale(1./10000);
  vtxProbability_sig[1] -> Scale(1./10000);
  vtxProbability_bkg[0] -> Scale(1./10000);
  vtxProbability_bkg[1] -> Scale(1./10000);

  totData = vtxProbability_sig[0]->GetSumOfWeights() + vtxProbability_bkg[0]->GetSumOfWeights();
  totMC   = vtxProbability_sig[1]->GetSumOfWeights() + vtxProbability_bkg[1]->GetSumOfWeights();
  vtxProbability_sig[1] -> DrawNormalized("histo", vtxProbability_sig[1] ->GetSumOfWeights()*totData/totMC); //normalize to data entries                                                                                                      
  vtxProbability_sig[0] -> Draw("esame");
  vtxProbability_bkg[1] -> DrawNormalized("histo same", vtxProbability_bkg[1] ->GetSumOfWeights()*totData/totMC);
  vtxProbability_bkg[0] -> Draw("esame");
  vtxProbability_sig[0] -> Draw("esame");
  legprob.Draw("same");
  latex->Draw("same");
  
  //--- save plots
  gSystem->mkdir(outdir.c_str(),true);
  gSystem->cd(outdir.c_str());
  cEventMva.SaveAs("event_mva_zmumu.png");
  cEventMva.SaveAs("event_mva_zmumu.pdf");
  cEventMva.SaveAs("event_mva_zmumu.C");
  cEventMva.SaveAs("event_mva_zmumu.root");
  cVtxProb.SaveAs("vtxprob_zmumu.png");
  cVtxProb.SaveAs("vtxprob_zmumu.pdf");
  cVtxProb.SaveAs("vtxprob_zmumu.C");
  cVtxProb.SaveAs("vtxprob_zmumu.root");
 
  //  gApplication->Run();


  string done;
  cout << "Done? " << endl;
  cin  >>  done ;

}
