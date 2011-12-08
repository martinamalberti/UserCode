#include "histoFunc.h"

// Template on MC: 1 per module
// Test on MC and data 
bool IsEtaGap(float eta){
  float feta = fabs(eta);
  if( fabs(feta - 0 )<2) return true;
  if( fabs(feta - 25)<2) return true;
  if( fabs(feta - 45)<2) return true;
  if( fabs(feta - 65)<2) return true;
  if( fabs(feta - 85)<2) return true;
  return false;
}

int templIndex(float eta){
    float feta = fabs(eta);
    if (feta <= 25)            {return 0;}
    if (feta> 25 && feta <= 45){return 1;}
    if (feta> 45 && feta <= 65){return 2;}
    if (feta> 65 && feta <= 85){return 3;}

    return -1;
}

void regionalEta()
{
  gROOT->SetStyle("Plain");

  float xtalWidth = 0.01745329;

  bool usePUweights = true;

  //--- weights for MC
  TFile weightsFile("weights/PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","READ"); // stessi pesi usati per analisi vertici Hgg  
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();


  // -- NTUPLES
  TChain *ntu_Data = new TChain("ntu");
  TChain *ntu_MC   = new TChain("ntu");


  //---- MC fall 2011
  ntu_MC->Add("/data2/calibrator/NTUPLES/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11_All.root");
   
  //---- DATA
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v6_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_data_20111122_158851_180363.root");


  std::cout << "     DATA: " << ntu_Data->GetEntries() << " entries in Data sample" << std::endl;
  std::cout << "     MC  : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;


  // observables
  int npu; 
  float EoP, scEta, scPhi, tkP;
  float scE3x3, scE5x5, scEne;  
  float charge, scLocalEta, scLocalPhi,crackCorr,localCorr; 

  float fbrem;
  float Etrue, dR;
  
  float sigP, kfP;
  // Set branch addresses for MC  
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_MC->SetBranchAddress("ele1_tkP", &tkP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_MC->SetBranchAddress("ele1_e5x5", &scE5x5);
  ntu_MC->SetBranchAddress("ele1_scE", &scEne);
  ntu_MC->SetBranchAddress("ele1_charge", &charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_MC->SetBranchAddress("ele1_fbrem",&fbrem); 

  ntu_MC->SetBranchAddress("ele1_sigmaP",&sigP); 
  ntu_MC->SetBranchAddress("ele1_KfP",&kfP); 

  // Set branch addresses for Data
  ntu_Data->SetBranchAddress("ele1_scEta", &scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_Data->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_Data->SetBranchAddress("ele1_tkP", &tkP);
  ntu_Data->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_Data->SetBranchAddress("ele1_e5x5", &scE5x5);
  ntu_Data->SetBranchAddress("ele1_scE", &scEne);
  ntu_Data->SetBranchAddress("ele1_charge", &charge);
  ntu_Data->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_Data->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_Data->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_Data->SetBranchAddress("ele1_fbrem",&fbrem); 

  ntu_Data->SetBranchAddress("ele1_sigmaP",&sigP); 
  ntu_Data->SetBranchAddress("ele1_KfP",&kfP); 


  

  unsigned int nBins = 85*2;
  std::cout << "nBins = " << nBins << std::endl;
  float etalimit = 1.44;
  float R9cut = 0.9;
  TProfile *pr = new TProfile("a","a",nBins,-1.*etalimit, etalimit);
  int rebin = 4;

  // MC
  TH1F** h_EoP_MC = new TH1F*[nBins];   
  TH1F** h_EoC_MC = new TH1F*[nBins];   
  

  for(unsigned int i = 0; i < nBins; ++i)
  {
    char histoName[80];
    sprintf(histoName, "EoP_MC_%d", i);
    h_EoP_MC[i] = new TH1F(histoName, histoName, 1200, 0., 3.);
    h_EoP_MC[i] -> SetFillColor(4);
    h_EoP_MC[i] -> SetFillStyle(3004);
    sprintf(histoName, "EoC_MC_%d", i);
    h_EoC_MC[i] = new TH1F(histoName, histoName, 1200, 0., 3.);
    h_EoC_MC[i] -> SetFillColor(3);
    h_EoC_MC[i] -> SetFillStyle(3004);
  }

  TH1F* h_template[4];
  for (unsigned int imod = 0; imod<4; imod++){
    char histoName[80];
    sprintf(histoName, "template_%d", imod);
    h_template[imod] = new TH1F(histoName, "", 1200, 0., 3.);  
  }
  std::cout << "Loop in MC events " << endl; 
  float ww = 1;
  // loop on MC, make refernce and fit dist
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%100000 == 0 ) std::cout << "reading MC saved entry " << entry << std::endl;
    // if (entry>1000) break;

    ntu_MC->GetEntry(entry);
      
    if (usePUweights) ww = w[npu];
    
    if( fabs(scEta) > etalimit ) continue;
    if (scE3x3/scEne <  R9cut ) continue; 
    
    //-- remove phi cracks 
    float phi = (scPhi+3.1415926536)/0.01745329;
    float modphi = (int)phi%20;
    if (fabs(modphi-10)<3.) continue;

    float fetaCry = fabs (scEta) / xtalWidth;
    int mod = templIndex(fetaCry);

    float var = EoP;
    
    // all the same reference
    //mod=0;
    // reference out of the gap
    if( IsEtaGap(fetaCry)==0 ) 
      h_template[mod] -> Fill(var*1,ww);

    // fill MC histos in eta bins
    int bin = pr->GetXaxis()->FindBin(scEta) - 1;
    //cout<<"bin"<<endl;
    h_EoP_MC[bin] -> Fill(var/crackCorr,ww);
    h_EoC_MC[bin] -> Fill(var,ww);
  } 

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////


 // Data
  TH1F** h_EoP_Data = new TH1F*[nBins];   
  TH1F** h_EoC_Data = new TH1F*[nBins];   
  for(unsigned int i = 0; i < nBins; ++i)
  {
    char histoName[80];
    sprintf(histoName, "EoP_Data_%d", i);
    h_EoP_Data[i] = new TH1F(histoName, histoName, 1200, 0., 3.);
    h_EoP_Data[i] -> SetFillColor(4);
    h_EoP_Data[i] -> SetFillStyle(3004);
    sprintf(histoName, "EoC_Data_%d", i);
    h_EoC_Data[i] = new TH1F(histoName, histoName, 1200, 0., 3.);
    h_EoC_Data[i] -> SetFillColor(3);
    h_EoC_Data[i] -> SetFillStyle(3004);
  }
  std::cout << "Loop in Data events " << endl; 
  // loop on Data
  for(int entry = 0; entry < ntu_Data->GetEntries(); ++entry) {
    if( entry%100000 == 0 ) std::cout << "reading data saved entry " << entry << std::endl;
    //    if (entry>1000) break;
    //cout<<"1111111111"<<endl;
    ntu_Data->GetEntry(entry);
    if( fabs(scEta) > etalimit ) continue;
    if (scE3x3/scEne <  R9cut ) continue; 
    
    //-- remove phi cracks
    float phi = (scPhi+3.1415926536)/0.01745329;
    float modphi = (int)phi%20;
    if (fabs(modphi-10)<3.) continue;

    float fetaCry = fabs (scEta) / xtalWidth;
    int mod = templIndex(fetaCry);
         
    float var = EoP;
    // fill Data histos in eta bins
    int bin = pr->GetXaxis()->FindBin(scEta) - 1;
    //cout<<"bin"<<endl;
    h_EoP_Data[bin] -> Fill(var/crackCorr);
    h_EoC_Data[bin] -> Fill(var);
  } 
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas *cc = new TCanvas("cc");
  cc->Divide(2,2); 
  histoFunc *templateHistoFunc[4]; 
  for (int mod=0; mod<4;mod++) {
    h_template[mod] -> Rebin(rebin);
    templateHistoFunc[mod] = new histoFunc(h_template[mod]);
    cc->cd(mod+1); 
    h_template[mod]->Draw(); 
  }
  /// Fitting  MC distribution ///////

  TGraphErrors* g_EoP_MC   = new TGraphErrors();g_EoP_MC->SetName("gEoP_MC");
  TGraphErrors* g_EoC_MC   = new TGraphErrors();g_EoC_MC->SetName("gEoC_MC");
  TGraphErrors* g_Rat_MC   = new TGraphErrors();g_Rat_MC->SetName("gCorr_Uncorr_MC");
  
  
  for(unsigned int i = 0; i < nBins; ++i)
  {
    h_EoP_MC[i] -> Rebin(rebin);    
    h_EoC_MC[i] -> Rebin(rebin);    

    float xval = pr->GetXaxis()->GetBinCenter(i+1);
    int ieta = fabs(xval)/0.01745329;
    int mod = 0;
    if (ieta>25) mod = 1;
    if (ieta>45) mod = 2;
    if (ieta>65) mod = 3;
  

    TF1 * templateFunc = new TF1("templateFunc", templateHistoFunc[mod], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    // uncorrected    
    double xNorm = h_EoP_MC[i]->GetEntries()/h_template[mod]->GetEntries() *
      h_EoP_MC[i]->GetBinWidth(1)/h_template[mod]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting uncorrected MC:  " << i << " " << mod; 
    int fStatus = h_EoP_MC[i] -> Fit("templateFunc", "MRQLNS+");
    g_EoP_MC -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_MC -> SetPointError(i, 2*etalimit/nBins, 10);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    // corrected    
    double xNorm = h_EoC_MC[i]->GetEntries()/h_template[mod]->GetEntries() *
      h_EoC_MC[i]->GetBinWidth(1)/h_template[mod]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting corrected MC " << i << " " << mod;
    fStatus = h_EoC_MC[i] -> Fit("templateFunc", "MRQLN+");
    g_EoC_MC -> SetPoint(i, xval , 1./templateFunc->GetParameter(1));
    g_EoC_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
     if(fStatus%10 == 4)  g_EoC_MC -> SetPointError(i, 2*etalimit/nBins, 10);
    cout << " ********** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio finalization
    rat /= templateFunc->GetParameter(1);
    era = rat*sqrt(era*era+templateFunc->GetParError(1)*templateFunc->GetParError(1)); 
    
    g_Rat_MC -> SetPoint(i,  xval , rat); 
    g_Rat_MC -> SetPointError(i,  0. , era); 
    g_Rat_MC->SetLineColor(kBlue+2); 

  }
 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //// fitting data distribution //////
   /// Fitting  MC distribution ///////

  TGraphErrors* g_EoP_Data   = new TGraphErrors();g_EoP_Data->SetName("gEoP_Data");
  TGraphErrors* g_EoC_Data   = new TGraphErrors();g_EoC_Data->SetName("gEoC_Data");
  TGraphErrors* g_Rat_Data   = new TGraphErrors();g_Rat_Data->SetName("gCorr_Uncorr_Data");
  
  
  for(unsigned int i = 0; i < nBins; ++i)
  {
    h_EoP_Data[i] -> Rebin(rebin);    
    h_EoC_Data[i] -> Rebin(rebin);    

    float xval = pr->GetXaxis()->GetBinCenter(i+1);
    int ieta = fabs(xval)/0.01745329;
    int mod = 0;
    if (ieta>25) mod = 1;
    if (ieta>45) mod = 2;
    if (ieta>65) mod = 3;
  

    TF1 * templateFunc = new TF1("templateFunc", templateHistoFunc[mod], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    // uncorrected    
    double xNorm = h_EoP_Data[i]->GetEntries()/h_template[mod]->GetEntries() *
      h_EoP_Data[i]->GetBinWidth(1)/h_template[mod]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting uncorrected Data:  " << i << " " << mod; 
    int fStatus = h_EoP_Data[i] -> Fit("templateFunc", "MRQLNS+");
    g_EoP_Data -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_Data -> SetPointError(i, 2*etalimit/nBins, 10);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    // corrected    
    double xNorm = h_EoC_Data[i]->GetEntries()/h_template[mod]->GetEntries() *
      h_EoC_Data[i]->GetBinWidth(1)/h_template[mod]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting corrected Data " << i << " " << mod;
    fStatus = h_EoC_Data[i] -> Fit("templateFunc", "MRQLN+");
    g_EoC_Data -> SetPoint(i, xval , 1./templateFunc->GetParameter(1));
    g_EoC_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
     if(fStatus%10 == 4)  g_EoC_Data -> SetPointError(i, 2*etalimit/nBins, 10);
    cout << " ********** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio finalization
    rat /= templateFunc->GetParameter(1);
    era = rat*sqrt(era*era+templateFunc->GetParError(1)*templateFunc->GetParError(1)); 
    
    g_Rat_Data -> SetPoint(i,  xval , rat); 
    g_Rat_Data -> SetPointError(i,  0. , era); 
    g_Rat_Data->SetLineColor(kBlue+2); 

  }
 

  /// Plotting plots


  TCanvas* c_g_fit_MC = new TCanvas("g_fit_MC", "g_fit_MC",100,100,700,500);
  c_g_fit_MC->cd(1);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(-1.*etalimit,0.93,etalimit,1.03);
  hPad->GetXaxis()->SetTitle("#eta_{SC} (deg)");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  g_EoP_MC -> SetMarkerStyle(20);
  g_EoP_MC -> SetMarkerSize(.7);
  g_EoP_MC -> SetMarkerColor(kRed+1); 
  g_EoP_MC -> Draw("PL");
  g_EoC_MC -> SetMarkerStyle(20);
  g_EoC_MC -> SetMarkerSize(.7);
  g_EoC_MC -> SetMarkerColor(kGreen+1); 
  g_EoC_MC -> Draw("PL");

  //g_Rat_MC->Draw("L"); 

  
  TLegend *tl_MC = new TLegend(0.80,0.85,1.01,1.01);
  tl_MC -> SetFillColor(0);
  //  tl_MC -> SetBorderSize(0); 
  tl_MC -> AddEntry(g_EoP_MC,"MC Uncorr","P");
  tl_MC -> AddEntry(g_EoC_MC,"MC corr","P");
  tl_MC -> Draw();
  
  //////////////////////////////
  TCanvas* c_g_fit_Data = new TCanvas("g_fit_Data", "g_fit_Data",100,100,700,500);
  c_g_fit_Data->cd(1);
  gPad->SetGrid();
  TH1F *hPad = (TH1F*)gPad->DrawFrame(-1.*etalimit,0.93,etalimit,1.03);
  hPad->GetXaxis()->SetTitle("#eta_{SC} (deg)");
  hPad->GetYaxis()->SetTitle("Relative E/p scale"); 
  g_EoP_Data -> SetMarkerStyle(20);
  g_EoP_Data -> SetMarkerSize(.7);
  g_EoP_Data -> SetMarkerColor(kRed+1); 
  g_EoP_Data -> Draw("PL");
  g_EoC_Data -> SetMarkerStyle(20);
  g_EoC_Data -> SetMarkerSize(.7);
  g_EoC_Data -> SetMarkerColor(kGreen+1); 
  g_EoC_Data -> Draw("PL");

  //g_Rat_Data->Draw("L"); 

  
  TLegend *tl_Data = new TLegend(0.80,0.85,1.01,1.01);
  tl_Data -> SetFillColor(0);
  //  tl_Data -> SetBorderSize(0); 
  tl_Data -> AddEntry(g_EoP_Data,"Data Uncorr","P");
  tl_Data -> AddEntry(g_EoC_Data,"Data corr","P");
  tl_Data -> Draw();
  


  /// outfile /////
  TFile outfile("regionalScale.root","recreate");

  g_EoP_MC->Write();
  g_EoC_MC->Write();
  g_Rat_MC->Write();

  g_EoP_Data->Write();
  g_EoC_Data->Write();
  g_Rat_Data->Write();

  for (unsigned int imod = 0; imod<4; imod++){
    h_template[imod]->Write();  
  }

}
