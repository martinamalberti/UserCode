#include "histoFunc.h"

// Template on MC: 1 per module
// Test on MC and data 
bool IsEtaGap(float eta){
  float feta = fabs(eta);
  if( fabs(feta - 0 )<3) return true;
  if( fabs(feta - 25)<3) return true;
  if( fabs(feta - 45)<3) return true;
  if( fabs(feta - 65)<3) return true;
  if( fabs(feta - 85)<3) return true;
  return false;
}

void EoPPhiFolded(int cryFold=20)
{

  TChain *ntu_MC = new TChain("ntu");
  TChain *ntu_Data = new TChain("ntu");

  ntu_MC->Add("/data2/calibrator/NTUPLES/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11_All.root");

  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_data_20111122_158851_180363.root"); 
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_data_20111122_158851_180363.root"); 
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_data_20111122_158851_180363.root"); 
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v6_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_data_20111122_158851_180363.root");

  std::cout << "     MC  : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;

  // observables
  float EoP, scEta, scPhi, tkP;
  float scE3x3, scE5x5, scEne;  
  float charge, scLocalEta, scLocalPhi,crackCorr,localCorr; 

  float fbrem;
  float Etrue, dR;
  
  float sigP, kfP;
  // Set branch addresses for MC  
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


  
  float def_error = 0.;
  unsigned int nBins = 100;
  std::cout << "nBins = " << nBins << std::endl;
  float etalimit = 0.8;
  float R9cut = 0.009;
  TProfile *pr = new TProfile("prPhiFolded","Correction profile Phi",nBins,0,cryFold,0.7,1.25);
  int rebin = 2;

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

  TH1F* h_template[2];
  h_template[0] = new TH1F("templateMC", "", 1200, 0., 3.);  
  h_template[1] = new TH1F("templateDATA", "", 1200, 0., 3.);  
  
  
  std::cout << "Loop in MC events " << endl; 
  // loop on MC, make refernce and fit dist
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%100000 == 0 ) std::cout << "reading MC saved entry " << entry << std::endl;
    //if (entry>500000) break;

    ntu_MC->GetEntry(entry);
    if( fabs(scEta) > etalimit ) continue;
    if (scE3x3/scEne <  R9cut ) continue; 
    float etacharge = scEta * charge; 
    if(etacharge >0) continue;

    //-- remove eta gaps
    float fetaCry = fabs (scEta) / 0.01745329;
    if( IsEtaGap(fetaCry) ) continue;
    
    float myphi = scPhi;
    if (scEta<0) myphi = -myphi;
    myphi = (myphi+3.1415926536)/0.01745329;
    int modphi = (int)myphi%20;
    //skip phi gap if folded inside a SM 
    if (cryFold<20. && fabs(modphi-10)<2.) continue;
    

    //float foldphi =  (int)myphi%cryFold + scLocalPhi + 0.5;
    float foldphi =  (int)myphi%cryFold+myphi- (int)myphi ; // non si capisce 
    float var = EoP;

    // reference out of the gap
    if ( fabs(modphi-10)>2. )        {
      h_template[0]-> Fill(var);
    }

    // fill MC histos in eta bins
    int bin = pr->GetXaxis()->FindBin(foldphi) - 1;
    //cout<<"bin"<<endl;
    h_EoP_MC[bin] -> Fill(var/crackCorr);
    h_EoC_MC[bin] -> Fill(var);
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
    //if (entry>50000) break;
    
    ntu_Data->GetEntry(entry);
     if( fabs(scEta) > etalimit ) continue;
    if (scE3x3/scEne <  R9cut ) continue; 
    float etacharge = scEta * charge; 
    if(etacharge >0) continue;

    //-- remove eta gaps
    float fetaCry = fabs (scEta) / 0.01745329;
    if( IsEtaGap(fetaCry) ) continue;
    
    float myphi = scPhi;
    if (scEta<0) myphi = -myphi;
    myphi = (myphi+3.1415926536)/0.01745329;
    int modphi = (int)myphi%20;
    //skip phi gap if folded inside a SM 
    if (cryFold<20. && fabs(modphi-10)<2.) continue;
    
    //float foldphi =  (int)myphi%cryFold + scLocalPhi + 0.5;
    float foldphi =  (int)myphi%cryFold+myphi- (int)myphi ; // non si capisce  
    float var = EoP;
    
    // reference out of the gap
    if ( fabs(modphi-10)>2. )        {
      h_template[1]-> Fill(var);
    }
    // fill Data histos in eta bins
    int bin = pr->GetXaxis()->FindBin(foldphi) - 1;
    //cout<<"bin"<<endl;
    h_EoP_Data[bin] -> Fill(var/crackCorr);
    h_EoC_Data[bin] -> Fill(var);
  } 

  //////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas *cc = new TCanvas("cc");

  histoFunc *templateHistoFunc[2]; 
  for (int jj = 0; jj < 2; jj++){
    h_template[jj]-> Rebin(rebin);
    templateHistoFunc[jj] = new histoFunc(h_template[jj]);
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
     
    TF1 *templateFunc = new TF1("templateFunc", templateHistoFunc[0], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    // uncorrected    
    double xNorm = h_EoP_MC[i]->GetEntries()/h_template[0]->GetEntries() *
                   h_EoP_MC[i]->GetBinWidth(1)/h_template[0]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting uncorrected MC:  " << i << " "; 
    int fStatus = h_EoP_MC[i] -> Fit("templateFunc", "MRQLNS+");
    g_EoP_MC -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_MC -> SetPointError(i, def_error, def_error);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    // corrected    
    double xNorm = h_EoC_MC[i]->GetEntries()/h_template[0]->GetEntries() *
                   h_EoC_MC[i]->GetBinWidth(1)/h_template[0]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting corrected MC " << i << " ";
    fStatus = h_EoC_MC[i] -> Fit("templateFunc", "MRQLN+");
    g_EoC_MC -> SetPoint(i, xval , 1./templateFunc->GetParameter(1));
    g_EoC_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoC_MC -> SetPointError(i,def_error, def_error);
    cout << " ********** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio finalization
    rat /= templateFunc->GetParameter(1);
    era = rat*sqrt(era*era+templateFunc->GetParError(1)*templateFunc->GetParError(1)); 
    
    g_Rat_MC -> SetPoint(i,  xval , rat); 
    g_Rat_MC -> SetPointError(i,  0. , era); 
    g_Rat_MC -> SetLineColor(kBlue+2); 

  }
 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //// fitting data distribution //////
  
  TGraphErrors* g_EoP_Data   = new TGraphErrors();g_EoP_Data->SetName("gEoP_Data");
  TGraphErrors* g_EoC_Data   = new TGraphErrors();g_EoC_Data->SetName("gEoC_Data");
  TGraphErrors* g_Rat_Data   = new TGraphErrors();g_Rat_Data->SetName("gCorr_Uncorr_Data");
    
  for(unsigned int i = 0; i < nBins; ++i)
  {
    h_EoP_Data[i] -> Rebin(rebin);    
    h_EoC_Data[i] -> Rebin(rebin);    

    float xval = pr->GetXaxis()->GetBinCenter(i+1);
  

    TF1 * templateFunc = new TF1("templateFunc", templateHistoFunc[1], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    // uncorrected    
    double xNorm = h_EoP_Data[i]->GetEntries()/h_template[1]->GetEntries() *
                   h_EoP_Data[i]->GetBinWidth(1)/h_template[1]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting uncorrected Data:  " << i <<" "; 
    int fStatus = h_EoP_Data[i] -> Fit("templateFunc", "MRQLNS+");
    g_EoP_Data -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_Data -> SetPointError(i,def_error , def_error);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    // corrected    
    double xNorm = h_EoC_Data[i]->GetEntries()/h_template[1]->GetEntries() *
      h_EoC_Data[i]->GetBinWidth(1)/h_template[1]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting corrected Data " << i << "  ";
    fStatus = h_EoC_Data[i] -> Fit("templateFunc", "MRQLN+");
    g_EoC_Data -> SetPoint(i, xval , 1./templateFunc->GetParameter(1));
    g_EoC_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
     if(fStatus%10 == 4)  g_EoC_Data -> SetPointError(i, def_error, def_error);
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
  TH1F *hPad = (TH1F*)gPad->DrawFrame(0,0.93,cryFold,1.03);
  hPad->GetXaxis()->SetTitle("#phi_{SC} (deg)");
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
  TH1F *hPad = (TH1F*)gPad->DrawFrame(0,0.93,cryFold,1.03);
  hPad->GetXaxis()->SetTitle("#phi_{SC} (deg)");
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
  TFile outfile("SaveMe.root","recreate");

  g_EoP_MC->Write();
  g_EoC_MC->Write();
  g_Rat_MC->Write();

  g_EoP_Data->Write();
  g_EoC_Data->Write();
  g_Rat_Data->Write();

  
  h_template[0]->Write();  
  h_template[1]->Write();  
  

}
