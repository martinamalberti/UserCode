// per compilare: g++ -Wall -o EoPPhiFolded `root-config --cflags --glibs` EoPPhiFolded.cpp

#include "../CommonTools/histoFunc.h"
#include "../CommonTools/TPileupReweighting.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

bool IsEtaGap(float eta){
  float feta = fabs(eta);
  if( fabs(feta - 0 )<3) return true;
  if( fabs(feta - 25)<3) return true;
  if( fabs(feta - 45)<3) return true;
  if( fabs(feta - 65)<3) return true;
  if( fabs(feta - 85)<3) return true;
  return false;
}

int main(int argc, char** argv){
  
  //--- folding
  int cryFold = 2;
  
  //---- output file to save graphs
  char outfilename[100];
  //sprintf(outfilename,"Graphs_EoP_vs_LocalPhi_Folded20_regression_etacharge_pos.root");
  sprintf(outfilename,"test_RelVal_V14.root");
    
  bool usePUweights = false;

  //---- PU weights for MC
  TPileupReweighting *puReweighting;
  puReweighting = new TPileupReweighting("../CommonTools/weights/PUweights_2011_0100_73500_WJetsToLL_Fall11_S6.root","hweights"); 
  
  //--- NTUPLES 
  TChain *ntu_MC = new TChain("ntu");
  TChain *ntu_Data = new TChain("ntu");
  
  //--- MC Fall 2011
  //  ntu_MC->Add("../NTUPLES/Fall11/WZAnalysis/WZAnalysis_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1.root");
  ntu_MC->Add("../NTUPLES/RelVal/WZAnalysis_RelValZee_V14.root");

  //---- DATA
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");

  std::cout << "     MC  : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;

  // observables
  float EoP, scEta, scPhi, tkP;
  float scE3x3, scE5x5, scE, scE_regression;  
  float charge, scLocalEta, scLocalPhi,crackCorr,localCorr; 
  int npu;
  float R9;
  int seedIphi;

  // Set branch addresses for MC 
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("ele1_seedIphi", &seedIphi);
  ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_MC->SetBranchAddress("ele1_tkP", &tkP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_MC->SetBranchAddress("ele1_e5x5", &scE5x5);
  ntu_MC->SetBranchAddress("ele1_scE", &scE);
  ntu_MC->SetBranchAddress("ele1_scE_regression", &scE_regression);
  ntu_MC->SetBranchAddress("ele1_charge", &charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 

  // Set branch addresses for Data
  ntu_Data->SetBranchAddress("ele1_scEta", &scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_Data->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_Data->SetBranchAddress("ele1_tkP", &tkP);
  ntu_Data->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_Data->SetBranchAddress("ele1_e5x5", &scE5x5);
  ntu_Data->SetBranchAddress("ele1_scE", &scE);
  ntu_Data->SetBranchAddress("ele1_scE_regression", &scE_regression);
  ntu_Data->SetBranchAddress("ele1_charge", &charge);
  ntu_Data->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_Data->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_Data->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 

  
  float def_error = 0.;
  unsigned int nBins = 10;
  std::cout << "nBins = " << nBins << std::endl;
  float etalimit = 0.8;
  float R9cut = 0.00;
  TProfile *pr = new TProfile("prPhiFolded","Correction profile Phi",nBins,0,cryFold,0.7,1.25);
  
  int rebin = 8;

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

  // -- define different templates for data/MC, corrected uncorrected
  // [0][0] --> mc uncorrected
  // [1][0] --> data uncorrected
  // [0][1] --> mc corrected
  // [1][1] --> data corrected

  TH1F* h_template[2][2];
  h_template[0][0] = new TH1F("templateMC", "", 1200, 0., 3.);  
  h_template[1][0] = new TH1F("templateDATA", "", 1200, 0., 3.);  
  h_template[0][1] = new TH1F("templateMCcorr", "", 1200, 0., 3.);  
  h_template[1][1] = new TH1F("templateDATAcorr", "", 1200, 0., 3.);    

  TH2F *htest = new TH2F("htest","htest", 360,1,361,20,0,20);
  TH1F *htest2 = new TH1F("htest2","htest2", 360,1,361);
  
  std::cout << "Loop in MC events " << std::endl; 

  //**** LOOP ON MC
  float ww = 1;
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%100000 == 0 ) std::cout << "reading MC saved entry " << entry << std::endl;
    //if (entry>10000) break;
    ntu_MC->GetEntry(entry);

    // -- PU weights
    if (usePUweights) ww = puReweighting->GetWeight(npu);
 
    //--- R9, eta, charge selections 
    // if( fabs(scEta) > etalimit ) continue;
    // if (scE3x3/scE <  R9cut ) continue; 
    // float etacharge = scEta * charge; 
    // if(etacharge <0) continue;

    //-- remove eta gaps
    float fetaCry = fabs (scEta) / 0.01745329;
    if( IsEtaGap(fetaCry) ) continue;
    
    float myphi = scPhi;
    if (scEta<0) myphi = -myphi; 
    myphi = (myphi+3.1415926536)/0.01745329;
    int modphi = (int)myphi%20;

    //---skip phi gap if folded inside a SM 
    if (cryFold<20. && fabs(modphi-10)<2.) continue;
    
    float foldphi = (int)myphi%cryFold + scLocalPhi + 0.5;
    // if (cryFold == 20) 
    //  foldphi =  (int)myphi%cryFold+myphi- (int)myphi ; // non si capisce perche' usare questo per fold=20


    // -- to correct every second crystal
    // if (scEta > 0){
    //   if (seedIphi%2==1) EoP*=0.998;
    //   if (seedIphi%2==0) EoP*=1.002;
    // } 
    // if (scEta < 0){
    //   if (seedIphi%2==0) EoP*=0.998;
    //   if (seedIphi%2==1) EoP*=1.002;
    // } 
    
    //if (int(myphi)%2==0) EoP*=0.998;
    //if (int(myphi)%2==1) EoP*=1.002;

    //--- reference out of the gap
    float correction = scE_regression/scE;
    if ( fabs(modphi-10)>2. )        {
      h_template[0][0]-> Fill(EoP,ww);
      h_template[0][1]-> Fill(EoP*correction,ww);
    }

    //--- fill MC histos in eta bins
    int bin = pr->GetXaxis()->FindBin(foldphi) - 1;
    if (fabs(scLocalPhi>=0.5) ) continue;
   
    // htest->Fill(seedIphi, myphi,ww);    
    // if (int(foldphi)%2==0) htest2->Fill(seedIphi,ww);
        
    
    h_EoP_MC[bin] -> Fill(EoP,ww);
    h_EoC_MC[bin] -> Fill(EoP*correction,ww);
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

  std::cout << "Loop in Data events " << std::endl; 

  // loop on Data
  for(int entry = 0; entry < ntu_Data->GetEntries(); ++entry) {
    if( entry%100000 == 0 ) std::cout << "reading data saved entry " << entry << std::endl;
    if (entry>50000) break;
    
    ntu_Data->GetEntry(entry);
    
    // -- R9, eta , charge selections
    if( fabs(scEta) > etalimit ) continue;
    if (scE3x3/scE <  R9cut ) continue; 
    float etacharge = scEta * charge; 
    if(etacharge <0) continue;

    //-- remove eta gap
    float fetaCry = fabs (scEta) / 0.01745329;
    if( IsEtaGap(fetaCry) ) continue;
    
    float myphi = scPhi;
    if (scEta<0) myphi = -myphi;
    myphi = (myphi+3.1415926536)/0.01745329;
    int modphi = (int)myphi%20;
  
    //-- skip phi gap if folded inside a SM 
    if (cryFold<20. && fabs(modphi-10)<2.) continue;
    
    float foldphi = (int)myphi%cryFold + scLocalPhi + 0.5;
    if (cryFold == 20) 
      foldphi =  (int)myphi%cryFold+myphi- (int)myphi ; // non si capisce perche' usare questo per fold=20
    
    //-- reference out of the gap
    float correction  = scE_regression/scE;
    if ( fabs(modphi-10)>2. )        {
      h_template[1][0]-> Fill(EoP);
      h_template[1][1]-> Fill(EoP*correction);
    }

    //-- fill Data histos in eta bins
    int bin = pr->GetXaxis()->FindBin(foldphi) - 1;
    h_EoP_Data[bin] -> Fill(EoP);
    h_EoC_Data[bin] -> Fill(EoP*correction);
  } 
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas *cc = new TCanvas("cc");

  histoFunc *templateHistoFunc[2][2]; 
  for (int jj = 0; jj < 2; jj++){
    for (int ll = 0; ll < 2; ll++){
      h_template[jj][ll]-> Rebin(rebin);
      templateHistoFunc[jj][ll] = new histoFunc(h_template[jj][ll]);
    }  
  }
    
  ///---- Fitting  MC distribution ///////
  TGraphErrors* g_EoP_MC   = new TGraphErrors();g_EoP_MC->SetName("gEoP_MC");
  TGraphErrors* g_EoC_MC   = new TGraphErrors();g_EoC_MC->SetName("gEoC_MC");
  TGraphErrors* g_Rat_MC   = new TGraphErrors();g_Rat_MC->SetName("gCorr_Uncorr_MC");
    
  for(unsigned int i = 0; i < nBins; ++i)
  {
    h_EoP_MC[i] -> Rebin(rebin);    
    h_EoC_MC[i] -> Rebin(rebin);    

    float xval = pr->GetXaxis()->GetBinCenter(i+1);
     
    TF1 *templateFunc = new TF1("templateFunc", templateHistoFunc[0][0], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    //--- uncorrected    
    double xNorm = h_EoP_MC[i]->Integral()/h_template[0][0]->Integral() *
                   h_EoP_MC[i]->GetBinWidth(1)/h_template[0][0]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting uncorrected MC:  " << i << " "; 
    int fStatus = h_EoP_MC[i] -> Fit("templateFunc", "MRQLNS+");
    g_EoP_MC -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_MC -> SetPointError(i, def_error, def_error);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    //cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    //--- corrected 
    templateFunc = new TF1("templateFunc", templateHistoFunc[0][1], 0.7, 1.3, 3, "histoFunc");
  
    xNorm = h_EoC_MC[i]->Integral()/h_template[0][1]->Integral() *
            h_EoC_MC[i]->GetBinWidth(1)/h_template[0][1]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting corrected MC " << i << " ";
    fStatus = h_EoC_MC[i] -> Fit("templateFunc", "MRQLN+");
    g_EoC_MC -> SetPoint(i, xval , 1./templateFunc->GetParameter(1));
    g_EoC_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoC_MC -> SetPointError(i,def_error, def_error);
    //cout << " ********** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio finalization
    rat /= templateFunc->GetParameter(1);
    era = rat*sqrt(era*era+templateFunc->GetParError(1)*templateFunc->GetParError(1)); 
    
    g_Rat_MC -> SetPoint(i,  xval , rat); 
    g_Rat_MC -> SetPointError(i,  0. , era); 
    g_Rat_MC -> SetLineColor(kBlue+2); 

  }
 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ////**** fitting data distribution //////
  
  TGraphErrors* g_EoP_Data   = new TGraphErrors();g_EoP_Data->SetName("gEoP_Data");
  TGraphErrors* g_EoC_Data   = new TGraphErrors();g_EoC_Data->SetName("gEoC_Data");
  TGraphErrors* g_Rat_Data   = new TGraphErrors();g_Rat_Data->SetName("gCorr_Uncorr_Data");
    
  for(unsigned int i = 0; i < nBins; ++i)
  {
    h_EoP_Data[i] -> Rebin(rebin);    
    h_EoC_Data[i] -> Rebin(rebin);    

    float xval = pr->GetXaxis()->GetBinCenter(i+1);
  

    TF1 * templateFunc = new TF1("templateFunc", templateHistoFunc[1][0], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    // --- uncorrected    
    double xNorm = h_EoP_Data[i]->Integral()/h_template[1][0]->Integral() *
                   h_EoP_Data[i]->GetBinWidth(1)/h_template[1][0]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting uncorrected Data:  " << i <<" "; 
    int fStatus = h_EoP_Data[i] -> Fit("templateFunc", "MRQLNS+");
    g_EoP_Data -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_Data -> SetPointError(i,def_error , def_error);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    //cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    // ---- corrected    
    templateFunc = new TF1("templateFunc", templateHistoFunc[1][1], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);
    xNorm = h_EoC_Data[i]->Integral()/h_template[1][1]->Integral() *
            h_EoC_Data[i]->GetBinWidth(1)/h_template[1][1]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting corrected Data " << i << "  ";
    fStatus = h_EoC_Data[i] -> Fit("templateFunc", "MRQLN+");
    g_EoC_Data -> SetPoint(i, xval , 1./templateFunc->GetParameter(1));
    g_EoC_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoC_Data -> SetPointError(i, def_error, def_error);
    //cout << " ********** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio finalization
    rat /= templateFunc->GetParameter(1);
    era = rat*sqrt(era*era+templateFunc->GetParError(1)*templateFunc->GetParError(1)); 
    
    g_Rat_Data -> SetPoint(i,  xval , rat); 
    g_Rat_Data -> SetPointError(i,  0. , era); 
    g_Rat_Data->SetLineColor(kBlue+2); 

  }

  g_EoP_MC-> SetMarkerStyle(20);
  g_EoC_MC-> SetMarkerStyle(20);
  g_EoP_Data-> SetMarkerStyle(20);
  g_EoC_Data-> SetMarkerStyle(20);
 
  g_EoP_MC-> SetMarkerColor(kRed);
  g_EoC_MC-> SetMarkerColor(kRed+2);
  g_EoP_Data-> SetMarkerColor(kGreen);
  g_EoC_Data-> SetMarkerColor(kGreen+2);

  /// outfile /////
  TFile outfile(outfilename,"recreate");

  g_EoP_MC->Write();
  g_EoC_MC->Write();
  g_Rat_MC->Write();

  g_EoP_Data->Write();
  g_EoC_Data->Write();
  g_Rat_Data->Write();
  
  h_template[0][0]->Write();  
  h_template[1][0]->Write();  
  h_template[0][1]->Write();  
  h_template[1][1]->Write();  
  
  htest->Write();
  htest2->Write();
  
}
