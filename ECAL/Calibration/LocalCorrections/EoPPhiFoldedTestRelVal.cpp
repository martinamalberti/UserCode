// per compilare: g++ -Wall -o EoPPhiFoldedTestRelVal `root-config --cflags --glibs` EoPPhiFoldedTestRelVal.cpp

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
  unsigned int nBins = 8;

  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"test_RelVal_V14A.root");
    
 
  //--- NTUPLES 
  TChain *ntu_MC = new TChain("ntu");
  TChain *ntu_Data = new TChain("ntu");
  
  //--- MC Fall 2011
  ntu_MC->Add("../NTUPLES/Fall11/WZAnalysis/WZAnalysis_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1.root");
  
  //---- DATA
  ntu_Data->Add("../NTUPLES/RelVal/WZAnalysis_RelValZee_V14A.root");

  std::cout << "     MC  : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;

  // observables
  float EoP, scEta, scPhi, tkP;
  float scE3x3, scE5x5, scE, scE_regression;  
  float charge, scLocalEta, scLocalPhi,crackCorr,localCorr; 
  int npu;
  float R9;
  int seedIphi;
  int isZ;

  float EoP2, scEta2, scPhi2, charge2;
  float scLocalPhi2;

  // Set branch addresses for MC 
  ntu_MC->SetBranchAddress("isZ", &isZ);
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

  ntu_MC->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_MC->SetBranchAddress("ele2_scPhi", &scPhi2);
  ntu_MC->SetBranchAddress("ele2_EOverP", &EoP2);
  ntu_MC->SetBranchAddress("ele2_scLocalPhi",&scLocalPhi2); 
  ntu_MC->SetBranchAddress("ele2_charge", &charge2);

  // Set branch addresses for Data
  ntu_Data->SetBranchAddress("isZ", &isZ);
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

  ntu_Data->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_Data->SetBranchAddress("ele2_scPhi", &scPhi2);
  ntu_Data->SetBranchAddress("ele2_EOverP", &EoP2);
  ntu_Data->SetBranchAddress("ele2_scLocalPhi",&scLocalPhi2); 
  ntu_Data->SetBranchAddress("ele2_charge", &charge2);

  float def_error = 1.;
  std::cout << "nBins = " << nBins << std::endl;
  float etalimit = 1.4442;
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
    //if (entry>1000000) break;
    ntu_MC->GetEntry(entry);

    // if (npu > 3) continue;

    if (isZ==0) continue;

    //--- R9, eta, charge selections 
    //if( fabs(scEta) > etalimit ) continue;
    //if( fabs(scEta2) > etalimit ) continue;
    // if (scE3x3/scE <  R9cut ) continue; 
    // float etacharge = scEta * charge; 
    // if(etacharge <0) continue;

    // fisrt electron
    //-- remove eta gaps
    float fetaCry = fabs (scEta) / 0.01745329;
        
    float myphi = scPhi;
    if (scEta<0) myphi = -myphi; 
    myphi = (myphi+3.1415926536)/0.01745329;
    int modphi = (int)myphi%20;
    
    //if ( fabs(scLocalPhi)<0.5 && IsEtaGap(fetaCry)==0 &&  fabs(modphi-10)>2 && fabs(scEta) < etalimit && (scEta * charge) > 0){
    if ( fabs(scLocalPhi)<0.5 && IsEtaGap(fetaCry)==0 &&  fabs(modphi-10)>2 && fabs(scEta) < etalimit ){
    
      float foldphi = (int)myphi%cryFold + scLocalPhi + 0.5;
      
      //--- reference 
      h_template[0][0]-> Fill(EoP,ww);
      
      
      //--- fill MC histos in eta bins
      int bin = pr->GetXaxis()->FindBin(foldphi) - 1;
            
      h_EoP_MC[bin] -> Fill(EoP,ww);
      h_EoC_MC[bin] -> Fill(EoP,ww);
    }

    
    // -- second electron
    //-- remove eta gaps
    float fetaCry2 = fabs (scEta2) / 0.01745329;
    
    float myphi2 = scPhi2;
    if (scEta2<0) myphi2 = -myphi2; 
    myphi2 = (myphi2+3.1415926536)/0.01745329;
    int modphi2 = (int)myphi2%20;
    
    if ( fabs(scLocalPhi2)< 0.5 && IsEtaGap(fetaCry2)==0 &&  fabs(modphi2-10)>2 && fabs(scEta2) < etalimit ){
      //      if ( fabs(scLocalPhi2)< 0.5 && IsEtaGap(fetaCry2)==0 &&  fabs(modphi2-10)>2 && fabs(scEta2) < etalimit && (scEta2 * charge2) >0 ){
      
      float foldphi2 = (int)myphi2%cryFold + scLocalPhi2 + 0.5;
      
      //--- reference 
      h_template[0][0]-> Fill(EoP2,ww);


      //--- fill MC histos in eta bins
      int bin2 = pr->GetXaxis()->FindBin(foldphi2) - 1;
           
      h_EoP_MC[bin2] -> Fill(EoP2,ww);
      h_EoC_MC[bin2] -> Fill(EoP2,ww);
    }

   
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
      
    ntu_Data->GetEntry(entry);
    
    if (isZ==0) continue;
 
    //--- R9, eta, charge selections 
    // if( fabs(scEta) > etalimit ) continue;
    // if (scE3x3/scE <  R9cut ) continue; 
    // float etacharge = scEta * charge; 
    // if(etacharge <0) continue;

    // fisrt electron
    //-- remove eta gaps
    float fetaCry = fabs (scEta) / 0.01745329;
        
    float myphi = scPhi;
    if (scEta<0) myphi = -myphi; 
    myphi = (myphi+3.1415926536)/0.01745329;
    int modphi = (int)myphi%20;
    

    if ( fabs(scLocalPhi)< 0.5 && IsEtaGap(fetaCry)==0 &&  fabs(modphi-10)>2 && fabs(scEta) < etalimit ){
      //if ( fabs(scLocalPhi)< 0.5 && IsEtaGap(fetaCry)==0 &&  fabs(modphi-10)>2 && fabs(scEta) < etalimit && (scEta * charge) > 0){

      float foldphi = (int)myphi%cryFold + scLocalPhi + 0.5;
      
      //--- reference
	h_template[1][0]-> Fill(EoP,ww);

      
      //--- fill MC histos in eta bins
      int bin = pr->GetXaxis()->FindBin(foldphi) - 1;
      if (fabs(scLocalPhi)>=0.5 ) continue;
      
      h_EoP_Data[bin] -> Fill(EoP,ww);
      h_EoC_Data[bin] -> Fill(EoP,ww);
    }


    // -- second electron
    //-- remove eta gaps
    float fetaCry2 = fabs (scEta2) / 0.01745329;
    
    float myphi2 = scPhi2;
    if (scEta2<0) myphi2 = -myphi2; 
    myphi2 = (myphi2+3.1415926536)/0.01745329;
    int modphi2 = (int)myphi2%20;
    
    if ( fabs(scLocalPhi2)< 0.5 && IsEtaGap(fetaCry2)==0 &&  fabs(modphi2-10)>2 && fabs(scEta2) < etalimit ){
      //if ( fabs(scLocalPhi2)< 0.5 && IsEtaGap(fetaCry2)==0 &&  fabs(modphi2-10)>2 && fabs(scEta2) < etalimit && (scEta2 * charge2) > 0){
      
      float foldphi2 = (int)myphi2%cryFold + scLocalPhi2 + 0.5;
      
      //--- reference
      h_template[1][0]-> Fill(EoP2,ww);
	

      //--- fill MC histos in eta bins
      int bin2 = pr->GetXaxis()->FindBin(foldphi2) - 1;
      if (fabs(scLocalPhi2)>=0.5 ) continue;
      
      h_EoP_Data[bin2] -> Fill(EoP2,ww);
      h_EoC_Data[bin2] -> Fill(EoP2,ww);
    }

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
 
    TFitResultPtr rp;
    int fStatus; 
    for (int trial=0;trial<5;trial++) {
      templateFunc -> SetParameter(1, 0.99 + trial*0.01 );
      rp = h_EoP_MC[i] -> Fit("templateFunc", "MQRLS+");
      fStatus = rp;
      if (fStatus < 2) break; 
    } 
        
    g_EoP_MC -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_MC -> SetPointError(i, 0, def_error);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    //cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    //--- corrected 
    templateFunc = new TF1("templateFunc", templateHistoFunc[0][0], 0.7, 1.3, 3, "histoFunc");
  
    xNorm = h_EoC_MC[i]->Integral()/h_template[0][0]->Integral() *
            h_EoC_MC[i]->GetBinWidth(1)/h_template[0][0]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);

    std::cout << "***** Fitting corrected MC " << i << " " << std::endl;
    for (int trial=0;trial<5;trial++) {
      templateFunc -> SetParameter(1, 0.99 + trial*0.01 );
      rp = h_EoC_MC[i] -> Fit("templateFunc", "MQRLS+");
      fStatus = rp;
      if (fStatus < 2) break; 
    } 

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
  

    TF1 * templateFunc = new TF1("templateFunc", templateHistoFunc[0][0], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    // --- uncorrected    
    double xNorm = h_EoP_Data[i]->Integral()/h_template[0][0]->Integral() *
                   h_EoP_Data[i]->GetBinWidth(1)/h_template[0][0]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting uncorrected Data:  " << i <<" "<< std::endl; 
    TFitResultPtr rp;
    int fStatus; 
    for (int trial=0;trial<5;trial++) {
      templateFunc -> SetParameter(1, 0.99 + trial*0.01 );
      rp = h_EoP_Data[i] -> Fit("templateFunc", "MQRLS+");
      fStatus = rp;
      if (fStatus < 2) break; 
    } 
    g_EoP_Data -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_Data -> SetPointError(i, 0, def_error);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    //cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    // ---- corrected    
    templateFunc = new TF1("templateFunc", templateHistoFunc[0][0], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);
    xNorm = h_EoC_Data[i]->Integral()/h_template[0][0]->Integral() *
            h_EoC_Data[i]->GetBinWidth(1)/h_template[0][0]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting corrected Data " << i << "  "<< std::endl;
    for (int trial=0;trial<5;trial++) {
      templateFunc -> SetParameter(1, 0.99 + trial*0.01 );
      rp = h_EoC_Data[i] -> Fit("templateFunc", "MQRLS+");
      fStatus = rp;
      if (fStatus < 2) break; 
    } 

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
  
  for(unsigned int i = 0; i < nBins; ++i)
    {
      h_EoP_MC[i] -> Write(); 
      h_EoP_Data[i] -> Write(); 
    }

  htest->Write();
  htest2->Write();
  
}
