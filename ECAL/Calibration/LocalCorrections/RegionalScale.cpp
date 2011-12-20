// per compilare: g++ -Wall -o RegionalScale `root-config --cflags --glibs` RegionalScale.cpp

#include "histoFunc.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraphErrors.h"
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

#define xtalWidth 0.01745329

using namespace std;


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




int main(int argc, char** argv)
{
  //---- output file to save graphs
  char outfilename[100];
  //sprintf(outfilename,"GraphsRegionalScaleEta_highR9_vtxZEtaPos_diffTemplates.root");
  sprintf(outfilename,"test.root");
  
  //---- variables for selections
  float etaMax = 1.44;
  float r9min = 0.90 ;
  float r9max = 9999 ;  
  
  bool usePUweights = true;

  //--- weights for MC
  TFile weightsFile("weights/PUweights_2011_0100_73500_WJetsToLL_Fall11_S6.root","READ"); 
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();



  // --- NTUPLES
  TChain *ntu_MC = new TChain("ntu");
  TChain *ntu_Data = new TChain("ntu");
    
  //---- MC fall 2011
  ntu_MC->Add("../NTUPLES/Fall11/WZAnalysis/WZAnalysis_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1.root");

  //---- DATA
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
  
  std::cout << "     MC    : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;


  //---- Observables
  int npu;
  float PV_z;
  float EoP, scEta, scPhi;
  float scE3x3, scE5x5, scE, scE_regression;  
  float charge, scLocalEta, scLocalPhi,crackCorr; 
  float R9;
  float fbrem ;

  //---- Set branch addresses for MC  
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("PV_z", &PV_z);
  ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_MC->SetBranchAddress("ele1_scE", &scE);
  ntu_MC->SetBranchAddress("ele1_scE_regression", &scE_regression);
  ntu_MC->SetBranchAddress("ele1_charge", &charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_MC->SetBranchAddress("ele1_fbrem",&fbrem); 
  
  //---- Set branch addresses for Data
  ntu_Data->SetBranchAddress("PV_z", &PV_z);
  ntu_Data->SetBranchAddress("ele1_scEta", &scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_Data->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_Data->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_Data->SetBranchAddress("ele1_scE", &scE);
  ntu_Data->SetBranchAddress("ele1_scE_regression", &scE_regression);
  ntu_Data->SetBranchAddress("ele1_charge", &charge);
  ntu_Data->SetBranchAddress("ele1_fbrem",&fbrem); 
  
  unsigned int nBins = 85*2;
  std::cout << "nBins = " << nBins << std::endl;
  
  TProfile *pr = new TProfile("a","a",nBins,-1.*etaMax, etaMax, 0, 2);
  TProfile *fbremMC   = new TProfile("fbremMC","fbremMC",nBins,-1.*etaMax, etaMax, 0, 2);
  TProfile *fbremData = new TProfile("fbremData","fbremData",nBins,-1.*etaMax, etaMax, 0, 2);
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

  const int Ntempl = 4;
  //---- book templates
  // [0] --> uncorrected
  // [1] --> corrected
  TH1F* h_template_MC[Ntempl][2];
  TH1F* h_template_Data[Ntempl][2];
  for(unsigned int i = 0; i < Ntempl; ++i){
    char histoName[100];
    sprintf(histoName, "template_MC_%d", i);
    h_template_MC[i][0] = new TH1F(histoName, "", 1200, 0., 3.);   
    sprintf(histoName, "template_DATA_%d", i);
    h_template_Data[i][0] = new TH1F(histoName, "", 1200, 0., 3.);   
    sprintf(histoName, "template_MCcorr_%d", i);
    h_template_MC[i][1] = new TH1F(histoName, "", 1200, 0., 3.);   
    sprintf(histoName, "template_DATAcorr_%d", i);
    h_template_Data[i][1] = new TH1F(histoName, "", 1200, 0., 3.);   
  }

 
  //******************************************************************************************
  //*************************************** MC  ********************************************** 
   
  std::cout << "Loop in MC events " << endl; 
  float ww = 1;
  // loop on MC, make refernce and fit dist
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%100000 == 0 ) std::cout << "reading MC saved entry " << entry << std::endl;
    // if (entry>1000) break;

    ntu_MC->GetEntry(entry);
      
    //     if (PV_z*scEta<0) continue;

 
    if (usePUweights) ww = w[npu];

    R9 = scE3x3/scE;

    //-- eta or R9 cuts
    if( fabs(scEta) > etaMax ) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 
    
    //-- remove phi cracks 
    float phi = (scPhi+3.1415926536)/0.01745329;
    float modphi = (int)phi%20;
    if (fabs(modphi-10)<3.) continue;

    float fetaCry = fabs (scEta) / xtalWidth;

    //-- fill templates for each mod
    int mod = templIndex(fetaCry);
    float correction = scE_regression/scE;
    // reference out of the gaps
    if( IsEtaGap(fetaCry)==0 ) {
      h_template_MC[mod][0]-> Fill(EoP,ww);
      h_template_MC[mod][1]-> Fill(EoP*correction,ww);
    }

    // fill MC histos in eta bins
    int bin = pr->GetXaxis()->FindBin(scEta) - 1;
    h_EoP_MC[bin] -> Fill(EoP,ww);
    h_EoC_MC[bin] -> Fill(EoP*correction,ww);
    
    fbremMC->Fill(scEta,fbrem,ww);

  } 

  
  //******************************************************************************************
  //*************************************** DATA ********************************************** 
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
  //---- loop on Data
  for(int entry = 0; entry < ntu_Data->GetEntries(); ++entry) {
    if( entry%100000 == 0 ) std::cout << "reading data saved entry " << entry << std::endl;
    //    if (entry>1000) break;
    ntu_Data->GetEntry(entry);
    
    //  if (PV_z*scEta<0) continue;

    R9 = scE3x3/scE;

    //-- eta or R9 cuts
    if ( fabs(scEta) > etaMax) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 

    //-- remove phi cracks
    float phi = (scPhi+3.1415926536)/0.01745329;
    float modphi = (int)phi%20;
    if (fabs(modphi-10)<3.) continue;

    float fetaCry = fabs (scEta) / xtalWidth;

    //-- fill template for each mod
    int mod = templIndex(fetaCry);
    float correction = scE_regression/scE;
     // reference out of the gaps
    if( IsEtaGap(fetaCry)==0 ) {
      h_template_Data[mod][0]-> Fill(EoP);
      h_template_Data[mod][1]-> Fill(EoP*correction);
    }
    
    // fill Data histos in eta bins
    int bin = pr->GetXaxis()->FindBin(scEta) - 1;
    h_EoP_Data[bin] -> Fill(EoP);
    h_EoC_Data[bin] -> Fill(EoP*correction);

    fbremData->Fill(scEta,fbrem);
  } 
  
  //************************************* FITTING ***************************************************//
  
  TGraphErrors* g_EoP_MC   = new TGraphErrors();g_EoP_MC->SetName("gEoP_MC");
  TGraphErrors* g_EoC_MC   = new TGraphErrors();g_EoC_MC->SetName("gEoC_MC");
  TGraphErrors* g_Rat_MC   = new TGraphErrors();g_Rat_MC->SetName("gCorr_Uncorr_MC");
  
  histoFunc *templateHistoFuncMC[4][2]; 
  histoFunc *templateHistoFuncData[4][2]; 
  for (int mod=0; mod<4;mod++) {
    for (int ll = 0; ll < 2; ll++){
      h_template_MC[mod][ll] -> Rebin(rebin);
      templateHistoFuncMC[mod][ll] = new histoFunc(h_template_MC[mod][ll]);
      h_template_Data[mod][ll] -> Rebin(rebin);
      templateHistoFuncData[mod][ll] = new histoFunc(h_template_Data[mod][ll]);
    }
  }
  
  TF1 * templateFunc;
  double xNorm;

  for(unsigned int i = 0; i < nBins; ++i) {
    
    h_EoP_MC[i] -> Rebin(rebin);    
    h_EoC_MC[i] -> Rebin(rebin);    
    
    float xval = pr->GetXaxis()->GetBinCenter(i+1);
    int ieta = fabs(xval)/xtalWidth;
    int mod = templIndex(ieta);
    
    //-- MC uncorrected    
    std::cout << "***** Fitting uncorrected MC:  " << i << " " << mod; 
    templateFunc = new TF1("templateFunc", templateHistoFuncMC[mod][0], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);
    
    xNorm = h_EoP_MC[i]->Integral()/h_template_MC[mod][0]->Integral() *
            h_EoP_MC[i]->GetBinWidth(1)/h_template_MC[mod][0]->GetBinWidth(1); 
    
    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    int fStatus = h_EoP_MC[i] -> Fit("templateFunc", "MRQLS+");
    g_EoP_MC -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_MC -> SetPointError(i, 2*etaMax/nBins, 10);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    //cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //---ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    //-- MC corrected    
    std::cout << "***** Fitting corrected MC:  " << i << " " << mod; 
    templateFunc = new TF1("templateFunc", templateHistoFuncMC[mod][1], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    xNorm = h_EoC_MC[i]->Integral()/h_template_MC[mod][1]->Integral() *
            h_EoC_MC[i]->GetBinWidth(1)/h_template_MC[mod][1]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    fStatus = h_EoC_MC[i] -> Fit("templateFunc", "MRQLS+");
    g_EoC_MC -> SetPoint(i, xval , 1./templateFunc->GetParameter(1));
    g_EoC_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoC_MC -> SetPointError(i, 2*etaMax/nBins, 10);
    //cout << " ********** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //--ratio finalization
    rat /= templateFunc->GetParameter(1);
    era = rat*sqrt(era*era+templateFunc->GetParError(1)*templateFunc->GetParError(1)); 
    
    g_Rat_MC -> SetPoint(i,  xval , rat); 
    g_Rat_MC -> SetPointError(i,  0. , era); 
    g_Rat_MC->SetLineColor(kBlue+2); 

  }
 
    
  TGraphErrors* g_EoP_Data   = new TGraphErrors();g_EoP_Data->SetName("gEoP_Data");
  TGraphErrors* g_EoC_Data   = new TGraphErrors();g_EoC_Data->SetName("gEoC_Data");
  TGraphErrors* g_Rat_Data   = new TGraphErrors();g_Rat_Data->SetName("gCorr_Uncorr_Data");
  
  
  for(unsigned int i = 0; i < nBins; ++i)
  {
    h_EoP_Data[i] -> Rebin(rebin);    
    h_EoC_Data[i] -> Rebin(rebin);    
    
    
    float xval = pr->GetXaxis()->GetBinCenter(i+1);
    int ieta = fabs(xval)/xtalWidth;
    int mod = templIndex(ieta);


    
    //-- DATA uncorrected    
    std::cout << "***** Fitting uncorrected DATA:  " << i << " " << mod; 
    templateFunc = new TF1("templateFunc", templateHistoFuncData[mod][0], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);
    xNorm = h_EoP_Data[i]->Integral()/h_template_Data[mod][0]->Integral() *
            h_EoP_Data[i]->GetBinWidth(1)/h_template_Data[mod][0]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    int fStatus = h_EoP_Data[i] -> Fit("templateFunc", "MRQLS+");
    g_EoP_Data -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
    if(fStatus%10 == 4)  g_EoP_Data -> SetPointError(i, 2*etaMax/nBins, 10);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    //cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //--ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    //--- corrected DATA    
    std::cout << "***** Fitting corrected DATA:  " << i << " " << mod; 
    templateFunc = new TF1("templateFunc", templateHistoFuncData[mod][1], 0.7, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    xNorm = h_EoC_Data[i]->Integral()/h_template_Data[mod][1]->Integral() *
            h_EoC_Data[i]->GetBinWidth(1)/h_template_Data[mod][1]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
    
    std::cout << "***** Fitting corrected Data " << i << " " << mod;
    fStatus = h_EoC_Data[i] -> Fit("templateFunc", "SMRQL+");
    g_EoC_Data -> SetPoint(i, xval , 1./templateFunc->GetParameter(1));
    g_EoC_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
    //if(fStatus%10 == 4)  g_EoC_Data -> SetPointError(i, 2*etaMax/nBins, 10);
    //cout << " ********** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //--ratio finalization
    rat /= templateFunc->GetParameter(1);
    era = rat*sqrt(era*era+templateFunc->GetParError(1)*templateFunc->GetParError(1)); 
    
    g_Rat_Data -> SetPoint(i,  xval , rat); 
    g_Rat_Data -> SetPointError(i,  0. , era); 
    g_Rat_Data->SetLineColor(kBlue+2); 
    
  }
 
  //out file
  TFile outfile(outfilename,"recreate");

  g_EoP_MC->Write();
  g_EoC_MC->Write();
  g_Rat_MC->Write();

  g_EoP_Data->Write();
  g_EoC_Data->Write();
  g_Rat_Data->Write();

  for (unsigned int imod = 0; imod<4; imod++){
    h_template_Data[imod][0]->Write();  
    h_template_Data[imod][1]->Write();
    h_template_MC[imod][1]->Write();  
    h_template_MC[imod][1]->Write();  
  }

  for (unsigned int i=0; i<nBins;i++ )
    {
      h_EoP_MC[i]->Write() ;
      h_EoC_MC[i]->Write() ;
      h_EoP_Data[i]->Write() ;
      h_EoC_Data[i]->Write() ;
 

    }

  fbremMC->Write();
  fbremData->Write();

}
