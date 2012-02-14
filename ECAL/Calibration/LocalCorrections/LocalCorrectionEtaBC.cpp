// per compilare: g++ -Wall -o LocalCorrectionEtaBC `root-config --cflags --glibs` LocalCorrectionEtaBC.cpp

#include "../CommonTools/histoFunc.h"
#include "../CommonTools/TPileupReweighting.h"
#include "../CommonTools/BasicClusterLocalContainementCorrection.C"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
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
#define PI        3.1415926536 

using std::vector;

//derive f(local eta) correction in MC, with 4 different f for each module.
//fold crystal in eta into one single module

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
  gROOT->ProcessLine("#include <vector>");

  //---- variables for selections
  float etaMax = 1.44;
  
  float r9min  = 0.0 ;
  float r9max  = 9999 ;  
  float bcNmin = 1;
  float bcNmax = 2;
  
  bool useW         = false;
  bool useZ         = true;
  bool usePUweights = true;
  
  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"GraphsLocalEta_singleBCcorrectionPol2_Zee_nBC1.root");

  //---- PU weights for MC
  TPileupReweighting *puReweighting;
  if (useW) puReweighting = new TPileupReweighting("../CommonTools/weights/PUweights_2011_0100_73500_WJetsToLL_Fall11_S6.root","hweights"); 
  if (useZ) puReweighting = new TPileupReweighting("../CommonTools/weights/PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","hweights"); 
    
  //---- NTUPLES
  TChain *ntu_MC = new TChain("ntu");
  TChain *ntu_Data = new TChain("ntu");
    
  //---- MC fall 2011  //---- MC fall 2011
  if (useW)
    ntu_MC->Add("../NTUPLES/Fall11/WZAnalysis/WZAnalysis_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1.root");
  if (useZ)
    ntu_MC->Add("../NTUPLES/Fall11/WZAnalysis/WZAnalysis_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1.root");
 
  //---- DATA
  if (useW){
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
  }
  if (useZ){
    // ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
    // ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
    // ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
    // ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
    // ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B_bc.root");
  }
 
  std::cout << "     MC    : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;

  //---- Observables
  int isZ;
  int npu;
  float PV_z;
  float ele1_EoP, ele1_scEta, ele1_scPhi;
  float ele1_scE3x3, ele1_scE, ele1_scE_regression;  
  float ele1_charge, ele1_scLocalEta, ele1_scLocalPhi,ele1_scLocalCorr; 
  float ele1_scCrackCorr, ele1_fCorr;
  float ele2_EoP, ele2_scEta, ele2_scPhi;
  float ele2_scE3x3, ele2_scE, ele2_scE_regression;  
  float ele2_charge, ele2_scLocalEta, ele2_scLocalPhi,ele2_scLocalCorr; 
  float ele2_scCrackCorr, ele2_fCorr;
  float R9;
  int ele1_bcN;
  int ele2_bcN;
  vector<float>   *ele1_bcE = 0;
  vector<float>   *ele1_bcEta= 0;
  vector<float>   *ele1_bcPhi= 0;
  vector<float>   *ele1_bcLocalEta= 0;  
  vector<float>   *ele1_bcLocalPhi= 0;
  vector<float>   *ele2_bcE = 0;
  vector<float>   *ele2_bcEta= 0;
  vector<float>   *ele2_bcPhi= 0;
  vector<float>   *ele2_bcLocalEta= 0;  
  vector<float>   *ele2_bcLocalPhi= 0;

  
  //---- Set branch addresses for MC  
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("isZ", &isZ);
  ntu_MC->SetBranchAddress("PV_z", &PV_z);
  ntu_MC->SetBranchAddress("ele1_scEta", &ele1_scEta);
  ntu_MC->SetBranchAddress("ele1_scPhi", &ele1_scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &ele1_EoP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &ele1_scE3x3);
  ntu_MC->SetBranchAddress("ele1_scE", &ele1_scE);
  ntu_MC->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression);
  ntu_MC->SetBranchAddress("ele1_charge", &ele1_charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&ele1_scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&ele1_scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_fCorrection",    &ele1_fCorr);
  ntu_MC->SetBranchAddress("ele1_scCrackCorr",    &ele1_scCrackCorr);
  ntu_MC->SetBranchAddress("ele1_scLocalContCorr",&ele1_scLocalCorr); 
  ntu_MC->SetBranchAddress("ele1_bcN", &ele1_bcN);
  ntu_MC->SetBranchAddress("ele1_bcE",            &ele1_bcE);
  ntu_MC->SetBranchAddress("ele1_bcEta",          &ele1_bcEta);
  ntu_MC->SetBranchAddress("ele1_bcPhi",          &ele1_bcPhi);
  ntu_MC->SetBranchAddress("ele1_bcLocalEta",     &ele1_bcLocalEta);
  ntu_MC->SetBranchAddress("ele1_bcLocalPhi",     &ele1_bcLocalPhi);
  ntu_MC->SetBranchAddress("ele2_scEta", &ele2_scEta);
  ntu_MC->SetBranchAddress("ele2_scPhi", &ele2_scPhi);
  ntu_MC->SetBranchAddress("ele2_EOverP", &ele2_EoP);
  ntu_MC->SetBranchAddress("ele2_e3x3", &ele2_scE3x3);
  ntu_MC->SetBranchAddress("ele2_scE", &ele2_scE);
  ntu_MC->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression);
  ntu_MC->SetBranchAddress("ele2_charge", &ele2_charge);
  ntu_MC->SetBranchAddress("ele2_scLocalPhi",&ele2_scLocalPhi); 
  ntu_MC->SetBranchAddress("ele2_scLocalEta",&ele2_scLocalEta); 
  ntu_MC->SetBranchAddress("ele2_fCorrection",    &ele2_fCorr);
  ntu_MC->SetBranchAddress("ele2_scCrackCorr",    &ele2_scCrackCorr);
  ntu_MC->SetBranchAddress("ele2_scLocalContCorr",&ele2_scLocalCorr); 
  ntu_MC->SetBranchAddress("ele2_bcN", &ele2_bcN);
  ntu_MC->SetBranchAddress("ele2_bcE",            &ele2_bcE);
  ntu_MC->SetBranchAddress("ele2_bcEta",          &ele2_bcEta);
  ntu_MC->SetBranchAddress("ele2_bcPhi",          &ele2_bcPhi);
  ntu_MC->SetBranchAddress("ele2_bcLocalEta",     &ele2_bcLocalEta);
  ntu_MC->SetBranchAddress("ele2_bcLocalPhi",     &ele2_bcLocalPhi);

  //---- Set branch addresses for Data
  ntu_Data->SetBranchAddress("isZ", &isZ);
  ntu_Data->SetBranchAddress("PV_z", &PV_z);
  ntu_Data->SetBranchAddress("ele1_scEta", &ele1_scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi", &ele1_scPhi);
  ntu_Data->SetBranchAddress("ele1_EOverP", &ele1_EoP);
  ntu_Data->SetBranchAddress("ele1_e3x3", &ele1_scE3x3);
  ntu_Data->SetBranchAddress("ele1_scE", &ele1_scE);
  ntu_Data->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression);
  ntu_Data->SetBranchAddress("ele1_charge", &ele1_charge);
  ntu_Data->SetBranchAddress("ele1_scLocalPhi",&ele1_scLocalPhi); 
  ntu_Data->SetBranchAddress("ele1_scLocalEta",&ele1_scLocalEta); 
  ntu_Data->SetBranchAddress("ele1_fCorrection",    &ele1_fCorr);
  ntu_Data->SetBranchAddress("ele1_scCrackCorr",    &ele1_scCrackCorr);
  ntu_Data->SetBranchAddress("ele1_scLocalContCorr",&ele1_scLocalCorr); 
  ntu_Data->SetBranchAddress("ele1_bcN", &ele1_bcN);
  ntu_Data->SetBranchAddress("ele1_bcE",            &ele1_bcE);
  ntu_Data->SetBranchAddress("ele1_bcEta",          &ele1_bcEta);
  ntu_Data->SetBranchAddress("ele1_bcPhi",          &ele1_bcPhi);
  ntu_Data->SetBranchAddress("ele1_bcLocalEta",     &ele1_bcLocalEta);
  ntu_Data->SetBranchAddress("ele1_bcLocalPhi",     &ele1_bcLocalPhi);
  ntu_Data->SetBranchAddress("ele2_scEta", &ele2_scEta);
  ntu_Data->SetBranchAddress("ele2_scPhi", &ele2_scPhi);
  ntu_Data->SetBranchAddress("ele2_EOverP", &ele2_EoP);
  ntu_Data->SetBranchAddress("ele2_e3x3", &ele2_scE3x3);
  ntu_Data->SetBranchAddress("ele2_scE", &ele2_scE);
  ntu_Data->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression);
  ntu_Data->SetBranchAddress("ele2_charge", &ele2_charge);
  ntu_Data->SetBranchAddress("ele2_scLocalPhi",&ele2_scLocalPhi); 
  ntu_Data->SetBranchAddress("ele2_scLocalEta",&ele2_scLocalEta); 
  ntu_Data->SetBranchAddress("ele2_fCorrection",    &ele2_fCorr);
  ntu_Data->SetBranchAddress("ele2_scCrackCorr",    &ele2_scCrackCorr);
  ntu_Data->SetBranchAddress("ele2_scLocalContCorr",&ele2_scLocalCorr); 
  ntu_Data->SetBranchAddress("ele2_bcN", &ele2_bcN);
  ntu_Data->SetBranchAddress("ele2_bcE",            &ele2_bcE);
  ntu_Data->SetBranchAddress("ele2_bcEta",          &ele2_bcEta);
  ntu_Data->SetBranchAddress("ele2_bcPhi",          &ele2_bcPhi);
  ntu_Data->SetBranchAddress("ele2_bcLocalEta",     &ele2_bcLocalEta);
  ntu_Data->SetBranchAddress("ele2_bcLocalPhi",     &ele2_bcLocalPhi);

  const unsigned int nBins = 20;
  const int Ntempl = 4;
  std::cout << "nBins = " << nBins << std::endl;
  
  // histogram definition
  TH1F* h_EoP_MC[nBins][Ntempl] ;   
  TH1F* h_EoC_MC[nBins][Ntempl] ;
  TH1F* h_EoP_Data[nBins][Ntempl];
  TH1F* h_EoC_Data[nBins][Ntempl] ;
  
  for(int mod=0; mod<Ntempl; mod++){
    for(unsigned int i = 0; i < nBins; ++i)
      {
	char histoName[80];
	sprintf(histoName, "EoP_MC_%d_mod%d", i,mod+1);
	h_EoP_MC[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoP_MC[i][mod] -> SetFillColor(4);
	h_EoP_MC[i][mod] -> SetFillStyle(3004);
	sprintf(histoName, "EoC_MC_%d_mod%d", i,mod+1);
	h_EoC_MC[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoC_MC[i][mod] -> SetFillColor(3);
	h_EoC_MC[i][mod] -> SetFillStyle(3004);
	sprintf(histoName, "EoP_Data_%d_mod%d", i,mod+1);
	h_EoP_Data[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoP_Data[i][mod] -> SetFillColor(4);
	h_EoP_Data[i][mod] -> SetFillStyle(3004);
	sprintf(histoName, "EoC_Data_%d_mod%d", i,mod+1);
	h_EoC_Data[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoC_Data[i][mod] -> SetFillColor(3);
	h_EoC_Data[i][mod] -> SetFillStyle(3004);
      }
  }

  TH1F* h_spreadEoP_MC[Ntempl] ;   
  TH1F* h_spreadEoC_MC[Ntempl] ;
  TH1F* h_spreadEoP_Data[Ntempl] ;   
  TH1F* h_spreadEoC_Data[Ntempl] ;
  TH1F* h_spreadRatio[Ntempl] ;
  TH1F* h_spreadRatioCorr[Ntempl] ;

  for(int mod=0; mod<Ntempl; mod++){
    char histoName[80];
    sprintf(histoName, "spreadEoP_MC_mod%d",mod+1);
    h_spreadEoP_MC[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadEoP_MC[mod] -> SetFillColor(4);
    h_spreadEoP_MC[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadEoC_MC_mod%d",mod+1);
    h_spreadEoC_MC[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadEoC_MC[mod] -> SetFillColor(3);
    h_spreadEoC_MC[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadEoP_Data_mod%d",mod+1);
    h_spreadEoP_Data[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadEoP_Data[mod] -> SetFillColor(4);
    h_spreadEoP_Data[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadEoC_Data_mod%d", mod+1);
    h_spreadEoC_Data[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadEoC_Data[mod] -> SetFillColor(3);
    h_spreadEoC_Data[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadRatio_mod%d", mod+1);
    h_spreadRatio[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadRatio[mod] -> SetFillColor(3);
    h_spreadRatio[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadRatioCorr_mod%d", mod+1);
    h_spreadRatioCorr[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadRatioCorr[mod] -> SetFillColor(3);
    h_spreadRatioCorr[mod] -> SetFillStyle(3004);

  }


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
  std::cout << "Loop on MC events ... " << std::endl; 
 
  float ww = 1 ;
  //---- loop on MC, make refernce and fit dist
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    
    if( entry%500000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //if (entry>500000) break;

    ntu_MC->GetEntry(entry);

    if (isZ==0) continue;

    R9 = ele1_scE3x3/ele1_scE;
    
    // -- PU weights
    if (usePUweights) ww = puReweighting->GetWeight(npu);
   
    //-- eta or R9 cuts
    if( fabs(ele1_scEta) > etaMax ) continue;
    if( fabs(ele2_scEta) > etaMax ) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 

    // -- first electron
    
    //-- remove phi cracks
    float phi = (ele1_scPhi+3.1415926536)/xtalWidth;
    float modphi = (int)phi%20;
 
    //-- remove eta gaps
    float fetaCry = fabs (ele1_scEta) / xtalWidth;


    if ( fabs(modphi-10)>2.  && IsEtaGap(fetaCry)==0 && ele1_bcN >= bcNmin && ele1_bcN < bcNmax){

      //-- fill template for each mod
      int mod = templIndex(fetaCry);
      // // corrected energy: use single BC correction
      // float scEc = 0;
      // for (int ibc = 0; ibc < ele1_bcN; ibc++){
      // 	scEc+= ele1_bcE->at(ibc)*BasicClusterLocalContainementCorrection(ele1_bcLocalEta->at(ibc),mod);
      // }
      // scEc*=ele1_fCorr*ele1_scCrackCorr;
      
      float correction = ele1_scE_regression/ele1_scE;
      //float correction = scEc/ele1_scE;
      h_template_MC[mod][0]-> Fill(ele1_EoP,ww);
      h_template_MC[mod][1]-> Fill(ele1_EoP*correction,ww);
      
      //-- fill MC histos in eta bins
      float locEta = ele1_scLocalEta+0.5; 
      if ( fabs(locEta-0.5) >= 0.5 ) continue;
      
      int bin = nBins * (locEta); 
      if (bin>nBins-1 || bin < 0 ) {
	std::cout << "Error in bins: " << bin << " " << ele1_scLocalEta << std::endl;
	continue;
      }
      
      h_EoP_MC[bin][mod] -> Fill(ele1_EoP,ww);
      h_EoC_MC[bin][mod] -> Fill(ele1_EoP*correction,ww);
    }


    // --- second electron
    
    //-- remove phi cracks
    phi = (ele2_scPhi+3.1415926536)/xtalWidth;
    modphi = (int)phi%20;
 
    //-- remove eta gaps
    fetaCry = fabs (ele2_scEta) / xtalWidth;

    if ( fabs(modphi-10)>2.  && IsEtaGap(fetaCry)==0 && ele2_bcN >= bcNmin && ele2_bcN < bcNmax){

      //-- fill template for each mod
      int mod = templIndex(fetaCry);
      // // corrected energy: use single BC correction
      // float scEc = 0;
      // for (int ibc = 0; ibc < ele2_bcN; ibc++){
      // 	scEc+= ele2_bcE->at(ibc)*BasicClusterLocalContainementCorrection(ele2_bcLocalEta->at(ibc),mod);
      // }
      // scEc*=ele2_fCorr*ele2_scCrackCorr;
      
      float correction = ele2_scE_regression/ele2_scE;
      //float correction = scEc/ele2_scE;
      h_template_MC[mod][0]-> Fill(ele2_EoP,ww);
      h_template_MC[mod][1]-> Fill(ele2_EoP*correction,ww);
      
      //-- fill MC histos in eta bins
      float locEta = ele2_scLocalEta+0.5; 
      if ( fabs(locEta-0.5) >= 0.5 ) continue;
      
      int bin = nBins * (locEta); 
      if (bin>nBins-1 || bin < 0 ) {
	std::cout << "Error in bins: " << bin << " " << ele2_scLocalEta << std::endl;
	continue;
      }
      
      h_EoP_MC[bin][mod] -> Fill(ele2_EoP,ww);
      h_EoC_MC[bin][mod] -> Fill(ele2_EoP*correction,ww);
    }

  }

  std::cout << "stat: "<< h_EoP_MC[0][10]->GetEntries()<< " " <<  h_EoC_MC[0][10]->GetEntries() << std::endl; 

  
  //******************************************************************************************
  //*************************************** DATA ********************************************** 
  std::cout << "Loop on Data events ..." << std::endl; 
  //--- loop on data
  for(int entry = 0; entry < ntu_Data->GetEntries(); ++entry) {
    if( entry%500000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //if (entry>500000) break;

    ntu_Data->GetEntry(entry);

    R9 = ele1_scE3x3/ele1_scE;

    //-- eta or R9 cuts
    if( fabs(ele1_scEta) > etaMax ) continue;
    if( fabs(ele2_scEta) > etaMax ) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 


    // -- first electron
    
    //-- remove phi cracks
    float phi = (ele1_scPhi+3.1415926536)/xtalWidth;
    float modphi = (int)phi%20;
 
    //-- remove eta gaps
    float fetaCry = fabs (ele1_scEta) / xtalWidth;


    if ( fabs(modphi-10)>2.  && IsEtaGap(fetaCry)==0 && ele1_bcN >= bcNmin && ele1_bcN < bcNmax){

      //-- fill template for each mod
      int mod = templIndex(fetaCry);
      // corrected energy: use single BC correction
      float scEc = 0;
      for (int ibc = 0; ibc < ele1_bcN; ibc++){
	scEc+= ele1_bcE->at(ibc)*BasicClusterLocalContainementCorrection(ele1_bcLocalEta->at(ibc),mod);
      }
      scEc*=ele1_fCorr*ele1_scCrackCorr;
      
      //float correction = scE_regression/scE;
      float correction = scEc/ele1_scE;
      h_template_Data[mod][0]-> Fill(ele1_EoP);
      h_template_Data[mod][1]-> Fill(ele1_EoP*correction);
      
      //-- fill MC histos in eta bins
      float locEta = ele1_scLocalEta+0.5; 
      if ( fabs(locEta-0.5) >= 0.5 ) continue;
      
      int bin = nBins * (locEta); 
      if (bin>nBins-1 || bin < 0 ) {
	std::cout << "Error in bins: " << bin << " " << ele1_scLocalEta << std::endl;
	continue;
      }
      
      h_EoP_Data[bin][mod] -> Fill(ele1_EoP);
      h_EoC_Data[bin][mod] -> Fill(ele1_EoP*correction);
    }


    // --- second electron
    
    //-- remove phi cracks
    phi = (ele2_scPhi+3.1415926536)/xtalWidth;
    modphi = (int)phi%20;
 
    //-- remove eta gaps
    fetaCry = fabs (ele2_scEta) / xtalWidth;

    if ( fabs(modphi-10)>2.  && IsEtaGap(fetaCry)==0 && ele2_bcN >= bcNmin && ele2_bcN < bcNmax){

      //-- fill template for each mod
      int mod = templIndex(fetaCry);
      // corrected energy: use single BC correction
      float scEc = 0;
      for (int ibc = 0; ibc < ele2_bcN; ibc++){
	scEc+= ele2_bcE->at(ibc)*BasicClusterLocalContainementCorrection(ele2_bcLocalEta->at(ibc),mod);
      }
      scEc*=ele2_fCorr*ele2_scCrackCorr;
      
      //float correction = scE_regression/scE;
      float correction = scEc/ele2_scE;
      h_template_Data[mod][0]-> Fill(ele2_EoP);
      h_template_Data[mod][1]-> Fill(ele2_EoP*correction);
      
      //-- fill MC histos in eta bins
      float locEta = ele2_scLocalEta+0.5; 
      if ( fabs(locEta-0.5) >= 0.5 ) continue;
      
      int bin = nBins * (locEta); 
      if (bin>nBins-1 || bin < 0 ) {
	std::cout << "Error in bins: " << bin << " " << ele2_scLocalEta << std::endl;
	continue;
      }
      
      h_EoP_Data[bin][mod] -> Fill(ele2_EoP);
      h_EoC_Data[bin][mod] -> Fill(ele2_EoP*correction);
    }



  }
  
  
  ///////////////****************** Fit the histograms and fill the graphs *************** ////////////////////////
  int rebin = 8;

  TGraphErrors* g_EoP_MC[Ntempl];
  TGraphErrors* g_EoC_MC[Ntempl];
  TGraphErrors* g_ratio_MC[Ntempl];
  
  TGraphErrors* g_EoP_Data[Ntempl];
  TGraphErrors* g_EoC_Data[Ntempl];
  TGraphErrors* g_ratio_Data[Ntempl];

  TGraphErrors* g_ratio_uncorr[Ntempl];
  TGraphErrors* g_ratio_corr[Ntempl];
 
  for (int mod=0; mod<4; mod++){
    char histoName[100];
    
    sprintf(histoName, "gEoP_MC_mod%d", mod+1);
    g_EoP_MC[mod]   = new TGraphErrors(); 
    g_EoP_MC[mod]->SetName(histoName);
    
    sprintf(histoName, "gEoC_MC_mod%d", mod+1);
    g_EoC_MC[mod]   = new TGraphErrors(); 
    g_EoC_MC[mod]->SetName(histoName);
    
    sprintf(histoName, "gRatio_MC_mod%d", mod+1);
    g_ratio_MC[mod]   = new TGraphErrors(); 
    g_ratio_MC[mod]->SetName(histoName);
    
    sprintf(histoName, "gEoP_Data_mod%d", mod+1);
    g_EoP_Data[mod]   = new TGraphErrors(); 
    g_EoP_Data[mod]->SetName(histoName);
    
    sprintf(histoName, "gEoC_Data_mod%d", mod+1);
    g_EoC_Data[mod]   = new TGraphErrors(); 
    g_EoC_Data[mod]->SetName(histoName);
    
    sprintf(histoName, "gRatio_Data_mod%d", mod+1);
    g_ratio_Data[mod]   = new TGraphErrors(); 
    g_ratio_Data[mod]->SetName(histoName);

    sprintf(histoName, "gRatio_uncorr_mod%d", mod+1);
    g_ratio_uncorr[mod]   = new TGraphErrors(); 
    g_ratio_uncorr[mod]->SetName(histoName);

    sprintf(histoName, "gRatio_corr_mod%d", mod+1);
    g_ratio_corr[mod]   = new TGraphErrors(); 
    g_ratio_corr[mod]->SetName(histoName);
 
  }
  
  

  //************************************* FITTING ***************************************************//
    
  for(int mod=0;mod<4;mod++){

    //--- define func template
    histoFunc *templateHistoFuncMC[2]; 
    histoFunc *templateHistoFuncData[2]; 
    for (int ll = 0; ll < 2; ll++){
      h_template_MC[mod][ll]   -> Rebin(rebin);
      h_template_Data[mod][ll] -> Rebin(rebin);
      templateHistoFuncMC[ll]   = new histoFunc(h_template_MC[mod][ll]);
      templateHistoFuncData[ll] = new histoFunc(h_template_Data[mod][ll]);
    } 

    for(unsigned int i = 0; i < nBins; ++i)
      {
	std::cout << "***** Fitting :  mod: " <<mod <<"  bin:  "<< i << std::endl; 
	
	h_EoP_MC[i][mod] -> Rebin(rebin);    
	h_EoC_MC[i][mod] -> Rebin(rebin);    
	h_EoP_Data[i][mod] -> Rebin(rebin);    
	h_EoC_Data[i][mod] -> Rebin(rebin);    


	//************************ MC ****************************************************************

	TF1 * templateFuncMC = new TF1("templateFuncMC", templateHistoFuncMC[0], 0.7, 1.3, 3, "histoFunc");
	templateFuncMC -> SetNpx(10000);

	float xval = (i+0.5)*1/(float)nBins - 0.5;

	//--- uncorrected MC   
	double xNorm = h_EoP_MC[i][mod]->Integral()/h_template_MC[mod][0]->Integral() *
	               h_EoP_MC[i][mod]->GetBinWidth(1)/h_template_MC[mod][0]->GetBinWidth(1); 

	templateFuncMC -> FixParameter(0, xNorm);
	templateFuncMC -> SetParameter(1, 1.05 );
	templateFuncMC -> FixParameter(2, 0.);
    

	TFitResultPtr	r = h_EoP_MC[i][mod] -> Fit("templateFuncMC", "SMRQL+");
	int fitStatus = r;
	
	g_EoP_MC[mod] -> SetPoint(i, xval , 1./templateFuncMC->GetParameter(1));
	g_EoP_MC[mod] -> SetPointError(i, 0., templateFuncMC->GetParError(1));
	h_spreadEoP_MC[mod] -> Fill(1./templateFuncMC->GetParameter(1));
    
//    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
	//    cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 
	
	float scaleMC  = 1./templateFuncMC->GetParameter(1);
	float escaleMC = templateFuncMC->GetParError(1)/scaleMC/scaleMC;


	//--- corrected MC  
	templateFuncMC = new TF1("templateFuncMC", templateHistoFuncMC[1], 0.7, 1.3, 3, "histoFunc");
	templateFuncMC -> SetNpx(10000);
	xNorm = h_EoC_MC[i][mod]->Integral()/h_template_MC[mod][1]->Integral() *
	        h_EoC_MC[i][mod]->GetBinWidth(1)/h_template_MC[mod][1]->GetBinWidth(1); 

	templateFuncMC -> FixParameter(0, xNorm);
	templateFuncMC -> SetParameter(1, 1.05 );
	templateFuncMC -> FixParameter(2, 0.);
    

	h_EoC_MC[i][mod] -> Fit("templateFuncMC", "MRQLS+");
	g_EoC_MC[mod] -> SetPoint(i, xval , 1./templateFuncMC->GetParameter(1));
	g_EoC_MC[mod] -> SetPointError(i, 0., templateFuncMC->GetParError(1));
	h_spreadEoC_MC[mod] -> Fill(1./templateFuncMC->GetParameter(1));
    

	float scaleMCcorr  = 1./templateFuncMC->GetParameter(1);
	float escaleMCcorr = templateFuncMC->GetParError(1)/scaleMCcorr/scaleMCcorr;




	//************************ DATA **********************************************
	TF1 * templateFuncData = new TF1("templateFuncData", templateHistoFuncData[0], 0.7, 1.3, 3, "histoFunc");
	templateFuncData -> SetNpx(10000);


	//--- uncorrected Data   
	xNorm = h_EoP_Data[i][mod]->Integral()/h_template_Data[mod][0]->Integral() *
	        h_EoP_Data[i][mod]->GetBinWidth(1)/h_template_Data[mod][0]->GetBinWidth(1); 

	templateFuncData -> FixParameter(0, xNorm);
	templateFuncData -> SetParameter(1, 1.05 );
	templateFuncData -> FixParameter(2, 0.);
    

	h_EoP_Data[i][mod] -> Fit("templateFuncData", "MRQLS+");
	g_EoP_Data[mod] -> SetPoint(i, xval , 1./templateFuncData->GetParameter(1));
	g_EoP_Data[mod] -> SetPointError(i, 0., templateFuncData->GetParError(1));
	h_spreadEoP_Data[mod] -> Fill(1./templateFuncData->GetParameter(1));
//    if ( templateFuncData->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
	//    cout  << " ***** " <<  1./templateFuncData->GetParameter(1) << " " << templateFuncData->GetParError(1) << endl; 

	float scaleDA  = 1./templateFuncData->GetParameter(1);
	float escaleDA = templateFuncData->GetParError(1)/scaleDA/scaleDA;

	//--- corrected Data 
	templateFuncData = new TF1("templateFuncData", templateHistoFuncData[1], 0.7, 1.3, 3, "histoFunc");
	templateFuncData -> SetNpx(10000);
	xNorm = h_EoC_Data[i][mod]->Integral()/h_template_Data[mod][1]->Integral() *
	        h_EoC_Data[i][mod]->GetBinWidth(1)/h_template_Data[mod][1]->GetBinWidth(1); 

	templateFuncData -> FixParameter(0, xNorm);
	templateFuncData -> SetParameter(1, 1.05 );
	templateFuncData -> FixParameter(2, 0.);
    
	h_EoC_Data[i][mod] -> Fit("templateFuncData", "MRQLS+");
	g_EoC_Data[mod] -> SetPoint(i, xval, 1./templateFuncData->GetParameter(1));
	g_EoC_Data[mod] -> SetPointError(i, 0., templateFuncData->GetParError(1));
	h_spreadEoC_Data[mod] -> Fill(1./templateFuncData->GetParameter(1));

	float scaleDAcorr  = 1./templateFuncData->GetParameter(1);
	float escaleDAcorr = templateFuncData->GetParError(1)/scaleDAcorr/scaleDAcorr;

	//--- ratio finalization MC corr/uncorr
	float ratioMC = scaleMCcorr/scaleMC;
	float eratioMC = ratioMC*sqrt(pow(escaleMC/scaleMC,2) + pow(escaleMCcorr/scaleMCcorr,2)); 
	g_ratio_MC[mod] -> SetPoint(i, xval, ratioMC);
	g_ratio_MC[mod] -> SetPointError(i, 0., eratioMC);

	//--- ratio finalization DATA corr/uncorr
	float ratioDA  = scaleDAcorr/scaleDA;
	float eratioDA = ratioDA*sqrt(pow(escaleDA/scaleDA,2) + pow(escaleDAcorr/scaleDAcorr,2)); 
	g_ratio_Data[mod] -> SetPoint(i, xval , ratioDA);
	g_ratio_Data[mod] -> SetPointError(i, 0., eratioDA);


	//--- ratio finalization data/MC uncorrected
	float ratioU = scaleDA/scaleMC;
	float eratioU = ratioU*sqrt(pow(escaleMC/scaleMC,2) + pow(escaleDA/scaleDA,2)); 
	g_ratio_uncorr[mod] -> SetPoint(i, xval , ratioU);
	g_ratio_uncorr[mod] -> SetPointError(i, 0., eratioU);
	h_spreadRatio[mod] -> Fill(ratioU);
    

	//--- ratio finalization data/MC corrected
	float ratioC = scaleDAcorr/scaleMCcorr;
	float eratioC = ratioC*sqrt(pow(escaleMCcorr/scaleMCcorr,2) + pow(escaleDAcorr/scaleDAcorr,2)); 
	g_ratio_corr[mod] -> SetPoint(i, xval , ratioC);
	g_ratio_corr[mod] -> SetPointError(i, 0., eratioC);
    	h_spreadRatioCorr[mod] -> Fill(ratioC);
      }

  }

  
  TFile fout(outfilename,"recreate");
  for(int mod=0;mod<4;mod++) {
    g_EoP_MC[mod]->Write();
    g_EoC_MC[mod]->Write();
    g_EoP_Data[mod]->Write();
    g_EoC_Data[mod]->Write();
    g_ratio_MC[mod]->Write();
    g_ratio_Data[mod]->Write();
    g_ratio_uncorr[mod]->Write();
    g_ratio_corr[mod]->Write();

    h_spreadEoP_MC[mod]->Write();
    h_spreadEoC_MC[mod]->Write();
    h_spreadEoP_Data[mod]->Write();
    h_spreadEoC_Data[mod]->Write();

    h_spreadRatio[mod]->Write();
    h_spreadRatioCorr[mod]->Write();

    for (int mod = 0; mod < 4; mod++){
      h_template_MC[mod][0]->Write();
      h_template_MC[mod][1]->Write();
      h_template_Data[mod][0]->Write();
      h_template_Data[mod][1]->Write();
    }
 
  }
  fout.Close();
}
