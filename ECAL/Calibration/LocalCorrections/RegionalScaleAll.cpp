// per compilare: g++ -Wall -o RegionalScaleAll `root-config --cflags --glibs` RegionalScaleAll.cpp

#include "histoFunc.h"
#include "TEndcapRings.h"
#include "TMomentumCalibration.h"

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
#define PI        3.1415926536 

using namespace std;


// Template on MC: 1 per module
// Test on MC and data 
bool IsEtaGap(float eta){
  float feta = fabs(eta);
  if( fabs(feta - 0 )<4) return true;
  if( fabs(feta - 25)<4) return true;
  if( fabs(feta - 45)<4) return true;
  if( fabs(feta - 65)<4) return true;
  if( fabs(feta - 85)<4) return true;
  if( fabs(feta) > 121 ) return true;
  return false;
}

int templIndex(float ieta){
    float feta = fabs(ieta);
    if (feta <= 25)               {return 0;}
    if (feta>  25 && feta <=  45) {return 1;}
    if (feta>  45 && feta <=  65) {return 2;}
    if (feta>  65 && feta <=  85) {return 3;}  
    if (feta>  85 && feta <=  98) {return 4;}
    if (feta>  98 && feta <= 108) {return 5;}
    if (feta> 108 && feta <= 118) {return 6;}
    if (feta> 118 )               {return 7;}

    return -1;
}



int main(int argc, char** argv)
{
  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"GraphsRegionalScaleEta_allR9_Pcalib_noES.root");
    
  //---- variables for selections
  float etaMax = 2.5;
  float r9min = 0.00 ;
  float r9max = 9999 ;  
  
  bool useZ              = false;
  bool useW              = true;
  bool applyPcalibration = true;
  bool usePUweights      = true;

  //--- weights for MC
  std::cout << "Reading PU weights from : weights/PUweights_2011_0100_73500_WJetsToLL_Fall11_S6.root" << std::endl;
  TFile weightsFile("weights/PUweights_2011_0100_73500_WJetsToLL_Fall11_S6.root","READ"); 
  //TFile weightsFile("weights/PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","READ"); 
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();

  //--- P calibration
  TMomentumCalibration *pCalib = new TMomentumCalibration(); 
  // for (int i = -125 ; i < 125 ; i++){
  //   cout << " ieta = " << i << "   p scale = " << pCalib->GetMomentumCalibration(i,0) << endl;;
  // }
  
  // --- NTUPLES
  TChain *ntu_MC = new TChain("ntu");
  TChain *ntu_Data = new TChain("ntu");
    
  //---- MC fall 2011
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
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
  }
  std::cout << "     MC    : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;


  //---- Observables
  int npu;
  float PV_z;
  float EoP, scEta, scPhi;
  float scE3x3, scE5x5, scE, scE_regression, esE, tkP;  
  float charge, scLocalEta, scLocalPhi,crackCorr; 
  float R9;
  float fbrem ;
  int iphi,ieta,ix,iy,iz; 

  //---- Set branch addresses for MC  
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("PV_z", &PV_z);
  ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_MC->SetBranchAddress("ele1_scE", &scE);
  ntu_MC->SetBranchAddress("ele1_es", &esE);
  ntu_MC->SetBranchAddress("ele1_scE_regression", &scE_regression);
  ntu_MC->SetBranchAddress("ele1_charge", &charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_MC->SetBranchAddress("ele1_fbrem",&fbrem); 
  ntu_MC->SetBranchAddress("ele1_seedIphi",&iphi); 
  ntu_MC->SetBranchAddress("ele1_seedIeta",&ieta); 
  ntu_MC->SetBranchAddress("ele1_seedIx",&ix); 
  ntu_MC->SetBranchAddress("ele1_seedIy",&iy); 
  ntu_MC->SetBranchAddress("ele1_seedZside",&iz); 
  ntu_MC->SetBranchAddress("ele1_tkP",&tkP); 


  //---- Set branch addresses for Data
  ntu_Data->SetBranchAddress("PV_z", &PV_z);
  ntu_Data->SetBranchAddress("ele1_scEta", &scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_Data->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_Data->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_Data->SetBranchAddress("ele1_scE", &scE);
  ntu_Data->SetBranchAddress("ele1_scE_regression", &scE_regression);
  ntu_Data->SetBranchAddress("ele1_es", &esE);
  ntu_Data->SetBranchAddress("ele1_charge", &charge);
  ntu_Data->SetBranchAddress("ele1_fbrem",&fbrem); 
  ntu_Data->SetBranchAddress("ele1_seedIphi",&iphi); 
  ntu_Data->SetBranchAddress("ele1_seedIeta",&ieta); 
  ntu_Data->SetBranchAddress("ele1_seedIx",&ix); 
  ntu_Data->SetBranchAddress("ele1_seedIy",&iy); 
  ntu_Data->SetBranchAddress("ele1_seedZside",&iz); 
  ntu_Data->SetBranchAddress("ele1_tkP",&tkP); 

  // Analysis bins
  // Int_t iEtaBins = 250; 
  Int_t iEtaBins = 236; // there are eight rings per endcap side w/o electrons (beyond TK acceptance)
  Int_t nBins = iEtaBins/1. ;
  std::cout << "nBins = " << nBins << std::endl;
  
  TProfile *pr = new TProfile("a","a",nBins,-1.*etaMax, etaMax, 0, 2);
  TProfile *fbremMC   = new TProfile("fbremMC","fbremMC",nBins,-1.*etaMax, etaMax, 0, 2);
  TProfile *fbremData = new TProfile("fbremData","fbremData",nBins,-1.*etaMax, etaMax, 0, 2);
  int rebin = 4;

  // MC histos
  TH1F** h_EoP_MC = new TH1F*[nBins];   
  TH1F** h_EoC_MC = new TH1F*[nBins];   
  TH1F** h_Eta = new TH1F*[nBins]; // used to map iEta (as defined for Barrel and Endcap geom) into eta 
  for(int i = 0; i < nBins; ++i) {
    char histoName[80];
    sprintf(histoName, "EoP_MC_%d", i);
    h_EoP_MC[i] = new TH1F(histoName, histoName, 1200, 0., 3.);
    h_EoP_MC[i] -> SetFillColor(4);
    h_EoP_MC[i] -> SetFillStyle(3004);
    sprintf(histoName, "EoC_MC_%d", i);
    h_EoC_MC[i] = new TH1F(histoName, histoName, 1200, 0., 3.);
    h_EoC_MC[i] -> SetFillColor(3);
    h_EoC_MC[i] -> SetFillStyle(3004);
    sprintf(histoName, "Eta_%d", i);
    h_Eta[i] = new TH1F(histoName, histoName, 1000, -250., 250.); 
  }
  // DATA histos
  TH1F** h_EoP_Data = new TH1F*[nBins];   
  TH1F** h_EoC_Data = new TH1F*[nBins];   
  for(int i = 0; i < nBins; ++i){
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

  //---- book templates
  // [0] --> uncorrected
  // [1] --> corrected (regression)
  std::vector<int> refId(nBins);
  const int Ntempl = 8;
  TH1F* h_template_MC[Ntempl][2];
  TH1F* h_template_Data[Ntempl][2];
  for(int i = 0; i < Ntempl; ++i){
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
  
  //--- Initialize endcap geometry
  TEndcapRings *eRings = new TEndcapRings(); 

 
  //******************************************************************************************
  //*************************************** MC  ********************************************** 
   
  std::cout << "Loop in MC events " << endl; 
  float ww = 1;
  //--------------- loop on MC, make reference and fit distributions
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%100000 == 0 ) std::cout << "reading MC saved entry " << entry << "\r" << std::flush;
    //if (entry>2000000) break;

    ntu_MC->GetEntry(entry);
      
    R9 = scE3x3/scE;
    
    //-- eta and R9 cuts
    if ( fabs(scEta) > etaMax ) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 
    
    //-- remove phi cracks 
        // if (abs(ieta)<86) {
    //   float phi = (scPhi+PI)/xtalWidth;
    //   float modphi = (int)phi%20;
    //   if (fabs(modphi-10)<3.) continue;
    // }

    // -- PU weights
    if (usePUweights) ww = w[npu];
   
    //--- set ieta for the Endcaps
    if (iz != 0) {
      ieta = eRings->GetEndcapIeta(ix,iy,iz); 
    }
    if (abs(ieta)>iEtaBins/2) continue; 
    
    //-- fill templates for each mod
    int mod = templIndex(ieta);

    if (iz!=0) EoP = (scE-esE)/(tkP-esE);

    float correction = scE_regression/scE;
    if (applyPcalibration) EoP*=1./pCalib->GetMomentumCalibration(ieta,0);


    // reference out of the gaps
    if ( (ieta>3 && ieta<22)  || (ieta>27 && ieta<43) || 
	 (ieta>47 && ieta<63) || (ieta>67 && ieta<83) || (ieta>89 && ieta <122) )  
      //if( IsEtaGap(ieta)==0 ) 
      {
	h_template_MC[mod][0]-> Fill(EoP,ww);
	h_template_MC[mod][1]-> Fill(EoP*correction,ww);
      }
  

    // fill MC histos in eta bins
    int bin = (ieta+iEtaBins/2) * ((float)nBins/(float)iEtaBins); 
    if (ieta>0) bin = bin-1; 
    if (bin>nBins-1 || bin < 0 ) {
      cout << "Error in bins with MC: "<< mod << " " << bin <<" "<< " " << ieta << " " << scEta << flush;
    }
    
    h_EoP_MC[bin] -> Fill(EoP,ww);
    h_EoC_MC[bin] -> Fill(EoP*correction,ww);
    
    fbremMC->Fill(scEta,fbrem,ww);

    refId.at(bin) = mod; 

  } 

  
  //******************************************************************************************
  //*************************************** DATA ********************************************** 


  std::cout << "Loop in Data events " << endl; 
  //---- loop on Data
  for(int entry = 0; entry < ntu_Data->GetEntries(); ++entry) {
    if( entry%100000 == 0 ) std::cout << "reading data saved entry "  << entry << "\r" << std::flush;
    //if (entry>2000000) break;
    
    ntu_Data->GetEntry(entry);
    
       
    R9 = scE3x3/scE;

    //-- eta or R9 cuts
    if ( fabs(scEta) > etaMax) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 

    //-- remove phi cracks
    // if (abs(ieta)<86) {
    //   float phi = (scPhi+PI)/xtalWidth;
    //   float modphi = (int)phi%20;
    //   if (fabs(modphi-10)<3.) continue;
    // }
   
    //--- set ieta for the Endcaps
    if (iz != 0) 
      ieta = eRings->GetEndcapIeta(ix,iy,iz); 
    if (abs(ieta)>iEtaBins/2) continue; 

    //-- fill template for each mod
    int mod = templIndex(ieta);
    float correction = scE_regression/scE;
    if (applyPcalibration) EoP*=1./pCalib->GetMomentumCalibration(ieta,1);

    // reference out of the gaps
    if ( (ieta>3 && ieta<22)  || (ieta>27 && ieta<43) || 
	 (ieta>47 && ieta<63) || (ieta>67 && ieta<83) || (ieta>89 && ieta <122) )  
      //    if( IsEtaGap(ieta)==0 ) 
      {
	h_template_Data[mod][0]-> Fill(EoP);
	h_template_Data[mod][1]-> Fill(EoP*correction);
      }


    
    // fill Data histos in eta bins
    int bin = (ieta+iEtaBins/2) * ((float)nBins/(float)iEtaBins); 
    if (ieta>0) bin = bin-1; 
    if (bin>nBins-1 || bin < 0 ) {
      cout << "Error in bins with Data: " << bin <<" "<< " " << " " << scEta << endl;
    }


 
    h_EoP_Data[bin] -> Fill(EoP);
    h_EoC_Data[bin] -> Fill(EoP*correction);

    
    h_Eta[bin] -> Fill((double)ieta); 

    fbremData->Fill(scEta,fbrem);
  } 
  
  //************************************* FITTING ***************************************************//
  
  TGraphErrors* g_EoP_MC   = new TGraphErrors();g_EoP_MC->SetName("gEoP_MC");
  TGraphErrors* g_EoC_MC   = new TGraphErrors();g_EoC_MC->SetName("gEoC_MC");
  TGraphErrors* g_Rat_MC   = new TGraphErrors();g_Rat_MC->SetName("gCorr_Uncorr_MC");
  
  histoFunc *templateHistoFuncMC[8][2]; 
  histoFunc *templateHistoFuncData[8][2]; 
  for (int mod=0; mod<8;mod++) {
    for (int ll = 0; ll < 2; ll++){
      h_template_MC[mod][ll] -> Rebin(rebin);
      templateHistoFuncMC[mod][ll] = new histoFunc(h_template_MC[mod][ll]);
      h_template_Data[mod][ll] -> Rebin(rebin);
      templateHistoFuncData[mod][ll] = new histoFunc(h_template_Data[mod][ll]);
    }
  }
  
  TF1 * templateFunc;
  double xNorm;

  for(int i = 0; i < nBins; ++i) {
    
    h_EoP_MC[i] -> Rebin(rebin);    
    h_EoC_MC[i] -> Rebin(rebin);    
    
    int mod = refId.at(i); 
    //float xval = pr->GetXaxis()->GetBinCenter(i+1);
    float xval = i * (iEtaBins/nBins) - iEtaBins/2 + 0.5;
    xval = h_Eta[i]->GetMean(); 

    if (xval == 0) {
      cout << "xval = "<< xval << endl;  
      continue;
    }

    //-- MC uncorrected    
    std::cout << "***** Fitting uncorrected MC:  " << i << " " << mod << endl; 
    templateFunc = new TF1("templateFunc", templateHistoFuncMC[mod][0], 0.6, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);
    
    xNorm = h_EoP_MC[i]->Integral()/h_template_MC[mod][0]->Integral() *
            h_EoP_MC[i]->GetBinWidth(1)/h_template_MC[mod][0]->GetBinWidth(1); 
    
    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
     
    TFitResultPtr rp;
    int fStatus; 
    for (int trial=0;trial<5;trial++) {
      rp = h_EoP_MC[i] -> Fit("templateFunc", "MQRLS+");
      fStatus = rp;
      if (fStatus < 2) break; 
    }
    g_EoP_MC -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    //if(fStatus%10 == 4)  g_EoP_MC -> SetPointError(i, 2*etaMax/nBins, 10);
    //if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    //cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //---ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    //-- MC corrected    
    std::cout << "***** Fitting corrected MC:  " << i << " " << mod << endl; 
    templateFunc = new TF1("templateFunc", templateHistoFuncMC[mod][1], 0.6, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    xNorm = h_EoC_MC[i]->Integral()/h_template_MC[mod][1]->Integral() *
            h_EoC_MC[i]->GetBinWidth(1)/h_template_MC[mod][1]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
   
    for (int trial=0;trial<5;trial++) {
      rp = h_EoC_MC[i] -> Fit("templateFunc", "MQRLS+");
      fStatus = rp;
      if (fStatus < 2) break; 
    }

    g_EoC_MC -> SetPoint(i, xval , 1./templateFunc->GetParameter(1));
    g_EoC_MC -> SetPointError(i, 0., templateFunc->GetParError(1));
    //if(fStatus%10 == 4)  g_EoC_MC -> SetPointError(i, 2*etaMax/nBins, 10);
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
  
  
  for(int i = 0; i < nBins; ++i)
  {
    h_EoP_Data[i] -> Rebin(rebin);    
    h_EoC_Data[i] -> Rebin(rebin);    
    
    int mod = refId.at(i); 
    //float xval = pr->GetXaxis()->GetBinCenter(i+1);
    float xval = i * (iEtaBins/nBins) - iEtaBins/2 + 0.5;
    xval = h_Eta[i]->GetMean(); 
    
    //-- DATA uncorrected    
    std::cout << "***** Fitting uncorrected DATA:  " << i << " " << mod << endl; 
    templateFunc = new TF1("templateFunc", templateHistoFuncData[mod][0], 0.6, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);
    xNorm = h_EoP_Data[i]->Integral()/h_template_Data[mod][0]->Integral() *
            h_EoP_Data[i]->GetBinWidth(1)/h_template_Data[mod][0]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);

    TFitResultPtr rp;
    int fStatus; 
    for (int trial=0;trial<5;trial++) {
      rp = h_EoP_Data[i] -> Fit("templateFunc", "MQRLS+");
      fStatus = rp;
      if (fStatus < 2) break; 
    }
    g_EoP_Data -> SetPoint(i,  xval , 1./templateFunc->GetParameter(1));
    g_EoP_Data -> SetPointError(i, 0., templateFunc->GetParError(1));
    //    if(fStatus%10 == 4)  g_EoP_Data -> SetPointError(i, 2*etaMax/nBins, 10);
    //    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
    //cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 

    //--ratio preparation
    float rat = templateFunc->GetParameter(1);
    float era = templateFunc->GetParError(1); 

    //--- corrected DATA    
    std::cout << "***** Fitting corrected DATA:  " << i << " " << mod << endl; 
    templateFunc = new TF1("templateFunc", templateHistoFuncData[mod][1], 0.6, 1.3, 3, "histoFunc");
    templateFunc -> SetNpx(10000);

    xNorm = h_EoC_Data[i]->Integral()/h_template_Data[mod][1]->Integral() *
            h_EoC_Data[i]->GetBinWidth(1)/h_template_Data[mod][1]->GetBinWidth(1); 

    templateFunc -> FixParameter(0, xNorm);
    templateFunc -> SetParameter(1, 1.05 );
    templateFunc -> FixParameter(2, 0.);
     
    for (int trial=0;trial<5;trial++) {
      rp = h_EoC_Data[i] -> Fit("templateFunc", "MQRLS+");
      fStatus = rp;
      if (fStatus < 2) break; 
    }
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

  for (int imod = 0; imod<4; imod++){
    h_template_Data[imod][0]->Write();  
    h_template_Data[imod][1]->Write();
    h_template_MC[imod][1]->Write();  
    h_template_MC[imod][1]->Write();  
  }

  // for (int i=0; i<nBins;i++ )
  //   {
  //     h_EoP_MC[i]->Write() ;
  //     h_EoC_MC[i]->Write() ;
  //     h_EoP_Data[i]->Write() ;
  //     h_EoC_Data[i]->Write() ;
 

  //   }

  fbremMC->Write();
  fbremData->Write();

}
