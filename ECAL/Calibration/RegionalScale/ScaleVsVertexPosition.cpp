// per compilare: g++ -Wall -o ScaleVsVertexPosition `root-config --cflags --glibs` ScaleVsVertexPosition.cpp

#include "../CommonTools/TMomentumCalibration.h"
#include "../CommonTools/TEndcapRings.h"
#include "../CommonTools/histoFunc.h"
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

bool IsEtaGap(float eta){
  float feta = fabs(eta);
  if( fabs(feta - 0 )<3) return true;
  if( fabs(feta - 25)<3) return true;
  if( fabs(feta - 45)<3) return true;
  if( fabs(feta - 65)<3) return true;
  if( fabs(feta - 85)<3) return true;
  return false;
}

int templIndex(float eta){
    float feta = fabs(eta);
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

int zVertexBin(float zvtx, int nBins){
  float bin = -1;
  float maxz = 10;
  if ( zvtx < -maxz) bin = 0;
  if ( zvtx >  maxz) bin = nBins-1;
  
  float dz = 2*maxz/(nBins-2);

  for (int i = 0 ; i < (nBins-2); i++ ){
    if (zvtx >= (-maxz + i*dz )  &&  zvtx < (-maxz + (i+1)*dz ) ) bin = i+1;
  }
  return (bin);
}


//**************  MAIN PROGRAM **************************************************************
int main(int argc, char** argv)
{
  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"outputFiles/RelativeScaleVSVertexPosition_mod4_EBM_noPcalibration.root");
  
  //--- select module   
  int selectedMod = 3;

  //--- P calibration
  TMomentumCalibration *pCalib = new TMomentumCalibration("./outputFiles/momentumScale_mod4_EBP.root"); 


  //---- variables for selection
  float etaMin = -1.4442;
  float etaMax = -1.2;
  
  bool usePUweights = true;
  bool applyPcalibration = false;

  //--- weights for MC
  TFile weightsFile("../CommonTools/weights/PUweights_2011_0100_73500_WJetsToLL_Fall11_S6.root","READ"); 
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();


  //----- NTUPLES--------------------
  TChain *ntu_DA = new TChain("ntu");
  TChain *ntu_MC = new TChain("ntu");

  //Data 
  ntu_DA->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
   
  // --- MC Fall 2011
  ntu_MC->Add("../NTUPLES/Fall11/WZAnalysis/WZAnalysis_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1.root");
  
  std::cout << "     DATA: " << ntu_DA->GetEntries() << " entries in Data sample" << std::endl;
  std::cout << "     MC  : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;

  // observables  
  int isW;
  float PV_z;
  float EoP, scEta, scPhi, mZ;
  float scEta2,scEne2,scPhi2;
  float scE3x3, scE5x5, scEne, scERaw, scEt;  
  float charge, scLocalEta, scLocalPhi,crackCorr,localCorr; 
  float pTK,pTK2; 
  float scEneReg,scEneReg2;
  int iphi,ieta,ix,iy,iz; 
  int iphi2,ieta2,ix2,iy2,iz2;
  int npu;
  
  // Set branch addresses for Data  
  ntu_DA->SetBranchAddress("isW", &isW);
  ntu_DA->SetBranchAddress("PV_z", &PV_z);
  ntu_DA->SetBranchAddress("ele1ele2_scM",     &mZ);
  ntu_DA->SetBranchAddress("ele1_scEta", &scEta);
  ntu_DA->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_DA->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_DA->SetBranchAddress("ele2_scPhi", &scPhi2);
  ntu_DA->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_DA->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_DA->SetBranchAddress("ele1_e5x5", &scE5x5);
  ntu_DA->SetBranchAddress("ele1_scERaw", &scERaw);
  ntu_DA->SetBranchAddress("ele1_scE", &scEne);
  ntu_DA->SetBranchAddress("ele2_scE", &scEne2);
  ntu_DA->SetBranchAddress("ele1_scE_regression", &scEneReg);
  ntu_DA->SetBranchAddress("ele2_scE_regression", &scEneReg2);
  ntu_DA->SetBranchAddress("ele1_charge", &charge);
  ntu_DA->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_DA->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_DA->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_DA->SetBranchAddress("ele1_scLocalContCorr",&localCorr); 
  ntu_DA->SetBranchAddress("ele1_scEt",&scEt); 
  ntu_DA->SetBranchAddress("ele1_tkP",&pTK); 
  ntu_DA->SetBranchAddress("ele2_tkP",&pTK2); 
  ntu_DA->SetBranchAddress("ele1_seedIphi",&iphi); 
  ntu_DA->SetBranchAddress("ele1_seedIeta",&ieta); 
  ntu_DA->SetBranchAddress("ele1_seedIx",&ix); 
  ntu_DA->SetBranchAddress("ele1_seedIy",&iy); 
  ntu_DA->SetBranchAddress("ele1_seedZside",&iz); 
  ntu_DA->SetBranchAddress("ele2_seedIphi",&iphi2); 
  ntu_DA->SetBranchAddress("ele2_seedIeta",&ieta2); 
  ntu_DA->SetBranchAddress("ele2_seedIx",&ix2); 
  ntu_DA->SetBranchAddress("ele2_seedIy",&iy2); 
  ntu_DA->SetBranchAddress("ele2_seedZside",&iz2); 

  // Set branch addresses for MC  
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("isW", &isW);
  ntu_MC->SetBranchAddress("PV_z", &PV_z);
  ntu_MC->SetBranchAddress("ele1ele2_scM",     &mZ);
  ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_MC->SetBranchAddress("ele1_e5x5", &scE5x5);
  ntu_MC->SetBranchAddress("ele1_scE", &scEne);
  ntu_MC->SetBranchAddress("ele2_scE", &scEne2);
  ntu_MC->SetBranchAddress("ele1_scERaw", &scERaw);
  ntu_MC->SetBranchAddress("ele1_charge", &charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_MC->SetBranchAddress("ele1_scEt",&scEt); 
  ntu_MC->SetBranchAddress("ele1_tkP",&pTK); 
  ntu_MC->SetBranchAddress("ele2_tkP",&pTK2); 
  ntu_MC->SetBranchAddress("ele1_seedIphi",&iphi); 
  ntu_MC->SetBranchAddress("ele1_seedIeta",&ieta); 
  ntu_MC->SetBranchAddress("ele1_seedIx",&ix); 
  ntu_MC->SetBranchAddress("ele1_seedIy",&iy); 
  ntu_MC->SetBranchAddress("ele1_seedZside",&iz); 
  
  // Analysis bins
  Int_t nBins = 16 ;
  std::cout << "nBins = " << nBins << std::endl;

  // histogram definition
  TH1F** h_EoP = new TH1F*[nBins];   
  TH1F** h_EoC = new TH1F*[nBins];   
  TF1** f_EoP = new TF1*[nBins];
  TF1** f_EoC = new TF1*[nBins];
  std::vector<int> refId(nBins);
  TH1F** h_zvtx = new TH1F*[nBins];

  for(Int_t i = 0; i < nBins; ++i)
  {
    char histoName[80];
    sprintf(histoName, "EoP_%d", i);
    h_EoP[i] = new TH1F(histoName, histoName, 1200, 0., 3.);
    h_EoP[i] -> SetFillColor(kRed+2);
    h_EoP[i] -> SetLineColor(kRed+2);
    h_EoP[i] -> SetFillStyle(3004);

    sprintf(histoName, "EoC_%d", i);
    h_EoC[i] = new TH1F(histoName, histoName, 1200, 0., 3.);
    h_EoC[i] -> SetFillColor(kGreen+2);
    h_EoC[i] -> SetLineColor(kGreen+2);
    h_EoC[i] -> SetFillStyle(3004);

    sprintf(histoName, "hz_%d", i);
    h_zvtx[i] = new TH1F(histoName, histoName, 1000, -50., 50.); 
  }

  TH1F* h_template[8];
  for (unsigned int imod = 0; imod<8; imod++){
    char histoName[80];
    sprintf(histoName, "template_%d", imod);
    h_template[imod] = new TH1F(histoName, "", 1200, 0., 3.);  
  }

  //**************************** loop on MC, make refernce and fit dist
  float ww = 1;
  float var = 0;
  std::cout << "Loop in MC events " << endl; 
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%10000 == 0 ) std::cout << "reading saved entry " << entry << "\r" << std::flush;
    //    if (entry>1000) break;

    ntu_MC->GetEntry(entry);

    //--- PU weights
    if (usePUweights) ww = w[npu];


    //--- only one module
    int mod  = templIndex(ieta);
    if ( mod != selectedMod ) continue;
    if ( scEta < etaMin || scEta > etaMax ) continue;

    //--- define bins in z vertex
    int bin = zVertexBin(PV_z, nBins);
    if (bin>nBins-1 || bin < 0 ) {
      cout << "Error in bins with MC: "<< bin << " " << PV_z << endl;
    }  
    
    if (applyPcalibration) EoP*=1./pCalib->GetMomentumCalibration(PV_z,0);
    var = EoP; 
   
    //--- reference out of the gap
    ieta = fabs(ieta);
    if ( (ieta>3 && ieta<22)  || (ieta>27 && ieta<43) || 
	 (ieta>47 && ieta<63) || (ieta>67 && ieta<83) || (ieta>89 && ieta <122) )  
      h_template[mod] ->  Fill(var,ww);
    
    h_EoP[bin] -> Fill(var,ww);  // This is MC

    refId.at(bin) = mod; 

  } 

  // loop on events in Data
  std::cout << "Loop in Data events " << endl; 
  for(int entry = 0; entry < ntu_DA->GetEntries(); ++entry) {
    if( entry%10000 == 0 ) std::cout << "reading saved entry " << entry << "\r" << std::flush;
    //    if (entry>1000) break;

    ntu_DA->GetEntry(entry);

    //--- only one module
    int mod  = templIndex(ieta);
    if ( mod != selectedMod ) continue;
    if ( scEta < etaMin || scEta > etaMax ) continue;


    //--- define bins in z vertex
    int bin = zVertexBin(PV_z,nBins);
    if (bin>nBins-1 || bin < 0 ) {
      cout << "Error in bins with MC: " << bin << " " << PV_z << endl;
    }  

    
    if (applyPcalibration) EoP*=1./pCalib->GetMomentumCalibration(PV_z,1);
    var  = EoP;
        
    (h_EoC[bin]) -> Fill(var);  // This is Data
    (h_zvtx[bin]) -> Fill( PV_z ); 
      
  }

  std::cout << "End loop: Analyze events " << endl; 
    
  // draw results
  TGraphErrors* g_EoP   = new TGraphErrors();
  TGraphErrors* g_EoC   = new TGraphErrors();
  TGraphErrors* g_Rat   = new TGraphErrors();
 
  int rebin = 6;
  
  histoFunc *templateHistoFunc[8]; 
  for (int mod=0; mod<8;mod++) {
    h_template[mod] -> Rebin(rebin*2);
    templateHistoFunc[mod] = new histoFunc(h_template[mod]);
  }



  for(unsigned int i = 0; i < nBins; ++i)
  {
    h_EoP[i] -> Rebin(rebin);    
    h_EoC[i] -> Rebin(rebin);    
    
    // define the fitting function
    // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )
    int mod = refId.at(i); 
    
    char funcName[50];
    sprintf(funcName,"f_EoP_%d_Ref_%d",i,mod);
    cout << funcName << endl; 
    f_EoP[i] = new TF1(funcName, templateHistoFunc[mod], 0.7, 1.3, 3, "histoFunc");
    f_EoP[i] -> SetParName(0,"Norm"); 
    f_EoP[i] -> SetParName(1,"Scale factor"); 
    f_EoP[i] -> SetLineWidth(1); 
    f_EoP[i] -> SetLineColor(kRed+2); 
    f_EoP[i] -> SetNpx(10000);

    // uncorrected    
    double xNorm = h_EoP[i]->Integral()/h_template[mod]->Integral() *
                   h_EoP[i]->GetBinWidth(1)/h_template[mod]->GetBinWidth(1); 

    f_EoP[i] -> FixParameter(0, xNorm);
    f_EoP[i] -> SetParameter(1, 1.05 );
    f_EoP[i] -> FixParameter(2, 0.);

    TFitResultPtr rp;
    int fStatus; 
    for (int trial=0;trial<5;trial++) {
      f_EoP[i] -> SetParameter(1, (0.99 + 0.01*i) );
      rp = h_EoP[i] -> Fit(funcName, "QMRL+");
      fStatus = rp;
      if (fStatus < 2) break; 
    }
    
    float zbin = h_zvtx[i]->GetMean(); 
    float ezbin = h_zvtx[i]->GetRMS(); 
    g_EoP -> SetPoint(i, zbin , 1./f_EoP[i]->GetParameter(1));
    g_EoP -> SetPointError(i, ezbin, f_EoP[i]->GetParError(1));
    //cout  << " ***** " <<  1./f_EoP[i]->GetParameter(1) << " " << f_EoP[i]->GetParError(1) << " " << fStatus <<  endl; 

    //ratio preparation
    float ymc  = 1./f_EoP[i]->GetParameter(1);
    float eymc = f_EoP[i]->GetParError(1); 

    // corrected    
    xNorm = h_EoC[i]->Integral()/h_template[mod]->Integral() *
            h_EoC[i]->GetBinWidth(1)/h_template[mod]->GetBinWidth(1); 

    sprintf(funcName,"f_EoC_%d_Ref_%d",i,mod);
    f_EoC[i] = new TF1(funcName, templateHistoFunc[mod], 0.7, 1.3, 3, "histoFunc");
    f_EoC[i] -> SetParName(0,"Norm"); 
    f_EoC[i] -> SetParName(1,"Scale factor"); 
    f_EoC[i] -> SetLineWidth(1); 
    f_EoC[i] -> SetLineColor(kGreen+2); 
    f_EoC[i] -> SetNpx(10000);

    f_EoC[i] -> FixParameter(0, xNorm);
    f_EoC[i] -> SetParameter(1, 1.05 );
    f_EoC[i] -> FixParameter(2, 0.);
    
    std::cout << "***** Re-Fitting ";
    for (int trial=0;trial<2;trial++) {
      f_EoC[i] -> SetParameter(1, (0.99 + 0.01*i) );
      rp = h_EoC[i] -> Fit(funcName, "QMRL+");
      fStatus = rp;
      if (fStatus < 2) break; 
    }

    g_EoC -> SetPoint(i, zbin, 1./f_EoC[i]->GetParameter(1));
    g_EoC -> SetPointError(i, ezbin, f_EoC[i]->GetParError(1));
    //    cout << " ********** " <<  1./f_EoC[i]->GetParameter(1) << " " << f_EoC[i]->GetParError(1) << endl; 

    float yda  = 1./f_EoC[i]->GetParameter(1);
    float eyda = f_EoC[i]->GetParError(1); 

    //ratio finalization: data/MC    
    float rat = yda/ymc;
    float era = rat*sqrt(eymc*eymc+eyda*eyda); 
    
    g_Rat -> SetPoint(i, zbin , rat); 
    g_Rat -> SetPointError(i,  0. , era); 
    g_Rat->SetLineColor(kBlue+2); 
    
  }

  TFile* o = new TFile(outfilename,"RECREATE");
  o -> cd();
  
  g_EoP -> Write("g_EoP");
  g_EoC -> Write("g_EoC");
  g_Rat -> Write("g_Rat");
  
  for (int imod = 0; imod<8; imod++){
    h_template[imod]->Write();  
  }
  
  for (int i=0; i<nBins;i++ )
    {
      h_EoP[i]->Write() ;
      h_EoC[i]->Write() ;
      h_zvtx[i]->Write();
    }

  o -> Close();
    

}
