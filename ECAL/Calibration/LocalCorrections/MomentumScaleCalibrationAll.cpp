// per compilare: g++ -Wall -o MomentumScaleCalibrationAll `root-config --cflags --glibs` MomentumScaleCalibrationAll.cpp

#include "./TEndcapRings.h"
#include "./histoFunc.h"
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

//**************  MAIN PROGRAM **************************************************************
int main(int argc, char** argv)
{
  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"outputFiles/momentumScale_scM.root");

  //---- variables for selection
  float r9min = 0.00 ;
  float r9max = 9999 ;  
  float etaMax  = 2.5;
  float eta2Max = 2.5;

  bool usePUweights = true;

  //--- weights for MC
  TFile weightsFile("weights/PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","READ"); 
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();


  //----- NTUPLES--------------------
  TChain *ntu_DA = new TChain("ntu");
  TChain *ntu_MC = new TChain("ntu");

  // Data 
  ntu_DA->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
 
  // --- MC Fall 2011
  ntu_MC->Add("../NTUPLES/Fall11/WZAnalysis/WZAnalysis_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1.root");


  std::cout << "     DATA: " << ntu_DA->GetEntries() << " entries in Data sample" << std::endl;
  std::cout << "     MC  : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;

  // observables  
    int isW;
  float EoP, scEta, scPhi, mZ;
  float scEta2,scEne2,scPhi2;
  float scE3x3, scE5x5, scEne, scERaw, zVtx, scEt;  
  float charge, scLocalEta, scLocalPhi,crackCorr,localCorr; 
  float pTK,pTK2; 
  float scEneReg,scEneReg2;
  int iphi,ieta,ix,iy,iz; 
  int iphi2,ieta2,ix2,iy2,iz2;
  int npu;
  
  // Set branch addresses for Data  
  ntu_DA->SetBranchAddress("isW", &isW);
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
  ntu_DA->SetBranchAddress("PV_z",&zVtx); 
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
  ntu_MC->SetBranchAddress("PV_z",&zVtx); 
  ntu_MC->SetBranchAddress("ele1_seedIphi",&iphi); 
  ntu_MC->SetBranchAddress("ele1_seedIeta",&ieta); 
  ntu_MC->SetBranchAddress("ele1_seedIx",&ix); 
  ntu_MC->SetBranchAddress("ele1_seedIy",&iy); 
  ntu_MC->SetBranchAddress("ele1_seedZside",&iz); 
  
  // Analysis bins
  //  Int_t iEtaBins = 250; 
  Int_t iEtaBins = 236; // there are eight rings per endcap side w/o electrons (beyond TK acceptance)
  Int_t nBins = iEtaBins/1. ;
  std::cout << "nBins = " << nBins << std::endl;

  // histogram definition
  TH1F** h_EoP = new TH1F*[nBins];   
  TH1F** h_EoC = new TH1F*[nBins];   
  TH1F** h_Eta = new TH1F*[nBins]; // used to map iEta (as defined for Barrel and Endcap geom) into eta 
  TF1** f_EoP = new TF1*[nBins];
  TF1** f_EoC = new TF1*[nBins];
  std::vector<int> refId(nBins);

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

    sprintf(histoName, "Eta_%d", i);
    h_Eta[i] = new TH1F(histoName, histoName, 1000, -250., 250.); 
  }

  TH1F* h_template[8];
  for (unsigned int imod = 0; imod<8; imod++){
    char histoName[80];
    sprintf(histoName, "template_%d", imod);
    h_template[imod] = new TH1F(histoName, "", 1200, 0., 3.);  
  }

  
  TH1F* h_eta_data = new TH1F("h_eta_data","h_eta_data",300,-3,3);
  TH1F* h_eta_mc   = new TH1F("h_eta_mc","h_eta_mc",300,-3,3);
  TH1F* h_et_data  = new TH1F("h_et_data","h_et_data",100,0,100);
  TH1F* h_et_mc    = new TH1F("h_et_mc","h_et_mc",100,0,100);


  // Initialize endcap geometry
  TEndcapRings *eRings = new TEndcapRings(); 

  //**************************** loop on MC, make refernce and fit dist
  float ww = 1;
  float var = 0;
  std::cout << "Loop in MC events " << endl; 
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%10000 == 0 ) std::cout << "reading saved entry " << entry << "\r" << std::flush;
    //    if (entry>1000) break;

    ntu_MC->GetEntry(entry);
    
    if (isW==1) continue;
    if( fabs(scEta2) > eta2Max ) continue;
    float R9 = scE3x3/scEne;
    if ( R9 < r9min || R9 > r9max ) continue; 

    //--- PU weights
    if (usePUweights) ww = w[npu];
 
    //--- set ieta for the Endcaps
    if (iz != 0) {
      ieta = eRings->GetEndcapIeta(ix,iy,iz); 
    }
    if (abs(ieta)>iEtaBins/2) continue; 
    
    // //--- cut phi cracks
    // if (abs(ieta)<86) {
    //   float phi = (scPhi+PI)/xtalWidth;
    //   float modphi = (int)phi%20;
    //   if (fabs(modphi-10)<3.) continue;
    // }

    var = mZ * sqrt(pTK/scEne)/91.19;    /// use the momentum for ele1
    
    int mod  = templIndex(ieta);

    // reference out of the gap
    if ( (ieta>3 && ieta<22)  || (ieta>27 && ieta<43) || 
	 (ieta>47 && ieta<63) || (ieta>67 && ieta<83) || (ieta>89 && ieta <122) )  
      //if( IsEtaGap(ieta)==0 )
      h_template[mod] ->  Fill(var,ww);

    // fill MC histos in eta bins
    int bin = (ieta+iEtaBins/2) * (nBins/iEtaBins); 
    if (ieta>0) bin = bin-1; 
    if (bin>nBins-1 || bin < 0 ) {
      cout << "Error in bins with MC: "<< mod << " " << bin <<" "<< " " << ieta << " " << scEta << endl;
    }
    
    h_EoP[bin] -> Fill(var,ww);  // This is MC
    
    h_eta_mc->Fill(scEta,ww);
    h_et_mc ->Fill(scEt,ww);
    
    refId.at(bin) = mod; 

  } 

  // loop on events in Data
  std::cout << "Loop in Data events " << endl; 
  for(int entry = 0; entry < ntu_DA->GetEntries(); ++entry) {
    if( entry%10000 == 0 ) std::cout << "reading saved entry " << entry << "\r" << std::flush;
    //    if (entry>1000) break;

    ntu_DA->GetEntry(entry);

    if (isW==1) continue;
    if ( fabs(scEta2) > eta2Max) continue;
    float R9 = scE3x3/scEne;
    if ( R9 < r9min || R9 > r9max ) continue; 



    //--- set ieta for the Endcaps
    if (iz != 0) 
      ieta = eRings->GetEndcapIeta(ix,iy,iz); 
    if (abs(ieta)>iEtaBins/2) continue; 

    // //--- cut phi cracks
    // if (abs(ieta)<86) {
    //   float phi = (scPhi+PI)/xtalWidth;
    //   float modphi = (int)phi%20;
    //   if (fabs(modphi-10)<3.) continue;
    // }

    var  = mZ * sqrt(pTK/scEne) / 91.19;    /// use the momentum for ele1
    //var*=var;
       
    // fill Data histos in eta bins
    int bin = (ieta+iEtaBins/2) * (nBins/iEtaBins); 
    if (ieta>0) bin = bin-1; 
    if (bin>nBins-1 || bin < 0 ) {
      cout << "Error in bins with Data: " << bin <<" "<< " " << " " << scEta << endl;
    }

    (h_EoC[bin]) -> Fill(var);  // This is Data
    (h_Eta[bin]) -> Fill((double)ieta); 


    //use also the other electron
    var = mZ * sqrt(pTK2/scEne2) / 91.19;    /// use the momentum for ele2

    if (iz2 != 0) 
      ieta2 = eRings->GetEndcapIeta(ix2,iy2,iz2); 
    if (abs(ieta2)>iEtaBins/2) continue; 

    // cut phi cracks in EB
    //    if (abs(ieta2)<86) {
    //     float phi = (scPhi2+3.1415926536)/0.01745329;
    //     float modphi = (int)phi%20;
    //     if (fabs(modphi-10)<3.) continue;
    //    }

    // fill Data histos in eta bins
    bin = (ieta2+iEtaBins/2) * (nBins/iEtaBins); 
    if (ieta2>0) bin = bin-1; 
    if (bin>nBins-1 || bin < 0 ) {
      cout << "Error in bins with Data: " << bin <<" "<< " " << " " << scEta << endl;
    }

    h_EoC[bin] -> Fill(var);  // This is Data
    h_Eta[bin] -> Fill((double)ieta2); 
    
    h_eta_data->Fill(scEta);
    h_et_data ->Fill(scEt);

  }

  std::cout << "End loop: Analyze events " << endl; 
    
  // draw results
  TGraphErrors* g_EoP   = new TGraphErrors();
  TGraphErrors* g_EoC   = new TGraphErrors();
  TGraphErrors* g_Rat   = new TGraphErrors();
 
  int rebin = 2;
  
  histoFunc *templateHistoFunc[8]; 
  for (int mod=0; mod<8;mod++) {
    h_template[mod] -> Rebin(rebin*4);
    templateHistoFunc[mod] = new histoFunc(h_template[mod]);
  }



  for(unsigned int i = 0; i < nBins; ++i)
  {
    h_EoP[i] -> Rebin(rebin*4);    
    h_EoC[i] -> Rebin(rebin*4);    
    
    // define the fitting function
    // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )
    int mod = refId.at(i); 
    
    char funcName[50];
    sprintf(funcName,"f_EoP_%d_Ref_%d",i,mod);
    cout << funcName << endl; 
    f_EoP[i] = new TF1(funcName, templateHistoFunc[mod], 0.6, 1.3, 3, "histoFunc");
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
      rp = h_EoP[i] -> Fit(funcName, "QMRL+");
      fStatus = rp;
      if (fStatus < 2) break; 
    }
    
    //    int bin = (ieta+iEtaBins/2) * (nBins/iEtaBins);  //questa da invertire   
    float flEta = i * (iEtaBins/nBins) - iEtaBins/2 + 0.5;
    flEta = h_Eta[i]->GetMean(); 
        
    // g_EoP -> SetPoint(i, flEta , 1./f_EoP[i]->GetParameter(1));
    // g_EoP -> SetPointError(i, 0., f_EoP[i]->GetParError(1));
    g_EoP -> SetPoint(i, flEta , pow(f_EoP[i]->GetParameter(1),2));
    g_EoP -> SetPointError(i, 0., 2*f_EoP[i]->GetParError(1));
    //cout  << " ***** " <<  1./f_EoP[i]->GetParameter(1) << " " << f_EoP[i]->GetParError(1) << " " << fStatus <<  endl; 

    //ratio preparation
    float rat = f_EoP[i]->GetParameter(1);
    float era = f_EoP[i]->GetParError(1); 

    // corrected    
    xNorm = h_EoC[i]->Integral()/h_template[mod]->Integral() *
            h_EoC[i]->GetBinWidth(1)/h_template[mod]->GetBinWidth(1); 

    sprintf(funcName,"f_EoC_%d_Ref_%d",i,mod);
    f_EoC[i] = new TF1(funcName, templateHistoFunc[mod], 0.6, 1.3, 3, "histoFunc");
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
      rp = h_EoC[i] -> Fit(funcName, "QMRL+");
      fStatus = rp;
      if (fStatus < 2) break; 
    }

    // g_EoC -> SetPoint(i, flEta, 1./f_EoC[i]->GetParameter(1));
    // g_EoC -> SetPointError(i, 0., f_EoC[i]->GetParError(1));
    g_EoC -> SetPoint(i, flEta, pow(f_EoC[i]->GetParameter(1),2));
    g_EoC -> SetPointError(i, 0., 2*f_EoC[i]->GetParError(1));
    //    cout << " ********** " <<  1./f_EoC[i]->GetParameter(1) << " " << f_EoC[i]->GetParError(1) << endl; 

    //ratio finalization
    rat /= f_EoC[i]->GetParameter(1);
    //    rat = 1+2.*(rat-1);
    era = rat*sqrt(era*era+f_EoC[i]->GetParError(1)*f_EoC[i]->GetParError(1)); 
    
    g_Rat -> SetPoint(i, flEta , rat); 
    g_Rat -> SetPointError(i,  0. , era); 
    g_Rat->SetLineColor(kBlue+2); 
    
  }

  TFile* o = new TFile(outfilename,"RECREATE");
  o -> cd();
  
  g_EoP -> Write("g_EoP");
  g_EoC -> Write("g_EoC");
  g_Rat -> Write("g_Rat");
  
  // for (int imod = 0; imod<8; imod++){
  //   h_template[imod]->Write();  
  //   h_template[imod]->Write();  
  // }
  
  // for (int i=0; i<nBins;i++ )
  //   {
  //     h_EoP[i]->Write() ;
  //     h_EoC[i]->Write() ;
  //   }

  h_eta_mc->Write();
  h_et_mc->Write();
  h_eta_data->Write();
  h_et_data->Write();

  o -> Close();
    

}
