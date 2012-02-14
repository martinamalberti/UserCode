// per compilare: g++ -Wall -o TestZpeak `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit -lMathMore TestZpeak.cpp

#include "../CommonTools/TPileupReweighting.h"

#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooExponential.h"
#include "RooAbsReal.h"
#include "RooMsgService.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>

#define  xtalWidth  0.01745329
#define  zmass      91.188


using namespace std;
using namespace RooFit;
using std::vector;


TF1 *f = new TF1("f","[0]+([1]*(x)-[2]*pow(x,2)-[3]*pow(x,3)-[4]*pow(x,4))",-1,1.);

//****Effective sigma from fitted func ****************************************************************
/*void evalEffectiveSigma( RooAddPdf* pdf , RooRealVar*  mass, double &ml, double &mh, double &sigmaeff){

  double testmass = 90.;
  double center = testmass-10.0;
  
  double cdfhi = 0;
  double cdflo = 0;
  double mlmin = 0.0;
  double mhmin = 0.0;
  double step  = 0.05;

  double minwidth = 999.0;

  RooAbsReal *cdf = 0;
  double mlow = -999;
  double mhigh = -999;
  
  for (int i=0; i<400; ++i) {
    if (i%100==0) cout <<  i << endl;
    mlow = center+i*step;
    mass->setRange("signal",60,mlow); 
    cdf = pdf->createIntegral(*mass, NormSet(*mass),Range("signal"));
    cdflo = cdf->getVal();
    for (int j=i; j<400; ++j) {
      mhigh = center+j*step;
      mass->setRange("signal",60,mhigh); 
      cdf   = pdf->createIntegral(*mass, NormSet(*mass),Range("signal"));
      cdfhi = cdf->getVal();
      if ( (cdfhi-cdflo)>0.684 ) {
	if ( (mhigh-mlow)<minwidth) {
	  minwidth = mhigh-mlow;
	  mlmin = mlow;
	  mhmin = mhigh;
	}
	break;
      }
    }
  }
 
  sigmaeff = minwidth/2.0;
  ml = mlmin;
  mh = mhmin;
  cout << "Mmin = " << mlmin << "  Mmax = " << mhmin << "  effective sigma = " << sigmaeff << endl;
  
  // return (sigmaeff);
  return;
}
*/
//*****************************************************************************************************


//****Effective sigma from histogram ****************************************************************
void evalEffectiveSigmaFromHisto( TH1F *h , double &ml, double &mh, double &sigmaeff){

  double testmass = 90.;
  double center = testmass-6.0;
  
  double cdfhi = 0;
  double cdflo = 0;
  double mlmin = 0.0;
  double mhmin = 0.0;
  double step  = 0.02;

  double minwidth = 999.0;

  double mlow = -999;
  double mhigh = -999;
  int binlow, binhigh;
  int nbins = h->GetNbinsX();

  for (int i=0; i<600; ++i) {
    if (i%100==0) cout <<  i << endl;
    mlow = center+i*step;
    binlow = h ->FindBin(mlow);
    cdflo = h->Integral(1, binlow)/h->Integral(1, nbins);
    for (int j=i; j<600; ++j) {
      mhigh = center+j*step;
      binhigh = h ->FindBin(mhigh);
      cdfhi = h->Integral(1, binhigh)/h->Integral(1, nbins);
      if ( (cdfhi-cdflo)>0.684 ) {
	if ( (mhigh-mlow)<minwidth) {
	  minwidth = mhigh-mlow;
	  mlmin = mlow;
	  mhmin = mhigh;
	}
	break;
      }
    }
  }
  
  sigmaeff = minwidth/2.0;
  ml = mlmin;
  mh = mhmin;
  cout << "Mmin = " << mlmin << "  Mmax = " << mhmin << "  effective sigma = " << sigmaeff << endl;
  
  // return (sigmaeff);
  return;
}

//*****************************************************************************************************

int modId(float eta){
  int ieta  = fabs(eta)/0.01745329;
  int im = 0;
  if (ieta>25) im = 1;
  if (ieta>45) im = 2;
  if (ieta>65) im = 3;
  return (im);
}

bool IsEtaGap(float eta){
  int ieta  = fabs(eta)/xtalWidth;
  float feta = fabs(ieta);
  if( fabs(feta - 0 )<2) return true;
  if( fabs(feta - 25)<2) return true;
  if( fabs(feta - 45)<2) return true;
  if( fabs(feta - 65)<2) return true;
  if( fabs(feta - 85)<2) return true;
  return false;
}


//----------------------------------------------------------------------------
double basicClusterCorrection (float localEta, int imod){
  
  // --- correction function derived from E/p data
  if (imod==0) f->SetParameters(1.00672,0.00804331 , 0.0775959 , 0.0272849 , -0.0124545);
  if (imod==1) f->SetParameters(1.00734,0.00995546 , 0.0951636 , 0.0317456 , -0.0535642);
  if (imod==2) f->SetParameters(1.00844,0.0147664 , 0.123315 , 0.0446714 , -0.120907);
  if (imod==3) f->SetParameters(1.0132,0.0117694 , 0.210148 , 0.0306087 , -0.305);

  double corr = f-> Eval(localEta);
  return(1./corr);
}
//-----------------------------------------------------------------------------

double singleBasicClusterCorrection (float localEta, int imod){
   
  // --- correction function derived from E/p data
  if (imod==0) f->SetParameters(1.00613,0.0049502 , 0.0677101 , 0.0148299 , 0.00523093);
  if (imod==1) f->SetParameters(1.00633,0.00403076 , 0.0689901 , 0.00242213 , 0.0239275);
  if (imod==2) f->SetParameters(1.00648,0.0111289 , 0.0878569 , 0.0299266 , -0.0597349);
  if (imod==3) f->SetParameters(1.01058,0.0135225 , 0.16731 , 0.0208505 , -0.21014);

  double corr = f-> Eval(localEta);
  return(1./corr);

}

int main(int argc, char** argv)
{

  RooMsgService::instance().Print() ;
  RooMsgService::instance().setStreamStatus(1,0);

  gROOT->ProcessLine("#include <vector>");

  //
  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);
  
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.03);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.03);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.03);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);
  gROOT->ForceStyle();

  //---- variables for selections
  float etaMin = 0.;
  float etaMax = 1.444;
  
  float r9min = 0.00 ;
  float r9max = 9999;
  float minNPV = 0;
  float maxNPV = 100;

  bool useData = true;
  bool useMC   = false;
  bool usePUweights = false;

  bool useRegression = false;
  bool usePUcleaning = false;
  bool useSingleBasicClusterCorrectionEta = true;
  bool useLocalCorrectionEta = false;
  bool useLocalContCorrection = false;
 
  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"dataZpeak_allR9_nBC1_singleBCcorr.root");

  //---- PU weights for MC
  TPileupReweighting *puReweighting = new TPileupReweighting("../CommonTools/weights/PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","hweights"); 
  
  //---- NTUPLES 
  TChain *ntu_Data = new TChain("ntu");
  if (useData){
    // ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
    // ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
    // ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
    // ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
    // ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B_bc.root");
  }
  if (useMC){
    ntu_Data->Add("../NTUPLES/Fall11/WZAnalysis/WZAnalysis_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1.root");
  }             
      
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;
  std::cout << std::endl;
  
  //---- Set branch addresses  
  int npu;
  int runId,isZ;
  int PV_n;
  int ele1_isEB,ele2_isEB, ele1_isEBEEGap,ele2_isEBEEGap ;
  float mZ, mZc;
  float ele1_scE,ele2_scE;
  float ele1_scERaw,ele2_scERaw;
  float ele1_scERaw_PUcleaned,ele2_scERaw_PUcleaned;
  float ele1_e3x3, ele2_e3x3;
  float ele1_scE_regression, ele2_scE_regression;
  float ele1_fCorr, ele2_fCorr;
  float ele1_fCorr_PUcleaned, ele2_fCorr_PUcleaned;
  float ele1_scCrackCorr, ele2_scCrackCorr;
  float ele1_scLocalContCorr, ele2_scLocalContCorr;
  float ele1_scEta,ele2_scEta;
  float ele1_scPhi,ele2_scPhi;
  float ele1_eta,ele2_eta;
  float ele1_phi,ele2_phi;
  float ele1_scLocalEta, ele2_scLocalEta;
  float ele1_scLocalPhi, ele2_scLocalPhi;
  int   ele1_charge , ele2_charge;
  int ele1_bcN, ele2_bcN;
  vector<float>   *ele1_bcE = 0;
  vector<float>   *ele1_bcEta= 0;
  vector<float>   *ele1_bcPhi= 0;
  vector<float>   *ele1_bcLocalEta= 0;  
  vector<float>   *ele1_bcLocalPhi= 0;
  vector<float>   *ele2_bcE= 0;
  vector<float>   *ele2_bcEta= 0;
  vector<float>   *ele2_bcPhi= 0;
  vector<float>   *ele2_bcLocalEta= 0;  
  vector<float>   *ele2_bcLocalPhi= 0;

  ntu_Data->SetBranchAddress("PV_n", &PV_n);
  ntu_Data->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_Data->SetBranchAddress("isZ",                 &isZ);
  ntu_Data->SetBranchAddress("ele1ele2_scM",        &mZ);
  ntu_Data->SetBranchAddress("ele1_isEB",           &ele1_isEB);
  ntu_Data->SetBranchAddress("ele1_isEBEEGap",      &ele1_isEBEEGap);
  ntu_Data->SetBranchAddress("ele1_e3x3",           &ele1_e3x3);
  ntu_Data->SetBranchAddress("ele1_scERaw",         &ele1_scERaw);
  ntu_Data->SetBranchAddress("ele1_scERaw_PUcleaned",&ele1_scERaw_PUcleaned);
  ntu_Data->SetBranchAddress("ele1_scE",            &ele1_scE);
  ntu_Data->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression);
  ntu_Data->SetBranchAddress("ele1_fCorrection",    &ele1_fCorr);
  ntu_Data->SetBranchAddress("ele1_fCorrection_PUcleaned", &ele1_fCorr_PUcleaned);
  ntu_Data->SetBranchAddress("ele1_scCrackCorr",    &ele1_scCrackCorr);
  ntu_Data->SetBranchAddress("ele1_scLocalContCorr",    &ele1_scLocalContCorr);
  ntu_Data->SetBranchAddress("ele1_scEta",          &ele1_scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi",          &ele1_scPhi);
  ntu_Data->SetBranchAddress("ele1_eta",            &ele1_eta);
  ntu_Data->SetBranchAddress("ele1_phi",            &ele1_phi);
  ntu_Data->SetBranchAddress("ele1_scLocalPhi",     &ele1_scLocalPhi);
  ntu_Data->SetBranchAddress("ele1_scLocalEta",     &ele1_scLocalEta); 
  ntu_Data->SetBranchAddress("ele1_charge",         &ele1_charge); 
  ntu_Data->SetBranchAddress("ele1_bcN",            &ele1_bcN);
  ntu_Data->SetBranchAddress("ele1_bcE",            &ele1_bcE);
  ntu_Data->SetBranchAddress("ele1_bcEta",          &ele1_bcEta);
  ntu_Data->SetBranchAddress("ele1_bcPhi",          &ele1_bcPhi);
  ntu_Data->SetBranchAddress("ele1_bcLocalEta",     &ele1_bcLocalEta);
  ntu_Data->SetBranchAddress("ele1_bcLocalPhi",     &ele1_bcLocalPhi);

  ntu_Data->SetBranchAddress("ele2_isEB",           &ele2_isEB);
  ntu_Data->SetBranchAddress("ele2_isEBEEGap",      &ele2_isEBEEGap);
  ntu_Data->SetBranchAddress("ele2_e3x3",           &ele2_e3x3);
  ntu_Data->SetBranchAddress("ele2_scERaw",         &ele2_scERaw);
  ntu_Data->SetBranchAddress("ele2_scERaw_PUcleaned",&ele2_scERaw_PUcleaned);
  ntu_Data->SetBranchAddress("ele2_scE",            &ele2_scE);
  ntu_Data->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression);
  ntu_Data->SetBranchAddress("ele2_fCorrection",    &ele2_fCorr);
  ntu_Data->SetBranchAddress("ele2_fCorrection_PUcleaned", &ele2_fCorr_PUcleaned);
  ntu_Data->SetBranchAddress("ele2_scCrackCorr",    &ele2_scCrackCorr);
  ntu_Data->SetBranchAddress("ele2_scLocalContCorr",    &ele2_scLocalContCorr);
  ntu_Data->SetBranchAddress("ele2_scEta",          &ele2_scEta);
  ntu_Data->SetBranchAddress("ele2_scPhi",          &ele2_scPhi);
  ntu_Data->SetBranchAddress("ele2_eta",            &ele2_eta);
  ntu_Data->SetBranchAddress("ele2_phi",            &ele2_phi);
  ntu_Data->SetBranchAddress("ele2_scLocalPhi",     &ele2_scLocalPhi);
  ntu_Data->SetBranchAddress("ele2_scLocalEta",     &ele2_scLocalEta); 
  ntu_Data->SetBranchAddress("ele2_charge",         &ele2_charge); 
  ntu_Data->SetBranchAddress("ele2_bcN",            &ele2_bcN);
  ntu_Data->SetBranchAddress("ele2_bcE",            &ele2_bcE);
  ntu_Data->SetBranchAddress("ele2_bcEta",          &ele2_bcEta);
  ntu_Data->SetBranchAddress("ele2_bcPhi",          &ele2_bcPhi);
  ntu_Data->SetBranchAddress("ele2_bcLocalEta",     &ele2_bcLocalEta);
  ntu_Data->SetBranchAddress("ele2_bcLocalPhi",     &ele2_bcLocalPhi);

  // Define histograms
  TH1F* h_mZ_EBEB = new TH1F("h_mZ_EBEB", "",2500,65.,115.);
  TH1F* h_mZ_EEEE = new TH1F("h_mZ_EEEE", "",2500,65.,115.);
  TH1F* h_mZ_EBEE = new TH1F("h_mZ_EBEE", "",2500,65.,115.);

  TH1F* h_mZc_EBEB = new TH1F("h_mZc_EBEB", "",2500,65.,115.);
  TH1F* h_mZc_EEEE = new TH1F("h_mZc_EEEE", "",2500,65.,115.);  
  TH1F* h_mZc_EBEE = new TH1F("h_mZc_EBEE", "",2500,65.,115.);
  
  TProfile *mZ_vs_localEta_EBEB = new TProfile("mZ_vs_localEta_EBEB","mZ_vs_localEta_EBEB",20,0,1,0.,2.);
  mZ_vs_localEta_EBEB->SetMarkerStyle(20);
  mZ_vs_localEta_EBEB->SetMarkerColor(kRed);
  mZ_vs_localEta_EBEB->SetLineColor(kRed);
  mZ_vs_localEta_EBEB-> GetXaxis()->SetTitle("#eta_{SC} (deg)");
  mZ_vs_localEta_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZc_vs_localEta_EBEB = new TProfile("mZc_vs_localEta_EBEB", "",20,0,1,0.,2.);
  mZc_vs_localEta_EBEB->SetMarkerStyle(20);
  mZc_vs_localEta_EBEB->SetMarkerColor(kGreen+2);
  mZc_vs_localEta_EBEB->SetLineColor(kGreen+2);
  mZc_vs_localEta_EBEB-> GetXaxis()->SetTitle("#eta_{SC} (deg)");
  mZc_vs_localEta_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZ_vs_localPhi_EBEB = new TProfile("mZ_vs_localPhi_EBEB", "",20,0,1,0.,2.);
  mZ_vs_localPhi_EBEB->SetMarkerStyle(20);
  mZ_vs_localPhi_EBEB->SetMarkerColor(kRed);
  mZ_vs_localPhi_EBEB->SetLineColor(kRed);
  mZ_vs_localPhi_EBEB-> GetXaxis()->SetTitle("#phi_{SC} (deg)");
  mZ_vs_localPhi_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZc_vs_localPhi_EBEB = new TProfile("mZc_vs_localPhi_EBEB", "",20,0,1,0.,2.);
  mZc_vs_localPhi_EBEB->SetMarkerStyle(20);
  mZc_vs_localPhi_EBEB->SetMarkerColor(kGreen+2);
  mZc_vs_localPhi_EBEB->SetLineColor(kGreen+2);
  mZc_vs_localPhi_EBEB-> GetXaxis()->SetTitle("#phi_{SC} (deg)");
  mZc_vs_localPhi_EBEB-> GetYaxis()->SetTitle("< M_{ee}> /M_{Z}");

  TProfile *mZ_vs_Eta_EBEB = new TProfile("mZ_vs_Eta_EBEB", "",50,-2.5,2.5,0.,2.);
  mZ_vs_Eta_EBEB->SetMarkerStyle(20);
  mZ_vs_Eta_EBEB->SetMarkerColor(kRed);
  mZ_vs_Eta_EBEB->SetLineColor(kRed);
  mZ_vs_Eta_EBEB-> GetXaxis()->SetTitle("#eta_{SC}");
  mZ_vs_Eta_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZc_vs_Eta_EBEB = new TProfile("mZc_vs_Eta_EBEB", "",50,-2.5,2.5,0.,2.);
  mZc_vs_Eta_EBEB->SetMarkerStyle(20);
  mZc_vs_Eta_EBEB->SetMarkerColor(kGreen+2);
  mZc_vs_Eta_EBEB->SetLineColor(kGreen+2);
  mZc_vs_Eta_EBEB-> GetXaxis()->SetTitle("#eta_{SC}");
  mZc_vs_Eta_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZ_vs_Eta_EEEE = new TProfile("mZ_vs_Eta_EEEE", "",50,-2.5,2.5,0.,2.);
  mZ_vs_Eta_EEEE->SetMarkerStyle(20);
  mZ_vs_Eta_EEEE->SetMarkerColor(kRed);
  mZ_vs_Eta_EEEE->SetLineColor(kRed);
  mZ_vs_Eta_EEEE-> GetXaxis()->SetTitle("#eta_{SC}");
  mZ_vs_Eta_EEEE-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZc_vs_Eta_EEEE = new TProfile("mZc_vs_Eta_EEEE", "",50,-2.5,2.5,0.,2.);
  mZc_vs_Eta_EEEE->SetMarkerStyle(20);
  mZc_vs_Eta_EEEE->SetMarkerColor(kGreen+2);
  mZc_vs_Eta_EEEE->SetLineColor(kGreen+2);
  mZc_vs_Eta_EEEE-> GetXaxis()->SetTitle("#eta_{SC}");
  mZc_vs_Eta_EEEE-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TH1F *hspread = new TH1F("hspread","spread",800,0.8,1.2);
  TH1F *hspreadc = new TH1F("hspreadc","spread",800,0.8,1.2);

  
  //---- Loop over entries
  int nEntries = ntu_Data -> GetEntriesFast();
  // nEntries = 20000;
  int ww = 1;
  for(int ientry = 0; ientry < nEntries; ++ientry) {
    
    if( ientry%1000 == 0 ) std::cout << ">>>>> Reading entry : " << ientry << "\r" << std::flush;
    
    ntu_Data -> GetEntry(ientry); 
    
    if( isZ == 0 ) continue;
       
    // -- PU weights
    if (usePUweights) ww = puReweighting->GetWeight(npu);
  
    //-- NPV selection
    if (PV_n <  minNPV || PV_n > maxNPV) continue;

    //--- eta cut 
    if ( ele1_isEBEEGap || ele2_isEBEEGap) continue;
    if ( fabs(ele1_scEta) < etaMin ) continue;
    if ( fabs(ele2_scEta) < etaMin ) continue;
    if ( fabs(ele1_scEta) > etaMax ) continue;
    if ( fabs(ele2_scEta) > etaMax ) continue;

    //-- nBC
    if ( ele1_bcN !=1 ) continue;
    if ( ele2_bcN !=1 ) continue;
    


    //--- remove eta gaps
    if ( IsEtaGap(ele1_scEta) || IsEtaGap(ele2_scEta) ) continue;
       
    //--- remove phi cracks
    float myphi = (ele1_scPhi+3.1415926536)/xtalWidth;
    float modphi = (int)myphi%20;
    if (fabs(modphi-10)<2.) continue;
    float myphi2 = (ele2_scPhi+3.1415926536)/xtalWidth;
    float modphi2 = (int)myphi2%20;
    if (fabs(modphi2-10)<2.) continue;
    
    //--- select on r9
    float ele1_R9 = ele1_e3x3/ele1_scE;
    float ele2_R9 = ele2_e3x3/ele2_scE;

    if ( ele1_R9 < r9min || ele1_R9 > r9max ) continue; 
    if ( ele2_R9 < r9min || ele2_R9 > r9max ) continue; 
 
    float ele1_scEc = 0;
    float ele2_scEc = 0;

    // corrected energy: PU cleaning (dynamic clustering)
    if (usePUcleaning){
      if ( ele1_scERaw_PUcleaned < 0) continue;
      if ( ele2_scERaw_PUcleaned < 0) continue;
      ele1_scEc= ele1_scERaw_PUcleaned*ele1_fCorr_PUcleaned*ele1_scCrackCorr;
      ele2_scEc= ele2_scERaw_PUcleaned*ele2_fCorr_PUcleaned*ele2_scCrackCorr;
    }

    // corrected energy: use regression
    if (useRegression){
      ele1_scEc = ele1_scE_regression;
      ele2_scEc = ele2_scE_regression;
    }
    
    // corrected energy: use single BC correction
    if (useLocalContCorrection){
      ele1_scEc = ele1_scE * ele1_scLocalContCorr;
      ele2_scEc = ele2_scE * ele2_scLocalContCorr;
    }

    // corrected energy: use single BC correction
    if (useLocalCorrectionEta){
      ele1_scEc = ele1_scE * basicClusterCorrection( ele1_scLocalEta , modId(ele1_scEta));
      ele2_scEc = ele2_scE * basicClusterCorrection( ele2_scLocalEta , modId(ele2_scEta));
    }
    
    // // corrected energy: use single BC correction
    if (useSingleBasicClusterCorrectionEta){
      for (int ibc = 0; ibc < ele1_bcN; ibc++){
    	ele1_scEc+= ele1_bcE->at(ibc)*singleBasicClusterCorrection(ele1_bcLocalEta->at(ibc),modId(ele1_scEta));
      }
      ele1_scEc*=ele1_fCorr*ele1_scCrackCorr;
      for (int ibc = 0; ibc < ele2_bcN; ibc++){
    	ele2_scEc+= ele2_bcE->at(ibc)*singleBasicClusterCorrection(ele2_bcLocalEta->at(ibc),modId(ele2_scEta));
      }
      ele2_scEc*=ele2_fCorr*ele2_scCrackCorr;
    }

    mZc = mZ*sqrt(ele1_scEc/ele1_scE)*sqrt(ele2_scEc/ele2_scE);
    
    //--- both in EB
    if( ( ele1_isEB == 1) && (ele2_isEB == 1) )      {
      h_mZ_EBEB      -> Fill(mZ,ww);
      h_mZc_EBEB     -> Fill(mZc,ww);
      mZ_vs_Eta_EBEB -> Fill(ele1_scEta, mZ/zmass);
      mZ_vs_Eta_EBEB -> Fill(ele2_scEta, mZ/zmass);
      mZc_vs_Eta_EBEB-> Fill(ele1_scEta, mZc/zmass);
      mZc_vs_Eta_EBEB-> Fill(ele2_scEta, mZc/zmass);
      
      mZ_vs_localEta_EBEB  -> Fill(ele1_scLocalEta+0.5,mZ/zmass);
      mZc_vs_localEta_EBEB -> Fill(ele1_scLocalEta+0.5,mZc/zmass);
      mZ_vs_localEta_EBEB  -> Fill(ele2_scLocalEta+0.5,mZ/zmass);
      mZc_vs_localEta_EBEB -> Fill(ele2_scLocalEta+0.5,mZc/zmass);
      
      if ( (ele1_scEta*ele1_charge) > 0  ){
	mZ_vs_localPhi_EBEB  -> Fill(ele1_scLocalPhi+0.5,mZ/zmass);
	mZc_vs_localPhi_EBEB -> Fill(ele1_scLocalPhi+0.5,mZc/zmass);
      }
      
      if ( (ele2_scEta*ele2_charge) > 0){
	mZ_vs_localPhi_EBEB  -> Fill(ele2_scLocalPhi+0.5,mZ/zmass);
	mZc_vs_localPhi_EBEB -> Fill(ele2_scLocalPhi+0.5,mZc/zmass);
      }
    }
    
    //--- one in EB, one in EE
    if( ( ele1_isEB && !ele2_isEB) ||  ( !ele1_isEB && ele2_isEB )  )
      //if( !isEB || !isEB2 )// at least one in EE
      {
	h_mZ_EBEE  -> Fill(mZ,ww);
	h_mZc_EBEE -> Fill(mZc,ww);
      }      

    //--- both in EE
    if( ( !ele1_isEB && !ele2_isEB)  ) {
      h_mZ_EEEE  -> Fill(mZ,ww);
      h_mZc_EEEE -> Fill(mZc,ww);
      
      mZ_vs_Eta_EBEB->Fill(ele1_scEta, mZ/zmass);
      mZ_vs_Eta_EBEB->Fill(ele2_scEta, mZ/zmass);
      
      mZc_vs_Eta_EBEB->Fill(ele1_scEta, mZc/zmass);
      mZc_vs_Eta_EBEB->Fill(ele2_scEta, mZc/zmass);
      
    }      
    
  }// end loop over entries
  
  
  for (int i =0 ; i< mZ_vs_localEta_EBEB->GetNbinsX(); i++){
    hspread  ->Fill(mZ_vs_localEta_EBEB->GetBinContent(i+1));
    hspreadc ->Fill(mZc_vs_localEta_EBEB->GetBinContent(i+1));
  }
  
  
  
  // //COMPUTE EFFECTIVE SIGMA FROM HISTOS
  // double mlEBEB = 0, mhEBEB = 0, effsigmaEBEB = 0;
  // double mlcEBEB = 0, mhcEBEB = 0, effsigmacEBEB = 0;
  // double mlEBEE = 0, mhEBEE = 0, effsigmaEBEE = 0;
  // double mlcEBEE = 0, mhcEBEE = 0, effsigmacEBEE = 0;
  // double mlEEEE = 0, mhEEEE = 0, effsigmaEEEE = 0;
  // double mlcEEEE = 0, mhcEEEE = 0, effsigmacEEEE = 0;

  // evalEffectiveSigmaFromHisto(h_mZ_EBEB, mlEBEB, mhEBEB, effsigmaEBEB);
  // evalEffectiveSigmaFromHisto(h_mZc_EBEB, mlcEBEB, mhcEBEB, effsigmacEBEB);
  // evalEffectiveSigmaFromHisto(h_mZ_EBEE, mlEBEE, mhEBEE, effsigmaEBEE);
  // evalEffectiveSigmaFromHisto(h_mZc_EBEE, mlcEBEE, mhcEBEE, effsigmacEBEE);
  // evalEffectiveSigmaFromHisto(h_mZ_EEEE, mlEEEE, mhEEEE, effsigmaEEEE);
  // evalEffectiveSigmaFromHisto(h_mZc_EEEE, mlcEEEE, mhcEEEE, effsigmacEEEE);
  

  // int nre = 25;
  // h_mZ_EBEB  -> Rebin(nre);
  // h_mZc_EBEB -> Rebin(nre);
  // h_mZ_EBEE  -> Rebin(nre);
  // h_mZc_EBEE -> Rebin(nre);
  // h_mZ_EEEE  -> Rebin(nre);
  // h_mZc_EEEE -> Rebin(nre);


  TLegend   *tl = new TLegend(0.60,0.15,0.89,0.35);
  tl-> SetFillColor(0);
  tl-> AddEntry(mZ_vs_localEta_EBEB,"uncorrected","PL");
  tl-> AddEntry(mZc_vs_localEta_EBEB,"corrected","PL");
    
  TCanvas* c1 = new TCanvas("c1","c1",100,100,700,500);
  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  mZ_vs_localEta_EBEB->GetYaxis()->SetRangeUser(0.96,1.02);
  //mZ_vs_localEta_EBEB->Draw("");
  //mZc_vs_localEta_EBEB->Draw("same");
  //tl-> Draw("same");

  TCanvas* c2 = new TCanvas("c2","c2",100,100,700,500);
  c2->cd();
  c2->SetGridx();
  c2->SetGridy();
  mZ_vs_localPhi_EBEB->GetYaxis()->SetRangeUser(0.96,1.02);
  //mZ_vs_localPhi_EBEB->Draw("");
  //mZc_vs_localPhi_EBEB->Draw("same");
  //tl-> Draw("same");

  TCanvas* c3 = new TCanvas("c3","c3",100,100,700,500);
  c3->cd();
  c3->SetGridx();
  c3->SetGridy();
  mZ_vs_Eta_EBEB->GetYaxis()->SetRangeUser(0.96,1.02);
  //mZ_vs_Eta_EBEB->Draw("");
  //mZc_vs_Eta_EBEB->Draw("same");
  //tl-> Draw("same");

  TCanvas* c4 = new TCanvas("c4","c4",100,100,700,500);
  c4->cd();
  c4->SetGridx();
  c4->SetGridy();
  mZ_vs_Eta_EEEE->GetYaxis()->SetRangeUser(0.96,1.02);
  //mZ_vs_Eta_EEEE->Draw("");
  //mZc_vs_Eta_EEEE->Draw("same");
  //tl-> Draw("same");
  /*

  //------------------------
  //---------------- fitting
  
  RooRealVar  mass("mass","M(e^{+}e^{-})", 70.0, 110.0,"GeV/c^{2}");
  mass.setBins(10000) ;


  // Parameters for Crystal Ball Lineshape 
  RooRealVar  dm("#Delta m", "offset", 0.0, -5.0, 5.0,"GeV/c^{2}"); 
  RooRealVar  sigma("#sigma_{CB}","sigmaCB", 1.5,0.5,7.5,"GeV/c^{2}"); 
  RooRealVar  alpha("#alpha","alpha", 1.607,0.6,2.0); 
  RooRealVar  n("n","n", 6., 0.5, 50.0); 
  //alpha.setConstant();
  //n.setConstant();

  // Parameters for Breit-Wigner Distribution
  RooRealVar  MZ("M_{Z}", "M_{Z}", 91.188, 80.0, 100.0,"GeV/c^{2}"); 
  RooRealVar  Gamma("#Gamma", "#Gamma", 2.45, 2.0,3.0,"GeV/c^{2}"); 
  MZ.setConstant();
  Gamma.setConstant();

  // Exponential Background
  RooRealVar  bkgshape("bkgshape", "Backgroung Shape", -0.1,-1.0,0.0, "1/GeV/c^{2}");
  RooRealVar  frac("frac", "Signal Fraction", 1.0,0.0,1.0);
  frac.setConstant();
  bkgshape.setConstant();
  
  // Crystal Ball Lineshape
  RooCBShape     cb("cb", "Crystal Ball Lineshape", mass, dm,sigma, alpha, n);
  // Breit-Wigner Distribution
  RooBreitWigner bw("bw","A Breit-Wigner Distribution",mass,MZ,Gamma);
  
  // Convolution p.d.f. using numeric convolution operator based on Fourier Transforms
  RooFFTConvPdf bwcb("bwcb","convolution", mass, bw, cb);
  
  // Background  p.d.f.
  RooExponential bkg("bkg", "Backgroung Distribution", mass, bkgshape);
  
  
  // --- EB-EB
  TCanvas *cEBEB = new TCanvas("cEBEB","cEBEB",700,700);
  cEBEB->SetLeftMargin(0.15);

  RooDataHist hdataEBEB("hdataEBEB","hdataEBEB",RooArgSet(mass), h_mZ_EBEB);
  RooAddPdf   modelEBEB("modelEBEB", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rEBEB = modelEBEB.fitTo(hdataEBEB, Range(70., 110.), Save(),Verbose(0) );
  RooPlot* plotEBEB = mass.frame(Range(70,110));
  hdataEBEB.plotOn(plotEBEB);
  modelEBEB.plotOn(plotEBEB);
  modelEBEB.plotOn(plotEBEB, LineColor(kRed));
  modelEBEB.paramOn(plotEBEB, Layout(0.16,0.46,0.89));
  modelEBEB.plotOn(plotEBEB, LineColor(kRed));
  modelEBEB.plotOn(plotEBEB, Range(mlEBEB,mhEBEB) , LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3005) , VLines(), MoveToBack() );
  plotEBEB->SetTitleOffset(1.7,"Y");
  plotEBEB->Draw();
  TPaveText *pEBEB = (TPaveText*)cEBEB->FindObject("modelEBEB_paramBox");
  pEBEB->SetTextSize(0.025);
  pEBEB->SetLineColor(kRed);

  //evalEffectiveSigma(&modelEBEB, &mass, mlEBEB, mhEBEB, effsigmaEBEB);


 
  RooDataHist hdatacEBEB("hdatacEBEB","hdatacEBEB",RooArgSet(mass), h_mZc_EBEB);
  RooAddPdf   modelcEBEB("modelcEBEB", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rcEBEB = modelcEBEB.fitTo(hdatacEBEB, Range(70., 110.), Save(),Verbose(0));
  RooPlot* plotcEBEB = mass.frame(Range(70,110));
  hdatacEBEB.plotOn(plotcEBEB);
  modelcEBEB.plotOn(plotcEBEB, LineColor(kGreen+2));
  modelcEBEB.paramOn(plotcEBEB, Layout(0.16,0.46,0.63));
  modelcEBEB.plotOn(plotcEBEB, LineColor(kGreen+2));
  modelcEBEB.plotOn(plotcEBEB, Range(mlcEBEB,mhcEBEB) , LineColor(kGreen+2), DrawOption("F"), FillColor(kGreen+2), FillStyle(3004) , VLines(), MoveToBack() );
  plotcEBEB->SetTitleOffset(1.7,"Y");
  plotcEBEB->Draw("same");
  TPaveText *pcEBEB = (TPaveText*)cEBEB->FindObject("modelcEBEB_paramBox");
  pcEBEB->SetTextSize(0.025);
  pcEBEB->SetLineColor(kGreen+2);
  //evalEffectiveSigma(&modelcEBEB, &mass, mlcEBEB, mhcEBEB, effsigmacEBEB);


  // --- EB-EE
  TCanvas *cEBEE = new TCanvas("cEBEE","cEBEE",700,700);
  cEBEE->SetLeftMargin(0.15);

  RooDataHist hdataEBEE("hdataEBEE","hdataEBEE",RooArgSet(mass), h_mZ_EBEE);
  RooAddPdf   modelEBEE("modelEBEE", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rEBEE = modelEBEE.fitTo(hdataEBEE, Range(70., 110.), Save() ,Verbose(0));
  RooPlot* plotEBEE = mass.frame(Range(70,110));
  hdataEBEE.plotOn(plotEBEE);
  modelEBEE.plotOn(plotEBEE);
  modelEBEE.plotOn(plotEBEE, LineColor(kRed));
  modelEBEE.paramOn(plotEBEE, Layout(0.16,0.46,0.89));
  modelEBEE.plotOn(plotEBEE, LineColor(kRed));
  modelEBEE.plotOn(plotEBEE, Range(mlEBEE,mhEBEE) , LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3005) , VLines(), MoveToBack() );
  plotEBEE->SetTitleOffset(1.7,"Y");
  plotEBEE->Draw();
  TPaveText *pEBEE = (TPaveText*)cEBEE->FindObject("modelEBEE_paramBox");
  pEBEE->SetTextSize(0.025);
  pEBEE->SetLineColor(kRed);
  //evalEffectiveSigma(&modelEBEE, &mass, mlEBEE, mhEBEE, effsigmaEBEE);

  RooDataHist hdatacEBEE("hdatacEBEE","hdatacEBEE",RooArgSet(mass), h_mZc_EBEE);
  RooAddPdf   modelcEBEE("modelcEBEE", "Signal + Background", bwcb, bkg, frac);
  RooFitResult *rcEBEE = modelcEBEE.fitTo(hdatacEBEE, Range(70., 110.), Save(),Verbose(0));
  RooPlot* plotcEBEE = mass.frame(Range(70,110));
  hdatacEBEE.plotOn(plotcEBEE);
  modelcEBEE.plotOn(plotcEBEE, LineColor(kGreen+2));
  modelcEBEE.paramOn(plotcEBEE, Layout(0.16,0.46,0.63));
  modelcEBEE.plotOn(plotcEBEE, LineColor(kGreen+2));
  modelcEBEE.plotOn(plotcEBEE, Range(mlcEBEE,mhcEBEE) , LineColor(kGreen+2), DrawOption("F"), FillColor(kGreen+2), FillStyle(3004) , VLines(), MoveToBack() );  
  plotcEBEE->SetTitleOffset(1.7,"Y");
  plotcEBEE->Draw("same");
  TPaveText *pcEBEE = (TPaveText*)cEBEE->FindObject("modelcEBEE_paramBox");
  pcEBEE->SetTextSize(0.025);
  pcEBEE->SetLineColor(kGreen+2);
  //evalEffectiveSigma(&modelcEBEE, &mass, mlcEBEE, mhcEBEE, effsigmacEBEE);
   
  // --- EE-EE
  TCanvas *cEEEE = new TCanvas("cEEEE","cEEEE",700,700);
  cEEEE->SetLeftMargin(0.15);

  RooDataHist hdataEEEE("hdataEEEE","hdataEEEE",RooArgSet(mass), h_mZ_EEEE);
  RooAddPdf   modelEEEE("modelEEEE", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rEEEE = modelEEEE.fitTo(hdataEEEE, Range(70., 110.), Save() ,Verbose(0));
  RooPlot* plotEEEE = mass.frame(Range(70,110));
  hdataEEEE.plotOn(plotEEEE);
  modelEEEE.plotOn(plotEEEE);
  modelEEEE.plotOn(plotEEEE, LineColor(kRed));
  modelEEEE.paramOn(plotEEEE, Layout(0.16,0.46,0.89));
  modelEEEE.plotOn(plotEEEE, LineColor(kRed));
  modelEEEE.plotOn(plotEEEE, Range(mlEEEE,mhEEEE) , LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3005) , VLines(), MoveToBack() );
  plotEEEE->SetTitleOffset(1.7,"Y");
  plotEEEE->Draw();
  TPaveText *pEEEE = (TPaveText*)cEEEE->FindObject("modelEEEE_paramBox");
  pEEEE->SetTextSize(0.025);
  pEEEE->SetLineColor(kRed);
  //evalEffectiveSigma(&modelEEEE, &mass, mlEEEE, mhEEEE, effsigmaEEEE);

  RooDataHist hdatacEEEE("hdatacEEEE","hdatacEEEE",RooArgSet(mass), h_mZc_EEEE);
  RooAddPdf   modelcEEEE("modelcEEEE", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rcEEEE = modelcEEEE.fitTo(hdatacEEEE, Range(70., 110.), Save(),Verbose(0));
  RooPlot* plotcEEEE = mass.frame(Range(70,110));
  hdatacEEEE.plotOn(plotcEEEE);
  modelcEEEE.plotOn(plotcEEEE, LineColor(kGreen+2));
  modelcEEEE.paramOn(plotcEEEE, Layout(0.16,0.46,0.63));
  modelcEEEE.plotOn(plotcEEEE, LineColor(kGreen+2));
  modelcEEEE.plotOn(plotcEEEE, Range(mlcEEEE,mhcEEEE) , LineColor(kGreen+2), DrawOption("F"), FillColor(kGreen+2), FillStyle(3004) , VLines(), MoveToBack() );
  plotcEEEE->SetTitleOffset(1.7,"Y");
  plotcEEEE->Draw("same");
  TPaveText *pcEEEE = (TPaveText*)cEEEE->FindObject("modelcEEEE_paramBox");
  pcEEEE->SetTextSize(0.025);
  pcEEEE->SetLineColor(kGreen+2);
  //evalEffectiveSigma(&modelcEEEE, &mass, mlcEEEE, mhcEEEE, effsigmacEEEE);



  cout << "EBEB uncorrected --> EffectiveSigma = " << effsigmaEBEB  << " GeV" << "  M_min = " << mlEBEB << "  M_max = " << mhEBEB<< endl;
  cout << "EBEB corrected   --> EffectiveSigma = " << effsigmacEBEB << " GeV" << "  M_min = " << mlcEBEB << "  M_max = " << mhcEBEB<< endl;
  if (effsigmaEBEB > effsigmacEBEB)  cout << "EBEB smearing    --> " << sqrt( effsigmaEBEB*effsigmaEBEB - effsigmacEBEB*effsigmacEBEB ) << endl;
  else cout <<  "EBEB smearing    -->  sigma_corr > sigma_corr !!! " << endl;

  cout << "EBEE uncorrected --> EffectiveSigma = " << effsigmaEBEE  << " GeV" << "  M_min = " << mlEBEE << "  M_max = " << mhEBEE<< endl;
  cout << "EBEE corrected   --> EffectiveSigma = " << effsigmacEBEE << " GeV" << "  M_min = " << mlcEBEE << "  M_max = " << mhcEBEE<< endl;
  if (effsigmaEBEE > effsigmacEBEE)  cout << "EBEE smearing    --> " << sqrt( effsigmaEBEE*effsigmaEBEE - effsigmacEBEE*effsigmacEBEE ) << endl;
  else cout <<  "EBEE smearing    -->  sigma_corr > sigma_corr !!! " << endl;
  cout << "EEEE uncorrected --> EffectiveSigma = " << effsigmaEEEE  << " GeV" << "  M_min = " << mlEEEE << "  M_max = " << mhEEEE<< endl;
  cout << "EEEE corrected   --> EffectiveSigma = " << effsigmacEEEE << " GeV" << "  M_min = " << mlcEEEE << "  M_max = " << mhcEEEE<< endl;
  if (effsigmaEEEE > effsigmacEEEE)  cout << "EEEE smearing    --> " << sqrt( effsigmaEEEE*effsigmaEEEE - effsigmacEEEE*effsigmacEEEE ) << endl;
  else cout <<  "EEEE smearing    -->  sigma_corr > sigma_corr !!! " << endl;

 // wring effectivesigma on canvas 
  char title1[100], title2[100];
  TLatex *latex1;
  TLatex *latex2;

  cEBEB->cd();
  sprintf (title1, "#sigma_{eff} = %.2f GeV",effsigmaEBEB);
  sprintf (title2, "#sigma_{eff} = %.2f GeV",effsigmacEBEB);
  latex1 = new TLatex(0.65,0.8,title1);
  latex1->SetNDC();
  latex1->SetTextFont(42);
  latex1->SetTextSize(0.03);
  latex1->SetTextColor(kRed);
  latex1->Draw("same");
  latex2 = new TLatex(0.65,0.75,title2);
  latex2->SetNDC();
  latex2->SetTextFont(42);
  latex2->SetTextSize(0.03);
  latex2->SetTextColor(kGreen+2);
  latex2->Draw("same");

  cEBEE->cd();
  sprintf (title1, "#sigma_{eff} = %.2f GeV",effsigmaEBEE);
  sprintf (title2, "#sigma_{eff} = %.2f GeV",effsigmacEBEE);
  latex1 = new TLatex(0.65,0.8,title1);
  latex1->SetNDC();
  latex1->SetTextFont(42);
  latex1->SetTextSize(0.03);
  latex1->SetTextColor(kRed);
  latex1->Draw("same");
  latex2 = new TLatex(0.65,0.75,title2);
  latex2->SetNDC();
  latex2->SetTextFont(42);
  latex2->SetTextSize(0.03);
  latex2->SetTextColor(kGreen+2);
  latex2->Draw("same");

  cEEEE->cd();
  sprintf (title1, "#sigma_{eff} = %.2f GeV",effsigmaEEEE);
  sprintf (title2, "#sigma_{eff} = %.2f GeV",effsigmacEEEE);
  latex1 = new TLatex(0.65,0.8,title1);
  latex1->SetNDC();
  latex1->SetTextFont(42);
  latex1->SetTextSize(0.03);
  latex1->SetTextColor(kRed);
  latex1->Draw("same");
  latex2 = new TLatex(0.65,0.75,title2);
  latex2->SetNDC();
  latex2->SetTextFont(42);
  latex2->SetTextSize(0.03);
  latex2->SetTextColor(kGreen+2);
  latex2->Draw("same");
  */
  // save in a file
  TFile *fout = new TFile(outfilename,"recreate");
  h_mZ_EBEB->Write();
  h_mZc_EBEB->Write();
  h_mZ_EBEE->Write();
  h_mZc_EBEE->Write();
  h_mZ_EEEE->Write();
  h_mZc_EEEE->Write();
  mZ_vs_localEta_EBEB->Write();
  mZc_vs_localEta_EBEB->Write();
  mZ_vs_localPhi_EBEB->Write();
  mZc_vs_localPhi_EBEB->Write();
  hspread->Write();
  hspreadc->Write();

  // plotEBEB->Write("plotEBEB");
  // plotcEBEB->Write("plotcEBEB");
  // plotEBEE->Write("plotEBEE");
  // plotcEBEE->Write("plotcEBEE");
  // plotEEEE->Write("plotEEEE");
  // plotcEEEE->Write("plotcEEEE");

  // cEBEB->Write("cEBEB");
  // cEBEB->Write("cEBEE");
  // cEBEE->Write("cEEEE");

  cout << "Closing file..." << endl;
  fout->Close();
  


}

