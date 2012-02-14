// per compilare: g++ -Wall -o TestZpeakVsNvtx `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit -lMathMore TestZpeakVsNvtx.cpp

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


//*** FIT BW (x) CB **************************************************************************************************************
void FitZpeak(TH1F *h , float minM, float maxM, float &deltaM, float &sigmaCB, float &deltaMerr, float &sigmaCBerr, RooPlot **plot, int kColor)
{
  RooRealVar  mass("mass","M(e^{+}e^{-})", minM, maxM,"GeV/c^{2}");
  mass.setBins(10000) ;

  // Parameters for Crystal Ball Lineshape 
  RooRealVar  dm("#Delta m", "offset", 0.0, -5.0, 5.0,"GeV/c^{2}"); 
  RooRealVar  sigma("#sigma_{CB}","sigmaCB", 1.5,0.5,7.5,"GeV/c^{2}"); 
  RooRealVar  alpha("#alpha","alpha", 1.65,0.6,4.0); 
  RooRealVar  n("n","n", 1.47, 0.5, 50.0); 
  //alpha.setConstant();// from MC
  //n.setConstant();// from MC

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
  
  // Fit
  RooDataHist hdata("hdata","hdata",RooArgSet(mass), h);
  RooAddPdf   model("model", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *r = model.fitTo(hdata, Range(minM, maxM), Save(),Verbose(0) );
  deltaM  = dm.getVal() ;
  sigmaCB = sigma.getVal() ;
  deltaMerr  = dm.getError() ;
  sigmaCBerr = sigma.getError() ;
 
  // plot
  (*plot) = mass.frame(Range(minM,maxM));
  hdata.plotOn(*plot);
  model.plotOn(*plot);
  model.plotOn(*plot, LineColor(kColor));
  model.paramOn(*plot, Layout(0.16,0.46,0.89));
  model.plotOn(*plot, LineColor(kColor));
  //model.plotOn(plot, Range(ml,mh) , LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3005) , VLines(), MoveToBack() );
  (*plot)->SetTitleOffset(1.7,"Y");
  (*plot)->Draw();
}
//*****************************************************************************************************

void SetGraphStyle (TGraphErrors **g, int kColor){
  (*g)-> SetLineColor(kColor);
  (*g)-> SetMarkerColor(kColor);
  (*g)-> SetMarkerStyle(20);
}


//***MOD**************************************************************************************************
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
  int maxNPV = 20; 

  float etaMin = 0.;
  float etaMax = 2.5;
  
  int useR9categorization = false;
  int numberOfHighR9electrons = 0;
  float r9min = 0.0;
  float r9max = 999;
  float r9threshold = 0.94;

  int useNchangeCategorization = true;
  int numberOfChangedElectrons = 2;

  bool useData = false;
  bool useMC   = true;
  bool usePUweights = true;

  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"outputFiles/mcZpeakVsNvtx_nChange2_PU.root");

  //---- PU weights for MC
  TPileupReweighting *puReweighting = new TPileupReweighting("../CommonTools/weights/PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","hweights"); 
  
  //---- NTUPLES 
  TChain *ntu_Data = new TChain("ntu");
  if (useData){
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
    ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
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
  float mZ, mZc, mZr;
  float ele1_EoP,ele2_EoP;
  float ele1_es,ele2_es;
  float ele1_scE,ele2_scE;
  float ele1_scERaw,ele2_scERaw;
  float ele1_scERaw_PUcleaned,ele2_scERaw_PUcleaned;
  float ele1_e3x3, ele2_e3x3;
  float ele1_scE_regression, ele2_scE_regression;
  float ele1_fCorr, ele2_fCorr;
  float ele1_fCorr_PUcleaned, ele2_fCorr_PUcleaned;
  float ele1_scCrackCorr, ele2_scCrackCorr;
  float ele1_scEta,ele2_scEta;
  float ele1_scPhi,ele2_scPhi;
  float ele1_eta,ele2_eta;
  float ele1_phi,ele2_phi;
  float ele1_scLocalEta, ele2_scLocalEta;
  float ele1_scLocalPhi, ele2_scLocalPhi;
  int   ele1_charge , ele2_charge;
  int ele1_bcN, ele2_bcN;
  
  ntu_Data->SetBranchAddress("PV_n", &PV_n);
  ntu_Data->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_Data->SetBranchAddress("isZ",                 &isZ);
  ntu_Data->SetBranchAddress("ele1ele2_scM",        &mZ);
  ntu_Data->SetBranchAddress("ele1_isEB",           &ele1_isEB);
  ntu_Data->SetBranchAddress("ele1_isEBEEGap",      &ele1_isEBEEGap);
  ntu_Data->SetBranchAddress("ele1_e3x3",           &ele1_e3x3);
  ntu_Data->SetBranchAddress("ele1_EOverP",         &ele1_EoP);
  ntu_Data->SetBranchAddress("ele1_es",             &ele1_es);
  ntu_Data->SetBranchAddress("ele1_scERaw",         &ele1_scERaw);
  ntu_Data->SetBranchAddress("ele1_scERaw_PUcleaned",&ele1_scERaw_PUcleaned);
  ntu_Data->SetBranchAddress("ele1_scE",            &ele1_scE);
  ntu_Data->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression);
  ntu_Data->SetBranchAddress("ele1_fCorrection",    &ele1_fCorr);
  ntu_Data->SetBranchAddress("ele1_fCorrection_PUcleaned", &ele1_fCorr_PUcleaned);
  ntu_Data->SetBranchAddress("ele1_scCrackCorr",    &ele1_scCrackCorr);
  ntu_Data->SetBranchAddress("ele1_scEta",          &ele1_scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi",          &ele1_scPhi);
  ntu_Data->SetBranchAddress("ele1_eta",            &ele1_eta);
  ntu_Data->SetBranchAddress("ele1_phi",            &ele1_phi);
  ntu_Data->SetBranchAddress("ele1_scLocalPhi",     &ele1_scLocalPhi);
  ntu_Data->SetBranchAddress("ele1_scLocalEta",     &ele1_scLocalEta); 
  ntu_Data->SetBranchAddress("ele1_charge",         &ele1_charge); 

  ntu_Data->SetBranchAddress("ele2_isEB",           &ele2_isEB);
  ntu_Data->SetBranchAddress("ele2_isEBEEGap",      &ele2_isEBEEGap);
  ntu_Data->SetBranchAddress("ele2_e3x3",           &ele2_e3x3);
  ntu_Data->SetBranchAddress("ele2_EOverP",         &ele2_EoP);
  ntu_Data->SetBranchAddress("ele2_es",             &ele2_es);
  ntu_Data->SetBranchAddress("ele2_scERaw",         &ele2_scERaw);
  ntu_Data->SetBranchAddress("ele2_scERaw_PUcleaned",&ele2_scERaw_PUcleaned);
  ntu_Data->SetBranchAddress("ele2_scE",            &ele2_scE);
  ntu_Data->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression);
  ntu_Data->SetBranchAddress("ele2_fCorrection",    &ele2_fCorr);
  ntu_Data->SetBranchAddress("ele2_fCorrection_PUcleaned", &ele2_fCorr_PUcleaned);
  ntu_Data->SetBranchAddress("ele2_scCrackCorr",    &ele2_scCrackCorr);
  ntu_Data->SetBranchAddress("ele2_scEta",          &ele2_scEta);
  ntu_Data->SetBranchAddress("ele2_scPhi",          &ele2_scPhi);
  ntu_Data->SetBranchAddress("ele2_eta",            &ele2_eta);
  ntu_Data->SetBranchAddress("ele2_phi",            &ele2_phi);
  ntu_Data->SetBranchAddress("ele2_scLocalPhi",     &ele2_scLocalPhi);
  ntu_Data->SetBranchAddress("ele2_scLocalEta",     &ele2_scLocalEta); 
  ntu_Data->SetBranchAddress("ele2_charge",         &ele2_charge); 

  //-- Define histograms
  TH1F* h_mZ_EBEB_all = new TH1F("h_mZ_EBEB_all","h_mZ_EBEB_all",2500,65.,115.);
  TH1F* h_mZc_EBEB_all = new TH1F("h_mZc_EBEB_all","h_mZc_EBEB_all",2500,65.,115.);
  TH1F* h_mZr_EBEB_all = new TH1F("h_mZr_EBEB_all","h_mZr_EBEB_all",2500,65.,115.);;
  TH1F* h_mZ_EBEB[50] ;
  TH1F* h_mZc_EBEB[50] ;
  TH1F* h_mZr_EBEB[50] ;
  TH1F* h_mZ_EEEE_all = new TH1F("h_mZ_EEEE_all","h_mZ_EEEE_all",2500,65.,115.);
  TH1F* h_mZc_EEEE_all = new TH1F("h_mZc_EEEE_all","h_mZc_EEEE_all",2500,65.,115.);
  TH1F* h_mZr_EEEE_all = new TH1F("h_mZr_EEEE_all","h_mZr_EEEE_all",2500,65.,115.);;
  TH1F* h_mZ_EEEE[50] ;
  TH1F* h_mZc_EEEE[50] ;
  TH1F* h_mZr_EEEE[50] ;
  
  char hname[100];
  for (int i = 0; i < 50; i++){
    sprintf(hname, "h_mZ_EBEB_%d",i);
    h_mZ_EBEB[i] = new TH1F(hname,hname,2500,65.,115.);
    sprintf(hname, "h_mZc_EBEB_%d",i);
    h_mZc_EBEB[i] = new TH1F(hname,hname,2500,65.,115.);
    sprintf(hname, "h_mZr_EBEB_%d",i);
    h_mZr_EBEB[i] = new TH1F(hname,hname,2500,65.,115.);
    sprintf(hname, "h_mZ_EEEE_%d",i);
    h_mZ_EEEE[i] = new TH1F(hname,hname,2500,65.,115.);
    sprintf(hname, "h_mZc_EEEE_%d",i);
    h_mZc_EEEE[i] = new TH1F(hname,hname,2500,65.,115.);
    sprintf(hname, "h_mZr_EEEE_%d",i);
    h_mZr_EEEE[i] = new TH1F(hname,hname,2500,65.,115.);
   }

  
  //---- Loop over entries
  int nEntries = ntu_Data -> GetEntriesFast();
  // nEntries = 20000;
  int ww = 1;
  int nChanged = -1;
  int nHighR9  = -1;
  int nHighR9c = -1;
  int nHighR9r = -1;

  for(int ientry = 0; ientry < nEntries; ++ientry) {
    
    if( ientry%1000 == 0 ) std::cout << ">>>>> Reading entry : " << ientry << "\r" << std::flush;
    
    ntu_Data -> GetEntry(ientry); 
    
    if( isZ == 0 ) continue;
       
    // -- PU weights
    if (usePUweights) ww = puReweighting->GetWeight(npu);
  
    //--- eta cut 
    if ( ele1_isEBEEGap || ele2_isEBEEGap) continue;
    if ( fabs(ele1_scEta) < etaMin ) continue;
    if ( fabs(ele2_scEta) < etaMin ) continue;
    if ( fabs(ele1_scEta) > etaMax ) continue;
    if ( fabs(ele2_scEta) > etaMax ) continue;

    //--- remove eta gaps
    if ( IsEtaGap(ele1_scEta) || IsEtaGap(ele2_scEta) ) continue;
       
    //--- remove phi cracks
    //     float myphi = (scPhi+3.1415926536)/xtalWidth;
    //     float modphi = (int)myphi%20;
    //     if (fabs(modphi-10)<2.) continue;
    //     float myphi2 = (scPhi2+3.1415926536)/xtalWidth;
    //     float modphi2 = (int)myphi2%20;
    //     if (fabs(modphi2-10)<2.) continue;
    

    //--- corrected energy: PU cleaning (dynamic clustering)
    if ( ele1_scERaw_PUcleaned < 0) continue;
    if ( ele2_scERaw_PUcleaned < 0) continue;
 
    float ele1_R9 = ele1_e3x3/ele1_scE;
    float ele2_R9 = ele2_e3x3/ele2_scE;
    float ele1_R9r, ele1_R9c;
    float ele2_R9r, ele2_R9c;

    float ele1_scEc = 0;
    float ele2_scEc = 0;
    float ele1_scEr = 0;
    float ele2_scEr = 0;

    // -- dynamic clustering
    ele1_scEc= ele1_scERaw_PUcleaned*ele1_fCorr_PUcleaned*ele1_scCrackCorr;
    ele2_scEc= ele2_scERaw_PUcleaned*ele2_fCorr_PUcleaned*ele2_scCrackCorr;
    if (ele1_isEB == 0) ele1_scEc= (ele1_scERaw_PUcleaned+ele1_es)*ele1_fCorr_PUcleaned*ele1_scCrackCorr;
    if (ele2_isEB == 0) ele2_scEc= (ele2_scERaw_PUcleaned+ele2_es)*ele2_fCorr_PUcleaned*ele2_scCrackCorr;
    ele1_R9c = ele1_e3x3/ele1_scEc;
    ele2_R9c = ele2_e3x3/ele2_scEc;
    
    //--- corrected energy: use regression
    ele1_scEr = ele1_scE_regression;
    ele2_scEr = ele2_scE_regression;
    ele1_R9r = ele1_e3x3/ele1_scEr;
    ele2_R9r = ele2_e3x3/ele2_scEr;
    
    //-- new mZ value
    mZc = mZ*sqrt(ele1_scEc/ele1_scE)*sqrt(ele2_scEc/ele2_scE);
    mZr = mZ*sqrt(ele1_scEr/ele1_scE)*sqrt(ele2_scEr/ele2_scE);
    
    //-- bin vtx
    int bin = (int)PV_n;
    if (PV_n>maxNPV) bin = maxNPV+1;

    //--- Categories
    // number of electrons cleaned by dynamic seed algo
    float ratio1 = (ele1_scERaw/ele1_scERaw_PUcleaned);
    float ratio2 = (ele2_scERaw/ele2_scERaw_PUcleaned);
    if (  ratio1==1 && ratio2==1 ) nChanged=0;
    if ( (ratio1!=1 && ratio2==1) || (ratio1==1 && ratio2!=1) ) nChanged=1;
    if (  ratio1!=1 && ratio2!=1 ) nChanged=2;
    
    // R9 categories
    if ( ele1_R9 > r9threshold  && ele2_R9 > r9threshold  ) nHighR9 = 2;
    if ((ele1_R9 > r9threshold  && ele2_R9 < r9threshold) || (ele1_R9 < r9threshold  && ele2_R9 > r9threshold)  ) nHighR9 = 1;
    if ( ele1_R9 < r9threshold  && ele2_R9 < r9threshold  ) nHighR9 = 0;
    if ( ele1_R9c > r9threshold  && ele2_R9c > r9threshold  ) nHighR9c = 2;
    if ((ele1_R9c > r9threshold  && ele2_R9c < r9threshold) || (ele1_R9c < r9threshold  && ele2_R9c > r9threshold)  ) nHighR9c = 1;
    if ( ele1_R9c < r9threshold  && ele2_R9c < r9threshold  ) nHighR9c = 0;
    if ( ele1_R9r > r9threshold  && ele2_R9r > r9threshold  ) nHighR9r = 2;
    if ((ele1_R9r > r9threshold  && ele2_R9r < r9threshold) || (ele1_R9r < r9threshold  && ele2_R9r > r9threshold)  ) nHighR9r = 1;
    if ( ele1_R9r < r9threshold  && ele2_R9r < r9threshold  ) nHighR9r = 0;

    //--- both in EB
    if( ( ele1_isEB == 1) && (ele2_isEB == 1) )  {
    
      if (useR9categorization && nHighR9 ==  numberOfHighR9electrons){
	h_mZ_EBEB[bin] -> Fill(mZ,ww);
	h_mZ_EBEB_all  -> Fill(mZ,ww);
      }
      
      if (useR9categorization && nHighR9c ==  numberOfHighR9electrons){
	h_mZc_EBEB[bin]-> Fill(mZc,ww);
	h_mZc_EBEB_all -> Fill(mZc,ww);
      }

      if (useR9categorization && nHighR9r ==  numberOfHighR9electrons){
	h_mZr_EBEB[bin]-> Fill(mZr,ww);
	h_mZr_EBEB_all -> Fill(mZr,ww);
      }      

      //...
      if ( useNchangeCategorization && nChanged == numberOfChangedElectrons && nHighR9 ==  numberOfHighR9electrons){
	h_mZ_EBEB[bin] -> Fill(mZ,ww);
	h_mZ_EBEB_all  -> Fill(mZ,ww);
      }
      if ( useNchangeCategorization && nChanged == numberOfChangedElectrons && nHighR9c ==  numberOfHighR9electrons){
	h_mZc_EBEB[bin]-> Fill(mZc,ww);
	h_mZc_EBEB_all -> Fill(mZc,ww);
      }
      if ( useNchangeCategorization && nChanged == numberOfChangedElectrons && nHighR9r ==  numberOfHighR9electrons){
	h_mZr_EBEB[bin]-> Fill(mZr,ww);
	h_mZr_EBEB_all -> Fill(mZr,ww);
      }

    }

    //--- both in EE
    if( ( ele1_isEB == 0) && (ele2_isEB == 0) )  {
    
      if (useR9categorization && nHighR9 ==  numberOfHighR9electrons){
	h_mZ_EEEE[bin] -> Fill(mZ,ww);
	h_mZ_EEEE_all  -> Fill(mZ,ww);
      }
      
      if (useR9categorization && nHighR9c ==  numberOfHighR9electrons){
	h_mZc_EEEE[bin]-> Fill(mZc,ww);
	h_mZc_EEEE_all -> Fill(mZc,ww);
      }

      if (useR9categorization && nHighR9r ==  numberOfHighR9electrons){
	h_mZr_EEEE[bin]-> Fill(mZr,ww);
	h_mZr_EEEE_all -> Fill(mZr,ww);
      }      
   
      //...
      if ( useNchangeCategorization && nChanged == numberOfChangedElectrons){
	h_mZ_EEEE[bin] -> Fill(mZ,ww);
	h_mZ_EEEE_all  -> Fill(mZ,ww);
	h_mZc_EEEE[bin]-> Fill(mZc,ww);
	h_mZc_EEEE_all -> Fill(mZc,ww);
	h_mZr_EEEE[bin]-> Fill(mZr,ww);
	h_mZr_EEEE_all -> Fill(mZr,ww);
      }
      
    }

    
  }// end loop over entries
  
  // FIT
  int rebin = 10;
  float minM = 70.;
  float maxM = 110.;

  // --- EB-EB and EE-EE graphs
  float sigmaCB , deltaM, sigmaCBerr, deltaMerr;

  TGraphErrors *gEBsigma   = new TGraphErrors();
  TGraphErrors *gEBsigmac  = new TGraphErrors();
  TGraphErrors *gEBsigmar  = new TGraphErrors();
  TGraphErrors *gEBdeltam  = new TGraphErrors();
  TGraphErrors *gEBdeltamc = new TGraphErrors();
  TGraphErrors *gEBdeltamr = new TGraphErrors();

  TGraphErrors *gEEsigma   = new TGraphErrors();
  TGraphErrors *gEEsigmac  = new TGraphErrors();
  TGraphErrors *gEEsigmar  = new TGraphErrors();
  TGraphErrors *gEEdeltam  = new TGraphErrors();
  TGraphErrors *gEEdeltamc = new TGraphErrors();
  TGraphErrors *gEEdeltamr = new TGraphErrors();

  // --- all vertices
  RooPlot *plotEBEB_all = new RooPlot();
  RooPlot *plotEBEBc_all= new RooPlot();
  RooPlot *plotEBEBr_all= new RooPlot();  
  RooPlot *plotEEEE_all = new RooPlot();
  RooPlot *plotEEEEc_all= new RooPlot();
  RooPlot *plotEEEEr_all= new RooPlot();  

   
  TCanvas *pippone = new TCanvas("pippone","pippone");
  h_mZ_EBEB_all->Rebin(10);
  h_mZc_EBEB_all->Rebin(10);
  h_mZr_EBEB_all->Rebin(10);
  h_mZ_EEEE_all->Rebin(10);
  h_mZc_EEEE_all->Rebin(10);
  h_mZr_EEEE_all->Rebin(10);
  

  if ( h_mZ_EBEB_all )    FitZpeak(h_mZ_EBEB_all, minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEBEB_all, 2);
  if ( h_mZc_EBEB_all )   FitZpeak(h_mZc_EBEB_all, minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEBEBc_all, 8);
  if ( h_mZr_EBEB_all )   FitZpeak(h_mZr_EBEB_all, minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEBEBr_all, 4);
  if ( h_mZ_EEEE_all )    FitZpeak(h_mZ_EEEE_all, minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEEEE_all, 2);
  if ( h_mZc_EEEE_all )   FitZpeak(h_mZc_EEEE_all, minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEEEEc_all, 8);
  if ( h_mZr_EEEE_all )   FitZpeak(h_mZr_EEEE_all, minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEEEEr_all, 4);

 

  //--- fit vs nvtx
  RooPlot *plotEBEB[100];
  RooPlot *plotEBEBc[100];
  RooPlot *plotEBEBr[100];
  RooPlot *plotEEEE[100];
  RooPlot *plotEEEEc[100];
  RooPlot *plotEEEEr[100];

  for (int i = 1; i < maxNPV+1; i++){
    plotEBEB[i]  = new RooPlot();
    plotEBEBc[i] = new RooPlot();
    plotEBEBr[i] = new RooPlot();
    plotEEEE[i]  = new RooPlot();
    plotEEEEc[i] = new RooPlot();
    plotEEEEr[i] = new RooPlot();
    
    if (i<15){
      h_mZ_EBEB[i]  -> Rebin(rebin);    
      h_mZc_EBEB[i] -> Rebin(rebin);    
      h_mZr_EBEB[i] -> Rebin(rebin);    
      h_mZ_EEEE[i]  -> Rebin(rebin);    
      h_mZc_EEEE[i] -> Rebin(rebin);    
      h_mZr_EEEE[i] -> Rebin(rebin);    
    }
    else {
      h_mZ_EBEB[i]  -> Rebin(2*rebin);    
      h_mZc_EBEB[i] -> Rebin(2*rebin);    
      h_mZr_EBEB[i] -> Rebin(2*rebin);    
      h_mZ_EEEE[i]  -> Rebin(2*rebin);    
      h_mZc_EEEE[i] -> Rebin(2*rebin);    
      h_mZr_EEEE[i] -> Rebin(2*rebin);   
    }
    
    // --- EB-EB uncorrected
    cout << "@@@@@@@@Fitting bin  " << i << " uncorrected " << endl;
    FitZpeak(h_mZ_EBEB[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEBEB[i], 2);
    gEBsigma  -> SetPoint(i, i, sigmaCB);
    gEBsigma  -> SetPointError(i, 0.5, sigmaCBerr);
    gEBdeltam -> SetPoint(i, i, deltaM);
    gEBdeltam -> SetPointError(i, 0.5, deltaMerr);
        
    // --- EB-EB corrected
    cout << "@@@@@@@@Fitting bin  " << i << " dyn. clust " << endl;
    FitZpeak(h_mZc_EBEB[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEBEBc[i], 8);
    gEBsigmac  -> SetPoint(i, i, sigmaCB);
    gEBsigmac  -> SetPointError(i, 0.5, sigmaCBerr);
    gEBdeltamc -> SetPoint(i, i, deltaM);
    gEBdeltamc -> SetPointError(i, 0.5, deltaMerr);

    
    // --- EB-EB corrected with regression
    cout << "@@@@@@@@Fitting bin  " << i << " regression " << endl;
    FitZpeak(h_mZr_EBEB[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEBEBr[i], 4);
    gEBsigmar  -> SetPoint(i, i, sigmaCB);
    gEBsigmar  -> SetPointError(i, 0.5, sigmaCBerr);
    gEBdeltamr -> SetPoint(i, i, deltaM);
    gEBdeltamr -> SetPointError(i, 0.5, deltaMerr);

    // --- EE-EE uncorrected
    FitZpeak(h_mZ_EEEE[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEEEE[i], 2 );
    gEEsigma  -> SetPoint(i, i, sigmaCB);
    gEEsigma  -> SetPointError(i, 0.5, sigmaCBerr);
    gEEdeltam -> SetPoint(i, i, deltaM);
    gEEdeltam -> SetPointError(i, 0.5, deltaMerr);
    
    // --- EE-EE corrected
    FitZpeak(h_mZc_EEEE[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEEEEc[i], 8);
    gEEsigmac  -> SetPoint(i, i, sigmaCB);
    gEEsigmac  -> SetPointError(i, 0.5, sigmaCBerr);
    gEEdeltamc -> SetPoint(i, i, deltaM);
    gEEdeltamc -> SetPointError(i, 0.5, deltaMerr);

    // --- EE-EE corrected with regression
    FitZpeak(h_mZr_EEEE[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, &plotEEEEr[i], 4);
    gEEsigmar  -> SetPoint(i, i, sigmaCB);
    gEEsigmar  -> SetPointError(i, 0.5, sigmaCBerr);
    gEEdeltamr -> SetPoint(i, i, deltaM);
    gEEdeltamr -> SetPointError(i, 0.5, deltaMerr);
    
  }
  

  char plotname[100];

  SetGraphStyle(&gEBsigma, 2);
  SetGraphStyle(&gEBsigmac, 8);
  SetGraphStyle(&gEBsigmar, 4);
  SetGraphStyle(&gEBdeltam, 2);
  SetGraphStyle(&gEBdeltamc, 8);
  SetGraphStyle(&gEBdeltamr, 4);
  SetGraphStyle(&gEEsigma, 2);
  SetGraphStyle(&gEEsigmac, 8);
  SetGraphStyle(&gEEsigmar, 4);
  SetGraphStyle(&gEEdeltam, 2);
  SetGraphStyle(&gEEdeltamc, 8);
  SetGraphStyle(&gEEdeltamr, 4);

  // save in a file
  TFile *fout = new TFile(outfilename,"recreate");

  gEBsigma->Write("gEBsigma"); 
  gEBsigmac->Write("gEBsigmac");
  gEBsigmar->Write("gEBsigmar");

  gEBdeltam->Write("gEBdeltam"); 
  gEBdeltamc->Write("gEBdeltamc");
  gEBdeltamr->Write("gEBdeltamr");

  gEEsigma->Write("gEEsigma"); 
  gEEsigmac->Write("gEEsigmac");
  gEEsigmar->Write("gEEsigmar");

  gEEdeltam->Write("gEEdeltam"); 
  gEEdeltamc->Write("gEEdeltamc");
  gEEdeltamr->Write("gEEdeltamr");


  h_mZ_EBEB_all  -> Write();
  h_mZc_EBEB_all -> Write();
  h_mZr_EBEB_all -> Write();
  h_mZ_EEEE_all  -> Write();
  h_mZc_EEEE_all -> Write();
  h_mZr_EEEE_all -> Write();

  plotEBEB_all ->Write("plotEBEB_all");
  plotEBEBc_all->Write("plotEBEBc_all");
  plotEBEBr_all->Write("plotEBEBr_all");
  plotEEEE_all ->Write("plotEEEE_all");
  plotEEEEc_all->Write("plotEEEEc_all");
  plotEEEEr_all->Write("plotEEEEr_all");

  cout << "++++++++++++ saving all" << endl;
  
  for (int i=1; i<maxNPV+1; i++){

    if ( h_mZ_EBEB[i] -> GetEntries()==0) continue; 
    if ( h_mZ_EEEE[i] -> GetEntries()==0) continue; 

    h_mZ_EBEB[i] -> Write();
    h_mZc_EBEB[i]-> Write();
    h_mZr_EBEB[i]-> Write();
    h_mZ_EEEE[i] -> Write();
    h_mZc_EEEE[i]-> Write();
    h_mZr_EEEE[i]-> Write();

    sprintf(plotname, "plotEBEB_%d",i);
    plotEBEB[i]->Write(plotname);

    sprintf(plotname, "plotEBEBc_%d",i);
    plotEBEBc[i]->Write(plotname);

    sprintf(plotname, "plotEBEBr_%d",i);
    plotEBEBr[i]->Write(plotname);

    sprintf(plotname, "plotEEEE_%d",i);
    plotEEEE[i]->Write(plotname);
    
    sprintf(plotname, "plotEEEEc_%d",i);
    plotEEEEc[i]->Write(plotname);

    sprintf(plotname, "plotEEEEr_%d",i);
    plotEEEEr[i]->Write(plotname);


    
  }

  cout << "Closing file..." << endl;
  fout->Close();
  


}

