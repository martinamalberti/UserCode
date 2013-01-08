/*
c++ `root-config --cflags --glibs`  -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit -lMathMore  -o 2DFitForSpin 2DFitForSpin.cc
*/

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include <fstream>

#include "TMath.h"

#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include <RooMCStudy.h>
#include "RooMsgService.h"
#include "RooCustomizer.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif 

#include "Math/QuantFuncMathCore.h"
#include "Math/DistFunc.h"

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

using namespace RooFit ;


 


//--- function to fix parameters in a pdf -------------------
void  setConstantParameters(RooAbsPdf *pdf, RooRealVar* var ){
  cout << "@@@ Fixing parameters of PDF: " << pdf->getTitle() << endl;
  TIterator* it = pdf->getParameters(*var)->createIterator();
  int npar = (pdf->getParameters(*var))->getSize();
  RooRealVar *tmpVar;
  for (int i = 0; i < npar; i++){
      tmpVar = (RooRealVar*)it->Next();
      tmpVar -> setConstant(true); 
  }
}
//------------------------------------------------------------

//--- median calculator---------------------------------------
float getMedian(vector<float> v){
  float median;
  int s = v.size();
  std::sort(v.begin(),v.end());
  if (s%2 != 0) median = v[int(s/2.)];
  else median = 0.5*(v[int(s/2.)] + v[int(s/2.)+1]  );
  return (median);

}
//------------------------------------------------------------


//--- significance calculator---------------------------------
float getSignificance (vector<float> v1, vector<float> v2 ){

  if ( v1.size()==0 || v2.size()==0 ) return(-1);

  TH1F htemp("htemp","htemp",5000,-500,500);
  for (int i=0; i<v2.size(); i++)
    htemp.Fill(v2[i]);
  
  float m = getMedian(v1);
  int a   = htemp.FindBin(m);  
  float pvalue = htemp.Integral(1,a)/ (htemp.Integral());
  float significance = ROOT::Math::normal_quantile_c(pvalue,1.) ;
  return (significance);
}
//------------------------------------------------------------

//--- significance calculator: use a gaus fit---------------------------------
float getSignificanceGausFit (vector<float> v1, vector<float> v2 ){

  if ( v1.size()==0 || v2.size()==0 ) return(-1);

  TH1F htemp("htemp","htemp",500,-500,500);
  htemp.Reset();
  for (int i=0; i<v2.size(); i++)
    htemp.Fill(v2[i]);
  double a = htemp.GetMean()-10*htemp.GetRMS();
  double b = htemp.GetMean()+10*htemp.GetRMS();
  TF1 gaus("gaus","gaus",a,b);
  htemp.Fit("gaus","Q");
  float m = getMedian(v1);
  double pvalue = gaus.Integral(a,m)/ (gaus.Integral(a,b));
  float significance = ROOT::Math::normal_quantile_c(pvalue,1.) ;
  return (significance);
}
//------------------------------------------------------------


//--- mean calculator ---------------------------------
float getMean (vector<float> v){

  float mean = 0;
  for (int i = 0; i < v.size(); i++){
    mean+=v[i];
  }
  mean/=v.size();
  return (mean);
}
//-----------------------------------------------------


//--- Compute conditional chi-squared statistics ------
float getTCC (vector<float> v){

  float m = getMean(v);
  float tcc = 0.;
  for (int i = 0; i < v.size(); i++){
    tcc+=pow(v[i]-m,2);
  }
  tcc/=m;
  return (tcc);
}
//------------------------------------------------------------


//--- print results for conditional chi-squared statistics ----
void printTCCResults (vector<float> v){
  float tcc = getTCC(v);
  float cdf = ROOT::Math::chisquared_cdf_c(tcc,v.size()-1);
  cout << "TCC = " << tcc << "    CDF = " << cdf << "  prob = " << 1 - cdf<< endl;
  
}
//------------------------------------------------------------

//------------------------------------------------------------
void fillFitResultsHistograms( double nsigtrue, double nbkgtrue, double nsig, double nbkg, double esig, double ebkg, TH1F *hsig, TH1F *hsigpull, TH1F *hbkg, TH1F *hbkgpull, int signalOnly  ){
  
  hsig -> Fill(nsig) ;
  hsigpull -> Fill((nsig-nsigtrue)/esig) ;
  if (!signalOnly){
    hbkg   -> Fill(nbkg) ;
    hbkgpull -> Fill((nbkg - nbkgtrue)/ebkg) ;
  }
}

// ******************* MAIN PROGRAM ********************************************
int main (int argc, char ** argv)
{

  // input arguments 
  string outfilename = "test.root"  ;       // output file name
  float lumi         = 30.; // integrated lumi
  int nToys          = 100; // number of toys
  int signalOnly     = 1; // toys with signal only
  int vbfOnly        = 1; // use only vbf as signal
  float bkgScale     = 1.; // scale bkg by this factor
 
  int opt= 0;
  static struct option long_options[] = {
    {"outfilename", required_argument , 0,  'o' },
    {"lumi",        required_argument , 0,  'l' },
    {"nToys",       required_argument , 0,  'n' },
    {"signalOnly",  required_argument , 0,  's' },
    {"vbfOnly",     required_argument , 0,  'v' },
    {"bkgScale",    required_argument , 0,  'b' },
    {0,             0,                  0,  0   }
  };
  
  int long_index =0;
  while ((opt = getopt_long(argc, argv,"o:l:n:s:v:b:", long_options, &long_index )) != -1) {
    switch (opt) {
    case 'o' : outfilename = optarg;
      break;
    case 'l' : lumi = atof(optarg);
      break;
    case 'n' : nToys= atoi(optarg); 
      break;
    case 's' : signalOnly = atoi(optarg);
      break;
    case 'v' : vbfOnly = atoi(optarg);
      break;
    case 'b' : bkgScale = atof(optarg);
      break;
    }
  }


  cout << "Outfile name    : " << outfilename.c_str() << endl;
  cout << "Integrated lumi : " << lumi  << endl;
  cout << "Number of toys  : " << nToys << endl;
  cout << "Signal only     : " << signalOnly << endl;
  cout << "VBF only        : " << vbfOnly << endl;
  cout << "Scaling bkg  by : " << bkgScale << endl;

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  gROOT->SetStyle ("Plain") ;	
  gStyle->SetOptStat ("mr") ;
  gStyle->SetOptFit (1111) ;
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.1);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.06);
  gStyle->SetStatW(0.3);


  //------------------------------------------------------------
  //--- Retrieve workspace from file to get pdf(mass)
  //TFile *f = new TFile("vbf_mvatag_104fb_v1/CMS-HGG_test.root") ;
  
  TFile *f = new TFile("massfact_spintTest_v2/CMS-HGG_test.root") ; 
  
  RooWorkspace* w = (RooWorkspace*) f->Get("cms_hgg_workspace") ;
  
  RooRealVar* l   = w->var("IntLumi") ;
  float lumifact = lumi/(l->getVal()/1000.);
  
  RooRealVar* m   = w->var("CMS_hgg_mass") ;

  // --- signal
  RooDataHist* h_sig_model_m   = (RooDataHist*)w->data("roohist_sig_vbf_mass_m125_cat4");
  float normvbf = h_sig_model_m->sum(false); 

  RooDataHist* ggh_sig_model_m = (RooDataHist*)w->data("roohist_sig_ggh_mass_m125_cat4");
  if (!vbfOnly) 
    h_sig_model_m->add(*ggh_sig_model_m);
  RooHistPdf  sig_model_m("sig_model_m","sig_model_m",*m,*h_sig_model_m,0);
  
  // --- background
  RooAbsPdf *bkg_model_m   = (RooAbsPdf*)w->pdf("data_pol_model_8TeV_cat4");
  w->var("pdf_data_pol_model_8TeV_cat4_norm")->setConstant() ;
  //setConstantParameters(bkg_model_m, m);

  float normggh = ggh_sig_model_m->sum(false); 
  float normsig = sig_model_m.createIntegral(*m,"all")->getVal(); 
  float normbkg = (w->var("pdf_data_pol_model_8TeV_cat4_norm"))->getVal(); 

  normbkg = normbkg*bkgScale; // bkg in CiC analysis, toght category

  cout << "number of vbf events      = " <<  normvbf << endl ;
  cout << "number of ggh events      = " <<  normggh << endl ;
  cout << "number of signal events      = " <<  normsig << endl ;
  cout << "number of background events  = " <<  normbkg << endl ;
  
  f->Close();
  //-------------------------------------------------------------

  
  //-------------------------------------------------------------
  //--- Retrieve workspace from file to get pdf(d) for the signal
  TFile *fld = new TFile("LDWorkspace.root") ;
  
  RooWorkspace* ldw = (RooWorkspace*) fld->Get("LDWorkspace") ;
  RooRealVar* d     = ldw->var("LD") ;
  //  RooAbsPdf *sig0_model_d = (RooAbsPdf*)ldw->pdf("sum_H0_DS");
  //  RooAbsPdf *sig2_model_d = (RooAbsPdf*)ldw->pdf("sum_H1_DS");
  //  setConstantParameters(sig0_model_d, d);
  //  setConstantParameters(sig2_model_d, d);
  
  RooAbsPdf *vbf_s0_model_d = (RooAbsPdf*)ldw->pdf("sum_H0_DS");
  RooAbsPdf *vbf_s2_model_d = (RooAbsPdf*)ldw->pdf("sum_H1_DS");
  setConstantParameters(vbf_s0_model_d, d);
  setConstantParameters(vbf_s2_model_d, d);
  fld->Close();
  //-------------------------------------------------------------

  
  //-------------------------------------------------------------
  //--- Read tree from file to build background pdf(d)
  TFile *fbkg = new TFile("massfact_spintTest_v2/histograms_CMS-HGG_test_spin.root") ;

  std::vector<TTree*> trees;
  trees.push_back((TTree*)fbkg->Get("Data"));
  trees.push_back((TTree*)fbkg->Get("ggh_m125_8TeV"));

  float weight;
  float MGamGam ;
  float deltaPhiJJ;
  float VBFSpin_Discriminant;
  int   category ;
   
  TH1F *hb   = new TH1F("hb","hb",10,-1,1);  
  TH1F *hggh = new TH1F("hggh","hggh",50,-1,1); 
  hb->GetXaxis()->SetTitle("LD");
  for (int i = 0; i < trees.size(); i++){
    trees[i]->SetBranchAddress("deltaPhiJJ",&deltaPhiJJ);
    trees[i]->SetBranchAddress("MGamGam",&MGamGam);
    trees[i]->SetBranchAddress("weight",&weight);
    trees[i]->SetBranchAddress("VBFSpin_Discriminant",&VBFSpin_Discriminant);
    trees[i]->SetBranchAddress("category",&category);
  }
  
  trees[0]->Project("hb","VBFSpin_Discriminant","category==4 && abs(MGamGam-125)>5");
  trees[1]->Project("hggh","VBFSpin_Discriminant","weight*(category==4)");

  hb ->SetDirectory(0);
  hggh ->SetDirectory(0);
  fbkg->Close();
  //-------------------------------------------------------------


  RooDataHist hggh_model_d("hggh_model_d","hggh_model_d",*d,Import(*hggh)) ;
  
  // fit with double gaussian
  RooRealVar  mean1("mean1","mean of gaussian",0.1,-1,1) ;
  RooRealVar  sigma1("sigma1","width of gaussian",0.1,0.,2.) ;
  RooGaussian gauss1("gauss1","gaussian PDF",*d,mean1,sigma1) ;  
  RooRealVar  mean2("mean2","mean of gaussian",0.4,-1,1) ;
  RooRealVar  sigma2("sigma2","width of gaussian",0.2,0.,2.) ;
  RooGaussian gauss2("gauss2","gaussian PDF",*d,mean2,sigma2) ;  
  RooRealVar  frac("frac","frac",0.5,0,1) ;
  RooAddPdf   ggh_model_d("ggh_model_d","ggh_model_d",RooArgList(gauss1,gauss2),frac) ;
  ggh_model_d.fitTo(hggh_model_d);
    
  RooPlot *ggframe = d->frame();
  hggh_model_d.plotOn(ggframe, RefreshNorm());
  ggh_model_d.plotOn(ggframe);

  RooAddPdf *sig0_model_d ;
  RooAddPdf *sig2_model_d ;
  RooRealVar  relnorm("relnorm","relnorm",0.25,0,1) ;

  if (vbfOnly){
    sig0_model_d = (RooAddPdf*)vbf_s0_model_d;
    sig2_model_d = (RooAddPdf*)vbf_s2_model_d ;
  }
  else{
    sig0_model_d = new RooAddPdf("sig0_model_d","vbf+hgg",RooArgList(*vbf_s0_model_d, ggh_model_d),relnorm) ;
    sig2_model_d = new RooAddPdf("sig2_model_d","vbf+hgg",RooArgList(*vbf_s2_model_d, ggh_model_d),relnorm) ;
    relnorm.setVal(normvbf/normsig);
    relnorm.setConstant(true);
  }
  setConstantParameters(sig0_model_d, d);
  setConstantParameters(sig2_model_d, d);

  //---------------------------------------------------------------
  // --- create 2D pdfs as product pdf(m,d)=pdf(m)x pdf(d)
 
  // --- signal
  cout << " *** Building signal model pdf(m,d)..."<<endl;
  RooProdPdf sig0_model ("sig0_model","sig_model_m*sig0_model_d",RooArgList(sig_model_m, *sig0_model_d)) ;
  RooProdPdf sig2_model ("sig2_model","sig_model_m*sig2_model_d",RooArgList(sig_model_m, *sig2_model_d)) ;
  //RooProdPdf sig0_model ("sig0_model","sig_model_m*sig0_model_d",RooArgList(sig_model_m)) ;
  //RooProdPdf sig2_model ("sig2_model","sig_model_m*sig2_model_d",RooArgList(sig_model_m)) ;

  // --- bkg
  cout << "*** Building background model pdf(m,d)..."<<endl;
  RooDataHist hbkg_model_d("hbkg_model_d","hbkg_model_d",*d,Import(*hb)) ;
  RooHistPdf  bkg_model_d("bkg_model_d","bkg_model_d",*d,hbkg_model_d,0);
  RooProdPdf bkg0_model ("bkg0_model","bkg_model_m*bkg_model_d",RooArgList(*bkg_model_m, bkg_model_d)) ;
  RooProdPdf bkg2_model ("bkg2_model","bkg_model_m*bkg_model_d",RooArgList(*bkg_model_m, bkg_model_d)) ;
  //RooProdPdf bkg0_model ("bkg0_model","bkg_model_m*bkg_model_d",RooArgList(*bkg_model_m)) ;
  //RooProdPdf bkg2_model ("bkg2_model","bkg_model_m*bkg_model_d",RooArgList(*bkg_model_m)) ;

  
  // --- Build signal+background 2D PDF ---
  RooRealVar nsig0("nsig0","#signal events",(normsig*lumifact),0.,50.) ;
  RooRealVar nbkg0("nbkg0","#background events",(normbkg*lumifact),100.,600.) ;
  RooRealVar nsig2("nsig2","#signal events",(normsig*lumifact),0.,50.) ;
  RooRealVar nbkg2("nbkg2","#background events",(normbkg*lumifact),100.,600.) ;

  RooAddPdf *model0 ;
  RooAddPdf *model2 ;

  if (signalOnly){
    normbkg = 0;
    model0 = new RooAddPdf("model0","sig_model+bkg_model",RooArgList(sig0_model),RooArgList(nsig0)) ;
    model2 = new RooAddPdf("model2","sig_model+bkg_model",RooArgList(sig2_model),RooArgList(nsig2)) ;
  }
  else {
    model0 = new RooAddPdf("model0","sig_model+bkg_model",RooArgList(sig0_model,bkg0_model),RooArgList(nsig0,nbkg0));
    model2 = new RooAddPdf("model2","sig_model+bkg_model",RooArgList(sig2_model,bkg2_model),RooArgList(nsig2,nbkg2)) ;
  }

  


  //---- Book histograms
  TH1F *hNsig0_m0     = new TH1F("hNsig0_m0","hNsig0_m0",100,0,100);
  TH1F *hNsig2_m0     = new TH1F("hNsig2_m0","hNsig2_m0",100,0,100);
  TH1F *hpullNsig0_m0 = new TH1F("hpullNsig0_m0","hpullNsig0_m0",50,-5,5);
  TH1F *hpullNsig2_m0 = new TH1F("hpullNsig2_m0","hpullNsig2_m0",50,-5,5);
  TH1F *hNbkg0_m0     = new TH1F("hNbkg0_m0","hNbkg0_m0",1000,0,1000);
  TH1F *hNbkg2_m0     = new TH1F("hNbkg2_m0","hNbkg2_m0",1000,0,1000);
  TH1F *hpullNbkg0_m0 = new TH1F("hpullNbkg0_m0","hpullNbkg0_m0",50,-5,5);
  TH1F *hpullNbkg2_m0 = new TH1F("hpullNbkg2_m0","hpullNbkg2_m0",50,-5,5);
  
  TH1F *hNsig0_m2     = new TH1F("hNsig0_m2","hNsig0_m2",100,0,100);
  TH1F *hNsig2_m2     = new TH1F("hNsig2_m2","hNsig2_m2",100,0,100);
  TH1F *hpullNsig0_m2 = new TH1F("hpullNsig0_m2","hpullNsig0_m2",50,-5,5);
  TH1F *hpullNsig2_m2 = new TH1F("hpullNsig2_m2","hpullNsig2_m2",50,-5,5);
  TH1F *hNbkg0_m2     = new TH1F("hNbkg0_m2","hNbkg0_m2",1000,0,1000);
  TH1F *hNbkg2_m2     = new TH1F("hNbkg2_m2","hNbkg2_m2",1000,0,1000);
  TH1F *hpullNbkg0_m2 = new TH1F("hpullNbkg0_m2","hpullNbkg0_m2",50,-5,5);
  TH1F *hpullNbkg2_m2 = new TH1F("hpullNbkg2_m2","hpullNbkg2_m2",50,-5,5);
  
  TH1F *hLLR0[5];
  TH1F *hLLR2[5];
  string hname0[5]={"hLLR0_nsigcut_0","hLLR0_nsigcut_1","hLLR0_nsigcut_2","hLLR0_nsigcut_3","hLLR0_nsigcut_4"} ;
  string hname2[5]={"hLLR2_nsigcut_0","hLLR2_nsigcut_1","hLLR2_nsigcut_2","hLLR2_nsigcut_3","hLLR2_nsigcut_4"} ;

  for (int i=0; i<5 ;i++){
    hLLR0[i]= new TH1F(hname0[i].c_str(),hname0[i].c_str(),200,-50,50);
    hLLR2[i]= new TH1F(hname2[i].c_str(),hname2[i].c_str(),200,-50,50);
    hLLR0[i]->GetXaxis()->SetTitle("LLR");
    hLLR2[i]->GetXaxis()->SetTitle("LLR");
    hLLR0[i]->SetLineColor(kGreen+3);
    hLLR2[i]->SetLineColor(kBlue);
    hLLR0[i]->SetFillColor(kGreen+3);
    hLLR2[i]->SetFillColor(kBlue);
    hLLR0[i]->SetFillStyle(3004);
    hLLR2[i]->SetFillStyle(3005);
  }

  vector<float> nmin; // limit over number of events
  nmin.push_back(0.0*normsig*lumifact);
  nmin.push_back(0.1*normsig*lumifact);
  nmin.push_back(0.3*normsig*lumifact);
  nmin.push_back(0.5*normsig*lumifact);
  nmin.push_back(0.7*normsig*lumifact);

  std::vector< std::vector<float>* > llr0;
  std::vector< std::vector<float>* > llr2;
  for (int j = 0; j<nmin.size();j++) {
    llr0.push_back(new std::vector<float>);
    llr2.push_back(new std::vector<float>);
  }


  // --- Throw toys to compute LLR 
  float nEventsToy = (normbkg+normsig)*lumifact;  
 
  RooMCStudy* MyToy0 = new RooMCStudy(*model0,RooArgSet(*m,*d),Extended(true),Silence()); 
  //RooMCStudy* MyToy0 = new RooMCStudy(*model0,RooArgSet(*m),Extended(true),Silence()); 
  //MyToy0->generateAndFit(nToys,nEventsToy,true);
  MyToy0->generate(nToys,nEventsToy,true);

  RooMCStudy* MyToy2 = new RooMCStudy(*model2,RooArgSet(*m,*d),Extended(true),Silence()); 
  //RooMCStudy* MyToy2 = new RooMCStudy(*model2,RooArgSet(*m),Extended(true),Silence()); 
  //MyToy2->generateAndFit(nToys,nEventsToy,true);
  MyToy2->generate(nToys,nEventsToy,true);
 
  RooDataSet* toySample;
  RooFitResult* fitResult0;
  RooFitResult* fitResult2;
  RooAbsReal* nll0, *nll2; 
  


  vector<float> s00, s02, s20,s22;
  vector<float> b00, b02, b20,b22;
  vector<float> n00, n02, n20,n22;


  for (unsigned int i = 0; i < nToys; i++) {
    
    cout << "  --->>>>>>>>>>>  Fitting toy " << i << "/" << nToys << endl;

    // -- LLR for hyp0
    toySample = (RooDataSet*)MyToy0->genData(i);
    fitResult0 = model0->fitTo(*toySample, Extended(true), Save(true), Verbose(false), PrintLevel(1));
    fitResult2 = model2->fitTo(*toySample, Extended(true), Save(true), Verbose(false), PrintLevel(1));

    fillFitResultsHistograms(normsig*lumifact, normbkg*lumifact, nsig0.getVal(),nbkg0.getVal(), nsig0.getError(), nbkg0.getError(), hNsig0_m0, hpullNsig0_m0, hNbkg0_m0, hpullNbkg0_m0, signalOnly  );
            
    fillFitResultsHistograms(normsig*lumifact, normbkg*lumifact, nsig2.getVal(),nbkg2.getVal(), nsig2.getError(), nbkg2.getError(), hNsig2_m0, hpullNsig2_m0, hNbkg2_m0, hpullNbkg2_m0, signalOnly  );

    s00.push_back(nsig0.getVal());
    s02.push_back(nsig2.getVal());
    b00.push_back(nbkg0.getVal());
    b02.push_back(nbkg2.getVal());
    n00.push_back(nsig0.getVal()+nbkg0.getVal());
    n02.push_back(nsig2.getVal()+nbkg2.getVal());
    
    // -- check fit quality and fill histos
    if ( fitResult0->covQual() == 3 && fitResult0->status() == 0  && fitResult2->covQual() == 3 && fitResult2->status() == 0) 
      {
	nll0 = model0->createNLL(*toySample, Extended(true), NumCPU(1), Verbose(false)) ;
	nll2 = model2->createNLL(*toySample, Extended(true), NumCPU(1), Verbose(false)) ;
	float llr = 2* (nll0->getVal() - nll2->getVal()) ;
	for (int j = 0; j < nmin.size(); j++) {
	  if (nsig0.getVal() >= nmin[j] && nsig2.getVal() >= nmin[j] ) {
	    hLLR0[j]->Fill( llr );
	    llr0[j]->push_back( llr );
	  }
	}
      }
    
    // -- LLR for hyp2
    toySample = (RooDataSet*)MyToy2->genData(i);
    fitResult0 = model0->fitTo(*toySample, Extended(true), Save(true), Verbose(false), PrintLevel(-1));
    fitResult2 = model2->fitTo(*toySample, Extended(true), Save(true), Verbose(false), PrintLevel(-1));
    
    fillFitResultsHistograms(normsig*lumifact, normbkg*lumifact, nsig0.getVal(),nbkg0.getVal(), nsig0.getError(), nbkg0.getError(), hNsig0_m2, hpullNsig0_m2, hNbkg0_m2, hpullNbkg0_m2, signalOnly  );
    
    fillFitResultsHistograms(normsig*lumifact, normbkg*lumifact, nsig2.getVal(),nbkg2.getVal(), nsig2.getError(), nbkg2.getError(), hNsig2_m2, hpullNsig2_m2, hNbkg2_m2, hpullNbkg2_m2, signalOnly  );

    s20.push_back(nsig0.getVal());
    s22.push_back(nsig2.getVal());
    b20.push_back(nbkg0.getVal());
    b22.push_back(nbkg2.getVal());
    n20.push_back(nsig0.getVal()+nbkg0.getVal());
    n22.push_back(nsig2.getVal()+nbkg2.getVal());
  
    // -- check fit quality and fill histos
    if ( fitResult0->covQual() == 3 && fitResult0->status() == 0  && fitResult2->covQual() == 3 && fitResult2->status() == 0)
      {
	nll0 = model0->createNLL(*toySample, Extended(true), NumCPU(1), Verbose(false)) ;
	nll2 = model2->createNLL(*toySample, Extended(true), NumCPU(1), Verbose(false)) ;
	float llr = 2* (nll0->getVal() - nll2->getVal()) ;
	for (int j = 0; j < nmin.size(); j++) {
	  if (nsig0.getVal() >= nmin[j] && nsig2.getVal() >= nmin[j] ) {
	    hLLR2[j]->Fill( llr );
	    llr2[j]->push_back( llr );
	  }
	}	
      }
  }



  //  // Plot x distribution of data and projection of gaussxy on x = Int(dy) gaussxy(x,y)
  RooPlot* xframe = m->frame(Title("X projection of pdf(m)*pdf(d)")) ;
  ((RooDataSet*)MyToy2->genData(nToys-1))->plotOn(xframe) ;
  model0->plotOn(xframe, LineColor(kRed),VisualizeError(*fitResult0), FillColor(kOrange)) ; 
  model0->plotOn(xframe,VisualizeError(*fitResult0,1,kFALSE),DrawOption("L"),LineWidth(2),LineColor(kRed),LineStyle(2));
  model0->plotOn(xframe, LineColor(kRed)) ; 
  model0->plotOn(xframe, Components(sig0_model), LineColor(kGreen+3));
  model2->plotOn(xframe, Components(sig2_model), LineColor(kBlue));
  //if (!signalOnly)  model0->plotOn(xframe, Components(bkg0_model) , LineColor(kRed), LineStyle(kDashed));
  //if (!signalOnly)  model2->plotOn(xframe, Components(bkg2_model) , LineColor(kBlue), LineStyle(kDashed));

  RooPlot* yframe = d->frame(Title("Y projection of pdf(m)*pdf(d)"), Bins(10),Range(-1,1)) ;
  ((RooDataSet*)MyToy2->genData(nToys-1))->plotOn(yframe) ;
  //  model0->plotOn(yframe, LineColor(kRed),VisualizeError(*fitResult0), FillColor(kOrange)) ; 
  model0->plotOn(yframe, LineColor(kRed)) ; 
  model0->plotOn(yframe, Components(sig0_model), LineColor(kGreen+3));
  model2->plotOn(yframe, Components(sig2_model), LineColor(kBlue));
  //if (!signalOnly)  model0->plotOn(yframe, Components(bkg0_model) , LineColor(kRed), LineStyle(kDashed));
  //if (!signalOnly)  model2->plotOn(yframe, Components(bkg2_model) , LineColor(kBlue), LineStyle(kDashed));

//   RooPlot* nsigpullframe0 =  MyToy0->plotPull(nsig0,Bins(30),FitGauss(kTRUE)) ;
//   RooPlot* nsigpullframe2 =  MyToy2->plotPull(nsig2,Bins(30),FitGauss(kTRUE)) ;
//   RooPlot* nsigframe0 =  MyToy0->plotParam(nsig0,Bins(100),Range(0.,100.)) ;
//   RooPlot* nsigframe2 =  MyToy2->plotParam(nsig2,Bins(100), Range(0.,100.)) ;
//   RooPlot* nbkgpullframe0, * nbkgpullframe2, * nbkgframe0, * nbkgframe2;
//   if (!signalOnly) {
//     nbkgpullframe0 =  MyToy0->plotPull(nbkg0,Bins(30),FitGauss(kTRUE)) ;
//     nbkgpullframe2 =  MyToy2->plotPull(nbkg2,Bins(30),FitGauss(kTRUE)) ;
//     nbkgframe0 =  MyToy0->plotParam(nbkg0,Bins(30)) ;
//     nbkgframe2 =  MyToy2->plotParam(nbkg2,Bins(30)) ;
//   }


  //--- compute significance
  for (int j = 0; j < nmin.size(); j++){
    cout << "sigma/sigma_SM > "<< nmin[j]/(normsig*lumifact)<< "  --> Expected significance = " << getSignificance(*llr0[j],*llr2[j])<< endl;
  }
  
  for (int j = 0; j < nmin.size(); j++){
    cout << "sigma/sigma_SM > "<< nmin[j]/(normsig*lumifact)<< "  --> Expected significance (gaus) = " << getSignificanceGausFit(*llr0[j],*llr2[j])<< endl;
  }


  
//   // test for poisson distribution
//   cout << "*** fit spin0 model on s0 dataset: " << endl;
//   printTCCResults(s00);
//   cout << "*** fit spin2 model on s0 dataset: " << endl;
//   printTCCResults(s02);

//   cout << "*** fit spin0 model on s2 dataset: " << endl;
//   printTCCResults(s20);
//   cout << "*** fit spin2 model on s2 dataset: " << endl;
//   printTCCResults(s22);





  //--- save histograms on file
  
  TFile* fout = new TFile(outfilename.c_str(), "RECREATE");
  xframe->Write();
  yframe->Write();
//   nsigpullframe0->Write();
//   nsigpullframe2->Write();
//   nsigframe0->Write();
//   nsigframe2->Write();
//   if (!signalOnly) {
//    nbkgpullframe0->Write(); 
//    nbkgpullframe2->Write();
//    nbkgframe0->Write();
//    nbkgframe2->Write();
//   }    
  for (int i=0; i<nmin.size(); i++){
    hLLR0[i]->Write();
    hLLR2[i]->Write();
  }
  hpullNsig0_m0 -> Write();
  hpullNsig2_m0 -> Write();
  hNsig0_m0 -> Write();
  hNsig2_m0 -> Write();
  hpullNbkg0_m0 -> Write();
  hpullNbkg2_m0 -> Write();
  hNbkg0_m0 -> Write();
  hNbkg2_m0 -> Write();
  hpullNsig0_m2 -> Write();
  hpullNsig2_m2 -> Write();
  hNsig0_m2 -> Write();
  hNsig2_m2 -> Write();
  hpullNbkg0_m2 -> Write();
  hpullNbkg2_m2 -> Write();
  hNbkg0_m2 -> Write();
  hNbkg2_m2 -> Write();
  hb->Write();
  //ggframe->Write();
  fout->Close();



  
  cout << "bye bye!" << endl;
  



}
