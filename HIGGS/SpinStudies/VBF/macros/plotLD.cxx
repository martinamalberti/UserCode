#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooDataHist.h"
#include "TCanvas.h"

using namespace RooFit ;

void plotLD(){
  
 
  RooWorkspace* ldw = (RooWorkspace*)_file0->Get("LDWorkspace") ;
  RooRealVar* d     = ldw->var("LD") ;
  RooAbsPdf *sig0   = (RooAbsPdf*)ldw->pdf("sum_H0_DS");
  RooAbsPdf *sig2   = (RooAbsPdf*)ldw->pdf("sum_H1_DS");


  RooDataSet *H0_DS = (RooDataSet*)ldw->data("H0_DS");
  RooDataSet *H2_DS = (RooDataSet*)ldw->data("H1_DS");

  float norm = H0_DS->numEntries()/H2_DS->numEntries();


  RooPlot* yframe = d->frame(Range(-0.5,1.)) ;
  H0_DS->plotOn(yframe,LineColor(kGreen+2),MarkerColor(kGreen+2),RefreshNorm(),Binning(50) );
  sig0->plotOn(yframe, LineColor(kGreen+2));
  sig0->plotOn(yframe, Components("gaus3_H0_DS"), LineStyle(kDashed), LineColor(kGreen+2));
  sig0->plotOn(yframe, Components("gaus4_H0_DS"), LineStyle(kDotted), LineColor(kGreen+2));

//   H2_DS->plotOn(yframe,LineColor(kBlue),MarkerColor(kBlue),Rescale(norm),Binning(50));
//   sig2->plotOn(yframe, LineColor(kBlue),Normalization(norm));
//   sig2->plotOn(yframe, Components("gaus3_H1_DS"), LineStyle(kDashed), LineColor(kBlue));
//   sig2->plotOn(yframe, Components("gaus4_H1_DS"), LineStyle(kDotted), LineColor(kBlue));
  yframe->Draw();



 
}
