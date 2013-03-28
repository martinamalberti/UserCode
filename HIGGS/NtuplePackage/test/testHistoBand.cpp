#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"


#include "TLine.h"
#include "TStyle.h"


int main(int argc, char** argv)
{
 
 gROOT->Reset();
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(111);
 
 TF1 gaussian("gaussian","-exp(-0.1*x) + exp(-0.2 * x)",0,50);
 
 TH1F histo("histo","histo",1000,0,50);
 histo.FillRandom("gaussian",100000);
 
 TCanvas cc("cc","cc",400,400);
 
 histo.SetLineColor(kRed);
 histo.Draw();
 
 std::cerr << "===== Get Neyman intervals ====" << std::endl;

 std::vector<double> band = getSigmaBands_FeldmanCousins (histo) ;
 
 
 std::cerr << "=======================" << std::endl;
 std::cerr << " " << band.at(0) << " <<  " << band.at(1) << " << " << band.at(2) << " << " << band.at(3) << " << " << band.at(4) << std::endl;
 std::cerr << "=======================" << std::endl;
 
 TLine* lVertLeft95 = new TLine(band.at(0),0,band.at(0),1000);
 lVertLeft95->SetLineColor(kBlue);
 lVertLeft95->SetLineWidth(2);
 lVertLeft95->SetLineStyle(5);
 
 TLine* lVertLeft68 = new TLine(band.at(1),0,band.at(1),1000);
 lVertLeft68->SetLineColor(kMagenta);
 lVertLeft68->SetLineWidth(2);
 lVertLeft68->SetLineStyle(5);
 
 TLine* lVertMiddle = new TLine(band.at(2),0,band.at(2),1000);
 lVertMiddle->SetLineColor(kGreen);
 lVertMiddle->SetLineWidth(2);
 lVertMiddle->SetLineStyle(5);
 
 TLine* lVertRight68 = new TLine(band.at(3),0,band.at(3),1000);
 lVertRight68->SetLineColor(kMagenta);
 lVertRight68->SetLineWidth(2);
 lVertRight68->SetLineStyle(5);

 TLine* lVertRight95 = new TLine(band.at(4),0,band.at(4),1000);
 lVertRight95->SetLineColor(kBlue);
 lVertRight95->SetLineWidth(2);
 lVertRight95->SetLineStyle(5);
 
 
 lVertLeft95->Draw();
 lVertLeft68->Draw();
 lVertMiddle->Draw();
 lVertRight68->Draw();
 lVertRight95->Draw();
 
 cc.SaveAs("exampleBand.png");
 
 
 
  
  return 0;
}
