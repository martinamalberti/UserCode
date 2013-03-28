#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"


#include "TPolyLine.h"
#include "TLine.h"
#include "TStyle.h"
#include "TRandom.h"


/**
return closed contour defined by the two set of 2D points, 
expected to be two non-crossing lines,
*/
TPolyLine makePoly (std::pair<std::vector<double>, std::vector<double> > lineAbove, std::pair<std::vector<double>, std::vector<double> > lineBelow) 
{
 std::vector<double> c_x ;
 std::vector<double> c_y ;
 for (int i = 0 ; i < lineAbove.first.size () ; ++i)
 {
  c_x.push_back (lineAbove.first.at (i)) ;      
  c_y.push_back (lineAbove.second.at (i)) ;          
 }
 for (int i = lineBelow.first.size () - 1 ; i >= 0 ; --i)
 {
  c_x.push_back (lineBelow.first.at (i)) ;      
  c_y.push_back (lineBelow.second.at (i)) ;          
 }
 c_x.push_back (c_x.at (0)) ;      
 c_y.push_back (c_y.at (0)) ; 
 TPolyLine pl (c_x.size (), &c_x.at (0), &c_y.at (0)) ;
 
 return pl ; 
}


int main(int argc, char** argv)
{
 
 gROOT->SetStyle ("Plain") ;	
 gStyle->SetOptStat ("mr") ;
 gStyle->SetOptFit (1111) ;
 gStyle->SetStatFont (42) ;
 gStyle->SetStatFontSize (0.1) ;
 gStyle->SetStatTextColor (1) ;
 gStyle->SetStatFormat ("6.4g") ;
 gStyle->SetStatBorderSize (1) ;
 gStyle->SetStatH (0.06) ;
 gStyle->SetStatW (0.3) ;
 gStyle->SetPalette (8) ;
 
 
 TRandom generator ;
 TH2F * h2_cava = new TH2F ("h2_cava","cava example", 100, 0, 100, 100, 0, 100) ;
 for (int i = 0 ; i < 1000000 ; ++i)
 {
  double y = generator.Uniform(0,100);
  double x = y ;
  x *= 1 + 0.2 * generator.Gaus () ;
  //          y *= 5 ;
  h2_cava->Fill (x,y) ;
 }
 
 
 
 h2_cava->FitSlicesY (0, 0, -1, 0, "QNR") ;
 TH1F* h2_cava_1 = (TH1F*) gDirectory->Get("h2_cava_1");
 h2_cava_1->SetMarkerSize(0.1);
 h2_cava_1->SetMarkerStyle(20);
 h2_cava_1->SetMarkerColor(kRed);
 
 
 //PG build bands with getBand and same Tails along x
 //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- --
 
 std::vector<std::vector<double> > out_sameTails_x = getBand_integrY (*h2_cava, getLimit_sameTails ()) ;
 
 std::pair<std::vector<double>, std::vector<double> > x_lineAbove_sameTails ;
 x_lineAbove_sameTails.first = out_sameTails_x.at (0) ;
 x_lineAbove_sameTails.second = out_sameTails_x.at (2) ;
 std::pair<std::vector<double>, std::vector<double> > x_lineBelow_sameTails ;
 x_lineBelow_sameTails.first = out_sameTails_x.at (0) ;
 x_lineBelow_sameTails.second = out_sameTails_x.at (1) ;
 
 TPolyLine x_cont_sameTails = makePoly (x_lineAbove_sameTails, x_lineBelow_sameTails) ;
 x_cont_sameTails.SetLineWidth (2.5) ;
 x_cont_sameTails.SetFillColor (kGreen) ;
 x_cont_sameTails.SetLineColor (kGreen) ;
 
 
 
 //PG build bands with getBand and Neyman intervals along x
 //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---
 
 std::vector<std::vector<double> > out_Neyman_x = getBand_integrY (*h2_cava, getLimit_FC ()) ;
 
 std::pair<std::vector<double>, std::vector<double> > x_lineAbove_Neyman ;
 x_lineAbove_Neyman.first = out_Neyman_x.at (0) ;
 x_lineAbove_Neyman.second = out_Neyman_x.at (2) ;
 std::pair<std::vector<double>, std::vector<double> > x_lineBelow_Neyman ;
 x_lineBelow_Neyman.first = out_Neyman_x.at (0) ;
 x_lineBelow_Neyman.second = out_Neyman_x.at (1) ;
 
 TPolyLine x_cont_Neyman = makePoly (x_lineAbove_Neyman, x_lineBelow_Neyman) ;
 x_cont_Neyman.SetLineWidth (2.5) ;
 x_cont_Neyman.SetFillColor (kBlue) ;
 x_cont_Neyman.SetLineColor (kBlue) ;
 
 
 TCanvas* cc = new TCanvas("cc","cc",800,800);
 
 h2_cava->Draw ("COLZ") ;
 h2_cava_1->Draw("Esame");
 x_cont_sameTails.Draw () ;
 x_cont_Neyman.Draw () ;
 cc->Print ("example2D.png") ;
 
 
 return 0;
}
