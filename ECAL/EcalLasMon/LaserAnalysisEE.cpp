// compile with:
// g++ `root-config --libs --cflags` -lTMVA LaserAnalysisEE.cpp -o LaserAnalysisEE

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "TDirectory.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TTimeStamp.h"

#include "./LaserDataAnalysis.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

using namespace TMVA;

using namespace std;

int main(int argc, char ** argv)
{
  int selected_fed = atoi(argv[1]);
  
  if ( selected_fed >= 610 &&  selected_fed <= 645){
    cout << ">>>>>  Fed  " << selected_fed << " is not in EE !!!"  << endl;
    cout << "       Exiting ...  "  << endl;
    return 0;
  }
    
  bool saveAllChannels = false;
  bool useRegression   = false;

  int * wl = NULL, nwl = 0;
  TChain * tx = new TChain("x");
  //tx->Add("/data2/EcalLaserMonitoringData/ntuples_2011_158851_178888/ntu_data_001*.root");
  tx->Add("/data2/EcalLaserMonitoringData/ntuples_2011_158851_178888/ntu_data_00175*.root");
  tx->Add("/data2/EcalLaserMonitoringData/ntuples_2011_158851_178888/ntu_data_00176*.root");
  tx->Add("/data2/EcalLaserMonitoringData/ntuples_2011_158851_178888/ntu_data_00177*.root");
  //tx->Add("/tmp/malberti/ntu_data_fed605.root"); // reduced ntuple
  
  init_ttree(tx, &x);
  tx->SetBranchStatus("*",0); //disable all branches
  // only read needed branches
  tx->SetBranchStatus("run",1);
  tx->SetBranchStatus("fed",1);
  tx->SetBranchStatus("seq",1);
  tx->SetBranchStatus("detId",1);
  tx->SetBranchStatus("ix",1);
  tx->SetBranchStatus("iy",1);
  tx->SetBranchStatus("iz",1);
  tx->SetBranchStatus("apdpnAB",1);
  tx->SetBranchStatus("apdpnA",1);
  tx->SetBranchStatus("apdpnB",1);
  tx->SetBranchStatus("time",1);
  tx->SetBranchStatus("eta",1);
  tx->SetBranchStatus("lumi",1);
  tx->SetBranchStatus("field", 1);
  tx->SetBranchStatus("wl",1);
  tx->SetBranchStatus("harness",1);
  tx->SetBranchStatus("elecId",1);
  tx->SetBranchStatus("tmax",1);
  tx->SetBranchStatus("l_ampli",1);
  tx->SetBranchStatus("l_rise_time",1);
  tx->SetBranchStatus("l_fwhm",1);
  tx->SetBranchStatus("l_prepulse",1);

  cout << " Number of entries " << tx->GetEntries() << endl;

  TGraph* g_vptpn_las[101][101][2];
  TGraph* g_vptpn_led[101][101][2];
  TGraph* g_ratio[101][101][2];
  TH1F*   h_ratio[101][101][2];
  
  char gname[101];
  int fed = 0;

  for (int k = 0; k < 2; k++){
    for (int i = 0; i < 101; i++){
      for (int j = 0; j < 101; j++){
	if (saveAllChannels){
	if (k==0) sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEM", i, j);
	else sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEP", i, j);
	g_vptpn_las[i][j][k] = new TGraph();
	g_vptpn_las[i][j][k] -> SetTitle(gname);
	g_vptpn_las[i][j][k] -> SetName(gname);
	g_vptpn_las[i][j][k] -> SetMarkerColor(kBlue);
	g_vptpn_las[i][j][k] -> SetMarkerStyle(20);
	g_vptpn_las[i][j][k] -> SetMarkerSize(0.5);
	
	if (k==0) sprintf(gname,"g_vptpn_led_ix%d_iy%d_EEM", i, j);
	else sprintf(gname,"g_vptpn_led_ix%d_iy%d_EEP", i, j);
	g_vptpn_led[i][j][k] = new TGraph();
	g_vptpn_led[i][j][k] -> SetTitle(gname);
	g_vptpn_led[i][j][k] -> SetName(gname);
	g_vptpn_led[i][j][k] -> SetMarkerColor(kCyan+2);
	g_vptpn_led[i][j][k] -> SetMarkerStyle(20);
	g_vptpn_led[i][j][k] -> SetMarkerSize(0.5);
	
	if (k==0) sprintf(gname,"g_ratio_ix%d_iy%d_EEM", i, j);
	else sprintf(gname,"g_ratio_ix%d_iy%d_EEP", i, j);
	g_ratio[i][j][k] = new TGraph();
	g_ratio[i][j][k] -> SetTitle(gname);
	g_ratio[i][j][k] -> SetName(gname);
	g_ratio[i][j][k] -> SetMarkerColor(kBlack);
	g_ratio[i][j][k] -> SetMarkerStyle(20);
	g_ratio[i][j][k] -> SetMarkerSize(0.5);
	}
	if (k==0) sprintf(gname,"h_ratio_ix%d_iy%d_EEM", i, j);
	else sprintf(gname,"h_ratio_ix%d_iy%d_EEP", i, j);
	h_ratio[i][j][k] = new TH1F(gname,gname,1000,0.9,1.1);
      }
    }
  }

  // Plot results:
  TTimeStamp dateMin(2011, 2,  01, 0, kTRUE, 0); 
  TTimeStamp dateMax(2011, 10, 31, 0, kTRUE, 0); 
    
  // profile averaging on one fed
  TProfile *p_fed_vptpn_las[18];
  TProfile *p_fed_vptpn_led[18];
  TProfile *p_fed_ratio[18];

  // profile averaging on one harness
  TProfile *p_vptpn_las[18][20];
  TProfile *p_vptpn_led[18][20];
  TProfile *p_ratio[18][20];
  TProfile *p_dT[18][20];
  TProfile *p_regressionOutput[18][20];

  // other quantities per harness 
  TProfile *p_ratiopn_las[18][20];
  TProfile *p_ratiopn_led[18][20];
  TProfile *p_tmax_las[18][20];
  TProfile *p_tmax_led[18][20];
  
  TProfile *p_matacqampli_las[18][20];
  TProfile *p_matacqrisetime_las[18][20];
  TProfile *p_matacqfwhm_las[18][20];
  TProfile *p_matacqprepulse_las[18][20];
  TProfile *p_matacqtmax_las[18][20];

  // TH2
  TH2F *h_ratio_vs_matacqfwhm[18][20];
  TH2F *h_ratio_vs_matacqampli[18][20];
  
  //TGraph2D
  TGraph2D *g2[18][20];
  
  for (int ism = 0; ism < 18; ism++){
    if (ism < 9 ) fed = 600 + ism+1;
    if (ism >= 9) fed = 636 + ism+1 ; // fix

    sprintf(gname,"p_vptpn_las_fed%d", fed);
    p_fed_vptpn_las[ism] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
    p_fed_vptpn_las[ism] ->SetLineColor(kBlue);
    p_fed_vptpn_las[ism] ->SetMarkerColor(kBlue);
    p_fed_vptpn_las[ism] ->SetMarkerStyle(20);
    p_fed_vptpn_las[ism] ->SetMarkerSize(0.5);
    
    sprintf(gname,"p_vptpn_led_fed%d", fed);
    p_fed_vptpn_led[ism] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
    p_fed_vptpn_led[ism] ->SetLineColor(kCyan+2);
    p_fed_vptpn_led[ism] ->SetMarkerColor(kCyan+2);
    p_fed_vptpn_led[ism] ->SetMarkerStyle(20);
    p_fed_vptpn_led[ism] ->SetMarkerSize(0.5);
    
    sprintf(gname,"p_ratio_fed%d", fed);
    p_fed_ratio[ism] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
    p_fed_ratio[ism] ->SetLineColor(kBlack);
    p_fed_ratio[ism] ->SetMarkerColor(kBlack);
    p_fed_ratio[ism] ->SetMarkerStyle(20);
    p_fed_ratio[ism] ->SetMarkerSize(0.5);   

    for (int ih = 0; ih < 20; ih++){
      
      sprintf(gname,"p_vptpn_las_fed%d_harness%d", fed, ih+1);
      p_vptpn_las[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_vptpn_las[ism][ih] ->SetLineColor(kBlue);
      p_vptpn_las[ism][ih] ->SetMarkerColor(kBlue);
      p_vptpn_las[ism][ih] ->SetMarkerStyle(20);
      p_vptpn_las[ism][ih] ->SetMarkerSize(0.5);

      sprintf(gname,"p_vptpn_led_fed%d_harness%d", fed, ih+1);
      p_vptpn_led[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_vptpn_led[ism][ih] ->SetLineColor(kCyan+2);
      p_vptpn_led[ism][ih] ->SetMarkerColor(kCyan+2);
      p_vptpn_led[ism][ih] ->SetMarkerStyle(20);
      p_vptpn_led[ism][ih] ->SetMarkerSize(0.5);

      sprintf(gname,"p_ratio_fed%d_harness%d", fed, ih+1);
      p_ratio[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_ratio[ism][ih] ->SetLineColor(kBlack);
      p_ratio[ism][ih] ->SetMarkerColor(kBlack);
      p_ratio[ism][ih] ->SetMarkerStyle(20);
      p_ratio[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_dT_fed%d_harness%d", fed, ih+1);
      p_dT[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, -5000,5000);
      p_dT[ism][ih] ->SetLineColor(kBlack);
      p_dT[ism][ih] ->SetMarkerColor(kBlack);
      p_dT[ism][ih] ->SetMarkerStyle(20);
      p_dT[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_ratiopn_las_fed%d_harness%d", fed, ih+1);
      p_ratiopn_las[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_ratiopn_las[ism][ih] ->SetLineColor(kBlue);
      p_ratiopn_las[ism][ih] ->SetMarkerColor(kBlue);
      p_ratiopn_las[ism][ih] ->SetMarkerStyle(20);
      p_ratiopn_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_ratiopn_led_fed%d_harness%d", fed, ih+1);
      p_ratiopn_led[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_ratiopn_led[ism][ih] ->SetLineColor(kCyan+2);
      p_ratiopn_led[ism][ih] ->SetMarkerColor(kCyan+2);
      p_ratiopn_led[ism][ih] ->SetMarkerStyle(20);
      p_ratiopn_led[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_tmax_las_fed%d_harness%d", fed, ih+1);
      p_tmax_las[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_tmax_las[ism][ih] ->SetLineColor(kBlue);
      p_tmax_las[ism][ih] ->SetMarkerColor(kBlue);
      p_tmax_las[ism][ih] ->SetMarkerStyle(20);
      p_tmax_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_tmax_led_fed%d_harness%d", fed, ih+1);
      p_tmax_led[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_tmax_led[ism][ih] ->SetLineColor(kCyan+2);
      p_tmax_led[ism][ih] ->SetMarkerColor(kCyan+2);
      p_tmax_led[ism][ih] ->SetMarkerStyle(20);
      p_tmax_led[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqampli_las_fed%d_harness%d", fed, ih+1);
      p_matacqampli_las[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_matacqampli_las[ism][ih] ->SetLineColor(kBlue);
      p_matacqampli_las[ism][ih] ->SetMarkerColor(kBlue);
      p_matacqampli_las[ism][ih] ->SetMarkerStyle(20);
      p_matacqampli_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqrisetime_las_fed%d_harness%d", fed, ih+1);
      p_matacqrisetime_las[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_matacqrisetime_las[ism][ih] ->SetLineColor(kBlue);
      p_matacqrisetime_las[ism][ih] ->SetMarkerColor(kBlue);
      p_matacqrisetime_las[ism][ih] ->SetMarkerStyle(20);
      p_matacqrisetime_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqfwhm_las_fed%d_harness%d", fed, ih+1);
      p_matacqfwhm_las[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_matacqfwhm_las[ism][ih] ->SetLineColor(kBlue);
      p_matacqfwhm_las[ism][ih] ->SetMarkerColor(kBlue);
      p_matacqfwhm_las[ism][ih] ->SetMarkerStyle(20);
      p_matacqfwhm_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqprepulse_las_fed%d_harness%d", fed, ih+1);
      p_matacqprepulse_las[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,100);
      p_matacqprepulse_las[ism][ih] ->SetLineColor(kBlue);
      p_matacqprepulse_las[ism][ih] ->SetMarkerColor(kBlue);
      p_matacqprepulse_las[ism][ih] ->SetMarkerStyle(20);
      p_matacqprepulse_las[ism][ih] ->SetMarkerSize(0.5);  

      sprintf(gname,"h_ratio_vs_matacqfwhm_fed%d_harness%d", fed, ih+1);
      h_ratio_vs_matacqfwhm[ism][ih] = new TH2F(gname,gname, 100,0,50, 200, 0.9,1.1);
      h_ratio_vs_matacqfwhm[ism][ih] ->GetXaxis()->SetTitle("fwhm (ns)");
      h_ratio_vs_matacqfwhm[ism][ih] ->GetYaxis()->SetTitle("laser/LED");

      sprintf(gname,"h_ratio_vs_matacqampli_fed%d_harness%d", fed, ih+1);
      h_ratio_vs_matacqampli[ism][ih] = new TH2F(gname,gname,500,0,1000, 200, 0.9,1.1);
      h_ratio_vs_matacqampli[ism][ih] ->GetXaxis()->SetTitle("amplitude (ADC)");
      h_ratio_vs_matacqampli[ism][ih] ->GetYaxis()->SetTitle("laser/LED");

      sprintf(gname,"g2_fed%d_harness%d", fed, ih+1);
      g2[ism][ih]= new TGraph2D();
      g2[ism][ih]->SetTitle(gname);
      g2[ism][ih]->SetName(gname);

      if (useRegression){
	sprintf(gname,"p_regressionOutput_fed%d_harness%d", fed, ih+1);
	p_regressionOutput[ism][ih] = new TProfile(gname,gname,10000,dateMin, dateMax, 0,10);
	p_regressionOutput[ism][ih] ->SetLineColor(kRed);
	p_regressionOutput[ism][ih] ->SetMarkerColor(kRed);
	p_regressionOutput[ism][ih] ->SetMarkerStyle(20);
	p_regressionOutput[ism][ih] ->SetMarkerSize(0.5);   
      }
    }
  }

  TH2F  *harnessMap[2];
  harnessMap[0]  = new TH2F("harnessMapEEM","harnessMapEEM",101,0,101,101,0,101);
  harnessMap[1]  = new TH2F("harnessMapEEP","harnessMapEEP",101,0,101,101,0,101);

  TH2F  *fedMap[2];
  fedMap[0]  = new TH2F("fedMapEEM","fedMapEEM",101,0,101,101,0,101);
  fedMap[1]  = new TH2F("fedMapEEP","fedMapEEP",101,0,101,101,0,101);
  
  // create the Reader object
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  float var1, var2, var3, var4, var5, var6;
  if (useRegression) {
    // create a set of variables and declare them to the reader
    // - the variable names must corresponds in name and type to 
    // those given in the weight file(s) that you use
    //reader->AddVariable( "qmax[0]", &var1 );
    //reader->AddVariable( "tmax[0]", &var2 );
    reader->AddVariable( "l_fwhm[0]", &var3);
    //reader->AddVariable( "l_prepulse[0]", &var4 );
    //reader->AddVariable( "l_width90[0]", &var5 );
    reader->AddVariable( "l_ampli[0]", &var6 );
    // book MVA
    reader->BookMVA( "BDTGmethod","weights/TMVARegression_BDTG.weights.xml" ); 
  }


  TTimeStamp tmin1(2011, 9, 1, 0, kTRUE, 0);
  TTimeStamp tmax1(2011, 9, 30, 0, kTRUE, 0);
  
  int evtlas[18][850] = {0};
  int evtled[18][850] = {0};
  int evt[18][850] = {0};
  int tAve;

  float apdpnlas0[18][850], apdpnled0[18][850], tmaxlas0[18][850], tmaxled0[18][850] ;
  float las_ampli0[18][850], las_risetime0[18][850], las_fwhm0[18][850] , las_prepulse0[18][850];

  float apdpntemp[18][850] = {-999.};
  int tledtemp[18][850] = {0};
  float apdpnlednew;

  float mva0[18][850]={0};

  float ratioTemp[101][101][2];
  double y1[101][101][2];
  double y2[101][101][2];
  double y3[101][101][2];
  
  for (int ientry = 0; ientry < tx->GetEntries(); ientry++){
    tx->GetEntry(ientry);
    
    if (ientry%10000000 == 0 ) cout << "Analyzing entry " << ientry << endl;
    
    // select only EE
    if ( isEB(x.fed)) continue;
    
    int harn       = x.harness;
    int fed        = x.fed;
    int ism        = iSM(fed);
    float apdpnlas = x.apdpnAB[0];
    float apdpnled = x.apdpnAB[1];
    int tlas       = x.time[0];
    int tled       = x.time[1];
    int xtal       = x.elecId;      
         
    int iX         = (x.ix) + 50 ;
    if (x.ix < 0) iX = iX + 1 ;
    int iY         = (x.iy) * (x.iz) - 1;
    int iZ         = x.iz;
    if (iZ<0)   iZ = 0;


    // september 2011
    if ( tlas < tmin1.GetSec() || tlas > tmax1.GetSec() ) continue;
    if ( tled < tmin1.GetSec() || tled > tmax1.GetSec() ) continue;

    // check only one fed
    if ( selected_fed!=-1 && fed != selected_fed ) continue;
    
    // check only one harness
    // if (harn !=6  ) continue;
    // if (fed  !=605) continue;
    
    if (apdpnlas<=0) continue;
    if (apdpnled<=0) continue;
   
    // normalize to first point 
    if (evtlas[ism][xtal] == 0) {
      apdpnlas0[ism][xtal]     = apdpnlas;
      tmaxlas0[ism][xtal]      = x.tmax[0];
      las_ampli0[ism][xtal]    = x.l_ampli[0];
      las_risetime0[ism][xtal] = x.l_rise_time[0];
      las_fwhm0[ism][xtal]     = x.l_fwhm[0];
      las_prepulse0[ism][xtal] = x.l_prepulse[0];
      harnessMap[iZ]->Fill(iX-1,iY, harn);
      fedMap[iZ]->Fill(iX-1,iY, fed);
    }

    if (evtled[ism][xtal] == 0) {
      apdpnled0[ism][xtal]     = apdpnled;
      tmaxled0[ism][xtal]      = x.tmax[1];
    }
    
    // average time between two measurements
    tAve = tlas + (tled - tlas)/ 2;
    if ( tlas > tled ) tAve = tled + (tlas - tled)/ 2;
    
    //extrapolate led measurement to laser time
    if ( evtled[ism][xtal] > 0 && apdpntemp[ism][xtal] > 0 && tledtemp[ism][xtal] > 0) 
      apdpnlednew = ( apdpnled - apdpntemp[ism][xtal] )/(tled - tledtemp[ism][xtal]) * (tlas-tled) + apdpnled;
    else apdpnlednew = apdpnled;
   
    apdpnlas/=apdpnlas0[ism][xtal];
    apdpnled/=apdpnled0[ism][xtal];

    float ratio = apdpnlas/apdpnled;

    // check of led quality
    y3[iX][iY][iZ] = apdpnled;

    // --- profile plots  (average on one FED)
    p_fed_vptpn_las[ism] -> Fill(tlas, apdpnlas);
    p_fed_vptpn_led[ism] -> Fill(tled, apdpnled);
    p_fed_ratio[ism]     -> Fill(tAve, ratioTemp[iX][iY][iZ]);
    
    // --- profile plots  (average on one harness)
    p_vptpn_las[ism][harn-1] -> Fill(tlas, apdpnlas);
    p_vptpn_led[ism][harn-1] -> Fill(tled, apdpnled);
    p_ratio[ism][harn-1]     -> Fill(tAve, ratioTemp[iX][iY][iZ]);
    //p_ratio[ism][harn-1]   -> Fill(tlas, apdpnlas/apdpnlednew);
      
    // single channel graphs
    if (saveAllChannels){
      g_vptpn_las[iX][iY][iZ] -> SetPoint(evtlas[ism][xtal], tlas , apdpnlas);
      g_vptpn_led[iX][iY][iZ] -> SetPoint(evtled[ism][xtal], tled , apdpnled);
      g_ratio[iX][iY][iZ]     -> SetPoint(evtled[ism][xtal], tAve , ratioTemp[iX][iY][iZ]);        
    }
    
    if ( tAve > tmin1.GetSec() && tAve < tmax1.GetSec() ){
      if ( fabs( 2*y2[iX][iY][iZ]-y1[iX][iY][iZ]-y3[iX][iY][iZ]) < 10 * fabs (y3[iX][iY][iZ]-y1[iX][iY][iZ]) ) {
	h_ratio[iX][iY][iZ] -> Fill(ratioTemp[iX][iY][iZ]);
      }
    }

    ratioTemp[iX][iY][iZ] = ratio; 
    y1[iX][iY][iZ] = y2[iX][iY][iZ];
    y2[iX][iY][iZ] = apdpnled;

    //-------------------------------------------------------------------------------------
    // --- other variables per harness 
    if (x.apdpnB[0]>0 ) p_ratiopn_las[ism][harn-1]->Fill(tlas, x.apdpnA[0]/x.apdpnB[0]);
    if (x.apdpnB[1]>0 ) p_ratiopn_led[ism][harn-1]->Fill(tled, x.apdpnA[1]/x.apdpnB[1]);
    
    p_tmax_las[ism][harn-1]->Fill(tlas, x.tmax[0]/tmaxlas0[ism][xtal]);
    p_tmax_led[ism][harn-1]->Fill(tled, x.tmax[1]/tmaxled0[ism][xtal]);
    p_matacqampli_las[ism][harn-1]    -> Fill(tlas, x.l_ampli[0]/las_ampli0[ism][xtal]);
    p_matacqrisetime_las[ism][harn-1] -> Fill(tlas, x.l_rise_time[0]/las_risetime0[ism][xtal]);
    p_matacqfwhm_las[ism][harn-1]     -> Fill(tlas, x.l_fwhm[0]/las_fwhm0[ism][xtal]);
        
    //      if ( tlas > tmin1.GetSec() && tlas < tmax1.GetSec() ){
    //       if ( fabs(ratio-1)<0.1 ){
    // 	h_ratio_vs_matacqfwhm[ism][harn-1] ->Fill( x.l_fwhm[0],ratio );
    //  	h_ratio_vs_matacqampli[ism][harn-1]->Fill( x.l_ampli[0],ratio );
    //  	g2[ism][harn-1]->SetPoint(evt[ism][harn-1],x.l_fwhm[0],x.l_ampli[0], ratio);
    //  	//g2[ism][harn-1]->Fill(x.l_fwhm[0],x.l_ampli[0],ratio);
    // 	evt[ism][harn-1]++;
    
    //       }
    //     }
        

    // ---  MVA regression plots
    if (useRegression){
      //var1 = x.qmax[0];
      //var2 = x.tmax[0];
      var3 = x.l_fwhm[0];
      //var4 = x.l_prepulse[0];
      //var5 = x.l_width90[0]; 
      var6 = x.l_ampli[0]; 
      
      double mva = (reader->EvaluateRegression( "BDTGmethod" ))[0];
      if (evtlas[ism][xtal] == 0){
	mva0[ism][xtal] = mva;
      }
      p_regressionOutput[ism][harn-1]->Fill(tlas, mva/mva0[ism][xtal])    ;
    }
    
    evtlas[ism][xtal]++;
    evtled[ism][xtal]++;
    
    apdpntemp[ism][xtal]= apdpnled;
    tledtemp[ism][xtal] = tled;
    
  }// end loop over entries

  cout << "pippone " << endl;
  
  TH1F *hmean = new TH1F("hmean","hmean",1000,0.5,1.5); 
  TH1F *hrms  = new TH1F("hrms","hrms",4000,0.,0.2); 

  TH2F *hmapMean[2];
  hmapMean[0] = new TH2F("hmapMeanEEM","hmapMeanEEM",101,0,101,101,0,101);
  hmapMean[1] = new TH2F("hmapMeanEEP","hmapMeanEEP",101,0,101,101,0,101);

  TH2F *hmapRMS[2];
  hmapRMS[0] = new TH2F("hmapRmsEEM","hmapRmsEEM",101,0,101,101,0,101);
  hmapRMS[1] = new TH2F("hmapRmsEEP","hmapRmsEEP",101,0,101,101,0,101);

  float mean = 0;
  float rms = 0;

  for (int k = 0; k < 2; k++){
    for (int i = 0; i < 101; i++){    
      for (int j = 0; j < 101; j++){
	if (h_ratio[i][j][k] -> GetEntries() == 0) continue;
	mean = h_ratio[i][j][k] -> GetMean();
	rms  = h_ratio[i][j][k] -> GetRMS();
	hmean->Fill(mean);
	hrms ->Fill(rms);
	hmapMean[k]->Fill(i-1,j,mean);
	hmapRMS[k]->Fill(i-1,j,rms);
      }
    }
  }

  // ****** SAVE GRAPHS ******
  cout << "Saving ..." << endl;

  char fname[100];
  if ( selected_fed!=-1) sprintf(fname,"EElaserAnalysis_fed%d_test2.root",selected_fed);
  else sprintf(fname,"EElaserAnalysis_all.root");
  
  TFile *fout = new TFile(fname,"recreate");
   
  if (saveAllChannels){
    TDirectory *mydir      =   fout -> mkdir("by_channel","by_channel");
    mydir->cd();
    for (int k = 0; k < 2; k++){
      for (int i = 0; i < 101; i++){    
	for (int j = 0; j < 101; j++){
	  if ( g_vptpn_las[i][j][k]-> GetN() > 0) g_vptpn_las[i][j][k]-> Write();
	  if ( g_vptpn_led[i][j][k]-> GetN() > 0) g_vptpn_led[i][j][k]-> Write();
	  if ( g_ratio[i][j][k]    -> GetN() > 0) g_ratio[i][j][k]-> Write();
	  if ( h_ratio[i][j][k]    -> GetEntries() > 0) h_ratio[i][j][k]    -> Write();
	}
      }
    }
  }
  
  fout->cd();
  
  for (int ism = 0; ism < 18; ism++){
    
    if ( p_fed_vptpn_las[ism]-> GetEntries() > 0)  	p_fed_vptpn_las[ism]->Write();
    if ( p_fed_vptpn_led[ism]-> GetEntries() > 0)  	p_fed_vptpn_led[ism]->Write();
    if ( p_fed_ratio[ism]    -> GetEntries() > 0)  	p_fed_ratio[ism]->Write();

    for (int ih = 0; ih < 20; ih++){    
      if ( p_vptpn_las[ism][ih]-> GetEntries() > 0)  	p_vptpn_las[ism][ih]->Write();
      if ( p_vptpn_led[ism][ih]-> GetEntries() > 0)  	p_vptpn_led[ism][ih]->Write();
      if ( p_ratio[ism][ih]    -> GetEntries() > 0)  	p_ratio[ism][ih]->Write();
      if ( p_dT[ism][ih]       -> GetEntries() > 0)  	p_dT[ism][ih]->Write();

      if ( p_ratiopn_las[ism][ih]-> GetEntries() > 0) p_ratiopn_las[ism][ih]->Write();
      if ( p_ratiopn_led[ism][ih]-> GetEntries() > 0) p_ratiopn_led[ism][ih]->Write();

      if ( p_tmax_las[ism][ih]-> GetEntries() > 0) p_tmax_las[ism][ih]->Write();
      if ( p_tmax_led[ism][ih]-> GetEntries() > 0) p_tmax_led[ism][ih]->Write();
      
      if ( p_matacqampli_las[ism][ih]   -> GetEntries() > 0) p_matacqampli_las[ism][ih]->Write();
      if ( p_matacqrisetime_las[ism][ih]-> GetEntries() > 0) p_matacqrisetime_las[ism][ih]->Write();
      if ( p_matacqfwhm_las[ism][ih]    -> GetEntries() > 0) p_matacqfwhm_las[ism][ih]->Write();
      if ( p_matacqprepulse_las[ism][ih]-> GetEntries() > 0) p_matacqprepulse_las[ism][ih]->Write();

      if ( h_ratio_vs_matacqfwhm[ism][ih] -> GetEntries() > 0 ) h_ratio_vs_matacqfwhm[ism][ih]->Write();
      if ( h_ratio_vs_matacqampli[ism][ih]-> GetEntries() > 0 ) h_ratio_vs_matacqampli[ism][ih]->Write();
      if ( g2[ism][ih]                    -> GetN() > 0       ) g2[ism][ih]->Write();
      
      if ( useRegression && p_regressionOutput[ism][ih]-> GetEntries() > 0) p_regressionOutput[ism][ih]->Write();
	 
    }
  }

  hmean->Write();
  hrms->Write();

  hmapMean[0]->Write();
  hmapMean[1]->Write();

  hmapRMS[0]->Write();
  hmapRMS[1]->Write();

  harnessMap[0]->Write();
  harnessMap[1]->Write();
  fedMap[0]->Write();
  fedMap[1]->Write();

  
}

