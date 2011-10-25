// compile with:
// g++ `root-config --libs --cflags` LaserAnalysisEE.cpp -o LaserAnalysisEE

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
#include "TProfile.h"
#include "TTimeStamp.h"


#include "./LaserDataAnalysis.h"

using namespace std;

int main(int argc, char ** argv)
{
  int selected_fed = atoi(argv[1]);
  
  if ( selected_fed >= 610 &&  selected_fed <= 645){
    cout << ">>>>>  Fed  " << selected_fed << " is not in EE !!!"  << endl;
    cout << "       Exiting ...  "  << endl;
    return 0;
  }
    
  bool saveAllChannels = false ;

  int * wl = NULL, nwl = 0;
  TChain * tx = new TChain("x");
  //tx->Add("/data2/EcalLaserMonitoringData/ntuples_2011_158851_178888/ntu_data_001*.root");
  tx->Add("/tmp/malberti/ntu_data_fed605.root");
  
  init_ttree(tx, &x);
  tx->SetBranchStatus("*",0); //disable all branches
  // only read needed branches
  tx->SetBranchStatus("run",1);
  tx->SetBranchStatus("fed",1);
  tx->SetBranchStatus("seq",1);
  tx->SetBranchStatus("detId",1);
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
  
  // 18 --> number of feds in EE
  // 850 : max number of channels in the fed
  TGraph* gvptpnlas[18][850];
  TGraph* gvptpnled[18][850];
  TGraph* gratio[18][850];
  char gname[100];
  int fed = 0;
  for (int ism = 0; ism < 18; ism++){
    if (ism < 9 ) fed = 600 + ism+1;
    if (ism >= 9) fed = 636 + ism+1 ; // fix
    cout << (ism+1) << "  "  << fed << endl;
    for (int ixtal = 0; ixtal < 850; ixtal++){
      sprintf(gname,"g_vptpnlas_fed%d_xtal%d", fed, ixtal);
      gvptpnlas[ism][ixtal] = new TGraph();
      gvptpnlas[ism][ixtal] -> SetTitle(gname);
      gvptpnlas[ism][ixtal] -> SetName(gname);
      gvptpnlas[ism][ixtal] -> SetMarkerColor(kBlue);
      gvptpnlas[ism][ixtal] -> SetMarkerStyle(20);
      gvptpnlas[ism][ixtal] -> SetMarkerSize(0.5);
      
      sprintf(gname,"g_vptpnled_fed%d_xtal%d", fed, ixtal);
      gvptpnled[ism][ixtal] = new TGraph();
      gvptpnled[ism][ixtal] -> SetTitle(gname);
      gvptpnled[ism][ixtal] -> SetName(gname);
      gvptpnled[ism][ixtal] -> SetMarkerColor(kCyan+2);
      gvptpnled[ism][ixtal] -> SetMarkerStyle(20);
      gvptpnled[ism][ixtal] -> SetMarkerSize(0.5);
      
      sprintf(gname,"g_ratio_fed%d_xtal%d", fed, ixtal);
      gratio[ism][ixtal] = new TGraph();
      gratio[ism][ixtal] -> SetTitle(gname);
      gratio[ism][ixtal] -> SetName(gname);
      gratio[ism][ixtal] -> SetMarkerColor(kBlack);
      gratio[ism][ixtal] -> SetMarkerStyle(20);
      gratio[ism][ixtal] -> SetMarkerSize(0.5);
    }
  }
  
  // Plot results:
  TTimeStamp dateMin(2011, 1,  01, 0, kTRUE, 0); 
  TTimeStamp dateMax(2011, 10, 31, 0, kTRUE, 0); 
  
  TTimeStamp tfake(2011, 9, 28, 0, kTRUE, 0);

  // profile averaging on one harness
  TProfile *p_vptpn_las[18][20];
  TProfile *p_vptpn_led[18][20];
  TProfile *p_ratio[18][20];

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
  
  for (int ism = 0; ism < 18; ism++){
    if (ism < 9 ) fed = 600 + ism+1;
    if (ism >= 9) fed = 636 + ism+1 ; // fix
    for (int ih = 0; ih < 20; ih++){
      
      sprintf(gname,"p_vptpn_las_fed%d_harness%d", fed, ih+1);
      p_vptpn_las[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_vptpn_las[ism][ih] ->SetLineColor(kBlue);
      p_vptpn_las[ism][ih] ->SetMarkerColor(kBlue);
      p_vptpn_las[ism][ih] ->SetMarkerStyle(20);
      p_vptpn_las[ism][ih] ->SetMarkerSize(0.5);

      sprintf(gname,"p_vptpn_led_fed%d_harness%d", fed, ih+1);
      p_vptpn_led[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_vptpn_led[ism][ih] ->SetLineColor(kCyan+2);
      p_vptpn_led[ism][ih] ->SetMarkerColor(kCyan+2);
      p_vptpn_led[ism][ih] ->SetMarkerStyle(20);
      p_vptpn_led[ism][ih] ->SetMarkerSize(0.5);

      sprintf(gname,"p_ratio_fed%d_harness%d", fed, ih+1);
      p_ratio[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_ratio[ism][ih] ->SetLineColor(kBlack);
      p_ratio[ism][ih] ->SetMarkerColor(kBlack);
      p_ratio[ism][ih] ->SetMarkerStyle(20);
      p_ratio[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_ratiopn_las_fed%d_harness%d", fed, ih+1);
      p_ratiopn_las[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_ratiopn_las[ism][ih] ->SetLineColor(kBlue);
      p_ratiopn_las[ism][ih] ->SetMarkerColor(kBlue);
      p_ratiopn_las[ism][ih] ->SetMarkerStyle(20);
      p_ratiopn_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_ratiopn_led_fed%d_harness%d", fed, ih+1);
      p_ratiopn_led[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_ratiopn_led[ism][ih] ->SetLineColor(kCyan+2);
      p_ratiopn_led[ism][ih] ->SetMarkerColor(kCyan+2);
      p_ratiopn_led[ism][ih] ->SetMarkerStyle(20);
      p_ratiopn_led[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_tmax_las_fed%d_harness%d", fed, ih+1);
      p_tmax_las[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_tmax_las[ism][ih] ->SetLineColor(kBlue);
      p_tmax_las[ism][ih] ->SetMarkerColor(kBlue);
      p_tmax_las[ism][ih] ->SetMarkerStyle(20);
      p_tmax_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_tmax_led_fed%d_harness%d", fed, ih+1);
      p_tmax_led[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_tmax_led[ism][ih] ->SetLineColor(kCyan+2);
      p_tmax_led[ism][ih] ->SetMarkerColor(kCyan+2);
      p_tmax_led[ism][ih] ->SetMarkerStyle(20);
      p_tmax_led[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqampli_las_fed%d_harness%d", fed, ih+1);
      p_matacqampli_las[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_matacqampli_las[ism][ih] ->SetLineColor(kBlue);
      p_matacqampli_las[ism][ih] ->SetMarkerColor(kBlue);
      p_matacqampli_las[ism][ih] ->SetMarkerStyle(20);
      p_matacqampli_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqrisetime_las_fed%d_harness%d", fed, ih+1);
      p_matacqrisetime_las[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_matacqrisetime_las[ism][ih] ->SetLineColor(kBlue);
      p_matacqrisetime_las[ism][ih] ->SetMarkerColor(kBlue);
      p_matacqrisetime_las[ism][ih] ->SetMarkerStyle(20);
      p_matacqrisetime_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqfwhm_las_fed%d_harness%d", fed, ih+1);
      p_matacqfwhm_las[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_matacqfwhm_las[ism][ih] ->SetLineColor(kBlue);
      p_matacqfwhm_las[ism][ih] ->SetMarkerColor(kBlue);
      p_matacqfwhm_las[ism][ih] ->SetMarkerStyle(20);
      p_matacqfwhm_las[ism][ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqprepulse_las_fed%d_harness%d", fed, ih+1);
      p_matacqprepulse_las[ism][ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_matacqprepulse_las[ism][ih] ->SetLineColor(kBlue);
      p_matacqprepulse_las[ism][ih] ->SetMarkerColor(kBlue);
      p_matacqprepulse_las[ism][ih] ->SetMarkerStyle(20);
      p_matacqprepulse_las[ism][ih] ->SetMarkerSize(0.5);  

    }
  }
    
  int evtlas[18][850] = {0};
  int evtled[18][850] = {0};
  int tAve;

  float apdpnlas0[18][850], apdpnled0[18][850], tmaxlas0[18][850], tmaxled0[18][850] ;
  float las_ampli0[18][850], las_risetime0[18][850], las_fwhm0[18][850] , las_prepulse0[18][850];

  float apdpntemp[18][850] = {-999.};
  int tledtemp[18][850] = {0};
  float apdpnlednew;

  for (int ientry = 0; ientry < tx->GetEntries(); ientry++){
    tx->GetEntry(ientry);
    
    if (ientry%1000000 == 0 ) cout << "Analyzing entry " << ientry << endl;
    
    // select only EE
    if ( isEB(x.fed)) continue;
    
    int harn    = x.harness;
    int fed     = x.fed;
    int ism     = iSM(fed);
    float apdpnlas = x.apdpnAB[0];
    float apdpnled = x.apdpnAB[1];
    int tlas       = x.time[0];
    int tled       = x.time[1];
    int xtal       = x.elecId;      

       
    // check only one fed
    if ( fed != selected_fed ) continue;
    
    // check only one harness
    if (harn !=6) continue;
 
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
    }

    if (evtled[ism][xtal] == 0) {
      apdpnled0[ism][xtal]     = apdpnled;
      tmaxled0[ism][xtal]      = x.tmax[1];
    }

    // average time between two measurements
    tAve = tlas + (tled - tlas)/ 2;
    if ( tlas > tled ) tAve = tled + (tlas - tled)/ 2;
    
    // extrapolate led measurement to laser time
    //     if ( apdpntemp[ism][xtal] > 0 && tledtemp[ism][xtal] > 0) 
    //       apdpnlednew = (apdpnled - apdpntemp[ism][xtal] )/(tled - tledtemp[ism][xtal]) * (tlas-tled) + apdpnled;
    //     else apdpnlednew = apdpnled;
    
    
    //apply a fake step
    //if ( tlas > int(tfake.GetSec())) apdpnlas+=0.005;


    // --- profile plots  (average on one harness)
    p_vptpn_las[ism][harn-1]->Fill(tlas, apdpnlas/apdpnlas0[ism][xtal]);
    p_vptpn_led[ism][harn-1]->Fill(tled, apdpnled/apdpnled0[ism][xtal]);
    p_ratio[ism][harn-1]   ->Fill(tAve, (apdpnlas/apdpnlas0[ism][xtal])/(apdpnled/apdpnled0[ism][xtal]));
    //p_ratio[ism][harn-1]   ->Fill(tlas, apdpnlas/apdpnlednew);
      
    p_ratiopn_las[ism][harn-1]->Fill(tlas, x.apdpnA[0]/x.apdpnB[0]);
    p_ratiopn_led[ism][harn-1]->Fill(tled, x.apdpnA[1]/x.apdpnB[1]);

    p_tmax_las[ism][harn-1]->Fill(tlas, x.tmax[0]/tmaxlas0[ism][xtal]);
    p_tmax_led[ism][harn-1]->Fill(tled, x.tmax[1]/tmaxled0[ism][xtal]);

    p_matacqampli_las[ism][harn-1]->Fill(tlas, x.l_ampli[0]/las_ampli0[ism][xtal]);

    p_matacqrisetime_las[ism][harn-1]->Fill(tlas, x.l_rise_time[0]/las_risetime0[ism][xtal]);

    p_matacqfwhm_las[ism][harn-1]->Fill(tlas, x.l_fwhm[0]/las_fwhm0[ism][xtal]);

    p_matacqprepulse_las[ism][harn-1]->Fill(tlas, x.l_prepulse[0]/las_prepulse0[ism][xtal]);

    // plot by channel
    if (saveAllChannels){
      gvptpnlas[ism][xtal]->SetPoint(evtlas[ism][xtal], tlas , apdpnlas);
      gvptpnled[ism][xtal]->SetPoint(evtled[ism][xtal], tled , apdpnled);
      gratio[ism][xtal]->SetPoint(evtled[ism][xtal], tAve , apdpnlas/apdpnled);
      //gratio[ism][xtal]->SetPoint(evtled[ism][xtal], tlas , apdpnlas/apdpnlednew);
    }
    
    evtlas[ism][xtal]++;
    evtled[ism][xtal]++;

    apdpntemp[ism][xtal]= apdpnled;
    tledtemp[ism][xtal] = tled;

    
  }// end loop over entries


//   // ratios : boh , da capire 
//   float y1, y2;
//   int t1, t2;

//   float las;
//   int t ;

//   for (int ism = 0; ism < 18; ism++){
//     for (int ih = 0; ih < 20; ih++){
//       if ( p_vptpnlas[ism][ih]-> GetEntries() < 1) continue;  
//       if ( p_vptpnled[ism][ih]-> GetEntries() < 1) continue;
//       for (int ibin = 1; ibin < p_vptpnlas[ism][ih]-> GetNbinsX(); ibin++){
// 	las = p_vptpnlas[ism][ih]-> GetBinContent(ibin);
// 	t   = p_vptpnlas[ism][ih]-> GetBinCenter(ibin);
	
// 	if (ibin > 1 && las > 0) { 

// 	y1 = p_vptpnled[ism][ih]-> GetBinContent(ibin);
// 	y2 = p_vptpnled[ism][ih]-> GetBinContent(ibin-1);
// 	t2 = p_vptpnled[ism][ih]-> GetBinCenter(ibin-1);
	
// 	t1 = p_vptpnled[ism][ih]-> GetBinCenter(ibin);
// 	if (ibin > 1) {

// 	}
	
// 	float apdpnlednew = (y1 - y2 )/(t1 - t2) * (tlas-tled) + apdpnled;
// //     else apdpnlednew = apdpnled;
// 	p_ratio[ism][ih]   -> Fill(tAve, apdpnlas/apdpnled);
//       }
//     }
//   }



  char fname[100];
  sprintf(fname,"testEElaserAnalysis_fed%d.root",selected_fed);

  TFile *fout = new TFile(fname,"recreate");
   
  if (saveAllChannels){
    TDirectory *mydir      =   fout -> mkdir("by_channel","by_channel");
    mydir->cd();
    for (int ism = 0; ism < 18; ism++){
      for (int ixtal = 0; ixtal < 850; ixtal++){
	if ( gvptpnlas[ism][ixtal]-> GetN() > 0) gvptpnlas[ism][ixtal]-> Write();
	if ( gvptpnled[ism][ixtal]-> GetN() > 0) gvptpnled[ism][ixtal]-> Write();
	if ( gratio[ism][ixtal]-> GetN()    > 0) gratio[ism][ixtal]-> Write();
      }
    }
  }
 
  fout->cd();
  for (int ism = 0; ism < 18; ism++){
    for (int ih = 0; ih < 20; ih++){
      if ( p_vptpn_las[ism][ih]-> GetEntries() > 0)  	p_vptpn_las[ism][ih]->Write();
      if ( p_vptpn_led[ism][ih]-> GetEntries() > 0)  	p_vptpn_led[ism][ih]->Write();
      if ( p_ratio[ism][ih]-> GetEntries()    > 0)  	p_ratio[ism][ih]->Write();
      
      if ( p_ratiopn_las[ism][ih]-> GetEntries() > 0) p_ratiopn_las[ism][ih]->Write();
      if ( p_ratiopn_led[ism][ih]-> GetEntries() > 0) p_ratiopn_led[ism][ih]->Write();

      if ( p_tmax_las[ism][ih]-> GetEntries() > 0) p_tmax_las[ism][ih]->Write();
      if ( p_tmax_led[ism][ih]-> GetEntries() > 0) p_tmax_led[ism][ih]->Write();
      
      if ( p_matacqampli_las[ism][ih]-> GetEntries() > 0) p_matacqampli_las[ism][ih]->Write();
      if ( p_matacqrisetime_las[ism][ih]-> GetEntries() > 0) p_matacqrisetime_las[ism][ih]->Write();
      if ( p_matacqfwhm_las[ism][ih]-> GetEntries() > 0) p_matacqfwhm_las[ism][ih]->Write();
      if ( p_matacqprepulse_las[ism][ih]-> GetEntries() > 0) p_matacqprepulse_las[ism][ih]->Write();
      
    }  
  }

  
}

