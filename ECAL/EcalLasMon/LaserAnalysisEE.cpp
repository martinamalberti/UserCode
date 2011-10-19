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
    

  int * wl = NULL, nwl = 0;
  TChain * tx = new TChain("x");
  tx->Add("Data/ntuples_2011_173559_178888/ntu_data_001*.root");
  
  init_ttree(tx, &x);
  tx->SetBranchStatus("*",0); //disable all branches
  // only read needed branches
  tx->SetBranchStatus("run",1);
  tx->SetBranchStatus("fed",1);
  tx->SetBranchStatus("seq",1);
  tx->SetBranchStatus("detId",1);
  tx->SetBranchStatus("apdpnAB",1);
  tx->SetBranchStatus("time",1);
  tx->SetBranchStatus("eta",1);
  tx->SetBranchStatus("lumi",1);
  tx->SetBranchStatus("field", 1);
  tx->SetBranchStatus("wl",1);
  tx->SetBranchStatus("harness",1);
  tx->SetBranchStatus("elecId",1);
  
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
  TTimeStamp dateMin(2011, 8,  01, 0, kTRUE, 0); 
  TTimeStamp dateMax(2011, 10, 31, 0, kTRUE, 0); 
  
  TTimeStamp tfake(2011, 9, 28, 0, kTRUE, 0);

  // profile averaging on one harness
  TProfile *p_vptpnlas[18][20];
  TProfile *p_vptpnled[18][20];
  TProfile *p_ratio[18][20];
  
  for (int ism = 0; ism < 18; ism++){
    if (ism < 9 ) fed = 600 + ism+1;
    if (ism >= 9) fed = 636 + ism+1 ; // fix
    for (int ih = 0; ih < 20; ih++){
      
      sprintf(gname,"p_vptpnlas_fed%d_harness%d", fed, ih+1);
      p_vptpnlas[ism][ih] = new TProfile(gname,gname,2000,dateMin, dateMax, 0,100);
      p_vptpnlas[ism][ih] ->SetLineColor(kBlue);
      p_vptpnlas[ism][ih] ->SetMarkerColor(kBlue);
      p_vptpnlas[ism][ih] ->SetMarkerStyle(20);
      p_vptpnlas[ism][ih] ->SetMarkerSize(0.5);

      sprintf(gname,"p_vptpnled_fed%d_harness%d", fed, ih+1);
      p_vptpnled[ism][ih] = new TProfile(gname,gname,2000,dateMin, dateMax, 0,100);
      p_vptpnled[ism][ih] ->SetLineColor(kCyan+2);
      p_vptpnled[ism][ih] ->SetMarkerColor(kCyan+2);
      p_vptpnled[ism][ih] ->SetMarkerStyle(20);
      p_vptpnled[ism][ih] ->SetMarkerSize(0.5);

      sprintf(gname,"p_ratio_fed%d_harness%d", fed, ih+1);
      p_ratio[ism][ih] = new TProfile(gname,gname,2000,dateMin, dateMax, 0,100);
      p_ratio[ism][ih] ->SetLineColor(kBlack);
      p_ratio[ism][ih] ->SetMarkerColor(kBlack);
      p_ratio[ism][ih] ->SetMarkerStyle(20);
      p_ratio[ism][ih] ->SetMarkerSize(0.5);   
    }
  }
    
  int evtlas[18][850] = {0};
  int evtled[18][850] = {0};
  int tAve;

  float apdpnlas0[18][850], apdpnled0[18][850];
  
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
    //if (harn !=4) continue;
 
    if (apdpnlas<=0) continue;
    if (apdpnled<=0) continue;
   
    // normalize to first point 
    if (evtlas[ism][xtal] == 0) apdpnlas0[ism][xtal] = apdpnlas;
    if (evtled[ism][xtal] == 0) apdpnled0[ism][xtal] = apdpnled;
    
    if (evtlas[ism][xtal] > 0) apdpnlas/=apdpnlas0[ism][xtal] ;
    if (evtled[ism][xtal] > 0) apdpnled/=apdpnled0[ism][xtal] ;

    // average time between two measurements
    tAve = tlas + (tled - tlas)/ 2;
    if ( tlas > tled ) tAve = tled + (tlas - tled)/ 2;
    
    // extrapolate led measurement to laser time
//     if ( apdpntemp[ism][xtal] > 0 && tledtemp[ism][xtal] > 0) 
//       apdpnlednew = (apdpnled - apdpntemp[ism][xtal] )/(tled - tledtemp[ism][xtal]) * (tlas-tled) + apdpnled;
//     else apdpnlednew = apdpnled;
   
    
    if ( tlas > int(tfake.GetSec())) apdpnlas+=0.005;


 
    p_vptpnlas[ism][harn-1]->Fill(tlas, apdpnlas);
    p_vptpnled[ism][harn-1]->Fill(tled, apdpnled);
    p_ratio[ism][harn-1]   ->Fill(tAve, apdpnlas/apdpnled);
    //p_ratio[ism][harn-1]   ->Fill(tlas, apdpnlas/apdpnlednew);
    
    gvptpnlas[ism][xtal]->SetPoint(evtlas[ism][xtal], tlas , apdpnlas);
    gvptpnled[ism][xtal]->SetPoint(evtled[ism][xtal], tled , apdpnled);
    gratio[ism][xtal]->SetPoint(evtled[ism][xtal], tAve , apdpnlas/apdpnled);
    //gratio[ism][xtal]->SetPoint(evtled[ism][xtal], tlas , apdpnlas/apdpnlednew);
 
    evtlas[ism][xtal]++;
    evtled[ism][xtal]++;

    apdpntemp[ism][xtal]= apdpnled;
    tledtemp[ism][xtal] = tled;

    
  }// end loop over entries


//   // ratios
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
  sprintf(fname,"EElaserAnalysis_fed%d.root",selected_fed);

  TFile *fout = new TFile(fname,"recreate");
    
  TDirectory *mydir      =   fout -> mkdir("by_channel","by_channel");
  mydir->cd();
  for (int ism = 0; ism < 18; ism++){
    for (int ixtal = 0; ixtal < 850; ixtal++){
      if ( gvptpnlas[ism][ixtal]-> GetN() > 0) gvptpnlas[ism][ixtal]-> Write();
      if ( gvptpnled[ism][ixtal]-> GetN() > 0) gvptpnled[ism][ixtal]-> Write();
      if ( gratio[ism][ixtal]-> GetN()    > 0) gratio[ism][ixtal]-> Write();
    }
  }

  fout->cd();
  for (int ism = 0; ism < 18; ism++){
    for (int ih = 0; ih < 20; ih++){
      if ( p_vptpnlas[ism][ih]-> GetEntries() > 0)  	p_vptpnlas[ism][ih]->Write();
      if ( p_vptpnled[ism][ih]-> GetEntries() > 0)  	p_vptpnled[ism][ih]->Write();
      if ( p_ratio[ism][ih]-> GetEntries()    > 0)  	p_ratio[ism][ih]->Write();
    }  
  }

  
}

