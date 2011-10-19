// compile with:
// g++ `root-config --libs --cflags` LaserAnalysisEB.cc -o LaserAnalysisEB

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>


#include "TChain.h"
#include "TGraph.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTimeStamp.h"


#include "./LaserDataAnalysis.h"

using namespace std;

int main(int argc, char ** argv)
{
  int selected_fed = atoi(argv[1]);
  
  if ( selected_fed < 610 ||  selected_fed > 645){
    cout << ">>>>>  Fed  " << selected_fed << " is not in EB !!!"  << endl;
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
  
  // 36 --> number of feds in EE
  // 1700 : max number of channels in the fed
  TGraph* gapdpn1[36][1700];
  TGraph* gapdpn2[36][1700];
  TGraph* gratio[36][1700];
  char gname[100];
  int fed = 0;
  for (int ism = 0; ism < 36; ism++){
    if (ism < 18 ) fed = 628 + ism;
    if (ism >= 18) fed = 628 - 36 + ism ; 
    cout << (ism+1) << "  "  << fed << endl;
    for (int ixtal = 0; ixtal < 1700; ixtal++){
      sprintf(gname,"g_apdpn1_fed%d_xtal%d", fed, ixtal);
      gapdpn1[ism][ixtal] = new TGraph();
      gapdpn1[ism][ixtal] -> SetTitle(gname);
      gapdpn1[ism][ixtal] -> SetName(gname);
      gapdpn1[ism][ixtal] -> SetMarkerColor(kBlue);
      gapdpn1[ism][ixtal] -> SetMarkerStyle(20);
      gapdpn1[ism][ixtal] -> SetMarkerSize(0.5);
      
      sprintf(gname,"g_apdpn2_fed%d_xtal%d", fed, ixtal);
      gapdpn2[ism][ixtal] = new TGraph();
      gapdpn2[ism][ixtal] -> SetTitle(gname);
      gapdpn2[ism][ixtal] -> SetName(gname);
      gapdpn2[ism][ixtal] -> SetMarkerColor(kRed+1);
      gapdpn2[ism][ixtal] -> SetMarkerStyle(20);
      gapdpn2[ism][ixtal] -> SetMarkerSize(0.5);
      
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
  TTimeStamp dateMin(2011, 8,  15, 0, kTRUE, 0); 
  TTimeStamp dateMax(2011, 10, 15, 0, kTRUE, 0); 
  
  // profile averaging on one harness
  TProfile *p_apdpn1[36][20];
  TProfile *p_apdpn2[36][20];
  TProfile *p_ratio[36][20];
  
  for (int ism = 0; ism < 36; ism++){
    if (ism < 18 ) fed = 628 + ism;
    if (ism >= 18) fed = 628 - 36 + ism ; 
    for (int ih = 0; ih < 20; ih++){
      
      sprintf(gname,"p_apdpn1_fed%d_harness%d", fed, ih+1);
      p_apdpn1[ism][ih] = new TProfile(gname,gname,2000,dateMin, dateMax, 0,100);
      p_apdpn1[ism][ih] ->SetLineColor(kBlue);
      p_apdpn1[ism][ih] ->SetMarkerColor(kBlue);
      p_apdpn1[ism][ih] ->SetMarkerStyle(20);
      p_apdpn1[ism][ih] ->SetMarkerSize(0.5);
      
      sprintf(gname,"p_apdpn2_fed%d_harness%d", fed, ih+1);
      p_apdpn2[ism][ih] = new TProfile(gname,gname,2000,dateMin, dateMax, 0,100);
      p_apdpn2[ism][ih] ->SetLineColor(kRed+1);
      p_apdpn2[ism][ih] ->SetMarkerColor(kRed+1);
      p_apdpn2[ism][ih] ->SetMarkerStyle(20);
      p_apdpn2[ism][ih] ->SetMarkerSize(0.5);

      sprintf(gname,"p_ratio_fed%d_harness%d", fed, ih+1);
      p_ratio[ism][ih] = new TProfile(gname,gname,2000,dateMin, dateMax, 0,100);
      p_ratio[ism][ih] ->SetLineColor(kBlack);
      p_ratio[ism][ih] ->SetMarkerColor(kBlack);
      p_ratio[ism][ih] ->SetMarkerStyle(20);
      p_ratio[ism][ih] ->SetMarkerSize(0.5);
    }
  }
    
  int evt1[36][1700] = {0};
  int evt2[36][1700] = {0};

  int tAve;

  float apdpn1ref[36][1700], apdpn2ref[36][1700];
  
  for (int ientry = 0; ientry < tx->GetEntries(); ientry++){
    tx->GetEntry(ientry);
    
    if (ientry%1000000 == 0 ) cout << "Analyzing entry " << ientry << endl;
    
    // select only EB
    if (!isEB(x.fed)) continue;
        
    int harn     = x.harness;
    int fed      = x.fed;
    int ism      = iSM(fed);
    float apdpn1 = x.apdpnAB[0];
    float apdpn2 = x.apdpnAB[1];
    int t1       = x.time[0];
    int t2       = x.time[1];
    int xtal     = x.elecId;      


    // check only one fed
    if ( fed != selected_fed ) continue;
     
    if (apdpn1<=0) continue;
    if (apdpn2<=0) continue;
   
    // normalize to first point 
    if (evt1[ism][xtal] == 0) apdpn1ref[ism][xtal] = apdpn1;
    if (evt2[ism][xtal] == 0) apdpn2ref[ism][xtal] = apdpn2;
    
    if (evt1[ism][xtal] > 0) apdpn1/=apdpn1ref[ism][xtal] ;
    if (evt2[ism][xtal] > 0) apdpn2/=apdpn2ref[ism][xtal] ;
    

    tAve = t1 + (t2 - t1)/ 2;
    if ( t1 > t2 ) tAve = t2 + (t1 - t2)/ 2;

    p_apdpn1[ism][harn-1]->Fill(t1, apdpn1);
    p_apdpn2[ism][harn-1]->Fill(t2, apdpn2);
    p_ratio[ism][harn-1]   ->Fill(tAve, apdpn1/apdpn2);
  
    gapdpn1[ism][xtal]->SetPoint(evt1[ism][xtal], t1 , apdpn1);
    
    gapdpn2[ism][xtal]->SetPoint(evt2[ism][xtal], t2 , apdpn2);

    gratio[ism][xtal]->SetPoint(evt2[ism][xtal], tAve , apdpn1/apdpn2);

    evt1[ism][xtal]++;    
    evt2[ism][xtal]++;
  }
  


//   // ratio plots
//   for (int ism = 0; ism < 36; ism++){
//     for (int ixtal = 0; ixtal < 1700; ixtal++){
//       if ( gapdpn1[ism][ixtal]-> GetN() == 0 ) continue();
//       if ( gapdpn2[ism][ixtal]-> GetN() == 0 ) continue();
      
//       for (int ibin = 1; ibin < gapdpn1[ism][ixtal]-> GetNbinsX()+1; ibin++){
// 	gratio[ism][ixtal]-> Fill();
//       }
//     }
//   }




  char fname[100];
  sprintf(fname,"EBlaserAnalysis_fed%d.root",selected_fed);

  TFile *fout = new TFile(fname,"recreate");
  
  for (int ism = 0; ism < 36; ism++){
    for (int ixtal = 0; ixtal < 1700; ixtal++){
      if ( gapdpn1[ism][ixtal]-> GetN() > 0) gapdpn1[ism][ixtal]-> Write();
      if ( gapdpn2[ism][ixtal]-> GetN() > 0) gapdpn2[ism][ixtal]-> Write();
      if ( gratio[ism][ixtal]-> GetN()  > 0) gratio[ism][ixtal]-> Write();
    }
    for (int ih = 0; ih < 20; ih++){
      if ( p_apdpn1[ism][ih]-> GetEntries() > 0)  	p_apdpn1[ism][ih]->Write();
      if ( p_apdpn2[ism][ih]-> GetEntries() > 0)  	p_apdpn2[ism][ih]->Write();
      if ( p_ratio[ism][ih]-> GetEntries()  > 0)  	p_ratio[ism][ih]->Write();
    }
  }


}
