#include "EcalAnalysis/SpikeStudies/interface/EcalTreeContent.h"

#include <iostream>
#include <fstream>
#include <string>
#include <boost/foreach.hpp>
#include <map>

#include "TMath.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TVector3.h"

int main (int argc, char** argv)
{
  

  TChain *chain = new TChain ("myanalysis/EcalAnalysisTree") ;
 
  // input files
  chain->Add("/tmp/malberti/SpikesCommissioning2010_GOODCOLLV8.root");
  int nEntries = chain->GetEntries () ;
  std::cout << "FOUND " << nEntries << " ENTRIES\n" ;    

  //ecalVariables variables
  unsigned int BX;
  unsigned int lumiId;
  unsigned int runId;
  unsigned int eventId;
  unsigned int eventNaiveId;

  int nEcalRecHits;
  float ecalRecHitType[1000];
  float ecalRecHitEnergy[1000];
  float ecalRecHitOutOfTimeEnergy[1000];
  float ecalRecHitIEta[1000];
  float ecalRecHitIPhi[1000];
  float ecalRecHitTime[1000];
  float ecalRecHitChi2[1000];
  float ecalRecHitOutOfTimeChi2[1000];
  int ecalRecHitRawId[1000];
  float ecalRecHitCoeff[1000];
  int ecalRecHitRecoFlag[1000];
  float ecalRecHitR9[1000];
  float ecalRecHitS4oS1[1000];
  float ecalRecHitIso03[1000][2];
  float ecalRecHitIso04[1000][2];
  int ecalDigis[1000][10];
  int ecalGainId[1000][10];
  float ecalRecHitMatrix[1000][5][5];
  float ecalRecHitMatrixFlag[1000][5][5];

  // L1 trigger variables
  int techL1Bit[64];
  int algoL1Bit[128];
  
  chain -> SetBranchAddress("BX",            &BX);
  chain -> SetBranchAddress("lumiId",        &lumiId);
  chain -> SetBranchAddress("runId",         &runId);
  chain -> SetBranchAddress("eventId",       &eventId);
  chain -> SetBranchAddress("eventNaiveId",  &eventNaiveId);

  chain -> SetBranchAddress("nEcalRecHits",             &nEcalRecHits);
  chain -> SetBranchAddress("ecalRecHitType",            ecalRecHitType);
  chain -> SetBranchAddress("ecalRecHitEnergy",          ecalRecHitEnergy);
  chain -> SetBranchAddress("ecalRecHitOutOfTimeEnergy", ecalRecHitOutOfTimeEnergy);
  chain -> SetBranchAddress("ecalRecHitIEta",            ecalRecHitIEta);
  chain -> SetBranchAddress("ecalRecHitIPhi",            ecalRecHitIPhi);
  chain -> SetBranchAddress("ecalRecHitTime",            ecalRecHitTime);
  chain -> SetBranchAddress("ecalRecHitChi2",            ecalRecHitChi2);
  chain -> SetBranchAddress("ecalRecHitOutOfTimeChi2",   ecalRecHitOutOfTimeChi2);
  chain -> SetBranchAddress("ecalRecHitRawId",           ecalRecHitRawId);
  chain -> SetBranchAddress("ecalRecHitCoeff",           ecalRecHitCoeff);
  chain -> SetBranchAddress("ecalRecHitRecoFlag",        ecalRecHitRecoFlag);
  chain -> SetBranchAddress("ecalRecHitR9",              ecalRecHitR9);
  chain -> SetBranchAddress("ecalRecHitS4oS1",           ecalRecHitS4oS1);
  chain -> SetBranchAddress("ecalRecHitIso03",           ecalRecHitIso03);
  chain -> SetBranchAddress("ecalRecHitIso04",           ecalRecHitIso04);
  chain -> SetBranchAddress("ecalDigis",                 ecalDigis);
  chain -> SetBranchAddress("ecalGainId",                ecalGainId);
  chain -> SetBranchAddress("ecalRecHitMatrix",          ecalRecHitMatrix);
  chain -> SetBranchAddress("ecalRecHitMatrixFlag",      ecalRecHitMatrixFlag);
  
  chain -> SetBranchAddress("techL1Bit",     techL1Bit);
  chain -> SetBranchAddress("algoL1Bit",     algoL1Bit);
 
  // output file
  std::string outputRootName = "spikeAnalysis.root" ;
  
  // output histos

  TH1F *hRun = new TH1F("hRun","hRun",300,132400,132700);

  TH1F *hS1oS9  = new TH1F("hS1oS9","hS1oS9",120,0,1.2); 
  hS1oS9 -> GetXaxis()-> SetTitle("S1/S9");

  TH1F *hS4oS1 = new TH1F("hS4oS1","hS4oS1",210,-0.1,2); 
  hS4oS1 -> GetXaxis()-> SetTitle("S4/S1");
  
  TH1F *hTime[2];
  hTime[0] = new TH1F("hTimeSpike","hTimeSpike",800,-100,100);
  hTime[0]->SetLineColor(2);
  hTime[0]-> GetXaxis()-> SetTitle("t(ns)");
  hTime[1] = new TH1F("hTimeNormal","hTimeNormal",800,-100,100);
  hTime[1]->SetLineColor(3);
  hTime[1]-> GetXaxis()-> SetTitle("t(ns)");

  TH1F *hChi2[2];
  hChi2[0] = new TH1F("hChi2Spike","hChi2Spike",100,0,100);
  hChi2[0] -> SetLineColor(2);
  hChi2[0] -> GetXaxis()-> SetTitle("chi^{2}/ndf");
  hChi2[1] = new TH1F("hChi2Normal","hChi2Normal",100,0,100);
  hChi2[1] ->SetLineColor(3);
  hChi2[1] -> GetXaxis()-> SetTitle("chi^{2}/ndf");

  TH1F *hOutOfTimeChi2[2];
  hOutOfTimeChi2[0] = new TH1F("hOutOfTimeChi2Spike","hOutOfTimeChi2Spike",100,0,100);
  hOutOfTimeChi2[0] -> SetLineColor(2);
  hOutOfTimeChi2[0] -> GetXaxis()-> SetTitle("chi^{2}/ndf");
  hOutOfTimeChi2[1] = new TH1F("hOutOfTimeChi2Normal","hOutOfTimeChi2Normal",100,0,100);
  hOutOfTimeChi2[1] ->SetLineColor(3);
  hOutOfTimeChi2[1] -> GetXaxis()-> SetTitle("chi^{2}/ndf");

  TH1F *hRatioEnergy[2];
  hRatioEnergy[0] = new TH1F("hRatioEnergySpike","hRatioEnergySpike",800,-100,100);
  hRatioEnergy[0]->SetLineColor(2);
  hRatioEnergy[0]-> GetXaxis()-> SetTitle("t(ns)");
  hRatioEnergy[1] = new TH1F("hRatioEnergyNormal","hRatioEnergyNormal",800,-100,100);
  hRatioEnergy[1]->SetLineColor(3);
  hRatioEnergy[1]-> GetXaxis()-> SetTitle("t(ns)");


  TH2F *hTime_vs_S1oS9 = new TH2F("hTime_vs_S1oS9","hTime_vs_S1oS9",150, 0.0,1.5,800,-100,100);
  hTime_vs_S1oS9->GetXaxis()->SetTitle("S1/S9");
  hTime_vs_S1oS9->GetYaxis()->SetTitle("time(ns)");

  TH2F *hTime_vs_S4oS1 = new TH2F("hTime_vs_S4oS1","hTime_vs_S4oS1",210,-0.1,2.0,800,-100,100);
  hTime_vs_S4oS1->GetXaxis()->SetTitle("S4/S1");
  hTime_vs_S4oS1->GetYaxis()->SetTitle("time(ns)");

  TH2F *hChi2_vs_S1oS9 = new TH2F("hChi2_vs_S1oS9","hChi2_vs_S1oS9",150, 0.0,1.5,100,0,100);
  hChi2_vs_S1oS9->GetXaxis()->SetTitle("S1/S9");
  hChi2_vs_S1oS9->GetYaxis()->SetTitle("chi^{2}");

  TH2F *hChi2_vs_S4oS1 = new TH2F("hChi2_vs_S4oS1","hChi2_vs_S4oS1",210,-0.1,2.0,100,0,100);
  hChi2_vs_S4oS1->GetXaxis()->SetTitle("S4/S1");
  hChi2_vs_S4oS1->GetYaxis()->SetTitle("chi^{2}");

  TH2F *hOutOfTimeChi2_vs_S1oS9 = new TH2F("hOutOfTimeChi2_vs_S1oS9","hOutOfTimeChi2_vs_S1oS9",150, 0.0,1.5,100,0,100);
  hOutOfTimeChi2_vs_S1oS9->GetXaxis()->SetTitle("S1/S9");
  hOutOfTimeChi2_vs_S1oS9->GetYaxis()->SetTitle("out-of-time chi^{2}");
  
  TH2F *hOutOfTimeChi2_vs_S4oS1 = new TH2F("hOutOfTimeChi2_vs_S4oS1","hOutOfTimeChi2_vs_S4oS1",210,-0.1,2.0,100,0,100);
  hOutOfTimeChi2_vs_S4oS1->GetXaxis()->SetTitle("S4/S1");
  hOutOfTimeChi2_vs_S4oS1->GetYaxis()->SetTitle("out-of-time chi^{2}");

  TH2F *hRatioEnergy_vs_S1oS9 = new TH2F("hRatioEnergy_vs_S1oS9","hRatioEnergy_vs_S1oS9",150, 0.0,1.5,100,0,100);
  hRatioEnergy_vs_S1oS9->GetXaxis()->SetTitle("S1/S9");
  hRatioEnergy_vs_S1oS9->GetYaxis()->SetTitle("energy/energyOutOfTime");

  TH2F *hRatioEnergy_vs_S4oS1 = new TH2F("hRatioEnergy_vs_S4oS1","hRatioEnergy_vs_S4oS1",210,-0.1,2.0,100,0,100);
  hRatioEnergy_vs_S4oS1->GetXaxis()->SetTitle("S4/S1");
  hRatioEnergy_vs_S4oS1->GetYaxis()->SetTitle("energy/energyOutOfTime");

  TH2F *hTime_vs_Chi2[2];
  hTime_vs_Chi2[0] = new TH2F("hTime_vs_Chi2_Spike","hTime_vs_Chi2_Spike",100, 0,100,800,-100,100);
  hTime_vs_Chi2[0] ->GetXaxis()->SetTitle("Chi^{2}");
  hTime_vs_Chi2[0] ->GetYaxis()->SetTitle("time(ns)");
  hTime_vs_Chi2[0] ->SetMarkerColor(2);
  hTime_vs_Chi2[1] = new TH2F("hTime_vs_Chi2_Normal","hTime_vs_Chi2_Normal",100, 0,100,800,-100,100);
  hTime_vs_Chi2[1] ->GetXaxis()->SetTitle("Chi^{2}");
  hTime_vs_Chi2[1] ->GetYaxis()->SetTitle("time(ns)");
  hTime_vs_Chi2[1] ->SetMarkerColor(3);

  TH2F *hTime_vs_OutOfTimeChi2[2];
  hTime_vs_OutOfTimeChi2[0] = new TH2F("hTime_vs_OutOfTimeChi2_Spike","hTime_vs_OutOfTimeChi2_Spike",100, 0,100,800,-100,100);
  hTime_vs_OutOfTimeChi2[0] ->GetXaxis()->SetTitle("out-of-time Chi^{2}");
  hTime_vs_OutOfTimeChi2[0] ->GetYaxis()->SetTitle("time(ns)");
  hTime_vs_OutOfTimeChi2[0] ->SetMarkerColor(2);
  hTime_vs_OutOfTimeChi2[1] = new TH2F("hTime_vs_OutOfTimeChi2_Normal","hTime_vs_OutOfTimeChi2_Normal",100, 0,100,800,-100,100);
  hTime_vs_OutOfTimeChi2[1] ->GetXaxis()->SetTitle("out-of-time Chi^{2}");
  hTime_vs_OutOfTimeChi2[1] ->GetYaxis()->SetTitle("time(ns)");
  hTime_vs_OutOfTimeChi2[1] ->SetMarkerColor(3);

  TH2F *hTime_vs_RatioEnergy[2];
  hTime_vs_RatioEnergy[0] = new TH2F("hTime_vs_RatioEnergy_Spike","hTime_vs_RatioEnergy_Spike",100, 0,100,800,-100,100);
  hTime_vs_RatioEnergy[0] ->GetXaxis()->SetTitle("energy/energyOutOfTime");
  hTime_vs_RatioEnergy[0] ->GetYaxis()->SetTitle("time(ns)");
  hTime_vs_RatioEnergy[0] ->SetMarkerColor(2);
  hTime_vs_RatioEnergy[1] = new TH2F("hTime_vs_RatioEnergy_Normal","hTime_vs_RatioEnergy_Normal",100, 0,100,800,-100,100);
  hTime_vs_RatioEnergy[1] ->GetXaxis()->SetTitle("energy/energyOutOfTime");
  hTime_vs_RatioEnergy[1] ->GetYaxis()->SetTitle("time(ns)");
  hTime_vs_RatioEnergy[1] ->SetMarkerColor(3);



  // loop over entries
  for (int entry = 0; entry < nEntries; ++entry)
    {
      chain->GetEntry (entry) ;
      if(entry%100000 == 0) std::cout << "event n. " << entry << std::endl;
      

      if ((techL1Bit[40]+techL1Bit[41])==0) continue; // c'e' gia'
             

      // select the good collisions runs + good lumi sections (no ecal time scans, etc...)
      bool goodrun = false;

      if (runId==132440 && lumiId >=86  && lumiId <= 138) goodrun = true;
      if (runId==132440 && lumiId >=141 && lumiId <= 401) goodrun = true;
   
      if (runId==132473 && lumiId >=1   && lumiId <= 29 ) goodrun = true;
      
      if (runId==132476 && lumiId >=23  && lumiId <= 28 ) goodrun = true;
      if (runId==132476 && lumiId >=54  && lumiId <= 57 ) goodrun = true;

      if (runId==132477 && lumiId >=34  && lumiId <= 35 ) goodrun = true;
      if (runId==132477 && lumiId >=63  && lumiId <= 64 ) goodrun = true;
      if (runId==132477 && lumiId >=90  && lumiId <= 93 ) goodrun = true;
      if (runId==132477 && lumiId >=118 && lumiId <= 121) goodrun = true;
      if (runId==132477 && lumiId >=148 && lumiId <= 149) goodrun = true;
      if (runId==132477 && lumiId >=176 && lumiId <= 179) goodrun = true;
      if (runId==132477 && lumiId >=225 && lumiId <= 236) goodrun = true;
      if (runId==132477 && lumiId >=368 && lumiId <= 384) goodrun = true;
      if (runId==132477 && lumiId >=517 && lumiId <= 520) goodrun = true;

      if (runId==132569 && lumiId >=222 && lumiId <= 224) goodrun = true;
      if (runId==132569 && lumiId >=310 && lumiId <= 310) goodrun = true;
      if (runId==132569 && lumiId >=411 && lumiId <= 419) goodrun = true;
      if (runId==132569 && lumiId >=529 && lumiId <= 582) goodrun = true;

      if (runId==132596 && lumiId >=383 && lumiId <= 383) goodrun = true;
      if (runId==132596 && lumiId >=447 && lumiId <= 453) goodrun = true;

      if (runId==132598 && lumiId >=80  && lumiId <= 82 ) goodrun = true;
      if (runId==132598 && lumiId >=174 && lumiId <= 188) goodrun = true;

      if (runId==132599 && lumiId >=1   && lumiId <= 74 ) goodrun = true;

      if (runId==132601 && lumiId >=261 && lumiId <= 1131)goodrun = true;

      if (runId==132602 && lumiId >=1   && lumiId <= 83)  goodrun = true;

      if (runId==132605 && lumiId >=446 && lumiId <= 622) goodrun = true;
      if (runId==132605 && lumiId >=624 && lumiId <= 829) goodrun = true;
      if (runId==132605 && lumiId >=831 && lumiId <= 968) goodrun = true;

      if (runId==132606 && lumiId >=1   && lumiId <= 37 ) goodrun = true;

      if (!goodrun) continue;      

      for (int ihit =0 ; ihit < nEcalRecHits; ihit++){

      if (ecalRecHitEnergy[ihit] < 3) continue;
      
      // check gain switch
      bool gainSwitch = false;
      for (int isample = 0; isample < 10 ; isample++){
        if ( (isample > 0 && ecalGainId[ihit][isample]!= ecalGainId[ihit][isample-1]) ||
             ecalGainId[ihit][isample] !=1 ) {
          gainSwitch = true;
          break;
        }
      }

      if (gainSwitch) continue;

      float ratioEnergy = ecalRecHitEnergy[ihit]/ecalRecHitOutOfTimeEnergy[ihit];

      hRun->Fill(runId);
      hS1oS9                 -> Fill(ecalRecHitR9[ihit] );
      hTime_vs_S1oS9         -> Fill(ecalRecHitR9[ihit],ecalRecHitTime[ihit] );
      hChi2_vs_S1oS9         -> Fill(ecalRecHitR9[ihit],ecalRecHitChi2[ihit] );
      hOutOfTimeChi2_vs_S1oS9-> Fill(ecalRecHitR9[ihit],ecalRecHitOutOfTimeChi2[ihit] );
      hRatioEnergy_vs_S1oS9  -> Fill(ecalRecHitR9[ihit],ratioEnergy );

      hS4oS1->Fill(ecalRecHitS4oS1[ihit]);
      hTime_vs_S4oS1         -> Fill(ecalRecHitS4oS1[ihit],ecalRecHitTime[ihit] );
      hChi2_vs_S4oS1         -> Fill(ecalRecHitS4oS1[ihit],ecalRecHitChi2[ihit] );
      hOutOfTimeChi2_vs_S4oS1-> Fill(ecalRecHitS4oS1[ihit],ecalRecHitOutOfTimeChi2[ihit] );
      hRatioEnergy_vs_S4oS1  -> Fill(ecalRecHitS4oS1[ihit],ratioEnergy );

      if (ecalRecHitR9[ihit] > 0.9 ) {
        hTime[0]                  -> Fill( ecalRecHitTime[ihit] );
	if (fabs(ecalRecHitTime[ihit]) >4 ) hChi2[0]                  -> Fill( ecalRecHitChi2[ihit]);
        hOutOfTimeChi2[0]         -> Fill( ecalRecHitOutOfTimeChi2[ihit]);
        hRatioEnergy[0]           -> Fill( ratioEnergy );
	hTime_vs_Chi2[0]          -> Fill( ecalRecHitChi2[ihit], ecalRecHitTime[ihit] );
	hTime_vs_OutOfTimeChi2[0] -> Fill( ecalRecHitOutOfTimeChi2[ihit], ecalRecHitTime[ihit] );
	hTime_vs_RatioEnergy[0]   -> Fill( ratioEnergy, ecalRecHitTime[ihit] );
      }

      if (ecalRecHitR9[ihit] < 0.9 ) {
	hTime[1]                  -> Fill( ecalRecHitTime[ihit] );
        if (fabs(ecalRecHitTime[ihit]) >4 ) hChi2[1]                  -> Fill( ecalRecHitChi2[ihit]);
        hOutOfTimeChi2[1]         -> Fill( ecalRecHitOutOfTimeChi2[ihit]);
        hRatioEnergy[1]           -> Fill( ratioEnergy );
	hTime_vs_Chi2[1]          -> Fill( ecalRecHitChi2[ihit], ecalRecHitTime[ihit] );
	hTime_vs_OutOfTimeChi2[1] -> Fill( ecalRecHitOutOfTimeChi2[ihit], ecalRecHitTime[ihit] );
	hTime_vs_RatioEnergy[1]   -> Fill( ratioEnergy, ecalRecHitTime[ihit] );
      }


    }


      



    } // loop over entries

  TFile saving (outputRootName.c_str (),"recreate") ;
  saving.cd () ;  
  
  // saving distributions
  hRun->Write();
  hS1oS9                  -> Write();
  hTime_vs_S1oS9          -> Write();    
  hChi2_vs_S1oS9          -> Write();   
  hOutOfTimeChi2_vs_S1oS9 -> Write();
  hRatioEnergy_vs_S1oS9   -> Write();

  hS4oS1 -> Write();
  hTime_vs_S4oS1          -> Write();    
  hChi2_vs_S4oS1          -> Write();   
  hOutOfTimeChi2_vs_S4oS1 -> Write();
  hRatioEnergy_vs_S4oS1   -> Write();


  for (int i = 0; i < 2 ; i++){
    hTime[i]                  -> Write();
    hChi2[i]                  -> Write();
    hOutOfTimeChi2[i]         -> Write();
    hRatioEnergy[i]           -> Write();
    hTime_vs_Chi2[i]          -> Write();
    hTime_vs_OutOfTimeChi2[i] -> Write();
    hTime_vs_RatioEnergy[i]   -> Write();
  }

  saving.Close () ;
 
  return 0 ;
}


