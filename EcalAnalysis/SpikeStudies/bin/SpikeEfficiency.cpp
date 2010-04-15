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
#include "TGraph.h"

int main (int argc, char** argv)
{
  

  TChain *chain = new TChain ("myanalysis/EcalAnalysisTree") ;
 
  // input files
  chain->Add("/tmp/malberti/SpikesCommissioning10_GOODCOLLV8.root");
  chain->Add("/tmp/malberti/SpikesCommissioning10_Apr1Skim_GOODCOLL-v1.root");
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
  std::string outputRootName = "spikeEfficiency_3GeV.root" ;
  
  // output histos

  TH1F *hRun = new TH1F("hRun","hRun",300,132400,132700);

  TH1F *hS1oS9  = new TH1F("hS1oS9","S1oS9 all Rec Hits",120,0,1.2); 
  hS1oS9 -> GetXaxis()-> SetTitle("S1/S9");

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
  hChi2[0] -> GetXaxis()-> SetTitle("chi^{2}");
  hChi2[1] = new TH1F("hChi2Normal","hChi2Normal",100,0,100);
  hChi2[1] ->SetLineColor(3);
  hChi2[1] -> GetXaxis()-> SetTitle("chi^{2}");

  TH1F *hOutOfTimeChi2[2];
  hOutOfTimeChi2[0] = new TH1F("hOutOfTimeChi2Spike","hOutOfTimeChi2Spike",100,0,100);
  hOutOfTimeChi2[0] -> SetLineColor(2);
  hOutOfTimeChi2[0] -> GetXaxis()-> SetTitle("OutOfTime chi^{2}");
  hOutOfTimeChi2[1] = new TH1F("hOutOfTimeChi2Normal","hOutOfTimeChi2Normal",100,0,100);
  hOutOfTimeChi2[1] ->SetLineColor(3);
  hOutOfTimeChi2[1] -> GetXaxis()-> SetTitle("OutOfTime chi^{2}");

  int nSpikeTot  = 0 ;
  int nNormalTot = 0 ;
 
  int nSpikeTotEB[2]  = {0} ;
  int nNormalTotEB[2] = {0} ;
 
  int n = 1000;

  int nSpikeChi2[1000] = {0.};
  int nNormalChi2[1000]= {0};

  int nSpikeChi2OutOfTime[1000] = {0};
  int nNormalChi2OutOfTime[1000]= {0};

  int nSpikeTime[1000] = {0};
  int nNormalTime[1000]= {0};


  int nSpikeChi2EB[1000][2] ;
  int nNormalChi2EB[1000][2];

  int nSpikeChi2OutOfTimeEB[1000][2] ;
  int nNormalChi2OutOfTimeEB[1000][2];

  int nSpikeTimeEB[1000][2] ;
  int nNormalTimeEB[1000][2];

  for (int i = 0; i < n ; i++){
    for (int j = 0; j < 2; j++){
      nSpikeChi2EB[i][j] = 0;
      nSpikeChi2OutOfTimeEB[i][j] = 0;
      nSpikeTimeEB[i][j] = 0;
      nNormalChi2EB[i][j] = 0;
      nNormalChi2OutOfTimeEB[i][j] = 0;
      nNormalTimeEB[i][j] = 0;
    }
  }

  //for time resolution parametrization
  float N = 31.;
  float ADCtoGeV = 0.03865 ;
  float rmsNoise = 1.1;
  float cterm = 0.7;

  // loop over entries
  for (int entry = 0; entry < nEntries; ++entry)
    {
      chain->GetEntry (entry) ;
      if(entry%500000 == 0) std::cout << "event n. " << entry << std::endl;
      

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

      hRun->Fill(runId);


      // c'e' un canale con problemi?
      if ( ecalRecHitIEta[ihit]==58 && ecalRecHitIPhi[ihit]==97 && runId==132569) continue; 


      if (ecalRecHitR9[ihit] > 0.95 ) {
	nSpikeTot++;
	if (ecalRecHitIEta[ihit] < 0)  nSpikeTotEB[0]++;
	if (ecalRecHitIEta[ihit] > 0)  nSpikeTotEB[1]++;
	hTime[0]->Fill(ecalRecHitTime[ihit]);
      }

      if (ecalRecHitR9[ihit] < 0.80 ) {
	nNormalTot++;
	if (ecalRecHitIEta[ihit] < 0)  nNormalTotEB[0]++;
	if (ecalRecHitIEta[ihit] > 0)  nNormalTotEB[1]++;
	hTime[1]->Fill(ecalRecHitTime[ihit]);
	//if (ecalRecHitTime[ihit]> 8) std::cout << ecalRecHitIEta[ihit] << "  " << ecalRecHitIPhi[ihit] << "  " << runId<< std::endl;
	
      }


      for (int i = 0 ; i < n ; i++){

	float chi2cut = (i+1) * 65./n ;

	// for the time cut , use parametrization by Giovanni et al.
	float A = ecalRecHitEnergy[ihit]/ADCtoGeV; 
	float sigmat = sqrt(pow(N/(A/rmsNoise),2) + cterm*cterm );
	float timecut = (i+1)*20./n;

	if ( ecalRecHitR9[ihit] > 0.95 ){

	  if (ecalRecHitChi2[ihit] < chi2cut) nSpikeChi2[i]++;
	  if (ecalRecHitChi2[ihit] < chi2cut  && ecalRecHitIEta[ihit] < 0) nSpikeChi2EB[i][0]++;
	  if (ecalRecHitChi2[ihit] < chi2cut  && ecalRecHitIEta[ihit] > 0) nSpikeChi2EB[i][1]++;
	
	  if (ecalRecHitOutOfTimeChi2[ihit] < chi2cut) nSpikeChi2OutOfTime[i]++;
	  if (ecalRecHitOutOfTimeChi2[ihit] < chi2cut  && ecalRecHitIEta[ihit] < 0) nSpikeChi2OutOfTimeEB[i][0]++;
	  if (ecalRecHitOutOfTimeChi2[ihit] < chi2cut  && ecalRecHitIEta[ihit] > 0) nSpikeChi2OutOfTimeEB[i][1]++;
	  
	  if (fabs(ecalRecHitTime[ihit]/sigmat) < timecut) nSpikeTime[i]++;
	  if (fabs(ecalRecHitTime[ihit]/sigmat ) < timecut && ecalRecHitIEta[ihit] < 0) nSpikeTimeEB[i][0]++;
	  if (fabs(ecalRecHitTime[ihit]/sigmat ) < timecut && ecalRecHitIEta[ihit] > 0) nSpikeTimeEB[i][1]++;

	}

	if ( ecalRecHitR9[ihit] < 0.80 ){
	  


	  if (ecalRecHitChi2[ihit] < chi2cut) nNormalChi2[i]++;
	  if (ecalRecHitChi2[ihit] < chi2cut  && ecalRecHitIEta[ihit] < 0) nNormalChi2EB[i][0]++;
	  if (ecalRecHitChi2[ihit] < chi2cut  && ecalRecHitIEta[ihit] > 0) nNormalChi2EB[i][1]++;
	
	  if (ecalRecHitOutOfTimeChi2[ihit] < chi2cut) nNormalChi2OutOfTime[i]++;
	  if (ecalRecHitOutOfTimeChi2[ihit] < chi2cut  && ecalRecHitIEta[ihit] < 0) nNormalChi2OutOfTimeEB[i][0]++;
	  if (ecalRecHitOutOfTimeChi2[ihit] < chi2cut  && ecalRecHitIEta[ihit] > 0) nNormalChi2OutOfTimeEB[i][1]++;
	  
	  if (fabs(ecalRecHitTime[ihit]/sigmat ) < timecut) nNormalTime[i]++;
	  if (fabs(ecalRecHitTime[ihit]/sigmat ) < timecut && ecalRecHitIEta[ihit] < 0) nNormalTimeEB[i][0]++;
	  if (fabs(ecalRecHitTime[ihit]/sigmat ) < timecut && ecalRecHitIEta[ihit] > 0) nNormalTimeEB[i][1]++;


	  
	}
      }


      }
    } // loop over entries

  TGraph *gTime = new TGraph();
  gTime->SetMarkerStyle(21);
  gTime->SetMarkerSize(0.6);
  gTime->SetMarkerColor(2);

  TGraph *gChi2 = new TGraph();
  gChi2->SetMarkerStyle(21);
  gChi2->SetMarkerSize(0.6);
  gChi2->SetMarkerColor(3);

  TGraph *gOutOfTimeChi2 = new TGraph();
  gOutOfTimeChi2->SetMarkerStyle(21);
  gOutOfTimeChi2->SetMarkerSize(0.6);
  gOutOfTimeChi2->SetMarkerColor(4);


  TGraph *gTimeEBP = new TGraph();
  gTimeEBP->SetMarkerStyle(20);
  gTimeEBP->SetMarkerSize(0.7);
  gTimeEBP->SetMarkerColor(2);

  TGraph *gChi2EBP = new TGraph();
  gChi2EBP->SetMarkerStyle(20);
  gChi2EBP->SetMarkerSize(0.7);
  gChi2EBP->SetMarkerColor(3);

  TGraph *gOutOfTimeChi2EBP = new TGraph();
  gOutOfTimeChi2EBP->SetMarkerStyle(20);
  gOutOfTimeChi2EBP->SetMarkerSize(0.7);
  gOutOfTimeChi2EBP->SetMarkerColor(4);
 
  TGraph *gTimeEBM = new TGraph();
  gTimeEBM->SetMarkerStyle(24);
  gTimeEBM->SetMarkerSize(0.7);
  gTimeEBM->SetMarkerColor(2);

  TGraph *gChi2EBM = new TGraph();
  gChi2EBM->SetMarkerStyle(24);
  gChi2EBM->SetMarkerSize(0.7);
  gChi2EBM->SetMarkerColor(3);

  TGraph *gOutOfTimeChi2EBM = new TGraph();
  gOutOfTimeChi2EBM->SetMarkerStyle(24);
  gOutOfTimeChi2EBM->SetMarkerSize(0.7);
  gOutOfTimeChi2EBM->SetMarkerColor(4);


  float effS, effN;

  for (int i = 0; i < n ; i++){

    effN = (float)nNormalTime[i]/(float)nNormalTot;
    effS = (float)nSpikeTime[i]/(float)nSpikeTot;
    gTime->SetPoint(i,effS,effN);

    effN = (float)nNormalChi2[i]/(float)nNormalTot;
    effS = (float)nSpikeChi2[i]/(float)nSpikeTot;
    gChi2->SetPoint(i,effS,effN);

    effN = (float)nNormalChi2OutOfTime[i]/(float)nNormalTot;
    effS = (float)nSpikeChi2OutOfTime[i]/(float)nSpikeTot;
    gOutOfTimeChi2->SetPoint(i,effS,effN);

    //EB-
    effN = (float)nNormalTimeEB[i][0]/(float)nNormalTotEB[0];
    effS = (float)nSpikeTimeEB[i][0]/(float)nSpikeTotEB[0];
    gTimeEBM->SetPoint(i,effS,effN);

    effN = (float)nNormalChi2EB[i][0]/(float)nNormalTotEB[0];
    effS = (float)nSpikeChi2EB[i][0]/(float)nSpikeTotEB[0];
    gChi2EBM->SetPoint(i,effS,effN);

    effN = (float)nNormalChi2OutOfTimeEB[i][0]/(float)nNormalTotEB[0];
    effS = (float)nSpikeChi2OutOfTimeEB[i][0]/(float)nSpikeTotEB[0];
    gOutOfTimeChi2EBM->SetPoint(i,effS,effN);

    //EB+
    effN = (float)nNormalTimeEB[i][1]/(float)nNormalTotEB[1];
    effS = (float)nSpikeTimeEB[i][1]/(float)nSpikeTotEB[1];
    gTimeEBP->SetPoint(i,effS,effN);

    effN = (float)nNormalChi2EB[i][1]/(float)nNormalTotEB[1];
    effS = (float)nSpikeChi2EB[i][1]/(float)nSpikeTotEB[1];
    gChi2EBP->SetPoint(i,effS,effN);

    effN = (float)nNormalChi2OutOfTimeEB[i][1]/(float)nNormalTotEB[1];
    effS = (float)nSpikeChi2OutOfTimeEB[i][1]/(float)nSpikeTotEB[1];
    gOutOfTimeChi2EBP->SetPoint(i,effS,effN);

  }




  TFile saving (outputRootName.c_str (),"recreate") ;
  saving.cd () ;  
  
  // saving distributions
  hRun->Write();

  gTime->Write("g_Time");
  gChi2->Write("g_Chi2");
  gOutOfTimeChi2->Write("g_OutOfTimeChi2");

  gTimeEBP->Write("g_Time_EBP");
  gChi2EBP->Write("g_Chi2_EBP");
  gOutOfTimeChi2EBP->Write("g_OutOfTimeChi2_EBP");

  gTimeEBM->Write("g_Time_EBM");
  gChi2EBM->Write("g_Chi2_EBM");
  gOutOfTimeChi2EBM->Write("g_OutOfTimeChi2_EBM");


  hTime[0]->Write();
  hTime[1]->Write();

  saving.Close () ;
 
  return 0 ;
}


