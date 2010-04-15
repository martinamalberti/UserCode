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
#include "TProfile.h"
#include "TFile.h"
#include "TVector3.h"

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
  // std::string outputRootName = "spikeAnalysis_EcalPhase-7ns.root" ;
  std::string outputRootName = "chi2Analysis.root" ;
  
  // output histos

  TH1F *hRun = new TH1F("hRun","hRun",300,132400,132700);

  TH1F *hChi2[11];
  TH1F *hChi2_EBp[11];
  TH1F *hChi2_EBm[11];

  TH1F *hChi2_Normal[11];
  TH1F *hChi2_EBp_Normal[11];
  TH1F *hChi2_EBm_Normal[11];

  TH1F *hChi2_Spike[11];
  TH1F *hChi2_EBp_Spike[11];
  TH1F *hChi2_EBm_Spike[11];

  TH1F *hOutOfTimeChi2[11];
  TH1F *hOutOfTimeChi2_Normal[11];
  TH1F *hOutOfTimeChi2_Spike[11];

  char hname[110];
  for (int i = 0; i < 11 ; i++){
    sprintf(hname,"hChi2_%d",i);
    hChi2[i] = new TH1F(hname, "chi^2{2}", 100,0,100);
    hChi2[i] ->SetLineColor(i+1);
    hChi2[i] ->GetXaxis()->SetTitle("chi^{2}");
 
    sprintf(hname,"hChi2_EBp_%d",i);
    hChi2_EBp[i] = new TH1F(hname, "chi^{2} (EB+)", 100,0,100);
    hChi2_EBp[i] ->SetLineColor(i+1);
    hChi2_EBp[i] ->GetXaxis()->SetTitle("chi^{2}");

    sprintf(hname,"hChi2_EBm_%d",i);
    hChi2_EBm[i] = new TH1F(hname,"chi^{2} (EB-)", 100,0,100);
    hChi2_EBm[i] ->SetLineColor(i+1);
    hChi2_EBm[i] ->GetXaxis()->SetTitle("chi^{2}");

    sprintf(hname,"hChi2_Normal_%d",i);
    hChi2_Normal[i] = new TH1F(hname, "chi^2{2}", 100,0,100);
    hChi2_Normal[i] ->SetLineColor(i+1);
    hChi2_Normal[i] ->GetXaxis()->SetTitle("chi^{2}");
 
    sprintf(hname,"hChi2_EBp_Normal_%d",i);
    hChi2_EBp_Normal[i] = new TH1F(hname, "chi^{2} (EB+)", 100,0,100);
    hChi2_EBp_Normal[i] ->SetLineColor(i+1);
    hChi2_EBp_Normal[i] ->GetXaxis()->SetTitle("chi^{2}");

    sprintf(hname,"hChi2_EBm_Normal_%d",i);
    hChi2_EBm_Normal[i] = new TH1F(hname,"chi^{2} (EB-)", 100,0,100);
    hChi2_EBm_Normal[i] ->SetLineColor(i+1);
    hChi2_EBm_Normal[i] ->GetXaxis()->SetTitle("chi^{2}");

    sprintf(hname,"hChi2_Spike_%d",i);
    hChi2_Spike[i] = new TH1F(hname, "chi^2{2}", 100,0,100);
    hChi2_Spike[i] ->SetLineColor(i+1);
    hChi2_Spike[i] ->GetXaxis()->SetTitle("chi^{2}");
 
    sprintf(hname,"hChi2_EBp_Spike_%d",i);
    hChi2_EBp_Spike[i] = new TH1F(hname, "chi^{2} (EB+)", 100,0,100);
    hChi2_EBp_Spike[i] ->SetLineColor(i+1);
    hChi2_EBp_Spike[i] ->GetXaxis()->SetTitle("chi^{2}");

    sprintf(hname,"hChi2_EBm_Spike_%d",i);
    hChi2_EBm_Spike[i] = new TH1F(hname,"chi^{2} (EB-)", 100,0,100);
    hChi2_EBm_Spike[i] ->SetLineColor(i+1);
    hChi2_EBm_Spike[i] ->GetXaxis()->SetTitle("chi^{2}");
 
    sprintf(hname,"hOutOfTimeChi2_%d",i);
    hOutOfTimeChi2[i] = new TH1F(hname, "OutOfTime chi^2{2}", 100,0,100);
    hOutOfTimeChi2[i] ->SetLineColor(i+1);
    hOutOfTimeChi2[i] ->GetXaxis()->SetTitle("OutOfTime chi^{2}");

    sprintf(hname,"hOutOfTimeChi2_Normal_%d",i);
    hOutOfTimeChi2_Normal[i] = new TH1F(hname, "OutOfTime chi^{2}", 100,0,100);
    hOutOfTimeChi2_Normal[i] ->SetLineColor(i+1);
    hOutOfTimeChi2_Normal[i] ->GetXaxis()->SetTitle("OutOfTime chi^{2}");
   
    sprintf(hname,"hOutOfTimeChi2_Spike_%d",i);
    hOutOfTimeChi2_Spike[i] = new TH1F(hname, "OutOfTime chi^{2}", 100,0,100);
    hOutOfTimeChi2_Spike[i] ->SetLineColor(i+1);
    hOutOfTimeChi2_Spike[i] ->GetXaxis()->SetTitle("OutOfTime chi^{2}");

  }

  TH2F *hChi2vsIeta = new TH2F("hChi2vsIeta","Chi^{2} vs Ieta",172,-85,85,100,0,100);
  TH2F *hTimevsIeta = new TH2F("hTimevsIeta","Time vs Ieta",172,-85,85,800,-100,100);
 
  TProfile *pChi2vsIeta = new TProfile("pChi2vsIeta","Chi^{2} vs Ieta",172,-85,85,0,100);
  TProfile *pTimevsIeta = new TProfile("pTimevsIeta","Time vs Ieta",172,-85,85,-100,100);

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

	//if (ecalRecHitEnergy[ihit] < 3) continue;
      
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

      if ( ecalRecHitEnergy[ihit] > 0.6 &&  ecalRecHitR9[ihit] < 0.9) {
	if (ecalRecHitChi2[ihit]!=64) hChi2vsIeta->Fill(ecalRecHitIEta[ihit],ecalRecHitChi2[ihit]);
	if (ecalRecHitChi2[ihit]!=64) pChi2vsIeta->Fill(ecalRecHitIEta[ihit],ecalRecHitChi2[ihit]);
	if (fabs(ecalRecHitTime[ihit])<4  ) hTimevsIeta->Fill(ecalRecHitIEta[ihit],ecalRecHitTime[ihit]); 
	if (fabs(ecalRecHitTime[ihit])<4  ) pTimevsIeta->Fill(ecalRecHitIEta[ihit],ecalRecHitTime[ihit]);
      }

      for (int i = 0; i < 11; i++){
	float eCut = (float)i;
	if ( ecalRecHitEnergy[ihit] > eCut) hChi2[i] ->Fill(ecalRecHitChi2[ihit]);
	if ( ecalRecHitEnergy[ihit] > eCut) hOutOfTimeChi2[i] ->Fill(ecalRecHitOutOfTimeChi2[ihit]);
	if ( ecalRecHitEnergy[ihit] > eCut &&  ecalRecHitIEta[ihit] > 0) hChi2_EBp[i] ->Fill(ecalRecHitChi2[ihit]);
	if ( ecalRecHitEnergy[ihit] > eCut &&  ecalRecHitIEta[ihit] < 0) hChi2_EBm[i] ->Fill(ecalRecHitChi2[ihit]);

	if ( ecalRecHitR9[ihit] < 0.9){
	
	  if ( ecalRecHitEnergy[ihit] > eCut) hChi2_Normal[i] ->Fill(ecalRecHitChi2[ihit]);
	  if ( ecalRecHitEnergy[ihit] > eCut) hOutOfTimeChi2_Normal[i] ->Fill(ecalRecHitOutOfTimeChi2[ihit]);
	  if ( ecalRecHitEnergy[ihit] > eCut &&  ecalRecHitIEta[ihit] > 0) hChi2_EBp_Normal[i] ->Fill(ecalRecHitChi2[ihit]);
	  if ( ecalRecHitEnergy[ihit] > eCut &&  ecalRecHitIEta[ihit] < 0) hChi2_EBm_Normal[i] ->Fill(ecalRecHitChi2[ihit]);
	  
	}

	if ( ecalRecHitR9[ihit] > 0.9){
	
	  if ( ecalRecHitEnergy[ihit] > eCut) hChi2_Spike[i] ->Fill(ecalRecHitChi2[ihit]);
	  if ( ecalRecHitEnergy[ihit] > eCut) hOutOfTimeChi2_Spike[i] ->Fill(ecalRecHitOutOfTimeChi2[ihit]);
	  if ( ecalRecHitEnergy[ihit] > eCut &&  ecalRecHitIEta[ihit] > 0) hChi2_EBp_Spike[i] ->Fill(ecalRecHitChi2[ihit]);
	  if ( ecalRecHitEnergy[ihit] > eCut &&  ecalRecHitIEta[ihit] < 0) hChi2_EBm_Spike[i] ->Fill(ecalRecHitChi2[ihit]);
	  
	}




      }

      }//end loop over rec hits
      
    } // loop over entries

  
  TFile saving (outputRootName.c_str (),"recreate") ;
  saving.cd () ;  
  
  // saving distributions
  hRun->Write();

  hChi2vsIeta->Write();
  hTimevsIeta->Write();
  pChi2vsIeta->Write();
  pTimevsIeta->Write();

  for (int i = 0; i < 11 ; i++){
    hChi2[i]     -> Write();
    hChi2_EBp[i] -> Write();
    hChi2_EBm[i] -> Write();
    
    hChi2_Normal[i]     -> Write();
    hChi2_EBp_Normal[i] -> Write();
    hChi2_EBm_Normal[i] -> Write();
    
    hChi2_Spike[i]     -> Write();
    hChi2_EBp_Spike[i] -> Write();
    hChi2_EBm_Spike[i] -> Write();
 
    hOutOfTimeChi2[i]     -> Write();
    hOutOfTimeChi2_Normal[i]     -> Write();
    hOutOfTimeChi2_Spike[i]     -> Write();
  }

  saving.Close () ;
 
  return 0 ;
}


