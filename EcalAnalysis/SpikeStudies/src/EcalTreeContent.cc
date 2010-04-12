#include "EcalAnalysis/SpikeStudies/interface/EcalTreeContent.h"

bool EcalTreeContent::ecalVariables       = true;
bool EcalTreeContent::preShowerVariables  = false;
bool EcalTreeContent::L1TriggerVariables  = true;



void setBranchAddresses(TTree* chain, EcalTreeContent& treeVars)
{
  chain -> SetBranchAddress("BX",            &treeVars.BX);
  chain -> SetBranchAddress("lumiId",        &treeVars.lumiId);
  chain -> SetBranchAddress("runId",         &treeVars.runId);
  chain -> SetBranchAddress("eventId",       &treeVars.eventId);
  chain -> SetBranchAddress("eventNaiveId",  &treeVars.eventNaiveId);


  // ECAL VARIABLES  
  if(EcalTreeContent::ecalVariables)
    {  
      //recHits variables
      chain -> SetBranchAddress("nEcalRecHits",             &treeVars.nEcalRecHits);
      chain -> SetBranchAddress("ecalRecHitType",            treeVars.ecalRecHitType);
      chain -> SetBranchAddress("ecalRecHitEnergy",          treeVars.ecalRecHitEnergy);
      chain -> SetBranchAddress("ecalRecHitOutOfTimeEnergy", treeVars.ecalRecHitOutOfTimeEnergy);
      chain -> SetBranchAddress("ecalRecHitIEta",            treeVars.ecalRecHitIEta);
      chain -> SetBranchAddress("ecalRecHitIPhi",            treeVars.ecalRecHitIPhi);
      chain -> SetBranchAddress("ecalRecHitTime",            treeVars.ecalRecHitTime);
      chain -> SetBranchAddress("ecalRecHitChi2",            treeVars.ecalRecHitChi2);
      chain -> SetBranchAddress("ecalRecHitOutOfTimeChi2",   treeVars.ecalRecHitOutOfTimeChi2);
      chain -> SetBranchAddress("ecalRecHitRawId",           treeVars.ecalRecHitRawId);
      chain -> SetBranchAddress("ecalRecHitCoeff",           treeVars.ecalRecHitCoeff);
      chain -> SetBranchAddress("ecalRecHitRecoFlag",        treeVars.ecalRecHitRecoFlag);
      chain -> SetBranchAddress("ecalRecHitR9",              treeVars.ecalRecHitR9);
      chain -> SetBranchAddress("ecalRecHitS4oS1",           treeVars.ecalRecHitS4oS1);
      chain -> SetBranchAddress("ecalRecHitIso03",           treeVars.ecalRecHitIso03);
      chain -> SetBranchAddress("ecalRecHitIso04",           treeVars.ecalRecHitIso04);
      chain -> SetBranchAddress("ecalDigis",                 treeVars.ecalDigis);
      chain -> SetBranchAddress("ecalGainId",                treeVars.ecalGainId);
      chain -> SetBranchAddress("ecalRecHitMatrix",          treeVars.ecalRecHitMatrix);
      chain -> SetBranchAddress("ecalRecHitMatrixFlag",      treeVars.ecalRecHitMatrixFlag);
    } // ECAL VARIABLES
  
  
  //PRESHOWER VARIABLES
  if(EcalTreeContent::preShowerVariables)
    {
      chain -> SetBranchAddress("nPreShowerRecHits",     &treeVars.nPreShowerRecHits);
      chain -> SetBranchAddress("preShowerRecHitType",   treeVars.preShowerRecHitType);
      chain -> SetBranchAddress("preShowerRecHitEnergy", treeVars.preShowerRecHitEnergy);
      chain -> SetBranchAddress("preShowerRecHitIx",     treeVars.preShowerRecHitIx);
      chain -> SetBranchAddress("preShowerRecHitIy",     treeVars.preShowerRecHitIy);
      chain -> SetBranchAddress("preShowerRecHitPlane",  treeVars.preShowerRecHitPlane);
      chain -> SetBranchAddress("preShowerRecHitStrip",  treeVars.preShowerRecHitStrip);
      chain -> SetBranchAddress("preShowerRecHitTime",   treeVars.preShowerRecHitTime);
      chain -> SetBranchAddress("preShowerRecHitRawId",  treeVars.preShowerRecHitRawId);

      chain -> SetBranchAddress("nPreShowerDigis",     &treeVars.nPreShowerDigis);
      chain -> SetBranchAddress("esDigiType",          treeVars.esDigiType);
      chain -> SetBranchAddress("esDigiIx",            treeVars.esDigiIx);
      chain -> SetBranchAddress("esDigiIy",            treeVars.esDigiIy);
      chain -> SetBranchAddress("esDigiPlane",         treeVars.esDigiPlane);
      chain -> SetBranchAddress("esDigiStrip",         treeVars.esDigiStrip);
      chain -> SetBranchAddress("esDigiRawId",         treeVars.esDigiRawId);
      chain -> SetBranchAddress("preShowerSample1",    treeVars.preShowerSample1);
      chain -> SetBranchAddress("preShowerSample2",    treeVars.preShowerSample2);
      chain -> SetBranchAddress("preShowerSample3",    treeVars.preShowerSample3);
    }  //PRESHOWER VARIABLES

  
  //L1 VARIABLES
  if(EcalTreeContent::L1TriggerVariables)
    {
      chain -> SetBranchAddress("techL1Bit",     treeVars.techL1Bit);
      chain -> SetBranchAddress("algoL1Bit",     treeVars.algoL1Bit);
    }  //L1 VARIABLES
}


 



void setBranches(TTree* chain, EcalTreeContent& treeVars)
{
  chain -> Branch("BX",            &treeVars.BX,                       "BX/i");
  chain -> Branch("lumiId",        &treeVars.lumiId,               "lumiId/i");
  chain -> Branch("runId",         &treeVars.runId,                 "runId/i");
  chain -> Branch("eventId",       &treeVars.eventId,             "eventId/i");
  chain -> Branch("eventNaiveId",  &treeVars.eventNaiveId,   "eventNaiveId/i");

  
  // ECAL VARIABLES  
  if(EcalTreeContent::ecalVariables)
  {    
    //recHits variables
    chain -> Branch("nEcalRecHits",             &treeVars.nEcalRecHits,              "nEcalRecHits/I");
    chain -> Branch("ecalRecHitType",            treeVars.ecalRecHitType,            "ecalRecHitType[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitEnergy",          treeVars.ecalRecHitEnergy,          "ecalRecHitEnergy[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitOutOfTimeEnergy", treeVars.ecalRecHitOutOfTimeEnergy, "ecalRecHitOutOfTimeEnergy[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitIEta",            treeVars.ecalRecHitIEta,        "ecalRecHitIEta[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitIPhi",            treeVars.ecalRecHitIPhi,        "ecalRecHitIPhi[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitTime",            treeVars.ecalRecHitTime,        "ecalRecHitTime[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitChi2",            treeVars.ecalRecHitChi2,         "ecalRecHitChi2[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitOutOfTimeChi2",   treeVars.ecalRecHitOutOfTimeChi2,"ecalRecHitOutOfTimeChi2[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitRawId",           treeVars.ecalRecHitRawId,      "ecalRecHitRawId[nEcalRecHits]/I");
    chain -> Branch("ecalRecHitCoeff",           treeVars.ecalRecHitCoeff,      "ecalRecHitCoeff[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitRecoFlag",        treeVars.ecalRecHitRecoFlag,   "ecalRecHitRecoFlag[nEcalRecHits]/I");
    chain -> Branch("ecalRecHitR9",              treeVars.ecalRecHitR9,         "ecalRecHitR9[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitS4oS1",           treeVars.ecalRecHitS4oS1,      "ecalRecHitS4oS1[nEcalRecHits]/F");
    chain -> Branch("ecalRecHitIso03",           treeVars.ecalRecHitIso03,      "ecalRecHitIso03[nEcalRecHits][2]/F");
    chain -> Branch("ecalRecHitIso04",           treeVars.ecalRecHitIso04,      "ecalRecHitIso04[nEcalRecHits][2]/F");
    chain -> Branch("ecalDigis",                 treeVars.ecalDigis,            "ecalDigis[nEcalRecHits][10]/I");
    chain -> Branch("ecalGainId",                treeVars.ecalGainId,           "ecalGainId[nEcalRecHits][10]/I");
    chain -> Branch("ecalRecHitMatrix",          treeVars.ecalRecHitMatrix,     "ecalRecHitMatrix[nEcalRecHits][5][5]/F");
    chain -> Branch("ecalRecHitMatrixFlag",      treeVars.ecalRecHitMatrixFlag, "ecalRecHitMatrixFlag[nEcalRecHits][5][5]/F");
  }
  
  //PRESHOWER VARIABLES
  if(EcalTreeContent::preShowerVariables)
    {    
      chain -> Branch("nPreShowerRecHits",     &treeVars.nPreShowerRecHits,       "nPreShowerRecHits/I");
      chain -> Branch("preShowerRecHitType",   treeVars.preShowerRecHitType,      "preShowerRecHitType[nPreShowerRecHits]/F");
      chain -> Branch("preShowerRecHitEnergy", treeVars.preShowerRecHitEnergy,    "preShowerRecHitEnergy[nPreShowerRecHits]/F");
      chain -> Branch("preShowerRecHitIx",     treeVars.preShowerRecHitIx,        "preShowerRecHitIx[nPreShowerRecHits]/F");
      chain -> Branch("preShowerRecHitIy",     treeVars.preShowerRecHitIy,        "preShowerRecHitIy[nPreShowerRecHits]/F");
      chain -> Branch("preShowerRecHitPlane",  treeVars.preShowerRecHitPlane,     "preShowerRecHitPlane[nPreShowerRecHits]/F");
      chain -> Branch("preShowerRecHitStrip",  treeVars.preShowerRecHitStrip,     "preShowerRecHitStrip[nPreShowerRecHits]/F");
      chain -> Branch("preShowerRecHitTime",   treeVars.preShowerRecHitTime,      "preShowerRecHitTime[nPreShowerRecHits]/F");
      chain -> Branch("preShowerRecHitRawId",  treeVars.preShowerRecHitRawId,     "preShowerRecHitRawId[nPreShowerRecHits]/I");
      
      chain -> Branch("nPreShowerDigis",     &treeVars.nPreShowerDigis,           "nPreShowerDigis/I");
      chain -> Branch("esDigiType",          treeVars.esDigiType,                 "esDigiType[nPreShowerDigis]/F");
      chain -> Branch("esDigiIx",            treeVars.esDigiIx,                   "esDigiIx[nPreShowerDigis]/F");
      chain -> Branch("esDigiIy",            treeVars.esDigiIy,                   "esDigiIy[nPreShowerDigis]/F");
      chain -> Branch("esDigiPlane",         treeVars.esDigiPlane,                "esDigiPlane[nPreShowerDigis]/F");
      chain -> Branch("esDigiStrip",         treeVars.esDigiStrip,                "esDigiStrip[nPreShowerDigis]/F");
      chain -> Branch("esDigiRawId",         treeVars.esDigiRawId,                "esDigiRawId[nPreShowerDigis]/I");
      chain -> Branch("preShowerSample1",    treeVars.preShowerSample1,           "preShowerSample1[nPreShowerDigis]/F");
      chain -> Branch("preShowerSample2",    treeVars.preShowerSample2,           "preShowerSample2[nPreShowerDigis]/F");
      chain -> Branch("preShowerSample3",    treeVars.preShowerSample3,           "preShowerSample3[nPreShowerDigis]/F");
    }  //PRESHOWER VARIABLES  


  //L1 VARIABLES
  if(EcalTreeContent::L1TriggerVariables)
    {
      chain -> Branch("techL1Bit",     treeVars.techL1Bit,    "techL1Bit[64]/I");
      chain -> Branch("algoL1Bit",     treeVars.algoL1Bit,    "algoL1Bit[128]/I");
    }  //L1 VARIABLES
  
  
}



void initializeBranches(TTree* chain, EcalTreeContent& treeVars)
{
  treeVars.BX = 0;
  treeVars.lumiId = 0;
  treeVars.runId = 0;
  treeVars.eventId = 0; 
  treeVars.eventNaiveId = 0; 
  
  
  // ECAL VARIABLES  
  if(EcalTreeContent::ecalVariables)
    {    

      //recHitsVariables
      for(int i = 0; i < 75850; ++i)
	{
	  treeVars.ecalRecHitType    [i] = -9999;
	  treeVars.ecalRecHitEnergy  [i] = -9999.;
	  treeVars.ecalRecHitOutOfTimeEnergy  [i] = -9999.;
	  treeVars.ecalRecHitIEta    [i] = -9999.;
	  treeVars.ecalRecHitIPhi    [i] = -9999.;
	  treeVars.ecalRecHitTime    [i] = -9999.;
          treeVars.ecalRecHitChi2    [i] = -9999;
          treeVars.ecalRecHitOutOfTimeChi2[i] = -9999;
	  treeVars.ecalRecHitRawId   [i] = -9999;
	  treeVars.ecalRecHitCoeff   [i] = -9999;
	  treeVars.ecalRecHitRecoFlag[i] = -9999;
	  treeVars.ecalRecHitR9      [i] = -9999;
	  treeVars.ecalRecHitS4oS1   [i] = -9999;
	  treeVars.ecalRecHitIso03   [i][0] = -9999;
	  treeVars.ecalRecHitIso03   [i][1] = -9999;
	  treeVars.ecalRecHitIso04   [i][0] = -9999;
	  treeVars.ecalRecHitIso04   [i][1] = -9999;

	  for (int xx = 0; xx < 5; ++xx)
	    for (int yy = 0; yy < 5; ++yy)
	      {
		treeVars.ecalRecHitMatrix[i][xx][yy] = -9999.;
		treeVars.ecalRecHitMatrixFlag[i][xx][yy] = -9999.;
	      }


	  for (int isample = 0; isample < 10 ; isample++){
	    treeVars.ecalDigis       [i][isample] = -9999;
	    treeVars.ecalGainId      [i][isample] = -9999;
	  }
	}
      
      treeVars.nEcalRecHits = 0;
    } // ECAL VARIABLES
  

  // PRESHOWER VARIABLES
  if(EcalTreeContent::preShowerVariables)
    {
      treeVars.nPreShowerRecHits = 0;
      treeVars.nPreShowerDigis   = 0;

      for(int i = 0; i < 144385; ++i)
	{
	  treeVars.preShowerRecHitType[i] = -9999.;
	  treeVars.preShowerRecHitEnergy[i] = -9999.;
	  treeVars.preShowerRecHitIx[i] = -9999.;
	  treeVars.preShowerRecHitIy[i] = -9999.;
	  treeVars.preShowerRecHitPlane[i] = -9999.;
	  treeVars.preShowerRecHitStrip[i] = -9999.;
	  treeVars.preShowerRecHitTime[i] = -9999.;
	  treeVars.preShowerRecHitRawId[i] = -9999;

	  treeVars.esDigiType[i] = -9999.;
	  treeVars.esDigiIx[i] = -9999.;
	  treeVars.esDigiIy[i] = -9999.;
	  treeVars.esDigiPlane[i] = -9999.;
	  treeVars.esDigiStrip[i] = -9999.;
	  treeVars.esDigiRawId[i] = -9999;
	  treeVars.preShowerSample1[i] = -9999.;
	  treeVars.preShowerSample2[i] = -9999.;
	  treeVars.preShowerSample3[i] = -9999.;
	}
    }
  
  
  //L1 VARIABLES
  if(EcalTreeContent::L1TriggerVariables)
    {
      for (int i = 0; i < 64 ; i++){
	treeVars.techL1Bit[i] = -9999;
      }
      
      for (int i = 0; i < 128 ; i++){
	treeVars.algoL1Bit[i] = -9999;
      }
    }  //L1 VARIABLES
  
  
  

}
