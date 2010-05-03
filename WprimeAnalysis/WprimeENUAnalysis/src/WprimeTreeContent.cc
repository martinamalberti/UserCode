#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeTreeContent.h"

bool WprimeTreeContent::electronVariables      = true;
bool WprimeTreeContent::metVariables           = true;
bool WprimeTreeContent::jetVariables           = true;
bool WprimeTreeContent::muonVariables          = true;
bool WprimeTreeContent::HLTrigVariables        = true;
bool WprimeTreeContent::L1TrigVariables        = true;
  
void setBranchAddresses(TTree* chain, WprimeTreeContent& treeVars)
{
  chain -> SetBranchAddress("BX",            &treeVars.BX);
  chain -> SetBranchAddress("lumiId",        &treeVars.lumiId);
  chain -> SetBranchAddress("runId",         &treeVars.runId);
  chain -> SetBranchAddress("eventId",       &treeVars.eventId);
  chain -> SetBranchAddress("eventNaiveId",  &treeVars.eventNaiveId);


  // ELECTRON VARIABLES  
  if(WprimeTreeContent::electronVariables)
    {  
      
      chain -> SetBranchAddress("nElectrons",             &treeVars.nElectrons);
      chain -> SetBranchAddress("elePx",                   treeVars.elePx);
      chain -> SetBranchAddress("elePy",                   treeVars.elePy);
      chain -> SetBranchAddress("elePz",                   treeVars.elePz);
      chain -> SetBranchAddress("eleE",                    treeVars.eleE);
      chain -> SetBranchAddress("eleEt",                   treeVars.eleEt);
      chain -> SetBranchAddress("eleEta",                  treeVars.eleEta);
      chain -> SetBranchAddress("elePhi",                  treeVars.elePhi);

      chain -> SetBranchAddress("eleId",                   treeVars.eleId);
      chain -> SetBranchAddress("eleSigmaIEtaIEta",        treeVars.eleSigmaIEtaIEta);
      chain -> SetBranchAddress("eleE1x5",                 treeVars.eleE1x5);
      chain -> SetBranchAddress("eleE2x5",                 treeVars.eleE2x5);
      chain -> SetBranchAddress("eleE5x5",                 treeVars.eleE5x5);
      chain -> SetBranchAddress("eleSeedSwissCross",       treeVars.eleSeedSwissCross);

      chain -> SetBranchAddress("eleCharge",               treeVars.eleCharge);
      chain -> SetBranchAddress("eleTrkIso",               treeVars.eleTrkIso);
      chain -> SetBranchAddress("eleEcalIso",              treeVars.eleEcalIso);
      chain -> SetBranchAddress("eleHcalIsoD1",            treeVars.eleHcalIsoD1);
      chain -> SetBranchAddress("eleHcalIsoD2",            treeVars.eleHcalIsoD2);

    } // ELECTRON VARIABLES
  
 
  // MET VARIABLES  
  if(WprimeTreeContent::metVariables)
    {  
      chain -> SetBranchAddress("Met",                   &treeVars.Met);
      chain -> SetBranchAddress("Mex",                   &treeVars.Mex);
      chain -> SetBranchAddress("Mey",                   &treeVars.Mey);
      chain -> SetBranchAddress("MetPhi",                &treeVars.MetPhi);
      
      chain -> SetBranchAddress("uncorrMet",             &treeVars.uncorrMet);
      chain -> SetBranchAddress("uncorrMex",             &treeVars.uncorrMex);
      chain -> SetBranchAddress("uncorrMey",             &treeVars.uncorrMey);
      chain -> SetBranchAddress("uncorrMetPhi",          &treeVars.uncorrMetPhi);
      
    } // MET VARIABLES
  

  // JET VARIABLES  
  if(WprimeTreeContent::jetVariables)
    {  
      
      chain -> SetBranchAddress("nJets",             &treeVars.nJets);
      chain -> SetBranchAddress("jetPx",              treeVars.jetPx);
      chain -> SetBranchAddress("jetPy",              treeVars.jetPy);
      chain -> SetBranchAddress("jetPz",              treeVars.jetPz);
      chain -> SetBranchAddress("jetPt",              treeVars.jetPt);
      chain -> SetBranchAddress("jetEta",             treeVars.jetEta);
      chain -> SetBranchAddress("jetPhi",             treeVars.jetPhi);
      chain -> SetBranchAddress("jetBdisc",           treeVars.jetBdisc);
      
    } // JET VARIABLES
  


  // MUON VARIABLES  
  if(WprimeTreeContent::muonVariables)
    {  
      
      chain -> SetBranchAddress("nMuons",             &treeVars.nMuons);
      chain -> SetBranchAddress("muonPx",              treeVars.muonPx);
      chain -> SetBranchAddress("muonPy",              treeVars.muonPy);
      chain -> SetBranchAddress("muonPz",              treeVars.muonPz);
      chain -> SetBranchAddress("muonPt",              treeVars.muonPt);
      chain -> SetBranchAddress("muonEta",             treeVars.muonEta);
      chain -> SetBranchAddress("muonPhi",             treeVars.muonPhi);
            
    } // MUON VARIABLES



  //L1 VARIABLES
  if(WprimeTreeContent::L1TrigVariables)
    {
      chain -> SetBranchAddress("techL1Bit",     treeVars.techL1Bit);
      chain -> SetBranchAddress("algoL1Bit",     treeVars.algoL1Bit);
    }  //L1 VARIABLES

  
  //HLT VARIABLES
  if(WprimeTreeContent::HLTrigVariables)
    {
      chain -> SetBranchAddress("HLT_Ele15_SW_L1R",           &treeVars.HLT_Ele15_SW_L1R);
      chain -> SetBranchAddress("HLT_Ele15_LooseTrackIso_L1R",&treeVars.HLT_Ele15_LooseTrackIso_L1R);
      chain -> SetBranchAddress("HLT_Photon15_L1R",           &treeVars.HLT_Photon15_L1R);
      chain -> SetBranchAddress("HLT_Photon25_L1R",           &treeVars.HLT_Photon25_L1R);

    }//HLT VARIABLES

}


 



void setBranches(TTree* chain, WprimeTreeContent& treeVars)
{
  chain -> Branch("BX",            &treeVars.BX,                       "BX/i");
  chain -> Branch("lumiId",        &treeVars.lumiId,               "lumiId/i");
  chain -> Branch("runId",         &treeVars.runId,                 "runId/i");
  chain -> Branch("eventId",       &treeVars.eventId,             "eventId/i");
  chain -> Branch("eventNaiveId",  &treeVars.eventNaiveId,   "eventNaiveId/i");

  
  // ELECTRON  VARIABLES  
  if(WprimeTreeContent::electronVariables)
    {
      
      chain -> Branch("nElectrons",        &treeVars.nElectrons,       "nElectrons/I");
      chain -> Branch("elePx",              treeVars.elePx,            "elePx[nElectrons]/F");
      chain -> Branch("elePy",              treeVars.elePy,            "elePy[nElectrons]/F");
      chain -> Branch("elePz",              treeVars.elePz,            "elePz[nElectrons]/F");
      chain -> Branch("eleE",               treeVars.eleE,             "eleE[nElectrons]/F");
      chain -> Branch("eleEt",              treeVars.eleEt,            "eleEt[nElectrons]/F");
      chain -> Branch("eleEta",             treeVars.eleEta,           "eleEta[nElectrons]/F");
      chain -> Branch("elePhi",             treeVars.elePhi,           "elePhi[nElectrons]/F");
      
      chain -> Branch("eleId",              treeVars.eleId,            "eleId[nElectrons]/I");
      chain -> Branch("eleSigmaIEtaIEta",   treeVars.eleSigmaIEtaIEta, "eleSigmaIEtaIEta[nElectrons]/F");
      chain -> Branch("eleE1x5",            treeVars.eleE1x5,          "eleE1xE5[nElectrons]/F");
      chain -> Branch("eleE2x5",            treeVars.eleE2x5,          "eleE2xE5[nElectrons]/F");
      chain -> Branch("eleE5x5",            treeVars.eleE5x5,          "eleE5xE5[nElectrons]/F");
      chain -> Branch("eleSeedSwissCross",  treeVars.eleSeedSwissCross,"eleSeedSwisscross[nElectrons]/F");
      
      chain -> Branch("eleCharge",          treeVars.eleCharge,       "eleCharge[nElectrons]/I");
      chain -> Branch("eleTrkIso",          treeVars.eleTrkIso,       "eleTrkIso[nElectrons]/F");
      chain -> Branch("eleEcalIso",         treeVars.eleEcalIso,      "eleEcalIso[nElectrons]/F");
      chain -> Branch("eleHcalIsoD1",       treeVars.eleHcalIsoD1,    "eleHcalIsoD1[nElectrons]/F");
      chain -> Branch("eleHcalIsoD2",       treeVars.eleHcalIsoD2,    "eleHcalIsoD2[nElectrons]/F");
  
    }
  
 

  // MET VARIABLES  
  if(WprimeTreeContent::metVariables)
    {  
      chain -> Branch("Met",          &treeVars.Met,         "Met/F");
      chain -> Branch("Mex",          &treeVars.Mex,         "Mex/F");
      chain -> Branch("Mey",          &treeVars.Mey,         "Mey/F");
      chain -> Branch("MetPhi",       &treeVars.MetPhi,      "MetPhi/F");
      
      chain -> Branch("uncorrMet",    &treeVars.uncorrMet,   "uncorrMet/F");
      chain -> Branch("uncorrMex",    &treeVars.uncorrMex,   "uncorrMex/F"); 
      chain -> Branch("uncorrMey",    &treeVars.uncorrMey,   "uncorrMey/F");
      chain -> Branch("uncorrMetPhi", &treeVars.uncorrMetPhi,"uncorrMetPhi/F");
      
    } // MET VARIABLES
  

  // JET VARIABLES  
  if(WprimeTreeContent::jetVariables)
    {  
      
      chain -> Branch("nJets",        &treeVars.nJets,       "nJets/I");
      chain -> Branch("jetPx",        treeVars.jetPx,        "jetPx[nJets]/F");
      chain -> Branch("jetPy",        treeVars.jetPy,        "jetPy[nJets]/F");
      chain -> Branch("jetPz",        treeVars.jetPz,        "jetPz[nJets]/F");
      chain -> Branch("jetPt",        treeVars.jetPt,        "jetPt[nJets]/F");
      chain -> Branch("jetEta",       treeVars.jetEta,       "jetEta[nJets]/F");
      chain -> Branch("jetPhi",       treeVars.jetPhi,       "jetPhi[nJets]/F");
      chain -> Branch("jetBdisc",     treeVars.jetBdisc,     "jetBdisc[nJets]/F");
      
    } // JET VARIABLES
  


  // MUON VARIABLES  
  if(WprimeTreeContent::muonVariables)
    {  
      
      chain -> Branch("nMuons",     &treeVars.nMuons,         "nMuons/I");
      chain -> Branch("muonPx",      treeVars.muonPx,         "muonPx[nMuons]/F");
      chain -> Branch("muonPy",      treeVars.muonPy,         "muonPy[nMuons]/F");
      chain -> Branch("muonPz",      treeVars.muonPz,         "muonPz[nMuons]/F");
      chain -> Branch("muonPt",      treeVars.muonPt,         "muonPt[nMuons]/F");
      chain -> Branch("muonEta",     treeVars.muonEta,        "muonEta[nMuons]/F");
      chain -> Branch("muonPhi",     treeVars.muonPhi,        "muonPhi[nMuons]/F");
            
    } // MUON VARIABLES



 

  //L1 VARIABLES
  if(WprimeTreeContent::L1TrigVariables)
    {
      chain -> Branch("techL1Bit",     treeVars.techL1Bit,    "techL1Bit[64]/I");
      chain -> Branch("algoL1Bit",     treeVars.algoL1Bit,    "algoL1Bit[128]/I");
    }  //L1 VARIABLES
  



  //HLT VARIABLES
  if(WprimeTreeContent::HLTrigVariables)
    {
      chain -> Branch("HLT_Ele15_SW_L1R",           &treeVars.HLT_Ele15_SW_L1R, "HLT_Ele15_SW_L1R/I");
      chain -> Branch("HLT_Ele15_LooseTrackIso_L1R",&treeVars.HLT_Ele15_LooseTrackIso_L1R,"HLT_Ele15_LooseTrackIso_L1R/I");
      chain -> Branch("HLT_Photon15_L1R",           &treeVars.HLT_Photon15_L1R, "HLT_Photon15_L1R/I");
      chain -> Branch("HLT_Photon25_L1R",           &treeVars.HLT_Photon25_L1R, "HLT_Photon25_L1R/I");

    }//HLT VARIABLES
  
}



void initializeBranches(TTree* chain, WprimeTreeContent& treeVars)
{
  treeVars.BX = 0;
  treeVars.lumiId = 0;
  treeVars.runId = 0;
  treeVars.eventId = 0; 
  treeVars.eventNaiveId = 0; 
  
  
  // ELECTRONS VARIABLES  
  if(WprimeTreeContent::electronVariables)
    {    
      for(int i = 0; i < MAXELECTRONS; ++i)
	{
	  treeVars.elePx[i] = -9999;
	  treeVars.elePy[i] = -9999;
	  treeVars.elePz[i] = -9999;
	  treeVars.eleE[i] = -9999;
	  treeVars.eleEt[i] = -9999;
	  treeVars.eleEta[i] = -9999;
	  treeVars.elePhi[i] = -9999;
	  treeVars.eleId[i] = -9999;
	  treeVars.eleSigmaIEtaIEta[i] = -9999;
	  treeVars.eleE1x5[i] = -9999;
	  treeVars.eleE2x5[i] = -9999;
	  treeVars.eleE5x5[i] = -9999;
	  treeVars.eleSeedSwissCross[i] = -9999;
 
	  treeVars.eleCharge[i] = -9999;
	  treeVars.eleTrkIso[i] = -9999;
	  treeVars.eleEcalIso[i] = -9999;
	  treeVars.eleHcalIsoD1[i] = -9999;
	  treeVars.eleHcalIsoD2[i] = -9999;

	}
      
      treeVars.nElectrons = 0;
    } // ELECTRONS VARIABLES
  



 // MET VARIABLES  
  if(WprimeTreeContent::metVariables)
    {  
      treeVars.Met = -9999;
      treeVars.Mex = -9999;
      treeVars.Mey = -9999;
      treeVars.MetPhi = -9999;
      
      treeVars.uncorrMet = -9999;
      treeVars.uncorrMex = -9999;
      treeVars.uncorrMey = -9999;
      treeVars.uncorrMetPhi = -9999;
      
    } // MET VARIABLES
  

  // JET VARIABLES  
  if(WprimeTreeContent::jetVariables)
    {  
      for(int i = 0; i < MAXJETS; ++i) {
	treeVars.jetPx[i] = -9999;
	treeVars.jetPy[i] = -9999;
	treeVars.jetPz[i] = -9999;
	treeVars.jetPt[i] = -9999;
	treeVars.jetEta[i] = -9999;
	treeVars.jetPhi[i] = -9999;
	treeVars.jetBdisc[i] = -9999;
      }
 
      treeVars.nJets = 0;
      
    } // JET VARIABLES
  


  // MUON VARIABLES  
  if(WprimeTreeContent::muonVariables)
    {  
      for(int i = 0; i < MAXMUONS; ++i) {
	treeVars.muonPx[i] = -9999;
	treeVars.muonPy[i] = -9999;
	treeVars.muonPz[i] = -9999;
	treeVars.muonPt[i] = -9999;
	treeVars.muonEta[i] = -9999;
	treeVars.muonPhi[i] = -9999;
      }
 
      treeVars.nMuons = 0;
     
            
    } // MUON VARIABLES





  
  
  //L1 VARIABLES
  if(WprimeTreeContent::L1TrigVariables)
    {
      for (int i = 0; i < 64 ; i++){
	treeVars.techL1Bit[i] = -9999;
      }
      
      for (int i = 0; i < 128 ; i++){
	treeVars.algoL1Bit[i] = -9999;
      }
    }  //L1 VARIABLES
  
  
    //HLT VARIABLES
  if(WprimeTreeContent::HLTrigVariables)
    {
      treeVars.HLT_Ele15_SW_L1R = 0;
      treeVars.HLT_Ele15_LooseTrackIso_L1R = 0; 

      treeVars.HLT_Photon15_L1R = 0; 
      treeVars.HLT_Photon25_L1R = 0; 


    }//HLT VARIABLES
  

}
