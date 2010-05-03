#ifndef WprimeTreeContent_h
#define WprimeTreeContent_h

#include "TChain.h" 

#define MAXELECTRONS 10
#define MAXMUONS 10
#define MAXJETS 100
#define MAXTECHL1BITS 64
#define MAXALGOL1BITS 128


struct WprimeTreeContent
{
  // Flags
  static bool electronVariables;
  static bool metVariables;
  static bool jetVariables;
  static bool muonVariables;
  static bool HLTrigVariables;
  static bool L1TrigVariables;
  
  unsigned int BX;
  unsigned int lumiId;
  unsigned int runId;
  unsigned int eventId;
  unsigned int eventNaiveId;
 

  // electron variables
  int nElectrons;
  float elePx[MAXELECTRONS];
  float elePy[MAXELECTRONS];
  float elePz[MAXELECTRONS];
  float eleE[MAXELECTRONS];
  float eleEt[MAXELECTRONS];
  float eleEta[MAXELECTRONS];
  float elePhi[MAXELECTRONS];
  int   eleId[MAXELECTRONS];
  float eleSigmaIEtaIEta[MAXELECTRONS];
  float eleE1x5[MAXELECTRONS];
  float eleE2x5[MAXELECTRONS];
  float eleE5x5[MAXELECTRONS];
  float eleSeedSwissCross[MAXELECTRONS];
  int   eleCharge[MAXELECTRONS];

  float eleTrkIso[MAXELECTRONS];
  float eleEcalIso[MAXELECTRONS];
  float eleHcalIsoD1[MAXELECTRONS];
  float eleHcalIsoD2[MAXELECTRONS];


  // MET VARIABLES
  float Met;
  float Mex;
  float Mey;
  float MetPhi;

  float uncorrMet;
  float uncorrMex;
  float uncorrMey;
  float uncorrMetPhi;

  // JET VARIABLES
  int nJets;
  float jetPx[MAXJETS];
  float jetPy[MAXJETS];
  float jetPz[MAXJETS];
  float jetPt[MAXJETS];
  float jetEta[MAXJETS];
  float jetPhi[MAXJETS];
  float jetBdisc[MAXJETS];


  // MUON VARIABLES
  int nMuons;
  float muonPx[MAXMUONS];
  float muonPy[MAXMUONS];
  float muonPz[MAXMUONS];
  float muonPt[MAXMUONS];
  float muonEta[MAXMUONS];
  float muonPhi[MAXMUONS];
  

  // HLT VARIABLES
  int HLT_Ele15_SW_L1R;
  int HLT_Ele15_LooseTrackIso_L1R;
  int HLT_Photon15_L1R;
  int HLT_Photon25_L1R;

  // L1 trigger variables
  int techL1Bit[MAXTECHL1BITS];
  int algoL1Bit[MAXALGOL1BITS];


};







// ------------------------------------------------------------------------
//! branch addresses settings

void setBranchAddresses(TTree* chain, WprimeTreeContent& treeVars);






// ------------------------------------------------------------------------
//! create branches for a tree

void setBranches(TTree* chain, WprimeTreeContent& treeVars);






// ------------------------------------------------------------------------
//! initialize branches

void initializeBranches(TTree* chain, WprimeTreeContent& treeVars);



#endif
