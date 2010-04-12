#ifndef EcalTreeContent_h
#define EcalTreeContent_h

#include "TChain.h" 

#define MAXPRESHOWERRECHITS 144385
#define MAXECALRECHITS 75850
#define MAXTECHL1BITS 64
#define MAXALGOL1BITS 128


struct EcalTreeContent
{
  // Flags
  static bool ecalVariables;
  static bool preShowerVariables;
  static bool L1TriggerVariables;
  
  unsigned int BX;
  unsigned int lumiId;
  unsigned int runId;
  unsigned int eventId;
  unsigned int eventNaiveId;
 
  
  //ecalVariables variables

  int nEcalRecHits;
  float ecalRecHitType[MAXECALRECHITS];
  float ecalRecHitEnergy[MAXECALRECHITS];
  float ecalRecHitOutOfTimeEnergy[MAXECALRECHITS];
  float ecalRecHitIEta[MAXECALRECHITS];
  float ecalRecHitIPhi[MAXECALRECHITS];
  float ecalRecHitTime[MAXECALRECHITS];
  float ecalRecHitChi2[MAXECALRECHITS];
  float ecalRecHitOutOfTimeChi2[MAXECALRECHITS];
  int ecalRecHitRawId[MAXECALRECHITS];
  float ecalRecHitCoeff[MAXECALRECHITS];
  int ecalRecHitRecoFlag[MAXECALRECHITS];
  float ecalRecHitR9[MAXECALRECHITS];
  float ecalRecHitS4oS1[MAXECALRECHITS];
  float ecalRecHitIso03[MAXECALRECHITS][2];
  float ecalRecHitIso04[MAXECALRECHITS][2];
  int ecalDigis[MAXECALRECHITS][10];
  int ecalGainId[MAXECALRECHITS][10];
  float ecalRecHitMatrix[MAXECALRECHITS][5][5];
  float ecalRecHitMatrixFlag[MAXECALRECHITS][5][5];


  //preshower variables
  int nPreShowerRecHits;
  float preShowerRecHitType[MAXPRESHOWERRECHITS];
  float preShowerRecHitEnergy[MAXPRESHOWERRECHITS];
  float preShowerRecHitIx[MAXPRESHOWERRECHITS];
  float preShowerRecHitIy[MAXPRESHOWERRECHITS];
  float preShowerRecHitPlane[MAXPRESHOWERRECHITS];
  float preShowerRecHitStrip[MAXPRESHOWERRECHITS];
  float preShowerRecHitTime[MAXPRESHOWERRECHITS];
  int preShowerRecHitRawId[MAXPRESHOWERRECHITS];

  int nPreShowerDigis; 
  float esDigiType[MAXPRESHOWERRECHITS];
  float esDigiIx[MAXPRESHOWERRECHITS];
  float esDigiIy[MAXPRESHOWERRECHITS];
  float esDigiPlane[MAXPRESHOWERRECHITS];
  float esDigiStrip[MAXPRESHOWERRECHITS];
  int esDigiRawId[MAXPRESHOWERRECHITS];
  float preShowerSample1[MAXPRESHOWERRECHITS];
  float preShowerSample2[MAXPRESHOWERRECHITS];
  float preShowerSample3[MAXPRESHOWERRECHITS];

  // L1 trigger variables
  int techL1Bit[MAXTECHL1BITS];
  int algoL1Bit[MAXALGOL1BITS];


};







// ------------------------------------------------------------------------
//! branch addresses settings

void setBranchAddresses(TTree* chain, EcalTreeContent& treeVars);






// ------------------------------------------------------------------------
//! create branches for a tree

void setBranches(TTree* chain, EcalTreeContent& treeVars);






// ------------------------------------------------------------------------
//! initialize branches

void initializeBranches(TTree* chain, EcalTreeContent& treeVars);



#endif
