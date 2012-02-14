// g++ `root-config --libs --cflags` PhotonFix.cc MakeNewTree.cpp -o MakeNewTree

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "TFile.h"
#include "TChain.h"

#include "./PhotonFix.h"

using namespace std;

int main(int argc, char ** argv)
{
 
  //Get old tree
  TChain *ntu = new TChain("ntu");
      
  //---- MC fall 2011
  ntu->Add("../NTUPLES/Fall11/WZAnalysis/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11_All.root");
   
  //---- DATA
//   ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
//   ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
//   ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
//   ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/SingleElectron_Run2011A-WElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
//   ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
  
  //---- Observables
  float pho1E, pho1Eta, pho1Phi, pho1R9;
  float pho2E, pho2Eta, pho2Phi, pho2R9;
  int ele1_isEB;
  int ele2_isEB;

  //---- Set branch addresses for MC  
  ntu->SetBranchAddress("ele1_isEB", &ele1_isEB);

  ntu->SetBranchAddress("ele1_ph_scEta", &pho1Eta);
  ntu->SetBranchAddress("ele1_ph_scPhi", &pho1Phi);
  ntu->SetBranchAddress("ele1_ph_R9", &pho1R9);
  ntu->SetBranchAddress("ele2_isEB", &ele2_isEB);
  ntu->SetBranchAddress("ele2_ph_E", &pho2E);
  ntu->SetBranchAddress("ele2_ph_scEta", &pho2Eta);
  ntu->SetBranchAddress("ele2_ph_scPhi", &pho2Phi);
  ntu->SetBranchAddress("ele2_ph_R9", &pho2R9);
  

  TTree *oldtree = (TTree*)ntu;

  Long64_t nentries = oldtree->GetEntries();
  cout<< "Number of entries in the tree : " << nentries << endl;
  
  //Create a new file + a clone of old tree in new file
  char fname[100];
  sprintf(fname,"WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11_All_new.root");
  TFile *newfile = new TFile(fname,"recreate");
  TTree *newtree = oldtree->CloneTree(0);
  
  float ele1_etaC_DK;
  float ele1_etaS_DK;
  float ele1_etaM_DK;
  float ele1_phiC_DK;
  float ele1_phiS_DK;
  float ele1_phiM_DK;

  float ele1_xC_DK;
  float ele1_xS_DK;
  float ele1_xM_DK;
  float ele1_yC_DK;
  float ele1_yS_DK;
  float ele1_yM_DK;

  float ele2_etaC_DK;
  float ele2_etaS_DK;
  float ele2_etaM_DK;
  float ele2_phiC_DK;
  float ele2_phiS_DK;
  float ele2_phiM_DK;

  float ele2_xC_DK;
  float ele2_xS_DK;
  float ele2_xM_DK;
  float ele2_yC_DK;
  float ele2_yS_DK;
  float ele2_yM_DK;

  newtree->Branch("ele1_etaC_DK",&ele1_etaC_DK,"ele1_etaC_DK/F");
  newtree->Branch("ele1_etaS_DK",&ele1_etaS_DK,"ele1_etaS_DK/F");
  newtree->Branch("ele1_etaM_DK",&ele1_etaM_DK,"ele1_etaM_DK/F");
  newtree->Branch("ele2_etaC_DK",&ele2_etaC_DK,"ele2_etaC_DK/F");
  newtree->Branch("ele2_etaS_DK",&ele2_etaS_DK,"ele2_etaS_DK/F");
  newtree->Branch("ele2_etaM_DK",&ele2_etaM_DK,"ele2_etaM_DK/F");

  newtree->Branch("ele1_phiC_DK",&ele1_phiC_DK,"ele1_phiC_DK/F");
  newtree->Branch("ele1_phiS_DK",&ele1_phiS_DK,"ele1_phiS_DK/F");
  newtree->Branch("ele1_phiM_DK",&ele1_phiM_DK,"ele1_phiM_DK/F");
  newtree->Branch("ele2_phiC_DK",&ele2_phiC_DK,"ele2_phiC_DK/F");
  newtree->Branch("ele2_phiS_DK",&ele2_phiS_DK,"ele2_phiS_DK/F");
  newtree->Branch("ele2_phiM_DK",&ele2_phiM_DK,"ele2_phiM_DK/F");

  newtree->Branch("ele1_xC_DK",&ele1_xC_DK,"ele1_xC_DK/F");
  newtree->Branch("ele1_xS_DK",&ele1_xS_DK,"ele1_xS_DK/F");
  newtree->Branch("ele1_xM_DK",&ele1_xM_DK,"ele1_xM_DK/F");
  newtree->Branch("ele2_xC_DK",&ele2_xC_DK,"ele2_xC_DK/F");
  newtree->Branch("ele2_xS_DK",&ele2_xS_DK,"ele2_xS_DK/F");
  newtree->Branch("ele2_xM_DK",&ele2_xM_DK,"ele2_xM_DK/F");

  newtree->Branch("ele1_yC_DK",&ele1_yC_DK,"ele1_yC_DK/F");
  newtree->Branch("ele1_yS_DK",&ele1_yS_DK,"ele1_yS_DK/F");
  newtree->Branch("ele1_yM_DK",&ele1_yM_DK,"ele1_yM_DK/F");
  newtree->Branch("ele2_yC_DK",&ele2_yC_DK,"ele2_yC_DK/F");
  newtree->Branch("ele2_yS_DK",&ele2_yS_DK,"ele2_yS_DK/F");
  newtree->Branch("ele2_yM_DK",&ele2_yM_DK,"ele2_yM_DK/F");


  PhotonFix::initialise("4_2e");



  for (int i=0; i<nentries; i++) {
    //for (int i=0; i<100; i++) {
    if (i%5000==0) cout << i << endl;
    oldtree->GetEntry(i);
    
    ele1_etaC_DK=-999;
    ele1_etaS_DK=-999;
    ele1_etaM_DK=-999;
    ele1_phiC_DK=-999;
    ele1_phiS_DK=-999;
    ele1_phiM_DK=-999;
    
    ele1_xC_DK=-999;
    ele1_xS_DK=-999;
    ele1_xM_DK=-999;
    ele1_yC_DK=-999;
    ele1_yS_DK=-999;
    ele1_yM_DK=-999;
    
    ele2_etaC_DK=-999;
    ele2_etaS_DK=-999;
    ele2_etaM_DK=-999;
    ele2_phiC_DK=-999;
    ele2_phiS_DK=-999;
    ele2_phiM_DK=-999;
    
    ele2_xC_DK=-999;
    ele2_xS_DK=-999;
    ele2_xM_DK=-999;
    ele2_yC_DK=-999;
    ele2_yS_DK=-999;
    ele2_yM_DK=-999;

    PhotonFix fix1 (pho1E, pho1Eta, pho1Phi, pho1R9);
    
    if (fabs(fix1.eta())<1.48){
      ele1_etaC_DK = fix1.etaC();
      ele1_etaS_DK = fix1.etaS();
      ele1_etaM_DK = fix1.etaM();
      ele1_phiC_DK = fix1.phiC();
      ele1_phiS_DK = fix1.phiS();
      ele1_phiM_DK = fix1.phiM();
    }

    else{
      ele1_xC_DK = fix1.xC();
      ele1_xS_DK = fix1.xS();
      ele1_xM_DK = fix1.xM();
      ele1_yC_DK = fix1.yC();
      ele1_yS_DK = fix1.yS();
      ele1_yM_DK = fix1.yM();
    }

    
    PhotonFix fix2 (pho2E, pho2Eta, pho2Phi, pho2R9);
   
    if (fabs(fix2.eta())<1.48){
      ele2_etaC_DK = fix2.etaC();
      ele2_etaS_DK = fix2.etaS();
      ele2_etaM_DK = fix2.etaM();
      ele2_phiC_DK = fix2.phiC();
      ele2_phiS_DK = fix2.phiS();
      ele2_phiM_DK = fix2.phiM();
    }
    else{
      ele2_xC_DK = fix2.xC();
      ele2_xS_DK = fix2.xS();
      ele2_xM_DK = fix2.xM();
      ele2_yC_DK = fix2.yC();
      ele2_yS_DK = fix2.yS();
      ele2_yM_DK = fix2.yM();
    }
    newtree->Fill();
    
  }
  
  newtree->Print();
  newtree->AutoSave();
  
  delete newfile;


}
