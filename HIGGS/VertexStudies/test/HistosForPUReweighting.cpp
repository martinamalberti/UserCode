///==== test program ====

#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "Math/GenVector/VectorUtil.h"
#include "TRandom3.h"
#include <time.h>
#include <sstream>
#include "MyTest.h"

#include <iostream>

#include "TClonesArray.h"
#include "TMatrix.h"

#define etaEB 1.444
#define etaEE 1.560


int main(int argc, char** argv)
{ 

  //Chain
  TChain* chain = new TChain("MiBiCommonNTTwoPhotons/SimpleNtuple");
  // Summer 12 DY sample
  chain->Add("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VERTEX/2012/MC/MiBiCommonNT_1*.root");
  //chain->Add("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VERTEX/2012/MC/MiBiCommonNT_2*.root");



  treeReader reader((TTree*)(chain));

  // histos
  TH1F *hgen = new TH1F("hgen","hgen",60,0,60);
  hgen  -> GetXaxis()->SetTitle("N_{gen} + 1");

  TH1F *hnpumc = new TH1F("hnpumc","hnpumc",60,0,60);
  hnpumc  -> GetXaxis()->SetTitle("N_{gen}");

  TH1F *htruenpumc = new TH1F("htruenpumc","htruenpumc",60,0,60);
  htruenpumc  -> GetXaxis()->SetTitle("N_{gen}");

  TH2F *hnvtxreco_vs_npugen     = new TH2F("hnvtxreco_vs_npugen","hnvtxreco_vs_npugen",60,0,60,60,0,60);
  hnvtxreco_vs_npugen -> GetXaxis()->SetTitle("N_{gen} + 1");
  hnvtxreco_vs_npugen -> GetYaxis()->SetTitle("N_{reco}");

  TH2F *hnpugen_vs_nvtxreco     = new TH2F("hnpugen_vs_nvtxreco","hnpugen_vs_nvtxreco",60,0,60,60,0,60);
  hnpugen_vs_nvtxreco -> GetXaxis()->SetTitle("N_{reco}"); 
  hnpugen_vs_nvtxreco -> GetYaxis()->SetTitle("N_{gen} + 1");
 


  int entryMax = reader.GetEntries();
  std::cout<< "Total number of entries : " << entryMax << std::endl;
  
  // estimate the response curve NPU vs NVXT
  for (int u = 0; u < entryMax; u++ )
    {
      if(u%10000 == 0) std::cout<<"reading event "<< u <<std::endl;
      reader.GetEntry(u);
      
     
      std::vector<int>* mc_PUit_NumInteractions  = reader.GetInt("mc_PUit_NumInteractions");
      int npu = mc_PUit_NumInteractions->at(0);

      std::vector<float>* mc_PUit_TrueNumInteractions  = reader.GetFloat("mc_PUit_TrueNumInteractions");
      int truenpu = mc_PUit_TrueNumInteractions->at(0);

      std::vector<float>* PV_z       = reader.GetFloat("PV_z");
      int nvtx = PV_z->size(); // number of sim PU vtx
      
      hgen->Fill(npu+1);
      hnpumc->Fill(npu);
      htruenpumc->Fill(truenpu);

      hnvtxreco_vs_npugen -> Fill(npu+1,nvtx);
      hnpugen_vs_nvtxreco -> Fill(nvtx,npu+1);
  
    }
  
  TFile ff("weights/HistosForPUReweighting_DYJetsToLL_Summer12_DR53X_PU-S10.root","recreate");

  hgen->Write();
  hnpumc->Write();
  htruenpumc->Write();
  hnvtxreco_vs_npugen -> Write();
  hnpugen_vs_nvtxreco -> Write();

  ff.Close();
  return 0;
  

}
