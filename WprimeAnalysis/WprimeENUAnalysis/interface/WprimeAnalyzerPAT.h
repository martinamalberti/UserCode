#ifndef WprimeAnalyzerPAT_h
#define WprimeAnalyzerPAT_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


// ROOT include
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"



//
// class decleration
//

class WprimeAnalyzerPAT : public edm::EDAnalyzer {
   public:
      explicit WprimeAnalyzerPAT(const edm::ParameterSet&);
      ~WprimeAnalyzerPAT();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


      // ----------member data ---------------------------
      edm::InputTag eleLabel_ ;
      edm::InputTag metLabel_ ;
      edm::InputTag jetLabel_ ;
      edm::InputTag muonLabel_ ;
      edm::InputTag hlTriggerResults_ ;
      edm::InputTag L1GlobalTrigger_  ;
      std::string electronID_;
      std::string btagAlgo_;
      edm::InputTag  TrkIsolationProducer_ ;
      edm::InputTag  EcalIsolationProducer_ ;
      edm::InputTag  HcalIsolationProducerDepth1_ ;
      edm::InputTag  HcalIsolationProducerDepth2_ ;

    
     
      
      std::vector<std::string>  hlNames_;  // name of each HLT algorithm


      TFile *thefile;
      int EvtCounter;

      TH1F *hNelectronCandidates;

      TH1F *hTrkIsolation;
      TH1F *hEcalIsolation;
      TH1F *hHcalIsolationDepth1;
      TH1F *hHcalIsolationDepth2;

      TH1F *hEMPlusHcalDepth1;


      TTree * tTreeUtilities_;

      std::vector<double> * MCelePx;
      std::vector<double> * MCelePy;
      std::vector<double> * MCelePz;
      std::vector<double> * MCeleEta;
      std::vector<double> * MCelePhi;
      std::vector<int>    * MCelePid;

      std::vector<double> * MCmuonPx;
      std::vector<double> * MCmuonPy;
      std::vector<double> * MCmuonPz;
      std::vector<double> * MCmuonEta;
      std::vector<double> * MCmuonPhi;
      std::vector<int>    * MCmuonPid;

      int NeleCand;
      std::vector<int>  * isMCmatched;
      std::vector<double> * elePx;
      std::vector<double> * elePy;
      std::vector<double> * elePz;
      std::vector<double> * eleE;
      std::vector<double> * eleEt;
      std::vector<double> * eleEta;
      std::vector<double> * elePhi;
      std::vector<double>  * eleCharge;
      std::vector<double>  * eleId;
      std::vector<double>  * eleSigmaIEtaIEta;
      std::vector<double>  * eleE1x5;
      std::vector<double>  * eleE2x5;
      std::vector<double>  * eleE5x5;

      std::vector<double> * eleTrkIsol;
      std::vector<double> * eleEcalIsol;
      std::vector<double> * eleHcalIsolD1;
      std::vector<double> * eleHcalIsolD2;


      double Met;
      double Mex;
      double Mey;
      double MetPhi;

      double uncorrMet;
      double uncorrMex;
      double uncorrMey;
      double uncorrMetPhi;

      int Njets;
      std::vector<double> * jetPx;
      std::vector<double> * jetPy;
      std::vector<double> * jetPz;
      std::vector<double> * jetPt;
      std::vector<double> * jetEta;
      std::vector<double> * jetPhi;
      std::vector<double> * jetBdisc;


      std::vector<double> * muonPt;
      std::vector<double> * muonEta;
      std::vector<double> * muonPhi;

      int HLTEle15;
      int HLTEle10;  
      int HLTEle20;

      int HLTLooseIsoEle15; 
    
      int HLTPhoton15;
      int HLTPhoton25;


};

#endif
