#ifndef WprimeAnalyzerPAT_h
#define WprimeAnalyzerPAT_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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
      edm::InputTag hlTriggerResults_ ;
      std::string electronID_;
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

      std::vector<double> * MCelePx_;
      std::vector<double> * MCelePy_;
      std::vector<double> * MCelePz_;
      std::vector<double> * MCeleEta_;
      std::vector<double> * MCelePhi_;
      std::vector<int> * MCelePid_;

      int NeleCand_;
      std::vector<int>  * isMCmatched_;
      std::vector<double> * elePx_;
      std::vector<double> * elePy_;
      std::vector<double> * elePz_;
      std::vector<double> * eleE_;
      std::vector<double> * eleEt_;
      std::vector<double> * eleEta_;
      std::vector<double> * elePhi_;
      std::vector<float>  * eleCharge_;
      std::vector<float>  * eleId_;
      std::vector<double> * eleTrkIsol_;
      std::vector<double> * eleEcalIsol_;
      std::vector<double> * eleHcalIsolD1_;
      std::vector<double> * eleHcalIsolD2_;


      double met_;
      double mex_;
      double mey_;
      double metPhi_;

      double uncorrMet_;
      double uncorrMex_;
      double uncorrMey_;
      double uncorrMetPhi_;

      int Njets_;
      std::vector<double> * jetPx_;
      std::vector<double> * jetPy_;
      std::vector<double> * jetPz_;
      std::vector<double> * jetPt_;
      std::vector<double> * jetEta_;
      std::vector<double> * jetPhi_;

      int HLTLooseIsoEle15_;
      int HLT_EM80_        ;
      int HLT_EM200_       ;
      int HLT_Photon15_       ;
      int HLT_Photon25_       ;


};

#endif
