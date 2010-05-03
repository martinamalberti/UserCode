#ifndef WprimeTreeHEEP_h
#define WprimeTreeHEEP_h

// system include files
#include <memory>
#include <vector>
#include <map>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"


#include "TFile.h"
#include "TTree.h"

#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeTreeContent.h"

using namespace cms ;
using namespace edm ;
using namespace std ;
using namespace reco;

class WprimeTreeHEEP : public edm::EDAnalyzer 
{
   
   public:

      explicit WprimeTreeHEEP (const edm::ParameterSet&) ;
      ~WprimeTreeHEEP () ;

      
   private:  
      virtual void beginJob ();
      virtual void analyze (const edm::Event&, const edm::EventSetup&) ;
      virtual void endJob () ;

    
   protected:
 
     
     
      void dumpL1Info(edm::Handle<L1GlobalTriggerReadoutRecord>  gtRecord,
		      WprimeTreeContent & myTreeVariables) ;

      void dumpHLTInfo(Handle<TriggerResults> hltresults,
		       WprimeTreeContent & myTreeVariables) ;

     
      //      void dumpElectronInfo(  View<pat::Electron>  electrons, WprimeTreeContent & myTreeVariables_) ;
      void dumpElectronInfo( View<pat::Electron>  electrons, const EcalRecHitCollection* theBarrelEcalRecHits,
			     const CaloTopology* theCaloTopology, WprimeTreeContent & myTreeVariables_) ;
      void dumpMetInfo     ( View<pat::MET>            mets, WprimeTreeContent & myTreeVariables_) ;
      void dumpJetInfo     ( View<pat::Jet>            jets, WprimeTreeContent & myTreeVariables_) ;
      void dumpMuonInfo    ( View<pat::Muon>          muons, WprimeTreeContent & myTreeVariables_) ;
      

        


      // ----------member data ---------------------------
      edm::InputTag recHitCollection_EB_;
      edm::InputTag eleLabel_ ;
      edm::InputTag metLabel_ ;
      edm::InputTag jetLabel_ ;
      edm::InputTag muonLabel_ ;
      std::string electronID_;
      std::string btagAlgo_;
      edm::InputTag  TrkIsolationProducer_ ;
      edm::InputTag  EcalIsolationProducer_ ;
      edm::InputTag  HcalIsolationProducerDepth1_ ;
      edm::InputTag  HcalIsolationProducerDepth2_ ;
      
      edm::InputTag L1InputTag_;
      edm::InputTag HLTInputTag_ ;
     
      int naiveId_;
    
      WprimeTreeContent myTreeVariables_ ;

      TTree* tree_ ;

     
} ;

#endif
