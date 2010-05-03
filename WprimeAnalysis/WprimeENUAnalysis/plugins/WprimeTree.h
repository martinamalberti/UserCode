#ifndef WprimeTree_h
#define WprimeTree_h

// -*- C++ -*-
//
// Package:   ESCosmicsTree
// Class:     ESCosmicsTree
//
/**\class ESCosmicsTree ESCosmicsTree.h

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Original Authors:  Federico DE GUIO Martina MALBERTI
//         Created:  Mo Jul 14 5:46:22 CEST 2008
// $Id: ESCosmicsTree.h,v 1.11 2009/08/31 16:35:59 abenagli Exp $
//
//

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

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"


#include "TFile.h"
#include "TTree.h"

#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeTreeContent.h"

using namespace cms ;
using namespace edm ;
using namespace std ;
using namespace reco;

class WprimeTree : public edm::EDAnalyzer 
{
   
   public:

      explicit WprimeTree (const edm::ParameterSet&) ;
      ~WprimeTree () ;

      
   private:  
      virtual void beginJob ();
      virtual void analyze (const edm::Event&, const edm::EventSetup&) ;
      virtual void endJob () ;

    
   protected:
 
     
     
      void dumpL1Info(edm::Handle<L1GlobalTriggerReadoutRecord>  gtRecord,
		      WprimeTreeContent & myTreeVariables) ;

      void dumpHLTInfo(Handle<TriggerResults> hltresults,
		       WprimeTreeContent & myTreeVariables) ;

     
      void dumpElectronInfo(  View<pat::Electron>  electrons, WprimeTreeContent & myTreeVariables_) ;
      void dumpMetInfo     (  View<pat::MET>            mets, WprimeTreeContent & myTreeVariables_) ;
      void dumpJetInfo     (  View<pat::Jet>            jets, WprimeTreeContent & myTreeVariables_) ;
      void dumpMuonInfo    (  View<pat::Muon>          muons, WprimeTreeContent & myTreeVariables_) ;
      

        


      // ----------member data ---------------------------
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
