#ifndef EcalTree_h
#define EcalTree_h

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
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "TFile.h"
#include "TTree.h"



#include "EcalAnalysis/SpikeStudies/interface/EcalTreeContent.h"

class EcalTree : public edm::EDAnalyzer 
{
   
   public:

      explicit EcalTree (const edm::ParameterSet&) ;
      ~EcalTree () ;

      
   private:  
      virtual void beginJob ();
      virtual void analyze (const edm::Event&, const edm::EventSetup&) ;
      virtual void endJob () ;

    
   protected:
 
     
      /*       ! dump  information */
      void dumpL1Info(edm::Handle<L1GlobalTriggerReadoutRecord>  gtRecord,
		      EcalTreeContent & myTreeVariables) ;

      //! dump  information
      void dumpBarrelInfo( const CaloTopology* theCaloTopology,
			   const CaloGeometry* theCaloGeometry,
			   const EBDigiCollection* theEcalBarrelDigis,
			   const EcalRecHitCollection* theBarrelEcalRecHits,
			   const reco::BasicClusterCollection* ebClusters,
			   EcalTreeContent & myTreeVariables_) ;
      

        
      edm::InputTag ebRecHitCollection_ ;
      edm::InputTag ebClusterCollection_ ;
      edm::InputTag ebDigiCollection_;
      edm::InputTag L1InputTag_;


      double radiusForIso_;
      double energyCutForIso_;
    
      int naiveId_;
    
      EcalTreeContent myTreeVariables_ ;

      TTree* tree_ ;

     
} ;

#endif
