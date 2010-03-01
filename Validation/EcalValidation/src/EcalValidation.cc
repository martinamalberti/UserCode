// -*- C++ -*-
//
// Package:    EcalValidation
// Class:      EcalValidation
// 
/**\class EcalValidation EcalValidation.cc Validation/EcalValidation/src/EcalValidation.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Federico Ferri
//         Created:  Fri Mar 21 18:06:59 CET 2008
// $Id: EcalValidation.cc,v 1.9 2010/01/10 12:27:08 ferriff Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"

#include "Validation/EcalValidation/interface/EcalValidation.h"

using namespace cms ;
using namespace edm ;
using namespace std ;
using namespace reco;

//
// constructors and destructor
//
EcalValidation::EcalValidation(const edm::ParameterSet& ps)
{
  //now do what ever initialization is needed
  recHitCollection_EB_       = ps.getParameter<edm::InputTag>("recHitCollection_EB");
  recHitCollection_EE_       = ps.getParameter<edm::InputTag>("recHitCollection_EE");
  basicClusterCollection_EB_ = ps.getParameter<edm::InputTag>("basicClusterCollection_EB");
  basicClusterCollection_EE_ = ps.getParameter<edm::InputTag>("basicClusterCollection_EE");
  superClusterCollection_EB_ = ps.getParameter<edm::InputTag>("superClusterCollection_EB");
  superClusterCollection_EE_ = ps.getParameter<edm::InputTag>("superClusterCollection_EE");
  esRecHitCollection_        = ps.getParameter<edm::InputTag>("recHitCollection_ES");
  esClusterCollectionX_      = ps.getParameter<edm::InputTag>("ClusterCollectionX_ES");
  esClusterCollectionY_      = ps.getParameter<edm::InputTag>("ClusterCollectionY_ES");

  ethrEB_                    = ps.getParameter<double>("ethrEB");
  ethrEE_                    = ps.getParameter<double>("ethrEE");

  edm::Service<TFileService> fs;
  TFileDirectory dirEB = fs->mkdir("EB");
  TFileDirectory dirEE = fs->mkdir("EE");
  TFileDirectory dirES = fs->mkdir("ES");
  
  //dirEB->cd();
  
  // RecHits ---------------------------------------------- 
  // ... barrel
  h_recHits_EB_size   = dirEB.make<TH1D>("h_recHits_EB_size", "h_recHitsEB_size", 1000, 0, 10000 );
  h_recHits_EB_energy = dirEB.make<TH1D>("h_recHits_EB_energy","h_recHitsEB_energy",2000,-50,350);
  h_recHits_EB_eta    = dirEB.make<TH1D>("h_recHits_EB_eta","h_recHits_EB_eta",300,-300.,300.);
  h_recHits_EB_phi    = dirEB.make<TH1D>("h_recHits_EB_phi","h_recHits_EB_phi",320,-3.2,3.2);
  h_recHits_EB_time   = dirEB.make<TH1D>("h_recHits_EB_time","h_recHitsEB_time",400,-100,100);
  h_recHits_EB_Chi2   = dirEB.make<TH1D>("h_recHits_EB_Chi2","h_recHitsEB_Chi2",1000,0,100);
  h_recHits_EB_OutOfTimeChi2 = dirEB.make<TH1D>("h_recHits_EB_OutOfTimeChi2","h_recHitsEB_OutOfTimeChi2",1000,0,100);  
  h_recHits_EB_occupancy = dirEB.make<TH2D>("h_recHits_EB_occupancy","h_recHitsEB_occupancy",172,-86.,86.,360,1.,361. );

  // ... endcap
  h_recHits_EE_eta = dirEE.make<TH1D>("h_recHits_EE_eta","h_recHits_EE_eta",300,-300.,300.);
  h_recHits_EE_phi = dirEE.make<TH1D>("h_recHits_EE_phi","h_recHits_EE_phi",320,-3.2,3.2);
  
  h_recHits_EEP_size   = dirEE.make<TH1D>("h_recHits_EEP_size","h_recHits_EEP_size",1000,0,10000);
  h_recHits_EEP_energy = dirEE.make<TH1D>("h_recHits_EEP_energy","h_recHits_EEP_energy",2000,-50,350);
  h_recHits_EEP_time   = dirEE.make<TH1D>("h_recHits_EEP_time","h_recHits_EEP_time",400,-100,100);
  h_recHits_EEP_Chi2   = dirEE.make<TH1D>("h_recHits_EEP_Chi2","h_recHits_EEP_Chi2",1000,0,100);
  h_recHits_EEP_OutOfTimeChi2 = dirEE.make<TH1D>("h_recHits_EEP_OutOfTimeChi2","h_recHits_EEP_OutOfTimeChi2",1000,0,100);  
  h_recHits_EEP_occupancy = dirEE.make<TH2D>("h_recHits_EEP_occupancy","h_recHits_EEP_occupancy",100,0.,100.,100,0.,100. );

  h_recHits_EEM_size   = dirEE.make<TH1D>("h_recHits_EEM_size","h_recHits_EEM_size",1000,0,10000);
  h_recHits_EEM_energy = dirEE.make<TH1D>("h_recHits_EEM_energy","h_recHits_EEM_energy",2000,-50,350);
  h_recHits_EEM_time   = dirEE.make<TH1D>("h_recHits_EEM_time","h_recHits_EEM_time",400,-100,100);
  h_recHits_EEM_Chi2   = dirEE.make<TH1D>("h_recHits_EEM_Chi2","h_recHits_EEM_Chi2",1000,0,100);
  h_recHits_EEM_OutOfTimeChi2 = dirEE.make<TH1D>("h_recHits_EEM_OutOfTimeChi2","h_recHits_EEM_OutOfTimeChi2",1000,0,100);  
  h_recHits_EEM_occupancy = dirEE.make<TH2D>("h_recHits_EEM_occupancy","h_recHits_EEM_occupancy",100,0.,100.,100,0.,100. );
  
  // Basic Clusters ----------------------------------------------
  // ... barrel
  h_basicClusters_EB_size   = dirEB.make<TH1D>("h_basicClusters_EB_size","h_basicClusters_EB_size",200,0.,200.);
  h_basicClusters_EB_nXtals = dirEB.make<TH1D>("h_basicClusters_EB_nXtals","h_basicClusters_EB_nXtals",400,0.,400.);
  h_basicClusters_EB_energy = dirEB.make<TH1D>("h_basicClusters_EB_energy","h_basicClusters_EB_energy",2000,0.,400.);
  h_basicClusters_EB_eta    = dirEB.make<TH1D>("h_basicClusters_EB_eta","h_basicClusters_EB_eta",300,-3.,3.);
  h_basicClusters_EB_phi    = dirEB.make<TH1D>("h_basicClusters_EB_phi","h_basicClusters_EB_phi",320,-3.2,3.2);
  h_basicClusters_EB_seedFlag = dirEB.make<TH1D>("h_basicClusters_EB_seedFlag","h_basicClusters_EB_seedFlag",20,0,20);
  // ... endcap
  h_basicClusters_EEP_size   = dirEE.make<TH1D>("h_basicClusters_EEP_size","h_basicClusters_EEP_size",200,0.,200.);
  h_basicClusters_EEP_nXtals = dirEE.make<TH1D>("h_basicClusters_EEP_nXtals","h_basicClusters_EEP_nXtals",400,0.,400.);
  h_basicClusters_EEP_energy = dirEE.make<TH1D>("h_basicClusters_EEP_energy","h_basicClusters_EEP_energy",2000,0.,400.);
  h_basicClusters_EEP_seedFlag = dirEE.make<TH1D>("h_basicClusters_EEP_seedFlag","h_basicClusters_EEP_seedFlag",20,0,20);

  h_basicClusters_EEM_size   = dirEE.make<TH1D>("h_basicClusters_EEM_size","h_basicClusters_EEM_size",200,0.,200.);
  h_basicClusters_EEM_nXtals = dirEE.make<TH1D>("h_basicClusters_EEM_nXtals","h_basicClusters_EEM_nXtals",400,0.,400.);
  h_basicClusters_EEM_energy = dirEE.make<TH1D>("h_basicClusters_EEM_energy","h_basicClusters_EEM_energy",2000,0.,400.);
  h_basicClusters_EEM_seedFlag = dirEE.make<TH1D>("h_basicClusters_EEM_seedFlag","h_basicClusters_EEM_seedFlag",20,0,20);

  h_basicClusters_EE_eta     = dirEE.make<TH1D>("h_basicClusters_EE_eta","h_basicClusters_EE_eta",300,-3.,3.);
  h_basicClusters_EE_phi     = dirEE.make<TH1D>("h_basicClusters_EE_phi","h_basicClusters_EE_phi",320,-3.2,3.2);

  // Super Clusters ----------------------------------------------
  // ... barrel
  h_superClusters_EB_size   = dirEB.make<TH1D>("h_superClusters_EB_size","h_superClusters_EB_size",200,0.,200.);
  h_superClusters_EB_nXtals = dirEB.make<TH1D>("h_superClusters_EB_nXtals","h_superClusters_EB_nXtals",400,0.,400.);
  h_superClusters_EB_energy = dirEB.make<TH1D>("h_superClusters_EB_energy","h_superClusters_EB_energy",2000,0.,400.);
  h_superClusters_EB_eta    = dirEB.make<TH1D>("h_superClusters_EB_eta","h_superClusters_EB_eta",300,-3.,3.);
  h_superClusters_EB_phi    = dirEB.make<TH1D>("h_superClusters_EB_phi","h_superClusters_EB_phi",320,-3.2,3.2);
  h_superClusters_EB_E1oE9  = dirEB.make<TH1D>("h_superClusters_EB_E1oE9","h_superClusters_EB_E1oE9",150,0,1.5);

  // ... endcap
  h_superClusters_EEP_size   = dirEE.make<TH1D>("h_superClusters_EEP_size","h_superClusters_EEP_size",200,0.,200.);
  h_superClusters_EEP_nXtals = dirEE.make<TH1D>("h_superClusters_EEP_nXtals","h_superClusters_EEP_nXtals",400,0.,400.);
  h_superClusters_EEP_energy = dirEE.make<TH1D>("h_superClusters_EEP_energy","h_superClusters_EEP_energy",2000,0.,400.);
  h_superClusters_EEP_E1oE9  = dirEE.make<TH1D>("h_superClusters_EEP_E1oE9","h_superClusters_EEP_E1oE9",150,0,1.5);  

  h_superClusters_EEM_size   = dirEE.make<TH1D>("h_superClusters_EEM_size","h_superClusters_EEM_size",200,0.,200.);
  h_superClusters_EEM_nXtals = dirEE.make<TH1D>("h_superClusters_EEM_nXtals","h_superClusters_EEM_nXtals",400,0.,400.);
  h_superClusters_EEM_energy = dirEE.make<TH1D>("h_superClusters_EEM_energy","h_superClusters_EEM_energy",2000,0.,400.);
  h_superClusters_EEM_E1oE9  = dirEE.make<TH1D>("h_superClusters_EEM_E1oE9","h_superClusters_EEM_E1oE9",150,0,1.5);  

  h_superClusters_EE_eta     = dirEE.make<TH1D>("h_superClusters_EE_eta","h_superClusters_EE_eta",300,-3.,3.);
  h_superClusters_EE_phi     = dirEE.make<TH1D>("h_superClusters_EE_phi","h_superClusters_EE_phi",320,-3.2,3.2);

  // preshower
  h_esRecHits_energy_F[0] = dirES.make<TH1D>("h_esRecHits_energy_F+","ES+F rec hit energy",1000,0.,0.01);
  h_esRecHits_energy_F[1] = dirES.make<TH1D>("h_esRecHits_energy_F-","ES-F rec hit energy",1000,0.,0.01);
  h_esRecHits_energy_R[0] = dirES.make<TH1D>("h_esRecHits_energy_R+","ES+R rec hit energy",1000,0.,0.01);
  h_esRecHits_energy_R[1] = dirES.make<TH1D>("h_esRecHits_energy_R-","ES-R rec hit energy",1000,0.,0.01);
  h_esClusters_energy_plane1 = dirES.make<TH1D>("h_esClusters_energy_plane1","h_esClusters_energy_plane1",1000,0.,0.01);
  h_esClusters_energy_plane2 = dirES.make<TH1D>("h_esClusters_energy_plane2","h_esClusters_energy_plane2",1000,0.,0.01);
  h_esClusters_energy_ratio  = dirES.make<TH1D>("h_esClusters_energy_ratio","h_esClusters_energy_ratio",100,0.,10.);


}


EcalValidation::~EcalValidation()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void EcalValidation::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  // calo geometry
  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);
  const CaloGeometry *geometry = pGeometry.product();

  // calo topology
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  const CaloTopology *topology = pTopology.product();

  // RecHits
  // ... barrel
  edm::Handle<EcalRecHitCollection> recHitsEB;
  ev.getByLabel( recHitCollection_EB_, recHitsEB );
  const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product () ;
  if ( ! recHitsEB.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsEB not found" << std::endl; 
  }
  
  for ( EcalRecHitCollection::const_iterator itr = theBarrelEcalRecHits->begin () ;
	itr != theBarrelEcalRecHits->end () ;++itr)
    {
      EBDetId ebid( itr -> id() );
      
      h_recHits_EB_energy        -> Fill( itr -> energy() );
      h_recHits_EB_time          -> Fill( itr -> time() );
      h_recHits_EB_Chi2          -> Fill( itr -> chi2() );
      h_recHits_EB_OutOfTimeChi2 -> Fill( itr -> outOfTimeChi2() );
      h_recHits_EB_occupancy     -> Fill( ebid.ieta() , ebid.iphi() );

      GlobalPoint mycell = geometry -> getPosition(DetId(itr->id()));
      if (  itr -> energy() > ethrEB_ ) {
	h_recHits_EB_eta            -> Fill( mycell.eta() );
	h_recHits_EB_phi            -> Fill( mycell.phi() );
      }
    }
  h_recHits_EB_size->Fill( recHitsEB->size() );


  // ... endcap
  edm::Handle<EcalRecHitCollection> recHitsEE;
  ev.getByLabel( recHitCollection_EE_, recHitsEE );
  const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;
  if ( ! recHitsEE.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsEE not found" << std::endl; 
  }
  
  int nHitsEEP = 0;
  int nHitsEEM = 0;
  
  for ( EcalRecHitCollection::const_iterator itr = theEndcapEcalRecHits->begin () ;
	itr != theEndcapEcalRecHits->end () ; ++itr)
    {
      EEDetId eeid( itr -> id() );
      
      GlobalPoint mycell = geometry->getPosition(itr->detid());
      
      if (  itr -> energy() > ethrEE_ ) {
	h_recHits_EE_eta             -> Fill( mycell.eta() );
	h_recHits_EE_phi             -> Fill( mycell.phi() );
      }

      if ( eeid.zside() > 0 ){
	h_recHits_EEP_energy        -> Fill( itr -> energy() );
	h_recHits_EEP_time          -> Fill( itr -> time() );
	h_recHits_EEP_Chi2          -> Fill( itr -> chi2() );
	h_recHits_EEP_OutOfTimeChi2 -> Fill( itr -> outOfTimeChi2() );
	h_recHits_EEP_occupancy   -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5 );
	nHitsEEP++;
      }

      if ( eeid.zside() < 0 ){
	h_recHits_EEM_energy        -> Fill( itr -> energy() );
	h_recHits_EEM_time          -> Fill( itr -> time() );
	h_recHits_EEM_Chi2          -> Fill( itr -> chi2() );
	h_recHits_EEM_OutOfTimeChi2 -> Fill( itr -> outOfTimeChi2() );
	h_recHits_EEM_occupancy     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5 );
	nHitsEEM++;
      }
    }
  
  h_recHits_EEP_size->Fill( nHitsEEP );
  h_recHits_EEM_size->Fill( nHitsEEM );
  
  
  // Basic Clusters
  // ... barrel
  edm::Handle<reco::BasicClusterCollection> basicClusters_EB_h;
  ev.getByLabel( basicClusterCollection_EB_, basicClusters_EB_h );
  if ( ! basicClusters_EB_h.isValid() ) {
    std::cerr << "EcalValidation::analyze --> basicClusters_EB_h not found" << std::endl; 
  }
  h_basicClusters_EB_size->Fill( basicClusters_EB_h->size() );
  for (unsigned int icl = 0; icl < basicClusters_EB_h->size(); ++icl) {
    h_basicClusters_EB_energy -> Fill( (*basicClusters_EB_h)[icl].energy() );
    h_basicClusters_EB_eta    -> Fill( (*basicClusters_EB_h)[icl].eta() );
    h_basicClusters_EB_phi    -> Fill( (*basicClusters_EB_h)[icl].phi() );
    h_basicClusters_EB_nXtals -> Fill( (*basicClusters_EB_h)[icl].hitsAndFractions().size() );
    
    EcalRecHitCollection::const_iterator seed_itEB = recHitsEB->find((*basicClusters_EB_h)[icl].seed());
    if ( seed_itEB != recHitsEB->end() ) continue;
    h_basicClusters_EB_seedFlag->Fill( seed_itEB->recoFlag());
  }
  

  // ... endcap
  edm::Handle<reco::BasicClusterCollection> basicClusters_EE_h;
  ev.getByLabel( basicClusterCollection_EE_, basicClusters_EE_h );
  if ( ! basicClusters_EE_h.isValid() ) {
    std::cerr << "EcalValidation::analyze --> basicClusters_EE_h not found" << std::endl; 
  }
 

  int nBasicClustersEEP = 0;
  int nBasicClustersEEM = 0;

  for (unsigned int icl = 0; icl < basicClusters_EE_h->size(); ++icl) {

    h_basicClusters_EE_eta    -> Fill( (*basicClusters_EE_h)[icl].eta() );
    h_basicClusters_EE_phi    -> Fill( (*basicClusters_EE_h)[icl].phi() );

    if ((*basicClusters_EE_h)[icl].z() > 0){
      h_basicClusters_EEP_energy -> Fill( (*basicClusters_EE_h)[icl].energy() );
      h_basicClusters_EEP_nXtals -> Fill( (*basicClusters_EE_h)[icl].hitsAndFractions().size() );
      EcalRecHitCollection::const_iterator seed_itEE = recHitsEE->find((*basicClusters_EE_h)[icl].seed());
      h_basicClusters_EEP_seedFlag -> Fill( seed_itEE->recoFlag());
      nBasicClustersEEP++;
    }

    if ((*basicClusters_EE_h)[icl].z() < 0){
      h_basicClusters_EEM_energy -> Fill( (*basicClusters_EE_h)[icl].energy() );
      h_basicClusters_EEM_nXtals -> Fill( (*basicClusters_EE_h)[icl].hitsAndFractions().size() );
      EcalRecHitCollection::const_iterator seed_itEE = recHitsEE->find((*basicClusters_EE_h)[icl].seed());
      h_basicClusters_EEM_seedFlag -> Fill( seed_itEE->recoFlag());
      nBasicClustersEEM++;
    }
  }
  h_basicClusters_EEP_size->Fill( nBasicClustersEEP );
  h_basicClusters_EEM_size->Fill( nBasicClustersEEM );


  // Super Clusters
  // ... barrel
  edm::Handle<reco::SuperClusterCollection> superClusters_EB_h;
  ev.getByLabel( superClusterCollection_EB_, superClusters_EB_h );
  const reco::SuperClusterCollection* theBarrelSuperClusters = superClusters_EB_h.product () ;
  if ( ! superClusters_EB_h.isValid() ) {
    std::cerr << "EcalValidation::analyze --> superClusters_EB_h not found" << std::endl; 
  }
 
  for (reco::SuperClusterCollection::const_iterator itSC = theBarrelSuperClusters->begin(); 
       itSC != theBarrelSuperClusters->end(); ++itSC ) {
    h_superClusters_EB_energy -> Fill( itSC -> energy() );
    h_superClusters_EB_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
    h_superClusters_EB_eta    -> Fill( itSC -> eta() );
    h_superClusters_EB_phi    -> Fill( itSC -> phi() );
    float E1oE9 = EcalClusterTools::eMax( *itSC, theBarrelEcalRecHits)/
                  EcalClusterTools::e3x3( *itSC, theBarrelEcalRecHits, topology );
    h_superClusters_EB_E1oE9  -> Fill( E1oE9 );
  }
  
  h_superClusters_EB_size->Fill( superClusters_EB_h->size() );

  
  // ... endcap
  edm::Handle<reco::SuperClusterCollection> superClusters_EE_h;
  ev.getByLabel( superClusterCollection_EE_, superClusters_EE_h );
  const reco::SuperClusterCollection* theEndcapSuperClusters = superClusters_EE_h.product () ;
  if ( ! superClusters_EE_h.isValid() ) {
    std::cerr << "EcalValidation::analyze --> superClusters_EE_h not found" << std::endl; 
  }

  int nSuperClustersEEP = 0;
  int nSuperClustersEEM = 0;

  for (reco::SuperClusterCollection::const_iterator itSC = theEndcapSuperClusters->begin(); 
       itSC != theEndcapSuperClusters->end(); ++itSC ) {
    h_superClusters_EE_eta    -> Fill( itSC -> eta() );
    h_superClusters_EE_phi    -> Fill( itSC -> phi() );
    
    float E1oE9 = EcalClusterTools::eMax( *itSC, theEndcapEcalRecHits)/
                  EcalClusterTools::e3x3( *itSC, theEndcapEcalRecHits, topology );

    if  ( itSC -> z() > 0 ){
      h_superClusters_EEP_energy -> Fill( itSC -> energy() );
      h_superClusters_EEP_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
      h_superClusters_EEP_E1oE9  -> Fill( E1oE9 );
      nSuperClustersEEP++;
    }

    if  ( itSC -> z() < 0 ){
      h_superClusters_EEM_energy -> Fill( itSC -> energy() );
      h_superClusters_EEM_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
      h_superClusters_EEM_E1oE9  -> Fill( E1oE9 );
      nSuperClustersEEM++;
    }
  }
  
  h_superClusters_EEP_size->Fill( nSuperClustersEEP );
  h_superClusters_EEM_size->Fill( nSuperClustersEEM );

  //--------------------------------------------------------
 
  //Preshower RecHits
  edm::Handle<ESRecHitCollection> recHitsES;
  ev.getByLabel (esRecHitCollection_, recHitsES) ;
  const ESRecHitCollection* thePreShowerRecHits = recHitsES.product () ;

  if ( ! recHitsES.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsES not found" << std::endl; 
  }

  for (ESRecHitCollection::const_iterator esItr = thePreShowerRecHits->begin();
       esItr != thePreShowerRecHits->end(); ++esItr)
    {
      ESDetId id = ESDetId(esItr->id());
      
      // front plane : id.plane()==1
      if ( id.plane()==1 && id.zside() > 0 )
	h_esRecHits_energy_F[0]->Fill( esItr->energy() );

      if ( id.plane()==1 && id.zside() < 0 )
	h_esRecHits_energy_F[1]->Fill( esItr->energy() );
      
      // rear plane : id.plane()==2
      if ( id.plane()==2 && id.zside() > 0 )
	h_esRecHits_energy_R[0]->Fill( esItr->energy() );

      if ( id.plane()==2 && id.zside() < 0 )
	h_esRecHits_energy_R[1]->Fill( esItr->energy() );
    
    }

  
  // ES clusters in X plane
  Handle<PreshowerClusterCollection> esClustersX;
  ev.getByLabel( esClusterCollectionX_, esClustersX);
  const PreshowerClusterCollection *ESclustersX = esClustersX.product();

  // ES clusters in Y plane
  Handle<PreshowerClusterCollection> esClustersY;
  ev.getByLabel( esClusterCollectionY_, esClustersY);
  const PreshowerClusterCollection *ESclustersY = esClustersY.product(); 
  
  // loop over all super clusters
  for (reco::SuperClusterCollection::const_iterator itSC = theEndcapSuperClusters->begin(); 
       itSC != theEndcapSuperClusters->end(); ++itSC ) {
    
    if ( fabs(itSC->eta()) < 1.65 || fabs(itSC->eta()) > 2.6 ) continue;

    // Loop over all ECAL Basic clusters in the supercluster
    for (CaloCluster_iterator ecalBasicCluster = itSC->clustersBegin(); ecalBasicCluster!= itSC->clustersEnd(); 
	 ecalBasicCluster++) {
      const CaloClusterPtr ecalBasicClusterPtr = *(ecalBasicCluster);
      
      float ESenergyPlane1 = -999.;
      float ESenergyPlane2 = -999.;
      
      for (PreshowerClusterCollection::const_iterator iESClus = ESclustersX->begin(); iESClus != ESclustersX->end(); 
	   ++iESClus) {
        const CaloClusterPtr preshBasicCluster = iESClus->basicCluster();
        const PreshowerCluster *esCluster = &*iESClus;
        if (preshBasicCluster == ecalBasicClusterPtr) {
	  ESenergyPlane1 = esCluster->energy();
	  h_esClusters_energy_plane1 ->Fill(esCluster->energy());
	}
      }  // end of x loop
      
      for (PreshowerClusterCollection::const_iterator iESClus = ESclustersY->begin(); iESClus != ESclustersY->end(); 
	   ++iESClus) {
        const CaloClusterPtr preshBasicCluster = iESClus->basicCluster();
        const PreshowerCluster *esCluster = &*iESClus;
        if (preshBasicCluster == ecalBasicClusterPtr) {
	  ESenergyPlane2 = esCluster->energy();
	  h_esClusters_energy_plane2 -> Fill(esCluster->energy());
	}
      } // end of y loop
      
      if ( ESenergyPlane1 != -999. && ESenergyPlane2 != -999. ) 
	h_esClusters_energy_ratio -> Fill(ESenergyPlane1/ESenergyPlane2);
      
      
    } // end loop over all basic clusters in the supercluster
  }// end loop over superclusters

}


// ------------ method called once each job just before starting event loop  ------------
        void 
EcalValidation::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EcalValidation::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EcalValidation);
