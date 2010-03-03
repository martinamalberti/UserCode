// -*- C++ -*-
//
// Package:    EcalValidation
// Class:      EcalValidation
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

  naiveId_ = 0;
 
  // histos 
  
  edm::Service<TFileService> fs;
  
//   TFileDirectory dirEB = fs->mkdir("EB");
//   TFileDirectory dirEE = fs->mkdir("EE");
//   TFileDirectory dirES = fs->mkdir("ES");
  
  h_numberOfEvents = fs->make<TH1D>("h_numberOfEvents","h_numberOfEvents",10,0,10);
  
  // RecHits ---------------------------------------------- 
  // ... barrel
  h_recHits_EB_size          = fs->make<TH1D>("h_recHits_EB_size", "h_recHitsEB_size", 1000, 0, 10000 );
  h_recHits_EB_energy        = fs->make<TH1D>("h_recHits_EB_energy","h_recHitsEB_energy",2000,-50,500);
  h_recHits_EB_energyMax     = fs->make<TH1D>("h_recHits_EB_energyMax","h_recHitsEB_energyMax",2000,-50,500);
  h_recHits_EB_time          = fs->make<TH1D>("h_recHits_EB_time","h_recHits_EB_time",400,-100,100);
  h_recHits_EB_Chi2          = fs->make<TH1D>("h_recHits_EB_Chi2","h_recHits_EB_Chi2",1000,0,100);
  h_recHits_EB_OutOfTimeChi2 = fs->make<TH1D>("h_recHits_EB_OutOfTimeChi2","h_recHits_EB_OutOfTimeChi2",1000,0,100);  
  h_recHits_EB_occupancy     = fs->make<TH2D>("h_recHits_EB_occupancy","h_recHits_EB_occupancy",172,-86.,86.,360,1.,361. );

  // ... endcap
  h_recHits_EEP_size          = fs->make<TH1D>("h_recHits_EEP_size","h_recHits_EEP_size",1000,0,10000);
  h_recHits_EEP_energy        = fs->make<TH1D>("h_recHits_EEP_energy","h_recHits_EEP_energy",2000,-50,500);
  h_recHits_EEP_energyMax     = fs->make<TH1D>("h_recHits_EEP_energyMax","h_recHitsEEP_energyMax",2000,-50,500);
  h_recHits_EEP_time          = fs->make<TH1D>("h_recHits_EEP_time","h_recHits_EEP_time",400,-100,100);
  h_recHits_EEP_Chi2          = fs->make<TH1D>("h_recHits_EEP_Chi2","h_recHits_EEP_Chi2",1000,0,100);
  h_recHits_EEP_OutOfTimeChi2 = fs->make<TH1D>("h_recHits_EEP_OutOfTimeChi2","h_recHits_EEP_OutOfTimeChi2",1000,0,100);  
  h_recHits_EEP_occupancy     = fs->make<TH2D>("h_recHits_EEP_occupancy","h_recHits_EEP_occupancy",100,0.,100.,100,0.,100. );

  h_recHits_EEM_size          = fs->make<TH1D>("h_recHits_EEM_size","h_recHits_EEM_size",1000,0,10000);
  h_recHits_EEM_energy        = fs->make<TH1D>("h_recHits_EEM_energy","h_recHits_EEM_energy",2000,-50,500);
  h_recHits_EEM_energyMax     = fs->make<TH1D>("h_recHits_EEM_energyMax","h_recHits_EEM_energyMax",2000,-50,500);
  h_recHits_EEM_time          = fs->make<TH1D>("h_recHits_EEM_time","h_recHits_EEM_time",400,-100,100);
  h_recHits_EEM_Chi2          = fs->make<TH1D>("h_recHits_EEM_Chi2","h_recHits_EEM_Chi2",1000,0,100);
  h_recHits_EEM_OutOfTimeChi2 = fs->make<TH1D>("h_recHits_EEM_OutOfTimeChi2","h_recHits_EEM_OutOfTimeChi2",1000,0,100);  
  h_recHits_EEM_occupancy     = fs->make<TH2D>("h_recHits_EEM_occupancy","h_recHits_EEM_occupancy",100,0.,100.,100,0.,100. );
  
  h_recHits_eta               = fs->make<TH1D>("h_recHits_eta","h_recHits_eta",300,-3.,3.);
  h_recHits_EB_phi            = fs->make<TH1D>("h_recHits_EB_phi","h_recHits_EB_phi",320,-3.2,3.2);
  h_recHits_EE_phi            = fs->make<TH1D>("h_recHits_EE_phi","h_recHits_EE_phi",320,-3.2,3.2);


  // Basic Clusters ----------------------------------------------
  // ... barrel
  h_basicClusters_EB_size    = fs->make<TH1D>("h_basicClusters_EB_size","h_basicClusters_EB_size",200,0.,200.);
  h_basicClusters_EB_nXtals  = fs->make<TH1D>("h_basicClusters_EB_nXtals","h_basicClusters_EB_nXtals",400,0.,400.);
  h_basicClusters_EB_energy  = fs->make<TH1D>("h_basicClusters_EB_energy","h_basicClusters_EB_energy",2000,0.,400.);

  // ... endcap
  h_basicClusters_EEP_size   = fs->make<TH1D>("h_basicClusters_EEP_size","h_basicClusters_EEP_size",200,0.,200.);
  h_basicClusters_EEP_nXtals = fs->make<TH1D>("h_basicClusters_EEP_nXtals","h_basicClusters_EEP_nXtals",400,0.,400.);
  h_basicClusters_EEP_energy = fs->make<TH1D>("h_basicClusters_EEP_energy","h_basicClusters_EEP_energy",2000,0.,400.);

  h_basicClusters_EEM_size   = fs->make<TH1D>("h_basicClusters_EEM_size","h_basicClusters_EEM_size",200,0.,200.);
  h_basicClusters_EEM_nXtals = fs->make<TH1D>("h_basicClusters_EEM_nXtals","h_basicClusters_EEM_nXtals",400,0.,400.);
  h_basicClusters_EEM_energy = fs->make<TH1D>("h_basicClusters_EEM_energy","h_basicClusters_EEM_energy",2000,0.,400.);
  
  h_basicClusters_eta        = fs->make<TH1D>("h_basicClusters_eta","h_basicClusters_eta",300,-3.,3.);
  h_basicClusters_EB_phi     = fs->make<TH1D>("h_basicClusters_EB_phi","h_basicClusters_EB_phi",320,-3.2,3.2);
  h_basicClusters_EE_phi     = fs->make<TH1D>("h_basicClusters_EE_phi","h_basicClusters_EE_phi",320,-3.2,3.2);

  // Super Clusters ----------------------------------------------
  // ... barrel
  h_superClusters_EB_size    = fs->make<TH1D>("h_superClusters_EB_size","h_superClusters_EB_size",200,0.,200.);
  h_superClusters_EB_nXtals  = fs->make<TH1D>("h_superClusters_EB_nXtals","h_superClusters_EB_nXtals",400,0.,400.);
  h_superClusters_EB_nBC     = fs->make<TH1D>("h_superClusters_EB_nBC","h_superClusters_EB_nBC",100,0.,100.);
  h_superClusters_EB_energy  = fs->make<TH1D>("h_superClusters_EB_energy","h_superClusters_EB_energy",2000,0.,400.);
  h_superClusters_EB_E1oE9   = fs->make<TH1D>("h_superClusters_EB_E1oE9","h_superClusters_EB_E1oE9",150,0,1.5);

  // ... endcap
  h_superClusters_EEP_size   = fs->make<TH1D>("h_superClusters_EEP_size","h_superClusters_EEP_size",200,0.,200.);
  h_superClusters_EEP_nXtals = fs->make<TH1D>("h_superClusters_EEP_nXtals","h_superClusters_EEP_nXtals",400,0.,400.);
  h_superClusters_EEP_nBC    = fs->make<TH1D>("h_superClusters_EEP_nBC","h_superClusters_EEP_nBC",100,0.,100.);
  h_superClusters_EEP_energy = fs->make<TH1D>("h_superClusters_EEP_energy","h_superClusters_EEP_energy",2000,0.,400.);
  h_superClusters_EEP_E1oE9  = fs->make<TH1D>("h_superClusters_EEP_E1oE9","h_superClusters_EEP_E1oE9",150,0,1.5);  

  h_superClusters_EEM_size   = fs->make<TH1D>("h_superClusters_EEM_size","h_superClusters_EEM_size",200,0.,200.);
  h_superClusters_EEM_nXtals = fs->make<TH1D>("h_superClusters_EEM_nXtals","h_superClusters_EEM_nXtals",400,0.,400.);
  h_superClusters_EEM_nBC    = fs->make<TH1D>("h_superClusters_EEM_nBC","h_superClusters_EEM_nBC",100,0.,100.);
  h_superClusters_EEM_energy = fs->make<TH1D>("h_superClusters_EEM_energy","h_superClusters_EEM_energy",2000,0.,400.);
  h_superClusters_EEM_E1oE9  = fs->make<TH1D>("h_superClusters_EEM_E1oE9","h_superClusters_EEM_E1oE9",150,0,1.5);  

  h_superClusters_eta        = fs->make<TH1D>("h_superClusters_eta","h_superClusters_eta",300,-3.,3.);
  h_superClusters_EB_phi     = fs->make<TH1D>("h_superClusters_EB_phi","h_superClusters_EB_phi",320,-3.2,3.2);
  h_superClusters_EE_phi     = fs->make<TH1D>("h_superClusters_EE_phi","h_superClusters_EE_phi",320,-3.2,3.2);

  // preshower
  h_recHits_ES_size          = fs->make<TH1D>("h_recHits_ES_size","h_recHits_ES_size",1000,0.,10000);
  h_recHits_ES_energy        = fs->make<TH1D>("h_recHits_ES_energy","h_recHits_ES_energy",1000,0.,0.01);
  h_recHits_ES_energyMax     = fs->make<TH1D>("h_recHits_ES_energyMax","h_recHits_ES_energyMax",1000,0.,0.01);
  h_recHits_ES_time          = fs->make<TH1D>("h_recHits_ES_time","h_recHits_ES_time",400,-100.,100.);
  h_recHits_ES_energy_F[0]   = fs->make<TH1D>("h_recHits_ES_energy_F+","h_recHits_ES_energy_F+",1000,0.,0.01);
  h_recHits_ES_energy_F[1]   = fs->make<TH1D>("h_recHits_ES_energy_F-","h_recHits_ES_energy_F-",1000,0.,0.01);
  h_recHits_ES_energy_R[0]   = fs->make<TH1D>("h_recHits_ES_energy_R+","h_recHits_ES_energy_R+",1000,0.,0.01);
  h_recHits_ES_energy_R[1]   = fs->make<TH1D>("h_recHits_ES_energy_R-","h_recHits_ES_energy_R-",1000,0.,0.01);
  h_esClusters_energy_plane1 = fs->make<TH1D>("h_esClusters_energy_plane1","h_esClusters_energy_plane1",1000,0.,0.01);
  h_esClusters_energy_plane2 = fs->make<TH1D>("h_esClusters_energy_plane2","h_esClusters_energy_plane2",1000,0.,0.01);
  h_esClusters_energy_ratio  = fs->make<TH1D>("h_esClusters_energy_ratio","h_esClusters_energy_ratio",100,0.,10.);


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

  naiveId_++;

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
  
  float maxRecHitEnergyEB = -999.;
  
  for ( EcalRecHitCollection::const_iterator itr = theBarrelEcalRecHits->begin () ;
	itr != theBarrelEcalRecHits->end () ;++itr)
    {
      
      EBDetId ebid( itr -> id() );
      
      if (itr -> energy() > maxRecHitEnergyEB ) maxRecHitEnergyEB = itr -> energy() ;
       
      h_recHits_EB_energy        -> Fill( itr -> energy() );
      if (  itr -> energy() > ethrEB_ )
	h_recHits_EB_time        -> Fill( itr -> time() );
      h_recHits_EB_Chi2          -> Fill( itr -> chi2() );
      h_recHits_EB_OutOfTimeChi2 -> Fill( itr -> outOfTimeChi2() );
      h_recHits_EB_occupancy     -> Fill( ebid.ieta() , ebid.iphi() );

      GlobalPoint mycell = geometry -> getPosition(DetId(itr->id()));
      //if (  itr -> energy() > ethrEB_ ) {
	h_recHits_eta            -> Fill( mycell.eta() );
	h_recHits_EB_phi         -> Fill( mycell.phi() );
	//}

      
    }
  h_recHits_EB_energyMax -> Fill( maxRecHitEnergyEB  );
  h_recHits_EB_size      -> Fill( recHitsEB->size() );


  // ... endcap
  edm::Handle<EcalRecHitCollection> recHitsEE;
  ev.getByLabel( recHitCollection_EE_, recHitsEE );
  const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;
  if ( ! recHitsEE.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsEE not found" << std::endl; 
  }
  
  int nHitsEEP = 0;
  int nHitsEEM = 0;

  float maxRecHitEnergyEEP = -999.;
  float maxRecHitEnergyEEM = -999.;

  for ( EcalRecHitCollection::const_iterator itr = theEndcapEcalRecHits->begin () ;
	itr != theEndcapEcalRecHits->end () ; ++itr)
    {
      EEDetId eeid( itr -> id() );
      
      GlobalPoint mycell = geometry->getPosition(itr->detid());
      //if (  itr -> energy() > ethrEE_ ) {
	h_recHits_eta        -> Fill( mycell.eta() );
	h_recHits_EE_phi     -> Fill( mycell.phi() );
	//}

      if ( eeid.zside() > 0 ){
	if (itr -> energy() > maxRecHitEnergyEEP ) maxRecHitEnergyEEP = itr -> energy() ;
	h_recHits_EEP_energy        -> Fill( itr -> energy() );
	if (  itr -> energy() > ethrEE_ ) 
	  h_recHits_EEP_time          -> Fill( itr -> time() );
	h_recHits_EEP_Chi2          -> Fill( itr -> chi2() );
	h_recHits_EEP_OutOfTimeChi2 -> Fill( itr -> outOfTimeChi2() );
	h_recHits_EEP_occupancy     -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5 );
	nHitsEEP++;
      }

      if ( eeid.zside() < 0 ){
	if (itr -> energy() > maxRecHitEnergyEEM ) maxRecHitEnergyEEM = itr -> energy() ;
	h_recHits_EEM_energy        -> Fill( itr -> energy() );
	if (  itr -> energy() > ethrEE_ ) 
	  h_recHits_EEM_time          -> Fill( itr -> time() );
	h_recHits_EEM_Chi2          -> Fill( itr -> chi2() );
	h_recHits_EEM_OutOfTimeChi2 -> Fill( itr -> outOfTimeChi2() );
	h_recHits_EEM_occupancy     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5 );
	nHitsEEM++;
      }
    }

  h_recHits_EEP_energyMax -> Fill( maxRecHitEnergyEEP );
  h_recHits_EEM_energyMax -> Fill( maxRecHitEnergyEEM );
  h_recHits_EEP_size      -> Fill( nHitsEEP );
  h_recHits_EEM_size      -> Fill( nHitsEEM );
  
  
  // Basic Clusters
  // ... barrel
  edm::Handle<reco::BasicClusterCollection> basicClusters_EB_h;
  ev.getByLabel( basicClusterCollection_EB_, basicClusters_EB_h );
  if ( ! basicClusters_EB_h.isValid() ) {
    std::cerr << "EcalValidation::analyze --> basicClusters_EB_h not found" << std::endl; 
  }
  
  h_basicClusters_EB_size     -> Fill( basicClusters_EB_h->size() );
  
  for (unsigned int icl = 0; icl < basicClusters_EB_h->size(); ++icl) {
    h_basicClusters_EB_nXtals -> Fill( (*basicClusters_EB_h)[icl].hitsAndFractions().size() );
    h_basicClusters_EB_energy -> Fill( (*basicClusters_EB_h)[icl].energy() );
    h_basicClusters_eta       -> Fill( (*basicClusters_EB_h)[icl].eta() );
    h_basicClusters_EB_phi    -> Fill( (*basicClusters_EB_h)[icl].phi() );
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

    h_basicClusters_eta       -> Fill( (*basicClusters_EE_h)[icl].eta() );
    h_basicClusters_EE_phi    -> Fill( (*basicClusters_EE_h)[icl].phi() );

    if ((*basicClusters_EE_h)[icl].z() > 0){
      h_basicClusters_EEP_nXtals -> Fill( (*basicClusters_EE_h)[icl].hitsAndFractions().size() );
      h_basicClusters_EEP_energy -> Fill( (*basicClusters_EE_h)[icl].energy() );
      nBasicClustersEEP++;
    }

    if ((*basicClusters_EE_h)[icl].z() < 0){
      h_basicClusters_EEM_nXtals -> Fill( (*basicClusters_EE_h)[icl].hitsAndFractions().size() );
      h_basicClusters_EEM_energy -> Fill( (*basicClusters_EE_h)[icl].energy() );
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

  h_superClusters_EB_size->Fill( superClusters_EB_h->size() );

  for (reco::SuperClusterCollection::const_iterator itSC = theBarrelSuperClusters->begin(); 
       itSC != theBarrelSuperClusters->end(); ++itSC ) {
    h_superClusters_EB_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
    h_superClusters_EB_nBC    -> Fill( itSC -> clustersSize());
    h_superClusters_EB_energy -> Fill( itSC -> energy() );
    h_superClusters_eta       -> Fill( itSC -> eta() );
    h_superClusters_EB_phi    -> Fill( itSC -> phi() );
    float E1oE9 = EcalClusterTools::eMax( *itSC, theBarrelEcalRecHits)/
                  EcalClusterTools::e3x3( *itSC, theBarrelEcalRecHits, topology );
    h_superClusters_EB_E1oE9  -> Fill( E1oE9 );
  }
 
  
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
    h_superClusters_eta       -> Fill( itSC -> eta() );
    h_superClusters_EE_phi    -> Fill( itSC -> phi() );
    
    float E1oE9 = EcalClusterTools::eMax( *itSC, theEndcapEcalRecHits)/
                  EcalClusterTools::e3x3( *itSC, theEndcapEcalRecHits, topology );

    if  ( itSC -> z() > 0 ){
      h_superClusters_EEP_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
      h_superClusters_EEP_nBC    -> Fill( itSC -> clustersSize() );      
      h_superClusters_EEP_energy -> Fill( itSC -> energy() );
      h_superClusters_EEP_E1oE9  -> Fill( E1oE9 );
      nSuperClustersEEP++;
    }

    if  ( itSC -> z() < 0 ){
      h_superClusters_EEM_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
      h_superClusters_EEM_nBC    -> Fill( itSC -> clustersSize() );      
      h_superClusters_EEM_energy -> Fill( itSC -> energy() );
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

  h_recHits_ES_size -> Fill( recHitsES->size());

  float maxRecHitEnergyES = -999.;

  for (ESRecHitCollection::const_iterator esItr = thePreShowerRecHits->begin(); esItr != thePreShowerRecHits->end(); ++esItr) 
    {
      
      h_recHits_ES_energy -> Fill(esItr->energy()); 
      h_recHits_ES_time   -> Fill(esItr->time()); 
      if (esItr -> energy() > maxRecHitEnergyES ) maxRecHitEnergyES = esItr -> energy() ;

      ESDetId id = ESDetId(esItr->id());
      // front plane : id.plane()==1
      if ( id.plane()==1 && id.zside() > 0 )
	h_recHits_ES_energy_F[0]->Fill( esItr->energy() );

      if ( id.plane()==1 && id.zside() < 0 )
	h_recHits_ES_energy_F[1]->Fill( esItr->energy() );
      
      // rear plane : id.plane()==2
      if ( id.plane()==2 && id.zside() > 0 )
	h_recHits_ES_energy_R[0]->Fill( esItr->energy() );

      if ( id.plane()==2 && id.zside() < 0 )
	h_recHits_ES_energy_R[1]->Fill( esItr->energy() );
      
    } // end loop over ES rec Hits

  h_recHits_ES_energyMax -> Fill(maxRecHitEnergyES ); 
  
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
  
  h_numberOfEvents ->Fill(0.,naiveId_);

}

//define this as a plug-in
DEFINE_FWK_MODULE(EcalValidation);
