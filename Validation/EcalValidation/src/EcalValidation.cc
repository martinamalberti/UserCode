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
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"

#include "Validation/EcalValidation/interface/EcalValidation.h"


#include "TVector3.h"

#include <iostream>
#include <cmath>
#include <fstream>

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

  
 

  ///for Pi0 barrel selection
  isMonEBpi0_ = ps.getUntrackedParameter<bool>("isMonEBpi0",false);
  
  selePtGamma_ = ps.getParameter<double> ("selePtGamma");  
  selePtPi0_ = ps.getParameter<double> ("selePtPi0");  
  seleMinvMaxPi0_ = ps.getParameter<double> ("seleMinvMaxPi0");  
  seleMinvMinPi0_ = ps.getParameter<double> ("seleMinvMinPi0");  
  seleS4S9Gamma_ = ps.getParameter<double> ("seleS4S9Gamma");  
  selePi0Iso_ = ps.getParameter<double> ("selePi0Iso");  
  ptMinForIsolation_ = ps.getParameter<double> ("ptMinForIsolation");
  selePi0BeltDR_ = ps.getParameter<double> ("selePi0BeltDR");  
  selePi0BeltDeta_ = ps.getParameter<double> ("selePi0BeltDeta");  


  ///for Pi0 endcap selection
  isMonEEpi0_ = ps.getUntrackedParameter<bool>("isMonEEpi0",false);
  
  selePtGamma_EE_ = ps.getParameter<double> ("selePtGamma_EE");  
  selePtPi0_EE_ = ps.getParameter<double> ("selePtPi0_EE");   
  seleS4S9Gamma_EE_ = ps.getParameter<double> ("seleS4S9Gamma_EE");  
  seleMinvMaxPi0_EE_ = ps.getParameter<double> ("seleMinvMaxPi0_EE");  
  seleMinvMinPi0_EE_ = ps.getParameter<double> ("seleMinvMinPi0_EE");  
  ptMinForIsolation_EE_ = ps.getParameter<double> ("ptMinForIsolation_EE");
  selePi0BeltDR_EE_ = ps.getParameter<double> ("selePi0BeltDR_EE");  
  selePi0BeltDeta_EE_ = ps.getParameter<double> ("selePi0BeltDeta_EE");  
  selePi0Iso_EE_ = ps.getParameter<double> ("selePi0Iso_EE");  
  
  ///Endcap regions for pi0:
  ///try to divide endcap region into 3 parts
  /// eta< 2 ; eta>2 && eta<2.5 ; eta>2.5; 
  region1_Pi0_EE_ = ps.getParameter<double> ("region1_Pi0_EE");
  selePtGammaPi0_EE_region1_ = ps.getParameter<double> ("selePtGammaPi0_EE_region1");  
  selePtPi0_EE_region1_ = ps.getParameter<double> ("selePtPi0_EE_region1");   

  region2_Pi0_EE_ = ps.getParameter<double> ("region2_Pi0_EE");
  selePtGammaPi0_EE_region2_ = ps.getParameter<double> ("selePtGammaPi0_EE_region2");  
  selePtPi0_EE_region2_ = ps.getParameter<double> ("selePtPi0_EE_region2");   

  selePtGammaPi0_EE_region3_ = ps.getParameter<double> ("selePtGammaPi0_EE_region3");  
  selePtPi0_EE_region3_ = ps.getParameter<double> ("selePtPi0_EE_region3"); 

  isMaskEB_ = ps.getUntrackedParameter<bool>("isMaskEB",false);
  isMaskEE_ = ps.getUntrackedParameter<bool>("isMaskEE",false);

  maskEBFile_=  ps.getUntrackedParameter<string>("maskEBFile","maskEB.txt");
  maskEEFile_=  ps.getUntrackedParameter<string>("maskEEFile","maskEE.txt");

  useRecoFlag_ = ps.getUntrackedParameter<bool>("useRecoFlag",false);

  clusSeedThr_ = ps.getParameter<double> ("clusSeedThr");
  clusSeedThr_EE_ = ps.getParameter<double> ("clusSeedThr_EE");
  clusEtaSize_ = ps.getParameter<int> ("clusEtaSize");
  clusPhiSize_ = ps.getParameter<int> ("clusPhiSize");
  
  seleXtalMinEnergy_ = ps.getParameter<double>("seleXtalMinEnergy");
  seleXtalMinEnergy_EE_ = ps.getParameter<double>("seleXtalMinEnergy_EE");
  
  ParameterLogWeighted_ = ps.getParameter<bool> ("ParameterLogWeighted");
  ParameterX0_ = ps.getParameter<double> ("ParameterX0");
  ParameterT0_barl_ = ps.getParameter<double> ("ParameterT0_barl");
  ParameterT0_endc_ = ps.getParameter<double> ("ParameterT0_endc");
  ParameterT0_endcPresh_ = ps.getParameter<double> ("ParameterT0_endcPresh");
  ParameterW0_ = ps.getParameter<double> ("ParameterW0");

  std::map<std::string,double> providedParameters;
  providedParameters.insert(std::make_pair("LogWeighted",ParameterLogWeighted_));
  providedParameters.insert(std::make_pair("X0",ParameterX0_));
  providedParameters.insert(std::make_pair("T0_barl",ParameterT0_barl_));
  providedParameters.insert(std::make_pair("T0_endc",ParameterT0_endc_));
  providedParameters.insert(std::make_pair("T0_endcPresh",ParameterT0_endcPresh_));
  providedParameters.insert(std::make_pair("W0",ParameterW0_));

  posCalculator_ = PositionCalc(providedParameters);



  naiveId_ = 0;
  
  // histos 
  
  edm::Service<TFileService> fs;
  
  //   TFileDirectory dirEB = fs->mkdir("EB");
  //   TFileDirectory dirEE = fs->mkdir("EE");
  //   TFileDirectory dirES = fs->mkdir("ES");

  TFileDirectory dirPi0 = fs->mkdir("Pi0");

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

  // Pi0 ----------------------------------------------
  // ... barrel
  h_Pi0_EB_mass   = dirPi0.make<TH1D>("h_Pi0_EB_mass","h_Pi0_EB_mass",100,0.,0.5);
  h_Pi0_EB_pt1    = dirPi0.make<TH1D>("h_Pi0_EB_pt1","h_Pi0_EB_pt1",100,0.,20.);
  h_Pi0_EB_pt2    = dirPi0.make<TH1D>("h_Pi0_EB_pt2","h_Pi0_EB_pt2",100,0.,20.);
  h_Pi0_EB_pt     = dirPi0.make<TH1D>("h_Pi0_EB_pt","h_Pi0_EB_pt",100,0.,20.);
  h_Pi0_EB_eta    = dirPi0.make<TH1D>("h_Pi0_EB_eta","h_Pi0_EB_eta",100,-1.5,1.5);
  h_Pi0_EB_phi    = dirPi0.make<TH1D>("h_Pi0_EB_phi","h_Pi0_EB_phi",100,-3.2,3.2);

  // ... endcap
  h_Pi0_EE_mass   = dirPi0.make<TH1D>("h_Pi0_EE_mass","h_Pi0_EE_mass",100,0.,0.5);
  h_Pi0_EE_pt1    = dirPi0.make<TH1D>("h_Pi0_EE_pt1","h_Pi0_EE_pt1",100,0.,20.);
  h_Pi0_EE_pt2    = dirPi0.make<TH1D>("h_Pi0_EE_pt2","h_Pi0_EE_pt2",100,0.,20.);
  h_Pi0_EE_pt     = dirPi0.make<TH1D>("h_Pi0_EE_pt","h_Pi0_EE_pt",100,0.,20.);
  h_Pi0_EE_eta    = dirPi0.make<TH1D>("h_Pi0_EE_eta","h_Pi0_EE_eta",100,-3.,3.);
  h_Pi0_EE_phi    = dirPi0.make<TH1D>("h_Pi0_EE_phi","h_Pi0_EE_phi",100,-3.2,3.2);

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


  // ---------- Do histos for pi0 peak

  doPi0Barrel(geometry, topology, recHitsEB);
  doPi0Endcap(geometry, topology, recHitsEE);

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


// ------------ Pi0 peak  ECAL BARREL ------------
void EcalValidation::doPi0Barrel ( const CaloGeometry *geometry,
				   const CaloTopology *topology,
				   edm::Handle<EcalRecHitCollection> recHitsEB ) 
{
  const CaloSubdetectorGeometry *geometry_EB = geometry->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
  const CaloSubdetectorGeometry *geometry_ES = geometry->getSubdetectorGeometry(DetId::Ecal,EcalPreshower);

  const CaloSubdetectorTopology *topology_EB = topology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);
  
  const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product () ;

  //---------------------Barrel--------------------------------//
  if (isMonEBpi0_ )
    {
      if (recHitsEB.isValid() && (recHitsEB->size() > 0))
	{
	  
	  vector<EBDetId> maskedIds;
	  maskedIds.clear();
	  
	  if(isMaskEB_)
	    {
	      //Read the mask
	      std::ifstream mask_file(maskEBFile_.c_str());
	      if (mask_file.is_open()) 
		{
		  
		  char buffer[100];
		  int ieta, iphi, ret;
		  while( ! mask_file.getline(buffer,100).eof() ){
		    ret=sscanf(buffer,"%d %d ", &ieta, &iphi);
		    if ((ret!=2)) 
		      {
			//cout<< " Mask file "<<maskEBFile_.c_str()<<" contains wrong info - must be two digits per line "<<endl;
			continue;
		      }
		    // cout <<" EB Mask: ieta iphi "<<ieta<<" "<<iphi<<endl;
		    
		    EBDetId maskId(ieta,iphi); 
		    maskedIds.push_back(maskId);
		    
		  } // end buffer reading
		  
		  mask_file.close();
		} else {
		//cout<< " EB Mask File "<<maskEBFile_.c_str()<<" was not opened properly. Status: "<<mask_file.is_open()<<endl; 
	      }//mask_file.isopen
    
	      //cout<< " EB Pi0: # of masked Xtals: "<<maskedIds.size()<<endl;
	      
	    }// isMaskEB

	  

	  std::vector<EcalRecHit> seeds;
	  seeds.clear();
	  
	  vector<EBDetId> usedXtals;
	  usedXtals.clear();
	  
	  std::vector<EBDetId> detIdEBRecHits;
	  detIdEBRecHits.clear(); 
	  
	  std::vector<EcalRecHit> EBRecHits; 
	  EBRecHits.clear();
	  
	  float etot =0;
	  EcalRecHitCollection::const_iterator itb;
	  
	  // Build Seeds
	  for(itb=recHitsEB->begin(); itb!=recHitsEB->end(); ++itb){
	    
	    if ( (useRecoFlag_) && (itb->recoFlag() != 0) ) continue;
	    
	    EBDetId id(itb->id());
	    double energy = itb->energy();
	    if( energy < seleXtalMinEnergy_) continue; 
	    
	    EBDetId det = itb->id();
	    
	    if(isMaskEB_)
	      {
		
		detIdEBRecHits.push_back(det);
		EBRecHits.push_back(*itb);
		
		bool  maskedId=false;
		std::vector<EBDetId>::const_iterator mIds;
		for(mIds=maskedIds.begin(); mIds!=maskedIds.end(); mIds++){
		  if(*mIds==id){
		    //cout<< " EB Seed was masked : E ieta iphi recoFlag "<<itb->energy()<<" "<<id.ieta()<<" "<<id.iphi()<<" "<<itb->recoFlag()<<endl;
		    maskedId = true;
		    break;
		  }
		}
		
		if ( (energy > clusSeedThr_) && (!(maskedId)) ) seeds.push_back(*itb);
		
		if(!(maskedId))
		  {
		    etot+= itb->energy();	 
		  }
		
	      } 
	    else
	      {
		
		detIdEBRecHits.push_back(det);
		EBRecHits.push_back(*itb);
		if ( (energy > clusSeedThr_) ) seeds.push_back(*itb);
		
		etot+= itb->energy();	 
		
	      }
	    
	  } // Eb rechits
	  
	  
	  int nClus;
	  vector<float> eClus;
	  vector<float> etClus;
	  vector<float> etaClus; 
	  vector<float> thetaClus;
	  vector<float> phiClus;
	  vector<EBDetId> max_hit;
	  
	  vector< vector<EcalRecHit> > RecHitsCluster;
	  vector< vector<EcalRecHit> > RecHitsCluster5x5;
	  vector<float> s4s9Clus;
	  vector<float> s9s25Clus;
	  
	  
	  nClus=0;
	  
	  // Make own simple clusters (3x3, 5x5 or clusPhiSize_ x clusEtaSize_)
	  sort(seeds.begin(), seeds.end(), ecalRecHitLess());
	  
	  //Cycle on seeds
	  for (std::vector<EcalRecHit>::iterator itseed=seeds.begin(); itseed!=seeds.end(); itseed++) {
	    EBDetId seed_id = itseed->id();
	    std::vector<EBDetId>::const_iterator usedIds;
	    
	    bool seedAlreadyUsed=false;
	    for(usedIds=usedXtals.begin(); usedIds!=usedXtals.end(); usedIds++){
	      if(*usedIds==seed_id){
		seedAlreadyUsed=true;
		//cout<< " Seed with energy "<<itseed->energy()<<" was used !"<<endl;
		break; 
	      }
	    }
	    if(seedAlreadyUsed)continue;
	    
	    std::vector<DetId> clus_v = topology_EB->getWindow(seed_id,clusEtaSize_,clusPhiSize_);       
	    std::vector<std::pair<DetId,float> > clus_used;
	    
	    vector<EcalRecHit> RecHitsInWindow;
	    vector<EcalRecHit> RecHitsInWindow5x5;
	    
	    double simple_energy = 0; 
	    
	    for (std::vector<DetId>::iterator det=clus_v.begin(); det!=clus_v.end(); det++) {
	      EBDetId EBdet = *det;
	      //      cout<<" det "<< EBdet<<" ieta "<<EBdet.ieta()<<" iphi "<<EBdet.iphi()<<endl;
	      bool  HitAlreadyUsed=false;
	      for(usedIds=usedXtals.begin(); usedIds!=usedXtals.end(); usedIds++){
		if(*usedIds==*det){
		  HitAlreadyUsed=true;
		  break;
		}
	      }
	      if(HitAlreadyUsed)continue;
	      
	      
	      std::vector<EBDetId>::iterator itdet = find( detIdEBRecHits.begin(),detIdEBRecHits.end(),EBdet);
	      if(itdet == detIdEBRecHits.end()) continue; 
	      
	      int nn = int(itdet - detIdEBRecHits.begin());
	      usedXtals.push_back(*det);
	      RecHitsInWindow.push_back(EBRecHits[nn]);
	      clus_used.push_back(std::make_pair(*det,1));
	      simple_energy = simple_energy + EBRecHits[nn].energy();
	      
	    }
	    
	    if(simple_energy <= 0) continue;
	    
	    math::XYZPoint clus_pos = posCalculator_.Calculate_Location(clus_used,theBarrelEcalRecHits,geometry_EB,geometry_ES);
	    
	    float theta_s = 2. * atan(exp(-clus_pos.eta()));
	    float et_s = simple_energy * sin(theta_s);
	    
	    //Compute S4/S9 variable
	    //We are not sure to have 9 RecHits so need to check eta and phi:
	    ///check s4s9
	    float s4s9_tmp[4];
	    for(int i=0;i<4;i++)s4s9_tmp[i]= 0;
	    
	    int seed_ieta = seed_id.ieta();
	    int seed_iphi = seed_id.iphi();
	    
	    convxtalid( seed_iphi,seed_ieta);
	    
	    float e3x3 = 0; 
	    float e5x5 = 0; 
	    
	    for(unsigned int j=0; j<RecHitsInWindow.size();j++){
	      EBDetId det = (EBDetId)RecHitsInWindow[j].id(); 
	      
	      int ieta = det.ieta();
	      int iphi = det.iphi();
  
	      convxtalid(iphi,ieta);
  
	      float en = RecHitsInWindow[j].energy(); 
  
	      int dx = diff_neta_s(seed_ieta,ieta);
	      int dy = diff_nphi_s(seed_iphi,iphi);
  
	      if(abs(dx)<=1 && abs(dy)<=1) {
		e3x3 += en; 
		if(dx <= 0 && dy <=0) s4s9_tmp[0] += en; 
		if(dx >= 0 && dy <=0) s4s9_tmp[1] += en; 
		if(dx <= 0 && dy >=0) s4s9_tmp[2] += en; 
		if(dx >= 0 && dy >=0) s4s9_tmp[3] += en; 
	      }

  
	    }
    
	    if(e3x3 <= 0) continue;

	    float s4s9_max = *max_element( s4s9_tmp,s4s9_tmp+4)/e3x3; 

	    ///calculate e5x5
	    std::vector<DetId> clus_v5x5 = topology_EB->getWindow(seed_id,5,5); 
	    for( std::vector<DetId>::const_iterator idItr = clus_v5x5.begin(); idItr != clus_v5x5.end(); idItr++){
	      EBDetId det = *idItr;
      
	      //inside collections
	      std::vector<EBDetId>::iterator itdet = find( detIdEBRecHits.begin(),detIdEBRecHits.end(),det);
	      if(itdet == detIdEBRecHits.end()) continue; 

	      int nn = int(itdet - detIdEBRecHits.begin());
      
	      RecHitsInWindow5x5.push_back(EBRecHits[nn]);
	      e5x5 += EBRecHits[nn].energy();
      
	    }



	    if(e5x5 <= 0) continue;

	    eClus.push_back(simple_energy);
	    etClus.push_back(et_s);
	    etaClus.push_back(clus_pos.eta());
	    thetaClus.push_back(theta_s);
	    phiClus.push_back(clus_pos.phi());
	    s4s9Clus.push_back(s4s9_max);
	    s9s25Clus.push_back(e3x3/e5x5);
	    RecHitsCluster.push_back(RecHitsInWindow);
	    RecHitsCluster5x5.push_back(RecHitsInWindow5x5);

	    nClus++;
    
	  }
    
	  // Selection, based on Simple clustering
	  //pi0 candidates
	  int npi0_s=0;

	  for(Int_t i=0 ; i<nClus ; i++){
	    for(Int_t j=i+1 ; j<nClus ; j++){
	      //cout<<" i "<<i<<"  etClus[i] "<<etClus[i]<<" j "<<j<<"  etClus[j] "<<etClus[j]<<endl;
	      if( etClus[i]>selePtGamma_ && etClus[j]>selePtGamma_ && s4s9Clus[i]>seleS4S9Gamma_ && s4s9Clus[j]>seleS4S9Gamma_){


		float p0x = etClus[i] * cos(phiClus[i]);
		float p1x = etClus[j] * cos(phiClus[j]);
		float p0y = etClus[i] * sin(phiClus[i]);
		float p1y = etClus[j] * sin(phiClus[j]);
		float p0z = eClus[i] * cos(thetaClus[i]);
		float p1z = eClus[j] * cos(thetaClus[j]);
      
      
		float pt_pair = sqrt( (p0x+p1x)*(p0x+p1x) + (p0y+p1y)*(p0y+p1y));

		if (pt_pair < selePtPi0_)continue;

		float m_inv = sqrt ( (eClus[i] + eClus[j])*(eClus[i] + eClus[j]) - (p0x+p1x)*(p0x+p1x) - (p0y+p1y)*(p0y+p1y) - (p0z+p1z)*(p0z+p1z) );  
		if ( (m_inv<seleMinvMaxPi0_) && (m_inv>seleMinvMinPi0_) ){

		  //New Loop on cluster to measure isolation:
		  vector<int> IsoClus;
		  IsoClus.clear();
		  float Iso = 0;
		  TVector3 pairVect = TVector3((p0x+p1x), (p0y+p1y), (p0z+p1z));
		  for(Int_t k=0 ; k<nClus ; k++){


		    if(etClus[k] < ptMinForIsolation_) continue;

		    if(k==i || k==j)continue;
		    TVector3 ClusVect = TVector3(etClus[k] *cos(phiClus[k]), etClus[k] * sin(phiClus[k]) , eClus[k] * cos(thetaClus[k]));

		    float dretacl = fabs(etaClus[k] - pairVect.Eta());
		    float drcl = ClusVect.DeltaR(pairVect);
		    //cout<< "   Iso: k, E, drclpi0, detaclpi0, dphiclpi0 "<<k<<" "<<eClus[k]<<" "<<drcl<<" "<<dretacl<<endl;
		    if((drcl<selePi0BeltDR_) && (dretacl<selePi0BeltDeta_) ){
		      //cout<< "   ... good iso cluster #: "<<k<<" etClus[k] "<<etClus[k] <<endl;
		      Iso = Iso + etClus[k];
		      IsoClus.push_back(k);
		    }
		  }

		  if(Iso/pt_pair<selePi0Iso_){

		    h_Pi0_EB_eta->Fill(pairVect.Eta());
		    h_Pi0_EB_phi->Fill(pairVect.Phi());
		    
		    int order[2] = {0,0};

		    if(etClus[i] > etClus[j]) {
		      order[0] = i;
		      order[1] = j;
		    }
		    if(etClus[j] > etClus[i]){
		      order[0] = j;
		      order[1] = i;
		    }
		    
		    h_Pi0_EB_mass->Fill(m_inv);
		    h_Pi0_EB_pt1 ->Fill(etClus[order[0]]);
		    h_Pi0_EB_pt2 ->Fill(etClus[order[1]]);
		    h_Pi0_EB_pt  ->Fill(pt_pair);
		    
		    npi0_s++;
		    
		  }
		}
	      }
	    } // End of the "j" loop over Simple Clusters
	  } // End of the "i" loop over Simple Clusters
	}
      
   }// End isMonEBpi0_

  //------------------ End of pi0 in EB --------------------------// 
}


// ------------ Pi0 peak  ECAL BARREL ------------

void EcalValidation::doPi0Endcap ( const CaloGeometry *geometry,
				   const CaloTopology *topology,
				   edm::Handle<EcalRecHitCollection> recHitsEE ) 
{
  const CaloSubdetectorGeometry *geometry_EE = geometry->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);
  const CaloSubdetectorGeometry *geometry_ES = geometry->getSubdetectorGeometry(DetId::Ecal,EcalPreshower);

  const CaloSubdetectorTopology *topology_EE = topology->getSubdetectorTopology(DetId::Ecal,EcalEndcap);
  
  const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;

  if (isMonEEpi0_ )
    {
      if (recHitsEE.isValid() && (recHitsEE->size() > 0))
	{
  	  vector<EEDetId> maskedIds;
	  maskedIds.clear();
	  
	  if(isMaskEE_)
	    {
	      //Read the mask
	      std::ifstream mask_file(maskEEFile_.c_str());
	      if (mask_file.is_open()) 
		{
		  char buffer[100];
		  int ix, iy, iz, ret;
		  while( ! mask_file.getline(buffer,100).eof() ){
		    ret=sscanf(buffer,"%d %d %d", &iz, &ix, &iy);
		    if ((ret!=3)) 
		      {
			//cout<< " Mask file "<<maskEEFile_.c_str()<<" contains wrong info - must be three digits per line "<<endl;
			continue;
		      }
		    //cout <<" EE Mask: iz ix iy "<<iz<<" "<<ix<<" "<<iy<<endl;
		    
		    EEDetId maskId(ix,iy,iz); 
		    maskedIds.push_back(maskId);
		    
		  } // end buffer reading
		  
		  mask_file.close();
		} else {
		//cout<< " EE Mask File "<<maskEEFile_.c_str()<<" was not opened properly. Status: "<<mask_file.is_open()<<endl; 
	      }//mask_file.isopen
	      
	      //cout<< " EE Pi0: # of masked Xtals: "<<maskedIds.size()<<endl;
	      
	    }// isMaskEE
  
  
	  std::vector<EcalRecHit> seeds_EE;
	  seeds_EE.clear();
  
	  vector<EEDetId> usedXtals_EE;
	  usedXtals_EE.clear();
  
	  std::vector<EEDetId> detIdEERecHits;
	  detIdEERecHits.clear(); 
  
	  std::vector<EcalRecHit> EERecHits; 
	  EERecHits.clear();
  
	  float etot_EE =0;
	  EERecHitCollection::const_iterator ite;
  
	  // Build Seeds
	  for(ite=recHitsEE->begin(); ite!=recHitsEE->end(); ++ite){
    
	    if ( (useRecoFlag_) && (ite->recoFlag() != 0) ) continue;
  
	    double energy = ite->energy();
	    if( energy < seleXtalMinEnergy_EE_) continue; 
  
	    EEDetId det = ite->id();
	    EEDetId id(ite->id());
    
	    if(isMaskEE_)
	      {
  
		detIdEERecHits.push_back(det);
		EERecHits.push_back(*ite);
  
		bool  maskedId=false;
		std::vector<EEDetId>::const_iterator mIds;
		for(mIds=maskedIds.begin(); mIds!=maskedIds.end(); mIds++){
		  if(*mIds==id){
		    //cout<< " EE Seed was masked : E ix iy iz "<<ite->energy()<<" "<<id.ix()<<" "<<id.iy()<<" "<<id.zside()<<endl;
		    maskedId = true;
		    break;
		  }
		}
    
		if ( (energy > clusSeedThr_EE_) && (!(maskedId)) ) seeds_EE.push_back(*ite);
  
		if(!(maskedId))
		  {
		    etot_EE+= ite->energy();	 
		  }
  
      
	      }
	    else
	      {
		detIdEERecHits.push_back(det);
		EERecHits.push_back(*ite);
  
		if ( (energy > clusSeedThr_EE_) ) seeds_EE.push_back(*ite);
      
		etot_EE+= ite->energy();	 
      
	      } //ismaskEE
	  } // Ee rechits
    
  
	  int nClus_EE;
	  vector<float> eClus_EE;
	  vector<float> etClus_EE;
	  vector<float> etaClus_EE; 
	  vector<float> thetaClus_EE;
	  vector<float> phiClus_EE;
	  vector< vector<EcalRecHit> > RecHitsCluster_EE;
	  vector< vector<EcalRecHit> > RecHitsCluster_EE_5x5;
	  vector<float> s4s9Clus_EE;
	  vector<float> s9s25Clus_EE;
	  vector<float> xClusEndCap;
	  vector<float> yClusEndCap;
	  vector<float> zClusEndCap;

  
  
	  nClus_EE=0;
  
	  // Make own simple clusters (3x3, 5x5 or clusPhiSize_ x clusEtaSize_)
	  sort(seeds_EE.begin(), seeds_EE.end(), ecalRecHitLess());
  
	  //Cycle on seeds_EE
	  for (std::vector<EcalRecHit>::iterator itseed=seeds_EE.begin(); itseed!=seeds_EE.end(); itseed++) {
	    EEDetId seed_id = itseed->id();
	    std::vector<EEDetId>::const_iterator usedIds;
  
	    bool seedAlreadyUsed=false;
	    for(usedIds=usedXtals_EE.begin(); usedIds!=usedXtals_EE.end(); usedIds++){
	      if(*usedIds==seed_id){
		seedAlreadyUsed=true;
		//cout<< " Seed with energy "<<itseed->energy()<<" was used !"<<endl;
		break; 
	      }
	    }
	    if(seedAlreadyUsed)continue;
  
	    std::vector<DetId> clus_v = topology_EE->getWindow(seed_id,clusEtaSize_,clusPhiSize_);       
	    std::vector<std::pair<DetId,float> > clus_used;
    
	    vector<EcalRecHit> RecHitsInWindow;
	    vector<EcalRecHit> RecHitsInWindow5x5;
  
	    double simple_energy = 0; 
    
	    for (std::vector<DetId>::iterator det=clus_v.begin(); det!=clus_v.end(); det++) {
	      EEDetId EEdet = *det;
	      //      cout<<" det "<< EBdet<<" ieta "<<EBdet.ieta()<<" iphi "<<EBdet.iphi()<<endl;
	      bool  HitAlreadyUsed=false;
	      for(usedIds=usedXtals_EE.begin(); usedIds!=usedXtals_EE.end(); usedIds++){
		if(*usedIds==*det){
		  HitAlreadyUsed=true;
		  break;
		}
	      }
	      if(HitAlreadyUsed)continue;
  
  
	      std::vector<EEDetId>::iterator itdet = find( detIdEERecHits.begin(),detIdEERecHits.end(),EEdet);
	      if(itdet == detIdEERecHits.end()) continue; 
  
	      int nn = int(itdet - detIdEERecHits.begin());
	      usedXtals_EE.push_back(*det);
	      RecHitsInWindow.push_back(EERecHits[nn]);
	      clus_used.push_back(std::make_pair(*det,1));
	      simple_energy = simple_energy + EERecHits[nn].energy();
  
	    }
    
	    if(simple_energy <= 0) continue;
  
	    math::XYZPoint clus_pos = posCalculator_.Calculate_Location(clus_used, theEndcapEcalRecHits, geometry_EE, geometry_ES);
  
	    float theta_s = 2. * atan(exp(-clus_pos.eta()));
	    float et_s = simple_energy * sin(theta_s);
    
	    //Compute S4/S9 variable
	    //We are not sure to have 9 RecHits so need to check eta and phi:
	    ///check s4s9
	    float s4s9_tmp[4];
	    for(int i=0;i<4;i++)s4s9_tmp[i]= 0;
    
	    int ixSeed = seed_id.ix();
	    int iySeed = seed_id.iy();
	    float e3x3 = 0; 
	    float e5x5 = 0;
    
	    for(unsigned int j=0; j<RecHitsInWindow.size();j++){
	      EEDetId det_this = (EEDetId)RecHitsInWindow[j].id(); 
  
	      int dx = ixSeed - det_this.ix();
	      int dy = iySeed - det_this.iy();
  
	      float en = RecHitsInWindow[j].energy(); 
  
	      if(abs(dx)<=1 && abs(dy)<=1) {
		e3x3 += en; 
		if(dx <= 0 && dy <=0) s4s9_tmp[0] += en; 
		if(dx >= 0 && dy <=0) s4s9_tmp[1] += en; 
		if(dx <= 0 && dy >=0) s4s9_tmp[2] += en; 
		if(dx >= 0 && dy >=0) s4s9_tmp[3] += en; 
	      }
  
  
	    }
    
	    if(e3x3 <= 0) continue;
  
	    float s4s9_max = *max_element( s4s9_tmp,s4s9_tmp+4)/e3x3; 
  
	    ///calculate e5x5
	    std::vector<DetId> clus_v5x5 = topology_EE->getWindow(seed_id,5,5); 
	    for( std::vector<DetId>::const_iterator idItr = clus_v5x5.begin(); idItr != clus_v5x5.end(); idItr++){
	      EEDetId det = *idItr;
      
	      //inside collections
	      std::vector<EEDetId>::iterator itdet = find( detIdEERecHits.begin(),detIdEERecHits.end(),det);
	      if(itdet == detIdEERecHits.end()) continue; 
  
	      int nn = int(itdet - detIdEERecHits.begin());
      
	      RecHitsInWindow5x5.push_back(EERecHits[nn]);
	      e5x5 += EERecHits[nn].energy();
      
	    }
  
  
  
	    if(e5x5 <= 0) continue;
  
  
  
	    xClusEndCap.push_back(clus_pos.x());
	    yClusEndCap.push_back(clus_pos.y());
	    zClusEndCap.push_back(clus_pos.z());

	    eClus_EE.push_back(simple_energy);
	    etClus_EE.push_back(et_s);
	    etaClus_EE.push_back(clus_pos.eta());
	    thetaClus_EE.push_back(theta_s);
	    phiClus_EE.push_back(clus_pos.phi());
	    s4s9Clus_EE.push_back(s4s9_max);
	    s9s25Clus_EE.push_back(e3x3/e5x5);
	    RecHitsCluster_EE.push_back(RecHitsInWindow);
	    RecHitsCluster_EE_5x5.push_back(RecHitsInWindow5x5);
  
	    nClus_EE++;
    
	  }
    
	  // Selection, based on Simple clustering
	  //pi0 candidates
	  int npi0_s=0;
  
	  for(Int_t i=0 ; i<nClus_EE ; i++){
	    for(Int_t j=i+1 ; j<nClus_EE ; j++){
	      //cout<<" i "<<i<<"  etClus_EE[i] "<<etClus_EE[i]<<" j "<<j<<"  etClus_EE[j] "<<etClus_EE[j]<<endl;
	      if( s4s9Clus_EE[i] < seleS4S9Gamma_EE_ || s4s9Clus_EE[j] < seleS4S9Gamma_EE_) continue ;
  
	      float p0x = etClus_EE[i] * cos(phiClus_EE[i]);
	      float p1x = etClus_EE[j] * cos(phiClus_EE[j]);
	      float p0y = etClus_EE[i] * sin(phiClus_EE[i]);
	      float p1y = etClus_EE[j] * sin(phiClus_EE[j]);
	      float p0z = eClus_EE[i] * cos(thetaClus_EE[i]);
	      float p1z = eClus_EE[j] * cos(thetaClus_EE[j]);
      
      
	      float pt_pair = sqrt( (p0x+p1x)*(p0x+p1x) + (p0y+p1y)*(p0y+p1y));
  
	      if (pt_pair < selePtPi0_)continue;
  
	      float m_inv = sqrt ( (eClus_EE[i] + eClus_EE[j])*(eClus_EE[i] + eClus_EE[j]) - (p0x+p1x)*(p0x+p1x) - (p0y+p1y)*(p0y+p1y) - (p0z+p1z)*(p0z+p1z) );  
        
	      if ( m_inv > seleMinvMaxPi0_EE_ || m_inv < seleMinvMinPi0_EE_) continue;
        
	      ////try different cut for different regions
	      TVector3 pairVect = TVector3((p0x+p1x), (p0y+p1y), (p0z+p1z));
	      float etapair = fabs(pairVect.Eta());
	      float ptmin = etClus_EE[i] < etClus_EE[j] ?  etClus_EE[i] : etClus_EE[j]; 
  
	      if(etapair <= region1_Pi0_EE_){
		if(ptmin < selePtGammaPi0_EE_region1_ || pt_pair < selePtPi0_EE_region1_) continue; 
	      }else if( etapair <= region2_Pi0_EE_){
		if(ptmin < selePtGammaPi0_EE_region2_ || pt_pair < selePtPi0_EE_region2_) continue;
	      }else{
		if(ptmin < selePtGammaPi0_EE_region3_ || pt_pair < selePtPi0_EE_region3_) continue;
	      }
  
  
	      //New Loop on cluster to measure isolation:
	      vector<int> IsoClus_EE;
	      IsoClus_EE.clear();
	      float Iso = 0;
	      for(Int_t k=0 ; k<nClus_EE ; k++){
  
  
		if(etClus_EE[k] < ptMinForIsolation_EE_) continue;
  
		if(k==i || k==j)continue;
		TVector3 ClusVect = TVector3(etClus_EE[k] *cos(phiClus_EE[k]), etClus_EE[k] * sin(phiClus_EE[k]) , eClus_EE[k] * cos(thetaClus_EE[k]));
  
		float dretacl = fabs(etaClus_EE[k] - pairVect.Eta());
		float drcl = ClusVect.DeltaR(pairVect);
		//cout<< "   Iso: k, E, drclpi0, detaclpi0, dphiclpi0 "<<k<<" "<<eClus_EE[k]<<" "<<drcl<<" "<<dretacl<<endl;
		if((drcl<selePi0BeltDR_EE_) && (dretacl<selePi0BeltDeta_EE_) ){
		  //cout<< "   ... good iso cluster #: "<<k<<" etClus_EE[k] "<<etClus_EE[k] <<endl;
		  Iso = Iso + etClus_EE[k];
		  IsoClus_EE.push_back(k);
		}
	      }
  
	      if(Iso/pt_pair<selePi0Iso_EE_){
		
		h_Pi0_EE_eta->Fill(pairVect.Eta());
		h_Pi0_EE_phi->Fill(pairVect.Phi());
    
		int order[2] = {0,0};
		
		if(etClus_EE[i] > etClus_EE[j]) {
		    order[0] = i;
		    order[1] = j;
		}
		if(etClus_EE[j] > etClus_EE[i]) {
		    order[0] = j;
		    order[1] = i;
		}
		                      
		h_Pi0_EE_mass->Fill(m_inv);
		h_Pi0_EE_pt1 ->Fill(etClus_EE[order[0]]);
		h_Pi0_EE_pt2 ->Fill(etClus_EE[order[1]]);
		h_Pi0_EE_pt  ->Fill(pt_pair);
  
		npi0_s++;
          
	      }
	    } // End of the "j" loop over Simple Clusters
	  } // End of the "i" loop over Simple Clusters
	  
	}//------------------ End of pi0 in EE --------------------------//
      
    }//End isMonEEpi0_

  
  
}



// ----------additional functions-------------------

void EcalValidation::convxtalid(Int_t &nphi,Int_t &neta)
{
  // Barrel only
  // Output nphi 0...359; neta 0...84; nside=+1 (for eta>0), or 0 (for eta<0).
  // neta will be [-85,-1] , or [0,84], the minus sign indicates the z<0 side.
  
  if(neta > 0) neta -= 1;
  if(nphi > 359) nphi=nphi-360;
  
} //end of convxtalid

int EcalValidation::diff_neta_s(Int_t neta1, Int_t neta2){
  Int_t mdiff;
  mdiff=(neta1-neta2);
  return mdiff;
}

// Calculate the distance in xtals taking into account the periodicity of the Barrel
int EcalValidation::diff_nphi_s(Int_t nphi1,Int_t nphi2) {
  Int_t mdiff;
  if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2))) {
    mdiff=nphi1-nphi2;
  }
  else {
    mdiff=360-abs(nphi1-nphi2);
    if(nphi1>nphi2) mdiff=-mdiff;
  }
  return mdiff;
}

//define this as a plug-in
DEFINE_FWK_MODULE(EcalValidation);
