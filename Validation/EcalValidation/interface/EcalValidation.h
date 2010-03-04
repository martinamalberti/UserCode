#ifndef EcalValidation_h
#define EcalValidation_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

// ROOT include
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"


// Less than operator for sorting EcalRecHits according to energy.
class ecalRecHitLess : public std::binary_function<EcalRecHit, EcalRecHit, bool>
{
public:
  bool operator()(EcalRecHit x, EcalRecHit y)
  {
    return (x.energy() > y.energy());
  }
};


class EcalValidation : public edm::EDAnalyzer {
  
      public:
         explicit EcalValidation(const edm::ParameterSet&);
	 ~EcalValidation();
  
  
      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
	 virtual void endJob() ;

	 // ----------additional functions-------------------
	 void convxtalid(int & , int &);
	 int diff_neta_s(int,int);
	 int diff_nphi_s(int,int);

      protected:
	 void doPi0Barrel(const CaloGeometry *theCaloGeometry,
			  const CaloTopology *theCaloTopology,
			  edm::Handle<EcalRecHitCollection> recHitsEB);

	 void doPi0Endcap(const CaloGeometry *theCaloGeometry,
			  const CaloTopology *theCaloTopology,
			  edm::Handle<EcalRecHitCollection> recHitsEE);
  

	 // ----------member data ---------------------------
	 edm::InputTag recHitCollection_EB_;
	 edm::InputTag recHitCollection_EE_;
	 edm::InputTag basicClusterCollection_EB_;
	 edm::InputTag basicClusterCollection_EE_;
	 edm::InputTag superClusterCollection_EB_;
	 edm::InputTag superClusterCollection_EE_;
	 edm::InputTag esRecHitCollection_;
	 edm::InputTag esClusterCollectionX_ ;
	 edm::InputTag esClusterCollectionY_ ;
	 
	 double ethrEB_;
	 double ethrEE_;
	 

	 // for Pi0
	 PositionCalc posCalculator_ ;
	 bool ParameterLogWeighted_;
	 double ParameterX0_;
	 double ParameterT0_barl_;
	 double ParameterT0_endc_;
	 double ParameterT0_endcPresh_;
	 double ParameterW0_;
	 
	 
	 double clusSeedThr_;
	 double clusSeedThr_EE_;
	 int clusEtaSize_;
	 int clusPhiSize_;
	 
	 double seleXtalMinEnergy_;
	 double seleXtalMinEnergy_EE_;
	 
	 bool isMaskEB_;
	 bool isMaskEE_;
	 
	 /// Files with xtals masked
	 std::string maskEBFile_;
	 std::string maskEEFile_;
	 
	 /// Use Reco Flag from RH
	 bool useRecoFlag_;
	 
	 // ... barrel
	 bool isMonEBpi0_;

	 double selePtGamma_;
	 double selePtPi0_;
	 double seleMinvMaxPi0_;
	 double seleMinvMinPi0_;
	 double seleS4S9Gamma_;
	 double selePi0BeltDR_;
	 double selePi0BeltDeta_;
	 double selePi0Iso_;
	 double ptMinForIsolation_;
	 
	 // ... endcap
	 bool isMonEEpi0_;
	 
	 double selePtGamma_EE_;
	 double selePtPi0_EE_;
	 double seleMinvMaxPi0_EE_;
	 double seleMinvMinPi0_EE_;
	 double seleS4S9Gamma_EE_;
	 double selePi0Iso_EE_;
	 double selePi0BeltDR_EE_;
	 double selePi0BeltDeta_EE_;
	 double ptMinForIsolation_EE_;
	 
	 double region1_Pi0_EE_;
	 double selePtGammaPi0_EE_region1_;
	 double selePtPi0_EE_region1_;
	 double region2_Pi0_EE_;
	 double selePtGammaPi0_EE_region2_;
	 double selePtPi0_EE_region2_;
	 double selePtGammaPi0_EE_region3_;
	 double selePtPi0_EE_region3_;



	 // ------------- HISTOGRAMS ------------------------------------

	 int naiveId_;
	 
	 TH1D *h_numberOfEvents;
	 
	 // RecHits ----------------------------------------------
	 // ... barrel 
	 TH1D *h_recHits_EB_size; 
	 TH1D *h_recHits_EB_energy;
	 TH1D *h_recHits_EB_energyMax;
	 TH1D *h_recHits_EB_time;
	 TH1D *h_recHits_EB_Chi2;
	 TH1D *h_recHits_EB_OutOfTimeChi2;
	 TH2D *h_recHits_EB_occupancy;
	 
	 // ... endcap
	 TH1D *h_recHits_EEP_size;
	 TH1D *h_recHits_EEP_energy;
	 TH1D *h_recHits_EEP_energyMax;
	 TH1D *h_recHits_EEP_time;
	 TH1D *h_recHits_EEP_Chi2;
	 TH1D *h_recHits_EEP_OutOfTimeChi2;
	 TH2D *h_recHits_EEP_occupancy;
	 
	 TH1D *h_recHits_EEM_size;
	 TH1D *h_recHits_EEM_energy;
	 TH1D *h_recHits_EEM_energyMax;
	 TH1D *h_recHits_EEM_time;
	 TH1D *h_recHits_EEM_Chi2;
	 TH1D *h_recHits_EEM_OutOfTimeChi2;
	 TH2D *h_recHits_EEM_occupancy;
	 
	 TH1D *h_recHits_EB_phi;
	 TH1D *h_recHits_EE_phi;
	 TH1D *h_recHits_eta;
	 
	 
	 // Basic Clusters ----------------------------------------------
	 
	 // ... barrel
	 TH1D *h_basicClusters_EB_size;
	 TH1D *h_basicClusters_EB_nXtals;
	 TH1D *h_basicClusters_EB_energy;
	 
	 // ... endcap
	 TH1D *h_basicClusters_EEP_size;
	 TH1D *h_basicClusters_EEP_nXtals;
	 TH1D *h_basicClusters_EEP_energy;
	 
	 TH1D *h_basicClusters_EEM_size;
	 TH1D *h_basicClusters_EEM_nXtals;
	 TH1D *h_basicClusters_EEM_energy;
	 
	 TH1D *h_basicClusters_eta;
	 TH1D *h_basicClusters_EB_phi;
	 TH1D *h_basicClusters_EE_phi;
	 
	 
	 // Super Clusters ----------------------------------------------
	 // ... barrel
	 TH1D *h_superClusters_EB_size;
	 TH1D *h_superClusters_EB_nXtals;
	 TH1D *h_superClusters_EB_nBC;
	 TH1D *h_superClusters_EB_energy;
	 TH1D *h_superClusters_EB_E1oE9;
	 
	 // ... endcap
	 TH1D *h_superClusters_EEP_size;
	 TH1D *h_superClusters_EEP_nXtals;
	 TH1D *h_superClusters_EEP_nBC;
	 TH1D *h_superClusters_EEP_energy;
	 TH1D *h_superClusters_EEP_E1oE9;
	 
	 TH1D *h_superClusters_EEM_size;
	 TH1D *h_superClusters_EEM_nXtals;
	 TH1D *h_superClusters_EEM_nBC;
	 TH1D *h_superClusters_EEM_energy;
	 TH1D *h_superClusters_EEM_E1oE9;
	 
	 TH1D *h_superClusters_eta;
	 TH1D *h_superClusters_EB_phi;
	 TH1D *h_superClusters_EE_phi;
	 
	 
	 // PRESHOWER ----------------------------------------------
	 
	 TH1D *h_recHits_ES_size;
	 TH1D *h_recHits_ES_energy;
	 TH1D *h_recHits_ES_energyMax;
	 TH1D *h_recHits_ES_time;
	 TH1D *h_recHits_ES_energy_F[2];
	 TH1D *h_recHits_ES_energy_R[2];
	 TH1D *h_esClusters_energy_plane1;
	 TH1D *h_esClusters_energy_plane2;
	 TH1D *h_esClusters_energy_ratio;
	 
	 // Pi0 peak ----------------------------------------------

	 TH1D *h_Pi0_EB_mass;
	 TH1D *h_Pi0_EB_pt1;
	 TH1D *h_Pi0_EB_pt2;
	 TH1D *h_Pi0_EB_pt;
	 TH1D *h_Pi0_EB_eta;
	 TH1D *h_Pi0_EB_phi;

	 TH1D *h_Pi0_EE_mass;
	 TH1D *h_Pi0_EE_pt1;
	 TH1D *h_Pi0_EE_pt2;
	 TH1D *h_Pi0_EE_pt;
	 TH1D *h_Pi0_EE_eta;
	 TH1D *h_Pi0_EE_phi;


};


#endif
