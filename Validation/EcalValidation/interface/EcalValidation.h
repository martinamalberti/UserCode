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

// ROOT include
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"



class EcalValidation : public edm::EDAnalyzer {
        public:
                explicit EcalValidation(const edm::ParameterSet&);
                ~EcalValidation();


        private:
                virtual void beginJob() ;
                virtual void analyze(const edm::Event&, const edm::EventSetup&);
                virtual void endJob() ;

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


		// PRESHOWER

		TH1D *h_recHits_ES_size;
		TH1D *h_recHits_ES_energy;
		TH1D *h_recHits_ES_energyMax;
		TH1D *h_recHits_ES_time;
		TH1D *h_recHits_ES_energy_F[2];
		TH1D *h_recHits_ES_energy_R[2];
		TH1D *h_esClusters_energy_plane1;
		TH1D *h_esClusters_energy_plane2;
		TH1D *h_esClusters_energy_ratio;


};


#endif
