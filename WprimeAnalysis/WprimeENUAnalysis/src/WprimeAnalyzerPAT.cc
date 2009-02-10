// -*- C++ -*-
//
// Package:    WprimeAnalyzerPAT
// Class:      WprimeAnalyzerPAT
// 
/**\class WprimeAnalyzerPAT WprimeAnalyzerPAT.cc WprimeAnalysis/WprimeENUAnalysis/src/WprimeAnalyzerPAT.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Martina Malberti
//         Created:  Wed Nov 26 14:30:50 CET 2008
// $Id$
//
//


// system include files
#include <memory>
#include<string>
#include "math.h"
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Flags.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <Math/GenVector/VectorUtil.h>
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "PhysicsTools/Utilities/interface/deltaR.h"
#include <CLHEP/Vector/LorentzVector.h>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeAnalyzerPAT.h"


#define PI 3.141592654
#define TWOPI 6.283185308

using namespace edm;
using namespace std;



WprimeAnalyzerPAT::WprimeAnalyzerPAT(const edm::ParameterSet& iConfig)

{
  eleLabel_      = iConfig.getParameter<edm::InputTag>("electronTag");
  metLabel_      = iConfig.getParameter<edm::InputTag>("metTag");
  jetLabel_      = iConfig.getParameter<edm::InputTag>("jetTag");
  hlTriggerResults_   = iConfig.getParameter<edm::InputTag>("TriggerResults");

  electronID_         = iConfig.getUntrackedParameter<std::string>("electronID") ;

  
  EvtCounter = 0 ;

}



WprimeAnalyzerPAT::~WprimeAnalyzerPAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void 
WprimeAnalyzerPAT::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;

  //---- Tree creation ----
  tTreeUtilities_ = fs->make<TTree>("tTreeUtilities","tTreeUtilities");
  MCelePx_  = new std::vector<double>; // MCtruth
  MCelePy_  = new std::vector<double>; // MCtruth
  MCelePz_  = new std::vector<double>; // MCtruth
  MCeleEta_ = new std::vector<double>;    // MCtruth
  MCelePhi_ = new std::vector<double>;    // MCtruth
  MCelePid_ = new std::vector<int>;    // MCtruth
  
  isMCmatched_  = new std::vector<int>; // MC truth matching 
  elePx_  = new std::vector<double>; // track momentum 
  elePy_  = new std::vector<double>; // track momentum 
  elePz_  = new std::vector<double>; // track momentum 
  eleE_ = new std::vector<double>;    // SC energy
  eleEt_ = new std::vector<double>;   // SC transverse energy
  eleEta_ = new std::vector<double>;  // SC pseudorapidity
  elePhi_ = new std::vector<double>;  // SC phi
  eleCharge_ = new std::vector<float>;  // electron charge
  eleId_ = new std::vector<float>;      // electron ID
  eleTrkIsol_ = new std::vector<double>; // track isolation
  eleEcalIsol_ = new std::vector<double>; // ecal isolation
  eleHcalIsolD1_ = new std::vector<double>; // hcal isolation
  eleHcalIsolD2_ = new std::vector<double>; // hcal isolation

  jetPx_  = new std::vector<double>; // 
  jetPy_  = new std::vector<double>; // 
  jetPz_  = new std::vector<double>; // 
  jetPt_  = new std::vector<double>; // 
  jetEta_  = new std::vector<double>; // 
  jetPhi_  = new std::vector<double>; // 

  tTreeUtilities_->Branch("MCelePx_","std::vector<double>",&MCelePx_);
  tTreeUtilities_->Branch("MCelePy_","std::vector<double>",&MCelePy_);
  tTreeUtilities_->Branch("MCelePz_","std::vector<double>",&MCelePz_);
  tTreeUtilities_->Branch("MCeleEta_","std::vector<double>",&MCeleEta_);
  tTreeUtilities_->Branch("MCelePhi_","std::vector<double>",&MCelePhi_);
  tTreeUtilities_->Branch("MCelePid_","std::vector<int>",&MCelePid_);

  tTreeUtilities_->Branch("NeleCand_",&NeleCand_,"NeleCand_/I");
  tTreeUtilities_->Branch("isMCmatched_","std::vector<int>",&isMCmatched_);
  tTreeUtilities_->Branch("elePx_","std::vector<double>",&elePx_);
  tTreeUtilities_->Branch("elePy_","std::vector<double>",&elePy_);
  tTreeUtilities_->Branch("elePz_","std::vector<double>",&elePz_);
  tTreeUtilities_->Branch("eleE_","std::vector<double>",&eleE_);
  tTreeUtilities_->Branch("eleEt_","std::vector<double>",&eleEt_);
  tTreeUtilities_->Branch("eleEta_","std::vector<double>",&eleEta_);
  tTreeUtilities_->Branch("elePhi_","std::vector<double>",&elePhi_);
  tTreeUtilities_->Branch("eleId_","std::vector<float>",&eleId_);
  tTreeUtilities_->Branch("eleCharge_","std::vector<float>",&eleCharge_);
  tTreeUtilities_->Branch("eleTrkIsol_","std::vector<double>",&eleTrkIsol_);
  tTreeUtilities_->Branch("eleEcalIsol_","std::vector<double>",&eleEcalIsol_);
  tTreeUtilities_->Branch("eleHcalIsolD1_","std::vector<double>",&eleHcalIsolD1_);
  tTreeUtilities_->Branch("eleHcalIsolD2_","std::vector<double>",&eleHcalIsolD2_);
  
  tTreeUtilities_->Branch("met_",&met_,"met_/D");
  tTreeUtilities_->Branch("mex_",&mex_,"mex_/D");
  tTreeUtilities_->Branch("mey_",&mey_,"mey_/D");
  tTreeUtilities_->Branch("metPhi_",&metPhi_,"metPhi_/D");

  tTreeUtilities_->Branch("uncorrMet_",&uncorrMet_,"uncorrMet_/D");
  tTreeUtilities_->Branch("uncorrMex_",&uncorrMex_,"uncorrMex_/D");
  tTreeUtilities_->Branch("uncorrMey_",&uncorrMey_,"uncorrMey_/D");
  tTreeUtilities_->Branch("uncorrMetPhi_",&uncorrMetPhi_,"uncorrMetPhi_/D");

  tTreeUtilities_->Branch("Njets_",&Njets_,"Njets_/I");
  tTreeUtilities_->Branch("jetPx_","std::vector<double>",&jetPx_);
  tTreeUtilities_->Branch("jetPy_","std::vector<double>",&jetPy_);
  tTreeUtilities_->Branch("jetPz_","std::vector<double>",&jetPz_);
  tTreeUtilities_->Branch("jetPt_","std::vector<double>",&jetPt_);
  tTreeUtilities_->Branch("jetEta_","std::vector<double>",&jetEta_);
  tTreeUtilities_->Branch("jetPhi_","std::vector<double>",&jetPhi_);

  tTreeUtilities_->Branch("HLTLooseIsoEle15_",&HLTLooseIsoEle15_,"HLTLooseIsoEle15_/I");			  
  tTreeUtilities_->Branch("HLT_EM80_",&HLT_EM80_,"HLT_EM80_/I");			  
  tTreeUtilities_->Branch("HLT_EM200_",&HLT_EM200_,"HLT_EM200_/I");			  
  tTreeUtilities_->Branch("HLT_Photon15_",&HLT_Photon15_,"HLT_Photon15_/I");
  tTreeUtilities_->Branch("HLT_Photon25_",&HLT_Photon25_,"HLT_Photon25_/I");
			  

}


// ------------ method called to for each event  ------------
void
WprimeAnalyzerPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     
  EvtCounter++;
 
  //---- clear the vectors ----
  MCelePx_ ->clear(); 
  MCelePy_ ->clear(); 
  MCelePz_ ->clear(); 
  MCeleEta_ ->clear(); 
  MCelePhi_ ->clear(); 
  MCelePid_ ->clear(); 
 
  NeleCand_    = 0;
  isMCmatched_ ->clear(); 
  elePx_ ->clear(); 
  elePy_ ->clear(); 
  elePz_ ->clear(); 
  eleE_  ->clear();
  eleEt_ ->clear();
  eleEta_->clear();
  elePhi_->clear();
  eleId_ ->clear();
  eleCharge_ ->clear();
  eleTrkIsol_ ->clear();
  eleEcalIsol_->clear();
  eleHcalIsolD1_->clear();
  eleHcalIsolD2_->clear();

  met_    = 0 ;
  mex_    = 0 ;
  mey_    = 0 ;
  metPhi_ = 0;

  uncorrMet_    = 0 ;
  uncorrMex_    = 0 ;
  uncorrMey_    = 0 ;
  uncorrMetPhi_ = 0;

  Njets_ = 0;
  jetPx_ ->clear(); 
  jetPy_ ->clear(); 
  jetPz_ ->clear(); 
  jetPt_ ->clear(); 
  jetEta_ ->clear(); 
  jetPhi_ ->clear(); 


  HLTLooseIsoEle15_ = 0;
  HLT_EM80_         = 0;
  HLT_EM200_        = 0;
  HLT_Photon15_     = 0;
  HLT_Photon25_     = 0;

  // Read in electrons
  Handle<View<pat::Electron> > electronHandle;
  iEvent.getByLabel(eleLabel_,electronHandle);
  View<pat::Electron> electrons = *electronHandle;

 
  
  NeleCand_ = electrons.size();
  // Loop over electrons
  for ( unsigned int i=0; i<electrons.size(); ++i ) {
  
    pat::Electron electron = electrons.at(i);

    const reco::Particle * mcMatch =  electron.genLepton();
    if (mcMatch) {
      isMCmatched_->push_back(1);
      MCelePx_  -> push_back( mcMatch->px() );
      MCelePy_  -> push_back( mcMatch->py() );
      MCelePz_  -> push_back( mcMatch->pz() );
      MCeleEta_ -> push_back( mcMatch->eta() );
      MCelePhi_ -> push_back( mcMatch->phi() );
      MCelePid_ -> push_back( mcMatch->pdgId() );
    }
    else {
      isMCmatched_->push_back(0);
      MCelePx_  -> push_back( -9999.);
      MCelePy_  -> push_back( -9999.);
      MCelePz_  -> push_back( -9999.);
      MCeleEta_ -> push_back( -9999.);
      MCelePhi_ -> push_back( -9999 );
      MCelePid_ -> push_back( -9999 );
    }


    reco::SuperClusterRef scRef = electron.superCluster();
    double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
    double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());

    elePx_ -> push_back( electron.trackMomentumAtVtx().x() );
    elePy_ -> push_back( electron.trackMomentumAtVtx().y() );
    elePz_ -> push_back( electron.trackMomentumAtVtx().z() );
    eleE_  -> push_back( scRef->energy());
    eleEt_ -> push_back( scRef->energy() * (Rt/R) );
    eleEta_-> push_back( scRef->eta());
    elePhi_-> push_back( scRef->phi());
    eleCharge_ ->push_back( electron.gsfTrack()->charge() );
    
    eleId_ ->push_back( electron.electronID(electronID_));

    // need to add isolations
 
    
    
    
  }// end loop over electron candidates


  //********* MET  from Calo Towers 
  // 
  edm::Handle<edm::View<pat::MET> > metHandle;
  iEvent.getByLabel(metLabel_,metHandle);
  View<pat::MET>  mets = *metHandle;
  for ( unsigned int i=0; i<mets.size(); ++i ) {
    pat::MET met = mets.at(i);
    // corrected MET
    met_ = met.et();
    mex_ = met.px();
    mey_ = met.py();
    metPhi_ = met.phi(); 

    // Uncorrected MET
    uncorrMet_ = met.uncorrectedPt(pat::MET::uncorrALL);
    uncorrMex_ = met.corEx(pat::MET::uncorrALL);
    uncorrMey_ = met.corEy(pat::MET::uncorrALL);
    uncorrMetPhi_ = met.uncorrectedPhi(pat::MET::uncorrALL); 

  }

  //********* JETS
  Handle<View<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetLabel_,jetHandle);
  View<pat::Jet> jets = *jetHandle;
   
  Njets_ = jets.size();
  for ( unsigned int i=0; i<jets.size(); ++i ) {
    pat::Jet jet = jets.at(i);
    jetPx_->push_back(jet.px());
    jetPy_->push_back(jet.py());
    jetPz_->push_back(jet.pz());
    jetPt_->push_back(jet.pt());
    jetEta_->push_back(jet.eta());
    jetPhi_->push_back(jet.phi());
  }
  
  
  //********** HLT INFO 
  edm::TriggerNames triggerNames_;
  Handle<TriggerResults> hltresults;
  iEvent.getByLabel(hlTriggerResults_,hltresults);
  if (&hltresults) 
    {
      int ntrigs=(*hltresults).size();
      
      if (ntrigs==0){std::cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;}
      triggerNames_.init(*hltresults);
      hlNames_=triggerNames_.triggerNames();
      
      const unsigned int n(hlNames_.size());
      for (unsigned int i = 0; i < n; ++i)
	{
	  if(hlNames_[i]=="HLT_LooseIsoEle15_LW_L1R" && hltresults->accept(i)) HLTLooseIsoEle15_ = 1;
	  if(hlNames_[i]=="HLT_EM80" && hltresults->accept(i)) HLT_EM80_ = 1;
	  if(hlNames_[i]=="HLT_EM200" && hltresults->accept(i)) HLT_EM200_ = 1;
	  if(hlNames_[i]=="HLT_Photon15_L1R" && hltresults->accept(i)) HLT_Photon15_ = 1;
	  if(hlNames_[i]=="HLT_Photon25_L1R" && hltresults->accept(i)) HLT_Photon25_ = 1;	  
	}
    }

  //tree fill
  tTreeUtilities_->Fill();
}




// ------------ method called once each job just after ending the event loop  ------------
void 
WprimeAnalyzerPAT::endJob() {
  
  delete MCelePx_ ; 
  delete MCelePy_ ; 
  delete MCelePz_ ; 
  delete MCeleEta_ ; 
  delete MCelePhi_ ; 
  delete MCelePid_ ; 
 
 
  delete isMCmatched_ ; 
  delete elePx_ ; 
  delete elePy_ ; 
  delete elePz_ ; 
  delete eleE_  ;
  delete eleEt_ ;
  delete eleEta_;
  delete elePhi_;
  delete eleId_ ;
  delete eleCharge_ ;
  delete eleTrkIsol_ ;
  delete eleEcalIsol_;
  delete eleHcalIsolD1_;
  delete eleHcalIsolD2_;

  delete jetPx_ ; 
  delete jetPy_ ; 
  delete jetPz_ ; 
  delete jetPt_ ; 
  delete jetEta_ ; 
  delete jetPhi_ ; 


  cout << "Analyzed " << EvtCounter <<" events"<< endl;
  //cout << "Writing information into file: " << thefile->GetName() << endl;



}

