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

#include <CLHEP/Vector/LorentzVector.h>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeAnalyzerPAT.h"


#define PI 3.141592654
#define TWOPI 6.283185308



using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;



WprimeAnalyzerPAT::WprimeAnalyzerPAT(const edm::ParameterSet& iConfig)

{
  eleLabel_      = iConfig.getParameter<edm::InputTag>("electronTag");
  metLabel_      = iConfig.getParameter<edm::InputTag>("metTag");
  jetLabel_      = iConfig.getParameter<edm::InputTag>("jetTag");
  muonLabel_     = iConfig.getParameter<edm::InputTag>("muonTag");

  hlTriggerResults_   = iConfig.getParameter<edm::InputTag>("TriggerResults");
  
  electronID_         = iConfig.getUntrackedParameter<std::string>("electronID") ;
  btagAlgo_           = iConfig.getUntrackedParameter<std::string>("btagAlgo") ;
   
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
 
  MCelePx    = new std::vector<double>; // MCtruth
  MCelePy    = new std::vector<double>; // MCtruth
  MCelePz    = new std::vector<double>; // MCtruth
  MCeleEta   = new std::vector<double>;    // MCtruth
  MCelePhi   = new std::vector<double>;    // MCtruth
  MCelePid   = new std::vector<int>;    // MCtruth
  
  MCmuonPx   = new std::vector<double>; // MCtruth
  MCmuonPy   = new std::vector<double>; // MCtruth
  MCmuonPz   = new std::vector<double>; // MCtruth
  MCmuonEta  = new std::vector<double>;    // MCtruth
  MCmuonPhi  = new std::vector<double>;    // MCtruth
  MCmuonPid  = new std::vector<int>;    // MCtruth

  isMCmatched = new std::vector<int>; // MC truth matching 
  elePx       = new std::vector<double>; // track momentum 
  elePy       = new std::vector<double>; // track momentum 
  elePz       = new std::vector<double>; // track momentum 
  eleE        = new std::vector<double>;    // SC energy
  eleEt       = new std::vector<double>;   // SC transverse energy
  eleEta      = new std::vector<double>;  // ele pseudorapidity
  elePhi      = new std::vector<double>;  // ele phi
  eleCharge   = new std::vector<double>;  // electron charge
  eleId       = new std::vector<double>;      // electron ID

  eleSigmaIEtaIEta = new std::vector<double>;    // sigmaIEtaIEta
  eleE1x5          = new std::vector<double>;
  eleE2x5          = new std::vector<double>;
  eleE5x5          = new std::vector<double>;

  eleTrkIsol       = new std::vector<double>; // track isolation
  eleEcalIsol      = new std::vector<double>; // ecal isolation
  eleHcalIsolD1    = new std::vector<double>; // hcal isolation
  eleHcalIsolD2    = new std::vector<double>; // hcal isolation

  jetPx     = new std::vector<double>; // 
  jetPy     = new std::vector<double>; // 
  jetPz     = new std::vector<double>; // 
  jetPt     = new std::vector<double>; // 
  jetEta    = new std::vector<double>; // 
  jetPhi    = new std::vector<double>; // 
  jetBdisc  = new std::vector<double>; // 

  muonPt    = new std::vector<double>; 
  muonEta   = new std::vector<double>; 
  muonPhi   = new std::vector<double>; 


  tTreeUtilities_->Branch("MCelePx","std::vector<double>",&MCelePx);
  tTreeUtilities_->Branch("MCelePy","std::vector<double>",&MCelePy);
  tTreeUtilities_->Branch("MCelePz","std::vector<double>",&MCelePz);
  tTreeUtilities_->Branch("MCeleEta","std::vector<double>",&MCeleEta);
  tTreeUtilities_->Branch("MCelePhi","std::vector<double>",&MCelePhi);
  tTreeUtilities_->Branch("MCelePid","std::vector<int>",&MCelePid);

  tTreeUtilities_->Branch("MCmuonPx","std::vector<double>",&MCmuonPx);
  tTreeUtilities_->Branch("MCmuonPy","std::vector<double>",&MCmuonPy);
  tTreeUtilities_->Branch("MCmuonPz","std::vector<double>",&MCmuonPz);
  tTreeUtilities_->Branch("MCmuonEta","std::vector<double>",&MCmuonEta);
  tTreeUtilities_->Branch("MCmuonPhi","std::vector<double>",&MCmuonPhi);
  tTreeUtilities_->Branch("MCmuonPid","std::vector<int>",&MCmuonPid);

  tTreeUtilities_->Branch("NeleCand",&NeleCand,"NeleCand/I");
  tTreeUtilities_->Branch("isMCmatched","std::vector<int>",&isMCmatched);
  tTreeUtilities_->Branch("elePx","std::vector<double>",&elePx);
  tTreeUtilities_->Branch("elePy","std::vector<double>",&elePy);
  tTreeUtilities_->Branch("elePz","std::vector<double>",&elePz);
  tTreeUtilities_->Branch("eleE","std::vector<double>",&eleE);
  tTreeUtilities_->Branch("eleEt","std::vector<double>",&eleEt);
  tTreeUtilities_->Branch("eleEta","std::vector<double>",&eleEta);
  tTreeUtilities_->Branch("elePhi","std::vector<double>",&elePhi);
  tTreeUtilities_->Branch("eleId","std::vector<double>",&eleId);
 
  tTreeUtilities_->Branch("eleSigmaIEtaIEta","std::vector<double>",&eleSigmaIEtaIEta);
  tTreeUtilities_->Branch("eleE1x5","std::vector<double>",&eleE1x5);
  tTreeUtilities_->Branch("eleE2x5","std::vector<double>",&eleE2x5);
  tTreeUtilities_->Branch("eleE5x5","std::vector<double>",&eleE5x5);
 
  tTreeUtilities_->Branch("eleCharge","std::vector<double>",&eleCharge);
  tTreeUtilities_->Branch("eleTrkIsol","std::vector<double>",&eleTrkIsol);
  tTreeUtilities_->Branch("eleEcalIsol","std::vector<double>",&eleEcalIsol);
  tTreeUtilities_->Branch("eleHcalIsolD1","std::vector<double>",&eleHcalIsolD1);
  tTreeUtilities_->Branch("eleHcalIsolD2","std::vector<double>",&eleHcalIsolD2);
  
  tTreeUtilities_->Branch("Met",&Met,"Met/D");
  tTreeUtilities_->Branch("Mex",&Mex,"Mex/D");
  tTreeUtilities_->Branch("Mey",&Mey,"Mey/D");
  tTreeUtilities_->Branch("MetPhi",&MetPhi,"MetPhi/D");

  tTreeUtilities_->Branch("uncorrMet",&uncorrMet,"uncorrMet/D");
  tTreeUtilities_->Branch("uncorrMex",&uncorrMex,"uncorrMex/D");
  tTreeUtilities_->Branch("uncorrMey",&uncorrMey,"uncorrMey/D");
  tTreeUtilities_->Branch("uncorrMetPhi",&uncorrMetPhi,"uncorrMetPhi/D");

  tTreeUtilities_->Branch("Njets",&Njets,"Njets/I");
  tTreeUtilities_->Branch("jetPx","std::vector<double>",&jetPx);
  tTreeUtilities_->Branch("jetPy","std::vector<double>",&jetPy);
  tTreeUtilities_->Branch("jetPz","std::vector<double>",&jetPz);
  tTreeUtilities_->Branch("jetPt","std::vector<double>",&jetPt);
  tTreeUtilities_->Branch("jetEta","std::vector<double>",&jetEta);
  tTreeUtilities_->Branch("jetPhi","std::vector<double>",&jetPhi);
  tTreeUtilities_->Branch("jetBdisc","std::vector<double>",&jetBdisc);

  tTreeUtilities_->Branch("muonPt","std::vector<double>",&muonPt);
  tTreeUtilities_->Branch("muonEta","std::vector<double>",&muonEta);
  tTreeUtilities_->Branch("muonPhi","std::vector<double>",&muonPhi);

  tTreeUtilities_->Branch("HLTLooseIsoEle15",&HLTLooseIsoEle15,"HLTLooseIsoEle15/I");		  
  tTreeUtilities_->Branch("HLTEle15",&HLTEle15,"HLTEle15/I");			  
  tTreeUtilities_->Branch("HLTEle10",&HLTEle10,"HLTEle10/I");			  
  tTreeUtilities_->Branch("HLTEle20",&HLTEle20,"HLTEle20/I");			  
  tTreeUtilities_->Branch("HLTPhoton15",&HLTPhoton15,"HLTPhoton15/I");
  tTreeUtilities_->Branch("HLTPhoton25",&HLTPhoton25,"HLTPhoton25/I");
			  

}


// ------------ method called to for each event  ------------
void
WprimeAnalyzerPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     
  EvtCounter++;
 
  //---- clear the vectors ----
  MCelePx  -> clear(); 
  MCelePy  -> clear(); 
  MCelePz  -> clear(); 
  MCeleEta -> clear(); 
  MCelePhi -> clear(); 
  MCelePid -> clear(); 

  MCmuonPx  -> clear();
  MCmuonPy  -> clear();
  MCmuonPz  -> clear();
  MCmuonEta -> clear();
  MCmuonPhi -> clear();
  MCmuonPid -> clear();

  NeleCand    =  0;
  isMCmatched -> clear(); 
  elePx  -> clear(); 
  elePy  -> clear(); 
  elePz  -> clear(); 
  eleE   -> clear();
  eleEt  -> clear();
  eleEta -> clear();
  elePhi -> clear();
  eleCharge ->clear();

  eleId  -> clear();
  eleSigmaIEtaIEta -> clear();
  eleE1x5 -> clear();
  eleE2x5 -> clear();
  eleE5x5 -> clear();

  eleTrkIsol    -> clear();
  eleEcalIsol   -> clear();
  eleHcalIsolD1 -> clear();
  eleHcalIsolD2 -> clear();

  Met    = 0 ;
  Mex    = 0 ;
  Mey    = 0 ;
  MetPhi = 0;

  uncorrMet    = 0 ;
  uncorrMex    = 0 ;
  uncorrMey    = 0 ;
  uncorrMetPhi = 0;

  Njets = 0;
  jetPx    -> clear(); 
  jetPy    -> clear(); 
  jetPz    -> clear(); 
  jetPt    -> clear(); 
  jetEta   -> clear(); 
  jetPhi   -> clear(); 
  jetBdisc -> clear(); 

  muonPt  -> clear(); 
  muonEta -> clear(); 
  muonPhi -> clear(); 
  
  HLTLooseIsoEle15 = 0;
  HLTEle15         = 0; 
  HLTEle10         = 0; 
  HLTEle20         = 0; 
  HLTPhoton15      = 0;
  HLTPhoton25      = 0;
  
  // some Gen Level infos
  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  for(size_t i = 0; i < genParticles->size(); ++ i) {

    const GenParticle & p = (*genParticles)[i];
    if ( abs(p.pdgId()) == 11 && abs(p.mother(0)->pdgId())== 24 && p.status()==3){
      MCelePx  -> push_back( p.px() );
      MCelePy  -> push_back( p.py() );
      MCelePz  -> push_back( p.pz() );
      MCeleEta -> push_back( p.eta() );
      MCelePhi -> push_back( p.phi() );
      MCelePid -> push_back( p.pdgId() );
    }
    if ( abs(p.pdgId()) == 13 && abs(p.mother(0)->pdgId())== 24 && p.status()==3){
      MCmuonPx  -> push_back( p.px() );
      MCmuonPy  -> push_back( p.py() );
      MCmuonPz  -> push_back( p.pz() );
      MCmuonEta -> push_back( p.eta() );
      MCmuonPhi -> push_back( p.phi() );
      MCmuonPid -> push_back( p.pdgId() );
    }
  }


  // *********** ELECTRONS 
  Handle<View<pat::Electron> > electronHandle;
  iEvent.getByLabel(eleLabel_,electronHandle);
  View<pat::Electron> electrons = *electronHandle;

 
  NeleCand = electrons.size();
  
  // Loop over electrons
  for ( unsigned int i=0; i<electrons.size(); ++i ) {
  
    pat::Electron electron = electrons.at(i);

    // keep only e/gamma electrons 
    //if ( electron.isEcalDriven() == false) continue;

    reco::SuperClusterRef scRef = electron.superCluster();
    double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
    double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());

    elePx  -> push_back( electron.trackMomentumAtVtx().x() );
    elePy  -> push_back( electron.trackMomentumAtVtx().y() );
    elePz  -> push_back( electron.trackMomentumAtVtx().z() );
    eleE   -> push_back( scRef->energy());
    eleEt  -> push_back( scRef->energy() * (Rt/R) );
    eleEta -> push_back( scRef->eta());
    elePhi -> push_back( scRef->phi());
    eleCharge -> push_back( electron.gsfTrack()->charge() );
    
    // electron ID
    eleId            -> push_back( electron.electronID(electronID_));
    eleSigmaIEtaIEta -> push_back( electron.scSigmaIEtaIEta());
    eleE1x5          -> push_back(electron.scE1x5());
    eleE2x5          -> push_back(electron.scE2x5Max());
    eleE5x5          -> push_back(electron.scE5x5());

    
    // electron isolation 
    eleTrkIsol    -> push_back(electron.dr03EcalRecHitSumEt()) ;
    eleEcalIsol   -> push_back(electron.dr03HcalDepth1TowerSumEt()) ;
    eleHcalIsolD1 -> push_back(electron.dr03HcalDepth2TowerSumEt()) ;
    eleHcalIsolD2 -> push_back(electron.dr03TkSumPt()) ;


  }// end loop over electron candidates


  //********* MET  
  edm::Handle<edm::View<pat::MET> > metHandle;
  iEvent.getByLabel(metLabel_,metHandle);
  View<pat::MET>  mets = *metHandle;
  for ( unsigned int i=0; i<mets.size(); ++i ) {
    pat::MET met = mets.at(i);
    // corrected MET
    Met    = met.et();
    Mex    = met.px();
    Mey    = met.py();
    MetPhi = met.phi(); 

    // Uncorrected MET
    uncorrMet    = met.uncorrectedPt(pat::MET::uncorrALL);
    uncorrMex    = met.corEx(pat::MET::uncorrALL);
    uncorrMey    = met.corEy(pat::MET::uncorrALL);
    uncorrMetPhi = met.uncorrectedPhi(pat::MET::uncorrALL); 

  }

  //********* JETS
  Handle<View<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetLabel_,jetHandle);
  View<pat::Jet> jets = *jetHandle;
   
  Njets = jets.size();
  for ( unsigned int i=0; i<jets.size(); ++i ) {
    pat::Jet jet = jets.at(i);
    jetPx    -> push_back(jet.px());
    jetPy    -> push_back(jet.py());
    jetPz    -> push_back(jet.pz());
    jetPt    -> push_back(jet.pt());
    jetEta   -> push_back(jet.eta());
    jetPhi   -> push_back(jet.phi());
    jetBdisc -> push_back(jet.bDiscriminator( btagAlgo_));
  }
   

  //**********  MUONS
  Handle<pat::MuonCollection> muons;
  iEvent.getByLabel( muonLabel_ , muons);
  for (pat::MuonCollection::const_iterator muon = muons->begin();  muon != muons->end(); ++muon){
    if(muon->isGlobalMuon()){
      muonPt  -> push_back(muon->pt());
      muonEta -> push_back(muon->eta());
      muonPhi -> push_back(muon->phi());
    }
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
	  if(hlNames_[i]=="HLT_Ele10_SW_L1R" && hltresults->accept(i)) HLTEle10 = 1;	  	  
	  if(hlNames_[i]=="HLT_Ele15_LW_L1R" && hltresults->accept(i)) HLTEle15 = 1;
	  if(hlNames_[i]=="HLT_Ele15_SW_LooseTrackIso_L1R" && hltresults->accept(i)) HLTLooseIsoEle15 = 1;
	  if(hlNames_[i]=="HLT_Ele20_SW_L1R" && hltresults->accept(i)) HLTEle20 = 1;
	  if(hlNames_[i]=="HLT_Photon15_L1R" && hltresults->accept(i)) HLTPhoton15 = 1;
	  if(hlNames_[i]=="HLT_Photon25_L1R" && hltresults->accept(i)) HLTPhoton25 = 1;	  
	}
    }

  //tree fill
  tTreeUtilities_->Fill();
}




// ------------ method called once each job just after ending the event loop  ------------
void 
WprimeAnalyzerPAT::endJob() {
  
  delete MCelePx  ; 
  delete MCelePy  ; 
  delete MCelePz  ; 
  delete MCeleEta ; 
  delete MCelePhi ; 
  delete MCelePid ; 

  delete MCmuonPx ; 
  delete MCmuonPy ; 
  delete MCmuonPz ; 
  delete MCmuonEta ; 
  delete MCmuonPhi ; 
  delete MCmuonPid ; 

  delete isMCmatched ; 
  delete elePx ; 
  delete elePy ; 
  delete elePz ; 
  delete eleE  ;
  delete eleEt ;
  delete eleEta;
  delete elePhi;
  delete eleId ;
  delete eleSigmaIEtaIEta;
  delete eleE1x5;
  delete eleE2x5;
  delete eleE5x5;
  delete eleCharge     ;
  delete eleTrkIsol    ;
  delete eleEcalIsol   ;
  delete eleHcalIsolD1 ;
  delete eleHcalIsolD2 ;

  delete jetPx    ; 
  delete jetPy    ; 
  delete jetPz    ; 
  delete jetPt    ; 
  delete jetEta   ; 
  delete jetPhi   ;  
  delete jetBdisc ; 
  
  delete muonPt  ;
  delete muonEta ;
  delete muonPhi ;
  
  cout << "Analyzed " << EvtCounter <<" events"<< endl;
  


}

