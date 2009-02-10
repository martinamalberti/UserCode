// -*- C++ -*-
//
// Package:    WprimeAnalyzer
// Class:      WprimeAnalyzer
// 
/**\class WprimeAnalyzer WprimeAnalyzer.cc WprimeAnalysis/WprimeENUAnalysis/src/WprimeAnalyzer.cc

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
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeAnalyzer.h"


#define PI 3.141592654
#define TWOPI 6.283185308

using namespace reco;
using namespace std;
using namespace edm;


WprimeAnalyzer::WprimeAnalyzer(const edm::ParameterSet& iConfig)

{
  electronCollection_ = iConfig.getParameter<edm::InputTag>("electrons");
  metCollection_      = iConfig.getParameter<edm::InputTag>("met");
  jetAlgo_            = iConfig.getParameter<edm::InputTag>("jetAlgo");
  hlTriggerResults_   = iConfig.getParameter<edm::InputTag>("TriggerResults");

  electronID_         = iConfig.getUntrackedParameter<std::string>("electronID") ;
  
  TrkIsolationProducer_        = iConfig.getParameter<edm::InputTag>("TrkIsolationProducer");
  EcalIsolationProducer_       = iConfig.getParameter<edm::InputTag>("EcalIsolationProducer");
  HcalIsolationProducerDepth1_ = iConfig.getParameter<edm::InputTag>("HcalIsolationProducerDepth1");
  HcalIsolationProducerDepth2_ = iConfig.getParameter<edm::InputTag>("HcalIsolationProducerDepth2");

//   std::string outputFilename = iConfig.getUntrackedParameter<std::string>("OutputFilename","dummy.root");
//   thefile = new TFile(outputFilename.c_str(),"recreate");
//   thefile->cd();
  
  EvtCounter = 0 ;

 

}


WprimeAnalyzer::~WprimeAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void 
WprimeAnalyzer::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;

  // initialize histograms
  //   hNelectronCandidates =  new TH1F("hNelectronCandidates","hNelectronCandidates",10,0,10);
  //   hTrkIsolation = new TH1F("hTrkIsolation","hTrkIsolation",200,0,100);
  //   hEcalIsolation = new TH1F("hEcalIsolation","hEcalIsolation",400,-100,100);
  //   hHcalIsolationDepth1 = new TH1F("hHcalIsolationDepth1","hHcalIsolationDepth1",400,-100,100);
  //   hHcalIsolationDepth2 = new TH1F("hHcalIsolationDepth2","hHcalIsolationDepth2",400,-100,100);
  //   hEMPlusHcalDepth1 = new TH1F("hEMPlusHcalDepth1","hEMPlusHcalDepth1",400,-100,100);
  
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
WprimeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  Handle<reco::GsfElectronCollection> gsfElectrons;
  iEvent.getByLabel(electronCollection_,gsfElectrons); 

   
  // Read eID results
  edm::Handle<edm::ValueMap<float> > eIDValueMap; 
  iEvent.getByLabel( electronID_ , eIDValueMap ); 
  const edm::ValueMap<float> & eIDmap = * eIDValueMap ;
  
  // Read eIsolation
  edm::Handle<edm::ValueMap<double> > eTrkIsoMap; 
  iEvent.getByLabel( TrkIsolationProducer_, eTrkIsoMap ); 
  
  edm::Handle<edm::ValueMap<double> > eEcalIsoMap; 
  iEvent.getByLabel( EcalIsolationProducer_, eEcalIsoMap ); 
  
  edm::Handle<edm::ValueMap<double> > eHcalDepth1IsoMap; 
  iEvent.getByLabel( HcalIsolationProducerDepth1_, eHcalDepth1IsoMap ); 
  
  edm::Handle<edm::ValueMap<double> > eHcalDepth2IsoMap; 
  iEvent.getByLabel( HcalIsolationProducerDepth2_, eHcalDepth2IsoMap ); 
  
  // Read MC truth association 
  Handle<GenParticleMatch> match;
  iEvent.getByLabel( "matchElectrons", match );
  
  int i = 0;
  NeleCand_ = gsfElectrons->size();
  // Loop over electrons
  for (GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin(); gsfIter!=gsfElectrons->end(); gsfIter++){
    
    edm::Ref<reco::GsfElectronCollection> electronRef(gsfElectrons,i);
    i++;
    edm::RefToBase<reco::Candidate> candRef = edm::RefToBase<reco::Candidate>(electronRef);
   
    GenParticleRef mcMatch = (*match)[candRef];
    cout<<"n ele:"<< i<<"     mc match:"<<mcMatch.isNonnull()<<endl;
    isMCmatched_->push_back(mcMatch.isNonnull());

    if (mcMatch.isNonnull()) {
      MCelePx_  -> push_back( mcMatch->px() );
      MCelePy_  -> push_back( mcMatch->py() );
      MCelePz_  -> push_back( mcMatch->pz() );
      MCeleEta_ -> push_back( mcMatch->eta() );
      MCelePhi_ -> push_back( mcMatch->phi() );
      MCelePid_ -> push_back( mcMatch->pdgId() );
    }
    else {
      MCelePx_  -> push_back( -9999.);
      MCelePy_  -> push_back( -9999.);
      MCelePz_  -> push_back( -9999.);
      MCeleEta_ -> push_back( -9999.);
      MCelePhi_ -> push_back( -9999 );
      MCelePid_ -> push_back( -9999 );
    }

    reco::SuperClusterRef scRef = (*gsfIter).superCluster();
    double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
    double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());

    elePx_ -> push_back( gsfIter->trackMomentumAtVtx().x() );
    elePy_ -> push_back( gsfIter->trackMomentumAtVtx().y() );
    elePz_ -> push_back( gsfIter->trackMomentumAtVtx().z() );
    eleE_  -> push_back( scRef->energy());
    eleEt_ -> push_back( scRef->energy() * (Rt/R) );
    eleEta_-> push_back( scRef->eta());
    elePhi_-> push_back( scRef->phi());
    eleId_ ->push_back( eIDmap[electronRef] );
    eleCharge_ ->push_back( (*gsfIter).gsfTrack()->charge() );
    eleTrkIsol_ ->push_back((*eTrkIsoMap)[ electronRef ]);
    eleEcalIsol_ ->push_back((*eEcalIsoMap)[ electronRef ]);
    eleHcalIsolD1_ ->push_back((*eHcalDepth1IsoMap)[ electronRef ]);
    eleHcalIsolD2_ ->push_back((*eHcalDepth2IsoMap)[ electronRef ]);

  }// end loop over electron candidates

 

  
  //********* RECO MET  from Calo Towers
  Handle<CaloMETCollection> metCollection;
  iEvent.getByLabel(metCollection_, metCollection) ;
  CaloMETCollection::const_iterator caloMET = metCollection->begin();
  met_ = caloMET->pt();
  mex_ = caloMET->px();
  mey_ = caloMET->py();
  metPhi_ = caloMET->phi(); 
  
  //********* RECO JETS 
  Handle<reco::CaloJetCollection> jets;
  iEvent.getByLabel(jetAlgo_,jets);
  reco::CaloJetCollection::const_iterator jet;
  
  Njets_ = jets->size();
  for( jet = jets->begin(); jet != jets->end(); ++ jet ) {
    jetPx_->push_back(jet->px());
    jetPy_->push_back(jet->py());
    jetPz_->push_back(jet->pz());
    jetPt_->push_back(jet->pt());
    jetEta_->push_back(jet->eta());
    jetPhi_->push_back(jet->phi());
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
WprimeAnalyzer::endJob() {
  
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

