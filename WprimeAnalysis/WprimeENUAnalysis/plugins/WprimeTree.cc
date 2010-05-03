// -*- C++ -*-
//
// Package:   WprimeTree
// Class:     WprimeTree
//
 
#include "WprimeAnalysis/WprimeENUAnalysis/plugins/WprimeTree.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Flags.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TSystem.h"
#include <memory>
#include <vector>
#include <iostream>
#include <iterator>

#define PI 3.141592654
#define TWOPI 6.283185308

using namespace cms ;
using namespace edm ;
using namespace std ;
using namespace reco;

WprimeTree::WprimeTree (const edm::ParameterSet& iConfig)
{
  eleLabel_      = iConfig.getParameter<edm::InputTag>("electronTag");
  metLabel_      = iConfig.getParameter<edm::InputTag>("metTag");
  jetLabel_      = iConfig.getParameter<edm::InputTag>("jetTag");
  muonLabel_     = iConfig.getParameter<edm::InputTag>("muonTag");

  HLTInputTag_   = iConfig.getParameter<edm::InputTag>("HLTInputTag");
  L1InputTag_    = iConfig.getParameter<edm::InputTag>("L1InputTag");

  electronID_    = iConfig.getUntrackedParameter<std::string>("electronID") ;
  btagAlgo_      = iConfig.getUntrackedParameter<std::string>("btagAlgo") ;

  naiveId_ = 0;


}


// -----------------------------------------------------------------------------------------


WprimeTree::~WprimeTree ()
{
}


// -----------------------------------------------------------------------------------------


void WprimeTree::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  naiveId_++;

  //*********** Algo and Technical L1 bits
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
   
  edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
  iEvent.getByLabel(L1InputTag_, gtRecord);

  //*********** HLT INFO
  edm::TriggerNames triggerNames_;
  Handle<TriggerResults> hltresults;
  iEvent.getByLabel(HLTInputTag_,hltresults);


  //************* ELECTRONS
  Handle<View<pat::Electron> > electronHandle;
  iEvent.getByLabel(eleLabel_,electronHandle);
  View<pat::Electron> electrons = *electronHandle;


  //************* MET 
  edm::Handle<edm::View<pat::MET> > metHandle;
  iEvent.getByLabel(metLabel_,metHandle);
  View<pat::MET>  mets = *metHandle;

  //************* JETS
  Handle<View<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetLabel_,jetHandle);
  View<pat::Jet> jets = *jetHandle;


  //************* MUONS
  // Handle<pat::MuonCollection> muons;
  //   iEvent.getByLabel( muonLabel_ , muons);
  Handle<View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muonLabel_,muonHandle);
  View<pat::Muon> muons = *muonHandle;
  
  
  //Fill Tree
  initializeBranches(tree_, myTreeVariables_);

  myTreeVariables_.lumiId       = iEvent.luminosityBlock();
  myTreeVariables_.BX           = iEvent.bunchCrossing();
  myTreeVariables_.runId        = iEvent.id ().run () ;
  myTreeVariables_.eventId      = iEvent.id ().event () ;
  myTreeVariables_.eventNaiveId = naiveId_ ;

  dumpL1Info(gtRecord, myTreeVariables_) ;
  dumpHLTInfo(hltresults, myTreeVariables_) ;
  dumpElectronInfo( electrons, myTreeVariables_) ;
  dumpMetInfo     (      mets, myTreeVariables_) ;
  dumpJetInfo     (      jets, myTreeVariables_) ;
  dumpMuonInfo    (     muons, myTreeVariables_) ;

  tree_ -> Fill();
}






// -----------------------------------------------------------------------------------------

void WprimeTree::endJob ()
{
  cout<< "Analyzed " <<  naiveId_ << " events" <<endl;
}

// ----------------------------------------------------------------------------------------

void WprimeTree::beginJob()
{
  //file to save output
  edm::Service<TFileService> fs;
  // Initialize Tree
  tree_ = fs->make<TTree>("WprimeAnalysisTree","WprimeAnalysisTree");
  setBranches (tree_, myTreeVariables_) ;
}






// -----------------------------------------------------------------------------------------

void WprimeTree::dumpElectronInfo ( View<pat::Electron> electrons,
				    WprimeTreeContent & myTreeVariables_)
{
  
  // Loop over electrons
  for ( unsigned int i=0; i<electrons.size(); ++i ) {
    
    pat::Electron electron = electrons.at(i);
    
    // keep only ecal driven electrons
    if (!electron.ecalDrivenSeed()) continue;

    reco::SuperClusterRef scRef = electron.superCluster();
    double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
    double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());

    myTreeVariables_.elePx[ myTreeVariables_.nElectrons ]  = electron.trackMomentumAtVtx().x() ;
    myTreeVariables_.elePy[ myTreeVariables_.nElectrons ]  = electron.trackMomentumAtVtx().y() ;
    myTreeVariables_.elePz[ myTreeVariables_.nElectrons ]  = electron.trackMomentumAtVtx().z() ;
    myTreeVariables_.eleE [ myTreeVariables_.nElectrons ]  = scRef->energy() ;
    myTreeVariables_.eleEt[ myTreeVariables_.nElectrons ]  = scRef->energy()*(Rt/R) ;
    myTreeVariables_.eleEta[ myTreeVariables_.nElectrons ] = scRef->eta();
    myTreeVariables_.elePhi[ myTreeVariables_.nElectrons ] = scRef->phi() ;

    myTreeVariables_.eleCharge[ myTreeVariables_.nElectrons ] = electron.gsfTrack()->charge();

    myTreeVariables_.eleId[ myTreeVariables_.nElectrons ] = int(electron.electronID(electronID_));

    
    myTreeVariables_.eleSigmaIEtaIEta[ myTreeVariables_.nElectrons ] =  electron.sigmaIetaIeta();
    myTreeVariables_.eleE1x5[ myTreeVariables_.nElectrons ] = electron.e1x5();
    myTreeVariables_.eleE2x5[ myTreeVariables_.nElectrons ] = electron.e2x5Max();
    myTreeVariables_.eleE5x5[ myTreeVariables_.nElectrons ] = electron.e5x5();

    myTreeVariables_.eleTrkIso   [ myTreeVariables_.nElectrons ] = electron.dr03TkSumPt();
    myTreeVariables_.eleEcalIso  [ myTreeVariables_.nElectrons ] = electron.dr03EcalRecHitSumEt();
    myTreeVariables_.eleHcalIsoD1[ myTreeVariables_.nElectrons ] = electron.dr03HcalDepth1TowerSumEt();
    myTreeVariables_.eleHcalIsoD2[ myTreeVariables_.nElectrons ] = electron.dr03HcalDepth2TowerSumEt();

    ++myTreeVariables_.nElectrons;

  }// end loop over electron candidates
  
  return ;
  
} // dumpElectronInfo  


// -----------------------------------------------------------------------------------------

// dump MET Info
void WprimeTree::dumpMetInfo ( View<pat::MET>  mets ,
			       WprimeTreeContent & myTreeVariables_)

{

  for ( unsigned int i=0; i<mets.size(); ++i ) {
    pat::MET met = mets.at(i);
    // corrected MET
    myTreeVariables_.Met = met.et();
    myTreeVariables_.Mex = met.px();
    myTreeVariables_.Mey = met.py();
    myTreeVariables_.MetPhi = met.phi();

    // Uncorrected MET
    myTreeVariables_.uncorrMet = met.uncorrectedPt(pat::MET::uncorrALL);
    myTreeVariables_.uncorrMex = met.corEx(pat::MET::uncorrALL);
    myTreeVariables_.uncorrMey = met.corEy(pat::MET::uncorrALL);
    myTreeVariables_.uncorrMetPhi = met.uncorrectedPhi(pat::MET::uncorrALL);
  }

  return ;
  
} // dumpMetInfo 


// -----------------------------------------------------------------------------------------

// dump JET Info
void WprimeTree::dumpJetInfo ( View<pat::Jet>  jets ,
			       WprimeTreeContent & myTreeVariables_)

{
  for ( unsigned int i=0; i<jets.size(); ++i ) {
    pat::Jet jet = jets.at(i);
    myTreeVariables_.jetPx[ myTreeVariables_.nJets ]    = jet.px();
    myTreeVariables_.jetPy[ myTreeVariables_.nJets ]    = jet.py();
    myTreeVariables_.jetPz[ myTreeVariables_.nJets ]    = jet.pz();
    myTreeVariables_.jetPt[ myTreeVariables_.nJets ]    = jet.pt();
    myTreeVariables_.jetEta[ myTreeVariables_.nJets ]   = jet.eta();
    myTreeVariables_.jetPhi[ myTreeVariables_.nJets ]   = jet.phi();
    myTreeVariables_.jetBdisc[ myTreeVariables_.nJets ] =jet.bDiscriminator( btagAlgo_);
  
    ++myTreeVariables_.nJets;
  }

  return ;
  
} // dumpJetInfo 

// -----------------------------------------------------------------------------------------

void WprimeTree::dumpMuonInfo ( View<pat::Muon>  muons ,
			       WprimeTreeContent & myTreeVariables_)

{
  for ( unsigned int i=0; i<muons.size(); ++i ) {
    pat::Muon muon = muons.at(i);
    myTreeVariables_.muonPx[ myTreeVariables_.nMuons ]    = muon.px();
    myTreeVariables_.muonPy[ myTreeVariables_.nMuons ]    = muon.py();
    myTreeVariables_.muonPz[ myTreeVariables_.nMuons ]    = muon.pz();
    myTreeVariables_.muonPt[ myTreeVariables_.nMuons ]    = muon.pt();
    myTreeVariables_.muonEta[ myTreeVariables_.nMuons ]   = muon.eta();
    myTreeVariables_.muonPhi[ myTreeVariables_.nMuons ]   = muon.phi();
      
    ++myTreeVariables_.nMuons;
  }

  return ;
  
} // dumpMuonInfo 


// -----------------------------------------------------------------------------------------

// dump L1Info
void WprimeTree::dumpL1Info ( edm::Handle<L1GlobalTriggerReadoutRecord>  gtRecord ,
			    WprimeTreeContent & myTreeVariables_)
{

  DecisionWord AlgoWord = gtRecord->decisionWord();
  TechnicalTriggerWord TechWord = gtRecord->technicalTriggerWord();
 
  // Loop over the technical bits
  for (unsigned int ibit = 0; ibit < TechWord.size(); ibit++) 
    {
      myTreeVariables_.techL1Bit[ibit] = TechWord[ibit];
    }
  
  // Loop over the algo bits
  for (unsigned int ibit = 0; ibit < AlgoWord.size(); ibit++) 
    {
      myTreeVariables_.algoL1Bit[ibit] = AlgoWord[ibit];
    }

  return ;

}// dumpL1Info  


// -----------------------------------------------------------------------------------------

// dump HLTInfo
void WprimeTree::dumpHLTInfo ( Handle<TriggerResults>  hltresults ,
			       WprimeTreeContent & myTreeVariables_)
{
  edm::TriggerNames triggerNames_;
  std::vector<std::string> hlNames_;

  //if (&hltresults)
  //  {
      int ntrigs=(*hltresults).size();

      if (ntrigs==0){std::cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;}
      triggerNames_.init(*hltresults);
      hlNames_=triggerNames_.triggerNames();

      const unsigned int n(hlNames_.size());
      for (unsigned int i = 0; i < n; ++i)
        {
	  if(hlNames_[i]=="HLT_Ele15_SW_L1R" && hltresults->accept(i)) myTreeVariables_.HLT_Ele15_SW_L1R = 1;
	  if(hlNames_[i]=="HLT_Ele15_LooseTrackIso_L1R" && hltresults->accept(i)) myTreeVariables_.HLT_Ele15_LooseTrackIso_L1R = 1;
	  if(hlNames_[i]=="HLT_Photon15_L1R" && hltresults->accept(i)) myTreeVariables_.HLT_Photon15_L1R = 1;
	  if(hlNames_[i]=="HLT_Photon25_L1R" && hltresults->accept(i)) myTreeVariables_.HLT_Photon25_L1R = 1;
        }
      //}




  return ;

}// dumpL1Info  


