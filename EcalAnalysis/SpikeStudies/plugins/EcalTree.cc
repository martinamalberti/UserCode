// -*- C++ -*-
//
// Package:   EcalTree
// Class:     EcalTree
//
 
#include "EcalAnalysis/SpikeStudies/plugins/EcalTree.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"

#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TSystem.h"
#include <memory>
#include <vector>
#include <iostream>
#include <iterator>

using namespace cms ;
using namespace edm ;
using namespace std ;
using namespace reco;

EcalTree::EcalTree (const edm::ParameterSet& iConfig)
{
  ebRecHitCollection_  = iConfig.getParameter<edm::InputTag> ("ebRecHitCollection");
  ebClusterCollection_ = iConfig.getParameter<edm::InputTag> ("ebClusterCollection");
  ebDigiCollection_    = iConfig.getParameter<edm::InputTag> ("ebDigiCollection");
  L1InputTag_          = iConfig.getParameter<edm::InputTag> ("L1InputTag");
 
  //radiusForIso_        = iConfig.getParameter<double> ("radiusForIso");
  energyCutForIso_     = iConfig.getUntrackedParameter<double> ("energyCutForIso",0.);
  
  naiveId_ = 0;

}


// -----------------------------------------------------------------------------------------


EcalTree::~EcalTree ()
{
}


// -----------------------------------------------------------------------------------------


void EcalTree::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  naiveId_++;

  // calo topology
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  const CaloTopology *topology = pTopology.product();

  // calo geometry
  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);
  const CaloGeometry *geometry = pGeometry.product();

  // Algo and Technical L1 bits
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* L1Menu = menuRcd.product();
  
  edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
  iEvent.getByLabel(L1InputTag_, gtRecord);
  
  // Ecal Basic Cluster Collection
  edm::Handle<BasicClusterCollection > pEBClusters;
  iEvent.getByLabel( ebClusterCollection_, pEBClusters );
  const reco::BasicClusterCollection* ebClusters = pEBClusters.product();
  
  // Ecal barrel RecHits 
  edm::Handle<EcalRecHitCollection> pBarrelEcalRecHits ;
  iEvent.getByLabel (ebRecHitCollection_, pBarrelEcalRecHits) ;
  const EcalRecHitCollection* theBarrelEcalRecHits = pBarrelEcalRecHits.product () ;   
  
  if (! (pBarrelEcalRecHits.isValid ()) )
    {
      LogWarning ("AnomalousChannelsAnalyzer") << ebRecHitCollection_ 
					       << " not available" ;
      return ;
    }
  
  //ECAL DIGIS
  edm::Handle<EBDigiCollection> ebDigis;
  iEvent.getByLabel (ebDigiCollection_, ebDigis) ;
  const EBDigiCollection* theEcalBarrelDigis = ebDigis.product () ;
  
  if (! (ebDigis.isValid ()) )
    {
      LogWarning ("EcalTree") << ebDigiCollection_
			      << " not available" ;
      return ;
    }

  //Fill Tree
  initializeBranches(tree_, myTreeVariables_);

  myTreeVariables_.lumiId       = iEvent.luminosityBlock();
  myTreeVariables_.BX           = iEvent.bunchCrossing();
  myTreeVariables_.runId        = iEvent.id ().run () ;
  myTreeVariables_.eventId      = iEvent.id ().event () ;
  myTreeVariables_.eventNaiveId = naiveId_ ;

  dumpL1Info(gtRecord, myTreeVariables_) ;
  dumpBarrelInfo(topology, geometry, theEcalBarrelDigis, theBarrelEcalRecHits, ebClusters, myTreeVariables_) ;

  tree_ -> Fill();
}

// -----------------------------------------------------------------------------------------

void EcalTree::endJob ()
{
  cout<< "Analyzed " <<  naiveId_ << " events" <<endl;
}

// ----------------------------------------------------------------------------------------

void EcalTree::beginJob()
{
  //file to save output
  edm::Service<TFileService> fs;
  // Initialize Tree
  tree_ = fs->make<TTree>("EcalAnalysisTree","EcalAnalysisTree");
  setBranches (tree_, myTreeVariables_) ;
}






// -----------------------------------------------------------------------------------------

void EcalTree::dumpBarrelInfo (	const CaloTopology * topology,
				const CaloGeometry * geometry,
				const EBDigiCollection* theEcalBarrelDigis,
				const EcalRecHitCollection* theBarrelEcalRecHits,
				const reco::BasicClusterCollection *ebClusters,
				EcalTreeContent & myTreeVariables_)
{
  
  EcalRecHitMetaCollection mhits(*theBarrelEcalRecHits);

  // cerco il BC corrispondente e trovo R9 = E1/E9
  float S1oS9=-9999.;
  float S4oS1=-9999.;

  for (reco::BasicClusterCollection::const_iterator it = ebClusters->begin(); it != ebClusters->end(); ++it ) 
    {

      //solo barrel: seed
      EBDetId ebid = EcalClusterTools::getMaximum( (*it).hitsAndFractions(), theBarrelEcalRecHits).first;
      EcalRecHitCollection::const_iterator thishitEB_seed = theBarrelEcalRecHits->find (ebid) ;
      myTreeVariables_.ecalRecHitMatrix[ myTreeVariables_.nEcalRecHits ][2][2] = thishitEB_seed -> energy();;
      myTreeVariables_.ecalRecHitMatrixFlag[ myTreeVariables_.nEcalRecHits ][2][2] = thishitEB_seed -> recoFlag();;
  
      for(int xx = 0; xx < 5; ++xx)
	for(int yy = 0; yy < 5; ++yy)
	  if(xx != 2 && yy != 2)
	    {
              myTreeVariables_.ecalRecHitMatrix[ myTreeVariables_.nEcalRecHits ][xx][yy] = 0.;
              myTreeVariables_.ecalRecHitMatrixFlag[ myTreeVariables_.nEcalRecHits ][xx][yy] = 0.;

	      std::vector<DetId> vector =  EcalClusterTools::matrixDetId(topology, ebid, xx-2, xx-2, yy-2, yy-2);
	      if(vector.size() == 0) continue;
	      EcalRecHitCollection::const_iterator iterator = theBarrelEcalRecHits->find (vector.at(0)) ;
	      if(iterator == theBarrelEcalRecHits->end()) continue;
	      myTreeVariables_.ecalRecHitMatrix[ myTreeVariables_.nEcalRecHits ][xx][yy] = iterator -> energy();
	      myTreeVariables_.ecalRecHitMatrixFlag[ myTreeVariables_.nEcalRecHits ][xx][yy] = iterator -> recoFlag();
	    }
      

//       std::vector<DetId> matrixDetId =  EcalClusterTools::matrixDetId(topology, ebid, -2, 2, -2, 2);
//       for (std::vector<DetId>::const_iterator matrixItr = matrixDetId.begin(); 
// 	   matrixItr != matrixDetId.end();
// 	   ++matrixItr)
// 	{
// 	  EBDetId dummy(*matrixItr);
// 	  cout << "ieta = " << dummy.ieta() << std::endl;
// 	  cout << "iphi = " << dummy.iphi() << std::endl;


// 	}
      
      S1oS9 = EcalClusterTools::eMax( *it, theBarrelEcalRecHits)/
	EcalClusterTools::e3x3( *it, theBarrelEcalRecHits, topology );
      
      float energyNeighbours = 
	EcalClusterTools::eTop( *it, theBarrelEcalRecHits, topology)+
	EcalClusterTools::eRight( *it, theBarrelEcalRecHits, topology)+
	EcalClusterTools::eBottom( *it, theBarrelEcalRecHits, topology)+
	EcalClusterTools::eLeft( *it, theBarrelEcalRecHits, topology);
      
      S4oS1 = energyNeighbours/EcalClusterTools::eMax( *it, theBarrelEcalRecHits);

      
      myTreeVariables_.ecalRecHitType     [ myTreeVariables_.nEcalRecHits ] = 0;
      myTreeVariables_.ecalRecHitEnergy   [ myTreeVariables_.nEcalRecHits ] = thishitEB_seed -> energy();
      myTreeVariables_.ecalRecHitOutOfTimeEnergy[ myTreeVariables_.nEcalRecHits ] = thishitEB_seed -> outOfTimeEnergy();
      myTreeVariables_.ecalRecHitIEta     [ myTreeVariables_.nEcalRecHits ] = ebid.ieta();
      myTreeVariables_.ecalRecHitIPhi     [ myTreeVariables_.nEcalRecHits ] = ebid.iphi();
      myTreeVariables_.ecalRecHitTime     [ myTreeVariables_.nEcalRecHits ] = thishitEB_seed -> time();
      myTreeVariables_.ecalRecHitChi2     [ myTreeVariables_.nEcalRecHits ] = thishitEB_seed->chi2() ;
      myTreeVariables_.ecalRecHitOutOfTimeChi2 [ myTreeVariables_.nEcalRecHits ] = thishitEB_seed->outOfTimeChi2();
      myTreeVariables_.ecalRecHitRawId    [ myTreeVariables_.nEcalRecHits ] = ebid.rawId();
      myTreeVariables_.ecalRecHitRecoFlag [ myTreeVariables_.nEcalRecHits ] = thishitEB_seed -> recoFlag();
      myTreeVariables_.ecalRecHitR9       [ myTreeVariables_.nEcalRecHits ] = S1oS9;
      myTreeVariables_.ecalRecHitS4oS1    [ myTreeVariables_.nEcalRecHits ] = S4oS1;
      
      
      // ecal activity in R = 0.3
      CaloConeSelector *sel03 = new CaloConeSelector(0.3, geometry , DetId::Ecal, EcalBarrel);
      
      std::auto_ptr<CaloRecHitMetaCollectionV> chosen1r03 = sel03->select(ebid.ieta(), ebid.iphi(), mhits);
      std::auto_ptr<CaloRecHitMetaCollectionV> chosen2r03 = sel03->select(-ebid.ieta(), -ebid.iphi(), mhits);
      
      float myIso03 = 0;
      for (CaloRecHitMetaCollectionV::const_iterator recIt = chosen1r03->begin(); 
	   recIt!= chosen1r03->end () ; ++recIt) {
	if ( recIt->energy() < energyCutForIso_) continue;  //dont fill if below E noise value
	if ( EBDetId(recIt->detid()).ieta() ==  ebid.ieta() && 
	     EBDetId(recIt->detid()).iphi() ==  ebid.iphi() ) continue;
	GlobalPoint mycell = geometry->getPosition(recIt->detid());
	myIso03 +=  recIt->energy()*sin(2*atan(exp(mycell.eta())));
      }
      
      myTreeVariables_.ecalRecHitIso03 [ myTreeVariables_.nEcalRecHits ][0] = myIso03  ;
      
  
      myIso03 = 0;
      for (CaloRecHitMetaCollectionV::const_iterator recIt = chosen2r03->begin(); 
	   recIt!= chosen2r03->end () ; ++recIt) {
	if ( recIt->energy() < energyCutForIso_) continue;  //dont fill if below E noise value	  
	if ( EBDetId(recIt->detid()).ieta() == -ebid.ieta() && 
	     EBDetId(recIt->detid()).iphi() == -ebid.iphi() ) continue;
	GlobalPoint mycell = geometry->getPosition(recIt->detid());
	myIso03 +=  recIt->energy()*sin(2*atan(exp(mycell.eta())));
      }
      
      myTreeVariables_.ecalRecHitIso03 [ myTreeVariables_.nEcalRecHits ][1] = myIso03 ;
  
      
      delete  sel03;
      
      
      
      // ecal activity in R = 0.4
      CaloConeSelector *sel04 = new CaloConeSelector(0.4, geometry , DetId::Ecal, EcalBarrel);
      
      std::auto_ptr<CaloRecHitMetaCollectionV> chosen1r04 = sel04->select(ebid.ieta(), ebid.iphi(), mhits);
      std::auto_ptr<CaloRecHitMetaCollectionV> chosen2r04 = sel04->select(-ebid.ieta(), -ebid.iphi(), mhits);
      
      float myIso04 = 0;
      for (CaloRecHitMetaCollectionV::const_iterator recIt = chosen1r04->begin(); 
	   recIt!= chosen1r04->end () ; ++recIt) {
	if ( recIt->energy() < energyCutForIso_) continue;  //dont fill if below E noise value
	if ( EBDetId(recIt->detid()).ieta() ==  ebid.ieta() && 
	     EBDetId(recIt->detid()).iphi() ==  ebid.iphi() ) continue;
	GlobalPoint mycell = geometry->getPosition(recIt->detid());
	myIso04 +=  recIt->energy()*sin(2*atan(exp(mycell.eta())));
      }
      
      myTreeVariables_.ecalRecHitIso04 [ myTreeVariables_.nEcalRecHits ][0] = myIso04  ;
      
      
      myIso04 = 0;
      for (CaloRecHitMetaCollectionV::const_iterator recIt = chosen2r04->begin(); 
	   recIt!= chosen2r04->end () ; ++recIt) {
	if ( recIt->energy() < energyCutForIso_) continue;  //dont fill if below E noise value	  
	if ( EBDetId(recIt->detid()).ieta() == -ebid.ieta() && 
	     EBDetId(recIt->detid()).iphi() == -ebid.iphi() ) continue;
	GlobalPoint mycell = geometry->getPosition(recIt->detid());
	myIso04 +=  recIt->energy()*sin(2*atan(exp(mycell.eta())));
      }
      
      myTreeVariables_.ecalRecHitIso04 [ myTreeVariables_.nEcalRecHits ][1] = myIso04  ;
      
      delete  sel04;
      
      
      
      
      
      
      // DIGIS
      EBDigiCollection::const_iterator digiItr = theEcalBarrelDigis->begin();
      while(digiItr != theEcalBarrelDigis->end() && ((*digiItr).id() != ebid))
	{
	  ++digiItr;
	}
      EcalDataFrame df = *digiItr;
      
      for (int i=0; i < df.size(); i++ ) {
	myTreeVariables_.ecalDigis       [ myTreeVariables_.nEcalRecHits ][i] = df.sample(i).adc();
	myTreeVariables_.ecalGainId      [ myTreeVariables_.nEcalRecHits ][i] = df.sample(i).gainId();
	
	//cout <<  i << "  " <<  myTreeVariables_.ecalDigis       [ myTreeVariables_.nEcalRecHits ][i] << endl;
	
      }
      
      
      ++myTreeVariables_.nEcalRecHits;
      
      return ;
      
    } // dumpBarrelInfo  
}
// -----------------------------------------------------------------------------------------


// dump L1Info
void EcalTree::dumpL1Info ( edm::Handle<L1GlobalTriggerReadoutRecord>  gtRecord ,
			    EcalTreeContent & myTreeVariables_)
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


