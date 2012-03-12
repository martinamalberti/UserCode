// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// to access recHits and BasicClusters
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

// to use the cluster tools
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalBasicClusterLocalContCorrection.h"



class testEcalClusterCorrections : public edm::EDAnalyzer {
  
 public:
  explicit testEcalClusterCorrections(const edm::ParameterSet &);
   ~testEcalClusterCorrections();
        
   
  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&);

    const CaloGeometry *caloGeometry;

    edm::InputTag barrelClusterCollection_;
    edm::ParameterSet pset;
    EcalBasicClusterLocalContCorrection* localCorr;
};



testEcalClusterCorrections::testEcalClusterCorrections(const edm::ParameterSet& ps)
{
  barrelClusterCollection_ = ps.getParameter<edm::InputTag>("barrelClusterCollection");

  pset  = ps;

 }



testEcalClusterCorrections::~testEcalClusterCorrections()
{
}



void testEcalClusterCorrections::analyze(const edm::Event& ev, const edm::EventSetup& es)
{
 
  localCorr = new EcalBasicClusterLocalContCorrection(pset, es);

 
  edm::Handle< reco::BasicClusterCollection > pEBClusters;
  ev.getByLabel( barrelClusterCollection_, pEBClusters );
  const reco::BasicClusterCollection *ebClusters = pEBClusters.product();
  
   
  
  std::cout << "========== BARREL ==========" << std::endl;
  
  for (reco::BasicClusterCollection::const_iterator it = ebClusters->begin(); it != ebClusters->end(); ++it ) {

    std::cout << "----- new cluster -----" << std::endl;
    std::cout<< "energy  = " << (*it).energy() << std::endl;

    std::cout<< "localContCorr = " << localCorr->getCorrectionVsLocalEta((*it)) <<  std::endl;
    std::cout<< "localContCorr = " << localCorr->getCorrectionVsLocalPhi((*it)) <<  std::endl;
     }
  
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(testEcalClusterCorrections);
