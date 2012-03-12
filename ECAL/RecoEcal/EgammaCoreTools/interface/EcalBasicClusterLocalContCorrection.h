#ifndef RecoEcal_EgammaCoreTools_EcalBasicClusterLocalContCorrection_h
#define RecoEcal_EgammaCoreTools_EcalBasicClusterLocalContCorrection_h

/** \class EcalBasicClusterLocalContCorrection
  *  Function to correct basic cluster energy for local containment
  *
  *  $Id: EcalBasicClusterLocalContCorrection.h
  *  $Date:
  *  $Revision:
  *  \author Martina Malberti, February 2012
  */

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "TF1.h"

class CaloGeometry;

class EcalBasicClusterLocalContCorrection  {

 public:

  EcalBasicClusterLocalContCorrection(const edm::ParameterSet &, const edm::EventSetup &);
  ~EcalBasicClusterLocalContCorrection();

  // -- to get the corrections
  float getCorrectionVsLocalEta(const reco::CaloCluster & ) ;
  float getCorrectionVsLocalPhi(const reco::CaloCluster & ) ;
  
 private:

  void setCorrectionParameters(const edm::ParameterSet &ps);
  void getGeometry( const edm::EventSetup &es);
  int  getEcalModule( DetId id );
  std::pair<double,double> getLocalPosition(const reco::CaloCluster & caloCluster) ;
  
  // -- geometry
  const CaloGeometry *caloGeometry;

  // -- correction parameters
  std::vector<double> paramsEta0_;
  std::vector<double> paramsEta1_;
  std::vector<double> paramsEta2_;
  std::vector<double> paramsEta3_;

  std::vector<double> paramsPhi0_;
  std::vector<double> paramsPhi1_;
  std::vector<double> paramsPhi2_;
  std::vector<double> paramsPhi3_;
  
  // -- functions
  TF1 *fEta[4];
  TF1 *fPhi[4];


};

#endif
