#include "RecoEcal/EgammaCoreTools/interface/EcalBasicClusterLocalContCorrection.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"

#include "TVector2.h"
#include "TMath.h"


#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <iostream>


//-----------------------------------------------------------------------------------------------------
EcalBasicClusterLocalContCorrection::EcalBasicClusterLocalContCorrection( const edm::ParameterSet & ps , const edm::EventSetup &es)
{  
  setCorrectionParameters(ps);

  fEta[0] = new TF1("fEta0","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);
  fEta[1] = new TF1("fEta1","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);
  fEta[2] = new TF1("fEta2","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);
  fEta[3] = new TF1("fEta3","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);

  fEta[0]->SetParameters(paramsEta0_[0],paramsEta0_[1],paramsEta0_[2]);
  fEta[1]->SetParameters(paramsEta1_[0],paramsEta1_[1],paramsEta1_[2]);
  fEta[2]->SetParameters(paramsEta2_[0],paramsEta2_[1],paramsEta2_[2]);
  fEta[3]->SetParameters(paramsEta3_[0],paramsEta3_[1],paramsEta3_[2]);

  fPhi[0] = new TF1("fPhi0","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);
  fPhi[1] = new TF1("fPhi1","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);
  fPhi[2] = new TF1("fPhi2","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);
  fPhi[3] = new TF1("fPhi3","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);

  fPhi[0]->SetParameters(paramsPhi0_[0],paramsPhi0_[1],paramsPhi0_[2]);
  fPhi[1]->SetParameters(paramsPhi1_[0],paramsPhi1_[1],paramsPhi1_[2]);
  fPhi[2]->SetParameters(paramsPhi2_[0],paramsPhi2_[1],paramsPhi2_[2]);
  fPhi[3]->SetParameters(paramsPhi3_[0],paramsPhi3_[1],paramsPhi3_[2]);

  
  getGeometry( es );


}
//------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------
void EcalBasicClusterLocalContCorrection::setCorrectionParameters(const edm::ParameterSet& ps) 
{

  // --- load parameters for corrections from config file
 
  const edm::ParameterSet params_ = ps.getParameter<edm::ParameterSet>("localCorrectionParameters") ;

  paramsEta0_   = params_.getParameter< std::vector<double> >("paramsEta0");
  paramsEta1_   = params_.getParameter< std::vector<double> >("paramsEta1");
  paramsEta2_   = params_.getParameter< std::vector<double> >("paramsEta2");
  paramsEta3_   = params_.getParameter< std::vector<double> >("paramsEta3");
  
  paramsPhi0_   = params_.getParameter< std::vector<double> >("paramsPhi0");
  paramsPhi1_   = params_.getParameter< std::vector<double> >("paramsPhi1");
  paramsPhi2_   = params_.getParameter< std::vector<double> >("paramsPhi2");
  paramsPhi3_   = params_.getParameter< std::vector<double> >("paramsPhi3");
  
}
//------------------------------------------------------------------------------------------------------


void EcalBasicClusterLocalContCorrection::getGeometry( const edm::EventSetup &es )
{
  edm::ESHandle<CaloGeometry> theCaloGeom;
  es.get<CaloGeometryRecord>().get(theCaloGeom);
  caloGeometry = theCaloGeom.product();

}



//------------------------------------------------------------------------------------------------------
float EcalBasicClusterLocalContCorrection::getCorrectionVsLocalEta(const reco::CaloCluster & basicCluster ) 
{

  //--- local cluster coordinates 
  std::pair<double,double> localPosition = getLocalPosition(basicCluster);
  float localEta = localPosition.first;

  //--- ecal module
  int imod    = getEcalModule(basicCluster.seed());

  //--- correction
  float corr = fEta[imod]-> Eval(localEta);

  return(1./corr);
 
}
//------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------
float EcalBasicClusterLocalContCorrection::getCorrectionVsLocalPhi(const reco::CaloCluster & basicCluster ) 
{

  //--- local cluster coordinates 
  std::pair<double,double> localPosition = getLocalPosition(basicCluster);
  float localPhi = localPosition.second;

  //--- ecal module
  int imod    = getEcalModule(basicCluster.seed());

  //--- correction
  float corr = fPhi[imod]-> Eval(localPhi);

  return(1./corr);
 
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
int EcalBasicClusterLocalContCorrection::getEcalModule( DetId id )
{
  int mod = 0;
  int ieta = (EBDetId(id)).ieta();
  
  if (fabs(ieta) <=25 ) mod = 0;
  if (fabs(ieta) <=45 ) mod = 1;
  if (fabs(ieta) <=65 ) mod = 2;
  if (fabs(ieta) <=85 ) mod = 3;

  return (mod);
 
}
//------------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------------
std::pair<double,double> EcalBasicClusterLocalContCorrection::getLocalPosition(const reco::CaloCluster & caloCluster)
{


  //--------------if barrel calculate local position wrt xtal center -------------------
  const CaloSubdetectorGeometry* geom = caloGeometry->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);//EcalBarrel = 1
  

  const math::XYZPoint position_ = caloCluster.position(); 
  double Theta = -position_.theta()+0.5*TMath::Pi();
  double Eta = position_.eta();
  double Phi = TVector2::Phi_mpi_pi(position_.phi());
  
  
  
  //Calculate expected depth of the maximum shower from energy (like in PositionCalc::Calculate_Location()):
  // The parameters X0 and T0 are hardcoded here because these values were used to calculate the corrections:
  const float X0 = 0.89; const float T0 = 7.4;
  double depth = X0 * (T0 + log(caloCluster.energy()));
  
  
  //search which crystal is closest to the cluster position and call it crystalseed:
  //std::vector<DetId> crystals_vector = *scRef.getHitsByDetId();   //deprecated
  std::vector< std::pair<DetId, float> > crystals_vector = caloCluster.hitsAndFractions();
  float dphimin=999.;
  float detamin=999.;
  int ietaclosest = 0;
  int iphiclosest = 0;
  
  
  for (unsigned int icry=0; icry!=crystals_vector.size(); ++icry) 
    {    
      
      EBDetId crystal(crystals_vector[icry].first);
      const CaloCellGeometry* cell=geom->getGeometry(crystal);// problema qui
      GlobalPoint center_pos = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);
      double EtaCentr = center_pos.eta();
      double PhiCentr = TVector2::Phi_mpi_pi(center_pos.phi());
      if (TMath::Abs(EtaCentr-Eta) < detamin) {
	detamin = TMath::Abs(EtaCentr-Eta); 
	ietaclosest = crystal.ieta();
      }
      if (TMath::Abs(TVector2::Phi_mpi_pi(PhiCentr-Phi)) < dphimin) {
	dphimin = TMath::Abs(TVector2::Phi_mpi_pi(PhiCentr-Phi)); 
	iphiclosest = crystal.iphi();
      }
      
    }
 
  
  EBDetId crystalseed(ietaclosest, iphiclosest);

  
  // Get center cell position from shower depth
  const CaloCellGeometry* cell=geom->getGeometry(crystalseed);
  GlobalPoint center_pos = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);
  
  //PHI
  double PhiCentr = TVector2::Phi_mpi_pi(center_pos.phi());
  double PhiWidth = (TMath::Pi()/180.);
  double PhiCry = (TVector2::Phi_mpi_pi(Phi-PhiCentr))/PhiWidth;
  if (PhiCry>0.5) PhiCry=0.5;
  if (PhiCry<-0.5) PhiCry=-0.5;
  //flip to take into account ECAL barrel symmetries:
  if (ietaclosest<0) PhiCry *= -1.;
  

  //ETA
  double ThetaCentr = -center_pos.theta()+0.5*TMath::Pi();
  double ThetaWidth = (TMath::Pi()/180.)*TMath::Cos(ThetaCentr);
  double EtaCry = (Theta-ThetaCentr)/ThetaWidth;    
  if (EtaCry>0.5) EtaCry=0.5;
  if (EtaCry<-0.5) EtaCry=-0.5;
  //flip to take into account ECAL barrel symmetries:
  if (ietaclosest<0) EtaCry *= -1.;
  
  //-------------- end calculate local position -------------
  
  std::pair<double,double> etaphi(EtaCry,PhiCry);
  return etaphi;
  
  
  
}




