
///==== include ====

#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "Math/GenVector/VectorUtil.h"
#include "TRandom3.h"
#include <time.h>
#include <sstream>
#include "MyTest.h"

#include <iostream>

#include "TClonesArray.h"
#include "TMatrix.h"


//---HGG----------------------------------------------------------------------------------------------
void hggSelection (std::vector<ROOT::Math::XYZTVector>* mcV1, 
		   std::vector<ROOT::Math::XYZTVector>* mcV2, 
		   std::vector<ROOT::Math::XYZTVector>* photons ,
		   std::vector<ROOT::Math::XYZTVector>* photons_SC,
		   std::vector<float>* photons_r9,
		   int& passSelection, int& i1, int& i2
		   )
{
  float photonsMC_eta[2], photonsMC_phi[2];
  
  if ( mcV1->size() != 1 ||  mcV2->size() != 1 ||  photons->size() == 0){
    passSelection = 0;
    i1 = -100;
    i2=-100;
    return;
  }
  
  
  // MC matching 
  if ( mcV1->at(0).E() > mcV2->at(0).E()) {
    photonsMC_eta[0] = mcV1->at(0).Eta();
    photonsMC_eta[1] = mcV2->at(0).Eta();
    photonsMC_phi[0] = mcV1->at(0).Phi();
    photonsMC_phi[1] = mcV2->at(0).Phi();
  }
  else      {
    photonsMC_eta[1] = mcV1->at(0).Eta();
    photonsMC_eta[0] = mcV2->at(0).Eta();
    photonsMC_phi[1] = mcV1->at(0).Phi();
    photonsMC_phi[0] = mcV2->at(0).Phi();
  }
  
  int index1 = -100, index2 = -100;
  double dR_1_min = 10000.;
  double dR_2_min = 10000.;
  for(unsigned int u=0; u < photons->size(); u++)
    {
      double dR_1 = deltaR( photonsMC_eta[0], photonsMC_phi[0], photons_SC->at(u).eta(), photons_SC->at(u).phi() );
      if (dR_1 < dR_1_min) { dR_1_min = dR_1; index1 = u; }
      
      double dR_2 = deltaR( photonsMC_eta[1], photonsMC_phi[1], photons_SC->at(u).eta(), photons_SC->at(u).phi() );
      if (dR_2 < dR_2_min) { dR_2_min = dR_2; index2 = u; }
      
    }

  
  if (photons->at(index1).E() > photons->at(index2).E())
    {
      i1 = index1;
      i2 = index2;
    }
  else
    {
      i1 = index2;
      i2 = index1;
    }
  


  bool is_mcmatched   = (dR_1_min < 0.15) && (dR_2_min < 0.15);
  //bool is_unconverted = (photons_r9->at(index1) > 0.93) &&  (photons_r9->at(index2) > 0.93);
  bool is_unconverted = (photons_r9->at(index1) > 0.) &&  (photons_r9->at(index2) > 0.);
  bool pass_kincuts   = (photons->at(index1).pt()> 40.) && ( photons->at(index2).pt() > 30.) ;
  
  if( is_mcmatched && is_unconverted && pass_kincuts ) {
    passSelection = 1;
    return;
  }


}
//-------------------------------------------------------------------------------------------------




//---ZEE----------------------------------------------------------------------------------------------
void zeeSelection (std::vector<ROOT::Math::XYZTVector>* electrons ,
		   std::vector<float>* eleid,
		   int& passSelection, int& i1, int& i2
		   )

{
  passSelection = 0;
  int index1 = -100, index2 = -100;
  int ngood = 0;
  for(unsigned int uu = 0; uu < eleid->size(); uu++)
    {
      if (  eleid->at(uu) == 7 && index1 < 0 )      {index1 = uu; ngood++;}
      else if ( eleid->at(uu) == 7 && index2 < 0 )  {index2 = uu; ngood++;}
      else if ( eleid->at(uu) == 7 )                { ngood++;}
    }
  
  if ( ngood != 2){
    passSelection = 0;
    i1 = -100;
    i2 = -100;
    return;
  }
  
  if (electrons->at(index1).E() > electrons->at(index2).E())
    {
      i1 = index1;
      i2 = index2;
    }
  else
    {
      i1 = index2;
      i2 = index1;
    }
  

  ROOT::Math::XYZTVector v = electrons->at(index1) + electrons->at(index2);
  bool pass_mcut      = (fabs( v.M()  - 91. ) < 8. ) ;
  //bool pass_kincuts   = (electrons->at(index1).pt()> 30.) && ( electrons->at(index2).pt() > 24.) ;  
  bool pass_kincuts   = (electrons->at(index1).pt()> 0.) && ( electrons->at(index2).pt() > 0.) ;
  
  if( pass_mcut && pass_kincuts ) {
    passSelection = 1;    
    return ;
  }

}



//---ZMM----------------------------------------------------------------------------------------------
void zmumuSelection (std::vector<ROOT::Math::XYZTVector>* muons ,
		     std::vector<int>* muons_global,
		     std::vector<int>* muons_tracker,
		     std::vector<float>* muons_tkIsoR03,
		     std::vector<float>* muons_normalizedChi2 , 
		     std::vector<int>* muons_numberOfValidMuonHits,
		     std::vector<int>* muons_numberOfValidPixelHits,
		     std::vector<float>* muons_dxy_PV,
		     std::vector<float>* muons_dz_PV,
		     int& passSelection, int& i1, int& i2)
{
  passSelection = 0;
  int index1 = -100, index2 = -100;
  int ngood = 0;

  for( unsigned int uu = 0; uu < muons->size(); uu++)
    {
      float relIso  = muons_tkIsoR03->at(uu)/muons->at(uu).pt();
      bool goodmuon =  (muons->at(uu).pt() > 10   && 
			muons_global->at(uu)==1   && 
			muons_tracker->at(uu)==1  && 
			muons_normalizedChi2->at(uu)<10   &&
			muons_numberOfValidMuonHits->at(uu)>0 &&
			muons_numberOfValidPixelHits->at(uu)>0 &&
			muons_dxy_PV->at(uu)<0.2 &&
			muons_dz_PV->at(uu)<0.5  &&
			relIso<0.10);
      if ( goodmuon && index1 < 0 )   { index1 = uu; ngood++; }
      else if (goodmuon && index2 < 0 )   { index2 = uu; ngood++; }
      else if (goodmuon)                 { ngood++;}
    }
    
  if ( ngood != 2){
    passSelection = 0;
    i1 = -100;
    i2 = -100;
    return;
  }
  
  if (muons->at(index1).E() > muons->at(index2).E())
    {
      i1 = index1;
      i2 = index2;
    }
  else
    {
      i1 = index2;
      i2 = index1;
    }
  

  ROOT::Math::XYZTVector v = muons->at(index1) + muons->at(index2);
  bool pass_mcut      = (fabs( v.M() - 91. ) < 8. ) ;
  //bool pass_kincuts   = (muons->at(index1).pt()> 30.) && ( muons->at(index2).pt() > 24.) ;
  bool pass_kincuts   = (muons->at(index1).pt()> 0.) && ( muons->at(index2).pt() > 0.) ;
  
  if( pass_mcut && pass_kincuts ) {
    passSelection = 1;    
    return;
  }

}
//-------------------------------------------------------------------------------------------------
