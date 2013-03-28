#ifndef ntpleUtils_h
#define ntpleUtils_h

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include "functional"
#include <utility>



#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TChain.h"
#include "TVector3.h"
#include "Math/Vector4D.h"
#include "ConfigParser.h"
#include "TVectorF.h"





/** get the number of events from a list of files */
std::map<int, int> GetTotalEvents(const std::string& histoName, const std::string& inputFileList);

/** fill a chain from a list of files */
bool FillChain(TChain& chain, const std::string& inputFileList);

/** get the parameters from a congiguration file */
int parseConfigFile (const TString& config) ;






/** compute delta phi */
double deltaPhi (const double& phi1, const double& phi2);

/** compute delta eta */
double deltaEta (const double& eta1, const double& eta2);

/** compute delta R */
double deltaR (const double& eta1, const double& phi1,
               const double& eta2, const double& phi2);





/** define the map of all possible matching */
template <class T1, class T2>
std::map<float, std::pair<int, int> > MatchingDRMap(const std::vector<T1>& collection1,
                                                    const std::vector<T2>& collection2)
{
  // define map of all possible DRs
  typedef std::map<float, std::pair<int, int> > DRMap;
  DRMap myDRMap;
  myDRMap.clear();  
  // compute all DRs 
  for(unsigned int it1 = 0; it1 <  collection1.size(); ++it1)
  {
    for(unsigned int it2 = 0; it2 < collection2.size(); ++it2)
    {
      float DR = deltaR((collection1.at(it1)).Eta(), (collection1.at(it1)).Phi(),
                        (collection2.at(it2)).Eta(), (collection2.at(it2)).Phi());
      std::pair<int, int> dummy(it1, it2);
      myDRMap[DR] = dummy;
    }
  }
  return myDRMap;
}

/** define the map of all possible matching */
template <class T1, class T2>
int GetMatching(const std::vector<T1>& collection1, //---- RECO
                const std::vector<T2>& collection2, //---- MC
                const float& DRMax,
                float ptRatioMin,
                float ptRatioMax,
                std::vector<int>* matchIt1 = 0) //---- index from RECO that matches with MC

{
  // get the DR map between two collection of particles
  typedef std::map<float, std::pair<int, int> > DRMap;
  DRMap myDRMap = MatchingDRMap(collection1, collection2);
  
    // define boolean flags to avoid double usage of the same particle
  std::vector<bool> isUsed1;
  for(unsigned int i = 0; i < collection1.size(); ++ i)
    isUsed1.push_back(false);
  
  std::vector<bool> isUsed2;
  for(unsigned int i = 0; i < collection2.size(); ++ i)
    isUsed2.push_back(false);
  
  
  // intialize the vector which will store the result of the matching
  if(matchIt1 != 0)
  {
    (*matchIt1).clear();
    for(unsigned int i = 0; i < collection2.size(); ++i)
      (*matchIt1).push_back(-1);
  }
  
  // loop over the DRmap to get the smallest DR matchings
  unsigned int nMatching = 0;
  
  for(DRMap::const_iterator mapIt = myDRMap.begin();
      mapIt != myDRMap.end(); ++mapIt)
  {
    int it1 = (mapIt -> second).first;
    int it2 = (mapIt -> second).second;
    
    float DR = mapIt -> first;
    
    TVector3 momentum1(collection1.at(it1).Px(),
                       collection1.at(it1).Py(),
                       collection1.at(it1).Pz());
    TVector3 momentum2(collection2.at(it2).Px(),
                       collection2.at(it2).Py(),
                       collection2.at(it2).Pz());
    
    if( 1. * collection1.at(it1).Pt() / (collection2.at(it2)).Pt() < ptRatioMin) 
      continue;
    if( 1. * collection1.at(it1).Pt() / (collection2.at(it2)).Pt() > ptRatioMax) 
      continue;
    
    if( (DR <= DRMax) && (isUsed1.at(it1) == false) && (isUsed2.at(it2) == false) )
    {
      isUsed1.at(it1) = true; 
      isUsed2.at(it2) = true; 
      ++nMatching;
      
      if(matchIt1 != 0)
        (*matchIt1).at(it2) = it1;
    }
    
    if( (nMatching == collection2.size()) || (DR > DRMax) )
      break;    
  }
  
  return nMatching;
}



/** Zeppenfeld Jet Veto */
int getZepp(std::vector<ROOT::Math::XYZTVector>& jets,
	   int q1,
	   int q2,
	   const double& EtMin,
	   const double& zeppMax,
	   const std::vector<int>* blacklist = 0);
	   

/** Central Jet Veto */
int getCJV(std::vector<ROOT::Math::XYZTVector>& jets,
	      int q1,
	      int q2,
	      const double& EtMin,
	      const std::vector<int>* blacklist = 0);


/** Jet Veto */
int getJV(std::vector<ROOT::Math::XYZTVector>& jets,
	      const double& EtMin,
	      const std::vector<int>* blacklist = 0);

	      
/** select jet pairs */
double SelectJets(std::vector<int>& it, std::vector<ROOT::Math::XYZTVector>& jets,
                  const std::string& method,
                  const double& etMin,
                  const std::vector<int>* blacklist = 0);

/** select leptons */
int SelectLepton(std::vector<ROOT::Math::XYZTVector>& leptons,
                 const std::string& method,
                 const double& ptMin,
                 const std::vector<int>* blacklist = 0);


/** select pair of objects with a given invariant mass*/
double SelectResonance(std::vector<int>& it, std::vector<ROOT::Math::XYZTVector>& objects,
		       const double& mass,
		       const double& ptMin,
		       const std::vector<int>* blacklist = 0);
 

/** select pair of objects with a given invariant mass and opposite charge*/
double SelectResonanceOppositeCharge(std::vector<int>& it,
		       std::vector<ROOT::Math::XYZTVector>& objects,
		       std::vector<float>& charge,	     
		       const double& mass,
		       const double& ptMin,
		       const std::vector<int>* blacklist = 0);


/** build combinations of n jets */
int Build4JetCombinations(std::vector<std::vector<int> >& comb, const int& nJets);

/** build combinations of n jets */
void Print4JetCombination(const std::vector<int>& combination);

/** build combinations (2 jets) of n jets */
int Build2JetCombinations(std::vector<std::vector<int> >& comb, const int& nJets);

/** smart profiling by double averaging */
TH1D * smartProfileX (TH2F * strip, double width) ; 

/** smart profiling by fixing gaussian parameters and
range from a first averaging */
TH1D * smartGausProfileX (TH2F * strip, double width) ; 

/** smart profiling by double averaging */
TH1D * smartProfileY (TH2F * strip, double width) ; 

/** smart profiling by fixing gaussian parameters and
range from a first averaging */
TH1D * smartGausProfileY (TH2F * strip, double width) ; 

/** 68% and 95% bands using Neyman intervals 
and Feldman-Cousins principle */
std::vector<double> getSigmaBands_FeldmanCousins (const TH1 & histo) ;
 
/** smart profiling Y using Neyman intervals 
and Feldman-Cousins principle */
std::vector<TH1D*> smartGausProfileY_FeldmanCousins (TH2F * strip, double width);
  
/** smart profiling X using Neyman intervals 
and Feldman-Cousins principle */
std::vector<TH1D*> smartGausProfileX_FeldmanCousins (TH2F * strip, double width);

/**
define a band in the histogram such that for any slice at fixed x-bin 
from the histogram in the Y a minimal and maximal points are determined
by the funcional getLimit
*/
template <class T> std::vector<std::vector<double> > //PG three vectors: x axis, lower line, upper line
getBand_integrY (TH2F & histo2D, T getLimit)
{
 int cut = 0 ; // minimum number of entries per fitted bin
 int nbins = histo2D.GetXaxis ()->GetNbins () ;
 int binmin = 1 ;
 int ngroup = 1 ; // bins per step
 int binmax = nbins ;
 
 std::vector<double> dummy ;
 std::vector<std::vector<double> > out (3, dummy) ;
 
 // loop over the 2D histo bins
 for (int bin=binmin ; bin <= binmax ; bin += ngroup) 
 {
  TH1D * h1_dummy = histo2D.ProjectionY ("_temp", bin, bin + ngroup - 1 , "e") ;
  if (h1_dummy == 0) continue ;
  int nentries = Int_t (h1_dummy->GetEntries ()) ;
  if (nentries == 0 || nentries < cut) {delete h1_dummy ; continue ;}   
  std::pair<double, double> limits = getLimit (*h1_dummy) ;
  out.at (0).push_back (histo2D.GetXaxis ()->GetBinCenter (bin+0.5)) ;
  out.at (1).push_back (limits.first) ;
  out.at (2).push_back (limits.second) ;
  delete h1_dummy ;
 } // loop over the bins
 return out ;
}


/**
define a band in the histogram such that for any slice at fixed y-bin 
from the histogram in the X a minimal and maximal points are determined
by the funcional getLimit
*/
template <class T> std::vector<std::vector<double> > //PG three vectors: x axis, lower line, upper line
getBand_integrX (TH2F & histo2D, T getLimit)
{
 int cut = 0 ; // minimum number of entries per fitted bin
 int nbins = histo2D.GetYaxis ()->GetNbins () ;
 int binmin = 1 ;
 int ngroup = 1 ; // bins per step
 int binmax = nbins ;
 
 std::vector<double> dummy ;
 std::vector<std::vector<double> > out (3, dummy) ;
 
 // loop over the 2D histo bins
 for (int bin=binmin ; bin <= binmax ; bin += ngroup) 
 {
  TH1D * h1_dummy = histo2D.ProjectionX ("_temp", bin, bin + ngroup - 1 , "e") ;
  if (h1_dummy == 0) continue ;
  int nentries = Int_t (h1_dummy->GetEntries ()) ;
  if (nentries == 0 || nentries < cut) {delete h1_dummy ; continue ;}   
  std::pair<double, double> limits = getLimit (*h1_dummy) ;
  out.at (0).push_back (histo2D.GetYaxis ()->GetBinCenter (bin+0.5)) ;
  out.at (1).push_back (limits.first) ;
  out.at (2).push_back (limits.second) ;
  delete h1_dummy ;
 } // loop over the bins
 return out ;
}



/**
get the points on the definition set of a TH1F such that the region contains 68% of the area
and the tails contain the same 15.8%
*/
struct getLimit_sameTails : public std::unary_function <const TH1D &, std::pair<double, double> > {
 std::pair<double, double> operator() (const TH1D & histo);
};

/**
get the points on the definition set of a TH1F such that the region contains 68% of the area
and the size of the region is minimal (Neyman intervals - Feldman Cousins)
*/
struct getLimit_FC : public std::unary_function <const TH1D &, std::pair<double, double> > {
 std::pair<double, double> operator() (const TH1D & histo);
};


/**
build a TGraphError starting from a TH1F
*/
TGraphErrors buildGEfromH (const TH1D & histo) ;

#endif
