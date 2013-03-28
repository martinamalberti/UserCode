#ifndef stdHisto_h
#define stdHisto_h

#include "treeReader.h"
#include "hFactory.h"

#include <iostream>
#include <vector>

#include "Math/GenVector/VectorUtil.h"



class stdHisto
{
 public:
  
  //!ctor
  stdHisto(const int& nStep,
           const std::string& outFileName,
           treeReader* reader = NULL);
  
  //!dtor
  ~stdHisto();
  
  
  //! methods
  
  // add histograms
  void Add1(const std::string& histoName,
            const int& nStep);
  
  void Add2(const std::string& histoName,
            const int& nStep);

  void Add1Float(const std::string& histoName,
           const int& nStep,
	   int nbins,
	   double min,
	   double max);

  // fill histograms
  void Fill1(const std::string& histoName,
             const std::string& branchName,
             const int& nStep,
             std::vector<int>* selectionIt = NULL);
  
  void Fill1(const std::vector<ROOT::Math::XYZTVector>& vet,
             const std::string& histoName,
             const int& step);
  
  void Fill1(const ROOT::Math::XYZTVector& p,
             const std::string& histoName,
             const int& step);
  
  void Fill2(const std::string& histoName,
             const std::string& branchName,
             const int& nStep,
             const int& it1, const int& it2);
  
  void Fill2(const ROOT::Math::XYZTVector& v1,
             const ROOT::Math::XYZTVector& v2, 
             const std::string& histoName,
             const int& step);
  
  void Fill1Float(const std::string& histoName,
             const std::string& branchName,
             const int& nStep,
             std::vector<int>* selectionIt = NULL);
  
  void Fill1Float(const std::string& histoName,
             const double& value,
             const int& nStep);
  
 private:
  
  int m_nStep;
  treeReader* m_reader; 
  std::string m_outFileName;
  
  hFactory* m_hFactory;
  
};

#endif
