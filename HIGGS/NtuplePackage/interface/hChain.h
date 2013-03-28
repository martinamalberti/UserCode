#ifndef hChain_h
#define hChain_h

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TNtuple.h"

#include <vector>


struct hChain 
{
  hChain (TString baseName, TString baseTitle, 
          int nbins, double min, double max, int NUM) ;          
  ~hChain () ;
  
  void SetColors (std::vector<int> colors) ;
  void Fill (int i, double val) ;
  void SetBinContent (int i, int bin, double val) ;
  void Print (bool isLog = false, int rebin = 1, TString altName = "default") ;
  void PrintEach (bool isLog = false, int rebin = 1) ;
  void Normalize (int index) ;
  void Scale (int index, double factor) ;
  void Write (TFile & outputFile) ;
  void Write (const std::string& dirName, TFile & outputFile) ;
  unsigned int Size() { return m_histos.size(); };
        
  private :
  
    TString m_baseName ;
    std::vector <TH1F*> m_histos ;
//    std::vector <TNtuple*> m_ntuples ;

    double findNMin () ;  
    double findNMax () ;
} ;

#endif
