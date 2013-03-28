#include "TFile.h"
#include "TH1.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>


class VertexIdEfficiencyReweighting {
 private:
  TH1F* h;
 public:
  VertexIdEfficiencyReweighting(std::string,std::string); 
  ~VertexIdEfficiencyReweighting();
  double GetWeight(float);
};

VertexIdEfficiencyReweighting::VertexIdEfficiencyReweighting(std::string fname, std::string hname) {
 
  TFile *f = TFile::Open(fname.c_str(), "READ");
  std::cout << "Reading data/MC scale factors from :"<<  fname.c_str()  << std::endl;
  h = (TH1F*)f->Get(hname.c_str());
  return; 
}

VertexIdEfficiencyReweighting::~VertexIdEfficiencyReweighting() { return;}

double VertexIdEfficiencyReweighting::GetWeight(float pt){
  float w = 1.;
  if (pt < 250.){
    float bin = h->FindBin(pt);
    w = h->GetBinContent(bin);
  }
  //std::cout << "vtxid SF = " << w << std::endl; 
  return (w);
}

