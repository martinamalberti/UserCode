#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

//---- CMSSW includes
#include "DataFormats/Math/interface/LorentzVector.h"


//---- root includes
#include "TH1.h"
#include "TH2.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define PI 3.14159265
#define TWOPI 6.28318530

using namespace std;

double deltaPhi(double phi1,double phi2) {
 double deltaphi = fabs(phi1-phi2);  
 if (deltaphi > TWOPI) deltaphi -= TWOPI;  
 if (deltaphi > PI ) deltaphi = TWOPI - deltaphi;  
 return deltaphi; 
}

//! main program
int main (int argc, char** argv)
{
  int mass = atoi(argv[1]);
  char outputFileName[100];
  sprintf(outputFileName,"outputHistos_WprimeM%d.root",mass);

  char inputFileName[100];
  sprintf(inputFileName,"/tmp/malberti/treeWprimeM%d.root",mass);

  std::string TreeName       = "myanalysis/tTreeUtilities";
  


  // TREE
  vector<double>  *MCelePx_;
  vector<double>  *MCelePy_;
  vector<double>  *MCelePz_;
  vector<double>  *MCeleEta_;
  vector<double>  *MCelePhi_;
  vector<int>     *MCelePid_;
  Int_t           NeleCand_;
  vector<int>     *isMCmatched_;
  vector<double>  *elePx_;
  vector<double>  *elePy_;
  vector<double>  *elePz_;
  vector<double>  *eleE_;
  vector<double>  *eleEt_;
  vector<double>  *eleEta_;
  vector<double>  *elePhi_;
  vector<float>   *eleId_;
  vector<float>   *eleCharge_;
  vector<double>  *eleTrkIsol_;
  vector<double>  *eleEcalIsol_;
  vector<double>  *eleHcalIsolD1_;
  vector<double>  *eleHcalIsolD2_;
  Double_t        uncorrMet_;
  Double_t        uncorrMex_;
  Double_t        uncorrMey_;
  Double_t        uncorrMetPhi_;
  Double_t        Met_;
  Double_t        Mex_;
  Double_t        Mey_;
  Double_t        MetPhi_;
  Int_t           Njets_;
  vector<double>  *jetPx_;
  vector<double>  *jetPy_;
  vector<double>  *jetPz_;
  vector<double>  *jetPt_;
  vector<double>  *jetEta_;
  vector<double>  *jetPhi_;
  Int_t           HLTEle15_;
  Int_t           HLTLooseIsoEle15_;
  Int_t           HLT_EM80_;
  Int_t           HLT_EM200_;
  Int_t           HLTPhoton15_;
  Int_t           HLTPhoton25_;
  

  MCelePx_  = new std::vector<double>; // MCtruth
  MCelePy_  = new std::vector<double>; // MCtruth
  MCelePz_  = new std::vector<double>; // MCtruth
  MCeleEta_ = new std::vector<double>;    // MCtruth
  MCelePhi_ = new std::vector<double>;    // MCtruth
  MCelePid_ = new std::vector<int>;    // MCtruth
  
  isMCmatched_  = new std::vector<int>; // MC truth matching 
  elePx_  = new std::vector<double>; // track momentum 
  elePy_  = new std::vector<double>; // track momentum 
  elePz_  = new std::vector<double>; // track momentum 
  eleE_ = new std::vector<double>;    // SC energy
  eleEt_ = new std::vector<double>;   // SC transverse energy
  eleEta_ = new std::vector<double>;  // SC pseudorapidity
  elePhi_ = new std::vector<double>;  // SC phi
  eleCharge_ = new std::vector<float>;  // electron charge
  eleId_ = new std::vector<float>;      // electron ID
  eleTrkIsol_ = new std::vector<double>; // track isolation
  eleEcalIsol_ = new std::vector<double>; // ecal isolation
  eleHcalIsolD1_ = new std::vector<double>; // hcal isolation
  eleHcalIsolD2_ = new std::vector<double>; // hcal isolation

  jetPx_  = new std::vector<double>; // 
  jetPy_  = new std::vector<double>; // 
  jetPz_  = new std::vector<double>; // 
  jetPt_  = new std::vector<double>; // 
  jetEta_  = new std::vector<double>; // 
  jetPhi_  = new std::vector<double>; // 

  // initialize tree
  TChain * fChain = new TChain (TreeName.c_str()) ;
  fChain->Add(inputFileName);
  fChain->SetBranchAddress("MCelePx_", &MCelePx_);
  fChain->SetBranchAddress("MCelePy_", &MCelePy_);
  fChain->SetBranchAddress("MCelePz_", &MCelePz_);
  fChain->SetBranchAddress("MCeleEta_", &MCeleEta_);
  fChain->SetBranchAddress("MCelePhi_", &MCelePhi_);
  fChain->SetBranchAddress("MCelePid_", &MCelePid_);
  fChain->SetBranchAddress("NeleCand_", &NeleCand_);
  fChain->SetBranchAddress("isMCmatched_", &isMCmatched_);
  fChain->SetBranchAddress("elePx_", &elePx_);
  fChain->SetBranchAddress("elePy_", &elePy_);
  fChain->SetBranchAddress("elePz_", &elePz_);
  fChain->SetBranchAddress("eleE_", &eleE_);
  fChain->SetBranchAddress("eleEt_", &eleEt_);
  fChain->SetBranchAddress("eleEta_", &eleEta_);
  fChain->SetBranchAddress("elePhi_", &elePhi_);
  fChain->SetBranchAddress("eleId_", &eleId_);
  fChain->SetBranchAddress("eleCharge_", &eleCharge_);
  fChain->SetBranchAddress("eleTrkIsol_", &eleTrkIsol_);
  fChain->SetBranchAddress("eleEcalIsol_", &eleEcalIsol_);
  fChain->SetBranchAddress("eleHcalIsolD1_", &eleHcalIsolD1_);
  fChain->SetBranchAddress("eleHcalIsolD2_", &eleHcalIsolD2_);
  fChain->SetBranchAddress("uncorrMet_",&uncorrMet_);
  fChain->SetBranchAddress("uncorrMex_",&uncorrMex_);
  fChain->SetBranchAddress("uncorrMey_",&uncorrMey_);
  fChain->SetBranchAddress("uncorrMetPhi_",&uncorrMetPhi_);
  fChain->SetBranchAddress("Met_", &Met_);
  fChain->SetBranchAddress("Mex_", &Mex_);
  fChain->SetBranchAddress("Mey_", &Mey_);
  fChain->SetBranchAddress("MetPhi_", &MetPhi_);
  fChain->SetBranchAddress("Njets_", &Njets_);
  fChain->SetBranchAddress("jetPx_", &jetPx_);
  fChain->SetBranchAddress("jetPy_", &jetPy_);
  fChain->SetBranchAddress("jetPz_", &jetPz_);
  fChain->SetBranchAddress("jetPt_", &jetPt_);
  fChain->SetBranchAddress("jetEta_", &jetEta_);
  fChain->SetBranchAddress("jetPhi_", &jetPhi_);
  fChain->SetBranchAddress("HLTEle15_", &HLTEle15_);
  fChain->SetBranchAddress("HLTLooseIsoEle15_", &HLTLooseIsoEle15_);
  fChain->SetBranchAddress("HLT_EM80_", &HLT_EM80_);
  fChain->SetBranchAddress("HLT_EM200_", &HLT_EM200_);
  fChain->SetBranchAddress("HLTPhoton15_", &HLTPhoton15_);
  fChain->SetBranchAddress("HLTPhoton25_", &HLTPhoton25_);
 
  int nEntries = fChain->GetEntries();
  std::cout << "------> Number of events "<< nEntries << " <------\n" ;

  // file to save histos
  TFile *fout = new TFile(outputFileName,"recreate");
  //INITIALIZING HISTOGRAMS
  TH1F *hNGoodElectrons = new TH1F("hNGoodElectrons","hNGoodElectrons",5,0,5);
  TH1F *het = new TH1F("het","het",1500,0,1500);
  TH1F *hmet = new TH1F("hmet","hmet",1500,0,1500);
  TH1F *hmetUncorr = new TH1F("hmetUncorr","hmetUncorr",1500,0,1500);
  TH1F *hmt = new TH1F("hmt","hmt",3000,0,3000);

  TH1F *hetWithJetVeto = new TH1F("hetWithJetVeto","hetWithJetVeto",1500,0,1500);
  TH1F *hmetWithJetVeto = new TH1F("hmetWithJetVeto","hmetWithJetVeto",1500,0,1500);
  TH1F *hmtWithJetVeto = new TH1F("hmtWithJetVeto","hmtWithJetVeto",3000,0,3000);

  TH1F *hetWithMetCut = new TH1F("hetWithMetCut","hetWithMetCut",1500,0,1500);
  TH1F *hmetWithMetCut = new TH1F("hmetWithMetCut","hmetWithMetCut",1500,0,1500);
  TH1F *hmtWithMetCut = new TH1F("hmtWithMetCut","hmtWithMetCut",3000,0,3000);

  TH1F *hjetPt = new TH1F("hjetPt","hjetPt",500,0,500);
  TH1F *hNjets = new TH1F("hNjets","hNjets",50,0,50);
  TH1F *hDphiJetEle = new TH1F("hDphiJetEle","hDphiJetEle",400,0,4);

  TH1F *hEtOverMet = new TH1F("hEtOverMet","hEtOverMet",1000,0,10);
  TH1F *hDphiMetEle = new TH1F("hDphiMetEle","hDphiMetEle",400,0,4);

  TH1F *hMHT = new TH1F("hMHT","hMHT",1000,0,1000); // Greg's suggestion

  
  TH1F *hR = new TH1F("hR","hR",100,0,10);

  float counter[10] = {0,0,0,0,0,0,0,0,0,0};
  float EtaCutEB    = 1.4442;
  float EtaCutEE    = 1.560;
  float combIsoCut  = 0.;
  float combinedIso = 0.;


  float lumi = 100;
  float xsec = 0;
  if (mass == 1000) xsec = 1.55  ; // Wprime , M = 1000
  if (mass == 1500) xsec = 0.24; // Wprime , M = 1500
  if (mass == 2000) xsec = 0.051 ; // Wprime , M = 2000

  float w = xsec*lumi/(float)nEntries;

  // START LOOP OVER ENTRIES
  for (int entry = 0 ; entry < nEntries ; ++entry) {
    fChain->GetEntry(entry) ;
    if(entry%10000==0) std::cout << "------> reading entry " << entry << " <------\n" ;

    counter[0]+=w;

    // HLT selection
    if (HLTEle15_==0) continue;
    counter[1]+=w;


    int nGoodElectrons = 0;
    int chosenEle = 0;

    //start loop over electron candidates 
    for (int i=0; i<NeleCand_; i++){

      // keep only electrons in ECAL fiducial volume
      if ( fabs(eleEta_->at(i)) > EtaCutEB &&  fabs(eleEta_->at(i)) < EtaCutEE) continue;

      // keep only electrons with ET > 30 GeV
      //if (eleEt_->at(i) < 30) continue ;
      
      // electron ID
      if (eleId_->at(i)==0) continue;

      // isolation
      if  ( fabs(eleEta_->at(i)) < EtaCutEB &&  eleTrkIsol_->at(i) > 7.5)  continue;
      if  ( fabs(eleEta_->at(i)) > EtaCutEE &&  eleTrkIsol_->at(i) > 15.)  continue;
      
      combinedIso = eleEcalIsol_->at(i) + eleHcalIsolD1_->at(i);
      if (fabs(eleEta_->at(i)) < EtaCutEB ) combIsoCut = 3 + 0.02*eleEt_->at(i);
      if (fabs(eleEta_->at(i)) > EtaCutEE && eleEt_->at(i)< 50) combIsoCut = 5.5;
      if (fabs(eleEta_->at(i)) > EtaCutEE && eleEt_->at(i)> 50) combIsoCut = 5.5 + 0.05*(eleEt_->at(i)-50);
      if (combinedIso > combIsoCut ) continue;
    
      if ( eleHcalIsolD2_->at(i) > 0.5 && fabs(eleEta_->at(i)) > EtaCutEE) continue;
      
      nGoodElectrons++;
      chosenEle = i;
      
    } // end loop over electron candidates     
    
    hNGoodElectrons->Fill( nGoodElectrons,w);

    if ( nGoodElectrons != 1 ) continue;
    counter[2]+=w;

    double et = eleEt_->at(chosenEle);
    double pt = sqrt(elePx_->at(chosenEle)*elePx_->at(chosenEle)+elePy_->at(chosenEle)*elePy_->at(chosenEle));

 
    double cphi = (elePx_->at(chosenEle)*Mex_ + elePy_->at(chosenEle)*Mey_)/(Met_*pt);
    double mt = sqrt(2 * et * Met_*(1-cphi));
   
    het  -> Fill(et,w);
    hmt  -> Fill(mt,w);
    hmet -> Fill(Met_,w);
    hmetUncorr -> Fill(uncorrMet_,w);
    
    
    // leading jet pT
    double maxPt = -1.;
    int jmax = 0;
    int njet = 0 ;
    double mh[3]={0.,0.,0.};
    for (int j=0; j<Njets_; j++){
      double deta = jetEta_->at(j) -  eleEta_->at(chosenEle);
      double dphi = deltaPhi(jetPhi_->at(j),elePhi_->at(chosenEle));
      double dR = sqrt(deta*deta+dphi*dphi);
      hR->Fill(dR);
      if (dR<0.1) continue;
      if (jetPt_->at(j) > maxPt) {
	jmax = j ; 
	maxPt = jetPt_->at(j);
      } 
      if (jetPt_->at(j) > 15 && fabs(jetEta_->at(j))<3 ) njet++;
      mh[0]+=jetPx_->at(j);
      mh[1]+=jetPy_->at(j);
      mh[2]+=jetPz_->at(j);
    }
    
    double mht = sqrt(mh[0]*mh[0]+mh[1]*mh[1]);

    hjetPt->Fill(maxPt,w);
    hNjets->Fill(njet,w);
    hDphiJetEle->Fill(deltaPhi(jetPhi_->at(jmax),elePhi_->at(chosenEle)),w);

    hEtOverMet -> Fill(et/Met_,w);
    hDphiMetEle-> Fill(deltaPhi(MetPhi_,elePhi_->at(chosenEle)),w);

    hMHT ->Fill(mht,w);

    if ( maxPt < 70 ) {
      counter[3]+=w;
      hetWithJetVeto ->Fill(et,w);
      hmetWithJetVeto ->Fill(Met_,w);
      hmtWithJetVeto ->Fill(mt,w);
    } 

    if (et/Met_>0.4 && et/Met_<1.5 && deltaPhi(MetPhi_,elePhi_->at(chosenEle))>2.5 ){
      counter[4]+=w;
      hetWithMetCut ->Fill(et,w);
      hmetWithMetCut ->Fill(Met_,w);
      hmtWithMetCut ->Fill(mt,w);
    }
   
  }// END LOOP OVER ENTRIES



  cout<<"Selection        # events          Efficiency"               <<endl;
  cout<<"No selections    " << float(counter[0]) << "  " << float(counter[0])/float(counter[0]) <<endl; 
  cout<<"HLT              " << float(counter[1]) << "  " << float(counter[1])/float(counter[0]) <<endl; 
  cout<<"elID+Isolation   " << float(counter[2]) << "  " << float(counter[2])/float(counter[1]) <<endl; 
  cout<<"Jet veto         " << float(counter[3]) << "  " << float(counter[3])/float(counter[2]) <<endl; 
  cout<<"Total eff.       " << float(counter[3]) << "  " << float(counter[3])/float(counter[0]) <<endl; 
  cout<<endl;

  int xbin ;
  int nbins = het->GetNbinsX();
  xbin = het->FindBin(200);
  cout<<"Number of events with ET > 200 GeV :"<< het->Integral(xbin,nbins)<<endl;
  xbin = het->FindBin(300);
  cout<<"Number of events with ET > 300 GeV :"<< het->Integral(xbin,nbins)<<endl;

  nbins = hmt->GetNbinsX();
  xbin = hmt->FindBin(200);
  cout<<"Number of events with MT > 200 GeV :"<< hmt->Integral(xbin,nbins)<<endl;
  xbin = hmt->FindBin(500);
  cout<<"Number of events with MT > 500 GeV :"<< hmt->Integral(xbin,nbins)<<endl;

  cout<<"With jet veto"<<endl;
  nbins = hetWithJetVeto->GetNbinsX();
  xbin = hetWithJetVeto->FindBin(200);
  cout<<"Number of events with ET > 200 GeV :"<< hetWithJetVeto->Integral(xbin,nbins)<<endl;
  xbin = hetWithJetVeto->FindBin(300);
  cout<<"Number of events with ET > 300 GeV :"<< hetWithJetVeto->Integral(xbin,nbins)<<endl;

  nbins = hmtWithJetVeto->GetNbinsX();
  xbin = hmtWithJetVeto->FindBin(200);
  cout<<"Number of events with MT > 200 GeV :"<< hmtWithJetVeto->Integral(xbin,nbins)<<endl;
  xbin = hmtWithJetVeto->FindBin(500);
  cout<<"Number of events with MT > 500 GeV :"<< hmtWithJetVeto->Integral(xbin,nbins)<<endl;

  cout<<"With MET cut"<<endl;
  nbins = hetWithMetCut->GetNbinsX();
  xbin = hetWithMetCut->FindBin(200);
  cout<<"Number of events with ET > 200 GeV :"<< hetWithMetCut->Integral(xbin,nbins)<<endl;
  xbin = hetWithMetCut->FindBin(300);
  cout<<"Number of events with ET > 300 GeV :"<< hetWithMetCut->Integral(xbin,nbins)<<endl;

  nbins = hmtWithMetCut->GetNbinsX();
  xbin = hmtWithMetCut->FindBin(200);
  cout<<"Number of events with MT > 200 GeV :"<< hmtWithMetCut->Integral(xbin,nbins)<<endl;
  xbin = hmtWithMetCut->FindBin(500);
  cout<<"Number of events with MT > 500 GeV :"<< hmtWithMetCut->Integral(xbin,nbins)<<endl;


  // write histos 
  fout->Write();
  fout->Close();

}
