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
 

  std::string TreeName       = "myanalysis/tTreeUtilities";


  // Reading input parameters

  //---- Input Variables ----
  std::string FileIn;
  std::string FileOut;
  float Xsec;
  float Lumi;
  int useJetVeto = 1;

  std::string inputName = argv[1];
  ifstream file(inputName.c_str());

  std::string variableNameFileIn  = "FileIn";
  std::string variableNameFileOut = "FileOut";
  std::string variableNameXsec    = "Xsec";
  std::string variableNameLumi    = "Lumi";
  std::string variableNameJetVeto = "JetVetoFlag";
  std::cerr << " Reading  " << inputName << " ... " << std::endl;

  file >>  variableNameFileIn  >> FileIn
       >>  variableNameFileOut >> FileOut
       >>  variableNameXsec    >> Xsec
       >>  variableNameLumi    >> Lumi
       >>  variableNameJetVeto >> useJetVeto;


  cout <<  variableNameFileIn  << "  "  << FileIn  << "  "
       <<  variableNameFileOut << "  "  << FileOut << "  "
       <<  variableNameXsec    << "  "  << Xsec    << "  "
       <<  variableNameLumi    << "  "  << Lumi    << "  "
       <<  variableNameJetVeto << "  "  << useJetVeto
       <<  endl;



  // TREE
  Int_t           NeleCand;
  vector<double>  *elePx;
  vector<double>  *elePy;
  vector<double>  *elePz;
  vector<double>  *eleE;
  vector<double>  *eleEt;
  vector<double>  *eleEta;
  vector<double>  *elePhi;
  vector<double>  *eleId;
  vector<double>  *eleTrkIsol;
  vector<double>  *eleEcalIsol;
  vector<double>  *eleHcalIsolD1;
  vector<double>  *eleHcalIsolD2;
  Double_t        uncorrMet;
  Double_t        uncorrMex;
  Double_t        uncorrMey;
  Double_t        uncorrMetPhi;
  Double_t        Met;
  Double_t        Mex;
  Double_t        Mey;
  Double_t        MetPhi;
  Int_t           Njets;
  vector<double>  *jetPx;
  vector<double>  *jetPy;
  vector<double>  *jetPz;
  vector<double>  *jetPt;
  vector<double>  *jetEta;
  vector<double>  *jetPhi;
  Int_t           HLTEle10;
  Int_t           HLTEle15;
  Int_t           HLTEle20;
  Int_t           HLTLooseIsoEle15;
  Int_t           HLTPhoton15;
  Int_t           HLTPhoton25;
  

  elePx  = new std::vector<double>; // track momentum 
  elePy  = new std::vector<double>; // track momentum 
  elePz  = new std::vector<double>; // track momentum 
  eleE   = new std::vector<double>;    // SC energy
  eleEt  = new std::vector<double>;   // SC transverse energy
  eleEta = new std::vector<double>;  // SC pseudorapidity
  elePhi = new std::vector<double>;  // SC phi
  eleId  = new std::vector<double>;      // electron ID
  eleTrkIsol    = new std::vector<double>; // track isolation
  eleEcalIsol   = new std::vector<double>; // ecal isolation
  eleHcalIsolD1 = new std::vector<double>; // hcal isolation
  eleHcalIsolD2 = new std::vector<double>; // hcal isolation

  jetPx  = new std::vector<double>; // 
  jetPy  = new std::vector<double>; // 
  jetPz  = new std::vector<double>; // 
  jetPt  = new std::vector<double>; // 
  jetEta = new std::vector<double>; // 
  jetPhi = new std::vector<double>; // 

  // initialize tree
  TChain * fChain = new TChain (TreeName.c_str()) ;

  fChain->Add(FileIn.c_str());
  int nEntries = fChain->GetEntries();
  std::cout << "------> Number of events "<< nEntries << " <------\n" ;
  float w = Xsec*Lumi/(float)nEntries;
  cout <<"xsec = " << Xsec  << "    N = " <<  nEntries << "     w =" << w << endl;

  fChain->SetBranchAddress("NeleCand", &NeleCand);
  fChain->SetBranchAddress("elePx", &elePx);
  fChain->SetBranchAddress("elePy", &elePy);
  fChain->SetBranchAddress("elePz", &elePz);
  fChain->SetBranchAddress("eleE", &eleE);
  fChain->SetBranchAddress("eleEt", &eleEt);
  fChain->SetBranchAddress("eleEta", &eleEta);
  fChain->SetBranchAddress("elePhi", &elePhi);
  fChain->SetBranchAddress("eleId", &eleId);
  fChain->SetBranchAddress("eleTrkIsol", &eleTrkIsol);
  fChain->SetBranchAddress("eleEcalIsol", &eleEcalIsol);
  fChain->SetBranchAddress("eleHcalIsolD1", &eleHcalIsolD1);
  fChain->SetBranchAddress("eleHcalIsolD2", &eleHcalIsolD2);
  fChain->SetBranchAddress("uncorrMet",&uncorrMet);
  fChain->SetBranchAddress("uncorrMex",&uncorrMex);
  fChain->SetBranchAddress("uncorrMey",&uncorrMey);
  fChain->SetBranchAddress("uncorrMetPhi",&uncorrMetPhi);
  fChain->SetBranchAddress("Met", &Met);
  fChain->SetBranchAddress("Mex", &Mex);
  fChain->SetBranchAddress("Mey", &Mey);
  fChain->SetBranchAddress("MetPhi", &MetPhi);
  fChain->SetBranchAddress("Njets", &Njets);
  fChain->SetBranchAddress("jetPx", &jetPx);
  fChain->SetBranchAddress("jetPy", &jetPy);
  fChain->SetBranchAddress("jetPz", &jetPz);
  fChain->SetBranchAddress("jetPt", &jetPt);
  fChain->SetBranchAddress("jetEta", &jetEta);
  fChain->SetBranchAddress("jetPhi", &jetPhi);
  fChain->SetBranchAddress("HLTEle10", &HLTEle10);
  fChain->SetBranchAddress("HLTEle15", &HLTEle15);
  fChain->SetBranchAddress("HLTEle20", &HLTEle20);
  fChain->SetBranchAddress("HLTLooseIsoEle15", &HLTLooseIsoEle15);
  fChain->SetBranchAddress("HLTPhoton15", &HLTPhoton15);
  fChain->SetBranchAddress("HLTPhoton25", &HLTPhoton25);
 
  
  // FILE to save histos
  TFile *fout = new TFile(FileOut.c_str(),"recreate");
  
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

  float counter[10] = {0,0,0,0,0,0,0,0,0,0};
  float EtaCutEB    = 1.4442;
  float EtaCutEE    = 1.560;
  float combIsoCut  = 0.;
  float combinedIso = 0.;

  
  // START LOOP OVER ENTRIES
  for (int entry = 0 ; entry < nEntries ; ++entry) {
    fChain->GetEntry(entry) ;
    if(entry%10000==0) std::cout << "------> reading entry " << entry << " <------\n" ;

    counter[0]+=w;

    // HLT selection
    if (HLTLooseIsoEle15==0) continue;
    counter[1]+=w;

    int nGoodElectrons = 0;
    int chosenEle = 0;
    
    //start loop over electron candidates 
    for (int i=0; i < NeleCand; i++){
    
      // keep only electrons in ECAL fiducial volume
      if ( fabs(eleEta->at(i)) > EtaCutEB &&  fabs(eleEta->at(i)) < EtaCutEE) continue;
     

      // keep only electrons with ET > 30 GeV
      if (eleEt->at(i) < 30) continue ;
      
      // electron ID
      if (eleId->at(i)==0) continue;
     
      // isolation
      if  ( fabs(eleEta->at(i)) < EtaCutEB &&  eleTrkIsol->at(i) > 7.5)  continue;
      if  ( fabs(eleEta->at(i)) > EtaCutEE &&  eleTrkIsol->at(i) > 15.)  continue;
  
      combinedIso = eleEcalIsol->at(i) + eleHcalIsolD1->at(i);
      if (fabs(eleEta->at(i)) < EtaCutEB ) combIsoCut = 2 + 0.03*eleEt->at(i);
      if (fabs(eleEta->at(i)) > EtaCutEE && eleEt->at(i)< 50) combIsoCut = 2.5;
      if (fabs(eleEta->at(i)) > EtaCutEE && eleEt->at(i)> 50) combIsoCut = 2.5 + 0.03*(eleEt->at(i)-50);
      if (combinedIso > combIsoCut ) continue;

      if ( eleHcalIsolD2->at(i) > 0.5 && fabs(eleEta->at(i)) > EtaCutEE) continue;

      nGoodElectrons++;
      chosenEle = i;
      
    } // end loop over electron candidates     
    

    hNGoodElectrons->Fill( nGoodElectrons,w );

    if ( nGoodElectrons != 1 ) continue;
    counter[2]+=w;

    double et = eleEt->at(chosenEle);
    double pt = sqrt(elePx->at(chosenEle)*elePx->at(chosenEle)+elePy->at(chosenEle)*elePy->at(chosenEle));
    double cphi = (elePx->at(chosenEle)*Mex + elePy->at(chosenEle)*Mey)/(Met*pt);
    double mt = sqrt(2 * et * Met*(1-cphi));
        
    het  -> Fill(et,w);
    hmt  -> Fill(mt,w);
    hmet -> Fill(Met,w);
    hmetUncorr -> Fill(uncorrMet,w);
    
    
    // leading jet pT
    double maxPt = -1.;
    int jmax = -1;
    int njet = 0 ;
    for (int j=0; j<Njets; j++){
      double deta = jetEta->at(j) -  eleEta->at(chosenEle);
      double dphi = deltaPhi(jetPhi->at(j),elePhi->at(chosenEle));
      double dR = sqrt(deta*deta+dphi*dphi);
      if (dR<0.1) continue;
      if (jetPt->at(j) > maxPt) {
	jmax = j ; 
	maxPt = jetPt->at(j);
      } 
    }

    if (jmax !=-1){
      hjetPt->Fill(maxPt,w);
      hDphiJetEle->Fill(deltaPhi(jetPhi->at(jmax),elePhi->at(chosenEle)),w);
    }
    
    hNjets->Fill(njet,w);
    if ( maxPt < 100 ) {
      counter[3]+=w;
      hetWithJetVeto ->Fill(et,w);
      hmetWithJetVeto ->Fill(Met,w);
      hmtWithJetVeto ->Fill(mt,w);
    }    
    

    hEtOverMet -> Fill(et/Met,w);
    hDphiMetEle-> Fill(deltaPhi(MetPhi,elePhi->at(chosenEle)),w);
    if (et/Met>0.4 && et/Met<1.5 && deltaPhi(MetPhi,elePhi->at(chosenEle))>2.5 ){
      counter[4]+=w;
      hetWithMetCut  -> Fill(et,w);
      hmetWithMetCut -> Fill(Met,w);
      hmtWithMetCut  -> Fill(mt,w);
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
