#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

//---- CMSSW includes
#include "DataFormats/Math/interface/LorentzVector.h"
#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeTreeContent.h"

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
 

  // Reading input parameters
  
  //---- Input Variables ----
  std::string FileIn;
  std::string FileOut;
  float Xsec;
  float Lumi;
  float GenFilterEff;
  int NeventsDataset;
  float MaxJetPt = 100;

  std::string inputName = argv[1];
  ifstream file(inputName.c_str());

  std::string variableNameFileIn    = "FileIn";
  std::string variableNameFileOut   = "FileOut";
  std::string variableNameXsec      = "Xsec";
  std::string variableNameLumi      = "Lumi";
  std::string variableNameFilterEff = "GenFilterEff";
  std::string variableNameNevents   = "NeventsDataset";
  std::string variableNameMaxJetPt  = "MaxJetPt";
  std::cerr << " Reading  " << inputName << " ... " << std::endl;

  file >>  variableNameFileIn    >> FileIn
       >>  variableNameFileOut   >> FileOut
       >>  variableNameXsec      >> Xsec
       >>  variableNameLumi      >> Lumi
       >>  variableNameFilterEff >> GenFilterEff
       >>  variableNameNevents   >> NeventsDataset
       >>  variableNameMaxJetPt  >> MaxJetPt;


  cout <<  variableNameFileIn    << "  "  << FileIn          << "  "
       <<  variableNameFileOut   << "  "  << FileOut         << "  "
       <<  variableNameXsec      << "  "  << Xsec            << "  "
       <<  variableNameLumi      << "  "  << Lumi            << "  "
       <<  variableNameFilterEff << "  "  << GenFilterEff    << "  "
       <<  variableNameNevents   << "  "  << NeventsDataset  << "  "
       <<  variableNameMaxJetPt  << "  "  << MaxJetPt
       <<  endl;

  
  // load the tree
  TChain *chain = new TChain ("myanalysis/WprimeAnalysisTree") ;
  WprimeTreeContent treeVars ;
  setBranchAddresses (chain, treeVars) ;
  
  chain->Add(FileIn.c_str());
  int nEntries = chain->GetEntries () ;
  std::cout << "FOUND " << nEntries << " ENTRIES\n" ;
  
  double w = Xsec*GenFilterEff*Lumi/(float)nEntries;
  cout << "Events weight = " << w << endl;

  // FILE to save histos
  TFile *fout = new TFile(FileOut.c_str(),"recreate");

  //INITIALIZING HISTOGRAMS
  TH1F *hNGoodElectrons = new TH1F("hNGoodElectrons","hNGoodElectrons",5,0,5);
  TH1F *het = new TH1F("het","het",1500,0,1500);
  TH1F *hmet = new TH1F("hmet","hmet",1500,0,1500);
  TH1F *hmetUncorr = new TH1F("hmetUncorr","hmetUncorr",1500,0,1500);
  TH1F *hmt = new TH1F("hmt","hmt",3000,0,3000);

  TH1F *hetWithJetVeto  = new TH1F("hetWithJetVeto","hetWithJetVeto",1500,0,1500);
  TH1F *hmetWithJetVeto = new TH1F("hmetWithJetVeto","hmetWithJetVeto",1500,0,1500);
  TH1F *hmtWithJetVeto  = new TH1F("hmtWithJetVeto","hmtWithJetVeto",3000,0,3000);
  
  TH1F *hetWithMetCut  = new TH1F("hetWithMetCut","hetWithMetCut",1500,0,1500);
  TH1F *hmetWithMetCut = new TH1F("hmetWithMetCut","hmetWithMetCut",1500,0,1500);
  TH1F *hmtWithMetCut  = new TH1F("hmtWithMetCut","hmtWithMetCut",3000,0,3000);
  
  TH1F *hjetPt      = new TH1F("hjetPt","hjetPt",510,-10,500);
  TH1F *hjetEta     = new TH1F("hjetEta","hjetEta",240,-6,6);
  TH1F *hNjets      = new TH1F("hNjets","hNjets",50,0,50);
  TH1F *hDphiJetEle = new TH1F("hDphiJetEle","hDphiJetEle",400,0,4);
  TH2F *h2JetEle    = new TH2F("h2JetEle","h2JetEle",400,0,4,500,0,500);

  TH1F *hjetPt_HighEtEle      = new TH1F("hjetPt_HighEtEle","hjetPt_HighEtEle",510,-10,500);
  TH1F *hjetEta_HighEtEle     = new TH1F("hjetEta_HighEtEle","hjetEta_HighEtEle",240,-6,6);
  TH1F *hNjets_HighEtEle      = new TH1F("hNjets_HighEtEle","hNjets_HighEtEle",50,0,50);
  TH1F *hDphiJetEle_HighEtEle = new TH1F("hDphiJetEle_HighEtEle","hDphiJetEle_HighEtEle",400,0,4);
  TH2F *h2JetEle_HighEtEle    = new TH2F("h2JetEle_HighEtEle","h2JetEle_HighEtEle",400,0,4,500,0,500);
  
  TH1F *hEtOverMet            = new TH1F("hEtOverMet","hEtOverMet",1000,0,10);
  TH1F *hDphiMetEle           = new TH1F("hDphiMetEle","hDphiMetEle",400,0,4);

  TH1F *hEtOverMet_HighEtEle  = new TH1F("hEtOverMet_HighEtEle","hEtOverMet_HighEtEle",1000,0,10);
  TH1F *hDphiMetEle_HighEtEle = new TH1F("hDphiMetEle_HightEtEle","hDphiMetEle_HighEtEle",400,0,4);
  
  float EtaCutEB    = 1.4442;
  float EtaCutEE    = 1.560;
  float combIsoCut  = 0.;
  float combinedIso = 0.;

  double counter[10] = {0,0,0,0,0,0,0,0,0,0};

  // loop over entries
  for (int entry = 0; entry < nEntries; ++entry) {
    
    chain->GetEntry (entry) ;
    if(entry%100000 == 0) std::cout << "event n. " << entry << std::endl;
    
    counter[0] = counter[0] + w;
    
    // HLT selection
    if (treeVars.HLT_Photon15_L1R==0) continue;
    counter[1] = counter[1] + w;
    
    
    int nGoodElectrons = 0;
    int chosenEle = 0;

    //start loop over electron candidates
    for (int i=0; i<treeVars.nElectrons; i++){
      
      float et  = treeVars.eleEt[i] ;
      float eta = treeVars.eleEta[i] ;
      
      // keep only electrons in ECAL fiducial volume
      if ( fabs(eta) > EtaCutEB && fabs(eta)  < EtaCutEE) continue;

      // keep only electrons with ET > 30 GeV
      
      if ( et < 30) continue ;
      
       if ( treeVars.eleId[i] != 0x0 ) continue;
      
      /*
      // electron ID
      if ( treeVars.eleId[i] == 0 ) continue;
      
      // isolation
      if  ( fabs(eta) < EtaCutEB &&  treeVars.eleTrkIso[i] > 7.5)  continue;
      if  ( fabs(eta) > EtaCutEE &&  treeVars.eleTrkIso[i] > 15.)  continue;
      
      combinedIso = treeVars.eleEcalIso[i] + treeVars.eleHcalIsoD1[i] ;
      if ( fabs(eta) < EtaCutEB ) combIsoCut = 2 + 0.03*et;
      if ( fabs(eta) > EtaCutEE && et < 50.) combIsoCut = 2.5;
      if ( fabs(eta) > EtaCutEE && et > 50.) combIsoCut = 2.5 + 0.03*(et-50);
      if ( combinedIso > combIsoCut ) continue;
      
      if ( treeVars.eleHcalIsoD2[i] > 0.5 && fabs(eta) > EtaCutEE) continue;
      */
      nGoodElectrons++;
      chosenEle = i;

    }// end loop over ele cand
  


    hNGoodElectrons->Fill( nGoodElectrons,w);
    
    if ( nGoodElectrons != 1 ) continue;
    counter[2] = counter[2] + w;
    
    double et = treeVars.eleEt[chosenEle];
    double pt = sqrt(treeVars.elePx[chosenEle]*treeVars.elePx[chosenEle]+
		     treeVars.elePy[chosenEle]*treeVars.elePy[chosenEle]);

    double met = treeVars.Met;
    
    double cphi = (treeVars.elePx[chosenEle]*treeVars.Mex + treeVars.elePy[chosenEle]*treeVars.Mey )
                  /(met*pt);
    double mt = sqrt(2*et*met*(1-cphi));
    

    // no selections
    het  -> Fill(et,w);
    hmt  -> Fill(mt,w);
    hmet -> Fill(met,w);
    hmetUncorr -> Fill(treeVars.uncorrMet,w);
  
  
    //jet variables: considero solo i jet entro |eta| < 3.
    double maxPt = -1.;
    int jmax = -1;
    int njet = 0 ;
    for (int j=0; j<treeVars.nJets; j++){
      double deta = treeVars.jetEta[j] -  treeVars.eleEta[chosenEle];
      double dphi = treeVars.jetPhi[j] -  treeVars.elePhi[chosenEle];
      double dR = sqrt(deta*deta+dphi*dphi);
      if (dR < 0.3) continue;
      if ( fabs( treeVars.jetEta[j] )<3 && treeVars.jetPt[j] > maxPt) {
	jmax  = j ;
        maxPt = treeVars.jetPt[j] ;
      }
      
      if ( treeVars.jetPt[j] > 15 && fabs( treeVars.jetEta[j] )<3 ) njet++;
    }
    
    hNjets -> Fill(njet,w);
    hjetPt -> Fill(maxPt,w); // riempio anche se non ci sono  jet
    
    if (jmax!=-1){
      hjetEta     -> Fill(treeVars.jetEta[jmax],w);
      hDphiJetEle -> Fill(deltaPhi(treeVars.jetPhi[jmax],treeVars.elePhi[chosenEle]),w);
      h2JetEle    -> Fill(deltaPhi(treeVars.jetPhi[jmax],treeVars.elePhi[chosenEle]),treeVars.jetPt[jmax],w);
    }
    
    hEtOverMet  -> Fill(et/met,w);
    hDphiMetEle -> Fill(deltaPhi(treeVars.MetPhi,treeVars.elePhi[chosenEle]),w);
    
    if (treeVars.eleEt[chosenEle] > 100 ){

      hNjets_HighEtEle -> Fill(njet,w);
      
      // guardo il pT dei jet opposti all'elettrone
      if ( jmax==-1 || (jmax!=-1 && deltaPhi(treeVars.jetPhi[jmax],treeVars.elePhi[chosenEle]) < 1.57) ) 
	hjetPt_HighEtEle -> Fill(-1.,w);
      
      if (jmax!=-1){
	if ( deltaPhi(treeVars.jetPhi[jmax],treeVars.elePhi[chosenEle]) > 1.57){    
	  hjetPt_HighEtEle -> Fill(maxPt,w);
	  hjetEta_HighEtEle-> Fill(treeVars.jetEta[jmax],w);
	}
	hDphiJetEle_HighEtEle -> Fill(deltaPhi(treeVars.jetPhi[jmax],treeVars.elePhi[chosenEle]),w);
	h2JetEle_HighEtEle    -> Fill(deltaPhi(treeVars.jetPhi[jmax],treeVars.elePhi[chosenEle]),
				      treeVars.jetPt[jmax],w);
      }
      
      hEtOverMet_HighEtEle  -> Fill(et/met,w);
      hDphiMetEle_HighEtEle -> Fill(deltaPhi(treeVars.MetPhi,treeVars.elePhi[chosenEle]),w);
    }
    


    //met based selections
    if (et/met>0.4 && et/met<1.5 && deltaPhi(treeVars.MetPhi,treeVars.elePhi[chosenEle])>2.5 ){
      counter[4] = counter[4]+ w;
      hetWithMetCut  -> Fill(et,w);
      hmetWithMetCut -> Fill(met,w);
      hmtWithMetCut  -> Fill(mt,w);
    }

    //jet based selections
    if (jmax!=-1 && (maxPt > MaxJetPt) && deltaPhi(treeVars.jetPhi[jmax],treeVars.elePhi[chosenEle])>1.57 ) continue;
    {
      counter[3] =counter[3]+ w;
      hetWithJetVeto  -> Fill(et,w);
      hmetWithJetVeto -> Fill(met,w);
      hmtWithJetVeto  -> Fill(mt,w);
    }
    
    
  } // end loop over entries
  
  int bin200 = het ->FindBin(200);
  int bin300 = het ->FindBin(300);
  int bin2 = het ->FindBin(1500);
  float nHighEt200 = het->Integral(bin200,bin2);
  float nHighEt300 = het->Integral(bin300,bin2);
  
  float nn = Xsec*Lumi;

  cout << "Total Number of Events in the dataset: " << nEntries << endl;
  cout << "Selection                      # events (L = " << Lumi << " pb^-1)"<<      "    Efficiency " <<endl;
  cout << "No selections                  " << counter[0] << "      " << 1.  <<endl;
  cout << "HLT Photon15                   " << counter[1] << "      " << float(counter[1])/float(counter[0]) <<endl;
  cout << "Electron ET > 30, eleId+iso.   " << counter[2] << "      " << float(counter[2])/float(counter[0]) <<endl;
  cout << endl;
  
  cout << "w/o Jet Veto                   " << counter[2] << "      " << float(counter[2])/float(counter[0]) <<endl;
  cout << "w/o Jet Veto, ET>200           " << nHighEt200 << "      " << nHighEt200/float(counter[0])  <<endl;
  cout << "w/o Jet Veto, ET>300           " << nHighEt300 << "      " << nHighEt300/float(counter[0])  <<endl;
  cout << endl;
 
  nHighEt200 = hetWithJetVeto->Integral(bin200,bin2);
  nHighEt300 = hetWithJetVeto->Integral(bin300,bin2);
 
  cout << "w/ Jet Veto                    " << counter[3] << "      " << float(counter[3])/float(counter[0]) <<endl;
  cout << "w/ Jet Veto, ET>200            " << nHighEt200 << "      " << nHighEt200/float(counter[0])  <<endl;
  cout << "w/ Jet Veto, ET>300            " << nHighEt300 << "      " << nHighEt300/float(counter[0])  <<endl;

 
  cout << endl;


  
  // write histos
  fout->Write();
  fout->Close();

  
}
