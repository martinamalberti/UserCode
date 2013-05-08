
///==== include ====

#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "readJSONFile.h"

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

#include "VertexAlgoParameters.h"
#include "HggVertexAnalyzer.h"
#include "HggVertexFromConversions.h"
#include "PhotonInfo.h"
#include "../src/selection.cc"
#include "../src/eleId95.cc"
#include "../src/BeamSpotReweighting.cc"
#include "../src/VertexIdEfficiencyReweighting.cc"

#include <iostream>

#include "TClonesArray.h"
#include "TMatrix.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

using namespace std;


void HZZ4LSelection( std::vector<ROOT::Math::XYZTVector>* leg, int& i11, int& i12, int &i21, int& i22){

  // select first pair of leptons (closest to MZ)
  float mindm = 99999999;
  float MZ = 91.188;
  float mass = 0 ;
  for (int i =0; i< leg->size(); i++){
    for (int j = i+1; j < leg->size(); j++){
      mass = (leg->at(i)+leg->at(j)).M();
      //cout << mass << "  " << i << "  " << j <<endl;
      if (fabs(mass-MZ) < mindm ){
	mindm = fabs(mass-MZ);
	i11=i;
	i12=j;
      }
    }
  }
  // select second pair --> two highest pt leptons
  float pt1 = -1.;
  float pt2 = -1.;
  for (int k = 0; k < leg->size(); k++){
    if (k == i11 ||  k == i12) continue;
    float pt = leg->at(k).Pt() ;
    //    cout << k << "  " << pt << endl; 
    if ( pt > pt1 ) {
      pt2 = pt1;
      i22 = i21;
      pt1 = pt;
      i21 = k;
    }
    else if (pt > pt2){
      pt2 = pt;
      i22 = k;
    }
  }
  

}



// --------- MAIN -------------------

int main(int argc, char** argv)
{ 
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>>>> VertexStudiesAnalysis::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
  
  // Parse the config file
  parseConfigFile (argv[1]) ;

  std::string inputFileList = gConfigParser -> readStringOption("Input::inputFileList");

  std::string treeName  = gConfigParser -> readStringOption("Input::treeName");
  
  std::string tmvaMethod    = gConfigParser -> readStringOption("Input::tmvaMethod");
  std::string tmvaWeights   = gConfigParser -> readStringOption("Input::tmvaWeights");
 
  std::string tmvaEventMethod    = gConfigParser -> readStringOption("Input::tmvaEventMethod");
  std::string tmvaEventWeights   = gConfigParser -> readStringOption("Input::tmvaEventWeights");
 
  std::string outputRootFilePath = gConfigParser -> readStringOption("Output::outputRootFilePath");
  std::string outputRootFileName = gConfigParser -> readStringOption("Output::outputRootFileName");  

  int entryMIN = gConfigParser -> readIntOption("Options::entryMIN");
  int entryMAX = gConfigParser -> readIntOption("Options::entryMAX");
 
  int isData   = gConfigParser -> readIntOption("Options::isData");

  double trackThr = gConfigParser -> readDoubleOption("Options::trackThr");

  int addConversionToMva = gConfigParser -> readIntOption("Options::addConversionToMva");
  int useMvaRanking      = gConfigParser -> readIntOption("Options::useMvaRanking");
  

  int useWeights      = gConfigParser -> readIntOption("Options::useWeights");
  int poissonWeights  = gConfigParser -> readIntOption("Options::poissonWeights");
  int nAvePU          = gConfigParser -> readIntOption("Options::nAvePU");
  
  std::string puweightsFileName = gConfigParser -> readStringOption("Options::puweightsFileName");  

  int useJSON      = gConfigParser -> readIntOption("Options::useJSON");
  std::string jsonFileName = gConfigParser -> readStringOption("Options::jsonFileName");  

  int doBSreweighting   = gConfigParser -> readIntOption("Options::doBSreweighting");
  float dzRightVertex   = gConfigParser -> readFloatOption("Options::dzRightVertex");

  //******* Get run/LS map from JSON file *******
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
  if ( isData && useJSON ) {
    std::cout << ">>> Getting GOOD  run/LS from JSON file" << std::endl;
    jsonMap = readJSONFile(jsonFileName);
    std::cout << std::endl;
    std::cout << std::endl;
  }

  
  //****** Get weights for MC ****** 
  float nmax;
  float w[100]={0.};
  TRandom *gRandom = new TRandom();
  
  if (useWeights){
    
    TFile weightsFile(puweightsFileName.c_str(),"READ");  
    TH1F* hweights;
     
    if ( poissonWeights ){
      std::cout << "N ave PU for Poisson PU reweighting : " << nAvePU << std::endl;
      char hname[100];
      sprintf(hname,"hwmc%d",nAvePU);
      hweights = (TH1F*)weightsFile.Get(hname);
    }
    else {
      hweights = (TH1F*)weightsFile.Get("hweights");
    }
    for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
      w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
    }
    hweights->GetXaxis()->SetRangeUser(5,60);
    nmax = hweights ->GetMaximum();
    std::cout << " Max weight " << nmax << std::endl;
    weightsFile.Close();
  }
  
  
  //****** Parameters for vertex finding algo ****** 
  VertexAlgoParameters vtxAlgoParams_;
  vtxAlgoParams_.rescaleTkPtByError = false;
  vtxAlgoParams_.trackCountThr      = 0.;
  vtxAlgoParams_.highPurityOnly     = false;
  vtxAlgoParams_.maxD0Signif        = 9999999.;
  vtxAlgoParams_.maxDzSignif        = 9999999.;
  vtxAlgoParams_.removeTracksInCone = 1;
  vtxAlgoParams_.coneSize           = 0.05;
  
  //--- sumpt2 ordering
  vector<string> ranksumpt2_;
  ranksumpt2_.push_back("logsumpt2");  

  //--- ranking product variables
  vector<string> rankVariables_;
  // variables order matters to resolve ties
  rankVariables_.push_back("ptbal"), rankVariables_.push_back("ptasym"),rankVariables_.push_back("logsumpt2");  
  
  //--- tmva variables
  vector<string> tmvaPerVtxVariables_;
  //original MVA 
  tmvaPerVtxVariables_.push_back("ptbal"), tmvaPerVtxVariables_.push_back("ptasym"),tmvaPerVtxVariables_.push_back("logsumpt2");  
  if( addConversionToMva ) {
    tmvaPerVtxVariables_.push_back("limPullToConv");
    tmvaPerVtxVariables_.push_back("nConv");
  } 
  
  //--- book per vertex TMVA
  TMVA::Reader *tmvaReader_ = new TMVA::Reader( "!Color:!Silent" );
  HggVertexAnalyzer::bookVariables( *tmvaReader_, tmvaPerVtxVariables_ );
  tmvaReader_->BookMVA( tmvaMethod.c_str(), tmvaWeights.c_str() );

  //--- book per event TMVA
  TMVA::Reader *tmvaPerEvtReader_ = new TMVA::Reader( "!Color:!Silent" );
  HggVertexAnalyzer::bookPerEventVariables( *tmvaPerEvtReader_, 3, true );
  tmvaPerEvtReader_->BookMVA( tmvaEventMethod.c_str(), tmvaEventWeights.c_str() );

  
  //****** BOOK OUTPUT HISTOGRAMS ******

  TH1F PtAll("PtAll","Pt of boson all",80,0,400);
  TH1F PtGood("PtGood","Pt of boson good",80,0,400);
  TH1F PtGood_BDT("PtGood_BDT","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK("PtGood_RANK","Pt of boson good (RANKING)",80,0,400);

  TH1F NvtAll("NvtAll","number of PV all",50,0,50);
  TH1F NvtGood("NvtGood","number of PV good",50,0,50);
  TH1F NvtGood_BDT("NvtGood_BDT","number of PV good (BDT)",50,0,50);
  TH1F NvtGood_RANK("NvtGood_RANK","number of PV good (RANKING)",50,0,50);

  TH1F NpuAll("NpuAll","number of PV all",50,0,50);
  TH1F NpuGood("NpuGood","number of PV good",50,0,50);
  TH1F NpuGood_BDT("NpuGood_BDT","number of PV good (BDT)",50,0,50);
  TH1F NpuGood_RANK("NpuGood_RANK","number of PV good (RANKING)",50,0,50);

  TH1F PtGood_matchedClosest("PtGood_matchedClosest","Pt of boson good",80,0,400);
  TH1F PtGood_BDT_matchedClosest("PtGood_BDT_matchedClosest","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK_matchedClosest("PtGood_RANK_matchedClosest","Pt of boson good (RANKING)",80,0,400);
  
  TH1F NvtGood_matchedClosest("NvtGood_matchedClosest","number of PV good",50,0,50);
  TH1F NvtGood_BDT_matchedClosest("NvtGood_BDT_matchedClosest","number of PV good (BDT)",50,0,50);
  TH1F NvtGood_RANK_matchedClosest("NvtGood_RANK_matchedClosest","number of PV good (RANKING)",50,0,50);

  TH1F NpuGood_matchedClosest("NpuGood_matchedClosest","number of PV good",50,0,50);
  TH1F NpuGood_BDT_matchedClosest("NpuGood_BDT_matchedClosest","number of PV good (BDT)",50,0,50);
  TH1F NpuGood_RANK_matchedClosest("NpuGood_RANK_matchedClosest","number of PV good (RANKING)",50,0,50);

  TH2F hdist("hdist"," hdist",80,0,200,400,-10,10);
  TH2F hdiff_dZ_muons("hdiff_dZ_muons","hdiff_dZ_muons",80,0,200,400,-10,10);
  TH2F hdiff_dZ_electrons("hdiff_dZ_electrons","hdiff_dZ_electrons",80,0,200,400,-10,10);

  TH1F BDToutput("BDToutput","BDT output",500,-1,1);
  TH1F BDToutput_sig("BDToutput_sig","BDT output - signal vertices",500,-1,1);
  TH1F BDToutput_bkg("BDToutput_bkg","BDT output - background",500,-1,1);
 
  TH1F perEventBDToutput("perEventBDToutput","BDT output",500,-1,1);
  TH1F perEventBDToutput_sig("perEventBDToutput_sig","BDT output - signal vertices",500,-1,1);
  TH1F perEventBDToutput_bkg("perEventBDToutput_bkg","BDT output - background",500,-1,1);

  TH1F ChosenVertex_BDT("ChosenVertex_BDT","index of chosen vertex (BDT)",50,0,50);
  TH1F ChosenVertexDz_BDT("ChosenVertexDz_BDT"," chosen vertex (z - z_{true}) (BDT)",100000,-50,50);
  TH2F ChosenVertexDz_BDT_vs_pt("ChosenVertexDz_BDT_vs_pt"," chosen vertex (z - z_{true}) (BDT)",200, 0, 200, 100000,-50,50);
  
  TH1F sumpt2("sumpt2","sumpt2",50,-5.,15.);
  TH1F sumpt2_sig("sumpt2_sig","sumpt2 - signal vertices",50,-5,15);
  TH1F sumpt2_bkg("sumpt2_bkg","sumpt2 - background",50,-5,15);

  TH1F ptbal("ptbal","ptbal",100,-50.,150.);
  TH1F ptbal_sig("ptbal_sig","ptbal - signal vertices",100,-50,150);
  TH1F ptbal_bkg("ptbal_bkg","ptbal - background",100,-50,150);
 
  TH1F ptasym("ptasym","ptasym",50,-1.,1.);
  TH1F ptasym_sig("ptasym_sig","ptasym - signal vertices",50,-1,1);
  TH1F ptasym_bkg("ptasym_bkg","ptasym - background",50,-1,1);

  TH1F hz("hz"," primary vertex z",10000,-50,50);

  TH1F pt2h("pt2h","pt2 H",500,0,500);
  TH1F pt2bkg("pt2bkg","pt2 bkg",500,0,500);

  float ww = 1;
  float mindz = dzRightVertex;
  float diff, bsweight; 

  //****** LOAD TREE ******
  TChain* chain = new TChain(treeName.c_str());
  FillChain(*chain, inputFileList.c_str());
  treeReader reader((TTree*)(chain));
  std::cout<<"found "<< reader.GetEntries() <<" entries"<<std::endl;

  //****** Start loop over entries ******
  int runId, lumiId;

  vector<int> runnum,lumisec,eventnum;

  for (int u = 0; u < reader.GetEntries(); u++ )
    {
      if(u == entryMAX) break;
      if(u < entryMIN)  continue;
      if(u%1 == 0) std::cout<<"reading event "<< u <<std::endl;
      reader.GetEntry(u);
      
      //*** filter bad runs/lumis
      runId = reader.GetInt("runId")->at(0);
      lumiId = reader.GetInt("lumiId")->at(0);
      
      bool skipEvent = false;
      if( isData && useJSON ){
	if(AcceptEventByRunAndLumiSection(runId,lumiId,jsonMap) == false)
	  skipEvent = true;
      }
      if( skipEvent == true ) continue;


      //*** pu weights
      std::vector<float>*PU_z ;
      std::vector<int>* mc_PUit_NumInteractions; 
      std::vector<float>* mc_PUit_TrueNumInteractions; 
      int npu ;
      float npuTrue;

      if ( !isData ){
	mc_PUit_NumInteractions  = reader.GetInt("mc_PUit_NumInteractions");
	npu = mc_PUit_NumInteractions->at(0);

      	mc_PUit_TrueNumInteractions  = reader.GetFloat("mc_PUit_TrueNumInteractions"); // needed for 2012 PU reweighting
	npuTrue = mc_PUit_TrueNumInteractions->at(0);

	//--- use weights 
	if (useWeights){
	  ww = w[int(npu)]; // observed
	}
      }

     
      //*** setup common branches ***
      std::vector<int>* PV_nTracks;
      std::vector<float>* PV_z;
      std::vector<float>* PV_z_muon;
      std::vector<float>* PV_d0;
      std::vector<ROOT::Math::XYZVector>* PVtracks;
      std::vector<int>* PVtracks_PVindex;
      std::vector<int>* tracks_PVindex;
      std::vector<float>* tracks_dz ; //?
      std::vector<float>* tracks_dz_PV ; //?
      std::vector<float>* tracks_dxy_PV ; //?
      
      float TrueVertex_Z;
      TLorentzVector Z1(0,0,0,0);
      TLorentzVector Z2(0,0,0,0);
      
      int i11=-1;
      int i12=-1;
      int i21=-1;
      int i22=-1;

      //*** selections for HZZ->4L ***
      std::vector<ROOT::Math::XYZTVector>* leg; // boson deacy products
      std::vector<ROOT::Math::XYZTVector>* muons = reader.Get4V("muons");
      std::vector<ROOT::Math::XYZTVector>* electrons = reader.Get4V("electrons");



      // count number of muons/ele with pt > 5 GeV
      int nMuons = 0;
      int nElectrons = 0;
      for (int i = 0 ; i < muons->size(); i++){
	if (muons->at(i).Pt()>5) nMuons++;
      }
      
      for (int i = 0 ; i < electrons->size(); i++){
	if (electrons->at(i).Pt()>5) nElectrons++;
      }

      cout << "number of muons     = " << nMuons << endl;
      cout << "number of electrons = " << nElectrons << endl;

      //      bool acceptEvent = muons->size()==4 || electrons->size()==4 || (muons->size()==2 && electrons->size()==2);
      bool is4e    = nElectrons >=4 && nMuons <= 2;
      bool is4mu   = nMuons >=4 && nElectrons <= 2;
      bool is2e2mu = nMuons == 2 && nElectrons == 2;
      bool acceptEvent = is4e || is4mu || is2e2mu;
      if (!acceptEvent) continue;
      
      if (is4e){
	std::vector<float>* electrons_dz_PV_noEle = reader.GetFloat("electrons_dz_PV_noEle");
	
	PV_nTracks       = reader.GetInt("PV_noEle_nTracks");
	PV_z             = reader.GetFloat("PV_noEle_z");
	PV_d0            = reader.GetFloat("PV_noEle_d0");
	PVtracks         = reader.Get3V("PVEleLessTracks");
	PVtracks_PVindex = reader.GetInt("PVEleLessTracks_PVindex");
	tracks_PVindex   = reader.GetInt("tracks_PVindex");
	tracks_dxy_PV    = reader.GetFloat("tracks_dxy_PV");
	tracks_dz_PV     = reader.GetFloat("tracks_dz_PV");
	tracks_dz        = reader.GetFloat("tracks_dz");
	leg              = reader.Get4V("electrons");

	//cout << " choosing Z1 and Z2" <<endl;
	HZZ4LSelection(leg,i11,i12,i21,i22);

	cout << i11 << "  " << i12 << "  " <<  i21 << "  "<< i22<<endl;

	//cout << " setting Z1 and Z2" <<endl;
	ROOT::Math::XYZTVector sum1 = leg->at(i11)+leg->at(i12); 
	ROOT::Math::XYZTVector sum2 = leg->at(i21)+leg->at(i22); 
	Z1.SetPxPyPzE(sum1.Px(), sum1.Py(), sum1.Pz(), sum1.E());
	Z2.SetPxPyPzE(sum2.Px(), sum2.Py(), sum2.Pz(), sum2.E());

	TrueVertex_Z = PV_z->at(0) + (electrons_dz_PV_noEle->at(i11) + electrons_dz_PV_noEle->at(i12) + electrons_dz_PV_noEle->at(i21) + electrons_dz_PV_noEle->at(i22))/4.;
      }

      if ( is4mu ){
	std::vector<float>* muons_dz_PV_noMuon = reader.GetFloat("muons_dz_PV_noMuon");
		
	PV_nTracks       = reader.GetInt("PV_noMuon_nTracks");
	PV_z             = reader.GetFloat("PV_noMuon_z");
	PV_d0            = reader.GetFloat("PV_noMuon_d0");
	PVtracks         = reader.Get3V("PVMuonLessTracks");
	PVtracks_PVindex = reader.GetInt("PVMuonLessTracks_PVindex");
	PV_z_muon        = reader.GetFloat("PV_z");
	tracks_PVindex   = reader.GetInt("tracks_PVindex");
	tracks_dxy_PV    = reader.GetFloat("tracks_dxy_PV");
	tracks_dz_PV     = reader.GetFloat("tracks_dz_PV");
	tracks_dz        = reader.GetFloat("tracks_dz");
	leg              = reader.Get4V("muons");  // use muon info for SC
	
	HZZ4LSelection(leg,i11,i12,i21,i22);
	ROOT::Math::XYZTVector sum1 = leg->at(i11)+leg->at(i12); 
	ROOT::Math::XYZTVector sum2 = leg->at(i21)+leg->at(i22); 
	Z1.SetPxPyPzE(sum1.Px(), sum1.Py(), sum1.Pz(), sum1.E());
	Z2.SetPxPyPzE(sum2.Px(), sum2.Py(), sum2.Pz(), sum2.E());
	TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(i11) + muons_dz_PV_noMuon->at(i12) + muons_dz_PV_noMuon->at(i21) + muons_dz_PV_noMuon->at(i22))/4.;
	
	//if ( fabs(muons_dz_PV_noMuon->at(indpho1) - muons_dz_PV_noMuon->at(indpho2))> 0.5) continue;
      }//4mu end


      if ( is2e2mu ){
	std::vector<float>* muons_dz_PV_noMuon = reader.GetFloat("muons_dz_PV_noMuon");
		
	PV_nTracks       = reader.GetInt("PV_noMuon_nTracks");
	PV_z             = reader.GetFloat("PV_noMuon_z");
	PV_d0            = reader.GetFloat("PV_noMuon_d0");
	PVtracks         = reader.Get3V("PVMuonLessTracks");
	PVtracks_PVindex = reader.GetInt("PVMuonLessTracks_PVindex");
	PV_z_muon        = reader.GetFloat("PV_z");
	tracks_PVindex   = reader.GetInt("tracks_PVindex");
	tracks_dxy_PV    = reader.GetFloat("tracks_dxy_PV");
	tracks_dz_PV     = reader.GetFloat("tracks_dz_PV");
	tracks_dz        = reader.GetFloat("tracks_dz");
	leg              = reader.Get4V("muons");  // use muon info for SC
	
	HZZ4LSelection(leg,i11,i12,i21,i22);
	ROOT::Math::XYZTVector sum1 = leg->at(i11)+leg->at(i12); 
	ROOT::Math::XYZTVector sum2 = leg->at(i21)+leg->at(i22); 
	Z1.SetPxPyPzE(sum1.Px(), sum1.Py(), sum1.Pz(), sum1.E());
	Z2.SetPxPyPzE(sum2.Px(), sum2.Py(), sum2.Pz(), sum2.E());
	TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(i11) + muons_dz_PV_noMuon->at(i12) + muons_dz_PV_noMuon->at(i21) + muons_dz_PV_noMuon->at(i22))/4.;
	
	//if ( fabs(muons_dz_PV_noMuon->at(indpho1) - muons_dz_PV_noMuon->at(indpho2))> 0.5) continue;
      }//4mu end

    

      // branches buffers
      int nvtx_;
      float  vtxx_[1000], vtxy_[1000], vtxz_[1000];
      int ntracks_;
      float tkpx_[5000], tkpy_[5000], tkpz_[5000], tkPtErr_[5000], tkWeight_[5000], tkd0_[5000], tkd0Err_[5000], tkdz_[5000], tkdzErr_[5000];
      int tkVtxId_[5000];
      bool tkIsHighPurity_[5000];
      
      float phocalox_[100], phocaloy_[100], phocaloz_[100], phoen_[100];
      
      // set variables
      
      // vertices 
      nvtx_    = (int) PV_z->size();
      for ( int iv = 0; iv < nvtx_; iv++){
	vtxx_[iv] =  0;
	vtxy_[iv] =  0;
	vtxz_[iv] =  PV_z->at(iv) ;
      }
      
      // tracks
      ntracks_ = PVtracks->size();
      for (int itrk = 0; itrk <ntracks_; itrk++ ){
	tkpx_[itrk]    = PVtracks->at(itrk).X();
	tkpy_[itrk]    = PVtracks->at(itrk).Y();
	tkpz_[itrk]    = PVtracks->at(itrk).Z();
	tkPtErr_[itrk] = 0;
	tkVtxId_[itrk] = PVtracks_PVindex->at(itrk);
	
	tkWeight_[itrk]= 1.;
	tkd0_[itrk]    = 0;
	tkd0Err_[itrk] = 0;
	tkdz_[itrk]    = 0;
	tkdzErr_[itrk] = 0;
	tkIsHighPurity_[itrk]= 1.;
      }
      
      
      //*** set vertex info
      TupleVertexInfo vinfo( nvtx_, vtxx_ , vtxy_, vtxz_, ntracks_, tkpx_, tkpy_, tkpz_, tkPtErr_, tkVtxId_, tkWeight_, tkd0_, tkd0Err_,tkdz_, tkdzErr_ , tkIsHighPurity_);
           
      //*** set photon info
      PhotonInfo pho1(0, Z1.Vect(),Z1.E());
      PhotonInfo pho2(1, Z2.Vect(),Z2.E());

      TLorentzVector Higgs = Z1 + Z2;
      

      //*** vertex analyzer
      HggVertexAnalyzer vAna(vtxAlgoParams_,nvtx_);
      vAna.setNConv(0);
      vAna.analyze(vinfo,pho1,pho2);
            
      //*** preselect vertices 
      std::vector<int> presel;
      for(int i=0; i<nvtx_; i++) {
	presel.push_back(i); 
	hz.Fill(PV_z->at(i));
      }
      vAna.preselection(presel);
      
      
      //*** Look if the H vertex matches one of the PV vertices
      float dmin = 10000;
      int iClosest = -1;
      for ( int uu = 0; uu < nvtx_; uu++){
	float distt = fabs( PV_z->at(uu) - TrueVertex_Z );
	if ( distt < dmin)   { dmin = distt; iClosest = uu; }
      }
      //*** NOW FILL HISTOGRAMS

      //*** BDT 
      vector<int> ranktmva = vAna.rank(*tmvaReader_,tmvaMethod);
      
      // -- BS reweighting 
      if (doBSreweighting){
        diff = PV_z->at(ranktmva[0])-TrueVertex_Z;
	bsweight = BSweight(diff);
	//cout << diff << "  "  << bsweight << endl;
	ww*=bsweight;
      }

      ChosenVertex_BDT.Fill(ranktmva[0],ww);
      ChosenVertexDz_BDT.Fill(TrueVertex_Z - PV_z->at(ranktmva[0]),ww);
      ChosenVertexDz_BDT_vs_pt.Fill(Higgs.Pt(),TrueVertex_Z - PV_z->at(ranktmva[0]),ww);
      float vtxmva, evtmva;

      // fill per vertex mva 
      for (int iv=0;iv<ranktmva.size();iv++) {
	float vtxmva = vAna.mva(ranktmva[iv]);
	if ( iClosest == ranktmva[iv] ) {
	  BDToutput_sig.Fill( vtxmva, ww );
	  sumpt2_sig.Fill(vAna.logsumpt2(ranktmva[iv]),ww ) ;
	  ptasym_sig.Fill(vAna.ptasym(ranktmva[iv]),ww );
	  ptbal_sig .Fill(vAna.ptbal(ranktmva[iv]),ww );
	}
	else {
	  BDToutput_bkg.Fill( vtxmva, ww );
	  sumpt2_bkg.Fill(vAna.logsumpt2(ranktmva[iv]),ww ) ;
	  ptasym_bkg.Fill(vAna.ptasym(ranktmva[iv]),ww );
	  ptbal_bkg .Fill(vAna.ptbal(ranktmva[iv]),ww );
	}
      }

      // fill per event mva
      evtmva = vAna.perEventMva(*tmvaPerEvtReader_, tmvaEventMethod.c_str(),ranktmva);
      perEventBDToutput.Fill( evtmva, ww);
      if (fabs( TrueVertex_Z - PV_z->at(ranktmva[0]) ) < mindz) {
      	perEventBDToutput_sig.Fill( evtmva, ww );
      }
      else{
	perEventBDToutput_bkg.Fill( evtmva, ww );
      }

      //-- matching 1cm
      if ( fabs( TrueVertex_Z - PV_z->at(ranktmva[0]) ) < mindz) {
	PtGood_BDT.Fill( Higgs.Pt(),ww );
       	NvtGood_BDT.Fill( nvtx_ ,ww);
	if (!isData) NpuGood_BDT.Fill(npu,ww);
      }
      //-- matching closest vtx
      if ( iClosest == ranktmva[0]){
	PtGood_BDT_matchedClosest.Fill( Higgs.Pt(),ww );
	NvtGood_BDT_matchedClosest.Fill( nvtx_,ww );
	if (!isData) NpuGood_BDT_matchedClosest.Fill(npu,ww);
      }
      
      // -- all vtx
      PtAll.Fill( Higgs.Pt(),ww );
      NvtAll.Fill( nvtx_,ww );
      if (!isData) NpuAll.Fill(npu,ww);
      
      
    }// end loop over entries

  std::cout << "END LOOP OVER ENTRIES" << std::endl;
  std::cout << "Saving histos on file ..." << std::endl;
  
  TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");

  //hAcceptedLumis -> Write();
  
  PtAll.Write();
  PtGood.Write();
  PtGood_BDT.Write();
  PtGood_RANK.Write();
    
  NvtAll.Write();
  NvtGood.Write();
  NvtGood_BDT.Write(); 
  NvtGood_RANK.Write();

  NpuAll.Write();
  NpuGood.Write();
  NpuGood_BDT.Write(); 
  NpuGood_RANK.Write();
    
  PtGood_matchedClosest.Write();
  PtGood_BDT_matchedClosest.Write();
  PtGood_RANK_matchedClosest.Write();
    
  NvtGood_matchedClosest.Write();
  NvtGood_BDT_matchedClosest.Write(); 
  NvtGood_RANK_matchedClosest.Write();

  NpuGood_matchedClosest.Write();
  NpuGood_BDT_matchedClosest.Write(); 
  NpuGood_RANK_matchedClosest.Write();

  hdist.Write();
  hdiff_dZ_muons.Write();
  hdiff_dZ_electrons.Write();
  
  BDToutput.Write();
  BDToutput_sig.Write();
  BDToutput_bkg.Write();
  perEventBDToutput.Write();
  perEventBDToutput_sig.Write();
  perEventBDToutput_bkg.Write();

  ChosenVertex_BDT.Write();
  ChosenVertexDz_BDT.Write();
  ChosenVertexDz_BDT_vs_pt.Write();
  hz.Write();

  sumpt2_sig.Write();
  ptbal_sig.Write();
  ptasym_sig.Write();
  sumpt2_bkg.Write();
  ptbal_bkg.Write();
  ptasym_bkg.Write();

  ff.Close();
  
  std::cout << "BYE BYE !!!! " << std::endl;

  return 0; 
  

}











