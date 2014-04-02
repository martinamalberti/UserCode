///==== include ====

#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "readJSONFile.h"

#include "TF1.h"
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



#define R_ECAL    129
#define Z_ENDCAP  317

#define etaEB   1.4442 
#define etaEE   1.566


using namespace std;




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
 
  int isZee    = gConfigParser -> readIntOption("Options::isZee");
  int isZmumu  = gConfigParser -> readIntOption("Options::isZmumu");
  int isHiggs  = gConfigParser -> readIntOption("Options::isHiggs");

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

  int applyVertexIdScaleFactor = gConfigParser -> readIntOption("Options::applyVertexIdScaleFactor");

//   int nVertexMin = gConfigParser -> readIntOption("Options::nVertexMin");
//   int nVertexMax = gConfigParser -> readIntOption("Options::nVertexMax");

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
  // my TMVA
  //tmvaPerVtxVariables_.push_back("logsumpt2"), tmvaPerVtxVariables_.push_back("ptbal"),tmvaPerVtxVariables_.push_back("ptasym");
  //my TMVA with 4 vars
  //tmvaPerVtxVariables_.push_back("logsumpt2"), tmvaPerVtxVariables_.push_back("ptbal"),tmvaPerVtxVariables_.push_back("ptasym"),tmvaPerVtxVariables_.push_back("nch");  
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

  // vertex probability formula
  TF1 *fprob = new TF1("fprob","1.-0.49*(x+1)",-2,2);

  
  //****** BOOK OUTPUT HISTOGRAMS ******

  TH1F PtAll("PtAll","Pt of boson all",80,0,400);
  TH1F PtGood("PtGood","Pt of boson good",80,0,400);
  TH1F PtGood_BDT("PtGood_BDT","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK("PtGood_RANK","Pt of boson good (RANKING)",80,0,400);

  TH1F PtAll_EBEB("PtAll_EBEB","Pt of boson all",80,0,400);
  TH1F PtGood_EBEB("PtGood_EBEB","Pt of boson good",80,0,400);
  TH1F PtGood_BDT_EBEB("PtGood_BDT_EBEB","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK_EBEB("PtGood_RANK_EBEB","Pt of boson good (RANKING)",80,0,400);

  TH1F PtAll_EBEE("PtAll_EBEE","Pt of boson all",80,0,400);
  TH1F PtGood_EBEE("PtGood_EBEE","Pt of boson good",80,0,400);
  TH1F PtGood_BDT_EBEE("PtGood_BDT_EBEE","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK_EBEE("PtGood_RANK_EBEE","Pt of boson good (RANKING)",80,0,400);

  TH1F PtAll_EEEE("PtAll_EEEE","Pt of boson all",80,0,400);
  TH1F PtGood_EEEE("PtGood_EEEE","Pt of boson good",80,0,400);
  TH1F PtGood_BDT_EEEE("PtGood_BDT_EEEE","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK_EEEE("PtGood_RANK_EEEE","Pt of boson good (RANKING)",80,0,400);

  TH1F EtaAll("EtaAll","Eta of max SC",50,-5,5);
  TH1F EtaGood("EtaGood","Eta of max SC good",50,-5,5);
  TH1F EtaGood_BDT("EtaGood_BDT","Eta of max SC good (BDT)",50,-5,5);
  TH1F EtaGood_RANK("EtaGood_RANK","Eta of max SC good (RANKING)",50,-5,5);
  
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

  TH1F EtaGood_matchedClosest("EtaGood_matchedClosest","Eta of max SC good",50,-5,5);
  TH1F EtaGood_BDT_matchedClosest("EtaGood_BDT_matchedClosest","Eta of max SC good (BDT)",50,-5,5);
  TH1F EtaGood_RANK_matchedClosest("EtaGood_RANK_matchedClosest","Eta of max SC good (RANKING)",50,-5,5);
  
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

  TH1F vtxProbability("vtxProbability","vertex probability",250,0,1);
  TH1F vtxProbability_sig("vtxProbability_sig","vertex probability- signal vertices",250,0,1);
  TH1F vtxProbability_bkg("vtxProbability_bkg","vertex probability - background",250,0,1);

  TH1F *perVertexBDToutput_vtxbin[3];
  TH1F *perVertexBDToutput_sig_vtxbin[3];
  TH1F *perVertexBDToutput_bkg_vtxbin[3];
  TH1F *perEventBDToutput_vtxbin[3];
  TH1F *perEventBDToutput_sig_vtxbin[3];
  TH1F *perEventBDToutput_bkg_vtxbin[3];
  for (int vtxbin = 0; vtxbin < 3; vtxbin++ ){
    perVertexBDToutput_vtxbin[vtxbin] = new TH1F(Form("perVertexBDToutput_vtxbin%d",vtxbin), Form("perVertexBDToutput_vtxbin%d",vtxbin),500,-1,1);
    perVertexBDToutput_sig_vtxbin[vtxbin] = new TH1F(Form("perVertexBDToutput_sig_vtxbin%d",vtxbin), Form("perVertexBDToutput_sig_vtxbin%d",vtxbin),500,-1,1);
    perVertexBDToutput_bkg_vtxbin[vtxbin] = new TH1F(Form("perVertexBDToutput_bkg_vtxbin%d",vtxbin), Form("perVertexBDToutput_bkg_vtxbin%d",vtxbin),500,-1,1);
    perEventBDToutput_vtxbin[vtxbin] = new TH1F(Form("perEventBDToutput_vtxbin%d",vtxbin), Form("perEventBDToutput_vtxbin%d",vtxbin),500,-1,1);
    perEventBDToutput_sig_vtxbin[vtxbin] = new TH1F(Form("perEventBDToutput_sig_vtxbin%d",vtxbin), Form("perEventBDToutput_sig_vtxbin%d",vtxbin),500,-1,1);
    perEventBDToutput_bkg_vtxbin[vtxbin] = new TH1F(Form("perEventBDToutput_bkg_vtxbin%d",vtxbin), Form("perEventBDToutput_bkg_vtxbin%d",vtxbin),500,-1,1);
  }
 
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


  //TH2F * hAcceptedLumis = new TH2F("hAcceptedLumis","hAcceptedLumis",20000, 160000, 180000, 10000, 0, 10000);

  VertexIdEfficiencyReweighting *myVtxIdEffRw = new VertexIdEfficiencyReweighting("/afs/cern.ch/work/m/malberti/HGG/LegacyPaper/CMSSW_5_3_10/src/HIGGS/VertexStudies/macros/vtxIdScaleFactorFromZmumu_2012_legacypaper_v0.root","hscaleFactor") ;

  float ww = 1;
  float r9cut = 0.93;
  float mindz = dzRightVertex;
  float diff, bsweight; 

  //****** LOAD TREE ******
  TChain* chain = new TChain(treeName.c_str());
  FillChain(*chain, inputFileList.c_str());
  treeReader reader((TTree*)(chain));
  std::cout<<"found "<< reader.GetEntries() <<" entries"<<std::endl;

  //****** Start loop over entries ******
  int runId, lumiId;

  for (int u = 0; u < reader.GetEntries(); u++ )
    {
      if(u == entryMAX) break;
      if(u < entryMIN)  continue;
      if(u%10000 == 0) std::cout<<"reading event "<< u <<std::endl;
      reader.GetEntry(u);
      
      //*** filter bad runs/lumis
      runId = reader.GetInt("runId")->at(0);
      lumiId = reader.GetInt("lumiId")->at(0);
      
      bool skipEvent = false;
      if( isData && useJSON ){
	if(AcceptEventByRunAndLumiSection(runId,lumiId,jsonMap) == false)
	  skipEvent = true;
	if ( runId==190949 || runId==191090 || runId==191112 || runId==191116 )  // skip low lumi runs
	  skipEvent = true;
      }
      if( skipEvent == true ) continue;
      //hAcceptedLumis -> Fill(runId, lumiId);

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

	//	if ( npuTrue<6 ) continue;  // skip low lumi runs

	//--- use weights 
	if (useWeights){
	  //float myrnd = gRandom->Uniform(0,nmax);
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
      std::vector<ROOT::Math::XYZTVector>* sc; // supercluster
      
      int accept = 0;
      int indpho1 = -100;
      int indpho2 = -100;
	
      ROOT::Math::XYZTVector sum2pho;
      float etaMaxSC ;
      float TrueVertex_Z;

      //*** selections for Hgg ***
      if (isHiggs){
	std::vector<ROOT::Math::XYZVector>* mc_H_vertex = reader.Get3V("mc_H_vertex");
	std::vector<ROOT::Math::XYZTVector>* mcV1 = reader.Get4V("mcV1");
	std::vector<ROOT::Math::XYZTVector>* mcV2 = reader.Get4V("mcV2");
	std::vector<ROOT::Math::XYZTVector>* photons = reader.Get4V("photons");
	std::vector<float>* photons_r9 = reader.GetFloat("photons_r9");
	
	if (mc_H_vertex->size() != 1) continue;
	
	PV_nTracks       = reader.GetInt("PV_nTracks");
	PV_z             = reader.GetFloat("PV_z");
	PV_d0            = reader.GetFloat("PV_d0");
	PVtracks         = reader.Get3V("PVtracks");
	PVtracks_PVindex = reader.GetInt("PVtracks_PVindex");
	tracks_PVindex   = reader.GetInt("tracks_PVindex");
	tracks_dxy_PV    = reader.GetFloat("tracks_dxy_PV");
	tracks_dz_PV     = reader.GetFloat("tracks_dz_PV");
	tracks_dz        = reader.GetFloat("tracks_dz");
	sc               = reader.Get4V("photons_SC");

	hggSelection(mcV1, mcV2, photons, sc, photons_r9, accept, indpho1, indpho2);
	if (!accept) continue;

	if ( photons_r9->at(indpho1) < r9cut ) continue;
	if ( photons_r9->at(indpho2) < r9cut ) continue;
	
	etaMaxSC = sc->at(indpho1).eta();
	sum2pho  = photons->at(indpho1)+ photons->at(indpho2);
	TrueVertex_Z = mc_H_vertex->at(0).Z();
      }// end Hgg selection


      //*** selections for Zee ***
      if (isZee){
	std::vector<ROOT::Math::XYZTVector>* electrons = reader.Get4V("electrons");
	// std::vector<float>* eleid = reader.GetFloat("simpleEleId95cIso");
	// no eleID95 available for data 2011 --> compute it by hand
	std::vector<float>* eleid = new std::vector<float>;
	eleid->clear();
	if ( electrons->size() < 2) continue; 
	
	for (unsigned int iele = 0; iele < electrons->size(); iele++){
	  
	  float pt = electrons->at(iele).pt();
	  float tkIso   = reader.GetFloat("electrons_tkIsoR03")->at(iele);
	  float emIso   = reader.GetFloat("electrons_emIsoR03")->at(iele);
	  float hadIso  = reader.GetFloat("electrons_hadIsoR03_depth1")->at(iele) + reader.GetFloat("electrons_hadIsoR03_depth2")->at(iele);
	  float combIso = tkIso + emIso + hadIso;
	
	  int isEB = reader.GetInt("electrons_isEB")->at(iele);
	  float sigmaIetaIeta = reader.GetFloat("electrons_sigmaIetaIeta")->at(iele);
	  float DetaIn        = reader.GetFloat("electrons_deltaEtaIn")->at(iele);
	  float DphiIn        = reader.GetFloat("electrons_deltaPhiIn")->at(iele);
	  float HOverE        = reader.GetFloat("electrons_hOverE")->at(iele);
	  int mishits         = reader.GetInt("electrons_mishits")->at(iele);

	  float id = eleId95 ( pt, tkIso, emIso, hadIso, combIso, isEB, sigmaIetaIeta, DetaIn, DphiIn, HOverE, mishits);	
	  id *= 7.; // to emulate simpleEleId95cIso

	  eleid->push_back( id );
	}

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
	//sc               = reader.Get4V("electrons_SC");
	sc               = reader.Get4V("electrons");
	
	zeeSelection(electrons, eleid, accept, indpho1, indpho2);
	if (!accept) continue;

	etaMaxSC = sc->at(indpho1).eta();
	sum2pho  = electrons->at(indpho1)+ electrons->at(indpho2);
	TrueVertex_Z = PV_z->at(0) + (electrons_dz_PV_noEle->at(indpho1) + electrons_dz_PV_noEle->at(indpho2))/2.;

	hdiff_dZ_electrons.Fill( sum2pho.pt(), electrons_dz_PV_noEle->at(indpho1) - electrons_dz_PV_noEle->at(indpho2) );
      }

      //*** selections for Zmumu
      if ( isZmumu ){
	std::vector<ROOT::Math::XYZTVector>* muons = reader.Get4V("muons");
	std::vector<int>* muons_global = reader.GetInt("muons_global");
	std::vector<int>* muons_tracker = reader.GetInt("muons_tracker");
	std::vector<float>* muons_tkIsoR03 = reader.GetFloat("muons_tkIsoR03");
	std::vector<float>* muons_normalizedChi2 = reader.GetFloat("muons_normalizedChi2");
	std::vector<int>* muons_numberOfValidMuonHits = reader.GetInt("muons_numberOfValidMuonHits");
	std::vector<int>* muons_numberOfValidPixelHits = reader.GetInt("muons_numberOfValidPixelHits");
	std::vector<float>* muons_dxy_PV = reader.GetFloat("muons_dxy_PV");
	std::vector<float>* muons_dz_PV = reader.GetFloat("muons_dz_PV");
	std::vector<float>* muons_dz_PV_noMuon = reader.GetFloat("muons_dz_PV_noMuon");
		
	PV_nTracks       = reader.GetInt("PV_noMuon_nTracks");
	PV_z             = reader.GetFloat("PV_noMuon_z");
	PV_d0            = reader.GetFloat("PV_noMuon_d0");
	PVtracks         = reader.Get3V("PVMuonLessTracks");
	PVtracks_PVindex = reader.GetInt("PVMuonLessTracks_PVindex");
	PV_z_muon        = reader.GetFloat("PV_z");


	//************ 

	tracks_PVindex   = reader.GetInt("tracks_PVindex");
	tracks_dxy_PV    = reader.GetFloat("tracks_dxy_PV");
	tracks_dz_PV     = reader.GetFloat("tracks_dz_PV");
	tracks_dz        = reader.GetFloat("tracks_dz");
	sc               = reader.Get4V("muons");  // use muon info for SC
	
	zmumuSelection(muons,muons_global,muons_tracker, muons_tkIsoR03, 
		       muons_normalizedChi2 , 
		       muons_numberOfValidMuonHits,
		       muons_numberOfValidPixelHits,
		       muons_dxy_PV,
		       muons_dz_PV,
		       accept, indpho1, indpho2);
	
	if (!accept) continue;
	
	etaMaxSC = muons->at(indpho1).eta();
	sum2pho  = muons->at(indpho1)+ muons->at(indpho2);
	TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(indpho1) + muons_dz_PV_noMuon->at(indpho2))/2.;
		
	hdiff_dZ_muons.Fill( sum2pho.pt(), muons_dz_PV_noMuon->at(indpho1) - muons_dz_PV_noMuon->at(indpho2) );
	
	if ( fabs(muons_dz_PV_noMuon->at(indpho1) - muons_dz_PV_noMuon->at(indpho2))> 0.5) continue;
      }//Zmumu end

    
      // branches buffers
      int nvtx_;
      float  vtxx_[1000], vtxy_[1000], vtxz_[1000];
      int ntracks_;
      float tkpx_[5000], tkpy_[5000], tkpz_[5000], tkPtErr_[5000], tkWeight_[5000], 
	tkd0_[5000], tkd0Err_[5000], tkdz_[5000], tkdzErr_[5000];
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
      
      // photons
      for (int ipho = 0 ; ipho < sc->size(); ipho++){
	float px = sc->at(ipho).X();
	float py = sc->at(ipho).Y();
	float pz = sc->at(ipho).Z();
	float pt = sqrt ( px*px+py*py ); 
	float theta = 2*atan(exp(-sc->at(ipho).eta()) );
	float tantheta = tan(theta);

	if ( fabs(sc->at(ipho).eta()) < etaEB ) {
	  phocalox_[ipho] = R_ECAL*px/pt;
	  phocaloy_[ipho] = R_ECAL*py/pt;
	  phocaloz_[ipho] = R_ECAL/tantheta;
	} 

	if ( fabs(sc->at(ipho).eta()) > etaEE ) {
	  float r_endcap  = fabs(Z_ENDCAP * tantheta);
	  phocalox_[ipho] = r_endcap * px/pt;
	  phocaloy_[ipho] = r_endcap * py/pt;
	  if (pz > 0) phocaloz_[ipho] = Z_ENDCAP;
	  else phocaloz_[ipho] = -Z_ENDCAP;
	} 
	
	phoen_[ipho] = sc->at(ipho).E();
	
      }


      float eta1 = sc->at(indpho1).eta();
      float eta2 = sc->at(indpho2).eta();

      if ( (fabs(eta1) > etaEB && fabs(eta1) < etaEE) || fabs(eta1) > 2.5) continue;
      if ( (fabs(eta2) > etaEB && fabs(eta2) < etaEE) || fabs(eta2) > 2.5) continue;
   
      //*** set vertex info
      TupleVertexInfo vinfo( nvtx_, vtxx_ , vtxy_, vtxz_, ntracks_, tkpx_, tkpy_, tkpz_, tkPtErr_, tkVtxId_, tkWeight_, tkd0_, tkd0Err_,tkdz_, tkdzErr_ , tkIsHighPurity_);
           
     

     

      //*** set photon info
      PhotonInfo pho1(indpho1, TVector3(phocalox_[indpho1],phocaloy_[indpho1],phocaloz_[indpho1]),phoen_[indpho1]); 
      PhotonInfo pho2(indpho2, TVector3(phocalox_[indpho2],phocaloy_[indpho2],phocaloz_[indpho2]),phoen_[indpho2]); 
     
      
      //*** vertex analyzer
      HggVertexAnalyzer vAna(vtxAlgoParams_,nvtx_);
      vAna.setNConv(0);
      vAna.analyze(vinfo,pho1,pho2);
            
      //*** preselect vertices 
      std::vector<int> presel;
      for(int i=0; i<nvtx_; i++) {
	presel.push_back(i); 
	hz.Fill(PV_z->at(i),ww);
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
      


      //*** SUMPT2 CRITERION
//       vector<int> ranksumpt2 = vAna.rankprod(ranksumpt2_);
//       //-- matching 1cm
//       if ( fabs( TrueVertex_Z - PV_z->at(ranksumpt2[0]) ) < mindz) {
// 	PtGood.Fill( sum2pho.pt(),ww );
// 	EtaGood.Fill( etaMaxSC ,ww);
// 	NvtGood.Fill( nvtx_ ,ww);
// 	if (!isData) NpuGood.Fill(npu,ww);
//       }
//       //-- matching closest vtx
//       if ( iClosest == ranksumpt2[0]){
// 	PtGood_matchedClosest.Fill( sum2pho.pt(),ww );
// 	EtaGood_matchedClosest.Fill( etaMaxSC ,ww);
// 	NvtGood_matchedClosest.Fill( nvtx_ ,ww);
// 	if (!isData) NpuGood_matchedClosest.Fill(npu,ww);
//       }        

   
      //*** RANKING PRODUCT
//       vector<int> rankprod = vAna.rankprod(rankVariables_);
//       //-- matching 1cm
//       if ( fabs( TrueVertex_Z - PV_z->at(rankprod[0]) ) < mindz) {
// 	PtGood_RANK.Fill( sum2pho.pt(),ww );
// 	EtaGood_RANK.Fill( etaMaxSC ,ww);
// 	NvtGood_RANK.Fill( nvtx_,ww );
// 	if (!isData) NpuGood_RANK.Fill(npu,ww);
//       }
//       //-- matching closest vtx
//       if ( iClosest == rankprod[0]){
// 	PtGood_RANK_matchedClosest.Fill( sum2pho.pt(),ww );
// 	EtaGood_RANK_matchedClosest.Fill( etaMaxSC ,ww);
// 	NvtGood_RANK_matchedClosest.Fill( nvtx_,ww );
// 	if (!isData) NpuGood_RANK_matchedClosest.Fill(npu,ww);
//       }
      
      int vtxbin=0;
      if (nvtx_ >= 10 && nvtx_< 21 ) vtxbin=1;
      if (nvtx_ >= 21 ) vtxbin=2;


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
      ChosenVertexDz_BDT_vs_pt.Fill(sum2pho.pt(),TrueVertex_Z - PV_z->at(ranktmva[0]),ww);

      float vtxmva, evtmva, vtxprob;
      // fill per vertex mva 
      for (int iv=0;iv<ranktmva.size();iv++) {
	float vtxmva = vAna.mva(ranktmva[iv]);
	if ( iClosest == ranktmva[iv] ) {
	  BDToutput_sig.Fill( vtxmva, ww );
	  sumpt2_sig.Fill(vAna.logsumpt2(ranktmva[iv]),ww ) ;
	  ptasym_sig.Fill(vAna.ptasym(ranktmva[iv]),ww );
	  ptbal_sig .Fill(vAna.ptbal(ranktmva[iv]),ww );
	  perVertexBDToutput_sig_vtxbin[vtxbin]->Fill( vtxmva, ww );
	}
	else {
	  BDToutput_bkg.Fill( vtxmva, ww );
	  sumpt2_bkg.Fill(vAna.logsumpt2(ranktmva[iv]),ww ) ;
	  ptasym_bkg.Fill(vAna.ptasym(ranktmva[iv]),ww );
	  ptbal_bkg .Fill(vAna.ptbal(ranktmva[iv]),ww );
	  perVertexBDToutput_bkg_vtxbin[vtxbin]->Fill( vtxmva, ww );
	}
      }
      // fill per event mva
      evtmva = vAna.perEventMva(*tmvaPerEvtReader_, tmvaEventMethod.c_str(),ranktmva);
      //      vtxprob = vAna.vertexProbability(evtmva);
      vtxprob = fprob->Eval(evtmva);
      float wid = 1;
      if (applyVertexIdScaleFactor) wid = myVtxIdEffRw->GetWeight(sum2pho.pt());
      perEventBDToutput.Fill( evtmva, ww*wid );
      vtxProbability.Fill( vtxprob, ww*wid );
      perEventBDToutput_vtxbin[vtxbin]->Fill( evtmva, ww*wid );
      if (fabs( TrueVertex_Z - PV_z->at(ranktmva[0]) ) < mindz) {
      	perEventBDToutput_sig.Fill( evtmva, ww*wid );
  	perEventBDToutput_sig_vtxbin[vtxbin]->Fill( evtmva, ww*wid );
	vtxProbability_sig.Fill( vtxprob, ww*wid );
      }
      else{
	perEventBDToutput_bkg.Fill( evtmva, ww*wid );
	perEventBDToutput_bkg_vtxbin[vtxbin]->Fill( evtmva, ww*wid );
	vtxProbability_bkg.Fill( vtxprob, ww*wid );
      }


      //-- matching 1cm
      if ( fabs( TrueVertex_Z - PV_z->at(ranktmva[0]) ) < mindz) {
	PtGood_BDT.Fill( sum2pho.pt(),ww );
	EtaGood_BDT.Fill( etaMaxSC ,ww);
	NvtGood_BDT.Fill( nvtx_ ,ww);
	if (!isData) NpuGood_BDT.Fill(npu,ww);
      }
      //-- matching closest vtx
      if ( iClosest == ranktmva[0]){
	PtGood_BDT_matchedClosest.Fill( sum2pho.pt(),ww );
	EtaGood_BDT_matchedClosest.Fill( etaMaxSC ,ww);
	NvtGood_BDT_matchedClosest.Fill( nvtx_,ww );
	if (!isData) NpuGood_BDT_matchedClosest.Fill(npu,ww);
      }
      

      // -- all vtx
      PtAll.Fill( sum2pho.pt(),ww );
      if (fabs(eta1) < etaEB && fabs(eta2) < etaEB)  PtAll_EBEB.Fill( sum2pho.pt(),ww );
      if (fabs(eta1) > etaEE && fabs(eta2) > etaEE)  PtAll_EEEE.Fill( sum2pho.pt(),ww );
      if ( (fabs(eta1) < etaEB && fabs(eta2) > etaEE) || (fabs(eta2) < etaEB && fabs(eta1) > etaEE))  PtAll_EBEE.Fill( sum2pho.pt(),ww );

      EtaAll.Fill( etaMaxSC ,ww);
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
  
  EtaAll.Write();
  EtaGood.Write();
  EtaGood_BDT.Write();
  EtaGood_RANK.Write();
  
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
  
  EtaGood_matchedClosest.Write();
  EtaGood_BDT_matchedClosest.Write();
  EtaGood_RANK_matchedClosest.Write();
  
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
  vtxProbability.Write();
  vtxProbability_sig.Write();
  vtxProbability_bkg.Write();
  for (int vtxbin=0; vtxbin<3; vtxbin++){
    perEventBDToutput_vtxbin[vtxbin]-> Write();
    perEventBDToutput_sig_vtxbin[vtxbin]-> Write();
    perEventBDToutput_bkg_vtxbin[vtxbin]-> Write();
    perVertexBDToutput_vtxbin[vtxbin]-> Write();
    perVertexBDToutput_sig_vtxbin[vtxbin]-> Write();
    perVertexBDToutput_bkg_vtxbin[vtxbin]-> Write();
  }

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











