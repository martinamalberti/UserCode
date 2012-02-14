
{
  gROOT->SetStyle("Plain");


  float etaMax = 1.444;

  bool usePUweights = true;

  //--- weights for MC
  TFile weightsFile("weights/PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","READ"); // stessi pesi usati per analisi vertici Hgg  
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();
    


  TChain *ntu_MC = new TChain("ntu");
  TChain *ntu_Data = new TChain("ntu");
    
  //---- MC fall 2011
  ntu_MC->Add("/data2/calibrator/NTUPLES/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11_All.root");
   
  //---- DATA
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v6_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_data_20111122_158851_180363.root");



 std::cout << "     MC    : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;

  //---- Observables
  int npu;
  int isEB;
  float EoP, scEta, scPhi;
  float scE3x3, scE5x5, scE, scEt;  
  float charge, scLocalEta, scLocalPhi,crackCorr,scLocalCorr; 
  float R9;

  //---- Set branch addresses for MC  
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("ele1_isEB", &isEB);
  ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_MC->SetBranchAddress("ele1_scE", &scE);
  ntu_MC->SetBranchAddress("ele1_scEt", &scEt);
  ntu_MC->SetBranchAddress("ele1_charge", &charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_MC->SetBranchAddress("ele1_scLocalContCorr",&scLocalCorr); 


  //---- Set branch addresses for Data
  ntu_Data->SetBranchAddress("ele1_isEB", &isEB);
  ntu_Data->SetBranchAddress("ele1_scEta", &scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_Data->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_Data->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_Data->SetBranchAddress("ele1_scE", &scE);
  ntu_Data->SetBranchAddress("ele1_scEt", &scEt);
  ntu_Data->SetBranchAddress("ele1_charge", &charge);
  ntu_Data->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_Data->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_Data->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_Data->SetBranchAddress("ele1_scLocalContCorr",&scLocalCorr); 



  TH1F *hr9EBData = new TH1F("hr9EBData","hr9EBdata",600,0,1.5);
  TH1F *hr9EBMC   = new TH1F("hr9EBMC","hr9EBmc",600,0,1.5);

  TH1F *hr9EEData = new TH1F("hr9EEData","hr9EEdata",600,0,1.5);
  TH1F *hr9EEMC   = new TH1F("hr9EEMC","hr9EEmc",600,0,1.5);

  TH1F *hscaledr9EBMC   = new TH1F("hscaledr9EBMC","hscaledr9EBmc",600,0,1.5);
  TH1F *hscaledr9EEMC   = new TH1F("hscaledr9EEMC","hscaledr9EEmc",600,0,1.5);


  TH1F *hscEtaEBData = new TH1F("hscEtaEBData","hscEtaEBData",600,-2.5,2.5);
  TH1F *hscEtaEBMC   = new TH1F("hscEtaEBMC","hscEtaEBMC",600,-2.5,2.5);

  TH1F *hscEtaEEData = new TH1F("hscEtaEEData","hscEtaEEData",600,-2.5,2.5);
  TH1F *hscEtaEEMC   = new TH1F("hscEtaEEMC","hscEtaEEMC",600,-2.5,2.5);

  TH1F *hscEtaData = new TH1F("hscEtaData","hscEtaData",600,-2.5,2.5);
  TH1F *hscEtaMC   = new TH1F("hscEtaMC","hscEtaMC",600,-2.5,2.5);

  TH1F *hscEtData = new TH1F("hscEtData","hscEtData",200,0,200);
  TH1F *hscEtMC   = new TH1F("hscEtMC","hscEtMC",200,0,200);

  TH1F *hscEtEBData = new TH1F("hscEtEBData","hscEtEBData",200,0,200);
  TH1F *hscEtEBMC   = new TH1F("hscEtEBMC","hscEtEBMC",200,0,200);

  TH1F *hscEtEEData = new TH1F("hscEtEEData","hscEtEEData",200,0,200);
  TH1F *hscEtEEMC   = new TH1F("hscEtEEMC","hscEtEEMC",200,0,200);


  //******************************************************************************************
  //*************************************** MC  ********************************************** 
  std::cout << "Loop on MC events ... " << endl; 
  float ww = 1 ;

  //---- loop on MC, make refernce and fit dist
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%200000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //if (entry>500000) break;

    ntu_MC->GetEntry(entry);

    if (usePUweights) ww = w[npu];

    R9 = scE3x3/scE;    

    hscEtaMC->Fill(scEta,ww);
    hscEtMC ->Fill(scEt,ww);

    if (isEB) {
      hscEtaEBMC->Fill(scEta,ww);
      hscEtEBMC ->Fill(scEt,ww);
      hr9EBMC   ->Fill(R9,ww);
      hscaledr9EBMC   ->Fill(R9*1.005,ww);
    }
    else{
      hscEtaEEMC->Fill(scEta,ww);
      hscEtEEMC ->Fill(scEt,ww);
      hr9EEMC   ->Fill(R9,ww);
      hscaledr9EEMC   ->Fill(R9*1.005,ww);
    }

  }


 //******************************************************************************************
  //*************************************** DATA ********************************************** 
  std::cout << "Loop on Data events ..." << endl; 
  //--- loop on data
  for(int entry = 0; entry < ntu_Data->GetEntries(); ++entry) {
    if( entry%200000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //if (entry>500000) break;

    ntu_Data->GetEntry(entry);


    R9 = scE3x3/scE;    

    hscEtaData->Fill(scEta);
    hscEtData ->Fill(scEt);

    if (isEB) {
      hscEtaEBData->Fill(scEta);
      hscEtEBData ->Fill(scEt);
      hr9EBData   ->Fill(R9);
    }
    else{
      hscEtaEEData->Fill(scEta);
      hscEtEEData ->Fill(scEt);
      hr9EEData   ->Fill(R9);
    }

  }


  hr9EBMC   -> Sumw2();
  hr9EBData -> Sumw2();

  hr9EEMC   -> Sumw2();
  hr9EEData -> Sumw2();

  hscaledr9EBMC-> Sumw2();
  hscaledr9EEMC-> Sumw2();

  hscEtaMC    -> Sumw2();
  hscEtaData  -> Sumw2();

  hscEtaEBMC    -> Sumw2();
  hscEtaEBData  -> Sumw2();

  hscEtaEEMC    -> Sumw2();
  hscEtaEEData  -> Sumw2();

  hscEtMC    -> Sumw2();
  hscEtData  -> Sumw2();

  hscEtEBMC    -> Sumw2();
  hscEtEBData  -> Sumw2();

  hscEtEEMC    -> Sumw2();
  hscEtEEData  -> Sumw2();




  TCanvas *c = new TCanvas("c","c",500,500);
  hr9EBMC  -> Draw("histo");
  hr9EBData->SetMarkerStyle(20);
  hr9EBData->SetMarkerSize(1);
  hr9EBData-> Draw("esame");


  float nEBMC = hr9EBMC->GetEntries();
  float nEBData = hr9EBData->GetEntries();
  float nEEMC = hr9EEMC->GetEntries();
  float nEEData = hr9EEData->GetEntries();

  TH1F * hr9ratioEB = (TH1F*)hr9EBData->Clone("hr9ratioEB");
  hr9ratioEB->Divide(hr9ratioEB,hr9EBMC,1./nEBData, 1./nEBMC);

  TH1F * hr9ratioEE = (TH1F*)hr9EEData->Clone("hr9ratioEE");
  hr9ratioEE->Divide(hr9ratioEE,hr9EEMC,1./nEEData, 1./nEEMC);

  TFile *fout = new TFile("r9weights.root","recreate");

  hscEtaMC    -> Write();
  hscEtaData  -> Write();

  hscEtaEBMC    -> Write();
  hscEtaEBData  -> Write();

  hscEtaEEMC    -> Write();
  hscEtaEEData  -> Write();

  hscEtMC    -> Write();
  hscEtData  -> Write();

  hscEtEBMC    -> Write();
  hscEtEBData  -> Write();

  hscEtEEMC    -> Write();
  hscEtEEData  -> Write();

  hr9EBMC    -> Write();
  hr9EBData  -> Write();
  hr9ratioEB -> Write();

  hr9EEMC    -> Write();
  hr9EEData  -> Write();
  hr9ratioEE -> Write();



}
