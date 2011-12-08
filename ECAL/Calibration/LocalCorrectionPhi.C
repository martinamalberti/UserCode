#include "histoFunc.h"

//derive f(local eta) correction in MC, with 4 different f for each module.
// foldi crystal in eta into one single module

bool IsEtaGap(float eta){
  float feta = fabs(eta);
  if( fabs(feta - 0 )<3) return true;
  if( fabs(feta - 25)<3) return true;
  if( fabs(feta - 45)<3) return true;
  if( fabs(feta - 65)<3) return true;
  if( fabs(feta - 85)<3) return true;
  return false;
}

int templIndex(float eta){
    float feta = fabs(eta);
    if (feta <= 25)            {return 0;}
    if (feta> 25 && feta <= 45){return 1;}
    if (feta> 45 && feta <= 65){return 2;}
    if (feta> 65 && feta <= 85){return 3;}

    return -1;
}

void LocalCorrectionPhi()
{

  gROOT->SetStyle("Plain");
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.05);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.05);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.05);
  gROOT->ForceStyle();

  float xtalWidth = 0.01745329;

  //---- variables for selections
  float etaMax = 1.44;
  
  float r9min = 0. ;
  float r9max = 999999. ;

  bool useOddCry  = false;
  bool useEvenCry = false;
  bool useEvenCry = false;
  bool usePUweights = true;

  //--- weights for MC
  TFile weightsFile("weights/PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","READ"); // stessi pesi usati per analisi vertici Hgg  
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();

  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"GraphsLocalPhi.root");
   
  // NTUPLES 
  TChain *ntu_MC = new TChain("ntu");
  TChain *ntu_Data = new TChain("ntu");
  
  //--- MC Fall 2011
  ntu_MC->Add("/data2/calibrator/NTUPLES/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11_All.root");
   
  //--- Data
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v6_data_20111122_158851_180363.root");
  ntu_Data->Add("/data2/calibrator/NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_data_20111122_158851_180363.root");
 


  std::cout << "     MC    : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;

  //---- observables
  int npu;
  float EoP, scEta, scPhi;
  float scE3x3, scE5x5, scEne;  
  float charge, scLocalEta, scLocalPhi,crackCorr,scLocalCorr; 

  float R9;
  
  //---- Set branch addresses for MC  
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_MC->SetBranchAddress("ele1_scE", &scEne);
  ntu_MC->SetBranchAddress("ele1_charge", &charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_MC->SetBranchAddress("ele1_scLocalContCorr",&scLocalCorr); 


  //---- Set branch addresses for Data
  ntu_Data->SetBranchAddress("ele1_scEta", &scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_Data->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_Data->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_Data->SetBranchAddress("ele1_scE", &scEne);
  ntu_Data->SetBranchAddress("ele1_charge", &charge);
  ntu_Data->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_Data->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_Data->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_Data->SetBranchAddress("ele1_scLocalContCorr",&scLocalCorr); 
  
  const unsigned int nBins = 20;
  const int Ntempl = 4;
  std::cout << "nBins = " << nBins << std::endl;
  
  // histogram definition
  TH1F* h_EoP_MC[nBins][Ntempl] ;   
  TH1F* h_EoC_MC[nBins][Ntempl] ;
  TH1F* h_EoP_Data[nBins][Ntempl];
  TH1F* h_EoC_Data[nBins][Ntempl] ;

  for(int mod=0; mod<Ntempl; mod++){
    for(unsigned int i = 0; i < nBins; ++i)
      {
	char histoName[80];

	sprintf(histoName, "EoP_MC_%d_mod%d", i,mod+1);
	h_EoP_MC[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoP_MC[i][mod] -> SetFillColor(4);
	h_EoP_MC[i][mod] -> SetFillStyle(3004);
	sprintf(histoName, "EoC_MC_%d_mod%d", i,mod+1);
	h_EoC_MC[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoC_MC[i][mod] -> SetFillColor(3);
	h_EoC_MC[i][mod] -> SetFillStyle(3004);

	sprintf(histoName, "EoP_Data_%d_mod%d", i,mod+1);
	h_EoP_Data[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoP_Data[i][mod] -> SetFillColor(4);
	h_EoP_Data[i][mod] -> SetFillStyle(3004);
	sprintf(histoName, "EoC_Data_%d_mod%d", i,mod+1);
	h_EoC_Data[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoC_Data[i][mod] -> SetFillColor(3);
	h_EoC_Data[i][mod] -> SetFillStyle(3004);
      }
  }

  
  //---- book templates
  TH1F* h_template_MC[Ntempl];
  TH1F* h_template_Data[Ntempl];
  for(unsigned int i = 0; i < Ntempl; ++i){
    char histoName[100];
    sprintf(histoName, "template_MC_%d", i);
    h_template_MC[i] = new TH1F(histoName, "", 1200, 0., 3.);   
    sprintf(histoName, "template_Data_%d", i);
    h_template_Data[i] = new TH1F(histoName, "", 1200, 0., 3.);   
  }

  
  TH2F* h_Corr_MC_lowR9 = new TH2F("h_Corr_MC_lowR9","Correction in MC",100,0,1,200,0.9,1.1);
  TH2F* h_Corr_MC_highR9  = new TH2F("h_Corr_MC_highR9","Correction in MC",100,0,1,200,0.9,1.1);
  TProfile* pr_Corr_MC_lowR9 =  new TProfile("pr_Corr_MC_lowR9","Correction profile MC",100,0,1,0.95,1.05);
  pr_Corr_MC_lowR9 -> SetMarkerColor(2); 
  pr_Corr_MC_lowR9 -> SetLineColor(2);
  TProfile* pr_Corr_MC_highR9 =  new TProfile("pr_Corr_MC_highR9","Correction profile MC",100,0,1,0.95,1.05);
  pr_Corr_MC_highR9 -> SetMarkerColor(2); 
  pr_Corr_MC_highR9 -> SetLineColor(2);

  TH2F* h_Corr_Data_lowR9  = new TH2F("h_Corr_Data_lowR9","Correction in Data",100,0,1,200,0.9,1.1);
  TH2F* h_Corr_Data_highR9 = new TH2F("h_Corr_Data_highR9","Correction in Data",100,0,1,200,0.9,1.1);
  TProfile* pr_Corr_Data_lowR9 =  new TProfile("pr_Corr_Data_lowR9","Correction profile Data",100,0,1,0.95,1.05);
  pr_Corr_Data_lowR9 -> SetMarkerColor(4); 
  pr_Corr_Data_lowR9 -> SetLineColor(4);
  TProfile* pr_Corr_Data_highR9 =  new TProfile("pr_Corr_Data_highR9","Correction profile Data",100,0,1,0.95,1.05);
  pr_Corr_Data_highR9 -> SetMarkerColor(4); 
  pr_Corr_Data_highR9 -> SetLineColor(4);
 


  //******************************************************************************************
  //*************************************** MC  ********************************************** 
  std::cout << "Loop over MC events ... " << endl; 
  float ww = 1 ;
 
  //---- loop on MC, make reference and fit dist
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%200000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //    if (entry>1000) break;

    ntu_MC->GetEntry(entry);

    if (usePUweights) ww = w[npu];

    R9 = scE3x3/scEne;

    //-- eta or R9 cuts
    if ( fabs(scEta) > etaMax ) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 
    //if ( (scEta*charge)>0 ) continue;
      
    //-- remove phi cracks
    float phi = (scPhi+3.1415926536)/xtalWidth;
    float modphi = (int)phi%20;
    if (fabs(modphi-10)<2.) continue;
    
    //-- use only even/odd crystals (for data)
    if ( useOddCry  && int(modphi)%2 == 0) continue; 
    if ( useEvenCry && int(modphi)%2 != 0) continue; 
 
    //-- remove gaps
    float fetaCry = fabs (scEta) / xtalWidth;
    if( IsEtaGap(fetaCry) ) continue;
    
    int mod = templIndex(fetaCry);
            
    //-- fill template for each mod
    h_template_MC[mod]-> Fill(EoP,ww);

    //-- fill MC histos in phi bins
    float locPhi = scLocalPhi+0.5; 
    if ( fabs(locPhi-0.5) >= 0.5 ) continue;
 
    int bin = nBins * locPhi; 
    if (bin>nBins-1 || bin < 0 ) {
      cout << "Error in bins: " << bin << " " << scLocalPhi << endl;
      continue;
    }
    
    float correction = scLocalCorr;
    correction = 1;
    crackCorr  = 1;
    
    (h_EoP_MC[bin][mod]) -> Fill(EoP/crackCorr,ww);
    (h_EoC_MC[bin][mod]) -> Fill(EoP/crackCorr*correction,ww);

    if (R9 > 0.94) {
      pr_Corr_MC_highR9->Fill(locPhi,correction);
      h_Corr_MC_highR9->Fill(locPhi,correction);
    }
    else {
      pr_Corr_MC_lowR9->Fill(locPhi,correction);
      h_Corr_MC_lowR9->Fill(locPhi,correction);

    }
  }
  cout << "stat: "<< h_EoP_MC[0][nBins/2]->GetEntries()<< " " <<  h_EoC_MC[0][nBins/2]->GetEntries() << endl; 

  
    
  //******************************************************************************************
  //*************************************** DATA ********************************************** 
  std::cout << "Loop on Data events ... " << endl; 
  //---- loop on data
  for(int entry = 0; entry < ntu_Data->GetEntries(); ++entry) {
    if( entry%200000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //    if (entry>1000) break;

    ntu_Data->GetEntry(entry);

    R9 = scE3x3/scEne;

    //-- eta or R9 cuts
    if ( fabs(scEta) > etaMax ) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 
    //if ( (scEta*charge)>0 ) continue;

    //-- remove phi cracks
    float phi = (scPhi+3.141592653)/xtalWidth;
    float modphi = (int)phi%20;
    if (fabs(modphi-10)<2.) continue;
        
    //-- use only even/odd crystals (for data)
    if ( useOddCry  && int(modphi)%2 == 0) continue; 
    if ( useEvenCry && int(modphi)%2 != 0) continue; 

    //-- remove eta gaps
    float fetaCry = fabs (scEta) / xtalWidth;
    if( IsEtaGap(fetaCry) ) continue;
    
    int mod = templIndex(fetaCry);

    //-- fill template for each mod
    h_template_Data[mod]-> Fill(EoP);

    //-- fill data histos in phi bins
    float locPhi = scLocalPhi + 0.5; 
    if ( fabs(locPhi-0.5) >= 0.5 ) continue;
     
    int bin = nBins * (locPhi); 
    if (bin>nBins-1 || bin < 0 ) {
      cout << "Error in bins: " << bin << " " << scLocalPhi << endl;
      continue;
    }

    float correction = scLocalCorr;
    correction = 1;
    crackCorr  = 1;

    (h_EoP_Data[bin][mod]) -> Fill(EoP/crackCorr);
    (h_EoC_Data[bin][mod]) -> Fill(EoP/crackCorr*correction);

    if (R9 > 0.94) {
      pr_Corr_Data_highR9->Fill(locPhi,correction);
      h_Corr_Data_highR9->Fill(locPhi,correction);
    }
    else {
      pr_Corr_Data_lowR9->Fill(locPhi,correction);
      h_Corr_Data_lowR9->Fill(locPhi,correction);
    }

  }
  
  cout << "stat: "<< h_EoP_Data[nBins/2][0]->GetEntries()<< " " <<  h_EoC_Data[nBins/2][0]->GetEntries() << endl; 
  

  ///
  TCanvas* ccorr_highR9 = new TCanvas("ccorrection_highR9","correction_highR9",100,100,500,700);
  ccorr_highR9 ->Divide(1,2);
  ccorr_highR9->cd(1);
  h_Corr_MC_highR9->Draw();  pr_Corr_MC_highR9->Draw("same");
  ccorr_highR9->cd(2);
  h_Corr_Data_highR9->Draw();  pr_Corr_Data_highR9->Draw("same");
  
  TCanvas* ccorr_lowR9 = new TCanvas("ccorrection_lowR9","correction_lowR9",100,100,500,700);
  ccorr_lowR9 ->Divide(1,2);
  ccorr_lowR9->cd(1);
  h_Corr_MC_lowR9->Draw();  pr_Corr_MC_lowR9->Draw("same");
  ccorr_lowR9->cd(2);
  h_Corr_Data_lowR9->Draw();  pr_Corr_Data_lowR9->Draw("same");
  

 
  ///////////////****************** Fit the histograms and fill the graphs *************** ////////////////////////

  int rebin = 4;
  
  TGraphErrors* g_EoP_MC[Ntempl];
  TGraphErrors* g_EoC_MC[Ntempl];
  TGraphErrors* g_ratio_MC[Ntempl];

  TGraphErrors* g_EoP_Data[Ntempl];
  TGraphErrors* g_EoC_Data[Ntempl];
  TGraphErrors* g_ratio_Data[Ntempl];

  TGraphErrors* g_ratio_uncorr[Ntempl];
  TGraphErrors* g_ratio_corr[Ntempl];

  TH1F* spread_EoP_MC[Ntempl];
  TH1F* spread_EoC_MC[Ntempl];
  
  TH1F* spread_EoP_Data[Ntempl];
  TH1F* spread_EoC_Data[Ntempl];
  
  for (int mod=0; mod<4; mod++){
    char histoName[100];
    
    sprintf(histoName, "gEoP_MC_mod%d", mod+1);
    g_EoP_MC[mod]   = new TGraphErrors(); 
    g_EoP_MC[mod]->SetName(histoName);
    
    sprintf(histoName, "gEoC_MC_mod%d", mod+1);
    g_EoC_MC[mod]   = new TGraphErrors(); 
    g_EoC_MC[mod]->SetName(histoName);
    
    sprintf(histoName, "gRatio_MC_mod%d", mod+1);
    g_ratio_MC[mod]   = new TGraphErrors(); 
    g_ratio_MC[mod]->SetName(histoName);

    sprintf(histoName, "gEoP_Data_mod%d", mod+1);
    g_EoP_Data[mod]   = new TGraphErrors(); 
    g_EoP_Data[mod]->SetName(histoName);
    
    sprintf(histoName, "gEoC_Data_mod%d", mod+1);
    g_EoC_Data[mod]   = new TGraphErrors(); 
    g_EoC_Data[mod]->SetName(histoName);
    
    sprintf(histoName, "gRatio_Data_mod%d", mod+1);
    g_ratio_Data[mod]   = new TGraphErrors(); 
    g_ratio_Data[mod]->SetName(histoName);

    sprintf(histoName, "gRatio_uncorr_mod%d", mod+1);
    g_ratio_uncorr[mod]   = new TGraphErrors(); 
    g_ratio_uncorr[mod]->SetName(histoName);

    sprintf(histoName, "gRatio_corr_mod%d", mod+1);
    g_ratio_corr[mod]   = new TGraphErrors(); 
    g_ratio_corr[mod]->SetName(histoName);

    if (mod<3) {
      h_template_MC[mod]   -> Rebin(rebin);
      h_template_Data[mod] -> Rebin(rebin);
    }
    else {
      h_template_MC[mod] -> Rebin(rebin*2);
      h_template_Data[mod] -> Rebin(rebin*2);    
    }

    sprintf(histoName, "spreadEoP_MC_mod%d", mod+1);
    spread_EoP_MC[mod]   = new TH1F(histoName, histoName,200,0.9,1.1); 

    sprintf(histoName, "spreadEoC_MC_mod%d", mod+1);
    spread_EoC_MC[mod]   = new TH1F(histoName, histoName,200,0.9,1.1); 

    sprintf(histoName, "spreadEoP_Data_mod%d", mod+1);
    spread_EoP_Data[mod]   = new TH1F(histoName, histoName,200,0.9,1.1); 

    sprintf(histoName, "spreadEoC_Data_mod%d", mod+1);
    spread_EoC_Data[mod]   = new TH1F(histoName, histoName,200,0.9,1.1); 
    
  }


  /////////////  MC ////////////////////
    
  for(int mod=0;mod<4;mod++){

    histoFunc *templateHistoFuncMC   = new histoFunc(h_template_MC[mod]);
    histoFunc *templateHistoFuncData = new histoFunc(h_template_Data[mod]);

    for(unsigned int i = 0; i < nBins; ++i)
      {
	std::cout << "***** Fitting :  mod: " <<mod <<"  bin:  "<< i << endl; 

	if (mod < 3){
	  h_EoP_MC[i][mod] -> Rebin(rebin);    
	  h_EoC_MC[i][mod] -> Rebin(rebin);    
	  h_EoP_Data[i][mod] -> Rebin(rebin);    
	  h_EoC_Data[i][mod] -> Rebin(rebin);    
	}
	else {
	  h_EoP_MC[i][mod] -> Rebin(rebin*2);    
	  h_EoC_MC[i][mod] -> Rebin(rebin*2);    
	  h_EoP_Data[i][mod] -> Rebin(rebin*2);    
	  h_EoC_Data[i][mod] -> Rebin(rebin*2);
	}

	//************************ MC ****************************************************************
	TF1 * templateFuncMC = new TF1("templateFuncMC", templateHistoFuncMC, 0.7, 1.3, 3, "histoFunc");
	templateFuncMC -> SetNpx(10000);

	float xval = (i+0.5)*1/(float)nBins - 0.5;

	//--- uncorrected MC   
	double xNorm = h_EoP_MC[i][mod]->GetEntries()/h_template_MC[mod]->GetEntries() *
	               h_EoP_MC[i][mod]->GetBinWidth(1)/h_template_MC[mod]->GetBinWidth(1); 

	templateFuncMC -> FixParameter(0, xNorm);
	templateFuncMC -> SetParameter(1, 1.05 );
	templateFuncMC -> FixParameter(2, 0.);
    

	h_EoP_MC[i][mod] -> Fit("templateFuncMC", "MRQLN+");
	g_EoP_MC[mod] -> SetPoint(i, xval , 1./templateFuncMC->GetParameter(1));
	g_EoP_MC[mod] -> SetPointError(i, 0., templateFuncMC->GetParError(1));
	spread_EoP_MC[mod] -> Fill(1./templateFuncMC->GetParameter(1));
	//    if ( templateFuncMC->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
	//    cout  << " ***** " <<  1./templateFuncMC->GetParameter(1) << " " << templateFuncMC->GetParError(1) << endl; 

	float scaleMC  = 1./templateFuncMC->GetParameter(1);
	float escaleMC = templateFuncMC->GetParError(1)/scaleMC/scaleMC;

	//--- corrected MC  
	xNorm = h_EoC_MC[i][mod]->GetEntries()/h_template_MC[mod]->GetEntries() *
	        h_EoC_MC[i][mod]->GetBinWidth(1)/h_template_MC[mod]->GetBinWidth(1); 

	templateFuncMC -> FixParameter(0, xNorm);
	templateFuncMC -> SetParameter(1, 1.05 );
	templateFuncMC -> FixParameter(2, 0.);

	h_EoC_MC[i][mod] -> Fit("templateFuncMC", "MRQLN+");
	g_EoC_MC[mod] -> SetPoint(i, xval , 1./templateFuncMC->GetParameter(1));
	g_EoC_MC[mod] -> SetPointError(i, 0., templateFuncMC->GetParError(1));
	spread_EoC_MC[mod] -> Fill(1./templateFuncMC->GetParameter(1));

    	float scaleMCcorr  = 1./templateFuncMC->GetParameter(1);
	float escaleMCcorr = templateFuncMC->GetParError(1)/scaleMCcorr/scaleMCcorr;




	//************************ DATA **********************************************
	TF1 * templateFuncData = new TF1("templateFuncData", templateHistoFuncData, 0.7, 1.3, 3, "histoFunc");
	templateFuncData -> SetNpx(10000);

	//--- uncorrected Data   
	xNorm = h_EoP_Data[i][mod]->GetEntries()/h_template_Data[mod]->GetEntries() *
	        h_EoP_Data[i][mod]->GetBinWidth(1)/h_template_Data[mod]->GetBinWidth(1); 
 
	templateFuncData -> FixParameter(0, xNorm);
	templateFuncData -> SetParameter(1, 1.05 );
	templateFuncData -> FixParameter(2, 0.);
    
	h_EoP_Data[i][mod] -> Fit("templateFuncData", "MRQLN+");
	g_EoP_Data[mod] -> SetPoint(i, xval , 1./templateFuncData->GetParameter(1));
	g_EoP_Data[mod] -> SetPointError(i, 0., templateFuncData->GetParError(1));
	spread_EoP_Data[mod] -> Fill(1./templateFuncData->GetParameter(1));
	//    if ( templateFuncData->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
	//    cout  << " ***** " <<  1./templateFuncData->GetParameter(1) << " " << templateFuncData->GetParError(1) << endl; 
	
	float scaleDA  = 1./templateFuncData->GetParameter(1);
	float escaleDA = templateFuncData->GetParError(1)/scaleDA/scaleDA;

	//--- corrected Data  
	xNorm = h_EoC_Data[i][mod]->GetEntries()/h_template_Data[mod]->GetEntries() *
	  h_EoC_Data[i][mod]->GetBinWidth(1)/h_template_Data[mod]->GetBinWidth(1); 

	templateFuncData -> FixParameter(0, xNorm);
	templateFuncData -> SetParameter(1, 1.05 );
	templateFuncData -> FixParameter(2, 0.);
    
	h_EoC_Data[i][mod] -> Fit("templateFuncData", "SNQR+");
	g_EoC_Data[mod] -> SetPoint(i, xval, 1./templateFuncData->GetParameter(1));
	g_EoC_Data[mod] -> SetPointError(i, 0., templateFuncData->GetParError(1));
	spread_EoC_Data[mod] -> Fill(1./templateFuncData->GetParameter(1));
    	
	float scaleDAcorr  = 1./templateFuncData->GetParameter(1);
	float escaleDAcorr = templateFuncData->GetParError(1)/scaleDAcorr/scaleDAcorr;
	
	
	
	//--- ratio finalization MC corr/uncorr
	float ratioMC = scaleMC/scaleMCcorr;
	float eratioMC = ratioMC*sqrt(pow(escaleMC/scaleMC,2) + pow(escaleMCcorr/scaleMCcorr,2)); 
	g_ratio_MC[mod] -> SetPoint(i, xval , ratioMC);
	g_ratio_MC[mod] -> SetPointError(i, 0., eratioMC);

	//--- ratio finalization DATA corr/uncorr
	float ratioDA  = scaleDA/scaleDAcorr;
	float eratioDA = ratioDA*sqrt(pow(escaleDA/scaleDA,2) + pow(escaleDAcorr/scaleDAcorr,2)); 
	g_ratio_Data[mod] -> SetPoint(i, xval, ratioDA);
	g_ratio_Data[mod] -> SetPointError(i, 0., eratioDA);


	//--- ratio finalization data/MC uncorrected
	float ratioU = scaleDA/scaleMC;
	float eratioU = ratioU*sqrt(pow(escaleMC/scaleMC,2) + pow(escaleDA/scaleDA,2)); 
	g_ratio_uncorr[mod] -> SetPoint(i, xval , ratioU);
	g_ratio_uncorr[mod] -> SetPointError(i, 0., eratioU);

	//--- ratio finalization data/MC corrected
	float ratioC = scaleDAcorr/scaleMCcorr;
	float eratioC = ratioC*sqrt(pow(escaleMCcorr/scaleMCcorr,2) + pow(escaleDAcorr/scaleDAcorr,2)); 
	g_ratio_corr[mod] -> SetPoint(i, xval , ratioC);
	g_ratio_corr[mod] -> SetPointError(i, 0., eratioC);
    
      }
  }
  
  TCanvas* c_g_fit[Ntempl];
  TLegend* tl[Ntempl];
  TLegend* tlr[Ntempl];

  for(int mod=0;mod<4;mod++) {
    char padName[100];
    sprintf(padName, "g_fit_mod%d", mod+1);
    c_g_fit[mod] = new TCanvas(padName,padName,100,100,700,600);
    c_g_fit[mod]->Divide(1,2);
    c_g_fit[mod]->cd(1);
    gPad->SetGrid();
    TH1F *hPad = (TH1F*)gPad->DrawFrame(-0.55,0.97,0.55,1.01);

    hPad->GetXaxis()->SetTitle("#phi_{SC} (deg)");
    hPad->GetYaxis()->SetTitle("Relative E/p scale");
    hPad->GetXaxis()->SetTitleOffset(0.8);
    hPad->GetYaxis()->SetTitleSize(0.05);
    hPad->GetYaxis()->SetLabelSize(0.05);
   
    if (mod > 1) hPad->GetYaxis()->SetRangeUser(0.95, 1.02);
 
    g_EoP_MC[mod] -> SetMarkerStyle(20);
    g_EoP_MC[mod] -> SetMarkerSize(1.);
    g_EoP_MC[mod] -> SetMarkerColor(kRed); 
    //    g_EoP_MC[mod] -> Draw("PL");
    g_EoC_MC[mod] -> SetMarkerStyle(20);
    g_EoC_MC[mod] -> SetMarkerSize(1.);
    g_EoC_MC[mod] -> SetMarkerColor(kRed+2); 
    g_EoC_MC[mod] -> Draw("PL");
  
    g_EoP_Data[mod] -> SetMarkerStyle(20);
    g_EoP_Data[mod] -> SetMarkerSize(1.);
    g_EoP_Data[mod] -> SetMarkerColor(kGreen); 
    //    g_EoP_Data[mod] -> Draw("PL");
    g_EoC_Data[mod] -> SetMarkerStyle(20);
    g_EoC_Data[mod] -> SetMarkerSize(1.);
    g_EoC_Data[mod] -> SetMarkerColor(kGreen+2); 
    g_EoC_Data[mod] -> Draw("PL");
  
    tl[mod] = new TLegend(0.60,0.15,0.89,0.45);
    tl[mod] -> SetFillColor(0);
    tl[mod] -> AddEntry(g_EoC_MC[mod],"MC","PL");
    tl[mod] -> AddEntry(g_EoC_Data[mod],"DATA","PL");
    //     tl[mod] -> AddEntry(g_EoP_MC[mod],"MC uncorrected","PL");
    //     tl[mod] -> AddEntry(g_EoC_MC[mod],"MC corrected","PL");
    //     tl[mod] -> AddEntry(g_EoP_Data[mod],"Data uncorrected","PL");
    //     tl[mod] -> AddEntry(g_EoC_Data[mod],"Data corrected","PL");
    tl[mod] -> Draw();

    c_g_fit[mod]->cd(2);
    gPad->SetGrid();
    TH1F *hPad2 = (TH1F*)gPad->DrawFrame(-0.55,0.98,0.55,1.02);
    hPad2->GetXaxis()->SetTitle("#phi_{SC} (deg)");
    hPad2->GetYaxis()->SetTitle("data/MC ratio");
    hPad2->GetXaxis()->SetTitleOffset(0.8);
    hPad2->GetYaxis()->SetTitleSize(0.05);
    hPad2->GetYaxis()->SetLabelSize(0.05);

    g_ratio_uncorr[mod]->SetLineColor(kBlue);
    g_ratio_uncorr[mod]->SetMarkerColor(kBlue);
    g_ratio_uncorr[mod]->SetMarkerStyle(20);
    g_ratio_uncorr[mod]->SetMarkerSize(0.7);
    //g_ratio_uncorr[mod]->Draw("PL");
    
    g_ratio_corr[mod]->SetLineColor(kBlue+3);
    g_ratio_corr[mod]->SetMarkerColor(kBlue+3);
    g_ratio_corr[mod]->SetMarkerStyle(20);
    g_ratio_corr[mod]->SetMarkerSize(0.7);
    g_ratio_corr[mod]->Draw("PL");
   
    tlr[mod] = new TLegend(0.60,0.15,0.89,0.45);
    tlr[mod] -> SetFillColor(0);
    tlr[mod] -> AddEntry(g_ratio_uncorr[mod],"uncorrected","PL");
    tlr[mod] -> AddEntry(g_ratio_corr[mod],"corrected","PL");
    //tlr[mod] -> Draw();

//     g_ratio_MC[mod]->SetLineColor(kBlue);
//     g_ratio_MC[mod]->SetMarkerColor(kBlue);
//     g_ratio_MC[mod]->SetMarkerStyle(20);
//     g_ratio_MC[mod]->SetMarkerSize(0.7);
//     g_ratio_MC[mod]->Draw("PL");
    
//     g_ratio_Data[mod]->SetLineColor(kBlue+3);
//     g_ratio_Data[mod]->SetMarkerColor(kBlue+3);
//     g_ratio_Data[mod]->SetMarkerStyle(20);
//     g_ratio_Data[mod]->SetMarkerSize(0.7);
//     g_ratio_Data[mod]->Draw("PL");
   
//     tlr[mod] = new TLegend(0.60,0.15,0.89,0.45);
//     tlr[mod] -> SetFillColor(0);
//     tlr[mod] -> AddEntry(g_ratio_MC[mod],"MC","PL");
//     tlr[mod] -> AddEntry(g_ratio_Data[mod],"DATA","PL");
//     tlr[mod] -> Draw();
  }

  
  TFile fout(outfilename,"recreate");
  for(int mod=0;mod<4;mod++) {
    g_EoP_MC[mod]->Write();
    g_EoC_MC[mod]->Write();
    g_EoP_Data[mod]->Write();
    g_EoC_Data[mod]->Write();
    g_ratio_MC[mod]->Write();
    g_ratio_Data[mod]->Write();
    g_ratio_uncorr[mod]->Write();
    g_ratio_corr[mod]->Write();
    spread_EoP_MC[mod]->Write();
    spread_EoC_MC[mod]->Write();
    spread_EoP_Data[mod]->Write();
    spread_EoC_Data[mod]->Write();
  }

  fout.Close();
  }
