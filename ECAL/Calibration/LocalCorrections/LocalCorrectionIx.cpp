// per compilare: g++ -Wall -o LocalCorrectionEE `root-config --cflags --glibs` LocalCorrectionEE.cpp

#include "histoFunc.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>


//derive f(local eta) correction in MC, with 4 different f for each module.
//fold crystal in eta into one single module


int EERing(flota eta){
  int myring = -1;
  if (fabs(eta) > 1.56 && fabs(eta)<=1.75) myring = 0;
  if (fabs(eta) > 1.75 && fabs(eta)<=2.00) myring = 1;
  if (fabs(eta) > 2.00 && fabs(eta)<=2.25) myring = 2;
  if (fabs(eta) > 2.25 && fabs(eta)<=2.50) myring = 3;
  return(myring);
}


int main(int argc, char** argv)
{

  gROOT->SetStyle("Plain");
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.04);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.04);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.04);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.04);
  gROOT->ForceStyle();

  float xtalWidth = 0.01745329;

  //---- variables for selections
  float etaMin = 1.56;
  float etaMax = 2.50;
  
  float r9min = 0.0 ;
  float r9max = 999999. ;  
  
  bool usePUweights = true;
  bool useR9weights = false;
  
  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"GraphsLocalEta_regression_EE.root");

  //--- weights for MC
  TFile weightsFile("weights/PUweights_2011_0100_73500_WJetsToLL_Fall11_S6.root","READ"); // stessi pesi usati per analisi vertici Hgg  
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();
 
  // --- NTUPLES
  TChain *ntu_MC = new TChain("ntu");
  TChain *ntu_Data = new TChain("ntu");
    
  //---- MC fall 2011
  ntu_MC->Add("../NTUPLES/Fall11/WZAnalysis/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11_All_new.root");
   
  //---- DATA
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B_new.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B_new.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_SingleElectron_Run2011A-WElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B_new.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/SingleElectron_Run2011A-WElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B_new.root");
  ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_SingleElectron_Run2011B-WElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B_new.root");
  
  std::cout << "     MC    : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;

  //---- Observables
  int npu;
  float PV_z;
  int isEEDeeGap;
  float EoP, scEta, scPhi;
  float scE3x3, scE5x5, scE, scE_regression;  
  float R9;
  int charge;
  float xC_DK, yC_DK;

  //---- Set branch addresses for MC  
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("PV_z", &PV_z);
  ntu_MC->SetBranchAddress("ele1_isEE", &isEEDeeGap);
  ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_MC->SetBranchAddress("ele1_scE",  &scE);
  ntu_MC->SetBranchAddress("ele1_xC_DK",&xC_DK); 
  ntu_MC->SetBranchAddress("ele1_yC_DK",&yC_DK); 


  //---- Set branch addresses for Data
  ntu_Data->SetBranchAddress("PV_z", &PV_z);
  ntu_Data->SetBranchAddress("ele1_isEE", &isEEDeeGap);
  ntu_Data->SetBranchAddress("ele1_scEta", &scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_Data->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_Data->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_Data->SetBranchAddress("ele1_scE", &scE);
  ntu_Data->SetBranchAddress("ele1_scE_regression", &scE_regression);
  ntu_Data->SetBranchAddress("ele1_charge", &charge);
  ntu_Data->SetBranchAddress("ele1_xC_DK",&xC_DK); 
  ntu_Data->SetBranchAddress("ele1_yC_DK",&yC_DK); 

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

  //-- book templates
  TH1F* h_template_MC[Ntempl];
  TH1F* h_template_Data[Ntempl];
  for(unsigned int i = 0; i < Ntempl; ++i){
    char histoName[100];
    sprintf(histoName, "template_MC_%d", i);
    h_template_MC[i] = new TH1F(histoName, "", 1200, 0., 3.);   
    sprintf(histoName, "template_Data_%d", i);
    h_template_Data[i] = new TH1F(histoName, "", 1200, 0., 3.);   
  }


  //******************************************************************************************
  //*************************************** MC  ********************************************** 
  std::cout << "Loop on MC events ... " << std::endl; 
 

  //---- loop on MC, make refernce and fit dist
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    
    if( entry%500000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //if (entry>500000) break;

    ntu_MC->GetEntry(entry);
    

    R9 = scE3x3/scE;
        
    float ww = 1 ;
    if (usePUweights) ww = w[npu];
    
    //-- eta or R9 cuts
    if ( fabs(scEta) < etaMin ) continue;
    if ( fabs(scEta) > etaMax ) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 

    //-- remove Dee gaps 
    if ( isEEDeeGap ) continue;

    //-- fill template for each mod
    int mod = EERing(scEta);
    if (mod == -1 ) continue;
    h_template_MC[mod]-> Fill(EoP,ww);
    
    //-- fill MC histos in x bins
    float locX = xC_DK/(xtalWidth/2)+0.5; 
    if ( locX < 0 || locX > 1) continue;
        
    int bin = nBins * (locX); 
    if (bin>nBins-1 || bin < 0 ) {
      std::cout << "Error in bins: " << bin << " " << scLocalEta << std::endl;
      continue;
    }

    float correction = 1;

    (h_EoP_MC[bin][mod]) -> Fill(EoP,ww);
    (h_EoC_MC[bin][mod]) -> Fill(EoP*correction,ww);

  }

  
  //******************************************************************************************
  //*************************************** DATA ********************************************** 
  std::cout << "Loop on Data events ..." << std::endl; 
  //--- loop on data
  for(int entry = 0; entry < ntu_Data->GetEntries(); ++entry) {
    if( entry%500000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //if (entry>500000) break;

    ntu_Data->GetEntry(entry);

    R9 = scE3x3/scE;

    //-- eta or R9 cuts
    if ( fabs(scEta) < etaMin ) continue;
    if ( fabs(scEta) > etaMax ) continue;
    if ( R9 < r9min || R9 > r9max ) continue; 
    
    //-- remove Dee gaps 
    if ( isEEDeeGap ) continue;

    //-- fill template for each mod
    int mod = EERing(scEta);
    if (mod == -1 ) continue;
    h_template_Data[mod]-> Fill(EoP);

    //-- fill data histos in X bins
    float locX = xC_DK/(xtalWidth/2)+0.5; 
    if ( locX < 0 || locX > 1) continue;
        
    int bin = nBins * (locX); 
    if (bin>nBins-1 || bin < 0 ) {
      std::cout << "Error in bins: " << bin << " " << scLocalEta << std::endl;
      continue;
    }

    float correction = scE_regression/scE;
    
    (h_EoP_Data[bin][mod]) -> Fill(EoP);
    (h_EoC_Data[bin][mod]) -> Fill(EoP*correction);
    
  }
  
  
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
    
    h_template_MC[mod]   -> Rebin(rebin);
    h_template_Data[mod] -> Rebin(rebin);
 
  }
  
  

  /////////////  MC ////////////////////
    
  for(int mod=0;mod<4;mod++){

    histoFunc *templateHistoFuncMC   = new histoFunc(h_template_MC[mod]);
    histoFunc *templateHistoFuncData = new histoFunc(h_template_Data[mod]);

    for(unsigned int i = 0; i < nBins; ++i)
      {
	std::cout << "***** Fitting :  mod: " <<mod <<"  bin:  "<< i << std::endl; 

	h_EoP_MC[i][mod] -> Rebin(rebin);    
	h_EoC_MC[i][mod] -> Rebin(rebin);    
	h_EoP_Data[i][mod] -> Rebin(rebin);    
	h_EoC_Data[i][mod] -> Rebin(rebin);    



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
    

	TFitResultPtr	r = h_EoP_MC[i][mod] -> Fit("templateFuncMC", "SMRQL+");
	int fitStatus = r;
	
	g_EoP_MC[mod] -> SetPoint(i, xval , 1./templateFuncMC->GetParameter(1));
	g_EoP_MC[mod] -> SetPointError(i, 0., templateFuncMC->GetParError(1));
	//    if ( templateFunc->GetParError(1) < 0.003) g_EoP -> SetPointError(i, 0., 0.003);
	//    cout  << " ***** " <<  1./templateFunc->GetParameter(1) << " " << templateFunc->GetParError(1) << endl; 
	
	float scaleMC  = 1./templateFuncMC->GetParameter(1);
	float escaleMC = templateFuncMC->GetParError(1)/scaleMC/scaleMC;


	//--- corrected MC  
	xNorm = h_EoC_MC[i][mod]->GetEntries()/h_template_MC[mod]->GetEntries() *
	        h_EoC_MC[i][mod]->GetBinWidth(1)/h_template_MC[mod]->GetBinWidth(1); 

	templateFuncMC -> FixParameter(0, xNorm);
	templateFuncMC -> SetParameter(1, 1.05 );
	templateFuncMC -> FixParameter(2, 0.);
    

	h_EoC_MC[i][mod] -> Fit("templateFuncMC", "MRQLS+");
	g_EoC_MC[mod] -> SetPoint(i, xval , 1./templateFuncMC->GetParameter(1));
	g_EoC_MC[mod] -> SetPointError(i, 0., templateFuncMC->GetParError(1));
    

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
    

	h_EoP_Data[i][mod] -> Fit("templateFuncData", "MRQLS+");
	g_EoP_Data[mod] -> SetPoint(i, xval , 1./templateFuncData->GetParameter(1));
	g_EoP_Data[mod] -> SetPointError(i, 0., templateFuncData->GetParError(1));
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
    

	h_EoC_Data[i][mod] -> Fit("templateFuncData", "MRQLS+");
	g_EoC_Data[mod] -> SetPoint(i, xval, 1./templateFuncData->GetParameter(1));
	g_EoC_Data[mod] -> SetPointError(i, 0., templateFuncData->GetParError(1));

	float scaleDAcorr  = 1./templateFuncData->GetParameter(1);
	float escaleDAcorr = templateFuncData->GetParError(1)/scaleDAcorr/scaleDAcorr;

	//--- ratio finalization MC corr/uncorr
	float ratioMC = scaleMC/scaleMCcorr;
	float eratioMC = ratioMC*sqrt(pow(escaleMC/scaleMC,2) + pow(escaleMCcorr/scaleMCcorr,2)); 
	g_ratio_MC[mod] -> SetPoint(i, xval, ratioMC);
	g_ratio_MC[mod] -> SetPointError(i, 0., eratioMC);

	//--- ratio finalization DATA corr/uncorr
	float ratioDA  = scaleDA/scaleDAcorr;
	float eratioDA = ratioDA*sqrt(pow(escaleDA/scaleDA,2) + pow(escaleDAcorr/scaleDAcorr,2)); 
	g_ratio_Data[mod] -> SetPoint(i, xval , ratioDA);
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

  /*
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
    TH1F *hPad = (TH1F*)gPad->DrawFrame(-0.55,0.96,0.55,1.03);
    hPad->GetXaxis()->SetTitle("#eta_{SC} (deg)");
    hPad->GetYaxis()->SetTitle("Relative E/p scale");
    hPad->GetXaxis()->SetTitleOffset(0.8);
    hPad->GetYaxis()->SetTitleSize(0.05);
    hPad->GetYaxis()->SetLabelSize(0.05);


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
    //g_EoP_Data[mod] -> Draw("PL");
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
    TH1F *hPad2 = (TH1F*)gPad->DrawFrame(-0.55,0.97,0.55,1.03);
    hPad2->GetXaxis()->SetTitle("#eta_{SC} (deg)");
    hPad2->GetYaxis()->SetTitle("data/MC ratio");
    hPad2->GetXaxis()->SetTitleOffset(0.8);
    hPad2->GetYaxis()->SetTitleSize(0.05);
    hPad2->GetYaxis()->SetLabelSize(0.05);
  
    g_ratio_uncorr[mod]->SetLineColor(kBlue);
    g_ratio_uncorr[mod]->SetMarkerColor(kBlue);
    g_ratio_uncorr[mod]->SetMarkerStyle(20);
    g_ratio_uncorr[mod]->SetMarkerSize(0.7);
    g_ratio_uncorr[mod]->Draw("PL");
    
    g_ratio_corr[mod]->SetLineColor(kBlue+3);
    g_ratio_corr[mod]->SetMarkerColor(kBlue+3);
    g_ratio_corr[mod]->SetMarkerStyle(20);
    g_ratio_corr[mod]->SetMarkerSize(0.7);
    g_ratio_corr[mod]->Draw("PL");
   
    tlr[mod] = new TLegend(0.60,0.15,0.89,0.35);
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

  */
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
    
  }
  fout.Close();
}
