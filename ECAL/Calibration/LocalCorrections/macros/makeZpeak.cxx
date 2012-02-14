typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > PtEtaPhiELorentzVector;

using namespace RooFit;

//****Effective sigma from fitted func ****************************************************************
void evalEffectiveSigma( RooAddPdf* pdf , RooRealVar*  mass, double &ml, double &mh, double &sigmaeff){

  double testmass = 90.;
  double center = testmass-10.0;
  
  double cdfhi = 0;
  double cdflo = 0;
  double mlmin = 0.0;
  double mhmin = 0.0;
  double step  = 0.05;

  double minwidth = 999.0;

  RooAbsReal *cdf = 0;
  double mlow = -999;
  double mhigh = -999;
  
  for (int i=0; i<400; ++i) {
    if (i%100==0) cout <<  i << endl;
    mlow = center+i*step;
    mass->setRange("signal",60,mlow); 
    cdf = pdf->createIntegral(*mass, NormSet(*mass),Range("signal"));
    cdflo = cdf->getVal();
    for (int j=i; j<400; ++j) {
      mhigh = center+j*step;
      mass->setRange("signal",60,mhigh); 
      cdf   = pdf->createIntegral(*mass, NormSet(*mass),Range("signal"));
      cdfhi = cdf->getVal();
      if ( (cdfhi-cdflo)>0.684 ) {
	if ( (mhigh-mlow)<minwidth) {
	  minwidth = mhigh-mlow;
	  mlmin = mlow;
	  mhmin = mhigh;
	}
	break;
      }
    }
  }
 
  sigmaeff = minwidth/2.0;
  ml = mlmin;
  mh = mhmin;
  cout << "Mmin = " << mlmin << "  Mmax = " << mhmin << "  effective sigma = " << sigmaeff << endl;
  
  // return (sigmaeff);
  return;
}
//*****************************************************************************************************


//****Effective sigma from histogram ****************************************************************
void evalEffectiveSigmaFromHisto( TH1F *h , double &ml, double &mh, double &sigmaeff){

  double testmass = 90.;
  double center = testmass-6.0;
  
  double cdfhi = 0;
  double cdflo = 0;
  double mlmin = 0.0;
  double mhmin = 0.0;
  double step  = 0.02;

  double minwidth = 999.0;

  double mlow = -999;
  double mhigh = -999;
  int binlow, binhigh;
  int nbins = h->GetNbinsX();

  for (int i=0; i<600; ++i) {
    if (i%100==0) cout <<  i << endl;
    mlow = center+i*step;
    binlow = h ->FindBin(mlow);
    cdflo = h->Integral(1, binlow)/h->Integral(1, nbins);
    for (int j=i; j<600; ++j) {
      mhigh = center+j*step;
      binhigh = h ->FindBin(mhigh);
      cdfhi = h->Integral(1, binhigh)/h->Integral(1, nbins);
      if ( (cdfhi-cdflo)>0.684 ) {
	if ( (mhigh-mlow)<minwidth) {
	  minwidth = mhigh-mlow;
	  mlmin = mlow;
	  mhmin = mhigh;
	}
	break;
      }
    }
  }
  
  sigmaeff = minwidth/2.0;
  ml = mlmin;
  mh = mhmin;
  cout << "Mmin = " << mlmin << "  Mmax = " << mhmin << "  effective sigma = " << sigmaeff << endl;
  
  // return (sigmaeff);
  return;
}
//*****************************************************************************************************

int modId(float eta){
  int ieta  = fabs(eta)/0.01745329;
  int im = 0;
  if (ieta>25) im = 1;
  if (ieta>45) im = 2;
  if (ieta>65) im = 3;
  return (im);
}

bool IsEtaGap(float eta){
  int ieta  = fabs(eta)/0.01745329;
  float feta = fabs(ieta);
  if( fabs(feta - 0 )<2) return true;
  if( fabs(feta - 25)<2) return true;
  if( fabs(feta - 45)<2) return true;
  if( fabs(feta - 65)<2) return true;
  if( fabs(feta - 85)<2) return true;
  return false;
}


void makeZPeak()
{
  RooMsgService::instance().Print() ;
  RooMsgService::instance().setStreamStatus(1,0);


  //
  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);
  
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.03);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.03);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.03);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);
  gROOT->ForceStyle();

  float xtalWidth = 0.01745329;
  float zmass = 91.187;

  //---- variables for selections
  float etaMin = 0.;
  float etaMax = 2.5;
  
  float r9min = 0. ;
  float r9max = 999.94;
  
  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"dataZpeak_allR9_regression.root");

  //----

  // Get trees
  TChain *ntu_Data = new TChain("ntu");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
  ntu_Data->Add("../NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");

    
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;
  std::cout << std::endl;
  
  // Set branch addresses  
  int runId,isZ;
  int PV_n;
  int isEB,isEB2, isEBEEGap,isEBEEGap2 ;
  float mZ, mZc;
  
  float p, p2;
  float scE,scE2;
  float scE3x3, scE3x3_2;
  float scE_regression, scE2_regression;
  float scERaw_PUcleaned, scERaw2_PUcleaned; 
  float fCorrPUcleaned, fCorrPUcleaned2;
  float crackCorr, crackCorr2;
  float scEta,scEta2;
  float scPhi,scPhi2;
  float eta,eta2;
  float phi,phi2;
  float scLocalEta, scLocalEta2;
  float scLocalPhi, scLocalPhi2;
  int charge , charge2;

  ntu_Data->SetBranchAddress("isZ",              &isZ);
  ntu_Data->SetBranchAddress("ele1_isEB",        &isEB);
  ntu_Data->SetBranchAddress("ele2_isEB",        &isEB2);
  ntu_Data->SetBranchAddress("ele1_isEBEEGap",   &isEBEEGap);
  ntu_Data->SetBranchAddress("ele2_isEBEEGap",   &isEBEEGap2);
  ntu_Data->SetBranchAddress("ele1_scE",         &scE);
  ntu_Data->SetBranchAddress("ele2_scE",         &scE2);
  ntu_Data->SetBranchAddress("ele1_scE_regression",  &scE_regression);
  ntu_Data->SetBranchAddress("ele2_scE_regression",  &scE2_regression);
  ntu_Data->SetBranchAddress("ele1_scERaw_PUcleaned",  &scERaw_PUcleaned);
  ntu_Data->SetBranchAddress("ele2_scERaw_PUcleaned",  &scERaw2_PUcleaned);
  ntu_Data->SetBranchAddress("ele1_fCorrection_PUcleaned", &fCorrPUcleaned);  
  ntu_Data->SetBranchAddress("ele2_fCorrection_PUcleaned", &fCorrPUcleaned2);  
  ntu_Data->SetBranchAddress("ele1_scEta",       &scEta);
  ntu_Data->SetBranchAddress("ele2_scEta",       &scEta2);
  ntu_Data->SetBranchAddress("ele1_scPhi",       &scPhi);
  ntu_Data->SetBranchAddress("ele2_scPhi",       &scPhi2);
  ntu_Data->SetBranchAddress("ele1_eta",         &eta);
  ntu_Data->SetBranchAddress("ele2_eta",         &eta2);
  ntu_Data->SetBranchAddress("ele1_phi",         &phi);
  ntu_Data->SetBranchAddress("ele2_phi",         &phi2);
  ntu_Data->SetBranchAddress("ele1_scLocalPhi",  &scLocalPhi); 
  ntu_Data->SetBranchAddress("ele1_scLocalEta",  &scLocalEta); 
  ntu_Data->SetBranchAddress("ele1_scCrackCorr", &crackCorr); 
  ntu_Data->SetBranchAddress("ele2_scLocalPhi",  &scLocalPhi2); 
  ntu_Data->SetBranchAddress("ele2_scLocalEta",  &scLocalEta2); 
  ntu_Data->SetBranchAddress("ele2_scCrackCorr", &crackCorr2); 
  ntu_Data->SetBranchAddress("ele1_charge", &charge); 
  ntu_Data->SetBranchAddress("ele2_charge", &charge2); 
  

  // Define histograms
  TH1F* h_mZ_EBEB = new TH1F("h_mZ_EBEB", "",2500,65.,115.);
  TH1F* h_mZ_EEEE = new TH1F("h_mZ_EEEE", "",2500,65.,115.);
  TH1F* h_mZ_EBEE = new TH1F("h_mZ_EBEE", "",2500,65.,115.);

  TH1F* h_mZc_EBEB = new TH1F("h_mZc_EBEB", "",2500,65.,115.);
  TH1F* h_mZc_EEEE = new TH1F("h_mZc_EEEE", "",2500,65.,115.);  
  TH1F* h_mZc_EBEE = new TH1F("h_mZc_EBEE", "",2500,65.,115.);
  
  TProfile *mZ_vs_localEta_EBEB = new TProfile("mZ_vs_localEta_EBEB", "",20,0,1,0.,2.);
  mZ_vs_localEta_EBEB->SetMarkerStyle(20);
  mZ_vs_localEta_EBEB->SetMarkerColor(kRed);
  mZ_vs_localEta_EBEB->SetLineColor(kRed);
  mZ_vs_localEta_EBEB-> GetXaxis()->SetTitle("#eta_{SC} (deg)");
  mZ_vs_localEta_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZc_vs_localEta_EBEB = new TProfile("mZc_vs_localEta_EBEB", "",20,0,1,0.,2.);
  mZc_vs_localEta_EBEB->SetMarkerStyle(20);
  mZc_vs_localEta_EBEB->SetMarkerColor(kGreen+2);
  mZc_vs_localEta_EBEB->SetLineColor(kGreen+2);
  mZc_vs_localEta_EBEB-> GetXaxis()->SetTitle("#eta_{SC} (deg)");
  mZc_vs_localEta_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZ_vs_localPhi_EBEB = new TProfile("mZ_vs_localPhi_EBEB", "",20,0,1,0.,2.);
  mZ_vs_localPhi_EBEB->SetMarkerStyle(20);
  mZ_vs_localPhi_EBEB->SetMarkerColor(kRed);
  mZ_vs_localPhi_EBEB->SetLineColor(kRed);
  mZ_vs_localPhi_EBEB-> GetXaxis()->SetTitle("#phi_{SC} (deg)");
  mZ_vs_localPhi_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZc_vs_localPhi_EBEB = new TProfile("mZc_vs_localPhi_EBEB", "",20,0,1,0.,2.);
  mZc_vs_localPhi_EBEB->SetMarkerStyle(20);
  mZc_vs_localPhi_EBEB->SetMarkerColor(kGreen+2);
  mZc_vs_localPhi_EBEB->SetLineColor(kGreen+2);
  mZc_vs_localPhi_EBEB-> GetXaxis()->SetTitle("#phi_{SC} (deg)");
  mZc_vs_localPhi_EBEB-> GetYaxis()->SetTitle("< M_{ee}> /M_{Z}");

  TProfile *mZ_vs_Eta_EBEB = new TProfile("mZ_vs_Eta_EBEB", "",50,-2.5,2.5,0.,2.);
  mZ_vs_Eta_EBEB->SetMarkerStyle(20);
  mZ_vs_Eta_EBEB->SetMarkerColor(kRed);
  mZ_vs_Eta_EBEB->SetLineColor(kRed);
  mZ_vs_Eta_EBEB-> GetXaxis()->SetTitle("#eta_{SC}");
  mZ_vs_Eta_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZc_vs_Eta_EBEB = new TProfile("mZc_vs_Eta_EBEB", "",50,-2.5,2.5,0.,2.);
  mZc_vs_Eta_EBEB->SetMarkerStyle(20);
  mZc_vs_Eta_EBEB->SetMarkerColor(kGreen+2);
  mZc_vs_Eta_EBEB->SetLineColor(kGreen+2);
  mZc_vs_Eta_EBEB-> GetXaxis()->SetTitle("#eta_{SC}");
  mZc_vs_Eta_EBEB-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZ_vs_Eta_EEEE = new TProfile("mZ_vs_Eta_EEEE", "",50,-2.5,2.5,0.,2.);
  mZ_vs_Eta_EEEE->SetMarkerStyle(20);
  mZ_vs_Eta_EEEE->SetMarkerColor(kRed);
  mZ_vs_Eta_EEEE->SetLineColor(kRed);
  mZ_vs_Eta_EEEE-> GetXaxis()->SetTitle("#eta_{SC}");
  mZ_vs_Eta_EEEE-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TProfile *mZc_vs_Eta_EEEE = new TProfile("mZc_vs_Eta_EEEE", "",50,-2.5,2.5,0.,2.);
  mZc_vs_Eta_EEEE->SetMarkerStyle(20);
  mZc_vs_Eta_EEEE->SetMarkerColor(kGreen+2);
  mZc_vs_Eta_EEEE->SetLineColor(kGreen+2);
  mZc_vs_Eta_EEEE-> GetXaxis()->SetTitle("#eta_{SC}");
  mZc_vs_Eta_EEEE-> GetYaxis()->SetTitle("< M_{ee} > / M_{Z}");

  TH1F *hspread = new TH1F("hspread","spread",800,0.8,1.2);
  TH1F *hspreadc = new TH1F("hspreadc","spread",800,0.8,1.2);

  
  // Loop over entries
  int nEntries = ntu_Data -> GetEntriesFast();
  //nEntries = 20000;
  for(int ientry = 0; ientry < nEntries; ++ientry)
  {
    if( ientry%10000 == 0 ) std::cout << ">>>>> Reading entry : " << ientry << "\r" << std::flush;

    ntu_Data -> GetEntry(ientry); 
    
    if( isZ == 0 ) continue;
    
    // eta cut 
    if ( isEBEEGap || isEBEEGap2) continue;
    if ( fabs(scEta)  < etaMin ) continue;
    if ( fabs(scEta2) < etaMin ) continue;
    if ( fabs(scEta)  > etaMax ) continue;
    if ( fabs(scEta2) > etaMax ) continue;

    // remove eta gaps
    //     if ( IsEtaGap(scEta) || IsEtaGap(scEta2) ) continue;
       
    // //    remove phi cracks
    //     float myphi = (scPhi+3.1415926536)/xtalWidth;
    //     float modphi = (int)myphi%20;
    //     if (fabs(modphi-10)<2.) continue;
    //     float myphi2 = (scPhi2+3.1415926536)/xtalWidth;
    //     float modphi2 = (int)myphi2%20;
    //     if (fabs(modphi2-10)<2.) continue;
    
    // select on r9
    float R9 = scE3x3/scE;
    float R9_2 = scE3x3_2/scE2;

    if ( R9   < r9min || R9 > r9max ) continue; 
    if ( R9_2 < r9min || R9_2 > r9max ) continue; 

    PtEtaPhiELorentzVector sc1(scE *sin(2*atan(exp(-1.*eta))), eta, phi, scE);
    PtEtaPhiELorentzVector sc2(scE2*sin(2*atan(exp(-1.*eta2))),eta2,phi2,scE2);
    mZ = (sc1+sc2).mass();
 
    //corrected energy : use DK corrections
    //float scEc  = scE*scLocalCorr_DK;
    //float scE2c = scE2*scLocalCorr_DK2;
    
    // corrected energy: use regression
    float scEc  = scE_regression;
    float scE2c = scE2_regression;

    // corrected energy: use PU cleaning
    // if (scERaw_PUcleaned< 0 || scERaw2_PUcleaned< 0) continue;
    // float scEc  = scERaw_PUcleaned*fCorrPUcleaned*crackCorr;
    // float scE2c = scERaw2_PUcleaned*fCorrPUcleaned2*crackCorr2;
    
    PtEtaPhiELorentzVector sc1c(scEc *sin(2*atan(exp(-1.*eta))), eta, phi, scEc);
    PtEtaPhiELorentzVector sc2c(scE2c*sin(2*atan(exp(-1.*eta2))),eta2,phi2,scE2c);
    mZc = (sc1c+sc2c).mass();
    
    // both in EB
    if( (isEB == 1) && (isEB2 == 1) )      {
      h_mZ_EBEB -> Fill(mZ);
      h_mZc_EBEB -> Fill(mZc);
	
      mZ_vs_Eta_EBEB->Fill(scEta, mZ/zmass);
      mZ_vs_Eta_EBEB->Fill(scEta2, mZ/zmass);
      
      mZc_vs_Eta_EBEB->Fill(scEta, mZc/zmass);
      mZc_vs_Eta_EBEB->Fill(scEta2, mZc/zmass);
      
      mZ_vs_localEta_EBEB  -> Fill(scLocalEta+0.5,mZ/zmass);
      mZc_vs_localEta_EBEB -> Fill(scLocalEta+0.5,mZc/zmass);
      mZ_vs_localEta_EBEB  -> Fill(scLocalEta2+0.5,mZ/zmass);
      mZc_vs_localEta_EBEB -> Fill(scLocalEta2+0.5,mZc/zmass);
      
      if ( (scEta*charge) > 0  ){
	mZ_vs_localPhi_EBEB  -> Fill(scLocalPhi+0.5,mZ/zmass);
	mZc_vs_localPhi_EBEB -> Fill(scLocalPhi+0.5,mZc/zmass);
      }
      
      if ( (scEta2*charge2) > 0){
	mZ_vs_localPhi_EBEB -> Fill(scLocalPhi2+0.5,mZ/zmass);
	mZc_vs_localPhi_EBEB -> Fill(scLocalPhi2+0.5,mZc/zmass);
      }
    }
    
    // one in EB, one in EE
    if( ( isEB && !isEB2) ||  ( !isEB && isEB2 )  )
      //if( !isEB || !isEB2 )// at least one in EE
      {
	h_mZ_EBEE  -> Fill(mZ);
	h_mZc_EBEE -> Fill(mZc);
      }      

    // both in EE
    if( ( !isEB && !isEB2)  ) {
      h_mZ_EEEE  -> Fill(mZ);
      h_mZc_EEEE -> Fill(mZc);
      
      mZ_vs_Eta_EBEB->Fill(scEta, mZ/zmass);
      mZ_vs_Eta_EBEB->Fill(scEta2, mZ/zmass);
      
      mZc_vs_Eta_EBEB->Fill(scEta, mZc/zmass);
      mZc_vs_Eta_EBEB->Fill(scEta2, mZc/zmass);
      
    }      
    
  }
  
  
  for (int i =0 ; i< mZ_vs_localEta_EBEB->GetNbinsX(); i++){
    hspread  ->Fill(mZ_vs_localEta_EBEB->GetBinContent(i+1));
    hspreadc ->Fill(mZc_vs_localEta_EBEB->GetBinContent(i+1));
  }
  
  
  
  //COMPUTE EFFECTIVE SIGMA FROM HISTOS
  double mlEBEB = 0, mhEBEB = 0, effsigmaEBEB = 0;
  double mlcEBEB = 0, mhcEBEB = 0, effsigmacEBEB = 0;
  double mlEBEE = 0, mhEBEE = 0, effsigmaEBEE = 0;
  double mlcEBEE = 0, mhcEBEE = 0, effsigmacEBEE = 0;
  double mlEEEE = 0, mhEEEE = 0, effsigmaEEEE = 0;
  double mlcEEEE = 0, mhcEEEE = 0, effsigmacEEEE = 0;

  evalEffectiveSigmaFromHisto(h_mZ_EBEB, mlEBEB, mhEBEB, effsigmaEBEB);
  evalEffectiveSigmaFromHisto(h_mZc_EBEB, mlcEBEB, mhcEBEB, effsigmacEBEB);
  evalEffectiveSigmaFromHisto(h_mZ_EBEE, mlEBEE, mhEBEE, effsigmaEBEE);
  evalEffectiveSigmaFromHisto(h_mZc_EBEE, mlcEBEE, mhcEBEE, effsigmacEBEE);
  evalEffectiveSigmaFromHisto(h_mZ_EEEE, mlEEEE, mhEEEE, effsigmaEEEE);
  evalEffectiveSigmaFromHisto(h_mZc_EEEE, mlcEEEE, mhcEEEE, effsigmacEEEE);
  

  int nre = 25;
  h_mZ_EBEB  -> Rebin(nre);
  h_mZc_EBEB -> Rebin(nre);
  h_mZ_EBEE  -> Rebin(nre);
  h_mZc_EBEE -> Rebin(nre);
  h_mZ_EEEE  -> Rebin(nre);
  h_mZc_EEEE -> Rebin(nre);


  TLegend   *tl = new TLegend(0.60,0.15,0.89,0.35);
  tl-> SetFillColor(0);
  tl-> AddEntry(mZ_vs_localEta_EBEB,"uncorrected","PL");
  tl-> AddEntry(mZc_vs_localEta_EBEB,"corrected","PL");
    
  TCanvas* c1 = new TCanvas("c1","c1",100,100,700,500);
  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  mZ_vs_localEta_EBEB->GetYaxis()->SetRangeUser(0.96,1.02);
  mZ_vs_localEta_EBEB->Draw("");
  mZc_vs_localEta_EBEB->Draw("same");
  tl-> Draw("same");

  TCanvas* c2 = new TCanvas("c2","c2",100,100,700,500);
  c2->cd();
  c2->SetGridx();
  c2->SetGridy();
  mZ_vs_localPhi_EBEB->GetYaxis()->SetRangeUser(0.96,1.02);
  mZ_vs_localPhi_EBEB->Draw("");
  mZc_vs_localPhi_EBEB->Draw("same");
  tl-> Draw("same");

  TCanvas* c3 = new TCanvas("c3","c3",100,100,700,500);
  c3->cd();
  c3->SetGridx();
  c3->SetGridy();
  mZ_vs_Eta_EBEB->GetYaxis()->SetRangeUser(0.96,1.02);
  mZ_vs_Eta_EBEB->Draw("");
  mZc_vs_Eta_EBEB->Draw("same");
  tl-> Draw("same");

  TCanvas* c4 = new TCanvas("c4","c4",100,100,700,500);
  c4->cd();
  c4->SetGridx();
  c4->SetGridy();
  mZ_vs_Eta_EEEE->GetYaxis()->SetRangeUser(0.96,1.02);
  mZ_vs_Eta_EEEE->Draw("");
  mZc_vs_Eta_EEEE->Draw("same");
  tl-> Draw("same");

  //------------------------
  //---------------- fitting
  
  RooRealVar  mass("mass","M(e^{+}e^{-})", 70.0, 110.0,"GeV/c^{2}");
  mass.setBins(10000) ;


  // Parameters for Crystal Ball Lineshape 
  RooRealVar  dm("#Delta m", "offset", 0.0, -5.0, 5.0,"GeV/c^{2}"); 
  RooRealVar  sigma("#sigma_{CB}","sigmaCB", 1.5,0.5,7.5,"GeV/c^{2}"); 
  RooRealVar  alpha("#alpha","alpha", 1.607,0.6,2.0); 
  RooRealVar  n("n","n", 6., 0.5, 50.0); 
  //alpha.setConstant();
  //n.setConstant();

  // Parameters for Breit-Wigner Distribution
  RooRealVar  MZ("M_{Z}", "M_{Z}", 91.188, 80.0, 100.0,"GeV/c^{2}"); 
  RooRealVar  Gamma("#Gamma", "#Gamma", 2.45, 2.0,3.0,"GeV/c^{2}"); 
  MZ.setConstant();
  Gamma.setConstant();

  // Exponential Background
  RooRealVar  bkgshape("bkgshape", "Backgroung Shape", -0.1,-1.0,0.0, "1/GeV/c^{2}");
  RooRealVar  frac("frac", "Signal Fraction", 1.0,0.0,1.0);
  frac.setConstant();
  bkgshape.setConstant();
  
  // Crystal Ball Lineshape
  RooCBShape     cb("cb", "Crystal Ball Lineshape", mass, dm,sigma, alpha, n);
  // Breit-Wigner Distribution
  RooBreitWigner bw("bw","A Breit-Wigner Distribution",mass,MZ,Gamma);
  
  // Convolution p.d.f. using numeric convolution operator based on Fourier Transforms
  RooFFTConvPdf bwcb("bwcb","convolution", mass, bw, cb);
  
  // Background  p.d.f.
  RooExponential bkg("bkg", "Backgroung Distribution", mass, bkgshape);
  
  
  // --- EB-EB
  TCanvas *cEBEB = new TCanvas("cEBEB","cEBEB",700,700);
  cEBEB->SetLeftMargin(0.15);

  RooDataHist hdataEBEB("hdataEBEB","hdataEBEB",RooArgSet(mass), h_mZ_EBEB);
  RooAddPdf   modelEBEB("modelEBEB", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rEBEB = modelEBEB.fitTo(hdataEBEB, Range(70., 110.), Save(),Verbose(0) );
  RooPlot* plotEBEB = mass.frame(Range(70,110));
  hdataEBEB.plotOn(plotEBEB);
  modelEBEB.plotOn(plotEBEB);
  modelEBEB.plotOn(plotEBEB, LineColor(kRed));
  modelEBEB.paramOn(plotEBEB, Layout(0.16,0.46,0.89));
  modelEBEB.plotOn(plotEBEB, LineColor(kRed));
  modelEBEB.plotOn(plotEBEB, Range(mlEBEB,mhEBEB) , LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3005) , VLines(), MoveToBack() );
  plotEBEB->SetTitleOffset(1.7,"Y");
  plotEBEB->Draw();
  TPaveText *pEBEB = (TPaveText*)cEBEB->FindObject("modelEBEB_paramBox");
  pEBEB->SetTextSize(0.025);
  pEBEB->SetLineColor(kRed);

  //evalEffectiveSigma(&modelEBEB, &mass, mlEBEB, mhEBEB, effsigmaEBEB);


 
  RooDataHist hdatacEBEB("hdatacEBEB","hdatacEBEB",RooArgSet(mass), h_mZc_EBEB);
  RooAddPdf   modelcEBEB("modelcEBEB", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rcEBEB = modelcEBEB.fitTo(hdatacEBEB, Range(70., 110.), Save(),Verbose(0));
  RooPlot* plotcEBEB = mass.frame(Range(70,110));
  hdatacEBEB.plotOn(plotcEBEB);
  modelcEBEB.plotOn(plotcEBEB, LineColor(kGreen+2));
  modelcEBEB.paramOn(plotcEBEB, Layout(0.16,0.46,0.63));
  modelcEBEB.plotOn(plotcEBEB, LineColor(kGreen+2));
  modelcEBEB.plotOn(plotcEBEB, Range(mlcEBEB,mhcEBEB) , LineColor(kGreen+2), DrawOption("F"), FillColor(kGreen+2), FillStyle(3004) , VLines(), MoveToBack() );
  plotcEBEB->SetTitleOffset(1.7,"Y");
  plotcEBEB->Draw("same");
  TPaveText *pcEBEB = (TPaveText*)cEBEB->FindObject("modelcEBEB_paramBox");
  pcEBEB->SetTextSize(0.025);
  pcEBEB->SetLineColor(kGreen+2);
  //evalEffectiveSigma(&modelcEBEB, &mass, mlcEBEB, mhcEBEB, effsigmacEBEB);


  // --- EB-EE
  TCanvas *cEBEE = new TCanvas("cEBEE","cEBEE",700,700);
  cEBEE->SetLeftMargin(0.15);

  RooDataHist hdataEBEE("hdataEBEE","hdataEBEE",RooArgSet(mass), h_mZ_EBEE);
  RooAddPdf   modelEBEE("modelEBEE", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rEBEE = modelEBEE.fitTo(hdataEBEE, Range(70., 110.), Save() ,Verbose(0));
  RooPlot* plotEBEE = mass.frame(Range(70,110));
  hdataEBEE.plotOn(plotEBEE);
  modelEBEE.plotOn(plotEBEE);
  modelEBEE.plotOn(plotEBEE, LineColor(kRed));
  modelEBEE.paramOn(plotEBEE, Layout(0.16,0.46,0.89));
  modelEBEE.plotOn(plotEBEE, LineColor(kRed));
  modelEBEE.plotOn(plotEBEE, Range(mlEBEE,mhEBEE) , LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3005) , VLines(), MoveToBack() );
  plotEBEE->SetTitleOffset(1.7,"Y");
  plotEBEE->Draw();
  TPaveText *pEBEE = (TPaveText*)cEBEE->FindObject("modelEBEE_paramBox");
  pEBEE->SetTextSize(0.025);
  pEBEE->SetLineColor(kRed);
  //evalEffectiveSigma(&modelEBEE, &mass, mlEBEE, mhEBEE, effsigmaEBEE);

  RooDataHist hdatacEBEE("hdatacEBEE","hdatacEBEE",RooArgSet(mass), h_mZc_EBEE);
  RooAddPdf   modelcEBEE("modelcEBEE", "Signal + Background", bwcb, bkg, frac);
  RooFitResult *rcEBEE = modelcEBEE.fitTo(hdatacEBEE, Range(70., 110.), Save(),Verbose(0));
  RooPlot* plotcEBEE = mass.frame(Range(70,110));
  hdatacEBEE.plotOn(plotcEBEE);
  modelcEBEE.plotOn(plotcEBEE, LineColor(kGreen+2));
  modelcEBEE.paramOn(plotcEBEE, Layout(0.16,0.46,0.63));
  modelcEBEE.plotOn(plotcEBEE, LineColor(kGreen+2));
  modelcEBEE.plotOn(plotcEBEE, Range(mlcEBEE,mhcEBEE) , LineColor(kGreen+2), DrawOption("F"), FillColor(kGreen+2), FillStyle(3004) , VLines(), MoveToBack() );  
  plotcEBEE->SetTitleOffset(1.7,"Y");
  plotcEBEE->Draw("same");
  TPaveText *pcEBEE = (TPaveText*)cEBEE->FindObject("modelcEBEE_paramBox");
  pcEBEE->SetTextSize(0.025);
  pcEBEE->SetLineColor(kGreen+2);
  //evalEffectiveSigma(&modelcEBEE, &mass, mlcEBEE, mhcEBEE, effsigmacEBEE);
   
  // --- EE-EE
  TCanvas *cEEEE = new TCanvas("cEEEE","cEEEE",700,700);
  cEEEE->SetLeftMargin(0.15);

  RooDataHist hdataEEEE("hdataEEEE","hdataEEEE",RooArgSet(mass), h_mZ_EEEE);
  RooAddPdf   modelEEEE("modelEEEE", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rEEEE = modelEEEE.fitTo(hdataEEEE, Range(70., 110.), Save() ,Verbose(0));
  RooPlot* plotEEEE = mass.frame(Range(70,110));
  hdataEEEE.plotOn(plotEEEE);
  modelEEEE.plotOn(plotEEEE);
  modelEEEE.plotOn(plotEEEE, LineColor(kRed));
  modelEEEE.paramOn(plotEEEE, Layout(0.16,0.46,0.89));
  modelEEEE.plotOn(plotEEEE, LineColor(kRed));
  modelEEEE.plotOn(plotEEEE, Range(mlEEEE,mhEEEE) , LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3005) , VLines(), MoveToBack() );
  plotEEEE->SetTitleOffset(1.7,"Y");
  plotEEEE->Draw();
  TPaveText *pEEEE = (TPaveText*)cEEEE->FindObject("modelEEEE_paramBox");
  pEEEE->SetTextSize(0.025);
  pEEEE->SetLineColor(kRed);
  //evalEffectiveSigma(&modelEEEE, &mass, mlEEEE, mhEEEE, effsigmaEEEE);

  RooDataHist hdatacEEEE("hdatacEEEE","hdatacEEEE",RooArgSet(mass), h_mZc_EEEE);
  RooAddPdf   modelcEEEE("modelcEEEE", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rcEEEE = modelcEEEE.fitTo(hdatacEEEE, Range(70., 110.), Save(),Verbose(0));
  RooPlot* plotcEEEE = mass.frame(Range(70,110));
  hdatacEEEE.plotOn(plotcEEEE);
  modelcEEEE.plotOn(plotcEEEE, LineColor(kGreen+2));
  modelcEEEE.paramOn(plotcEEEE, Layout(0.16,0.46,0.63));
  modelcEEEE.plotOn(plotcEEEE, LineColor(kGreen+2));
  modelcEEEE.plotOn(plotcEEEE, Range(mlcEEEE,mhcEEEE) , LineColor(kGreen+2), DrawOption("F"), FillColor(kGreen+2), FillStyle(3004) , VLines(), MoveToBack() );
  plotcEEEE->SetTitleOffset(1.7,"Y");
  plotcEEEE->Draw("same");
  TPaveText *pcEEEE = (TPaveText*)cEEEE->FindObject("modelcEEEE_paramBox");
  pcEEEE->SetTextSize(0.025);
  pcEEEE->SetLineColor(kGreen+2);
  //evalEffectiveSigma(&modelcEEEE, &mass, mlcEEEE, mhcEEEE, effsigmacEEEE);



  cout << "EBEB uncorrected --> EffectiveSigma = " << effsigmaEBEB  << " GeV" << "  M_min = " << mlEBEB << "  M_max = " << mhEBEB<< endl;
  cout << "EBEB corrected   --> EffectiveSigma = " << effsigmacEBEB << " GeV" << "  M_min = " << mlcEBEB << "  M_max = " << mhcEBEB<< endl;
  if (effsigmaEBEB > effsigmacEBEB)  cout << "EBEB smearing    --> " << sqrt( effsigmaEBEB*effsigmaEBEB - effsigmacEBEB*effsigmacEBEB ) << endl;
  else cout <<  "EBEB smearing    -->  sigma_corr > sigma_corr !!! " << endl;

  cout << "EBEE uncorrected --> EffectiveSigma = " << effsigmaEBEE  << " GeV" << "  M_min = " << mlEBEE << "  M_max = " << mhEBEE<< endl;
  cout << "EBEE corrected   --> EffectiveSigma = " << effsigmacEBEE << " GeV" << "  M_min = " << mlcEBEE << "  M_max = " << mhcEBEE<< endl;
  if (effsigmaEBEE > effsigmacEBEE)  cout << "EBEE smearing    --> " << sqrt( effsigmaEBEE*effsigmaEBEE - effsigmacEBEE*effsigmacEBEE ) << endl;
  else cout <<  "EBEE smearing    -->  sigma_corr > sigma_corr !!! " << endl;
  cout << "EEEE uncorrected --> EffectiveSigma = " << effsigmaEEEE  << " GeV" << "  M_min = " << mlEEEE << "  M_max = " << mhEEEE<< endl;
  cout << "EEEE corrected   --> EffectiveSigma = " << effsigmacEEEE << " GeV" << "  M_min = " << mlcEEEE << "  M_max = " << mhcEEEE<< endl;
  if (effsigmaEEEE > effsigmacEEEE)  cout << "EEEE smearing    --> " << sqrt( effsigmaEEEE*effsigmaEEEE - effsigmacEEEE*effsigmacEEEE ) << endl;
  else cout <<  "EEEE smearing    -->  sigma_corr > sigma_corr !!! " << endl;



  // wring effectivesigma on canvas 
  char title1[100], title2[100];
  TLatex *latex1;
  TLatex *latex2;

  cEBEB->cd();
  sprintf (title1, "#sigma_{eff} = %.2f GeV",effsigmaEBEB);
  sprintf (title2, "#sigma_{eff} = %.2f GeV",effsigmacEBEB);
  latex1 = new TLatex(0.65,0.8,title1);
  latex1->SetNDC();
  latex1->SetTextFont(42);
  latex1->SetTextSize(0.03);
  latex1->SetTextColor(kRed);
  latex1->Draw("same");
  latex2 = new TLatex(0.65,0.75,title2);
  latex2->SetNDC();
  latex2->SetTextFont(42);
  latex2->SetTextSize(0.03);
  latex2->SetTextColor(kGreen+2);
  latex2->Draw("same");

  cEBEE->cd();
  sprintf (title1, "#sigma_{eff} = %.2f GeV",effsigmaEBEE);
  sprintf (title2, "#sigma_{eff} = %.2f GeV",effsigmacEBEE);
  latex1 = new TLatex(0.65,0.8,title1);
  latex1->SetNDC();
  latex1->SetTextFont(42);
  latex1->SetTextSize(0.03);
  latex1->SetTextColor(kRed);
  latex1->Draw("same");
  latex2 = new TLatex(0.65,0.75,title2);
  latex2->SetNDC();
  latex2->SetTextFont(42);
  latex2->SetTextSize(0.03);
  latex2->SetTextColor(kGreen+2);
  latex2->Draw("same");

  cEEEE->cd();
  sprintf (title1, "#sigma_{eff} = %.2f GeV",effsigmaEEEE);
  sprintf (title2, "#sigma_{eff} = %.2f GeV",effsigmacEEEE);
  latex1 = new TLatex(0.65,0.8,title1);
  latex1->SetNDC();
  latex1->SetTextFont(42);
  latex1->SetTextSize(0.03);
  latex1->SetTextColor(kRed);
  latex1->Draw("same");
  latex2 = new TLatex(0.65,0.75,title2);
  latex2->SetNDC();
  latex2->SetTextFont(42);
  latex2->SetTextSize(0.03);
  latex2->SetTextColor(kGreen+2);
  latex2->Draw("same");


  // save in a file
  TFile *fout = new TFile(outfilename,"recreate");
  h_mZ_EBEB->Write();
  h_mZc_EBEB->Write();
  h_mZ_EBEE->Write();
  h_mZc_EBEE->Write();
  h_mZ_EEEE->Write();
  h_mZc_EEEE->Write();
  mZ_vs_localEta_EBEB->Write();
  mZc_vs_localEta_EBEB->Write();
  mZ_vs_localPhi_EBEB->Write();
  mZc_vs_localPhi_EBEB->Write();
  hspread->Write();
  hspreadc->Write();

  plotEBEB->Write("plotEBEB");
  plotcEBEB->Write("plotcEBEB");
  plotEBEE->Write("plotEBEE");
  plotcEBEE->Write("plotcEBEE");
  plotEEEE->Write("plotEEEE");
  plotcEEEE->Write("plotcEEEE");

  cout << "Closing file..." << endl;
  fout->Close();



}

