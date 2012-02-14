using namespace RooFit;

//****Effective sigma from histogram ****************************************************************
void evalEffectiveSigmaFromHisto( TH1F *h , double &ml, double &mh, double &sigmaeff){

  cout << "computing effective sigma..."<< endl;

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
    //    if (i%100==0) cout <<  i << endl;
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
  //cout << "Mmin = " << mlmin << "  Mmax = " << mhmin << "  effective sigma = " << sigmaeff << endl;
  
  // return (sigmaeff);
  return;
}

void DrawZpeak()
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
  

  TFile *_file0 = TFile::Open("dataZpeak_allR9_nBC1_singleBCcorr.root");
  gDirectory.ls();
  
  //COMPUTE EFFECTIVE SIGMA FROM HISTOS
  double mlEBEB = 0, mhEBEB = 0, effsigmaEBEB = 0;
  double mlcEBEB = 0, mhcEBEB = 0, effsigmacEBEB = 0;
  double mlEBEE = 0, mhEBEE = 0, effsigmaEBEE = 0;
  double mlcEBEE = 0, mhcEBEE = 0, effsigmacEBEE = 0;
  double mlEEEE = 0, mhEEEE = 0, effsigmaEEEE = 0;
  double mlcEEEE = 0, mhcEEEE = 0, effsigmacEEEE = 0;
  //evalEffectiveSigmaFromHisto(h_mZ_EBEB, mlEBEB, mhEBEB, effsigmaEBEB);
  //evalEffectiveSigmaFromHisto(h_mZc_EBEB, mlcEBEB, mhcEBEB, effsigmacEBEB);
  // evalEffectiveSigmaFromHisto(h_mZ_EBEE, mlEBEE, mhEBEE, effsigmaEBEE);
  // evalEffectiveSigmaFromHisto(h_mZc_EBEE, mlcEBEE, mhcEBEE, effsigmacEBEE);
  // evalEffectiveSigmaFromHisto(h_mZ_EEEE, mlEEEE, mhEEEE, effsigmaEEEE);
  // evalEffectiveSigmaFromHisto(h_mZc_EEEE, mlcEEEE, mhcEEEE, effsigmacEEEE);

  int nre = 10;
  h_mZ_EBEB  -> Rebin(nre);
  h_mZc_EBEB -> Rebin(nre);
 
  //------------------------
  //---------------- fitting
  float minM = 70.0;
  float maxM = 110.0;

  RooRealVar  mass("mass","M(e^{+}e^{-})", minM, maxM,"GeV/c^{2}");
  mass.setBins(10000) ;


  // Parameters for Crystal Ball Lineshape 
  RooRealVar  dm("#Delta m", "offset", 0.0, -5.0, 5.0,"GeV/c^{2}"); 
  RooRealVar  sigma("#sigma_{CB}","sigmaCB", 1.5,0.5,7.5,"GeV/c^{2}"); 
  RooRealVar  alpha("#alpha","alpha", 1.65,0.6,4.0); 
  RooRealVar  n("n","n", 1.47, 0.5, 50.0); 
  //alpha.setConstant();// from MC
  //n.setConstant();// from MC

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
  TCanvas *cEBEB = new TCanvas("cEBEB","cEBEB",450,600); 
  
  cEBEB->SetLeftMargin(0.15);  
  RooDataHist hdataEBEB("hdataEBEB","hdataEBEB",RooArgSet(mass), h_mZ_EBEB);
  RooAddPdf   modelEBEB("modelEBEB", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rEBEB = modelEBEB.fitTo(hdataEBEB, Range(minM, maxM), Save(),Verbose(0) );
  RooPlot* plotEBEB = mass.frame(Range(minM,maxM));
  hdataEBEB.plotOn(plotEBEB);
  modelEBEB.paramOn(plotEBEB, Layout(0.16,0.46,0.89));
  modelEBEB.plotOn(plotEBEB, LineColor(kRed+1),LineWidth(2) );
  //modelEBEB.plotOn(plotEBEB, Range(mlEBEB,mhEBEB) , LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3005) , VLines(), MoveToBack() );
  plotEBEB->SetTitleOffset(1.7,"Y");
  plotEBEB->Draw();
  TPaveText *pEBEB = (TPaveText*)cEBEB->FindObject("modelEBEB_paramBox");
  pEBEB->SetTextSize(0.025);
  pEBEB->SetLineColor(kRed+1);

  TCanvas *cEBEBc = new TCanvas("cEBEBc","cEBEBc",450,600);
  cEBEBc->SetLeftMargin(0.15);
  RooDataHist hdatacEBEB("hdatacEBEB","hdatacEBEB",RooArgSet(mass), h_mZc_EBEB);
  RooAddPdf   modelcEBEB("modelcEBEB", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rcEBEB = modelcEBEB.fitTo(hdatacEBEB, Range(minM, maxM), Save(),Verbose(0));
  RooPlot* plotcEBEB = mass.frame(Range(minM,maxM));
  hdatacEBEB.plotOn(plotcEBEB);
  modelcEBEB.paramOn(plotcEBEB, Layout(0.16,0.46,0.63));
  modelcEBEB.plotOn(plotcEBEB, LineColor(kGreen+2),LineWidth(2));
  //modelcEBEB.plotOn(plotcEBEB, Range(mlcEBEB,mhcEBEB) , LineColor(kGreen+2), DrawOption("F"), FillColor(kGreen+2), FillStyle(3004) , VLines(), MoveToBack() );
  plotcEBEB->SetTitleOffset(1.7,"Y");
  plotcEBEB->Draw("");
  TPaveText *pcEBEB = (TPaveText*)cEBEBc->FindObject("modelcEBEB_paramBox");
  pcEBEB ->SetFillStyle(1);
  pcEBEB ->SetFillColor(0);
  pcEBEB->SetTextSize(0.025);
  pcEBEB->SetLineColor(kGreen+2);
  
  /*
  // --- EB-EE
  TCanvas *cEBEE = new TCanvas("cEBEE","cEBEE",700,700);
  cEBEE->SetLeftMargin(0.15);

  RooDataHist hdataEBEE("hdataEBEE","hdataEBEE",RooArgSet(mass), h_mZ_EBEE);
  RooAddPdf   modelEBEE("modelEBEE", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *rEBEE = modelEBEE.fitTo(hdataEBEE, Range(minM, maxM), Save() ,Verbose(0));
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
  RooFitResult *rcEBEE = modelcEBEE.fitTo(hdatacEBEE, Range(minM, maxM), Save(),Verbose(0));
  RooPlot* plotcEBEE = mass.frame(Range(minM, maxM));
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
  RooFitResult *rEEEE = modelEEEE.fitTo(hdataEEEE, Range(minM, maxM), Save() ,Verbose(0));
  RooPlot* plotEEEE = mass.frame(Range(minM,maxM));
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
  RooFitResult *rcEEEE = modelcEEEE.fitTo(hdatacEEEE, Range(minM, maxM), Save(),Verbose(0));
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
  */


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

  
  /*
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
  /*
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
  */

}
