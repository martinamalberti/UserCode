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

void FitZpeak(TH1F *h , float minM, float maxM, float &deltaM, float &sigmaCB, float &deltaMerr, float &sigmaCBerr, RooPlot *plot)
{
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
  
  // Fit
  RooDataHist hdata("hdata","hdata",RooArgSet(mass), h);
  RooAddPdf   model("model", "Signal + Background", bwcb, bkg, frac);  
  RooFitResult *r = model.fitTo(hdata, Range(minM, maxM), Save(),Verbose(0) );
  cout <<  "DeltaM = " << dm.getVal() << endl;
  cout <<  "SigmaCB= " << sigma.getVal()  << endl;
  // cout <<  "DeltaM = " << cb.getVal("dm")  << endl;
  // cout <<  "SigmaCB= " << cb.getVal("sigma")  << endl;
  deltaM  = dm.getVal() ;
  sigmaCB = sigma.getVal() ;
  deltaMerr  = dm.getError() ;
  sigmaCBerr = sigma.getError() ;
 
  //  RooPlot* plot = mass.frame(Range(minM,maxM));
  plot = mass.frame(Range(minM,maxM));
  hdata.plotOn(plot);
  model.plotOn(plot);
  model.plotOn(plot, LineColor(kRed));
  model.paramOn(plot, Layout(0.16,0.46,0.89));
  model.plotOn(plot, LineColor(kRed));
  //model.plotOn(plot, Range(ml,mh) , LineColor(kRed), DrawOption("F"), FillColor(kRed), FillStyle(3005) , VLines(), MoveToBack() );
  plot->SetTitleOffset(1.7,"Y");
  // if (drawSameOption==0)  plot->Draw();
  //if (drawSameOption==1)  plot->Draw("same");
  //TPaveText *p = (TPaveText*)c->FindObject("model_paramBox");
  //p->SetTextSize(0.025);
  //p->SetLineColor(kRed);
  
}

void DrawZpeakVsNvtx()
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
  

  TFile *f = TFile::Open("outputFiles/dataZpeakVsNvtx_allR9_PU.root");
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
  
  TH1F *h_mZ_EBEB[50];
  TH1F *h_mZc_EBEB[50];
  TH1F *h_mZr_EBEB[50];
  int counter = 0;

  for (int i=0; i<50; i++){
    char hname[100];
    sprintf(hname, "h_mZ_EBEB_%d",i);
    h_mZ_EBEB[i] = (TH1F*)f->Get(hname);
    if (h_mZ_EBEB[i]){
      h_mZ_EBEB[i] -> Rebin(nre);
      counter++; 
    }

    sprintf(hname, "h_mZc_EBEB_%d",i);
    h_mZc_EBEB[i] = (TH1F*)f->Get(hname);
    if (h_mZc_EBEB[i]) h_mZc_EBEB[i] -> Rebin(nre);

    sprintf(hname, "h_mZr_EBEB_%d",i);
    h_mZr_EBEB[i] = (TH1F*)f->Get(hname);
    if (h_mZr_EBEB[i]) h_mZr_EBEB[i] -> Rebin(nre);
  }

  cout << counter << endl;
  
  //------------------------
  //---------------- fitting
  float minM = 70.0;
  float maxM = 110.0;
  float MZ = 91.188;

  // --- EB-EB
  TGraphErrors *gEBsigma   = new TGraphErrors();
  TGraphErrors *gEBsigmac  = new TGraphErrors();
  TGraphErrors *gEBsigmar  = new TGraphErrors();
  TGraphErrors *gEBdeltam  = new TGraphErrors();
  TGraphErrors *gEBdeltamc = new TGraphErrors();
  TGraphErrors *gEBdeltamr = new TGraphErrors();

  float sigmaCB , deltaM, sigmaCBerr, deltaMerr;
  TCanvas *cEBEB[100];
  TCanvas *cEBEBc[100];
  TCanvas *cEBEBr[100];
  RooPlot *plotEBEB[100];
  RooPlot *plotEBEBc[100];
  RooPlot *plotEBEBr[100];

  char cname[100];
  for (int i = 0; i < 22; i++){
 
    // sprintf(cname, "cEBEB_%d",i);
    // cEBEB[i] = new TCanvas(cname,cname,700,700);
    // cEBEB[i]->SetLeftMargin(0.15);
    //if ( h_mZ_EBEB[i] )    FitZpeak(h_mZ_EBEB[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, cEBEB[i], 0);
    if ( h_mZ_EBEB[i] )    FitZpeak(h_mZ_EBEB[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, plotEBEB[i]);

    // sigmaCB/=(MZ+deltaM);
    // sigmaCBerr/=(MZ+deltaM);
    gsigma  -> SetPoint(i, i, sigmaCB);
    gsigma  -> SetPointError(i, 0, sigmaCBerr);
    gdeltam -> SetPoint(i, i, deltaM);
    gdeltam -> SetPointError(i, 0, deltaMerr);
    
    //sprintf(cname, "cEBEBc_%d",i);
    //cEBEBc[i] = new TCanvas(cname,cname,700,700);
    //cEBEBc[i]->SetLeftMargin(0.15);
    //if ( h_mZc_EBEB[i] )  FitZpeak(h_mZc_EBEB[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, cEBEBc[i], 0);
    if ( h_mZc_EBEB[i] )  FitZpeak(h_mZc_EBEB[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, plotEBEBc[i]);
    
    // sigmaCB/=(MZ+deltaM);
    // sigmaCBerr/=(MZ+deltaM);
    gsigmac  -> SetPoint(i, i, sigmaCB);
    gsigmac  -> SetPointError(i, 0, sigmaCBerr);
    gdeltamc -> SetPoint(i, i, deltaM);
    gdeltamc -> SetPointError(i, 0, deltaMerr);

    //sprintf(cname, "cEBEBr_%d",i);
    //cEBEBr[i] = new TCanvas(cname,cname,700,700);
    //cEBEBr[i]->SetLeftMargin(0.15);
    //if ( h_mZr_EBEB[i] )  FitZpeak(h_mZr_EBEB[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, cEBEBr[i], 0);
    if ( h_mZr_EBEB[i] )  FitZpeak(h_mZr_EBEB[i], minM, maxM, deltaM, sigmaCB, deltaMerr, sigmaCBerr, plotEBEBr[i]);
    
    // sigmaCB/=(MZ+deltaM);
    // sigmaCBerr/=(MZ+deltaM);
    gsigmar  -> SetPoint(i, i, sigmaCB);
    gsigmar  -> SetPointError(i, 0, sigmaCBerr);
    gdeltamr -> SetPoint(i, i, deltaM);
    gdeltamr -> SetPointError(i, 0, deltaMerr);


    /*
    // writing effectivesigma on canvas 
    char title1[100], title2[100];
    TLatex *latex1;
    TLatex *latex2;
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
    */
  }

  gsigma->SetMarkerStyle(20);
  gsigma->SetMarkerColor(2);
  gdeltam->SetMarkerStyle(20);
  gdeltam->SetMarkerColor(2);

  gsigmac->SetMarkerStyle(20);
  gsigmac->SetMarkerColor(4);
  gdeltamc->SetMarkerStyle(20);
  gdeltamc->SetMarkerColor(4);

  gsigmar->SetMarkerStyle(20);
  gsigmar->SetMarkerColor(3);
  gdeltamr->SetMarkerStyle(20);
  gdeltamr->SetMarkerColor(3);


  TCanvas *c = new TCanvas("c","c",500,500);
  gsigma-> GetYaxis()->SetRangeUser(0, 3);
  gsigma-> GetYaxis()->SetTitle("#sigma_{CB} (GeV)");
  gsigma-> GetXaxis()->SetTitle("number of vertices");
  gsigma->Draw("ap");
  gsigmac->Draw("psame");
  gsigmar->Draw("psame");

  TCanvas *cc = new TCanvas("cc","cc",500,500);
  gdeltam -> GetYaxis()->SetRangeUser(-2, 2);
  gdeltam -> GetYaxis()->SetTitle("#Delta M (GeV)");
  gdeltam -> GetXaxis()->SetTitle("number of vertices");
  gdeltam->Draw("ap");
  gdeltamc->Draw("psame");
  gdeltamr->Draw("psame");


  TFile *fout = new TFile("out.root","recreate");
  gsigma->Write();
  gsigmac->Write();
  gsigmar->Write();
  gdeltam ->Write();
  gdeltamc ->Write();
  gdeltamr ->Write();

  for (int i = 0; i < 4 ; i++){
    plotEBEB[i]->Write();
    plotEBEBc[i]->Write();
    plotEBEBr[i]->Write();

  }

  cout << "Closing file..." << endl;
  fout->Close();
  

}
