//
// Macro to produce ECAL Pi0 plots
//

// int Wait() {
//      cout << " Continue [<RET>|q]?  ";
//      char x;
//      x = getchar();
//      if ((x == 'q') || (x == 'Q')) return 1;
//      return 0;
// }


void DrawValidationPlotsPi0(Char_t* infile1 = 0, 
		     Char_t* infile2 = 0, 
		     Char_t* fileType = "png", 
		     Char_t* dirName = ".")
{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 
  gStyle->SetOptStat(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  if (!infile1 || !infile2) {
    cout << " No input file specified !" << endl;
    return;
  }
  
  cout << "Producing validation plots for: " << infile1 << " and " << infile2 << endl;

  TFile* f[2];
  f[0] = new TFile(infile1);
  f[1] = new TFile(infile2);

  char name[500];

  int nHists = 100;
  
  TCanvas* c[100];
  
  char cname[100];

  for (int i=0; i<nHists+1; i++) {
    int x = (i%3)*600;     //int x = (i%3)*600;
    int y = (i/3)*100;     //int y = (i/3)*200;
    sprintf(cname,"c%i",i);
    c[i] =  new TCanvas(cname,cname,x,y,600,400);
  }
  

  float minvMax = 0.350;
  float minvMin = 0.06;

  TH1D *h_numberOfEvents[2];

  TH1D *h_Pi0_EB_mass[2];
  TH1D *h_Pi0_EE_mass[2];

  for (int i = 0; i < 2 ; i++) {
    
    // number of events (needed for normalization)
    h_numberOfEvents[i] = (TH1D*)f[i]->Get("ecalvalidation/h_numberOfEvents") ; 

    h_Pi0_EB_mass[i]    = (TH1D*)f[i]->Get("ecalvalidation/Pi0/h_Pi0_EB_mass") ; 
    h_Pi0_EE_mass[i]    = (TH1D*)f[i]->Get("ecalvalidation/Pi0/h_Pi0_EE_mass") ; 

  }


  float nEvents1 = h_numberOfEvents[0]->GetBinContent(1);
  float nEvents2 = h_numberOfEvents[1]->GetBinContent(1);
  float s        = nEvents1/nEvents2;
  float maxY;

  // pi0 mass
  c[0]->cd();
  h_Pi0_EB_mass[0]->SetTitle("Pi0 peak (EB)");
  h_Pi0_EB_mass[0]->GetXaxis()->SetTitle("mass (GeV)");
  h_Pi0_EB_mass[0]->GetXaxis()->SetRangeUser(minvMin,minvMax);
  maxY = max(h_Pi0_EB_mass[0]->GetMaximum(), h_Pi0_EB_mass[1]->GetMaximum()*s);
  h_Pi0_EB_mass[0]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_Pi0_EB_mass[0]->SetLineColor(4);
  h_Pi0_EB_mass[0]->Draw(); 
  h_Pi0_EB_mass[1]->Scale(s);
  h_Pi0_EB_mass[1]->SetLineColor(2);
  h_Pi0_EB_mass[1]->Draw("same"); 
  sprintf(name,"%s/pi0_EB_mass.%s",dirName,fileType); 
  c[0]->Print(name);

  
  c[1]->cd();
  h_Pi0_EE_mass[0]->SetTitle("Pi0 peak (EE)");
  h_Pi0_EE_mass[0]->GetXaxis()->SetTitle("mass (GeV)");
  h_Pi0_EE_mass[0]->GetXaxis()->SetRangeUser(minvMin,minvMax);
  maxY = max(h_Pi0_EE_mass[0]->GetMaximum(), h_Pi0_EE_mass[1]->GetMaximum()*s);
  h_Pi0_EE_mass[0]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_Pi0_EE_mass[0]->SetLineColor(4);
  h_Pi0_EE_mass[0]->Draw(); 
  h_Pi0_EE_mass[1]->Scale(s);
  h_Pi0_EE_mass[1]->SetLineColor(2);
  h_Pi0_EE_mass[1]->Draw("same"); 
  sprintf(name,"%s/pi0_EE_mass.%s",dirName,fileType); 
  c[1]->Print(name);

}
