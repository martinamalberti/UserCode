//
// Macro to produce ECAL cosmic plots
//

// int Wait() {
//      cout << " Continue [<RET>|q]?  ";
//      char x;
//      x = getchar();
//      if ((x == 'q') || (x == 'Q')) return 1;
//      return 0;
// }

#include <algorithm>

void DrawValidationPlots(Char_t* infile1 = 0, 
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
  
  TCanvas* cRecHits[100];
  TCanvas* cSC[100];
  TCanvas* cES[100];
 
  char cname[100];

  for (int i=0; i<nHists+1; i++) {
    int x = (i%3)*600;     //int x = (i%3)*600;
    int y = (i/3)*100;     //int y = (i/3)*200;
    sprintf(cname,"cRecHits%i",i);
    cRecHits[i] =  new TCanvas(cname,cname,x,y,600,400);
    sprintf(cname,"cSC%i",i);
    cSC[i] =  new TCanvas(cname,cname,x,y,600,400);
    sprintf(cname,"cES%i",i);
    cES[i] =  new TCanvas(cname,cname,x,y,600,400);
  }
  

  float nMaxHits = 2000;
  float eMinHits = -1.0;
  float eMaxHits = 10;
  float nMaxBC   = 100;
  float eMaxBC   = 200;
  float nMaxSC   = 100;
  float eMaxSC   = 200;
  float nBCMax    = 10;
  float nSCMax    = 50;
  float nXtalsMax = 50;
 

  TH1D *h_numberOfEvents[2];


  // RecHits ----------------------------------------------

  // Rec Hits multiplicity
  TH1D *h_recHits_EB_size[2];
  TH1D *h_recHits_EEP_size[2];
  TH1D *h_recHits_EEM_size[2];
  TH1D *h_recHits_ES_size[2];
 
  // Rec Hits energy
  TH1D *h_recHits_EB_energy[2];
  TH1D *h_recHits_EEP_energy[2];
  TH1D *h_recHits_EEM_energy[2];
  TH1D *h_recHits_ES_energy[2];

  // Rec Hits energy MAX
  TH1D *h_recHits_EB_energyMax[2];
  TH1D *h_recHits_EEP_energyMax[2];
  TH1D *h_recHits_EEM_energyMax[2];
  TH1D *h_recHits_ES_energyMax[2];

  // Rec Hits Time
  TH1D *h_recHits_EB_time[2];
  TH1D *h_recHits_EEP_time[2];
  TH1D *h_recHits_EEM_time[2];
  TH1D *h_recHits_ES_time[2];

  // Rec Hits PHI
  TH1D *h_recHits_EB_phi[2];
  TH1D *h_recHits_EE_phi[2];

  // Rec Hits ETA
  TH1D *h_recHits_eta[2];

  
  // Clusters ----------------------------------------------
  
  // number of SC per event
  TH1D *h_superClusters_EB_size[2];
  TH1D *h_superClusters_EEP_size[2];
  TH1D *h_superClusters_EEM_size[2];

  // number of xtals in SC
  TH1D *h_superClusters_EB_nXtals[2];
  TH1D *h_superClusters_EEP_nXtals[2];
  TH1D *h_superClusters_EEM_nXtals[2];
 
  // number of BC in SC
  TH1D *h_superClusters_EB_nBC[2];
  TH1D *h_superClusters_EEP_nBC[2];
  TH1D *h_superClusters_EEM_nBC[2];

  // SC energy
  TH1D *h_superClusters_EB_energy[2];
  TH1D *h_superClusters_EEP_energy[2];
  TH1D *h_superClusters_EEM_energy[2];

  // SC eta
  TH1D *h_superClusters_eta[2];

  // SC phi
  TH1D *h_superClusters_EB_phi[2];
  TH1D *h_superClusters_EE_phi[2];
  
  // ES
  TH1D *h_esClusters_energy_plane1[2];
  TH1D *h_esClusters_energy_plane2[2];
  TH1D *h_esClusters_energy_ratio[2];

  for (int i = 0; i < 2 ; i++) {
    
    // number of events (needed for normalization)
    h_numberOfEvents[i]           = (TH1D*)f[i]->Get("ecalvalidation/h_numberOfEvents") ; 

    // rec hits
    h_recHits_EB_size[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EB_size") ;
    h_recHits_EEP_size[i]    = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EEP_size") ;
    h_recHits_EEM_size[i]    = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EEM_size") ;
    h_recHits_ES_size[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_ES_size") ;
 
    h_recHits_EB_energy[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EB_energy") ;
    h_recHits_EEP_energy[i]    = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EEP_energy") ;
    h_recHits_EEM_energy[i]    = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EEM_energy") ;
    h_recHits_ES_energy[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_ES_energy") ;
  
    h_recHits_EB_energyMax[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EB_energyMax") ;
    h_recHits_EEP_energyMax[i]    = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EEP_energyMax") ;
    h_recHits_EEM_energyMax[i]    = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EEM_energyMax") ;
    h_recHits_ES_energyMax[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_ES_energyMax") ;
  
    h_recHits_EB_time[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EB_time") ;
    h_recHits_EEP_time[i]    = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EEP_time") ;
    h_recHits_EEM_time[i]    = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EEM_time") ;
    h_recHits_ES_time[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_ES_time") ;

    h_recHits_eta[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_eta") ;
    h_recHits_EB_phi[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EB_phi") ;
    h_recHits_EE_phi[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_recHits_EE_phi") ;


    // Super Clusters ----------------------------------------------
    h_superClusters_EB_size[i]    = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EB_size") ;
    h_superClusters_EEP_size[i]   = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EEP_size") ;
    h_superClusters_EEM_size[i]   = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EEM_size") ;
    
    h_superClusters_EB_nXtals[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EB_nXtals") ;
    h_superClusters_EEP_nXtals[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EEP_nXtals") ;
    h_superClusters_EEM_nXtals[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EEM_nXtals") ;
   
    h_superClusters_EB_nBC[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EB_nBC") ;
    h_superClusters_EEP_nBC[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EEP_nBC") ;
    h_superClusters_EEM_nBC[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EEM_nBC") ;
   
    h_superClusters_EB_energy[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EB_energy") ;
    h_superClusters_EEP_energy[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EEP_energy") ;
    h_superClusters_EEM_energy[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EEM_energy") ;
 
    h_superClusters_eta[i]       = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_eta") ;
    h_superClusters_EB_phi[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EB_phi") ;
    h_superClusters_EE_phi[i]     = (TH1D*)f[i]->Get("ecalvalidation/h_superClusters_EE_phi") ;

   
    h_esClusters_energy_plane1[i] = (TH1D*)f[i]->Get("ecalvalidation/h_esClusters_energy_plane1") ;
    h_esClusters_energy_plane2[i] = (TH1D*)f[i]->Get("ecalvalidation/h_esClusters_energy_plane2") ;
    h_esClusters_energy_ratio[i]  = (TH1D*)f[i]->Get("ecalvalidation/h_esClusters_energy_ratio") ;

  }


  float nEvents1 = h_numberOfEvents[0]->GetBinContent(1);
  float nEvents2 = h_numberOfEvents[1]->GetBinContent(1);
  float s        = nEvents1/nEvents2;
  float maxY;

  int i1 = 0;
  int i2 = 1;
  int icolor1 = 4;
  int icolor2 = 2;

 
  if ( nEvents2  > nEvents1) {
    i1 = 1;
    i2 = 0;
    s  = nEvents2/nEvents1;
    icolor1 = 2;
    icolor2 = 4;
   }


  // Rec Hits number
  cRecHits[0]->cd();
  cRecHits[0]->SetLogy();
  h_recHits_EB_size[i1]->SetTitle("Number of RecHits (EB)");
  h_recHits_EB_size[i1]->GetXaxis()->SetTitle("Number of RecHits");
  h_recHits_EB_size[i1]->GetXaxis()->SetRangeUser(0,nMaxHits);
  maxY = max(h_recHits_EB_size[i1]->GetMaximum(), h_recHits_EB_size[i2]->GetMaximum()*s);
  h_recHits_EB_size[i1]->GetYaxis()->SetRangeUser(0.1,maxY*10);
  h_recHits_EB_size[i1]->SetLineColor(icolor1);
  h_recHits_EB_size[i1]->Draw(); 
  h_recHits_EB_size[i2]->Scale(s);
  h_recHits_EB_size[i2]->SetLineColor(icolor2);
  h_recHits_EB_size[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EB_size.%s",dirName,fileType); 
  cRecHits[0]->Print(name);
  
  cRecHits[1]->cd();
  cRecHits[1]->SetLogy();
  h_recHits_EEP_size[i1]->SetTitle("Number of RecHits (EE+)");
  h_recHits_EEP_size[i1]->GetXaxis()->SetTitle("Number of RecHits");
  h_recHits_EEP_size[i1]->GetXaxis()->SetRangeUser(0,nMaxHits);
  maxY = max(h_recHits_EEP_size[i1]->GetMaximum(), h_recHits_EEP_size[i2]->GetMaximum()*s);
  h_recHits_EEP_size[i1]->GetYaxis()->SetRangeUser(0.1,maxY*10);
  h_recHits_EEP_size[i1]->SetLineColor(icolor1);
  h_recHits_EEP_size[i1]->Draw(); 
  h_recHits_EEP_size[i2]->Scale(s);
  h_recHits_EEP_size[i2]->SetLineColor(icolor2);
  h_recHits_EEP_size[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EEP_size.%s",dirName,fileType); 
  cRecHits[1]->Print(name);

  cRecHits[2]->cd();
  cRecHits[2]->SetLogy();
  h_recHits_EEM_size[i1]->SetTitle("Number of RecHits (EE-)");
  h_recHits_EEM_size[i1]->GetXaxis()->SetTitle("Number of RecHits");
  h_recHits_EEM_size[i1]->GetXaxis()->SetRangeUser(0,nMaxHits);
  maxY = max(h_recHits_EEM_size[i1]->GetMaximum(), h_recHits_EEM_size[i2]->GetMaximum()*s);
  h_recHits_EEM_size[i1]->GetYaxis()->SetRangeUser(0.1,maxY*10);
  h_recHits_EEM_size[i1]->SetLineColor(icolor1);
  h_recHits_EEM_size[i1]->Draw(); 
  h_recHits_EEM_size[i2]->Scale(s);
  h_recHits_EEM_size[i2]->SetLineColor(icolor2);
  h_recHits_EEM_size[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EEM_size.%s",dirName,fileType); 
  cRecHits[2]->Print(name);

  cRecHits[3]->cd();
  cRecHits[3]->SetLogy();
  h_recHits_ES_size[i1]->SetTitle("Number of RecHits (ES)");
  h_recHits_ES_size[i1]->GetXaxis()->SetTitle("Number of RecHits");
  h_recHits_ES_size[i1]->GetXaxis()->SetRangeUser(0,nMaxHits);
  maxY = max(h_recHits_ES_size[i1]->GetMaximum(), h_recHits_ES_size[i2]->GetMaximum()*s);
  h_recHits_ES_size[i1]->GetYaxis()->SetRangeUser(0.1,maxY*10);
  h_recHits_ES_size[i1]->SetLineColor(icolor1);
  h_recHits_ES_size[i1]->Draw(); 
  h_recHits_ES_size[i2]->Scale(s);
  h_recHits_ES_size[i2]->SetLineColor(icolor2);
  h_recHits_ES_size[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_ES_size.%s",dirName,fileType); 
  cRecHits[3]->Print(name);


  // Rec Hits energy
  cRecHits[4]->cd();
  cRecHits[4]->SetLogy();
  h_recHits_EB_energy[i1]->SetTitle("Rec Hits Energy (EB)");
  h_recHits_EB_energy[i1]->GetXaxis()->SetTitle("energy (GeV)");
  h_recHits_EB_energy[i1]->GetXaxis()->SetRangeUser(eMinHits,eMaxHits);
  h_recHits_EB_energy[i1]->SetLineColor(icolor1);
  h_recHits_EB_energy[i1]->Draw(); 
  h_recHits_EB_energy[i2]->Scale(s);
  h_recHits_EB_energy[i2]->SetLineColor(icolor2);
  h_recHits_EB_energy[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EB_energy.%s",dirName,fileType); 
  cRecHits[4]->Print(name);
  
  cRecHits[5]->cd();
  cRecHits[5]->SetLogy();
  h_recHits_EEP_energy[i1]->SetTitle("Rec Hits Energy (EE+)");
  h_recHits_EEP_energy[i1]->GetXaxis()->SetTitle("energy (GeV)");
  h_recHits_EEP_energy[i1]->GetXaxis()->SetRangeUser(eMinHits,eMaxHits);
  h_recHits_EEP_energy[i1]->SetLineColor(icolor1);
  h_recHits_EEP_energy[i1]->Draw(); 
  h_recHits_EEP_energy[i2]->Scale(s);
  h_recHits_EEP_energy[i2]->SetLineColor(icolor2);
  h_recHits_EEP_energy[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EEP_energy.%s",dirName,fileType); 
  cRecHits[5]->Print(name);

  cRecHits[6]->cd();
  cRecHits[6]->SetLogy();
  h_recHits_EEM_energy[i1]->SetTitle("Rec Hits Energy (EE-)");
  h_recHits_EEM_energy[i1]->GetXaxis()->SetTitle("energy (GeV)");
  h_recHits_EEM_energy[i1]->GetXaxis()->SetRangeUser(eMinHits,eMaxHits);
  h_recHits_EEM_energy[i1]->SetLineColor(icolor1);
  h_recHits_EEM_energy[i1]->Draw(); 
  h_recHits_EEM_energy[i2]->Scale(s);
  h_recHits_EEM_energy[i2]->SetLineColor(icolor2);
  h_recHits_EEM_energy[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EEM_energy.%s",dirName,fileType); 
  cRecHits[6]->Print(name);

  cRecHits[7]->cd();
  cRecHits[7]->SetLogy();
  h_recHits_ES_energy[i1]->SetTitle("Rec Hits Energy(ES)");
  h_recHits_ES_energy[i1]->GetXaxis()->SetTitle("energy (GeV)");
  //h_recHits_ES_energy[i1]->GetXaxis()->SetRangeUser(eMinHits,eMaxHits);
  h_recHits_ES_energy[i1]->SetLineColor(icolor1);
  h_recHits_ES_energy[i1]->Draw(); 
  h_recHits_ES_energy[i2]->Scale(s);
  h_recHits_ES_energy[i2]->SetLineColor(icolor2);
  h_recHits_ES_energy[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_ES_energy.%s",dirName,fileType); 
  cRecHits[7]->Print(name);
 

  // Rec Hits energy MAX
  cRecHits[8]->cd();
  cRecHits[8]->SetLogy();
  h_recHits_EB_energyMax[i1]->SetTitle("Rec Hits Max Energy (EB)");
  h_recHits_EB_energyMax[i1]->GetXaxis()->SetTitle("energy (GeV)");
  h_recHits_EB_energyMax[i1]->GetXaxis()->SetRangeUser(0,150);
  h_recHits_EB_energyMax[i1]->SetLineColor(icolor1);
  h_recHits_EB_energyMax[i1]->Draw(); 
  h_recHits_EB_energyMax[i2]->Scale(s);
  h_recHits_EB_energyMax[i2]->SetLineColor(icolor2);
  h_recHits_EB_energyMax[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EB_energyMax.%s",dirName,fileType); 
  cRecHits[8]->Print(name);
  
  cRecHits[9]->cd();
  cRecHits[9]->SetLogy();
  h_recHits_EEP_energyMax[i1]->SetTitle("Rec Hits Max Energy (EE+)");
  h_recHits_EEP_energyMax[i1]->GetXaxis()->SetTitle("energy (GeV)");
  h_recHits_EEP_energyMax[i1]->GetXaxis()->SetRangeUser(0,150);
  h_recHits_EEP_energyMax[i1]->SetLineColor(icolor1);
  h_recHits_EEP_energyMax[i1]->Draw(); 
  h_recHits_EEP_energyMax[i2]->Scale(s);
  h_recHits_EEP_energyMax[i2]->SetLineColor(icolor2);
  h_recHits_EEP_energyMax[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EEP_energyMax.%s",dirName,fileType); 
  cRecHits[9]->Print(name);

  cRecHits[10]->cd();
  cRecHits[10]->SetLogy();
  h_recHits_EEM_energyMax[i1]->SetTitle("Rec Hits Max Energy (EE-)");
  h_recHits_EEM_energyMax[i1]->GetXaxis()->SetTitle("energy (GeV)");
  h_recHits_EEM_energyMax[i1]->GetXaxis()->SetRangeUser(0,150);
  h_recHits_EEM_energyMax[i1]->SetLineColor(icolor1);
  h_recHits_EEM_energyMax[i1]->Draw(); 
  h_recHits_EEM_energyMax[i2]->Scale(s);
  h_recHits_EEM_energyMax[i2]->SetLineColor(icolor2);
  h_recHits_EEM_energyMax[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EEM_energyMax.%s",dirName,fileType); 
  cRecHits[10]->Print(name);

  cRecHits[11]->cd();
  cRecHits[11]->SetLogy();
  h_recHits_ES_energyMax[i1]->SetTitle("Rec Hits Max Energy(ES)");
  h_recHits_ES_energyMax[i1]->GetXaxis()->SetTitle("energy (GeV)");
  //h_recHits_ES_energyMax[i1]->GetXaxis()->SetRangeUser(eMinHits,eMaxHits);
  h_recHits_ES_energyMax[i1]->SetLineColor(icolor1);
  h_recHits_ES_energyMax[i1]->Draw(); 
  h_recHits_ES_energyMax[i2]->Scale(s);
  h_recHits_ES_energyMax[i2]->SetLineColor(icolor2);
  h_recHits_ES_energyMax[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_ES_energyMax.%s",dirName,fileType); 
  cRecHits[11]->Print(name);

  // rec hits timing
  cRecHits[12]->cd();
  h_recHits_EB_time[i1]->SetTitle("Rec Hits Time (EB)");
  h_recHits_EB_time[i1]->GetXaxis()->SetTitle("time (ns)");
  maxY = max(h_recHits_EB_time[i1]->GetMaximum(), h_recHits_EB_time[i2]->GetMaximum()*s);
  h_recHits_EB_time[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_recHits_EB_time[i1]->SetLineColor(icolor1);
  h_recHits_EB_time[i1]->Draw(); 
  h_recHits_EB_time[i2]->Scale(s);
  h_recHits_EB_time[i2]->SetLineColor(icolor2);
  h_recHits_EB_time[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EB_time.%s",dirName,fileType); 
  cRecHits[12]->Print(name);
  
  cRecHits[13]->cd();
  h_recHits_EEP_time[i1]->SetTitle("Rec Hits Time (EE+)");
  h_recHits_EEP_time[i1]->GetXaxis()->SetTitle("time (ns)");
  maxY = max(h_recHits_EEP_time[i1]->GetMaximum(), h_recHits_EEP_time[i2]->GetMaximum()*s);
  h_recHits_EEP_time[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_recHits_EEP_time[i1]->SetLineColor(icolor1);
  h_recHits_EEP_time[i1]->Draw(); 
  h_recHits_EEP_time[i2]->Scale(s);
  h_recHits_EEP_time[i2]->SetLineColor(icolor2);
  h_recHits_EEP_time[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EEP_time.%s",dirName,fileType); 
  cRecHits[13]->Print(name);

  cRecHits[14]->cd();
  h_recHits_EEM_time[i1]->SetTitle("Rec Hits Time (EE-)");
  h_recHits_EEM_time[i1]->GetXaxis()->SetTitle("time (ns)");
  maxY = max(h_recHits_EEM_time[i1]->GetMaximum(), h_recHits_EEM_time[i2]->GetMaximum()*s);
  h_recHits_EEM_time[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_recHits_EEM_time[i1]->SetLineColor(icolor1);
  h_recHits_EEM_time[i1]->Draw(); 
  h_recHits_EEM_time[i2]->Scale(s);
  h_recHits_EEM_time[i2]->SetLineColor(icolor2);
  h_recHits_EEM_time[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EEM_time.%s",dirName,fileType); 
  cRecHits[14]->Print(name);

  cRecHits[15]->cd();
  h_recHits_ES_time[i1]->SetTitle("Rec Hits Time (ES)");
  h_recHits_ES_time[i1]->GetXaxis()->SetTitle("time (ns)");
  maxY = max(h_recHits_ES_time[i1]->GetMaximum(), h_recHits_ES_time[i2]->GetMaximum()*s);
  h_recHits_ES_time[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_recHits_ES_time[i1]->SetLineColor(icolor1);
  h_recHits_ES_time[i1]->Draw(); 
  h_recHits_ES_time[i2]->Scale(s);
  h_recHits_ES_time[i2]->SetLineColor(icolor2);
  h_recHits_ES_time[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_ES_time.%s",dirName,fileType); 
  cRecHits[15]->Print(name);


  // rec hits eta/phi
  cRecHits[16]->cd();
  cRecHits[16]->SetLogy();
  h_recHits_eta[i1]->SetTitle("Rec Hits Eta (EE+EB)");
  h_recHits_eta[i1]->GetXaxis()->SetTitle("#eta");
  h_recHits_eta[i1]->SetLineColor(icolor1);
  h_recHits_eta[i1]->SetLineColor(icolor1);
  h_recHits_eta[i1]->Draw(); 
  h_recHits_eta[i2]->Scale(s);
  h_recHits_eta[i2]->SetLineColor(icolor2);
  h_recHits_eta[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_eta.%s",dirName,fileType); 
  cRecHits[16]->Print(name);

  cRecHits[17]->cd();
  cRecHits[17]->SetLogy();
  h_recHits_EB_phi[i1]->SetTitle("Rec Hits Phi (EB)");
  h_recHits_EB_phi[i1]->GetXaxis()->SetTitle("#phi");
  h_recHits_EB_phi[i1]->SetLineColor(icolor1);
  h_recHits_EB_phi[i1]->Draw(); 
  h_recHits_EB_phi[i2]->Scale(s);
  h_recHits_EB_phi[i2]->SetLineColor(icolor2);
  h_recHits_EB_phi[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EB_phi.%s",dirName,fileType); 
  cRecHits[17]->Print(name);

  cRecHits[18]->cd();
  cRecHits[18]->SetLogy();
  h_recHits_EE_phi[i1]->SetTitle("Rec Hits Phi (EE)");
  h_recHits_EE_phi[i1]->GetXaxis()->SetTitle("#phi");
  h_recHits_EE_phi[i1]->SetLineColor(icolor1);
  h_recHits_EE_phi[i1]->Draw(); 
  h_recHits_EE_phi[i2]->Scale(s);
  h_recHits_EE_phi[i2]->SetLineColor(icolor2);
  h_recHits_EE_phi[i2]->Draw("same"); 
  sprintf(name,"%s/recHits_EE_phi.%s",dirName,fileType); 
  cRecHits[18]->Print(name);


  // number of SC
  cSC[0]->cd();
  cSC[0]->SetLogy();
  h_superClusters_EB_size[i1]->SetTitle("Number of SuperClusters (EB)");
  h_superClusters_EB_size[i1]->GetXaxis()->SetTitle("number of superclusters");
  h_superClusters_EB_size[i1]->GetXaxis()->SetRangeUser(0,nSCMax);
  maxY = max(h_superClusters_EB_size[i1]->GetMaximum(), h_superClusters_EB_size[i2]->GetMaximum()*s);
  h_superClusters_EB_size[i1]->GetYaxis()->SetRangeUser(0.01,maxY*10);
  h_superClusters_EB_size[i1]->SetLineColor(icolor1);
  h_superClusters_EB_size[i1]->Draw(); 
  h_superClusters_EB_size[i2]->Scale(s);
  h_superClusters_EB_size[i2]->SetLineColor(icolor2);
  h_superClusters_EB_size[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EB_size.%s",dirName,fileType); 
  cSC[0]->Print(name);

  cSC[1]->cd();
  cSC[1]->SetLogy();
  h_superClusters_EEP_size[i1]->SetTitle("Number of SuperClusters (EE+)");
  h_superClusters_EEP_size[i1]->GetXaxis()->SetTitle("number of superclusters");
  h_superClusters_EEP_size[i1]->GetXaxis()->SetRangeUser(0,nSCMax);
  maxY = max(h_superClusters_EEP_size[i1]->GetMaximum(), h_superClusters_EEP_size[i2]->GetMaximum()*s);
  h_superClusters_EEP_size[i1]->GetYaxis()->SetRangeUser(0.01,maxY*10);
  h_superClusters_EEP_size[i1]->SetLineColor(icolor1);
  h_superClusters_EEP_size[i1]->Draw(); 
  h_superClusters_EEP_size[i2]->Scale(s);
  h_superClusters_EEP_size[i2]->SetLineColor(icolor2);
  h_superClusters_EEP_size[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EEP_size.%s",dirName,fileType); 
  cSC[1]->Print(name);
  
  cSC[2]->cd();
  cSC[2]->SetLogy();
  h_superClusters_EEM_size[i1]->SetTitle("Number of SuperClusters (EE-)");
  h_superClusters_EEM_size[i1]->GetXaxis()->SetTitle("number of superclusters");
  h_superClusters_EEM_size[i1]->GetXaxis()->SetRangeUser(0,nSCMax);
  maxY = max(h_superClusters_EEM_size[i1]->GetMaximum(), h_superClusters_EEM_size[i2]->GetMaximum()*s);
  h_superClusters_EEM_size[i1]->GetYaxis()->SetRangeUser(0.01,maxY*10);
  h_superClusters_EEM_size[i1]->SetLineColor(icolor1);
  h_superClusters_EEM_size[i1]->Draw(); 
  h_superClusters_EEM_size[i2]->Scale(s);
  h_superClusters_EEM_size[i2]->SetLineColor(icolor2);
  h_superClusters_EEM_size[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EEM_size.%s",dirName,fileType); 
  cSC[2]->Print(name);



  // number of crystals/SC
  cSC[3]->cd();
  h_superClusters_EB_nXtals[i1]->SetTitle("Number of crystals in SC (EB)");
  h_superClusters_EB_nXtals[i1]->GetXaxis()->SetTitle("number of crystals");
  h_superClusters_EB_nXtals[i1]->GetXaxis()->SetRangeUser(0,nXtalsMax);
  maxY = max(h_superClusters_EB_nXtals[i1]->GetMaximum(), h_superClusters_EB_nXtals[i2]->GetMaximum()*s);
  h_superClusters_EB_nXtals[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EB_nXtals[i1]->SetLineColor(icolor1);
  h_superClusters_EB_nXtals[i1]->Draw(); 
  h_superClusters_EB_nXtals[i2]->Scale(s);
  h_superClusters_EB_nXtals[i2]->SetLineColor(icolor2);
  h_superClusters_EB_nXtals[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EB_nXtals.%s",dirName,fileType); 
  cSC[3]->Print(name);

  cSC[4]->cd();
  h_superClusters_EEP_nXtals[i1]->SetTitle("Number of crystals in SC (EE+)");
  h_superClusters_EEP_nXtals[i1]->GetXaxis()->SetTitle("number of crystals");
  h_superClusters_EEP_nXtals[i1]->GetXaxis()->SetRangeUser(0,nXtalsMax);
  maxY = max(h_superClusters_EEP_nXtals[i1]->GetMaximum(), h_superClusters_EEP_nXtals[i2]->GetMaximum()*s);
  h_superClusters_EEP_nXtals[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EEP_nXtals[i1]->SetLineColor(icolor1);
  h_superClusters_EEP_nXtals[i1]->Draw(); 
  h_superClusters_EEP_nXtals[i2]->Scale(s);
  h_superClusters_EEP_nXtals[i2]->SetLineColor(icolor2);
  h_superClusters_EEP_nXtals[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EEP_nXtals.%s",dirName,fileType); 
  cSC[4]->Print(name);
  
  cSC[5]->cd();
  h_superClusters_EEM_nXtals[i1]->SetTitle("Number of crystals in SC (EE-)");
  h_superClusters_EEM_nXtals[i1]->GetXaxis()->SetTitle("number of crystals");
  h_superClusters_EEM_nXtals[i1]->GetXaxis()->SetRangeUser(0,nXtalsMax);
  maxY = max(h_superClusters_EEM_nXtals[i1]->GetMaximum(), h_superClusters_EEM_nXtals[i2]->GetMaximum()*s);
  h_superClusters_EEM_nXtals[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EEM_nXtals[i1]->SetLineColor(icolor1);
  h_superClusters_EEM_nXtals[i1]->Draw(); 
  h_superClusters_EEM_nXtals[i2]->Scale(s);
  h_superClusters_EEM_nXtals[i2]->SetLineColor(icolor2);
  h_superClusters_EEM_nXtals[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EEM_nXtals.%s",dirName,fileType); 
  cSC[5]->Print(name);

  // number of BC per SC
  cSC[6]->cd();
  h_superClusters_EB_nBC[i1]->SetTitle("Number of Basic Clusters in SC (EB)");
  h_superClusters_EB_nBC[i1]->GetXaxis()->SetTitle("number of BC");
  h_superClusters_EB_nBC[i1]->GetXaxis()->SetRangeUser(0,nBCMax);
  maxY = max(h_superClusters_EB_nBC[i1]->GetMaximum(), h_superClusters_EB_nBC[i2]->GetMaximum()*s);
  h_superClusters_EB_nBC[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EB_nBC[i1]->SetLineColor(icolor1);
  h_superClusters_EB_nBC[i1]->Draw(); 
  h_superClusters_EB_nBC[i2]->Scale(s);
  h_superClusters_EB_nBC[i2]->SetLineColor(icolor2);
  h_superClusters_EB_nBC[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EB_nBC.%s",dirName,fileType); 
  cSC[6]->Print(name);

  cSC[7]->cd();
  h_superClusters_EEP_nBC[i1]->SetTitle("Number of Basic Clusters in SC (EE+)");
  h_superClusters_EEP_nBC[i1]->GetXaxis()->SetTitle("number of BC");
  h_superClusters_EEP_nBC[i1]->GetXaxis()->SetRangeUser(0,nBCMax);
  maxY = max(h_superClusters_EEP_nBC[i1]->GetMaximum(), h_superClusters_EEP_nBC[i2]->GetMaximum()*s);
  h_superClusters_EEP_nBC[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EEP_nBC[i1]->SetLineColor(icolor1);
  h_superClusters_EEP_nBC[i1]->Draw(); 
  h_superClusters_EEP_nBC[i2]->Scale(s);
  h_superClusters_EEP_nBC[i2]->SetLineColor(icolor2);
  h_superClusters_EEP_nBC[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EEP_nBC.%s",dirName,fileType); 
  cSC[7]->Print(name);
  
  cSC[8]->cd();
  h_superClusters_EEM_nBC[i1]->SetTitle("Number of Basic Clusters in SC (EE-)");
  h_superClusters_EEM_nBC[i1]->GetXaxis()->SetTitle("number of BC");
  h_superClusters_EEM_nBC[i1]->GetXaxis()->SetRangeUser(0,nBCMax);
  maxY = max(h_superClusters_EEM_nBC[i1]->GetMaximum(), h_superClusters_EEM_nBC[i2]->GetMaximum()*s);
  h_superClusters_EEM_nBC[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EEM_nBC[i1]->SetLineColor(icolor1);
  h_superClusters_EEM_nBC[i1]->Draw(); 
  h_superClusters_EEM_nBC[i2]->Scale(s);
  h_superClusters_EEM_nBC[i2]->SetLineColor(icolor2);
  h_superClusters_EEM_nBC[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EEM_nBC.%s",dirName,fileType); 
  cSC[8]->Print(name);


 
  // SC eta/phi
  cSC[9]->cd();
  h_superClusters_eta[i1]->SetTitle("SuperClusters eta (EB+EE)");
  h_superClusters_eta[i1]->GetXaxis()->SetTitle("#eta");
  h_superClusters_eta[i1]->SetLineColor(icolor1);
  maxY = max(h_superClusters_eta[i1]->GetMaximum(), h_superClusters_eta[i2]->GetMaximum()*s);
  h_superClusters_eta[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_eta[i1]->Draw(); 
  h_superClusters_eta[i2]->Scale(s);
  h_superClusters_eta[i2]->SetLineColor(icolor2);
  h_superClusters_eta[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_eta.%s",dirName,fileType); 
  cSC[9]->Print(name);

  cSC[10]->cd();
  h_superClusters_EB_phi[i1]->SetTitle("SuperClusters phi (EB)");
  h_superClusters_EB_phi[i1]->GetXaxis()->SetTitle("#phi");
  maxY = max(h_superClusters_EB_phi[i1]->GetMaximum(), h_superClusters_EB_phi[i2]->GetMaximum()*s);
  h_superClusters_EB_phi[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EB_phi[i1]->SetLineColor(icolor1);
  h_superClusters_EB_phi[i1]->Draw(); 
  h_superClusters_EB_phi[i2]->Scale(s);
  h_superClusters_EB_phi[i2]->SetLineColor(icolor2);
  h_superClusters_EB_phi[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EB_phi.%s",dirName,fileType); 
  cSC[10]->Print(name);

  cSC[11]->cd();
  h_superClusters_EE_phi[i1]->SetTitle("SuperClusters phi (EE)");
  h_superClusters_EE_phi[i1]->GetXaxis()->SetTitle("#phi");
  maxY = max(h_superClusters_EE_phi[i1]->GetMaximum(), h_superClusters_EE_phi[i2]->GetMaximum()*s);
  h_superClusters_EE_phi[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EE_phi[i1]->SetLineColor(icolor1);
  h_superClusters_EE_phi[i1]->Draw(); 
  h_superClusters_EE_phi[i2]->Scale(s);
  h_superClusters_EE_phi[i2]->SetLineColor(icolor2);
  h_superClusters_EE_phi[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EE_phi.%s",dirName,fileType); 
  cSC[11]->Print(name);


  // SC energy
  cSC[12]->cd();
  h_superClusters_EB_energy[i1]->SetTitle("Supercluster Energy (EB)");
  h_superClusters_EB_energy[i1]->GetXaxis()->SetTitle("energy (GeV)");
  h_superClusters_EB_energy[i1]->GetXaxis()->SetRangeUser(0,eMaxSC);
  maxY = max(h_superClusters_EB_energy[i1]->GetMaximum(), h_superClusters_EB_energy[i2]->GetMaximum()*s);
  h_superClusters_EB_energy[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EB_energy[i1]->SetLineColor(icolor1);
  h_superClusters_EB_energy[i1]->Draw(); 
  h_superClusters_EB_energy[i2]->Scale(s);
  h_superClusters_EB_energy[i2]->SetLineColor(icolor2);
  h_superClusters_EB_energy[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EB_energy.%s",dirName,fileType); 
  cSC[12]->Print(name);

  cSC[13]->cd();
  h_superClusters_EEP_energy[i1]->SetTitle("Supercluster Energy (EE+)");
  h_superClusters_EEP_energy[i1]->GetXaxis()->SetTitle("energy (GeV)");
  h_superClusters_EEP_energy[i1]->GetXaxis()->SetRangeUser(0,eMaxSC);
  maxY = max(h_superClusters_EEP_energy[i1]->GetMaximum(), h_superClusters_EEP_energy[i2]->GetMaximum()*s);
  h_superClusters_EEP_energy[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EEP_energy[i1]->SetLineColor(icolor1);
  h_superClusters_EEP_energy[i1]->Draw(); 
  h_superClusters_EEP_energy[i2]->Scale(s);
  h_superClusters_EEP_energy[i2]->SetLineColor(icolor2);
  h_superClusters_EEP_energy[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EEP_energy.%s",dirName,fileType); 
  cSC[13]->Print(name);
  
  cSC[14]->cd();
  h_superClusters_EEM_energy[i1]->SetTitle("Supercluster Energy (EE-)");
  h_superClusters_EEM_energy[i1]->GetXaxis()->SetTitle("energy (GeV)");
  h_superClusters_EEM_energy[i1]->GetXaxis()->SetRangeUser(0,eMaxSC);
  maxY = max(h_superClusters_EEM_energy[i1]->GetMaximum(), h_superClusters_EEM_energy[i2]->GetMaximum()*s);
  h_superClusters_EEM_energy[i1]->GetYaxis()->SetRangeUser(0,maxY*1.1);
  h_superClusters_EEM_energy[i1]->SetLineColor(icolor1);
  h_superClusters_EEM_energy[i1]->Draw(); 
  h_superClusters_EEM_energy[i2]->Scale(s);
  h_superClusters_EEM_energy[i2]->SetLineColor(icolor2);
  h_superClusters_EEM_energy[i2]->Draw("same"); 
  sprintf(name,"%s/superClusters_EEM_energy.%s",dirName,fileType); 
  cSC[14]->Print(name);



  //ES
  cES[0]->cd();
  cES[0]->SetLogy(); 
  h_esClusters_energy_plane1[i1] -> SetTitle("Cluster Energy (plane 1)");
  h_esClusters_energy_plane1[i1] -> GetXaxis() -> SetTitle("energy (GeV)");
  h_esClusters_energy_plane1[i1] -> SetLineColor(icolor1);
  h_esClusters_energy_plane1[i1] -> Draw();
  h_esClusters_energy_plane1[i2] -> Scale(s);
  h_esClusters_energy_plane1[i2] -> SetLineColor(icolor2);
  h_esClusters_energy_plane1[i2] -> Draw("same");
  sprintf(name,"%s/clusters_ES_plane1_energy.%s",dirName,fileType);
  cES[0]->Print(name);

  cES[1]->cd();
  cES[1]->SetLogy(); 
  h_esClusters_energy_plane2[i1] -> SetTitle("Cluster Energy (plane 2)");
  h_esClusters_energy_plane2[i1] -> GetXaxis() -> SetTitle("energy (GeV)");
  h_esClusters_energy_plane2[i1] -> SetLineColor(icolor1);
  h_esClusters_energy_plane2[i1] -> Draw();
  h_esClusters_energy_plane2[i2] -> Scale(s);
  h_esClusters_energy_plane2[i2] -> SetLineColor(icolor2);
  h_esClusters_energy_plane2[i2] -> Draw("same");
  sprintf(name,"%s/clusters_ES_plane2_energy.%s",dirName,fileType);
  cES[i2]->Print(name);

  cES[2]->cd();
  cES[2]->SetLogy(); 
  h_esClusters_energy_ratio[i1] -> SetTitle("Cluster Energy ratio plane1/plane2)");
  h_esClusters_energy_ratio[i1] -> GetXaxis() -> SetTitle("ratio");
  h_esClusters_energy_ratio[i1] -> SetLineColor(icolor1);
  h_esClusters_energy_ratio[i1] -> Draw();
  h_esClusters_energy_ratio[i2] -> Scale(s);
  h_esClusters_energy_ratio[i2] -> SetLineColor(icolor2);
  h_esClusters_energy_ratio[i2] -> Draw("same");
  sprintf(name,"%s/clusters_ES_ratio_energy.%s",dirName,fileType);
  cES[2]->Print(name);

}
