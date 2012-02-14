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
  
  TFile *f[2];
  f[0] = TFile::Open("outputFiles/GraphsLocalEta_regression_nBC1.root");
  f[1] = TFile::Open("outputFiles/GraphsLocalEta_regression_nBCgt1.root");

  const int Ntempl = 4;
  
  TGraphErrors* g_EoP_MC[Ntempl][2];
  TGraphErrors* g_EoC_MC[Ntempl][2];
   
  TGraphErrors* g_EoP_Data[Ntempl][2];
  TGraphErrors* g_EoC_Data[Ntempl][2];

  for (int ifile = 0; ifile < 2; ifile++){ 
    for (int mod=0; mod<Ntempl; mod++){
      char histoName[100];
      
      sprintf(histoName, "gEoP_MC_mod%d", mod+1);
      g_EoP_MC[mod][ifile]   = (TGraphErrors*)f[ifile]->Get(histoName); 
      
      sprintf(histoName, "gEoC_MC_mod%d", mod+1);
      g_EoC_MC[mod][ifile]   = (TGraphErrors*)f[ifile]->Get(histoName); 
      
      sprintf(histoName, "gEoP_Data_mod%d", mod+1);
      g_EoP_Data[mod][ifile]   = (TGraphErrors*)f[ifile]->Get(histoName); 
      
      sprintf(histoName, "gEoC_Data_mod%d", mod+1);
      g_EoC_Data[mod][ifile]   = (TGraphErrors*)f[ifile]->Get(histoName); 

      int mycolor = 2;
      if (ifile == 1) mycolor = 4; 
      //-- data uncorrected
      g_EoP_Data[mod][ifile] -> SetMarkerStyle(20);
      g_EoP_Data[mod][ifile] -> SetMarkerSize(1.);
      g_EoP_Data[mod][ifile] -> SetMarkerColor(mycolor); 
      g_EoP_Data[mod][ifile] -> SetLineColor(mycolor); 

      g_EoP_Data[mod][ifile] -> GetXaxis()->SetTitle("#eta_{SC} (deg)");
      g_EoP_Data[mod][ifile] -> GetXaxis()->SetLabelSize(0.04);
      g_EoP_Data[mod][ifile] -> GetXaxis()->SetTitleOffset(0.8);
      g_EoP_Data[mod][ifile] -> GetYaxis()->SetTitle("Relative E/p scale");
      g_EoP_Data[mod][ifile] -> GetYaxis()->SetLabelSize(0.04);
      g_EoP_Data[mod][ifile] -> GetYaxis()->SetTitleOffset(1.3);

    }
  }

  TLegend *leg = new TLegend(0.12,0.70,0.35,0.87);
  leg-> SetFillColor(0);
  leg->AddEntry(g_EoP_Data[0][0],"n(BC)=1","PL");
  leg->AddEntry(g_EoP_Data[0][1],"n(BC)>1","PL");
  

  TCanvas *c[4];
  
  c[0] = new TCanvas("cmod1","cmod1");
  c[0]->SetGridx();
  c[0]->SetGridy();
  g_EoP_Data[0][0] -> GetYaxis()->SetRangeUser(0.97,1.03);
  g_EoP_Data[0][0] -> Draw("apl");
  g_EoP_Data[0][1] -> Draw("plsame");
  leg->Draw("same");

  c[1] = new TCanvas("cmod2","cmod2");
  c[1]->SetGridx();
  c[1]->SetGridy();
  g_EoP_Data[1][0] -> GetYaxis()->SetRangeUser(0.97,1.03);
  g_EoP_Data[1][0] -> Draw("apl");
  g_EoP_Data[1][1] -> Draw("plsame");
  leg->Draw("same");

  c[2] = new TCanvas("cmod3","cmod3");
  c[2]->SetGridx();
  c[2]->SetGridy();
  g_EoP_Data[2][0] -> GetYaxis()->SetRangeUser(0.97,1.03);
  g_EoP_Data[2][0] -> Draw("apl");
  g_EoP_Data[2][1] -> Draw("plsame");
  leg->Draw("same");

  c[3] = new TCanvas("cmod4","cmod4");
  c[3]->SetGridx();
  c[3]->SetGridy();
  g_EoP_Data[3][0] -> GetYaxis()->SetRangeUser(0.97,1.03);
  g_EoP_Data[3][0] -> Draw("apl");
  g_EoP_Data[3][1] -> Draw("plsame");
  leg->Draw("same");

  
  g_EoP_Data[0][0] -> Draw("ap");
  g_EoP_Data[1][0] -> Draw("psame");
  g_EoP_Data[2][0] -> Draw("psame");
  g_EoP_Data[3][0] -> Draw("psame");


  //define local corrections
  TF1 * fData[4];
  for(int ii=0;ii<4; ii++){
    char fName[80];
    sprintf(fName, "fData_%d", ii);
    //fData[ii] = new TF1(fName,"[0]+([1]*(x-0.0)-[2]*pow(x-0.0,2)-[3]*pow(x-0.0,3)-[4]*pow(x-0.0,4))",-1,1.);
    fData[ii] = new TF1(fName,"[0]+([1]*(x-0.0)-[2]*pow(x-0.0,2))",-1,1.);
    g_EoP_Data[ii][0] -> Fit(fData[ii],"");
  }
  
  // cout <<  "f->SetParameters(" << fData[0]->GetParameter(0) << "," << fData[0]->GetParameter(1) << " , " << fData[0]->GetParameter(2) << " , " << fData[0]->GetParameter(3) <<  " , " << fData[0]->GetParameter(4) << ");" << endl ;
  // cout <<  "f->SetParameters(" << fData[1]->GetParameter(0) << "," << fData[1]->GetParameter(1) << " , " << fData[1]->GetParameter(2) << " , " << fData[1]->GetParameter(3) <<  " , " << fData[1]->GetParameter(4) << ");" << endl ;
  // cout <<  "f->SetParameters(" << fData[2]->GetParameter(0) << "," << fData[2]->GetParameter(1) << " , " << fData[2]->GetParameter(2) << " , " << fData[2]->GetParameter(3) <<  " , " << fData[2]->GetParameter(4) << ");" << endl ;
  // cout <<  "f->SetParameters(" << fData[3]->GetParameter(0) << "," << fData[3]->GetParameter(1) << " , " << fData[3]->GetParameter(2) << " , " << fData[3]->GetParameter(3) <<  " , " << fData[3]->GetParameter(4) << ");" << endl ;

  cout <<  "f->SetParameters(" << fData[0]->GetParameter(0) << "," << fData[0]->GetParameter(1) << " , " << fData[0]->GetParameter(2) << ");" << endl ;
  cout <<  "f->SetParameters(" << fData[1]->GetParameter(0) << "," << fData[1]->GetParameter(1) << " , " << fData[1]->GetParameter(2) << ");" << endl ;
  cout <<  "f->SetParameters(" << fData[2]->GetParameter(0) << "," << fData[2]->GetParameter(1) << " , " << fData[2]->GetParameter(2) << ");" << endl ;
  cout <<  "f->SetParameters(" << fData[3]->GetParameter(0) << "," << fData[3]->GetParameter(1) << " , " << fData[3]->GetParameter(2) << ");" << endl ;





}
