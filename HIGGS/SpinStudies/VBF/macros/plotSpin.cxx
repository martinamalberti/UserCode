
void plotLLR(TCanvas *c){
  
  TH1F *hLLR0 = (TH1F*)_file0->Get("hLLR0_nsigcut_4"); 
  TH1F *hLLR2 = (TH1F*)_file0->Get("hLLR2_nsigcut_4"); 

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  hLLR0->Rebin(2);
  hLLR2->Rebin(2);
  hLLR2->GetXaxis()->SetTitle("LLR");
  hLLR2->GetYaxis()->SetTitle("a.u.");
  hLLR2->GetYaxis()->SetTitleOffset(1.3);
  hLLR2->GetXaxis()->SetRangeUser(-30,30);
  hLLR2->SetMaximum(1.3*hLLR2->GetMaximum());
  hLLR2->DrawNormalized();
  hLLR0->DrawNormalized("same");
  
  TLegend *l = new TLegend(0.68,0.7,0.89,0.89);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->AddEntry(hLLR0,"S = 0","F");
  l->AddEntry(hLLR2,"S = 2","F");
  l->Draw("same");
}

void plotNsig(TH1F *h, TCanvas *c){
  gStyle->SetOptStat(1110);
  h->GetXaxis()->SetTitle("n_{sig}");
  h->GetXaxis()->SetRangeUser(0,50);
  h->Draw();
  c->Print("",".gif");
}

void plotNbkg(TH1F *h, TCanvas *c){
  gStyle->SetOptStat(1110);
  h->GetXaxis()->SetTitle("n_{bkg}");
  h->GetXaxis()->SetRangeUser(200,500);
  h->Draw();
  c->Print("",".gif");
}

void plotSpin(){
  
  //-- LLR
  TCanvas *cLLR = new TCanvas("cLLR","cLLR",600,600);
  plotLLR(cLLR);
  cLLR->Print("",".gif");

  //-- nsig fitted
  TH1F *hnsig0_m0 = (TH1F*)_file0->Get("hNsig0_m0"); 
  TH1F *hnsig2_m0 = (TH1F*)_file0->Get("hNsig2_m0"); 
  TH1F *hnsig0_m2 = (TH1F*)_file0->Get("hNsig0_m2"); 
  TH1F *hnsig2_m2 = (TH1F*)_file0->Get("hNsig2_m2"); 

  TCanvas *cnsig0_m0 = new TCanvas("cnsig0_m0","cnsig0_m0",600,600);
  plotNsig(hnsig0_m0,cnsig0_m0);
 
  TCanvas *cnsig2_m0 = new TCanvas("cnsig2_m0","cnsig2_m0",600,600);
  plotNsig(hnsig2_m0,cnsig2_m0);

  TCanvas *cnsig0_m2 = new TCanvas("cnsig0_m2","cnsig0_m2",600,600);
  plotNsig(hnsig0_m2,cnsig0_m2);
 
  TCanvas *cnsig2_m2 = new TCanvas("cnsig2_m2","cnsig2_m2",600,600);
  plotNsig(hnsig2_m2,cnsig2_m2);

  //-- nbkg fitted
  TH1F *hnbkg0_m0 = (TH1F*)_file0->Get("hNbkg0_m0"); 
  TH1F *hnbkg2_m0 = (TH1F*)_file0->Get("hNbkg2_m0"); 
  TH1F *hnbkg0_m2 = (TH1F*)_file0->Get("hNbkg0_m2"); 
  TH1F *hnbkg2_m2 = (TH1F*)_file0->Get("hNbkg2_m2"); 

  TCanvas *cnbkg0_m0 = new TCanvas("cnbkg0_m0","cnbkg0_m0",600,600);
  plotNbkg(hnbkg0_m0,cnbkg0_m0);
 
  TCanvas *cnbkg2_m0 = new TCanvas("cnbkg2_m0","cnbkg2_m0",600,600);
  plotNbkg(hnbkg2_m0,cnbkg2_m0);

  TCanvas *cnbkg0_m2 = new TCanvas("cnbkg0_m2","cnbkg0_m2",600,600);
  plotNbkg(hnbkg0_m2,cnbkg0_m2);
 
  TCanvas *cnbkg2_m2 = new TCanvas("cnbkg2_m2","cnbkg2_m2",600,600);
  plotNbkg(hnbkg2_m2,cnbkg2_m2);
  
}
