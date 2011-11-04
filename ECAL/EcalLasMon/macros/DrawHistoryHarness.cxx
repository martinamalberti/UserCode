{
 
  gStyle->SetOptStat(0);
  gStyle->SetNdivisions(1020,"X");

  // sel fed
  int selected_fed = 603;
  
  // time interval
  TTimeStamp dateMin(2011, 8, 31, 0, kTRUE, 0); 
  TTimeStamp dateMax(2011, 10, 1, 0, kTRUE, 0);   

  // input file
  TFile *f = TFile::Open("EElaserAnalysis_all.root");

  char gname[100];
  int fed;
  // plots
  // profile averaging on one harness
  TProfile *las[18][20];
  TProfile *led[18][20];
  TProfile *ratio[18][20];

  float ymin, ymax;
 
  TCanvas *c[18][20];
  char cname[100];


  for (int ism = 0; ism < 18; ism++){
    if (ism < 9 ) fed = 600 + ism+1;
    if (ism >= 9) fed = 636 + ism+1 ; 

    if (fed != selected_fed) continue;
    
    for (int ih = 0; ih < 20; ih++){
      
      sprintf(gname,"p_vptpn_las_fed%d_harness%d", fed, ih+1);
      las[ism][ih] = (TProfile*)f->Get(gname);
      
      sprintf(gname,"p_vptpn_led_fed%d_harness%d", fed, ih+1);
      led[ism][ih] = (TProfile*)f->Get(gname);
      
      sprintf(gname,"p_ratio_fed%d_harness%d", fed, ih+1);
      ratio[ism][ih] = (TProfile*)f->Get(gname);

      // draw
       
      if (las[ism][ih] != NULL && led[ism][ih] != NULL){
	sprintf(cname,"c_fed%d_harness%d",fed, ih+1);
	c[ism][ih] = new TCanvas(cname,cname,1200,400);
	c[ism][ih]->SetGridx();
	c[ism][ih]->SetGridy();
	
	TLegend *leg = new TLegend(0.12,0.15,0.3,0.3);
	leg->SetFillColor(0);
	leg->AddEntry(las[ism][ih], "blue laser", "LP");
	leg->AddEntry(led[ism][ih], "blue led", "LP");
	leg->AddEntry(ratio[ism][ih],  "ratio", "LP");
	//       ymin = las[ism][ih]->GetMinimum(0.8);
	//       ymax = las[ism][ih]->GetMaximum(1.05);
	ymin = 0.85;
	ymax = 1.05;
	las[ism][ih]->GetXaxis()->SetTitle("date");
	las[ism][ih]->GetYaxis()->SetTitle("VPT/PN");
	las[ism][ih]->GetXaxis()->SetTimeFormat("%d/%m%F1970-01-01 00:00:00");
	las[ism][ih]->GetXaxis()->SetTimeDisplay(1);
	las[ism][ih]->GetXaxis()->SetRangeUser(dateMin, dateMax);
	las[ism][ih]->GetXaxis()->SetNdivisions(1020);
	las[ism][ih]->GetYaxis()->SetRangeUser(ymin,ymax);  
	las[ism][ih]->Draw();
	led[ism][ih]->Draw("same");
	ratio[ism][ih]->Draw("same");
	leg->Draw("same");
	gPad->Update();
	
	sprintf(cname,"fed%d_harness%d_Sept2011.gif",fed, ih+1);
	//c[ism][ih]->Print(cname);
      }
    }
    
  }




}
