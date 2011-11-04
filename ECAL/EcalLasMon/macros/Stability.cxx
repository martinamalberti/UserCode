//****  LIST OF GOOD HARNESSES                                                                                                                                                 
bool goodHarness(int fed, int harness){

  bool is_good = true;

  // bad harnesses September 2011                                                                                                                                              
  if ( fed == 601 && harness == 12  ) is_good = false; // bad laser and led. PN?
  if ( fed == 601 && harness == 13  ) is_good = false; // bad laser and led. PN?
  if ( fed == 602 && harness == 6   ) is_good = false; // bad led, laser ok 
  if ( fed == 603 && harness == 1   ) is_good = false; // bad led, laser ok 
  if ( fed == 603 && harness == 4   ) is_good = false; // 
  if ( fed == 608 && harness == 19  ) is_good = false;

  if ( fed == 646 && harness == 11  ) is_good = false;
  if ( fed == 646 && harness == 9   ) is_good = false;
  if ( fed == 649 && harness == 4   ) is_good = false;
  if ( fed == 650 && harness == 8   ) is_good = false;
  if ( fed == 651 && harness == 10  ) is_good = false;
  if ( fed == 651 && harness == 11  ) is_good = false;
  if ( fed == 652 && harness == 16  ) is_good = false;
  if ( fed == 653 && harness == 19  ) is_good = false;
  if ( fed == 654                   ) is_good = false;

  return (is_good);

}



void Stability(){
 

  //TFile *f = TFile::Open("EElaserAnalysis_all_September2011.root");
  TFile *f = TFile::Open("EElaserAnalysis_all.root");

  TH1F *hstabilityEEP  = new TH1F("hstabilityEEP","hstabilityEEP",1000,0.,0.1);
  TH1F *hstabilityEEM  = new TH1F("hstabilityEEM","hstabilityEEM",1000,0.,0.1);
 
  float meanEEP, rmsEEP;
  float meanEEM, rmsEEM;
  int fed, harness;

  for (int i = 0; i < 101; i++){
    for (int j = 0; j < 101; j++){
      meanEEP = hmapMeanEEP   -> GetCellContent(i,j);
      rmsEEP  = hmapRmsEEP    -> GetCellContent(i,j);
      fed     = fedMapEEP     -> GetCellContent(i,j);
      harness = harnessMapEEP -> GetCellContent(i,j);
      if (goodHarness(fed,harness) && rmsEEP> 0  )
	hstabilityEEP ->Fill(rmsEEP);
      else{
// 	hmapMeanEEP   -> SetCellContent(i,j, 0);
// 	hmapRmsEEP   -> SetCellContent(i,j, 0);
      }
      meanEEM = hmapMeanEEM   -> GetCellContent(i,j);
      rmsEEM  = hmapRmsEEM    -> GetCellContent(i,j);
      fed     = fedMapEEM     -> GetCellContent(i,j);
      harness = harnessMapEEM -> GetCellContent(i,j);
      if (goodHarness(fed,harness)&& rmsEEM> 0)
      	hstabilityEEM ->Fill(rmsEEM);
      else{
//  	hmapMeanEEM   -> SetCellContent(i,j, 0);
//  	hmapRmsEEM   -> SetCellContent(i,j, 0);
      }
    }
  }
  
  TCanvas *cstabEEP = new TCanvas("cstabEEP","cstabEEP");
  hstabilityEEP->Draw();

  TCanvas *cmeanEEP = new TCanvas("cmeanEEP","cmeanEEP");
  hmapMeanEEP->GetZaxis()->SetRangeUser(0.96,1.01);
  hmapMeanEEP->Draw("colz");
  drawEELines();

  TCanvas *crmsEEP = new TCanvas("crmsEEP","crmsEEP");
  hmapRmsEEP->GetZaxis()->SetRangeUser(0.,0.01);
  hmapRmsEEP->Draw("colz");
  drawEELines();

  TCanvas *cstabEEM = new TCanvas("cstabEEM","cstabEEM");
  hstabilityEEM->Draw();

  TCanvas *cmeanEEM = new TCanvas("cmeanEEM","cmeanEEM");
  hmapMeanEEM->GetZaxis()->SetRangeUser(0.96,1.01);
  hmapMeanEEM->Draw("colz");
  drawEELines();
 
  TCanvas *crmsEEM = new TCanvas("crmsEEM","crmsEEM");
  hmapRmsEEM->GetZaxis()->SetRangeUser(0.,0.01);
  hmapRmsEEM->Draw("colz");
  drawEELines();

  new TCanvas(); 
  harnessMapEEP->Draw("colz");
  new TCanvas(); 
  fedMapEEP->GetZaxis()->SetRangeUser(646,654);
  fedMapEEP->Draw("colz");

  new TCanvas(); 
  harnessMapEEM->Draw("colz");
  new TCanvas(); 
  fedMapEEM->GetZaxis()->SetRangeUser(601,609);
  fedMapEEM->Draw("colz");

}



void drawEELines() {

  int ixSectorsEE[202] = {61, 61, 60, 60, 59, 59, 58, 58, 57, 57, 55, 55, 45, 45, 43, 43, 42, 42, 41, 41, 40, 40, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 45, 45, 55, 55, 57, 57, 58, 58, 59, 59, 60, 60, 61, 61, 0,100,100, 97, 97, 95, 95, 92, 92, 87, 87, 85, 85, 80, 80, 75, 75, 65, 65, 60, 60, 40, 40, 35, 35, 25, 25, 20, 20, 15, 15, 13, 13,  8,  8,  5,  5,  3,  3,  0,  0,  3,  3,  5,  5,  8,  8, 13, 13, 15, 15, 20, 20, 25, 25, 35, 35, 40, 40, 60, 60, 65, 65, 75, 75, 80, 80, 85, 85, 87, 87, 92, 92, 95, 95, 97, 97,100,100,  0, 61, 65, 65, 70, 70, 80, 80, 90, 90, 92,  0, 61, 65, 65, 90, 90, 97,  0, 57, 60, 60, 65, 65, 70, 70, 75, 75, 80, 80,  0, 50, 50,  0, 43, 40, 40, 35, 35, 30, 30, 25, 25, 20, 20,  0, 39, 35, 35, 10, 10,  3,  0, 39, 35, 35, 30, 30, 20, 20, 10, 10,  8,  0, 45, 45, 40, 40, 35, 35,  0, 55, 55, 60, 60, 65, 65};
 
  int iySectorsEE[202] = {50, 55, 55, 57, 57, 58, 58, 59, 59, 60, 60, 61, 61, 60, 60, 59, 59, 58, 58, 57, 57, 55, 55, 45, 45, 43, 43, 42, 42, 41, 41, 40, 40, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 45, 45, 50,  0, 50, 60, 60, 65, 65, 75, 75, 80, 80, 85, 85, 87, 87, 92, 92, 95, 95, 97, 97,100,100, 97, 97, 95, 95, 92, 92, 87, 87, 85, 85, 80, 80, 75, 75, 65, 65, 60, 60, 40, 40, 35, 35, 25, 25, 20, 20, 15, 15, 13, 13,  8,  8,  5,  5,  3,  3,  0,  0,  3,  3,  5,  5,  8,  8, 13, 13, 15, 15, 20, 20, 25, 25, 35, 35, 40, 40, 50,  0, 45, 45, 40, 40, 35, 35, 30, 30, 25, 25,  0, 50, 50, 55, 55, 60, 60,  0, 60, 60, 65, 65, 70, 70, 75, 75, 85, 85, 87,  0, 61,100,  0, 60, 60, 65, 65, 70, 70, 75, 75, 85, 85, 87,  0, 50, 50, 55, 55, 60, 60,  0, 45, 45, 40, 40, 35, 35, 30, 30, 25, 25,  0, 39, 30, 30, 15, 15,  5,  0, 39, 30, 30, 15, 15,  5};

 TLine l;
 l.SetLineWidth(1);
 for ( int i=0; i<201; i=i+1) {
   if ( (ixSectorsEE[i]!=0 || iySectorsEE[i]!=0) && 
   	(ixSectorsEE[i+1]!=0 || iySectorsEE[i+1]!=0) ) {
     l.DrawLine(ixSectorsEE[i], iySectorsEE[i], 
		ixSectorsEE[i+1], iySectorsEE[i+1]);
   }
 }


}
