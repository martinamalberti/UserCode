/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVARegressionApplication                                          *
 *                                                                                *
 * This macro provides a simple example on how to use the trained regression MVAs *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TStopwatch.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTimeStamp.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

#define NWL 3 // if you change this, mind to change the Branch declaration of TTree...

#define ETA_NBIN 100
#define ETA_MIN   -3.
#define ETA_MAX   +3.
   struct ntu_xtals {
    int   run;
    int   seq;
    int   fed;
    int   ix;
    int   iy;
    int   iz;
    int   detId;
    int   elecId;
    int   harness;
    int   side;
    float eta;
    float phi;
    float field;
    float alpha;
    int   wl[NWL];
    int   time[NWL];
    int   nevt[NWL];
    int   laser_power[NWL];
    int   vinj[NWL];
    int   nxtal25[NWL];
    int   npn[NWL];
    int   vfe_mode[NWL];
    int   mem_gain[NWL];
    int   numfill[NWL];
    float lumi[NWL];
    float qmax[NWL];
    float qmaxS[NWL];
    float tmax[NWL];
    float apdpnA[NWL];
    float apdpnAS[NWL];
    float apdpnB[NWL];
    float apdpnBS[NWL];
    float apdpnAB[NWL];
    float apdpnABS[NWL];
    float ped[NWL];
    float pedS[NWL];
    float ped_tp[NWL];
    float ped_tpS[NWL];
    float ped_laser[NWL];
    float ped_laserS[NWL];
    float tp[NWL];
    float tpS[NWL];
    float corrwidth[NWL];
    float apdapdA[NWL];
    float apdapdAS[NWL];
    float apdapdB[NWL];
    float apdapdBS[NWL];
    float l_ampli[NWL];
    float l_rise_time[NWL];
    float l_fwhm[NWL];
    float l_width90[NWL];
    float l_width95[NWL];
    float l_prepulse[NWL];
    float l_nphot0[NWL];
    float l_nphot1[NWL];
    float l_nphot2[NWL];
    float l_nphot3[NWL];
    float l_nphot4[NWL];
    float pnA_qmax[NWL];
    float pnA_qmaxS[NWL];
    float pnA_pnpnB[NWL];
    float pnA_pnpnBS[NWL];
    float pnA_corrwidth[NWL];
    float pnA_ped[NWL];
    float pnA_pedS[NWL];
    float pnA_tp_ped[NWL];
    float pnA_tp_pedS[NWL];
    float pnA_l_ped[NWL];
    float pnA_l_pedS[NWL];
    float pnA_tp[NWL];
    float pnA_tpS[NWL];
    float pnB_qmax[NWL];
    float pnB_qmaxS[NWL];
    float pnB_pnpnB[NWL];
    float pnB_pnpnBS[NWL];
    float pnB_corrwidth[NWL];
    float pnB_ped[NWL];
    float pnB_pedS[NWL];
    float pnB_tp_ped[NWL];
    float pnB_tp_pedS[NWL];
    float pnB_l_ped[NWL];
    float pnB_l_pedS[NWL];
    float pnB_tp[NWL];
    float pnB_tpS[NWL];
};


struct ntu_xtals x;

void init_ttree(TTree * t, struct ntu_xtals * x)
{
    t->SetBranchAddress("run", &x->run);
    t->SetBranchAddress("seq", &x->seq);
    t->SetBranchAddress("fed", &x->fed);
    t->SetBranchAddress("ix", &x->ix);
    t->SetBranchAddress("iy", &x->iy);
    t->SetBranchAddress("iz", &x->iz);
    t->SetBranchAddress("detId", &x->detId);
    t->SetBranchAddress("elecId", &x->elecId);
    t->SetBranchAddress("harness", &x->harness);
    t->SetBranchAddress("side", &x->side);
    t->SetBranchAddress("eta", &x->eta);
    t->SetBranchAddress("phi", &x->phi);
    t->SetBranchAddress("field", &x->field);
    t->SetBranchAddress("alpha", &x->alpha);
    t->SetBranchAddress("wl", &x->wl);
    t->SetBranchAddress("time", &x->time);
    t->SetBranchAddress("nevt", &x->nevt);
    t->SetBranchAddress("laser_power", &x->laser_power);
    t->SetBranchAddress("vinj", &x->vinj);
    t->SetBranchAddress("nxtal25", &x->nxtal25);
    t->SetBranchAddress("npn", &x->npn);
    t->SetBranchAddress("vfe_mode", &x->vfe_mode);
    t->SetBranchAddress("mem_gain", &x->mem_gain);
    t->SetBranchAddress("numfill", &x->numfill);
    t->SetBranchAddress("lumi", &x->lumi);
    t->SetBranchAddress("qmax", &x->qmax);
    t->SetBranchAddress("qmaxS", &x->qmaxS);
    t->SetBranchAddress("tmax", &x->tmax);
    t->SetBranchAddress("apdpnA", &x->apdpnA);
    t->SetBranchAddress("apdpnAS", &x->apdpnAS);
    t->SetBranchAddress("apdpnB", &x->apdpnB);
    t->SetBranchAddress("apdpnBS", &x->apdpnBS);
    t->SetBranchAddress("apdpnAB", &x->apdpnAB);
    t->SetBranchAddress("apdpnABS", &x->apdpnABS);
    t->SetBranchAddress("ped", &x->ped);
    t->SetBranchAddress("pedS", &x->pedS);
    t->SetBranchAddress("ped_tp", &x->ped_tp);
    t->SetBranchAddress("ped_tpS", &x->ped_tpS);
    t->SetBranchAddress("ped_laser", &x->ped_laser);
    t->SetBranchAddress("ped_laserS", &x->ped_laserS);
    t->SetBranchAddress("tp", &x->tp);
    t->SetBranchAddress("tpS", &x->tpS);
    t->SetBranchAddress("corrwidth", &x->corrwidth);
    t->SetBranchAddress("apdapdA", &x->apdapdA);
    t->SetBranchAddress("apdapdAS", &x->apdapdAS);
    t->SetBranchAddress("apdapdB", &x->apdapdB);
    t->SetBranchAddress("apdapdBS", &x->apdapdBS);
    t->SetBranchAddress("l_ampli", &x->l_ampli);
    t->SetBranchAddress("l_rise_time", &x->l_rise_time);
    t->SetBranchAddress("l_fwhm", &x->l_fwhm);
    t->SetBranchAddress("l_width90", &x->l_width90);
    t->SetBranchAddress("l_width95", &x->l_width95);
    t->SetBranchAddress("l_prepulse", &x->l_prepulse);
    t->SetBranchAddress("l_nphot0", &x->l_nphot0);
    t->SetBranchAddress("l_nphot1", &x->l_nphot1);
    t->SetBranchAddress("l_nphot2", &x->l_nphot2);
    t->SetBranchAddress("l_nphot3", &x->l_nphot3);
    t->SetBranchAddress("l_nphot4", &x->l_nphot4);
    t->SetBranchAddress("pnA_qmax", &x->pnA_qmax);
    t->SetBranchAddress("pnA_qmaxS", &x->pnA_qmaxS);
    t->SetBranchAddress("pnA_pnpnB", &x->pnA_pnpnB);
    t->SetBranchAddress("pnA_pnpnBS", &x->pnA_pnpnBS);
    t->SetBranchAddress("pnA_corrwidth", &x->pnA_corrwidth);
    t->SetBranchAddress("pnA_ped", &x->pnA_ped);
    t->SetBranchAddress("pnA_pedS", &x->pnA_pedS);
    t->SetBranchAddress("pnA_tp_ped", &x->pnA_tp_ped);
    t->SetBranchAddress("pnA_tp_pedS", &x->pnA_tp_pedS);
    t->SetBranchAddress("pnA_l_ped", &x->pnA_l_ped);
    t->SetBranchAddress("pnA_l_pedS", &x->pnA_l_pedS);
    t->SetBranchAddress("pnA_tp", &x->pnA_tp);
    t->SetBranchAddress("pnA_tpS", &x->pnA_tpS);
    t->SetBranchAddress("pnB_qmax", &x->pnB_qmax);
    t->SetBranchAddress("pnB_qmaxS", &x->pnB_qmaxS);
    t->SetBranchAddress("pnB_pnpnB", &x->pnB_pnpnB);
    t->SetBranchAddress("pnB_pnpnBS", &x->pnB_pnpnBS);
    t->SetBranchAddress("pnB_corrwidth", &x->pnB_corrwidth);
    t->SetBranchAddress("pnB_ped", &x->pnB_ped);
    t->SetBranchAddress("pnB_pedS", &x->pnB_pedS);
    t->SetBranchAddress("pnB_tp_ped", &x->pnB_tp_ped);
    t->SetBranchAddress("pnB_tp_pedS", &x->pnB_tp_pedS);
    t->SetBranchAddress("pnB_l_ped", &x->pnB_l_ped);
    t->SetBranchAddress("pnB_l_pedS", &x->pnB_l_pedS);
    t->SetBranchAddress("pnB_tp", &x->pnB_tp);
    t->SetBranchAddress("pnB_tpS", &x->pnB_tpS);
}

int ebwl[] = { 440, 800 };
int eewl[] = { 440, 455, 617 };

int isEB(int ifed)
{
    return (ifed >= 610 && ifed <= 645);
}

using namespace TMVA;

void TMVARegressionApplication( TString myMethodList = "" ) 
{

  gROOT->SetStyle("Plain");

   //---------------------------------------------------------------
   // default MVA methods to be trained + tested

   // this loads the library
   TMVA::Tools::Instance();

   std::map<std::string,int> Use;

   Use["PDERS"]           = 0;
   Use["PDERSkNN"]        = 0; 
   Use["PDEFoam"]         = 0;
   // ---
   Use["KNN"]             = 0;
   // ---
   Use["LD"]		  = 0;
   // ---
   Use["FDA_GA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   // ---
   Use["MLP"]             = 0; 
   // ---
   Use["SVM"]             = 0;
   // ---
   Use["BDT"]             = 0;
   Use["BDTG"]            = 1;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegressionApplication" << std::endl;

   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   //
   // create the Reader object
   //
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // create a set of variables and declare them to the reader
   // - the variable names must corresponds in name and type to 
   // those given in the weight file(s) that you use
   Float_t var1, var2, var3, var4, var5;
   reader->AddVariable( "qmax[0]", &var1 );
   reader->AddVariable( "tmax[0]", &var2 );
   reader->AddVariable( "l_fwhm[0]", &var3);
   reader->AddVariable( "l_prepulse[0]", &var4 );
   reader->AddVariable( "l_width90[0]", &var5 );

   //Spectator variables declared in the training have to be added to the reader, too
   //  Float_t spec1;
   //   reader->AddSpectator( "spec1:=eventn.tStamp",  &spec1 );

   //
   // book the MVA methods
   //
   TString dir    = "weights/";
   TString prefix = "TMVARegression";

   // book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = it->first + " method";
	 TString weightfile = dir + prefix + "_" + TString(it->first) + ".weights.xml";
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   // example how to use your own method as plugin
   if (Use["Plugin"]) {
      // the weight file contains a line 
      // Method         : MethodName::InstanceName

      // if MethodName is not a known TMVA method, it is assumed to be
      // a user implemented method which has to be loaded via the
      // plugin mechanism
      
      // for user implemented methods the line in the weight file can be
      // Method         : PluginName::InstanceName
      // where PluginName can be anything

      // before usage the plugin has to be defined, which can happen
      // either through the following line in .rootrc:
      // # plugin handler          plugin       class            library        constructor format
      // Plugin.TMVA@@MethodBase:  PluginName   MethodClassName  UserPackage    "MethodName(DataSet&,TString)"
      //  
      // or by telling the global plugin manager directly
      gPluginMgr->AddHandler("TMVA@@MethodBase", "PluginName", "MethodClassName", "UserPackage", "MethodName(DataSet&,TString)");
      // the class is then looked for in libUserPackage.so

      // now the method can be booked like any other
      reader->BookMVA( "User method", dir + prefix + "_User.weights.txt" );
   }

   // book output histograms
   TH1* h_regOutput = new TH1F("h_regOutput","regression output", 1000, -2, 2 );
   h_regOutput->SetFillStyle(3005); 
   h_regOutput->SetFillColor(kGreen);;
  
   TGraph *g_regOutput_vs_time = new TGraph();
   g_regOutput_vs_time -> SetName("g_regOutput_vs_time");
   g_regOutput_vs_time -> SetTitle("g_regOutput_vs_time");

   TGraph *g_regOutput_vs_apd = new TGraph();
   g_regOutput_vs_apd -> SetName("g_regOutput_vs_apd");
   g_regOutput_vs_apd -> SetTitle("g_regOutput_vs_apd");

   TGraph *g_regOutput_vs_width = new TGraph();
   g_regOutput_vs_width -> SetName("g_regOutput_vs_width");
   g_regOutput_vs_width -> SetTitle("g_regOutput_vs_width");
 
   
   TH1 *hLas = new TH1F("hLas","",500,-2,2.); 
   hLas->SetFillStyle(3004); 
   hLas->SetFillColor(kRed); 

   TGraph *gLas[1700];  
   for (int iApd = 0;iApd<1700;iApd++){
     gLas[iApd] = new TGraph(); 
     gLas[iApd]->SetMarkerColor(kRed); 
     char gName[80];
     sprintf(gName,"ApdPN_%04d",iApd);
     gLas[iApd]->SetName(gName);
     gLas[iApd]->SetTitle(gName);
   }
   
 
   TChain * theTree = new TChain("x");
   theTree->Add("/data2/EcalLaserMonitoringData/ntuples_2011_158851_178888/ntu_data_001*.root");

   init_ttree(theTree, &x);
   theTree->SetBranchStatus("*",0); //disable all branches
   // only read needed branches
   theTree->SetBranchStatus("run",1);
   theTree->SetBranchStatus("fed",1);
   theTree->SetBranchStatus("field",1);
   theTree->SetBranchStatus("apdpnAB",1);
   theTree->SetBranchStatus("time",1);
   theTree->SetBranchStatus("harness",1);
   theTree->SetBranchStatus("qmax",1);
   theTree->SetBranchStatus("tmax",1);
   theTree->SetBranchStatus("l_fwhm",1);
   theTree->SetBranchStatus("l_prepulse",1);
   theTree->SetBranchStatus("l_width90",1);
   
   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   int evt[1700]={0};
   
   for (Long64_t ientry=0; ientry<theTree->GetEntries();ientry++) {
     
     if (ientry%1000000 == 0) {
       std::cout << "--- ... Processing entry: " << ientry << std::endl;
     }
          
     theTree->GetEntry(ientry);
     
     //      if (timeStamp > 1302e+6) continue; 
     if (x.run < 161200) continue;
     if (x.field < 0.5) continue;
     if (x.fed != 605) continue;
     if (x.harness != 6) continue;
     if (x.apdpnAB[0]<=0) continue;
     int t = x.time[0];

     // 
     // retrieve the MVA target values (regression outputs) and fill into histograms
     // NOTE: EvaluateRegression(..) returns a vector for multi-target regression
     // 
     var1 = x.qmax[0];
     var2 = x.tmax[0];
     var3 = x.l_fwhm[0];
     var4 = x.l_prepulse[0];
     var5 = x.l_width90[0]; 
 
     iApd = x.elecId;     
 
     double val = (reader->EvaluateRegression( "BDTG method" ))[0];
     g_regOutput_vs_time->SetPoint(evt[iApd],t,val);
     h_regOutput->Fill( x.apdpnAB[0] - val );    
     
     gLas[iApd] -> SetPoint(evt[iApd],t,x.apdpnAB[0]); 
     hLas -> Fill( x.apdpnAB[0] - 1 ); 
     
     evt[iApd]++;

   }  // end loop over entries

   


   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
    
   // Plot results:
   TTimeStamp dateMin(2011, 2,  20, 0, kTRUE, 0); 
   TTimeStamp dateMax(2011, 10, 30, 0, kTRUE, 0); 
   
   nhists=1;

   TH1F *hC; 
   TCanvas *c = new TCanvas("c","c");
   gPad->SetGridx(); gPad->SetGridy();
   hC = (TH1F*)gPad->DrawFrame(dateMin, 0, dateMax, 3); 
   hC->GetXaxis()->SetTimeDisplay(1); 
   hC->GetXaxis()->SetTitle("date"); 
   hC->GetYaxis()->SetTitle("APD/PN");
   for (int iApd=0;iApd<1700;iApd++)
     if (gLas[iApd]->GetN() > 0) gLas[iApd]->Draw("PSAME");
   g_regOutput_vs_time->Draw("PSAME");
   

     
   //save histograms
   
   TFile *target  = new TFile( "TMVARegApp.root","RECREATE" );
   
   h_regOutput->Write();
   g_regOutput_vs_time->Write();

   for (int iApd=0;iApd<1700;iApd++)
     if (gLas[iApd]->GetN() > 0)  gLas[iApd]->Write();
   hLas->Write(); 
   
   target->Close();
   
   std::cout << "--- Created root file: \"" << target->GetName() 
	     << "\" containing the MVA output histograms" << std::endl;
   
   delete reader;
   
   std::cout << "==> TMVARegressionApplication is done!" << std::endl << std::endl;

   
}
