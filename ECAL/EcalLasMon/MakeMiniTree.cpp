
// compile with:
// g++ `root-config --libs --cflags` MakeMiniTree.cpp -o MakeMiniTree

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "TDirectory.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TTimeStamp.h"


#include "./LaserDataAnalysis.h"

using namespace std;

int main(int argc, char ** argv)
{
  // select fed
  int selected_fed = atoi(argv[1]);
 
  //Get old tree 
  TChain *tx = new TChain("x");
  tx->Add("/data2/EcalLaserMonitoringData/ntuples_2011_158851_178888/ntu_data_001*.root");
 
  init_ttree(tx, &x);
  
  TTree *oldtree = (TTree*)tx;

  Long64_t nentries = oldtree->GetEntries();
  cout<< "Number of entries in the tree : " << nentries << endl;
  
  //Create a new file + a clone of old tree in new file
  char fname[100];
  sprintf(fname,"/tmp/malberti/ntu_data_fed%d.root",selected_fed);
  TFile *newfile = new TFile(fname,"recreate");
  TTree *newtree = oldtree->CloneTree(0);
  
  for (int i=0; i<nentries; i++) {
    if (i%10000000==0) cout << i << endl;
    oldtree->GetEntry(i);
    if ( x.fed == selected_fed ) newtree->Fill();
  }
  
  newtree->Print();
  newtree->AutoSave();
  
  delete newfile;


}
