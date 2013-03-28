#include "ntpleUtils.h"



std::map<int, int> GetTotalEvents(const std::string& histoName, const std::string& inputFileList)
{
  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;
  std::map<int, int> totalEvents;
  
  if(!inFile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return totalEvents;
  }
  
  while(1)
  {
    inFile >> buffer;
    if(!inFile.good()) break;

    TFile* f = new TFile(buffer.c_str(), "READ");
    TH1F* histo = (TH1F*)(f -> Get(histoName.c_str()));
    
    for(int bin = 1; bin <= histo -> GetNbinsX(); ++bin)
      totalEvents[bin] += int(histo -> GetBinContent(bin));
    
    f -> Close();
    
    delete f;
  }

  return totalEvents;
}

//  ------------------------------------------------------------


bool FillChain(TChain& chain, const std::string& inputFileList)
{
  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;

  if(!inFile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return false;
  }
  
  while(1)
  {
    inFile >> buffer;
    if(!inFile.good()) break;
    chain.Add(buffer.c_str());
  }

  return true;
}

//  ------------------------------------------------------------

int parseConfigFile (const TString& config)
{
  std::cout << ">>> Parsing " << config << " file" << std::endl ;
  
  if (gConfigParser) return 1 ;
  gConfigParser = new ConfigParser();
  
  if( !(gConfigParser -> init(config)) )
  {
    std::cout << ">>> parseConfigFile::Could not open configuration file "
              << config << std::endl;
     return -1;
  }
  
  gConfigParser -> print();
  
  return 0 ;
}

//  ------------------------------------------------------------






double deltaPhi(const double& phi1, const double& phi2)
{ 
  double deltaphi = fabs(phi1 - phi2);
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308 - deltaphi;
  return deltaphi;
}

//  ------------------------------------------------------------

double deltaEta(const double& eta1, const double& eta2)
{ 
  double deltaeta = fabs(eta1 - eta2);
  return deltaeta;
}

//  ------------------------------------------------------------

double deltaR(const double& eta1, const double& phi1,
              const double& eta2, const double& phi2)
{ 
 double deltaphi = deltaPhi(phi1, phi2);
 double deltaeta = deltaEta(eta1, eta2);
 double deltar = sqrt(deltaphi*deltaphi + deltaeta*deltaeta);
 return deltar;
}

//  ------------------------------------------------------------



int getCJV(std::vector<ROOT::Math::XYZTVector>& jets,
	      int q1,
	      int q2,
	      const double& EtMin,
	      const std::vector<int>* blacklist){
 
 int CJV = 0;
 double etaMin = jets.at(q1).Eta();
 double etaMax = jets.at(q2).Eta();
 
 if (etaMax < etaMin) std::swap(etaMin,etaMax);
 
 for(unsigned int i = 0; i < jets.size(); ++i)
 {
  if (i==q1 || i==q2) continue;
  if (jets.at(i).Et() < EtMin) continue;
  
  bool skipJet = false;
  if(blacklist){
   for(unsigned int kk = 0; kk < blacklist -> size(); ++kk) {
    if(blacklist -> at(kk) == static_cast<int>(i)) skipJet = true;
   }
  }
  
  if (skipJet) continue;
   
  if(jets.at(i).Eta() > etaMax || jets.at(i).Eta() < etaMin) continue;
  
  CJV++;
 } 
 
 return CJV;
 
}


//  ------------------------------------------------------------


int getJV(std::vector<ROOT::Math::XYZTVector>& jets,
	      const double& EtMin,
	      const std::vector<int>* blacklist){
 
 int JV = 0;
  
 for(unsigned int i = 0; i < jets.size(); ++i)
 {
  if (jets.at(i).Et() < EtMin) continue;
  
  bool skipJet = false;
  if(blacklist){
   for(unsigned int kk = 0; kk < blacklist -> size(); ++kk) {
    if(blacklist -> at(kk) == static_cast<int>(i)) skipJet = true;
   }
  }
  
  if (skipJet) continue;
  JV++;  
 } 
 
 return JV;
 
}

//  ------------------------------------------------------------

int getZepp(std::vector<ROOT::Math::XYZTVector>& jets,
	  int q1,
	  int q2,
	  const double& EtMin,
	  const double& zeppMax, 
	  const std::vector<int>* blacklist){
 
 int nJet = 0;
 
 double etaMin = jets.at(q1).Eta();
 double etaMax = jets.at(q2).Eta();

 if (etaMax < etaMin) std::swap(etaMin,etaMax);

 double etaMean = (etaMax + etaMin) / 2.;
 double dEta = (etaMax - etaMin);
 
 
 for(unsigned int i = 0; i < jets.size(); ++i)
 {
  if (i==q1 || i==q2) continue;
  if (jets.at(i).Et() < EtMin) continue;
  
  bool skipJet = false;
  if(blacklist){
   for(unsigned int kk = 0; kk < blacklist -> size(); ++kk) {
    if(blacklist -> at(kk) == static_cast<int>(i)) skipJet = true;
   }
  }
  
  if (fabs((jets.at(i).Eta() - etaMean) / dEta) > zeppMax) continue;
  
  if (skipJet) continue;
  nJet++;  
 } 
 
 return nJet;
 
}


//  ------------------------------------------------------------



double SelectJets(std::vector<int>& it, std::vector<ROOT::Math::XYZTVector>& jets,
                  const std::string& method,
                  const double& etMin,
                  const std::vector<int>* blacklist)
{
  // initialize vector with result
  it.clear();
  it.push_back(-1);
  it.push_back(-1);
  
  
  
  // initialize the selection variable
  double maxDeta = -999999.;
  double tempDeta = 0.;
  
  double maxMJJ = -999999.;
  double tempMJJ = 0.;
  
  double maxPt = -999999.;
  double tempPt = 0.;
  
  double maxSumPt = -999999.;
  double tempSumPt = 0.;
  
  
  
  // loop over 1st jet
  for(unsigned int i = 0; i < jets.size(); ++i)
  {
    if(jets.at(i).Et() < etMin) continue;
    
    bool skipJet1 = false;
    if(blacklist)
      for(unsigned int kk = 0; kk < blacklist -> size(); ++kk)
        if(blacklist -> at(kk) == static_cast<int>(i)) skipJet1 = true;
    if(skipJet1) continue;
    
    
    
    // loop over 2nd jet
    for(unsigned int j = i+1; j < jets.size(); ++j)
    {
      if(jets.at(j).Et() < etMin) continue;
      
      bool skipJet2 = false;
      if(blacklist)
        for(unsigned int kk = 0; kk < blacklist -> size(); ++kk)
          if(blacklist -> at(kk) == static_cast<int>(j)) skipJet2 = true;
      if(skipJet2) continue;
      
      
      
      // -------------------------------------
      // select jets with different techniques
      // -------------------------------------
      
      if(method == "maxDeta")
      {
        tempDeta = deltaEta(jets.at(i).Eta(), jets.at(j).Eta());
        if(tempDeta > maxDeta)
        {
          maxDeta = tempDeta;
          
    	    it.at(0) = i;
	        it.at(1) = j;
        }
      }
      
      // -------------------------------------
      
      if(method == "maxMJJ")
      {
        tempMJJ = (jets.at(i) + jets.at(j)).mass();
        if(tempMJJ > maxMJJ)
        {
          maxMJJ = tempMJJ;
          
    	    it.at(0) = i;
	        it.at(1) = j;
        }
      }
      
      // -------------------------------------
      else if(method == "maxPt")
      {
        tempPt = sqrt( (jets.at(i) + jets.at(j)).perp2() );
        if(tempPt > maxPt)
        {
          maxPt = tempPt;
          
	        it.at(0) = i;
	        it.at(1) = j;
        }
      }
      
      // -------------------------------------
      
      else if(method == "maxSumPt")
      {
        tempSumPt = sqrt(jets.at(i).perp2()) + sqrt(jets.at(j).perp2());
        if(tempSumPt > maxSumPt)
        {
          maxSumPt = tempSumPt;
          
   	      it.at(0) = i;
	        it.at(1) = j;
        }
      }
      
      // -------------------------------------
      
      
      
    } // loop over 2nd jet
  } // loop over 1st jet
  
  
  
  if(method == "maxMJJ")
    return maxMJJ;
  
  else if(method == "maxPt")
    return maxPt;
  
  else if(method == "maxSumPt")
    return maxSumPt;
  
  else return -1.;
}

//  ------------------------------------------------------------

int SelectLepton(std::vector<ROOT::Math::XYZTVector>& leptons,
                 const std::string& method,
                 const double& ptMin,
                 const std::vector<int>* blacklist)
{
  // initialize variable with result
  int it = -1;
  
  
  
  // initialize the selection variable
  double maxPt = -999999.;
  double tempPt = 0.;
  
  
  
  // loop over leptons
  for(unsigned int i = 0; i < leptons.size(); ++i)
  {
    if( sqrt(leptons.at(i).Perp2()) < ptMin ) continue;
    
    bool skipLep = false;
    if(blacklist)
      for(unsigned int kk = 0; kk < blacklist -> size(); ++kk)
        if(blacklist -> at(kk) == static_cast<int>(i)) skipLep = true;
    if(skipLep) continue;
    
    
    
    // -------------------------------------
    // select jets with different techniques
    // -------------------------------------
    
    if(method == "maxPt")
    {
      tempPt = sqrt(leptons.at(i).perp2());
      if(tempPt > maxPt)
      {
        maxPt = tempPt;
        
        it = i;
      }
    }
    
    // -------------------------------------
    
    
    
  } // loop over leptons
  
  
  
  if(method == "maxPt")
    return it;
  
  else return -1;
}

//  ------------------------------------------------------------






int Build4JetCombinations(std::vector<std::vector<int> >& combinations, const int& nJets)
{
  combinations.clear();
  
  std::vector<int> vi;
  for(int i = 0; i < nJets; ++i)
 	  vi.push_back(i);
  
  std::vector<int> buffer;
  buffer.push_back(0);
  buffer.push_back(1);
  buffer.push_back(2);
  buffer.push_back(3);
  
  combinations.push_back(buffer);
  

  std::vector<int> oldCombination = buffer;
  while( next_permutation(vi.begin(), vi.end()) )      
  {
    if( (vi.at(0) < vi.at(1)) && (vi.at(2) < vi.at(3)) )
    {
      buffer.at(0) = vi.at(0);
      buffer.at(1) = vi.at(1);
      buffer.at(2) = vi.at(2);
      buffer.at(3) = vi.at(3);                  
      
      if(buffer == oldCombination) continue;
      
      combinations.push_back(buffer);
      oldCombination = buffer;
    }  
  }
  
  return combinations.size();
}

//  ------------------------------------------------------------

void Print4JetCombination(const std::vector<int>& combination)
{
  std::cout << "(" << combination.at(0) << "," << combination.at(1) << ")";
  std::cout << "   ---   ";
  std::cout << "(" << combination.at(2) << "," << combination.at(3) << ")";  
  std::cout << std::endl;
}

//  ------------------------------------------------------------



int Build2JetCombinations(std::vector<std::vector<int> >& combinations, const int& nJets)
{
 combinations.clear();
 
 std::vector<int> vi;
 for(int i = 0; i < nJets; ++i)
  vi.push_back(i);
 
 std::vector<int> buffer;
 buffer.push_back(0);
 buffer.push_back(1);
 
 combinations.push_back(buffer);
 
 std::vector<int> oldCombination = buffer;
 while( next_permutation(vi.begin(), vi.end()) )      
 {
  if(vi.at(0) < vi.at(1))
  {
   buffer.at(0) = vi.at(0);
   buffer.at(1) = vi.at(1);

   if(buffer == oldCombination) continue;
   
   combinations.push_back(buffer);
   oldCombination = buffer;
  }  
 }
 
 return combinations.size();
}


//  ------------------------------------------------------------


double SelectResonance(std::vector<int>& it, std::vector<ROOT::Math::XYZTVector>& objects,
		       const double& mass,
		       const double& ptMin,
		       const std::vector<int>* blacklist){
 // initialize vector with result
 it.clear();
 it.push_back(-1);
 it.push_back(-1);
 
 double minDMass = 999999.;
 double tempDMass = 0;
 
 // loop over 1st object
 for(unsigned int i = 0; i < objects.size(); ++i){
  
  if(objects.at(i).Pt() < ptMin) continue;
  
  bool skipObj1 = false;
  if(blacklist)
   for(unsigned int kk = 0; kk < blacklist -> size(); ++kk)
    if(blacklist -> at(kk) == static_cast<int>(i)) skipObj1 = true;
    if(skipObj1) continue;
    
    // loop over 2nd object
    for(unsigned int j = i+1; j < objects.size(); ++j)
    {
     if(objects.at(j).Pt() < ptMin) continue;
     
     bool skipObj2 = false;
     if(blacklist)
      for(unsigned int kk = 0; kk < blacklist -> size(); ++kk)
       if(blacklist -> at(kk) == static_cast<int>(j)) skipObj2 = true;
       if(skipObj2) continue;
       
       tempDMass = fabs((objects.at(i) + objects.at(j)).mass() - mass);
      if(tempDMass < minDMass)
      {
       minDMass = tempDMass;
       it.at(0) = i;
       it.at(1) = j;
      }
    } // loop over 2nd object
 } // loop over 1st object
 return minDMass;
}
//  ------------------------------------------------------------


double SelectResonanceOppositeCharge(std::vector<int>& it,
		       std::vector<ROOT::Math::XYZTVector>& objects,
		       std::vector<float>& charge,	     
		       const double& mass,
		       const double& ptMin,
		       const std::vector<int>* blacklist){
 // initialize vector with result
 it.clear();
 it.push_back(-1);
 it.push_back(-1);
 
 double minDMass = 999999.;
 double tempDMass = 0;
 
 // loop over 1st object
 for(unsigned int i = 0; i < objects.size(); ++i){
  if(objects.at(i).Pt() < ptMin) continue;
  double charge1 = charge.at(i);
  
  bool skipObj1 = false;
  if(blacklist)
   for(unsigned int kk = 0; kk < blacklist -> size(); ++kk)
    if(blacklist -> at(kk) == static_cast<int>(i)) skipObj1 = true;
    if(skipObj1) continue;
    
    // loop over 2nd object
    for(unsigned int j = i+1; j < objects.size(); ++j)
    {
     if(objects.at(j).Pt() < ptMin) continue;
     double charge2 = charge.at(j);     
     if(charge1 * charge2 > 0) continue;
     
     bool skipObj2 = false;
     
     if(blacklist)
      for(unsigned int kk = 0; kk < blacklist -> size(); ++kk)
       if(blacklist -> at(kk) == static_cast<int>(j)) skipObj2 = true;
       if(skipObj2) continue;
       
       tempDMass = fabs((objects.at(i) + objects.at(j)).mass() - mass);
      if(tempDMass < minDMass)
      {
       minDMass = tempDMass;
       it.at(0) = i;
       it.at(1) = j;
      }
    } // loop over 2nd object
 } // loop over 1st object
 return minDMass;
}
//  ------------------------------------------------------------



// -------------------------------------------------------------


TH1D * smartProfileX (TH2F * strip, double width){
  TProfile * stripProfile = strip->ProfileX () ;

  // (from FitSlices of TH2.h)

  double xmin = stripProfile->GetXaxis ()->GetXmin () ;
  double xmax = stripProfile->GetXaxis ()->GetXmax () ;
  int profileBins = stripProfile->GetNbinsX () ;

  std::string name = strip->GetName () ;
  name += "_smart_X" ; 
  TH1D * prof = new TH1D(name.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;
   
  int cut = 0 ; // minimum number of entries per fitted bin
  int nbins = strip->GetXaxis ()->GetNbins () ;
  int binmin = 1 ;
  int ngroup = 1 ; // bins per step
  int binmax = nbins ;

  // loop over the strip bins
  for (int bin=binmin ; bin<=binmax ; bin += ngroup) 
    {
      TH1D *hpy = strip->ProjectionY ("_temp",bin,bin+ngroup-1,"e") ;
      if (hpy == 0) continue ;
      int nentries = Int_t (hpy->GetEntries ()) ;
      if (nentries == 0 || nentries < cut) {delete hpy ; continue ;} 
 
      Int_t biny = bin + ngroup/2 ;
      
      hpy->GetXaxis ()->SetRangeUser ( hpy->GetMean () - width * hpy->GetRMS (), hpy->GetMean () + width * hpy->GetRMS ()) ;         
      prof->Fill (strip->GetXaxis ()->GetBinCenter (biny), hpy->GetMean ()) ;       
      prof->SetBinError (biny,hpy->GetRMS()) ;
      
      delete hpy ;
    } // loop over the bins

  delete stripProfile ;
  return prof ;
}


// -------------------------------------------------------------

TH1D * smartGausProfileX (TH2F * strip, double width){
  TProfile * stripProfile = strip->ProfileX () ;

  // (from FitSlices of TH2.h)

  double xmin = stripProfile->GetXaxis ()->GetXmin () ;
  double xmax = stripProfile->GetXaxis ()->GetXmax () ;
  int profileBins = stripProfile->GetNbinsX () ;

  std::string name = strip->GetName () ;
  name += "_smartGaus_X" ; 
  TH1D * prof = new TH1D(name.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;
   
  int cut = 0 ; // minimum number of entries per fitted bin
  int nbins = strip->GetXaxis ()->GetNbins () ;
  int binmin = 1 ;
  int ngroup = 1 ; // bins per step
  int binmax = nbins ;

  // loop over the strip bins
  for (int bin=binmin ; bin<=binmax ; bin += ngroup) 
    {
      TH1D *hpy = strip->ProjectionY ("_temp",bin,bin+ngroup-1,"e") ;
      if (hpy == 0) continue ;
      int nentries = Int_t (hpy->GetEntries ()) ;
      if (nentries == 0 || nentries < cut) {delete hpy ; continue ;} 
 
      Int_t biny = bin + ngroup/2 ;

      TF1 * gaussian = new TF1 ("gaussian","gaus", hpy->GetMean () - width * hpy->GetRMS (), hpy->GetMean () + width * hpy->GetRMS ()) ; 
      gaussian->SetParameter (1,hpy->GetMean ()) ;
      gaussian->SetParameter (2,hpy->GetRMS ()) ;
      hpy->Fit ("gaussian","RQL") ;           

//       hpy->GetXaxis ()->SetRangeUser ( hpy->GetMean () - width * hpy->GetRMS (), hpy->GetMean () + width * hpy->GetRMS ()) ;         
      prof->Fill (strip->GetXaxis ()->GetBinCenter (biny), gaussian->GetParameter (1)) ;       
      prof->SetBinError (biny,gaussian->GetParameter (2)) ;
      
      delete gaussian ;
      delete hpy ;
    } // loop over the bins

  delete stripProfile ;
  return prof ;
}


// -------------------------------------------------------------


TH1D * smartProfileY (TH2F * strip, double width){
 TProfile * stripProfile = strip->ProfileY () ;
 
 // (from FitSlices of TH2.h)
 
 double xmin = stripProfile->GetXaxis ()->GetXmin () ;
 double xmax = stripProfile->GetXaxis ()->GetXmax () ;
 int profileBins = stripProfile->GetNbinsX () ;
 
 std::string name = strip->GetName () ;
 name += "_smart_Y" ; 
 TH1D * prof = new TH1D(name.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;
 
 int cut = 0 ; // minimum number of entries per fitted bin
 int nbins = strip->GetYaxis ()->GetNbins () ;
 int binmin = 1 ;
 int ngroup = 1 ; // bins per step
 int binmax = nbins ;
 
 // loop over the strip bins
 for (int bin=binmin ; bin<=binmax ; bin += ngroup) 
 {
  TH1D *hpx = strip->ProjectionX ("_temp",bin,bin+ngroup-1,"e") ;
  if (hpx == 0) continue ;
  int nentries = Int_t (hpx->GetEntries ()) ;
  if (nentries == 0 || nentries < cut) {delete hpx ; continue ;} 
  
  Int_t biny = bin + ngroup/2 ;
  
  hpx->GetXaxis ()->SetRangeUser ( hpx->GetMean () - width * hpx->GetRMS (), hpx->GetMean () + width * hpx->GetRMS ()) ;         
  prof->Fill (strip->GetYaxis ()->GetBinCenter (biny), hpx->GetMean ()) ;       
  prof->SetBinError (biny,hpx->GetRMS()) ;
  
  delete hpx ;
 } // loop over the bins
 
 delete stripProfile ;
 return prof ;
}


// -------------------------------------------------------------

TH1D * smartGausProfileY (TH2F * strip, double width){
 TProfile * stripProfile = strip->ProfileY () ;
 
 // (from FitSlices of TH2.h)
 
 double xmin = stripProfile->GetXaxis ()->GetXmin () ;
 double xmax = stripProfile->GetXaxis ()->GetXmax () ;
 int profileBins = stripProfile->GetNbinsX () ;
 
 std::string name = strip->GetName () ;
 name += "_smartGaus_Y" ; 
 TH1D * prof = new TH1D(name.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;
 
 int cut = 0 ; // minimum number of entries per fitted bin
 int nbins = strip->GetYaxis ()->GetNbins () ;
 int binmin = 1 ;
 int ngroup = 1 ; // bins per step
 int binmax = nbins ;
 
 // loop over the strip bins
 for (int bin=binmin ; bin<=binmax ; bin += ngroup) 
 {
  TH1D *hpx = strip->ProjectionX ("_temp",bin,bin+ngroup-1,"e") ;
  if (hpx == 0) continue ;
  int nentries = Int_t (hpx->GetEntries ()) ;
  if (nentries == 0 || nentries < cut) {delete hpx ; continue ;} 
  
  Int_t biny = bin + ngroup/2 ;
  
  TF1 * gaussian = new TF1 ("gaussian","gaus", hpx->GetMean () - width * hpx->GetRMS (), hpx->GetMean () + width * hpx->GetRMS ()) ; 
  gaussian->SetParameter (1,hpx->GetMean ()) ;
  gaussian->SetParameter (2,hpx->GetRMS ()) ;
  hpx->Fit ("gaussian","RQL") ;           
  
  //       hpy->GetXaxis ()->SetRangeUser ( hpy->GetMean () - width * hpy->GetRMS (), hpy->GetMean () + width * hpy->GetRMS ()) ;         
  prof->Fill (strip->GetYaxis ()->GetBinCenter (biny), gaussian->GetParameter (1)) ;       
  prof->SetBinError (biny,gaussian->GetParameter (2)) ;
  
  delete gaussian ;
  delete hpx ;
 } // loop over the bins
 
 delete stripProfile ;
 return prof ;
}

// -------------------------------------------------------------


std::vector<double> getSigmaBands_FeldmanCousins (const TH1 & histo)
{
 ///==== AM Get Neyman Intervals, MPV, +/- 1 sigma ====
 
 int StepY = 1000;
 
 std::vector<double> result (5, 0.) ;
 
 double maxY = histo.GetMaximum();
 int totBins = histo.GetNbinsX();
 int iBinCenter = histo.GetMaximumBin();
 int iBinMin = iBinCenter;
 int iBinMax = iBinMin;
 int iBinMin_cycle = iBinCenter;
 int iBinMax_cycle = iBinMin_cycle;
 double totalEntries = histo.GetEffectiveEntries();
 double integral = histo.GetBinContent(iBinCenter);
 double area = integral / totalEntries;
 double DeltaY = maxY / StepY;
 double PositionY = maxY;
 
 double ValueSx = maxY;
 double ValueDx = maxY;
 
 bool doContinue68 = true; 
 int increasePerformance = 0;
 
 if (totalEntries != 0) { ///==== if 0 entries return 0,0,0,0,0
  while (doContinue68 && PositionY>0){ /// && PositionY>0 test for infinite-loop stop
   PositionY = PositionY - DeltaY; 
   ///==== look left ====
   for (int iBinSx = (iBinCenter-iBinMin_cycle); iBinSx < iBinCenter; iBinSx++){
    ValueSx = histo.GetBinContent(iBinCenter-iBinSx);
    if (ValueSx <= PositionY) {
     iBinMin = iBinCenter-iBinSx;
     break;
    }
    else {
     if (iBinSx == (iBinCenter-1)){
      iBinMin = iBinSx;
     }
    }
   }
   ///==== look right ====
   for (int iBinDx = (iBinMax_cycle-iBinCenter); iBinDx <= (totBins - iBinCenter); iBinDx++){
    ValueDx = histo.GetBinContent(iBinCenter+iBinDx);
    if (ValueDx <= PositionY) {
     iBinMax = iBinCenter+iBinDx;
     break;
    }
    else {
     if (iBinDx == ((totBins - iBinCenter)-1)){
      iBinMax = iBinDx;
     }
    }
   }
   
   integral = histo.Integral(iBinMin,iBinMax);
   area = integral / totalEntries;
//    std::cerr << " area 68[" << iBinMin << ":" << iBinCenter << ":" << iBinMax << "] [" << PositionY << ":" << maxY << "] = " << area << " = " << integral << " / " << totalEntries << std::endl;
   if (area > 0.68) {
    if (increasePerformance == 0){
     increasePerformance = 1;
     PositionY = PositionY + DeltaY;
     DeltaY  = DeltaY / 10.;
    }
    else {
     doContinue68 = false;
    }
   }
   else {
    iBinMin_cycle = iBinMin;
    iBinMax_cycle = iBinMax;
   }
  }
  
  ///=== 68% and mean ===
  result.at(1) = histo.GetBinCenter(iBinMin);
  result.at(2) = histo.GetBinCenter(iBinCenter);
  result.at(3) = histo.GetBinCenter(iBinMax);
  ///====================
  
  
  bool doContinue95 = true; 
  increasePerformance = 0;
  
  while (doContinue95 && PositionY>0){ /// && PositionY>0 test for infinite-loop stop
   PositionY = PositionY - DeltaY;
   ///==== look left ====
   for (int iBinSx = (iBinCenter-iBinMin_cycle); iBinSx < iBinCenter; iBinSx++){
    ValueSx = histo.GetBinContent(iBinCenter-iBinSx);
//     std::cerr << " ValueSx[" << iBinSx << "] = " << ValueSx << " " ;
    if (ValueSx <= PositionY) {
//      std::cerr << "found SX " << std::endl;
     iBinMin = iBinCenter-iBinSx;
     break;
    }
    else {
     if (iBinSx == (iBinCenter-1)){
      iBinMin = iBinSx;
     }
    }
   }
   ///==== look right ====
   for (int iBinDx = (iBinMax_cycle-iBinCenter); iBinDx <= (totBins - iBinCenter); iBinDx++){
    ValueDx = histo.GetBinContent(iBinCenter+iBinDx);
//     std::cerr << " ValueDx[" << iBinDx << "] = " << ValueDx << " " ;
    if (ValueDx <= PositionY) {
//      std::cerr << "found DX " << std::endl;
     iBinMax = iBinCenter+iBinDx;
     break;
    }
    else {
     if (iBinDx == ((totBins - iBinCenter)-1)){
      iBinMax = iBinDx;
     }
    }
   }
   
   integral = histo.Integral(iBinMin,iBinMax);
   area = integral / totalEntries;
//    std::cerr << " area 95[" << iBinMin << ":" << iBinCenter << ":" << iBinMax << "] [" << PositionY << ":" << maxY << "] = " << area << " = " << integral << " / " << totalEntries << "    SX : " << iBinCenter-iBinMin_cycle << " : " << iBinCenter << " DX : " << (iBinMax_cycle-iBinCenter) << " : " << (totBins - iBinCenter) << std::endl;
   if (area > 0.95) {
    if (increasePerformance == 0){
     increasePerformance = 1;
     PositionY = PositionY + DeltaY;
     DeltaY  = DeltaY / 10.;
    }
    else {
     doContinue95 = false;
    }
   }
   else {
    iBinMin_cycle = iBinMin;
    iBinMax_cycle = iBinMax;
   }
  }
  
  ///=== 95% ===
  result.at(0) = histo.GetBinCenter(iBinMin);
  result.at(4) = histo.GetBinCenter(iBinMax);
  ///===========
 }
 
 return result ;
}



// -------------------------------------------------------------


std::vector<TH1D*> smartGausProfileY_FeldmanCousins (TH2F * strip, double width){
 TProfile * stripProfile = strip->ProfileY () ;
 
 double xmin = stripProfile->GetXaxis ()->GetXmin () ;
 double xmax = stripProfile->GetXaxis ()->GetXmax () ;
 int profileBins = stripProfile->GetNbinsX () ;
 
 std::string nameLow = strip->GetName () ;
 nameLow += "_smartGaus_FeldmanCousins_Low_Y" ; 
 TH1D * prof_Low = new TH1D(nameLow.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;

 std::string nameMean = strip->GetName () ;
 nameMean += "_smartGaus_FeldmanCousins_Mean_Y" ; 
 TH1D * prof_Mean = new TH1D(nameMean.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;
 
 std::string nameHigh = strip->GetName () ;
 nameHigh += "_smartGaus_FeldmanCousins_High_Y" ; 
 TH1D * prof_High = new TH1D(nameHigh.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;
 
 int cut = 0 ; // minimum number of entries per fitted bin
 int nbins = strip->GetYaxis ()->GetNbins () ;
 int binmin = 1 ;
 int ngroup = 1 ; // bins per step
 int binmax = nbins ;
 
 // loop over the strip bins
 for (int bin=binmin ; bin<=binmax ; bin += ngroup) 
 {
  TH1D *hpx = strip->ProjectionX ("_temp",bin,bin+ngroup-1,"e") ;
  if (hpx == 0) continue ;
  int nentries = Int_t (hpx->GetEntries ()) ;
  if (nentries == 0 || nentries < cut) {delete hpx ; continue ;} 
  
  Int_t biny = bin + ngroup/2 ;
  
  std::vector<double> band = getSigmaBands_FeldmanCousins(*hpx);

  prof_Low->Fill (strip->GetYaxis ()->GetBinCenter (biny), band.at(1)) ;       
  prof_Mean->Fill (strip->GetYaxis ()->GetBinCenter (biny), band.at(2)) ;       
  prof_High->Fill (strip->GetYaxis ()->GetBinCenter (biny), band.at(3)) ;       
  
  delete hpx ;
 } // loop over the bins
 
 delete stripProfile ;
 
 std::vector<TH1D*> result;
 result.push_back(prof_Low);
 result.push_back(prof_Mean);
 result.push_back(prof_High);
 
 return result ;
}


// -------------------------------------------------------------



std::vector<TH1D*> smartProfileX_FeldmanCousins (TH2F * strip, double width){
 TProfile * stripProfile = strip->ProfileX () ;
 
 // (from FitSlices of TH2.h)
 
 double xmin = stripProfile->GetXaxis ()->GetXmin () ;
 double xmax = stripProfile->GetXaxis ()->GetXmax () ;
 int profileBins = stripProfile->GetNbinsX () ;
 
 std::string nameLow = strip->GetName () ;
 nameLow += "_smartGaus_FeldmanCousins_Low_X" ; 
 TH1D * prof_Low = new TH1D(nameLow.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;
 
 std::string nameMean = strip->GetName () ;
 nameMean += "_smartGaus_FeldmanCousins_Mean_X" ; 
 TH1D * prof_Mean = new TH1D(nameMean.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;
 
 std::string nameHigh = strip->GetName () ;
 nameHigh += "_smartGaus_FeldmanCousins_High_X" ; 
 TH1D * prof_High = new TH1D(nameHigh.c_str (),strip->GetTitle (),profileBins,xmin,xmax) ;
 
 int cut = 0 ; // minimum number of entries per fitted bin
 int nbins = strip->GetXaxis ()->GetNbins () ;
 int binmin = 1 ;
 int ngroup = 1 ; // bins per step
 int binmax = nbins ;
 
 // loop over the strip bins
 for (int bin=binmin ; bin<=binmax ; bin += ngroup) 
 {
  TH1D *hpy = strip->ProjectionY ("_temp",bin,bin+ngroup-1,"e") ;
  if (hpy == 0) continue ;
  int nentries = Int_t (hpy->GetEntries ()) ;
  if (nentries == 0 || nentries < cut) {delete hpy ; continue ;} 
  
  Int_t biny = bin + ngroup/2 ;
  
  std::vector<double> band = getSigmaBands_FeldmanCousins(*hpy);
  
  prof_Low->Fill (strip->GetXaxis ()->GetBinCenter (biny), band.at(1)) ;       
  prof_Mean->Fill (strip->GetXaxis ()->GetBinCenter (biny), band.at(2)) ;       
  prof_High->Fill (strip->GetXaxis ()->GetBinCenter (biny), band.at(3)) ;       
  
  delete hpy ;
 } // loop over the bins
 
 delete stripProfile ;

 std::vector<TH1D*> result;
 result.push_back(prof_Low);
 result.push_back(prof_Mean);
 result.push_back(prof_High);
 
 return result ;
}


// -------------------------------------------------------------



/**
get the points on the definition set of a TH1F such that the region contains 68% of the area
and the tails contain the same 15.8%
*/
std::pair<double, double> getLimit_sameTails::operator() (const TH1D & histo) 
{
 double threshold = histo.GetEntries () * 0.158655 ;
 
 int entries1s_low = 0;
 double banda1s_low = 0.;
 
 // left tail
 for (int bandaBin = 1 ; bandaBin <= histo.GetNbinsX () ; ++bandaBin)
 {
  entries1s_low += histo.GetBinContent (bandaBin) ;
  if ( entries1s_low >= threshold)
  {
   banda1s_low = histo.GetBinCenter (bandaBin) ;
   break ;
  }
 } // left tail
 
 int entries1s_high = 0;
 double banda1s_high = 0.;
 
 // right tail
 for (int bandaBin = histo.GetNbinsX () ; bandaBin > 0 ; --bandaBin)
 {
  entries1s_high += histo.GetBinContent (bandaBin) ;
  if ( entries1s_high >= threshold)
  {
   banda1s_high = histo.GetBinCenter (bandaBin) ;
   break ;
  }
 } // right tail
 return std::pair<double, double> (banda1s_low, banda1s_high) ;  
}


// -------------------------------------------------------------


/**
get the points on the definition set of a TH1F such that the region contains 68% of the area
and the size of the region is minimal (Neyman intervals - Feldman Cousins)
*/

std::pair<double, double> getLimit_FC::operator() (const TH1D & histo) 
{
 std::vector<double> band = getSigmaBands_FeldmanCousins (histo);
 return std::pair<double, double> (band.at(1), band.at(3)) ;  
}


// -------------------------------------------------------------

/**
build a TGraphError starting from a TH1F
*/
TGraphErrors buildGEfromH (const TH1D & histo) 
{
 TVectorF xV(histo.GetNbinsX());
 TVectorF yV(histo.GetNbinsX());
 TVectorF errxV(histo.GetNbinsX());
 TVectorF erryV(histo.GetNbinsX());
 for (int iBin = 0; iBin<histo.GetNbinsX(); iBin++){
  xV[iBin] = histo.GetBinCenter(iBin);
  yV[iBin] = histo.GetBinContent(iBin);
  errxV[iBin] = histo.GetBinWidth(iBin);
  erryV[iBin] = histo.GetBinError(iBin);
 }  
 TGraphErrors g (xV, yV, errxV, erryV) ;
 return g ;
}



// -------------------------------------------------------------



