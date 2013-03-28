#include "stdHisto.h"



//!ctor
stdHisto::stdHisto(const int& nStep,
                   const std::string& outFileName,
                   treeReader* reader):
 m_nStep(nStep),
 m_reader(reader),
 m_outFileName(outFileName),
 m_hFactory(new hFactory(m_outFileName, true))
{}

// ------------------------------------------------






//!dtor
stdHisto::~stdHisto()
{
  delete m_hFactory;
}

// ------------------------------------------------






//! add histograms
void stdHisto::Add1(const std::string& histoName,
                    const int& nStep)
{
  m_hFactory -> add_h1(histoName+"_n",      histoName+"_n",       100,   0.,  100., nStep);
  m_hFactory -> add_h1(histoName+"_energy", histoName+"_energy", 3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_et",     histoName+"_et",     3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_p",      histoName+"_p",      3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_pt",     histoName+"_pt",     3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_pl",     histoName+"_pl",     3000,   0., 3000., nStep);  
  m_hFactory -> add_h1(histoName+"_eta",    histoName+"_eta",     400, -10.,   10., nStep);
  m_hFactory -> add_h1(histoName+"_absEta", histoName+"_absEta",  200,   0.,   10., nStep);
  m_hFactory -> add_h1(histoName+"_phi",    histoName+"_phi",     200,   0.,    5., nStep);
}

// ------------------------------------------------






void stdHisto::Add2(const std::string& histoName,
                    const int& nStep)
{
  m_hFactory -> add_h1(histoName+"_mass",   histoName+"_energy", 3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_energy", histoName+"_energy", 3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_et",     histoName+"_et",     3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_p",      histoName+"_p",      3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_pt",     histoName+"_pt",     3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_pl",     histoName+"_pl",     3000,   0., 3000., nStep);  
  m_hFactory -> add_h1(histoName+"_eta",    histoName+"_eta",     400, -10.,   10., nStep);
  m_hFactory -> add_h1(histoName+"_absEta", histoName+"_absEta",  200,   0.,   10., nStep);
  m_hFactory -> add_h1(histoName+"_phi",    histoName+"_phi",     200,   0.,    5., nStep);
  m_hFactory -> add_h1(histoName+"_Deta",   histoName+"_Deta",    200,   0.,   10., nStep);
  m_hFactory -> add_h1(histoName+"_Dphi",   histoName+"_Dphi",    200,   0.,    5., nStep);
  m_hFactory -> add_h1(histoName+"_DR",     histoName+"_DR",      200,   0.,   10., nStep);
  
  m_hFactory -> add_h1(histoName+"_max_energy", histoName+"_max_energy", 3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_max_et",     histoName+"_max_et",     3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_max_p",      histoName+"_max_p",      3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_max_pt",     histoName+"_max_pt",     3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_max_pl",     histoName+"_max_pl",     3000,   0., 3000., nStep);  
  m_hFactory -> add_h1(histoName+"_max_eta",    histoName+"_max_eta",     400, -10.,   10., nStep);
  m_hFactory -> add_h1(histoName+"_max_absEta", histoName+"_max_absEta",  200,   0.,   10., nStep);
  
  m_hFactory -> add_h1(histoName+"_min_energy", histoName+"_min_energy", 3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_min_et",     histoName+"_min_et",     3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_min_p",      histoName+"_min_p",      3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_min_pt",     histoName+"_min_pt",     3000,   0., 3000., nStep);
  m_hFactory -> add_h1(histoName+"_min_pl",     histoName+"_min_pl",     3000,   0., 3000., nStep);  
  m_hFactory -> add_h1(histoName+"_min_eta",    histoName+"_min_eta",     400, -10.,   10., nStep);
  m_hFactory -> add_h1(histoName+"_min_absEta", histoName+"_min_absEta",  200,   0.,   10., nStep);

}

// ------------------------------------------------

void stdHisto::Add1Float(const std::string& histoName,
                    const int& nStep,
		    int nbins,
		    double min,
		    double max
		    )
{
 m_hFactory -> add_h1(histoName+"_value",      histoName+"_value",       nbins,   min,  max, nStep);
}

// ------------------------------------------------





//! fill histograms
void stdHisto::Fill1(const std::string& histoName,
                     const std::string& branchName,
                     const int& step,
                     std::vector<int>* selectionIt)
{
  int n = 0;
  for(unsigned int it = 0; it < m_reader -> Get4V(branchName) -> size(); ++it)
  {
    if(selectionIt != NULL)
      if( (selectionIt -> at(it)) != 1 )
        continue;
    ++n;
    m_hFactory -> Fill( histoName+"_energy", step, m_reader->Get4V(branchName)->at(it).E());
    m_hFactory -> Fill( histoName+"_et",     step, m_reader->Get4V(branchName)->at(it).Et());
    m_hFactory -> Fill( histoName+"_p",      step, sqrt(m_reader->Get4V(branchName)->at(it).P()));    
    m_hFactory -> Fill( histoName+"_pt",     step, m_reader->Get4V(branchName)->at(it).pt());        
    m_hFactory -> Fill( histoName+"_pl",     step, m_reader->Get4V(branchName)->at(it).pz());            
    m_hFactory -> Fill( histoName+"_eta",    step, m_reader->Get4V(branchName)->at(it).eta() );
    m_hFactory -> Fill( histoName+"_absEta", step, fabs(m_reader->Get4V(branchName)->at(it).eta()) );
    m_hFactory -> Fill( histoName+"_phi",    step, m_reader->Get4V(branchName)->at(it).phi() );
  }
  
  m_hFactory -> Fill( histoName+"_n", step, n ) ;
}

// ------------------------------------------------




//! fill histograms
void stdHisto::Fill1(const std::vector<ROOT::Math::XYZTVector>& vet,
                     const std::string& histoName,
                     const int& step)
{
  for(unsigned int it = 0; it < vet.size(); ++it)
  {
    m_hFactory -> Fill( histoName+"_energy", step, (vet.at(it)).E());
    m_hFactory -> Fill( histoName+"_et",     step, (vet.at(it)).Et());
    m_hFactory -> Fill( histoName+"_p",      step, (vet.at(it)).P());
    m_hFactory -> Fill( histoName+"_pt",     step, (vet.at(it)).pt());
    m_hFactory -> Fill( histoName+"_pl",     step, (vet.at(it)).pz());
    m_hFactory -> Fill( histoName+"_eta",    step, (vet.at(it)).eta());
    m_hFactory -> Fill( histoName+"_absEta", step, (vet.at(it)).eta());
    m_hFactory -> Fill( histoName+"_phi",    step, (vet.at(it)).phi());
  }
  
  m_hFactory -> Fill( histoName+"_n", step, vet.size() );
}

// ------------------------------------------------



//! fill histograms
void stdHisto::Fill1(const ROOT::Math::XYZTVector& p,
                     const std::string& histoName,
                     const int& step)
{
  m_hFactory -> Fill( histoName+"_energy", step, p.E());
  m_hFactory -> Fill( histoName+"_et",     step, p.Et());
  m_hFactory -> Fill( histoName+"_p",      step, p.P());
  m_hFactory -> Fill( histoName+"_pt",     step, p.pt());
  m_hFactory -> Fill( histoName+"_pl",     step, p.pz());
  m_hFactory -> Fill( histoName+"_eta",    step, p.eta());
  m_hFactory -> Fill( histoName+"_absEta", step, p.eta());
  m_hFactory -> Fill( histoName+"_phi",    step, p.phi());
  m_hFactory -> Fill( histoName+"_n", step, 1 );
}

// ------------------------------------------------






//! fill histograms
void stdHisto::Fill2(const std::string& histoName,
                     const std::string& branchName,
                     const int& step,
                     const int& it1, const int& it2)
{
  ROOT::Math::XYZTVector v1 = m_reader -> Get4V(branchName) -> at(it1);
  ROOT::Math::XYZTVector v2 = m_reader -> Get4V(branchName) -> at(it2);
  ROOT::Math::XYZTVector vSum = v1 + v2;

  m_hFactory -> Fill( histoName+"_mass",   step, vSum.mass() );
  m_hFactory -> Fill( histoName+"_energy", step, vSum.E() );
  m_hFactory -> Fill( histoName+"_et",     step, vSum.Et() );
  m_hFactory -> Fill( histoName+"_p",      step, vSum.P() );    
  m_hFactory -> Fill( histoName+"_pt",     step, vSum.pt() );        
  m_hFactory -> Fill( histoName+"_pl",     step, vSum.pz() );            
  m_hFactory -> Fill( histoName+"_eta",    step, vSum.eta() );
  m_hFactory -> Fill( histoName+"_absEta", step, fabs(vSum.eta()) );
  m_hFactory -> Fill( histoName+"_phi",    step, vSum.phi() );
  m_hFactory -> Fill( histoName+"_Deta",   step, fabs(v1.eta() - v2.eta()) );
  m_hFactory -> Fill( histoName+"_Dphi",   step, ROOT::Math::VectorUtil::DeltaPhi(v1, v2) );
  m_hFactory -> Fill( histoName+"_DR",     step, ROOT::Math::VectorUtil::DeltaR(v1, v2) );
  
  m_hFactory -> Fill( histoName+"_max_energy", step, std::max(v1.E(), v2.E()) );
  m_hFactory -> Fill( histoName+"_max_et",     step, std::max(v1.Et(), v2.Et()) );
  m_hFactory -> Fill( histoName+"_max_p",      step, std::max(v1.P(), v2.P()) );    
  m_hFactory -> Fill( histoName+"_max_pt",     step, std::max(v1.pt(), v2.pt()) );        
  m_hFactory -> Fill( histoName+"_max_pl",     step, std::max(v1.pz(), v2.pz()) );            
  m_hFactory -> Fill( histoName+"_max_eta",    step, std::max(v1.eta(), v2.eta()) );
  m_hFactory -> Fill( histoName+"_max_absEta", step, std::max(fabs(v1.eta()), fabs(v2.eta())) );

  m_hFactory -> Fill( histoName+"_min_energy", step, std::min(v1.E(), v2.E()) );
  m_hFactory -> Fill( histoName+"_min_et",     step, std::min(v1.Et(), v2.Et()) );
  m_hFactory -> Fill( histoName+"_min_p",      step, std::min(v1.P(), v2.P()) );    
  m_hFactory -> Fill( histoName+"_min_pt",     step, std::min(v1.pt(), v2.pt()) );        
  m_hFactory -> Fill( histoName+"_min_pl",     step, std::min(v1.pz(), v2.pz()) );            
  m_hFactory -> Fill( histoName+"_min_eta",    step, std::min(v1.eta(), v2.eta()) );
  m_hFactory -> Fill( histoName+"_min_absEta", step, std::min(fabs(v1.eta()), fabs(v2.eta())) );

}

// ------------------------------------------------



//! fill histograms
void stdHisto::Fill2(const ROOT::Math::XYZTVector& v1,
                     const ROOT::Math::XYZTVector& v2, 
                     const std::string& histoName,
                     const int& step)
{
  ROOT::Math::XYZTVector vSum = v1 + v2;

  m_hFactory -> Fill( histoName+"_mass",   step, vSum.mass() );
  m_hFactory -> Fill( histoName+"_energy", step, vSum.E() );
  m_hFactory -> Fill( histoName+"_et",     step, vSum.Et() );
  m_hFactory -> Fill( histoName+"_p",      step, vSum.P() );    
  m_hFactory -> Fill( histoName+"_pt",     step, vSum.pt() );        
  m_hFactory -> Fill( histoName+"_pl",     step, vSum.pz() );            
  m_hFactory -> Fill( histoName+"_eta",    step, vSum.eta() );
  m_hFactory -> Fill( histoName+"_absEta", step, vSum.eta() );
  m_hFactory -> Fill( histoName+"_phi",    step, vSum.phi() );
  m_hFactory -> Fill( histoName+"_Deta",   step, fabs(v1.eta() - v2.eta()) );
  m_hFactory -> Fill( histoName+"_Dphi",   step, ROOT::Math::VectorUtil::DeltaPhi(v1, v2) );
  m_hFactory -> Fill( histoName+"_DR",     step, ROOT::Math::VectorUtil::DeltaR(v1, v2) );
  
  m_hFactory -> Fill( histoName+"_max_energy", step, std::max(v1.E(), v2.E()) );
  m_hFactory -> Fill( histoName+"_max_et",     step, std::max(v1.Et(), v2.Et()) );
  m_hFactory -> Fill( histoName+"_max_p",      step, std::max(v1.P(), v2.P()) );    
  m_hFactory -> Fill( histoName+"_max_pt",     step, std::max(v1.pt(), v2.pt()) );        
  m_hFactory -> Fill( histoName+"_max_pl",     step, std::max(v1.pz(), v2.pz()) );            
  m_hFactory -> Fill( histoName+"_max_eta",    step, std::max(v1.eta(), v2.eta()) );
  m_hFactory -> Fill( histoName+"_max_absEta", step, std::max(fabs(v1.eta()), fabs(v2.eta())) );

  m_hFactory -> Fill( histoName+"_min_energy", step, std::min(v1.E(), v2.E()) );
  m_hFactory -> Fill( histoName+"_min_et",     step, std::min(v1.Et(), v2.Et()) );
  m_hFactory -> Fill( histoName+"_min_p",      step, std::min(v1.P(), v2.P()) );    
  m_hFactory -> Fill( histoName+"_min_pt",     step, std::min(v1.pt(), v2.pt()) );        
  m_hFactory -> Fill( histoName+"_min_pl",     step, std::min(v1.pz(), v2.pz()) );            
  m_hFactory -> Fill( histoName+"_min_eta",    step, std::min(v1.eta(), v2.eta()) );
  m_hFactory -> Fill( histoName+"_min_absEta", step, std::min(fabs(v1.eta()), fabs(v2.eta())) );

}

// ------------------------------------------------


void stdHisto::Fill1Float(const std::string& histoName,
                     const std::string& branchName,
                     const int& step,
                     std::vector<int>* selectionIt)
{
  int n = 0;
  for(unsigned int it = 0; it < m_reader -> GetFloat(branchName) -> size(); ++it)
  {
    if(selectionIt != NULL)
      if( (selectionIt -> at(it)) != 1 )
        continue;
    ++n;
    m_hFactory -> Fill( histoName+"_value", step, m_reader->GetFloat(branchName)->at(it));
  }
}

// ------------------------------------------------


void stdHisto::Fill1Float(const std::string& histoName,
                     const double& value,
                     const int& step)
{
 m_hFactory -> Fill( histoName+"_value", step, value); 
}

// ------------------------------------------------

