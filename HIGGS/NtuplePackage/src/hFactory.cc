#include "hFactory.h"

#include <iostream>
#include <string>

hFactory::hFactory (const std::string& fileName, 
                    bool print) :
 m_fileName (fileName) ,
 m_print (print)                     
{}


//PG --------------------------------------------------------   


hFactory::~hFactory () 
{
//  if (m_H1content.size () && m_print)
//    for (std::map <std::string, hChain *>::iterator mapIt = m_H1content.begin () ;
//         mapIt != m_H1content.end () ;
//         ++mapIt)
//      {
//        mapIt->second->Print (0,4) ; 
//        //PG FIXME generalizzare le opzioni di stampa    
//      }
  if (m_H1content.size () && m_fileName != "no")
    {
      TFile histosFile (m_fileName.c_str (),"update") ;
      histosFile.cd () ;
      for (std::map <std::string, hChain *>::iterator mapIt = m_H1content.begin () ;
           mapIt != m_H1content.end () ;
           ++mapIt)
        {
          mapIt->second->Write (mapIt->first, histosFile) ; 
        }
      histosFile.Close () ;
    }  
}


//PG --------------------------------------------------------   


void
hFactory::add_h1 (const std::string& baseName, const std::string& baseTitle, 
                  int nbins, double min, double max, int NUM) 
{
  if (m_H1content.find (baseName) != m_H1content.end ())
    {
      std::cerr << "ERROR : histograms series " << baseName
                << " already existing, NOT replaced" << std::endl ;
      return ;                
    }
  hChain * dummy = new hChain (baseName, baseTitle, 
                               nbins, min, max, NUM) ;
  m_H1content[baseName] = dummy ;     
  return ;                          
}


//PG --------------------------------------------------------   


void 
hFactory::Fill (const std::string & name, int i, double val) 
  {
    if (m_H1content.find (name) != m_H1content.end ())
      m_H1content[name]->Fill (i,val) ;
    return ;
  }


//PG --------------------------------------------------------   


hChain * 
hFactory::operator[] (const std::string& name)
  {
    if (m_H1content.find (name) == m_H1content.end ())
      {
        std::cerr << "ERROR : histograms series " << name
                  << " not existing, returning the first one" << std::endl ;
        return m_H1content.begin ()->second ;                
      }
    return m_H1content[name] ;
  
  }


//PG --------------------------------------------------------   


void 
hFactory::Print (int isLog, int rebin)
{
  if (!m_H1content.size ()) return ;
  for (std::map <std::string, hChain *>::iterator mapIt = m_H1content.begin () ;
       mapIt != m_H1content.end () ;
       ++mapIt)
    {
      mapIt->second->Print (isLog, rebin) ; 
    }
}
