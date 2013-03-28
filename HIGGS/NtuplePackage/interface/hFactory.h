#ifndef hFactory_h
#define hFactory_h

#include <map>
#include <string>
#include <functional>
#include "hChain.h"
#include "h2Chain.h"

class hFactory
{
  public :
    hFactory (const std::string& fileName = "no", bool print = true) ;
    ~hFactory () ;

    void add_h1 (const std::string& baseName, const std::string& baseTitle, 
                 int nbins, double min, double max, int NUM) ;

    template <class UnaryFunction>
    void applyToAll (UnaryFunction function) 
      {
        for (std::map <std::string, hChain *>::iterator mapIt = m_H1content.begin () ;
             mapIt != m_H1content.end () ;
             ++mapIt)
          {
            function (mapIt->second) ;
          }
      }

    void Fill (const std::string & name, int i, double val) ;
    hChain * operator[] (const std::string& name) ;
    void Print (int isLog = 0, int rebin = 1) ;

  private :
  
    std::string m_fileName ;
    bool m_print ;
    std::map <std::string, hChain *> m_H1content ;

} ;

#endif
