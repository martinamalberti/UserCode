/*                                                                                                                                                                                 
source ~/Desktop/setupRoot5.25.sh                                                                                                                                                  
rootcint -f dict.cpp -c LinkDef.h                                                                                                                                                  
c++ -o testReader `root-config --cflags --glibs` -lGenVector testReader.cpp dict.cpp treeReader.cc                                                                                 
c++ -Wall -o analysis `root-config --cflags --glibs` -lGenVector analysis.cpp dict.cpp treeReader.cc hFactory.cc hChain.cc stdHisto.cc
*/

#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"






int main(int argc, char** argv)
{
  //Check if all nedeed arguments to parse are there                                                                                                                               
  if(argc != 2)
  {
    std::cerr << ">>>>> analysis.cpp::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
  
  
  
  // Parse the config file                                                                                                                                                          
  parseConfigFile (argv[1]) ;
  
  std::string treeName      = gConfigParser -> readStringOption("Input::treeName");
  std::string inputFileList = gConfigParser -> readStringOption("Input::inputFileList");
  
  int entryMAX = gConfigParser -> readIntOption("Input::entryMAX");
  int entryMODULO = gConfigParser -> readIntOption("Input::entryMODULO");
  
  
  
  // Open ntple
  TChain* chain = new TChain(treeName.c_str());
  if(!FillChain(*chain, inputFileList.c_str())) return 1;
  treeReader reader((TTree*)(chain));
  
    
  
  for(int entry = 0 ; entry < chain -> GetEntries() ; ++entry)
  {
    reader.GetEntry(entry);
    if((entry%entryMODULO) == 0) std::cout << ">>>>> analysis::GetEntry " << entry << std::endl;   
    if(entry == entryMAX) break;
    
  } //PG loop over the events
  
  
  
  return 0;
}
