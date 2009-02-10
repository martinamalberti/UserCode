#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeAnalyzer.h"
#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeAnalyzerPAT.h"
#include "WprimeAnalysis/WprimeENUAnalysis/interface/ElectronDuplicateRemover.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(WprimeAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(WprimeAnalyzerPAT);
DEFINE_ANOTHER_FWK_MODULE(ElectronDuplicateRemover);
