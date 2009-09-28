import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
    )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
           'rfio:/castor/cern.ch/user/m/malberti/CMSSW31X/WPRIME/PAT/pat_Wenu_10.root',
           'rfio:/castor/cern.ch/user/m/malberti/CMSSW31X/WPRIME/PAT/pat_Wenu_11.root',
           'rfio:/castor/cern.ch/user/m/malberti/CMSSW31X/WPRIME/PAT/pat_Wenu_12.root'
           )
)

#Analysis
process.myanalysis = cms.EDAnalyzer('WprimeAnalyzerPAT',
                                   electronTag = cms.InputTag("cleanLayer1Electrons"),
                                   jetTag      = cms.InputTag("cleanLayer1Jets"),
                                   metTag      = cms.InputTag("layer1METs"),
                                   muonTag     = cms.InputTag("cleanLayer1Muons"),
                                   electronID  = cms.untracked.string("eidRobustHighEnergy"),
                                   btagAlgo    = cms.untracked.string("jetBProbabilityBJetTags"),
                                    TriggerResults = cms.InputTag("TriggerResults","","HLT"),
                                   )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("/tmp/malberti/treeWenu31X.root")
)




process.p = cms.Path( process.myanalysis) 



## process.load("Configuration.EventContent.EventContent_cff")
## process.out = cms.OutputModule("PoolOutputModule",
##             process.FEVTSIMEventContent,
##             fileName = cms.untracked.string('file:/tmp/malberti/testPat.root')
## )

## process.out.outputCommands.append('drop *_*_*_*')
## process.out.outputCommands.append('keep *_*_*_myprocess')
## process.out.outputCommands.append('keep *_*_*_HLT')

## process.outpath = cms.EndPath(process.wprimePATSequence *  process.out)

