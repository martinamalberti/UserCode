import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
    )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V9::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )

process.source = cms.Source("PoolSource",
       #fileNames = cms.untracked.vstring('file:/tmp/malberti/pythia_Wenu_1000_1100_cfi_py_RAW2DIGI_RECO_recoSW217_rawSW217.root')
       fileNames = cms.untracked.vstring('file:/tmp/malberti/pythia_Wenu_200_250_cfi_py_RAW2DIGI_RECO_10.root')
)



process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")
process.load("WprimeAnalysis.WprimeENUAnalysis.WprimePATConfig_cfi")

#Analysis
process.myanalysis = cms.EDAnalyzer('WprimeAnalyzerPAT',
                                   electronTag = cms.InputTag("allLayer1Electrons"),
                                   jetTag      = cms.InputTag("allLayer1Jets"),
                                   metTag      = cms.InputTag("allLayer1METs"),
                                   electronID = cms.untracked.string("eidRobustHighEnergy"),
                                   TriggerResults = cms.InputTag("TriggerResults","","HLT"),
                                   )


process.TFileService = cms.Service("TFileService",
    #fileName = cms.string("testRelValPt1000.root")
    fileName = cms.string("testWenu1000.root")
)

## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)


#process.p = cms.Path( process.patLayer0 * process.patLayer1 * process.myanalysis)
process.p = cms.Path( process.wprimePATSequence * process.myanalysis) 



## process.load("Configuration.EventContent.EventContent_cff")
## process.out = cms.OutputModule("PoolOutputModule",
##             process.FEVTSIMEventContent,
##             fileName = cms.untracked.string('file:/tmp/malberti/testPat.root')
## )

## process.out.outputCommands.append('drop *_*_*_*')
## process.out.outputCommands.append('keep *_*_*_myprocess')

## process.outpath = cms.EndPath(process.wprimePATSequence *  process.out)

