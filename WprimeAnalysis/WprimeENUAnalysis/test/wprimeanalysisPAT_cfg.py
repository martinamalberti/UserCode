import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
#process.MessageLogger.cerr.threshold = 'WARNING'
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
    input = cms.untracked.int32(10)
    )

process.source = cms.Source("PoolSource",
       #fileNames = cms.untracked.vstring('/store/relval/CMSSW_2_1_8/RelValWE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/06763DA3-8C82-DD11-9682-000423D99658.root')
       fileNames = cms.untracked.vstring('/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/0AF3FBC9-AA99-DD11-9E92-000423D98920.root',
        '/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/24C672DF-A499-DD11-8F6C-000423D98EC4.root',
        '/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/52ABBD98-9799-DD11-8D22-001617C3B73A.root',
        '/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/AC5E763A-A899-DD11-8B94-001617C3B79A.root',
        '/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/BCFFE400-9A99-DD11-A23B-000423D98B28.root',
        '/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/C21A3BD8-A099-DD11-B725-0019DB29C620.root',
        '/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/C2FD6F4A-FD99-DD11-88D7-000423D996B4.root',
        '/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/D42C1BA8-9F99-DD11-987E-001617C3B778.root',
        '/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/DE84CC0E-9999-DD11-9C87-000423D6A6F4.root',
        '/store/relval/CMSSW_2_1_10/RelValSingleElectronPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/F8B5BF99-A299-DD11-AD4A-000423D6CAF2.root')
)

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

process.myanalysis = cms.EDAnalyzer('WprimeAnalyzerPAT',
                                   electronTag = cms.InputTag("selectedLayer1Electrons"),
                                   jetTag      = cms.InputTag("selectedLayer1Jets"),
                                   metTag      = cms.InputTag("selectedLayer1METs"),
                                   electronID = cms.untracked.string("eidRobustHighEnergy"),
                                   TriggerResults = cms.InputTag("TriggerResults","","HLT"),
                                   )

# HEEP Electron ID and Isolation


#MATCHING with MC truth


process.TFileService = cms.Service("TFileService",
    fileName = cms.string("testRelValPt1000.root")
)

## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)

process.p = cms.Path( process.patLayer0 + process.patLayer1 + process.myanalysis)


## process.load("Configuration.EventContent.EventContent_cff")
## process.out = cms.OutputModule("PoolOutputModule",
##             process.FEVTSIMEventContent,
##             fileName = cms.untracked.string('file:/tmp/malberti/eleIso.root')
## )

## process.out.outputCommands.append('drop *_*_*_*')
## process.out.outputCommands.append('keep *_*_*_myprocess')

## process.outpath = cms.EndPath(process.patLayer0 + process.patLayer1 + process.out)
