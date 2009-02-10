import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.load("Configuration.StandardSequences.Geometry_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
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

# HEEP Electron ID and Isolation
process.load("WprimeAnalysis.WprimeENUAnalysis.electronCleaning_cff")
from WprimeAnalysis.WprimeENUAnalysis.electronCleaning_cff import *

process.load("WprimeAnalysis.WprimeENUAnalysis.HighEnergyElectronID_cff")
from WprimeAnalysis.WprimeENUAnalysis.HighEnergyElectronID_cff import *

process.load("WprimeAnalysis.WprimeENUAnalysis.HighEnergyElectronIsolation_cfi")
from WprimeAnalysis.WprimeENUAnalysis.HighEnergyElectronIsolation_cfi import *
eleIsoDepositHcalDepth1FromTowers.src = cms.InputTag("theGsfElectrons")
eleIsoDepositHcalDepth2FromTowers.src = cms.InputTag("theGsfElectrons")
electronTrackIsolation.electronProducer = cms.InputTag("theGsfElectrons")
egammaEcalRecHitIsolation.emObjectProducer = cms.InputTag("theGsfElectrons")


process.matchElectrons = cms.EDProducer( "MCTruthDeltaRMatcherNew",
        src = cms.InputTag("theGsfElectrons"),
#        src = cms.InputTag("pixelMatchGsfElectrons"),
        matched = cms.InputTag("genParticles"),
        distMin = cms.double(0.15),
        matchPDGId = cms.vint32(11)
)

process.myanalysis = cms.EDAnalyzer("WprimeAnalyzer",
        TriggerResults = cms.InputTag("TriggerResults::HLT"),    
        electrons = cms.InputTag("theGsfElectrons"),
#        electrons = cms.InputTag("pixelMatchGsfElectrons"),
        electronID = cms.untracked.string("eidRobustHighEnergy"),        
        TrkIsolationProducer = cms.InputTag("electronTrackIsolation"),
        EcalIsolationProducer = cms.InputTag("egammaEcalRecHitIsolation"),
        HcalIsolationProducerDepth1 = cms.InputTag("eleIsoFromDepsHcalDepth1FromTowers"),
        HcalIsolationProducerDepth2 = cms.InputTag("eleIsoFromDepsHcalDepth2FromTowers"),
        met = cms.InputTag("met"),
        jetAlgo = cms.InputTag("iterativeCone5CaloJets")                
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("testRelValPt1000.root")
)

process.p = cms.Path(electronCleaning+process.eIDSequence+process.eIsolationSequence+process.matchElectrons+process.myanalysis)

#process.p = cms.Path(electronCleaning * process.eIDSequence*process.eIsolationSequence*process.matchElectrons)

## process.load("Configuration.EventContent.EventContent_cff")
## process.out = cms.OutputModule("PoolOutputModule",
##     process.FEVTSIMEventContent,
##     fileName = cms.untracked.string('file:/tmp/malberti/eleIso.root')
## )

## process.out.outputCommands.append('drop *_*_*_*')
## process.out.outputCommands.append('keep *_pixelMatchGsfElectrons_*_*')
## process.out.outputCommands.append('keep *_*_*_myprocess')

## process.outpath = cms.EndPath(process.out)
