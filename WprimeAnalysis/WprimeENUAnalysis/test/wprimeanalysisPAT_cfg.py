import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START3X_V26A::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )

process.source = cms.Source("PoolSource",
       duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),    
       fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_5_7/RelValWE/GEN-SIM-RECO/START3X_V26-v1/0012/EC9F278B-6949-DF11-99F2-003048678B0E.root',
        #'/store/relval/CMSSW_3_5_7/RelValWE/GEN-SIM-RECO/START3X_V26-v1/0012/EAB2D968-4649-DF11-9FDB-003048679166.root',
        #'/store/relval/CMSSW_3_5_7/RelValWE/GEN-SIM-RECO/START3X_V26-v1/0012/E2E0DBBC-4749-DF11-B841-003048678FC4.root',
        #'/store/relval/CMSSW_3_5_7/RelValWE/GEN-SIM-RECO/START3X_V26-v1/0012/60FB455C-4649-DF11-AF23-003048678B0C.root',
        #'/store/relval/CMSSW_3_5_7/RelValWE/GEN-SIM-RECO/START3X_V26-v1/0012/2275C211-4749-DF11-A7D8-003048678FFE.root'
        
                                                )
)


#no configuration of the pat is necessary for us at the moment
process.load("PhysicsTools.PatAlgos.patSequences_cff");

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                          #eleLabel = cms.InputTag("allLayer1Electrons"),
                                          eleLabel   = cms.InputTag("patElectrons"),
                                          barrelCuts = cms.PSet(heepBarrelCuts),
                                          endcapCuts = cms.PSet(heepEndcapCuts)
                                          )


#Analysis
process.myanalysis = cms.EDAnalyzer('WprimeTreeHEEP',
                                    recHitCollection_EB       = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                                    electronTag = cms.InputTag("heepPatElectrons"),
                                    jetTag      = cms.InputTag("patJets"),
                                    metTag      = cms.InputTag("patMETs"),
                                    muonTag     = cms.InputTag("patMuons"),
                                    electronID  = cms.untracked.string("eidRobustHighEnergy"),
                                    btagAlgo    = cms.untracked.string("jetBProbabilityBJetTags"),
                                    HLTInputTag = cms.InputTag("TriggerResults","","HLT"),
                                    L1InputTag =  cms.InputTag("gtDigis")
                                    )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("testHEEP.root")
)


process.p = cms.Path( process.patDefaultSequence *
                      process.heepPatElectrons* #heepifies the pat electrons (resets energy to ecal energy and adds heep id)
                      process.myanalysis) 



## process.load("Configuration.EventContent.EventContent_cff")
## process.out = cms.OutputModule("PoolOutputModule",
##             process.FEVTSIMEventContent,
##             fileName = cms.untracked.string('file:/tmp/malberti/testPat.root')
## )

## process.out.outputCommands.append('drop *_*_*_*')
## process.out.outputCommands.append('keep *_*_*_myprocess')
## process.out.outputCommands.append('keep *_*_*_HLT')

## process.outpath = cms.EndPath(process.wprimePATSequence *  process.out)

