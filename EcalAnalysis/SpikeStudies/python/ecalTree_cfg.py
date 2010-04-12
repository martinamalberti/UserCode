import FWCore.ParameterSet.Config as cms

process = cms.Process("ECALANALYSIS")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)


# Geometry
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

# Global Tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.GlobalTag.globaltag = 'GR10_P_V4::All'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
      #'file:/tmp/malberti/skimMonoPhoton_evt20176336_123596_ReReco.root'
#    'rfio:/castor/cern.ch/user/m/malberti/Collisions09/MonsterRecHitsSkim/ecalAnomalousRecHitsSkim_6.root'   
     #'file:/tmp/malberti/ecalAnomalousRecHitsSkim_99.root'

    #/MinimumBias/Commissioning10-GOODCOLL-v7/RAW-RECO
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/FCCA95EA-4A3C-DF11-B108-00E0817918BF.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/F82798F3-4A3C-DF11-9694-001A649747B0.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/F4C2295F-433C-DF11-89C2-00E0817917A7.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/E6C5424A-4A3C-DF11-8F0C-0030486361DC.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/D823C7A8-443C-DF11-82EF-00E08179174D.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/BCF1B595-4D3C-DF11-A3AA-001A64789358.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/AA56E655-413C-DF11-9F06-003048D46296.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/8C78DC01-463C-DF11-BF2C-00E0817917A7.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/84D87194-4D3C-DF11-B686-00E0817918C5.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/5E973302-443C-DF11-8FDF-00E081791815.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/5474B830-4A3C-DF11-90DE-00E0817918C5.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/3AA920FD-4A3C-DF11-896F-003048D47A7E.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/2EE0FC21-4C3C-DF11-A18A-00E08178C0FB.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/0C061554-493C-DF11-AFA0-00E0817918BF.root',
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v7/000/132/440/0636C91D-4C3C-DF11-9A45-001A649747B0.root'
    )
)


process.myanalysis = cms.EDAnalyzer('EcalTree',
    ebRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    ebClusterCollection = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters","RECO"),
    #ebClusterCollection = cms.InputTag("multi5x5BasicClusters","multi5x5BarrelBasicClusters","EGammaCommissioning"),
    ebDigiCollection = cms.InputTag("ecalDigis","ebDigis"),
    L1InputTag =  cms.InputTag("gtDigis")                              
)

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string("treeAnomalousChannels_test.root")
    fileName = cms.string("SpikesCommissioning2010.root")
)


## process.output = cms.OutputModule("PoolOutputModule",
##     outputCommands = cms.untracked.vstring('keep *','drop *_MEtoEDMConverter_*_*'),
##     fileName = cms.untracked.string("/tmp/malberti/mytest.root"),
##       dataset = cms.untracked.PSet(
##       dataTier = cms.untracked.string('RAW-RECO'),
##       filterName = cms.untracked.string('test')
##       ),
##       SelectEvents = cms.untracked.PSet(
##       SelectEvents = cms.vstring('p')
##     )
   
##       )
                                  
## process.outpath = cms.EndPath(process.output)


process.p = cms.Path(
    process.ecalDigis
    *process.myanalysis
    )
