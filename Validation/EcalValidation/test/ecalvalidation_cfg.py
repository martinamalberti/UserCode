import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

# Geometry
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(

    '/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0013/F6E3FFAC-6113-DF11-AF06-001BFCDBD15E.root',
    '/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0012/BA183BD4-3813-DF11-B015-00304867902C.root',
    '/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0012/A6D53A27-3813-DF11-B460-0026189437FA.root'
    
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


process.myvalidation = cms.EDFilter("EcalValidation",
    superClusterCollection_EB = cms.InputTag("correctedHybridSuperClusters"),
    superClusterCollection_EE = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
    basicClusterCollection_EE = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"),
    recHitCollection_EE       = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    recHitCollection_EB       = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    basicClusterCollection_EB = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    recHitCollection_ES       = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
    ClusterCollectionX_ES      = cms.InputTag("multi5x5SuperClustersWithPreshower","preshowerXClusters"),
    ClusterCollectionY_ES      = cms.InputTag("multi5x5SuperClustersWithPreshower","preshowerYClusters"),
                                    
    ethrEB = cms.double(-999.),
    ethrEE = cms.double(-999.)
                                    
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_validation.root')
)

process.p = cms.Path(process.myvalidation)

