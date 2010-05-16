import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# Geometry
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_1.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_2.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_3.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_4.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_5.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_6.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_7.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_8.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_9.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_10.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_11.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_12.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_13.root",
    "rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03_newMask/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_newMask_14.root",
        
    #"rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_1.root",
    #"rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_2.root",
    #"rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_4.root",
    #"rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_5.root",
    #"rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_6.root",
    #"rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_7.root",
    #"rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_9.root",
    #"rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_10.root",
    #"rfio:/castor/cern.ch/user/d/deguio/REDIGI_minBias7TeV_Spring10-START3X_V25B-v1_newNoise_360_v03/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_360_v03_11.root",
            
    )
)



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# filter on PhysDeclared bit
process.skimming = cms.EDFilter("PhysDecl",
   applyfilter = cms.untracked.bool(True)
)

# filter on bit 40 || 41 nad !(bit36 || bit37 || bit38 || bit39)
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')


#Good Vertex Filter (see GOODCOLL skim)
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), 
  filter = cms.bool(True)   
)

# FilterOutScraping
process.noscraping = cms.EDFilter("FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.25)
)

process.load("Validation.EcalValidation.ecalvalidation_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('EcalValidation.root')
)

process.p = cms.Path(
    #process.skimming*
    process.hltLevel1GTSeed*
    #process.noscraping*
    process.primaryVertexFilter*
    process.ecalvalidation
    )

