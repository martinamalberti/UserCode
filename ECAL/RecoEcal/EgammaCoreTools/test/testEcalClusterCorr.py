from FWCore.ParameterSet.Config import *

process = Process("test")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1


process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.source = Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data2/amassiro/CMSSWRoot/WJetFall11/WJetFall11_Example.root'
    )
)

process.maxEvents = untracked.PSet( input = untracked.int32( 10 ) )

from RecoEcal.EgammaCoreTools.ecalBasicClusterLocalContCorrectionParameters_cfi import *

process.testEcalClusterCorrections = EDAnalyzer("testEcalClusterCorrections",
    barrelClusterCollection = InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    localCorrectionParameters = ecalBasicClusterLocalContCorrectionParameters
                                                
)

process.p1 = Path( process.testEcalClusterCorrections )
