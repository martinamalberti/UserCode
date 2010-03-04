import FWCore.ParameterSet.Config as cms

ecalvalidation = cms.EDFilter("EcalValidation",
    superClusterCollection_EB = cms.InputTag("correctedHybridSuperClusters"),
    superClusterCollection_EE = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
    basicClusterCollection_EE = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"),
    recHitCollection_EE       = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    recHitCollection_EB       = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    basicClusterCollection_EB = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    recHitCollection_ES       = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
    ClusterCollectionX_ES     = cms.InputTag("multi5x5SuperClustersWithPreshower","preshowerXClusters"),
    ClusterCollectionY_ES     = cms.InputTag("multi5x5SuperClustersWithPreshower","preshowerYClusters"),
                                    
    ethrEB = cms.double(0.800),
    ethrEE = cms.double(1.200),
                                    
    # for pi0

    isMonEBpi0 = cms.untracked.bool(True),
    selePtGamma = cms.double(1 ),
    selePtPi0 = cms.double( 2. ),
    #seleMinvMaxPi0 = cms.double( 0.22 ),
    seleMinvMaxPi0 = cms.double( 0.35 ),
    seleMinvMinPi0 = cms.double( 0.06 ),
    seleS4S9Gamma = cms.double( 0.83 ),
    selePi0Iso = cms.double( 0.5 ),
    ptMinForIsolation = cms.double( 1 ),
    selePi0BeltDR = cms.double( 0.2 ),
    selePi0BeltDeta = cms.double( 0.05 ),

    isMonEEpi0 = cms.untracked.bool(True),
    selePtGamma_EE = cms.double( 0.8 ),
    selePtPi0_EE = cms.double( 3.0 ),
    seleS4S9Gamma_EE = cms.double( 0.9 ),
    seleMinvMaxPi0_EE = cms.double( 0.3 ),
    seleMinvMinPi0_EE = cms.double( 0.05 ),
    ptMinForIsolation_EE = cms.double( 0.5 ),
    selePi0Iso_EE = cms.double( 0.5 ),
    selePi0BeltDR_EE  = cms.double( 0.2 ),
    selePi0BeltDeta_EE  = cms.double( 0.05 ),

    region1_Pi0_EE = cms.double(2),
    selePtGammaPi0_EE_region1 = cms.double(0.7),
    selePtPi0_EE_region1 = cms.double(3),

    region2_Pi0_EE = cms.double(2.5),
    selePtGammaPi0_EE_region2 = cms.double(0.5),
    selePtPi0_EE_region2 = cms.double(2),

    selePtGammaPi0_EE_region3 = cms.double(0.3),
    selePtPi0_EE_region3 = cms.double(1.2),

    clusSeedThr = cms.double( 0.5 ),
    clusSeedThr_EE = cms.double( 1.0 ),
    clusEtaSize = cms.int32( 3 ),
    clusPhiSize = cms.int32( 3 ),
    seleXtalMinEnergy = cms.double( -0.15 ),
    seleXtalMinEnergy_EE = cms.double( -0.75 ),
                                    
    isMaskEB = cms.untracked.bool(True),
    isMaskEE = cms.untracked.bool(False),

    maskEBFile = cms.untracked.string('maskEB.txt'),
    maskEEFile = cms.untracked.string('maskEE.txt'),

    useRecoFlag = cms.untracked.bool(False),

    ParameterLogWeighted = cms.bool( True ),
    ParameterX0 = cms.double( 0.89 ),
    ParameterT0_barl = cms.double( 5.7 ),
    ParameterT0_endc = cms.double( 3.1 ),
    ParameterT0_endcPresh = cms.double( 1.2 ),
    ParameterW0 = cms.double( 4.2 )
                                 
)

