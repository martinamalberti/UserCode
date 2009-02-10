#need the EleIsoHcalFromTowersExtractorBlock to configure the CandIsoDepositProducer
#the defaults are sensible
from RecoEgamma.EgammaIsolationAlgos.eleHcalExtractorBlocks_cff import *


# Hcal isolation Depth1
eleIsoDepositHcalDepth1FromTowers = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pixelMatchGsfElectrons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(EleIsoHcalFromTowersExtractorBlock)
)
eleIsoDepositHcalDepth1FromTowers.hcalDepth = cms.int32(1)

eleIsoFromDepsHcalDepth1FromTowers = cms.EDFilter("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("eleIsoDepositHcalDepth1FromTowers"), #the input isodeposits
        deltaR = cms.double(0.3), #outer cone size
        weight = cms.string('1'),
        vetos = cms.vstring('0.1'), # the veto (here vetoing inner cone of 0.1)
        skipDefaultVeto = cms.bool(True), #default veto is depreciated and usually skipped, using the vetoes specifed in veto above
        mode = cms.string('sum') #sum the Ets
    ))
)

# Hcal isolation Depth2
eleIsoDepositHcalDepth2FromTowers = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pixelMatchGsfElectrons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(EleIsoHcalFromTowersExtractorBlock)
)
eleIsoDepositHcalDepth2FromTowers.hcalDepth = cms.int32(2)

eleIsoFromDepsHcalDepth2FromTowers = cms.EDFilter("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("eleIsoDepositHcalDepth2FromTowers"),
        deltaR = cms.double(0.3),
        weight = cms.string('1'),
        vetos = cms.vstring('0.0'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


#define a sequence, produce the IsoDeposits first and then make the isolations from them
eleHcalDepthSplitIsolSequence = cms.Sequence(eleIsoDepositHcalDepth1FromTowers*eleIsoFromDepsHcalDepth1FromTowers*eleIsoDepositHcalDepth2FromTowers*eleIsoFromDepsHcalDepth2FromTowers)  


#Track isolation 
electronTrackIsolation = cms.EDProducer("EgammaElectronTkIsolationProducer",
    absolut = cms.bool(True),
    trackProducer = cms.InputTag("generalTracks"),
    intRadius = cms.double(0.02),
    electronProducer = cms.InputTag("pixelMatchGsfElectrons"),
    extRadius = cms.double(0.2),
    ptMin = cms.double(1.5),
    maxVtxDist = cms.double(0.1),
    BeamspotProducer = cms.InputTag("offlineBeamSpot"),
    maxVtxDistXY     = cms.double(99999.)
)

#ECAL isolation (based on RecHits)
egammaEcalRecHitIsolation = cms.EDProducer("EgammaEcalRecHitIsolationProducer",
    
    ecalBarrelRecHitProducer = cms.InputTag("ecalRecHit"),
    ecalEndcapRecHitCollection = cms.InputTag("EcalRecHitsEE"),
    ecalEndcapRecHitProducer = cms.InputTag("ecalRecHit"),
    ecalBarrelRecHitCollection = cms.InputTag("EcalRecHitsEB"),

    intRadiusBarrel = cms.double(0.045),
    intRadiusEndcap = cms.double(0.07),
    extRadius = cms.double(0.3),
    etMinBarrel = cms.double(-9999),
    eMinBarrel = cms.double(0.08),
    etMinEndcap = cms.double(-9999),
    eMinEndcap = cms.double(0.3),
    jurassicWidth = cms.double(0.02),

    useIsolEt = cms.bool(True), # choose if ET or E isolation
    subtract  = cms.bool(False), # subtract SC hits
    tryBoth   = cms.bool(True), # barrel and endcaps
                                           
    emObjectProducer = cms.InputTag("pixelMatchGsfElectrons")
)

eIsolationSequence = cms.Sequence(eleHcalDepthSplitIsolSequence+electronTrackIsolation+egammaEcalRecHitIsolation)
