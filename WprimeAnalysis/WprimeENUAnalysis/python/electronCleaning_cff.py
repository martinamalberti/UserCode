import FWCore.ParameterSet.Config as cms

# To remove electron duplicates 
theGsfElectrons = cms.EDFilter("ElectronDuplicateRemover",
    src = cms.untracked.string('pixelMatchGsfElectrons'),
    #ptMin = cms.untracked.double(20.0),
    #EndcapMinEta = cms.untracked.double(1.56),
    #ptMax = cms.untracked.double(10000.0),
    #BarrelMaxEta = cms.untracked.double(1.4442),
    #EndcapMaxEta = cms.untracked.double(3.0)

    ptMin = cms.untracked.double(0.0),
    EndcapMinEta = cms.untracked.double(1.4),
    ptMax = cms.untracked.double(10000.0),
    BarrelMaxEta = cms.untracked.double(1.5),
    EndcapMaxEta = cms.untracked.double(3.0)                         
)

## gsfSelection = cms.EDFilter("GsfElectronSelector",
##      src = cms.InputTag("gsfElectrons"),
##      cut = cms.string('et > 0.0')
## )

## theGsfElectrons = cms.EDProducer("PixelMatchGsfElectronShallowCloneProducer",
##     src = cms.InputTag("gsfElectrons")
## )


electronCleaning = cms.Sequence(theGsfElectrons)
