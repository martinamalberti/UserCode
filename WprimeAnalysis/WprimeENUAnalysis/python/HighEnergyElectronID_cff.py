import FWCore.ParameterSet.Config as cms

from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustHighEnergy = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustHighEnergy.src = cms.InputTag("theGsfElectrons")
eidRobustHighEnergy.robustEleIDCuts.barrel = [0.050, 0.011, 0.090, 0.005]
eidRobustHighEnergy.robustEleIDCuts.endcap = [0.100, 0.0275, 0.090, 0.007]

eIDSequence = cms.Sequence(eidRobustHighEnergy)

