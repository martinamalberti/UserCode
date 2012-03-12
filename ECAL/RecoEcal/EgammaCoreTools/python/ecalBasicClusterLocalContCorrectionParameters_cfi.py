import FWCore.ParameterSet.Config as cms

ecalBasicClusterLocalContCorrectionParameters = cms.PSet(
    
 paramsEta0 = cms.vdouble(1.00603,0.00300789 , 0.0667232),
 paramsEta1 = cms.vdouble(1.00655,0.00386189 , 0.073931),
 paramsEta2 = cms.vdouble(1.00634,0.00631341 , 0.0764134),
 paramsEta3 = cms.vdouble(1.00957,0.0113306 , 0.123808),

 paramsPhi0 = cms.vdouble(1.00403,-0.0012733 , 0.042925),
 paramsPhi1 = cms.vdouble(1.00394,-0.00137567 , 0.0416698),
 paramsPhi2 = cms.vdouble(1.00298,-0.00111589 , 0.0320377),
 paramsPhi3 = cms.vdouble(1.00269,-0.00153347 , 0.0296769),

)
