setenv THISDIR `pwd`
setenv H2GGLOBEDIR /afs/cern.ch/work/m/malberti/HGG/LegacyPaper/CMSSW_5_3_10/src/h2gglobe

echo $THISDIR
echo $H2GGLOBEDIR

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/../VertexAnalysis/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${H2GGLOBEDIR}/VertexAnalysis/lib

if (${?DYLD_LIBRARY_PATH}) then
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}/lib
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
#setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}/../VertexAnalysis/lib
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${H2GGLOBEDIR}/VertexAnalysis/lib
endif

setenv PATH ${PATH}:${THISDIR}/bin
setenv PATH ${PATH}:${THISDIR}/../NtuplePackage/bin
#setenv PATH ${PATH}:${THISDIR}/../VertexAnalysis/bin
setenv PATH ${PATH}:${H2GGLOBEDIR}/VertexAnalysis/bin

setenv NTUPLEPKGINCLUDE ${THISDIR}/../NtuplePackage/interface
#setenv VERTEXANALYSISINCLUDE ${THISDIR}/../VertexAnalysis/interface
setenv VERTEXANALYSISINCLUDE ${H2GGLOBEDIR}/VertexAnalysis/interface

setenv NTUPLEPKGLIB ${THISDIR}/../NtuplePackage/lib
#setenv VERTEXANALYSISLIB ${THISDIR}/../VertexAnalysis/lib
setenv VERTEXANALYSISLIB ${H2GGLOBEDIR}/VertexAnalysis/lib
