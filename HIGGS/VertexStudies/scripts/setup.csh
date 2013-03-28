setenv THISDIR `pwd`

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/../VertexAnalysis/lib

if (${?DYLD_LIBRARY_PATH}) then
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}/lib
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}/../VertexAnalysis/lib
endif

setenv PATH ${PATH}:${THISDIR}/bin
setenv PATH ${PATH}:${THISDIR}/../NtuplePackage/bin
setenv PATH ${PATH}:${THISDIR}/../VertexAnalysis/bin

setenv NTUPLEPKGINCLUDE ${THISDIR}/../NtuplePackage/interface
setenv VERTEXANALYSISINCLUDE ${THISDIR}/../VertexAnalysis/interface
setenv NTUPLEPKGLIB ${THISDIR}/../NtuplePackage/lib
setenv VERTEXANALYSISLIB ${THISDIR}/../VertexAnalysis/lib
