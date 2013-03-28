export THISDIR=`pwd`

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THISDIR}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THISDIR}/../VertexAnalysis/lib

if [ -n "${DYLD_LIBRARY_PATH}" ] ; then
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${THISDIR}/lib
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${THISDIR}/../VertexAnalysis/lib
fi

export PATH=${PATH}:${THISDIR}/bin
export PATH=${PATH}:${THISDIR}/../NtuplePackage/bin
export PATH=${PATH}:${THISDIR}/../VertexAnalysis/bin

export NTUPLEPKGINCLUDE=${THISDIR}/../NtuplePackage/interface
export VERTEXANALYSISINCLUDE=${THISDIR}/../VertexAnalysis/interface
export NTUPLEPKGLIB=${THISDIR}/../NtuplePackage/lib
export VERTEXANALYSISLIB=${THISDIR}/../VertexAnalysis/lib

