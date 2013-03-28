HDR   = ../interface/
OBJ   = ../obj/
LIB   = ../lib/
BIN   =  ../bin/
PRG  =  ../test/


HdrSuf  =  .h
SrcSuf  =  .cc
ObjSuf  =  .o
PrgSuf  =  .cpp
BinSuf  =  .exe
LibSuf  =  .so


HDRS  =  $(wildcard $(HDR)*$(HdrSuf))
SRCS  =  $(wildcard *$(SrcSuf))
_OBJS  =  $(patsubst %$(SrcSuf), %$(ObjSuf), $(SRCS))
OBJS   =  $(addprefix $(OBJ),$(_OBJS))
PRGS  =  $(wildcard $(PRG)*$(PrgSuf))
_BINS    =  $(wildcard $(PRG)*$(PrgSuf))
__BINS   =  $(_BINS:$(PrgSuf)=$(BinSuf))
___BINS  =  $(notdir $(__BINS))
BINS	 =  $(addprefix $(BIN),${___BINS})

LINKDEF   =  $(wildcard ${HDR}*LinkDef.h)
DICTHDRS  =  $(patsubst $(LINKDEF),,$(HDRS)) $(LINKDEF)


ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CLHEPCFLAGS   =  $(shell clhep-config --include)
CLHEPLIBS     =  $(shell clhep-config --libs)

NTUPLEPACKAGEFLAGS = $(shell ntuplepkg-config --include)
NTUPLEPACKAGELIBS = $(shell ntuplepkg-config --libs)

VERTEXANALYSISFLAGS = $(shell vtxpkg-config --include)
VERTEXANALYSISLIBS = $(shell vtxpkg-config --libs)

ARCHL = -m32
ARCH  =  $(shell root-config --arch)
ifeq ($(ARCH),linuxx8664gcc)
ARCHL = -m64
endif

CXX  =  c++
CXXFLAGS  =  -Wall -O -fPIC -I$(HDR) $(ROOTCFLAGS) $(NTUPLEPACKAGEFLAGS) $(VERTEXANALYSISFLAGS)
CPP  =  c++
CPPFLAGS  = -Wall $(ARCHL) -I$(HDR) $(ROOTCFLAGS) $(NTUPLEPACKAGEFLAGS) $(VERTEXANALYSISFLAGS)


LD       =  c++
LDFLAGS  =  -rdynamic -shared $(ARCHL)
SONAME	 =  libMyTest.so
SOFLAGS  =  -Wl,-soname,$(SONAME)

GLIBS   =  -lm -ldl -rdynamic $(ROOTGLIBS) -lGenVector  -lTMVA $(NTUPLEPACKAGELIBS) $(VERTEXANALYSISLIBS)

#################################################
#if mac
ARCH  =  $(shell root-config --arch)
ifeq ($(ARCH),macosx64)
LibSuf  =  .dylib

CPPFLAGS  =  -Wall -W -Woverloaded-virtual -O2 $(ARCHL) -pipe -I$(HDR) $(ROOTCFLAGS) $(NTUPLEPACKAGEFLAGS) $(VERTEXANALYSISFLAGS)

CXXFLAGS  =  -Wall -W -Woverloaded-virtual -O2 $(ARCHL) -pipe -I$(HDR) $(ROOTCFLAGS) $(NTUPLEPACKAGEFLAGS) $(VERTEXANALYSISFLAGS)

LDFLAGS  =  -dynamiclib -shared -single_module -undefined dynamic_lookup $(ARCHL)
SONAME	 =  libMyTest.dylib
SOFLAGS  =
endif
#################################################





.PHONY: all exe test clean


all: $(LIB)$(SONAME)


exe: $(BINS) 


test:
	@echo "HDRS = $(HDRS)"
	@echo "DICTHDRS = $(DICTHDRS)"
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "PRGS = $(PRGS)"
	@echo "BINS = $(BINS)"


$(OBJ)%$(ObjSuf): %$(SrcSuf) $(HDRS)
	$(CXX) -c $(CXXFLAGS) -o $@ $< 


$(OBJ)mydict.o: $(OBJ)mydict.cc $(OBJ)mydict.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<


$(OBJ)mydict.cc: $(DICTHDRS)
	@echo "Generating dictionary for  ..."
	rootcint -f $(OBJ)mydict.cc -c -p ${CXXFLAGS} $(DICTHDRS)

$(LIB)$(SONAME): $(OBJS) $(OBJ)mydict.o
	@echo "Linking $(SONAME):"
	$(LD) $(LDFLAGS) $(OBJS) $(OBJ)mydict.o $(SOFLAGS) -o $(LIB)$(SONAME)

$(BIN)%$(BinSuf): $(PRG)%$(PrgSuf) $(HDRS) $(LIB)$(SONAME)
	$(CPP) $(CPPFLAGS) -L$(LIB) $(GLIBS) -lMyTest  -o $@ $<


clean:
	rm -f $(OBJ)*$(ObjSuf) $(LIB)*$(LibSuf) $(OBJ)mydict*
