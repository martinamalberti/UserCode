/********************************************************************
* ../obj/mydict.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error ../obj/mydict.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtablemydict();
extern void G__cpp_setup_inheritancemydict();
extern void G__cpp_setup_typetablemydict();
extern void G__cpp_setup_memvarmydict();
extern void G__cpp_setup_globalmydict();
extern void G__cpp_setup_memfuncmydict();
extern void G__cpp_setup_funcmydict();
extern void G__set_cpp_environmentmydict();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "../interface/MyTest.h"
#include "../interface/readJSONFile.h"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__mydictLN_vectorlEcharcOallocatorlEchargRsPgR;
extern G__linked_taginfo G__mydictLN_vectorlEfloatcOallocatorlEfloatgRsPgR;
extern G__linked_taginfo G__mydictLN_vectorlEdoublecOallocatorlEdoublegRsPgR;
extern G__linked_taginfo G__mydictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR;
extern G__linked_taginfo G__mydictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_vectorlEintcOallocatorlEintgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRcOallocatorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRcOallocatorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRcOallocatorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRcOallocatorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_maplEstringcOvectorlEdoublecOallocatorlEdoublegRsPgRmUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOvectorlEdoublecOallocatorlEdoublegRsPgRmUgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_maplEstringcOvectorlEfloatcOallocatorlEfloatgRsPgRmUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOvectorlEfloatcOallocatorlEfloatgRsPgRmUgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_maplEstringcOvectorlEintcOallocatorlEintgRsPgRmUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOvectorlEintcOallocatorlEintgRsPgRmUgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_maplEstringcOvectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRcOallocatorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgRsPgRmUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOvectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRcOallocatorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgRsPgRmUgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_maplEstringcOvectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRcOallocatorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgRsPgRmUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOvectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRcOallocatorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgRsPgRmUgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_maplEstringcOTClonesArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTClonesArraymUgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_TMatrixTSymlEdoublegR;
extern G__linked_taginfo G__mydictLN_TMatrixTlEdoublegR;
extern G__linked_taginfo G__mydictLN_TVectorTlEdoublegR;
extern G__linked_taginfo G__mydictLN_maplETStringcOTMVAcLcLTypescLcLEMVAcOlesslETStringgRcOallocatorlEpairlEconstsPTStringcOTMVAcLcLTypescLcLEMVAgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_TVectorTlEfloatgR;
extern G__linked_taginfo G__mydictLN_vectorlEfloatmUcOallocatorlEfloatmUgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlEfloatmUcOallocatorlEfloatmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_vectorlETMVAcLcLVariableInfocOallocatorlETMVAcLcLVariableInfogRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlETMVAcLcLVariableInfocOallocatorlETMVAcLcLVariableInfogRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_vectorlETStringcOallocatorlETStringgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_vectorlETMVAcLcLClassInfomUcOallocatorlETMVAcLcLClassInfomUgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlETMVAcLcLClassInfomUcOallocatorlETMVAcLcLClassInfomUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_pairlEcharcOunsignedsPintgR;
extern G__linked_taginfo G__mydictLN_vectorlEpairlEcharcOunsignedsPintgRcOallocatorlEpairlEcharcOunsignedsPintgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlEpairlEcharcOunsignedsPintgRcOallocatorlEpairlEcharcOunsignedsPintgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_vectorlETMVAcLcLEventmUcOallocatorlETMVAcLcLEventmUgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlETMVAcLcLEventmUcOallocatorlETMVAcLcLEventmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_vectorlETMatrixTSymlEdoublegRmUcOallocatorlETMatrixTSymlEdoublegRmUgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlETMatrixTSymlEdoublegRmUcOallocatorlETMatrixTSymlEdoublegRmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_vectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_maplETMVAcLcLEMsgTypecOstringcOlesslETMVAcLcLEMsgTypegRcOallocatorlEpairlEconstsPTMVAcLcLEMsgTypecOstringgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_vectorlETMVAcLcLTreeInfocOallocatorlETMVAcLcLTreeInfogRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlETMVAcLcLTreeInfocOallocatorlETMVAcLcLTreeInfogRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_maplETStringcOvectorlETMVAcLcLTreeInfocOallocatorlETMVAcLcLTreeInfogRsPgRcOlesslETStringgRcOallocatorlEpairlEconstsPTStringcOvectorlETMVAcLcLTreeInfocOallocatorlETMVAcLcLTreeInfogRsPgRsPgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_maplEstringcOboolcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOboolgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_vectorlEstringcOallocatorlEstringgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlEstringcOallocatorlEstringgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_maplETStringcOTMVAcLcLIMethodmUcOlesslETStringgRcOallocatorlEpairlEconstsPTStringcOTMVAcLcLIMethodmUgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_vectorlEpairlEintcOintgRcOallocatorlEpairlEintcOintgRsPgRsPgR;
extern G__linked_taginfo G__mydictLN_reverse_iteratorlEvectorlEpairlEintcOintgRcOallocatorlEpairlEintcOintgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__mydictLN_maplEintcOvectorlEpairlEintcOintgRcOallocatorlEpairlEintcOintgRsPgRsPgRcOlesslEintgRcOallocatorlEpairlEconstsPintcOvectorlEpairlEintcOintgRcOallocatorlEpairlEintcOintgRsPgRsPgRsPgRsPgRsPgR;

/* STUB derived class for protected member access */
