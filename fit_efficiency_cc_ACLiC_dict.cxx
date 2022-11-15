// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedImdiottidIalicedIhypertriton_kinkdIfit_efficiency_cc_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/mdiotti/alice/hypertriton_kink/./fit_efficiency.cc"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_fit_efficiency_cc_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./fit_efficiency.cc",
nullptr
    };
    static const char* includePaths[] = {
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-04-patches-alice2-4/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/O2Physics/master-local2/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/MCStepLogger/v0.5.0-16/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/KFParticle/v1.1-5-2/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/O2/dev-local2/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/O2/dev-local2/include/GPU",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ONNXRuntime/v1.12.1-alice1-2/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/fastjet/v3.4.0_1.045-alice1-9/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/JAliEn-ROOT/0.6.8-5/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/Common-O2/v1.6.1-1/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/Monitoring/v3.14.1-3/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/HepMC3/3.2.5-12/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/Vc/1.4.1-6/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/FairRoot/v18.4.9-7/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/FairMQ/v1.4.55-2/include/fairmq",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/FairMQ/v1.4.55-2/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/asio/v1.19.1-6/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/GEANT3/v4-1-16/include/TGeant3",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/GEANT4_VMC/v6-1-p2-3/include/g4root",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/GEANT4_VMC/v6-1-p2-3/include/geant4vmc",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/VMC/v2-0-16/include/vmc",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/TBB/v2021.5.0-6/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/libInfoLogger/v2.5.3-1/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ms_gsl/4.0.0-3/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/FairLogger/v1.11.1-2/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/fmt/9.1.0-1/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/pythia/v8304-18/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/boost/v1.75.0-17/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/OpenSSL/v1.1.1m-1/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-04-patches-alice2-4/etc/",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-04-patches-alice2-4/etc//cling",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-04-patches-alice2-4/etc//cling/plugins/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-04-patches-alice2-4/include/",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-04-patches-alice2-4/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-04-patches-alice2-4/include/",
"/home/mdiotti/alice/hypertriton_kink/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "fit_efficiency_cc_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace o2{namespace dataformats{class __attribute__((annotate("$clingAutoload$ReconstructionDataFormats/VtxTrackIndex.h")))  __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  VtxTrackIndex;}}
using GIndex __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  = o2::dataformats::VtxTrackIndex;
namespace o2{namespace dataformats{class __attribute__((annotate("$clingAutoload$ReconstructionDataFormats/V0.h")))  __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  V0;}}
using V0 __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  = o2::dataformats::V0;
namespace o2{template <class _T> class __attribute__((annotate("$clingAutoload$SimulationDataFormat/MCTrack.h")))  __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  MCTrackT;
}
namespace o2{using MCTrack __attribute__((annotate("$clingAutoload$SimulationDataFormat/MCTrack.h")))  __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  = MCTrackT<float>;}
using MCTrack __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  = o2::MCTrack;
namespace o2{namespace dataformats{template <typename FirstEntry = int, typename NElem = int> class __attribute__((annotate("$clingAutoload$CommonDataFormat/RangeReference.h")))  __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  RangeReference;
}}
using RRef __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  = o2::dataformats::RangeReference<int, int>;
namespace o2{namespace its{class __attribute__((annotate("$clingAutoload$DataFormatsITS/TrackITS.h")))  __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  TrackITS;}}
using TrackITS __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  = o2::its::TrackITS;
namespace o2{class __attribute__((annotate("$clingAutoload$SimulationDataFormat/MCCompLabel.h")))  __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  MCCompLabel;}
namespace o2{namespace dataformats{template <typename TruthElement> class __attribute__((annotate("$clingAutoload$SimulationDataFormat/MCTruthContainer.h")))  __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  MCTruthContainer;
}}
using MCLabCont __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
namespace o2{namespace itsmft{class __attribute__((annotate("$clingAutoload$DataFormatsITSMFT/CompCluster.h")))  __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  CompClusterExt;}}
using CompClusterExt __attribute__((annotate("$clingAutoload$./fit_efficiency.cc")))  = o2::itsmft::CompClusterExt;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "fit_efficiency_cc_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./fit_efficiency.cc"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"CompClusterExt", payloadCode, "@",
"FITTEROPTION", payloadCode, "@",
"GIndex", payloadCode, "@",
"MCLabCont", payloadCode, "@",
"MCTrack", payloadCode, "@",
"RRef", payloadCode, "@",
"TrackITS", payloadCode, "@",
"V0", payloadCode, "@",
"Vec3", payloadCode, "@",
"calcRadius", payloadCode, "@",
"chiLabel", payloadCode, "@",
"fit_efficiency", payloadCode, "@",
"fontSize", payloadCode, "@",
"hypLabel", payloadCode, "@",
"hypMassTh", payloadCode, "@",
"hypPDG", payloadCode, "@",
"markerSize", payloadCode, "@",
"nBins", payloadCode, "@",
"pi0Mass", payloadCode, "@",
"tritonMass", payloadCode, "@",
"tritonPDG", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("fit_efficiency_cc_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_fit_efficiency_cc_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_fit_efficiency_cc_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_fit_efficiency_cc_ACLiC_dict() {
  TriggerDictionaryInitialization_fit_efficiency_cc_ACLiC_dict_Impl();
}
