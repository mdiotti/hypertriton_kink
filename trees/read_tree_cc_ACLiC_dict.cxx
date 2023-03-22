// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedImdiottidIalicedIcodicidIread_tree_cc_ACLiC_dict
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
#include "/home/mdiotti/alice/codici/./read_tree.cc"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_read_tree_cc_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./read_tree.cc",
nullptr
    };
    static const char* includePaths[] = {
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-10-alice5-5/include",
"/home/mdiotti/alice/codici",
"/home/mdiotti/Desktop/Repositories/Utils/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/O2/str_tr-local1/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ms_gsl/4.0.0-4/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/O2/str_tr-local1/include/GPU",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/O2Physics/master-local5/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/MCStepLogger/v0.5.0-27/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/KFParticle/v1.1-5-13/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/O2/str_tr-local1/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/O2/str_tr-local1/include/GPU",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ONNXRuntime/v1.12.1-alice1-9/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/fastjet/v3.4.0_1.045-alice1-15/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/JAliEn-ROOT/0.6.9-6/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/Common-O2/v1.6.1-7/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/Monitoring/v3.15.1-3/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/HepMC3/3.2.5-23/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/Vc/1.4.1-7/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/FairRoot/v18.4.9-20/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/FairMQ/v1.4.56-3/include/fairmq",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/FairMQ/v1.4.56-3/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/asio/v1.19.1-7/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/GEANT3/v4-2-3/include/TGeant3",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/GEANT4_VMC/v6-1-p3-3/include/g4root",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/GEANT4_VMC/v6-1-p3-3/include/geant4vmc",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/VMC/v2-0-27/include/vmc",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/TBB/v2021.5.0-7/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/libInfoLogger/v2.5.3-7/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ms_gsl/4.0.0-4/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/FairLogger/v1.11.1-3/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/fmt/9.1.0-2/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/pythia/v8304-24/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/boost/v1.75.0-23/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/OpenSSL/v1.1.1m-2/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-10-alice5-5/etc/",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-10-alice5-5/etc//cling",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-10-alice5-5/etc//cling/plugins/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-10-alice5-5/include/",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-10-alice5-5/include",
"/home/mdiotti/alice/codici",
"/home/mdiotti/Desktop/Repositories/Utils/include",
"/home/mdiotti/alice/sw/ubuntu2004_x86-64/ROOT/v6-26-10-alice5-5/include/",
"/home/mdiotti/alice/codici/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "read_tree_cc_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace o2{namespace dataformats{class __attribute__((annotate("$clingAutoload$ReconstructionDataFormats/VtxTrackIndex.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  VtxTrackIndex;}}
using GIndex __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::dataformats::VtxTrackIndex;
namespace o2{namespace dataformats{class __attribute__((annotate("$clingAutoload$ReconstructionDataFormats/V0.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  V0;}}
using V0 __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::dataformats::V0;
namespace o2{template <class _T> class __attribute__((annotate("$clingAutoload$SimulationDataFormat/MCTrack.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  MCTrackT;
}
namespace o2{using MCTrack __attribute__((annotate("$clingAutoload$SimulationDataFormat/MCTrack.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = MCTrackT<float>;}
using MCTrack __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::MCTrack;
namespace o2{namespace dataformats{template <typename FirstEntry = int, typename NElem = int> class __attribute__((annotate("$clingAutoload$CommonDataFormat/RangeReference.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  RangeReference;
}}
using RRef __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::dataformats::RangeReference<int, int>;
namespace o2{namespace its{class __attribute__((annotate("$clingAutoload$DataFormatsITS/TrackITS.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  TrackITS;}}
using TrackITS __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::its::TrackITS;
namespace o2{namespace dataformats{class __attribute__((annotate("$clingAutoload$ReconstructionDataFormats/TrackTPCITS.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  TrackTPCITS;}}
using TrackTPCITS __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::dataformats::TrackTPCITS;
namespace o2{namespace dataformats{struct __attribute__((annotate("$clingAutoload$ReconstructionDataFormats/StrangeTrack.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  StrangeTrack;}}
using StrangeTrack __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::dataformats::StrangeTrack;
namespace o2{namespace dataformats{struct __attribute__((annotate("$clingAutoload$ReconstructionDataFormats/KinkTrack.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  KinkTrack;}}
using KinkTrack __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::dataformats::KinkTrack;
namespace o2{class __attribute__((annotate("$clingAutoload$SimulationDataFormat/MCCompLabel.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  MCCompLabel;}
namespace o2{namespace dataformats{template <typename T> class __attribute__((annotate("$clingAutoload$ITStracking/IOUtils.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  MCTruthContainer;
}}
using MCLabCont __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
namespace o2{namespace itsmft{class __attribute__((annotate("$clingAutoload$DataFormatsITSMFT/CompCluster.h")))  __attribute__((annotate("$clingAutoload$./read_tree.cc")))  CompClusterExt;}}
using CompClusterExt __attribute__((annotate("$clingAutoload$./read_tree.cc")))  = o2::itsmft::CompClusterExt;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "read_tree_cc_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./read_tree.cc"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"CompClusterExt", payloadCode, "@",
"GIndex", payloadCode, "@",
"KinkTrack", payloadCode, "@",
"MCLabCont", payloadCode, "@",
"MCTrack", payloadCode, "@",
"RRef", payloadCode, "@",
"StrangeTrack", payloadCode, "@",
"TrackITS", payloadCode, "@",
"TrackTPCITS", payloadCode, "@",
"V0", payloadCode, "@",
"VBracket", payloadCode, "@",
"Vec3", payloadCode, "@",
"calcDecayLenght", payloadCode, "@",
"calcRadius", payloadCode, "@",
"hypLabel", payloadCode, "@",
"hypPDG", payloadCode, "@",
"nBins", payloadCode, "@",
"pi0Mass", payloadCode, "@",
"ptLabel", payloadCode, "@",
"read_tree", payloadCode, "@",
"tritonMass", payloadCode, "@",
"tritonPDG", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("read_tree_cc_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_read_tree_cc_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_read_tree_cc_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_read_tree_cc_ACLiC_dict() {
  TriggerDictionaryInitialization_read_tree_cc_ACLiC_dict_Impl();
}
