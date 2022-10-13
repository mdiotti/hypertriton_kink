#if !defined(CLING) || defined(ROOTCLING)
#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "DataFormatsITS/TrackITS.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"

#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"

#endif

#include <iostream>
using namespace std;

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using Vec3 = ROOT::Math::SVector<double, 3>;

using TrackITS = o2::its::TrackITS;

const int hypPDG = 1010010030;
const int tritonPDG = 1000010030;
const int pi0PDG = 111;
const int tf_min = 1;
const int tf_max = 40;
int tf_lenght = tf_max - tf_min + 1;

void check_generated_hyp(TString path, bool verbose = true)
{

    for (int tf = tf_min; tf < tf_max; tf++)
    {
        LOG(info) << "Processing TF " << tf;
        int counter_hyp = 0;
        int counter_hyp_no_kink = 0;
        int counter_pi0 = 0;

        TString tf_string = Form("%d", tf);
        TString tf_path = path + "tf" + tf_string;

        auto fMCTracks = TFile::Open(tf_path + "/sgn_" + tf_string + "_Kine.root");

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");

        // Tracks
        std::vector<MCTrack> *MCtracks = nullptr;

        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);

        // mc start

        std::vector<std::vector<MCTrack>> mcTracksMatrix;
        auto nev = treeMCTracks->GetEntriesFast();
        unsigned int nTracks[nev];

        mcTracksMatrix.resize(nev);
        for (int n = 0; n < nev; n++) // fill mcTracksMatrix
        {
            treeMCTracks->GetEvent(n);
            unsigned int size = MCtracks->size();
            nTracks[n] = size;

            mcTracksMatrix[n].resize(size);

            for (unsigned int mcI{0}; mcI < size; ++mcI)
            {
                mcTracksMatrix[n][mcI] = MCtracks->at(mcI);
            }
        }

        for (int n = 0; n < nev; n++)
        {
            for (unsigned int mcI{0}; mcI < nTracks[n]; mcI++)
            {
                auto mcTrack = mcTracksMatrix[n][mcI];

                if (abs(mcTrack.GetPdgCode()) == hypPDG)
                {
                    counter_hyp++;
                    int firstDauID = mcTrack.getFirstDaughterTrackId();
                    int nDau = mcTrack.getLastDaughterTrackId();
                    bool hasTriton = false;
                    for (int iDau = firstDauID; iDau < nDau; iDau++)
                    {
                        auto dauTrack = mcTracksMatrix[n][iDau];
                        if (abs(dauTrack.GetPdgCode()) == tritonPDG)
                        {
                            hasTriton = true;
                            break;
                        }
                    }

                    if (!hasTriton)
                    {
                        counter_hyp_no_kink++;
                        LOG(info) << "hyp without triton";
                    }
                    
                        if(verbose) cout << "Particle PDGs in " <<n <<"x" <<mcI <<" event: " << endl;
                        for (int iDau = firstDauID; iDau < nDau; iDau++)
                        {
                            auto dauTrack = mcTracksMatrix[n][iDau];
                            
                                if (abs(dauTrack.GetPdgCode()) == pi0PDG)
                                    counter_pi0++;

                                if(verbose) cout << dauTrack.GetPdgCode() << endl;
                            
                        }
                    
                }
            }
        }
        LOG(info) << "TF " << tf << " has " << counter_hyp << " hyps and " << counter_hyp_no_kink << " without triton and " << counter_pi0 << " with pi0";
        LOG(info) << "---------------------------------------";
    }
}