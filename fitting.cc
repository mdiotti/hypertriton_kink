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

#include "DetectorsVertexing/DCAFitterN.h"
#include "DataFormatsParameters/GRPObject.h"

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
TString chiLabel = "{chi}^2";
int nBins = 100;
double min_bins = 0;
double max_bins = 10;

string FITTEROPTION = "DCA"; // "DCA_false" or "KFParticle"

void fitting(TString path, TString filename, int tf_max = 40)
{
    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;

    TH1F *daughter_chi = new TH1F("Daughters chi2", "Daughters " + chiLabel + ";" + chiLabel + ";counts", nBins, min_bins, max_bins);
    TH1F *nondaughter_chi = new TH1F("Non-daughters chi2", "Non-daughters " + chiLabel + ";" + chiLabel + ";counts", nBins, min_bins, max_bins);

    for (int tf = tf_min; tf < tf_max; tf++)
    {
        TString tf_string = Form("%d", tf);
        TString tf_path = path + "tf" + tf_string;

        auto fITS = TFile::Open(tf_path + "/o2trac_its.root");
        auto fITSTPC = TFile::Open(tf_path + "/o2match_itstpc.root");
        auto fMCTracks = TFile::Open(tf_path + "/sgn_" + tf_string + "_Kine.root");

        TString string_to_convert =  tf_path + "/o2sim_grp.root";
        std::string path_string(string_to_convert.Data());
        const auto grp = o2::parameters::GRPObject::loadFrom(path_string);

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");

        // Tracks
        std::vector<MCTrack> *MCtracks = nullptr;
        std::vector<TrackITS> *ITStracks = nullptr;
        std::vector<o2::dataformats::TrackTPCITS> *ITSTPCtracks = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;

        // Branches
        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeITSTPC->SetBranchAddress("TPCITS", &ITSTPCtracks);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
        // std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}};

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

        for (int event = 0; event < treeITS->GetEntriesFast(); event++)
        {
            if (!treeITS->GetEvent(event) || !treeITSTPC->GetEvent(event))
                continue;

            for (unsigned int iTrack{0}; iTrack < labITSvec->size(); ++iTrack)
            {
                auto lab = labITSvec->at(iTrack);
                int trackID, evID, srcID;
                bool fake;
                lab.get(trackID, evID, srcID, fake);
                if (!lab.isNoise() && lab.isValid())
                {
                    auto MCTrack = mcTracksMatrix[evID][trackID];
                    if (abs(MCTrack.GetPdgCode()) == hypPDG)
                    {
                        auto hypITSTrack = ITStracks->at(iTrack);
                        if (hypITSTrack.getNumberOfClusters() == 3) // 3 clusters recontruction doesn't work well
                            continue;

                        int firstDauID = MCTrack.getFirstDaughterTrackId();
                        int nDau = MCTrack.getLastDaughterTrackId();
                        int tritID = 0;
                        for (int iDau = firstDauID; iDau < nDau; iDau++)
                        {
                            if (mcTracksMatrix[evID][iDau].GetPdgCode() == tritonPDG)
                            {
                                tritID = iDau;
                                break;
                            }
                        }

                        if (tritID == 0)
                            continue; // if no triton daughter, improves speed

                        for (unsigned int jTrack{0}; jTrack < labITSTPCvec->size(); ++jTrack)
                        {
                            bool isDaughter = false;

                            auto tritlab = labITSTPCvec->at(jTrack);
                            int trittrackID, tritevID, tritsrcID;
                            bool tritfake;
                            tritlab.get(trittrackID, tritevID, tritsrcID, tritfake);
                            if (!tritlab.isNoise() && tritlab.isValid())
                            {
                                auto tritMCTrack = mcTracksMatrix[tritevID][trittrackID];
                                if (abs(tritMCTrack.GetPdgCode()) == tritonPDG)
                                {
                                    if (trittrackID == tritID && tritevID == evID)
                                        isDaughter = true;

                                    auto tritITSTPCtrack = ITSTPCtracks->at(jTrack);
                                    if (!tritfake && !fake)
                                    {
                                        // Fitting start

                                        if (FITTEROPTION == "KFParticle")
                                        {
                                        }

                                        if (FITTEROPTION == "DCA")
                                        {
                                            try
                                            {
                                                DCAFitter2 ft2;
                                                hypITSTrack.checkCovariance();
                                                tritITSTPCtrack.checkCovariance();
                                                ft2.setMaxChi2(5);
                                                ft2.setBz(grp->getNominalL3Field());
                                                ft2.process(hypITSTrack, tritITSTPCtrack);
                                                ft2.propagateTracksToVertex();
                                            }
                                            catch (std::runtime_error &e)
                                            {
                                                continue;
                                            }
                                        }
                                    }
                                } // Fitting end
                            }
                        }
                    }
                }
            }
        } // end of event loop
    }     // end of tf loop

    auto fFile = TFile(filename, "recreate");

    fFile.Close();
} // end of fitting function
