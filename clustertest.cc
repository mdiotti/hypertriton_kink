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
#include "DataFormatsITSMFT/TrkClusRef.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsITSMFT/CompCluster.h"

#include "DetectorsVertexing/DCAFitterN.h"
#include "DataFormatsParameters/GRPObject.h"

#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"

#endif

#include <fstream>
#include <iostream>
using namespace std;

using namespace o2;
using namespace vertexing;

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using Vec3 = ROOT::Math::SVector<double, 3>;

using TrackITS = o2::its::TrackITS;
using MCLabCont = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
using CompClusterExt = o2::itsmft::CompClusterExt;

const int hypPDG = 1010010030;
const int tritonPDG = 1000010030;
const double tritonMass = 2.808921;
const double pi0Mass = 0.1349766;
const double hypMassTh = 2.99131;
TString chiLabel = "#chi^{2}";
TString hypLabel = "M_{^{3}_{#Lambda}H} (GeV/c^{2})";
int nBins = 100;
const double fontSize = 0.055;
const double markerSize = 4;

string FITTEROPTION = "DCA"; // "DCA_false" or "KFParticle"

double calcRadius(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        //        if (dauTrack.GetPdgCode() == dauPDG)
        if (abs(dauTrack.GetPdgCode()) == dauPDG)
        {
            auto radius = TMath::Sqrt((dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) * (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) + (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) * (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()));
            return radius;
        }
    }
    return -1;
}

void clustertest(TString path, TString filename, int tf_max = 40)
{

    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;
    int min_r = 0;
    int max_r = 50;

    TH1F *gen_r = new TH1F("Topology gen r", "Topology gen r;r (cm);counts", nBins, min_r, max_r);
    TH1F *rec_r = new TH1F("Topology rec r", "Topology rec r;r (cm);counts", nBins, min_r, max_r);
    TH1F *true_r = new TH1F("Topology true r", "Topology true r;r (cm);counts", nBins, min_r, max_r);
    TH1F *fit_r = new TH1F("Topology fit r", "Topology fit r;r (cm);counts", nBins, min_r, max_r);
    TH1F *true_fit_r = new TH1F("Topology true fit r", "Topology true fit r;r (cm);counts", nBins, min_r, max_r);

    TH1F *p_res_fake = new TH1F("p_res_fake", "p_res_fake;p_{rec} - p_{gen} (GeV/c);counts", nBins, -5, 5);
    TH1F *p_res_true = new TH1F("p_res", "p_res;p_{rec} - p_{gen} (GeV/c);counts", nBins, -5, 5);

    int clusterMother = 0;
    int clusterMotherTot = 0;

    for (int tf = tf_min; tf < tf_max; tf++)
    {
        TString tf_string = Form("%d", tf);
        TString tf_path = path + "tf" + tf_string;

        auto fITS = TFile::Open(tf_path + "/o2trac_its.root");
        auto fITSTPC = TFile::Open(tf_path + "/o2match_itstpc.root");
        auto fMCTracks = TFile::Open(tf_path + "/sgn_" + tf_string + "_Kine.root");
        auto fClusITS = TFile::Open(tf_path + "/o2clus_its.root");

        TString string_to_convert = tf_path + "/o2sim_grp.root";
        std::string path_string(string_to_convert.Data());
        const auto grp = o2::parameters::GRPObject::loadFrom(path_string);

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");

        // Tracks
        std::vector<MCTrack> *MCtracks = nullptr;
        std::vector<TrackITS> *ITStracks = nullptr;
        std::vector<o2::dataformats::TrackTPCITS> *ITSTPCtracks = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
        MCLabCont *clusLabArr = nullptr;

        // Clusters
        std::vector<CompClusterExt> *ITSclus = nullptr;
        std::vector<int> *ITSTrackClusIdx = nullptr;

        // Branches
        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeITSTPC->SetBranchAddress("TPCITS", &ITSTPCtracks);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);

        treeITSclus->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);
        treeITSclus->SetBranchAddress("ITSClusterComp", &ITSclus);
        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);

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
            if (!treeITS->GetEvent(event) || !treeITSTPC->GetEvent(event) || !treeITSclus->GetEvent(event))
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

                        int firstDauID = MCTrack.getFirstDaughterTrackId();
                        int nDau = MCTrack.getLastDaughterTrackId();
                        int tritID = 0;

                        for (int iDau = firstDauID; iDau <= nDau; iDau++)
                        {
                            if (mcTracksMatrix[evID][iDau].GetPdgCode() == tritonPDG)
                            {
                                tritID = iDau;
                                break;
                            }
                        }

                        if (tritID == 0)
                            continue; // if no triton daughter, improves speed
                        if (fake)
                            continue;

                        cout << "trackID: " << trackID << " evID: " << evID << endl;
                        cout << "track fake: " << fake << endl;
                        cout << "particle PDG " << MCTrack.GetPdgCode() << endl;

                        int firstClus = hypITSTrack.getFirstClusterEntry();
                        int nClus = hypITSTrack.getNumberOfClusters();
                        int lastClus = firstClus + nClus;

                        cout << "firstClus: " << firstClus << ", lastClus: " << lastClus << endl;
                        for (int iClus = firstClus; iClus < lastClus; iClus++)
                        {
                            cout << "\nprocessing cluster " << iClus << endl;
                            auto &labCls = (clusLabArr->getLabels(ITSTrackClusIdx->at(iClus)))[0];

                            int clustertrackID, clusterevID, clustersrcID;
                            bool clusterfake;
                            labCls.get(clustertrackID, clusterevID, clustersrcID, clusterfake);
                            cout << "clusterfake: " << clusterfake << endl;
                            cout << "clustertrackID: " << clustertrackID << ", clusterevID: " << clusterevID << endl;
                            auto clusterMCTrack = mcTracksMatrix[clusterevID][clustertrackID];
                            cout << "clusterMCTrack PDG " << clusterMCTrack.GetPdgCode() << endl;
                        }
                    }
                }
            }
        }
    }
}
