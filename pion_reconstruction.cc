#if !defined(CLING) || defined(ROOTCLING)
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "DataFormatsITS/TrackITS.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"

#include "DetectorsVertexing/DCAFitterN.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "ITStracking/IOUtils.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"

#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITSBase/GeometryTGeo.h"

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

using namespace o2;
using namespace vertexing;

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using Vec3 = ROOT::Math::SVector<double, 3>;

using TrackITS = o2::its::TrackITS;
using ITSCluster = o2::BaseCluster<float>;
using TrackTPCITS = o2::dataformats::TrackTPCITS;

using MCLabCont = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
using CompClusterExt = o2::itsmft::CompClusterExt;

const int hypPDG = 1010010030;
const int tritonPDG = 1000010030;
const double tritonMass = 2.808921;
const double pi0Mass = 0.1349766;
const double hypMassTh = 2.99131;
TString ptLabel = "#it{p}_{T} (GeV/#it{c})";
TString resLabel = "Resolution: (gen-rec)/gen";
TString chiLabel = "#chi^{2}";
TString hypLabel = "M_{^{3}_{#Lambda}H} (GeV/c^{2})";

const int electronPDG = 11;
const int muonPDG = 13;
const int kaonPDG = 321;
const int pionPDG = 211;
const int protonPDG = 2212;

const double fontSize = 0.055;
const int nBins = 100;
string FITTEROPTION = "DCA"; // "DCA_false" or "KFParticle"

std::array<int, 2> matchCompLabelToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix,
                                      o2::MCCompLabel compLabel)
{
    std::array<int, 2> compRef = {-1, -1};
    int trackID, evID, srcID;
    bool fake;
    compLabel.get(trackID, evID, srcID, fake);
    if (compLabel.isValid())
    {
        compRef = {evID, trackID};
    }
    return compRef;
}

double calcRadius(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        if (abs(dauTrack.GetPdgCode()) == dauPDG)
        {
            auto radius = TMath::Sqrt((dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) * (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) + (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) * (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()));
            return radius;
        }
    }
    return -1;
}

double calcDecayLenght(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        if (abs(dauTrack.GetPdgCode()) == dauPDG)
        {
            auto dl = TMath::Sqrt((dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) * (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) + (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) * (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) + (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ()) * (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ()));
            return dl;
        }
    }
    return -1;
}

void pion_reconstruction(TString path, TString filename, bool propagation = false, int tf_max = 80)
{
    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;

    // Geometry
    TString string_to_convert = path + "tf1/o2sim_geometry.root";
    std::string path_string(string_to_convert.Data());
    o2::base::GeometryManager::loadGeometry(path_string);

    // GRP
    TString string_to_convert2 = path + "tf1/o2sim_grp.root";
    std::string path_string_grp(string_to_convert2.Data());
    const auto grp = o2::parameters::GRPObject::loadFrom(path_string_grp);

    // Load propagator
    o2::base::Propagator::initFieldFromGRP(grp);

    // Matching ITS tracks to MC tracks and V0
    std::array<int, 2> ITSref = {-1, 1};
    o2::its::TrackITS ITStrack;
    std::array<std::array<int, 2>, 7> clsRef;

    // Load Geometry
    auto gman = o2::its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

    for (int tf = tf_min; tf <= tf_max; tf++)
    {
        LOG(info) << "Processing TF " << tf;
        TString tf_string = Form("%d", tf);
        TString tf_path = path + "tf" + tf_string;

        auto fITS = TFile::Open(tf_path + "/o2trac_its.root");
        auto fITSTPC = TFile::Open(tf_path + "/o2match_itstpc.root");
        auto fMCTracks = TFile::Open(tf_path + "/sgn_" + tf_string + "_Kine.root");
        auto fClusITS = TFile::Open(tf_path + "/o2clus_its.root");

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");

        // Tracks
        std::vector<MCTrack> *MCtracks = nullptr;
        std::vector<TrackITS> *ITStracks = nullptr;
        std::vector<TrackTPCITS> *ITSTPCtracks = nullptr;

        // Clusters
        std::vector<CompClusterExt> *ITSclus = nullptr;
        o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
        std::vector<int> *ITSTrackClusIdx = nullptr;
        std::vector<unsigned char> *ITSpatt = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;

        // Branches
        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeITSTPC->SetBranchAddress("TPCITS", &ITSTPCtracks);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);

        treeITSclus->SetBranchAddress("ITSClusterComp", &ITSclus);
        treeITSclus->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);
        treeITSclus->SetBranchAddress("ITSClusterPatt", &ITSpatt);
        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);

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
                auto mcTrack = MCtracks->at(mcI);
                mcTracksMatrix[n][mcI] = mcTrack;
            }
        }
        // mc end

        for (int event = 0; event < treeITS->GetEntriesFast(); event++)
        {
            if (!treeITS->GetEvent(event) || !treeITSTPC->GetEvent(event) || !treeITSclus->GetEvent(event))
                continue;

            for (unsigned int iTrack{0}; iTrack < labITSvec->size(); ++iTrack)
            {
                double hypgenPabs = 0;
                TVector3 hypgenP;
                auto lab = labITSvec->at(iTrack);
                int trackID, evID, srcID;
                bool fake;
                lab.get(trackID, evID, srcID, fake);
                if (!lab.isNoise() && lab.isValid())
                {
                    auto mcTrack = mcTracksMatrix[evID][trackID];
                    if (abs(mcTrack.GetPdgCode()) == hypPDG)
                    {
                        // hypertriton histos fill
                        auto hypITSTrack = ITStracks->at(iTrack);

                        // topology histos fill

                        int firstDauID = mcTrack.getFirstDaughterTrackId();
                        int nDau = mcTrack.getLastDaughterTrackId();
                        int tritID = -10;
                        for (int iDau = firstDauID; iDau <= nDau; iDau++)
                        {
                            if (abs(mcTracksMatrix[evID][iDau].GetPdgCode()) == tritonPDG)
                            {
                                tritID = iDau;
                                break;
                            }
                        }

                        if (tritID == -10)
                            continue; // if no triton daughter, improves speed

                        int nLayers = hypITSTrack.getNumberOfClusters();
                        if (nLayers == 3) // 3 clusters recontruction doesn't work well
                            continue;

                        auto dautherTrack = mcTracksMatrix[evID][tritID];
                        for (unsigned int jTrack{0}; jTrack < labITSTPCvec->size(); ++jTrack)
                        {
                            bool isDaughter = false;
                            auto tritlab = labITSTPCvec->at(jTrack);
                            int trittrackID, tritevID, tritsrcID;
                            bool tritfake;
                            tritlab.get(trittrackID, tritevID, tritsrcID, tritfake);
                            if (!tritlab.isNoise() && tritlab.isValid())
                            {
                                auto tritmcTrack = mcTracksMatrix[tritevID][trittrackID];
                                if (abs(tritmcTrack.GetPdgCode()) == tritonPDG)
                                {
                                    int motherID = tritmcTrack.getMotherTrackId();
                                    auto motherTrack = mcTracksMatrix[evID][motherID];

                                    if (abs(motherTrack.GetPdgCode()) != hypPDG)
                                        continue;

                                    if (evID == tritevID && tritID == trittrackID)
                                        isDaughter = true;

                                    // if (!isDaughter)
                                    //     continue;

                                    auto tritITSTPCtrack = ITSTPCtracks->at(jTrack);
                                    float hypMass = 0;
                                    float hypMassProp = 0;
                                    double chi2 = 10000;
                                    double chi2Prop = 10000;
                                    bool hypTrue = false;
                                    bool tritTrue = false;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (unsigned int iTrack{0}; iTrack < labITSTPCvec->size(); ++iTrack) // triton histos fill
        {
            auto lab = labITSTPCvec->at(iTrack);
            int trackID, evID, srcID;
            bool fake;
            lab.get(trackID, evID, srcID, fake);
            if (!lab.isNoise() && lab.isValid())
            {
                auto mcTrack = mcTracksMatrix[evID][trackID];
                if (abs(mcTrack.GetPdgCode()) == tritonPDG)
                {
                    int motherID = mcTrack.getMotherTrackId();
                    auto motherTrack = mcTracksMatrix[evID][motherID];
                    if (abs(motherTrack.GetPdgCode()) != hypPDG)
                        continue;
                }
            }
        }

    } // tf loop end

    auto fFile = TFile(filename, "recreate");

    fFile.Close();
}
