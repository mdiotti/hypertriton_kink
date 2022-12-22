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

int particle_output(int PDG)
{
    if (abs(PDG) == electronPDG)
        return 1;
    if (abs(PDG) == muonPDG)
        return 2;
    if (abs(PDG) == pionPDG)
        return 3;
    if (abs(PDG) == kaonPDG)
        return 4;
    if (abs(PDG) == protonPDG)
        return 5;
    if (abs(PDG) == hypPDG)
        return 6;
    if (abs(PDG) == tritonPDG)
        return 7;
    return 0;
}

int NFake(TrackITS track)
{
    int nFake{0};
    int firstClus = track.getFirstClusterLayer();
    for (int iCl = firstClus; iCl < track.getNClusters() + firstClus; ++iCl)
    {
        if (track.hasHitOnLayer(iCl) && track.isFakeOnLayer(iCl))
        {
            ++nFake;
        }
    }
    return nFake;
}

double_t getTrackClusChi2(o2::track::TrackParCov tritonTrack, const ITSCluster &clus)
{
    auto corrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrLUT;
    auto propInstance = o2::base::Propagator::Instance();
    auto gman = o2::its::GeometryTGeo::Instance();
    float alpha = gman->getSensorRefAlpha(clus.getSensorID()), x = clus.getX();
    int layer{gman->getLayer(clus.getSensorID())};

    if (!tritonTrack.rotate(alpha))
        return -1;

    if (!propInstance->propagateToX(tritonTrack, x, propInstance->getNominalBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, corrType))
        return -1;

    auto chi2 = tritonTrack.getPredictedChi2(clus);
    return chi2;
}

double_t getTrackClusChi2(o2::track::TrackParCov tritonTrack, const ITSCluster &clus, o2::track::TrackParCov &propagatedTrack)
{
    auto corrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrLUT;
    auto propInstance = o2::base::Propagator::Instance();
    auto gman = o2::its::GeometryTGeo::Instance();
    float alpha = gman->getSensorRefAlpha(clus.getSensorID()), x = clus.getX();
    int layer{gman->getLayer(clus.getSensorID())};

    if (!tritonTrack.rotate(alpha))
        return -1;

    if (!propInstance->propagateToX(tritonTrack, x, propInstance->getNominalBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, corrType))
        return -1;

    propagatedTrack = tritonTrack;
    auto chi2 = tritonTrack.getPredictedChi2(clus);
    return chi2;
}

bool updateTrack(const ITSCluster &clus, o2::track::TrackParCov &track, bool verbose = false)
{
    auto corrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrLUT;
    auto propInstance = o2::base::Propagator::Instance();
    auto gman = o2::its::GeometryTGeo::Instance();
    float alpha = gman->getSensorRefAlpha(clus.getSensorID()), x = clus.getX();
    int layer{gman->getLayer(clus.getSensorID())};

    if (verbose)
        LOG(info) << "Momentum before update: " << track.getP();
    if (!track.rotate(alpha))
        return false;

    if (!propInstance->propagateToX(track, x, propInstance->getNominalBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, corrType))
        return false;

    if (!track.update(clus))
        return false;

    if (verbose)
        LOG(info) << "Momentum after update: " << track.getP();

    return true;
}

bool createTrack(o2::track::TrackParCov &track, const ITSCluster *clus, int nClus, bool verbose = false)
{
    for (int i = 0; i < nClus; i++)
    {
        if (!updateTrack(clus[i], track, verbose))
        {
            return false;
        }
    }
    return true;
}

void topology(TString path, TString filename, bool propagation = false, int tf_max = 80)
{
    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;

    // counts histogram of fake hyp
    TH1F *hCount = new TH1F("Counts", "Hypertriton Fake Track Belonging; 1: daughter triton, 2: non daughter triton, 3: other particle;counts", 3, 0.5, 3.5); // 1: daughter triton, 2: non daughter triton, 3: other particle

    // define hypertriton track histograms
    TH1F *hist_gen_pt = new TH1F("Hyp Gen pt", "Hypertriton Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_rec_pt = new TH1F("Hyp Rec pt", "Hypertriton Reconstructed p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_fake_pt = new TH1F("Hyp True pt", "Hypertriton True p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_ris_3 = new TH1F("Hyp Res 3 hit", "Hypertriton Resolution p_{T} 3 hit;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *hist_ris_4 = new TH1F("Hyp Res 4 hit", "Hypertriton Resolution p_{T} 4 hit;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *res_fake_3 = new TH1F("Hyp Res Fake 3 hit", "Hypertriton Resolution p_{T} 3 hit;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *res_fake_4 = new TH1F("Hyp Res Fake 4 hit", "Hypertriton Resolution p_{T} 4 hit;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *hist_gen_ct = new TH1F("Hyp Gen ct", "Hypertriton Generated c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *hist_rec_ct = new TH1F("Hyp Rec ct", "Hypertriton Reconstructed c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *hist_fake_ct = new TH1F("Hyp True ct", "Hypertriton True c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *hist_gen_r = new TH1F("Hyp Gen r", "Hypertriton Generated r; r (cm);counts", nBins, 0, 50);
    TH1F *hist_rec_r = new TH1F("Hyp Rec r", "Hypertriton Reconstructed r; r (cm);counts", nBins, 0, 50);
    TH1F *hist_fake_r = new TH1F("Hyp True r", "Hypertriton True r; r (cm);counts", nBins, 0, 50);

    // define triton track histograms
    TH1F *hist_gen_pt_trit = new TH1F("Trit Gen pt", "Triton Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_rec_pt_trit = new TH1F("Trit Rec pt", "Triton Reconstructed p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_fake_pt_trit = new TH1F("Trit True pt", "Triton True p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_ris_trit = new TH1F("Trit Res", "Triton Resolution p_{T};" + resLabel + ";counts", nBins, -0.25, 0.25);
    TH1F *hist_gen_ct_trit = new TH1F("Trit Gen ct", "Triton Generated c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *hist_rec_ct_trit = new TH1F("Trit Rec ct", "Triton Reconstructed c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *hist_fake_ct_trit = new TH1F("Trit True ct", "Triton True c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *hist_gen_r_trit = new TH1F("Trit Gen r", "Triton Generated r; r (cm);counts", nBins, 0, 50);
    TH1F *hist_rec_r_trit = new TH1F("Trit Rec r", "Triton Reconstructed r; r (cm);counts", nBins, 0, 50);
    TH1F *hist_fake_r_trit = new TH1F("Trit True r", "Triton True r; r (cm);counts", nBins, 0, 50);

    // define topology histograms
    TH1F *hist_gen_pt_top = new TH1F("Top Gen pt", "Topology Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_rec_pt_top = new TH1F("Top Rec pt", "Topology Reconstructed p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_fake_pt_top = new TH1F("Top True pt", "Topology True p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_ris_3_top = new TH1F("Top Res 3 hit", "Topology Resolution p_{T} 3 hit;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *hist_ris_4_top = new TH1F("Top Res 4 hit", "Topology Resolution p_{T} 4 hit;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *hist_gen_ct_top = new TH1F("Top Gen ct", "Topology Generated c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *hist_rec_ct_top = new TH1F("Top Rec ct", "Topology Reconstructed c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *hist_fake_ct_top = new TH1F("Top True ct", "Topology True c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *hist_gen_r_top = new TH1F("Top Gen r", "Topology Generated r; r (cm);counts", nBins, 0, 50);
    TH1F *hist_rec_r_top = new TH1F("Top Rec r", "Topology Reconstructed r; r (cm);counts", nBins, 0, 50);
    TH1F *hist_fake_r_top = new TH1F("Top True r", "Topology True r; r (cm);counts", nBins, 0, 50);

    // define fit histograms
    TH1F *chi_squared = new TH1F("chi2", "" + chiLabel + ";" + chiLabel + ";counts", nBins, 0, 1);
    TH1F *chi_squared_fake = new TH1F("chi2_fake", "Fake " + chiLabel + ";" + chiLabel + ";counts", nBins, 0, 1);
    TH1F *inv_mass = new TH1F("Invariant mass", "Invariant mass;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *inv_mass_fake = new TH1F("Invariant mass fake", "Invariant mass fake;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *fit_ct = new TH1F("Topology fit ct", "Topology fit c_{t};c_{t} (cm);counts", nBins, 0, 50);
    TH1F *fit_r = new TH1F("Topology fit r", "Topology fit r; r (cm);counts", nBins, 0, 50);
    TH1F *fit_pt = new TH1F("Topology fit pT", "Topology fit p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *fit_rec_r = new TH1F("Topology fit rec r", "Topology fit rec r; r (cm);counts", nBins, 0, 50);
    TH1F *fit_r_res = new TH1F("Topology fit r res", "Fittted r resolution; (r_{rec} - r_{gen})/r_{gen};counts", nBins, -1, 1);

    // define cluster histograms
    TH1F *n_cluster_fake = new TH1F("n_cluster_fake", "Number of fake clusters per track;N_{fake cluster};counts", 7, -0.5, 6.5);
    TH1F *fake_layer = new TH1F("fake_layer", "Layer of the fake cluster;layer;counts", 7, -0.5, 6.5);
    TH1F *nontriton_cluster = new TH1F("nontriton_cluster", "Type of non-triton fake clusters; 1:electron 2:muon 3:pion 4:kaon 5:proton 6:hypertriton ;counts", 6, 0.5, 6.5);

    // define propagation histograms
    TH1F *chi_cluster_mother = new TH1F("chi_cluster_mother", "Chi2 of the cluster mother;#chi^{2};counts", nBins, 0, 10);
    TH1F *chi_cluster_daughter = new TH1F("chi_cluster_daughter", "Chi2 of the cluster daughter;#chi^{2};counts", nBins, 0, 10);
    TH1F *chi_cluster_other = new TH1F("chi_cluster_other", "Chi2 of the cluster other;#chi^{2};counts", nBins, 0, 10);
    /*TH1F *inv_mass_propagation = new TH1F("Invariant mass propagation", "Invariant mass propagation;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *chi2_propagation = new TH1F("chi2 propagation", "Chi2 propagation", nBins, 0, 1);
    TH1F *fit_ct_propagation = new TH1F("Topology fit ct propagation", "Topology fit c_{t} propagation;c_{t} (cm);counts", nBins, 0, 50);
    TH1F *fit_r_propagation = new TH1F("Topology fit r propagation", "Topology fit r propagation; r (cm);counts", nBins, 0, 50);
    TH1F *fit_pt_propagation = new TH1F("Topology fit pT propagation", "Topology fit p_{T} propagation;" + ptLabel + ";counts", nBins, 0, 12);*/
    TH1F *counts_propagation = new TH1F("Counts propagation", "Counts propagation; 1: Hyp non-prop 2: Hyp prop 3: Trit non-prop 4: Trit prop;counts", 4, 0.5, 4.5);
    TH1F *p_res_before = new TH1F("p_res_before", "p_{T} resolution before propagation;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *p_res_after = new TH1F("p_res_after", "p_{T} resolution after propagation;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *p_res_nonpropagated = new TH1F("p_res_nonpropagated", "p_{T} resolution non-propagated;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *p_res_5hit_before = new TH1F("p_res_5hit_before", "p_{T} resolution 5 hit before propagation;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *p_res_5hit_after = new TH1F("p_res_5hit_after", "p_{T} resolution 5 hit after propagation;" + resLabel + ";counts", nBins, -1, 1);
    TH1F *p_res_mother = new TH1F("p_res_mother", "p_{T} resolution when the cluster is mother's; resolution [GeV];counts", nBins, -1, 1);
    TH1F *p_res_daughter = new TH1F("p_res_daughter", "p_{T} resolution when the cluster is daughter's;resolution [GeV];counts", nBins, -1, 1);

    int nonPropagated = 0;
    int propagated = 0;

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
                if (abs(mcTrack.GetPdgCode()) == hypPDG)
                {
                    double rGen = calcRadius(MCtracks, mcTrack, tritonPDG);
                    double dl = calcDecayLenght(MCtracks, mcTrack, tritonPDG);
                    double ct = hypMassTh * dl / mcTrack.GetP();
                    hist_gen_pt->Fill(mcTrack.GetPt());
                    hist_gen_ct->Fill(ct);
                    hist_gen_r->Fill(rGen);
                    int firstDauID = mcTrack.getFirstDaughterTrackId();
                    int nDau = mcTrack.getLastDaughterTrackId();
                    for (int iDau = firstDauID; iDau <= nDau; iDau++)
                    {
                        auto dauTrack = MCtracks->at(iDau);
                        if (abs(dauTrack.GetPdgCode()) == tritonPDG)
                        {
                            hist_gen_pt_top->Fill(mcTrack.GetPt());
                            hist_gen_ct_top->Fill(ct);
                            hist_gen_r_top->Fill(rGen);
                            hist_gen_pt_trit->Fill(dauTrack.GetPt());
                            hist_gen_ct_trit->Fill(ct);
                            hist_gen_r_trit->Fill(rGen);
                        }
                    }
                }
            }
        }
        // mc end

        for (int event = 0; event < treeITS->GetEntriesFast(); event++)
        {
            if (!treeITS->GetEvent(event) || !treeITSTPC->GetEvent(event) || !treeITSclus->GetEvent(event))
                continue;

            o2::itsmft::TopologyDictionary mdict;
            std::vector<ITSCluster> ITSClusXYZ;
            if (propagation)
            {
                ITSClusXYZ.reserve(ITSclus->size());
                gsl::span<const unsigned char> spanPatt = *ITSpatt;
                auto pattIt = spanPatt.begin();
                mdict.readFromFile(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(o2::detectors::DetID::ITS, "utils/ITSdictionary.bin"));
                o2::its::ioutils::convertCompactClusters(*ITSclus, pattIt, ITSClusXYZ, &mdict);
            }

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
                        hist_rec_pt->Fill(mcTrack.GetPt());
                        auto rad = calcRadius(&mcTracksMatrix[evID], mcTrack, tritonPDG);
                        double dl = calcDecayLenght(&mcTracksMatrix[evID], mcTrack, tritonPDG);
                        double ct = hypMassTh * dl / mcTrack.GetP();
                        hist_rec_ct->Fill(ct);
                        hist_rec_r->Fill(rad);

                        if (!fake)
                        {
                            hist_fake_pt->Fill(mcTrack.GetPt());
                            hist_fake_ct->Fill(ct);
                            hist_fake_r->Fill(rad);

                            if (hypITSTrack.getNumberOfClusters() >= 4)
                                hist_ris_4->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                            if (hypITSTrack.getNumberOfClusters() == 3)
                                hist_ris_3->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                        }
                        else
                        {
                            if (hypITSTrack.getNumberOfClusters() >= 4)
                                res_fake_4->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                            if (hypITSTrack.getNumberOfClusters() == 3)
                                res_fake_3->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                        }

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

                                    hist_rec_pt_top->Fill(mcTrack.GetPt());
                                    hist_rec_ct_top->Fill(ct);
                                    hist_rec_r_top->Fill(rad);
                                    if (!tritfake)
                                    {
                                        hist_fake_pt_top->Fill(mcTrack.GetPt());
                                        hist_fake_ct_top->Fill(ct);
                                        hist_fake_r_top->Fill(rad);
                                    }

                                    if (!isDaughter)
                                        continue;

                                    auto tritITSTPCtrack = ITSTPCtracks->at(jTrack);
                                    o2::dataformats::TrackTPCITS newTritonTrack;
                                    bool hasPropagate = false;
                                    float hypMass = 0;
                                    float hypMassProp = 0;
                                    double chi2 = 10000;
                                    double chi2Prop = 10000;
                                    bool hypTrue = false;
                                    bool tritTrue = false;

                                    double chi2_clus = 1000;
                                    double chi2_clus_prop = 1000;
                                    o2::track::TrackParCov newHypTrack = hypITSTrack;
                                    newHypTrack.resetCovariance();
                                    o2::track::TrackParCov newTritTrack = tritITSTPCtrack;
                                    newTritTrack.resetCovariance();

                                    bool isBetterPropagated = false;

                                    if (propagation)
                                    {
                                        auto firstCls = hypITSTrack.getFirstClusterEntry(); // ultimo cluster della madre poichè in ordine inverso
                                        auto nCls = hypITSTrack.getNumberOfClusters();
                                        auto lastCls = firstCls + nCls - 1;
                                        auto &clusXYZ = ITSClusXYZ[(*ITSTrackClusIdx)[firstCls]];
                                        auto &clus = (*ITSclus)[(*ITSTrackClusIdx)[firstCls]];
                                        auto layer = gman->getLayer(clus.getSensorID());
                                        auto &labCls = (clusLabArr->getLabels(ITSTrackClusIdx->at(firstCls)))[0];
                                        clsRef[layer] = matchCompLabelToMC(mcTracksMatrix, labCls);

                                        if (clsRef[layer][0] <= -1 || clsRef[layer][1] <= -1)
                                            continue;

                                        auto MCTrack = mcTracksMatrix[clsRef[layer][0]][clsRef[layer][1]];
                                        int PDG = abs(MCTrack.GetPdgCode());

                                        if (PDG == hypPDG)
                                            hypTrue = true;
                                        else if (PDG == tritonPDG)
                                            tritTrue = true;

                                        double p_before = newTritTrack.getPt();
                                        updateTrack(clusXYZ, newTritTrack);
                                        double p_after = newTritTrack.getPt();
                                        if (hypTrue)
                                            p_res_mother->Fill((p_after - p_before));
                                        else if (tritTrue)
                                            p_res_daughter->Fill((p_after - p_before));

                                        for (int icl = lastCls; icl > firstCls; icl--) // creation of the new track
                                        {
                                            int i = -(icl - lastCls);
                                            auto &clusXYZToAdd = ITSClusXYZ[(*ITSTrackClusIdx)[icl]];
                                        }

                                        chi2_clus_prop = getTrackClusChi2(tritITSTPCtrack, clusXYZ, newTritonTrack);
                                        chi2_clus = getTrackClusChi2(newHypTrack, clusXYZ);

                                        if (chi2_clus_prop != -1)
                                            hasPropagate = true;

                                        if (hypTrue)
                                            chi_cluster_mother->Fill(chi2_clus_prop);
                                        else if (tritTrue)
                                            chi_cluster_daughter->Fill(chi2_clus_prop);
                                        else
                                            chi_cluster_other->Fill(chi2_clus_prop);

                                        // 1: Hyp best non-propagated 2: Hyp Best propagated 3: Trit best non-propagated 4: Trit best propagated
                                        if (!hasPropagate)
                                            counts_propagation->Fill(5);

                                        if (hasPropagate)
                                        {
                                            if (chi2_clus <= chi2_clus_prop)
                                            {
                                                if (hypTrue)
                                                    counts_propagation->Fill(1);
                                                else if (tritTrue)
                                                    counts_propagation->Fill(3);
                                            }
                                            else
                                            {
                                                isBetterPropagated = true;
                                                if (hypTrue)
                                                    counts_propagation->Fill(2);
                                                else if (tritTrue)
                                                    counts_propagation->Fill(4);
                                            }
                                        }

                                        if (isBetterPropagated)
                                        {
                                            p_res_before->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                                            p_res_after->Fill((mcTrack.GetPt() - newHypTrack.getPt()) / mcTrack.GetPt());
                                            if (nLayers >= 5)
                                            {
                                                p_res_5hit_before->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                                                p_res_5hit_after->Fill((mcTrack.GetPt() - newHypTrack.getPt()) / mcTrack.GetPt());
                                            }
                                        }
                                        else
                                            p_res_nonpropagated->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                                    }

                                    if (FITTEROPTION == "DCA")
                                    {
                                        try
                                        {
                                            auto bz = o2::base::Propagator::Instance()->getNominalBz();
                                            DCAFitter2 ft2;
                                            auto hypITSTrackParamOut = hypITSTrack.getParamOut();
                                            hypITSTrackParamOut.checkCovariance();
                                            tritITSTPCtrack.checkCovariance();
                                            // ft2.setMaxChi2(5);
                                            ft2.setMaxChi2(50);
                                            ft2.setUseAbsDCA(true);
                                            ft2.setBz(bz);
                                            ft2.process(hypITSTrackParamOut, tritITSTPCtrack);
                                            ft2.propagateTracksToVertex();
                                            if (ft2.isPropagateTracksToVertexDone() == true)
                                            {
                                                auto hypTrackDCA = ft2.getTrack(0);
                                                auto tritTrackDCA = ft2.getTrack(1);

                                                std::array<float, 3> hypP = {0, 0, 0};
                                                std::array<float, 3> tritP = {0, 0, 0};
                                                float hypPabs = 0;
                                                float tritPabs = 0;
                                                float etaHyp = 0;
                                                float phiHyp = 0;
                                                float etaTrit = 0;
                                                float phiTrit = 0;

                                                hypTrackDCA.getPxPyPzGlo(hypP);
                                                hypPabs = hypTrackDCA.getP();
                                                etaHyp = hypITSTrackParamOut.getEta();
                                                phiHyp = hypITSTrackParamOut.getPhi();

                                                tritTrackDCA.getPxPyPzGlo(tritP);
                                                tritPabs = tritTrackDCA.getP();
                                                etaTrit = tritITSTPCtrack.getEta();
                                                phiTrit = tritITSTPCtrack.getPhi();

                                                if (ft2.getChi2AtPCACandidate() < 0)
                                                    continue;

                                                chi2 = ft2.getChi2AtPCACandidate();
                                                if (!tritfake && !fake)
                                                    chi_squared->Fill(chi2);
                                                else
                                                    chi_squared_fake->Fill(chi2);

                                                std::array<float, 3> R = ft2.getPCACandidatePos();
                                                double recDl = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
                                                double recCt = hypMassTh * recDl / hypPabs;
                                                double recR = sqrt(R[0] * R[0] + R[1] * R[1]);

                                                float tritE = sqrt(tritPabs * tritPabs + tritonMass * tritonMass);
                                                std::array<float, 3> piP = {hypP[0] - tritP[0], hypP[1] - tritP[1], hypP[2] - tritP[2]};
                                                float piPabs = sqrt(piP[0] * piP[0] + piP[1] * piP[1] + piP[2] * piP[2]);
                                                float piE = sqrt(pi0Mass * pi0Mass + piPabs * piPabs);
                                                float hypE = piE + tritE;
                                                hypMass = sqrt(hypE * hypE - hypPabs * hypPabs);
                                                double pt = sqrt(hypP[0] * hypP[0] + hypP[1] * hypP[1]);

                                                if (fake)
                                                {
                                                    auto firstClus = hypITSTrack.getFirstClusterEntry();
                                                    auto ncl = hypITSTrack.getNumberOfClusters();

                                                    int nFake = NFake(hypITSTrack);
                                                    int nFakeFound = 0;
                                                    n_cluster_fake->Fill(nFake);
                                                    bool secondbin = false;
                                                    bool firstbin = false;
                                                    bool thirdbin = false;

                                                    for (int icl = 0; icl < ncl; icl++)
                                                    {
                                                        auto &labCls = (clusLabArr->getLabels(ITSTrackClusIdx->at(firstClus + icl)))[0];
                                                        auto &clus = (*ITSclus)[(*ITSTrackClusIdx)[firstClus + icl]];
                                                        auto layer = gman->getLayer(clus.getSensorID());
                                                        clsRef[layer] = matchCompLabelToMC(mcTracksMatrix, labCls);
                                                        if (clsRef[layer][0] > -1 && clsRef[layer][1] > -1)
                                                        {
                                                            auto MCTrack = mcTracksMatrix[clsRef[layer][0]][clsRef[layer][1]];
                                                            int PDG = MCTrack.GetPdgCode();

                                                            if (hypITSTrack.isFakeOnLayer(layer))
                                                            {
                                                                fake_layer->Fill(layer);

                                                                if (abs(PDG) == tritonPDG)
                                                                {
                                                                    if (clsRef[layer][0] == evID && clsRef[layer][1] == tritID)
                                                                    {
                                                                        nFakeFound++;
                                                                        if (nFakeFound == nFake)
                                                                        {
                                                                            firstbin = true;
                                                                        }
                                                                    }
                                                                    else
                                                                    {
                                                                        secondbin = true;
                                                                    }
                                                                }
                                                                else
                                                                {
                                                                    thirdbin = true;
                                                                    nontriton_cluster->Fill(particle_output(PDG));
                                                                }
                                                            }
                                                        }
                                                    }

                                                    if (firstbin)
                                                        hCount->Fill(1);
                                                    else if (secondbin)
                                                        hCount->Fill(2);
                                                    else if (thirdbin)
                                                        hCount->Fill(3);
                                                    else
                                                        hCount->Fill(4);
                                                } // end of fake

                                                if (recR < 18 || std::abs(etaHyp - etaTrit) > 0.03 || std::abs(phiHyp - phiTrit) > 0.03)
                                                    continue;

                                                if (!tritfake && !fake)
                                                {
                                                    inv_mass->Fill(hypMass);

                                                    fit_ct->Fill(recCt);
                                                    fit_pt->Fill(pt);
                                                    fit_r->Fill(rad);
                                                    fit_rec_r->Fill(recR);

                                                    double resR = (recR - rad) / rad;
                                                    fit_r_res->Fill(resR);
                                                }
                                                else
                                                    inv_mass_fake->Fill(hypMass);
                                            }
                                        }
                                        catch (std::runtime_error &e)
                                        {
                                            continue;
                                        }
                                    }
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

                    auto tritITSTPCTrack = ITSTPCtracks->at(iTrack);
                    double recPt = mcTrack.GetPt();
                    hist_rec_pt_trit->Fill(recPt);

                    auto tritR = calcRadius(&mcTracksMatrix[evID], motherTrack, tritonPDG);
                    auto tritDl = calcDecayLenght(&mcTracksMatrix[evID], motherTrack, tritonPDG);
                    double tritCt = hypMassTh * tritDl / mcTrack.GetP();
                    hist_rec_ct_trit->Fill(tritCt);
                    hist_rec_r_trit->Fill(tritR);

                    if (!fake)
                    {
                        hist_fake_pt_trit->Fill(recPt);
                        hist_fake_ct_trit->Fill(tritCt);
                        hist_fake_r_trit->Fill(tritR);

                        hist_ris_trit->Fill((mcTrack.GetPt() - tritITSTPCTrack.getPt()) / mcTrack.GetPt());
                    }
                }
            }
        }

    } // tf loop end

    /*
    Non propagated: 3258
    Propagated: 10659
    */

    auto fFile = TFile(filename, "recreate");

    auto effDir = fFile.mkdir("efficiencies");
    auto fitDir = fFile.mkdir("fit");
    auto clusterDir = fFile.mkdir("cluster");
    auto propDir = fFile.mkdir("propagation");

    // efficieny directory
    effDir->cd();
    hist_gen_pt->Write();
    hist_rec_pt->Write();
    hist_fake_pt->Write();
    hist_ris_3->Write();
    hist_ris_4->Write();
    res_fake_3->Write();
    res_fake_4->Write();
    hist_gen_ct->Write();
    hist_rec_ct->Write();
    hist_fake_ct->Write();
    hist_gen_r->Write();
    hist_rec_r->Write();
    hist_fake_r->Write();

    TH1F *eff_hist_pt = (TH1F *)hist_rec_pt->Clone("Hyp Eff pt");
    eff_hist_pt->GetYaxis()->SetTitle("Efficiency");
    eff_hist_pt->SetTitle("Hypertriton p_{T} Efficiency");
    eff_hist_pt->Divide(hist_gen_pt);
    eff_hist_pt->Write();

    TH1F *pur_hist_pt = (TH1F *)hist_fake_pt->Clone("Hyp Pur pt");
    pur_hist_pt->GetYaxis()->SetTitle("Purity");
    pur_hist_pt->SetTitle("Hypertriton p_{T} Purity");
    pur_hist_pt->Divide(hist_rec_pt);
    pur_hist_pt->Write();

    TH1F *eff_hist_r = (TH1F *)hist_rec_ct->Clone("Hyp Eff ct");
    eff_hist_r->GetYaxis()->SetTitle("Efficiency");
    eff_hist_r->SetTitle("Hypertriton ct Efficiency");
    eff_hist_r->Divide(hist_gen_ct);
    eff_hist_r->Write();

    TH1F *pur_hist_r = (TH1F *)hist_fake_ct->Clone("Hyp Pur ct");
    pur_hist_r->GetYaxis()->SetTitle("Purity");
    pur_hist_r->SetTitle("Hypertriton ct Purity");
    pur_hist_r->Divide(hist_rec_ct);
    pur_hist_r->Write();

    TH1F *eff_r = (TH1F *)hist_rec_r->Clone("Hyp Eff r");
    eff_r->GetYaxis()->SetTitle("Efficiency");
    eff_r->SetTitle("Hypertriton r Efficiency");
    eff_r->Divide(hist_gen_r);
    eff_r->Write();

    TH1F *pur_r = (TH1F *)hist_fake_r->Clone("Hyp Pur r");
    pur_r->GetYaxis()->SetTitle("Purity");
    pur_r->SetTitle("Hypertriton r Purity");
    pur_r->Divide(hist_rec_r);
    pur_r->Write();

    TH1F *eff_pur = (TH1F *)hist_fake_r->Clone("Hyp Eff Pur");
    eff_pur->GetYaxis()->SetTitle("Efficiency * Purity");
    eff_pur->SetTitle("Hypertriton Efficiency * Purity");
    eff_pur->Divide(hist_gen_r);
    eff_pur->Write();

    hist_gen_pt_trit->Write();
    hist_rec_pt_trit->Write();
    hist_fake_pt_trit->Write();
    hist_ris_trit->Write();
    hist_gen_ct_trit->Write();
    hist_rec_ct_trit->Write();
    hist_fake_ct_trit->Write();
    hist_gen_r_trit->Write();
    hist_rec_r_trit->Write();
    hist_fake_r_trit->Write();

    TH1F *eff_hist_pt_trit = (TH1F *)hist_rec_pt_trit->Clone("Trit Eff pt");
    eff_hist_pt_trit->GetYaxis()->SetTitle("Efficiency");
    eff_hist_pt_trit->SetTitle("Triton p_{T} Efficiency");
    eff_hist_pt_trit->Divide(hist_gen_pt_trit);
    eff_hist_pt_trit->Write();

    TH1F *pur_hist_pt_trit = (TH1F *)hist_fake_pt_trit->Clone("Trit Pur pt");
    pur_hist_pt_trit->GetYaxis()->SetTitle("Purity");
    pur_hist_pt_trit->SetTitle("Triton p_{T} Purity");
    pur_hist_pt_trit->Divide(hist_rec_pt_trit);
    pur_hist_pt_trit->Write();

    TH1F *eff_hist_r_trit = (TH1F *)hist_rec_ct_trit->Clone("Trit Eff ct");
    eff_hist_r_trit->GetYaxis()->SetTitle("Efficiency");
    eff_hist_r_trit->SetTitle("Triton ct Efficiency");
    eff_hist_r_trit->Divide(hist_gen_ct_trit);
    eff_hist_r_trit->Write();

    TH1F *pur_hist_r_trit = (TH1F *)hist_fake_ct_trit->Clone("Trit Pur ct");
    pur_hist_r_trit->GetYaxis()->SetTitle("Purity");
    pur_hist_r_trit->SetTitle("Triton ct Purity");
    pur_hist_r_trit->Divide(hist_rec_ct_trit);
    pur_hist_r_trit->Write();

    TH1F *eff_r_trit = (TH1F *)hist_rec_r_trit->Clone("Trit Eff r");
    eff_r_trit->GetYaxis()->SetTitle("Efficiency");
    eff_r_trit->SetTitle("Triton r Efficiency");
    eff_r_trit->Divide(hist_gen_r_trit);
    eff_r_trit->Write();

    TH1F *pur_r_trit = (TH1F *)hist_fake_r_trit->Clone("Trit Pur r");
    pur_r_trit->GetYaxis()->SetTitle("Purity");
    pur_r_trit->SetTitle("Triton r Purity");
    pur_r_trit->Divide(hist_rec_r_trit);
    pur_r_trit->Write();

    hist_gen_pt_top->Write();
    hist_rec_pt_top->Write();
    hist_fake_pt_top->Write();
    hist_ris_3_top->Write();
    hist_ris_4_top->Write();
    hist_gen_ct_top->Write();
    hist_rec_ct_top->Write();
    hist_fake_ct_top->Write();
    hist_gen_r_top->Write();
    hist_rec_r_top->Write();
    hist_fake_r_top->Write();

    TH1F *eff_hist_pt_top = (TH1F *)hist_rec_pt_top->Clone("Top Eff pt");
    eff_hist_pt_top->GetYaxis()->SetTitle("Efficiency");
    eff_hist_pt_top->SetTitle("Topopolgy p_{T} Efficiency");
    eff_hist_pt_top->Divide(hist_gen_pt_top);
    eff_hist_pt_top->Write();

    TH1F *pur_hist_pt_top = (TH1F *)hist_fake_pt_top->Clone("Top Pur pt");
    pur_hist_pt_top->GetYaxis()->SetTitle("Purity");
    pur_hist_pt_top->SetTitle("Topopolgy p_{T} Purity");
    pur_hist_pt_top->Divide(hist_rec_pt_top);
    pur_hist_pt_top->Write();

    TH1F *eff_hist_r_top = (TH1F *)hist_rec_ct_top->Clone("Top Eff ct");
    eff_hist_r_top->GetYaxis()->SetTitle("Efficiency");
    eff_hist_r_top->SetTitle("Topopolgy ct Efficiency");
    eff_hist_r_top->Divide(hist_gen_ct_top);
    eff_hist_r_top->Write();

    TH1F *pur_hist_r_top = (TH1F *)hist_fake_ct_top->Clone("Top Pur ct");
    pur_hist_r_top->GetYaxis()->SetTitle("Purity");
    pur_hist_r_top->SetTitle("Topopolgy ct Purity");
    pur_hist_r_top->Divide(hist_rec_ct_top);
    pur_hist_r_top->Write();

    TH1F *eff_r_top = (TH1F *)hist_ris_3_top->Clone("Top Eff r");
    eff_r_top->GetYaxis()->SetTitle("Efficiency");
    eff_r_top->SetTitle("Topopolgy r Efficiency");
    eff_r_top->Divide(hist_gen_ct_top);
    eff_r_top->Write();

    TH1F *pur_r_top = (TH1F *)hist_ris_4_top->Clone("Top Pur r");
    pur_r_top->GetYaxis()->SetTitle("Purity");
    pur_r_top->SetTitle("Topopolgy r Purity");
    pur_r_top->Divide(hist_rec_ct_top);
    pur_r_top->Write();

    TH1F *eff_pur_top = (TH1F *)hist_fake_r_top->Clone("Top Eff Pur");
    eff_pur_top->GetYaxis()->SetTitle("Efficiency * Purity");
    eff_pur_top->SetTitle("Topopolgy Efficiency * Purity");
    eff_pur_top->Divide(hist_gen_r_top);
    eff_pur_top->Write();

    // fit directory
    fitDir->cd();
    chi_squared->Write();
    chi_squared_fake->Write();
    inv_mass->Write();
    inv_mass_fake->Write();
    fit_ct->Write();

    TH1F *eff_fit_ct = (TH1F *)fit_ct->Clone("Top Eff fit ct");
    eff_fit_ct->GetXaxis()->SetTitleSize(fontSize);
    eff_fit_ct->GetYaxis()->SetTitleSize(fontSize);
    eff_fit_ct->GetYaxis()->SetTitle("Efficiency");
    eff_fit_ct->SetTitle("Topology Fit ct Efficiency");
    eff_fit_ct->Divide(hist_gen_ct);
    eff_fit_ct->Write();

    fit_pt->Write();

    TH1F *eff_fit_pt = (TH1F *)fit_pt->Clone("Top Eff fit pt");
    eff_fit_pt->GetXaxis()->SetTitleSize(fontSize);
    eff_fit_pt->GetYaxis()->SetTitleSize(fontSize);
    eff_fit_pt->GetYaxis()->SetTitle("Efficiency");
    eff_fit_pt->SetTitle("Topology Fit p_{T} Efficiency");
    eff_fit_pt->Divide(hist_gen_pt);
    eff_fit_pt->Write();

    fit_r->Write();

    TH1F *eff_fit_r = (TH1F *)fit_r->Clone("Top Eff fit r");
    eff_fit_r->GetXaxis()->SetTitleSize(fontSize);
    eff_fit_r->GetYaxis()->SetTitleSize(fontSize);
    eff_fit_r->GetYaxis()->SetTitle("Efficiency");
    eff_fit_r->SetTitle("Topology Fit r Efficiency");
    eff_fit_r->Divide(hist_gen_r);
    eff_fit_r->Write();

    fit_rec_r->Write();
    fit_r_res->Write();

    // cluster directory
    clusterDir->cd();
    hCount->Write();
    n_cluster_fake->Write();
    fake_layer->Write();
    nontriton_cluster->Write();

    // chi2 directory
    propDir->cd();
    chi_cluster_mother->Write();
    chi_cluster_daughter->Write();
    chi_cluster_other->Write();
    /*inv_mass_propagation->Write();
    chi2_propagation->Write();
    fit_ct_propagation->Write();
    fit_pt_propagation->Write();
    fit_r_propagation->Write()*/
    counts_propagation->Write();
    p_res_before->Write();
    p_res_after->Write();
    p_res_nonpropagated->Write();
    p_res_5hit_before->Write();
    p_res_5hit_after->Write();
    p_res_mother->Write();
    p_res_daughter->Write();

    fFile.Close();
}
