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

void topology(TString path, TString filename, int tf_max = 80)
{
    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;

    // counts histogram of fake hyp
    TH1F *hCount = new TH1F("Counts", "Hypertriton Fake Track Belonging; 1: daughter triton, 2: non daughter triton, 3: other particle;counts", 3, 0.5, 3.5); // 1: daughter triton, 2: non daughter triton, 3: other particle

    // define hypertriton track histograms
    TH1F *hist_gen_pt = new TH1F("Hyp Gen pt", "Hypertriton Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_rec_pt = new TH1F("Hyp Rec pt", "Hypertriton Reconstructed p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_fake_pt = new TH1F("Hyp True pt", "Hypertriton True p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_ris_3 = new TH1F("Hyp Res 3 hit", "Hypertriton Resolution p_{T} 3 hit;" + resLabel + ";counts", 30, -1, 1);
    TH1F *hist_ris_4 = new TH1F("Hyp Res 4 hit", "Hypertriton Resolution p_{T} 4 hit;" + resLabel + ";counts", 30, -1, 1);
    TH1F *hist_gen_ct = new TH1F("Hyp Gen ct", "Hypertriton Generated c_{t};c_{t} (cm);counts", 50, 0, 50);
    TH1F *hist_rec_ct = new TH1F("Hyp Rec ct", "Hypertriton Reconstructed c_{t};c_{t} (cm);counts", 50, 0, 50);
    TH1F *hist_fake_ct = new TH1F("Hyp True ct", "Hypertriton True c_{t};c_{t} (cm);counts", 50, 0, 50);

    // define triton track histograms
    TH1F *hist_gen_pt_trit = new TH1F("Trit Gen pt", "Triton Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_rec_pt_trit = new TH1F("Trit Rec pt", "Triton Reconstructed p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_fake_pt_trit = new TH1F("Trit True pt", "Triton True p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_ris_trit = new TH1F("Trit Res", "Triton Resolution p_{T};" + resLabel + ";counts", 30, -0.25, 0.25);
    TH1F *hist_gen_ct_trit = new TH1F("Trit Gen r", "Triton Generated c_{t};c_{t} (cm);counts", 50, 0, 50);
    TH1F *hist_rec_ct_trit = new TH1F("Trit Rec r", "Triton Reconstructed c_{t};c_{t} (cm);counts", 50, 0, 50);
    TH1F *hist_fake_ct_trit = new TH1F("Trit True r", "Triton True c_{t};c_{t} (cm);counts", 50, 0, 50);

    // define topology histograms
    TH1F *hist_gen_pt_top = new TH1F("Top Gen pt", "Topology Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_rec_pt_top = new TH1F("Top Rec pt", "Topology Reconstructed p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_fake_pt_top = new TH1F("Top True pt", "Topology True p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *hist_ris_3_top = new TH1F("Top Res 3 hit", "Topology Resolution p_{T} 3 hit;" + resLabel + ";counts", 30, -1, 1);
    TH1F *hist_ris_4_top = new TH1F("Top Res 4 hit", "Topology Resolution p_{T} 4 hit;" + resLabel + ";counts", 30, -1, 1);
    TH1F *hist_gen_ct_top = new TH1F("Top Gen c_{t}", "Topology Generated c_{t};c_{t} (cm);counts", 50, 0, 50);
    TH1F *hist_rec_ct_top = new TH1F("Top Rec c_{t}", "Topology Reconstructed c_{t};c_{t} (cm);counts", 50, 0, 50);
    TH1F *hist_fake_ct_top = new TH1F("Top True c_{t}", "Topology True c_{t};c_{t} (cm);counts", 50, 0, 50);

    // define fit histograms
    TH1F *chi_squared = new TH1F("chi2", "" + chiLabel + ";" + chiLabel + ";counts", nBins, 0, 1);
    TH1F *chi_squared_fake = new TH1F("chi2_fake", "Fake " + chiLabel + ";" + chiLabel + ";counts", nBins, 0, 1);
    TH1F *inv_mass = new TH1F("Invariant mass", "Invariant mass;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *inv_mass_fake = new TH1F("Invariant mass fake", "Invariant mass fake;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *fit_ct = new TH1F("Topology fit ct", "Topology fit c_{t};c_{t} (cm);counts", 50, 0, 50);
    TH1F *fit_pt = new TH1F("Topology fit pT", "Topology fit p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *inv_mass_hight_layer = new TH1F("Invariant mass high layer", "Invariant mass first layer != 0;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *inv_mass_very_hight_layer = new TH1F("Invariant mass very high layer", "Invariant mass first layer >=2;" + hypLabel + ";counts", nBins, 2.9, 4);

    // define cluster histograms
    TH1F *n_cluster_fake = new TH1F("n_cluster_fake", "Number of fake clusters per track;N_{fake cluster};counts", 7, -0.5, 6.5);
    TH1F *fake_layer = new TH1F("fake_layer", "Layer of the fake cluster;layer;counts", 7, -0.5, 6.5);
    TH1F *nontriton_cluster = new TH1F("nontriton_cluster", "Type of non-triton fake clusters; 1:electron 2:muon 3:pion 4:kaon 5:proton 6:hypertriton ;counts", 6, 0.5, 6.5);

    int clusterMotherFake = 0;
    int allClusterMotherFake = 0;

    // Geometry
    TString string_to_convert = path + "tf1/o2sim_geometry.root";
    std::string path_string(string_to_convert.Data());
    o2::base::GeometryManager::loadGeometry(path_string);

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

        // GRP
        TString string_to_convert2 = tf_path + "/o2sim_grp.root";
        std::string path_string_grp(string_to_convert2.Data());
        const auto grp = o2::parameters::GRPObject::loadFrom(path_string_grp);

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");

        // Tracks
        std::vector<MCTrack> *MCtracks = nullptr;
        std::vector<TrackITS> *ITStracks = nullptr;
        std::vector<o2::dataformats::TrackTPCITS> *ITSTPCtracks = nullptr;

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
                    double dl = calcDecayLenght(MCtracks, mcTrack, tritonPDG);
                    // double dl = calcDecayLenght(&mcTracksMatrix[n], mcTrack, tritonPDG);
                    double ct = hypMassTh * dl / mcTrack.GetP();
                    hist_gen_pt->Fill(mcTrack.GetPt());
                    hist_gen_ct->Fill(ct);
                    int firstDauID = mcTrack.getFirstDaughterTrackId();
                    int nDau = mcTrack.getLastDaughterTrackId();
                    for (int iDau = firstDauID; iDau <= nDau; iDau++)
                    {
                        auto dauTrack = MCtracks->at(iDau);
                        // auto dauTrack = mcTracksMatrix[n][iDau];
                        if (abs(dauTrack.GetPdgCode()) == tritonPDG)
                        {
                            hist_gen_pt_top->Fill(mcTrack.GetPt());
                            hist_gen_ct_top->Fill(ct);
                            hist_gen_pt_trit->Fill(dauTrack.GetPt());
                            hist_gen_ct_trit->Fill(ct);
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

                        if (!fake)
                        {
                            hist_fake_pt->Fill(mcTrack.GetPt());
                            hist_fake_ct->Fill(ct);

                            if (hypITSTrack.getNumberOfClusters() >= 4)
                                hist_ris_4->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                            if (hypITSTrack.getNumberOfClusters() == 3)
                                hist_ris_3->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
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

                                    if (!isDaughter)
                                        continue;

                                    auto tritITSTPCtrack = ITSTPCtracks->at(jTrack);

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
                                                etaHyp = hypITSTrack.getEta();
                                                phiHyp = hypITSTrack.getPhi();

                                                tritTrackDCA.getPxPyPzGlo(tritP);
                                                tritPabs = tritTrackDCA.getP();
                                                etaTrit = tritITSTPCtrack.getEta();
                                                phiTrit = tritITSTPCtrack.getPhi();

                                                if (ft2.getChi2AtPCACandidate() < 0)
                                                    continue;

                                                double chi2 = ft2.getChi2AtPCACandidate();
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
                                                float hypMass = sqrt(hypE * hypE - hypPabs * hypPabs);
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

                                                if (recR < 18)
                                                    continue;
                                                if (std::abs(etaHyp - etaTrit) > 0.03 || std::abs(phiHyp - phiTrit) > 0.03)
                                                    continue;

                                                if (hypITSTrack.getFirstClusterLayer() != 0)
                                                inv_mass_hight_layer->Fill(hypMass);
                                                if (hypITSTrack.getFirstClusterLayer() >=2)
                                                inv_mass_very_hight_layer->Fill(hypMass);
                                                
                                                if (!tritfake && !fake)
                                                {
                                                    inv_mass->Fill(hypMass);
                                                }
                                                else
                                                    inv_mass_fake->Fill(hypMass);
                                                fit_ct->Fill(recCt);
                                                fit_pt->Fill(pt);
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

                        auto tritDl = calcDecayLenght(&mcTracksMatrix[evID], motherTrack, tritonPDG);
                        double tritCt = hypMassTh * tritDl / mcTrack.GetP();
                        hist_rec_ct_trit->Fill(tritCt);

                        if (!fake)
                        {
                            hist_fake_pt_trit->Fill(recPt);
                            hist_fake_ct_trit->Fill(tritCt);

                            hist_ris_trit->Fill((mcTrack.GetPt() - tritITSTPCTrack.getPt()) / mcTrack.GetPt());
                        }
                    }
                }
            }
        }
    } // tf loop end
    auto fFile = TFile(filename, "recreate");

    auto effDir = fFile.mkdir("efficiencies");
    auto fitDir = fFile.mkdir("fit");
    auto clusterDir = fFile.mkdir("cluster");

    // efficieny directory
    effDir->cd();
    hist_gen_pt->Write();
    hist_rec_pt->Write();
    hist_fake_pt->Write();
    hist_ris_3->Write();
    hist_ris_4->Write();
    hist_gen_ct->Write();
    hist_rec_ct->Write();
    hist_fake_ct->Write();

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

    TH1F *eff_pur = (TH1F *)hist_fake_ct->Clone("Hyp Eff Pur");
    eff_pur->GetYaxis()->SetTitle("Efficiency * Purity");
    eff_pur->SetTitle("Hypertriton Efficiency * Purity");
    eff_pur->Divide(hist_gen_ct);
    eff_pur->Write();

    hist_gen_pt_trit->Write();
    hist_rec_pt_trit->Write();
    hist_fake_pt_trit->Write();
    hist_ris_trit->Write();
    hist_gen_ct_trit->Write();
    hist_rec_ct_trit->Write();
    hist_fake_ct_trit->Write();

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

    hist_gen_pt_top->Write();
    hist_rec_pt_top->Write();
    hist_fake_pt_top->Write();
    hist_ris_3_top->Write();
    hist_ris_4_top->Write();
    hist_gen_ct_top->Write();
    hist_rec_ct_top->Write();
    hist_fake_ct_top->Write();

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

    TH1F *eff_pur_top = (TH1F *)hist_fake_ct_top->Clone("Top Eff Pur");
    eff_pur_top->GetYaxis()->SetTitle("Efficiency * Purity");
    eff_pur_top->SetTitle("Topopolgy Efficiency * Purity");
    eff_pur_top->Divide(hist_gen_ct_top);
    eff_pur_top->Write();

    // fit directory
    fitDir->cd();
    chi_squared->Write();
    chi_squared_fake->Write();
    inv_mass->Write();
    inv_mass_fake->Write();
    inv_mass_hight_layer->Write();
    inv_mass_very_hight_layer->Write();
    fit_ct->Write();

    TH1F *eff_fit_r = (TH1F *)fit_ct->Clone("Top Eff fit ct");
    eff_fit_r->GetXaxis()->SetTitleSize(fontSize);
    eff_fit_r->GetYaxis()->SetTitleSize(fontSize);
    eff_fit_r->GetYaxis()->SetTitle("Efficiency");
    eff_fit_r->SetTitle("Topology Fit ct Efficiency");
    eff_fit_r->Divide(hist_gen_ct);
    eff_fit_r->Write();

    fit_pt->Write();

    TH1F *eff_fit_pt = (TH1F *)fit_pt->Clone("Top Eff fit pt");
    eff_fit_pt->GetXaxis()->SetTitleSize(fontSize);
    eff_fit_pt->GetYaxis()->SetTitleSize(fontSize);
    eff_fit_pt->GetYaxis()->SetTitle("Efficiency");
    eff_fit_pt->SetTitle("Topology Fit p_{T} Efficiency");
    eff_fit_pt->Divide(hist_gen_pt);
    eff_fit_pt->Write();

    // cluster directory
    clusterDir->cd();
    hCount->Write();
    n_cluster_fake->Write();
    fake_layer->Write();
    nontriton_cluster->Write();

    fFile.Close();
}
