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

#include "StrangenessTracking/StrangenessTracker.h"

#include "ReconstructionDataFormats/StrangeTrack.h"
#include "ReconstructionDataFormats/KinkTrack.h"

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

using StrangeTrack = o2::dataformats::StrangeTrack;
using KinkTrack = o2::dataformats::KinkTrack;

using MCLabCont = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
using CompClusterExt = o2::itsmft::CompClusterExt;

TString ptLabel = "#it{p}_{T} (GeV/#it{c})";
TString hypLabel = "M_{^{3}_{#Lambda}H} (GeV/c^{2})";
const int nBins = 100;
const int hypPDG = 1010010030;
const int tritonPDG = 1000010030;
const double tritonMass = 2.808921;
const double pi0Mass = 0.1349766;

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

void tree_test(TString path, TString filename, int tf_max = 80)
{
    TH1F *gen_mother_pt = new TH1F("Hyp Gen pt", "Hypertriton Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *gen_daughter_pt = new TH1F("Trit Gen pt", "Triton Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *gen_decay_lenght = new TH1F("Decay Lenght", "Decay Lenght (cm);cm;counts", nBins, 0, 50);
    TH1F *mother_pt = new TH1F("Hyp pt", "Hypertriton p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *daughter_pt = new TH1F("Trit pt", "Triton p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *decay_lenght = new TH1F("Gen Decay Lenght", "Decay Lenght (cm);cm;counts", nBins, 0, 50);
    TH1F *rec_r = new TH1F("Rec r", "Rec r (cm); r [cm];counts", nBins, 0, 50);

    TH1F *inv_mass = new TH1F("Invariant mass", "Invariant mass;" + hypLabel + ";counts", nBins, 2.9, 4);

    TH1F *r_res = new TH1F("r resolution", "r resolution; (r_{rec} - r_{gen})/r_{gen};counts", nBins, -5, 5);

    for (int tf = 1; tf <= tf_max; tf++)
    {
        LOG(info) << "Processing TF " << tf;

        TString tf_string = Form("%d", tf);
        TString tf_path = path + "tf" + tf_string;

        auto fITS = TFile::Open(tf_path + "/o2trac_its.root");
        auto fITSTPC = TFile::Open(tf_path + "/o2match_itstpc.root");
        auto fMCTracks = TFile::Open(tf_path + "/sgn_" + tf_string + "_Kine.root");
        auto fClusITS = TFile::Open(tf_path + "/o2clus_its.root");
        auto fStrangeTrack = TFile::Open(tf_path + "/o2_strange_tracks.root");

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");
        auto treeStrangeTrack = (TTree *)fStrangeTrack->Get("o2sim");

        // Tracks
        std::vector<MCTrack> *MCtracks = nullptr;
        std::vector<TrackITS> *ITStracks = nullptr;
        std::vector<TrackTPCITS> *ITSTPCtracks = nullptr;
        std::vector<StrangeTrack> *strangeTracks = nullptr;
        std::vector<KinkTrack> *kinkTracks = nullptr;

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
        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);
        treeStrangeTrack->SetBranchAddress("StrangeTracks", &strangeTracks);
        treeStrangeTrack->SetBranchAddress("KinkTracks", &kinkTracks);

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

                    int firstDauID = mcTrack.getFirstDaughterTrackId();
                    int nDau = mcTrack.getLastDaughterTrackId();
                    for (int iDau = firstDauID; iDau <= nDau; iDau++)
                    {
                        auto dauTrack = MCtracks->at(iDau);
                        if (abs(dauTrack.GetPdgCode()) == tritonPDG)
                        {
                            gen_mother_pt->Fill(mcTrack.GetPt());
                            gen_daughter_pt->Fill(dauTrack.GetPt());
                            gen_decay_lenght->Fill(dl);
                        }
                    }
                }
            }
        }

        // mc end

        for (int event = 0; event < treeStrangeTrack->GetEntriesFast(); event++)
        {
            if (!treeStrangeTrack->GetEvent(event) || !treeITS->GetEvent(event) || !treeITSTPC->GetEvent(event) || !treeMCTracks->GetEvent(event))
                continue;

            unsigned int size = kinkTracks->size();
            for (unsigned int mcI{0}; mcI < size; ++mcI)
            {
                auto kinkTrack = kinkTracks->at(mcI);
                if (kinkTrack.mITSRef == -1 || kinkTrack.mDecayRef == -1)
                    continue;
                auto mother = kinkTrack.mMother;
                auto daughter = kinkTrack.mDaughter;
                auto decayVtx = kinkTrack.mDecayVtx;
                double decayL = sqrt(decayVtx[0] * decayVtx[0] + decayVtx[1] * decayVtx[1] + decayVtx[2] * decayVtx[2]);
                double recR = sqrt(decayVtx[0] * decayVtx[0] + decayVtx[1] * decayVtx[1]);

                int ITSRef = kinkTrack.mITSRef;
                auto ITStrack = ITStracks->at(ITSRef);
                auto ITSLab = labITSvec->at(ITSRef);

                int trackID, evID, srcID;
                bool fake;
                ITSLab.get(trackID, evID, srcID, fake);

                int decayRef = kinkTrack.mDecayRef;
                auto daughterTrack = ITSTPCtracks->at(decayRef);
                auto daughterLab = labITSTPCvec->at(decayRef);

                int dauTrackID, dauEvID, dauSrcID;
                bool dauFake;
                daughterLab.get(dauTrackID, dauEvID, dauSrcID, dauFake);

                if (!fake && !dauFake)
                {
                    auto motherTrackMC = mcTracksMatrix[evID][trackID];
                    auto daughterTrackMC = mcTracksMatrix[dauEvID][dauTrackID];

                    int firstDauID = motherTrackMC.getFirstDaughterTrackId();
                    int nDau = motherTrackMC.getLastDaughterTrackId();
                    int recTritID = -10;
                    for (int iDau = firstDauID; iDau <= nDau; iDau++)
                    {
                        if (abs(mcTracksMatrix[evID][iDau].GetPdgCode()) == tritonPDG)
                        {
                            recTritID = iDau;
                            break;
                        }
                    }

                    if (abs(motherTrackMC.GetPdgCode()) == hypPDG && abs(daughterTrackMC.GetPdgCode()) == tritonPDG && recTritID == dauTrackID)
                    {
                        double genR = calcRadius(&mcTracksMatrix[evID], motherTrackMC, tritonPDG);
                        double resR = (recR - genR) / genR;
                        r_res->Fill(resR);
                        mother_pt->Fill(mother.getPt());
                        daughter_pt->Fill(daughter.getPt());
                        decay_lenght->Fill(decayL);
                        rec_r->Fill(recR);
                        double tritPabs = daughter.getP();
                        double hypPabs = mother.getP();
                        std::array<float, 3> hypP = {0, 0, 0};
                        std::array<float, 3> tritP = {0, 0, 0};
                        mother.getPxPyPzGlo(hypP);
                        daughter.getPxPyPzGlo(tritP);
                        float tritE = sqrt(tritPabs * tritPabs + tritonMass * tritonMass);
                        TVector3 piP = {hypP[0] - tritP[0], hypP[1] - tritP[1], hypP[2] - tritP[2]};
                        float piPabs = sqrt(piP[0] * piP[0] + piP[1] * piP[1] + piP[2] * piP[2]);
                        float piE = sqrt(pi0Mass * pi0Mass + piPabs * piPabs);
                        float hypE = piE + tritE;
                        double hypMass = sqrt(hypE * hypE - hypPabs * hypPabs);

                        inv_mass->Fill(hypMass);
                    }
                }
            }
        }
    }

    auto fFile = TFile::Open(filename, "RECREATE");
    gen_mother_pt->Write();
    gen_daughter_pt->Write();
    gen_decay_lenght->Write();

    mother_pt->Write();
    daughter_pt->Write();
    decay_lenght->Write();
    rec_r->Write();

    inv_mass->Write();

    r_res->Write();

    fFile->Close();
}
