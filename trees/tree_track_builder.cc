#if !defined(CLING) || defined(ROOTCLING)

#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ITSMFTSimulation/Hit.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITStracking/IOUtils.h"
#include "DataFormatsParameters/GRPObject.h"

#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"

#include <gsl/gsl>
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TSystemDirectory.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "CommonDataFormat/RangeReference.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "StrangenessTracking/StrangenessTracker.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "StrangenessTracking/StrangenessTracker.h"
#include "GPUCommonArray.h"
#include "DetectorsBase/Propagator.h"

#endif

using namespace o2;
using namespace o2::framework;

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using Cascade = o2::dataformats::Cascade;

using MCTrack = o2::MCTrack;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using CompClusterExt = o2::itsmft::CompClusterExt;
using ITSCluster = o2::BaseCluster<float>;
using Vec3 = ROOT::Math::SVector<double, 3>;
using KinkTrack = o2::dataformats::KinkTrack;

const int motherPDG = 1010010030;        // hypertriton
const int firstDaughterPDG = 1000020030; // triton
const int secondDaughterPDG = 211;       // pi0

const float hypMass = 2.99131;
const float tritonMass = 2.80839;

// std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"ITS-TPC", labITSTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"TPC-TRD", labTPCTRDvec}, {"ITS-TPC-TRD", labITSTPCTRDvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}, {"TPC-TRD-TOF", labTPCTRDTOFvec}, {"ITS-TPC-TRD-TOF", labITSTPCTRDTOFvec}};

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

// KinkTrackTreeBuilder
void tree_track_builder(std::string path, std::string outSuffix = "")
{

    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 2) == "tf")
        {
            // if (stoi(file.substr(2)) > 1)
            if (stoi(file.substr(2)) > 80)
                continue;

            dirs.push_back(path + "/" + file);
            auto innerdir = (TSystemDirectory *)fileObj;
            auto innerfiles = innerdir->GetListOfFiles();
            for (auto innerfileObj : *innerfiles)
            {
                TString innerfile = ((TSystemFile *)innerfileObj)->GetName();
                if (innerfile.EndsWith("Kine.root") && innerfile.Contains("sgn"))
                {
                    kine_files.push_back(innerfile);
                }
            }
        }
    }

    TFile outFile = TFile(Form("TrackedKinkTreeTracks%s.root", outSuffix.data()), "recreate");
    TTree *outTree = new TTree("KinkTree", "KinkTree");

    KinkTrack gKinkTrack;
    outTree->Branch("kinkTrack", &gKinkTrack);
    // outTree->Branch("gX", &x);

    // create MC tree for efficiency calculation
    TTree *mcTree = new TTree("MCTree", "MCTree");
    /* MCTrack mcMother, mcDaughter;
       float mcMotherPt, mcDaughterPt, mcDL, mcCt;
       mcTree->Branch("mcMotherPt", &mcMotherPt);
       mcTree->Branch("mcDaughterPt", &mcDaughterPt);
       mcTree->Branch("mcDL", &mcDL);
       mcTree->Branch("mcCt", &mcCt);

       mcTree->Branch("mcMotherTrack", &mcMother);
       mcTree->Branch("mcDaughterTrack", &mcDaughter);
       */

    MCTrack gMCMother;
    MCTrack gMCDaughter;
    bool gIsDaughter;
    
    mcTree->Branch("mcMother", &gMCMother);
    mcTree->Branch("mcDaughter", &gMCDaughter);
    mcTree->Branch("isDaughter", &gIsDaughter);

    // Geometry
    o2::base::GeometryManager::loadGeometry(dirs[0] + "/o2sim_geometry.root");
    auto gman = o2::its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

    const auto grp = o2::parameters::GRPObject::loadFrom(dirs[0] + "/" + "o2sim_grp.root");
    o2::base::Propagator::initFieldFromGRP(grp);
    auto propagator = o2::base::Propagator::Instance();

    int counter = 0;
    for (unsigned int i = 0; i < dirs.size(); i++)
    {
        auto &dir = dirs[i];
        auto &kine_file = kine_files[i];
        LOG(info) << "Processing " << dir;
        // Files
        auto fMCTracks = TFile::Open((TString(dir + "/") + kine_file));
        auto fStrangeTracks = TFile::Open((dir + "/o2_strange_tracks.root").data());

        auto fITSTPC = TFile::Open((dir + "/o2match_itstpc.root").data());
        auto fTPCTOF = TFile::Open((dir + "/o2match_tof_tpc.root").data());
        auto fTPCTRD = TFile::Open((dir + "/trdmatches_tpc.root").data());
        auto fITSTPCTOF = TFile::Open((dir + "/o2match_tof_itstpc.root").data());
        auto fITS = TFile::Open((dir + "/o2trac_its.root").data());
        auto fClusITS = TFile::Open((dir + "/o2clus_its.root").data());
        auto fTPC = TFile::Open((dir + "/tpctracks.root").data());
        auto fITSTPCTRD = TFile::Open((dir + "/trdmatches_itstpc.root").data());
        auto fTPCTRDTOF = TFile::Open((dir + "/o2match_tof_tpctrd.root").data());
        auto fITSTPCTRDTOF = TFile::Open((dir + "/o2match_tof_itstpctrd.root").data());

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeStrangeTracks = (TTree *)fStrangeTracks->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSTPCTOF = (TTree *)fITSTPCTOF->Get("matchTOF");
        auto treeTPCTOF = (TTree *)fTPCTOF->Get("matchTOF");
        auto treeTPCTRD = (TTree *)fTPCTRD->Get("tracksTRD");
        auto treeITSTPCTRD = (TTree *)fITSTPCTRD->Get("tracksTRD");
        auto treeTPCTRDTOF = (TTree *)fTPCTRDTOF->Get("matchTOF");
        auto treeITSTPCTRDTOF = (TTree *)fITSTPCTRDTOF->Get("matchTOF");

        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");
        auto treeTPC = (TTree *)fTPC->Get("tpcrec");

        // MC Tracks
        std::vector<o2::MCTrack> *MCtracks = nullptr;
        std::vector<o2::itsmft::Hit> *ITSHits = nullptr;

        // Primary Vertex
        std::vector<o2::dataformats::PrimaryVertex> *primVertices = nullptr;

        // Hypertracks
        std::vector<KinkTrack> *kinkTrackVec = nullptr;
        std::vector<o2::track::TrackParametrizationWithError<float>> *hypertrackVec = nullptr;
        std::vector<o2::strangeness_tracking::ClusAttachments> *nAttachments = nullptr;

        // ITS tracks
        std::vector<o2::its::TrackITS> *ITStracks = nullptr;

        // Tracks
        std::vector<KinkTrack> *kinkTracks = nullptr;

        // Labels
        std::vector<o2::MCEventLabel> *pvMcArr = nullptr;

        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTRDvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTRDvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTRDTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTRDTOFvec = nullptr;

        // Clusters
        std::vector<CompClusterExt> *ITSclus = nullptr;
        o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
        std::vector<int> *ITSTrackClusIdx = nullptr;
        std::vector<unsigned char> *ITSpatt = nullptr;

        // Setting branches
        // treeStrangeTracks->SetBranchAddress("ClusUpdates", &nAttachments);
        treeStrangeTracks->SetBranchAddress("KinkTracks", &kinkTracks);
        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);

        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
        treeTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTOFvec);

        treeTPCTRD->SetBranchAddress("labels", &labTPCTRDvec);
        treeITSTPCTRD->SetBranchAddress("labelsTRD", &labITSTPCTRDvec);

        treeTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTRDTOFvec);
        treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);
        treeITSTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTRDTOFvec);

        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);
        treeITSclus->SetBranchAddress("ITSClusterComp", &ITSclus);
        treeITSclus->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);

        // define detector map
        std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"ITS-TPC", labITSTPCvec}, {"TPC", labTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"TPC-TRD", labTPCTRDvec}, {"ITS-TPC-TRD", labITSTPCTRDvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}, {"TPC-TRD-TOF", labTPCTRDTOFvec}, {"ITS-TPC-TRD-TOF", labITSTPCTRDTOFvec}};

        // fill MC matrix
        std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;

        auto nev = treeMCTracks->GetEntriesFast();
        mcTracksMatrix.resize(nev);
        for (int n = 0; n < nev; n++) // fill mcTracksMatrix
        {
            treeMCTracks->GetEvent(n);
            unsigned int size = MCtracks->size();
            mcTracksMatrix[n].resize(size);
            for (unsigned int mcI{0}; mcI < size; ++mcI)
            {
                auto mcTrack = MCtracks->at(mcI);
                mcTracksMatrix[n][mcI] = mcTrack;
            }
        }

        for (int event = 0; event < treeStrangeTracks->GetEntriesFast(); event++)
        {
            if (!treeStrangeTracks->GetEvent(event) || !treeMCTracks->GetEvent(event) || !treeITS->GetEvent(event) || !treeITSTPC->GetEvent(event) || !treeTPCTOF->GetEvent(event) || !treeTPCTRD->GetEvent(event) || !treeITSTPCTRD->GetEvent(event) || !treeITSTPCTOF->GetEvent(event) || !treeTPCTRDTOF->GetEvent(event) || !treeITSTPCTRDTOF->GetEvent(event) || !treeTPC->GetEvent(event))
                continue;

            unsigned int size = kinkTracks->size();
            for (unsigned int mcI{0}; mcI < size; ++mcI)
            {
                auto kinkTrack = kinkTracks->at(mcI);

                TString source = kinkTrack.mTrackIdx.getSourceName();

                int ITSRef = kinkTrack.mITSRef;
                int decayRef = kinkTrack.mTrackIdx.getIndex();

                // if(kinkTrack.mTrackIdx.getSourceName()=="TPC") continue; //exclude TPC only tracks (why?)
                // without TPC: 33k topology, 16k reconstructed
                // with TPC: 40k topology, 18k reconstructed

                auto ITSLab = labITSvec->at(ITSRef);
                auto labTrackType = map[kinkTrack.mTrackIdx.getSourceName()];
                auto decayLab = labTrackType->at(decayRef);

                if (!ITSLab.isValid())
                    continue;

                if (!decayLab.isValid())
                    continue;

                int trackID, evID, srcID;
                bool fake;
                ITSLab.get(trackID, evID, srcID, fake);

                int decayTrackID, decayEvID, decaySrcID;
                bool decayFake;
                decayLab.get(decayTrackID, decayEvID, decaySrcID, decayFake);

                auto motherTrackMC = mcTracksMatrix[evID][trackID];
                auto daughterTrackMC = mcTracksMatrix[decayEvID][decayTrackID];

                bool isHyp = false;
                if (abs(motherTrackMC.GetPdgCode()) == motherPDG)
                    isHyp = true;

                int firstDauID = motherTrackMC.getFirstDaughterTrackId();
                int nDau = motherTrackMC.getLastDaughterTrackId();
                int tritID = -10;
                if (isHyp)
                {
                    for (int iDau = firstDauID; iDau <= nDau; iDau++)
                    {
                        if (abs(mcTracksMatrix[evID][iDau].GetPdgCode()) == firstDaughterPDG)
                        {
                            tritID = iDau;
                            break;
                        }
                    }
                }

                bool isDaughter = false;
                if (tritID == decayTrackID && evID == decayEvID)
                    isDaughter = true;

                gMCMother = motherTrackMC;
                gMCDaughter = daughterTrackMC;
                gKinkTrack = kinkTrack;
                gIsDaughter = isDaughter;

                mcTree->Fill();
                outTree->Fill();
            }
        }

    } // end of files loop

    outFile.cd();
    outTree->Write();
    mcTree->Write();
}

std::array<int, 2> matchV0ToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, V0 &v0)
{
    std::array<int, 2> motherVec{-1, -1};
    std::array<std::array<int, 2>, 2> v0DauRefs;

    for (unsigned int iV0 = 0; iV0 < 2; iV0++)
    {
        v0DauRefs[iV0] = {-1, -1};
        if (map[v0.getProngID(iV0).getSourceName()])
        {
            auto labTrackType = map[v0.getProngID(iV0).getSourceName()];
            auto lab = labTrackType->at(v0.getProngID(iV0).getIndex());

            int trackID, evID, srcID;
            bool fake;
            lab.get(trackID, evID, srcID, fake);
            if (lab.isValid())
            {
                v0DauRefs[iV0] = {lab.getEventID(), lab.getTrackID()};
            }
        }
    }

    if (v0DauRefs[0][1] == -1 || v0DauRefs[1][1] == -1)
        return motherVec;

    auto &dau1MC = mcTracksMatrix[v0DauRefs[0][0]][v0DauRefs[0][1]];
    auto &dau2MC = mcTracksMatrix[v0DauRefs[1][0]][v0DauRefs[1][1]];

    if (!(std::abs(dau1MC.GetPdgCode()) == firstDaughterPDG && std::abs(dau2MC.GetPdgCode()) == secondDaughterPDG) && !(std::abs(dau1MC.GetPdgCode()) == secondDaughterPDG && std::abs(dau2MC.GetPdgCode()) == firstDaughterPDG))
        return motherVec;

    if (!dau1MC.isSecondary() || !dau2MC.isSecondary() || dau1MC.getMotherTrackId() != dau2MC.getMotherTrackId())
        return motherVec;
    auto v0MC = mcTracksMatrix[v0DauRefs[0][0]][dau1MC.getMotherTrackId()];

    if (std::abs(v0MC.GetPdgCode()) != motherPDG)
        return motherVec;

    motherVec = {v0DauRefs[0][0], dau1MC.getMotherTrackId()};
    return motherVec;
}
