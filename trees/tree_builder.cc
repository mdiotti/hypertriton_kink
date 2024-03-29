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
#include "StrangenessTracking/StrangenessTracker.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "StrangenessTracking/StrangenessTracker.h"
#include "GPUCommonArray.h"
#include "DetectorsBase/Propagator.h"

#include <algorithm>

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
const int firstDaughterPDG = 1000010030; // triton
const int secondDaughterPDG = 211;       // pi0

const float hypMass = 2.99131;
const float tritonMass = 2.80839;

// std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"ITS-TPC", labITSTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"TPC-TRD", labTPCTRDvec}, {"ITS-TPC-TRD", labITSTPCTRDvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}, {"TPC-TRD-TOF", labTPCTRDTOFvec}, {"ITS-TPC-TRD-TOF", labITSTPCTRDTOFvec}};

double calcDecayLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
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

double getPt(std::array<float, 3> momentum)
{
    return TMath::Sqrt(momentum[0] * momentum[0] + momentum[1] * momentum[1]);
}

void tree_builder(std::string path, std::string outSuffix = "", bool KFP = true)
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

    std::map<std::string, int> detectorMapN{{"ITS", 1}, {"ITS-TPC", 2}, {"TPC", 3}, {"TPC-TOF", 4}, {"TPC-TRD", 5}, {"ITS-TPC-TRD", 6}, {"ITS-TPC-TOF", 7}, {"TPC-TRD-TOF", 8}, {"ITS-TPC-TRD-TOF", 9}};
    //ITS,ITS-TPC,TPC-TOF,TPC-TRD,ITS-TPC-TRD,ITS-TPC-TOF,TPC-TRD-TOF,ITS-TPC-TRD-TOF

    TString name;
    if(KFP) name = "TrackedKFPKinkTree%s.root";
    else name = "TrackedKinkTree%s.root";
    TFile outFile = TFile(Form(name, outSuffix.data()), "recreate");

    TTree *outTree = new TTree("KinkTree", "KinkTree");
    float rMotherPt, rDaughterPt, rDecayLength, rAngle, Mass, rRadius, Chi2Match, Chi2DCA, gMotherPt, gDaughterPt, gDecayLength, gRadius;
    int NLayers, Detector;
    bool isTopology, isHyp, isTriton, motherFake, daughterFake, usingKFP;

    outTree->Branch("rMotherPt", &rMotherPt);
    outTree->Branch("rDaughterPt", &rDaughterPt);
    outTree->Branch("rDecayLength", &rDecayLength);
    outTree->Branch("rRadius", &rRadius);
    outTree->Branch("rAngle", &rAngle);
    outTree->Branch("Mass", &Mass);
    outTree->Branch("isTopology", &isTopology);
    outTree->Branch("isHyp", &isHyp);
    outTree->Branch("isTriton", &isTriton);
    outTree->Branch("Chi2Match", &Chi2Match);
    outTree->Branch("Chi2DCA", &Chi2DCA);
    outTree->Branch("NLayers", &NLayers);
    outTree->Branch("gMotherPt", &gMotherPt);
    outTree->Branch("gDaughterPt", &gDaughterPt);
    outTree->Branch("gDecayLength", &gDecayLength);
    outTree->Branch("gRadius", &gRadius);
    outTree->Branch("motherFake", &motherFake);
    outTree->Branch("daughterFake", &daughterFake);
    outTree->Branch("Detector", &Detector);
    outTree->Branch("usingKFP", &usingKFP);

    // create MC tree for efficiency calculation
    TTree *mcTree = new TTree("MCTree", "MCTree");
    float mcMotherPt, mcDaughterPt, mcDecayLength, mcRadius;

    mcTree->Branch("mcMotherPt", &mcMotherPt);
    mcTree->Branch("mcDaughterPt", &mcDaughterPt);
    mcTree->Branch("mcDecayLength", &mcDecayLength);
    mcTree->Branch("mcRadius", &mcRadius);

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
        // LOG(info) << "Processing " << dir;
        counter++;
        LOG(info) << "Processing TF" << counter;
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

                if (abs(mcTrack.GetPdgCode()) == motherPDG)
                {
                    int firstDauID = mcTrack.getFirstDaughterTrackId();
                    int nDau = mcTrack.getLastDaughterTrackId();
                    for (int iDau = firstDauID; iDau <= nDau; iDau++)
                    {
                        auto dauTrack = MCtracks->at(iDau);
                        if (abs(dauTrack.GetPdgCode()) == firstDaughterPDG)
                        {
                            float dl = calcDecayLength(MCtracks, mcTrack, firstDaughterPDG);
                            float r = calcRadius(MCtracks, mcTrack, firstDaughterPDG);

                            mcMotherPt = mcTrack.GetPt();
                            mcDaughterPt = dauTrack.GetPt();
                            mcDecayLength = dl;
                            mcRadius = r;
                            mcTree->Fill();
                        }
                    }
                }
            }
        }

        for (int event = 0; event < treeStrangeTracks->GetEntriesFast(); event++)
        {
            if (!treeStrangeTracks->GetEvent(event) || !treeMCTracks->GetEvent(event) || !treeITS->GetEvent(event) || !treeITSTPC->GetEvent(event) || !treeTPCTOF->GetEvent(event) || !treeTPCTRD->GetEvent(event) || !treeITSTPCTRD->GetEvent(event) || !treeITSTPCTOF->GetEvent(event) || !treeTPCTRDTOF->GetEvent(event) || !treeITSTPCTRDTOF->GetEvent(event) || !treeTPC->GetEvent(event))
                continue;

            unsigned int size = kinkTracks->size();
            for (unsigned int mcI{0}; mcI < size; ++mcI)
            {

                // setting default values
                rMotherPt = -1, rDaughterPt = -1, rDecayLength = -1, Mass = -1, rMotherPt = -1, rDaughterPt = -1, rAngle = -1, gRadius = -1, Chi2Match = -1, Chi2DCA = -1, gMotherPt = -1, gDaughterPt = -1;
                isTopology = false, isHyp = false, isTriton = false, motherFake = false, daughterFake = false;
                Detector = -1;

                auto kinkTrack = kinkTracks->at(mcI);
                if (kinkTrack.mITSRef == -1)
                    continue;

                if (kinkTrack.mKFused != KFP) {
                    LOG(info) << "KF tag mismatch!";
                    continue;
                }

                TString source = kinkTrack.mTrackIdx.getSourceName();

                int ITSRef = kinkTrack.mITSRef;
                int decayRef = kinkTrack.mTrackIdx.getIndex();
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

                if (abs(motherTrackMC.GetPdgCode()) == motherPDG)
                    isHyp = true;

                if (abs(daughterTrackMC.GetPdgCode()) == firstDaughterPDG)
                    isTriton = true;

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

                if (tritID == decayTrackID && evID == decayEvID)
                    isTopology = true;

                
                std::array<float, 3> hypP = kinkTrack.mMotherP;
                std::array<float, 3> tritP = kinkTrack.mDaughterP;
                float motherPabs = sqrt(hypP[0] * hypP[0] + hypP[1] * hypP[1] + hypP[2] * hypP[2]);
                float daughterPabs = sqrt(tritP[0] * tritP[0] + tritP[1] * tritP[1] + tritP[2] * tritP[2]);
                float scalarProduct = hypP[0] * tritP[0] + hypP[1] * tritP[1] + hypP[2] * tritP[2];
                float cosTheta = scalarProduct / (motherPabs * daughterPabs);
                float angle = acos(cosTheta) / TMath::Pi() * 180;

                std::array<float, 3> Vertex = {0, 0, 0};

                rMotherPt = getPt(hypP);
                rDaughterPt = getPt(tritP);
                Mass = kinkTrack.mMasses[0];
                Vertex = kinkTrack.mDecayVtx;
                rRadius = sqrt(Vertex[0] * Vertex[0] + Vertex[1] * Vertex[1]);
                rDecayLength = sqrt(Vertex[0] * Vertex[0] + Vertex[1] * Vertex[1] + Vertex[2] * Vertex[2]);
                Chi2Match = kinkTrack.mChi2Match;
                Chi2DCA = kinkTrack.mChi2Vertex;
                NLayers = kinkTrack.mNLayers;
                gMotherPt = motherTrackMC.GetPt();
                gDaughterPt = daughterTrackMC.GetPt();
                gRadius = calcRadius(&mcTracksMatrix[evID], motherTrackMC, firstDaughterPDG);
                gDecayLength = calcDecayLength(&mcTracksMatrix[evID], motherTrackMC, firstDaughterPDG);
                motherFake = fake;
                daughterFake = decayFake;
                Detector = detectorMapN[kinkTrack.mTrackIdx.getSourceName()];
                usingKFP = kinkTrack.mKFused;
                rAngle = angle;

                outTree->Fill();
            }
        }

    } // end of files loop

    outFile.cd();
    outTree->Write();
    mcTree->Write();
}
