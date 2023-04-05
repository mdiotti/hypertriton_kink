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
using ITSCluster = o2::BaseCluster<double>;
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

double calcDecayLenght(MCTrack *motherTrack, MCTrack *dauTrack)
{
    return TMath::Sqrt((dauTrack->GetStartVertexCoordinatesX() - motherTrack->GetStartVertexCoordinatesX()) * (dauTrack->GetStartVertexCoordinatesX() - motherTrack->GetStartVertexCoordinatesX()) + (dauTrack->GetStartVertexCoordinatesY() - motherTrack->GetStartVertexCoordinatesY()) * (dauTrack->GetStartVertexCoordinatesY() - motherTrack->GetStartVertexCoordinatesY()) + (dauTrack->GetStartVertexCoordinatesZ() - motherTrack->GetStartVertexCoordinatesZ()) * (dauTrack->GetStartVertexCoordinatesZ() - motherTrack->GetStartVertexCoordinatesZ()));
}

double calcRadius(MCTrack *motherTrack, MCTrack *dauTrack)
{
    return TMath::Sqrt((dauTrack->GetStartVertexCoordinatesX() - motherTrack->GetStartVertexCoordinatesX()) * (dauTrack->GetStartVertexCoordinatesX() - motherTrack->GetStartVertexCoordinatesX()) + (dauTrack->GetStartVertexCoordinatesY() - motherTrack->GetStartVertexCoordinatesY()) * (dauTrack->GetStartVertexCoordinatesY() - motherTrack->GetStartVertexCoordinatesY()));
}

// ITS,ITS-TPC,TPC-TOF,TPC-TRD,ITS-TPC-TRD,ITS-TPC-TOF,TPC-TRD-TOF,ITS-TPC-TRD-TOF

void read_tree(TString filename)
{
    TH1F *gen_mother_pt = new TH1F("Hyp Gen pt", "Hypertriton Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *gen_daughter_pt = new TH1F("Trit Gen pt", "Triton Generated p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *gen_decay_lenght = new TH1F("Gen Decay Lenght", "Decay Lenght (cm);cm;counts", nBins, 0, 50);
    TH1F *gen_r = new TH1F("Gen r", "Gen r (cm); r [cm];counts", nBins, 0, 50);
    TH1F *mother_pt = new TH1F("Hyp Rec pt", "Hypertriton p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *daughter_pt = new TH1F("Trit Red pt", "Triton p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *decay_lenght = new TH1F("Red Decay Lenght", "Decay Lenght (cm);cm;counts", nBins, 0, 50);
    TH1F *rec_r = new TH1F("Rec r", "Rec r (cm); r [cm];counts", nBins, 0, 50);

    TH1F *inv_mass = new TH1F("Invariant mass", "Invariant mass;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *inv_mass_MC = new TH1F("Invariant mass MC", "Invariant mass MC;" + hypLabel + ";counts", nBins, 2.9, 4);

    TH1F *r_res = new TH1F("r resolution", "r resolution; (r_{rec} - r_{gen})/r_{gen};counts", nBins, -5, 5);
    TH1F *r_res_MC = new TH1F("r resolution MC", "r resolution MC; (r_{rec} - r_{gen})/r_{gen};counts", nBins, -5, 5);

    TH1F *m_differences = new TH1F("m differences", "m differences; m_{pxyz} - m_{|p|};counts", nBins, -0.1, 0.1);
    TH1F *m_differences_MC = new TH1F("m differences MC", "m differences MC; m_{pxyz} - m_{|p|};counts", nBins, -0.1, 0.1);

    TH1F *clus_size = new TH1F("Cluster size", "Cluster size; Cluster size;counts", nBins, 0, 10);
    TH1F *clus_size_MC = new TH1F("Cluster size MC", "Cluster size MC; Cluster size;counts", nBins, 0, 10);
    TH1F *chi_all = new TH1F("Chi2 all", "Chi2 all; Chi2;counts", nBins, 0, 50);

    TH1F *chi2 = new TH1F("Chi2", "Chi2; Chi2;counts", nBins, 0, 1);
    TH1F *chi2_MC = new TH1F("Chi2 MC", "Chi2 MC; Chi2;counts", nBins, 0, 1);

    TH1F *true_mother_pt = new TH1F("True Hyp pt", "True Hypertriton p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *true_daughter_pt = new TH1F("True Trit pt", "True Triton p_{T};" + ptLabel + ";counts", nBins, 0, 12);
    TH1F *true_decay_lenght = new TH1F("True Decay Lenght", "True Decay Lenght (cm);cm;counts", nBins, 0, 50);
    TH1F *true_r = new TH1F("True r", "True r (cm); r [cm];counts", nBins, 0, 50);

    TH1F *lose_point = new TH1F("Lose point", "Lose point; Lose point;counts", 3, 0.5, 3.5);
    TH1F *lose_point_final = new TH1F("Lose point final", "Lose point final; Lose point;counts", 3, 0.5, 3.5);

    // Teoric:  ITS,ITS-TPC,TPC-TOF,TPC-TRD,ITS-TPC-TRD,ITS-TPC-TOF,TPC-TRD-TOF,ITS-TPC-TRD-TOF
    // Found:   ITS-TPC,TPC-TOF,TPC-TRD,ITS-TPC-TRD,ITS-TPC-TOF

    // Files
    auto fTrees = TFile::Open("TrackedKinkTreeTracks.root");

    // Trees
    auto treeMC = (TTree *)fTrees->Get("MCTree");
    auto treeKinkTracks = (TTree *)fTrees->Get("KinkTree");

    // Tracks
    MCTrack *MCMotherTracks = nullptr;
    MCTrack *MCDaughterTracks = nullptr;
    KinkTrack *kinkTracks = nullptr;

    Bool_t *isDaughterVec = nullptr;

    // Branches
    treeMC->SetBranchAddress("mcMother", &MCMotherTracks);
    treeMC->SetBranchAddress("mcDaughter", &MCDaughterTracks);
    treeMC->SetBranchAddress("isDaughter", &isDaughterVec);
    treeKinkTracks->SetBranchAddress("kinkTrack", &kinkTracks);

    for (int event = 0; event < treeKinkTracks->GetEntriesFast(); event++)
    {
        if (!treeKinkTracks->GetEvent(event) || !treeMC->GetEvent(event))
            continue;

        auto motherTrackMC = MCMotherTracks;
        auto daughterTrackMC = MCDaughterTracks;
        auto kinkTrack = kinkTracks;
        auto isDaughter = isDaughterVec;
        if(isDaughter) LOG(info) <<isDaughter;

        double genR = calcRadius(motherTrackMC, daughterTrackMC);
        double genL = calcDecayLenght(motherTrackMC, daughterTrackMC);

        gen_mother_pt->Fill(motherTrackMC->GetPt());
        gen_daughter_pt->Fill(daughterTrackMC->GetPt());
        gen_decay_lenght->Fill(genL);
        gen_r->Fill(genR);

        // if(kinkTrack.mTrackIdx.getSourceName()=="TPC") continue;
        // without TPC: 33k topology, 16k reconstructed
        // with TPC: 40k topology, 18k reconstructed

        bool isHyp = false;
        if (abs(motherTrackMC->GetPdgCode()) == hypPDG)
            isHyp = true;

        auto mother = kinkTrack->mMother;
        auto daughter = kinkTrack->mDaughter;
        auto decayVtx = kinkTrack->mDecayVtx;
        double decayL = sqrt(decayVtx[0] * decayVtx[0] + decayVtx[1] * decayVtx[1] + decayVtx[2] * decayVtx[2]);
        double recR = sqrt(decayVtx[0] * decayVtx[0] + decayVtx[1] * decayVtx[1]);

        int nLayers = kinkTrack->mNLayers;
        if (nLayers == 3) // 3 clusters recontruction doesn't work well
            continue;

        float etaHyp = 0;
        float phiHyp = 0;
        float etaTrit = 0;
        float phiTrit = 0;

        etaHyp = mother.getEta();
        phiHyp = mother.getPhi();

        etaTrit = daughter.getEta();
        phiTrit = daughter.getPhi();

        if (recR < 18 || std::abs(etaHyp - etaTrit) > 0.03 || std::abs(phiHyp - phiTrit) > 0.03)
            continue;

        double resR = (recR - genR) / genR;
        r_res->Fill(resR);

        mother_pt->Fill(mother.getPt());
        daughter_pt->Fill(daughter.getPt());
        decay_lenght->Fill(decayL);
        rec_r->Fill(recR);

        double tritPabs = daughter.getP();
        double hypPabs = mother.getP();
        double tritE = sqrt(tritPabs * tritPabs + tritonMass * tritonMass);

        std::array<float, 3> hypP = {0, 0, 0};
        std::array<float, 3> tritP = {0, 0, 0};
        mother.getPxPyPzGlo(hypP);
        daughter.getPxPyPzGlo(tritP);
        TVector3 piP = {hypP[0] - tritP[0], hypP[1] - tritP[1], hypP[2] - tritP[2]};
        double piPabs = sqrt(piP[0] * piP[0] + piP[1] * piP[1] + piP[2] * piP[2]);
        double piE = sqrt(pi0Mass * pi0Mass + piPabs * piPabs);
        double hypE = piE + tritE;
        double hypMass = sqrt(hypE * hypE - hypPabs * hypPabs);

        double hypMassKink = kinkTrack->mMasses[0];

        double mDiff = hypMass - hypMassKink;

        if (hypE < tritE)
            continue;

        float clusSize = kinkTrack->mITSClusSize;
        float chi2k = kinkTrack->mChi2Match;

        chi_all->Fill(chi2k);

        if (isHyp && isDaughter)
        {
            r_res_MC->Fill(resR);
            inv_mass_MC->Fill(hypMass);
            m_differences_MC->Fill(mDiff);
            clus_size_MC->Fill(clusSize);
            chi2_MC->Fill(chi2k);
            true_mother_pt->Fill(hypPabs);
            true_daughter_pt->Fill(tritPabs);
            true_decay_lenght->Fill(decayL);
            true_r->Fill(recR);
        }

        inv_mass->Fill(hypMass);
        m_differences->Fill(mDiff);
        clus_size->Fill(clusSize);
        chi2->Fill(chi2k);
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
    inv_mass_MC->Write();

    r_res->Write();
    r_res_MC->Write();

    m_differences->Write();
    m_differences_MC->Write();

    clus_size->Write();
    clus_size_MC->Write();

    chi2->Write();
    chi2_MC->Write();
    chi_all->Write();

    true_mother_pt->Write();
    true_daughter_pt->Write();
    true_decay_lenght->Write();
    true_r->Write();

    lose_point->Write();
    lose_point_final->Write();

    fFile->Close();
}
