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
TString ptLabel = "#it{p}_{T} (GeV/#it{c})";
TString resLabel = "Resolution: (gen-rec)/gen";

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

void efficiency_functions(TString path, TString filename, int tf_max = 40)
{
    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;

    // define hypertriton track histograms
    TH1F *hist_gen_pt = new TH1F("Hyp Gen pt", "Hypertriton Generated p_{T};" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_rec_pt = new TH1F("Hyp Rec pt", "Hypertriton Reconstructed p_{T};" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_fake_pt = new TH1F("Hyp True pt", "Hypertriton True p_{T};" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_ris_3 = new TH1F("Hyp Res 3 hit", "Hypertriton Resolution p_{T} 3 hit;" + resLabel + ";counts", 30, -1, 1);
    TH1F *hist_ris_4 = new TH1F("Hyp Res 4 hit", "Hypertriton Resolution p_{T} 4 hit;" + resLabel + ";counts", 30, -1, 1);
    TH1F *hist_gen_r = new TH1F("Hyp Gen r", "Hypertriton Generated Radius;Radius (cm);counts", 50, 0, 50);
    TH1F *hist_rec_r = new TH1F("Hyp Rec r", "Hypertriton Reconstructed Radius;Radius (cm);counts", 50, 0, 50);
    TH1F *hist_fake_r = new TH1F("Hyp True r", "Hypertriton True Radius;Radius (cm);counts", 50, 0, 50);

    // define triton track histograms
    TH1F *hist_gen_pt_trit = new TH1F("Trit Gen pt", "Triton Generated p_{T};" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_rec_pt_trit = new TH1F("Trit Rec pt", "Triton Reconstructed p_{T};" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_fake_pt_trit = new TH1F("Trit True pt", "Triton True p_{T};" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_ris_trit = new TH1F("Trit Res", "Triton Resolution p_{T} 3 hit;" + resLabel + ";counts", 30, -0.25, 0.25);
    TH1F *hist_gen_r_trit = new TH1F("Trit Gen r", "Triton Generated Radius;Radius (cm);counts", 50, 0, 50);
    TH1F *hist_rec_r_trit = new TH1F("Trit Rec r", "Triton Reconstructed Radius;Radius (cm);counts", 50, 0, 50);
    TH1F *hist_fake_r_trit = new TH1F("Trit True r", "Triton True Radius;Radius (cm);counts", 50, 0, 50);

    // define topology histograms
    TH1F *hist_gen_pt_top = new TH1F("Top Gen pt", "Topology Generated p_{T};" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_rec_pt_top = new TH1F("Top Rec pt", "Topology Reconstructed p_{T};" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_fake_pt_top = new TH1F("Top True pt", "Topology True p_{T};" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_ris_3_top = new TH1F("Top Res 3 hit", "Topology Resolution p_{T} 3 hit;" + resLabel + ";counts", 30, -1, 1);
    TH1F *hist_ris_4_top = new TH1F("Top Res 4 hit", "Topology Resolution p_{T} 4 hit;" + resLabel + ";counts", 30, -1, 1);
    TH1F *hist_gen_r_top = new TH1F("Top Gen r", "Topology Generated Radius;Radius (cm);counts", 50, 0, 50);
    TH1F *hist_rec_r_top = new TH1F("Top Rec r", "Topology Reconstructed Radius;Radius (cm);counts", 50, 0, 50);
    TH1F *hist_fake_r_top = new TH1F("Top True r", "Topology True Radius;Radius (cm);counts", 50, 0, 50);

    for (int tf = tf_min; tf < tf_max; tf++)
    {
        TString tf_string = Form("%d", tf);
        TString tf_path = path + "tf" + tf_string;

        auto fITS = TFile::Open(tf_path + "/o2trac_its.root");
        auto fITSTPC = TFile::Open(tf_path + "/o2match_itstpc.root");
        auto fMCTracks = TFile::Open(tf_path + "/sgn_" + tf_string + "_Kine.root");

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

        for (int n = 0; n < nev; n++) // fill histos
        {
            for (unsigned int mcI{0}; mcI < nTracks[n]; mcI++)
            {
                auto mcTrack = mcTracksMatrix[n][mcI];

                if (abs(mcTrack.GetPdgCode()) == hypPDG)
                {
                    double radius = calcRadius(&mcTracksMatrix[n], mcTrack, tritonPDG);
                    hist_gen_pt->Fill(mcTrack.GetPt());
                    hist_gen_r->Fill(radius);
                    int firstDauID = mcTrack.getFirstDaughterTrackId();
                    int nDau = mcTrack.getLastDaughterTrackId();
                    for (int iDau = firstDauID; iDau < nDau; iDau++)
                    {
                        auto dauTrack = mcTracksMatrix[n][iDau];
                        if (abs(dauTrack.GetPdgCode()) == tritonPDG)
                        {
                            hist_gen_pt_top->Fill(mcTrack.GetPt());
                            hist_gen_r_top->Fill(radius);
                            hist_gen_pt_trit->Fill(dauTrack.GetPt());
                            hist_gen_r_trit->Fill(radius);
                        }
                    }
                }
            }
        }
        // mc end

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
                    auto mcTrack = mcTracksMatrix[evID][trackID];
                    if (abs(mcTrack.GetPdgCode()) == hypPDG)
                    {
                        // hypertriton histos fill
                        auto hypITSTrack = ITStracks->at(iTrack);
                        hist_rec_pt->Fill(mcTrack.GetPt());
                        auto radius = calcRadius(&mcTracksMatrix[evID], mcTrack, tritonPDG);
                        hist_rec_r->Fill(radius);

                        if (!fake)
                        {
                            hist_fake_pt->Fill(mcTrack.GetPt());
                            hist_fake_r->Fill(radius);

                            if (hypITSTrack.getNumberOfClusters() >= 4)
                                hist_ris_4->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                            if (hypITSTrack.getNumberOfClusters() == 3)
                                hist_ris_3->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                        }

                        // topology histos fill

                        int firstDauID = mcTrack.getFirstDaughterTrackId();
                        int nDau = mcTrack.getLastDaughterTrackId();
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
                            continue; // if no triton daughter exits the loop
                        auto dautherTrack = mcTracksMatrix[evID][tritID];
                        for (unsigned int jTrack{0}; jTrack < labITSTPCvec->size(); ++jTrack)
                        {
                            auto tritlab = labITSTPCvec->at(jTrack);
                            int trittrackID, tritevID, tritsrcID;
                            bool tritfake;
                            tritlab.get(trittrackID, tritevID, tritsrcID, tritfake);
                            if (!tritlab.isNoise() && tritlab.isValid())
                            {
                                auto tritmcTrack = mcTracksMatrix[tritevID][trittrackID];
                                if (abs(tritmcTrack.GetPdgCode()) == tritonPDG && trittrackID == tritID && tritevID == evID)
                                {
                                    hist_rec_pt_top->Fill(mcTrack.GetPt());
                                    auto radiusTop = calcRadius(&mcTracksMatrix[evID], mcTrack, tritonPDG);
                                    hist_rec_r_top->Fill(radiusTop);
                                    if (!tritfake && !fake)
                                    {
                                        hist_fake_pt_top->Fill(mcTrack.GetPt());
                                        hist_fake_r_top->Fill(radiusTop);

                                        if (hypITSTrack.getNumberOfClusters() >= 4)
                                            hist_ris_4_top->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                                        if (hypITSTrack.getNumberOfClusters() == 3)
                                            hist_ris_3_top->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned int iTrack{0}; iTrack < labITSTPCvec->size(); ++iTrack)
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
                        hist_rec_pt_trit->Fill(mcTrack.GetPt());

                        auto radius = calcRadius(&mcTracksMatrix[evID], motherTrack, tritonPDG);
                        hist_rec_r_trit->Fill(radius);

                        if (!fake)
                        {
                            hist_fake_pt_trit->Fill(mcTrack.GetPt());
                            hist_fake_r_trit->Fill(radius);

                            hist_ris_trit->Fill((mcTrack.GetPt() - tritITSTPCTrack.getPt()) / mcTrack.GetPt());
                        }
                    }
                }
            }

        } // event loop
    }     // end of tf loop

    auto fFile = TFile(filename, "recreate");
    hist_gen_pt->Write();
    hist_rec_pt->Write();
    hist_fake_pt->Write();
    hist_ris_3->Write();
    hist_ris_4->Write();
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

    TH1F *eff_hist_r = (TH1F *)hist_rec_r->Clone("Hyp Eff r");
    eff_hist_r->GetYaxis()->SetTitle("Efficiency");
    eff_hist_r->SetTitle("Hypertriton Radius Efficiency");
    eff_hist_r->Divide(hist_gen_r);
    eff_hist_r->Write();

    TH1F *pur_hist_r = (TH1F *)hist_fake_r->Clone("Hyp Pur r");
    pur_hist_r->GetYaxis()->SetTitle("Purity");
    pur_hist_r->SetTitle("Hypertriton Radius Purity");
    pur_hist_r->Divide(hist_rec_r);
    pur_hist_r->Write();

    hist_gen_pt_trit->Write();
    hist_rec_pt_trit->Write();
    hist_fake_pt_trit->Write();
    hist_ris_trit->Write();
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

    TH1F *eff_hist_r_trit = (TH1F *)hist_rec_r_trit->Clone("Trit Eff r");
    eff_hist_r_trit->GetYaxis()->SetTitle("Efficiency");
    eff_hist_r_trit->SetTitle("Triton Radius Efficiency");
    eff_hist_r_trit->Divide(hist_gen_r_trit);
    eff_hist_r_trit->Write();

    TH1F *pur_hist_r_trit = (TH1F *)hist_fake_r_trit->Clone("Trit Pur r");
    pur_hist_r_trit->GetYaxis()->SetTitle("Purity");
    pur_hist_r_trit->SetTitle("Triton Radius Purity");
    pur_hist_r_trit->Divide(hist_rec_r_trit);
    pur_hist_r_trit->Write();

    hist_gen_pt_top->Write();
    hist_rec_pt_top->Write();
    hist_fake_pt_top->Write();
    hist_ris_3_top->Write();
    hist_ris_4_top->Write();
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

    TH1F *eff_hist_r_top = (TH1F *)hist_rec_r_top->Clone("Top Eff r");
    eff_hist_r_top->GetYaxis()->SetTitle("Efficiency");
    eff_hist_r_top->SetTitle("Topopolgy Radius Efficiency");
    eff_hist_r_top->Divide(hist_gen_r_top);
    eff_hist_r_top->Write();

    TH1F *pur_hist_r_top = (TH1F *)hist_fake_r_top->Clone("Top Pur r");
    pur_hist_r_top->GetYaxis()->SetTitle("Purity");
    pur_hist_r_top->SetTitle("Topopolgy Radius Purity");
    pur_hist_r_top->Divide(hist_rec_r_top);
    pur_hist_r_top->Write();

    fFile.Close();
} // end of efficiency functions