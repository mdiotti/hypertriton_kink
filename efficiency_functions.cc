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

const int tf_min = 1;
const int tf_max = 40;
int tf_lenght = tf_max - tf_min + 1;

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
            auto decLength = TMath::Sqrt((dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) * (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) + (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) * (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()));
            return decLength;
        }
    }
    return -1;
}

void efficiency_functions(TString path, TString filename)
{

    // define hypertriton track histograms
    TH1F *hist_gen_pt = new TH1F("gen_pt", ";" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_rec_pt = new TH1F("rec_pt", ";" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_fake_pt = new TH1F("fake_pt", ";" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_ris_pt = new TH1F("res_pt", "gen_pt - rec_pt;" + ptLabel + ";counts", 30, -5, 5);
    TH1F *hist_ris_pt_perc = new TH1F("res_pt_perc", "Resolution;resolution;counts", 30, -1, 1);
    TH1F *hist_gen_r = new TH1F("gen_r", ";Radius (cm);counts", 50, 0, 50);
    TH1F *hist_rec_r = new TH1F("rec_r", ";Radius (cm);counts", 50, 0, 50);
    TH1F *hist_fake_r = new TH1F("fake_r", ";Radius (cm);counts", 50, 0, 50);

    // define triton track histograms
    TH1F *hist_gen_pt_trit = new TH1F("gen_pt_trit", ";" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_rec_pt_trit = new TH1F("rec_pt_trit", ";" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_fake_pt_trit = new TH1F("fake_pt_trit", ";" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_ris_pt_trit = new TH1F("res_pt_trit", "gen_pt - rec_pt;" + ptLabel + ";counts", 30, -5, 5);
    TH1F *hist_ris_pt_perc_trit = new TH1F("res_pt_perc_trit", "Resolution;resolution;counts", 30, -1, 1);
    TH1F *hist_gen_r_trit = new TH1F("gen_r_trit", ";Radius (cm);counts", 50, 0, 50);
    TH1F *hist_rec_r_trit = new TH1F("rec_r_trit", ";Radius (cm);counts", 50, 0, 50);
    TH1F *hist_fake_r_trit = new TH1F("fake_r_trit", ";Radius (cm);counts", 50, 0, 50);

    // define topology histograms
    TH1F *hist_gen_pt_top = new TH1F("gen_pt_top", ";" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_rec_pt_top = new TH1F("rec_pt_top", ";" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_fake_pt_top = new TH1F("fake_pt_top", ";" + ptLabel + ";counts", 30, 1, 10);
    TH1F *hist_ris_pt_top = new TH1F("res_pt_top", "gen_pt - rec_pt;" + ptLabel + ";counts", 30, -5, 5);
    TH1F *hist_ris_pt_perc_top = new TH1F("res_pt_perc_top", "Resolution;resolution;counts", 30, -1, 1);
    TH1F *hist_gen_r_top = new TH1F("gen_r_top", ";Radius (cm);counts", 50, 0, 50);
    TH1F *hist_rec_r_top = new TH1F("rec_r_top", ";Radius (cm);counts", 50, 0, 50);
    TH1F *hist_fake_r_top = new TH1F("fake_r_top", ";Radius (cm);counts", 50, 0, 50);

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

        std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}};

        // mc start

        std::vector<std::vector<MCTrack>> mcTracksMatrix;
        auto nev = treeMCTracks->GetEntriesFast();

        mcTracksMatrix.resize(nev);
        for (int n = 0; n < nev; n++)
        { // loop over MC events
            treeMCTracks->GetEvent(n);

            mcTracksMatrix[n].resize(MCtracks->size());

            for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI)
            {
                auto mcTrack = MCtracks->at(mcI);
                mcTracksMatrix[n][mcI] = mcTrack;

                if (abs(mcTrack.GetPdgCode()) == tritonPDG)
                {
                    int motherID = mcTrack.getMotherTrackId();
                    auto motherTrack = mcTracksMatrix[n][motherID];
                    if (abs(motherTrack.GetPdgCode()) == hypPDG)
                    {
                        hist_gen_pt_trit->Fill(motherTrack.GetPt());
                        hist_gen_r_trit->Fill(calcRadius(&mcTracksMatrix[n], motherTrack, tritonPDG));
                    }
                }
                else if (abs(mcTrack.GetPdgCode()) == hypPDG)
                {
                    hist_gen_pt->Fill(mcTrack.GetPt());
                    hist_gen_r->Fill(calcRadius(&mcTracksMatrix[n], mcTrack, tritonPDG));

                    int dauID = mcTrack.getFirstDaughterTrackId();
                    auto dauTrack = mcTracksMatrix[n][dauID];
                    if(abs(dauTrack.GetPdgCode()) == tritonPDG){
                        hist_gen_pt_top->Fill(mcTrack.GetPt());
                        hist_gen_r_top->Fill(calcRadius(&mcTracksMatrix[n], mcTrack, tritonPDG));
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
                        auto decLength = calcRadius(&mcTracksMatrix[evID], mcTrack, tritonPDG);
                        hist_rec_r->Fill(decLength);

                        if (!fake)
                        {
                            hist_fake_pt->Fill(mcTrack.GetPt());
                            hist_fake_r->Fill(decLength);

                            hist_ris_pt->Fill(mcTrack.GetPt() - hypITSTrack.getPt());
                            hist_ris_pt_perc->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
                        }

                        // topology histos fill
                        int dauID = mcTrack.getFirstDaughterTrackId();
                        auto dautherTrack = mcTracksMatrix[evID][dauID];

                        for (unsigned int jTrack{0}; jTrack < labITSTPCvec->size(); ++jTrack)
                        {
                            auto tritlab = labITSTPCvec->at(jTrack);
                            int trittrackID, tritevID, tritsrcID;
                            bool tritfake;
                            tritlab.get(trittrackID, tritevID, tritsrcID, tritfake);
                            if (!tritlab.isNoise() && tritlab.isValid())
                            {
                                auto tritmcTrack = mcTracksMatrix[tritevID][trittrackID];
                                if (abs(tritmcTrack.GetPdgCode()) == tritonPDG && trittrackID == dauID && tritevID == evID)
                                {
                                    hist_rec_pt_top->Fill(mcTrack.GetPt());
                                    auto decLengthTop = calcRadius(&mcTracksMatrix[evID], mcTrack, tritonPDG);
                                    hist_rec_r_top->Fill(decLengthTop);
                                    if (!tritfake && !fake)
                                    {
                                        hist_fake_pt_top->Fill(mcTrack.GetPt());
                                        hist_fake_r_top->Fill(decLengthTop);

                                        hist_ris_pt_top->Fill(mcTrack.GetPt() - hypITSTrack.getPt());
                                        hist_ris_pt_perc_top->Fill((mcTrack.GetPt() - hypITSTrack.getPt()) / mcTrack.GetPt());
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

                        auto decLength = calcRadius(&mcTracksMatrix[evID], mcTrack, tritonPDG);
                        hist_rec_r_trit->Fill(decLength);

                        if (!fake)
                        {
                            hist_fake_pt_trit->Fill(mcTrack.GetPt());
                            hist_fake_r_trit->Fill(decLength);

                            hist_ris_pt_trit->Fill(mcTrack.GetPt() - tritITSTPCTrack.getPt());
                            hist_ris_pt_perc_trit->Fill((mcTrack.GetPt() - tritITSTPCTrack.getPt()) / mcTrack.GetPt());
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
    hist_ris_pt->Write();
    hist_ris_pt_perc->Write();
    hist_gen_r->Write();
    hist_rec_r->Write();
    hist_fake_r->Write();

    TH1F *eff_hist_pt = (TH1F *)hist_rec_pt->Clone("hist_eff_pt");
    eff_hist_pt->GetYaxis()->SetTitle("Hypertriton Efficiency");
    eff_hist_pt->Divide(hist_gen_pt);
    eff_hist_pt->Write();

    TH1F *pur_hist_pt = (TH1F *)hist_fake_pt->Clone("hist_pur_pt");
    pur_hist_pt->GetYaxis()->SetTitle("Hypertriton Purity");
    pur_hist_pt->Divide(hist_rec_pt);
    pur_hist_pt->Write();

    TH1F *eff_hist_r = (TH1F *)hist_rec_r->Clone("hist_eff_r");
    eff_hist_r->GetYaxis()->SetTitle("Hypertriton Efficiency");
    eff_hist_r->Divide(hist_gen_r);
    eff_hist_r->Write();

    TH1F *pur_hist_r = (TH1F *)hist_fake_r->Clone("hist_pur_r");
    pur_hist_r->GetYaxis()->SetTitle("Hypertriton Purity");
    pur_hist_r->Divide(hist_rec_r);
    pur_hist_r->Write();

    hist_gen_pt_trit->Write();
    hist_rec_pt_trit->Write();
    hist_fake_pt_trit->Write();
    hist_ris_pt_trit->Write();
    hist_ris_pt_perc_trit->Write();
    hist_gen_r_trit->Write();
    hist_rec_r_trit->Write();
    hist_fake_r_trit->Write();

    TH1F *eff_hist_pt_trit = (TH1F *)hist_rec_pt_trit->Clone("hist_eff_pt_trit");
    eff_hist_pt_trit->GetYaxis()->SetTitle("Triton Efficiency");
    eff_hist_pt_trit->Divide(hist_gen_pt_trit);
    eff_hist_pt_trit->Write();

    TH1F *pur_hist_pt_trit = (TH1F *)hist_fake_pt_trit->Clone("hist_pur_pt_trit");
    pur_hist_pt_trit->GetYaxis()->SetTitle("Triton Purity");
    pur_hist_pt_trit->Divide(hist_rec_pt_trit);
    pur_hist_pt_trit->Write();

    TH1F *eff_hist_dl_trit = (TH1F *)hist_rec_r_trit->Clone("hist_eff_r_trit");
    eff_hist_dl_trit->GetYaxis()->SetTitle("Triton Efficiency");
    eff_hist_dl_trit->Divide(hist_gen_r_trit);
    eff_hist_dl_trit->Write();

    TH1F *pur_hist_dl_trit = (TH1F *)hist_fake_r_trit->Clone("hist_pur_r_trit");
    pur_hist_dl_trit->GetYaxis()->SetTitle("Triton Purity");
    pur_hist_dl_trit->Divide(hist_rec_r_trit);
    pur_hist_dl_trit->Write();

    hist_gen_pt_top->Write();
    hist_rec_pt_top->Write();
    hist_fake_pt_top->Write();
    hist_ris_pt_top->Write();
    hist_ris_pt_perc_top->Write();
    hist_gen_r_top->Write();
    hist_rec_r_top->Write();
    hist_fake_r_top->Write();

    TH1F *eff_hist_pt_top = (TH1F *)hist_rec_pt_top->Clone("hist_eff_pt_top");
    eff_hist_pt_top->GetYaxis()->SetTitle("Top Efficiency");
    eff_hist_pt_top->Divide(hist_gen_pt_top);
    eff_hist_pt_top->Write();

    TH1F *pur_hist_pt_top = (TH1F *)hist_fake_pt_top->Clone("hist_pur_pt_top");
    pur_hist_pt_top->GetYaxis()->SetTitle("Top Purity");
    pur_hist_pt_top->Divide(hist_rec_pt_top);
    pur_hist_pt_top->Write();

    TH1F *eff_hist_r_top = (TH1F *)hist_rec_r_top->Clone("hist_eff_r_top");
    eff_hist_r_top->GetYaxis()->SetTitle("Top Efficiency");
    eff_hist_r_top->Divide(hist_gen_r_top);
    eff_hist_r_top->Write();

    fFile.Close();
} // end of efficiency functions