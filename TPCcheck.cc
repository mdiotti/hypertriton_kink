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
#include "DataFormatsTPC/TrackTPC.h"

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
using TrackTPC = o2::tpc::TrackTPC;

const int hypPDG = 1010010030;
const int tritonPDG = 1000010030;
const int pi0PDG = 111;
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

double calcDecayLenght(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        //        if (dauTrack.GetPdgCode() == dauPDG)
        if (abs(dauTrack.GetPdgCode()) == dauPDG)
        {
            auto dl = TMath::Sqrt((dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) * (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) + (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) * (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) + (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ()) * (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ()));
            return dl;
        }
    }
    return -1;
}

void TPCcheck(TString path, TString filename, int tf_max = 80)
{
    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;

    TH1F *trit_gen_pt = new TH1F("Trit gen pt", "Trit gen pt;p_{T} (GeV/c);counts", nBins, 0, 12);
    TH1F *trit_gen_r = new TH1F("Trit gen r", "Trit gen r;r (cm);counts", nBins, 0, 50);
    TH1F *trit_rec_pt = new TH1F("Trit rec pt", "Trit rec pt;p_{T} (GeV/c);counts", nBins, 0, 12);
    TH1F *trit_rec_r = new TH1F("Trit rec r", "Trit rec r;r (cm);counts", nBins, 0, 50);
    TH1F *trit_true_pt = new TH1F("Trit true pt", "Trit true pt;p_{T} (GeV/c);counts", nBins, 0, 12);
    TH1F *trit_true_r = new TH1F("Trit true r", "Trit true r;r (cm);counts", nBins, 0, 50);

    for (int tf = tf_min; tf <= tf_max; tf++)
    {
        LOG(info) << "Processing TF " << tf;
        TString tf_string = Form("%d", tf);
        TString tf_path = path + "tf" + tf_string;

        auto fITS = TFile::Open(tf_path + "/o2trac_its.root");
        auto fITSTPC = TFile::Open(tf_path + "/o2match_itstpc.root");
        auto fTPC = TFile::Open(tf_path + "/tpctracks.root");
        auto fMCTracks = TFile::Open(tf_path + "/sgn_" + tf_string + "_Kine.root");

        TString string_to_convert = tf_path + "/o2sim_grp.root";
        std::string path_string(string_to_convert.Data());
        const auto grp = o2::parameters::GRPObject::loadFrom(path_string);

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeTPC = (TTree *)fTPC->Get("tpcrec");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");

        // Tracks
        std::vector<MCTrack> *MCtracks = nullptr;
        std::vector<TrackTPC> *TPCtracks = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labTPCvec = nullptr;

        // Branches
        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeTPC->SetBranchAddress("TPCTracks", &TPCtracks);

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
                auto mcTrack = MCtracks->at(mcI);
                if (abs(mcTrack.GetPdgCode()) == tritonPDG)
                {
                    auto hypID = mcTrack.getMotherTrackId();
                    // auto hypTrack = mcTracksMatrix[n][hypID];
                    auto hypTrack = MCtracks->at(hypID);

                    if (abs(hypTrack.GetPdgCode()) != hypPDG)
                        continue;

                    trit_gen_pt->Fill(mcTrack.GetPt());
                    // double rGen = calcRadius(&mcTracksMatrix[n], hypTrack, tritonPDG);
                    double rGen = calcRadius(MCtracks, hypTrack, tritonPDG);
                    trit_gen_r->Fill(rGen);
                }
            }
        }

        for (int event = 0; event < treeTPC->GetEntriesFast(); event++)
        {
            if (!treeTPC->GetEvent(event))
                continue;

            for (unsigned int iTrack{0}; iTrack < labTPCvec->size(); ++iTrack)
            {
                auto lab = labTPCvec->at(iTrack);
                int trackID, evID, srcID;
                bool fake;
                lab.get(trackID, evID, srcID, fake);
                if (!lab.isNoise() && lab.isValid())
                {
                    auto MCTrack = mcTracksMatrix[evID][trackID];
                    if (abs(MCTrack.GetPdgCode()) == tritonPDG)
                    {
                        auto hypID = MCTrack.getMotherTrackId();
                        auto hypTrack = mcTracksMatrix[evID][hypID];
                        if (abs(hypTrack.GetPdgCode()) != hypPDG)
                            continue;

                        auto tritTPCTrack = TPCtracks->at(iTrack);
                        double recR = calcRadius(&mcTracksMatrix[evID], hypTrack, tritonPDG);
                        double tritPt = MCTrack.GetPt();
                        trit_rec_pt->Fill(tritPt);
                        trit_rec_r->Fill(recR);

                        if (!fake)
                        {
                            trit_true_pt->Fill(tritPt);
                            trit_true_r->Fill(recR);
                        }
                    }
                }
            }
        }
    }

    auto fFile = TFile(filename, "recreate");
    trit_gen_pt->Write();
    trit_rec_pt->Write();
    trit_true_pt->Write();
    trit_gen_r->Write();
    trit_rec_r->Write();
    trit_true_r->Write();

    TH1F *eff_pt = (TH1F *)trit_rec_pt->Clone("Trit Eff p");
    eff_pt->GetYaxis()->SetTitle("Efficiency");
    eff_pt->SetTitle("Triton p_{T} Efficiency");
    eff_pt->Divide(trit_gen_pt);
    eff_pt->Write();

    TH1F *eff_r = (TH1F *)trit_rec_r->Clone("Trit Eff r");
    eff_r->GetYaxis()->SetTitle("Efficiency");
    eff_r->SetTitle("Triton r Efficiency");
    eff_r->Divide(trit_gen_r);
    eff_r->Write();

    TH1F *pur_pt = (TH1F *)trit_true_pt->Clone("Trit Pur p");
    pur_pt->GetYaxis()->SetTitle("Purity");
    pur_pt->SetTitle("Triton p_{T} Purity");
    pur_pt->Divide(trit_rec_pt);
    pur_pt->Write();

    TH1F *pur_r = (TH1F *)trit_true_r->Clone("Trit Pur r");
    pur_r->GetYaxis()->SetTitle("Purity");
    pur_r->SetTitle("Triton r Purity");
    pur_r->Divide(trit_rec_r);
    pur_r->Write();

    fFile.Close();
}