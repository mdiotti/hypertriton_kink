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

const int hypPDG = 1010010030;
const int tritonPDG = 1000010030;
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

void fit_efficiency(TString path, TString filename, int tf_max = 40)
{

    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;
    int min_r = 0;
    int max_r = 50;

    TH1F *gen_r = new TH1F("Topology gen r", "Topology gen r;r (cm);counts", nBins, min_r, max_r);
    TH1F *rec_r = new TH1F("Topology rec r", "Topology rec r;r (cm);counts", nBins, min_r, max_r);
    TH1F *true_r = new TH1F("Topology true r", "Topology true r;r (cm);counts", nBins, min_r, max_r);
    TH1F *fit_r = new TH1F("Topology fit r", "Topology fit r;r (cm);counts", nBins, min_r, max_r);
    TH1F *true_fit_r = new TH1F("Topology true fit r", "Topology true fit r;r (cm);counts", nBins, min_r, max_r);

    for (int tf = tf_min; tf < tf_max; tf++)
    {
        TString tf_string = Form("%d", tf);
        TString tf_path = path + "tf" + tf_string;

        auto fITS = TFile::Open(tf_path + "/o2trac_its.root");
        auto fITSTPC = TFile::Open(tf_path + "/o2match_itstpc.root");
        auto fMCTracks = TFile::Open(tf_path + "/sgn_" + tf_string + "_Kine.root");
        auto fClusITS = TFile::Open(tf_path + "/o2clus_its.root");

        TString string_to_convert = tf_path + "/o2sim_grp.root";
        std::string path_string(string_to_convert.Data());
        const auto grp = o2::parameters::GRPObject::loadFrom(path_string);

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");

        // Tracks
        std::vector<MCTrack> *MCtracks = nullptr;
        std::vector<TrackITS> *ITStracks = nullptr;
        std::vector<o2::dataformats::TrackTPCITS> *ITSTPCtracks = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;

        // Branches
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

        for (int n = 0; n < nev; n++) // fill Gen histos
        {
            treeMCTracks->GetEvent(n);
            unsigned int size = MCtracks->size();

            for (unsigned int mcI{0}; mcI < size; ++mcI)
            {
                auto mcTrack = mcTracksMatrix[n][mcI];

                if (abs(mcTrack.GetPdgCode()) == hypPDG)
                {
                    auto mcTrack = MCtracks->at(mcI);
                    double rGen = calcRadius(&mcTracksMatrix[n], mcTrack, tritonPDG);
                    int firstDauID = mcTrack.getFirstDaughterTrackId();
                    int nDau = mcTrack.getLastDaughterTrackId();
                    for (int iDau = firstDauID; iDau <= nDau; iDau++)
                    {
                        auto dauTrack = mcTracksMatrix[n][iDau];
                        if (abs(dauTrack.GetPdgCode()) == tritonPDG)
                        {
                            gen_r->Fill(rGen);
                        }
                    }
                }
            }
        }

        for (int event = 0; event < treeITS->GetEntriesFast(); event++)
        {
            if (!treeITS->GetEvent(event) || !treeITSTPC->GetEvent(event) || !treeITSclus->GetEvent(event))
                continue;

            for (unsigned int iTrack{0}; iTrack < labITSvec->size(); ++iTrack)
            {
                auto lab = labITSvec->at(iTrack);
                int trackID, evID, srcID;
                bool fake;
                lab.get(trackID, evID, srcID, fake);
                if (!lab.isNoise() && lab.isValid())
                {
                    double hypgenPabs = 0;
                    TVector3 hypgenP;
                    auto MCTrack = mcTracksMatrix[evID][trackID];
                    if (abs(MCTrack.GetPdgCode()) == hypPDG)
                    {
                        double genR = calcRadius(&mcTracksMatrix[evID], MCTrack, tritonPDG);

                        hypgenPabs = MCTrack.GetP();
                        double hypgenE = sqrt(hypgenPabs * hypgenPabs + hypMassTh * hypMassTh);
                        MCTrack.GetMomentum(hypgenP);
                        auto hypITSTrack = ITStracks->at(iTrack);
                        int nLayers = hypITSTrack.getNumberOfClusters();

                        if (nLayers == 3) // 3 clusters recontruction doesn't work well
                            continue;

                        int firstDauID = MCTrack.getFirstDaughterTrackId();
                        int nDau = MCTrack.getLastDaughterTrackId();
                        int tritID = 0;

                        double tritgenPabs = 0;
                        TVector3 tritgenP;
                        for (int iDau = firstDauID; iDau <= nDau; iDau++)
                        {
                            if (mcTracksMatrix[evID][iDau].GetPdgCode() == tritonPDG)
                            {
                                tritID = iDau;
                                tritgenPabs = mcTracksMatrix[evID][iDau].GetP();
                                mcTracksMatrix[evID][iDau].GetMomentum(tritgenP);
                                break;
                            }
                        }

                        if (tritID == 0)
                            continue; // if no triton daughter, improves speed

                        bool counted = false;
                        bool countedTrue = false;

                        for (unsigned int jTrack{0}; jTrack < labITSTPCvec->size(); ++jTrack)
                        {
                            bool isDaughter = false;

                            auto tritlab = labITSTPCvec->at(jTrack);
                            int trittrackID, tritevID, tritsrcID;
                            bool tritfake;
                            tritlab.get(trittrackID, tritevID, tritsrcID, tritfake);
                            if (!tritlab.isNoise() && tritlab.isValid())
                            {
                                auto tritMCTrack = mcTracksMatrix[tritevID][trittrackID];
                                if (abs(tritMCTrack.GetPdgCode()) == tritonPDG)
                                {
                                    if (!counted)
                                    {
                                        rec_r->Fill(genR);
                                        counted = true;
                                    }

                                    if (trittrackID == tritID && tritevID == evID)
                                        isDaughter = true;

                                    bool isSignal = false;

                                    if (isDaughter && !tritfake && !fake)
                                        isSignal = true;

                                    auto tritITSTPCtrack = ITSTPCtracks->at(jTrack);
                                    // if (!tritfake && !fake)
                                    if (true) // fits also the fake tracks
                                    {

                                        if (!tritfake && !fake && !countedTrue)
                                        {
                                            true_r->Fill(genR);
                                            countedTrue = true;
                                        }
                                        // Fitting start

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

                                                    std::array<float, 3> R = ft2.getPCACandidatePos();
                                                    double recR = sqrt(R[0] * R[0] + R[1] * R[1]);

                                                    if (recR < 18)
                                                        continue;
                                                    if (std::abs(etaHyp - etaTrit) > 0.03 || std::abs(phiHyp - phiTrit) > 0.03)
                                                        continue;

                                                    float tritE = sqrt(tritPabs * tritPabs + tritonMass * tritonMass);
                                                    std::array<float, 3> piP = {hypP[0] - tritP[0], hypP[1] - tritP[1], hypP[2] - tritP[2]};
                                                    float piPabs = sqrt(piP[0] * piP[0] + piP[1] * piP[1] + piP[2] * piP[2]);
                                                    float piE = sqrt(pi0Mass * pi0Mass + piPabs * piPabs);
                                                    float hypE = piE + tritE;
                                                    float hypMass = sqrt(hypE * hypE - hypPabs * hypPabs);
                                                    double pt = sqrt(hypP[0] * hypP[0] + hypP[1] * hypP[1]);

                                                    if (hypPabs < tritPabs)
                                                        continue;

                                                    double rRec = calcRadius(&mcTracksMatrix[evID], MCTrack, tritonPDG);
                                                    fit_r->Fill(rRec);
                                                    if (!tritfake && !fake)
                                                        true_fit_r->Fill(rRec);
                                                }
                                            }
                                            catch (std::runtime_error &e)
                                            {
                                                continue;
                                            }
                                        }
                                    }
                                } // Fitting end
                            }
                        }
                    }
                }
            }
        }
    }

    auto fFile = TFile(filename, "recreate");

    gen_r->Write();
    rec_r->Write();
    true_r->Write();
    fit_r->Write();
    true_fit_r->Write();

    TH1F *eff_r = (TH1F *)rec_r->Clone("Top Eff r");
    eff_r->GetXaxis()->SetTitleSize(fontSize);
    eff_r->GetYaxis()->SetTitleSize(fontSize);
    eff_r->GetYaxis()->SetTitle("Efficiency");
    eff_r->SetTitle("Topology r Efficiency");
    eff_r->Divide(gen_r);
    eff_r->Write();

    TH1F *eff_true_r = (TH1F *)true_r->Clone("Top Eff true r");
    eff_true_r->GetXaxis()->SetTitleSize(fontSize);
    eff_true_r->GetYaxis()->SetTitleSize(fontSize);
    eff_true_r->GetYaxis()->SetTitle("Efficiency");
    eff_true_r->SetTitle("Topology true r Efficiency");
    eff_true_r->Divide(gen_r);
    eff_true_r->Write();

    TH1F *eff_fit_r = (TH1F *)fit_r->Clone("Top Eff fit r");
    eff_fit_r->GetXaxis()->SetTitleSize(fontSize);
    eff_fit_r->GetYaxis()->SetTitleSize(fontSize);
    eff_fit_r->GetYaxis()->SetTitle("Efficiency");
    eff_fit_r->SetTitle("Topology Fit r Efficiency");
    eff_fit_r->Divide(gen_r);
    eff_fit_r->Write();

    TH1F *eff_true_fit_r = (TH1F *)true_fit_r->Clone("Top Eff true fit r");
    eff_true_fit_r->GetXaxis()->SetTitleSize(fontSize);
    eff_true_fit_r->GetYaxis()->SetTitleSize(fontSize);
    eff_true_fit_r->GetYaxis()->SetTitle("Efficiency");
    eff_true_fit_r->SetTitle("Topology Fit True r Efficiency");
    eff_true_fit_r->Divide(gen_r);
    eff_true_fit_r->Write();

    fFile.Close();
}