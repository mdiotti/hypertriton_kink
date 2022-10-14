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

TString p_filename = "Momentum distributions.txt";

const int hypPDG = 1010010030;
const int tritonPDG = 1000010030;
const int pi0PDG = 111;
const double tritonMass = 2.808921;
const double pi0Mass = 0.1349766;
TString chiLabel = "#chi^{2}";
TString hypLabel = "M_{^{3}_{#Lambda}H} (GeV/c^{2})";
int nBins = 100;
double min_bins = 0;
double min_r = 0;
double res_bin_lim = 0.25;
double eta_bin_lim = 0.1;
double phi_bin_lim = 0.1;

double lim0 = 0.01;
TString lim0Label = "0.01";
double lim = 0.02;
TString limLabel = Form("%f", lim);
double lim2 = 0.05;
TString limLabel2 = Form("%f", lim2);
double lim3 = 0.1;
TString limLabel3 = Form("%f", lim3);

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

void fitting(TString path, TString filename, int tf_max = 40, bool cut = true)
{
    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;
    if (cut)
    {
        min_r = 18;
        eta_bin_lim = 0.03;
        phi_bin_lim = 0.03;
    }

    TH1F *daughter_chi = new TH1F("Daughters chi2", "Daughters " + chiLabel + ";" + chiLabel + ";counts", nBins, min_bins, 2);
    TH1F *nondaughter_chi = new TH1F("Non-daughters chi2", "Non-daughters " + chiLabel + ";" + chiLabel + ";counts", nBins, 0, 50);
    TH1F *nondaughter_chi_normalized = new TH1F("Non-daughters chi2 normalized", chiLabel + " normalized;" + chiLabel + ";counts", nBins, min_bins, 2);
    TH1F *resolution = new TH1F("Resolution", "Resolution;#Delta r;counts", nBins, -res_bin_lim, res_bin_lim);
    TH1F *daughter_radius = new TH1F("Daughter radius", "Daughter radius;Rrec(cm);counts", nBins, min_r, 50);
    TH1F *nondaughter_radius = new TH1F("Non-daughter radius", "Non-daughter radius;Rrec(cm);counts", nBins, min_r, 50);
    TH1F *inv_mass = new TH1F("Invariant mass", "Invariant mass;" + hypLabel + ";counts", nBins, 2.5, 3.5);
    TH1F *inv_mass_daughter = new TH1F("Invariant mass daughter", "Invariant mass daughter;" + hypLabel + ";counts", nBins, 2.5, 3.5);
    TH1F *inv_mass_nondaughter = new TH1F("Invariant mass non-daughter", "Invariant mass non-daughter;" + hypLabel + ";counts", nBins, 2.5, 3.5);
    TH1F *pi0_p_resolution = new TH1F("Pi0 p resolution", "Pi0 p resolution;Resolution;counts", nBins, -30, 30);
    TH1F *hyp_p_resolution = new TH1F("Hyp p resolution", "Hyp p resolution;Resolution;counts", nBins, -1, 1);
    TH1F *triton_p_resolution = new TH1F("Triton p resolution", "Triton p resolution;Resolution;counts", nBins, -1, 1);
    TH1F *pi0_partial_resolution0 = new TH1F("Pi0 partial resolution0", "#pi^{0} p resolution with hyp and trit res < " + limLabel + ";Resolution;counts", nBins, -2, 2);
    TH1F *pi0_partial_resolution = new TH1F("Pi0 partial resolution", "#pi^{0} p resolution with hyp and trit res < " + limLabel + ";Resolution;counts", nBins, -2, 2);
    TH1F *pi0_partial_resolution2 = new TH1F("Pi0 p resolution2", "#pi^{0} p resolution with hyp and trit res < " + limLabel2 + ";Resolution;counts", nBins, -2, 2);
    TH1F *pi0_partial_resolution3 = new TH1F("Pi0 p resolution3", "#pi^{0} p resolution with hyp and trit res < " + limLabel3 + ";Resolution;counts", nBins, -2, 2);

    TH1F *pi0_resolution_zoomed = new TH1F("Pi0 resolution zoommato", "#pi^{0} p resolution zoomed;Resolution;counts", nBins, -2, 2);

    TH2F *resolution_vs_chi = new TH2F("Resolution vs chi2", "Resolution vs " + chiLabel + ";Resolution;" + chiLabel, nBins, -res_bin_lim, res_bin_lim, nBins, min_bins, 2);
    TH2F *eta_vs_phi = new TH2F("Eta vs Phi daughter", "Eta vs Phi daughter;#eta;#phi", nBins, -eta_bin_lim, eta_bin_lim, nBins, -phi_bin_lim, phi_bin_lim);
    TH2F *eta_vs_phi_nondaughter = new TH2F("Eta vs Phi non-daughter", "Eta vs Phi non-daughter;#eta;#phi", nBins, -eta_bin_lim, eta_bin_lim, nBins, -phi_bin_lim, phi_bin_lim);
    TH2F *resolution_vs_rrec = new TH2F("Resolution vs Rrec", "Resolution vs Rrec;Resolution;Rrec(cm)", nBins, -res_bin_lim, res_bin_lim, nBins, min_r, 50);
    TH2F *pi0_p_resolution_vs_pgen = new TH2F("Pi0 p resolution vs pgen", "#pi^{0} p resolution vs pgen;Resolution;Pgen(GeV/c)", nBins, -30, 30, nBins, 0, 30);
    TH2F *hyp_p_resolution_vs_pgen = new TH2F("Hyp p resolution vs pgen", "Hyp p resolution vs pgen;Resolution;Pgen(GeV/c)", nBins, -1, 1, nBins, 0, 30);
    TH2F *triton_p_resolution_vs_pgen = new TH2F("Triton p resolution vs pgen", "Triton p resolution vs pgen;Resolution;Pgen(GeV/c)", nBins, -1, 1, nBins, 0, 30);

     ofstream oFile(p_filename);

    for (int tf = tf_min; tf < tf_max; tf++)
    {
        TString tf_string = Form("%d", tf);
        TString tf_path = path + "tf" + tf_string;

        auto fITS = TFile::Open(tf_path + "/o2trac_its.root");
        auto fITSTPC = TFile::Open(tf_path + "/o2match_itstpc.root");
        auto fMCTracks = TFile::Open(tf_path + "/sgn_" + tf_string + "_Kine.root");

        TString string_to_convert = tf_path + "/o2sim_grp.root";
        std::string path_string(string_to_convert.Data());
        const auto grp = o2::parameters::GRPObject::loadFrom(path_string);

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
                    double hypgenPabs = 0;
                    TVector3 hypgenP;
                    auto MCTrack = mcTracksMatrix[evID][trackID];
                    if (abs(MCTrack.GetPdgCode()) == hypPDG)
                    {
                        hypgenPabs = MCTrack.GetP();
                        MCTrack.GetMomentum(hypgenP);
                        auto hypITSTrack = ITStracks->at(iTrack);
                        if (hypITSTrack.getNumberOfClusters() == 3) // 3 clusters recontruction doesn't work well
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

                        double pi0genPabs = 0;
                        TVector3 pi0genP;
                        for (int iDau = firstDauID; iDau <= nDau; iDau++)
                        {
                            if (iDau == tritID)
                                continue;

                            if (abs(mcTracksMatrix[evID][iDau].GetPdgCode()) == pi0PDG)
                            {
                                pi0genPabs = mcTracksMatrix[evID][iDau].GetP();
                                mcTracksMatrix[evID][iDau].GetMomentum(pi0genP);
                            }
                        }

                        double genR = calcRadius(&mcTracksMatrix[evID], MCTrack, tritonPDG);

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
                                    if (trittrackID == tritID && tritevID == evID)
                                        isDaughter = true;

                                    auto tritITSTPCtrack = ITSTPCtracks->at(jTrack);
                                    if (!tritfake && !fake)
                                    {
                                        // Fitting start

                                        if (FITTEROPTION == "KFParticle")
                                        {
                                        }

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

                                                    if (cut)
                                                    {
                                                        if (recR < 18)
                                                            continue;
                                                        if (std::abs(etaHyp - etaTrit) > 0.03 || std::abs(phiHyp - phiTrit) > 0.03)
                                                            continue;
                                                    }

                                                    double res = (genR - recR) / genR;

                                                    float tritE = sqrt(tritPabs * tritPabs + tritonMass * tritonMass);
                                                    std::array<float, 3> piP = {hypP[0] - tritP[0], hypP[1] - tritP[1], hypP[2] - tritP[2]};
                                                    float piPabs = sqrt(piP[0] * piP[0] + piP[1] * piP[1] + piP[2] * piP[2]);
                                                    float piE = sqrt(pi0Mass * pi0Mass + piPabs * piPabs);
                                                    float hypE = piE + tritE;

                                                    float hypMass = sqrt(hypE * hypE - hypPabs * hypPabs);

                                                    if (isDaughter)
                                                    {
                                                        daughter_chi->Fill(chi2);
                                                        resolution->Fill(res);
                                                        resolution_vs_chi->Fill(res, chi2);
                                                        eta_vs_phi->Fill((etaHyp - etaTrit), (phiHyp - phiTrit));
                                                        daughter_radius->Fill(recR);
                                                        inv_mass_daughter->Fill(hypMass);
                                                    }
                                                    else
                                                    {
                                                        nondaughter_chi->Fill(chi2);
                                                        nondaughter_chi_normalized->Fill(chi2);
                                                        eta_vs_phi_nondaughter->Fill((etaHyp - etaTrit), (phiHyp - phiTrit));
                                                        nondaughter_radius->Fill(recR);
                                                        inv_mass_nondaughter->Fill(hypMass);
                                                    }
                                                    resolution_vs_rrec->Fill(res, recR);
                                                    inv_mass->Fill(hypMass);
                                                    double pi0_p_res = (pi0genPabs - piPabs) / pi0genPabs;
                                                    double trit_p_res = (tritgenPabs - tritPabs) / tritgenPabs;
                                                    double hyp_p_res = (hypgenPabs - hypPabs) / hypgenPabs;

                                                    pi0_p_resolution_vs_pgen->Fill(pi0_p_res, pi0genPabs);
                                                    triton_p_resolution_vs_pgen->Fill(trit_p_res, tritgenPabs);
                                                    hyp_p_resolution_vs_pgen->Fill(hyp_p_res, hypgenPabs);

                                                    TString pi0String = "pi0 resolution = "+ std::to_string(pi0_p_res) + " Generated: px= " + std::to_string(pi0genP[0]) + " py= " + std::to_string(pi0genP[1]) + " pz= " + std::to_string(pi0genP[2]) + " Reconstructed: px= " + std::to_string(piP[0]) + " py= " + std::to_string(piP[1]) + " pz= " + std::to_string(piP[2]);
                                                    TString tritString = "Triton resolution = "+ std::to_string(trit_p_res) +" Generated: px= " + std::to_string(tritgenP[0]) + " py= " + std::to_string(tritgenP[1]) + " pz= " + std::to_string(tritgenP[2]) + " Reconstructed: px= " + std::to_string(tritP[0]) + " py= " + std::to_string(tritP[1]) + " pz= " + std::to_string(tritP[2]);
                                                    TString hypString = "Hyp resolution = "+ std::to_string(hyp_p_res) +" Generated: px= " + std::to_string(hypgenP[0]) + " py= " + std::to_string(hypgenP[1]) + " pz= " + std::to_string(hypgenP[2]) + " Reconstructed: px= " + std::to_string(hypP[0]) + " py= " + std::to_string(hypP[1]) + " pz= " + std::to_string(hypP[2]);

                                                    if (abs(trit_p_res) <= lim0 && abs(hyp_p_res) <= lim0)
                                                        pi0_partial_resolution0->Fill(pi0_p_res);

                                                    if (abs(trit_p_res) <= lim && abs(hyp_p_res) <= lim)
                                                        pi0_partial_resolution->Fill(pi0_p_res);

                                                    if (abs(trit_p_res) <= lim2 && abs(hyp_p_res) <= lim2)
                                                        pi0_partial_resolution2->Fill(pi0_p_res);

                                                    if (abs(trit_p_res) <= lim3 && abs(hyp_p_res) <= lim3)
                                                    {
                                                        pi0_partial_resolution3->Fill(pi0_p_res);

                                                        oFile << hypString << endl
                                                              << tritString << endl
                                                              << pi0String << endl
                                                              << endl;
                                                    }

                                                    // if (abs(pi0_p_res) >= 5)
                                                    //  cout << "pi0 res >= 5; trit_p_res: " << trit_p_res << " hyp_p_res: " << hyp_p_res << " pi0_p_res: " << pi0_p_res << endl;

                                                    pi0_p_resolution->Fill(pi0_p_res);
                                                    triton_p_resolution->Fill(trit_p_res);
                                                    hyp_p_resolution->Fill(hyp_p_res);

                                                    pi0_resolution_zoomed->Fill(pi0_p_res);
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
        } // end of event loop
    }     // end of tf loop

    oFile.close();

    auto fFile = TFile(filename, "recreate");
    daughter_chi->Write();
    nondaughter_chi->Write();
    resolution->Write();
    resolution_vs_chi->Write();
    eta_vs_phi->Write();
    eta_vs_phi_nondaughter->Write();
    resolution_vs_rrec->Write();

    // inv_mass->Fit("gaus");
    inv_mass->Write();

    inv_mass_daughter->Write();
    inv_mass_nondaughter->Write();
    pi0_p_resolution->Write();
    triton_p_resolution->Write();
    hyp_p_resolution->Write();
    pi0_partial_resolution0->Write();
    pi0_partial_resolution->Write();
    pi0_partial_resolution2->Write();
    pi0_partial_resolution3->Write();

    pi0_p_resolution_vs_pgen->Write();
    triton_p_resolution_vs_pgen->Write();
    hyp_p_resolution_vs_pgen->Write();

    auto c = new TCanvas("c", "c", 800, 600);
    nondaughter_chi_normalized->SetLineColor(kRed);
    nondaughter_chi_normalized->DrawNormalized("P");
    daughter_chi->DrawNormalized("sameP");
    c->Write();
    /*
        auto c1 = new TCanvas("c", "c", 800, 600);
        nondaughter_radius->SetLineColor(kRed);
        nondaughter_radius->SetTitle("Radius normalized");
        nondaughter_radius->DrawNormalized("P");
        daughter_radius->DrawNormalized("sameP");
        c1->Write();
    */

    fFile.Close();
} // end of fitting function
