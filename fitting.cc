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
    TH1F *pi0_p_resolution = new TH1F("Pi0 P resolution", "Pi0 P resolution;Resolution;counts", nBins, -50, 50);

    TH2F *resolution_vs_chi = new TH2F("Resolution vs chi2", "Resolution vs " + chiLabel + ";Resolution;" + chiLabel, nBins, -res_bin_lim, res_bin_lim, nBins, min_bins, 2);
    TH2F *eta_vs_phi = new TH2F("Eta vs Phi daughter", "Eta vs Phi daughter;#eta;#phi", nBins, -eta_bin_lim, eta_bin_lim, nBins, -phi_bin_lim, phi_bin_lim);
    TH2F *eta_vs_phi_nondaughter = new TH2F("Eta vs Phi non-daughter", "Eta vs Phi non-daughter;#eta;#phi", nBins, -eta_bin_lim, eta_bin_lim, nBins, -phi_bin_lim, phi_bin_lim);
    TH2F *resolution_vs_rrec = new TH2F("Resolution vs Rrec", "Resolution vs Rrec;Resolution;Rrec(cm)", nBins, -res_bin_lim, res_bin_lim, nBins, min_r, 50);

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
                    auto MCTrack = mcTracksMatrix[evID][trackID];
                    if (abs(MCTrack.GetPdgCode()) == hypPDG)
                    {
                        auto hypITSTrack = ITStracks->at(iTrack);
                        if (hypITSTrack.getNumberOfClusters() == 3) // 3 clusters recontruction doesn't work well
                            continue;

                        int firstDauID = MCTrack.getFirstDaughterTrackId();
                        int nDau = MCTrack.getLastDaughterTrackId();
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
                            continue; // if no triton daughter, improves speed

                        double pi0genPabs = 0;
                        for (int iDau = firstDauID; iDau < nDau; iDau++)
                        {
                            if (iDau == tritID)
                                continue;

                            if (abs(mcTracksMatrix[evID][iDau].GetPdgCode()) == pi0PDG)
                            {
                                cout << "pi0 found" << endl;
                                pi0genPabs = mcTracksMatrix[evID][iDau].GetP();
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
                                                    auto hypITSTrack = ft2.getTrack(1);
                                                    auto tritITSTPStrack = ft2.getTrack(0);

                                                    std::array<float, 3> hypP = {0, 0, 0};
                                                    std::array<float, 3> tritP = {0, 0, 0};
                                                    float hypPabs = 0;
                                                    float tritPabs = 0;
                                                    float etaHyp = 0;
                                                    float phiHyp = 0;
                                                    float etaTrit = 0;
                                                    float phiTrit = 0;

                                                    hypITSTrack.getPxPyPzGlo(hypP);
                                                    hypPabs = hypITSTrack.getP();
                                                    etaHyp = hypITSTrack.getEta();
                                                    phiHyp = hypITSTrack.getPhi();

                                                    tritITSTPStrack.getPxPyPzGlo(tritP);
                                                    tritPabs = tritITSTPStrack.getP();
                                                    etaTrit = tritITSTPStrack.getEta();
                                                    phiTrit = tritITSTPStrack.getPhi();

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
                                                    double p_res = (pi0genPabs - piPabs) / pi0genPabs;
                                                    if (pi0genPabs == 0)
                                                        p_res = -10;
                                                    pi0_p_resolution->Fill(p_res);
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

    auto c = new TCanvas("c", "c", 800, 600);
    nondaughter_chi_normalized->SetLineColor(kRed);
    nondaughter_chi_normalized->DrawNormalized("P");
    daughter_chi->DrawNormalized("sameP");
    c->Write();

    auto c1 = new TCanvas("c", "c", 800, 600);
    nondaughter_radius->SetLineColor(kRed);
    nondaughter_radius->SetTitle("Radius normalized");
    nondaughter_radius->DrawNormalized("P");
    daughter_radius->DrawNormalized("sameP");
    c1->Write();

    fFile.Close();
} // end of fitting function
