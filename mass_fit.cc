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
const int pi0PDG = 111;
const double tritonMass = 2.808921;
const double pi0Mass = 0.1349766;
const double hypMassTh = 2.99131;
TString chiLabel = "#chi^{2}";
TString hypLabel = "M_{^{3}_{#Lambda}H} (GeV/c^{2})";
int nBins = 100;
double min_bins = 0;
double min_r = 0;
double res_bin_lim = 0.25;
double eta_bin_lim = 0.1;
double phi_bin_lim = 0.1;

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

void mass_fit(TString path, TString filename, int tf_max = 40)
{
    bool cut = true;
    const int tf_min = 1;
    int tf_lenght = tf_max - tf_min + 1;
    if (cut)
    {
        min_r = 18;
        eta_bin_lim = 0.03;
        phi_bin_lim = 0.03;
    }

    int counts = 0;
    int b4chi = 0;
    int b4cut = 0;
    int b4p = 0;
    int final = 977;

    TH1F *chi_squared = new TH1F("chi2", "" + chiLabel + ";" + chiLabel + ";counts", nBins, min_bins, 1);
    TH1F *tot_chi_squared = new TH1F("tot_chi2", "Signal and Background " + chiLabel + ";" + chiLabel + ";counts", nBins, min_bins, 1);
    TH1F *bkg_chi_squared = new TH1F("bkg_chi2", "Background " + chiLabel + ";" + chiLabel + ";counts", nBins, min_bins, 50);
    TH1F *resolution = new TH1F("Radius Resolution", "Resolution;#Delta r;counts", nBins, -0.1, 0.1);
    TH1F *resolution_bkg = new TH1F("Radius Resolution Background", "Resolution;#Delta r;counts", nBins, -0.1, 0.1);
    TH1F *resolution_dl = new TH1F("Decay Length Resolution", "Resolution;(L_{gen} - L_{rec})/L_{gen};counts", nBins, -0.25, 0.25);
    TH1F *radius = new TH1F("Radius", "Radius;Rrec(cm);counts", nBins, min_r, 40);
    TH1F *inv_mass = new TH1F("Invariant mass", "Invariant mass;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *inv_mass_pi = new TH1F("Invariant mass pi", "#pi^{0} Invariant mass; m_{#pi^{0}};counts", nBins, 0, 0.3);
    TH1F *bkg_inv_mass = new TH1F("Background Invariant mass", "Background Invariant mass;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *tot_inv_mass = new TH1F("Total Invariant mass", "Invariant mass: Signal and Background;" + hypLabel + ";counts", nBins, 2.9, 4);

    TH1F *pi0_resolution = new TH1F("Pi0 p resolution", "#pi^{0} p resolution;Resolution;counts", nBins, -10, 10);
    TH1F *triton_resolution = new TH1F("Triton p resolution", "Triton p resolution;Resolution;counts", nBins, -5, 5);
    TH1F *hyp_resolution = new TH1F("Hyp p resolution", "Hyp p resolution;Resolution;counts", nBins, -10, 10);
    TH1F *pi0_resolution_normalized = new TH1F("Pi0 p resolution normalized", "#pi^{0} p resolution normalized;Resolution;counts", nBins, -10, 10);
    TH1F *triton_resolution_normalized = new TH1F("Triton p resolution normalized", "Triton p resolution normalized;Resolution;counts", nBins, -1, 1);
    TH1F *hyp_resolution_normalized = new TH1F("Hyp p resolution normalized", "Hyp p resolution normalized;Resolution;counts", nBins, -5, 5);

    TH1F *top_gen_pt = new TH1F("Hyp gen pt", "Hyp gen pt;p_{T} (GeV/c);counts", nBins, 0, 12);
    TH1F *top_rec_pt = new TH1F("Hyp rec pt", "Hyp rec pt;p_{T} (GeV/c);counts", nBins, 0, 12);
    TH1F *top_fitted_pt = new TH1F("Hyp fitted pt", "Hyp fitted pt;p_{T} (GeV/c);counts", nBins, 0, 12);
    TH1F *top_gen_r = new TH1F("Hyp gen r", "Hyp gen r;r (cm);counts", nBins, 15, 50);
    TH1F *top_rec_r = new TH1F("Hyp rec r", "Hyp rec r;r (cm);counts", nBins, 15, 50);
    TH1F *top_fitted_r = new TH1F("Hyp fitted r", "Hyp fitted r;r (cm);counts", nBins, 15, 50);

    TH2F *hyp_res_decay = new TH2F("Hyp p resolution vs decay radius", "Hyp p resolution vs decay radius;Radius (cm) ;Resolution (GeV/c)", nBins, 15, 40, nBins, -6, 6);
    TH2F *trit_res_decay = new TH2F("Triton p resolution vs decay radius", "Triton p resolution vs decay radius;Radius (cm) ;Resolution (GeV/c)", nBins, 15, 40, nBins, -2, 2);
    TH2F *hyp_res_layers = new TH2F("Hyp p resolution vs layers", "Hyp p resolution vs layers;Layers;Resolution (GeV/c)", 5, 2.5, 7.5, nBins, -6, 6);
    TH2F *trit_res_layers = new TH2F("Triton p resolution vs layers", "Triton p resolution vs layers;Layers;Resolution (GeV/c)", 5, 2.5, 7.5, nBins, -2, 2);
    TH2F *p_vs_e = new TH2F("p_vs_e", "(p_{rec} - p_{gen}) vs (E_{rec} - E_{gen}); (p_{rec} - p_{gen}) (GeV/c); (E_{rec} - E_{gen}) (GeV/c);counts", nBins, -1, 1, nBins, -1, 1);
    TH2F *mass_vs_p = new TH2F("mass_vs_p", "Mass vs p_{gen};p_{gen} (GeV/c);" + hypLabel + ";counts", nBins, 0, 16, nBins, 2.9, 3.7);

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
                            top_gen_pt->Fill(mcTrack.GetPt());
                            top_gen_r->Fill(rGen);
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

                        double pi0genPabs = 0;
                        TVector3 pi0genP;
                        for (int iDau = firstDauID; iDau <= nDau; iDau++)
                        {
                            if (abs(mcTracksMatrix[evID][iDau].GetPdgCode()) == pi0PDG)
                            {
                                pi0genPabs = mcTracksMatrix[evID][iDau].GetP();
                                mcTracksMatrix[evID][iDau].GetMomentum(pi0genP);
                            }
                        }

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

                                    bool isSignal = false;

                                    if (isDaughter && !tritfake && !fake)
                                        isSignal = true;

                                    if (isSignal)
                                        counts++;

                                    auto tritITSTPCtrack = ITSTPCtracks->at(jTrack);
                                    if (!tritfake && !fake)
                                    {
                                        top_rec_pt->Fill(MCTrack.GetPt());
                                        top_rec_r->Fill(genR);
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

                                                    if (isSignal)
                                                        b4chi++;

                                                    if (ft2.getChi2AtPCACandidate() < 0)
                                                        continue;

                                                    double chi2 = ft2.getChi2AtPCACandidate();

                                                    std::array<float, 3> R = ft2.getPCACandidatePos();
                                                    double recR = sqrt(R[0] * R[0] + R[1] * R[1]);
                                                    double genDL = calcDecayLenght(&mcTracksMatrix[evID], MCTrack, tritonPDG);
                                                    double recDL = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);

                                                    if (isSignal)
                                                        b4cut++;

                                                    if (cut)
                                                    {
                                                        if (recR < 18)
                                                            continue;
                                                        if (std::abs(etaHyp - etaTrit) > 0.03 || std::abs(phiHyp - phiTrit) > 0.03)
                                                            continue;
                                                    }

                                                    double res = (genR - recR) / genR;
                                                    double resDL = (genDL - recDL) / genDL;

                                                    float tritE = sqrt(tritPabs * tritPabs + tritonMass * tritonMass);
                                                    std::array<float, 3> piP = {hypP[0] - tritP[0], hypP[1] - tritP[1], hypP[2] - tritP[2]};
                                                    float piPabs = sqrt(piP[0] * piP[0] + piP[1] * piP[1] + piP[2] * piP[2]);
                                                    float piE = sqrt(pi0Mass * pi0Mass + piPabs * piPabs);
                                                    float hypE = piE + tritE;

                                                    float hypMass = sqrt(hypE * hypE - hypPabs * hypPabs);
                                                    tot_inv_mass->Fill(hypMass);
                                                    if (!isDaughter)
                                                    {
                                                        bkg_inv_mass->Fill(hypMass);
                                                        bkg_chi_squared->Fill(chi2);
                                                        tot_chi_squared->Fill(chi2);
                                                        resolution_bkg->Fill(res);
                                                        continue;
                                                    }

                                                    if (isSignal)
                                                        b4p++;

                                                    if (hypPabs < tritPabs)
                                                        continue;

                                                    double pt = sqrt(hypP[0] * hypP[0] + hypP[1] * hypP[1]);

                                                    top_fitted_pt->Fill(pt);
                                                    double rRec = calcRadius(&mcTracksMatrix[evID], MCTrack, tritonPDG);
                                                    top_fitted_r->Fill(rRec);

                                                    float hypEFound = sqrt(hypPabs * hypPabs + hypMassTh * hypMassTh);
                                                    float piEFound = hypEFound - tritE;
                                                    float piMassFound = sqrt(piEFound * piEFound - piPabs * piPabs);
                                                    if (piEFound * piEFound - piPabs * piPabs > 0)
                                                        inv_mass_pi->Fill(piMassFound);

                                                    chi_squared->Fill(chi2);
                                                    resolution->Fill(res);
                                                    resolution_dl->Fill(resDL);
                                                    radius->Fill(recR);
                                                    inv_mass->Fill(hypMass);

                                                    double pi0_p_res = (pi0genPabs - piPabs);
                                                    double trit_p_res = (tritgenPabs - tritPabs);
                                                    double hyp_p_res = (hypgenPabs - hypPabs);

                                                    pi0_resolution->Fill(pi0_p_res);
                                                    triton_resolution->Fill(trit_p_res);
                                                    hyp_resolution->Fill(hyp_p_res);

                                                    pi0_resolution_normalized->Fill(pi0_p_res / pi0genPabs);
                                                    triton_resolution_normalized->Fill(trit_p_res / tritgenPabs);
                                                    hyp_resolution_normalized->Fill(hyp_p_res / hypgenPabs);

                                                    hyp_res_decay->Fill(recR, hyp_p_res);
                                                    trit_res_decay->Fill(recR, trit_p_res);
                                                    hyp_res_layers->Fill(nLayers, hyp_p_res);
                                                    trit_res_layers->Fill(nLayers, trit_p_res);

                                                    p_vs_e->Fill(hypPabs - hypgenPabs, hypE - hypgenE);
                                                    mass_vs_p->Fill(hypgenPabs, hypMass);
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

    cout << "counts = " << counts << endl;
    cout << "b4chi = " << b4chi << endl;
    cout << "b4cut = " << b4cut << endl;
    cout << "b4p = " << b4p << endl;
    cout << "final = " << final << endl;
    /*
    efficiency = 0.010438
    counts = 1715
    startfit = 1715
    b4chi = 1681
    b4cut = 1679
    b4p = 1523
    final = 977
    */

    inv_mass->GetXaxis()->SetTitleSize(fontSize);
    inv_mass->GetYaxis()->SetTitleSize(fontSize);

    bkg_inv_mass->GetXaxis()->SetTitleSize(fontSize);
    bkg_inv_mass->GetYaxis()->SetTitleSize(fontSize);

    tot_inv_mass->GetXaxis()->SetTitleSize(fontSize);
    tot_inv_mass->GetYaxis()->SetTitleSize(fontSize);

    top_rec_r->GetXaxis()->SetTitleSize(fontSize);
    top_rec_r->GetYaxis()->SetTitleSize(fontSize);

    top_rec_pt->GetXaxis()->SetTitleSize(fontSize);
    top_rec_pt->GetYaxis()->SetTitleSize(fontSize);

    inv_mass_pi->GetXaxis()->SetTitleSize(fontSize);
    inv_mass_pi->GetYaxis()->SetTitleSize(fontSize);

    chi_squared->GetXaxis()->SetTitleSize(fontSize);
    chi_squared->GetYaxis()->SetTitleSize(fontSize);

    resolution->GetXaxis()->SetTitleSize(fontSize);
    resolution->GetYaxis()->SetTitleSize(fontSize);

    resolution_dl->GetXaxis()->SetTitleSize(fontSize);
    resolution_dl->GetYaxis()->SetTitleSize(fontSize);

    radius->GetXaxis()->SetTitleSize(fontSize);
    radius->GetYaxis()->SetTitleSize(fontSize);

    pi0_resolution->GetXaxis()->SetTitleSize(fontSize);
    pi0_resolution->GetYaxis()->SetTitleSize(fontSize);

    triton_resolution->GetXaxis()->SetTitleSize(fontSize);
    triton_resolution->GetYaxis()->SetTitleSize(fontSize);

    hyp_resolution->GetXaxis()->SetTitleSize(fontSize);
    hyp_resolution->GetYaxis()->SetTitleSize(fontSize);

    pi0_resolution_normalized->GetXaxis()->SetTitleSize(fontSize);
    pi0_resolution_normalized->GetYaxis()->SetTitleSize(fontSize);

    triton_resolution_normalized->GetXaxis()->SetTitleSize(fontSize);
    triton_resolution_normalized->GetYaxis()->SetTitleSize(fontSize);

    hyp_resolution_normalized->GetXaxis()->SetTitleSize(fontSize);
    hyp_resolution_normalized->GetYaxis()->SetTitleSize(fontSize);

    hyp_res_decay->GetXaxis()->SetTitleSize(fontSize);
    hyp_res_decay->GetYaxis()->SetTitleSize(fontSize);

    trit_res_decay->GetXaxis()->SetTitleSize(fontSize);
    trit_res_decay->GetYaxis()->SetTitleSize(fontSize);

    hyp_res_layers->GetXaxis()->SetTitleSize(fontSize);
    hyp_res_layers->GetYaxis()->SetTitleSize(fontSize);

    trit_res_layers->GetXaxis()->SetTitleSize(fontSize);
    trit_res_layers->GetYaxis()->SetTitleSize(fontSize);

    p_vs_e->GetXaxis()->SetTitleSize(fontSize);
    p_vs_e->GetYaxis()->SetTitleSize(fontSize);

    auto fFile = TFile(filename, "recreate");
    chi_squared->Write();
    bkg_chi_squared->Write();
    resolution->Write();
    resolution_dl->Write();
    radius->Write();

    inv_mass->Write();
    inv_mass_pi->Write();
    bkg_inv_mass->Write();
    tot_inv_mass->Write();

    pi0_resolution->Write();
    hyp_resolution->Write();
    triton_resolution->Write();

    pi0_resolution_normalized->Write();
    hyp_resolution_normalized->Write();
    triton_resolution_normalized->Write();

    hyp_res_decay->Write();
    trit_res_decay->Write();

    hyp_res_layers->Write();
    trit_res_layers->Write();

    top_gen_pt->Write();
    top_rec_pt->Write();
    top_fitted_pt->Write();
    top_gen_r->Write();
    top_rec_r->Write();
    top_fitted_r->Write();


    TH1F *eff_r = (TH1F *)top_rec_r->Clone("Hyp Eff r");
    eff_r->GetYaxis()->SetTitle("Efficiency");
    eff_r->SetTitle("Hypertriton r Efficiency");
    eff_r->Divide(top_gen_r);
    eff_r->Write();

    TH1F *eff_pt = (TH1F *)top_rec_pt->Clone("Hyp Eff pt");
    eff_pt->GetYaxis()->SetTitle("Efficiency");
    eff_pt->SetTitle("Hypertriton pt Efficiency");
    eff_pt->Divide(top_gen_pt);
    eff_pt->Write();

    TH1F *eff_fit_pt = (TH1F *)top_fitted_pt->Clone("Hyp Fit Eff pt");
    eff_fit_pt->GetYaxis()->SetTitle("Fitting Efficiency");
    eff_fit_pt->SetTitle("Hypertriton pt Fitting Efficiency");
    eff_fit_pt->Divide(top_rec_pt);
    eff_fit_pt->Write();

    TH1F *eff_fit_r = (TH1F *)top_fitted_r->Clone("Hyp Fit Eff r");
    eff_fit_r->GetYaxis()->SetTitle("Fitting Efficiency");
    eff_fit_r->SetTitle("Hypertriton r Fitting Efficiency");
    eff_fit_r->Divide(top_rec_r);
    eff_fit_r->Write();

    p_vs_e->Write();
    mass_vs_p->Write();

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    chi_squared->SetMarkerSize(markerSize);
    chi_squared->GetXaxis()->SetTitleSize(fontSize);
    chi_squared->GetYaxis()->SetTitleSize(fontSize);
    chi_squared->Draw("EP");
    tot_chi_squared->SetMarkerSize(markerSize);
    tot_chi_squared->GetXaxis()->SetTitleSize(fontSize);
    tot_chi_squared->GetYaxis()->SetTitleSize(fontSize);
    tot_chi_squared->SetLineColor(kRed);
    tot_chi_squared->Draw("sameEP");
    c1->Write();

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    inv_mass->SetMarkerSize(markerSize);
    inv_mass->GetXaxis()->SetTitleSize(fontSize);
    inv_mass->GetYaxis()->SetTitleSize(fontSize);
    inv_mass->Draw("EP");
    bkg_inv_mass->SetMarkerSize(markerSize);
    bkg_inv_mass->GetXaxis()->SetTitleSize(fontSize);
    bkg_inv_mass->GetYaxis()->SetTitleSize(fontSize);
    // bkg_inv_mass->SetLineColor(kRed);
    bkg_inv_mass->SetLineColor(kGreen);
    bkg_inv_mass->Draw("sameEP");
    c2->Write();

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    resolution->SetMarkerSize(markerSize);
    resolution->GetXaxis()->SetTitleSize(fontSize);
    resolution->GetYaxis()->SetTitleSize(fontSize);
    resolution->Draw("EP");
    resolution_bkg->SetMarkerSize(markerSize);
    resolution_bkg->GetXaxis()->SetTitleSize(fontSize);
    resolution_bkg->GetYaxis()->SetTitleSize(fontSize);
    resolution_bkg->SetLineColor(kRed);
    resolution_bkg->Draw("sameEP");
    c3->Write();

    fFile.Close();

} // end of fitting function
