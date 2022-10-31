
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

const double fontSize = 0.04;

double lim1 = 1;
TString limLabel1 = Form("%f", lim1);

double lim2 = 2;
TString limLabel2 = Form("%f", lim2);

double lim3 = 4;
TString limLabel3 = Form("%f", lim3);

double lim4 = 6;
TString limLabel4 = Form("%f", lim4);

double lim5 = 8;
TString limLabel5 = Form("%f", lim5);

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

void mass_fit(TString path, TString filename, int tf_max = 40, bool partial = false)
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

    TH1F *chi_squared = new TH1F("chi2", chiLabel + ";" + chiLabel + ";counts", nBins, min_bins, 1);
    TH1F *bkg_chi_squared = new TH1F("bkg_chi2", "Background " + chiLabel + ";" + chiLabel + ";counts", nBins, min_bins, 1);
    TH1F *resolution = new TH1F("Radius Resolution", "Resolution;#Delta r;counts", nBins, -res_bin_lim, res_bin_lim);
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
    TH1F *hyp_gen_p = new TH1F("Hyp gen p", "Hyp gen p;p (GeV/c);counts", nBins, 0, 10);
    TH1F *hyp_rec_p = new TH1F("Hyp rec p", "Hyp rec p;p (GeV/c);counts", nBins, 0, 10);
    TH1F *hyp_gen_r = new TH1F("Hyp gen r", "Hyp gen r;r (cm);counts", nBins, 15, 40);
    TH1F *hyp_rec_r = new TH1F("Hyp rec r", "Hyp rec r;r (cm);counts", nBins, 15, 40);

    TH1F *hyp_rel_lim1 = new TH1F("Hyp p relative resolution " + limLabel1, "Hyp p relative resolution when #pi^{0} resolution is <-" + limLabel1 + ";Resolution;counts", nBins, -1.5, 1.5);
    TH1F *hyp_rel_lim2 = new TH1F("Hyp p relative resolution " + limLabel2, "Hyp p relative resolution when #pi^{0} resolution is <-" + limLabel2 + ";Resolution;counts", nBins, -1.5, 1.5);
    TH1F *hyp_rel_lim3 = new TH1F("Hyp p relative resolution " + limLabel3, "Hyp p relative resolution when #pi^{0} resolution is <-" + limLabel3 + ";Resolution;counts", nBins, -1.5, 1.5);
    TH1F *hyp_rel_lim4 = new TH1F("Hyp p relative resolution " + limLabel4, "Hyp p relative resolution when #pi^{0} resolution is <-" + limLabel4 + ";Resolution;counts", nBins, -1.5, 1.5);
    TH1F *hyp_rel_lim5 = new TH1F("Hyp p relative resolution " + limLabel5, "Hyp p relative resolution when #pi^{0} resolution is <-" + limLabel5 + ";Resolution;counts", nBins, -1.5, 1.5);

    TH1F *hyp_lim1 = new TH1F("Hyp p resolution " + limLabel1, "Hyp p resolution when #pi^{0} resolution is <-" + limLabel1 + ";Resolution;counts", nBins, -15, 15);
    TH1F *hyp_lim2 = new TH1F("Hyp p resolution " + limLabel2, "Hyp p resolution when #pi^{0} resolution is <-" + limLabel2 + ";Resolution;counts", nBins, -15, 15);
    TH1F *hyp_lim3 = new TH1F("Hyp p resolution " + limLabel3, "Hyp p resolution when #pi^{0} resolution is <-" + limLabel3 + ";Resolution;counts", nBins, -15, 15);
    TH1F *hyp_lim4 = new TH1F("Hyp p resolution " + limLabel4, "Hyp p resolution when #pi^{0} resolution is <-" + limLabel4 + ";Resolution;counts", nBins, -15, 15);
    TH1F *hyp_lim5 = new TH1F("Hyp p resolution " + limLabel5, "Hyp p resolution when #pi^{0} resolution is <-" + limLabel5 + ";Resolution;counts", nBins, -15, 15);

    TH1F *trit_lim1 = new TH1F("Triton p resolution " + limLabel1, "Triton p resolution when #pi^{0} resolution is <-" + limLabel1 + ";Resolution;counts", nBins, -5, 5);
    TH1F *trit_lim2 = new TH1F("Triton p resolution " + limLabel2, "Triton p resolution when #pi^{0} resolution is <-" + limLabel2 + ";Resolution;counts", nBins, -5, 5);
    TH1F *trit_lim3 = new TH1F("Triton p resolution " + limLabel3, "Triton p resolution when #pi^{0} resolution is <-" + limLabel3 + ";Resolution;counts", nBins, -5, 5);
    TH1F *trit_lim4 = new TH1F("Triton p resolution " + limLabel4, "Triton p resolution when #pi^{0} resolution is <-" + limLabel4 + ";Resolution;counts", nBins, -5, 5);
    TH1F *trit_lim5 = new TH1F("Triton p resolution " + limLabel5, "Triton p resolution when #pi^{0} resolution is <-" + limLabel5 + ";Resolution;counts", nBins, -5, 5);

    TH1F *pi0_lim1 = new TH1F("Pi0 p resolution " + limLabel1, "#pi^{0} p resolution when #pi^{0} resolution is <-" + limLabel1 + ";Resolution;counts", nBins, 0, 25);
    TH1F *pi0_lim2 = new TH1F("Pi0 p resolution " + limLabel2, "#pi^{0} p resolution when #pi^{0} resolution is <-" + limLabel2 + ";Resolution;counts", nBins, 0, 25);
    TH1F *pi0_lim3 = new TH1F("Pi0 p resolution " + limLabel3, "#pi^{0} p resolution when #pi^{0} resolution is <-" + limLabel3 + ";Resolution;counts", nBins, 0, 25);
    TH1F *pi0_lim4 = new TH1F("Pi0 p resolution " + limLabel4, "#pi^{0} p resolution when #pi^{0} resolution is <-" + limLabel4 + ";Resolution;counts", nBins, 0, 25);
    TH1F *pi0_lim5 = new TH1F("Pi0 p resolution " + limLabel5, "#pi^{0} p resolution when #pi^{0} resolution is <-" + limLabel5 + ";Resolution;counts", nBins, 0, 25);

    TH2F *hyp_res_decay = new TH2F("Hyp p resolution vs decay radius", "Hyp p resolution vs decay radius;Radius (cm) ;Resolution (GeV/c)", nBins, 15, 40, nBins, -6, 6);
    TH2F *trit_res_decay = new TH2F("Triton p resolution vs decay radius", "Triton p resolution vs decay radius;Radius (cm) ;Resolution (GeV/c)", nBins, 15, 40, nBins, -2, 2);
    TH2F *hyp_res_layers = new TH2F("Hyp p resolution vs layers", "Hyp p resolution vs layers;Layers;Resolution (GeV/c)", 5, 2.5, 7.5, nBins, -6, 6);
    TH2F *trit_res_layers = new TH2F("Triton p resolution vs layers", "Triton p resolution vs layers;Layers;Resolution (GeV/c)", 5, 2.5, 7.5, nBins, -2, 2);
    TH2F *p_vs_e = new TH2F("p_vs_e", "(p_{rec} - p_{gen}) vs (E_{rec} - E_{gen}); (p_{rec} - p_{gen}) (GeV/c); (E_{rec} - E_{gen}) (GeV/c);counts", nBins, -1, 1, nBins, -1, 1);
    TH2F *mass_vs_p = new TH2F("mass_vs_p", "Mass vs p_{gen};p_{gen} (GeV/c);" + hypLabel + ";counts", nBins, 0, 16, nBins, 2.9, 3.7);

    double genEntries = 0;

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

        for (int n = 0; n < nev; n++) // fill Gen histos
        {
            treeMCTracks->GetEvent(n);
            unsigned int size = MCtracks->size();

            for (unsigned int mcI{0}; mcI < size; ++mcI)
            {
                auto mcTrack = MCtracks->at(mcI);
                hyp_gen_p->Fill(mcTrack.GetP());
                double rGen = calcRadius(&mcTracksMatrix[n], mcTrack, tritonPDG);
                hyp_gen_r->Fill(rGen);
                genEntries++;
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

                                    //if (!isDaughter)
                                    //    continue;

                                    auto tritITSTPCtrack = ITSTPCtracks->at(jTrack);
                                    if (!tritfake && !fake)
                                    {
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

                                                    tot_inv_mass->Fill(hypMass);
                                                    if (!isDaughter){
                                                        bkg_inv_mass->Fill(hypMass);
                                                        bkg_chi_squared->Fill(chi2);
                                                        continue;
                                                    }

                                                    if (hypPabs < tritPabs)
                                                        continue;

                                                    hyp_rec_p->Fill(hypPabs);
                                                    double rRec = calcRadius(&mcTracksMatrix[evID], MCTrack, tritonPDG);
                                                    hyp_rec_r->Fill(rRec);

                                                    float hypEFound = sqrt(hypPabs * hypPabs + hypMassTh * hypMassTh);
                                                    float piEFound = hypEFound - tritE;
                                                    float piMassFound = sqrt(piEFound * piEFound - piPabs * piPabs);
                                                    if (piEFound * piEFound - piPabs * piPabs > 0)
                                                        inv_mass_pi->Fill(piMassFound);

                                                    chi_squared->Fill(chi2);
                                                    resolution->Fill(res);
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

                                                    if (partial)
                                                    {

                                                        if (pi0_p_res < -lim1)
                                                        {
                                                            hyp_lim1->Fill(hyp_p_res);
                                                            trit_lim1->Fill(trit_p_res);
                                                            pi0_lim1->Fill(pi0_p_res);
                                                            hyp_rel_lim1->Fill(hyp_p_res / hypgenPabs);
                                                        }
                                                        if (pi0_p_res < -lim2)
                                                        {
                                                            hyp_lim2->Fill(hyp_p_res);
                                                            trit_lim2->Fill(trit_p_res);
                                                            pi0_lim2->Fill(pi0_p_res);
                                                            hyp_rel_lim2->Fill(hyp_p_res / hypgenPabs);
                                                        }
                                                        if (pi0_p_res < -lim3)
                                                        {
                                                            hyp_lim3->Fill(hyp_p_res);
                                                            trit_lim3->Fill(trit_p_res);
                                                            pi0_lim3->Fill(pi0_p_res);
                                                            hyp_rel_lim3->Fill(hyp_p_res / hypgenPabs);
                                                        }
                                                        if (pi0_p_res < -lim4)
                                                        {
                                                            hyp_lim4->Fill(hyp_p_res);
                                                            trit_lim4->Fill(trit_p_res);
                                                            pi0_lim4->Fill(pi0_p_res);
                                                            hyp_rel_lim4->Fill(hyp_p_res / hypgenPabs);
                                                        }
                                                        if (pi0_p_res < -lim5)
                                                        {
                                                            hyp_lim5->Fill(hyp_p_res);
                                                            trit_lim5->Fill(trit_p_res);
                                                            pi0_lim5->Fill(pi0_p_res);
                                                            hyp_rel_lim5->Fill(hyp_p_res / hypgenPabs);
                                                        }
                                                    }
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

    double recEntries = hyp_rec_p->GetEntries();
    double eff = recEntries / genEntries;

    cout << "efficiency = " << eff << endl;
    // efficiency = 1.56696e-06

    inv_mass->GetXaxis()->SetTitleSize(fontSize);
    inv_mass->GetYaxis()->SetTitleSize(fontSize);
    bkg_inv_mass->GetXaxis()->SetTitleSize(fontSize);
    bkg_inv_mass->GetYaxis()->SetTitleSize(fontSize);
    tot_inv_mass->GetXaxis()->SetTitleSize(fontSize);
    tot_inv_mass->GetYaxis()->SetTitleSize(fontSize);

    auto fFile = TFile(filename, "recreate");
    chi_squared->Write();
    bkg_chi_squared->Write();
    resolution->Write();
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

    TH1F *eff_p = (TH1F *)hyp_rec_p->Clone("Hyp Eff p");
    eff_p->GetYaxis()->SetTitle("Efficiency");
    eff_p->SetTitle("Hypertriton p Efficiency");
    eff_p->Divide(hyp_gen_p);
    eff_p->Write();

    TH1F *eff_r = (TH1F *)hyp_rec_r->Clone("Hyp Eff r");
    eff_r->GetYaxis()->SetTitle("Efficiency");
    eff_r->SetTitle("Hypertriton r Efficiency");
    eff_r->Divide(hyp_gen_r);
    eff_r->Write();

    p_vs_e->Write();
    mass_vs_p->Write();

    fFile.Close();

    if (partial)
    {
        auto fFile2 = TFile("partial_hyp.root", "recreate");
        hyp_lim1->Write();
        hyp_lim2->Write();
        hyp_lim3->Write();
        hyp_lim4->Write();
        hyp_lim5->Write();
        fFile2.Close();

        auto fFile3 = TFile("partial_trit.root", "recreate");
        trit_lim1->Write();
        trit_lim2->Write();
        trit_lim3->Write();
        trit_lim4->Write();
        trit_lim5->Write();
        fFile3.Close();

        auto fFile4 = TFile("partial_pi0.root", "recreate");
        pi0_lim1->Write();
        pi0_lim2->Write();
        pi0_lim3->Write();
        pi0_lim4->Write();
        pi0_lim5->Write();
        fFile4.Close();

        auto fFile5 = TFile("partial_hyp_rel.root", "recreate");
        hyp_rel_lim1->Write();
        hyp_rel_lim2->Write();
        hyp_rel_lim3->Write();
        hyp_rel_lim4->Write();
        hyp_rel_lim5->Write();
        fFile5.Close();
    }

} // end of fitting function
