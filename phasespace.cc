#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "TGenPhaseSpace.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRandom3.h"

#include "DetectorsVertexing/DCAFitterN.h"
#include "DataFormatsParameters/GRPObject.h"

#include "SimulationDataFormat/MCTrack.h"
#include "DataFormatsITS/TrackITS.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"

using namespace std;
using namespace o2;
using namespace vertexing;

const double hypMass = 2.99131;
const double piMass = 0.1349766;
const double tritonMass = 2.808921;

const double fontSize = 0.045;

const int hypPDG = 1010010030;
const int tritonPDG = 1000010030;
const int pi0PDG = 111;
TString chiLabel = "#chi^{2}";
TString hypLabel = "M_{^{3}_{#Lambda}H} (GeV/c^{2})";
int nBins = 100;
double min_bins = 0;

const double piSmearing = 0.05;     // 5%
const double tritonSmearing = 0.05; // 5%
const double hypSmearing = 0.2;     // 20%

const double hypPRange = 10;
const double etaRange = 1;

// string FITTEROPTION = "DCA"; // "DCA_false" or "KFParticle"

void phasespace(TString filename, int nEvents = 100000, int seed = 0)
{
    gRandom->SetSeed(seed);

    if (!gROOT->GetClass("TGenPhaseSpace"))
        gSystem->Load("libPhysics");

    TH1F *h_piP = new TH1F("piP", "#pi^{0}  p;p (GeV/c);counts", 100, 0, 1.5);
    TH1F *h_tritonP = new TH1F("tritonP", "^{3}H p;p (GeV/c);counts", 100, 0, 6);

    TH1F *inv_mass = new TH1F("Invariant mass", "Invariant mass;" + hypLabel + ";counts", nBins, 2.9, 4);
    TH1F *inv_mass_pi = new TH1F("Invariant mass pi", "#pi^{0} Invariant mass;M_{#pi^{0}};counts", nBins, 0, 0.05);

    TH1F *kink_angle = new TH1F("Kink angle", "Kink angle;#theta_{kink} (rad);counts", nBins, 0, 10);

    //TH1F *inv_mass_corrected = new TH1F("Invariant mass corrected", "Invariant mass corrected;" + hypLabel + ";counts", nBins, 2.9, 4);

    TH1F *pi0_resolution = new TH1F("Pi0 p resolution", "#pi^{0} p resolution;Resolution;counts", nBins, -10, 10);
    TH1F *triton_resolution = new TH1F("Triton p resolution", "Triton p resolution;Resolution;counts", nBins, -5, 5);
    TH1F *hyp_resolution = new TH1F("Hyp p resolution", "Hyp p resolution;Resolution;counts", nBins, -10, 10);

    TH2F *mass_vs_p = new TH2F("mass_vs_p", "Mass vs p_{gen};p_{gen} (GeV/c);" + hypLabel + ";counts", nBins, 0, 16, nBins, 2.9, 3.7);
    TH2F *p_vs_e = new TH2F("p_vs_e", "(p_{rec} - p_{gen}) vs (E_{rec} - E_{gen}); (p_{rec} - p_{gen}) (GeV/c); (E_{rec} - E_{gen}) (GeV/c);counts", nBins, -1, 1, nBins, -1, 1);

    for (Int_t n = 0; n < nEvents; n++)
    {
        double pT = gRandom->Uniform(0, hypPRange);
        double eta = gRandom->Uniform(-etaRange, etaRange);
        double pi = gRandom->Uniform(0, TMath::Pi() * 2);

        TLorentzVector hypGen = TLorentzVector(0, 0, 0, hypMass);
        hypGen.SetPtEtaPhiM(pT, eta, pi, hypMass);

        double pGen = hypGen.P();
        double EGen = hypGen.E();

        TGenPhaseSpace event;
        event.SetDecay(hypGen, 2, new Double_t[2]{piMass, tritonMass});

        Double_t weight = event.Generate();

        TLorentzVector *GenPi = event.GetDecay(0);
        TLorentzVector *LorentzTriton = event.GetDecay(1);

        double tritonPGen = LorentzTriton->P();
        double piPGen = GenPi->P();

        double HypPtSmeared = gRandom->Gaus(pT, hypSmearing * pT);
        double TritonPtSmeared = gRandom->Gaus(LorentzTriton->Pt(), tritonSmearing * LorentzTriton->Pt());

        LorentzTriton->SetPtEtaPhiM(TritonPtSmeared, LorentzTriton->Eta(), LorentzTriton->Phi(), tritonMass);
        double tritonPabs = LorentzTriton->P();

        TLorentzVector HypRec = TLorentzVector(0, 0, 0, hypMass);
        HypRec.SetPtEtaPhiM(HypPtSmeared, eta, pi, hypMass);
        double hypPabs = HypRec.P();

       // bool corrected = false;
/*
        while ((hypPabs * hypPabs + hypMass * hypMass) < (tritonPabs * tritonPabs + tritonMass * tritonMass))
        {
            double pTrec = HypRec.Pt();
            double deltapT = abs(gRandom->Gaus(0, hypSmearing * pTrec));
            double pTrecSmeared = pTrec + deltapT;
            HypRec.SetPtEtaPhiM(pTrecSmeared, HypRec.Eta(), HypRec.Phi(), hypMass);
            hypPabs = HypRec.P();
            corrected = true;
        }
*/

        TLorentzVector LorentzPi = TLorentzVector(HypRec.Px() - LorentzTriton->Px(), HypRec.Py() - LorentzTriton->Py(), HypRec.Pz() - LorentzTriton->Pz(), piMass);
        LorentzPi.SetPtEtaPhiM(LorentzPi.Pt(), LorentzPi.Eta(), LorentzPi.Phi(), piMass);
        double tritonE = LorentzTriton->E();
        double piPabs = LorentzPi.P();
        double piE = LorentzPi.E();

        double hypE = tritonE + piE;
        double hypRecM = sqrt(hypE * hypE - hypPabs * hypPabs);

        float hypEFound = sqrt(hypPabs * hypPabs + hypMass * hypMass);
        float piEFound = hypEFound - tritonE;
        float piMassFound = (piEFound * piEFound - piPabs * piPabs);
        if (piE * piE - piPabs * piPabs > 0)
            inv_mass_pi->Fill(piMassFound);

        double angle = HypRec.Angle(LorentzTriton->Vect())/TMath::Pi()*180;

        h_piP->Fill(piPabs);
        h_tritonP->Fill(tritonPabs);
        //if (!corrected)
           inv_mass->Fill(hypRecM);
        //else
        //   inv_mass_corrected->Fill(hypRecM);
        mass_vs_p->Fill(pGen, hypRecM);
        p_vs_e->Fill(hypPabs - pGen, hypE - EGen);
        pi0_resolution->Fill(piPGen - piPabs);
        triton_resolution->Fill(tritonPGen - tritonPabs);
        hyp_resolution->Fill(pGen - hypPabs);
        kink_angle->Fill(angle);
    }

    h_piP->GetXaxis()->SetTitleSize(fontSize);
    h_piP->GetYaxis()->SetTitleSize(fontSize);

    h_tritonP->GetXaxis()->SetTitleSize(fontSize);
    h_tritonP->GetYaxis()->SetTitleSize(fontSize);

    inv_mass->GetXaxis()->SetTitleSize(fontSize);
    inv_mass->GetYaxis()->SetTitleSize(fontSize);

    inv_mass_pi->GetXaxis()->SetTitleSize(fontSize);
    inv_mass_pi->GetYaxis()->SetTitleSize(fontSize);

    pi0_resolution->GetXaxis()->SetTitleSize(fontSize);
    pi0_resolution->GetYaxis()->SetTitleSize(fontSize);

    triton_resolution->GetXaxis()->SetTitleSize(fontSize);
    triton_resolution->GetYaxis()->SetTitleSize(fontSize);

    hyp_resolution->GetXaxis()->SetTitleSize(fontSize);
    hyp_resolution->GetYaxis()->SetTitleSize(fontSize);

    mass_vs_p->GetXaxis()->SetTitleSize(fontSize);
    mass_vs_p->GetYaxis()->SetTitleSize(fontSize);

    p_vs_e->GetXaxis()->SetTitleSize(fontSize);
    p_vs_e->GetYaxis()->SetTitleSize(fontSize);

    auto fFile = TFile(filename, "recreate");

    h_piP->Write();
    h_tritonP->Write();
    inv_mass->Write();
    //inv_mass_corrected->Write();
    inv_mass_pi->Write();
    pi0_resolution->Write();
    triton_resolution->Write();
    hyp_resolution->Write();
    mass_vs_p->Write();
    p_vs_e->Write();
    kink_angle->Write();

    double markerSize = 0.8;
/*
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    inv_mass->SetMarkerSize(markerSize);
    inv_mass->GetXaxis()->SetTitleSize(fontSize);
    inv_mass->GetYaxis()->SetTitleSize(fontSize);
    inv_mass->DrawNormalized("EP");
    inv_mass_corrected->SetMarkerSize(markerSize);
    inv_mass_corrected->GetXaxis()->SetTitleSize(fontSize);
    inv_mass_corrected->GetYaxis()->SetTitleSize(fontSize);
    inv_mass_corrected->SetLineColor(kGreen);
    inv_mass_corrected->DrawNormalized("sameEP");
    c1->Write();
*/
    fFile.Close();
}