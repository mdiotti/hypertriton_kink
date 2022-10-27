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

    TH1F *inv_mass = new TH1F("Invariant mass", "Invariant mass;" + hypLabel + ";counts", nBins, 2.9, 3.7);
    TH1F *inv_mass_pi = new TH1F("Invariant mass pi", "Invariant mass;M_{#pi^{0}};counts", nBins, 0, 0.3);

    TH2F *mass_vs_p = new TH2F("mass_vs_p", "Mass vs p;p_{gen} (GeV/c);" + hypLabel + ";counts", nBins, 0, 16, nBins, 2.9, 3.7);

    for (Int_t n = 0; n < nEvents; n++)
    {
        double pT = gRandom->Uniform(0, hypPRange);
        double eta = gRandom->Uniform(-etaRange, etaRange);
        double pi = gRandom->Uniform(0, TMath::Pi() * 2);

        TLorentzVector hypGen = TLorentzVector(0, 0, 0, hypMass);
        hypGen.SetPtEtaPhiM(pT, eta, pi, hypMass);

        // double Energy = hypGen.E();
        double pGen = hypGen.P();

        TGenPhaseSpace event;
        event.SetDecay(hypGen, 2, new Double_t[2]{piMass, tritonMass});

        Double_t weight = event.Generate();

        // TLorentzVector *LorentzPi = event.GetDecay(0);
        TLorentzVector *LorentzTriton = event.GetDecay(1);

        double HypPtSmeared = gRandom->Gaus(pT, hypSmearing * pT);
        double TritonPtSmeared = gRandom->Gaus(LorentzTriton->Pt(), tritonSmearing * LorentzTriton->Pt());

        LorentzTriton->SetPtEtaPhiM(TritonPtSmeared, LorentzTriton->Eta(), LorentzTriton->Phi(), tritonMass);

        TLorentzVector HypRec = TLorentzVector(0, 0, 0, hypMass);
        HypRec.SetPtEtaPhiM(HypPtSmeared, eta, pi, hypMass);

        TLorentzVector LorentzPi = TLorentzVector(HypRec.Px() - LorentzTriton->Px(), HypRec.Py() - LorentzTriton->Py(), HypRec.Pz() - LorentzTriton->Pz(), piMass);
        LorentzPi.SetPtEtaPhiM(LorentzPi.Pt(), LorentzPi.Eta(), LorentzPi.Phi(), piMass);
        double tritonPabs = LorentzTriton->P();
        double tritonE = LorentzTriton->E();
        double piPabs = LorentzPi.P();
        double piE = LorentzPi.E();

        double hypPabs = HypRec.P();
        double hypE = tritonE + piE;

        double hypRecM = sqrt(hypE * hypE - hypPabs * hypPabs);

        if (hypPabs < tritonPabs)
            continue;

        float hypEFound = sqrt(hypPabs * hypPabs + hypMass * hypMass);
        float piEFound = hypEFound - tritonE;
        float piMassFound = (piEFound * piEFound - piPabs * piPabs);
        if (piE * piE - piPabs * piPabs > 0)
            inv_mass_pi->Fill(piMassFound);

        h_piP->Fill(piPabs);
        h_tritonP->Fill(tritonPabs);
        inv_mass->Fill(hypRecM);
        mass_vs_p->Fill(pGen, hypRecM);
    }

    auto fFile = TFile(filename, "recreate");

    h_piP->Write();
    h_tritonP->Write();
    inv_mass->Write();
    inv_mass_pi->Write();
    mass_vs_p->Write();

    fFile.Close();
}