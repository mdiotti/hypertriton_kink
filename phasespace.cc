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

const double hypPRange = 3;

// string FITTEROPTION = "DCA"; // "DCA_false" or "KFParticle"

void phasespace(TString filename, int nEvents = 100000, int seed = 0)
{
    gRandom->SetSeed(seed);

    if (!gROOT->GetClass("TGenPhaseSpace"))
        gSystem->Load("libPhysics");

    TH1F *h_piP = new TH1F("piP", "#pi^{0}  p;p (GeV/c);counts", 100, 0, 1.5);
    TH1F *h_tritonP = new TH1F("tritonP", "^{3}H p;p (GeV/c);counts", 100, 0, 6);

    TH1F *inv_mass = new TH1F("Invariant mass", "Invariant mass;" + hypLabel + ";counts", nBins, 2.8, 4.5);
    TH1F *inv_mass_pi = new TH1F("Invariant mass pi", "Invariant mass;M_{#pi^{0}};counts", nBins, 0, 0.3);

    TH2F *mass_vs_p = new TH2F("mass_vs_p", "Mass vs p;p_{gen} (GeV/c);" + hypLabel + ";counts", nBins, 0, 6, nBins, 2.8, 4.5);

    for (Int_t n = 0; n < nEvents; n++)
    {
        double px = gRandom->Uniform(-hypPRange, hypPRange);
        double py = gRandom->Uniform(-hypPRange, hypPRange);
        double pz = gRandom->Uniform(-hypPRange, hypPRange);
        double pGen = sqrt(px * px + py * py + pz * pz);
        double Energy = sqrt(pGen * pGen + hypMass * hypMass);

        TLorentzVector hypGen = TLorentzVector(px, py, pz, Energy);

        TGenPhaseSpace event;
        event.SetDecay(hypGen, 2, new Double_t[2]{piMass, tritonMass});

        Double_t weight = event.Generate();

        TLorentzVector *LorentzPi = event.GetDecay(0);
        TLorentzVector *LorentzTriton = event.GetDecay(1);

        std::array<double, 3> tritonP = {gRandom->Gaus(LorentzTriton->Px(), tritonSmearing * LorentzTriton->Px()),
                                         gRandom->Gaus(LorentzTriton->Py(), tritonSmearing * LorentzTriton->Py()),
                                         gRandom->Gaus(LorentzTriton->Pz(), tritonSmearing * LorentzTriton->Pz())};
        std::array<double, 3> hypP = {gRandom->Gaus(hypGen.Px(), hypSmearing * hypGen.Px()),
                                      gRandom->Gaus(hypGen.Py(), hypSmearing * hypGen.Py()),
                                      gRandom->Gaus(hypGen.Pz(), hypSmearing * hypGen.Pz())};
        std::array<double, 3> piP = {hypP[0] - tritonP[0], hypP[1] - tritonP[1], hypP[2] - tritonP[2]};

        double piPabs = sqrt(piP[0] * piP[0] + piP[1] * piP[1] + piP[2] * piP[2]);
        double tritonPabs = sqrt(tritonP[0] * tritonP[0] + tritonP[1] * tritonP[1] + tritonP[2] * tritonP[2]);
        double hypPabs = sqrt(hypP[0] * hypP[0] + hypP[1] * hypP[1] + hypP[2] * hypP[2]);

        double piE = sqrt(piPabs * piPabs + piMass * piMass);
        double tritonE = sqrt(tritonPabs * tritonPabs + tritonMass * tritonMass);
        double hypE = piE + tritonE;

        if (hypPabs < tritonPabs)
            continue;

        double hypRecM = sqrt(hypE * hypE - hypPabs * hypPabs);
        float hypEFound = sqrt(hypPabs * hypPabs + hypMass * hypMass);
        float piEFound = hypEFound - tritonE;
        float piMassFound = sqrt(piEFound * piEFound - piPabs * piPabs);

        if (piEFound * piEFound - piPabs * piPabs > 0)
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