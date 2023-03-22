import uproot
import pandas as pd
import ROOT

def fill_th1_hist(df, var,h):
    for var_val in df[var]:
        h.Fill(var_val)


kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC  = ROOT.TColor.GetColor("#ff7f00")

mc_mother_pt = ROOT.TH1F("mc_mother_pt", "mc_mother_pt", 100, 0, 12)
mc_daughter_pt = ROOT.TH1F("mc_daughter_pt", "mc_daughter_pt", 100, 0, 12)
mc_dl = ROOT.TH1F("mc_dl", "mc_dl", 100, 0, 50)
mc_radius = ROOT.TH1F("mc_radius", "mc_radius", 100, 0, 50)

mother_pt = ROOT.TH1F("mother_pt", "mother_pt", 100, 0, 12)
daughter_pt = ROOT.TH1F("daughter_pt", "daughter_pt", 100, 0, 12)
dl = ROOT.TH1F("dl", "dl", 100, 0, 50)
radius = ROOT.TH1F("radius", "radius", 100, 0, 50)
mass_rec = ROOT.TH1F("mass", "mass", 100, 2.9, 4)
chi2 = ROOT.TH1F("chi2", "chi2", 100, 0, 100)

top_mother_pt = ROOT.TH1F("top_mother_pt", "top_mother_pt", 100, 0, 12)
top_daughter_pt = ROOT.TH1F("top_daughter_pt", "top_daughter_pt", 100, 0, 12)
top_dl = ROOT.TH1F("top_dl", "top_dl", 100, 0, 50)
top_radius = ROOT.TH1F("top_radius", "top_radius", 100, 0, 50)
top_mass = ROOT.TH1F("top_mass", "top_mass", 100, 2.9, 4)
top_chi2 = ROOT.TH1F("top_chi2", "top_chi2", 100, 0, 100)

top_gen_pt = ROOT.TH1F("top_gen_pt", "top_gen_pt", 100, 0, 12)

df_mc = uproot.open("TrackedKinkTree.root")["MCTree"].arrays(library="pd")

df_rec = uproot.open("TrackedKinkTree.root")["KinkTree"].arrays(library="pd").query("isHyp == 1").query("gNLayers > 3")

df_rec.drop("gVertex[0]",1)
df_rec.drop("gVertex[1]",1)
df_rec.drop("gVertex[2]",1)

df_topology = df_rec.query("isTopology == 1")

fill_th1_hist(df_mc, "mcMotherPt", mc_mother_pt)
fill_th1_hist(df_mc, "mcDaughterPt", mc_daughter_pt)
fill_th1_hist(df_mc, "mcDecayLength", mc_dl)
fill_th1_hist(df_mc, "mcRadius", mc_radius)

fill_th1_hist(df_rec, "gMotherPt", mother_pt)
fill_th1_hist(df_rec, "gDaughterPt", daughter_pt)
fill_th1_hist(df_rec, "gDecayLength", dl)
fill_th1_hist(df_rec, "gRadius", radius)
fill_th1_hist(df_rec, "gMass", mass_rec)
fill_th1_hist(df_rec, "gChi2", chi2)

fill_th1_hist(df_topology, "gMotherPt", top_mother_pt)
fill_th1_hist(df_topology, "gDaughterPt", top_daughter_pt)
fill_th1_hist(df_topology, "gDecayLength", top_dl)
fill_th1_hist(df_topology, "gRadius", top_radius)
fill_th1_hist(df_topology, "gMass", top_mass)
fill_th1_hist(df_topology, "gChi2", top_chi2)
fill_th1_hist(df_topology, "gGeneratedMotherPt", top_gen_pt)

### sum w2 for all the histograms

if (0):
    mc_mother_pt.Sumw2()
    mc_daughter_pt.Sumw2()
    mc_dl.Sumw2()
    mc_radius.Sumw2()

    mother_pt.Sumw2()
    daughter_pt.Sumw2()
    dl.Sumw2()
    radius.Sumw2()
    mass_rec.Sumw2()

    top_mother_pt.Sumw2()
    top_daughter_pt.Sumw2()
    top_dl.Sumw2()
    top_radius.Sumw2()
    top_mass.Sumw2()

#mass_rec.SetStats(0)
#top_mass.SetStats(0)
mass_rec.GetYaxis().SetTitle("Counts")
mass_rec.GetXaxis().SetTitle("#it{M}(^{3}He + #pi^{-} and c.c.)   (GeV/#it{c}^{2})")
top_mass.GetYaxis().SetTitle("Counts")
top_mass.GetXaxis().SetTitle("#it{M}(^{3}He + #pi^{-} and c.c.)   (GeV/#it{c}^{2})")

if(0):
    mass = ROOT.RooRealVar('m', '#it{M}(^{3}He + #pi^{-} and c.c.)', 2.96, 3.04, 'GeV/c^{2}')
    mu = ROOT.RooRealVar('mu', 'hypernucl mass', 2.96, 3.04, 'GeV/c^{2}')
    sigma = ROOT.RooRealVar('sigma', 'hypernucl width', 0.001, 0.004, 'GeV/c^{2}')
    a1 = ROOT.RooRealVar('a1', 'a1', 0, 3.)
    a2 = ROOT.RooRealVar('a2', 'a2', 0.3, 0.8)

    c0 = ROOT.RooRealVar('c0', 'constant c0', -100., 100)
    c1 = ROOT.RooRealVar('c1', 'constant c1', -100., 100)
    c2 = ROOT.RooRealVar('c2', 'constant c2', -100., 100)

    n1 = ROOT.RooRealVar('n1', 'n1', 1, 10.)
    n2 = ROOT.RooRealVar('n2', 'n2', 1, 10.)
    signal = ROOT.RooDSCBShape('cb', 'cb', mass, mu, sigma, a1, n1, a2, n2)

    background1 = ROOT.RooChebychev('bkg', 'pol3 bkg', mass, ROOT.RooArgList(c0,c1, c2))
    background2 = ROOT.RooChebychev('bkg', 'pol1 bkg', mass, ROOT.RooArgList(c0))

    n = ROOT.RooRealVar('n', 'n const', 0.01, 1)
    # define the fit funciton and perform the actual fit
    fit_function1 = ROOT.RooAddPdf('total_pdf1', 'signal + background 1', ROOT.RooArgList(signal, background1), ROOT.RooArgList(n))
    fit_function2 = ROOT.RooAddPdf('total_pdf2', 'signal + background 2', ROOT.RooArgList(signal, background2), ROOT.RooArgList(n))


outfile = ROOT.TFile("histos.root", "recreate")

#_, _, frame1 = fit_and_plot(mass_v0, mass, fit_function1, signal, background1, sigma, mu, n)
#_, _, frame2 = fit_and_plot(mass_tracked, mass, fit_function2, signal, background2, sigma, mu, n)


## rescale background and signal
# currently, one signal event per bkg is generated. to be rescaled: 0.8*1e-7
yield_signal = 0.8*1e-7

 
#bkg_dataset = background1.generate(ROOT.RooArgSet(mass), 1e6)
#signal_dataset = signal.generate(ROOT.RooArgSet(mass), 1e3)
#scaled_bkg_histo = bkg_dataset.createHistogram("scaled_bkg_histo", mass, ROOT.RooFit.Binning(100, 2.96, 3.04))
#scaled_signal_histo = signal_dataset.createHistogram("scaled_signal_histo", mass, ROOT.RooFit.Binning(100, 2.96, 3.04))

## sum the scaled background and signal
#sum_histo = scaled_bkg_histo.Clone("sum_histo")
#sum_histo.Add(scaled_signal_histo)
#sum_histo.Write()

mc_mother_pt.Write()
mc_daughter_pt.Write()
mc_dl.Write()
mc_radius.Write()

mother_pt.Write()
daughter_pt.Write()
dl.Write()
radius.Write()
chi2.Write()

top_mother_pt.Write()
top_daughter_pt.Write()
top_dl.Write()
top_radius.Write()
top_chi2.Write()



if(0):
    cv1 = ROOT.TCanvas()
    #frame1.Draw()
    leg1 = ROOT.TLegend(0.52,0.52,0.93,0.76)
    leg1.AddEntry("data","{}_{#Lambda}^{3}H + {}_{#bar{#Lambda}}^{3}#bar{H}", "PE")
    leg1.AddEntry("fit_func","Signal + Background", "L")
    leg1.AddEntry("bkg","Background", "L")
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(42)
    leg1.Draw()
    cv1.Write()

if(0):
    cv2 = ROOT.TCanvas()
    #frame2.Draw()
    cv2.Write()

mass_rec.Write()
top_mass.Write()

if(0):
    pinfo = ROOT.TPaveText(0.7,0.7,0.9,0.9,"NDC")
    pinfo.AddText("This thesis, MC simulation")
    pinfo.AddText("pp, #sqrt{#it{s}_{NN}} = 13.6 TeV")
    pinfo.SetBorderSize(0)
    pinfo.SetFillStyle(0)
    pinfo.SetTextAlign(11)
    pinfo.SetTextFont(42)


## efficiency plots
eff_top = top_gen_pt.Clone("eff_top_pt")
eff_top.Divide(mc_mother_pt)
eff_top.SetStats(0)
eff_top.GetYaxis().SetTitle("Topology efficiency")
eff_top.GetXaxis().SetTitle("#it{p}_{T}   (GeV/#it{c})")
eff_top.SetTitle("")
eff_top.Write()


eff_radius_top = top_radius.Clone("eff_radius_top")
eff_radius_top.Divide(mc_radius)
eff_radius_top.SetStats(0)
eff_radius_top.GetYaxis().SetTitle("Topology efficiency")
eff_radius_top.GetXaxis().SetTitle("#it{R}   (cm)")
eff_radius_top.SetTitle("")
eff_radius_top.Write()
