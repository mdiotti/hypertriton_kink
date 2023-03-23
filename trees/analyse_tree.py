import uproot
import pandas as pd
import ROOT

def fill_th1_hist(df, var,h):
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th2_hist(df, var1, var2, h):
    for ind in range(len(df)):
        h.Fill(df.iloc[ind][var1], df.iloc[ind][var2])


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
mass = ROOT.TH1F("mass", "mass", 100, 2.9, 4)
chi2 = ROOT.TH1F("chi2Matching", "chi2Matching", 100, 0, 100000)
chi2DCA = ROOT.TH1F("chi2DCA", "chi2DCA", 100, 0, 0.5)

chi2_bkg = ROOT.TH1F("chi2Matching_bkg", "chi2Matching_bkg", 100, 0, 100000)
chi2DCA_bkg = ROOT.TH1F("chi2DCA_bkg", "chi2DCA_bkg", 100, 0, 0.5)

radius_bkg = ROOT.TH1F("radius_bkg", "radius_bkg", 100, 0, 50)
mass_bkg = ROOT.TH1F("mass_bkg", "mass_bkg", 100, 2.9, 4)

gen_pt = ROOT.TH1F("gen_pt", "gen_pt", 100, 0, 12)
gen_r = ROOT.TH1F("gen_r", "gen_r", 100, 0, 50)

res_vs_pgen = ROOT.TH2F("res_vs_pgen", "res_vs_pgen", 100, 0, 12, 100, -3, 3)
res_vs_pgen_dau = ROOT.TH2F("res_vs_pgen_dau", "res_vs_pgen_dau", 100, 0, 12, 100, -3, 3)
res_vs_radius = ROOT.TH2F("res_vs_radius", "res_vs_radius", 100, 0, 50, 100, -3, 3)

df_mc = uproot.open("TrackedKinkTree.root")["MCTree"].arrays(library="pd")

df_rec = uproot.open("TrackedKinkTree.root")["KinkTree"].arrays(library="pd").query("NLayers > 3").query("motherFake == 0").query("daughterFake == 0")
df_topology = df_rec.query("isTopology == 1").query("isHyp == 1")
df_bkg = df_rec.query("isTopology == 0")


res = ROOT.TH1F("res_mother", "res_mother", 100, -3, 3)
res_dau = ROOT.TH1F("res_dau", "res_dau", 100, -3, 3)
radius_res = ROOT.TH1F("radius_res", "radius_res", 100, -3, 3)

df_topology["res"] = (df_topology["rMotherPt"] - df_topology["gMotherPt"])/df_topology["gMotherPt"]
df_topology["res_dau"] = (df_topology["rDaughterPt"] - df_topology["gDaughterPt"])/df_topology["gDaughterPt"]
df_topology["radius_res"] = (df_topology["rRadius"] - df_topology["gRadius"])/df_topology["gRadius"]

fill_th1_hist(df_topology, "res", res)
fill_th1_hist(df_topology, "res_dau", res_dau)
fill_th1_hist(df_topology, "radius_res", radius_res)

fill_th1_hist(df_mc, "mcMotherPt", mc_mother_pt)
fill_th1_hist(df_mc, "mcDaughterPt", mc_daughter_pt)
fill_th1_hist(df_mc, "mcDecayLength", mc_dl)
fill_th1_hist(df_mc, "mcRadius", mc_radius)

fill_th1_hist(df_topology, "rMotherPt", mother_pt)
fill_th1_hist(df_topology, "rDaughterPt", daughter_pt)
fill_th1_hist(df_topology, "rDecayLength", dl)
fill_th1_hist(df_topology, "rRadius", radius)
fill_th1_hist(df_topology, "Mass", mass)
fill_th1_hist(df_topology, "Chi2Match", chi2)
fill_th1_hist(df_topology, "Chi2DCA", chi2DCA)
fill_th1_hist(df_topology, "gMotherPt", gen_pt)
fill_th1_hist(df_topology, "gRadius", gen_r)

fill_th1_hist(df_bkg, "Chi2Match", chi2_bkg)
fill_th1_hist(df_bkg, "Chi2DCA", chi2DCA_bkg)
fill_th1_hist(df_bkg, "rRadius", radius_bkg)
fill_th1_hist(df_bkg, "Mass", mass_bkg)

fill_th2_hist(df_topology, "gMotherPt", "res", res_vs_pgen)
fill_th2_hist(df_topology, "gDaughterPt", "res_dau", res_vs_pgen_dau)
fill_th2_hist(df_topology, "gRadius", "radius_res", res_vs_radius)


### sum w2 for all the histograms

if (False):
    mc_mother_pt.Sumw2()
    mc_daughter_pt.Sumw2()
    mc_dl.Sumw2()
    mc_radius.Sumw2()

    mother_pt.Sumw2()
    daughter_pt.Sumw2()
    dl.Sumw2()
    radius.Sumw2()
    mass_rec.Sumw2()


mass.GetYaxis().SetTitle("Counts")
mass.GetXaxis().SetTitle("#it{M}(^{3}He + #pi^{-} and c.c.)   (GeV/#it{c}^{2})")

outfile = ROOT.TFile("histos.root", "recreate")

th2Dir = outfile.mkdir("th2")
th2Dir.cd()
res_vs_pgen.Write()
res_vs_pgen_dau.Write()
res_vs_radius.Write()

genDir = outfile.mkdir("gen")
genDir.cd()
mc_mother_pt.Write()
mc_daughter_pt.Write()
mc_dl.Write()
mc_radius.Write()

recDir = outfile.mkdir("rec")
recDir.cd()
mother_pt.Write()
daughter_pt.Write()
dl.Write()
radius.Write()
mass.Write()

chi2Dir = outfile.mkdir("chi2")
chi2Dir.cd()
chi2.Write()
chi2DCA.Write()
chi2_bkg.Write()
chi2DCA_bkg.Write()

resDir = outfile.mkdir("resolutions")
resDir.cd()
res_dau.Write()
radius_res.Write()
res.Write()


normalisedDir = outfile.mkdir("normalised")
normalisedDir.cd()

cv1 = ROOT.TCanvas()
cv1.SetName("chi2DCA")
cv1.SetLogy()
chi2DCA.SetLineColor(ROOT.kRed)
chi2DCA.SetTitle("Chi2DCA: signal in red, background in blue")
chi2DCA.SetMarkerStyle(2)
chi2DCA.DrawNormalized("EP")
chi2DCA_bkg.SetLineColor(ROOT.kBlue)
chi2DCA_bkg.SetMarkerStyle(2)
chi2DCA_bkg.DrawNormalized("EPsame")
cv1.Write()

cv2 = ROOT.TCanvas()
cv2.SetName("chi2Matching")
cv2.SetLogy()
chi2.SetLineColor(ROOT.kRed)
chi2.SetTitle("Chi2Matching: signal in red, background in blue")
chi2.SetMarkerStyle(2)
chi2.DrawNormalized("EP")
chi2_bkg.SetLineColor(ROOT.kBlue)
chi2_bkg.SetMarkerStyle(2)
chi2_bkg.DrawNormalized("EPsame")
cv2.Write()

cv3 = ROOT.TCanvas()
cv3.SetName("radius")
radius.SetLineColor(ROOT.kRed)
radius.SetTitle("Radius: signal in red, background in blue")
radius.SetMarkerStyle(2)
radius.DrawNormalized("EP")
radius_bkg.SetLineColor(ROOT.kBlue)
radius_bkg.SetMarkerStyle(2)
radius_bkg.DrawNormalized("EPsame")
cv3.Write()

cv4 = ROOT.TCanvas()
cv4.SetName("mass")
mass.SetLineColor(ROOT.kRed)
mass.SetTitle("Mass: signal in red, background in blue")
mass.SetMarkerStyle(2)
mass.DrawNormalized("EP")
mass_bkg.SetLineColor(ROOT.kBlue)
mass_bkg.SetMarkerStyle(2)
mass_bkg.DrawNormalized("EPsame")
cv4.Write()


## efficiency plots

effDir = outfile.mkdir("efficiencies")
effDir.cd()

eff_top = gen_pt.Clone("eff_pt")
eff_top.Divide(mc_mother_pt)
eff_top.SetStats(0)
eff_top.GetYaxis().SetTitle("Topology efficiency")
eff_top.GetXaxis().SetTitle("#it{p}_{T}   (GeV/#it{c})")
eff_top.SetTitle("")
eff_top.Write()


eff_radius_top = gen_r.Clone("eff_radius_top")
eff_radius_top.Divide(mc_radius)
eff_radius_top.SetStats(0)
eff_radius_top.GetYaxis().SetTitle("Topology efficiency")
eff_radius_top.GetXaxis().SetTitle("#it{R}   (cm)")
eff_radius_top.SetTitle("")
eff_radius_top.Write()
