import uproot
import pandas as pd
import numpy as np
import ROOT

def fill_th1_hist(df, var,h):
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th2_hist(df, var1, var2, h):
    for ind in range(len(df)):
        h.Fill(df.iloc[ind][var1], df.iloc[ind][var2])

usingKFP = False


kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC  = ROOT.TColor.GetColor("#ff7f00")

mc_mother_pt = ROOT.TH1F("mc_mother_pt", "mc_mother_pt", 100, 0, 12)
mc_daughter_pt = ROOT.TH1F("mc_daughter_pt", "mc_daughter_pt", 100, 0, 12)
mc_dl = ROOT.TH1F("mc_dl", "mc_dl", 100, 0, 50)
mc_radius = ROOT.TH1F("mc_radius", "mc_radius", 100, 0, 50)

mother_pt = ROOT.TH1F("mother_pt", "mother_pt", 100, 0, 12)
daughter_pt = ROOT.TH1F("daughter_pt", "daughter_pt", 100, 0, 12)
dl = ROOT.TH1F("dl", "dl", 100, 15, 60)
radius = ROOT.TH1F("radius", "radius", 100, 15, 60)
mass = ROOT.TH1F("mass", "mass", 100, 2.9, 3.2)
chi2 = ROOT.TH1F("chi2Matching", "chi2Matching", 100, 0, 100000)
chi2DCA = ROOT.TH1F("chi2DCA", "chi2DCA", 100, 0, 0.1)

chi2_bkg = ROOT.TH1F("chi2Matching_bkg", "chi2Matching_bkg", 100, 0, 100000)
chi2DCA_bkg = ROOT.TH1F("chi2DCA_bkg", "chi2DCA_bkg", 100, 0, 0.1)

radius_bkg = ROOT.TH1F("radius_bkg", "radius_bkg", 100, 15, 60)
mass_bkg = ROOT.TH1F("mass_bkg", "mass_bkg", 100, 2.9, 3.2)

gen_pt = ROOT.TH1F("gen_pt", "gen_pt", 100, 0, 12)
gen_r = ROOT.TH1F("gen_r", "gen_r", 100, 0, 50)

angle = ROOT.TH1F("angle", "angle", 100, 0, 5)
angle_bkg = ROOT.TH1F("angle_bkg", "angle_bkg", 100, 0, 5)

angle_vs_pgen = ROOT.TH2F("angle_vs_pgen", "angle_vs_pgen", 100, 0, 12, 100, 0, 5)
angle_vs_pgen_bkg = ROOT.TH2F("angle_vs_pgen_bkg", "angle_vs_pgen_bkg", 100, 0, 12, 100, 0, 5)

detector = ROOT.TH1F("Daughter detector", "Daughter detector; Detector[ITS(1),ITS-TPC(2),TPC(3),TPC-TOF(4),TPC-TRD(5),ITS-TPC-TRD(6),ITS-TPC-TOF(7),TPC-TRD-TOF(8),ITS-TPC-TRD-TOF(9)]; counts", 9, 0.5, 9.5)
detector_bkg = ROOT.TH1F("Daughter detector bkg", "Daughter detector bkg; Detector[ITS(1),ITS-TPC(2),TPC(3),TPC-TOF(4),TPC-TRD(5),ITS-TPC-TRD(6),ITS-TPC-TOF(7),TPC-TRD-TOF(8),ITS-TPC-TRD-TOF(9)]; counts", 9, 0.5, 9.5)

res_vs_pgen = ROOT.TH2F("res_vs_pgen", "res_vs_pgen", 100, 0, 12, 100, -3, 3)
res_vs_pgen_dau = ROOT.TH2F("res_vs_pgen_dau", "res_vs_pgen_dau", 100, 0, 12, 100, -3, 3)
res_vs_radius = ROOT.TH2F("res_vs_radius", "res_vs_radius", 100, 0, 50, 100, -3, 3)

res_vs_pgen_4L = ROOT.TH2F("res_vs_pgen_4L", "res_vs_pgen_4L", 100, 0, 12, 100, -1, 1)
res_vs_pgen_5L = ROOT.TH2F("res_vs_pgen_5L", "res_vs_pgen_5L", 100, 0, 12, 100, -1, 1)
res_vs_pgen_6L = ROOT.TH2F("res_vs_pgen_6L", "res_vs_pgen_6L", 100, 0, 12, 100, -1, 1)
res_vs_pgen_7L = ROOT.TH2F("res_vs_pgen_7L", "res_vs_pgen_7L", 100, 0, 12, 100, -1, 1)

if(usingKFP):
    name = "TrackedKFPKinkTree.root"
else:
    name = "TrackedKinkTree.root"
df_mc = uproot.open(name)["MCTree"].arrays(library="pd")

df_rec = uproot.open(name)["KinkTree"].arrays(library="pd").query("NLayers > 3").query("motherFake == 0").query("daughterFake == 0").query("Detector != 3")
df_topology = df_rec.query("isTopology == 1")
df_bkg = df_rec.query("isTopology == 0").query("isHyp == 0").query("isTriton == 0")

print(df_topology)

if(usingKFP != df_rec["usingKFP"].all()):
    print("ERROR: KFP flag is not the same in the tree")

mass_4Layers = ROOT.TH1F("mass_4Layers", "mass_4Layers", 50, 2.9, 3.2)
mass_5Layers = ROOT.TH1F("mass_5Layers", "mass_5Layers", 50, 2.9, 3.2)
mass_6Layers = ROOT.TH1F("mass_6Layers", "mass_6Layers", 50, 2.9, 3.2)
mass_7Layers = ROOT.TH1F("mass_7Layers", "mass_7Layers", 50, 2.9, 3.2)

resolution_4Layers = ROOT.TH1F("resolution_4Layers", "resolution_4Layers", 100, -1, 1)
resolution_5Layers = ROOT.TH1F("resolution_5Layers", "resolution_5Layers", 100, -1, 1)
resolution_6Layers = ROOT.TH1F("resolution_6Layers", "resolution_6Layers", 100, -1, 1)
resolution_7Layers = ROOT.TH1F("resolution_7Layers", "resolution_7Layers", 100, -1, 1)


res = ROOT.TH1F("res_mother", "res_mother", 100, -3, 3)
res_dau = ROOT.TH1F("res_dau", "res_dau", 100, -3, 3)
radius_res = ROOT.TH1F("radius_res", "radius_res", 100, -3, 3)

df_topology["res"] = (df_topology["rMotherPt"] - df_topology["gMotherPt"])/df_topology["gMotherPt"]
df_topology["res_dau"] = (df_topology["rDaughterPt"] - df_topology["gDaughterPt"])/df_topology["gDaughterPt"]
df_topology["radius_res"] = (df_topology["rRadius"] - df_topology["gRadius"])/df_topology["gRadius"]

df_4L = df_topology.query("NLayers == 4")
df_5L = df_topology.query("NLayers == 5")
df_6L = df_topology.query("NLayers == 6")
df_7L = df_topology.query("NLayers == 7")


fill_th1_hist(df_4L, "Mass", mass_4Layers)
fill_th1_hist(df_5L, "Mass", mass_5Layers)
fill_th1_hist(df_6L, "Mass", mass_6Layers)
fill_th1_hist(df_7L, "Mass", mass_7Layers)

fill_th1_hist(df_4L, "res", resolution_4Layers)
fill_th1_hist(df_5L, "res", resolution_5Layers)
fill_th1_hist(df_6L, "res", resolution_6Layers)
fill_th1_hist(df_7L, "res", resolution_7Layers)

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
fill_th1_hist(df_topology, "Detector", detector)

fill_th1_hist(df_bkg, "Chi2Match", chi2_bkg)
fill_th1_hist(df_bkg, "Chi2DCA", chi2DCA_bkg)
fill_th1_hist(df_bkg, "rRadius", radius_bkg)
fill_th1_hist(df_bkg, "Mass", mass_bkg)
fill_th1_hist(df_bkg, "Detector", detector_bkg)

fill_th1_hist(df_topology, "rAngle", angle)
fill_th1_hist(df_bkg, "rAngle", angle_bkg)

fill_th2_hist(df_topology, "gMotherPt", "rAngle", angle_vs_pgen)
fill_th2_hist(df_bkg, "gMotherPt", "rAngle", angle_vs_pgen_bkg)

fill_th2_hist(df_topology, "gMotherPt", "res", res_vs_pgen)
fill_th2_hist(df_topology, "gDaughterPt", "res_dau", res_vs_pgen_dau)
fill_th2_hist(df_topology, "gRadius", "radius_res", res_vs_radius)

fill_th2_hist(df_4L, "gMotherPt", "res", res_vs_pgen_4L)
fill_th2_hist(df_5L, "gMotherPt", "res", res_vs_pgen_5L)
fill_th2_hist(df_6L, "gMotherPt", "res", res_vs_pgen_6L)
fill_th2_hist(df_7L, "gMotherPt", "res", res_vs_pgen_7L)

mean_4L_array = np.array([])
mean_5L_array = np.array([])
mean_6L_array = np.array([])
mean_7L_array = np.array([])

std_dev_4L_array = np.array([])
std_dev_5L_array = np.array([])
std_dev_6L_array = np.array([])
std_dev_7L_array = np.array([])

std_dev_4L_error = np.array([])
std_dev_5L_error = np.array([])
std_dev_6L_error = np.array([])
std_dev_7L_error = np.array([])

pt_range = np.array([])
bins_step = 10
midstep = int(bins_step/2)

if (usingKFP):
    fitfile = ROOT.TFile("histos_KFP_fit.root", "recreate")
else:
    fitfile = ROOT.TFile("histos_fit.root", "recreate")

fitRange = 0.8
minEntries = 0

for idx in range(0, res_vs_pgen.GetNbinsX(), bins_step):
    pt = res_vs_pgen.GetXaxis().GetBinCenter(idx + midstep)
    pt_range = np.append(pt_range, pt)

    proj4 = res_vs_pgen_4L.ProjectionY("res_4L", idx, idx + bins_step -1)
    if(proj4.GetEntries() > minEntries):
        canv4 = ROOT.TCanvas()
        func4 = ROOT.TF1("func4", "gaus", -fitRange, fitRange)
        proj4.Fit(func4)
        std_dev_4L_error = np.append(std_dev_4L_error, func4.GetParError(2))
        std_dev_4L_array = np.append(std_dev_4L_array, func4.GetParameter(2))
        mean_4L_array = np.append(mean_4L_array, abs(res_vs_pgen_4L.ProjectionY("res_4L", idx, idx + bins_step -1).GetMean()))
        func4.Draw("same")
        canv4.Write()
    else:
        std_dev_4L_array = np.append(std_dev_4L_array, proj4.GetStdDev())
        std_dev_4L_error = np.append(std_dev_4L_error, 0)
        mean_4L_array = np.append(mean_4L_array, abs(res_vs_pgen_4L.ProjectionY("res_4L", idx, idx + bins_step -1).GetMean()))
        
    proj5 = res_vs_pgen_5L.ProjectionY("res_5L", idx, idx + bins_step -1)
    if(proj5.GetEntries() > minEntries):
        canv5 = ROOT.TCanvas()
        func5 = ROOT.TF1("func5", "gaus", -fitRange, fitRange)
        proj5.Fit(func5)
        std_dev_5L_error = np.append(std_dev_5L_error, func5.GetParError(2))
        std_dev_5L_array = np.append(std_dev_5L_array, func5.GetParameter(2))
        mean_5L_array = np.append(mean_5L_array, abs(res_vs_pgen_5L.ProjectionY("res_5L", idx, idx + bins_step -1).GetMean()))
        func5.Draw("same")
        canv5.Write()
    else:
        std_dev_5L_array = np.append(std_dev_5L_array, proj5.GetStdDev())
        std_dev_5L_error = np.append(std_dev_5L_error, 0)
        mean_5L_array = np.append(mean_5L_array, abs(res_vs_pgen_5L.ProjectionY("res_5L", idx, idx + bins_step -1).GetMean()))
    
    proj6 = res_vs_pgen_6L.ProjectionY("res_6L", idx, idx + bins_step -1)
    if(proj6.GetEntries() > minEntries):
        canv6 = ROOT.TCanvas()
        func6 = ROOT.TF1("func6", "gaus", -fitRange, fitRange)
        proj6.Fit(func6)
        std_dev_6L_error = np.append(std_dev_6L_error, func6.GetParError(2))
        std_dev_6L_array = np.append(std_dev_6L_array, func6.GetParameter(2))
        mean_6L_array = np.append(mean_6L_array, abs(res_vs_pgen_6L.ProjectionY("res_6L", idx, idx + bins_step -1).GetMean()))
        func6.Draw("same")
        canv6.Write()
    else:
        std_dev_6L_array = np.append(std_dev_6L_array, proj6.GetStdDev())
        std_dev_6L_error = np.append(std_dev_6L_error, 0)
        mean_6L_array = np.append(mean_6L_array, abs(res_vs_pgen_6L.ProjectionY("res_6L", idx, idx + bins_step -1).GetMean()))

    proj7 = res_vs_pgen_7L.ProjectionY("res_7L", idx, idx + bins_step -1)
    if(proj7.GetEntries() > minEntries):
        canv7 = ROOT.TCanvas()
        func7 = ROOT.TF1("func7", "gaus", -fitRange, fitRange)
        proj7.Fit(func7)
        std_dev_7L_error = np.append(std_dev_7L_error, func7.GetParError(2))
        std_dev_7L_array = np.append(std_dev_7L_array, func7.GetParameter(2))
        mean_7L_array = np.append(mean_7L_array, abs(res_vs_pgen_7L.ProjectionY("res_7L", idx, idx + bins_step -1).GetMean()))
        func7.Draw("same")
        canv7.Write()
    else:
        std_dev_7L_array = np.append(std_dev_7L_array, proj7.GetStdDev())
        std_dev_7L_error = np.append(std_dev_7L_error, 0)
        mean_7L_array = np.append(mean_7L_array, abs(res_vs_pgen_7L.ProjectionY("res_7L", idx, idx + bins_step -1).GetMean()))

ptStep = (pt_range[1] - pt_range[0])/2
ptStepArray = np.array([ptStep]*len(pt_range))

histo_4L = ROOT.TH1F("4 Layers mean", "4 Layers mean", len(pt_range), pt_range[0], pt_range[-1])
histo_5L = ROOT.TH1F("5 Layers mean", "5 Layers mean", len(pt_range), pt_range[0], pt_range[-1])
histo_6L = ROOT.TH1F("6 Layers mean", "6 Layers mean", len(pt_range), pt_range[0], pt_range[-1])
histo_7L = ROOT.TH1F("7 Layers mean", "7 Layers mean", len(pt_range), pt_range[0], pt_range[-1])

graph4 = ROOT.TGraphErrors(len(pt_range), pt_range, mean_4L_array, ptStepArray, std_dev_4L_error)
graph5 = ROOT.TGraphErrors(len(pt_range), pt_range, mean_5L_array, ptStepArray, std_dev_5L_error)
graph6 = ROOT.TGraphErrors(len(pt_range), pt_range, mean_6L_array, ptStepArray, std_dev_6L_error)
graph7 = ROOT.TGraphErrors(len(pt_range), pt_range, mean_7L_array, ptStepArray, std_dev_7L_error)

graph4.SetMarkerStyle(4)
graph4.SetMarkerColor(ROOT.kRed)
graph4.SetLineColor(ROOT.kRed)
graph4.SetDrawOption("EP")
graph4.SetTitle("4 Layers")
graph5.SetMarkerStyle(4)
graph5.SetMarkerColor(ROOT.kBlue)
graph5.SetLineColor(ROOT.kBlue)
graph5.SetDrawOption("EP")
graph5.SetTitle("5 Layers")
graph6.SetMarkerStyle(4)
graph6.SetMarkerColor(ROOT.kGreen)
graph6.SetLineColor(ROOT.kGreen)
graph6.SetDrawOption("EP")
graph6.SetTitle("6 Layers")
graph7.SetMarkerStyle(4)
graph7.SetMarkerColor(ROOT.kMagenta +4)
graph7.SetLineColor(ROOT.kMagenta +4)
graph7.SetDrawOption("EP")
graph7.SetTitle("7 Layers")


multigraph = ROOT.TMultiGraph("Standard Deviation", "Standard Deviation vs p_{T}; p_{T} (GeV/c); #sigma_{#it{M}} (GeV/c)")
multigraph.Add(graph4)
multigraph.Add(graph5)
multigraph.Add(graph6)
multigraph.Add(graph7)

multigraph.GetYaxis().SetRangeUser(0, 0.1)

for idx in range(0, len(pt_range)):
    histo_4L.SetBinContent(idx, mean_4L_array[idx])
    histo_5L.SetBinContent(idx, mean_5L_array[idx])
    histo_6L.SetBinContent(idx, mean_6L_array[idx])
    histo_7L.SetBinContent(idx, mean_7L_array[idx])


mass.GetYaxis().SetTitle("Counts")
mass.GetXaxis().SetTitle("#it{M}(^{3}He + #pi^{-} and c.c.)   (GeV/#it{c}^{2})")

if (usingKFP):
    outfile = ROOT.TFile("histosKFP.root", "recreate")
else: 
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

cv0 = ROOT.TCanvas()
cv0.SetName("Mass vs NLayers")
mass_6Layers.SetStats(0)
mass_6Layers.SetMarkerColor(ROOT.kGreen)
mass_6Layers.SetLineColor(ROOT.kGreen)
mass_6Layers.SetTitle("Mass vs NLayers")
mass_6Layers.DrawNormalized("line")
mass_4Layers.SetStats(0)
mass_4Layers.SetMarkerColor(ROOT.kRed)
mass_4Layers.SetLineColor(ROOT.kRed)
mass_4Layers.DrawNormalized("sameline")
mass_5Layers.SetMarkerColor(ROOT.kBlue)
mass_5Layers.SetLineColor(ROOT.kBlue)
mass_5Layers.DrawNormalized("sameline")
mass_7Layers.SetMarkerColor(ROOT.kMagenta +4)
mass_7Layers.SetLineColor(ROOT.kMagenta +4)
mass_7Layers.DrawNormalized("sameline")

legend0 = ROOT.TLegend(0.8,0.7,0.9,0.9)
legend0.AddEntry(mass_4Layers,"4 Layers","l")
legend0.AddEntry(mass_5Layers,"5 Layers","l")
legend0.AddEntry(mass_6Layers,"6 Layers","l")
legend0.AddEntry(mass_7Layers,"7 Layers","l")
legend0.Draw()
cv0.Write()

cv = ROOT.TCanvas()
cv.SetName("Pt Resolution vs NLayers")
resolution_7Layers.SetTitle("Pt Resolution vs NLayers")
resolution_7Layers.SetStats(0)
resolution_7Layers.SetMarkerColor(ROOT.kMagenta +4)
resolution_7Layers.SetLineColor(ROOT.kMagenta +4)
resolution_7Layers.DrawNormalized("line")
resolution_6Layers.SetStats(0)
resolution_6Layers.SetTitle("Pt Resolution vs NLayers")
resolution_6Layers.SetMarkerColor(ROOT.kGreen)
resolution_6Layers.SetLineColor(ROOT.kGreen)
resolution_6Layers.DrawNormalized("sameline")
resolution_5Layers.SetMarkerColor(ROOT.kBlue)
resolution_5Layers.SetLineColor(ROOT.kBlue)
resolution_5Layers.DrawNormalized("sameline")
resolution_4Layers.SetMarkerColor(ROOT.kRed)
resolution_4Layers.SetLineColor(ROOT.kRed)
resolution_4Layers.DrawNormalized("sameline")

legend = ROOT.TLegend(0.8,0.67,0.9,0.87)
legend.AddEntry(resolution_4Layers,"4 Layers","l")
legend.AddEntry(resolution_5Layers,"5 Layers","l")
legend.AddEntry(resolution_6Layers,"6 Layers","l")
legend.AddEntry(resolution_7Layers,"7 Layers","l")
legend.Draw()
cv.Write()

res_vs_pgen_4L.Write()
res_vs_pgen_5L.Write()
res_vs_pgen_6L.Write()
res_vs_pgen_7L.Write()

cvv = ROOT.TCanvas()
cvv.SetName("Mean Pt Resolution vs Generated Pt")
histo_4L.SetStats(0)
histo_4L.SetMarkerColor(ROOT.kRed)
histo_4L.SetMarkerStyle(4)
histo_4L.SetLineColor(ROOT.kRed)
histo_4L.SetTitle("Mean Pt Resolution vs Generated Pt")
histo_4L.GetYaxis().SetRangeUser(0,0.1)
histo_4L.Draw("P")
histo_5L.SetMarkerColor(ROOT.kBlue) 
histo_5L.SetMarkerStyle(4)
histo_5L.SetLineColor(ROOT.kBlue)
histo_5L.Draw("Psame")
histo_6L.SetMarkerColor(ROOT.kGreen)
histo_6L.SetMarkerStyle(4)
histo_6L.SetLineColor(ROOT.kGreen)
histo_6L.Draw("Psame")
histo_7L.SetMarkerColor(ROOT.kMagenta +4)
histo_7L.SetMarkerStyle(4)
histo_7L.SetLineColor(ROOT.kMagenta +4)
histo_7L.Draw("Psame")

legend2 = ROOT.TLegend(0.8,0.7,0.9,0.9)
legend2.AddEntry(histo_4L,"4 Layers","p")
legend2.AddEntry(histo_5L,"5 Layers","p")
legend2.AddEntry(histo_6L,"6 Layers","p")
legend2.AddEntry(histo_7L,"7 Layers","p")
legend2.Draw()
cvv.Write()


cvMG = ROOT.TCanvas("Standard Deviation")
multigraph.Draw()
cvMG.BuildLegend(0.8,0.67,0.9,0.87)
cvMG.Write()


detectorDir = outfile.mkdir("detector")
detectorDir.cd()
detector.Write()
detector_bkg.Write()


normalisedDir = outfile.mkdir("normalised")
normalisedDir.cd()

cv1 = ROOT.TCanvas()
cv1.SetName("chi2DCA")
cv1.SetLogy()
chi2DCA.SetStats(0)
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
chi2.SetStats(0)
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
radius.SetStats(0)
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
mass.SetStats(0)
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

angleDir = outfile.mkdir("angles")
angleDir.cd()
angle.Write()
angle_bkg.Write()
angle_vs_pgen.Write()
angle_vs_pgen_bkg.Write()

fitDir = outfile.mkdir("fit")
fitDir.cd()
canv = ROOT.TCanvas()
funcGaus = ROOT.TF1("funcGaus","gaus",2.9,3.1)
funcGaus.SetParLimits(1,2.9,3)
funcGaus.SetParLimits(2,0.01,0.1)
mass.Fit("funcGaus","R")
mass.Draw()
funcGaus.Draw("same")
canv.Write()
