import sys, os, ROOT, argparse

ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetErrorX(0)

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--inputDir", dest="inputDir", help="Path to input directory"    , default="NULL", type=str) 
parser.add_argument("--year"    , dest="year"    , help="What year are we looking at", default="NULL", type=str) 

arg = parser.parse_args()

OPTIONSMAP = {"h_njets_nvtx_1l_HT300_ge3j_ge1b_Mbl"              : {"X" : {"rebin" : 2,  "min" : 0, "max" : 100, "title" : "NVtx"}, "Y" : {"title" : "N_{J}"}},
              "h_njets_nvtx_1l_HT300_ge3j_ge1b_Mbl_noPuWeight"   : {"X" : {"rebin" : 2,  "min" : 0, "max" : 100, "title" : "NVtx"}, "Y" : {"title" : "N_{J}"}},
}

def doOptions(histo, histoName):

    is1D = "TH1" in histo.ClassName()

    for axis, options in OPTIONSMAP[histoName].iteritems():

        if axis == "X":
            if "rebin" in options:
                if is1D: histo.Rebin(options["rebin"])
                else: histo.RebinX(options["rebin"])
            if "min" in options and "max" in options: histo.GetXaxis().SetRangeUser(options["min"],options["max"])
            if "title" in options: histo.GetXaxis().SetTitle(options["title"])
        if axis == "Y":
            if "rebin" in options:
                if is1D: histo.Rebin(options["rebin"])
                else: histo.RebinY(options["rebin"])
            if "min" in options and "max" in options: histo.GetYaxis().SetRangeUser(options["min"],options["max"])
            if "title" in options: histo.GetYaxis().SetTitle(options["title"])
        if axis == "Z":
            if "min" in options and "max" in options: histo.GetZaxis().SetRangeUser(options["min"],options["max"])

def prettyHisto(histo,magicFactor=1.0,magicFactor2=1.0):
    histo.GetYaxis().SetLabelSize(0.045); histo.GetYaxis().SetTitleSize(0.06); histo.GetYaxis().SetTitleOffset(0.8)
    histo.GetXaxis().SetLabelSize(0.045); histo.GetXaxis().SetTitleSize(0.06); histo.GetXaxis().SetTitleOffset(0.8)
    histo.GetZaxis().SetLabelSize(0.045); histo.GetZaxis().SetTitleSize(0.06)

def fillMap(inRootFile, theMap):

    keyName = inRootFile.split("_")[-1].split(".root")[0]

    if ".root" not in inRootFile: return
    histoFile = ROOT.TFile.Open(inRootFile, "READ")
    for hkey in histoFile.GetListOfKeys():
        if "TH" not in hkey.GetClassName(): continue

        name = hkey.GetName()
        if name not in OPTIONSMAP: continue

        histo = hkey.ReadObj()
        histo.SetDirectory(0)

        histo.Sumw2()
        
        if name in theMap[keyName].keys(): theMap[keyName][name].Add(histo)
        else: theMap[keyName][name] = histo

if __name__ == '__main__':

    XCANVAS = 2400; YCANVAS = 2400

    inputDir = arg.inputDir
    year     = arg.year

    stub = ""
    if inputDir == "NULL": quit()
    else: stub = inputDir.split("/")[-1]

    outpath = "./plots/AppendixC/%s/"%(stub)
    if not os.path.exists(outpath): os.makedirs(outpath)

    mapPFAhistos = {"TT" : {}, "Data" : {}}

    ttbarFile = inputDir + "/%s_TT.root"%(year)
    dataFile  = inputDir + "/%s_Data.root"%(year)

    fillMap(ttbarFile, mapPFAhistos)
    fillMap(dataFile,  mapPFAhistos)

    # Save the final histograms
    magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}

    c1 = ROOT.TCanvas("Data_2D", "Data_2D", XCANVAS, YCANVAS); c1.cd(); c1.SetLogz()

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])

    nominalName = "h_njets_nvtx_1l_HT300_ge3j_ge1b_Mbl"; noPuWeightName = "h_njets_nvtx_1l_HT300_ge3j_ge1b_Mbl_noPuWeight" 
    ttNominal    = mapPFAhistos["TT"][nominalName];    prettyHisto(ttNominal,0.875,1.1);    doOptions(ttNominal, nominalName);       ttNominal.SetTitle("NJets v Number of Vertices");    ttNominal.SetContour(255)
    ttNoPuWeight = mapPFAhistos["TT"][noPuWeightName]; prettyHisto(ttNoPuWeight,0.875,1.1); doOptions(ttNoPuWeight, noPuWeightName); ttNoPuWeight.SetTitle("NJets v Number of Vertices"); ttNoPuWeight.SetContour(255)
    data         = mapPFAhistos["Data"][nominalName];  prettyHisto(data,0.875,1.1);         doOptions(data, nominalName);            data.SetTitle("NJets v Number of Vertices");         data.SetContour(255)

    ttNominalProf    = ttNominal.ProfileX("ttNominal", 1, -1, "");       ttNominalProf.SetLineColor(ROOT.kBlue);   ttNominalProf.SetMarkerColor(ROOT.kBlue);    ttNominalProf.SetLineWidth(3);    ttNominalProf.SetTitle("Comparison of Profile with and without PU Weights")
    ttNoPuWeightProf = ttNoPuWeight.ProfileX("ttNoPuWeight", 1, -1, ""); ttNoPuWeightProf.SetLineColor(ROOT.kRed); ttNoPuWeightProf.SetMarkerColor(ROOT.kRed); ttNoPuWeightProf.SetLineWidth(3); ttNoPuWeightProf.SetTitle("Comparison of Profile with and without PU Weights")
    dataProf         = data.ProfileX("data", 1, -1, "");                 dataProf.SetLineColor(ROOT.kBlue);        data.SetMarkerColor(ROOT.kCyan);             dataProf.SetLineWidth(3);         dataProf.SetTitle("Comparison of Profile with and without PU Weights")

    data.Draw("COLZ")
    dataProf.Draw("HIST SAME")

    c1.SaveAs("%s/Data_NJets_vs_NVtx.pdf"%(outpath))

    c3 = ROOT.TCanvas("TT_2D_nom", "TT_2D_nom", XCANVAS, YCANVAS); c3.cd(); c3.SetLogz()

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])

    ttNominal.Draw("COLZ")
    ttNominalProf.Draw("HIST SAME")

    c3.SaveAs("%s/TT_Nominal_NJets_vs_NVtx.pdf"%(outpath))

    c4 = ROOT.TCanvas("TT_2D_noPU", "TT_2D_noPU", XCANVAS, YCANVAS); c4.cd(); c4.SetLogz()

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])

    ttNoPuWeight.Draw("COLZ")
    ttNoPuWeightProf.Draw("HIST SAME")

    c4.SaveAs("%s/TT_NoPuWeight_NJets_vs_NVtx.pdf"%(outpath))

    c2 = ROOT.TCanvas("_1D", "_1D", XCANVAS, YCANVAS); c2.cd()

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])

    dataProf.SetLineColor(ROOT.kCyan)

    ttNominalProf.SetMinimum(3)
    ttNominalProf.SetMaximum(9)

    iamLegend = ROOT.TLegend(0.12, 0.82, 0.55, 0.92)
    iamLegend.AddEntry(ttNominalProf, "with PU weight")
    iamLegend.AddEntry(ttNoPuWeightProf, "without PU weight")
    iamLegend.AddEntry(dataProf, "Data")

    ttNominalProf.Draw("HIST"); ttNoPuWeightProf.Draw("HIST SAME"); dataProf.Draw("HIST SAME"); iamLegend.Draw("SAME")

    c2.SaveAs("%s/Profiles_Comparison.pdf"%(outpath))
