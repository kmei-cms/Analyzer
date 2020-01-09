import sys, os, ROOT, argparse
from collections import defaultdict

ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--inputDir", dest="inputDir", help="Path to input", default="NULL", type=str) 

arg = parser.parse_args()

OPTIONSMAP = {"h_njets_1l_HT300_ge7j_ge1b_Mbl_d1" : {"X" : {"min" : 7, "max" : 15, "title" : "N_{J} D1"}},
              "h_njets_1l_HT300_ge7j_ge1b_Mbl_d2" : {"X" : {"min" : 7, "max" : 15, "title" : "N_{J} D2"}},
              "h_njets_1l_HT300_ge7j_ge1b_Mbl_d3" : {"X" : {"min" : 7, "max" : 15, "title" : "N_{J} D3"}},
              "h_njets_1l_HT300_ge7j_ge1b_Mbl_d4" : {"X" : {"min" : 7, "max" : 15, "title" : "N_{J} D4"}},
              "h_njets_1l_HT300_ge7j_ge1b_Mbl"    : {"X" : {"min" : 7, "max" : 15, "title" : "N_{J}"}}
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
    histo.GetYaxis().SetLabelSize(magicFactor*0.055); histo.GetYaxis().SetTitleSize(magicFactor*0.08); histo.GetYaxis().SetTitleOffset(0.7/magicFactor)
    histo.GetXaxis().SetLabelSize(magicFactor*0.055); histo.GetXaxis().SetTitleSize(magicFactor*0.08); histo.GetXaxis().SetTitleOffset(0.8/magicFactor2)
    histo.GetZaxis().SetLabelSize(magicFactor*0.055); histo.GetZaxis().SetTitleSize(magicFactor*0.06)

def fillMap(inRootFile, theMap):

    if ".root" not in inRootFile: return
    histoFile = ROOT.TFile.Open(inRootFile, "READ")
    for hkey in histoFile.GetListOfKeys():
        if "TH" not in hkey.GetClassName(): continue

        if hkey.GetName() == "EventCounter" or hkey.GetName().find("njets") == -1: continue

        histo = hkey.ReadObj()
        histo.SetDirectory(0)

        histo.Sumw2()
        
        theMap.setdefault(hkey.GetName(), histo)

if __name__ == '__main__':

    XCANVAS = 2400; YCANVAS = 2400

    if arg.inputDir == "NULL": quit()
    stub = arg.inputDir.split("condor/")[-1]

    inRootFile = arg.inputDir + "/2017_MC.root"
       
    outpath = "./plots/%s/"%(stub)
    if not os.path.exists(outpath): os.makedirs(outpath)

    mapPFAhistos = {}

    fillMap(inRootFile, mapPFAhistos)

    # Save the final histograms

    njetsD1 = mapPFAhistos["h_njets_1l_HT300_ge7j_ge1b_Mbl_d1"]; prettyHisto(njetsD1)
    njetsD2 = mapPFAhistos["h_njets_1l_HT300_ge7j_ge1b_Mbl_d2"]; prettyHisto(njetsD2)
    njetsD3 = mapPFAhistos["h_njets_1l_HT300_ge7j_ge1b_Mbl_d3"]; prettyHisto(njetsD3)
    njetsD4 = mapPFAhistos["h_njets_1l_HT300_ge7j_ge1b_Mbl_d4"]; prettyHisto(njetsD4)

    njets   = mapPFAhistos["h_njets_1l_HT300_ge7j_ge1b_Mbl"]; prettyHisto(njets)

    XMin = 0; XMax = 1
    YMin = 0; YMax = 1

    njetsD1.SetTitle(""); njetsD1.Scale(1./njetsD1.Integral()); doOptions(njetsD1, "h_njets_1l_HT300_ge7j_ge1b_Mbl_d1")
    njetsD2.SetTitle(""); njetsD2.Scale(1./njetsD2.Integral()); doOptions(njetsD2, "h_njets_1l_HT300_ge7j_ge1b_Mbl_d2")
    njetsD3.SetTitle(""); njetsD3.Scale(1./njetsD3.Integral()); doOptions(njetsD3, "h_njets_1l_HT300_ge7j_ge1b_Mbl_d3")
    njetsD4.SetTitle(""); njetsD4.Scale(1./njetsD4.Integral()); doOptions(njetsD4, "h_njets_1l_HT300_ge7j_ge1b_Mbl_d4")
    njets.SetTitle("");   njets.Scale(1./njets.Integral());     doOptions(njets,   "h_njets_1l_HT300_ge7j_ge1b_Mbl")

    njetsD1.Divide(njets); njetsD1.SetMinimum(0.50); njetsD1.SetMaximum(1.50); njetsD1.GetYaxis().SetNdivisions(308)
    njetsD2.Divide(njets); njetsD2.SetMinimum(0.50); njetsD2.SetMaximum(1.50); njetsD2.GetYaxis().SetNdivisions(308)
    njetsD3.Divide(njets); njetsD3.SetMinimum(0.50); njetsD3.SetMaximum(1.50); njetsD3.GetYaxis().SetNdivisions(308)
    njetsD4.Divide(njets); njetsD4.SetMinimum(0.50); njetsD4.SetMaximum(1.50); njetsD4.GetYaxis().SetNdivisions(308)

    njetsD1.SetMarkerColor(ROOT.kBlack);   njetsD1.SetLineColor(ROOT.kBlack);   njetsD1.SetMarkerSize(4); njetsD1.SetMarkerStyle(20); njetsD1.SetLineWidth(3)
    njetsD2.SetMarkerColor(ROOT.kRed);     njetsD2.SetLineColor(ROOT.kRed);     njetsD2.SetMarkerSize(4); njetsD2.SetMarkerStyle(20); njetsD2.SetLineWidth(3)
    njetsD3.SetMarkerColor(ROOT.kBlue);    njetsD3.SetLineColor(ROOT.kBlue);    njetsD3.SetMarkerSize(4); njetsD3.SetMarkerStyle(20); njetsD3.SetLineWidth(3)
    njetsD4.SetMarkerColor(ROOT.kGreen+2); njetsD4.SetLineColor(ROOT.kGreen+2); njetsD4.SetMarkerSize(4); njetsD4.SetMarkerStyle(20); njetsD4.SetLineWidth(3)

    mvaBins = ["D1", "D2", "D3", "D4"]

    for mva in mvaBins:

        c1 = ROOT.TCanvas("njets%s"%(mva), "njets%s"%(mva), XCANVAS, YCANVAS); 
        c1.cd(); ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)

        ROOT.gPad.SetGridy(); ROOT.gPad.SetGridx()
        ROOT.gPad.SetTopMargin(0.03)
        ROOT.gPad.SetLeftMargin(0.11)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetRightMargin(0.04)

        if   mva == "D1": njetsD1.Draw("L")
        elif mva == "D2": njetsD2.Draw("L")
        elif mva == "D3": njetsD3.Draw("L")
        elif mva == "D4": njetsD4.Draw("L")

        c1.SaveAs("%s/njets%s_Total_Ratio.pdf"%(outpath,mva))
