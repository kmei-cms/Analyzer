import sys, os, ROOT, argparse
from collections import defaultdict

ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetErrorX(0)

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--inputDir", dest="inputDir", help="Path to input", default="NULL", type=str) 

arg = parser.parse_args()

OPTIONSMAP = {"h_njets" : {"X" : {"min" : 7, "max" : 15, "title" : "N_{J}"}}}

def makeBandGraph(histoUp, histoDown, color = ROOT.kRed):
    
    npoints = histoUp.GetNbinsX() / 2

    graphUp = ROOT.TGraph(npoints)
    graphDown = ROOT.TGraph(npoints)
    graphBand = ROOT.TGraph(2*npoints)

    offset = 7
    for iPoint in xrange(0, npoints+1):

        graphUp.SetPoint(iPoint, histoUp.GetBinCenter(iPoint+offset+1), histoUp.GetBinContent(iPoint+offset+1))
        graphDown.SetPoint(iPoint, histoUp.GetBinCenter(iPoint+offset+1), histoDown.GetBinContent(iPoint+offset+1))

        graphBand.SetPoint(iPoint, histoUp.GetBinCenter(iPoint+offset+1), histoUp.GetBinContent(iPoint+offset+1))
        graphBand.SetPoint(npoints+iPoint, histoUp.GetBinCenter(npoints-iPoint+offset+2) , histoDown.GetBinContent(npoints-iPoint+offset+2))

    graphUp.SetLineColor(color); graphUp.SetLineWidth(2)
    graphDown.SetLineColor(color); graphDown.SetLineWidth(2)
    graphBand.SetFillColorAlpha(color, 0.35)

    return graphUp, graphDown, graphBand

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

        if hkey.GetName() == "EventCounter": continue

        keyStubs = hkey.GetName().split("_")

        mva = keyStubs[-2]
        variation = keyStubs[-1]

        # Little corner case...
        if mva == "hemcorr":
            mva = variation
            variation = "hemcorr"

        histo = hkey.ReadObj()
        histo.SetDirectory(0)

        histo.Sumw2()
        
        theMap.setdefault(mva, {}).setdefault(variation, histo)

if __name__ == '__main__':

    XCANVAS = 2400; YCANVAS = 2400

    if arg.inputDir == "NULL": quit()
    stub = arg.inputDir.split("condor/")[-1].split(".root")[0].split("_")[-1]

    inRootFile = arg.inputDir
       
    outpath = "./plots/NJetsStudy/%s/"%(stub)
    if not os.path.exists(outpath): os.makedirs(outpath)

    mapPFAhistos = {}

    fillMap(inRootFile, mapPFAhistos)

    mvas = ["d1", "d2", "d3", "d4"]
    # Save the final histograms
    for mva in mvas:

        hemCorr = mapPFAhistos[mva]["hemcorr"]; prettyHisto(hemCorr)
        JECUp = mapPFAhistos[mva]["JECup"]; prettyHisto(JECUp)
        JECDown = mapPFAhistos[mva]["JECdown"]; prettyHisto(JECDown)
        JERUp = mapPFAhistos[mva]["JERup"]; prettyHisto(JERUp)
        JERDown = mapPFAhistos[mva]["JERdown"]; prettyHisto(JERDown)
        totUp = mapPFAhistos[mva]["up"]; prettyHisto(totUp)
        totDown = mapPFAhistos[mva]["down"]; prettyHisto(totDown)

        XMin = 0; XMax = 1
        YMin = 0; YMax = 1

        JECUp.SetTitle("");   JECUp.Scale(1./JECUp.Integral())
        JECDown.SetTitle(""); JECDown.Scale(1./JECDown.Integral())
        JERUp.SetTitle("");   JERUp.Scale(1./JERUp.Integral())
        JERDown.SetTitle(""); JERDown.Scale(1./JERDown.Integral())
        totUp.SetTitle("");   totUp.Scale(1./totUp.Integral())
        totDown.SetTitle(""); totDown.Scale(1./totDown.Integral())

        doOptions(hemCorr, "h_njets"); doOptions(JECUp, "h_njets"); doOptions(JECDown, "h_njets")

        jecup, jecdown, jecband = makeBandGraph(JECUp, JECDown)
        jerup, jerdown, jerband = makeBandGraph(JERUp, JERDown, ROOT.kBlue)
        totup, totdown, totband = makeBandGraph(totUp, totDown, ROOT.kGreen+2)

        hemCorr.SetTitle(""); hemCorr.SetMarkerColor(ROOT.kBlack); hemCorr.SetLineColor(ROOT.kBlack); hemCorr.SetMarkerSize(4); hemCorr.SetLineWidth(3); hemCorr.SetMarkerStyle(20); hemCorr.Scale(1./hemCorr.Integral())

        hemCorr.SetNdivisions(16)

        c1 = ROOT.TCanvas("%s"%(mva), "%s"%(mva), XCANVAS, YCANVAS); 
        c1.cd(); ROOT.gPad.SetLogy(); ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)

        ROOT.gPad.SetGridy(); ROOT.gPad.SetGridx()
        ROOT.gPad.SetTopMargin(0.03)
        ROOT.gPad.SetLeftMargin(0.11)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetRightMargin(0.04)

        hemCorr.Draw()
        jecband.Draw("F")
        jecup.Draw("L")
        jecdown.Draw("L")

        jerband.Draw("F")
        jerup.Draw("L")
        jerdown.Draw("L")

        c1.SaveAs("%s/%s_JECJER.pdf"%(outpath,mva))

        c2 = ROOT.TCanvas("%s_tot"%(mva), "%s_tot"%(mva), XCANVAS, YCANVAS); 
        c2.cd(); ROOT.gPad.SetLogy(); ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)

        ROOT.gPad.SetGridy(); ROOT.gPad.SetGridx()
        ROOT.gPad.SetTopMargin(0.03)
        ROOT.gPad.SetLeftMargin(0.11)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetRightMargin(0.04)

        hemCorr.Draw()
        totband.Draw("F")
        totup.Draw("L")
        totdown.Draw("L")
        hemCorr.Draw("EP SAME")

        c2.SaveAs("%s/%s_total.pdf"%(outpath,mva))
