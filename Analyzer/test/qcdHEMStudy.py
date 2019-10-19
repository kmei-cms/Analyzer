import sys, os, ROOT, argparse
from collections import defaultdict

ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetErrorX(0)

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--tag"     , dest="tag"     , help="Unique tag for output"    , type=str, default="")
parser.add_argument("--preHEM", dest="preHEM", help="Path to QCD pre HEM", default="NULL", type=str) 
parser.add_argument("--postHEM", dest="postHEM", help="Path to QCD post HEM", default="NULL", type=str) 
parser.add_argument("--ratio", dest="ratio", help="Draw ratio", action="store_true", default=False) 

arg = parser.parse_args()

optionsMap = {"h_njets"        : {"X" : {"rebin" : 1, "min" : 4.5, "max" : 15.5, "title" : "N_{J}"}},
              "h_ht"           : {"X" : {"rebin" : 20, "min" : 0, "max" : 5000, "title" : "H_{T} [GeV]"}},
              "h_deepESM"      : {"X" : {"rebin" : 5, "min" : 0, "max" : 1, "title" : "DeepESM"}},
}

def doOptions(histo, histoName, theMap):

    is1D = "TH1" in histo.ClassName()

    for axis, options in theMap[histoName].iteritems():

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
    histo.GetYaxis().SetLabelSize(magicFactor*0.055); histo.GetYaxis().SetTitleSize(magicFactor*0.08); histo.GetYaxis().SetTitleOffset(0.6/magicFactor)
    histo.GetXaxis().SetLabelSize(magicFactor*0.055); histo.GetXaxis().SetTitleSize(magicFactor*0.08); histo.GetXaxis().SetTitleOffset(0.7/magicFactor2)
    histo.GetZaxis().SetLabelSize(magicFactor*0.055); histo.GetZaxis().SetTitleSize(magicFactor*0.06)

def fillMap(inRootDir, theMap):

    for histoFile in os.listdir(inRootDir):

        if ".root" not in histoFile: continue
        histoFile = ROOT.TFile.Open(inRootDir + "/" + histoFile, "READ")

        for hkey in histoFile.GetListOfKeys():

            if "TH" not in hkey.GetClassName(): continue
            if hkey.GetName() == "EventCounter": continue
            if "passQCDCR" not in hkey.GetName(): continue

            vetoTag = ""
            if "passHEMVeto" in hkey.GetName(): vetoTag = "veto"
            else: vetoTag = "noveto"

            name = hkey.GetName()
            name = name.replace("_passHEMVeto","").replace("_passQCDCR", "")
            histo = hkey.ReadObj()
            histo.SetDirectory(0)

            histo.Sumw2()
            
            theMap.setdefault(name, {}).setdefault(vetoTag, histo)

if __name__ == '__main__':

    ratio = arg.ratio
    XCANVAS = 2400; YCANVAS = 2400

    if arg.preHEM != "NULL" and arg.postHEM != "NULL": stub = arg.preHEM.split("condor/")[-1]
    else: quit()

    tag = arg.tag

    preHEMDir = arg.preHEM
    postHEMDir = arg.postHEM
       
    outpath = "./plots/HEMStudy/%s/%s"%(stub,tag)
    if not os.path.exists(outpath): os.makedirs(outpath)

    mapPreHEMhistos = {}
    mapPostHEMhistos = {}

    fillMap(preHEMDir, mapPreHEMhistos)
    fillMap(postHEMDir, mapPostHEMhistos)

    # Save the final histograms
    for name in mapPreHEMhistos:
        for vetoTag in ["veto", "noveto"]:

            magicMargins = {"T" : 0.02625, "B" : 0.13375, "L" : 0.11, "R" : 0.12}

            preHEM  = mapPreHEMhistos[name][vetoTag].Clone()
            postHEM = mapPostHEMhistos[name][vetoTag].Clone()

            XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1 
            YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
            PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)

            c1 = ROOT.TCanvas("%s_%s"%(name,vetoTag), "%s_%s"%(name,vetoTag), XCANVAS, YCANVAS); 
            c1.Divide(1,2); c1.cd(1); ROOT.gPad.SetLogy(); ROOT.gPad.SetLogz(); ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)

            ROOT.gPad.SetGridy(); ROOT.gPad.SetGridx()
            ROOT.gPad.SetTopMargin(0.03)
            ROOT.gPad.SetLeftMargin(0.11)
            ROOT.gPad.SetBottomMargin(0.02)
            ROOT.gPad.SetRightMargin(0.04)

            ratio = ROOT.TH1F("ratio_%s_%s"%(name,vetoTag), "ratio_%s_%s"%(name,vetoTag), preHEM.GetNbinsX(), preHEM.GetXaxis().GetXmin(), preHEM.GetXaxis().GetXmax()) 

            prettyHisto(preHEM); prettyHisto(postHEM); prettyHisto(ratio,PadFactor)

            theName = preHEM.GetName().replace("_passQCDCR", "").replace("_passHEMVeto", "")
            if theName in optionsMap:

                doOptions(preHEM, theName, optionsMap)
                doOptions(postHEM, theName, optionsMap)
                doOptions(ratio, theName, optionsMap)

            else: break

            postHEM.SetMarkerColor(ROOT.kRed); postHEM.SetLineColor(ROOT.kRed); postHEM.SetMarkerSize(3); postHEM.SetLineWidth(2); postHEM.SetMarkerStyle(20)
            preHEM.SetMarkerColor(ROOT.kBlack);    preHEM.SetLineColor(ROOT.kBlack);    preHEM.SetMarkerSize(3);  preHEM.SetLineWidth(2);  preHEM.SetMarkerStyle(20)
            postHEM.SetTitle(""); preHEM.SetTitle("")

            preHEM.Scale(1./preHEM.Integral()); postHEM.Scale(1./postHEM.Integral())

            #preHEM.GetYaxis().SetRangeUser(0.001,1.1*preHEM.GetMaximum())

            preHEM.GetXaxis().SetLabelSize(0)

            preHEM.Draw("EP")
            postHEM.Draw("EP SAME")

            c1.cd(2)

            ROOT.gPad.SetGridy()
            ROOT.gPad.SetTopMargin(0.10)
            ROOT.gPad.SetBottomMargin(0.30)
            ROOT.gPad.SetRightMargin(0.04)
            ROOT.gPad.SetLeftMargin(0.11)
            ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)

            ratio.SetTitle("")
            ratio.Divide(preHEM,postHEM)

            ratio.GetYaxis().SetRangeUser(0.50,1.5)
            ratio.GetYaxis().SetNdivisions(-304)
            ratio.GetYaxis().SetTitle("HEM / Nominal")
            ratio.GetYaxis().SetTitleSize(ratio.GetYaxis().GetTitleSize()*0.6)
            ratio.GetXaxis().SetTitleSize(ratio.GetXaxis().GetTitleSize()*0.8)
            ratio.GetYaxis().SetTitleOffset(0.5)
            ratio.GetXaxis().SetTitleOffset(0.9)
            ratio.SetMarkerStyle(20); ratio.SetMarkerSize(3); ratio.SetMarkerColor(ROOT.kBlue+2)
            ratio.SetLineWidth(2); ratio.SetLineColor(ROOT.kBlue+2)

            ratio.Draw("EP")
            c1.SaveAs("%s/%s_%s.pdf"%(outpath,name,vetoTag))
