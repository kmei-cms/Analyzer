import sys, os, ROOT, argparse

ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.3f")
ROOT.gStyle.SetFrameLineWidth(2)

parser = argparse.ArgumentParser()
parser.add_argument("--noVeto", dest="noVeto", help="Path to no veto input dir", default="NULL", type=str) 
parser.add_argument("--partialVeto", dest="partialVeto", help="Path to partialVeto input dir", default="NULL", type=str) 
parser.add_argument("--fullVeto", dest="fullVeto", help="Path to fullVeto input dir", default="NULL", type=str) 
arg = parser.parse_args()

OPTIONSMAP = {"Passbaseline"           : {"X" : {"title" : "N_{J}"}, "Y" : {"title" : "MVA bin"}, "Z" : {"min" : 0.001, "max" : 4e5}},
              "Passbaseline_Ratio"     : {"X" : {"title" : "N_{J}"}, "Y" : {"title" : "MVA bin"}, "Z" : {"min" : 0.0,   "max" : 0.75}},
              "Passbaseline_Raw"       : {"X" : {"title" : "N_{J}"}, "Y" : {"title" : "MVA bin"}, "Z" : {"min" : 0.001, "max" : 4e5}},
              "Passbaseline_Raw_Ratio" : {"X" : {"title" : "N_{J}"}, "Y" : {"title" : "MVA bin"}, "Z" : {"min" : 0.0,   "max" : 0.75}}}

def doOptions(histo, histoName):

    is1D = "TH1" in histo.ClassName()

    for axis, options in OPTIONSMAP[histoName].iteritems():

        if axis == "X":
            histo.GetXaxis().SetNdivisions(10)
            histo.GetXaxis().SetTickLength(0)
            if "rebin" in options:
                if is1D: histo.Rebin(options["rebin"])
                else: histo.RebinX(options["rebin"])
            if "min" in options and "max" in options: histo.GetXaxis().SetRangeUser(options["min"],options["max"])
            if "title" in options: histo.GetXaxis().SetTitle(options["title"])
        if axis == "Y":
            histo.GetYaxis().SetNdivisions(4)
            histo.GetYaxis().SetTickLength(0)
            if "rebin" in options:
                if is1D: histo.Rebin(options["rebin"])
                else: histo.RebinY(options["rebin"])
            if "min" in options and "max" in options: histo.SetMinimum(options["min"]); histo.SetMaximum(options["max"])
            if "title" in options: histo.GetYaxis().SetTitle(options["title"])
        if axis == "Z":
            if "min" in options and "max" in options: histo.SetMinimum(options["min"]); histo.SetMaximum(options["max"])

def prettyHisto(histo,magicFactor=1.0,magicFactor2=1.0):
    histo.GetYaxis().SetLabelSize(0.055); histo.GetYaxis().SetTitleSize(0.06); histo.GetYaxis().SetTitleOffset(0.5)
    histo.GetXaxis().SetLabelSize(0.055); histo.GetXaxis().SetTitleSize(0.06); histo.GetXaxis().SetTitleOffset(1.0)
    histo.GetZaxis().SetLabelSize(0.055); histo.GetZaxis().SetTitleSize(0.06)

def fillMap(inRootFile, tag):

    if ".root" not in inRootFile: return
    histoFile = ROOT.TFile.Open(inRootFile, "READ")
    for hkey in histoFile.GetListOfKeys():

        name = hkey.GetName()
        spec = name.split("_")[-1].split(".")[0]
        
        if "TH" not in hkey.GetClassName(): continue
        if "Passbaseline" not in name: continue

        histo = hkey.ReadObj()
        histo.SetDirectory(0)
        ROOT.SetOwnership(histo, False)

        histo.Sumw2()
        
        if name in MAPPFAHISTOS[tag].keys(): MAPPFAHISTOS[tag][name].Add(histo)
        else: MAPPFAHISTOS[tag][name] = histo

if __name__ == '__main__':

    XCANVAS = 2400; YCANVAS = 1400

    noVeto = arg.noVeto; partialVeto = arg.partialVeto; fullVeto = arg.fullVeto

    if arg.noVeto == "NULL" or arg.partialVeto == "NULL" or arg.fullVeto == "NULL": quit()

    stub = noVeto.split("_")[-1].split(".")[0]

    outpath = "./plots/HEMVetoStudy/"
    if not os.path.exists(outpath): os.makedirs(outpath)

    MAPPFAHISTOS = {"noVeto" : {}, "partialVeto" : {}, "fullVeto" : {}}

    fillMap(arg.noVeto,      "noVeto")
    fillMap(arg.partialVeto, "partialVeto")
    fillMap(arg.fullVeto,    "fullVeto")

    names = ["Passbaseline", "Passbaseline_Raw"]

    # Save the final histograms
    for name in names:
        magicMargins = {"T" : 0.03, "L" : 0.07, "B" : 0.13, "R" : 0.12}

        c1 = ROOT.TCanvas("%s_noVeto"%(name), "%s_noVeto"%(name), XCANVAS, YCANVAS); 
        ROOT.gPad.SetLogz()

        ROOT.gPad.SetTopMargin(magicMargins["T"])
        ROOT.gPad.SetLeftMargin(magicMargins["L"])
        ROOT.gPad.SetBottomMargin(magicMargins["B"])
        ROOT.gPad.SetRightMargin(magicMargins["R"])

        c2 = ROOT.TCanvas("%s_partialVeto"%(name), "%s_partialVeto"%(name), XCANVAS, YCANVAS); 
        ROOT.gPad.SetLogz()

        ROOT.gPad.SetTopMargin(magicMargins["T"])
        ROOT.gPad.SetLeftMargin(magicMargins["L"])
        ROOT.gPad.SetBottomMargin(magicMargins["B"])
        ROOT.gPad.SetRightMargin(magicMargins["R"])

        c3 = ROOT.TCanvas("%s_fullVeto"%(name), "%s_fullVeto"%(name), XCANVAS, YCANVAS); 
        ROOT.gPad.SetLogz()

        ROOT.gPad.SetTopMargin(magicMargins["T"])
        ROOT.gPad.SetLeftMargin(magicMargins["L"])
        ROOT.gPad.SetBottomMargin(magicMargins["B"])
        ROOT.gPad.SetRightMargin(magicMargins["R"])

        c4 = ROOT.TCanvas("%s_fullVetoRatio"%(name), "%s_fullVetoRatio"%(name), XCANVAS, YCANVAS); 

        ROOT.gPad.SetTopMargin(magicMargins["T"])
        ROOT.gPad.SetLeftMargin(magicMargins["L"])
        ROOT.gPad.SetBottomMargin(magicMargins["B"])
        ROOT.gPad.SetRightMargin(magicMargins["R"])

        c5 = ROOT.TCanvas("%s_partialVetoRatio"%(name), "%s_partialVetoRatio"%(name), XCANVAS, YCANVAS); 

        ROOT.gPad.SetTopMargin(magicMargins["T"])
        ROOT.gPad.SetLeftMargin(magicMargins["L"])
        ROOT.gPad.SetBottomMargin(magicMargins["B"])
        ROOT.gPad.SetRightMargin(magicMargins["R"])

        noVeto = MAPPFAHISTOS["noVeto"][name]; partialVeto = MAPPFAHISTOS["partialVeto"][name]; fullVeto = MAPPFAHISTOS["fullVeto"][name]
        fullVetoRatio    = noVeto.Clone("fullVetoRatio");    fullVetoRatio.Add(fullVeto, -1.0);       fullVetoRatio.Divide(noVeto)
        partialVetoRatio = noVeto.Clone("partialVetoRatio"); partialVetoRatio.Add(partialVeto, -1.0); partialVetoRatio.Divide(noVeto)
        noVeto.SetContour(255); partialVeto.SetContour(255); fullVeto.SetContour(255); fullVetoRatio.SetContour(255); partialVetoRatio.SetContour(255)

        prettyHisto(noVeto); prettyHisto(partialVeto); prettyHisto(fullVeto)
        prettyHisto(fullVetoRatio); prettyHisto(partialVetoRatio)

        theName      = noVeto.GetName()
        theNameRatio = theName + "_Ratio"

        if theName in OPTIONSMAP: doOptions(noVeto, theName); doOptions(partialVeto, theName); doOptions(fullVeto, theName); doOptions(fullVetoRatio, theNameRatio); doOptions(partialVetoRatio, theNameRatio)

        noVeto.SetTitle(""); partialVeto.SetTitle(""); fullVeto.SetTitle(""); fullVetoRatio.SetTitle(""); partialVetoRatio.SetTitle("")

        c1.cd(); noVeto.Draw("COLZ TEXT E");           c1.SaveAs("%s/%s_%s_noVeto.pdf"%(outpath,stub,name))
        c2.cd(); partialVeto.Draw("COLZ TEXT E");      c2.SaveAs("%s/%s_%s_partialVeto.pdf"%(outpath,stub,name))
        c3.cd(); fullVeto.Draw("COLZ TEXT E");         c3.SaveAs("%s/%s_%s_fullVeto.pdf"%(outpath,stub,name))
        c4.cd(); fullVetoRatio.Draw("COLZ TEXT E");    c4.SaveAs("%s/%s_%s_fullVetoRatio.pdf"%(outpath,stub,name))
        c5.cd(); partialVetoRatio.Draw("COLZ TEXT E"); c5.SaveAs("%s/%s_%s_partialVetoRatio.pdf"%(outpath,stub,name))
