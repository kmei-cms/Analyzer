import sys, os, ROOT, argparse

ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetErrorX(0)

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--data2018A", dest="data2018A", help="Path to 2018A data dir", default="NULL", type=str) 
parser.add_argument("--data2018D", dest="data2018D", help="Path to 2018D data dir", default="NULL", type=str) 
parser.add_argument("--ratio", dest="ratio", help="Draw ratio", action="store_true", default=False) 

arg = parser.parse_args()

OPTIONSMAP = {"h_jet_pt"       : {"X" : {"rebin" : 12, "min" : 0, "max" : 1000, "title" : "Jet p_{T} [GeV]"}},
              "h_jet_eta"      : {"X" : {"rebin" : 12, "title" : "#eta"}},
              "h_jet_etaphi"   : {"X" : {"rebin" : 30, "title" : "#eta"}, "Y" : {"rebin" : 30, "title" : "#phi"}},
              "h_jet_pteta"    : {"X" : {"rebin" : 3,  "min" : 0, "max" : 200, "title" : "p_{T} [GeV]"}, "Y" : {"rebin" : 30, "title" : "#eta"}},
              "h_jet_ptem"     : {"X" : {"rebin" : 3,  "min" : 0, "max" : 200, "title" : "p_{T} [GeV]"}, "Y" : {"rebin" : 30, "title" : "Neutral EM Fraction"}},
              "h_jet_pthad"    : {"X" : {"rebin" : 3,  "min" : 0, "max" : 200, "title" : "p_{T} [GeV]"}, "Y" : {"rebin" : 30, "title" : "Neutral Hadron Fraction"}},
              "h_jet_etaem"    : {"X" : {"rebin" : 30, "title" : "#eta"}, "Y" : {"rebin" : 30, "title" : "Neutral EM Fraction"}},
              "h_jet_etahad"   : {"X" : {"rebin" : 30, "title" : "#eta"}, "Y" : {"rebin" : 30, "title" : "Neutral Hadron Fraction"}},
              "h_jet30_pt"     : {"X" : {"rebin" : 12, "min" : 0, "max" : 1000, "title" : "Jet p_{T} [GeV]"}},
              "h_jet30_eta"    : {"X" : {"rebin" : 12, "title" : "#eta"}},
              "h_jet30_etaphi" : {"X" : {"rebin" : 30, "title" : "#eta"}, "Y" : {"rebin" : 30, "title" : "#phi"}},
              "h_jet30_pteta"  : {"X" : {"rebin" : 3,  "min" : 0, "max" : 200, "title" : "p_{T} [GeV]"}, "Y" : {"rebin" : 30, "title" : "#eta"}},
              "h_jet30_ptem"   : {"X" : {"rebin" : 3,  "min" : 0, "max" : 200, "title" : "p_{T} [GeV]"}, "Y" : {"rebin" : 30, "title" : "Neutral EM Fraction"}},
              "h_jet30_pthad"  : {"X" : {"rebin" : 3,  "min" : 0, "max" : 200, "title" : "p_{T} [GeV]"}, "Y" : {"rebin" : 30, "title" : "Neutral Hadron Fraction"}},
              "h_jet30_etaem"  : {"X" : {"rebin" : 30, "title" : "#eta"}, "Y" : {"rebin" : 30, "title" : "Neutral EM Fraction"}},
              "h_jet30_etahad" : {"X" : {"rebin" : 30, "title" : "#eta"}, "Y" : {"rebin" : 30, "title" : "Neutral Hadron Fraction"}},
              "h_fwm2_top6"    : {"X" : {"rebin" : 20, "title" : "Fox-Wolfram 2nd Moment"}},
              "h_fwm3_top6"    : {"X" : {"rebin" : 20, "title" : "Fox-Wolfram 3rd Moment"}},
              "h_fwm4_top6"    : {"X" : {"rebin" : 20, "title" : "Fox-Wolfram 4th Moment"}},
              "h_fwm5_top6"    : {"X" : {"rebin" : 20, "title" : "Fox-Wolfram 5th Moment"}},
              "h_jmt_ev0_top6" : {"X" : {"min" : 0.0, "max" : 1.0, "rebin" : 20, "title" : "Jet Momentum Tensor Eigenvalue 0"}},
              "h_jmt_ev1_top6" : {"X" : {"min" : 0.0, "max" : 1.0, "rebin" : 20, "title" : "Jet Momentum Tensor Eigenvalue 1"}},
              "h_jmt_ev2_top6" : {"X" : {"min" : 0.0, "max" : 1.0, "rebin" : 20, "title" : "Jet Momentum Tensor Eigenvalue 2"}},
              "h_beta_z"       : {"X" : {"rebin" : 20, "title" : "#beta_{z}"}},
              "h_beta_z_pt20"  : {"X" : {"rebin" : 20, "title" : "#beta_{z}"}},
              "h_beta_z_nvtx"       : {"X" : {"rebin" : 20, "title" : "#beta_{z}"}, "Y" : {"title" : "Number of Vertices"}},
              "h_beta_z_pt20_nvtx"  : {"X" : {"rebin" : 20, "title" : "#beta_{z}"}, "Y" : {"title" : "Number of Vertices"}},
              "h_deepESM"  : {"X" : {"rebin" : 20, "title" : "DeepESM"}},
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
    histo.GetXaxis().SetLabelSize(magicFactor*0.055); histo.GetXaxis().SetTitleSize(magicFactor*0.08); histo.GetXaxis().SetTitleOffset(1.0/magicFactor2)
    histo.GetZaxis().SetLabelSize(magicFactor*0.055); histo.GetZaxis().SetTitleSize(magicFactor*0.06)

def fillMap(inRootFile, theMap):

    if ".root" not in inRootFile: return
    histoFile = ROOT.TFile.Open(inRootFile, "READ")
    for hkey in histoFile.GetListOfKeys():
        if "TH" not in hkey.GetClassName(): continue

        if hkey.GetName() == "EventCounter": continue

        # Will be either 2018D or 2018A
        keyName = hkey.GetName().split("_")[-1]

        #if keyName == "2018D" and "2018pre" in inRootFile: continue
        #if keyName == "2018A" and "2018post" in inRootFile: continue

        name = hkey.GetName()
        name = name.replace("_2018A","").replace("_2018D", "")
        histo = hkey.ReadObj()
        histo.SetDirectory(0)

        histo.Sumw2()
        
        if name in theMap[keyName].keys(): theMap[keyName][name].Add(histo)
        else: theMap[keyName][name] = histo

if __name__ == '__main__':

    ratio = arg.ratio
    XCANVAS = 2400; YCANVAS = 2400

    if arg.data2018D == "NULL" or arg.data2018A == "NULL": quit()
    if arg.data2018A != "NULL": stub = arg.data2018A.split("condor/")[-1].split(".root")[0].split("_")[-1]

    inRootFile2018A = arg.data2018A
    inRootFile2018D = arg.data2018D
       
    outpath = "./plots/EENoiseStudy/%s/"%(stub)
    if not os.path.exists(outpath): os.makedirs(outpath)

    mapPFAhistos = {"2018A" : {}, "2018D" : {}}

    fillMap(inRootFile2018A, mapPFAhistos)
    fillMap(inRootFile2018D, mapPFAhistos)

    # Save the final histograms
    for name in mapPFAhistos.values()[0].keys():

        magicMargins = {"T" : 0.02625, "B" : 0.15, "L" : 0.11, "R" : 0.12}


        if "TH2" in mapPFAhistos.values()[0][name].ClassName():
            
            c1 = ROOT.TCanvas("%s_2018A"%(name), "%s_2018D"%(name), XCANVAS, int(0.8*YCANVAS)); c1.cd()

            ROOT.gPad.SetTopMargin(magicMargins["T"])
            ROOT.gPad.SetBottomMargin(magicMargins["B"])
            ROOT.gPad.SetLeftMargin(magicMargins["L"])
            ROOT.gPad.SetRightMargin(magicMargins["R"])

            data2018A = mapPFAhistos["2018A"][name]
            theName = data2018A.GetName().replace("_2018A", "")
            prettyHisto(data2018A,0.875,1.1)
            doOptions(data2018A, theName)

            data2018A.SetTitle(""); data2018A.SetContour(255); 

            data2018A.Draw("COLZ")

            c1.SaveAs("%s/%s_2018A.pdf"%(outpath,name))

            c2 = ROOT.TCanvas("%s_2018D"%(name), "%s_2018D"%(name), XCANVAS, int(0.8*YCANVAS)); c2.cd()

            ROOT.gPad.SetTopMargin(magicMargins["T"])
            ROOT.gPad.SetBottomMargin(magicMargins["B"])
            ROOT.gPad.SetLeftMargin(magicMargins["L"])
            ROOT.gPad.SetRightMargin(magicMargins["R"])

            data2018D = mapPFAhistos["2018D"][name]
            prettyHisto(data2018D,0.875,1.1)
            doOptions(data2018D, theName)

            data2018D.SetTitle(""); data2018D.SetContour(255)

            data2018D.Draw("COLZ")

            c2.SaveAs("%s/%s_2018D.pdf"%(outpath,name))

            if data2018A.Integral() != 0:

                c1 = ROOT.TCanvas("%s_ratio"%(name), "%s_ratio"%(name), XCANVAS, int(0.8*YCANVAS)); c1.cd()

                ROOT.gPad.SetTopMargin(magicMargins["T"])
                ROOT.gPad.SetBottomMargin(magicMargins["B"])
                ROOT.gPad.SetLeftMargin(magicMargins["L"])
                ROOT.gPad.SetRightMargin(magicMargins["R"])

                data2018A.Scale(1./data2018A.Integral())
                data2018D.Scale(1./data2018D.Integral())

                data2018D.Divide(data2018A)

                data2018D.GetZaxis().SetRangeUser(0.5,3.0)
                
                data2018D.Draw("COLZ")
                c1.SaveAs("%s/%s_ratio.pdf"%(outpath,name))

        elif "TH1" in mapPFAhistos.values()[0][name].ClassName():

            if ratio:
                XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1 
                YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
                PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)

                c1 = ROOT.TCanvas("%s"%(name), "%s"%(name), XCANVAS, YCANVAS); 
                c1.Divide(1,2); c1.cd(1); ROOT.gPad.SetLogy(); ROOT.gPad.SetLogz(); ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)

                ROOT.gPad.SetGridy(); ROOT.gPad.SetGridx()
                ROOT.gPad.SetTopMargin(0.03)
                ROOT.gPad.SetLeftMargin(0.11)
                ROOT.gPad.SetBottomMargin(0.02)
                ROOT.gPad.SetRightMargin(0.04)

                data2018A = mapPFAhistos["2018A"][name]; data2018D = mapPFAhistos["2018D"][name]
                ratio = ROOT.TH1F("ratio_%s"%(name), "ratio_%s"%(name), data2018A.GetNbinsX(), data2018A.GetXaxis().GetXmin(), data2018A.GetXaxis().GetXmax()) 

                prettyHisto(data2018A); prettyHisto(data2018D); prettyHisto(ratio,PadFactor)

                theName = data2018A.GetName().replace("_2018A", "")
                if theName in OPTIONSMAP:
                    doOptions(data2018A, theName); doOptions(data2018D, theName); doOptions(ratio, theName)

                data2018A.SetMarkerColor(ROOT.kBlack); data2018A.SetLineColor(ROOT.kBlack); data2018A.SetMarkerSize(3); data2018A.SetLineWidth(2); data2018A.SetMarkerStyle(20); data2018A.Scale(1./data2018A.Integral())
                data2018D.SetMarkerColor(ROOT.kRed);   data2018D.SetLineColor(ROOT.kRed);   data2018D.SetMarkerSize(3); data2018D.SetLineWidth(2); data2018D.SetMarkerStyle(20); data2018D.Scale(1./data2018D.Integral())
                data2018A.SetTitle(""); data2018D.SetTitle("")

                data2018A.GetXaxis().SetLabelSize(0)

                data2018A.Draw("EP"); data2018D.Draw("EP SAME")

                c1.cd(2)

                ROOT.gPad.SetGridy()
                ROOT.gPad.SetTopMargin(0.10)
                ROOT.gPad.SetBottomMargin(0.35)
                ROOT.gPad.SetRightMargin(0.04)
                ROOT.gPad.SetLeftMargin(0.11)
                ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)

                ratio.SetTitle("")
                ratio.Divide(data2018D,data2018A)

                ratio.GetYaxis().SetRangeUser(0.25,1.75)
                ratio.GetYaxis().SetNdivisions(604)
                ratio.GetYaxis().SetTitle("2018D / 2018A")
                ratio.GetYaxis().SetTitleSize(ratio.GetYaxis().GetTitleSize()*0.6)
                ratio.GetXaxis().SetTitleSize(ratio.GetXaxis().GetTitleSize()*0.8)
                ratio.GetYaxis().SetTitleOffset(0.5)
                ratio.GetXaxis().SetTitleOffset(1.08)
                ratio.SetMarkerStyle(20); ratio.SetMarkerSize(3); ratio.SetMarkerColor(ROOT.kBlue+2)
                ratio.SetLineWidth(2); ratio.SetLineColor(ROOT.kBlue+2)

                ratio.Draw("EP")
                c1.SaveAs("%s/%s.pdf"%(outpath,name))
