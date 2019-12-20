import sys, os, ROOT, argparse

ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)

parser = argparse.ArgumentParser()
parser.add_argument("--tt2016", dest="tt2016", help="Path to 2016 fit input dir", default="NULL", type=str) 
parser.add_argument("--tt2017", dest="tt2017", help="Path to 2017 fit input dir", default="NULL", type=str) 
parser.add_argument("--tt2018pre", dest="tt2018pre", help="Path to 2018pre fit input dir", default="NULL", type=str) 
parser.add_argument("--tt2018post", dest="tt2018post", help="Path to 2018post fit input dir", default="NULL", type=str) 
arg = parser.parse_args()

OPTIONSMAP = {"D1_TT_h_njets_pt30_1l" : {"X" : {"title" : "N_{J} (Shifted)"}, "Y" : {"min" : 2e-4, "max" : 2}},
              "D2_TT_h_njets_pt30_1l" : {"X" : {"title" : "N_{J} (Shifted)"}, "Y" : {"min" : 2e-4, "max" : 2}},
              "D3_TT_h_njets_pt30_1l" : {"X" : {"title" : "N_{J} (Shifted)"}, "Y" : {"min" : 2e-4, "max" : 2}},
              "D4_TT_h_njets_pt30_1l" : {"X" : {"title" : "N_{J} (Shifted)"}, "Y" : {"min" : 2e-4, "max" : 2}}
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
            if "min" in options and "max" in options: histo.SetMinimum(options["min"]); histo.SetMaximum(options["max"])
            if "title" in options: histo.GetYaxis().SetTitle(options["title"])
        if axis == "Z":
            if "min" in options and "max" in options: histo.GetZaxis().SetRangeUser(options["min"],options["max"])

def prettyHisto(histo,magicFactor=1.0,magicFactor2=1.0):
    histo.GetYaxis().SetLabelSize(0.055); histo.GetYaxis().SetTitleSize(0.06); histo.GetYaxis().SetTitleOffset(0.7)
    histo.GetXaxis().SetLabelSize(0.055); histo.GetXaxis().SetTitleSize(0.06); histo.GetXaxis().SetTitleOffset(1.0)
    histo.GetZaxis().SetLabelSize(0.055); histo.GetZaxis().SetTitleSize(0.06)

def fillMap(inRootFile, year):

    if ".root" not in inRootFile: return
    histoFile = ROOT.TFile.Open(inRootFile, "READ")
    for hkey in histoFile.GetListOfKeys():

        name = hkey.GetName()
        spec = name.split("_")[-1].split(".")[0]
        
        if "TH" not in hkey.GetClassName(): continue
        if name == "EventCounter": continue
        if "_TT_" not in name or spec != "1l": continue

        histo = hkey.ReadObj()
        histo.SetDirectory(0)
        ROOT.SetOwnership(histo, False)

        histo.Sumw2()
        
        if name in MAPPFAHISTOS[year].keys(): MAPPFAHISTOS[year][name].Add(histo)
        else: MAPPFAHISTOS[year][name] = histo

if __name__ == '__main__':

    XCANVAS = 2400; YCANVAS = 2400

    tt2016 = arg.tt2016; tt2017 = arg.tt2017; tt2018pre = arg.tt2018pre; tt2018post = arg.tt2018post

    if arg.tt2016 == "NULL" or arg.tt2017 == "NULL" or arg.tt2018pre == "NULL" or arg.tt2018post == "NULL": quit()

    outpath = "./plots/NJetsTTbar/"
    if not os.path.exists(outpath): os.makedirs(outpath)

    MAPPFAHISTOS = {"2016" : {}, "2017" : {}, "2018pre" : {}, "2018post" : {}}

    fillMap(arg.tt2016,     "2016")
    fillMap(arg.tt2017,     "2017")
    fillMap(arg.tt2018pre,  "2018pre")
    fillMap(arg.tt2018post, "2018post")

    # Save the final histograms
    for name in MAPPFAHISTOS.values()[0].keys():

        c1 = ROOT.TCanvas("%s"%(name), "%s"%(name), XCANVAS, YCANVAS); 
        c1.cd(); ROOT.gPad.SetLogy(); ROOT.gPad.SetLogz()

        ROOT.gPad.SetGridy(); ROOT.gPad.SetGridx()
        ROOT.gPad.SetTopMargin(0.03)
        ROOT.gPad.SetLeftMargin(0.11)
        ROOT.gPad.SetBottomMargin(0.13)
        ROOT.gPad.SetRightMargin(0.04)

        tt2016 = MAPPFAHISTOS["2016"][name]; tt2017 = MAPPFAHISTOS["2017"][name]; tt2018pre = MAPPFAHISTOS["2018pre"][name]; tt2018post = MAPPFAHISTOS["2018post"][name]

        prettyHisto(tt2016); prettyHisto(tt2017); prettyHisto(tt2018pre); prettyHisto(tt2018post)

        theName = tt2016.GetName().replace("_2016", "")

        tt2016.SetMarkerColor(ROOT.kBlack);       tt2016.SetLineColor(ROOT.kBlack);       tt2016.SetMarkerSize(4);     tt2016.SetLineWidth(4);     tt2016.SetMarkerStyle(20);     tt2016.Scale(1./tt2016.Integral())
        tt2017.SetMarkerColor(ROOT.kBlue);        tt2017.SetLineColor(ROOT.kBlue);        tt2017.SetMarkerSize(4);     tt2017.SetLineWidth(4);     tt2017.SetMarkerStyle(20);     tt2017.Scale(1./tt2017.Integral())
        tt2018pre.SetMarkerColor(ROOT.kRed);      tt2018pre.SetLineColor(ROOT.kRed);      tt2018pre.SetMarkerSize(4);  tt2018pre.SetLineWidth(4);  tt2018pre.SetMarkerStyle(20);  tt2018pre.Scale(1./tt2018pre.Integral())
        tt2018post.SetMarkerColor(ROOT.kGreen+2); tt2018post.SetLineColor(ROOT.kGreen+2); tt2018post.SetMarkerSize(4); tt2018post.SetLineWidth(4); tt2018post.SetMarkerStyle(20); tt2018post.Scale(1./tt2018post.Integral())

        if theName in OPTIONSMAP: doOptions(tt2016, theName); doOptions(tt2017, theName); doOptions(tt2018pre, theName); doOptions(tt2018post, theName)

        tt2016.SetTitle(""); tt2017.SetTitle(""); tt2018pre.SetTitle(""); tt2018post.SetTitle("")

        iamLegend = ROOT.TLegend(0.70, 0.6, 0.9, 0.95) 
        iamLegend.SetLineWidth(0)
        iamLegend.AddEntry(tt2016, "2016")
        iamLegend.AddEntry(tt2017, "2017")
        iamLegend.AddEntry(tt2018pre, "2018pre")
        iamLegend.AddEntry(tt2018post, "2018post")

        tt2016.Draw("EP"); #tt2016.Draw("HIST SAME")
        tt2017.Draw("EP SAME"); #tt2017.Draw("HIST SAME")
        tt2018pre.Draw("EP SAME"); #tt2018pre.Draw("HIST SAME")
        tt2018post.Draw("EP SAME"); #tt2018post.Draw("HIST SAME")
        iamLegend.Draw("SAME")

        c1.SaveAs("%s/%s.pdf"%(outpath,name))
