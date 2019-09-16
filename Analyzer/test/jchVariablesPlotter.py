import os
import sys
import ROOT
import math
import argparse
import collections

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetLineWidth(3)
ROOT.gStyle.SetFrameLineWidth(3)

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--inputDir"  , dest="inputDir"  , help="the input directory", required=True)
parser.add_argument("--tag"       , dest="tag"       , help="output name tag", required=True)
parser.add_argument("--units"     , dest="units"     , help="x axis units", default=None)
parser.add_argument("--setXTitle" , dest="xTitle" , nargs="+", help="Set x-axis title", default=None, type=str)
parser.add_argument("--setYTitle" , dest="yTitle" , nargs="+", help="Set y-axis title", default=None, type=str)
parser.add_argument("--setLogX"   , dest="setLogX" , help="make x-axis log scale", action="store_true", default=False)
parser.add_argument("--setLogY"   , dest="setLogY" , help="make y-axis log scale", action="store_true", default=False)
parser.add_argument("--setLogZ"   , dest="setLogZ" , help="make z-axis log scale", action="store_true", default=False)
parser.add_argument("--setXRange" , dest="xRange" , nargs="+", help="set x range", default=None, type=float)
parser.add_argument("--setYRange" , dest="yRange" , nargs="+", help="set y range", default=None, type=float)
parser.add_argument("--rebinX"    , dest="rebinX" , help="rebin histo in X", default=None, type=int)
parser.add_argument("--histoName" , dest="histoName" , help="histo name to grab",required=True,type=str)
parser.add_argument("--legend"    , dest="legend" , help="draw a legend",action="store_true", default=False)
parser.add_argument("--Project3D" , dest="project3D", nargs="+", help="See ROOT doc.", default=None)
parser.add_argument("--ProjectionX" , dest="projectionX", help="See ROOT doc.", type=int, default=None)
parser.add_argument("--dimension" , dest="dimension", help="histo dim", default=1, type=int)
parser.add_argument("--normalize" , dest="normalize" , help="normalize the hist", action="store_true", default=False)
parser.add_argument("--same"      , dest="same"     , help="Draw same", action="store_true", default=False)
arg = parser.parse_args()

HISTONAME = arg.histoName

HISTOCONF = [["2016_TT.root", ROOT.TColor.GetColor("#7a0019"), "TTbar"],
             ["2016_RPV_2t6j_mStop-300.root", ROOT.kRed, "RPV 2t+6j, m_{#tilde{t}} = 300 GeV"],
             ["2016_RPV_2t6j_mStop-350.root", ROOT.kRed, "RPV 2t+6j, m_{#tilde{t}} = 350 GeV"],
             ["2016_RPV_2t6j_mStop-400.root", ROOT.kRed, "RPV 2t+6j, m_{#tilde{t}} = 400 GeV"],
             ["2016_RPV_2t6j_mStop-450.root", ROOT.kRed, "RPV 2t+6j, m_{#tilde{t}} = 450 GeV"],
             ["2016_RPV_2t6j_mStop-500.root", ROOT.kGray+2, "RPV 2t+6j, m_{#tilde{t}} = 500 GeV"],
             ["2016_RPV_2t6j_mStop-700.root", ROOT.kGreen+1, "RPV 2t+6j, m_{#tilde{t}} = 700 GeV"],
             ["2016_RPV_2t6j_mStop-900.root", ROOT.kBlue, "RPV 2t+6j, m_{#tilde{t}} = 900 GeV"]]

#HISTOCONF = [["2016_TT.root", ROOT.TColor.GetColor("#7a0019"), "TTbar"],
#                ["2016_StealthSYY_2t6j_mStop-300.root", ROOT.kRed, "SYY 2t+6j, m_{#tilde{t}} = 300 GeV"],
#                ["2016_StealthSYY_2t6j_mStop-500.root", ROOT.kGray+2, "SYY 2t+6j, m_{#tilde{t}} = 500 GeV"],
#                ["2016_StealthSYY_2t6j_mStop-700.root", ROOT.kGreen+1, "SYY 2t+6j, m_{#tilde{t}} = 700 GeV"],
#                ["2016_StealthSYY_2t6j_mStop-900.root", ROOT.kBlue, "SYY 2t+6j, m_{#tilde{t}} = 900 GeV"]]

def getVarHisto(datasetFile):

    f = ROOT.TFile.Open(datasetFile, "READ")
    h_varHisto = f.Get(HISTONAME).Clone()
    h_varHisto.SetDirectory(0)
    f.Close()

    h_varHisto.SetTitle("")
    return h_varHisto 

def printEffs(aHisto):

    total = aHisto.Integral()

    running = total 
    for abin in xrange(aHisto.GetNbinsX()):

        if abin == 0: continue

        running -= aHisto.GetBinContent(abin)

        print "Cut at %f has efficiency %f"%(aHisto.GetBinCenter(abin),running/total)

def prettyHisto(histoConf):

    if arg.project3D != None:

        histoConf[-1].GetZaxis().FindBin(arg.project3D[0])
        histoConf[-1].GetZaxis().FindBin(arg.project3D[0])

        h2 = 0
        if len(arg.project3D) == 0:
            h2 = histoConf[-1].Project3D("xy")
        elif len(arg.project3D) == 1:
            histoConf[-1].GetZaxis().SetRange(histoConf[-1].GetZaxis().FindBin(arg.project3D[0]),histoConf[-1].GetZaxis().FindBin(arg.project3D[0]))
            histoConf[-1].GetZaxis().SetBit(ROOT.TAxis.kAxisRange)
            h2 = histoConf[-1].Project3D("xy")
        elif len(arg.project3D) == 2:
            histoConf[-1].GetZaxis().SetRange(histoConf[-1].GetZaxis().FindBin(arg.project3D[0]),histoConf[-1].GetZaxis().FindBin(arg.project3D[1]))
            histoConf[-1].GetZaxis().SetBit(ROOT.TAxis.kAxisRange)
            h2 = histoConf[-1].Project3D("xy")    

        h2.SetName("")
        h2.SetTitle("")
        histoConf[-1] = h2

    if arg.projectionX != None:

        hProjX = histoConf[-1].ProjectionX(histoConf[-1].GetName()+"_proj_%s"%(histoConf[0]), histoConf[-1].GetYaxis().FindBin(arg.projectionX), histoConf[-1].GetYaxis().FindBin(arg.projectionX))

        histoConf[-1] = hProjX

    # Do all the things to the histoConf[-1]
    if arg.rebinX and arg.dimension == 1:
        histoConf[-1].Rebin(arg.rebinX)

    if arg.dimension == 1:
        histoConf[-1].SetLineColor(histoConf[1])
        histoConf[-1].SetLineWidth(0)
        histoConf[-1].SetMarkerSize(0)
        if not "TT" in histoConf[2]:
            histoConf[-1].SetLineStyle(7)
            histoConf[-1].SetLineWidth(5)
        if "TT" in histoConf[2]: histoConf[-1].SetFillColorAlpha(histoConf[1], 0.3)
    if arg.dimension == 2:
        histoConf[-1].SetContour(255)

    if arg.xRange:
        histoConf[-1].GetXaxis().SetRangeUser(arg.xRange[0], arg.xRange[1])
    if arg.yRange:
        if arg.dimension == 1:
            histoConf[-1].SetMinimum(arg.yRange[0])
            histoConf[-1].SetMaximum(arg.yRange[1])
        else:
            histoConf[-1].GetYaxis().SetRangeUser(arg.yRange[0], arg.yRange[1])
    if arg.units != None:
        if not arg.normalize:
            histoConf[-1].GetYaxis().SetTitle("Events / %.3g %s"%(histoConf[-1].GetBinWidth(1), arg.units))
        else:
            histoConf[-1].GetYaxis().SetTitle("Fraction of Events / %.3g %s"%(histoConf[-1].GetBinWidth(1), arg.units))
    elif arg.yTitle != None:
        histoConf[-1].GetYaxis().SetTitle(" ".join(arg.yTitle))

    histoConf[-1].GetXaxis().SetTitle(" ".join(arg.xTitle))
    histoConf[-1].GetXaxis().SetTitleSize(XTitleSize)
    histoConf[-1].GetXaxis().SetLabelSize(XLabelSize)
    histoConf[-1].GetYaxis().SetTitleSize(YTitleSize)
    histoConf[-1].GetYaxis().SetLabelSize(YLabelSize)
    histoConf[-1].GetYaxis().SetTitleOffset(YTitleOffset)
    histoConf[-1].GetXaxis().SetTitleOffset(XTitleOffset)

    if arg.normalize:
        if histoConf[-1].Integral() > 0: histoConf[-1].Scale(1./histoConf[-1].Integral())

if __name__ == '__main__':

    XTitleSize = 0.035;  XLabelSize = 0.03;  XTitleOffset = 1.0 
    YTitleSize = 0.035;  YLabelSize = 0.03;  YTitleOffset = 1.3 

    StackXMin = 0;    StackXMax = 1; RatioXMin = 0; RatioXMax = 1
    StackYMin = 0.20; StackYMax = 1; RatioYMin = 0; RatioYMax = 0.20

    PadFactor = (StackYMax-StackYMin) / (RatioYMax-RatioYMin)
    
    DataMarkSize = 1.3

    outputDir = "/uscms/home/jhiltb/nobackup/susy/CMSSW_9_3_3/src/Analyzer/Analyzer/test/plots"
    inputDir = arg.inputDir

    iamLegend = 0
    if arg.legend:
        iamLegend = ROOT.TLegend(0.68, 0.76, 0.97, 0.98, "", "trNDC")
        iamLegend.SetTextSize(0.023)

    for histoConf in HISTOCONF:
    
        process = histoConf[2]
        varHisto = getVarHisto(inputDir+"/"+histoConf[0])
        histoConf.append(varHisto)
        prettyHisto(histoConf)

    canvas = ROOT.TCanvas("c1", "c1", 1200, 1200)
    ROOT.gPad.Clear()

    if arg.setLogX: ROOT.gPad.SetLogx()        
    if arg.setLogY: ROOT.gPad.SetLogy()        
    if arg.setLogZ: ROOT.gPad.SetLogz()        

    ROOT.gPad.SetPad(StackXMin, StackYMin, StackXMax, StackYMax)
    ROOT.gPad.SetTopMargin(0.01)
    ROOT.gPad.SetBottomMargin(0.08)
    if arg.dimension == 2:
        ROOT.gPad.SetRightMargin(0.15)
    else:
        ROOT.gPad.SetRightMargin(0.025)
    ROOT.gPad.SetLeftMargin(0.09)
    ROOT.gPad.SetTicks()

    canvas.cd()
    if arg.same:
        for histoConf in HISTOCONF:

            if arg.legend:
                if histoConf[-1].Integral() > 0:
                    if "TT" not in histoConf[0]:
                        iamLegend.AddEntry(histoConf[-1], histoConf[2], "L")
                    else:
                        iamLegend.AddEntry(histoConf[-1], histoConf[2], "F")

            if arg.dimension == 2: histoConf[-1].Draw("COLZ SAME")
            else: histoConf[-1].Draw("HIST SAME")

        if arg.legend: iamLegend.Draw("SAME")
        canvas.SaveAs("%s/StackPlot_%s_%s.pdf"%(outputDir,HISTONAME,arg.tag))
    else:
        for histoConf in HISTOCONF:

            if arg.legend:
                if histoConf[-1].Integral() > 0:
                    if "TT" not in histoConf[0]:
                        iamLegend.AddEntry(histoConf[-1], histoConf[2], "L")
                    else:
                        iamLegend.AddEntry(histoConf[-1], histoConf[2], "F")

            if arg.dimension == 2: histoConf[-1].Draw("COLZ")
            else: histoConf[-1].Draw("HIST")

            if arg.legend: iamLegend.Draw("SAME")
            canvas.SaveAs("%s/StackPlot_%s_%s_%s.pdf"%(outputDir,HISTONAME,histoConf[0].split(".root")[0],arg.tag))
            canvas.Clear()
