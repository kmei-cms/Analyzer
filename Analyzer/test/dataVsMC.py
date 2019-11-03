import ROOT
ROOT.gROOT.SetBatch(True)
import math
import DataSetInfo as info
import optparse

def getData(path, scale=1.0, year = "2016"):
    Data = [
        info.DataSetInfo(basedir=path, fileName=year+"_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]
    
    bgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_Triboson.root",        sys= 0.3, label="Triboson",        scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_Diboson.root",         sys= 0.3, label="Diboson",         scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_DYJetsToLL_M-50.root", sys= 0.3, label="DYJetsToLL_M-50", scale=scale),        
        info.DataSetInfo(basedir=path, fileName=year+"_TTX.root",             sys= 0.3, label="TTX",             scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",           sys= 0.3, label="WJets",           scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              sys= 0.3, label="ST",              scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",             sys= 0.5, label="QCD",             scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_TT.root",              sys= 0.2, label="T#bar{T}",        scale=scale),
    ]

    sgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-300.root", sys=-1.0, label="RPV 300", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-350.root", sys=-1.0, label="RPV 350", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-400.root", sys=-1.0, label="RPV 400", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-450.root", sys=-1.0, label="RPV 450", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-500.root", sys=-1.0, label="RPV 500", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-550.root", sys=-1.0, label="RPV 550", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-600.root", sys=-1.0, label="RPV 600", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-650.root", sys=-1.0, label="RPV 650", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-700.root", sys=-1.0, label="RPV 700", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-750.root", sys=-1.0, label="RPV 750", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-800.root", sys=-1.0, label="RPV 800", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-850.root", sys=-1.0, label="RPV 850", scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_RPV_2t6j_mStop-900.root", sys=-1.0, label="RPV 900", scale=scale),
    ]
    return Data, sgData, bgData

def plotRatio(data, histoName, year="2016", xTitle="", yTitle="", rebinx=-1.0, rebiny=-1.0, xmin=None, xmax=None):
    hData = data[0][0].getHisto(histoName, rebinx, rebiny, xmin, xmax)
    hMC = data[1][0].getHisto(histoName, rebinx, rebiny, xmin, xmax)
    for i in range(len(data[1])):
        if i == 0: continue
        h = data[1][i].getHisto(histoName, rebinx, rebiny, xmin, xmax)
        hMC.Add(h)

    c1 = ROOT.TCanvas( "c", "c", 1200, 1200)
    ROOT.gPad.Clear()
    ROOT.gStyle.SetOptStat("")
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.16)
    ROOT.gPad.SetTicks(1,1)
    ROOT.gStyle.SetPaintTextFormat("4.2f")
    #ROOT.gPad.SetLogy()

    ratio = hData.Clone()
    ratio.Divide(hMC)
    ratio.SetTitle(histoName)
    ratio.GetXaxis().SetTitle(xTitle)                                                                                                                              
    ratio.GetYaxis().SetTitle(yTitle)
    ratio.SetMaximum(1.5)
    ratio.SetMinimum(0.5)
    ratio.GetXaxis().SetTitleSize(0.045)
    ratio.GetYaxis().SetTitleSize(0.045)
    ratio.GetYaxis().SetLabelSize(0.036)
    ratio.GetYaxis().SetTitleOffset(0.9)
    ratio.Draw("colz text")
    
    c1.SaveAs("outputPlots/FullRun2_"+year+"/"+histoName+".pdf")
    del c1                

def plotStack(data, histoName, outputPath="./", xTitle="", yTitle="", rebinx=-1.0, rebiny=-1.0, xmin=None, xmax=None):
    hMC = data[1][0].getHisto(histoName, rebinx, rebiny, xmin, xmax)
    for i in range(len(data[1])):
        if i == 0: continue
        h = data[1][i].getHisto(histoName, rebinx, rebiny, xmin, xmax)
        hMC.Add(h)

    c1 = ROOT.TCanvas( "c", "c", 1200, 1200)
    ROOT.gPad.Clear()
    ROOT.gStyle.SetOptStat("")
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.16)
    ROOT.gPad.SetTicks(1,1)
    ROOT.gStyle.SetPaintTextFormat("4.2f")
    #ROOT.gPad.SetLogy()

    hMC.SetTitle(histoName)
    hMC.GetXaxis().SetTitle(xTitle)                                                                                                                              
    hMC.GetYaxis().SetTitle(yTitle)
    hMC.SetMaximum(1.5)
    hMC.SetMinimum(0.5)
    hMC.GetXaxis().SetTitleSize(0.045)
    hMC.GetYaxis().SetTitleSize(0.045)
    hMC.GetYaxis().SetLabelSize(0.036)
    hMC.GetYaxis().SetTitleOffset(0.9)
    hMC.Draw("colz text")
    
    c1.SaveAs(outputPath+"/"+histoName+".pdf")
    del c1                

def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-y', dest='year', type='string', default='2016', help="Can pass in the run year")
    options, args = parser.parse_args()
    
    cuts = ["1l_HT300_ge7j_ge1b_Mbl"]
    Data, sgData, bgData = getData("condor/Analyze1Lep_"+options.year+"_v1.0/", 1.0, options.year)
    
    for cut in cuts:
        #plotRatio((Data, bgData), "h_jEta_jPhi_"+cut, options.year, "#eta_{j}", "#phi_{j}", 8, 12, -2.4, 2.4)
        #plotRatio((Data, bgData), "h_lEta_lPhi_"+cut, options.year, "#eta_{l}", "#phi_{l}", 8, 12, -2.4, 2.4)
        plotStack((Data, bgData), "h_jEta_jPhi_"+cut, options.year, "#eta_{j}", "#phi_{j}", 8, 12, -2.4, 2.4)

if __name__ == '__main__':
    main()
