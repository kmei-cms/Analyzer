import ROOT
import math
import DataSetInfo as info

def getSig(sgData, bgData, histoName):
    # Loop over all background and get their histograms
    bgHistos = []
    for sample in bgData:
        h = sample.getHisto(histoName)
        bgHistos.append( (h, sample.sys) )

    # Loop over all signal models and calculate sig. for each one 
    sigVec = []
    for sample in sgData:
        hSg = sample.getHisto(histoName)
            
        sigTot2 = 0.0
        for nJetBin in range(7, hSg.GetNbinsX()):            
            nSg = hSg.GetBinContent(nJetBin)
            nBg = 0.0
            
            # Loop over all background histos and calculate sigma
            sigma2 = 0.0
            for t in bgHistos:
                n = t[0].GetBinContent(nJetBin)
                nBg += n
                sys = ((t[1])*n)**2
                sigma2 += n + sys

            # Calculate and fill sig.
            sig = 0.0
            if(sigma2 > 0 and nSg >= 1.0 and nBg >= 1.0): 
                sig = nSg/(math.sqrt(sigma2))
            sigTot2 += sig**2

        sigTot = math.sqrt(sigTot2)
        sigVec.append(sigTot)
        #print sample.label, sigTot

    return sigVec

def getData(path, scale=1.0, year = "2018"):
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
    return sgData, bgData

def main():
    XMin = 0.0; XMax = 1.0; RatioXMin = 0.0; RatioXMax = 1.0
    YMin = 0.3; YMax = 1.0; RatioYMin = 0.0; RatioYMax = 0.3
    PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)

    XTitleSize = 0.045;  XLabelSize = 0.036;  XTitleOffset = 4.0 
    YTitleSize = 0.045;  YLabelSize = 0.036;  YTitleOffset = 0.9

    c1 = ROOT.TCanvas( "c", "c", 1200, 1200)
    ROOT.gPad.Clear()
    #c1.Divide(1,2)
    #c1.cd(1)
    #ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)
    #ROOT.gStyle.SetStatY(0.5)
    #ROOT.gStyle.SetStatX(0.85)
    ROOT.gStyle.SetOptStat("")
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetTopMargin(0.08)
    #ROOT.gPad.SetBottomMargin(0.01)
    ROOT.gPad.SetBottomMargin(0.16)
    ROOT.gPad.SetTicks(1,1)
    #ROOT.gPad.SetLogy()

    iamLegend = ROOT.TLegend(0.55, 0.65, 0.85, 0.88)
    iamLegend.SetFillStyle(0)
    iamLegend.SetBorderSize(0)
    iamLegend.SetLineWidth(1)
    iamLegend.SetNColumns(1)
    iamLegend.SetTextFont(42)

    #Get the sig. for 2018pre
    sgData2018pre, bgData2018pre = getData("condor/Analyze1Lep_2018pre_v1.0/", 1.0, "2018pre")
    sigVec2018pre = getSig(sgData2018pre, bgData2018pre, "h_njets_1l_HT300_ge7j_ge1b_Mbl")
    sig2018pre = ROOT.TH1D("hist","hist", 13, 300, 900)
    sig2018pre.SetTitle("RPV: Simple Significance")
    #sig2018pre.GetXaxis().SetTitle("")
    sig2018pre.GetXaxis().SetTitle("Stop Mass [GeV]")                                                                                                                              
    sig2018pre.GetYaxis().SetTitle("Significance")
    sig2018pre.SetLineColor(ROOT.kCyan)
    sig2018pre.SetMaximum(50)
    sig2018pre.SetMinimum(0.05)
    sig2018pre.GetXaxis().SetTitleSize(XTitleSize)
    sig2018pre.GetYaxis().SetTitleSize(YTitleSize)
    #sig2018pre.GetXaxis().SetLabelSize(0)
    sig2018pre.GetYaxis().SetLabelSize(YLabelSize)
    #sig2018pre.GetXaxis().SetTitleOffset(XTitleOffset)
    sig2018pre.GetYaxis().SetTitleOffset(YTitleOffset)
    for i in range(1, len(sigVec2018pre)+1):
        sig2018pre.SetBinContent(i, sigVec2018pre[i-1])

    #Get the sig. for 2018post
    sgData2018post, bgData2018post = getData("condor/Analyze1Lep_2018post_v1.0/", 1.0, "2018post")
    sigVec2018post = getSig(sgData2018post, bgData2018post, "h_njets_1l_HT300_ge7j_ge1b_Mbl")
    sig2018post = ROOT.TH1D("hist","hist", 13, 300, 900)
    sig2018post.SetLineColor(ROOT.kRed)
    for i in range(1, len(sigVec2018post)+1):
        sig2018post.SetBinContent(i, sigVec2018post[i-1])

    #Get the sig. for 2017
    sgData2017, bgData2017 = getData("condor/Analyze1Lep_2017_v1.0/", 1.0, "2017")
    sigVec2017 = getSig(sgData2017, bgData2017, "h_njets_1l_HT300_ge7j_ge1b_Mbl")
    sig2017 = ROOT.TH1D("hist","hist", 13, 300, 900)
    sig2017.SetLineColor(ROOT.kBlue)
    for i in range(1, len(sigVec2017)+1):
        sig2017.SetBinContent(i, sigVec2017[i-1])

    #Get the sig. for 2016
    sgData2016, bgData2016 = getData("condor/Analyze1Lep_2016_v1.0/", 1.0, "2016")
    sigVec2016 = getSig(sgData2016, bgData2016, "h_njets_1l_HT300_ge7j_ge1b_Mbl")
    sig2016 = ROOT.TH1D("hist","hist", 13, 300, 900)
    sig2016.SetLineColor(ROOT.kGreen + 2)
    for i in range(1, len(sigVec2016)+1):
        sig2016.SetBinContent(i, sigVec2016[i-1])

    #Get total sig. for all 4 periods
    sigTot = sig2018pre.Clone()
    sigTot.SetLineColor(ROOT.kBlack)
    for i in range(1, len(sigVec2016)+1):
        sigTot.SetBinContent(i, math.sqrt(sigVec2018pre[i-1]**2 + sigVec2018post[i-1]**2 + sigVec2017[i-1]**2 + sigVec2016[i-1]**2))

    sig2018pre.Draw("hist")
    sig2018post.Draw("hist same")
    sig2017.Draw("hist same")
    sig2016.Draw("hist same")
    sigTot.Draw("hist same")

    iamLegend.AddEntry(sig2018pre,  "2018pre", "L")
    iamLegend.AddEntry(sig2018post, "2018post", "L")
    iamLegend.AddEntry(sig2017, "2017", "L")
    iamLegend.AddEntry(sig2016, "2016", "L")
    iamLegend.AddEntry(sigTot, "Total", "L")
    iamLegend.Draw("same")

    ##Make Ratio Plot
    #c1.cd(2)
    #ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)
    #ROOT.gPad.SetLeftMargin(0.12)
    #ROOT.gPad.SetRightMargin(0.15)
    #ROOT.gPad.SetTopMargin(0.01)
    #ROOT.gPad.SetBottomMargin(0.26)
    #ROOT.gPad.SetTicks(1,1)
    #ROOT.gPad.SetGridy()
    #
    #ratio = sigHEMFullYear.Clone()
    #ratio.Divide(sigHEMPartYear)
    #ratio.SetTitle("")
    #ratio.GetXaxis().SetTitle("Stop Mass [GeV]")
    #ratio.GetYaxis().SetTitle("Rato of 2018 Options")
    #ratio.SetMaximum(1.25)
    #ratio.SetMinimum(0.05)
    #ratio.GetXaxis().SetTitleSize(PadFactor*XTitleSize); 
    #ratio.GetYaxis().SetTitleSize(0.8*PadFactor*YTitleSize)
    #ratio.GetXaxis().SetLabelSize(PadFactor*XLabelSize); 
    #ratio.GetYaxis().SetLabelSize(PadFactor*YLabelSize)
    #ratio.GetXaxis().SetTitleOffset(0.65*XTitleOffset/PadFactor)
    #ratio.GetYaxis().SetTitleOffset(1.3*YTitleOffset/PadFactor); 
    #
    #ratio.Draw()

    c1.SaveAs("sig.pdf")
    del c1                

if __name__ == '__main__':
    main()
