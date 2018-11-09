import ROOT
import math

class DataSetInfo:
    def __init__(self, path, filename, sys, dataType):
        self.path = path
        self.filename = filename
        self.sys = sys
        self.dataType = dataType
        self.file = ROOT.TFile.Open(path+filename)

    def getHisto(self, name):
        return self.file.Get(name)
        
    def __del__(self):
        self.file.Close()

def main():
    path = "condor/Analyze1Lep_Kerasv1.2.0/"
    histoNames = ["blind_ht_1l_5to6j_ge1b", "blind_ht_1l_ge7j_ge1b"]

    for histoName in histoNames:
        c = ROOT.TCanvas( "c", "c", 0, 0, 800, 800)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTicks(1,1)
        ROOT.TH1.SetDefaultSumw2()
        ROOT.gStyle.SetOptFit();

        data = DataSetInfo(path=path, filename="Data.root",  sys=-1, dataType="Data")
        bg   = DataSetInfo(path=path, filename="AllBG.root", sys=-1, dataType="Background")
        dHist = data.getHisto(histoName)
        bHist = bg.getHisto(histoName)
        dHist.Rebin(10)
        bHist.Rebin(10)
        dHist.Scale(1.0/dHist.Integral())
        bHist.Scale(1.0/bHist.Integral())

        print dHist.Integral(), bHist.Integral()
        
        ratio = dHist.Clone("Ratio")
        ratio.Divide(bHist)
        ratio.SetMaximum(2.0)
        ratio.SetMinimum(0.0)
        ratio.SetLineColor(ROOT.kBlack)
        ratio.SetMarkerColor(ROOT.kBlack)
        ratio.SetMarkerStyle(21)
        ratio.SetTitle(histoName+" Ratio")
        ratio.GetXaxis().SetTitle("H_{T} [GeV]")
        ratio.GetYaxis().SetTitle("Data/BG")
        ratio.Draw("ep")
        
        line = ROOT.TF1("1" ,"1" ,-2000,20000)
        line.SetLineColor(ROOT.kBlack)
        line.Draw("same")
        
        fit = ROOT.TF1("fit", "[0]*x + [1]", 300.0, 3000.0, 2)
        fit.SetParLimits(0,  -1.0, 1.0)
        fit.SetParLimits(1,  -2.0, 2.0)
        fit.SetParameter(0, -1.0)
        fit.SetParameter(1,  0.0)        
        fit.SetLineColor(ROOT.kRed)
        ratio.Fit(fit, "", "", 300.0, 3000.0)
        fit.Draw("same")
        
        c.SaveAs(histoName+"_ratio.png")
        del c                

if __name__ == '__main__':
    main()
