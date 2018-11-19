import ROOT
import math
import numpy as np
import DataSetInfo as info

def main():
    #path = "condor/Analyze1Lep_Kerasv1.2.0/"
    path = "condor/Analyze1Lep_Kerasv1.2.0_MCTrigger_bTag_leptonWeight_ht300/"
    histoNames = [
        "blind_ht_1l_4j_ge1b", "blind_ht_1l_5j_ge1b", "blind_ht_1l_6j_ge1b", "blind_ht_1l_ge7j_ge1b",
        #"blind_deepESM_1l_4j_ge1b", "blind_deepESM_1l_5j_ge1b", "blind_deepESM_1l_6j_ge1b", "blind_deepESM_1l_ge7j_ge1b",
                  ]
    nJets = []
    slope = []
    yInt = []

    nJet = 4
    for histoName in histoNames:
        c = ROOT.TCanvas( "c", "c", 0, 0, 800, 800)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTicks(1,1)
        ROOT.TH1.SetDefaultSumw2()
        ROOT.gStyle.SetOptFit();

        data = info.DataSetInfo(basedir=path, fileName="Data.root",   label="Data")
        bg   = info.DataSetInfo(basedir=path, fileName="AllBG.root",  label="Background")
        dHist = data.getHisto(histoName)
        bHist = bg.getHisto(histoName)
        dHist.Rebin(10)
        bHist.Rebin(10)
        dHist.Scale(1.0/dHist.Integral())
        bHist.Scale(1.0/bHist.Integral())

        ratio = dHist.Clone("Ratio")
        ratio.Divide(bHist)
        ratio.SetMaximum(2.0)
        ratio.SetMinimum(0.0)
        ratio.SetLineColor(ROOT.kBlack)
        ratio.SetMarkerColor(ROOT.kBlack)
        ratio.SetMarkerStyle(21)
        ratio.SetTitle(histoName+" Ratio")
        ratio.GetYaxis().SetTitle("Data/BG")
        ratio.Draw("ep")
        
        line = ROOT.TF1("1" ,"1" ,-2000,20000)
        line.SetLineColor(ROOT.kBlack)
        line.Draw("same")
        
        ratio.GetXaxis().SetTitle("H_{T} [GeV]")
        fit = ROOT.TF1("fit", "[0]*x + [1]", 300.0, 3000.0, 2)
        fit.SetParLimits(0,  -1.0, 1.0)
        fit.SetParLimits(1,  -2.0, 2.0)
        fit.SetParameter(0, -1.0)
        fit.SetParameter(1,  0.0)        
        fit.SetLineColor(ROOT.kRed)
        ratio.Fit(fit, "", "", 300.0, 3000.0)
        fit.Draw("same")

        #ratio.GetXaxis().SetTitle("DeepESM Score")
        #fit = ROOT.TF1("fit", "[0]*x + [1]", 0.0, 1.0, 2)
        #fit.SetParLimits(0,  -1.0, 1.0)
        #fit.SetParLimits(1,  -2.0, 2.0)
        #fit.SetParameter(0, -1.0)
        #fit.SetParameter(1,  0.0)        
        #fit.SetLineColor(ROOT.kRed)
        #ratio.Fit(fit, "", "", 0.0, 1.0)
        #fit.Draw("same")
        
        nJets.append(nJet)
        slope.append(fit.GetParameter(0))
        yInt.append(fit.GetParameter(1))

        c.SaveAs(histoName+"_ratio.png")
        del c           
        nJet+=1

    #Fit slope of all the fits
    print nJets
    print slope
    print yInt

    nJets = np.array(nJets)
    slope = np.array(slope)
    yInt = np.array(yInt)

    s = ROOT.TGraph(nJets.size, nJets.astype(np.double), slope.astype(np.double))
    y = ROOT.TGraph(nJets.size, nJets.astype(np.double), yInt.astype(np.double))

    List = [("slope",s, -0.001, 0.001),("yInt",y, 1, 1.25)]
    for l in List:
        print np.max(l[2])
        c = ROOT.TCanvas( "c", "c", 0, 0, 800, 800)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTicks(1,1)
        ROOT.TH1.SetDefaultSumw2()
        ROOT.gStyle.SetOptFit();
        
        l[1].SetMaximum(l[3])
        l[1].SetMinimum(l[2])
        l[1].SetLineColor(ROOT.kBlack)
        l[1].SetMarkerColor(ROOT.kBlack)
        l[1].SetMarkerStyle(21)
        l[1].SetTitle(l[0]+" Fits")
        l[1].GetXaxis().SetTitle("NJet")
        l[1].GetYaxis().SetTitle(l[0])
        l[1].Draw("AP")
        
        fit = ROOT.TF1("fit", "[0]*x + [1]", 0.0, 10.0, 3)
        fit.SetParLimits(0,  -10.0, 10.0)
        fit.SetParLimits(1,  -10.0, 10.0)
        fit.SetParameter(0, -1.0)
        fit.SetParameter(1,  0.0)        
        fit.SetLineColor(ROOT.kRed)
        l[1].Fit(fit, "", "", 0.0, 10.0)
        fit.Draw("same")
        
        c.SaveAs(l[0]+"_fits.png")
        del c           

if __name__ == '__main__':
    main()
