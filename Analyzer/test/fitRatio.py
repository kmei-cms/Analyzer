import ROOT
import math
import numpy as np
import DataSetInfo as info

def main():
    #path = "condor/Analyze1Lep_Kerasv1.2.0/"
    path = "condor/Analyze1Lep_Kerasv1.2.0_MCTrigger_bTag_leptonWeight_ht300/"
    histoNames = [
        #"blind_ht_1l_4j_ge1b", "blind_ht_1l_5j_ge1b", "blind_ht_1l_6j_ge1b", "blind_ht_1l_ge7j_ge1b",
        "blind_ht_1l_5j_ge1b", "blind_ht_1l_6j_ge1b", "blind_ht_1l_ge7j_ge1b",
        #"blind_deepESM_1l_4j_ge1b", "blind_deepESM_1l_5j_ge1b", "blind_deepESM_1l_6j_ge1b", "blind_deepESM_1l_ge7j_ge1b",
                  ]
    nJets = []
    norm = []
    expo = []

    nJet = 5
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
        fit = ROOT.TF1("fit", "[0]*exp([1]*x)", 300.0, 3000.0, 2)
        #fit.SetParLimits(0,  -1.0, 1.0)
        #fit.SetParLimits(1,  -2.0, 2.0)
        #fit.SetParameter(0, -1.0)
        #fit.SetParameter(1,  0.0)        
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
        norm.append(fit.GetParameter(0))
        expo.append(fit.GetParameter(1))

        c.SaveAs(histoName+"_ratio.png")
        del c           
        nJet+=1

    ####################################################
    #Fit results of all the fits
    ####################################################
    nJets = np.array(nJets)
    norm = np.array(norm)
    expo = 1000*np.array(expo)
    print nJets
    print norm
    print expo
    
    s = ROOT.TGraph(nJets.size, nJets.astype(np.double), norm.astype(np.double))
    y = ROOT.TGraph(nJets.size, nJets.astype(np.double), expo.astype(np.double))
    
    List = [("Norm Term",s, 0.0, 2.0, -5, 5, -1),
            #("expo",y, -0.0005, 0.0, -0.002, 0.0, -0.0001)]
            ("Exp Term",y, -1, 0.0, -5, 5, -1)]
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
        fit.SetParLimits(0,  l[4], l[5])
        fit.SetParLimits(1,  -5.0, 5.0)
        fit.SetParameter(0, l[6])
        fit.SetParameter(1,  0.0)        
        fit.SetLineColor(ROOT.kRed)
        l[1].Fit(fit, "M", "", 0.0, 10.0)
        fit.Draw("same")
        
        c.SaveAs(l[0]+"_fits.png")
        del c           

    ####################################################
    #Plot Scale Factor per nJet
    ####################################################
    c = ROOT.TCanvas( "c", "c", 0, 0, 800, 800)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)
    ROOT.TH1.SetDefaultSumw2()
    ROOT.gStyle.SetOptStat(0);    

    leg = ROOT.TLegend(0.6, 0.65, 0.9, 0.88)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(1)
    leg.SetTextFont(42)
        
    h = ROOT.TH1F("Ht Scale Factor", "Ht Scale Factor", 1,0,5000)
    h.GetXaxis().SetTitle("Ht [GeV]")
    h.GetYaxis().SetTitle("Scale Factor")
    h.Set
    h.Draw()

    dummyList = []
    for nJet in [(ROOT.kCyan+1, "7") , (ROOT.kCyan+2, "8"), (ROOT.kCyan+3, "9"), (ROOT.kCyan+4, "10"), (ROOT.kRed, "11"), (ROOT.kRed+1, "12"), (ROOT.kRed+2, "13"), (ROOT.kRed+3, "14")]:        
        sf = ROOT.TF1("sf"+nJet[1], "(0.06146*"+nJet[1]+"+0.7908)*exp((-0.06063*"+nJet[1]+"+0.1018)*(x/1000))", 0.0, 5000)
        sf.SetLineColor(nJet[0])
        #sf.SetMarkerColor(nJet[0])
        sf.DrawCopy("l same")

        hd = ROOT.TH1F("Ht Scale Factor"+nJet[1], "Ht Scale Factor"+nJet[1], 1, 0, 5000)
        hd.Draw("P same")
        hd.SetLineColor(nJet[0])
        hd.SetLineWidth(3)
        leg.AddEntry(hd, "NJet "+nJet[1], "l")
        dummyList.append(hd)

    leg.Draw()
    c.SaveAs("sf.png")
        
if __name__ == '__main__':
    main()
