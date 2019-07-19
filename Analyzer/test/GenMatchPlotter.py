import ROOT
import math
import numpy as np

def GenWRecoplot(path, DR, resPart, filename, cut):

    c = ROOT.TCanvas( "c", "c", 0, 0, 1000, 1000)
    ROOT.TH1.AddDirectory(0)
   # ROOT.gPad.SetLeftMargin(0.1)
   # ROOT.gPad.SetRightMargin(0.1)
   # ROOT.gPad.SetTopMargin(0.08)
   # ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)
    ROOT.TH1.SetDefaultSumw2()
    file = ROOT.TFile.Open("condor/"+path+"/"+filename+".root")
    AllGenFile = ROOT.TFile.Open("condor/2016_RT_DRle5_MC_Samples/"+filename+".root")
    AllGen = AllGenFile.Get("h_"+resPart+"GenMass_"+cut)
    AllGen.SetLineColor(ROOT.kGreen)
    AllGen.Rebin(4)
    AllGen.SetTitle(resPart+", "+filename+ ", DeltaR < "+DR)
    AllGen.GetXaxis().SetTitle("Mass [GeV]")
    AllGen.SetStats(0)
    Gen = file.Get("h_"+resPart+"GenMass_")
    Gen.SetLineColor(ROOT.kRed)
    Gen.GetXaxis().SetTitle("Mass [GeV]")
    Gen.GetYaxis().SetTitle("Events")
    Gen.SetStats(0)
    Gen.Rebin(4)
#    Gen.SetTitle(resPart+" RPV 550,  DeltaR = "+DR)
    legend = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legend.AddEntry(Gen, "Matched Gen Particles", "l")
    legend.AddEntry(AllGen, "All Gen Particles", "l")
    Reco = file.Get("h_"+resPart+"Mass_")
    Reco.SetFillColor(ROOT.kBlue)
    Reco.SetTitle(resPart+", "+filename+ ", DeltaR < "+DR)
    Reco.GetXaxis().SetTitle("Mass [GeV]")
#    Reco.GetYaxis().SetTitle("Events")
    Reco.SetStats(0)
    Reco.Rebin(4)
    legend.AddEntry(Reco, "Reconstructed jets and leptons", "l")
    if resPart == "Nlino1" or resPart == "Nlino2" or resPart == "Single1" or resPart == "Single2":
        Gen.GetXaxis().SetRangeUser(0,500)
        Reco.GetXaxis().SetRangeUser(0,500)
        AllGen.GetXaxis().SetRangeUser(0,500)
    if filename == "2016_StealthSYY_2t6j_mStop-900" and (resPart == "Stop1" or resPart == "Stop2" or resPart == "StopMT2"):
        Reco.GetXaxis().SetRangeUser(0,1500)
        Gen.GetXaxis().SetRangeUser(0,1500)
        AllGen.GetXaxis().SetRangeUser(0,1500)
    max = Reco.GetMaximum(10000)
   # Reco.SetMaximum(max*2)
    Reco.GetYaxis().SetRangeUser(0,max*2)
    Reco.Draw("h")
    AllGen.Draw("same h")
    Gen.Draw("same h")
   # Gen.Draw("same l")
    legend.Draw("same")

    c.SaveAs("ResTest_"+filename+"_"+resPart+"_DR"+DR+"_"+cut+".pdf")
    del c

def GenPlot(path, DR, resPart, filename):
    c = ROOT.TCanvas( "c", "c", 0, 0, 1000, 1000)
    file = ROOT.TFile.Open("condor/"+path+"/"+filename+".root")
    Gen = file.Get("h_"+resPart+"GenMass_")
    Gen.SetLineColor(ROOT.kRed)
    Gen.SetTitle("Gen Particles")
    Gen.GetXaxis().SetTitle("Mass [GeV]")
    Gen.GetYaxis().SetTitle("Events")
    Gen.SetStats(0)
    Gen.Rebin(2)
#    Gen.SetTitle(resPart+" RPV 550,  DeltaR = "+DR)
    legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
    legend.AddEntry(Gen, "Gen Particles", "l")
    if resPart == "Nlino1" or resPart == "Nlino2" or resPart == "Single1" or resPart == "Single2":
        Gen.GetXaxis().SetRangeUser(0,500)
    Gen.Draw()
 #   legend.Draw("same")
    c.SaveAs("ResTest_Gen_"+filename+"_"+resPart+"_DR"+DR+".pdf")
    del c


def main():

    pathList = ["2016_RT_DRle1_MC_Samples", "2016_RT_DRle2_MC_Samples", "2016_RT_DRle5_MC_Samples"]
    DRList = ["1", "2"]
    resPartList = ["Stop1", "Stop2", "Nlino1", "Nlino2", "Single1", "Single2", "StopMT2"]
    fileNames = ["2016_RPV_2t6j_mStop-350", "2016_RPV_2t6j_mStop-550", "2016_StealthSYY_2t6j_mStop-900"]
    cuts = ["", "_0l_HT500_ge2b_ge2t", "_1l_base" , "_2l_ge1b_BothMblge25le250_HTge200"]

    for r in resPartList:
        for f in fileNames:
            for d in DRList:
                for c in cuts:
                    GenWRecoplot("2016_RT_DRle"+d+"_MC_Samples", d, r, f, c)
   # for f in fileNames:
       # for d in DRList:
           # MT2Plot("2016_RT_DRle"+d+"_MC_Samples", d, "Stop", f)
                
           # GenPlot("2016_RT_DRle5_MC_Samples", "5", r, f)
               
  # plot("2016_RT_DRle02_MC_Samples", "01", "Nlino1")

if __name__ == '__main__':
    main()
