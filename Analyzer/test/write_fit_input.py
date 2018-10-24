import ROOT
ROOT.gStyle.SetOptStat(0)
from math import sqrt
import random

class DataSetInfo:
    def __init__(self, basedir, fileName, label):
        self.basedir = basedir
        self.fileName = fileName
        self.label = label
        self.file = ROOT.TFile.Open(basedir+fileName)

    def __del__(self):
        self.file.Close()

class WriteNJetPlots:
    def __init__(self):    
        self.colors = [ROOT.kCyan+1, ROOT.kMagenta+1, ROOT.kYellow+1,
                       ROOT.kRed+1, ROOT.kGreen+1, ROOT.kBlue+1,
                       ROOT.kBlack]
        self.signalcolors = [ROOT.kViolet-6, ROOT.kTeal-6, ROOT.kPink-6, ROOT.kAzure-6]
        
    def calcUnc(self,a,unc_a,b,unc_b):
        return  a/b * sqrt( (unc_a/a)**2+(unc_b/b)**2)
    
    def writeHistos(self, dataList, basename, cut):
        histos = []
        for dsi in dataList:
            h = dsi.file.Get(basename+"_"+cut)
            h.SetName(h.GetName()+"_"+dsi.label)
            h.SetTitle(h.GetTitle()+"_"+dsi.label)
            h.Write()
            histos.append(h)
        return histos

    def makePseudoData(self, histos, signalhistos, sgDataList, jettype, cut):
        #make some pseudo_data
        mynewh = ROOT.TH1D(basename+"_"+cut+"_pseudodata", basename+"_"+cut+"_pseudodata",
                           histos[0].GetNbinsX(), histos[0].GetBinLowEdge(1), histos[0].GetBinLowEdge(1)+histos[0].GetNbinsX())
        
        #make some pseudo_data with signal
        pseudodataS_histos = []
        for dsi in sgDataList:
            sig = dsi.label
            pseudodataS_histos.append( ROOT.TH1D(basename+"_"+cut+"_pseudodataS_"+sig, basename+"_"+cut+"_pseudodataS_"+sig,
                                                 histos[0].GetNbinsX(), histos[0].GetBinLowEdge(1), histos[0].GetBinLowEdge(1)+histos[0].GetNbinsX())
                                       )
        for bin in range(histos[0].GetNbinsX()):
            sumdata = 0
            #sumunc = 0
            for h in histos:
                sumdata += h.GetBinContent(bin+1)
                #sumunc += h.GetBinError(bin+1)
            #sumdata = sumdata + sumunc*2*(random.random()-0.5)
            mynewh.SetBinContent(bin+1, int(round(sumdata)))
        
            for i, sigh in enumerate(signalhistos):
                sumsig = sigh.GetBinContent(bin+1)
                #sumsigunc = sigh.GetBinError(bin+1)
                sumdatasigh = sumdata + sumsig #+ sumsigunc*2*(random.random()-0.5)
                pseudodataS_histos[i].SetBinContent(bin+1, int(round(sumdatasigh)))
        
        mynewh.Write()
        for h in pseudodataS_histos:
            h.Write()

class MakeDataCard:
    def __init__(self, nMVABins, outFile, bgDataList, sgDataList):
        self.nMVABins = nMVABins
        self.outFile = outFile
        self.bgDataList = bgDataList
        self.sgDataList = sgDataList

    def printDataCard(self):
        print "Signal Region 1 Datacard -- signal category"
        print ""
        print "imax * number of bins"
        print "jmax * number of processes minus 1"
        print "kmax * number of nuisance parameters"
        print ""
        print "-------------------------------------------------------------------------------------------------------------------------------------------"
        print ""
        for i in range(self.nMVABins):
            print "shapes data_obs    sigD{0}    {1} wspace:data_obs_D{0}".format(i+1, self.outFile)
            print "shapes bkg_tt      sigD{0}    {1} wspace:background_tt_D{0}".format(i+1, self.outFile)
            print "shapes bkg_other   sigD{0}    {1} wspace:background_other_D{0}".format(i+1, self.outFile)
            print "shapes signal      sigD{0}    {1} wspace:signal_D{0}".format(i+1, self.outFile)
            print ""
        print ""
        print "-------------------------------------------------------------------------------------------------------------------------------------------"
        bins        = "bin         "
        observation = "observation "
        for i in range(self.nMVABins):
            bins        += " sigD{0} ".format(i+1)
            observation += "  -1   "
        print bins
        print observation
        print "-------------------------------------------------------------------------------------------------------------------------------------------"
        print "# background rate must be taken from _norm param x 1"
        print "bin                 sigD1       sigD1       sigD1        sigD2       sigD2       sigD2        sigD3       sigD3       sigD3        sigD4       sigD4       sigD4"
        print "process             signal      bkg_tt      bkg_other    signal      bkg_tt      bkg_other    signal      bkg_tt      bkg_other    signal      bkg_tt      bkg_other"
        print "process             0           1           2            0           1           2            0           1           2            0           1           2"
        print "rate                46.6614     1           7289.53      133.046     1           5997.22      289.856     1           6255.23      986.308     1           5385.85"
        print "-------------------------------------------------------------------------------------------------------------------------------------------"
        print "# Normal uncertainties in the signal region"
        print "lumi_13TeV         lnN    1.05     -    -    1.05     -    -    1.05     -    -    1.05     -    -"
        print "-------------------------------------------------------------------------------------------------------------------------------------------"
        for i in range(self.nMVABins):
            print "ttBkgRateD{0} rateParam sigD{0} bkg_tt {1}:wspace".format(i+1, self.outFile)
        for i in range(self.nMVABins):
            print "np_tt_D{0} param 1.0 1.0".format(i+1)
        print "#np_tt   param 1.0 1.0"

if __name__ == "__main__":
    # Where the root files are stored
    basedir = "condor/rootfiles/"
    cutlevels = [ "1l_deepESMbin1", "1l_deepESMbin2","1l_deepESMbin3","1l_deepESMbin4"]
    jettypes = ["pt30"]#, "pt45"]
    outputfile = ROOT.TFile.Open("njets_for_Aron.root","RECREATE")

    # I hadd my ttbar files into TT.root, and I hadd all other backgrounds into BG_noTT.root
    bgDataList = [
        DataSetInfo(basedir, "TT.root",      "TT"),
        DataSetInfo(basedir, "BG_noTT.root", "other"),
    ]

    sgDataList = [
        DataSetInfo(basedir, "stealth_stop_350_SYY.root", "SYY_350"),
        DataSetInfo(basedir, "stealth_stop_450_SYY.root", "SYY_450"),
        DataSetInfo(basedir, "stealth_stop_550_SYY.root", "SYY_550"),
        DataSetInfo(basedir, "stealth_stop_650_SYY.root", "SYY_650"),
        DataSetInfo(basedir, "stealth_stop_750_SYY.root", "SYY_750"),
        DataSetInfo(basedir, "stealth_stop_850_SYY.root", "SYY_850"),
        DataSetInfo(basedir, "rpv_stop_350.root",         "RPV_350"),
        DataSetInfo(basedir, "rpv_stop_450.root",         "RPV_450"),
        DataSetInfo(basedir, "rpv_stop_550.root",         "RPV_550"),
        DataSetInfo(basedir, "rpv_stop_650.root",         "RPV_650"),
        DataSetInfo(basedir, "rpv_stop_750.root",         "RPV_750"),
        DataSetInfo(basedir, "rpv_stop_850.root",         "RPV_850"),
    ]

    #Write all histos to outputfile
    outputfile.cd()
    wp = WriteNJetPlots()
    for jettype in jettypes:
        basename = "h_njets_" + jettype
        for cut in cutlevels:
            histos = wp.writeHistos(bgDataList, basename, cut)
            signalhistos = wp.writeHistos(sgDataList, basename, cut)
            wp.makePseudoData(histos, signalhistos, sgDataList, jettype, cut)
    
    #Close outfile
    outputfile.Close()
    
    #Make data card
    md = MakeDataCard(nMVABins=4, outFile="multiF_ESM_ws.root", bgDataList=bgDataList, sgDataList=sgDataList)
    md.printDataCard()
