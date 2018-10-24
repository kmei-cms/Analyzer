import ROOT
ROOT.gStyle.SetOptStat(0)
from math import sqrt
import random

class MakeNJetPlots:
    def __init__(self):    
        self.colors = [ROOT.kCyan+1, ROOT.kMagenta+1, ROOT.kYellow+1,
                       ROOT.kRed+1, ROOT.kGreen+1, ROOT.kBlue+1,
                       ROOT.kBlack]
        self.signalcolors = [ROOT.kViolet-6, ROOT.kTeal-6, ROOT.kPink-6, ROOT.kAzure-6]
        
    def calcUnc(self,a,unc_a,b,unc_b):
        return  a/b * sqrt( (unc_a/a)**2+(unc_b/b)**2)
    
    def makeNJplot(self,inputfiles, labels, signalfiles, signallabels, jettype, cut, outputfile):
        basename = "h_njets_" + jettype
        histos = [inputfile.Get(basename+"_"+cut) for inputfile in inputfiles]
        # scale histos
        for i,h in enumerate(histos):
            #h.Scale(35900.)
            h.SetName(h.GetName()+"_"+labels[i])
            h.SetTitle(h.GetTitle()+"_"+labels[i])
            h.Write()
    
        signalhistos = []
        if signalfiles:
            signalhistos = [inputfile.Get(basename+"_"+cut) for inputfile in signalfiles]
            for i,h in enumerate(signalhistos):
                #h.Scale(35900.)
                h.SetName(h.GetName()+"_"+signallabels[i])
                h.SetTitle(h.GetTitle()+"_"+signallabels[i])
                h.Write()
                
        #make some pseudo_data
        mynewh = ROOT.TH1D(basename+"_"+cut+"_pseudodata", basename+"_"+cut+"_pseudodata",
                           histos[0].GetNbinsX(), histos[0].GetBinLowEdge(1), histos[0].GetBinLowEdge(1)+histos[0].GetNbinsX())
        #make some pseudo_data with signal
        pseudodataS_histos = []
        for sig in signallabels:
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
    def __init__(self, nMVABins, outFile):
        self.nMVABins = nMVABins
        self.outFile = outFile

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
    # I hadd my ttbar files into TT.root, and I hadd all other backgrounds into BG_noTT.root
    inputfilenames = [basedir + "TT.root", basedir + "BG_noTT.root"]
    inputfiles = [ROOT.TFile.Open(inputfilename) for inputfilename in inputfilenames]
    labels = ["TT", "other"]

    signalfilenames = [basedir+"stealth_stop_350_SYY.root", 
                       basedir+"stealth_stop_450_SYY.root", 
                       basedir+"stealth_stop_550_SYY.root", 
                       basedir+"stealth_stop_650_SYY.root", 
                       basedir+"stealth_stop_750_SYY.root", 
                       basedir+"stealth_stop_850_SYY.root", 
                       basedir+"rpv_stop_350.root", 
                       basedir+"rpv_stop_450.root", 
                       basedir+"rpv_stop_550.root", 
                       basedir+"rpv_stop_650.root", 
                       basedir+"rpv_stop_750.root", 
                       basedir+"rpv_stop_850.root", 
                       ]
    signalfiles = [ROOT.TFile.Open(signalfilename) for signalfilename in signalfilenames]
    signallabels = ["SYY_350","SYY_450","SYY_550","SYY_650","SYY_750","SYY_850", 
                    "RPV_350","RPV_450","RPV_550","RPV_650","RPV_750","RPV_850",]

    outputfile = ROOT.TFile.Open("njets_for_Aron.root","RECREATE")

    cutlevels = [
        "1l_deepESMbin1", "1l_deepESMbin2","1l_deepESMbin3","1l_deepESMbin4",                 
        ]
    jettypes = ["pt30"]#, "pt45"]

    mp = MakeNJetPlots()
    for jettype in jettypes:
        for cut in cutlevels:
            mp.makeNJplot(inputfiles, labels, signalfiles, signallabels, jettype, cut, outputfile)

    outputfile.Close()
    for inputfile in inputfiles:
        inputfile.Close()

    #Make data card
    md = MakeDataCard(nMVABins = 4, outFile = "multiF_ESM_ws.root")
    md.printDataCard()
