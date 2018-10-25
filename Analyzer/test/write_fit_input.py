import ROOT
ROOT.gStyle.SetOptStat(0)
from math import sqrt
import random

colors = [ROOT.kCyan+1, ROOT.kMagenta+1, ROOT.kYellow+1,
          ROOT.kRed+1, ROOT.kGreen+1, ROOT.kBlue+1,
          ROOT.kBlack]
signalcolors = [ROOT.kViolet-6, ROOT.kTeal-6, ROOT.kPink-6, ROOT.kAzure-6]

class DataSetInfo:
    def __init__(self, basedir, fileName, label, processName, process, rate, lumiSys):
        self.basedir = basedir
        self.fileName = fileName
        self.label = label
        self.file = ROOT.TFile.Open(basedir+fileName)
        self.processName = processName
        self.process = process
        self.rate = rate
        self.lumiSys = lumiSys

    def getHisto(self, name):
        return self.file.Get(name)

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
    
    def writeHistos(self, data, basename, cut):
        histos = []
        for key, dsi in data.iteritems():
            h = dsi.getHisto(basename+"_"+cut)
            h.SetName(h.GetName()+"_"+dsi.label)
            h.SetTitle(h.GetTitle()+"_"+dsi.label)
            h.Write()
            histos.append(h)
        return histos

    def makePseudoData(self, histos, signalhistos, sgData, jettype, cut):
        #make some pseudo_data
        mynewh = ROOT.TH1D(basename+"_"+cut+"_pseudodata", basename+"_"+cut+"_pseudodata",
                           histos[0].GetNbinsX(), histos[0].GetBinLowEdge(1), histos[0].GetBinLowEdge(1)+histos[0].GetNbinsX())
        
        #make some pseudo_data with signal
        pseudodataS_histos = []
        for key, dsi in sgData.iteritems():
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
    def __init__(self, outFile, bgData, otData, sgData, basename, cutlevels):
        self.nMVABins = len(cutlevels)
        self.outFile = outFile
        self.bgData = bgData
        self.otData = otData
        self.sgData = sgData
        self.basename = basename
        self.cutlevels = cutlevels
        bgMVAHistos = self.getMVAHists(bgData, cutlevels, basename)
        otMVAHistos = self.getMVAHists(otData, cutlevels, basename)
        sgMVAHistos = self.getMVAHists(sgData, cutlevels, basename)        
        self.dataSets = [(sgMVAHistos, sgData), (bgMVAHistos, bgData), (otMVAHistos, otData)]

    def getMVAHists(self, data, cutlevels, basename):
        histos = []
        for cut in self.cutlevels:
            h = data.getHisto(self.basename+"_"+cut)
            histos.append(h)
        return histos

    def writeDataCard(self, name):
        f = open(name, 'w')

        f.write( "Signal Region 1 Datacard -- signal category" )
        f.write( "\n" )
        f.write( "imax * number of bins\n" )
        f.write( "jmax * number of processes minus 1\n" )
        f.write( "kmax * number of nuisance parameters\n" )
        f.write( "\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "\n" )
        for i in range(self.nMVABins):
            f.write( "shapes data_obs    sigD{0}    {1} wspace:data_obs_D{0}\n".format(i+1, self.outFile) )
            f.write( "shapes bkg_tt      sigD{0}    {1} wspace:background_tt_D{0}\n".format(i+1, self.outFile) )
            f.write( "shapes bkg_other   sigD{0}    {1} wspace:background_other_D{0}\n".format(i+1, self.outFile) )
            f.write( "shapes signal      sigD{0}    {1} wspace:signal_D{0}\n".format(i+1, self.outFile) )
            f.write( "\n" )
        f.write( "\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        bins        = "bin         "
        observation = "observation "
        for i in range(self.nMVABins):
            bins        += " sigD{0} ".format(i+1)
            observation += "  -1   "
        f.write( bins+"\n" )
        f.write( observation+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "# background rate must be taken from _norm param x 1\n" )
        bin      = "bin              "
        process1 = "process          "
        process2 = "process          "
        rate     = "rate             "
        for i in range(self.nMVABins):
            for d in self.dataSets:
                r = 1
                if(d[1].rate):
                    r = d[0][i].Integral()
                bin      += "{0: <16} ".format("sigD"+str(i+1))
                process1 += "{0: <16} ".format(d[1].processName)
                process2 += "{0: <16} ".format(d[1].process)
                rate     += "{0: <16} ".format(r)
        f.write( bin+"\n" )
        f.write( process1+"\n" )
        f.write( process2+"\n" )
        f.write( rate+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "# Normal uncertainties in the signal region\n" )
        lumiSys = "lumi_13TeV      lnN  "
        for i in range(self.nMVABins):
            for d in self.dataSets:
                lumiSys += "{0: <5} ".format(d[1].lumiSys)
        f.write( lumiSys+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        for i in range(self.nMVABins):
            f.write( "ttBkgRateD{0} rateParam sigD{0} bkg_tt {1}:wspace\n".format(i+1, self.outFile) )
        for i in range(self.nMVABins):
            f.write( "np_tt_D{0} param 1.0 1.0\n".format(i+1) )
        f.write( "#np_tt   param 1.0 1.0\n" )

if __name__ == "__main__":
    # Where the root files are stored
    basedir = "condor/rootfiles/"
    cutlevels = [ "1l_deepESMbin1", "1l_deepESMbin2","1l_deepESMbin3","1l_deepESMbin4"]
    jettypes = ["pt30"]#, "pt45"]
    outputfile = ROOT.TFile.Open("njets_for_Aron.root","RECREATE")
    outputDataCard = "dataCard.txt"

    # I hadd my ttbar files into TT.root, and I hadd all other backgrounds into BG_noTT.root
    bgData = {
        "TT"    : DataSetInfo(basedir=basedir, fileName="TT.root",      label="TT",    processName="bkg_tt",    process="1", rate=False, lumiSys="-"),
        "other" : DataSetInfo(basedir=basedir, fileName="BG_noTT.root", label="other", processName="bkg_other", process="2", rate=True,  lumiSys="-"),
    }

    sgData = {
        "SYY_350" : DataSetInfo(basedir=basedir, fileName="stealth_stop_350_SYY.root", label="SYY_350", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_450" : DataSetInfo(basedir=basedir, fileName="stealth_stop_450_SYY.root", label="SYY_450", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_550" : DataSetInfo(basedir=basedir, fileName="stealth_stop_550_SYY.root", label="SYY_550", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_650" : DataSetInfo(basedir=basedir, fileName="stealth_stop_650_SYY.root", label="SYY_650", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_750" : DataSetInfo(basedir=basedir, fileName="stealth_stop_750_SYY.root", label="SYY_750", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "SYY_850" : DataSetInfo(basedir=basedir, fileName="stealth_stop_850_SYY.root", label="SYY_850", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_350" : DataSetInfo(basedir=basedir, fileName="rpv_stop_350.root",         label="RPV_350", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_450" : DataSetInfo(basedir=basedir, fileName="rpv_stop_450.root",         label="RPV_450", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_550" : DataSetInfo(basedir=basedir, fileName="rpv_stop_550.root",         label="RPV_550", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_650" : DataSetInfo(basedir=basedir, fileName="rpv_stop_650.root",         label="RPV_650", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_750" : DataSetInfo(basedir=basedir, fileName="rpv_stop_750.root",         label="RPV_750", processName="signal", process="0", rate=True, lumiSys="1.05"),
        "RPV_850" : DataSetInfo(basedir=basedir, fileName="rpv_stop_850.root",         label="RPV_850", processName="signal", process="0", rate=True, lumiSys="1.05"),
    }

    #Write all histos to outputfile
    outputfile.cd()
    wp = WriteNJetPlots()
    for jettype in jettypes:
        basename = "h_njets_" + jettype
        for cut in cutlevels:
            histos = wp.writeHistos(bgData, basename, cut)
            signalhistos = wp.writeHistos(sgData, basename, cut)
            wp.makePseudoData(histos, signalhistos, sgData, jettype, cut)
    
    #Close outfile
    outputfile.Close()
    
    #Make data card
    md = MakeDataCard(outFile="multiF_ESM_ws.root", bgData=bgData["TT"], otData=bgData["other"], sgData=sgData["RPV_550"], basename="h_njets_pt30", cutlevels=cutlevels)
    md.writeDataCard(outputDataCard)

    import os
    os.system("cat "+outputDataCard)
