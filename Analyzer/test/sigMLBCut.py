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
    #path = "condor/Analyze1Lep_Kerasv1.2.0_noMCTrigger_PU_bTag_leptonWeight/"
    #histoName = "h_mbl_1l_ge7j_ge1b_noMbl"
    histoName = "h_njets_mbl_1l_ge7j_ge1b_noMbl"
    rebinVal = 5

    bgData = {
        "TT"              : DataSetInfo(path=path, filename="TT.root",              sys=0.3, dataType="Background"),
        #"TTJets"          : DataSetInfo(path=path, filename="TTJets.root",          sys=0.3, dataType="Background"),
        "QCD"             : DataSetInfo(path=path, filename="QCD.root",             sys=1.0, dataType="Background"),
        #"DYJetsToLL_M-50" : DataSetInfo(path=path, filename="DYJetsToLL_M-50.root", sys=1.0, dataType="Background"),
        #"Rare"            : DataSetInfo(path=path, filename="Rare.root",            sys=1.0, dataType="Background"),
        #"Diboson"         : DataSetInfo(path=path, filename="Diboson.root",         sys=1.0, dataType="Background"),
        #"ST"              : DataSetInfo(path=path, filename="ST.root",              sys=1.0, dataType="Background"),
        #"WJetsToLNu"      : DataSetInfo(path=path, filename="WJetsToLNu.root",      sys=1.0, dataType="Background"),
    }

    sgData = {
        "rpv_stop_350" : DataSetInfo(path=path, filename="rpv_stop_350.root", sys=-1.0, dataType="Signal"),
        #"rpv_stop_450" : DataSetInfo(path=path, filename="rpv_stop_450.root", sys=-1.0, dataType="Signal"),
        #"rpv_stop_550" : DataSetInfo(path=path, filename="rpv_stop_550.root", sys=-1.0, dataType="Signal"),
        #"rpv_stop_650" : DataSetInfo(path=path, filename="rpv_stop_650.root", sys=-1.0, dataType="Signal"),
        #"rpv_stop_750" : DataSetInfo(path=path, filename="rpv_stop_750.root", sys=-1.0, dataType="Signal"),
        #"rpv_stop_850" : DataSetInfo(path=path, filename="rpv_stop_850.root", sys=-1.0, dataType="Signal"),
    }

    # Loop over all background and get their histograms
    bgHistos = {}
    for key in bgData:
        h = bgData[key].getHisto(histoName)
        h.RebinY(rebinVal)
        bgHistos[key] = {"hist" : h}

    # Loop over all signal models and calculate sig. for each one 
    for key in sgData:
        hSg = sgData[key].getHisto(histoName)
        hSg.RebinY(rebinVal)
        nBins = hSg.GetNbinsY()
        low = hSg.GetYaxis().GetXmin()
        high = hSg.GetYaxis().GetXmax()
        binToGeV = (high - low)/nBins
        sigDic = {}
            
        # Find all possible MLB cut values
        for nJetBin in range(6, hSg.GetNbinsX()+1):            
            c2 = ROOT.TCanvas( "c", "c", 0, 0, 800, 800)
            ROOT.gStyle.SetStatY(0.5)
            ROOT.gStyle.SetStatX(0.85)
            ROOT.gStyle.SetPalette(ROOT.kRainBow)
            ROOT.gPad.SetLeftMargin(0.12)
            ROOT.gPad.SetRightMargin(0.15)
            ROOT.gPad.SetTopMargin(0.08)
            ROOT.gPad.SetBottomMargin(0.12)
            ROOT.gPad.SetTicks(1,1)
            sig2DHisto = ROOT.TH2D(key+"hist"+str(nJetBin),key+"hist"+str(nJetBin), nBins, low, high, nBins, low, high)
            for lowBin in range(1, nBins+1):       
                for highBin in range(nBins+1):
                    if(lowBin < highBin):
                        nSg = hSg.Integral(nJetBin, nJetBin, lowBin, highBin)
                        
                        # Loop over all background histos and calculate sigma
                        sigma2 = 0.0
                        for k in bgHistos:
                            n = bgHistos[k]["hist"].Integral(nJetBin, nJetBin, lowBin, highBin)
                            sys = ((bgData[k].sys)*n)**2
                            sigma2 += n + sys
            
                        # Calculate and fill sig.
                        sig = 0.0
                        if(sigma2 != 0): sig = nSg/math.sqrt(sigma2)
                        sig2DHisto.SetBinContent(lowBin, highBin, sig)
                        #print lowBin, highBin, nJetBin
                        #sigDic[lowBin] = {highBin : {nJetBin : sig}}
                        #sigDic[lowBin,highBin,nJetBin] = sig

            # Make the plot nice
            sig2DHisto.SetTitle(key+" nJetBin "+str(nJetBin)+" Significance")
            sig2DHisto.GetXaxis().SetTitle("Lower Cut [GeV]")
            sig2DHisto.GetYaxis().SetTitle("Higher Cut [GeV]")
            sig2DHisto.Draw("colz")
            c2.SaveAs(key+"_sig_nJetBin"+str(nJetBin)+".png")
            del c2                

if __name__ == '__main__':
    main()
