import ROOT
import math
import DataSetInfo as info
from ROOT import TFile

def main():
    path = "condor/hadd_2016_MC_4thSlides_25.07.2019/"
    histoName = "h_njets_dR_bjet1_bjet2_0l_HT500_ge2b_ge2t" # 2D histo from analyzer / x axes is deltaR, y axes is njets
    #rebinVal = 20 # 20

    bgData = {
        "TT"     : info.DataSetInfo(basedir=path, fileName="2016_TT.root",  sys=0.2, label="Background"),
        "QCD"    : info.DataSetInfo(basedir=path, fileName="2016_QCD.root", sys=0.2, label="Background"),
        "Others" : info.DataSetInfo(basedir=path, fileName="others.root",   sys=0.2, label="Background"),       
        }

    sgData = {
        "rpv_stop_300"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-300.root",        sys=-1.0, label="Signal"),
        "rpv_stop_350"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-350.root",        sys=-1.0, label="Signal"),
        "rpv_stop_400"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-400.root",        sys=-1.0, label="Signal"),
        "rpv_stop_450"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-450.root",        sys=-1.0, label="Signal"),
        "rpv_stop_500"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-500.root",        sys=-1.0, label="Signal"),
        "rpv_stop_550"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-550.root",        sys=-1.0, label="Signal"),
        "rpv_stop_600"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-600.root",        sys=-1.0, label="Signal"),        
        "rpv_stop_650"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-650.root",        sys=-1.0, label="Signal"),
        "rpv_stop_700"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-700.root",        sys=-1.0, label="Signal"),
        "rpv_stop_750"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-750.root",        sys=-1.0, label="Signal"),
        "rpv_stop_800"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-800.root",        sys=-1.0, label="Signal"),
        "rpv_stop_850"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-850.root",        sys=-1.0, label="Signal"),
        "rpv_stop_900"         : info.DataSetInfo(basedir=path, fileName="2016_RPV_2t6j_mStop-900.root",        sys=-1.0, label="Signal"),    


        "stealth_stop_300_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-300.root", sys=-1.0, label="Signal"),
        "stealth_stop_350_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-350.root", sys=-1.0, label="Signal"),
        "stealth_stop_400_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-400.root", sys=-1.0, label="Signal"),
        "stealth_stop_450_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-450.root", sys=-1.0, label="Signal"),
        "stealth_stop_500_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-500.root", sys=-1.0, label="Signal"),
        "stealth_stop_550_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-550.root", sys=-1.0, label="Signal"),
        "stealth_stop_600_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-600.root", sys=-1.0, label="Signal"),
        "stealth_stop_650_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-650.root", sys=-1.0, label="Signal"),
        "stealth_stop_700_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-700.root", sys=-1.0, label="Signal"),
        "stealth_stop_750_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-750.root", sys=-1.0, label="Signal"),
        "stealth_stop_800_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-800.root", sys=-1.0, label="Signal"),
        "stealth_stop_850_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-850.root", sys=-1.0, label="Signal"),
        "stealth_stop_900_SYY" : info.DataSetInfo(basedir=path, fileName="2016_StealthSYY_2t6j_mStop-900.root", sys=-1.0, label="Signal"),
        

    }

    # Loop over all background and get their histograms
    bgHistos = {}
    for key in bgData:
        h = bgData[key].getHisto(histoName)
        h.RebinX(rebinVal)
        bgHistos[key] = {"hist" : h}

    # Loop over all signal models and calculate sig. for each one 
    for key in sgData:
        hSg = sgData[key].getHisto(histoName)
        hSg.RebinX(rebinVal)
        nBins = hSg.GetNbinsX()
        low = hSg.GetXaxis().GetXmin()
        high = hSg.GetXaxis().GetXmax()
        binToGeV = (high - low)/nBins
        sigDic = {}
            
        # Find all possible MLB cut values
        c1 = ROOT.TCanvas( "c", "c", 0, 0, 800, 800)
        ROOT.gStyle.SetStatY(0.5)
        ROOT.gStyle.SetStatX(0.85)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTicks(1,1)
        sig1DHisto = ROOT.TH1D(key+"hist", key+"hist", nBins, low, high)         

        for lowBin in range(1, nBins+1):       
            sigTot2=0.0
            for nJetBin in range(0, hSg.GetNbinsY()+1):            
                highBin=hSg.GetNbinsX()
                #nSg = hSg.Integral(nJetBin, nJetBin, lowBin, highBin)
                nSg = hSg.Integral(lowBin, highBin, nJetBin, nJetBin) # x axes is deltaR, y axes is njets
 
                # Loop over all background histos and calculate sigma
                sigma2 = 0.0
                for k in bgHistos:
                    #n = bgHistos[k]["hist"].Integral(nJetBin, nJetBin, lowBin, highBin)
                    n = bgHistos[k]["hist"].Integral(lowBin, highBin, nJetBin, nJetBin)
                    sys = ((bgData[k].sys)*n)**2
                    sigma2 += n + sys
            
                # Calculate and fill sig.
                sig = 0.0
                if(sigma2 > 0.0): 
                    sig = nSg/math.sqrt(sigma2)
                sigTot2 += sig**2

            sig1DHisto.SetBinContent(lowBin, math.sqrt(sigTot2))

        # Make the plot nice
        sig1DHisto.SetTitle(key+" Significance")
        sig1DHisto.GetXaxis().SetTitle("#DeltaR_{bj1-bj2}")
        sig1DHisto.GetYaxis().SetTitle("Significance")
        #binmax = sig1DHisto.GetMaximumBin() # for the best cut bin
        #sig1DHisto.GetBinCenter(binmax)
        sig1DHisto.Draw("")
        c1.SaveAs(key+"_sig.pdf")
        del c1                

if __name__ == '__main__':
    main()
