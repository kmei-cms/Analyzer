import ROOT
import plot
import math
import sys
import argparse
import os

ROOT.gROOT.SetBatch(True)

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--year", dest="year", help="Which period of Run2", default="NULL", type=str) 
parser.add_argument("--fitdir", dest="fitdir", help="Directory containing fits", default="NULL", type=str) 
parser.add_argument("--fittag", dest="fittag", help="Unique tag for fit results", default="NULL", type=str) 

arg = parser.parse_args()

if arg.fitdir == "NULL":
    print "Must provide directory containing fits!"
    print "Exiting..."
    quit()

if arg.fittag == "NULL":
    print "Must provide tag for fit results!"
    print "Exiting..."
    quit()

if not os.path.exists(arg.fitdir):
    print "Fit directory \"%s\" does not exist!"%(arg.fitdir)
    print "Exiting..."
    quit()

year = arg.year
fitdir = arg.fitdir
fittag = arg.fittag

inputfilename = "/uscms_data/d3/jhiltb/susy/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/"
if year == "2016":
    inputfilename += "Keras_2016_v1.1" 
elif year == "2017":
    inputfilename += "Keras_2017_v1.1" 
elif year == "2018pre":
    inputfilename += "Keras_2018pre_v1.0" 
elif year == "2018post":
    inputfilename += "Keras_2018post_v1.0" 
else:
    print "Year \"%s\" is not valid!"%(year)
    print "Exiting..."
    quit()
 
inputfilename += "/njets_for_Aron.root"

print "Opening root file \"%s\""%(inputfilename)
inputfile = ROOT.TFile.Open(inputfilename)

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()

mvas = ["D1","D2","D3","D4"]
base = "TT_h_njets_pt30_1l"
systs = ["", 
         "_JECUp",
         "_JECDown",
         "_JERUp",
         "_JERDown",
         "_btgUp",
         "_btgDown",
         "_lepUp",
         "_lepDown",
         "_pdfUp",
         "_pdfDown",
         "_sclUp",
         "_sclDown",
         "_htUp",
         "_htDown",
         "_puUp",
         "_puDown"
]
if "2017" in year or "2018" in year:
    systs.extend([
         "_isrUp",
         "_isrDown",
         "_fsrUp",
         "_fsrDown",
         ])
colors = [ROOT.kBlack,
          ROOT.kCyan+1,
          ROOT.kCyan+2,
          ROOT.kBlue-7,
          ROOT.kBlue-5,
          ROOT.kRed,
          ROOT.kRed+2,
          ROOT.kGreen-6,
          ROOT.kGreen+2,
          ROOT.kBlue-7,
          ROOT.kBlue-5,
          ROOT.kRed,
          ROOT.kRed+2,
          ROOT.kGreen-6,
          ROOT.kGreen+2,
          ROOT.kBlue-7,
          ROOT.kBlue-5,
          ROOT.kGreen-6,
          ROOT.kGreen+2,
          ROOT.kBlue-7,
          ROOT.kBlue-5]

histos = {}
for mva in mvas: 
    # get all histograms
    #print mva
    print base
    histnames = [mva+"_"+base+syst for syst in systs]
    #print histnames
    Dhistos = [inputfile.Get(histname) for histname in histnames]
    #print Dhistos
    #print "contents: " , Dhistos[-2].GetBinContent(1), " ", Dhistos[-2].GetBinContent(2), " ", Dhistos[-2].GetBinContent(3), " ", Dhistos[-2].GetBinContent(4), " "
    histos[mva] = Dhistos

extra_histos = {}
extra_systs = ["_isrUp",
               "_isrDown",
               "_fsrUp",
               "_fsrDown",
               ]
extra_base = "_h_njets_pt30_1l"
if "2016" in year:
    # deal with the isr and fsr uncertainties
    for mva in mvas:
        extra_histnames = [mva+"_TT"+syst+extra_base for syst in extra_systs]
        print extra_histnames
        extra_Dhistos = [inputfile.Get(histname) for histname in extra_histnames]
        extra_histos[mva] = extra_Dhistos
        histos[mva].extend(extra_Dhistos)
    # now append to other lists
    systs = systs+extra_systs
 
#print systs
#print histos   
# make the summed histograms
sumhistos = []
for i in range(len(systs)):
    #print "making sum for " , systs[i]
    basename = base+systs[i]
    h1 = histos["D1"][i].Clone("Sum"+systs[i])
    h1.Add(histos["D2"][i])
    h1.Add(histos["D3"][i])
    h1.Add(histos["D4"][i])
    sumhistos.append(h1)
    #print "contents: " , h1.GetBinContent(1), " ", h1.GetBinContent(2), " ", h1.GetBinContent(3), " ", h1.GetBinContent(4), " "
    #print "contents: " , histos["D1"][i].GetBinContent(1), " ", histos["D1"][i].GetBinContent(2), " ", histos["D1"][i].GetBinContent(3), " ", histos["D1"][i].GetBinContent(4), " "
    #print "contents: " , histos["D2"][i].GetBinContent(1), " ", histos["D2"][i].GetBinContent(2), " ", histos["D2"][i].GetBinContent(3), " ", histos["D2"][i].GetBinContent(4), " "
    #print "contents: " , histos["D3"][i].GetBinContent(1), " ", histos["D3"][i].GetBinContent(2), " ", histos["D3"][i].GetBinContent(3), " ", histos["D3"][i].GetBinContent(4), " "
    #print "contents: " , histos["D4"][i].GetBinContent(1), " ", histos["D4"][i].GetBinContent(2), " ", histos["D4"][i].GetBinContent(3), " ", histos["D4"][i].GetBinContent(4), " "


sumhistos_norm = []
for sumh in sumhistos:
    newh = sumh.Clone(sumh.GetName()+"_norm")
    newh.Scale(1./newh.Integral())
    sumhistos_norm.append(newh)


outputfile = ROOT.TFile.Open("ttbar_systematics_%s.root"%year, "RECREATE")
outputfile.cd()

for mva in mvas: 
    # make a plot
    print "Making plot for ", mva, " and syst "
    for start in range(1,len(systs),2):        
        # Compare shape in given mva bin to total shape
        ratios = []
        for i in [0, start, start+1]:
            ratio = histos[mva][i].Clone(mva+"_ratio"+systs[i])
            ratio.Scale(1./ratio.Integral())
            ratio.Divide(sumhistos_norm[i])
            ratios.append(ratio)

        # Also make plots of difference wrt nominal
        divratios = []
        for index, i in enumerate([0, start, start+1]):
            divratio = ratios[index].Clone(mva+"_ratio_div"+systs[i])
            divratio.Divide(ratios[0])
            divratios.append(divratio)
            #if i > 4 :
            #    divratio.Write(mva+"_ratio_div"+systs[i])

        # diffratios = []
        # for ratio in ratios:
        #     diffratio = ratio.Clone(mva+"_ratio_diff"+systs[i])
        #     diffratio.Add(ratios[0], -1)
        #     diffratios.append(diffratio)

        # Also make symmetrized versions of the lepton and btag and pdf uncertainties and isr
        if start > 4 and start <= len(systs)-5:
            print "making final histo for systematic: ", systs[start]
            # take the average of Up and 1/Down
            mytemph = divratios[0].Clone(mva+systs[start].replace("Up","temp"))
            #for bin in range(mytemph.GetNbinsX()):
            #    mytemph.SetBinContent(bin+1, 1.)
            #    mytemph.SetBinError(bin+1, 0.)
            #mytemph.Divide(divratios[2])
            #mytemph.Add(divratios[1])
            #mytemph.Scale(0.5)
            myh = ROOT.TH1D(mva+systs[start].replace("Up",""), mva+systs[start].replace("Up",""),
                           divratios[0].GetNbinsX(), 0, divratios[0].GetNbinsX())
            #for bin in range(mytemph.GetNbinsX()):
            #    if mytemph.GetBinContent(bin+1) > 0:
            #        myh.SetBinContent(bin+1, mytemph.GetBinContent(bin+1))
            #    else: # if bin is empty, grab the previous one
            #        myh.SetBinContent(bin+1, mytemph.GetBinContent(bin))
            # Make more robust in case Up/Down go in same direction
            for bin in range(mytemph.GetNbinsX()):
                upval = divratios[1].GetBinContent(bin+1)
                downval = divratios[2].GetBinContent(bin+1)
                if upval > 1 and downval < 1:
                    if downval > 0:
                        downval = 1 / downval
                    else:
                        downval = upval
                elif upval < 1 and downval > 1: 
                    downval = 1 / downval
                average = (upval+downval)/2
                if average > 0:
                    myh.SetBinContent(bin+1, average)
                else:
                    myh.SetBinContent(bin+1, myh.GetBinContent(bin))
            # 2016 PDF unc in bin D3 is subject to very large stat fluctuation, set by hand to value from bin before
            #if "2016" in fitversion and mva is "D3":
            #    if "pdf" in systs[start]:
            #        myh.SetBinContent(8, myh.GetBinContent(7))
            # 2016 pileup unc in bin D4 is subject to very large stat fluctuation, set by hand to value from bin before
            if "2016" in year and mva in ["D4"]:
                if "pu" in systs[start]:
                    myh.SetBinContent(7, myh.GetBinContent(6))
                    myh.SetBinContent(8, myh.GetBinContent(6))
            myh.Write()


        plot.makeplot([histos[mva][0], histos[mva][start], histos[mva][start+1]], [systs[s] for s in [0, start, start+1]], 
                      "N_{j}-7", mva+"_"+systs[start]+"_njets", plotdir="./", linear=False, legendColumns=1, 
                      ymin=0.5, ymax=1e5, ylabel="",
                      colors=[colors[s] for s in [0, start, start+1]], norm=False, drawstyle="hist")


        plot.makeplot(ratios, [systs[s] for s in [0, start, start+1]], "N_{j}-7", mva+"_"+systs[start]+"_comp", plotdir="./", linear=True, legendColumns=1, 
                      ymin=0.5, ymax=1.5, ylabel="",
                      colors=[colors[s] for s in [0, start, start+1]], norm=False, drawstyle="comp")


        plot.makeplot(divratios, [systs[s] for s in [0, start, start+1]], "N_{j}-7", mva+"_"+systs[start]+"_comp_div", plotdir="./", linear=True, legendColumns=1, 
                      ymin=0.9, ymax=1.1, ylabel="", dropzeroes=False,
                      colors=[colors[s] for s in [0, start, start+1]], norm=False, drawstyle="hist")

        # plot.makeplot(diffratios, [systs[s] for s in [0, start, start+1]], "N_{j}-7", mva+"_"+systs[start]+"_comp_diff", plotdir="./", linear=True, legendColumns=1, 
        #               ymin=-0.1, ymax=0.1, dropzeroes=False, ylabel="",
        #               colors=[colors[s] for s in [0, start, start+1]], norm=False, drawstyle="hist")

        plot.makeplot([sumhistos_norm[start].Clone(), histos[mva][0].Clone(), histos[mva][start].Clone(), histos[mva][start+1].Clone()], 
                      ["Nom "+systs[start]]+[systs[s] for s in [0, start, start+1]], 
                      "N_{j}-7", mva+"_"+systs[start]+"_njets_norm", plotdir="./", linear=False, legendColumns=1, 
                      ymin=1e-5, ymax=1, ylabel="",
                      colors=[ROOT.kGray+1]+[colors[s] for s in [0, start, start+1]], norm=True, drawstyle="hist")


#for 2016, also copy over isr from the 2017 file -- Should be replaced by the dedicated systematic samples
#f_2017 = ROOT.TFile.Open("FitInputs/ttbar_systematics_2017_freezing.root")
#h_isr_2017 = [f_2017.Get("%s_isr"%mva) for mva in mvas]
#outputfile.cd()
#for h in h_isr_2017:
#    h.Write()

# Add the QCD CR systematic and the new HT systematics
#if "2016" in fitversion:
#    #f_CR_2016 = ROOT.TFile.Open("/uscms/homes/k/kmei91/public/forNadja/FinalSystematicFiles/qcdCR_shape_systematic_2016.root")
#    f_CR_2016 = ROOT.TFile.Open("/uscms/homes/k/kmei91/public/forNadja/RatioOfRatiosSystematic_12BinInclusive/2016_qcdCR_qcdCRErr_Systematics.root")
#    h_CR_2016 = [f_CR_2016.Get("%s_qcdCR"%mva) for mva in mvas]
#    herr_CR_2016 = [f_CR_2016.Get("%s_qcdCRErr"%mva) for mva in mvas]
#    outputfile.cd()
#    for h in h_CR_2016:
#        newh = ROOT.TH1D(h.GetName(), h.GetTitle(), 6, 0, 6)
#        for bin in range(h.GetNbinsX()):
#            newh.SetBinContent(bin+1, h.GetBinContent(bin+1) if h.GetBinContent(bin+1) > 0 else h.GetBinContent(bin))
#        newh.Write()
#    for h in herr_CR_2016:
#        newh = ROOT.TH1D(h.GetName(), h.GetTitle(), 6, 0, 6)
#        for bin in range(h.GetNbinsX()):
#            newh.SetBinContent(bin+1, h.GetBinContent(bin+1) if h.GetBinContent(bin+1) > 0 else h.GetBinContent(bin))
#        newh.Write()
#
#    f_HT_2016 = ROOT.TFile.Open("/uscms/homes/k/kmei91/public/forNadja/RatioOfRatiosSystematic_12BinInclusive/2016_htRatioSyst.root")
#    h_HTtail_2016 = [f_HT_2016.Get("%s_httail"%mva) for mva in mvas]
#    h_HTnjet_2016 = [f_HT_2016.Get("%s_htnjet"%mva) for mva in mvas]
#    outputfile.cd()
#    for i, h in enumerate(h_HTtail_2016):
#        newh = ROOT.TH1D(h.GetName(), h.GetTitle(), 6, 0, 6)
#        for bin in range(h.GetNbinsX()):
#            newh.SetBinContent(bin+1, h.GetBinContent(bin+1) if h.GetBinContent(bin+1) > 0 else h.GetBinContent(bin))
#        newh.Write()
#    for h in h_HTnjet_2016:
#        newh = ROOT.TH1D(h.GetName(), h.GetTitle(), 6, 0, 6)
#        for bin in range(h.GetNbinsX()):
#            newh.SetBinContent(bin+1, h.GetBinContent(bin+1) if h.GetBinContent(bin+1) > 0 else h.GetBinContent(bin))
#        newh.Write()
#
#        
#if "2017" in fitversion:
#    f_CR_2017 = ROOT.TFile.Open("/uscms/homes/k/kmei91/public/forNadja/RatioOfRatiosSystematic_12BinInclusive/2017_qcdCR_qcdCRErr_Systematics.root")
#    #f_CR_2017 = ROOT.TFile.Open("/uscms_data/d1/owenl/stealth-rpv/mva-shape-syst/dratio-rw5-fab5.root")
#    #print ["h_njets_dratio_mva%s_ab_rw5_fab"%mva for mva in ["1","2","3","4"]]
#    #h_CR_2017 = [f_CR_2017.Get("h_njets_dratio_mva%s_ab_rw5_fab5"%mva) for mva in ["1","2","3","4"]]
#    #f_CR_2017 = ROOT.TFile.Open("/uscms/homes/k/kmei91/public/forNadja/FinalSystematicFiles/qcdCR_shape_systematic_2017.root")
#    h_CR_2017 = [f_CR_2017.Get("%s_qcdCR"%mva) for mva in mvas]
#    herr_CR_2017 = [f_CR_2017.Get("%s_qcdCRErr"%mva) for mva in mvas]
#    outputfile.cd()
#    for i,h in enumerate(h_CR_2017):
#        #print h
#        newh = ROOT.TH1D("D%s_qcdCR"%(i+1), "D%s_qcdCR"%(i+1), 6, 0, 6)
#        for bin in range(h.GetNbinsX()):
#            newh.SetBinContent(bin+1, h.GetBinContent(bin+1) if h.GetBinContent(bin+1) > 0 else h.GetBinContent(bin))
#            # if i == 2 and bin >= 5:
#            #     newh.SetBinContent(bin+1, 1.)
#            # elif i == 3 and bin >=4:
#            #     newh.SetBinContent(bin+1, 1.)
#            # else:
#            #     newh.SetBinContent(bin+1, h.GetBinContent(bin+1) if h.GetBinContent(bin+1) > 0 else h.GetBinContent(bin))
#        newh.Write()
#    for i,h in enumerate(herr_CR_2017):
#        newh = ROOT.TH1D("D%s_qcdCRErr"%(i+1), "D%s_qcdCRErr"%(i+1), 6, 0, 6)
#        for bin in range(h.GetNbinsX()):
#            newh.SetBinContent(bin+1, h.GetBinContent(bin+1) if h.GetBinContent(bin+1) > 0 else h.GetBinContent(bin))
#        newh.Write()
#
#    f_HT_2017 = ROOT.TFile.Open("/uscms/homes/k/kmei91/public/forNadja/RatioOfRatiosSystematic_12BinInclusive/2017_htRatioSyst.root")
#    h_HTtail_2017 = [f_HT_2017.Get("%s_httail"%mva) for mva in mvas]
#    h_HTnjet_2017 = [f_HT_2017.Get("%s_htnjet"%mva) for mva in mvas]
#    outputfile.cd()
#    for h in h_HTtail_2017:
#        newh = ROOT.TH1D(h.GetName(), h.GetTitle(), 6, 0, 6)
#        for bin in range(h.GetNbinsX()):
#            newh.SetBinContent(bin+1, h.GetBinContent(bin+1) if h.GetBinContent(bin+1) > 0 else h.GetBinContent(bin))
#        newh.Write()
#    for h in h_HTnjet_2017:
#        newh = ROOT.TH1D(h.GetName(), h.GetTitle(), 6, 0, 6)
#        for bin in range(h.GetNbinsX()):
#            newh.SetBinContent(bin+1, h.GetBinContent(bin+1) if h.GetBinContent(bin+1) > 0 else h.GetBinContent(bin))
#        newh.Write()


# Add the systematic from Owen
# if "2016" in fitversion:
#     f_owen = ROOT.TFile.Open("/uscms_data/d1/owenl/stealth-rpv/signal-syst-shape-2016_RPV_2t6j_mStop-450.root")
#     h_owen = [f_owen.Get("h_njets_double_ratio_%s_rebin"%mva.replace("D","d")).Clone("%s_ttHTshape"%mva) for mva in mvas]
#     outputfile.cd()
#     for h in h_owen:
#         h.Write()

#*(1) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of nominal samples*
#Results are here:
fname_nom_shared = "%s/%s_nom_shared/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_nom_shared = ROOT.TFile.Open(fname_nom_shared)
histos_nom_shared = {}
for mva in mvas:
    h = f_nom_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_nom_shared[mva] = h

#*(2) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of nominal samples*
#Results are here:
fname_nom_sep = "%s/%s_nom_sep/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_nom_sep = ROOT.TFile.Open(fname_nom_sep)
histos_nom_sep = {}
for mva in mvas:
    h = f_nom_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_nom_sep[mva] = h


#*(3) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of JEC UP samples*
#Results are here:
fname_JECUp_shared = "%s/%s_JECUp_shared/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_JECUp_shared = ROOT.TFile.Open(fname_JECUp_shared)
histos_JECUp_shared = {}
for mva in mvas:
    h = f_JECUp_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_JECUp_shared[mva] = h

#*(4) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of JEC UP samples*
#Results are here:
fname_JECUp_sep = "%s/%s_JECUp_sep/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_JECUp_sep = ROOT.TFile.Open(fname_JECUp_sep)
histos_JECUp_sep = {}
for mva in mvas:
    h = f_JECUp_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_JECUp_sep[mva] = h

#*(5) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of JEC Down samples*
#Results are here:
fname_JECDown_shared = "%s/%s_JECDown_shared/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_JECDown_shared = ROOT.TFile.Open(fname_JECDown_shared)
histos_JECDown_shared = {}
for mva in mvas:
    h = f_JECDown_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_JECDown_shared[mva] = h

#*(6) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of JEC Down samples*
#Results are here:
fname_JECDown_sep = "%s/%s_JECDown_sep/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_JECDown_sep = ROOT.TFile.Open(fname_JECDown_sep)
histos_JECDown_sep = {}
for mva in mvas:
    h = f_JECDown_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_JECDown_sep[mva] = h




#*(7) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of JER UP samples*
#Results are here:
fname_JERUp_shared = "%s/%s_JERUp_shared/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_JERUp_shared = ROOT.TFile.Open(fname_JERUp_shared)
histos_JERUp_shared = {}
for mva in mvas:
    h = f_JERUp_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_JERUp_shared[mva] = h

#*(8) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of JER UP samples*
#Results are here:
fname_JERUp_sep = "%s/%s_JERUp_sep/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_JERUp_sep = ROOT.TFile.Open(fname_JERUp_sep)
histos_JERUp_sep = {}
for mva in mvas:
    h = f_JERUp_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_JERUp_sep[mva] = h

#*(9) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of JER Down samples*
#Results are here:
fname_JERDown_shared = "%s/%s_JERDown_shared/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_JERDown_shared = ROOT.TFile.Open(fname_JERDown_shared)
histos_JERDown_shared = {}
for mva in mvas:
    h = f_JERDown_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_JERDown_shared[mva] = h

#*(10) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of JER Down samples*
#Results are here:
fname_JERDown_sep = "%s/%s_JERDown_sep/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
f_JERDown_sep = ROOT.TFile.Open(fname_JERDown_sep)
histos_JERDown_sep = {}
for mva in mvas:
    h = f_JERDown_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_JERDown_sep[mva] = h



histos_FSRUp_shared = {}
histos_FSRUp_sep = {}
histos_FSRDown_shared = {}
histos_FSRDown_sep = {}
if "201" in year:
    #*(11) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of FSR UP samples*
    #Results are here:
    fname_FSRUp_shared = "%s/%s_FSRUp_shared/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
    f_FSRUp_shared = ROOT.TFile.Open(fname_FSRUp_shared)
    for mva in mvas:
        h = f_FSRUp_shared.Get("shapes_fit_b/%s/TT" % mva)
        histos_FSRUp_shared[mva] = h
    
    #*(12) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of FSR UP samples*
    #Results are here:
    fname_FSRUp_sep = "%s/%s_FSRUp_sep/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
    f_FSRUp_sep = ROOT.TFile.Open(fname_FSRUp_sep)
    for mva in mvas:
        h = f_FSRUp_sep.Get("shapes_fit_b/%s/TT" % mva)
        histos_FSRUp_sep[mva] = h
    
    #*(13) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of FSR Down samples*
    #Results are here:
    fname_FSRDown_shared = "%s/%s_FSRDown_shared/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
    f_FSRDown_shared = ROOT.TFile.Open(fname_FSRDown_shared)
    for mva in mvas:
        h = f_FSRDown_shared.Get("shapes_fit_b/%s/TT" % mva)
        histos_FSRDown_shared[mva] = h
    
    #*(14) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of FSR Down samples*
    #Results are here:
    fname_FSRDown_sep = "%s/%s_FSRDown_sep/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
    f_FSRDown_sep = ROOT.TFile.Open(fname_FSRDown_sep)
    for mva in mvas:
        h = f_FSRDown_sep.Get("shapes_fit_b/%s/TT" % mva)
        histos_FSRDown_sep[mva] = h

histos_ISRUp_shared = {}
histos_ISRUp_sep = {}
histos_ISRDown_shared = {}
histos_ISRDown_sep = {}
if "201" in year:
    #*(15) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of ISR UP samples*
    #Results are here:
    fname_ISRUp_shared = "%s/%s_ISRUp_shared/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
    f_ISRUp_shared = ROOT.TFile.Open(fname_ISRUp_shared)
    for mva in mvas:
        h = f_ISRUp_shared.Get("shapes_fit_b/%s/TT" % mva)
        histos_ISRUp_shared[mva] = h
    
    #*(16) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of ISR UP samples*
    #Results are here:
    fname_ISRUp_sep = "%s/%s_ISRUp_sep/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
    f_ISRUp_sep = ROOT.TFile.Open(fname_ISRUp_sep)
    for mva in mvas:
        h = f_ISRUp_sep.Get("shapes_fit_b/%s/TT" % mva)
        histos_ISRUp_sep[mva] = h
    
    #*(17) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of ISR Down samples*
    #Results are here:
    fname_ISRDown_shared = "%s/%s_ISRDown_shared/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
    f_ISRDown_shared = ROOT.TFile.Open(fname_ISRDown_shared)
    for mva in mvas:
        h = f_ISRDown_shared.Get("shapes_fit_b/%s/TT" % mva)
        histos_ISRDown_shared[mva] = h
    
    #*(18) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of ISR Down samples*
    #Results are here:
    fname_ISRDown_sep = "%s/%s_ISRDown_sep/fitDiagnostics%sRPV550.root"%(fitdir,fittag,year)
    f_ISRDown_sep = ROOT.TFile.Open(fname_ISRDown_sep)
    for mva in mvas:
        h = f_ISRDown_sep.Get("shapes_fit_b/%s/TT" % mva)
        histos_ISRDown_sep[mva] = h
    


#For each of the above files you can
#- Navigate to the TDirectoryFile called shapes_fit_s
#- Navigate to each TDirectoryFile called sigD1, sigD2, sigD3, sigD4 for the corresponding MVA bin
#- Extract information from bins 1 though 8 of the histogram called bkg_tt

# fname_rebin = "/uscms_data/d2/soha/stealth/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/Keras_V1.2.4/njets_rebin_for_Aron.root"
# inputfile_rebin = ROOT.TFile.Open(fname_rebin)
inputfile_rebin=inputfile # normal file now has correct bin edges
histos_rebin = histos
# for mva in mvas: 
#     # get all histograms
#     histnames = [mva+"_"+base+syst for syst in systs]
#     Dhistos = [inputfile_rebin.Get(histname) for histname in histnames]
#     histos_rebin[mva] = Dhistos

colors2 = [ROOT.kBlack,
           ROOT.kCyan+1,
           #ROOT.kCyan+3,
           ROOT.kRed,
           ]
colors3 = [ROOT.kGray+1,
           ROOT.kCyan+1,
           ROOT.kCyan+3,
           ROOT.kRed,
           ]

for mva in mvas: 
    # first make plot of the actual distributions
    # for nominal case
    plot.makeplot([histos_nom_sep[mva].Clone(), histos_nom_shared[mva].Clone(), histos_rebin[mva][0].Clone()], ["Separate", "Shared", "MC"], 
                  "N_{j}-7", mva+"_nom_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                  ymin=0.01, ymax=1e5, 
                  ylabel="",
                  colors=colors2, norm=False, drawstyle="lastP")

    # for syst JEC Up case
    plot.makeplot([histos_nom_shared[mva].Clone(), histos_JECUp_sep[mva].Clone(), histos_JECUp_shared[mva].Clone(), histos_rebin[mva][1].Clone()], 
                  ["Nominal","Separate", "Shared", "MC"], 
                  "N_{j}-7", mva+"_JECUp_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                  ymin=0.01, ymax=1e5, 
                  ylabel="",
                  colors=colors3, norm=False, drawstyle="lastP")
    
    # for JEC Down case
    plot.makeplot([histos_nom_shared[mva].Clone(), histos_JECDown_sep[mva].Clone(), histos_JECDown_shared[mva].Clone(), histos_rebin[mva][2].Clone()], 
                  ["Nominal","Separate", "Shared", "MC"], 
                  "N_{j}-7", mva+"_JECDown_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                  ymin=0.01, ymax=1e5, 
                  ylabel="",
                  colors=colors3, norm=False, drawstyle="lastP")
    
    # for JER Up case
    plot.makeplot([histos_nom_shared[mva].Clone(), histos_JERUp_sep[mva].Clone(), histos_JERUp_shared[mva].Clone(), histos_rebin[mva][3].Clone()], 
                  ["Nominal","Separate", "Shared", "MC"], 
                  "N_{j}-7", mva+"_JERUp_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                  ymin=0.01, ymax=1e5, 
                  ylabel="",
                  colors=colors3, norm=False, drawstyle="lastP")
    
    # for JER Down case
    plot.makeplot([histos_nom_shared[mva].Clone(), histos_JERDown_sep[mva].Clone(), histos_JERDown_shared[mva].Clone(), histos_rebin[mva][4].Clone()], 
                  ["Nominal","Separate", "Shared", "MC"], 
                  "N_{j}-7", mva+"_JERDown_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                  ymin=0.01, ymax=1e5, 
                  ylabel="",
                  colors=colors3, norm=False, drawstyle="lastP")
   
    if "201" in year:
        # for FSR Up case
        index_FSRUp = systs.index("_fsrUp")
        plot.makeplot([histos_nom_shared[mva].Clone(), histos_FSRUp_sep[mva].Clone(), histos_FSRUp_shared[mva].Clone(), histos_rebin[mva][index_FSRUp].Clone()], 
                      ["Nominal","Separate", "Shared", "MC"], 
                      "N_{j}-7", mva+"_FSRUp_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                      ymin=0.01, ymax=1e5, 
                      ylabel="",
                      colors=colors3, norm=False, drawstyle="lastP")
        
        # for FSR Down case
        index_FSRDown = systs.index("_fsrDown")
        plot.makeplot([histos_nom_shared[mva].Clone(), histos_FSRDown_sep[mva].Clone(), histos_FSRDown_shared[mva].Clone(), histos_rebin[mva][index_FSRDown].Clone()], 
                      ["Nominal","Separate", "Shared", "MC"], 
                      "N_{j}-7", mva+"_FSRDown_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                      ymin=0.01, ymax=1e5, 
                      ylabel="",
                      colors=colors3, norm=False, drawstyle="lastP")
        
    if "201" in year: 
        # for ISR Up case
        index_ISRUp = systs.index("_isrUp")
        plot.makeplot([histos_nom_shared[mva].Clone(), histos_ISRUp_sep[mva].Clone(), histos_ISRUp_shared[mva].Clone(), histos_rebin[mva][index_ISRUp].Clone()], 
                      ["Nominal","Separate", "Shared", "MC"], 
                      "N_{j}-7", mva+"_ISRUp_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                      ymin=0.01, ymax=1e5, 
                      ylabel="",
                      colors=colors3, norm=False, drawstyle="lastP")
        
        # for ISR Down case
        index_ISRDown = systs.index("_isrDown")
        plot.makeplot([histos_nom_shared[mva].Clone(), histos_ISRDown_sep[mva].Clone(), histos_ISRDown_shared[mva].Clone(), histos_rebin[mva][index_ISRDown].Clone()], 
                      ["Nominal","Separate", "Shared", "MC"], 
                      "N_{j}-7", mva+"_ISRDown_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                      ymin=0.01, ymax=1e5, 
                      ylabel="",
                      colors=colors3, norm=False, drawstyle="lastP")
        

    # make a plot
    print "Making plot for ", mva, " and syst JEC with fits"
    # Compare shape in given mva bin to total shape when fitting all mva bins together
    # first for the nominal case
    ratio_nom = histos_nom_sep[mva].Clone(mva+"_ratio_nom_fit")
    #ratio_nom.Scale(1./ratio_nom.Integral())
    ratio_nom.Divide(histos_nom_shared[mva])

    # then for JEC Up case
    ratio_JECUp = histos_JECUp_sep[mva].Clone(mva+"_ratio_JECUp_fit")
    #ratio_sys.Scale(1./ratio_sys.Integral())
    ratio_JECUp.Divide(histos_JECUp_shared[mva])
    # For the JEC Down case
    ratio_JECDown = histos_JECDown_sep[mva].Clone(mva+"_ratio_JECDown_fit")
    ratio_JECDown.Divide(histos_JECDown_shared[mva])
    
    # For the JER Up case
    ratio_JERUp = histos_JERUp_sep[mva].Clone(mva+"_ratio_JERUp_fit")
    ratio_JERUp.Divide(histos_JERUp_shared[mva])
    # For the JER Down case
    ratio_JERDown = histos_JERDown_sep[mva].Clone(mva+"_ratio_JERDown_fit")
    ratio_JERDown.Divide(histos_JERDown_shared[mva])
    
    ratio_FSRUp = None
    ratio_FSRDown = None
    if "201" in year:
        # For the FSR Up case
        ratio_FSRUp = histos_FSRUp_sep[mva].Clone(mva+"_ratio_FSRUp_fit")
        ratio_FSRUp.Divide(histos_FSRUp_shared[mva])
        # For the FSR Down case
        ratio_FSRDown = histos_FSRDown_sep[mva].Clone(mva+"_ratio_FSRDown_fit")
        ratio_FSRDown.Divide(histos_FSRDown_shared[mva])

    ratio_ISRUp = None
    ratio_ISRDown = None
    if "201" in year:
        # For the ISR Up case
        ratio_ISRUp = histos_ISRUp_sep[mva].Clone(mva+"_ratio_ISRUp_fit")
        ratio_ISRUp.Divide(histos_ISRUp_shared[mva])
        # For the ISR Down case
        ratio_ISRDown = histos_ISRDown_sep[mva].Clone(mva+"_ratio_ISRDown_fit")
        ratio_ISRDown.Divide(histos_ISRDown_shared[mva])


    # Also make plots of difference wrt nominal
    divratio_JECUp = ratio_JECUp.Clone(mva+"_ratio_JECUp_div_fit")
    divratio_JECUp.Divide(ratio_nom)
    divratio_JECDown = ratio_JECDown.Clone(mva+"_ratio_JECDown_div_fit")
    divratio_JECDown.Divide(ratio_nom)
    divratio_JERUp = ratio_JERUp.Clone(mva+"_ratio_JERUp_div_fit")
    divratio_JERUp.Divide(ratio_nom)
    divratio_JERDown = ratio_JERDown.Clone(mva+"_ratio_JERDown_div_fit")
    divratio_JERDown.Divide(ratio_nom)
    divratio_FSRUp = None
    divratio_FSRDown = None
    if "201" in year:
        divratio_FSRUp = ratio_FSRUp.Clone(mva+"_ratio_FSRUp_div_fit")
        divratio_FSRUp.Divide(ratio_nom)
        divratio_FSRDown = ratio_FSRDown.Clone(mva+"_ratio_FSRDown_div_fit")
        divratio_FSRDown.Divide(ratio_nom)
    divratio_ISRUp = None
    divratio_ISRDown = None
    if "201" in year:
        divratio_ISRUp = ratio_ISRUp.Clone(mva+"_ratio_ISRUp_div_fit")
        divratio_ISRUp.Divide(ratio_nom)
        divratio_ISRDown = ratio_ISRDown.Clone(mva+"_ratio_ISRDown_div_fit")
        divratio_ISRDown.Divide(ratio_nom)

    outputfile.cd()
    # Rather than put JECUp and JERUp, construct maximum of the Up and Down variations, but preserve the "sign"
    mva_JEC = None
    if "2016" in year:
        # Use JEC Down
        mva_JEC = divratio_JECDown.Clone(mva+"_JEC")
        for i in range(divratio_JECDown.GetNbinsX()):
            #print i
            up_value = divratio_JECUp.GetBinContent(i+1)
            #print "up: " , up_value
            if up_value>1:
                up_value = 1/up_value
            down_value = divratio_JECDown.GetBinContent(i+1)
            #print "down: " , down_value
            if down_value>1:
                down_value = 1/down_value
            if up_value < down_value:
                #print "changing value"
                if divratio_JECUp.GetBinContent(i+1) < 1:
                    mva_JEC.SetBinContent(i+1,up_value)
                    #print "to ", up_value
                else:
                    mva_JEC.SetBinContent(i+1,1/up_value)
                    #print "to " , 1/up_value
    else:
        mva_JEC = divratio_JECUp.Clone(mva+"_JEC")
        for i in range(divratio_JECUp.GetNbinsX()):
            #print i
            up_value = divratio_JECUp.GetBinContent(i+1)
            #print "up: " , up_value
            if up_value<1:
                up_value = 1/up_value
            down_value = divratio_JECDown.GetBinContent(i+1)
            #print "down: " , down_value
            if down_value<1:
                down_value = 1/down_value
            if down_value > up_value:
                #print "changing value"
                if divratio_JECUp.GetBinContent(i+1) > 1:
                    mva_JEC.SetBinContent(i+1,down_value)
                    #print "to ", down_value
                else:
                    mva_JEC.SetBinContent(i+1,1/down_value)
                    #print "to " , 1/down_value

    mva_JER = divratio_JERUp.Clone(mva+"_JER")  # should be updated to use both if fit succeeds
    for i in range(divratio_JERUp.GetNbinsX()):
        up_value = divratio_JERUp.GetBinContent(i+1)
        if up_value<1:
            up_value = 1/up_value
        down_value = divratio_JERDown.GetBinContent(i+1)
        if down_value<1:
            down_value = 1/down_value
        if down_value > up_value:
            if divratio_JERUp.GetBinContent(i+1) > 1:
                mva_JER.SetBinContent(i+1,down_value)
            else:
                mva_JER.SetBinContent(i+1,1/down_value)

    mva_FSR = None
    if "201" in year:
        mva_FSR = divratio_FSRUp.Clone(mva+"_FSR")
        for i in range(divratio_FSRUp.GetNbinsX()):
            print "bin ", i
            up_value = divratio_FSRUp.GetBinContent(i+1)
            print "Up: ", up_value
            if up_value<1:
                up_value = 1/up_value
                print "inverting up, ", up_value
            down_value = divratio_FSRDown.GetBinContent(i+1)
            print "Down: " , down_value
            if down_value<1:
                down_value = 1/down_value
                print "inverting down, ", down_value
            if down_value > up_value:
                if divratio_FSRUp.GetBinContent(i+1) > 1:
                    mva_FSR.SetBinContent(i+1,down_value)
                    print "writing down value", down_value
                else:
                    mva_FSR.SetBinContent(i+1,1/down_value)
                    print "writing 1/down", 1/down_value
    
    mva_ISR = None
    if "201" in year:
        mva_ISR = divratio_ISRUp.Clone(mva+"_ISR")
        for i in range(divratio_ISRUp.GetNbinsX()):
            print "bin ", i
            up_value = divratio_ISRUp.GetBinContent(i+1)
            print "Up: ", up_value
            if up_value<1:
                up_value = 1/up_value
                print "inverting up, ", up_value
            down_value = divratio_ISRDown.GetBinContent(i+1)
            print "Down: " , down_value
            if down_value<1:
                down_value = 1/down_value
                print "inverting down, ", down_value
            if down_value > up_value:
                if divratio_ISRUp.GetBinContent(i+1) > 1:
                    mva_ISR.SetBinContent(i+1,down_value)
                    print "writing down value", down_value
                else:
                    mva_ISR.SetBinContent(i+1,1/down_value)
                    print "writing 1/down", 1/down_value
    if "2016" in year:
        for bin in range(mva_ISR.GetNbinsX()):
            old_value = mva_ISR.GetBinContent(bin+1)
            diff_value = abs(old_value-1)/math.sqrt(2)
            new_value = 1 + diff_value if (old_value > 1) else 1 - diff_value
            mva_ISR.SetBinContent(bin+1, new_value)
        for bin in range(mva_FSR.GetNbinsX()):
            old_value = mva_FSR.GetBinContent(bin+1)
            diff_value = abs(old_value-1)/math.sqrt(2)
            new_value = 1 + diff_value if (old_value > 1) else 1 - diff_value
            mva_FSR.SetBinContent(bin+1, new_value)

    divratio_JECUp.Write(mva+"_JECUp")
    divratio_JECDown.Write(mva+"_JECDown")
    divratio_JERUp.Write(mva+"_JERUp")
    divratio_JERDown.Write(mva+"_JERDown")
    #mva_JEC.Write(mva+"_JEC")
    #mva_JER.Write(mva+"_JER")
    if "201" in year:
        mva_FSR.Write(mva+"_FSR")
    if "201" in year:
        mva_ISR.Write(mva+"_ISR")


    # Also write the shape differences for the no-systematics case to the file
    ratio_nom.Write(mva+"_nom")

    #diffratio = ratio_sys.Clone(mva+"_ratio_diff_fit")
    #diffratio.Add(ratio_nom, -1)
    

    plot.makeplot([ratio_nom, ratio_JECUp, ratio_JECDown], ["Nominal", "JEC Up", "JEC Down"], "N_{j}-7", mva+"_fit"+"_comp_JEC", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.7, ymax=1.3, ylabel="",
    colors=colors, norm=False, drawstyle="hist")

    plot.makeplot([ratio_nom, ratio_JERUp, ratio_JERDown], ["Nominal", "JER Up", "JER Down"], "N_{j}-7", mva+"_fit"+"_comp_JER", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.7, ymax=1.3, ylabel="",
    colors=colors, norm=False, drawstyle="hist")
    #plot.makeplot([ratio_nom, ratio_JERDown], ["Nominal", "JER Down"], "N_{j}-7", mva+"_fit"+"_comp_JER", plotdir="./", linear=True, legendColumns=1, 
    #ymin=0.7, ymax=1.3, ylabel="",
    #colors=colors, norm=False, drawstyle="hist")

    if "201" in year:
        plot.makeplot([ratio_nom, ratio_FSRUp, ratio_FSRDown], ["Nominal", "FSR Up", "FSR Down"], "N_{j}-7", mva+"_fit"+"_comp_FSR", plotdir="./", linear=True, legendColumns=1, 
                      ymin=0.7, ymax=1.3, ylabel="",
                      colors=colors, norm=False, drawstyle="hist")
    if "201" in year:
        plot.makeplot([ratio_nom, ratio_ISRUp, ratio_ISRDown], ["Nominal", "ISR Up", "ISR Down"], "N_{j}-7", mva+"_fit"+"_comp_ISR", plotdir="./", linear=True, legendColumns=1, 
                      ymin=0.7, ymax=1.3, ylabel="",
                      colors=colors, norm=False, drawstyle="hist")

    
    plot.makeplot([divratio_JECUp, divratio_JECDown], ["JEC Up","JEC Down"], "N_{j}-7", mva+"_fit_comp_div_JEC", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
    colors=colors, norm=False, drawstyle="hist")

    plot.makeplot([divratio_JERUp, divratio_JERDown], ["JER Up","JER Down"], "N_{j}-7", mva+"_fit_comp_div_JER", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
    colors=colors, norm=False, drawstyle="hist")
    #plot.makeplot([divratio_JERDown], ["JER Down"], "N_{j}-7", mva+"_fit_comp_div_JER", plotdir="./", linear=True, legendColumns=1, 
    #ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
    #colors=colors, norm=False, drawstyle="hist")

    if "201" in year:
        plot.makeplot([divratio_FSRUp, divratio_FSRDown], ["FSR Up","FSR Down"], "N_{j}-7", mva+"_fit_comp_div_FSR", plotdir="./", linear=True, legendColumns=1, 
                      ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
                      colors=colors, norm=False, drawstyle="hist")
    if "201" in year:
        plot.makeplot([divratio_ISRUp, divratio_ISRDown], ["ISR Up","ISR Down"], "N_{j}-7", mva+"_fit_comp_div_ISR", plotdir="./", linear=True, legendColumns=1, 
                      ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
                      colors=colors, norm=False, drawstyle="hist")

    plot.makeplot([ratio_nom], ["Nominal"], "N_{j}-7", mva+"_fit"+"_nom", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.7, ymax=1.3, ylabel="",
    colors=colors, norm=False, drawstyle="hist")


    #plot.makeplot([diffratio], ["JEC Up"], "N_{j}-7", mva+"_fit_comp_diff", plotdir="./", linear=True, legendColumns=1, 
    #ymin=-0.2, ymax=0.2, dropzeroes=False, ylabel="",
    #colors=colors, norm=False, drawstyle="hist")



outputfile.Close()


