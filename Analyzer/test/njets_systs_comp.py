import ROOT
import plot

version = "Keras_V1.2.5_v2"
#version = "Keras_V3.0.1_v2"

fitversion = "FS"
#fitversion = "FS2017"

#inputfilename = "~cmadrid/nobackup/ana/SUSY/Stealth/AnaNTuples/CMSSW_9_3_3/src/Analyzer/Analyzer/test/FitInput/Keras_V1.2.4/njets_for_Aron.root"
inputfilename = "~cmadrid/nobackup/ana/SUSY/Stealth/AnaNTuples/CMSSW_9_3_3/src/Analyzer/Analyzer/test/FitInput/%s/njets_for_Aron.root"%version
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
         #"_isrUp",
         #"_isrDown",
         #"_fsrUp",
         #"_fsrDown",
         #"_sclUp",
         #"_sclDown",
         "_htUp",
         "_htDown"
         ]
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
          ROOT.kBlue-5]

histos = {}
for mva in mvas: 
    # get all histograms
    #print mva
    histnames = [mva+"_"+base+syst for syst in systs]
    #print histnames
    Dhistos = [inputfile.Get(histname) for histname in histnames]
    #print Dhistos
    #print "contents: " , Dhistos[-2].GetBinContent(1), " ", Dhistos[-2].GetBinContent(2), " ", Dhistos[-2].GetBinContent(3), " ", Dhistos[-2].GetBinContent(4), " "
    histos[mva] = Dhistos

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


outputfile = ROOT.TFile.Open("ttbar_systematics.root", "RECREATE")
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
        if start > 4: # and start <= 12:
            # take the average of Up and 1/Down
            mytemph = divratios[0].Clone(mva+systs[start].replace("Up","temp"))
            for bin in range(mytemph.GetNbinsX()):
                mytemph.SetBinContent(bin+1, 1.)
                mytemph.SetBinError(bin+1, 0.)
            mytemph.Divide(divratios[2])
            mytemph.Add(divratios[1])
            mytemph.Scale(0.5)
            myh = ROOT.TH1D(mva+systs[start].replace("Up",""), mva+systs[start].replace("Up",""),
                           divratios[0].GetNbinsX(), 0, divratios[0].GetNbinsX())
            for bin in range(mytemph.GetNbinsX()):
                if mytemph.GetBinContent(bin+1) > 0:
                    myh.SetBinContent(bin+1, mytemph.GetBinContent(bin+1))
                else: # if bin is empty, grab the previous one
                    myh.SetBinContent(bin+1, mytemph.GetBinContent(bin))
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

        plot.makeplot([sumhistos_norm[start], histos[mva][0].Clone(), histos[mva][start], histos[mva][start+1]], ["Nom "+systs[start]]+[systs[s] for s in [0, start, start+1]], 
                      "N_{j}-7", mva+"_"+systs[start]+"_njets_norm", plotdir="./", linear=False, legendColumns=1, 
                      ymin=1e-5, ymax=1, ylabel="",
                      colors=[ROOT.kGray+1]+[colors[s] for s in [0, start, start+1]], norm=True, drawstyle="hist")


#for 2016, also copy over isr from the 2017 file
f_2017 = ROOT.TFile.Open("FitInputs/ttbar_systematics_2017_freezing.root")
h_isr_2017 = [f_2017.Get("%s_isr"%mva) for mva in mvas]
outputfile.cd()
for h in h_isr_2017:
    h.Write()
# outputfile.Close()
# import sys
# sys.exit()

# Look at results from Aron's fits -- update with my own background-only fits

#*(1) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of nominal samples*
#Results are here:
#fname_nom_shared = "/uscms_data/d2/soha/stealth/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/fit_results_v1_Dec19_2018/fitDiagnostics.root"
fname_nom_shared = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_nom_shared/fitDiagnostics.root"%fitversion
f_nom_shared = ROOT.TFile.Open(fname_nom_shared)
histos_nom_shared = {}
for mva in mvas:
    h = f_nom_shared.Get("shapes_fit_b/%s/TT" % mva)
    #h = f_nom_shared.Get("shapes_fit_s/sig%s/bkg_tt" % mva)
    histos_nom_shared[mva] = h

#*(2) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of nominal samples*
#Results are here:
#fname_nom_sep = "/uscms_data/d2/soha/stealth/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/fit_results_v2_Dec19_2018/fitDiagnostics.root"
fname_nom_sep = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_nom_sep/fitDiagnostics.root"%fitversion
f_nom_sep = ROOT.TFile.Open(fname_nom_sep)
histos_nom_sep = {}
for mva in mvas:
    h = f_nom_sep.Get("shapes_fit_b/%s/TT" % mva)
    #h = f_nom_sep.Get("shapes_fit_s/sig%s/bkg_tt" % mva)
    histos_nom_sep[mva] = h

#*(3) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of JEC UP samples*
#Results are here:
#fname_sys_shared = "/uscms_data/d2/soha/stealth/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/fit_results_v1_Dec20_2018/fitDiagnostics_JECUp.root"
fname_sys_shared = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_JECUp_shared/fitDiagnostics.root"%fitversion
f_sys_shared = ROOT.TFile.Open(fname_sys_shared)
histos_sys_shared = {}
for mva in mvas:
    h = f_sys_shared.Get("shapes_fit_b/%s/TT" % mva)
    #h = f_sys_shared.Get("shapes_fit_s/sig%s/bkg_tt" % mva)
    histos_sys_shared[mva] = h

#*(4) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of JEC UP samples*
#Results are here:
#fname_sys_sep = "/uscms_data/d2/soha/stealth/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/fit_results_v2_Dec20_2018/fitDiagnostics_JECUp.root"
fname_sys_sep = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_JECUp_sep/fitDiagnostics.root"%fitversion
f_sys_sep = ROOT.TFile.Open(fname_sys_sep)
histos_sys_sep = {}
for mva in mvas:
    h = f_sys_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_sys_sep[mva] = h

#*(5) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of JEC Down samples*
#Results are here:
fname_JECDown_shared = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_JECDown_shared/fitDiagnostics.root"%fitversion
f_JECDown_shared = ROOT.TFile.Open(fname_JECDown_shared)
histos_JECDown_shared = {}
for mva in mvas:
    h = f_JECDown_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_JECDown_shared[mva] = h

#*(6) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of JEC Down samples*
#Results are here:
fname_JECDown_sep = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_JECDown_sep/fitDiagnostics.root"%fitversion
f_JECDown_sep = ROOT.TFile.Open(fname_JECDown_sep)
histos_JECDown_sep = {}
for mva in mvas:
    h = f_JECDown_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_JECDown_sep[mva] = h




#*(7) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of JER UP samples*
#Results are here:
fname_JERUp_shared = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_JERUp_shared/fitDiagnostics.root"%fitversion
f_JERUp_shared = ROOT.TFile.Open(fname_JERUp_shared)
histos_JERUp_shared = {}
for mva in mvas:
    h = f_JERUp_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_JERUp_shared[mva] = h

#*(8) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of JER UP samples*
#Results are here:
# fname_JERUp_sep = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_JERUp_sep/fitDiagnostics.root"%fitversion
# f_JERUp_sep = ROOT.TFile.Open(fname_JERUp_sep)
# histos_JERUp_sep = {}
# for mva in mvas:
#     h = f_JERUp_sep.Get("shapes_fit_b/%s/TT" % mva)
#     histos_JERUp_sep[mva] = h

#*(9) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of JER Down samples*
#Results are here:
fname_JERDown_shared = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_JERDown_shared/fitDiagnostics.root"%fitversion
f_JERDown_shared = ROOT.TFile.Open(fname_JERDown_shared)
histos_JERDown_shared = {}
for mva in mvas:
    h = f_JERDown_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_JERDown_shared[mva] = h

#*(10) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of JER Down samples*
#Results are here:
fname_JERDown_sep = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_JERDown_sep/fitDiagnostics.root"%fitversion
f_JERDown_sep = ROOT.TFile.Open(fname_JERDown_sep)
histos_JERDown_sep = {}
for mva in mvas:
    h = f_JERDown_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_JERDown_sep[mva] = h


#*(11) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of FSR UP samples*
#Results are here:
fname_FSRUp_shared = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_fsrUp_shared/fitDiagnostics.root"%("freezing2017")
f_FSRUp_shared = ROOT.TFile.Open(fname_FSRUp_shared)
histos_FSRUp_shared = {}
for mva in mvas:
    h = f_FSRUp_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_FSRUp_shared[mva] = h

#*(12) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of FSR UP samples*
#Results are here:
fname_FSRUp_sep = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_fsrUp_sep/fitDiagnostics.root"%("freezing2017")
f_FSRUp_sep = ROOT.TFile.Open(fname_FSRUp_sep)
histos_FSRUp_sep = {}
for mva in mvas:
    h = f_FSRUp_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_FSRUp_sep[mva] = h

#*(13) Fit with a0, a1, a2 shared across MVA bins, fitting to pseudo-data composed of FSR Down samples*
#Results are here:
fname_FSRDown_shared = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_fsrDown_shared/fitDiagnostics.root"%("freezing2017")
f_FSRDown_shared = ROOT.TFile.Open(fname_FSRDown_shared)
histos_FSRDown_shared = {}
for mva in mvas:
    h = f_FSRDown_shared.Get("shapes_fit_b/%s/TT" % mva)
    histos_FSRDown_shared[mva] = h

#*(14) Fit with separate a0, a1, a2 for each MVA bin, fitting to pseudo-data composed of FSR Down samples*
#Results are here:
fname_FSRDown_sep = "~/nobackup/StealthRPV/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/%s_bkgonly_fsrDown_sep/fitDiagnostics.root"%("freezing2017")
f_FSRDown_sep = ROOT.TFile.Open(fname_FSRDown_sep)
histos_FSRDown_sep = {}
for mva in mvas:
    h = f_FSRDown_sep.Get("shapes_fit_b/%s/TT" % mva)
    histos_FSRDown_sep[mva] = h



#For each of the above files you can
#- Navigate to the TDirectoryFile called shapes_fit_s
#- Navigate to each TDirectoryFile called sigD1, sigD2, sigD3, sigD4 for the corresponding MVA bin
#- Extract information from bins 1 though 8 of the histogram called bkg_tt

# fname_rebin = "/uscms_data/d2/soha/stealth/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/Keras_V1.2.4/njets_rebin_for_Aron.root"
# inputfile_rebin = ROOT.TFile.Open(fname_rebin)
inputfile_rebin=inputfile # normal file now has correct bin edges
histos_rebin = {}
for mva in mvas: 
    # get all histograms
    histnames = [mva+"_"+base+syst for syst in systs]
    Dhistos = [inputfile_rebin.Get(histname) for histname in histnames]
    histos_rebin[mva] = Dhistos

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
    plot.makeplot([histos_nom_sep[mva].Clone(), histos_nom_shared[mva].Clone(), histos_rebin[mva][0].Clone()], ["Separate", "Shared", "MC"], "N_{j}-7", mva+"_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
    ymin=0.01, ymax=1e5, 
    ylabel="",
    colors=colors2, norm=False, drawstyle="lastP")

    # for syst JEC Up case
    plot.makeplot([histos_sys_sep[mva].Clone(), histos_sys_shared[mva].Clone(), histos_rebin[mva][1].Clone()], ["Separate", "Shared", "MC"], "N_{j}-7", mva+"_sys_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
    ymin=0.01, ymax=1e5, 
    ylabel="",
    colors=colors2, norm=False, drawstyle="lastP")
    
    # for JEC Down case
    plot.makeplot([histos_nom_shared[mva].Clone(), histos_JECDown_sep[mva].Clone(), histos_JECDown_shared[mva].Clone(), histos_rebin[mva][2].Clone()], 
                  ["Nominal","Separate", "Shared", "MC"], 
                  "N_{j}-7", mva+"_JECDown_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                  ymin=0.01, ymax=1e5, 
                  ylabel="",
                  colors=colors3, norm=False, drawstyle="lastP")
    
    # for JER Up case
    # plot.makeplot([histos_nom_shared[mva].Clone(), histos_JERUp_sep[mva].Clone(), histos_JERUp_shared[mva].Clone(), histos_rebin[mva][3].Clone()], 
    #               ["Nominal","Separate", "Shared", "MC"], 
    #               "N_{j}-7", mva+"_JERUp_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
    #               ymin=0.01, ymax=1e5, 
    #               ylabel="",
    #               colors=colors3, norm=False, drawstyle="lastP")
    
    # for JER Down case
    plot.makeplot([histos_nom_shared[mva].Clone(), histos_JERDown_sep[mva].Clone(), histos_JERDown_shared[mva].Clone(), histos_rebin[mva][4].Clone()], 
                  ["Nominal","Separate", "Shared", "MC"], 
                  "N_{j}-7", mva+"_JERDown_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                  ymin=0.01, ymax=1e5, 
                  ylabel="",
                  colors=colors3, norm=False, drawstyle="lastP")
   
    # for FSR Up case
    plot.makeplot([histos_nom_shared[mva].Clone(), histos_FSRUp_sep[mva].Clone(), histos_FSRUp_shared[mva].Clone(), histos_rebin[mva][3].Clone()], 
                  ["Nominal","Separate", "Shared", "MC"], 
                  "N_{j}-7", mva+"_FSRUp_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
                  ymin=0.01, ymax=1e5, 
                  ylabel="",
                  colors=colors3, norm=False, drawstyle="lastP")
    
    # for FSR Down case
    plot.makeplot([histos_nom_shared[mva].Clone(), histos_FSRDown_sep[mva].Clone(), histos_FSRDown_shared[mva].Clone(), histos_rebin[mva][4].Clone()], 
                  ["Nominal","Separate", "Shared", "MC"], 
                  "N_{j}-7", mva+"_FSRDown_fit"+"_njets", plotdir="./", linear=False, legendColumns=1, 
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
    ratio_sys = histos_sys_sep[mva].Clone(mva+"_ratio_sys_fit")
    #ratio_sys.Scale(1./ratio_sys.Integral())
    ratio_sys.Divide(histos_sys_shared[mva])
    # For the JEC Down case
    ratio_JECDown = histos_JECDown_sep[mva].Clone(mva+"_ratio_JECDown_fit")
    ratio_JECDown.Divide(histos_JECDown_shared[mva])
    # For the JER Up case
    # ratio_JERUp = histos_JERUp_sep[mva].Clone(mva+"_ratio_JERUp_fit")
    # ratio_JERUp.Divide(histos_JERUp_shared[mva])
    # For the JER Down case
    ratio_JERDown = histos_JERDown_sep[mva].Clone(mva+"_ratio_JERDown_fit")
    ratio_JERDown.Divide(histos_JERDown_shared[mva])
    # For the FSR Up case
    ratio_FSRUp = histos_FSRUp_sep[mva].Clone(mva+"_ratio_FSRUp_fit")
    ratio_FSRUp.Divide(histos_FSRUp_shared[mva])
    # For the JER Down case
    ratio_FSRDown = histos_FSRDown_sep[mva].Clone(mva+"_ratio_FSRDown_fit")
    ratio_FSRDown.Divide(histos_FSRDown_shared[mva])


    # Also make plots of difference wrt nominal
    divratio = ratio_sys.Clone(mva+"_ratio_div_fit")
    divratio.Divide(ratio_nom)
    divratio_JECDown = ratio_JECDown.Clone(mva+"_ratio_JECDown_div_fit")
    divratio_JECDown.Divide(ratio_nom)
    #divratio_JERUp = ratio_JERUp.Clone(mva+"_ratio_JERUp_div_fit")
    #divratio_JERUp.Divide(ratio_nom)
    divratio_JERDown = ratio_JERDown.Clone(mva+"_ratio_JERDown_div_fit")
    divratio_JERDown.Divide(ratio_nom)
    divratio_FSRUp = ratio_FSRUp.Clone(mva+"_ratio_FSRUp_div_fit")
    divratio_FSRUp.Divide(ratio_nom)
    divratio_FSRDown = ratio_FSRDown.Clone(mva+"_ratio_FSRDown_div_fit")
    divratio_FSRDown.Divide(ratio_nom)

    outputfile.cd()
    # Rather than put JECUp and JERUp, construct maximum of the Up and Down variations, but preserve the "sign"
    mva_JEC = divratio.Clone(mva+"_JEC")
    for i in range(divratio.GetNbinsX()):
        #print i
        up_value = divratio.GetBinContent(i+1)
        #print "up: " , up_value
        if up_value<1:
            up_value = 1/up_value
        down_value = divratio_JECDown.GetBinContent(i+1)
        #print "down: " , down_value
        if down_value<1:
            down_value = 1/down_value
        if down_value > up_value:
            #print "changing value"
            if divratio.GetBinContent(i+1) > 1:
                mva_JEC.SetBinContent(i+1,down_value)
                #print "to ", down_value
            else:
                mva_JEC.SetBinContent(i+1,1/down_value)
                #print "to " , 1/down_value

    mva_JER = divratio_JERDown.Clone(mva+"_JER")
    # for i in range(divratio_JERUp.GetNbinsX()):
    #     up_value = divratio_JERUp.GetBinContent(i+1)
    #     if up_value<1:
    #         up_value = 1/up_value
    #     down_value = divratio_JERDown.GetBinContent(i+1)
    #     if down_value<1:
    #         down_value = 1/down_value
    #     if down_value < up_value:
    #         if divratio_JERDown.GetBinContent(i+1) > 1:
    #             mva_JER.SetBinContent(i+1,up_value)
    #         else:
    #             mva_JER.SetBinContent(i+1,1/up_value)

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

        
    #divratio.Write(mva+"_JEC")
    #divratio_JERUp.Write(mva+"_JER")
    mva_JEC.Write(mva+"_JEC")
    mva_JER.Write(mva+"_JER")
    mva_FSR.Write(mva+"_FSR")


    # Also write the shape differences for the no-systematics case to the file
    ratio_nom.Write(mva+"_nom")

    #diffratio = ratio_sys.Clone(mva+"_ratio_diff_fit")
    #diffratio.Add(ratio_nom, -1)
    

    plot.makeplot([ratio_nom, ratio_sys, ratio_JECDown], ["Nominal", "JEC Up", "JEC Down"], "N_{j}-7", mva+"_fit"+"_comp", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.7, ymax=1.3, ylabel="",
    colors=colors, norm=False, drawstyle="hist")

    #plot.makeplot([ratio_nom, ratio_JERUp, ratio_JERDown], ["Nominal", "JER Up", "JER Down"], "N_{j}-7", mva+"_fit"+"_comp_JER", plotdir="./", linear=True, legendColumns=1, 
    #ymin=0.7, ymax=1.3, ylabel="",
    #colors=colors, norm=False, drawstyle="hist")
    plot.makeplot([ratio_nom, ratio_JERDown], ["Nominal", "JER Down"], "N_{j}-7", mva+"_fit"+"_comp_JER", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.7, ymax=1.3, ylabel="",
    colors=colors, norm=False, drawstyle="hist")

    plot.makeplot([ratio_nom, ratio_FSRUp, ratio_FSRDown], ["Nominal", "FSR Up", "FSR Down"], "N_{j}-7", mva+"_fit"+"_comp_FSR", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.7, ymax=1.3, ylabel="",
    colors=colors, norm=False, drawstyle="hist")

    
    plot.makeplot([divratio, divratio_JECDown], ["JEC Up","JEC Down"], "N_{j}-7", mva+"_fit_comp_div", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
    colors=colors, norm=False, drawstyle="hist")

    #plot.makeplot([divratio_JERUp, divratio_JERDown], ["JER Up","JER Down"], "N_{j}-7", mva+"_fit_comp_div_JER", plotdir="./", linear=True, legendColumns=1, 
    #ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
    #colors=colors, norm=False, drawstyle="hist")
    plot.makeplot([divratio_JERDown], ["JER Down"], "N_{j}-7", mva+"_fit_comp_div_JER", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
    colors=colors, norm=False, drawstyle="hist")

    plot.makeplot([divratio_FSRUp, divratio_FSRDown], ["FSR Up","FSR Down"], "N_{j}-7", mva+"_fit_comp_div_FSR", plotdir="./", linear=True, legendColumns=1, 
    ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
    colors=colors, norm=False, drawstyle="hist")


    #plot.makeplot([diffratio], ["JEC Up"], "N_{j}-7", mva+"_fit_comp_diff", plotdir="./", linear=True, legendColumns=1, 
    #ymin=-0.2, ymax=0.2, dropzeroes=False, ylabel="",
    #colors=colors, norm=False, drawstyle="hist")



outputfile.Close()


# Also look at some other systematics
        # KEY: TH1Dh_njets_1l_ge7j_ge1b;1h_njets_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_1l_ge7j_ge1b_NVtx0_20;1h_njets_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_1l_ge7j_ge1b_NVtx20_40;1h_njets_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_1l_ge7j_ge1b_NVtx40_60;1h_njets_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_1l_ge7j_ge1b_NVtx60_80;1h_njets_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_1l_ge7j_ge1b_d1;1h_njets_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_1l_ge7j_ge1b_d2;1h_njets_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_1l_ge7j_ge1b_d3;1h_njets_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_1l_ge7j_ge1b_d4;1h_njets_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_1l_ge7j_ge1b_noMbl;1h_njets_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b;1h_njets_FSRDown_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b_NVtx0_20;1h_njets_FSRDown_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b_NVtx20_40;1h_njets_FSRDown_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b_NVtx40_60;1h_njets_FSRDown_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b_NVtx60_80;1h_njets_FSRDown_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b_d1;1h_njets_FSRDown_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b_d2;1h_njets_FSRDown_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b_d3;1h_njets_FSRDown_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b_d4;1h_njets_FSRDown_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_FSRDown_1l_ge7j_ge1b_noMbl;1h_njets_FSRDown_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b;1h_njets_FSRUp_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b_NVtx0_20;1h_njets_FSRUp_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b_NVtx20_40;1h_njets_FSRUp_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b_NVtx40_60;1h_njets_FSRUp_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b_NVtx60_80;1h_njets_FSRUp_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b_d1;1h_njets_FSRUp_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b_d2;1h_njets_FSRUp_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b_d3;1h_njets_FSRUp_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b_d4;1h_njets_FSRUp_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_FSRUp_1l_ge7j_ge1b_noMbl;1h_njets_FSRUp_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b;1h_njets_ISRDown_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b_NVtx0_20;1h_njets_ISRDown_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b_NVtx20_40;1h_njets_ISRDown_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b_NVtx40_60;1h_njets_ISRDown_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b_NVtx60_80;1h_njets_ISRDown_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b_d1;1h_njets_ISRDown_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b_d2;1h_njets_ISRDown_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b_d3;1h_njets_ISRDown_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b_d4;1h_njets_ISRDown_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_ISRDown_1l_ge7j_ge1b_noMbl;1h_njets_ISRDown_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b;1h_njets_ISRUp_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b_NVtx0_20;1h_njets_ISRUp_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b_NVtx20_40;1h_njets_ISRUp_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b_NVtx40_60;1h_njets_ISRUp_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b_NVtx60_80;1h_njets_ISRUp_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b_d1;1h_njets_ISRUp_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b_d2;1h_njets_ISRUp_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b_d3;1h_njets_ISRUp_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b_d4;1h_njets_ISRUp_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_ISRUp_1l_ge7j_ge1b_noMbl;1h_njets_ISRUp_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b;1h_njets_PDFDown_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b_NVtx0_20;1h_njets_PDFDown_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b_NVtx20_40;1h_njets_PDFDown_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b_NVtx40_60;1h_njets_PDFDown_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b_NVtx60_80;1h_njets_PDFDown_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b_d1;1h_njets_PDFDown_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b_d2;1h_njets_PDFDown_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b_d3;1h_njets_PDFDown_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b_d4;1h_njets_PDFDown_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_PDFDown_1l_ge7j_ge1b_noMbl;1h_njets_PDFDown_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b;1h_njets_PDFUp_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b_NVtx0_20;1h_njets_PDFUp_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b_NVtx20_40;1h_njets_PDFUp_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b_NVtx40_60;1h_njets_PDFUp_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b_NVtx60_80;1h_njets_PDFUp_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b_d1;1h_njets_PDFUp_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b_d2;1h_njets_PDFUp_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b_d3;1h_njets_PDFUp_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b_d4;1h_njets_PDFUp_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_PDFUp_1l_ge7j_ge1b_noMbl;1h_njets_PDFUp_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b;1h_njets_pileupDown_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b_NVtx0_20;1h_njets_pileupDown_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b_NVtx20_40;1h_njets_pileupDown_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b_NVtx40_60;1h_njets_pileupDown_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b_NVtx60_80;1h_njets_pileupDown_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b_d1;1h_njets_pileupDown_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b_d2;1h_njets_pileupDown_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b_d3;1h_njets_pileupDown_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b_d4;1h_njets_pileupDown_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_pileupDown_1l_ge7j_ge1b_noMbl;1h_njets_pileupDown_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b;1h_njets_pileupUp_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b_NVtx0_20;1h_njets_pileupUp_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b_NVtx20_40;1h_njets_pileupUp_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b_NVtx40_60;1h_njets_pileupUp_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b_NVtx60_80;1h_njets_pileupUp_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b_d1;1h_njets_pileupUp_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b_d2;1h_njets_pileupUp_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b_d3;1h_njets_pileupUp_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b_d4;1h_njets_pileupUp_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_pileupUp_1l_ge7j_ge1b_noMbl;1h_njets_pileupUp_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b;1h_njets_pileup_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b_NVtx0_20;1h_njets_pileup_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b_NVtx20_40;1h_njets_pileup_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b_NVtx40_60;1h_njets_pileup_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b_NVtx60_80;1h_njets_pileup_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b_d1;1h_njets_pileup_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b_d2;1h_njets_pileup_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b_d3;1h_njets_pileup_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b_d4;1h_njets_pileup_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_pileup_1l_ge7j_ge1b_noMbl;1h_njets_pileup_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b;1h_njets_scaleDown_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b_NVtx0_20;1h_njets_scaleDown_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b_NVtx20_40;1h_njets_scaleDown_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b_NVtx40_60;1h_njets_scaleDown_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b_NVtx60_80;1h_njets_scaleDown_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b_d1;1h_njets_scaleDown_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b_d2;1h_njets_scaleDown_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b_d3;1h_njets_scaleDown_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b_d4;1h_njets_scaleDown_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_scaleDown_1l_ge7j_ge1b_noMbl;1h_njets_scaleDown_1l_ge7j_ge1b_noMbl
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b;1h_njets_scaleUp_1l_ge7j_ge1b
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b_NVtx0_20;1h_njets_scaleUp_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b_NVtx20_40;1h_njets_scaleUp_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b_NVtx40_60;1h_njets_scaleUp_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b_NVtx60_80;1h_njets_scaleUp_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b_d1;1h_njets_scaleUp_1l_ge7j_ge1b_d1
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b_d2;1h_njets_scaleUp_1l_ge7j_ge1b_d2
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b_d3;1h_njets_scaleUp_1l_ge7j_ge1b_d3
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b_d4;1h_njets_scaleUp_1l_ge7j_ge1b_d4
        # KEY: TH1Dh_njets_scaleUp_1l_ge7j_ge1b_noMbl;1h_njets_scaleUp_1l_ge7j_ge1b_noMbl
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b;1h_njets_NVtx_1l_ge7j_ge1b
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b_NVtx0_20;1h_njets_NVtx_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b_NVtx20_40;1h_njets_NVtx_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b_NVtx40_60;1h_njets_NVtx_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b_NVtx60_80;1h_njets_NVtx_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b_d1;1h_njets_NVtx_1l_ge7j_ge1b_d1
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b_d2;1h_njets_NVtx_1l_ge7j_ge1b_d2
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b_d3;1h_njets_NVtx_1l_ge7j_ge1b_d3
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b_d4;1h_njets_NVtx_1l_ge7j_ge1b_d4
        # KEY: TH2Dh_njets_NVtx_1l_ge7j_ge1b_noMbl;1h_njets_NVtx_1l_ge7j_ge1b_noMbl
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b;1h_njets_TrueNumInteractions_1l_ge7j_ge1b
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b_NVtx0_20;1h_njets_TrueNumInteractions_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b_NVtx20_40;1h_njets_TrueNumInteractions_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b_NVtx40_60;1h_njets_TrueNumInteractions_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b_NVtx60_80;1h_njets_TrueNumInteractions_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b_d1;1h_njets_TrueNumInteractions_1l_ge7j_ge1b_d1
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b_d2;1h_njets_TrueNumInteractions_1l_ge7j_ge1b_d2
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b_d3;1h_njets_TrueNumInteractions_1l_ge7j_ge1b_d3
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b_d4;1h_njets_TrueNumInteractions_1l_ge7j_ge1b_d4
        # KEY: TH2Dh_njets_TrueNumInteractions_1l_ge7j_ge1b_noMbl;1h_njets_TrueNumInteractions_1l_ge7j_ge1b_noMbl
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b;1h_njets_deepESM_1l_ge7j_ge1b
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b_NVtx0_20;1h_njets_deepESM_1l_ge7j_ge1b_NVtx0_20
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b_NVtx20_40;1h_njets_deepESM_1l_ge7j_ge1b_NVtx20_40
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b_NVtx40_60;1h_njets_deepESM_1l_ge7j_ge1b_NVtx40_60
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b_NVtx60_80;1h_njets_deepESM_1l_ge7j_ge1b_NVtx60_80
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b_d1;1h_njets_deepESM_1l_ge7j_ge1b_d1
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b_d2;1h_njets_deepESM_1l_ge7j_ge1b_d2
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b_d3;1h_njets_deepESM_1l_ge7j_ge1b_d3
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b_d4;1h_njets_deepESM_1l_ge7j_ge1b_d4
        # KEY: TH2Dh_njets_deepESM_1l_ge7j_ge1b_noMbl;1h_njets_deepESM_1l_ge7j_ge1b_noMbl

