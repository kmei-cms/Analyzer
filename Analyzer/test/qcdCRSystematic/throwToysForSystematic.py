#This is the scripts that will throw n (default n = 100 ) number of toys using an additional multiplicative weight derived from using the data from the control region and the Powheg TT MC in the signal region

#Last updated September 23rd, 2019 by Kelvin Mei

#!/bin/python
import ROOT
import copy
import os.path
import array
import numpy as np

ROOT.gROOT.SetBatch( True )

from optparse import OptionParser

parser  = OptionParser()

parser.add_option('--log', action='store_true',
                    dest='log',
                    default = False, help = 'Make all plots log on the y scale' )

parser.add_option('--year', action='store_true',
                    dest='year',
                    default = False, help = 'Choose between 2016 and 2017 (true for 2017)' )

(options,args) = parser.parse_args()

def main():

    # Make sure this is set so errors are properly propagated
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()

    #Define some useful variables here

    yearTag = "2016"
    if options.year :
        yearTag = "2017"
    njetNameArray                       = [ 7, 8, 9, 10, 11, 12 ] #Changed with the move from 14 to 12 jets
    mvaCutArray                         = []
    mvaBinArray                         = [ "D1", "D2", "D3", "D4" ]
    seedNumber                          = 689300203
    
    #This is the location of the signal miniTrres
    inputDir                            = "../../haddfiles/SignalMiniTrees/"
    inputFile                           = ROOT.TFile.Open( inputDir+yearTag+"_TT.root" )
    inputTree                           = inputFile.Get( "myMiniTree" )

    #This corresponds to the binning of the MVA histograms
    #CHECK THIS TO MAKE SURE THIS IS RIGHT IF THE MVA HISTOGRAM BINNING CHANGES (can also streamline)
    for i in range( 24 ):
        mvaCutArray.append( float(i-1)*0.05 )
   
    #Define numbenr of toys that you want to throw here. The output file will have all
    # of the D1-D4 histograms for each toy
    #For 2016, 100 toys takes about 30 - 60 minutes running locally on a 2015 Macbook
    nToys                               = 100

    #Get all the information for the ratio of ratios and place them in a dictionary for easy access. The SF dictionary will hold the central value and the Err dictionary will hold the error bar
    ratioOfRatiosSF                     = {}
    ratioOfRatiosErr                    = {}
    inputFileName                       = "ratioOfRatios_"+yearTag+".root"
    inputRatioOfRatiosFile              = ROOT.TFile.Open( inputFileName )
        
    for njet in njetNameArray :
        if njet > 11 : #Decision was made to do the mva plots only up to njet = 11 for statistics reasons
            continue
        ratioOfRatioHisto                       = inputRatioOfRatiosFile.Get( "h_ratio_DataCR_over_TTSR_"+str(njet) )
        ratioOfRatiosSF[njet]           = {}
        ratioOfRatiosErr[njet]          = {}
        
        for binIter in range( ratioOfRatioHisto.GetNbinsX() ) :
           ratioOfRatiosSF[njet][mvaCutArray[binIter]] = ratioOfRatioHisto.GetBinContent(binIter+1)
           ratioOfRatiosErr[njet][mvaCutArray[binIter]] = ratioOfRatioHisto.GetBinError(binIter+1)
    
    inputRatioOfRatiosFile.Close()

    #Define two dictionaries to hold plots to validate procedure and the first attempt at calculating the systematic. The 1D histogram dictionary will hold the mean and standard deviation values (just for sanity checks) and the 2D histogram dictionary will hold the plot showing the extra chosen multiplicative factor.
    histoDict                           = {}
    histo2dDict                         = {}

    #Make plots 
    for njet in njetNameArray :
        
        #Histogram filled with all 100 toys with the MVA score
        histoDict[njet]                 = ROOT.TH1D( "h_mva_100toys_"+str(njet)+"j", "h_mva_100toys_"+str(njet)+"j", 20, 0, 1.0 )
        
        histo2dDict[njet]               = ROOT.TH2D( "h_mva_eventweight_"+str(njet)+"j", "h_mva_eventweight_"+str(njet)+"j", 20, 0, 1.0, 200, 0.0, 2.0 )
    
    histoDict["std"]                    = {}

    for mva in mvaBinArray :
        histoDict[mva]                  = ROOT.TH1D( "h_"+mva+"_qcdCR_meanValue", "h_"+mva+"_qcdCR_meanValue", 6, 7, 13 )
        histoDict["std"][mva]           = ROOT.TH1D( "h_"+mva+"_qcdCR_stdDev", "h_"+mva+"_qcdCR_stdDev", 6, 7, 13 )
    
    #Define bin edges based on the function. "signal" means derived from tt bar signal region.
    binEdges                            = binEdgeDictReturn( "Signal", yearTag )

    #Define some more dictionaries for the second attempt to calculate the systematic
    #   totalEventArray holds all the final number of events in each njet bin and mva bin (24)
    #       This is useful for using numpy functions
    #   toyHistoDict holds all the D1-D4 njet histograms for each toy trial

    totalEventArray                         = {}
    toyHistoDict                            = {}
    for i in mvaBinArray :
        totalEventArray[i]                  = {}
        for j in njetNameArray :
            totalEventArray[i][j]           = []
            for k in xrange( 0, nToys ) :
                totalEventArray[i][j].append(0.0)

    for i in xrange( 0, nToys ) :
        toyHistoDict[i]                     = {}
        for mva in mvaBinArray :
            toyHistoDict[i][mva]            = ROOT.TH1D( "h_"+mva+"_toy"+str(i), "h_"+mva+"_toy"+str(i), 6, 7, 13 )

    for entry in xrange( 0, inputTree.GetEntries() ):

        inputTree.GetEntry( entry )

        #Out of habit, initialize variables (and also to clear the previous value)
        scaleFactorMean                 = 0.0
        scaleFactorSigma                = 1.0
        mvaBin                          = ""

        if inputTree.NGoodJets_pt30 < 7 : #Is here only because I had a set of trees that looked at njet = 5 and greater in order to do HT weight fits
            continue

        if inputTree.NGoodJets_pt30 < 12 :
            
            for njet in njetNameArray :    
                #Determine which MVA bin this event falls into
                for binEdgeIter in range( len( mvaBinArray ) ) :
                    if inputTree.deepESM_val < binEdges[njet][binEdgeIter] :
                        mvaBin          = mvaBinArray[binEdgeIter]
                        break

                #Determine what the new event based weight is 
                if inputTree.NGoodJets_pt30 == njet:
                    for mvaBinIter in xrange( 0, len( mvaCutArray ) ) :
                        if inputTree.deepESM_val < mvaCutArray[mvaBinIter] :
                            scaleFactorMean     = ratioOfRatiosSF[njet][mvaCutArray[mvaBinIter]]
                            scaleFactorSigma    = ratioOfRatiosErr[njet][mvaCutArray[mvaBinIter]]
                            break

                    for i in xrange( 0, nToys ) :
                        myRandom                        = ROOT.TRandom3(seedNumber)
                        seedNumber                      += 1
                        newEventWeightSF                = myRandom.Gaus( scaleFactorMean, scaleFactorSigma )
                        histo2dDict[njet].Fill( inputTree.deepESM_val, newEventWeightSF )
                        histoDict[njet].Fill( inputTree.deepESM_val, newEventWeightSF )
                        histoDict[mvaBin].Fill( inputTree.NGoodJets_pt30, inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF )
                        totalEventArray[mvaBin][inputTree.NGoodJets_pt30][i] += inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF
                        toyHistoDict[i][mvaBin].Fill( inputTree.NGoodJets_pt30, inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF )
                    
                    break
        
        else : #For all events with more than 11 jets
            
            #Determine which MVA bin this event falls into - use njets >= 12 bin edges
            for binEdgeIter in range( len( mvaBinArray ) ) :
                if inputTree.deepESM_val < binEdges[12][binEdgeIter] :
                    mvaBin          = mvaBinArray[binEdgeIter]
                    break
            #Determine what the new event based weight is using njets >= 11 ratio of ratios histogram
            for mvaBinIter in xrange( 0, len( mvaCutArray ) ) :
                if inputTree.deepESM_val < mvaCutArray[mvaBinIter] :
                    scaleFactorMean     = ratioOfRatiosSF[11][mvaCutArray[mvaBinIter]]
                    scaleFactorSigma    = ratioOfRatiosErr[11][mvaCutArray[mvaBinIter]]
                    break

            for i in xrange( 0, nToys ) :
                myRandom                        = ROOT.TRandom3(seedNumber)
                seedNumber                      += 1
                newEventWeightSF                = myRandom.Gaus( scaleFactorMean, scaleFactorSigma )
                histo2dDict[12].Fill( inputTree.deepESM_val, newEventWeightSF )
                histoDict[12].Fill( inputTree.deepESM_val, newEventWeightSF )
                histoDict[mvaBin].Fill( 12, inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF )
                totalEventArray[mvaBin][12][i] += inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF
                toyHistoDict[i][mvaBin].Fill( 12, inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF )


    #This is to make plots of the RMS of the 100 toys using numpy functions
    rmsHistoDict                        = {}
    for mva in mvaBinArray:
        rmsHistoDict[mva]               = {}
        for njet in njetNameArray:
            y = np.array( totalEventArray[mva][njet] )

            print mva, histoDict[mva].GetBinLowEdge( njet-6 ), np.mean(y)
            rmsHistoDict[mva][njet]     = ROOT.TH1D( "h_"+mva+"_rms_"+str(njet)+"j", "h_"+mva+"_rms_"+str(njet)+"j", 100, -10.0, 10.0 )
            for iteration in totalEventArray[mva][njet] :
                rmsHistoDict[mva][njet].Fill( iteration - np.mean(y) )
            histoDict["std"][mva].SetBinContent( njet-6, np.std(y) )
   
    outputFileName = yearTag+"_evtWghtSyst_allToys.root"
    outputFile  = ROOT.TFile.Open( outputFileName, "RECREATE" )
    
    for njet in njetNameArray :
        histoDict[njet].Write()
        histo2dDict[njet].Write()
    
    for mva in mvaBinArray :
        histoDict[mva].Write()
        histoDict["std"][mva].Write()
    
        for njet in njetNameArray:
            rmsHistoDict[mva][njet].Write() 
    
    for i in xrange(0, nToys):
        for mva in mvaBinArray :
            toyHistoDict[i][mva].Write()

    outputFile.Close()

def binEdgeDictReturn( region, year ) :
    
    binEdges = {}
   
    # These are the bin edges derived from the QCD control region using data from the Single Muon data set
    if region == "QCDCR" and year == "2017" :
        binEdges[7]     = [ 0.391, 0.771, 0.869, 1.000 ]
        binEdges[8]     = [ 0.413, 0.808, 0.892, 1.000 ]
        binEdges[9]     = [ 0.432, 0.832, 0.913, 1.000 ]
        binEdges[10]    = [ 0.437, 0.863, 0.932, 1.000 ]
        binEdges[11]    = [ 0.446, 0.880, 0.942, 1.000 ]
        binEdges[12]    = [ 0.471, 0.848, 0.959, 1.000 ]
        binEdges[13]    = [ 0.471, 0.848, 0.959, 1.000 ]
        binEdges[14]    = [ 0.471, 0.848, 0.959, 1.000 ]
                
    if region == "QCDCR" and year == "2016" :
        binEdges[7]     = [ 0.408, 0.770, 0.885, 1.000 ]
        binEdges[8]     = [ 0.429, 0.814, 0.908, 1.000 ]
        binEdges[9]     = [ 0.455, 0.848, 0.933, 1.000 ]
        binEdges[10]    = [ 0.468, 0.862, 0.937, 1.000 ]
        binEdges[11]    = [ 0.506, 0.892, 0.948, 1.000 ]
        binEdges[12]    = [ 0.506, 0.892, 0.952, 1.000 ]
        binEdges[13]    = [ 0.506, 0.892, 0.952, 1.000 ]
        binEdges[14]    = [ 0.506, 0.892, 0.952, 1.000 ]
            
    # These bin edges are derived from the signal region using TT bar MC
    if region == "Signal" and year == "2017" :
        binEdges[7]     = [ 0.346, 0.715, 0.833, 1.000 ]
        binEdges[8]     = [ 0.363, 0.751, 0.863, 1.000 ]
        binEdges[9]     = [ 0.376, 0.781, 0.886, 1.000 ]
        binEdges[10]    = [ 0.385, 0.796, 0.896, 1.000 ]
        binEdges[11]    = [ 0.391, 0.796, 0.900, 1.000 ]
        binEdges[12]    = [ 0.391, 0.796, 0.900, 1.000 ]
        binEdges[13]    = [ 0.391, 0.796, 0.900, 1.000 ]
        binEdges[14]    = [ 0.391, 0.796, 0.900, 1.000 ]
            
    if region == "Signal" and year == "2016" :
        binEdges[7]     = [ 0.349, 0.678, 0.833, 1.000 ]
        binEdges[8]     = [ 0.357, 0.737, 0.874, 1.000 ]
        binEdges[9]     = [ 0.369, 0.783, 0.899, 1.000 ]
        binEdges[10]    = [ 0.383, 0.810, 0.913, 1.000 ]
        binEdges[11]    = [ 0.410, 0.848, 0.937, 1.000 ]
        binEdges[12]    = [ 0.436, 0.853, 0.944, 1.000 ]
        binEdges[13]    = [ 0.436, 0.853, 0.944, 1.000 ]
        binEdges[14]    = [ 0.436, 0.853, 0.944, 1.000 ]

    return binEdges

if __name__ == '__main__' :
    main()
