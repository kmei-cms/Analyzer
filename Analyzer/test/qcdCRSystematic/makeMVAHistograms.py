#This script takes the output of MakeMiniTrees and generates the first step to the QCD control region systematic
#Edited: Kelvin Mei on September 4th, 2019

import math
import ROOT
import array
import random
import copy
import os

ROOT.gROOT.SetBatch(True)

from optparse import OptionParser

parser = OptionParser()

#This is a boolean that is true for 2017 and false for 2016. This can be changed into a string to accompany all three years.
parser.add_option( '--year', action = 'store_true',
                    dest = 'year',
                    default = False, help = 'Use 2017 data file instead of 2016' )

#The signal option uses the bin edges derived in the TT bar MC that will also be used in the signal region. This is the option you will use if you want to make the QCD control region systematic
parser.add_option( '--signal', action = 'store_true',
                    dest = 'signal',
                    default = False, help = 'Use bin edges derived in signal TT bar region')

#The data option is used when taking the bin edges derived using binEdgeAnalyzer.py on the Single Muon data set. This was used for studies where you derive the bin edges on QCD and see how the fit function we used works.
parser.add_option( '--data', action = 'store_true',
                    dest = 'data',
                    default = False, help = 'Use bin edges derived with Single Muon data in QCD control region')

#Use the signal region baseline mini trees to make TT MVA histograms
parser.add_option( '--ttbar', action = 'store_true',
                    dest = 'ttbar',
                    default = False, help = 'Run over the signal region baseline tt bar mini trees to make the proper histograms')


(options, args) = parser.parse_args()

def main() :

    #Define variables needed for the analyzer here.
    ROOT.TH1.SetDefaultSumw2()

    #This is the locaiton of the miniTrees
    inputDir            = '../../haddfiles/QCDCRMiniTrees/'

    if options.ttbar :
        inputDir        = '../../haddfiles/SignalMiniTrees/'

    if not os.path.exists(inputDir):
        print "Directory to miniTrees do not exist"

    outputDir           = 'mvaPlots'

    if not os.path.exists(outputDir):
        print "Output directory does not exist, so making output directory", outputDir
        os.makedirs(outputDir)
    
    #Utilize this for the year option and to swap between different year files
    yearTag        = '2016'
    if( options.year ) :
        yearTag    = '2017'
    
    #Originally, this code was used to derive plots also in the QCD_Mu samples and up to 14 njets bins

    #mainBkgdArray       = [ "QCD", "QCD_Mu" ]
    mainBkgdArray       = [ "QCD" ]
    dataArray           = [ "Data_SingleMuon" ]
    #njetArray           = [ 7, 8, 9, 10, 11, 12, 13, 14 ]
    njetArray           = [ 7, 8, 9, 10, 11 ] #Because of statistics, use the values for Nj = 11 for Nj = 12 and up
    
    if options.ttbar :
        mainBkgdArray   = [ "TT" ]
        dataArray       = []
    
    #Choose the correct version of the bin edges (can add more options for more different bin edges)
    if options.data and options.signal :
        print "Cannot have two options for bin edges"
        return 0

    if not options.data and not options.signal :
        print "Must choose a set of bin edges"
        return 0
    
    #Define bin edges here
    binEdges            = {}

    if( options.year ):

        if( options.data ) :

            #2017 Data in control reigon derived bin edges - last updated in May 2019 and have to be updated with htDerivedweight applied to the QCDCR
            binEdges[7]     = [ 0.391, 0.771, 0.869, 1.000 ]
            binEdges[8]     = [ 0.413, 0.808, 0.892, 1.000 ]
            binEdges[9]     = [ 0.432, 0.832, 0.913, 1.000 ]
            binEdges[10]    = [ 0.437, 0.863, 0.932, 1.000 ]
            binEdges[11]    = [ 0.446, 0.880, 0.942, 1.000 ]
            binEdges[12]    = [ 0.471, 0.848, 0.959, 1.000 ]
            binEdges[13]    = [ 0.471, 0.848, 0.959, 1.000 ]
            binEdges[14]    = [ 0.471, 0.848, 0.959, 1.000 ]
            
        elif( options.signal ) :
            #2017 TT MC in signal region derived bin edges - HT scale factor added - checked to be up to date as of September 2019
            binEdges[7]     = [ 0.346, 0.715, 0.833, 1.000 ]
            binEdges[8]     = [ 0.363, 0.751, 0.863, 1.000 ]
            binEdges[9]     = [ 0.376, 0.781, 0.886, 1.000 ]
            binEdges[10]    = [ 0.385, 0.796, 0.896, 1.000 ]
            binEdges[11]    = [ 0.391, 0.796, 0.900, 1.000 ]
            binEdges[12]    = [ 0.391, 0.706, 0.900, 1.000 ]
            
    else :
        if( options.data ) :
            #2016 Data in control reigon derived bin edges - last updated in May 2019 and have to be updated with htDerivedweight applied to the QCDCR
            binEdges[7]     = [ 0.408, 0.770, 0.885, 1.000 ]
            binEdges[8]     = [ 0.429, 0.814, 0.908, 1.000 ]
            binEdges[9]     = [ 0.455, 0.848, 0.933, 1.000 ]
            binEdges[10]    = [ 0.468, 0.862, 0.937, 1.000 ]
            binEdges[11]    = [ 0.506, 0.892, 0.948, 1.000 ]
            binEdges[12]    = [ 0.506, 0.892, 0.952, 1.000 ]
            binEdges[13]    = [ 0.506, 0.892, 0.952, 1.000 ]
            binEdges[14]    = [ 0.506, 0.892, 0.952, 1.000 ]
 
        elif( options.signal ) :
            #2016 TT MC in signal region derived bin edges - HT scale factor added - checked to be up to date as of September 2019
            binEdges[7]     = [ 0.349, 0.678, 0.833, 1.000 ]
            binEdges[8]     = [ 0.357, 0.737, 0.874, 1.000 ]
            binEdges[9]     = [ 0.369, 0.783, 0.899, 1.000 ]
            binEdges[10]    = [ 0.383, 0.810, 0.913, 1.000 ]
            binEdges[11]    = [ 0.410, 0.848, 0.937, 1.000 ]
            binEdges[12]    = [ 0.436, 0.853, 0.944, 1.000 ]
    
    #Create a dictionary to hold onto all the histograms that will be used to write to an output file
    histogramDict                   = {}
    
    #Make all the input MVA histograms for MC
    for mainBkgd in mainBkgdArray :
        
        inputFileName   = inputDir+yearTag+"_"+mainBkgd+".root"
        inputFile       = ROOT.TFile.Open( inputFileName )
        inputTree       = ROOT.TTree()
        inputTree       = inputFile.Get( "myMiniTree" )
       
        histogramDict[mainBkgd]             = {}
        histogramDict[mainBkgd]["All"]      = ROOT.TH1D( "h_mva_all", "h_mva_all", 24, -0.1, 1.1 )
        
        for njet in njetArray :
            histogramDict[mainBkgd][njet]   = ROOT.TH1D( "h_mva_"+str(njet)+"j", "h_mva_"+str(njet)+"j", 24, -0.1, 1.1 )
        
        for entry in xrange( 0, inputTree.GetEntries() ) :
            inputTree.GetEntry( entry )

            if inputTree.NGoodJets_pt30 < 7 :
                continue
            
            histogramDict[mainBkgd]["All"].Fill( inputTree.deepESM_val, inputTree.Weight*inputTree.Lumi*inputTree.puWeightCorr*inputTree.htDerivedweight )

            if inputTree.NGoodJets_pt30 < 12 :
                for njet in njetArray :
                    if inputTree.NGoodJets_pt30 == njet :
                        histogramDict[mainBkgd][njet].Fill( inputTree.deepESM_val, inputTree.Weight*inputTree.Lumi*inputTree.puWeightCorr*inputTree.htDerivedweight )
                        break

            else:
                histogramDict[mainBkgd][11].Fill( inputTree.deepESM_val, inputTree.Weight*inputTree.Lumi*inputTree.puWeightCorr*inputTree.htDerivedweight )

        outputFile  = ROOT.TFile( outputDir+"/"+mainBkgd+"_ratio_"+yearTag+".root", "RECREATE" )

        histogramDict[mainBkgd]["All"].Write()
        for njet in njetArray :
            histogramDict[mainBkgd][njet].Write()
        outputFile.Close()
        inputFile.Close()
    
    #Make all the input MVA histograms for data
    for data in dataArray :

        inputFileName   = inputDir+yearTag+"_"+data+".root"
        inputFile       = ROOT.TFile.Open( inputFileName )
        inputTree       = ROOT.TTree()
        inputTree       = inputFile.Get( "myMiniTree" )
        
        histogramDict[data]             = {}
        histogramDict[data]["All"]      = ROOT.TH1D( "h_mva_all", "h_mva_all", 24, -0.1, 1.1 )
        
        for njet in njetArray :
            histogramDict[data][njet]   = ROOT.TH1D( "h_mva_"+str(njet)+"j", "h_mva_"+str(njet)+"j", 24, -0.1, 1.1 )
        
        for entry in xrange( 0, inputTree.GetEntries() ) :
            inputTree.GetEntry( entry )

            histogramDict[data]["All"].Fill( inputTree.deepESM_val, 1.0 )

            if inputTree.NGoodJets_pt30 < 12 :
                for njet in njetArray :
                    if inputTree.NGoodJets_pt30 == njet :
                        histogramDict[data][njet].Fill( inputTree.deepESM_val, 1.0 )
                        break
            else:
                histogramDict[data][11].Fill( inputTree.deepESM_val, 1.0 )
        
        outputFile  = ROOT.TFile( outputDir+"/"+data+"_ratio_"+yearTag+".root", "RECREATE" )

        histogramDict[data]["All"].Write()
        for njet in njetArray :
            histogramDict[data][njet].Write()
        outputFile.Close()
        inputFile.Close()

if __name__ == '__main__' :
    main()
