#This script takes the mini tree outputs for the MC and data in both the signal and control region and makes histograms of njets in different categories for comparison.
#Comparisons include:
#   - Using CR derived bin edges (from SingleMuon dataset) and using SR derived bin edges (from TT MC)
#   - Looking at events with just muons in the signal region, just electrons in the signal region, just muons in the control region, or both leptons in the signal region
#   - Between the different years [2016,2017] - 2018 can be added of course

## Last edited by Kelvin Mei - September 24th, 2019

#!/bin/python
import ROOT
import copy
import os
import array

ROOT.gROOT.SetBatch( True )

def main() :

    #Set the DefaultSumw2 for good error propagation
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()
    
    # MC Input Files - make sure you have [2016,2017]_Other.root where Other is the histogram add (ROOT command hadd) of all the other backgrounds (you can break it down further if you want)
    categoryList      = [ "QCD", "TT", "Other", "Data" ]
    
    # Output directory check
    outputDir       = "njetRootHistos"
    if not os.path.exists(outputDir):
        print "Output directory does not exist - making directory",outputDir
        os.makedirs(outputDir)
    
    #Make histograms for each year, MVA bin, region, and type of lepton
    yearList        = [ "2016", "2017" ] # year of data
    regionList      = [ "CR", "SR" ] # pass CR baseline or pass SR baseline
    binEdgeList     = [ "CRBE", "SRBE" ] # Using control region derived bin edges or signal region bin edges
    leptonList      = [ "el", "mu", "both" ] #Which lepton category
    mvaBinList      = [ "D1", "D2", "D3", "D4", "All" ] #Which MVA bin

    for category in categoryList: #A separate file will be made for each MC breakdown

        #Define dictionary with all the histograms
        njetHistoDict   = {}
    
        for region in regionList :
            njetHistoDict[region] = {}
            
            for year in yearList :
                njetHistoDict[region][year] = {}
    
                for binEdge in binEdgeList:
    
                    #We will never look at signal region using control region Single Muon data derived bin edges
                    if binEdge == "CRBE" and region == "SR" :
                        continue

                    njetHistoDict[region][year][binEdge] = {}
                
                    for lepton in leptonList :
                        
                        #Ignore anything with electrons in the control region (of course)
                        if lepton == "el" and region == "CR" :
                            continue

                        if lepton == "both" and region == "CR" :
                            continue

                        else :
                            njetHistoDict[region][year][binEdge][lepton] = {}
        
                            for mvaBin in mvaBinList :
                                njetHistoDict[region][year][binEdge][lepton][mvaBin] = ROOT.TH1D( "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin, "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin, 6, 7, 13 )
    
        #Get the input file mini trees
        for region in regionList :
            inputDir    = ""
            if region == "CR" :
                inputDir = "../../haddfiles/QCDCRMiniTrees/"
            if region == "SR" :
                inputDir    = "../../haddfiles/SignalMiniTrees/"

            if not os.path.exists(inputDir):
                print "Error: you do not have the correct path/directory to the miniTrees for the ",region
            for year in yearList :
                
                for binEdge in binEdgeList :
                    
                    if binEdge == "CRBE" and region == "SR" :
                        continue
                    
                    print "Initializing and making histograms for",category,"passing the",region,"baseline in the year",year,"using bin edges derived from the",binEdge[:2]
                    
                    
                    binEdges        = binEdgeDictReturn( binEdge, year )
                    inputFile       = ROOT.TFile.Open( inputDir+year+"_"+category+".root")
                    myTree          = inputFile.Get( "myMiniTree" )

                    for entry in xrange( 0, myTree.GetEntries() ) : 
                        myTree.GetEntry(entry)
                        eventWeight = 1.0
                        
                        if region == "CR" :
                            #Legacy requirement where the control region mini trees had events with less than 7 jets to see trends in the bulk
                            if myTree.NGoodJets_pt30 < 7:
                                continue

                            if category != "Data" :
                                eventweight = myTree.totalEventWeightNIM*myTree.Lumi
                            else : 
                                eventweight = 1.0
                        
                            njetHistoDict[region][year][binEdge]["mu"]["All"].Fill( myTree.NGoodJets_pt30, eventweight ) 
                            
                            if myTree.NGoodJets_pt30 < 13 :
                                for itBinEdge in xrange( 0, len( binEdges[ myTree.NGoodJets_pt30 ] ) ) :
                                    if myTree.deepESM_val < binEdges[ myTree.NGoodJets_pt30 ][ itBinEdge ] :
                                        njetHistoDict[region][year][binEdge]["mu"][mvaBinList[itBinEdge]].Fill( myTree.NGoodJets_pt30, eventweight )
                                        break

                            else :
                                for itBinEdge in xrange( 0, len( binEdges[ 12 ] ) ):
                                    if myTree.deepESM_val < binEdges[ 12 ][ itBinEdge ] :
                                        njetHistoDict[region][year][binEdge]["mu"][mvaBinList[itBinEdge]].Fill( 12, eventweight )
                                        break
                        
                        elif region == "SR" :
                            
                            if myTree.NGoodJets_pt30 < 7:
                                continue
                            
                            if category != "Data" :
                                eventweight = myTree.totalEventWeight*myTree.Lumi
                            else :
                                eventweight = 1.0

                            njetHistoDict[region][year][binEdge]["both"]["All"].Fill( myTree.NGoodJets_pt30, eventweight )
                            if myTree.NGoodMuons == 1:
                                njetHistoDict[region][year][binEdge]["mu"]["All"].Fill( myTree.NGoodJets_pt30, eventweight )
                            if myTree.NGoodElectrons == 1:
                                njetHistoDict[region][year][binEdge]["el"]["All"].Fill( myTree.NGoodJets_pt30, eventweight )
                                
                            if myTree.NGoodJets_pt30 < 13 :
                                for itBinEdge in xrange( 0, len( binEdges[ myTree.NGoodJets_pt30 ] ) ) :
                                    if myTree.deepESM_val < binEdges[ myTree.NGoodJets_pt30 ][ itBinEdge ] :
                                        njetHistoDict[region][year][binEdge]["both"][mvaBinList[itBinEdge]].Fill( myTree.NGoodJets_pt30, eventweight )
                                        if myTree.NGoodMuons == 1:
                                            njetHistoDict[region][year][binEdge]["mu"][mvaBinList[itBinEdge]].Fill( myTree.NGoodJets_pt30, eventweight )
                                        if myTree.NGoodElectrons == 1:
                                            njetHistoDict[region][year][binEdge]["el"][mvaBinList[itBinEdge]].Fill( myTree.NGoodJets_pt30, eventweight )
                                        break

                            else :
                                for itBinEdge in xrange( 0, len( binEdges[ 12 ] ) ) :
                                    if myTree.deepESM_val < binEdges[ 12 ][ itBinEdge ] :
                                        njetHistoDict[region][year][binEdge]["both"][mvaBinList[itBinEdge]].Fill( 12, eventweight )
                                        if myTree.NGoodMuons == 1:
                                            njetHistoDict[region][year][binEdge]["mu"][mvaBinList[itBinEdge]].Fill( 12, eventweight )
                                        if myTree.NGoodElectrons == 1:
                                            njetHistoDict[region][year][binEdge]["el"][mvaBinList[itBinEdge]].Fill( 12, eventweight )
                                        break
                        
        
        outputFile = ROOT.TFile( outputDir+"/"+category+".root", "RECREATE")
        for region in regionList :
            for year in yearList :
                for binEdge in binEdgeList:
                    if binEdge == "CRBE" and region == "SR" :
                        continue
                    for lepton in leptonList :
                        if lepton == "el" and region == "CR" :
                            continue
                        if lepton == "both" and region == "CR" :
                            continue
                        else :
                            for mvaBin in mvaBinList :
                                njetHistoDict[region][year][binEdge][lepton][mvaBin].Write()
        
        outputFile.Close()
    
def binEdgeDictReturn( region, year ) :
    
    binEdges = {}
   
    # These are the bin edges derived from the QCD control region using data from the Single Muon data set
    if region == "CRBE" and year == "2017" :
        binEdges[7]     = [ 0.391, 0.771, 0.869, 1.000 ]
        binEdges[8]     = [ 0.413, 0.808, 0.892, 1.000 ]
        binEdges[9]     = [ 0.432, 0.832, 0.913, 1.000 ]
        binEdges[10]    = [ 0.437, 0.863, 0.932, 1.000 ]
        binEdges[11]    = [ 0.446, 0.880, 0.942, 1.000 ]
        binEdges[12]    = [ 0.471, 0.848, 0.959, 1.000 ]
        binEdges[13]    = [ 0.471, 0.848, 0.959, 1.000 ]
        binEdges[14]    = [ 0.471, 0.848, 0.959, 1.000 ]
                
    if region == "CRBE" and year == "2016" :
        binEdges[7]     = [ 0.408, 0.770, 0.885, 1.000 ]
        binEdges[8]     = [ 0.429, 0.814, 0.908, 1.000 ]
        binEdges[9]     = [ 0.455, 0.848, 0.933, 1.000 ]
        binEdges[10]    = [ 0.468, 0.862, 0.937, 1.000 ]
        binEdges[11]    = [ 0.506, 0.892, 0.948, 1.000 ]
        binEdges[12]    = [ 0.506, 0.892, 0.952, 1.000 ]
        binEdges[13]    = [ 0.506, 0.892, 0.952, 1.000 ]
        binEdges[14]    = [ 0.506, 0.892, 0.952, 1.000 ]
            
    # These bin edges are derived from the signal region using TT bar MC
    if region == "SRBE" and year == "2017" :
        binEdges[7]     = [ 0.346, 0.715, 0.833, 1.000 ]
        binEdges[8]     = [ 0.363, 0.751, 0.863, 1.000 ]
        binEdges[9]     = [ 0.376, 0.781, 0.886, 1.000 ]
        binEdges[10]    = [ 0.385, 0.796, 0.896, 1.000 ]
        binEdges[11]    = [ 0.391, 0.796, 0.900, 1.000 ]
        binEdges[12]    = [ 0.391, 0.796, 0.900, 1.000 ]
        binEdges[13]    = [ 0.391, 0.796, 0.900, 1.000 ]
        binEdges[14]    = [ 0.391, 0.796, 0.900, 1.000 ]
            
    if region == "SRBE" and year == "2016" :
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
