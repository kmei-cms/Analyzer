#All this script does is take the histograms produced from makeQcdCrComparisonHistos.py and takes the data and subtract out the TT and Other MC component to look at what the "pure QCD" data looks like

## Last edited by Kelvin Mei on September 24th, 2019

#!/bin/python
import ROOT
import copy
import os
import array

ROOT.gROOT.SetBatch( True )

def main() :
    
    #Set the defaultsumw2
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()
    
    #Make histograms for each year, MVA bin, region, and type of lepton
    yearList        = [ "2016" , "2017" ]
    regionList      = [ "CR", "SR" ]
    binEdgeList     = [ "CRBE", "SRBE" ]
    leptonList      = [ "el", "mu", "both" ]
    mvaBinList      = [ "D1", "D2", "D3", "D4", "All" ]

    inputDir                    = "njetRootHistos"
    if not os.path.exists(inputDir):
        print "Input directory",inputDir,"does not exist."

    inputFileDataName           = inputDir+"/Data.root"
    inputFileMcSubtractName     = inputDir+"/Other.root"
    inputFileMcSubtract2Name    = inputDir+"/TT.root"

    outputFile = ROOT.TFile( inputDir+"/Data_QcdOnly.root", "RECREATE")

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
                            
                            dataQcdOnlyHist     = ROOT.TH1D()
                            inputFileData       = ROOT.TFile.Open( inputFileDataName )
                            dataAllHist         = inputFileData.Get( "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin )
                            
                            dataQcdOnlyHist     = copy.deepcopy( dataAllHist.Clone() )
                            inputFileData.Close()
                            
                            inputFileMC         = ROOT.TFile.Open( inputFileMcSubtractName )
                            mcSubtractNameHist  = inputFileMC.Get( "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin ) 
                            dataQcdOnlyHist.Add( mcSubtractNameHist, -1 )
                            inputFileMC.Close()
                            
                            inputFileMC2        = ROOT.TFile.Open( inputFileMcSubtract2Name )
                            mcSubtract2NameHist = inputFileMC2.Get( "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin ) 
                            dataQcdOnlyHist.Add( mcSubtract2NameHist, -1 )
                            inputFileMC2.Close()

                            outputFile.cd()
                            dataQcdOnlyHist.Write()

    outputFile.Close()

if __name__ == '__main__' :
    main()
