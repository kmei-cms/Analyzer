# This script plots all the MVA distribution histograms per njet and then creates a file with the ratios calculated as (MVA for all events with 7 and greater jets) / (MVA for njets = 7,8,9,10,11) 
# Last edited by Kelvin Mei on September 19th, 2019

#!/bin/python
import ROOT
import copy
import os.path
import array

ROOT.gROOT.SetBatch( True )

from optparse import OptionParser

parser  = OptionParser()

parser.add_option('--log', action='store_true',
                    dest='log',
                    default = False, help = 'Make all plots log on the y scale' )

parser.add_option('--year', action='store_true',
                    dest='year',
                    default = False, help = 'Run over 2017 instead of 2016' )

(options,args) = parser.parse_args()

def main():

    yearTag = "2016"
    
    if options.year :
        yearTag = "2017"

    inputDir            = "mvaPlots" #default output directory for MakeMVAHistograms.py

    if not os.path.exists(inputDir):
        print "Input directory does not exist. Please run makeMVAHistograms.py first"

    outputDir           = "mvaRatioPlots"
    if not os.path.exists(outputDir):
        print "Output directory does not exist, so making output directory", outputDir
        os.makedirs(outputDir)

    categoryNameArray   = [ "TT", "QCD", "Data_SingleMuon" ]
    njetNameArray       = [ "7", "8", "9", "10", "11" ]

    histoDict           = {}
    
    #Create dictionary with all the histograms to plot later
    for categoryName in categoryNameArray:
        histoDict[categoryName]             = {}
        inputFileName                       = inputDir+"/"+categoryName+"_ratio_"+yearTag+".root"
        myInputFile                         = ROOT.TFile.Open( inputFileName )
        myHisto                             = myInputFile.Get( "h_mva_all" )
        histoDict[categoryName]["All"]      = copy.deepcopy( myHisto )

        for njetName in njetNameArray :
            myHisto                         = myInputFile.Get( "h_mva_"+njetName+"j" )
            histoDict[categoryName][njetName]    = copy.deepcopy( myHisto )
    
    #This function plots everything and writes it to an output file
    plotRatioPlots( yearTag, outputDir, histoDict, categoryNameArray, njetNameArray, "All" )

def plotRatioPlots( yearTag, outputDir, histoDict, categoryNameArray, njetNameArray, numeratorName ) :

    for categoryName in categoryNameArray:
        
        for njetName in njetNameArray :

            #Define canvas and pad
            c1              = ROOT.TCanvas( "c1", "c1", 0, 0, 700, 700 )
            p1              = ROOT.TPad( "p1", "p1", 0, 0.41, 1.0, 1.0 )
            p1.SetLeftMargin( 0.15 )
            p1.SetBottomMargin( 0 )
            if( options.log ):
                p1.SetLogy()
            p1.Draw()
            p2              = ROOT.TPad( "p2", "p2", 0, 0, 1.0, 0.4 )
            p2.SetTopMargin( 0 )
            p2.SetLeftMargin( 0.15 )
            p2.SetBottomMargin( 0.3 )
            p2.Draw()
            c1.Update()
        
            p1.cd()

            #Create legend
            leg1            = ROOT.TLegend( 0.62, 0.70, 0.89, 0.88 )
            leg1.SetFillStyle( 0 )
            leg1.SetBorderSize( 0 )
            leg1.SetLineWidth( 1 )
            leg1.SetNColumns( 1 )
            leg1.SetTextFont( 42 )
            leg1.SetTextSize( 0.055 )
               
            h_numerator = histoDict[categoryName][numeratorName]
            h_numerator.SetStats(0)
            h_numerator.SetLineWidth(2)
            h_numerator.SetLineColor(ROOT.kRed)
            
            if h_numerator.Integral() != 0 :
                h_numerator.Scale( 1.0 / h_numerator.Integral() )
            else :
                print "Numerator histogram has no entries and cannot be normalized"
            
            h_njets = histoDict[categoryName][njetName]
            h_njets.SetStats(0)
            h_njets.SetLineWidth(2)
            h_njets.SetLineColor(ROOT.kBlue)
                
            if h_njets.Integral() != 0 :
                h_njets.Scale( 1.0 / h_njets.Integral() )
            else :
                print "Denominator histogram has no entries and cannot be normalized"
            
            h_numerator.SetTitle( "MVA Shape for "+categoryName+" for N_{j} = "+njetName )
                
            leg1.AddEntry(h_numerator, "All N_{j}","l")
            leg1.AddEntry(h_njets, "N_{j} = "+njetName,"l")
               
            h_numerator.Draw("HIST")
            h_njets.Draw("SAME E")
    
            leg1.Draw()
    
            h_ratio = h_numerator.Clone( "h_ratio_"+categoryName+"_"+njetName )
            h_ratio.SetLineColor(ROOT.kBlue)
            h_ratio.SetTitle("")
            h_ratio.GetXaxis().SetTickLength( 0.1 )
            h_ratio.GetYaxis().SetNdivisions( 12, 5, 0 )

            h_ratio.GetXaxis().SetTitle("MVA Score")
            h_ratio.GetYaxis().SetTitle("Ratio: All NJets / NJets = "+njetName ) 
            h_ratio.GetXaxis().SetTitleSize(0.075*.8)
            h_ratio.GetYaxis().SetTitleSize(0.075*.8)
            h_ratio.GetYaxis().SetTitleOffset( 0.9 )
            h_ratio.GetXaxis().SetLabelSize(0.08*.8)
            h_ratio.GetYaxis().SetLabelSize(0.08*.8)
            h_ratio.Divide( h_ratio, h_njets ) 
            h_ratio.SetMinimum(0.0)
            h_ratio.SetMaximum(1.55)
            h_ratio.SetStats(ROOT.kFALSE)
            
            p2.cd()
            
            ROOT.gPad.SetGridy(1)
            ROOT.gPad.SetTicks( 1,1) 
            h_ratio.Draw("PE")
            
            line = ROOT.TLine( -0.1, 1.0, 1.1, 1.0 )
            line.SetLineWidth( 2 )
            line.SetLineColor( ROOT.kRed )
            line.Draw("SAME")
            c1.Update()
            
            c1.SaveAs(outputDir+"/MVAratio_"+yearTag+"_"+categoryName+"_"+njetName+".pdf")
            c1.Clear()

            outputFileName = outputDir+"/"+yearTag+"_"+categoryName+"_"+njetName+"_ratio.root"
            outputFile     = ROOT.TFile( outputFileName, "RECREATE" )
            h_ratio.Write()
            outputFile.Close()

if __name__ == '__main__' :
    main()
