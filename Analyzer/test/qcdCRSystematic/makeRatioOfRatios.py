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
                    default = False, help = 'Make all plots norm on the y scale' )

(options,args) = parser.parse_args()

def main():

    ROOT.TH1.SetDefaultSumw2()

    yearTag = "2016"
    if options.year :
        yearTag = "2017"

    #Define the new ratio of ratios.
    njetNameArray       = [ "7", "8", "9", "10", "11" ]
    histoNameArray      = []

    histoDict           = {}

    #Name the input file after what you hadded all the root files in mvaRatioPlots.
    inputFileName                       = "ratios_"+yearTag+".root"
    inputFile                         = ROOT.TFile.Open( inputFileName )
        
    for njetName in njetNameArray :
        myTTHisto                       = inputFile.Get( "h_ratio_TT_"+njetName )
        myDataHisto                     = inputFile.Get( "h_ratio_Data_SingleMuon_"+njetName )
        myRatioOfRatioHisto             = myDataHisto.Clone("h_ratio_DataCR_over_TTSR_"+njetName )
        myRatioOfRatioHisto.Divide( myRatioOfRatioHisto, myTTHisto, 1, 1, "" )
        histoDict[njetName]             = copy.deepcopy( myRatioOfRatioHisto )
    
    inputFile.Close()

    #Write all histograms to output file
    outputFileName = "ratioOfRatios_"+yearTag+".root"
    outputFile  = ROOT.TFile.Open( outputFileName, "RECREATE" )
    for njetName in njetNameArray :
        histoDict[njetName].Write()
    outputFile.Close()

    #Check output pdf directory
    outputDir = "ratioPdfs"
    
    if not os.path.exists(outputDir) :
        print "Output directory does not exist, so making output directory", outputDir
        os.makedirs(outputDir)

    plotAllPlots( yearTag, outputDir, histoDict, njetNameArray )
    plotIndividualPlots( yearTag, outputDir, histoDict, njetNameArray )

def plotAllPlots( yearTag, outputDir, histoDict, njetNameArray ) :

    c1              = ROOT.TCanvas( "c1", "c1", 0, 0, 1000, 1000 )
    p1              = ROOT.TPad( "p1", "p1", 0, 0, 1.0, 1.0 )
    p1.SetLeftMargin( 0.15 )
    if( options.log ):
        p1.SetLogy()
    p1.Draw()
    c1.Update()

    p1.cd()
    leg1            = ROOT.TLegend( 0.30, 0.78, 0.89, 0.88 )
    leg1.SetFillStyle( 0 )
    leg1.SetBorderSize( 0 )
    leg1.SetLineWidth( 1 )
    leg1.SetNColumns( 2 )
    leg1.SetTextFont( 42 )
    leg1.SetTextSize( 0.035 )
    
    for njetName in njetNameArray :
           
        h1 = histoDict[njetName]
        h1.SetStats(0)
        h1.SetLineWidth(2)
        h1.SetLineColor(decideColor(njetName))
            
        h1.SetTitle( "Data CR / Signal MC Ratios Per NJet" )
        h1.GetXaxis().SetTitleSize( 0.04 )
        h1.GetYaxis().SetTitleSize( 0.04 )
        h1.GetXaxis().SetLabelSize( 0.04 )
        h1.GetYaxis().SetLabelSize( 0.04 )
        h1.GetYaxis().SetTitle("")
        h1.GetYaxis().SetTitleOffset( 1.2 )
        h1.SetMaximum( 2.0 )
        h1.SetMinimum( 0.0 )

        leg1.AddEntry(h1, "N_{j} = "+njetName,"l")
        if njetName == "7" :   
            h1.Draw("HIST")
        else :
            h1.Draw("HIST SAME" )


    leg1.Draw()
    c1.SaveAs("ratioPdfs/"+yearTag+"_ratioOfRatioComp.pdf")
    c1.Clear()

def plotIndividualPlots( yearTag, histoDict, njetNameArray ) :

    c1              = ROOT.TCanvas( "c1", "c1", 0, 0, 1000, 1000 )
    p1              = ROOT.TPad( "p1", "p1", 0, 0, 1.0, 1.0 )
    p1.SetLeftMargin( 0.15 )
    if( options.log ):
        p1.SetLogy()
    p1.Draw()
    c1.Update()

    p1.cd()
    for njetName in njetNameArray :
        leg1            = ROOT.TLegend( 0.30, 0.78, 0.89, 0.88 )
        leg1.SetFillStyle( 0 )
        leg1.SetBorderSize( 0 )
        leg1.SetLineWidth( 1 )
        leg1.SetNColumns( 2 )
        leg1.SetTextFont( 42 )
        leg1.SetTextSize( 0.035 )
           
        h1 = histoDict[njetName]
        h1.SetStats(0)
        h1.SetLineWidth(2)
        h1.SetLineColor(decideColor(njetName))
            
        h1.SetTitle( "Data CR / Signal TT MC Ratio for N_{j} = "+njetName )
        h1.GetXaxis().SetTitleSize( 0.04 )
        h1.GetYaxis().SetTitleSize( 0.04 )
        h1.GetXaxis().SetLabelSize( 0.04 )
        h1.GetYaxis().SetLabelSize( 0.04 )
        h1.GetYaxis().SetTitle("")
        h1.GetYaxis().SetTitleOffset( 1.2 )
        h1.SetMaximum( 2.0 )
        h1.SetMinimum( 0.0 )

        leg1.AddEntry(h1, "N_{j} = "+njetName,"l")
        h1.Draw("HIST E")

        leg1.Draw()
        c1.SaveAs(outputDir+"/"+yearTag+"_"+njetName+"j_ratioOfRatioComp.pdf")
        c1.Clear()

def decideColor( njetName ) :
    if njetName == "7" :
        #return ROOT.kRed
        return ROOT.kBlue
    if njetName == "8" :
        return ROOT.kRed
    if njetName == "9" :
        #return ROOT.kRed
        return ROOT.kGreen+2
    if njetName == "10" :
        return ROOT.kMagenta
    if njetName == "11" :
        return ROOT.kOrange+2

if __name__ == '__main__' :
    main()
