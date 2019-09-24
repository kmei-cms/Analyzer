#This script takes the output of the 100 toys and calculating the final systematic 

#!/bin/python
import ROOT
import copy
import os.path
import array
import math
import numpy as np

ROOT.gROOT.SetBatch( True )

from optparse import OptionParser

parser  = OptionParser()

parser.add_option('--year', action='store_true',
                    dest='year',
                    default = False, help = 'Make all plots norm on the y scale' )

(options,args) = parser.parse_args()

def main():

    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()

    yearTag = "2016"
    if options.year :
        yearTag = "2017"

    mvaBinArray         = [ "D1", "D2", "D3", "D4" ]
    njetNameArray       = [ "7", "8", "9", "10", "11", "12" ]
    nToys               = 100

    #Get the ttbar nominal shape systematic from the njets_for_Aron.root file by taking the actual njets shape and doing the division
    ttFileName                          = yearTag+"_njets_for_Aron.root"
    ttFile                              = ROOT.TFile.Open( ttFileName )

    ttHistoAll                          = ROOT.TH1D()
    ttHisto                             = {}

    #Make output directories for plots
    if not os.path.exists("toyPdfs"):
        print "Making output directory toyPdfs"
        os.makedirs("toyPdfs")
    
    if not os.path.exists("njetPdfs"):
        print "Making output directory njetPdfs"
        os.makedirs("njetPdfs")

    for mva in mvaBinArray:
        tempHisto                       = ttFile.Get( mva+"_TT_h_njets_pt30_1l" )
        ttHisto[mva]                    = copy.deepcopy( tempHisto.Clone( mva+"_TT" ) )

        ttHisto[mva].Scale( 1.0 / ttHisto[mva].Integral() )
        
        if mva == "D1" :
            ttHistoAll                  = copy.deepcopy( tempHisto.Clone( "All_TT" ) )
        else :
            ttHistoAll.Add( tempHisto )

    ttHistoAll.Scale( 1.0 / ttHistoAll.Integral() )

    for mva in mvaBinArray :
        ttHisto[mva].Divide( ttHisto[mva], ttHistoAll, 1, 1, "" )
    
    ttFile.Close()

    #Calculate the error based on all the toys - input file is the output of throwToysForSystematic.py
    toyFileName                         = yearTag+"_evtWghtSyst_AllToys.root"
    toyFile                             = ROOT.TFile.Open( toyFileName )

    toyHisto                            = {} #Holds all the histograms for each toy

    contHisto                           = {}
    errHisto                            = {}

    for mva in mvaBinArray :
        contHisto[mva]                  = {}
        errHisto[mva]                   = {}

        for njet in njetNameArray :
            contHisto[mva][njet]        = []
            errHisto[mva][njet]         = []


    #Everything had to be rebinned because we were using the shifted histograms from njets_for_Aron.root (0 to 6)file and it would not divide properly with the not shifted histograms from throwToysForSystematic.py (7 to 12) - hence the seemingly redundant histograms. This can be fixed, but was done in a bind and worked.
    for i in range( nToys ) :
        toyHisto[i]                     = {}
        toyHisto[i]["All"]              = ROOT.TH1D( "h_All_shift_toy"+str(i), "h_All_shift_toy"+str(i), 6, 0, 6 )
        for mva in mvaBinArray :
            toyHisto[i][mva]            = ROOT.TH1D( "h_"+mva+"_shift_toy"+str(i), "h_"+mva+"_shift_toy"+str(i), 6, 0, 6 )

        for mva in mvaBinArray :
            tempHisto                   = toyFile.Get( "h_"+mva+"_toy"+str(i) )
            for binIter in range( 0, tempHisto.GetNbinsX() + 1 ) :
                toyHisto[i][mva].SetBinContent( binIter, tempHisto.GetBinContent( binIter ) )
                toyHisto[i][mva].SetBinError( binIter, tempHisto.GetBinError( binIter ) )

        for mva in mvaBinArray :
            toyHisto[i]["All"].Add( toyHisto[i][mva] )
            toyHisto[i][mva].Scale( 1.0 / toyHisto[i][mva].Integral() )

        toyHisto[i]["All"].Scale( 1.0 / toyHisto[i]["All"].Integral() )

        for mva in mvaBinArray :
            #Divide each njets shape for each MVA bin by thte njets shape for all the MVA bins
            toyHisto[i][mva].Divide( toyHisto[i][mva], toyHisto[i]["All"], 1, 1, "" )
            #Divide the resulting ratio with the normalized shape from ttbar nominal (like all systematics)
            toyHisto[i][mva].Divide( toyHisto[i][mva], ttHisto[mva], 1, 1, "" )

    
    for mva in mvaBinArray :
        for njet in njetNameArray :
            for i in range( nToys ) :
                contHisto[mva][njet].append( toyHisto[i][mva].GetBinContent( int(njet) - 6 ) )
                errHisto[mva][njet].append( toyHisto[i][mva].GetBinError( int(njet) - 6 ) )
   
    #The next few lines of code were just added to visualize the spread of the value for number of events per MVA bin and per njet bin (may not be useful to see anymore)
    contHistoArray                      = {}
    errHistoArray                       = {}

    for mva in mvaBinArray :
        contHistoArray[mva]             = {}
        errHistoArray[mva]              = {}
        for njet in njetNameArray :
            contHistoArray[mva][njet]   = ROOT.TH1D( "h_"+mva+"_"+njet+"j", "h_"+mva+"_"+njet+"j", 100, np.mean( contHisto[mva][njet] ) - 3*np.std( contHisto[mva]["12"] ), np.mean( contHisto[mva][njet] ) + 3*np.std( contHisto[mva]["12"] ) )
            errHistoArray[mva][njet]    = ROOT.TH1D( "h_"+mva+"_"+njet+"j_err", "h_"+mva+"_"+njet+"j_err", 300, np.mean( errHisto[mva][njet] ) - .5, np.mean( errHisto[mva][njet] ) + .5 )

    for mva in mvaBinArray :
        for njet in njetNameArray :
            for i in range( len(contHisto[mva][njet]) ) :
                contHistoArray[mva][njet].Fill( contHisto[mva][njet][i] )
                errHistoArray[mva][njet].Fill( errHisto[mva][njet][i] )

    outputFileName                  = yearTag+"_toys_StatErrPlusFullDev.root"
    outputFile                      = ROOT.TFile( outputFileName, "RECREATE" ) 
    for mva in mvaBinArray:
        qcdCRHisto                  = ROOT.TH1D( mva+"_qcdCR", mva+"_qcdCR", 6, 0, 6 )
        qcdCRErrHisto               = ROOT.TH1D( mva+"_qcdCRErr", mva+"_qcdCRErr", 6, 0, 6 )
        
        for binIter in xrange( 1, 7 ):
            errValue                = np.std( contHisto[mva][str(binIter+6)] ) #Base error is the standard deviation for each MVA bin each njet bin for the 100 toys
            
            #errValueAdded           = 0.0 #No full deviation from 1
            errValueAdded           = ( np.mean( contHisto[mva][str(binIter+6)] ) - 1 ) #full deviation from one - will be squared in any meaningful calculation below

            #print mva,binIter+6,np.mean( errHisto[mva][str(binIter+6)] ), np.std( contHisto[mva][str(binIter+6)] ) #Just used to compare the size of the average error in each toy with the standard deviation of the toy values
            qcdCRHisto.SetBinContent( binIter, np.mean( contHisto[mva][str(binIter+6)] ) )
            qcdCRHisto.SetBinError( binIter, math.sqrt( errValue**2 + errValueAdded**2 ) )
            qcdCRErrHisto.SetBinContent( binIter, 1 + math.sqrt( errValue**2 + errValueAdded**2 )/ qcdCRHisto.GetBinContent( binIter ) )
            
            #contHistoArray[mva][str(binIter+6)].Write()
            #errHistoArray[mva][str(binIter+6)].Write()

        qcdCRHisto.Write()
        qcdCRErrHisto.Write()

        #Plot error histograms
        c1              = ROOT.TCanvas( "c1", "c1", 0, 0, 1000, 1000 )
        p1              = ROOT.TPad( "p1", "p1", 0, 0, 1.0, 1.0 )
        p1.SetLeftMargin( 0.15 )
        p1.SetLogy(False)
        p1.Draw()
        p1.SetGridy(1)
        c1.Update()
    
        p1.cd()
        leg1            = ROOT.TLegend( 0.30, 0.78, 0.89, 0.88 )
        leg1.SetFillStyle( 0 )
        leg1.SetBorderSize( 0 )
        leg1.SetLineWidth( 1 )
        leg1.SetNColumns( 2 )
        leg1.SetTextFont( 42 )
        leg1.SetTextSize( 0.035 )
           
        qcdCRHisto.SetStats(0)
        qcdCRHisto.SetLineWidth(2)
        qcdCRHisto.SetLineColor(decideColorMVA(mva))
            
        qcdCRHisto.SetTitle( "QCD CR Nominal R Values "+mva+" in "+yearTag )
        qcdCRHisto.GetXaxis().SetTitleSize( 0.04 )
        qcdCRHisto.GetYaxis().SetTitleSize( 0.04 )
        qcdCRHisto.GetXaxis().SetLabelSize( 0.04 )
        qcdCRHisto.GetYaxis().SetLabelSize( 0.04 )
        qcdCRHisto.GetYaxis().SetTitle("")
        qcdCRHisto.GetYaxis().SetTitleOffset( 1.2 )
        qcdCRHisto.GetYaxis().SetRangeUser( 0.5, 1.5 )

        qcdCRHisto.Draw("HIST E")
        c1.SaveAs("njetPdfs/"+yearTag+"_ratio_"+mva+".pdf")
       
        for mva in mvaBinArray :
            for njet in njetNameArray :
                c2              = ROOT.TCanvas( "c2", "c2", 0, 0, 1000, 1000 )
                p2              = ROOT.TPad( "p2", "p2", 0, 0, 1.0, 1.0 )
                p2.SetLeftMargin( 0.15 )
                p2.SetLogy(False)
                p2.Draw()
                p2.SetGridy(1)
                c2.Update()
            
                p2.cd()
                leg1            = ROOT.TLegend( 0.30, 0.78, 0.89, 0.88 )
                leg1.SetFillStyle( 0 )
                leg1.SetBorderSize( 0 )
                leg1.SetLineWidth( 1 )
                leg1.SetNColumns( 2 )
                leg1.SetTextFont( 42 )
                leg1.SetTextSize( 0.035 )
        
                tempHisto       = contHistoArray[mva][njet]
                tempHisto.SetLineWidth(2)
                tempHisto.SetLineColor(decideColorMVA(mva))
                    
                tempHisto.SetTitle( "Nominal R Values for 100 Toys in "+mva+" for "+njet+"j in "+yearTag )
                tempHisto.GetXaxis().SetTitleSize( 0.04 )
                tempHisto.GetYaxis().SetTitleSize( 0.04 )
                tempHisto.GetXaxis().SetLabelSize( 0.04 )
                tempHisto.GetYaxis().SetLabelSize( 0.04 )
                tempHisto.GetYaxis().SetTitle("")
                tempHisto.GetYaxis().SetTitleOffset( 1.2 )
        
                tempHisto.Draw("HIST E")
                c2.SaveAs("toyPdfs/"+yearTag+"_toyR_"+mva+"_"+njet+"j.pdf")
                
                tempHisto2       = errHistoArray[mva][njet]
                tempHisto2.SetLineWidth(2)
                tempHisto2.SetLineColor(decideColorMVA(mva))
                    
                tempHisto2.SetTitle( "Nominal Error Bar Values for 100 Toys in "+mva+" for "+njet+"j in "+yearTag )
                tempHisto2.GetXaxis().SetTitleSize( 0.04 )
                tempHisto2.GetYaxis().SetTitleSize( 0.04 )
                tempHisto2.GetXaxis().SetLabelSize( 0.04 )
                tempHisto2.GetYaxis().SetLabelSize( 0.04 )
                tempHisto2.GetYaxis().SetTitle("")
                tempHisto2.GetYaxis().SetTitleOffset( 1.2 )
        
                tempHisto2.Draw("HIST E")
                c2.SaveAs("toyPdfs/"+yearTag+"_toyRprime_"+mva+"_"+njet+"j.pdf")

def decideColorMVA( mva ) :
    if mva == "D1" :
        return ROOT.kBlue
    if mva == "D2" :
        return ROOT.kRed
    if mva == "D3" :
        return ROOT.kGreen+2
    if mva == "D4" :
        return ROOT.kMagenta

if __name__ == '__main__' :
    main()
