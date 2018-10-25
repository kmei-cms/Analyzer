import ROOT
import math
import array
import copy

ROOT.gROOT.SetBatch( True )

def main():

    myGraphs                = []
    myMassPts               = [ "350", "450", "550", "650", "750", "850" ]
#    myMassPts               = [ "550" ]
    myLegendTags            = []
    
    myTreeDirectory         = "MiniTrees_18Oct2018_V14_DeepESM_v2.0"
    TTBkgdFileName          = "TT.root"
    QCDBkgdFileName         = "QCD.root"
    myBkgdEffCutOffs        = [ i/1000.0 + 0.01 for i in range(190) ]
    
    totBKGDHisto        = ROOT.TH1F( "h_njets_total", "h_njets_total", 8, 7, 15 )
    totBKGDHistoNorm    = ROOT.TH1F( "h_njets_total_norm", "h_njets_total_norm", 8, 7, 15 )
        
    myTTBkgdFile        = ROOT.TFile( myTreeDirectory+"/"+TTBkgdFileName )
    myTTBkgdMiniTree    = myTTBkgdFile.Get( "myMiniTree" )
    myTTBkgdDESM        = []
    myTTBkgdNJets       = []
    
    myTTBkgdDESM_7      = []
    myTTBkgdDESM_8      = []
    myTTBkgdDESM_9      = []
    myTTBkgdDESM_10     = []
    myTTBkgdDESM_11     = []
    myTTBkgdDESM_12     = []
    myTTBkgdDESM_13     = []
    myTTBkgdDESM_14     = []
    
    myTTBkgdEventWeight = 0.0
    
    nEntries            = myTTBkgdMiniTree.GetEntries()
    
    for myEntryIt in xrange( 0, nEntries ):
        myTTBkgdMiniTree.GetEntry( myEntryIt )
    
        if myEntryIt is 0 :
            myTTBkgdEventWeight = myTTBkgdMiniTree.Weight * 80000.
    
        myTTBkgdDESM.append( myTTBkgdMiniTree.deepESM_val )
        myTTBkgdNJets.append( myTTBkgdMiniTree.NGoodJets_pt30 )
        totBKGDHisto.Fill( myTTBkgdMiniTree.NGoodJets_pt30 , myTTBkgdEventWeight )
    
        if myTTBkgdMiniTree.NGoodJets_pt30 is 7:
            myTTBkgdDESM_7.append( myTTBkgdMiniTree.deepESM_val )
        
        elif myTTBkgdMiniTree.NGoodJets_pt30 is 8:
            myTTBkgdDESM_8.append( myTTBkgdMiniTree.deepESM_val )
        
        elif myTTBkgdMiniTree.NGoodJets_pt30 is 9:
            myTTBkgdDESM_9.append( myTTBkgdMiniTree.deepESM_val )
    
        elif myTTBkgdMiniTree.NGoodJets_pt30 is 10:
            myTTBkgdDESM_10.append( myTTBkgdMiniTree.deepESM_val )
        
        elif myTTBkgdMiniTree.NGoodJets_pt30 is 11:
            myTTBkgdDESM_11.append( myTTBkgdMiniTree.deepESM_val )
        
        elif myTTBkgdMiniTree.NGoodJets_pt30 is 12:
            myTTBkgdDESM_12.append( myTTBkgdMiniTree.deepESM_val )
        
        elif myTTBkgdMiniTree.NGoodJets_pt30 is 13:
            myTTBkgdDESM_13.append( myTTBkgdMiniTree.deepESM_val )
        
        elif myTTBkgdMiniTree.NGoodJets_pt30 > 14:
            myTTBkgdDESM_14.append( myTTBkgdMiniTree.deepESM_val )
    
    myTTBkgdFile.Close()
    
    myQCDBkgdFile        = ROOT.TFile( myTreeDirectory+"/"+QCDBkgdFileName )
    myQCDBkgdMiniTree    = myQCDBkgdFile.Get( "myMiniTree" )
    myQCDBkgdDESM        = []
    myQCDBkgdNJets       = []

    myQCDBkgdDESM_7      = []
    myQCDBkgdDESM_8      = []
    myQCDBkgdDESM_9      = []
    myQCDBkgdDESM_10     = []
    myQCDBkgdDESM_11     = []
    myQCDBkgdDESM_12     = []
    myQCDBkgdDESM_13     = []
    myQCDBkgdDESM_14     = []
    
    myQCDBkgdEventWeight = 0.0
    
    nEntries            = myQCDBkgdMiniTree.GetEntries()
    
    for myEntryIt in xrange( 0, nEntries ):
        myQCDBkgdMiniTree.GetEntry( myEntryIt )
    
        if myEntryIt is 0 :
            myQCDBkgdEventWeight = myQCDBkgdMiniTree.Weight * 80000.
    
        myQCDBkgdDESM.append( myQCDBkgdMiniTree.deepESM_val )
        myQCDBkgdNJets.append( myQCDBkgdMiniTree.NGoodJets_pt30 )
        totBKGDHisto.Fill( myQCDBkgdMiniTree.NGoodJets_pt30 , myQCDBkgdEventWeight )
    
        if myQCDBkgdMiniTree.NGoodJets_pt30 is 7:
            myQCDBkgdDESM_7.append( myQCDBkgdMiniTree.deepESM_val )
        
        elif myQCDBkgdMiniTree.NGoodJets_pt30 is 8:
            myQCDBkgdDESM_8.append( myQCDBkgdMiniTree.deepESM_val )
        
        elif myQCDBkgdMiniTree.NGoodJets_pt30 is 9:
            myQCDBkgdDESM_9.append( myQCDBkgdMiniTree.deepESM_val )
    
        elif myQCDBkgdMiniTree.NGoodJets_pt30 is 10:
            myQCDBkgdDESM_10.append( myQCDBkgdMiniTree.deepESM_val )
        
        elif myQCDBkgdMiniTree.NGoodJets_pt30 is 11:
            myQCDBkgdDESM_11.append( myQCDBkgdMiniTree.deepESM_val )
        
        elif myQCDBkgdMiniTree.NGoodJets_pt30 is 12:
            myQCDBkgdDESM_12.append( myQCDBkgdMiniTree.deepESM_val )
        
        elif myQCDBkgdMiniTree.NGoodJets_pt30 is 13:
            myQCDBkgdDESM_13.append( myQCDBkgdMiniTree.deepESM_val )
        
        elif myQCDBkgdMiniTree.NGoodJets_pt30 > 14:
            myQCDBkgdDESM_14.append( myQCDBkgdMiniTree.deepESM_val )
    
    myQCDBkgdFile.Close()
   
    totWeightedTTEvents_7 = float( len(myTTBkgdDESM_7) ) * myTTBkgdEventWeight
    totWeightedTTEvents_8 = float( len(myTTBkgdDESM_8) ) * myTTBkgdEventWeight
    totWeightedTTEvents_9 = float( len(myTTBkgdDESM_9) ) * myTTBkgdEventWeight
    totWeightedTTEvents_10 = float( len(myTTBkgdDESM_10) ) * myTTBkgdEventWeight
    totWeightedTTEvents_11 = float( len(myTTBkgdDESM_11) ) * myTTBkgdEventWeight
    totWeightedTTEvents_12 = float( len(myTTBkgdDESM_12) ) * myTTBkgdEventWeight
    totWeightedTTEvents_13 = float( len(myTTBkgdDESM_13) ) * myTTBkgdEventWeight
    totWeightedTTEvents_14 = float( len(myTTBkgdDESM_14) ) * myTTBkgdEventWeight
    
    totWeightedQCDEvents_7 = float( len(myQCDBkgdDESM_7) ) * myQCDBkgdEventWeight
    totWeightedQCDEvents_8 = float( len(myQCDBkgdDESM_8) ) * myQCDBkgdEventWeight
    totWeightedQCDEvents_9 = float( len(myQCDBkgdDESM_9) ) * myQCDBkgdEventWeight
    totWeightedQCDEvents_10 = float( len(myQCDBkgdDESM_10) ) * myQCDBkgdEventWeight
    totWeightedQCDEvents_11 = float( len(myQCDBkgdDESM_11) ) * myQCDBkgdEventWeight
    totWeightedQCDEvents_12 = float( len(myQCDBkgdDESM_12) ) * myQCDBkgdEventWeight
    totWeightedQCDEvents_13 = float( len(myQCDBkgdDESM_13) ) * myQCDBkgdEventWeight
    totWeightedQCDEvents_14 = float( len(myQCDBkgdDESM_14) ) * myQCDBkgdEventWeight

    totWeightedBKGDEvents_7 = totWeightedTTEvents_7 + totWeightedQCDEvents_7
    totWeightedBKGDEvents_8 = totWeightedTTEvents_8 + totWeightedQCDEvents_8
    totWeightedBKGDEvents_9 = totWeightedTTEvents_9 + totWeightedQCDEvents_9
    totWeightedBKGDEvents_10 = totWeightedTTEvents_10 + totWeightedQCDEvents_10
    totWeightedBKGDEvents_11 = totWeightedTTEvents_11 + totWeightedQCDEvents_11
    totWeightedBKGDEvents_12 = totWeightedTTEvents_12 + totWeightedQCDEvents_12
    totWeightedBKGDEvents_13 = totWeightedTTEvents_13 + totWeightedQCDEvents_13
    totWeightedBKGDEvents_14 = totWeightedTTEvents_14 + totWeightedQCDEvents_14

    totWeightedBKGDEvents = totWeightedBKGDEvents_7 + totWeightedBKGDEvents_8 + totWeightedBKGDEvents_9 + totWeightedBKGDEvents_10 + totWeightedBKGDEvents_11 + totWeightedBKGDEvents_12 + totWeightedBKGDEvents_13 + totWeightedBKGDEvents_14

    totWeightedBKGDRatio_7 = totWeightedBKGDEvents_7 / totWeightedBKGDEvents
    totWeightedBKGDRatio_8 = totWeightedBKGDEvents_8 / totWeightedBKGDEvents
    totWeightedBKGDRatio_9 = totWeightedBKGDEvents_9 / totWeightedBKGDEvents
    totWeightedBKGDRatio_10 = totWeightedBKGDEvents_10 / totWeightedBKGDEvents
    totWeightedBKGDRatio_11 = totWeightedBKGDEvents_11 / totWeightedBKGDEvents
    totWeightedBKGDRatio_12 = totWeightedBKGDEvents_12 / totWeightedBKGDEvents
    totWeightedBKGDRatio_13 = totWeightedBKGDEvents_13 / totWeightedBKGDEvents
    totWeightedBKGDRatio_14 = totWeightedBKGDEvents_14 / totWeightedBKGDEvents
    
    for myBkgdEntry in myTTBkgdNJets :
        totBKGDHistoNorm.Fill( myBkgdEntry, myTTBkgdEventWeight/totWeightedBKGDEvents )
    for myBkgdEntry in myQCDBkgdNJets :
        totBKGDHistoNorm.Fill( myBkgdEntry, myQCDBkgdEventWeight/totWeightedBKGDEvents )

    totBKGDHistoNorm.SetLineColor(ROOT.kRed)

    canv0 = ROOT.TCanvas( "canv", "histograms", 0, 0, 600, 600) 
    canv0.SetLogy()
    pad0  = ROOT.TPad( "pad", "pad", 0, 0, 1.0, 1.0 )
    pad0.SetLogy()
    pad0.Draw()
    totBKGDHisto.Draw()
    canv0.Update()
    canv0.SaveAs("Histograms/tot_bkgd_histo.pdf")

    totBKGDHistoNorm.Draw()
    canv0.Update()
    canv0.SaveAs("Histograms/tot_bkgd_histo_norm.pdf")

    del canv0, totBKGDHisto 
    myBin1Histos        = []   
    myBin2Histos        = []


    print "MassPt,BkgdEffCutOff,NJets,Nsig,NTTbkgd,NQCDbkgd"
    graphMax = 0.0

    for myMassPt in myMassPts :
    
        #Define all the directories and locations of files
        mySignalFileName    = "rpv_stop_"+myMassPt+".root"
        myLegendTags.append( decideLegendTag( myMassPt ) )
       
        #Get the DeepESM branch and distribution
        mySignalFile        = ROOT.TFile( myTreeDirectory+"/"+mySignalFileName )
        mySignalMiniTree    = mySignalFile.Get( "myMiniTree" )
        mySignalDESM        = []
        mySignalNJets       = []

        mySignalDESM_7      = []
        mySignalDESM_8      = []
        mySignalDESM_9      = []
        mySignalDESM_10     = []
        mySignalDESM_11     = []
        mySignalDESM_12     = []
        mySignalDESM_13     = []
        mySignalDESM_14     = []
        
        mySignalEventWeight = 0.0
        
        nEntries            = mySignalMiniTree.GetEntries()
        
        for myEntryIt in xrange( 0, nEntries ):
            mySignalMiniTree.GetEntry( myEntryIt )
        
            if myEntryIt is 0 :
                mySignalEventWeight = mySignalMiniTree.Weight * 80000.0
        
            mySignalDESM.append( mySignalMiniTree.deepESM_val )
            mySignalNJets.append( mySignalMiniTree.NGoodJets_pt30 )
        
            if mySignalMiniTree.NGoodJets_pt30 is 7:
                mySignalDESM_7.append( mySignalMiniTree.deepESM_val )
            
            elif mySignalMiniTree.NGoodJets_pt30 is 8:
                mySignalDESM_8.append( mySignalMiniTree.deepESM_val )
            
            elif mySignalMiniTree.NGoodJets_pt30 is 9:
                mySignalDESM_9.append( mySignalMiniTree.deepESM_val )
        
            elif mySignalMiniTree.NGoodJets_pt30 is 10:
                mySignalDESM_10.append( mySignalMiniTree.deepESM_val )
            
            elif mySignalMiniTree.NGoodJets_pt30 is 11:
                mySignalDESM_11.append( mySignalMiniTree.deepESM_val )
            
            elif mySignalMiniTree.NGoodJets_pt30 is 12:
                mySignalDESM_12.append( mySignalMiniTree.deepESM_val )
            
            elif mySignalMiniTree.NGoodJets_pt30 is 13:
                mySignalDESM_13.append( mySignalMiniTree.deepESM_val )
            
            elif mySignalMiniTree.NGoodJets_pt30 > 14:
                mySignalDESM_14.append( mySignalMiniTree.deepESM_val )
        
        mySignalFile.Close()
        
        #Calculate the significance
      
        myXArray      = array.array('d')
        myYArray      = array.array('d')
        
    #    for i in range(9) :
    #       myXArray.append(i)
    #       myYArray.append(-1*i)
        myDESMCut_7  = []
        myDESMCut_8  = []
        myDESMCut_9  = []
        myDESMCut_10 = []
        myDESMCut_11 = []
        myDESMCut_12 = []
        myDESMCut_13 = []
        myDESMCut_14 = []
       
        myTempDESMCut_7 = 0.999            
        myTempDESMCut_8 = 0.999            
        myTempDESMCut_9 = 0.999            
        myTempDESMCut_10 = 0.999            
        myTempDESMCut_11 = 0.999            
        myTempDESMCut_12 = 0.999            
        myTempDESMCut_13 = 0.999            
        myTempDESMCut_14 = 0.999            

        for myBkgdEffCutOff in myBkgdEffCutOffs:
           
            if 'myNJetsHisto' in locals() :
               del myNJetsHisto
               del myNJetsHistoBin1
            myHistoName  = "h_njets_"+str(myBkgdEffCutOff)+"_mass"+str(myMassPt)
            myHistoName1 = myHistoName+"_bin1"
            myHistoName2 = myHistoName+"_bin2"
            myNJetsHistoBin1 = ROOT.TH1F( myHistoName1, myHistoName1, 8, 7, 15 )
            myNJetsHisto     = ROOT.TH1F( myHistoName2, myHistoName2, 8, 7, 15 )
            
            foundMyCut_7    = False
            foundMyCut_8    = False
            foundMyCut_9    = False
            foundMyCut_10   = False
            foundMyCut_11   = False
            foundMyCut_12   = False
            foundMyCut_13   = False
            foundMyCut_14   = False

            while( not foundMyCut_7 ) :
                myBackgroundEff_7 = ( float( sum( i > myTempDESMCut_7 for i in myTTBkgdDESM_7 ) )*myTTBkgdEventWeight + float( sum( j > myTempDESMCut_7 for j in myQCDBkgdDESM_7 ) )* myQCDBkgdEventWeight )/ float( totWeightedTTEvents_7 + totWeightedQCDEvents_7 )
                if myBackgroundEff_7 > myBkgdEffCutOff :
                    myDESMCut_7.append( myTempDESMCut_7 )
                    foundMyCut_7 = True
                else :
                    myTempDESMCut_7 = myTempDESMCut_7 - .001
            
            while( not foundMyCut_8 ) :
                myBackgroundEff_8 = ( float( sum( i > myTempDESMCut_8 for i in myTTBkgdDESM_8 ) )*myTTBkgdEventWeight + float( sum( j > myTempDESMCut_8 for j in myQCDBkgdDESM_8 ) )* myQCDBkgdEventWeight )/ float( totWeightedTTEvents_8 + totWeightedQCDEvents_8 )
                if myBackgroundEff_8 > myBkgdEffCutOff :
                    myDESMCut_8.append( myTempDESMCut_8 )
                    foundMyCut_8 = True
                else :
                    myTempDESMCut_8 = myTempDESMCut_8 - .001
            
            while( not foundMyCut_9 ) :
                myBackgroundEff_9 = ( float( sum( i > myTempDESMCut_9 for i in myTTBkgdDESM_9 ) )*myTTBkgdEventWeight + float( sum( j > myTempDESMCut_9 for j in myQCDBkgdDESM_9 ) )* myQCDBkgdEventWeight )/ float( totWeightedTTEvents_9 + totWeightedQCDEvents_9 )
                if myBackgroundEff_9 > myBkgdEffCutOff :
                    myDESMCut_9.append( myTempDESMCut_9 )
                    foundMyCut_9 = True
                else :
                    myTempDESMCut_9 = myTempDESMCut_9 - .001
            
            while( not foundMyCut_10 ) :
                myBackgroundEff_10 = ( float( sum( i > myTempDESMCut_10 for i in myTTBkgdDESM_10 ) )*myTTBkgdEventWeight + float( sum( j > myTempDESMCut_10 for j in myQCDBkgdDESM_10 ) )* myQCDBkgdEventWeight )/ float( totWeightedTTEvents_10 + totWeightedQCDEvents_10 )
                if myBackgroundEff_10 > myBkgdEffCutOff :
                    myDESMCut_10.append( myTempDESMCut_10 )
                    foundMyCut_10 = True
                else :
                    myTempDESMCut_10 = myTempDESMCut_10 - .001
            
            while( not foundMyCut_11 ) :
                myBackgroundEff_11 = ( float( sum( i > myTempDESMCut_11 for i in myTTBkgdDESM_11 ) )*myTTBkgdEventWeight + float( sum( j > myTempDESMCut_11 for j in myQCDBkgdDESM_11 ) )* myQCDBkgdEventWeight )/ float( totWeightedTTEvents_11 + totWeightedQCDEvents_11 )
                if myBackgroundEff_11 > myBkgdEffCutOff :
                    myDESMCut_11.append( myTempDESMCut_11 )
                    foundMyCut_11 = True
                else :
                    myTempDESMCut_11 = myTempDESMCut_11 - .001
            
            while( not foundMyCut_12 ) :
                myBackgroundEff_12 = ( float( sum( i > myTempDESMCut_12 for i in myTTBkgdDESM_12 ) )*myTTBkgdEventWeight + float( sum( j > myTempDESMCut_12 for j in myQCDBkgdDESM_12 ) )* myQCDBkgdEventWeight )/ float( totWeightedTTEvents_12 + totWeightedQCDEvents_12 )
                if myBackgroundEff_12 > myBkgdEffCutOff :
                    myDESMCut_12.append( myTempDESMCut_12 )
                    foundMyCut_12 = True
                else :
                    myTempDESMCut_12 = myTempDESMCut_12 - .001
            
            while( not foundMyCut_13 ) :
                myBackgroundEff_13 = ( float( sum( i > myTempDESMCut_13 for i in myTTBkgdDESM_13 ) )*myTTBkgdEventWeight + float( sum( j > myTempDESMCut_13 for j in myQCDBkgdDESM_13 ) )* myQCDBkgdEventWeight )/ float( totWeightedTTEvents_13 + totWeightedQCDEvents_13 )
                if myBackgroundEff_13 > myBkgdEffCutOff :
                    myDESMCut_13.append( myTempDESMCut_13 )
                    foundMyCut_13 = True
                else :
                    myTempDESMCut_13 = myTempDESMCut_13 - .001
            
            while( not foundMyCut_14 ) :
                myBackgroundEff_14 = ( float( sum( i > myTempDESMCut_14 for i in myTTBkgdDESM_14 ) )*myTTBkgdEventWeight + float( sum( j > myTempDESMCut_14 for j in myQCDBkgdDESM_14 ) )* myQCDBkgdEventWeight )/ float( totWeightedTTEvents_14 + totWeightedQCDEvents_14 )
                if myBackgroundEff_14 > myBkgdEffCutOff :
                    myDESMCut_14.append( myTempDESMCut_14 )
                    foundMyCut_14 = True
                else :
                    myTempDESMCut_14 = myTempDESMCut_14 - .001
            
            #Make Histograms Here
            for i in myTTBkgdDESM_7 :
                if i > myDESMCut_7[-1] :
                    myNJetsHisto.Fill( 7, myTTBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 7, myTTBkgdEventWeight )
            
            for i in myTTBkgdDESM_8 :
                if i > myDESMCut_8[-1] :
                    myNJetsHisto.Fill( 8, myTTBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 8, myTTBkgdEventWeight )
            
            for i in myTTBkgdDESM_9 :
                if i > myDESMCut_9[-1] :
                    myNJetsHisto.Fill( 9, myTTBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 9, myTTBkgdEventWeight )
            
            for i in myTTBkgdDESM_10 :
                if i > myDESMCut_10[-1] :
                    myNJetsHisto.Fill( 10, myTTBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 10, myTTBkgdEventWeight )
            
            for i in myTTBkgdDESM_11 :
                if i > myDESMCut_11[-1] :
                    myNJetsHisto.Fill( 11, myTTBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 11, myTTBkgdEventWeight )
            
            for i in myTTBkgdDESM_12 :
                if i > myDESMCut_12[-1] :
                    myNJetsHisto.Fill( 12, myTTBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 12, myTTBkgdEventWeight )
            
            for i in myTTBkgdDESM_13 :
                if i > myDESMCut_13[-1] :
                    myNJetsHisto.Fill( 13, myTTBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 13, myTTBkgdEventWeight )
            
            for i in myTTBkgdDESM_14 :
                if i > myDESMCut_14[-1] :
                    myNJetsHisto.Fill( 14, myTTBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 14, myTTBkgdEventWeight )
            
            for i in myQCDBkgdDESM_7 :
                if i > myDESMCut_7[-1] :
                    myNJetsHisto.Fill( 7, myQCDBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 7, myQCDBkgdEventWeight )
            
            for i in myQCDBkgdDESM_8 :
                if i > myDESMCut_8[-1] :
                    myNJetsHisto.Fill( 8, myQCDBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 8, myQCDBkgdEventWeight )
            
            for i in myQCDBkgdDESM_9 :
                if i > myDESMCut_9[-1] :
                    myNJetsHisto.Fill( 9, myQCDBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 9, myQCDBkgdEventWeight )
            
            for i in myQCDBkgdDESM_10 :
                if i > myDESMCut_10[-1] :
                    myNJetsHisto.Fill( 10, myQCDBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 10, myQCDBkgdEventWeight )
            
            for i in myQCDBkgdDESM_11 :
                if i > myDESMCut_11[-1] :
                    myNJetsHisto.Fill( 11, myQCDBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 11, myQCDBkgdEventWeight )
            
            for i in myQCDBkgdDESM_12 :
                if i > myDESMCut_12[-1] :
                    myNJetsHisto.Fill( 12, myQCDBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 12, myQCDBkgdEventWeight )
            
            for i in myQCDBkgdDESM_13 :
                if i > myDESMCut_13[-1] :
                    myNJetsHisto.Fill( 13, myQCDBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 13, myQCDBkgdEventWeight )
            
            for i in myQCDBkgdDESM_14 :
                if i > myDESMCut_14[-1] :
                    myNJetsHisto.Fill( 14, myQCDBkgdEventWeight )
                else :
                    myNJetsHistoBin1.Fill( 14, myQCDBkgdEventWeight )

            myBin1Histos.append( copy.deepcopy(myNJetsHistoBin1) )
            myBin2Histos.append( copy.deepcopy(myNJetsHisto) )

            leg1 = ROOT.TLegend( 0.7, 0.7, 0.95, 0.88 )

            canv1 = ROOT.TCanvas( "canv", "histograms", 0, 0, 600, 600) 
            canv1.SetLogy()
            pad1  = ROOT.TPad( "pad", "pad", 0, 0, 1.0, 1.0 )
            pad1.SetLogy()
            pad1.SetLeftMargin(0.20)
            pad1.Draw()
            myNJetsHisto.Scale( 1.0/myNJetsHisto.Integral() )
            myNJetsHisto.GetXaxis().SetTitle('Number of Jets')
            #myNJetsHisto.GetYaxis().SetTitle('Normalized Number of Events')
            myNJetsHisto.SetTitle('NJets Distribution with mass '+str(myMassPt)+' and bkgd eff. '+str(myBkgdEffCutOff))
            myNJetsHisto.SetLineColor(ROOT.kGreen)
            myNJetsHisto.SetMinimum( 0.00001 )
            myNJetsHisto.SetMaximum( 1.0 )
            myNJetsHisto.SetStats( 0 )
            myNJetsHisto.Draw()
            leg1.AddEntry(myNJetsHisto,'Bin 2','l')
            myNJetsHistoBin1.Scale( 1.0/myNJetsHistoBin1.Integral() )
            myNJetsHistoBin1.SetLineColor(ROOT.kMagenta)
            leg1.AddEntry(myNJetsHistoBin1,'Bin 1','l')
            myNJetsHistoBin1.SetStats( 0 )
            myNJetsHistoBin1.Draw('SAME')
            leg1.AddEntry(totBKGDHistoNorm,'No DeepESM Cut','l')
            totBKGDHistoNorm.SetStats( 0 ) 
            totBKGDHistoNorm.Draw('SAME')
            leg1.Draw()
            canv1.Update()
            canv1.SaveAs("Histograms/cut"+str(myBkgdEffCutOff)+"_m"+str(myMassPt)+"_bkgd_histo.pdf")

            print("%s,%.2f,%i,%.2f,%.2f,%.2f" % ( myMassPt, myBkgdEffCutOff, 7, sum( i > myDESMCut_7[-1] for i in mySignalDESM_7)*mySignalEventWeight, sum( i > myDESMCut_7[-1] for i in myTTBkgdDESM_7)*myTTBkgdEventWeight, sum(i > myDESMCut_7[-1] for i in myQCDBkgdDESM_7)*myQCDBkgdEventWeight ) )
            print("%s,%.2f,%i,%.2f,%.2f,%.2f" % ( myMassPt, myBkgdEffCutOff, 8, sum( i > myDESMCut_8[-1] for i in mySignalDESM_8)*mySignalEventWeight, sum( i > myDESMCut_8[-1] for i in myTTBkgdDESM_8)*myTTBkgdEventWeight, sum(i > myDESMCut_8[-1] for i in myQCDBkgdDESM_8)*myQCDBkgdEventWeight ) )
            print("%s,%.2f,%i,%.2f,%.2f,%.2f" % ( myMassPt, myBkgdEffCutOff, 9, sum( i > myDESMCut_9[-1] for i in mySignalDESM_9)*mySignalEventWeight, sum( i > myDESMCut_9[-1] for i in myTTBkgdDESM_9)*myTTBkgdEventWeight, sum(i > myDESMCut_9[-1] for i in myQCDBkgdDESM_9)*myQCDBkgdEventWeight ) )
            print("%s,%.2f,%i,%.2f,%.2f,%.2f" % ( myMassPt, myBkgdEffCutOff, 10, sum( i > myDESMCut_10[-1] for i in mySignalDESM_10)*mySignalEventWeight, sum( i > myDESMCut_10[-1] for i in myTTBkgdDESM_10)*myTTBkgdEventWeight, sum(i > myDESMCut_10[-1] for i in myQCDBkgdDESM_10)*myQCDBkgdEventWeight ) )
            print("%s,%.2f,%i,%.2f,%.2f,%.2f" % ( myMassPt, myBkgdEffCutOff, 11, sum( i > myDESMCut_11[-1] for i in mySignalDESM_11)*mySignalEventWeight, sum( i > myDESMCut_11[-1] for i in myTTBkgdDESM_11)*myTTBkgdEventWeight, sum(i > myDESMCut_11[-1] for i in myQCDBkgdDESM_11)*myQCDBkgdEventWeight ) )
            print("%s,%.2f,%i,%.2f,%.2f,%.2f" % ( myMassPt, myBkgdEffCutOff, 12, sum( i > myDESMCut_12[-1] for i in mySignalDESM_12)*mySignalEventWeight, sum( i > myDESMCut_12[-1] for i in myTTBkgdDESM_12)*myTTBkgdEventWeight, sum(i > myDESMCut_12[-1] for i in myQCDBkgdDESM_12)*myQCDBkgdEventWeight ) )
            print("%s,%.2f,%i,%.2f,%.2f,%.2f" % ( myMassPt, myBkgdEffCutOff, 13, sum( i > myDESMCut_13[-1] for i in mySignalDESM_13)*mySignalEventWeight, sum( i > myDESMCut_13[-1] for i in myTTBkgdDESM_13)*myTTBkgdEventWeight, sum(i > myDESMCut_13[-1] for i in myQCDBkgdDESM_13)*myQCDBkgdEventWeight ) )
            print("%s,%.2f,%i,%.2f,%.2f,%.2f" % ( myMassPt, myBkgdEffCutOff, 14, sum( i > myDESMCut_14[-1] for i in mySignalDESM_14)*mySignalEventWeight, sum( i > myDESMCut_14[-1] for i in myTTBkgdDESM_14)*myTTBkgdEventWeight, sum(i > myDESMCut_14[-1] for i in myQCDBkgdDESM_14)*myQCDBkgdEventWeight ) )

            sig7 = calcPerBinSignificance( myDESMCut_7[-1], mySignalEventWeight, mySignalDESM_7, myTTBkgdEventWeight, myTTBkgdDESM_7, myQCDBkgdEventWeight, myQCDBkgdDESM_7 )
            sig8 = calcPerBinSignificance( myDESMCut_8[-1], mySignalEventWeight, mySignalDESM_8, myTTBkgdEventWeight, myTTBkgdDESM_8, myQCDBkgdEventWeight, myQCDBkgdDESM_8 )
            sig9 = calcPerBinSignificance( myDESMCut_9[-1], mySignalEventWeight, mySignalDESM_9, myTTBkgdEventWeight, myTTBkgdDESM_9, myQCDBkgdEventWeight, myQCDBkgdDESM_9 )
            sig10 = calcPerBinSignificance( myDESMCut_10[-1], mySignalEventWeight, mySignalDESM_10, myTTBkgdEventWeight, myTTBkgdDESM_10, myQCDBkgdEventWeight, myQCDBkgdDESM_10 )
            sig11 = calcPerBinSignificance( myDESMCut_11[-1], mySignalEventWeight, mySignalDESM_11, myTTBkgdEventWeight, myTTBkgdDESM_11, myQCDBkgdEventWeight, myQCDBkgdDESM_11 )
            sig12 = calcPerBinSignificance( myDESMCut_12[-1], mySignalEventWeight, mySignalDESM_12, myTTBkgdEventWeight, myTTBkgdDESM_12, myQCDBkgdEventWeight, myQCDBkgdDESM_12 )
            sig13 = calcPerBinSignificance( myDESMCut_13[-1], mySignalEventWeight, mySignalDESM_13, myTTBkgdEventWeight, myTTBkgdDESM_13, myQCDBkgdEventWeight, myQCDBkgdDESM_13 )
            sig14 = calcPerBinSignificance( myDESMCut_14[-1], mySignalEventWeight, mySignalDESM_14, myTTBkgdEventWeight, myTTBkgdDESM_14, myQCDBkgdEventWeight, myQCDBkgdDESM_14 )
        
            #print sig7, sig8, sig9, sig10, sig11, sig12

            myTotSigValue = float( calcTotalSignificance( sig7, sig8, sig9, sig10, sig11, sig12, sig13, sig14 ) )
            if myTotSigValue > graphMax :
                graphMax = myTotSigValue
            myXArray.append( float(myBkgdEffCutOff) )
            myYArray.append( myTotSigValue )

        graphMax = 1.1* graphMax
        myGraph = ROOT.TGraph(len(myXArray), myXArray, myYArray)
        myGraph.SetTitle( 'Significance' )
        myGraph.GetXaxis().SetTitle( 'Bkgd Eff (#epsilon)' )
        myGraph.GetYaxis().SetTitle( 'Significance' )
        myGraph.SetLineColor( decideColor( myMassPt ) )
        myGraph.SetLineWidth( 2 )
        myGraphs.append( copy.deepcopy(myGraph) )
    
    
    cs1     = ROOT.TCanvas( "cs1", "cs1", 700, 700 )
    pad1    = ROOT.TPad( "pad1", "pad1", 0, 0, 1.0, 1.0 )
    pad1.SetLeftMargin(0.15)
    pad1.Draw()
    leg1    = ROOT.TLegend( 0.7, 0.7, 0.95, 0.95 )
    pad1.cd()
    
    count   = 0

    for graph in myGraphs :
        
        if count is 0 :
            graph.SetMinimum( 0.0 )
            graph.SetMaximum( graphMax )
            graph.Draw()
        
        else :
            graph.Draw('SAME')
        
        leg1.AddEntry( graph, myLegendTags[count], 'l' )
        count += 1
        
    leg1.Draw()
    cs1.SaveAs("temp_eff_all_NEW_FAST.pdf")

def decideColor( filename ) :
    if "350" in filename :
        return ROOT.kRed
    elif "450" in filename :
        return ROOT.kOrange
    elif "550" in filename :
        return ROOT.kGreen
    elif "650" in filename :
        return ROOT.kCyan
    elif "750" in filename :
        return ROOT.kBlue
    elif "850" in filename :
        return ROOT.kViolet

def decideLegendTag( filename ) :
    if "350" in filename :
        return "m_{S} = 350 GeV"
    elif "450" in filename :
        return "m_{S} = 450 GeV"
    elif "550" in filename :
        return "m_{S} = 550 GeV"
    elif "650" in filename :
        return "m_{S} = 650 GeV"
    elif "750" in filename :
        return "m_{S} = 750 GeV"
    elif "850" in filename :
        return "m_{S} = 850 GeV"

def calcPerBinSignificance( myCut, mySignalWeight, mySignalList, myTTBackgroundWeight, myTTBackgroundList, myQCDBackgroundWeight, myQCDBackgroundList ) :
    numSigEvents = float(sum( i > myCut for i in mySignalList ) )* mySignalWeight 
    numTTBkgdEvents = float( sum( i > myCut for i in myTTBackgroundList ) ) * myTTBackgroundWeight
    numQCDBkgdEvents = float( sum( i > myCut for i in myQCDBackgroundList ) ) * myQCDBackgroundWeight
    significance = numSigEvents / ( numTTBkgdEvents + numQCDBkgdEvents + numQCDBkgdEvents*numQCDBkgdEvents + ( 0.30 * numTTBkgdEvents ) * ( 0.30 * numTTBkgdEvents ) )
    return significance

def calcTotalSignificance( sig7, sig8, sig9, sig10, sig11, sig12, sig13, sig14 ) :
    return math.sqrt( sig7*sig7 + sig8*sig8 + sig9*sig9 + sig10*sig10 + sig11*sig11 + sig12*sig12 + sig13*sig13 + sig14*sig14 )

def calcTotalSignificanceNo12( sig7, sig8, sig9, sig10, sig11, sig12 ) :
    return math.sqrt( sig7*sig7 + sig8*sig8 + sig9*sig9 + sig10*sig10 + sig11*sig11 )


if __name__=='__main__':
    main()
