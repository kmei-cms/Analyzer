#define CalculateHtSF_cxx
#include "Analyzer/Analyzer/include/CalculateHtSF.h"
#include "Analyzer/Analyzer/include/Histo.h"

#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

CalculateHtSF::CalculateHtSF()
{
}

void CalculateHtSF::InitHistos(std::string filetag)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_Histos.emplace_back(new Histo1D(filetag, 100, 0, 5, "htDerivedweightUncor", {"passMadHT"}, {"Lumi", "Weight"}));
}

void CalculateHtSF::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    const auto& runtype = tr.getVar<std::string>("runtype");
    const auto& filetag = tr.getVar<std::string>("filetag");

    InitHistos(filetag);

    while( tr.getNextEvent() )
    {
        //------------------------------------
        //-- Print Event Number
        //------------------------------------
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        //-----------------------------------
        //-- Make sure you are running over MC
        //-- For now (soon will run on data)
        //-----------------------------------
        if( runtype != "MC" ) 
        {
            std::cerr<<"Please run over an MC file since these scale factors should not be applied to data!!"<<std::endl;
            break;
        }
               
        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        for(auto& h : my_Histos)
        {
            h->Fill(tr);
        }
    }//END of while tr.getNextEvent loop   
}//END of function
      
void CalculateHtSF::WriteHistos( TFile* outfile ) 
{
    outfile->cd();
    
    for(const auto& h : my_Histos)
    {
        h->Write();
    }
}
