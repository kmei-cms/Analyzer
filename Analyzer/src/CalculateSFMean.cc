#define CalculateSFMean_cxx
#include "Analyzer/Analyzer/include/CalculateSFMean.h"
#include "Analyzer/Analyzer/include/Histo.h"

#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

CalculateSFMean::CalculateSFMean()
{
}

void CalculateSFMean::InitHistos(std::string filetag)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_Histos.emplace_back(new Histo1D(filetag+"_ht",       100, 0, 5, "htDerivedweightUncor", {"passMadHT"}, {"Lumi", "Weight"}));
    my_Histos.emplace_back(new Histo1D(filetag+"_sclUp",    100, 0, 5, "scaleWeightUpUncor",   {"passMadHT"}, {"Lumi", "Weight"}));
    my_Histos.emplace_back(new Histo1D(filetag+"_sclDown",  100, 0, 5, "scaleWeightDownUncor", {"passMadHT"}, {"Lumi", "Weight"}));
    my_Histos.emplace_back(new Histo1D(filetag+"_pdf_Up",   100, 0, 5, "PDFweightUpUncor",     {"passMadHT"}, {"Lumi", "Weight"}));
    my_Histos.emplace_back(new Histo1D(filetag+"_pdf_Down", 100, 0, 5, "PDFweightDownUncor",   {"passMadHT"}, {"Lumi", "Weight"}));    
}

void CalculateSFMean::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
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
      
void CalculateSFMean::WriteHistos( TFile* outfile ) 
{
    outfile->cd();
    
    for(const auto& h : my_Histos)
    {
        h->Write();
    }
}
