#define AnalyzeEventShape_cxx
#include "Analyzer/Analyzer/include/AnalyzeEventShape.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>


// includes for the event shapes
#include "Framework/Framework/include/bdt_350to650_fwm10_jmtev_top6.h"
#include "Framework/Framework/include/EventShapeVariables.h"
#include "Framework/Framework/src/get_cmframe_jets.c"
//#include "Framework/Framework/include/fisher_350to650_fwm10_jmtev_top6.h"

AnalyzeEventShape::AnalyzeEventShape()
{
    printf("\n\n In AnalyzeEventShape constructor.\n\n") ; fflush( stdout ) ;
    InitHistos();
}

void AnalyzeEventShape::InitHistos()
{

    printf("\n\n In AnalyzeEventShape::InitHistos\n\n" ) ;

    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("h_met", std::make_shared<TH1D>("h_met","h_met", 20, 0, 200));
    my_histos.emplace("h_ht", std::make_shared<TH1D>("h_ht","h_ht", 60, 0, 3000));


}

void AnalyzeEventShape::Loop(NTupleReader& tr, double weight, int maxevents, std::string filetag, bool isQuiet)
{

    printf("\n\n In AnalyzeEventShape::Loop\n\n" ) ; fflush( stdout ) ;

    int total_events = tr.getNEntries() ;

    if ( maxevents > 0 ) {
       printf("\n\n  maxevents set to %d\n", maxevents ) ; fflush( stdout ) ;
    } else {
       printf("\n\n  Will run over all events in tree : %d\n", total_events ) ; fflush( stdout ) ;
    }

    if ( maxevents > 0 ) total_events = maxevents ;

    int modnum = total_events / 100 ;
    if ( total_events > 0 && modnum == 0 ) modnum = 1 ;


    while( tr.getNextEvent() )
    {
        const std::string& runtype = tr.getVar<std::string>("runtype");
        const double& MET     = tr.getVar<double>("MET");
        const double& HT      = tr.getVar<double>("HT");
        double eventweight(1.) ;
        if ( runtype != "Data" ) {
           const double& Weight  = tr.getVar<double>("Weight");
           eventweight = Weight ;
        }

        if ( tr.getEvtNum() % modnum == 0 ) {
           printf("   %9d / %9d  (%.3f)\n", tr.getEvtNum(), total_events, ((1.*tr.getEvtNum())/(1.*total_events)) ) ; fflush( stdout ) ;
        }


        if (maxevents > 0 && tr.getEvtNum() >= maxevents) break;        
        
        my_histos["h_met"]->Fill(MET, eventweight);
        my_histos["h_ht"]->Fill(HT, eventweight);

    } // end of event loop

}

bool AnalyzeEventShape::PassTriggerGeneral(std::vector<std::string> &mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    bool passTrigger = false;
    for(unsigned int i=0; i<TriggerNames.size(); ++i)
    {
        if(TriggerPass.at(i) != 1)
            continue;
        std::string trigname = TriggerNames.at(i);
        if( std::any_of(mytriggers.begin(), mytriggers.end(), [&] (std::string s) { return trigname.find(s)!=std::string::npos; }) )
        {
            passTrigger = true;
            break;
        }
    }
    return passTrigger;

}


bool AnalyzeEventShape::PassTriggerAllHad(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {
        //"HLT_PFHT1050", // 2017 trigger
        //"HLT_PFHT900"
        //"HLT_PFHT380_SixPFJet32_DoublePFBTagCSV", // 2017 trigger
        //"HLT_PFHT430_SixPFJet40_PFBTagCSV", // 2017 trigger
        "HLT_PFHT450_SixJet40_BTagCSV",
            "HLT_PFHT400_SixJet30_DoubleBTagCSV",            
            };
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

bool AnalyzeEventShape::PassTriggerMuon(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {"HLT_IsoMu24","HLT_IsoTkMu24_v"};
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

bool AnalyzeEventShape::PassTriggerElectron(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {"HLT_Ele27_WPTight_Gsf"};
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

void AnalyzeEventShape::WriteHistos()
{
    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }
    
}
