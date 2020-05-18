#define AnalyzeSignalKM_cxx
#include "Analyzer/Analyzer/include/AnalyzeSignalKM.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <TFile.h>

AnalyzeSignalKM::AnalyzeSignalKM()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeSignalKM::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    std::vector<std::string> nomWeightList      { "noWghts", "allWghts" };
    std::vector<std::string> pdfWeightList      { "noPdf", "pdfUp", "pdfDown" };
    std::vector<std::string> sclWeightList      { "noScl", "sclUp", "sclDown" };
    std::vector<std::string> isrWeightList      { "noIsr", "isrUp", "isrDown" };
    std::vector<std::string> fsrWeightList      { "noFsr", "fsrUp", "fsrDown" };

    //Define 1D histograms
    for( std::string nomWeight : nomWeightList ) {
        for( std::string pdfWeight : pdfWeightList ) {
            for( std::string sclWeight : sclWeightList ) {
                for( std::string isrWeight : isrWeightList ) {
                    for( std::string fsrWeight : fsrWeightList ) {
                        my_histos.emplace("h_njets_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_njets_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_njets_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 22, 3.5, 25.5 ) );
                        my_histos.emplace("h_njets_12binIncl_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_njets_12binIncl_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_njets_12binIncl_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 6, 7, 13) );
    
                        my_histos.emplace("h_ht_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_ht_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_ht_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0.0, 3000.0 ) );
                        my_histos.emplace("h_ht_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_ht_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_ht_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0.0, 3000.0 ) );
                        my_histos.emplace("h_ht_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_ht_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_ht_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0.0, 3000.0 ) );
                        
                        my_histos.emplace("h_jetpt_1_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_1_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_1_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_1_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_1_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_1_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_1_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_1_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_1_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );

                        my_histos.emplace("h_jetpt_2_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_2_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_2_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_2_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_2_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_2_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_2_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_2_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_2_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                    
                        my_histos.emplace("h_jetpt_3_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_3_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_3_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_3_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_3_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_3_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_3_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_3_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_3_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                    
                        my_histos.emplace("h_jetpt_4_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_4_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_4_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_4_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_4_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_4_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_4_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_4_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_4_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                    
                        my_histos.emplace("h_jetpt_5_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_5_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_5_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_5_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_5_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_5_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_5_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_5_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_5_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                    
                        my_histos.emplace("h_jetpt_6_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_6_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_6_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_6_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_6_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_6_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_6_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_6_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_6_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                    
                        my_histos.emplace("h_jetpt_7_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_7_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_7_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_7_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_7_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_7_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_7_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_7_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_7_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                    
                    
                        my_histos.emplace("h_nb_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_nb_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_nb_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 7, -0.5, 6 ) );
                        my_histos.emplace("h_nb_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_nb_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_nb_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 7, -0.5, 6 ) );
                        my_histos.emplace("h_nb_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_nb_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_nb_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 7, -0.5, 6 ) );
                    
                        my_histos.emplace("h_mva_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_mva_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_mva_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 15, 0.0, 1.0 ) );
                        my_histos.emplace("h_mva_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_mva_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_mva_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 15, 0.0, 1.0 ) );
                        my_histos.emplace("h_mva_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_mva_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_mva_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 15, 0.0, 1.0 ) );
                        
                        my_histos.emplace("h_nvtx_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_nvtx_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_nvtx_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0, 60 ) );
                        my_histos.emplace("h_nvtx_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_nvtx_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_nvtx_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0, 60 ) );
                        my_histos.emplace("h_nvtx_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_nvtx_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_nvtx_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0, 60 ) );
                    
                        my_histos.emplace("h_met_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_met_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_met_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0, 3000.0 ) );
                        my_histos.emplace("h_met_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_met_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_met_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0, 3000.0 ) );
                        my_histos.emplace("h_met_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_met_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_met_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0, 3000.0 ) );
                        
                        my_histos.emplace("h_mbl_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_mbl_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_mbl_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 20, 200 ) );
                        my_histos.emplace("h_mbl_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_mbl_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_mbl_ge4j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 20, 200 ) );
                        my_histos.emplace("h_mbl_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_mbl_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_mbl_ge7j_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 20, 200 ) );
                    }
                }
            }
        }
    }
}

//Put everything you want to do per event here.
void AnalyzeSignalKM::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );
        if( !isQuiet ) std::cout<<weight<<std::endl;

        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        const auto& runtype         = tr.getVar<std::string>("runtype");     
        const auto& NGoodBJets_pt30 = tr.getVar<int>("NGoodBJets_pt30");
        const auto& Mbl             = tr.getVar<double>("Mbl");
        
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30");
        const auto& NGoodJets_pt30  = tr.getVar<int>("NGoodJets_pt30");
        
        const auto& nVtx            = tr.getVar<int>("NVtx");
        const auto& passMadHT       = tr.getVar<bool>("passMadHT");
        const auto& passBaseline    = tr.getVar<bool>("passBaselineGoodOffline1l");
        const auto& met             = tr.getVar<double>("MET");

        const auto& passTrigger     = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC   = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodMuons      = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons  = tr.getVar<int>("NGoodElectrons");

        const auto& Jet_pt_1        = tr.getVar<double>("Jet_pt_1");
        const auto& Jet_pt_2        = tr.getVar<double>("Jet_pt_2");
        const auto& Jet_pt_3        = tr.getVar<double>("Jet_pt_3");
        const auto& Jet_pt_4        = tr.getVar<double>("Jet_pt_4");
        const auto& Jet_pt_5        = tr.getVar<double>("Jet_pt_5");
        const auto& Jet_pt_6        = tr.getVar<double>("Jet_pt_6");
        const auto& Jet_pt_7        = tr.getVar<double>("Jet_pt_7");

        const auto& deepESM_val     = tr.getVar<double>("deepESM_val");
       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------
        // -- Define weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
        // ------------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double htDerivedScaleFactor = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;

        double sclUpFactor          = 1.0;
        double sclDownFactor        = 1.0;
        double pdfUpFactor          = 1.0;
        double pdfDownFactor        = 1.0;
        double isrUpFactor          = 1.0;
        double isrDownFactor        = 1.0;
        double fsrUpFactor          = 1.0;
        double fsrDownFactor        = 1.0;

        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            weight = eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;

            const auto& sclUp           = tr.getVar<double>("scaleWeightUp");
            const auto& sclDown         = tr.getVar<double>("scaleWeightDown");
            const auto& pdfUp           = tr.getVar<double>("PDFweightUp");
            const auto& pdfDown         = tr.getVar<double>("PDFweightDown");
            const auto& isrUp           = tr.getVar<double>("PSweight_ISRUp");
            const auto& isrDown         = tr.getVar<double>("PSweight_ISRDown");
            const auto& fsrUp           = tr.getVar<double>("PSweight_FSRUp");
            const auto& fsrDown         = tr.getVar<double>("PSweight_FSRDown");

            sclUpFactor                 = sclUp;
            sclDownFactor               = sclDown;
            pdfUpFactor                 = pdfUp;
            pdfDownFactor               = pdfDown;
            isrUpFactor                 = isrUp;
            isrDownFactor               = isrDown;
            fsrUpFactor                 = fsrUp;
            fsrDownFactor               = fsrDown;
        }
    
        std::vector<std::string> pdfWeightList      { "noPdf", "pdfUp", "pdfDown" };
        std::vector<std::string> sclWeightList      { "noScl", "sclUp", "sclDown" };
        std::vector<std::string> isrWeightList      { "noIsr", "isrUp", "isrDown" };
        std::vector<std::string> fsrWeightList      { "noFsr", "fsrUp", "fsrDown" };

        std::map<std::string, double> extraWeightDict  = {};
        double extraPdfWeight                          = 1.0;
        double extraSclWeight                          = 1.0;
        double extraIsrWeight                          = 1.0;
        double extraFsrWeight                          = 1.0;

        for( std::string pdfWeight : pdfWeightList ) {
           
            extraPdfWeight                                 = 1.0;
            if( pdfWeight == "pdfUp" )      extraPdfWeight = pdfUpFactor;
            if( pdfWeight == "pdfDown" )    extraPdfWeight = pdfDownFactor;
            
            for( std::string sclWeight : sclWeightList ) {

                extraSclWeight                                     = 1.0;
                if( sclWeight == "sclUp" )          extraSclWeight = sclUpFactor;
                if( sclWeight == "sclDown" )        extraSclWeight = sclDownFactor;
                
                for( std::string isrWeight : isrWeightList ) {
                    
                    extraIsrWeight                                     = 1.0;
                    if( isrWeight == "isrUp" )          extraIsrWeight = isrUpFactor;
                    if( isrWeight == "isrDown" )        extraIsrWeight = isrDownFactor;
                    
                    for( std::string fsrWeight : fsrWeightList ) {
                
                        extraFsrWeight                                     = 1.0;
                        if( fsrWeight == "fsrUp" )          extraFsrWeight = fsrUpFactor;
                        if( fsrWeight == "fsrDown" )        extraFsrWeight = fsrDownFactor;
                        
                        extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight] = extraPdfWeight*extraSclWeight*extraIsrWeight*extraFsrWeight; 

                        //std::cout<<pdfWeight<<" "<<sclWeight<<" "<<isrWeight<<" "<<fsrWeight<<" "<<extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]<<std::endl;
                    }
                }
            }
        }
        
        //Make cuts and fill histograms here
        if( passBaseline && ( ( NGoodElectrons == 0 && NGoodMuons == 1 ) || ( NGoodElectrons == 1 && NGoodMuons == 0 ) ) && passTrigger && passTriggerMC ) {
            for( std::string pdfWeight : pdfWeightList ) {
                for( std::string sclWeight : sclWeightList ) {
                    for( std::string isrWeight : isrWeightList ) {
                        for( std::string fsrWeight : fsrWeightList ) {
                            
                            my_histos["h_njets_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodJets_pt30, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_ht_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_trigger_pt30, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_1_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_2_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_3_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_4_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_5_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_6_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_7_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_nb_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodBJets_pt30, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_mva_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_met_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( met, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_nvtx_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_mbl_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Mbl, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            
                            if( NGoodJets_pt30 >= 4 ) {
                                my_histos["h_ht_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_trigger_pt30, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_1_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_2_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_3_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_4_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_5_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_6_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_7_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_nb_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodBJets_pt30, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_mva_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_met_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( met, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_nvtx_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_mbl_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Mbl, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                
                                if( NGoodJets_pt30 >= 7 ) {
                                    if( NGoodJets_pt30 < 12.5 ) {
                                        my_histos["h_njets_12binIncl_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodJets_pt30, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]); }
                                    else {
                                        my_histos["h_njets_12binIncl_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( 12, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]); }
                                    
                                    my_histos["h_ht_ge7j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_trigger_pt30, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_1_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_2_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_3_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_4_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_5_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_6_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_7_ge4j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_nb_ge7j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodBJets_pt30, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_mva_ge7j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_met_ge7j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( met, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_nvtx_ge7j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_mbl_ge7j_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Mbl, weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                }
                            }
                            
                            my_histos["h_njets_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodJets_pt30, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_ht_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_trigger_pt30, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_1_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_2_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_3_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_4_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_5_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_6_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_jetpt_7_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_nb_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodBJets_pt30, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_mva_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_met_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( met, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_nvtx_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            my_histos["h_mbl_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Mbl, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                            
                            if( NGoodJets_pt30 >= 4 ) {

                                my_histos["h_ht_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_trigger_pt30, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_1_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_2_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_3_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_4_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_5_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_6_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_jetpt_7_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_nb_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodBJets_pt30, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_mva_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_met_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( met, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_nvtx_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                my_histos["h_mbl_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Mbl, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                
                                if( NGoodJets_pt30 >= 7 ) {
                                    if( NGoodJets_pt30 < 12.5 ) {
                                        my_histos["h_njets_12binIncl_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodJets_pt30, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]); }
                                    else {
                                        my_histos["h_njets_12binIncl_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( 12, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]); }
                                    
                                    my_histos["h_ht_ge7j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_trigger_pt30, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_1_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_2_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_3_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_4_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_5_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_6_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_jetpt_7_ge4j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_nb_ge7j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodBJets_pt30, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_mva_ge7j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_met_ge7j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( met, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_nvtx_ge7j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                    my_histos["h_mbl_ge7j_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Mbl, eventweight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]);
                                }
                            }
                        }
                    }
                }
            }
        }
    } 
}

void AnalyzeSignalKM::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }

    for (const auto &p : my_3d_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }
    
}
