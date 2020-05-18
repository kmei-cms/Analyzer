#define NonIsoMuonBStudy_cxx
#include "Analyzer/Analyzer/include/NonIsoMuonBStudy.h"
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

NonIsoMuonBStudy::NonIsoMuonBStudy()
{
    InitHistos();
}

//Define all your histograms here. 
void NonIsoMuonBStudy::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    std::vector<std::string> njetTags   = {"1", "2", "3", "4", "5", "6", "7" };
    std::vector<std::string> nbTags     = { "0b", "1b", "2b", "ge0b", "ge1b", "ge2b" };
    std::vector<std::string> n7cutTags  = { "n7cutTag", "n5cutTag", "non7cutTag" };
    //Define 1D histograms
   
    //Make Histograms for All Input Variables - Regular
    
    for( std::string nbTag : nbTags ) {
        for( std::string n7cutTag : n7cutTags ) {
            my_histos.emplace( "h_mva_"+n7cutTag+"_"+nbTag, std::make_shared<TH1D>( ( "h_mva_"+n7cutTag+"_"+nbTag ).c_str(), ( "h_mva_"+n7cutTag+"_"+nbTag ).c_str(), 10, 0.0, 1.0 ) );
            my_histos.emplace( "h_nb_"+n7cutTag+"_"+nbTag, std::make_shared<TH1D>( ( "h_nb_"+n7cutTag+"_"+nbTag ).c_str(), ( "h_nb_"+n7cutTag+"_"+nbTag ).c_str(), 6, -0.5, 5.5 ) );
            my_histos.emplace( "h_njets_"+n7cutTag+"_"+nbTag, std::make_shared<TH1D>( ( "h_njets_"+n7cutTag+"_"+nbTag ).c_str(), ( "h_njets_"+n7cutTag+"_"+nbTag ).c_str(), 6, 6.5, 12.5 ) );
            my_histos.emplace( "h_rho_"+n7cutTag+"_"+nbTag, std::make_shared<TH1D>( ( "h_rho_"+n7cutTag+"_"+nbTag ).c_str(), ( "h_rho_"+n7cutTag+"_"+nbTag ).c_str(), 50, 0.0, 100.0 ) );
            my_histos.emplace( "h_ht_"+n7cutTag+"_"+nbTag, std::make_shared<TH1D>( ( "h_ht_"+n7cutTag+"_"+nbTag ).c_str(), ( "h_ht_"+n7cutTag+"_"+nbTag ).c_str(), 60, 0.0, 3000.0 ) );
            my_histos.emplace( "h_leppt_"+n7cutTag+"_"+nbTag, std::make_shared<TH1D>( ( "h_leppt_"+n7cutTag+"_"+nbTag ).c_str(), ( "h_leppt_"+n7cutTag+"_"+nbTag ).c_str(), 50, 0.0, 1000.0 ) );
    
            for( std::string njetTag : njetTags ) {
                my_histos.emplace( "h_jetpt_"+njetTag+"_"+n7cutTag+"_"+nbTag, std::make_shared<TH1D>( ( "h_jetpt_"+njetTag+"_"+n7cutTag+"_"+nbTag ).c_str(), ( "h_jetpt_"+njetTag+"_"+n7cutTag+"_"+nbTag ).c_str(), 50, 0, 1000 ) );
            }
        }
    }
}

//Put everything you want to do per event here.
void NonIsoMuonBStudy::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );
        if( !isQuiet ) { std::cout<<weight<<std::endl; }
        
        const auto& lumi             = tr.getVar<double>("Lumi");
        const auto& runtype          = tr.getVar<std::string>("runtype");     
        //const auto& filetag          = tr.getVar<std::string>("filetag");     
        
        const auto& passBaseline1l_NIM     = tr.getVar<bool>("passBaseline1l_NonIsoMuon");
        const auto& passBaseline1l_NIM_noNJets     = tr.getVar<bool>("passBaseline1l_NonIsoMuon_noNJets");

        double totalWeightNIM       = 1.0;

        if(runtype == "MC")
        {
            const auto& totalEventWeightNIM = tr.getVar<double>("totalEventWeightNIM");
         
            totalWeightNIM = totalEventWeightNIM*lumi;
        }

        const auto& NNonIsoMuonJets_pt30    = tr.getVar<int>( "NNonIsoMuonJets_pt30" );
        const auto& NBNonIsoMuonJets_pt30   = tr.getVar<int>( "NBNonIsoMuonJets_pt30" );
        const auto& deepESM_valNonIsoMuon   = tr.getVar<double>( "deepESM_valNonIsoMuon" );
        const auto& HT_NonIsoMuon_pt30      = tr.getVar<double>( "HT_NonIsoMuon_pt30" );
        const auto& fixedGridRhoFastjetAll  = tr.getVar<double>( "fixedGridRhoFastjetAll" );
        
        const auto& GoodNonIsoMuons         = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodNonIsoMuons" );
        const auto& JetNonIsoMuons_pt_1     = tr.getVar<double>( "JetNonIsoMuons_pt_1" );
        const auto& JetNonIsoMuons_pt_2     = tr.getVar<double>( "JetNonIsoMuons_pt_2" );
        const auto& JetNonIsoMuons_pt_3     = tr.getVar<double>( "JetNonIsoMuons_pt_3" );
        const auto& JetNonIsoMuons_pt_4     = tr.getVar<double>( "JetNonIsoMuons_pt_4" );
        const auto& JetNonIsoMuons_pt_5     = tr.getVar<double>( "JetNonIsoMuons_pt_5" );
        const auto& JetNonIsoMuons_pt_6     = tr.getVar<double>( "JetNonIsoMuons_pt_6" );
        const auto& JetNonIsoMuons_pt_7     = tr.getVar<double>( "JetNonIsoMuons_pt_7" );

        std::vector<std::string> nbTags     = { "0b", "1b", "2b", "ge0b", "ge1b", "ge2b" };
        std::vector<std::string> n7cutTags  = { "n7cutTag", "non7cutTag" };

        const std::map<std::string, bool> cut_map {
            { "n5cutTag_ge0b"           , passBaseline1l_NIM_noNJets && NNonIsoMuonJets_pt30 >= 5 },
            { "n5cutTag_ge1b"           , NBNonIsoMuonJets_pt30 > 0 && passBaseline1l_NIM_noNJets && NNonIsoMuonJets_pt30 >= 5 },
            { "n5cutTag_ge2b"           , NBNonIsoMuonJets_pt30 > 1 && passBaseline1l_NIM_noNJets && NNonIsoMuonJets_pt30 >= 5 },
            { "n5cutTag_0b"             , NBNonIsoMuonJets_pt30 == 0 && passBaseline1l_NIM_noNJets && NNonIsoMuonJets_pt30 >= 5 },
            { "n5cutTag_1b"             , NBNonIsoMuonJets_pt30 == 1 && passBaseline1l_NIM_noNJets && NNonIsoMuonJets_pt30 >= 5 },
            { "n5cutTag_2b"             , NBNonIsoMuonJets_pt30 == 2 && passBaseline1l_NIM_noNJets && NNonIsoMuonJets_pt30 >= 5 },
            { "n7cutTag_ge0b"             , passBaseline1l_NIM },
            { "n7cutTag_0b"               , NBNonIsoMuonJets_pt30 == 0 && passBaseline1l_NIM },
            { "n7cutTag_1b"               , NBNonIsoMuonJets_pt30 == 1 && passBaseline1l_NIM },
            { "n7cutTag_2b"               , NBNonIsoMuonJets_pt30 == 2 && passBaseline1l_NIM },
            { "n7cutTag_ge1b"             , NBNonIsoMuonJets_pt30 > 0 && passBaseline1l_NIM },
            { "n7cutTag_ge2b"             , NBNonIsoMuonJets_pt30 > 1 && passBaseline1l_NIM },
            { "non7cutTag_ge0b"           , passBaseline1l_NIM_noNJets },
            { "non7cutTag_0b"             , NBNonIsoMuonJets_pt30 == 0 && passBaseline1l_NIM_noNJets },
            { "non7cutTag_1b"             , NBNonIsoMuonJets_pt30 == 1 && passBaseline1l_NIM_noNJets },
            { "non7cutTag_2b"             , NBNonIsoMuonJets_pt30 == 2 && passBaseline1l_NIM_noNJets },
            { "non7cutTag_ge1b"           , NBNonIsoMuonJets_pt30 > 0 && passBaseline1l_NIM_noNJets },
            { "non7cutTag_ge2b"           , NBNonIsoMuonJets_pt30 > 1 && passBaseline1l_NIM_noNJets }
        };

        for( auto & kv : cut_map ) {
            if( kv.second ) {
                my_histos["h_njets_"+kv.first]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                my_histos["h_mva_"+kv.first]->Fill( deepESM_valNonIsoMuon, totalWeightNIM );
                my_histos["h_nb_"+kv.first]->Fill( NBNonIsoMuonJets_pt30, totalWeightNIM );
                my_histos["h_ht_"+kv.first]->Fill( HT_NonIsoMuon_pt30, totalWeightNIM );
                my_histos["h_rho_"+kv.first]->Fill( fixedGridRhoFastjetAll, totalWeightNIM );
                my_histos["h_leppt_"+kv.first]->Fill( GoodNonIsoMuons.at(0).second.Pt(), totalWeightNIM );
                my_histos["h_jetpt_1_"+kv.first]->Fill( JetNonIsoMuons_pt_1, totalWeightNIM );
                my_histos["h_jetpt_2_"+kv.first]->Fill( JetNonIsoMuons_pt_2, totalWeightNIM );
                my_histos["h_jetpt_3_"+kv.first]->Fill( JetNonIsoMuons_pt_3, totalWeightNIM );
                my_histos["h_jetpt_4_"+kv.first]->Fill( JetNonIsoMuons_pt_4, totalWeightNIM );
                my_histos["h_jetpt_5_"+kv.first]->Fill( JetNonIsoMuons_pt_5, totalWeightNIM );
                my_histos["h_jetpt_6_"+kv.first]->Fill( JetNonIsoMuons_pt_6, totalWeightNIM );
                my_histos["h_jetpt_7_"+kv.first]->Fill( JetNonIsoMuons_pt_7, totalWeightNIM );

            }
        }
    } 
}

void NonIsoMuonBStudy::WriteHistos(TFile* outfile)
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
