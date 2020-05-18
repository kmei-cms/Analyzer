#define ExtraChecks_cxx
#include "Analyzer/Analyzer/include/ExtraChecks.h"
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

ExtraChecks::ExtraChecks()
{
    InitHistos();
}

//Define all your histograms here. 
void ExtraChecks::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    std::vector<std::string> njetTags = {"1", "2", "3", "4", "5", "6", "7" };
    //Define 1D histograms
    
    //Make Histograms for All Input Variables - Regular
    for( std::string njetTag : njetTags ) {
        my_histos.emplace( "h_jet_pt_"+njetTag, std::make_shared<TH1D>( ( "h_jet_pt_"+njetTag ).c_str(), ( "h_jet_pt_"+njetTag ).c_str(), 50, 0, 1000 ) );
        my_histos.emplace( "h_jet_eta_"+njetTag, std::make_shared<TH1D>( ( "h_jet_eta_"+njetTag ).c_str(), ( "h_jet_eta_"+njetTag ).c_str(), 50, -5, 5 ) );
        my_histos.emplace( "h_jet_phi_"+njetTag, std::make_shared<TH1D>( ( "h_jet_phi_"+njetTag ).c_str(), ( "h_jet_phi_"+njetTag ).c_str(), 70, -3.5, 3.5 ) );
        my_histos.emplace( "h_jet_m_"+njetTag, std::make_shared<TH1D>( ( "h_jet_m_"+njetTag ).c_str(), ( "h_jet_m_"+njetTag ).c_str(), 30, 0, 300.0 ) );
    }
    
    my_histos.emplace( "h_lep_pt", std::make_shared<TH1D>( "h_lep_pt", "h_lep_pt", 50, 0, 1000 ) );
    my_histos.emplace( "h_lep_eta", std::make_shared<TH1D>( "h_lep_eta", "h_lep_eta", 50, -5, 5 ) );
    my_histos.emplace( "h_lep_phi", std::make_shared<TH1D>( "h_lep_phi", "h_lep_phi", 70, -3.5, 3.5 ) );
    my_histos.emplace( "h_lep_m", std::make_shared<TH1D>( "h_lep_m", "h_lep_m", 30, 0, 300.0 ) );
        
    my_histos.emplace( "h_fwm2", std::make_shared<TH1D>( "h_fwm2", "h_fwm2", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_fwm3", std::make_shared<TH1D>( "h_fwm3", "h_fwm3", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_fwm4", std::make_shared<TH1D>( "h_fwm4", "h_fwm4", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_fwm5", std::make_shared<TH1D>( "h_fwm5", "h_fwm5", 25, 0.0, 1.0 ) );
    
    my_histos.emplace( "h_jmt0", std::make_shared<TH1D>( "h_jmt0", "h_jmt0", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_jmt1", std::make_shared<TH1D>( "h_jmt1", "h_jmt1", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_jmt2", std::make_shared<TH1D>( "h_jmt2", "h_jmt2", 25, 0.0, 1.0 ) );
    
    for( std::string njetTag : njetTags ) {
        my_histos.emplace( "h_noHT_jet_pt_"+njetTag, std::make_shared<TH1D>( ( "h_noHT_jet_pt_"+njetTag ).c_str(), ( "h_noHT_jet_pt_"+njetTag ).c_str(), 50, 0, 1000 ) );
        my_histos.emplace( "h_noHT_jet_eta_"+njetTag, std::make_shared<TH1D>( ( "h_noHT_jet_eta_"+njetTag ).c_str(), ( "h_noHT_jet_eta_"+njetTag ).c_str(), 50, -5, 5 ) );
        my_histos.emplace( "h_noHT_jet_phi_"+njetTag, std::make_shared<TH1D>( ( "h_noHT_jet_phi_"+njetTag ).c_str(), ( "h_noHT_jet_phi_"+njetTag ).c_str(), 70, -3.5, 3.5 ) );
        my_histos.emplace( "h_noHT_jet_m_"+njetTag, std::make_shared<TH1D>( ( "h_noHT_jet_m_"+njetTag ).c_str(), ( "h_noHT_jet_m_"+njetTag ).c_str(), 30, 0, 300.0 ) );
    }
        
    my_histos.emplace( "h_noHT_lep_pt", std::make_shared<TH1D>( "h_noHT_lep_pt", "h_noHT_lep_pt", 50, 0, 1000 ) );
    my_histos.emplace( "h_noHT_lep_eta", std::make_shared<TH1D>( "h_noHT_lep_eta", "h_noHT_lep_eta", 50, -5, 5 ) );
    my_histos.emplace( "h_noHT_lep_phi", std::make_shared<TH1D>( "h_noHT_lep_phi", "h_noHT_lep_phi", 70, -3.5, 3.5 ) );
    my_histos.emplace( "h_noHT_lep_m", std::make_shared<TH1D>( "h_noHT_lep_m", "h_noHT_lep_m", 30, 0, 300.0 ) );

    my_histos.emplace( "h_noHT_fwm2", std::make_shared<TH1D>( "h_noHT_fwm2", "h_noHT_fwm2", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noHT_fwm3", std::make_shared<TH1D>( "h_noHT_fwm3", "h_noHT_fwm3", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noHT_fwm4", std::make_shared<TH1D>( "h_noHT_fwm4", "h_noHT_fwm4", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noHT_fwm5", std::make_shared<TH1D>( "h_noHT_fwm5", "h_noHT_fwm5", 25, 0.0, 1.0 ) );
    
    my_histos.emplace( "h_noHT_jmt0", std::make_shared<TH1D>( "h_noHT_jmt0", "h_noHT_jmt0", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noHT_jmt1", std::make_shared<TH1D>( "h_noHT_jmt1", "h_noHT_jmt1", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noHT_jmt2", std::make_shared<TH1D>( "h_noHT_jmt2", "h_noHT_jmt2", 25, 0.0, 1.0 ) );
    
    for( std::string njetTag : njetTags ) {
        my_histos.emplace( "h_noSF_jet_pt_"+njetTag, std::make_shared<TH1D>( ( "h_noSF_jet_pt_"+njetTag ).c_str(), ( "h_noSF_jet_pt_"+njetTag ).c_str(), 50, 0, 1000 ) );
        my_histos.emplace( "h_noSF_jet_eta_"+njetTag, std::make_shared<TH1D>( ( "h_noSF_jet_eta_"+njetTag ).c_str(), ( "h_noSF_jet_eta_"+njetTag ).c_str(), 50, -5, 5 ) );
        my_histos.emplace( "h_noSF_jet_phi_"+njetTag, std::make_shared<TH1D>( ( "h_noSF_jet_phi_"+njetTag ).c_str(), ( "h_noSF_jet_phi_"+njetTag ).c_str(), 70, -3.5, 3.5 ) );
        my_histos.emplace( "h_noSF_jet_m_"+njetTag, std::make_shared<TH1D>( ( "h_noSF_jet_m_"+njetTag ).c_str(), ( "h_noSF_jet_m_"+njetTag ).c_str(), 30, 0, 300.0 ) );
    }
        
    my_histos.emplace( "h_noSF_lep_pt", std::make_shared<TH1D>( "h_noSF_lep_pt", "h_noSF_lep_pt", 50, 0, 1000 ) );
    my_histos.emplace( "h_noSF_lep_eta", std::make_shared<TH1D>( "h_noSF_lep_eta", "h_noSF_lep_eta", 50, -5, 5 ) );
    my_histos.emplace( "h_noSF_lep_phi", std::make_shared<TH1D>( "h_noSF_lep_phi", "h_noSF_lep_phi", 70, -3.5, 3.5 ) );
    my_histos.emplace( "h_noSF_lep_m", std::make_shared<TH1D>( "h_noSF_lep_m", "h_noSF_lep_m", 30, 0, 300.0 ) );

    my_histos.emplace( "h_noSF_fwm2", std::make_shared<TH1D>( "h_noSF_fwm2", "h_noSF_fwm2", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noSF_fwm3", std::make_shared<TH1D>( "h_noSF_fwm3", "h_noSF_fwm3", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noSF_fwm4", std::make_shared<TH1D>( "h_noSF_fwm4", "h_noSF_fwm4", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noSF_fwm5", std::make_shared<TH1D>( "h_noSF_fwm5", "h_noSF_fwm5", 25, 0.0, 1.0 ) );
    
    my_histos.emplace( "h_noSF_jmt0", std::make_shared<TH1D>( "h_noSF_jmt0", "h_noSF_jmt0", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noSF_jmt1", std::make_shared<TH1D>( "h_noSF_jmt1", "h_noSF_jmt1", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_noSF_jmt2", std::make_shared<TH1D>( "h_noSF_jmt2", "h_noSF_jmt2", 25, 0.0, 1.0 ) );
    
    //Make Histograms for All Input Variables - Non Iso Muon
    for( std::string njetTag : njetTags ) {
        my_histos.emplace( "h_jet_pt_"+njetTag+"_nim", std::make_shared<TH1D>( ( "h_jet_pt_"+njetTag+"_nim" ).c_str(), ( "h_jet_pt_"+njetTag+"_nim" ).c_str(), 50, 0, 1000 ) );
        my_histos.emplace( "h_jet_eta_"+njetTag+"_nim", std::make_shared<TH1D>( ( "h_jet_eta_"+njetTag+"_nim" ).c_str(), ( "h_jet_eta_"+njetTag+"_nim" ).c_str(), 50, -5, 5 ) );
        my_histos.emplace( "h_jet_phi_"+njetTag+"_nim", std::make_shared<TH1D>( ( "h_jet_phi_"+njetTag+"_nim" ).c_str(), ( "h_jet_phi_"+njetTag+"_nim" ).c_str(), 70, -3.5, 3.5 ) );
        my_histos.emplace( "h_jet_m_"+njetTag+"_nim", std::make_shared<TH1D>( ( "h_jet_m_"+njetTag+"_nim" ).c_str(), ( "h_jet_m_"+njetTag+"_nim" ).c_str(), 30, 0, 300.0 ) );
    }
        
    my_histos.emplace( "h_lep_pt_nim", std::make_shared<TH1D>( "h_lep_pt_nim", "h_lep_pt_nim", 50, 0, 1000 ) );
    my_histos.emplace( "h_lep_eta_nim", std::make_shared<TH1D>( "h_lep_eta_nim", "h_lep_eta_nim", 50, -5, 5 ) );
    my_histos.emplace( "h_lep_phi_nim", std::make_shared<TH1D>( "h_lep_phi_nim", "h_lep_phi_nim", 70, -3.5, 3.5 ) );
    my_histos.emplace( "h_lep_m_nim", std::make_shared<TH1D>( "h_lep_m_nim", "h_lep_m_nim", 30, 0, 300.0 ) );

    my_histos.emplace( "h_fwm2_nim", std::make_shared<TH1D>( "h_fwm2_nim", "h_fwm2_nim", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_fwm3_nim", std::make_shared<TH1D>( "h_fwm3_nim", "h_fwm3_nim", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_fwm4_nim", std::make_shared<TH1D>( "h_fwm4_nim", "h_fwm4_nim", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_fwm5_nim", std::make_shared<TH1D>( "h_fwm5_nim", "h_fwm5_nim", 25, 0.0, 1.0 ) );
    
    my_histos.emplace( "h_jmt0_nim", std::make_shared<TH1D>( "h_jmt0_nim", "h_jmt0_nim", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_jmt1_nim", std::make_shared<TH1D>( "h_jmt1_nim", "h_jmt1_nim", 25, 0.0, 1.0 ) );
    my_histos.emplace( "h_jmt2_nim", std::make_shared<TH1D>( "h_jmt2_nim", "h_jmt2_nim", 25, 0.0, 1.0 ) );

    //Make Historgrams for New TTX Baseline
    my_histos.emplace( "h_njets_ttx", std::make_shared<TH1D>( "h_njets_ttx", "h_njets_ttx", 20, 0, 20 ) ); 
    my_histos.emplace( "h_ht_ttx", std::make_shared<TH1D>( "h_ht_ttx", "h_ht_ttx", 40, 0, 4000 ) ); 
    my_histos.emplace( "h_deepESM_ttx", std::make_shared<TH1D>( "h_deepESM_ttx", "h_deepESM_ttx", 15, 0, 1.0 ) ); 

}

//Put everything you want to do per event here.
void ExtraChecks::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
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
        if( tr.getEvtNum() && 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );
        
        const auto& lumi             = tr.getVar<double>("Lumi");
        const auto& runtype          = tr.getVar<std::string>("runtype");     
        const auto& filetag          = tr.getVar<std::string>("filetag");     
        
        const auto& passBaseline1l   = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline1l_NIM     = tr.getVar<bool>("passBaseline1l_NonIsoMuon");

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;
        if( !isQuiet ) { std::cout<<weight<<std::endl; }

        // ------------------------
        // -- Define weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
        // ------------------------
        double totalWeight          = 1.0;
        double nohtweight           = 1.0;
        double totalWeightNIM       = 1.0;
        double eventweight          = 1.0;

        if(runtype == "MC")
        {
            const auto& weight = tr.getVar<double>("Weight");
            const auto& totalEventWeight = tr.getVar<double>("totalEventWeight");
            const auto& totalEventWeightNIM = tr.getVar<double>("totalEventWeightNIM");
            const auto& htDerivedweight = tr.getVar<double>("htDerivedweight");
         
            eventweight = weight*lumi;
            totalWeight = totalEventWeight*lumi;
            nohtweight = totalEventWeight*lumi/htDerivedweight;
            totalWeightNIM = totalEventWeightNIM*lumi;
        }

        //Declare branches for the variables in MakeMVAVariables.h
        const auto& Jet_pt_1    = tr.getVar<double>("Jet_pt_1");
        const auto& Jet_eta_1    = tr.getVar<double>("Jet_eta_1");
        const auto& Jet_phi_1    = tr.getVar<double>("Jet_phi_1");
        const auto& Jet_m_1    = tr.getVar<double>("Jet_m_1");
        const auto& Jet_pt_2    = tr.getVar<double>("Jet_pt_2");
        const auto& Jet_eta_2    = tr.getVar<double>("Jet_eta_2");
        const auto& Jet_phi_2    = tr.getVar<double>("Jet_phi_2");
        const auto& Jet_m_2    = tr.getVar<double>("Jet_m_2");
        const auto& Jet_pt_3    = tr.getVar<double>("Jet_pt_3");
        const auto& Jet_eta_3    = tr.getVar<double>("Jet_eta_3");
        const auto& Jet_phi_3    = tr.getVar<double>("Jet_phi_3");
        const auto& Jet_m_3    = tr.getVar<double>("Jet_m_3");
        const auto& Jet_pt_4    = tr.getVar<double>("Jet_pt_4");
        const auto& Jet_eta_4    = tr.getVar<double>("Jet_eta_4");
        const auto& Jet_phi_4    = tr.getVar<double>("Jet_phi_4");
        const auto& Jet_m_4    = tr.getVar<double>("Jet_m_4");
        const auto& Jet_pt_5    = tr.getVar<double>("Jet_pt_5");
        const auto& Jet_eta_5    = tr.getVar<double>("Jet_eta_5");
        const auto& Jet_phi_5    = tr.getVar<double>("Jet_phi_5");
        const auto& Jet_m_5    = tr.getVar<double>("Jet_m_5");
        const auto& Jet_pt_6    = tr.getVar<double>("Jet_pt_6");
        const auto& Jet_eta_6    = tr.getVar<double>("Jet_eta_6");
        const auto& Jet_phi_6    = tr.getVar<double>("Jet_phi_6");
        const auto& Jet_m_6    = tr.getVar<double>("Jet_m_6");
        const auto& Jet_pt_7    = tr.getVar<double>("Jet_pt_7");
        const auto& Jet_eta_7    = tr.getVar<double>("Jet_eta_7");
        const auto& Jet_phi_7    = tr.getVar<double>("Jet_phi_7");
        const auto& Jet_m_7    = tr.getVar<double>("Jet_m_7");
        const auto& GoodLeptons_pt_1    = tr.getVar<double>("GoodLeptons_pt_1");
        const auto& GoodLeptons_eta_1    = tr.getVar<double>("GoodLeptons_eta_1");
        const auto& GoodLeptons_phi_1    = tr.getVar<double>("GoodLeptons_phi_1");
        const auto& GoodLeptons_m_1    = tr.getVar<double>("GoodLeptons_m_1");
        const auto& fwm2_top6    = tr.getVar<double>("fwm2_top6");
        const auto& fwm3_top6    = tr.getVar<double>("fwm3_top6");
        const auto& fwm4_top6    = tr.getVar<double>("fwm4_top6");
        const auto& fwm5_top6    = tr.getVar<double>("fwm5_top6");
        const auto& jmt_ev0_top6    = tr.getVar<double>("jmt_ev0_top6");
        const auto& jmt_ev1_top6    = tr.getVar<double>("jmt_ev1_top6");
        const auto& jmt_ev2_top6    = tr.getVar<double>("jmt_ev2_top6");


        if( passBaseline1l ) {
            my_histos["h_jet_pt_1"]->Fill( Jet_pt_1, totalWeight );
            my_histos["h_jet_eta_1"]->Fill( Jet_eta_1, totalWeight );
            my_histos["h_jet_phi_1"]->Fill( Jet_phi_1, totalWeight );
            my_histos["h_jet_m_1"]->Fill( Jet_m_1, totalWeight );
            my_histos["h_jet_pt_2"]->Fill( Jet_pt_2, totalWeight );
            my_histos["h_jet_eta_2"]->Fill( Jet_eta_2, totalWeight );
            my_histos["h_jet_phi_2"]->Fill( Jet_phi_2, totalWeight );
            my_histos["h_jet_m_2"]->Fill( Jet_m_2, totalWeight );
            my_histos["h_jet_pt_3"]->Fill( Jet_pt_3, totalWeight );
            my_histos["h_jet_eta_3"]->Fill( Jet_eta_3, totalWeight );
            my_histos["h_jet_phi_3"]->Fill( Jet_phi_3, totalWeight );
            my_histos["h_jet_m_3"]->Fill( Jet_m_3, totalWeight );
            my_histos["h_jet_pt_4"]->Fill( Jet_pt_4, totalWeight );
            my_histos["h_jet_eta_4"]->Fill( Jet_eta_4, totalWeight );
            my_histos["h_jet_phi_4"]->Fill( Jet_phi_4, totalWeight );
            my_histos["h_jet_m_4"]->Fill( Jet_m_4, totalWeight );
            my_histos["h_jet_pt_5"]->Fill( Jet_pt_5, totalWeight );
            my_histos["h_jet_eta_5"]->Fill( Jet_eta_5, totalWeight );
            my_histos["h_jet_phi_5"]->Fill( Jet_phi_5, totalWeight );
            my_histos["h_jet_m_5"]->Fill( Jet_m_5, totalWeight );
            my_histos["h_jet_pt_6"]->Fill( Jet_pt_6, totalWeight );
            my_histos["h_jet_eta_6"]->Fill( Jet_eta_6, totalWeight );
            my_histos["h_jet_phi_6"]->Fill( Jet_phi_6, totalWeight );
            my_histos["h_jet_m_6"]->Fill( Jet_m_6, totalWeight );
            my_histos["h_jet_pt_7"]->Fill( Jet_pt_7, totalWeight );
            my_histos["h_jet_eta_7"]->Fill( Jet_eta_7, totalWeight );
            my_histos["h_jet_phi_7"]->Fill( Jet_phi_7, totalWeight );
            my_histos["h_jet_m_7"]->Fill( Jet_m_7, totalWeight );
            my_histos["h_lep_pt"]->Fill( GoodLeptons_pt_1, totalWeight );
            my_histos["h_lep_eta"]->Fill( GoodLeptons_eta_1, totalWeight );
            my_histos["h_lep_phi"]->Fill( GoodLeptons_phi_1, totalWeight );
            my_histos["h_lep_m"]->Fill( GoodLeptons_m_1, totalWeight );
            my_histos["h_fwm2"]->Fill( fwm2_top6, totalWeight );
            my_histos["h_fwm3"]->Fill( fwm3_top6, totalWeight );
            my_histos["h_fwm4"]->Fill( fwm4_top6, totalWeight );
            my_histos["h_fwm5"]->Fill( fwm5_top6, totalWeight );
            my_histos["h_jmt0"]->Fill( jmt_ev0_top6, totalWeight );
            my_histos["h_jmt1"]->Fill( jmt_ev1_top6, totalWeight );
            my_histos["h_jmt2"]->Fill( jmt_ev2_top6, totalWeight );
            
            my_histos["h_noHT_jet_pt_1"]->Fill( Jet_pt_1, nohtweight );
            my_histos["h_noHT_jet_eta_1"]->Fill( Jet_eta_1, nohtweight );
            my_histos["h_noHT_jet_phi_1"]->Fill( Jet_phi_1, nohtweight );
            my_histos["h_noHT_jet_m_1"]->Fill( Jet_m_1, nohtweight );
            my_histos["h_noHT_jet_pt_2"]->Fill( Jet_pt_2, nohtweight );
            my_histos["h_noHT_jet_eta_2"]->Fill( Jet_eta_2, nohtweight );
            my_histos["h_noHT_jet_phi_2"]->Fill( Jet_phi_2, nohtweight );
            my_histos["h_noHT_jet_m_2"]->Fill( Jet_m_2, nohtweight );
            my_histos["h_noHT_jet_pt_3"]->Fill( Jet_pt_3, nohtweight );
            my_histos["h_noHT_jet_eta_3"]->Fill( Jet_eta_3, nohtweight );
            my_histos["h_noHT_jet_phi_3"]->Fill( Jet_phi_3, nohtweight );
            my_histos["h_noHT_jet_m_3"]->Fill( Jet_m_3, nohtweight );
            my_histos["h_noHT_jet_pt_4"]->Fill( Jet_pt_4, nohtweight );
            my_histos["h_noHT_jet_eta_4"]->Fill( Jet_eta_4, nohtweight );
            my_histos["h_noHT_jet_phi_4"]->Fill( Jet_phi_4, nohtweight );
            my_histos["h_noHT_jet_m_4"]->Fill( Jet_m_4, nohtweight );
            my_histos["h_noHT_jet_pt_5"]->Fill( Jet_pt_5, nohtweight );
            my_histos["h_noHT_jet_eta_5"]->Fill( Jet_eta_5, nohtweight );
            my_histos["h_noHT_jet_phi_5"]->Fill( Jet_phi_5, nohtweight );
            my_histos["h_noHT_jet_m_5"]->Fill( Jet_m_5, nohtweight );
            my_histos["h_noHT_jet_pt_6"]->Fill( Jet_pt_6, nohtweight );
            my_histos["h_noHT_jet_eta_6"]->Fill( Jet_eta_6, nohtweight );
            my_histos["h_noHT_jet_phi_6"]->Fill( Jet_phi_6, nohtweight );
            my_histos["h_noHT_jet_m_6"]->Fill( Jet_m_6, nohtweight );
            my_histos["h_noHT_jet_pt_7"]->Fill( Jet_pt_7, nohtweight );
            my_histos["h_noHT_jet_eta_7"]->Fill( Jet_eta_7, nohtweight );
            my_histos["h_noHT_jet_phi_7"]->Fill( Jet_phi_7, nohtweight );
            my_histos["h_noHT_jet_m_7"]->Fill( Jet_m_7, nohtweight );
            my_histos["h_noHT_lep_pt"]->Fill( GoodLeptons_pt_1, nohtweight );
            my_histos["h_noHT_lep_eta"]->Fill( GoodLeptons_eta_1, nohtweight );
            my_histos["h_noHT_lep_phi"]->Fill( GoodLeptons_phi_1, nohtweight );
            my_histos["h_noHT_lep_m"]->Fill( GoodLeptons_m_1, nohtweight );
            my_histos["h_noHT_fwm2"]->Fill( fwm2_top6, nohtweight );
            my_histos["h_noHT_fwm3"]->Fill( fwm3_top6, nohtweight );
            my_histos["h_noHT_fwm4"]->Fill( fwm4_top6, nohtweight );
            my_histos["h_noHT_fwm5"]->Fill( fwm5_top6, nohtweight );
            my_histos["h_noHT_jmt0"]->Fill( jmt_ev0_top6, nohtweight );
            my_histos["h_noHT_jmt1"]->Fill( jmt_ev1_top6, nohtweight );
            my_histos["h_noHT_jmt2"]->Fill( jmt_ev2_top6, nohtweight );
            
            my_histos["h_noSF_jet_pt_1"]->Fill( Jet_pt_1, eventweight );
            my_histos["h_noSF_jet_eta_1"]->Fill( Jet_eta_1, eventweight );
            my_histos["h_noSF_jet_phi_1"]->Fill( Jet_phi_1, eventweight );
            my_histos["h_noSF_jet_m_1"]->Fill( Jet_m_1, eventweight );
            my_histos["h_noSF_jet_pt_2"]->Fill( Jet_pt_2, eventweight );
            my_histos["h_noSF_jet_eta_2"]->Fill( Jet_eta_2, eventweight );
            my_histos["h_noSF_jet_phi_2"]->Fill( Jet_phi_2, eventweight );
            my_histos["h_noSF_jet_m_2"]->Fill( Jet_m_2, eventweight );
            my_histos["h_noSF_jet_pt_3"]->Fill( Jet_pt_3, eventweight );
            my_histos["h_noSF_jet_eta_3"]->Fill( Jet_eta_3, eventweight );
            my_histos["h_noSF_jet_phi_3"]->Fill( Jet_phi_3, eventweight );
            my_histos["h_noSF_jet_m_3"]->Fill( Jet_m_3, eventweight );
            my_histos["h_noSF_jet_pt_4"]->Fill( Jet_pt_4, eventweight );
            my_histos["h_noSF_jet_eta_4"]->Fill( Jet_eta_4, eventweight );
            my_histos["h_noSF_jet_phi_4"]->Fill( Jet_phi_4, eventweight );
            my_histos["h_noSF_jet_m_4"]->Fill( Jet_m_4, eventweight );
            my_histos["h_noSF_jet_pt_5"]->Fill( Jet_pt_5, eventweight );
            my_histos["h_noSF_jet_eta_5"]->Fill( Jet_eta_5, eventweight );
            my_histos["h_noSF_jet_phi_5"]->Fill( Jet_phi_5, eventweight );
            my_histos["h_noSF_jet_m_5"]->Fill( Jet_m_5, eventweight );
            my_histos["h_noSF_jet_pt_6"]->Fill( Jet_pt_6, eventweight );
            my_histos["h_noSF_jet_eta_6"]->Fill( Jet_eta_6, eventweight );
            my_histos["h_noSF_jet_phi_6"]->Fill( Jet_phi_6, eventweight );
            my_histos["h_noSF_jet_m_6"]->Fill( Jet_m_6, eventweight );
            my_histos["h_noSF_jet_pt_7"]->Fill( Jet_pt_7, eventweight );
            my_histos["h_noSF_jet_eta_7"]->Fill( Jet_eta_7, eventweight );
            my_histos["h_noSF_jet_phi_7"]->Fill( Jet_phi_7, eventweight );
            my_histos["h_noSF_jet_m_7"]->Fill( Jet_m_7, eventweight );
            my_histos["h_noSF_lep_pt"]->Fill( GoodLeptons_pt_1, eventweight );
            my_histos["h_noSF_lep_eta"]->Fill( GoodLeptons_eta_1, eventweight );
            my_histos["h_noSF_lep_phi"]->Fill( GoodLeptons_phi_1, eventweight );
            my_histos["h_noSF_lep_m"]->Fill( GoodLeptons_m_1, eventweight );
            my_histos["h_noSF_fwm2"]->Fill( fwm2_top6, eventweight );
            my_histos["h_noSF_fwm3"]->Fill( fwm3_top6, eventweight );
            my_histos["h_noSF_fwm4"]->Fill( fwm4_top6, eventweight );
            my_histos["h_noSF_fwm5"]->Fill( fwm5_top6, eventweight );
            my_histos["h_noSF_jmt0"]->Fill( jmt_ev0_top6, eventweight );
            my_histos["h_noSF_jmt1"]->Fill( jmt_ev1_top6, eventweight );
            my_histos["h_noSF_jmt2"]->Fill( jmt_ev2_top6, eventweight );
        }
        
        const auto& JetNIM_pt_1    = tr.getVar<double>("JetNonIsoMuons_pt_1");
        const auto& JetNIM_eta_1    = tr.getVar<double>("JetNonIsoMuons_eta_1");
        const auto& JetNIM_phi_1    = tr.getVar<double>("JetNonIsoMuons_phi_1");
        const auto& JetNIM_m_1    = tr.getVar<double>("JetNonIsoMuons_m_1");
        const auto& JetNIM_pt_2    = tr.getVar<double>("JetNonIsoMuons_pt_2");
        const auto& JetNIM_eta_2    = tr.getVar<double>("JetNonIsoMuons_eta_2");
        const auto& JetNIM_phi_2    = tr.getVar<double>("JetNonIsoMuons_phi_2");
        const auto& JetNIM_m_2    = tr.getVar<double>("JetNonIsoMuons_m_2");
        const auto& JetNIM_pt_3    = tr.getVar<double>("JetNonIsoMuons_pt_3");
        const auto& JetNIM_eta_3    = tr.getVar<double>("JetNonIsoMuons_eta_3");
        const auto& JetNIM_phi_3    = tr.getVar<double>("JetNonIsoMuons_phi_3");
        const auto& JetNIM_m_3    = tr.getVar<double>("JetNonIsoMuons_m_3");
        const auto& JetNIM_pt_4    = tr.getVar<double>("JetNonIsoMuons_pt_4");
        const auto& JetNIM_eta_4    = tr.getVar<double>("JetNonIsoMuons_eta_4");
        const auto& JetNIM_phi_4    = tr.getVar<double>("JetNonIsoMuons_phi_4");
        const auto& JetNIM_m_4    = tr.getVar<double>("JetNonIsoMuons_m_4");
        const auto& JetNIM_pt_5    = tr.getVar<double>("JetNonIsoMuons_pt_5");
        const auto& JetNIM_eta_5    = tr.getVar<double>("JetNonIsoMuons_eta_5");
        const auto& JetNIM_phi_5    = tr.getVar<double>("JetNonIsoMuons_phi_5");
        const auto& JetNIM_m_5    = tr.getVar<double>("JetNonIsoMuons_m_5");
        const auto& JetNIM_pt_6    = tr.getVar<double>("JetNonIsoMuons_pt_6");
        const auto& JetNIM_eta_6    = tr.getVar<double>("JetNonIsoMuons_eta_6");
        const auto& JetNIM_phi_6    = tr.getVar<double>("JetNonIsoMuons_phi_6");
        const auto& JetNIM_m_6    = tr.getVar<double>("JetNonIsoMuons_m_6");
        const auto& JetNIM_pt_7    = tr.getVar<double>("JetNonIsoMuons_pt_7");
        const auto& JetNIM_eta_7    = tr.getVar<double>("JetNonIsoMuons_eta_7");
        const auto& JetNIM_phi_7    = tr.getVar<double>("JetNonIsoMuons_phi_7");
        const auto& JetNIM_m_7    = tr.getVar<double>("JetNonIsoMuons_m_7");
        const auto& GoodNonIsoMuons_pt_1    = tr.getVar<double>("GoodNonIsoMuons_pt_1");
        const auto& GoodNonIsoMuons_eta_1    = tr.getVar<double>("GoodNonIsoMuons_eta_1");
        const auto& GoodNonIsoMuons_phi_1    = tr.getVar<double>("GoodNonIsoMuons_phi_1");
        const auto& GoodNonIsoMuons_m_1    = tr.getVar<double>("GoodNonIsoMuons_m_1");
        const auto& NonIsoMuons_fwm2_top6    = tr.getVar<double>("NonIsoMuons_fwm2_top6");
        const auto& NonIsoMuons_fwm3_top6    = tr.getVar<double>("NonIsoMuons_fwm3_top6");
        const auto& NonIsoMuons_fwm4_top6    = tr.getVar<double>("NonIsoMuons_fwm4_top6");
        const auto& NonIsoMuons_fwm5_top6    = tr.getVar<double>("NonIsoMuons_fwm5_top6");
        const auto& NonIsoMuons_jmt_ev0_top6    = tr.getVar<double>("NonIsoMuons_jmt_ev0_top6");
        const auto& NonIsoMuons_jmt_ev1_top6    = tr.getVar<double>("NonIsoMuons_jmt_ev1_top6");
        const auto& NonIsoMuons_jmt_ev2_top6    = tr.getVar<double>("NonIsoMuons_jmt_ev2_top6");
        
        if( passBaseline1l_NIM ) {
            my_histos["h_jet_pt_1_nim"]->Fill( JetNIM_pt_1, totalWeightNIM );
            my_histos["h_jet_eta_1_nim"]->Fill( JetNIM_eta_1, totalWeightNIM );
            my_histos["h_jet_phi_1_nim"]->Fill( JetNIM_phi_1, totalWeightNIM );
            my_histos["h_jet_m_1_nim"]->Fill( JetNIM_m_1, totalWeightNIM );
            my_histos["h_jet_pt_2_nim"]->Fill( JetNIM_pt_2, totalWeightNIM );
            my_histos["h_jet_eta_2_nim"]->Fill( JetNIM_eta_2, totalWeightNIM );
            my_histos["h_jet_phi_2_nim"]->Fill( JetNIM_phi_2, totalWeightNIM );
            my_histos["h_jet_m_2_nim"]->Fill( JetNIM_m_2, totalWeightNIM );
            my_histos["h_jet_pt_3_nim"]->Fill( JetNIM_pt_3, totalWeightNIM );
            my_histos["h_jet_eta_3_nim"]->Fill( JetNIM_eta_3, totalWeightNIM );
            my_histos["h_jet_phi_3_nim"]->Fill( JetNIM_phi_3, totalWeightNIM );
            my_histos["h_jet_m_3_nim"]->Fill( JetNIM_m_3, totalWeightNIM );
            my_histos["h_jet_pt_4_nim"]->Fill( JetNIM_pt_4, totalWeightNIM );
            my_histos["h_jet_eta_4_nim"]->Fill( JetNIM_eta_4, totalWeightNIM );
            my_histos["h_jet_phi_4_nim"]->Fill( JetNIM_phi_4, totalWeightNIM );
            my_histos["h_jet_m_4_nim"]->Fill( JetNIM_m_4, totalWeightNIM );
            my_histos["h_jet_pt_5_nim"]->Fill( JetNIM_pt_5, totalWeightNIM );
            my_histos["h_jet_eta_5_nim"]->Fill( JetNIM_eta_5, totalWeightNIM );
            my_histos["h_jet_phi_5_nim"]->Fill( JetNIM_phi_5, totalWeightNIM );
            my_histos["h_jet_m_5_nim"]->Fill( JetNIM_m_5, totalWeightNIM );
            my_histos["h_jet_pt_6_nim"]->Fill( JetNIM_pt_6, totalWeightNIM );
            my_histos["h_jet_eta_6_nim"]->Fill( JetNIM_eta_6, totalWeightNIM );
            my_histos["h_jet_phi_6_nim"]->Fill( JetNIM_phi_6, totalWeightNIM );
            my_histos["h_jet_m_6_nim"]->Fill( JetNIM_m_6, totalWeightNIM );
            my_histos["h_jet_pt_7_nim"]->Fill( JetNIM_pt_7, totalWeightNIM );
            my_histos["h_jet_eta_7_nim"]->Fill( JetNIM_eta_7, totalWeightNIM );
            my_histos["h_jet_phi_7_nim"]->Fill( JetNIM_phi_7, totalWeightNIM );
            my_histos["h_jet_m_7_nim"]->Fill( JetNIM_m_7, totalWeightNIM );
            my_histos["h_lep_pt_nim"]->Fill( GoodNonIsoMuons_pt_1, totalWeightNIM );
            my_histos["h_lep_eta_nim"]->Fill( GoodNonIsoMuons_eta_1, totalWeightNIM );
            my_histos["h_lep_phi_nim"]->Fill( GoodNonIsoMuons_phi_1, totalWeightNIM );
            my_histos["h_lep_m_nim"]->Fill( GoodNonIsoMuons_m_1, totalWeightNIM );
            my_histos["h_fwm2_nim"]->Fill( NonIsoMuons_fwm2_top6, totalWeightNIM );
            my_histos["h_fwm3_nim"]->Fill( NonIsoMuons_fwm3_top6, totalWeightNIM );
            my_histos["h_fwm4_nim"]->Fill( NonIsoMuons_fwm4_top6, totalWeightNIM );
            my_histos["h_fwm5_nim"]->Fill( NonIsoMuons_fwm5_top6, totalWeightNIM );
            my_histos["h_jmt0_nim"]->Fill( NonIsoMuons_jmt_ev0_top6, totalWeightNIM );
            my_histos["h_jmt1_nim"]->Fill( NonIsoMuons_jmt_ev1_top6, totalWeightNIM );
            my_histos["h_jmt2_nim"]->Fill( NonIsoMuons_jmt_ev2_top6, totalWeightNIM );
        }

        //Making Hannjeorg region
        const auto& passBaselineGoodOffline1l       = tr.getVar<bool>("passBaselineGoodOffline1l");
        const auto& NGoodLeptons                    = tr.getVar<int>("NGoodLeptons");
        const auto& allElectrons                    = tr.getVec<TLorentzVector>("Electrons");
        const auto& allElectrons_passIso            = tr.getVec<bool>("Electrons_passIso");
        const auto& allElectrons_charge             = tr.getVec<int>("Electrons_charge");
        const auto& allElectrons_tightID            = tr.getVec<bool>("Electrons_tightID");
        
        const auto& allMuons                        = tr.getVec<TLorentzVector>("Muons");
        const auto& allMuons_passIso                = tr.getVec<bool>("Muons_passIso");
        const auto& allMuons_charge                 = tr.getVec<int>("Muons_charge");
        const auto& allMuons_mediumID               = tr.getVec<bool>("Muons_mediumID");
        
        double  nPosEl_pt15                         = 0.0;
        double  nNegEl_pt15                         = 0.0;
        double  nPosMu_pt15                         = 0.0;
        double  nNegMu_pt15                         = 0.0;
        
        for( unsigned int iel = 0; iel < allElectrons.size(); ++iel ) {
            TLorentzVector lvel = allElectrons.at(iel);
            if( abs( lvel.Eta() ) < 2.4 && lvel.Pt() > 15.0 && allElectrons_passIso.at(iel) && allElectrons_tightID.at(iel) ) {
                if( allElectrons_charge.at(iel) > 0 ) nPosEl_pt15++;
                else if( allElectrons_charge.at(iel) < 0 ) nNegEl_pt15++;
            }
        }
        for( unsigned int imu = 0; imu < allMuons.size(); ++imu ) {
            TLorentzVector lvel = allMuons.at(imu);
            if( abs( lvel.Eta() ) < 2.4 && lvel.Pt() > 15.0 && allMuons_passIso.at(imu) && allMuons_mediumID.at(imu) ) {
                if( allMuons_charge.at(imu) > 0 ) nPosMu_pt15++;
                else if( allMuons_charge.at(imu) < 0 ) nNegMu_pt15++;
            }
        }

        const auto& passTrigger                     = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC                   = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodJets_pt30                  = tr.getVar<int>("NGoodJets_pt30");
        const auto& HT_trigger_pt30                 = tr.getVar<double>("HT_trigger_pt30");
        const auto& deepESM_val                     = tr.getVar<double>("deepESM_val");

        bool goodMuonEvent                          = (runtype != "Data" || filetag.find("Data_SingleElectron") !=     std::string::npos);
        bool goodElectronEvent                      = (runtype != "Data" || filetag.find("Data_SingleElectron") !=     std::string::npos);

        double totalLeptons                         = nPosEl_pt15 + nNegEl_pt15 + nPosMu_pt15 + nNegMu_pt15;
        if( passBaselineGoodOffline1l && NGoodLeptons >= 1 && totalLeptons >= 2 && goodMuonEvent && goodElectronEvent && NGoodJets_pt30 >= 7 && passTrigger && passTriggerMC ) {

            if( passBaseline1l ) continue;
            if( totalLeptons == 2 ) {
                if( ( ( nPosEl_pt15 + nPosMu_pt15 ) == 2 ) || ( nNegEl_pt15 + nNegMu_pt15 ) == 2 )
                {
                    my_histos["h_njets_ttx"]->Fill( NGoodJets_pt30, eventweight);
                    my_histos["h_ht_ttx"]->Fill( HT_trigger_pt30, eventweight);
                    my_histos["h_deepESM_ttx"]->Fill( deepESM_val, eventweight);
                }
            }
            else {
                my_histos["h_njets_ttx"]->Fill( NGoodJets_pt30, eventweight);
                my_histos["h_ht_ttx"]->Fill( HT_trigger_pt30, eventweight);
                my_histos["h_deepESM_ttx"]->Fill( deepESM_val, eventweight);
            }
        }
    } 
}

void ExtraChecks::WriteHistos(TFile* outfile)
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
