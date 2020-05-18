#define MakeQCDChecks_cxx
#include "Analyzer/Analyzer/include/MakeQCDChecks.h"
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

MakeQCDChecks::MakeQCDChecks()
{
    InitHistos();
}

//Define all your histograms here. 
void MakeQCDChecks::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    std::vector<std::string> qcdSamples  =  { "200to300", "300to500", "500to700", "700to1000", "1000to1500", "1500to2000", "2000toInf", "All" };

    for( std::string qcdSampleString : qcdSamples ) {
        my_histos.emplace( "h_cutflow_"+qcdSampleString, std::make_shared<TH1D>( ( "h_cutflow_"+qcdSampleString ).c_str(), ( "h_cutflow_"+qcdSampleString ).c_str(), 11, 0, 11 ) );
        
        my_histos.emplace( "h_npv_noTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_npv_noTrig_"+qcdSampleString ).c_str(), ( "h_npv_noTrig_"+qcdSampleString ).c_str(), 80, 0, 80 ) );
        my_histos.emplace( "h_npv_passTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_npv_passTrig_"+qcdSampleString ).c_str(), ( "h_npv_passTrig_"+qcdSampleString ).c_str(), 80, 0, 80 ) );
        my_histos.emplace( "h_npv_noTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_npv_noTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_npv_noTrig_passBaseline_"+qcdSampleString ).c_str(), 80, 0, 80 ) );
        my_histos.emplace( "h_npv_passTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_npv_passTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_npv_passTrig_passBaseline_"+qcdSampleString ).c_str(), 80, 0, 80 ) );

    }
    for( std::string qcdSampleString : qcdSamples ) {
        my_histos.emplace( "h_madHT_"+qcdSampleString, std::make_shared<TH1D>( ( "h_madHT_"+qcdSampleString ).c_str(), ( "h_madHT_"+qcdSampleString ).c_str(), 100, 0, 5000 ) );
        my_histos.emplace( "h_njets_noTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_njets_noTrig_"+qcdSampleString ).c_str(), ( "h_njets_noTrig_"+qcdSampleString ).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_njets_passTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_njets_passTrig_"+qcdSampleString ).c_str(), ( "h_njets_passTrig_"+qcdSampleString ).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_ht_noTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_ht_noTrig_"+qcdSampleString ).c_str(), ( "h_ht_noTrig_"+qcdSampleString ).c_str(), 50, 0, 5000 ) );
        my_histos.emplace( "h_ht_passTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_ht_passTrig_"+qcdSampleString ).c_str(), ( "h_ht_passTrig_"+qcdSampleString ).c_str(), 50, 0, 5000 ) );
        
        my_histos.emplace( "h_tw_madHT_"+qcdSampleString, std::make_shared<TH1D>( ( "h_tw_madHT_"+qcdSampleString ).c_str(), ( "h_tw_madHT_"+qcdSampleString ).c_str(), 100, 0, 5000 ) );
        my_histos.emplace( "h_tw_njets_noTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_tw_njets_noTrig_"+qcdSampleString ).c_str(), ( "h_tw_njets_noTrig_"+qcdSampleString ).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_tw_njets_passTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_tw_njets_passTrig_"+qcdSampleString ).c_str(), ( "h_tw_njets_passTrig_"+qcdSampleString ).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_tw_ht_noTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_tw_ht_noTrig_"+qcdSampleString ).c_str(), ( "h_tw_ht_noTrig_"+qcdSampleString ).c_str(), 50, 0, 5000 ) );
        my_histos.emplace( "h_tw_ht_passTrig_"+qcdSampleString, std::make_shared<TH1D>( ( "h_tw_ht_passTrig_"+qcdSampleString ).c_str(), ( "h_tw_ht_passTrig_"+qcdSampleString ).c_str(), 50, 0, 5000 ) );
        
        my_histos.emplace( "h_njets_noTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_njets_noTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_njets_noTrig_passBaseline_"+qcdSampleString ).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_njets_passTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_njets_passTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_njets_passTrig_passBaseline_"+qcdSampleString ).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_ht_noTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_ht_noTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_ht_noTrig_passBaseline_"+qcdSampleString ).c_str(), 50, 0, 5000 ) );
        my_histos.emplace( "h_ht_passTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_ht_passTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_ht_passTrig_passBaseline_"+qcdSampleString ).c_str(), 50, 0, 5000 ) );
        
        my_histos.emplace( "h_tw_njets_noTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_tw_njets_noTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_tw_njets_noTrig_passBaseline_"+qcdSampleString ).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_tw_njets_passTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_tw_njets_passTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_tw_njets_passTrig_passBaseline_"+qcdSampleString ).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_tw_ht_noTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_tw_ht_noTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_tw_ht_noTrig_passBaseline_"+qcdSampleString ).c_str(), 50, 0, 5000 ) );
        my_histos.emplace( "h_tw_ht_passTrig_passBaseline_"+qcdSampleString, std::make_shared<TH1D>( ( "h_tw_ht_passTrig_passBaseline_"+qcdSampleString ).c_str(), ( "h_tw_ht_passTrig_passBaseline_"+qcdSampleString ).c_str(), 50, 0, 5000 ) );
    }

}

//Put everything you want to do per event here.
void MakeQCDChecks::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
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
        
        const auto& lumi                    = tr.getVar<double>("Lumi");
        const auto& runtype                 = tr.getVar<std::string>("runtype");     
        const auto& filetag                 = tr.getVar<std::string>("filetag");
        
        const auto& passMadHT               = tr.getVar<bool>("passMadHT");
        const auto& NNonIsoMuonJets_pt30    = tr.getVar<int>("NNonIsoMuonJets_pt30");
        const auto& HT_NonIsoMuon_pt30      = tr.getVar<double>("HT_NonIsoMuon_pt30");
        const auto& passBaseline            = tr.getVar<bool>("passBaseline1l_NonIsoMuon_noTrigger");

        const auto& TriggerNames            = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass             = tr.getVec<int>("TriggerPass");
       
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------
        // -- Define weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
        // ------------------------
        double eventweight          = 1.0;
        double totalWeightNIM       = 1.0;
        double madHT                = 1.0;

        if( isQuiet ) {
            if( weight < -9999 ) continue;
        }

        if(runtype == "MC")
        {
            const auto& madht                   = tr.getVar<double>("madHT");
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& totalEventWeightNIM = tr.getVar<double>("totalEventWeightNIM");
            
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
            eventweight = lumi*Weight;
            totalWeightNIM = totalEventWeightNIM*lumi;
            madHT = madht;
        }

        std::vector<std::string> myNonIsoMuonTriggers { "HLT_Mu50_v", "HLT_TkMu50_v" };
        bool passNonIsoMuonTriggers = passTriggerGeneral( myNonIsoMuonTriggers, TriggerNames, TriggerPass );
       
        const std::map<std::string, bool> cut_map {
            { "200to300"                    , filetag.find("200to300") != std::string::npos },
            { "300to500"                    , filetag.find("300to500") != std::string::npos },
            { "500to700"                    , filetag.find("500to700") != std::string::npos },
            { "700to1000"                   , filetag.find("700to1000") != std::string::npos },
            { "1000to1500"                  , filetag.find("1000to1500") != std::string::npos },
            { "1500to2000"                  , filetag.find("1500to2000") != std::string::npos },
            { "2000toInf"                   , filetag.find("2000toInf") != std::string::npos },
            { "All"                         , true }
        };

        const auto& correct2018Split        = tr.getVar<bool>( "correct2018Split" );
        const auto& passHEMVeto             = tr.getVar<bool>( "passHEMVeto" );
        const auto& passMETFilters          = tr.getVar<bool>( "passMETFilters" );
        const auto& NNonIsoMuons            = tr.getVar<int>( "NNonIsoMuons" );
        const auto& NGoodElectrons          = tr.getVar<int>( "NGoodElectrons" );
        const auto& NGoodMuons              = tr.getVar<int>( "NGoodMuons" );
        const auto& JetID                   = tr.getVar<bool>( "JetID" );
        const auto& NVtx                    = tr.getVar<int>( "NVtx" );

        bool passHTcut                      = ( HT_NonIsoMuon_pt30 > 300 );
        bool correctMuonDataset             = ( runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos );
        bool passNonIsoMuonCut              = ( NNonIsoMuons == 1 );
        bool requireNoGoodMuonsCut          = ( NGoodMuons == 0 );
        bool requireNoGoodElectronsCut      = ( NGoodElectrons == 0 );
        bool require7JetCut                 = ( NNonIsoMuonJets_pt30 >= 7 );

        if( correct2018Split && correctMuonDataset ) {
            for( auto& kv : cut_map ) {
                if( kv.second ) {
                    my_histos["h_cutflow_"+kv.first]->Fill( 0.5, eventweight );
                    if( JetID ) my_histos["h_cutflow_"+kv.first]->Fill( 1.5, eventweight );
                    if( JetID && passHTcut ) my_histos["h_cutflow_"+kv.first]->Fill( 2.5, eventweight );
                    if( JetID && passHTcut && passNonIsoMuonTriggers ) my_histos["h_cutflow_"+kv.first]->Fill( 3.5, eventweight );
                    if( JetID && passHTcut && passNonIsoMuonTriggers && passNonIsoMuonCut ) my_histos["h_cutflow_"+kv.first]->Fill( 4.5, eventweight );
                    if( JetID && passHTcut && passNonIsoMuonTriggers && passNonIsoMuonCut && requireNoGoodMuonsCut ) my_histos["h_cutflow_"+kv.first]->Fill( 5.5, eventweight );
                    if( JetID && passHTcut && passNonIsoMuonTriggers && passNonIsoMuonCut && requireNoGoodMuonsCut && requireNoGoodElectronsCut ) my_histos["h_cutflow_"+kv.first]->Fill( 6.5, eventweight );
                    if( JetID && passHTcut && passNonIsoMuonTriggers && passNonIsoMuonCut && requireNoGoodMuonsCut && requireNoGoodElectronsCut && require7JetCut ) my_histos["h_cutflow_"+kv.first]->Fill( 7.5, eventweight );
                    if( JetID && passHTcut && passNonIsoMuonTriggers && passNonIsoMuonCut && requireNoGoodMuonsCut && requireNoGoodElectronsCut && require7JetCut ) my_histos["h_cutflow_"+kv.first]->Fill( 8.5, eventweight );
                    if( JetID && passHTcut && passNonIsoMuonTriggers && passNonIsoMuonCut && requireNoGoodMuonsCut && requireNoGoodElectronsCut && require7JetCut && passHEMVeto ) my_histos["h_cutflow_"+kv.first]->Fill( 9.5, eventweight );
                    if( JetID && passHTcut && passNonIsoMuonTriggers && passNonIsoMuonCut && requireNoGoodMuonsCut && requireNoGoodElectronsCut && require7JetCut && passHEMVeto && passMETFilters ) my_histos["h_cutflow_"+kv.first]->Fill( 10.5, eventweight );
                    
                    my_histos["h_npv_noTrig_"+kv.first]->Fill( NVtx, eventweight );
                    if( passBaseline ) my_histos["h_npv_noTrig_passBaseline_"+kv.first]->Fill( NVtx, eventweight );
                    if( passNonIsoMuonTriggers ) my_histos["h_npv_passTrig_"+kv.first]->Fill( NVtx, eventweight );
                    if( passNonIsoMuonTriggers && passBaseline ) my_histos["h_npv_passTrig_passBaseline_"+kv.first]->Fill( NVtx, eventweight ); 

                    my_histos["h_madHT_"+kv.first]->Fill( madHT, eventweight ); 
                    my_histos["h_tw_madHT_"+kv.first]->Fill( madHT, totalWeightNIM ); 
    
                    my_histos["h_njets_noTrig_"+kv.first]->Fill( NNonIsoMuonJets_pt30, eventweight ); 
                    my_histos["h_tw_njets_noTrig_"+kv.first]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM ); 
                    my_histos["h_ht_noTrig_"+kv.first]->Fill( HT_NonIsoMuon_pt30, eventweight ); 
                    my_histos["h_tw_ht_noTrig_"+kv.first]->Fill( HT_NonIsoMuon_pt30, totalWeightNIM ); 
                    if( passBaseline ) {
                        my_histos["h_njets_noTrig_passBaseline_"+kv.first]->Fill( NNonIsoMuonJets_pt30, eventweight ); 
                        my_histos["h_tw_njets_noTrig_passBaseline_"+kv.first]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM ); 
                        my_histos["h_ht_noTrig_passBaseline_"+kv.first]->Fill( HT_NonIsoMuon_pt30, eventweight ); 
                        my_histos["h_tw_ht_noTrig_passBaseline_"+kv.first]->Fill( HT_NonIsoMuon_pt30, totalWeightNIM ); 
                    }
    
                if( passNonIsoMuonTriggers ) {
                        my_histos["h_njets_passTrig_"+kv.first]->Fill( NNonIsoMuonJets_pt30, eventweight );
                        my_histos["h_tw_njets_passTrig_"+kv.first]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                        my_histos["h_ht_passTrig_"+kv.first]->Fill( HT_NonIsoMuon_pt30, eventweight );
                        my_histos["h_tw_ht_passTrig_"+kv.first]->Fill( HT_NonIsoMuon_pt30, totalWeightNIM );
                        if( passBaseline ) {
                            my_histos["h_njets_passTrig_passBaseline_"+kv.first]->Fill( NNonIsoMuonJets_pt30, eventweight );
                            my_histos["h_tw_njets_passTrig_passBaseline_"+kv.first]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                            my_histos["h_ht_passTrig_passBaseline_"+kv.first]->Fill( HT_NonIsoMuon_pt30, eventweight );
                            my_histos["h_tw_ht_passTrig_passBaseline_"+kv.first]->Fill( HT_NonIsoMuon_pt30, totalWeightNIM );
                        }
                    }
                }
            }
        }
    }
}

void MakeQCDChecks::WriteHistos(TFile* outfile)
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


//This is just the general trigger passing function
bool MakeQCDChecks::passTriggerGeneral( std::vector<std::string>& myTriggerVector, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass ) {
    bool passTrigger = false;
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) { 
        if( TriggerPass.at(i) != 1 ) continue;
        std::string trigname = TriggerNames.at(i);
    
        //Now comes the fun bit logic
        if( std::any_of( myTriggerVector.begin(), myTriggerVector.end(), [&] (std::string s) { return trigname.find(s) != std::string::npos; } ) ) { 
            passTrigger = true;
            break;
        }   
    }   
    return passTrigger;
}
