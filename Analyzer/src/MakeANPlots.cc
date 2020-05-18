#define MakeANPlots_cxx
#include "Analyzer/Analyzer/include/MakeANPlots.h"
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

MakeANPlots::MakeANPlots()
{
    InitHistos();
}

//Define all your histograms here. 
void MakeANPlots::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    std::vector<std::string> njetTags = {"7", "8", "9", "10", "11", "12" };
    //Define 1D histograms
    
    //Section 6.1 - looser njets distribution
    my_histos.emplace( "h_njets_all", std::make_shared<TH1D>( "h_njets_all", "h_njets_all", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_all_noHT", std::make_shared<TH1D>( "h_njets_all_noHT", "h_njets_all_noHT", 20, 0, 20 ) );
    my_histos.emplace( "h_mva_all", std::make_shared<TH1D>( "h_mva_all", "h_mva_all", 20, 0.0, 1.0 ) );

    //Section 6.2 - MVA and HT output per Njet
    for( std::string njetTag : njetTags ) {
        my_histos.emplace( "h_mva_"+njetTag+"j", std::make_shared<TH1D>( ( "h_mva_"+njetTag+"j" ).c_str(), ( "h_mva_"+njetTag+"j" ).c_str(), 20, 0.0, 1.0 ) );
        my_histos.emplace( "h_ht_"+njetTag+"j", std::make_shared<TH1D>( ( "h_ht_"+njetTag+"j" ).c_str(), ( "h_ht_"+njetTag+"j" ).c_str(), 30, 0.0, 3000.0 ) );
        
        my_histos.emplace( "h_mva_"+njetTag+"j_noHT", std::make_shared<TH1D>( ( "h_mva_"+njetTag+"j_noHT" ).c_str(), ( "h_mva_"+njetTag+"j_noHT" ).c_str(), 20, 0.0, 1.0 ) );
        my_histos.emplace( "h_ht_"+njetTag+"j_noHT", std::make_shared<TH1D>( ( "h_ht_"+njetTag+"j_noHT" ).c_str(), ( "h_ht_"+njetTag+"j_noHT" ).c_str(), 30, 0.0, 3000.0 ) );
    }

    //Section 6.3 - QCD CR Njets related plots
    my_histos.emplace( "h_njets_qcdcr_all", std::make_shared<TH1D>( "h_njets_qcdcr_all", "h_njets_qcdcr_all", 20, 0, 20 ) );
    my_histos.emplace( "h_mva_qcdcr_all", std::make_shared<TH1D>( "h_mva_qcdcr_all", "h_mva_qcdcr_all", 20, 0.0, 1.0 ) );
    my_histos.emplace( "h_njets_qcdcr_ge7j", std::make_shared<TH1D>( "h_njets_qcdcr_ge7j", "h_njets_qcdcr_ge7j", 20, 0, 20 ) );
    my_histos.emplace( "h_mva_qcdcr_ge7j", std::make_shared<TH1D>( "h_mva_qcdcr_ge7j", "h_mva_qcdcr_ge7j", 20, 0.0, 1.0 ) );
    my_histos.emplace( "h_njets_qcdcr_all_noHT", std::make_shared<TH1D>( "h_njets_qcdcr_all_noHT", "h_njets_qcdcr_all_noHT", 20, 0, 20 ) );
    my_histos.emplace( "h_mva_qcdcr_all_noHT", std::make_shared<TH1D>( "h_mva_qcdcr_all_noHT", "h_mva_qcdcr_all_noHT", 20, 0.0, 1.0 ) );
    my_histos.emplace( "h_njets_qcdcr_ge7j_noHT", std::make_shared<TH1D>( "h_njets_qcdcr_ge7j_noHT", "h_njets_qcdcr_ge7j_noHT", 20, 0, 20 ) );
    my_histos.emplace( "h_mva_qcdcr_ge7j_noHT", std::make_shared<TH1D>( "h_mva_qcdcr_ge7j_noHT", "h_mva_qcdcr_ge7j_noHT", 20, 0.0, 1.0 ) );

    //Section 6.4 - Alternate tt bar dominated regions
    my_histos.emplace( "h_njets_ge2b_le7j", std::make_shared<TH1D>( "h_njets_2b_le7j", "h_njets_2b_le7j", 20, 0, 20 ) );
    my_histos.emplace( "h_mbl_ge2b_le7j", std::make_shared<TH1D>( "h_mbl_2b_le7j", "h_mbl_2b_le7j", 30, 0, 300 ) );
    my_histos.emplace( "h_ht_ge2b_le7j", std::make_shared<TH1D>( "h_ht_2b_le7j", "h_ht_2b_le7j", 50, 0, 5000 ) );
    my_histos.emplace( "h_nb_ge2b_le7j", std::make_shared<TH1D>( "h_nb_2b_le7j", "h_nb_2b_le7j", 10, 0, 10 ) );
    
    //Section 6.5 - Exactly 7 jet region
    my_histos.emplace( "h_mbl_7j", std::make_shared<TH1D>( "h_mbl_7j", "h_mbl_7j", 30, 0, 300 ) );
    my_histos.emplace( "h_mvabin_7j", std::make_shared<TH1D>( "h_mvabin_7j", "h_mvabin_7j", 4, 0, 4 ) );
    my_histos.emplace( "h_nb_7j", std::make_shared<TH1D>( "h_nb_7j", "h_nb_7j", 10, 0, 10 ) );
    my_histos.emplace( "h_jetpt_7j", std::make_shared<TH1D>( "h_jetpt_7j", "h_jetpt_7j", 50, 0, 1000 ) );
    my_histos.emplace( "h_jeteta_7j", std::make_shared<TH1D>( "h_jeteta_7j", "h_jeteta_7j", 30, -3.0, 3.0 ) );
    
    my_histos.emplace( "h_ht_7j_D1", std::make_shared<TH1D>( "h_ht_7j_D1", "h_ht_7j_D1", 30, 0, 3000 ) );
    my_histos.emplace( "h_njets_7j_D1", std::make_shared<TH1D>( "h_njets_7j_D1", "h_njets_7j_D1", 20, 0, 20 ) );
    my_histos.emplace( "h_mbl_7j_D1", std::make_shared<TH1D>( "h_mbl_7j_D1", "h_mbl_7j_D1", 30, 0, 300 ) );
    my_histos.emplace( "h_nb_7j_D1", std::make_shared<TH1D>( "h_nb_7j_D1", "h_nb_7j_D1", 10, 0, 10 ) );
    my_histos.emplace( "h_jetpt_7j_D1", std::make_shared<TH1D>( "h_jetpt_7j_D1", "h_nb_7j_D1", 50, 0, 1000 ) );
    
    my_histos.emplace( "h_ht_ge7j_D1", std::make_shared<TH1D>( "h_ht_ge7j_D1", "h_ht_ge7j_D1", 30, 0, 3000 ) );
    my_histos.emplace( "h_njets_ge7j_D1", std::make_shared<TH1D>( "h_njets_ge7j_D1", "h_njets_ge7j_D1", 20, 0, 20 ) );
    my_histos.emplace( "h_mbl_ge7j_D1", std::make_shared<TH1D>( "h_mbl_ge7j_D1", "h_mbl_ge7j_D1", 30, 0, 300 ) );
    my_histos.emplace( "h_nb_ge7j_D1", std::make_shared<TH1D>( "h_nb_ge7j_D1", "h_nb_ge7j_D1", 10, 0, 10 ) );
    my_histos.emplace( "h_jetpt_ge7j_D1", std::make_shared<TH1D>( "h_jetpt_ge7j_D1", "h_nb_ge7j_D1", 50, 0, 1000 ) );

    //Section 6.6 - Full baseline selection
    my_histos.emplace( "h_mbl_ge7j", std::make_shared<TH1D>( "h_mbl_ge7j", "h_mbl_ge7j", 30, 0, 300 ) );
    my_histos.emplace( "h_mbl_ge7j_noMblCut", std::make_shared<TH1D>( "h_mbl_ge7j_noMblCut", "h_mbl_ge7j_noMblCut", 30, 0, 300 ) );
    my_histos.emplace( "h_mva_ge7j", std::make_shared<TH1D>( "h_mva_ge7j", "h_mva_ge7j", 20, 0, 1.0 ) );
    my_histos.emplace( "h_mvabin_ge7j", std::make_shared<TH1D>( "h_mvabin_ge7j", "h_mvabin_ge7j", 4, 0, 4 ) );
    my_histos.emplace( "h_ht_ge7j", std::make_shared<TH1D>( "h_ht_ge7j", "h_ht_ge7j", 50, 0, 5000 ) );
    
    my_histos.emplace( "h_njets_ge7j", std::make_shared<TH1D>( "h_njets_ge7j", "h_njets_ge7j", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_ge7j_D2", std::make_shared<TH1D>( "h_njets_ge7j_D2", "h_njets_ge7j_D2", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_ge7j_D3", std::make_shared<TH1D>( "h_njets_ge7j_D3", "h_njets_ge7j_D3", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_ge7j_D4", std::make_shared<TH1D>( "h_njets_ge7j_D4", "h_njets_ge7j_D4", 20, 0, 20 ) );
    
    my_histos.emplace( "h_njets_qcdcr_ge7j_D1", std::make_shared<TH1D>( "h_njets_qcdcr_ge7j_D1", "h_njets_qcdcr_ge7j_D1", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_qcdcr_ge7j_D2", std::make_shared<TH1D>( "h_njets_qcdcr_ge7j_D2", "h_njets_qcdcr_ge7j_D2", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_qcdcr_ge7j_D3", std::make_shared<TH1D>( "h_njets_qcdcr_ge7j_D3", "h_njets_qcdcr_ge7j_D3", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_qcdcr_ge7j_D4", std::make_shared<TH1D>( "h_njets_qcdcr_ge7j_D4", "h_njets_qcdcr_ge7j_D4", 20, 0, 20 ) );

    //Appendix I - Looser Selection
    my_histos.emplace( "h_njets_0b", std::make_shared<TH1D>( "h_njets_0b", "h_njets_0b", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_1l", std::make_shared<TH1D>( "h_njets_1l", "h_njets_1l", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_0b_1l", std::make_shared<TH1D>( "h_njets_0b_1l", "h_njets_0b_1l", 20, 0, 20 ) );
    my_histos.emplace( "h_njets_ge2b_1l", std::make_shared<TH1D>( "h_njets_ge2b_1l", "h_njets_ge2b_1l", 20, 0, 20 ) );
}

//Put everything you want to do per event here.
void MakeANPlots::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
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
        //const auto& GoodLeptons      = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& Mbl             = tr.getVar<double>("Mbl");
        
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30");
        const auto& NGoodJets_pt30  = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30 = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NNonIsoMuonJets_pt30  = tr.getVar<int>("NNonIsoMuonJets_pt30");
        
        //const auto& nVtx                = tr.getVar<int>("NVtx");
        //const auto& passBaseline        = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaselineNIM     = tr.getVar<bool>("passBaseline1l_NonIsoMuon");
        const auto& JetID               = tr.getVar<bool>("JetID");
        const auto& passHEMVeto         = tr.getVar<bool>("passHEMVeto");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& correct2018Split    = tr.getVar<bool>("correct2018Split");
        const auto& passMETFilters      = tr.getVar<bool>("passMETFilters");
        //const auto& met                 = tr.getVar<double>("MET");

        const auto& passTrigger     = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC   = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodMuons      = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons  = tr.getVar<int>("NGoodElectrons");

        const auto& deepESM_val     = tr.getVar<double>("deepESM_val");
        const auto& deepESM_valNonIsoMuon     = tr.getVar<double>("deepESM_valNonIsoMuon");

        const auto& deepESM_bin1    = tr.getVar<bool>("deepESM_bin1");
        const auto& deepESM_bin2    = tr.getVar<bool>("deepESM_bin2");
        const auto& deepESM_bin3    = tr.getVar<bool>("deepESM_bin3");
        const auto& deepESM_bin4    = tr.getVar<bool>("deepESM_bin4");

        const auto& goodJets        = tr.getVec<bool>("GoodJets_pt30");
        const auto& Jets            = tr.getVec<TLorentzVector>("Jets");

        const auto& Jet_pt_1        = tr.getVar<double>("Jet_pt_1");
        const auto& Jet_eta_1       = tr.getVar<double>("Jet_eta_1");

       
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;
        if( !isQuiet ) std::cout<<weight<<std::endl;

        // ------------------------
        // -- Define weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
        // ------------------------
        double nohtweight           = 1.0;
        double qcdcrhtweight        = 1.0;
        double totalWeight          = 1.0;
        double totalWeightNIM       = 1.0;

        if(runtype == "MC")
        {
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& totalEventWeight = tr.getVar<double>("totalEventWeight");
            const auto& totalEventWeightNIM = tr.getVar<double>("totalEventWeightNIM");
            const auto& htDerivedweight = tr.getVar<double>("htDerivedweight");
            
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
            totalWeight = totalEventWeight*lumi;
            nohtweight = totalEventWeight*lumi/htDerivedweight;
            qcdcrhtweight = totalEventWeightNIM*lumi*htDerivedweight;
            totalWeightNIM = totalEventWeightNIM*lumi;

        }

        std::vector<int> nJets       = { 7, 8, 9, 10, 11, 12 };
        std::vector<bool> nJetsBool  = { ( NGoodJets_pt30 == 7 ), ( NGoodJets_pt30 == 8 ),
                                    ( NGoodJets_pt30 == 9 ), ( NGoodJets_pt30 == 10 ),
                                    ( NGoodJets_pt30 == 11 ), ( NGoodJets_pt30 >= 12 ) };
        
        double leadingJetPt  = 0.0;
        for( unsigned int itJet = 0; itJet < Jets.size(); itJet++ ) {
            if( goodJets[itJet] ) {
                leadingJetPt = Jets.at(itJet).Pt();
                break;
            }
        }

        bool passBaselineOffline    = JetID && passHEMVeto && correct2018Split && passMETFilters && HT_trigger_pt30 > 300 && passMadHT  && NGoodJets_pt30 >= 4 && passTriggerMC && passTrigger;
        bool passMuonBaseline       = ( runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos ) && NGoodMuons == 1 && NGoodElectrons == 0;
        bool passElectronBaseline   = ( runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos ) && NGoodElectrons == 1 && NGoodMuons == 0;
        bool passMblCut             = 50 < Mbl && Mbl < 250;

        if( passBaselineOffline ) {
            if( NGoodBJets_pt30 == 0 ) {
                my_histos["h_njets_0b"]->Fill( NGoodJets_pt30, nohtweight );
            }
        }

        if( passBaselineOffline && ( passMuonBaseline || passElectronBaseline )) {
            
            my_histos["h_njets_1l"]->Fill( NGoodJets_pt30, nohtweight );
            
            if( NGoodBJets_pt30 == 0 ) {
                my_histos["h_njets_0b_1l"]->Fill( NGoodJets_pt30, nohtweight );
            }

            if( NGoodBJets_pt30 >= 2 ) {
                my_histos["h_njets_ge2b_1l"]->Fill( NGoodJets_pt30, nohtweight );
            }

            if( NGoodJets_pt30 <= 7 && NGoodBJets_pt30 >= 2 ) {
                my_histos["h_njets_ge2b_le7j"]->Fill( NGoodJets_pt30, nohtweight );
                my_histos["h_ht_ge2b_le7j"]->Fill( HT_trigger_pt30, nohtweight );
                my_histos["h_mbl_ge2b_le7j"]->Fill( Mbl, nohtweight );
                my_histos["h_nb_ge2b_le7j"]->Fill( NGoodBJets_pt30, nohtweight );
            }
        }

        //Make cuts and fill histograms here
        if( passBaselineOffline && ( passMuonBaseline || passElectronBaseline ) && NGoodBJets_pt30 >= 1) {

            my_histos["h_mbl_ge7j_noMblCut"]->Fill( Mbl, totalWeight );
            if( passMblCut ) {
                my_histos["h_njets_all"]->Fill( NGoodJets_pt30, totalWeight );
                my_histos["h_njets_all_noHT"]->Fill( NGoodJets_pt30, nohtweight );
                my_histos["h_mva_all"]->Fill( deepESM_val, totalWeight );
    
                for( unsigned i = 0; i < nJets.size(); i++ ) {
                    if( nJetsBool[i] ) {
                        my_histos["h_mva_"+std::to_string(nJets[i])+"j"]->Fill( deepESM_val, totalWeight );
                        my_histos["h_ht_"+std::to_string(nJets[i])+"j"]->Fill( HT_trigger_pt30, totalWeight );
                        my_histos["h_mva_"+std::to_string(nJets[i])+"j_noHT"]->Fill( deepESM_val, nohtweight );
                        my_histos["h_ht_"+std::to_string(nJets[i])+"j_noHT"]->Fill( HT_trigger_pt30, nohtweight );
                    }
                }
    
                if( nJetsBool[0] ) {
                    if( deepESM_bin1 ) {
                        my_histos["h_mvabin_7j"]->Fill( 0.5, totalWeight );
                        my_histos["h_ht_7j_D1"]->Fill( HT_trigger_pt30, totalWeight );
                        my_histos["h_mbl_7j_D1"]->Fill( Mbl, totalWeight );
                        my_histos["h_njets_7j_D1"]->Fill( NGoodJets_pt30, totalWeight );
                        my_histos["h_nb_7j_D1"]->Fill( NGoodBJets_pt30, totalWeight );
                        my_histos["h_jetpt_7j_D1"]->Fill( leadingJetPt, totalWeight );
                    }
                    if( deepESM_bin2 ) my_histos["h_mvabin_7j"]->Fill( 1.5, totalWeight );
                    if( deepESM_bin3 ) my_histos["h_mvabin_7j"]->Fill( 2.5, totalWeight );
                    if( deepESM_bin4 ) my_histos["h_mvabin_7j"]->Fill( 3.5, totalWeight );
                    
                    my_histos["h_mbl_7j"]->Fill( Mbl, totalWeight );
                    my_histos["h_nb_7j"]->Fill( NGoodBJets_pt30, totalWeight );
                    my_histos["h_jetpt_7j"]->Fill( Jet_pt_1, totalWeight );
                    my_histos["h_jeteta_7j"]->Fill( Jet_eta_1, totalWeight );
                }
    
                if( NGoodJets_pt30 >= 7 ) {
                    my_histos["h_njets_ge7j"]->Fill( NGoodJets_pt30, totalWeight );
                    if( deepESM_bin1 ) {
                        my_histos["h_ht_ge7j_D1"]->Fill( HT_trigger_pt30, totalWeight );
                        my_histos["h_njets_ge7j_D1"]->Fill( NGoodJets_pt30, totalWeight );
                        my_histos["h_mbl_ge7j_D1"]->Fill( Mbl, totalWeight );
                        my_histos["h_nb_ge7j_D1"]->Fill( NGoodBJets_pt30, totalWeight );
                        my_histos["h_jetpt_ge7j_D1"]->Fill( leadingJetPt, totalWeight );
                        my_histos["h_mvabin_ge7j"]->Fill( 0.5, totalWeight );
                    }
                    if( deepESM_bin2 ) {
                        my_histos["h_mvabin_ge7j"]->Fill( 1.5, totalWeight );
                        my_histos["h_njets_ge7j_D2"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                    if( deepESM_bin3 ) {
                        my_histos["h_mvabin_ge7j"]->Fill( 2.5, totalWeight );
                        my_histos["h_njets_ge7j_D3"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                    if( deepESM_bin4 ) {
                        my_histos["h_mvabin_ge7j"]->Fill( 3.5, totalWeight );
                        my_histos["h_njets_ge7j_D4"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                    
                    my_histos["h_mbl_ge7j"]->Fill( Mbl, totalWeight );
                    my_histos["h_mva_ge7j"]->Fill( deepESM_val, totalWeight );
                    my_histos["h_ht_ge7j"]->Fill( HT_trigger_pt30, totalWeight );
                }
            }
        }

        if( passBaselineNIM ) {
            if( filetag.find("Data_SingleElectron") != std::string::npos) continue;
            
            my_histos["h_njets_qcdcr_all"]->Fill( NNonIsoMuonJets_pt30, qcdcrhtweight );
            my_histos["h_mva_qcdcr_all"]->Fill( deepESM_valNonIsoMuon, qcdcrhtweight );
            my_histos["h_njets_qcdcr_all_noHT"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
            my_histos["h_mva_qcdcr_all_noHT"]->Fill( deepESM_valNonIsoMuon, totalWeightNIM );
            
            if( NNonIsoMuonJets_pt30 >= 7 ) {
                my_histos["h_njets_qcdcr_ge7j"]->Fill( NNonIsoMuonJets_pt30, qcdcrhtweight );
                my_histos["h_mva_qcdcr_ge7j"]->Fill( deepESM_valNonIsoMuon, qcdcrhtweight );
                my_histos["h_njets_qcdcr_ge7j_noHT"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                my_histos["h_mva_qcdcr_ge7j_noHT"]->Fill( deepESM_valNonIsoMuon, totalWeightNIM );
                
            }
        }
    } 
}

void MakeANPlots::WriteHistos(TFile* outfile)
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
