#define ExtraChecksv2_cxx
#include "Analyzer/Analyzer/include/ExtraChecksv2.h"
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

ExtraChecksv2::ExtraChecksv2()
{
    InitHistos();
}

//Define all your histograms here. 
void ExtraChecksv2::InitHistos()
{
    std::vector<std::string> leptonTags { "nonIsoMuon", "isoMuon", "isoElectron", "isoAll" };
    std::vector<std::string> mvaTags    { "D1", "D2", "D3", "D4" };

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    for( std::string leptonTag : leptonTags ) {
        for( std::string mvaTag : mvaTags ) {
            my_histos.emplace( "h_njets_"+mvaTag+"_"+leptonTag, std::make_shared<TH1D>( ( "h_njets_"+mvaTag+"_"+leptonTag ).c_str(), ( "h_njets_"+mvaTag+"_"+leptonTag ).c_str(), 6, 7, 13 ) );
            my_histos.emplace( "h_njets_"+mvaTag+"_"+leptonTag+"_pt55", std::make_shared<TH1D>( ( "h_njets_"+mvaTag+"_"+leptonTag+"_pt55" ).c_str(), ( "h_njets_"+mvaTag+"_"+leptonTag+"_pt55" ).c_str(), 6, 7, 13 ) );
        }
    }

    my_histos.emplace( "h_nEvtsPerMvaBin", std::make_shared<TH1D>( "h_nEvtsPerMvaBin", "h_nEvtsPerMvaBin", 4, 0, 4 ) );
    my_histos.emplace( "h_nEvtsPerMvaBin_NIM", std::make_shared<TH1D>( "h_nEvtsPerMvaBin_NIM", "h_nEvtsPerMvaBin_NIM", 4, 0, 4 ) );

}

//Put everything you want to do per event here.
void ExtraChecksv2::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    std::vector<std::string> leptonTags { "nonIsoMuon", "isoMuon", "isoElectron", "isoAll" };
    std::vector<std::string> mvaTags    { "D1", "D2", "D3", "D4" };
    
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
        
        const auto& passBaseline1l          = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline1l_NIM      = tr.getVar<bool>("passBaseline1l_NonIsoMuon");

        //const auto& NGoodJets_pt30          = tr.getVar<int>("NGoodJets_pt30");
        //const auto& NNonIsoMuonJets_pt30    = tr.getVar<int>("NNonIsoMuonJets_pt30");

        //const auto& GoodLeptons             = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        //const auto& GoodNonIsoMuons         = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodNonIsoMuons");

        const auto& deepESM_bin1            = tr.getVar<bool>("deepESM_bin1");
        const auto& deepESM_bin2            = tr.getVar<bool>("deepESM_bin2");
        const auto& deepESM_bin3            = tr.getVar<bool>("deepESM_bin3");
        const auto& deepESM_bin4            = tr.getVar<bool>("deepESM_bin4");
        
        const auto& deepESM_binNonIsoMuon1  = tr.getVar<bool>("deepESM_binNonIsoMuon1");
        const auto& deepESM_binNonIsoMuon2  = tr.getVar<bool>("deepESM_binNonIsoMuon2");
        const auto& deepESM_binNonIsoMuon3  = tr.getVar<bool>("deepESM_binNonIsoMuon3");
        const auto& deepESM_binNonIsoMuon4  = tr.getVar<bool>("deepESM_binNonIsoMuon4");

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;
        if( !isQuiet ) { std::cout<<weight<<std::endl; }

        // ------------------------
        // -- Define weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
        // ------------------------
        double totalWeight          = 1.0;
        double totalWeightNIM       = 1.0;

        if(runtype == "MC") {
            const auto& totalEventWeight = tr.getVar<double>("totalEventWeight");
            const auto& totalEventWeightNIM = tr.getVar<double>("totalEventWeightNIM");
         
            totalWeight = totalEventWeight*lumi;
            totalWeightNIM = totalEventWeightNIM*lumi;
        }

        if( passBaseline1l ) {
            if( deepESM_bin1 ) my_histos[ "h_nEvtsPerMvaBin" ]->Fill( 0.5, totalWeight );
            if( deepESM_bin2 ) my_histos[ "h_nEvtsPerMvaBin" ]->Fill( 1.5, totalWeight );
            if( deepESM_bin3 ) my_histos[ "h_nEvtsPerMvaBin" ]->Fill( 2.5, totalWeight );
            if( deepESM_bin4 ) my_histos[ "h_nEvtsPerMvaBin" ]->Fill( 3.5, totalWeight );
        }
        
        if( passBaseline1l_NIM ) {
            if( deepESM_binNonIsoMuon1 ) my_histos[ "h_nEvtsPerMvaBin_NIM" ]->Fill( 0.5, totalWeightNIM );
            if( deepESM_binNonIsoMuon2 ) my_histos[ "h_nEvtsPerMvaBin_NIM" ]->Fill( 1.5, totalWeightNIM );
            if( deepESM_binNonIsoMuon3 ) my_histos[ "h_nEvtsPerMvaBin_NIM" ]->Fill( 2.5, totalWeightNIM );
            if( deepESM_binNonIsoMuon4 ) my_histos[ "h_nEvtsPerMvaBin_NIM" ]->Fill( 3.5, totalWeightNIM );
        }
/*
        if( passBaseline1l_NIM ) {
            if( NNonIsoMuonJets_pt30 < 12 ) {
                if( deepESM_binNonIsoMuon1 ) {
                    my_histos["h_njets_D1_nonIsoMuon"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                    if( GoodNonIsoMuons[0].second.Pt() > 55.0 ) my_histos["h_njets_D1_nonIsoMuon_pt55"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                }
                if( deepESM_binNonIsoMuon2 ) {
                    my_histos["h_njets_D2_nonIsoMuon"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                    if( GoodNonIsoMuons[0].second.Pt() > 55.0 ) my_histos["h_njets_D2_nonIsoMuon_pt55"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                }
                if( deepESM_binNonIsoMuon3 ) {
                    my_histos["h_njets_D3_nonIsoMuon"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                    if( GoodNonIsoMuons[0].second.Pt() > 55.0 ) my_histos["h_njets_D3_nonIsoMuon_pt55"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                }
                if( deepESM_binNonIsoMuon4 ) {
                    my_histos["h_njets_D4_nonIsoMuon"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                    if( GoodNonIsoMuons[0].second.Pt() > 55.0 ) my_histos["h_njets_D4_nonIsoMuon_pt55"]->Fill( NNonIsoMuonJets_pt30, totalWeightNIM );
                }
            }
            else {
                if( deepESM_binNonIsoMuon1 ) {
                    my_histos["h_njets_D1_nonIsoMuon"]->Fill( 12, totalWeightNIM );
                    if( GoodNonIsoMuons[0].second.Pt() > 55.0 ) my_histos["h_njets_D1_nonIsoMuon_pt55"]->Fill( 12, totalWeightNIM );
                }
                if( deepESM_binNonIsoMuon2 ) {
                    my_histos["h_njets_D2_nonIsoMuon"]->Fill( 12, totalWeightNIM );
                    if( GoodNonIsoMuons[0].second.Pt() > 55.0 ) my_histos["h_njets_D2_nonIsoMuon_pt55"]->Fill( 12, totalWeightNIM );
                }
                if( deepESM_binNonIsoMuon3 ) {
                    my_histos["h_njets_D3_nonIsoMuon"]->Fill( 12, totalWeightNIM );
                    if( GoodNonIsoMuons[0].second.Pt() > 55.0 ) my_histos["h_njets_D3_nonIsoMuon_pt55"]->Fill( 12, totalWeightNIM );
                }
                if( deepESM_binNonIsoMuon4 ) {
                    my_histos["h_njets_D4_nonIsoMuon"]->Fill( 12, totalWeightNIM );
                    if( GoodNonIsoMuons[0].second.Pt() > 55.0 ) my_histos["h_njets_D4_nonIsoMuon_pt55"]->Fill( 12, totalWeightNIM );
                }
            }
        }
        
        if( passBaseline1l ) {
            if( NGoodJets_pt30 < 12 ) {
                if( deepESM_bin1 ) {
                    my_histos["h_njets_D1_isoAll"]->Fill( NGoodJets_pt30, totalWeight );
                    if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D1_isoAll_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                }
                if( deepESM_bin2 ) {
                    my_histos["h_njets_D2_isoAll"]->Fill( NGoodJets_pt30, totalWeight );
                    if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D2_isoAll_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                }
                if( deepESM_bin3 ) {
                    my_histos["h_njets_D3_isoAll"]->Fill( NGoodJets_pt30, totalWeight );
                    if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D3_isoAll_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                }
                if( deepESM_bin4 ) {
                    my_histos["h_njets_D4_isoAll"]->Fill( NGoodJets_pt30, totalWeight );
                    if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D4_isoAll_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                }
            }
            else {
                if( deepESM_bin1 ) {
                    my_histos["h_njets_D1_isoAll"]->Fill( 12, totalWeight );
                    if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D1_isoAll_pt55"]->Fill( 12, totalWeight );
                }
                if( deepESM_bin2 ) {
                    my_histos["h_njets_D2_isoAll"]->Fill( 12, totalWeight );
                    if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D2_isoAll_pt55"]->Fill( 12, totalWeight );
                }
                if( deepESM_bin3 ) {
                    my_histos["h_njets_D3_isoAll"]->Fill( 12, totalWeight );
                    if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D3_isoAll_pt55"]->Fill( 12, totalWeight );
                }
                if( deepESM_bin4 ) {
                    my_histos["h_njets_D4_isoAll"]->Fill( 12, totalWeight );
                    if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D4_isoAll_pt55"]->Fill( 12, totalWeight );
                }
            }

            if( GoodLeptons[0].first == "e" ) {
                if( NGoodJets_pt30 < 12 ) {
                    if( deepESM_bin1 ) {
                        my_histos["h_njets_D1_isoElectron"]->Fill( NGoodJets_pt30, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D1_isoElectron_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                    if( deepESM_bin2 ) {
                        my_histos["h_njets_D2_isoElectron"]->Fill( NGoodJets_pt30, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D2_isoElectron_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                    if( deepESM_bin3 ) {
                        my_histos["h_njets_D3_isoElectron"]->Fill( NGoodJets_pt30, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D3_isoElectron_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                    if( deepESM_bin4 ) {
                        my_histos["h_njets_D4_isoElectron"]->Fill( NGoodJets_pt30, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D4_isoElectron_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                }
                else {
                    if( deepESM_bin1 ) {
                        my_histos["h_njets_D1_isoElectron"]->Fill( 12, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D1_isoElectron_pt55"]->Fill( 12, totalWeight );
                    }
                    if( deepESM_bin2 ) {
                        my_histos["h_njets_D2_isoElectron"]->Fill( 12, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D2_isoElectron_pt55"]->Fill( 12, totalWeight );
                    }
                    if( deepESM_bin3 ) {
                        my_histos["h_njets_D3_isoElectron"]->Fill( 12, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D3_isoElectron_pt55"]->Fill( 12, totalWeight );
                    }
                    if( deepESM_bin4 ) {
                        my_histos["h_njets_D4_isoElectron"]->Fill( 12, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D4_isoElectron_pt55"]->Fill( 12, totalWeight );
                    }
                }
            }
            
            if( GoodLeptons[0].first == "m" ) {
                if( NGoodJets_pt30 < 12 ) {
                    if( deepESM_bin1 ) {
                        my_histos["h_njets_D1_isoMuon"]->Fill( NGoodJets_pt30, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D1_isoMuon_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                    if( deepESM_bin2 ) {
                        my_histos["h_njets_D2_isoMuon"]->Fill( NGoodJets_pt30, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D2_isoMuon_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                    if( deepESM_bin3 ) {
                        my_histos["h_njets_D3_isoMuon"]->Fill( NGoodJets_pt30, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D3_isoMuon_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                    if( deepESM_bin4 ) {
                        my_histos["h_njets_D4_isoMuon"]->Fill( NGoodJets_pt30, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D4_isoMuon_pt55"]->Fill( NGoodJets_pt30, totalWeight );
                    }
                }
                else {
                    if( deepESM_bin1 ) {
                        my_histos["h_njets_D1_isoMuon"]->Fill( 12, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D1_isoMuon_pt55"]->Fill( 12, totalWeight );
                    }
                    if( deepESM_bin2 ) {
                        my_histos["h_njets_D2_isoMuon"]->Fill( 12, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D2_isoMuon_pt55"]->Fill( 12, totalWeight );
                    }
                    if( deepESM_bin3 ) {
                        my_histos["h_njets_D3_isoMuon"]->Fill( 12, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D3_isoMuon_pt55"]->Fill( 12, totalWeight );
                    }
                    if( deepESM_bin4 ) {
                        my_histos["h_njets_D4_isoMuon"]->Fill( 12, totalWeight );
                        if( GoodLeptons[0].second.Pt() > 55.0 ) my_histos["h_njets_D4_isoMuon_pt55"]->Fill( 12, totalWeight );
                    }
                }
            }
        }
*/
    } 
}

void ExtraChecksv2::WriteHistos(TFile* outfile)
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
