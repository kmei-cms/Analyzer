#define AnalyzeLepTrigger_NIM_cxx
#include "Analyzer/Analyzer/include/AnalyzeLepTrigger_NIM.h"
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

AnalyzeLepTrigger_NIM::AnalyzeLepTrigger_NIM()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeLepTrigger_NIM::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define some strings that are used for different scenarios that we want to calculate trigger efficiencies for
    std::vector<std::string> trigTags       { "trig", "noTrig", "justTrig" }; //Trig ( inclusive adding TkMu100) ; noTrig ( nominal ) ; justTrig ( Tk100 but not Mu50 )
    std::vector<std::string> mvaTags        { "All", "D1", "D2", "D3", "D4" }; // NN Bin
 
    for( std::string trigTag : trigTags ) {
        for( std::string mvaTag : mvaTags ) {
            my_histos.emplace( "h_mva_"+trigTag+"_"+mvaTag, std::make_shared<TH1D>( ( "h_mva_"+trigTag+"_"+mvaTag ).c_str(), ( "h_mva_"+trigTag+"_"+mvaTag ).c_str(), 20, 0, 1) );
            my_histos.emplace( "h_njets_"+trigTag+"_"+mvaTag, std::make_shared<TH1D>( ( "h_njets_"+trigTag+"_"+mvaTag ).c_str(), ( "h_njets_"+trigTag+"_"+mvaTag ).c_str(), 6, 7, 13) );
            my_histos.emplace( "h_muonPt_"+trigTag+"_"+mvaTag, std::make_shared<TH1D>( ( "h_muonPt_"+trigTag+"_"+mvaTag ).c_str(), ( "h_muonPt_"+trigTag+"_"+mvaTag ).c_str(), 60, 0, 3000) );
            my_histos.emplace( "h_muonEta_"+trigTag+"_"+mvaTag, std::make_shared<TH1D>( ( "h_muonEta_"+trigTag+"_"+mvaTag ).c_str(), ( "h_muonEta_"+trigTag+"_"+mvaTag ).c_str(), 26, -2.6, 2.6) );
            my_histos.emplace( "h_muonPhi_"+trigTag+"_"+mvaTag, std::make_shared<TH1D>( ( "h_muonPhi_"+trigTag+"_"+mvaTag ).c_str(), ( "h_muonPhi_"+trigTag+"_"+mvaTag ).c_str(), 30, -3.0, 3.0 ) );
        }//End of mvaTags loop
    }//End of trigTags loop
}

//Put everything you want to do per event here.
void AnalyzeLepTrigger_NIM::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter            = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //Define useful variables here
        const auto& runtype                 = tr.getVar<std::string>("runtype");
        const auto& NonIsoMuons             = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodNonIsoMuons");
        
        const auto& passBaseline            = tr.getVar<bool>("passBaseline1l_NonIsoMuon");
        const auto& passBaselineNoTrig      = tr.getVar<bool>("passBaseline1l_NonIsoMuon_noTrigger");

        const auto& deepESM_bin1            = tr.getVar<bool>("deepESM_binNonIsoMuon1");
        const auto& deepESM_bin2            = tr.getVar<bool>("deepESM_binNonIsoMuon2");
        const auto& deepESM_bin3            = tr.getVar<bool>("deepESM_binNonIsoMuon3");
        const auto& deepESM_bin4            = tr.getVar<bool>("deepESM_binNonIsoMuon4");

        const auto& deepESM_valNonIsoMuon   = tr.getVar<double>("deepESM_valNonIsoMuon");
        const auto& NNonIsoMuonJets         = tr.getVar<int>("NNonIsoMuonJets_pt30");
    
        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        //--------------------------------------------------
        //-- Print list of triggers (only if you want to see them)
        //--------------------------------------------------
        const auto& TriggerNames        = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass         = tr.getVec<int>("TriggerPass");
        if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

        // ------------------------
        // -- Define theweight
        // ------------------------
        double eventweight          = 1.0;
        
        if(runtype == "MC") {
            double mcweight         = tr.getVar<double>("totalEventWeightNIM");
            double lumi             = tr.getVar<double>("Lumi");
            eventweight             = mcweight*lumi;
        }

        std::vector<std::string> newTriggers    = { "HLT_TkMu100" };
        bool passNewTrigger                     = PassTriggerGeneral( newTriggers, TriggerNames, TriggerPass );

        if( passBaseline ) {
            my_histos["h_mva_noTrig_All"]->Fill( deepESM_valNonIsoMuon, eventweight );
            my_histos["h_njets_noTrig_All"]->Fill( NNonIsoMuonJets, eventweight );
            my_histos["h_muonPt_noTrig_All"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
            my_histos["h_muonEta_noTrig_All"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
            my_histos["h_muonPhi_noTrig_All"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            if( deepESM_bin1 ) {
                my_histos["h_mva_noTrig_D1"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_noTrig_D1"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_noTrig_D1"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_noTrig_D1"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_noTrig_D1"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
            else if( deepESM_bin2 ) {
                my_histos["h_mva_noTrig_D2"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_noTrig_D2"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_noTrig_D2"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_noTrig_D2"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_noTrig_D2"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
            else if( deepESM_bin3 ) {
                my_histos["h_mva_noTrig_D3"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_noTrig_D3"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_noTrig_D3"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_noTrig_D3"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_noTrig_D3"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
            else if( deepESM_bin4 ) {
                my_histos["h_mva_noTrig_D4"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_noTrig_D4"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_noTrig_D4"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_noTrig_D4"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_noTrig_D4"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
        }

        if( ( passNewTrigger && passBaselineNoTrig ) || passBaseline ) {
            my_histos["h_mva_trig_All"]->Fill( deepESM_valNonIsoMuon, eventweight );
            my_histos["h_njets_trig_All"]->Fill( NNonIsoMuonJets, eventweight );
            my_histos["h_muonPt_trig_All"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
            my_histos["h_muonEta_trig_All"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
            my_histos["h_muonPhi_trig_All"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            if( deepESM_bin1 ) {
                my_histos["h_mva_trig_D1"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_trig_D1"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_trig_D1"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_trig_D1"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_trig_D1"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
            else if( deepESM_bin2 ) {
                my_histos["h_mva_trig_D2"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_trig_D2"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_trig_D2"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_trig_D2"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_trig_D2"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
            else if( deepESM_bin3 ) {
                my_histos["h_mva_trig_D3"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_trig_D3"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_trig_D3"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_trig_D3"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_trig_D3"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
            else if( deepESM_bin4 ) {
                my_histos["h_mva_trig_D4"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_trig_D4"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_trig_D4"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_trig_D4"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_trig_D4"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
        }
        
        if( ( passNewTrigger && passBaselineNoTrig ) && !passBaseline ) {
            my_histos["h_mva_justTrig_All"]->Fill( deepESM_valNonIsoMuon, eventweight );
            my_histos["h_njets_justTrig_All"]->Fill( NNonIsoMuonJets, eventweight );
            my_histos["h_muonPt_justTrig_All"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
            my_histos["h_muonEta_justTrig_All"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
            my_histos["h_muonPhi_justTrig_All"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            if( deepESM_bin1 ) {
                my_histos["h_mva_justTrig_D1"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_justTrig_D1"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_justTrig_D1"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_justTrig_D1"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_justTrig_D1"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
            else if( deepESM_bin2 ) {
                my_histos["h_mva_justTrig_D2"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_justTrig_D2"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_justTrig_D2"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_justTrig_D2"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_justTrig_D2"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
            else if( deepESM_bin3 ) {
                my_histos["h_mva_justTrig_D3"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_justTrig_D3"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_justTrig_D3"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_justTrig_D3"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_justTrig_D3"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
            else if( deepESM_bin4 ) {
                my_histos["h_mva_justTrig_D4"]->Fill( deepESM_valNonIsoMuon, eventweight );
                my_histos["h_njets_justTrig_D4"]->Fill( NNonIsoMuonJets, eventweight );
                my_histos["h_muonPt_justTrig_D4"]->Fill( NonIsoMuons[0].second.Pt(), eventweight );
                my_histos["h_muonEta_justTrig_D4"]->Fill( NonIsoMuons[0].second.Eta(), eventweight );
                my_histos["h_muonPhi_justTrig_D4"]->Fill( NonIsoMuons[0].second.Phi(), eventweight );
            }
        }
    } 
}

void AnalyzeLepTrigger_NIM::WriteHistos(TFile* outfile)
{
    outfile->cd();

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

bool AnalyzeLepTrigger_NIM::containsGoodLepton( const std::vector<TLorentzVector>& leptons, const std::vector<bool>& goodLeptons, double ptThreshold, double etaSelection) { 
    //Require a good muon in the single muon dataset
    for( unsigned int iLep = 0; iLep < leptons.size(); ++iLep ) {
        if( !goodLeptons.at( iLep ) ) continue; 
    
        TLorentzVector myLepton = leptons.at( iLep );
    
        if( myLepton.Pt() >= ptThreshold && std::fabs( myLepton.Eta() ) < etaSelection ) return true;
    }

    return false;
}

int AnalyzeLepTrigger_NIM::goodLeptonIndex( const std::vector<TLorentzVector>& leptons, const std::vector<bool>& goodLeptons) {
    for( unsigned int iLep = 0; iLep < leptons.size(); ++iLep ) {
        if( !goodLeptons.at( iLep ) ) continue;

        return iLep;
    }

    return -1;
}

void AnalyzeLepTrigger_NIM::fillHistos( const std::map<std::string, bool>& cutMap, bool passLeptonTriggers, const TLorentzVector& lepton, double theWeight ) {
    for( auto& kv : cutMap ) {
        if( kv.second ) {
            my_histos["h_trig_den_"+kv.first+"_wLepPtBin"]->Fill( lepton.Pt(), theWeight );
            my_histos["h_trig_den_"+kv.first+"_wLepEtaBin"]->Fill( lepton.Eta(), theWeight );
            my_2d_histos["h2_trig_den_"+kv.first+"_wLepPtLepEtaBin"]->Fill( lepton.Pt(), lepton.Eta(), theWeight );

            if( passLeptonTriggers ) {
                my_histos["h_trig_num_"+kv.first+"_wLepPtBin"]->Fill( lepton.Pt(), theWeight );
                my_histos["h_trig_num_"+kv.first+"_wLepEtaBin"]->Fill( lepton.Eta(), theWeight );
                my_2d_histos["h2_trig_num_"+kv.first+"_wLepPtLepEtaBin"]->Fill( lepton.Pt(), lepton.Eta(), theWeight );
            }
        }
    }
}

//Use this function to print out the entire trigger list in an ntuple (useful when trying to figure out which triggers to use)
void AnalyzeLepTrigger_NIM::printTriggerList( const std::vector<std::string>& TriggerNames ) {
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) {
        std::string myString = TriggerNames.at(i);
        printf("%s\n", myString.c_str());
    }
}

bool AnalyzeLepTrigger_NIM::PassTriggerGeneral(std::vector<std::string>& mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass) {   
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

