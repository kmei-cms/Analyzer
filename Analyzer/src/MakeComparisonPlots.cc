#define MakeComparisonPlots_cxx
#include "Analyzer/Analyzer/include/MakeComparisonPlots.h"
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

MakeComparisonPlots::MakeComparisonPlots()
{
    InitHistos();
}

//Define all your histograms here. 
void MakeComparisonPlots::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define some strings that are used for different scenarios that we want to calculate trigger efficiencies for
    std::vector<std::string> effTags        { "den", "num" }; 
    std::vector<std::string> lepTags        { "el", "mu" }; //Electron, muon, non-iso muon
    std::vector<std::string> ptTags         { "pt35", "pt40" }; //Pt threshold of other lepton
    std::vector<std::string> trigTags       { "trig", "noTrig" }; //Require other lepton trigger? Used for studies, but conclusion was that it does not affect statistics much so it's okay
    std::vector<std::string> nJetCutTags    { "1", "2", "3", "4", "5", "6" }; //Cutting at 6 jets really killed statistics, but the 1 jet cut region is not representative of our signal region, so look at all intermediate steps as well
 
    //Define binning for the histograms
    const Int_t nPtBins = 6;
    Double_t ptBinEdges[ nPtBins + 1 ] =    { 30, 40, 50, 60, 70, 120, 200 };
    const Int_t nEtaBins = 10;
    Double_t etaBinEdges[ nEtaBins + 1 ] =  { -2.4, -2.1, -1.8, -1.4, -1.0, 0, 1.0, 1.4, 1.8, 2.1, 2.4 };

    std::vector<std::string> ptBins         { "30to40", "40to50", "50to60", "60to70", "70to120", "120toInf" };
    std::vector<std::string> etaBins        { "0p0to1p4", "1p4to1p8", "1p8to2p1", "2p1to2p4" };

    for( std::string ptBin : ptBins ) {
        for( std::string etaBin : etaBins ) {
            my_histos.emplace( "h_neutEMfrac_allJets_pt"+ptBin+"_eta"+etaBin, std::make_shared<TH1D>( ( "h_neutEMfrac_allJets_pt"+ptBin+"_eta"+etaBin ).c_str(), ( "h_neutEMfrac_allJets_pt"+ptBin+"_eta"+etaBin ).c_str(), 25, 0.0, 1.0 ) );
            my_histos.emplace( "h_neutHadfrac_allJets_pt"+ptBin+"_eta"+etaBin, std::make_shared<TH1D>( ( "h_neutHadfrac_allJets_pt"+ptBin+"_eta"+etaBin ).c_str(), ( "h_neutHadfrac_allJets_pt"+ptBin+"_eta"+etaBin ).c_str(), 25, 0.0, 1.0 ) );
            my_histos.emplace( "h_neutEMfrac_allJets_pt"+ptBin+"_eta"+etaBin+"_10j", std::make_shared<TH1D>( ( "h_neutEMfrac_allJets_pt"+ptBin+"_eta"+etaBin+"_10j" ).c_str(), ( "h_neutEMfrac_allJets_pt"+ptBin+"_eta"+etaBin+"_10j" ).c_str(), 25, 0.0, 1.0 ) );
            my_histos.emplace( "h_neutHadfrac_allJets_pt"+ptBin+"_eta"+etaBin+"_10j", std::make_shared<TH1D>( ( "h_neutHadfrac_allJets_pt"+ptBin+"_eta"+etaBin+"_10j" ).c_str(), ( "h_neutHadfrac_allJets_pt"+ptBin+"_eta"+etaBin+"_10j" ).c_str(), 25, 0.0, 1.0 ) );
            
            my_histos.emplace( "h_neutEMfrac_leadingJet_pt"+ptBin+"_eta"+etaBin, std::make_shared<TH1D>( ( "h_neutEMfrac_leadingJet_pt"+ptBin+"_eta"+etaBin ).c_str(), ( "h_neutEMfrac_leadingJet_pt"+ptBin+"_eta"+etaBin ).c_str(), 25, 0.0, 1.0 ) );
            my_histos.emplace( "h_neutHadfrac_leadingJet_pt"+ptBin+"_eta"+etaBin, std::make_shared<TH1D>( ( "h_neutHadfrac_leadingJet_pt"+ptBin+"_eta"+etaBin ).c_str(), ( "h_neutHadfrac_leadingJet_pt"+ptBin+"_eta"+etaBin ).c_str(), 25, 0.0, 1.0 ) );
            my_histos.emplace( "h_neutEMfrac_leadingJet_pt"+ptBin+"_eta"+etaBin+"_10j", std::make_shared<TH1D>( ( "h_neutEMfrac_leadingJet_pt"+ptBin+"_eta"+etaBin+"_10j" ).c_str(), ( "h_neutEMfrac_leadingJet_pt"+ptBin+"_eta"+etaBin+"_10j" ).c_str(), 25, 0.0, 1.0 ) );
            my_histos.emplace( "h_neutHadfrac_leadingJet_pt"+ptBin+"_eta"+etaBin+"_10j", std::make_shared<TH1D>( ( "h_neutHadfrac_leadingJet_pt"+ptBin+"_eta"+etaBin+"_10j" ).c_str(), ( "h_neutHadfrac_leadingJet_pt"+ptBin+"_eta"+etaBin+"_10j" ).c_str(), 25, 0.0, 1.0 ) );

            my_histos.emplace( "h_njets_leadingJetPt"+ptBin+"_eta"+etaBin, std::make_shared<TH1D>( ( "h_njets_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), ( "h_njets_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), 8, 7, 15 ) );
            my_histos.emplace( "h_mva_leadingJetPt"+ptBin+"_eta"+etaBin, std::make_shared<TH1D>( ( "h_mva_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), ( "h_mva_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), 25, 0.0, 1.0 ) );
            my_histos.emplace( "h_HT_leadingJetPt"+ptBin+"_eta"+etaBin, std::make_shared<TH1D>( ( "h_HT_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), ( "h_HT_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), 60, 0.0, 3000.0 ) );
            my_histos.emplace( "h_lepPt_leadingJetPt"+ptBin+"_eta"+etaBin, std::make_shared<TH1D>( ( "h_lepPt_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), ( "h_lepPt_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), nPtBins, ptBinEdges ) );
            my_histos.emplace( "h_lepEta_leadingJetPt"+ptBin+"_eta"+etaBin, std::make_shared<TH1D>( ( "h_lepEta_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), ( "h_lepEta_leadingJetPt"+ptBin+"_eta"+etaBin ).c_str(), nEtaBins, etaBinEdges ) );
        }
    }
}

//Put everything you want to do per event here.
void MakeComparisonPlots::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //Define useful variables here
        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& runYear             = tr.getVar<std::string>("runYear");
        const auto& filetag             = tr.getVar<std::string>("filetag");

        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30");

        const auto& Muons               = tr.getVec<TLorentzVector>("Muons");
        const auto& Electrons           = tr.getVec<TLorentzVector>("Electrons");
        const auto& Jets                = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");

        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passBaseline        = tr.getVar<bool>("passBaselineGoodOffline1l");

        const auto& GoodMuons           = tr.getVec<bool>("GoodMuons");
        const auto& GoodElectrons       = tr.getVec<bool>("GoodElectrons");
        const auto& GoodJets_pt30       = tr.getVec<bool>("GoodJets_pt30");

       
        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        //--------------------------------------------------
        //-- Print list of triggers (only if you want to see them)
        //--------------------------------------------------
        //if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

        // ------------------------
        // -- Define theweight
        // ------------------------
        double theweight            = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;

        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            const auto& puWeight = tr.getVar<double>("puWeightCorr");

            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }

            const auto& bTagScaleFactor   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            const auto& prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");

            theweight = lumi*Weight*puWeight*leptonScaleFactor*bTagScaleFactor*prefiringScaleFactor;
            
        }

        if( passBaseline ) {
            
            int leadingJetIt    = 0;
            for( int i = 0; i < Jets.size(); i++ ) {
                if( !GoodJets_pt30.at(i) ) continue;
                
                leadingJetIt = i;
                break;
            }
        
            bool leadJet_pt30to40   = ( Jets.at(leadingJetIt).Pt() < 40.0 );
            bool leadJet_pt40to50   = ( Jets.at(leadingJetIt).Pt() >= 40.0 ) && ( Jets.at(leadingJetIt).Pt() <= 50.0 );
            bool leadJet_pt50to60   = ( Jets.at(leadingJetIt).Pt() >= 50.0 ) && ( Jets.at(leadingJetIt).Pt() <= 60.0 );
            bool leadJet_pt60to70   = ( Jets.at(leadingJetIt).Pt() >= 60.0 ) && ( Jets.at(leadingJetIt).Pt() <= 70.0 );
            bool leadJet_pt70to120   = ( Jets.at(leadingJetIt).Pt() >= 70.0 ) && ( Jets.at(leadingJetIt).Pt() <= 120.0 );
            bool leadJet_pt120toInf   = ( Jets.at(leadingJetIt).Pt() >= 120.0 );

            bool leadJet_absEta0p0to1p4 = ( std::fabs( Jets.at(leadingJetIt).Eta() ) < 1.4 );
            bool leadJet_absEta1p4to1p8 = ( std::fabs( Jets.at(leadingJetIt).Eta() ) >= 1.4 ) && ( std::fabs( Jets.at(leadingJetIt).Eta() ) < 1.8 ) ;
            bool leadJet_absEta1p8to2p1 = ( std::fabs( Jets.at(leadingJetIt).Eta() ) >= 1.8 ) && ( std::fabs( Jets.at(leadingJetIt).Eta() ) < 2.1 ) ;
            bool leadJet_absEta2p1to2p4 = ( std::fabs( Jets.at(leadingJetIt).Eta() ) >= 2.1 ) && ( std::fabs( Jets.at(leadingJetIt).Eta() ) < 2.4 ) ;

            const std::map<std::string, bool> cut_map {
                { "leadingJet_pt30to40_eta0p0to1p4",  leadJet_pt30to40 && leadJet_absEta0p0to1p4 }, 
                { "leadingJet_pt30to40_eta1p4to1p8",  leadJet_pt30to40 && leadJet_absEta1p4to1p8 }, 
                { "leadingJet_pt30to40_eta1p8to2p1",  leadJet_pt30to40 && leadJet_absEta1p8to2p1 }, 
                { "leadingJet_pt30to40_eta2p1to2p4",  leadJet_pt30to40 && leadJet_absEta2p1to2p4 }, 
                { "leadingJet_pt40to50_eta0p0to1p4", leadJet_pt40to50 && leadJet_absEta0p0to1p4 }, 
                { "leadingJet_pt40to50_eta1p4to1p8", leadJet_pt40to50 && leadJet_absEta1p4to1p8 }, 
                { "leadingJet_pt40to50_eta1p8to2p1", leadJet_pt40to50 && leadJet_absEta1p8to2p1 }, 
                { "leadingJet_pt40to50_eta2p1to2p4", leadJet_pt40to50 && leadJet_absEta2p1to2p4 }, 
                { "leadingJet_pt50to60_eta0p0to1p4", leadJet_pt50to60 && leadJet_absEta0p0to1p4 }, 
                { "leadingJet_pt50to60_eta1p4to1p8", leadJet_pt50to60 && leadJet_absEta1p4to1p8 }, 
                { "leadingJet_pt50to60_eta1p8to2p1", leadJet_pt50to60 && leadJet_absEta1p8to2p1 }, 
                { "leadingJet_pt50to60_eta2p1to2p4", leadJet_pt50to60 && leadJet_absEta2p1to2p4 }, 
                { "leadingJet_pt60to70_eta0p0to1p4", leadJet_pt60to70 && leadJet_absEta0p0to1p4 }, 
                { "leadingJet_pt60to70_eta1p4to1p8", leadJet_pt60to70 && leadJet_absEta1p4to1p8 }, 
                { "leadingJet_pt60to70_eta1p8to2p1", leadJet_pt60to70 && leadJet_absEta1p8to2p1 }, 
                { "leadingJet_pt60to70_eta2p1to2p4", leadJet_pt60to70 && leadJet_absEta2p1to2p4 }, 
                { "leadingJet_pt70to120_eta0p0to1p4", leadJet_pt70to120 && leadJet_absEta0p0to1p4 }, 
                { "leadingJet_pt70to120_eta1p4to1p8", leadJet_pt70to120 && leadJet_absEta1p4to1p8 }, 
                { "leadingJet_pt70to120_eta1p8to2p1", leadJet_pt70to120 && leadJet_absEta1p8to2p1 }, 
                { "leadingJet_pt70to120_eta2p1to2p4", leadJet_pt70to120 && leadJet_absEta2p1to2p4 }, 
                { "leadingJet_pt120toInf_eta0p0to1p4", leadJet_pt120toInf && leadJet_absEta0p0to1p4 }, 
                { "leadingJet_pt120toInf_eta1p4to1p8", leadJet_pt120toInf && leadJet_absEta1p4to1p8 }, 
                { "leadingJet_pt120toInf_eta1p8to2p1", leadJet_pt120toInf && leadJet_absEta1p8to2p1 }, 
                { "leadingJet_pt120toInf_eta2p1to2p4", leadJet_pt120toInf && leadJet_absEta2p1to2p4 }
            
            };
   
            const auto& Jets_neutEMfrac         = tr.getVec<double>("Jets_neutralEmEnergyFraction");
            const auto& Jets_neutHadfrac        = tr.getVec<double>("Jets_neutralHadronEnergyFraction");
            const auto& deepESM_val             = tr.getVar<double>("deepESM_val");
            const auto& HT_trigger_pt30         = tr.getVar<double>("HT_trigger_pt30");
            
            bool pass10JetCut                   = ( NGoodJets_pt30 >= 10 );

            for( auto& kv : cut_map ) {
                if( kv.second ) {
                    my_histos["h_neutEMfrac_"+kv.first ]->Fill( Jets_neutEMfrac.at(leadingJetIt), theweight );
                    my_histos["h_neutHadfrac_"+kv.first ]->Fill( Jets_neutHadfrac.at(leadingJetIt), theweight );
                    my_histos["h_njets_"+kv.first ]->Fill( NGoodJets_pt30, theweight );
                    my_histos["h_mva_"+kv.first ]->Fill( deepESM_val, theweight );
                    my_histos["h_HT_"+kv.first ]->Fill( HT_trigger_pt30, theweight );
                    //my_histos["h_lepPt_"+kv.first ]->Fill( GoodLeptons[0].second().Pt(), theweight );
                    //my_histos["h_lepEta_"+kv.first ]->Fill( GoodLeptons[0].second().Eta(), theweight );
                    if( pass10JetCut ) {
                        my_histos["h_neutEMfrac_"+kv.first+"_10j" ]->Fill( Jets_neutEMfrac.at(leadingJetIt), theweight );
                        my_histos["h_neutHadfrac_"+kv.first+"_10j" ]->Fill( Jets_neutHadfrac.at(leadingJetIt), theweight );
                    }
                }
            }
            for( int i = 0; i < Jets.size(); i++ ) {
                
                if( !GoodJets_pt30.at(i) ) continue;

                bool jet_pt30to40   = ( Jets.at(i).Pt() < 40.0 );
                bool jet_pt40to50   = ( Jets.at(i).Pt() >= 40.0 ) && ( Jets.at(i).Pt() <= 50.0 );
                bool jet_pt50to60   = ( Jets.at(i).Pt() >= 50.0 ) && ( Jets.at(i).Pt() <= 60.0 );
                bool jet_pt60to70   = ( Jets.at(i).Pt() >= 60.0 ) && ( Jets.at(i).Pt() <= 70.0 );
                bool jet_pt70to120   = ( Jets.at(i).Pt() >= 70.0 ) && ( Jets.at(i).Pt() <= 120.0 );
                bool jet_pt120toInf   = ( Jets.at(i).Pt() >= 120.0 );
    
                bool jet_absEta0p0to1p4 = ( std::fabs( Jets.at(i).Eta() ) < 1.4 );
                bool jet_absEta1p4to1p8 = ( std::fabs( Jets.at(i).Eta() ) >= 1.4 ) && ( std::fabs( Jets.at(i).Eta() ) < 1.8 ) ;
                bool jet_absEta1p8to2p1 = ( std::fabs( Jets.at(i).Eta() ) >= 1.8 ) && ( std::fabs( Jets.at(i).Eta() ) < 2.1 ) ;
                bool jet_absEta2p1to2p4 = ( std::fabs( Jets.at(i).Eta() ) >= 2.1 ) && ( std::fabs( Jets.at(i).Eta() ) < 2.4 ) ;

                const std::map<std::string, bool> cut_map_jet {
                    { "allJets_pt30to40_eta0p0to1p4",  jet_pt30to40 && jet_absEta0p0to1p4 }, 
                    { "allJets_pt30to40_eta1p4to1p8",  jet_pt30to40 && jet_absEta1p4to1p8 }, 
                    { "allJets_pt30to40_eta1p8to2p1",  jet_pt30to40 && jet_absEta1p8to2p1 }, 
                    { "allJets_pt30to40_eta2p1to2p4",  jet_pt30to40 && jet_absEta2p1to2p4 }, 
                    { "allJets_pt40to50_eta0p0to1p4", jet_pt40to50 && jet_absEta0p0to1p4 }, 
                    { "allJets_pt40to50_eta1p4to1p8", jet_pt40to50 && jet_absEta1p4to1p8 }, 
                    { "allJets_pt40to50_eta1p8to2p1", jet_pt40to50 && jet_absEta1p8to2p1 }, 
                    { "allJets_pt40to50_eta2p1to2p4", jet_pt40to50 && jet_absEta2p1to2p4 }, 
                    { "allJets_pt50to60_eta0p0to1p4", jet_pt50to60 && jet_absEta0p0to1p4 }, 
                    { "allJets_pt50to60_eta1p4to1p8", jet_pt50to60 && jet_absEta1p4to1p8 }, 
                    { "allJets_pt50to60_eta1p8to2p1", jet_pt50to60 && jet_absEta1p8to2p1 }, 
                    { "allJets_pt50to60_eta2p1to2p4", jet_pt50to60 && jet_absEta2p1to2p4 }, 
                    { "allJets_pt60to70_eta0p0to1p4", jet_pt60to70 && jet_absEta0p0to1p4 }, 
                    { "allJets_pt60to70_eta1p4to1p8", jet_pt60to70 && jet_absEta1p4to1p8 }, 
                    { "allJets_pt60to70_eta1p8to2p1", jet_pt60to70 && jet_absEta1p8to2p1 }, 
                    { "allJets_pt60to70_eta2p1to2p4", jet_pt60to70 && jet_absEta2p1to2p4 }, 
                    { "allJets_pt70to120_eta0p0to1p4", jet_pt70to120 && jet_absEta0p0to1p4 }, 
                    { "allJets_pt70to120_eta1p4to1p8", jet_pt70to120 && jet_absEta1p4to1p8 }, 
                    { "allJets_pt70to120_eta1p8to2p1", jet_pt70to120 && jet_absEta1p8to2p1 }, 
                    { "allJets_pt70to120_eta2p1to2p4", jet_pt70to120 && jet_absEta2p1to2p4 }, 
                    { "allJets_pt120toInf_eta0p0to1p4", jet_pt120toInf && jet_absEta0p0to1p4 }, 
                    { "allJets_pt120toInf_eta1p4to1p8", jet_pt120toInf && jet_absEta1p4to1p8 }, 
                    { "allJets_pt120toInf_eta1p8to2p1", jet_pt120toInf && jet_absEta1p8to2p1 }, 
                    { "allJets_pt120toInf_eta2p1to2p4", jet_pt120toInf && jet_absEta2p1to2p4 }, 

                
                };
                for( auto& kv : cut_map_jet ) {
                    if( kv.second ) {
                        my_histos["h_neutEMfrac_"+kv.first ]->Fill( Jets_neutEMfrac.at(i), theweight );
                        my_histos["h_neutHadfrac_"+kv.first ]->Fill( Jets_neutHadfrac.at(i), theweight );
                        if( pass10JetCut ) {
                            my_histos["h_neutEMfrac_"+kv.first+"_10j" ]->Fill( Jets_neutEMfrac.at(i), theweight );
                            my_histos["h_neutHadfrac_"+kv.first+"_10j" ]->Fill( Jets_neutHadfrac.at(i), theweight );
                        }
                    }
                }             
            }
        }
    }
}

void MakeComparisonPlots::WriteHistos(TFile* outfile)
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

