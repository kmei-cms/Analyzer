#define AnalyzeHadTrigger_cxx
#include "Analyzer/Analyzer/include/AnalyzeHadTrigger.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <algorithm>

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "Framework/Framework/include/SetUpTopTagger.h"

AnalyzeHadTrigger::AnalyzeHadTrigger()
{
}


void AnalyzeHadTrigger::InitHistos(NTupleReader &tr)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    std::vector<std::string> num_tags       { "OLS", "OLSwIsoMu", "OLSwIsoMuwSingleBTag", "OLSwIsoMuwDoubleBTag", "OLSwIsoMuwHadTrig", "OLSwIsoMuwCont", "OLSwIsoMuwContwHadTrig", "OLSwRef", "OLSwRefwCont","OLSwIsoMuwContwHadTrigwHT", "OLSwRefwContwHT" };
    //OLS is offline selection, IsoMu is HLT_IsoMu27, Had Trigger is all hadronic triggers, Cont is control path HLT_PFHT430/380_SixPFJet40/32, Ref is HLT_PFHT430/270
    std::vector<std::string> eff_tags       { "l", "t" }; //Loose or tight

    for( std::string numTag : num_tags ) {
        for( std::string effTag : eff_tags ) {
            my_histos.emplace("h_trig_"+numTag+"_"+effTag+"_w_HT_bin", std::make_shared<TH1D>( ( "h_trig_"+numTag+"_"+effTag+"_w_HT_bin" ).c_str(), ( "h_trig_"+numTag+"_"+effTag+"_w_HT_bin" ).c_str(), 40, 500, 2500 ) );
            my_histos.emplace("h_trig_"+numTag+"_"+effTag+"_w_Btag_bin", std::make_shared<TH1D>( ( "h_trig_"+numTag+"_"+effTag+"_w_Btag_bin" ).c_str(), ( "h_trig_"+numTag+"_"+effTag+"_w_Btag_bin" ).c_str(), 6, 2, 8 ) );
            my_histos.emplace("h_trig_"+numTag+"_"+effTag+"_w_JetSixPt_bin", std::make_shared<TH1D>( ( "h_trig_"+numTag+"_"+effTag+"_w_JetSixPt_bin" ).c_str(), ( "h_trig_"+numTag+"_"+effTag+"_w_JetSixPt_bin" ).c_str(), 34, 30, 200 ) );
        }//END of eff_tags
    }//END of num_tags
}//END of Init_Histos

void AnalyzeHadTrigger::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& allJets             = tr.getVec<TLorentzVector>("Jets");
        const auto& allJetsCSV          = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");
        
        const auto& TriggerNames        = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass         = tr.getVec<int>("TriggerPass");

        const auto& HT                  = tr.getVar<double>("HT");
        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        
        const auto& Jets                = tr.getVec<TLorentzVector>("GoodJets_pt45");
        const auto& NJets               = tr.getVar<int>("NGoodJets_pt45");
        const auto& NBJets              = tr.getVar<int>("NGoodBJets_pt45");
        const auto& NBJets_tight        = tr.getVar<int>("NGoodBJets_pt45_tight");
        

        const auto& HT_trigger          = tr.getVar<double>("HT_trigger");

        //------------------------------------
        //-- Print Event Number
        //------------------------------------

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );
        if( tr.getEvtNum() == 1) printTriggerList( TriggerNames );

        //------------------------------------
        //-- Get the Proper Event Weight
        //------------------------------------

        double eventweight = 1.0; //For data

        //------------------------------------
        //-- Calculate the hadronic trigger variables by myself
        //-----------------------------------
        
        std::vector<double> myJetPts_loose;
        std::vector<double> myJetBTagCSVs_loose;
        std::vector<double> myJetPts_tight;
        std::vector<double> myJetBTagCSVs_tight;

        double myHT_pt30 = 0.0;
        double myHT_pt40 = 0.0;
        double myBTagCSVCut = 0.8484;

        myJetPts_loose.clear(); myJetBTagCSVs_loose.clear();
        myJetPts_tight.clear(); myJetBTagCSVs_tight.clear();

        for( unsigned int i = 0; i < allJets.size(); ++i ) {
            
            TLorentzVector myJetLV = allJets.at(i);
            
            if( passLooseSelection( myJetLV ) ) {
                myJetPts_loose.push_back( myJetLV.Pt() );
                myJetBTagCSVs_loose.push_back( allJetsCSV.at(i) );
                myHT_pt30 += myJetLV.Pt();
            }

            if( passTightSelection( myJetLV ) ) {
                myJetPts_tight.push_back( myJetLV.Pt() );
                myJetBTagCSVs_tight.push_back( allJetsCSV.at(i) );
                myHT_pt40 += myJetLV.Pt();
            }

        }

        bool passLoose = true;
        bool passTight = true;
        bool passMuon  = (NGoodMuons > 0);
        
        int  nBTagJets_loose = 0;
        int  nBTagJets_tight = 0;

        if( myJetPts_loose.size() < 6 ) passLoose = false;

        if( passLoose ) {
            
            for( double itJetCSV : myJetBTagCSVs_loose ) {
                if( itJetCSV > myBTagCSVCut ) nBTagJets_loose++;
            }

            if( nBTagJets_loose   <  2.0 ) passLoose = false;
            if( myJetPts_loose[0] < 70.0 ) passLoose = false;
            if( myJetPts_loose[1] < 55.0 ) passLoose = false;
            if( myJetPts_loose[2] < 40.0 ) passLoose = false;
            if( myJetPts_loose[3] < 35.0 ) passLoose = false;
            if( myHT_pt30         < 450. ) passLoose = false;
        }
        
        if( myJetPts_tight.size() < 6 ) passTight = false;

        if( passTight ) {
            
            for( double itJetCSV : myJetBTagCSVs_tight ) {
                if( itJetCSV > myBTagCSVCut ) nBTagJets_tight++;
            }

            if( nBTagJets_tight   <  2.0 ) passLoose = false;
            if( myJetPts_tight[0] < 70.0 ) passLoose = false;
            if( myJetPts_tight[1] < 55.0 ) passLoose = false;
            if( myHT_pt30         < 500. ) passLoose = false;
        }
        
        //---------------------------------------------------
        //-- Decide Whether the Event Passes the Denominator 
        //---------------------------------------------------
        
        std::vector<std::string> myIsoMuTriggers        { "HLT_IsoMu22", "HLT_IsoMu22_eta2p1", "HLT_IsoMu24", "HLT_IsoMu27" };
        std::vector<std::string> mySingleBTagTriggers   { "HLT_PFHT450_SixJet40_BTagCSV_p056"};
        std::vector<std::string> myDoubleBTagTriggers   { "HLT_PFHT400_SixJet30_DoubleBTagCSV_p056" };
        std::vector<std::string> myHadronicTriggers     { "HLT_PFHT450_SixJet40_BTagCSV_p056", "HLT_PFHT400_SixJet30_DoubleBTagCSV_p056" };
        std::vector<std::string> myControlTriggers      { "HLT_PFHT450_SixJet40", "HLT_PFHT400_SixJet30" }; //TODO: NEED TO FIND THE CORRECT CONTROL TRIGGERS
        std::vector<std::string> myReferenceTriggers    { "HLT_PFHT900", "HLT_PFHT1050" };
        std::vector<std::string> myHTTriggers           { "HLT_AK8PFJet450" };

        //Check if it passes the offline selection
        if( passLoose ) {     
            my_histos["h_trig_OLS_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
            my_histos["h_trig_OLS_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
            my_histos["h_trig_OLS_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );

            //Check to see if it passes the offline selection and IsoMu triggers
            if( passTriggerGeneral( myIsoMuTriggers, TriggerNames, TriggerPass ) && NGoodMuons > 0 ) {
                my_histos["h_trig_OLSwIsoMu_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                my_histos["h_trig_OLSwIsoMu_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                my_histos["h_trig_OLSwIsoMu_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                if( passTriggerGeneral( mySingleBTagTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwIsoMuwSingleBTag_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                    my_histos["h_trig_OLSwIsoMuwSingleBTag_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                    my_histos["h_trig_OLSwIsoMuwSingleBTag_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                }
                if( passTriggerGeneral( myDoubleBTagTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwIsoMuwDoubleBTag_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                    my_histos["h_trig_OLSwIsoMuwDoubleBTag_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                    my_histos["h_trig_OLSwIsoMuwDoubleBTag_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                }
                if( passTriggerGeneral( myHadronicTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwIsoMuwHadTrig_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                    my_histos["h_trig_OLSwIsoMuwHadTrig_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                    my_histos["h_trig_OLSwIsoMuwHadTrig_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                }
                if( passTriggerGeneral( myControlTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwIsoMuwCont_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                    my_histos["h_trig_OLSwIsoMuwCont_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                    my_histos["h_trig_OLSwIsoMuwCont_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );

                    if( passTriggerGeneral( myHadronicTriggers, TriggerNames, TriggerPass ) ) {
                        my_histos["h_trig_OLSwIsoMuwContwHadTrig_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrig_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrig_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                    }
                    else if( passTriggerGeneral( myHTTriggers, TriggerNames, TriggerPass ) ) {
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                    }
                }
            }//End of passing IsoMu triggers

            if( passTriggerGeneral( myReferenceTriggers, TriggerNames, TriggerPass ) ) {
                my_histos["h_trig_OLSwRef_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                my_histos["h_trig_OLSwRef_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                my_histos["h_trig_OLSwRef_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );

                if( passTriggerGeneral( myControlTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwRefwCont_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                    my_histos["h_trig_OLSwRefwCont_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                    my_histos["h_trig_OLSwRefwCont_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                }//END of passing Control triggers

                else if( passTriggerGeneral( myHTTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwRefwContwHT_l_w_Btag_bin"]->Fill( nBTagJets_loose, eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_l_w_HT_bin"]->Fill( myHT_pt30, eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_l_w_JetSixPt_bin"]->Fill( myJetPts_loose[5], eventweight );
                }

            }//END of passing Reference triggers
        }//End of passLoose

        if( passTight ) {   
            my_histos["h_trig_OLS_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
            my_histos["h_trig_OLS_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
            my_histos["h_trig_OLS_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
            //Check to see if it passes the offline selection and IsoMu triggers
            if( passTriggerGeneral( myIsoMuTriggers, TriggerNames, TriggerPass ) && NGoodMuons > 0 ) {
                my_histos["h_trig_OLSwIsoMu_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                my_histos["h_trig_OLSwIsoMu_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                my_histos["h_trig_OLSwIsoMu_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                if( passTriggerGeneral( mySingleBTagTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwIsoMuwSingleBTag_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                    my_histos["h_trig_OLSwIsoMuwSingleBTag_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                    my_histos["h_trig_OLSwIsoMuwSingleBTag_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                }
                if( passTriggerGeneral( myDoubleBTagTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwIsoMuwDoubleBTag_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                    my_histos["h_trig_OLSwIsoMuwDoubleBTag_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                    my_histos["h_trig_OLSwIsoMuwDoubleBTag_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                }
                if( passTriggerGeneral( myHadronicTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwIsoMuwHadTrig_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                    my_histos["h_trig_OLSwIsoMuwHadTrig_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                    my_histos["h_trig_OLSwIsoMuwHadTrig_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                }
                if( passTriggerGeneral( myControlTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwIsoMuwCont_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                    my_histos["h_trig_OLSwIsoMuwCont_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                    my_histos["h_trig_OLSwIsoMuwCont_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                    if( passTriggerGeneral( myHadronicTriggers, TriggerNames, TriggerPass ) ) {
                        my_histos["h_trig_OLSwIsoMuwContwHadTrig_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrig_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrig_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                    }
                    else if( passTriggerGeneral( myHTTriggers, TriggerNames, TriggerPass ) ) {
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                        my_histos["h_trig_OLSwIsoMuwContwHadTrigwHT_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                    }
                }
            }//End of passing IsoMu triggers

            if( passTriggerGeneral( myReferenceTriggers, TriggerNames, TriggerPass ) ) {
                my_histos["h_trig_OLSwRef_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                my_histos["h_trig_OLSwRef_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                my_histos["h_trig_OLSwRef_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );

                if( passTriggerGeneral( myControlTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwRefwCont_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                    my_histos["h_trig_OLSwRefwCont_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                    my_histos["h_trig_OLSwRefwCont_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                }//END of passing Control triggers
                else if( passTriggerGeneral( myHTTriggers, TriggerNames, TriggerPass ) ) {
                    my_histos["h_trig_OLSwRefwContwHT_t_w_Btag_bin"]->Fill( nBTagJets_tight, eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_t_w_HT_bin"]->Fill( myHT_pt40, eventweight );
                    my_histos["h_trig_OLSwRefwContwHT_t_w_JetSixPt_bin"]->Fill( myJetPts_tight[5], eventweight );
                }
            }//END of passing Reference triggers
        }//End of passTight
    }//END of while tr.getNextEvent loop   
}//END of function

void AnalyzeHadTrigger::WriteHistos( TFile* outfile ) 
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

bool AnalyzeHadTrigger::passTriggerGeneral( std::vector<std::string> &mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass )
{
    bool passTrigger = false;
    for( unsigned int i = 0; i < TriggerNames.size(); ++i ) {
         
        if( TriggerPass.at(i) != 1 ) continue;       
        std::string trigname = TriggerNames.at(i);
        
        if( std::any_of( mytriggers.begin(), mytriggers.end(), [&] (std::string s) { return trigname.find(s) != std::string::npos; } ) ) {
                    
            passTrigger = true;
            break;
        }
    }
    return passTrigger;
}

bool AnalyzeHadTrigger::passLooseSelection( TLorentzVector& lv ) {
    bool passLoose = true;
    if( lv.Pt() < 30 ) passLoose = false;
    if( std::fabs( lv.Eta() ) < 2.4 ) passLoose = false;
}

bool AnalyzeHadTrigger::passTightSelection( TLorentzVector& lv ) {
    bool passLoose = true;
    if( lv.Pt() < 40 ) passLoose = false;
    if( std::fabs( lv.Eta() ) < 2.4 ) passLoose = false;
}

void AnalyzeHadTrigger::printTriggerList( const std::vector<std::string>& TriggerNames ) {
    for( int i = 0; i < TriggerNames.size(); i++ ) {
        std::string myString = TriggerNames.at(i);
        printf("%s\n",myString.c_str());
    }
}

void AnalyzeHadTrigger::doesTriggerExist( std::vector<std::string> &myTrigTest, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass ) {
    if( passTriggerGeneral( myTrigTest, TriggerNames, TriggerPass) ) printf("%s\n","This Trigger Works"); 
}
