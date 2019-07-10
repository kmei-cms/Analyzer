#define TwoLepAnalyzer_cxx
#include "Analyzer/Analyzer/include/TwoLepAnalyzer.h"
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

TwoLepAnalyzer::TwoLepAnalyzer() : inithisto(false) 
{
}


/// Define histos
void TwoLepAnalyzer::InitHistos(const std::map<std::string, bool>& cutmap)
{
    	TH1::SetDefaultSumw2();
    	TH2::SetDefaultSumw2();

	my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;


    	for (const auto& cutVar : cutmap) 
        {
		my_histos.emplace( "h_ntops_"+cutVar.first, std::make_shared<TH1D> ( ("h_ntops_"+cutVar.first).c_str(), ("h_ntops_"+cutVar.first).c_str(), 10, 0, 10 ) );
        	my_histos.emplace( "h_ht_"+cutVar.first, std::make_shared<TH1D> ( ("h_ht_"+cutVar.first).c_str(), ("h_ht_"+cutVar.first).c_str(), 60, 0, 3000 ) );
        	my_histos.emplace( "h_njets_"+cutVar.first, std::make_shared<TH1D> ( ("h_njets_"+cutVar.first).c_str(), ("h_njets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        	my_histos.emplace( "h_nbjets_"+cutVar.first, std::make_shared<TH1D> ( ("h_nbjets_"+cutVar.first).c_str(), ("h_nbjets_"+cutVar.first).c_str(), 15, 0, 15 ) );
		my_2d_histos.emplace( "h_njets_MVA_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_MVA_"+cutVar.first).c_str(), ("h_njets_MVA_"+cutVar.first).c_str(), 8, 7, 15, 50, 0, 1.0 ) );
                my_histos.emplace( "h_jet_pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_jet_pt_"+cutVar.first).c_str(), ("h_jet_pt_"+cutVar.first).c_str(), 150, 0, 1500 ) );
                my_histos.emplace( "h_LepLepDeltaR_"+cutVar.first, std::make_shared<TH1D> ( ("h_LepLepDeltaR_"+cutVar.first).c_str(), ("h_LepLepDeltaR_"+cutVar.first).c_str(), 50, 0, 10 ) );
                my_histos.emplace( "h_lepton_pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_lepton_pt_"+cutVar.first).c_str(), ("h_lepton_pt_"+cutVar.first).c_str(), 100, 0, 700 ) );
                my_histos.emplace( "h_lepton_mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_lepton_mass_"+cutVar.first).c_str(), ("h_lepton_mass_"+cutVar.first).c_str(), 50, 0, 500 ) );
                my_histos.emplace( "h_Mbl_"+cutVar.first, std::make_shared<TH1D> ( ("h_Mbl_"+cutVar.first).c_str(), ("h_Mbl_"+cutVar.first).c_str(), 100, 0, 400 ) );
                my_histos.emplace( "h_2b_nonbMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_2b_nonbMass_"+cutVar.first).c_str(), ("h_2b_nonbMass_"+cutVar.first).c_str(), 100, 0, 3000) );
                my_histos.emplace( "h_jet_M_"+cutVar.first, std::make_shared<TH1D> ( ("h_jet_M_"+cutVar.first).c_str(), ("h_jet_M_"+cutVar.first).c_str(), 100, 0, 500 ) );
                my_histos.emplace( "h_Mbl1_"+cutVar.first, std::make_shared<TH1D> ( ("h_Mbl1_"+cutVar.first).c_str(), ("h_Mbl1_"+cutVar.first).c_str(), 100, 0, 500 ) );
                my_histos.emplace( "h_Mbl2_"+cutVar.first, std::make_shared<TH1D> ( ("h_Mbl2_"+cutVar.first).c_str(), ("h_Mbl2_"+cutVar.first).c_str(), 100, 0, 500 ) );
                my_histos.emplace( "h_njetsPlus1_"+cutVar.first, std::make_shared<TH1D> ( ("h_njetsPlus1_"+cutVar.first).c_str(), ("h_njetsPlus1_"+cutVar.first).c_str(), 20, 1, 21 ) );
	}

	//Define TEfficiencies if you are doing trigger studies (for proper error bars) or cut flow charts.
    	my_efficiencies.emplace("event_sel_weight", std::make_shared<TEfficiency>("event_sel_weight","event_sel_weight",9,0,9));
}


/// Put everything you want to do per event 
void TwoLepAnalyzer::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& eventCounter        = tr.getVar<int>("eventCounter");        
        const auto& runtype             = tr.getVar<std::string>("runtype");     
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");

        const auto& JetID               = tr.getVar<bool>("JetID");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& passTriggerMC       = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodBJets_pt30     = tr.getVar<int>("NGoodBJets_pt30");
        const auto& Mbl                 = tr.getVar<double>("Mbl");
        const auto& HT_trigger_pt30     = tr.getVar<double>("HT_trigger_pt30");
        const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30");
        const auto& GoodBJets_pt30      = tr.getVec<bool>("GoodBJets_pt30");
        const auto& deepESM_val         = tr.getVar<double>("deepESM_val");
        
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passBaseline        = tr.getVar<bool>("passBaseline1l_Good");
        const auto& ntops               = tr.getVar<int>("ntops");
        const auto& GoodLeptonsCharge   = tr.getVec<int>("GoodLeptonsCharge");
        const auto& Jets                = tr.getVec<TLorentzVector>("Jets");
        const auto& passTrigger         = tr.getVar<bool>("passTrigger");
        const auto& passBlind           = tr.getVar<bool>("passBlindLep_Good");

        const auto& deepESM_bin1        = tr.getVar<bool>("deepESM_bin1");
        const auto& deepESM_bin2        = tr.getVar<bool>("deepESM_bin2");
        const auto& deepESM_bin3        = tr.getVar<bool>("deepESM_bin3");
        const auto& deepESM_bin4        = tr.getVar<bool>("deepESM_bin4");
        const auto& deepESM_binNum      = tr.getVar<int>("deepESM_binNum");
        const auto& onZ                 = tr.getVar<bool>("onZ");
        const auto& GoodJets_pt30       = tr.getVec<bool>("GoodJets_pt30");

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double htDerivedScaleFactor = 1.0;
        double topPtScaleFactor     = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            bTagScaleFactor   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
            topPtScaleFactor = tr.getVar<double>("topPtScaleFactor");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor = tr.getVar<double>("puWeightCorr");
            
            weight *= eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;
        }
        // Create vector of only bjet TLorentzVectors
        std::vector<TLorentzVector> GoodBJetsVec_pt30;
        for (int j=0; j < Jets.size(); j++)
        {
            if (GoodBJets_pt30.at(j)) GoodBJetsVec_pt30.push_back(Jets.at(j));
        }
        std::vector<TLorentzVector> GoodJetsVec_pt30;
        for (int j=0; j < Jets.size(); j++)
        {
            if (GoodJets_pt30.at(j)) GoodJetsVec_pt30.push_back(Jets.at(j));
        }
        // two lepton M(b,l) definition
        TLorentzVector lep1, lep2;
        int  used_jet1, used_jet2;
        double Mbl1=-1, Mbl2=-1;
        double mbl1diff=999, mbl2diff=999;
        if (NGoodLeptons==2)
        {
            if (GoodLeptons.at(0).second.Pt() >= GoodLeptons.at(1).second.Pt())
            {
                lep1 = GoodLeptons.at(0).second;
                lep2 = GoodLeptons.at(1).second;
            }
            else
            {
                lep1 = GoodLeptons.at(1).second;
                lep2 = GoodLeptons.at(0).second;
            }
            if (NGoodBJets_pt30 >= 2)
            {
                for (int b=0; b < NGoodBJets_pt30; b++)
                {
                    double  mbl1 = (lep1 + GoodBJetsVec_pt30.at(b)).M();
                    if( abs(mbl1 - 105) < mbl1diff)
                    {
                        Mbl1 = mbl1;
                        mbl1diff = abs(mbl1 - 105);
                        used_jet1 = b;
                    }
                    double mbl2 = (lep2 + GoodBJetsVec_pt30.at(b)).M();
                    if (abs(mbl2 - 105) < mbl2diff)
                    {
                        Mbl2 = mbl2;
                        mbl2diff = abs(mbl2 - 105);
                        used_jet2 = b;
                    }
                }                
                if (used_jet1 == used_jet2)
                {
                    if (mbl1diff > mbl2diff)
                    {
                        for (int b=0; b < NGoodBJets_pt30; b++)
                        {
                            double mbl1 = (lep1 + GoodBJetsVec_pt30.at(b)).M();
                            if( abs(mbl1 - 105) < mbl1diff && b != used_jet2)
                            {                            
                                Mbl1 = mbl1;
                                mbl1diff = abs(mbl1 - 105);
                            }
                        }
                    }
                    else 
                    {
                        for (int b=0; b < NGoodBJets_pt30; b++)
                        {
                            double mbl2 = (lep2+ GoodBJetsVec_pt30.at(b)).M();
                            if (abs(mbl2 - 105) < mbl2diff && b != used_jet1)
                            {
                                Mbl2 = mbl2;
                                mbl2diff = abs(mbl2 - 105);
                            }
                        }
                    }
                }                        
            }
            else if (NGoodBJets_pt30 == 1)
            {
                if ( abs((lep1 + GoodBJetsVec_pt30.at(0)).M() - 105) <= abs((lep2 + GoodBJetsVec_pt30.at(0)).M() - 105) )
                {
                    Mbl1 = (lep1 + GoodBJetsVec_pt30.at(0)).M();
                    for (int j = 0; j < NGoodJets_pt30; j++)
                    {
                        double mbl2 = (lep2 + GoodJetsVec_pt30.at(j)).M();
                        if (abs(mbl2 - 105) < mbl2diff && !GoodBJets_pt30.at(j))
                        {
                            Mbl2 = mbl2;
                            mbl2diff = abs(mbl2 - 105);
                        }
                    }
                }
                else
                { 
                    Mbl2 = (lep2 + GoodBJetsVec_pt30.at(0)).M();
                    for (int j = 0; j < NGoodJets_pt30; j++)
                    {
                        double mbl1 = (lep1 + GoodJetsVec_pt30.at(j)).M();
                        if (abs(mbl1 - 105) < mbl1diff && !GoodBJets_pt30.at(j))
                        {
                            Mbl1 = mbl1;
                            mbl1diff = abs(mbl1 - 105);
                        }
                    }
                }
            }
            else            
            {
                for (int j=0; j < NGoodJets_pt30; j++)
                {
                    double  mbl1 = (lep1 + GoodJetsVec_pt30.at(j)).M();
                    if( abs(mbl1 - 105) < mbl1diff)
                    {
                        Mbl1 = mbl1;
                        mbl1diff = abs(mbl1 - 105);
                        used_jet1 = j;
                    }
                    double mbl2 = (lep2 + GoodJetsVec_pt30.at(j)).M();
                    if( abs(mbl2 - 105) < mbl2diff)
                    {
                        Mbl2 = mbl2;
                        mbl2diff = abs(mbl2 - 105);
                        used_jet2 = j;
                    } 
                }
                if (used_jet1 == used_jet2)
                {
                    if (mbl1diff >  mbl2diff)
                    {
                        for (int j=0; j < NGoodJets_pt30; j++)
                        {
                            double mbl1 = (lep1 +GoodJetsVec_pt30.at(j)).M();
                            if (abs(mbl1 - 105) < mbl1diff && j != used_jet2)
                            {
                                Mbl1 = mbl1;
                                mbl1diff = abs(mbl1 - 105);
                            }
                        }
                    }
                    else
                    {
                        for (int j=0; j < NGoodJets_pt30; j++)
                        {
                            double mbl2 = (lep2 + GoodJetsVec_pt30.at(j)).M();
                            if (abs(mbl2-105) < mbl2diff && j != used_jet1)
                            {
                                Mbl2 = mbl2;
                                mbl2diff = abs(mbl2-105);
                            }
                        }
                    }
                }
            }
                
        }
        //Some other plotting variables
        TLorentzVector NonBJetsSum;
        if (NGoodLeptons==2 && NGoodBJets_pt30==2) 
        {
            for (unsigned int i = 0; i < Jets.size(); ++i) 
            {
                if(!GoodBJets_pt30.at(i))  NonBJetsSum += Jets.at(i);
            }
        }
        TLorentzVector lepton_sum;
        for(int l = 0; l < GoodLeptons.size(); l++)
        {
            lepton_sum += GoodLeptons.at(l).second;
        }  
        // cutmap booleans
        bool pass_general =  passMadHT && passBlind;
        bool pass_jgeneral = passMadHT && passBlind && JetID;
        bool pass_1l = NGoodLeptons==1;
        bool is_muon =  pass_1l ? GoodLeptons[0].first=="m" : false;
        bool is_elec =  pass_1l ?  GoodLeptons[0].first=="e" : false; 
        bool pass_2l = NGoodLeptons==2;
        bool pass_2l_opc =  pass_2l ? GoodLeptonsCharge[0]!=GoodLeptonsCharge[1] : false;
        bool is_ee = false, is_em = false, is_me = false, is_mm = false;
        double twolepDeltaR = -1;
        if(pass_2l_opc) {
            is_ee = GoodLeptons[0].first=="e" && GoodLeptons[1].first=="e";
            is_em = GoodLeptons[0].first=="e" && GoodLeptons[1].first=="m";
            is_me = GoodLeptons[0].first=="m" && GoodLeptons[1].first=="e";
            is_mm = GoodLeptons[0].first=="m" && GoodLeptons[1].first=="m";
            twolepDeltaR = GoodLeptons[0].second.DeltaR(GoodLeptons[1].second);
        }
        bool pass_Mbl1 = pass_2l_opc ? (25 <= Mbl1 && Mbl1 <= 250) : false;
        bool pass_Mbl2 =  pass_2l_opc ? (25 <= Mbl2 && Mbl2 <= 250) : false;
        bool pass_bothMbl = pass_Mbl1 && pass_Mbl2;
        // cutmap
        const std::map<std::string, bool>& cutmap
	{
            {"", true},
            {"_2l_", pass_general && pass_2l_opc},
                /*  //assorted cuts and combinations to turn on or off
            {"_1l", pass_general &&  pass_1l},
            {"_ge6j",pass_general && JetID &&  NGoodJets_pt30 >= 6},
            {"_1e", pass_general && is_elec},
            {"_1m", pass_general && is_muon},
            {"_2e", pass_general && is_ee},
            {"_1e1m", pass_general && ( is_em || is_me )},
            {"_2m", pass_general &&  is_mm},
            {"_baseline1l", pass_general &&  passBaseline},
            {"_Mblge50le250", pass_general &&  50 <= Mbl && 250 >= Mbl},
            {"_2l_Mblge50le250", pass_general && pass_2l_opc &&  50 <= Mbl && 250 >= Mbl},
            {"_2e_Mblge50le250", pass_general && is_ee  &&  50 <= Mbl && 250 >= Mbl},
            {"_2m_Mblge50le250", pass_general && is_mm  &&  50 <= Mbl && 250 >= Mbl},
            {"_1e1m_Mblge50le250", pass_general && (is_me || is_em)  &&  50 <= Mbl && 250 >= Mbl},
            {"_2l_Mblge50le250_HTge300", pass_general && pass_2l_opc &&  50 <= Mbl && 250 >= Mbl && HT_trigger_pt30 >=300 },
            {"_2l_0b" , pass_jgeneral && pass_2l_opc && NGoodBJets_pt30  == 0},
            {"_ge1b", pass_general && JetID &&  NGoodBJets_pt30 >= 1},
            {"_ge2b", pass_general && JetID && NGoodBJets_pt30 >= 2},
            {"_2e_ge1b", pass_jgeneral && is_ee && JetID && NGoodBJets_pt30 >=1},
            {"_2e_ge2b", pass_jgeneral && is_ee && JetID && NGoodBJets_pt30 >=2},
            {"_2m_ge1b", pass_jgeneral && is_mm  && JetID && NGoodBJets_pt30 >=1},
            {"_2m_ge2b", pass_jgeneral && is_mm && JetID && NGoodBJets_pt30 >=2},
            {"_1e1m_ge1b", pass_jgeneral && (is_me || is_em)  && NGoodBJets_pt30 >=1},
            {"_1e1m_ge2b", pass_jgeneral && (is_me || is_em) && NGoodBJets_pt30 >=2},
            {"_2l_ge1b_HTge300", pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 1 && HT_trigger_pt30 >= 300},
            {"_2l_ge1b_Mblge50le250", pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 1 && 50 <= Mbl && 250 >= Mbl},
            {"_2l_ge1b_Mblge50le250_HTge300", pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 1 && 50 <= Mbl && 250 >= Mbl && HT_trigger_pt30 >= 300},
            {"_2l_ge2b_HTge300", pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 2 && HT_trigger_pt30 >= 300},
            {"_2l_ge2b_Mblge50le250", pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 2 && 50 <= Mbl && 250 >= Mbl},
            {"_2l_ge1b", pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 1},
            {"_2l_ge2b", pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 2},
            {"_2l_ge1b_d1", pass_jgeneral && pass_2l_opc  && NGoodBJets_pt30 >= 1 && deepESM_bin1},
            {"_2l_ge1b_d2", pass_jgeneral && pass_2l_opc  && NGoodBJets_pt30 >= 1 && deepESM_bin2},
            {"_2l_ge1b_d3", pass_jgeneral && pass_2l_opc  && NGoodBJets_pt30 >= 1 && deepESM_bin3},
            {"_2l_ge1b_d4", pass_jgeneral && pass_2l_opc  && NGoodBJets_pt30 >= 1 && deepESM_bin4},
            {"_2l_1j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 1 && NGoodBJets_pt30 >=1},
            {"_2l_2j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 2 && NGoodBJets_pt30 >=1},
            {"_2l_3j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 3 && NGoodBJets_pt30 >=1},
            {"_2l_4j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 4 && NGoodBJets_pt30 >=1},
            {"_2l_5j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 5 && NGoodBJets_pt30 >=1},
            {"_2l_6j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 6 && NGoodBJets_pt30 >=1},
            {"_2l_7j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 7 && NGoodBJets_pt30 >=1},
            {"_2l_8j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 8 && NGoodBJets_pt30 >=1},
            {"_2l_9j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 9 && NGoodBJets_pt30 >=1},
            {"_2l_10j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 10 && NGoodBJets_pt30 >=1},
            {"_2l_11j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 11 && NGoodBJets_pt30 >=1},
            {"_2l_12j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 12 && NGoodBJets_pt30 >=1},
            {"_2l_13j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 13 && NGoodBJets_pt30 >=1},
            {"_2l_14j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 14 && NGoodBJets_pt30 >=1},
            {"_2l_15j_ge1b", pass_jgeneral && pass_2l_opc && NGoodJets_pt30 == 15 && NGoodBJets_pt30 >=1},
            {"_2l_ge1b_ge7j_Mblge50le250_HTge300", pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 1 && NGoodJets_pt30 >= 7 && 50 <= Mbl && 250 >= Mbl && HT_trigger_pt30 >= 300},
            {"_2l_1t", pass_general && pass_2l_opc && ntops == 1},
            {"_2l_ge1t", pass_general && pass_2l_opc && ntops >= 1},
            {"_2l_onZ", pass_general && pass_2l_opc && onZ},
            {"_2e_onZ", pass_general && onZ && is_ee},
            {"_2m_onZ", pass_general && onZ && is_mm},
            {"_2l_onZ_HTge300_Mblge50le250", pass_general && pass_2l_opc && onZ && HT_trigger_pt30 >= 300 && 50 <= Mbl && 250 >= Mbl},
            {"_2l_offZ_HTge300", pass_general && pass_2l_opc && !onZ &&  HT_trigger_pt30 >= 300},
            {"_2l_offZ_Mblge50le250", pass_general && pass_2l_opc && !onZ && 50 <= Mbl && 250 >= Mbl},
            {"_2l_offZ_HTge300_Mblge50le250", pass_general && pass_2l_opc && !onZ && HT_trigger_pt30 >= 300 && 50 <= Mbl && 250 >= Mbl},
            {"_1e1m_offZ", pass_general && (is_em || is_me) && !onZ},*/
            {"_2l_offZ", pass_general && pass_2l_opc && !onZ},
            //njetcuts
            {"_2l_offZ_ge1j", pass_jgeneral && pass_2l_opc && !onZ && NGoodJets_pt30 >= 1},
            {"_2l_offZ_ge2j", pass_jgeneral && pass_2l_opc && !onZ && NGoodJets_pt30 >= 2},
            {"_2l_offZ_ge3j", pass_jgeneral && pass_2l_opc && !onZ && NGoodJets_pt30 >= 3},
            {"_2l_offZ_ge4j", pass_jgeneral && pass_2l_opc && !onZ && NGoodJets_pt30 >= 4},
            {"_2l_offZ_ge5j", pass_jgeneral && pass_2l_opc && !onZ && NGoodJets_pt30 >= 5},
            {"_2l_offZ_ge6j", pass_jgeneral && pass_2l_opc && !onZ && NGoodJets_pt30 >= 6},
            {"_2l_offZ_ge7j", pass_jgeneral && pass_2l_opc && !onZ && NGoodJets_pt30 >= 7},
                // bcuts
            {"_2l_offZ_ge1b", pass_jgeneral && pass_2l_opc && !onZ && NGoodBJets_pt30 >= 1},
            {"_2l_offZ_ge2b", pass_jgeneral && pass_2l_opc && !onZ && NGoodBJets_pt30 >= 2},
                //Adjusted Mbl cuts
            {"_2l_offZ_Mbl1ge25le250", pass_general && pass_2l_opc && !onZ && pass_Mbl1},
            {"_2l_offZ_Mbl2ge25le250", pass_general && pass_2l_opc && !onZ && pass_Mbl2},
            {"_2l_offZ_BothMblge25le250", pass_general && pass_2l_opc && !onZ && pass_bothMbl},
                // HT cuts
            {"_2l_HTge200", pass_general && pass_2l_opc && HT_trigger_pt30 >= 200},
            {"_2l_offZ_HTge200", pass_general && pass_2l_opc && !onZ && HT_trigger_pt30 >=200},
            // combo cuts
            {"_2l_offZ_ge1b_ge4j", pass_jgeneral && pass_2l_opc && !onZ && NGoodBJets_pt30 >= 1 && NGoodJets_pt30 >= 4},
            {"_2l_offZ_ge2b_ge4j", pass_jgeneral && pass_2l_opc && !onZ && NGoodBJets_pt30 >= 2 && NGoodJets_pt30 >= 4},
            {"_2l_offZ_ge1b_BothMblge25le250", pass_jgeneral && pass_2l_opc && !onZ && NGoodBJets_pt30 >= 1 && pass_bothMbl},
            {"_2l_offZ_ge2b_BothMblge25le250", pass_jgeneral && pass_2l_opc && !onZ && NGoodBJets_pt30  >= 2 && pass_bothMbl},
            {"_2l_offZ_ge1b_ge4j_BothMblge25le250", pass_jgeneral && pass_2l_opc && !onZ && NGoodBJets_pt30 >= 1 && NGoodJets_pt30 >= 4 && pass_bothMbl},
            {"_2l_offZ_ge1b_BothMblge25le250_HTge200", pass_jgeneral && pass_2l_opc && !onZ && NGoodBJets_pt30 >= 1 && pass_bothMbl && HT_trigger_pt30 >= 200},
            {"_2l_offZ_ge1b_BothMblge25le250_HTge200_ge4j", pass_jgeneral && pass_2l_opc && !onZ && NGoodBJets_pt30 >= 1 && pass_bothMbl && HT_trigger_pt30 >= 200 && NGoodJets_pt30 >= 4}
       	};
                
	if (!inithisto) {
		InitHistos(cutmap);
		inithisto = true;
	}
 
	my_histos["EventCounter" ]->Fill( eventCounter );
        if( !passMadHT ) continue; //Make sure not to double count DY events

	for (const auto& cutVar: cutmap) 
        {
		if (cutVar.second) {
			my_histos["h_ntops_"+cutVar.first]->Fill( ntops, weight );
			my_histos["h_ht_"+cutVar.first]->Fill( HT_trigger_pt30, weight );		
			my_histos["h_njets_"+cutVar.first]->Fill( NGoodJets_pt30, weight );
			my_histos["h_nbjets_"+cutVar.first]->Fill( NGoodBJets_pt30, weight );
			my_2d_histos["h_njets_MVA_"+cutVar.first]->Fill( NGoodJets_pt30, deepESM_val, weight );
                        my_histos["h_LepLepDeltaR_"+cutVar.first]->Fill( twolepDeltaR, weight );
                        my_histos["h_Mbl_"+cutVar.first]->Fill( Mbl, weight );
                        my_histos["h_2b_nonbMass_"+cutVar.first]->Fill( NonBJetsSum.M(), weight);
                        my_histos["h_lepton_mass_"+cutVar.first]->Fill(lepton_sum.M(), weight);
                        my_histos["h_Mbl1_"+cutVar.first]->Fill(Mbl1, weight);
                        my_histos["h_Mbl2_"+cutVar.first]->Fill(Mbl2, weight);
                        my_histos["h_njetsPlus1_"+cutVar.first]->Fill( NGoodJets_pt30, weight );
                        for(int j = 0; j < NGoodJets_pt30; j++) {
                            my_histos["h_jet_pt_"+cutVar.first]->Fill(GoodJetsVec_pt30.at(j).Pt(), weight);
                            my_histos["h_jet_M_"+cutVar.first]->Fill(GoodJetsVec_pt30.at(j).M(), weight);
                        }   
                        for(int l = 0; l < GoodLeptons.size(); l++) {
                            my_histos["h_lepton_pt_"+cutVar.first]->Fill(GoodLeptons.at(l).second.Pt(), weight);
                           
                        }
		}
	}

        
        // Example Fill event selection efficiencies
        my_efficiencies["event_sel_weight"]->SetUseWeightedEvents();
        // my_efficiencies["event_sel_weight"]->FillWeighted(true,eventweight,0);
        // my_efficiencies["event_sel_weight"]->FillWeighted(true && pass_jgeneral,eventweight,1);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && pass_jgeneral && pass_2l_opc,eventweight,0);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 1,eventweight,1);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 1 && NGoodJets_pt30 >= 4,eventweight,2);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && pass_jgeneral && pass_2l_opc && NGoodBJets_pt30 >= 1 && NGoodJets_pt30 >= 4 && pass_bothMbl ,eventweight,3);
    }
   
}


void TwoLepAnalyzer::WriteHistos(TFile* outfile)
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
