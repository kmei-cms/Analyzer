#define Analyze1Lep_cxx
#include "Analyzer/Analyzer/include/Analyze1Lep.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

Analyze1Lep::Analyze1Lep() : initHistos(false)
{
}

void Analyze1Lep::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    int fB = 200;

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("fwm2_top6", std::make_shared<TH1D>("fwm2_top6","fwm2_top6", 50, 0, 1 ) );
    my_histos.emplace("fwm3_top6", std::make_shared<TH1D>("fwm3_top6","fwm3_top6", 50, 0, 1 ) );
    my_histos.emplace("fwm4_top6", std::make_shared<TH1D>("fwm4_top6","fwm4_top6", 50, 0, 1 ) );
    my_histos.emplace("fwm5_top6", std::make_shared<TH1D>("fwm5_top6","fwm5_top6", 50, 0, 1 ) );
    my_histos.emplace("fwm6_top6", std::make_shared<TH1D>("fwm6_top6","fwm6_top6", 50, 0, 1 ) );
    my_histos.emplace("fwm7_top6", std::make_shared<TH1D>("fwm7_top6","fwm7_top6", 50, 0, 1 ) );
    my_histos.emplace("fwm8_top6", std::make_shared<TH1D>("fwm8_top6","fwm8_top6", 50, 0, 1 ) );
    my_histos.emplace("fwm9_top6", std::make_shared<TH1D>("fwm9_top6","fwm9_top6", 50, 0, 1 ) );
    my_histos.emplace("fwm10_top6", std::make_shared<TH1D>("fwm10_top6","fwm10_top6", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev0_top6", std::make_shared<TH1D>("jmt_ev0_top6","jmt_ev0_top6", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev1_top6", std::make_shared<TH1D>("jmt_ev1_top6","jmt_ev1_top6", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev2_top6", std::make_shared<TH1D>("jmt_ev2_top6","jmt_ev2_top6", 50, 0, 1 ) );

    for(unsigned int i = 1; i <= 7 ; i++) //Bad hard code
    {
        my_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b",  std::make_shared<TH1D>(("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 150, 0, 1500 ));
        my_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH1D>(("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 100, -6, 6 ));
        my_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH1D>(("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 80, -4, 4 ));
        my_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b",   std::make_shared<TH1D>(("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 20, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), fB, 0.0, 1.0, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), fB, 0.0, 1.0, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), fB, 0.0, 1.0, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), fB, 0.0, 1.0, 40, 0, 200));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 40, 0, 200));
    }

    for(auto& mycut : cutMap)
    {
        my_histos.emplace("h_njets"+mycut.first, std::make_shared<TH1D>(("h_njets"+mycut.first).c_str(),("h_njets"+mycut.first).c_str(), 20, 0, 20));
        my_histos.emplace("h_ntops"+mycut.first, std::make_shared<TH1D>(("h_ntops"+mycut.first).c_str(),("h_ntops"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("h_nb"+mycut.first,    std::make_shared<TH1D>(("h_nb"+mycut.first).c_str(),("h_nb"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("h_fisher"+mycut.first,std::make_shared<TH1D>(("h_fisher"+mycut.first).c_str(),("h_fisher"+mycut.first).c_str(), fB, -0.5, 0.5));
        my_histos.emplace("h_deepESM"+mycut.first,std::make_shared<TH1D>(("h_deepESM"+mycut.first).c_str(),("h_deepESM"+mycut.first).c_str(), fB, 0.0, 1.0));
        my_histos.emplace("h_ht"+mycut.first,std::make_shared<TH1D>(("h_ht"+mycut.first).c_str(),("h_ht"+mycut.first).c_str(), 300, 0.0, 3000));
        my_histos.emplace("h_mbl"+mycut.first,std::make_shared<TH1D>(("h_mbl"+mycut.first).c_str(),("h_mbl"+mycut.first).c_str(), 300, 0.0, 300));
        my_histos.emplace("h_lPt"+mycut.first,std::make_shared<TH1D>(("h_lPt"+mycut.first).c_str(),("h_lPt"+mycut.first).c_str(), 200, 0.0, 2000));
        my_histos.emplace("h_lEta"+mycut.first,std::make_shared<TH1D>(("h_lEta"+mycut.first).c_str(),("h_lEta"+mycut.first).c_str(), 100, -6.0, 6.0));

        my_2d_histos.emplace("h_njets_bdt"+mycut.first, std::make_shared<TH2D>(("h_njets_bdt"+mycut.first).c_str(),("h_njets_bdt"+mycut.first).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_fisher"+mycut.first, std::make_shared<TH2D>(("h_njets_fisher"+mycut.first).c_str(),("h_njets_fisher"+mycut.first).c_str(), 15, 0, 15, fB, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_deepESM"+mycut.first, std::make_shared<TH2D>(("h_njets_deepESM"+mycut.first).c_str(),("h_njets_deepESM"+mycut.first).c_str(), 15, 0, 15, fB, 0.0, 1.0));
        my_2d_histos.emplace("h_njets_mbl"+mycut.first, std::make_shared<TH2D>(("h_njets_mbl"+mycut.first).c_str(),("h_njets_mbl"+mycut.first).c_str(), 15, 0, 15, 300, 0, 300));
        my_2d_histos.emplace("h_ht_deepESM"+mycut.first, std::make_shared<TH2D>(("h_ht_deepESM"+mycut.first).c_str(),("h_ht_deepESM"+mycut.first).c_str(), 300, 0, 3000, fB, 0.0, 1.0));

        my_histos.emplace("blind_njets"+mycut.first, std::make_shared<TH1D>(("blind_njets"+mycut.first).c_str(),("blind_njets"+mycut.first).c_str(), 20, 0, 20));
        my_histos.emplace("blind_ntops"+mycut.first, std::make_shared<TH1D>(("blind_ntops"+mycut.first).c_str(),("blind_ntops"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("blind_nb"+mycut.first,    std::make_shared<TH1D>(("blind_nb"+mycut.first).c_str(),("blind_nb"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("blind_bdt"+mycut.first,   std::make_shared<TH1D>(("blind_bdt"+mycut.first).c_str(),("blind_bdt"+mycut.first).c_str(), 40, -0.5, 0.5));
        my_histos.emplace("blind_fisher"+mycut.first,std::make_shared<TH1D>(("blind_fisher"+mycut.first).c_str(),("blind_fisher"+mycut.first).c_str(), fB, -0.5, 0.5));
        my_histos.emplace("blind_deepESM"+mycut.first,std::make_shared<TH1D>(("blind_deepESM"+mycut.first).c_str(),("blind_deepESM"+mycut.first).c_str(), fB, 0.0, 1.0));
        my_histos.emplace("blind_ht"+mycut.first,std::make_shared<TH1D>(("blind_ht"+mycut.first).c_str(),("blind_ht"+mycut.first).c_str(), 300, 0.0, 3000));
        my_histos.emplace("blind_mbl"+mycut.first,std::make_shared<TH1D>(("blind_mbl"+mycut.first).c_str(),("blind_mbl"+mycut.first).c_str(), 300, 0.0, 300));
        my_histos.emplace("blind_lPt"+mycut.first,std::make_shared<TH1D>(("blind_lPt"+mycut.first).c_str(),("blind_lPt"+mycut.first).c_str(), 200, 0.0, 2000));
        my_histos.emplace("blind_lEta"+mycut.first,std::make_shared<TH1D>(("blind_lEta"+mycut.first).c_str(),("blind_lEta"+mycut.first).c_str(), 100, -6.0, 6.0));
    
        my_2d_histos.emplace("blind_njets_bdt"+mycut.first, std::make_shared<TH2D>(("blind_njets_bdt"+mycut.first).c_str(),("blind_njets_bdt"+mycut.first).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("blind_njets_fisher"+mycut.first, std::make_shared<TH2D>(("blind_njets_fisher"+mycut.first).c_str(),("blind_njets_fisher"+mycut.first).c_str(), 15, 0, 15, fB, -0.5, 0.5));
        my_2d_histos.emplace("blind_njets_deepESM"+mycut.first, std::make_shared<TH2D>(("blind_njets_deepESM"+mycut.first).c_str(),("blind_njets_deepESM"+mycut.first).c_str(), 15, 0, 15, fB, 0.0, 1.0));
        my_2d_histos.emplace("blind_njets_mbl"+mycut.first, std::make_shared<TH2D>(("blind_njets_mbl"+mycut.first).c_str(),("blind_njets_mbl"+mycut.first).c_str(), 15, 0, 15, 300, 0, 300));
        my_2d_histos.emplace("blind_ht_deepESM"+mycut.first, std::make_shared<TH2D>(("blind_ht_deepESM"+mycut.first).c_str(),("blind_ht_deepESM"+mycut.first).c_str(), 300, 0, 3000, fB, 0.0, 1.0));
    }
}

void Analyze1Lep::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& MET                  = tr.getVar<double>("MET");
        const auto& ntops                = tr.getVar<int>("ntops");
        const auto& ntops_3jet           = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet           = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet           = tr.getVar<int>("ntops_1jet");
        const auto& runtype              = tr.getVar<std::string>("runtype");     
        const auto& filetag              = tr.getVar<std::string>("filetag");
        const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodJets_pt45       = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets_pt30      = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodBJets_pt45      = tr.getVar<int>("NGoodBJets_pt45");
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons");
        const auto& GoodLeptons          = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& fisher_val           = tr.getVar<double>("fisher_val");
        const auto& HT_trigger           = tr.getVar<double>("HT_trigger");
        const auto& JetID                = tr.getVar<bool>("JetID");
        const auto& passTrigger          = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC        = tr.getVar<bool>("passTriggerMC");
        const auto& passMadHT            = tr.getVar<bool>("passMadHT");
              auto  passBaseline1l_Good  = tr.getVar<bool>("passBaseline1l_Good");
        const auto& Mbl                  = tr.getVar<double>("Mbl");
                    passBaseline1l_Good  = passBaseline1l_Good && Mbl>30 && Mbl<180;
        const auto& passBlind            = tr.getVar<bool>("passBlindLep_Good");            
        const auto& deepESM_val          = tr.getVar<double>("deepESM_val");
        const auto& deepESM_bin1         = tr.getVar<bool>("deepESM_bin1");
        const auto& deepESM_bin2         = tr.getVar<bool>("deepESM_bin2");
        const auto& deepESM_bin3         = tr.getVar<bool>("deepESM_bin3");
        const auto& deepESM_bin4         = tr.getVar<bool>("deepESM_bin4");
        const auto& fwm2_top6            = tr.getVar<double>("fwm2_top6");
        const auto& fwm3_top6            = tr.getVar<double>("fwm3_top6");
        const auto& fwm4_top6            = tr.getVar<double>("fwm4_top6");
        const auto& fwm5_top6            = tr.getVar<double>("fwm5_top6");
        const auto& fwm6_top6            = tr.getVar<double>("fwm6_top6");
        const auto& fwm7_top6            = tr.getVar<double>("fwm7_top6");
        const auto& fwm8_top6            = tr.getVar<double>("fwm8_top6");
        const auto& fwm9_top6            = tr.getVar<double>("fwm9_top6");
        const auto& fwm10_top6           = tr.getVar<double>("fwm10_top6");
        const auto& jmt_ev0_top6         = tr.getVar<double>("jmt_ev0_top6");
        const auto& jmt_ev1_top6         = tr.getVar<double>("jmt_ev1_top6");
        const auto& jmt_ev2_top6         = tr.getVar<double>("jmt_ev2_top6");
        const auto& Jets_cm_top6         = tr.getVec<TLorentzVector>("Jets_cm_top6");

        // ------------------------
        // -- Print event number
        // ------------------------       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if(tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() );

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight = 1.0;
        double leptonweight = 1.0;
        double PileupWeight = 1.0;
        double bTagWeight = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            double eventweight = 1.0;        
            const auto& Weight  = tr.getVar<double>("Weight");
            double lumi = 35900; // Lumi for 2016
            eventweight = lumi*Weight;
            
            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonweight = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }

            PileupWeight = tr.getVar<double>("_PUweightFactor");
            bTagWeight   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");

            weight *= eventweight*leptonweight*bTagWeight;
        }

        // -------------------------------
        // -- Define cuts
        // -------------------------------
        bool pass_0l              = NGoodLeptons == 0;
        bool pass_1l              = NGoodLeptons == 1;
        bool pass_njet_pt30       = NGoodJets_pt30 >= 7;
        bool pass_njet_pt30_1btag = NGoodBJets_pt30 >= 1;
        bool pass_njet_pt30_2btag = NGoodBJets_pt30 >= 2;
        bool pass_1e              = false;
        bool pass_1m              = false;
        if(pass_1l)
        { 
            if(GoodLeptons[0].first == "e") pass_1e = true;
            if(GoodLeptons[0].first == "m") pass_1m = true;
        }
        bool pass_5to6njet_pt30 = (NGoodJets_pt30 == 5 || NGoodJets_pt30 == 6);

        // ---------------------------
        // --    1 Top selection
        // ---------------------------
        //exactly selections
        bool pass_1t    = ntops==1;
        bool pass_1t_d1 = pass_1t && deepESM_bin1, pass_1t_d2 = pass_1t && deepESM_bin2, pass_1t_d3 = pass_1t && deepESM_bin3, pass_1t_d4 = pass_1t && deepESM_bin4;

        bool pass_1t1 = pass_1t && ntops_1jet==1;
        bool pass_1t2 = pass_1t && ntops_2jet==1;
        bool pass_1t3 = pass_1t && ntops_3jet==1;
        bool pass_1t2or3 = pass_1t2 or pass_1t3;

        bool pass_1t1_d1 = pass_1t1 && deepESM_bin1, pass_1t1_d2 = pass_1t1 && deepESM_bin2, pass_1t1_d3 = pass_1t1 && deepESM_bin3, pass_1t1_d4 = pass_1t1 && deepESM_bin4;
        bool pass_1t2_d1 = pass_1t2 && deepESM_bin1, pass_1t2_d2 = pass_1t2 && deepESM_bin2, pass_1t2_d3 = pass_1t2 && deepESM_bin3, pass_1t2_d4 = pass_1t2 && deepESM_bin4;
        bool pass_1t3_d1 = pass_1t3 && deepESM_bin1, pass_1t3_d2 = pass_1t3 && deepESM_bin2, pass_1t3_d3 = pass_1t3 && deepESM_bin3, pass_1t3_d4 = pass_1t3 && deepESM_bin4;
        bool pass_1t2or3_d1 = pass_1t2or3 && deepESM_bin1, pass_1t2or3_d2 = pass_1t2or3 && deepESM_bin2, pass_1t2or3_d3 = pass_1t2or3 && deepESM_bin3, pass_1t2or3_d4 = pass_1t2or3 && deepESM_bin4;

        //at least selections
        bool pass_ge1t    = ntops>=1;
        bool pass_ge1t_d1 = pass_ge1t && deepESM_bin1, pass_ge1t_d2 = pass_ge1t && deepESM_bin2, pass_ge1t_d3 = pass_ge1t && deepESM_bin3, pass_ge1t_d4 = pass_ge1t && deepESM_bin4;

        bool pass_ge1t1 = pass_ge1t && ntops_1jet>=1;
        bool pass_ge1t2 = pass_ge1t && ntops_2jet>=1;
        bool pass_ge1t3 = pass_ge1t && ntops_3jet>=1;

        bool pass_ge1t1_d1 = pass_ge1t1 && deepESM_bin1, pass_ge1t1_d2 = pass_ge1t1 && deepESM_bin2, pass_ge1t1_d3 = pass_ge1t1 && deepESM_bin3, pass_ge1t1_d4 = pass_ge1t1 && deepESM_bin4;
        bool pass_ge1t2_d1 = pass_ge1t2 && deepESM_bin1, pass_ge1t2_d2 = pass_ge1t2 && deepESM_bin2, pass_ge1t2_d3 = pass_ge1t2 && deepESM_bin3, pass_ge1t2_d4 = pass_ge1t2 && deepESM_bin4;
        bool pass_ge1t3_d1 = pass_ge1t3 && deepESM_bin1, pass_ge1t3_d2 = pass_ge1t3 && deepESM_bin2, pass_ge1t3_d3 = pass_ge1t3 && deepESM_bin3, pass_ge1t3_d4 = pass_ge1t3 && deepESM_bin4;
        // ---------------------------
        // --    2 Top selection
        // ---------------------------

        // exactly selections
        bool pass_2t    = ntops==2;
        bool pass_2t_d1 = pass_2t && deepESM_bin1, pass_2t_d2 = pass_2t && deepESM_bin2, pass_2t_d3 = pass_2t && deepESM_bin3, pass_2t_d4 = pass_2t && deepESM_bin4;

        bool pass_2t11 = pass_2t && ntops_1jet==2;
        bool pass_2t12 = pass_2t && ntops_1jet==1 && ntops_2jet==1;
        bool pass_2t13 = pass_2t && ntops_1jet==1 && ntops_3jet==1;
        bool pass_2t22 = pass_2t && ntops_2jet==2;
        bool pass_2t23 = pass_2t && ntops_2jet==1 && ntops_3jet==1;
        bool pass_2t33 = pass_2t && ntops_3jet==2;

        bool pass_2t11_d1 = pass_2t11 && deepESM_bin1, pass_2t11_d2 = pass_2t11 && deepESM_bin2, pass_2t11_d3 = pass_2t11 && deepESM_bin3, pass_2t11_d4 = pass_2t11 && deepESM_bin4;
        bool pass_2t12_d1 = pass_2t12 && deepESM_bin1, pass_2t12_d2 = pass_2t12 && deepESM_bin2, pass_2t12_d3 = pass_2t12 && deepESM_bin3, pass_2t12_d4 = pass_2t12 && deepESM_bin4;
        bool pass_2t13_d1 = pass_2t13 && deepESM_bin1, pass_2t13_d2 = pass_2t13 && deepESM_bin2, pass_2t13_d3 = pass_2t13 && deepESM_bin3, pass_2t13_d4 = pass_2t13 && deepESM_bin4;
        bool pass_2t22_d1 = pass_2t22 && deepESM_bin1, pass_2t22_d2 = pass_2t22 && deepESM_bin2, pass_2t22_d3 = pass_2t22 && deepESM_bin3, pass_2t22_d4 = pass_2t22 && deepESM_bin4;
        bool pass_2t23_d1 = pass_2t23 && deepESM_bin1, pass_2t23_d2 = pass_2t23 && deepESM_bin2, pass_2t23_d3 = pass_2t23 && deepESM_bin3, pass_2t23_d4 = pass_2t23 && deepESM_bin4;
        bool pass_2t33_d1 = pass_2t33 && deepESM_bin1, pass_2t33_d2 = pass_2t33 && deepESM_bin2, pass_2t33_d3 = pass_2t33 && deepESM_bin3, pass_2t33_d4 = pass_2t33 && deepESM_bin4;

        // at least selections
        bool pass_ge2t    = ntops>=2;
        bool pass_ge2t_d1 = pass_ge2t && deepESM_bin1, pass_ge2t_d2 = pass_ge2t && deepESM_bin2, pass_ge2t_d3 = pass_ge2t && deepESM_bin3, pass_ge2t_d4 = pass_ge2t && deepESM_bin4;

        bool pass_ge2t11 = pass_ge2t && ntops_1jet>=2;
        bool pass_ge2t12 = pass_ge2t && ntops_1jet>=1 && ntops_2jet>=1;
        bool pass_ge2t13 = pass_ge2t && ntops_1jet>=1 && ntops_3jet>=1;
        bool pass_ge2t22 = pass_ge2t && ntops_2jet>=2;
        bool pass_ge2t23 = pass_ge2t && ntops_2jet>=1 && ntops_3jet>=1;
        bool pass_ge2t33 = pass_ge2t && ntops_3jet>=2;

        //Merge or no Merge selections
        bool pass_Merge = ntops_1jet >= 1, pass_noMerge = ntops_1jet == 0;
        bool pass_ge2t11or12or13 = pass_Merge   && (pass_ge2t11 or pass_ge2t12 or pass_ge2t13);
        bool pass_ge2t22or23or33 = pass_noMerge && (pass_ge2t22 or pass_ge2t23 or pass_ge2t33);

        bool pass_ge2t11_d1 = pass_ge2t11 && deepESM_bin1, pass_ge2t11_d2 = pass_ge2t11 && deepESM_bin2, pass_ge2t11_d3 = pass_ge2t11 && deepESM_bin3, pass_ge2t11_d4 = pass_ge2t11 && deepESM_bin4;
        bool pass_ge2t12_d1 = pass_ge2t12 && deepESM_bin1, pass_ge2t12_d2 = pass_ge2t12 && deepESM_bin2, pass_ge2t12_d3 = pass_ge2t12 && deepESM_bin3, pass_ge2t12_d4 = pass_ge2t12 && deepESM_bin4;
        bool pass_ge2t13_d1 = pass_ge2t13 && deepESM_bin1, pass_ge2t13_d2 = pass_ge2t13 && deepESM_bin2, pass_ge2t13_d3 = pass_ge2t13 && deepESM_bin3, pass_ge2t13_d4 = pass_ge2t13 && deepESM_bin4;
        bool pass_ge2t22_d1 = pass_ge2t22 && deepESM_bin1, pass_ge2t22_d2 = pass_ge2t22 && deepESM_bin2, pass_ge2t22_d3 = pass_ge2t22 && deepESM_bin3, pass_ge2t22_d4 = pass_ge2t22 && deepESM_bin4;
        bool pass_ge2t23_d1 = pass_ge2t23 && deepESM_bin1, pass_ge2t23_d2 = pass_ge2t23 && deepESM_bin2, pass_ge2t23_d3 = pass_ge2t23 && deepESM_bin3, pass_ge2t23_d4 = pass_ge2t23 && deepESM_bin4;
        bool pass_ge2t33_d1 = pass_ge2t33 && deepESM_bin1, pass_ge2t33_d2 = pass_ge2t33 && deepESM_bin2, pass_ge2t33_d3 = pass_ge2t33 && deepESM_bin3, pass_ge2t33_d4 = pass_ge2t33 && deepESM_bin4;

        bool pass_ge2t11or12or13_d1 = pass_ge2t11or12or13 && deepESM_bin1, pass_ge2t11or12or13_d2 = pass_ge2t11or12or13 && deepESM_bin2, pass_ge2t11or12or13_d3 = pass_ge2t11or12or13 && deepESM_bin3, pass_ge2t11or12or13_d4 = pass_ge2t11or12or13 && deepESM_bin4;
        bool pass_ge2t22or23or33_d1 = pass_ge2t22or23or33 && deepESM_bin1, pass_ge2t22or23or33_d2 = pass_ge2t22or23or33 && deepESM_bin2, pass_ge2t22or23or33_d3 = pass_ge2t22or23or33 && deepESM_bin3, pass_ge2t22or23or33_d4 = pass_ge2t22or23or33 && deepESM_bin4;

        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map_1l 
        {
            {""                                , true                                                                     },
            {"_1l"                             , pass_1l                                                                  },
            {"_1l_ge7j"                        , pass_1l && pass_njet_pt30 && JetID                                       },
            {"_1l_ge2b"                        , pass_1l && pass_njet_pt30_2btag && JetID                                 },
            {"_1l_1t"                          , pass_1l && pass_1t                                                       },
            {"_1l_2t"                          , pass_1l && pass_2t                                                       },
            {"_1l_ge1t"                        , pass_1l && pass_ge1t                                                     },
            {"_1l_ge2t"                        , pass_1l && pass_ge2t                                                     },
            {"_1l_5to6j_ge1b"                  , pass_1l && pass_5to6njet_pt30 && pass_njet_pt30_1btag && JetID           },
            {"_1l_ge7j_ge2b"                   , pass_1l && pass_njet_pt30 && pass_njet_pt30_2btag && JetID               },
            {"_1l_ge7j_ge1b_noMbl"             , pass_1l && pass_njet_pt30 && pass_njet_pt30_1btag && JetID               },
            {"_1l_ge7j_ge1b"                   , passBaseline1l_Good                                                      },                         
            {"_1e_ge7j_ge1b"                   , passBaseline1l_Good && pass_1e                                           },                         
            {"_1m_ge7j_ge1b"                   , passBaseline1l_Good && pass_1m                                           },                         
            {"_1l_ge7j_ge1b_d1"                , passBaseline1l_Good && deepESM_bin1                                      },                         
            {"_1l_ge7j_ge1b_d2"                , passBaseline1l_Good && deepESM_bin2                                      },                         
            {"_1l_ge7j_ge1b_d3"                , passBaseline1l_Good && deepESM_bin3                                      },                         
            {"_1l_ge7j_ge1b_d4"                , passBaseline1l_Good && deepESM_bin4                                      },                         
            {"_1l_6j_ge1b"                     , passBaseline1l_Good && NGoodJets_pt30 == 6                               },
            {"_1l_7j_ge1b"                     , passBaseline1l_Good && NGoodJets_pt30 == 7                               },
            {"_1l_8j_ge1b"                     , passBaseline1l_Good && NGoodJets_pt30 == 8                               },
            {"_1l_9j_ge1b"                     , passBaseline1l_Good && NGoodJets_pt30 == 9                               },
            {"_1l_10j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 10                              },
            {"_1l_11j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 11                              },
            {"_1l_12j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 12                              },
            {"_1l_13j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 13                              },
            {"_1l_14j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 14                              },
            {"_1l_15j_ge1b"                    , passBaseline1l_Good && NGoodJets_pt30 == 15                              },
                                                 
            {"_1l_ge7j_ge1b_1t"                , passBaseline1l_Good && pass_1t                                           },
            {"_1l_ge7j_ge1b_1t_d1"             , passBaseline1l_Good && pass_1t_d1                                        },
            {"_1l_ge7j_ge1b_1t_d2"             , passBaseline1l_Good && pass_1t_d2                                        }, 
            {"_1l_ge7j_ge1b_1t_d3"             , passBaseline1l_Good && pass_1t_d3                                        }, 
            {"_1l_ge7j_ge1b_1t_d4"             , passBaseline1l_Good && pass_1t_d4                                        },
            {"_1l_ge7j_ge1b_ge1t"              , passBaseline1l_Good && pass_ge1t                                         },
            {"_1l_ge7j_ge1b_ge1t_d1"           , passBaseline1l_Good && pass_ge1t_d1                                      },
            {"_1l_ge7j_ge1b_ge1t_d2"           , passBaseline1l_Good && pass_ge1t_d2                                      }, 
            {"_1l_ge7j_ge1b_ge1t_d3"           , passBaseline1l_Good && pass_ge1t_d3                                      }, 
            {"_1l_ge7j_ge1b_ge1t_d4"           , passBaseline1l_Good && pass_ge1t_d4                                      },
                                                 
            {"_1l_ge7j_ge1b_1t1"               , passBaseline1l_Good && pass_1t1                                          },
            {"_1l_ge7j_ge1b_1t2"               , passBaseline1l_Good && pass_1t2                                          },
            {"_1l_ge7j_ge1b_1t3"               , passBaseline1l_Good && pass_1t3                                          },
            {"_1l_ge7j_ge1b_1t2or3"            , passBaseline1l_Good && pass_1t2or3                                       },
            {"_1l_ge7j_ge1b_1t1_d1"            , passBaseline1l_Good && pass_1t1_d1                                       },
            {"_1l_ge7j_ge1b_1t1_d2"            , passBaseline1l_Good && pass_1t1_d2                                       },
            {"_1l_ge7j_ge1b_1t1_d3"            , passBaseline1l_Good && pass_1t1_d3                                       },
            {"_1l_ge7j_ge1b_1t1_d4"            , passBaseline1l_Good && pass_1t1_d4                                       },
            {"_1l_ge7j_ge1b_1t2_d1"            , passBaseline1l_Good && pass_1t2_d1                                       },
            {"_1l_ge7j_ge1b_1t2_d2"            , passBaseline1l_Good && pass_1t2_d2                                       },
            {"_1l_ge7j_ge1b_1t2_d3"            , passBaseline1l_Good && pass_1t2_d3                                       },
            {"_1l_ge7j_ge1b_1t2_d4"            , passBaseline1l_Good && pass_1t2_d4                                       },
            {"_1l_ge7j_ge1b_1t3_d1"            , passBaseline1l_Good && pass_1t3_d1                                       },
            {"_1l_ge7j_ge1b_1t3_d2"            , passBaseline1l_Good && pass_1t3_d2                                       },
            {"_1l_ge7j_ge1b_1t3_d3"            , passBaseline1l_Good && pass_1t3_d3                                       },
            {"_1l_ge7j_ge1b_1t3_d4"            , passBaseline1l_Good && pass_1t3_d4                                       },
            {"_1l_ge7j_ge1b_1t2or3_d1"         , passBaseline1l_Good && pass_1t2or3_d1                                    },
            {"_1l_ge7j_ge1b_1t2or3_d2"         , passBaseline1l_Good && pass_1t2or3_d2                                    },
            {"_1l_ge7j_ge1b_1t2or3_d3"         , passBaseline1l_Good && pass_1t2or3_d3                                    },
            {"_1l_ge7j_ge1b_1t2or3_d4"         , passBaseline1l_Good && pass_1t2or3_d4                                    },
            {"_1l_ge7j_ge1b_ge1t1"             , passBaseline1l_Good && pass_ge1t1                                        },
            {"_1l_ge7j_ge1b_ge1t2"             , passBaseline1l_Good && pass_ge1t2                                        },
            {"_1l_ge7j_ge1b_ge1t3"             , passBaseline1l_Good && pass_ge1t3                                        },
            {"_1l_ge7j_ge1b_ge1t1_d1"          , passBaseline1l_Good && pass_ge1t1_d1                                     },
            {"_1l_ge7j_ge1b_ge1t1_d2"          , passBaseline1l_Good && pass_ge1t1_d2                                     },
            {"_1l_ge7j_ge1b_ge1t1_d3"          , passBaseline1l_Good && pass_ge1t1_d3                                     },
            {"_1l_ge7j_ge1b_ge1t1_d4"          , passBaseline1l_Good && pass_ge1t1_d4                                     },
            {"_1l_ge7j_ge1b_ge1t2_d1"          , passBaseline1l_Good && pass_ge1t2_d1                                     },
            {"_1l_ge7j_ge1b_ge1t2_d2"          , passBaseline1l_Good && pass_ge1t2_d2                                     },
            {"_1l_ge7j_ge1b_ge1t2_d3"          , passBaseline1l_Good && pass_ge1t2_d3                                     },
            {"_1l_ge7j_ge1b_ge1t2_d4"          , passBaseline1l_Good && pass_ge1t2_d4                                     },
            {"_1l_ge7j_ge1b_ge1t3_d1"          , passBaseline1l_Good && pass_ge1t3_d1                                     },
            {"_1l_ge7j_ge1b_ge1t3_d2"          , passBaseline1l_Good && pass_ge1t3_d2                                     },
            {"_1l_ge7j_ge1b_ge1t3_d3"          , passBaseline1l_Good && pass_ge1t3_d3                                     },
            {"_1l_ge7j_ge1b_ge1t3_d4"          , passBaseline1l_Good && pass_ge1t3_d4                                     },
                                                 
            {"_1l_ge7j_ge1b_2t"                , passBaseline1l_Good && pass_2t                                           },
            {"_1l_ge7j_ge1b_2t_d1"             , passBaseline1l_Good && pass_2t_d1                                        },
            {"_1l_ge7j_ge1b_2t_d2"             , passBaseline1l_Good && pass_2t_d2                                        }, 
            {"_1l_ge7j_ge1b_2t_d3"             , passBaseline1l_Good && pass_2t_d3                                        }, 
            {"_1l_ge7j_ge1b_2t_d4"             , passBaseline1l_Good && pass_2t_d4                                        },
            {"_1l_ge7j_ge1b_ge2t"              , passBaseline1l_Good && pass_ge2t                                         },
            {"_1l_ge7j_ge1b_ge2t_d1"           , passBaseline1l_Good && pass_ge2t_d1                                      },
            {"_1l_ge7j_ge1b_ge2t_d2"           , passBaseline1l_Good && pass_ge2t_d2                                      }, 
            {"_1l_ge7j_ge1b_ge2t_d3"           , passBaseline1l_Good && pass_ge2t_d3                                      }, 
            {"_1l_ge7j_ge1b_ge2t_d4"           , passBaseline1l_Good && pass_ge2t_d4                                      },
                                                 
            //{"_1l_ge7j_ge1b_2t11"              , passBaseline1l_Good && pass_2t11                                         },
            //{"_1l_ge7j_ge1b_2t12"              , passBaseline1l_Good && pass_2t12                                         },
            //{"_1l_ge7j_ge1b_2t13"              , passBaseline1l_Good && pass_2t13                                         },
            //{"_1l_ge7j_ge1b_2t22"              , passBaseline1l_Good && pass_2t22                                         },
            //{"_1l_ge7j_ge1b_2t23"              , passBaseline1l_Good && pass_2t23                                         },
            //{"_1l_ge7j_ge1b_2t33"              , passBaseline1l_Good && pass_2t33                                         },
            //{"_1l_ge7j_ge1b_2t11_d1"           , passBaseline1l_Good && pass_2t11_d1                                      },
            //{"_1l_ge7j_ge1b_2t11_d2"           , passBaseline1l_Good && pass_2t11_d2                                      },
            //{"_1l_ge7j_ge1b_2t11_d3"           , passBaseline1l_Good && pass_2t11_d3                                      },
            //{"_1l_ge7j_ge1b_2t11_d4"           , passBaseline1l_Good && pass_2t11_d4                                      },
            //{"_1l_ge7j_ge1b_2t12_d1"           , passBaseline1l_Good && pass_2t12_d1                                      },
            //{"_1l_ge7j_ge1b_2t12_d2"           , passBaseline1l_Good && pass_2t12_d2                                      },
            //{"_1l_ge7j_ge1b_2t12_d3"           , passBaseline1l_Good && pass_2t12_d3                                      },
            //{"_1l_ge7j_ge1b_2t12_d4"           , passBaseline1l_Good && pass_2t12_d4                                      },
            //{"_1l_ge7j_ge1b_2t13_d1"           , passBaseline1l_Good && pass_2t13_d1                                      },
            //{"_1l_ge7j_ge1b_2t13_d2"           , passBaseline1l_Good && pass_2t13_d2                                      },
            //{"_1l_ge7j_ge1b_2t13_d3"           , passBaseline1l_Good && pass_2t13_d3                                      },
            //{"_1l_ge7j_ge1b_2t13_d4"           , passBaseline1l_Good && pass_2t13_d4                                      },
            //{"_1l_ge7j_ge1b_2t22_d1"           , passBaseline1l_Good && pass_2t22_d1                                      },
            //{"_1l_ge7j_ge1b_2t22_d2"           , passBaseline1l_Good && pass_2t22_d2                                      },
            //{"_1l_ge7j_ge1b_2t22_d3"           , passBaseline1l_Good && pass_2t22_d3                                      },
            //{"_1l_ge7j_ge1b_2t22_d4"           , passBaseline1l_Good && pass_2t22_d4                                      },
            //{"_1l_ge7j_ge1b_2t23_d1"           , passBaseline1l_Good && pass_2t23_d1                                      },
            //{"_1l_ge7j_ge1b_2t23_d2"           , passBaseline1l_Good && pass_2t23_d2                                      },
            //{"_1l_ge7j_ge1b_2t23_d3"           , passBaseline1l_Good && pass_2t23_d3                                      },
            //{"_1l_ge7j_ge1b_2t23_d4"           , passBaseline1l_Good && pass_2t23_d4                                      },
            //{"_1l_ge7j_ge1b_2t33_d1"           , passBaseline1l_Good && pass_2t33_d1                                      },
            //{"_1l_ge7j_ge1b_2t33_d2"           , passBaseline1l_Good && pass_2t33_d2                                      },
            //{"_1l_ge7j_ge1b_2t33_d3"           , passBaseline1l_Good && pass_2t33_d3                                      },
            //{"_1l_ge7j_ge1b_2t33_d4"           , passBaseline1l_Good && pass_2t33_d4                                      },
            //{"_1l_ge7j_ge1b_ge2t11"            , passBaseline1l_Good && pass_ge2t11                                       },
            //{"_1l_ge7j_ge1b_ge2t12"            , passBaseline1l_Good && pass_ge2t12                                       },
            //{"_1l_ge7j_ge1b_ge2t13"            , passBaseline1l_Good && pass_ge2t13                                       },
            //{"_1l_ge7j_ge1b_ge2t11or12or13"    , passBaseline1l_Good && pass_ge2t11or12or13                               },
            //{"_1l_ge7j_ge1b_ge2t22"            , passBaseline1l_Good && pass_ge2t22                                       },
            //{"_1l_ge7j_ge1b_ge2t23"            , passBaseline1l_Good && pass_ge2t23                                       },
            //{"_1l_ge7j_ge1b_ge2t33"            , passBaseline1l_Good && pass_ge2t33                                       },
            //{"_1l_ge7j_ge1b_ge2t22or23or33"    , passBaseline1l_Good && pass_ge2t22or23or33                               },
            //{"_1l_ge7j_ge1b_ge2t11_d1"         , passBaseline1l_Good && pass_ge2t11_d1                                    },
            //{"_1l_ge7j_ge1b_ge2t11_d2"         , passBaseline1l_Good && pass_ge2t11_d2                                    },
            //{"_1l_ge7j_ge1b_ge2t11_d3"         , passBaseline1l_Good && pass_ge2t11_d3                                    },
            //{"_1l_ge7j_ge1b_ge2t11_d4"         , passBaseline1l_Good && pass_ge2t11_d4                                    },
            //{"_1l_ge7j_ge1b_ge2t12_d1"         , passBaseline1l_Good && pass_ge2t12_d1                                    },
            //{"_1l_ge7j_ge1b_ge2t12_d2"         , passBaseline1l_Good && pass_ge2t12_d2                                    },
            //{"_1l_ge7j_ge1b_ge2t12_d3"         , passBaseline1l_Good && pass_ge2t12_d3                                    },
            //{"_1l_ge7j_ge1b_ge2t12_d4"         , passBaseline1l_Good && pass_ge2t12_d4                                    },
            //{"_1l_ge7j_ge1b_ge2t13_d1"         , passBaseline1l_Good && pass_ge2t13_d1                                    },
            //{"_1l_ge7j_ge1b_ge2t13_d2"         , passBaseline1l_Good && pass_ge2t13_d2                                    },
            //{"_1l_ge7j_ge1b_ge2t13_d3"         , passBaseline1l_Good && pass_ge2t13_d3                                    },
            //{"_1l_ge7j_ge1b_ge2t13_d4"         , passBaseline1l_Good && pass_ge2t13_d4                                    },
            //{"_1l_ge7j_ge1b_ge2t11or12or13_d1" , passBaseline1l_Good && pass_ge2t11or12or13_d1                            },
            //{"_1l_ge7j_ge1b_ge2t11or12or13_d2" , passBaseline1l_Good && pass_ge2t11or12or13_d2                            },
            //{"_1l_ge7j_ge1b_ge2t11or12or13_d3" , passBaseline1l_Good && pass_ge2t11or12or13_d3                            },
            //{"_1l_ge7j_ge1b_ge2t11or12or13_d4" , passBaseline1l_Good && pass_ge2t11or12or13_d4                            },
            //{"_1l_ge7j_ge1b_ge2t22_d1"         , passBaseline1l_Good && pass_ge2t22_d1                                    },
            //{"_1l_ge7j_ge1b_ge2t22_d2"         , passBaseline1l_Good && pass_ge2t22_d2                                    },
            //{"_1l_ge7j_ge1b_ge2t22_d3"         , passBaseline1l_Good && pass_ge2t22_d3                                    },
            //{"_1l_ge7j_ge1b_ge2t22_d4"         , passBaseline1l_Good && pass_ge2t22_d4                                    },
            //{"_1l_ge7j_ge1b_ge2t23_d1"         , passBaseline1l_Good && pass_ge2t23_d1                                    },
            //{"_1l_ge7j_ge1b_ge2t23_d2"         , passBaseline1l_Good && pass_ge2t23_d2                                    },
            //{"_1l_ge7j_ge1b_ge2t23_d3"         , passBaseline1l_Good && pass_ge2t23_d3                                    },
            //{"_1l_ge7j_ge1b_ge2t23_d4"         , passBaseline1l_Good && pass_ge2t23_d4                                    },
            //{"_1l_ge7j_ge1b_ge2t33_d1"         , passBaseline1l_Good && pass_ge2t33_d1                                    },
            //{"_1l_ge7j_ge1b_ge2t33_d2"         , passBaseline1l_Good && pass_ge2t33_d2                                    },
            //{"_1l_ge7j_ge1b_ge2t33_d3"         , passBaseline1l_Good && pass_ge2t33_d3                                    },
            //{"_1l_ge7j_ge1b_ge2t33_d4"         , passBaseline1l_Good && pass_ge2t33_d4                                    },
            //{"_1l_ge7j_ge1b_ge2t22or23or33_d1" , passBaseline1l_Good && pass_ge2t22or23or33_d1                            },
            //{"_1l_ge7j_ge1b_ge2t22or23or33_d2" , passBaseline1l_Good && pass_ge2t22or23or33_d2                            },
            //{"_1l_ge7j_ge1b_ge2t22or23or33_d3" , passBaseline1l_Good && pass_ge2t22or23or33_d3                            },
            //{"_1l_ge7j_ge1b_ge2t22or23or33_d4" , passBaseline1l_Good && pass_ge2t22or23or33_d4                            },
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map_1l);
            initHistos = true;
        }

        // Global cuts
        if ( !(passTriggerMC && passTrigger && passMadHT && passBlind) ) continue;

        for(auto& kv : cut_map_1l)
        {
            if(kv.second)
            {
                my_histos["h_njets"   +kv.first]->Fill(NGoodJets_pt30, weight);
                my_histos["h_ntops"   +kv.first]->Fill(ntops, weight);
                my_histos["h_nb"      +kv.first]->Fill(NGoodBJets_pt30, weight);
                my_histos["h_fisher"  +kv.first]->Fill(fisher_val, weight);
                my_histos["h_deepESM" +kv.first]->Fill(deepESM_val, weight);
                my_histos["h_ht"      +kv.first]->Fill(HT_trigger, weight);
                my_histos["h_mbl"     +kv.first]->Fill(Mbl, weight);
                for(const auto l : GoodLeptons)
                {
                    my_histos["h_lPt"+kv.first]->Fill(l.second.Pt(), weight);
                    my_histos["h_lEta"+kv.first]->Fill(l.second.Eta(), weight);
                }

                my_2d_histos["h_njets_fisher"+kv.first]->Fill(NGoodJets_pt30, fisher_val, weight);
                my_2d_histos["h_njets_deepESM"+kv.first]->Fill(NGoodJets_pt30, deepESM_val, weight);
                my_2d_histos["h_njets_mbl"+kv.first]->Fill(NGoodJets_pt30, Mbl, weight);
                my_2d_histos["h_ht_deepESM"+kv.first]->Fill(HT_trigger, deepESM_val, weight);

                if ( NGoodJets_pt30 <= 7 )
                {
                    my_histos["blind_njets"   +kv.first]->Fill(NGoodJets_pt30, weight);
                    my_histos["blind_ntops"   +kv.first]->Fill(ntops, weight);
                    my_histos["blind_nb"      +kv.first]->Fill(NGoodBJets_pt30, weight);
                    my_histos["blind_fisher"  +kv.first]->Fill(fisher_val, weight);
                    my_histos["blind_deepESM" +kv.first]->Fill(deepESM_val, weight);
                    my_histos["blind_ht"      +kv.first]->Fill(HT_trigger, weight);
                    my_histos["blind_mbl"     +kv.first]->Fill(Mbl, weight);
                    for(const auto l : GoodLeptons)
                    {
                        my_histos["blind_lPt"+kv.first]->Fill(l.second.Pt(), weight);
                        my_histos["blind_lEta"+kv.first]->Fill(l.second.Eta(), weight);
                    }
                
                    my_2d_histos["blind_njets_fisher"+kv.first]->Fill(NGoodJets_pt30, fisher_val, weight);
                    my_2d_histos["blind_njets_deepESM"+kv.first]->Fill(NGoodJets_pt30, deepESM_val, weight);
                    my_2d_histos["blind_njets_mbl"+kv.first]->Fill(NGoodJets_pt30, Mbl, weight);
                    my_2d_histos["blind_ht_deepESM"+kv.first]->Fill(HT_trigger, deepESM_val, weight);
                }
            }
        }

        // Fill histos for deepESM and Fisher training
        if(passBaseline1l_Good)
        {
            my_histos["fwm2_top6"]->Fill(fwm2_top6, weight);
            my_histos["fwm3_top6"]->Fill(fwm3_top6, weight);
            my_histos["fwm4_top6"]->Fill(fwm4_top6, weight);
            my_histos["fwm5_top6"]->Fill(fwm5_top6, weight);
            my_histos["fwm6_top6"]->Fill(fwm6_top6, weight);
            my_histos["fwm7_top6"]->Fill(fwm7_top6, weight);
            my_histos["fwm8_top6"]->Fill(fwm8_top6, weight);
            my_histos["fwm9_top6"]->Fill(fwm9_top6, weight);
            my_histos["fwm10_top6"]->Fill(fwm10_top6, weight);
            my_histos["jmt_ev0_top6"]->Fill(jmt_ev0_top6, weight);
            my_histos["jmt_ev1_top6"]->Fill(jmt_ev1_top6, weight);
            my_histos["jmt_ev2_top6"]->Fill(jmt_ev2_top6, weight);

            for(unsigned int i = 0; i < Jets_cm_top6.size(); i++)
            {
                my_histos["Jet_cm_pt_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Pt()), weight);
                my_histos["Jet_cm_eta_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Eta()), weight);
                my_histos["Jet_cm_phi_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Phi()), weight);
                my_histos["Jet_cm_m_"+std::to_string(i+1)+"_1l_ge7j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).M()), weight);

                my_2d_histos["Jet_cm_pt_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, Jets_cm_top6.at(i).Pt(), weight);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, Jets_cm_top6.at(i).Eta(), weight);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, Jets_cm_top6.at(i).Phi(), weight);
                my_2d_histos["Jet_cm_m_"+std::to_string(i+1)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, Jets_cm_top6.at(i).M(), weight);

                my_2d_histos["Jet_cm_pt_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Pt(), weight);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Eta(), weight);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).Phi(), weight);
                my_2d_histos["Jet_cm_m_"+std::to_string(i+1)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, Jets_cm_top6.at(i).M(), weight);
            }
        }
    } // end of event loop
}

void Analyze1Lep::WriteHistos(TFile* outfile)
{
    outfile->cd();
    
    for(const auto& p : my_histos) 
    {
        p.second->Write();
    }
    
    for(const auto& p : my_2d_histos) 
    {
        p.second->Write();
    }

    for(const auto& p : my_tp_histos) 
    {
        p.second->Write();
    }
    
    for(const auto& p : my_efficiencies) 
    {
        p.second->Write();
    }    
}
