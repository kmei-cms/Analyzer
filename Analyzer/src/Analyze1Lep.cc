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

    int fB = 200;

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("h_met",      std::make_shared<TH1D>("h_met",     "h_met",      20,  0,    200  ) );
    my_histos.emplace("h_bdt",      std::make_shared<TH1D>("h_bdt",     "h_bdt",      40, -0.5,    0.5) );
    my_histos.emplace("h_fisher",   std::make_shared<TH1D>("h_fisher",  "h_fisher",   fB, -0.5,    0.5) );
    my_histos.emplace("h_deepESM",  std::make_shared<TH1D>("h_deepESM", "h_deepESM",  fB,  0,      1) );
    my_histos.emplace("h_njets",    std::make_shared<TH1D>("h_njets",   "h_njets",    20,  0,     20  ) );
    my_histos.emplace("h_nb",       std::make_shared<TH1D>("h_nb",      "h_nb",       10,  0,     10  ) );
    my_histos.emplace("h_ntops",    std::make_shared<TH1D>("h_ntops",   "h_ntops",    10,  0,     10  ) );
    my_histos.emplace("h_ntops_j1", std::make_shared<TH1D>("h_ntops_j1","h_ntops_j1", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j2", std::make_shared<TH1D>("h_ntops_j2","h_ntops_j2", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j3", std::make_shared<TH1D>("h_ntops_j3","h_ntops_j3", 10,  0,     10  ) );

    my_histos.emplace("blind_met",      std::make_shared<TH1D>("blind_met",     "blind_met",      20,  0,    200  ) );
    my_histos.emplace("blind_bdt",      std::make_shared<TH1D>("blind_bdt",     "blind_bdt",      40, -0.5,    0.5) );
    my_histos.emplace("blind_fisher",   std::make_shared<TH1D>("blind_fisher",  "blind_fisher",   fB, -0.5,    0.5) );
    my_histos.emplace("blind_deepESM",  std::make_shared<TH1D>("blind_deepESM", "blind_deepESM",  fB,  0,      1) );
    my_histos.emplace("blind_njets",    std::make_shared<TH1D>("blind_njets",   "blind_njets",    20,  0,     20  ) );
    my_histos.emplace("blind_nb",       std::make_shared<TH1D>("blind_nb",      "blind_nb",       10,  0,     10  ) );
    my_histos.emplace("blind_ntops",    std::make_shared<TH1D>("blind_ntops",   "blind_ntops",    10,  0,     10  ) );
    my_histos.emplace("blind_ntops_j1", std::make_shared<TH1D>("blind_ntops_j1","blind_ntops_j1", 10,  0,     10  ) );
    my_histos.emplace("blind_ntops_j2", std::make_shared<TH1D>("blind_ntops_j2","blind_ntops_j2", 10,  0,     10  ) );
    my_histos.emplace("blind_ntops_j3", std::make_shared<TH1D>("blind_ntops_j3","blind_ntops_j3", 10,  0,     10  ) );

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

    for(auto& mycut : cutMap)
    {
        my_histos.emplace("h_njets_1l_"+mycut.first, std::make_shared<TH1D>(("h_njets_1l_"+mycut.first).c_str(),("h_njets_1l_"+mycut.first).c_str(), 20, 0, 20));
        my_histos.emplace("h_ntops_1l_"+mycut.first, std::make_shared<TH1D>(("h_ntops_1l_"+mycut.first).c_str(),("h_ntops_1l_"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("h_nb_1l_"+mycut.first,    std::make_shared<TH1D>(("h_nb_1l_"+mycut.first).c_str(),("h_nb_1l_"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("h_bdt_1l_"+mycut.first,   std::make_shared<TH1D>(("h_bdt_1l_"+mycut.first).c_str(),("h_bdt_1l_"+mycut.first).c_str(), 40, -0.5, 0.5));
        my_histos.emplace("h_fisher_1l_"+mycut.first,std::make_shared<TH1D>(("h_fisher_1l_"+mycut.first).c_str(),("h_fisher_1l_"+mycut.first).c_str(), fB, -0.5, 0.5));
        my_histos.emplace("h_deepESM_1l_"+mycut.first,std::make_shared<TH1D>(("h_deepESM_1l_"+mycut.first).c_str(),("h_deepESM_1l_"+mycut.first).c_str(), fB, 0.0, 1.0));
    
        my_2d_histos.emplace("h_njets_bdt_1l_"+mycut.first, std::make_shared<TH2D>(("h_njets_bdt_1l_"+mycut.first).c_str(),("h_njets_bdt_1l_"+mycut.first).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_fisher_1l_"+mycut.first, std::make_shared<TH2D>(("h_njets_fisher_1l_"+mycut.first).c_str(),("h_njets_fisher_1l_"+mycut.first).c_str(), 15, 0, 15, fB, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_deepESM_1l_"+mycut.first, std::make_shared<TH2D>(("h_njets_deepESM_1l_"+mycut.first).c_str(),("h_njets_deepESM_1l_"+mycut.first).c_str(), 15, 0, 15, fB, 0.0, 1.0));
        my_tp_histos.emplace("hTp_njets_fisher_1l_"+mycut.first, std::make_shared<TProfile>(("hTp_njets_fisher_1l_"+mycut.first).c_str(),("hTp_njets_fisher_1l_"+mycut.first).c_str(), 15, 0, 15, -0.5, 0.5));
        my_tp_histos.emplace("hTp_njets_deepESM_1l_"+mycut.first, std::make_shared<TProfile>(("hTp_njets_deepESM_1l_"+mycut.first).c_str(),("hTp_njets_deepESM_1l_"+mycut.first).c_str(), 15, 0, 15, 0.0, 1.0));

        my_histos.emplace("blind_njets_1l_"+mycut.first, std::make_shared<TH1D>(("blind_njets_1l_"+mycut.first).c_str(),("blind_njets_1l_"+mycut.first).c_str(), 20, 0, 20));
        my_histos.emplace("blind_ntops_1l_"+mycut.first, std::make_shared<TH1D>(("blind_ntops_1l_"+mycut.first).c_str(),("blind_ntops_1l_"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("blind_nb_1l_"+mycut.first,    std::make_shared<TH1D>(("blind_nb_1l_"+mycut.first).c_str(),("blind_nb_1l_"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("blind_bdt_1l_"+mycut.first,   std::make_shared<TH1D>(("blind_bdt_1l_"+mycut.first).c_str(),("blind_bdt_1l_"+mycut.first).c_str(), 40, -0.5, 0.5));
        my_histos.emplace("blind_fisher_1l_"+mycut.first,std::make_shared<TH1D>(("blind_fisher_1l_"+mycut.first).c_str(),("blind_fisher_1l_"+mycut.first).c_str(), fB, -0.5, 0.5));
        my_histos.emplace("blind_deepESM_1l_"+mycut.first,std::make_shared<TH1D>(("blind_deepESM_1l_"+mycut.first).c_str(),("blind_deepESM_1l_"+mycut.first).c_str(), fB, 0.0, 1.0));
    
        my_2d_histos.emplace("blind_njets_bdt_1l_"+mycut.first, std::make_shared<TH2D>(("blind_njets_bdt_1l_"+mycut.first).c_str(),("blind_njets_bdt_1l_"+mycut.first).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("blind_njets_fisher_1l_"+mycut.first, std::make_shared<TH2D>(("blind_njets_fisher_1l_"+mycut.first).c_str(),("blind_njets_fisher_1l_"+mycut.first).c_str(), 15, 0, 15, fB, -0.5, 0.5));
        my_2d_histos.emplace("blind_njets_deepESM_1l_"+mycut.first, std::make_shared<TH2D>(("blind_njets_deepESM_1l_"+mycut.first).c_str(),("blind_njets_deepESM_1l_"+mycut.first).c_str(), 15, 0, 15, fB, 0.0, 1.0));
        my_tp_histos.emplace("blindTp_njets_fisher_1l_"+mycut.first, std::make_shared<TProfile>(("blindTp_njets_fisher_1l_"+mycut.first).c_str(),("blindTp_njets_fisher_1l_"+mycut.first).c_str(), 15, 0, 15, -0.5, 0.5));
        my_tp_histos.emplace("blindTp_njets_deepESM_1l_"+mycut.first, std::make_shared<TProfile>(("blindTp_njets_deepESM_1l_"+mycut.first).c_str(),("blindTp_njets_deepESM_1l_"+mycut.first).c_str(), 15, 0, 15, 0.0, 1.0));
    }

    // Cut flows
    //my_efficiencies.emplace("event_sel",       std::make_shared<TEfficiency>("event_sel","Event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    //my_efficiencies.emplace("event_sel_total", std::make_shared<TEfficiency>("event_sel_total","Total event selection efficiency;Cut;#epsilon",8,0,8));
    //
    //my_efficiencies.emplace("event_sel_1l",       std::make_shared<TEfficiency>("event_sel_1l","0 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    //my_efficiencies.emplace("event_sel_total_1l", std::make_shared<TEfficiency>("event_sel_total_1l","Total 0 lepton event selection efficiency;Cut;#epsilon",8,0,8));
    //
    //my_efficiencies.emplace("event_sel_1l",       std::make_shared<TEfficiency>("event_sel_1l","1 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    //my_efficiencies.emplace("event_sel_total_1l", std::make_shared<TEfficiency>("event_sel_total_1l","Total 1 lepton event selection efficiency;Cut;#epsilon",8,0,8));
    //
    //my_efficiencies.emplace("event_sel_2l",       std::make_shared<TEfficiency>("event_sel_2l","2 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    //my_efficiencies.emplace("event_sel_total_2l", std::make_shared<TEfficiency>("event_sel_total_2l","Total 2 lepton event selection efficiency;Cut;#epsilon",8,0,8));

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
        const auto& NJets_pt30           = tr.getVar<int>("NJets_pt30");
        const auto& NJets_pt45           = tr.getVar<int>("NJets_pt45");
        const auto& NBJets               = tr.getVar<int>("NBJets");
        const auto& NBJets_pt45          = tr.getVar<int>("NBJets_pt45");
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons");
        const auto& fisher_bin1          = tr.getVar<bool>("fisher_bin1");
        const auto& fisher_bin2          = tr.getVar<bool>("fisher_bin2");
        const auto& fisher_bin3          = tr.getVar<bool>("fisher_bin3");
        const auto& fisher_bin4          = tr.getVar<bool>("fisher_bin4");
        const auto& eventshape_bdt_val   = tr.getVar<double>("eventshape_bdt_val");
        const auto& fisher_val           = tr.getVar<double>("fisher_val");
        const auto& passTrigger          = tr.getVar<bool>("passTrigger");
        const auto& passMadHT            = tr.getVar<bool>("passMadHT");
        const auto& passBaseline1l_Good  = tr.getVar<bool>("passBaseline1l_Good");
        const auto& deepESM_val          = tr.getVar<double>("deepEventShape_val");
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
        // ------------------------
        // -- Print event number
        // -----------------------        
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------
        // -- Define weight
        // -----------------------

        double eventweight = 1.0;        
        // Weight from samples.cc
        //eventweight = weight;

        if(runtype == "MC"){
            const auto& Weight  = tr.getVar<double>("Weight");
            double lumi = 35900; // Lumi for 2016
            // Weight from NTuples            
            eventweight = lumi*Weight;
        }

        // -------------------------------
        // -- Define cuts
        // -------------------------------

        // Global cuts
        if ( !(passTrigger && passMadHT) ) continue;
        
        bool pass_1l              = NGoodLeptons==1;
        bool pass_njet_pt45       = NJets_pt45>=6;
        bool pass_njet_pt45_1btag = NBJets_pt45 >= 1;
        bool pass_njet_pt45_2btag = NBJets_pt45 >= 2;

        // ---------------------------
        // --    1 Top selection
        // ---------------------------
        //exactly selections
        bool pass_1t    = ntops==1;
        bool pass_1t_f1 = pass_1t && fisher_bin1, pass_1t_f2 = pass_1t && fisher_bin2, pass_1t_f3 = pass_1t && fisher_bin3, pass_1t_f4 = pass_1t && fisher_bin4;

        bool pass_1t1 = pass_1t && ntops_1jet==1;
        bool pass_1t2 = pass_1t && ntops_2jet==1;
        bool pass_1t3 = pass_1t && ntops_3jet==1;
        bool pass_1t2or3 = pass_1t2 or pass_1t3;

        bool pass_1t1_f1 = pass_1t1 && fisher_bin1, pass_1t1_f2 = pass_1t1 && fisher_bin2, pass_1t1_f3 = pass_1t1 && fisher_bin3, pass_1t1_f4 = pass_1t1 && fisher_bin4;
        bool pass_1t2_f1 = pass_1t2 && fisher_bin1, pass_1t2_f2 = pass_1t2 && fisher_bin2, pass_1t2_f3 = pass_1t2 && fisher_bin3, pass_1t2_f4 = pass_1t2 && fisher_bin4;
        bool pass_1t3_f1 = pass_1t3 && fisher_bin1, pass_1t3_f2 = pass_1t3 && fisher_bin2, pass_1t3_f3 = pass_1t3 && fisher_bin3, pass_1t3_f4 = pass_1t3 && fisher_bin4;
        bool pass_1t2or3_f1 = pass_1t2or3 && fisher_bin1, pass_1t2or3_f2 = pass_1t2or3 && fisher_bin2, pass_1t2or3_f3 = pass_1t2or3 && fisher_bin3, pass_1t2or3_f4 = pass_1t2or3 && fisher_bin4;

        //at least selections
        bool pass_ge1t    = ntops>=1;
        bool pass_ge1t_f1 = pass_ge1t && fisher_bin1, pass_ge1t_f2 = pass_ge1t && fisher_bin2, pass_ge1t_f3 = pass_ge1t && fisher_bin3, pass_ge1t_f4 = pass_ge1t && fisher_bin4;

        bool pass_ge1t1 = pass_ge1t && ntops_1jet>=1;
        bool pass_ge1t2 = pass_ge1t && ntops_2jet>=1;
        bool pass_ge1t3 = pass_ge1t && ntops_3jet>=1;

        bool pass_ge1t1_f1 = pass_ge1t1 && fisher_bin1, pass_ge1t1_f2 = pass_ge1t1 && fisher_bin2, pass_ge1t1_f3 = pass_ge1t1 && fisher_bin3, pass_ge1t1_f4 = pass_ge1t1 && fisher_bin4;
        bool pass_ge1t2_f1 = pass_ge1t2 && fisher_bin1, pass_ge1t2_f2 = pass_ge1t2 && fisher_bin2, pass_ge1t2_f3 = pass_ge1t2 && fisher_bin3, pass_ge1t2_f4 = pass_ge1t2 && fisher_bin4;
        bool pass_ge1t3_f1 = pass_ge1t3 && fisher_bin1, pass_ge1t3_f2 = pass_ge1t3 && fisher_bin2, pass_ge1t3_f3 = pass_ge1t3 && fisher_bin3, pass_ge1t3_f4 = pass_ge1t3 && fisher_bin4;
        // ---------------------------
        // --    2 Top selection
        // ---------------------------

        // exactly selections
        bool pass_2t    = ntops==2;
        bool pass_2t_f1 = pass_2t && fisher_bin1, pass_2t_f2 = pass_2t && fisher_bin2, pass_2t_f3 = pass_2t && fisher_bin3, pass_2t_f4 = pass_2t && fisher_bin4;

        bool pass_2t11 = pass_2t && ntops_1jet==2;
        bool pass_2t12 = pass_2t && ntops_1jet==1 && ntops_2jet==1;
        bool pass_2t13 = pass_2t && ntops_1jet==1 && ntops_3jet==1;
        bool pass_2t22 = pass_2t && ntops_2jet==2;
        bool pass_2t23 = pass_2t && ntops_2jet==1 && ntops_3jet==1;
        bool pass_2t33 = pass_2t && ntops_3jet==2;

        bool pass_2t11_f1 = pass_2t11 && fisher_bin1, pass_2t11_f2 = pass_2t11 && fisher_bin2, pass_2t11_f3 = pass_2t11 && fisher_bin3, pass_2t11_f4 = pass_2t11 && fisher_bin4;
        bool pass_2t12_f1 = pass_2t12 && fisher_bin1, pass_2t12_f2 = pass_2t12 && fisher_bin2, pass_2t12_f3 = pass_2t12 && fisher_bin3, pass_2t12_f4 = pass_2t12 && fisher_bin4;
        bool pass_2t13_f1 = pass_2t13 && fisher_bin1, pass_2t13_f2 = pass_2t13 && fisher_bin2, pass_2t13_f3 = pass_2t13 && fisher_bin3, pass_2t13_f4 = pass_2t13 && fisher_bin4;
        bool pass_2t22_f1 = pass_2t22 && fisher_bin1, pass_2t22_f2 = pass_2t22 && fisher_bin2, pass_2t22_f3 = pass_2t22 && fisher_bin3, pass_2t22_f4 = pass_2t22 && fisher_bin4;
        bool pass_2t23_f1 = pass_2t23 && fisher_bin1, pass_2t23_f2 = pass_2t23 && fisher_bin2, pass_2t23_f3 = pass_2t23 && fisher_bin3, pass_2t23_f4 = pass_2t23 && fisher_bin4;
        bool pass_2t33_f1 = pass_2t33 && fisher_bin1, pass_2t33_f2 = pass_2t33 && fisher_bin2, pass_2t33_f3 = pass_2t33 && fisher_bin3, pass_2t33_f4 = pass_2t33 && fisher_bin4;

        // at least selections
        bool pass_ge2t    = ntops>=2;
        bool pass_ge2t_f1 = pass_ge2t && fisher_bin1, pass_ge2t_f2 = pass_ge2t && fisher_bin2, pass_ge2t_f3 = pass_ge2t && fisher_bin3, pass_ge2t_f4 = pass_ge2t && fisher_bin4;

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

        bool pass_ge2t11_f1 = pass_ge2t11 && fisher_bin1, pass_ge2t11_f2 = pass_ge2t11 && fisher_bin2, pass_ge2t11_f3 = pass_ge2t11 && fisher_bin3, pass_ge2t11_f4 = pass_ge2t11 && fisher_bin4;
        bool pass_ge2t12_f1 = pass_ge2t12 && fisher_bin1, pass_ge2t12_f2 = pass_ge2t12 && fisher_bin2, pass_ge2t12_f3 = pass_ge2t12 && fisher_bin3, pass_ge2t12_f4 = pass_ge2t12 && fisher_bin4;
        bool pass_ge2t13_f1 = pass_ge2t13 && fisher_bin1, pass_ge2t13_f2 = pass_ge2t13 && fisher_bin2, pass_ge2t13_f3 = pass_ge2t13 && fisher_bin3, pass_ge2t13_f4 = pass_ge2t13 && fisher_bin4;
        bool pass_ge2t22_f1 = pass_ge2t22 && fisher_bin1, pass_ge2t22_f2 = pass_ge2t22 && fisher_bin2, pass_ge2t22_f3 = pass_ge2t22 && fisher_bin3, pass_ge2t22_f4 = pass_ge2t22 && fisher_bin4;
        bool pass_ge2t23_f1 = pass_ge2t23 && fisher_bin1, pass_ge2t23_f2 = pass_ge2t23 && fisher_bin2, pass_ge2t23_f3 = pass_ge2t23 && fisher_bin3, pass_ge2t23_f4 = pass_ge2t23 && fisher_bin4;
        bool pass_ge2t33_f1 = pass_ge2t33 && fisher_bin1, pass_ge2t33_f2 = pass_ge2t33 && fisher_bin2, pass_ge2t33_f3 = pass_ge2t33 && fisher_bin3, pass_ge2t33_f4 = pass_ge2t33 && fisher_bin4;

        bool pass_ge2t11or12or13_f1 = pass_ge2t11or12or13 && fisher_bin1, pass_ge2t11or12or13_f2 = pass_ge2t11or12or13 && fisher_bin2, pass_ge2t11or12or13_f3 = pass_ge2t11or12or13 && fisher_bin3, pass_ge2t11or12or13_f4 = pass_ge2t11or12or13 && fisher_bin4;
        bool pass_ge2t22or23or33_f1 = pass_ge2t22or23or33 && fisher_bin1, pass_ge2t22or23or33_f2 = pass_ge2t22or23or33 && fisher_bin2, pass_ge2t22or23or33_f3 = pass_ge2t22or23or33 && fisher_bin3, pass_ge2t22or23or33_f4 = pass_ge2t22or23or33 && fisher_bin4;

        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map_1l 
        {
            {""                            , pass_1l                                                                  },
            {"ge6j"                        , pass_1l && pass_njet_pt45                                                },
            {"ge2b"                        , pass_1l && pass_njet_pt45_2btag                                          },
            {"1t"                          , pass_1l && pass_1t                                                       },
            {"2t"                          , pass_1l && pass_2t                                                       },
            {"ge1t"                        , pass_1l && pass_ge1t                                                     },
            {"ge2t"                        , pass_1l && pass_ge2t                                                     },
            {"ge6j_ge2b"                   , pass_1l && pass_njet_pt45 && pass_njet_pt45_2btag                        },
            {"ge6j_ge1b"                   , passBaseline1l_Good                                                      },                         
            {"ge6j_ge1b_f1"                , passBaseline1l_Good && fisher_bin1                                       },                         
            {"ge6j_ge1b_f2"                , passBaseline1l_Good && fisher_bin2                                       },                         
            {"ge6j_ge1b_f3"                , passBaseline1l_Good && fisher_bin3                                       },                         
            {"ge6j_ge1b_f4"                , passBaseline1l_Good && fisher_bin4                                       },                         
            {"6j_ge1b"                     , passBaseline1l_Good && NJets_pt45 == 6                                   },
            {"7j_ge1b"                     , passBaseline1l_Good && NJets_pt45 == 7                                   },
            {"8j_ge1b"                     , passBaseline1l_Good && NJets_pt45 == 8                                   },
            {"9j_ge1b"                     , passBaseline1l_Good && NJets_pt45 == 9                                   },
            {"10j_ge1b"                    , passBaseline1l_Good && NJets_pt45 == 10                                  },
            {"11j_ge1b"                    , passBaseline1l_Good && NJets_pt45 == 11                                  },
            {"12j_ge1b"                    , passBaseline1l_Good && NJets_pt45 == 12                                  },
            {"13j_ge1b"                    , passBaseline1l_Good && NJets_pt45 == 13                                  },
            {"14j_ge1b"                    , passBaseline1l_Good && NJets_pt45 == 14                                  },
            {"15j_ge1b"                    , passBaseline1l_Good && NJets_pt45 == 15                                  },
                                                 
            {"ge6j_ge1b_1t"                , passBaseline1l_Good && pass_1t                                           },
            {"ge6j_ge1b_1t_f1"             , passBaseline1l_Good && pass_1t_f1                                        },
            {"ge6j_ge1b_1t_f2"             , passBaseline1l_Good && pass_1t_f2                                        }, 
            {"ge6j_ge1b_1t_f3"             , passBaseline1l_Good && pass_1t_f3                                        }, 
            {"ge6j_ge1b_1t_f4"             , passBaseline1l_Good && pass_1t_f4                                        },
            {"ge6j_ge1b_ge1t"              , passBaseline1l_Good && pass_ge1t                                         },
            {"ge6j_ge1b_ge1t_f1"           , passBaseline1l_Good && pass_ge1t_f1                                      },
            {"ge6j_ge1b_ge1t_f2"           , passBaseline1l_Good && pass_ge1t_f2                                      }, 
            {"ge6j_ge1b_ge1t_f3"           , passBaseline1l_Good && pass_ge1t_f3                                      }, 
            {"ge6j_ge1b_ge1t_f4"           , passBaseline1l_Good && pass_ge1t_f4                                      },
                                                 
            {"ge6j_ge1b_1t1"               , passBaseline1l_Good && pass_1t1                                          },
            {"ge6j_ge1b_1t2"               , passBaseline1l_Good && pass_1t2                                          },
            {"ge6j_ge1b_1t3"               , passBaseline1l_Good && pass_1t3                                          },
            {"ge6j_ge1b_1t2or3"            , passBaseline1l_Good && pass_1t2or3                                       },
            {"ge6j_ge1b_1t1_f1"            , passBaseline1l_Good && pass_1t1_f1                                       },
            {"ge6j_ge1b_1t1_f2"            , passBaseline1l_Good && pass_1t1_f2                                       },
            {"ge6j_ge1b_1t1_f3"            , passBaseline1l_Good && pass_1t1_f3                                       },
            {"ge6j_ge1b_1t1_f4"            , passBaseline1l_Good && pass_1t1_f4                                       },
            {"ge6j_ge1b_1t2_f1"            , passBaseline1l_Good && pass_1t2_f1                                       },
            {"ge6j_ge1b_1t2_f2"            , passBaseline1l_Good && pass_1t2_f2                                       },
            {"ge6j_ge1b_1t2_f3"            , passBaseline1l_Good && pass_1t2_f3                                       },
            {"ge6j_ge1b_1t2_f4"            , passBaseline1l_Good && pass_1t2_f4                                       },
            {"ge6j_ge1b_1t3_f1"            , passBaseline1l_Good && pass_1t3_f1                                       },
            {"ge6j_ge1b_1t3_f2"            , passBaseline1l_Good && pass_1t3_f2                                       },
            {"ge6j_ge1b_1t3_f3"            , passBaseline1l_Good && pass_1t3_f3                                       },
            {"ge6j_ge1b_1t3_f4"            , passBaseline1l_Good && pass_1t3_f4                                       },
            {"ge6j_ge1b_1t2or3_f1"         , passBaseline1l_Good && pass_1t2or3_f1                                    },
            {"ge6j_ge1b_1t2or3_f2"         , passBaseline1l_Good && pass_1t2or3_f2                                    },
            {"ge6j_ge1b_1t2or3_f3"         , passBaseline1l_Good && pass_1t2or3_f3                                    },
            {"ge6j_ge1b_1t2or3_f4"         , passBaseline1l_Good && pass_1t2or3_f4                                    },
            {"ge6j_ge1b_ge1t1"             , passBaseline1l_Good && pass_ge1t1                                        },
            {"ge6j_ge1b_ge1t2"             , passBaseline1l_Good && pass_ge1t2                                        },
            {"ge6j_ge1b_ge1t3"             , passBaseline1l_Good && pass_ge1t3                                        },
            {"ge6j_ge1b_ge1t1_f1"          , passBaseline1l_Good && pass_ge1t1_f1                                     },
            {"ge6j_ge1b_ge1t1_f2"          , passBaseline1l_Good && pass_ge1t1_f2                                     },
            {"ge6j_ge1b_ge1t1_f3"          , passBaseline1l_Good && pass_ge1t1_f3                                     },
            {"ge6j_ge1b_ge1t1_f4"          , passBaseline1l_Good && pass_ge1t1_f4                                     },
            {"ge6j_ge1b_ge1t2_f1"          , passBaseline1l_Good && pass_ge1t2_f1                                     },
            {"ge6j_ge1b_ge1t2_f2"          , passBaseline1l_Good && pass_ge1t2_f2                                     },
            {"ge6j_ge1b_ge1t2_f3"          , passBaseline1l_Good && pass_ge1t2_f3                                     },
            {"ge6j_ge1b_ge1t2_f4"          , passBaseline1l_Good && pass_ge1t2_f4                                     },
            {"ge6j_ge1b_ge1t3_f1"          , passBaseline1l_Good && pass_ge1t3_f1                                     },
            {"ge6j_ge1b_ge1t3_f2"          , passBaseline1l_Good && pass_ge1t3_f2                                     },
            {"ge6j_ge1b_ge1t3_f3"          , passBaseline1l_Good && pass_ge1t3_f3                                     },
            {"ge6j_ge1b_ge1t3_f4"          , passBaseline1l_Good && pass_ge1t3_f4                                     },
                                                 
            {"ge6j_ge1b_2t"                , passBaseline1l_Good && pass_2t                                           },
            {"ge6j_ge1b_2t_f1"             , passBaseline1l_Good && pass_2t_f1                                        },
            {"ge6j_ge1b_2t_f2"             , passBaseline1l_Good && pass_2t_f2                                        }, 
            {"ge6j_ge1b_2t_f3"             , passBaseline1l_Good && pass_2t_f3                                        }, 
            {"ge6j_ge1b_2t_f4"             , passBaseline1l_Good && pass_2t_f4                                        },
            {"ge6j_ge1b_ge2t"              , passBaseline1l_Good && pass_ge2t                                         },
            {"ge6j_ge1b_ge2t_f1"           , passBaseline1l_Good && pass_ge2t_f1                                      },
            {"ge6j_ge1b_ge2t_f2"           , passBaseline1l_Good && pass_ge2t_f2                                      }, 
            {"ge6j_ge1b_ge2t_f3"           , passBaseline1l_Good && pass_ge2t_f3                                      }, 
            {"ge6j_ge1b_ge2t_f4"           , passBaseline1l_Good && pass_ge2t_f4                                      },
                                                 
            {"ge6j_ge1b_2t11"              , passBaseline1l_Good && pass_2t11                                         },
            {"ge6j_ge1b_2t12"              , passBaseline1l_Good && pass_2t12                                         },
            {"ge6j_ge1b_2t13"              , passBaseline1l_Good && pass_2t13                                         },
            {"ge6j_ge1b_2t22"              , passBaseline1l_Good && pass_2t22                                         },
            {"ge6j_ge1b_2t23"              , passBaseline1l_Good && pass_2t23                                         },
            {"ge6j_ge1b_2t33"              , passBaseline1l_Good && pass_2t33                                         },
            {"ge6j_ge1b_2t11_f1"           , passBaseline1l_Good && pass_2t11_f1                                      },
            {"ge6j_ge1b_2t11_f2"           , passBaseline1l_Good && pass_2t11_f2                                      },
            {"ge6j_ge1b_2t11_f3"           , passBaseline1l_Good && pass_2t11_f3                                      },
            {"ge6j_ge1b_2t11_f4"           , passBaseline1l_Good && pass_2t11_f4                                      },
            {"ge6j_ge1b_2t12_f1"           , passBaseline1l_Good && pass_2t12_f1                                      },
            {"ge6j_ge1b_2t12_f2"           , passBaseline1l_Good && pass_2t12_f2                                      },
            {"ge6j_ge1b_2t12_f3"           , passBaseline1l_Good && pass_2t12_f3                                      },
            {"ge6j_ge1b_2t12_f4"           , passBaseline1l_Good && pass_2t12_f4                                      },
            {"ge6j_ge1b_2t13_f1"           , passBaseline1l_Good && pass_2t13_f1                                      },
            {"ge6j_ge1b_2t13_f2"           , passBaseline1l_Good && pass_2t13_f2                                      },
            {"ge6j_ge1b_2t13_f3"           , passBaseline1l_Good && pass_2t13_f3                                      },
            {"ge6j_ge1b_2t13_f4"           , passBaseline1l_Good && pass_2t13_f4                                      },
            {"ge6j_ge1b_2t22_f1"           , passBaseline1l_Good && pass_2t22_f1                                      },
            {"ge6j_ge1b_2t22_f2"           , passBaseline1l_Good && pass_2t22_f2                                      },
            {"ge6j_ge1b_2t22_f3"           , passBaseline1l_Good && pass_2t22_f3                                      },
            {"ge6j_ge1b_2t22_f4"           , passBaseline1l_Good && pass_2t22_f4                                      },
            {"ge6j_ge1b_2t23_f1"           , passBaseline1l_Good && pass_2t23_f1                                      },
            {"ge6j_ge1b_2t23_f2"           , passBaseline1l_Good && pass_2t23_f2                                      },
            {"ge6j_ge1b_2t23_f3"           , passBaseline1l_Good && pass_2t23_f3                                      },
            {"ge6j_ge1b_2t23_f4"           , passBaseline1l_Good && pass_2t23_f4                                      },
            {"ge6j_ge1b_2t33_f1"           , passBaseline1l_Good && pass_2t33_f1                                      },
            {"ge6j_ge1b_2t33_f2"           , passBaseline1l_Good && pass_2t33_f2                                      },
            {"ge6j_ge1b_2t33_f3"           , passBaseline1l_Good && pass_2t33_f3                                      },
            {"ge6j_ge1b_2t33_f4"           , passBaseline1l_Good && pass_2t33_f4                                      },
            {"ge6j_ge1b_ge2t11"            , passBaseline1l_Good && pass_ge2t11                                       },
            {"ge6j_ge1b_ge2t12"            , passBaseline1l_Good && pass_ge2t12                                       },
            {"ge6j_ge1b_ge2t13"            , passBaseline1l_Good && pass_ge2t13                                       },
            {"ge6j_ge1b_ge2t11or12or13"    , passBaseline1l_Good && pass_ge2t11or12or13                               },
            {"ge6j_ge1b_ge2t22"            , passBaseline1l_Good && pass_ge2t22                                       },
            {"ge6j_ge1b_ge2t23"            , passBaseline1l_Good && pass_ge2t23                                       },
            {"ge6j_ge1b_ge2t33"            , passBaseline1l_Good && pass_ge2t33                                       },
            {"ge6j_ge1b_ge2t22or23or33"    , passBaseline1l_Good && pass_ge2t22or23or33                               },
            {"ge6j_ge1b_ge2t11_f1"         , passBaseline1l_Good && pass_ge2t11_f1                                    },
            {"ge6j_ge1b_ge2t11_f2"         , passBaseline1l_Good && pass_ge2t11_f2                                    },
            {"ge6j_ge1b_ge2t11_f3"         , passBaseline1l_Good && pass_ge2t11_f3                                    },
            {"ge6j_ge1b_ge2t11_f4"         , passBaseline1l_Good && pass_ge2t11_f4                                    },
            {"ge6j_ge1b_ge2t12_f1"         , passBaseline1l_Good && pass_ge2t12_f1                                    },
            {"ge6j_ge1b_ge2t12_f2"         , passBaseline1l_Good && pass_ge2t12_f2                                    },
            {"ge6j_ge1b_ge2t12_f3"         , passBaseline1l_Good && pass_ge2t12_f3                                    },
            {"ge6j_ge1b_ge2t12_f4"         , passBaseline1l_Good && pass_ge2t12_f4                                    },
            {"ge6j_ge1b_ge2t13_f1"         , passBaseline1l_Good && pass_ge2t13_f1                                    },
            {"ge6j_ge1b_ge2t13_f2"         , passBaseline1l_Good && pass_ge2t13_f2                                    },
            {"ge6j_ge1b_ge2t13_f3"         , passBaseline1l_Good && pass_ge2t13_f3                                    },
            {"ge6j_ge1b_ge2t13_f4"         , passBaseline1l_Good && pass_ge2t13_f4                                    },
            {"ge6j_ge1b_ge2t11or12or13_f1" , passBaseline1l_Good && pass_ge2t11or12or13_f1                            },
            {"ge6j_ge1b_ge2t11or12or13_f2" , passBaseline1l_Good && pass_ge2t11or12or13_f2                            },
            {"ge6j_ge1b_ge2t11or12or13_f3" , passBaseline1l_Good && pass_ge2t11or12or13_f3                            },
            {"ge6j_ge1b_ge2t11or12or13_f4" , passBaseline1l_Good && pass_ge2t11or12or13_f4                            },
            {"ge6j_ge1b_ge2t22_f1"         , passBaseline1l_Good && pass_ge2t22_f1                                    },
            {"ge6j_ge1b_ge2t22_f2"         , passBaseline1l_Good && pass_ge2t22_f2                                    },
            {"ge6j_ge1b_ge2t22_f3"         , passBaseline1l_Good && pass_ge2t22_f3                                    },
            {"ge6j_ge1b_ge2t22_f4"         , passBaseline1l_Good && pass_ge2t22_f4                                    },
            {"ge6j_ge1b_ge2t23_f1"         , passBaseline1l_Good && pass_ge2t23_f1                                    },
            {"ge6j_ge1b_ge2t23_f2"         , passBaseline1l_Good && pass_ge2t23_f2                                    },
            {"ge6j_ge1b_ge2t23_f3"         , passBaseline1l_Good && pass_ge2t23_f3                                    },
            {"ge6j_ge1b_ge2t23_f4"         , passBaseline1l_Good && pass_ge2t23_f4                                    },
            {"ge6j_ge1b_ge2t33_f1"         , passBaseline1l_Good && pass_ge2t33_f1                                    },
            {"ge6j_ge1b_ge2t33_f2"         , passBaseline1l_Good && pass_ge2t33_f2                                    },
            {"ge6j_ge1b_ge2t33_f3"         , passBaseline1l_Good && pass_ge2t33_f3                                    },
            {"ge6j_ge1b_ge2t33_f4"         , passBaseline1l_Good && pass_ge2t33_f4                                    },
            {"ge6j_ge1b_ge2t22or23or33_f1" , passBaseline1l_Good && pass_ge2t22or23or33_f1                            },
            {"ge6j_ge1b_ge2t22or23or33_f2" , passBaseline1l_Good && pass_ge2t22or23or33_f2                            },
            {"ge6j_ge1b_ge2t22or23or33_f3" , passBaseline1l_Good && pass_ge2t22or23or33_f3                            },
            {"ge6j_ge1b_ge2t22or23or33_f4" , passBaseline1l_Good && pass_ge2t22or23or33_f4                            },
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map_1l);
            initHistos = true;
        }

        for(auto& kv : cut_map_1l)
        {
            if(kv.second)
            {
                my_histos["h_njets_1l_"   +kv.first]->Fill(NJets_pt30, eventweight);
                my_histos["h_ntops_1l_"   +kv.first]->Fill(ntops, eventweight);
                my_histos["h_nb_1l_"      +kv.first]->Fill(NBJets, eventweight);
                my_histos["h_bdt_1l_"     +kv.first]->Fill(eventshape_bdt_val, eventweight);
                my_histos["h_fisher_1l_"  +kv.first]->Fill(fisher_val, eventweight);
                my_histos["h_deepESM_1l_"  +kv.first]->Fill(deepESM_val, eventweight);

                my_2d_histos["h_njets_bdt_1l_"+kv.first]->Fill(NJets_pt30, eventshape_bdt_val, eventweight);
                my_2d_histos["h_njets_fisher_1l_"+kv.first]->Fill(NJets_pt30, fisher_val, eventweight);
                my_2d_histos["h_njets_deepESM_1l_"+kv.first]->Fill(NJets_pt30, deepESM_val, eventweight);
                my_tp_histos["hTp_njets_fisher_1l_"+kv.first]->Fill(NJets_pt30, fisher_val, eventweight);
                my_tp_histos["hTp_njets_deepESM_1l_"+kv.first]->Fill(NJets_pt30, deepESM_val, eventweight);

                if ( NJets_pt30 < 9 )
                {
                    my_histos["blind_njets_1l_"   +kv.first]->Fill(NJets_pt30, eventweight);
                    my_histos["blind_ntops_1l_"   +kv.first]->Fill(ntops, eventweight);
                    my_histos["blind_nb_1l_"      +kv.first]->Fill(NBJets, eventweight);
                    my_histos["blind_bdt_1l_"     +kv.first]->Fill(eventshape_bdt_val, eventweight);
                    my_histos["blind_fisher_1l_"  +kv.first]->Fill(fisher_val, eventweight);
                    my_histos["blind_deepESM_1l_"  +kv.first]->Fill(deepESM_val, eventweight);
                
                    my_2d_histos["blind_njets_bdt_1l_"+kv.first]->Fill(NJets_pt30, eventshape_bdt_val, eventweight);
                    my_2d_histos["blind_njets_fisher_1l_"+kv.first]->Fill(NJets_pt30, fisher_val, eventweight);
                    my_2d_histos["blind_njets_deepESM_1l_"+kv.first]->Fill(NJets_pt30, deepESM_val, eventweight);
                    my_tp_histos["blindTp_njets_fisher_1l_"+kv.first]->Fill(NJets_pt30, fisher_val, eventweight);
                    my_tp_histos["blindTp_njets_deepESM_1l_"+kv.first]->Fill(NJets_pt30, deepESM_val, eventweight);
                }
            }
        }

        // Fill event selection efficiencies
        //my_efficiencies["event_sel_total"]->Fill(true,0);
        //my_efficiencies["event_sel_total"]->Fill(HT_trigger>500,1);
        //my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 ,2);
        //my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 ,3);
        //my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 ,4);
        //my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 ,5);
        //my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 && ntops>1 ,6);
        //my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 && ntops>1 && NJets_pt30>=8 ,7);
        //
        //my_efficiencies["event_sel"]->Fill(true,0);
        //my_efficiencies["event_sel"]->Fill(HT_trigger>500,1);
        //if(HT_trigger>500)
        //{
        //    my_efficiencies["event_sel"]->Fill(NJets_pt45>=6,2);
        //    if (NJets_pt45>=6)
        //    {
        //        my_efficiencies["event_sel"]->Fill(NBJets_pt45>0,3);
        //        if (NBJets_pt45>0)
        //        {
        //            my_efficiencies["event_sel"]->Fill(ntops>0,4);
        //            if (ntops>0)
        //            {
        //                my_efficiencies["event_sel"]->Fill(NBJets_pt45>1,5);
        //                if (NBJets_pt45>1)
        //                {
        //                    my_efficiencies["event_sel"]->Fill(ntops>1,6);
        //                    if (ntops>1)
        //                    {
        //                        my_efficiencies["event_sel"]->Fill(NJets_pt30>=8,7);
        //                    }
        //                }
        //            }
        //        }
        //    }
        //}
        
        // No local cuts applied here
        my_histos["h_met"     ]->Fill(MET, eventweight);
        my_histos["h_bdt"     ]->Fill(eventshape_bdt_val, eventweight);
        my_histos["h_fisher"  ]->Fill(fisher_val, eventweight);
        my_histos["h_deepESM"  ]->Fill(deepESM_val, eventweight);
        my_histos["h_njets"   ]->Fill(NJets_pt30, eventweight);
        my_histos["h_nb"      ]->Fill(NBJets, eventweight);
        my_histos["h_ntops"   ]->Fill(ntops, eventweight);        
        my_histos["h_ntops_j1"]->Fill(ntops_1jet, eventweight);
        my_histos["h_ntops_j2"]->Fill(ntops_2jet, eventweight);
        my_histos["h_ntops_j3"]->Fill(ntops_3jet, eventweight);        

        // Do the blind cuts, for now, for data MC plots
        if ( NJets_pt30 < 9 )
        {
            my_histos["blind_met"     ]->Fill(MET, eventweight);
            my_histos["blind_bdt"     ]->Fill(eventshape_bdt_val, eventweight);
            my_histos["blind_fisher"  ]->Fill(fisher_val, eventweight);
            my_histos["blind_deepESM"  ]->Fill(deepESM_val, eventweight);
            my_histos["blind_njets"   ]->Fill(NJets_pt30, eventweight);
            my_histos["blind_nb"      ]->Fill(NBJets, eventweight);
            my_histos["blind_ntops"   ]->Fill(ntops, eventweight);        
            my_histos["blind_ntops_j1"]->Fill(ntops_1jet, eventweight);
            my_histos["blind_ntops_j2"]->Fill(ntops_2jet, eventweight);
            my_histos["blind_ntops_j3"]->Fill(ntops_3jet, eventweight);        
        }

        // Fill histos for deepESM and Fisher training
        if(passBaseline1l_Good)
        {
            my_histos["fwm2_top6"]->Fill(fwm2_top6, eventweight);
            my_histos["fwm3_top6"]->Fill(fwm3_top6, eventweight);
            my_histos["fwm4_top6"]->Fill(fwm4_top6, eventweight);
            my_histos["fwm5_top6"]->Fill(fwm5_top6, eventweight);
            my_histos["fwm6_top6"]->Fill(fwm6_top6, eventweight);
            my_histos["fwm7_top6"]->Fill(fwm7_top6, eventweight);
            my_histos["fwm8_top6"]->Fill(fwm8_top6, eventweight);
            my_histos["fwm9_top6"]->Fill(fwm9_top6, eventweight);
            my_histos["fwm10_top6"]->Fill(fwm10_top6, eventweight);
            my_histos["jmt_ev0_top6"]->Fill(jmt_ev0_top6, eventweight);
            my_histos["jmt_ev1_top6"]->Fill(jmt_ev1_top6, eventweight);
            my_histos["jmt_ev2_top6"]->Fill(jmt_ev2_top6, eventweight);
        }

    } // end of event loop

}

void Analyze1Lep::WriteHistos(TFile* outfile)
{
    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }

    for (const auto &p : my_tp_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }
    
}
