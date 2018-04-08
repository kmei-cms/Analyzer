#define Analyze0Lep_cxx
#include "Analyzer/Analyzer/include/Analyze0Lep.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

void Analyze0Lep::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("h_met",      std::make_shared<TH1D>("h_met",     "h_met",      20,  0,    200  ) );
    my_histos.emplace("h_ht",       std::make_shared<TH1D>("h_ht",      "h_ht",       60,  0,   3000  ) );
    my_histos.emplace("h_bdt",      std::make_shared<TH1D>("h_bdt",     "h_bdt",      40, -0.5,    0.5) );
    my_histos.emplace("h_fisher",   std::make_shared<TH1D>("h_fisher",  "h_fisher",   40, -0.5,    0.5) );
    my_histos.emplace("h_njets",    std::make_shared<TH1D>("h_njets",   "h_njets",    20,  0,     20  ) );
    my_histos.emplace("h_nb",       std::make_shared<TH1D>("h_nb",      "h_nb",       10,  0,     10  ) );
    my_histos.emplace("h_ntops",    std::make_shared<TH1D>("h_ntops",   "h_ntops",    10,  0,     10  ) );
    my_histos.emplace("h_ntops_j1", std::make_shared<TH1D>("h_ntops_j1","h_ntops_j1", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j2", std::make_shared<TH1D>("h_ntops_j2","h_ntops_j2", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j3", std::make_shared<TH1D>("h_ntops_j3","h_ntops_j3", 10,  0,     10  ) );
    
    std::vector<std::string> mycuts_0l 
    {
        ""                     ,
        "g6j"                  ,
        "HT500"                ,
        "g2b"                  ,
        "1t"                   ,
        "2t"                   ,
        "g6j_HT500"            ,
        "g6j_HT500_g1b"        ,
        "g6j_HT500_g2b"        ,

        "g6j_HT500_g2b_1t"     ,
        "g6j_HT500_g2b_1t_f1"  , "g6j_HT500_g2b_1t_f2"  , "g6j_HT500_g2b_1t_f3"  , "g6j_HT500_g2b_1t_f4"  ,

        "g6j_HT500_g2b_1t1"    , "g6j_HT500_g2b_1t2"    , "g6j_HT500_g2b_1t3"    ,
        "g6j_HT500_g2b_1t1_f1" , "g6j_HT500_g2b_1t1_f2" , "g6j_HT500_g2b_1t1_f3" , "g6j_HT500_g2b_1t1_f4" ,
        "g6j_HT500_g2b_1t2_f1" , "g6j_HT500_g2b_1t2_f2" , "g6j_HT500_g2b_1t2_f3" , "g6j_HT500_g2b_1t2_f4" ,
        "g6j_HT500_g2b_1t3_f1" , "g6j_HT500_g2b_1t3_f2" , "g6j_HT500_g2b_1t3_f3" , "g6j_HT500_g2b_1t3_f4" ,

        "g6j_HT500_g2b_2t"     ,
        "g6j_HT500_g2b_2t_f1"  , "g6j_HT500_g2b_2t_f2"  , "g6j_HT500_g2b_2t_f3"  , "g6j_HT500_g2b_2t_f4"  ,

        "g6j_HT500_g2b_2t11"   , "g6j_HT500_g2b_2t12"   , "g6j_HT500_g2b_2t13"   , "g6j_HT500_g2b_2t22"   , "g6j_HT500_g2b_2t23", "g6j_HT500_g2b_2t33",
        "g6j_HT500_g2b_2t11_f1", "g6j_HT500_g2b_2t11_f2", "g6j_HT500_g2b_2t11_f3", "g6j_HT500_g2b_2t11_f4",
        "g6j_HT500_g2b_2t12_f1", "g6j_HT500_g2b_2t12_f2", "g6j_HT500_g2b_2t12_f3", "g6j_HT500_g2b_2t12_f4",
        "g6j_HT500_g2b_2t13_f1", "g6j_HT500_g2b_2t13_f2", "g6j_HT500_g2b_2t13_f3", "g6j_HT500_g2b_2t13_f4",
        "g6j_HT500_g2b_2t22_f1", "g6j_HT500_g2b_2t22_f2", "g6j_HT500_g2b_2t22_f3", "g6j_HT500_g2b_2t22_f4",
        "g6j_HT500_g2b_2t23_f1", "g6j_HT500_g2b_2t23_f2", "g6j_HT500_g2b_2t23_f3", "g6j_HT500_g2b_2t23_f4",
        "g6j_HT500_g2b_2t33_f1", "g6j_HT500_g2b_2t33_f2", "g6j_HT500_g2b_2t33_f3", "g6j_HT500_g2b_2t33_f4",
    };

    for(std::string mycut : mycuts_0l)
    {
        my_histos.emplace("h_njets_0l_"+mycut, std::make_shared<TH1D>(("h_njets_0l_"+mycut).c_str(),("h_njets_0l_"+mycut).c_str(), 20, 0, 20));
        my_histos.emplace("h_ntops_0l_"+mycut, std::make_shared<TH1D>(("h_ntops_0l_"+mycut).c_str(),("h_ntops_0l_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_nb_0l_"+mycut,    std::make_shared<TH1D>(("h_nb_0l_"+mycut).c_str(),("h_nb_0l_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_0l_"+mycut,    std::make_shared<TH1D>(("h_HT_0l_"+mycut).c_str(),("h_HT_0l_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_bdt_0l_"+mycut,   std::make_shared<TH1D>(("h_bdt_0l_"+mycut).c_str(),("h_bdt_0l_"+mycut).c_str(), 40, -0.5, 0.5));
        my_histos.emplace("h_fisher_0l_"+mycut,std::make_shared<TH1D>(("h_fisher_0l_"+mycut).c_str(),("h_fisher_0l_"+mycut).c_str(), 40, -0.5, 0.5));

        my_2d_histos.emplace("h_njets_bdt_0l_"+mycut, std::make_shared<TH2D>(("h_njets_bdt_0l_"+mycut).c_str(),("h_njets_bdt_0l_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_fisher_0l_"+mycut, std::make_shared<TH2D>(("h_njets_fisher_0l_"+mycut).c_str(),("h_njets_fisher_0l_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
    }

    // Cut flows
    my_efficiencies.emplace("event_sel",       std::make_shared<TEfficiency>("event_sel","Event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total", std::make_shared<TEfficiency>("event_sel_total","Total event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_0l",       std::make_shared<TEfficiency>("event_sel_0l","0 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_0l", std::make_shared<TEfficiency>("event_sel_total_0l","Total 0 lepton event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_1l",       std::make_shared<TEfficiency>("event_sel_1l","1 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_1l", std::make_shared<TEfficiency>("event_sel_total_1l","Total 1 lepton event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_2l",       std::make_shared<TEfficiency>("event_sel_2l","2 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_2l", std::make_shared<TEfficiency>("event_sel_total_2l","Total 2 lepton event selection efficiency;Cut;#epsilon",8,0,8));

}

void Analyze0Lep::Loop(NTupleReader& tr, double weight, int maxevents, std::string filetag, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& MET                = tr.getVar<double>("MET");
        const auto& HT                 = tr.getVar<double>("HT");
        const auto& ntops              = tr.getVar<int>("ntops");
        const auto& ntops_3jet         = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet         = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet         = tr.getVar<int>("ntops_1jet");
        const auto& runtype            = tr.getVar<std::string>("runtype");     
        const auto& TriggerNames       = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass        = tr.getVec<int>("TriggerPass");
        const auto& Jets               = tr.getVec<TLorentzVector>("Jets");
        const auto& NBJets             = tr.getVar<int>("NBJets");
        const auto& NBJets_pt45        = tr.getVar<int>("NBJets_pt45");
        const auto& NGoodLeptons       = tr.getVar<int>("NGoodLeptons");
        const auto& HT_trigger         = tr.getVar<double>("HT_trigger");
        const auto& fisher_bin1        = tr.getVar<bool>("fisher_bin1");
        const auto& fisher_bin2        = tr.getVar<bool>("fisher_bin2");
        const auto& fisher_bin3        = tr.getVar<bool>("fisher_bin3");
        const auto& fisher_bin4        = tr.getVar<bool>("fisher_bin4");
        //const auto& bdt_bin1           = tr.getVar<bool>("bdt_bin1");
        //const auto& bdt_bin2           = tr.getVar<bool>("bdt_bin2");
        //const auto& bdt_bin3           = tr.getVar<bool>("bdt_bin3");
        //const auto& bdt_bin4           = tr.getVar<bool>("bdt_bin4");
        const auto& eventshape_bdt_val = tr.getVar<double>("eventshape_bdt_val");
        const auto& fisher_val         = tr.getVar<double>("fisher_val");

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

        // ------------------------
        // -- MC dependent stuff
        // -----------------------

        if(runtype == "MC"){
            const auto& madHT   = tr.getVar<double>("madHT");
            const auto& Weight  = tr.getVar<double>("Weight");
            double lumi = 35900; // Lumi for 2016

            // Exclude events with MadGraph HT > 100 from the DY inclusive sample
            if(filetag == "DYJetsToLL_M-50_Incl" && madHT > 100) continue;

            // Weight from NTuples            
            eventweight = lumi*Weight;
        }
        
        // ------------------------------
        // -- Trigger for data
        // ------------------------------
        
        bool passTriggerAllHad   = PassTriggerAllHad(TriggerNames, TriggerPass);
        bool passTriggerMuon     = PassTriggerMuon(TriggerNames, TriggerPass);
        bool passTriggerElectron = PassTriggerElectron(TriggerNames, TriggerPass);
        if (runtype == "Data")
        {
            if (filetag == "Data_JetHT" && !passTriggerAllHad) continue;
            if (filetag == "Data_SingleMuon" && !passTriggerMuon) continue;
            if (filetag == "Data_SingleElectron" && !passTriggerElectron) continue;
        }
        
        // -------------------------------
        // -- Basic event selection stuff
        // -------------------------------
        
        // Count jets
        int rec_njet_pt45(0) ;
        int rec_njet_pt30(0) ;
        for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) 
        {
            TLorentzVector jlv( Jets.at(rji) ) ;
            if (abs(jlv.Eta()) > 2.4) continue;
            if ( jlv.Pt() > 30 )
            { 
                rec_njet_pt30++;
            }
            if ( jlv.Pt() > 45 ) 
            {
                rec_njet_pt45++ ;
            }
        }         
        
        // Define selections
        bool passBaseline0l = NGoodLeptons==0 && rec_njet_pt45>=6 && HT_trigger > 500 && NBJets_pt45 >= 2;
        if(runtype == "Data" && filetag == "Data_JetHT")
        {
            passBaseline0l = passBaseline0l && passTriggerAllHad;
        }
        
        bool pass_0l              = NGoodLeptons==0;
        bool pass_njet_pt45       = rec_njet_pt45>=6;
        bool pass_HT_trigger      = HT_trigger > 500;
        bool pass_njet_pt45_1btag = NBJets_pt45 >= 1;
        bool pass_njet_pt45_2btag = NBJets_pt45 >= 2;

        // 1 Top selection
        bool pass_1t    = ntops>=1;
        bool pass_1t_f1 = pass_1t && fisher_bin1, pass_1t_f2 = pass_1t && fisher_bin2, pass_1t_f3 = pass_1t && fisher_bin3, pass_1t_f4 = pass_1t && fisher_bin4;

        bool pass_1t1 = pass_1t && ntops_1jet>=1;
        bool pass_1t2 = pass_1t && ntops_2jet>=1;
        bool pass_1t3 = pass_1t && ntops_3jet>=1;

        bool pass_1t1_f1 = pass_1t1 && fisher_bin1, pass_1t1_f2 = pass_1t1 && fisher_bin2, pass_1t1_f3 = pass_1t1 && fisher_bin3, pass_1t1_f4 = pass_1t1 && fisher_bin4;
        bool pass_1t2_f1 = pass_1t2 && fisher_bin1, pass_1t2_f2 = pass_1t2 && fisher_bin2, pass_1t2_f3 = pass_1t2 && fisher_bin3, pass_1t2_f4 = pass_1t2 && fisher_bin4;
        bool pass_1t3_f1 = pass_1t3 && fisher_bin1, pass_1t3_f2 = pass_1t3 && fisher_bin2, pass_1t3_f3 = pass_1t3 && fisher_bin3, pass_1t3_f4 = pass_1t3 && fisher_bin4;

        // 2 Top selection
        bool pass_2t    = ntops>=2;
        bool pass_2t_f1 = pass_2t && fisher_bin1, pass_2t_f2 = pass_2t && fisher_bin2, pass_2t_f3 = pass_2t && fisher_bin3, pass_2t_f4 = pass_2t && fisher_bin4;

        bool pass_2t11 = pass_2t && ntops_1jet>=2;
        bool pass_2t12 = pass_2t && ntops_1jet>=1 && ntops_2jet>=1;
        bool pass_2t13 = pass_2t && ntops_1jet>=1 && ntops_3jet>=1;
        bool pass_2t22 = pass_2t && ntops_2jet>=2;
        bool pass_2t23 = pass_2t && ntops_2jet>=1 && ntops_3jet>=1;
        bool pass_2t33 = pass_2t && ntops_3jet>=2;

        bool pass_2t11_f1 = pass_2t11 && fisher_bin1, pass_2t11_f2 = pass_2t11 && fisher_bin2, pass_2t11_f3 = pass_2t11 && fisher_bin3, pass_2t11_f4 = pass_2t11 && fisher_bin4;
        bool pass_2t12_f1 = pass_2t12 && fisher_bin1, pass_2t12_f2 = pass_2t12 && fisher_bin2, pass_2t12_f3 = pass_2t12 && fisher_bin3, pass_2t12_f4 = pass_2t12 && fisher_bin4;
        bool pass_2t13_f1 = pass_2t13 && fisher_bin1, pass_2t13_f2 = pass_2t13 && fisher_bin2, pass_2t13_f3 = pass_2t13 && fisher_bin3, pass_2t13_f4 = pass_2t13 && fisher_bin4;
        bool pass_2t22_f1 = pass_2t22 && fisher_bin1, pass_2t22_f2 = pass_2t22 && fisher_bin2, pass_2t22_f3 = pass_2t22 && fisher_bin3, pass_2t22_f4 = pass_2t22 && fisher_bin4;
        bool pass_2t23_f1 = pass_2t23 && fisher_bin1, pass_2t23_f2 = pass_2t23 && fisher_bin2, pass_2t23_f3 = pass_2t23 && fisher_bin3, pass_2t23_f4 = pass_2t23 && fisher_bin4;
        bool pass_2t33_f1 = pass_2t33 && fisher_bin1, pass_2t33_f2 = pass_2t33 && fisher_bin2, pass_2t33_f3 = pass_2t33 && fisher_bin3, pass_2t33_f4 = pass_2t33 && fisher_bin4;

        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map_0l {
            {""                     , pass_0l                                                             },
            {"g6j"                  , pass_0l && pass_njet_pt45                                           },
            {"HT500"                , pass_0l && pass_HT_trigger                                          },
            {"g2b"                  , pass_0l && pass_njet_pt45_2btag                                     },
            {"1t"                   , pass_0l && pass_1t                                                  },
            {"2t"                   , pass_0l && pass_2t                                                  },
            {"g6j_HT500"            , pass_0l && pass_njet_pt45 && pass_HT_trigger                        },
            {"g6j_HT500_g1b"        , pass_0l && pass_njet_pt45 && pass_HT_trigger && pass_njet_pt45_1btag},
            {"g6j_HT500_g2b"        , passBaseline0l                                                      },                         

            {"g6j_HT500_g2b_1t"     , passBaseline0l && pass_1t                                           },
            {"g6j_HT500_g2b_1t_f1"  , passBaseline0l && pass_1t_f1                                        },
            {"g6j_HT500_g2b_1t_f2"  , passBaseline0l && pass_1t_f2                                        }, 
            {"g6j_HT500_g2b_1t_f3"  , passBaseline0l && pass_1t_f3                                        }, 
            {"g6j_HT500_g2b_1t_f4"  , passBaseline0l && pass_1t_f4                                        },

            {"g6j_HT500_g2b_1t1"    , passBaseline0l && pass_1t1                                          },
            {"g6j_HT500_g2b_1t2"    , passBaseline0l && pass_1t2                                          },
            {"g6j_HT500_g2b_1t3"    , passBaseline0l && pass_1t3                                          },
            {"g6j_HT500_g2b_1t1_f1" , passBaseline0l && pass_1t1_f1                                       },
            {"g6j_HT500_g2b_1t1_f2" , passBaseline0l && pass_1t1_f2                                       },
            {"g6j_HT500_g2b_1t1_f3" , passBaseline0l && pass_1t1_f3                                       },
            {"g6j_HT500_g2b_1t1_f4" , passBaseline0l && pass_1t1_f4                                       },
            {"g6j_HT500_g2b_1t2_f1" , passBaseline0l && pass_1t2_f1                                       },
            {"g6j_HT500_g2b_1t2_f2" , passBaseline0l && pass_1t2_f2                                       },
            {"g6j_HT500_g2b_1t2_f3" , passBaseline0l && pass_1t2_f3                                       },
            {"g6j_HT500_g2b_1t2_f4" , passBaseline0l && pass_1t2_f4                                       },
            {"g6j_HT500_g2b_1t3_f1" , passBaseline0l && pass_1t3_f1                                       },
            {"g6j_HT500_g2b_1t3_f2" , passBaseline0l && pass_1t3_f2                                       },
            {"g6j_HT500_g2b_1t3_f3" , passBaseline0l && pass_1t3_f3                                       },
            {"g6j_HT500_g2b_1t3_f4" , passBaseline0l && pass_1t3_f4                                       },

            {"g6j_HT500_g2b_2t"     , passBaseline0l && pass_2t                                           },
            {"g6j_HT500_g2b_2t_f1"  , passBaseline0l && pass_2t_f1                                        },
            {"g6j_HT500_g2b_2t_f2"  , passBaseline0l && pass_2t_f2                                        }, 
            {"g6j_HT500_g2b_2t_f3"  , passBaseline0l && pass_2t_f3                                        }, 
            {"g6j_HT500_g2b_2t_f4"  , passBaseline0l && pass_2t_f4                                        },

            {"g6j_HT500_g2b_2t11"   , passBaseline0l && pass_2t11                                         },
            {"g6j_HT500_g2b_2t12"   , passBaseline0l && pass_2t12                                         },
            {"g6j_HT500_g2b_2t13"   , passBaseline0l && pass_2t13                                         },
            {"g6j_HT500_g2b_2t22"   , passBaseline0l && pass_2t22                                         },
            {"g6j_HT500_g2b_2t23"   , passBaseline0l && pass_2t23                                         },
            {"g6j_HT500_g2b_2t33"   , passBaseline0l && pass_2t33                                         },
            {"g6j_HT500_g2b_2t11_f1", passBaseline0l && pass_2t11_f1                                      },
            {"g6j_HT500_g2b_2t11_f2", passBaseline0l && pass_2t11_f2                                      },
            {"g6j_HT500_g2b_2t11_f3", passBaseline0l && pass_2t11_f3                                      },
            {"g6j_HT500_g2b_2t11_f4", passBaseline0l && pass_2t11_f4                                      },
            {"g6j_HT500_g2b_2t12_f1", passBaseline0l && pass_2t12_f1                                      },
            {"g6j_HT500_g2b_2t12_f2", passBaseline0l && pass_2t12_f2                                      },
            {"g6j_HT500_g2b_2t12_f3", passBaseline0l && pass_2t12_f3                                      },
            {"g6j_HT500_g2b_2t12_f4", passBaseline0l && pass_2t12_f4                                      },
            {"g6j_HT500_g2b_2t13_f1", passBaseline0l && pass_2t13_f1                                      },
            {"g6j_HT500_g2b_2t13_f2", passBaseline0l && pass_2t13_f2                                      },
            {"g6j_HT500_g2b_2t13_f3", passBaseline0l && pass_2t13_f3                                      },
            {"g6j_HT500_g2b_2t13_f4", passBaseline0l && pass_2t13_f4                                      },
            {"g6j_HT500_g2b_2t22_f1", passBaseline0l && pass_2t22_f1                                      },
            {"g6j_HT500_g2b_2t22_f2", passBaseline0l && pass_2t22_f2                                      },
            {"g6j_HT500_g2b_2t22_f3", passBaseline0l && pass_2t22_f3                                      },
            {"g6j_HT500_g2b_2t22_f4", passBaseline0l && pass_2t22_f4                                      },
            {"g6j_HT500_g2b_2t23_f1", passBaseline0l && pass_2t23_f1                                      },
            {"g6j_HT500_g2b_2t23_f2", passBaseline0l && pass_2t23_f2                                      },
            {"g6j_HT500_g2b_2t23_f3", passBaseline0l && pass_2t23_f3                                      },
            {"g6j_HT500_g2b_2t23_f4", passBaseline0l && pass_2t23_f4                                      },
            {"g6j_HT500_g2b_2t33_f1", passBaseline0l && pass_2t33_f1                                      },
            {"g6j_HT500_g2b_2t33_f2", passBaseline0l && pass_2t33_f2                                      },
            {"g6j_HT500_g2b_2t33_f3", passBaseline0l && pass_2t33_f3                                      },
            {"g6j_HT500_g2b_2t33_f4", passBaseline0l && pass_2t33_f4                                      },
        };
        
        for(auto& kv : cut_map_0l)
        {
            if(kv.second)
            {
                my_histos["h_njets_0l_"   +kv.first]->Fill(rec_njet_pt30, eventweight);
                my_histos["h_ntops_0l_"   +kv.first]->Fill(ntops, eventweight);
                my_histos["h_nb_0l_"      +kv.first]->Fill(NBJets, eventweight);
                my_histos["h_HT_0l_"      +kv.first]->Fill(HT_trigger, eventweight);
                my_histos["h_bdt_0l_"     +kv.first]->Fill(eventshape_bdt_val, eventweight);
                my_histos["h_fisher_0l_"  +kv.first]->Fill(fisher_val, eventweight);

                my_2d_histos["h_njets_bdt_0l_"+kv.first]->Fill(rec_njet_pt30, eventshape_bdt_val, eventweight);
                my_2d_histos["h_njets_fisher_0l_"+kv.first]->Fill(rec_njet_pt30, fisher_val, eventweight);
            }
        }

        // Fill event selection efficiencies
        my_efficiencies["event_sel_total"]->Fill(true,0);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500,1);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 ,2);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && NBJets_pt45>0 ,3);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && NBJets_pt45>0 && ntops>0 ,4);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 ,5);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 && ntops>1 ,6);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 && ntops>1 && rec_njet_pt30>=8 ,7);
        
        my_efficiencies["event_sel"]->Fill(true,0);
        my_efficiencies["event_sel"]->Fill(HT_trigger>500,1);
        if(HT_trigger>500)
        {
            my_efficiencies["event_sel"]->Fill(rec_njet_pt45>=6,2);
            if (rec_njet_pt45>=6)
            {
                my_efficiencies["event_sel"]->Fill(NBJets_pt45>0,3);
                if (NBJets_pt45>0)
                {
                    my_efficiencies["event_sel"]->Fill(ntops>0,4);
                    if (ntops>0)
                    {
                        my_efficiencies["event_sel"]->Fill(NBJets_pt45>1,5);
                        if (NBJets_pt45>1)
                        {
                            my_efficiencies["event_sel"]->Fill(ntops>1,6);
                            if (ntops>1)
                            {
                                my_efficiencies["event_sel"]->Fill(rec_njet_pt30>=8,7);
                            }
                        }
                    }
                }
            }
        }
        
        // Not cuts applied here
        my_histos["h_met"     ]->Fill(MET, eventweight);
        my_histos["h_ht"      ]->Fill(HT, eventweight);
        my_histos["h_bdt"     ]->Fill(eventshape_bdt_val, eventweight);
        my_histos["h_fisher"  ]->Fill(fisher_val, eventweight);
        my_histos["h_njets"   ]->Fill(rec_njet_pt30, eventweight);
        my_histos["h_nb"      ]->Fill(NBJets, eventweight);
        my_histos["h_ntops"   ]->Fill(ntops, eventweight);        
        my_histos["h_ntops_j1"]->Fill(ntops_1jet, eventweight);
        my_histos["h_ntops_j2"]->Fill(ntops_2jet, eventweight);
        my_histos["h_ntops_j3"]->Fill(ntops_3jet, eventweight);        

    } // end of event loop

}

bool Analyze0Lep::PassTriggerGeneral(std::vector<std::string> &mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
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


bool Analyze0Lep::PassTriggerAllHad(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {
        //"HLT_PFHT1050", // 2017 trigger
        //"HLT_PFHT900"
        //"HLT_PFHT380_SixPFJet32_DoublePFBTagCSV", // 2017 trigger
        //"HLT_PFHT430_SixPFJet40_PFBTagCSV", // 2017 trigger
        "HLT_PFHT450_SixJet40_BTagCSV",
            "HLT_PFHT400_SixJet30_DoubleBTagCSV",            
            };
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

bool Analyze0Lep::PassTriggerMuon(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {"HLT_IsoMu24","HLT_IsoTkMu24_v"};
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

bool Analyze0Lep::PassTriggerElectron(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {"HLT_Ele27_WPTight_Gsf"};
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

void Analyze0Lep::WriteHistos()
{
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
