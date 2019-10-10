#define Semra_Analyzer_cxx
#include "Analyzer/Analyzer/include/Semra_Analyzer.h"
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
#include <TDirectory.h>
#include <TH1F.h>

Semra_Analyzer::Semra_Analyzer() : inithisto(false) // define inithisto variable
{
}

// -------------------
// -- Define histos
// -------------------
void Semra_Analyzer::InitHistos(const std::map<std::string, bool>& cutmap) // define variable map
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;
    
    for (const auto& cutVar : cutmap) 
    {
        my_histos.emplace( "h_ntops_"+cutVar.first, std::make_shared<TH1D> ( ("h_ntops_"+cutVar.first).c_str(), ("h_ntops_"+cutVar.first).c_str(), 10, 0, 10 ) );
        my_histos.emplace( "h_njets_"+cutVar.first, std::make_shared<TH1D> ( ("h_njets_"+cutVar.first).c_str(), ("h_njets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_nbjets_"+cutVar.first, std::make_shared<TH1D> ( ("h_nbjets_"+cutVar.first).c_str(), ("h_nbjets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_ht_"+cutVar.first, std::make_shared<TH1D> ( ("h_ht_"+cutVar.first).c_str(), ("h_ht_"+cutVar.first).c_str(), 60, 0, 3000 ) );
        my_histos.emplace( "h_met_"+cutVar.first, std::make_shared<TH1D> ( ("h_met_"+cutVar.first).c_str(), ("h_met_"+cutVar.first).c_str(), 200, 0, 2000 ) ) ;
        my_histos.emplace( "h_jetsPt_"+cutVar.first, std::make_shared<TH1D> ( ("h_jetsPt_"+cutVar.first).c_str(), ("h_jetsPt_"+cutVar.first).c_str(), 1000, 0, 2000 ) );
        my_histos.emplace( "h_jetsMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_jetsMass_"+cutVar.first).c_str(), ("h_jetsMass_"+cutVar.first).c_str(), 1000, 0, 500) );
        my_histos.emplace( "h_bjetsPt_"+cutVar.first, std::make_shared<TH1D> ( ("h_bjetsPt_"+cutVar.first).c_str(), ("h_bjetsPt_"+cutVar.first).c_str(), 1000, 0, 2000 ) );
        my_histos.emplace( "h_bjetsMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_bjetsMass_"+cutVar.first).c_str(), ("h_bjetsMass_"+cutVar.first).c_str(), 1000, 0, 500 ) );
        my_histos.emplace( "h_topsMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_topsMass_"+cutVar.first).c_str(), ("h_topsMass_"+cutVar.first).c_str(), 1000, 0, 500 ) );
        my_histos.emplace( "h_topsEta_"+cutVar.first, std::make_shared<TH1D> ( ("h_topsEta_"+cutVar.first).c_str(), ("h_topsEta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_topsPt_"+cutVar.first, std::make_shared<TH1D> ( ("h_topsPt_"+cutVar.first).c_str(), ("h_topsPt_"+cutVar.first).c_str(), 1000, 0, 2000 ) );
        my_histos.emplace( "h_bestTopMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_bestTopMass_"+cutVar.first).c_str(), ("h_bestTopMass_"+cutVar.first).c_str(), 1000, 0, 500 ) );
        my_histos.emplace( "h_bestTopEta_"+cutVar.first, std::make_shared<TH1D> ( ("h_bestTopEta_"+cutVar.first).c_str(), ("h_bestTopEta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_bestTopPt_"+cutVar.first, std::make_shared<TH1D> ( ("h_bestTopPt_"+cutVar.first).c_str(), ("h_bestTopPt_"+cutVar.first).c_str(), 1000, 0, 2000 ) );
        my_histos.emplace( "h_dR_bjet1_bjet2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_bjet1_bjet2_"+cutVar.first).c_str(), ("h_dR_bjet1_bjet2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dR_top1_top2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_top1_top2_"+cutVar.first).c_str(), ("h_dR_top1_top2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dR_tops_bjets_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_tops_bjets_"+cutVar.first).c_str(), ("h_dR_tops_bjets_"+cutVar.first).c_str(), 50, 0, 10 ) );

        // stop MT2 hemispheres  
        my_histos.emplace( "h_MT2_"+cutVar.first, std::make_shared<TH1D> ( ("h_MT2_"+cutVar.first).c_str(), ("h_MT2_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_"+cutVar.first).c_str(), ("h_stop1Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Eta_"+cutVar.first).c_str(), ("h_stop1Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop1Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Phi_"+cutVar.first).c_str(), ("h_stop1Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop1Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Pt_"+cutVar.first).c_str(), ("h_stop1Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) ); 
        my_histos.emplace( "h_stop2Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_"+cutVar.first).c_str(), ("h_stop2Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop2Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Eta_"+cutVar.first).c_str(), ("h_stop2Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop2Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Phi_"+cutVar.first).c_str(), ("h_stop2Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop2Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Pt_"+cutVar.first).c_str(), ("h_stop2Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_dR_stop1stop2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_stop1stop2_"+cutVar.first).c_str(), ("h_dR_stop1stop2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dPhi_stop1stop2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dPhi_stop1stop2_"+cutVar.first).c_str(), ("h_dPhi_stop1stop2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_stop1Mass_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_PtRank_"+cutVar.first).c_str(), ("h_stop1Mass_PtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop2Mass_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_PtRank_"+cutVar.first).c_str(), ("h_stop2Mass_PtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Mass_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_MassRank_"+cutVar.first).c_str(), ("h_stop1Mass_MassRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop2Mass_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_MassRank_"+cutVar.first).c_str(), ("h_stop2Mass_MassRank_"+cutVar.first).c_str(), 500, 0, 1500) );

        my_2d_histos.emplace( "h_Mass_stop1vsstop2_"+cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2_"+cutVar.first).c_str(), ("h_Mass_stop1vsstop2_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Eta_stop1vsstop2_"+cutVar.first, std::make_shared<TH2D>( ("h_Eta_stop1vsstop2_"+cutVar.first).c_str(), ("h_Eta_stop1vsstop2_"+cutVar.first).c_str(), 100, -6, 6, 100, -6, 6 ) );
        my_2d_histos.emplace( "h_Phi_stop1vsstop2_"+cutVar.first, std::make_shared<TH2D>( ("h_Phi_stop1vsstop2_"+cutVar.first).c_str(), ("h_Phi_stop1vsstop2_"+cutVar.first).c_str(), 80, -4, 4, 80, -4, 4 ) );
        my_2d_histos.emplace( "h_Pt_stop1vsstop2_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt_stop1vsstop2_"+cutVar.first).c_str(), ("h_Pt_stop1vsstop2_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000 ) );
        my_2d_histos.emplace( "h_NJetsVsMT2_"+cutVar.first, std::make_shared<TH2D> ( ("h_NJetsVsMT2_"+cutVar.first).c_str(), ("h_NJetsVsMT2_"+cutVar.first).c_str(), 500, 0, 1500, 20, 0, 20 ) );  
        my_2d_histos.emplace( "h_MT2Vsstop1Mass_"+cutVar.first, std::make_shared<TH2D> ( ("h_MT2Vsstop1Mass_"+cutVar.first).c_str(), ("h_MT2Vsstop1Mass_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_MT2Vsstop2Mass_"+cutVar.first, std::make_shared<TH2D> ( ("h_MT2Vsstop2Mass_"+cutVar.first).c_str(), ("h_MT2Vsstop2Mass_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_NJetsVSstop1Mass_"+cutVar.first, std::make_shared<TH2D> ( ("h_NJetsVSstop1Mass_"+cutVar.first).c_str(), ("h_NJetsVSstop1Mass_"+cutVar.first).c_str(), 500, 0, 1500, 20, 0, 20 ) );
        my_2d_histos.emplace( "h_NJetsVSstop2Mass_"+cutVar.first, std::make_shared<TH2D> ( ("h_NJetsVSstop2Mass_"+cutVar.first).c_str(), ("h_NJetsVSstop2Mass_"+cutVar.first).c_str(), 500, 0, 1500, 20, 0, 20 ) );

        my_2d_histos.emplace( "h_njets_MVA_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_MVA_"+cutVar.first).c_str(), ("h_njets_MVA_"+cutVar.first).c_str(), 8, 7, 15, 50, 0, 1.0 ) );
        my_2d_histos.emplace( "h_njets_dR_bjets_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_dR_bjets_"+cutVar.first).c_str(), ("h_njets_dR_bjets_"+cutVar.first).c_str(), 1000, 0, 10, 20, 0, 20 ) ); // for cut optimization of dR_bjets cut
    }

    // cut flow absolute numbers 
    my_histos.emplace( "h_cutFlow_absolute", std::make_shared<TH1D>("h_cutFlow_absolute", "h_cutFlow_absolute", 9,0,9));    

}

// ---------------------------------------------
// -- Put everything you want to do per event 
// ---------------------------------------------
void Semra_Analyzer::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter    = tr.getVar<int>("eventCounter");
        
        //-------------------------
        // -- Print Event Number 
        //-------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype            = tr.getVar<std::string>("runtype");     
        const auto& filetag            = tr.getVar<std::string>("filetag");
        const auto& JetID              = tr.getVar<bool>("JetID");
        const auto& NGoodLeptons       = tr.getVar<int>("NGoodLeptons");
        const auto& passTriggerMC      = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodBJets_pt45    = tr.getVar<int>("NGoodBJets_pt45");
        const auto& HT_trigger_pt45    = tr.getVar<double>("HT_trigger_pt45");
        const auto& NGoodJets_pt45     = tr.getVar<int>("NGoodJets_pt45");
        const auto& passMadHT          = tr.getVar<bool>("passMadHT");
        const auto& MET                = tr.getVar<double>("MET");       
        const auto& Jets               = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt45      = tr.getVec<bool>("GoodJets_pt45");
        const auto& GoodBJets_pt45     = tr.getVec<bool>("GoodBJets_pt45");
        const auto& dR_bjets           = tr.getVar<double>("dR_bjets");               
 
        // ------------------------------
        // -- Define Top Tag variables
        // ------------------------------
        const auto& deepESM_val        = tr.getVar<double>("deepESM_val");
        const auto& ntops              = tr.getVar<int>("ntops");
        const auto& ntops_1jet         = tr.getVar<int>("ntops_1jet");
        const auto& ntops_2jet         = tr.getVar<int>("ntops_2jet");
        const auto& ntops_3jet         = tr.getVar<int>("ntops_3jet");
        const auto& topsMass           = tr.getVec<double>("topsMass");
        const auto& topsEta            = tr.getVec<double>("topsEta");
        const auto& topsPt             = tr.getVec<double>("topsPt");
        const auto& bestTopMass        = tr.getVar<double>("bestTopMass");
        const auto& bestTopEta         = tr.getVar<double>("bestTopEta");
        const auto& bestTopPt          = tr.getVar<double>("bestTopPt");
        const auto& dR_top1_top2       = tr.getVar<double>("dR_top1_top2");
        const auto& topsLV             = tr.getVec<TLorentzVector>("topsLV");
        const auto& passBaseline0l     = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passMETFilters     = tr.getVar<bool>("passMETFilters");
        const bool pass_general        = JetID && passMETFilters && passMadHT;
        const bool pass_0l             = NGoodLeptons==0;  
        const bool pass_HT500          = HT_trigger_pt45 > 500;
        const bool pass_ge2b           = NGoodBJets_pt45 >= 2;
        const bool pass_ge6j           = NGoodJets_pt45 >= 6;
        const bool pass_ge2t           = ntops >= 2;
        const bool pass_ge2t1j         = ntops >= 2 && ntops_3jet == 0 && ntops_2jet==0;
        const bool pass_ge2t3j         = ntops >= 2 && ntops_1jet == 0 && ntops_2jet==0;
        const bool pass_ge2t1j3j       = ntops >= 2 && ntops_1jet >= 1 && ntops_3jet >= 1 && ntops_2jet==0;
        const bool pass_ge1dRbjets     = dR_bjets >= 1.0;       
    
        // -------------------------------
        // -- MT2 hemispheres variables
        // -------------------------------
        const auto& MT2                = tr.getVar<double>("MT2_0l"); 
        const auto& stop1Mass          = tr.getVar<double>("stop1Mass_0l");
        const auto& stop1Eta           = tr.getVar<double>("stop1Eta_0l");
        const auto& stop1Phi           = tr.getVar<double>("stop1Phi_0l");
        const auto& stop1Pt            = tr.getVar<double>("stop1Pt_0l");
        const auto& stop2Mass          = tr.getVar<double>("stop2Mass_0l");
        const auto& stop2Eta           = tr.getVar<double>("stop2Eta_0l");
        const auto& stop2Phi           = tr.getVar<double>("stop2Phi_0l");
        const auto& stop2Pt            = tr.getVar<double>("stop2Pt_0l");
        const auto& dR_stop1stop2      = tr.getVar<double>("dR_stop1stop2_0l");
        const auto& dPhi_stop1stop2    = tr.getVar<double>("dPhi_stop1stop2_0l");
        const auto& stop1Mass_PtRank   = tr.getVar<double>("stop1Mass_PtRank");
        const auto& stop2Mass_PtRank   = tr.getVar<double>("stop2Mass_PtRank");
        const auto& stop1Mass_MassRank = tr.getVar<double>("stop1Mass_MassRank");
        const auto& stop2Mass_MassRank = tr.getVar<double>("stop2Mass_MassRank");       
 
        // -------------------
        // -- Define weight
        // -------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double bTagScaleFactor      = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight = tr.getVar<double>("Weight");
            const auto& lumi   = tr.getVar<double>("Lumi");
            eventweight        = lumi*Weight;
        
            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");
        
            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        // ---------------------------------------------
        // -- Calculate DeltaR between tops and bjets
        // ---------------------------------------------
        std::vector<TLorentzVector> bjets;
        for(int ijet = 0; ijet < Jets.size(); ijet++) {
            if(!GoodBJets_pt45[ijet]) continue;
            bjets.push_back(Jets.at(ijet));        
        }
        std::vector<double> dR_top_bjet;
        for (double t = 0; t < topsLV.size(); t++) {
            for (double b = 0; b < bjets.size(); b++) {
                double deltaR = topsLV.at(t).DeltaR(bjets.at(b));
                dR_top_bjet.push_back(deltaR);
            }
        }
   
        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            {"",                                       true                                                                },
            {"0l",                                     pass_general && pass_0l                                             },
            {"0l_HT500",                               pass_general && pass_0l && pass_HT500                               },           
            {"0l_HT500_ge2b",                          pass_general && pass_0l && pass_HT500 && pass_ge2b                  },     
            {"0l_HT500_ge2b_ge6j",                     pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j     },
            
            // >= 2 tops
            {"0l_HT500_ge2b_ge6j_ge2t",                passBaseline0l && pass_ge2t     },
            {"0l_HT500_ge2b_ge6j_ge2t1j",              passBaseline0l && pass_ge2t1j   },
            {"0l_HT500_ge2b_ge6j_ge2t3j",              passBaseline0l && pass_ge2t3j   },
            {"0l_HT500_ge2b_ge6j_ge2t1j3j",            passBaseline0l && pass_ge2t1j3j }, 
            
            // dR_bjets >= 1
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets",     passBaseline0l && pass_ge2t && pass_ge1dRbjets     },
            {"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets   },
            {"0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets   },
            {"0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjets", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets },
            
            // stopMass1 > 200
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_stopMass1g200",    passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop1Mass > 200 },
            {"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets_stopMass1g200",  passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop1Mass > 200 },
            {"0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets_stopMass1g200",  passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop1Mass > 200 },
            {"0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjetsstopMass1g200", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop1Mass > 200 },

            // stopMass2 <= 200
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_stopMass1le200",    passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop1Mass <= 200 },
            {"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets_stopMass1le200",  passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop1Mass <= 200 },
            {"0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets_stopMass1le200",  passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop1Mass <= 200 },
            {"0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjetsstopMass1le200", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop1Mass <= 200 },

            // stopMass2 > 200
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_stopMass2g200",     passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop2Mass > 200 },
            {"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets_stopMass2g200",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop2Mass > 200 },
            {"0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets_stopMass2g200",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop2Mass > 200 },
            {"0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjetsstopMass2g200",  passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop2Mass > 200 },

            // stopMass2 <= 200
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_stopMass2le200",    passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop2Mass <= 200 },
            {"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets_stopMass2le200",  passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop2Mass <= 200 },
            {"0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets_stopMass2le200",  passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop2Mass <= 200 },
            {"0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjetsstopMass2le200", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop2Mass <= 200 },

            // NJets cuts for stop hemispheres
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_Njet7",             passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 7  },
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_Njet8",             passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 8  },
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_Njet9",             passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 9  },
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_Njet10",            passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 10 },
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_Njet11",            passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 11 },
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_Njet12",            passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 12 },
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_Njet13",            passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 13 },
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_Njet14",            passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 14 },
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets_Njet15",            passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 15 },
        };

        if (!inithisto) {
            InitHistos(cutmap);
            inithisto = true;
        }
    
        my_histos["EventCounter"]->Fill( eventCounter );
    
        // --------------------------------
        // -- Fill the cutmap histograms
        // --------------------------------     
        for (const auto& cutVar: cutmap) {        
            if (cutVar.second) {
                my_histos["h_ntops_"+cutVar.first]->Fill( ntops, weight );
                my_histos["h_njets_"+cutVar.first]->Fill( NGoodJets_pt45, weight );
                my_histos["h_nbjets_"+cutVar.first]->Fill( NGoodBJets_pt45, weight );
                my_histos["h_ht_"+cutVar.first]->Fill( HT_trigger_pt45, weight );
                my_histos["h_met_"+cutVar.first]->Fill( MET, weight );
                
                // -----------------------------
                // -- jets & bjets mass & pT
                // -----------------------------
                for(int ijet = 0; ijet < Jets.size(); ijet++) {
                    if(!GoodJets_pt45[ijet]) continue;
                    my_histos["h_jetsPt_"+cutVar.first]->Fill(Jets.at(ijet).Pt(), weight);
                    my_histos["h_jetsMass_"+cutVar.first]->Fill(Jets.at(ijet).M(), weight);                 
        
                    if(!GoodBJets_pt45[ijet]) continue;
                    const TLorentzVector& bjet = Jets.at(ijet);                     
                    my_histos["h_bjetsPt_"+cutVar.first]->Fill(bjet.Pt(), weight);
                    my_histos["h_bjetsMass_"+cutVar.first]->Fill(bjet.M(), weight);
                }
            
                // --------------------------
                // -- tops mass & eta & pT
                // --------------------------
                for (int itops = 0; itops < topsMass.size(); itops++) {
                    my_histos["h_topsMass_"+cutVar.first]->Fill( topsMass.at(itops), weight );
                }
            
                for (int itops = 0; itops < topsEta.size(); itops++) {
                    my_histos["h_topsEta_"+cutVar.first]->Fill( topsEta.at(itops), weight );
                }
            
                for (int itops = 0; itops < topsPt.size(); itops++) {
                    my_histos["h_topsPt_"+cutVar.first]->Fill( topsPt.at(itops), weight );
                }       
                     
                my_histos["h_bestTopMass_"+cutVar.first]->Fill( bestTopMass, weight );
                my_histos["h_bestTopEta_"+cutVar.first]->Fill( bestTopEta, weight );
                my_histos["h_bestTopPt_"+cutVar.first]->Fill( bestTopPt, weight );
                my_histos["h_dR_bjet1_bjet2_"+cutVar.first]->Fill( dR_bjets, weight );
                my_histos["h_dR_top1_top2_"+cutVar.first]->Fill( dR_top1_top2, weight );
        
                // ---------------------------------
                // -- deltaR between top and bjet
                // ---------------------------------
                for (int idR = 0; idR < dR_top_bjet.size(); idR++) {
                    my_histos["h_dR_tops_bjets_"+cutVar.first]->Fill( dR_top_bjet.at(idR), weight );        
                }
              
                // --------------------------
                // -- stop MT2 hemispheres 
                // --------------------------
                my_histos["h_MT2_"+cutVar.first]->Fill( MT2, weight );
                my_histos["h_stop1Mass_"+cutVar.first]->Fill( stop1Mass, weight );
                my_histos["h_stop1Eta_"+cutVar.first]->Fill( stop1Eta, weight );
                my_histos["h_stop1Phi_"+cutVar.first]->Fill( stop1Phi, weight );
                my_histos["h_stop1Pt_"+cutVar.first]->Fill( stop1Pt, weight );
                my_histos["h_stop2Mass_"+cutVar.first]->Fill( stop2Mass, weight );
                my_histos["h_stop2Eta_"+cutVar.first]->Fill( stop2Eta, weight );
                my_histos["h_stop2Phi_"+cutVar.first]->Fill( stop2Phi, weight );
                my_histos["h_stop2Pt_"+cutVar.first]->Fill( stop2Pt, weight );
                my_histos["h_dR_stop1stop2_"+cutVar.first]->Fill( dR_stop1stop2, weight );
                my_histos["h_dPhi_stop1stop2_"+cutVar.first]->Fill( dPhi_stop1stop2, weight );
                my_histos["h_stop1Mass_PtRank_"+cutVar.first]->Fill( stop1Mass_PtRank, weight );
                my_histos["h_stop2Mass_PtRank_"+cutVar.first]->Fill( stop2Mass_PtRank, weight );
                my_histos["h_stop1Mass_MassRank_"+cutVar.first]->Fill( stop1Mass_MassRank, weight );
                my_histos["h_stop2Mass_MassRank_"+cutVar.first]->Fill( stop2Mass_MassRank, weight );
                my_2d_histos["h_Mass_stop1vsstop2_"+cutVar.first]->Fill( stop1Mass, stop2Mass, weight );
                my_2d_histos["h_Mass_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{1} [GeV]");
                my_2d_histos["h_Mass_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{2} [GeV]");
                my_2d_histos["h_Eta_stop1vsstop2_"+cutVar.first]->Fill( stop1Eta, stop2Eta, weight );
                my_2d_histos["h_Eta_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("#eta_{#tildet}_{1}");
                my_2d_histos["h_Eta_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("#eta_{#tildet}_{2}");
                my_2d_histos["h_Phi_stop1vsstop2_"+cutVar.first]->Fill( stop1Phi, stop2Phi, weight );
                my_2d_histos["h_Phi_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("#phi_{#tildet}_{1}");
                my_2d_histos["h_Phi_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("#phi_{#tildet}_{2}");
                my_2d_histos["h_Pt_stop1vsstop2_"+cutVar.first]->Fill( stop1Pt, stop2Pt, weight );
                my_2d_histos["h_Pt_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("pT_{#tildet}_{1}");
                my_2d_histos["h_Pt_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("pT_{#tildet}_{2}");
                my_2d_histos["h_NJetsVsMT2_"+cutVar.first]->Fill( NGoodJets_pt45, MT2, weight );
                my_2d_histos["h_NJetsVsMT2_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_NJetsVsMT2_"+cutVar.first]->GetYaxis()->SetTitle("MT2");
                my_2d_histos["h_MT2Vsstop1Mass_"+cutVar.first]->Fill( MT2, stop1Mass, weight );
                my_2d_histos["h_MT2Vsstop1Mass_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_MT2Vsstop1Mass_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{1} [GeV]");
                my_2d_histos["h_MT2Vsstop2Mass_"+cutVar.first]->Fill( MT2, stop2Mass, weight );
                my_2d_histos["h_MT2Vsstop2Mass_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_MT2Vsstop2Mass_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{2} [GeV]");
                my_2d_histos["h_NJetsVSstop1Mass_"+cutVar.first]->Fill( NGoodJets_pt45, stop1Mass, weight );
                my_2d_histos["h_NJetsVSstop1Mass_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_NJetsVSstop1Mass_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{1}");
                my_2d_histos["h_NJetsVSstop2Mass_"+cutVar.first]->Fill( NGoodJets_pt45, stop2Mass, weight );
                my_2d_histos["h_NJetsVSstop2Mass_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_NJetsVSstop2Mass_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{2}");

                my_2d_histos["h_njets_MVA_"+cutVar.first]->Fill( NGoodJets_pt45, deepESM_val, weight );
                my_2d_histos["h_njets_MVA_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_njets_MVA_"+cutVar.first]->GetYaxis()->SetTitle("MVA");
                my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->Fill( dR_bjets, NGoodJets_pt45, weight );
                my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaR_{bjets}");
                my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->GetYaxis()->SetTitle("N_{J}");
            }
        }

        // -------------------------------
        // -- Cut flow absolute numbers
        // -------------------------------
        if(true) my_histos["h_cutFlow_absolute"]->Fill(0.5, weight);
        if(true && pass_general) my_histos["h_cutFlow_absolute"]->Fill(1.5, weight);  
        if(true && pass_general && pass_0l) my_histos["h_cutFlow_absolute"]->Fill(2.5, weight);
        if(true && pass_general && pass_0l && pass_HT500) my_histos["h_cutFlow_absolute"]->Fill(3.5, weight);
        if(true && pass_general && pass_0l && pass_HT500 && pass_ge2b) my_histos["h_cutFlow_absolute"]->Fill(4.5, weight);
        if(true && pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j) my_histos["h_cutFlow_absolute"]->Fill(5.5, weight);
        if(true && pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t) my_histos["h_cutFlow_absolute"]->Fill(6.5, weight);
        if(true && pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t && pass_ge1dRbjets) my_histos["h_cutFlow_absolute"]->Fill(7.5, weight);        
    
    } 
}

void Semra_Analyzer::WriteHistos(TFile* outfile)
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
