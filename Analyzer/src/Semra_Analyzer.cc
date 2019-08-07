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


/// Define histos
void Semra_Analyzer::InitHistos(const std::map<std::string, bool>& cutmap) // define variable map
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    for (const auto& cutVar : cutmap) {

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
        my_2d_histos.emplace( "h_njets_MVA_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_MVA_"+cutVar.first).c_str(), ("h_njets_MVA_"+cutVar.first).c_str(), 8, 7, 15, 50, 0, 1.0 ) );
        my_2d_histos.emplace( "h_njets_dR_bjet1_bjet2_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_dR_bjet1_bjet2_"+cutVar.first).c_str(), ("h_njets_dR_bjet1_bjet2_"+cutVar.first).c_str(), 1000, 0, 10, 20, 0, 20 ) ); // for cut optimization of dR_bjet1_bjet2 cut

    }

    // cut flow absolute numbers 
    my_absoluteCutFlow.emplace("h_cutFlow_absolute_ge2t_ge1dRbjets", std::make_shared<TH1F>("h_cutFlow_absolute_ge2t_ge1dRbjets", "h_cutFlow_absolute_ge2t_ge1dRbjets", 9,0,9));

    // Define TEfficiencies if you are doing trigger studies (for proper error bars) or cut flow charts.
    my_efficiencies.emplace("event_sel_weight_ge2t_ge1dRbjets", std::make_shared<TEfficiency>("event_sel_weight_ge2t_ge1dRbjets","event_sel_weight_ge2t_ge1dRbjets",9,0,9));

}


/// Put everything you want to do per event 
void Semra_Analyzer::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{

 
    while( tr.getNextEvent() )
    {
        const auto& eventCounter    = tr.getVar<int>("eventCounter");

        //--------------------------------------------------
        // -- Print Event Number 
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype         = tr.getVar<std::string>("runtype");     
        const auto& filetag         = tr.getVar<std::string>("filetag");
        const auto& GoodLeptons     = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& JetID           = tr.getVar<bool>("JetID");
        const auto& NGoodLeptons    = tr.getVar<int>("NGoodLeptons");
        const auto& passTriggerMC   = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodBJets_pt45 = tr.getVar<int>("NGoodBJets_pt45");
        const auto& Mbl             = tr.getVar<double>("Mbl");
        const auto& HT_trigger_pt45 = tr.getVar<double>("HT_trigger_pt45");
        const auto& NGoodJets_pt45  = tr.getVar<int>("NGoodJets_pt45");
        const auto& passMadHT       = tr.getVar<bool>("passMadHT");
        const auto& passBlindHad    = tr.getVar<bool>("passBlindHad");
        const auto& passBaseline    = tr.getVar<bool>("passBaseline1l_Good");
        const auto& MET             = tr.getVar<double>("MET");       
        const auto& Jets            = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt45   = tr.getVec<bool>("GoodJets_pt45");
        const auto& GoodBJets_pt45  = tr.getVec<bool>("GoodBJets_pt45");

        // -----------------------------------------------
        // -- SEMRA / Define Top Tag variables
        // -----------------------------------------------
        const auto& deepESM_val     = tr.getVar<double>("deepESM_val");
        const bool& pass_0l         = NGoodLeptons==0;  
        const auto& ntops           = tr.getVar<int>("ntops");
        const auto& ntops_1jet      = tr.getVar<int>("ntops_1jet");
        const auto& ntops_2jet      = tr.getVar<int>("ntops_2jet");
        const auto& ntops_3jet      = tr.getVar<int>("ntops_3jet");
        const auto& topsMass        = tr.getVec<double>("topsMass");
        const auto& topsEta         = tr.getVec<double>("topsEta");
        const auto& topsPt          = tr.getVec<double>("topsPt");
        const auto& bestTopMass     = tr.getVar<double>("bestTopMass");
        const auto& bestTopEta      = tr.getVar<double>("bestTopEta");
        const auto& bestTopPt       = tr.getVar<double>("bestTopPt");
        const auto& dR_top1_top2    = tr.getVar<double>("dR_top1_top2");
        const auto& topsLV          = tr.getVec<TLorentzVector>("topsLV");
        const auto& pass_HT500      = HT_trigger_pt45 > 500;
        const auto& pass_ge6j       = NGoodJets_pt45 >= 6;
        const auto& pass_ge1b       = NGoodBJets_pt45 >= 1;
        const auto& pass_ge2b       = NGoodBJets_pt45 >= 2;
        
        const auto& pass_ge2t       = ntops >= 2;
        const auto& pass_ge2t1j     = ntops >= 2 && ntops_3jet == 0 && ntops_2jet==0;
        const auto& pass_ge2t3j     = ntops >= 2 && ntops_1jet == 0 && ntops_2jet==0;
        const auto& pass_ge2t1j3j   = ntops >= 2 && ntops_1jet >= 1 && ntops_3jet >= 1 && ntops_2jet==0;

        const auto& pass_general    = JetID && passMadHT;

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
            const auto& Weight = tr.getVar<double>("Weight");
            const auto& lumi   = tr.getVar<double>("Lumi");
            eventweight        = lumi*Weight;
            
            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
            topPtScaleFactor     = tr.getVar<double>("topPtScaleFactor");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");
            
            weight *= eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        // --------------------------------------
        // -- Calculate DeltaR between 2 bjets
        // --------------------------------------
        double dR_bjet1_bjet2 = -1; 
        if (NGoodBJets_pt45 == 2) {
                std::vector<TLorentzVector> bjets;
                for(int ijet = 0; ijet < Jets.size(); ijet++) {
                        if(!GoodBJets_pt45[ijet]) continue;
                        bjets.push_back(Jets.at(ijet));                
                }
                dR_bjet1_bjet2 = bjets[0].DeltaR(bjets[1]);
        }
        double pass_ge1dRbjets = (dR_bjet1_bjet2 >= 1.0);

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
                {"", true                                                                                               },
                {"0l",                                     pass_general && pass_0l                                                      },
                {"0l_HT500",                               pass_general && pass_0l && pass_HT500                                        },           
                {"0l_HT500_ge2b",                          pass_general && pass_0l && pass_HT500 && pass_ge2b                           },     
                {"0l_HT500_ge2b_ge6j",                     pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j              },

                //
                {"0l_ge1dRbjets",                          pass_general && pass_0l && pass_ge1dRbjets                                  },

                // >= 2 tops
                {"0l_HT500_ge2b_ge2t",                     pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge2t                 },
                {"0l_HT500_ge2b_ge2t1j",                   pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge2t1j               },
                {"0l_HT500_ge2b_ge2t3j",                   pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge2t3j               },
                {"0l_HT500_ge2b_ge2t1j3j",                 pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge2t1j3j             },

                {"0l_HT500_ge2b_ge6j_ge2t",                pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t    },
                {"0l_HT500_ge2b_ge6j_ge2t1j",              pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t1j  },
                {"0l_HT500_ge2b_ge6j_ge2t3j",              pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t3j  },
                {"0l_HT500_ge2b_ge6j_ge2t1j3j",            pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t1j3j}, 

                // dR_bjet1_bjet2 >= 1
                {"0l_HT500_ge2b_ge2t_ge1dRbjets",          pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge2t && pass_ge1dRbjets     },
                {"0l_HT500_ge2b_ge2t1j_ge1dRbjets",        pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge2t1j && pass_ge1dRbjets   },
                {"0l_HT500_ge2b_ge2t3j_ge1dRbjets",        pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge2t3j && pass_ge1dRbjets   },
                {"0l_HT500_ge2b_ge2t1j3j_ge1dRbjets",      pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge2t1j3j && pass_ge1dRbjets },
                
                {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets",     pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t && pass_ge1dRbjets     },
                {"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets",   pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t1j && pass_ge1dRbjets   },
                {"0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets",   pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t3j && pass_ge1dRbjets   },
                {"0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjets", pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t1j3j && pass_ge1dRbjets },      

        };

        if (!inithisto) {
                InitHistos(cutmap);
                inithisto = true;
        }

        my_histos["EventCounter"]->Fill( eventCounter );

        if( !passMadHT ) continue; //Make sure not to double count DY events


        // --------------------------------
        // -- Fill the cutmap histograms
        // --------------------------------     
        for (const auto& cutVar: cutmap) {

                if (cutVar.second) {
                        my_histos["h_ntops_"+cutVar.first]->Fill( ntops, weight );
                        my_histos["h_njets_"+cutVar.first]->Fill( NGoodJets_pt45, weight );
                        my_histos["h_nbjets_"+cutVar.first]->Fill( NGoodBJets_pt45, weight );
                        my_histos["h_ht_"+cutVar.first]->Fill( HT_trigger_pt45 );
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
                        my_histos["h_dR_bjet1_bjet2_"+cutVar.first]->Fill( dR_bjet1_bjet2, weight );
                        my_histos["h_dR_top1_top2_"+cutVar.first]->Fill( dR_top1_top2, weight );
        
                        // ---------------------------------
                        // -- deltaR between top and bjet
                        // ---------------------------------
                        for (int idR = 0; idR < dR_top_bjet.size(); idR++) {
                                my_histos["h_dR_tops_bjets_"+cutVar.first]->Fill( dR_top_bjet.at(idR), weight );        
                        }
                        
                        my_2d_histos["h_njets_MVA_"+cutVar.first]->Fill( NGoodJets_pt45, deepESM_val, weight );
                        my_2d_histos["h_njets_dR_bjet1_bjet2_"+cutVar.first]->Fill( dR_bjet1_bjet2, NGoodJets_pt45, weight );
                }
        }

        // -------------------------------
        // -- Cut flow absolute numbers
        // -------------------------------
        if (true) {
            my_absoluteCutFlow["h_cutFlow_absolute_ge2t_ge1dRbjets"]->AddBinContent(1, weight);
        } if (pass_general) {
            my_absoluteCutFlow["h_cutFlow_absolute_ge2t_ge1dRbjets"]->AddBinContent(2, weight);  
        } if (pass_general && pass_0l) {
            my_absoluteCutFlow["h_cutFlow_absolute_ge2t_ge1dRbjets"]->AddBinContent(3, weight);
        } if (pass_general && pass_0l && pass_HT500) {
            my_absoluteCutFlow["h_cutFlow_absolute_ge2t_ge1dRbjets"]->AddBinContent(4, weight);
        } if (pass_general && pass_0l && pass_HT500 && pass_ge2b) {
            my_absoluteCutFlow["h_cutFlow_absolute_ge2t_ge1dRbjets"]->AddBinContent(5, weight);
        } if (pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j) {
            my_absoluteCutFlow["h_cutFlow_absolute_ge2t_ge1dRbjets"]->AddBinContent(6, weight);
        } if (pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t) {
            my_absoluteCutFlow["h_cutFlow_absolute_ge2t_ge1dRbjets"]->AddBinContent(7, weight);
        } if (pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t && pass_ge1dRbjets) {
            my_absoluteCutFlow["h_cutFlow_absolute_ge2t_ge1dRbjets"]->AddBinContent(8, weight);
        }

        // --------------------------------------------
        // -- Cut flow (event selection efficiencies)
        // --------------------------------------------
        my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->SetUseWeightedEvents();
        my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->FillWeighted(true, weight, 0);
        my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->FillWeighted(true && pass_general, weight, 1);
        my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->FillWeighted(true && pass_general && pass_0l, weight, 2);
        my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->FillWeighted(true && pass_general && pass_0l && pass_HT500, weight, 3);
        my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->FillWeighted(true && pass_general && pass_0l && pass_HT500 && pass_ge2b , weight, 4);
        my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->FillWeighted(true && pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j, weight, 5);
        my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->FillWeighted(true && pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t, weight, 6);
        my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->FillWeighted(true && pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t && pass_ge1dRbjets, weight, 7);
    
    } 
}


void Semra_Analyzer::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
 
    // Calculate the efficiency from 2 internal histograms which are called total & passed 
    TH1* total_ge2t_ge1dRbjets = my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->GetCopyTotalHisto();
    total_ge2t_ge1dRbjets->SetDirectory(outfile);
    total_ge2t_ge1dRbjets->Write();

    TH1* passed_ge2t_ge1dRbjets = my_efficiencies["event_sel_weight_ge2t_ge1dRbjets"]->GetCopyPassedHisto();
    passed_ge2t_ge1dRbjets->SetDirectory(outfile);
    passed_ge2t_ge1dRbjets->Write();  
    
    for (const auto &p : my_absoluteCutFlow) {
        p.second->Write();
    }
 
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }

    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }
    
}
