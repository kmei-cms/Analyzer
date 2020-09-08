#define ISRJets_Analyzer_cxx
#include "Analyzer/Analyzer/include/ISRJets_Analyzer.h"
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

ISRJets_Analyzer::ISRJets_Analyzer() : inithisto(false) // define inithisto variable
{
}

// -------------------
// -- Define histos
// -------------------
void ISRJets_Analyzer::InitHistos(const std::map<std::string, bool>& cutmap) // define variable map
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;
    
    for (const auto& cutVar : cutmap) 
    {  
        
        // -----------------------
        // -- ISR gen variables
        // -----------------------
        my_histos.emplace( "h_nGenISR_"+cutVar.first, std::make_shared<TH1D> ( ("h_nGenISR_"+cutVar.first).c_str(), ("h_nGenISR_"+cutVar.first).c_str(), 20, 0, 20 ) );
        
        my_histos.emplace( "h_GenISR_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_Mass_"+cutVar.first).c_str(), ("h_GenISR_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_GenISR_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_Pt_"+cutVar.first).c_str(), ("h_GenISR_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_GenISR_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_Phi_"+cutVar.first).c_str(), ("h_GenISR_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_GenISR_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_Eta_"+cutVar.first).c_str(), ("h_GenISR_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );

        my_histos.emplace( "h_nRecoISR_"+cutVar.first, std::make_shared<TH1D> ( ("h_nRecoISR_"+cutVar.first).c_str(), ("h_nRecoISR_"+cutVar.first).c_str(), 20, 0, 20 ) );

        // --------------------------------
        // -- ISR gen matching variables
        // --------------------------------
        // ISR gen matching deltaR
        my_histos.emplace( "h_GM_ISRmatching_allDR_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_allDR_"+cutVar.first).c_str(), ("h_GM_ISRmatching_allDR_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_bestDR_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_bestDR_"+cutVar.first).c_str(), ("h_GM_ISRmatching_bestDR_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first).c_str(), ("h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first).c_str(), ("h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first).c_str(), 1000, 0, 10 ) );        
        // ISR gen matching pt ratio
        my_histos.emplace( "h_GM_ISRmatching_allPtRatio_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_allPtRatio_"+cutVar.first).c_str(), ("h_GM_ISRmatching_allPtRatio_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_bestPtRatio_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_bestPtRatio_"+cutVar.first).c_str(), ("h_GM_ISRmatching_bestPtRatio_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first).c_str(), ("h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first).c_str(), ("h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first).c_str(), 1000, 0, 10 ) ); 

        // -----------------------------------
        // -- ISR and NON ISR Jet variables
        // -----------------------------------
        // ISR Jets - matches cutting on dr and ptr
        my_histos.emplace( "h_nISRJets_drPTR_"+cutVar.first, std::make_shared<TH1D> ( ("h_nISRJets_drPTR_"+cutVar.first).c_str(), ("h_nISRJets_drPTR_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_ISRJets_drPTR_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_drPTR_Mass_"+cutVar.first).c_str(), ("h_ISRJets_drPTR_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_ISRJets_drPTR_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_drPTR_Pt_"+cutVar.first).c_str(), ("h_ISRJets_drPTR_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_ISRJets_drPTR_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_drPTR_Phi_"+cutVar.first).c_str(), ("h_ISRJets_drPTR_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );       
        my_histos.emplace( "h_ISRJets_drPTR_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_drPTR_Eta_"+cutVar.first).c_str(), ("h_ISRJets_drPTR_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) ); 
        // ISR Jets - matches cutting on dr
        my_histos.emplace( "h_nISRJets_dr_"+cutVar.first, std::make_shared<TH1D> ( ("h_nISRJets_dr_"+cutVar.first).c_str(), ("h_nISRJets_dr_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_ISRJets_dr_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_dr_Mass_"+cutVar.first).c_str(), ("h_ISRJets_dr_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_ISRJets_dr_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_dr_Pt_"+cutVar.first).c_str(), ("h_ISRJets_dr_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_ISRJets_dr_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_dr_Phi_"+cutVar.first).c_str(), ("h_ISRJets_dr_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_ISRJets_dr_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_dr_Eta_"+cutVar.first).c_str(), ("h_ISRJets_dr_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        // Jets_ISRMask - From Ntuple to check
        my_histos.emplace( "h_nJets_ISRMask_"+cutVar.first, std::make_shared<TH1D> ( ("h_nJets_ISRMask_"+cutVar.first).c_str(), ("h_nJets_ISRMask_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_Jets_ISRMask_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Jets_ISRMask_Mass_"+cutVar.first).c_str(), ("h_Jets_ISRMask_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_Jets_ISRMask_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_Jets_ISRMask_Pt_"+cutVar.first).c_str(), ("h_Jets_ISRMask_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_Jets_ISRMask_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_Jets_ISRMask_Phi_"+cutVar.first).c_str(), ("h_Jets_ISRMask_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_Jets_ISRMask_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_Jets_ISRMask_Eta_"+cutVar.first).c_str(), ("h_Jets_ISRMask_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        // NON ISR Jets
        my_histos.emplace( "h_nNonISRJets_"+cutVar.first, std::make_shared<TH1D> ( ("h_nNonISRJets_"+cutVar.first).c_str(), ("h_nNonISRJets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_NonISRJets_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_Mass_"+cutVar.first).c_str(), ("h_NonISRJets_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_NonISRJets_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_Pt_"+cutVar.first).c_str(), ("h_NonISRJets_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_NonISRJets_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_Phi_"+cutVar.first).c_str(), ("h_NonISRJetsPhi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_NonISRJets_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_Eta_"+cutVar.first).c_str(), ("h_NonISRJets_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        // DeltaR between ISR and Non ISR jets
        my_histos.emplace( "h_dR_ISRJetsDrPtr_NonISRJets_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_ISRJetsDrPtr_NonISRJets_"+cutVar.first).c_str(), ("h_dR_ISRJetsDrPtr_NonISRJets_"+cutVar.first).c_str(), 50, 0, 10 ) );

        // -------------------
        // -- 2D histograms
        // -------------------
        my_2d_histos.emplace( "h_GM_all_PtRvsDR_"+cutVar.first, std::make_shared<TH2D>( ("h_GM_all_PtRvsDR_"+cutVar.first).c_str(), ("h_GM_all_PtRvsDR_"+cutVar.first).c_str(), 1000, 0, 10, 1000, 0, 10) );
        my_2d_histos.emplace( "h_GM_best_PtRvsDR_"+cutVar.first, std::make_shared<TH2D>( ("h_GM_best_PtRvsDR_"+cutVar.first).c_str(), ("h_GM_best_PtRvsDR_"+cutVar.first).c_str(), 1000, 0, 10, 1000, 0, 10) );
        my_2d_histos.emplace( "h_GM_CutOnDR_PtRvsDR_"+cutVar.first, std::make_shared<TH2D>( ("h_GM_CutOnDR_PtRvsDR_"+cutVar.first).c_str(), ("h_GM_CutOnDR_PtRvsDR_"+cutVar.first).c_str(), 1000, 0, 10, 1000, 0, 10) );
        my_2d_histos.emplace( "h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first, std::make_shared<TH2D>( ("h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first).c_str(), ("h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first).c_str(), 1000, 0, 10, 1000, 0, 10) );

    }
}

// ---------------------------------------------
// -- Put everything you want to do per event 
// ---------------------------------------------
void ISRJets_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter    = tr.getVar<int>("eventCounter");
        
        //-------------------------
        // -- Print Event Number 
        //-------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype         = tr.getVar<std::string>("runtype");     
        const auto& Jets            = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt20   = tr.getVec<bool>("GoodJets_pt20");
        const auto& GenParticles    = tr.getVec<TLorentzVector>("GenParticles");
        const auto& Jets_ISRMask    = tr.getVec<bool>("Jets_ISRMask");
        const auto& passBaseline0l  = tr.getVar<bool>("passBaseline0l_Good");
        const auto& ntops           = tr.getVar<int>("ntops");
        const auto& ntops_1jet      = tr.getVar<int>("ntops_1jet"); // merged
        const auto& ntops_2jet      = tr.getVar<int>("ntops_2jet");
        const auto& ntops_3jet      = tr.getVar<int>("ntops_3jet"); // resolved 
        const auto& dR_bjets        = tr.getVar<double>("dR_bjets");
        const bool  pass_ge2t       = ntops >= 2;
        const bool  pass_ge2t1j     = ntops >= 2 && ntops_3jet == 0 && ntops_2jet==0;
        const bool  pass_ge2t3j     = ntops >= 2 && ntops_1jet == 0 && ntops_2jet==0;
        const bool  pass_ge2t1j3j   = ntops >= 2 && ntops_1jet >= 1 && ntops_3jet >= 1 && ntops_2jet==0;
        const bool  pass_ge1dRbjets = dR_bjets >= 1.0;       

        // ----------------------
        // -- ISR gen matching
        // ----------------------
        const auto& GenISR                                  = tr.getVec<bool>("GenISR");
        const auto& nGenISR                                 = tr.getVar<int>("nGenISR");
        const auto& nRecoISR                                = tr.getVar<int>("nRecoISR");
        const auto& ISRmatched_dr_ptr                       = tr.getVec<bool>("ISRmatched_dr_ptr");
        const auto& ISRmatched_dr                           = tr.getVec<bool>("ISRmatched_dr");
        const auto& GM_ISRmatching_allDR                    = tr.getVec<double>("GM_ISRmatching_allDR");
        const auto& GM_ISRmatching_bestDR                   = tr.getVec<double>("GM_ISRmatching_bestDR");
        const auto& GM_ISRmatching_justCutOnDR_DR           = tr.getVec<double>("GM_ISRmatching_justCutOnDR_DR");
        const auto& GM_ISRmatching_justCutOnPtRatio_DR      = tr.getVec<double>("GM_ISRmatching_justCutOnPtRatio_DR");
        const auto& GM_ISRmatching_allPtRatio               = tr.getVec<double>("GM_ISRmatching_allPtRatio");
        const auto& GM_ISRmatching_bestPtRatio              = tr.getVec<double>("GM_ISRmatching_bestPtRatio");
        const auto& GM_ISRmatching_justCutOnDR_PtRatio      = tr.getVec<double>("GM_ISRmatching_justCutOnDR_PtRatio");
        const auto& GM_ISRmatching_justCutOnPtRatio_PtRatio = tr.getVec<double>("GM_ISRmatching_justCutOnPtRatio_PtRatio");

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
            const auto& Weight      = tr.getVar<double>("Weight");
            const auto& lumi        = tr.getVar<double>("Lumi");
            eventweight             = lumi*Weight;
        
            bTagScaleFactor         = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor    = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor           = tr.getVar<double>("puWeightCorr");
        
            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        // --------------------------------------
        // DeltaR between ISR and Non ISR jets
        // --------------------------------------
        std::vector<TLorentzVector> ISRJets_dr_ptr;
        std::vector<TLorentzVector> NonISRJets;
        for(unsigned int j = 0; j < Jets.size(); j++) 
        {
            if(!GoodJets_pt20[j]) continue;
            if (ISRmatched_dr_ptr[j])  
            { 
                ISRJets_dr_ptr.push_back(Jets.at(j));
            } else 
            {
                NonISRJets.push_back(Jets.at(j));
            }
        }
        
        std::vector<double> dR_ISRJetsDrPtr_NonISRJets;
        for (unsigned int i = 0; i < ISRJets_dr_ptr.size(); i++) 
        {
            for (unsigned int n = 0; n < NonISRJets.size(); n++) 
            {
                double deltaR = ISRJets_dr_ptr.at(i).DeltaR(NonISRJets.at(n));
                dR_ISRJetsDrPtr_NonISRJets.push_back(deltaR);
            }
        }


        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            {"",                                       true                                               },
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets",     passBaseline0l && pass_ge2t && pass_ge1dRbjets     },
            {"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets   },
            {"0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets   },
            {"0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjets", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets },
    
        };

        if (!inithisto) 
        {
            InitHistos(cutmap);
            inithisto = true;
        }

        my_histos["EventCounter"]->Fill( eventCounter );

        // --------------------------------
        // -- Fill the cutmap histograms
        // --------------------------------     
        for (const auto& cutVar: cutmap) 
        { 
               
            if (cutVar.second) 
            {

                // -----------------------
                // -- ISR gen variables
                // -----------------------
                my_histos["h_nGenISR_"+cutVar.first]->Fill( nGenISR, weight );
                my_histos["h_nGenISR_"+cutVar.first]->GetXaxis()->SetTitle("nGenISR");
                my_histos["h_nGenISR_"+cutVar.first]->GetYaxis()->SetTitle("Events");                

                for (unsigned int g = 0; g < GenParticles.size(); g++)
                {
                    if (!GoodJets_pt20[g]) continue;
                    if (GenISR[g])
                    {  
                        my_histos["h_GenISR_Mass_"+cutVar.first]->Fill( GenParticles[g].M(), weight );
                        my_histos["h_GenISR_Mass_"+cutVar.first]->GetXaxis()->SetTitle("GenISR Mass");
                        my_histos["h_GenISR_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_GenISR_Pt_"+cutVar.first]->Fill( GenParticles[g].Pt(), weight );
                        my_histos["h_GenISR_Pt_"+cutVar.first]->GetXaxis()->SetTitle("GenISR Pt");
                        my_histos["h_GenISR_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_GenISR_Phi_"+cutVar.first]->Fill( GenParticles[g].Phi(), weight );
                        my_histos["h_GenISR_Phi_"+cutVar.first]->GetXaxis()->SetTitle("GenISR #phi");
                        my_histos["h_GenISR_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_GenISR_Eta_"+cutVar.first]->Fill( GenParticles[g].Eta(), weight );
                        my_histos["h_GenISR_Eta_"+cutVar.first]->GetXaxis()->SetTitle("GenISR #eta");
                        my_histos["h_GenISR_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    }
                }
                my_histos["h_nRecoISR_"+cutVar.first]->Fill( nRecoISR, weight );
                my_histos["h_nRecoISR_"+cutVar.first]->GetXaxis()->SetTitle("nRecoISR");
                my_histos["h_nRecoISR_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                // -------------------------------
                // -- ISR gen matching variables
                // -------------------------------
                // --------------------
                // all possible matches
                // --------------------
                for (double dr = 0.0; dr < GM_ISRmatching_allDR.size(); dr++)
                {
                    my_histos["h_GM_ISRmatching_allDR_"+cutVar.first]->Fill( GM_ISRmatching_allDR.at(dr), weight );
                    my_histos["h_GM_ISRmatching_allDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching all #DeltaR");
                    my_histos["h_GM_ISRmatching_allDR_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                    // 2D histogram - dr vs pt ratio
                    my_2d_histos["h_GM_all_PtRvsDR_"+cutVar.first]->Fill(GM_ISRmatching_allPtRatio.at(dr), GM_ISRmatching_allDR.at(dr), weight);
                    my_2d_histos["h_GM_all_PtRvsDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching all (1 - (Pt^{Reco}/Pt^{Gen}))");
                    my_2d_histos["h_GM_all_PtRvsDR_"+cutVar.first]->GetYaxis()->SetTitle("ISRmatching all #DeltaR");
                }

                for (double pt = 0.0; pt < GM_ISRmatching_allPtRatio.size(); pt++)
                {
                    my_histos["h_GM_ISRmatching_allPtRatio_"+cutVar.first]->Fill( GM_ISRmatching_allPtRatio.at(pt), weight );
                    my_histos["h_GM_ISRmatching_allPtRatio_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching all (1 - (Pt^{Reco}/Pt^{Gen}))");
                    my_histos["h_GM_ISRmatching_allPtRatio_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }

                // ------------
                // best matches
                // ------------
                for (double dr = 0.0; dr < GM_ISRmatching_bestDR.size(); dr++)
                {
                    my_histos["h_GM_ISRmatching_bestDR_"+cutVar.first]->Fill( GM_ISRmatching_bestDR.at(dr), weight );
                    my_histos["h_GM_ISRmatching_bestDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching best #DeltaR");
                    my_histos["h_GM_ISRmatching_bestDR_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                    // 2D histogram - dr vs pt ratio
                    my_2d_histos["h_GM_best_PtRvsDR_"+cutVar.first]->Fill(GM_ISRmatching_bestPtRatio.at(dr), GM_ISRmatching_bestDR.at(dr), weight);
                    my_2d_histos["h_GM_best_PtRvsDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching best (1 - (Pt^{Reco}/Pt^{Gen}))");
                    my_2d_histos["h_GM_best_PtRvsDR_"+cutVar.first]->GetYaxis()->SetTitle("ISRmatching best #DeltaR");
                }

                for (double pt = 0.0; pt < GM_ISRmatching_bestPtRatio.size(); pt++)
                {
                    my_histos["h_GM_ISRmatching_bestPtRatio_"+cutVar.first]->Fill( GM_ISRmatching_bestPtRatio.at(pt), weight );
                    my_histos["h_GM_ISRmatching_bestPtRatio_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching best (1 - (Pt^{Reco}/Pt^{Gen}))");
                    my_histos["h_GM_ISRmatching_bestPtRatio_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }

                // ------------------------------
                // matches with cutting on deltaR
                // ------------------------------
                for (double dr = 0.0; dr < GM_ISRmatching_justCutOnDR_DR.size(); dr++)
                {
                    my_histos["h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first]->Fill( GM_ISRmatching_justCutOnDR_DR.at(dr), weight );
                    my_histos["h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnDR #DeltaR");
                    my_histos["h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                    // 2D histogram - dr vs pt ratio
                    my_2d_histos["h_GM_CutOnDR_PtRvsDR_"+cutVar.first]->Fill(GM_ISRmatching_justCutOnDR_PtRatio.at(dr), GM_ISRmatching_justCutOnDR_DR.at(dr), weight);
                    my_2d_histos["h_GM_CutOnDR_PtRvsDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnDR (1 - (Pt^{Reco}/Pt^{Gen}))");
                    my_2d_histos["h_GM_CutOnDR_PtRvsDR_"+cutVar.first]->GetYaxis()->SetTitle("ISRmatching CutOnDR #DeltaR");   
                }

                for (double pt = 0.0; pt < GM_ISRmatching_justCutOnDR_PtRatio.size(); pt++)
                {
                    my_histos["h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first]->Fill( GM_ISRmatching_justCutOnDR_PtRatio.at(pt), weight );
                    my_histos["h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnDR (1 - (Pt^{Reco}/Pt^{Gen}))");
                    my_histos["h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }

                // --------------------------------
                // matches with cutting on pt ratio
                // --------------------------------
                for (double dr = 0.0; dr < GM_ISRmatching_justCutOnPtRatio_DR.size(); dr++)
                {
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first]->Fill( GM_ISRmatching_justCutOnPtRatio_DR.at(dr), weight );
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnPtRatio #DeltaR");
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first]->GetYaxis()->SetTitle("Events");                
                    // 2D histogram - dr vs pt ratio
                    my_2d_histos["h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first]->Fill(GM_ISRmatching_justCutOnPtRatio_PtRatio.at(dr), GM_ISRmatching_justCutOnPtRatio_DR.at(dr), weight);
                    my_2d_histos["h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnPtRatio (1 - (Pt^{Reco}/Pt^{Gen}))");
                    my_2d_histos["h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first]->GetYaxis()->SetTitle("ISRmatching CutOnPtRatio #DeltaR");
                }

                for (double pt = 0.0; pt < GM_ISRmatching_justCutOnPtRatio_PtRatio.size(); pt++)
                {
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first]->Fill( GM_ISRmatching_justCutOnPtRatio_PtRatio.at(pt), weight );
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnPtRatio (1 - (Pt^{Reco}/Pt^{Gen}))");
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }            

                // -----------------------------------
                // -- ISR and NON ISR Jet variables
                // ------------------------------------ 
                int nISRJets_dr_ptr = 0;
                int nISRJets_dr     = 0;
                int nJets_ISRMask   = 0;
                int nNonISRJets     = 0;
                for (unsigned int j = 0; j < Jets.size(); j++)
                {
                    if (!GoodJets_pt20[j]) continue;
                    
                    // ISR Jets - matches cutting on dr and ptr
                    if (ISRmatched_dr_ptr[j])
                    {
                        my_histos["h_ISRJets_drPTR_Mass_"+cutVar.first]->Fill( Jets[j].M(), weight );
                        my_histos["h_ISRJets_drPTR_Mass_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_histos["h_ISRJets_drPTR_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_drPTR_Pt_"+cutVar.first]->Fill( Jets[j].Pt(), weight );
                        my_histos["h_ISRJets_drPTR_Pt_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_histos["h_ISRJets_drPTR_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_drPTR_Phi_"+cutVar.first]->Fill( Jets[j].Phi(), weight );
                        my_histos["h_ISRJets_drPTR_Phi_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets #phi");
                        my_histos["h_ISRJets_drPTR_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_drPTR_Eta_"+cutVar.first]->Fill( Jets[j].Eta(), weight );
                        my_histos["h_ISRJets_drPTR_Eta_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets #eta");
                        my_histos["h_ISRJets_drPTR_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        nISRJets_dr_ptr++;
                    }
            
                    // ISR Jets - matches cutting on dr 
                    if(ISRmatched_dr[j])
                    {
                        my_histos["h_ISRJets_dr_Mass_"+cutVar.first]->Fill( Jets[j].M(), weight );
                        my_histos["h_ISRJets_dr_Mass_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_histos["h_ISRJets_dr_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_dr_Pt_"+cutVar.first]->Fill( Jets[j].Pt(), weight );
                        my_histos["h_ISRJets_dr_Pt_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_histos["h_ISRJets_dr_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_dr_Phi_"+cutVar.first]->Fill( Jets[j].Phi(), weight );
                        my_histos["h_ISRJets_dr_Phi_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets #phi");
                        my_histos["h_ISRJets_dr_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_dr_Eta_"+cutVar.first]->Fill( Jets[j].Eta(), weight );
                        my_histos["h_ISRJets_dr_Eta_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets #eta");
                        my_histos["h_ISRJets_dr_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");        
                        nISRJets_dr++;
                    }        

                    // From Ntuple to check
                    if (Jets_ISRMask[j])
                    {                
                        my_histos["h_Jets_ISRMask_Mass_"+cutVar.first]->Fill( Jets[j].M(), weight );
                        my_histos["h_Jets_ISRMask_Mass_"+cutVar.first]->GetXaxis()->SetTitle("Jets_ISRMask Mass");
                        my_histos["h_Jets_ISRMask_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_Jets_ISRMask_Pt_"+cutVar.first]->Fill( Jets[j].Pt(), weight );
                        my_histos["h_Jets_ISRMask_Pt_"+cutVar.first]->GetXaxis()->SetTitle("Jets_ISRMask Pt");
                        my_histos["h_Jets_ISRMask_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_Jets_ISRMask_Phi_"+cutVar.first]->Fill( Jets[j].Phi(), weight );
                        my_histos["h_Jets_ISRMask_Phi_"+cutVar.first]->GetXaxis()->SetTitle("Jets_ISRMask #phi");
                        my_histos["h_Jets_ISRMask_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_Jets_ISRMask_Eta_"+cutVar.first]->Fill( Jets[j].Eta(), weight );
                        my_histos["h_Jets_ISRMask_Eta_"+cutVar.first]->GetXaxis()->SetTitle("Jets_ISRMask #eta");
                        my_histos["h_Jets_ISRMask_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        nJets_ISRMask++;

                    // Non ISR Jets
                    } else
                    {
                        my_histos["h_NonISRJets_Mass_"+cutVar.first]->Fill( Jets[j].M(), weight );
                        my_histos["h_NonISRJets_Mass_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_histos["h_NonISRJets_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_Pt_"+cutVar.first]->Fill( Jets[j].Pt(), weight );
                        my_histos["h_NonISRJets_Pt_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_histos["h_NonISRJets_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_Phi_"+cutVar.first]->Fill( Jets[j].Phi(), weight );
                        my_histos["h_NonISRJets_Phi_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets #phi");
                        my_histos["h_NonISRJets_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_Eta_"+cutVar.first]->Fill( Jets[j].Eta(), weight );
                        my_histos["h_NonISRJets_Eta_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets #eta");
                        nNonISRJets++;
                    }
                }
                my_histos["h_nISRJets_drPTR_"+cutVar.first]->Fill( nISRJets_dr_ptr, weight );
                my_histos["h_nISRJets_drPTR_"+cutVar.first]->GetXaxis()->SetTitle("nISRJets");
                my_histos["h_nISRJets_drPTR_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                my_histos["h_nISRJets_dr_"+cutVar.first]->Fill( nISRJets_dr, weight );
                my_histos["h_nISRJets_dr_"+cutVar.first]->GetXaxis()->SetTitle("nISRJets");
                my_histos["h_nISRJets_dr_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                
                my_histos["h_nJets_ISRMask_"+cutVar.first]->Fill( nJets_ISRMask, weight );
                my_histos["h_nJets_ISRMask_"+cutVar.first]->GetXaxis()->SetTitle("nJets_ISRMask");
                my_histos["h_nJets_ISRMask_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                my_histos["h_nNonISRJets_"+cutVar.first]->Fill( nNonISRJets, weight );
                my_histos["h_nNonISRJets_"+cutVar.first]->GetXaxis()->SetTitle("nNonISRJets");
                my_histos["h_nNonISRJets_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                // -----------------------------------
                // DeltaR between ISR and Non ISR jets
                // ----------------------------------- 
                for (unsigned int i = 0; i < dR_ISRJetsDrPtr_NonISRJets.size(); i++) 
                {
                    my_histos["h_dR_ISRJetsDrPtr_NonISRJets_"+cutVar.first]->Fill( dR_ISRJetsDrPtr_NonISRJets.at(i), weight ); 
                }
                
            }
        }
    }    
}

void ISRJets_Analyzer::WriteHistos(TFile* outfile)
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
