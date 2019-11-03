#define AnalyzeEENoise_cxx
#include "Analyzer/Analyzer/include/AnalyzeEENoise.h"
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

AnalyzeEENoise::AnalyzeEENoise()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeEENoise::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    my_histos.emplace( "h_fwm2_top6_2018A",std::make_shared<TH1D>( "h_fwm2_top6_2018A", "h_fwm2_top6_2018A", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm3_top6_2018A",std::make_shared<TH1D>( "h_fwm3_top6_2018A", "h_fwm3_top6_2018A", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm4_top6_2018A",std::make_shared<TH1D>( "h_fwm4_top6_2018A", "h_fwm4_top6_2018A", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm5_top6_2018A",std::make_shared<TH1D>( "h_fwm5_top6_2018A", "h_fwm5_top6_2018A", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev0_top6_2018A",std::make_shared<TH1D>( "h_jmt_ev0_top6_2018A", "h_jmt_ev0_top6_2018A", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev1_top6_2018A",std::make_shared<TH1D>( "h_jmt_ev1_top6_2018A", "h_jmt_ev1_top6_2018A", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev2_top6_2018A",std::make_shared<TH1D>( "h_jmt_ev2_top6_2018A", "h_jmt_ev2_top6_2018A", 1440, 0, 1) ) ;
    my_histos.emplace( "h_beta_z_2018A", std::make_shared<TH1D>( "h_beta_z_2018A", "h_beta_z_2018A", 1440, -1, 1) ) ;
    my_histos.emplace( "h_beta_z_pt20_2018A", std::make_shared<TH1D>( "h_beta_z_pt20_2018A", "h_beta_z_pt20_2018A", 1440, -1, 1) ) ;

    my_2d_histos.emplace( "h_beta_z_nvtx_2018A", std::make_shared<TH2D>( "h_beta_z_nvtx_2018A", "h_beta_z_nvtx_2018A", 1440, -1, 1, 121, -0.5, 120.5) ) ;
    my_2d_histos.emplace( "h_beta_z_pt20_nvtx_2018A", std::make_shared<TH2D>( "h_beta_z_pt20_nvtx_2018A", "h_beta_z_pt20_nvtx_2018A", 1440, -1, 1, 121, -0.5, 120.5) ) ;

    my_histos.emplace( "h_deepESM_2018A", std::make_shared<TH1D>( "h_deepESM_2018A", "h_deepESM_val_2018A", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jet_pt_2018A", std::make_shared<TH1D>( "h_jet_pt_2018A", "h_jet_pt_2018A", 1440, 0, 1440) );
    my_histos.emplace( "h_jet_eta_2018A", std::make_shared<TH1D>( "h_jet_eta_2018A", "h_jet_eta_2018A", 1440, -5, 5) );
    my_2d_histos.emplace( "h_jet_etaphi_2018A", std::make_shared<TH2D>( "h_jet_etaphi_2018A", "h_jet_etaphi_2018A", 1440, -5, 5, 1440, -3.14, 3.14) );
    my_2d_histos.emplace( "h_jet_pteta_2018A", std::make_shared<TH2D>( "h_jet_pteta_2018A", "h_jet_pteta_2018A", 1440, 0, 1440, 1440, -5, 5) );
    my_2d_histos.emplace( "h_jet_ptem_2018A", std::make_shared<TH2D>( "h_jet_ptem_2018A", "h_jet_ptem_2018A", 1440, 0, 1440, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet_pthad_2018A", std::make_shared<TH2D>( "h_jet_pthad_2018A", "h_jet_pthad_2018A", 1440, 0, 1440, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet_etaem_2018A", std::make_shared<TH2D>( "h_jet_etaem_2018A", "h_jet_etaem_2018A", 1440, -5, 5, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet_etahad_2018A", std::make_shared<TH2D>( "h_jet_etahad_2018A", "h_jet_etahad_2018A", 1440, -5, 5, 1440, 0, 1) );
    my_histos.emplace( "h_jet30_pt_2018A", std::make_shared<TH1D>( "h_jet30_pt_2018A", "h_jet30_pt_2018A", 1440, 0, 1440) );
    my_histos.emplace( "h_jet30_eta_2018A", std::make_shared<TH1D>( "h_jet30_eta_2018A", "h_jet30_eta_2018A", 1440, -5, 5) );
    my_2d_histos.emplace( "h_jet30_etaphi_2018A", std::make_shared<TH2D>( "h_jet30_etaphi_2018A", "h_jet30_etaphi_2018A", 1440, -5, 5, 1440, -3.14, 3.14) );
    my_2d_histos.emplace( "h_jet30_pteta_2018A", std::make_shared<TH2D>( "h_jet30_pteta_2018A", "h_jet30_pteta_2018A", 1440, 0, 1440, 1440, -5, 5) );
    my_2d_histos.emplace( "h_jet30_ptem_2018A", std::make_shared<TH2D>( "h_jet30_ptem_2018A", "h_jet30_ptem_2018A", 1440, 0, 1440, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet30_pthad_2018A", std::make_shared<TH2D>( "h_jet30_pthad_2018A", "h_jet30_pthad_2018A", 1440, 0, 1440, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet30_etaem_2018A", std::make_shared<TH2D>( "h_jet30_etaem_2018A", "h_jet30_etaem_2018A", 1440, -5, 5, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet30_etahad_2018A", std::make_shared<TH2D>( "h_jet30_etahad_2018A", "h_jet30_etahad_2018A", 1440, -5, 5, 1440, 0, 1) );

    my_histos.emplace( "h_fwm2_top6_2018D",std::make_shared<TH1D>( "h_fwm2_top6_2018D", "h_fwm2_top6_2018D", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm3_top6_2018D",std::make_shared<TH1D>( "h_fwm3_top6_2018D", "h_fwm3_top6_2018D", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm4_top6_2018D",std::make_shared<TH1D>( "h_fwm4_top6_2018D", "h_fwm4_top6_2018D", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm5_top6_2018D",std::make_shared<TH1D>( "h_fwm5_top6_2018D", "h_fwm5_top6_2018D", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev0_top6_2018D",std::make_shared<TH1D>( "h_jmt_ev0_top6_2018D", "h_jmt_ev0_top6_2018D", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev1_top6_2018D",std::make_shared<TH1D>( "h_jmt_ev1_top6_2018D", "h_jmt_ev1_top6_2018D", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev2_top6_2018D",std::make_shared<TH1D>( "h_jmt_ev2_top6_2018D", "h_jmt_ev2_top6_2018D", 1440, 0, 1) ) ;
    my_histos.emplace( "h_beta_z_2018D", std::make_shared<TH1D>( "h_beta_z_2018D", "h_beta_z_2018D", 1440, -1, 1) ) ;
    my_histos.emplace( "h_beta_z_pt20_2018D", std::make_shared<TH1D>( "h_beta_z_pt20_2018D", "h_beta_z_pt20_2018D", 1440, -1, 1) ) ;

    my_2d_histos.emplace( "h_beta_z_nvtx_2018D", std::make_shared<TH2D>( "h_beta_z_nvtx_2018D", "h_beta_z_nvtx_2018D", 1440, -1, 1, 121, -0.5, 120.5) ) ;
    my_2d_histos.emplace( "h_beta_z_pt20_nvtx_2018D", std::make_shared<TH2D>( "h_beta_z_pt20_nvtx_2018D", "h_beta_z_pt20_nvtx_2018D", 1440, -1, 1, 121, -0.5, 120.5) ) ;

    my_histos.emplace( "h_deepESM_2018D", std::make_shared<TH1D>( "h_deepESM_2018D", "h_deepESM_2018D", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jet_pt_2018D", std::make_shared<TH1D>( "h_jet_pt_2018D", "h_jet_pt_2018D", 1440, 0, 1440) );
    my_histos.emplace( "h_jet_eta_2018D", std::make_shared<TH1D>( "h_jet_eta_2018D", "h_jet_eta_2018D", 1440, -5, 5) );
    my_2d_histos.emplace( "h_jet_etaphi_2018D", std::make_shared<TH2D>( "h_jet_etaphi_2018D", "h_jet_etaphi_2018D", 1440, -5, 5, 1440, -3.14, 3.14) );
    my_2d_histos.emplace( "h_jet_pteta_2018D", std::make_shared<TH2D>( "h_jet_pteta_2018D", "h_jet_pteta_2018D", 1440, 0, 1440, 1440, -5, 5) );
    my_2d_histos.emplace( "h_jet_ptem_2018D", std::make_shared<TH2D>( "h_jet_ptem_2018D", "h_jet_ptem_2018D", 1440, 0, 1440, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet_pthad_2018D", std::make_shared<TH2D>( "h_jet_pthad_2018D", "h_jet_pthad_2018D", 1440, 0, 1440, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet_etaem_2018D", std::make_shared<TH2D>( "h_jet_etaem_2018D", "h_jet_etaem_2018D", 1440, -5, 5, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet_etahad_2018D", std::make_shared<TH2D>( "h_jet_etahad_2018D", "h_jet_etahad_2018D", 1440, -5, 5, 1440, 0, 1) );
    my_histos.emplace( "h_jet30_pt_2018D", std::make_shared<TH1D>( "h_jet30_pt_2018D", "h_jet30_pt_2018D", 1440, 0, 1440) );
    my_histos.emplace( "h_jet30_eta_2018D", std::make_shared<TH1D>( "h_jet30_eta_2018D", "h_jet30_eta_2018D", 1440, -5, 5) );
    my_2d_histos.emplace( "h_jet30_etaphi_2018D", std::make_shared<TH2D>( "h_jet30_etaphi_2018D", "h_jet30_etaphi_2018D", 1440, -5, 5, 1440, -3.14, 3.14) );
    my_2d_histos.emplace( "h_jet30_pteta_2018D", std::make_shared<TH2D>( "h_jet30_pteta_2018D", "h_jet30_pteta_2018D", 1440, 0, 1440, 1440, -5, 5) );
    my_2d_histos.emplace( "h_jet30_ptem_2018D", std::make_shared<TH2D>( "h_jet30_ptem_2018D", "h_jet30_ptem_2018D", 1440, 0, 1440, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet30_pthad_2018D", std::make_shared<TH2D>( "h_jet30_pthad_2018D", "h_jet30_pthad_2018D", 1440, 0, 1440, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet30_etaem_2018D", std::make_shared<TH2D>( "h_jet30_etaem_2018D", "h_jet30_etaem_2018D", 1440, -5, 5, 1440, 0, 1) );
    my_2d_histos.emplace( "h_jet30_etahad_2018D", std::make_shared<TH2D>( "h_jet30_etahad_2018D", "h_jet30_etahad_2018D", 1440, -5, 5, 1440, 0, 1) );
}

//Put everything you want to do per event here.
void AnalyzeEENoise::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {

        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );
       
        const auto& RunNum              = tr.getVar<unsigned int>("RunNum");
        const auto& Jets                = tr.getVec<TLorentzVector>("Jets");
        const auto& NVtx                = tr.getVar<int>("NVtx");
        const auto& fwm2_top6           = tr.getVar<double>("fwm2_top6");
        const auto& fwm3_top6           = tr.getVar<double>("fwm3_top6");
        const auto& fwm4_top6           = tr.getVar<double>("fwm4_top6");
        const auto& fwm5_top6           = tr.getVar<double>("fwm5_top6");
        const auto& jmt_ev0_top6        = tr.getVar<double>("jmt_ev0_top6");
        const auto& jmt_ev1_top6        = tr.getVar<double>("jmt_ev1_top6");
        const auto& jmt_ev2_top6        = tr.getVar<double>("jmt_ev2_top6");
        const auto& event_beta_z        = tr.getVar<double>("event_beta_z");
        const auto& event_beta_z_pt20   = tr.getVar<double>("event_beta_z_pt20");
        const auto& passBaseline        = tr.getVar<bool>("passBaseline1l_Good");
        const auto& jetNeutralEM        = tr.getVec<double>("Jets_neutralEmEnergyFraction");
        const auto& jetNeutralHad       = tr.getVec<double>("Jets_neutralHadronEnergyFraction");
        const auto& nnDisc              = tr.getVar<double>("deepESM_val");
     
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        //Make cuts and fill histograms here
        if ( passBaseline ) {

            if ( RunNum >= 315252 && RunNum <= 316995 ) {

                for (unsigned int iJet = 0; iJet < Jets.size(); iJet++) {
                    my_2d_histos["h_jet_pteta_2018A"]->Fill(Jets[iJet].Pt(), Jets[iJet].Eta(), 1.0);
                    my_2d_histos["h_jet_etaphi_2018A"]->Fill(Jets[iJet].Eta(), Jets[iJet].Phi(), 1.0);
                    my_2d_histos["h_jet_ptem_2018A"]->Fill(Jets[iJet].Pt(), jetNeutralEM[iJet], 1.0);
                    my_2d_histos["h_jet_pthad_2018A"]->Fill(Jets[iJet].Pt(), jetNeutralHad[iJet], 1.0);
                    my_2d_histos["h_jet_etaem_2018A"]->Fill(Jets[iJet].Eta(), jetNeutralEM[iJet], 1.0);
                    my_2d_histos["h_jet_etahad_2018A"]->Fill(Jets[iJet].Eta(), jetNeutralHad[iJet], 1.0);

                    my_histos["h_jet_pt_2018A"]->Fill(Jets[iJet].Pt(), 1.0);
                    my_histos["h_jet_eta_2018A"]->Fill(Jets[iJet].Eta(), 1.0);

                    if (Jets[iJet].Pt() > 30) {
                        my_2d_histos["h_jet30_pteta_2018A"]->Fill(Jets[iJet].Pt(), Jets[iJet].Eta(), 1.0);
                        my_2d_histos["h_jet30_etaphi_2018A"]->Fill(Jets[iJet].Eta(), Jets[iJet].Phi(), 1.0);
                        my_2d_histos["h_jet30_ptem_2018A"]->Fill(Jets[iJet].Pt(), jetNeutralEM[iJet], 1.0);
                        my_2d_histos["h_jet30_pthad_2018A"]->Fill(Jets[iJet].Pt(), jetNeutralHad[iJet], 1.0);
                        my_2d_histos["h_jet30_etaem_2018A"]->Fill(Jets[iJet].Eta(), jetNeutralEM[iJet], 1.0);
                        my_2d_histos["h_jet30_etahad_2018A"]->Fill(Jets[iJet].Eta(), jetNeutralHad[iJet], 1.0);

                        my_histos["h_jet30_pt_2018A"]->Fill(Jets[iJet].Pt(), 1.0);
                        my_histos["h_jet30_eta_2018A"]->Fill(Jets[iJet].Eta(), 1.0);
                    }

                }
               
                my_histos["h_fwm2_top6_2018A"]->Fill(fwm2_top6, 1.0);
                my_histos["h_fwm3_top6_2018A"]->Fill(fwm3_top6, 1.0);
                my_histos["h_fwm4_top6_2018A"]->Fill(fwm4_top6, 1.0);
                my_histos["h_fwm5_top6_2018A"]->Fill(fwm5_top6, 1.0);
                my_histos["h_jmt_ev0_top6_2018A"]->Fill(jmt_ev0_top6, 1.0);
                my_histos["h_jmt_ev1_top6_2018A"]->Fill(jmt_ev1_top6, 1.0);
                my_histos["h_jmt_ev2_top6_2018A"]->Fill(jmt_ev2_top6, 1.0);
                my_histos["h_beta_z_2018A"]->Fill(event_beta_z, 1.0);
                my_histos["h_beta_z_pt20_2018A"]->Fill(event_beta_z_pt20, 1.0);

                my_2d_histos["h_beta_z_nvtx_2018A"]->Fill(event_beta_z, NVtx, 1.0);
                my_2d_histos["h_beta_z_pt20_nvtx_2018A"]->Fill(event_beta_z_pt20, NVtx, 1.0);

                my_histos["h_deepESM_2018A"]->Fill(nnDisc, 1.0);

            } else if ( RunNum >= 320673 && RunNum <= 325175 ) {

                for (unsigned int iJet = 0; iJet < Jets.size(); iJet++) {
                    my_2d_histos["h_jet_pteta_2018D"]->Fill(Jets[iJet].Pt(), Jets[iJet].Eta(), 1.0);
                    my_2d_histos["h_jet_etaphi_2018D"]->Fill(Jets[iJet].Eta(), Jets[iJet].Phi(), 1.0);
                    my_2d_histos["h_jet_ptem_2018D"]->Fill(Jets[iJet].Pt(), jetNeutralEM[iJet], 1.0);
                    my_2d_histos["h_jet_pthad_2018D"]->Fill(Jets[iJet].Pt(), jetNeutralHad[iJet], 1.0);
                    my_2d_histos["h_jet_etaem_2018D"]->Fill(Jets[iJet].Eta(), jetNeutralEM[iJet], 1.0);
                    my_2d_histos["h_jet_etahad_2018D"]->Fill(Jets[iJet].Eta(), jetNeutralHad[iJet], 1.0);

                    my_histos["h_jet_pt_2018D"]->Fill(Jets[iJet].Pt(), 1.0);
                    my_histos["h_jet_eta_2018D"]->Fill(Jets[iJet].Eta(), 1.0);

                    if (Jets[iJet].Pt() > 30) {
                        my_2d_histos["h_jet30_pteta_2018D"]->Fill(Jets[iJet].Pt(), Jets[iJet].Eta(), 1.0);
                        my_2d_histos["h_jet30_etaphi_2018D"]->Fill(Jets[iJet].Eta(), Jets[iJet].Phi(), 1.0);
                        my_2d_histos["h_jet30_ptem_2018D"]->Fill(Jets[iJet].Pt(), jetNeutralEM[iJet], 1.0);
                        my_2d_histos["h_jet30_pthad_2018D"]->Fill(Jets[iJet].Pt(), jetNeutralHad[iJet], 1.0);
                        my_2d_histos["h_jet30_etaem_2018D"]->Fill(Jets[iJet].Eta(), jetNeutralEM[iJet], 1.0);
                        my_2d_histos["h_jet30_etahad_2018D"]->Fill(Jets[iJet].Eta(), jetNeutralHad[iJet], 1.0);

                        my_histos["h_jet30_pt_2018D"]->Fill(Jets[iJet].Pt(), 1.0);
                        my_histos["h_jet30_eta_2018D"]->Fill(Jets[iJet].Eta(), 1.0);
                    }
                }

                my_histos["h_fwm2_top6_2018D"]->Fill(fwm2_top6, 1.0);
                my_histos["h_fwm3_top6_2018D"]->Fill(fwm3_top6, 1.0);
                my_histos["h_fwm4_top6_2018D"]->Fill(fwm4_top6, 1.0);
                my_histos["h_fwm5_top6_2018D"]->Fill(fwm5_top6, 1.0);
                my_histos["h_jmt_ev0_top6_2018D"]->Fill(jmt_ev0_top6, 1.0);
                my_histos["h_jmt_ev1_top6_2018D"]->Fill(jmt_ev1_top6, 1.0);
                my_histos["h_jmt_ev2_top6_2018D"]->Fill(jmt_ev2_top6, 1.0);
                my_histos["h_beta_z_2018D"]->Fill(event_beta_z, 1.0);
                my_histos["h_beta_z_pt20_2018D"]->Fill(event_beta_z_pt20, 1.0);

                my_2d_histos["h_beta_z_nvtx_2018D"]->Fill(event_beta_z, NVtx, 1.0);
                my_2d_histos["h_beta_z_pt20_nvtx_2018D"]->Fill(event_beta_z_pt20, NVtx, 1.0);

                my_histos["h_deepESM_2018D"]->Fill(nnDisc, 1.0);

            }
        }
    } 
}

void AnalyzeEENoise::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
}
