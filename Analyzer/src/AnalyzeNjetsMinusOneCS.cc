#define AnalyzeNjetsMinusOneCS_cxx
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Analyzer/Analyzer/include/AnalyzeNjetsMinusOneCS.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Framework/Framework/include/MakeMVAVariables.h"


#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

char pname[100] ;
const char* mcname( int pdgid ) ;

AnalyzeNjetsMinusOneCS::AnalyzeNjetsMinusOneCS() : initHistos(false)
{
   printf("\n\n AnalyzeNjetsMinusOneCS::AnalyzeNjetsMinusOneCS : Creating instance of DeepEventShape for alternative calculation.\n\n") ;
   des_ = new DeepEventShape( "DeepEventShape_alt.cfg", "Info", "_alt1" ) ;
}

void AnalyzeNjetsMinusOneCS::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    int fB = 200;

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("h_met",      std::make_shared<TH1D>("h_met",     "h_met",      20,  0,    200  ) );
    my_histos.emplace("h_bdt",      std::make_shared<TH1D>("h_bdt",     "h_bdt",      40, -0.5,    0.5) );
    my_histos.emplace("h_fisher",   std::make_shared<TH1D>("h_fisher",  "h_fisher",   fB, -0.5,    0.5) );
    my_histos.emplace("h_deepESM",  std::make_shared<TH1D>("h_deepESM", "h_deepESM",  fB,  0,      1) );
    my_histos.emplace("h_njets",    std::make_shared<TH1D>("h_njets",   "h_njets",    20,  0,     20  ) );
    my_histos.emplace("h_nb",       std::make_shared<TH1D>("h_nb",      "h_nb",       10,  0,     10  ) );
    for(unsigned int i = 1; i <= 7 ; i++) //Bad hard code
    {
        my_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_1l_ge6j_ge1b",  std::make_shared<TH1D>(("Jet_cm_pt_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(), 150, 0, 1500 ));
        my_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_1l_ge6j_ge1b", std::make_shared<TH1D>(("Jet_cm_eta_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(), 100, -6, 6 ));
        my_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_1l_ge6j_ge1b", std::make_shared<TH1D>(("Jet_cm_phi_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(), 80, -4, 4 ));
        my_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_1l_ge6j_ge1b",   std::make_shared<TH1D>(("Jet_cm_m_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(), 20, 0, 200 ));

        my_histos.emplace("Jet_pt_"+std::to_string(i)+"_1l_ge6j_ge1b",  std::make_shared<TH1D>(("Jet_pt_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(),("Jet_pt_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(), 150, 0, 1500 ));
        my_histos.emplace("Jet_eta_"+std::to_string(i)+"_1l_ge6j_ge1b", std::make_shared<TH1D>(("Jet_eta_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(),("Jet_eta_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(), 100, -6, 6 ));
        my_histos.emplace("Jet_phi_"+std::to_string(i)+"_1l_ge6j_ge1b", std::make_shared<TH1D>(("Jet_phi_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(),("Jet_phi_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(), 80, -4, 4 ));
        my_histos.emplace("Jet_m_"+std::to_string(i)+"_1l_ge6j_ge1b",   std::make_shared<TH1D>(("Jet_m_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(),("Jet_m_"+std::to_string(i)+"_1l_ge6j_ge1b").c_str(), 20, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b").c_str(), fB, 0.0, 1.0, 150, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b").c_str(), fB, 0.0, 1.0, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b").c_str(), fB, 0.0, 1.0, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge6j_ge1b").c_str(), fB, 0.0, 1.0, 20, 0, 200));
    }
} // InitHistos

void AnalyzeNjetsMinusOneCS::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{

    if ( !isQuiet ) {
       printf("\n\n AnalyzeNjetsMinusOneCS::Loop :  Starting loop over events.\n\n") ;
    }

    while( tr.getNextEvent() )
    {
        const auto& MET                  = tr.getVar<double>("MET");
        const auto& runtype              = tr.getVar<std::string>("runtype");     
        const auto& filetag              = tr.getVar<std::string>("filetag");
        const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodJets_pt45       = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets           = tr.getVar<int>("NGoodBJets");
        const auto& NGoodBJets_pt45      = tr.getVar<int>("NGoodBJets_pt45");
        const auto& NGoodBJets_pt30      = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons");

        const auto& eventshape_bdt_val   = tr.getVar<double>("eventshape_bdt_val");
        const auto& fisher_val           = tr.getVar<double>("fisher_val");
        const auto& passTrigger          = tr.getVar<bool>("passTrigger");
        const auto& passMadHT            = tr.getVar<bool>("passMadHT");
              auto  passBaseline1l_Good  = tr.getVar<bool>("passBaseline1l_Good");
        const auto& Mbl                  = tr.getVar<double>("Mbl");
                    passBaseline1l_Good  = passBaseline1l_Good && Mbl>30 && Mbl<180;
        const auto& Photons              = tr.getVec<TLorentzVector>("Photons");
        const auto& GoodPhotons          = tr.getVec<bool>("GoodPhotons");
        const auto& NGoodPhotons         = tr.getVar<int>("NGoodPhotons");
        const auto& passBaseline1g_Good  = tr.getVar<bool>("passBaseline1photon_Good"); 

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
        const auto& Jets_top6            = tr.getVec<TLorentzVector>("Jets_top6");

        const auto& GenParticles           = tr.getVec<TLorentzVector>("GenParticles") ;
        const auto& GenParticles_PdgId     = tr.getVec<int>("GenParticles_PdgId") ;
        const auto& GenParticles_ParentId  = tr.getVec<int>("GenParticles_ParentId") ;
        const auto& GenParticles_ParentIdx = tr.getVec<int>("GenParticles_ParentIdx") ;
        const auto& GenParticles_Status    = tr.getVec<int>("GenParticles_Status") ;

        const auto& GenJets           = tr.getVec<TLorentzVector>("GenJets") ;

        const auto& GoodJets_pt30 = tr.getVec<bool>("GoodJets_pt30") ;
        const auto& Jets = tr.getVec<TLorentzVector>("Jets") ;


      //----- baseline selection

        if ( NGoodLeptons < 1 ) continue ;
        if ( NGoodJets_pt30 < 7 ) continue ;
        if ( NGoodBJets_pt30 < 1 ) continue ;
        if ( Mbl < 30 || Mbl > 180 ) continue ;



        // ------------------------
        // -- Print event number
        // -----------------------        
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        if ( !isQuiet ) {
           printf("\n\n\n\n =========== Count %9lld : Run %9u , Lumi %9u , Event %9llu  ==========================================\n",
             tr.getEvtNum(),
             tr.getVar<unsigned int>("RunNum"),
             tr.getVar<unsigned int>("LumiBlockNum"),
             tr.getVar<long int>("EvtNum")
             ) ;
        }



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



        int top1_gpi(-1) ;
        int top2_gpi(-1) ;
        int w1_gpi(-1) ;
        int w2_gpi(-1) ;
        int w1_lep_gpi(-1) ;
        int w2_lep_gpi(-1) ;
        int b1_gpi(-1) ;
        int b2_gpi(-1) ;
        int w1_qd1_gpi(-1) ;
        int w1_qd2_gpi(-1) ;
        int w2_qd1_gpi(-1) ;
        int w2_qd2_gpi(-1) ;

        std::vector<int> observable_top_daughters ;

        for ( unsigned int gpi=0; gpi < GenParticles.size() ; gpi++ ) {

           int spdgid = GenParticles_PdgId.at(gpi) ;
           int pdgid = abs( GenParticles_PdgId.at(gpi) ) ;
           int smomid = GenParticles_ParentId.at(gpi) ;
           int momid = abs( GenParticles_ParentId.at(gpi) ) ;
           int momidx = GenParticles_ParentIdx.at(gpi) ;
           int status = GenParticles_Status.at(gpi) ;

           TLorentzVector gplv( GenParticles.at(gpi) ) ;

           if ( spdgid ==  6  && top1_gpi < 0 ) top1_gpi = gpi ;
           if ( spdgid == -6  && top2_gpi < 0 ) top2_gpi = gpi ;

           if ( spdgid ==  24 && smomid == 6 && w1_gpi < 0 ) { w1_gpi = gpi ; }
           if ( spdgid == -24 && smomid ==-6 && w2_gpi < 0 ) { w2_gpi = gpi ; }

           if ( spdgid ==  5 && b1_gpi < 0 && GenParticles.at(gpi).Pt() > 1. ) { b1_gpi = gpi ; observable_top_daughters.push_back( gpi ) ; }
           if ( spdgid == -5 && b2_gpi < 0 && GenParticles.at(gpi).Pt() > 1. ) { b2_gpi = gpi ; observable_top_daughters.push_back( gpi ) ; }

           if ( ( spdgid==-11 || spdgid==-13 ) && smomid ==  24 && GenParticles.at(gpi).Pt() > 10. && w1_lep_gpi < 0 )  { w1_lep_gpi = gpi ; observable_top_daughters.push_back( gpi ) ; }
           if ( ( spdgid== 11 || spdgid== 13 ) && smomid == -24 && GenParticles.at(gpi).Pt() > 10. && w2_lep_gpi < 0 )  { w2_lep_gpi = gpi ; observable_top_daughters.push_back( gpi ) ; }

           if ( pdgid < 6 && smomid ==  24 && GenParticles.at(gpi).Pt() > 1. ) {
              if ( w1_qd1_gpi < 0 ) {
                 w1_qd1_gpi = gpi ; observable_top_daughters.push_back( gpi ) ;
              } else if ( w1_qd2_gpi < 0 ) {
                 w1_qd2_gpi = gpi ; observable_top_daughters.push_back( gpi ) ;
              }
           }

           if ( pdgid < 6 && smomid == -24 && GenParticles.at(gpi).Pt() > 1. ) {
              if ( w2_qd1_gpi < 0 ) {
                 w2_qd1_gpi = gpi ; observable_top_daughters.push_back( gpi ) ;
              } else if ( w2_qd2_gpi < 0 )  {
                 w2_qd2_gpi = gpi ; observable_top_daughters.push_back( gpi ) ;
              }
           }

        } // gpi


        if ( !isQuiet ) {
           printf("\n\n GenParticles:\n") ;
        }

        for ( unsigned int gpi=0; gpi < GenParticles.size() ; gpi++ ) {

            int spdgid = GenParticles_PdgId.at(gpi) ;
            int pdgid = abs( GenParticles_PdgId.at(gpi) ) ;
            int smomid = GenParticles_ParentId.at(gpi) ;
            int momid = abs( GenParticles_ParentId.at(gpi) ) ;
            int momidx = GenParticles_ParentIdx.at(gpi) ;

            TLorentzVector gplv( GenParticles.at(gpi) ) ;


            char pname[100] ;
            char mname[100] ;
            sprintf( pname, "%s", mcname( GenParticles_PdgId.at(gpi) ) ) ;
            sprintf( mname, "%s", mcname( GenParticles_ParentId.at(gpi) ) ) ;
            double eta = 99. ;
            if ( GenParticles.at(gpi).Pt() > 0 ) eta = GenParticles.at(gpi).Eta() ;
            double phi = GenParticles.at(gpi).Phi() ;
            double pt = GenParticles.at(gpi).Pt() ;
            if ( !isQuiet ) {
               printf("  %3u :  ID=%9d %10s : MomID=%9d %10s MomIdx=%3d status=%2d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f",
                    gpi,
                    GenParticles_PdgId.at(gpi), pname,
                    GenParticles_ParentId.at(gpi), mname, GenParticles_ParentIdx.at(gpi),
                    GenParticles_Status.at(gpi),
                    GenParticles.at(gpi).Pt(),
                    eta,
                    phi
                    ) ;
               if ( gpi == top1_gpi ) printf("  == top1 " ) ;
               if ( gpi == top2_gpi ) printf("  == top2 " ) ;
               if ( gpi  == w1_gpi  ) printf("  == W1 " ) ;
               if ( gpi  == w2_gpi  ) printf("  == W2 " ) ;
               if ( gpi  == w1_lep_gpi  ) printf("  == W1 lep " ) ;
               if ( gpi  == w2_lep_gpi  ) printf("  == W2 lep " ) ;
               if ( gpi  == b1_gpi ) printf("  == b1 " ) ;
               if ( gpi  == b2_gpi ) printf("  == b2 " ) ;
               if ( gpi  == w1_qd1_gpi ) printf("  == W1 quark daughter 1 " ) ;
               if ( gpi  == w1_qd2_gpi ) printf("  == W1 quark daughter 2 " ) ;
               if ( gpi  == w2_qd1_gpi ) printf("  == W2 quark daughter 1 " ) ;
               if ( gpi  == w2_qd2_gpi ) printf("  == W2 quark daughter 2 " ) ;
               printf("\n") ;
            }
         } // gpi

         if ( !isQuiet ) {
            printf( "\n" ) ;
            printf( "   top1 = %3d ,  w1 = %3d \n", top1_gpi, w1_gpi ) ;
            printf( "   top2 = %3d ,  w2 = %3d \n", top2_gpi, w2_gpi ) ;
         }




        //-------- go through gen jets and check for top daughters.


         if ( !isQuiet ) {
            printf("\n\n  GenJets: %lu\n", GenJets.size() ) ;
         }
         for ( unsigned int gji=0; gji < GenJets.size() ; gji++ ) {

            TLorentzVector gjlv( GenJets.at(gji) ) ;

            double smallest_dr(999999.) ;
            int smallest_dr_tdgpi(-1) ;

            for ( int tdi=0; tdi < observable_top_daughters.size() ; tdi++ ) {

               TLorentzVector tdlv( GenParticles.at( observable_top_daughters[tdi] ) ) ;

               double dr = gjlv.DeltaR( tdlv ) ;

               if ( dr < smallest_dr ) {
                  smallest_dr = dr ;
                  smallest_dr_tdgpi = observable_top_daughters[tdi] ;
               }

            } // tdi

            if ( !isQuiet ) {
               printf("  %3d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f   ",
                 gji,
                   GenJets.at(gji).Pt(),
                   GenJets.at(gji).Eta(),
                   GenJets.at(gji).Phi()
               ) ;
               if ( smallest_dr < 0.15 ) {
                  printf("   matches top daughter : DR = %7.3f ,   ", smallest_dr ) ;
                  int gpi = smallest_dr_tdgpi ;
                  if ( gpi  == w1_lep_gpi  ) printf("  == W1 lep " ) ;
                  if ( gpi  == w2_lep_gpi  ) printf("  == W2 lep " ) ;
                  if ( gpi  == b1_gpi ) printf("  == b1 " ) ;
                  if ( gpi  == b2_gpi ) printf("  == b2 " ) ;
                  if ( gpi  == w1_qd1_gpi ) printf("  == W1 quark daughter 1 " ) ;
                  if ( gpi  == w1_qd2_gpi ) printf("  == W1 quark daughter 2 " ) ;
                  if ( gpi  == w2_qd1_gpi ) printf("  == W2 quark daughter 1 " ) ;
                  if ( gpi  == w2_qd2_gpi ) printf("  == W2 quark daughter 2 " ) ;
                  char pname[100] ;
                  sprintf( pname, "%s", mcname( GenParticles_PdgId.at(smallest_dr_tdgpi) ) ) ;
                  printf( " (%s) ", pname ) ;
               }
               printf("\n") ;
            }

         } // gji





    //---------  find which reconstructed jets are associated to top daughters.

      std::vector<int> nontop_rji ;
      std::vector<int> topdau_rji ;
      std::vector<double> topdau_match_dr ;
      std::vector<int> topdau_match_gpi ;
      int         b1_rji(-1) ;
      int         b2_rji(-1) ;

      for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) {

         if ( !GoodJets_pt30[rji] ) continue ;

         TLorentzVector jlv( Jets[rji] ) ;

         double smallest_dr(999999.) ;
         int smallest_dr_tdgpi(-1) ;

         for ( int tdi=0; tdi < observable_top_daughters.size() ; tdi++ ) {

            TLorentzVector tdlv( GenParticles.at( observable_top_daughters[tdi] ) ) ;

            double dr = jlv.DeltaR( tdlv ) ;

            if ( dr < smallest_dr ) {
               smallest_dr = dr ;
               smallest_dr_tdgpi = observable_top_daughters[tdi] ;
            }

         } // tdi

         if ( smallest_dr > 0.2 ) {
            nontop_rji.push_back( rji ) ;
         }
         if ( smallest_dr < 0.15 ) {
            topdau_rji.push_back( rji ) ;
            topdau_match_dr.push_back( smallest_dr ) ;
            topdau_match_gpi.push_back( smallest_dr_tdgpi ) ;
            if ( smallest_dr_tdgpi == b1_gpi ) b1_rji = rji ;
            if ( smallest_dr_tdgpi == b2_gpi ) b2_rji = rji ;
         }

      } // rji






     //---------------- Look at all nontop jet pairs

     
      if ( !isQuiet ) {
         printf("\n\n Non-top jet pairs:\n") ;
      }

      double smallest_mjj_val(99999999.) ;
      int    smallest_mjj_j1(-1) ;
      int    smallest_mjj_j2(-1) ;
      double smallest_mjj_dr(9999999.) ;
      TLorentzVector smallest_mjj_tlv ;

      if ( nontop_rji.size() > 1 ) {

         for ( int j1=0; j1<(nontop_rji.size()-1); j1++ ) {
            TLorentzVector j1tlv( Jets.at(nontop_rji[j1]) ) ;
            if ( !isQuiet ) printf("  j1 = %2d : Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", nontop_rji[j1],
                j1tlv.Pt(), j1tlv.Eta(), j1tlv.Phi() ) ;
            for ( int j2=j1+1; j2<nontop_rji.size(); j2++ ) {
               TLorentzVector j2tlv( Jets.at(nontop_rji[j2]) ) ;
               TLorentzVector jjtlv = j1tlv + j2tlv ;
               if ( !isQuiet ) printf("  j2 = %2d : Pt = %7.1f , Eta = %6.3f, Phi = %6.3f , mjj = %7.1f , dRjj = %7.3f\n", nontop_rji[j2],
                   j2tlv.Pt(), j2tlv.Eta(), j2tlv.Phi(), jjtlv.M(), j1tlv.DeltaR( j2tlv ) ) ;
               if ( jjtlv.M() < smallest_mjj_val ) {
                  smallest_mjj_val = jjtlv.M() ;
                  smallest_mjj_dr = j1tlv.DeltaR( j2tlv ) ;
                  smallest_mjj_j1 = nontop_rji[j1] ;
                  smallest_mjj_j2 = nontop_rji[j2] ;
                  smallest_mjj_tlv = jjtlv ;
               }
            } // j2
         } // ji

         if ( !isQuiet ) {
            printf("\n\n Pair of reconstructed nontop jets with smallest inv. mass:  Mjj = %7.1f ,  dRjj = %7.3f, j1 = %d , j2 = %d\n",
                smallest_mjj_val, smallest_mjj_dr, smallest_mjj_j1, smallest_mjj_j2 ) ;
            printf("     vector sum of jet pair:  Pt = %7.1f,  Eta = %7.3f\n", smallest_mjj_tlv.Pt(), smallest_mjj_tlv.Eta() ) ;
         }

         //h_mjj_nontop_smallest_mjj -> Fill( smallest_mjj_val ) ;
         //h_drjj_nontop_smallest_mjj -> Fill( smallest_mjj_dr ) ;
         //h_drjj_vs_mjj_nontop_smallest_mjj -> Fill( smallest_mjj_val, smallest_mjj_dr ) ;


         //------- how close is the sum of the two jets to other jets?

         double smallest_dr_wrt_jj_pair(9999999.) ;
         for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) {

            if ( !GoodJets_pt30[rji] ) continue ;

            if ( rji == smallest_mjj_j1 ) continue ;
            if ( rji == smallest_mjj_j2 ) continue ;

            TLorentzVector jlv( Jets.at(rji) ) ;

            double dr = jlv.DeltaR( smallest_mjj_tlv ) ;
            if ( dr < smallest_dr_wrt_jj_pair ) {
               smallest_dr_wrt_jj_pair = dr ;
            }

         } // rji
         if ( !isQuiet ) {
            printf("     smallest DR of jet pair vector sum with another jet:  DR = %7.3f\n\n", smallest_dr_wrt_jj_pair ) ;
         }
         //h_jjsum_smallest_dr_with_another_jet -> Fill( smallest_dr_wrt_jj_pair ) ;
         //h_jjsum_pt -> Fill( smallest_mjj_tlv.Pt() ) ;
         //h_jjsum_eta -> Fill( smallest_mjj_tlv.Eta() ) ;


      } // more than 1 nontop jet?



      if ( smallest_mjj_j1 < 0 || smallest_mjj_j2 < 0 ) {
         if ( !isQuiet ) {
            printf("\n\n *** don't have a smallest mjj pair.  Skipping the rest of this event.\n\n") ;
         }
         continue ;
      }


      TLorentzVector j1_tlv( Jets[ smallest_mjj_j1 ] ) ;
      TLorentzVector j2_tlv( Jets[ smallest_mjj_j2 ] ) ;



  /// h_pt2_over_pt1_nontop_smallest_mjj -> Fill(  (j2_tlv.Pt()) / (j1_tlv.Pt()) ) ;
  /// h_pt2_over_pt1_vs_drjj_nontop_smallest_mjj -> Fill( smallest_mjj_dr , (j2_tlv.Pt()) / (j1_tlv.Pt()) ) ;

  /// h_j1_pt -> Fill( j1_tlv.Pt() ) ;
  /// h_j2_pt -> Fill( j2_tlv.Pt() ) ;

  /// h_j1_m -> Fill( j1_tlv.M() ) ;
  /// h_j2_m -> Fill( j2_tlv.M() ) ;

  /// h_j1_m_vs_pt -> Fill( j1_tlv.Pt(), j1_tlv.M() ) ;
  /// h_j2_m_vs_pt -> Fill( j2_tlv.Pt(), j2_tlv.M() ) ;

  /// h_j1_m_over_pt -> Fill( j1_tlv.M() / j1_tlv.Pt() ) ;
  /// h_j2_m_over_pt -> Fill( j2_tlv.M() / j2_tlv.Pt() ) ;

  /// h_j2pt_vs_j1pt -> Fill( j1_tlv.Pt(), j2_tlv.Pt() ) ;
  /// h_j2eta_vs_j1eta -> Fill( j1_tlv.Eta(), j2_tlv.Eta() ) ;
  /// h_j2phi_vs_j1phi -> Fill( j1_tlv.Phi(), j2_tlv.Phi() ) ;

      double dphi = j1_tlv.Phi() - j2_tlv.Phi() ;
      if ( dphi > 3.1415926 ) dphi = dphi - 2*3.1415926 ;
      if ( dphi <-3.1415926 ) dphi = dphi + 2*3.1415926 ;

  /// h_j1j2_deta_vs_dphi    -> Fill( dphi, j1_tlv.Eta() - j2_tlv.Eta() ) ;
  /// h_j1j2_deta_vs_dphi_b2 -> Fill( dphi, j1_tlv.Eta() - j2_tlv.Eta() ) ;
  /// h_j1j2_deta_vs_dphi_b3 -> Fill( dphi, j1_tlv.Eta() - j2_tlv.Eta() ) ;
  /// h_j1j2_deta_vs_dphi_b4 -> Fill( dphi, j1_tlv.Eta() - j2_tlv.Eta() ) ;








      if ( !isQuiet ) {
         printf("\n\n  Reconstructed Jets: %lu\n", NGoodJets_pt30 ) ;
      }

      for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) {

         if ( !GoodJets_pt30[rji] ) continue ;

         TLorentzVector jlv( Jets[rji] ) ;

         if ( !isQuiet ) {
            printf("  %3d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f   ",
              rji,
                Jets[rji].Pt(),
                Jets[rji].Eta(),
                Jets[rji].Phi()
            ) ;
            for ( int tdmi=0; tdmi<topdau_rji.size(); tdmi++ ) {
               if ( topdau_rji[tdmi] == rji ) {
                  printf("   matches top daughter : DR = %7.3f ,   ", topdau_match_dr[tdmi] ) ;
                  int gpi = topdau_match_gpi[tdmi] ;
                  if ( gpi  == w1_lep_gpi  ) printf("  == W1 lep " ) ;
                  if ( gpi  == w2_lep_gpi  ) printf("  == W2 lep " ) ;
                  if ( gpi  == b1_gpi ) printf("  == b1 " ) ;
                  if ( gpi  == b2_gpi ) printf("  == b2 " ) ;
                  if ( gpi  == w1_qd1_gpi ) printf("  == W1 quark daughter 1 " ) ;
                  if ( gpi  == w1_qd2_gpi ) printf("  == W1 quark daughter 2 " ) ;
                  if ( gpi  == w2_qd1_gpi ) printf("  == W2 quark daughter 1 " ) ;
                  if ( gpi  == w2_qd2_gpi ) printf("  == W2 quark daughter 2 " ) ;
                  char pname[100] ;
                  sprintf( pname, "%s", mcname( GenParticles_PdgId.at(gpi) ) ) ;
                  printf( " (%s) ", pname ) ;
               }
            } // tdmi
            printf("\n") ;
         }


      } // rji


      for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) {

         if ( !GoodJets_pt30[rji] ) continue ;

         //h_csv_all -> Fill( Jets_bDiscriminatorCSV[ rji ] ) ;

      } // rji


      if ( b1_rji >= 0 ) {
         //h_csv_topb -> Fill( Jets_bDiscriminatorCSV[ b1_rji ] ) ;
      }
      if ( b2_rji >= 0 ) {
         //h_csv_topb -> Fill( Jets_bDiscriminatorCSV[ b2_rji ] ) ;
      }

      for ( int ntji=0; ntji<nontop_rji.size(); ntji++ ) {
         //h_csv_nontop -> Fill( Jets_bDiscriminatorCSV[ nontop_rji[ntji] ] ) ;
      } // ntji









     //------------ Now rank the reconstructed jets in pt order, excluding the pair of jets and including the vector sum of the pair

      std::vector<int> good_jets_rji_ptrank_order ;
      std::vector<int> good_jets_rji_ptrank_order_using_mjj_sum ;

      for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) {
         if ( !GoodJets_pt30[rji] ) continue ;
         good_jets_rji_ptrank_order.push_back( rji ) ;
      } // rji

      bool mjj_pair_in_list(false) ;
      for ( unsigned int jet_pt_rank=1; jet_pt_rank <= (good_jets_rji_ptrank_order.size()-1); jet_pt_rank++ ) {

         double highest_pt_val(0.) ;
         double highest_pt_rji(-1) ;

         for ( unsigned int goodjetind=0; goodjetind < good_jets_rji_ptrank_order.size() ; goodjetind++ ) {

            int rji = good_jets_rji_ptrank_order[goodjetind] ;
            if ( rji == smallest_mjj_j1 ) continue ;
            if ( rji == smallest_mjj_j2 ) continue ;

            bool skip_this_one(false) ;
            for ( int i=0; i<good_jets_rji_ptrank_order_using_mjj_sum.size(); i++ ) {
               if ( rji == good_jets_rji_ptrank_order_using_mjj_sum[i] ) {
                  skip_this_one = true ;
                  break ;
               }
            } // i
            if ( skip_this_one ) continue ;

            TLorentzVector jtlv( Jets[ rji ] ) ;

            if ( jtlv.Pt() > highest_pt_val ) {
               highest_pt_val = jtlv.Pt() ;
               highest_pt_rji = rji ;
            }

         } // rji
         if ( !mjj_pair_in_list && smallest_mjj_tlv.Pt() > highest_pt_val ) {
            highest_pt_val = smallest_mjj_tlv.Pt() ;
            highest_pt_rji = -2 ;
            mjj_pair_in_list = true ;
         }
         good_jets_rji_ptrank_order_using_mjj_sum.push_back( highest_pt_rji ) ;

      } // jet_pt_rank

      if ( !isQuiet ) {
         printf("\n\n rjis in Pt rank order using sum of rji = %d and rji = %d:\n", smallest_mjj_j1, smallest_mjj_j2 ) ;
      }
      for ( unsigned int i=0; i<good_jets_rji_ptrank_order_using_mjj_sum.size(); i++ ) {
         int rji = good_jets_rji_ptrank_order_using_mjj_sum[i] ;
         TLorentzVector jtlv ;
         if ( rji >= 0 ) {
            jtlv = Jets[rji] ;
         } else {
            jtlv = smallest_mjj_tlv ;
         }
         if ( rji >= 0 ) {
            for ( int tdmi=0; tdmi<topdau_rji.size(); tdmi++ ) {
               if ( topdau_rji[tdmi] == rji ) {
                  ////h_ptrank_topdau -> Fill( i ) ;
               }
            } //tdmi
            for ( int ntji=0; ntji<nontop_rji.size(); ntji++ ) {
               if ( nontop_rji[ntji] == rji ) {
                  ////h_ptrank_nontop -> Fill( i ) ;
               }
            } // ntji
         } else if ( rji == -2 ) {
            ////h_ptrank_smallest_mjj_sum -> Fill( i ) ;
         }
         if ( !isQuiet ) printf("  Rank = %d , rji = %2d : Pt = %7.1f\n", i, rji, jtlv.Pt() ) ;
      } // i

      if ( !isQuiet ) {
         printf("\n deepESM_val = %7.3f\n\n", deepESM_val ) ;
      }





     //--- Try the code for using an alternative jet collection.


      std::string alt_jet_postfix("_alt1") ;

      auto* alt_Jets = new std::vector<TLorentzVector>(); 
      auto* alt_GoodJets = new std::vector<bool>() ;
      int alt_NGoodJets(0) ; 

      for ( unsigned int rji=0; rji<Jets.size(); rji++ ) {
         if ( !GoodJets_pt30[rji] ) continue ;
         alt_Jets->push_back( Jets.at(rji) ) ;
         alt_GoodJets->push_back( true ) ;
         alt_NGoodJets ++ ;
      } // rji
      TLorentzVector test_fakejet_tlv ;
      test_fakejet_tlv.SetPtEtaPhiM( 100., 0.5, 2.0, 5.0 ) ;
      alt_Jets->push_back( test_fakejet_tlv ) ;
      alt_GoodJets->push_back( true ) ;
      alt_NGoodJets ++ ;

      if ( !isQuiet ) {
         printf( "\n\n About to test registering alternative Jet collection vars.\n\n" ) ;
      }
      tr.registerDerivedVec( "Jets"+alt_jet_postfix, alt_Jets ) ;
      tr.registerDerivedVec( "GoodJets"+alt_jet_postfix, alt_GoodJets ) ;
      tr.registerDerivedVar( "NGoodJets"+alt_jet_postfix, alt_NGoodJets ) ;
      if ( !isQuiet ) {
         printf( "\n\n Done registering alternative Jet collection vars.\n\n" ) ;
         printf( "  Calling MakeMVAVariables for this alternative jet collection.\n\n") ;
      }

      MakeMVAVariables mmv( true, alt_jet_postfix, alt_jet_postfix ) ;
      mmv(tr) ;

      if ( !isQuiet ) {
         printf( "  back from MakeMVAVariables for this alternative jet collection.\n\n") ;
         printf( "  trying DeepEventShape calculator with this alternative jet collection.\n\n") ;
      }
      (*des_)(tr) ;
      

      const auto& deepESM_val_alt1 = tr.getVar<double>("deepESM_val_alt1") ;
      if ( !isQuiet ) {
         printf( "     deepESM_val_alt1 = %7.3f\n", deepESM_val_alt1 ) ;
      }
      








        // -------------------------------
        // -- Define cuts
        // -------------------------------

        bool pass_0l              = NGoodLeptons == 0;
        bool pass_1l              = NGoodLeptons == 1;
        bool pass_njet_pt45       = NGoodJets_pt45 >= 7;
        bool pass_njet_pt45_1btag = NGoodBJets_pt45 >= 1;
        bool pass_njet_pt45_2btag = NGoodBJets_pt45 >= 2;




        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map_1l 
        {
            {"1l"                             , pass_1l                                                                  },
            {"1l_ge6j"                        , pass_1l && pass_njet_pt45                                                },
            {"1l_ge2b"                        , pass_1l && pass_njet_pt45_2btag                                          },
            {"1l_ge6j_ge2b"                   , pass_1l && pass_njet_pt45 && pass_njet_pt45_2btag                        },
            {"1l_ge6j_ge1b"                   , passBaseline1l_Good                                                      },                         
            {"1l_ge6j_ge1b_ge5esm"            , passBaseline1l_Good && deepESM_val >= 0.5                                },                         
            {"1l_ge6j_ge1b_ge5-6esm"          , passBaseline1l_Good && deepESM_val >= 0.5 && deepESM_val < 0.6           },                         
            {"1l_ge6j_ge1b_ge6-7esm"          , passBaseline1l_Good && deepESM_val >= 0.6 && deepESM_val < 0.7           },                         
            {"1l_ge6j_ge1b_ge7-8esm"          , passBaseline1l_Good && deepESM_val >= 0.7 && deepESM_val < 0.8           },                         
            {"1l_ge6j_ge1b_ge8-95esm"         , passBaseline1l_Good && deepESM_val >= 0.8 && deepESM_val < 0.95          },                         
            {"1l_ge6j_ge1b_ge8esm"            , passBaseline1l_Good && deepESM_val >= 0.8                                },                         
            {"1l_ge6j_ge1b_ge95esm"           , passBaseline1l_Good && deepESM_val >= 0.95                               },                         
            {"1l_ge6j_ge1b_l8esm"             , passBaseline1l_Good && deepESM_val <  0.8                                },                         
            {"1l_ge6j_ge1b_d1"                , passBaseline1l_Good && deepESM_bin1                                      },                         
            {"1l_ge6j_ge1b_d2"                , passBaseline1l_Good && deepESM_bin2                                      },                         
            {"1l_ge6j_ge1b_d3"                , passBaseline1l_Good && deepESM_bin3                                      },                         
            {"1l_ge6j_ge1b_d4"                , passBaseline1l_Good && deepESM_bin4                                      },                         
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map_1l);
            initHistos = true;
        }

        // Global cuts
        if ( !(passTrigger && passMadHT) ) continue;

    /// for(auto& kv : cut_map_1l)
    /// {
    ///     if(kv.second)
    ///     {
    ///         my_histos["h_njets_"   +kv.first]->Fill(NGoodJets_pt30, eventweight);
    ///         my_histos["h_ntops_"   +kv.first]->Fill(ntops, eventweight);
    ///         my_histos["h_nb_"      +kv.first]->Fill(NGoodBJets, eventweight);
    ///         my_histos["h_bdt_"     +kv.first]->Fill(eventshape_bdt_val, eventweight);
    ///         my_histos["h_fisher_"  +kv.first]->Fill(fisher_val, eventweight);
    ///         my_histos["h_deepESM_" +kv.first]->Fill(deepESM_val, eventweight);
    ///     }
    /// }

        
        // No local cuts applied here
        my_histos["h_met"     ]->Fill(MET, eventweight);
        my_histos["h_bdt"     ]->Fill(eventshape_bdt_val, eventweight);
        my_histos["h_fisher"  ]->Fill(fisher_val, eventweight);
        my_histos["h_deepESM" ]->Fill(deepESM_val, eventweight);
        my_histos["h_njets"   ]->Fill(NGoodJets_pt30, eventweight);
        my_histos["h_nb"      ]->Fill(NGoodBJets, eventweight);


        // Fill histos for deepESM and Fisher training
        if(passBaseline1l_Good)
        {
            for(unsigned int i = 0; i < Jets_cm_top6.size(); i++)
            {
                my_histos["Jet_cm_pt_"+std::to_string(i+1)+"_1l_ge6j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Pt()), eventweight);
                my_histos["Jet_cm_eta_"+std::to_string(i+1)+"_1l_ge6j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Eta()), eventweight);
                my_histos["Jet_cm_phi_"+std::to_string(i+1)+"_1l_ge6j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).Phi()), eventweight);
                my_histos["Jet_cm_m_"+std::to_string(i+1)+"_1l_ge6j_ge1b"]->Fill(static_cast<double>(Jets_cm_top6.at(i).M()), eventweight);

                my_histos["Jet_pt_"+std::to_string(i+1)+"_1l_ge6j_ge1b"]->Fill(static_cast<double>(Jets_top6.at(i).Pt()), eventweight);
                my_histos["Jet_eta_"+std::to_string(i+1)+"_1l_ge6j_ge1b"]->Fill(static_cast<double>(Jets_top6.at(i).Eta()), eventweight);
                my_histos["Jet_phi_"+std::to_string(i+1)+"_1l_ge6j_ge1b"]->Fill(static_cast<double>(Jets_top6.at(i).Phi()), eventweight);
                my_histos["Jet_m_"+std::to_string(i+1)+"_1l_ge6j_ge1b"]->Fill(static_cast<double>(Jets_top6.at(i).M()), eventweight);

                my_2d_histos["Jet_cm_pt_"+std::to_string(i+1)+"_deepESM_1l_ge6j_ge1b"]->Fill(static_cast<double>(deepESM_val, Jets_cm_top6.at(i).Pt()), eventweight);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i+1)+"_deepESM_1l_ge6j_ge1b"]->Fill(static_cast<double>(deepESM_val, Jets_cm_top6.at(i).Eta()), eventweight);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i+1)+"_deepESM_1l_ge6j_ge1b"]->Fill(static_cast<double>(deepESM_val, Jets_cm_top6.at(i).Phi()), eventweight);
                my_2d_histos["Jet_cm_m_"+std::to_string(i+1)+"_deepESM_1l_ge6j_ge1b"]->Fill(static_cast<double>(deepESM_val, Jets_cm_top6.at(i).M()), eventweight);
            }
        }

        if (!isQuiet) {
           printf("\n\n =============== End of processing for event ==========================\n\n") ;
        }

    } // end of event loop

} // Loop

void AnalyzeNjetsMinusOneCS::WriteHistos(TFile* outfile)
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


//---------------------

const char* mcname( int pdgid ) {

   sprintf( pname, "" ) ;

   if ( pdgid == 1 ) sprintf( pname, "d" ) ;
   if ( pdgid == 2 ) sprintf( pname, "u" ) ;
   if ( pdgid == 3 ) sprintf( pname, "s" ) ;
   if ( pdgid == 4 ) sprintf( pname, "c" ) ;
   if ( pdgid == 5 ) sprintf( pname, "b" ) ;
   if ( pdgid == 6 ) sprintf( pname, "t" ) ;

   if ( pdgid == -1 ) sprintf( pname, "d-bar" ) ;
   if ( pdgid == -2 ) sprintf( pname, "u-bar" ) ;
   if ( pdgid == -3 ) sprintf( pname, "s-bar" ) ;
   if ( pdgid == -4 ) sprintf( pname, "c-bar" ) ;
   if ( pdgid == -5 ) sprintf( pname, "b-bar" ) ;
   if ( pdgid == -6 ) sprintf( pname, "t-bar" ) ;

   if ( pdgid == 11 ) sprintf( pname, "e-" ) ;
   if ( pdgid == 12 ) sprintf( pname, "nu_e" ) ;
   if ( pdgid == 13 ) sprintf( pname, "mu-" ) ;
   if ( pdgid == 14 ) sprintf( pname, "nu_mu" ) ;
   if ( pdgid == 15 ) sprintf( pname, "tau-" ) ;
   if ( pdgid == 16 ) sprintf( pname, "nu_tau" ) ;

   if ( pdgid == -11 ) sprintf( pname, "e+" ) ;
   if ( pdgid == -12 ) sprintf( pname, "nu_e-bar" ) ;
   if ( pdgid == -13 ) sprintf( pname, "mu+" ) ;
   if ( pdgid == -14 ) sprintf( pname, "nu_mu-bar" ) ;
   if ( pdgid == -15 ) sprintf( pname, "tau+" ) ;
   if ( pdgid == -16 ) sprintf( pname, "nu_tau-bar" ) ;

   if ( pdgid == 21 ) sprintf( pname, "gluon" ) ;
   if ( pdgid == 22 ) sprintf( pname, "photon" ) ;
   if ( pdgid == 23 ) sprintf( pname, "Z0" ) ;
   if ( pdgid == 24 ) sprintf( pname, "W+" ) ;
   if ( pdgid ==-24 ) sprintf( pname, "W-" ) ;
   if ( pdgid == 25 ) sprintf( pname, "h" ) ;
   if ( pdgid == 35 ) sprintf( pname, "H" ) ;
   if ( pdgid == 36 ) sprintf( pname, "a" ) ;


   if ( pdgid == 1000001 ) sprintf( pname, "~dL" ) ;
   if ( pdgid == 1000002 ) sprintf( pname, "~uL" ) ;
   if ( pdgid == 1000003 ) sprintf( pname, "~sL" ) ;
   if ( pdgid == 1000004 ) sprintf( pname, "~cL" ) ;
   if ( pdgid == 1000005 ) sprintf( pname, "~b1" ) ;
   if ( pdgid == 1000006 ) sprintf( pname, "~t1" ) ;

   if ( pdgid == -1000001 ) sprintf( pname, "~dL*" ) ;
   if ( pdgid == -1000002 ) sprintf( pname, "~uL*" ) ;
   if ( pdgid == -1000003 ) sprintf( pname, "~sL*" ) ;
   if ( pdgid == -1000004 ) sprintf( pname, "~cL*" ) ;
   if ( pdgid == -1000005 ) sprintf( pname, "~b1*" ) ;
   if ( pdgid == -1000006 ) sprintf( pname, "~t1*" ) ;

   if ( pdgid == 1000022 ) sprintf( pname, "~chi01" ) ;

   return pname ;


} // mcname




