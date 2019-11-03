#define AnalyzeNjetsMinusOneCSFillDijetHists_cxx
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Analyzer/Analyzer/include/AnalyzeNjetsMinusOneCSFillDijetHists.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/DeepEventShape.h"


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

//===========================================================================================================================


AnalyzeNjetsMinusOneCSFillDijetHists::AnalyzeNjetsMinusOneCSFillDijetHists() : initHistos(false)
{
}

//===========================================================================================================================


void AnalyzeNjetsMinusOneCSFillDijetHists::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //---- Declare all your histograms here, that way we can fill them for multiple chains.
    char hname[1000] ;

    sprintf( hname, "h_mct__cvs_nontop" ) ;
    my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 60, -0.1, 1.1 ) ) ;

    sprintf( hname, "h_mct__cvs_topb" ) ;
    my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 60, -0.1, 1.1 ) ) ;

    sprintf( hname, "h_mct__cvs_all" ) ;
    my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 60, -0.1, 1.1 ) ) ;



    for ( int nj=7; nj<=15; nj++ ) {

       sprintf( hname, "h_mct_njet%02d__mjj_nontop_smallest_mjj", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 400. ) ) ;

       sprintf( hname, "h_mct_njet%02d__drjj_nontop_smallest_mjj", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 6. ) ) ;

       sprintf( hname, "h_mct_njet%02d__drjj_vs_mjj_nontop_smallest_mjj", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  40, 0., 400., 40, 0., 6. ) ) ;

       sprintf( hname, "h_mct_njet%02d__jjsum_smallest_dr_with_another_jet", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 6. ) ) ;

       sprintf( hname, "h_mct_njet%02d__jjsum_pt", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 600. ) ) ;

       sprintf( hname, "h_mct_njet%02d__jjsum_eta", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, -5., 5. ) ) ;

       sprintf( hname, "h_mct_njet%02d__pt2_over_pt1_nontop_smallest_mjj", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, -0.1, 1.1 ) ) ;

       sprintf( hname, "h_mct_njet%02d__j1_pt", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 400. ) ) ;

       sprintf( hname, "h_mct_njet%02d__j2_pt", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 400. ) ) ;

       sprintf( hname, "h_mct_njet%02d__j1_m", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 200. ) ) ;

       sprintf( hname, "h_mct_njet%02d__j2_m", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 200. ) ) ;

       sprintf( hname, "h_mct_njet%02d__j1_m_vs_pt", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  40, 0., 400., 40, 0., 100. ) ) ;

       sprintf( hname, "h_mct_njet%02d__j2_m_vs_pt", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  40, 0., 400., 40, 0., 100. ) ) ;

       sprintf( hname, "h_mct_njet%02d__j1_m_over_pt", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 0.4 ) ) ;

       sprintf( hname, "h_mct_njet%02d__j2_m_over_pt", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 0.4 ) ) ;

       sprintf( hname, "h_mct_njet%02d__j2pt_vs_j1pt", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  40, 0., 400., 40, 0., 400. ) ) ;

       sprintf( hname, "h_mct_njet%02d__j2eta_vs_j1eta", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  40, -3., 3., 40, -3., 3. ) ) ;

       sprintf( hname, "h_mct_njet%02d__j2phi_vs_j1phi", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  40, -3.2, 3.2, 40, -3.2, 3.2 ) ) ;

       sprintf( hname, "h_mct_njet%02d__deta_vs_dphi", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  40, -2., 2., 40, -2., 2. ) ) ;

       sprintf( hname, "h_mct_njet%02d__deta_vs_dphi_b2", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  80, -2., 2., 80, -2., 2. ) ) ;

       sprintf( hname, "h_mct_njet%02d__deta_vs_dphi_b3", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  80, -4., 4., 80, -4., 4. ) ) ;

       sprintf( hname, "h_mct_njet%02d__deta_vs_dphi_b4", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  80, -8., 8., 80, -8., 8. ) ) ;


       sprintf( hname, "h_mct_njet%02d__ptrank_nontop", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 21, -0.5, 20.5 ) ) ;
       
       sprintf( hname, "h_mct_njet%02d__ptrank_topdau", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 21, -0.5, 20.5 ) ) ;
       
       sprintf( hname, "h_mct_njet%02d__ptrank_smallest_mjj_sum", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 21, -0.5, 20.5 ) ) ;



       sprintf( hname, "h_mct_njet%02d__pt1ratio_vs_pt0", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  160, 0., 400.,   160, -0.2, 4.0 ) ) ;
       
       sprintf( hname, "h_mct_njet%02d__pt1_vs_pt0", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  160, 0., 200.,   160, 0., 200. ) ) ;

       sprintf( hname, "h_mct_njet%02d__pt2_vs_pt1", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  160, 0., 200.,   160, 0., 200. ) ) ;

       sprintf( hname, "h_mct_njet%02d__drjj_vs_pt0", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  160, 0., 400.,   160, 0., 6. ) ) ;

       sprintf( hname, "h_mct_njet%02d__pt1ratio_vs_drjj", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  160, 0., 6., 160, -0.2, 4.0 ) ) ;

       sprintf( hname, "h_mct_njet%02d__pt1ratio_vs_drjj_pt0gt100", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname,  160, 0., 6., 160, -0.2, 4.0 ) ) ;


       
       sprintf( hname, "h_mct_njet%02d__mva_val_for_true_event", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, -0.1, 1.1 ) ) ;

       sprintf( hname, "h_mct_njet%02d__mva_val_with_dijet_merged", nj ) ;
       my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, -0.1, 1.1 ) ) ;

       sprintf( hname, "h_mct_njet%02d__mva_val_djm_vs_true", nj ) ;
       my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname, 80, -0.1, 1.1,  80, -0.1, 1.1 ) ) ;

    

    } // nj

} // InitHistos

//===========================================================================================================================

void AnalyzeNjetsMinusOneCSFillDijetHists::Loop(NTupleReader& tr, double, int maxevents, bool isQuiet)
{

    char hname[1000] ;

    if ( !isQuiet ) {
       printf("\n\n AnalyzeNjetsMinusOneCSFillDijetHists::Loop :  Starting loop over events.\n\n") ;
    }

    while( tr.getNextEvent() )
    {

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos();
            initHistos = true;
        }

        const auto& MET                  = tr.getVar<double>("MET");
        const auto& runtype              = tr.getVar<std::string>("runtype");     
        const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30      = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons");

        const auto& passTrigger          = tr.getVar<bool>("passTrigger");
        const auto& passMadHT            = tr.getVar<bool>("passMadHT");
              auto  passBaseline1l_Good  = tr.getVar<bool>("passBaseline1l_Good");
        const auto& Mbl                  = tr.getVar<double>("Mbl");
                    passBaseline1l_Good  = passBaseline1l_Good && Mbl>30 && Mbl<180;
        const auto& deepESM_val          = tr.getVar<double>("deepESM_val");

        const auto& GenParticles           = tr.getVec<TLorentzVector>("GenParticles") ;
        const auto& GenParticles_PdgId     = tr.getVec<int>("GenParticles_PdgId") ;
        const auto& GenParticles_ParentId  = tr.getVec<int>("GenParticles_ParentId") ;
        const auto& GenParticles_ParentIdx = tr.getVec<int>("GenParticles_ParentIdx") ;
        const auto& GenParticles_Status    = tr.getVec<int>("GenParticles_Status") ;

        const auto& GenJets           = tr.getVec<TLorentzVector>("GenJets") ;

        const auto& GoodJets_pt30 = tr.getVec<bool>("GoodJets_pt30") ;
        const auto& Jets = tr.getVec<TLorentzVector>("Jets") ;
        const auto& Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV") ;

        // ------------------------
        // -- Print event number
        // -----------------------        
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 100 == 0 ) { printf("  Event %i\n", tr.getEvtNum() ) ; fflush( stdout ) ; }


        if ( !isQuiet ) {
           printf("\n\n\n\n =========== Count %9d : Run %9u , Lumi %9u , Event %9ld  ==========================================\n",
             tr.getEvtNum(),
             tr.getVar<unsigned int>("RunNum"),
             tr.getVar<unsigned int>("LumiBlockNum"),
             tr.getVar<long int>("EvtNum")
             ) ;
        }

        // Global cuts
        if ( !(passTrigger && passMadHT) ) continue;


      //----- baseline selection

        if ( NGoodLeptons < 1 ) continue ;
        if ( NGoodJets_pt30 < 7 ) continue ;
        if ( NGoodBJets_pt30 < 1 ) continue ;
        if ( Mbl < 30 || Mbl > 180 ) continue ;


        if ( NGoodJets_pt30 > 15 ) continue ; // don't have histograms for > 15 jets.



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





      //--- Collect relevant GenParticle indices.

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

            TLorentzVector gplv( GenParticles.at(gpi) ) ;
            char pname[100] ;
            char mname[100] ;
            sprintf( pname, "%s", mcname( GenParticles_PdgId.at(gpi) ) ) ;
            sprintf( mname, "%s", mcname( GenParticles_ParentId.at(gpi) ) ) ;
            double eta = 99. ;
            if ( GenParticles.at(gpi).Pt() > 0 ) eta = GenParticles.at(gpi).Eta() ;
            double phi = GenParticles.at(gpi).Phi() ;
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
               if ( int(gpi) == top1_gpi ) printf("  == top1 " ) ;
               if ( int(gpi) == top2_gpi ) printf("  == top2 " ) ;
               if ( int(gpi) == w1_gpi  ) printf("  == W1 " ) ;
               if ( int(gpi) == w2_gpi  ) printf("  == W2 " ) ;
               if ( int(gpi) == w1_lep_gpi  ) printf("  == W1 lep " ) ;
               if ( int(gpi) == w2_lep_gpi  ) printf("  == W2 lep " ) ;
               if ( int(gpi) == b1_gpi ) printf("  == b1 " ) ;
               if ( int(gpi) == b2_gpi ) printf("  == b2 " ) ;
               if ( int(gpi) == w1_qd1_gpi ) printf("  == W1 quark daughter 1 " ) ;
               if ( int(gpi) == w1_qd2_gpi ) printf("  == W1 quark daughter 2 " ) ;
               if ( int(gpi) == w2_qd1_gpi ) printf("  == W2 quark daughter 1 " ) ;
               if ( int(gpi) == w2_qd2_gpi ) printf("  == W2 quark daughter 2 " ) ;
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

            for ( unsigned int tdi=0; tdi < observable_top_daughters.size() ; tdi++ ) {

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

         for ( unsigned int tdi=0; tdi < observable_top_daughters.size() ; tdi++ ) {

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

         for ( unsigned int j1=0; j1<(nontop_rji.size()-1); j1++ ) {
            TLorentzVector j1tlv( Jets.at(nontop_rji[j1]) ) ;
            if ( !isQuiet ) printf("  j1 = %2d : Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", nontop_rji[j1],
                j1tlv.Pt(), j1tlv.Eta(), j1tlv.Phi() ) ;
            for ( unsigned int j2=j1+1; j2<nontop_rji.size(); j2++ ) {
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

         sprintf( hname, "h_mct_njet%02d__mjj_nontop_smallest_mjj", NGoodJets_pt30 ) ;
         my_histos[ hname ] -> Fill( smallest_mjj_val, eventweight ) ;

         sprintf( hname, "h_mct_njet%02d__drjj_nontop_smallest_mjj", NGoodJets_pt30 ) ;
         my_histos[ hname ] -> Fill( smallest_mjj_dr, eventweight ) ;

         sprintf( hname, "h_mct_njet%02d__drjj_vs_mjj_nontop_smallest_mjj", NGoodJets_pt30 ) ;
         my_2d_histos[ hname ] -> Fill( smallest_mjj_val, smallest_mjj_dr, eventweight ) ;


         //------- how close is the sum of the two jets to other jets?

         double smallest_dr_wrt_jj_pair(9999999.) ;
         for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) {

            if ( !GoodJets_pt30[rji] ) continue ;

            if ( int(rji) == smallest_mjj_j1 ) continue ;
            if ( int(rji) == smallest_mjj_j2 ) continue ;

            TLorentzVector jlv( Jets.at(rji) ) ;

            double dr = jlv.DeltaR( smallest_mjj_tlv ) ;
            if ( dr < smallest_dr_wrt_jj_pair ) {
               smallest_dr_wrt_jj_pair = dr ;
            }

         } // rji
         if ( !isQuiet ) {
            printf("     smallest DR of jet pair vector sum with another jet:  DR = %7.3f\n\n", smallest_dr_wrt_jj_pair ) ;
         }

         sprintf( hname, "h_mct_njet%02d__jjsum_smallest_dr_with_another_jet", NGoodJets_pt30 ) ;
         my_histos[ hname ] -> Fill( smallest_dr_wrt_jj_pair, eventweight ) ;

         sprintf( hname, "h_mct_njet%02d__jjsum_pt", NGoodJets_pt30 ) ;
         my_histos[ hname ] -> Fill( smallest_mjj_tlv.Pt(), eventweight ) ;

         sprintf( hname, "h_mct_njet%02d__jjsum_eta", NGoodJets_pt30 ) ;
         my_histos[ hname ] -> Fill( smallest_mjj_tlv.Eta(), eventweight ) ;


      } // more than 1 nontop jet?



      if ( smallest_mjj_j1 < 0 || smallest_mjj_j2 < 0 ) {
         if ( !isQuiet ) {
            printf("\n\n *** don't have a smallest mjj pair.  Skipping the rest of this event.\n\n") ;
         }
         continue ;
      }


      TLorentzVector j1_tlv( Jets[ smallest_mjj_j1 ] ) ;
      TLorentzVector j2_tlv( Jets[ smallest_mjj_j2 ] ) ;


      sprintf( hname, "h_mct_njet%02d__pt2_over_pt1_nontop_smallest_mjj", NGoodJets_pt30 ) ;
      my_histos[ hname ] -> Fill( (j2_tlv.Pt()) / (j1_tlv.Pt()), eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j1_pt", NGoodJets_pt30 ) ;
      my_histos[ hname ] -> Fill( j1_tlv.Pt() , eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j2_pt", NGoodJets_pt30 ) ;
      my_histos[ hname ] -> Fill( j2_tlv.Pt() , eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j1_m", NGoodJets_pt30 ) ;
      my_histos[ hname ] -> Fill( j1_tlv.M() , eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j2_m", NGoodJets_pt30 ) ;
      my_histos[ hname ] -> Fill( j2_tlv.M() , eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j1_m_vs_pt", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( j1_tlv.Pt(), j1_tlv.M(), eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j2_m_vs_pt", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( j2_tlv.Pt(), j2_tlv.M(), eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j1_m_over_pt", NGoodJets_pt30 ) ;
      my_histos[ hname ] -> Fill( j1_tlv.M() / j1_tlv.Pt() , eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j2_m_over_pt", NGoodJets_pt30 ) ;
      my_histos[ hname ] -> Fill( j2_tlv.M() / j2_tlv.Pt() , eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j2pt_vs_j1pt", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( j1_tlv.Pt(), j2_tlv.Pt(), eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j2eta_vs_j1eta", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( j1_tlv.Eta(), j2_tlv.Eta(), eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__j2phi_vs_j1phi", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( j1_tlv.Phi(), j2_tlv.Phi(), eventweight ) ;


      double dphi = j1_tlv.Phi() - j2_tlv.Phi() ;
      if ( dphi > 3.1415926 ) dphi = dphi - 2*3.1415926 ;
      if ( dphi <-3.1415926 ) dphi = dphi + 2*3.1415926 ;

      sprintf( hname, "h_mct_njet%02d__deta_vs_dphi", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( dphi, j1_tlv.Eta() - j2_tlv.Eta(), eventweight ) ;


      sprintf( hname, "h_mct_njet%02d__deta_vs_dphi_b2", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( dphi, j1_tlv.Eta() - j2_tlv.Eta(), eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__deta_vs_dphi_b3", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( dphi, j1_tlv.Eta() - j2_tlv.Eta(), eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__deta_vs_dphi_b4", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( dphi, j1_tlv.Eta() - j2_tlv.Eta(), eventweight ) ;


      double pt0 = smallest_mjj_tlv.Pt() ;
      double pt1 = j1_tlv.Pt() ;
      double pt2 = j2_tlv.Pt() ;
      double pt1ratio = -9. ;
      if ( (pt0 - 30.) > pt0/2. ) {
         pt1ratio = ( pt1 - pt0/2. ) / ( pt0 - 30. - pt0/2. ) ;
      }

      sprintf( hname, "h_mct_njet%02d__pt1ratio_vs_pt0", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( pt0, pt1ratio, eventweight ) ;
      
      sprintf( hname, "h_mct_njet%02d__pt1_vs_pt0", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( pt0, pt1, eventweight ) ;
      
      sprintf( hname, "h_mct_njet%02d__pt2_vs_pt1", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( pt1, pt2, eventweight ) ;
      
      sprintf( hname, "h_mct_njet%02d__drjj_vs_pt0", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( pt0, smallest_mjj_dr , eventweight ) ;
      
      sprintf( hname, "h_mct_njet%02d__pt1ratio_vs_drjj", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( smallest_mjj_dr, pt1ratio, eventweight ) ;
      
      if ( pt0 > 100 ) {
         sprintf( hname, "h_mct_njet%02d__pt1ratio_vs_drjj_pt0gt100", NGoodJets_pt30 ) ;
         my_2d_histos[ hname ] -> Fill( smallest_mjj_dr, pt1ratio, eventweight ) ;
      }
      






      if ( !isQuiet ) {
         printf("\n\n  Reconstructed Jets: %d\n", NGoodJets_pt30 ) ;
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
            for ( unsigned int tdmi=0; tdmi<topdau_rji.size(); tdmi++ ) {
               if ( topdau_rji[tdmi] == int(rji) ) {
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

         sprintf( hname, "h_mct__cvs_all" ) ;
         my_histos[ hname ] -> Fill( Jets_bDiscriminatorCSV[ rji ], eventweight ) ;

      } // rji


      if ( b1_rji >= 0 ) {
         sprintf( hname, "h_mct__cvs_topb" ) ;
         my_histos[ hname ] -> Fill( Jets_bDiscriminatorCSV[ b1_rji ], eventweight ) ;

      }
      if ( b2_rji >= 0 ) {
         sprintf( hname, "h_mct__cvs_topb" ) ;
         my_histos[ hname ] -> Fill( Jets_bDiscriminatorCSV[ b2_rji ], eventweight ) ;
      }

      for ( unsigned int ntji=0; ntji<nontop_rji.size(); ntji++ ) {
         sprintf( hname, "h_mct__cvs_nontop" ) ;
         my_histos[ hname ] -> Fill( Jets_bDiscriminatorCSV[ nontop_rji[ntji] ], eventweight ) ;
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
            for ( unsigned int i=0; i<good_jets_rji_ptrank_order_using_mjj_sum.size(); i++ ) {
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
            for ( unsigned int tdmi=0; tdmi<topdau_rji.size(); tdmi++ ) {
               if ( topdau_rji[tdmi] == rji ) {
                  sprintf( hname, "h_mct_njet%02d__ptrank_topdau", NGoodJets_pt30 ) ;
                  my_histos[ hname ] -> Fill( i, eventweight ) ;
               }
            } //tdmi
            for ( unsigned int ntji=0; ntji<nontop_rji.size(); ntji++ ) {
               if ( nontop_rji[ntji] == rji ) {
                  sprintf( hname, "h_mct_njet%02d__ptrank_nontop", NGoodJets_pt30 ) ;
                  my_histos[ hname ] -> Fill( i, eventweight ) ;
               }
            } // ntji
         } else if ( rji == -2 ) {
            sprintf( hname, "h_mct_njet%02d__ptrank_smallest_mjj_sum", NGoodJets_pt30 ) ;
            my_histos[ hname ] -> Fill( i, eventweight ) ;
         }
         if ( !isQuiet ) printf("  Rank = %d , rji = %2d : Pt = %7.1f\n", i, rji, jtlv.Pt() ) ;
      } // i

      if ( !isQuiet ) {
         printf("\n deepESM_val = %7.3f\n\n", deepESM_val ) ;
      }

      sprintf( hname, "h_mct_njet%02d__mva_val_for_true_event", NGoodJets_pt30 ) ;
      my_histos[ hname ] -> Fill( deepESM_val, eventweight ) ;



             //--- Recompute the MVA with these all jets except replace the pair with the 4-vector sum.

        const auto& GoodJets = tr.getVec<bool>("GoodJets");
        const auto& GoodLeptons = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& METPhi = tr.getVar<double>("METPhi");

               auto* altJets = new std::vector<TLorentzVector> ;
               auto* altGoodJets = new std::vector<bool> ;
               auto* altGoodJets_pt30 = new std::vector<bool> ;
               auto* altGoodLeptons = new std::vector<std::pair<std::string, TLorentzVector>>(GoodLeptons); // is this necessary???
               int   altNGoodJets(0) ;
               int   altNGoodJets_pt30(0) ;

               for ( unsigned int oji=0; oji<Jets.size(); oji++ ) {

                  ////////////////////////if ( !GoodJets[oji] ) continue ; ////-------- Don't do this!
                  if ( int(oji) == smallest_mjj_j1 ) continue ;
                  if ( int(oji) == smallest_mjj_j2 ) continue ;

                  altJets->push_back( Jets[oji] ) ;
                  altGoodJets->push_back( GoodJets[oji] ) ;
                  altGoodJets_pt30->push_back( GoodJets_pt30[oji] ) ;
                  if ( GoodJets[oji] ) altNGoodJets ++ ;
                  if ( GoodJets_pt30[oji] ) altNGoodJets_pt30 ++ ;

               } // oji

               altJets->push_back( smallest_mjj_tlv ) ;
               altGoodJets->push_back( true ) ;
               altGoodJets_pt30->push_back( true ) ;
               altNGoodJets ++ ;
               altNGoodJets_pt30 ++ ;


               NTupleReader newtr;
               newtr.registerDerivedVec("Jets", altJets);
               newtr.registerDerivedVec("GoodJets", altGoodJets);
               newtr.registerDerivedVar("NGoodJets", altNGoodJets);
               newtr.registerDerivedVec("GoodJets_pt30", altGoodJets_pt30);
               newtr.registerDerivedVar("NGoodJets_pt30", altNGoodJets_pt30);
               newtr.registerDerivedVec("GoodLeptons", altGoodLeptons);
               newtr.registerDerivedVar("NGoodLeptons", NGoodLeptons);
               newtr.registerDerivedVar("MET", MET);
               newtr.registerDerivedVar("METPhi", METPhi);
           
               MakeMVAVariables makeMVAVariables(false,"","GoodJets_pt30",false,false);
               makeMVAVariables(newtr);
               DeepEventShape deepEventShape("DeepEventShape.cfg", "keras_frozen.pb", "Info", false);
               deepEventShape(newtr);
               const auto& alt_deepESM_val = newtr.getVar<double>("deepESM_val");
               if ( !isQuiet ) {
                  printf("\n\n   Alternate MVA evaluation: altNGoodJets = %2d  , deepESM_val = %7.3f\n", altNGoodJets, alt_deepESM_val ) ;
                  for ( unsigned int ji=0; ji<altJets->size(); ji++ ) {
                     printf("   jet %2d : Pt = %7.1f  Eta = %7.3f  Phi = %7.3f   %s\n",
                         ji, altJets->at(ji).Pt(), altJets->at(ji).Eta(), altJets->at(ji).Phi(), (altGoodJets->at(ji)?"good":"BAD") ) ;
                  } // ji
                  printf("\n\n") ;
               }

      sprintf( hname, "h_mct_njet%02d__mva_val_with_dijet_merged", NGoodJets_pt30 ) ;
      my_histos[ hname ] -> Fill( alt_deepESM_val, eventweight ) ;

      sprintf( hname, "h_mct_njet%02d__mva_val_djm_vs_true", NGoodJets_pt30 ) ;
      my_2d_histos[ hname ] -> Fill( deepESM_val, alt_deepESM_val, eventweight ) ;



        if (!isQuiet) {
           printf("\n\n =============== End of processing for event ==========================\n\n") ;
        }

    } // end of event loop

} // Loop

//===========================================================================================================================


void AnalyzeNjetsMinusOneCSFillDijetHists::WriteHistos(TFile* outfile)
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

} // WriteHistos

//===========================================================================================================================



const char* mcname( int pdgid ) {
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




