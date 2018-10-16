#define AnalyzeNjetsMinusOneCSJetReplacement_cxx
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Analyzer/Analyzer/include/AnalyzeNjetsMinusOneCSJetReplacement.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/DeepEventShape.h"
#include "Framework/Framework/include/histio.h"


#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

//===========================================================================================================================

TH1* get_hist_1d( const char* hname ) {
   TH1* rp = (TH1*) gDirectory -> FindObject( hname ) ;
   if ( rp == 0x0 ) {
      printf("\n\n *** get_hist_1d : can't find %s\n\n", hname ) ;
      gSystem -> Exit(-1) ;
   }
   return rp ;
}

//===========================================================================================================================

TH2* get_hist_2d( const char* hname ) {
   TH2* rp = (TH2*) gDirectory -> FindObject( hname ) ;
   if ( rp == 0x0 ) {
      printf("\n\n *** get_hist_1d : can't find %s\n\n", hname ) ;
      gSystem -> Exit(-1) ;
   }
   return rp ;
}

//===========================================================================================================================

double gen_drjj( int nj ) {

   char hname[1000] ;
   sprintf( hname, "h_mct_njet%02d__drjj_nontop_smallest_mjj", nj ) ;
   TH1* hp = get_hist_1d( hname ) ;
   return hp -> GetRandom() ;

} // gen_drjj

//=====================================================================

double gen_alpha() {

   return gRandom->Uniform( -3.1415926, 3.1415926 ) ;

} // gen_alpha

//=====================================================================

double gen_beta( int nj, double pt0 ) {

   char hname[1000] ;
   sprintf( hname, "h_mct_njet%02d__pt1ratio_vs_pt0", nj ) ;
   TH2* hp2d = get_hist_2d( hname ) ;
   int nbinsx = hp2d -> GetNbinsX() ;
   TAxis* xaxis = hp2d -> GetXaxis() ;
   int xbin = xaxis -> FindBin( pt0 ) ;
   int xfirst(1) ;
   int xlast(nbinsx) ;
   int nbuf = 20 ;
   if ( xbin > nbuf ) xfirst = xbin - nbuf ;
   if ( xbin < (nbinsx-nbuf) ) xlast = xbin + nbuf ;
   if ( xfirst < 1 ) xfirst = 1 ;
   if ( xlast > nbinsx ) xlast = nbinsx ;
   TH1D* projy = hp2d -> ProjectionY( "gen_beta_py", xfirst, xlast ) ;
   return projy->GetRandom() ;

} // gen_beta

//=====================================================================

bool   calc_replacement_jet_p3( const TLorentzVector& j0tlv, TLorentzVector& j1tlv, TLorentzVector& j2tlv, double dr_jj, double alpha, double beta, bool isQuiet ) {

      double pt0 = j0tlv.Pt() ;

      double dphi_jj = dr_jj * cos( alpha ) ;
      double deta_jj = dr_jj * sin( alpha ) ;
      if ( !isQuiet ) printf("\n\n calc_replacement_jet_p3 : dphi_jj = %7.3f, deta_jj = %7.3f\n\n", dphi_jj, deta_jj ) ;

      double pt1 = beta * ( pt0 - 30. - pt0/2. ) + pt0 / 2. ;
      if ( !isQuiet ) printf(" calc_replacement_jet_p3 : pt0 = %7.2f  ,  pt1 = %7.2f\n", pt0, pt1 ) ;
      if ( pt1 < 0 ) return false ;

      double cosdpjj = cos( dphi_jj ) ;
      double sindpjj = sin( dphi_jj ) ;

      //--- quadratic solution with + is always the right one.
      double sqrt_arg = cosdpjj*cosdpjj + pow(pt0/pt1,2) - 1. ;
      if ( sqrt_arg < 0 ) {
         if ( !isQuiet ) printf(" calc_replacement_jet_p3 : ***  negative arg to sqrt : cosdpjj*cosdpjj + pow(pt0/pt1,2) - 1. = %7.3f + %7.3f - 1. = %7.3f\n",
            cosdpjj*cosdpjj, pow(pt0/pt1,2), sqrt_arg ) ;
         return false ;
      }

      double pt2 = -1. * pt1 * cosdpjj + pt1 * sqrt( sqrt_arg ) ;

      if ( !isQuiet ) printf(" calc_replacement_jet_p3 :   pt2 = %7.2f \n", pt2 ) ;

      double dphi1 = atan2( sindpjj , cosdpjj + pt1 / pt2 ) ;

      double dphi2 = dphi_jj - dphi1 ;

      if ( !isQuiet ) printf(" calc_replacement_jet_p3 :  dphi1 = %7.3f \n", dphi1 ) ;
      if ( !isQuiet ) printf(" calc_replacement_jet_p3 :  dphi2 = %7.3f \n", dphi2 ) ;


      double eta0 = j0tlv.Eta() ;
      double phi0 = j0tlv.Phi() ;

      //--------
      //double eta1_start = eta0  ;
      //double eta1_end = eta0 + deta_jj ;
      //--------
      //double eta1_start = eta0 - deta_jj*0.2  ;
      //double eta1_end = eta0 + 1.2*deta_jj ;
      //--------
      double eta1_start = eta0 - 1.5  ;
      double eta1_end = eta0 + 1.5 ;
      double pz0 = pt0 * sinh( eta0 ) ;
      if ( !isQuiet ) printf(" calc_replacement_jet_p3 : pz0 = %7.2f : eta0 = %7.3f, deta_jj = %7.3f, eta_start = %7.3f, eta_end = %7.3f\n", pz0, eta0, deta_jj, eta1_start, eta1_end ) ;

      double scan_eta1[500] ;
      double scan_sumpz_minus_pz0[500] ;
      double previous_sumpz_minus_pz0(0.) ;
      double previous_eta1(0.) ;

      int nscan(100) ;
      double best_eta1(0) ;
      int npoints_gr(0) ;
      for ( int i=0; i<=nscan; i++ ) {

         double this_eta1 = eta1_start + i * ( (eta1_end - eta1_start) / nscan ) ;
         double this_eta2 = this_eta1 - deta_jj ;

         double pz1 = pt1 * sinh( this_eta1 ) ;
         double pz2 = pt2 * sinh( this_eta2 ) ;

         double this_sumpz_minus_pz0 = pz1+pz2 - pz0 ;

         scan_eta1[i] = this_eta1 ;
         scan_sumpz_minus_pz0[i] = this_sumpz_minus_pz0 ;
         npoints_gr ++ ;

         if ( !isQuiet ) printf(" calc_replacement_jet_p3 :      %3d : eta1 = %7.3f , eta2 = %7.3f :  pz1 = %7.2f  pz2 = %7.2f  pz1+pz2 = %7.2f  pz0 = %7.2f  pz1+pz2-pz0 = %7.2f",
             i, this_eta1, this_eta2, pz1, pz2, pz1+pz2, pz0, pz1+pz2-pz0 ) ;

        //-- look for where the sign changes.
         if ( i > 0 ) {
            if ( this_sumpz_minus_pz0 * previous_sumpz_minus_pz0 < 0. ) {
               double delta_dpz = this_sumpz_minus_pz0 - previous_sumpz_minus_pz0 ;
               double delta_eta1 = this_eta1 - previous_eta1 ;
               double deta_over_dpz = 0 ;
               if ( delta_dpz != 0 ) deta_over_dpz = delta_eta1 / delta_dpz ;
               best_eta1 = previous_eta1 - deta_over_dpz * (previous_sumpz_minus_pz0) ;
               if ( !isQuiet ) printf("    estimated best eta1 = %7.3f ", best_eta1 ) ;
               break ;
            }
         }
         if ( !isQuiet ) printf("\n") ;

         previous_sumpz_minus_pz0 = this_sumpz_minus_pz0 ;
         previous_eta1 = this_eta1 ;

      } // i

      if ( !isQuiet ) printf("\n") ;

      double eta1 = best_eta1 ;
      double eta2 = eta1 - deta_jj ;

      double phi1 = phi0 + dphi1 ;
      double phi2 = phi0 - dphi2 ;

     //---- do some checks with 4 vectors.

      double m1 = 0.2*pt1 ; //********* fix this.
      double m2 = 0.2*pt2 ; //********* fix this.

      if ( pt1 > pt2 ) {
         j1tlv.SetPtEtaPhiM( pt1, eta1, phi1, m1 ) ;
         j2tlv.SetPtEtaPhiM( pt2, eta2, phi2, m2 ) ;
      } else {
         if ( !isQuiet ) printf(" calc_replacement_jet_p3 :  j1 pt < j2 pt ???\n") ;
         j1tlv.SetPtEtaPhiM( pt2, eta2, phi2, m2 ) ;
         j2tlv.SetPtEtaPhiM( pt1, eta1, phi1, m1 ) ;
      }

      if ( !isQuiet ) printf( " calc_replacement_jet_p3 :   j1   pt = %7.1f , eta = %7.3f , phi = %7.3f\n", j1tlv.Pt(), j1tlv.Eta(), j1tlv.Phi() ) ;
      if ( !isQuiet ) printf( " calc_replacement_jet_p3 :   j2   pt = %7.1f , eta = %7.3f , phi = %7.3f\n", j2tlv.Pt(), j2tlv.Eta(), j2tlv.Phi() ) ;

      TLorentzVector j1j2tlv = j1tlv + j2tlv ;
      TLorentzVector pconschecktlv = j0tlv - j1tlv - j2tlv ;
      double deltap = pconschecktlv.P() ;
      if ( !isQuiet ) {
         printf( " calc_replacement_jet_p3 :   j1+j2 px = %7.1f  j1+j2 py = %7.1f  j1+j2 pz = %7.1f\n", j1j2tlv.Px(), j1j2tlv.Py(), j1j2tlv.Pz() ) ;
         printf( " calc_replacement_jet_p3 :   j0    px = %7.1f  j0    py = %7.1f  j0    pz = %7.1f\n", j0tlv.Px(), j0tlv.Py(), j0tlv.Pz() ) ;
         printf( " calc_replacement_jet_p3 :   delta P = %7.3f\n", deltap ) ;
         printf( " calc_replacement_jet_p3 :   j1, j2 deltaR = %7.3f  ,  generated deltaR = %7.3f\n", j1tlv.DeltaR( j2tlv ), dr_jj ) ;
      }
      if ( deltap > 1.0 ) {
         if ( !isQuiet ) printf("\n\n *** calc_replacement_jet_p3 failed.  deltap = %8.3f\n\n", deltap ) ;
         return false ;
      }
      if ( fabs( j1tlv.DeltaR( j2tlv ) - dr_jj ) > 0.1 ) {
         if ( !isQuiet ) printf("\n\n *** calc_replacement_jet_p3 failed.  deltaR does not agree.\n\n" ) ;
         return false ;
      }

      return true ;


} // calc_replacement_jet_p3

//=====================================================================









AnalyzeNjetsMinusOneCSJetReplacement::AnalyzeNjetsMinusOneCSJetReplacement() : initHistos(false)
{
   char inhistfile[10000] ;
   sprintf( inhistfile, "njet-control-input-hists-ttjets.root" ) ;
   printf("\n\n") ;
   printf("   pwd : " ) ;
   gDirectory -> pwd() ;
   printf("\n") ;
   printf("   reading in input histograms from %s\n", inhistfile ) ;

   loadHist( inhistfile ) ; // NFG

   printf("   done reading input histograms.\n" ) ;
   ////gDirectory -> ls() ;
   printf("\n\n") ;
}

//===========================================================================================================================


void AnalyzeNjetsMinusOneCSJetReplacement::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    int fB = 200;

    //---- Declare all your histograms here, that way we can fill them for multiple chains.

    char hname[1000] ;

   for ( int nj=7; nj<=15; nj++ ) {

      sprintf( hname, "h_cs_for_njet%02d__pt2_vs_pt1", nj ) ;
      my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname, 160, 0., 200., 160, 0., 200. ) ) ;

      sprintf( hname, "h_cs_for_njet%02d__drjj", nj ) ;
      my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 6. ) ) ;

      sprintf( hname, "h_cs_for_njet%02d__mjj", nj ) ;
      my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 400. ) ) ;

      sprintf( hname, "h_cs_for_njet%02d__generated_drjj", nj ) ;
      my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 6. ) ) ;

      sprintf( hname, "h_cs_for_njet%02d__generated_alpha", nj ) ;
      my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, -3.2, 3.2 ) ) ;

      sprintf( hname, "h_cs_for_njet%02d__generated_beta", nj ) ;
      my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, 0., 4. ) ) ;

      sprintf( hname, "h_cs_for_njet%02d__drjj_vs_mjj", nj ) ;
      my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname, 40, 0., 400., 40, 0., 6. ) ) ;


      sprintf( hname, "h_cs_for_njet%02d__true_njm1_mva_val", nj ) ;
      my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, -0.1, 1.1 ) ) ;

      sprintf( hname, "h_cs_for_njet%02d__alt_mva_val", nj ) ;
      my_histos.emplace( hname, std::make_shared<TH1D>( hname, hname, 80, -0.1, 1.1 ) ) ;

      sprintf( hname, "h_cs_for_njet%02d__alt_mva_val_vs_true_njm1_mva_val", nj ) ;
      my_2d_histos.emplace( hname, std::make_shared<TH2D>( hname, hname, 80, -0.1, 1.1,  80, -0.1, 1.1 ) ) ;

   } // nj



} // InitHistos

//===========================================================================================================================

void AnalyzeNjetsMinusOneCSJetReplacement::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{

    char hname[1000] ;

    if ( !isQuiet ) {
       printf("\n\n AnalyzeNjetsMinusOneCSJetReplacement::Loop :  Starting loop over events.\n\n") ;
    }

    while( tr.getNextEvent() )
    {

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos();
            initHistos = true;
        }

        const auto& passTrigger          = tr.getVar<bool>("passTrigger");
        const auto& passMadHT            = tr.getVar<bool>("passMadHT");
        const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30      = tr.getVar<int>("NGoodBJets_pt30");
        const auto& Mbl                  = tr.getVar<double>("Mbl");
        const auto& Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV") ;
        const auto& GoodJets_pt30 = tr.getVec<bool>("GoodJets_pt30");


        const auto& Jets = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets = tr.getVec<bool>("GoodJets");
        const auto& NGoodJets = tr.getVar<int>("NGoodJets");
        const auto& GoodLeptons = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& NGoodLeptons = tr.getVar<int>("NGoodLeptons");
        const auto& MET = tr.getVar<double>("MET"); 
        const auto& METPhi = tr.getVar<double>("METPhi");
        const auto& runtype = tr.getVar<std::string>("runtype");

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


        double eventweight = 1.0;        
        // Weight from samples.cc
        //eventweight = weight;

        if(runtype == "MC"){
            const auto& Weight  = tr.getVar<double>("Weight");
            double lumi = 35900; // Lumi for 2016
            // Weight from NTuples            
            eventweight = lumi*Weight;
        }



        // Global cuts
        if ( !(passTrigger && passMadHT) ) continue;


      //----- baseline selection

        if ( NGoodLeptons < 1 ) continue ;
        if ( NGoodJets_pt30 < 6 ) continue ;
        if ( NGoodBJets_pt30 < 1 ) continue ;
        if ( Mbl < 30 || Mbl > 180 ) continue ;


        if ( NGoodJets_pt30 > 15 ) continue ; // don't have histograms for > 15 jets.

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


     //----- Count the number of jets in the event.


      if ( !isQuiet ) {
         printf("\n\n Number of jets (pt>30, |eta|<2.4) : %2d\n\n", NGoodJets_pt30 ) ;
      }

      if ( NGoodJets_pt30 < 6 ) continue ;
      if ( NGoodJets_pt30 > 15 ) continue ; // for now, only have hists for up to 15.


     //------- Loop over jets.

      for ( int rji=0; rji<Jets.size(); rji++ ) {

         TLorentzVector j0tlv( Jets.at(rji) ) ;

         if ( !GoodJets_pt30[rji] ) {
            if ( !isQuiet ) printf("   jet %2d is not good.  Pt=%7.1f , Eta = %7.3f\n", rji, j0tlv.Pt(), j0tlv.Eta() ) ;
            continue ;
         }

         if ( Jets_bDiscriminatorCSV.at(rji) > 0.7 ) {
            if ( !isQuiet ) printf("   jet %2d is b tagged.  CSV = %7.3f\n", rji, Jets_bDiscriminatorCSV.at(rji) ) ;
            continue ;
         }

         bool is_lepjet(false) ;
         double dr(99.) ;
         if ( NGoodLeptons > 0 ) {
            dr = GoodLeptons[0].second.DeltaR( j0tlv ) ;
         }
         if ( dr < 0.4 ) {
            is_lepjet = true ;
            if ( !isQuiet ) printf("   jet %2d is close to lepton.  DeltaR = %7.3f\n", rji, dr ) ;
         }

         if ( is_lepjet ) continue ;

         if ( j0tlv.Pt() < 60 ) continue ;

         if ( !isQuiet ) printf("   jet %2d is ok.\n", rji ) ;


            bool keep_going(true) ;
            int ntry(0) ;

            while ( keep_going ) {

               ntry ++ ;
               if ( ntry > 100 ) {
                  if ( !isQuiet ) printf("\n\n *** Tried %d times with this jet and haven't succeeded.  Giving up.\n\n", ntry ) ;
                  keep_going = false ;
                  continue ;
               }

              //--- generate values

               double pt0 = j0tlv.Pt() ;

               double drjj = gen_drjj( NGoodJets_pt30+1 ) ;
               if ( !isQuiet ) { printf( "  jet %2d : generated drjj = %7.3f\n", rji, drjj ) ; }

               double alpha = gen_alpha() ;
               if ( !isQuiet ) { printf( "  jet %2d : generated alpha = %7.3f\n", rji, alpha ) ; }

               double beta = gen_beta( NGoodJets_pt30+1, pt0 ) ;
               if ( !isQuiet ) { printf( "  jet %2d : generated beta = %7.3f\n", rji, beta ) ; }

               sprintf( hname, "h_cs_for_njet%02d__generated_drjj", NGoodJets_pt30+1 ) ;
               my_histos[ hname ] -> Fill( drjj , eventweight ) ;

               sprintf( hname, "h_cs_for_njet%02d__generated_alpha", NGoodJets_pt30+1 ) ;
               my_histos[ hname ] -> Fill( alpha , eventweight ) ;

               sprintf( hname, "h_cs_for_njet%02d__generated_beta", NGoodJets_pt30+1 ) ;
               my_histos[ hname ] -> Fill( beta , eventweight ) ;


              //--- Calculate replacement jet momentum vectors.
               TLorentzVector j1tlv ;
               TLorentzVector j2tlv ;
               bool calc_ok = calc_replacement_jet_p3( j0tlv, j1tlv, j2tlv, drjj, alpha, beta, true ) ;

               if ( ! calc_ok ) {
                  if ( !isQuiet ) printf(" jet %2d : calc_replacement_jet_p3 returned not ok.  Try again.\n\n", rji ) ;
                  continue ;
               }

               if ( !isQuiet ) {
                  printf( "\n" ) ;
                  printf( " jet %2d main loop :   j1   pt = %7.1f , eta = %7.3f , phi = %7.3f\n", rji, j1tlv.Pt(), j1tlv.Eta(), j1tlv.Phi() ) ;
                  printf( " jet %2d main loop :   j2   pt = %7.1f , eta = %7.3f , phi = %7.3f\n", rji, j2tlv.Pt(), j2tlv.Eta(), j2tlv.Phi() ) ;
               }

               if ( j1tlv.Pt() < 30 ) continue ;
               if ( j2tlv.Pt() < 30 ) continue ;
               if ( fabs( j1tlv.Eta() ) > 2.4 ) continue ;
               if ( fabs( j2tlv.Eta() ) > 2.4 ) continue ;

               if ( !isQuiet ) printf(" jet %2d : j1 and j2 passed pt and eta cuts.\n", rji ) ;

              //-- Make sure jets aren't too close to another jet in the event.

               bool too_close(false) ;
               for ( int oji=0; oji<Jets.size(); oji++ ) {

                  if ( !GoodJets[oji] ) continue ;
                  if ( oji == rji ) continue ;

                  TLorentzVector ojtlv( Jets.at(oji) ) ;

                  if ( j1tlv.DeltaR( ojtlv ) < 0.4 ) {
                     if ( !isQuiet ) printf("    j1 too close to jet ind=%2d :  DeltaR = %7.3f\n", oji, j1tlv.DeltaR( ojtlv ) ) ;
                     too_close = true ;
                     break ;
                  }
                  if ( j2tlv.DeltaR( ojtlv ) < 0.4 ) {
                     if ( !isQuiet ) printf("    j2 too close to jet ind=%2d :  DeltaR = %7.3f\n", oji, j2tlv.DeltaR( ojtlv ) ) ;
                     too_close = true ;
                     break ;
                  }


               } // oji
               if ( too_close ) {
                  if ( !isQuiet ) printf("   one of j1 and j2 too close to something else.  Rejecting this pair.\n") ;
                  continue ;
               }

              //-- generate mass values
               {
                  double gr ;
                  gr = gRandom -> Gaus( 0.17, 0.039 ) ;
                  double m1 = gr * j1tlv.Pt() ;
                  gr = gRandom -> Gaus( 0.17, 0.039 ) ;
                  double m2 = gr * j2tlv.Pt() ;
                  j1tlv.SetPtEtaPhiM( j1tlv.Pt(), j1tlv.Eta(), j1tlv.Phi(), m1 ) ;
                  j2tlv.SetPtEtaPhiM( j2tlv.Pt(), j2tlv.Eta(), j2tlv.Phi(), m2 ) ;
               }

               TLorentzVector j1j2tlv = j1tlv + j2tlv ;

               sprintf( hname, "h_cs_for_njet%02d__mjj", NGoodJets_pt30+1 ) ;
               my_histos[ hname ] -> Fill( j1j2tlv.M() , eventweight ) ;

               sprintf( hname, "h_cs_for_njet%02d__drjj", NGoodJets_pt30+1 ) ;
               my_histos[ hname ] -> Fill( j1tlv.DeltaR( j2tlv ) , eventweight ) ;

               sprintf( hname, "h_cs_for_njet%02d__pt2_vs_pt1", NGoodJets_pt30+1 ) ;
               my_2d_histos[ hname ] -> Fill( j1tlv.Pt(), j2tlv.Pt() , eventweight ) ;

               sprintf( hname, "h_cs_for_njet%02d__drjj_vs_mjj", NGoodJets_pt30+1 ) ;
               my_2d_histos[ hname ] -> Fill( j1j2tlv.M(), j1tlv.DeltaR( j2tlv ) , eventweight ) ;



               const auto& deepESM_val = tr.getVar<double>("deepESM_val");
               if ( !isQuiet ) {
                  printf("\n\n  deepESM_val with nominal set of jets:  %7.3f\n", deepESM_val ) ;
                  for ( unsigned int ji=0; ji<Jets.size(); ji++ ) {
                     printf("   jet %2d : Pt = %7.1f  Eta = %7.3f  Phi = %7.3f   %s\n",
                      ji, Jets.at(ji).Pt(), Jets.at(ji).Eta(), Jets.at(ji).Phi(), (GoodJets.at(ji)?"good":"BAD") ) ;
                  } // ji
               }


               sprintf( hname, "h_cs_for_njet%02d__true_njm1_mva_val", NGoodJets_pt30+1 ) ;
               my_histos[ hname ] ->  Fill( deepESM_val , eventweight ) ;


             //--- Recompute the MVA with these jets.

  //////       auto* altJets = new std::vector<TLorentzVector> ;
  //////       auto* altGoodJets = new std::vector<bool> ;
  //////       auto* altGoodLeptons = new std::vector<std::pair<std::string, TLorentzVector>>(GoodLeptons); // is this necessary???
  //////       int   altNGoodJets(0) ;

  //////       for ( int oji=0; oji<Jets.size(); oji++ ) {

  //////          if ( !GoodJets[oji] ) continue ;
  //////          if ( oji == rji ) continue ;

  //////          altJets->push_back( Jets[oji] ) ;
  //////          altGoodJets->push_back( GoodJets[oji] ) ;
  //////          altNGoodJets ++ ;

  //////       } // oji

  //////       altJets->push_back( j1tlv ) ;
  //////       altGoodJets->push_back( true ) ;
  //////       altNGoodJets ++ ;

  //////       altJets->push_back( j2tlv ) ;
  //////       altGoodJets->push_back( true ) ;
  //////       altNGoodJets ++ ;


  //////       NTupleReader newtr;
  //////       newtr.registerDerivedVec("Jets", altJets);
  //////       newtr.registerDerivedVec("GoodJets", altGoodJets);
  //////       newtr.registerDerivedVar("NGoodJets", altNGoodJets);
  //////       newtr.registerDerivedVec("GoodLeptons", altGoodLeptons);
  //////       newtr.registerDerivedVar("NGoodLeptons", NGoodLeptons);
  //////       newtr.registerDerivedVar("MET", MET);
  //////       newtr.registerDerivedVar("METPhi", METPhi);
  //////   
  //////       MakeMVAVariables makeMVAVariables;
  //////       makeMVAVariables(newtr);
  //////       DeepEventShape deepEventShape;
  //////       deepEventShape(newtr);
  //////       const auto& alt_deepESM_val = newtr.getVar<double>("deepESM_val");
  //////       if ( !isQuiet ) {
  //////          printf("\n\n   Alternate MVA evaluation: altNGoodJets = %2d  , deepESM_val = %7.3f\n", altNGoodJets, alt_deepESM_val ) ;
  //////          for ( unsigned int ji=0; ji<altJets->size(); ji++ ) {
  //////             printf("   jet %2d : Pt = %7.1f  Eta = %7.3f  Phi = %7.3f   %s\n",
  //////                 ji, altJets->at(ji).Pt(), altJets->at(ji).Eta(), altJets->at(ji).Phi(), (altGoodJets->at(ji)?"good":"BAD") ) ;
  //////          } // ji
  //////          printf("\n\n") ;
  //////       }


  //////       sprintf( hname, "h_cs_for_njet%02d__alt_mva_val", NGoodJets_pt30+1 ) ;
  //////       my_histos[ hname ] -> Fill( alt_deepESM_val , eventweight ) ;

  //////       sprintf( hname, "h_cs_for_njet%02d__alt_mva_val_vs_true_njm1_mva_val", NGoodJets_pt30+1 ) ;
  //////       my_2d_histos[ hname ] -> Fill( deepESM_val, alt_deepESM_val , eventweight ) ;

  //////       delete altJets ;
  //////       delete altGoodJets ;
  //////       //delete altGoodLeptons ;


               keep_going = false ;

            } // keep going?




      } // rji

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




//    //-------------------------------------------------------------------------------------
//    ///// Owen try something like this instead
//    
//    int index = -1;
//    for(auto jet : Jets)
//    {
//        index++;           
//        if(!GoodJets[index]) continue;
//        auto* newJets = new std::vector<TLorentzVector>(Jets);
//        auto* newGoodJets = new std::vector<bool>(GoodJets);
//        auto* GoodLeptons_ = new std::vector<std::pair<std::string, TLorentzVector>>(GoodLeptons);
//          
//        newJets->push_back(jet);
//        newGoodJets->push_back(true);
//        int newNGoodJets = NGoodJets + 1;
//          
//        NTupleReader newtr;
//        newtr.registerDerivedVec("Jets", newJets);
//        newtr.registerDerivedVec("GoodJets", newGoodJets);
//        newtr.registerDerivedVar("NGoodJets", newNGoodJets);
//        newtr.registerDerivedVec("GoodLeptons", GoodLeptons_);
//        newtr.registerDerivedVar("NGoodLeptons", NGoodLeptons);
//        newtr.registerDerivedVar("MET", MET);
//        newtr.registerDerivedVar("METPhi", METPhi);
//    
//        MakeMVAVariables makeMVAVariables;
//        makeMVAVariables(newtr);
//        DeepEventShape deepEventShape;
//        deepEventShape(newtr);
//        const auto& deepESM_val = newtr.getVar<double>("deepESM_val");
//        ////////////std::cout<<index<<" "<<newNGoodJets<<" "<<deepESM_val<<std::endl;
//        printf("   Testing new MVA evaluation: index = %2d : newNGoodJets  %2d , deepESM_val = %7.3f\n", index, newNGoodJets, deepESM_val ) ;
//        for ( unsigned int ji=0; ji<newJets->size(); ji++ ) {
//           printf("   jet %2d : Pt = %7.1f  Eta = %7.3f  Phi = %7.3f   %s\n",
//               ji, newJets->at(ji).Pt(), newJets->at(ji).Eta(), newJets->at(ji).Phi(), (newGoodJets->at(ji)?"good":"BAD") ) ;
//        } // ji
//        printf("\n\n") ;
//    }





      


        if (!isQuiet) {
           printf("\n\n =============== End of processing for event ==========================\n\n") ;
        }

    } // end of event loop





} // Loop

//===========================================================================================================================


void AnalyzeNjetsMinusOneCSJetReplacement::WriteHistos(TFile* outfile)
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







