#define ExploreEventSelection_cxx
#include "Analyzer/Analyzer/include/ExploreEventSelection.h"
#include "Framework/Framework/include/Utility.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "Framework/Framework/include/SetUpTopTagger.h"

// includes for the event shapes
#include "Framework/Framework/include/bdt_350to650_fwm10_jmtev_top6.h"
#include "Framework/Framework/include/EventShapeVariables.h"
#include "Framework/Framework/src/get_cmframe_jets.c"
//#include "Framework/Framework/include/fisher_350to650_fwm10_jmtev_top6.h"

void ExploreEventSelection::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("h_met", new TH1D("h_met","h_met", 20, 0, 200));
    my_histos.emplace("h_ht", new TH1D("h_ht","h_ht", 60, 0, 3000));
    my_histos.emplace("h_ntops", new TH1D("h_ntops","h_ntops", 5, 0, 5));
    my_histos.emplace("h_mbl_2l_test", new TH1D("h_mbl_2l_test","h_mbl_2l_test", 50, 0, 200));
    
    //Kelvin added btag type histograms
    std::vector<std::string> btag_types_tags        { "medium" , "tight" };
    std::vector<std::string> pt_threshold_tags      { "pt30" , "pt40" , "pt45" };
    std::vector<std::string> nleptons_tags          { "0l" , "1l" , "2l" };
    std::vector<std::string> n_med_btags            { "0", "1", "2", "3", "4", "5", "6" };
    std::vector<std::string> n_tight_btags          { "0", "1", "2", "3", "4", "5", "6" };
    std::vector<std::string> n_subjets_top_tags     { "1", "2", "3" };

    for( std::string btagTag : btag_types_tags ) {
        for( std::string ptTag : pt_threshold_tags ) {
            for( std::string nLepTag : nleptons_tags ) {
                my_histos.emplace( "h_nbtags_"+btagTag+"_"+ptTag+"_"+nLepTag, new TH1D( ( "h_nbtags_"+btagTag+"_"+ptTag+"_"+nLepTag ).c_str(), ( "h_nbtags_"+btagTag+"_"+ptTag+"_"+nLepTag ).c_str(), 8, 0, 8 ) );
            }//END of nleptons_tags
        }//END of ptTag
    }//END of btagTag

    for( std::string nMedTag : n_med_btags ) {
        for( std::string nTightTag: n_tight_btags ) {

            if( nTightTag == nMedTag ) break; //Does not make sense to require more tight tags than medium tags

            for( std::string ptTag : pt_threshold_tags ) {
                for( std::string nLepTag : nleptons_tags ) {
                    my_histos.emplace( "h_ntops_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag_"+nTightTag+"tTag", new TH1D( ("h_ntops_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag_"+nTightTag+"tTag").c_str(), ("h_ntops_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag_"+nTightTag+"tTag").c_str(), 5, 0, 5 ) );
                    my_histos.emplace( "h_njets_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag_"+nTightTag+"tTag", new TH1D( ("h_njets_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag_"+nTightTag+"tTag").c_str(), ("h_njets_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag_"+nTightTag+"tTag").c_str(), 15, 0, 15 ) );
                    for( std::string nSubJetsTag : n_subjets_top_tags ) {
                        my_histos.emplace( "h_ntops_"+nSubJetsTag+"j_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag_"+nTightTag+"tTag", new TH1D( ( "h_ntops_"+nSubJetsTag+"j_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag_"+nTightTag+"tTag" ).c_str(), ("h_ntops_"+nSubJetsTag+"j_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag_"+nTightTag+"tTag" ).c_str(), 5, 0, 5 ) );
                    }//END of n_subjets_top_tags
                }//END of nleptons_tags
            }//END of pt_threshold_tags

        }//END of n_tight_btags
    }//END of n_med_btags

    //Histograms for scanning the number of tight tags given a set number of medium tags
    for( std::string nMedTag : n_med_btags ) {
        for( std::string ptTag : pt_threshold_tags ) {
            for( std::string nLepTag : nleptons_tags ) {
                my_efficiencies.emplace( "event_sel_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag", new TEfficiency( ( "event_sel_"+ptTag+"_"+nLepTag+"_"+nMedTag+"mTag" ).c_str(), "Efficiency wrt previous cut;Tight btag cut;#epsilon", 9, 0, 9 ) );
            }//END of nleptons_tags
        }//END of pt_threshold_tags
    }//END of n_med_btags   

    for( std::string ptTag : pt_threshold_tags ) {
        for( std::string nLepTag : nleptons_tags ) {
            my_2d_histos.emplace( "h_nbtags_"+ptTag+"_"+nLepTag, new TH2D( ( "h_nbtags_"+ptTag+"_"+nLepTag ).c_str(), ( "h_nbtags_"+ptTag+"_"+nLepTag ).c_str(), 8, 0, 8, 8, 0, 8 ) );
        }
    }//END of pt_threshold_tags
}

void ExploreEventSelection::Loop(double weight, int maxevents, std::string type, std::string filetag, bool isQuiet)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes_total = 0, nbytes = 0;

   // make a toptagger object
   TopTagger tt;
   tt.setCfgFile("TopTagger.cfg");

   // Set up Event shape BDT
   std::vector<std::string> inputVarNames_top6 ;
   std::vector<double> bdtInputVals_top6 ;

   {
       std::string vname ;
       vname = "fwm2_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm3_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm4_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm5_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm6_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm7_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm8_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm9_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm10_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "jmt_ev0_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "jmt_ev1_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "jmt_ev2_top6" ; inputVarNames_top6.push_back( vname ) ;
       
       for ( unsigned int i=0; i < inputVarNames_top6.size() ; i++ ) {
           bdtInputVals_top6.push_back( 0.5 ) ; //--- load vector with dummy values.
       } // i
       
   }
   ReadBDT_350to650_fwm10_jmtev_top6 eventshapeBDT( inputVarNames_top6 ) ;
   //ReadFisher_350to650_fwm10_jmtev_top6 read_fisher_350to650_fwm10_jmtev_top6( inputVarNames_top6 ) ;


   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if(maxevents != -1 && jentry >= maxevents) break;

      nbytes = fChain->GetEntry(jentry);   
      nbytes_total += nbytes;

      if ( jentry % 10000 == 0 ) printf("  Event %9llu\n", jentry ) ;

      // Exclude events with MadGraph HT > 100 from the DY inclusive sample
      if(filetag == "DYJetsToLL_M-50_Incl" && madHT > 100) continue;

      // Make sure event weight is not 0 for data
      double eventweight = 1.;
      if(type != "Data")
          eventweight = Weight;
      
      // -----------------
      // check for number of hadronic tops at gen level
      // -----------------
      int nhadWs = 0;
      std::vector<TLorentzVector> hadtops;
      std::vector<TLorentzVector> hadWs;
      std::vector<int> hadtops_idx;
      std::vector<std::vector<const TLorentzVector*> > hadtopdaughters;
      std::vector<TLorentzVector> neutralinos;
      std::vector<TLorentzVector> singlets;
      std::vector<TLorentzVector> singlinos;
      if(type != "Data")
      {
          for ( unsigned int gpi=0; gpi < GenParticles->size() ; gpi++ ) 
          {
              int pdgid = abs( GenParticles_PdgId->at(gpi) ) ;
              int momid = abs( GenParticles_ParentId->at(gpi) ) ;
              int momidx = GenParticles_ParentIdx->at(gpi);
              int status = GenParticles_Status->at(gpi);
              if(pdgid == 1000022 && (status==22 || status == 52))
              {
                  neutralinos.push_back(GenParticles->at(gpi));
              }
              if(pdgid == 5000001 && (status == 22 || status == 52))
              {
                  singlinos.push_back(GenParticles->at(gpi));
              }
              if(pdgid == 5000002 && (status == 22 || status == 52))
              {
                  singlets.push_back(GenParticles->at(gpi));
              }
              if(status == 23 && momid == 24 && pdgid < 6)
              {
                  // Should be the quarks from W decay
                  nhadWs++;
                  // find the top
                  int Wmotherid = GenParticles_ParentId->at(momidx);
                  if (abs(Wmotherid) == 6){
                      int Wmotheridx = GenParticles_ParentIdx->at(momidx);
                      std::vector<int>::iterator found = std::find(hadtops_idx.begin(), hadtops_idx.end(), Wmotheridx);
                      if (found != hadtops_idx.end())
                      {
                          // already found before
                          // std::cout << "Found this top before: " << *found << std::endl;
                          int position = distance(hadtops_idx.begin(),found);
                          // add the daughter to the list
                          hadtopdaughters[position].push_back(&(GenParticles->at(gpi)));
                      } else
                      {
                          // not yet found
                          hadtops_idx.push_back(Wmotheridx);
                          hadtops.push_back(GenParticles->at(Wmotheridx));
                          hadWs.push_back(GenParticles->at(momidx));
                          std::vector<const TLorentzVector*> daughters;
                          daughters.push_back(&(GenParticles->at(gpi)));
                          hadtopdaughters.push_back(daughters);
                          //std::cout << "Found a new top at idx " << Wmotheridx << std::endl;
                      }
                  }
              } 
          }
          // Now check the b quarks (we only want the ones associated with a hadronic W decay for now)
          for ( unsigned int gpi=0; gpi < GenParticles->size() ; gpi++ ) 
          {
              int pdgid = abs( GenParticles_PdgId->at(gpi) ) ;
              int momid = abs( GenParticles_ParentId->at(gpi) ) ;
              int momidx = GenParticles_ParentIdx->at(gpi);
              int status = GenParticles_Status->at(gpi);
              
              if(status == 23 && momid == 6 && pdgid == 5)
              {
                  // found a b quark from top decay, need to add this to the list of daughters
                  std::vector<int>::iterator found = std::find(hadtops_idx.begin(), hadtops_idx.end(), momidx);
                  if (found != hadtops_idx.end())
                  {
                      // already found
                      int position = distance(hadtops_idx.begin(),found);
                      hadtopdaughters[position].push_back(&(GenParticles->at(gpi)));
                      //std::cout << "(b) Found this top before: " << *found << std::endl;
                  } 
                  //else
                  //{
                  // not yet found
                  //std::cout << "(b) Found a new leptonic top at idx " << momidx << std::endl;
                  //}
              }
          }
      }

      // ------------------------------
      // -- Trigger for data
      // ------------------------------
      
      bool passTriggerAllHad = PassTriggerAllHad();
      bool passTriggerMuon = PassTriggerMuon();
      bool passTriggerElectron = PassTriggerElectron();
      if (type == "Data")
      {
          if (filetag == "Data_JetHT" && !passTriggerAllHad) continue;
          if (filetag == "Data_SingleMuon" && !passTriggerMuon) continue;
          if (filetag == "Data_SingleElectron" && !passTriggerElectron) continue;
      }

      // ------------------
      // --- TOP TAGGER ---
      // ------------------
      
      // setup variables needed for top tagger
      SetUpTopTagger st(*static_cast<const NtupleClass*> (this) , hadtops, hadtopdaughters);
      std::vector<Constituent> constituents = st.getConstituents();
      
      // run the top tagger
      tt.runTagger(constituents);

      // retrieve the top tagger results object
      const TopTaggerResults& ttr = tt.getResults();

      // get reconstructed top
      const std::vector<TopObject*>& tops = ttr.getTops();
      my_histos["h_ntops"]->Fill(tops.size(), eventweight);

      // get set of all constituents (i.e. AK4 and AK8 jets) used in one of the tops
      std::set<Constituent const *> usedConstituents = ttr.getUsedConstituents();

      // count number of tops per type
      int ntops_3jet=0;
      int ntops_2jet=0;
      int ntops_1jet=0;
      for (const TopObject* top : tops)
      {
          if(top->getNConstituents() == 3 )
          {
              ntops_3jet++;
          }
          else if(top->getNConstituents() == 2 )
          {
              ntops_2jet++;
          }
          else if(top->getNConstituents() == 1 )
          {
              ntops_1jet++;
          }
      }


      // -------------------------------
      // -- Event shape BDT
      // -------------------------------

      std::vector<math::RThetaPhiVector> cm_frame_jets ;
      get_cmframe_jets( Jets, cm_frame_jets, 6 ) ;
      EventShapeVariables esv_top6( cm_frame_jets ) ;
      TVectorD eigen_vals_norm_top6 = esv_top6.getEigenValues() ;

      {
          int vi(0) ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(2) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(3) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(4) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(5) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(6) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(7) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(8) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(9) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(10) ; vi++ ;
          bdtInputVals_top6.at(vi) = eigen_vals_norm_top6[0] ; vi++ ;
          bdtInputVals_top6.at(vi) = eigen_vals_norm_top6[1] ; vi++ ;
          bdtInputVals_top6.at(vi) = eigen_vals_norm_top6[2] ; vi++ ;
      }

      double eventshape_bdt_val = eventshapeBDT.GetMvaValue( bdtInputVals_top6 ) ;
      //double fisher_val = read_fisher_350to650_fwm10_jmtev_top6.GetMvaValue( bdtInputVals_top6 ) ;



      // -------------------------------
      // -- Basic event selection stuff
      // -------------------------------

      // Count jets & bjets
      int rec_njet_pt45(0) ;
      int rec_njet_pt30(0) ;
      int rec_njet_pt30_btag(0) ;
      int rec_njet_pt45_btag(0) ;
      double HT_trigger = 0.0;
      std::vector<TLorentzVector> rec_bjets_pt30;
      for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {
          TLorentzVector jlv( Jets->at(rji) ) ;
          if (abs(jlv.Eta()) > 2.4) continue;
          if ( jlv.Pt() > 30 )
          { 
              rec_njet_pt30++;
              if ( Jets_bDiscriminatorCSV->at(rji) > 0.8484) 
              {
                  rec_njet_pt30_btag++;
                  rec_bjets_pt30.push_back(jlv);
              }
          }
          if (jlv.Pt() > 40)
              HT_trigger += jlv.Pt();
          if ( jlv.Pt() > 45 ) 
          {
              rec_njet_pt45++ ;
              if ( Jets_bDiscriminatorCSV->at(rji) > 0.8484) 
                  rec_njet_pt45_btag++;
          }
      } 

      // Count leptons > 30 GeV
      std::vector<TLorentzVector> rec_muon_pt30;
      std::vector<int> rec_charge_muon_pt30;
      for (unsigned int imu = 0; imu < Muons->size(); ++imu)
      {
          TLorentzVector lvmu(Muons->at(imu));
          if( abs(lvmu.Eta()) < 2.4 && lvmu.Pt() > 30 && Muons_passIso->at(imu))
          {
              rec_muon_pt30.push_back(lvmu);
              rec_charge_muon_pt30.push_back(Muons_charge->at(imu));
          }
      }
      std::vector<TLorentzVector> rec_electron_pt30;
      std::vector<int> rec_charge_electron_pt30;
      for (unsigned int iel = 0; iel < Electrons->size(); ++iel)
      {
          TLorentzVector lvel(Electrons->at(iel));
          if( abs(lvel.Eta()) < 2.4 && lvel.Pt() > 30 && Electrons_tightID->at(iel) && Electrons_passIso->at(iel))
          {
              rec_electron_pt30.push_back(lvel);
              rec_charge_electron_pt30.push_back(Electrons_charge->at(iel));
          }
      }
      int nleptons = rec_muon_pt30.size() + rec_electron_pt30.size();
      bool onZ = false;
      bool passMbl_2l = false;
      if ( nleptons == 2 )
      {
          if ( (rec_muon_pt30.size() == 2) && (rec_charge_muon_pt30[0] != rec_charge_muon_pt30[1]) )
          {
              double mll = (rec_muon_pt30[0] + rec_muon_pt30[1]).M();
              if( mll > 81 && mll < 101)
                  onZ = true; 
              // check whether a bl pair passes the M(b,l) cut
              for (TLorentzVector myb : rec_bjets_pt30)
              {
                  double mass_bl_1 = (rec_muon_pt30[0] + myb).M();
                  if(mass_bl_1 < 180 && mass_bl_1 > 30)
                      passMbl_2l = true;
                  double mass_bl_2 = (rec_muon_pt30[1] + myb).M();
                  if(mass_bl_2 < 180 && mass_bl_2 > 30)
                      passMbl_2l = true;
                  my_histos["h_mbl_2l_test"]->Fill(mass_bl_1,eventweight);
                  my_histos["h_mbl_2l_test"]->Fill(mass_bl_2,eventweight);
              }
          } 
          else if ( (rec_electron_pt30.size() == 2) && (rec_charge_electron_pt30[0] != rec_charge_electron_pt30[1]) )
          {
              double mll = (rec_electron_pt30[0] + rec_electron_pt30[1]).M();
              if( mll > 81 && mll < 101)
                  onZ = true;  
              // check whether a bl pair passes the M(b,l) cut
              for (TLorentzVector myb : rec_bjets_pt30)
              {
                  double mass_bl_1 = (rec_electron_pt30[0] + myb).M();
                  if(mass_bl_1 < 180 && mass_bl_1 > 30)
                      passMbl_2l = true;
                  double mass_bl_2 = (rec_electron_pt30[1] + myb).M();
                  if(mass_bl_2 < 180 && mass_bl_2 > 30)
                      passMbl_2l = true;
                  my_histos["h_mbl_2l_test"]->Fill(mass_bl_1,eventweight);
                  my_histos["h_mbl_2l_test"]->Fill(mass_bl_2,eventweight);
              }
          }
      }


      bool passBaseline0l = nleptons==0 && rec_njet_pt45>=6 && HT_trigger > 500 && rec_njet_pt45_btag >= 1;
      bool passBaseline1l = nleptons==1 && rec_njet_pt30>=6 ;
      bool passBaseline1mu = rec_muon_pt30.size()==1 && rec_njet_pt30>=6 ;
      bool passBaseline1el = rec_electron_pt30.size()==1 && rec_njet_pt30>=6 ;
      bool passBaseline2l = nleptons==2;
      if(type == "Data")
      {
          passBaseline0l = passBaseline0l && passTriggerAllHad && (filetag == "Data_JetHT");
          if (rec_muon_pt30.size() > 0)
          {
              passBaseline1l = passBaseline1l && passTriggerMuon && (filetag == "Data_SingleMuon");
              passBaseline2l = passBaseline2l && passTriggerMuon && (filetag == "Data_SingleMuon");
          } 
          else if (rec_electron_pt30.size() > 0)
          {
              passBaseline1l = passBaseline1l && passTriggerElectron && (filetag == "Data_SingleElectron");
              passBaseline2l = passBaseline2l && passTriggerElectron && (filetag == "Data_SingleElectron");
          }
      }
      bool pass_g1b = rec_njet_pt30_btag >= 1;
      bool pass_0t = tops.size()==0, pass_1t = tops.size()==1, pass_2t = tops.size()==2;
      bool pass_1t1 = tops.size()==1 && ntops_1jet==1, pass_1t2 = tops.size()==1 && ntops_2jet==1, pass_1t3 = tops.size()==1 && ntops_3jet==1;
      double mbl = -1;
      TLorentzVector used_bjet;
      double mblmet = -1;
      TLorentzVector metlv;
      metlv.SetPtEtaPhiM(MET, 0, METPhi, 0);
      if(nleptons == 1 && pass_g1b)
      {
          TLorentzVector mylepton = (rec_electron_pt30.size() == 1) ? rec_electron_pt30[0] : rec_muon_pt30[0];
          bool passMtop = false;
          //std::cout << "found lepton and " << rec_njet_pt30_btag << " bjets" << std::endl;
          for (TLorentzVector myb : rec_bjets_pt30)
          {
              double mass_bl = (mylepton + myb).M();
              double mass_blmet = (mylepton + myb + metlv).M();
              //std::cout << "mbl and mblmet are " << mass_bl << " and " << mass_blmet << std::endl;
              if (mbl == -1)
              {
                  mbl = mass_bl;
                  mblmet = mass_blmet;
                  used_bjet = myb;
              }
              else if( abs(mass_bl-172.5) < abs(mbl-172.5))
              {
                  mbl = mass_bl;
                  mblmet = mass_blmet;
                  used_bjet = myb;
              }
          }
      }
      bool pass_mbl = mbl > 30 && mbl < 180;
      // Now check that used_bjet isn't also used for the top tagger
      if (pass_mbl)
      {
          int top_type1_to_remove = 0;
          int top_type2_to_remove = 0;
          int top_type3_to_remove = 0;
          for(const TopObject* mytop : tops)
          {
              const std::vector<Constituent const *> mytop_constituents = mytop->getConstituents();
              bool usedup = false;
              for(const Constituent* c: mytop_constituents)
              {
                  if(c->p() == used_bjet)
                  {
                      usedup = true;
                      //std::cout << "Already used this b for the leptonic top" << std::endl;
                  }
              }
              if (usedup)
              {
                  if (mytop_constituents.size() == 1)
                      top_type1_to_remove++;
                  else if(mytop_constituents.size() == 2)
                      top_type2_to_remove++;
                  else if(mytop_constituents.size() == 3)
                      top_type3_to_remove++;
              }
          }
          int ntops_to_remove = top_type1_to_remove + top_type2_to_remove + top_type3_to_remove;
          //std::cout << "Old top counting: " << pass_0t << " " << pass_1t << " " << pass_2t << " " << pass_1t1 << " " << pass_1t2 << " " << pass_1t3 << std::endl;
          if (ntops_to_remove > 0)
          {
              pass_0t = (tops.size() - ntops_to_remove) == 0;
              pass_1t = (tops.size() - ntops_to_remove) == 1;
              pass_2t = (tops.size() - ntops_to_remove) == 2;
              pass_1t1 = (ntops_1jet - top_type1_to_remove) == 1;
              pass_1t2 = (ntops_2jet - top_type2_to_remove) == 1;
              pass_1t3 = (ntops_3jet - top_type3_to_remove) == 1;
          }
          //std::cout << "New top counting: " << pass_0t << " " << pass_1t << " " << pass_2t << " " << pass_1t1 << " " << pass_1t2 << " " << pass_1t3 << std::endl;
      }

      // Basic event selection stuff
      //
      //   Impose a pT / eta cut ( pT > 20, 30, 40 ) / ( eta < 2.4 )
      //   Apply preselection of njets >= 6 and HT > 500 GeV
      //   Three different CSV btag requirements:
      //       Loose    - 0.5438 ( 10% )
      //       Medium   - 0.8484 (  1% )
      //       Tight    - 0.9535 ( .1% )
      
      int rec_njet_KM_pt45(0);
      int rec_njet_KM_pt40(0);
      int rec_njet_KM_pt30(0);

      int rec_njet_passMediumBtag_pt30(0);
      int rec_njet_passTightBtag_pt30(0);
      
      int rec_njet_passMediumBtag_pt40(0);
      int rec_njet_passTightBtag_pt40(0);

      int rec_njet_passMediumBtag_pt45(0);
      int rec_njet_passTightBtag_pt45(0);

      double HT_pt40(0.0);
      bool pass0lPreselection(true);
      bool pass1lPreselection(true);
      bool pass2lPreselection(true);

      for( unsigned int itJet = 0; itJet < Jets->size(); itJet++ ) {

          TLorentzVector jlv ( Jets->at(itJet) );

          if( std::fabs( jlv.Eta() ) > 2.4 ) continue;

          bool passMediumBtag   = ( Jets_bDiscriminatorCSV->at(itJet) > 0.8484 );
          bool passTightBtag    = ( Jets_bDiscriminatorCSV->at(itJet) > 0.9535 );

          bool passPt30         = ( jlv.Pt() > 30.0 );
          bool passPt40         = ( jlv.Pt() > 40.0 );
          bool passPt45         = ( jlv.Pt() > 35.0 );

          if( passPt30 ) {
            
            rec_njet_KM_pt30++;

            if( passTightBtag ) {
                rec_njet_passTightBtag_pt30++;
                rec_njet_passMediumBtag_pt30++;
            }

            else if( passMediumBtag ) {
                rec_njet_passMediumBtag_pt30++;
            }

          }//END of passPt30
          
          if( passPt40 ) {
            
            HT_pt40 += jlv.Pt();
            rec_njet_KM_pt40++;

            if( passTightBtag ) {
                rec_njet_passTightBtag_pt40++;
                rec_njet_passMediumBtag_pt40++;
            }

            else if( passMediumBtag ) {
                rec_njet_passMediumBtag_pt40++;
            }

          }//END of passPt30
          
          if( passPt45 ) {
            
            rec_njet_KM_pt45++;

            if( passTightBtag ) {
                rec_njet_passTightBtag_pt45++;
                rec_njet_passMediumBtag_pt45++;
            }

            else if( passMediumBtag ) {
                rec_njet_passMediumBtag_pt45++;
            }

          }//END of passPt45

      }//END of jet loop
      //Now filling the histograms
      
      std::vector<std::string> n_medium_btags   { "0", "1", "2", "3", "4", "5", "6" };
      std::vector<std::string> n_tight_btags    { "0", "1", "2", "3", "4", "5", "6" };

      if( passBaseline0l ) {
          
          my_histos["h_nbtags_medium_pt30_0l"]->Fill(rec_njet_passMediumBtag_pt30, eventweight);
          my_histos["h_nbtags_tight_pt30_0l"]->Fill(rec_njet_passTightBtag_pt30, eventweight);
          my_histos["h_nbtags_medium_pt40_0l"]->Fill(rec_njet_passMediumBtag_pt40, eventweight);
          my_histos["h_nbtags_tight_pt40_0l"]->Fill(rec_njet_passTightBtag_pt40, eventweight);
          my_histos["h_nbtags_medium_pt45_0l"]->Fill(rec_njet_passMediumBtag_pt45, eventweight);
          my_histos["h_nbtags_tight_pt45_0l"]->Fill(rec_njet_passTightBtag_pt45, eventweight);

          my_2d_histos["h_nbtags_pt30_0l"]->Fill(rec_njet_passMediumBtag_pt30, rec_njet_passTightBtag_pt30, eventweight);
          my_2d_histos["h_nbtags_pt40_0l"]->Fill(rec_njet_passMediumBtag_pt40, rec_njet_passTightBtag_pt40, eventweight);
          my_2d_histos["h_nbtags_pt45_0l"]->Fill(rec_njet_passMediumBtag_pt45, rec_njet_passTightBtag_pt45, eventweight);
          
          for( int i = 0; i < 6; i++ ) {
            
            //Event Efficiency Plots
            
            if( rec_njet_passMediumBtag_pt30 > i) {
                my_efficiencies["event_sel_pt30_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > -1, 0);
                my_efficiencies["event_sel_pt30_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 0, 1);
                my_efficiencies["event_sel_pt30_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 1, 2);
                my_efficiencies["event_sel_pt30_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 2, 3);
                my_efficiencies["event_sel_pt30_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 3, 4);
                my_efficiencies["event_sel_pt30_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 4, 5);
                my_efficiencies["event_sel_pt30_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 5, 6);
                my_efficiencies["event_sel_pt30_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 6, 7);
                my_efficiencies["event_sel_pt30_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 7, 8);

               for( int j = 0; j < i; ++j ) {
                    
                    if( rec_njet_passTightBtag_pt30 > j ) {

                       int totalTops = ntops_1jet + ntops_2jet +ntops_3jet;
                       my_histos["h_ntops_pt30_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill( totalTops, eventweight);
                       my_histos["h_ntops_1j_pt30_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_1jet, eventweight);
                       my_histos["h_ntops_2j_pt30_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_2jet, eventweight);
                       my_histos["h_ntops_3j_pt30_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_3jet, eventweight);
                       my_histos["h_njets_pt30_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(rec_njet_passTightBtag_pt30, eventweight);
                    }
                }
            }
        
            if( rec_njet_passMediumBtag_pt40 > i) {
                my_efficiencies["event_sel_pt40_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > -1, 0);
                my_efficiencies["event_sel_pt40_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 0, 1);
                my_efficiencies["event_sel_pt40_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 1, 2);
                my_efficiencies["event_sel_pt40_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 2, 3);
                my_efficiencies["event_sel_pt40_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 3, 4);
                my_efficiencies["event_sel_pt40_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 4, 5);
                my_efficiencies["event_sel_pt40_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 5, 6);
                my_efficiencies["event_sel_pt40_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 7, 8);
                for( int j = 0; j < i; ++j ) {
                    if( rec_njet_passTightBtag_pt40 > j ) {
                       int totalTops = ntops_1jet + ntops_2jet +ntops_3jet;
                       my_histos["h_ntops_pt40_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill( totalTops, eventweight);
                       my_histos["h_ntops_1j_pt40_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_1jet, eventweight);
                       my_histos["h_ntops_2j_pt40_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_2jet, eventweight);
                       my_histos["h_ntops_3j_pt40_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_3jet, eventweight);
                       my_histos["h_njets_pt40_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(rec_njet_passTightBtag_pt40, eventweight);
                    }
                }
            }
            if( rec_njet_passMediumBtag_pt45 > i) {
                my_efficiencies["event_sel_pt45_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > -1, 0);
                my_efficiencies["event_sel_pt45_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 0, 1);
                my_efficiencies["event_sel_pt45_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 1, 2);
                my_efficiencies["event_sel_pt45_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 2, 3);
                my_efficiencies["event_sel_pt45_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 3, 4);
                my_efficiencies["event_sel_pt45_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 4, 5);
                my_efficiencies["event_sel_pt45_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 5, 6);
                my_efficiencies["event_sel_pt45_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 6, 7);
                my_efficiencies["event_sel_pt45_0l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 7, 8);
              for( int j = 0; j < i; ++j ) {
                    if( rec_njet_passTightBtag_pt45 > j ) {
                       int totalTops = ntops_1jet + ntops_2jet +ntops_3jet;
                       my_histos["h_ntops_pt45_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill( totalTops, eventweight);
                       my_histos["h_ntops_1j_pt45_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_1jet, eventweight);
                       my_histos["h_ntops_2j_pt45_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_2jet, eventweight);
                       my_histos["h_ntops_3j_pt45_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_3jet, eventweight);
                       my_histos["h_njets_pt45_0l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(rec_njet_passTightBtag_pt45, eventweight);                       
                    }
                }
            }
          }
      }

      if( passBaseline1l ) {
          
          my_histos["h_nbtags_medium_pt30_1l"]->Fill(rec_njet_passMediumBtag_pt30, eventweight);
          my_histos["h_nbtags_tight_pt30_1l"]->Fill(rec_njet_passTightBtag_pt30, eventweight);
          my_histos["h_nbtags_medium_pt40_1l"]->Fill(rec_njet_passMediumBtag_pt40, eventweight);
          my_histos["h_nbtags_tight_pt40_1l"]->Fill(rec_njet_passTightBtag_pt40, eventweight);
          my_histos["h_nbtags_medium_pt45_1l"]->Fill(rec_njet_passMediumBtag_pt45, eventweight);
          my_histos["h_nbtags_tight_pt45_1l"]->Fill(rec_njet_passTightBtag_pt45, eventweight);

          my_2d_histos["h_nbtags_pt30_1l"]->Fill(rec_njet_passMediumBtag_pt30, rec_njet_passTightBtag_pt30, eventweight);
          my_2d_histos["h_nbtags_pt40_1l"]->Fill(rec_njet_passMediumBtag_pt40, rec_njet_passTightBtag_pt40, eventweight);
          my_2d_histos["h_nbtags_pt45_1l"]->Fill(rec_njet_passMediumBtag_pt45, rec_njet_passTightBtag_pt45, eventweight);
          
          for( int i = 0; i < 6; i++ ) {
            
            //Event Efficiency Plots
            
            if( rec_njet_passMediumBtag_pt30 > i) {
                my_efficiencies["event_sel_pt30_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > -1, 0);
                my_efficiencies["event_sel_pt30_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 0, 1);
                my_efficiencies["event_sel_pt30_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 1, 2);
                my_efficiencies["event_sel_pt30_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 2, 3);
                my_efficiencies["event_sel_pt30_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 3, 4);
                my_efficiencies["event_sel_pt30_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 4, 5);
                my_efficiencies["event_sel_pt30_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 5, 6);
                my_efficiencies["event_sel_pt30_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 6, 7);
                my_efficiencies["event_sel_pt30_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 7, 8);

               for( int j = 0; j < i; ++j ) {
                    
                    if( rec_njet_passTightBtag_pt30 > j ) {

                       int totalTops = ntops_1jet + ntops_2jet +ntops_3jet;
                       my_histos["h_ntops_pt30_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill( totalTops, eventweight);
                       my_histos["h_ntops_1j_pt30_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_1jet, eventweight);
                       my_histos["h_ntops_2j_pt30_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_2jet, eventweight);
                       my_histos["h_ntops_3j_pt30_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_3jet, eventweight);
                       my_histos["h_njets_pt30_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(rec_njet_passTightBtag_pt30, eventweight);
                    }
                }
            }
        
            if( rec_njet_passMediumBtag_pt40 > i) {
                my_efficiencies["event_sel_pt40_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > -1, 0);
                my_efficiencies["event_sel_pt40_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 0, 1);
                my_efficiencies["event_sel_pt40_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 1, 2);
                my_efficiencies["event_sel_pt40_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 2, 3);
                my_efficiencies["event_sel_pt40_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 3, 4);
                my_efficiencies["event_sel_pt40_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 4, 5);
                my_efficiencies["event_sel_pt40_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 5, 6);
                my_efficiencies["event_sel_pt40_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 7, 8);
                for( int j = 0; j < i; ++j ) {
                    if( rec_njet_passTightBtag_pt40 > j ) {
                       int totalTops = ntops_1jet + ntops_2jet +ntops_3jet;
                       my_histos["h_ntops_pt40_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill( totalTops, eventweight);
                       my_histos["h_ntops_1j_pt40_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_1jet, eventweight);
                       my_histos["h_ntops_2j_pt40_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_2jet, eventweight);
                       my_histos["h_ntops_3j_pt40_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_3jet, eventweight);
                       my_histos["h_njets_pt40_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(rec_njet_passTightBtag_pt40, eventweight);
                    }
                }
            }
            if( rec_njet_passMediumBtag_pt45 > i) {
                my_efficiencies["event_sel_pt45_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > -1, 0);
                my_efficiencies["event_sel_pt45_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 0, 1);
                my_efficiencies["event_sel_pt45_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 1, 2);
                my_efficiencies["event_sel_pt45_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 2, 3);
                my_efficiencies["event_sel_pt45_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 3, 4);
                my_efficiencies["event_sel_pt45_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 4, 5);
                my_efficiencies["event_sel_pt45_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 5, 6);
                my_efficiencies["event_sel_pt45_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 6, 7);
                my_efficiencies["event_sel_pt45_1l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 7, 8);
              for( int j = 0; j < i; ++j ) {
                    if( rec_njet_passTightBtag_pt45 > j ) {
                       int totalTops = ntops_1jet + ntops_2jet +ntops_3jet;
                       my_histos["h_ntops_pt45_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill( totalTops, eventweight);
                       my_histos["h_ntops_1j_pt45_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_1jet, eventweight);
                       my_histos["h_ntops_2j_pt45_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_2jet, eventweight);
                       my_histos["h_ntops_3j_pt45_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_3jet, eventweight);
                       my_histos["h_njets_pt45_1l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(rec_njet_passTightBtag_pt45, eventweight);                       
                    }
                }
            }
          }
      }

      if( passBaseline2l ) {
          
          my_histos["h_nbtags_medium_pt30_2l"]->Fill(rec_njet_passMediumBtag_pt30, eventweight);
          my_histos["h_nbtags_tight_pt30_2l"]->Fill(rec_njet_passTightBtag_pt30, eventweight);
          my_histos["h_nbtags_medium_pt40_2l"]->Fill(rec_njet_passMediumBtag_pt40, eventweight);
          my_histos["h_nbtags_tight_pt40_2l"]->Fill(rec_njet_passTightBtag_pt40, eventweight);
          my_histos["h_nbtags_medium_pt45_2l"]->Fill(rec_njet_passMediumBtag_pt45, eventweight);
          my_histos["h_nbtags_tight_pt45_2l"]->Fill(rec_njet_passTightBtag_pt45, eventweight);

          my_2d_histos["h_nbtags_pt30_2l"]->Fill(rec_njet_passMediumBtag_pt30, rec_njet_passTightBtag_pt30, eventweight);
          my_2d_histos["h_nbtags_pt40_2l"]->Fill(rec_njet_passMediumBtag_pt40, rec_njet_passTightBtag_pt40, eventweight);
          my_2d_histos["h_nbtags_pt45_2l"]->Fill(rec_njet_passMediumBtag_pt45, rec_njet_passTightBtag_pt45, eventweight);
          
          for( int i = 0; i < 6; i++ ) {
            
            //Event Efficiency Plots
            
            if( rec_njet_passMediumBtag_pt30 > i) {
                my_efficiencies["event_sel_pt30_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > -1, 0);
                my_efficiencies["event_sel_pt30_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 0, 1);
                my_efficiencies["event_sel_pt30_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 1, 2);
                my_efficiencies["event_sel_pt30_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 2, 3);
                my_efficiencies["event_sel_pt30_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 3, 4);
                my_efficiencies["event_sel_pt30_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 4, 5);
                my_efficiencies["event_sel_pt30_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 5, 6);
                my_efficiencies["event_sel_pt30_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 6, 7);
                my_efficiencies["event_sel_pt30_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt30 > 7, 8);

               for( int j = 0; j < i; ++j ) {
                    
                    if( rec_njet_passTightBtag_pt30 > j ) {

                       int totalTops = ntops_1jet + ntops_2jet +ntops_3jet;
                       my_histos["h_ntops_pt30_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill( totalTops, eventweight);
                       my_histos["h_ntops_1j_pt30_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_1jet, eventweight);
                       my_histos["h_ntops_2j_pt30_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_2jet, eventweight);
                       my_histos["h_ntops_3j_pt30_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_3jet, eventweight);
                       my_histos["h_njets_pt30_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(rec_njet_passTightBtag_pt30, eventweight);
                    }
                }
            }
        
            if( rec_njet_passMediumBtag_pt40 > i) {
                my_efficiencies["event_sel_pt40_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > -1, 0);
                my_efficiencies["event_sel_pt40_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 0, 1);
                my_efficiencies["event_sel_pt40_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 1, 2);
                my_efficiencies["event_sel_pt40_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 2, 3);
                my_efficiencies["event_sel_pt40_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 3, 4);
                my_efficiencies["event_sel_pt40_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 4, 5);
                my_efficiencies["event_sel_pt40_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 5, 6);
                my_efficiencies["event_sel_pt40_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt40 > 7, 8);
                for( int j = 0; j < i; ++j ) {
                    if( rec_njet_passTightBtag_pt40 > j ) {
                       int totalTops = ntops_1jet + ntops_2jet +ntops_3jet;
                       my_histos["h_ntops_pt40_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill( totalTops, eventweight);
                       my_histos["h_ntops_1j_pt40_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_1jet, eventweight);
                       my_histos["h_ntops_2j_pt40_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_2jet, eventweight);
                       my_histos["h_ntops_3j_pt40_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_3jet, eventweight);
                       my_histos["h_njets_pt40_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(rec_njet_passTightBtag_pt40, eventweight);
                    }
                }
            }
            if( rec_njet_passMediumBtag_pt45 > i) {
                my_efficiencies["event_sel_pt45_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > -1, 0);
                my_efficiencies["event_sel_pt45_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 0, 1);
                my_efficiencies["event_sel_pt45_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 1, 2);
                my_efficiencies["event_sel_pt45_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 2, 3);
                my_efficiencies["event_sel_pt45_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 3, 4);
                my_efficiencies["event_sel_pt45_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 4, 5);
                my_efficiencies["event_sel_pt45_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 5, 6);
                my_efficiencies["event_sel_pt45_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 6, 7);
                my_efficiencies["event_sel_pt45_2l_"+n_medium_btags[i]+"mTag"]->Fill(rec_njet_passTightBtag_pt45 > 7, 8);
              for( int j = 0; j < i; ++j ) {
                    if( rec_njet_passTightBtag_pt45 > j ) {
                       int totalTops = ntops_1jet + ntops_2jet +ntops_3jet;
                       my_histos["h_ntops_pt45_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill( totalTops, eventweight);
                       my_histos["h_ntops_1j_pt45_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_1jet, eventweight);
                       my_histos["h_ntops_2j_pt45_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_2jet, eventweight);
                       my_histos["h_ntops_3j_pt45_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(ntops_3jet, eventweight);
                       my_histos["h_njets_pt45_2l_"+n_medium_btags[i]+"mTag_"+n_tight_btags[j]+"tTag"]->Fill(rec_njet_passTightBtag_pt45, eventweight);                       
                    }
                }
            }
          }
      }

      my_histos["h_met"]->Fill(MET, eventweight);
      my_histos["h_ht"]->Fill(HT, eventweight);


   } // end of event loop

}

bool ExploreEventSelection::PassTriggerGeneral(std::vector<std::string> &mytriggers)
{
    bool passTrigger = false;
    for(unsigned int i=0; i<TriggerNames->size(); ++i)
    {
        if(TriggerPass->at(i) != 1)
            continue;
        std::string trigname = TriggerNames->at(i);
        if( std::any_of(mytriggers.begin(), mytriggers.end(), [&] (std::string s) { return trigname.find(s)!=std::string::npos; }) )
        {
            passTrigger = true;
            break;
        }
    }
    return passTrigger;

}


bool ExploreEventSelection::PassTriggerAllHad()
{
    std::vector<std::string> mytriggers {
        //"HLT_PFHT1050", // 2017 trigger
        //"HLT_PFHT900"
            //"HLT_PFHT380_SixPFJet32_DoublePFBTagCSV", // 2017 trigger
            //"HLT_PFHT430_SixPFJet40_PFBTagCSV", // 2017 trigger
            "HLT_PFHT450_SixJet40_BTagCSV",
            "HLT_PFHT400_SixJet30_DoubleBTagCSV",            
            };
    return PassTriggerGeneral(mytriggers);
}

bool ExploreEventSelection::PassTriggerMuon()
{
    std::vector<std::string> mytriggers {"HLT_IsoMu24","HLT_IsoTkMu24_v"};
    return PassTriggerGeneral(mytriggers);
}

bool ExploreEventSelection::PassTriggerElectron()
{
    std::vector<std::string> mytriggers {"HLT_Ele27_WPTight_Gsf"};
    return PassTriggerGeneral(mytriggers);
}

void ExploreEventSelection::WriteHistos()
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
