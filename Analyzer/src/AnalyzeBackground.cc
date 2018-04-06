#define AnalyzeBackground_cxx
#include "Analyzer/Analyzer/include/AnalyzeBackground.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
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
#include "Framework/Framework/src/fisher_350to650_fwm10_jmtev_top6.c"
#include "Framework/Framework/src/fisher_350to650_fwm6_jmtev_top6_gt_v2.c"

void AnalyzeBackground::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    std::vector<std::string> jettypes {"pt30", "pt45"};
    for(std::string jettype : jettypes)
    {
        std::string base = "h_njets_" + jettype;
        my_histos.emplace(base,new TH1D(base.c_str(),base.c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l",new TH1D( (base+"_0l").c_str(),(base+"_0l").c_str(),17,0,17));
        my_histos.emplace(base + "_1l",new TH1D((base+"_1l").c_str(),(base+"_1l").c_str(),15,0,15));
        my_histos.emplace(base + "_2l",new TH1D((base+"_2l").c_str(),(base+"_2l").c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l_g1b",new TH1D((base+"_0l_g1b").c_str(),(base+"_0l_g1b").c_str(),17,0,17));
        my_histos.emplace(base + "_1l_g1b",new TH1D((base+"_1l_g1b").c_str(),(base+"_1l_g1b").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_g1b",new TH1D((base+"_2l_g1b").c_str(),(base+"_2l_g1b").c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l_g1b_ht500",new TH1D((base+"_0l_g1b_ht500").c_str(),(base+"_0l_g1b_ht500").c_str(),17,0,17));
        my_histos.emplace(base + "_1l_g1b_ht500",new TH1D((base+"_1l_g1b_ht500").c_str(),(base+"_1l_g1b_ht500").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_g1b_ht500",new TH1D((base+"_2l_g1b_ht500").c_str(),(base+"_2l_g1b_ht500").c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l_g1b_g1t",new TH1D((base+"_0l_g1b_g1t").c_str(),(base+"_0l_g1b_g1t").c_str(),17,0,17));
        my_histos.emplace(base + "_1l_g1b_g1t",new TH1D((base+"_1l_g1b_g1t").c_str(),(base+"_1l_g1b_g1t").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_g1b_g1t",new TH1D((base+"_2l_g1b_g1t").c_str(),(base+"_2l_g1b_g1t").c_str(),15,0,15));
        
        my_histos.emplace(base + "_0l_g1b_g1t_ht500",new TH1D((base+"_0l_g1b_g1t_ht500").c_str(),(base+"_0l_g1b_g1t_ht500").c_str(),17,0,17));
        my_histos.emplace(base + "_1l_g1b_g1t_ht500",new TH1D((base+"_1l_g1b_g1t_ht500").c_str(),(base+"_1l_g1b_g1t_ht500").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_g1b_g1t_ht500",new TH1D((base+"_2l_g1b_g1t_ht500").c_str(),(base+"_2l_g1b_g1t_ht500").c_str(),15,0,15));

        my_histos.emplace(base + "_0l_g1b_2t_ht500",new TH1D((base+"_0l_g1b_2t_ht500").c_str(),(base+"_0l_g1b_2t_ht500").c_str(),17,0,17));

        my_histos.emplace(base + "_1l_g1b_mbl",new TH1D((base+"_1l_g1b_mbl").c_str(),(base+"_1l_g1b_mbl").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt1",new TH1D((base+"_1l_g1b_mbl_bdt1").c_str(),(base+"_1l_g1b_mbl_bdt1").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt2",new TH1D((base+"_1l_g1b_mbl_bdt2").c_str(),(base+"_1l_g1b_mbl_bdt2").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt3",new TH1D((base+"_1l_g1b_mbl_bdt3").c_str(),(base+"_1l_g1b_mbl_bdt3").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt4",new TH1D((base+"_1l_g1b_mbl_bdt4").c_str(),(base+"_1l_g1b_mbl_bdt4").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher1",new TH1D((base+"_1l_g1b_mbl_fisher1").c_str(),(base+"_1l_g1b_mbl_fisher1").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher2",new TH1D((base+"_1l_g1b_mbl_fisher2").c_str(),(base+"_1l_g1b_mbl_fisher2").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher3",new TH1D((base+"_1l_g1b_mbl_fisher3").c_str(),(base+"_1l_g1b_mbl_fisher3").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher4",new TH1D((base+"_1l_g1b_mbl_fisher4").c_str(),(base+"_1l_g1b_mbl_fisher4").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_g1t",new TH1D((base+"_1l_g1b_mbl_g1t").c_str(),(base+"_1l_g1b_mbl_g1t").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_0t",new TH1D((base+"_1l_g1b_mbl_0t").c_str(),(base+"_1l_g1b_mbl_0t").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_1t1",new TH1D((base+"_1l_g1b_mbl_1t1").c_str(),(base+"_1l_g1b_mbl_1t1").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_1t2",new TH1D((base+"_1l_g1b_mbl_1t2").c_str(),(base+"_1l_g1b_mbl_1t2").c_str(),15,0,15));
        my_histos.emplace(base + "_1l_g1b_mbl_1t3",new TH1D((base+"_1l_g1b_mbl_1t3").c_str(),(base+"_1l_g1b_mbl_1t3").c_str(),15,0,15));
        
        // For Z->ll control region
        my_histos.emplace(base + "_2l_onZ",new TH1D((base+"_2l_onZ").c_str(),(base+"_2l_onZ").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b",new TH1D((base+"_2l_onZ_g1b").c_str(),(base+"_2l_onZ_g1b").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl",new TH1D((base+"_2l_onZ_g1b_nombl").c_str(),(base+"_2l_onZ_g1b_nombl").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_g1t",new TH1D((base+"_2l_onZ_g1b_g1t").c_str(),(base+"_2l_onZ_g1b_g1t").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt1",new TH1D((base+"_2l_onZ_g1b_nombl_bdt1").c_str(),(base+"_2l_onZ_g1b_nombl_bdt1").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt2",new TH1D((base+"_2l_onZ_g1b_nombl_bdt2").c_str(),(base+"_2l_onZ_g1b_nombl_bdt2").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt3",new TH1D((base+"_2l_onZ_g1b_nombl_bdt3").c_str(),(base+"_2l_onZ_g1b_nombl_bdt3").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt4",new TH1D((base+"_2l_onZ_g1b_nombl_bdt4").c_str(),(base+"_2l_onZ_g1b_nombl_bdt4").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher1",new TH1D((base+"_2l_onZ_g1b_fisher1").c_str(),(base+"_2l_onZ_g1b_fisher1").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher2",new TH1D((base+"_2l_onZ_g1b_fisher2").c_str(),(base+"_2l_onZ_g1b_fisher2").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher3",new TH1D((base+"_2l_onZ_g1b_fisher3").c_str(),(base+"_2l_onZ_g1b_fisher3").c_str(),15,0,15));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher4",new TH1D((base+"_2l_onZ_g1b_fisher4").c_str(),(base+"_2l_onZ_g1b_fisher4").c_str(),15,0,15));
    }
}

void AnalyzeBackground::Loop(NTupleReader& tr, double weight, int maxevents, std::string filetag, bool isQuiet)
{
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
    ReadFisher_350to650_fwm10_jmtev_top6 read_fisher_350to650_fwm10_jmtev_top6( inputVarNames_top6 ) ;

    std::vector<std::string> inputVarNames_top6_fwm6 ;
    std::vector<double> bdtInputVals_top6_fwm6 ;

    {
        std::string vname ;
        vname = "fwm2_top6" ; inputVarNames_top6_fwm6.push_back( vname ) ;
        vname = "fwm3_top6" ; inputVarNames_top6_fwm6.push_back( vname ) ;
        vname = "fwm4_top6" ; inputVarNames_top6_fwm6.push_back( vname ) ;
        vname = "fwm5_top6" ; inputVarNames_top6_fwm6.push_back( vname ) ;
        vname = "fwm6_top6" ; inputVarNames_top6_fwm6.push_back( vname ) ;
        vname = "jmt_ev0_top6" ; inputVarNames_top6_fwm6.push_back( vname ) ;
        vname = "jmt_ev1_top6" ; inputVarNames_top6_fwm6.push_back( vname ) ;
        vname = "jmt_ev2_top6" ; inputVarNames_top6_fwm6.push_back( vname ) ;

        for ( unsigned int i=0; i < inputVarNames_top6_fwm6.size() ; i++ ) {
            bdtInputVals_top6_fwm6.push_back( 0.5 ) ; //--- load vector with dummy values.
        } // i

    }

    ReadFisherG_350to650_fwm6_jmtev_top6_gt_v2 read_fisher_350to650_fwm6_jmtev_top6_gt_v2( inputVarNames_top6_fwm6 ) ;

    while(tr.getNextEvent())
    {
        const double& madHT   = tr.getVar<double>("madHT");
        const double& Weight  = tr.getVar<double>("Weight");
        const int& ntops_3jet = tr.getVar<int>("ntops_3jet");
        const int& ntops_2jet = tr.getVar<int>("ntops_2jet");
        const int& ntops_1jet = tr.getVar<int>("ntops_1jet");
        const std::string& runtype = tr.getVar<std::string>("runtype");
        const std::vector<TLorentzVector>& Jets         = tr.getVec<TLorentzVector>("Jets");
        const std::vector<double>&      Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const std::vector<std::string>& TriggerNames           = tr.getVec<std::string>("TriggerNames");
        const std::vector<int>&         TriggerPass            = tr.getVec<int>("TriggerPass");
        const TopTaggerResults* ttr         = tr.getVar<TopTaggerResults*>("ttr");
        const std::vector<TopObject*>& tops = ttr->getTops(); 

        const std::vector<TLorentzVector>& BJets = tr.getVec<TLorentzVector>("BJets");
        const int& NBJets = tr.getVar<int>("NBJets");
        const std::vector<TLorentzVector>& BJets_pt45 = tr.getVec<TLorentzVector>("BJets_pt45");
        const int& NBJets_pt45 = tr.getVar<int>("NBJets_pt45");
        const std::vector<TLorentzVector>& GoodMuons = tr.getVec<TLorentzVector>("GoodMuons");
        const std::vector<int>& GoodMuonsCharge      = tr.getVec<int>("GoodMuonsCharge");
        const std::vector<TLorentzVector>& GoodElectrons = tr.getVec<TLorentzVector>("GoodElectrons");
        const std::vector<int>& GoodElectronsCharge      = tr.getVec<int>("GoodElectronsCharge");
        const std::vector<TLorentzVector>& GoodLeptons = tr.getVec<TLorentzVector>("GoodLeptons");
        const int& NGoodLeptons = tr.getVar<int>("NGoodLeptons");
        const double& HT_trigger = tr.getVar<double>("HT_trigger");
        const double& Mbl        = tr.getVar<double>("Mbl");

        if (maxevents > 0 && tr.getEvtNum() >= maxevents) break;
        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;
        
        // Exclude events with MadGraph HT > 100 from the DY inclusive sample
        if(filetag == "DYJetsToLL_M-50_Incl" && madHT > 100) continue;

        // Make sure event weight is not 0 for data
        double eventweight = 1.;
        if(runtype != "Data")
            eventweight = Weight;
        
        // ------------------------------
        // -- Trigger for data
        // ------------------------------
      
        bool passTriggerAllHad = PassTriggerAllHad(TriggerNames, TriggerPass);
        bool passTriggerMuon = PassTriggerMuon(TriggerNames, TriggerPass);
        bool passTriggerElectron = PassTriggerElectron(TriggerNames, TriggerPass);
        if (runtype == "Data")
        {
            if (filetag == "Data_JetHT")
            {
                if(!passTriggerAllHad) continue;
                else
                {
                    passTriggerMuon = false;
                    passTriggerElectron = false;
                }
            }
            if (filetag == "Data_SingleMuon")
            { 
                if(!passTriggerMuon) continue;
                else
                {
                    passTriggerAllHad = false;
                    passTriggerElectron = false;
                }
            }
            if (filetag == "Data_SingleElectron")
            {
                if(!passTriggerElectron) continue;
                else
                {
                    passTriggerAllHad = false;
                    passTriggerMuon = false;
                }
            }
        }

        // -------------------------------
        // -- Event shape BDT
        // -------------------------------

        std::vector<math::RThetaPhiVector> cm_frame_jets ;
        get_cmframe_jets( &Jets, cm_frame_jets, 6 ) ;
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
        double fisher_val = read_fisher_350to650_fwm10_jmtev_top6.GetMvaValue( bdtInputVals_top6 ) ;

        {
            int vi(0) ;
            bdtInputVals_top6_fwm6.at(vi) = esv_top6.getFWmoment(2) ; vi++ ;
            bdtInputVals_top6_fwm6.at(vi) = esv_top6.getFWmoment(3) ; vi++ ;
            bdtInputVals_top6_fwm6.at(vi) = esv_top6.getFWmoment(4) ; vi++ ;
            bdtInputVals_top6_fwm6.at(vi) = esv_top6.getFWmoment(5) ; vi++ ;
            bdtInputVals_top6_fwm6.at(vi) = esv_top6.getFWmoment(6) ; vi++ ;
            bdtInputVals_top6_fwm6.at(vi) = eigen_vals_norm_top6[0] ; vi++ ;
            bdtInputVals_top6_fwm6.at(vi) = eigen_vals_norm_top6[1] ; vi++ ;
            bdtInputVals_top6_fwm6.at(vi) = eigen_vals_norm_top6[2] ; vi++ ;
        }

        double fisher_val_v2 = read_fisher_350to650_fwm6_jmtev_top6_gt_v2.GetMvaValue( bdtInputVals_top6_fwm6 ) ;



        bool bdt_bin1 = eventshape_bdt_val > -1.   && eventshape_bdt_val <= -0.04;
        bool bdt_bin2 = eventshape_bdt_val > -0.04 && eventshape_bdt_val <= 0;
        bool bdt_bin3 = eventshape_bdt_val > 0     && eventshape_bdt_val <= 0.04;
        bool bdt_bin4 = eventshape_bdt_val > 0.04  && eventshape_bdt_val <= 1;

        bool fisher_bin1 = fisher_val_v2 > -1.    && fisher_val_v2 <= -0.035;
        bool fisher_bin2 = fisher_val_v2 > -0.035 && fisher_val_v2 <= 0.03;
        bool fisher_bin3 = fisher_val_v2 > 0.03   && fisher_val_v2 <= 0.095;
        bool fisher_bin4 = fisher_val_v2 > 0.095  && fisher_val_v2 <= 1;

        // -------------------------------
        // -- Basic event selection stuff
        // -------------------------------

        // Check whether event would pass the trigger requirement
        bool passTrigger = true;
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
        if ( !( HT_trigger>500 && rec_njet_pt45>=6 ) ) 
        {
            passTrigger = false;
        }

        bool passTrigger0l = false;
        bool passTrigger1l = false;
        bool passTrigger2l = false;
        if(runtype == "Data")
        {
            if (GoodMuons.size() > 0)
            {
                passTrigger1l = passTriggerMuon && (filetag == "Data_SingleMuon");
                passTrigger2l = passTriggerMuon && (filetag == "Data_SingleMuon");
            } 
            else if (GoodElectrons.size() > 0)
            {
                passTrigger1l = passTriggerElectron && (filetag == "Data_SingleElectron");
                passTrigger2l = passTriggerElectron && (filetag == "Data_SingleElectron");
            }
            else
            {
                passTrigger0l = passTriggerAllHad && (filetag == "Data_JetHT");
            }
        } 
        else 
        {
            passTrigger0l = true;
            passTrigger1l = true;
            passTrigger2l = true;
        }


        //int nleptons = GoodMuons.size() + GoodElectrons.size();
        bool passBaseline0l = NGoodLeptons==0 && rec_njet_pt45>=6 && HT_trigger > 500 && NBJets_pt45 >= 1;
        bool passBaseline1l = NGoodLeptons==1 && rec_njet_pt30>=6 ;
        bool passBaseline2l = NGoodLeptons==2;

        bool passNtop = tops.size() >= 1;
        bool passNb = NBJets >= 1;
        bool onZ = false;
        bool passMbl_2l = false;
        if ( (GoodMuons.size() == 2) && (GoodMuonsCharge[0] != GoodMuonsCharge[1]) )
        {
            double mll = (GoodMuons[0] + GoodMuons[1]).M();
            if( mll > 81.2 && mll < 101.2)
                onZ = true;          

            // check whether a bl pair passes the M(b,l) cut
            for (TLorentzVector myb : BJets)
            {
                double mass_bl_1 = (GoodMuons[0] + myb).M();
                if(mass_bl_1 < 180 && mass_bl_1 > 30)
                    passMbl_2l = true;
                double mass_bl_2 = (GoodMuons[1] + myb).M();
                if(mass_bl_2 < 180 && mass_bl_2 > 30)
                    passMbl_2l = true;
            }
        } 
        else if ( (GoodElectrons.size() == 2) && (GoodElectronsCharge[0] != GoodElectronsCharge[1]) )
        {
            double mll = (GoodElectrons[0] + GoodElectrons[1]).M();
            if( mll > 81.2 && mll < 101.2)
                onZ = true;  
            // check whether a bl pair passes the M(b,l) cut
            for (TLorentzVector myb : BJets)
            {
                double mass_bl_1 = (GoodElectrons[0] + myb).M();
                if(mass_bl_1 < 180 && mass_bl_1 > 30)
                    passMbl_2l = true;
                double mass_bl_2 = (GoodElectrons[1] + myb).M();
                if(mass_bl_2 < 180 && mass_bl_2 > 30)
                    passMbl_2l = true;
            }
        }
        
        // -------------------------
        // -- Check r(j) behavior --
        // -------------------------

        int njets_rj = 0;
        std::vector<std::string> jettypes {"pt30", "pt45"};
        for(std::string jettype : jettypes)
        {
            if(jettype == "pt30")
                njets_rj = rec_njet_pt30;
            else if(jettype == "pt45")
                njets_rj = rec_njet_pt45;
            std::string base = "h_njets_" + jettype;

            my_histos[base]->Fill(njets_rj, eventweight);
            if(NGoodLeptons == 0 && passTrigger0l)
                my_histos[base+"_0l"]->Fill(njets_rj, eventweight);
            else if(NGoodLeptons == 1 && passTrigger1l)
                my_histos[base+"_1l"]->Fill(njets_rj, eventweight);
            else if(NGoodLeptons == 2 && passTrigger2l)
                my_histos[base+"_2l"]->Fill(njets_rj, eventweight);
            if (passNb)
            {
                if(NGoodLeptons == 0 && passTrigger0l)
                    my_histos[base+"_0l_g1b"]->Fill(njets_rj, eventweight);
                else if(NGoodLeptons == 1 && passTrigger1l)
                    my_histos[base+"_1l_g1b"]->Fill(njets_rj, eventweight);
                else if(NGoodLeptons == 2 && passTrigger2l)
                    my_histos[base+"_2l_g1b"]->Fill(njets_rj, eventweight);
              
                if (HT_trigger > 500)
                {
                    if(NGoodLeptons == 0 && passTrigger0l)
                        my_histos[base+"_0l_g1b_ht500"]->Fill(njets_rj, eventweight);
                    else if(NGoodLeptons == 1 && passTrigger1l)
                        my_histos[base+"_1l_g1b_ht500"]->Fill(njets_rj, eventweight);
                    else if(NGoodLeptons == 2 && passTrigger2l)
                        my_histos[base+"_2l_g1b_ht500"]->Fill(njets_rj, eventweight);
                }
              
                if (passNtop)
                {
                    if(NGoodLeptons == 0 && passTrigger0l)
                        my_histos[base+"_0l_g1b_g1t"]->Fill(njets_rj, eventweight);
                    else if(NGoodLeptons == 1 && passTrigger1l)
                        my_histos[base+"_1l_g1b_g1t"]->Fill(njets_rj, eventweight);
                    else if(NGoodLeptons == 2 && passTrigger2l)
                        my_histos[base+"_2l_g1b_g1t"]->Fill(njets_rj, eventweight);
                  
                    if (HT_trigger > 500)
                    {
                        if(NGoodLeptons == 0 && passTrigger0l)
                        {
                            my_histos[base+"_0l_g1b_g1t_ht500"]->Fill(njets_rj, eventweight);
                            if(tops.size() == 2)
                                my_histos[base+"_0l_g1b_2t_ht500"]->Fill(njets_rj, eventweight);
                        }
                        else if(NGoodLeptons == 1 && passTrigger1l)
                            my_histos[base+"_1l_g1b_g1t_ht500"]->Fill(njets_rj, eventweight);
                        else if(NGoodLeptons == 2 && passTrigger2l)
                            my_histos[base+"_2l_g1b_g1t_ht500"]->Fill(njets_rj, eventweight);
                    }
                }
            }

            // Dedicated 1l region
            if(NGoodLeptons == 1 && passNb && passTrigger1l)
            {
                bool passMtop = Mbl < 180 && Mbl > 30;
                
                if(passMtop)
                {
                    my_histos[base+"_1l_g1b_mbl"]->Fill(njets_rj, eventweight);

                    if(bdt_bin1)
                        my_histos[base+"_1l_g1b_mbl_bdt1"]->Fill(njets_rj, eventweight);
                    if(bdt_bin2)
                        my_histos[base+"_1l_g1b_mbl_bdt2"]->Fill(njets_rj, eventweight);
                    if(bdt_bin3)
                        my_histos[base+"_1l_g1b_mbl_bdt3"]->Fill(njets_rj, eventweight);
                    if(bdt_bin4)
                        my_histos[base+"_1l_g1b_mbl_bdt4"]->Fill(njets_rj, eventweight);

                    if(fisher_bin1)
                        my_histos[base+"_1l_g1b_mbl_fisher1"]->Fill(njets_rj, eventweight);
                    if(fisher_bin2)
                        my_histos[base+"_1l_g1b_mbl_fisher2"]->Fill(njets_rj, eventweight);
                    if(fisher_bin3)
                        my_histos[base+"_1l_g1b_mbl_fisher3"]->Fill(njets_rj, eventweight);
                    if(fisher_bin4)
                        my_histos[base+"_1l_g1b_mbl_fisher4"]->Fill(njets_rj, eventweight);

                    if(passNtop)
                    {
                        my_histos[base+"_1l_g1b_mbl_g1t"]->Fill(njets_rj, eventweight);

                        if(tops.size() == 1 && ntops_1jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t1"]->Fill(njets_rj, eventweight);
                        if(tops.size() == 1 && ntops_2jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t2"]->Fill(njets_rj, eventweight);
                        if(tops.size() == 1 && ntops_3jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t3"]->Fill(njets_rj, eventweight);
                    }
                    else 
                    {
                        my_histos[base+"_1l_g1b_mbl_0t"]->Fill(njets_rj, eventweight);
                    }
                }
            }
          
            // Now for the Z->ll region
            // h_njets_2l_onZ_g1b_g1t
            if (onZ && passTrigger2l)
            {
                my_histos[base+"_2l_onZ"]->Fill(njets_rj, eventweight);
                if(passNb)
                {
                    my_histos[base+"_2l_onZ_g1b"]->Fill(njets_rj, eventweight);

                    if(fisher_bin1)
                        my_histos[base+"_2l_onZ_g1b_fisher1"]->Fill(njets_rj, eventweight);
                    if(fisher_bin2)
                        my_histos[base+"_2l_onZ_g1b_fisher2"]->Fill(njets_rj, eventweight);
                    if(fisher_bin3)
                        my_histos[base+"_2l_onZ_g1b_fisher3"]->Fill(njets_rj, eventweight);
                    if(fisher_bin4)
                        my_histos[base+"_2l_onZ_g1b_fisher4"]->Fill(njets_rj, eventweight);
                  
                    if(!passMbl_2l)
                    {
                        my_histos[base+"_2l_onZ_g1b_nombl"]->Fill(njets_rj, eventweight);
                        if(bdt_bin1)
                            my_histos[base+"_2l_onZ_g1b_nombl_bdt1"]->Fill(njets_rj, eventweight);
                        if(bdt_bin2)
                            my_histos[base+"_2l_onZ_g1b_nombl_bdt2"]->Fill(njets_rj, eventweight);
                        if(bdt_bin3)
                            my_histos[base+"_2l_onZ_g1b_nombl_bdt3"]->Fill(njets_rj, eventweight);
                        if(bdt_bin4)
                            my_histos[base+"_2l_onZ_g1b_nombl_bdt4"]->Fill(njets_rj, eventweight);
                    }
                  
                    if(passNtop)
                    {
                        my_histos[base+"_2l_onZ_g1b_g1t"]->Fill(njets_rj, eventweight);
                    }
                }
            }
        }
    }
}


bool AnalyzeBackground::PassTriggerGeneral(std::vector<std::string> &mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
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

bool AnalyzeBackground::PassTriggerAllHad(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
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

bool AnalyzeBackground::PassTriggerMuon(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {"HLT_IsoMu24","HLT_IsoTkMu24_v"};
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

bool AnalyzeBackground::PassTriggerElectron(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
{
    std::vector<std::string> mytriggers {"HLT_Ele27_WPTight_Gsf"};
    return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
}

void AnalyzeBackground::WriteHistos()
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

//  LocalWords:  GoodMuonsCharge
