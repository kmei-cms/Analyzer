#define AnalyzeEventSelection_cxx
#include "Analyzer/Analyzer/include/AnalyzeEventSelection.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

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

AnalyzeEventSelection::AnalyzeEventSelection()
{
    InitHistos();
}

void AnalyzeEventSelection::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("h_met", std::make_shared<TH1D>("h_met","h_met", 20, 0, 200));
    my_histos.emplace("h_ht", std::make_shared<TH1D>("h_ht","h_ht", 60, 0, 3000));
    my_histos.emplace("h_ntops", std::make_shared<TH1D>("h_ntops","h_ntops", 5, 0, 5));
    
    // 0 lepton plots
    // "6j" is the control region only. Only look at data for that region
    // Add cuts on Owen's BDT 
    std::vector<std::string> mycuts_0l {"g6j_HT500_g1b", "g6j_HT500_g1b_0t", "g6j_HT500_g1b_1t", "g6j_HT500_g1b_2t",
            "6j_HT500_g1b_0t", "6j_HT500_g1b_1t1", "6j_HT500_g1b_1t2","6j_HT500_g1b_1t3","6j_HT500_g1b_2t", 
            "g7j_HT500_g1b_0t", "g7j_HT500_g1b_1t1", "g7j_HT500_g1b_1t2", "g7j_HT500_g1b_1t3", "g7j_HT500_g1b_2t"};
    for(std::string mycut : mycuts_0l)
    {
        my_histos.emplace("h_njets_0l_"+mycut, std::make_shared<TH1D>(("h_njets_0l_"+mycut).c_str(),("h_njets_0l_"+mycut).c_str(), 19, 0, 19));
        my_histos.emplace("h_ntops_0l_"+mycut, std::make_shared<TH1D>(("h_ntops_0l_"+mycut).c_str(),("h_ntops_0l_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_0l_"+mycut, std::make_shared<TH1D>(("h_nb_0l_"+mycut).c_str(),("h_nb_0l_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_0l_"+mycut, std::make_shared<TH1D>(("h_HT_0l_"+mycut).c_str(),("h_HT_0l_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_bdt_0l_"+mycut, std::make_shared<TH1D>(("h_bdt_0l_"+mycut).c_str(),("h_bdt_0l_"+mycut).c_str(), 40, -0.5, 0.5));
        my_histos.emplace("h_fisher_0l_"+mycut, std::make_shared<TH1D>(("h_fisher_0l_"+mycut).c_str(),("h_fisher_0l_"+mycut).c_str(), 50, -0.25, 0.25));

        my_2d_histos.emplace("h_njets_bdt_0l_"+mycut, std::make_shared<TH2D>(("h_njets_bdt_0l_"+mycut).c_str(),("h_njets_bdt_0l_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_fisher_0l_"+mycut, std::make_shared<TH2D>(("h_njets_fisher_0l_"+mycut).c_str(),("h_njets_fisher_0l_"+mycut).c_str(), 15, 0, 15, 50, -0.25, 0.25));
    }

    // 1 lepton plots
    // attempt to have "6j" be the control region, also include 6j in case that works better
    // TODO: add BDT cuts
    std::vector<std::string> mycuts_1l {"g6j_g1b", "g6j_g1b_mbl", "g6j_g1b_mbl_0t", "g6j_g1b_mbl_1t1", "g6j_g1b_mbl_1t2", "g6j_g1b_mbl_1t3",
            "6j_g1b", "6j_g1b_mbl", "6j_g1b_mbl_0t", "6j_g1b_mbl_1t1", "6j_g1b_mbl_1t2", "6j_g1b_mbl_1t3",
            "g7j_g1b", "g7j_g1b_mbl", "g7j_g1b_mbl_0t", "g7j_g1b_mbl_1t1", "g7j_g1b_mbl_1t2", "g7j_g1b_mbl_1t3",
            };
    std::vector<std::string> mycuts_1mu {
        "6j_g1b", "6j_g1b_mbl", "6j_g1b_mbl_0t", "6j_g1b_mbl_1t1", "6j_g1b_mbl_1t2", "6j_g1b_mbl_1t3",
            };
    std::vector<std::string> mycuts_1el {
        "6j_g1b", "6j_g1b_mbl", "6j_g1b_mbl_0t", "6j_g1b_mbl_1t1", "6j_g1b_mbl_1t2", "6j_g1b_mbl_1t3",
            };
    for(std::string mycut : mycuts_1l)
    {
        my_histos.emplace("h_njets_1l_"+mycut, std::make_shared<TH1D>(("h_njets_1l_"+mycut).c_str(),("h_njets_1l_"+mycut).c_str(), 15, 0, 15));
        my_histos.emplace("h_ntops_1l_"+mycut, std::make_shared<TH1D>(("h_ntops_1l_"+mycut).c_str(),("h_ntops_1l_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_1l_"+mycut, std::make_shared<TH1D>(("h_nb_1l_"+mycut).c_str(),("h_nb_1l_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_1l_"+mycut, std::make_shared<TH1D>(("h_HT_1l_"+mycut).c_str(),("h_HT_1l_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_mbl_1l_"+mycut, std::make_shared<TH1D>(("h_mbl_1l_"+mycut).c_str(),("h_mbl_1l_"+mycut).c_str(), 30, 0, 300));
        my_histos.emplace("h_bdt_1l_"+mycut, std::make_shared<TH1D>(("h_bdt_1l_"+mycut).c_str(),("h_bdt_1l_"+mycut).c_str(), 40, -0.5, 0.5));
        my_histos.emplace("h_fisher_1l_"+mycut, std::make_shared<TH1D>(("h_fisher_1l_"+mycut).c_str(),("h_fisher_1l_"+mycut).c_str(), 50, -0.25, 0.25));

        my_2d_histos.emplace("h_njets_bdt_1l_"+mycut, std::make_shared<TH2D>(("h_njets_bdt_1l_"+mycut).c_str(),("h_njets_bdt_1l_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_fisher_1l_"+mycut, std::make_shared<TH2D>(("h_njets_fisher_1l_"+mycut).c_str(),("h_njets_fisher_1l_"+mycut).c_str(), 15, 0, 15, 50, -0.25, 0.25));
    }
    for(std::string mycut : mycuts_1mu)
    {
        my_histos.emplace("h_mupt_1mu_"+mycut, std::make_shared<TH1D>(("h_mupt_1mu_"+mycut).c_str(),("h_mupt_1mu_"+mycut).c_str(), 50, 0, 500));
        my_histos.emplace("h_njets_1mu_"+mycut, std::make_shared<TH1D>(("h_njets_1mu_"+mycut).c_str(),("h_njets_1mu_"+mycut).c_str(), 15, 0, 15));
        my_histos.emplace("h_ntops_1mu_"+mycut, std::make_shared<TH1D>(("h_ntops_1mu_"+mycut).c_str(),("h_ntops_1mu_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_1mu_"+mycut, std::make_shared<TH1D>(("h_nb_1mu_"+mycut).c_str(),("h_nb_1mu_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_1mu_"+mycut, std::make_shared<TH1D>(("h_HT_1mu_"+mycut).c_str(),("h_HT_1mu_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_mbl_1mu_"+mycut, std::make_shared<TH1D>(("h_mbl_1mu_"+mycut).c_str(),("h_mbl_1mu_"+mycut).c_str(), 30, 0, 300));
        my_histos.emplace("h_bdt_1mu_"+mycut, std::make_shared<TH1D>(("h_bdt_1mu_"+mycut).c_str(),("h_bdt_1mu_"+mycut).c_str(), 40, -0.5, 0.5));
        my_histos.emplace("h_fisher_1mu_"+mycut, std::make_shared<TH1D>(("h_fisher_1mu_"+mycut).c_str(),("h_fisher_1mu_"+mycut).c_str(), 50, -0.25, 0.25));

        my_2d_histos.emplace("h_njets_bdt_1mu_"+mycut, std::make_shared<TH2D>(("h_njets_bdt_1mu_"+mycut).c_str(),("h_njets_bdt_1mu_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_fisher_1mu_"+mycut, std::make_shared<TH2D>(("h_njets_fisher_1mu_"+mycut).c_str(),("h_njets_fisher_1mu_"+mycut).c_str(), 15, 0, 15, 50, -0.25, 0.25));
    }
    for(std::string mycut : mycuts_1el)
    {
        my_histos.emplace("h_elpt_1el_"+mycut, std::make_shared<TH1D>(("h_elpt_1el_"+mycut).c_str(),("h_elpt_1el_"+mycut).c_str(), 50, 0, 500));
        my_histos.emplace("h_njets_1el_"+mycut, std::make_shared<TH1D>(("h_njets_1el_"+mycut).c_str(),("h_njets_1el_"+mycut).c_str(), 15, 0, 15));
        my_histos.emplace("h_ntops_1el_"+mycut, std::make_shared<TH1D>(("h_ntops_1el_"+mycut).c_str(),("h_ntops_1el_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_1el_"+mycut, std::make_shared<TH1D>(("h_nb_1el_"+mycut).c_str(),("h_nb_1el_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_1el_"+mycut, std::make_shared<TH1D>(("h_HT_1el_"+mycut).c_str(),("h_HT_1el_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_mbl_1el_"+mycut, std::make_shared<TH1D>(("h_mbl_1el_"+mycut).c_str(),("h_mbl_1el_"+mycut).c_str(), 30, 0, 300));
        my_histos.emplace("h_bdt_1el_"+mycut, std::make_shared<TH1D>(("h_bdt_1el_"+mycut).c_str(),("h_bdt_1el_"+mycut).c_str(), 40, -0.5, 0.5));
        my_histos.emplace("h_fisher_1el_"+mycut, std::make_shared<TH1D>(("h_fisher_1el_"+mycut).c_str(),("h_fisher_1el_"+mycut).c_str(), 50, -0.25, 0.25));

        my_2d_histos.emplace("h_njets_bdt_1el_"+mycut, std::make_shared<TH2D>(("h_njets_bdt_1el_"+mycut).c_str(),("h_njets_bdt_1el_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_fisher_1el_"+mycut, std::make_shared<TH2D>(("h_njets_fisher_1el_"+mycut).c_str(),("h_njets_fisher_1el_"+mycut).c_str(), 15, 0, 15, 50, -0.25, 0.25));
    }

    
    // 2 lepton plots
    // onZ control region only for now, also check for 1top to see if a fake top changes any behavior
    std::vector<std::string> mycuts_2l {"onZ_g1b", "onZ_g1b_g1t", "2b",
            "onZ_g1b_bdt1","onZ_g1b_bdt2","onZ_g1b_bdt3","onZ_g1b_bdt4"};
    for(std::string mycut : mycuts_2l)
    {
        my_histos.emplace("h_njets_2l_"+mycut, std::make_shared<TH1D>(("h_njets_2l_"+mycut).c_str(),("h_njets_2l_"+mycut).c_str(), 15, 0, 15));
        my_histos.emplace("h_ntops_2l_"+mycut, std::make_shared<TH1D>(("h_ntops_2l_"+mycut).c_str(),("h_ntops_2l_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_2l_"+mycut, std::make_shared<TH1D>(("h_nb_2l_"+mycut).c_str(),("h_nb_2l_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_2l_"+mycut, std::make_shared<TH1D>(("h_HT_2l_"+mycut).c_str(),("h_HT_2l_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_bdt_2l_"+mycut, std::make_shared<TH1D>(("h_bdt_2l_"+mycut).c_str(),("h_bdt_2l_"+mycut).c_str(), 40, -0.5, 0.5));
        my_histos.emplace("h_fisher_2l_"+mycut, std::make_shared<TH1D>(("h_fisher_2l_"+mycut).c_str(),("h_fisher_2l_"+mycut).c_str(), 50, -0.25, 0.25));

        my_2d_histos.emplace("h_njets_bdt_2l_"+mycut, std::make_shared<TH2D>(("h_njets_bdt_2l_"+mycut).c_str(),("h_njets_bdt_2l_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
        my_2d_histos.emplace("h_njets_fisher_2l_"+mycut, std::make_shared<TH2D>(("h_njets_fisher_2l_"+mycut).c_str(),("h_njets_fisher_2l_"+mycut).c_str(), 15, 0, 15, 50, -0.25, 0.25));

    }

    // Cut flows
    my_efficiencies.emplace("event_sel", std::make_shared<TEfficiency>("event_sel","Event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total", std::make_shared<TEfficiency>("event_sel_total","Total event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_0l", std::make_shared<TEfficiency>("event_sel_0l","0 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_0l", std::make_shared<TEfficiency>("event_sel_total_0l","Total 0 lepton event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_1l", std::make_shared<TEfficiency>("event_sel_1l","1 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_1l", std::make_shared<TEfficiency>("event_sel_total_1l","Total 1 lepton event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_2l", std::make_shared<TEfficiency>("event_sel_2l","2 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_2l", std::make_shared<TEfficiency>("event_sel_total_2l","Total 2 lepton event selection efficiency;Cut;#epsilon",8,0,8));

}

void AnalyzeEventSelection::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& MET          = tr.getVar<double>("MET");
        const auto& METPhi       = tr.getVar<double>("METPhi");
        const auto& HT           = tr.getVar<double>("HT");
        const auto& ntops        = tr.getVar<int>("ntops");
        const auto& ntops_3jet   = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet   = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet   = tr.getVar<int>("ntops_1jet");
        const auto& runtype      = tr.getVar<std::string>("runtype");
        const auto& filetag      = tr.getVar<std::string>("filetag");
        const auto& passTrigger  = tr.getVar<bool>("passTrigger");
        const auto* ttr          = tr.getVar<TopTaggerResults*>("ttr");
        const std::vector<TopObject*>& tops = ttr->getTops();

        const auto& NJets_pt30          = tr.getVar<int>("NJets_pt30");
        const auto& NJets_pt45          = tr.getVar<int>("NJets_pt45");
        const auto& BJets_pt30          = tr.getVec<TLorentzVector>("BJets_pt30");
        const auto& NBJets_pt30         = tr.getVar<int>("NBJets_pt30");
        const auto& BJets_pt45          = tr.getVec<TLorentzVector>("BJets_pt45");
        const auto& NBJets_pt45         = tr.getVar<int>("NBJets_pt45");
        const auto& GoodMuons           = tr.getVec<TLorentzVector>("GoodMuons");
        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");
        const auto& GoodMuonsCharge     = tr.getVec<int>("GoodMuonsCharge");
        const auto& GoodElectrons       = tr.getVec<TLorentzVector>("GoodElectrons");
        const auto& NGoodElectrons      = tr.getVar<int>("NGoodElectrons");
        const auto& GoodElectronsCharge = tr.getVec<int>("GoodElectronsCharge");
        const auto& GoodLeptons         = tr.getVec<TLorentzVector>("GoodLeptons");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& HT_trigger          = tr.getVar<double>("HT_trigger");
        const auto& Mbl                 = tr.getVar<double>("Mbl");
        const auto& used_bjet_for_mbl   = tr.getVar<TLorentzVector>("used_bjet_for_mbl");
        
        const auto& fisher_bin1 = tr.getVar<bool>("fisher_bin1");
        const auto& fisher_bin2 = tr.getVar<bool>("fisher_bin2");
        const auto& fisher_bin3 = tr.getVar<bool>("fisher_bin3");
        const auto& fisher_bin4 = tr.getVar<bool>("fisher_bin4");
        const auto& bdt_bin1    = tr.getVar<bool>("bdt_bin1");
        const auto& bdt_bin2    = tr.getVar<bool>("bdt_bin2");
        const auto& bdt_bin3    = tr.getVar<bool>("bdt_bin3");
        const auto& bdt_bin4    = tr.getVar<bool>("bdt_bin4");
        const auto& eventshape_bdt_val = tr.getVar<double>("eventshape_bdt_val");
        const auto& fisher_val         = tr.getVar<double>("fisher_val");


        const auto& passBaseline0l        = tr.getVar<bool>("passBaseline0l");
        const auto& passBaseline1l        = tr.getVar<bool>("passBaseline1l");
        const auto& passBaseline1mu       = tr.getVar<bool>("passBaseline1mu");
        const auto& passBaseline1el       = tr.getVar<bool>("passBaseline1el");
        const auto& passBaseline2l        = tr.getVar<bool>("passBaseline2l");
        const auto& passBaseline2lonZ     = tr.getVar<bool>("passBaseline2lonZ");
       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;
        
        // Exclude events with MadGraph HT > 100 from the DY inclusive sample
        if(filetag == "DYJetsToLL_M-50_Incl")
        {
            const double& madHT   = tr.getVar<double>("madHT");
            if(madHT > 100) continue;
        } 
        
        // Make sure event weight is not 0 for data
        double eventweight = 1.;
        if(runtype != "Data")
            eventweight = tr.getVar<double>("Weight");
        
        
        // ------------------
        // --- TOP TAGGER ---
        // ------------------
        my_histos["h_ntops"]->Fill(ntops, eventweight);
                
        // -------------------------------
        // -- Basic event selection stuff
        // -------------------------------
                
        bool pass_0t = ntops==0, pass_1t = ntops==1, pass_2t = ntops==2;
        bool pass_1t1 = ntops==1 && ntops_1jet==1, pass_1t2 = ntops==1 && ntops_2jet==1, pass_1t3 = ntops==1 && ntops_3jet==1;

        bool pass_mbl = Mbl > 30 && Mbl < 180;
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
                    if(c->p() == used_bjet_for_mbl)
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
                pass_0t = (ntops - ntops_to_remove) == 0;
                pass_1t = (ntops - ntops_to_remove) == 1;
                pass_2t = (ntops - ntops_to_remove) == 2;
                pass_1t1 = (ntops_1jet - top_type1_to_remove) == 1;
                pass_1t2 = (ntops_2jet - top_type2_to_remove) == 1;
                pass_1t3 = (ntops_3jet - top_type3_to_remove) == 1;
            }
            //std::cout << "New top counting: " << pass_0t << " " << pass_1t << " " << pass_2t << " " << pass_1t1 << " " << pass_1t2 << " " << pass_1t3 << std::endl;
        }
        
        const std::map<std::string, bool> cut_map_0l {
            {"g6j_HT500_g1b",     passBaseline0l},
            {"g6j_HT500_g1b_0t",  passBaseline0l && pass_0t},
            {"g6j_HT500_g1b_1t",  passBaseline0l && pass_1t},
            {"g6j_HT500_g1b_2t",  passBaseline0l && pass_2t},
            {"6j_HT500_g1b_0t",   passBaseline0l && NJets_pt30==6 && pass_0t},
            {"6j_HT500_g1b_1t1",  passBaseline0l && NJets_pt30==6 && pass_1t1},
            {"6j_HT500_g1b_1t2",  passBaseline0l && NJets_pt30==6 && pass_1t2},
            {"6j_HT500_g1b_1t3",  passBaseline0l && NJets_pt30==6 && pass_1t3},
            {"6j_HT500_g1b_2t",   passBaseline0l && NJets_pt30==6 && pass_2t},
            {"g7j_HT500_g1b_0t",  passBaseline0l && NJets_pt30>=7 && pass_0t},
            {"g7j_HT500_g1b_1t1", passBaseline0l && NJets_pt30>=7 && pass_1t1},
            {"g7j_HT500_g1b_1t2", passBaseline0l && NJets_pt30>=7 && pass_1t2},
            {"g7j_HT500_g1b_1t3", passBaseline0l && NJets_pt30>=7 && pass_1t3},
            {"g7j_HT500_g1b_2t",  passBaseline0l && NJets_pt30>=7 && pass_2t}
        };
        
        for(auto& kv : cut_map_0l)
        {
            if(kv.second)
            {
                my_histos["h_njets_0l_"+kv.first]->Fill(NJets_pt30, eventweight);
                my_histos["h_ntops_0l_"+kv.first]->Fill(ntops, eventweight);
                my_histos["h_nb_0l_"+kv.first]->Fill(NBJets_pt30, eventweight);
                my_histos["h_HT_0l_"+kv.first]->Fill(HT_trigger, eventweight);
                my_histos["h_bdt_0l_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
                my_histos["h_fisher_0l_"+kv.first]->Fill(fisher_val, eventweight);
                my_2d_histos["h_njets_bdt_0l_"+kv.first]->Fill(NJets_pt30, eventshape_bdt_val, eventweight);
                my_2d_histos["h_njets_fisher_0l_"+kv.first]->Fill(NJets_pt30, fisher_val, eventweight);
            }
        }
        
        const std::map<std::string, bool> cut_map_1l {
            {"g6j_g1b",         passBaseline1l},
            {"g6j_g1b_mbl",     passBaseline1l && pass_mbl},
            {"g6j_g1b_mbl_0t",  passBaseline1l && pass_mbl && pass_0t},
            {"g6j_g1b_mbl_1t1", passBaseline1l && pass_mbl && pass_1t1},
            {"g6j_g1b_mbl_1t2", passBaseline1l && pass_mbl && pass_1t2},
            {"g6j_g1b_mbl_1t3", passBaseline1l && pass_mbl && pass_1t3},
            {"6j_g1b",          passBaseline1l && NJets_pt30==6},
            {"6j_g1b_mbl",      passBaseline1l && NJets_pt30==6 && pass_mbl},
            {"6j_g1b_mbl_0t",   passBaseline1l && NJets_pt30==6 && pass_mbl && pass_0t},
            {"6j_g1b_mbl_1t1",  passBaseline1l && NJets_pt30==6 && pass_mbl && pass_1t1},
            {"6j_g1b_mbl_1t2",  passBaseline1l && NJets_pt30==6 && pass_mbl && pass_1t2},
            {"6j_g1b_mbl_1t3",  passBaseline1l && NJets_pt30==6 && pass_mbl && pass_1t3},
            {"g7j_g1b",         passBaseline1l && NJets_pt30>=7},
            {"g7j_g1b_mbl",     passBaseline1l && NJets_pt30>=7 && pass_mbl},
            {"g7j_g1b_mbl_0t",  passBaseline1l && NJets_pt30>=7 && pass_mbl && pass_0t},
            {"g7j_g1b_mbl_1t1", passBaseline1l && NJets_pt30>=7 && pass_mbl && pass_1t1},
            {"g7j_g1b_mbl_1t2", passBaseline1l && NJets_pt30>=7 && pass_mbl && pass_1t2},
            {"g7j_g1b_mbl_1t3", passBaseline1l && NJets_pt30>=7 && pass_mbl && pass_1t3},
                };
        const std::map<std::string, bool> cut_map_1mu {
            {"6j_g1b",         passBaseline1mu && NJets_pt30==6 },
            {"6j_g1b_mbl",     passBaseline1mu && NJets_pt30==6 && pass_mbl},
            {"6j_g1b_mbl_0t",  passBaseline1mu && NJets_pt30==6 && pass_mbl && pass_0t},
            {"6j_g1b_mbl_1t1", passBaseline1mu && NJets_pt30==6 && pass_mbl && pass_1t1},
            {"6j_g1b_mbl_1t2", passBaseline1mu && NJets_pt30==6 && pass_mbl && pass_1t2},
            {"6j_g1b_mbl_1t3", passBaseline1mu && NJets_pt30==6 && pass_mbl && pass_1t3},
                };
        const std::map<std::string, bool> cut_map_1el {
            {"6j_g1b",         passBaseline1el && NJets_pt30==6},
            {"6j_g1b_mbl",     passBaseline1el && NJets_pt30==6 && pass_mbl},
            {"6j_g1b_mbl_0t",  passBaseline1el && NJets_pt30==6 && pass_mbl && pass_0t},
            {"6j_g1b_mbl_1t1", passBaseline1el && NJets_pt30==6 && pass_mbl && pass_1t1},
            {"6j_g1b_mbl_1t2", passBaseline1el && NJets_pt30==6 && pass_mbl && pass_1t2},
            {"6j_g1b_mbl_1t3", passBaseline1el && NJets_pt30==6 && pass_mbl && pass_1t3},
                };
                
        for(auto& kv : cut_map_1l)
        {
            if(kv.second)
            {
                my_histos["h_njets_1l_"+kv.first]->Fill(NJets_pt30, eventweight);
                my_histos["h_ntops_1l_"+kv.first]->Fill(ntops, eventweight);
                my_histos["h_nb_1l_"+kv.first]->Fill(NBJets_pt30, eventweight);
                my_histos["h_HT_1l_"+kv.first]->Fill(HT_trigger, eventweight);
                my_histos["h_bdt_1l_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
                my_histos["h_fisher_1l_"+kv.first]->Fill(fisher_val, eventweight);
                my_histos["h_mbl_1l_"+kv.first]->Fill(Mbl, eventweight);
                my_2d_histos["h_njets_bdt_1l_"+kv.first]->Fill(NJets_pt30, eventshape_bdt_val, eventweight);
                my_2d_histos["h_njets_fisher_1l_"+kv.first]->Fill(NJets_pt30, fisher_val, eventweight);
            }
        }
        for(auto& kv : cut_map_1mu)
        {
            if(kv.second)
            {
                my_histos["h_mupt_1mu_"+kv.first]->Fill(GoodMuons[0].Pt(), eventweight);
                my_histos["h_njets_1mu_"+kv.first]->Fill(NJets_pt30, eventweight);
                my_histos["h_ntops_1mu_"+kv.first]->Fill(ntops, eventweight);
                my_histos["h_nb_1mu_"+kv.first]->Fill(NBJets_pt30, eventweight);
                my_histos["h_HT_1mu_"+kv.first]->Fill(HT_trigger, eventweight);
                my_histos["h_bdt_1mu_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
                my_histos["h_fisher_1mu_"+kv.first]->Fill(fisher_val, eventweight);
                my_histos["h_mbl_1mu_"+kv.first]->Fill(Mbl, eventweight);
                my_2d_histos["h_njets_bdt_1mu_"+kv.first]->Fill(NJets_pt30, eventshape_bdt_val, eventweight);
                my_2d_histos["h_njets_fisher_1mu_"+kv.first]->Fill(NJets_pt30, fisher_val, eventweight);
            }
        }
        for(auto& kv : cut_map_1el)
        {
            if(kv.second)
            {
                my_histos["h_elpt_1el_"+kv.first]->Fill(GoodElectrons[0].Pt(), eventweight);
                my_histos["h_njets_1el_"+kv.first]->Fill(NJets_pt30, eventweight);
                my_histos["h_ntops_1el_"+kv.first]->Fill(ntops, eventweight);
                my_histos["h_nb_1el_"+kv.first]->Fill(NBJets_pt30, eventweight);
                my_histos["h_HT_1el_"+kv.first]->Fill(HT_trigger, eventweight);
                my_histos["h_bdt_1el_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
                my_histos["h_fisher_1el_"+kv.first]->Fill(fisher_val, eventweight);
                my_histos["h_mbl_1el_"+kv.first]->Fill(Mbl, eventweight);
                my_2d_histos["h_njets_bdt_1el_"+kv.first]->Fill(NJets_pt30, eventshape_bdt_val, eventweight);
                my_2d_histos["h_njets_fisher_1el_"+kv.first]->Fill(NJets_pt30, fisher_val, eventweight);
            }
        }
        
        const std::map<std::string, bool> cut_map_2l {
            {"onZ_g1b",            passBaseline2lonZ},
            {"onZ_g1b_bdt1",       passBaseline2lonZ && bdt_bin1},
            {"onZ_g1b_bdt2",       passBaseline2lonZ && bdt_bin2},
            {"onZ_g1b_bdt3",       passBaseline2lonZ && bdt_bin3},
            {"onZ_g1b_bdt4",       passBaseline2lonZ && bdt_bin4},
            {"onZ_g1b_g1t",        passBaseline2lonZ && pass_1t}, 
            {"2b",                 passBaseline2l && NBJets_pt30 == 2} 
        };
        
        for(auto& kv : cut_map_2l)
        {
            if(kv.second)
            {
                my_histos["h_njets_2l_"+kv.first]->Fill(NJets_pt30, eventweight);
                my_histos["h_ntops_2l_"+kv.first]->Fill(ntops, eventweight);
                my_histos["h_nb_2l_"+kv.first]->Fill(NBJets_pt30, eventweight);
                my_histos["h_HT_2l_"+kv.first]->Fill(HT_trigger, eventweight);
                my_histos["h_bdt_2l_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
                my_histos["h_fisher_2l_"+kv.first]->Fill(fisher_val, eventweight);
                my_2d_histos["h_njets_bdt_2l_"+kv.first]->Fill(NJets_pt30, eventshape_bdt_val, eventweight);
                my_2d_histos["h_njets_fisher_2l_"+kv.first]->Fill(NJets_pt30, fisher_val, eventweight);
            }
        }
        
        // Fill event selection efficiencies
        my_efficiencies["event_sel_total"]->Fill(true,0);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500,1);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 ,2);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 ,3);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 ,4);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 ,5);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 && ntops>1 ,6);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 && ntops>1 && NJets_pt30>=8 ,7);
        
        my_efficiencies["event_sel"]->Fill(true,0);
        my_efficiencies["event_sel"]->Fill(HT_trigger>500,1);
        if(HT_trigger>500)
        {
            my_efficiencies["event_sel"]->Fill(NJets_pt45>=6,2);
            if (NJets_pt45>=6)
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
                                my_efficiencies["event_sel"]->Fill(NJets_pt30>=8,7);
                            }
                        }
                    }
                }
            }
        }
        
        my_histos["h_met"]->Fill(MET, eventweight);
        my_histos["h_ht"]->Fill(HT, eventweight);


    } // end of event loop

}


void AnalyzeEventSelection::WriteHistos(TFile* outfile)
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
