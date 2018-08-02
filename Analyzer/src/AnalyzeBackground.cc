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

AnalyzeBackground::AnalyzeBackground()
{
    InitHistos();
}

void AnalyzeBackground::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    std::vector<std::string> jettypes {"pt30", "pt45"};
    double min_x = 5.5;
    int nbins_0l = 11;
    double max_x_0l = 16.5;
    int nbins_1l = 9;
    double max_x_1l = 14.5;
    for(std::string jettype : jettypes)
    {
        std::string base = "h_njets_" + jettype;
        my_histos.emplace(base,std::make_shared<TH1D>(base.c_str(),base.c_str(),nbins_1l,min_x,max_x_1l));
        
        my_histos.emplace(base + "_0l",std::make_shared<TH1D>( (base+"_0l").c_str(),(base+"_0l").c_str(),nbins_0l,min_x,max_x_0l));
        my_histos.emplace(base + "_1l",std::make_shared<TH1D>((base+"_1l").c_str(),(base+"_1l").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l",std::make_shared<TH1D>((base+"_2l").c_str(),(base+"_2l").c_str(),nbins_1l,min_x,max_x_1l));
        
        my_histos.emplace(base + "_0l_g1b",std::make_shared<TH1D>((base+"_0l_g1b").c_str(),(base+"_0l_g1b").c_str(),nbins_0l,min_x,max_x_0l));
        my_histos.emplace(base + "_1l_g1b",std::make_shared<TH1D>((base+"_1l_g1b").c_str(),(base+"_1l_g1b").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_g1b",std::make_shared<TH1D>((base+"_2l_g1b").c_str(),(base+"_2l_g1b").c_str(),nbins_1l,min_x,max_x_1l));
        
        my_histos.emplace(base + "_0l_g1b_ht500",std::make_shared<TH1D>((base+"_0l_g1b_ht500").c_str(),(base+"_0l_g1b_ht500").c_str(),nbins_0l,min_x,max_x_0l));
        my_histos.emplace(base + "_1l_g1b_ht500",std::make_shared<TH1D>((base+"_1l_g1b_ht500").c_str(),(base+"_1l_g1b_ht500").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_g1b_ht500",std::make_shared<TH1D>((base+"_2l_g1b_ht500").c_str(),(base+"_2l_g1b_ht500").c_str(),nbins_1l,min_x,max_x_1l));
        
        my_histos.emplace(base + "_0l_g1b_g1t",std::make_shared<TH1D>((base+"_0l_g1b_g1t").c_str(),(base+"_0l_g1b_g1t").c_str(),nbins_0l,min_x,max_x_0l));
        my_histos.emplace(base + "_1l_g1b_g1t",std::make_shared<TH1D>((base+"_1l_g1b_g1t").c_str(),(base+"_1l_g1b_g1t").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_g1b_g1t",std::make_shared<TH1D>((base+"_2l_g1b_g1t").c_str(),(base+"_2l_g1b_g1t").c_str(),nbins_1l,min_x,max_x_1l));
        
        my_histos.emplace(base + "_0l_g1b_g1t_ht500",std::make_shared<TH1D>((base+"_0l_g1b_g1t_ht500").c_str(),(base+"_0l_g1b_g1t_ht500").c_str(),nbins_0l,min_x,max_x_0l));
        my_histos.emplace(base + "_1l_g1b_g1t_ht500",std::make_shared<TH1D>((base+"_1l_g1b_g1t_ht500").c_str(),(base+"_1l_g1b_g1t_ht500").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_g1b_g1t_ht500",std::make_shared<TH1D>((base+"_2l_g1b_g1t_ht500").c_str(),(base+"_2l_g1b_g1t_ht500").c_str(),nbins_1l,min_x,max_x_1l));

        my_histos.emplace(base + "_0l_g1b_2t_ht500",std::make_shared<TH1D>((base+"_0l_g1b_2t_ht500").c_str(),(base+"_0l_g1b_2t_ht500").c_str(),nbins_0l,min_x,max_x_0l));

        my_histos.emplace(base + "_1l_g1b_mbl",std::make_shared<TH1D>((base+"_1l_g1b_mbl").c_str(),(base+"_1l_g1b_mbl").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt1",std::make_shared<TH1D>((base+"_1l_g1b_mbl_bdt1").c_str(),(base+"_1l_g1b_mbl_bdt1").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt2",std::make_shared<TH1D>((base+"_1l_g1b_mbl_bdt2").c_str(),(base+"_1l_g1b_mbl_bdt2").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt3",std::make_shared<TH1D>((base+"_1l_g1b_mbl_bdt3").c_str(),(base+"_1l_g1b_mbl_bdt3").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_bdt4",std::make_shared<TH1D>((base+"_1l_g1b_mbl_bdt4").c_str(),(base+"_1l_g1b_mbl_bdt4").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher1",std::make_shared<TH1D>((base+"_1l_g1b_mbl_fisher1").c_str(),(base+"_1l_g1b_mbl_fisher1").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher2",std::make_shared<TH1D>((base+"_1l_g1b_mbl_fisher2").c_str(),(base+"_1l_g1b_mbl_fisher2").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher3",std::make_shared<TH1D>((base+"_1l_g1b_mbl_fisher3").c_str(),(base+"_1l_g1b_mbl_fisher3").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_fisher4",std::make_shared<TH1D>((base+"_1l_g1b_mbl_fisher4").c_str(),(base+"_1l_g1b_mbl_fisher4").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_g1t",std::make_shared<TH1D>((base+"_1l_g1b_mbl_g1t").c_str(),(base+"_1l_g1b_mbl_g1t").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_g1t_fisher1",std::make_shared<TH1D>((base+"_1l_g1b_mbl_g1t_fisher1").c_str(),(base+"_1l_g1b_mbl_g1t_fisher1").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_g1t_fisher2",std::make_shared<TH1D>((base+"_1l_g1b_mbl_g1t_fisher2").c_str(),(base+"_1l_g1b_mbl_g1t_fisher2").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_g1t_fisher3",std::make_shared<TH1D>((base+"_1l_g1b_mbl_g1t_fisher3").c_str(),(base+"_1l_g1b_mbl_g1t_fisher3").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_g1t_fisher4",std::make_shared<TH1D>((base+"_1l_g1b_mbl_g1t_fisher4").c_str(),(base+"_1l_g1b_mbl_g1t_fisher4").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_0t",std::make_shared<TH1D>((base+"_1l_g1b_mbl_0t").c_str(),(base+"_1l_g1b_mbl_0t").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_1t1",std::make_shared<TH1D>((base+"_1l_g1b_mbl_1t1").c_str(),(base+"_1l_g1b_mbl_1t1").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_1t2",std::make_shared<TH1D>((base+"_1l_g1b_mbl_1t2").c_str(),(base+"_1l_g1b_mbl_1t2").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_1t3",std::make_shared<TH1D>((base+"_1l_g1b_mbl_1t3").c_str(),(base+"_1l_g1b_mbl_1t3").c_str(),nbins_1l,min_x,max_x_1l));
        
        // For Z->ll control region
        my_histos.emplace(base + "_2l_onZ",std::make_shared<TH1D>((base+"_2l_onZ").c_str(),(base+"_2l_onZ").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b",std::make_shared<TH1D>((base+"_2l_onZ_g1b").c_str(),(base+"_2l_onZ_g1b").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_g1t",std::make_shared<TH1D>((base+"_2l_onZ_g1b_g1t").c_str(),(base+"_2l_onZ_g1b_g1t").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher1",std::make_shared<TH1D>((base+"_2l_onZ_g1b_fisher1").c_str(),(base+"_2l_onZ_g1b_fisher1").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher2",std::make_shared<TH1D>((base+"_2l_onZ_g1b_fisher2").c_str(),(base+"_2l_onZ_g1b_fisher2").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher3",std::make_shared<TH1D>((base+"_2l_onZ_g1b_fisher3").c_str(),(base+"_2l_onZ_g1b_fisher3").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher4",std::make_shared<TH1D>((base+"_2l_onZ_g1b_fisher4").c_str(),(base+"_2l_onZ_g1b_fisher4").c_str(),nbins_1l,min_x,max_x_1l));
    }
}

void AnalyzeBackground::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while(tr.getNextEvent())
    {
        const auto& ntops                   = tr.getVar<int>("ntops");
        const auto& ntops_3jet              = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet              = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet              = tr.getVar<int>("ntops_1jet");
        const auto& runtype                 = tr.getVar<std::string>("runtype");
        const auto& filetag                 = tr.getVar<std::string>("filetag");
        const auto& JetID                   = tr.getVar<bool>("JetID");
        const auto& passTrigger             = tr.getVar<bool>("passTrigger");

        const auto& NJets_pt30          = tr.getVar<int>("NJets_pt30");
        const auto& NJets_pt45          = tr.getVar<int>("NJets_pt45");
        const auto& BJets_pt30          = tr.getVec<TLorentzVector>("BJets_pt30");
        const auto& NBJets_pt30         = tr.getVar<int>("NBJets_pt30");
        const auto& BJets_pt45          = tr.getVec<TLorentzVector>("BJets_pt45");
        const auto& NBJets_pt45         = tr.getVar<int>("NBJets_pt45");
        const auto& GoodLeptons         = tr.getVec<TLorentzVector>("GoodLeptons");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& HT_trigger          = tr.getVar<double>("HT_trigger");
        const auto& Mbl                 = tr.getVar<double>("Mbl");
        const auto& onZ                 = tr.getVar<bool>("onZ");

        const auto& fisher_bin1 = tr.getVar<bool>("fisher_bin1");
        const auto& fisher_bin2 = tr.getVar<bool>("fisher_bin2");
        const auto& fisher_bin3 = tr.getVar<bool>("fisher_bin3");
        const auto& fisher_bin4 = tr.getVar<bool>("fisher_bin4");
        const auto& bdt_bin1    = tr.getVar<bool>("bdt_bin1");
        const auto& bdt_bin2    = tr.getVar<bool>("bdt_bin2");
        const auto& bdt_bin3    = tr.getVar<bool>("bdt_bin3");
        const auto& bdt_bin4    = tr.getVar<bool>("bdt_bin4");

        const auto& passBaseline0l    = tr.getVar<bool>("passBaseline0l");
        const auto& passBaseline1l    = tr.getVar<bool>("passBaseline1l");
        const auto& passBaseline2lonZ = tr.getVar<bool>("passBaseline2lonZ");

        if (maxevents > 0 && tr.getEvtNum() >= maxevents) break;        

        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;
        
        // Exclude events with MadGraph HT > 100 from the DY inclusive sample
        if(filetag == "DYJetsToLL_M-50_Incl")
        {
            const double& madHT = tr.getVar<double>("madHT");
            if (madHT > 100) continue;
        }

        // Make sure event weight is not 0 for data
        double eventweight = 1.;
        if(runtype != "Data")
        {
            eventweight = tr.getVar<double>("Weight");
        }        

        // Exclude events with bad jets
        if(!JetID) continue;


        // -------------------------------
        // -- Basic event selection stuff
        // -------------------------------

        bool passNtop = ntops >= 1;
        bool passNb = NBJets_pt30 >= 1;
        
        // -------------------------
        // -- Check r(j) behavior --
        // -------------------------
        // temp
        bool passTrigger0l = true;
        bool passTrigger1l = true;
        bool passTrigger2l = true;

        int njets_rj = 0;
        std::vector<std::string> jettypes {"pt30", "pt45"};
        for(std::string jettype : jettypes)
        {
            if(jettype == "pt30")
                njets_rj = NJets_pt30;
            else if(jettype == "pt45")
                njets_rj = NJets_pt45;
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
                            if(ntops == 2)
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
            if(passBaseline1l)
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

                        if(ntops == 1 && ntops_1jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t1"]->Fill(njets_rj, eventweight);
                        if(ntops == 1 && ntops_2jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t2"]->Fill(njets_rj, eventweight);
                        if(ntops == 1 && ntops_3jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t3"]->Fill(njets_rj, eventweight);

                        if(fisher_bin1)
                            my_histos[base+"_1l_g1b_mbl_g1t_fisher1"]->Fill(njets_rj, eventweight);
                        if(fisher_bin2)
                            my_histos[base+"_1l_g1b_mbl_g1t_fisher2"]->Fill(njets_rj, eventweight);
                        if(fisher_bin3)
                            my_histos[base+"_1l_g1b_mbl_g1t_fisher3"]->Fill(njets_rj, eventweight);
                        if(fisher_bin4)
                            my_histos[base+"_1l_g1b_mbl_g1t_fisher4"]->Fill(njets_rj, eventweight);

                    }
                    else 
                    {
                        my_histos[base+"_1l_g1b_mbl_0t"]->Fill(njets_rj, eventweight);
                    }
                }
            }
          
            // Now for the Z->ll region
            if (passBaseline2lonZ)
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
                                    
                    if(passNtop)
                    {
                        my_histos[base+"_2l_onZ_g1b_g1t"]->Fill(njets_rj, eventweight);
                    }
                }
            }
        }
    }
}



void AnalyzeBackground::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }
    
}

