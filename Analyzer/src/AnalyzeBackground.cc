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
        my_histos.emplace(base + "_1l_g1b_mbl_0t",std::make_shared<TH1D>((base+"_1l_g1b_mbl_0t").c_str(),(base+"_1l_g1b_mbl_0t").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_1t1",std::make_shared<TH1D>((base+"_1l_g1b_mbl_1t1").c_str(),(base+"_1l_g1b_mbl_1t1").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_1t2",std::make_shared<TH1D>((base+"_1l_g1b_mbl_1t2").c_str(),(base+"_1l_g1b_mbl_1t2").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_1l_g1b_mbl_1t3",std::make_shared<TH1D>((base+"_1l_g1b_mbl_1t3").c_str(),(base+"_1l_g1b_mbl_1t3").c_str(),nbins_1l,min_x,max_x_1l));
        
        // For Z->ll control region
        my_histos.emplace(base + "_2l_onZ",std::make_shared<TH1D>((base+"_2l_onZ").c_str(),(base+"_2l_onZ").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b",std::make_shared<TH1D>((base+"_2l_onZ_g1b").c_str(),(base+"_2l_onZ_g1b").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl",std::make_shared<TH1D>((base+"_2l_onZ_g1b_nombl").c_str(),(base+"_2l_onZ_g1b_nombl").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_g1t",std::make_shared<TH1D>((base+"_2l_onZ_g1b_g1t").c_str(),(base+"_2l_onZ_g1b_g1t").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt1",std::make_shared<TH1D>((base+"_2l_onZ_g1b_nombl_bdt1").c_str(),(base+"_2l_onZ_g1b_nombl_bdt1").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt2",std::make_shared<TH1D>((base+"_2l_onZ_g1b_nombl_bdt2").c_str(),(base+"_2l_onZ_g1b_nombl_bdt2").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt3",std::make_shared<TH1D>((base+"_2l_onZ_g1b_nombl_bdt3").c_str(),(base+"_2l_onZ_g1b_nombl_bdt3").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_nombl_bdt4",std::make_shared<TH1D>((base+"_2l_onZ_g1b_nombl_bdt4").c_str(),(base+"_2l_onZ_g1b_nombl_bdt4").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher1",std::make_shared<TH1D>((base+"_2l_onZ_g1b_fisher1").c_str(),(base+"_2l_onZ_g1b_fisher1").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher2",std::make_shared<TH1D>((base+"_2l_onZ_g1b_fisher2").c_str(),(base+"_2l_onZ_g1b_fisher2").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher3",std::make_shared<TH1D>((base+"_2l_onZ_g1b_fisher3").c_str(),(base+"_2l_onZ_g1b_fisher3").c_str(),nbins_1l,min_x,max_x_1l));
        my_histos.emplace(base + "_2l_onZ_g1b_fisher4",std::make_shared<TH1D>((base+"_2l_onZ_g1b_fisher4").c_str(),(base+"_2l_onZ_g1b_fisher4").c_str(),nbins_1l,min_x,max_x_1l));
    }
}

void AnalyzeBackground::Loop(NTupleReader& tr, double weight, int maxevents, std::string filetag, bool isQuiet)
{
    while(tr.getNextEvent())
    {
        const auto& ntops                   = tr.getVar<int>("ntops");
        const auto& ntops_3jet              = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet              = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet              = tr.getVar<int>("ntops_1jet");
        const auto& runtype                 = tr.getVar<std::string>("runtype");
        const auto& JetID                   = tr.getVar<bool>("JetID");
        const auto& TriggerNames            = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass             = tr.getVec<int>("TriggerPass");

        const auto& NJets_pt30          = tr.getVar<int>("NJets_pt30");
        const auto& NJets_pt45          = tr.getVar<int>("NJets_pt45");
        const auto& BJets_pt30          = tr.getVec<TLorentzVector>("BJets_pt30");
        const auto& NBJets_pt30         = tr.getVar<int>("NBJets_pt30");
        const auto& BJets_pt45          = tr.getVec<TLorentzVector>("BJets_pt45");
        const auto& NBJets_pt45         = tr.getVar<int>("NBJets_pt45");
        const auto& GoodMuons           = tr.getVec<TLorentzVector>("GoodMuons");
        const auto& GoodMuonsCharge     = tr.getVec<int>("GoodMuonsCharge");
        const auto& GoodElectrons       = tr.getVec<TLorentzVector>("GoodElectrons");
        const auto& GoodElectronsCharge = tr.getVec<int>("GoodElectronsCharge");
        const auto& GoodLeptons         = tr.getVec<TLorentzVector>("GoodLeptons");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& HT_trigger          = tr.getVar<double>("HT_trigger");
        const auto& Mbl                 = tr.getVar<double>("Mbl");

        const auto& fisher_bin1 = tr.getVar<bool>("fisher_bin1");
        const auto& fisher_bin2 = tr.getVar<bool>("fisher_bin2");
        const auto& fisher_bin3 = tr.getVar<bool>("fisher_bin3");
        const auto& fisher_bin4 = tr.getVar<bool>("fisher_bin4");
        const auto& bdt_bin1    = tr.getVar<bool>("bdt_bin1");
        const auto& bdt_bin2    = tr.getVar<bool>("bdt_bin2");
        const auto& bdt_bin3    = tr.getVar<bool>("bdt_bin3");
        const auto& bdt_bin4    = tr.getVar<bool>("bdt_bin4");

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
        // -- Basic event selection stuff
        // -------------------------------

        // Check whether event would pass the trigger requirement
        bool passTrigger = true;
        if ( !( HT_trigger>500 && NJets_pt45>=6 ) ) 
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

        bool passBaseline0l = NGoodLeptons==0 && NJets_pt45>=6 && HT_trigger > 500 && NBJets_pt45 >= 1;
        bool passBaseline1l = NGoodLeptons==1 && NJets_pt30>=6 ;
        bool passBaseline2l = NGoodLeptons==2;

        bool passNtop = ntops >= 1;
        bool passNb = NBJets_pt30 >= 1;
        bool onZ = false;
        bool passMbl_2l = false;
        if ( (GoodMuons.size() == 2) && (GoodMuonsCharge[0] != GoodMuonsCharge[1]) )
        {
            double mll = (GoodMuons[0] + GoodMuons[1]).M();
            if( mll > 81.2 && mll < 101.2)
                onZ = true;          

            // check whether a bl pair passes the M(b,l) cut
            for (TLorentzVector myb : BJets_pt30)
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
            for (TLorentzVector myb : BJets_pt30)
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

                        if(ntops == 1 && ntops_1jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t1"]->Fill(njets_rj, eventweight);
                        if(ntops == 1 && ntops_2jet == 1)
                            my_histos[base+"_1l_g1b_mbl_1t2"]->Fill(njets_rj, eventweight);
                        if(ntops == 1 && ntops_3jet == 1)
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

