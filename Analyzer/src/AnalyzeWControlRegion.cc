#define AnalyzeWControlRegion_cxx
#include "Analyzer/Analyzer/include/AnalyzeWControlRegion.h"
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

AnalyzeWControlRegion::AnalyzeWControlRegion() : initHistos(false)
{
}

void AnalyzeWControlRegion::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains

    for(auto& mycut : cutMap)
    {
        my_histos.emplace("h_njets_"+mycut.first,std::make_shared<TH1D>(("h_njets_"+mycut.first).c_str(),("h_njets_"+mycut.first).c_str(),15,0,15));
        my_histos.emplace("h_mT_"+mycut.first,std::make_shared<TH1D>(("h_mT_"+mycut.first).c_str(),("h_mT_"+mycut.first).c_str(),40,0,200));
        my_histos.emplace("h_HT_"+mycut.first,std::make_shared<TH1D>(("h_HT_"+mycut.first).c_str(),("h_HT_"+mycut.first).c_str(),30,0,1500));
    }   
}

void AnalyzeWControlRegion::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while(tr.getNextEvent())
    {
        //const auto& ntops                   = tr.getVar<int>("ntops");
        //const auto& ntops_3jet              = tr.getVar<int>("ntops_3jet");
        //const auto& ntops_2jet              = tr.getVar<int>("ntops_2jet");
        //const auto& ntops_1jet              = tr.getVar<int>("ntops_1jet");
        const auto& runtype                 = tr.getVar<std::string>("runtype");
        const auto& filetag                 = tr.getVar<std::string>("filetag");
        const auto& JetID                   = tr.getVar<bool>("JetID");
        const auto& passTrigger             = tr.getVar<bool>("passTrigger");

        //const auto& NJets_pt30          = tr.getVar<int>("NJets_pt30");
        //const auto& NJets_pt45          = tr.getVar<int>("NJets_pt45");
        //const auto& BJets_pt30          = tr.getVec<TLorentzVector>("BJets_pt30");
        //const auto& NBJets_pt30         = tr.getVar<int>("NBJets_pt30");
        //const auto& BJets_pt45          = tr.getVec<TLorentzVector>("BJets_pt45");
        //const auto& NBJets_pt45         = tr.getVar<int>("NBJets_pt45");
        const auto& NGoodBJets_pt30     = tr.getVar<int>("NGoodBJets_pt30");
        const auto& GoodLeptons         = tr.getVec<TLorentzVector>("GoodLeptons");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& GoodMuons           = tr.getVec<TLorentzVector>("GoodMuons");
        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");
        const auto& GoodElectrons       = tr.getVec<TLorentzVector>("GoodElectrons");
        const auto& NGoodElectrons      = tr.getVar<int>("NGoodElectrons");
        //const auto& HT_trigger          = tr.getVar<double>("HT_trigger");
        const auto& Mbl                 = tr.getVar<double>("Mbl");
        const auto& MET                 = tr.getVar<double>("MET");
        const auto& METPhi                 = tr.getVar<double>("METPhi");

        //const auto& passBaseline0l    = tr.getVar<bool>("passBaseline0l");
        //const auto& passBaseline1l    = tr.getVar<bool>("passBaseline1l");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");

        if (maxevents > 0 && tr.getEvtNum() >= maxevents) break;        

        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;
        
        // Make sure event weight is not 0 for data
        double eventweight = 1.;
        if(runtype != "Data")
        {
            eventweight = tr.getVar<double>("Weight");
        }        

        // -------------------------------
        // -- Basic event selection stuff
        // -------------------------------

        // basic event cleaning
        bool goodEvent = JetID && passMadHT && passTrigger;
        // Exclude events with bad jets
        if(!goodEvent) continue;

        // lepton cuts
        bool pass_1mu = (runtype != "Data" || filetag == "Data_SingleMuon")     && NGoodMuons == 1;
        bool pass_1el = (runtype != "Data" || filetag == "Data_SingleElectron") && NGoodElectrons == 1;
        bool pass_1l = pass_1mu || pass_1el;
        // compute mT 
        TLorentzVector metLV;
        metLV.SetPtEtaPhiE(MET,0,METPhi,MET);
        double mT = utility::calcMT(GoodLeptons[0], metLV);
        bool pass_mT = mT > 30 && mT < 100;

        // jet cuts
        bool pass_0b = NGoodBJets_pt30 == 0;


        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map 
        {
            {"1l"                                 , goodEvent && pass_1l                                                             },
            {"1l_mT"                              , goodEvent && pass_1l && pass_mT                                                  },
            {"1l_mT_0b"                           , goodEvent && pass_1l && pass_mT && pass_0b                                       },
        };
        
        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map);
            initHistos = true;
        }


    }
}



void AnalyzeWControlRegion::WriteHistos(TFile* outfile)
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

