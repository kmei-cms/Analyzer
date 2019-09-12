#define AnalyzeTopTagger_cxx
#include "Analyzer/Analyzer/include/AnalyzeTopTagger.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

// Mandatory includes to use top tagger
#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "Framework/Framework/include/SetUpTopTagger.h"

AnalyzeTopTagger::AnalyzeTopTagger() : hists("histos") 
{
    InitHistos();
}

void AnalyzeTopTagger::InitHistos()
{
    TH1::SetDefaultSumw2();
}

void AnalyzeTopTagger::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    TRandom3 rand(123);
   
    while(tr.getNextEvent())
    {
        const auto& eventCounter           = tr.getVar<int>("eventCounter");

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );        

        const auto& MET                    = tr.getVar<double>("MET");
        const auto& HT                     = tr.getVar<double>("HT");
        const auto& NJets                  = tr.getVar<int>("NJets");
        const auto& ntops_3jet             = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet             = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet             = tr.getVar<int>("ntops_1jet");
        const auto& Jets                   = tr.getVec<TLorentzVector>("Jets");
        const auto& Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const auto* ttr                    = tr.getVar<TopTaggerResults*>("ttr");
        const auto& hadtopdaughters        = tr.getVec<std::vector<const TLorentzVector*>>("hadtopdaughters");
        const auto& hadtops                = tr.getVec<TLorentzVector>("hadtops");
        const auto& neutralinos            = tr.getVec<TLorentzVector>("neutralinos");
        const auto& singlinos              = tr.getVec<TLorentzVector>("singlinos");
        const auto& singlets               = tr.getVec<TLorentzVector>("singlets");
        const auto& hadtops_idx            = tr.getVec<int>("hadtops_idx");

        const auto& runtype                = tr.getVar<std::string>("runtype");
        const auto& filetag                = tr.getVar<std::string>("filetag");
        const auto& NGoodLeptons           = tr.getVar<int>("NGoodLeptons");
        const auto& dR_bjets               = tr.getVar<double>("dR_bjets");
        const auto& GoodJets_pt45          = tr.getVec<bool>("GoodJets_pt45");
        const auto& ntops                  = tr.getVar<int>("ntops"); 
        const bool  passBaseline0l         = tr.getVar<bool>("passBaseline0l_Good");     
        const bool  pass_ge1dRbjets        = (dR_bjets >= 1.0); 

        // -------------------
        // -- Define weight
        // -------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double bTagScaleFactor      = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight = tr.getVar<double>("Weight");
            const auto& lumi   = tr.getVar<double>("Lumi");
            eventweight        = lumi*Weight;

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");

            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        // -----------------------------------------------------
        // -- Create variables for setStealthStopVar function
        // -----------------------------------------------------
        auto& goodjets_pt45 = tr.createDerivedVec<TLorentzVector>("GoodJets_pt45_tlv");
        for (int i = 0; i < Jets.size(); ++i )
        {
            if (!GoodJets_pt45[i]) continue;
            goodjets_pt45.emplace_back(Jets.at(i));            
        }

        // -----------------------------------------
        // -- Fill the fake rate study histograms 
        // -----------------------------------------
        std::vector<std::pair<std::string, bool>> ttbarCuts =
        {   
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
        };
        hists.fillWithCutFlow(ttbarCuts, tr, weight, &rand); 
    }
}

void AnalyzeTopTagger::WriteHistos(TFile* outfile)
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
   
    hists.save(outfile);
    outfile->Write();
    outfile->Close();
}
