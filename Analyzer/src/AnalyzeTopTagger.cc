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

AnalyzeTopTagger::AnalyzeTopTagger() : hists("histos"), histNjet7("Njet7"), histNjet8("Njet8"), histNjet9("Njet9"), histNjet10("Njet10"), histNjet11("Njet11"), histNjet12("Njet12"), histNjet13("Njet13"), histNjet14("Njet14"), histNjet15("Njet15") 
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

        const auto& runtype                = tr.getVar<std::string>("runtype");
        const auto& filetag                = tr.getVar<std::string>("filetag");
        const auto& NGoodLeptons           = tr.getVar<int>("NGoodLeptons");
        const auto& dR_bjets               = tr.getVar<double>("dR_bjets");
        const auto& Jets                   = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt45          = tr.getVec<bool>("GoodJets_pt45");
        const auto& NGoodJets_pt45         = tr.getVar<int>("NGoodJets_pt45");
        const auto& GoodBJets_pt45         = tr.getVec<bool>("GoodBJets_pt45");
        const auto& NGoodBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45");
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
        // ----------------
        // -- baseline cuts
        // ----------------
        std::vector<std::pair<std::string, bool>> ttbarCuts =
        {   
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
        };
        hists.fillWithCutFlow(ttbarCuts, tr, weight, &rand);

        // -----------------------------
        // -- baseline cuts + Njets == 7
        // -----------------------------
        std::vector<std::pair<std::string, bool>> Njets7 =
        {
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
            {"Njet7"         , NGoodJets_pt45 == 7},
        };
        histNjet7.fillWithCutFlow(Njets7, tr, weight, &rand);

        // -----------------------------
        // -- baseline cuts + Njets == 8
        // -----------------------------
        std::vector<std::pair<std::string, bool>> Njets8 =
        {
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets}, 
            {"Njet8"         , NGoodJets_pt45 == 8},
        };
        histNjet8.fillWithCutFlow(Njets8, tr, weight, &rand);

        // -----------------------------
        // -- baseline cuts + Njets == 9
        // -----------------------------
        std::vector<std::pair<std::string, bool>> Njets9 =
        {
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
            {"Njet9"         , NGoodJets_pt45 == 9},
        };
        histNjet9.fillWithCutFlow(Njets9, tr, weight, &rand);        
        
        // ------------------------------
        // -- baseline cuts + Njets == 10
        // ------------------------------
        std::vector<std::pair<std::string, bool>> Njets10 =
        {
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
            {"Njet10"        , NGoodJets_pt45 == 10},
        };
        histNjet10.fillWithCutFlow(Njets10, tr, weight, &rand);

        // ------------------------------
        // -- baseline cuts + Njets == 11
        // ------------------------------    
        std::vector<std::pair<std::string, bool>> Njets11 =
        {
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
            {"Njet11"        , NGoodJets_pt45 == 11},
        };
        histNjet11.fillWithCutFlow(Njets11, tr, weight, &rand);

        // ------------------------------
        // -- baseline cuts + Njets == 12
        // ------------------------------
        std::vector<std::pair<std::string, bool>> Njets12 =
        {
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
            {"Njet12"        , NGoodJets_pt45 == 12},
        };
        histNjet12.fillWithCutFlow(Njets12, tr, weight, &rand);

        // ------------------------------
        // -- baseline cuts + Njets == 13
        // ------------------------------
        std::vector<std::pair<std::string, bool>> Njets13 =
        {
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
            {"Njet13"        , NGoodJets_pt45 == 13},
        };
        histNjet13.fillWithCutFlow(Njets13, tr, weight, &rand);   
 
        // ------------------------------    
        // -- baseline cuts + Njets == 14
        // ------------------------------
        std::vector<std::pair<std::string, bool>> Njets14 =
        {
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
            {"Njet10"        , NGoodJets_pt45 == 14},
        };
        histNjet14.fillWithCutFlow(Njets14, tr, weight, &rand);
        
        // ------------------------------
        // -- baseline cuts + Njets == 15
        // ------------------------------
        std::vector<std::pair<std::string, bool>> Njets14 =
        {
            {"passBaseline0l", passBaseline0l && pass_ge1dRbjets},
            {"Njet15"        , NGoodJets_pt45 == 15},
        };
        histNjet14.fillWithCutFlow(Njets14, tr, weight, &rand);

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
