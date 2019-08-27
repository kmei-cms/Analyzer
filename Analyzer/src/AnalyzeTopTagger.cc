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

//mandatory includes to use top tagger
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
    TRandom3* rand = new TRandom3(123);
   
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
        const auto& GoodLeptons            = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons"); 
        const auto& NGoodBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45");
        const auto& GoodBJets_pt45         = tr.getVec<bool>("GoodBJets_pt45"); 
        const auto& pass_baseline_0l       = tr.getVar<bool>("pass_baseline_0l");  
     
        // -------------------
        // -- Define weight
        // -------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double htDerivedScaleFactor = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
 
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight = tr.getVar<double>("Weight");
            const auto& lumi   = tr.getVar<double>("Lumi");
            eventweight        = lumi*Weight;

            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor        = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");

            weight *= eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        // --------------------------------------
        // -- Calculate DeltaR between 2 bjets 
        // --------------------------------------
        double dR_bjet1_bjet2 = -1;
        if (NGoodBJets_pt45 == 2) {
            std::vector<TLorentzVector> bjets;
            for(int ijet = 0; ijet < Jets.size(); ijet++) {
                if(!GoodBJets_pt45[ijet]) continue;
                bjets.push_back(Jets.at(ijet));
            }
            dR_bjet1_bjet2 = bjets[0].DeltaR(bjets[1]);
        }
        double pass_ge1dRbjets = (dR_bjet1_bjet2 >= 1.0); 


        std::vector<std::pair<std::string, bool>> ttbarCuts =
        {   
            {"pass_baseline_0l", pass_baseline_0l && pass_ge1dRbjets},
        };
        hists.fillWithCutFlow(ttbarCuts, tr, weight, rand); 

    }
    
    delete rand;

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
