#ifndef Confg_h
#define Confg_h

#include "SusyAnaTools/Tools/NTupleReader.h"

#include "Framework/Framework/include/PrepNTupleVars.h"
#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Photon.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/DeepEventShape.h"
#include "Framework/Framework/include/BTagCorrector.h"
#include "Framework/Framework/include/ScaleFactors.h"
#include "Framework/Framework/include/PartialUnBlinding.h"
#include "Framework/Framework/include/StopGenMatch.h"
#include "Framework/Framework/include/MegaJetCombine.h"
#include "Framework/Framework/include/MakeMT2Hemispheres.h"

class Config
{
private:
    void registerModules(NTupleReader& tr, const std::vector<std::string>&& modules) const
    {
        const auto& runtype               = tr.getVar<std::string>("runtype");
        const auto& runYear               = tr.getVar<std::string>("runYear");
        const auto& DeepESMCfg            = tr.getVar<std::string>("DeepESMCfg");
        const auto& DeepESMCfg_NonIsoMuon = tr.getVar<std::string>("DeepESMCfg_NonIsoMuon");
        const auto& ModelFile             = tr.getVar<std::string>("ModelFile");
        const auto& leptonFileName        = tr.getVar<std::string>("leptonFileName");
        const auto& puFileName            = tr.getVar<std::string>("puFileName");
        const auto& meanFileName          = tr.getVar<std::string>("meanFileName");
        const auto& bjetFileName          = tr.getVar<std::string>("bjetFileName");
        const auto& bjetCSVFileName       = tr.getVar<std::string>("bjetCSVFileName");
        const auto& filetag               = tr.getVar<std::string>("filetag");
        
        for(const auto& module : modules)
        {
            if     (module=="PartialUnBlinding")           tr.emplaceModule<PartialUnBlinding>();
            else if(module=="PrepNTupleVars")              tr.emplaceModule<PrepNTupleVars>();
            else if(module=="RunTopTagger")                tr.emplaceModule<RunTopTagger>();
            else if(module=="Muon")                        tr.emplaceModule<Muon>();
            else if(module=="Electron")                    tr.emplaceModule<Electron>();
            else if(module=="Photon")                      tr.emplaceModule<Photon>();
            else if(module=="Jet")                         tr.emplaceModule<Jet>();
            else if(module=="BJet")                        tr.emplaceModule<BJet>();
            else if(module=="CommonVariables")             tr.emplaceModule<CommonVariables>();
            else if(module=="MakeMVAVariables")            tr.emplaceModule<MakeMVAVariables>();
            else if(module=="MakeMVAVariables_NonIsoMuon") tr.emplaceModule<MakeMVAVariables>(false, "", "NonIsoMuonJets_pt30");
            else if(module=="Baseline")                    tr.emplaceModule<Baseline>();
            else if(module=="StopGenMatch")                tr.emplaceModule<StopGenMatch>();
            else if(module=="MegaJetCombine")              tr.emplaceModule<MegaJetCombine>();
            else if(module=="MakeMT2Hemispheres_0l")       tr.emplaceModule<MakeMT2Hemispheres>("GoodJets_pt45", "NGoodJets_pt45", "_0l");
            else if(module=="MakeMT2Hemispheres_1l")       tr.emplaceModule<MakeMT2Hemispheres>("GoodJets_pt30", "NGoodJets_pt30", "_1l");
            else if(module=="DeepEventShape")              tr.emplaceModule<DeepEventShape>(DeepESMCfg, ModelFile);
            else if(module=="DeepEventShape_NonIsoMuon")   tr.emplaceModule<DeepEventShape>(DeepESMCfg_NonIsoMuon, ModelFile);
            if(runtype == "MC")
            {
                if     (module=="ScaleFactors")  tr.emplaceModule<ScaleFactors>(runYear, leptonFileName, puFileName, meanFileName);
                else if(module=="BTagCorrector")
                {
                    auto& bTagCorrector = tr.emplaceModule<BTagCorrector>(bjetFileName, "", bjetCSVFileName, filetag);
                    bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets", "GoodJets_pt30", "Jets_bJetTagDeepCSVtotb", "Jets_partonFlavor");
                }
            }
        }
    }

public:
    Config() 
    {
    }

    void setUp(NTupleReader& tr) const
    {
        //Get and make needed info
        const auto& filetag = tr.getVar<std::string>("filetag");
        const auto& analyzer = tr.getVar<std::string>("analyzer");
        const bool isSignal = (filetag.find("_stop") != std::string::npos || filetag.find("_mStop") != std::string::npos) ? true : false;

        std::string runYear, puFileName, DeepESMCfg, DeepESMCfg_NonIsoMuon, ModelFile, leptonFileName, bjetFileName, bjetCSVFileName, meanFileName;
        double Lumi, deepCSV_WP_loose, deepCSV_WP_medium, deepCSV_WP_tight;
        bool blind;

        if(filetag.find("2016") != std::string::npos)
        {
            runYear = "2016";
            Lumi = 35900.0;
            deepCSV_WP_loose  = 0.2217;
            deepCSV_WP_medium = 0.6321;
            deepCSV_WP_tight  = 0.8953;            
            puFileName = "PileupHistograms_0121_69p2mb_pm4p6.root";
            DeepESMCfg = "DeepEventShape_2016.cfg";
            DeepESMCfg_NonIsoMuon = "DeepEventShape_NonIsoMuon_2016.cfg";
            ModelFile = "keras_frozen_2016.pb";
            leptonFileName = "allInOne_leptonSF_2016.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_2016LegacySF_WP_V1.csv";
            meanFileName = "allInOne_SFMean.root";
            blind = false;
        }
        else if(filetag.find("2017") != std::string::npos)
        { 
            runYear = "2017";
            Lumi = 41525.0;
            deepCSV_WP_loose  = 0.1522;
            deepCSV_WP_medium = 0.4941;       
            deepCSV_WP_tight  = 0.8001;
            puFileName = "pu_ratio.root";
            DeepESMCfg = "DeepEventShape_2017.cfg";
            DeepESMCfg_NonIsoMuon = "DeepEventShape_NonIsoMuon_2017.cfg";
            ModelFile = "keras_frozen_2017.pb";
            leptonFileName = "allInOne_leptonSF_2017.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_94XSF_WP_V4_B_F.csv";
            meanFileName = "allInOne_SFMean.root";
            blind = false;
        }
        else if(filetag.find("2018pre") != std::string::npos) 
        {
            runYear = "2018pre";
            Lumi = 21071.0;
            deepCSV_WP_loose  = 0.1241;
            deepCSV_WP_medium = 0.4184;       
            deepCSV_WP_tight  = 0.7527;
            puFileName = "PileupHistograms_2018_69mb_pm5.root";
            DeepESMCfg = "DeepEventShape_2018pre.cfg";
            DeepESMCfg_NonIsoMuon = "DeepEventShape_NonIsoMuon_2018pre.cfg";
            ModelFile = "keras_frozen_2018pre.pb";
            leptonFileName = "allInOne_leptonSF_2018.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_102XSF_WP_V1.csv";
            meanFileName = "allInOne_SFMean.root";
            blind = true;
        }
        else if(filetag.find("2018post") != std::string::npos) 
        {
            runYear = "2018post";
            Lumi = 38654.0;
            deepCSV_WP_loose  = 0.1241;
            deepCSV_WP_medium = 0.4184;       
            deepCSV_WP_tight  = 0.7527;
            puFileName = "PileupHistograms_2018_69mb_pm5.root";
            DeepESMCfg = "DeepEventShape_2018post.cfg";
            DeepESMCfg_NonIsoMuon = "DeepEventShape_NonIsoMuon_2018post.cfg";
            ModelFile = "keras_frozen_2018post.pb";
            leptonFileName = "allInOne_leptonSF_2018.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_102XSF_WP_V1.csv";
            meanFileName = "allInOne_SFMean.root";
            blind = true;
        }

        tr.registerDerivedVar("runYear",runYear);
        tr.registerDerivedVar("Lumi",Lumi);
        tr.registerDerivedVar("deepCSV_WP_loose",deepCSV_WP_loose);
        tr.registerDerivedVar("deepCSV_WP_medium",deepCSV_WP_medium);
        tr.registerDerivedVar("deepCSV_WP_tight",deepCSV_WP_tight);
        tr.registerDerivedVar("isSignal",isSignal);
        tr.registerDerivedVar("DeepESMCfg",DeepESMCfg);
        tr.registerDerivedVar("DeepESMCfg_NonIsoMuon",DeepESMCfg_NonIsoMuon);
        tr.registerDerivedVar("ModelFile",ModelFile);        
        tr.registerDerivedVar("puFileName",puFileName);
        tr.registerDerivedVar("leptonFileName",leptonFileName);        
        tr.registerDerivedVar("bjetFileName",bjetFileName);        
        tr.registerDerivedVar("bjetCSVFileName",bjetCSVFileName);        
        tr.registerDerivedVar("meanFileName",meanFileName);        
        tr.registerDerivedVar("etaCut",2.4); 
        tr.registerDerivedVar("blind",blind);

        //Register Modules that are needed for each Analyzer
        if(analyzer=="Analyze0Lep" || analyzer=="Semra_Analyzer" || analyzer=="AnalyzeTopTagger")
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "RunTopTagger",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "MakeMVAVariables",
                "Baseline",
                "DeepEventShape",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="MakeNJetDists")
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "MakeMVAVariables",
                "Baseline",
                "DeepEventShape",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="CalculateBTagSF")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "Baseline",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="TwoLepAnalyzer" || analyzer=="Make2LInputTrees")
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "MakeMVAVariables",
                "Baseline",
                "DeepEventShape",
                "StopGenMatch",
                "MegaJetCombine",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="StealthHemispheres")
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "RunTopTagger",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "MakeMVAVariables",
                "Baseline",
                "DeepEventShape",
                "MakeMT2Hemispheres_0l"
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="AnalyzeLepTrigger" || analyzer=="CalculateSFMean")
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if (analyzer=="AnalyzeTest")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "PartialUnBlinding",
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "MakeMVAVariables",
                "Baseline",
                "DeepEventShape",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="Analyze1Lep")
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "RunTopTagger",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "MakeMVAVariables",
                "Baseline",
                "DeepEventShape",
                "MakeMVAVariables_NonIsoMuon",
                "DeepEventShape_NonIsoMuon",
                //"MakeMT2Hemispheres_1l",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "MakeMVAVariables",
                "Baseline",
                "DeepEventShape",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
    }
};

#endif
