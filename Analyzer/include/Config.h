#ifndef Confg_h
#define Confg_h

#include "SusyAnaTools/Tools/BTagCalibrationStandalone.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
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
#include "Framework/Framework/include/ScaleFactors.h"
#include "Framework/Framework/include/PartialUnBlinding.h"

class Config
{
public:
    Config()
    {
    }

    void registerModules(NTupleReader& tr) const
    {
        //Get and make needed info
        const auto& runtype = tr.getVar<std::string>("runtype");
        const auto& filetag = tr.getVar<std::string>("filetag");
        const auto& analyzer = tr.getVar<std::string>("analyzer");

        std::string runYear, puFileName, DeepESMCfg, ModelFile, leptonFileName, bjetFileName, bjetCSVFileName, meanFileName;
        double Lumi, deepCSV_WP_loose, deepCSV_WP_medium, deepCSV_WP_tight;
        if(filetag.find("2016") != std::string::npos)
        {
            runYear = "2016";
            Lumi = 35900.0;
            deepCSV_WP_loose  = 0.2217;
            deepCSV_WP_medium = 0.6321;
            deepCSV_WP_tight  = 0.8953;            
            puFileName = "PileupHistograms_0121_69p2mb_pm4p6.root";
            DeepESMCfg = "DeepEventShape_2016.cfg";
            ModelFile = "keras_frozen_2016.pb";
            leptonFileName = "allInOne_leptonSF_2016.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_2016LegacySF_V1.csv";
            meanFileName = "allInOne_SFMean.root";
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
            ModelFile = "keras_frozen_2017.pb";
            leptonFileName = "allInOne_leptonSF_2017.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_94XSF_V4_B_F.csv";
            meanFileName = "allInOne_SFMean.root";
        }
        else if(filetag.find("2018") != std::string::npos) 
        {
            runYear = "2018";
            Lumi = 59740.0;
            deepCSV_WP_loose  = 0.1241;
            deepCSV_WP_medium = 0.4184;       
            deepCSV_WP_tight  = 0.7527;
            puFileName = "pu_ratio.root";
            DeepESMCfg = "DeepEventShape_2018.cfg";
            ModelFile = "keras_frozen_2018.pb";
            leptonFileName = "allInOne_leptonSF_2017.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_102XSF_V1.csv";
            meanFileName = "allInOne_SFMean.root";
        }

        const bool isSignal = (filetag.find("_stop") != std::string::npos || filetag.find("_mStop") != std::string::npos) ? true : false;
        const bool doQCDCR = false;//bool to determine to use qcd control region

        tr.registerDerivedVar("runYear",runYear);
        tr.registerDerivedVar("etaCut",2.4); 
        tr.registerDerivedVar("Lumi",Lumi);
        tr.registerDerivedVar("deepCSV_WP_loose",deepCSV_WP_loose);
        tr.registerDerivedVar("deepCSV_WP_medium",deepCSV_WP_medium);
        tr.registerDerivedVar("deepCSV_WP_tight",deepCSV_WP_tight);
        tr.registerDerivedVar("isSignal",isSignal);
        tr.registerDerivedVar("DeepESMCfg",DeepESMCfg);
        tr.registerDerivedVar("ModelFile",ModelFile);        
        tr.registerDerivedVar("puFileName",puFileName);
        tr.registerDerivedVar("leptonFileName",leptonFileName);        
        tr.registerDerivedVar("bjetFileName",bjetFileName);        
        tr.registerDerivedVar("bjetCSVFileName",bjetCSVFileName);        
        tr.registerDerivedVar("meanFileName",meanFileName);        
        tr.registerDerivedVar("blind",true);
        tr.registerDerivedVar("doQCDCR",doQCDCR);

        //Create all possible modules
        PartialUnBlinding partUnBlind;
        PrepNTupleVars prep;
        RunTopTagger rtt;
        Muon muon;
        Electron electron;
        Photon photon;
        Jet jet;
        BJet bjet;
        CommonVariables commonVariables;
        MakeMVAVariables makeMVAVariables;
        Baseline baseline;
        DeepEventShape deepEventShape(DeepESMCfg,ModelFile);
        BTagCorrectorTemplate<double>* bTagCorrector = nullptr;
        ScaleFactors* scaleFactors = nullptr;
        if( runtype == "MC" && analyzer != "CalculateBTagSF")
        {
            bTagCorrector = new BTagCorrectorTemplate<double>(bjetFileName, "", bjetCSVFileName, false, filetag);
            bTagCorrector->SetVarNames("GenParticles_PdgId", "Jets", "Jets_bJetTagDeepCSVtotb", "Jets_partonFlavor");
            scaleFactors = new ScaleFactors( runYear, leptonFileName, puFileName, meanFileName );
        }

        //Register Modules that are needed for each Analyzer
        if(analyzer=="Analyze1Lep" || analyzer=="Analyze0Lep" || analyzer=="Semra_Analyzer")
        {
            tr.registerFunction(partUnBlind);
            tr.registerFunction(prep);                   
            tr.registerFunction(rtt);
            tr.registerFunction(muon);
            tr.registerFunction(electron);
            tr.registerFunction(photon);
            tr.registerFunction(jet);
            tr.registerFunction(bjet);
            tr.registerFunction(commonVariables);
            tr.registerFunction(makeMVAVariables);
            tr.registerFunction(baseline);
            tr.registerFunction(deepEventShape);        
            if( runtype == "MC")
            {
                tr.registerFunction(*bTagCorrector);
                tr.registerFunction(*scaleFactors);
            }
        }
        else if(analyzer=="MakeNJetDists")
        {
            tr.registerFunction(partUnBlind);
            tr.registerFunction(prep);                   
            tr.registerFunction(muon);
            tr.registerFunction(electron);
            tr.registerFunction(photon);
            tr.registerFunction(jet);
            tr.registerFunction(bjet);
            tr.registerFunction(commonVariables);
            tr.registerFunction(makeMVAVariables);
            tr.registerFunction(baseline);
            tr.registerFunction(deepEventShape);        
            if( runtype == "MC")
            {
                tr.registerFunction(*bTagCorrector);
                tr.registerFunction(*scaleFactors);
            }            
        }
        else if(analyzer=="CalculateBTagSF")
        {
            tr.registerFunction(partUnBlind);
            tr.registerFunction(prep);                   
            tr.registerFunction(muon);
            tr.registerFunction(electron);
            tr.registerFunction(photon);
            tr.registerFunction(jet);
            tr.registerFunction(bjet);
            tr.registerFunction(commonVariables);
            tr.registerFunction(makeMVAVariables);
            tr.registerFunction(baseline);
        }
        else
        {
            tr.registerFunction(partUnBlind);
            tr.registerFunction(prep);                   
            tr.registerFunction(muon);
            tr.registerFunction(electron);
            tr.registerFunction(photon);
            tr.registerFunction(jet);
            tr.registerFunction(bjet);
            tr.registerFunction(commonVariables);
            tr.registerFunction(makeMVAVariables);
            tr.registerFunction(baseline);
            tr.registerFunction(deepEventShape);        
            if( runtype == "MC")
            {
                tr.registerFunction(*bTagCorrector);
                tr.registerFunction(*scaleFactors);
            }            
        }
    }
};

#endif
