#ifndef Confg_h
#define Confg_h

#include "SusyAnaTools/Tools/BTagCalibrationStandalone.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "SusyAnaTools/Tools/PileupWeights.h"
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
        //Get needed info
        const auto& runtype = tr.getVar<std::string>("runtype");
        const auto& runYear = tr.getVar<std::string>("runYear");
        const auto& filetag = tr.getVar<std::string>("filetag");
        const auto& DeepESMCfg = tr.getVar<std::string>("DeepESMCfg");
        const auto& ModelFile = tr.getVar<std::string>("ModelFile");
        const auto& analyzer = tr.getVar<std::string>("analyzer");

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
        Pileup_SysTemplate<double>* pileup = nullptr;
        ScaleFactors* scaleFactors = nullptr;
        if( runtype == "MC" )
        {
            bTagCorrector = new BTagCorrectorTemplate<double>("allInOne_BTagEff.root","", false, filetag);
            bTagCorrector->SetVarNames("GenParticles_PdgId", "Jets", "Jets_bJetTagDeepCSVtotb", "Jets_partonFlavor");
            pileup = new Pileup_SysTemplate<double>("PileupHistograms_0121_69p2mb_pm4p6.root");
            const std::string scaleFactorHistoFileName = "allInOne_leptonSF_"+runYear+".root";
            std::string puFileName;
            if     (runYear == "2016") puFileName = "PileupHistograms_0121_69p2mb_pm4p6.root";
            else if(runYear == "2017") puFileName = "pu_ratio.root";
            else if(runYear == "2018") puFileName = "pu_ratio.root";
            scaleFactors = new ScaleFactors( scaleFactorHistoFileName, puFileName );
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
