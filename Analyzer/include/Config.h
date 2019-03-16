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
#include "Framework/Framework/include/RunFisher.h"
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

    void registerModules(NTupleReader& tr)
    {
        const auto& runtype = tr.getVar<std::string>("runtype");
        const auto& runYear = tr.getVar<std::string>("runYear");
        const auto& filetag = tr.getVar<std::string>("filetag");
        const auto& DeepESMCfg = tr.getVar<std::string>("DeepESMCfg");
        const auto& ModelFile = tr.getVar<std::string>("ModelFile");

        // Define classes/functions that add variables on the fly
        PartialUnBlinding partUnBlind;
        PrepNTupleVars prep;
        RunTopTagger rtt;
        Muon muon;
        Electron electron;
        Photon photon;
        Jet jet;
        BJet bjet;
        RunFisher runFisher("v3");
        CommonVariables commonVariables;
        MakeMVAVariables makeMVAVariables;
        Baseline baseline;
        DeepEventShape deepEventShape(DeepESMCfg,ModelFile);
        
        // Register classes/functions that add variables on the fly
        tr.registerFunction(partUnBlind);
        tr.registerFunction(prep);
        tr.registerFunction(rtt);
        tr.registerFunction(muon);
        tr.registerFunction(electron);
        tr.registerFunction(photon);
        tr.registerFunction(jet);
        tr.registerFunction(bjet);
        tr.registerFunction(runFisher);
        tr.registerFunction(commonVariables);
        tr.registerFunction(makeMVAVariables);
        tr.registerFunction(baseline);
        tr.registerFunction(deepEventShape);
        
        if( runtype == "MC" ) 
        {
            //std::string eosPath =   "root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/";
            BTagCorrectorTemplate<double> bTagCorrector("allInOne_BTagEff.root","", false, filetag);
            bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets", "Jets_bDiscriminatorCSV", "Jets_partonFlavor");
            Pileup_SysTemplate<double> pileup("PileupHistograms_0121_69p2mb_pm4p6.root");
            std::string scaleFactorHistoFileName = (runYear == "2017") ? "allInOne_leptonSF_2017.root" : "allInOne_leptonSF_2016.root";
            ScaleFactors scaleFactors( scaleFactorHistoFileName );

            tr.registerFunction(bTagCorrector);
            tr.registerFunction(pileup);
            tr.registerFunction(scaleFactors);
        }        
    }
};

#endif
