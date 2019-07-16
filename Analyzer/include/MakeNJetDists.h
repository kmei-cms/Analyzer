#ifndef MakeNJetDists_h
#define MakeNJetDists_h

#include "Analyzer/Analyzer/include/AnalyzeBase.h"
#include "Analyzer/Analyzer/include/Histo.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "SusyAnaTools/Tools/PileupWeights.h"

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

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <iostream>

class MakeNJetDists : public AnalyzeBase
{
private:
    std::vector<std::pair<std::string, std::string>> myVarSuffixPairs;

public:
    void InitHistos(const std::string& runtype)
    {
        TH1::SetDefaultSumw2();
        TH2::SetDefaultSumw2();

        my_Histos.emplace_back(new Histo1D("EventCounter", 2,  -1.1, 1.1, "eventCounter", {}, {}));

        std::vector<std::string> weightVec;
        if( runtype == "MC" )
        {
            myVarSuffixPairs = {{"",""}, {"JECup","_JECUp"}, {"JECdown","_JECDown"}, {"JERup","_JERUp"}, {"JERdown","_JERDown"}};
            weightVec = {"Lumi", "Weight"};

            //--------------------------------------------------------------------------------
            // Plots that are made only once and for MC only
            //--------------------------------------------------------------------------------

            //Normal nJet plots
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_btgUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_btgDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_lepUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_lepDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_isrUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_isrDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_fsrUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_fsrDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_isr2Up",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_isr2Down", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown_2"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_fsr2Up",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_fsr2Down", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_pdfUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_pdfDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightDown"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_htUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp", "puWeightCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_htDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "puWeightCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_puUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysUpCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_puDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysDownCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_sclUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_sclDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightDown"}));

            //Shifted nJet plots
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_btgUp",   6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_btgDown", 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_lepUp",   6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_lepDown", 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isrUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isrDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsrUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsrDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isr2Up",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isr2Down", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown_2"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsr2Up",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsr2Down", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_pdfUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_pdfDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightDown"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_htUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp", "puWeightCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_htDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "puWeightCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_puUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysUpCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_puDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysDownCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_sclUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_sclDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightDown"}));
       
            for(int i = 0; i < 4; i++)
            {
                std::string index = std::to_string(i+1);

                //Normal nJet plots
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_btgUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_btgDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_lepUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_lepDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_isrUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_isrDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown"}));

                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_fsrUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_fsrDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown"}));

                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_isr2Up",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_isr2Down", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown_2"}));

                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_fsr2Up",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_fsr2Down", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_pdfUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightUp"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_pdfDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightDown"}));

                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_htUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp", "puWeightCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_htDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "puWeightCorr", "prefiringScaleFactor"}));
                
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_puUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysUpCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_puDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysDownCorr", "prefiringScaleFactor"}));

                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_sclUp",   6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightUp"}));
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_sclDown", 6,  7.0, 13.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightDown"}));

                //Shifted nJet plots
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_btgUp",   6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_btgDown", 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_lepUp",   6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_lepDown", 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isrUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isrDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown"}));
            
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsrUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsrDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown"}));
            
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isr2Up",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isr2Down", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown_2"}));
        
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsr2Up",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsr2Down", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_pdfUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightUp"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_pdfDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightDown"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_htUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp","puWeightCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_htDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "puWeightCorr", "prefiringScaleFactor"}));
                
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_puUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight","puSysUpCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_puDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysDownCorr", "prefiringScaleFactor"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_sclUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightUp"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_sclDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightDown"}));
            }            
        }
        else
        {
            myVarSuffixPairs = {{"",""}};
            weightVec = {};
        }    

        //--------------------------------------------------------------------------------
        // Plots that are made with JEC/R variation and Data
        //--------------------------------------------------------------------------------
        for(const auto& pair : myVarSuffixPairs)
        {
            const std::string& s = pair.first;
            const std::string& n = pair.second;

            std::vector<std::string> weightVecAll;
            if( runtype == "MC" )
            {
                weightVecAll = {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central"+s, "totGoodElectronSF"+s, "totGoodMuonSF"+s, "htDerivedweight"+s, "puWeightCorr"+s, "prefiringScaleFactor"+s};
            }
            else
            {
                weightVecAll = {};
            }

            //-----------------------------------------------------------------
            // NJet plots
            //-----------------------------------------------------------------       
            my_Histos.emplace_back(new Histo1D("h_njets"+n,         6, 7.0, 13.0, "NGoodJets_pt30_inclusive"+s,       {}, weightVec));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted"+n,  6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift"+s, {}, weightVec));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l"+n,        6, 7.0, 13.0, "NGoodJets_pt30_inclusive"+s,       {"passBaseline1l_Good"+s}, weightVecAll));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l"+n, 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift"+s, {"passBaseline1l_Good"+s}, weightVecAll));
            for(int i = 0; i < 4; i++)
            {
                std::string index = std::to_string(i+1);
                my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+n,        6, 7.0, 13.0, "NGoodJets_pt30_inclusive"+s,       {"passBaseline1l_Good"+s,"deepESM_bin"+index+s}, weightVecAll));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+n, 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift"+s, {"passBaseline1l_Good"+s,"deepESM_bin"+index+s}, weightVecAll));
            }
        }
    
    }//END of init histos

    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false)
    {
        const auto& runtype = tr.getVar<std::string>("runtype");
        const auto& filetag = tr.getVar<std::string>("filetag");
        const auto& runYear = tr.getVar<std::string>("runYear");
        const auto& DeepESMCfg = tr.getVar<std::string>("DeepESMCfg");
        const auto& ModelFile = tr.getVar<std::string>("ModelFile");

        //-------------------------------------
        //-- Initialize histograms to be filled
        //-------------------------------------
        InitHistos(runtype);

        for(const auto& pair : myVarSuffixPairs)
        {
            const std::string& myVarSuffix = pair.first;
            if(myVarSuffix == "") continue;
            //RunTopTagger rtt("TopTagger.cfg", myVarSuffix);
            Muon muon(myVarSuffix);
            Electron electron(myVarSuffix);
            Photon photon(myVarSuffix);
            Jet jet(myVarSuffix);
            BJet bjet(myVarSuffix);
            CommonVariables commonVariables(myVarSuffix);
            MakeMVAVariables makeMVAVariables(false, myVarSuffix);
            Baseline baseline(myVarSuffix);
            DeepEventShape deepEventShape(DeepESMCfg, ModelFile, "Info", true, myVarSuffix);
            BTagCorrectorTemplate<double> bTagCorrector("allInOne_BTagEff.root","", false, filetag);
            bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets"+myVarSuffix, "Jets"+myVarSuffix+"_bJetTagDeepCSVtotb", "Jets"+myVarSuffix+"_partonFlavor", myVarSuffix);
            //Pileup_SysTemplate<double> pileup("PileupHistograms_0121_69p2mb_pm4p6.root");
            std::string scaleFactorHistoFileName = (runYear == "2017") ? "allInOne_leptonSF_2017.root" : "allInOne_leptonSF_2016.root";
            std::string puRootFileName = ( runYear == "2017" ) ? "pu_ratio.root" : "PileupHistograms_0121_69p2mb_pm4p6.root";
            ScaleFactors scaleFactors( scaleFactorHistoFileName, puRootFileName, "allInOne_SFMean.root", myVarSuffix );
        
            //tr.registerFunction(rtt);
            tr.registerFunction(muon);
            tr.registerFunction(electron);
            tr.registerFunction(photon);
            tr.registerFunction(jet);
            tr.registerFunction(bjet);
            tr.registerFunction(commonVariables);
            tr.registerFunction(makeMVAVariables);
            tr.registerFunction(baseline);
            tr.registerFunction(deepEventShape);
            tr.registerFunction(bTagCorrector);
            //tr.registerFunction(pileup);
            tr.registerFunction(scaleFactors);
        }

        while( tr.getNextEvent() )
        {
            //------------------------------------
            //-- Print Event Number
            //------------------------------------
            const bool breakLoop = printEventNum(maxevents, tr.getEvtNum());
            if(breakLoop) break;

            //-----------------------------------
            //-- Fill Histograms Below
            //-----------------------------------
            Fill(tr);
        }//END of while tr.getNextEvent loop   
    }//END of function
};

#endif
