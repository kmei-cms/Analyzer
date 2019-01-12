#define MakeNJetDists_cxx
#include "Analyzer/Analyzer/include/MakeNJetDists.h"
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
#include <TRandom3.h>
#include <iostream>

MakeNJetDists::MakeNJetDists()
{
}

void MakeNJetDists::InitHistos(const std::string& runtype)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    std::vector<std::string> weightVec;
    if( runtype == "MC" )
    {
        myVarSuffixPairs = {{"",""}, {"JECup","_JECUp"}, {"JECdown","_JECDown"}, {"JERup","_JERUp"}, {"JERdown","_JERDown"}};
        weightVec = {"Lumi", "Weight"};

        //--------------------------------------------------------------------------------
        // Plots that are made only once and for MC only
        //--------------------------------------------------------------------------------

        //Normal nJet plots
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_btgUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor"}));
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_btgDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor"}));

        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_lepUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "prefiringScaleFactor"}));
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_lepDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "prefiringScaleFactor"}));

        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_isrUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRUp"}));
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_isrDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRDown"}));

        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_fsrUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRUp"}));
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_fsrDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRDown"}));

        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_isr2Up",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_isr2Down", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRDown_2"}));

        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_fsr2Up",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_fsr2Down", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_pdfUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PDFweightUp"}));
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_pdfDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PDFweightDown"}));

        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_htUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp", "prefiringScaleFactor"}));
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_htDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "prefiringScaleFactor"}));

        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_sclUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "scaleWeightUp"}));
        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_sclDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "scaleWeightDown"}));

        //Shifted nJet plots
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_btgUp",   8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor"}));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_btgDown", 8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor"}));

        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_lepUp",   8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "prefiringScaleFactor"}));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_lepDown", 8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "prefiringScaleFactor"}));

        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isrUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRUp"}));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isrDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRDown"}));

        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsrUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRUp"}));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsrDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRDown"}));

        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isr2Up",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isr2Down", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRDown_2"}));

        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsr2Up",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsr2Down", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_pdfUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PDFweightUp"}));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_pdfDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PDFweightDown"}));

        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_htUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp", "prefiringScaleFactor"}));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_htDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "prefiringScaleFactor"}));

        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_sclUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "scaleWeightUp"}));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_sclDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "scaleWeightDown"}));
       
        for(int i = 0; i < 4; i++)
        {
            std::string index = std::to_string(i+1);

            //Normal nJet plots
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_btgUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_btgDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_lepUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_lepDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_isrUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_isrDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRDown"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_fsrUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_fsrDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRDown"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_isr2Up",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_isr2Down", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRDown_2"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_fsr2Up",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_fsr2Down", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_pdfUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PDFweightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_pdfDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PDFweightDown"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_htUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_htDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_sclUp",   8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "scaleWeightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+"_sclDown", 8,  7.0, 15.0, "NGoodJets_pt30_inclusive", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "scaleWeightDown"}));

            //Shifted nJet plots
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_btgUp",   8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_btgDown", 8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_lepUp",   8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_lepDown", 8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isrUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isrDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRDown"}));
            
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsrUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsrDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRDown"}));
            
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isr2Up",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isr2Down", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_ISRDown_2"}));
        
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsr2Up",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsr2Down", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_pdfUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PDFweightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_pdfDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "PDFweightDown"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_htUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_htDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_sclUp",   8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "scaleWeightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_sclDown", 8,  0.0, 8.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "prefiringScaleFactor", "scaleWeightDown"}));
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
            weightVecAll = {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central"+s, "totGoodElectronSF"+s, "totGoodMuonSF"+s, "htDerivedweight"+s, "prefiringScaleFactor"+s};
        }
        else
        {
            weightVecAll = {};
        }

        //-----------------------------------------------------------------
        // NJet plots
        //-----------------------------------------------------------------       
        my_Histos.emplace_back(new Histo1D("h_njets"+n,         8, 7.0, 15.0, "NGoodJets_pt30_inclusive"+s,       {}, weightVec));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted"+n,  8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift"+s, {}, weightVec));

        my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l"+n,        8, 7.0, 15.0, "NGoodJets_pt30_inclusive"+s,       {"passBaseline1l_Good"+s}, weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l"+n, 8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift"+s, {"passBaseline1l_Good"+s}, weightVecAll));
        for(int i = 0; i < 4; i++)
        {
            std::string index = std::to_string(i+1);
            my_Histos.emplace_back(new Histo1D("h_njets_pt30_1l_D"+index+n,        8, 7.0, 15.0, "NGoodJets_pt30_inclusive"+s,       {"passBaseline1l_Good"+s,"deepESM_bin"+index+s}, weightVecAll));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+n, 8, 0.0,  8.0, "NGoodJets_pt30_inclusive_shift"+s, {"passBaseline1l_Good"+s,"deepESM_bin"+index+s}, weightVecAll));
        }
    }
    
}//END of init histos

void MakeNJetDists::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    const auto& runtype = tr.getVar<std::string>("runtype");
    const auto& filetag = tr.getVar<std::string>("filetag");
    const auto& runYear = tr.getVar<std::string>("runYear");
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
        DeepEventShape deepEventShape("DeepEventShape.cfg", "Info", true, myVarSuffix);
        BTagCorrectorTemplate<double> bTagCorrector("allInOne_BTagEff.root","", false, filetag);
        bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets"+myVarSuffix, "Jets"+myVarSuffix+"_bDiscriminatorCSV", "Jets"+myVarSuffix+"_partonFlavor", myVarSuffix);
        //Pileup_SysTemplate<double> pileup("PileupHistograms_0121_69p2mb_pm4p6.root");
        std::string scaleFactorHistoFileName = (runYear == "2017") ? "allInOne_leptonSF_2017.root" : "allInOne_leptonSF_Moriond17.root";
        ScaleFactors scaleFactors( scaleFactorHistoFileName, "allInOne_HtSFDist_2016.root", myVarSuffix );
        
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
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        for(auto& h : my_Histos)
        {
            h->Fill(tr);
        }
    }//END of while tr.getNextEvent loop   
}//END of function
      
void MakeNJetDists::WriteHistos( TFile* outfile ) 
{
    outfile->cd();
    
    for(const auto& h : my_Histos)
    {
        h->Write();
    }
}
