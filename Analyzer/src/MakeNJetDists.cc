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
    InitHistos();
}

void MakeNJetDists::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    myVarSuffixs = {"", "JECup", "JECdown", "JERup", "JERdown"};
    for(const auto& s : myVarSuffixs)
    {
        my_Histos.emplace_back(new Histo1D<int>("h_njets"+s, 8, 6.5, 14.5, "NGoodJets_pt30_inclusive"+s));
        std::vector<std::string> weightVec = {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central"+s, "totGoodElectronSF"+s, "totGoodMuonSF"+s, "htDerivedweight"+s};

        my_Histos.emplace_back(new Histo1D<int>("h_njets_pt30_1l"+s, 8, 6.5, 14.5, "NGoodJets_pt30_inclusive"+s, {"passBaseline1l_Good"+s}, weightVec));
        for(int i = 0; i < 4; i++)
        {
            std::string index = std::to_string(i+1);
            my_Histos.emplace_back(new Histo1D<int>("h_njets_pt30_1l_deepESMbin"+index+s, 8, 6.5, 14.5, "NGoodJets_pt30_inclusive"+s, {"passBaseline1l_Good"+s,"deepESM_bin"+index+s}, weightVec));
        }        
    }
    
}//END of init histos

void MakeNJetDists::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    const auto& runtype = tr.getVar<std::string>("runtype");
    const auto& filetag = tr.getVar<std::string>("filetag");

//    for(const auto& myVarSuffix : myVarSuffixs)
//    {
//        if(myVarSuffix == "") continue;
//        //RunTopTagger rtt("TopTagger.cfg", myVarSuffix);
//        Muon muon(myVarSuffix);
//        Electron electron(myVarSuffix);
//        Photon photon(myVarSuffix);
//        Jet jet(myVarSuffix);
//        BJet bjet(myVarSuffix);
//        CommonVariables commonVariables(myVarSuffix);
//        MakeMVAVariables makeMVAVariables(false, myVarSuffix);
//        Baseline baseline(myVarSuffix);
//        DeepEventShape deepEventShape("DeepEventShape.cfg", "Info", true, myVarSuffix);
//        BTagCorrectorTemplate<double> bTagCorrector("allInOne_BTagEff.root","", false, filetag);
//        bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets"+myVarSuffix, "Jets"+myVarSuffix+"_bDiscriminatorCSV", "Jets"+myVarSuffix+"_partonFlavor", myVarSuffix);
//        //Pileup_SysTemplate<double> pileup("PileupHistograms_0121_69p2mb_pm4p6.root");
//        ScaleFactors scaleFactors("allInOne_leptonSF_Moriond17.root", myVarSuffix);
//        
//        //tr.registerFunction(rtt);
//        tr.registerFunction(muon);
//        tr.registerFunction(electron);
//        tr.registerFunction(photon);
//        tr.registerFunction(jet);
//        tr.registerFunction(bjet);
//        tr.registerFunction(commonVariables);
//        tr.registerFunction(makeMVAVariables);
//        tr.registerFunction(baseline);
//        tr.registerFunction(deepEventShape);
//        tr.registerFunction(bTagCorrector);
//        //tr.registerFunction(pileup);
//        tr.registerFunction(scaleFactors);
//    }

    while( tr.getNextEvent() )
    {
        //------------------------------------
        //-- Print Event Number
        //------------------------------------
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        
        //-----------------------------------
        //-- Make sure you are running over MC
        //-- Doesn't really make sense to run 
        //--   on Data (also not all variables
        //--   are there
        //-----------------------------------
        if( runtype != "MC" ) 
        {
            std::cerr<<"Please run over an MC file since these scale factors should not be applied to data!!"<<std::endl;
            break;
        }
               
        //-------------------------------------
        //-- Make sure we do not double DY events 
        //------------------------------------
        if( !passMadHT ) continue; 
        
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
