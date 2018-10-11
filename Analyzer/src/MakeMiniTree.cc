#define MakeMiniTree_cxx
#include "Analyzer/Analyzer/include/MakeMiniTree.h"
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
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "Framework/Framework/include/SetUpTopTagger.h"

#include "SusyAnaTools/Tools/MiniTupleMaker.h"
#include "SusyAnaTools/Tools/PileupWeights.h"
#include "SusyAnaTools/Tools/BTagCalibrationStandalone.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"

MakeMiniTree::MakeMiniTree()
{
    InitHistos();
}


void MakeMiniTree::InitHistos()
{

}//END of init histos

void MakeMiniTree::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& filetag             = tr.getVar<std::string>("filetag");

        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passTriggerMC       = tr.getVar<bool>("passTriggerMC");

        const auto& passBaseline0l      = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBaseline1l      = tr.getVar<bool>("passBaseline1l_Good");

        const auto& NJets_pt30          = tr.getVar<int>("NGoodJets_pt30");
        const auto& NJets_pt45          = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodElectrons      = tr.getVar<int>("NGoodElectrons");
        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");
        
        const auto& HT                  = tr.getVar<double>("HT");
        const auto& MadHT               = tr.getVar<double>("madHT");

        //------------------------------------
        //-- Print Event Number
        //------------------------------------

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        //-----------------------------------
        //-- Make sure you are running over MC
        //-- Doesn't really make sense to run 
        //--   on Data (also not all variables
        //--   are there
        //-----------------------------------
        
        if( runtype != "MC" ) {
            std::cerr<<"Please run over an MC file since these scale factors should not be applied to data!!"<<std::endl;
            break;
        }

        const auto& scaleWeight         = tr.getVar<double>("scaleWeightNom");
        const auto& scaleWeightUp       = tr.getVar<double>("scaleWeightUp");
        const auto& scaleWeightDown     = tr.getVar<double>("scaleWeightDown");
        
        const auto& PDFWeight           = tr.getVar<double>("PDFweightNom");
        const auto& PDFWeightUp         = tr.getVar<double>("PDFweightUp");
        const auto& PDFWeightDown       = tr.getVar<double>("PDFweightDown");
        
        const auto& PileupWeight        = tr.getVar<double>("_PUweightFactor");
        const auto& PileupWeightUp      = tr.getVar<double>("_PUSysUp");
        const auto& PileupWeightDown    = tr.getVar<double>("_PUSysDown");

        const auto& eleLepWeight        = tr.getVar<double>("totGoodElectronSF");
        const auto& eleLepWeightUp      = tr.getVar<double>("totGoodElectronSF_Up");
        const auto& eleLepWeightDown    = tr.getVar<double>("totGoodElectronSF_Down");
        
        const auto& muLepWeight         = tr.getVar<double>("totGoodMuonSF");
        const auto& muLepWeightUp       = tr.getVar<double>("totGoodMuonSF_Up");
        const auto& muLepWeightDown     = tr.getVar<double>("totGoodMuonSF_Down");

        //const auto& bTagWeight          = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
        //const auto& bTagWeightUp        = tr.getVar<double>("bTagSF_EventWeightSimple_Up");
        //const auto& bTagWeightDown      = tr.getVar<double>("bTagSF_EventWeightSimple_Down");
        
        const auto& deepESMValue        = tr.getVar<double>("deepESM_val");
        const auto& deepESMbin1         = tr.getVar<bool>("deepESM_bin1");
        const auto& deepESMbin2         = tr.getVar<bool>("deepESM_bin2");
        const auto& deepESMbin3         = tr.getVar<bool>("deepESM_bin3");
        const auto& deepESMbin4         = tr.getVar<bool>("deepESM_bin4");
        
        //-----------------------------------
        //  Initialize the tree
        //-----------------------------------
        
        std::set<std::string> variables = {
            "Weight",
            "deepESM_val",
            "NGoodJets_pt30",
            "NGoodBJets_pt30",
            "NGoodBJets_pt30_tight",
            "Mbl",
            "ntops_1jet",
            "ntops_2jet",
            "ntops_3jet",
        };

        if( tr.isFirstEvent() ) {
            std::string myTreeName = "myMiniTree";
            myTree = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
            myMiniTuple = new MiniTupleMaker( myTree );
            myMiniTuple->setTupleVars(variables);
            myMiniTuple->initBranches(tr);
        }

        //-------------------------------------
        //-- Make sure we do not double DY events 
        //------------------------------------
        
        if( !passMadHT ) continue; 
        
        //------------------------------------
        //-- Get the Proper Event Weight
        //------------------------------------

        double eventweight          = 1.0;

        double Lumi = 35900;
        const auto& Weight = tr.getVar<double>("Weight");
            
        eventweight         = Lumi*Weight;

        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        
        if( passBaseline1l && passTriggerMC ) {

            myMiniTuple->fill();
        }

    }//END of while tr.getNextEvent loop   
}//END of function
      
void MakeMiniTree::WriteHistos( TFile* outfile ) 
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }

    myTree->Write();
    delete myTree;    
    delete myMiniTuple;

}

