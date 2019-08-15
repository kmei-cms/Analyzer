#define MakeMiniTree_cxx
#include "Analyzer/Analyzer/include/MakeMiniTree.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

#include "SusyAnaTools/Tools/MiniTupleMaker.h"

MakeMiniTree::MakeMiniTree()
{
    InitHistos();
}


void MakeMiniTree::InitHistos()
{
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ); 
}//END of init histos

void MakeMiniTree::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );
       
        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& passBaseline1l      = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline1lNIM   = tr.getVar<bool>("passBaseline1l_NonIsoMuon");

        const auto& doQCDCR             = tr.getVar<bool>("doQCDCR");

        //------------------------------------
        //-- Print Event Number
        //------------------------------------

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        //-----------------------------------
        //  Initialize the tree
        //-----------------------------------
        
        std::set<std::string> variables = {
            "Lumi",
            "deepESM_val",
            "filetag",
            "NGoodJets_pt30",
            "NGoodBJets_pt30",
            "Mbl",
            "HT_trigger_pt30",
            "NGoodElectrons",
            "NGoodMuons",
            "MET",
            "NVtx",
            "bTagSF_EventWeightSimple_Central",
            "htDerivedweight",
            "htDerivedweightMG",
            "prefiringScaleFactor",
            "totGoodMuonSF",
            "totGoodElectronSF",
            "puWeightCorr",
            "Weight",
            "totalEventWeight",
            "totalEventWeightMG",
            "totalEventWeightNIM",
        };
        if( runtype != "MC" ) {
            variables = {
            "Lumi",
            "deepESM_val",
            "filetag",
            "NGoodJets_pt30",
            "NGoodBJets_pt30",
            "Mbl",
            "NVtx",
            "HT_trigger_pt30",
            "NGoodElectrons",
            "NGoodMuons",
            "MET",
            };
        }
        if( tr.isFirstEvent() ) {
            std::string myTreeName = "myMiniTree";
            myTree = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
            myMiniTuple = new MiniTupleMaker( myTree );
            myMiniTuple->setTupleVars(variables);
            myMiniTuple->initBranches(tr);
        }
        
        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        if( (doQCDCR && passBaseline1lNIM) || (!doQCDCR && passBaseline1l) ) {
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

