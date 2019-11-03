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
        const auto& filetag             = tr.getVar<std::string>("filetag");

        const auto& passBaseline1l      = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline1l_NIM  = tr.getVar<bool>("passBaseline1l_NonIsoMuon");

        //------------------------------------
        //-- Print Event Number
        //------------------------------------

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        //-----------------------------------
        //  Initialize the tree
        //-----------------------------------
        
        std::set<std::string> variables = {
            //General variables
            "Lumi",
            "NVtx",
            "MET",
            "htDerivedweight",
            "bTagSF_EventWeightSimple_Central",
            "htDerivedweightFlat2000",
            "htDerivedweightNJet7",
            "prefiringScaleFactor",
            "totGoodMuonSF",
            "totGoodElectronSF",
            "puWeightCorr",
            "Weight",
            //Variables for Signal selection
            //"deepESM_val",
            //"NGoodJets_pt30",
            //"NGoodBJets_pt30",
            //"HT_trigger_pt30",
            //"NGoodElectrons",
            //"NGoodMuons",
            //"passBaseline1l_Good",
            //"totalEventWeight",
            //Variables for QCD CR selection
            "deepESM_valNonIsoMuon",
            "NNonIsoMuonJets_pt30",
            "HT_NonIsoMuon_pt30",
            "NNonIsoMuons",
            "passBaseline1l_NonIsoMuon",
            "totalEventWeightNIM",
        };
        if( runtype != "MC" ) {
            variables = {
            //General variables
            "Lumi",
            "NVtx",
            "MET",
            //Variables for Signal selection
            //"deepESM_val",
            //"NGoodJets_pt30",
            //"NGoodBJets_pt30",
            //"HT_trigger_pt30",
            //"NGoodElectrons",
            //"NGoodMuons",
            //"passBaseline1l_Good",
            //Variables for QCD CR selection
            "deepESM_valNonIsoMuon",
            "NNonIsoMuonJets_pt30",
            "HT_NonIsoMuon_pt30",
            "NNonIsoMuons",
            "passBaseline1l_NonIsoMuon"
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
        if( passBaseline1l_NIM ) {
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

