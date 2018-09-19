#define MakeNJetDists_cxx
#include "Analyzer/Analyzer/include/MakeNJetDists.h"
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

#include "SusyAnaTools/Tools/BTagCalibrationStandalone.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"

MakeNJetDists::MakeNJetDists()
{
    InitHistos();
}


void MakeNJetDists::InitHistos()
{ 
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains

    std::vector<std::string> nlepTags   { "0l", "1l" };
    std::vector<std::string> sfTags     { "std", "btg", "qcd", "pdf", "pup", "lep", "all"}; 
    std::vector<std::string> uncertTags { "central", "up", "down" }; 

    //Histogram with no scale factors applied
    
    for( std::string nlepTag : nlepTags ) {
        for( std::string sfTag: sfTags ) {
            for( std::string uncertTag : uncertTags ) {
                
                if( nlepTag == "0l" && sfTag == "lep" ) continue; //If it is the hadronic channel, there is no lepton scale factor

                if( ( sfTag == "all" || sfTag == "std" ) && uncertTag != "central" ) continue; //No fluctuations with all scale factors implemented and when no scale factors are applied

                my_histos.emplace( "h_njets_"+sfTag+"_"+nlepTag+"_"+uncertTag, std::make_shared<TH1D>( ( "h_njets_"+sfTag+"_"+nlepTag+"_"+uncertTag ).c_str() , ( "h_njets_"+sfTag+"_"+nlepTag+"_"+uncertTag ).c_str() , 16, 0, 16 ) );
            }
        }
    }

//    my_histos.emplace( "h_temp", std::make_shared<TH1D>( "h_temp", "htemp", 30, 0, 2 ) );
//    my_histos.emplace( "h_temp1", std::make_shared<TH1D>( "h_temp1", "htemp1", 30, -1, 1 ) );
//    my_histos.emplace( "h_temp2", std::make_shared<TH1D>( "h_temp2", "htemp2", 30, 0, 2 ) );
//    my_histos.emplace( "h_temp3", std::make_shared<TH1D>( "h_temp3", "htemp3", 30, 0, 2 ) );

}//END of init histos

void MakeNJetDists::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& filetag             = tr.getVar<std::string>("filetag");

        const auto& passMadHT           = tr.getVar<bool>("passMadHT");

        const auto& passBaseline0l      = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBaseline1l      = tr.getVar<bool>("passBaseline1l_Good");

        const auto& NJets_pt30          = tr.getVar<int>("NGoodJets_pt30");
        const auto& NJets_pt45          = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodElectrons      = tr.getVar<int>("NGoodElectrons");
        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");

        //------------------------------------
        //-- Print Event Number
        //------------------------------------

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1024 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

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
        
        const auto& PileupWeight        = tr.getVar<float>("_PUweightFactor");
        const auto& PileupWeightUp      = tr.getVar<float>("_PUSysUp");
        const auto& PileupWeightDown    = tr.getVar<float>("_PUSysDown");

        const auto& eleLepWeight        = tr.getVar<double>("totGoodElectronSF");
        const auto& eleLepWeightUp      = tr.getVar<double>("totGoodElectronSF_Up");
        const auto& eleLepWeightDown    = tr.getVar<double>("totGoodElectronSF_Down");
        
        const auto& muLepWeight         = tr.getVar<double>("totGoodMuonSF");
        const auto& muLepWeightUp       = tr.getVar<double>("totGoodMuonSF_Up");
        const auto& muLepWeightDown     = tr.getVar<double>("totGoodMuonSF_Down");

        const auto& bTagWeight          = tr.getVar<float>("bTagSF_EventWeightSimple_Central");
        const auto& bTagWeightUp        = tr.getVar<float>("bTagSF_EventWeightSimple_Up");
        const auto& bTagWeightDown      = tr.getVar<float>("bTagSF_EventWeightSimple_Down");

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
        //--
        //-----------------------------------


        //------------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        
        if( passBaseline0l ) {
            my_histos["h_njets_std_0l_central"]->Fill( NJets_pt45, eventweight );
            my_histos["h_njets_btg_0l_central"]->Fill( NJets_pt45, eventweight*bTagWeight );
            my_histos["h_njets_qcd_0l_central"]->Fill( NJets_pt45, eventweight*scaleWeight );
            my_histos["h_njets_pdf_0l_central"]->Fill( NJets_pt45, eventweight*PDFWeight );
            my_histos["h_njets_pup_0l_central"]->Fill( NJets_pt45, eventweight*PileupWeight );
            my_histos["h_njets_all_0l_central"]->Fill( NJets_pt45, eventweight*bTagWeight*scaleWeight*PDFWeight*PileupWeight );
            
            my_histos["h_njets_btg_0l_up"]->Fill( NJets_pt45, eventweight*bTagWeightUp );
            my_histos["h_njets_qcd_0l_up"]->Fill( NJets_pt45, eventweight*scaleWeightUp );
            my_histos["h_njets_pdf_0l_up"]->Fill( NJets_pt45, eventweight*PDFWeightUp );
            my_histos["h_njets_pup_0l_up"]->Fill( NJets_pt45, eventweight*PileupWeightUp );

            my_histos["h_njets_btg_0l_down"]->Fill( NJets_pt45, eventweight*bTagWeightDown );
            my_histos["h_njets_qcd_0l_down"]->Fill( NJets_pt45, eventweight*scaleWeightDown );
            my_histos["h_njets_pdf_0l_down"]->Fill( NJets_pt45, eventweight*PDFWeightDown );
            my_histos["h_njets_pup_0l_down"]->Fill( NJets_pt45, eventweight*PileupWeightDown );
        }

        double lepWeight        = 0.0;
        double lepWeightUp      = 0.0;
        double lepWeightDown    = 0.0;

        if( NGoodElectrons == 1 ) {
            lepWeight       = eleLepWeight;
            lepWeightUp     = eleLepWeightUp;
            lepWeightDown   = eleLepWeightDown;
        }

        if( NGoodMuons == 1 ) {
            lepWeight       = muLepWeight;
            lepWeightUp     = muLepWeightUp;
            lepWeightDown   = muLepWeightDown;
        }

        if( passBaseline1l ) {
            my_histos["h_njets_std_1l_central"]->Fill( NJets_pt30, eventweight );
            my_histos["h_njets_btg_1l_central"]->Fill( NJets_pt30, eventweight*bTagWeight );
            my_histos["h_njets_qcd_1l_central"]->Fill( NJets_pt30, eventweight*scaleWeight );
            my_histos["h_njets_pdf_1l_central"]->Fill( NJets_pt30, eventweight*PDFWeight );
            my_histos["h_njets_pup_1l_central"]->Fill( NJets_pt30, eventweight*PileupWeight );
            my_histos["h_njets_lep_1l_central"]->Fill( NJets_pt30, eventweight*lepWeight );
            my_histos["h_njets_all_1l_central"]->Fill( NJets_pt30, eventweight*bTagWeight*scaleWeight*PDFWeight*PileupWeight );
            
            my_histos["h_njets_btg_1l_up"]->Fill( NJets_pt30, eventweight*bTagWeightUp );
            my_histos["h_njets_qcd_1l_up"]->Fill( NJets_pt30, eventweight*scaleWeightUp );
            my_histos["h_njets_pdf_1l_up"]->Fill( NJets_pt30, eventweight*PDFWeightUp );
            my_histos["h_njets_pup_1l_up"]->Fill( NJets_pt30, eventweight*PileupWeightUp );
            my_histos["h_njets_lep_1l_up"]->Fill( NJets_pt30, eventweight*lepWeightUp );
            
            my_histos["h_njets_btg_1l_down"]->Fill( NJets_pt30, eventweight*bTagWeightDown );
            my_histos["h_njets_qcd_1l_down"]->Fill( NJets_pt30, eventweight*scaleWeightDown );
            my_histos["h_njets_pdf_1l_down"]->Fill( NJets_pt30, eventweight*PDFWeightDown );
            my_histos["h_njets_pup_1l_down"]->Fill( NJets_pt30, eventweight*PileupWeightDown );
            my_histos["h_njets_lep_1l_down"]->Fill( NJets_pt30, eventweight*lepWeightDown );
        }
    }//END of while tr.getNextEvent loop   
}//END of function
      
void MakeNJetDists::WriteHistos( TFile* outfile ) 
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }   
}

