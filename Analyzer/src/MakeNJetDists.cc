#define MakeNJetDists_cxx
#include "Analyzer/Analyzer/include/MakeNJetDists.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/DeepEventShape.h"

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

    // Declare all your histograms here, that way we can fill them for multiple chains
    std::vector<std::string> ptTags      { "pt30", "pt45" };
    std::vector<std::string> nlepTags    { "0l", "1l" };
    std::vector<std::string> deepESMTags { "0", "1", "2" , "3", "4" };
    std::vector<std::string> sfTags      { "std", "qcd", "pdf", "pup", "lep", "all"}; 
    std::vector<std::string> uncertTags  { "central", "up", "down" };

    //Histogram with no scale factors applied
        for( auto& ptTag : ptTags ) 
    {
        for( auto& nlepTag : nlepTags ) 
        {
            for( auto& deepESMTag : deepESMTags ) 
            {
                for( auto& sfTag: sfTags ) 
                {
                    for( auto& uncertTag : uncertTags ) 
                    {                        
                        if( nlepTag == "0l" && sfTag == "lep"  ) continue; //If it is the hadronic channel, there is no lepton scale factor
                        if( nlepTag == "0l" && ptTag == "pt30" ) continue; //For the zero lepton channel, the pt threshold is 45 GeV;
                        if( nlepTag == "1l" && ptTag == "pt45" ) continue; //For the one lepton channel, the pt threshold is 30 GeV;
                        if( ( sfTag == "all" || sfTag == "std" ) && uncertTag != "central" ) continue; //No fluctuations with all scale factors implemented and when no scale factors are applied
                        
                        std::string myHistoName = "h_njets_"+ptTag+"_"+nlepTag;
                        std::string myHTHistoName = "h_ht_"+ptTag+"_"+nlepTag;
                        std::string myMadHTHistoName = "h_madht_"+ptTag+"_"+nlepTag;
                        
                        if( deepESMTag != "0" ) myHistoName = myHistoName+"_deepESMbin"+deepESMTag;
                        if( sfTag != "std" ) myHistoName = myHistoName+"_"+sfTag;
                        if( uncertTag != "central" ) myHistoName = myHistoName+"_"+uncertTag;

                        if( my_histos.find(myHistoName) == my_histos.end() )      my_histos.emplace( myHistoName,      std::make_shared<TH1D>( ( myHistoName ).c_str() ,      ( myHistoName ).c_str() ,       8, 6.5, 14.5 ) );
                        if( my_histos.find(myHTHistoName) == my_histos.end() )    my_histos.emplace( myHTHistoName,    std::make_shared<TH1D>( ( myHTHistoName ).c_str() ,    ( myHTHistoName ).c_str() ,    28, 200, 2500 ) );
                        if( my_histos.find(myMadHTHistoName) == my_histos.end() ) my_histos.emplace( myMadHTHistoName, std::make_shared<TH1D>( ( myMadHTHistoName ).c_str() , ( myMadHTHistoName ).c_str() , 28, 200, 2500 ) );
                    }
                }
            }
        }
    }
}//END of init histos

void MakeNJetDists::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    std::string myVarSuffix = "JECup";

    //Jet jet(myVarSuffix);
    //BJet bjet(myVarSuffix);
    //CommonVariables commonVariables(myVarSuffix);
    //MakeMVAVariables makeMVAVariables(false, myVarSuffix, true);
    ////Baseline baseline(myVarSuffix);
    //DeepEventShape deepEventShape(myVarSuffix);
    //
    //tr.registerFunction(jet);
    //tr.registerFunction(bjet);
    //tr.registerFunction(commonVariables);
    //tr.registerFunction(makeMVAVariables);
    ////tr.registerFunction(baseline);
    //tr.registerFunction(deepEventShape);

    while( tr.getNextEvent() )
    {
        //------------------------------------
        //-- Print Event Number
        //------------------------------------
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passBaseline0l      = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBaseline1l      = tr.getVar<bool>("passBaseline1l_Good");
        const auto& NJets_pt30          = tr.getVar<int>("NGoodJets_pt30");
        const auto& NJets_pt45          = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodElectrons      = tr.getVar<int>("NGoodElectrons");
        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");
        const auto& HT                  = tr.getVar<double>("HT");
        const auto& MadHT               = tr.getVar<double>("madHT");
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
        std::vector<bool> deepESMbins   = {deepESMbin1, deepESMbin2, deepESMbin3, deepESMbin4};

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
        
        //------------------------------------
        //-- Get the Proper Event Weight
        //------------------------------------
        double eventweight = 1.0;
        double Lumi = 35900;
        const auto& Weight = tr.getVar<double>("Weight");
        eventweight        = Lumi*Weight;

        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        if( passBaseline0l ) 
        {
            my_histos["h_njets_pt45_0l"]->Fill( NJets_pt45, eventweight );
            my_histos["h_njets_pt45_0l_pdf"]->Fill( NJets_pt45, eventweight*PDFWeight );
            my_histos["h_njets_pt45_0l_pup"]->Fill( NJets_pt45, eventweight*PileupWeight );
            
            my_histos["h_njets_pt45_0l_qcd_up"]->Fill( NJets_pt45, eventweight*scaleWeightUp );
            my_histos["h_njets_pt45_0l_pdf_up"]->Fill( NJets_pt45, eventweight*PDFWeightUp );
            my_histos["h_njets_pt45_0l_pup_up"]->Fill( NJets_pt45, eventweight*PileupWeightUp );

            my_histos["h_njets_pt45_0l_qcd_down"]->Fill( NJets_pt45, eventweight*scaleWeightDown );
            my_histos["h_njets_pt45_0l_pdf_down"]->Fill( NJets_pt45, eventweight*PDFWeightDown );
            my_histos["h_njets_pt45_0l_pup_down"]->Fill( NJets_pt45, eventweight*PileupWeightDown );

            for(int index = 0; index < deepESMbins.size(); index++)
            {
                if( deepESMbins.at(index) ) 
                {
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)       ]->Fill( NJets_pt45, eventweight );
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)+"_qcd"]->Fill( NJets_pt45, eventweight*scaleWeight );
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)+"_pdf"]->Fill( NJets_pt45, eventweight*PDFWeight );
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)+"_pup"]->Fill( NJets_pt45, eventweight*PileupWeight );
                
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)+"_qcd_up"]->Fill( NJets_pt45, eventweight*scaleWeightUp );
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)+"_pdf_up"]->Fill( NJets_pt45, eventweight*PDFWeightUp );
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)+"_pup_up"]->Fill( NJets_pt45, eventweight*PileupWeightUp );
    
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)+"_qcd_down"]->Fill( NJets_pt45, eventweight*scaleWeightDown );
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)+"_pdf_down"]->Fill( NJets_pt45, eventweight*PDFWeightDown );
                    my_histos["h_njets_pt45_0l_deepESMbin"+std::to_string(index+1)+"_pup_down"]->Fill( NJets_pt45, eventweight*PileupWeightDown );
                }
            }
        }

        double lepWeight, lepWeightUp, lepWeightDown = 0.0;
        if( NGoodElectrons == 1 ) 
        {
            lepWeight       = eleLepWeight;
            lepWeightUp     = eleLepWeightUp;
            lepWeightDown   = eleLepWeightDown;
        }
        if( NGoodMuons == 1 ) 
        {
            lepWeight       = muLepWeight;
            lepWeightUp     = muLepWeightUp;
            lepWeightDown   = muLepWeightDown;
        }

        if( passBaseline1l ) 
        {            
            my_histos["h_njets_pt30_1l"]->Fill( NJets_pt30, eventweight );
            my_histos["h_ht_pt30_1l"]->Fill( HT, eventweight );
            my_histos["h_madht_pt30_1l"]->Fill( MadHT, eventweight );
            my_histos["h_njets_pt30_1l_qcd"]->Fill( NJets_pt30, eventweight*scaleWeight );
            my_histos["h_njets_pt30_1l_pdf"]->Fill( NJets_pt30, eventweight*PDFWeight );
            my_histos["h_njets_pt30_1l_pup"]->Fill( NJets_pt30, eventweight*PileupWeight );
            my_histos["h_njets_pt30_1l_lep"]->Fill( NJets_pt30, eventweight*lepWeight );
            
            my_histos["h_njets_pt30_1l_qcd_up"]->Fill( NJets_pt30, eventweight*scaleWeightUp );
            my_histos["h_njets_pt30_1l_pdf_up"]->Fill( NJets_pt30, eventweight*PDFWeightUp );
            my_histos["h_njets_pt30_1l_pup_up"]->Fill( NJets_pt30, eventweight*PileupWeightUp );
            my_histos["h_njets_pt30_1l_lep_up"]->Fill( NJets_pt30, eventweight*lepWeightUp );
            
            my_histos["h_njets_pt30_1l_qcd_down"]->Fill( NJets_pt30, eventweight*scaleWeightDown );
            my_histos["h_njets_pt30_1l_pdf_down"]->Fill( NJets_pt30, eventweight*PDFWeightDown );
            my_histos["h_njets_pt30_1l_pup_down"]->Fill( NJets_pt30, eventweight*PileupWeightDown );
            my_histos["h_njets_pt30_1l_lep_down"]->Fill( NJets_pt30, eventweight*lepWeightDown );

            for(int index = 0; index < deepESMbins.size(); index++)
            {
                if( deepESMbins.at(index) )
                {
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)       ]->Fill( NJets_pt30, eventweight );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_qcd"]->Fill( NJets_pt30, eventweight*scaleWeight );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_pdf"]->Fill( NJets_pt30, eventweight*PDFWeight );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_pup"]->Fill( NJets_pt30, eventweight*PileupWeight );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_lep"]->Fill( NJets_pt30, eventweight*lepWeight );
                
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_qcd_up"]->Fill( NJets_pt30, eventweight*scaleWeightUp );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_pdf_up"]->Fill( NJets_pt30, eventweight*PDFWeightUp );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_pup_up"]->Fill( NJets_pt30, eventweight*PileupWeightUp );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_lep_up"]->Fill( NJets_pt30, eventweight*lepWeightUp );
                
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_qcd_down"]->Fill( NJets_pt30, eventweight*scaleWeightDown );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_pdf_down"]->Fill( NJets_pt30, eventweight*PDFWeightDown );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_pup_down"]->Fill( NJets_pt30, eventweight*PileupWeightDown );
                    my_histos["h_njets_pt30_1l_deepESMbin"+std::to_string(index+1)+"_lep_down"]->Fill( NJets_pt30, eventweight*lepWeightDown );                
                }
            }
        }
    }//END of while tr.getNextEvent loop   
}//END of function
      
void MakeNJetDists::WriteHistos( TFile* outfile ) 
{
    outfile->cd();

    for(const auto& p : my_histos) 
    {
        p.second->Write();
    }
    
    for(const auto& p : my_2d_histos) 
    {
        p.second->Write();
    }
    
    for(const auto& p : my_efficiencies) 
    {
        p.second->Write();
    }
}
