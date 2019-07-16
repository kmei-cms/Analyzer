#define AnalyzeLepTrigger_cxx
#include "Analyzer/Analyzer/include/AnalyzeLepTrigger.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <TFile.h>

AnalyzeLepTrigger::AnalyzeLepTrigger()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeLepTrigger::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define some strings that are used for different scenarios that we want to calculate trigger efficiencies for
    std::vector<std::string> effTags        { "den", "num" }; 
    std::vector<std::string> lepTags        { "el", "mu" }; //Electron, muon, non-iso muon
    std::vector<std::string> ptTags         { "pt35", "pt40" }; //Pt threshold of other lepton
    std::vector<std::string> trigTags       { "trig", "noTrig" }; //Require other lepton trigger? Used for studies, but conclusion was that it does not affect statistics much so it's okay
    std::vector<std::string> nJetCutTags    { "1", "2", "3", "4", "5", "6" }; //Cutting at 6 jets really killed statistics, but the 1 jet cut region is not representative of our signal region, so look at all intermediate steps as well
   
 
    //Define binning for the histograms
    const Int_t nPtBins = 6;
    Double_t ptBinEdges[ nPtBins + 1 ] = { 30, 40, 50, 60, 70, 120, 200 };
    const Int_t nEtaBins = 6;
    Double_t etaBinEdges[ nEtaBins + 1 ] = { -2.4, -1.8, -1.0, 0, 1.0, 1.8, 2.4 };

    for( std::string effTag : effTags ) {
        for( std::string lepTag : lepTags ) {
            for( std::string ptTag : ptTags ) {
                if( lepTag == "mu" && ptTag == "pt35" ) continue;
                for( std::string trigTag : trigTags ) {
                    for( std::string nJetCutTag : nJetCutTags ) {
                        
                        //Define 1D histograms - use these to get event yields ( a little more convenient than when extracting it from the TEfficiency )
                        my_histos.emplace( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtBin", std::make_shared<TH1D>( ( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtBin" ).c_str(), ( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtBin" ).c_str(), nPtBins, ptBinEdges ) );
                        my_histos.emplace( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepEtaBin", std::make_shared<TH1D>( ( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepEtaBin" ).c_str(), ( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepEtaBin" ).c_str(), nEtaBins, etaBinEdges ) );

                        //Define 2D histograms - use these for the final plots ( i.e., it is easier to get these values and use these to create TEfficiencies later on )
                        my_2d_histos.emplace( "h2_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin", std::make_shared<TH2D>( ( "h2_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin" ).c_str(), ( "h2_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin" ).c_str(), nPtBins, ptBinEdges, nEtaBins, etaBinEdges ) );

                        //Define TEfficincies - use these for sanity checks
//                        my_efficiencies.emplace( "h_eff_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtBin", std::make_shared<TEfficiency>( ( "h_eff_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtBin" ).c_str(), ( "h_eff_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtBin" ).c_str(), nPtBins, ptBinEdges ) );
//                        my_efficiencies.emplace( "h_eff_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepEtaBin", std::make_shared<TEfficiency>( ( "h_eff_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepEtaBin" ).c_str(), ( "h_eff_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepEtaBin" ).c_str(), nEtaBins, etaBinEdges ) );
//                        my_efficiencies.emplace( "h2_eff_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin", std::make_shared<TEfficiency>( ( "h2_eff_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin" ).c_str(), ( "h2_eff_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin" ).c_str(), nPtBins, ptBinEdges, nEtaBins, etaBinEdges ) );
                    }//End of nJetCutTags loop
                }//End of trigTags loop
            }//End of ptTags loop
        }//End of lepTags loop
    }//End of effTags loop

}

//Put everything you want to do per event here.
void AnalyzeLepTrigger::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );


        //Define useful variables here
        const auto& TriggerNames        = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass         = tr.getVec<int>("TriggerPass");

        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& runYear             = tr.getVar<std::string>("runYear");
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");

        const auto& JetID               = tr.getVar<bool>("JetID");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& passTriggerMC       = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodBJets_pt30     = tr.getVar<int>("NGoodBJets_pt30");
        const auto& Mbl                 = tr.getVar<double>("Mbl");
        const auto& HT_trigger_pt30     = tr.getVar<double>("HT_trigger_pt30");
        const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30");

        const auto& Muons               = tr.getVec<TLorentzVector>("Muons");
        const auto& Electrons           = tr.getVec<TLorentzVector>("Electrons");
        
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passBaseline        = tr.getVar<bool>("passBaseline1l_Good");

        const auto& GoodMuons           = tr.getVec<bool>("GoodMuons");
        const auto& GoodElectrons       = tr.getVec<bool>("GoodElectrons");
       
        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;


        //--------------------------------------------------
        //-- Print list of triggers (only if you want to see them)
        //--------------------------------------------------

        //if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double htDerivedScaleFactor = 1.0;
        double topPtScaleFactor     = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        
        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            bTagScaleFactor   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
            topPtScaleFactor = tr.getVar<double>("topPtScaleFactor");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor = tr.getVar<double>("puWeightCorr");
            
            weight *= eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;
        }
       
        //Define trigger cuts
        std::vector<std::string> myMuonTriggers             { "HLT_IsoMu24_v", "HLT_IsoTkMu24_v", "HLT_Mu50_v", "HLT_TkMu50_v" };
        std::vector<std::string> myElectronTriggers         { "HLT_Ele27_WPTight_Gsf_v", "HLT_Photon175_v", "HLT_Ele115_CaloIdVT_GsfTrkIdT_v" };


        if( runYear == "2017" ) {
            std::vector<std::string> myMuonTriggers         { "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_TkMu50" };
            std::vector<std::string> myElectronTriggers     { "HLT_Ele35_WPTight_Gsf", "HLT_Photon200", "HLT_Ele115CaloIdVT_GsfTrkIdT" };
        }

        bool passMuonTriggers       = passTriggerGeneral( myMuonTriggers, TriggerNames, TriggerPass );
        bool passElectronTriggers   = passTriggerGeneral( myElectronTriggers, TriggerNames, TriggerPass );
        
        //Define offline baseline for one lepton analysis
        bool passOfflineBaseline    = ( JetID && passMadHT && NGoodBJets_pt30 >= 1 && Mbl > 50 && Mbl < 250 && HT_trigger_pt30 > 300 );

        bool pass6JetCut            = ( NGoodJets_pt30 >= 6 );
        bool pass5JetCut            = ( NGoodJets_pt30 >= 5 );
        bool pass4JetCut            = ( NGoodJets_pt30 >= 4 );
        bool pass3JetCut            = ( NGoodJets_pt30 >= 3 );
        bool pass2JetCut            = ( NGoodJets_pt30 >= 2 );
        bool pass1JetCut            = ( NGoodJets_pt30 >= 1 ); //A little redundant because we require one good b jet but for symmetry reasons

        bool atLeastOneGoodEl       = ( std::any_of( GoodElectrons.begin(), GoodElectrons.end(), [] ( bool boolEl ) { return boolEl; } ) );
        bool atLeastOneGoodMu       = ( std::any_of( GoodMuons.begin(), GoodMuons.end(), [] ( bool boolMu ) { return boolMu; } ) );


        //------------------------------------------------------
        //-- Do the Electron Trigger Efficiency on the Single Muon Dataset
        //------------------------------------------------------
        
        if( filetag == "2016_Data_SingleMuon" || filetag == "2017_Data_SingleMuon" ) {
            
            if( !passOfflineBaseline || !atLeastOneGoodMu ) continue;

            bool foundMuonPt35      = false;
            bool foundMuonPt40      = false;
            
            //Require a good muon in the single muon dataset
            for( unsigned int itMu = 0; itMu < Muons.size(); ++itMu ) {
                if( !GoodMuons.at( itMu ) ) continue; 

                TLorentzVector myMuon = Muons.at( itMu );

                if( myMuon.Pt() >= 35 && std::fabs( myMuon.Eta() ) < 2.4 )       foundMuonPt35 = true;
                if( myMuon.Pt() >= 40 && std::fabs( myMuon.Eta() ) < 2.4 )       foundMuonPt40 = true;
            }

            //Look at the first good electron
            int myGoodElectronIndex     = -1;
            for( unsigned int itEl = 0; itEl < Electrons.size(); ++itEl ) {
                if( !GoodElectrons.at( itEl ) ) continue;

                myGoodElectronIndex = itEl;
                break;
            }
       
            if( myGoodElectronIndex != -1 ) {
                const std::map<std::string, bool> cut_map_elTriggers {
                    {   "el_pt35_trig_1jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass1JetCut }, 
                    {   "el_pt35_trig_2jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass2JetCut }, 
                    {   "el_pt35_trig_3jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass3JetCut }, 
                    {   "el_pt35_trig_4jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass4JetCut }, 
                    {   "el_pt35_trig_5jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass5JetCut }, 
                    {   "el_pt35_trig_6jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass6JetCut }, 
                    
                    {   "el_pt35_noTrig_1jCut",       passOfflineBaseline && foundMuonPt35 && pass1JetCut }, 
                    {   "el_pt35_noTrig_2jCut",       passOfflineBaseline && foundMuonPt35 && pass2JetCut }, 
                    {   "el_pt35_noTrig_3jCut",       passOfflineBaseline && foundMuonPt35 && pass3JetCut }, 
                    {   "el_pt35_noTrig_4jCut",       passOfflineBaseline && foundMuonPt35 && pass4JetCut }, 
                    {   "el_pt35_noTrig_5jCut",       passOfflineBaseline && foundMuonPt35 && pass5JetCut }, 
                    {   "el_pt35_noTrig_6jCut",       passOfflineBaseline && foundMuonPt35 && pass6JetCut }, 
                    
                    {   "el_pt40_trig_1jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass1JetCut }, 
                    {   "el_pt40_trig_2jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass2JetCut }, 
                    {   "el_pt40_trig_3jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass3JetCut }, 
                    {   "el_pt40_trig_4jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass4JetCut }, 
                    {   "el_pt40_trig_5jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass5JetCut }, 
                    {   "el_pt40_trig_6jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass6JetCut }, 
                    
                    {   "el_pt40_noTrig_1jCut",       passOfflineBaseline && foundMuonPt40 && pass1JetCut }, 
                    {   "el_pt40_noTrig_2jCut",       passOfflineBaseline && foundMuonPt40 && pass2JetCut }, 
                    {   "el_pt40_noTrig_3jCut",       passOfflineBaseline && foundMuonPt40 && pass3JetCut }, 
                    {   "el_pt40_noTrig_4jCut",       passOfflineBaseline && foundMuonPt40 && pass4JetCut }, 
                    {   "el_pt40_noTrig_5jCut",       passOfflineBaseline && foundMuonPt40 && pass5JetCut }, 
                    {   "el_pt40_noTrig_6jCut",       passOfflineBaseline && foundMuonPt40 && pass6JetCut } 
                };
    
                for( auto& kv : cut_map_elTriggers ) {
                    if( kv.second ) {
                        my_histos["h_trig_den_"+kv.first+"_wLepPtBin"]->Fill( Electrons.at( myGoodElectronIndex ).Pt(), eventweight );
                        my_histos["h_trig_den_"+kv.first+"_wLepEtaBin"]->Fill( Electrons.at( myGoodElectronIndex ).Eta(), eventweight );
                        my_2d_histos["h2_trig_den_"+kv.first+"_wLepPtLepEtaBin"]->Fill( Electrons.at( myGoodElectronIndex ).Pt(), Electrons.at( myGoodElectronIndex ).Eta(), eventweight );
    
                        if( passElectronTriggers ) {
                            my_histos["h_trig_num_"+kv.first+"_wLepPtBin"]->Fill( Electrons.at( myGoodElectronIndex ).Pt(), eventweight );
                            my_histos["h_trig_num_"+kv.first+"_wLepEtaBin"]->Fill( Electrons.at( myGoodElectronIndex ).Eta(), eventweight );
                            my_2d_histos["h2_trig_num_"+kv.first+"_wLepPtLepEtaBin"]->Fill( Electrons.at( myGoodElectronIndex ).Pt(), Electrons.at( myGoodElectronIndex ).Eta(), eventweight );
                        }
                    }
                }
            }
        }
        
        //------------------------------------------------------
        //-- Do the Muon Trigger Efficiency on the Single Electron Dataset
        //------------------------------------------------------
        
        if( filetag == "2016_Data_SingleElectron" || filetag == "2017_Data_SingleElectron" ) {
            
            if( !( passOfflineBaseline || passNonIsoMuonOfflineBaseline )  || !atLeastOneGoodEl ) continue;

            bool foundElectronPt40      = false;
            
            //Require a good electron in the single electron dataset
            for( unsigned int itEl = 0; itEl < Electrons.size(); ++itEl ) {
                if( !GoodElectrons.at( itEl ) ) continue; 

                TLorentzVector myElectron = Electrons.at( itEl );

                if( myElectron.Pt() >= 40 && std::fabs( myElectron.Eta() ) < 2.4 )       foundElectronPt40 = true;
            }

            //Look at the first good muon
            int myGoodMuonIndex     = -1;
            for( unsigned int itMu = 0; itMu < Muons.size(); ++itMu ) {
                if( !GoodMuons.at( itMu ) ) continue;

                myGoodMuonIndex = itMu;
                std::cout<<"Iso Muon Pt: "<<Muons.at(itMu).Pt()<<"; Eta: "<<Muons.at(itMu).Eta()<<std::endl;
                break;
            }
            
            if( myGoodMuonIndex != -1 ) {
        
                const std::map<std::string, bool> cut_map_muTriggers {
                    
                    {   "mu_pt40_trig_1jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass1JetCut }, 
                    {   "mu_pt40_trig_2jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass2JetCut }, 
                    {   "mu_pt40_trig_3jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass3JetCut }, 
                    {   "mu_pt40_trig_4jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass4JetCut }, 
                    {   "mu_pt40_trig_5jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass5JetCut }, 
                    {   "mu_pt40_trig_6jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass6JetCut }, 
                    
                    {   "mu_pt40_noTrig_1jCut",       passOfflineBaseline && foundElectronPt40 && pass1JetCut }, 
                    {   "mu_pt40_noTrig_2jCut",       passOfflineBaseline && foundElectronPt40 && pass2JetCut }, 
                    {   "mu_pt40_noTrig_3jCut",       passOfflineBaseline && foundElectronPt40 && pass3JetCut }, 
                    {   "mu_pt40_noTrig_4jCut",       passOfflineBaseline && foundElectronPt40 && pass4JetCut }, 
                    {   "mu_pt40_noTrig_5jCut",       passOfflineBaseline && foundElectronPt40 && pass5JetCut }, 
                    {   "mu_pt40_noTrig_6jCut",       passOfflineBaseline && foundElectronPt40 && pass6JetCut }, 
                    
                };
            
                for( auto& kv : cut_map_muTriggers ) {
                    if( kv.second ) {
                        
                        my_histos["h_trig_den_"+kv.first+"_wLepPtBin"]->Fill( Muons.at( myGoodMuonIndex ).Pt(), eventweight );
                        my_histos["h_trig_den_"+kv.first+"_wLepEtaBin"]->Fill( Muons.at( myGoodMuonIndex ).Eta(), eventweight );
                        my_2d_histos["h2_trig_den_"+kv.first+"_wLepPtLepEtaBin"]->Fill( Muons.at( myGoodMuonIndex ).Pt(), Muons.at( myGoodMuonIndex ).Eta(), eventweight );
                        
                        if( passMuonTriggers ) {
                            my_histos["h_trig_num_"+kv.first+"_wLepPtBin"]->Fill( Muons.at( myGoodMuonIndex ).Pt(), eventweight );
                            my_histos["h_trig_num_"+kv.first+"_wLepEtaBin"]->Fill( Muons.at( myGoodMuonIndex ).Eta(), eventweight );
                            my_2d_histos["h2_trig_num_"+kv.first+"_wLepPtLepEtaBin"]->Fill( Muons.at( myGoodMuonIndex ).Pt(), Muons.at( myGoodMuonIndex ).Eta(), eventweight );
                        }
                    }
                }
            }
        }
        
        //------------------------------------------------------
        //-- Do both trigger efficiencies on the MC
        //------------------------------------------------------
        
        if( runtype == "MC" ) {
          
            if( atLeastOneGoodEl ) {
                //Do the muon studies
                bool foundElectronPt40      = false;
                
                //Require a good electron in the single electron dataset
                for( unsigned int itEl = 0; itEl < Electrons.size(); ++itEl ) {
                    if( !GoodElectrons.at( itEl ) ) continue; 
    
                    TLorentzVector myElectron = Electrons.at( itEl );
    
                    if( myElectron.Pt() >= 40 && std::fabs( myElectron.Eta() ) < 2.4 )       foundElectronPt40 = true;
                }
    
                //Look at the first good muon
                int myGoodMuonIndex     = -1;
                for( unsigned int itMu = 0; itMu < Muons.size(); ++itMu ) {
                    if( !GoodMuons.at( itMu ) ) continue;
    
                    myGoodMuonIndex = itMu;
                    std::cout<<"MC Iso Muon Pt: "<<Muons.at(itMu).Pt()<<"; Eta: "<<Muons.at(itMu).Eta()<<std::endl;
                    break;
                }
                
                if( myGoodMuonIndex != -1 ) { 
                    const std::map<std::string, bool> cut_map_muTriggers {
                        
                        {   "mu_pt40_trig_1jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass1JetCut }, 
                        {   "mu_pt40_trig_2jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass2JetCut }, 
                        {   "mu_pt40_trig_3jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass3JetCut }, 
                        {   "mu_pt40_trig_4jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass4JetCut }, 
                        {   "mu_pt40_trig_5jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass5JetCut }, 
                        {   "mu_pt40_trig_6jCut",       passOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass6JetCut }, 
                        
                        {   "mu_pt40_noTrig_1jCut",       passOfflineBaseline && foundElectronPt40 && pass1JetCut }, 
                        {   "mu_pt40_noTrig_2jCut",       passOfflineBaseline && foundElectronPt40 && pass2JetCut }, 
                        {   "mu_pt40_noTrig_3jCut",       passOfflineBaseline && foundElectronPt40 && pass3JetCut }, 
                        {   "mu_pt40_noTrig_4jCut",       passOfflineBaseline && foundElectronPt40 && pass4JetCut }, 
                        {   "mu_pt40_noTrig_5jCut",       passOfflineBaseline && foundElectronPt40 && pass5JetCut }, 
                        {   "mu_pt40_noTrig_6jCut",       passOfflineBaseline && foundElectronPt40 && pass6JetCut }, 
                        
                    };
                
                    for( auto& kv : cut_map_muTriggers ) {
                        if( kv.second ) {
                            
                            my_histos["h_trig_den_"+kv.first+"_wLepPtBin"]->Fill( Muons.at( myGoodMuonIndex ).Pt(), eventweight );
                            my_histos["h_trig_den_"+kv.first+"_wLepEtaBin"]->Fill( Muons.at( myGoodMuonIndex ).Eta(), eventweight );
                            my_2d_histos["h2_trig_den_"+kv.first+"_wLepPtLepEtaBin"]->Fill( Muons.at( myGoodMuonIndex ).Pt(), Muons.at( myGoodMuonIndex ).Eta(), eventweight );
                            
                            if( passMuonTriggers ) {
                                my_histos["h_trig_num_"+kv.first+"_wLepPtBin"]->Fill( Muons.at( myGoodMuonIndex ).Pt(), eventweight );
                                my_histos["h_trig_num_"+kv.first+"_wLepEtaBin"]->Fill( Muons.at( myGoodMuonIndex ).Eta(), eventweight );
                                my_2d_histos["h2_trig_num_"+kv.first+"_wLepPtLepEtaBin"]->Fill( Muons.at( myGoodMuonIndex ).Pt(), Muons.at( myGoodMuonIndex ).Eta(), eventweight );
                            }
                        }
                    }
                }
                
            }
            
            if( passOfflineBaseline && atLeastOneGoodMu ) {
    
                bool foundMuonPt35      = false;
                bool foundMuonPt40      = false;
                
                //Require a good muon in the single muon dataset
                for( unsigned int itMu = 0; itMu < Muons.size(); ++itMu ) {
                    if( !GoodMuons.at( itMu ) ) continue; 
    
                    TLorentzVector myMuon = Muons.at( itMu );
    
                    if( myMuon.Pt() >= 35 && std::fabs( myMuon.Eta() ) < 2.4 )       foundMuonPt35 = true;
                    if( myMuon.Pt() >= 40 && std::fabs( myMuon.Eta() ) < 2.4 )       foundMuonPt40 = true;
                }
    
                //Look at the first good electron
                int myGoodElectronIndex     = -1;
                for( unsigned int itEl = 0; itEl < Electrons.size(); ++itEl ) {
                    if( !GoodElectrons.at( itEl ) ) continue;
    
                    myGoodElectronIndex = itEl;
                    break;
                }
            
                if( myGoodElectronIndex != -1 ) {
                    const std::map<std::string, bool> cut_map_elTriggers {
                        {   "el_pt35_trig_1jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass1JetCut }, 
                        {   "el_pt35_trig_2jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass2JetCut }, 
                        {   "el_pt35_trig_3jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass3JetCut }, 
                        {   "el_pt35_trig_4jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass4JetCut }, 
                        {   "el_pt35_trig_5jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass5JetCut }, 
                        {   "el_pt35_trig_6jCut",       passOfflineBaseline && foundMuonPt35 && passMuonTriggers && pass6JetCut }, 
                        
                        {   "el_pt35_noTrig_1jCut",       passOfflineBaseline && foundMuonPt35 && pass1JetCut }, 
                        {   "el_pt35_noTrig_2jCut",       passOfflineBaseline && foundMuonPt35 && pass2JetCut }, 
                        {   "el_pt35_noTrig_3jCut",       passOfflineBaseline && foundMuonPt35 && pass3JetCut }, 
                        {   "el_pt35_noTrig_4jCut",       passOfflineBaseline && foundMuonPt35 && pass4JetCut }, 
                        {   "el_pt35_noTrig_5jCut",       passOfflineBaseline && foundMuonPt35 && pass5JetCut }, 
                        {   "el_pt35_noTrig_6jCut",       passOfflineBaseline && foundMuonPt35 && pass6JetCut }, 
                        
                        {   "el_pt40_trig_1jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass1JetCut }, 
                        {   "el_pt40_trig_2jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass2JetCut }, 
                        {   "el_pt40_trig_3jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass3JetCut }, 
                        {   "el_pt40_trig_4jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass4JetCut }, 
                        {   "el_pt40_trig_5jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass5JetCut }, 
                        {   "el_pt40_trig_6jCut",       passOfflineBaseline && foundMuonPt40 && passMuonTriggers && pass6JetCut }, 
                        
                        {   "el_pt40_noTrig_1jCut",       passOfflineBaseline && foundMuonPt40 && pass1JetCut }, 
                        {   "el_pt40_noTrig_2jCut",       passOfflineBaseline && foundMuonPt40 && pass2JetCut }, 
                        {   "el_pt40_noTrig_3jCut",       passOfflineBaseline && foundMuonPt40 && pass3JetCut }, 
                        {   "el_pt40_noTrig_4jCut",       passOfflineBaseline && foundMuonPt40 && pass4JetCut }, 
                        {   "el_pt40_noTrig_5jCut",       passOfflineBaseline && foundMuonPt40 && pass5JetCut }, 
                        {   "el_pt40_noTrig_6jCut",       passOfflineBaseline && foundMuonPt40 && pass6JetCut } 
                    };
        
                    for( auto& kv : cut_map_elTriggers ) {
                        if( kv.second ) {
                            std::cout<<kv.first<<std::endl;
                            my_histos["h_trig_den_"+kv.first+"_wLepPtBin"]->Fill( Electrons.at( myGoodElectronIndex ).Pt(), eventweight );
                            my_histos["h_trig_den_"+kv.first+"_wLepEtaBin"]->Fill( Electrons.at( myGoodElectronIndex ).Eta(), eventweight );
                            my_2d_histos["h2_trig_den_"+kv.first+"_wLepPtLepEtaBin"]->Fill( Electrons.at( myGoodElectronIndex ).Pt(), Electrons.at( myGoodElectronIndex ).Eta(), eventweight );
        
                            if( passElectronTriggers ) {
                                my_histos["h_trig_num_"+kv.first+"_wLepPtBin"]->Fill( Electrons.at( myGoodElectronIndex ).Pt(), eventweight );
                                my_histos["h_trig_num_"+kv.first+"_wLepEtaBin"]->Fill( Electrons.at( myGoodElectronIndex ).Eta(), eventweight );
                                my_2d_histos["h2_trig_num_"+kv.first+"_wLepPtLepEtaBin"]->Fill( Electrons.at( myGoodElectronIndex ).Pt(), Electrons.at( myGoodElectronIndex ).Eta(), eventweight );
                            }
                        }
                    }
                }
            }
        }
    } 
}

void AnalyzeLepTrigger::WriteHistos(TFile* outfile)
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
    
}

//This is just the general trigger passing function
bool AnalyzeLepTrigger::passTriggerGeneral( std::vector<std::string>& myTriggerVector, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass ) {
    bool passTrigger = false;
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) {
        if( TriggerPass.at(i) != 1 ) continue;
        std::string trigname = TriggerNames.at(i);
        
        //Now comes the fun bit logic
        if( std::any_of( myTriggerVector.begin(), myTriggerVector.end(), [&] (std::string s) { return trigname.find(s) != std::string::npos; } ) ) {
            passTrigger = true;
            break;
        }
    }
    return passTrigger;
}

//Use this function to print out the entire trigger list in an ntuple (useful when trying to figure out which triggers to use)
void AnalyzeLepTrigger::printTriggerList( const std::vector<std::string>& TriggerNames ) {
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) {
        std::string myString = TriggerNames.at(i);
        printf("%s\n", myString.c_str());
    }
}

//Use this function if you want to make sure a specific trigger exists (sometimes the version number can cause some issues)
void AnalyzeLepTrigger::doesTriggerExist( std::vector<std::string>& myTrigTestVector, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass ) {
    if( passTriggerGeneral( myTrigTestVector, TriggerNames, TriggerPass ) ) printf( "%s\n", "This trigger works and the string comparison comes out valid." );
}
