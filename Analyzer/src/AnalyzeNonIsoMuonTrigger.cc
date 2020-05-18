#define AnalyzeNonIsoMuonTrigger_cxx
#include "Analyzer/Analyzer/include/AnalyzeNonIsoMuonTrigger.h"
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

AnalyzeNonIsoMuonTrigger::AnalyzeNonIsoMuonTrigger()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeNonIsoMuonTrigger::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define some strings that are used for different scenarios that we want to calculate trigger efficiencies for
    std::vector<std::string> effTags        { "den", "num" }; 
    std::vector<std::string> trigTags       { "trig", "noTrig" }; //Require other lepton trigger? Used for studies, but conclusion was that it does not affect statistics much so it's okay
    std::vector<std::string> nJetCutTags    { "1", "2", "3", "4", "5", "6", "7" }; //Cutting at 6 jets really killed statistics, but the 1 jet cut region is not representative of our signal region, so look at all intermediate steps as well
    std::vector<std::string> htWeightTags   { "htWeight", "noHtWeight" };
   
 
    //Define binning for the histograms
    const Int_t nPtBins = 2;
    Double_t ptBinEdges[ nPtBins + 1 ] = { 50, 85, 200 };
    const Int_t nEtaBins = 2;
    Double_t etaBinEdges[ nEtaBins + 1 ] = { 0, 1.4, 2.4 };

    for( std::string effTag : effTags ) {
        for( std::string trigTag : trigTags ) {
            for( std::string nJetCutTag : nJetCutTags ) {
                for( std::string htWeightTag : htWeightTags ) {
                    //Define 1D histograms - use these to get event yields ( a little more convenient than when extracting it from the TEfficiency )
                    my_histos.emplace( "h_trig_"+effTag+"_nimu_pt40_"+trigTag+"_"+nJetCutTag+"jCut_"+htWeightTag+"_wLepPtBin", std::make_shared<TH1D>( ( "h_trig_"+effTag+"_nimu_pt40_"+trigTag+"_"+nJetCutTag+"jCut_"+htWeightTag+"_wLepPtBin" ).c_str(), ( "h_trig_"+effTag+"_nimu_pt40_"+trigTag+"_"+nJetCutTag+"jCut_"+htWeightTag+"_wLepPtBin" ).c_str(), nPtBins, ptBinEdges ) );
                    my_histos.emplace( "h_trig_"+effTag+"_nimu_pt40_"+trigTag+"_"+nJetCutTag+"jCut_"+htWeightTag+"_wLepEtaBin", std::make_shared<TH1D>( ( "h_trig_"+effTag+"_nimu_pt40_"+trigTag+"_"+nJetCutTag+"jCut_"+htWeightTag+"_wLepEtaBin" ).c_str(), ( "h_trig_"+effTag+"_nimu_pt40_"+trigTag+"_"+nJetCutTag+"jCut_"+htWeightTag+"_wLepEtaBin" ).c_str(), nEtaBins, etaBinEdges ) );
    
                    //Define 2D histograms - use these for the final plots ( i.e., it is easier to get these values and use these to create TEfficiencies later on )
                    my_2d_histos.emplace( "h2_trig_"+effTag+"_nimu_pt40_"+trigTag+"_"+nJetCutTag+"jCut_"+htWeightTag+"_wLepPtLepEtaBin", std::make_shared<TH2D>( ( "h2_trig_"+effTag+"_nimu_pt40_"+trigTag+"_"+nJetCutTag+"jCut_"+htWeightTag+"_wLepPtLepEtaBin" ).c_str(), ( "h2_trig_"+effTag+"_nimu_pt40_"+trigTag+"_"+nJetCutTag+"jCut_"+htWeightTag+"_wLepPtLepEtaBin" ).c_str(), nPtBins, ptBinEdges, nEtaBins, etaBinEdges ) );
                }//End of htWeightTags loop
            }//End of nJetCutTags loop
        }//End of trigTags loop
    }//End of effTags loop
}

//Put everything you want to do per event here.
void AnalyzeNonIsoMuonTrigger::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter            = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //Define useful variables here
        const auto& TriggerNames            = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass             = tr.getVec<int>("TriggerPass");

        const auto& runtype                 = tr.getVar<std::string>("runtype");
        const auto& runYear                 = tr.getVar<std::string>("runYear");
        const auto& passMadHT               = tr.getVar<bool>("passMadHT");

        const auto& JetID                   = tr.getVar<bool>("JetID");
        const auto& HT_NonIsoMuon_pt30      = tr.getVar<double>("HT_NonIsoMuon_pt30");
        const auto& NNonIsoMuonJets_pt30    = tr.getVar<int>("NNonIsoMuonJets_pt30");

        const auto& Muons                   = tr.getVec<TLorentzVector>("Muons");
        const auto& Electrons               = tr.getVec<TLorentzVector>("Electrons");
        
        const auto& GoodElectrons           = tr.getVec<bool>("GoodElectrons");
        const auto& NonIsoMuons             = tr.getVec<bool>("NonIsoMuons");

        const auto& passHEMVeto             = tr.getVar<bool>("passHEMVeto");
        const auto& correct2018Split        = tr.getVar<bool>("correct2018Split");
        const auto& passMETFilters          = tr.getVar<bool>("passMETFilters");
       
        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        //--------------------------------------------------
        //-- Print list of triggers (only if you want to see them)
        //--------------------------------------------------

        //if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

        // ------------------------
        // -- Define weight
        // ------------------------
        double eventweight          = 1.0;
        double hteventweight        = 1.0;
        
        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight              = tr.getVar<double>("Weight");
            const auto& Lumi                = tr.getVar<double>("Lumi");
            const auto& totNonIsoMuonSF     = tr.getVar<double>("totNonIsoMuonSF");
            const auto& puWeight            = tr.getVar<double>("puWeightCorr");
            const auto& htWeight            = tr.getVar<double>("htDerivedweight");

            eventweight = Weight*Lumi*totNonIsoMuonSF*puWeight;
            hteventweight = eventweight*htWeight;
        }
       
        //Define trigger cuts
        std::vector<std::string> myNonIsoMuonTriggers       { "HLT_Mu50_v", "HLT_TkMu50_v" };
        std::vector<std::string> myElectronTriggers         { "HLT_Ele27_WPTight_Gsf_v", "HLT_Photon175_v", "HLT_Ele115_CaloIdVT_GsfTrkIdT_v" };

        if( ( runYear == "2017" ) || ( runYear == "2018" ) ) {
            std::vector<std::string> myElectronTriggers     { "HLT_Ele35_WPTight_Gsf", "HLT_Photon200", "HLT_Ele115CaloIdVT_GsfTrkIdT" };
        }

        bool passNonIsoMuonTriggers = passTriggerGeneral( myNonIsoMuonTriggers, TriggerNames, TriggerPass );
        bool passElectronTriggers   = passTriggerGeneral( myElectronTriggers, TriggerNames, TriggerPass );
        
        bool pass7JetCut            = ( NNonIsoMuonJets_pt30 >= 7 );
        bool pass6JetCut            = ( NNonIsoMuonJets_pt30 >= 6 );
        bool pass5JetCut            = ( NNonIsoMuonJets_pt30 >= 5 );
        bool pass4JetCut            = ( NNonIsoMuonJets_pt30 >= 4 );
        bool pass3JetCut            = ( NNonIsoMuonJets_pt30 >= 3 );
        bool pass2JetCut            = ( NNonIsoMuonJets_pt30 >= 2 );
        bool pass1JetCut            = ( NNonIsoMuonJets_pt30 >= 1 ); //A little redundant because we require one good b jet but for symmetry reasons

        bool atLeastOneGoodEl       = ( std::any_of( GoodElectrons.begin(), GoodElectrons.end(), [] ( bool boolEl ) { return boolEl; } ) );

        //Define offfline baseline for the qcd control region
        bool passNonIsoMuonOfflineBaseline = ( JetID && passMadHT && HT_NonIsoMuon_pt30 > 300 && passHEMVeto && correct2018Split && passMETFilters );

        if( !passNonIsoMuonOfflineBaseline || !atLeastOneGoodEl ) continue;

        bool foundElectronPt40      = false;
            
        //Require a good electron in the single electron dataset
        for( unsigned int itEl = 0; itEl < Electrons.size(); ++itEl ) {
            if( !GoodElectrons.at( itEl ) ) continue; 

            TLorentzVector myElectron = Electrons.at( itEl );
            
            if( myElectron.Pt() >= 40 && std::fabs( myElectron.Eta() ) < 2.4 )       foundElectronPt40 = true;
        }

        //Look at the first non iso muon
        int myNonIsoMuonIndex = -1;
        for( unsigned int itNIMu = 0; itNIMu < Muons.size(); ++itNIMu ) {
            if( !NonIsoMuons.at( itNIMu ) ) continue;

            myNonIsoMuonIndex = itNIMu;
            //std::cout<<"Non Iso Muon Pt: "<<Muons.at(itNIMu).Pt()<<"; Eta: "<<Muons.at(itNIMu).Eta()<<std::endl;
            break;
        }

        if( myNonIsoMuonIndex != -1 ) {
                    
            const std::map<std::string, bool> cut_map_Triggers {
    
                {   "pt40_trig_1jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass1JetCut }, 
                {   "pt40_trig_2jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass2JetCut }, 
                {   "pt40_trig_3jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass3JetCut }, 
                {   "pt40_trig_4jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass4JetCut }, 
                {   "pt40_trig_5jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass5JetCut }, 
                {   "pt40_trig_6jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass6JetCut }, 
                {   "pt40_trig_7jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && passElectronTriggers && pass7JetCut }, 
                    
                {   "pt40_noTrig_1jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && pass1JetCut }, 
                {   "pt40_noTrig_2jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && pass2JetCut }, 
                {   "pt40_noTrig_3jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && pass3JetCut }, 
                {   "pt40_noTrig_4jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && pass4JetCut }, 
                {   "pt40_noTrig_5jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && pass5JetCut }, 
                {   "pt40_noTrig_6jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && pass6JetCut },
                {   "pt40_noTrig_7jCut",       passNonIsoMuonOfflineBaseline && foundElectronPt40 && pass7JetCut }
                
            };
                
            for( auto& kv : cut_map_Triggers ) {
                if( kv.second ) {
                    my_histos["h_trig_den_nimu_"+kv.first+"_htWeight_wLepPtBin"]->Fill( Muons.at( myNonIsoMuonIndex ).Pt(), hteventweight );
                    my_histos["h_trig_den_nimu_"+kv.first+"_htWeight_wLepEtaBin"]->Fill( std::fabs( Muons.at( myNonIsoMuonIndex ).Eta() ), hteventweight );
                    my_2d_histos["h2_trig_den_nimu_"+kv.first+"_htWeight_wLepPtLepEtaBin"]->Fill( Muons.at( myNonIsoMuonIndex ).Pt(), std::fabs( Muons.at( myNonIsoMuonIndex ).Eta() ), hteventweight );
                    
                    my_histos["h_trig_den_nimu_"+kv.first+"_noHtWeight_wLepPtBin"]->Fill( Muons.at( myNonIsoMuonIndex ).Pt(), eventweight );
                    my_histos["h_trig_den_nimu_"+kv.first+"_noHtWeight_wLepEtaBin"]->Fill( std::fabs( Muons.at( myNonIsoMuonIndex ).Eta() ), eventweight );
                    my_2d_histos["h2_trig_den_nimu_"+kv.first+"_noHtWeight_wLepPtLepEtaBin"]->Fill( Muons.at( myNonIsoMuonIndex ).Pt(), std::fabs( Muons.at( myNonIsoMuonIndex ).Eta() ), eventweight );
                    
                    if( passNonIsoMuonTriggers ) {
                        my_histos["h_trig_num_nimu_"+kv.first+"_htWeight_wLepPtBin"]->Fill( Muons.at( myNonIsoMuonIndex ).Pt(), hteventweight );
                        my_histos["h_trig_num_nimu_"+kv.first+"_htWeight_wLepEtaBin"]->Fill( std::fabs( Muons.at( myNonIsoMuonIndex ).Eta() ), hteventweight );
                        my_2d_histos["h2_trig_num_nimu_"+kv.first+"_htWeight_wLepPtLepEtaBin"]->Fill( Muons.at( myNonIsoMuonIndex ).Pt(), std::fabs( Muons.at( myNonIsoMuonIndex ).Eta() ), hteventweight );
                        
                        my_histos["h_trig_num_nimu_"+kv.first+"_noHtWeight_wLepPtBin"]->Fill( Muons.at( myNonIsoMuonIndex ).Pt(), eventweight );
                        my_histos["h_trig_num_nimu_"+kv.first+"_noHtWeight_wLepEtaBin"]->Fill( std::fabs( Muons.at( myNonIsoMuonIndex ).Eta() ), eventweight );
                        my_2d_histos["h2_trig_num_nimu_"+kv.first+"_noHtWeight_wLepPtLepEtaBin"]->Fill( Muons.at( myNonIsoMuonIndex ).Pt(), std::fabs( Muons.at( myNonIsoMuonIndex ).Eta() ), eventweight );
                    }
                }
            }
        }
    } 
}

void AnalyzeNonIsoMuonTrigger::WriteHistos(TFile* outfile)
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
bool AnalyzeNonIsoMuonTrigger::passTriggerGeneral( std::vector<std::string>& myTriggerVector, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass ) {
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
void AnalyzeNonIsoMuonTrigger::printTriggerList( const std::vector<std::string>& TriggerNames ) {
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) {
        std::string myString = TriggerNames.at(i);
        printf("%s\n", myString.c_str());
    }
}

//Use this function if you want to make sure a specific trigger exists (sometimes the version number can cause some issues)
void AnalyzeNonIsoMuonTrigger::doesTriggerExist( std::vector<std::string>& myTrigTestVector, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass ) {
    if( passTriggerGeneral( myTrigTestVector, TriggerNames, TriggerPass ) ) printf( "%s\n", "This trigger works and the string comparison comes out valid." );
}
