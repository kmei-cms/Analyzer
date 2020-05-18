#define MakeNewQCDSystematic_cxx
#include "Analyzer/Analyzer/include/MakeNewQCDSystematic.h"
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

std::vector<std::string> leptonTags { "nonIsoMuon", "isoMuon", "isoElectron" };
std::vector<std::string> njetTags   { "7", "8", "9", "10", "11" };
std::vector<std::string> mvaTags    { "D1", "D2", "D3", "D4" };
std::vector<std::string> nBinTags   { "N1", "N2", "N3", "N4" };
std::vector<std::string> nomWeights { "noWghts", "allWghts" };
std::vector<std::string> pdfWeights { "noPdf", "pdfUp", "pdfDown" };
std::vector<std::string> sclWeights { "noScl", "sclUp", "sclDown" };
std::vector<std::string> isrWeights { "noIsr", "isrUp", "isrDown" };
std::vector<std::string> fsrWeights { "noFsr", "fsrUp", "fsrDown" };

MakeNewQCDSystematic::MakeNewQCDSystematic()
{
    InitHistos();
}

//Define all your histograms here. 
void MakeNewQCDSystematic::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define 1D histograms
    
    //Make plots of mva shapes per njet for electrons, muons, and non isolated muons
   
    for( std::string leptonTag : leptonTags ) {
        
        for( std::string njetTag : njetTags ) {
            my_histos.emplace( "h_mva_"+njetTag+"j_"+leptonTag, std::make_shared<TH1D>( ( "h_mva_"+njetTag+"j_"+leptonTag ).c_str() , ( "h_mva_"+njetTag+"j_"+leptonTag ).c_str(), 10, 0.0, 1.0 ) );
        }

        for( std::string mvaTag : mvaTags ) {
            my_histos.emplace( "h_njets_"+mvaTag+"_"+leptonTag, std::make_shared<TH1D>( ( "h_njets_"+mvaTag+"_"+leptonTag ).c_str(), ( "h_njets_"+mvaTag+"_"+leptonTag ).c_str(), 6, 7, 13 ) );
        }

        for( std::string nBinTag : nBinTags ) {
            my_histos.emplace( "h_njets_"+nBinTag+"_"+leptonTag, std::make_shared<TH1D>( ( "h_njets_"+nBinTag+"_"+leptonTag ).c_str(), ( "h_njets_"+nBinTag+"_"+leptonTag ).c_str(), 6, 7, 13 ) );
        }
    }
    
    for( std::string nomWeight : nomWeights ) {
        for( std::string pdfWeight : pdfWeights ) {
            for( std::string sclWeight : sclWeights ) {
                for( std::string isrWeight : isrWeights ) {
                    for( std::string fsrWeight : fsrWeights ) {
                        
                        my_histos.emplace("h_jetpt_1_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_1_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_1_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_2_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_2_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_2_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_3_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_3_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_3_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_4_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_4_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_4_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_5_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_5_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_5_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_6_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_6_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_6_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_jetpt_7_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_jetpt_7_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_jetpt_7_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                    
                        my_histos.emplace("h_njets_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_njets_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_njets_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 6, 7, 13) );
                        my_histos.emplace("h_ht_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_ht_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_ht_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0.0, 3000.0 ) );
                        my_histos.emplace("h_nb_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_nb_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_nb_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 7, -0.5, 6 ) );
                        my_histos.emplace("h_mva_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_mva_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_mva_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 15, 0.0, 1.0 ) );
                        my_histos.emplace("h_nvtx_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_nvtx_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_nvtx_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0, 60 ) );
                        my_histos.emplace("h_mbl_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_mbl_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_mbl_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 20, 200 ) );
                        
                        my_histos.emplace("h_n_jetpt_1_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_jetpt_1_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_jetpt_1_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_n_jetpt_2_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_jetpt_2_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_jetpt_2_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_n_jetpt_3_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_jetpt_3_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_jetpt_3_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_n_jetpt_4_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_jetpt_4_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_jetpt_4_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_n_jetpt_5_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_jetpt_5_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_jetpt_5_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_n_jetpt_6_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_jetpt_6_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_jetpt_6_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                        my_histos.emplace("h_n_jetpt_7_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_jetpt_7_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_jetpt_7_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 50, 0.0, 1000.0 ) );
                    
                        my_histos.emplace("h_n_njets_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_njets_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_njets_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 6, 7, 13) );
                        my_histos.emplace("h_n_ht_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_ht_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_ht_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0.0, 3000.0 ) );
                        my_histos.emplace("h_n_mva_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_mva_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_mva_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 15, 0.0, 1.0 ) );
                        my_histos.emplace("h_n_nvtx_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight, std::make_shared<TH1D>( ( "h_n_nvtx_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), ( "h_n_nvtx_"+nomWeight+"_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight ).c_str(), 60, 0, 60 ) );
                    }
                }
            }
        }
    }
}

//Put everything you want to do per event here.
void MakeNewQCDSystematic::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );
        if( !isQuiet ) std::cout<<weight<<std::endl;

        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        const auto& runtype                 = tr.getVar<std::string>("runtype");     
        const auto& passMadHT               = tr.getVar<bool>("passMadHT");
        const auto& lumi                    = tr.getVar<double>("Lumi");

        const auto& NGoodJets_pt30          = tr.getVar<int>("NGoodJets_pt30");
        const auto& NNonIsoMuonJets_pt30    = tr.getVar<int>("NNonIsoMuonJets_pt30");
        
        const auto& passBaseline            = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaselineNIM         = tr.getVar<bool>("passBaseline1l_NonIsoMuon");
        const auto& NGoodMuons              = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons          = tr.getVar<int>("NGoodElectrons");
        
        const auto& deepESM_val             = tr.getVar<double>("deepESM_val");
        const auto& deepESM_valNonIsoMuon   = tr.getVar<double>("deepESM_valNonIsoMuon");

        const auto& deepESM_bin1            = tr.getVar<bool>("deepESM_bin1");
        const auto& deepESM_bin2            = tr.getVar<bool>("deepESM_bin2");
        const auto& deepESM_bin3            = tr.getVar<bool>("deepESM_bin3");
        const auto& deepESM_bin4            = tr.getVar<bool>("deepESM_bin4");
        
        const auto& deepESM_binNonIsoMuon1  = tr.getVar<bool>("deepESM_binNonIsoMuon1");
        const auto& deepESM_binNonIsoMuon2  = tr.getVar<bool>("deepESM_binNonIsoMuon2");
        const auto& deepESM_binNonIsoMuon3  = tr.getVar<bool>("deepESM_binNonIsoMuon3");
        const auto& deepESM_binNonIsoMuon4  = tr.getVar<bool>("deepESM_binNonIsoMuon4");

        const auto& NGoodBJets_pt30         = tr.getVar<int>("NGoodBJets_pt30");
        const auto& Mbl                     = tr.getVar<double>("Mbl");
        const auto& HT_trigger_pt30         = tr.getVar<double>("HT_trigger_pt30");        
        const auto& HT_NonIsoMuon_pt30      = tr.getVar<double>("HT_NonIsoMuon_pt30");        
        const auto& nVtx                    = tr.getVar<int>("NVtx");
        
        const auto& Jet_pt_1                = tr.getVar<double>("Jet_pt_1");
        const auto& Jet_pt_2                = tr.getVar<double>("Jet_pt_2");
        const auto& Jet_pt_3                = tr.getVar<double>("Jet_pt_3");
        const auto& Jet_pt_4                = tr.getVar<double>("Jet_pt_4");
        const auto& Jet_pt_5                = tr.getVar<double>("Jet_pt_5");
        const auto& Jet_pt_6                = tr.getVar<double>("Jet_pt_6");
        const auto& Jet_pt_7                = tr.getVar<double>("Jet_pt_7");
      
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------
        // -- Define weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
        // ------------------------
         double eventweight                 = 1.0;
         double totalEventWeight            = 1.0;
         double totalEventWeightNIM         = 1.0;
         double Lumi                        = 1.0;

         double sclUpFactor                 = 1.0;
         double sclDownFactor               = 1.0;
         double pdfUpFactor                 = 1.0;
         double pdfDownFactor               = 1.0;
         double isrUpFactor                 = 1.0;
         double isrDownFactor               = 1.0;
         double fsrUpFactor                 = 1.0;
         double fsrDownFactor               = 1.0;

        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            const auto& totalEventWeight_temp        = tr.getVar<double>("totalEventWeight");
            const auto& totalEventWeightNIM_temp     = tr.getVar<double>("totalEventWeightNIM");
            totalEventWeight                         = totalEventWeight_temp;
            totalEventWeightNIM                      = totalEventWeightNIM_temp;
            
            // Define Lumi weight*extraWeightDict[pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]
            const auto& Weight  = tr.getVar<double>("Weight");
            eventweight = lumi*Weight;

            pdfUpFactor                     = tr.getVar<double>("PDFweightUp");
            pdfDownFactor                   = tr.getVar<double>("PDFweightDown");
            sclUpFactor                     = tr.getVar<double>("scaleWeightUp");
            sclDownFactor                   = tr.getVar<double>("scaleWeightDown");
            isrUpFactor                     = tr.getVar<double>("PSweight_ISRUp");
            isrDownFactor                   = tr.getVar<double>("PSweight_ISRDown");
            fsrUpFactor                     = tr.getVar<double>("PSweight_FSRUp");
            fsrDownFactor                   = tr.getVar<double>("PSweight_FSRDown");
            Lumi                            = lumi;
        }

        std::map<std::string, double> extraWeightDict                   = {};
        double extraPdfWeight                                           = 1.0;
        double extraSclWeight                                           = 1.0;
        double extraIsrWeight                                           = 1.0;
        double extraFsrWeight                                           = 1.0;

        for( std::string pdfWeight : pdfWeights ) {
            
            extraPdfWeight                                              = 1.0;
            if( pdfWeight == "pdfUp" ) extraPdfWeight                   = pdfUpFactor;
            if( pdfWeight == "pdfDown" ) extraPdfWeight                 = pdfDownFactor;

            for( std::string sclWeight : sclWeights ) {
            
                extraSclWeight                                          = 1.0;
                if( sclWeight == "sclUp" ) extraSclWeight               = sclUpFactor;
                if( sclWeight == "sclDown" ) extraSclWeight             = sclDownFactor;

                for( std::string isrWeight : isrWeights ) {
                    
                    extraIsrWeight                                      = 1.0;
                    if( isrWeight == "isrUp" ) extraIsrWeight           = isrUpFactor;
                    if( isrWeight == "isrDown" ) extraIsrWeight         = isrDownFactor;

                    for( std::string fsrWeight : fsrWeights ) {
                        extraFsrWeight                                  = 1.0;
                        if( fsrWeight == "fsrUp" ) extraFsrWeight       = fsrUpFactor;
                        if( fsrWeight == "fsrDown" ) extraFsrWeight     = fsrDownFactor;

                        extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] = extraPdfWeight*extraSclWeight*extraIsrWeight*extraFsrWeight;
                    }
                }
            }
        }

        if( passBaselineNIM ) {
            
            for( std::string pdfWeight : pdfWeights ) {
                for( std::string sclWeight : sclWeights ) {
                    for( std::string isrWeight : isrWeights ) {
                        for( std::string fsrWeight : fsrWeights ) {
                            my_histos[ "h_n_njets_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NNonIsoMuonJets_pt30, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_ht_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_NonIsoMuon_pt30, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_nvtx_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_mva_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            
                            my_histos[ "h_n_jetpt_1_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_2_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_3_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_4_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_5_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_6_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_7_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, Lumi*totalEventWeightNIM*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            
                            my_histos[ "h_n_njets_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NNonIsoMuonJets_pt30, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_ht_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_NonIsoMuon_pt30, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_nvtx_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_mva_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            
                            my_histos[ "h_n_jetpt_1_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_2_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_3_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_4_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_5_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_6_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_n_jetpt_7_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                        }
                    }
                }
            }

            for( std::string njetTag : njetTags ) {
                if( std::stoi( njetTag ) == NNonIsoMuonJets_pt30 ) {
                    my_histos[ "h_mva_"+njetTag+"j_nonIsoMuon" ]->Fill( deepESM_valNonIsoMuon, Lumi*totalEventWeightNIM );
                }
                else if( NNonIsoMuonJets_pt30 > 11 ) {
                    my_histos[ "h_mva_11j_nonIsoMuon" ]->Fill( deepESM_valNonIsoMuon, Lumi*totalEventWeightNIM );
                }
            }
           
            if( NNonIsoMuonJets_pt30 < 12 ) {
                if( deepESM_binNonIsoMuon1 ) my_histos[ "h_njets_D1_nonIsoMuon" ]->Fill( NNonIsoMuonJets_pt30, Lumi*totalEventWeightNIM );
                if( deepESM_binNonIsoMuon2 ) my_histos[ "h_njets_D2_nonIsoMuon" ]->Fill( NNonIsoMuonJets_pt30, Lumi*totalEventWeightNIM );
                if( deepESM_binNonIsoMuon3 ) my_histos[ "h_njets_D3_nonIsoMuon" ]->Fill( NNonIsoMuonJets_pt30, Lumi*totalEventWeightNIM );
                if( deepESM_binNonIsoMuon4 ) my_histos[ "h_njets_D4_nonIsoMuon" ]->Fill( NNonIsoMuonJets_pt30, Lumi*totalEventWeightNIM );
                
                if( deepESM_valNonIsoMuon < 0.25 ) my_histos[ "h_njets_N1_nonIsoMuon" ]->Fill( NNonIsoMuonJets_pt30, Lumi*totalEventWeightNIM );
                if( deepESM_valNonIsoMuon >= 0.25 && deepESM_valNonIsoMuon < 0.50 ) my_histos[ "h_njets_N2_nonIsoMuon" ]->Fill( NNonIsoMuonJets_pt30, Lumi*totalEventWeightNIM );
                if( deepESM_valNonIsoMuon >= 0.50 && deepESM_valNonIsoMuon < 0.75 ) my_histos[ "h_njets_N3_nonIsoMuon" ]->Fill( NNonIsoMuonJets_pt30, Lumi*totalEventWeightNIM );
                if( deepESM_valNonIsoMuon >= 0.75 && deepESM_valNonIsoMuon < 1.00 ) my_histos[ "h_njets_N4_nonIsoMuon" ]->Fill( NNonIsoMuonJets_pt30, Lumi*totalEventWeightNIM );
            }
            else {
                if( deepESM_binNonIsoMuon1 ) my_histos[ "h_njets_D1_nonIsoMuon" ]->Fill( 12, Lumi*totalEventWeightNIM );
                if( deepESM_binNonIsoMuon2 ) my_histos[ "h_njets_D2_nonIsoMuon" ]->Fill( 12, Lumi*totalEventWeightNIM );
                if( deepESM_binNonIsoMuon3 ) my_histos[ "h_njets_D3_nonIsoMuon" ]->Fill( 12, Lumi*totalEventWeightNIM );
                if( deepESM_binNonIsoMuon4 ) my_histos[ "h_njets_D4_nonIsoMuon" ]->Fill( 12, Lumi*totalEventWeightNIM );
                
                if( deepESM_valNonIsoMuon < 0.25 ) my_histos[ "h_njets_N1_nonIsoMuon" ]->Fill( 12, Lumi*totalEventWeightNIM );
                if( deepESM_valNonIsoMuon >= 0.25 && deepESM_valNonIsoMuon < 0.50 ) my_histos[ "h_njets_N2_nonIsoMuon" ]->Fill( 12, Lumi*totalEventWeightNIM );
                if( deepESM_valNonIsoMuon >= 0.50 && deepESM_valNonIsoMuon < 0.75 ) my_histos[ "h_njets_N3_nonIsoMuon" ]->Fill( 12, Lumi*totalEventWeightNIM );
                if( deepESM_valNonIsoMuon >= 0.75 && deepESM_valNonIsoMuon < 1.00 ) my_histos[ "h_njets_N4_nonIsoMuon" ]->Fill( 12, Lumi*totalEventWeightNIM );
            }
        }

        if( passBaseline ) {

            for( std::string pdfWeight : pdfWeights ) {
                for( std::string sclWeight : sclWeights ) {
                    for( std::string isrWeight : isrWeights ) {
                        for( std::string fsrWeight : fsrWeights ) {
                            my_histos[ "h_njets_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodJets_pt30, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_ht_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_trigger_pt30, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_nb_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodBJets_pt30, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_nvtx_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_mbl_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Mbl, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_mva_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            
                            my_histos[ "h_jetpt_1_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_2_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_3_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_4_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_5_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_6_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_7_allWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, Lumi*totalEventWeight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            
                            my_histos[ "h_njets_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodJets_pt30, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_ht_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( HT_trigger_pt30, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_nb_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( NGoodBJets_pt30, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_nvtx_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( nVtx, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_mbl_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Mbl, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_mva_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( deepESM_val, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            
                            my_histos[ "h_jetpt_1_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_1, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_2_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_2, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_3_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_3, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_4_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_4, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_5_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_5, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_6_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_6, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                            my_histos[ "h_jetpt_7_noWghts_"+pdfWeight+"_"+sclWeight+"_"+isrWeight+"_"+fsrWeight]->Fill( Jet_pt_7, Lumi*eventweight*extraWeightDict[pdfWeight+sclWeight+isrWeight+fsrWeight] );
                        }
                    }
                }
            }

            for( std::string njetTag : njetTags ) {
                if( NGoodMuons == 1 && std::stoi( njetTag ) == NGoodJets_pt30 ) {
                    my_histos[ "h_mva_"+njetTag+"j_isoMuon" ]->Fill( deepESM_val, Lumi*totalEventWeight );
                }
                else if( NGoodMuons == 1 && NGoodJets_pt30 > 11 ) {
                    my_histos[ "h_mva_11j_isoMuon" ]->Fill( deepESM_val, Lumi*totalEventWeight );
                }
                else if( NGoodElectrons == 1 && std::stoi( njetTag ) == NGoodJets_pt30 ) {
                    my_histos[ "h_mva_"+njetTag+"j_isoElectron" ]->Fill( deepESM_val, Lumi*totalEventWeight );
                }
                else if( NGoodElectrons == 1 && NGoodJets_pt30 > 11) {
                    my_histos[ "h_mva_11j_isoElectron" ]->Fill( deepESM_val, Lumi*totalEventWeight );
                }
            }

            if( NGoodMuons == 1 ){
                if( NGoodJets_pt30 < 12 ) {
                    if( deepESM_bin1 ) my_histos[ "h_njets_D1_isoMuon" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_bin2 ) my_histos[ "h_njets_D2_isoMuon" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_bin3 ) my_histos[ "h_njets_D3_isoMuon" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_bin4 ) my_histos[ "h_njets_D4_isoMuon" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                
                    if( deepESM_val < 0.25 ) my_histos[ "h_njets_N1_isoMuon" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.25 && deepESM_val < 0.50 ) my_histos[ "h_njets_N2_isoMuon" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.50 && deepESM_val < 0.75 ) my_histos[ "h_njets_N3_isoMuon" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.75 && deepESM_val < 1.00 ) my_histos[ "h_njets_N4_isoMuon" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                }
                else {
                    if( deepESM_bin1 ) my_histos[ "h_njets_D1_isoMuon" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_bin2 ) my_histos[ "h_njets_D2_isoMuon" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_bin3 ) my_histos[ "h_njets_D3_isoMuon" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_bin4 ) my_histos[ "h_njets_D4_isoMuon" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_val < 0.25 ) my_histos[ "h_njets_N1_isoMuon" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.25 && deepESM_val < 0.50 ) my_histos[ "h_njets_N2_isoMuon" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.50 && deepESM_val < 0.75 ) my_histos[ "h_njets_N3_isoMuon" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.75 && deepESM_val < 1.00 ) my_histos[ "h_njets_N4_isoMuon" ]->Fill( 12, Lumi*totalEventWeight );
                }
            }
            
            if( NGoodElectrons == 1 ){
                if( NGoodJets_pt30 < 12 ) {
                    if( deepESM_bin1 ) my_histos[ "h_njets_D1_isoElectron" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_bin2 ) my_histos[ "h_njets_D2_isoElectron" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_bin3 ) my_histos[ "h_njets_D3_isoElectron" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_bin4 ) my_histos[ "h_njets_D4_isoElectron" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    
                    if( deepESM_val < 0.25 ) my_histos[ "h_njets_N1_isoElectron" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.25 && deepESM_val < 0.50 ) my_histos[ "h_njets_N2_isoElectron" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.50 && deepESM_val < 0.75 ) my_histos[ "h_njets_N3_isoElectron" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.75 && deepESM_val < 1.00 ) my_histos[ "h_njets_N4_isoElectron" ]->Fill( NGoodJets_pt30, Lumi*totalEventWeight );
                }
                else {
                    if( deepESM_bin1 ) my_histos[ "h_njets_D1_isoElectron" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_bin2 ) my_histos[ "h_njets_D2_isoElectron" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_bin3 ) my_histos[ "h_njets_D3_isoElectron" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_bin4 ) my_histos[ "h_njets_D4_isoElectron" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_val < 0.25 ) my_histos[ "h_njets_N1_isoElectron" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.25 && deepESM_val < 0.50 ) my_histos[ "h_njets_N2_isoElectron" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.50 && deepESM_val < 0.75 ) my_histos[ "h_njets_N3_isoElectron" ]->Fill( 12, Lumi*totalEventWeight );
                    if( deepESM_val >= 0.75 && deepESM_val < 1.00 ) my_histos[ "h_njets_N4_isoElectron" ]->Fill( 12, Lumi*totalEventWeight );
                }
            }
        }
    } 
}

void MakeNewQCDSystematic::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }

    for (const auto &p : my_3d_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }
    
}
