#include "Analyzer/Analyzer/include/StealthHemispheres.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/interface/TTException.h"
#include "Framework/Framework/include/SetUpTopTagger.h"


#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"

StealthHemispheres::StealthHemispheres() : inithisto(false)
{
}

void StealthHemispheres::InitHistos(const std::map<std::string, bool>& cutmap) 
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    for (const auto& cutVar : cutmap)
    {
        // 1 Lepton Case
        my_histos.emplace( "h_MT2_"+cutVar.first, std::make_shared<TH1D> ( ("h_MT2_"+cutVar.first).c_str(), ("h_MT2_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_"+cutVar.first).c_str(), ("h_stop1Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Eta_"+cutVar.first).c_str(), ("h_stop1Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop1Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Phi_"+cutVar.first).c_str(), ("h_stop1Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop1Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Pt_"+cutVar.first).c_str(), ("h_stop1Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_stop2Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_"+cutVar.first).c_str(), ("h_stop2Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop2Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Eta_"+cutVar.first).c_str(), ("h_stop2Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop2Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Phi_"+cutVar.first).c_str(), ("h_stop2Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop2Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Pt_"+cutVar.first).c_str(), ("h_stop2Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_dR_stop1stop2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_stop1stop2_"+cutVar.first).c_str(), ("h_dR_stop1stop2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dPhi_stop1stop2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dPhi_stop1stop2_"+cutVar.first).c_str(), ("h_dPhi_stop1stop2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_2d_histos.emplace( "h_Mass_stop1vsstop2_"+cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2_"+cutVar.first).c_str(), ("h_Mass_stop1vsstop2_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Eta_stop1vsstop2_"+cutVar.first, std::make_shared<TH2D>( ("h_Eta_stop1vsstop2_"+cutVar.first).c_str(), ("h_Eta_stop1vsstop2_"+cutVar.first).c_str(), 100, -6, 6, 100, -6, 6 ) );
        my_2d_histos.emplace( "h_Phi_stop1vsstop2_"+cutVar.first, std::make_shared<TH2D>( ("h_Phi_stop1vsstop2_"+cutVar.first).c_str(), ("h_Phi_stop1vsstop2_"+cutVar.first).c_str(), 80, -4, 4, 80, -4, 4 ) );
        my_2d_histos.emplace( "h_Pt_stop1vsstop2_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt_stop1vsstop2_"+cutVar.first).c_str(), ("h_Pt_stop1vsstop2_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000 ) );
        my_2d_histos.emplace( "h_MT2VsNJets_"+cutVar.first, std::make_shared<TH2D> ( ("h_MT2VsNJets_"+cutVar.first).c_str(), ("h_MT2VsNJets_"+cutVar.first).c_str(), 500, 0, 1500, 20, 0, 20 ) );
        my_2d_histos.emplace( "h_MT2Vsstop1Mass_"+cutVar.first, std::make_shared<TH2D> ( ("h_MT2Vsstop1Mass_"+cutVar.first).c_str(), ("h_MT2Vsstop1Mass_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_MT2Vsstop2Mass_"+cutVar.first, std::make_shared<TH2D> ( ("h_MT2Vsstop2Mass_"+cutVar.first).c_str(), ("h_MT2Vsstop2Mass_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_stop1MassVsNJets_"+cutVar.first, std::make_shared<TH2D> ( ("h_stop1MassVsNJets_"+cutVar.first).c_str(), ("h_stop1MassVsNJets_"+cutVar.first).c_str(), 500, 0, 1500, 20, 0, 20 ) );
        my_2d_histos.emplace( "h_stop2MassVsNJets_"+cutVar.first, std::make_shared<TH2D> ( ("h_stop2MassVsNJets_"+cutVar.first).c_str(), ("h_stop2MassVsNJets_"+cutVar.first).c_str(), 500, 0, 1500, 20, 0, 20 ) );

    }
}

void StealthHemispheres::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {

        const auto& eventCounter    = tr.getVar<int>("eventCounter");
        
        // ------------------------
        // -- Print event number   
        // ------------------------     
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& NGoodJets_pt45      = tr.getVar<int>("NGoodJets_pt45");
        const auto& passBaseline1l_Good = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline2l_Good = tr.getVar<bool>("passBaseline2l_Good");       
 
        // -------------------------------
        // -- MT2 hemispheres variables
        // -------------------------------
        const auto& MT2                = tr.getVar<double>("MT2_1l");
        const auto& stop1Mass          = tr.getVar<double>("stop1Mass_1l");
        const auto& stop1Eta           = tr.getVar<double>("stop1Eta_1l");
        const auto& stop1Phi           = tr.getVar<double>("stop1Phi_1l");
        const auto& stop1Pt            = tr.getVar<double>("stop1Pt_1l");
        const auto& stop2Mass          = tr.getVar<double>("stop2Mass_1l");
        const auto& stop2Eta           = tr.getVar<double>("stop2Eta_1l");
        const auto& stop2Phi           = tr.getVar<double>("stop2Phi_1l");
        const auto& stop2Pt            = tr.getVar<double>("stop2Pt_1l");
        const auto& dR_stop1stop2      = tr.getVar<double>("dR_stop1stop2_1l");
        const auto& dPhi_stop1stop2    = tr.getVar<double>("dPhi_stop1stop2_1l");

        // -------------------
        // -- Define weight
        // -------------------
        double weight=1.0, weightNoHT=1.0, weightQCDCR=1.0, weightNoBTag=1.0;
        double eventweight=1.0, leptonweight=1.0, bTagWeight=1.0, prefiringScaleFactor=1.0, pileupWeight=1.0, htDerivedweight=1.0;
        double weightNoLepton=1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight       = tr.getVar<double>("Weight");
            const auto& lumi         = tr.getVar<double>("Lumi");
            eventweight              = lumi*Weight;
            
            const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
            const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
            const auto& muNonIso     = tr.getVar<double>("totNonIsoMuonSF");
            leptonweight             = eleLepWeight*muLepWeight;
            
            pileupWeight             = tr.getVar<double>("puWeightCorr");
            bTagWeight               = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedweight          = tr.getVar<double>("htDerivedweight");
            prefiringScaleFactor     = tr.getVar<double>("prefiringScaleFactor");
            
            weightQCDCR    *= eventweight*muNonIso*prefiringScaleFactor*pileupWeight;
            weightNoHT     *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight;
            weightNoLepton *= eventweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight;
            weightNoBTag   *= eventweight*leptonweight*prefiringScaleFactor*pileupWeight*htDerivedweight;
            weight         *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight;
        }

        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            {"baseline_1l", passBaseline1l_Good},
            {"baseline_2l", passBaseline2l_Good},
        };

        if (!inithisto) 
        {
            InitHistos(cutmap);
            inithisto = true;
        }
       
        my_histos["EventCounter"]->Fill( eventCounter );
 
        for (const auto& cutVar: cutmap) 
        {
            if (cutVar.second) 
            {
                // 1 Lepton Case
                my_histos["h_MT2_"+cutVar.first]->Fill( MT2, weight );
                my_histos["h_stop1Mass_"+cutVar.first]->Fill( stop1Mass, weight );
                my_histos["h_stop1Eta_"+cutVar.first]->Fill( stop1Eta, weight );
                my_histos["h_stop1Phi_"+cutVar.first]->Fill( stop1Phi, weight );
                my_histos["h_stop1Pt_"+cutVar.first]->Fill( stop1Pt, weight );
                my_histos["h_stop2Mass_"+cutVar.first]->Fill( stop2Mass, weight );
                my_histos["h_stop2Eta_"+cutVar.first]->Fill( stop2Eta, weight );
                my_histos["h_stop2Phi_"+cutVar.first]->Fill( stop2Phi, weight );
                my_histos["h_stop2Pt_"+cutVar.first]->Fill( stop2Pt, weight );
                my_histos["h_dR_stop1stop2_"+cutVar.first]->Fill( dR_stop1stop2, weight );
                my_histos["h_dPhi_stop1stop2_"+cutVar.first]->Fill( dPhi_stop1stop2, weight );
                my_2d_histos["h_Mass_stop1vsstop2_"+cutVar.first]->Fill( stop1Mass, stop2Mass, weight );
                my_2d_histos["h_Mass_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{1} [GeV]");
                my_2d_histos["h_Mass_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{2} [GeV]");
                my_2d_histos["h_Eta_stop1vsstop2_"+cutVar.first]->Fill( stop1Eta, stop2Eta, weight );
                my_2d_histos["h_Eta_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("#eta_{#tildet}_{1}");
                my_2d_histos["h_Eta_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("#eta_{#tildet}_{2}");
                my_2d_histos["h_Phi_stop1vsstop2_"+cutVar.first]->Fill( stop1Phi, stop2Phi, weight );
                my_2d_histos["h_Phi_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("#phi_{#tildet}_{1}");
                my_2d_histos["h_Phi_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("#phi_{#tildet}_{2}");
                my_2d_histos["h_Pt_stop1vsstop2_"+cutVar.first]->Fill( stop1Pt, stop2Pt, weight );
                my_2d_histos["h_Pt_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("pT_{#tildet}_{1}");
                my_2d_histos["h_Pt_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("pT_{#tildet}_{2}");
                my_2d_histos["h_MT2VsNJets_"+cutVar.first]->Fill( MT2,  NGoodJets_pt45, weight );
                my_2d_histos["h_MT2VsNJets_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_MT2VsNJets_"+cutVar.first]->GetYaxis()->SetTitle("N_{J}");
                my_2d_histos["h_MT2Vsstop1Mass_"+cutVar.first]->Fill( MT2, stop1Mass, weight );
                my_2d_histos["h_MT2Vsstop1Mass_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_MT2Vsstop1Mass_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{1} [GeV]");
                my_2d_histos["h_MT2Vsstop2Mass_"+cutVar.first]->Fill( MT2, stop2Mass, weight );
                my_2d_histos["h_MT2Vsstop2Mass_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_MT2Vsstop2Mass_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{2} [GeV]");
                my_2d_histos["h_stop1MassVsNJets_"+cutVar.first]->Fill( stop1Mass, NGoodJets_pt45, weight );
                my_2d_histos["h_stop1MassVsNJets_"+cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{1}");
                my_2d_histos["h_stop1MassVsNJets_"+cutVar.first]->GetYaxis()->SetTitle("N_{J}");
                my_2d_histos["h_stop2MassVsNJets_"+cutVar.first]->Fill( stop2Mass, NGoodJets_pt45, weight );
                my_2d_histos["h_stop2MassVsNJets_"+cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{2}");
                my_2d_histos["h_stop2MassVsNJets_"+cutVar.first]->GetYaxis()->SetTitle("N_{J}");
         
            }
        }
    } 
} 

void StealthHemispheres::WriteHistos(TFile* outfile)
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
