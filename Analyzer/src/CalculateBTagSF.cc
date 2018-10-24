#define CalculateBTagSF_cxx
#include "Analyzer/Analyzer/include/CalculateBTagSF.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

CalculateBTagSF::CalculateBTagSF() : initHistos(false)
{
}


void CalculateBTagSF::InitHistos(std::string histoFileTag)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    const int nPtBins   = 17;
    const int nEtaBins  = 3;

    const double ptBins[]  = { 20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1000 };
    const double etaBins[] = { 0.0, 0.8, 1.6, 2.4 };

    // Declare all your histograms here, that way we can fill them for multiple chains

    std::cout<<histoFileTag<<std::endl;

    my_2d_histos.emplace( "n_eff_b_"+histoFileTag,     std::make_shared<TH2D>( ("n_eff_b_"+histoFileTag ).c_str(),     ( "n_eff_b_Efficiency_"+histoFileTag ).c_str(),     nPtBins, ptBins, nEtaBins, etaBins ) );
    my_2d_histos.emplace( "n_eff_c_"+histoFileTag,     std::make_shared<TH2D>( ("n_eff_c_"+histoFileTag ).c_str(),     ( "n_eff_c_Efficiency_"+histoFileTag ).c_str(),     nPtBins, ptBins, nEtaBins, etaBins ) );
    my_2d_histos.emplace( "n_eff_udsg_"+histoFileTag,  std::make_shared<TH2D>( ("n_eff_udsg_"+histoFileTag ).c_str(),  ( "n_eff_udsg_Efficiency_"+histoFileTag ).c_str(),  nPtBins, ptBins, nEtaBins, etaBins ) );
    my_2d_histos.emplace( "d_eff_b_"+histoFileTag,     std::make_shared<TH2D>( ("d_eff_b_"+histoFileTag ).c_str(),     ( "d_eff_b_Efficiency_"+histoFileTag ).c_str(),     nPtBins, ptBins, nEtaBins, etaBins ) );
    my_2d_histos.emplace( "d_eff_c_"+histoFileTag,     std::make_shared<TH2D>( ("d_eff_c_"+histoFileTag ).c_str(),     ( "d_eff_c_Efficiency_"+histoFileTag ).c_str(),     nPtBins, ptBins, nEtaBins, etaBins ) );
    my_2d_histos.emplace( "d_eff_udsg_"+histoFileTag,  std::make_shared<TH2D>( ("d_eff_udsg_"+histoFileTag ).c_str(),  ( "d_eff_udsg_Efficiency_"+histoFileTag ).c_str(),  nPtBins, ptBins, nEtaBins, etaBins ) );

    my_2d_histos["n_eff_b_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_eff_b_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );
    
    my_2d_histos["n_eff_c_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_eff_c_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );
    
    my_2d_histos["n_eff_udsg_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_eff_udsg_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );

/*    std::vector<std::string> genericVectorOfTags    { "myTag" };

    for( std::string nJetsTag : njets_fisherDist_tags ) {
        my_histos.emplace( "h_fishDist_"+nJetsTag+"j_"+jetPtTag+"_"+fishBinTag, std::make_shared<TH1D>( ( "h_fishDist_"+nJetsTag+"j_"+jetPtTag+"_"+fishBinTag ).c_str(), ( "h_fishDist_"+nJetsTag+"j_"+jetPtTag+"_"+fishBinTag ).c_str(), 120, -0.6, 0.6 ) );
    }//END of njets_fisherDist_tags*/

}//END of init histos

void CalculateBTagSF::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while( tr.getNextEvent() )
    {
        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& TriggerNames        = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass         = tr.getVec<int>("TriggerPass");

        const auto& NJets_pt30          = tr.getVar<int>("NGoodJets_pt30");
        const auto& NJets_pt45          = tr.getVar<int>("NGoodJets_pt45");

        const auto& HT_trigger          = tr.getVar<double>("HT_trigger");

        const auto& passBlindHad        = tr.getVar<bool>("passBlindHad");
        const auto& passBlindLep        = tr.getVar<bool>("passBlindLep");
        const auto& passTrigger         = tr.getVar<bool>("passTrigger");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        
        const auto& Jets                = tr.getVec<TLorentzVector>("Jets");
        const auto& recoJetsBtag        = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const auto& recoJetsFlavor      = tr.getVec<int>("Jets_hadronFlavor");

        //-----------------------------------
        //-- Initialize Histograms
        //----------------------------------
        
        if( !initHistos ) {
            InitHistos(filetag);
            initHistos = true;
        }
        //------------------------------------
        //-- Print Event Number
        //------------------------------------

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        //------------------------------------
        //-- Get the Proper Event Weight
        //------------------------------------

        double eventweight          = 1.0; //For data
        double PUaddweight          = 1.0, PU_JPaddweight           = 1.0;
        double scaleSFaddweight     = 1.0, scaleSFaddweight_u       = 1.0, scaleSFaddweight_d       = 1.0;
        double pdfSFaddweight       = 1.0, pdfSFaddweight_u         = 1.0, pdfSFaddweight_d         = 1.0;
        double Lumi = 35900;

        if( passMadHT )
        {
            const auto& Weight = tr.getVar<double>("Weight");
            eventweight         = Lumi*Weight;
            
            for( unsigned int ij = 0; ij < Jets.size(); ++ij ) {
                
                if( !JetPassCuts( Jets.at( ij ) ) ) continue;
                int myJetFlavor = std::abs( recoJetsFlavor.at( ij ) );
                double myJetCSV = recoJetsBtag.at( ij );
    
                if( myJetFlavor == 5 ) {
                    my_2d_histos["d_eff_b_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight );
                    if( myJetCSV > 0.8484 ) { my_2d_histos["n_eff_b_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight ); }
                }
                else if( myJetFlavor == 4 ) {
                    my_2d_histos["d_eff_c_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight );
                    if( myJetCSV > 0.8484 ) { my_2d_histos["n_eff_c_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight ); }
                }
                else if( myJetFlavor < 4 || myJetFlavor == 21 ) {
                    my_2d_histos["d_eff_udsg_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight );
                    if( myJetCSV > 0.8484 ) { my_2d_histos["n_eff_udsg_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight ); }
                }
            }
        }
        else continue;
    }//END of while tr.getNextEvent loop   
}//END of function
      
void CalculateBTagSF::WriteHistos( TFile* outfile ) 
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

bool CalculateBTagSF::JetPassCuts( const TLorentzVector& jet )
{
    const double minAbsEta = 0.0;
    const double maxAbsEta = 2.4;
    const double minPt     = 30.0;
    const double maxPt     = 999.9;

    return ( minAbsEta == -1 || std::fabs( jet.Eta() ) >= minAbsEta )
        && ( maxAbsEta == -1 || std::fabs( jet.Eta() ) <  maxAbsEta )
        && ( minPt == -1     || jet.Pt() >= minPt )
        && ( maxPt == -1     || jet.Pt() <  maxPt );

}
