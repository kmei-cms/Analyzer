#define AnalyzeWControlRegion_cxx
#include "Analyzer/Analyzer/include/AnalyzeWControlRegion.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <iostream>

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "Framework/Framework/include/SetUpTopTagger.h"

AnalyzeWControlRegion::AnalyzeWControlRegion() : initHistos(false)
{
}

void AnalyzeWControlRegion::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains

    for(auto& mycut : cutMap)
    {
        my_histos.emplace("h_njets_"+mycut.first,std::make_shared<TH1D>(("h_njets_"+mycut.first).c_str(),("h_njets_"+mycut.first).c_str(),15,0,15));
        my_histos.emplace("h_nbjets_"+mycut.first,std::make_shared<TH1D>(("h_nbjets_"+mycut.first).c_str(),("h_nbjets_"+mycut.first).c_str(),5,0,5));
        my_histos.emplace("h_met_"+mycut.first,std::make_shared<TH1D>(("h_met_"+mycut.first).c_str(),("h_met_"+mycut.first).c_str(),50,0,500));
        my_histos.emplace("h_mT_"+mycut.first,std::make_shared<TH1D>(("h_mT_"+mycut.first).c_str(),("h_mT_"+mycut.first).c_str(),40,0,200));
        my_histos.emplace("h_HT_"+mycut.first,std::make_shared<TH1D>(("h_HT_"+mycut.first).c_str(),("h_HT_"+mycut.first).c_str(),30,0,1500));
        my_histos.emplace("h_dphi_"+mycut.first,std::make_shared<TH1D>(("h_dphi_"+mycut.first).c_str(),("h_dphi_"+mycut.first).c_str(),50,0,5));
        my_histos.emplace("h_Mbl_maxpT_"+mycut.first,std::make_shared<TH1D>(("h_Mbl_maxpT_"+mycut.first).c_str(),("h_Mbl_maxpT_"+mycut.first).c_str(),30,0,300));
        my_histos.emplace("h_Mbl_maxCSV_"+mycut.first,std::make_shared<TH1D>(("h_Mbl_maxCSV_"+mycut.first).c_str(),("h_Mbl_maxCSV_"+mycut.first).c_str(),30,0,300));
        my_histos.emplace("h_Mbl_all_"+mycut.first,std::make_shared<TH1D>(("h_Mbl_all_"+mycut.first).c_str(),("h_Mbl_all_"+mycut.first).c_str(),30,0,300));

        my_histos.emplace("h_lepton_pT_"+mycut.first,std::make_shared<TH1D>(("h_lepton_pT_"+mycut.first).c_str(),("h_lepton_pT_"+mycut.first).c_str(),50,0,500));
        my_histos.emplace("h_muon_pT_"+mycut.first,std::make_shared<TH1D>(("h_muon_pT_"+mycut.first).c_str(),("h_muon_pT_"+mycut.first).c_str(),50,0,500));
        my_histos.emplace("h_electron_pT_"+mycut.first,std::make_shared<TH1D>(("h_electron_pT_"+mycut.first).c_str(),("h_electron_pT_"+mycut.first).c_str(),50,0,500));
        my_histos.emplace("h_lepton_eta_"+mycut.first,std::make_shared<TH1D>(("h_lepton_eta_"+mycut.first).c_str(),("h_lepton_eta_"+mycut.first).c_str(),50,-2.5,2.5));
        my_histos.emplace("h_muon_eta_"+mycut.first,std::make_shared<TH1D>(("h_muon_eta_"+mycut.first).c_str(),("h_muon_eta_"+mycut.first).c_str(),50,-2.5,2.5));
        my_histos.emplace("h_electron_eta_"+mycut.first,std::make_shared<TH1D>(("h_electron_eta_"+mycut.first).c_str(),("h_electron_eta_"+mycut.first).c_str(),50,-2.5,2.5));
        my_histos.emplace("h_lepton_phi_"+mycut.first,std::make_shared<TH1D>(("h_lepton_phi_"+mycut.first).c_str(),("h_lepton_phi_"+mycut.first).c_str(),40,-4,4));
        my_histos.emplace("h_muon_phi_"+mycut.first,std::make_shared<TH1D>(("h_muon_phi_"+mycut.first).c_str(),("h_muon_phi_"+mycut.first).c_str(),40,-4,4));
        my_histos.emplace("h_electron_phi_"+mycut.first,std::make_shared<TH1D>(("h_electron_phi_"+mycut.first).c_str(),("h_electron_phi_"+mycut.first).c_str(),40,-4,4));
    }   
}

void AnalyzeWControlRegion::Loop(NTupleReader& tr, double weight, int maxevents, bool isQuiet)
{
    while(tr.getNextEvent())
    {
        //const auto& ntops                   = tr.getVar<int>("ntops");
        //const auto& ntops_3jet              = tr.getVar<int>("ntops_3jet");
        //const auto& ntops_2jet              = tr.getVar<int>("ntops_2jet");
        //const auto& ntops_1jet              = tr.getVar<int>("ntops_1jet");
        const auto& runtype                 = tr.getVar<std::string>("runtype");
        const auto& filetag                 = tr.getVar<std::string>("filetag");
        const auto& JetID                   = tr.getVar<bool>("JetID");
        const auto& passTrigger             = tr.getVar<bool>("passTrigger");

        const auto& NGoodJets_pt30          = tr.getVar<int>("NGoodJets_pt30");
        const auto& GoodJets_pt30           = tr.getVec<TLorentzVector>("GoodJets_pt30");
        const auto& GoodJets_pt30_CSV       = tr.getVec<double>("GoodJets_pt30_bDiscriminatorCSV");
        //const auto& NJets_pt45          = tr.getVar<int>("NJets_pt45");
        //const auto& BJets_pt30          = tr.getVec<TLorentzVector>("BJets_pt30");
        //const auto& NBJets_pt30         = tr.getVar<int>("NBJets_pt30");
        //const auto& BJets_pt45          = tr.getVec<TLorentzVector>("BJets_pt45");
        //const auto& NBJets_pt45         = tr.getVar<int>("NBJets_pt45");
        const auto& NGoodBJets_pt30     = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodBJets_pt30_loose     = tr.getVar<int>("NGoodBJets_pt30_loose");
        const auto& GoodLeptons         = tr.getVec<TLorentzVector>("GoodLeptons");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& GoodMuons           = tr.getVec<TLorentzVector>("GoodMuons");
        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");
        const auto& GoodMuonsMTW        = tr.getVec<double>("GoodMuonsMTW");
        const auto& GoodElectrons       = tr.getVec<TLorentzVector>("GoodElectrons");
        const auto& NGoodElectrons      = tr.getVar<int>("NGoodElectrons");
        const auto& GoodElectronsMTW    = tr.getVec<double>("GoodElectronsMTW");
        //const auto& HT_trigger          = tr.getVar<double>("HT_trigger");
        const auto& Mbl                 = tr.getVar<double>("Mbl");
        const auto& HT                  = tr.getVar<double>("HT");
        const auto& MET                 = tr.getVar<double>("MET");
        const auto& METPhi              = tr.getVar<double>("METPhi");

        //const auto& passBaseline0l    = tr.getVar<bool>("passBaseline0l");
        //const auto& passBaseline1l    = tr.getVar<bool>("passBaseline1l");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");

        if (maxevents > 0 && tr.getEvtNum() >= maxevents) break;        

        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;
        
        // Make sure event weight is not 0 for data
        double eventweight = 1.;
        if(runtype != "Data")
        {
            eventweight = tr.getVar<double>("Weight");
        }        

        // -------------------------------
        // -- Basic event selection stuff
        // -------------------------------

        // basic event cleaning
        bool goodEvent = JetID && passMadHT && passTrigger;
        // Exclude events with bad jets
        if(!goodEvent) continue;

        // lepton cuts
        bool pass_1mu = (runtype != "Data" || filetag == "Data_SingleMuon")     && NGoodMuons == 1;
        bool pass_1el = (runtype != "Data" || filetag == "Data_SingleElectron") && NGoodElectrons == 1;
        bool pass_1l = pass_1mu || pass_1el;
        // compute mT 
        TLorentzVector metLV;
        metLV.SetPtEtaPhiE(MET,0,METPhi,MET);
        double mT = -1;
        //std::cout << "MET , METPhi: " << MET << ", " << METPhi << std::endl;
        //std::cout << "Nleptons: " << NGoodLeptons << std::endl;
        if (NGoodLeptons > 0)
        {
            //std::cout << "Found lepton, pt " << GoodLeptons[0].Pt() << std::endl; 
            mT = utility::calcMT(GoodLeptons[0], metLV);
        }
        bool pass_mT = mT > 30 && mT < 100;

        // jet cuts
        bool pass_0b = NGoodBJets_pt30 == 0;
        // compute dphi with highest pT jet
        double dphi = -1;
        if (NGoodJets_pt30 > 0)
        {
            dphi = utility::calcDPhi(METPhi,GoodJets_pt30[0].Phi());
        }
        std::vector<double> Mbl_all;
        double Mbl_maxCSV = -1;
        double Mbl_maxpT = -1;
        if (pass_1l && NGoodJets_pt30 > 0)
        {
            Mbl_maxpT = (GoodLeptons[0]+GoodJets_pt30[0]).M();
            int maxCSV_index = 0;
            for (int i=0; i<NGoodJets_pt30; ++i)
            {
                if(GoodJets_pt30_CSV[i] > GoodJets_pt30_CSV[maxCSV_index])
                    maxCSV_index = i;
                double temp = (GoodLeptons[0]+GoodJets_pt30[i]).M();
                Mbl_all.push_back(temp);
            }
            Mbl_maxCSV = (GoodLeptons[0]+GoodJets_pt30[maxCSV_index]).M();
        }

        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map 
        {
            {"1l"                                 , goodEvent && pass_1l                                                             },
            {"1l_mT"                              , goodEvent && pass_1l && pass_mT                                                  },
            {"1l_mT_0b"                           , goodEvent && pass_1l && pass_mT && pass_0b                                       },
        };
        
        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map);
            initHistos = true;
        }

        // Fill Histograms
        for(auto& kv : cut_map)
        {
            if(kv.second)
            {
                my_histos["h_njets_"   +kv.first]->Fill(NGoodJets_pt30, eventweight);
                my_histos["h_nbjets_"  +kv.first]->Fill(NGoodBJets_pt30, eventweight);
                my_histos["h_met_"     +kv.first]->Fill(MET, eventweight);
                my_histos["h_mT_"      +kv.first]->Fill(mT, eventweight);
                my_histos["h_HT_"      +kv.first]->Fill(HT, eventweight);

                my_histos["h_lepton_pT_"   +kv.first]->Fill(GoodLeptons[0].Pt(), eventweight);
                my_histos["h_lepton_eta_"  +kv.first]->Fill(GoodLeptons[0].Eta(), eventweight);
                my_histos["h_lepton_phi_"  +kv.first]->Fill(GoodLeptons[0].Phi(), eventweight);

                if (NGoodMuons > 0)
                {
                    my_histos["h_muon_pT_"   +kv.first]->Fill(GoodMuons[0].Pt(), eventweight);
                    my_histos["h_muon_eta_"  +kv.first]->Fill(GoodMuons[0].Eta(), eventweight);
                    my_histos["h_muon_phi_"  +kv.first]->Fill(GoodMuons[0].Phi(), eventweight);
                } 
                if (NGoodElectrons > 0)
                {
                    my_histos["h_electron_pT_"   +kv.first]->Fill(GoodElectrons[0].Pt(), eventweight);
                    my_histos["h_electron_eta_"  +kv.first]->Fill(GoodElectrons[0].Eta(), eventweight);
                    my_histos["h_electron_phi_"  +kv.first]->Fill(GoodElectrons[0].Phi(), eventweight);
                }
            }
        }
    }
}



void AnalyzeWControlRegion::WriteHistos(TFile* outfile)
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

