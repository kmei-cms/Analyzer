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
        my_histos.emplace("h_nbjets_loose_"+mycut.first,std::make_shared<TH1D>(("h_nbjets_loose_"+mycut.first).c_str(),("h_nbjets_loose_"+mycut.first).c_str(),5,0,5));
        my_histos.emplace("h_met_"+mycut.first,std::make_shared<TH1D>(("h_met_"+mycut.first).c_str(),("h_met_"+mycut.first).c_str(),50,0,500));
        my_histos.emplace("h_mT_"+mycut.first,std::make_shared<TH1D>(("h_mT_"+mycut.first).c_str(),("h_mT_"+mycut.first).c_str(),40,0,200));
        my_histos.emplace("h_HT_"+mycut.first,std::make_shared<TH1D>(("h_HT_"+mycut.first).c_str(),("h_HT_"+mycut.first).c_str(),30,0,3000));
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
        const auto& runtype                 = tr.getVar<std::string>("runtype");
        const auto& filetag                 = tr.getVar<std::string>("filetag");
        const auto& JetID                   = tr.getVar<bool>("JetID");
        const auto& passTrigger             = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMuon         = tr.getVar<bool>("passTriggerMuon");
        const auto& passTriggerElectron     = tr.getVar<bool>("passTriggerElectron");

        const auto& Jets                    = tr.getVec<TLorentzVector>("Jets");
        const auto& Jets_CSV                = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const auto& NGoodJets_pt30          = tr.getVar<int>("NGoodJets_pt30");
        const auto& GoodJets_pt30           = tr.getVec<bool>("GoodJets_pt30");
        const auto& NGoodBJets_pt30         = tr.getVar<int>("NGoodBJets_pt30");
        const auto& GoodBJets_pt30          = tr.getVec<bool>("GoodBJets_pt30");
        const auto& NGoodBJets_pt30_loose   = tr.getVar<int>("NGoodBJets_pt30_loose");
        const auto& GoodBJets_pt30_loose    = tr.getVec<bool>("GoodBJets_pt30_loose");

        const auto& GoodLeptons         = tr.getVec< std::pair<std::string,TLorentzVector> >("GoodLeptons");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        
        const auto& Muons               = tr.getVec<TLorentzVector>("Muons");
        const auto& MuonsMTW            = tr.getVec<double>("MuonsMTW");
        const auto& GoodMuons           = tr.getVec<bool>("GoodMuons");
        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");
        const auto& Electrons           = tr.getVec<TLorentzVector>("Electrons");
        const auto& ElectronsMTW        = tr.getVec<double>("ElectronsMTW");
        const auto& GoodElectrons       = tr.getVec<bool>("GoodElectrons");
        const auto& NGoodElectrons      = tr.getVar<int>("NGoodElectrons");
        
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
        bool pass_1mu = (runtype != "Data" || filetag == "Data_SingleMuon")     && NGoodMuons == 1 && passTriggerMuon;
        bool pass_1el = (runtype != "Data" || filetag == "Data_SingleElectron") && NGoodElectrons == 1 && passTriggerElectron;
        bool pass_1l = pass_1mu || pass_1el;
        // compute mT 
        TLorentzVector metLV;
        metLV.SetPtEtaPhiE(MET,0,METPhi,MET);
        double mT = -1;
        //std::cout << "MET , METPhi: " << MET << ", " << METPhi << std::endl;
        //std::cout << "Nleptons: " << NGoodLeptons << std::endl;
        if (NGoodLeptons > 0)
        {
            //std::cout << "Found lepton, pt " << (GoodLeptons[0].second).Pt() << std::endl; 
            mT = utility::calcMT(GoodLeptons[0].second, metLV);
        }
        //std::cout << "mT: " << mT << std::endl;
        bool pass_mT = mT > 30 && mT < 100;

        // jet cuts
        bool pass_0b = NGoodBJets_pt30 == 0;
        bool pass_0b_loose = NGoodBJets_pt30_loose == 0;
        // compute dphi with highest pT jet
        double dphi = -1;
        if (NGoodJets_pt30 > 0)
        {
            for (int i = 0; i<Jets.size() ; ++i)
            {
                if(GoodJets_pt30[i])
                {
                    dphi = utility::calcDPhi(METPhi,Jets[i].Phi());
                    break;
                }
            }
        }
        bool pass_dphi = dphi > 0.5; 

        std::vector<double> Mbl_all;
        double Mbl_maxCSV = -1;
        double Mbl_maxpT = -1;
        if (pass_1l && NGoodJets_pt30 > 0)
        {
            TLorentzVector Jet_maxpT;
            for (int i = 0; i<Jets.size() ; ++i)
            {
                if(GoodJets_pt30[i])
                {
                    Jet_maxpT = Jets[i];
                    break;
                }
            }
            Mbl_maxpT = (GoodLeptons[0].second+Jet_maxpT).M();

            int maxCSV_index = -1;
            for (int i=0; i<Jets.size(); ++i)
            {
                if(!GoodJets_pt30[i]) continue;
                if(maxCSV_index == -1)
                    maxCSV_index = i;
                else if(Jets_CSV[i] > Jets_CSV[maxCSV_index])
                    maxCSV_index = i;

                double temp = (GoodLeptons[0].second+Jets[i]).M();
                Mbl_all.push_back(temp);
            }
            Mbl_maxCSV = (GoodLeptons[0].second+Jets[maxCSV_index]).M();
        }
        bool pass_Mbl_maxpT = Mbl_maxpT > 30 && Mbl_maxpT < 180;
        bool pass_Mbl_maxCSV = Mbl_maxCSV > 30 && Mbl_maxCSV < 180;
        bool pass_Mbl_all = false;
        for(double mbl : Mbl_all)
        {
            if (mbl > 30 && mbl < 180)
            {
                pass_Mbl_all = true;
                break;
            }
        }

        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map 
        {
            {"1l"                                 , goodEvent && pass_1l                                                             },
            {"1l_met"                             , goodEvent && pass_1l && MET>30                                                   },
            {"1l_mT"                              , goodEvent && pass_1l && pass_mT                                                  },
            {"1l_mT_met"                          , goodEvent && pass_1l && pass_mT && MET>30                                        },
                //{"1l_mT_0b"                           , goodEvent && pass_1l && pass_mT && pass_0b                                       },
                //{"1l_mT_0b_dphi"                      , goodEvent && pass_1l && pass_mT && pass_0b && pass_dphi                          },
                //{"1l_mT_0b_dphi_mbl_maxpT"            , goodEvent && pass_1l && pass_mT && pass_0b && pass_dphi && !pass_Mbl_maxpT        },
                //{"1l_mT_0b_dphi_mbl_maxCSV"           , goodEvent && pass_1l && pass_mT && pass_0b && pass_dphi && !pass_Mbl_maxCSV       },
                //{"1l_mT_0b_dphi_mbl_all"              , goodEvent && pass_1l && pass_mT && pass_0b && pass_dphi && !pass_Mbl_all          },
            {"1l_mT_0b_loose"                           , goodEvent && pass_1l && pass_mT && pass_0b_loose                                       },
            {"1l_mT_0b_loose_met"                       , goodEvent && pass_1l && pass_mT && pass_0b_loose && MET>30                             },
            {"1l_mT_0b_loose_met_dphi"                  , goodEvent && pass_1l && pass_mT && pass_0b_loose && pass_dphi && MET>30                         },
            {"1l_mT_0b_loose_met_dphi_mbl_maxpT"        , goodEvent && pass_1l && pass_mT && pass_0b_loose && pass_dphi && MET>30 && !pass_Mbl_maxpT        },
            {"1l_mT_0b_loose_met_dphi_mbl_maxCSV"       , goodEvent && pass_1l && pass_mT && pass_0b_loose && pass_dphi && MET>30 && !pass_Mbl_maxCSV       },
            {"1l_mT_0b_loose_met_dphi_mbl_all"          , goodEvent && pass_1l && pass_mT && pass_0b_loose && pass_dphi && MET>30 && !pass_Mbl_all          },
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
                my_histos["h_nbjets_loose_"  +kv.first]->Fill(NGoodBJets_pt30_loose, eventweight);
                my_histos["h_met_"     +kv.first]->Fill(MET, eventweight);
                my_histos["h_mT_"      +kv.first]->Fill(mT, eventweight);
                my_histos["h_HT_"      +kv.first]->Fill(HT, eventweight);
                my_histos["h_dphi_"    +kv.first]->Fill(dphi, eventweight);

                my_histos["h_lepton_pT_"   +kv.first]->Fill(GoodLeptons[0].second.Pt(), eventweight);
                my_histos["h_lepton_eta_"  +kv.first]->Fill(GoodLeptons[0].second.Eta(), eventweight);
                my_histos["h_lepton_phi_"  +kv.first]->Fill(GoodLeptons[0].second.Phi(), eventweight);

                my_histos["h_Mbl_maxCSV_"  +kv.first]->Fill(Mbl_maxCSV, eventweight);
                my_histos["h_Mbl_maxpT_"  +kv.first]->Fill(Mbl_maxpT, eventweight);
                for (double mbl: Mbl_all)
                    my_histos["h_Mbl_all_"  +kv.first]->Fill(mbl, eventweight);

                if (NGoodMuons > 0)
                {
                    for (int m=0; m<Muons.size(); ++m)
                    {
                        if(GoodMuons[m])
                        {
                            my_histos["h_muon_pT_"   +kv.first]->Fill(Muons[m].Pt(), eventweight);
                            my_histos["h_muon_eta_"  +kv.first]->Fill(Muons[m].Eta(), eventweight);
                            my_histos["h_muon_phi_"  +kv.first]->Fill(Muons[m].Phi(), eventweight);
                            break;
                        }
                    }
                } 
                if (NGoodElectrons > 0)
                {
                    for (int m=0; m<Electrons.size(); ++m)
                    {
                        if(GoodElectrons[m])
                        {
                            my_histos["h_electron_pT_"   +kv.first]->Fill(Electrons[m].Pt(), eventweight);
                            my_histos["h_electron_eta_"  +kv.first]->Fill(Electrons[m].Eta(), eventweight);
                            my_histos["h_electron_phi_"  +kv.first]->Fill(Electrons[m].Phi(), eventweight);
                            break;
                        }
                    }
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

