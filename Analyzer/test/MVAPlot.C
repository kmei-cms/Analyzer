#include "Analyzer/Analyzer/test/HistInfoCollection.h"
#include "plotter.h"
#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TF1.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <map>

Double_t gaufun(Double_t* x, Double_t* par) {
    //par[0]=mean, par[1]=sigma, par[2]=amplitude  
    return par[2]*TMath::Gaus(x[0],par[0],par[1]);
}

std::unique_ptr<TH1> AddedHistos(const std::string& histName, std::vector<histInfo>& bgEntries)
{
    std::unique_ptr<TH1> hTotal; 
    bool firstPass = true;
    for(auto& entry : bgEntries)
    {
        entry.histName = histName;
        entry.retrieveHistogram();
        if(firstPass)
        {
            hTotal.reset( static_cast<TH1*>(entry.h->Clone()) );
            firstPass = false;
        }
        else 
        {
            hTotal->Add(entry.h.get());
        }
    }
    return hTotal;
}

template<typename H> std::unique_ptr<H> GetHisto(const std::string& histName, histInfo& entry)
{
    std::unique_ptr<H> h; 
    entry.histName = histName;
    entry.retrieveHistogram();
    h.reset( static_cast<H*>(entry.h->Clone()) );
    h->SetFillColor(entry.color);
    return h;    
}

void MVAPlot()
{
    TH1::AddDirectory(false);
    //std::string path = "condor/oldTest/Analyze1Lep_Kerasv1.2.5";
    //td::string year = "2016";
    std::string path = "condor/oldTest/Analyze1Lep_Kerasv3.0.1";
    std::string year = "2017";

    std::vector<histInfo> bgEntries = {
        //{"DYJetsToLL_M-50", path + "/"+year+"_DYJetsToLL_M-50.root", "hist", kOrange + 2 },        
        //{"Rare",            path + "/"+year+"_Rare.root",            "hist", kCyan + 1   },
        //{"Diboson",         path + "/"+year+"_Diboson.root",         "hist", kMagenta + 1},
        //{"WJetsToLNu",      path + "/"+year+"_WJetsToLNu.root",      "hist", kYellow + 1 },
        //{"ST",              path + "/"+year+"_ST.root",              "hist", kRed + 1    },
        //{"QCD",             path + "/"+year+"_QCD.root",             "hist", kGreen + 1  },
        {"TT",        path + "/"+year+"_TT.root",              "hist", kBlue       },
    };

    //vector summarizing signal histograms to include in the plot
    std::vector<histInfo> sigEntries;
    if(year == "2016")
    {        
        sigEntries = {
            {"RPV 850", path + "/"+year+"_rpv_stop_850.root",               "hist", kRed      },
            {"SYY 650", path + "/"+year+"_stealth_stop_650_SYY.root",       "hist", kMagenta  },
            {"RPV 350", path + "/"+year+"_rpv_stop_350.root",               "hist", kCyan     },        
        };
    }
    else
    {
        sigEntries = {        
            {"RPV 850", path + "/"+year+"_RPV_2t6j_mStop-850.root",         "hist", kRed     },
            {"SYY 650", path + "/"+year+"_StealthSYY_2t6j_mStop-650.root",  "hist", kMagenta },
            {"RPV 350", path + "/"+year+"_RPV_2t6j_mStop-350.root",         "hist", kCyan    },        
        };
    }

    std::vector<std::string> histNameVec = {
        "h_deepESM_1l_6j_ge1b","h_deepESM_1l_7j_ge1b",
        "h_deepESM_1l_8j_ge1b","h_deepESM_1l_9j_ge1b",
        "h_deepESM_1l_10j_ge1b","h_deepESM_1l_11j_ge1b",
        "h_deepESM_1l_12j_ge1b","h_deepESM_1l_13j_ge1b",
        "h_deepESM_1l_14j_ge1b","h_deepESM_1l_15j_ge1b"
    };

    //Make NJet plots of MVA distributions and get info for correction
    TCanvas* c = new TCanvas("c", "c", 800, 800);
    c->cd();

    //std::vector<double> means; double wSum = 0; double nEvts = 0;
    //for(const auto& name : histNameVec)
    //{
    //    std::unique_ptr<TH1> histo = AddedHistos(name, bgEntries);
    //    double mean = histo->GetMean();
    //    double nE   = histo->GetEntries();
    //    wSum += nE*mean;
    //    nEvts += nE;
    //    means.push_back(mean);
    //    histo->Scale(1.0 / histo->Integral() );
    //    histo->Draw("hist");        
    //    c->Print(("all_MVA_"+name+".png").c_str());    
    //}
    //
    ////Print out NJet corrections for MVA training
    //double wAvg = wSum/nEvts;
    //std::cout<<"Weighted Avg: "<<wAvg<<std::endl;
    //std::cout<<"------------Corrections----------"<<std::endl;
    //int i = 5;
    //for(const auto& m : means)
    //{
    //    i++;
    //    std::cout<<"{"<<i<<", "<<wAvg - m<<"},"<<std::endl;
    //}
    //
    //Get the baseline distribution and define the 4 MVA bins
    //std::string baselineName = "h_fisher_0l_ge6j_HT500_ge2b";
    //std::string baselineName = "h_deepESM_1l_ge6j_ge1b";
    //std::unique_ptr<TH1> baseline = AddedHistos(baselineName, bgEntries);
    //std::unique_ptr<TH1> integralBaseline;
    //integralBaseline.reset(static_cast<TH1*>(baseline->Clone("Integral")));
    //std::cout<<"Might have to tune the ratio bounds until you get a value"<<std::endl;
    //for(int i = 0; i < baseline->GetNbinsX(); i++)
    //{
    //    double integral = baseline->Integral(0,i);
    //    double iMax = baseline->Integral(0,baseline->GetNbinsX());
    //    double ratio = integral/iMax;
    //    //std::cout<<integral<<" "<<iMax<<" "<<ratio<<std::endl;
    //    if(abs(ratio - 0.25) < 0.005) std::cout<<"1/4 of Events: "<<baseline->GetBinCenter(i)<<std::endl;
    //    if(abs(ratio - 0.50) < 0.005) std::cout<<"2/4 of Events: "<<baseline->GetBinCenter(i)<<std::endl;
    //    if(abs(ratio - 0.75) < 0.005) std::cout<<"3/4 of Events: "<<baseline->GetBinCenter(i)<<std::endl;
    //    integralBaseline->SetBinContent(i, ratio);
    //}
    //integralBaseline->Draw("hist");
    //c->Print(("Integral_"+baselineName+".png").c_str());    
    
    //Make plots to show NJet dependance 
    std::vector<std::string> variables = {"deepESM"};
    //std::vector<std::string> variables = {"njets_deepESM", "njets_fisher"};
    //std::vector<std::string> variables = {
    //    //"Jet_cm_pt_1_deepESM", "Jet_cm_pt_2_deepESM", "Jet_cm_pt_3_deepESM", "Jet_cm_pt_4_deepESM", "Jet_cm_pt_5_deepESM", "Jet_cm_pt_6_deepESM", "Jet_cm_pt_7_deepESM",
    //    //"Jet_cm_eta_1_deepESM", "Jet_cm_eta_2_deepESM", "Jet_cm_eta_3_deepESM", "Jet_cm_eta_4_deepESM", "Jet_cm_eta_5_deepESM", "Jet_cm_eta_6_deepESM", "Jet_cm_eta_7_deepESM",
    //    //"Jet_cm_phi_1_deepESM", "Jet_cm_phi_2_deepESM", "Jet_cm_phi_3_deepESM", "Jet_cm_phi_4_deepESM", "Jet_cm_phi_5_deepESM", "Jet_cm_phi_6_deepESM", "Jet_cm_phi_7_deepESM",
    //    //"Jet_cm_m_1_deepESM", "Jet_cm_m_2_deepESM", "Jet_cm_m_3_deepESM", "Jet_cm_m_4_deepESM", "Jet_cm_m_5_deepESM", "Jet_cm_m_6_deepESM", "Jet_cm_m_7_deepESM",
    //
    //    "Jet_cm_pt_1_njets", "Jet_cm_pt_2_njets", "Jet_cm_pt_3_njets", "Jet_cm_pt_4_njets", "Jet_cm_pt_5_njets", "Jet_cm_pt_6_njets", "Jet_cm_pt_7_njets",
    //    "Jet_cm_eta_1_njets", "Jet_cm_eta_2_njets", "Jet_cm_eta_3_njets", "Jet_cm_eta_4_njets", "Jet_cm_eta_5_njets", "Jet_cm_eta_6_njets", "Jet_cm_eta_7_njets",
    //    "Jet_cm_phi_1_njets", "Jet_cm_phi_2_njets", "Jet_cm_phi_3_njets", "Jet_cm_phi_4_njets", "Jet_cm_phi_5_njets", "Jet_cm_phi_6_njets", "Jet_cm_phi_7_njets",
    //    "Jet_cm_m_1_njets", "Jet_cm_m_2_njets", "Jet_cm_m_3_njets", "Jet_cm_m_4_njets", "Jet_cm_m_5_njets", "Jet_cm_m_6_njets", "Jet_cm_m_7_njets",
    //};
    //std::vector<std::string> variables = {
    //    "Jet_cm_pt_1", "Jet_cm_pt_2", "Jet_cm_pt_3", "Jet_cm_pt_4", "Jet_cm_pt_5", "Jet_cm_pt_6", "Jet_cm_pt_7",
    //    "Jet_cm_eta_1", "Jet_cm_eta_2", "Jet_cm_eta_3", "Jet_cm_eta_4", "Jet_cm_eta_5", "Jet_cm_eta_6", "Jet_cm_eta_7",
    //    "Jet_cm_phi_1", "Jet_cm_phi_2", "Jet_cm_phi_3", "Jet_cm_phi_4", "Jet_cm_phi_5", "Jet_cm_phi_6", "Jet_cm_phi_7",
    //    "Jet_cm_m_1", "Jet_cm_m_2", "Jet_cm_m_3", "Jet_cm_m_4", "Jet_cm_m_5", "Jet_cm_m_6", "Jet_cm_m_7",
    //    "fwm2_top6","fwm3_top6","fwm4_top6","fwm5_top6","jmt_ev0_top6","jmt_ev1_top6","jmt_ev2_top6",
    //};
    //std::vector<std::string> variables = {"deepESM"};
    std::vector<histInfo> entries = sigEntries; 
    entries.push_back(bgEntries[0]);
    for(const auto& v : variables)
    {
        for(auto& entry : entries)
        {
            //std::unique_ptr<TH1F> dummy(new TH1F("Dummy","dummy", 10, 0, 20));
            //dummy->SetMaximum(1.1);
            //dummy->SetMinimum(-0.1);
            //dummy->Draw();
            std::unique_ptr<TH2> h2d = GetHisto<TH2>("h_njets_"+v+"_1l_ge7j_ge1b", entry);
            //h2d->GetYaxis()->SetRangeUser(-10, 10);
            //h2d->GetXaxis()->SetLimits(0, 1);
            h2d->Draw("COLZ");
            //h2d->RebinX(2);
            h2d->RebinY(10);
            h2d->GetYaxis()->SetTitle(v.c_str());
            h2d->GetXaxis()->SetTitle("N_{J}");
            h2d->SetTitle( (entry.legEntry+": "+entry.histName).c_str() );
            //std::unique_ptr<TProfile> hprofile = GetHisto<TProfile>("hTp_njets_"+v+"_1l_ge6j_ge1b", entry);
            std::unique_ptr<TProfile> hprofile( h2d->ProfileX(("hTp_"+v+"_1l_ge7j_ge1b").c_str()) );
            hprofile->Draw("same");
            hprofile->SetLineColor(kRed);
            hprofile->SetMarkerColor(kRed);
        
            gStyle->SetPalette(kRainBow);
            gStyle->SetStatY(0.89);
            gStyle->SetStatX(0.35);
            //gStyle->SetStatW(0.3);
            //gStyle->SetStatH(0.2); 
        
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.15);
            gPad->SetTopMargin(0.08);
            gPad->SetBottomMargin(0.12);
            gPad->SetTicks(1,1);
            gPad->SetLogz();
            c->Print((entry.legEntry+"_"+v+"_nJetDependance_"+year+".pdf").c_str());            
        }
    
        //TLegend *leg = new TLegend(0.20, 0.76, 0.89, 0.88);
        //leg->SetFillStyle(0);
        //leg->SetBorderSize(0);
        //leg->SetLineWidth(1);
        //leg->SetNColumns(3);
        //leg->SetTextFont(42);
        //
        //double norm = 1.0;
        //std::unique_ptr<TH1> h1 = GetHisto<TH1>(v+"_1l_ge7j_ge1b", bgEntries[0]);
        ////h1->Rebin(10);
        //h1->Scale(norm / h1->Integral() );
        //h1->SetStats(0);
        //h1->GetXaxis()->SetTitle((v).c_str());
        //h1->GetYaxis()->SetTitle("A.U.");            
        //double h1max = h1->GetBinContent(h1->GetMaximumBin()) * 1.3;
        //h1->Draw("hist");
        //leg->AddEntry(h1.get(), ("TT_"+v).c_str(), "F");
        //
        //std::vector<std::unique_ptr<TH1>> hVec;
        //for(auto& entry : sigEntries)
        //{        
        //    hVec.emplace_back( GetHisto<TH1>(v+"_1l_ge7j_ge1b", entry) );            
        //    //hVec.back()->Rebin(10);
        //    hVec.back()->Scale(norm / hVec.back()->Integral() );
        //    double h2max = hVec.back()->GetBinContent(hVec.back()->GetMaximumBin()) * 1.3;
        //    hVec.back()->SetFillStyle(3005);
        //    hVec.back()->Draw("hist same");
        //    leg->AddEntry(hVec.back().get(), (entry.legEntry+"_"+v).c_str(), "F");
        //    h1max = std::max(h1max, h2max);
        //
        //    //std::unique_ptr<TH1> h3 = GetHisto<TH1>("h_"v+"_1l_ge6j_ge1b", bgEntries[1]);
        //    //h3->Rebin(1);
        //    //h3->Scale(norm / h3->Integral() );
        //    //double h3max = h3->GetBinContent(h3->GetMaximumBin()) * 1.3;
        //    //h3->Draw("PEX0 same");
        //    //leg->AddEntry(h3.get(), ("Data_"+v).c_str(), "F");
        //    //h1->SetMaximum( std::max(std::max(h1max, h2max), h3max) );
        //}
        //
        //leg->Draw();
        //h1->SetMaximum( h1max );
        ////h1->Draw("hist same");
        //c->Print((year+"_"+v+"_.png").c_str());
    }
    //
    ////Make input variable plots
    ////std::vector<std::string> trainingVars = {"fwm2_top6","fwm3_top6","fwm4_top6","fwm5_top6","fwm6_top6","fwm7_top6","fwm8_top6","fwm9_top6","fwm10_top6",
    ////                                         "jmt_ev0_top6","jmt_ev1_top6","jmt_ev2_top6", "h_BestComboMass_1l_ge6j_ge1b"};
    //std::vector<std::string> trainingVars = {"h_BestComboMass_1l_ge6j_ge1b"};
    //for(const auto& v : trainingVars)
    //{
    //    TLegend *leg = new TLegend(0.65, 0.80, 0.89, 0.88);
    //    leg->SetFillStyle(0);
    //    leg->SetBorderSize(0);
    //    leg->SetLineWidth(1);
    //    leg->SetNColumns(1);
    //    leg->SetTextFont(42);
    //
    //    double norm = 1.0;
    //    std::unique_ptr<TH1> hBG = GetHisto<TH2>(v, bgEntries[0]);
    //    hBG->Draw("hist");
    //    hBG->SetStats(0);
    //    hBG->SetFillColor(kBlue - 9);
    //    hBG->Scale(norm / hBG->Integral() );
    //    hBG->SetMaximum( hBG->GetBinContent(hBG->GetMaximumBin()) * 1.3);
    //    hBG->GetYaxis()->SetTitle("A.U.");
    //    hBG->GetXaxis()->SetTitle(v.c_str());
    //    hBG->SetTitle( ("Input variable: "+v).c_str() );
    //    leg->AddEntry(hBG.get() , "TT", "F");
    //    std::unique_ptr<TH1> hSig = AddedHistos(v, sigEntries);
    //    hSig->Scale(norm / hSig->Integral() );
    //    hSig->Draw("hist same");
    //    hSig->SetLineColor(kRed);
    //    hSig->SetMarkerColor(kRed);
    //    hSig->SetFillColor(kRed);
    //    hSig->SetFillStyle(3005);
    //    leg->AddEntry(hSig.get() , "Signal", "F");
    //    leg->Draw();
    //
    //    gStyle->SetStatY(0.89);
    //    gStyle->SetStatX(0.35);
    //    //gStyle->SetStatW(0.3);
    //    //gStyle->SetStatH(0.2); 
    //
    //    gPad->SetLeftMargin(0.12);
    //    gPad->SetRightMargin(0.15);
    //    gPad->SetTopMargin(0.08);
    //    gPad->SetBottomMargin(0.12);
    //    gPad->SetTicks(1,1);
    //    c->Print((v+".png").c_str());                    
    //}

    delete c;
}

int main()
{
    MVAPlot();
}
