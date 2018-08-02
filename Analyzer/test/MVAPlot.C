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
    return h;    
}

void MVAPlot()
{
    TH1::AddDirectory(false);

    //std::string path = "lep0Ana_TestFisherV8-May-30-2018";
    //std::string path = "deepESM_v1";
    //std::string path = "deepESM_v2";
    //std::string path = "deepESM_GRfalse_1Layer";
    std::string path = "deepESM_GRtrue_1Layer";
    
    std::vector<histInfo> bgEntries = {
        //{"DYJetsToLL_M-50", "condor/output-files/" + path + "/DYJetsToLL_M-50/DYJetsToLL_M-50.root", "hist", kBlack      },        
        //{"Rare",            "condor/output-files/" + path + "/Rare/Rare.root",                       "hist", kCyan + 1   },
        //{"Diboson",         "condor/output-files/" + path + "/Diboson/Diboson.root",                 "hist", kMagenta + 1},
        //{"WJetsToLNu",      "condor/output-files/" + path + "/WJetsToLNu/WJetsToLNu.root",           "hist", kYellow + 1 },
        //{"ST",              "condor/output-files/" + path + "/ST/ST.root",                           "hist", kRed + 1    },
        //{"QCD",             "condor/output-files/" + path + "/QCD/QCD.root",                         "hist", kGreen + 1  },
        {"TT",        "condor/output-files/" + path + "/TT/TT.root",                           "hist", kBlue - 7   },
    };

    //vector summarizing signal histograms to include in the plot
    std::vector<histInfo> sigEntries = {
        {"RPV_350", "condor/output-files/" + path + "/AllSignal/MyAnalysis_rpv_stop_350_0.root",         "hist", kMagenta + 2},
        {"SYY_650", "condor/output-files/" + path + "/AllSignal/MyAnalysis_stealth_stop_650_SYY_0.root", "hist", kGreen + 3  },
        {"RPV_850", "condor/output-files/" + path + "/AllSignal/MyAnalysis_rpv_stop_850_0.root",         "hist", kRed + 1    },
    };


    //std::vector<std::string> histNameVec = {
    //    "h_fisher_0l_6j_HT500_ge2b","h_fisher_0l_7j_HT500_ge2b",
    //    "h_fisher_0l_8j_HT500_ge2b","h_fisher_0l_9j_HT500_ge2b",
    //    "h_fisher_0l_10j_HT500_ge2b","h_fisher_0l_11j_HT500_ge2b",
    //    "h_fisher_0l_12j_HT500_ge2b","h_fisher_0l_13j_HT500_ge2b",
    //    "h_fisher_0l_14j_HT500_ge2b","h_fisher_0l_15j_HT500_ge2b"
    //};
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

    std::vector<double> means; double wSum = 0; double nEvts = 0;
    for(const auto& name : histNameVec)
    {
        std::unique_ptr<TH1> histo = AddedHistos(name, bgEntries);
        double mean = histo->GetMean();
        double nE   = histo->GetEntries();
        wSum += nE*mean;
        nEvts += nE;
        means.push_back(mean);
        histo->Draw("hist");        
        c->Print(("all_MVA_"+name+".png").c_str());    
    }
    
    //Print out NJet corrections for MVA training
    double wAvg = wSum/nEvts;
    std::cout<<"Weighted Avg: "<<wAvg<<std::endl;
    std::cout<<"------------Corrections----------"<<std::endl;
    int i = 5;
    for(const auto& m : means)
    {
        i++;
        std::cout<<"{"<<i<<", "<<wAvg - m<<"},"<<std::endl;
    }
    
    //Get the baseline distribution and define the 4 MVA bins
    //std::string baselineName = "h_fisher_0l_ge6j_HT500_ge2b";
    std::string baselineName = "h_deepESM_1l_ge6j_ge1b";
    std::unique_ptr<TH1> baseline = AddedHistos(baselineName, bgEntries);
    std::unique_ptr<TH1> integralBaseline;
    integralBaseline.reset(static_cast<TH1*>(baseline->Clone("Integral")));
    std::cout<<"Might have to tune the ratio bounds until you get a value"<<std::endl;
    for(int i = 0; i < baseline->GetNbinsX(); i++)
    {
        double integral = baseline->Integral(0,i);
        double iMax = baseline->Integral(0,baseline->GetNbinsX());
        double ratio = integral/iMax;
        if(0.240 < ratio && ratio < 0.260) std::cout<<"1/4 of Events: "<<baseline->GetBinCenter(i)<<std::endl;
        if(0.495 < ratio && ratio < 0.505) std::cout<<"1/2 of Events: "<<baseline->GetBinCenter(i)<<std::endl;
        if(0.745 < ratio && ratio < 0.755) std::cout<<"3/4 of Events: "<<baseline->GetBinCenter(i)<<std::endl;
        integralBaseline->SetBinContent(i, ratio);
    }
    integralBaseline->Draw("hist");
    c->Print(("Integral_"+baselineName+".png").c_str());    

    //Make plots to show NJet dependance 
    std::vector<std::string> variables = {"deepESM", "fisher"};
    //std::vector<std::string> variables = {"deepESM"};
    std::vector<histInfo> entries = sigEntries; 
    entries.push_back(bgEntries[0]);
    for(const auto& v : variables)
    {
        for(auto& entry : entries)
        {
            std::unique_ptr<TH2> h2d = GetHisto<TH2>("h_njets_"+v+"_1l_ge6j_ge1b", entry);
            h2d->Draw("COLZ");
            h2d->RebinY(4);
            h2d->GetYaxis()->SetTitle((v+" Discriminator").c_str());
            h2d->GetXaxis()->SetTitle("N_{J}");
            h2d->SetTitle( (entry.legEntry+": "+entry.histName).c_str() );
            std::unique_ptr<TProfile> hprofile = GetHisto<TProfile>("hTp_njets_"+v+"_1l_ge6j_ge1b", entry);
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
            c->Print((entry.legEntry+"_"+v+"_nJetDependance.png").c_str());            
        }

        for(auto& entry : sigEntries)
        {
            TLegend *leg = new TLegend(0.20, 0.76, 0.89, 0.88);
            leg->SetFillStyle(0);
            leg->SetBorderSize(0);
            leg->SetLineWidth(1);
            leg->SetNColumns(3);
            leg->SetTextFont(42);

            double norm = 1.0;
            std::unique_ptr<TH1> h1 = GetHisto<TH1>("h_"+v+"_1l_ge6j_ge1b", entry);
            h1->SetStats(0);
            h1->GetXaxis()->SetTitle((v+" Discriminator").c_str());
            h1->GetYaxis()->SetTitle("Events");            
            h1->Scale(norm / h1->Integral() );
            h1->SetMaximum( h1->GetBinContent(h1->GetMaximumBin()) * 1.3);
            h1->Draw("hist");
            leg->AddEntry(h1.get(), (entry.legEntry+"_"+v).c_str(), "F");

            std::unique_ptr<TH1> h2 = GetHisto<TH1>("h_"+v+"_1l_ge6j_ge1b", bgEntries[0]);
            h2->Scale(norm / h2->Integral() );
            h2->Draw("hist same");
            leg->AddEntry(h2.get(), ("TT_"+v).c_str(), "F");
            leg->Draw();
            c->Print((entry.legEntry+"_"+v+"_.png").c_str());            
        }
    }

    //Make Roc Plots
    //std::unique_ptr<TH1> h = GetHisto<TH1>("h_deepESM_1l_ge6j_ge1b", sigEntries[0]);
    //rocInfo deepESMRocInfo = { makeFisherVec(std::move(h)), sigEntries[0].legEntry, sigEntries[0].color};


    delete c;
}

int main()
{
    MVAPlot();
}
