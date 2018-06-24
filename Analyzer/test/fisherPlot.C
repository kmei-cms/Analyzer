#include "Analyzer/Analyzer/test/HistInfoCollection.h"
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

void fisherPlot()
{
    TH1::AddDirectory(false);

    //std::string path = "lep0Ana_TestFisherV8-May-30-2018";
    std::string path = "lep0Ana_TestFisherV8_Corr-June-3-2018";
    
    std::vector<histInfo> bgEntries = {
        {"DYJetsToLL_M-50", "condor/output-files/" + path + "/DYJetsToLL_M-50/DYJetsToLL_M-50.root", "hist", kBlack      },        
        {"Rare",            "condor/output-files/" + path + "/Rare/Rare.root",                       "hist", kCyan + 1   },
        {"Diboson",         "condor/output-files/" + path + "/Diboson/Diboson.root",                 "hist", kMagenta + 1},
        {"WJetsToLNu",      "condor/output-files/" + path + "/WJetsToLNu/WJetsToLNu.root",           "hist", kYellow + 1 },
        {"ST",              "condor/output-files/" + path + "/ST/ST.root",                           "hist", kRed + 1    },
        {"QCD",             "condor/output-files/" + path + "/QCD/QCD.root",                         "hist", kGreen + 1  },
        {"T#bar{T}",        "condor/output-files/" + path + "/TT/TT.root",                           "hist", kBlue - 7   },
    };

    std::vector<std::string> histNameVec = {
        "h_fisher_0l_6j_HT500_ge2b","h_fisher_0l_7j_HT500_ge2b",
        "h_fisher_0l_8j_HT500_ge2b","h_fisher_0l_9j_HT500_ge2b",
        "h_fisher_0l_10j_HT500_ge2b","h_fisher_0l_11j_HT500_ge2b",
        "h_fisher_0l_12j_HT500_ge2b","h_fisher_0l_13j_HT500_ge2b",
        "h_fisher_0l_14j_HT500_ge2b","h_fisher_0l_15j_HT500_ge2b"
    };

    //Make NJet plots of Fisher distributions and get info for correction
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
        c->Print(("all_fisher_"+name+".png").c_str());    
    }

    //Print out NJet corrections for fisher training
    double wAvg = wSum/nEvts;
    std::cout<<"Weighted Avg: "<<wAvg<<std::endl;
    std::cout<<"------------Corrections----------"<<std::endl;
    int i = 5;
    for(const auto& m : means)
    {
        i++;
        std::cout<<"{"<<i<<", "<<wAvg - m<<"},"<<std::endl;
    }

    //Get the baseline distribution and define the 4 fisher bins
    std::string baselineName = "h_fisher_0l_ge6j_HT500_ge2b";
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

    delete c;
}

int main()
{
    fisherPlot();
}
