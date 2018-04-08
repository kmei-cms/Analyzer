#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>

//This is a helper function which will keep the plot from overlapping with the legend
void smartMax(const TH1 * const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax, const bool error)
{
    const bool isLog = p->GetLogy();
    double min = 9e99;
    double max = -9e99;
    double pThreshMax = -9e99;
    int threshold = static_cast<int>(h->GetNbinsX()*(l->GetX1() - p->GetLeftMargin())/((1 - p->GetRightMargin()) - p->GetLeftMargin()));

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        double bin = 0.0;
        if(error) bin = h->GetBinContent(i) + h->GetBinError(i);
        else      bin = h->GetBinContent(i);
        if(bin > max) max = bin;
        else if(bin > 1e-10 && bin < min) min = bin;
        if(i >= threshold && bin > pThreshMax) pThreshMax = bin;
    }

    gpThreshMax = std::max(gpThreshMax, pThreshMax);
    gmax = std::max(gmax, max);
    gmin = std::min(gmin, min);
}

//Class to hold TH1* with various helper functions 
class histInfo
{
public:
    std::string legEntry, histFile, histName, drawOptions;
    int color, rebin;
    std::shared_ptr<TH1> h;

    //helper function to get histogram from file and configure its optional settings
    void retrieveHistogram()
    {
        //Open the file for this histogram
        TFile *f = TFile::Open(histFile.c_str());

        //check that the file was opened successfully
        if(!f)
        {
            printf("File \"%s\" could not be opened!!!\n", histFile.c_str());
            h = nullptr;
            return;
        }

        //get the histogram from the file
        h.reset(static_cast<TH1*>(f->Get(histName.c_str())));

        //with the histogram retrieved, close the file
        f->Close();
        delete f;

        //check that the histogram was retireved from the file successfully
        if(!h)
        {
            printf("Histogram \"%s\" could not be found in file \"%s\"!!!\n", histName.c_str(), histFile.c_str());
            return;
        }

        //set the histogram color
        h->SetLineColor(color);
        h->SetLineWidth(3);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(20);

        // rebin the histogram if desired
        if(rebin >0) h->Rebin(rebin);
    }

    //helper function for axes
    void setupAxes()
    {
        h->SetStats(0);
        h->SetTitle(0);
        h->GetYaxis()->SetTitleOffset(1.2);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetXaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetLabelSize(0.045);
        if(h->GetXaxis()->GetNdivisions() % 100 > 5) h->GetXaxis()->SetNdivisions(6, 5, 0);
    }

    //wrapper to draw histogram
    void draw(const std::string& additionalOptions = "", bool noSame = false) const
    {
        h->Draw(((noSame?"":"same " + drawOptions + " " + additionalOptions)).c_str());
    }

    void setFillColor(int newColor = -1)
    {
        if(newColor >= 0) h->SetFillColor(newColor);
        else              h->SetFillColor(color);
    }

    void setLineStyle(int style = 1)
    {
        h->SetLineStyle(style);
    }

    histInfo(const std::string& legEntry, const std::string& histFile, const std::string& drawOptions, const int color) : legEntry(legEntry), histFile(histFile), histName(""), drawOptions(drawOptions), color(color), rebin(-1), h(nullptr)
    {
    }

    histInfo(TH1* h) : legEntry(h->GetName()), histFile(""), histName(h->GetName()), drawOptions(""), color(0), rebin(0), h(h)
    {
    }

    ~histInfo()
    {
    }
};

class Plotter
{
private:
    //entry for data
    std::vector<histInfo> data_;
    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries_;
    //vector summarizing signal histograms to include in the plot
    std::vector<histInfo> sigEntries_;
    
public:
    Plotter(std::vector<histInfo>&& data, std::vector<histInfo>&& bgEntries, std::vector<histInfo>&& sigEntries) : data_(data), bgEntries_(bgEntries), sigEntries_(sigEntries) {}

    void plot(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
    {
        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        //switch to the canvas to ensure it is the active object
        c->cd();

        //Create TLegend
        TLegend *leg = new TLegend(0.20, 0.76, 0.89, 0.88);
        //TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(3);
        leg->SetTextFont(42);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        //Create the THStack for the background ... warning, THStacks are terrible and must be filled "backwards"
        THStack *bgStack = new THStack();
        //Make seperate histogram from sum of BG histograms because I don't know how to make a THStack give me this 
        TH1* hbgSum = nullptr;
        for(int iBG = bgEntries_.size() - 1; iBG >= 0; --iBG)
        {
            //Get new histogram
            bgEntries_[iBG].histName = histName;
            bgEntries_[iBG].rebin = rebin;
            bgEntries_[iBG].retrieveHistogram();

            bgStack->Add(bgEntries_[iBG].h.get(), bgEntries_[iBG].drawOptions.c_str());
            if(!hbgSum) hbgSum = static_cast<TH1*>(bgEntries_[iBG].h->Clone());
            else        hbgSum->Add(bgEntries_[iBG].h.get());
        }

        //data
        //get new histogram from file
        for(auto& entry : data_)
        {
            //get new histogram
            entry.histName = histName;
            entry.rebin = rebin;
            entry.retrieveHistogram();

            //add histograms to TLegend
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), entry.drawOptions.c_str());
            smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, true);
        }

        //background
        for(auto& entry : bgEntries_)
        {
            //set fill color so BG will have solid fill
            entry.setFillColor();

            //add histograms to TLegend
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "F");
        }
        smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);

        //signal 
        for(auto& entry : sigEntries_)
        {
            //get new histogram
            entry.histName = histName;
            entry.rebin = rebin;
            entry.retrieveHistogram();
            entry.setLineStyle(2);

            //add histograms to TLegend
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L");
            smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        }

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.12);

        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hbgSum->GetBinLowEdge(1), hbgSum->GetBinLowEdge(hbgSum->GetNbinsX()) + hbgSum->GetBinWidth(hbgSum->GetNbinsX())));
        dummy.setupAxes();
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        dummy.h->SetTitle(histName.c_str());
        //Set the y-range of the histogram
        if(isLogY)
        {
            double locMin = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin = legSpan + log10(locMin);
            if(log10(lmax) > legMin)
            {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
            }
            dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
        }
        else
        {
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
            dummy.h->GetYaxis()->SetRangeUser(0.0, max*1.2);
        }
        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //plot background stack
        bgStack->Draw("same");

        //plot signal histograms
        for(const auto& entry : sigEntries_)
        {
            entry.draw();
        }

        //plot data histogram
        for(const auto& entry : data_)
        {
            entry.draw();
        }

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);

        TLatex mark;
        mark.SetNDC(true);

        //Draw CMS mark
        mark.SetTextAlign(11);
        mark.SetTextSize(0.050);
        mark.SetTextFont(61);
        //mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");
        mark.SetTextSize(0.040);
        mark.SetTextFont(52);
        //mark.DrawLatex(gPad->GetLeftMargin() + 0.11, 1 - (gPad->GetTopMargin() - 0.017), "Stealth 2018");

        //Draw lumistamp
        mark.SetTextFont(42);
        mark.SetTextAlign(31);
        //mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);

        //save new plot to file
        //c->Print(("outputPlots/" + histName + ".pdf").c_str());
        c->Print(("outputPlots/" + histName + ".png").c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
        delete bgStack;
        delete hbgSum;
    }
};

int main()
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry root file                 draw options  draw color
    std::vector<histInfo> data = {
        //{"All BG", "allBG.root"            , "PEX0",       kBlack}
    };

    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries = {
        {"T#bar{T}",        "condor/output-files/TT/TT.root",                           "hist", kBlue - 7   },
        {"QCD",             "condor/output-files/QCD/QCD.root",                         "hist", kGreen + 1  },
        {"ST",              "condor/output-files/ST/ST.root",                           "hist", kRed + 1    },
        {"WJetsToLNu",      "condor/output-files/WJetsToLNu/WJetsToLNu.root",           "hist", kYellow + 1 },
        {"Diboson",         "condor/output-files/Diboson/Diboson.root",                 "hist", kMagenta + 1},
        {"Rare",            "condor/output-files/Rare/Rare.root",                       "hist", kCyan + 1   },
        {"DYJetsToLL_M-50", "condor/output-files/DYJetsToLL_M-50/DYJetsToLL_M-50.root", "hist", kOrange + 2 },
    };

    //vector summarizing signal histograms to include in the plot
    std::vector<histInfo> sigEntries = {
        {"RPV 350", "condor/output-files/AllSignal/MyAnalysis_rpv_stop_350_0.root",         "hist", kMagenta + 2},
        {"SYY 650", "condor/output-files/AllSignal/MyAnalysis_stealth_stop_650_SYY_0.root", "hist", kGreen + 3  },
    };

    //make plotter object with the required sources for histograms specified
    Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));

    std::vector<std::string> mycuts_0l 
    {
        ""                     ,
        "g6j"                  ,
        "HT500"                ,
        "g2b"                  ,
        "1t"                   ,
        "2t"                   ,
        "g6j_HT500"            ,
        "g6j_HT500_g1b"        ,
        "g6j_HT500_g2b"        ,

        "g6j_HT500_g2b_1t"     ,
        "g6j_HT500_g2b_1t_f1"  , "g6j_HT500_g2b_1t_f2"  , "g6j_HT500_g2b_1t_f3"  , "g6j_HT500_g2b_1t_f4"  ,

        "g6j_HT500_g2b_1t1"    , "g6j_HT500_g2b_1t2"    , "g6j_HT500_g2b_1t3"    ,
        "g6j_HT500_g2b_1t1_f1" , "g6j_HT500_g2b_1t1_f2" , "g6j_HT500_g2b_1t1_f3" , "g6j_HT500_g2b_1t1_f4" ,
        "g6j_HT500_g2b_1t2_f1" , "g6j_HT500_g2b_1t2_f2" , "g6j_HT500_g2b_1t2_f3" , "g6j_HT500_g2b_1t2_f4" ,
        "g6j_HT500_g2b_1t3_f1" , "g6j_HT500_g2b_1t3_f2" , "g6j_HT500_g2b_1t3_f3" , "g6j_HT500_g2b_1t3_f4" ,

        "g6j_HT500_g2b_2t"     ,
        "g6j_HT500_g2b_2t_f1"  , "g6j_HT500_g2b_2t_f2"  , "g6j_HT500_g2b_2t_f3"  , "g6j_HT500_g2b_2t_f4"  ,

        "g6j_HT500_g2b_2t11"   , "g6j_HT500_g2b_2t12"   , "g6j_HT500_g2b_2t13"   , "g6j_HT500_g2b_2t22"   , "g6j_HT500_g2b_2t23", "g6j_HT500_g2b_2t33",
        "g6j_HT500_g2b_2t11_f1", "g6j_HT500_g2b_2t11_f2", "g6j_HT500_g2b_2t11_f3", "g6j_HT500_g2b_2t11_f4",
        "g6j_HT500_g2b_2t12_f1", "g6j_HT500_g2b_2t12_f2", "g6j_HT500_g2b_2t12_f3", "g6j_HT500_g2b_2t12_f4",
        "g6j_HT500_g2b_2t13_f1", "g6j_HT500_g2b_2t13_f2", "g6j_HT500_g2b_2t13_f3", "g6j_HT500_g2b_2t13_f4",
        "g6j_HT500_g2b_2t22_f1", "g6j_HT500_g2b_2t22_f2", "g6j_HT500_g2b_2t22_f3", "g6j_HT500_g2b_2t22_f4",
        "g6j_HT500_g2b_2t23_f1", "g6j_HT500_g2b_2t23_f2", "g6j_HT500_g2b_2t23_f3", "g6j_HT500_g2b_2t23_f4",
        "g6j_HT500_g2b_2t33_f1", "g6j_HT500_g2b_2t33_f2", "g6j_HT500_g2b_2t33_f3", "g6j_HT500_g2b_2t33_f4",
    };

    for(std::string mycut : mycuts_0l)
    {
        plt.plot( "h_njets_0l_"+mycut, "N_{J}" , "Events", true);
        plt.plot( "h_ntops_0l_"+mycut, "N_{T}" , "Events", true);
        plt.plot( "h_nb_0l_"   +mycut, "N_{B}" , "Events", true);        
        plt.plot( "h_HT_0l_"   +mycut, "H_{T}" , "Events", true);        
    }

    plt.plot("h_met"     , "MET"   , "Events", true);
    plt.plot("h_ht"      , "H_{T}" , "Events", true);
    plt.plot("h_bdt"     , "bdt"   , "Events", true);
    plt.plot("h_fisher"  , "fisher", "Events", true);
    plt.plot("h_njets"   , "N_{J}" , "Events", true);
    plt.plot("h_nb"      , "N_{B}" , "Events", true);
    plt.plot("h_ntops"   , "N_{T}" , "Events", true);
    plt.plot("h_ntops_j1", "N_{1T}", "Events", true);
    plt.plot("h_ntops_j2", "N_{2T}", "Events", true);
    plt.plot("h_ntops_j3", "N_{3T}", "Events", true);
   
    //plt.plot("HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
}
