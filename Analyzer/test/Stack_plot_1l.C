#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <algorithm>


// -----------------------------------------------------------------------------
// smartMax FUNCTION: to keep the plot from overlapping with the legend 
// -----------------------------------------------------------------------------


void smartMax(const TH1 * const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax, const bool error)
{
    //const bool isLog  = p->GetLogy();
    double min    = 9e99;
    double max    = -9e99;
    double pThreshMax = -9e99;
    int threshold     = static_cast<int>(h->GetNbinsX()*(l->GetX1() - p->GetLeftMargin())/((1 - p->GetRightMargin()) - p->GetLeftMargin()));

    for(int i = 1; i <= h->GetNbinsX(); ++i) {

        double bin = 0.0;
        
        if(error) 
            bin = h->GetBinContent(i) + h->GetBinError(i);
        else      
            bin = h->GetBinContent(i);
        if(bin > max) 
            max = bin;
        else if(bin > 1e-10 && bin < min) 
            min = bin;
        if(i >= threshold && bin > pThreshMax) 
            pThreshMax = bin;
    }

    gpThreshMax = std::max(gpThreshMax, pThreshMax);
    gmax    = std::max(gmax, max);
    gmin    = std::min(gmin, min);

}


// -----------------------------------------------------------------------------
// histInfo CLASS: to hold TH1* with various helper functions 
// -----------------------------------------------------------------------------
class histInfo
{
public:
    std::string legName, legEntry, histFile, histName, drawOptions;
    int color, rebin;
    double nEvents;
    std::shared_ptr<TH1> h;
    bool drawHisto; // define for empty data plot

    // -----------------------------------------------------------------------------------------------
    // retrieveHistogram FUNCTION:  to get histogram from file and configure its optional settings
    // -----------------------------------------------------------------------------------------------
    void retrieveHistogram()
    {
        if(drawHisto) // define for empty data plot
        {
            //Open the file for this histogram
            TFile *f = TFile::Open(histFile.c_str());

            if(!f) {
                printf("File \"%s\" could not be opened!!!\n", histFile.c_str());
                h = nullptr;
                return;
            }

            //get histogram & close file
            h.reset(static_cast<TH1*>(f->Get(histName.c_str())));
            f->Close();
            delete f;

            if(!h) {
                printf("Histogram \"%s\" could not be found in file \"%s\"!!!\n", histName.c_str(), histFile.c_str());
                return;
            }

            h->SetLineColor(color);
            h->SetLineWidth(3);
            h->SetMarkerColor(color);
            h->SetMarkerStyle(20);

            //Get number of events and save it as a member variable
            nEvents  = h->Integral();
            legEntry = legName + " ("+std::to_string(nEvents)+") ";     

            // rebin the histogram if desired
            if(rebin >0) 
                h->Rebin(rebin);
        } 
    }

    // --------------------------------------------
    // setupAxes FUNCTION: to help for axes
    // -------------------------------------------- 
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
        if(h->GetXaxis()->GetNdivisions() % 100 > 5) 
            h->GetXaxis()->SetNdivisions(6, 5, 0);
    }

    // ------------------------------------------------
    // draw FUNCTION: wrapper to draw histogram
    // ------------------------------------------------
    void draw(const std::string& additionalOptions = "", bool noSame = false) const 
    {
        h->Draw(((noSame?"":"same " + drawOptions + " " + additionalOptions)).c_str());
    }

    // ------------------------------------------------
    // clone FUNCTION: wrapper to clone histogram
    // ------------------------------------------------
    TH1F* clone(const std::string& newname = "0")
    {
        return (TH1F*)h->Clone(newname.c_str());
    }

    //------------------------------------------------
    // getHisto FUNCTION: get the histo from the class
    // -----------------------------------------------
    std::shared_ptr<TH1> getHisto()
    {
        return h;
    }

    // -------------------------------
    // setFillColor FUNCTION:
    // -------------------------------
    void setFillColor(int newColor = -1) {
        if(newColor >= 0) 
            h->SetFillColor(newColor);
        else          
            h->SetFillColor(color);
    }

    histInfo(const std::string& legName, const std::string& histFile, const std::string& drawOptions, const int color, const bool drawHisto = true) : legName(legName), legEntry(legName), histFile(histFile), histName(""), drawOptions(drawOptions), color(color), rebin(-1), nEvents(-1), h(nullptr), drawHisto(drawHisto)
    {
    }

    histInfo(TH1* h) : legName(h->GetName()), legEntry(h->GetName()), histFile(""), histName(h->GetName()), drawOptions(""), color(0), rebin(0), nEvents(-1), h(h)
    {
    }

    ~histInfo()
    {
    }
};


// -----------------------------------------------------------------------------
// Plotter CLASS: 
// -----------------------------------------------------------------------------
class Plotter
{

private:
    //entry for data
    histInfo data_;
    //vector summarizing background & signal histograms to include in the plot
    std::vector<histInfo> bgEntries_;
    std::vector<histInfo> sigEntries_;

    static bool compareNEvents(histInfo h1, histInfo h2) 
    {
        const double nE1 = h1.nEvents;
        const double nE2 = h2.nEvents;  
        return (nE1 > nE2);
    }
    
public:
    //Plotter(histInfo&& data, std::vector<histInfo>&& bgEntries, std::vector<histInfo>&& sigEntries) : data_(data), bgEntries_(bgEntries), sigEntries_(sigEntries) {}
    Plotter(histInfo&& data, std::vector<histInfo>&& bgEntries) : data_(data), bgEntries_(bgEntries) {}

    // -------------------------------
    // plot FUNCTION:
    // -------------------------------   
    void plot(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const std::string& cutlabel = "", const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 39500) // lumi 2016 = 39500, lumi 2017= 41500 
    {
        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        double XMin = 0;    double XMax = 1; double RatioXMin = 0; double RatioXMax = 1;
        double YMin = 0.20; double YMax = 1; double RatioYMin = 0; double RatioYMax = 0.20;

        double XTitleSize = 0.04;  double XLabelSize = 0.036;  double XTitleOffset = 4.0;
        double YTitleSize = 0.04;  double YLabelSize = 0.036;  double YTitleOffset = 0.9;

        double PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin);

        //create the canvas & TLegend for the plot
        TCanvas *c = new TCanvas("c1", "c1", 1200, 1300);
        c->Divide(1,2);
        c->cd(1);
        gPad->Clear();
        gPad->SetPad(XMin, YMin, XMax, YMax);
        TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88); 
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42); 

        //get maximum from histos and fill TLegend
        double min  = 0.0;
        double max  = 0.0;
        double lmax = 0.0;

        //Create the THStack for the background 
        THStack *bgStack = new THStack();
        //Make seperate histogram from sum of BG histograms  
        TH1* hbgSum      = nullptr;

        //sort the bgEntries by number of events
        for(auto& h : bgEntries_)
        {
            h.histName = histName;
            h.rebin    = rebin;
            h.retrieveHistogram();
        }
        std::sort(bgEntries_.begin(), bgEntries_.end(), compareNEvents);

        for(int iBG = bgEntries_.size() - 1; iBG >= 0; --iBG) {
            bgStack->Add(bgEntries_[iBG].h.get(), bgEntries_[iBG].drawOptions.c_str());

            if(!hbgSum)
                hbgSum = static_cast<TH1*>(bgEntries_[iBG].h->Clone());
            else
            hbgSum->Add(bgEntries_[iBG].h.get());
        }

        //data / get new histogram from file
        data_.histName = histName;
        data_.rebin    = rebin;
        data_.retrieveHistogram();
        
        if(data_.drawHisto) // define for empty data plot  
            leg->AddEntry(data_.h.get(), data_.legEntry.c_str(), data_.drawOptions.c_str());
        smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, true);
    

        //background
        for(auto& entry : bgEntries_) {
            entry.setFillColor();

            //add histograms to TLegend
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "F");
        }
    
        smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);


        //signal / get new histogram
        //for(auto& entry : sigEntries_) {
        //    entry.histName = histName;
        //    entry.rebin    = rebin;
        //    entry.retrieveHistogram();

        //    //add histograms to TLegend
        //    leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L"); 
        //    smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        //}
    

        //Set Canvas margin 
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.01);


        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hbgSum->GetBinLowEdge(1), hbgSum->GetBinLowEdge(hbgSum->GetNbinsX()) + hbgSum->GetBinWidth(hbgSum->GetNbinsX())));
        dummy.setupAxes();
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());


        //Set the y-range of the histogram
        if(isLogY) {

            double locMin  = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin  = legSpan + log10(locMin);

            if(log10(lmax) > legMin) {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
        }

        dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
    
        } else {
        
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
                dummy.h->GetYaxis()->SetRangeUser(0.0, max*1.2);
        }


        //set x-axis range
        if(xmin < xmax) 
            dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        dummy.draw(); 

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //plot background stack
        bgStack->Draw("same");

        TH1F* dataMc = data_.clone("aclone");
        dataMc->Divide(data_.getHisto().get(), hbgSum);
        dataMc->SetMinimum(0.5);
        dataMc->SetMaximum(1.5);
        dataMc->GetYaxis()->SetNdivisions(305);
        dataMc->SetTitle("");
        dataMc->SetLineWidth(3);
        dataMc->SetMarkerSize(1.3);
        dataMc->SetMarkerStyle(20);
        dataMc->SetMarkerColor(dataMc->GetLineColor());
        dataMc->GetYaxis()->SetTitle("Data / MC");
        dataMc->GetXaxis()->SetTitle(xAxisLabel.c_str());
        dataMc->GetXaxis()->SetTitleSize(PadFactor*XTitleSize); dataMc->GetYaxis()->SetTitleSize(0.8*PadFactor*YTitleSize);
        dataMc->GetXaxis()->SetLabelSize(PadFactor*XLabelSize); dataMc->GetYaxis()->SetLabelSize(PadFactor*YLabelSize);
        dataMc->GetYaxis()->SetTitleOffset(1.4*YTitleOffset/PadFactor); dataMc->GetXaxis()->SetTitleOffset(XTitleOffset/PadFactor);

        //plot signal histograms
        //for(const auto& entry : sigEntries_) {
        //    entry.draw();
        //}
    
        //plot data histogram
        if(data_.drawHisto) // define for empty data plot
            data_.draw();

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");


        //Draw CMS & lumi lables & cut labels
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);
    
        TLatex mark;
        mark.SetNDC(true);
        
        mark.SetTextAlign(11);
        mark.SetTextSize(0.050);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");
        
        mark.SetTextSize(0.030); 
        mark.SetTextFont(52);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.11, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary"); 
        
        mark.SetTextFont(42); 
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);

        mark.SetTextAlign(11);
        mark.SetTextFont(42);
        mark.SetTextSize(0.030);
        mark.DrawLatex(0.51, 0.89, cutlabel.c_str()); 

        // Now time for the ratio
        c->cd(2);

        gPad->Clear();
        gPad->SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax);
        gPad->SetGridy();

        //Set Canvas margin 
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.01);
        gPad->SetBottomMargin(0.35);

        dataMc->Draw("E1 P");

        //Calculate simple significance      
        //const double totBG = hbgSum->Integral();
        //const double nSig = sigEntries_.at(0).nEvents;          
        //const double sig = nSig / sqrt( totBG + pow ( 0.2*totBG, 2) ) ;

        //TLatex significance;  
        //significance.SetNDC(true);
        //significance.SetTextAlign(11);
        //significance.SetTextFont(52);
        //significance.SetTextSize(0.030);
        //significance.DrawLatex(0.15, 0.89, ("Significance = "+std::to_string(sig)).c_str());

        //save new plot to file
        //c->Print((histName + ".png").c_str());
        c->Print(("plots/" + histName + ".pdf").c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
        delete bgStack;
        delete hbgSum;
    }
};


// -----------------------------------------------------------------------------
// Main FUNCTION 
/// -----------------------------------------------------------------------------

int main()
{
    //gROOT->SetBatch(1);
    gStyle->SetOptStat("");
    // entry for data
    // this uses the initializer syntax to initialize the histInfo object
    // 'leg entry'  'root file'     'draw options'  'draw color'
    histInfo data = {"Data", "condor/2016Data_hadd/2016_Data_SingleLepton.root", "PEX0", kBlack, true};

    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries = {

        {"T#bar{T}",        "condor/2016Bkgd_hadd/2016_TT.root",              "hist", kBlue - 7   },
        {"WJetsToLNu",      "condor/2016Bkgd_hadd/2016_WJetsToLNu.root",      "hist", kYellow + 1 },
        {"DYJetsToLL_M-50", "condor/2016Bkgd_hadd/2016_DYJetsToLL_M-50.root", "hist", kOrange + 2 },
        {"QCD",             "condor/2016Bkgd_hadd/2016_QCD.root",             "hist", kGreen + 1  },
        {"ST",              "condor/2016Bkgd_hadd/2016_ST.root",              "hist", kRed + 1    },
        {"Diboson",         "condor/2016Bkgd_hadd/2016_Diboson.root",         "hist", kMagenta + 1},
        {"TTX",             "condor/2016Bkgd_hadd/2016_TTX.root",             "hist", kCyan + 1   },
        {"Triboson",        "condor/2016Bkgd_hadd/2016_Triboson.root",        "hist", kGray       },
    };

    //vector summarizing signal histograms to include in the plot
    //std::vector<histInfo> sigEntries = { 
 
    //    {"RPV m_{#tildet} = 550", "condor/hadd_2016_MC_4thSlides_25.07.2019/2016_RPV_2t6j_mStop-550.root",        "hist", kOrange - 3}, 
    //    {"RPV m_{#tildet} = 350", "condor/hadd_2016_MC_4thSlides_25.07.2019/2016_RPV_2t6j_mStop-350.root",        "hist", kGreen + 3 },
    //    {"SYY m_{#tildet} = 900", "condor/hadd_2016_MC_4thSlides_25.07.2019/2016_StealthSYY_2t6j_mStop-900.root", "hist", kBlue + 1  },
    //};

    //make plotter object with the required sources for histograms specified
    //Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));
    Plotter plt(std::move(data), std::move(bgEntries));

    std::vector<std::string> cut {

        "1l_7j_ge1b", 
        "1l_ge2b" 
    };

    for (const auto& cutlabel : cut) {
        plt.plot( "h_ht_"+cutlabel, "H_{T} [GeV]",    "Events", true, cutlabel, 999.9, -999.9, 4); // for log scale 
        plt.plot( "h_lPt_"+cutlabel, "Lepton p_{T} [GeV]",    "Events", true, cutlabel, 999.9, -999.9, 4); // for log scale 
    }
}
