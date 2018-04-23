#ifndef PLOTTER_H
#define PLOTTER_H

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
#include <iostream>
#include <map>

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
    void setupAxes(double xOffset, double yOffset, double xTitle, double yTitle, double xLabel, double yLabel)
    {
        h->SetStats(0);
        h->SetTitle(0);
        h->GetXaxis()->SetTitleOffset(xOffset);
        h->GetYaxis()->SetTitleOffset(yOffset);
        h->GetXaxis()->SetTitleSize(xTitle);
        h->GetYaxis()->SetTitleSize(yTitle);
        h->GetXaxis()->SetLabelSize(xLabel);
        h->GetYaxis()->SetLabelSize(yLabel);
        if(h->GetXaxis()->GetNdivisions() % 100 > 5) h->GetXaxis()->SetNdivisions(6, 5, 0);
    }

    //helper function for pads
    void setupPad(double left, double right, double top, double bottom)
    {
        gPad->SetLeftMargin(left);
        gPad->SetRightMargin(right);
        gPad->SetTopMargin(top);
        gPad->SetBottomMargin(bottom);
        gPad->SetTicks(1,1);
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

    histInfo(TH1* h) : legEntry(h->GetName()), histFile(""), histName(h->GetName()), drawOptions(""), color(kWhite), rebin(0), h(h)
    {
    }

    ~histInfo()
    {
    }
};

class HistInfoCollection
{
private:
    void makeYieldMap(std::map<std::string, double>& yieldMap, const std::vector<std::vector<histInfo>*>& dataSets, const std::shared_ptr<TH1>& hbgSum, const std::string& histType, const int min, const int max)
    {
        int index = -1;
        double yield = 0;
        for(const auto* set : dataSets)
        {
            for(const auto& entry : *set)
            {
                if( entry.histName.find(histType) != std::string::npos )
                {
                    index++;
                    if(index==0)
                    {
                        yield = 0;
                        for(int i = min; i<=max; i++)
                        {
                            yield+=hbgSum->GetBinContent(i + 1);
                            //std::cout<<hbgSum->GetBinContent(i)<<std::endl;
                        }
                        yieldMap.insert ( std::pair<std::string,double>( "AllBG", yield ) );
                    }
                    yield = 0;
                    for(int i = min; i<=max; i++)
                    {
                        yield+=entry.h->GetBinContent(i + 1);
                    }
                    yieldMap.insert ( std::pair<std::string,double>( entry.legEntry, yield ) );
                }
            }
        }

        //for(const auto& entry : yieldMap)
        //{
        //    std::cout<<entry.first<<"  "<<entry.second<<std::endl;
        //}
    }

public:
    std::vector<histInfo> dataVec_, bgVec_, sigVec_;

    HistInfoCollection(std::vector<histInfo>&& dataVec, std::vector<histInfo>&& bgVec, std::vector<histInfo>&& sigVec) : dataVec_(dataVec), bgVec_(bgVec), sigVec_(sigVec) {}
    HistInfoCollection(std::vector<histInfo>& dataVec, std::vector<histInfo>& bgVec, std::vector<histInfo>& sigVec) : dataVec_(dataVec), bgVec_(bgVec), sigVec_(sigVec) {}

    void setUpBG(const std::string& histName, int rebin, THStack* bgStack, std::shared_ptr<TH1>& hbgSum, const bool& setFill = true)
    {
        bool firstPass = true;
        for(auto& entry : bgVec_)
        {
            entry.histName = histName;
            entry.rebin = rebin;
            entry.retrieveHistogram();
    
            bgStack->Add(entry.h.get(), entry.drawOptions.c_str());
            if(firstPass) 
            {
                hbgSum.reset( static_cast<TH1*>(entry.h->Clone()) );
                firstPass = false;
            }
            else 
            {
                hbgSum->Add(entry.h.get());
            }
            
            if(setFill)
            {
                entry.setFillColor();
            }
        }
    }

    void setUpSignal(const std::string& histName, int rebin)
    {
        for(auto& entry : sigVec_)
        {
            entry.histName = histName;
            entry.rebin = rebin;
            entry.retrieveHistogram();
            entry.setLineStyle(2);
    
        }
    }

    void setUpData(const std::string& histName, int rebin)
    {
        for(auto& entry : dataVec_)
        {
            entry.histName = histName;
            entry.rebin = rebin;
            entry.retrieveHistogram();
        }
    }

    std::map<std::string,double> computeYields(const std::string& histName, const std::string& histType, const int min, const int max, const int rebin = -1)
    {
        std::map<std::string,double> yieldMap;
        std::vector<std::vector<histInfo>*> dataSets  = { &dataVec_, &bgVec_, &sigVec_ };
        std::shared_ptr<TH1> hbgSum(nullptr);
        THStack* bgStack = new THStack();

        setUpBG(histName, rebin, bgStack, hbgSum, false);
        setUpSignal(histName, rebin);
        setUpData(histName, rebin);
        
        makeYieldMap(yieldMap, dataSets, hbgSum, histType, min, max);

        delete bgStack;

        return yieldMap;
    }

    std::map<std::string,double> computeYields(const std::string& histName, const std::string& histType, const std::shared_ptr<TH1>& hbgSum, const int min, const int max)
    {
        std::map<std::string,double> yieldMap;
        std::vector<std::vector<histInfo>*> dataSets  = { &dataVec_, &bgVec_, &sigVec_ };
        
        makeYieldMap(yieldMap, dataSets, hbgSum, histType, min, max);

        return yieldMap;
    }

};

class Plotter
{
private:
    //Collection of histInfo
    HistInfoCollection hc_;
    std::shared_ptr<TH1> hbgSum_;

public:
    Plotter(HistInfoCollection&& hc) : hc_(hc), hbgSum_(nullptr) {}

    //This is a helper function which will keep the plot from overlapping with the legend
    void smartMax(const TH1* const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax, const bool error)
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

    void plotStack(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
    {
        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        //switch to the canvas to ensure it is the active object
        c->cd();

        // Upper plot will be in pad1: TPad(x1, y1, x2, y2)
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        //pad1->SetBottomMargin(0); // Upper and lower plot are joined
        //pad1->SetGridy();         // Horizontal grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad

        //Create TLegend
        TLegend *leg = new TLegend(0.20, 0.76, 0.89, 0.88);
        //TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(3);
        leg->SetTextFont(42);

        // ------------------------
        // -  Setup plots
        // ------------------------
        THStack *bgStack = new THStack();
        hc_.setUpBG(histName, rebin, bgStack, hbgSum_);
        hc_.setUpSignal(histName, rebin);
        hc_.setUpData(histName, rebin);

        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hbgSum_->GetBinLowEdge(1), hbgSum_->GetBinLowEdge(hbgSum_->GetNbinsX()) + hbgSum_->GetBinWidth(hbgSum_->GetNbinsX())));

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        // -----------------------
        // -  Plot Background
        // -----------------------
        for(auto& entry : hc_.bgVec_)
        {
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "F");
        }
        smartMax(hbgSum_.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        bgStack->Draw("same");

        // -------------------------
        // -   Plot Signal
        // -------------------------
        for(const auto& entry : hc_.sigVec_)
        {
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L");
            smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
            entry.draw();
        }

        // ------------------------
        // -  Plot Data
        // ------------------------
        for(const auto& entry : hc_.dataVec_)
        {
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), entry.drawOptions.c_str());
            smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, true);
            entry.draw();
        }

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        setupDummy(dummy, leg, histName, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, min, max, lmax);
        dummy.setupPad(0.12, 0.06, 0.08, 0.0);
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        //drawLables(lumi);

        //Compute and draw yields for njets min to max
        drawYields(histName,"njets",12,20);

        // lower plot will be in pad2
        c->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1, 0.3);
        //pad2->SetTopMargin(0);
        //pad2->SetBottomMargin(0.2);
        pad2->SetGridy(); // Horizontal grid
        pad2->Draw();
        pad2->cd();       // pad2 becomes the current pad        

        //make ratio dummy
        histInfo ratioDummy(new TH1D("rdummy", "rdummy", 1000, hc_.dataVec_[0].h->GetBinLowEdge(1), hc_.dataVec_[0].h->GetBinLowEdge(hc_.dataVec_[0].h->GetNbinsX()) + hc_.dataVec_[0].h->GetBinWidth(hc_.dataVec_[0].h->GetNbinsX())));
        ratioDummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        //ratioDummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        ratioDummy.h->GetYaxis()->SetTitle("Data / BG");
        ratioDummy.h->GetXaxis()->SetTickLength(0.1);
        ratioDummy.h->GetYaxis()->SetTickLength(0.045);
        ratioDummy.setupAxes(1.2, 0.4, 0.15, 0.15, 0.13, 0.13);
        ratioDummy.h->GetYaxis()->SetNdivisions(6, 5, 0);
        ratioDummy.h->GetXaxis()->SetRangeUser(xmin, xmax);
        ratioDummy.h->GetYaxis()->SetRangeUser(0.5, 1.5);
        ratioDummy.h->SetStats(0);
        //ratioDummy.h->SetMinimum(0.5);
        //ratioDummy.h->SetMaximum(1.5);

        //Make ratio histogram for data / background.
        histInfo ratio((TH1*)hc_.dataVec_[0].h->Clone());

        // set pad margins: setupPad(left, right, top, bottom)
        ratio.setupPad(0.12, 0.06, 0.0, 0.40);
        
        ratio.drawOptions = "ep";
        ratio.color = kBlack;

        //ratio.h->SetLineColor(kBlack);
        //ratio.h->Sumw2();
        //ratio.h->SetStats(0);
        ratio.h->Divide(hbgSum_.get());
        ratio.h->SetMarkerStyle(21);

        ratioDummy.draw();
        ratio.draw("same");

        //save new plot to file
        c->Print(("outputPlots/" + histName + ".pdf").c_str());
        c->Print(("outputPlots/" + histName + ".png").c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
        delete bgStack;
    }

    void plotNormFisher(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
    {
        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        //switch to the canvas to ensure it is the active object
        c->cd();

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.12);
        gPad->SetTicks(1,1);

        //Create TLegend
        TLegend *leg = new TLegend(0.20, 0.76, 0.89, 0.88);
        //TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(3);
        leg->SetTextFont(42);

        // ------------------------
        // -  Setup plots
        // ------------------------
        THStack *bgStack = new THStack();
        hc_.setUpBG(histName, rebin, bgStack, hbgSum_, false);
        hc_.setUpSignal(histName, rebin);
        hc_.setUpData(histName, rebin);

        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hbgSum_->GetBinLowEdge(1), hbgSum_->GetBinLowEdge(hbgSum_->GetNbinsX()) + hbgSum_->GetBinWidth(hbgSum_->GetNbinsX())));

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        // -----------------------
        // -  Plot Background
        // -----------------------
        for(auto& entry : hc_.bgVec_)
        {
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L");
            entry.h->Scale( 1.0/entry.h->Integral( 0, entry.h->GetNbinsX() + 1 ) );
            smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);            
            entry.draw();
        }
        hbgSum_->SetLineColor(kBlack);
        hbgSum_->SetMarkerColor(kBlack);
        //std::cout<<hbgSum_->GetNbinsX() + 1<<std::endl;
        leg->AddEntry(hbgSum_.get(), "AllBG", "L");
        hbgSum_->Scale( 1.0/hbgSum_->Integral( 0, hbgSum_->GetNbinsX() + 1 ) );
        smartMax(hbgSum_.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        hbgSum_->Draw("hist same");

        // -------------------------
        // -   Plot Signal
        // -------------------------
        for(const auto& entry : hc_.sigVec_)
        {
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L");
            entry.h->Scale( 1.0/entry.h->Integral( 0, entry.h->GetNbinsX() + 1 ) );
            smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
            entry.draw();
        }

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        setupDummy(dummy, leg, histName, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, min, max, lmax);
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        //drawLables(lumi);

        //save new plot to file
        c->Print(("outputPlots/fisherNorm_" + histName + ".pdf").c_str());
        c->Print(("outputPlots/fisherNorm_" + histName + ".png").c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
        delete bgStack;
    }

    void plotFisher(const std::vector<std::string>& histNameVec, const std::string& histTitle, const std::string& xAxisLabel, 
                    const std::string& yAxisLabel = "Events",    const bool isLogY = false,    const int fixedJetBin = -1,  
                    const double xmin = 999.9, const double xmax = -999.9, int rebin = -1,  double lumi = 36100)
    {
        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        //switch to the canvas to ensure it is the active object
        c->cd();

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.12);
        gPad->SetTicks(1,1);

        //Create TLegend
        TLegend *leg = new TLegend(0.20, 0.76, 0.89, 0.88);
        //TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(3);
        leg->SetTextFont(42);
        
        //Setup color and name for fisher plots
        std::vector<int> color = {kRed,kBlue,kGreen+2,kMagenta};
        std::vector<std::string> fisherNames = {"Fisher Bin 1","Fisher Bin 2","Fisher Bin 3","Fisher Bin 4"};

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        // -----------------------
        // -  Background
        // -----------------------
        
        std::vector<TH1*> hbgSumVec;
        int index = -1;
        double scale = 1;
        int fixedBin;
        for(auto& histName : histNameVec)
        {
            index++;
            TH1* hbgSum = nullptr;
            for(auto& entry : hc_.bgVec_)
            {
                //Get new histogram
                entry.histName = histName;
                entry.rebin = rebin;
                entry.retrieveHistogram();
            
                if(!hbgSum) hbgSum = static_cast<TH1*>(entry.h->Clone());
                else        hbgSum->Add(entry.h.get());
            
            }
            //std::cout<<"Index = "<<index<<"  histName = "<<histName<<std::endl;
            hbgSum->SetLineColor(color[index]);
            hbgSum->SetMarkerColor(color[index]);
            leg->AddEntry(hbgSum,(fisherNames[index]).c_str(), "l");

            if(fixedJetBin == -1 && index == 0) 
            {
                fixedBin = hbgSum->GetMaximumBin();
            }
            else if (index == 0)
            {
                fixedBin = fixedJetBin + 1;
            }

            double numEvents = hbgSum->GetBinContent( fixedBin );
            if(index == 0)
            {
                scale = numEvents;
            }
            else
            {
                hbgSum->Scale(scale/numEvents);
            }            

            smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);
            hbgSumVec.push_back(hbgSum);
        }
        
        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hbgSumVec[0]->GetBinLowEdge(1), hbgSumVec[0]->GetBinLowEdge(hbgSumVec[0]->GetNbinsX()) + hbgSumVec[0]->GetBinWidth(hbgSumVec[0]->GetNbinsX())));
        setupDummy(dummy, leg, histTitle, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, min, max, lmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //draw all fisher plots
        for(auto& h : hbgSumVec)
        {
            h->Draw("same E");
        }

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        //drawLables(lumi);

        //save new plot to file
        c->Print( ("outputPlots/fisher_" + histTitle + ".pdf").c_str() );
        c->Print( ("outputPlots/fisher_" + histTitle + ".png").c_str() );

        //clean up dynamic memory
        delete c;
        delete leg;
        for (auto* obj :  hbgSumVec)
        {
            delete obj;
        }
    }

    void setupDummy(histInfo dummy, TLegend *leg, const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel, const bool isLogY, 
                    const double xmin, const double xmax, double min, double max, double lmax)
    {
        dummy.setupAxes(1.2, 1.1, 0.045, 0.045, 0.045, 0.045);
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
    }

    void drawLables(double lumi)
    {
        //Draw CMS and lumi lables
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);

        TLatex mark;
        mark.SetNDC(true);

        //Draw CMS mark
        mark.SetTextAlign(11);
        mark.SetTextSize(0.050);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");
        mark.SetTextSize(0.040);
        mark.SetTextFont(52);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.11, 1 - (gPad->GetTopMargin() - 0.017), "Stealth 2018");

        //Draw lumistamp
        mark.SetTextFont(42);
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);        
    }    

    void drawYields(std::string histName, std::string histType, int min, int max)
    {
        const auto& yieldMap = hc_.computeYields(histName, histType, hbgSum_,  min, max);
        
        if(hc_.bgVec_[0].histName.find( histType ) != std::string::npos)
        {
            auto data   = yieldMap.find("Data_JetHT");
            auto allbg  = yieldMap.find("AllBG");
            auto qcd    = yieldMap.find("QCD");
            auto ttbar  = yieldMap.find("T#bar{T}");
            auto rpv350 = yieldMap.find("RPV 350");
            auto syy650 = yieldMap.find("SYY 650");
    
            char stamp_data[128];
            char stamp_allbg[128];
            char stamp_qcd[128];
            char stamp_ttbar[128];
            char stamp_rpv350[128];
            char stamp_syy650[128];
            sprintf(stamp_data,  "%-15.15s %8.0lf" , (data->first).c_str()  , data->second  );
            sprintf(stamp_allbg, "%-15.15s %8.0lf" , (allbg->first).c_str() , allbg->second );
            sprintf(stamp_qcd,   "%-15.15s %8.0lf" , (qcd->first  ).c_str() , qcd->second   );
            sprintf(stamp_ttbar, "%-15.15s %8.0lf" , "TTbar"                , ttbar->second );
            sprintf(stamp_rpv350,"%-15.15s %8.0lf" , (rpv350->first).c_str(), rpv350->second);
            sprintf(stamp_syy650,"%-15.15s %8.0lf" , (syy650->first).c_str(), syy650->second);
    
            TLatex mark;
            mark.SetNDC(true);
    
            //Draw CMS mark
            mark.SetTextAlign(11);
            mark.SetTextSize(0.030);
            mark.SetTextFont(61);
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18       ), "Njets>=12 Yield");
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.03), stamp_data       );
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.06), stamp_allbg      );
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.09), stamp_qcd        );
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.12), stamp_ttbar      );
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.15), stamp_rpv350     );
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.18), stamp_syy650     );
        }
    
    }

};

#endif
