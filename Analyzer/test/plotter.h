#ifndef PLOTTER_H
#define PLOTTER_H

#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TF1.h"
#include "HistInfoCollection.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <map>

class rocInfo
{
public:
    std::vector<double> rocVec;
    std::string legEntry;
    int color;
    bool firstOnly;
};

class Plotter
{
private:
    //Collection of histInfo
    HistInfoCollection hc_;
    std::shared_ptr<TH1> hbgSum_;
    //std::vector<HistInfoCollection> chc_;
    std::map< std::string, HistInfoCollection > mhc_;
    std::vector<std::shared_ptr<TH1>> hbgSumVec_;

public:
    Plotter(HistInfoCollection&& hc) : hc_(hc) {}
    //Plotter(std::vector<HistInfoCollection>&& chc) : chc_(chc) {}
    Plotter(std::map< std::string, HistInfoCollection >&& mhc) : mhc_(mhc) {}

    void plotStack(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, int rebin = -1, const double xmin = 999.9, const double xmax = -999.9, double lumi = 36100)
    {
        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        //switch to the canvas to ensure it is the active object
        c->cd();

        // Upper plot will be in pad1: TPad(x1, y1, x2, y2)
        TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.31, 1.0, 1.0);
        pad1->SetBottomMargin(0); // Upper and lower plot are joined
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
        hc_.setUpBG(histName, rebin, bgStack, hbgSum_, true, false);
        hc_.setUpSignal(histName, rebin);
        hc_.setUpData(histName, rebin);

        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hbgSum_->GetBinLowEdge(1), hbgSum_->GetBinLowEdge(hbgSum_->GetNbinsX()) + hbgSum_->GetBinWidth(hbgSum_->GetNbinsX())));

        //draw dummy axes
        dummy.setupAxes(1.2, 0.4, 0.15, 0.15, 0.13, 0.13);
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);
        gPad->SetTicks(1,1);

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

        //Draw CMS and lumi lables
        //drawLables(lumi);
            
        //Compute and draw yields for njets min to max
        //drawYields(histName,"njets",12,20);
        //drawYields(histName, "h", 0, hc_.dataVec_[0].h->GetNbinsX());

        //Draw dummy hist again to get axes on top of histograms
        setupDummy(dummy, leg, histName, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, min, max, lmax);
        //dummy.setupPad(0.12, 0.06, 0.08, 0.0);
        dummy.draw("AXIS");
                        
        // lower plot will be in pad2
        c->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        //pad2->SetGridy(); // Horizontal grid
        pad2->Draw();
        pad2->cd();       // pad2 becomes the current pad        
        gPad->SetTicks(1,1);
        
        //make ratio dummy
        histInfo ratioDummy(new TH1D("rdummy", "rdummy", 1000, hc_.dataVec_[0].h->GetBinLowEdge(1), 
                                     hc_.dataVec_[0].h->GetBinLowEdge(hc_.dataVec_[0].h->GetNbinsX()) + hc_.dataVec_[0].h->GetBinWidth(hc_.dataVec_[0].h->GetNbinsX())));
        ratioDummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        //ratioDummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        ratioDummy.h->GetYaxis()->SetTitle("Data / BG");
        ratioDummy.h->GetXaxis()->SetTickLength(0.1);
        ratioDummy.h->GetYaxis()->SetTickLength(0.045);
        ratioDummy.setupAxes(1.2, 0.5, 0.1, 0.1, 0.1, 0.1);
        ratioDummy.h->GetYaxis()->SetNdivisions(6, 5, 0);
        ratioDummy.h->GetXaxis()->SetRangeUser(xmin, xmax);
        ratioDummy.h->GetYaxis()->SetRangeUser(0, 2.1);
        ratioDummy.h->SetStats(0);
        //ratioDummy.h->SetMinimum(0.5);
        //ratioDummy.h->SetMaximum(1.5);

        //Make ratio histogram for data / background.
        histInfo ratio((TH1*)hc_.dataVec_[0].h->Clone());
            
        // set pad margins: setupPad(left, right, top, bottom)
        //ratio.setupPad(0.12, 0.06, 0.0, 0.40);
        
        ratio.drawOptions = "ep";
        ratio.color = kBlack;
        
        //ratio.h->SetLineColor(kBlack);
        //ratio.h->Sumw2();
        //ratio.h->SetStats(0);
        ratio.h->Divide(hbgSum_.get());
        ratio.h->SetMarkerStyle(21);
        
        ratioDummy.draw();
        ratio.draw("same");

        TF1* line = new TF1("1" ,"1" ,-2000,20000);
        line->SetLineColor(kRed);
        line->Draw("same");
        ratio.draw("same");

        //save new plot to file
        c->Print(("outputPlots/" + histName + ".pdf").c_str());
        //c->Print(("outputPlots/" + histName + ".png").c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
        delete bgStack;
        delete line;
    }

    void plotNormFisher(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, int rebin = -1, const double xmin = 999.9, const double xmax = -999.9, double lumi = 36100)
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
            entry.h->Scale( 1.0/entry.h->Integral() );
            smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);            
            entry.draw();
        }
        hbgSum_->SetLineColor(kBlack);
        hbgSum_->SetMarkerColor(kBlack);
        //std::cout<<hbgSum_->GetNbinsX() + 1<<std::endl;
        leg->AddEntry(hbgSum_.get(), "AllBG", "L");
        hbgSum_->Scale( 1.0/hbgSum_->Integral() );
        smartMax(hbgSum_.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        hbgSum_->Draw("hist same");

        // -------------------------
        // -   Plot Signal
        // -------------------------
        for(const auto& entry : hc_.sigVec_)
        {
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L");
            entry.h->Scale( 1.0/entry.h->Integral() );
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
        //c->Print(("outputPlots/fisherNorm_" + histName + ".png").c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
        delete bgStack;
    }

    void plotRocFisher(std::string histCut, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool firstOnly = false, int rebin = -1, const double xmin = 999.9, const double xmax = -999.9, double lumi = 36100)
    {
        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        c->cd();

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.12);
        gPad->SetTicks(1,1);

        //Create TLegend
        TLegend *leg = new TLegend(0.13, 0.75, 0.88, 0.9);
        //TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(3);
        leg->SetTextFont(42);
        leg->SetTextSize(0.02);

        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, 0, 1));
        dummy.draw();

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 1.0;
        double lmax = 1.0;

        // --------------------------
        // -  Make Roc Info and Plot
        // --------------------------
        std::vector<rocInfo> bgSumRocInfoVec;
        std::vector<TGraph*> graphVec;
        std::string histName;
        for(auto& mhc : mhc_)
        {
            if(mhc.first == "fisher")
                histName = "h_"+mhc.first+"_1l_"+histCut;
            else
                histName = "h_deepESM_1l_"+histCut;
            std::cout<<histName<<std::endl;
            THStack *bgStack = new THStack();
            std::shared_ptr<TH1> hbgSum;
            mhc.second.setUpBG(histName, rebin, bgStack, hbgSum, false);
            delete bgStack;
            mhc.second.setUpSignal(histName, rebin);
            rocInfo bgSumRocInfo = { makeFisherVec(hbgSum), "AllBG", mhc.second.bgVec_[0].color };
            std::vector<rocInfo> rocBgVec  = makeRocVec(mhc.second.bgVec_);
            std::vector<rocInfo> rocSigVec = makeRocVec(mhc.second.sigVec_);
            if(firstOnly) rocBgVec.emplace(rocBgVec.begin(), bgSumRocInfo);
            int lineStyle = (mhc.first == "Test") ?  kSolid : kDashed;
            int markStyle = (mhc.first == "Test") ?  kFullCircle : kFullSquare;
            drawRocCurve(mhc.first, graphVec, rocBgVec, rocSigVec, firstOnly, leg, lineStyle, markStyle);
            //std::cout<<histName<<" "<<mhc.first<<std::endl;
        }

        TF1* line1 = new TF1( "line1","1",0,1);
        line1->SetLineColor(kBlack);
        line1->Draw("same");
        TF1* line2 = new TF1( "line2","x",0,1);
        line2->SetLineColor(kBlack);
        line2->SetLineStyle(kDotted);
        line2->Draw("same");
        
        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        setupDummy(dummy, leg, histName, xAxisLabel, yAxisLabel, false, xmin, xmax, min, max, lmax);
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        //drawLables(lumi);

        //save new plot to file
        if(firstOnly) 
        {
            c->Print(("outputPlots/fisherRocCompare_" + histName + ".pdf").c_str());
            //c->Print(("outputPlots/fisherRocCompare_" + histName + ".png").c_str());
        }
        else
        {
            c->Print(("outputPlots/fisherRoc_" + histName + ".pdf").c_str());
            //c->Print(("outputPlots/fisherRoc_" + histName + ".png").c_str());
        }
        //c->Print("test.pdf");

        //clean up dynamic memory
        delete c;
        delete leg;
        for(auto* g : graphVec) delete g;
    }

    void plotFisher(const std::vector<std::string>& histNameVec, const std::string& histTitle, const std::string& xAxisLabel, 
                    const std::string& yAxisLabel = "Events",    const bool isLogY = false,    const int fixedJetBin = -1, const std::string name = "Fisher",
                    int rebin = -1, const double xmin = 999.9, const double xmax = -999.9, double lumi = 36100)
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
        std::vector<std::string> fisherNames = {name+" Bin 1",name+" Bin 2",name+" Bin 3",name+" Bin 4", name+" Bin 5"};

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
        c->Print( ("outputPlots/"+name+"_" + histTitle + ".pdf").c_str() );
        //c->Print( ("outputPlots/"+name+"_" + histTitle + ".png").c_str() );

        //clean up dynamic memory
        delete c;
        delete leg;
        for (auto* obj :  hbgSumVec)
        {
            delete obj;
        }
    }

    void plotRatioFisher(const std::vector<std::string>& histNameVec, const std::string& histTitle, const std::string& xAxisLabel, 
                         const std::string& yAxisLabel = "Events",    const bool isLogY = false,  const int start = 6, const std::string name = "Fisher",
                         int rebin = -1, const double xmin = 999.9, const double xmax = -999.9, double lumi = 36100)
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
        std::vector<std::string> fisherNames = {name+" Bin 1",name+" Bin 2",name+" Bin 3",name+" Bin 4", name+" Bin 5"};

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        // -----------------------
        // -  Background
        // -----------------------
        
        std::vector<TH1*> hbgRatioSumVec;
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

            TH1* hbgRatioSum = static_cast<TH1*>(new TH1F(hbgSum->GetName(),hbgSum->GetName(),hbgSum->GetNbinsX(),0, hbgSum->GetNbinsX()));
            hbgRatioSum->SetLineColor(color[index]);
            hbgRatioSum->SetLineWidth(2.0);
            hbgRatioSum->SetMarkerColor(color[index]);
            leg->AddEntry(hbgRatioSum,(fisherNames[index]).c_str(), "l");
            for(int i = 0; i < hbgSum->GetNbinsX(); i++)
            {
                auto getErrorRatioSq = [](int i, TH1* h){return pow(h->GetBinError(i)/h->GetBinContent(i),2);};
                double val = hbgSum->GetBinContent(i+1) / hbgSum->GetBinContent(i);
                double error = val*sqrt( getErrorRatioSq(i+1, hbgSum) + getErrorRatioSq(i, hbgSum) );
                if(hbgSum->GetBinContent(i) == 0 || hbgSum->GetBinContent(i+1) == 0 || i <= (start+1)) 
                {
                    hbgRatioSum->SetBinContent(i, 0);
                }
                else
                {
                    hbgRatioSum->SetBinContent(i, hbgSum->GetBinContent(i+1)/hbgSum->GetBinContent(i));
                    hbgRatioSum->SetBinError(i, error);
                }
            }

            smartMax(hbgRatioSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);
            hbgRatioSumVec.push_back(hbgRatioSum);
        }
        
        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hbgRatioSumVec[0]->GetBinLowEdge(1), hbgRatioSumVec[0]->GetBinLowEdge(hbgRatioSumVec[0]->GetNbinsX()) + hbgRatioSumVec[0]->GetBinWidth(hbgRatioSumVec[0]->GetNbinsX())));
        setupDummy(dummy, leg, histTitle, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, min, max, lmax);

        //draw dummy axes
        dummy.draw();
        dummy.h->SetMaximum(0.5);
        dummy.h->SetMinimum(0.0);

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //draw all fisher plots
        for(auto& h : hbgRatioSumVec)
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
        c->Print( ("outputPlots/"+name+"Ratio_" + histTitle + ".pdf").c_str() );
        //c->Print( ("outputPlots/"+name+"Ratio_" + histTitle + ".png").c_str() );

        //clean up dynamic memory
        delete c;
        delete leg;
        for (auto* obj :  hbgRatioSumVec)
        {
            delete obj;
        }
    }

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

    std::vector<double> makeFisherVec(std::shared_ptr<TH1> h)
    {
        h->Scale( 1.0 / h->Integral() );
        std::vector<double> v;
        for(int ii = 0; ii <= h->GetNbinsX(); ii++)
        {
            double val = h->Integral( ii, h->GetNbinsX());
            v.push_back( val );
        }
        return v;
    }

    std::vector<rocInfo> makeRocVec(const std::vector<histInfo>& vec)
    {
        std::vector<rocInfo> rocVec;
        for(const auto& entry : vec)
        {
            std::vector<double> v = makeFisherVec(entry.h);
            rocVec.push_back( {v, entry.legEntry, entry.color} );
        }
        return rocVec;    
}

    void drawRocCurve(const std::string& fType, std::vector<TGraph*>& graphVec, const std::vector<rocInfo>& rocBgVec, const std::vector<rocInfo>& rocSigVec, const bool firstOnly, TLegend* leg, int lineStyle, int markStyle)
    {
        int index = 0;
        for(const auto& mBg : rocBgVec)
        {
            index++;
            if(index > 2) break;
            for(const auto& mSig : rocSigVec)
            {
                int n = mBg.rocVec.size();
                double x[n], y[n];
                for(int i = 0; i < n; i++)
                {
                    x[i] = mBg.rocVec[i];
                    y[i] = mSig.rocVec[i];
                    //std::cout<<mBg.legEntry<<" "<<x[i]<<" "<<mSig.legEntry<<" "<<y[i]<<std::endl;
                }
                TGraph* g = new TGraph (n, x, y);
                g->SetLineWidth(2);
                g->SetLineStyle(lineStyle);
                g->SetLineColor( mBg.color );
                g->SetMarkerSize(0.7);
                g->SetMarkerStyle(markStyle);
                g->SetMarkerColor( mSig.color );
                g->Draw("same LP");                
                leg->AddEntry(g, (fType + " " + mBg.legEntry + " vs " + mSig.legEntry).c_str(), "LP");
                graphVec.push_back(g);
                //if(firstOnly) break; 
            }
            if(firstOnly) break; 
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
        //mark.SetTextSize(0.040);
        mark.SetTextFont(52);
        //mark.DrawLatex(gPad->GetLeftMargin() + 0.1, 1 - (gPad->GetTopMargin() - 0.017), "Stealth 2018");
        mark.DrawLatex(gPad->GetLeftMargin() + 0.08, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary");

        //Draw lumistamp
        mark.SetTextFont(42);
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);        
    }    

    void drawYields(std::string histName, std::string histType, int min, int max)
    {
        const auto& yieldMap = hc_.computeYields(histName, histType, hbgSum_,  min, max);
        
        int index = 0;
        for(const auto& pair : yieldMap)
        {
            char stamp[128];
            sprintf(stamp, "%-15.15s %8.0lf", pair.first.c_str(), pair.second);
            TLatex mark;
            mark.SetNDC(true);
            
            //Draw CMS mark
            mark.SetTextAlign(11);
            mark.SetTextSize(0.030);
            mark.SetTextFont(61);
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.03*index), stamp);
            index++;
        }

        //if(hc_.bgVec_[0].histName.find( histType ) != std::string::npos)
        //{
        //    auto data   = yieldMap.find("Data_JetHT");
        //    auto allbg  = yieldMap.find("AllBG");
        //    auto qcd    = yieldMap.find("QCD");
        //    auto ttbar  = yieldMap.find("T#bar{T}");
        //    auto rpv350 = yieldMap.find("RPV 350");
        //    auto syy650 = yieldMap.find("SYY 650");
        //
        //    char stamp_data[128];
        //    char stamp_allbg[128];
        //    char stamp_qcd[128];
        //    char stamp_ttbar[128];
        //    char stamp_rpv350[128];
        //    char stamp_syy650[128];
        //    sprintf(stamp_data,  "%-15.15s %8.0lf" , (data->first).c_str()  , data->second  );
        //    sprintf(stamp_allbg, "%-15.15s %8.0lf" , (allbg->first).c_str() , allbg->second );
        //    sprintf(stamp_qcd,   "%-15.15s %8.0lf" , (qcd->first  ).c_str() , qcd->second   );
        //    sprintf(stamp_ttbar, "%-15.15s %8.0lf" , "TTbar"                , ttbar->second );
        //    sprintf(stamp_rpv350,"%-15.15s %8.0lf" , (rpv350->first).c_str(), rpv350->second);
        //    sprintf(stamp_syy650,"%-15.15s %8.0lf" , (syy650->first).c_str(), syy650->second);
        //
        //    TLatex mark;
        //    mark.SetNDC(true);
        //
        //    //Draw CMS mark
        //    mark.SetTextAlign(11);
        //    mark.SetTextSize(0.030);
        //    mark.SetTextFont(61);
        //    mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18       ), "Njets>=12 Yield");
        //    mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.03), stamp_data       );
        //    mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.06), stamp_allbg      );
        //    mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.09), stamp_qcd        );
        //    mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.12), stamp_ttbar      );
        //    mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.15), stamp_rpv350     );
        //    mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.18), stamp_syy650     );
        //}    
    }
};

#endif
