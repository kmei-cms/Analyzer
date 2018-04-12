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

class HistInfoCollection
{
public:
    std::vector<histInfo> dataVec_, bgVec_, sigVec_;

    HistInfoCollection(std::vector<histInfo>&& dataVec, std::vector<histInfo>&& bgVec, std::vector<histInfo>&& sigVec) : dataVec_(dataVec), bgVec_(bgVec), sigVec_(sigVec)
    {
    }

    void setUpBG(const std::string& histName, int rebin, THStack* bgStack, std::shared_ptr<TH1>& hbgSum)
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
    
            entry.setFillColor();
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

    std::map<std::string,double> computeYields(const std::string& histName, const std::shared_ptr<TH1>& hbgSum, const int min, const int max)
    {
        std::map<std::string,double> yieldMap;
        std::vector<std::vector<histInfo>> dataSets  = {dataVec_,bgVec_,sigVec_};
        int index = -1;
        double yield = 0;

        for(const auto& set : dataSets)
        {
            for(const auto& entry : set)
            {
                if( entry.histName.find(histName) != std::string::npos )
                {
                    index++;
                    if(index==0)
                    {
                        yield = 0;
                        for(int i = min; i<=max; i++)
                        {
                            yield+=hbgSum->GetBinContent(i);
                            //std::cout<<hbgSum->GetBinContent(i)<<std::endl;
                        }
                        yieldMap.insert ( std::pair<std::string,double>( "AllBG", yield ) );
                    }
                    yield = 0;
                    for(int i = min; i<=max; i++)
                    {
                        yield+=entry.h.get()->GetBinContent(i);
                    }
                    yieldMap.insert ( std::pair<std::string,double>( entry.legEntry, yield ) );
                }
            }
        }

        //for(const auto& entry : yieldMap)
        //{
        //    std::cout<<entry.first<<"  "<<entry.second<<std::endl;
        //}
        
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
    Plotter(HistInfoCollection&& hc) : hc_(hc), hbgSum_(nullptr)
    {
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

    void plot(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
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
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        //drawLables(lumi);

        //Compute and draw yields for njets min to max
        drawYields("njets",12,20);

        //save new plot to file
        c->Print(("outputPlots/" + histName + ".pdf").c_str());
        c->Print(("outputPlots/" + histName + ".png").c_str());

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

    void drawYields(std::string histName, int min, int max)
    {
        const auto& yieldMap = hc_.computeYields(histName, hbgSum_, min, max);
        
        if(hc_.bgVec_[0].histName.find(histName) != std::string::npos)
        {
            auto allbg  = yieldMap.find("AllBG");
            auto qcd    = yieldMap.find("QCD");
            auto ttbar  = yieldMap.find("T#bar{T}");
            auto rpv350 = yieldMap.find("RPV 350");
            auto syy650 = yieldMap.find("SYY 650");
    
            char stamp_allbg[128];
            char stamp_qcd[128];
            char stamp_ttbar[128];
            char stamp_rpv350[128];
            char stamp_syy650[128];
            sprintf(stamp_allbg, "%-15s %i" , (allbg->first).c_str() , int(allbg->second ));
            sprintf(stamp_qcd,   "%-15s %i" , (qcd->first  ).c_str() , int(qcd->second   ));
            sprintf(stamp_ttbar, "%-15s %i" , "TTbar"                , int(ttbar->second ));
            sprintf(stamp_rpv350,"%-15s %i" , (rpv350->first).c_str(), int(rpv350->second));
            sprintf(stamp_syy650,"%-15s %i" , (syy650->first).c_str(), int(syy650->second));
    
            TLatex mark;
            mark.SetNDC(true);
    
            //Draw CMS mark
            mark.SetTextAlign(11);
            mark.SetTextSize(0.030);
            mark.SetTextFont(61);
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18       ), "Njets>=12 Yield");
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.03), stamp_allbg      );
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.06), stamp_qcd        );
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.09), stamp_ttbar      );
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.12), stamp_rpv350     );
            mark.DrawLatex(1 - (gPad->GetLeftMargin() + 0.25), 1 - (gPad->GetTopMargin() + 0.18 + 0.15), stamp_syy650     );
        }
    
    }

};

#endif
