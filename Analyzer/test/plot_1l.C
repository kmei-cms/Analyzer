#include "Analyzer/Analyzer/test/plotter.h"
#include <iostream>

class FisherHolder
{
public:
    std::vector<std::string> cutNames_;
    std::string plotName_;
};

void setHistInfo(const std::string& path, std::vector<histInfo>& data, std::vector<histInfo>& bg, std::vector<histInfo>& sig, int color = 0)
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry root file                 draw options  draw color
    data = {
        {"Data_JetHT", "condor/output-files/" + path + "/Data_JetHT/Data_JetHT.root", "PEX0", kBlack},
    };
    
    std::vector<int> bgColor = {kRed, kBlack, kBlue};
    std::vector<int> sigColor = {kBlack, kBlue, kRed};

    bg = {
        {"T#bar{T}",        "condor/output-files/" + path + "/TT/TT.root",                           "hist", bgColor[color]        },
    };
    sig = {
        {"RPV 850", "condor/output-files/" + path + "/AllSignal/MyAnalysis_rpv_stop_850_0.root",         "hist", sigColor[color]        },
        //{"SYY 650", "condor/output-files/" + path + "/AllSignal/MyAnalysis_stealth_stop_650_SYY_0.root", "hist", kBlue      },
        {"RPV 350", "condor/output-files/" + path + "/AllSignal/MyAnalysis_rpv_stop_350_0.root",         "hist", kGreen       },
    };
}

int main()
{
    TH1::AddDirectory(false);

    std::string pathGRfalse = "deepESM_GRfalse_1Layer_8Vars";
    //std::string pathGRfalse = "deepESM_GRfalse_1Layer_9Vars_nJets";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_8Vars";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_9Vars_nJets";
    //std::string pathGRtrue = "deepESM_GRtrue15-2_4Vars_1Layer";
    std::string pathGRtrue = "deepESM_GRtrue15-20_4Vars_1Layer";
    std::string pathFisher = "deepESM_GRfalse_1Layer";

    std::vector<histInfo> data_GRfalse, bg_GRfalse, sig_GRfalse;
    std::vector<histInfo> data_GRtrue, bg_GRtrue, sig_GRtrue;
    std::vector<histInfo> data_fisher, bg_fisher, sig_fisher;

    setHistInfo(pathGRfalse, data_GRfalse, bg_GRfalse, sig_GRfalse, 0);
    setHistInfo(pathGRtrue, data_GRtrue, bg_GRtrue, sig_GRtrue, 1);
    setHistInfo(pathFisher, data_fisher, bg_fisher, sig_fisher, 2);

    //make histInfoCollection
    HistInfoCollection histInfoCollection_GRfalse(data_GRfalse, bg_GRfalse, sig_GRfalse);
    HistInfoCollection histInfoCollection_GRtrue(data_GRtrue, bg_GRtrue, sig_GRtrue);
    HistInfoCollection histInfoCollection_fisher(data_fisher, bg_fisher, sig_fisher);

    // vector of histInfoCollection for Roc Curves
    std::map< std::string, HistInfoCollection > rocMap = { {"NN", histInfoCollection_GRfalse},
                                                           {"NN GR", histInfoCollection_GRtrue},
                                                           {"fisher", histInfoCollection_fisher},
    };

    //make plotter object with the required sources for histograms specified
    //Plotter pltRoc( std::move(rocMapTest) );
    Plotter pltRocCompare( std::move(rocMap) );
    Plotter plt( std::move(histInfoCollection_GRtrue) );
    //Plotter plt( std::move(histInfoCollection_GRfalse) );
    //Plotter plt( std::move(histInfoCollection_fisher) );

    // --------------------
    // - Make stack plots
    // --------------------

    std::vector<std::string> mycuts_1l 
    {
        //"ge6j"                        
        //"ge2b"                        
        //"1t"                          
        //"2t"                          
        //"ge1t"                        
        //"ge2t"                        
        //"ge6j_ge2b"                   
        "ge6j_ge1b"                                            
        //"ge6j_ge1b_1t"                
        //"ge6j_ge1b_ge1t"                                                    
        //"ge6j_ge1b_1t1"               
        //"ge6j_ge1b_1t2"               
        //"ge6j_ge1b_1t3"               
        //"ge6j_ge1b_1t2or3"            
        //"ge6j_ge1b_ge1t1"             
        //"ge6j_ge1b_ge1t2"             
        //"ge6j_ge1b_ge1t3"             
   };

    for(std::string mycut : mycuts_1l)
    {
        //plt.plotStack( "h_njets_1l_"+mycut, "N_{J}" , "Events", true);
        //plt.plotStack( "h_ntops_1l_"+mycut, "N_{T}" , "Events", true);
        //plt.plotStack( "h_nb_1l_"   +mycut, "N_{B}" , "Events", true);        
        //plt.plotStack( "h_fisher_1l_"+mycut, "fisher value" , "Events", true, 4);        
        //
        //// Make Normalized fisher
        //pltSkim.plotNormFisher("h_fisher_1l_"+mycut, "fisher value" , "Events", false, 4);
        //
        // - Make  Roc Curve
        //pltRoc.plotRocFisher("h_deepESM_1l_"+mycut,"Background","Signal", false);
        pltRocCompare.plotRocFisher("h_deepESM_1l_"+mycut,"Background","Signal", true);
        //        
        ////Need these until we un blind
        //plt.plotStack( "blind_njets_1l_"+mycut, "N_{J}" , "Events", true);
        //plt.plotStack( "blind_ntops_1l_"+mycut, "N_{T}" , "Events", true);
        //plt.plotStack( "blind_nb_1l_"   +mycut, "N_{B}" , "Events", true);        
        //plt.plotStack( "blind_fisher_1l_"+mycut, "fisher value" , "Events", true, 4);        
    }
    
    //plt.plotStack("h_met"     , "MET"   , "Events", true);
    //plt.plotStack("h_bdt"     , "bdt"   , "Events", true);
    //plt.plotStack("h_fisher"  , "fisher", "Events", true, 4);
    //plt.plotStack("h_njets"   , "N_{J}" , "Events", true);
    //plt.plotStack("h_nb"      , "N_{B}" , "Events", true);
    //plt.plotStack("h_ntops"   , "N_{T}" , "Events", true);
    //plt.plotStack("h_ntops_j1", "N_{1T}", "Events", true);
    //plt.plotStack("h_ntops_j2", "N_{2T}", "Events", true);
    //plt.plotStack("h_ntops_j3", "N_{3T}", "Events", true);

    //plt.plotStack("h_met"     , "MET"   , "Events", true);
    
    //plt.plotStack("blind_met"     , "MET"   , "Events", true);
    //plt.plotStack("blind_ht"      , "H_{T}" , "Events", true);
    //plt.plotStack("blind_bdt"     , "bdt"   , "Events", true);
    //plt.plotStack("blind_fisher"  , "fisher", "Events", true);
    //plt.plotStack("blind_njets"   , "N_{J}" , "Events", true);
    //plt.plotStack("blind_nb"      , "N_{B}" , "Events", true);
    //plt.plotStack("blind_ntops"   , "N_{T}" , "Events", true);
    //plt.plotStack("blind_ntops_j1", "N_{1T}", "Events", true);
    //plt.plotStack("blind_ntops_j2", "N_{2T}", "Events", true);
    //plt.plotStack("blind_ntops_j3", "N_{3T}", "Events", true);

    // --------------------
    // - Make fisher plots
    // --------------------

    std::vector< FisherHolder > fisherHolder
    {
        { {"h_njets_1l_ge6j_ge1b_d1"  , "h_njets_1l_ge6j_ge1b_d2"  , "h_njets_1l_ge6j_ge1b_d3"  , "h_njets_1l_ge6j_ge1b_d4"  } , "njets_1l_ge6j_ge1b"  },
        //{ {"h_njets_1l_ge6j_ge1b_f1"  , "h_njets_1l_ge6j_ge1b_f2"  , "h_njets_1l_ge6j_ge1b_f3"  , "h_njets_1l_ge6j_ge1b_f4"  } , "njets_1l_ge6j_ge1b"  },
    };
    
    for (auto& f : fisherHolder)
    {
        plt.plotFisher(f.cutNames_,  f.plotName_, "N_{J}", "Events", true, 9, "DeepESM");
        plt.plotRatioFisher(f.cutNames_,  f.plotName_, "N_{J}", "N_{J+1} / N_{J}", false, 6, "DeepESM");
        //plt.plotFisher(f.cutNames_,  f.plotName_, "N_{J}", "Events", true, 9, "Fisher");
        //plt.plotRatioFisher(f.cutNames_,  f.plotName_, "N_{J}", "N_{J+1} / N_{J}", false, 6, "Fisher");
    }
    
    // --------------------
    // - Compute Yields
    // --------------------

    //for(std::string mycut : mycuts_1l)
    //{
    //    const auto& yieldMap = histInfoCollection.computeYields("h_njets_1l_"+mycut,"njets",12,20);
    //    //const auto& yieldMap = histInfoCollection.computeYields("h_njets_1l_"+mycut,"njets",6,8);
    //    auto allbg  = yieldMap.find("AllBG"); auto data_JetHT = yieldMap.find("Data_JetHT"); auto ttbar = yieldMap.find("T#bar{T}"); auto QCD = yieldMap.find("QCD");
    //    printf("%45.45s:  %s: %12.1lf   %s: %12.1lf   %s: %12.1lf   %s: %12.1lf\n",("h_njets_1l_"+mycut).c_str(), (data_JetHT->first).c_str(), data_JetHT->second, (allbg->first).c_str(), allbg->second, "ttbar", ttbar->second, (QCD->first).c_str(), QCD->second );
    //}       
}
