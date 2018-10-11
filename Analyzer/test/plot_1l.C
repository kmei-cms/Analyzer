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
        //{"Data_JetHT", "condor/output-files/" + path + "/Data_JetHT/Data_JetHT.root", "PEX0", kBlack},
        {"Data_SingleLepton" , "condor/output-files/" + path + "/Data_SingleLepton/Data_SingleLepton.root", "PEX0", kBlack},
        //{"Data 1 #gamma", "condor/output-files/" + path + "/Data_SinglePhoton/Data_SinglePhoton.root", "PEX0", kBlack}, 
    };
    
    std::vector<int> bgColor = {kRed, kBlack, kBlue};
    std::vector<int> sigColor = {kBlack, kBlue, kRed};

    bg = {
        {"DYJetsToLL_M-50", "condor/output-files/" + path + "/DYJetsToLL_M-50/DYJetsToLL_M-50.root", "hist", kOrange + 2 },        
        {"Rare",            "condor/output-files/" + path + "/Rare/Rare.root",                       "hist", kCyan + 1   },
        {"Diboson",         "condor/output-files/" + path + "/Diboson/Diboson.root",                 "hist", kMagenta + 1},
        {"WJetsToLNu",      "condor/output-files/" + path + "/WJetsToLNu/WJetsToLNu.root",           "hist", kYellow + 1 },
        {"ST",              "condor/output-files/" + path + "/ST/ST.root",                           "hist", kRed + 1    },
        {"QCD",             "condor/output-files/" + path + "/QCD/QCD.root",                         "hist", kGreen + 1  },
        {"T#bar{T}",        "condor/output-files/" + path + "/TT/TT.root",                           "hist", kBlue - 7   },
    };
    //bg = {
    //    //{"T#bar{T}",   "condor/output-files/" + path + "/TT/TT.root", "hist", bgColor[color]        },
    //    {"Rare",           "condor/output-files/" + path + "/Rare/Rare.root",             "hist", kOrange + 2 },  
    //    {"Diboson",        "condor/output-files/" + path + "/Diboson/Diboson.root",       "hist", kCyan + 1   },
    //    {"ST",             "condor/output-files/" + path + "/ST/ST.root",                 "hist", kMagenta + 1},  
    //    {"T#bar{T}",       "condor/output-files/" + path + "/TT/TT.root",                 "hist", kYellow + 1 },
    //    {"T#bar{T}#gamma", "condor/output-files/" + path + "/TTGJets/TTGJets.root",       "hist", kBlue + 1   },
    //    {"W(l#nu) + jets", "condor/output-files/" + path + "/WJetsToLNu/WJetsToLNu.root", "hist", kRed + 1    },
    //    {"QCD",            "condor/output-files/" + path + "/QCD/QCD.root",               "hist", kBlue - 7   },  
    //    {"#gamma + jets",  "condor/output-files/" + path + "/SP/SP.root",                 "hist", kGreen + 1  },  
    //};
    sig = {
        //{"RPV 850", "condor/output-files/" + path + "/AllSignal/MyAnalysis_rpv_stop_850_0.root",         "hist", sigColor[color]        },
        {"RPV 850", "condor/output-files/" + path + "/AllSignal/MyAnalysis_rpv_stop_850_0.root",         "hist", kRed        },
        {"SYY 650", "condor/output-files/" + path + "/AllSignal/MyAnalysis_stealth_stop_650_SYY_0.root", "hist", kMagenta    },
        //{"RPV 450", "condor/output-files/" + path + "/AllSignal/MyAnalysis_rpv_stop_450_0.root",         "hist", kRed       },
        {"RPV 350", "condor/output-files/" + path + "/AllSignal/MyAnalysis_rpv_stop_350_0.root",         "hist", kCyan       },
    };
}

int main()
{
    TH1::AddDirectory(false);

    //std::string pathGRfalse = "deepESM_GRfalse_1Layer_8Vars";
    std::string pathGRfalse = "deepESM_GRfalse_1Layer_24Vars";
    //std::string pathGRfalse = "deepESM_GRfalse_1Layer_9Vars_nJets";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_8Vars";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_9Vars_nJets";
    //std::string pathGRtrue = "deepESM_GRtrue15-2_4Vars_1Layer";
    //std::string pathGRtrue = "deepESM_GRtrue15-20_4Vars_1Layer";
    //std::string pathGRtrue = "deepESM_7Vars_GRtrue_HighWeight";
    //std::string pathGRtrue = "deepESM_7Vars_GRtrue_HighWeight_noMETMassVar";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_24Vars";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_24Vars_v2";
    std::string pathGRtrue = "deepESM_MyCodeTrue_3Layer_28Vars";
    std::string pathFisher = "oldTest/DeepESMTests_JoesCode/deepESM_v1";
    //std::string pathPhoton = "photonCR";
    std::string pathPhoton = "photonCR_Barrel";

    std::vector<histInfo> data_GRfalse, bg_GRfalse, sig_GRfalse;
    std::vector<histInfo> data_GRtrue, bg_GRtrue, sig_GRtrue;
    std::vector<histInfo> data_fisher, bg_fisher, sig_fisher;
    std::vector<histInfo> data_photon, bg_photon, sig_photon;

    setHistInfo(pathGRfalse, data_GRfalse, bg_GRfalse, sig_GRfalse, 0);
    setHistInfo(pathGRtrue, data_GRtrue, bg_GRtrue, sig_GRtrue, 1);
    setHistInfo(pathFisher, data_fisher, bg_fisher, sig_fisher, 2);
    setHistInfo(pathPhoton, data_photon, bg_photon, sig_photon, 0);

    //make histInfoCollection
    HistInfoCollection histInfoCollection_GRfalse(data_GRfalse, bg_GRfalse, sig_GRfalse);
    HistInfoCollection histInfoCollection_GRtrue(data_GRtrue, bg_GRtrue, sig_GRtrue);
    HistInfoCollection histInfoCollection_fisher(data_fisher, bg_fisher, sig_fisher);
    HistInfoCollection histInfoCollection_photon(data_photon, bg_photon, sig_photon);

    // vector of histInfoCollection for Roc Curves
    std::map< std::string, HistInfoCollection > rocMap = { //{"NN", histInfoCollection_GRfalse},
                                                           {"NN GR", histInfoCollection_GRtrue},
                                                           {"fisher", histInfoCollection_fisher},
    };

    //make plotter object with the required sources for histograms specified
    //Plotter pltRoc( std::move(rocMapTest) );
    Plotter pltRocCompare( std::move(rocMap) );
    Plotter plt( std::move(histInfoCollection_GRtrue) );
    //Plotter plt( std::move(histInfoCollection_GRfalse) );
    //Plotter plt( std::move(histInfoCollection_fisher) );
    //Plotter plt( std::move(histInfoCollection_photon) );

    // --------------------
    // - Make stack plots
    // --------------------

    std::vector<std::string> mycuts_1l = {
        //"1l_ge6j",                        
        //"1l_ge2b",                       
        //"1l_1t",                          
        //"1l_2t",                          
        //"1l_ge1t",                        
        //"1l_ge2t",                        
        "1l_ge6j_ge1b",
        //"1l_ge6j_ge2b",                 
        //"1l_ge6j_ge1b_1t",
        //"1l_ge6j_ge1b_ge1t",                                                    
        //"1l_ge6j_ge1b_1t1",               
        //"1l_ge6j_ge1b_1t2",               
        //"1l_ge6j_ge1b_1t3",               
        //"1l_ge6j_ge1b_1t2or3",            
        //"1l_ge6j_ge1b_ge1t1",             
        //"1l_ge6j_ge1b_ge1t2",             
        //"1l_ge6j_ge1b_ge1t3",             
        //"1l_ge6j_ge1b_ge5-6esm", 
        //"1l_ge6j_ge1b_ge6-7esm", 
        //"1l_ge6j_ge1b_ge7-8esm", 
        //"1l_ge6j_ge1b_ge8-95esm",
        //"1l_ge6j_ge1b_ge95esm",
        //"1l_ge6j_ge1b_ge8esm",
        //"1l_ge6j_ge1b_1t_ge95esm",
        //"1l_ge6j_ge1b_1t_ge8esm",
        //"0l",
        //"0l_1g",
        //"0l_ge7j_1g",
    };

    for(std::string mycut : mycuts_1l)
    {
        plt.plotStack( "h_njets_"+mycut, "N_{J}" , "Events", true);
        plt.plotStack( "h_deepESM_"+mycut, "DeepESM" , "Events", true, 10);
        //plt.plotStack( "h_ntops_"+mycut, "N_{T}" , "Events", true);
        //plt.plotStack( "h_nb_"   +mycut, "N_{B}" , "Events", true);        
        //plt.plotStack( "h_fisher_"+mycut, "fisher value" , "Events", true, 4);        
        //plt.plotStack( "h_photonPt_"+mycut, "P_{T}^{#gamma}" , "Events", false, 10);
        //
        //// Make Normalized fisher
        //pltSkim.plotNormFisher("h_fisher_1l_"+mycut, "fisher value" , "Events", false, 4);
        //plt.plotNormFisher("h_BestComboMass_1l_"+mycut, "Average BestCombo Mass [GeV]" , "Events", false, 4);
        //plt.plotNormFisher("h_BestComboPt_1l_"+mycut, "Average BestCombo P_{T} [GeV]" , "Events", false, 4);
        //plt.plotNormFisher("h_BestComboMassDiff_1l_"+mycut, "BestCombo Mass Diff [GeV]" , "Events", false, 2);
        //plt.plotNormFisher("h_BestComboMassDiffAbs_1l_"+mycut, "BestCombo Mass Abs(Diff) [GeV]" , "Events", false, 2);
        //plt.plotNormFisher("h_BestComboRelDiff_1l_"+mycut, "BestCombo Rel Diff" , "Events", false, 2);
        //plt.plotNormFisher("h_BestComboRelDiffAbs_1l_"+mycut, "BestCombo Rel Abs(Diff)" , "Events", false, 2);
        //
        // - Make  Roc Curve
        //pltRoc.plotRocFisher("h_deepESM_1l_"+mycut,"Background","Signal", false);
        //pltRocCompare.plotRocFisher(mycut,"Background","Signal", true);
        //        
        ////Need these until we un blind
        plt.plotStack( "blind_njets_"+mycut, "N_{J}" , "Events", true);
        plt.plotStack( "blind_deepESM_"+mycut, "DeepESM" , "Events", true, 10);
        //plt.plotStack( "blind_ntops_"+mycut, "N_{T}" , "Events", true);
        //plt.plotStack( "blind_nb_"   +mycut, "N_{B}" , "Events", true);        
        //plt.plotStack( "blind_fisher_"+mycut, "fisher value" , "Events", true, 10);        
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

    //std::vector< FisherHolder > fisherHolder
    //{
    //    //{ {"h_njets_1l_ge6j_ge1b_d1"  , "h_njets_1l_ge6j_ge1b_d2"  , "h_njets_1l_ge6j_ge1b_d3"  , "h_njets_1l_ge6j_ge1b_d4"  } , "njets_1l_ge6j_ge1b"  },
    //    //{ {"h_njets_1l_ge6j_ge1b_ge6-7esm", "h_njets_1l_ge6j_ge1b_ge7-8esm", "h_njets_1l_ge6j_ge1b_ge8-95esm", "h_njets_1l_ge6j_ge1b_ge95esm" }, "njets_1l_ge6j_ge1b_esmBins"},
    //    //{ {"h_njets_1l_ge6j_ge1b_ge5-6esm", "h_njets_1l_ge6j_ge1b_ge6-7esm", "h_njets_1l_ge6j_ge1b_ge7-8esm", "h_njets_1l_ge6j_ge1b_ge8-95esm", "h_njets_1l_ge6j_ge1b_ge95esm" }, "njets_1l_ge6j_ge1b_esmBins"},
    //    { {"h_njets_1l_ge6j_ge1b_f1"  , "h_njets_1l_ge6j_ge1b_f2"  , "h_njets_1l_ge6j_ge1b_f3"  , "h_njets_1l_ge6j_ge1b_f4"  } , "njets_1l_ge6j_ge1b"  },
    //};
    //
    //for (auto& f : fisherHolder)
    //{
    //    //plt.plotFisher(f.cutNames_,  f.plotName_, "N_{J}", "Events", true, 9, "DeepESM");
    //    //plt.plotRatioFisher(f.cutNames_,  f.plotName_, "N_{J}", "N_{J+1} / N_{J}", false, 6, "DeepESM");
    //    plt.plotFisher(f.cutNames_,  f.plotName_, "N_{J}", "Events", true, 9, "Fisher");
    //    plt.plotRatioFisher(f.cutNames_,  f.plotName_, "N_{J}", "N_{J+1} / N_{J}", false, 6, "Fisher");
    //}
    
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
