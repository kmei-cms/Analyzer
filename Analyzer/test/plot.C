#include "Analyzer/Analyzer/test/plotter.h"

class FisherHolder
{
public:
    std::vector<std::string> cutNames_;
    std::string plotName_;
};

int main()
{
    std::string path = "";
    //std::string path = "atleast2Tops";
    //std::string path = "exactly2Tops";

    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry root file                 draw options  draw color
    std::vector<histInfo> data = {
        //{"All BG", "allBG.root"            , "PEX0",       kBlack}
    };

    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries = {
        {"DYJetsToLL_M-50", "condor/output-files/" + path + "/DYJetsToLL_M-50/DYJetsToLL_M-50.root", "hist", kOrange + 2 },        
        {"Rare",            "condor/output-files/" + path + "/Rare/Rare.root",                       "hist", kCyan + 1   },
        {"Diboson",         "condor/output-files/" + path + "/Diboson/Diboson.root",                 "hist", kMagenta + 1},
        {"WJetsToLNu",      "condor/output-files/" + path + "/WJetsToLNu/WJetsToLNu.root",           "hist", kYellow + 1 },
        {"ST",              "condor/output-files/" + path + "/ST/ST.root",                           "hist", kRed + 1    },
        {"QCD",             "condor/output-files/" + path + "/QCD/QCD.root",                         "hist", kGreen + 1  },
        {"T#bar{T}",        "condor/output-files/" + path + "/TT/TT.root",                           "hist", kBlue - 7   },
    };

    //vector summarizing signal histograms to include in the plot
    std::vector<histInfo> sigEntries = {
        {"RPV 350", "condor/output-files/" + path + "/AllSignal/MyAnalysis_rpv_stop_350_0.root",         "hist", kMagenta + 2},
        {"SYY 650", "condor/output-files/" + path + "/AllSignal/MyAnalysis_stealth_stop_650_SYY_0.root", "hist", kGreen + 3  },
    };

    //make plotter object with the required sources for histograms specified
    Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));


    // --------------------
    // - Make stack plots
    // --------------------

    std::vector<std::string> mycuts_0l 
    {
        ""                      ,
        "ge6j"                   ,
        "HT500"                 ,
        "ge2b"                   ,
        "1t"                    ,
        "2t"                    ,
        "ge1t"                   ,
        "ge2t"                   ,
        "ge6j_HT500"             ,
        "ge6j_HT500_ge1b"         ,
        "ge6j_HT500_ge2b"         ,

        "ge6j_HT500_ge2b_1t"      ,
        "ge6j_HT500_ge2b_1t_f1"   , "ge6j_HT500_ge2b_1t_f2"   , "ge6j_HT500_ge2b_1t_f3"   , "ge6j_HT500_ge2b_1t_f4"   ,
        "ge6j_HT500_ge2b_ge1t"     ,
        "ge6j_HT500_ge2b_ge1t_f1"  , "ge6j_HT500_ge2b_ge1t_f2"  , "ge6j_HT500_ge2b_ge1t_f3"  , "ge6j_HT500_ge2b_ge1t_f4"  ,

        "ge6j_HT500_ge2b_1t1"     , "ge6j_HT500_ge2b_1t2"     , "ge6j_HT500_ge2b_1t3"     ,
        "ge6j_HT500_ge2b_1t1_f1"  , "ge6j_HT500_ge2b_1t1_f2"  , "ge6j_HT500_ge2b_1t1_f3"  , "ge6j_HT500_ge2b_1t1_f4"  ,
        "ge6j_HT500_ge2b_1t2_f1"  , "ge6j_HT500_ge2b_1t2_f2"  , "ge6j_HT500_ge2b_1t2_f3"  , "ge6j_HT500_ge2b_1t2_f4"  ,
        "ge6j_HT500_ge2b_1t3_f1"  , "ge6j_HT500_ge2b_1t3_f2"  , "ge6j_HT500_ge2b_1t3_f3"  , "ge6j_HT500_ge2b_1t3_f4"  ,
        "ge6j_HT500_ge2b_ge1t1"    , "ge6j_HT500_ge2b_ge1t2"    , "ge6j_HT500_ge2b_ge1t3"    ,
        "ge6j_HT500_ge2b_ge1t1_f1" , "ge6j_HT500_ge2b_ge1t1_f2" , "ge6j_HT500_ge2b_ge1t1_f3" , "ge6j_HT500_ge2b_ge1t1_f4" ,
        "ge6j_HT500_ge2b_ge1t2_f1" , "ge6j_HT500_ge2b_ge1t2_f2" , "ge6j_HT500_ge2b_ge1t2_f3" , "ge6j_HT500_ge2b_ge1t2_f4" ,
        "ge6j_HT500_ge2b_ge1t3_f1" , "ge6j_HT500_ge2b_ge1t3_f2" , "ge6j_HT500_ge2b_ge1t3_f3" , "ge6j_HT500_ge2b_ge1t3_f4" ,

        "ge6j_HT500_ge2b_2t"      ,
        "ge6j_HT500_ge2b_2t_f1"   , "ge6j_HT500_ge2b_2t_f2"   , "ge6j_HT500_ge2b_2t_f3"   , "ge6j_HT500_ge2b_2t_f4"   ,
        "ge6j_HT500_ge2b_ge2t"     ,
        "ge6j_HT500_ge2b_ge2t_f1"  , "ge6j_HT500_ge2b_ge2t_f2"  , "ge6j_HT500_ge2b_ge2t_f3"  , "ge6j_HT500_ge2b_ge2t_f4"  ,

        "ge6j_HT500_ge2b_2t11"    , "ge6j_HT500_ge2b_2t12"    , "ge6j_HT500_ge2b_2t13"    , "ge6j_HT500_ge2b_2t22"    , "ge6j_HT500_ge2b_2t23" , "ge6j_HT500_ge2b_2t33" ,
        "ge6j_HT500_ge2b_2t11_f1" , "ge6j_HT500_ge2b_2t11_f2" , "ge6j_HT500_ge2b_2t11_f3" , "ge6j_HT500_ge2b_2t11_f4" ,
        "ge6j_HT500_ge2b_2t12_f1" , "ge6j_HT500_ge2b_2t12_f2" , "ge6j_HT500_ge2b_2t12_f3" , "ge6j_HT500_ge2b_2t12_f4" ,
        "ge6j_HT500_ge2b_2t13_f1" , "ge6j_HT500_ge2b_2t13_f2" , "ge6j_HT500_ge2b_2t13_f3" , "ge6j_HT500_ge2b_2t13_f4" ,
        "ge6j_HT500_ge2b_2t22_f1" , "ge6j_HT500_ge2b_2t22_f2" , "ge6j_HT500_ge2b_2t22_f3" , "ge6j_HT500_ge2b_2t22_f4" ,
        "ge6j_HT500_ge2b_2t23_f1" , "ge6j_HT500_ge2b_2t23_f2" , "ge6j_HT500_ge2b_2t23_f3" , "ge6j_HT500_ge2b_2t23_f4" ,
        "ge6j_HT500_ge2b_2t33_f1" , "ge6j_HT500_ge2b_2t33_f2" , "ge6j_HT500_ge2b_2t33_f3" , "ge6j_HT500_ge2b_2t33_f4" ,
        "ge6j_HT500_ge2b_ge2t11"   , "ge6j_HT500_ge2b_ge2t12"   , "ge6j_HT500_ge2b_ge2t13"   , "ge6j_HT500_ge2b_ge2t22"   , "ge6j_HT500_ge2b_ge2t23", "ge6j_HT500_ge2b_ge2t33",
        "ge6j_HT500_ge2b_ge2t11_f1", "ge6j_HT500_ge2b_ge2t11_f2", "ge6j_HT500_ge2b_ge2t11_f3", "ge6j_HT500_ge2b_ge2t11_f4",
        "ge6j_HT500_ge2b_ge2t12_f1", "ge6j_HT500_ge2b_ge2t12_f2", "ge6j_HT500_ge2b_ge2t12_f3", "ge6j_HT500_ge2b_ge2t12_f4",
        "ge6j_HT500_ge2b_ge2t13_f1", "ge6j_HT500_ge2b_ge2t13_f2", "ge6j_HT500_ge2b_ge2t13_f3", "ge6j_HT500_ge2b_ge2t13_f4",
        "ge6j_HT500_ge2b_ge2t22_f1", "ge6j_HT500_ge2b_ge2t22_f2", "ge6j_HT500_ge2b_ge2t22_f3", "ge6j_HT500_ge2b_ge2t22_f4",
        "ge6j_HT500_ge2b_ge2t23_f1", "ge6j_HT500_ge2b_ge2t23_f2", "ge6j_HT500_ge2b_ge2t23_f3", "ge6j_HT500_ge2b_ge2t23_f4",
        "ge6j_HT500_ge2b_ge2t33_f1", "ge6j_HT500_ge2b_ge2t33_f2", "ge6j_HT500_ge2b_ge2t33_f3", "ge6j_HT500_ge2b_ge2t33_f4",
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

    // --------------------
    // - Make fisher plots
    // --------------------

    std::vector< FisherHolder > fisherHolder
    {
        { {"h_njets_0l_ge6j_HT500_ge2b_1t_f1"  , "h_njets_0l_ge6j_HT500_ge2b_1t_f2"  , "h_njets_0l_ge6j_HT500_ge2b_1t_f3"  , "h_njets_0l_ge6j_HT500_ge2b_1t_f4"  } , "njets_0l_ge6j_HT500_ge2b_1t"  },
        { {"h_njets_0l_ge6j_HT500_ge2b_1t1_f1" , "h_njets_0l_ge6j_HT500_ge2b_1t1_f2" , "h_njets_0l_ge6j_HT500_ge2b_1t1_f3" , "h_njets_0l_ge6j_HT500_ge2b_1t1_f4" } , "njets_0l_ge6j_HT500_ge2b_1t1" },
        { {"h_njets_0l_ge6j_HT500_ge2b_1t2_f1" , "h_njets_0l_ge6j_HT500_ge2b_1t2_f2" , "h_njets_0l_ge6j_HT500_ge2b_1t2_f3" , "h_njets_0l_ge6j_HT500_ge2b_1t2_f4" } , "njets_0l_ge6j_HT500_ge2b_1t2" },
        { {"h_njets_0l_ge6j_HT500_ge2b_1t3_f1" , "h_njets_0l_ge6j_HT500_ge2b_1t3_f2" , "h_njets_0l_ge6j_HT500_ge2b_1t3_f3" , "h_njets_0l_ge6j_HT500_ge2b_1t3_f4" } , "njets_0l_ge6j_HT500_ge2b_1t3" },
        { {"h_njets_0l_ge6j_HT500_ge2b_2t_f1"  , "h_njets_0l_ge6j_HT500_ge2b_2t_f2"  , "h_njets_0l_ge6j_HT500_ge2b_2t_f3"  , "h_njets_0l_ge6j_HT500_ge2b_2t_f4"  } , "njets_0l_ge6j_HT500_ge2b_2t"  },
        { {"h_njets_0l_ge6j_HT500_ge2b_2t11_f1", "h_njets_0l_ge6j_HT500_ge2b_2t11_f2", "h_njets_0l_ge6j_HT500_ge2b_2t11_f3", "h_njets_0l_ge6j_HT500_ge2b_2t11_f4"} , "njets_0l_ge6j_HT500_ge2b_2t11"},
        { {"h_njets_0l_ge6j_HT500_ge2b_2t12_f1", "h_njets_0l_ge6j_HT500_ge2b_2t12_f2", "h_njets_0l_ge6j_HT500_ge2b_2t12_f3", "h_njets_0l_ge6j_HT500_ge2b_2t12_f4"} , "njets_0l_ge6j_HT500_ge2b_2t12"},
        { {"h_njets_0l_ge6j_HT500_ge2b_2t13_f1", "h_njets_0l_ge6j_HT500_ge2b_2t13_f2", "h_njets_0l_ge6j_HT500_ge2b_2t13_f3", "h_njets_0l_ge6j_HT500_ge2b_2t13_f4"} , "njets_0l_ge6j_HT500_ge2b_2t13"},
        { {"h_njets_0l_ge6j_HT500_ge2b_2t22_f1", "h_njets_0l_ge6j_HT500_ge2b_2t22_f2", "h_njets_0l_ge6j_HT500_ge2b_2t22_f3", "h_njets_0l_ge6j_HT500_ge2b_2t22_f4"} , "njets_0l_ge6j_HT500_ge2b_2t22"},
        { {"h_njets_0l_ge6j_HT500_ge2b_2t23_f1", "h_njets_0l_ge6j_HT500_ge2b_2t23_f2", "h_njets_0l_ge6j_HT500_ge2b_2t23_f3", "h_njets_0l_ge6j_HT500_ge2b_2t23_f4"} , "njets_0l_ge6j_HT500_ge2b_2t23"},
        { {"h_njets_0l_ge6j_HT500_ge2b_2t33_f1", "h_njets_0l_ge6j_HT500_ge2b_2t33_f2", "h_njets_0l_ge6j_HT500_ge2b_2t33_f3", "h_njets_0l_ge6j_HT500_ge2b_2t33_f4"} , "njets_0l_ge6j_HT500_ge2b_2t33"},


        { {"h_njets_0l_ge6j_HT500_ge2b_ge1t_f1"  , "h_njets_0l_ge6j_HT500_ge2b_ge1t_f2"  , "h_njets_0l_ge6j_HT500_ge2b_ge1t_f3"  , "h_njets_0l_ge6j_HT500_ge2b_ge1t_f4"  } , "njets_0l_ge6j_HT500_ge2b_ge1t"  },
        { {"h_njets_0l_ge6j_HT500_ge2b_ge1t1_f1" , "h_njets_0l_ge6j_HT500_ge2b_ge1t1_f2" , "h_njets_0l_ge6j_HT500_ge2b_ge1t1_f3" , "h_njets_0l_ge6j_HT500_ge2b_ge1t1_f4" } , "njets_0l_ge6j_HT500_ge2b_ge1t1" },
        { {"h_njets_0l_ge6j_HT500_ge2b_ge1t2_f1" , "h_njets_0l_ge6j_HT500_ge2b_ge1t2_f2" , "h_njets_0l_ge6j_HT500_ge2b_ge1t2_f3" , "h_njets_0l_ge6j_HT500_ge2b_ge1t2_f4" } , "njets_0l_ge6j_HT500_ge2b_ge1t2" },
        { {"h_njets_0l_ge6j_HT500_ge2b_ge1t3_f1" , "h_njets_0l_ge6j_HT500_ge2b_ge1t3_f2" , "h_njets_0l_ge6j_HT500_ge2b_ge1t3_f3" , "h_njets_0l_ge6j_HT500_ge2b_ge1t3_f4" } , "njets_0l_ge6j_HT500_ge2b_ge1t3" },
        { {"h_njets_0l_ge6j_HT500_ge2b_ge2t_f1"  , "h_njets_0l_ge6j_HT500_ge2b_ge2t_f2"  , "h_njets_0l_ge6j_HT500_ge2b_ge2t_f3"  , "h_njets_0l_ge6j_HT500_ge2b_ge2t_f4"  } , "njets_0l_ge6j_HT500_ge2b_ge2t"  },
        { {"h_njets_0l_ge6j_HT500_ge2b_ge2t11_f1", "h_njets_0l_ge6j_HT500_ge2b_ge2t11_f2", "h_njets_0l_ge6j_HT500_ge2b_ge2t11_f3", "h_njets_0l_ge6j_HT500_ge2b_ge2t11_f4"} , "njets_0l_ge6j_HT500_ge2b_ge2t11"},
        { {"h_njets_0l_ge6j_HT500_ge2b_ge2t12_f1", "h_njets_0l_ge6j_HT500_ge2b_ge2t12_f2", "h_njets_0l_ge6j_HT500_ge2b_ge2t12_f3", "h_njets_0l_ge6j_HT500_ge2b_ge2t12_f4"} , "njets_0l_ge6j_HT500_ge2b_ge2t12"},
        { {"h_njets_0l_ge6j_HT500_ge2b_ge2t13_f1", "h_njets_0l_ge6j_HT500_ge2b_ge2t13_f2", "h_njets_0l_ge6j_HT500_ge2b_ge2t13_f3", "h_njets_0l_ge6j_HT500_ge2b_ge2t13_f4"} , "njets_0l_ge6j_HT500_ge2b_ge2t13"},
        { {"h_njets_0l_ge6j_HT500_ge2b_ge2t22_f1", "h_njets_0l_ge6j_HT500_ge2b_ge2t22_f2", "h_njets_0l_ge6j_HT500_ge2b_ge2t22_f3", "h_njets_0l_ge6j_HT500_ge2b_ge2t22_f4"} , "njets_0l_ge6j_HT500_ge2b_ge2t22"},
        { {"h_njets_0l_ge6j_HT500_ge2b_ge2t23_f1", "h_njets_0l_ge6j_HT500_ge2b_ge2t23_f2", "h_njets_0l_ge6j_HT500_ge2b_ge2t23_f3", "h_njets_0l_ge6j_HT500_ge2b_ge2t23_f4"} , "njets_0l_ge6j_HT500_ge2b_ge2t23"},
        { {"h_njets_0l_ge6j_HT500_ge2b_ge2t33_f1", "h_njets_0l_ge6j_HT500_ge2b_ge2t33_f2", "h_njets_0l_ge6j_HT500_ge2b_ge2t33_f3", "h_njets_0l_ge6j_HT500_ge2b_ge2t33_f4"} , "njets_0l_ge6j_HT500_ge2b_ge2t33"},

    };
    
    for (auto& f : fisherHolder)
    {
        plt.plotFisher(f.cutNames_,  f.plotName_, "N_{J}", "Events", true);
    }
}
