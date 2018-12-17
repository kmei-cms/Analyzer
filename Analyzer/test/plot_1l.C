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
        {"Data_SingleLepton" , path + "/Data.root", "PEX0", kBlack},
        //{"Data 1 #gamma", "condor/output-files/" + path + "/Data_SinglePhoton/Data_SinglePhoton.root", "PEX0", kBlack}, 
    };
    
    std::vector<int> bgColor = {kRed, kBlack, kBlue};
    std::vector<int> sigColor = {kBlack, kBlue, kRed};

    bg = {
        {"DYJetsToLL_M-50", path + "/DYJetsToLL_M-50.root", "hist", kOrange + 2 },        
        {"Rare",            path + "/Rare.root",            "hist", kCyan + 1   },
        {"Diboson",         path + "/Diboson.root",         "hist", kMagenta + 1},
        {"WJetsToLNu",      path + "/WJetsToLNu.root",      "hist", kYellow + 1 },
        {"ST",              path + "/ST.root",              "hist", kRed + 1    },
        {"QCD",             path + "/QCD.root",             "hist", kGreen + 1  },
        {"T#bar{T}",        path + "/TT.root",              "hist", kBlue       },
        //{"T#bar{T}",        path + "/TTJets.root",              "hist", kBlue       },
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
        //{"RPV 850", path + "/AllSignal/MyAnalysis_rpv_stop_850_0.root",         "hist", sigColor[color]        },
        {"RPV 850", path + "/rpv_stop_850.root",         "hist", kRed        },
        {"SYY 650", path + "/stealth_stop_650_SYY.root", "hist", kMagenta    },
        //{"RPV 450", path + "/rpv_stop_450.root",         "hist", kRed       },
        {"RPV 350", path + "/rpv_stop_350.root",         "hist", kCyan       },
    };
}

int main()
{
    TH1::AddDirectory(false);

    //std::string pathGRfalse = "deepESM_GRfalse_1Layer_8Vars";
    //std::string pathGRfalse = "deepESM_GRfalse_1Layer_24Vars";
    //std::string pathGRfalse = "deepESM_GRfalse_1Layer_9Vars_nJets";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_8Vars";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_9Vars_nJets";
    //std::string pathGRtrue = "deepESM_GRtrue15-2_4Vars_1Layer";
    //std::string pathGRtrue = "deepESM_GRtrue15-20_4Vars_1Layer";
    //std::string pathGRtrue = "deepESM_7Vars_GRtrue_HighWeight";
    //std::string pathGRtrue = "deepESM_7Vars_GRtrue_HighWeight_noMETMassVar";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_24Vars";
    //std::string pathGRtrue = "deepESM_GRtrue_1Layer_24Vars_v2";

    //std::string pathGRtrue = "deepESM_MyCodeTrue_3Layer_28Vars";
    //std::string pathGRtrue = "deepESM_perNJet_1Layer";
    //std::string pathGRtrue = "Training_V2";
    std::string pathGRtrue = "condor/Analyze1Lep_Kerasv1.2.0";

    //std::string pathFisher = "oldTest/DeepESMTests_JoesCode/deepESM_v1";
    //std::string pathPhoton = "photonCR";
    //std::string pathPhoton = "photonCR_Barrel";

    std::vector<histInfo> data_GRfalse, bg_GRfalse, sig_GRfalse;
    std::vector<histInfo> data_GRtrue, bg_GRtrue, sig_GRtrue;
    std::vector<histInfo> data_fisher, bg_fisher, sig_fisher;
    std::vector<histInfo> data_photon, bg_photon, sig_photon;

    //setHistInfo(pathGRfalse, data_GRfalse, bg_GRfalse, sig_GRfalse, 0);
    setHistInfo(pathGRtrue, data_GRtrue, bg_GRtrue, sig_GRtrue, 1);
    //setHistInfo(pathFisher, data_fisher, bg_fisher, sig_fisher, 2);
    //setHistInfo(pathPhoton, data_photon, bg_photon, sig_photon, 0);

    //make histInfoCollection
    //HistInfoCollection histInfoCollection_GRfalse(data_GRfalse, bg_GRfalse, sig_GRfalse);
    HistInfoCollection histInfoCollection_GRtrue(data_GRtrue, bg_GRtrue, sig_GRtrue);
    //HistInfoCollection histInfoCollection_fisher(data_fisher, bg_fisher, sig_fisher);
    //HistInfoCollection histInfoCollection_photon(data_photon, bg_photon, sig_photon);

    // vector of histInfoCollection for Roc Curves
    std::map< std::string, HistInfoCollection > rocMap = { //{"NN", histInfoCollection_GRfalse},
                                                           {"NN GR", histInfoCollection_GRtrue},
                                                           //{"fisher", histInfoCollection_fisher},
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
        "",
        "_1l",
        "_1l_ge7j",                        
        "_1l_ge2b",                       
        //"_1l_1t",                          
        //"_1l_2t",                          
        //"_1l_ge1t",                        
        //"_1l_ge2t",                        
        "_1l_5to6j_ge1b",
        //"_1l_7j_ge1b", 
        //"_1l_8j_ge1b", 
        //"_1l_9j_ge1b", 
        //"_1l_10j_ge1b", 
        //"_1l_11j_ge1b", 
        //"_1l_12j_ge1b", 
        //"_1l_13j_ge1b", 
        //"_1l_14j_ge1b", 
        //"_1l_15j_ge1b",
        "_1l_ge7j_ge1b_noMbl",
        "_1l_ge7j_ge1b",
        "_1e_ge7j_ge1b",
        "_1m_ge7j_ge1b",
        "_1l_ge7j_ge1b_d1",
        "_1l_ge7j_ge1b_d2",
        "_1l_ge7j_ge1b_d3",
        "_1l_ge7j_ge1b_d4",
        //"_1l_ge7j_ge2b",                 
        "_1l_ge7j_ge1b_1t",
        //"_1l_ge7j_ge1b_ge1t",                                                    
        //"_1l_ge7j_ge1b_1t1",               
        //"_1l_ge7j_ge1b_1t2",               
        //"_1l_ge7j_ge1b_1t3",               
        //"_1l_ge7j_ge1b_1t2or3",            
        //"_1l_ge7j_ge1b_ge1t1",             
        //"_1l_ge7j_ge1b_ge1t2",             
        //"_1l_ge7j_ge1b_ge1t3",             
    };

    for(std::string mycut : mycuts_1l)
    {
        //plt.plotStack( "h_njets"+mycut, "N_{J}" , "Events", true);
        //plt.plotStack( "h_deepESM"+mycut, "DeepESM" , "Events", true, 10);
        plt.plotStack( "h_mbl"+mycut, "M(l,b) [GeV]",        "Events", true, 10);
        plt.plotStack( "h_lPt"+mycut, "Lepton P_{T} [GeV]",  "Events", true, 2, false, 0, 1000);
        plt.plotStack( "h_lEta"+mycut, "Lepton Eta",         "Events", true, 2);
        //plt.plotStack( "h_ntops"+mycut, "N_{T}" , "Events", true);
        //plt.plotStack( "h_nb"   +mycut, "N_{B}" , "Events", true);        
        //plt.plotStack( "h_fisher"+mycut, "fisher value" , "Events", true, 4);        
        //plt.plotStack( "h_photonPt"+mycut, "P_{T}^{#gamma}" , "Events", false, 10);
        //
        //// Make Normalized fisher
        //pltSkim.plotNormFisher("h_fisher_1l"+mycut, "fisher value" , "Events", false, 4);
        //plt.plotNormFisher("h_deepESM"+mycut, "DeepESM" , "Events", false, 10);
        //plt.plotNormFisher("h_njets"+mycut, "DeepESM" , "Events", true);
        //plt.plotNormFisher("h_BestComboMass_1l"+mycut, "Average BestCombo Mass [GeV]" , "Events", false, 4);
        //plt.plotNormFisher("h_BestComboPt_1l"+mycut, "Average BestCombo P_{T} [GeV]" , "Events", false, 4);
        //plt.plotNormFisher("h_BestComboMassDiff_1l"+mycut, "BestCombo Mass Diff [GeV]" , "Events", false, 2);
        //plt.plotNormFisher("h_BestComboMassDiffAbs_1l"+mycut, "BestCombo Mass Abs(Diff) [GeV]" , "Events", false, 2);
        //plt.plotNormFisher("h_BestComboRelDiff_1l"+mycut, "BestCombo Rel Diff" , "Events", false, 2);
        //plt.plotNormFisher("h_BestComboRelDiffAbs_1l"+mycut, "BestCombo Rel Abs(Diff)" , "Events", false, 2);
        //
        // - Make  Roc Curve
        //pltRoc.plotRocFisher("h_deepESM_1l"+mycut,"Background","Signal", false);
        //pltRocCompare.plotRocFisher(mycut,"Background","Signal", true);
        //        
        ////Need these until we un blind
        plt.plotStack( "blind_njets"+mycut,   "N_{J}",               "Events", true, -1, false);
        plt.plotStack( "blind_deepESM"+mycut, "DeepESM",             "Events", true, 10, false);
        plt.plotStack( "blind_mbl"+mycut,     "M(l,b) [GeV]",        "Events", true, 10, false);
        plt.plotStack( "blind_ht"+mycut,      "H_{T} [GeV]",         "Events", true, 10, false);
        plt.plotStack( "blind_ntops"+mycut,   "N_{T}",               "Events", true, -1, false);
        plt.plotStack( "blind_nb"   +mycut,   "N_{B}",               "Events", true, -1, false);        
        plt.plotStack( "blind_lPt"+mycut,     "Lepton P_{T} [GeV]",  "Events", true,  2, false, 0, 1000);
        plt.plotStack( "blind_lEta"+mycut,    "Lepton Eta",          "Events", true,  2);
        //plt.plotStack( "blind_fisher"+mycut, "fisher value" , "Events", true, 10);        
    }

    // --------------------
    // - Make fisher plots
    // --------------------

    std::vector< FisherHolder > fisherHolder
    {
        { {"h_njets_1l_ge7j_ge1b_d1"  , "h_njets_1l_ge7j_ge1b_d2"  , "h_njets_1l_ge7j_ge1b_d3"  , "h_njets_1l_ge7j_ge1b_d4"  } , "njets_1l_ge7j_ge1b"  },
        //{ {"h_njets_1l_ge7j_ge1b_ge6-7esm", "h_njets_1l_ge7j_ge1b_ge7-8esm", "h_njets_1l_ge7j_ge1b_ge8-95esm", "h_njets_1l_ge7j_ge1b_ge95esm" }, "njets_1l_ge7j_ge1b_esmBins"},
        //{ {"h_njets_1l_ge7j_ge1b_ge5-6esm", "h_njets_1l_ge7j_ge1b_ge6-7esm", "h_njets_1l_ge7j_ge1b_ge7-8esm", "h_njets_1l_ge7j_ge1b_ge8-95esm", "h_njets_1l_ge7j_ge1b_ge95esm" }, "njets_1l_ge7j_ge1b_esmBins"},
        //{ {"h_njets_1l_ge7j_ge1b_f1"  , "h_njets_1l_ge7j_ge1b_f2"  , "h_njets_1l_ge7j_ge1b_f3"  , "h_njets_1l_ge7j_ge1b_f4"  } , "njets_1l_ge7j_ge1b"  },
    };
    
    for (auto& f : fisherHolder)
    {
        //plt.plotFisher(f.cutNames_,  f.plotName_, "N_{J}", "Events", true, 9, "DeepESM");
        //plt.plotRatioFisher(f.cutNames_,  f.plotName_, "N_{J}", "N_{J+1} / N_{J}", false, 6, "DeepESM");
        plt.plotFisher(f.cutNames_,  f.plotName_, "N_{J}", "Events", true, 9, "Fisher");
        plt.plotRatioFisher(f.cutNames_,  f.plotName_, "N_{J}", "N_{J+1} / N_{J}", false, 6, "Fisher");
    }
    
    // --------------------
    // - Compute Yields
    // --------------------

    //for(std::string mycut : mycuts_1l)
    //{
    //    const auto& yieldMap = histInfoCollection.computeYields("h_njets_1l"+mycut,"njets",12,20);
    //    //const auto& yieldMap = histInfoCollection.computeYields("h_njets_1l"+mycut,"njets",6,8);
    //    auto allbg  = yieldMap.find("AllBG"); auto data_JetHT = yieldMap.find("Data_JetHT"); auto ttbar = yieldMap.find("T#bar{T}"); auto QCD = yieldMap.find("QCD");
    //    printf("%45.45s:  %s: %12.1lf   %s: %12.1lf   %s: %12.1lf   %s: %12.1lf\n",("h_njets_1l"+mycut).c_str(), (data_JetHT->first).c_str(), data_JetHT->second, (allbg->first).c_str(), allbg->second, "ttbar", ttbar->second, (QCD->first).c_str(), QCD->second );
    //}       
}
